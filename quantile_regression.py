import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
from sklearn.utils import shuffle
from numpy.random import choice
import statsmodels.api as sm
import statsmodels.formula.api as smf
import argparse

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from file_utils import *
from df_utils import *


def scatter_plot(x, y, title=None, fontsize=12, dotsize=3, color=None,
                 xlabel=None, ylabel=None, xlim=[-0.05, 1], ylim=[-0.05, 1],
                 outfname=None, save_svg=False, save_jpg=False, transparent=True):
    plt.figure(figsize=(6,6))
    if color is None:
        plt.scatter(x, y, s=dotsize)
    else:
        plt.scatter(x, y, s=dotsize, c=color)
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    if xlabel != None:
        plt.xlabel(xlabel, fontsize=fontsize)
    if ylabel != None:
        plt.ylabel(ylabel, fontsize=fontsize)
    if title != None:
        plt.title(title)
    sns.despine()
    if save_svg:
        if outfname[-4:] != '.svg':
            outfname += '.svg'
        plt.savefig(outfname, format='svg', transparent=transparent)
    if save_jpg:
        if outfname[-4:] != '.jpg':
            outfname += '.jpg'
        plt.savefig(outfname, format='jpg', transparent=transparent)
    plt.show()
    return

parser = argparse.ArgumentParser()
parser.add_argument('--q', type=float)
parser.add_argument('--balance', action='store_true')
parser.add_argument('--outdir')
parser.add_argument('--rf_pred_path')
parser.add_argument('--sdir')
parser.add_argument('--hierarchy_path')
parser.add_argument('--predicted_sim_path')
args = parser.parse_args()

outdir = args.outdir
rf_pred = np.load(args.rf_pred_path)
ytrue = np.load('{}/y_test_genepair.all.npy'.format(args.sdir))

bin_size = 0.01
bin_start = np.arange(0, 1, bin_size, dtype=float)
bin_end = np.arange(bin_size, 1+bin_size, bin_size, dtype=float)
# # start <= x < end
# # Add a small value to the end position of last bin to include everything
bin_end[-1] += 0.01

bin_start = np.round(bin_start, 3)
bin_end = np.round(bin_end, 3)
bins = list(zip(bin_start, bin_end))

N = 1000000 # maximum number of training samples
n = int(np.floor(N / len(bins))) # number of training samples per resnik bin
print('Each bin has {} samples'.format(n))

sample_index = np.asarray([], dtype=int)
for b in bins:
     ridx = np.where((ytrue >= b[0]) & (ytrue < b[1]))[0]
     if len(ridx) == 0:
         continue
     if len(ridx) < n:
         rand_idx = choice(ridx, n, replace=True)
     else:
         rand_idx = choice(ridx, n, replace=False)
     sample_index = np.concatenate((sample_index, rand_idx))
# np.save('{}/sample_idx.balanced.npy'.format(outdir), sample_index)
# # Shuffle training data
np.random.shuffle(sample_index)
np.save('{}/sample_idx.balanced.shuffled.npy'.format(outdir), sample_index)

data = pd.DataFrame()
if args.balance:
    data['rf_pred'] = rf_pred[sample_index]
    data['ytrue'] = ytrue[sample_index]
    outprefix = '{}/{}_balanced'.format(outdir, args.q)
else:
    data['rf_pred'] = rf_pred
    data['ytrue'] = ytrue
    outprefix = '{}/{}'.format(outdir, args.q)

mod = smf.quantreg('ytrue ~ rf_pred', data)
res = mod.fit(q=args.q)
with open('{}.summary.txt'.format(outprefix), 'w') as outf:
    outf.write(res.summary().as_text())
save_obj(res, '{}.summary.pkl'.format(outprefix))

res_rf_pred = res.params['rf_pred'] * rf_pred + res.params['Intercept']
np.save('{}.res_rf_pred.npy'.format(outprefix), res_rf_pred)

scatter_plot(rf_pred, res_rf_pred, dotsize=1,
                  ylabel='Predicted similarity after Quantile regression',
                  xlabel='Predicted similarity',
                  outfname='{}.rfPred_VS_resRfPred.jpg'.format(outprefix), save_jpg=True)

scatter_plot(res_rf_pred, ytrue, dotsize=1,
                  xlabel='Predicted similarity after Quantile regression', 
                  ylabel='Resnik similarity', 
                  outfname='{}.resRfPred_VS_ytrue.jpg'.format(outprefix), save_jpg=True)

plt.figure(figsize=(6,6))
sns.kdeplot(x=res_rf_pred, y=ytrue, shade=True)
sns.despine()
plt.xlabel('Predicted similarity after Quantile regression')
plt.ylabel('Resnik similarity')
plt.savefig('{}.resRfPred_VS_ytrue.kde.png'.format(outprefix), format='png')


res_rf_pred = np.load('{}/0.2_balanced.res_rf_pred.npy'.format(outdir))
res = load_obj('{}/0.2_balanced.summary.pkl'.format(outdir))

def predict(res, x):
    return res.params['rf_pred'] * x + res.params['Intercept']

def resnik_to_nm(resnik):
    power = -3.968*resnik + 4.326
    return 10**power

ont = pd.read_table(args.hierarchy_path)
ont.columns = ['cluster_id', 'size', 'genes']

ddot = pd.read_table(args.predicted_sim_path,
                     header=None)
ddot.columns = ['geneA', 'geneB', 'sim']
ddot['lt_sim'] = predict(res, ddot['sim'].values)
ddot['nm'] = resnik_to_nm(ddot['lt_sim'].values)
median_nm = [] # Median
for comp, row in ont.iterrows():
    comp_genes = row['genes'].split(',')[:-1]
    comp_ddot = ddot[ddot['geneA'].isin(comp_genes) & ddot['geneB'].isin(comp_genes)]
    median_nm.append(np.median(comp_ddot['nm'].values))
ont['median_recal_nm'] = median_nm
sphere_constant = 35/18
ont['estimated_size'] = sphere_constant*ont['median_recal_nm']
std_error = 0.1251855862826624
t = 2.306 # 95% interval with df = 8
ont['log_estimated_size'] = np.log10(ont['estimated_size'])
for idx, row in ont.iterrows():
    ont.at[idx, 'lower_bound'] = 10**(row['log_estimated_size'] - t*std_error)
    ont.at[idx, 'upper_bound'] = 10**(row['log_estimated_size'] + t*std_error)
output_path = args.hierarchy_path.strip('termStats') + 'calibrated.termStats'
ont.to_csv(output_path, index=False, sep='\t')

print('=== finished! ===')

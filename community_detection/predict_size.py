import pandas as pd
import numpy as np
import sys
import os
import argparse
import statsmodels.api as sm
import statsmodels.formula.api as smf
from numpy.random import choice
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from music_utils import *


parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--n_samples', default=1000000, type=int,
                    help='Number of samples to use for fitting linear model.')
parser.add_argument('--q_step', type=float, default=0.1, help='Step for scanning best quantile to use.')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology.')
args = parser.parse_args()

outprefix = args.outprefix
q_step = args.q_step
# Check if all required files exist
if not os.path.exists('{}.louvain.termStats'.format(outprefix)):
    raise ValueError('{}.louvain.termStats file does not exist.'.format(outprefix))
if not os.path.exists('{}_avgPred_ytrue.csv'.format(outprefix)):
    raise ValueError('{}_avgPred_ytrue.csv file does not exist. Please run random_forest_output.py script with the same outprefix first.'.format(outprefix))

qr_dir = '{}_qr'.format(outprefix)
if not os.path.exists(qr_dir):
    os.mkdir(qr_dir)
df = pd.read_csv('{}_avgPred_ytrue.csv'.format(outprefix))

# Balance samples used for quantile regression
bin_size = 0.01
bin_start = np.arange(0, 1, bin_size, dtype=float)
bin_end = np.arange(bin_size, 1+bin_size, bin_size, dtype=float)
# start <= x < end
# Add a small value to the end position of last bin to include everything
bin_end[-1] += 0.01

bin_start = np.round(bin_start, 3)
bin_end = np.round(bin_end, 3)
bins = list(zip(bin_start, bin_end))

N = args.n_samples
n = int(np.floor(N / len(bins))) 

sample_index = np.asarray([], dtype=int)
for b in bins:
    ridx = df[(df['ytrue'] >= b[0]) & (df['ytrue'] < b[1])].index.values
    if len(ridx) == 0:
        continue
    if len(ridx) < n:
        rand_idx = choice(ridx, n, replace=True)
    else:
        rand_idx = choice(ridx, n, replace=False)
    sample_index = np.concatenate((sample_index, rand_idx))
print('Sampled {} data points for calibration'.format(len(sample_index)))
# Shuffle training data
np.random.shuffle(sample_index)
np.save('{}/sample_idx.balanced.shuffled.npy'.format(qr_dir), sample_index, allow_pickle=True)

data = df.loc[sample_index]
print('Start fitting linear model...')
q_values = np.round(np.arange(q_step, 1, q_step), 3)
q_stats = np.asarray([-10]*len(q_values), dtype=float)
for i, q in enumerate(q_values):
    print('testing quantile {} ...'.format(q))
    mod = smf.quantreg('ytrue ~ avg_pred', data)
    res = mod.fit(q=q)
    q_stats[i] = res.prsquared
stats_df = pd.DataFrame()
stats_df['q'] = q_values
stats_df['prsquared'] = q_stats
best_q = stats_df.loc[stats_df['prsquared'].idxmax()]['q']
print('Best quantile for linear fit: {}'.format(best_q))
res = mod.fit(q=best_q)
res_rf_pred = res.params['avg_pred'] * df['avg_pred'].values + res.params['Intercept']
df['linear_adjusted'] = res_rf_pred
df['nm'] = [scaled_P_to_nm(x) for x in res_rf_pred]

ont = pd.read_table('{}.louvain.termStats'.format(outprefix), header=None, index_col=0)
ont.columns = ['Number_of_proteins', 'Proteins']
median_nm = []
for comp, row in ont.iterrows():
    comp_genes = row['Proteins'].split(',')[:-1]
    comp_df = df[df['geneA'].isin(comp_genes) & df['geneB'].isin(comp_genes)]
    median_nm.append(np.median(comp_df['nm'].values))
ont['median_recal_nm'] = median_nm
cal_constant = 1.7905959998061936
ont['Estimated_size_in_nm'] = cal_constant*ont['median_recal_nm']

# Propagate to ensure parent larger than child
cmd = '{}/ontologyTermStats {}.louvain.ddot descendents > {}.louvain.descendents'.format(args.path_to_alignOntology.rstrip('/'), 
                                                                                         outprefix, outprefix)
os.system(cmd)

descendent = pd.read_table('{}.louvain.descendents'.format(outprefix), header=None, index_col=0)
descendent.columns = ['child']
for idx, row in ont.iterrows():
    child_idx = descendent.at[idx, 'child'].split(',')[:-1]
    ont.at[idx, 'Estimated_size_in_nm'] = max(ont.loc[child_idx]['Estimated_size_in_nm'].values)
ont.to_csv('{}.louvain.termStats'.format(outprefix), sep='\t')
print('=== finished predict_size ===')
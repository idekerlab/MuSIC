import numpy as np
import pandas as pd
import sys
import os
from glob import glob
import argparse
from music_utils import *

parser = argparse.ArgumentParser(description='Organize output from random forest regressors.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
args = parser.parse_args()

if args.outprefix is None:
    raise ValueError('Please specify --outprefix: Prefix for files generated. E.g., /path/to/output/directory/fileIdentifier')
outprefix = args.outprefix
# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)

folds = [x.split('.')[-2].split('_')[-1] for x in glob('{}_train_test_data/X_test_genepair_*.npy'.format(outprefix))]
X_test_genepair = np.load('{}_train_test_data/X_test_genepair_{}.npy'.format(outprefix, folds[0]), allow_pickle=True)
y_test_genepair = np.load('{}_train_test_data/y_test_genepair_{}.npy'.format(outprefix, folds[0]), allow_pickle=True)
# Average all random forest output for a protein pair as the final protein-protein proximity
pred_files = glob('{}_yPred_output/*.fold_{}.npy'.format(outprefix, folds[0]))
if len(pred_files) == 0:
    raise ValueError('Fold {} does not have ypred files.'.format(fold))
fold_ypred = []
for f in pred_files:
    fold_ypred.append(np.load(f, allow_pickle=True))
fold_ypred = np.asarray(fold_ypred)
avg_ypred = fold_ypred.mean(axis=0)

for i in range(1, len(folds)):
    fold_X_test = np.load('{}_train_test_data/X_test_genepair_{}.npy'.format(outprefix, folds[i]), allow_pickle=True)
    fold_y_test = np.load('{}_train_test_data/y_test_genepair_{}.npy'.format(outprefix, folds[i]), allow_pickle=True)
    X_test_genepair = np.concatenate((X_test_genepair, fold_X_test))
    y_test_genepair = np.concatenate((y_test_genepair, fold_y_test))
    # Average all random forest output for a protein pair as the final protein-protein proximity
    pred_files = glob('{}_yPred_output/*.fold_{}.npy'.format(outprefix, folds[i]))
    if len(pred_files) == 0:
        raise ValueError('Fold {} does not have ypred files.'.format(fold))
    fold_ypred = []
    for f in pred_files:
        fold_ypred.append(np.load(f, allow_pickle=True))
    fold_ypred = np.asarray(fold_ypred)
    avg_ypred = np.concatenate((avg_ypred, fold_ypred.mean(axis=0)))

df = pd.DataFrame(X_test_genepair, columns=['geneA', 'geneB'])
df['avg_pred'] = avg_ypred
df['ytrue'] = y_test_genepair
df.to_csv('{}_avgPred_ytrue.csv'.format(outprefix), index=False)

# Prepare input for community detection
rest_genepair = np.load('{}_rest_genepair.npy'.format(outprefix), allow_pickle=True)
X_test_genepair = np.concatenate((X_test_genepair, rest_genepair))
rest_pred_files = glob('{}_yPred_output/*.restGP.npy'.format(outprefix))
rest_ypred = []
for f in rest_pred_files:
    rest_ypred.append(np.load(f, allow_pickle=True))
rest_ypred = np.asarray(rest_ypred)
avg_ypred = np.concatenate((avg_ypred, rest_ypred.mean(axis=0)))
df = pd.DataFrame(X_test_genepair, columns=['geneA', 'geneB'])
df['avg_pred'] = avg_ypred
df.to_csv('{}_predicted_proximity.txt'.format(outprefix), sep='\t', index=False, header=False)

print('=== finished! ===')
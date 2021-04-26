from sklearn.ensemble import RandomForestRegressor
import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr
import os
import argparse
from music_utils import *

parser = argparse.ArgumentParser(description='Generate train/test samples for random forest regressors.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--fold', type=int, help='Specify which fold of k-fold cross validation to train.')
parser.add_argument('--emd_label', nargs='+', help='Label for each embedding file.')
parser.add_argument('--train_set', nargs='+', help='For each embedding data, specify which training set to use. In order with <emd_label>. Default to 1 for each emd_label.')
parser.add_argument('--n_estimators', type=int, default=1000, help='The number of trees in the forest.')
parser.add_argument('--max_depth', type=int, default=30, help='The maximum depth of the tree.')
parser.add_argument('--n_jobs', type=int, default=8, help='The number of jobs to run in parallel.')
args = parser.parse_args()

if args.outprefix is None:
    raise ValueError('Please specify --outprefix: Prefix for files generated. E.g., /path/to/output/directory/fileIdentifier')
outprefix = args.outprefix
# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)

fold = args.fold
max_depth = args.max_depth
n_estimators = args.n_estimators
n_jobs = args.n_jobs
emd_label = args.emd_label
if args.train_set is None:
    train_set = [1] * len(emd_label)
else:
    train_set = args.train_set
    if len(train_set) != len(emd_label):
        raise ValueError('Entered train_set does not match emd_label.')
# Check if train_set for each emd_label exists
tmp_label = []
for i in range(len(emd_label)):
    tmp_label.append('{}_{}'.format(emd_label[i], train_set[i]))
    if not os.path.exists('{}_{}/training_set_{}'.format(outprefix, emd_label[i], train_set[i])):
        raise ValueError('{} does not have training set {}'.format(emd_label[i], train_set[i]))

outfname = '{}.{}.RF_maxDep_{}_nEst_{}.fold_{}.pkl'.format(outprefix, '_'.join(x for x in tmp_label), 
                                                           max_depth, n_estimators, fold)

print('Trained model will be saved in {}'.format(outfname))
print('Start training random forest model...')
# Build Random Forest Model
rf = RandomForestRegressor(max_depth=max_depth,
                           n_estimators=n_estimators,
                           n_jobs=n_jobs)
print(rf)

sample_dir = '{}_train_test_data'.format(outprefix)
train_idx = np.load('{}/train_idx_{}.balanced.shuffled.npy'.format(sample_dir, fold))
# Load y_train data
y_train = np.load('{}/y_train_genepair_{}.npy'.format(sample_dir, fold))[train_idx]
# Load training data for each data modality
X_train = np.load('{}_{}/training_set_{}/X_train_{}.by_gp.npy'.format(outprefix, 
                                                                      emd_label[0], 
                                                                      train_set[0], 
                                                                      fold))[train_idx]
emd_dim = [X_train.shape[1]]
for i in range(1, len(emd_label)):
    emd_X_train = np.load('{}_{}/training_set_{}/X_train_{}.by_gp.npy'.format(outprefix, 
                                                                              emd_label[i], 
                                                                              train_set[i], 
                                                                              fold))[train_idx]
    emd_dim.append(emd_X_train.shape[1])
    X_train = np.concatenate((X_train, emd_X_train), axis=1)
print('... loaded training data')
if len(set(emd_dim)) > 1:
    print('NOTE: we recommend using same number of features for each data modality.')
print('Start training random forest regressor on {} samples with {} total features...'.format(X_train[0], X_train[1]))


rf.fit(X_train, y_train)
print('... finished training')

save_obj(rf, outfname)
print('... saved trained random forest regressor')

# Predict test data
print('\nStart predicting test data...')
X_test = np.load('{}_{}/training_set_{}/X_test_{}.npy'.format(outprefix,
                                                              emd_label[0], 
                                                              train_set[0], fold))
for i in range(1, len(emd_label)):
    emd_X_test = np.load('{}_{}/training_set_{}/X_test_{}.npy'.format(outprefix, 
                                                                      emd_label[i], 
                                                                      train_set[i], fold))
    X_test = np.concatenate((X_test, emd_X_test), axis=1)
print('... loaded {} test samples'.format(X_test.shape[0]))
y_pred = rf.predict(X_test)
pred_dir = '{}_yPred_output'.format(outprefix)
if not os.path.exists(pred_dir):
    os.mkdir(pred_dir)
np.save('{}/{}.npy'.format(pred_dir, outfname.split('/')[-1]), y_pred)
print('... finished predicting test samples')
y_test = np.load('{}/y_test_genepair_{}.npy'.format(sample_dir, fold))
print('Pearson r: {}'.format(pearsonr(y_test, y_pred)[0]))

# Predict gene pairs without specific GO annotations
print('\nStart predicting gene pairs without specific GO annotations...')
X_rest = np.load('{}_{}/training_set_{}/rest_genepair.npy'.format(outprefix,
                                                                  emd_label[0], 
                                                                  train_set[0]))
for i in range(1, len(emd_label)):
    emd_X_rest = np.load('{}_{}/training_set_{}/rest_genepair.npy'.format(outprefix, 
                                                                          emd_label[i], 
                                                                          train_set[i]))
    X_rest = np.concatenate((X_rest, emd_X_rest), axis=1)
y_rest = rf.predict(X_rest)
np.save('{}/{}.restGP.npy'.format(pred_dir, outfname.split('/')[-1]), y_rest)
print('=== finished! ===')
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import pandas as pd
import sys
import os
import argparse
from file_utils import *


def train_RF(arg_max_features, batch, img_fold, img_label, ppi_label, max_depth, n_estimators, sdir, outdir):
    # format max_features
    if arg_max_features == 'auto':
        max_features = 'auto'
    else:
        max_features = int(arg_max_features)

    # format batch
    print('... train on batch {}'.format(batch))

    # format output filename
    if img_fold == -1:
        outfname = '{}/{}.{}.RF_maxDep_{}_nEst_{}_maxFeat_{}_batch_{}.pkl'.format(outdir, img_label,
                                                                                  ppi_label, max_depth,
                                                                                  n_estimators, max_features,
                                                                                  batch)
    else:
        outfname = '{}/{}.{}.RF_maxDep_{}_nEst_{}_maxFeat_{}_batch_{}.fold_{}.pkl'.format(outdir, img_label,
                                                                                          ppi_label,
                                                                                          max_depth,
                                                                                          n_estimators,
                                                                                          max_features,
                                                                                          batch, img_fold)
    print(outfname)

    # Build Random Forest Model
    rf = RandomForestRegressor(max_depth=max_depth,
                               n_estimators=n_estimators,
                               max_features=max_features,
                               n_jobs=24)
    print(rf)

    train_idx = np.load('{}/train_idx_{}.balanced.shuffled.npy'.format(sdir, batch))
    # Load y_train data
    y_train = np.load('{}/y_train_genepair_{}.npy'.format(sdir, batch))[train_idx]
    # Load Image training data
    if img_fold == -1:
        img_X_train = np.load('{}/X_train_{}.by_gp.npy'.format(sdir, batch))[train_idx]
    else:
        img_X_train = np.load('{}/X_train_{}.by_gp.fold_{}.npy'.format(sdir, batch, img_fold))[train_idx]
    # Load PPI training data
    ppi_X_train = np.load('{}/{}_X_train_{}.by_gp.npy'.format(sdir, ppi_label, batch))[train_idx]
    # Concatenate Image and PPI data to generate training data
    X_train = np.concatenate((img_X_train, ppi_X_train), axis=1)

    print('... loaded data')
    print(X_train.shape)
    print(y_train.shape)

    print('... start training Random Forest model')
    rf.fit(X_train, y_train)
    print('... finished training')

    save_obj(rf, outfname)

    print('=== finished! ===')


parser = argparse.ArgumentParser(description='Arguments for Random Forest Regressor')
parser.add_argument('--sdir', help='directory with train and test data')
parser.add_argument('--n_estimators', type=int)
parser.add_argument('--max_depth', type=int)
parser.add_argument('--max_features')
parser.add_argument('--n2v_prefix', help='the label for node2vec feature')
parser.add_argument('--densenet_prefix', help='the label for image feature')
parser.add_argument('--outputdir', help='output directory')
args = parser.parse_args()

for batch in range(1, 6):
    for fold in range(1,7):
        train_RF(args.max_features, batch, fold, args.densenet_prefix, args.n2v_prefix, args.max_depth, args.n_estimators, args.sdir, args.outputdir)


import argparse
import pandas as pd
import numpy as np
import sys
import os
from file_utils import *
from df_utils import *

def format_x_and_y_data(workdir):

    with open('{}/hek_ccNoHPA_resnik.602.scaled.pkl'.format(workdir), 'rb') as f:
        resnik = pickle.load(f)
    genes = resnik.index.values

    if not os.path.exists('{}/train_test_data'.format(workdir)):
        os.mkdir('{}/train_test_data'.format(workdir))
    outputdir = '{}/train_test_data'.format(workdir)

    # make gene pair files
    X = []
    y = []
    for i in range(len(genes) - 1):
        for j in range(i + 1, len(genes)):
            X.append(tuple([genes[i], genes[j]]))
            y.append(resnik.at[genes[i], genes[j]])
    X = np.asarray(X)
    y = np.asarray(y)
    np.save('{}/X_genepairs.npy'.format(outputdir), X)
    np.save('{}/y_genepairs.npy'.format(outputdir), y)

    # split gene pairs into 5 batches
    kf = KFold(n_splits=5, random_state=10, shuffle=True)
    batch = 0
    y_test_all = list()
    for train_index, test_index in kf.split(X):
        batch += 1
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        np.save('{}/X_train_genepair_{}.npy'.format(outputdir, batch), X_train)
        np.save('{}/X_test_genepair_{}.npy'.format(outputdir, batch), X_test)
        np.save('{}/y_train_genepair_{}.npy'.format(outputdir, batch), y_train)
        np.save('{}/y_test_genepair_{}.npy'.format(outputdir, batch), y_test)
        y_test_all.extend(y_test)
    np.save('{}/y_test_genepair.all.npy'.format(outputdir), y_test_all)

    for batch in range(1, 6):
        y_train = np.load('{}/y_train_genepair_{}.npy'.format(outputdir, batch))

        # Balance each batch by bins
        bin_size = 0.01
        bin_start = np.arange(0, 1, bin_size, dtype=float)
        bin_end = np.arange(bin_size, 1 + bin_size, bin_size, dtype=float)
        # start <= x < end
        # Add a small value to the end position of last bin to include everything
        bin_end[-1] += 0.01

        bin_start = np.round(bin_start, 3)
        bin_end = np.round(bin_end, 3)
        bins = list(zip(bin_start, bin_end))

        N = 1000000  # maximum number of training samples
        n = int(np.floor(N / len(bins)))  # number of training samples per resnik bin
        print('Each bin has {} samples'.format(n))

        sample_index = np.asarray([], dtype=int)
        for b in bins:
            ridx = np.where((y_train >= b[0]) & (y_train < b[1]))[0]
            if len(ridx) == 0:
                continue
            if len(ridx) < n:
                rand_idx = choice(ridx, n, replace=True)
            else:
                rand_idx = choice(ridx, n, replace=False)
            sample_index = np.concatenate((sample_index, rand_idx))
        print('Sampled {} training samples for batch {}'.format(len(sample_index), batch))
        np.save('{}/train_idx_{}.balanced.npy'.format(outputdir, batch), sample_index)
        # Shuffle training data
        np.random.shuffle(sample_index)
        np.save('{}/train_idx_{}.balanced.shuffled.npy'.format(outputdir, batch), sample_index)

import numpy as np
import sys
import os
from numpy.random import choice
import argparse
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from music_utils import *

def balance_sample(sample_dir, batch, N):
    '''
    Balance the training samples for fold of K-fold cross validation.
    
    Args:
        sample_dir: Path to folder containing training and testing protein pairs.
        batch: The fold to balance.
        N: Maximum number of training samples for each fold.
    '''
    y_train = np.load('{}/y_train_genepair_{}.npy'.format(sample_dir, batch))

    bin_size = 0.01
    bin_start = np.arange(0, 1, bin_size, dtype=float)
    bin_end = np.arange(bin_size, 1+bin_size, bin_size, dtype=float)
    # start <= x < end
    # Add a small value to the end position of last bin to include everything
    bin_end[-1] += 0.01
    bin_start = np.round(bin_start, 3)
    bin_end = np.round(bin_end, 3)
    bins = list(zip(bin_start, bin_end))
    
    # number of training samples per bin
    n = int(np.floor(N / len(bins)))

    sample_index = np.asarray([], dtype=int)
    for b in bins:
        ridx = np.where((y_train >= b[0]) & (y_train < b[1]))[0]
        if len(ridx) == 0:
            continue
        if len(ridx) < n:
            rand_idx = choice(ridx, n, replace=True) # Upsample
        else:
            rand_idx = choice(ridx, n, replace=False) # Downsample
        sample_index = np.concatenate((sample_index, rand_idx))
    print('Sampled {} training samples for fold {}'.format(len(sample_index), batch))
    # Shuffle training data
    np.random.shuffle(sample_index)
    np.save('{}/train_idx_{}.balanced.shuffled.npy'.format(sample_dir, batch), sample_index)
    return
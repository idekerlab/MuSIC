import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from sklearn.model_selection import KFold
from numpy.random import choice
from random import randint
import random
import argparse
cdDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cdDir + '/random_forest')
from music_utils import *
from format_emd_data import *
from balance_sample import *

parser = argparse.ArgumentParser(description='Generate train/test samples for random forest regressors.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--protein_file', help='Path to the file containing list of proteins to analyze. E.g., /Examples/MuSIC_proteins.txt')
parser.add_argument('--emd_files', nargs='+', help='Path to each embedding file.')
parser.add_argument('--emd_label', nargs='+', help='Label for each embedding file. In order with <emd_files>.')
parser.add_argument('--num_set', nargs='+', help='Number of training sets for each embedding file. In order with <emd_files>. Default to auto')
parser.add_argument('--n_samples', type=int, default=1000000, 
                    help='Maximum number of samples to train/test random forest regressor.')
parser.add_argument('--k', type=int, default=5,
                    help='Specify k for k-fold cross-validation.')
args = parser.parse_args()

if args.outprefix is None:
    raise ValueError('Please specify --outprefix: Prefix for files generated. E.g., /path/to/output/directory/fileIdentifier')
outprefix = args.outprefix
# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Check if protein_file exists
if not os.path.exists(args.protein_file):
    raise ValueError('File path entered for protein_file does not exist!')
proteins = [x.rstrip('\n') for x in open(args.protein_file, 'r').readlines()]

# Sanity check on entered embedding files
if args.emd_files is None:
    raise ValueError('Please enter path to each embedding file.')
if args.emd_label is None:
    raise ValueError('please specify a unique identifier for each entered embedding file. E.g., IF_emd, APMS_emd')
if len(args.emd_files) != len(args.emd_label):
    raise ValueError('Number of elements entered for <emd_files> and <emd_label> does not match.')
if args.num_set != None:
    if len(args.emd_files) != len(args.num_set):
        raise ValueError('Number of elements entered for <emd_files> and <num_set> does not match.')
        
# Check if calibraiton output files exist
if not os.path.exists('{}.annot_proteins.npy'.format(outprefix)):
    raise ValueError('Missing annot_proteins.npy file. Please run calibrate_pairwise_distance.py first using the same outprefix.')
if not os.path.exists('{}.not_annot_proteins.npy'.format(outprefix)):
    raise ValueError('Missing not_annot_proteins.npy file. Please run calibrate_pairwise_distance.py first using the same outprefix.')
if not os.path.exists('{}.calibrated_distance.csv'.format(outprefix)):
    raise ValueError('Missing calibrated_distance.csv file. Please run calibrate_pairwise_distance.py first using the same outprefix.')
     

# Format each embedding file into specified training sets and calculate pairwise protein similarity.
for i in range(len(args.emd_files)):
    if args.num_set is None:
        print('Formatting {} using file {} ...'.format(args.emd_label[i], args.emd_files[i]))
        format_emd(args.emd_files[i], args.emd_label[i], proteins, outprefix)
    else:
        print('Formatting {} into {} training sets using file {} ...'.format(args.emd_label[i], 
                                                                            args.num_set[i],
                                                                            args.emd_files[i]))
        format_emd(args.emd_files[i], args.emd_label[i], proteins, outprefix, num_set=args.num_set[i])
    print('\n')

# Create directory to store training and testing samples
sample_dir = '{}_train_test_data'.format(outprefix)
if not os.path.exists(sample_dir):
    os.mkdir(sample_dir)
# Scale protein-protein proximity into [0, 1]
df = pd.read_csv('{}.calibrated_distance.csv'.format(outprefix))
df['scaled_P'] = df['P'] - df['P'].min()
max_p = df['scaled_P'].max()
df['scaled_P'] = df['scaled_P'] / max_p
# Format train-test-split for k-fold cross validation
print('Split protein pairs into training and testing sets for {}-fold cross validation...'.format(args.k))
X = df[['geneA', 'geneB']].values
y = df['scaled_P'].values
kf = KFold(n_splits=args.k, shuffle=True)
batch = 0
for train_index, test_index in kf.split(X):
    batch += 1
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    np.save('{}/X_train_genepair_{}.npy'.format(sample_dir, batch), X_train, allow_pickle=True)
    np.save('{}/X_test_genepair_{}.npy'.format(sample_dir, batch), X_test, allow_pickle=True)
    np.save('{}/y_train_genepair_{}.npy'.format(sample_dir, batch), y_train, allow_pickle=True)
    np.save('{}/y_test_genepair_{}.npy'.format(sample_dir, batch), y_test, allow_pickle=True)
# Balance training sample for each fold
for fold in range(1, args.k+1):
    balance_sample(sample_dir, fold, args.n_samples)
    
# Format gene pairs that lack specific annotation in Gene Ontology
rest_gp = []
annot_proteins = np.load('{}.annot_proteins.npy'.format(outprefix), allow_pickle=True)
unannot_proteins = np.load('{}.not_annot_proteins.npy'.format(outprefix), allow_pickle=True)
for i in range(len(unannot_proteins)-1):
    for j in range(i+1, len(unannot_proteins)):
        rest_gp.append([unannot_proteins[i], unannot_proteins[j]])
for i in range(len(unannot_proteins)):
    for j in range(len(annot_proteins)):
        rest_gp.append([unannot_proteins[i], annot_proteins[j]])
rest_gp = np.asarray(rest_gp)
np.save('{}_rest_genepair.npy'.format(outprefix), rest_gp, allow_pickle=True)

# Generate protein pairwise features input to random forest
print('Start generating protein pairwise features input to random forest...')
for i in range(len(args.emd_files)):
    avail_train_set = glob('{}_{}/training_set_*'.format(outprefix, args.emd_label[i]))
    avail_train_set = [x.split('_')[-1] for x in avail_train_set]
    for fold in range(1, args.k+1):
        for train_set in avail_train_set:
            get_emd_X(outprefix, fold, args.emd_files[i], args.emd_label[i], train_set, rest_gp)

print('=== finished! ===')
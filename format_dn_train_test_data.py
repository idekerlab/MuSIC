import argparse
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
from file_utils import *
from df_utils import *

parser = argparse.ArgumentParser(description='Format training and testing data for Node2vec')
parser.add_argument('--prefix', help='Image prefix')
parser.add_argument('--fold', help='The fold of image ID to use')
parser.add_argument('--batch', help='Data batch to use')
parser.add_argument('--emdfile', help='path to feature embedding file')
parser.add_argument('--sdir', help='directory for X train and test data')
parser.add_argument('--simdir', help='directory for similarity files')
parser.add_argument('--nodeiddir', help='directory for gname to node id files')
parser.add_argument('--outputdir', help='output directory')
args = parser.parse_args()

prefix = args.prefix
fold = args.fold
batch = args.batch

sdir = args.sdir
X_train_gp = np.load('{}/X_train_genepair_{}.npy'.format(sdir, batch))
X_test_gp = np.load('{}/X_test_genepair_{}.npy'.format(sdir, batch))

sdir = args.sdir
simdir = args.simdir
outputdir = args.outputdir
nodeiddir = args.nodeiddir
cosine = load_obj('{}/{}.cosine.scaled.pkl'.format(simdir, prefix))
pearson = load_obj('{}/{}.pearson.scaled.pkl'.format(simdir, prefix))
spearman = load_obj('{}/{}.spearman.scaled.pkl'.format(simdir, prefix))
kendall = load_obj('{}/{}.kendall.scaled.pkl'.format(simdir, prefix))
manhattan = load_obj('{}/{}.manhattan.scaled.pkl'.format(simdir, prefix))
euclidean = load_obj('{}/{}.euclidean.scaled.pkl'.format(simdir, prefix))
imgID_to_gname = load_obj('{}/imgID_to_gname.dict.pkl'.format(nodeiddir))
imgID = [x for x in imgID_to_gname]
gname = [imgID_to_gname[x] for x in imgID_to_gname]

emd = load_obj(args.emdfile).loc[imgID]
emd.index = gname

X_train = []
for gp in X_train_gp:
    ga, gb = gp
    feature = []
    feature.append(cosine.at[ga, gb])
    feature.append(pearson.at[ga, gb])
    feature.append(spearman.at[ga, gb])
    feature.append(kendall.at[ga, gb])
    feature.append(manhattan.at[ga, gb])
    feature.append(euclidean.at[ga, gb])
    feature += list(np.abs(emd.loc[ga].values - emd.loc[gb].values))
    X_train.append(feature)
X_train = np.asarray(X_train)
np.save('{}/X_train_{}.by_gp.fold_{}.npy'.format(outputdir, batch, fold), X_train)
print('... finished formatting and saving training data')

X_test = []
for gp in X_test_gp:
    ga, gb = gp
    feature = []
    feature.append(cosine.at[ga, gb])
    feature.append(pearson.at[ga, gb])
    feature.append(spearman.at[ga, gb])
    feature.append(kendall.at[ga, gb])
    feature.append(manhattan.at[ga, gb])
    feature.append(euclidean.at[ga, gb])
    feature += list(np.abs(emd.loc[ga].values - emd.loc[gb].values))
    X_test.append(feature)
X_test = np.asarray(X_test)
np.save('{}/X_test_{}.fold_{}.npy'.format(outputdir, batch, fold), X_test)
print('... finished formatting and saving testing data')
print('=== finished! ===')

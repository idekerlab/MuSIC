import argparse
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
from file_utils import *
from df_utils import *

parser = argparse.ArgumentParser(description='Format data for genes not in GO')
parser.add_argument('--prefix', help='Image prefix')
parser.add_argument('--fold', help='The fold of image ID to use')
parser.add_argument('--emdfile', help='path to feature embedding file')
parser.add_argument('--sdir', help='directory for X train and test data')
parser.add_argument('--simdir', help='directory for similarity files')
parser.add_argument('--nodeiddir', help='directory for gname to node id files')
parser.add_argument('--outputdir', help='output directory')
args = parser.parse_args()

prefix = args.prefix
fold = args.fold

sdir = args.sdir
rest_gp = np.load('{}/rest_genepair.npy'.format(sdir))

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

gp_data = []
for gp in rest_gp:
    ga, gb = gp
    feature = []
    feature.append(cosine.at[ga, gb])
    feature.append(pearson.at[ga, gb])
    feature.append(spearman.at[ga, gb])
    feature.append(kendall.at[ga, gb])
    feature.append(manhattan.at[ga, gb])
    feature.append(euclidean.at[ga, gb])
    feature += list(np.abs(emd.loc[ga].values - emd.loc[gb].values))
    gp_data.append(feature)
gp_data = np.asarray(gp_data)
np.save('{}/rest_genepair.fold_{}.npy'.format(outputdir, fold), gp_data)

print('... finished formatting and saving rest genepair data')
print('=== finished! ===')

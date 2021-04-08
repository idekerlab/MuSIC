import argparse
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
from file_utils import *
from df_utils import *

def format_densenet_data(prefix, batch, fold, emdfile, workdir):

    outputdir = '{}/train_test_data'.format(workdir)
    sdir = '{}/train_test_data'.format(workdir)
    simdir = '{}/densenet_emb_sim/fold_{}'.format(workdir, fold)

    X_train_gp = np.load('{}/X_train_genepair_{}.npy'.format(sdir, batch))
    X_test_gp = np.load('{}/X_test_genepair_{}.npy'.format(sdir, batch))

    cosine = load_obj('{}/{}.cosine.scaled.pkl'.format(simdir, prefix))
    pearson = load_obj('{}/{}.pearson.scaled.pkl'.format(simdir, prefix))
    spearman = load_obj('{}/{}.spearman.scaled.pkl'.format(simdir, prefix))
    kendall = load_obj('{}/{}.kendall.scaled.pkl'.format(simdir, prefix))
    manhattan = load_obj('{}/{}.manhattan.scaled.pkl'.format(simdir, prefix))
    euclidean = load_obj('{}/{}.euclidean.scaled.pkl'.format(simdir, prefix))
    imgID_to_gname = load_obj('{}/imgID_to_gname.dict.pkl'.format(simdir))
    imgID = [x for x in imgID_to_gname]
    gname = [imgID_to_gname[x] for x in imgID_to_gname]

    emd = load_obj(emdfile).loc[imgID]
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


def format_dn_restGP_data(fold, rest_gp_input, prefix, workdir, emdfile):

    outputdir = '{}/train_test_data'.format(workdir)
    simdir = '{}/densenet_emb_sim/fold_{}'.format(workdir, fold)

    rest_gp = np.load(rest_gp_input)

    cosine = load_obj('{}/{}.cosine.scaled.pkl'.format(simdir, prefix))
    pearson = load_obj('{}/{}.pearson.scaled.pkl'.format(simdir, prefix))
    spearman = load_obj('{}/{}.spearman.scaled.pkl'.format(simdir, prefix))
    kendall = load_obj('{}/{}.kendall.scaled.pkl'.format(simdir, prefix))
    manhattan = load_obj('{}/{}.manhattan.scaled.pkl'.format(simdir, prefix))
    euclidean = load_obj('{}/{}.euclidean.scaled.pkl'.format(simdir, prefix))
    imgID_to_gname = load_obj('{}/imgID_to_gname.dict.pkl'.format(simdir))
    imgID = [x for x in imgID_to_gname]
    gname = [imgID_to_gname[x] for x in imgID_to_gname]

    emd = load_obj(emdfile).loc[imgID]
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

import argparse
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
from file_utils import *
from df_utils import *

def format_n2v_data(prefix, batch, emdfile, workdir):

    outputdir = '{}/train_test_data'.format(workdir)
    sdir = '{}/train_test_data'.format(workdir)
    simdir = '{}/node2vec_emb_sim'.format(workdir)

    X_train_gp = np.load('{}/X_train_genepair_{}.npy'.format(sdir, batch))
    X_test_gp = np.load('{}/X_test_genepair_{}.npy'.format(sdir, batch))

    cosine = load_obj('{}/{}.cosine.scaled.pkl'.format(simdir, prefix))
    pearson = load_obj('{}/{}.pearson.scaled.pkl'.format(simdir, prefix))
    spearman = load_obj('{}/{}.spearman.scaled.pkl'.format(simdir, prefix))
    kendall = load_obj('{}/{}.kendall.scaled.pkl'.format(simdir, prefix))
    manhattan = load_obj('{}/{}.manhattan.scaled.pkl'.format(simdir, prefix))
    euclidean = load_obj('{}/{}.euclidean.scaled.pkl'.format(simdir, prefix))

    emd = load_obj(emdfile)

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
    np.save('{}/{}_X_train_{}.by_gp.npy'.format(outputdir, prefix, batch), X_train)
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
    np.save('{}/{}_X_test_{}.npy'.format(outputdir, prefix, batch), X_test)
    print('... finished formatting and saving testing data')
    print('=== finished! ===')


def format_n2v_restGP_data(rest_gp_input, prefix, workdir, emdfile):

    if not os.path.exists('{}/train_test_data'.format(workdir)):
        os.mkdir('{}/train_test_data'.format(workdir))
    outputdir = '{}/train_test_data'.format(workdir)

    rest_gp = np.load(rest_gp_input)
    simdir = workdir + '/node2vec_emb_sim'

    cosine = load_obj('{}/{}.cosine.scaled.pkl'.format(simdir, prefix))
    pearson = load_obj('{}/{}.pearson.scaled.pkl'.format(simdir, prefix))
    spearman = load_obj('{}/{}.spearman.scaled.pkl'.format(simdir, prefix))
    kendall = load_obj('{}/{}.kendall.scaled.pkl'.format(simdir, prefix))
    manhattan = load_obj('{}/{}.manhattan.scaled.pkl'.format(simdir, prefix))
    euclidean = load_obj('{}/{}.euclidean.scaled.pkl'.format(simdir, prefix))

    emd = load_obj(emdfile)

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
    np.save('{}/rest_genepair_data.npy'.format(outputdir), gp_data)

    print('... finished formatting and saving rest data')
    print('=== finished! ===')

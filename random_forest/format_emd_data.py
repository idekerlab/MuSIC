import pandas as pd
import numpy as np
import os
from random import randint
import random
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from music_utils import *

def format_emd(emd_file, emd_label, proteins, outprefix, num_set='auto'):
    '''
    Format embedding file into different training sets and calculate pairwise protein similarity for each set.
    
    Args:
        emd_file: path to the embedding file. 
                  Format requirement - comma separated file. 
                                       First column is index; 
                                       second column is protein name; 
                                       following columns are embedding values.
                                       Please refer to /Examples/IF_image_embedding.csv as guideline.
        emd_label: label of the embedding file. E.g., IF_emd, APMS_emd
        num_set: {'auto', 1, 2, 3, ...} number of training sets to generate. 
                 'auto': num_set set to the maximum number of embeddings a protein has in the emd_file.
        outprefix: full path to the folder where results will be saved in with unique file identifier. 
    '''
    emd = pd.read_csv(emd_file, index_col=0, header=None)

    # Filter for proteins of interest
    emd = emd[emd[1].isin(proteins)]

    # Determine number of training sets to generate
    count_num_emd = emd[1].value_counts()
    print('Number of embeddings per protein ranges from {} to {}.'.format(count_num_emd.min(), count_num_emd.max()))
    if num_set == 'auto':
        num_set = count_num_emd.max()
    else:
        num_set = int(num_set)
    print('Format {} into {} training sets...'.format(emd_label, num_set))
    
    # NOTE:
    # Original MuSIC study hard-coded to ensure each embedding is used at least once.
    # In consideration of generalizability, we have removed this constraint.
    if not os.path.exists('{}_{}'.format(outprefix, emd_label)):
        os.mkdir('{}_{}'.format(outprefix, emd_label))
    for i in range(1, num_set+1):
        sample_emd_idx = []
        for g in proteins:
            tmp_emdid = list(emd[emd[1] == g].index.values)
            gsample = tmp_emdid*int(np.ceil(num_set/len(tmp_emdid)))
            random.shuffle(gsample)
            sample_emd_idx.append(gsample[:num_set])
    for i in range(1, num_set+1):
        print('Start calculating pairwise protein similarity for {} with training set {} ...'.format(emd_label, i))
        outdir = '{}_{}/training_set_{}'.format(outprefix, emd_label, i)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        set_emdid = [x[i-1] for x in sample_emd_idx]
        perSet_gname_to_emdID = dict(zip(proteins, set_emdid))
        perSet_emdID_to_gname = dict(zip(set_emdid, proteins))
        save_obj(perSet_gname_to_emdID, '{}/gname_to_emdID.dict.pkl'.format(outdir))
        save_obj(perSet_emdID_to_gname, '{}/emdID_to_gname.dict.pkl'.format(outdir))

        # Compute protein pairwise similarity: cosine, Manhattan, Pearson, Spearman, Kendall, Euclidean
        df = emd.loc[set_emdid]
        df.set_index(1, inplace=True)
        cos_sim = cosine_similarity_scaled(df)
        save_obj(cos_sim, '{}/cosine.scaled.pkl'.format(outdir))
        mht_sim = manhattan_similarity(df)
        save_obj(mht_sim, '{}/manhattan.scaled.pkl'.format(outdir))
        euc_sim = euclidean_similarity(df)
        save_obj(euc_sim, '{}/euclidean.scaled.pkl'.format(outdir))
        pcorr = pearson_scaled(df)
        save_obj(pcorr, '{}/pearson.scaled.pkl'.format(outdir))
        scorr = spearman_scaled(df)
        save_obj(scorr, '{}/spearman.scaled.pkl'.format(outdir))
        kcorr = kendall_scaled(df)
        save_obj(kcorr, '{}/kendall.scaled.pkl'.format(outdir))
        
    # Assess distinctness among different training sets generated
    if num_set > 1:
        print('Assess distinctness among different training sets generated (Jaccard)...')
        index = ['train_set_{}'.format(x) for x in range(1, num_set+1)]
        ji_df = pd.DataFrame(1, index=index, columns=index, dtype=float)
        for i in range(1, num_set+1):
            outdir_i = '{}_{}/training_set_{}'.format(outprefix, emd_label, i)
            emdid_set_i = set(list(load_obj('{}/emdID_to_gname.dict.pkl'.format(outdir_i)).keys()))
            for j in range(1, num_set+1):
                outdir_j = '{}_{}/training_set_{}'.format(outprefix, emd_label, j)
                emdid_set_j = set(list(load_obj('{}/emdID_to_gname.dict.pkl'.format(outdir_j)).keys()))
                tmp_ji = jaccard(emdid_set_i, emdid_set_j)
                ji_df.at['train_set_{}'.format(i), 'train_set_{}'.format(j)] = tmp_ji
                ji_df.at['train_set_{}'.format(j), 'train_set_{}'.format(i)] = tmp_ji
        ji_df.to_csv('{}_{}/training_set_distinctness.txt'.format(outprefix, emd_label), sep='\t')
        sns.heatmap(ji_df)
        plt.savefig('{}_{}/training_set_distinctness.png'.format(outprefix, emd_label), 
                    format='png', transparent=True, bbox_inches='tight')
    return

def get_emd_X(outprefix, fold, emdfile, emd_label, train_set, rest_gp):
    ''' 
    Generate features input into random forest regressor.
    
    Args:
        outprefix: prefix of output files
        fold: fold number of k-fold cross validation
        emdfile: path to the file containing all embedding features
        emd_label: unique label for the embedding file
        train_set: the number of training set to use for feature formatting
        rest_gp: gene pairs that lack specific annotations in Gene Ontology
    '''
    sample_dir = '{}_train_test_data'.format(outprefix)
    X_train_gp = np.load('{}/X_train_genepair_{}.npy'.format(sample_dir, fold), allow_pickle=True)
    X_test_gp = np.load('{}/X_test_genepair_{}.npy'.format(sample_dir, fold), allow_pickle=True)

    workdir = '{}_{}/training_set_{}'.format(outprefix, emd_label, train_set)
    cosine = load_obj('{}/cosine.scaled.pkl'.format(workdir))
    pearson = load_obj('{}/pearson.scaled.pkl'.format(workdir))
    spearman = load_obj('{}/spearman.scaled.pkl'.format(workdir))
    kendall = load_obj('{}/kendall.scaled.pkl'.format(workdir))
    manhattan = load_obj('{}/manhattan.scaled.pkl'.format(workdir))
    euclidean = load_obj('{}/euclidean.scaled.pkl'.format(workdir))
    emdID_to_gname = load_obj('{}/emdID_to_gname.dict.pkl'.format(workdir))
    emdID = [x for x in emdID_to_gname]
    gname = [emdID_to_gname[x] for x in emdID_to_gname]

    emd = pd.read_csv(emdfile, index_col=0, header=None).loc[emdID]
    emd.set_index(1, inplace=True)

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
    np.save('{}/X_train_{}.by_gp.npy'.format(workdir, fold), X_train, allow_pickle=True)
    print('... finished formatting and saving training data for {} (fold {}, train set {})'.format(emd_label, 
                                                                                                   fold, train_set))

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
    np.save('{}/X_test_{}.npy'.format(workdir, fold), X_test, allow_pickle=True)
    print('... finished formatting and saving testing data for {} (fold {}, train set {})'.format(emd_label, 
                                                                                                  fold, train_set))
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
    np.save('{}/rest_genepair.npy'.format(workdir), gp_data, allow_pickle=True)
    print('... finished formatting and saving rest genepair data (train set {})'.format(train_set))

    return
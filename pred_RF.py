from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
import argparse
from file_utils import *

parser = argparse.ArgumentParser(description='Arguments for Random Forest Predictor')
parser.add_argument('--modelfname_prefix', help='Prefix filepath to random forest models')
parser.add_argument('--sdir', help='directory for X train and test data')
parser.add_argument('--outputdir', help='output directory')
parser.add_argument('--rest_gp_input', help='Full path to file with gene pairs not in Gene Ontology')
args = parser.parse_args()


for batch in range(1, 6):
    for fold in range(1,7):
        modelfname = args.modelfname_prefix + '_batch_' + str(batch) + '.fold_' + str(fold) + '.pkl'
        pred_RF(args.modelfname, args.sdir, args.outputdir)
        pred_RF_restGP(args.modelfname, args.sdir, args.outputdir, args.rest_gp_input)


#concatenate batches for each fold, then average across folds
pred_dir = '{}/yPred_output'.format(workdir)
model_prefix = model_fname.split('/')[-1].rstrip('.pkl')
for fold in range(1, 7):
    ypred_all = np.asarray([])
    for batch in range(1, 6):
        ypred = np.load('{}/yPred_{}_batch_{}.fold_{}.npy'.format(pred_dir, model_prefix, batch, fold))
        ypred_all = np.concatenate((ypred_all, ypred))
    np.save('{}/allyPred_{}.fold_{}.npy'.format(pred_dir, model_prefix, fold), ypred_all)

pred_fnames = glob('{}/allyPred_{}.*'.format(pred_dir, model_prefix))
allpred = []
for file in pred_fnames:
    allpred.append(list(np.load(file)))
allpred = np.asarray(allpred)
print(allpred.shape)
avg_allpred = allpred.mean(axis=0)
print(avg_allpred.shape)
np.save('{}/avg_allyPred_{}.npy'.format(pred_dir, model_prefix),
        avg_allpred)

restGP_fnames = glob('{}/restGP_{}*'.format(pred_dir, model_prefix))
restGP = []
for file in restGP_fnames:
    restGP.append(list(np.load(file)))
restGP = np.asarray(restGP)
print(restGP.shape)
avg_restGP = restGP.mean(axis=0)
print(avg_restGP.shape)
np.save('{}/avg_restGP_{}.npy'.format(pred_dir, model_prefix),
        avg_restGP)



def pred_RF(model_fname, sdir, outdir):
    model = load_obj(model_fname)
    batch = model_fname.split('_batch_')[-1].split('.')[0]
    if 'fold' in model_fname:
        fold = int(model_fname.split('fold_')[-1].rstrip('.pkl'))
    else:
        fold = -1
    img_label = model_fname.split('/')[-1].split('.')[0]
    ppi_label = model_fname.split('/')[-1].split('.')[1]
    print('... test on batch {} and fold {} for {} (image) and {} (PPI)'.format(batch, fold, img_label, ppi_label))

    # Load Image testing data
    if fold == -1:
        img_X_test = np.load('{}/X_test_{}.npy'.format(sdir, batch))
    else:
        img_X_test = np.load('{}/X_test_{}.fold_{}.npy'.format(sdir, batch, fold))
    # Load PPI testing data
    ppi_X_test = np.load('{}/{}_X_test_{}.npy'.format(sdir, ppi_label, batch))
    # Concatenate Image and PPI data to generate testing data
    X_test = np.concatenate((img_X_test, ppi_X_test), axis=1)
    print('... loaded test data')

    y_pred = model.predict(X_test)

    np.save('{}/yPred_{}.npy'.format(outdir, model_fname.split('/')[-1].rstrip('.pkl')), y_pred)
    print('=== finished! ===')

def pred_RF_restGP(model_fname, sdir, outdir, rest_gp_input):
    model = load_obj(model_fname)
    fold = int(model_fname.split('fold_')[-1].replace('.pkl', ''))
    img_label = model_fname.split('/')[-1].split('.')[0]
    ppi_label = model_fname.split('/')[-1].split('.')[1]
    print('... test on fold {} for {} (image) and {} (PPI)'.format(fold, img_label, ppi_label))

    # Load Image testing data
    img_X_test = np.load('{}/rest_genepair.fold_{}.npy'.format(sdir, fold))

    # Load PPI testing data
    ppi_X_test = np.load(rest_gp_input)
    # Concatenate Image and PPI data to generate testing data
    X_test = np.concatenate((img_X_test, ppi_X_test), axis=1)
    print('... loaded data of rest_genepair')

    y_pred = model.predict(X_test)

    np.save('{}/restGP_{}.npy'.format(outdir, model_fname.split('/')[-1].rstrip('.pkl')), y_pred)
    print('=== finished! ===')
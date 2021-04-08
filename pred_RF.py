from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import numpy as np
import sys
import os
from glob import glob
from tqdm import tqdm
import argparse
from file_utils import *



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

def pred_RF_restGP(model_fname, sdir, outdir):
    model = load_obj(model_fname)
    fold = int(model_fname.split('fold_')[-1].replace('.pkl', ''))
    img_label = model_fname.split('/')[-1].split('.')[0]
    ppi_label = model_fname.split('/')[-1].split('.')[1]
    print('... test on fold {} for {} (image) and {} (PPI)'.format(fold, img_label, ppi_label))

    # Load Image testing data
    img_X_test = np.load('{}/rest_genepair.fold_{}.npy'.format(sdir, fold))

    # Load PPI testing data
    ppi_X_test = np.load('{}/rest_genepair_data.npy'.format(sdir))
    # Concatenate Image and PPI data to generate testing data
    X_test = np.concatenate((img_X_test, ppi_X_test), axis=1)
    print('... loaded data of rest_genepair')

    y_pred = model.predict(X_test)

    np.save('{}/restGP_{}.npy'.format(outdir, model_fname.split('/')[-1].rstrip('.pkl')), y_pred)
    print('=== finished! ===')


parser = argparse.ArgumentParser(description='Arguments for Random Forest Predictor')
parser.add_argument('--modelfname_prefix', help='Prefix filepath to random forest models')
parser.add_argument('--sdir', help='directory for X train and test data')
parser.add_argument('--outputdir', help='output directory')
parser.add_argument('--rest_gp_input', help='Full path to file with gene pairs not in Gene Ontology (Example file rest_genepair.npy in Examples folder)')
args = parser.parse_args()


for batch in range(1, 6):
    for fold in range(1,7):
        modelfname = args.modelfname_prefix + '_batch_' + str(batch) + '.fold_' + str(fold) + '.pkl'
        pred_RF(modelfname, args.sdir, args.outputdir)
        pred_RF_restGP(modelfname, args.sdir, args.outputdir)


#concatenate batches for each fold, then average across folds
model_prefix = args.modelfname_prefix.split('/')[-1].rstrip('.pkl')
for fold in range(1, 7):
    ypred_all = np.asarray([])
    for batch in range(1, 6):
        ypred = np.load('{}/yPred_{}_batch_{}.fold_{}.npy'.format(args.outputdir, model_prefix, batch, fold))
        ypred_all = np.concatenate((ypred_all, ypred))
    np.save('{}/allyPred_{}.fold_{}.npy'.format(args.outputdir, model_prefix, fold), ypred_all)

pred_fnames = glob('{}/allyPred_{}.*'.format(args.outputdir, model_prefix))
allpred = []
for file in pred_fnames:
    allpred.append(list(np.load(file)))
allpred = np.asarray(allpred)
avg_allpred = allpred.mean(axis=0)
np.save('{}/avg_allyPred_{}.npy'.format(args.outputdir, model_prefix),
        avg_allpred)

restGP_fnames = glob('{}/restGP_{}*'.format(args.outputdir, model_prefix))
restGP = []
for file in restGP_fnames:
    restGP.append(list(np.load(file)))
restGP = np.asarray(restGP)
avg_restGP = restGP.mean(axis=0)
np.save('{}/avg_restGP_{}.npy'.format(args.outputdir, model_prefix),
        avg_restGP)

avg_ypred = np.load('{}/avg_allyPred_{}.npy'.format(args.outputdir, model_prefix))
avg_restgp = np.load('{}/avg_restGP_{}.npy'.format(args.outputdir, model_prefix))
gp = np.empty(0)
for batch in range(1, 6):
    gp_bybatch = np.load('{}/X_test_genepair_{}.npy'.format(args.sdir, batch))
    if batch == 1:
        gp = gp_bybatch
    else:
        gp = np.concatenate((gp, gp_bybatch))
restgp = np.load(args.rest_gp_input)

gp_all = np.concatenate([gp, restgp])
pred_all = np.concatenate([avg_ypred, avg_restgp])

df = pd.DataFrame(gp_all)
df.columns = ['geneA', 'geneB']
df['weight'] = pred_all

save_obj(df, '{}/predicted_resnik_sim.pkl'.format(args.outputdir))
df.to_csv('{}/predicted_resnik_sim.ddot'.format(args.outputdir), sep='\t', header=False, index=False)
import argparse
import pandas as pd
import numpy as np
import sys
import os
from file_utils import *
from df_utils import *

def format_GO_data(go_input, workdir, gene_list):

    # Load CC GO without HPA evidence
    go = pd.read_table(go_input,
                       header=None, index_col=0)
    go.columns = ['tsize', 'genes']

    # Load 661 HEK293 genes
    hek_genes = np.load(gene_list)

    # Remove genes directly assigned to root (suggests not well known genes)
    root = 'GO:0005575'
    root_genes = go.loc[root]['genes'].split(',')[:-1]
    nonRootGenes = []
    for idx, row in go.iterrows():
        if idx == root:
            continue
        nonRootGenes += row['genes'].split(',')[:-1]
    nonRootGenes = set(nonRootGenes)

    hek_x_cc = set(hek_genes).intersection(set(nonRootGenes))
    hek_x_cc = np.asarray(list(hek_x_cc))
    np.save('{}/hek293_x_cc_genes.602.npy'.format(workdir), hek_x_cc)

    # Sort GO by term size !!!
    go.sort_values('tsize', inplace=True)
    go.head()

    hek_resnik = pd.DataFrame(-1, index=hek_x_cc, columns=hek_x_cc, dtype=float)
    hek_x_cc_geneset = set(hek_x_cc)
    count = 0  # for process tracking
    N = len(nonRootGenes)
    for idx, row in go.iterrows():
        count += 1
        if count % 10 == 0:
            print('... processed {}%'.format(count * 100 / go.shape[0]))
        rsim = -np.log10(row['tsize'] / N)
        genes = set(row['genes'].split(',')[:-1]).intersection(hek_x_cc_geneset)
        if len(genes) > 1:
            genes = list(genes)
            for i in range(len(genes) - 1):
                ga = genes[i]
                for j in range(i + 1, len(genes)):
                    gb = genes[j]
                    if hek_resnik.at[ga, gb] == -1:
                        hek_resnik.at[ga, gb] = rsim
                        hek_resnik.at[gb, ga] = rsim

    # Assign genes never appeared in the same term similarity of 0
    hek_resnik.replace(-1, 0, inplace=True)
    # Fill in the diagonal with maximum resnik value
    max_resnik = hek_resnik.max().max()
    np.fill_diagonal(hek_resnik.values, max_resnik)
    # Scale resnik into [0, 1]
    hek_resnik_scaled = hek_resnik / max_resnik
    # Save scaled file
    with open('{}/hek_ccNoHPA_resnik.602.scaled.pkl'.format(workdir), 'wb') as f:
        pickle.dump(hek_resnik_scaled, f)


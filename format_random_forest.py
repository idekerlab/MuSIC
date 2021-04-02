import os
import pickle
from random import randint
from df_utils import *
from file_utils import *
from format_GO_data import *
from format_dn_data import *
from format_n2v_data import *
from format_x_and_y_data import *
import random
import argparse
from sklearn.model_selection import train_test_split, KFold
from numpy.random import choice

parser = argparse.ArgumentParser(description='Format files for random forest analysis.')
parser.add_argument('--go_input', help='Full path to file containing list of Gene Ontology terms size and genes.')
parser.add_argument ('--gene_list', help='Full path to list of genes to analyze')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--n2v_prefix', help='node2vec embedding prefix E.g. dim1024_p_2_q_1')
parser.add_argument('--densenet_prefix', help='densenet embedding prefix E.g. densenet_raw_1024D')
parser.add_argument('--n2v_emdfile', help='Full path to node2vec feature embedding file')
parser.add_argument('--densenet_emdfile', help='Full path to densenet feature embedding file')
parser.add_argument('--rest_gp_input', help='Full path to file with gene pairs not in the Gene Ontology')

args = parser.parse_args()

workdir = args.outprefix

format_GO_data(args.go_input, workdir, args.gene_list)
format_x_and_y_data(workdir)

for batch in range(1, 6):
    format_n2v_data(args.n2v_prefix, batch, args.sdir, workdir)
format_n2v_restGP_data(args.rest_gp_input, args.simdir, args.n2v_prefix, workdir, args.n2v_emdfile)

for batch in range(1, 6):
    for fold in range(1, 7):
        format_densenet_data(args.densenet_prefix, batch, fold, args.densenet_emdfile, workdir)
for fold in range(1, 7):
    format_n2v_restGP_data(args.rest_gp_input, args.densenet_prefix, workdir, args.densenet_emdfile)

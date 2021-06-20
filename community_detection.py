import pandas as pd
import numpy as np
import sys
import os
import argparse
from music_utils import *


parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--path_to_clixo', help='Full path to CliXO folder.')
parser.add_argument('--clixo_i', help='Path to input similarity network for CliXO.')
parser.add_argument('--clixo_a', type=float, default=0.1, help='CliXO parameter alpha.')
parser.add_argument('--clixo_b', type=float, default=0.5, help='CliXO parameter beta.')
parser.add_argument('--clixo_m', type=float, default=0, help='CliXO parameter m.')
parser.add_argument('--clixo_z', type=float, default=0, help='CliXO parameter z.')
parser.add_argument('--clixo_s', type=float, default=0, help='CliXO parameter s.')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology folder.')
parser.add_argument('--minSystemSize', type=int, default=4, 
                    help='Minimum number of proteins requiring each system to have.')
parser.add_argument('--ci_thre', type=float, default=0.75, help='Containment index threshold')
parser.add_argument('--ji_thre', type=float, default=0.9, 
                    help='Jaccard index threshold for merging similar clusters')
parser.add_argument('--niter', type=int, default=1000, help='Number of iterations Louvain clustering will run to select partition with the best modularity.')
parser.add_argument('--min_diff', type=int, default=1, help='Minimum difference in number of proteins for every parent-child pair.')
parser.add_argument('--keep_all_files', help='Keep all intermediate output.', action='store_true')
parser.add_argument('--predict_nm_size', help='Predict size in nm for each system.', action='store_true')
parser.add_argument('--n_samples', default=1000000, type=int,
                    help='Maximum number of samples to use for fitting linear model.')
parser.add_argument('--q_step', type=float, default=0.1, help='Step for scanning best quantile to use.')
args = parser.parse_args()

outprefix = args.outprefix

# Check CliXO
if args.path_to_clixo is None:
    path_to_clixo = os.path.dirname(os.path.realpath(__file__)) + '/CliXO-1.0/clixo'
    if not os.path.exists(path_to_clixo):
        raise ValueError('Please first install CliXO following installation instructions.')
else:
    path_to_clixo = args.path_to_clixo.rstrip('/') + '/clixo'
# Check if clixo path is valid
if not os.access(path_to_clixo, os.X_OK):
    raise ValueError('{} is not executable!'.format(path_to_clixo))

# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)
# Check if input to clixo exists
if not os.path.exists(args.clixo_i):
    raise ValueError('{} does not exist!'.format(args.clixo_i))

# Check alignOntology
if args.path_to_alignOntology is None:
    path_to_alignOntology = os.path.dirname(os.path.realpath(__file__)) + '/ddot/ddot/alignOntology'
else:
    path_to_alignOntology = args.path_to_alignOntology.rstrip('/')
if not os.path.exists(path_to_alignOntology):
    raise ValueError('Please first install DDOT following installation instructions.')
# Check if ontologyTermStats is valid
if not os.access('{}/ontologyTermStats'.format(path_to_alignOntology), os.X_OK):
    raise ValueError('{}/ontologyTermStats not executable!'.format(path_to_alignOntology))

print('CliXO parameters: a={}, b={}, m={}, z={}, s={}'.format(args.clixo_a, args.clixo_b, args.clixo_m, args.clixo_z, args.clixo_s))
cdDir = os.path.dirname(os.path.realpath(__file__))

with open('{}.sh'.format(outprefix), 'w') as scriptfile:
    scriptfile.write('#!/bin/bash\n')
    scriptfile.write('{} -i {} -a {} -b {} -m {} -z {} -s {} > {}.clixoOutFile\n'.format(path_to_clixo, 
                                                                                         args.clixo_i, 
                                                                                         args.clixo_a, 
                                                                                         args.clixo_b, 
                                                                                         args.clixo_m, 
                                                                                         args.clixo_z,
                                                                                         args.clixo_s,
                                                                                         outprefix))
    scriptfile.write('python {}/community_detection/process_clixoOutFile.py --outprefix {} --path_to_alignOntology {}\n'.format(cdDir, outprefix, path_to_alignOntology))
    scriptfile.write('python {}/community_detection/trim_clixo.py --outprefix {} --minSystemSize {}\n'.format(cdDir, outprefix, args.minSystemSize))
    scriptfile.write('python {}/community_detection/mature_hier_structure.py --outprefix {} --ci_thre {} --ji_thre {} --minSystemSize {} --path_to_alignOntology {} --min_diff {}\n'.format(cdDir, outprefix, args.ci_thre, args.ji_thre, args.minSystemSize, path_to_alignOntology, args.min_diff))
    scriptfile.write('python {}/community_detection/louvain_partition.py --outprefix {} --clixo_i {} --niter {}\n'.format(cdDir, outprefix, args.clixo_i, args.niter))
    scriptfile.write('python {}/community_detection/cap_louvain.py --outprefix {} --path_to_alignOntology {} --minSystemSize {} --niter {}\n'.format(cdDir, outprefix, path_to_alignOntology, args.minSystemSize, args.niter))
    if args.predict_nm_size:
        scriptfile.write('python {}/community_detection/predict_size.py --outprefix {} --path_to_alignOntology {} --n_samples {} --q_step {}\n'.format(cdDir, outprefix, path_to_alignOntology, args.n_samples, args.q_step))
    if not args.keep_all_files:
        scriptfile.write('python {}/community_detection/clean_up.py --outprefix {}\n'.format(cdDir, outprefix))


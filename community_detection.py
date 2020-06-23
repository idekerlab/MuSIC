import pandas as pd
import numpy as np
import sys
import os
import argparse
from file_utils import *


parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--path_to_clixo', help='Full path to CliXO folder.')
parser.add_argument('--clixo_i', help='Path to input similarity network for CliXO.')
parser.add_argument('--clixo_a', type=float, default=0.1, help='CliXO parameter alpha.')
parser.add_argument('--clixo_b', type=float, default=0.5, help='CliXO parameter beta.')
parser.add_argument('--clixo_m', type=float, default=0, help='CliXO parameter m.')
parser.add_argument('--clixo_z', type=float, default=0, help='CliXO parameter z.')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology folder.')
parser.add_argument('--minSystemSize', type=int, default=2, 
                    help='Minimum number of proteins requiring each system to have.')
parser.add_argument('--ci_thre', type=float, default=0.75, help='Containment index threshold')
parser.add_argument('--ji_thre', type=float, default=0.9, 
                    help='Jaccard index threshold for merging similar clusters')
args = parser.parse_args()

outprefix = args.outprefix
path_to_clixo = args.path_to_clixo.rstrip('/') + '/clixo'
cdDir = os.getcwd()

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

# Check if ontologyTermStats is valid
if not os.access('{}/ontologyTermStats'.format(args.path_to_alignOntology.rstrip('/')), os.X_OK):
    raise ValueError('{}/ontologyTermStats not executable!'.format(args.path_to_alignOntology.rstrip('/')))

print('CliXO parameters: a={}, b={}, m={}, z={}'.format(args.clixo_a, args.clixo_b, args.clixo_m, args.clixo_z))



with open('{}.sh'.format(outprefix), 'w') as scriptfile:
    scriptfile.write('#!/bin/bash\n')
    scriptfile.write('{} -i {} -a {} -b {} -m {} -z {} > {}.clixoOutFile\n'.format(path_to_clixo, 
                                                                                   args.clixo_i, 
                                                                                   args.clixo_a, 
                                                                                   args.clixo_b, 
                                                                                   args.clixo_m, 
                                                                                   args.clixo_z, 
                                                                                   outprefix))
    scriptfile.write('python {}/community_detection/process_clixoOutFile.py --outprefix {} --path_to_alignOntology {}\n'.format(cdDir, outprefix, args.path_to_alignOntology))
    scriptfile.write('python {}/community_detection/trim_clixo.py --outprefix {} --minSystemSize {}\n'.format(cdDir, outprefix, args.minSystemSize))
    scriptfile.write('python {}/community_detection/mature_hier_structure.py --outprefix {} --ci_thre {} --ji_thre {} --minSystemSize {} --path_to_alignOntology {}\n'.format(cdDir, outprefix, args.ci_thre, args.ji_thre, args.minSystemSize, args.path_to_alignOntology))
    scriptfile.write('python {}/community_detection/louvain_partition.py --outprefix {} --clixo_i {}\n'.format(cdDir, outprefix, args.clixo_i))
    scriptfile.write('python {}/community_detection/cap_louvain.py --outprefix {} --path_to_alignOntology {} --minSystemSize {}\n'.format(cdDir, outprefix, args.path_to_alignOntology, args.minSystemSize))


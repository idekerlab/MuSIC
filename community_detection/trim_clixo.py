import pandas as pd
import numpy as np
import sys
import os
import ddot
import argparse
from ddot import Ontology

parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--minSystemSize', type=int, default=4, 
                    help='Minimum number of proteins requiring each system to have.')
args = parser.parse_args()

minSystemSize = args.minSystemSize
print('... start trimming clixo to keep only systems with at least {} proteins'.format(minSystemSize))

f = '{}.clixoOutFile.oaInput'.format(args.outprefix)
ont = Ontology.from_table(f, clixo_format=True)
tsize = list(zip(ont.terms, ont.term_sizes))
delete_term = [x[0] for x in tsize if x[1] < minSystemSize]
trimmed_ont = ont.delete(to_delete=delete_term)
trimmed_ont = trimmed_ont.collapse_ontology(method='python')
trimmed_ont.to_table(f+'.mintsize{}'.format(minSystemSize), clixo_format=True)
print('=== finished trim_clixo.py ===')
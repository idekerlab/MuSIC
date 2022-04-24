import pandas as pd
import numpy as np
import sys
import os
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.stats import hypergeom
from music_utils import *

parser = argparse.ArgumentParser(description='Analyze each system in the given hierarchy.')
parser.add_argument('--hier_fname', help='Full path of the input termStats file.')
parser.add_argument('--genes', help='Path to the file storing all genes in the hierarchy. One gene per line.')
parser.add_argument('--outprefix', help='Prefix of output file path. Format: path_to_outdir/file_identifier')
parser.add_argument('--w_root', action='store_true', help='Do analysis for root term.')
parser.add_argument('--eval_matrix', help='Path to evaluation matrix file.')
parser.add_argument('--minTermSize', default=3, type=int, 
                    help='Do analysis for all terms with term size >= minTermSize')
args = parser.parse_args()

f = args.hier_fname
minTermSize = args.minTermSize
outprefix = args.outprefix
hiergenes = [x.rstrip() for x in open(args.genes).readlines()]
root_size = len(hiergenes)

# read in hierarchy termStats file and perform sanity check
if not os.path.exists(f):
    raise ValueError('Input termStats file does not exist!')
if os.path.getsize(f) == 0:
    print('=== No term in hierarchy! ===')
    sys.exit()
df = pd.read_table(f, header=None)
df.columns = ['term', 'tsize', 'genes']
if args.w_root:
    df = df[df['tsize'] >= minTermSize]
else:
    df = df[(df['tsize'] >= minTermSize) & (df['tsize'] < root_size)]
if df.shape[0] == 0:
    print('=== No system left after size filter ===')
    sys.exit()
df.set_index('term', inplace=True)


# Load evaluation edge matrix and perform sanity check
eval_m = pd.read_table(args.eval_matrix, sep=',', index_col=0)
if set(eval_m.values.flatten()) != set([0, 1]):
    raise ValueError('Please only have binary values (0 and 1) in evaluation matrix.')
if not check_symmetric(eval_m):
    raise ValueError('Please ensure evaluation matrix is symmetric.')
if np.any(eval_m.index != eval_m.columns):
    raise ValueError('Index and column of evalation matrix do not agree.')

eval_genes = np.asarray(list(set(eval_m.index.values).intersection(set(hiergenes))))
print('{} genes in hierarchy has observation in evaluation dataset.'.format(len(eval_genes))) 
eval_total = upper_tri_values(eval_m).sum()
eval_M = len(upper_tri_values(eval_m))


eval_count = []
eval_pvalue = [] # Note: calculate pvalue based on genes with observation only

for idx, row in df.iterrows():
    termgenes = row['genes'].split(',')[:-1]
    eval_x_termgenes = list(set(termgenes).intersection(set(eval_genes)))
    if len(eval_x_termgenes) < 2:
        eval_count.append(0)
        eval_pvalue.append(1)
    else:
        eval_count_value = 0
        for i in range(len(eval_x_termgenes)-1):
            ga = eval_x_termgenes[i]
            for j in range(i+1, len(eval_x_termgenes)):
                gb = eval_x_termgenes[j]
                if eval_m.at[ga, gb] == 1:
                    eval_count_value += 1
        eval_count.append(eval_count_value)
        eval_pvalue.append(hypergeom.sf(eval_count_value-1, eval_M, eval_total, num_comb(len(eval_x_termgenes))))

df['eval_count'] = eval_count
df['eval_pvalue'] = eval_pvalue
df['eval_adjp'] = multipletests(eval_pvalue, method='fdr_bh')[1]
df.drop(columns='genes', inplace=True)

root_file_label = 'noRoot'
if args.w_root:
    root_file_label = 'wRoot'
df.to_csv('{}.eval_{}.csv'.format(outprefix, root_file_label), sep='\t')
print('=== finished ===')

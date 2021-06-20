import numpy as np
import pandas as pd
import os
import sys
from tqdm import tqdm
import argparse
from scipy.stats import linregress
from music_utils import *

parser = argparse.ArgumentParser(description='Calibrate protein-protein distance and proximity based on Gene Ontology.')
parser.add_argument('--protein_file', help='Path to the file containing list of proteins to analyze. E.g., /Examples/MuSIC_proteins.txt')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g., /path/to/output/directory/fileIdentifier')
parser.add_argument('--C_file', help='Path to precomputed matrix containing the size (number of proteins) of smallest Gene Ontology component shared by the gene pair.')
args = parser.parse_args()

if args.outprefix is None:
    raise ValueError('Please specify --outprefix: Prefix for files generated. E.g., /path/to/output/directory/fileIdentifier')
outprefix = args.outprefix
# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
# Check if protein_file exists
if not os.path.exists(args.protein_file):
    raise ValueError('File path entered for protein_file does not exist!')
    
# If C_file is specified, check if path entered is valid.
if args.C_file != None:
    if not os.path.exists(args.C_file):
        raise ValueError('File path entered for C_file does not exist!')

music_dir = os.path.dirname(os.path.realpath(__file__))
# Get Calibration Function
calibration_data = pd.read_table('{}/data/calibration.txt'.format(music_dir))
slope, intercept, r_value, _, _ = linregress(np.log10(calibration_data['C']), 
                                             np.log10(calibration_data['D']))
print('Calibration Function:\nlog10(D) = {} * log10(C) + {}'.format(round(slope, 2), round(intercept, 2)))
print('R2 = {}\n'.format(round(r_value**2, 2)))

proteins = [x.rstrip('\n') for x in open(args.protein_file, 'r').readlines()]
print('Analyzing {} proteins:'.format(len(proteins)))
print(', '.join(x for x in proteins))

# Load GO CC without HPA evidence
go = pd.read_table('{}/data/GO_CC_human_no_hpa.txt'.format(music_dir), header=None, index_col=0)
go.columns = ['tsize', 'genes']
go_proteins = [x.rstrip('\n') for x in open('{}/data/GO_annotated_proteins.txt'.format(music_dir), 'r').readlines()]

annot_proteins = set(proteins).intersection(set(go_proteins))
rest_proteins = set(proteins) - annot_proteins
annot_proteins = np.asarray(list(annot_proteins))
rest_proteins = np.asarray(list(rest_proteins))
print('\n{} proteins entered have specific annotation in GO other than root.'.format(len(annot_proteins)))
np.save('{}.annot_proteins.npy'.format(outprefix), annot_proteins)
np.save('{}.not_annot_proteins.npy'.format(outprefix), rest_proteins)

if args.C_file is None:
    print('\nComputing protein pairwise matrix for C...')
    go.sort_values('tsize', inplace=True)
    pairwise_C = pd.DataFrame(-1, index=annot_proteins, columns=annot_proteins, dtype=int)
    for idx in tqdm(go.index):
        row = go.loc[idx]
        glist = list(set(row['genes'].split(',')[:-1]).intersection(set(annot_proteins)))
        for i in range(len(glist)-1):
            ga = glist[i]
            for j in range(i+1, len(glist)):
                gb = glist[j]
                if pairwise_C.at[ga, gb] == -1:
                    pairwise_C.at[ga, gb] = row['tsize']
                    pairwise_C.at[gb, ga] = row['tsize']
        if np.all(upper_tri_values(pairwise_C) != -1):
            break # early termination
    pairwise_C.replace(18535, len(go_proteins), inplace=True)
    pairwise_C.replace(-1, len(go_proteins), inplace=True)
else:
    pairwise_C = pd.read_csv(args.C_file, index_col=0).loc[annot_proteins][annot_proteins]

# Compute D and P from C
print('Computing calibrated protein-protein distance (D, nm) and protein-protein proximity (P) from C...')
geneA = []
geneB = []
C_list = []
for i in range(len(annot_proteins)-1):
    for j in range(i+1, len(annot_proteins)):
        geneA.append(annot_proteins[i])
        geneB.append(annot_proteins[j])
        C_list.append(pairwise_C.at[annot_proteins[i], annot_proteins[j]])
df = pd.DataFrame()
df['geneA'] = geneA
df['geneB'] = geneB
df['C'] = C_list
df['log10D'] = slope*np.log10(df['C']) + intercept
df['D'] = 10**df['log10D']
df['P'] = -df['log10D']
df.to_csv('{}.calibrated_distance.csv'.format(outprefix), index=False)
print('=== finished! ===')
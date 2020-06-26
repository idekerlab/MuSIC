import pandas as pd
import numpy as np
import sys
import os
import igraph as ig
import louvain
import argparse
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from file_utils import *


parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--clixo_i', help='Path to input similarity network for CliXO.')
parser.add_argument('--niter', type=int, default=1000, help='Number of iterations Louvain clustering will run to select partition with the best modularity.')
args = parser.parse_args()

print('... start community detection with louvain algorithm')
outprefix = args.outprefix
niter = args.niter
df = pd.read_table(args.clixo_i, header=None)
df.columns = ['geneA', 'geneB', 'weight']
hiergenes = list(set(df['geneA'].values).union(set(df['geneB'].values)))

G = ig.Graph()
G.add_vertices(hiergenes)
edges = list(zip(df['geneA'].values, df['geneB'].values))
G.add_edges(edges)
G.es['weight'] = df['weight'].values
if G.is_directed():
    raise ValueError('Graph should not be directed!')
if not G.is_weighted():
    raise ValueError('Graph should be weighted!')
# Partition niter times, find the partition with best modularity
modularity_list = []
membership_list = []
for i in range(niter):
    partition = louvain.find_partition(G, louvain.ModularityVertexPartition, weights='weight')
    modularity_list.append(partition.modularity)
    membership_list.append(partition.membership)
modularity_list = np.asarray(modularity_list)
membership_list = np.asarray(membership_list)
# save niter iteration results
np.save('{}.mvp_partition_modularity.{}.npy'.format(outprefix, niter), modularity_list)
np.save('{}.mvp_partition_membership.{}.npy'.format(outprefix, niter), membership_list)
# Get node ids
node_idx = []
node_name = []
for v in G.vs():
    node_idx.append(v.index)
    node_name.append(v['name'])
name_to_idx = dict(zip(node_name, node_idx))
idx_to_name = dict(zip(node_idx, node_name))
save_obj(G, '{}.G.pkl'.format(outprefix))
save_obj(name_to_idx, '{}.node_name_to_idx.dict.pkl'.format(outprefix))
save_obj(idx_to_name, '{}.node_idx_to_name.dict.pkl'.format(outprefix))
# Find the best modularity
best_idx = np.where(modularity_list == modularity_list.max())[0][0]
best_partition = membership_list[best_idx]
unique_cluster = list(set(best_partition))
print('Louvain best modularity partition gave {} clusters.'.format(len(unique_cluster)))
print('=== finished louvain_partition.py ===')
import pandas as pd
import numpy as np
import sys
import os
import argparse
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from music_utils import *


def contain_frac(query_list, target_list):
    query = set(query_list)
    target = set(target_list)
    return len(query.intersection(target)) / len(query)

def get_orphan_genes(ts_path, root_gene_list, minSystemSize=3):
    ts_df = pd.read_table(ts_path, header=None)
    ts_df.columns = ['term', 'tsize', 'genes']
    ts_df = ts_df[(ts_df['tsize'] >= minSystemSize) & (ts_df['tsize'] < len(root_gene_list))]
    non_root_gene = []
    for idx, row in ts_df.iterrows():
        non_root_gene += row['genes'].split(',')[:-1]
    non_root_geneset = set(non_root_gene)
    orphan_genes = list(set(root_gene_list) - non_root_geneset)
    return orphan_genes



parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology folder.')
parser.add_argument('--minSystemSize', type=int, default=2, 
                    help='Minimum number of proteins requiring each system to have.')
parser.add_argument('--niter', type=int, default=1000, help='Number of iterations Louvain clustering will run to select partition with the best modularity.')
args = parser.parse_args()

outprefix = args.outprefix
niter = args.niter

ddotfname = '{}.ddot'.format(outprefix)
# Load processed ddot file
ont_edge = pd.read_table(ddotfname, header=None, dtype=str)
ont_edge.columns = ['parent', 'child', 'type']
# Load termStats file
ont_ts = pd.read_table('{}.termStats'.format(outprefix), header=None, dtype=str)
ont_ts.columns = ['comp', 'tsize', 'genes']
ont_ts.set_index('comp', inplace=True)
# Load louvain outputs
idx_to_name = load_obj('{}.node_idx_to_name.dict.pkl'.format(outprefix))
partition_list = np.load('{}.mvp_partition_membership.{}.npy'.format(outprefix, niter))
modularity_list = np.load('{}.mvp_partition_modularity.{}.npy'.format(outprefix, niter))

hiergenes = list(set(ont_edge[ont_edge['type'] == 'gene']['child'].values))
best_idx = np.where(modularity_list == modularity_list.max())[0][0]
best_partition = partition_list[best_idx]
unique_partition = list(set(best_partition))
louvain_cluster = []
for c in unique_partition:
    cluster_idx = np.where(best_partition == c)[0]
    cluster_gene = [idx_to_name[x] for x in cluster_idx]
    louvain_cluster.append(cluster_gene)

root = ont_ts[ont_ts['tsize'] == str(len(hiergenes))].index[0]
first_layer_term = ont_edge[(ont_edge['parent'] == root) & (ont_edge['type'] == 'default')]['child'].values

df = pd.DataFrame(index=first_layer_term)
for fl in first_layer_term:
    fl_genes = ont_ts.loc[fl]['genes'].split(',')[:-1]
    for i, cluster in enumerate(louvain_cluster):
        df.at[fl, 'Cluster_{}'.format(i)] = contain_frac(fl_genes, cluster)
# Delete all root->node edges
drop_idx = ont_edge[ont_edge['parent'] == root].index.values
ont_edge.drop(index=drop_idx, inplace=True)

contain_thre = 0.5
# Add first layer term to louvain
child_term = []
for c in df.columns:
    c_child_list = []
    for comp, row in df.iterrows():
        if row[c] >= contain_thre:
            c_child_list.append(comp)
    child_term.append(c_child_list)
# If a first layer term did not pass contain_thre, add to louvain with largest contain
for comp, row in df.iterrows():
    if np.any(row.values >= contain_thre):
        continue
    lc_idx_list = np.where(row.values == row.values.max())[0]
    for lc_idx in lc_idx_list:
        child_term[lc_idx].append(comp)
# Add parent child relationships
add_parent = []
add_child = []
add_type = []
# Add root->louvain_cluster edges
for c in df.columns:
    add_parent.append(root)
    add_child.append(c)
    add_type.append('default')
# Add louvain_cluster->first_layer_term edges
for idx, c in enumerate(df.columns):
    if str(idx) != c.split('_')[-1]:
        raise ValueError('index position and column name disagree!')
    c_child_list = child_term[idx]
    for comp in c_child_list:
        add_parent.append(c)
        add_child.append(comp)
        add_type.append('default')

orphan_genes = get_orphan_genes('{}.termStats'.format(outprefix), hiergenes, minSystemSize=args.minSystemSize)
# Add louvain_cluster->gene edges
for idx, c in enumerate(df.columns):
    if str(idx) != c.split('_')[-1]:
        raise ValueError('index position and column name disagree!')
    c_only_gene = list(set(louvain_cluster[idx]).intersection(set(orphan_genes)))
    add_parent += [c]*len(c_only_gene)
    add_child += c_only_gene
    add_type += ['gene']*len(c_only_gene)

add_edge = pd.DataFrame()
add_edge['parent'] = add_parent
add_edge['child'] = add_child
add_edge['type'] = add_type
edge_df = pd.concat([ont_edge, add_edge])
edge_df.to_csv('{}.louvain.ddot'.format(outprefix), header=False, index=False, sep='\t')

cmd = '{}/ontologyTermStats {} genes > {}'.format(args.path_to_alignOntology.rstrip('/'), 
                                                  '{}.louvain.ddot'.format(outprefix), 
                                                  '{}.louvain.termStats'.format(outprefix))
os.system(cmd)

print('=== finished cap_louvain ===')
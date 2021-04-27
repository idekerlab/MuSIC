import pandas as pd
import numpy as np
import sys
import os
import argparse
import networkx as nx
cdDir = '/'.join(x for x in os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
sys.path.append(cdDir)
from music_utils import *

def to_pandas_dataframe(G):
    e = G.edges(data=True)
    df = pd.DataFrame()
    df['source'] = [x[0] for x in e]
    df['target'] = [x[1] for x in e]
    df['type'] = [x[2]['type'] for x in e]
    return df

def get_termStats(G, hiergeneset):
    clusters = list(set(list(G.nodes())) - hiergeneset)
    tsize_list = []
    cgene_list = []
    descendent_list = []
    for c in clusters:
        infoset = nx.descendants(G, c)
        cgeneset = infoset.intersection(hiergeneset)
        tsize_list.append(len(cgeneset))
        cgene_list.append(list(cgeneset))
        descendent_list.append(list(infoset - cgeneset))
    df = pd.DataFrame(index=clusters)
    df['tsize'] = tsize_list
    df['genes'] = cgene_list
    df['descendent'] = descendent_list
    return df

def jaccard(A, B):
    if type(A) != set:
        A = set(A)
    if type(B) != set:
        B = set(B)
    return len(A.intersection(B)) / len(A.union(B))

def clean_shortcut(G):
    edge_df = to_pandas_dataframe(G)
    edge_df.columns = ['parent', 'child', 'type']
    for idx, row in edge_df.iterrows():
        if len(list(nx.all_simple_paths(G, row['parent'], row['child']))) > 1:
            G.remove_edge(row['parent'], row['child'])
    return

def reorganize(G, hiergeneset, ci_thre, cluster_weight):
    iterate = True
    n_iter = 1
    while iterate:
        clear = True
        print('... starting iteration {}'.format(n_iter))
        ts_df = get_termStats(G, hiergeneset)
        ts_df.sort_values('tsize', ascending=False, inplace=True)
        for comp, row in ts_df.iterrows():
            tmp = ts_df[ts_df['tsize'] < row['tsize']]
            if tmp.shape[0] == 0:
                continue
            comp_geneset = set(row['genes'])
            descendent = row['descendent']
            for tmp_comp, tmp_row in tmp.iterrows():
                if tmp_comp in descendent:
                    continue
                tmp_comp_geneset = set(tmp_row['genes'])
                # Check if satisfy ci_thre
                if len(comp_geneset.intersection(tmp_comp_geneset))/tmp_row['tsize'] >= ci_thre:
                    # Check if child having higher weight than parent
                    if cluster_weight[comp] < cluster_weight[tmp_comp]:
                        G.add_edge(comp, tmp_comp, type='default')
                        clear = False
                        descendent += tmp_row['descendent']
        # Further clean up using networkx to remove shortcut edges
        clean_shortcut(G)
        # Update variables
        n_iter += 1
        if clear:
            iterate = False
    if n_iter == 2:
        modified = False
    else:
        modified = True
    return modified
    
def merge_parent_child(G, hiergeneset, ji_thre):
    # Delete child term if highly similar with parent term
    # One parent-child relationship at a time to avoid complicacies involved in potential long tail
    print('... start removing highly similar parent-child relationship')
    similar = True
    merged = False
    while similar:
        clear = True
        edge_df = to_pandas_dataframe(G)
        ts_df = get_termStats(G, hiergeneset)
        default_edge = edge_df[edge_df['type'] == 'default']
        for idx, row in default_edge.iterrows():
            if jaccard(ts_df.loc[row['source']]['genes'], ts_df.loc[row['target']]['genes']) >= ji_thre:
                print('# Cluster pair {}->{} failed Jaccard, removing cluster {}'.format(row['source'], row['target'], 
                                                                                         row['target']))
                clear = False
                merged = True
                parents = edge_df[edge_df['target'] == row['target']]['source'].values
                children = edge_df[edge_df['source'] == row['target']]['target'].values
                # Remove all parent->node edges
                for pnode in parents:
                    G.remove_edge(pnode, row['target'])
                for child_node in children:
                    etype = G[row['target']][child_node]['type']
                    # Remove all node->child edges
                    G.remove_edge(row['target'], child_node)
                    # Add all parent->child edges
                    for pnode in parents:
                        G.add_edge(pnode, child_node, type=etype)
                # Remove target node
                G.remove_node(row['target'])
                break
        if clear:
            similar = False
    # Clean up shortcuts introduced during node deleteing process
    clean_shortcut(G)
    return merged

def collapse_redundant(G, hiergeneset, cluster_weight, min_diff):
    # Delete child term if highly similar with parent term
    # One parent-child relationship at a time to avoid complicacies involved in potential long tail
    print('... start removing highly redundant systems')
    while True:
        edge_df = to_pandas_dataframe(G)
        ts_df = get_termStats(G, hiergeneset)
        default_edge = edge_df[edge_df['type'] == 'default']
        to_collapse = []
        for idx, row in default_edge.iterrows():
            parentSys, childSys, _ = row.values
            if parentSys not in cluster_weight:
                cluster_weight[parentSys] = 0
            if ts_df.loc[parentSys]['tsize'] - ts_df.loc[childSys]['tsize'] < min_diff:
                to_collapse.append([parentSys, childSys, cluster_weight[parentSys]])
        if len(to_collapse) == 0:
            return
        to_collapse = pd.DataFrame(to_collapse, columns=['parent', 'child', 'weight'])
        cidx = to_collapse['weight'].idxmin()
        deleteSys = to_collapse.loc[cidx]['child']
        print('# Cluster pair {}->{} highly redundant, removing cluster {}'.format(to_collapse.loc[cidx]['parent'], 
                                                                                   to_collapse.loc[cidx]['child'], 
                                                                                   deleteSys))
        parents = edge_df[edge_df['target'] == deleteSys]['source'].values
        children = edge_df[edge_df['source'] == deleteSys]['target'].values
        # Remove all parent->node edges
        for pnode in parents:
            G.remove_edge(pnode, deleteSys)
        for child_node in children:
            etype = G[deleteSys][child_node]['type']
            # Remove all node->child edges
            G.remove_edge(deleteSys, child_node)
            # Add all parent->child edges
            for pnode in parents:
                G.add_edge(pnode, child_node, type=etype)
        # Remove target node
        G.remove_node(deleteSys)

parser = argparse.ArgumentParser()
parser.add_argument('--outprefix', help='output_dir/file_prefix for the output file')
parser.add_argument('--ci_thre', type=float, default=0.75, help='Containment index threshold')
parser.add_argument('--ji_thre', type=float, default=0.9, 
                    help='Jaccard index threshold for merging similar clusters')
parser.add_argument('--minSystemSize', type=int, default=4, 
                    help='Minimum number of proteins requiring each system to have.')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology.')
parser.add_argument('--min_diff', type=int, default=1, help='Minimum difference in number of proteins for every parent-child pair.')
args = parser.parse_args()

outprefix = args.outprefix
minSystemSize = args.minSystemSize

# Get the weight threshold at which each system was found
f = '{}.clixoOutFile.sed'.format(outprefix)
ont_edge = pd.read_table(f, header=None)
ont_edge.columns = ['parent', 'child', 'type', 'alpha']
cluster_weight = dict(list(set(list(zip(np.asarray(ont_edge['parent'].values, dtype=str), 
                                        ont_edge['alpha'].values)))))
save_obj(cluster_weight, f+'.clusterWeightDict.pkl')

input_ddot = '{}.clixoOutFile.oaInput.mintsize{}'.format(outprefix, minSystemSize)
ci_thre = args.ci_thre
ji_thre = args.ji_thre
print('Containment index threshold: {}'.format(ci_thre))
print('Jaccard index threshold: {}'.format(ji_thre))

edge_df = pd.read_table(input_ddot, header=None, dtype=str)
edge_df.columns = ['parent', 'child', 'type']
input_n = edge_df.shape[0]
edge_df.drop_duplicates(inplace=True)
print('# {} duplicate edges removed'.format(input_n - edge_df.shape[0]))

hiergeneset = set(edge_df[edge_df['type'] == 'gene']['child'].values)
G = nx.from_pandas_dataframe(edge_df, source='parent', target='child', edge_attr='type', create_using=nx.DiGraph())
if not nx.is_directed_acyclic_graph(G):
    raise ValueError('Input hierarchy is not DAG!')

while True:
    modified = reorganize(G, hiergeneset, ci_thre, cluster_weight)
    merged = merge_parent_child(G, hiergeneset, ji_thre)
    if not modified and not merged:
        break

collapse_redundant(G, hiergeneset, cluster_weight, args.min_diff)
# Output as ddot edge file
clean_shortcut(G)
edge_df = to_pandas_dataframe(G)
edge_df.to_csv('{}.ddot'.format(outprefix), header=False, index=False, sep='\t')
cmd = '{}/ontologyTermStats {} genes > {}'.format(args.path_to_alignOntology.rstrip('/'), 
                                                  '{}.ddot'.format(outprefix), 
                                                  '{}.termStats'.format(outprefix))
os.system(cmd)
# Annotate cluster weight
edge_df['weight'] = [cluster_weight[x] for x in edge_df['source'].values]
edge_df.to_csv('{}.wClusterWeight.ddot'.format(outprefix), header=False, index=False, sep='\t')
print('=== finished mature_hier_structure.py ====')

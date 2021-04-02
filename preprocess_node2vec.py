import sys
import networkx as nx
from ndex.networkn import NdexGraph
import pickle
import os
import argparse
from file_utils import *


parser = argparse.ArgumentParser(description='Download BioPlex v2 and prepare for node2vec network embedding.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
args = parser.parse_args()

outprefix = args.outprefix

# Check if output directory exists, make directory if not exist
outdir = '/'.join(x for x in outprefix.split('/')[:-1])
if not os.path.exists(outdir):
    os.mkdir(outdir)

if not os.path.exists('{}/node2vec_emb_sim'.format(outprefix)):
    os.mkdir('{}/node2vec_emb_sim'.format(outprefix))
outputdir = '{}/node2vec_emb_sim'.format(outprefix)

bioplex = nx.Graph(NdexGraph(server='http://www.ndexbio.org', uuid='98ba6a19-586e-11e7-8f50-0ac135e8bacf'))
node_id = []
node_name = []
for node in bioplex.nodes(data=True):
    node_id.append(node[0])
    node_name.append(node[1]['name'])
idx_to_name = dict(zip(node_id, node_name))
name_to_idx = dict(zip(node_name, node_id))

with open('{}/node_idx_to_name.dict.pkl'.format(outputdir), 'wb') as f:
    pickle.dump(idx_to_name, f)
with open('{}/node_name_to_idx.dict.pkl'.format(outputdir), 'wb') as f:
    pickle.dump(name_to_idx, f)

#save edge file from bioplex w/ index as node names
f = open('{}/bioplex_edges.txt'.format(outprefix), "w")
for edge in bioplex.edges():
    f.write(str(edge[0]) + "\t" + str(edge[1]) + "\n")
f.close()



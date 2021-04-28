import sys
import os
import shutil
import argparse

parser = argparse.ArgumentParser(description='Generate bash file for community detection pipeline.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
args = parser.parse_args()
outprefix = args.outprefix

os.remove('{}.clixoOutFile'.format(outprefix))
os.remove('{}.clixoOutFile.oaInput'.format(outprefix))
os.remove('{}.clixoOutFile.oaInput.mintsize4'.format(outprefix))
os.remove('{}.clixoOutFile.sed'.format(outprefix))
os.remove('{}.clixoOutFile.sed.clusterWeightDict.pkl'.format(outprefix))
os.remove('{}.clixoOutFile.termStats'.format(outprefix))
os.remove('{}.ddot'.format(outprefix))
os.remove('{}.G.pkl'.format(outprefix))
os.remove('{}.mvp_partition_membership.1000.npy'.format(outprefix))
os.remove('{}.mvp_partition_modularity.1000.npy'.format(outprefix))
os.remove('{}.node_idx_to_name.dict.pkl'.format(outprefix))
os.remove('{}.node_name_to_idx.dict.pkl'.format(outprefix))
os.remove('{}.termStats'.format(outprefix))
os.remove('{}.wClusterWeight.ddot'.format(outprefix))
shutil.rmtree('{}_qr'.format(outprefix))
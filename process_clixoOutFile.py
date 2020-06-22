import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Process CliXO output file.')
parser.add_argument('--outprefix', help='Prefix for files generated. E.g. /path/to/output/directory/fileIdentifier')
parser.add_argument('--path_to_alignOntology', help='Full path to alignOntology.')
args = parser.parse_args()
outprefix = args.outprefix

print('... start processing CliXO output file')
clixo_outfname = '{}.clixoOutFile'.format(outprefix)
# == sed clixo output ==
sed_fname = clixo_outfname + '.sed'
cmd = 'sed \'/^#/ d\' ' + clixo_outfname + ' > ' + sed_fname
os.system(cmd)
print('... generated sed file with commented lines in CliXO output file deleted')
# == generate oa_input file ==
oa_fname = clixo_outfname + '.oaInput'
cmd = 'awk \'{print $1 "\t" $2 "\t" $3}\' ' + sed_fname + ' > ' + oa_fname
os.system(cmd)
print('... generated oa_input files with the fourth score column from sed files deleted')
# == generate termStats file ===
ts_fname = clixo_outfname + '.termStats'
cmd = '{}/ontologyTermStats {} genes > {}'.format(args.path_to_alignOntology.rstrip('/'), oa_fname, ts_fname)
os.system(cmd)
print('... generated termStat files')

print('=== finished process_clixo_output.py ===')
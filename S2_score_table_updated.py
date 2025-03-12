#! /usr/bin/env python
#
# This work is copyrhight 2024 Ansgar Gruber and Marta Vohnoutová, University of South Bohemia,
# and licenced under CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International) licence,
# available at: https://creativecommons.org/licenses/by-sa/4.0/
#
# Main changes are the adjustments of package requirements. 
# For updates see our:
# Repository on GitHub: https://github.com/ASAFind/ASAFind-2
# ASAFind web service: https://asafind.jcu.cz/
# Contacts:
# Ansgar Gruber <ansgar.gruber@paru.cas.cz>
# Marta Vohnoutová <mvohnoutova@prf.jcu.cz>
#
#
# This work is derivative of Score_Table.py, which is copyright 2020 Cedar McKay and Gabrielle Rocap,
# University of Washington, licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License,
# available at: http://creativecommons.org/licenses/by-sa/3.0/
# Score_Table.py is accessible via the following url:
# https://bitbucket.org/rocaplab/asafind/src/30f925e2684e1f9ee08df7af24df98bf5ee5fdb3/ASAFind2.py
#
#


VERSION = '1.0_beta5'
PROJ = 'Score_Table'

import sys
import os.path
import inspect
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
# The FreqTable import has been replaced with dictionary-based frequency handling.
import math
import pickle
from fill_constants import *

if sys.version_info < (3, 10):
    print(inspect.cleandoc('''\n\nWARNING: This version of Score_Table was developed on Python 3.6
        and tested on python 3.10 It is not tested on Python version:
        {}.{}\n\n'''.format(sys.version_info[0], sys.version_info[1])))

#### Collect Input ####
parser = argparse.ArgumentParser(description=f'''
{PROJ} {VERSION}

Takes a Fasta and generates a scoring table based on the aa frequencies at each position.
The output of this script is a tab delimited table.
Python >= 3.6 required.''',
formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', '--fasta_file',
                    help='Specify the input fasta FILE. Must not contain gaps.', required=True)

parser.add_argument('-o', '--out_file', help='Specify the path and name of the output file you'
                    'wish to create. Default will be the same as the fasta_file, but with a '
                    '".tab" suffix.')

parser.add_argument('--version', action='version', version=f'{PROJ} {VERSION}')

args = parser.parse_args()

# Figure out some names and paths
fasta_file = os.path.abspath(args.fasta_file)
(fasta_file_path, fasta_file_name) = os.path.split(fasta_file)
(fasta_file_base_name, fasta_file_ext) = os.path.splitext(fasta_file_name)
# Figure out what our out_file is.
if args.out_file:
    out_file = os.path.abspath(args.out_file)
else:
    out_file = os.path.abspath(os.path.join(fasta_file_path, fasta_file_base_name + '.tab'))


#####open log file####
orig_stdout = sys.stdout
log_file=open(f'{temp_file}log_S2_{fasta_file_name}.txt', 'w')
sys.stdout = log_file
print(f'Log file for S2 program {fasta_file_name}')
###################

#### Main ####
##############

# Test input files to make sure they are parseable
try:
    records_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
except ValueError as e:
    raise Exception(f'The input fasta could not be parsed. The error is: {e}')
# Check if all sequences are equal length, ungapped, no illegal characters
valid_letters = 'ACDEFGHIKLMNPQRSTVWY'
lengths = []
for gene in records_dict.keys():
    seq = records_dict[gene].seq.upper()
    for aa in seq:
        if aa == '-':
            raise ValueError(f'No gaps permitted in alignment, found gap in {gene}')
        if aa not in valid_letters:
            raise ValueError(f'Letters must be one of {valid_letters}, found: "{aa}" in {gene}')
    lengths.append(len(seq))
if len(set(lengths)) > 1:
    raise Exception(f'Sequences must be of equal length')

# Now get to the good stuff

# set up aa list and generic expected frequency table
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
random_freq = [0.05] * 20
random_freq_dict = dict(zip(aas, random_freq))
random_freq_table = random_freq_dict  # Using a simple dictionary for frequencies

# read alignment file into alignment object
transit_align = AlignIO.read(fasta_file, 'fasta')
num_taxa = len(transit_align)
align_length = transit_align.get_alignment_length()

#print(f'Transit alignment {transit_align} has {num_taxa} taxa and {align_length} positions') #Marta: I added this line to check the alignment

# calculate e_n for this sample size
e_n = (1/(math.log(2))) * (19/(2*num_taxa))

# get alignment info and PSSM objects
Alignment_Info = SummaryInfo(transit_align)
PSSM = Alignment_Info.pos_specific_score_matrix()

keys = PSSM.pssm[0][1].keys()
missing_keys=[]
for k in aas:
    if k not in keys:
        missing_keys.append(k)

if len(missing_keys) > 0:    
    for i in PSSM.pssm:
        i[1].update({k: 0 for k in missing_keys})
print(f'{missing_keys} were added to the PSSM') #Marta: I added this line to check the PSSM


# make a list with the information content at each position in it
Pos_Info = []
for pos in range(align_length):
    InfContent = Alignment_Info.information_content(start=pos, end=pos+1,
		e_freq_table=random_freq_table)
    InfContent = InfContent - e_n
    Pos_Info.append(InfContent)

#print(f'Positional Information {Pos_Info}') #Marta:

# set up the matrix
rows = align_length + 1     # extra for header
cols = len(aas) + 1         # extra for X row
weight_matrix = [[0 for i in range(cols)] for j in range(rows)]

# go through each position and calculate information for each aa at that position
max_score = 0
weight_matrix[0] = aas
for pos in range(align_length):
    weights = []
    for aa in aas: #Marta: I have to change this to aas_abbr
        try:
            weight = (PSSM[pos][aa]/num_taxa)*Pos_Info[pos]
            weights.append(weight)
        except KeyError:
            print(f'KeyError at position {pos} and aa {aa}')
            #print(PSSM)
    weight_matrix[pos + 1] = weights
    max_score = max_score + max(weights)

weight_matrix_transposed = [list(x) for x in zip(*weight_matrix)]

# add the X row, X contains no information, so all 0
X_row = [0] * align_length
X_row.insert(0, 'X')
weight_matrix_transposed.append(X_row)

# print some useful output about this matrix
print(f'''The maximum possible sequence score for this scoring table based on {num_taxa}
sequences is {max_score}''')

# write matrix (2Dlist) to string
scoring_table_string = ''
for row in weight_matrix_transposed:
    for elem in row:
        scoring_table_string = scoring_table_string + str(elem) + '\t'
    scoring_table_string = scoring_table_string + '\n'

# write string to file
out_handle = open(out_file, 'w')
out_handle.write(scoring_table_string)
out_handle.close()  

# save to pickle after targetp step
with open(out_file+'.pkl', 'wb') as f:
    pickle.dump(scoring_table_string, f, pickle.HIGHEST_PROTOCOL)
print(type(scoring_table_string))

print('Scoring table saved to', out_file)
print('Scoring table in pickle saved to', out_file+'.pkl')


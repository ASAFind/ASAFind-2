#! /usr/bin/env python
#
# This work is copyrhight 2025 Ansgar Gruber and Marta Vohnoutová, University of South Bohemia,
# and licenced under CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International) licence,
# available at: https://creativecommons.org/licenses/by-sa/4.0/
# 
#
# For updates see our repository on GitHub: https://github.com/ASAFind/ASAFind-2
# ASAFind web service: https://asafind.jcu.cz/
# Publication: https://doi.org/10.1111/tpj.70138 (please cite if you publish results of the scripts, or derivative work)
# 
# Contacts:
# Ansgar Gruber <ansgar.gruber@paru.cas.cz>
# Marta Vohnoutová <mvohnoutova@prf.jcu.cz>
#
#
# This work is derivative of ASAFind2.py version 2.0_beta19, which is copyright 2020 Cedar McKay and Gabrielle Rocap,
# University of Washington, licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License,
# available at: http://creativecommons.org/licenses/by-sa/3.0/
# ASAFind2.py version 2.0_beta19 is accessible via the following url:
# https://bitbucket.org/rocaplab/asafind/src/30f925e2684e1f9ee08df7af24df98bf5ee5fdb3/ASAFind2.py
# Main change is the inclusion of targeting predictions for the periplastidic compartment
#
#

VERSION = '2.0'
PROJ = 'ASAFind'

import sys
import os.path
import inspect
from collections import defaultdict
import argparse


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from fill_constants import *
import pickle
#import pickle5 as pickle

# Check Python version
if sys.version_info < (3, 10):
    print(inspect.cleandoc('''\n\nWARNING: This version of ASAFind was developed on Python 3.10, it is not tested on Python version {}.{}.\n\n'''.format(sys.version_info[0], sys.version_info[1])))

#### Collect Input ####
#######################

parser = argparse.ArgumentParser(description=f'''
{PROJ} {VERSION}

Takes a Fasta and companion SignalP (versions 3.0, 4.0, 4.1, 5.0) or TargetP 2.0 short format file
as input, with the complete SignalP or TargetP header (two lines starting with '#'). Some versions of SignalP
truncate the sequence names. SignalP-3.0 to 20 characters, and 4.0, 4.1 to 58 characters. Therefore,
ASAFind only considers the first corresponding characters of the fasta name (and the first 90 in the
case of SignalP 5.0 and TargetP 2.0), which must be unique within the file. Parts of the fasta name
after that character are ignored. Additionally, the fasta name may not contain a '-' or '|'. This
requirement is because SignalP converts characters in sequence names (e.g. '-' is changed to '_').
ASAFind requires at least 7 aa upstream and 22 aa downstream of the cleavage site suggested by
SignalP. The output of this script is a tab delimited table. 
Python >= 3.6 required.''',
formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-f', '--fasta_file', help='Specify the input fasta FILE.',
                    required=True)

parser.add_argument('-p', '--signalp_file',
                    help='Specify the input SignalP or TargetP FILE.', required=True)

parser.add_argument('-s', '--simple_score_cutoff', default=None, type=float,
                    help='Optionally, specify an explicit score cutoff, rather than '
                    'using ASAFind\'s default algorithm, not compatible with option -v1. The score given here will not be normalized and therefore should be obtained from a distribution of normalized scores.')

parser.add_argument('-t', '--score_table_file',
                    help='Optionally, specify a custom scoring table. The scoring table will be normalized with the maximum score, which allows for processing of non-normalized as well as normalized scoring tables.')

parser.add_argument('-o', '--out_file', help='Specify the path and name of the output file '
                    'you wish to create. Default will be the same as the fasta_file, but with a '
                    '".tab" suffix.')

parser.add_argument('-w', '--web_output', action='store_true',
                    help='Format output for web display. This is mostly useful when called by a '
                    'web app.')

parser.add_argument('-v1', '--reproduce_ASAFind_1', action='store_true',
                    help='Reproduce ASAFind 1.x scores and results (non-normalized scores, if no custom scoring table is specified, the original default scoring table generated without small sample size correction will be used, not compatible with option -s).')

parser.add_argument('-ppc', '--include_ppc_prediction', action='store_true',
                    help='Include prediction of proteins that might be targeted to the periplastidic compartment.')

parser.add_argument('-s_ppc', '--score_cutoff_ppc', default=None, type=float,
                    help='Optionally, specify an explicit score cutoff for the ppc protein prediction, if given, ppc protein prediction will be included. The score given here will not be normalized and therefore should be obtained from a distribution of normalized scores.')

parser.add_argument('-t_ppc', '--score_table_file_ppc',
                    help='Optionally, specify a custom scoring table for the ppc protein prediction, if given, ppc protein prediction will be included. The scoring table will be normalized with the maximum score, which allows for processing of non-normalized as well as normalized scoring tables.')

parser.add_argument('-v', '--version', action='version', version=f'{PROJ} {VERSION}')

args = parser.parse_args()


#### Variables and Names ####
#############################

# Figure out some names and paths
fasta_file = os.path.abspath(args.fasta_file)
signalp_file = os.path.abspath(args.signalp_file)

(fasta_file_path, fasta_file_name) = os.path.split(fasta_file)
(fasta_file_base_name, fasta_file_ext) = os.path.splitext(fasta_file_name)

if args.score_table_file:
    custom_scoring_table = open(args.score_table_file, 'r').read()
else:
    custom_scoring_table = None

if args.score_table_file_ppc:
    custom_scoring_table_ppc = open(args.score_table_file_ppc, 'r').read()
else:
    custom_scoring_table_ppc = None
    
simple_score_cutoff = args.simple_score_cutoff
score_cutoff_ppc = args.score_cutoff_ppc
web_output = args.web_output
reproduce_ASAFind_1 = args.reproduce_ASAFind_1
include_ppc_prediction = args.include_ppc_prediction

#####open log file####
orig_stdout = sys.stdout
log_file=open(f'{temp_file}log_S1_{fasta_file_name}.txt', 'w')
sys.stdout = log_file
print(f'Log file for S1 program {fasta_file_name}')
###################


# Figure out what our out_file is.
if args.out_file:
    out_file = os.path.abspath(args.out_file)
else:
    out_file = os.path.abspath(os.path.join(fasta_file_path, fasta_file_base_name + '.tab'))

# Check compatibility of options

if simple_score_cutoff and reproduce_ASAFind_1:
    raise Exception('options -s and -v1 are not campatible')

if args.score_cutoff_ppc:
    include_ppc_prediction = True

if args.score_table_file_ppc:
    include_ppc_prediction = True

if args.score_table_file and args.simple_score_cutoff is False:
    print ('You specified a custom scoring table, but no custom scoring cut-off. Please be aware that a custom scoring table most likely works best with a custom score cut-off.')

if args.score_table_file_ppc and args.score_cutoff_ppc is False:
    print ('You specified a custom scoring table for the ppc protein prediction, but no custom scoring cut-off. Please be aware that a custom scoring table most likely works best with a custom score cut-off.')

#### Functions ####
###################

# modified by Marta
def clean_aa_sequence(seq_record):
    '''Takes a Bio.Seq object. Returns a Bio.Seq.Seq object. Checks the validity of a protein
    aa sequence. 1) Must be entirely composed of ACDEFGHIKLMNPQRSTVWYBZJUOX* 2) Replace BZJUO
    with X'''
    valid_letters = 'ACDEFGHIKLMNPQRSTVWYBZJUOX*'
    replace_these_letters = 'BZJUO-'
    #seq = seq_record.upper()
    seq_mutable = str(seq_record.upper())
    #seq_mutable = MutableSeq(str(seq))
    #print(f'Seq mutable {seq_mutable}')
    i = 0
    for aa in seq_mutable:
        if aa not in valid_letters + replace_these_letters:
            raise ValueError(f'Letters must be one of {valid_letters}, found character: "{aa}"')
        if aa in replace_these_letters:
            if 0 <= i < len(seq_mutable):
                seq_mutable = seq_mutable[:i] + 'X' + seq_mutable[i+1:]
            else:
                raise IndexError("Index out of range")
            snippet = seq_mutable[max(0,i-10):i+10]
            #print(f'WARNING: Replaced "{aa}" in "{snippet}" with an "X" in {seq_record.name}')
        i += 1
    #print((seq_mutable))
    return  seq_mutable
            

def parse_signalp(p_file):
    '''Processes a SignalP or TargetP output file, and returns a parsed dictionary of results

    Args: Path to a SignalP 3.0, SignalP 4.0, SignalP 4.1, SignalP 5.0, or TargetP 2.0 format
    file, short format, header required.

    Returns: Dict with parsed keys. Each value is a tuple (score, signalp-Y/N, position). Smean
    and D don't have a position.

    '''
    d = {}  # results dict
    p_handle = open(p_file, 'r')
    seen = set()

    # Account for different column order. Thanks guys.
    line = next(p_handle)  # First line
    if not line.startswith('#'):
        print('Could not find a header. Are you sure this is a valid input file?')
        raise ValueError('The SignalP or TargetP file is an unrecognized format')
    if line.startswith('# SignalP-NN'):  # or len(line.split()) == 21:
        d['signalp_version'] = 'SignalP-3.0'
        d['character_limit'] = 20
        organism = line.split()[2]
        if organism != 'euk':
            print(inspect.cleandoc(f'''######################################## WARNING WARNING
                WARNING: ASAFind is designed for SignalP 3 Eukaryotes predictions, your input
                file used {organism} predictions. We will proceed, but you may want to double
                check your SignalP settings.
                ########################################'''))
        Cmax_score = 1
        Cmax_position = 2
        Ymax_score = 4
        Ymax_position = 5
        Smax_score = 7
        Smax_position = 8
        Smean_score = 10
        D_score = 12
        D_conclusion = 13

    elif line.startswith('# SignalP-4'):  # or len(line.split()) == 12:
        if line.startswith('# SignalP-4'):
            d['signalp_version'] = line.split()[1]  # Right now this can be 4.0 or 4.1
        else:
            d['signalp_version'] = 'SignalP-4.x'
        d['character_limit'] = 58
        organism = line.split()[2]
        if organism != 'euk':
            print(inspect.cleandoc(f'''######################################## WARNING WARNING
                WARNING: ASAFind is designed for SignalP 4 Eukaryotes predictions, your input
                file used {organism} predictions. We will proceed, but you may want to double
                check your SignalP settings.
                ########################################'''))
        Cmax_score = 1
        Cmax_position = 2
        Ymax_score = 3
        Ymax_position = 4
        Smax_score = 5
        Smax_position = 6
        Smean_score = 7
        D_score = 8
        D_conclusion = 9

    elif line.startswith('# SignalP-5.0'):  # or len(line.split()) in (4, 6, 10):
        d['signalp_version'] = 'SignalP-5.0'
        d['character_limit'] = 90
        organism = line.split('\t')[1].split(':')[1].strip()
        if organism != 'Eukarya':
            print(inspect.cleandoc(f'''######################################## WARNING WARNING
                WARNING: ASAFind is designed for SignalP 5 Eukarya predictions, your input file
                used {organism} predictions. We will proceed, but you may want to double check
                your SignalP settings.
                ########################################'''))
            Position = 6  # the other organism settings have more columns before the CS position
        else:
            Position = 4
        SP = 1

    elif line.startswith('# TargetP-2.0'):  # or len(line.split()) in (5, 7, 11):
        d['signalp_version'] = 'TargetP-2.0'
        d['character_limit'] = 90
        organism = line.split('\t')[1].split(':')[1].strip()
        if organism != 'Non-Plant':
            print(inspect.cleandoc(f'''######################################## WARNING WARNING
                WARNING: ASAFind is designed for TargetP 2.0 "Non-Plant" predictions, your input
                file used {organism} predictions. We will proceed, but you may want to double
                check your TargetP settings.
                ########################################'''))
            Position = 7  # the other organism settings have more columns before the CS position
        else:
            Position = 5

        SP = 1

    else:
        raise ValueError('The SignalP or TargetP file is an unrecognized format')

    p_handle.seek(0)  # In case we didn't have a header, we need to go back to beginning of file.
    for line in p_handle:
        if not line.startswith('#') and len(line) > 1:
            atoms = line.split()
            id = atoms[0][:d['character_limit']]
            if id in seen:
                raise ValueError(inspect.cleandoc(f'''The SignalP or TargetP file has entries not
                    unique within first {d['character_limit']} characters: {id}'''))

            elif d['signalp_version'] in ('SignalP-3.0', 'SignalP-4.0', 'SignalP-4.1'):
                seen.add(id)
                d[id] = {'Cmax': {'score': atoms[Cmax_score], 'position': int(atoms[Cmax_position])},
                         'Ymax': {'score': atoms[Ymax_score], 'position': int(atoms[Ymax_position])},
                         'Smax': {'score': atoms[Smax_score], 'position': int(atoms[Smax_position])},
                         'Smean': {'score': atoms[Smean_score]}, 'D': {'score': atoms[D_score],
                         'conclusion': atoms[D_conclusion]}}

            elif d['signalp_version'] in ('SignalP-5.0', 'TargetP-2.0'):
                seen.add(id)
                d[id] = {'SP': atoms[SP], 'Position': atoms[Position:]}

            else:
                raise ValueError(inspect.cleandoc(f'''The SignalP or TargetP file is an unrecognized
                    format. Did not recognize the version: {d['signalp_version']}'''))
    return d


def parse_score_table(scoring_table_string=None, reproduce_ASAFind_1=None):
    '''Takes a scoring table, and returns a 2D dict (like a lookup table), default tables are used if no custom table is specified, ASAFind 1 and ASAFind 2 default tables are stored

    The input table should be in the following format:

    A   0.151009115 0.085083382 0.583560009
    C   0.031243265 0.026179502 0.063661092
    D   0.005207211 0   0

    For example, the score weight given to a 'C' in the second position is 0.026179502

    Args:
        Scoring table as a string. If none, the default table will be used.

    Returns:
        A 2D dict describing the table. Access the scoring weight to a 'C' in the second
        position like this: d['C'][2]. Note position is not zero based, to match SignalP
        results.
    '''

    if not scoring_table_string:
        if reproduce_ASAFind_1:
            with open(f'{where_are_programs}/reproduce_ASAFind_1.tab.pkl', 'rb') as f:
                scoring_table_string = pickle.load(f)
        else:   # newly added from diatom_scoring_matrix.fasta counted by S2 program
            with open(f'{where_are_programs}/diatom_output.tab.pkl', 'rb') as f:
                scoring_table_string = pickle.load(f)	

    d = defaultdict(dict)  # results dict
    for line in scoring_table_string.split('\n'):
        atoms = line.split()

        #Ignore blank lines
        if len(atoms) == 0: 
        	continue
        # Test input
        if len(atoms) < 26:
            raise ValueError(f'Expected 26 columns, not {len(atoms)}')
        elif len(atoms) > 26:
            print(f'''WARNING Expected 26 columns, not {len(atoms)}. Will proceed, but ASAFind will
                only consider the first 26 columns.''')
        aa = atoms.pop(0)  # amino acid short code
        i = 1  # positions are 1 based
        while len(atoms) > 0:
            try:
                atom = atoms.pop(0)
                score = float(atom)
            except ValueError:
                print(f'Expected a number, but found "{atom}" instead.')
                raise()
            d[aa][i] = score
            i += 1

    return d

def parse_score_table_ppc(scoring_table_string_ppc=None):
    '''Takes a scoring table, and returns a 2D dict (like a lookup table), default ppc protein prediction table is used if no custom table is specified

    The input table should be in the following format:

    A   0.151009115 0.085083382 0.583560009
    C   0.031243265 0.026179502 0.063661092
    D   0.005207211 0   0

    For example, the score weight given to a 'C' in the second position is 0.026179502

    Args:
        Scoring table as a string. If none, the default table will be used.

    Returns:
        A 2D dict describing the table. Access the scoring weight to a 'C' in the second
        position like this: d['C'][2]. Note position is not zero based, to match SignalP
        results.
    '''

    if not scoring_table_string_ppc:
        with open('ppc_output.tab.pkl', 'rb') as f:
            scoring_table_string_ppc = pickle.load(f)
    d = defaultdict(dict)  # results dict
    for line in scoring_table_string_ppc.split('\n'):
        atoms = line.split()

        #Ignore blank lines
        if len(atoms) == 0: 
        	continue
        # Test input
        if len(atoms) < 26:
            raise ValueError(f'Expected 26 columns, not {len(atoms)}')
        elif len(atoms) > 26:
            print(f'''WARNING Expected 26 columns, not {len(atoms)}. Will proceed, but ASAFind will
                only consider the first 26 columns.''')
        aa = atoms.pop(0)  # amino acid short code
        i = 1  # positions are 1 based
        while len(atoms) > 0:
            try:
                atom = atoms.pop(0)
                score = float(atom)
            except ValueError:
                print(f'Expected a number, but found "{atom}" instead.')
                raise()
            d[aa][i] = score
            i += 1

    return d

def get_max_tscore(score_dict):
    '''Takes a scoring table dictionary and returns the max possible transit peptide score
    
    Args:
        A  2D dict representation of a scoring table.

    Returns:
    	float of max possible score
    '''
    max_tscore = 0
    for pos in range(6,26):
        scores = []
        for aa in score_dict.keys():
            scores.append(score_dict[aa][pos])
        max_pos = max(scores)    
        max_tscore = max_tscore + max_pos
    return max_tscore

def get_max_cscore(score_dict):
    '''Takes a scoring table dictionary and returns the max possible cleavage site score
    
    Args:
        A  2D dict representation of a scoring table.

    Returns:
    	float of max possible score
    '''
    max_cscore = 0
    for pos in range(1,26):
        scores = []
        for aa in score_dict.keys():
            scores.append(score_dict[aa][pos])
        max_pos = max(scores)    
        max_cscore = max_cscore + max_pos
    return max_cscore

def score_peptide(seq, table=None):
    '''Takes a 25mer and calculates the cleavage site score and transit peptide score.

    Args:
        A 25 aa long string, and a 2D dict representation of a scoring table.

    Returns:
        A tuple of scores. First is cleavage site score, 2nd is transit peptide score.
    '''

    if not isinstance(table, dict):
        raise TypeError('Expected a scoring table as a 2D dict')
    if len(seq) != 25:
        #print(f'{str(seq)}, lenght seq {len(seq)} type seq {type(seq)}')  #marta added
        raise ValueError('Expected a 25mer')
    seq = seq.upper()
    cleavage_score = 0
    transit_score = 0

    i = 1  # 1 based counting conforming to signalP
    for aa in seq:
        try:                    # Marta added
            cleavage_score += table[aa][i]
        except KeyError:
            print(f'KeyError!!!!!!!: {aa} {i} {table}')  # Marta added
            pass

        if i > 5:  # transit score is calculated omitting 1st 5 aa.
            try:
                transit_score += table[aa][i]
            except KeyError:
                print(f'KeyError!!i>5!!!!!: {aa} {i} {table}')  # Marta added
                pass
        i += 1

    return (cleavage_score, transit_score)


#### Main ####
##############


if web_output:
    print('<PRE>')

# Test input files to make sure they are parseable, valid, no dups.
try:
    records_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
except ValueError as e:
    raise ValueError(f'The input fasta could not be parsed. The error is: {e}')


try:
    signalps = parse_signalp(signalp_file)  # It parses.That means no duplicate entries. This is a dict

    records_keys_short = [x[:signalps['character_limit']] for x in list(records_dict.keys())]
    problem_names = [x for x in records_keys_short if any([s in x for s in ('-', '|')])]
    if problem_names:
        print('Problem names: ', problem_names)
        raise ValueError('''SignalP changes "-" and "|" characters in names to "_", therefore
            ASAFind does not allow those characters either. Therefore, these files cannot be analysed, please change the sequence names so that
they match between the files.''')

    if not len(records_keys_short) == len(set(records_keys_short)):
        dups = tuple(set([x for x in records_keys_short if records_keys_short.count(x) > 1]))
        raise ValueError(inspect.cleandoc(f'''The input fasta has records that are not unique to the
            first {signalps['character_limit']} characters.  Therefore, these files cannot be analysed, please change the sequence names so that
each name in the FASTA file is unique, and that the names match between the files. The offending records begin: {dups}'''))

    if not set(signalps.keys()) >= set(records_keys_short):
        missing_entries = [x for x in set(records_keys_short) - set(signalps.keys())]
        raise ValueError(inspect.cleandoc(f'''SignalP/TargetP file does not have an entry (unique to first
            {signalps['character_limit']} characters) for every sequence in fasta file. Therefore, these files cannot be analysed, please change the SignalP/TargetP file so that
it contains exactly one prediction for each sequence in the FASTA file, with matching sequence names between the files. Missing
            entry for {missing_entries}'''))
except Exception as e:
    print(f'The SignalP/TargetP file could not be validated. The error is:\n{e}')
    raise(e)


# Input files passed our tests. Onward!
out_handle = open(out_file, 'w')
records = SeqIO.parse(fasta_file, 'fasta')
scoring_table = parse_score_table(custom_scoring_table, reproduce_ASAFind_1)
scoring_table_ppc = parse_score_table_ppc(custom_scoring_table_ppc)
maximum_tscore = get_max_tscore(scoring_table)
maximum_cscore = get_max_cscore(scoring_table)
if include_ppc_prediction:
    maximum_tscore_ppc = get_max_tscore(scoring_table_ppc)
    maximum_cscore_ppc = get_max_cscore(scoring_table_ppc)
    if not score_cutoff_ppc:
        score_cutoff_ppc = 0.49

if reproduce_ASAFind_1:
    default_score_cutoff = 2
else:
    default_score_cutoff = (2/maximum_tscore)

# Initialize stats
total_proteins = 0
signalp_positive = 0
chloroplast_targeted = 0
chloroplast_targeted_high = 0
chloroplast_targeted_low = 0
ppc_targeted = 0
targetp_mtp = 0
targetp_other = 0
skipped_proteins = 0


results = []

for record in records:
    #print(f'Record {record}')
    total_proteins += 1
    if include_ppc_prediction:
        best_window = {'window': 0, 'score_25mer': 0, 'score_20mer': 0, 'score_20mer_ppc': 0, 'peptide_25mer': 0}
    else:
        best_window = {'window': 0, 'score_25mer': 0, 'score_20mer': 0, 'peptide_25mer': 0}

    # Get a single parsed signalp result for the record. Some version of SignalP only uses 20
    # chars of the name
    signalp = signalps[record.name[0:signalps['character_limit']]] 
 

    if signalps['signalp_version'] in ('SignalP-5.0', 'TargetP-2.0'):
        signalp_marker = signalp['SP']
        try:
            if signalp['Position']:
                #print(signalp['Position']) # Marta added
                #signalp_parsed_position=int(row.split('\t')[4].split('-')[0][-2:])
                signalp_parsed_position = int(signalp['Position'][2].split('-')[1].split('.')[0])
                #print(signalp_parsed_position)
            else:
                signalp_parsed_position = 0  # 5.0 leaves this column blank if no prediction.
        except ValueError:
#            print(f'Value Error in predicted cleavage position of {record.name}.')
            pass # ValueError here does not disturb the calculations, program should be able to complete (tested)
        except IndexError:  # Marta added
#            print(f'Index Error in predicted cleavage position of {record.name}.') 
            pass # IndexError here does not disturb the calculations, program should be able to complete (tested)

    else:  # signalp_version is not 5.0
        signalp_marker = signalp['D']['conclusion']
        signalp_parsed_position = signalp['Ymax']['position']

    # There are 3 scenarios to contemplate. 1) Positive signal, but too short, 2) Positive signal,
    # and long enough, 3) negative signal.
    # First, tackle scenario 1.

    if signalp_marker in ('SP', 'SP(Sec/SPI)', 'Y'):  # SP(Sec/SPI) Equiv to old D conclusion == 'Y'.
                        # LIPO(Sec/SPII), TAT(Tat/SPI), OTHER, N,  all equiv to NO
        signalp_positive += 1

        if len(record.seq[signalp_parsed_position-8: signalp_parsed_position+21]) < 29:
            results.append((record.name,
                           signalp_marker,
                           'NA',
                           'NA',
                           'NA',
                           'NA',
                           'No prediction, sequence too short'))

            skipped_proteins += 1

        else:  # Not too short, so tackle scenario 2

            sliding_window_range = (-2, -1, 0, 1, 2)  # Hardcode since no plans to change

            for window in sliding_window_range:
                # subtract 1 to convert to 0 based positions
                # 25mer starts 5 before cleavage position.
                start = signalp_parsed_position + window - 5 - 1
                stop = start + 25
                #print(f'Record {record[start:stop]}')
                peptide_25mer = clean_aa_sequence(record.seq[start:stop])
                this_window_score = score_peptide(peptide_25mer, scoring_table)

                # Score returned is a tuple. 0 is 25mer score and 1 is 20mer score
                if not reproduce_ASAFind_1:
                    this_window_score = ((this_window_score[0]/maximum_cscore), (this_window_score[1]/maximum_tscore))
                
                if include_ppc_prediction:
                    this_window_score_ppc = score_peptide(peptide_25mer, scoring_table_ppc)
                    this_window_score_ppc = ((this_window_score_ppc[0]/maximum_cscore_ppc), (this_window_score_ppc[1]/maximum_tscore_ppc))
                    if this_window_score[0] > best_window['score_25mer']:
                        best_window = {'window': window, 'score_25mer': this_window_score[0],
                                       'score_20mer': this_window_score[1],
                                       'score_20mer_ppc': this_window_score_ppc[1],
                                       'peptide_25mer': peptide_25mer}
    
                    if window == 0:  # Store the original predicted window when we run across it.
                        window_0 = {'window': window, 'score_25mer': this_window_score[0],
                                    'score_20mer': this_window_score[1],
                                    'score_20mer_ppc': this_window_score_ppc[1],
                                    'peptide_25mer': peptide_25mer}
                else:    
                    if this_window_score[0] > best_window['score_25mer']:
                        best_window = {'window': window, 'score_25mer': this_window_score[0],
                                       'score_20mer': this_window_score[1],
                                       'peptide_25mer': peptide_25mer}
    
                    if window == 0:  # Store the original predicted window when we run across it.
                        window_0 = {'window': window, 'score_25mer': this_window_score[0],
                                    'score_20mer': this_window_score[1],
                                    'peptide_25mer': peptide_25mer}

            # simple_score_cutoff algorithm
            if simple_score_cutoff:
                if best_window['score_20mer'] > simple_score_cutoff:
                    best_window['prediction'] = 'Plastid'
                    chloroplast_targeted += 1
                    chloroplast_targeted_high += 1
                elif include_ppc_prediction and best_window['score_20mer_ppc'] > score_cutoff_ppc:
                    best_window['prediction'] = 'PPC'
                    ppc_targeted += 1
                elif include_ppc_prediction:
                    best_window['prediction'] = 'Not plastid, not ppc, signal peptide identified'
                else:
                    best_window['prediction'] = 'Not plastid, signal peptide identified'

            # normal algorithm
            else:
                # Check if 1st of cleaved peptide of best window starts with F,W,Y, L
                plus1_aa = best_window['peptide_25mer'][5].upper()
                if plus1_aa in 'FWYL':
                    best_window['plus1_aa'] = plus1_aa

                    if best_window['window'] == 0 and best_window['score_20mer'] > default_score_cutoff:
                        best_window['prediction'] = 'Plastid, high confidence'
                        chloroplast_targeted += 1
                        chloroplast_targeted_high += 1
                    else:
                        best_window['prediction'] = 'Plastid, low confidence'
                        chloroplast_targeted += 1
                        chloroplast_targeted_low += 1
                elif include_ppc_prediction and best_window['score_20mer_ppc'] > score_cutoff_ppc:
                    best_window['plus1_aa'] = 'No'
                    best_window['prediction'] = 'PPC'
                    ppc_targeted += 1
                elif include_ppc_prediction:
                    best_window['plus1_aa'] = 'No'
                    best_window['prediction'] = 'Not plastid, not ppc, signal peptide identified'
                else:
                    best_window['plus1_aa'] = 'No'
                    best_window['prediction'] = 'Not plastid, signal peptide identified'

            # Calculate a few values that would be too cumbersome to put right into the results list
            window_0['position'] = signalp_parsed_position
            best_window['position'] = signalp_parsed_position + best_window['window']
            # -1 business is to convert to 0 based counting
            cleavage_seq_start = record.seq[best_window['position']-3-1:best_window['position']-1]
            cleavage_seq_end = record.seq[best_window['position']-1:best_window['position']+3-1]

            if include_ppc_prediction:
                results.append((record.name,
                                signalp_marker,
                                best_window['position'],
                                f'{cleavage_seq_start}-{cleavage_seq_end}',
                                best_window['window'],
                                best_window['score_20mer'],
                                best_window['score_20mer_ppc'],
                                best_window['prediction']))
            else:
                results.append((record.name,
                                signalp_marker,
                                best_window['position'],
                                f'{cleavage_seq_start}-{cleavage_seq_end}',
                                best_window['window'],
                                best_window['score_20mer'],
                                best_window['prediction']))
    # Finally, deal with scenario 3
    elif signalp_marker in ['mTP', 'luTP', 'LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'OTHER', 'N', 'noTP']:
        if signalps['signalp_version'] == 'TargetP-2.0' and signalp_marker == 'mTP': # count mitochondrial predictions only if TargetP is used
            targetp_mtp += 1
        if signalps['signalp_version'] == 'TargetP-2.0' and signalp_marker in ['OTHER', 'noTP']: # count "other" predictions only if TargetP is used
            targetp_other += 1
        if include_ppc_prediction:
            results.append((record.name,
                            signalp_marker,
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            signalps['signalp_version'] + ': ' + signalp_marker))
        else:
            results.append((record.name,
                            signalp_marker,
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            signalps['signalp_version'] + ': ' + signalp_marker))
    else:  # Didn't recognize the contents of the 'conclusion' column.
        if signalps['signalp_version'] in ('SignalP-5.0', 'TargetP-2.0'):
            raise ValueError(f'Did not recognize the contents of the Prediction column: "{signalp_marker}" in {record.name}.')
        else:
            raise ValueError(f'Did not recognize the contents of the "D Conclusion" column: "{signalp_marker}" in {record.name}.')


# Construct file output
if simple_score_cutoff:
    prediction_header = f'ASAFind {VERSION} Prediction, score threshold = {simple_score_cutoff}, no FWYL check'
    prediction_cutoff_summary = f'{simple_score_cutoff} (custom threshold), no FWYL check'
else:
    if reproduce_ASAFind_1:
        prediction_header = f'ASAFind {VERSION} Prediction, reproducing ASAFind 1, score threshold = {default_score_cutoff}, with FWYL check'
        prediction_cutoff_summary = f'{default_score_cutoff} (non-normalized default threshold), with an FWYL check'
    else:
        prediction_header = f'ASAFind {VERSION} Prediction, score threshold = {default_score_cutoff}, with FWYL check'
        prediction_cutoff_summary = f'{default_score_cutoff} (normalized default threshold), with an FWYL check'
    
if custom_scoring_table:
    score_header = 'ASAFind 20aa score, custom scoring table'
    scoring_table_summary = args.score_table_file
else:
    if reproduce_ASAFind_1:
        score_header = 'ASAFind 20aa score, default scoring table, no small sample correction'
        scoring_table_summary = 'the default scoring table  without small sample correction'
    else:
        score_header = 'ASAFind 20aa score, default scoring table with small sample correction'
        scoring_table_summary = 'the default scoring table with small sample correction'
        
if include_ppc_prediction:
    prediction_header += ', PPC protein prediction included'
    prediction_cutoff_summary += ', PPC protein prediction included'
    if args.score_cutoff_ppc:
        prediction_header += f', custom ppc score cut-off: {score_cutoff_ppc}'
        prediction_cutoff_summary += f', custom ppc score cut-off: {score_cutoff_ppc}'
    else:
        prediction_header += f', default ppc score cut-off: {score_cutoff_ppc}'
        prediction_cutoff_summary += f', default ppc score cut-off: {score_cutoff_ppc}'
    if custom_scoring_table_ppc:
        score_header_ppc = 'PPC 20aa score, custom scoring table'
        scoring_table_summary_ppc = args.score_table_file_ppc
    else:
        score_header_ppc = 'PPC 20aa score, default scoring table'
        scoring_table_summary_ppc = 'default PPC scoring table'
    headers = ('Identifier',
                signalps['signalp_version'],
                'ASAFind cleavage position',
                'ASAFind cleavage site sequence',
                f'ASAFind/{signalps["signalp_version"]} cleavage site offset',
                f'{score_header}',
                f'{score_header_ppc}',
                f'{prediction_header}')
else:
    headers = ('Identifier',
                signalps['signalp_version'],
                'ASAFind cleavage position',
                'ASAFind cleavage site sequence',
                f'ASAFind/{signalps["signalp_version"]} cleavage site offset',
                f'{score_header}',
                f'{prediction_header}')

out_handle.write('\t'.join(headers) + '\n')

for result in results:
    # convert every item in list to string, and round floats to 9 decimal places.
    result = [f'{x:.9f}' if type(x) is float else str(x) for x in result]
    out_handle.write('\t'.join(result) + '\n')


# Construct terminal output
chloroplast_not_targeted = signalp_positive - chloroplast_targeted - ppc_targeted - skipped_proteins

print(f'''########################################
This is ASAFind {VERSION}.
We detected a {signalps['signalp_version']} input file.

Proteins were scored against {scoring_table_summary}
The maximum possible non-normalized transit peptide score for this table is {maximum_tscore}''')
if reproduce_ASAFind_1:
    print('Scores were not normalized to the maximum score')
else:
    print('Scores were normalized to the maximum score')

print(f'''The protein scoring threshold for high confidence predictions was {prediction_cutoff_summary}''')


print(f'''You submitted {total_proteins} proteins
{total_proteins-signalp_positive} of your proteins had no signal peptide detected''')
if signalps['signalp_version'] == 'TargetP-2.0':
    print(f'''    {signalps['signalp_version']} predicted a mitochondrial transit peptide in {targetp_mtp} proteins, and location = "OTHER" for {targetp_other} proteins''')
print(f'''{signalp_positive} of your proteins have a signal peptide''')
print(f'''    {chloroplast_targeted} of these were predicted to go to the plastid''')
if not simple_score_cutoff:
    print(f'''        {chloroplast_targeted_high} of these were predicted with high confidence
        {chloroplast_targeted_low} of these were predicted with with low confidence''')
if include_ppc_prediction:
    print(f'''    {ppc_targeted} of these were predicted to go to the periplastidic compartment''')
if include_ppc_prediction:
    print(f'''    {chloroplast_not_targeted} of these were predicted not to go to the plastid or PPC''')
else:
    print(f'''    {chloroplast_not_targeted} of these were predicted not to go to the plastid''')
if skipped_proteins:
    print(f'    {skipped_proteins} of the proteins could not be analyzed because the length is too short')



print('''\nCitation: If you use ASAFind in your research please cite our publication (Gruber et al., 2025, doi:10.1111/tpj.70138) as well as the appropriate publications for SignalP or TargetP, you used ''')

if 'SignalP-5' in signalps['signalp_version']:
    print('SignalP 5.0, (Almagro Armenteros et al. 2019, doi: 10.1038/s41587-019-0036-z)')
elif 'SignalP-4' in signalps['signalp_version']:
    print('SignalP 4, (Petersen et al., 2011, doi: 10.1038/nmeth.1701)')
elif 'SignalP-3' in signalps['signalp_version']:
    print('SignalP 3.0, (Bendtsen et al., 2004, doi: 10.1016/j.jmb.2004.05.028)')
elif 'TargetP-2' in signalps['signalp_version']:
    print('TargetP 2.0, (Almagro Armenteros et al., 2019, doi: 10.26508/lsa.201900429)')

if web_output:
    print('</PRE>')
else:
    print(f'\nWrote results to {out_file}')
####close log file####
sys.stdout = orig_stdout
log_file.close()
######end#####

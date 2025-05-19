#! /usr/bin/env python
#
# This work is copyrhight 2025 Ansgar Gruber and Marta Vohnoutová, University of South Bohemia,
# and licenced under CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International) licence,
# available at: https://creativecommons.org/licenses/by-sa/4.0/
#
# For updates see our repository on GitHub: https://github.com/ASAFind/ASAFind-2
# ASAFind web service: https://asafind.jcu.cz/
# Publication: https://doi.org/10.1111/tpj.70138 (please cite if you publish results of the scripts, or derivative work)
# 
# Contacts:
# Ansgar Gruber <ansgar.gruber@paru.cas.cz>
# Marta Vohnoutová <mvohnoutova@prf.jcu.cz>
#



VERSION = '2.0'
PROJ = 'ASAFind_graphics'

#import part
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import pickle
#import pickle5 as pickle

pd.options.mode.chained_assignment = None  # default='warn'
import os, sys
import os.path
import numpy as np
import copy
import csv
import shutil


from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio import Entrez
import logomaker as lm
import inspect
from collections import defaultdict
import argparse

import math

from pathlib import Path

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib import font_manager, cm
import matplotlib.ticker as ticker

import logomaker as lm
from PIL import Image, ImageEnhance
from copy import copy
import subprocess

from fill_constants import *


# imports for logomaker images

print('In S0 program inside')

if sys.version_info < (3, 10):
    print(inspect.cleandoc('''\n\nWARNING: This version of ASAFind was developed on Python 3.10,
        it is not tested on Python version:
        {}.{}\n\n'''.format(sys.version_info[0], sys.version_info[1])))

#### Collect Input ####
#######################
parser = argparse.ArgumentParser(description=f'''
{PROJ} {VERSION}

Takes a Fasta and companion TargetP v.2.0 short format file
as input, with the complete TargetP header (two lines starting with '#'). Some versions of SignalP
truncate the sequence names. SignalP-3.0 to 20 characters, and 4.0, 4.1 to 58 characters. Therefore,
ASAFind only considers the first corresponding characters of the fasta name (and the first 90 in the
case of TargetP 2.0), which must be unique within the file. Parts of the fasta name
after that character are ignored. Additionally, the fasta name may not contain a '-' or '|'. This
requirement is because SignalP converts characters in sequence names (e.g. '-' is changed to '_').
ASAFind requires at least 7 aa upstream and 22 aa downstream of the cleavage site suggested by 
SignalP. The output of this script is a tab delimited table. 
Python >= 3.10 required.''',
formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-f', '--fasta_file', help='Specify the input fasta FILE.',
                    required=True)

parser.add_argument('-p', '--signalp_file',
                    help='Specify the input TargetP FILE.', required=True)

parser.add_argument('-s', '--simple_score_cutoff', default=None, type=float,
                    help='Optionally, specify an explicit score cutoff, rather than '
                    'using ASAFind\'s default algorithm, not compatible with option -v1. The score given here will not be normalized and therefore should be obtained form a distribution of normalized scores.')

parser.add_argument('-t', '--fasta_file_with_motifs',
                    help='Optionally, specify a custom scoring table. The scoring table will be normalized with the maximum score, which allows for processing of non-normalized as well as normalized scoring tables.',
                    required=False)

# this file will be automaticall in temp directory as a part of final output file
#parser.add_argument('-o', '--out_file', help='Specify the path and name of the output file '
                    #'you wish to create. Default will be the same as the fasta_file, but with a '
                    #'".tab" suffix.',
                    #required=True)

parser.add_argument('-w', '--web_output', action='store_true',
                    help='Format output for web display. This is mostly useful when called by a '
                    'web app.')

parser.add_argument('-v1', '--reproduce_ASAFind_1', action='store_true',
                    help='Reproduce ASAFind 1.x scores and results (non-normalized scores, if no custom scoring table is specified, the original default scoring table generated without small sample size correction will be used, not compatible with option -s).')

parser.add_argument('-ppc', '--include_ppc_prediction', action='store_true',
                    help='Include prediction of proteins that might be targeted to the periplastidic compartment.')

parser.add_argument('-s_ppc', '--score_cutoff_ppc', default=None, type=float,
                    help='Optionally, specify an explicit score cutoff for the ppc protein prediction, if given, ppc protein prediction will be included. The score given here will not be normalized and therefore should be obtained form a distribution of normalized scores.')

parser.add_argument('-t_ppc', '--fasta_file_with_motifs_ppc',
                    help='Optionally, specify a custom scoring table for the ppc protein prediction, if given, ppc protein prediction will be included. The scoring table will be normalized with the maximum score, which allows for processing of non-normalized as well as normalized scoring tables.')

parser.add_argument('-l', '--logomaker', action='store_true',
                    help='If defined - the program will generate also the logomaker pictures in .png and .svg formats. They will be include into the output compressed package.')

parser.add_argument('-my_org', '--my_organism',
                    help='Specify the name of organism.',
                    required=False)

parser.add_argument('-v', '--version', action='version', version=f'{PROJ} {VERSION}')


args = parser.parse_args()
my_args = vars(parser.parse_args())
#print(my_args)

####create temp directory for temporary files and output directory for output files
where_are_programs = os.getcwd() + '/'
temp_path = f'{os.getcwd()}/temp/'  # check it
output_path = f'{os.getcwd()}/output/' 


def recreate_directory(directory_path):
    # Check if the directory exists
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)  # Remove the entire directory
    os.makedirs(directory_path)  # Recreate the directory
    print(f"Recreated the directory: {directory_path}")

recreate_directory(temp_path)
recreate_directory(output_path)
#######################import of paths and constants #####################
#########################constants#####use before Logomaker###########################################
allowed_graphics=['Plastid, high confidence', 'Plastid, low confidence','PPC','Not plastid, not ppc, signal peptide identified','Not plastid, signal peptide identified']
#########################

fasta_all_path = args.fasta_file

(fasta_file_path, fasta_file_name) = os.path.split(fasta_all_path)
(signalp_file_path, signalp_file_name) = os.path.split(signalp_file)

if args.my_organism:  # organism name arg is not obligatory
    my_organism = args.my_organism
else:
    my_organism = "my_organism"

my_logomaker = args.logomaker

output_all_path = f'{output_path}output_{fasta_file_name}.zip'     # name of fasta file is part of the output.zip file name - Marta added
S1_out_file  = f'{temp_path}{fasta_file_name}.tab'  # name of the output file in temp directory

#####open log file####
orig_stdout = sys.stdout
log_file=open(f'{temp_file}log_S0_{fasta_file_name}.txt', 'w')
sys.stdout = log_file
print(f'Log file for S0 program {fasta_file_name}')
###################



#### Variables and Names ####
#############################

# Figure out some names and paths
if args.logomaker:
    logomaker_yes_no = True
else:
    logomaker_yes_no = False


fasta_file = os.path.abspath(args.fasta_file)
#signalp_file = os.path.abspath(args.signalp_file)

(fasta_file_path, fasta_file_name) = os.path.split(fasta_file)
(fasta_file_base_name, fasta_file_ext) = os.path.splitext(fasta_file_name)

#########define variables from args#####????#####    
simple_score_cutoff = args.simple_score_cutoff
score_cutoff_ppc = args.score_cutoff_ppc
include_ppc_prediction = args.include_ppc_prediction
web_output = args.web_output
reproduce_ASAFind_1 = args.reproduce_ASAFind_1
include_ppc_prediction = args.include_ppc_prediction

###############################################################
# Check compatibility of options

if simple_score_cutoff and reproduce_ASAFind_1:
    raise Exception('options -s and -v1 are not campatible')

if args.score_cutoff_ppc:
    include_ppc_prediction = True

if args.fasta_file_with_motifs:
    fasta_file_with_motifs = os.path.abspath(args.fasta_file_with_motifs)
    (fasta_file_with_motifs_path, fasta_file_with_motifs_name) = os.path.split(fasta_file_with_motifs)

if args.fasta_file_with_motifs_ppc:
    fasta_file_with_motifs_ppc = os.path.abspath(args.fasta_file_with_motifs_ppc)
    (fasta_file_with_motifs_ppc_path, fasta_file_with_motifs_ppc_name) = os.path.split(fasta_file_with_motifs_ppc)

if args.fasta_file_with_motifs and args.simple_score_cutoff is False:
    print ('You specified a custom scoring table, but no custom scoring cut-off. Please be aware that a custom scoring table most likely works best with a custom score cut-off.')

if args.fasta_file_with_motifs_ppc and args.score_cutoff_ppc is False:
    print ('You specified a custom scoring table for the ppc protein prediction, but no custom scoring cut-off. Please be aware that a custom scoring table most likely works best with a custom score cut-off.')

if args.include_ppc_prediction:
    ppc=True
else:
    ppc=False

###############################################################

print('fasta',args.fasta_file)
print('signalp file',args.signalp_file)
print('simple score cutoff',args.simple_score_cutoff)
print('fasta file with motifs',args.fasta_file_with_motifs)
print('web output',args.web_output)
print('reproduce ASAFind 1',args.reproduce_ASAFind_1)
print('include ppc prediction',args.include_ppc_prediction)
print('score cutoff ppc',args.score_cutoff_ppc,type(args.score_cutoff_ppc))
print('fasta file with motifs ppc',args.fasta_file_with_motifs_ppc)
print('logomaker',args.logomaker)
print('my organism',args.my_organism)
print('-s parameter',args.simple_score_cutoff,type(args.simple_score_cutoff))

###############################################################
# before logomaker we need to check the fasta file or the pandas table

# desired_columns_logo2 = list('ACDEFGHIKLMNPQRSTVWY')
scoring_background_matrixes = []  # list of scoring matrixes for background
fasta_background_files=[] # list of fasta files for background
if not args.include_ppc_prediction and not args.fasta_file_with_motifs and not args.fasta_file_with_motifs_ppc:
    scoring_background_matrixes.append('diatom_output.tab.pkl') 
    fasta_background_files.append(f'{diatom_fasta_file}')
if args.include_ppc_prediction and not args.fasta_file_with_motifs and not args.fasta_file_with_motifs_ppc:
    scoring_background_matrixes.append('diatom_output.tab.pkl')
    scoring_background_matrixes.append('ppc_output.tab.pkl')  
    fasta_background_files.append(f'{diatom_fasta_file}')
    fasta_background_files.append(f'{ppc_fasta_file}')  
if not args.include_ppc_prediction and args.fasta_file_with_motifs and not args.fasta_file_with_motifs_ppc:
    cmd = ['python', 'S2_score_table_updated.py']  # Base command to call the second program S2_ASAFind.py in case of user own scoring matrixes
    cmd.extend(['-f', args.fasta_file_with_motifs])  
    cmd.extend(['-o', f'{temp_path}{fasta_file_with_motifs_name}.tab'])  # out_file is the name of the output file in temp directory
    custom_scoring_table = True
    # Print the command to debug (optional)
    print("Subprocess S2 command:", cmd)
    # Call the second Python program
    subprocess.run(cmd)
    scoring_background_matrixes.append(f'{temp_path}{fasta_file_with_motifs_name}.tab.pkl')
    fasta_background_files.append(f'{fasta_file_with_motifs}')
if args.include_ppc_prediction and not args.fasta_file_with_motifs and args.fasta_file_with_motifs_ppc:
    scoring_background_matrixes.append('diatom_output.tab.pkl')
    fasta_background_files.append(f'{diatom_fasta_file}')
    cmd = ['python', 'S2_score_table_updated.py']  # Base command to call the second program S2_ASAFind.py in case of user own scoring matrixes
    cmd.extend(['-f', args.fasta_file_with_motifs_ppc])  
    cmd.extend(['-o', f'{temp_path}{fasta_file_with_motifs_ppc_name}.tab'])  # out_file is the name of the output file in temp directory
    custom_scoring_table_ppc = True
    # Print the command to debug (optional)
    print("Subprocess S2 command:", cmd)
    # Call the second Python program
    subprocess.run(cmd)
    scoring_background_matrixes.append(f'{temp_path}{fasta_file_with_motifs_ppc_name}.tab.pkl')
    fasta_background_files.append(f'{fasta_file_with_motifs_ppc}')
if args.include_ppc_prediction and args.fasta_file_with_motifs and args.fasta_file_with_motifs_ppc:
    cmd = ['python', 'S2_score_table_updated.py']  # Base command to call the second program S2_ASAFind.py in case of user own scoring matrixes
    cmd.extend(['-f', args.fasta_file_with_motifs])  
    cmd.extend(['-o', f'{temp_path}{fasta_file_with_motifs_name}.tab'])  # out_file is the name of the output file in temp directory
    custom_scoring_table = True
    # Print the command to debug (optional)
    print("Subprocess S2 command:", cmd)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("STDOUT:", result.stdout.decode())
    print("STDERR:", result.stderr.decode())

    # Call the second Python program
    subprocess.run(cmd)
    scoring_background_matrixes.append(f'{temp_path}{fasta_file_with_motifs_name}.tab.pkl')
    fasta_background_files.append(f'{fasta_file_with_motifs}')
    cmd = ['python', 'S2_score_table_updated.py']  # Base command to call the second program S2_ASAFind.py in case of user own scoring matrixes
    cmd.extend(['-f', args.fasta_file_with_motifs_ppc])  
    cmd.extend(['-o', f'{temp_path}{fasta_file_with_motifs_ppc_name}.tab'])  # out_file is the name of the output file in temp directory
    custom_scoring_table_ppc = True
    # Print the command to debug (optional)
    print("Subprocess S2 command:", cmd)
    # Call the second Python program
    subprocess.run(cmd)
    scoring_background_matrixes.append(f'{temp_path}{fasta_file_with_motifs_ppc_name}.tab.pkl')
    fasta_background_files.append(f'{fasta_file_with_motifs_ppc}')

###############################################################
#### Functions ####
# from S0 call S1 with parameters
##################
# ! python3 S1_ASAFind.py -f /home/ubuntu/asafind_doc/data/examples/example_5.fasta  -ppc -p /home/ubuntu/asafind_doc/data/examples/fasta_example5.fsa.targetp2 > temp/log.txt
# Construct the arguments for the second Python program
# cmd = ['python', 'ASAFind2.py']  # Base command to call the second program S1_ASAFind.py
cmd = ['python', 'S1_ASAFind.py']  # Base command to call the second program S1_ASAFind.py
# Add arguments to the command only if they are not None or True
if args.fasta_file:
    cmd.extend(['--fasta_file', fasta_file])
if args.signalp_file:
    cmd.extend(['--signalp_file', args.signalp_file])
if args.simple_score_cutoff:
    cmd.extend(['--simple_score_cutoff', str(args.simple_score_cutoff)])
if args.fasta_file_with_motifs:
    cmd.extend(['--score_table_file', f'{temp_path}{fasta_file_with_motifs_name}.tab'])  # score_table_file is the name of the output file in temp directory
if args.fasta_file_with_motifs_ppc:
    cmd.extend(['--score_table_file_ppc', f'{temp_path}{fasta_file_with_motifs_ppc_name}.tab'])  # score_table_file_ppc is the name of the output file in temp directory
if args.web_output:
    cmd.extend(['--web_output', args.web_output])
if args.reproduce_ASAFind_1:
    cmd.extend(['--reproduce_ASAFind_1', args.reproduce_ASAFind_1])
if include_ppc_prediction:
    cmd.extend(['--include_ppc_prediction'])
if args.score_cutoff_ppc:
    cmd.extend(['--score_cutoff_ppc'])

cmd.extend(['--out_file', S1_out_file])  # S1_out_file is the name of the output file in temp directory

# Print the command to debug (optional)
print("Subprocess command:", cmd)

# Call the second Python program
subprocess.run(cmd)

#########create pandas tables ############################
# from S1
def clean_aa_sequence(seq_record):
    '''Takes a Bio.Seq object. Returns a Bio.Seq.Seq object. Checks the validity of a protein
    aa sequence. 1) Must be entirely composed of ACDEFGHIKLMNPQRSTVWYBZJUOX* 2) Replace BZJUO
    with X'''
    valid_letters = 'ACDEFGHIKLMNPQRSTVWYBZJUOX*'    
    replace_these_letters = 'BZJUO-'
    seq = seq_record.upper()
    seq_mutable = MutableSeq(str(seq))
    i = 0
    for aa in seq:
        if aa not in valid_letters + replace_these_letters:
            raise ValueError(f'Letters must be one of {valid_letters}, found character: "{aa}"')
        if aa in replace_these_letters:
            seq_mutable[i] = 'X'
            snippet = seq[max(0,i-10):i+10]
            #print(f'WARNING: Replaced "{aa}" in "{snippet}" with an "X" in {seq_record.name}')
        i += 1
    #print(str(seq_mutable))
    return  str(seq_mutable)

print(S1_out_file)
##########df_signalp_table#######################
df_signalp_file = pd.read_csv(f'{S1_out_file}', sep='\t') # targetp output recalculated by S1
df_signalp_file.columns.values[0] = 'Hit ID'
df_signalp_file.columns.values[2] = 'Cleavage position'
df_signalp_file.columns.values[6] = 'PPC prediction'

##########df_cleavage_table#######################
columns=['Organism','Protein ID','Hit ID','Full sequence','Motif','Cleavage position','PPC prediction','TargetP-2.0']
df_cleavage_table=pd.DataFrame(columns=columns)
protein_id = 0             # in case of unknown fasta file no parsing of header !!! check parsing
for rec in SeqIO.parse(f"{fasta_all_path}", "fasta"): 
    protein_id += 1
    #print(rec.seq)
    df_cleavage_table.loc[len(df_cleavage_table.index)]= [my_organism,str(protein_id),rec.id[:90],clean_aa_sequence(rec.seq) ,None,None,None, None] 

##map df_signalp_file to df_cleavage_table##Cleavage position, PPC prediction, TargetP-2.0 ############################
# Create a dictionary mapping from Hit ID to Cleavage position in the second table (df2)
cleavage_position_map = df_signalp_file.set_index('Hit ID')['Cleavage position'].to_dict()
# Map the "Cleavage position" from df2 to df1 based on "Hit ID"
df_cleavage_table['Cleavage position'] = df_cleavage_table['Hit ID'].map(cleavage_position_map)

# Create a dictionary mapping from Hit ID to Cleavage position in the second table (df2)
cleavage_position_map = df_signalp_file.set_index('Hit ID')['PPC prediction'].to_dict()
# Map the "Cleavage position" from df2 to df1 based on "Hit ID"
df_cleavage_table['PPC prediction'] = df_cleavage_table['Hit ID'].map(cleavage_position_map)

# Create a dictionary mapping from Hit ID to Cleavage position in the second table (df2)
cleavage_position_map = df_signalp_file.set_index('Hit ID')['TargetP-2.0'].to_dict()
# Map the "Cleavage position" from df2 to df1 based on "Hit ID"
df_cleavage_table['TargetP-2.0'] = df_cleavage_table['Hit ID'].map(cleavage_position_map)

# get rid of no signal
# from here proceed only if at least one cleavage position is found  ### Marta added
if not df_cleavage_table['Cleavage position'].isna().all():

    # Count non-None/NaN values in column 'Cleavage position'
    count_non_none = df_cleavage_table['Cleavage position'].notna().sum()


    
    print(f"The FASTA file {fasta_file_name} contains {count_non_none} sequence(s) for which a predicted cleavage site was found in \n{signalp_file_name}, therefore, if choosen by a user, graphical output is generated.")
    
    df_cleavage_table = df_cleavage_table[df_cleavage_table['Cleavage position'].notna()]  # get rid of no signal
    df_cleavage_table['Cleavage position'] = df_cleavage_table['Cleavage position'].apply(int)

    #df_cleavage_table['Motif']=df_cleavage_table.apply(lambda row: row['Full sequence'][row['Cleavage position']-5:row['Cleavage position']+20],axis=1)
    df_cleavage_table["Motif"] = df_cleavage_table.apply(lambda row: row["Full sequence"][row["Cleavage position"] - 6:row["Cleavage position"] +19], axis=1)
    #get rid of query Protein ID > 4 but keep unique Hit ID

    count_pr={}
    accepted_hits=[]
    df_cleavage_table  = df_cleavage_table .reset_index()

    for index, row in df_cleavage_table.iterrows():
        #print(row['Protein ID'], row['Hit ID'])
        try:
            count_pr[row['Protein ID']] += 1
        except KeyError:
            count_pr[row['Protein ID']] = 1
                
        if count_pr[row['Protein ID']] < 4 and row['Hit ID'] in accepted_hits:      # Erase not unique Hit ID        
            df_cleavage_table = df_cleavage_table.drop(df_cleavage_table[((df_cleavage_table['Protein ID'] ==  row['Protein ID']) & (df_cleavage_table['Hit ID'] == row['Hit ID']))].index)
            count_pr[row['Protein ID']] -= 1
                
        if count_pr[row['Protein ID']] < 4 and row['Hit ID'] not in accepted_hits:      # Hit ID already not accepted
            accepted_hits.append(row['Hit ID'])
                
        if count_pr[row['Protein ID']] > 3:      # Erase in all cases
            df_cleavage_table = df_cleavage_table.drop(df_cleavage_table[((df_cleavage_table['Protein ID'] ==  row['Protein ID']) & (df_cleavage_table['Hit ID'] == row['Hit ID']))].index)
            count_pr[row['Protein ID']] -= 1
            
    # get rid of duplicities in Hit ID
    df_cleavage_table=df_cleavage_table.drop_duplicates(subset='Hit ID', keep="first") # just to be sure

    # get rid of duplicities in motifs
    df_cleavage_table=df_cleavage_table.drop_duplicates(subset='Motif', keep="first")

    # get rid of shorter motifs
    df_cleavage_table=df_cleavage_table[df_cleavage_table['Motif'].str.len() == 25]

    arr_motifs = df_cleavage_table["Motif"].to_numpy()
    
    # Sequence logo
    l_motifs_foreground=[str(i) for i in arr_motifs]

    print(df_cleavage_table)
    ################ logomaker#####only if switch --logomaker or -l is added ################
    if logomaker_yes_no:  # argument -l activated Marta added
        # create color scheme
        color_scheme = {
            'ACFGILMPVWY' : [0, 0, 0], # hydrophobic
            'NQST' : [0, 0.7, 0],  # hydrophylic
            'HKR' : [0, 0, 1],  # basic
            'DE' : [1, 0, 0]  # acidic
        }

        # Function to determine the color based on the amino acid
        def get_color_for_amino_acid(amino_acid):
            for key, color in color_scheme.items():
                if amino_acid in key:
                    return color
            return [0, 0, 0]  # Default to black if not in any group

        def what_color(amino_acit):
            for i in color_scheme.keys():
                if amino_acit in i:
                    color=color_scheme[i]
                    return color
            return [0.5, 0.5, 0.5]    
        
        # for bottom headline
        def resize_width_only(image_path, new_width):
            # Open the image
            image = Image.open(image_path)
            
            # Get the original height
            original_width, original_height = image.size
            
            # Resize the image, keeping the original height
            resized_image = image.resize((new_width, original_height))
            
            return resized_image

        def get_concat(all_images, im_width, im_height, bottom_image_path):
            # Load all images and resize them to have the same width
            images = [Image.open(all_images[i]).resize((im_width, im_height)) for i in range(len(all_images))]
            
            # Load and resize the additional image for the bottom
            #bottom_image = Image.open(bottom_image_path).resize((im_width, im_height))
            bottom_image = resize_width_only(bottom_image_path, im_width)
            bottom_image_height = bottom_image.height
            
            # Calculate total height to accommodate the additional image
            total_height = (len(images) * im_height) + bottom_image_height
            
            # Create a new blank image with the updated height
            dst = Image.new('RGB', (im_width, total_height))
            
            # Paste each resized image in the list
            for i in range(len(images)):
                dst.paste(images[i], (0, i * im_height))
            
            # Paste the resized additional image at the bottom
            dst.paste(bottom_image, (0, len(images) * im_height))
            
            return dst

        l_motifs_foreground=list(df_cleavage_table['Motif'])   # foreground motifs

        # Define a function to format y-axis as percentages
        def to_percent(y, position):
            """Convert normalized y-axis values to percentage format."""
            return f"{y:.0f}%"  # Display as whole percentage

        print(scoring_background_matrixes)

        for file_index in range(len(scoring_background_matrixes)):
            # do not edit PPC either or
            with open(scoring_background_matrixes[file_index], 'rb') as background:  # background the same for all
                    my_diatom_pattern_scoring_matrix1 = pickle.load(background)
            
            my_diatom_pattern_scoring_matrix1=my_diatom_pattern_scoring_matrix1.split('\n')
            
            background_dict = {}
            
            for i in my_diatom_pattern_scoring_matrix1:
                if len(i) > 0:
                    l = i.split('\t')
                    l = [_ for _ in l if _!='']
                    
                    background_dict[l[0]] = [float(j) for j in l[1:]]
            
            background_info_df = pd.DataFrame(background_dict)
            
            # count once used for all 
            input_fasta_background=fasta_background_files[file_index] # background
            l_motifs_background=[]
            fasta_sequences = SeqIO.parse(open(f'{input_fasta_background}'),'fasta')
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                l_motifs_background.append(clean_aa_sequence(sequence))
            
            background_counts_mat = lm.alignment_to_matrix(l_motifs_background)
            
            bc1=background_counts_mat/(background_counts_mat.loc[0].sum()/100) # correction of height of graph

            ### complete the missing columns (if any) in the bc1 matrix Marta added 
            # Find missing columns
            missing_columns = [col for col in desired_columns_logo2 if col not in bc1.columns]
            # Add missing columns filled with zeros
            for col in missing_columns:
                bc1[col] = 0

            # Ensure the column order matches the desired list (optional)
            bc1 = bc1[desired_columns_logo2]
            ###########
            #print(f'background_counts_mat k logo2 {background_counts_mat}')
            #print(f'background_info_df k logo1 {background_info_df}')
            #print(f'l_motifs_background k logo2 {l_motifs_background}')

            if file_index == 0:
                my_background='Plastid scoring'
            else:
                my_background='PPC scoring'
            # whole cycles - both graphs
            
            for protein_iloc in range(len(l_motifs_foreground)):        

                if 'X' not in df_cleavage_table.iloc[protein_iloc]['Motif'] and '*' not in df_cleavage_table.iloc[protein_iloc]['Motif']:  # !!!! get rid of illegal X to logomaker
                    # load ww information matrix - skip zero
                    rel_position = [i for i in range(-5,0)]+[i for i in range(1,21)]
                
                    #fig = plt.figure(111)
                    #ax = fig.add_subplot(111)
                
                    # create Logo object
                    #plt.subplot(2, 1, 1)  # 2 rows, 1 column, plot 1
                    ww_logo1 = lm.Logo(background_info_df,
                                font_name='Arial',
                                color_scheme='lightgray',
                                vpad=.1,
                                width=.8)
                    
                    for seq_no in range(len(l_motifs_foreground[protein_iloc])):
                        ww_logo1.style_single_glyph(c=l_motifs_foreground[protein_iloc][seq_no], p=seq_no, color=what_color(l_motifs_foreground[protein_iloc][seq_no]))
                
                    # style using Logo methods
                    ww_logo1.ax.set_xlabel("Amino acids")
                    #plt.xticks(ticks = tickvalues ,labels = labellist, rotation = 'vertical')
                    plt.xticks(ticks= range(len(l_motifs_foreground[protein_iloc])),labels=list(l_motifs_foreground[protein_iloc]))
                    #ww_logo1.style_xticks(anchor=0, spacing=1, rotation=45)
                    ww_logo1.highlight_position(p=5, color='gold', alpha=.5)   # in 2.0 version changed to 5
                    
                    # style using Axes methods
                    try:
                        ww_logo1.ax.set_title(f"{df_cleavage_table.iloc[protein_iloc]['Hit ID']} : {df_cleavage_table.iloc[protein_iloc]['PPC prediction']}")
                    except KeyError:
                        ww_logo1.ax.set_title(f"Protein")
                        
                    ww_logo1.ax.set_ylabel('bits')
                
                    #ww_logo1.ax.set_xlim([-1, len(counts_mat)])
                    try:
                        cleavage = df_cleavage_table.iloc[protein_iloc]['Cleavage position']
                    except KeyError:
                        cleavage=0
                        
                    # rainbow_text(0.0, -100.00, words, colors, size=12)
                        
                    sec_x_ticks=[i for i in range(cleavage-5,cleavage+20)]
                    
                    color = 'black'
                    ww_logo1.ax2 = ww_logo1.ax.twiny() 
                    
                    ww_logo1.ax2.xaxis.set_ticks_position("bottom")
                    ww_logo1.ax2.xaxis.set_label_position("bottom")
                    
                    ww_logo1.ax2.tick_params(axis='x', labelcolor=color)
                    ww_logo1.ax2.set_xticks(ticks= range(0,len(sec_x_ticks)),labels=sec_x_ticks)
                    ww_logo1.ax2.set_xlim(left=-0.5, right=len(sec_x_ticks)-0.5)
                
                    ww_logo1.ax2.spines["bottom"].set_position(("axes", -0.25))
                
                    ww_logo1.ax.text(0.95, 0.95, f"{my_background}", fontsize=12, color='black', ha='right', va='top', transform=ww_logo1.ax.transAxes)
                    
                    try:
                        ww_logo1.ax2.set_xlabel(f"..position reative to predicted cleavage site is {df_cleavage_table.iloc[protein_iloc]['Cleavage position']}..")
                    except KeyError:
                        ww_logo1.ax2.set_xlabel(f"..position reative to predicted cleavage site is 0..")  
                
                    for tick in  ww_logo1.ax.get_xticklabels():     
                        amino_acid = tick.get_text()
                        color = get_color_for_amino_acid(amino_acid)
                        tick.set_color(color)


                    for tick in ww_logo1.ax2.get_xticklabels():
                        tick_label = tick.get_text()
                        if tick_label.strip():  # Check if tick label is not empty
                            if int(tick_label) == df_cleavage_table.iloc[protein_iloc]['Cleavage position']:
                                tick.set_color('r')
                                tick.set_fontsize(16)
                                tick.set_fontweight
                    
                        
                    #ww_logo1.ax.set_xlim([-1, len(counts_mat)])
                    #plt_logo1 = plt
                    if file_index == 0:
                        try:
                            plt.savefig(f"{temp_path}logo1_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png", bbox_inches='tight')
                        except:
                            plt.savefig(f'{temp_path}logo1_{my_organism}_some_protein.png', bbox_inches='tight')
                        
                        
                        try:
                            plt.savefig(f"{temp_path}logo1_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.svg", bbox_inches='tight', format="svg")
                        except:
                            plt.savefig(f'{temp_path}logo1_{my_organism}_some_protein.svg', bbox_inches='tight',format="svg")
                    else:
                        try:
                            plt.savefig(f"{temp_path}logo1_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png", bbox_inches='tight')
                        except:
                            plt.savefig(f'{temp_path}logo1_PPC_{my_organism}_some_protein.png', bbox_inches='tight')
                        
                        
                        try:
                            plt.savefig(f"{temp_path}logo1_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.svg", bbox_inches='tight', format="svg")
                        except:
                            plt.savefig(f'{temp_path}logo1_PPC_{my_organism}_some_protein.svg', bbox_inches='tight',format="svg")

                    plt.close('all')
                    #plt.show()
                    
                    # LOGO2 ######################################
                    # load ww information matrix - frequency
                    #fig = plt.figure(111)
                    #ax = fig.add_subplot(111)
                    #plt.subplot(2, 1, 2)  # 2 rows, 1 column, plot 1
                    # create Logo object
                    # Normalize the counts to a 0-100% scale for logo2
                    #max_count = max(counts)  # Find the maximum value in counts for normalization
                    #normalized_counts = [100 * count / max_count for count in counts]  # Scale to 100%
                    
                    ww_logo2 = lm.Logo(bc1,
                                        font_name='Arial',
                                        color_scheme='lightgray',
                                        vpad=.1,
                                        width=.8)
                    
                    for seq_no in range(len(l_motifs_foreground[protein_iloc])):
                        ww_logo2.style_single_glyph(c=l_motifs_foreground[protein_iloc][seq_no], p=seq_no, color=what_color(l_motifs_foreground[protein_iloc][seq_no]))
                
                    # style using Logo methods
                    ww_logo2.ax.set_xlabel("Amino acids")
                    #plt.xticks(ticks = tickvalues ,labels = labellist, rotation = 'vertical')
                    plt.xticks(ticks= range(len(l_motifs_foreground[protein_iloc])),labels=list(l_motifs_foreground[protein_iloc]))
                
                    
                    #ww_logo2.style_xticks(anchor=0, spacing=1, rotation=45)
                    ww_logo2.highlight_position(p=5, color='gold', alpha=.5)   # in 2.0 version changed to 5
                    
                    # style using Axes methods
                    try:
                        ww_logo2.ax.set_title(f"{df_cleavage_table.iloc[protein_iloc]['Hit ID']} : {df_cleavage_table.iloc[protein_iloc]['PPC prediction']}\n{my_background}")
                    except KeyError:
                        ww_logo2.ax.set_title(f"Protein")
                
                    # logo2 y axis
                    # Format y-axis as percentages
                
                    #ww_logo2.ax.set_xlim([-1, len(counts_mat)])
                    try:
                        cleavage = df_cleavage_table.iloc[protein_iloc]['Cleavage position']
                    except KeyError:
                        cleavage=0
                        
                    # rainbow_text(0.0, -100.00, words, colors, size=12)
                        
                    sec_x_ticks=[i for i in range(cleavage-5,cleavage+20)]
                    
                
                    color = 'black'
                    ww_logo2.ax2 = ww_logo2.ax.twiny() 
                    
                    ww_logo2.ax2.xaxis.set_ticks_position("bottom")
                    ww_logo2.ax2.xaxis.set_label_position("bottom")
                    
                    ww_logo2.ax2.tick_params(axis='x', labelcolor=color)
                    ww_logo2.ax2.set_xticks(ticks= range(0,len(sec_x_ticks)),labels=sec_x_ticks)
                    ww_logo2.ax2.set_xlim(left=-0.5, right=len(sec_x_ticks)-0.5)
                    
                    ww_logo2.ax2.spines["bottom"].set_position(("axes", -0.25))
                    
                    # set position ax2 axis
                    #ww_logo2.ax2.set_xticklabels(tick_function(new_tick_locations))
                
                    # y ticks
                    ww_logo2.ax.set_ylabel('Frequency')
                    ww_logo2.ax.set_yticks([0, 25, 50, 75, 100])  # Define fixed tick positions
                
                    
                    # Normalize the y-axis to represent 0-100%
                    ww_logo2.ax.set_ylim(0, 100)  # Set y-axis from 0 to 100%
                    # Apply the percentage formatter to the y-axis
                    ww_logo2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(to_percent))
                    # Apply the percentage formatter to display y-axis in percentages
                        
                    try:
                        ww_logo2.ax2.set_xlabel(f"..position in sequence..cleavage position is {df_cleavage_table.iloc[protein_iloc]['Cleavage position']}..")
                    except KeyError:
                        ww_logo2.ax2.set_xlabel(f"..position in sequence..cleavage position is 0..")  
                
                    for tick in  ww_logo2.ax.get_xticklabels():     
                        amino_acid = tick.get_text()
                        color = get_color_for_amino_acid(amino_acid)
                        tick.set_color(color)
                        
                    for tick in  ww_logo2.ax2.get_xticklabels():
                        tick_label = tick.get_text()
                        if tick_label.strip():  # Check if tick label is not empty
                            if int(tick_label) == df_cleavage_table.iloc[protein_iloc]['Cleavage position']:
                                tick.set_color('r')
                                tick.set_fontsize(16)
                                tick.set_fontweight
                        
                            
                    #bottom_title="colour code for observed residues NQST (hydrophilic), HKR (basic), DE (acidic), ACFGILMPVWY (other)"
                    
                    #plt.title(bottom_title, y=-0.65)
                    
                    #print("{where_are_data}{organism}_logomaker/logo2_{organism}_{df_cleavage_table.iloc[protein_iloc]['Protein ID']}_{df_cleavage_table.iloc[protein_iloc]['Sequence order']}.png")
                    plt_logo2 = plt 
                    if file_index == 0:
                        try:
                            plt.savefig(f"{temp_path}logo2_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png", bbox_inches='tight')
                        except:
                            plt.savefig(f'{temp_path}logo2_{my_organism}_some_protein.png', bbox_inches='tight')
                        
                        
                        try:
                            plt.savefig(f"{temp_path}logo2_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.svg", bbox_inches='tight', format="svg")
                        except:
                            plt.savefig(f'{temp_path}logo2_{my_organism}_some_protein.svg', bbox_inches='tight',format="svg")
                    else:
                        try:
                            plt.savefig(f"{temp_path}logo2_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png", bbox_inches='tight')
                        except:
                            plt.savefig(f'{temp_path}logo2_PPC_{my_organism}_some_protein.png', bbox_inches='tight')
                        
                        
                        try:
                            plt.savefig(f"{temp_path}logo2_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.svg", bbox_inches='tight', format="svg")
                        except:
                            plt.savefig(f'{temp_path}logo2_PPC_{my_organism}_some_protein.svg', bbox_inches='tight',format="svg")
                        
                    
                    # put two images together
                    if file_index == 0:
                        all_images = [f"{temp_path}logo1_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png",
                                    f"{temp_path}logo2_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png"]
                    else:
                        all_images = [f"{temp_path}logo1_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png",
                                    f"{temp_path}logo2_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png"]
                    
                    
                    im = Image.open(all_images[0])
                    im_width=im.size[0]
                    im_height=im.size[1]
                    im.close()
                    bottom_image_path = f"colour_code.png"
                    
                    #result= get_concat(all_images,im_width,im_height).save(f"{where_are_data}{organism}_logomaker/logo_{organism}_{df_cleavage_table.iloc[protein_iloc]['Protein ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png")
                    result= get_concat(all_images,im_width,im_height,  bottom_image_path)
                    
                    #result.show()  # To display the image
                    if file_index == 0:
                        result.save(f"{temp_path}logo_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Protein ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png")  # To save the final combined image
                    else:
                        result.save(f"{temp_path}logo_PPC_{my_organism}_{df_cleavage_table.iloc[protein_iloc]['Protein ID']}_{df_cleavage_table.iloc[protein_iloc]['Hit ID']}.png")  # To save the final combined image
                    
                    plt.close('all')
                    #break
                    #plt.show()
                else:
                    log_file = open('temp/log.txt','a')
                    log_file.write(f"Not acceptable marks 'X' or '*' in Motif {df_cleavage_table.iloc[protein_iloc]['Motif']}\n")
                    log_file.close()

else:
    
    print(f"The FASTA file {fasta_file_name} does not contain any sequence for which a predicted cleavage site was found in \n{signalp_file_name}, therefore, no graphical output is generated.")
    


####### output zip file ############################
# zip files for user download

# importing required modules 
from zipfile import ZipFile 
import os 

def get_all_file_paths(directory): 

    # initializing empty file paths list 
    file_paths = [] 
  
    # crawling through directory and subdirectories 
    for root, directories, files in os.walk(directory): 
        for filename in files: 
            # join the two strings in order to form the full filepath. 
            filepath = os.path.join(root, filename) 
            file_paths.append(filepath) 
  
    # returning all file paths 
    return file_paths         
  


# path to folder which needs to be zipped 
directory = f"{temp_path}"
    
# calling function to get all file paths in the directory 
file_paths = get_all_file_paths(directory) 
    
# printing the list of all files to be zipped  
# writing files to a zipfile 

with ZipFile(f'{output_all_path}','a') as zip: 
        
    # writing each file one by one 
    for file in file_paths: 
            
        filename = os.path.basename(file)            
        zip.write(file, arcname=filename) 
  
print('All files zipped successfully!')         

######remove temp files#####################
shutil.rmtree(f'{temp_path}')

####close log file####
sys.stdout = orig_stdout
log_file.close()
######end#####
# the end





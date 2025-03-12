# define constants for S1_ASAFind
# 
# Marta Vohnoutova
import os
aminoacids_abbr='ACDEFGHIKLMNPQRSTVWYX'
my_home = os.path.expanduser('~')

# edit this part
# no need now
#Entrez_email = "mvohnoutova@jcu.cz" # Always tell NCBI who you are
#root_dir=f'{my_home}/asafind_line_command/'
#where_are_programs=f'{my_home}/asafind_line_command/'
root_dir=os.getcwd()
where_are_programs=os.getcwd()
organism="my_organism"
temp_file=f'{root_dir}/temp/'
diatom_fasta_file=f'{where_are_programs}/diatom_scoring_matrix.fasta'
ppc_fasta_file=f'{where_are_programs}/ppc_scoring_matrix.fasta'

# no need, targetp file is obligatory input for S1 program
#targetp_path = os.path.expanduser('~') +"/asafind_doc/TargetP/targetp-2.0/bin/targetp"  # check it




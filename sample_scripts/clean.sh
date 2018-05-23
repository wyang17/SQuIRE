######################################################
# Bash script clean.sh to run the Clean step of SQuIRE with input arguments
# from input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=3G
#$ -l h_vmem=3G

#$ -l h_fsize=500G
#$ -m e
#$ -M squire@email.com


#Load arguments
echo 'Loading arguments'
argument_file=$1
. $argument_file

# Set up environment and modules for SQuIRE
echo 'Setting up environment'
source activate $virtual_env



# Run SQuIRE Build
echo 'Running Clean'

if [ -z $repeatmasker_file ]
then
  squire Clean --build $build --fetch_folder $fetch_folder --clean_folder $clean_folder --extra $extra $verbosity 

elif
then
  squire Clean --rmsk $repeatmasker_file --fetch_folder $fetch_folder --clean_folder $clean_folder --extra $non_reference $verbosity 
fi

echo 'Clean Complete on' `date`

# clean.sh

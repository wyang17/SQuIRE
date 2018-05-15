######################################################
# Bash script fetch.sh to run the Fetch step of SQuIRE with input arguments
# Input file: arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=7G
#$ -l h_vmem=7G
#$ -pe local 8
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
echo 'Loading modules'
pthreads=8

# Run SQuIRE Fetch
echo 'Running Fetch'
squire Fetch --build $build --fetch_folder $fetch_folder --fasta --rmsk  --chrom_info  --index  --gene --pthreads $pthreads $verbosity

echo 'Fetch Complete on' `date`

# fetch.sh

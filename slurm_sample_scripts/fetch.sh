#!/bin/bash -l
######################################################
# Bash script fetch.sh to run the Fetch step of SQuIRE with input arguments
# Input file: arguments.sh
# Last update: 2018_01_12
# Tested on version 0.9.9.92
# cpacyna
# ---
# Run with:
# sbatch -D . --export=argument_file='arguments.sh' fetch.sh
######################################################


#SBATCH
#SBATCH --job-name=fetch
#SBATCH --time=24:0:0
#SBATCH --partition=parallel
# number of tasks (processes) per node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=end
#SBATCH --mail-user=EMAIL@EMAIL.EDU



#Load arguments
echo 'Loading arguments'
pwd

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


######################################################
# Bash script count.sh to run individual squire Draw jobs
# using input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=8G
#$ -l h_vmem=8G
#$ -pe local 2
#$ -l h_fsize=500G
#$ -m e
#$ -M squire@email.com

#Load arguments
echo 'Loading arguments'
basefile=$1
argument_file=$2
. $argument_file
pthreads=8

# Set up environment and modules for SQuIRE
source activate $virtual_env

# Run SQuIRE Draw
echo 'Running Draw'
if squire Draw --map_folder $map_folder --draw_folder $draw_folder --name $basefile --normlib $normlib --pthreads $pthreads --strandedness $strandedness $verbosity --build $build --fetch_folder $fetch_folder
then
  echo 'Draw is complete'
fi

# draw.sh

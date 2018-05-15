######################################################
# Bash script call.sh to run squire Call job
# using input file arguments.sh
# Last update: 2018_02_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_fsize=500G
#$ -pe local 8
#$ -m e
#$ -M squire@email.com

#Load arguments
echo 'Loading arguments'
argument_file=$1
. $argument_file
pthreads=8

# Set up environment and modules for SQuIRE
source activate $virtual_env

# Run SQuIRE Call
echo 'Running Call'
if squire Call --group1 $group1 --group2 $group2 --condition1 $condition1 --condition2 $condition2 --projectname $projectname --pthreads $pthreads --output_format $output_format  --call_folder $call_folder $verbosity
then
  echo 'squire Call is complete'
fi

#call.sh

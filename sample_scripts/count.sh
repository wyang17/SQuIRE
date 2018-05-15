######################################################
# Bash script count.sh to run individual squire Count jobs
# using input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=36G
#$ -l h_vmem=36G
#$ -l h_fsize=500G
#$ -m e
#$ -M squire@email.com

#Load arguments
echo 'Loading arguments'
name=$1
argument_file=$2
. $argument_file

# Set up environment and modules for SQuIRE
echo 'Setting up environment'
source activate $virtual_env


# Run SQuIRE Count
echo 'Running Count'
if [ -z $tempfolder ]
then
  temp_folder=$count_folder
fi

if [ -z $EM ]
then
  EM="auto"
fi

squire Count --map_folder $map_folder --clean_folder $clean_folder --count_folder $count_folder --temp_folder $temp_folder --name $name --build $build --strandedness $strandedness --EM $EM $verbosity


echo 'Count Complete on' `date`

# count.sh

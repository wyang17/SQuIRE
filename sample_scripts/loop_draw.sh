######################################################
# Bash script count.sh to run individual squire Count jobs
# using input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=2G
#$ -l h_vmem=2G

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

# Loop through
for file in $map_folder/*.bam
do
basefile=$(basename $file)
basename=${basefile//.bam/}
qsub -cwd -V draw.sh $basefile $argument_file
done

#loop_draw.sh

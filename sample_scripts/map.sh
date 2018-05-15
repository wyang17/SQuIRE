######################################################
# Bash script map.sh to run individual squire Map jobs
# using input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l h_fsize=500G
#$ -pe local 8
#$ -m e
#$ -M cpacyna@jhu.edu

#Load arguments
echo 'Loading arguments'
r1file=$1
r2file=$2
name=$3
argument_file=$4
pthreads=8
. $argument_file

# Set up environment and modules for SQuIRE
echo 'Setting up environment'
source activate $virtual_env

# Run SQuIRE Map
echo 'Running Map'

if [ -z $non_reference ]
then
  non_reference=False
fi

if [[ $r2file != 'False' ]]
then

  if squire Map --read1 $r1file --read2 $r2file --map_folder $map_folder --read_length $read_length --index_folder $fetch_folder --extra $non_reference --pthreads $pthreads --verbosity
  then
    echo $name >> success_map_$projectname.txt
  else
    echo $name >> fail_map_$projectname.txt
  fi

elif [[ $r2file = 'False' ]]
then

  if squire Map --read1 $r1file  --map_folder $map_folder --read_length $read_length --index_folder $fetch_folder --extra $non_reference --pthreads $pthreads --verbosity
  then
    echo $name >> success_map_$projectname.txt
  else
    echo $name >> fail_map_$projectname.txt
  fi

fi

echo 'Map Complete on' `date`

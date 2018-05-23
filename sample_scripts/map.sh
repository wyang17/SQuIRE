######################################################
# Bash script map.sh to run individual squire Map jobs
# using input file arguments.sh
# Last update: 2018_05_21
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
basename=${r1file//}
basename=$2
argument_file=$3
r2file=$4
pthreads=8
. $argument_file

# Set up environment and modules for SQuIRE
echo 'Setting up environment'
source activate $virtual_env

# Run SQuIRE Map
echo 'Running Map'

if [[ $r2file != 'False' ]]
then

  if [ -z $non_reference ]; then
    if squire Map --read1 $r1file --read2 $r2file --map_folder $map_folder --read_length $read_length --fetch_folder $fetch_folder --pthreads $pthreads --build $build --name $basename $verbosity
    then
      echo $basename >> success_map_$projectname.txt
    else
      echo $basename >> fail_map_$projectname.txt
    fi
  else
    if squire Map --read1 $r1file --read2 $r2file --map_folder $map_folder --read_length $read_length --fetch_folder $fetch_folder --extra $non_reference --pthreads $pthreads --build $build --name $basename $verbosity
    then
      echo $basename >> success_map_$projectname.txt
    else
      echo $basename >> fail_map_$projectname.txt
    fi
  fi

elif [[ $r2file = 'False' ]]
then
  if [ -z $non_reference ]; then
    if squire Map --read1 $r1file --map_folder $map_folder --read_length $read_length --fetch_folder $fetch_folder --pthreads $pthreads --build $build --name $basename $verbosity
    then
      echo $basename >> success_map_$projectname.txt
    else
      echo $basename >> fail_map_$projectname.txt
    fi
  else
    if squire Map --read1 $r1file --map_folder $map_folder --read_length $read_length --fetch_folder $fetch_folder --extra $non_reference --pthreads $pthreads --build $build --name $basename $verbosity
    then
      echo $basename >> success_map_$projectname.txt
    else
      echo $basename >> fail_map_$projectname.txt
    fi
  fi

fi

echo 'Map Complete on' `date`

# map.sh

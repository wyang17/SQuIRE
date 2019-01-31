######################################################
# Bash script loop_count.sh to loop through fastq files and send Count jobs
# using input file arguments.sh
# Last update: 2018_01_12
# cpacyna
######################################################

#! /bin/bash
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l h_fsize=500G
#$ -m e
#$ -M cpacyna@jhu.edu


#Load arguments
echo 'Loading arguments'
argument_file=$1
. $argument_file

# Loop through read 1 fastq files
for samplename in $samplenames
do
  if [[ `stat -c %s $count_folder/${samplename}_counts.txt` < 1000 ]]
  then
    qsub -cwd -V count.sh $samplename $argument_file
  else
    echo "Count is already complete for $file"
  fi
done

# loop_count.sh


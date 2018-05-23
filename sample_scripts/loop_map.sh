######################################################
# Bash script loop_map.sh to loop through fastq files and send Map jobs
# using input file arguments.sh
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


# Loop through read 1 fastq files
for samplename in ${samplenames//,/ }
do
for fastq in $fastq_folder/${samplename}*$r1suffix

do
  read1+=${fastq}
done


  # Check if the data is paired or unpaired
  if [ $r2suffix != 'False' ]
  then
  for fastq in $fastq_folder/${samplename}*$r2suffix
    do
    read2+=${fastq}
    done

    # Check if map has already been comleted
    if [[ `stat -c %s $map_folder/${samplename}.bam` < 1000 ]]
    then
      # Send paired map.sh job
      qsub -cwd -V map.sh $read1 $samplename $argument_file $read2
    fi

  else
    # Check if map has already been comleted
    if [[ `stat -c %s $map_folder/$base.bam` < 1000 ]]
    then
      # Send unpaired map.sh job
      read2=False
      qsub -cwd -V map.sh $read1 $samplename $argument_file $read2
    fi
  fi
done


# loop_map.sh

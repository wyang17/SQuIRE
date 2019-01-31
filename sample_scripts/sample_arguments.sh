
############ sample_arguments.sh ############
#Sample arguments file for squire

# Assumes you have four SRR fastq files in directory '/cpacyna/drosophila/fastq'
# Those fastq files should be named in the form 'SRR4016864_1.fastq.gz' and
# 'SRR4016864_2.fastq.gz


# Bash script to save and import arguments for the entire SQuIRE project

## Data information REQUIRED
virtual_env='squire'
fastq_folder=/cpacyna/drosophila/fastq
#Location of fastq files
samplenames=SRR4016864,SRR4016865,SRR4016861,SRR4016860
#Comma separated list of sample names found in Fastq files
build=dm6
#UCSC designation for genome build, eg. 'hg38'; if using custom repeatmasker file (see below), enter custom build name
strandedness=2
#'0' if unstranded, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard
read_length=50
#read length of FASTQ data
projectname=drosrepo
#
r1suffix=_1.fastq.gz
# Common suffix of all r1 fastq files (e.g. for project1_file1_r1.fastq.gz, r1suffix='_r1').  If unpaired data, r1suffix='.fastq'.
r2suffix=_2.fastq.gz
# Common suffix of all r2 fastq files (e.g. for project1_file1_r2.fastq, r2suffix='_r2').  If unpaired data, r2suffix='False'

#ADVANCED
repeatmasker_file=
#filepath of repeatmasker file if using RepeatMasker software output (Leave blank if you want repeatmasker file from UCSC)
non_reference=
#filepath of table of non-reference TEs (eg polymorphic or plasmid TEs) of interest
EM=auto
#desired number of EM iterations other than auto
temp_folder='$TMPDIR'
#Location or variable (such as $TMPDIR) to store intermediate files
group1=SRR4016864,SRR4016865
#Name of basenames of samples in group 1
group2=SRR4016861,SRR4016860
#Name of basenames of samples in group 2
condition1=treated
#Name of condition for group 1 in squire Call
condition2=control
#Name of condition for group 2 in squire Call
output_format=pdf
#Desired output of figures as html or pdf
#normlib="--normlib"
#Uncomment to normalize bedgraphs by library size (default = False) for squire Draw.
verbosity="--verbosity"
#report progress of SQuIRE script
trim3=0

fetch_folder=/cpacyna/squire_fetch
clean_folder=/cpacyna/squire_clean
map_folder=/cpacyna/drosophila/squire_map
count_folder=/cpacyna/drosophila/squire_count
draw_folder=/cpacyna/drosophila/squire_draw
call_folder=/cpacyna/drosophila/squire_call_repo
seek_folder=/cpacyna/drosophila/squire_seek
download_folder=/cpacyna/drosophila/downloads

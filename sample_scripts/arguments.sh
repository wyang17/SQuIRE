############ arguments.sh ############

# Bash script to save and import arguments for the entire SQuIRE project



virtual_env='squire' # Name of virtual environment


## Data information REQUIRED
fastq_folder=
#Location of fastq files

samplenames=
#Space separated list of sample names found in Fastq files

build=
#UCSC designation for genome build, eg. 'hg38'; if using custom repeatmasker file (see below), enter custom build name

strandedness=
#'0' if unstranded, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard

read_length=
#read length of FASTQ data

projectname=
#

r1suffix=
# Common suffix of all r1 fastq files (e.g. for project1_file1_r1.fastq.gz, r1suffix='_r1').  If unpaired data, r1suffix='.fastq'.

r2suffix=
# Common suffix of all r2 fastq files (e.g. for project1_file1_r2.fastq, r2suffix='_r2').  If unpaired data, r2suffix='False'


#ADVANCED

repeatmasker_file=
#filepath of repeatmasker file if using RepeatMasker software output (Leave blank if you want repeatmasker file from UCSC)

non_reference=
#filepath of table of non-reference TEs (eg polymorphic or plasmid TEs) of interest

EM=
#desired number of EM iterations other than auto

temp_folder=
#Location or variable (such as $TMPDIR) to store intermediate files

group1=
#Name of basenames of samples in group 1

group2=
#Name of basenames of samples in group 2

condition1=
#Name of condition for group 1 in squire Call

condition2=
#Name of condition for group 2 in squire Call

output_format=
#Desired output of figures as html or pdf

#normlib="--normlib"
#Uncomment to normalize bedgraphs by library size (default = False) for squire Draw.

verbosity="--verbosity"
#report progress of SQuIRE script



fetch_folder=squire_fetch
clean_folder=squire_clean
map_folder=squire_map
count_folder=squire_count
draw_folder=squire_draw
call_folder=squire_call
seek_folder=squire_seek

#!/bin/env python

#################### MODULES ###################
import sys
import os
import shutil
import subprocess
from subprocess import *
import argparse #module that passes command-line arguments into script
from pkg_resources import get_distribution
__version__ = get_distribution("SQuIRE").version
script_folder=os.path.dirname(os.path.realpath(__file__))
currentWorkingDirectory = os.getcwd()
sys.path.append(currentWorkingDirectory)
sys.path.append(script_folder)
## import the processes to be called
import Build as s1
import Fetch as s2
import Clean as s3
import Map as s4
import Count as s5
import Call as s6
import Draw as s7
import Seek as s8


#from squire import __version__
##################################

def main():
    ## create the top level parser
    parser = argparse.ArgumentParser()
    parser._positionals.title = "SQuIRE Steps"
    parser.add_argument('--version', action="version", version=__version__, help="print SQuIRE version number")
    subparsers = parser.add_subparsers()
    # create subparser for Download Step b, "Build"
    parser1 = subparsers.add_parser("Build", help = "Installs required software")
    parser1._optionals.title = "Arguments"
    parser1.add_argument("-b","--build_folder", help = "Destination folder for downloaded UCSC file(s) (optional; default='squire_build')", type=str, default="squire_build", metavar = "<folder>")
    parser1.add_argument("-s","--software", help = "Install required SQuIRE software and add to PATH - specify 'all' or provide comma-separated list (no spaces) of: STAR,bedtools,samtools,stringtie (optional; default = False)" , type=str, metavar = "<software>", default=False)
    parser1.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default=False)
    parser1.set_defaults(func=s1.main)

    ## create subparser for Download Step 1, "Fetch"
    parser2 = subparsers.add_parser("Fetch", help ="Downloads input files from UCSC")
    parser2._optionals.title = "Arguments"
    parser2.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg37' (required)", type=str, required = True, metavar = "<build>")
    parser2.add_argument("-o","--fetch_folder", help = "Destination folder for downloaded UCSC file(s) (optional; default='squire_fetch')", type=str, default="squire_fetch", metavar = "<folder>")
    parser2.add_argument("-f","--fasta", help = "Download chromosome fasta files for build chromosomes (optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-c","--chrom_info", help = "Download chrom_info.txt file with lengths of each chromosome (optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-r","--rmsk", help = "Download Repeatmasker file (optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-g","--gene", help = "Download UCSC gene annotation(optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-x","--index", help = "Create STAR index, WARNING will take a lot of time and memory (optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
    parser2.add_argument("-k","--keep", help = "Keep downloaded compressed files (optional; default=False)", action = "store_true", default=False)
    parser2.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default=False)
    parser2.set_defaults(func=s2.main)



    ## create subparser for Step1, "Clean"
    parser3 = subparsers.add_parser("Clean", help = "Filters Repeatmasker file for Repeats of interest, collapses overlapping repeats, and returns as BED file.")
    parser3._optionals.title = "Arguments"
    parser3.add_argument("-r","--rmsk", help = "Repeatmasker file (optional; will search 'squire_fetch' folder for rmsk.txt or .out file by default)", type=str, metavar = "<rmsk.txt or file.out>")
    parser3.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg37' (optional; will be basename of rmsk.txt file by default)", type=str, metavar = "<build>")
    parser3.add_argument("-i","--fetch_folder", help = "Destination folder for downloaded UCSC file(s) (optional; default='squire_fetch')", type=str, default="squire_fetch", metavar = "<folder>")
    parser3.add_argument("-o","--clean_folder", help = "Destination folder for output BED file (optional; default = 'squire_clean')", type=str, default = "squire_clean", metavar = "<folder>")
    parser3.add_argument("-c","--repclass", help = "Comma-separated list of desired repeat class/classes, aka superfamily, eg DNA, LTR. Column 12 in repeatmasker file. Can use UNIX wildcard patterns. (optional; default=False)", type=str, metavar = "<classes>")
    parser3.add_argument("-f","--family", help = "Comma-separated list of desired repeat family/families, eg 'ERV1,ERVK,ERVL. Column 13 in repeatmasker file. Can use UNIX wildcard patterns.  (optional; default=False)", type=str, metavar = "<subfamilies>")
    parser3.add_argument("-s","--subfamily", help = "Comma-separated list of desired repeat subfamilies, eg 'L1HS,AluYb'. Column 11 in repeatmasker file. Can use UNIX wildcard patterns.  (optional; default=False)", type=str, metavar = "<families>")
    parser3.add_argument("-e","--extra", help = "Filepath of extra file containing non-reference repeat sequences. Columns should be chr, start, stop, strand, subfamily, and sequence (optional)", type=str, metavar = "<file>", default=False)
    parser3.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default = False)
    parser3.set_defaults(func=s3.main)


    ## create subparser for Step2, 'Map'
    parser4 = subparsers.add_parser('Map', help='Aligns RNAseq reads to STAR index allowing for multiple alignments')
    parser4._optionals.title = "Arguments"
    parser4.add_argument("-1","--read1", help = "RNASeq data fastq file(s); read1 if providing paired end data. If more than one file, separate with commas, no spaces. Can be gzipped.", type = str, metavar = "<file_1.fastq or file_1.fastq.gz>")
    parser4.add_argument("-2","--read2", help = "RNASeq data read2 fastq file(s). if more than one file, separate with commas, no spaces. Can be gzipped.  (optional, can skip or enter 'False' if data is unpaired)", type = str, metavar = "<file_2.fastq or file_2.fastq.gz>")
    parser4.add_argument("-o","--map_folder", help = "Location of SQuIRE Map outputs (optional, default = 'squire_map')", type = str, metavar = "<folder>", default = "squire_map")
    parser4.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch'",type = str, metavar = "<folder>",default="squire_fetch")
    parser4.add_argument("-r","--read_length", help = "Read length (if trim3 selected, after trimming; required).", type = int, metavar = "<int>", required=True)
    parser4.add_argument("-n","--name", help = "Common basename for input files (optional; uses basename of read1 as default)", type = str, metavar = "<str>",default=False)
    parser4.add_argument("-3","--trim3", help = "Trim <int> bases from right end of each read before alignment (optional; default=0).", type = int, default = 0, metavar = "<int>")
    parser4.add_argument("-e","--extra", help = "Filepath of text file containing non-reference repeat sequence and genome information", type=str, metavar = "<file.txt>")
    parser4.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required if more than 1 build in clean_folder)", type=str, metavar = "<build>",default=False)        
    # parser.add_argument("-m","--mask", help = "Separate reads from bamfile that map to plasmid or transgene into another file (optional; default=False)", action = "store_true", default = False)
    parser4.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
    parser4.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

    parser4.set_defaults(func=s4.main)

    ## create subparser for Step3, 'Count'
    parser5 = subparsers.add_parser('Count', help = "Quantifies RNAseq reads aligning to TEs and genes")
    parser5._optionals.title = "Arguments"
    parser5.add_argument("-m","--map_folder", help = "Folder location of outputs from SQuIRE Map (optional, default = 'squire_map')", type = str, metavar = "<folder>",default="squire_map")
    parser5.add_argument("-c","--clean_folder", help = "Folder location of outputs from SQuIRE Clean (optional, default = 'squire_clean')", type = str, metavar = "<folder>",default = "squire_clean")
    parser5.add_argument("-o","--count_folder", help = "Destination folder for output files(optional, default = 'squire_count')", type = str, metavar = "<folder>", default="squire_count")
    parser5.add_argument("-t","--tempfolder", help = "Folder for tempfiles (optional; default=count_folder')", type = str, metavar = "<folder>", default=False)
    parser5.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch)'",type = str, metavar = "<folder>",default="squire_fetch")
    parser5.add_argument("-r","--read_length", help = "Read length (if trim3 selected, after trimming; required).", type = int, metavar = "<int>", required=True)
    parser5.add_argument("-n","--name", help = "Common basename for input files (required if more than one bam file in map_folder)", type = str, metavar = "<str>",default=False)
    parser5.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required if more than 1 build in clean_folder)", type=str, metavar = "<build>",default=False)
    parser5.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
    parser5.add_argument("-s","--strandedness", help = " '0' if unstranded eg Standard Illumina, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard SOLiD (optional,default=0)", type = int, metavar = "<int>", default = 0)
    parser5.add_argument("-e","--EM", help = "Run estimation-maximization on TE counts given number of times (optional, specify 0 if no EM desired; default=auto)", type=str, default = "auto")
    parser5.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

## set which program to be associated with this parser
    parser5.set_defaults(func=s5.main)

    parser6 = subparsers.add_parser("Call",help = """Performs differential expression analysis on TEs and genes""")
    parser6._optionals.title = "Arguments"
    parser6.add_argument("-1","--group1", help = "List of basenames for group1 (Treatment) samples, can also provide string pattern common to all group1 basenames",required = True, type = str, metavar = "<str1,str2> or <*str*>")
    parser6.add_argument("-2","--group2", help = "List of basenames for group2 (Control) samples, can also provide string pattern common to all group2 basenames",required = True, type = str, metavar = "<str1,str2> or <*str*>")
    parser6.add_argument("-A","--condition1", help = "Name of condition for group1",required = True, type = str, metavar = "<str>")
    parser6.add_argument("-B","--condition2", help = "Name of condition for group2",required = True, type = str, metavar = "<str>")
    parser.add_argument("-i","--count_folder", help = "Folder location of outputs from SQuIRE Count (optional, default = 'squire_count')", type = str, metavar = "<folder>",default="squire_count")
    parser6.add_argument("-o","--call_folder", help = "Destination folder for output files (optional; default='squire_call')", type = str, metavar = "<folder>", default="squire_call")
    parser6.add_argument("-s","--subfamily", help = "Compare TE counts by subfamily. Otherwise, compares TEs at locus level (optional; default=False)", action = "store_true", default = False)
    parser6.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
    parser6.add_argument("-N","--projectname", help = "Basename for project", type = str, metavar = "<str>",default=False)
    parser6.add_argument("-f","--output_format", help = "Output figures as html or pdf", type = str, metavar = "<str>",default="html")
    parser6.add_argument("-t","--table_only", help = "Output count table only, don't want to perform differential expression with DESeq2", action = "store_true", default = False)
    parser6.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

    parser6.set_defaults(func=s6.main)

    parser7 = subparsers.add_parser('Draw', help  = """Makes bedgraphs and bedwigs from RNAseq data""")
    parser7._optionals.title = "Arguments"
    parser7.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch')",type = str, metavar = "<folder>",default="squire_fetch")
    parser7.add_argument("-m","--map_folder", help = "Folder location of outputs from SQuIRE Map (optional, default = 'squire_map')", type = str, metavar = "<folder>", default="squire_map")
    parser7.add_argument("-o","--draw_folder", help = "Destination folder for output files (optional; default='squire_draw')", type = str, metavar = "<folder>", default="squire_draw")
    parser7.add_argument("-n","--name", help = "Basename for bam file (required if more than one bam file in map_folder)", type = str, metavar = "<str>",default=False)
    parser7.add_argument("-s","--strandedness", help = " '0' if unstranded, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard  (optional,default=1)", type = int, metavar = "<int>", default = False)
    parser7.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required)", type=str, metavar = "<build>",default=False,required=True)
    parser7.add_argument("-l","--normlib", help = "Normalize bedgraphs by library size (optional; default=False)", action = "store_true", default = False)
    parser7.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
    parser7.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)
    parser7.set_defaults(func=s7.main)


    parser8 = subparsers.add_parser("Seek", help = """Retrieves sequences from chromosome fasta files designated by BED file coordinates""")
    parser8._optionals.title = "Arguments"
    parser8.add_argument("-i","--infile", help = """Repeat genomic coordinates, can be TE_ID, bedfile, or gff (required)""", type=argparse.FileType('r'), metavar = "<file.bed>", required=True)
    parser8.add_argument("-o","--outfile", help = """Repeat sequences output file (FASTA), can use "-" for stdout (required)""", type = argparse.FileType('w'), metavar = "<file.fa>", required=True)
    parser8.add_argument("-g","--genome", help = "Genome build's fasta chromosomes - .fa file or .chromFa folder (required)", type = str, metavar="<file.fa or folder.chromFa>", required=True)
    parser8.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default = False)

    parser8.set_defaults(func=s8.main)


    ## parse the args and call the specific program
    subargs,extra_args = parser.parse_known_args()
    subargs.func(args = subargs)

    # print help usage if no arguments are supplied
    if len(sys.argv)==1 and not ext_args:
        parser.print_help()
        sys.exit(1)

if __name__=="__main__":
    main()

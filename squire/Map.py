#!/usr/bin/env python
# -*- coding: utf-8 -*-
############MODULES#########################
from __future__ import print_function,division
import sys
import os
import errno
import argparse #module that passes command-line arguments into script
from datetime import datetime
import operator #for doing operations on tuple
from operator import itemgetter
import subprocess as sp
from subprocess import Popen, PIPE,STDOUT
import io
import tempfile
#for creating interval from start
from collections import defaultdict #for dictionary
import glob
import re
from six import itervalues
import textwrap
import shutil

####### FUNCTIONS ##############################
def isempty(filepath):
    if os.path.getsize(filepath) == 0:
        raise Exception(filepath + " is empty")

def make_dir(path):
    try:
        original_umask = os.umask(0)
        os.makedirs(path, 0770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    finally:
        os.umask(original_umask)

def get_basename(filepath):
        filename = os.path.basename(filepath)
        filebase = os.path.splitext(filename)[0]
        return filebase

def make_tempfile(basename,step,outfolder):
    tmpfile = tempfile.NamedTemporaryFile(delete=False, dir = outfolder, prefix= basename + "_" + step +  ".tmp")
    tmpname = tmpfile.name
    tmpfile.close()
    return tmpname

def rev_comp(sequence):
    rev_seq = sequence[::-1]
    new_seq=""
    for base in rev_seq:
        if base == "A":
            newbase = "T"
        elif base == "T":
            newbase = "A"
        elif base == "C":
            newbase = "G"
        elif base == "G":
            newbase = "C"
        new_seq += newbase
    return new_seq
def rename_file(oldname,newname):
    shutil.move(oldname, newname)


def combine_files(file1,file2,outfile,debug):
    catcommand_list = ["cat", file1, file2, ">", outfile] #combines multi_aligned reads
    catcommand = " ".join(catcommand_list)
    sp.check_call(["/bin/sh","-c",catcommand])

    if not debug:
        os.unlink(file1)
        os.unlink(file2)

def find_file(folder,pattern,base, wildpos, needed):
    foundfile=False
    if wildpos == 1:
        file_list=glob.glob(folder + "/" + "*" + pattern)
    elif wildpos ==2:
        file_list=glob.glob(folder + "/" + pattern + "*")
    if len(file_list)>1: #if more than one file in folder
        if not base:
            raise Exception("More than 1 " + pattern + " file")
        for i in file_list:
            if base in i:
                foundfile = i        
    elif len(file_list) == 0:
        foundfile = False  
    else:
        foundfile = file_list[0]
    if not foundfile:
        if needed:
            raise Exception("No " + pattern + " file")
        else:
            foundfile = False             
    return foundfile


def align_paired(fastq1,fastq2,pthreads,trim3,index,outfile,gtf,gzip,prefix,read_length,extra_fa):
        ##### ALIGN FASTQ FILE(S) TO GENOME OR REPCHR ########
        #-p16: allows hyperhreading over 16 cores
        #-t: outputs time of alignment
        #--tryhard Puts in maximal effort in finding valid alignments for paired end reads
        #-a:  reports all valid alignments for reads
        #-3 trim3: trims user-specified bases from 3' end of FASTQ sequences (useful for if sequencing read > subsequence length)
        gtf_option = []
        gzip_option = []
        extra_option = []
        if gtf:
            gtf_option = ["--sjdbGTFfile", gtf, "--sjdbOverhang",str(read_length-1), "--twopassMode", "Basic"]
        if gzip:
            gzip_option = ["""--readFilesCommand""", "zcat"]
        if extra_fa:
            extra_option=["""--genomeFastaFiles""",extra_fa]


        add_options = gtf_option + gzip_option + extra_option

        multi_align = ["""--outFilterMultimapNmax""", "100", """--winAnchorMultimapNmax""", "100", "--alignEndsType","EndToEnd" ,"--alignEndsProtrude","100 DiscordantPair"]
        trim = ["""--clip3pNbases""", str(trim3)]
        single_reads = ["""--outFilterScoreMinOverLread""", "0.4", """--outFilterMatchNminOverLread""", """0.4"""]
        #single_reads=[]
        discordant = ["--chimSegmentMin", str(read_length)]
        #discordant = []
        inputs = ["""--genomeDir""", index,"""--readFilesIn""",fastq1, fastq2]
        outputs = [ """--outFileNamePrefix""", prefix, """--outSAMtype""", "BAM Unsorted", "--outSAMattributes", "All","--outSAMstrandField", "intronMotif", "--outSAMattrIHstart", "0"]

        STARcommand_list = ["STAR","""--runThreadN""",str(pthreads)] + trim + multi_align + single_reads + discordant + inputs + outputs + add_options
        STARcommand=" ".join(STARcommand_list)
        sp.check_call(["/bin/sh", "-c", STARcommand])

        STAR_output = prefix + "Aligned.out.bam"

        sortcommand_list = ["samtools", "sort", "-@",str(pthreads), STAR_output, prefix]
        sortcommand = " ".join(sortcommand_list)
        sp.check_call(["/bin/sh", "-c", sortcommand])

        indexcommand_list = ["samtools", "index", outfile]
        indexcommand = " ".join(indexcommand_list)
        sp.check_call(["/bin/sh", "-c", indexcommand])   

        os.unlink(STAR_output)
        os.unlink(prefix + "Log.out")
        os.unlink(prefix + "Log.progress.out")
        rename_file(prefix + "Log.final.out",prefix +".log")

def align_unpaired(fastq,pthreads,trim3,index,outfile,gtf,gzip,prefix,read_length,extra_fa):
        ##### ALIGN FASTQ FILE(S) TO GENOME OR REPCHR ########
        #-p16: allows hyperhreading over 16 cores
        #-t: outputs time of alignment
        #-v3: allows maximum of 3 mismatches to account for population variants, increases stringency of Tag finding
        #-a -m1:  reports all valid alignments for reads with only 1 reportable alignment
        #-3 trim3: trims user-specified bases from 3' end of FASTQ sequences (useful for if sequencing read > subsequence length)

        gtf_option = []
        gzip_option = []
        extra_option = []
        if gtf:
            gtf_option = ["--sjdbGTFfile", gtf, "--sjdbOverhang",str(read_length-1), "--twopassMode", "Basic"]
        if gzip:
            gzip_option = ["""--readFilesCommand""", "zcat"]
        if extra_fa:
            extra_option=["""--genomeFastaFiles""",extra_fa]
        add_options = gtf_option + gzip_option + extra_option

        multi_align = ["""--outFilterMultimapNmax""", "100", """--winAnchorMultimapNmax""", "100"]
        trim = ["""--clip3pNbases""", str(trim3)]
        inputs = ["""--genomeDir""", index,"""--readFilesIn""",fastq]
        outputs = [ """--outFileNamePrefix""", prefix, """--outSAMtype""", "BAM Unsorted", "--outSAMattributes", "All","--outSAMstrandField", "intronMotif", "--outSAMattrIHstart", "0"]

        STARcommand_list = ["STAR","""--runThreadN""",str(pthreads)] + trim + multi_align + inputs + outputs + add_options
        STARcommand=" ".join(STARcommand_list)
        sp.check_call(["/bin/sh", "-c", STARcommand])
        STAR_output = prefix + "Aligned.out.bam"

        sortcommand_list = ["samtools", "sort", "-@",str(pthreads), STAR_output, prefix]
        sortcommand = " ".join(sortcommand_list)
        sp.check_call(["/bin/sh", "-c", sortcommand])


        indexcommand_list = ["samtools", "index", outfile]
        indexcommand = " ".join(indexcommand_list)
        sp.check_call(["/bin/sh", "-c", indexcommand])        

        os.unlink(STAR_output)
        # os.unlink(prefix + "Log.out")
        # os.unlink(prefix + "Log.progress.out")
        rename_file(prefix + "Log.final.out",prefix +".log")

def get_header(bamfile,headerfile):
    samtoolscommand_list = ["samtools","view","-H", bamfile, ">",headerfile]
    samtoolscommand = " ".join(samtoolscommand_list)
    sp.check_call(["/bin/sh", "-c", samtoolscommand])

def mask_reads(infile,extra,chrom_list,basename,outfolder,pthreads,debug):
    read_dict={}
    sam_temp=make_tempfile(basename,"sam_temp",outfolder)
    ectopic_alignments = make_tempfile(basename,"ectopic",outfolder)
    ectopic_reads = infile.replace(".bam","_ectopic.bam")
    nonectopic_reads = infile.replace(".bam","_masked.bam")
    with open(extra,'r') as nonreftable:
        for line in nonreftable:
            line=line.rstrip()
            line=line.split("\t")
            chrom=line[0]
            strand=line[3]
            TE_type=line[5].lower()
            if strand=="Strand":
                continue
            if TE_type=="plasmid":
                chrom_list.append(chrom)
            elif TE_type=="transgene":
                chrom_list.append(chrom)

        get_header(infile,ectopic_reads)
        get_header(infile,nonectopic_reads)

    for chrom in chrom_list:
        search="""'"""+ '$3 ~ /' + chrom + '/' + """'"""
        dupe_command_list = ["samtools","view",infile,chrom, ">", ectopic_alignments] #skips lines if the read has already appeared in the file
        dupe_command = " ".join(dupe_command_list)
        sp.check_call(["/bin/sh", "-c", dupe_command])

        awkcommand_list = ["samtools","view",  infile,  ">",sam_temp] #writes lines in combined_tempfile that are not in unique_tempfile2 -> duplicates
        awkcommand = " ".join(awkcommand_list)
        sp.check_call(["/bin/sh","-c",awkcommand])


        awkcommand_list = ["awk", """'FNR==NR{a[$1]++;next}a[$1]'""", ectopic_alignments, sam_temp,  ">>", ectopic_reads] #writes lines in combined_tempfile that are not in unique_tempfile2 -> duplicates
        awkcommand = " ".join(awkcommand_list)
        sp.check_call(["/bin/sh","-c",awkcommand])

        awkcommand_list = ["awk", """'FNR==NR{a[$1]++;next}!a[$1]'""",  ectopic_alignments, sam_temp,  ">", nonectopic_reads] #writes lines in combined_tempfile that are not in unique_tempfile2 -> duplicates
        awkcommand = " ".join(awkcommand_list)
        sp.check_call(["/bin/sh","-c",awkcommand])
        if not debug:
            os.unlink(ectopic_alignments)



def main(**kwargs):
    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = """Aligns RNAseq reads to STAR index allowing for multiple alignments""")
        parser._optionals.title = "Arguments"
        parser.add_argument("-1","--read1", help = "RNASeq data fastq file; read1 if providing paired end data. If more than one file, separate with commas, no spaces. Can be gzipped. (Required for single-end data; optional for paired-end)", type = str, metavar = "<file_1.fastq or file_1.fastq.gz>")
        parser.add_argument("-2","--read2", help = "RNASeq data read2 fastq file. if more than one file, separate with commas, no spaces. Can be gzipped.  (optional, can skip or enter 'False' if data is unpaired)", type = str, metavar = "<file_2.fastq or file_2.fastq.gz>")
        parser.add_argument("-o","--map_folder", help = "Location of SQuIRE Map outputs (optional, default = 'squire_map')", type = str, metavar = "<folder>", default = "squire_map")
        parser.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch'",type = str, metavar = "<folder>",default="squire_fetch")
        parser.add_argument("-r","--read_length", help = "Read length (if trim3 selected, after trimming; required).", type = int, metavar = "<int>", required=True)
        parser.add_argument("-n","--name", help = "Common basename for input files (optional; uses basename of read1 as default)", type = str, metavar = "<str>",default=False)
        parser.add_argument("-3","--trim3", help = "Trim <int> bases from right end of each read before alignment (optional; default=0).", type = int, default = 0, metavar = "<int>")
        parser.add_argument("-e","--extra", help = "Filepath of text file containing non-reference repeat sequence and genome information", type=str, metavar = "<file.txt>")
        parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required if more than 1 build in clean_folder)", type=str, metavar = "<build>",default=False)
        parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)
        args,extra_args = parser.parse_known_args()

########## I/O #########
    ###### ARGUMENTS ######
    read1=args.read1
    read2=args.read2
    outfolder = args.map_folder
    read_length = args.read_length
    fetch_folder=args.fetch_folder    
    #index = args.index
    basename = args.name
    trim3 = args.trim3
    extra=args.extra
    build=args.build
    #gtf=args.gtf
    # mask=args.mask
    pthreads = args.pthreads
    verbosity=args.verbosity

    ######### START TIMING SCRIPT ############
    if verbosity:
        startTime = datetime.now()
        print("start time is:" + str(startTime) + '\n', file = sys.stderr)# Prints start time
        print(os.path.basename(__file__) + '\n', file = sys.stderr) #prints script name to std err
        print("Script Arguments" + '\n' + "=================", file = sys.stderr)
        args_dict = vars(args)
        for option,arg in args_dict.iteritems():
            print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
        print("\n", file = sys.stderr)

    #### SET DEFAULTS #####

    if not read1 and not read2:
        raise Exception("read1 or read2 must be provided")
    if read2:
        if read2.lower()=="false":
            read2=False
    debug=True 

    ### CHECK INPUTS#####
    index = find_file(fetch_folder,"_STAR",build, 1,True)

    gtf = find_file("squire_fetch","_refGene.gtf",build, 1,True)

    if not basename:
        basename = get_basename(read1)
    make_dir(outfolder)
    outfile = outfolder + "/" + basename + ".bam"

    prefix = outfolder + "/" + basename
    if ".gz" in read1:
        gzip=True
    else:
        gzip = False

    extra_fapath=False

    if extra:
        extra_fapath = outfolder + "/" + get_basename(extra) +  ".fa"
        extra_fa=open(extra_fapath,'wb')
        if verbosity:
            print("Making fasta file from extra file" + "\n", file = sys.stderr)

        previous_chrom=0 #This is needed to avoid reopening chromosome sequence files, which would make the script run time a lot longer.
        buffer_sequence = "N" * 200
        chrom_dict = {}
        seq_dict=defaultdict(str)
        maskchrom_list=[]

        with open(extra,'r') as extra_file:
            nonref_types=["polymorphism","novel","plasmid","transgene"]
            for line in extra_file:
                line = line.rstrip()
                line=line.split("\t")
                chrom=line[0]                    
                start = line[1]
                stop = line[2]
                strand = line[3]
                if strand.lower()=="strand":
                    continue
                taxo = line[4]
                TE_type=line[5].lower()                
                if TE_type not in nonref_types:
                    raise Exception('TE type needs to be "polymorphism","novel","plasmid",or "transgene"')
                chrom = chrom + "_" + TE_type                    
                if not chrom.startswith("chr"):
                    chrom="chr"+chrom 
                #chrom=chrom + "_" + TE_type
                if "plasmid" in TE_type: #if plasmid
                    score = "999"
                    maskchrom_list.append(chrom)
                elif "transgene" in TE_type:
                    score="999"
                    maskchrom_list.append(chrom)
                else: #if insertion polymorphism
                    score = "1000"
                left_flankseq = line[6]
                right_flankseq = line[7]
                TEsequence = line[8]

                sequence=left_flankseq + TEsequence + right_flankseq

                seq_dict[chrom] += sequence + buffer_sequence

        for chrom, sequence in seq_dict.iteritems():
            extra_fa.writelines(">" + chrom + "\n")
            new_chromseq = textwrap.fill(seq_dict[chrom],50)
            extra_fa.writelines(new_chromseq + "\n")

        extra_fa.close()
    else:
        extra_fa=None
    if read1 and not read2: #if single-end
        if read1.endswith(","):
            read1=read1[:-1]
        if verbosity:
            print("Aligning FastQ files " + str(datetime.now()) + "\n",file = sys.stderr)
        align_unpaired(read1,pthreads,trim3,index,outfile,gtf,gzip,prefix, read_length,extra_fapath)

    if read1 and read2:
        if read1.endswith(","):
            read1=read1[:-1]
        if read2.endswith(","):
            read2=read2[:-1]            
        if verbosity:
            print("Aligning FastQ files for Read1 and Read2 " + str(datetime.now()) + "\n",file = sys.stderr)
        align_paired(read1,read2,pthreads,trim3,index,outfile,gtf,gzip,prefix, read_length,extra_fapath)


    # if mask:
    #     mask_reads(outfile,chrom_list,basename,outfolder,pthreads,debug)


    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs" + "\n",file = sys.stderr)

        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
        print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)

###################
if __name__ == "__main__":
    main()

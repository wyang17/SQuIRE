#!/bin/env python

#################### MODULES ###################
from __future__ import print_function,division
import sys
import os
import errno # error code module
import os.path
from datetime import datetime
import argparse #module that passes command-line arguments into script
from pyfaidx import Fasta #Pyfasta module flattens fasta data without spaces or headers so fasta doesn't need to be read into memory
import glob
import tempfile
import subprocess as sp
import re
from subprocess import Popen, PIPE,STDOUT


###################################################
def make_dir(path):
    try:
        original_umask = os.umask(0)
        os.makedirs(path, 0770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    finally:
        os.umask(original_umask)

def isempty(filepath):
    if os.path.getsize(filepath) == 0:
        raise Exception(filepath + " is empty")

def basename(filepath):
        filename = os.path.basename(filepath)
        filebase = os.path.splitext(filename)[0]
        return filebase

class bed(object):
    def __init__(self, line):
        self.chromosome = line[0] # chr = first tab/first in list
        self.start = int(line[1])
        self.end = int(line[2])
        self.name=line[3]
        self.score=float(line[4])
        self.strand = line[5]

class gtf(object):
    def __init__(self,line):
        self.chromosome = line[0] # chr = first tab/first in list
        self.source = (line[1])
        self.feature = (line[2])
        self.start=int(line[3])
        self.end = int(line[4])
        self.score=float(line[5])
        self.strand = str(line[6])
        self.frame = line[7]
        self.attributes = line[8]

def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = "Retrieves sequences from chromosome fasta files")
        parser._optionals.title = "Arguments"
        parser.add_argument("-i","--infile", help = """Repeat genomic coordinates, can be TE_ID, bedfile, or gff (required)""", type=argparse.FileType('r'), metavar = "<file.bed>", required=True)
        parser.add_argument("-o","--outfile", help = """Repeat sequences output file (FASTA), can use "-" for stdout (required)""", type = argparse.FileType('w'), metavar = "<file.fa>", required=True)
        parser.add_argument("-g","--genome", help = "Genome build's fasta chromosomes - .fa file or .chromFa folder (required)", type = str, metavar="<file.fa or folder.chromFa>", required=True)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default = False)

        args,extra_args = parser.parse_known_args()

    ########## PARSE ARGUMENTS  #########
    infile = args.infile
    outfile = args.outfile #if outfile not given, give basename of infile to outfile with .seq extension
    genome = args.genome
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

    ####### CHECK ARGUMENTS AND SET DEFAULTS ###########

    ########### REFERENCES #################
    required_columns = 6 ###For checking if infile is BED format
    previous_chromosome=0 #This is needed to avoid reopening chromosome sequence files, which would make the script run time a lot longer.

    if os.path.isfile(genome): #if genome is file
        chromosome_infile = Fasta(genome)

    ### START FOR LOOP ####

    for line in infile:
        line = line.rstrip() #removes white space at end of line
        if line.startswith("track"):
            continue

        line = line.split("\t")  # returns list of items that were separated by tab in original file

    #########CHECK FILE FORMAT ###########
        column_count = len(line)
        if column_count == 1:
            line=line.split("|")
            bedline=bed(line)
            chromosome = bedline.chromosome # chr = first tab/first in list
            repstart = bedline.start
            repstop = bedline.end
            name=bedline.name
            strand = bedline.strand
            score=bedline.score
            header = str(chromosome) + ":" + str(repstart) + "-" + str(repstop) + "/" + str(strand )+ "/" + str(name)

        elif column_count > 1:
            if re.match("\d+", line[1]):
                bedline = bed(line)
                chromosome = bedline.chromosome # chr = first tab/first in list
                repstart = bedline.start
                repstop = bedline.end
                name=bedline.name
                strand = bedline.strand
                score=bedline.score
                header = str(chromosome) + ":" + str(repstart) + "-" + str(repstop) + "/" + str(strand )+ "/" + str(name)
            else:
                gtfline = gtf(chromosome)
                chromosome = gtfline.chromosome # chr = first tab/first in list
                repstart = gtf.start
                repstop = gtfline.end
                name=gtfline.feature
                strand = gtfline.strand
                score=bedline.score
                header = str(chromosome) + ":" + str(repstart) + "-" + str(repstop) + "/" + str(strand )+ "/" + str(name)




        ######## FETCH SEQUENCES ###############

        if (chromosome != previous_chromosome):  #only reopens new chromosome file if a new chr is reached in coordinates file
            print("Opening " + chromosome + "file" + '\n',file=sys.stderr)
            previous_chromosome = chromosome
            chromstart=0
            if os.path.isdir(genome): #if genome is folder
                chrom_infile = genome + '/'  + chromosome + '.fa'
                chromosome_infile = Fasta(chrom_infile)


        plus_strand_sequence = chromosome_infile[chromosome][repstart:repstop]
        if strand == '-':
            desired_sequence = -(plus_strand_sequence) #if negative strand, give reverse complement in fasta file
        else:
            desired_sequence = plus_strand_sequence

#        fix_BED.writelines(str(chromosome) + "\t" + str(repstart) + "\t" + str(repstop) + "\t" + str(TE_ID) + "\t" + str(score) + "\t" + str(strand) + "\t" + str(repstart) + "\t" + str(repstop) + "\t" + str(RGB) + "\n")

        #FASTA id for each repeat sequence is first 6 columns of the BED file
        outfile.writelines('>' + header + '\n' + str(desired_sequence) + '\n')


    print("Finished writing " + str(outfile) + '\n',file=sys.stderr)


    if verbosity:
        print("Finished writing RepChr FASTA file" + "\n", file = sys.stderr)

    ###### I/O ###############
    infile.close()
    outfile.close()
    ###### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing: " + outfile.name + '\n', file = sys.stderr)

        endTime = datetime.now()

        print('end time is: '+ str(endTime) + '\n', file = sys.stderr)
        print('it took: ' + str(endTime-startTime) + '\n', file = sys.stderr)

###################
if __name__ == "__main__":
    main()

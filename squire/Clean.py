#!/usr/bin/env python

############MODULES #################
from __future__ import print_function,division
import sys
import os
import errno
import argparse #module that passes command-line arguments into script
from datetime import datetime
import subprocess as sp
from pyfaidx import Fasta
import tempfile
from collections import defaultdict
import glob
import re
from fnmatch import fnmatch, fnmatchcase
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

def get_basename(filepath):
        filename = os.path.basename(filepath)
        filebase = os.path.splitext(filename)[0]
        return filebase

def is_list(value):
    if "," in value:
        return value.split(",")
    else:
        return [value]

def uniqify_name(name_list, idfun=None):
   # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in name_list:
        marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    name=",".join(result)
    return name

def dict_col(group,colno,column_list,entries_list,name_list): #add group column and names to be put into dictionary
            column_list.append(colno)
            entries_list.append(is_list(group))
            for string in is_list(group):
                if "*" in string:
                    no_wildcards = string.count("*")
                    regex_tuple = (colno,string)
                    regex_list.append(regex_tuple)
            if "*" in group:
                group = group.replace("*","")
            name_list.append(group)

def split_subF(subF):
    if "," in subF:
        subfamily_list = subF.split(",")
        return subfamily_list
    else:
        return [subF]


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
        if not foundfile:
            if needed:
                raise Exception("No " + pattern + " file")
            else:
                foundfile = False
    elif len(file_list) == 0:
        foundfile = False
    else:
        foundfile = file_list[0]
    return foundfile


def main(**kwargs):


    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = "Filters Repeatmasker file for Repeats of interest, collapses overlapping repeats, returns BED file and count of subfamily copies.")
        parser._optionals.title = "Arguments"
        parser.add_argument("-r","--rmsk", help = "Repeatmasker file (optional; will search 'squire_fetch' folder for rmsk.txt or .out file by default)", type=str, metavar = "<rmsk.txt or file.out>")
        parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (optional; will be basename of rmsk.txt file by default)", type=str, metavar = "<build>")
        parser.add_argument("-o","--clean_folder", help = "Destination folder for output BED file (optional; default = 'squire_clean')", type=str, default = "squire_clean", metavar = "<folder>")
        parser.add_argument("-c","--repclass", help = "Comma-separated list of desired repeat class/classes, aka superfamily, eg DNA, LTR. Column 12 in repeatmasker file. Can use UNIX wildcard patterns. (optional; default=False)", type=str, metavar = "<classes>")
        parser.add_argument("-f","--family", help = "Comma-separated list of desired repeat family/families, eg 'ERV1,ERVK,ERVL. Column 13 in repeatmasker file. Can use UNIX wildcard patterns.  (optional; default=False)", type=str, metavar = "<families>")
        parser.add_argument("-s","--subfamily", help = "Comma-separated list of desired repeat subfamilies, eg 'L1HS,AluYb'. Column 11 in repeatmasker file. Can use UNIX wildcard patterns.  (optional; default=False)", type=str, metavar = "<subfamilies>")
        parser.add_argument("-e","--extra", help = "Filepath of extra file containing non-reference repeat sequences. Columns should be chr, start, stop, strand, subfamily, and sequence (optional)", type=str, metavar = "<file>", default=False)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default = False)

        args,extra_args = parser.parse_known_args()

    ### DEFINE ARGUMENTS ####
    rmsk = args.rmsk
    build = args.build
    outfolder = args.clean_folder
    repclass = args.repclass
    family = args.family
    subfamily = args.subfamily
    verbosity = args.verbosity
    extra = args.extra

    #### DEFINE DEFAULTS #######
    RM_out = False
    print_unique=True
    if not rmsk:  #if infile not given, find *rmsk.txt file in squire_fetch folder and open as infile
        rmsk=find_file(clean_folder,".bed",build,1,True)
    elif ".out" in os.path.basename(rmsk):
        RM_out = True

    if not build:
        rmsk_name=os.path.basename(rmsk)
        if RM_out:
            build=rmsk_name.replace(".out","")
        else:
            build=rmsk_name.replace("_rmsk.txt","")

    make_dir(outfolder)

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

    ##############CREATE DICTIONARY OF COLUMNS AND COLUMN ENTRIES############################
    column_list = []
    entries_list = []
    name_list = [build]
    regex_list = []
    unwantedChr = ["hap", "M", "alt"]
    unwantedClass = ["Simple_repeat", "Satellite", "Unknown", "Low_complexity"]


    get_all=True
    if repclass:
        dict_col(repclass,11,column_list,entries_list,name_list)
        get_all=False
        for i in is_list(repclass):
            if i in unwantedClass:
                unwantedClass.remove(i)
    if family:
        dict_col(family,12,column_list,entries_list,name_list)
        get_all=False
    if subfamily:
        dict_col(subfamily,10,column_list,entries_list,name_list)
        get_all=False
    column_dict = dict(zip(column_list,entries_list))
    if repclass or family or subfamily:
        outfile_name = "_".join(name_list)
    elif get_all:
        outfile_name = build + "_all"

    if extra:
        outfile_name = outfile_name + "_" + get_basename(extra)

    if repclass or family or subfamily or get_all:
            outfilepath = outfolder + "/" + outfile_name + ".bed"
            outfile = open(outfilepath,'w')
    if print_unique:
        unique_filepath = outfolder + "/" + outfile_name + "_copies.txt"
        unique_file = open(unique_filepath,'w')
        unique_file.writelines("Subfamily:Family:Class" + "\t" + "Copies" + "\t" + "Tot.Length" + "\t" "Avg.Length" + "\n") #add header

    #############CREATE BED FILE FROM REPEATMASKER FILE BASED ON SUPPLIED COLUMN INFO#############



    class line_class(object):
        def __init__(self,line,RM_out):

            self.line = line.rstrip() #removes white space at end of line

            if RM_out:
                self.line_split = self.line.split()
                #print(self.line_split)
                self.line_split.insert(0,'0') #insert 0 into beginning of list so column numbers are same between UCSC track and RM output
                #print(self.line_split)
                self.classfamily = str(self.line_split[11])
                self.classfamily_split = self.classfamily.split('/')
                #print(self.classfamily_split)
                if len(self.classfamily_split) == 2:
                    self.line_split[11:12] = self.classfamily_split[0],self.classfamily_split[1]
                elif len(self.classfamily_split) == 1:
                    self.line_split[11:12] = self.classfamily_split[0],self.classfamily_split[0]

                if len(self.line_split)>17:
                    self.line_split=self.line_split[0:17]

            else:
                self.line_split = self.line.split('\t')

            self.listoflists=[]
            for i in self.line_split:
                self.listoflists.append([(i)])

            self.milliDiv = str(self.line_split[2])
            self.chrom = str(self.line_split[5])
            self.start = int(self.line_split[6])
            self.stop = int(self.line_split[7])
            self.subfamily = str(self.line_split[10])
            self.family = str(self.line_split[12])
            self.repclass = str(self.line_split[11])
            self.name=self.subfamily + ":" + self.family + ":" + self.repclass
            self.strand = str(self.line_split[9])
            self.taxo = self.subfamily + ":" + self.family + ":" + self.repclass
            self.length = int(self.stop) - int(self.start)


            if self.strand == "+":
                self.score=self.milliDiv
                self.RGB = "120,91,12" #if TE is on + strand, RGB track will be brown
            else:
                self.strand = "-"
                self.score=self.milliDiv
                self.RGB = "94,137,255" #if TE is on - strand, RGB track will be blue

            self.TE_ID = str(self.chrom) + '|' + str(self.start) + '|' + str(self.stop) + '|' + str(self.name) + '|' + str(self.score) + '|' + str(self.strand)
            self.bed9_list = [self.chrom, str(self.start), str(self.stop), self.TE_ID, self.score, self.strand, str(self.start), str(self.stop), self.RGB]
            self.bed9 = "\t".join(self.bed9_list)
            self.id_list = [self.name, self.strand]
            self.id = ".".join(self.id_list)
            self.name_list = [self.id]

    class collapsed_line(object):
        def __init__(self,line1,line2):
            self.listoflists = []
            #Turn columns in collapsed line into lists
            for i in range(len(line1.line_split)):
                self.listoflists.append(line1.listoflists[i])

                self.listoflists[i].append(line2.line_split[i])

            self.line_split = self.listoflists
            self.start = str(min(int(line1.start), int(line2.start)))
            self.stop = str(max(int(line1.stop), int(line2.stop)))
            self.chrom = line1.chrom
            self.name_list = line1.name_list
            self.name_list.append(line2.id)
            self.name = uniqify_name(self.name_list)
            self.strand = "*"
            self.score = "998"
            self.RGB = "107,114,134"
            self.TE_ID = str(self.chrom) + '|' + str(self.start) + '|' + str(self.stop) + '|' + str(self.name) + '|' + str(self.score) + '|' + str(self.strand)
            self.bed9_list = [self.chrom, str(self.start), str(self.stop), self.TE_ID, self.score, self.strand, str(self.start), str(self.stop), self.RGB]
            self.bed9 = "\t".join(self.bed9_list)

    subfamily_countdict = defaultdict(int)
    subfamily_length=defaultdict(int)

    with open(rmsk) as infile:
        if verbosity:
            if repclass or family or subfamily or get_all:
                print("Adding repeatmasker repeats to BED file" + str(datetime.now()) + "\n", file = sys.stderr)

        for line in infile:
            line_split=line.split()
            if not line_split:
                continue
            elif line_split[0].isdigit()==False:
                continue
            doprint=False
            check_line = line_class(line,RM_out)
            #########FILTER OUT UNWANTED CHROMOSOMES AND CLASSES AND OVERLAPPING TEs #############
            if any(x in check_line.chrom for x in unwantedChr):
                continue

            unwantedPresent=False
            for column in check_line.line_split:
                if any(x in column for x in unwantedClass):
                    unwantedPresent=True

            if unwantedPresent:
                continue

            prev = line_class(line,RM_out)
            subfamily_countdict[prev.taxo] +=1 #adds class, family and subfamily to dictionary, skips overlaps
            subfamily_length[prev.taxo] += prev.length

            #######FILTER FOR DESIRED COLUMN ENTRIES ##############

            for key,listvalue in column_dict.iteritems():
                if any(x in prev.listoflists[int(key)] for x in listvalue):
                    doprint = True

            if regex_list:
                for regex_tuple in regex_list:
                    if not doprint: #check for wildcard if doprint is still not true
                        column_no = regex_tuple[0]
                        pattern = regex_tuple[1]
                        column_entrylist = prev.listoflists[column_no]
                        column_entry = "".join(column_entrylist)
                        if fnmatch(column_entry,pattern):
                            doprint = True
            if get_all:
                doprint = True
            ####### WRITE DESIRED LINES AS BED FILE ##############
            if repclass or family or subfamily or get_all:
                if doprint:
                    outfile.writelines(prev.bed9 + "\n")



    if extra:
        if verbosity:
            print("Adding extra repeats to BED file" + "\n", file = sys.stderr)
        chromdict = defaultdict(int)

        RGB = "107,114,134"

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
                if not chrom.startswith("chr"):
                    chrom="chr"+chrom 
                if "plasmid" in TE_type: #if plasmid
                    score = "999"
                elif "transgene" in TE_type:
                    score="999"
                else: #if insertion polymorphism
                    score = "1000"
                left_flank = line[6]
                right_flank = line[7]
                TEsequence = line[8]
                TE_length = int(stop) - int(start)
                TEseq_length = len(TEsequence)                
                left_flanklength=len(left_flank)
                right_flanklength=len(right_flank)
                seq_length=TEseq_length+left_flanklength + right_flanklength
                newstart = chromdict[chrom] + int(left_flanklength)
                newstop = newstart + TEseq_length
                chromdict[chrom] += seq_length+200
                subfamily_countdict[taxo] +=1
                subfamily_length[taxo] += TE_length
                TE_ID = chrom + "|" + start + "|" + stop + "|" + taxo + "|"  + score + "|" + strand
                bed9_list = [chrom, str(newstart), str(newstop), TE_ID, score, strand, str(start), str(stop), RGB]
                bed9 = "\t".join(bed9_list)
                outfile.writelines(bed9 + "\n")

    if verbosity:
        if repclass or family or subfamily or get_all:
            print("Finished writing bedfile " + outfile.name + " " + str(datetime.now())  + "\n",file = sys.stderr)

    #### Print unique file
    if print_unique:
        if verbosity:
            print("Writing subfamily counts and length "+ str(datetime.now())  + "\n",file = sys.stderr)
        ### FIND REPEAT SUBFAMILIES ####
        for taxo, count in sorted(subfamily_countdict.items()):
            taxo_list = taxo.split("|")
            taxo_tab = taxo.replace("|","\t")
            taxo_length = subfamily_length[taxo]
            taxo_avg = int(taxo_length)/int(count)
            unique_file.writelines(taxo_tab + "\t" + str(count) + "\t" + str(taxo_length) + "\t" + "{0:.2f}".format(taxo_avg) + "\n")

        if verbosity:
            print("Finished writing subfamily counts and length " + unique_file.name + " " + str(datetime.now())  + "\n",file = sys.stderr)


    ######### I/O ########
    if repclass or family or subfamily or get_all:
        outfile.close()
    if print_unique:
        unique_file.close()

    ###### STOP TIMING SCRIPT #######################
    if verbosity:

        endTime = datetime.now()

        print('end time is: '+ str(endTime) + '\n', file = sys.stderr)
        print('it took: ' + str(endTime-startTime) + '\n', file = sys.stderr)

#########################################################
if __name__ == "__main__":
    main()


########################################

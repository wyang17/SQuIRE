#!/usr/bin/env python

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
import shutil


def find_file(folder,pattern,base, wildpos):
    foundfile=False
    needed=False
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


def make_tempfile(basename, step, outfolder):
    tmpfile = tempfile.NamedTemporaryFile(delete=False, dir = outfolder, prefix= basename + "_" + step +  ".tmp")
    tmpname = tmpfile.name
    tmpfile.close()
    return tmpname

def filter_files(file_in,file_out, string, column):
    command = "'$" + str(column) + "==" + '"' + string + '"'+ "'"
    pastecommandlist = ["awk", "-v", "OFS='\\t'",command,file_in, ">", file_out]
    pastecommand = " ".join(pastecommandlist)
    sp.check_call(["/bin/sh","-c",pastecommand])

def rename_file(oldname,newname):
    shutil.move(oldname, newname)

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


def sort_coord(infile, outfile,chrcol,startcol):
    chrfieldsort = "-k" + str(chrcol) + "," + str(chrcol)
    startfieldsort = "-k" + str(startcol) + "," + str(startcol) + "n"
    sort_command_list = ["sort",chrfieldsort,startfieldsort, infile, ">", outfile]
    sort_command = " ".join(sort_command_list)
    sp.check_call(["/bin/sh", "-c", sort_command])
    os.unlink(infile)

def bedgraph(infile,strandedness,outfolder,basename,normlib,pthreads,bedgraph_list):
    if strandedness==1:
        stranded_yesno= "Stranded"
        plus_bedgraph_unique=outfolder + "/" + basename + "Signal.Unique.str2.out.bg"
        minus_bedgraph_unique = outfolder + "/" + basename + "Signal.Unique.str1.out.bg"
        plus_bedgraph_multi=outfolder + "/" + basename + "Signal.UniqueMultiple.str2.out.bg"
        minus_bedgraph_multi = outfolder + "/" + basename + "Signal.UniqueMultiple.str1.out.bg"
    elif strandedness==2:
        stranded_yesno= "Stranded"
        plus_bedgraph_multi=outfolder + "/" + basename + "Signal.UniqueMultiple.str1.out.bg"
        minus_bedgraph_multi = outfolder + "/" + basename + "Signal.UniqueMultiple.str2.out.bg"
        plus_bedgraph_unique=outfolder + "/" + basename + "Signal.Unique.str1.out.bg"
        minus_bedgraph_unique = outfolder + "/" + basename + "Signal.Unique.str2.out.bg"
    else:
        stranded_yesno="Unstranded"
        bedgraph_unique = outfolder + "/" + basename + "Signal.Unique.str1.out.bg"
        bedgraph_multi = outfolder + "/" + basename + "Signal.UniqueMultiple.str1.out.bg"

    inputs = ["""--inputBAMfile""", infile]
    outputs = ["""--outWigType""", "bedGraph", """--outWigStrand""", stranded_yesno, """--outFileNamePrefix""", outfolder + "/" + basename]

    if not normlib:
        normalization=["""--outWigNorm""", "None"]
    else:
        normalization=["""--outWigNorm""", "RPM"]
    STARcommand_list = ["STAR","""--runMode""","inputAlignmentsFromBAM","""--runThreadN""",str(pthreads)] + inputs + outputs + normalization
    STARcommand=" ".join(STARcommand_list)
    sp.check_call(["/bin/sh", "-c", STARcommand])

    if strandedness !=0:
        sort_coord(plus_bedgraph_unique,outfolder + "/" + basename + "_plus_unique.bedgraph",1,2)
        sort_coord(minus_bedgraph_unique,outfolder + "/" + basename + "_minus_unique.bedgraph",1,2)
        sort_coord(plus_bedgraph_multi,outfolder + "/" + basename + "_plus_multi.bedgraph",1,2)
        sort_coord(minus_bedgraph_multi,outfolder + "/" + basename + "_minus_multi.bedgraph",1,2)

        bedgraph_list += [outfolder + "/" + basename + "_plus_unique.bedgraph",outfolder + "/" + basename + "_minus_unique.bedgraph",outfolder + "/" + basename + "_plus_multi.bedgraph",outfolder + "/" + basename + "_minus_multi.bedgraph"]
    else:
        sort_coord(bedgraph_unique,outfolder + "/" + basename + "_unique.bedgraph",1,2)
        sort_coord(bedgraph_multi,outfolder + "/" + basename + "_multi.bedgraph",1,2)
        bedgraph_list += [outfolder + "/" + basename + "_unique.bedgraph",outfolder + "/" + basename + "_multi.bedgraph"]

def make_bigwig(chrominfo,bedgraph_list):
    for bedgraph in bedgraph_list:
        outfile=bedgraph + ".bw"
        igvcommand_list = ["bedGraphToBigWig",bedgraph, chrominfo,outfile]
        igvcommand=" ".join(igvcommand_list)
        sp.check_call(["/bin/sh", "-c", igvcommand])

def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = """Makes unique and multi bedgraph files""")
        parser._optionals.title = "Arguments"
        parser.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch'",type = str, metavar = "<folder>",default="squire_fetch")
        parser.add_argument("-m","--map_folder", help = "Folder location of outputs from SQuIRE Map (optional, default = 'squire_map')", type = str, metavar = "<folder>", default="squire_map")
        parser.add_argument("-o","--draw_folder", help = "Destination folder for output files (optional; default='squire_draw')", type = str, metavar = "<folder>", default="squire_draw")
        parser.add_argument("-n","--name", help = "Basename for bam file (required if more than one bam file in map_folder)", type = str, metavar = "<str>",default=False)
        parser.add_argument("-s","--strandedness", help = " '0' if unstranded, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard  (optional,default=1)", type = int, metavar = "<int>", default = False)
        parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required)", type=str, metavar = "<build>",default=False,required=True)
        parser.add_argument("-l","--normlib", help = "Normalize bedgraphs by library size (optional; default=False)", action = "store_true", default = False)
        parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

        args,extra_args = parser.parse_known_args()
########## I/O #########
    ###### ARGUMENTS ######
    fetch_folder=args.fetch_folder
    map_folder = args.map_folder
    outfolder=args.draw_folder
    basename = args.name
    verbosity=args.verbosity
    build=args.build
    pthreads = args.pthreads
    strandedness=args.strandedness
    normlib=args.normlib
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
    make_dir(outfolder)
    infile = find_file(map_folder,".bam",basename, 1)
    if not basename:
        basename = get_basename(infile)


    if verbosity:
        print("Making unique and total bedgraphs "+ str(datetime.now())  + "\n",file = sys.stderr)
    chrominfo = find_file(fetch_folder,"_chromInfo.txt",build,1)
    bedgraph_list=[]
    bedgraph(infile,strandedness,outfolder,basename,normlib,pthreads,bedgraph_list)
    if verbosity:
        print("Making unique and total bigwigs "+ str(datetime.now())  + "\n",file = sys.stderr)
    make_bigwig(chrominfo,bedgraph_list)
    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
        print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)

###################
if __name__ == "__main__":
    main()

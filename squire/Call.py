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
from collections import defaultdict #for dictionary
import glob
import re
from six import itervalues
import call_deseq2
import call_deseq2_prefilter
import shutil


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

def get_groupfiles(group,gene_files,subF_files,TE_files,subfamily):
    if "*" not in group:
        if "," in group:
            group_list = group.split(",")    
        else:
            group_list=[group]    
        for sample in group_list: 
            if subfamily:
                subF_files.append(find_file(count_folder,"_subFcounts.txt",sample,1,True))
            else:
                TE_files.append(find_file(count_folder,"_TEcounts.txt",sample,1,True))
            gene_files.append(glob.glob(count_folder + "/" + sample + "_refGenecounts.txt")[0])
    elif "*" in group:
        if subfamily:
            subF_files+=(glob.glob(count_folder + "/" + group + "_subFcounts.txt"))
        else:
            TE_files+=(glob.glob(count_folder + "/" + group + "_TEcounts.txt"))
        gene_files+=(glob.glob(count_folder + "/" + group + "_refGenecounts.txt"))
        group_list=[get_basename(gene_file).replace("_refGenecounts","") for gene_file in (glob.glob(count_folder + "/" + group + "_refGenecounts.txt"))]
    return group_list  

def create_count_dict(infilepath,count_dict,stringtie_list):
    name=get_basename(infilepath).replace("_refGenecounts","")
    with open(infilepath,'r') as infile:
        header = infile.readline().rstrip()
        for line in infile:
            line = line.rstrip()
            line = line.split("\t")
            chrom = line[0]
            start=line[1]
            stop = line[2]
            gene_ID = line[3]
            fpkm=line[4]
            strand = line[5]
            count = line[6]        
            if (gene_ID,strand) in stringtie_list:
                continue
            if (gene_ID,strand) not in count_dict:
                count_dict[(gene_ID,strand)] = {name:count}
            else:
                count_dict[(gene_ID,strand)][name]=count

#subF_file_header.writelines("Sample" + "\t" + "aligned_libsize" + "\t" + "Subfamily:Family:Class" + "\t" + "copies"  + "\t" + "EM_iteration" + "\t" + "uniq_counts" + "\t" + "tot_counts_preEM" + "\t" + "tot_counts_postEM" + "\t" + "tot_reads" + "\t" + "avg_conf"  + "\n")

def create_TE_dict(infilepath,sample_count_dict,threshold):
    conf_dict={}
    count_dict={}
    with open(infilepath,'r') as infile:
        header = infile.readline().rstrip()
        for line in infile:            
            line = line.rstrip()
            line = line.split("\t")
            if "milliDiv" in line[12]:
                continue
            TE_ID = line[3]
            strand = line[5]
            milliDiv = int(line[12])
            count = str(int(float(line[15])))
            conf = float(line[17])
            sample = line[6]                       
            if (TE_ID,strand) not in sample_count_dict:
                sample_count_dict[(TE_ID,strand)] = {sample:count}
                conf_dict[(TE_ID,strand)] = [conf]
                count_dict[(TE_ID,strand)] = [count]
            else:
                sample_count_dict[(TE_ID,strand)][sample]=count
                conf_dict[(TE_ID,strand)].append(conf)
                count_dict[(TE_ID,strand)].append(count)
    
    for TE_tuple,conf_list in conf_dict.iteritems():
        mean_conf=sum(conf_list)/len(conf_list)
        if mean_conf <= threshold:
            sample_count_dict.pop(TE_tuple, None)
    for TE_tuple,TEcount_list in count_dict.iteritems():
        TEcount_list=TEcount_list = [int(i) for i in TEcount_list]
        mean_count=sum(TEcount_list)/len(TEcount_list)
        if mean_count <= 5:
            sample_count_dict.pop(TE_tuple, None)        

def create_subfamily_dict(infilepath,count_dict):
    TE_classes=["LTR","LINE","SINE","Retroposon","DNA","RC"] 
    with open(infilepath,'r') as infile:
        for line in infile:
            line = line.rstrip()
            line = line.split("\t")
            taxo = line[2]
            count=line[6]                  
            if any(x in taxo for x in TE_classes):                
                if count=="tot_counts":
                    continue
                else:
                    count = str(int(round(float(line[5]))))
                sample = line[0]
                if taxo not in count_dict:
                    count_dict[taxo] = {sample:count}
                else:   
                    count_dict[taxo][sample]=count

def combinefiles(infile,catfile):
    with open(catfile, 'a') as outFile:
        with open(infile, 'rb') as inFile:
            shutil.copyfileobj(inFile, outFile)                     

def create_rscript(count_table,coldata,outfolder,output_format,projectname,verbosity,pthreads,prefilter,condition1,condition2,label_no):
    r_script = make_tempfile(projectname,"R_script",outfolder)
    outfolder=os.path.abspath(outfolder)
    count_table = os.path.abspath(count_table)
    coldata=os.path.abspath(coldata)
    if prefilter:
        call_deseq2_prefilter.write_Rscript(r_script)
    else:
        call_deseq2.write_Rscript(r_script)
    #outfile = open(outfolder + "/" + projectname + "call_results.txt","w")

    if verbosity:
        print("Creating DESeq2 results"+ str(datetime.now())  + "\n",file = sys.stderr)


    Rcommandlist = ["Rscript", r_script, count_table,coldata,outfolder,projectname,pthreads,condition1,condition2,str(label_no)]
    Rcommand = " ".join(Rcommandlist)
    sp.check_call(["/bin/sh","-c",Rcommand])

    # if output_format=="html":
    #     render_command="rmarkdown::render('" + r_script + "')"
    # elif output_format == "pdf":
    #     render_command="rmarkdown::render('" + r_script + "', 'pdf_document')"
    # Rcommandlist = ["R","-e", render_command]
    # Rcommand = " ".join(Rcommandlist)
    # sp.check_call(["/bin/sh","-c",Rcommand])

    os.unlink(r_script)

def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = """Performs differential expression analysis on TEs and genes""")
        parser._optionals.title = "Arguments"
        parser.add_argument("-1","--group1", help = "List of basenames for group1 (Treatment) samples, can also provide string pattern common to all group1 basenames with * ",required = True, type = str, metavar = "<str1,str2> or <*str*>")
        parser.add_argument("-2","--group2", help = "List of basenames for group2 (Control) samples, can also provide string pattern common to all group2 basenames with * ",required = True, type = str, metavar = "<str1,str2> or <*str*>")
        parser.add_argument("-A","--condition1", help = "Name of condition for group1",required = True, type = str, metavar = "<str>")
        parser.add_argument("-B","--condition2", help = "Name of condition for group2",required = True, type = str, metavar = "<str>")
        parser.add_argument("-i","--count_folder", help = "Folder location of outputs from SQuIRE Count (optional, default = 'squire_count')", type = str, metavar = "<folder>",default="squire_count")
        parser.add_argument("-o","--call_folder", help = "Destination folder for output files (optional; default='squire_call')", type = str, metavar = "<folder>", default="squire_call")
        parser.add_argument("-s","--subfamily", help = "Compare TE counts by subfamily. Otherwise, compares TEs at locus level (optional; default=False)", action = "store_true", default = False)
        parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
        parser.add_argument("-N","--projectname", help = "Basename for project, default='SQuIRE'",type = str, metavar = "<str>",default="SQuIRE")
        parser.add_argument("-f","--output_format", help = "Output figures as html or pdf", type = str, metavar = "<str>",default="html")
        parser.add_argument("-t","--table_only", help = "Output count table only, don't want to perform differential expression with DESeq2", action = "store_true", default = False)
        #parser.add_argument("-c","--cluster", help = "Want to cluster samples by gene and TE expression", action = "store_true", default = False)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

        args,extra_args = parser.parse_known_args()
########## I/O #########
    ###### ARGUMENTS ######

    group1 = args.group1
    group2 = args.group2
    condition1=args.condition1
    condition2 = args.condition2
    count_folder = args.count_folder
    outfolder=args.call_folder
    verbosity=args.verbosity
    projectname = args.projectname
    subfamily=args.subfamily
    output_format = args.output_format
    pthreads= args.pthreads
    table_only=args.table_only
    debug = True
    label_no=20
    threshold=0
    ######### Call TIMING SCRIPT ############
    if verbosity:
        CallTime = datetime.now()
        print("Script start time is:" + str(CallTime) + '\n', file = sys.stderr)# Prints Call time
        print("Script Arguments" + '\n' + "=================", file = sys.stderr)
        args_dict = vars(args)
        for option,arg in args_dict.iteritems():
            print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
        print("\n", file = sys.stderr)
    if os.path.isfile(outfolder):
        raise Exception (outfolder + " exists as a file" )


    make_dir(outfolder)
    gene_files = []  
    subF_files=[]   
    TE_files=[]    


    group1_list=get_groupfiles(group1,gene_files,subF_files,TE_files,subfamily)
    group2_list=get_groupfiles(group2,gene_files,subF_files,TE_files,subfamily)

    count_dict = {}
    gene_list=set()

    TE_dict={}

    if subfamily: 
        for subF in subF_files:
            combinefiles(subF,subF_combo)       
        create_subfamily_dict(subF_combo,TE_dict)
    else:
        for TE in TE_files:
            combinefiles(TE,TE_combo)        
        create_TE_dict(TE_combo,TE_dict,threshold)

    for genefile in gene_files:        
        create_count_dict(genefile,count_dict,gene_list)

    coldata=outfolder + "/" + projectname + "_coldata.txt"
    with open(coldata,'w') as datafile:
        datafile.writelines("sample" + "\t" + "condition" + "\n")
        for group1_sample in group1_list:
            datafile.writelines(group1_sample + "\t" + condition1 + "\n")
        for group2_sample in group2_list:
            datafile.writelines(group2_sample + "\t" + condition2 + "\n")
    if subfamily:
        counttable = outfolder + "/" + projectname + "_gene_subF_counttable.txt"
    else:
        counttable = outfolder + "/" + projectname + "_gene_TE_counttable.txt"        

    with open(counttable,'w') as DEfile:
        sample_list = group1_list + group2_list
        header_list = ["gene_id"] + sample_list
        header = "\t".join(header_list)
        DEfile.writelines(header + "\n")
        for gene_key,sample_dict in count_dict.iteritems():
            if type(gene_key) is tuple:
                gene=",".join(gene_key)
            else:
                gene=gene_key
            count_list = []
            for sample in sample_list:
                if sample in sample_dict:
                    count_list.append(str(sample_dict[sample]))
                else:
                    count_list.append("0")
            countline = "\t".join(count_list)
            DEfile.writelines(gene + "\t" + countline + "\n")
        for TE_key,sample_dict in TE_dict.iteritems():
            TE_out=",".join(TE_key)
            count_list = []
            for sample in sample_list:
                if sample in sample_dict:
                    count_list.append(str(sample_dict[sample]))
                else:
                    count_list.append("0")
            countline = "\t".join(count_list)
            DEfile.writelines(TE_out + "\t" + countline + "\n")

    prefilter = True
    if not table_only:
        create_rscript(counttable,coldata,outfolder,output_format,projectname,verbosity,str(pthreads),prefilter,condition1,condition2,label_no)

    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
        print('it took: ' + str(endTime-CallTime) + "\n", file = sys.stderr)

###################
if __name__ == "__main__":
    main()

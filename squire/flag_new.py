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
#for creating interval from StringTiet
from collections import defaultdict #for dictionary
import glob
import re
from six import itervalues
import shutil
from math import ceil

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


def isempty(filepath):
    if os.path.getsize(filepath) == 0:
        raise Exception(filepath + " is empty")


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


def isempty(filepath):
    if os.path.getsize(filepath) == 0:
        raise Exception(filepath + " is empty")


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


def Stringtie(bamfile,outfolder,basename,strandedness,read_length,pthreads,gtf, verbosity,outgtf):
    ###Stringtie parameters
    extra_files=True
    if strandedness ==1:
        stringtie_strand = "--rf"
    elif strandedness == 2:
        stringtie_strand = "--fr"
    else:
        stringtie_strand = ""
    if gtf:
        inputs = ["-G", gtf, bamfile]
        pct_max_fpkm=0.1
        flanklength = 10
        flankdepth = 1
        read_gap = 50
        min_tx_length=200
        max_multi_pct = 1
        min_coverage = 2.5
        TEoptions = [stringtie_strand, "-f",str(pct_max_fpkm),"-m", str(min_tx_length), "-a", str(flanklength), "-j", str(flankdepth), "-g", str(read_gap), "-M", str(max_multi_pct), "-c", str(min_coverage),"-e"]
    else:
        inputs = [bamfile]
        pct_max_fpkm=0.1
        flanklength = 10
        flankdepth = .1
        read_gap = 50
        min_tx_length=200
        max_multi_pct = 1.0
        min_coverage = 2.5
        TEoptions = [stringtie_strand,"-l",basename, "-f",str(pct_max_fpkm),"-m", str(min_tx_length), "-a", str(flanklength), "-j", str(flankdepth), "-g", str(read_gap), "-M", str(max_multi_pct), "-c", str(min_coverage), "-t"]
    runoptions = ["-p", str(pthreads), ]
    if verbosity:
        if gtf:
            print("Running Guided Stringtie on each bamfile " + basename + " " + str(datetime.now()) + "\n",file = sys.stderr)
        else:
            print("Running Unguided Stringtie on each bamfile " + basename + " " + str(datetime.now()) + "\n",file = sys.stderr)            
    outputs=["-o", outgtf]
    if extra_files:
        out_abund = outgtf.replace("outgtf","outabund")
        outputs= outputs + ["-A", out_abund]
    StringTiecommand_list = ["stringtie"] + runoptions + TEoptions + outputs + inputs
    StringTiecommand=" ".join(StringTiecommand_list)
    sp.check_call(["/bin/sh", "-c", StringTiecommand])


def sort_coord(infile, outfile,chrcol,startcol):
    chrfieldsort = "-k" + str(chrcol) + "," + str(chrcol)
    startfieldsort = "-k" + str(startcol) + "," + str(startcol) + "n"
    sort_command_list = ["sort",chrfieldsort,startfieldsort, infile, ">", outfile]
    sort_command = " ".join(sort_command_list)
    sp.check_call(["/bin/sh", "-c", sort_command])

def sort_counts(tempfile,headerfile,countsfile, field,debug):
    sorted_countsfile = tempfile + ".sorted"
    field_command = str(field) + "," + str(field) + "rn"
    sort_command_list = ["sort","-k",field_command, tempfile, ">", sorted_countsfile]
    sort_command = " ".join(sort_command_list)
    sp.check_call(["/bin/sh", "-c", sort_command])
    #add header to countsfile
    catcommand_list = ["cat", headerfile, sorted_countsfile, ">",countsfile ] #combines multi_aligned reads
    catcommand = " ".join(catcommand_list)
    sp.check_call(["/bin/sh","-c",catcommand])
    if not debug:
        os.unlink(sorted_countsfile)
        os.unlink(tempfile)
        os.unlink(headerfile)


def combine_files(file1,file2,outfile,debug):
    catcommand_list = ["cat", file1, file2, ">", outfile] #combines multi_aligned reads
    catcommand = " ".join(catcommand_list)
    sp.check_call(["/bin/sh","-c",catcommand])
    if not debug:
        os.unlink(file1)
        os.unlink(file2)


def combinefiles(infile,catfile):
    with open(catfile, 'a') as outFile:
        with open(infile, 'rb') as inFile:
            shutil.copyfileobj(inFile, outFile) 


def filter_TE(infile,TE_strandfile,tx_strandfile,strandedness):
    TE_classes=["LTR","LINE","SINE","Retroposon","DNA","RC"]
    TEout = open(TE_strandfile,'w')
    txout=open(tx_strandfile,'w')
    with open(infile,'r') as filterin:
        header=filterin.readline()
        for line in filterin:
            line = line.rstrip()
            line = line.split("\t")
            TE_class=line[11].split(":")[2]
            if TE_class in TE_classes:
                TE_ID = line[3]
                chrom = line[0]
                strand  = line[5]
                start = line[1]
                stop = line[2]
                fpkm = float(line[4])
                sample=line[6]
                libsize=line[7]
                TE_info=line[8:14]
                subfamily=line[11]            
                milliDiv=line[12]
                unique_count=float(line[14])
                total_count = float(line[15])
                total_reads=int(line[16])
                score = float(line[17])
                if strand==".":
                    outstrand = line[13]
                else:
                    outstrand = strand
                if total_count >= 1:
                    outline = [chrom, start,stop, TE_ID,fpkm,outstrand,unique_count,total_count,total_reads,score]
                    stringout = [str(i) for i in outline ]
                    txout.writelines("\t".join(stringout) + "\n")
                    outline = [TE_info[0], TE_info[1],TE_info[2], TE_ID,fpkm,strand,unique_count,total_count,total_reads,score]
                    stringout = [str(i) for i in outline ]
                    TEout.writelines("\t".join(stringout) + "\n")                    
        TEout.close()
        txout.close()


class gene_info(object):
    def __init__(self,line):
        self.Gene_ID = line[0]
        self.Gene_name = line[1]
        self.chrom = line[2]
        self.strand  = line[3]
        self.start = str(int(line[4])-1) #changes from 1 base to 0-base
        self.stop = int(line[5])
        self.coverage = float(line[6])
        self.fpkm = float(line[7])
        self.tpm = float(line[8])
        self.counts=0
        self.tx_IDs=set()
        self.tx_ID_string=",".join(self.tx_IDs)
        self.flagout=[self.Gene_ID,self.fpkm,self.counts]
        self.countsout=[self.chrom,self.start,self.stop,self.Gene_ID,int(self.counts),self.strand,self.tx_ID_string]
    def add_counts(self,counts):
        self.counts += counts        
    def add_tx(self,txID):
        self.tx_IDs.add(txID)
        self.tx_ID_string=",".join(self.tx_IDs)
        self.flagout = [self.Gene_ID,self.fpkm,self.counts]
        self.countsout = [self.chrom,self.start,self.stop,self.Gene_ID,int(round(self.counts)),self.strand,self.tx_ID_string]


def filter_abund(infile,gene_dict,notinref_dict):
    with open(infile,'r') as filterin:        
        for line in filterin:
            line = line.rstrip()
            line = line.split("\t")
            if "Gene" in line[0] and "TPM" in line[-1]:
                continue
            gene_data=gene_info(line)          
            if not notinref_dict:
                gene_dict[(gene_data.Gene_ID,gene_data.strand)] = gene_data
            else:
                if gene_data.Gene_ID in notinref_dict:
                    gene_dict[(gene_data.Gene_ID,gene_data.strand)] = gene_data

            
class gtfline(object):
    def __init__(self,line):
        self.line=line
        self.chrom = line[0]
        self.source=line[1]
        self.category=line[2]
        self.start = (int(line[3])-1)
        self.stop = int(line[4])
        self.score=(line[5])
        self.strand  = line[6]
        self.frame=line[7]
        self.attributes=line[8].split("; ")
        for attribute_pair in self.attributes:
            self.attribute = attribute_pair.replace(" ","").split('"')
            if self.attribute[0]=="FPKM":
                self.fpkm=float(self.attribute[1])
            elif self.attribute[0]=="TPM":
                self.tpm=float(self.attribute[1])
            elif self.attribute[0]=="gene_id":
                self.Gene_ID=self.attribute[1]
            elif self.attribute[0]=="cov" :
                self.coverage=float(self.attribute[1])
            elif self.attribute[0]== "transcript_id":
                self.transcript_id=self.attribute[1]
    def replace_geneid(self,newgeneid):
        newgeneid=[str(x) for x in newgeneid]        
        newgeneid=",".join(newgeneid)
        self.attributes[0] = "gene_id" + " " + '"' + newgeneid + '"'
        attributesout = self.attributes= "; ".join(self.attributes)
        gtfout = [self.chrom,self.source,self.category,self.start+1,self.stop,self.score,self.strand,self.frame,attributesout]
        self.gtfout = [str(i) for i in gtfout]


def filter_tx(infile,gene_bed,genefile,tx_gene_dict,gene_dict,read_length,genebeddict):    
    if genefile:
        geneout = open(genefile,'w')
        with open(gene_bed,'r') as genein:
            for line in genein:
                line=line.rstrip()
                line=line.split("\t")            
                chrom=line[0]
                start=line[1]
                stop=line[2]
                ref_transcript=line[3]
                strand = line[5]
                genebeddict[ref_transcript]=line 
    with open(infile,'r') as filterin:   
        for line in filterin:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            line = line.split("\t")
            gtf_line = gtfline(line[0:9])        
            if len(line) == 9:
                if gtf_line.category=="exon":  
                    transcribed_length=int(gtf_line.stop) - int(gtf_line.start)                        
                    counts = gtf_line.coverage*transcribed_length/int(read_length)
                    if counts > 0:
                        gene_dict[(gtf_line.Gene_ID,gtf_line.strand)].add_counts(counts)  
                        gene_dict[(gtf_line.Gene_ID,gtf_line.strand)].add_tx(gtf_line.transcript_id)    
                        tx_gene_dict[gtf_line.transcript_id].add((gtf_line.Gene_ID,gtf_line.strand))
                if gtf_line.category=="transcript" and  len(line) == 9 and genefile:
                    if float(gtf_line.fpkm) > 0 or float(gtf_line.tpm) > 0:
                        bedline=genebeddict[gtf_line.transcript_id]
                        geneout.writelines("\t".join(bedline) + "\n")
            else:
                ref_line=gtfline(line[9:18])
                tx_gene_dict[gtf_line.transcript_id].add((ref_line.Gene_ID,ref_line.strand))
                gene_dict[(ref_line.Gene_ID,ref_line.strand)].add_tx(gtf_line.transcript_id)  
    if genefile:
        geneout.close()   



def gtf_to_bed(gtf,bed,debug):    
    #convert gtf to genepred
    genepred=gtf.replace("gtf","genepred")
    gtftogenepredcommand_list = ["gtfToGenePred",gtf,genepred] 
    gtftogenepredcommand=" ".join(gtftogenepredcommand_list)
    sp.check_call(["/bin/sh", "-c", gtftogenepredcommand]) 
    #convert genepred to bed
    genepredtobedcommand_list = ["genePredToBed ",genepred,bed] 
    genepredtobedcommand=" ".join(genepredtobedcommand_list)
    sp.check_call(["/bin/sh", "-c", genepredtobedcommand])     
    if not debug:
        os.unlink(genepred)  


def closest_sense(TEfile,genefile,outfile,strandedness):
    ######## closest WITH BED FILE #########################
    #identify nearest gene that is upstream of TE with respect to reference gnome
    if strandedness == 0:
        closest_list = ["bedtools", "closest", "-a",TEfile,"-b",genefile,"-D","b","-t","all","-s",">",outfile]
        closest_command = " ".join(closest_list)
        sp.check_call(["/bin/sh", "-c", closest_command])
    else:
        closest_list = ["bedtools", "closest", "-a",TEfile,"-b",genefile,"-D","b","-t","all","-s",">",outfile]
        closest_command = " ".join(closest_list)
        sp.check_call(["/bin/sh", "-c", closest_command])


def closest_antisense(TEfile,genefile,outfile,strandedness):
    ######## closest WITH BED FILE #########################
    #identify nearest gene that is upstream of TE with respect to reference gnome
    if strandedness == 0:
        closest_list = ["bedtools", "closest", "-a",TEfile,"-b",genefile,"-D","b","-t","all","-S",">",outfile]
        closest_command = " ".join(closest_list)
        sp.check_call(["/bin/sh", "-c", closest_command])
    else:
        closest_list = ["bedtools", "closest", "-a",TEfile,"-b",genefile,"-D","b","-t","all","-S",">",outfile]
        closest_command = " ".join(closest_list)
        sp.check_call(["/bin/sh", "-c", closest_command])    


def intersect_wawb(TEfile,transcriptfile,outfile,strandedness):
    ######## intersect WITH BED FILE #########################
    #keep read if 50% of read overlaps with TE range
    if strandedness == 0:
        intersect_list = ["bedtools", "intersect", "-a",TEfile,"-b",transcriptfile,"-sorted", "-wa","-wb",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])
    else:
        intersect_list = ["bedtools", "intersect", "-a",TEfile,"-b",transcriptfile,"-sorted","-s","-wa", "-wb",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])

def intersect_not(TEfile,transcriptfile,outfile,strandedness):
    ######## intersect WITH BED FILE #########################
    #keep read if 50% of read overlaps with TE range
    if strandedness == 0:
        intersect_list = ["bedtools", "intersect", "-a",TEfile,"-b",transcriptfile,"-sorted", "-v",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])
    else:
        intersect_list = ["bedtools", "intersect", "-a",TEfile,"-b",transcriptfile,"-sorted","-s","-v",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])


def intersect_loj(file1,file2,outfile,strandedness):
    ######## intersect WITH BED FILE #########################
    #keep read if 50% of read overlaps with TE range
    if strandedness == 0:
        intersect_list = ["bedtools", "intersect", "-a",file1,"-b",file2,"-sorted", "-loj",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])
    else:
        intersect_list = ["bedtools", "intersect", "-a",file1,"-b",file2,"-sorted","-s","-loj",">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])


def merge(infile,outfile,distance,strandedness, function):
    if strandedness == 0:     
        intersect_list = ["bedtools", "merge", "-d",str(distance),"-i",infile]
        intersect_list += function
        intersect_list += ["|","awk","-v", "OFS='\\t'", """'{ print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS "." }'"""]
        intersect_list += [">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])
    else:     
    #only merge on same strand   
        intersect_list = ["bedtools", "merge", "-d",str(distance),"-i",infile, "-s"]
        intersect_list+= function
        intersect_list += ["|","awk","-v", "OFS='\\t'","""'{ print $1 OFS $2 OFS $3 OFS $5 OFS $6 OFS $4 }'"""]
        intersect_list += [">",outfile]
        intersect_command = " ".join(intersect_list)
        sp.check_call(["/bin/sh", "-c", intersect_command])        


# def make_bed(infile,outfile,strandedness):
#     if strandedness >0:
#         replacecommand_list = ["awk","-v", "OFS='\\t'", """'{ print $1 OFS $2 OFS $3 OFS $5 OFS $6 OFS $4 }'""",infile, ">" ,outfile]
#         replacecommand=" ".join(replacecommand_list)
#         sp.check_call(["/bin/sh","-c",replacecommand])
#     else:
#         replacecommand_list = ["awk","-v", "OFS='\\t'", """'{ print $1 OFS $2 OFS $3 OFS $5 OFS $6 OFS  "." }'""",infile, ">" ,outfile]
#         replacecommand=" ".join(replacecommand_list)
#         sp.check_call(["/bin/sh","-c",replacecommand])


def get_gtf_table(outfolder,outgtf,filename):
    gtf_table = make_tempfile(filename,"gtf_table",outfolder)
    gtf_list = open(gtf_table,'w')
    gtf_path=outgtf
    gtf_list.writelines(filename + "\t" + gtf_path + "\n")
    gtf_list.close()
    return  gtf_table

def find_splice_sites(gtf_file,start_dict,stop_dict):       
    with open(gtf_file,'r') as infile:
        for line in infile:
            line=line.rstrip().split("\t")
            gtf_line = gtfline(line)
            if gtf_line.category=="exon":
                start_dict[(gtf_line.chrom, str(gtf_line.start),gtf_line.strand)].add(gtf_line.Gene_ID)
                stop_dict[(gtf_line.chrom,str(gtf_line.stop),gtf_line.strand)].add(gtf_line.Gene_ID)


def filter_splice(infile,splice_start,splice_stop, start_dict,stop_dict):
    motif_dict={"0":"non-canonical","1":"GT/AG","2":"CT/AC","3":"GC/AG", "4":"CT/GC", "5":"AT/AC", "6": "GT/AT"}
    startout = open(splice_start,'w')
    stopout = open(splice_stop,'w')
    with open(infile,'r') as filterin:
        for line in filterin:
            line = line.rstrip()
            line = line.split("\t")
            chrom = line[0]
            strand_info  = line[3]
            if strand_info=="0":
                strand = "."
            elif strand_info =="1":
                strand="+"
            elif strand_info =="2":
                strand = "-"
            start = str(int(line[1])-1)
            stop = str(int(line[2]))
            motif = motif_dict[line[4]]
            if line[5]=="0":
                start_annot="novel"
                stop_annot = "novel"
            elif line[5]=="1":
                if (chrom,start,strand) in stop_dict:
                    stop_annot = stop_dict[(chrom,start,strand)] 
                    stop_annot=",".join([str(i) for i in stop_annot])
                else:
                    stop_annot = "novel"
                if (chrom,stop,strand) in start_dict:
                    start_annot = start_dict[(chrom,stop,strand)]
                    start_annot=",".join([str(i) for i in start_annot])
                else:
                    start_annot = "novel"
            unique = line[6]
            if int(unique) < 1:
                continue
            multi = line[7]
            overhang = line[8]
            name= start_annot  + "-" + stop_annot 
            score = motif + "/" + unique + "/" + multi + "/" + overhang
            splice_ID = start + "-" + stop + "|" + name + "|" + score  + "|" +  strand               
            stopline = "\t".join([chrom, start,start, splice_ID, score, strand])
            startline = "\t".join([chrom, stop,stop, splice_ID, score, strand])
            startout.writelines(startline + "\n")
            stopout.writelines(stopline + "\n")
    startout.close()
    stopout.close()


def find_spliced_TEs(spliced_file,splicedTE_dict,splice_site_dict,merged_dict):    
    with open(spliced_file,'r') as infile:
        for line in infile:
            line=line.rstrip()
            line=line.split("\t")
            TE_line=line[0:10]
            splice_line=line[10:16]
            TE_start=int(TE_line[1])
            TE_stop = int(TE_line[2])
            TE_ID=TE_line[3]
            # if strandedness==0:
            #     TE_strand="."
            # else:
            #     TE_strand=TE_line[5]
            TE_strand=TE_line[5]
            splice_ID=splice_line[3]            
            intron_start = int(splice_ID.split("|")[0].split("-")[0])
            intron_stop = int(splice_ID.split("|")[0].split("-")[1])      
            splice_site_dict[splice_ID].add((TE_ID,TE_strand))
            if TE_start < intron_start and TE_stop > intron_stop:
                splicedTE_dict[(TE_ID,TE_strand)].add((splice_ID,False))
            else:
                splicedTE_dict[(TE_ID,TE_strand)].add((splice_ID,True))  

def get_exons(gtf,exonfile):
    exonfile_presort=exonfile+"pre"
    awk_list = ["awk","-v", "OFS='\\t'", """'$3=="exon"'""",gtf,">",exonfile_presort]
    awk_command = " ".join(awk_list)
    sp.check_call(["/bin/sh", "-c", awk_command])  
    sort_coord(exonfile_presort,exonfile,1,4)
    os.unlink(exonfile_presort)

def get_exp_TE(TEfile,expfile):
    awk_list = ["awk","-v", "OFS='\\t'", """'$5 > 0.01 && $8 > 5'""",TEfile,">",expfile]
    awk_command = " ".join(awk_list)
    sp.check_call(["/bin/sh", "-c", awk_command]) 


def tx_loop(outdis,gene_dict,tx_gene_dict,splicedTE_dict,splice_site_dict,merged_dict,basename,strandedness,tempflag):
    flagout=open(tempflag,'w')
    merged_out = ["contigID","TE_subfamilies","merged_counts", "merged_pos"]  #4                 
    flag_out = ["transcript_label","spliced"] #2
    TE_out=["tx_chr","tx_start","tx_stop","TE_ID","TE_fpkm","tx_strand","unique_counts", "TE_counts","aligned_reads","TE_score"] #10
    rel_gene_out = ["sense_dis"," sense_pos"," sense_init", "sense_fpkm", "antisense_dis","antisense_pos","antisense_init", "antisense_fpkm",] #8
    header_list = TE_out + rel_gene_out +  merged_out + flag_out 9D
    header = "\t".join(header_list)
    flagout.writelines(header + "\n")
    TE_dict={}
    extra_list=[]
    with open(outdis,'r') as infile:        
        for line in infile:
            line=line.rstrip().split("\t")
            if line[10]==".":
                if (line[3],line[5]) not in extra_list:
                    extra_list.append((line[3],line[5]))
                    extra_out = line[0:10] + ["."]*14
                    flagout.writelines("\t".join(extra_out) + "\n")
                continue        
            tx_class = tx_line(line,tx_gene_dict,gene_dict)                
            if (tx_class.TE_ID,tx_class.TE_strand) in merged_dict:
                merged_dict[(tx_class.TE_ID,tx_class.TE_strand)].add_tx_info(tx_class)
            if (tx_class.TE_ID,tx_class.TE_strand) in splicedTE_dict:
                splice_dict=defaultdict(set)
                other_TEs=[]
                for splice_site_tuple in splicedTE_dict[(tx_class.TE_ID,tx_class.TE_strand)]:
                    splice_ID=splice_site_tuple[0]
                    spliced_TEs = splice_site_dict[splice_ID]
                    other_TEs += [TE_tuple for TE_tuple in spliced_TEs if TE_tuple != (tx_class.TE_ID,tx_class.TE_strand)]
                    splice_dict[splice_site_tuple]=other_TEs
            else:
                spliced = False
            if (tx_class.TE_ID,tx_class.TE_strand) not in TE_dict:
                TE_dict[(tx_class.TE_ID,tx_class.TE_strand)] = TE_allinfo(line,spliced,tx_class,gene_dict,strandedness)
            else:
                TE_dict[(tx_class.TE_ID,tx_class.TE_strand)].add_transcript(tx_class,gene_dict)
    for TE,allinfo in TE_dict.iteritems():
        if TE in merged_dict:
            merged = merged_dict[TE]
        else:
            merged = False
        allinfo.label_transcript(tx_gene_dict,gene_dict,merged,filename,flagout)
    flagout.close()

class bedline(object):
    def __init(self,line):
        self.tx_chr = line[0]
        self.tx_start = int(line[1])
        self.tx_stop = int(line[2])
        self.tx_ID = line[3]
        self.tx_score=line[4]
        self.tx_strand = line[5]
        self.tx_CDS_start = int(line[6])
        self.tx_CDS_stop = int(line[7])
        self.tx_rgb = (line[8])
        self.tx_exon_count = int(line[9])
        self.tx_exon_sizes = line[10].split(",")
        self.tx_exon_starts = line[11].split(",")
        self.tx_distance = int(line[12]) 

def filter_bed(bedfile,outfile):
    with open(bedfile,'r') as infile:
        for line in infile:
            line=line.rstrip().split("\t")
            bed_line=bedline(line)
            if bed_line.tx_exon_count == 1: 
                pass

class tx_line(object):
    def __init__(self,line,tx_gene_dict,gene_dict):
        self.TE_chr = line[0]
        self.TE_start = int(line[1])
        self.TE_stop = int(line[2])
        self.TE_ID = line[3]
        self.TE_fpkm = float(line[4])
        self.TE_strand = line[5]
        self.unique_count = float(line[6])
        self.total_count = float(line[7])
        self.total_reads = float(line[8])
        self.score = float(line[9])
        self.TE_out=line[0:10]
        self.tx_chr = line[10]
        self.tx_start = int(line[11])
        self.tx_stop = int(line[12])
        self.tx_ID = line[13]
        self.tx_score=line[14]
        self.tx_strand = line[15]
        self.tx_CDS_start = int(line[16])
        self.tx_CDS_stop = int(line[17])
        self.tx_rgb = (line[18])
        self.tx_exon_count = int(line[19])
        self.tx_exon_sizes = line[20].split(",")
        self.tx_exon_starts = line[21].split(",")
        self.tx_distance = int(line[22])        
        self.fpkm_dict = {}
        self.gene_ID_tuples=tx_gene_dict[self.tx_ID]
        self.gene_IDs=set()
        for gene_ID_tuple in self.gene_ID_tuples:            
            self.gene_IDs.add(gene_ID_tuple[0])
            gene_flagout=gene_dict[gene_ID_tuple].flagout
            self.fpkm_dict[gene_ID_tuple]=gene_flagout[1]
        self.max_gene_tuple=max(self.fpkm_dict, key=self.fpkm_dict.get)            
        self.max_gene_ID = self.max_gene_tuple[0]
        self.max_fpkm=self.fpkm_dict[self.max_gene_tuple] 
        self.init = False
        self.add=False
        if self.tx_strand == self.TE_strand:
            self.rel_strand = "sense"            
        else:
            self.rel_strand = "antisense"   
            #if TE is unit and tx within 1kb of TE and NM nor NR are in transcript ID, skip 
        if self.tx_distance != 0 or (self.TE_start - self.tx_start) > 1000 or (self.tx_stop - self.TE_stop) > 1000 or self.tx_exon_count > 1: #if transcript extends beyond TE
            self.label_position_to_gene()
            self.label_start(gene_dict)        
            self.add = True
    def label_position_to_gene(self):
        if self.tx_distance < 0:
            self.tx_pos = ("upstream")
        elif self.tx_distance > 0:
            self.tx_pos = ("downstream")
        elif (self.tx_CDS_stop - self.tx_CDS_start) > 0  and self.TE_stop < int(self.tx_CDS_start):
            if self.tx_strand == "+":
                self.tx_pos = ("5'UTR")
            else:
                self.tx_pos = ("3'UTR")
        elif (self.tx_CDS_stop - self.tx_CDS_start) > 0  and self.TE_start >= int(self.tx_CDS_stop):
            if self.tx_strand == "+":
                self.tx_pos = ("5'UTR")
            else:
                self.tx_pos = ("3'UTR")
        else:
            for i in range(int(self.tx_exon_count)):
                if (self.TE_start - self.tx_start) >= int(self.tx_exon_starts[i]) and (self.TE_start - self.tx_start) <= (int(self.tx_exon_starts[i]) + int(self.tx_exon_sizes[i])):
                    self.tx_pos = ("exon")
                    break
                elif (self.TE_stop - self.tx_start) >= int(self.tx_exon_starts[i]) and (self.TE_stop - self.tx_start) <= (int(self.tx_exon_starts[i]) + int(self.tx_exon_sizes[i])):
                    self.tx_pos = ("exon")  
                    break 
                elif (self.TE_start - self.tx_start) <= int(self.tx_exon_starts[i]) and (self.TE_stop - self.tx_start) >= (int(self.tx_exon_starts[i]) + int(self.tx_exon_sizes[i])):
                    self.tx_pos = ("exon")   
                    break
                else:
                    self.tx_pos = ("intron")  
    def label_start(self,gene_dict):
        self.annot_start = int(self.TE_ID.split("|")[1])
        self.annot_stop = int(self.TE_ID.split("|")[2])
        for gene_ID_tuple in self.gene_ID_tuples:
            genedata=gene_dict[gene_ID_tuple]
            if self.tx_strand=="+" or self.tx_strand == "." :
                if genedata.start in range(self.annot_start-100,self.annot_stop+100):
                    self.init = True or self.init
            elif self.tx_strand == "-" or self.tx_strand == "." :
                if genedata.stop in range(self.annot_start-100,self.annot_stop+100):
                    self.init = True or self.init


def get_mergedinfo(mergedTE,merged_dict,merged_and_ref,merged_and_noref):
    merged_temp={}
    with open(mergedTE,'r') as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            merged_info = mergedTE_info(line)
            if (merged_info.contig_chr,merged_info.contig_start,merged_info.contig_stop,merged_info.contig_strand) not in merged_temp:
                merged_temp[(merged_info.contig_chr,merged_info.contig_start,merged_info.contig_stop,merged_info.contig_strand)] = merged_info
            else:
                merged_temp[(merged_info.contig_chr,merged_info.contig_start,merged_info.contig_stop,merged_info.contig_strand)].add_TE(line)
    with open(merged_and_ref,'r') as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            merged_chr=line[0]
            merged_start=int(line[1])
            merged_stop=int(line[2])
            merged_strand=line[5]
            exon_info=gtfline(line[6:15])
            if (merged_chr,merged_start,merged_stop,merged_strand) in merged_temp:
                merged_temp[(merged_chr,merged_start,merged_stop,merged_strand)].add_exon_info(exon_info.Gene_ID)
    with open(merged_and_noref,'r') as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            merged_chr=line[0]
            merged_start=int(line[1])
            merged_stop=int(line[2])
            merged_strand=line[5]
            exon_info=gtfline(line[6:15])
            if (merged_chr,merged_start,merged_stop,merged_strand) in merged_temp:
                merged_temp[(merged_chr,merged_start,merged_stop,merged_strand)].add_exon_info(exon_info.Gene_ID)  
    for mergedinfo in merged_temp.itervalues():
        for TE_ID in mergedinfo.contig_TE_IDs:
            merged_dict[(TE_ID,mergedinfo.contig_strand)]=mergedinfo



def label_files(file_in,file_out, string,debug):
    command = "'{print $0," + '"' + string + '"' + "}'"
    pastecommandlist = ["awk", "-v", "OFS='\\t'",command,file_in, ">", file_out]
    pastecommand = " ".join(pastecommandlist)
    sp.check_call(["/bin/sh","-c",pastecommand])
    if not debug:
        os.unlink(file_in)

def simplify_info(geneinfo_dict):
    genevaluelist=[]
    if not geneinfo_dict:
        genevaluestring="."
    else:
        for gene,values in geneinfo_dict.iteritems():            
            stringvalues = [str(i) for i in values]
            commavalues = ",".join(stringvalues)
            genevalues=gene + ":" +  commavalues
            genevaluelist.append(genevalues)
        genevaluestring = ";".join(genevaluelist)
    return genevaluestring


def simplify_info_bool(geneinfo_dict):
    genevaluelist=[]
    if not geneinfo_dict:
        genevaluestring="False"
    for gene,values in geneinfo_dict.iteritems():            
        stringvalues = str(values)            
        genevalues=gene + ":" + stringvalues
        genevaluelist.append(genevalues)
    genevaluestring = ";".join(genevaluelist)
    return genevaluestring     


class mergedTE_info(object):    
    def __init__(self,line):
        self.merged1k_set=set()
        self.TE_chr = line[0]
        self.TE_start = int(line[1])
        self.TE_stop = int(line[2])
        self.TE_ID = line[3]
        self.TE_fpkm = float(line[4])
        self.TE_strand = line[5]
        self.unique_count = float(line[6])
        self.total_count = float(line[7])
        self.total_reads = float(line[8])
        self.score = float(line[9])
        self.contig_chr = line[10]
        self.contig_start = int(line[11])
        self.contig_stop = int(line[12])        
        self.contig_cov= float(line[13])
        self.contig_ID = line[14]        
        self.contig_TE_IDs = [self.TE_ID]
        self.contig_strand = line[15]             
        self.merged_1k_chr = line[16]
        self.merged_1k_start = int(line[17])
        self.merged_1k_stop = int(line[18])
        self.merged_1k_TE_IDs=[self.TE_ID]
        self.merged_1k_counts= self.total_count
        self.merged_1k_strand = line[21]     
        self.merged1k_set.add((self.merged_1k_chr,self.merged_1k_start,self.merged_1k_stop,self.merged_1k_strand))
        self.sense_pos=defaultdict(set)   
        self.exon_genes=set()      
        splice_joined=False  
    def add_TE(self,line):
        self.merged1k_set.add((line[16],line[17],line[18],line[21]))
        self.contig_TE_IDs.append(line[3])        
    def add_exon_info(self,gene_id):
        self.sense_pos[gene_id].add("exon")
        self.exon_genes.add(gene_id)
    def add_tx_info(self,tx):
        exon = ["exon","3'UTR","5'UTR"]   
        if tx.add:
            if tx.rel_strand=="sense" :
                for gene_ID in tx.gene_IDs:
                    self.sense_pos[gene_ID].add(tx.tx_pos)
                    if tx.tx_pos in exon:
                        self.exon_genes.add(gene_ID)
    def check_unit(self):        
        self.family_set=set()
        self.taxo_set=set()        
        self.order_set = set()
        self.strand_set=set() 
        sorted_TE_IDs = set(sorted(self.contig_TE_IDs, key=lambda TE_ID: int(TE_ID.split("|")[1])))
        self.TE_count = len(sorted_TE_IDs) 
        TE_coverage=0       
        for TE_ID in self.contig_TE_IDs:
            TE_info=TE_ID.split("|")
            TE_start=TE_info[1]
            TE_stop = TE_info[2]
            TE_taxo = TE_info[3]
            TE_family=TE_taxo.split(":")[1]
            TE_order=TE_taxo.split(":")[2]
            TE_strand=TE_info[5]
            self.taxo_set.add(TE_taxo)
            self.family_set.add(TE_family)
            self.order_set.add(TE_order)
            self.strand_set.add(TE_strand) 
            TE_length=TE_stop - TE_start
            TE_coverage += TE_length
        self.contig_length = self.contig_stop - self.contig_start            
        self.pct_TE_cov = TE_coverage / self.contig_length * 100
        if len(self.merged1k_set) > 1:
            self.merged_type = "multi_TE"
        else:
            if len(self.family_set) > 1 or len(self.order_set) > 1: #if different families or orders
                self.merged_type = "multi_TE"
            elif len(self.taxo_set) > 1: #if same family but different subfamilies
                if len(self.strand_set) > 1: #if different strand orientation, not same TE
                    self.merged_type = "multi_TE"
                elif "LINE" in self.order_set or "LTR" in self.order_set: #give more flexibity for LINES and LTRs
                    self.merged_type = "ITL"
                else:
                    self.merged_type = "multi_TE"
            else:
                self.merged_type = "ITL"      
        #set TE info 
        self.taxo_set_out = ",".join([str(i) for i in self.taxo_set])        
        if len(self.contig_TE_IDs)>1:
            self.mergedID = self.contig_TE_IDs[0] +":" + self.contig_TE_IDs[-1]  + ":" + self.contig_strand
        else:
            self.mergedID = self.TE_ID
        self.contig_TE_start=self.contig_TE_IDs[0].split("|")[1]
        self.contig_TE_stop=self.contig_TE_IDs[-1].split("|")[2]
        self.merged_positions = simplify_info(self.sense_pos)         
        self.merged_out = [self.contig_start,self.contig_stop,self.contig_ID,self.taxo_set_out, self.TE_count, self.pct_TE_cov, self.contig_counts,self.merged_positions]

class merged_exon(object):
    def __init__(self,line):
        self.chrom = line[0]
        self.start = (int(line[1]))
        self.stop = int(line[2])
        self.name=(line[4])
        self.strand  = line[5]

def replace_geneid(self,newgeneid):
    newgeneid=[str(x) for x in newgeneid]        
    newgeneid=",".join(newgeneid)
    self.attributes[0] = "gene_id" + " " + '"' + newgeneid + '"'
    attributesout = self.attributes= "; ".join(self.attributes)
    gtfout = [self.chrom,self.source,self.category,self.start+1,self.stop,self.score,self.strand,self.frame,attributesout]
    self.gtfout = [str(i) for i in gtfout]


def eval_intron(bedgraph,all_exons_merged,strand, filename,tempfolder,splice_start_sorted,splice_stop_sorted,filtered_bedgraph,filtered_bed, intron_exp,interexon_exp, intron_exp_spliced_start, intron_exp_spliced_stop):
    chrom_dict={}
    with open(chrominfo,'r') as chromfile:
        for line in chromfile:
            line=line.rstrip().split("\t")
            chrom_dict[line[0]]=int(line[1])
    exon_tuple={}
    with open(all_exons_merged,'r') as infile:
        prev_chrom=False
        for line in infile:
            line=line.rstrip().split("\t")
            exonline=merged_exon(line)
            if exonline.strand == strand or exonline.strand==".":
                if not prev_chrom:
                    prev_chrom=exonline.chrom
                    exon_tuple[exonline.chrom]=[(0,exonline.start)]
                    prev_start=exonline.stop                
                elif exonline.chrom != prev_chrom:
                    exon_tuple[prev_chrom].append((prev_start,chrom_dict[prev_chrom]))
                    exon_tuple[exonline.chrom]=[(0,exonline.start)]
                    prev_start=exonline.stop   
                    prev_chrom=exonline.chrom
                else:
                    exon_tuple[exonline.chrom].append((prev_start,exonline.start))
                    prev_start = exonline.stop
        exon_tuple[prev_chrom].append((prev_start,chrom_dict[prev_chrom]))
    # for chrom in exon_tuple.keys():
    #     exon_tuple[chrom]=sorted(exon_tuple[chrom])
    block_dict={}
    outfile=open(intron_exp,'w')
    exon_count=0
    block=False
    with open(bedgraph,'r') as infile:
        for line in infile:
            line=(line.rstrip().split("\t"))
            bg_line=bgline(line)
            if bg_line.chrom not in exon_tuple:
                continue
            if not block:
                block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)] = interexon(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand) 
                block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]
            elif bg_line.chrom != block.chrom:
                exon_count=0
                block.write_block(outfile)
                block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)] = interexon(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand) 
                block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]
            while bg_line.keep == True:
                if bg_line.stop >= block.start: 
                    if bg_line.start < block.stop and bg_line.stop <= block.stop: #all bedgraph line included
                        length_start=max(bg_line.start,block.start)
                        length=bg_line.stop - length_start
                        total_count=bg_line.count * length
                        block.add_block(total_count)
                        bg_line.keep = False
                    elif bg_line.start < block.stop and bg_line.stop > block.stop: #partial bedgraph line included
                        length_start=max(bg_line.start,block.start)
                        length=block.stop- length_start
                        total_count = bg_line.count * length                    
                        block.add_block(total_count)
                        block.write_block(outfile)
                        exon_count += 1 
                        if exon_count < len(exon_tuple[bg_line.chrom]):
                            block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)] = interexon(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand) 
                            block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]                          
                    else:   
                        block.write_block(outfile)
                        bg_line.keep=False      
                        exon_count += 1  
                        block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)] = interexon(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand) 
                        block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]
                else: #if bedgraph line after block
                    bg_line.keep=False
        block.write_block(outfile) 
    outfile.close()
    outfile=open(filtered_bedgraph,'w')
    exon_count=0   
    #Check if line > mean
    with open(bedgraph,'r') as infile:
        for line in infile:
            line=(line.rstrip().split("\t"))
            bg_line=bgline(line)
            line.insert(3,strand)
            if bg_line.chrom not in exon_tuple:
                continue
            if not block:
                block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]
            elif bg_line.chrom != block.chrom:
                exon_count=0            
                block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]         
            if bg_line.stop >= block.start: 
                if bg_line.start < block.stop and bg_line.stop <= block.stop: #all bedgraph line included
                    if bg_line.count > 2 * block.mean:
                        outfile.writelines("\t".join(line)+ "\n")
                elif bg_line.start < block.stop and bg_line.stop > block.stop: # bedgraph line partially included in exon
                    outfile.writelines("\t".join(line)+ "\n")
                    if exon_count < len(exon_tuple[bg_line.chrom]):
                        block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]                    
                else:   
                    outfile.writelines("\t".join(line)+ "\n")      
                    exon_count += 1  
                    block=block_dict[(bg_line.chrom,exon_tuple[bg_line.chrom][exon_count][0],exon_tuple[bg_line.chrom][exon_count][1],strand)]
            else: #if bedgraph line after block
                outfile.writelines("\t".join(line)+ "\n")
    outfile.close()
    merge_bed(interexon_exp, filtered_bed , strand,2500)
    intersect_wawb(filtered_bed, splice_start_sorted,intron_exp_spliced_start,strandedness) #intersect TEs with intron
    intersect_wawb(filtered_bed, splice_stop_sorted,intron_exp_spliced_stop,strandedness)   


def merge_bed(bedgraph, bed, strand, distance):
    outfile=open(bed,'w')
    prev=False
    with open(bedgraph) as infile:
        for line in infile:
            line=line.rstrip()
            line=line.split("\t")
            if not prev:
                prev=contig(line,strand)
            else:
                current=contig(line,strand)
                if (current.start - prev.stop ) < distance:
                    prev = prev.add_contig(current)
                else:
                    prev.write_contig(outfile)
                    prev=current
        prev.write_contig(outfile)             
    outfile.close()


class contig(object):
    def __init__(self,line,strand):
        self.chrom = line[0]
        self.start = int(line[1])
        self.stop = int(line[2])
        self.count = float(line[3]) * (self.stop - self.start)
        self.cov = float(line[3]) 
        self.strand=strand
    def add_contig(self,contig):
        self.stop = contig.stop
        self.count +=contig.count
        self.cov +=contig.cov
    def write_contig(self,outfile):
        self.name = self.chrom + ":" + str(self.start) + "-" + str(self.stop) 
        self.bedline=[self.chrom,str(self.start),str(self.stop), "{0:.2f}".format(self.cov),self.name, self.strand]
        outfile.writelines("\t".join(self.bedline) + "\n")
    def add_splice(self,line):
        self.splice_chr = line[5]
        self.splice_start=int(line[6])

            

class bgline(object):
    def __init__(self,line):
        self.chrom = line[0]
        self.start = int(line[1])
        self.stop = int(line[2])
        self.count = float(line[3])
        self.keep = True


class interexon(object):
    def __init__(self,chrom,start,stop,strand):
        self.chrom = chrom
        self.start=start
        self.count=0
        self.stop = stop
        self.strand = strand
        self.name=self.chrom + ":" + str(self.start) + "-" + str(self.stop)
        self.mean=0
    def add_block(self,count):
        self.count += count
    def write_block(self,outfile):        
        self.length=self.stop- self.start
        if self.length > 0:
            self.mean = self.count/self.length            
            block_out=[self.chrom,str(self.start),str(self.stop),self.name,str(self.mean),self.strand]
            block_line="\t".join(block_out)
            outfile.writelines(block_line + "\n")
            

class TE_allinfo(object):    
    def __init__(self,line,spliced,tx,gene_dict,strandedness):
        #initiate variables        
        self.TE_ID = tx.TE_ID
        self.TE_start=int(self.TE_ID.split("|")[1])
        self.TE_stop=int(self.TE_ID.split("|")[2])  
        self.tx_start = tx.TE_start
        self.tx_stop = tx.TE_stop
        self.biflanked=False
        self.stop_flanked=False
        self.start_flanked=False
        self.merged_init = False
        if self.tx_stop > self.TE_stop:
            self.length_stop = self.TE_stop        
        if self.tx_start < self.TE_start:
            self.length_start = self.TE_start
        self.pct_length = (self.tx_stop-self.tx_start)/(self.TE_stop - self.TE_start)
        if strandedness > 0:
            self.TE_strand = tx.TE_strand
        else:
            self.TE_strand = "."
        self.TE_fpkm=tx.TE_fpkm
        self.TE_counts = tx.total_count
        self.TE_out = tx.TE_out
        self.exon_tx=[]
        self.sense_pos=defaultdict(set)
        self.sense_dis=defaultdict(set)
        self.sense_fpkm =defaultdict(set)
        self.sense_init = defaultdict(lambda: False)
        self.antisense_pos=defaultdict(set)
        self.antisense_dis=defaultdict(set)
        self.antisense_fpkm =defaultdict(set)
        self.antisense_init =  defaultdict(lambda: False)
        self.splice_tf=False
        self.splice_outs=set()
        self.spliced=spliced
        if self.spliced:
            for splice_site_tuple,spliced_TEs in self.spliced.iteritems():
                splice_type=splice_site_tuple[1]
                self.splice_tf = self.splice_tf or splice_type
                splice_site=splice_site_tuple[0] 
                spliced_TE_IDs=set([TE_tuple[0] for TE_tuple in spliced_TEs])
                if not spliced_TEs:
                    spliced_TE_out="NA"
                else:
                    spliced_TE_out=",".join(spliced_TE_IDs)
                splice_out=splice_site + ":" + spliced_TE_out         
                self.splice_outs.add(splice_out)
            self.spliced_out=[";".join([str(i) for i in self.splice_outs])]
        else:            
            self.splice_tf=False
            self.spliced_TEs="NA"
            self.spliced_out=["False"]
        self.add_transcript(tx,gene_dict)
    def add_transcript(self,tx,gene_dict):
        exon = ["exon","3'UTR","5'UTR"]   
        if tx.add:
            if tx.rel_strand=="sense" :
                if tx.tx_pos in exon:
                    self.exon_tx.append(tx.tx_ID)
                for gene_ID in tx.gene_IDs:
                    self.sense_pos[gene_ID].add(tx.tx_pos)                    
                    self.sense_dis[gene_ID].add(tx.tx_distance)
                    self.sense_fpkm[gene_ID].add(tx.fpkm_dict[(gene_ID,tx.tx_strand)])
                    self.sense_init[gene_ID] = self.sense_init[gene_ID] or tx.init                    
            elif tx.rel_strand == "antisense":
                for gene_ID in tx.gene_IDs:      
                    self.antisense_pos[gene_ID].add(tx.tx_pos)
                    self.antisense_dis[gene_ID].add(tx.tx_distance)
                    self.antisense_fpkm[gene_ID].add(tx.fpkm_dict[(gene_ID,tx.tx_strand)])
                    self.antisense_init[gene_ID] = self.antisense_init[gene_ID] or tx.init   
    def label_transcript(self,tx_gene_dict,gene_dict,merged,basename,outfile):         
        self.mergedinfo = merged  
        if self.mergedinfo:
            self.mergedinfo.check_unit()   
            self.merged_tx_start=int(self.mergedinfo.contig_start)
            self.merged_tx_stop=int(self.mergedinfo.contig_stop)
            self.merged_TE_start=int(self.mergedinfo.contig_TE_start)
            self.merged_TE_stop=int(self.mergedinfo.contig_TE_stop)
            if self.mergedinfo.contig_strand=="+" or self.mergedinfo.contig_strand == "." :
                if self.TE_ID == self.mergedinfo.contig_TE_IDs[0] and self.merged_tx_start >= self.merged_TE_start:
                    self.merged_init = True or self.merged_init
            if self.mergedinfo.contig_strand=="-" or self.mergedinfo.contig_strand == "." :
                if self.TE_ID == self.mergedinfo.contig_TE_IDs[-1]and self.merged_tx_stop <= self.merged_TE_stop:
                    self.merged_init = True or self.merged_init
            if self.merged_tx_stop > self.merged_TE_stop:
                self.length_stop = self.TE_stop
                self.stop_flanked=True
            else:
                self.length_stop = self.tx_stop
            if self.merged_tx_start < self.merged_TE_start:
                self.length_start = self.TE_start
                self.start_flanked=True
            else:
                self.length_start = self.tx_start
            if self.start_flanked and self.stop_flanked:
                self.biflanked = True  
        if self.exon_tx: 
            if not any(basename in tx_ID for tx_ID in self.exon_tx):
                self.transcript_label = "RefGene"            
            else:
                self.transcript_label = "StringTie"
        elif self.mergedinfo and self.mergedinfo.exon_genes:
            if not any(basename in gene_ID for gene_ID in self.mergedinfo.exon_genes):
                self.transcript_label = "RefGene"            
            else:
                self.transcript_label = "StringTie"          
        elif self.mergedinfo and self.mergedinfo.merged_type == "multi_TE":
            self.transcript_label = "SQuIRE"
        elif self.mergedinfo and self.biflanked==True:
            self.transcript_label = "SQuIRE"
        elif self.mergedinfo and self.spliced and self.splice_tf:
            self.transcript_label="SQuIRE"
        elif self.mergedinfo and self.TE_fpkm > 0.01 and self.TE_counts > 10:
            self.transcript_label = "ITL"
        else:
            self.transcript_label = "low_exp"
        flag_out = [self.transcript_label]
        gene_out = []
        gene_out.append(simplify_info(self.sense_dis))
        gene_out.append(simplify_info(self.sense_pos))
        gene_out.append(simplify_info_bool(self.sense_init))        
        gene_out.append(simplify_info(self.sense_fpkm))
        gene_out.append(simplify_info(self.antisense_dis))
        gene_out.append(simplify_info(self.antisense_pos))
        gene_out.append(simplify_info_bool(self.antisense_init))          
        gene_out.append(simplify_info(self.antisense_fpkm))        
        if self.mergedinfo:
            self.merged_out=self.mergedinfo.merged_out
        else:
            self.merged_out = ["."]*6     
        outline =  self.TE_out  + gene_out + self.merged_out + [str(self.merged_init)] + self.spliced_out + flag_out 
        stringout = [str(i) for i in outline]
        outfile.writelines("\t".join(stringout) + "\n")


def compare_exons (exonfile,inrefdict):
    with open(exonfile,'r') as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            norefline = gtfline(line[0:9])
            refline=gtfline(line[9:18])
            if norefline.category== "exon":
                if norefline.transcript_id not in inrefdict:
                    inrefdict[norefline.transcript_id]=set([refline.Gene_ID])
                else:
                    inrefdict[norefline.transcript_id].add(refline.Gene_ID)


def contrast_exons (exonfile,notinrefdict):
    with open(exonfile,'r') as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            norefline = gtfline(line[0:9])
            if norefline.category== "exon":
                notinrefdict.add(norefline.Gene_ID)


def find_inref(gtffile,inrefdict,notinrefdict, inref):
    inref_outfile=open(inref,'w')
    with open(gtffile,'r') as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            gtf_line = gtfline(line[0:9])
            if gtf_line.Gene_ID not in notinrefdict: #unguided transcript shares all of its exons with reference
                gtf_line.replace_geneid(inrefdict[gtf_line.transcript_id])
                inref_outfile.writelines("\t".join(gtf_line.gtfout + line[9:18]) + "\n")
    inref_outfile.close()


def find_notinref(gtffile,notinrefdict, notinref):
    notinref_outfile=open(notinref,'w')
    with open(gtffile,'r') as infile:
        for line in infile:
            if line.startswith("#"):
                continue            
            line = line.rstrip().split("\t")
            gtf_line = gtfline(line)
            if gtf_line.Gene_ID in notinrefdict:
                notinref_outfile.writelines("\t".join(line) + "\n")            
    notinref_outfile.close()


def filter_single_exon_genes(ingtf,outgtf_multiexon):
    exoncount=defaultdict(int)
    outfile=open(outgtf_multiexon,'w')
    with open(ingtf,'r') as infile:
        for line in infile:
            line=line.rstrip().split("\t")
            gtf_line=gtfline(line)
            if gtf_line.category =="exon":
                exoncount[gtf_line.Gene_ID]+=1
    with open(ingtf,'r') as infile:
        for line in infile:
            line=line.rstrip().split("\t")
            gtf_line=gtfline(line)
            if gtf_line.category =="exon":
                if exoncount[gtf_line.Gene_ID] > 1:
                    outfile.writelines("\t".join(line) + "\n")


def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = """Assembles mapped reads into transcripts, merges nearby TEs, and flags transcript type""")
        parser._optionals.title = "Arguments"
        parser.add_argument("-m","--map_folder", help = "Folder location(s) of outputs from SQuIRE Map (optional, default = 'squire_map')", type = str, metavar = "<folder>",default="squire_map")
        parser.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch'",type = str, metavar = "<folder>",default="squire_fetch")
        parser.add_argument("-c","--count_folder", help = "Folder location(s) of outputs from SQuIRE Count (optional, default = 'squire_count')", type = str, metavar = "<folder>",default="squire_count")
        parser.add_argument("-d","--draw_folder", help = "Folder location(s) of outputs from SQuIRE Count (optional, default = 'squire_draw')", type = str, metavar = "<folder>",default="squire_draw")
        parser.add_argument("-o","--flag_folder", help = "Destination folder for output files (optional; default='squire_flag')", type = str, metavar = "<folder>", default="squire_flag")
        parser.add_argument("-r","--read_length", help = "Read length", type = int, metavar = "<int>",required=True)
        parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required if more than 1 build in fetch_folder)", type=str, metavar = "<build>",default=False)
        parser.add_argument("-s","--strandedness", help = " '0' if unstranded, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard  (optional,default=1)", type = int, metavar = "<file>", default = False)
        parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
        parser.add_argument("-n","--name", help = "Basename for bam file (required if more than one bam file in map_folder)", type = str, metavar = "<str>",default=False)
        parser.add_argument("-N","--projectname", help = "Basename for project", type = str, metavar = "<str>",default=False)
        parser.add_argument("-t","--tempfolder", help = "Folder for tempfiles (optional; default=outfolder')", type = str, metavar = "<folder>", default=False)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)

        args,extra_args = parser.parse_known_args()
########## I/O #########
    ###### ARGUMENTS ######

    map_folder = args.map_folder
    count_folder = args.count_folder
    draw_folder=args.draw_folder
    outfolder=args.flag_folder
    read_length = args.read_length
    fetch_folder = args.fetch_folder    
    verbosity=args.verbosity
    pthreads = args.pthreads
    filename=args.name
    tempfolder=args.tempfolder
    projectname=args.projectname
    strandedness=args.strandedness
    build = args.build

######### StringTieT TIMING SCRIPT ############
if verbosity:
    StringTieTime = datetime.now()
    print("Script start time is:" + str(StringTieTime) + '\n', file = sys.stderr)# Prints StringTiet time
    print(os.path.basename(__file__) + '\n', file = sys.stderr) #prints script name to std err
    print("Script Arguments" + '\n' + "=================", file = sys.stderr)
    args_dict = vars(args)
    for option,arg in args_dict.iteritems():
        print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
    print("\n", file = sys.stderr)
make_dir(outfolder)
debug = True
if debug:
    debug = True

if not tempfolder:
    tempfolder=outfolder

##### GET INPUTS ####
bamfile = find_file(map_folder, ".bam",filename,1)
if not bamfile:
    raise Exception ("Bamfile with " + filename + " not found in " + map_folder + "\n")
countfile = find_file(count_folder,"_counts.txt",filename,1)
if not countfile:
    raise Exception ("countfile with " + filename + " not found in" + count_folder + "\n") 
plus_bedgraph = find_file(draw_folder,"_plus_multi.bedgraph",filename,1)
minus_bedgraph= find_file(draw_folder,"_minus_multi.bedgraph",filename,1)
unstranded_bedgraph=False
if not plus_bedgraph and not minus_bedgraph:
    unstranded_bedgraph = find_file(draw_folder,"_multi.bedgraph",filename,1)

if not unstranded_bedgraph:
    raise Exception ("bedgraphs with " + filename + " not found in" + draw_folder + "\n")   

if not plus_bedgraph:
    raise Exception ("countfile with " + filename + " not found in" + count_folder + "\n")     

gtf = find_file(fetch_folder,"_refGene.gtf",build,1)
outgtf_ref = find_file(count_folder,".gtf",filename,1)
if not gtf:
    raise Exception ("gtffile with " + filename + " not found in" + fetch_folder + "\n")    


genebed = find_file(fetch_folder,"_refGene.bed",build,1)

if not genebed:
    raise Exception ("bedfile with " + filename + " not found in"  + fetch_folder + "\n")    

splice_file=find_file(map_folder, "SJ.out.tab",filename,1)

if not splice_file:
    raise Exception ("SJ.out.tab file with " + filename + " not found in"  + map_folder+ "\n")  

chrominfo = find_file(fetch_folder,"_chromInfo.txt",build,1)     

####PREPARE STRINGTIE#####

#unguided stringtie
outgtf_noref =make_tempfile(filename, "outgtf_noref_sorted", tempfolder)
abund_noref =outgtf_noref.replace("outgtf","outabund")
outgtf_noref_temp =  make_tempfile(filename, "outgtf_noref", tempfolder)
abund_noref_temp = outgtf_noref_temp.replace("outgtf","outabund")

Stringtie(bamfile,outfolder,filename,strandedness,read_length,pthreads,False, verbosity,outgtf_noref_temp)
sort_coord(outgtf_noref_temp,outgtf_noref,1,4)
sort_coord(abund_noref_temp,abund_noref,3,5)    


if verbosity:
    print("Comparing guided and unguided transcripts for " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)

#filter_tx - find genes with any exons that don't overlap with reference
exon_ref = make_tempfile(filename, "exon_ref", tempfolder)
exon_noref = make_tempfile(filename, "exon_noref", tempfolder)
exons_notinref = make_tempfile(filename, "exons_notinref", tempfolder)
exons_inref = make_tempfile(filename, "exons_inref", tempfolder)

gtf_notinref = make_tempfile(filename, "gtf_notinref", tempfolder)
gtfvsref = make_tempfile(filename, "gtfvsref", tempfolder)
gtf_inref = make_tempfile(filename, "gtf_inref", tempfolder)

get_exons(outgtf_ref,exon_ref)
get_exons(outgtf_noref,exon_noref)


intersect_wawb(outgtf_noref,exon_ref,gtfvsref,strandedness)
intersect_not(outgtf_noref,exon_ref,exons_notinref,strandedness)
inrefdict={}
notinref_dict = set()
compare_exons(gtfvsref,inrefdict)
contrast_exons(exons_notinref,notinref_dict)

find_inref(gtfvsref,inrefdict,notinref_dict,gtf_inref)
find_notinref(outgtf_noref,notinref_dict,gtf_notinref)

gene_dict={}

genebeddict={}
if verbosity:
    print("Getting gene counts and transcript information for " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)


#filter_abund - add to gene_dict guided
filter_abund(abund_ref,gene_dict,False)
filter_abund(abund_noref,gene_dict,notinref_dict)


exp_bed = make_tempfile(filename, "exp_bed", tempfolder)
exp_bed_sorted = make_tempfile(filename, "exp_bed_sorted", tempfolder)

tx_gene_dict=defaultdict(set)
    #calculcate coverage and match transcript ids to gene ids, output exp_bed
filter_tx(outgtf_ref,genebed, exp_bed,tx_gene_dict,gene_dict,read_length,genebeddict)
filter_tx(gtf_notinref,genebed, False,tx_gene_dict,gene_dict,read_length,genebeddict)
filter_tx(gtf_inref,genebed, False,tx_gene_dict,gene_dict,read_length,genebeddict)

#gtf to bed 
outbed_noref= make_tempfile(filename, "outbed_noref", tempfolder)
gtf_to_bed(outgtf_noref,outbed_noref,debug)
outgtf_multiexon = make_tempfile(filename, "outgtf_multiexon", tempfolder)
filter_single_exon_genes(exon_noref,outgtf_multiexon)


all_exons = make_tempfile(filename, "all_exons", tempfolder)
all_exons_sorted = make_tempfile(filename, "all_exons", tempfolder)
all_exons_merged = make_tempfile(filename, "all_exons_merged", tempfolder)

combine_files(exon_ref,outgtf_multiexon,all_exons,debug)
sort_coord(all_exons,all_exons_sorted,1,4)

splice_start = make_tempfile(filename, "splice_start", tempfolder)
splice_stop = make_tempfile(filename, "splice_stop", tempfolder)
spliced_TE = make_tempfile(filename, "spliced_TE", tempfolder)
spliced_TE_start = make_tempfile(filename, "spliced_TE_start", tempfolder)
spliced_TE_stop = make_tempfile(filename, "spliced_TE_stop", tempfolder)
splice_start_sorted = make_tempfile(filename, "spliced_start_sorted", tempfolder)
splice_stop_sorted = make_tempfile(filename, "spliced_stop_sorted", tempfolder)
spliced_TE_sorted = make_tempfile(filename, "spliced_TE_sorted", tempfolder)
start_dict=defaultdict(set)
stop_dict=defaultdict(set)
find_splice_sites(gtf,start_dict,stop_dict)
filter_splice(splice_file,splice_start,splice_stop,start_dict,stop_dict) #convert coordinates of introns into bed format
sort_coord(splice_start,splice_start_sorted,1,2)     
sort_coord(splice_stop,splice_stop_sorted,1,2)      
spliced_TE_sorted = outfolder + "/" + filename + "_spliced_TEs.txt"



merge(all_exons_sorted,all_exons_merged,0,strandedness,["-c 7,3,5","-o distinct,distinct,sum"])
interexon_exp_bed= make_tempfile(filename, "interexon_exp" + "bed", tempfolder)
interexon_exp_bed_sorted= make_tempfile(filename, "interexon_exp_sorted" + "bed", tempfolder)
if strandedness > 0:
strand="+"
filtered_bedgraph_plus =make_tempfile(filename, "filtered_bedgraph" + "_" + strand, tempfolder)
filtered_bed_plus =make_tempfile(filename, "filtered_bed" + "_" + strand, tempfolder)
intron_exp_plus = make_tempfile(filename, "intron_exp" + "_" + strand, tempfolder)
interexon_exp_plus= make_tempfile(filename, "interexon_exp" + "_" + strand, tempfolder)
intron_exp_spliced_start_plus= make_tempfile(filename, "intron_exp_spliced_start" + "_" + strand, tempfolder)
intron_exp_spliced_stop_plus= make_tempfile(filename, "intron_exp_spliced_stop" + "_" + strand, tempfolder)
eval_intron(plus_bedgraph,all_exons_merged,strand, filename,tempfolder,splice_start_sorted,splice_stop_sorted, filtered_bedgraph_plus,filtered_bed_plus, intron_exp_plus,interexon_exp_plus, intron_exp_spliced_start_plus, intron_exp_spliced_stop_plus)
strand="-"
filtered_bedgraph_minus =make_tempfile(filename, "filtered_bedgraph" + "_" + strand, tempfolder)
filtered_bed_minus =make_tempfile(filename, "filtered_bed" + "_" + strand, tempfolder)
intron_exp_minus = make_tempfile(filename, "intron_exp" + "_" + strand, tempfolder)
interexon_exp_minus= make_tempfile(filename, "interexon_exp" + "_" + strand, tempfolder)
intron_exp_spliced_start_minus= make_tempfile(filename, "intron_exp_spliced_start" + "_" + strand, tempfolder)
intron_exp_spliced_stop_minus= make_tempfile(filename, "intron_exp_spliced_stop" + "_" + strand, tempfolder)
eval_intron(minus_bedgraph,all_exons_merged,strand, filename,tempfolder,splice_start_sorted,splice_stop_sorted, filtered_bedgraph_minus,filtered_bed_minus, intron_exp_minus,interexon_exp_minus, intron_exp_spliced_start_minus, intron_exp_spliced_stop_minus)
combine_files(filtered_bed_plus,filtered_bed_minus,interexon_exp_bed,debug) 
sort_coord(interexon_exp_bed,interexon_exp_bed_sorted,1,2)
else:
    strand="."
    eval_intron(unstranded_bedgraph,all_exons_merged, strand, filename,tempfolder,splice_start_sorted,splice_stop_sorted)    


    ##### PREPARE TES #####
        #filter_TE

simple_TE = make_tempfile(filename, "simple_TE", tempfolder)
simple_TE_sorted = make_tempfile(filename, "simple_TE_sorted", tempfolder)
simple_tx = make_tempfile(filename, "simple_tx", tempfolder)
simple_tx_sorted = make_tempfile(filename, "simple_tx_sorted", tempfolder)    
exp_TE = make_tempfile(filename, "exp_TE", tempfolder)
exp_TE_sorted = make_tempfile(filename, "exp_TE_sorted", tempfolder)   

if verbosity:
    print("Merging expressed TEs" + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)

intersect_wawb(simple_TE_sorted, splice_start_sorted,spliced_TE_start,strandedness) #intersect TEs with intron
intersect_wawb(simple_TE_sorted, splice_stop_sorted,spliced_TE_stop,strandedness)   
combine_files(spliced_TE_start,spliced_TE_stop,spliced_TE,debug) 
sort_coord(spliced_TE,spliced_TE_sorted,1,2)
splicedTE_dict=defaultdict(set)
splice_site_dict=defaultdict(set)
find_spliced_TEs(spliced_TE_sorted,splicedTE_dict,splice_site_dict,merged_dict)
filter_TE(countfile,simple_TE,simple_tx,strandedness)
get_exp_TE(simple_tx,exp_TE)
sort_coord(exp_TE,exp_TE_sorted,1,2)
sort_coord(simple_TE,simple_TE_sorted,1,2)
sort_coord(simple_tx,simple_tx_sorted,1,2)
sort_coord(exp_bed,exp_bed_sorted,1,2)
#merge TE by 1000
if verbosity:
    print("Merging nearby expressed TEs for " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)


mergedTE_1kb = make_tempfile( filename, "mergedTE" + "_1kb", tempfolder)

merge(exp_TE_sorted,mergedTE_1kb,1000,strandedness,["-c 6,4,8","-o distinct,collapse,sum"])

mergedTE_5kb = make_tempfile( filename, "mergedTE" + "_5kb", tempfolder)

merge(exp_TE_sorted,mergedTE_5kb,5000,strandedness,["-c 6,4,8","-o distinct,collapse,sum"])    


if verbosity:
    print("Comparing merged TE locations with expressed transcript exons" + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)

    #intersect merge with filtered TE
TE_and_merged = make_tempfile( filename, "TE_and_merged", tempfolder)
merged1andcontig = make_tempfile( filename, "merged1andcontig", tempfolder)
merged_and_ref = make_tempfile( filename, "merged_and_ref", tempfolder)
merged_and_noref = make_tempfile( filename, "merged_and_noref", tempfolder)    
intersect_wawb(interexon_exp_bed_sorted,mergedTE_1kb,merged1andcontig,strandedness)
intersect_wawb(simple_tx_sorted,merged1andcontig,TE_and_merged,strandedness)
intersect_wawb(interexon_exp_bed_sorted,exon_ref,merged_and_ref,strandedness)
intersect_wawb(interexon_exp_bed_sorted,exon_ref,merged_and_noref,strandedness)
merged_dict={}
get_mergedinfo(TE_and_merged,merged_dict,merged_and_ref,merged_and_noref)

if verbosity:
    print("Comparing TE locations with spliced intron locations " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)

    #find spliced TEs


####COMBINE INFO ####
#get closest sense 
#merge single exon unguided transcripts
if verbosity:
    print("Finding closest sense gene " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)
sensedis_noref = make_tempfile(filename,"sensedis_noref",tempfolder)
antisensedis_noref = make_tempfile(filename,"antisensedis_noref",tempfolder)
noref_dis_unlabeled = make_tempfile(filename,"norefdis",tempfolder) 
noref_dis = make_tempfile(filename,"norefdis",tempfolder) 

closest_sense(simple_tx_sorted,outbed_noref,sensedis_noref,strandedness)
closest_antisense(simple_tx_sorted,outbed_noref,antisensedis_noref,strandedness)
combine_files(sensedis_noref,antisensedis_noref,noref_dis,debug) 
#label_files(noref_dis_unlabeled,noref_dis,"noref",debug) 
#get closest antisense     
if verbosity:
    print("Finding closest antisense gene " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)

sensedis_ref = make_tempfile(filename,"sensedis_ref",tempfolder)
antisensedis_ref = make_tempfile(filename,"antisensedis_ref",tempfolder)
ref_dis_unlabeled = make_tempfile(filename,"refdis",tempfolder) 
ref_dis = make_tempfile(filename,"refdis",tempfolder) 

closest_sense(simple_tx_sorted,exp_bed_sorted,sensedis_ref,strandedness)
closest_antisense(simple_tx_sorted,exp_bed_sorted,antisensedis_ref,strandedness)
combine_files(sensedis_ref,antisensedis_ref,ref_dis,debug) 
#label_files(ref_dis_unlabeled,ref_dis,"ref",debug) 

outdis = make_tempfile(filename,"outdis",tempfolder)
combine_files(ref_dis,noref_dis,outdis,debug)

# outdis_sorted = outdis + "_sorted"
# sort_coord(outdis,outdis_sorted,1,2)
#for each TE_gene intersection:
if verbosity:
    print("Flagging TEs  for " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)    

tempflag = make_tempfile(filename,"flag",tempfolder)

tx_loop(outdis,gene_dict,tx_gene_dict,splicedTE_dict,splice_site_dict,merged_dict,filename,strandedness,tempflag)
flag_file =  outfolder + "/" + filename + "_flag.txt"
### Write counts ###
if verbosity:
    print("Writing gene counts  for " + filename + " " + str(datetime.now()) + "\n",file = sys.stderr)    
gene_file=outfolder + "/" + filename + "_gene_count.txt"
gene_count = open(gene_file,'w')
count_header = ["chr","start","stop","gene_ID","counts","strand","tx_IDs"]
gene_count.writelines("\t".join(count_header) + "\n")
for gene,geneinfo in gene_dict.iteritems():
    stringout = [str(i) for i in geneinfo.countsout]
    gene_count.writelines("\t".join(stringout) + "\n")


gene_count.close()

if not debug:
    os.unlink(outgtf_ref_temp)
    os.unlink(outgtf_noref_temp)    
    os.unlink(abund_ref_temp)
    os.unlink(abund_noref_temp)
    os.unlink(exon_ref)
    os.unlink(exon_noref)
    os.unlink(exon_notinref)
    os.unlink(inref)
    os.unlink(gtf_notinref)
    os.unlink(gtfvsref)
    os.unlink(gtf_inref)
    os.unlink(exp_bed)
    os.unlink(exp_bed_sorted)
    os.unlink(outbed_noref)
    os.unlink(simple_TE)
    os.unlink(simple_TE_sorted)
    os.unlink(simple_tx)
    os.unlink(simple_tx_sorted)
    os.unlink(exp_TE)
    os.unlink(exp_TE_sorted)
    os.unlink(mergedTE_1kb)
    os.unlink(mergedTE_1kb_bed)
    os.unlink(TE_and_merged)
    os.unlink(merged_and_ref)
    os.unlink(merged_and_noref)
    os.unlink(splice_bed)    
    os.unlink(spliced_TE)
    os.unlink(splice_bed_sorted)
    os.unlink(sensedis_noref)
    os.unlink(antisensedis_noref)
    os.unlink(noref_dis_unlabeled)
    os.unlink(noref_dis)
    os.unlink(sensedis_ref)
    os.unlink(antisense_ref)
    os.unlink(ref_dis_unlabeled)
    os.unlink(ref_dis)
    os.unlink(outdis)    
    os.unlink(tempflag)


####### STOP TIMING SCRIPT #######################
if verbosity:
    print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

    endTime = datetime.now()
    print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
    print('it took: ' + str(endTime-StringTieTime) + "\n", file = sys.stderr)

###################
if __name__ == "__main__":
main()









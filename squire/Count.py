#!/usr/bin/env python

############MODULES#########################
from __future__ import print_function,division
import sys
import os
import errno
import re
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

###################################################
#########INITIATE VARIABLES #####################################################

RepCalc_dict = {}
subF_reads = defaultdict(int)

####### FUNCTIONS ##############################
def isempty(filepath):
	if os.path.getsize(filepath) == 0:
		raise Exception(filepath + " is empty")

def get_basename(filepath):
		filename = os.path.basename(filepath)
		filebase = os.path.splitext(filename)[0]
		return filebase

def make_dir(path):
	try:
		original_umask = os.umask(0)
		os.makedirs(path, 0770)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
	finally:
		os.umask(original_umask)

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


def rename_file(oldname,newname):
	shutil.move(oldname, newname)

	####create tempfiles ###
def make_tempfile(basename, step, outfolder):
	tmpfile = tempfile.NamedTemporaryFile(delete=False, dir = outfolder, prefix= basename + "_" + step +  ".tmp")
	tmpname = tmpfile.name
	tmpfile.close()
	return tmpname

def getlibsize(logfile, infile,multi_bed,uniq_bed,paired_end,debug):
	if logfile:
		STAR_logfile=open(logfile,'r')
		for line in STAR_logfile:
			line = line.strip()
			unique_string = """Uniquely mapped reads number"""
			multi_string = """Number of reads mapped to multiple loci"""
			if unique_string in line:
				unique_libsize  = int(re.search("\d+",line).group(0))
			elif multi_string in line:
				multi_libsize =int(re.search("\d+",line).group(0))
		libsize = (unique_libsize + multi_libsize)
		STAR_logfile.close()
	else:
		count_temp = infile + "libsize"
		linecountcommandlist = ["samtools", "view", infile, "|", "cut", "-f1", "|", "sort", "-k1,1", "|" , "uniq","|", "wc -l", ">", count_temp]
		linecountcommand = " ".join(linecountcommandlist)
		sp.check_call(["/bin/bash","-c",linecountcommand])
		with open(count_temp, 'r') as count_file:
			first_line = count_file.readline()
			first_line_split = first_line.split()
			libsize = int(first_line_split[0])
		if paired_end:
			libsize = libsize
		if not debug:
			os.unlink(count_temp)
	return libsize


def getlinecount(first_file,name):
	count_temp = first_file +"_" + name +  ".libsize"
	linecountcommandlist = ["wc","-l",first_file,">", count_temp]
	linecountcommand = " ".join(linecountcommandlist)
	sp.check_call(["/bin/bash","-c",linecountcommand])

	with open(count_temp, 'r') as count_file:
		first_line = count_file.readline()
		first_line_split = first_line.split()
		libsize = first_line_split[0]
		return int(libsize)

	# os.unlink(count_temp)


def Stringtie(bamfile,outfolder,basename,strandedness,pthreads,gtf, verbosity,outgtf):
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
        max_multi_pct = .95
        min_coverage = 2.5
        TEoptions = [stringtie_strand, "-f",str(pct_max_fpkm),"-m", str(min_tx_length), "-a", str(flanklength), "-j", str(flankdepth), "-g", str(read_gap), "-M", str(max_multi_pct), "-c", str(min_coverage), "-e"]
    else:
        inputs = [bamfile]
        pct_max_fpkm=0.1
        flanklength = 10
        flankdepth = .1
        read_gap = 50
        min_tx_length=200
        max_multi_pct = 1.0
        min_coverage = 1.5
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
    sp.check_call(["/bin/bash", "-c", StringTiecommand])

			
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

def filter_tx(infile,gene_dict,read_length,genecounts):    
	with open(infile,'r') as filterin:
		header=filterin.readline()
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
			else:
				ref_line=gtfline(line[9:18])
				gene_dict[(ref_line.Gene_ID,ref_line.strand)].add_tx(gtf_line.transcript_id)
	with open(genecounts,'w') as outfile:
		for genestrand,geneinfo in gene_dict.iteritems():
			outline="\t".join(geneinfo.countsout)
			outfile.writelines(outline+"\n")

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
		self.countsout=[self.chrom,self.start,self.stop,self.Gene_ID,self.fpkm,self.strand,int(round(self.counts)),self.tx_ID_string]
		self.countsout = [str(i) for i in self.countsout]
	def add_counts(self,counts):
		self.counts += counts        
	def add_tx(self,txID):
		self.tx_IDs.add(txID)
		self.tx_ID_string=",".join(self.tx_IDs)
		self.flagout = [self.Gene_ID,self.fpkm,self.counts]
		self.countsout = [self.chrom,self.start,self.stop,self.Gene_ID,self.fpkm,self.strand,int(round(self.counts)),self.tx_ID_string]
		self.countsout = [str(i) for i in self.countsout]

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


def intersect(bamfile,bedfile,out_bed):
	######## INTERSECT WITH BED FILE #########################
	intersect_list = ["bedtools", "intersect", "-a",bamfile,"-b",bedfile,"-wo", "-bed",">",out_bed]
	intersect_command = " ".join(intersect_list)
	sp.check_call(["/bin/bash", "-c", intersect_command])

def intersect_flank(bamfile,bedfile,out_bed,debug):
	######## INTERSECT WITH BED FILE #########################
	#keep read if 50% of read overlaps with TE range
	intersect_list = ["bedtools", "intersect", "-a",bamfile,"-b",bedfile,"-wo", "-bed","-f", ".5",">",out_bed]
	intersect_command = " ".join(intersect_list)
	sp.check_call(["/bin/bash", "-c", intersect_command])

def label_files(file_in,file_out, string,debug):
	command = "'{print $0," + '"' + string + '"' + "}'"
	pastecommandlist = ["awk", "-v", "OFS='\\t'",command,file_in, ">", file_out]
	pastecommand = " ".join(pastecommandlist)
	sp.check_call(["/bin/bash","-c",pastecommand])
	if not debug:
		os.unlink(file_in)
def combine_files(file1,file2,outfile,debug):
	catcommand_list = ["cat", file1, file2, ">", outfile] #combines multi_aligned reads
	catcommand = " ".join(catcommand_list)
	sp.check_call(["/bin/bash","-c",catcommand])

	if not debug:
		os.unlink(file1)
		os.unlink(file2)

def sort_temp(tempfile, field,sorted_tempfile,debug):
	field_command = str(field) + "," + str(field)
	sort_command_list = ["sort","-k",field_command, tempfile, ">", sorted_tempfile]
	sort_command = " ".join(sort_command_list)
	sp.check_call(["/bin/bash", "-c", sort_command])
	if not debug:
		os.unlink(tempfile)

def get_header(bamfile,headerfile):
	samtoolscommand_list = ["samtools","view","-H", bamfile, ">",headerfile]
	samtoolscommand = " ".join(samtoolscommand_list)
	sp.check_call(["/bin/bash", "-c", samtoolscommand])

def is_paired(bamfile,basename,tempfolder,debug):
	bam_temp = make_tempfile(basename,"bam_header",tempfolder)
	get_header(bamfile,bam_temp)
	with open(bam_temp,'r') as header:
		for line in header:
			if line.startswith("@CO"):
				fastq=re.search("--readFilesIn(.+)--outFileNamePrefix",line).group(1)
				fastq_list = fastq.split()
				if len(fastq_list) > 1:
					paired = True
				else:
					paired = False
	if not debug:
		os.unlink(bam_temp)
	return paired


def find_properpair(paired_bam, proper,nonproper):
	##### FILTER INTO CONCORDANT AND DISCORDANT/SINGLE READS ####
	#-b: output in BAM format
	#-h: keep header
	#-S: input is SAM File
	#-F4: skip unmapped reads (bit flag = 4)
	#-f2 = keep proper pair
	#-F2 = discard proper pair
	samtoolscommand_list = ["samtools","view","-bf2", "-o", proper, paired_bam]
	samtoolscommand = " ".join(samtoolscommand_list)
	sp.check_call(["/bin/bash", "-c", samtoolscommand])
	samtoolscommand_list = ["samtools","view","-bF2", "-o", nonproper, paired_bam]
	samtoolscommand = " ".join(samtoolscommand_list)
	sp.check_call(["/bin/bash", "-c", samtoolscommand])

def split_paired(paired_bed, paired_bed1, paired_bed2,debug):
	#separate read 1 and read2 into separate files
	awkcommand_list = ["awk","'$4 ~ v'","v='/1'", paired_bed,">", paired_bed1]
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash", "-c", awkcommand])
	awkcommand_list = ["awk","'$4 ~ v'","v='/2'", paired_bed,">", paired_bed2]
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash", "-c", awkcommand])
	if not debug:
		os.unlink(paired_bed)

def reduce_reads(read_file,new_readfile,debug):
	#Find reads aligned to same position but different TE_IDs (overlapping flanks) and merge
	prev = False
	with open(read_file,'r') as infile:
		with open(new_readfile,'w') as outfile:
			for line in infile:
				if not prev:
					prev=bedline(line)
					prev.TE_ID = prev.line_split[15]
					prev_TE_ID = prev.TE_ID
					continue
				else:
					current = bedline(line)
					current.TE_ID = current.line_split[15]
					if current.Read_ID == prev.Read_ID and current.Read_chr == prev.Read_chr and current.Read_geno_start==prev.Read_geno_start and current.Read_geno_stop == prev.Read_geno_stop and current.Read_strand == prev.Read_strand:
						if current.TE_ID != prev.TE_ID:
							prev_TE_ID = prev_TE_ID + "&" + current.TE_ID
					else:
						prev.line_split[15]  = prev_TE_ID
						prev.line = "\t".join(prev.line_split)
						outfile.writelines(prev.line + "\n")
						prev= current
						prev_TE_ID = current.TE_ID
	#end of loop
			prev.line_split[15]  = prev_TE_ID
			prev.line = "\t".join(prev.line_split)
			outfile.writelines(prev.line + "\n")
	if not debug:
		os.unlink(read_file)


def get_coords(file_in,read_end,strandedness, file_out,debug):
	####Get genomic coordinates from bed file
	temp_file_coords = file_in + "_temp_coords"
	temp_file_chr = file_in + "_temp_chr"
	temp_file_plus = file_in + "_temp_plus"
	temp_file_minus = file_in + "_temp_minus"
	temp_file_new = file_in + "_temp_new"
	coords_commandlist = ["awk", "-v", "OFS='\\t'","""'{print $1 OFS $19-$14+$2 OFS $19-$14+$3 OFS $4 OFS $5 OFS "orig_"$6 OFS $16 OFS $23}'""",file_in, ">", temp_file_coords]
	coords_command = " ".join(coords_commandlist)
	sp.check_call(["/bin/bash","-c",coords_command])
	remove_underscore_command_list = ["awk","-v", "OFS='\\t'", """'{ gsub(/_polymorphism/,"",$1); gsub(/_novel/,"",$1);print $0 }'""", temp_file_coords, ">", temp_file_chr]
	remove_underscore_command = " ".join(remove_underscore_command_list)
	sp.check_call(["/bin/bash","-c",remove_underscore_command])
	if not debug:
		os.unlink(temp_file_coords)
		os.unlink(file_in)
	if strandedness==0:
		strandedness=1 #change strandedness just so paired-end reads are switched to the same strand
	if strandedness == read_end: #switch strand
		plus_command_list = ["awk","-v", "OFS='\\t'", """'{ gsub(/orig_\+/,"new_-",$6); print $0 }'""", temp_file_chr, ">", temp_file_plus]
		plus_command = " ".join(plus_command_list)
		sp.check_call(["/bin/bash","-c",plus_command])
		minus_command_list = ["awk", "-v", "OFS='\\t'","""'{ gsub(/orig_\-/,"new_+",$6); print $0 }'""", temp_file_plus, ">", temp_file_minus]
		minus_command = " ".join(minus_command_list)
		sp.check_call(["/bin/bash","-c",minus_command])
		new_command_list = ["awk", "-v", "OFS='\\t'","""'{ gsub("new_","",$6); print $0 }'""", temp_file_minus, ">", file_out]
		new_command = " ".join(new_command_list)
		sp.check_call(["/bin/bash","-c",new_command])
		if not debug:
			os.unlink(temp_file_chr)
			os.unlink(temp_file_plus)
			os.unlink(temp_file_minus)
	else: #keep strand
		new_command_list = ["awk","-v", "OFS='\\t'", """'{ gsub("orig_","",$6); print $0 }'""", temp_file_chr, ">", file_out]
		new_command = " ".join(new_command_list)
		sp.check_call(["/bin/bash","-c",new_command])
		if not debug:
			os.unlink(temp_file_chr)

def fix_paired(file1,file2,fixed_file1,fixed_file2, debug): #remove "/1" or "/2"
	remove1_command_list = ["sed", """'s@/1@@g'""", file1, ">", fixed_file1]
	remove1_command = " ".join(remove1_command_list)
	sp.check_call(["/bin/bash","-c",remove1_command])
	remove2_command_list = ["sed", """'s@/2@@g'""", file2, ">", fixed_file2]
	remove2_command = " ".join(remove2_command_list)
	sp.check_call(["/bin/bash","-c",remove2_command])
	if not debug:
		os.unlink(file1)
		os.unlink(file2)

def find_uniq(combined_tempfile, first_tempfile,unique_tempfile, multi_tempfile,debug):
	##### SEPARATE UNIQUELY ALIGNED AND MULTI-ALIGNED READS ########
	dupe_tempfile = combined_tempfile + "_dupe"
	dupe_tempfile1 = combined_tempfile + "_dupe1"
	dupe_command_list = ["awk","'!a[$4]++'", combined_tempfile, ">", first_tempfile] #skips lines if the read has already appeared in the file
	dupe_command = " ".join(dupe_command_list)
	sp.check_call(["/bin/bash", "-c", dupe_command])
	awkcommand_list = ["awk", """'FNR==NR{a[$0]++;next}!a[$0]'""", first_tempfile, combined_tempfile,  ">", dupe_tempfile] #writes lines in combined_tempfile that are not in unique_tempfile2 -> duplicates
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""", dupe_tempfile, first_tempfile,  ">", unique_tempfile] #writes lines where read is in unique2 but not multi file -> truly unique
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$4]++;next}a[$4]{print $0}'""", dupe_tempfile, first_tempfile,  ">", dupe_tempfile1] #writes lines in read is in unique2 and multi file -> gets first appearance of multi-aligned reads
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	#delete unneeded tempfiles
	catcommand_list = ["cat", dupe_tempfile, dupe_tempfile1, ">", multi_tempfile ] #combines multi_aligned reads
	catcommand = " ".join(catcommand_list)
	sp.check_call(["/bin/bash","-c",catcommand])
	if not debug:
		os.unlink(dupe_tempfile)
		os.unlink(dupe_tempfile1)
		os.unlink(combined_tempfile)
		os.unlink(first_tempfile)

def match_reads(R1, R2, strandedness, matched_file,unmatched_file1,unmatched_file2,debug):
	#match  read1 and read2 if within 2kb of each other on same strand
	#add rough location to read_ID to reduce combinations for join
	rounded_1_v1 = R1 + "_rounded_v1"
	rounded_2_v1 = R2 + "_rounded_v1"
	newread_1_v1 = R1 + "_newread_v1"
	newread_2_v1 = R2 + "_newread_v1"
	rounded_1_v2 = R1 + "_rounded_v2"
	rounded_2_v2 = R2 + "_rounded_v2"
	newread_1_v2 = R1 + "_newread_v2"
	newread_2_v2 = R2 + "_newread_v2"
	matched_file_v1 = matched_file + "_v1"
	matched_file_v2 = matched_file + "_v2"
	matched_file_10k_v1 = matched_file + "_10k_v1"
	matched_file_10k_v2 = matched_file + "_10k_v2"
	unmatched_file1_v1 = unmatched_file1 + "_v1"
	unmatched_file2_v1 = unmatched_file2 + "_v1"
	roundcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $2/10000}'""", """OFMT='%.f'""", R1, ">", rounded_1_v1]
	roundcommand=" ".join(roundcommand_list)
	sp.check_call(["/bin/bash","-c",roundcommand])
	roundcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $2/10000}'""", """OFMT='%.f'""", R2, ">", rounded_2_v1]
	roundcommand=" ".join(roundcommand_list)
	sp.check_call(["/bin/bash","-c",roundcommand])
	#create new read to join on that is read/chro
	newreadcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $4 "/" $1 "/" $11 "/" $6}'""", rounded_1_v1,"|", "sort -k12", ">", newread_1_v1]
	newreadcommand=" ".join(newreadcommand_list)
	sp.check_call(["/bin/bash","-c",newreadcommand])
	newreadcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $4 "/" $1 "/" $11 "/" $6}'""",  rounded_2_v1,"|", "sort -k12", ">", newread_2_v1]
	newreadcommand=" ".join(newreadcommand_list)
	sp.check_call(["/bin/bash","-c",newreadcommand])
	#use join not awk because awk only takes 1st hit with shared value to find match
	joincommand_list = ["join", "-j", "12", "-t", "$'\\t'", "-o", "1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10", newread_1_v1, newread_2_v1, ">" , matched_file_10k_v1]
	joincommand=" ".join(joincommand_list)
	sp.check_call(["/bin/bash","-c",joincommand])
	pos_strand_2 =  """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="+" && $5!=1000 && $15!=1000)""" #insert size < 500 & end of read1 will be after beginning of read 2 & start of read1 will be after beginning of read2
	minus_strand_2 = """($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="-" && $5!=1000 && $15!=1000)""" #insert size < 500 & end of read2 will be after beginning of read 1 & start of read2 will be after beginning of read1
	poly_pos_2 = """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="+" && $5==1000 && $15 !=1000) || ($13 -  $2 <= 500 && $13-  $2 >= 0 && $12 >= $2 && $6=="+" && $5!=1000 && $15==1000)||($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="+" && $5==1000 && $15==1000)"""
	poly_minus_2 =  """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="-" && $5!=1000 && $15 ==1000) || ($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="-" && $5==1000 && $15!=1000)||($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="-"  && $5==1000 && $15==1000)"""#if plus strand: pos_strand 2 before insert, minus_strand 2 after; if minus strand: minus strand_2 before insert, pos strand 2 after
	pos_strand_1 =  """($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="+" && $5!=1000)""" #insert size < 500 & end of read2 will be after beginning of read 1 & start of read2 will be after beginning of read1
	minus_strand_1 = """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="-" && $5!=1000)""" #insert size < 500 & end of read1 will be after beginning of read 2 & start of read1 will be after beginning of read2
	poly_minus_1 = """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="-" && $5==1000 && $15 !=1000) || ($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="-" && $5!=1000 && $15==1000)||($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="-" && $5==1000 && $15==1000)"""
	poly_pos_1 =  """($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12 && $6=="+" && $5!=1000 && $15 ==1000) || ($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="+" && $5==1000 && $15!=1000)||($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 && $6=="+"  && $5==1000 && $15==1000)"""#if plus strand: pos_strand 2 before insert, minus_strand 2 after; if minus strand: minus strand_2 before insert, pos strand 2 after
	unstranded = """($13 -  $2 <= 500 && $13 -  $2 >= 0 && $12 >= $2 ) || ($3 -$12 <= 500 && $3 -$12 >= 0 && $2 >= $12)"""
	awk_inout_v1 = [matched_file_10k_v1,  ">",  matched_file_v1]
	if strandedness==1: #first strand RNA synthesis (Illumina, dUTP, NSR, NNSR)
		awkcommand_list2 = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'","""'(""",pos_strand_2,"""||""",minus_strand_2,"||",poly_pos_2,"||",poly_minus_2,"){print $0}'"""]  #find pairs that match TE_ID and strand
		awkcommand2 = " ".join(awkcommand_list2 + awk_inout_v1)
		sp.check_call(["/bin/bash","-c",awkcommand2])
	if strandedness==2: #second strand (Ligation, standard Solid)
		awkcommand_list1 = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'","""'(""",pos_strand_1,"""||""", minus_strand_1 ,"||",poly_pos_1,"||",poly_minus_1,""") {print $0}'"""] #find pairs that match TE_ID and strand
		awkcommand1 = " ".join(awkcommand_list1 + awk_inout_v1)
		sp.check_call(["/bin/bash","-c",awkcommand1])
	if strandedness==0:
		awkcommand_list0 = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'","""'(""", unstranded,""") {print $0}'"""] #find pairs that match TE_ID and strand
		awkcommand0 = " ".join(awkcommand_list0 + awk_inout_v1)
		sp.check_call(["/bin/bash","-c",awkcommand0])
	awkcommand_list = ["awk","-v", "OFS='\\t'","-v", "FS='\\t'", """'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""", matched_file_v1,R1,   ">", unmatched_file1_v1] #writes lines in read1 that is not in matched file -> unmatched
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'","""'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""", matched_file_v1, R2,  ">", unmatched_file2_v1] #writes lines in read2 that is not in matched file -> unmatched
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	roundcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $2/100000}'""", """OFMT='%.f'""", unmatched_file1_v1, ">", rounded_1_v2]
	roundcommand=" ".join(roundcommand_list)
	sp.check_call(["/bin/bash","-c",roundcommand])
	roundcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $2/100000}'""", """OFMT='%.f'""", unmatched_file2_v1, ">", rounded_2_v2]
	roundcommand=" ".join(roundcommand_list)
	sp.check_call(["/bin/bash","-c",roundcommand])
	newreadcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $4 "/" $1 "/" $11 "/" $6}'""", rounded_1_v2,"|", "sort -k12", ">", newread_1_v2]
	newreadcommand=" ".join(newreadcommand_list)
	sp.check_call(["/bin/bash","-c",newreadcommand])
	newreadcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'", """'{print $0, $4 "/" $1 "/" $11 "/" $6}'""",  rounded_2_v2,"|", "sort -k12", ">", newread_2_v2]
	newreadcommand=" ".join(newreadcommand_list)
	sp.check_call(["/bin/bash","-c",newreadcommand])
	#use join not awk because awk only takes 1st hit with shared value to find match
	joincommand_list = ["join", "-j", "12", "-t", "$'\\t'", "-o", "1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10", newread_1_v2, newread_2_v2, ">" , matched_file_10k_v2]
	joincommand=" ".join(joincommand_list)
	sp.check_call(["/bin/bash","-c",joincommand])
	awk_inout_v2 = [matched_file_10k_v2,  ">",  matched_file_v2]
	if strandedness==1:
		awkcommand2 = " ".join(awkcommand_list2+ awk_inout_v2)
		sp.check_call(["/bin/bash","-c",awkcommand2])
	if strandedness==2:
		awkcommand1 = " ".join(awkcommand_list1+ awk_inout_v2)
		sp.check_call(["/bin/bash","-c",awkcommand1])
	if strandedness==0:
		awkcommand0 = " ".join(awkcommand_list0+ awk_inout_v2)
		sp.check_call(["/bin/bash","-c",awkcommand0])
	catcommand_list = ["cat", matched_file_v1, matched_file_v2, ">", matched_file ] #combines multi_aligned reads
	catcommand = " ".join(catcommand_list)
	sp.check_call(["/bin/bash","-c",catcommand])
	awkcommand_list = ["awk","-v", "OFS='\\t'","-v", "FS='\\t'", """'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""", matched_file,R1,   ">", unmatched_file1] #writes lines in read1 that is not in matched file -> unmatched
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", "-v", "OFS='\\t'","-v", "FS='\\t'","""'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""", matched_file, R2,  ">", unmatched_file2] #writes lines in read2 that is not in matched file -> unmatched
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	if not debug:
		os.unlink(rounded_1_v1)
		os.unlink(rounded_2_v1)
		os.unlink(newread_1_v1)
		os.unlink(newread_2_v1)
		os.unlink(rounded_1_v2)
		os.unlink(rounded_2_v2)
		os.unlink(newread_1_v2)
		os.unlink(newread_2_v2)
		os.unlink(matched_file_10k_v1)
		os.unlink(matched_file_10k_v2)
		os.unlink(matched_file_v1)
		os.unlink(matched_file_v2)
		os.unlink(unmatched_file1_v1)
		os.unlink(unmatched_file2_v1)
		os.unlink(R1)
		os.unlink(R2)



def merge_coords(paired_file,merged_paired,debug): #combine coordinates for paired reads
	outfile = open(merged_paired,'w')
	with open(paired_file,'r') as infile:
		for line in infile:
			line = line.rstrip()
			line_split = line.split("\t")
			chrom = line_split[0]
			R1_start = line_split[1]
			R1_end = line_split[2]
			R1_score = line_split[4]
			R2_score = line_split[14]
			R2_start = line_split[11]
			R2_end = line_split[12]
			new_read = "paired"
			R1_proper = line_split[7]
			R2_proper = line_split[17]
			R1_uniq = line_split[9]
			R2_uniq = line_split[19]
			new_start = str(min(int(R1_start), int(R2_start)))
			new_end = str(max(int(R1_end),int(R2_end)))
			read_ID = line_split[3]
			TE_ID_1 = line_split[6]
			TE_ID_2 = line_split[16]
			if TE_ID_1 != TE_ID_2:
				new_TE_ID = TE_ID_1 + "&" + TE_ID_2
				new_score=R1_score + "&" + R2_score
			else:
				new_TE_ID = TE_ID_1
				new_score = R1_score
			strand = line_split[5]
			new_uniq="R1" + "_" + R1_start + "_" + R1_uniq + ":" + "R2" + "_" + R2_start + "_" + R2_uniq
			new_proper = "R1" + "_" + R1_proper + ":" + "R2" + "_" + R2_proper
			insert_size=abs(int(new_end) - int(new_start))
			new_line = "\t".join([chrom, new_start,new_end,read_ID,new_score,strand,new_TE_ID,new_proper,new_read,new_uniq])
			outfile.writelines(new_line + "\n")
	outfile.close()
	infile.close()
	if not debug:
		os.unlink(paired_file)


def find_proper(single_bed,nonproper_bed, proper_bed,debug):
	#separate read 1 and read2 into separate files
	awkcommand_list = ["awk","'$8 ~ v'","v='nonproper'", single_bed,">", nonproper_bed]
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash", "-c", awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$0]++;next}!a[$0]{print $0}'""", nonproper_bed, single_bed, ">", proper_bed] #writes lines that are in single_bed and not in nonproper bed (all proper alignments)
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	if not debug:
		os.unlink(single_bed)

def remove_repeat_reads(paired_bed,unpaired_bed,new_unpaired_bed, debug): #removes reads from unpaired files that are already present in paired files
	removecommandlist = ["awk", "-v", "OFS='\\t'","'FNR==NR{a[$4]++;next}!a[$4]{print $0}'",paired_bed, unpaired_bed, ">", new_unpaired_bed]
	removecommand = " ".join(removecommandlist)
	sp.check_call(["/bin/bash","-c",removecommand])
	if not debug:
		os.unlink(unpaired_bed)

def find_paired_uniq(multi_tempfile,paired_uniq_tempfile,new_multi_tempfile, new_uniq_tempfile,debug):
	##### Find reads which exclusively align as paired + uniq -> means one read is unique and the other is multi-aligned to nearby TEs
	paired_uniq_tempfile_1 = paired_uniq_tempfile + "_1"
	new_multi_tempfile_1 = new_multi_tempfile + "_0"
	awkcommand_list = ["awk", """'$9=="paired" && $10 ~/uniq/'""", multi_tempfile,  ">", paired_uniq_tempfile_1] #writes lines labeled with "paired" and uniq
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$0]++;next}!a[$0]{print $0}'""", paired_uniq_tempfile_1, multi_tempfile, ">", new_multi_tempfile_1] #writes lines that are not paired and uniq
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$4]++;next}!a[$4]{print $0}'""",  new_multi_tempfile_1, paired_uniq_tempfile_1,  ">", paired_uniq_tempfile] #writes lines in read is in unique2 and multi file -> gets first appearance of multi-aligned reads
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	awkcommand_list = ["awk", """'FNR==NR{a[$0]++;next}!a[$0]{print $0}'""", paired_uniq_tempfile, multi_tempfile, ">", new_multi_tempfile] #writes lines that are not paired and uniq
	awkcommand = " ".join(awkcommand_list)
	sp.check_call(["/bin/bash","-c",awkcommand])
	new_multi_bed = open(new_multi_tempfile, 'a')
	new_uniq_bed = open(new_uniq_tempfile, 'a')

	def eval_paired_uniq(startlist,stoplist,beddict,uniqfile,multifile):
		if len(beddict)==1: #if all alignments of read are to same TE_ID
			TE_ID=beddict.keys()[0]
			newbedline = beddict[TE_ID]
			newbedline.Read_geno_start =str(max(start_list))
			newbedline.Read_geno_stop = str(min(stop_list))
			newbedline.write_bedline(uniqfile)
		else:
			for TE_ID, newbedline in beddict.iteritems():
				newbedline.write_bedline(multifile)
	paired_uniq_bed = open(paired_uniq_tempfile,'r')
	paired_uniq_bed.seek(0)
	prev_ID = False
	bed_dict={}
	start_list=[]
	stop_list=[]
	for line in paired_uniq_bed:
		bed_line = bedline(line)
		if not prev_ID: #if first line
			prev_ID = bed_line.Read_ID
			bed_dict={bed_line.TE_ID:bed_line}
			start_list = [int(bed_line.Read_geno_start)]
			stop_list = [int(bed_line.Read_geno_stop)]
		elif prev_ID == bed_line.Read_ID:
			bed_dict[bed_line.TE_ID] = bed_line
			start_list.append(int(bed_line.Read_geno_start))
			stop_list.append(int(bed_line.Read_geno_stop))
		else: #if new read, evaluate previous read and replace info
			eval_paired_uniq(start_list,stop_list,bed_dict,new_uniq_bed,new_multi_bed)
			#start new ID
			prev_ID = bed_line.Read_ID
			bed_dict={bed_line.TE_ID:bed_line}
			start_list = [int(bed_line.Read_geno_start)]
			stop_list = [int(bed_line.Read_geno_stop)]
	#End of loop
	eval_paired_uniq(start_list,stop_list,bed_dict,new_uniq_bed,new_multi_bed)
	if not debug: #delete unneeded tempfiles
		os.unlink(paired_uniq_tempfile_1)
		os.unlink(new_multi_tempfile_1)
		os.unlink(multi_tempfile)
		os.unlink(paired_uniq_tempfile)

def get_subF(TE_ID):
	TE_ID_components = TE_ID.split("|")
	TE_ID_subfamily = TE_ID_components[3]
	return TE_ID_subfamily

def get_strand(TE_ID):
	TE_ID_components = TE_ID.split("|")
	TE_ID_strand = TE_ID_components[5]
	return TE_ID_strand

def get_chr(TE_ID):
	TE_ID_components = TE_ID.split("|")
	TE_ID_strand = TE_ID_components[0]
	return TE_ID_strand

def split_subF(subF):
	if "overlap" in subF:
		subfamily_string=subF.replace("overlap:","") #remove overlap
		subfamily_string=subfamily_string.replace(".","") #remove period
		subfamily_string=subfamily_string.replace("+","") #remove strand
		subfamily_string=subfamily_string.replace("-","") #remove strand
		subfamily_list = subfamily_string.split(",")
		subfamily_list = set(subfamily_list)
		return subfamily_list
	else:
		return [subF]

def get_length(TE_ID):
	TE_ID_components = TE_ID.split("|")
	TE_ID_stop = TE_ID_components[2]
	TE_ID_start = TE_ID_components[1]
	TE_length = int(TE_ID_stop) - int(TE_ID_start)
	return TE_length

#### DEFINE CLASS TO CALCULATE COUNTS ###############

class RepCalc(object):
	def __init__(self, TE_ID):
		self.TE_ID = TE_ID
		self.uniqcounts = 0
		self.uniq_plus = 0
		self.uniq_minus = 0
		self.multilist = []
		self.multi_plus=0
		self.multi_minus = 0
		self.multi_u_plus=0
		self.multi_u_minus = 0
		self.multimax_plus=0
		self.multimax_minus = 0
		self.counts_plus = 0
		self.counts_minus = 0
		self.multi_counts = 0
		self.counts_tot = 0
		self.total_reads_tot = 0
		self.uniq_reads_plus = 0
		self.uniq_reads_minus = 0
		self.total_1_plus = 0
		self.total_1_minus = 0
		self.multi_u_reads_plus = 0
		self.multi_u_reads_minus = 0
		self.multi_reads_plus = 0
		self.multi_reads_minus = 0
		self.multi_max_reads_plus = 0
		self.multi_max_reads_minus = 0
		self.total_reads_plus = 0
		self.total_reads_minus = 0
		self.uniq_fragment_plus = 0
		self.uniq_fragment_minus = 0
		self.uniq_plus_perTagkb = 0
		self.uniq_minus_perTagkb = 0
		self.start_plus=False
		self.stop_plus=False
		self.start_minus=False
		self.stop_minus=False
		self.start_tot=False
		self.stop_tot=False
		self.efflength_plus=0
		self.efflength_minus=0
		self.length_plus=0
		self.length_minus=0
		self.length_tot=0
		self.fpkm_tot=0
		self.uniq_starts_plus=set()
		self.uniq_starts_minus=set()
		self.TEstrand = get_strand(TE_ID)
		self.read_list_plus = []
		self.read_list_minus = []
		self.total_fragment_plus = 0
		self.total_fragment_minus = 0
		self.min_fragment_plus = False
		self.min_fragment_minus = False
		self.new_counts_plus = 0
		self.new_counts_minus = 0
		self.fpkm=0
		# self.same_strandcount=0
		# self.opp_strandcount = 0
	def add_uniq(self,strand,count):
			if strand == "+":
				self.uniq_plus += count
				# self.uniq_reads_plus += count

			elif strand == "-":
				self.uniq_minus += count
				# self.uniq_reads_minus += count
	def add_multi(self,strand,read_fraction):
		count = read_fraction
		readcount=1
		if strand == "+":
			self.multi_plus += count
			# self.multi_reads_plus += readcount
		else:
			self.multi_minus += count
			# self.multi_reads_minus += readcount
	def add_multi_u(self,strand,read_fraction):
		count = read_fraction
		readcount=1
		if strand == "+":
			self.multi_u_plus += count
			# self.multi_u_reads_plus += readcount
		else:
			self.multi_u_minus += count
			# self.multi_u_reads_minus += readcount
	def add_read(self,read_ID,strand):
		if strand == "+":
			self.read_list_plus.append(read_ID)
		elif strand == "-":
			self.read_list_minus.append(read_ID)
	def get_fragment(self,start,stop, strand):
		start=int(start)
		stop = int(stop)
		length = stop - start
		if strand == "+":
			self.total_fragment_plus += length
			if not self.start_plus or start < self.start_plus:
				self.start_plus = start
			if not self.stop_plus or stop > self.stop_plus:
				self.stop_plus = stop
			if not self.min_fragment_plus:
				self.min_fragment_plus = length
			else:
				self.min_fragment_plus = min(length,self.min_fragment_plus)
			self.total_reads_plus +=1
			self.avg_fraglength_plus  = self.total_fragment_plus/self.total_reads_plus
		else:
			self.total_fragment_minus += length
			if not self.start_minus or start < self.start_minus:
				self.start_minus = start
			if not self.stop_minus or stop > self.stop_minus:
				self.stop_minus = stop
			if not self.min_fragment_minus:
				self.min_fragment_minus = length
			else:
				self.min_fragment_minus = min(length,self.min_fragment_minus)
			self.total_reads_minus +=1
			self.avg_fraglength_minus = self.total_fragment_minus/self.total_reads_minus
		if self.start_plus and self.start_minus: #minimum comparing False and a number gives False
			self.start_tot = min(self.start_plus, self.start_minus)
		elif self.start_plus:
			self.start_tot = self.start_plus
		elif self.start_minus:
			self.start_tot = self.start_minus
		if self.stop_plus and self.stop_minus: #minimum comparing False and a number gives False
			self.stop_tot = max(self.stop_plus, self.stop_minus)
			self.avg_fraglength_tot = (self.avg_fraglength_plus + self.avg_fraglength_minus)/2
		elif self.stop_plus:
			self.stop_tot = self.stop_plus
			self.avg_fraglength_tot = self.avg_fraglength_plus
		elif self.stop_minus:
			self.stop_tot = self.stop_minus
			self.avg_fraglength_tot = self.avg_fraglength_minus


	def get_uniqfragment(self,start,stop, strand, read_end, read_uniq):
		start=int(start)
		stop = int(stop)
		if read_end != "paired":
			if strand == "+":
				self.uniq_fragment_plus += 1
				if start not in self.uniq_starts_plus:
					self.uniq_starts_plus.add(start)
			else:
				self.uniq_fragment_minus += 1   #add fragment length to count
				if start not in self.uniq_starts_minus:
						self.uniq_starts_minus.add(start)
		else: #if reads are paired
			uniq_list = read_uniq.split(":")
			R1_uniq_list = uniq_list[0].split("_")
			R2_uniq_list = uniq_list[1].split("_")
			R1_start = R1_uniq_list[1]
			R2_start = R2_uniq_list[1]
			R1_uniq=R1_uniq_list[2]
			R2_uniq=R2_uniq_list[2]
			# Effective uniquely alignable length
			if R1_uniq == "uniq":
				if strand == "+":
					self.uniq_fragment_plus += 1   #add to fragment count
					self.uniq_starts_plus.add(R1_start)
				if strand == "-":
					self.uniq_fragment_minus += 1  #add to fragment count
					self.uniq_starts_minus.add(R1_start)
			if R2_uniq == "uniq":
				if strand == "+":
					self.uniq_fragment_plus += 1   #add to fragment count
					self.uniq_starts_plus.add(R2_start)
				if strand == "-":
					self.uniq_fragment_minus += 1  #add to fragment count
					self.uniq_starts_minus.add(R2_start)
		### Transcript start stop
		if strand == "+":
			if not self.start_plus or start < self.start_plus:
				self.start_plus = start
			if not self.stop_plus or stop > self.stop_plus:
				self.stop_plus = stop
		else:
			if not self.start_minus or start < self.start_minus:
				self.start_minus = start
			if not self.stop_minus or stop > self.stop_minus:
				self.stop_minus = stop
		if self.start_plus and self.start_minus: #minimum comparing False and a number gives False
			self.start_tot = min(self.start_plus, self.start_minus)
		elif self.start_plus:
			self.start_tot = self.start_plus
		elif self.start_minus:
			self.start_tot = self.start_minus
		self.stop_tot = max(self.stop_plus, self.stop_minus)
	def calcuniqRep(self):
		self.efflength_plus = len(self.uniq_starts_plus)
		self.efflength_minus = len(self.uniq_starts_minus)
		if self.efflength_plus:
			self.uniq_plus_perTagkb = self.uniq_fragment_plus/(int(self.efflength_plus)/1000)
		if self.efflength_minus:
			self.uniq_minus_perTagkb = self.uniq_fragment_minus/(int(self.efflength_minus)/1000)
		self.uniqcounts=self.uniq_plus + self.uniq_minus
	def calcmultiRep(self,iteration):
		self.length_plus = self.stop_plus - self.start_plus
		self.length_minus = self.stop_minus - self.start_minus
		self.length_tot = self.stop_tot - self.start_tot
		if iteration ==1:
			self.total_1_plus = self.uniq_plus + self.multi_plus + self.multi_u_plus
			self.total_1_minus = self.uniq_minus + self.multi_minus + self.multi_u_minus
			self.counts_tot = self.total_1_plus + self.total_1_minus
		tot_change = 0
		changed_strands =0
		pct_change = 0
		###old  likelihood
		if self.multi_plus:
			avg_fraglength  = self.avg_fraglength_plus
			min_fraglength = self.min_fragment_plus
			if iteration > 1:
				self.old_counts_plus = self.uniq_plus + self.multi_plus + self.multi_u_plus
				self.oldmulti_plus_perkb = self.old_counts_plus/((self.length_plus - avg_fraglength +1)/1000)
				if self.oldmulti_plus_perkb <0:
					self.oldmulti_plus_perkb = self.old_counts_plus/((self.length_plus - min_fraglength +1)/1000)
				self.multi_plus = self.multimax_plus
				self.new_counts_plus = self.uniq_plus + self.multimax_plus + self.multi_u_plus
				if self.new_counts_plus < 1:
					self.newmulti_plus_perkb = 0
				else:
					self.newmulti_plus_perkb = self.new_counts_plus/((self.length_plus - avg_fraglength +1)/1000)
					if self.newmulti_plus_perkb  < 0:
						self.newmulti_plus_perkb = self.new_counts_plus/((self.length_plus - min_fraglength +1)/1000)
				tot_change += abs(self.old_counts_plus - self.new_counts_plus)
				changed_strands +=1
				pct_change = tot_change/self.old_counts_plus *100
			else:
				self.oldmulti_plus_perkb = self.uniq_plus_perTagkb
				self.old_counts_plus = self.uniq_plus + self.multi_plus + self.multi_u_plus
				if self.old_counts_plus < 1:
					self.newmulti_plus_perkb = 0
				else:
					self.newmulti_plus_perkb = self.old_counts_plus/((self.length_plus - avg_fraglength +1)/1000)
					if self.newmulti_plus_perkb  < 0:
						self.newmulti_plus_perkb = self.old_counts_plus/((self.length_plus - min_fraglength +1)/1000)
		if self.multi_minus:
			avg_fraglength = self.avg_fraglength_minus
			min_fraglength = self.min_fragment_minus
			if iteration > 1:
				self.old_counts_minus =self.uniq_minus + self.multi_minus + self.multi_u_minus
				self.oldmulti_minus_perkb = self.old_counts_minus/((self.length_minus - avg_fraglength +1)/1000)
				if self.oldmulti_minus_perkb <0:
					self.oldmulti_minus_perkb = self.old_counts_minus/((self.length_minus - min_fraglength +1)/1000)
				self.multi_minus = self.multimax_minus
				self.new_counts_minus = self.uniq_minus + self.multimax_minus + self.multi_u_minus
				if self.new_counts_minus < 1:
					self.newmulti_minus_perkb = 0
				else:
					self.newmulti_minus_perkb = self.new_counts_minus/((self.length_minus - avg_fraglength +1)/1000)
					if self.newmulti_minus_perkb  < 0:
						self.newmulti_minus_perkb = self.new_counts_minus/((self.length_minus - min_fraglength +1)/1000)
				tot_change +=abs(self.old_counts_minus - self.new_counts_minus)
				pct_change = tot_change/self.old_counts_minus*100
				changed_strands +=1
			else:
				self.oldmulti_minus_perkb = self.uniq_minus_perTagkb
				self.old_counts_minus = self.uniq_minus + self.multi_minus + self.multi_u_minus
				if self.old_counts_minus < 1:
					self.newmulti_minus_perkb = 0
				else:
					self.newmulti_minus_perkb = self.old_counts_minus/((self.length_minus - avg_fraglength +1)/1000)
					if self.newmulti_minus_perkb  < 0:
						self.newmulti_minus_perkb = self.old_counts_minus/((self.length_minus - min_fraglength +1)/1000)
		#reset multimax
		self.multimax_plus = 0
		self.multimax_minus = 0
		if iteration > 1:
			self.counts_tot = self.new_counts_plus + self.new_counts_minus
		if changed_strands > 0:
			return (tot_change/changed_strands)
		else:
			return 0

	def add_multimax(self,strand,read_fraction):
		count = read_fraction
		if strand == "+":
			self.multimax_plus += count
			# self.multi_max_reads_plus +=1
		else:
			self.multimax_minus += count
			# self.multi_max_reads_minus +=1
	def calc_transcript_coords(self,read_locdict):
		for read in self.read_list_plus:
			start=int(read_locdict[read][(self.TE_ID,"+")][1])
			stop = int(read_locdict[read][(self.TE_ID,"+")][2])
			if not self.start_plus or start < self.start_plus:
				self.start_plus = start
			if not self.stop_plus or stop > self.stop_plus:
				self.stop_plus = stop
		for read in self.read_list_minus:
			start=int(read_locdict[read][(self.TE_ID,"-")][1])
			stop = int(read_locdict[read][(self.TE_ID,"-")][2])
			if not self.start_minus or start < self.start_minus:
				self.start_minus = start
			if not self.stop_minus or stop > self.stop_minus:
				self.stop_minus = stop
		if self.start_plus and self.start_minus: #minimum comparing False and a number gives False
			self.start_tot = min(self.start_plus, self.start_minus)
		elif self.start_plus:
			self.start_tot = self.start_plus
		elif self.start_minus:
			self.start_tot = self.start_minus
		if self.stop_plus and self.stop_minus: #minimum comparing False and a number gives False
			self.stop_tot = max(self.stop_plus, self.stop_minus)
		elif self.stop_plus:
			self.stop_tot = self.stop_plus
		elif self.stop_minus:
			self.stop_tot = self.stop_minus
	def calc_total_reads(self):
		self.total_reads_plus = len(self.read_list_plus)
		self.total_reads_minus = len(self.read_list_minus)
	def writeRep(self,aligned_libsize, counts_file,basename,strandedness,iteration):
		if iteration ==0:
			self.total_1_plus = self.uniq_plus + self.multi_plus + self.multi_u_plus
			self.total_1_minus = self.uniq_minus + self.multi_minus + self.multi_u_minus
		self.length_plus = self.stop_plus - self.start_plus
		self.length_minus = self.stop_minus - self.start_minus
		self.length_tot = self.stop_tot - self.start_tot
		self.counts_plus = self.uniq_plus + self.multi_plus + self.multi_u_plus
		self.counts_minus =self.uniq_minus + self.multi_minus + self.multi_u_minus
		self.total_1 = self.total_1_plus + self.total_1_minus
		self.uniq_tot = self.uniq_plus + self.uniq_minus
		self.counts_tot = self.counts_plus + self.counts_minus
		self.total_reads_tot = self.total_reads_plus + self.total_reads_minus
		self.TE_ID_tab = self.TE_ID.split("|")
		self.TE_chr =self.TE_ID_tab[0]
		if strandedness> 0:
			if self.counts_plus > 0 and self.total_reads_plus > 0 and self.length_plus > 0:
				outline = squire_bed(self.TE_chr,self.start_plus,self.stop_plus, self.avg_fraglength_plus,"+",self.TE_ID,self.uniq_plus,self.counts_plus,self.total_reads_plus,basename,aligned_libsize)
				counts_file.writelines(outline.out_line + "\n")
				self.fpkm += outline.fpkm

			if self.counts_minus > 0 and self.total_reads_minus > 0 and self.length_minus > 0:
				outline = squire_bed(self.TE_chr,self.start_minus,self.stop_minus, self.avg_fraglength_minus,"-",self.TE_ID,self.uniq_minus,self.counts_minus,self.total_reads_minus,basename,aligned_libsize)
				counts_file.writelines(outline.out_line + "\n")
				self.fpkm += outline.fpkm

		else:
			if self.counts_tot > 0 and self.total_reads_tot > 0 and self.length_tot > 0:
				outline = squire_bed(self.TE_chr,self.start_tot,self.stop_tot,self.avg_fraglength_tot, ".",self.TE_ID,self.uniq_tot, self.counts_tot,self.total_reads_tot,basename,aligned_libsize)
				counts_file.writelines(outline.out_line + "\n")
				self.fpkm += outline.fpkm
		

class squire_bed(object):
	def __init__(self,chrom, start, stop, avg_fraglength,strand, TE_ID,uniq_counts,total_counts, reads, basename, aligned_libsize):
		self.seqname = chrom
		self.source = "SQuIRE"
		self.feature = "TE"
		self.start = str(start)
		self.end = str(stop)
		self.length = stop - start
		self.score = TE_ID.split("|")[4]
		self.bed = TE_ID.split("|")
		self.strand = strand
		self.uniq_counts = str(uniq_counts)
		self.total_counts = "{0:.2f}".format(total_counts)
		self.reads = str(reads)
		self.conf =  "{0:.2f}".format(total_counts/reads * 100)
		self.fpkm =  (total_counts/((self.length /1000)*(int(aligned_libsize)/1000000)))
		self.aligned_libsize = str(aligned_libsize)
		self.chrom = self.bed[0]
		self.out_list = [self.chrom, self.start, self.end, TE_ID,"{0:.2f}".format(self.fpkm),strand, basename, self.aligned_libsize] + self.bed + [self.uniq_counts,self.total_counts, self.reads, self.conf ]
		self.out_line = "\t".join(self.out_list)

class subfamily(object):
	def __init__(self,subF,multi_reads):
		self.subF = subF
		self.uniq = 0
		self.unique_copies = 0
		self.total_counts_pre = 0
		self.total_counts = 0
		self.total_reads = 0
		self.multi_reads=multi_reads
		self.fpkm=0
		# self.same=0
		# self.opp = 0
	def add_TE_count(self,RepClass,strandedness):
		self.uniq +=RepClass.uniqcounts
		self.total_counts+=RepClass.counts_tot
		self.total_counts_pre +=RepClass.total_1
		self.fpkm += RepClass.fpkm
		# if strandedness > 0:
		# 	self.same +=RepClass.same_strandcount
		# 	self.opp += RepClass.opp_strandcount
		# else:
		# 	self.same= "NA"
		# 	self.opp="NA"
	def add_copy_info(self,line):
		self.line_copies = line[1]
		self.line_length = line[2]

	def write_subfamily(self,outfile,basename,aligned_libsize,iteration):
		self.total_reads = self.multi_reads + self.uniq
		if self.total_reads >  0:
			self.conf = "{0:.2f}".format(self.total_counts/self.total_reads * 100) #confidence = total counts assigned to subfamily divided by total reads
		else:
			self.conf = "NA"
		outfile.writelines(basename + "\t" + str(aligned_libsize) + "\t" + self.subF + "\t" + self.line_copies + "\t" + str(self.fpkm)  + "\t" + str(self.uniq)  + "\t" + "{0:.2f}".format(self.total_counts) + "\t" + str(self.total_reads) + "\t" + str(self.conf) + "\n")




class bedline(object):
	def __init__(self,line):
		self.line = line.rstrip() #removes white space at end of line
		self.line_split = self.line.split('\t')  # returns list of items that were separated by tab in original file
		### Read variables #####
		col_no = len(self.line_split)
		self.Read_chr = self.line_split[0]
		self.Read_geno_start = self.line_split[1]
		self.Read_geno_stop = self.line_split[2]
		self.Read_name = self.line_split[3]
		self.Read_ID = re.split("[ #/]", self.Read_name)
		self.Read_ID=self.Read_ID[0]
		self.Read_score = self.line_split[4]
		self.Read_strand = self.line_split[5]
		if col_no >=7:
			self.TE_ID = self.line_split[6]
			self.TE_ID_list = self.TE_ID.split('&')
		else:
			self.TE_ID = False
		if col_no >=8:
			self.proper = self.line_split[7]
		else:
			self.proper = False
		if col_no >=9:
			self.Read_end = self.line_split[8]
		else:
			self.Read_end = False
		if col_no >=10:
			self.uniq = self.line_split[9]
		else: self.uniq = False
			#### TE variables #####
		self.Read_length = int(self.Read_geno_stop) - int(self.Read_geno_start)


	def write_bedline(self,outfile):
		outline = [self.Read_chr, self.Read_geno_start, self.Read_geno_stop, self.Read_name, self.Read_score,self.Read_strand]
		if self.TE_ID:
			outline.append(self.TE_ID)
		if self.proper:
			outline.append(self.proper)
		if self.Read_end:
			outline.append(self.Read_end)
		if self.uniq:
			outline.append(self.uniq)
		outline = "\t".join(outline)
		outfile.writelines(outline + "\n")

def uniquecount(tempBED,RepCalc_dict,read_locdict):
	unique_fragsum=0
	tempBED.seek(0)
	unique_linecount=0
	for line in tempBED:
		bed_line = bedline(line)
		if bed_line.Read_length < 25:
			continue
		unique_fragsum += bed_line.Read_length
		unique_linecount+=1
		if len(bed_line.TE_ID_list) == 1:
		### Convert Read coordinates from Rep chrom coordinates to genomic coordinates ##
		#if RNAseq data was aligned to whole genome, Read_start-Seq_start=0, so the result will be the same
			TE_id_list = bed_line.TE_ID.split('|')
			TE_start = TE_id_list[1]
			TE_stop = TE_id_list[2]
			if bed_line.TE_ID not in RepCalc_dict:
				RepCalc_dict[bed_line.TE_ID] = RepCalc(bed_line.TE_ID)
			if "uniq" in bed_line.uniq:
				RepCalc_dict[bed_line.TE_ID].add_uniq(bed_line.Read_strand,1)
				RepCalc_dict[bed_line.TE_ID].get_uniqfragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand,bed_line.Read_end, bed_line.uniq)
				RepCalc_dict[bed_line.TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
				RepCalc_dict[bed_line.TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
			else:
				RepCalc_dict[bed_line.TE_ID].add_multi_u(bed_line.Read_strand,1)
				RepCalc_dict[bed_line.TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
				RepCalc_dict[bed_line.TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
		else: #if two TE_IDs
			for TE_ID in bed_line.TE_ID_list:
				TE_id_list = TE_ID.split('|')
				TE_start = TE_id_list[1]
				TE_stop = TE_id_list[2]
				if TE_ID not in RepCalc_dict:
					RepCalc_dict[TE_ID] = RepCalc(TE_ID)
				if "uniq" in bed_line.uniq: #if read is unique when aligning single end, otherwise only uniquely aligned because paired
					RepCalc_dict[TE_ID].add_uniq(bed_line.Read_strand,1)
					RepCalc_dict[TE_ID].get_uniqfragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand,bed_line.Read_end, bed_line.uniq)
					RepCalc_dict[TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
					RepCalc_dict[TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
				else:
					RepCalc_dict[TE_ID].add_multi_u(bed_line.Read_strand,1)
					RepCalc_dict[TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
					RepCalc_dict[TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
	unique_fragavg=unique_fragsum/int(unique_linecount)
	return unique_fragavg

def multicount(tempBED,RepCalc_dict, multidict,read_locdict):
	tempBED.seek(0)
	for line in tempBED:
		adj = False
		bed_line = bedline(line)
		if bed_line.Read_length < 25:
			continue
		if len(bed_line.TE_ID_list) == 1: #if read not aligned to more than one TE_ID at same position
			if bed_line.TE_ID not in RepCalc_dict:
				RepCalc_dict[bed_line.TE_ID] = RepCalc(bed_line.TE_ID) #Initiate Repeat class object, Add RepeatTagNo, Repeat Total Tag Length to Repeat class object
			if bed_line.Read_ID not in multidict: #if TE_ID not in multi dictionary:
				multidict[bed_line.Read_ID]={bed_line.TE_ID:bed_line.Read_strand} #For each Read_ID, include TE it's aligned to and strand it's on
				# locdict[bed_line.Read_ID] = {bed_line.TE_ID:[bed_line.Read_chr,str(bed_line.Read_geno_start),str(bed_line.Read_geno_stop),adj]} #For each Read_Id, include location of alignment
				RepCalc_dict[bed_line.TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
				RepCalc_dict[bed_line.TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
			else: #if TE_ID in dictionary:
				multidict[bed_line.Read_ID][bed_line.TE_ID] = bed_line.Read_strand
				# locdict[bed_line.Read_ID][bed_line.TE_ID] = [bed_line.Read_chr,str(bed_line.Read_geno_start),str(bed_line.Read_geno_stop),adj]
				RepCalc_dict[bed_line.TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
				RepCalc_dict[bed_line.TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
		else: ##if read aligned to more than one TE_ID at same position
			adj= True
			for TE_ID in bed_line.TE_ID_list:
				TE_id_list = TE_ID.split('|')
				TE_start = TE_id_list[1]
				TE_stop = TE_id_list[2]

				if TE_ID not in RepCalc_dict:
					RepCalc_dict[TE_ID] = RepCalc(TE_ID)
				if bed_line.Read_ID not in multidict: #if TE_ID not in dictionary:
					multidict[bed_line.Read_ID]={TE_ID:bed_line.Read_strand} #Initiate Repeat class object, Add RepeatTagNo, Repeat Total Tag Length to Repeat class object
					# locdict[bed_line.Read_ID] = {TE_ID:[bed_line.Read_chr,str(bed_line.Read_geno_start),str(bed_line.Read_geno_stop),adj]}
					RepCalc_dict[TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)
					RepCalc_dict[TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
				else: #if TE_ID in dictionary:
					multidict[bed_line.Read_ID][TE_ID] = bed_line.Read_strand
					# locdict[bed_line.Read_ID][TE_ID] = [bed_line.Read_chr,str(bed_line.Read_geno_start),str(bed_line.Read_geno_stop),adj]
					RepCalc_dict[TE_ID].add_read(bed_line.Read_ID,bed_line.Read_strand)
					RepCalc_dict[TE_ID].get_fragment(bed_line.Read_geno_start,bed_line.Read_geno_stop,bed_line.Read_strand)


def comparedict(read_multidict, RepCalc_dict):
	for Read_ID,TE_dict in read_multidict.iteritems():
		setfraction(Read_ID,TE_dict,RepCalc_dict)


def setfraction(Read_ID,TE_dict,RepCalc_dict): #compare multi with unique
	ID_TagKb_dict = {}
	UnTagged_TEs = 0
	read_subF = []
	read_sum=0
	### Use count from  unique count based on strand
	for TE_ID,strand in TE_dict.iteritems():
		#adj=locdict[TE_ID][3]
		#if TE_ID is untagged, would not have any unique reads
		if strand == "+":
			if RepCalc_dict[TE_ID].uniq_plus_perTagkb==0:
				read_fraction = 1/len(TE_dict)
				RepCalc_dict[TE_ID].add_multi(strand,read_fraction)
				#RepCalc_dict[TE_ID].get_fragment(locdict[TE_ID][1],locdict[TE_ID][2],strand)
				UnTagged_TEs +=1
				read_subF.append(get_subF(TE_ID))
				read_sum+=read_fraction
			elif RepCalc_dict[TE_ID].uniq_plus_perTagkb:
		#if TE is tagged, evaluate likelihood of contribution by length of uniq Tags on appropriate strand per Tagkb
				ID_TagKb_dict[TE_ID] = RepCalc_dict[TE_ID].uniq_plus_perTagkb
		if strand == "-": #if TE_ID is untagged, would not have any unique reads
			if RepCalc_dict[TE_ID].uniq_minus_perTagkb==0:
				read_fraction = 1/len(TE_dict)
				read_sum+=read_fraction
				RepCalc_dict[TE_ID].add_multi(strand,read_fraction)
				#RepCalc_dict[TE_ID].get_fragment(locdict[TE_ID][1],locdict[TE_ID][2],strand)
				UnTagged_TEs +=1
				read_subF.append(get_subF(TE_ID))
			elif RepCalc_dict[TE_ID].uniq_minus_perTagkb:
		#if TE is tagged, evaluate likelihood of contribution by length of uniq Tags on appropriate strand per Tagkb
				ID_TagKb_dict[TE_ID] = RepCalc_dict[TE_ID].uniq_minus_perTagkb
	TagKb_sum = sum(itervalues(ID_TagKb_dict))
	#print("TagKb_sum " + Read_ID + " " + str(TagKb_sum),file = sys.stderr)
	remainder_fraction = (len(TE_dict) - UnTagged_TEs)/len(TE_dict) #fraction of TEs in TE_dict that are tagged
	if TagKb_sum > 0: #if some unique elements can be assigned
		for TE_ID, uniq_sum in ID_TagKb_dict.iteritems():
			# adj=locdict[TE_ID][3]
			strand = TE_dict[TE_ID]
			#print(TE_ID + " " + str(uniq_sum),file = sys.stderr)
			read_fraction = (uniq_sum/TagKb_sum) * remainder_fraction #defines read_fraction by TE_IDs contribution to total uniqperTagKb_sum
			read_sum+=read_fraction
			RepCalc_dict[TE_ID].add_multi(strand,read_fraction)
			# RepCalc_dict[TE_ID].get_fragment(locdict[TE_ID][1],locdict[TE_ID][2],strand)
			read_subF.append(get_subF(TE_ID))
	### If no unique counts for any element in new dict, read fraction is 1/number of elements per read
	else:
		for TE_ID, uniq_sum in ID_TagKb_dict.iteritems():
			strand = TE_dict[TE_ID]
			#adj=locdict[TE_ID][3]
			read_fraction = 1/len(TE_dict) #defines read_fraction by TE_IDs contribution to total uniqperTagKb_sum
			read_sum+=read_fraction
			RepCalc_dict[TE_ID].add_multi(strand,read_fraction)
			# RepCalc_dict[TE_ID].get_fragment(locdict[TE_ID][1],locdict[TE_ID][2],strand)
			read_subF.append(get_subF(TE_ID))
	unique_subF = set(read_subF) #add 1 read for each uniquely represented subfamily
	if read_sum > 1.01:
		print("read_sum is greater than 1 for " + Read_ID + "\n",file = sys.stderr)
	for subfamily in unique_subF:
		subF_list=split_subF(subfamily)
		for subF in subF_list:
			subF_reads[subF] +=1

def estdict(read_multidict, RepCalc_dict):
	changed_likelihood_sum = 0
	read_no=0
	changed_reads = 0
	for Read_ID,TE_dict in read_multidict.iteritems():
		changed_likelihood = maxfraction(Read_ID,TE_dict,RepCalc_dict)
		changed_likelihood_sum +=changed_likelihood
		read_no +=1
		if changed_likelihood >= 0.1:
			changed_reads +=1
	if read_no > 0:
		avg_changed_likelihood = changed_likelihood_sum/read_no
	else:
		avg_changed_likelihood = 0
	print("Expectation maximization changed the average likelihood by: " + str(avg_changed_likelihood) + " " + str(datetime.now())  + "\n",file = sys.stderr)
	print("Number of reads with average TE likelihoods changed by at least 0.1: " + str(changed_reads) + " " + str(datetime.now())  + "\n",file = sys.stderr)
	return avg_changed_likelihood

def maxfraction(Read_ID,TE_dict,RepCalc_dict): #calculate new likelihoods and compare with previous
	ID_TagKb_dict = {}
	UnTagged_TEs = 0
	read_subF = []
	oldTagKb_sum=0
	newTagKb_sum=0
	changed_likelihood=0
	changed_TEs = 0
	read_sum=0
	TEcount=0
	### Use count from  unique count based on strand
	for TE_ID,strand in TE_dict.iteritems():
		if strand == "+":
			ID_TagKb_dict[(TE_ID,strand)] = (RepCalc_dict[TE_ID].oldmulti_plus_perkb,RepCalc_dict[TE_ID].newmulti_plus_perkb)
			oldTagKb_sum += RepCalc_dict[TE_ID].oldmulti_plus_perkb
			newTagKb_sum += RepCalc_dict[TE_ID].newmulti_plus_perkb
			if RepCalc_dict[TE_ID].newmulti_plus_perkb > 0:
				TEcount+=1
		elif strand == "-":
			ID_TagKb_dict[(TE_ID,strand)] = (RepCalc_dict[TE_ID].oldmulti_minus_perkb,RepCalc_dict[TE_ID].newmulti_minus_perkb)
			oldTagKb_sum += RepCalc_dict[TE_ID].oldmulti_minus_perkb
			newTagKb_sum += RepCalc_dict[TE_ID].newmulti_minus_perkb
			if RepCalc_dict[TE_ID].newmulti_minus_perkb > 0:
				TEcount+=1

	for TE_ID_tuple, multi_sum_tuple in ID_TagKb_dict.iteritems():
		TE_ID = TE_ID_tuple[0]
		strand = TE_ID_tuple[1]
		old_multi = multi_sum_tuple[0]
		new_multi = multi_sum_tuple[1]
		 #defines read_fraction by TE_IDs contribution to total multipekb_sum
		if newTagKb_sum == 0:
			newread_fraction = 1/len(TE_dict)
		else:
			newread_fraction = (new_multi/newTagKb_sum)
			if newread_fraction > 1:
				raise Exception("Fraction is greater than 1:" + TE_ID + " " + Read_ID + " " + str(TEcount) + " " + str(len(TE_dict)) + " " + str(newread_fraction) + " " + str(new_multi) + " " + str(newTagKb_sum) + '\n')
		read_sum+=newread_fraction
		RepCalc_dict[TE_ID].add_multimax(strand,newread_fraction)
		if oldTagKb_sum == 0:
			oldread_fraction =  1/len(TE_dict)
		else:
			oldread_fraction = (old_multi/oldTagKb_sum)

		changed_likelihood += abs(newread_fraction-oldread_fraction)
		changed_TEs +=1

	if read_sum > 1.01:
		print("read_sum is greater than 1 for " + Read_ID + "\n",file = sys.stderr)
		print("newTagKb_sum is " + str(newTagKb_sum) + "\n",file = sys.stderr)
		print("oldTagKb_sum is " + str(oldTagKb_sum) + "\n",file = sys.stderr)
		print("read_sum is " + str(read_sum) + "\n",file = sys.stderr)
	if changed_TEs == 0:
		avg_change = 0
	else:
		avg_change = changed_likelihood/changed_TEs
	return avg_change

def sort_coord(infile, outfile,chrcol,startcol,debug):
	chrfieldsort = "-k" + str(chrcol) + "," + str(chrcol)
	startfieldsort = "-k" + str(startcol) + "," + str(startcol) + "n"
	sort_command_list = ["sort",chrfieldsort,startfieldsort, infile, ">", outfile]
	sort_command = " ".join(sort_command_list)
	sp.check_call(["/bin/bash", "-c", sort_command])
	if not debug:
		os.unlink(infile)

def sort_coord_header(infile, outfile,chrcol,startcol,debug):
	chrfieldsort = "-k" + str(chrcol) + "," + str(chrcol)
	startfieldsort = "-k" + str(startcol) + "," + str(startcol) + "n"
	sort_command_list = ["(head -n 2", infile,"&& tail -n +3", infile, "|", "sort",chrfieldsort,startfieldsort, ")", ">", outfile]
	sort_command = " ".join(sort_command_list)
	sp.check_call(["/bin/bash", "-c", sort_command])
	if not debug:
		os.unlink(infile)

def sort_counts(tempfile,headerfile,countsfile, field,debug):
	sorted_countsfile = tempfile + ".sorted"
	field_command = str(field) + "," + str(field) + "rn"
	sort_command_list = ["sort","-k",field_command, tempfile, ">", sorted_countsfile]
	sort_command = " ".join(sort_command_list)
	sp.check_call(["/bin/bash", "-c", sort_command])
	catcommand_list = ["cat", headerfile, sorted_countsfile, ">",countsfile ] #combines multi_aligned reads
	catcommand = " ".join(catcommand_list)
	sp.check_call(["/bin/bash","-c",catcommand])
	if not debug:
		os.unlink(sorted_countsfile)
		os.unlink(tempfile)
		os.unlink(headerfile)

def bedgraph(infile,strandedness,outfolder,basename):
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
	normalization=["""--outWigNorm""", "None"]
	STARcommand_list = ["STAR","""--runMode""","inputAlignmentsFromBAM"] + inputs + outputs + normalization
	STARcommand=" ".join(STARcommand_list)
	sp.check_call(["/bin/bash", "-c", STARcommand])
	if strandedness !=0:
		rename_file(plus_bedgraph_unique,outfolder + "/" + basename + "_plus_unique.bedgraph")
		rename_file(minus_bedgraph_unique,outfolder + "/" + basename + "_minus_unique.bedgraph")
		rename_file(plus_bedgraph_multi,outfolder + "/" + basename + "_plus_multi.bedgraph")
		rename_file(minus_bedgraph_multi,outfolder + "/" + basename + "_minus_multi.bedgraph")
	else:
		rename_file(bedgraph_unique,outfolder + "/" + basename + "_unique.bedgraph")
		rename_file(bedgraph_multi,outfolder + "/" + basename + "_multi.bedgraph")

def main(**kwargs):
	######## ARGUMENTS ###########
	#check if already args is provided, i.e. main() is called from the top level script
	args = kwargs.get('args', None)
	if args is None: ## i.e. standalone script called from command line in normal way
		parser = argparse.ArgumentParser(description = """Quantifies RNAseq reads aligning to TEs. Outputs TE count file and subfamily count file""")
		parser._optionals.title = "Arguments"
		parser.add_argument("-m","--map_folder", help = "Folder location of outputs from SQuIRE Map (optional, default = 'squire_map')", type = str, metavar = "<folder>",default="squire_map")
		parser.add_argument("-c","--clean_folder", help = "Folder location of outputs from SQuIRE Clean (optional, default = 'squire_clean')", type = str, metavar = "<folder>",default = "squire_clean")
		parser.add_argument("-o","--count_folder", help = "Destination folder for output files(optional, default = 'squire_count')", type = str, metavar = "<folder>", default="squire_count")
		parser.add_argument("-t","--tempfolder", help = "Folder for tempfiles (optional; default=count_folder')", type = str, metavar = "<folder>", default=False)
		parser.add_argument("-f","--fetch_folder", help = "Folder location of outputs from SQuIRE Fetch (optional, default = 'squire_fetch')",type = str, metavar = "<folder>",default="squire_fetch")
		parser.add_argument("-r","--read_length", help = "Read length (if trim3 selected, after trimming; required).", type = int, metavar = "<int>", required=True)
		parser.add_argument("-n","--name", help = "Common basename for input files (required if more than one bam file in map_folder)", type = str, metavar = "<str>",default=False)
		parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required if more than 1 build in clean_folder)", type=str, metavar = "<build>",default=False)
		parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
		parser.add_argument("-s","--strandedness", help = " '0' if unstranded eg Standard Illumina, 1 if first-strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard SOLiD (optional,default=0)", type = int, metavar = "<int>", default = 0)
		parser.add_argument("-e","--EM", help = "Run estimation-maximization on TE counts given number of times (optional, specify 0 if no EM desired; default=auto)", type=str, default = "auto")
		parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)", action = "store_true", default = False)
		args = parser.parse_args()
########## I/O #########
	###### ARGUMENTS ######
	map_folder = args.map_folder
	clean_folder = args.clean_folder
	count_folder=args.count_folder
	fetch_folder=args.fetch_folder
	read_length=args.read_length
	tempfolder=args.tempfolder
	basename = args.name
	pthreads=args.pthreads
	strandedness=args.strandedness
	EM=args.EM
	build = args.build
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

	debug = False
	######## FIND INPUTS ############
	outfolder = count_folder
	make_dir(outfolder) #Create outfolder if doesn't exist
	if not os.path.isdir(clean_folder):
		raise Exception(clean_folder + " is not a folder" )
	if not os.path.isdir(map_folder):
		raise Exception(map_folder + " is not a folder" )

	logfile = find_file(map_folder,".log",basename, 1,False)
	bamfile = find_file(map_folder,".bam",basename,1,True)
	if not bamfile:
		if basename:
			raise Exception("Cannot find bamfile matching " + basename )
		else:
			raise Exception("Cannot find bamfile in map_folder" )
	if not basename:
		basename = get_basename(bamfile)
	rmsk_bed=find_file(clean_folder,".bed",build,1,True)
	copies = find_file(clean_folder,"_copies.txt",build,1,True)
	if not rmsk_bed:
		if build:
			raise Exception("Cannot find bedfile matching " + build )
		else:
			raise Exception("Cannot find bedfile in clean_folder" )
	if not copies:
		if build:
			raise Exception("Cannot find copies.txt file matching " + build )
		else:
			raise Exception("Cannot find copies.txt file in clean_folder")

	if not tempfolder:
		tempfolder = outfolder

	paired_end = is_paired(bamfile,basename,tempfolder,debug)

	if verbosity:
		print("Quantifying Gene expression "+ str(datetime.now())  + "\n",file = sys.stderr)
	ingtf = find_file(fetch_folder,".gtf",build,1,True)
	outgtf_ref =outfolder + "/" + basename + ".gtf"
	abund_ref =outgtf_ref.replace(".gtf","_abund.txt")
	outgtf_ref_temp =  make_tempfile(basename, "outgtf_ref", tempfolder)
	abund_ref_temp = outgtf_ref_temp.replace("outgtf","outabund")
	Stringtie(bamfile,outfolder,basename,strandedness,pthreads,ingtf, verbosity,outgtf_ref_temp) 
	sort_coord(outgtf_ref_temp,outgtf_ref,1,4,debug)
	sort_coord_header(abund_ref_temp,abund_ref,3,5,debug)	       
	genename_dict={}
	filter_abund(abund_ref,genename_dict,False)
	genecounts=outfolder + "/" + basename + "_refGenecounts.txt"
	filter_tx(outgtf_ref, genename_dict,read_length,genecounts)


	#### OPEN OUTPUTS & WRITE HEADER INFORMATION#############
	if verbosity:
		print("Creating temporary files"+ str(datetime.now())  + "\n",file = sys.stderr)
	counts_temp = tempfile.NamedTemporaryFile(delete=False, dir = tempfolder, prefix="count" +  ".tmp")
	countsfilepath = outfolder + "/" + basename + "_TEcounts.txt"
	counts_file_header = open(countsfilepath +".header",'w')

	counts_file_header.writelines("tx_chr" + "\t"  + "tx_start"  + "\t" + "tx_stop"  + "\t" + "TE_ID" + "\t" + "fpkm"  + "\t" + "tx_strand" + "\t" + "Sample" + "\t" + "alignedsize" + "\t" + "TE_chr" + "\t" + "TE_start" + "\t" + "TE_stop" + "\t" + "TE_name" + "\t" + "milliDiv" + "\t" + "TE_strand" + "\t" + "uniq_counts" + "\t" + "tot_counts"  + "\t" + "tot_reads" +"\t" + "score" + "\n" )

	counts_file_header.close()


#####CREATE TEMPFILES #######
	if verbosity:
		print("Creating unique and multiple alignment bedfiles "+ str(datetime.now())  + "\n",file = sys.stderr)

	if not paired_end:
		single_bam = bamfile
		if verbosity:
			print("Intersecting bam file with TE bedfile "+ str(datetime.now())  + "\n",file = sys.stderr)
		#intersect bam files with TE bed files
		single_bed_tempfile1 = make_tempfile(basename,"single_bed_1", tempfolder)
		intersect_flank(single_bam, rmsk_bed, single_bed_tempfile1,debug)
		if verbosity:
			print("Combining adjacent TEs with same read alignment "+ str(datetime.now())  + "\n",file = sys.stderr)
		#reduce reads   #Find reads aligned to same position but different TE_IDs (overlapping flanks) and merge
		single_reduced_tempfile1 = make_tempfile(basename,"single_reduced_1", tempfolder)
		single_reduced_tempfile1_sorted =single_reduced_tempfile1  + "_sorted"
		sort_coord(single_bed_tempfile1,single_reduced_tempfile1_sorted,1,2,debug)
		reduce_reads(single_reduced_tempfile1_sorted, single_reduced_tempfile1,debug)
		if verbosity:
			print("Getting genomic coordinates of read"+ str(datetime.now())  + "\n",file = sys.stderr)
		#get genomic coordinates and RNA strand for all alignments
		single_coords_tempfile1= make_tempfile(basename,"single_coords_1", tempfolder)
		get_coords(single_reduced_tempfile1,1,strandedness,single_coords_tempfile1,debug)

		# os.unlink(single_bed_tempfile1)

		if verbosity:
			print("Identifying and labeling unique and multi reads"+ str(datetime.now())  + "\n",file = sys.stderr)
		single_labeled_tempfile1 = make_tempfile(basename,"single_labeled_1", tempfolder)
		single_labeled_tempfile2 = make_tempfile(basename,"single_labeled_2", tempfolder)
		label_files(single_coords_tempfile1,single_labeled_tempfile1,"single",debug)
		label_files(single_labeled_tempfile1,single_labeled_tempfile2,"R1",debug)

		#find unique single alignments
		first_tempfile1 = make_tempfile(basename,"first_1", tempfolder)
		unique_tempfile1 = make_tempfile(basename,"unique_1", tempfolder)
		multi_tempfile1 = make_tempfile(basename,"multi_1", tempfolder)

		find_uniq(single_labeled_tempfile2,first_tempfile1,unique_tempfile1, multi_tempfile1,debug)


		#label uniq, multi, or single
		multi_bed = make_tempfile(basename,"multi_bed", tempfolder)
		unique_bed = make_tempfile(basename,"unique_bed", tempfolder)

		label_files(unique_tempfile1, unique_bed, "uniq",debug)
		label_files(multi_tempfile1, multi_bed, "multi",debug)

		aligned_libsize = getlibsize(logfile, bamfile,multi_bed,unique_bed,paired_end,debug)

	if paired_end:
		#intersect bam files with TE bed files
		if verbosity:
			print("Identifying properly paired reads "+ str(datetime.now())  + "\n",file = sys.stderr)
		paired_bam = bamfile
		proper_bam = make_tempfile(basename,"proper_bam", tempfolder)
		nonproper_bam = make_tempfile(basename,"nonproper_bam", tempfolder)
		find_properpair(paired_bam, proper_bam,nonproper_bam)

		if verbosity:
			print("Intersecting bam files with TE bedfile "+ str(datetime.now())  + "\n",file = sys.stderr)

		proper_bed = make_tempfile(basename,"proper_bed", tempfolder)
		nonproper_bed = make_tempfile(basename,"nonproper_bed", tempfolder)
		intersect_flank(proper_bam, rmsk_bed, proper_bed,debug)
		intersect_flank(nonproper_bam, rmsk_bed, nonproper_bed,debug)

		proper_labeled_tempfile = make_tempfile(basename,"proper_labeled", tempfolder)
		nonproper_labeled_tempfile = make_tempfile(basename,"nonproper_labeled", tempfolder)
		label_files(proper_bed,proper_labeled_tempfile,"proper",debug)
		label_files(nonproper_bed,nonproper_labeled_tempfile,"nonproper",debug)

		proper_nonproper_labeled_tempfile = make_tempfile(basename,"proper_nonproper_labeled", tempfolder)
		combine_files(proper_labeled_tempfile,nonproper_labeled_tempfile,proper_nonproper_labeled_tempfile,debug)

		if verbosity:
			print("Splitting into read1 and read 2 "+ str(datetime.now())  + "\n",file = sys.stderr)
		paired_bed_tempfile1 = make_tempfile(basename,"paired_1.bed",tempfolder)
		paired_bed_tempfile2 = make_tempfile(basename,"paired_2.bed",tempfolder)
		split_paired(proper_nonproper_labeled_tempfile,paired_bed_tempfile1,paired_bed_tempfile2,debug)
		paired_bed_tempfile1_sorted = paired_bed_tempfile1 + "_sorted"
		paired_bed_tempfile2_sorted = paired_bed_tempfile2 + "_sorted"
		if not debug:
			os.unlink(proper_bam)
			os.unlink(nonproper_bam)
		#reduce reads   #Find reads aligned to same position but different TE_IDs (overlapping flanks) and merge
		if verbosity:
			print("Combining adjacent TEs with same read alignment "+ str(datetime.now())  + "\n",file = sys.stderr)
		paired_reduced_tempfile1 = make_tempfile(basename,"paired_reduced_1", tempfolder)
		paired_reduced_tempfile2 = make_tempfile(basename,"paired_reduced_2", tempfolder)
		sort_coord(paired_bed_tempfile1,paired_bed_tempfile1_sorted,1,2,debug)
		sort_coord(paired_bed_tempfile2,paired_bed_tempfile2_sorted,1,2,debug)

		reduce_reads(paired_bed_tempfile1_sorted, paired_reduced_tempfile1,debug)
		reduce_reads(paired_bed_tempfile2_sorted, paired_reduced_tempfile2,debug)

		#get genomic coordinates and RNA strand for all alignments
		if verbosity:
			print("Getting genomic coordinates of read"+ str(datetime.now())  + "\n",file = sys.stderr)
		paired_coords_tempfile1= make_tempfile(basename,"paired_coords_1", tempfolder)
		paired_coords_tempfile2= make_tempfile(basename,"paired_coords_2", tempfolder)
		get_coords(paired_reduced_tempfile1,1,strandedness,paired_coords_tempfile1,debug)
		get_coords(paired_reduced_tempfile2,2,strandedness,paired_coords_tempfile2,debug)

		paired_labeled_tempfile1 = make_tempfile(basename,"paired_labeled_1", tempfolder)
		paired_labeled_tempfile2 = make_tempfile(basename,"paired_labeled_2", tempfolder)
		label_files(paired_coords_tempfile1,paired_labeled_tempfile1,"R1",debug)
		label_files(paired_coords_tempfile2,paired_labeled_tempfile2,"R2",debug)

		#remove /1 and /2 from read ID column
		paired_fixed_tempfile1 = make_tempfile(basename,"paired_fixed_1", tempfolder)
		paired_fixed_tempfile2 = make_tempfile(basename,"paired_fixed_2", tempfolder)
		fix_paired(paired_labeled_tempfile1,paired_labeled_tempfile2, paired_fixed_tempfile1,paired_fixed_tempfile2,debug)

		#find unique single alignments
		if verbosity:
			print("Identifying and labeling unique and multi reads"+ str(datetime.now())  + "\n",file = sys.stderr)
		first_tempfile1 = make_tempfile(basename,"first_1", tempfolder)
		unique_tempfile1 = make_tempfile(basename,"unique_1", tempfolder)
		multi_tempfile1 = make_tempfile(basename,"multi_1", tempfolder)


		first_tempfile2 = make_tempfile(basename,"first_2", tempfolder)
		unique_tempfile2 = make_tempfile(basename,"unique_2", tempfolder)
		multi_tempfile2 = make_tempfile(basename,"multi_2", tempfolder)

		find_uniq(paired_fixed_tempfile1,first_tempfile1,unique_tempfile1, multi_tempfile1,debug)
		find_uniq(paired_fixed_tempfile2,first_tempfile2, unique_tempfile2, multi_tempfile2,debug)

		#label uniq, multi, or paired

		unique_tempfile1_labeled = make_tempfile(basename,"unique_labeled_1", tempfolder)
		multi_tempfile1_labeled = make_tempfile(basename,"multi_labeled_1", tempfolder)

		unique_tempfile2_labeled = make_tempfile(basename,"unique_labeled_2", tempfolder)
		multi_tempfile2_labeled = make_tempfile(basename,"multi_labeled_2", tempfolder)

		label_files(unique_tempfile1, unique_tempfile1_labeled, "uniq",debug)
		label_files(unique_tempfile2, unique_tempfile2_labeled, "uniq",debug)
		label_files(multi_tempfile1, multi_tempfile1_labeled, "multi",debug)
		label_files(multi_tempfile2, multi_tempfile2_labeled, "multi",debug)

		paired_tempfile1_ulabeled = make_tempfile(basename,"paired_ulabeled_1", tempfolder)
		paired_tempfile2_ulabeled = make_tempfile(basename,"paired_ulabeled_2", tempfolder)
		combine_files(unique_tempfile1_labeled,multi_tempfile1_labeled, paired_tempfile1_ulabeled,debug)
		combine_files(unique_tempfile2_labeled,multi_tempfile2_labeled, paired_tempfile2_ulabeled,debug)

		#combine pairs
		if verbosity:
			print("Matching paired-end mates and merging coordinates"+ str(datetime.now())  + "\n",file = sys.stderr)
		paired_unmatched1= make_tempfile(basename,"paired_unmatched_1", tempfolder)
		paired_unmatched2 = make_tempfile(basename,"paired_unmatched_2", tempfolder)
		paired_matched_tempfile = make_tempfile(basename,"paired_matched", tempfolder)
		match_reads(paired_tempfile1_ulabeled,paired_tempfile2_ulabeled,strandedness,paired_matched_tempfile,paired_unmatched1, paired_unmatched2,debug) #match pairs between paired files

		#sort matched
		matched_tempfile_sorted = make_tempfile(basename,"paired_matched_sorted", tempfolder)
		sort_temp(paired_matched_tempfile,4,matched_tempfile_sorted,debug)

		#combine start and stop of paired reads
		matched_bed = make_tempfile(basename,"matched_bed", tempfolder)
		merge_coords(matched_tempfile_sorted,matched_bed,debug)

		# os.unlink(matched_tempfile_sorted)

		combined_unmatched = make_tempfile(basename,"combined_unmatched", tempfolder)
		combine_files(paired_unmatched1, paired_unmatched2, combined_unmatched,debug)

		#Find single reads that are matched outside of TE but are still proper pair
		if verbosity:
			print("Adding properly paired reads that have mates outside of TE into matched file"+ str(datetime.now())  + "\n",file = sys.stderr)
		proper_single = make_tempfile(basename,"proper_single",tempfolder)
		combined_unmatched2 = make_tempfile(basename,"combined_unmatched2", tempfolder)

		find_proper(combined_unmatched,combined_unmatched2,proper_single,debug) #outputs only nonproper pairs in combined_unmatched2, and single reads that are part of proper pairs in proper_single

		combined_matched = make_tempfile(basename,"combined_matched", tempfolder)
		combine_files(matched_bed,proper_single,combined_matched,debug)
		# os.unlink(proper_single)
		###Remove single alignments of reads that have paired matches using other valid alignments
		if verbosity:
			print("Removing single-end reads that have matching paired-end mates at other alignment locations"+ str(datetime.now())  + "\n",file = sys.stderr)
		only_unmatched = make_tempfile(basename,"only_unmatched", tempfolder)
		remove_repeat_reads(combined_matched,combined_unmatched2,only_unmatched,debug)

		#combine matched and unmatched alignments
		combined_bed = make_tempfile(basename,"combined_bed", tempfolder)
		combine_files(combined_matched,only_unmatched,combined_bed,debug)

		if verbosity:
			print("Identifying and labeling unique and multi fragments"+ str(datetime.now())  + "\n",file = sys.stderr)
		first_tempfile = make_tempfile(basename,"first", tempfolder)
		paired_uniq_tempfile = make_tempfile(basename,"paired_uniq_tempfile",tempfolder)
		unique_bed = make_tempfile(basename,"unique_bed", tempfolder)
		multi_bed_pre = make_tempfile(basename,"multi_bed_pre", tempfolder)
		find_uniq(combined_bed,first_tempfile,unique_bed,multi_bed_pre,debug)

		#find unique pairs in multi_bed
		multi_bed = make_tempfile(basename,"multi_bed", tempfolder)
		if verbosity:
			print("Identifying multi read pairs with one end unique"+ str(datetime.now())  + "\n",file = sys.stderr)
		find_paired_uniq(multi_bed_pre,paired_uniq_tempfile,multi_bed,unique_bed,debug)

		aligned_libsize = getlibsize(logfile, bamfile,multi_bed,unique_bed,paired_end,debug)

	######## COUNT READ(S) #########################
	read_multidict={}  #dictionary to store TE_IDs for each read alignment
	read_locdict = {}  #dictionary to store genomic location of each alignment

	if verbosity:
		print("counting unique alignments "+ str(datetime.now())  + "\n",file = sys.stderr)

	unique_bedfile = open(unique_bed,'r')
	avg_fraglength=uniquecount(unique_bedfile,RepCalc_dict,read_locdict)

	if verbosity:
		print("counting multi alignments "+ str(datetime.now()) + "\n",file = sys.stderr)

	multi_bedfile = open(multi_bed, 'r')
	multicount(multi_bedfile,RepCalc_dict,read_multidict,read_locdict)

	unique_bedfile.close()
	multi_bedfile.close()

	if verbosity:
		print("Adding Tag information to aligned TEs "+ str(datetime.now())  + "\n",file = sys.stderr)

	for TE_ID,RepClass in RepCalc_dict.iteritems():
		RepClass.calcuniqRep()

	if verbosity:
				print("Calculating multialignment assignments "+ str(datetime.now())  + "\n",file = sys.stderr)

	comparedict(read_multidict,RepCalc_dict)
	iteration=0
	if EM == "auto":
			notconverged=True
			prev_read_change=1
			prev_count_change = 0
			max_count_change = 0
			while notconverged:
				iteration +=1
				changed_count = 0
				total_TE =0
				total_TE_0 = 0
				total_TE_1 = 0
				total_TE_10 = 0
				avg_changed_count_pct =0
				max_count_change=0
				total_TE_10_1pct =0
				if verbosity:
					print("Running expectation-maximization calculation for iteration:" + str(iteration) + " " + str(datetime.now())  + "\n",file = sys.stderr)

				for TE_ID,RepClass in RepCalc_dict.iteritems():
					TE_changecount = RepClass.calcmultiRep(iteration)
					max_count_change = max(TE_changecount,max_count_change)
					changed_count +=TE_changecount
					total_TE +=1
					if TE_changecount > 0:
						total_TE_0 +=1
					if TE_changecount >= 1:
						total_TE_1 += 1
					if TE_changecount >= 1 and RepClass.counts_tot >= 10:
						total_TE_10 += 1
					if TE_changecount >= 1 and RepClass.counts_tot >= 10 and (TE_changecount/RepClass.counts_tot) > 0.01:
						total_TE_10_1pct += 1
					avg_changed_count_pct = changed_count/total_TE

				if verbosity:
					print("Average change in TE count:" + str(avg_changed_count_pct) + " " + str(datetime.now())  + "\n",file = sys.stderr)
					print("Max change in TE count:" + str(max_count_change) + " " + str(datetime.now())  + "\n",file = sys.stderr)
					print("Number changed TE:" + str(total_TE_0) + " " + str(datetime.now())  + "\n",file = sys.stderr)
					print("Number TEs changed by at least 1 count:" + str(total_TE_1) + " " + str(datetime.now())  + "\n",file = sys.stderr)
					print("Number TEs changed by at least 1 count with at least 10 counts:" + str(total_TE_10) + " " + str(datetime.now())  + "\n",file = sys.stderr)
					print("Number TEs changed by at least 1 count with at least 10 counts and > 1pct total count:" + str(total_TE_10_1pct) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				new_read_change = estdict(read_multidict,RepCalc_dict)
				if total_TE_10_1pct == 0 and iteration > 1:
					notconverged = False
				else:
					prev_read_change = new_read_change
			if verbosity:
				print("Finished running expectation-maximization calculation after iteration:" + str(iteration) + " " + str(datetime.now())  + "\n",file = sys.stderr)

	elif int(EM) > 0:
		notconverged=True
		prev_read_change=1
		prev_count_change = 0
		max_count_change = 0
		while iteration < int(EM):
			iteration +=1
			changed_count = 0
			total_TE =0
			total_TE_0 = 0
			total_TE_1 = 0
			total_TE_10 = 0
			avg_changed_count_pct =0
			max_count_change=0
			total_TE_10_1pct =0
			if verbosity:
				print("Running expectation-maximization calculation for iteration:" + str(iteration) + " " + str(datetime.now())  + "\n",file = sys.stderr)

			for TE_ID,RepClass in RepCalc_dict.iteritems():
				TE_changecount = RepClass.calcmultiRep(iteration)
				max_count_change = max(TE_changecount,max_count_change)
				changed_count +=TE_changecount
				total_TE +=1
				if TE_changecount > 0:
					total_TE_0 +=1
				if TE_changecount >= 1:
					total_TE_1 += 1
				if TE_changecount >= 1 and RepClass.counts_tot >= 10:
					total_TE_10 += 1
				avg_changed_count_pct = changed_count/total_TE
				if TE_changecount >= 1 and RepClass.counts_tot >= 10 and (TE_changecount/RepClass.counts_tot) > 0.01:
					total_TE_10_1pct += 1
			if verbosity:
				print("Average change in TE count:" + str(avg_changed_count_pct) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				print("Max change in TE count:" + str(max_count_change) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				print("Number changed TE:" + str(total_TE_0) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				print("Number TEs changed by at least 1 count:" + str(total_TE_1) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				print("Number TEs changed by at least 1 count with at least 10 counts:" + str(total_TE_10) + " " + str(datetime.now())  + "\n",file = sys.stderr)
				print("Number TEs changed by at least 1 count with at least 10 counts and > 1pct total count:" + str(total_TE_10_1pct) + " " + str(datetime.now())  + "\n",file = sys.stderr)
			new_read_change = estdict(read_multidict,RepCalc_dict)

	if verbosity:
				print("Writing counts "+ str(datetime.now())  + "\n",file = sys.stderr)

	read_multidict.clear()


	if copies:
		temp_subF = tempfile.NamedTemporaryFile(delete=False, dir = tempfolder, prefix="count" +  ".SFtmp")
		subF_filepath = outfolder + "/" + basename + "_subFcounts.txt"
		subF_file_header = open(subF_filepath + ".header",'w')
#		subF_file_header.writelines("Sample" + "\t" + "aligned_libsize" + "\t" + "Subfamily:Family:Class" + "\t" + "copies" + "\t" + "exp_copies" + "\t" + "uniq_counts" + "\t" + "tot_counts" + "\t" + "avg_conf"  + "\t" + "tot_sense" + "\t" + "tot_antisense" + "\n")
		subF_file_header.writelines("Sample" + "\t" + "aligned_libsize" + "\t" + "Subfamily:Family:Class" + "\t" + "copies"  + "\t" + "fpkm" + "\t" + "uniq_counts" + "\t" + "tot_counts" + "\t" + "tot_reads" + "\t" + "score"  + "\n")

		subF_file_header.close()
		subF_dict = {}

	for TE_ID,RepClass in RepCalc_dict.iteritems(): #for each TE_ID
		RepClass.writeRep(aligned_libsize,counts_temp,basename,strandedness,iteration)
	##########Sort by highest total counts before writing and add to subF dictionary

		if copies:
			subF = get_subF(TE_ID)
			subF_list = split_subF(subF)
			if subF not in subF_dict:
				subF_dict[subF]=subfamily(subF,subF_reads[subF])
				subF_dict[subF].add_TE_count(RepClass,strandedness)
			else:
				subF_dict[subF].add_TE_count(RepClass,strandedness)

	#Close dictionaries from memory
	counts_temp.close()

	sort_counts(counts_temp.name,counts_file_header.name,countsfilepath,5,debug) #sort on 5th column (fpkm)
	read_locdict.clear()
	RepCalc_dict.clear()

	if copies:
		if verbosity:
				print("Writing subfamily counts "+ str(datetime.now())  + "\n",file = sys.stderr)
		with open(copies,'r') as copiesfile: #copiesfile is sorted
			copiesfile.readline() #skip header
			for line in copiesfile:
				line = line.rstrip()
				line_tabs = line.split("\t")
				line_subF = line_tabs[0]
				if line_subF in subF_dict:
					subF_dict[line_subF].add_copy_info(line_tabs)
				#### Write lines ###
				# temp_subF.writelines(basename + "\t" + str(aligned_libsize) + "\t" + line_subF + "\t" + line_copies + "\t" + str(uniq) + "\t" + "{0:.2f}".format(multi) + "\t" + str(subF_conf) + "\t" + "{0:.2f}".format(sense_reads) + "\t"  +  "{0:.2f}".format(antisense_reads) + "\n")
					subF_dict[line_subF].write_subfamily(temp_subF,basename,aligned_libsize,iteration)
				else:
					subF_dict[line_subF]=subfamily(line_subF,0)
					subF_dict[line_subF].add_copy_info(line_tabs)
					subF_dict[line_subF].write_subfamily(temp_subF,basename,aligned_libsize,iteration)

		copiesfile.close()
		temp_subF.close()

		sort_counts(temp_subF.name,subF_file_header.name,subF_filepath,6,debug) #Sort by 7th field (multi)

		if not debug:
			os.unlink(unique_bed)
			os.unlink(multi_bed)

	####### STOP TIMING SCRIPT #######################
	if verbosity:
		print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

		endTime = datetime.now()
		print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
		print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)

###################
if __name__ == "__main__":
	main()

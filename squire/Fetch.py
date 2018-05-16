#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################### MODULES ###################
from __future__ import print_function
import sys
import os
import errno
import argparse #module that passes command-line arguments into script
import subprocess
import glob
import urllib
import urllib2
import tarfile
import gzip
from datetime import datetime
import subprocess as sp
import zipfile
from urllib2 import urlopen
import re
import shutil
import tempfile
import pkg_resources
import warnings
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

def decompress(compressed, decompressed):     #Function for decompressing gzip files
    inF = gzip.open(compressed, 'rb')
    outF = file(decompressed, 'wb')
    for line in inF:
        outF.write(line)
    outF.close()

def unzip(compressed, decompressed): #unzip .zip files
    zip_ref = zipfile.ZipFile(compressed, 'r')
    zip_ref.extractall(decompressed)
    zip_ref.close()

def failed_dl(filepath):          # If the path created by previous steps is empty, break
    with open(filepath) as downloadedfile:
        for i,line in enumerate(downloadedfile):
            if "not found" in line.lower():
                return True
                break
            elif i>10: #if past line 10
                return False
                break

def gtf_to_bed(gtf,bed):    
    #convert gtf to genepred
    genepred=gtf.replace(".gtf",".genepred")
    gtftogenepredcommand_list = ["gtfToGenePred",gtf,genepred] 
    gtftogenepredcommand=" ".join(gtftogenepredcommand_list)
    sp.check_call(["/bin/sh", "-c", gtftogenepredcommand]) 
    #convert genepred to bed
    genepredtobedcommand_list = ["genePredToBed ",genepred,bed] 
    genepredtobedcommand=" ".join(genepredtobedcommand_list)
    sp.check_call(["/bin/sh", "-c", genepredtobedcommand])     

    if not debug:
        os.unlink(genepred)  


def make_tempfile(step, outfolder):
    tmpfile = tempfile.NamedTemporaryFile(delete=False, dir = outfolder, prefix= step +  ".tmp")
    tmpname = tmpfile.name
    tmpfile.close()
    return tmpname

def find_files(folder,pattern, wildpos):
    if wildpos == 1:
        file_list=glob.glob(folder + "/" + "*" + pattern)
    elif wildpos ==2:
        file_list=glob.glob(folder + "/" + pattern + "*")
    if len(file_list) == 0:
        raise Exception("No files found in folder; please give specific " + pattern + " file")
    else:
        return file_list
        
def fix_gtf(infile,outfile):
    outgtf=open(outfile,'w')
    with open(infile,'r') as gtf:
        for line in gtf:
            line = line.rstrip()
            line=line.split()
            attributes=line[8:]
            attribute_col = " ".join(attributes)
            gtf_cols = "\t".join(line[:8])
            outgtf.writelines(gtf_cols + "\t" + attribute_col + "\n")

    outgtf.close()

def sort_coord(infile, outfile,chrcol,startcol):
    chrfieldsort = "-k" + str(chrcol) + "," + str(chrcol)
    startfieldsort = "-k" + str(startcol) + "," + str(startcol) + "n"
    sort_command_list = ["sort",chrfieldsort,startfieldsort, infile, ">", outfile]
    sort_command = " ".join(sort_command_list)
    sp.check_call(["/bin/sh", "-c", sort_command])

def get_basename(filepath):
        filename = os.path.basename(filepath)
        filebase = os.path.splitext(filename)[0]
        return filebase

def gtf_to_bed(gtf,bed):    
    #convert gtf to genepred
    genepred=gtf.replace("gtf","genepred")
    gtftogenepredcommand_list = ["gtfToGenePred",gtf,genepred] 
    gtftogenepredcommand=" ".join(gtftogenepredcommand_list)
    sp.check_call(["/bin/sh", "-c", gtftogenepredcommand]) 
    #convert genepred to bed
    genepredtobedcommand_list = ["genePredToBed ",genepred,bed] 
    genepredtobedcommand=" ".join(genepredtobedcommand_list)
    sp.check_call(["/bin/sh", "-c", genepredtobedcommand])     
    os.unlink(genepred)  


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)  # if no arguments, the below parser statements will be printed
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = "Downloads input files from UCSC")
        parser._optionals.title = "Arguments"
        parser.add_argument("-b","--build", help = "UCSC designation for genome build, eg. 'hg38' (required)", type=str, required = True, metavar = "<build>")
        parser.add_argument("-o","--fetch_folder", help = "Destination folder for downloaded UCSC file(s) (optional; default='squire_fetch')", type=str, default="squire_fetch", metavar = "<folder>")
        parser.add_argument("-f","--fasta", help = "Download chromosome fasta files for build chromosomes (optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-c","--chrom_info", help = "Download chrom_info.txt file with lengths of each chromosome (optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-r","--rmsk", help = "Download Repeatmasker file (optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-g","--gene", help = "Download UCSC gene annotation(optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-x","--index", help = "Create STAR index, WARNING will take a lot of time and memory (optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-p","--pthreads", help = "Launch <int> parallel threads(optional; default='1')", type = int, metavar = "<int>", default=1)
        parser.add_argument("-k","--keep", help = "Keep downloaded compressed files (optional; default=False)", action = "store_true", default=False)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default=False)
        args,extra_args = parser.parse_known_args()

    ###### I/O ############
    build=args.build
    outfolder=args.fetch_folder
    fasta = args.fasta
    chrom_info=args.chrom_info
    rmsk = args.rmsk
    keep = args.keep
    gene = args.gene
    index=args.index
    pthreads=args.pthreads
    verbosity = args.verbosity

    ######### START TIMING SCRIPT ############
    if verbosity:
        startTime = datetime.now()
        print("start time is:" + str(startTime) + '\n', file = sys.stderr)# Prints start time
        print(os.path.basename(__file__) + '\n', file = sys.stderr) #prints script name to std err
        print("Script Arguments" + '\n' + "=================", file = sys.stderr)  #
        args_dict = vars(args)
        for option,arg in args_dict.iteritems():
            print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
        print("\n", file = sys.stderr)


    ######## CHECK IF FOLDER DOESN'T EXIST, OTHERWISE CREATE #########

    make_dir(outfolder)

    ####### DOWNLOAD CHROMOSOME FASTA FILES #########
    if fasta:
        if verbosity:
            print("Downloading Compressed Chromosome files..." + "\n", file = sys.stderr)

        chrom_loc1 = "http://hgdownload.cse.ucsc.edu/goldenPath"  + "/" + build + "/" + "bigZips" + "/" + "chromFa.tar.gz"  # Different file types depending on size/format of chromosome data
        chrom_loc2 = "http://hgdownload.cse.ucsc.edu/goldenPath"  + "/" + build + "/" + "bigZips" + "/" + build + ".chromFa.tar.gz"
        chrom_loc3 = "http://hgdownload.cse.ucsc.edu/goldenPath"  + "/" + build + "/" + "bigZips" + "/" + build + ".fa.gz"
        chrom_loc4 = "http://hgdownload.cse.ucsc.edu/goldenPath"  + "/" + build + "/" + "bigZips" + "/" + "chromFa.zip"
        chrom_basename = outfolder + "/" + build
        chrom_outfolder= chrom_basename + ".chromFa"
        #Download chromosome fasta files

        chrom_name_compressed = chrom_basename + "chromFa.tar.gz"
        urllib.urlretrieve(chrom_loc1, filename=chrom_name_compressed)
        df_fail1=failed_dl(chrom_name_compressed)
        if df_fail1:
            os.unlink(chrom_name_compressed)
            chrom_name_compressed = chrom_basename + "chromFa.tar.gz"
            urllib.urlretrieve(chrom_loc2, filename=chrom_name_compressed)
            df_fail2=failed_dl(chrom_name_compressed)
            if df_fail2:
                os.unlink(chrom_name_compressed)
                chrom_name_compressed = chrom_basename + ".fa.gz"
                urllib.urlretrieve(chrom_loc3, filename=chrom_name_compressed)
                df_fail3=failed_dl(chrom_name_compressed)
                if df_fail3:
                    os.unlink(chrom_name_compressed)
                    chrom_name_compressed = outfolder + "/" + "chromFa.zip"
                    urllib.urlretrieve(chrom_loc4, filename=chrom_name_compressed)
                    df_fail4=failed_dl(chrom_name_compressed)
                    if df_fail4:
                        os.unlink(chrom_name_compressed)
                        raise Exception("Was not able to download chromosome file from UCSC" + "\n", file = sys.stderr)


        if verbosity:
            print("Finished Downloading Compressed Chromosome folder, Decompressing..." + "\n", file = sys.stderr)

        #Unzip
        if "tar.gz" in chrom_name_compressed:
            chrom_name = chrom_outfolder
            with tarfile.TarFile.open(chrom_name_compressed, 'r') as tarredgzippedFile:
                tarredgzippedFile.extractall(path=chrom_name)
        elif "fa.gz" in chrom_name_compressed:
            chrom_name = chrom_outfolder + "/" + build + ".fa"
            decompress(compressed = chrom_name_compressed, decompressed = chrom_name)
        elif "chromFa.zip" in chrom_name_compressed:
            chrom_name = chrom_outfolder
            unzip(chrom_name_compressed,chrom_name)

        if verbosity:
            print("Finished Decompressing Chromosome folder" + "\n", file = sys.stderr)

        #Removes compressed file
        if keep == False:
            if verbosity:
               print("Deleting Compressed Chromosome folder", file=sys.stderr)
            os.remove(chrom_name_compressed)
        #filter for fasta files and filter out unwanted chromosomes
        if os.path.isdir(chrom_name): # if genome is folder and not file
            fasta_folder = chrom_name + "/" + "chroms"
            if not os.path.isdir(fasta_folder): #if chromFa folder does not have "chroms" subdirectory
                fasta_folder = chrom_name #then fasta files are in chromFa folder

            unwanted_folder = chrom_name + "/" + "unwanted"  # create unwanted folder
            make_dir(unwanted_folder)

            file_list=os.listdir(fasta_folder) #list all files in unwanted folder (previously fasta folder)

            unwantedChr = ["hap", "M", "alt"]

            for i in file_list:    # Cleans up unwanted characters from the files before
                i=i.rstrip()
                i_file = fasta_folder + "/" + i
                wanted_file = chrom_outfolder + "/" + i
                unwanted_file = unwanted_folder + "/" + i
                basename = os.path.splitext(i)[0]
                extension = os.path.splitext(i)[1]

                #Filter out folders, non-fasta files, unwanted chromosome fasta files
                if any(x in basename for x in unwantedChr):
                    os.rename(i_file,unwanted_file) # move unwanted chromosome files to chrom.Fa/unwanted folder
                    continue
                if os.path.isdir(i):
                    continue
                if i_file != wanted_file:
                    os.rename(i_file,wanted_file) # move wanted chromosome files to chromFa folder

            if "chroms" in fasta_folder:
                os.rmdir(fasta_folder)

            if verbosity:
                print("Chromosome fasta files are in" + chrom_outfolder + "\n", file = sys.stderr)

    ####### DOWNLOAD CHROM_INFO FILE ##########
    if chrom_info:
        if verbosity:
            print("Downloading Chrom_info file..." + "\n", file = sys.stderr)

        chrom_info_loc = "http://hgdownload.cse.ucsc.edu/goldenPath"   + "/" + build  + "/" + "database" + "/"+ "chromInfo.txt.gz"
        chrom_info_name = outfolder + "/" + build + "_chromInfo.txt"
        chrom_info_name_compressed = chrom_info_name + ".gz"
        #Downloads Chromosome info file
        urllib.urlretrieve(chrom_info_loc, filename=chrom_info_name_compressed)

        if verbosity:
            print("Finished Downloading Chrom_info file, Decompressing..." + "\n", file = sys.stderr)

        #Decompresses chromosome info file
        decompress(compressed = chrom_info_name_compressed, decompressed = chrom_info_name)

        if verbosity:
            print("Finished Decompressing Chrom_info file: "  + "\t" +  chrom_info_name + "\n", file = sys.stderr)
        #Deletes compressed chromosome info file
        if keep == False:
            if verbosity:
               print("Deleting Compressed Chrom_info file" + "\n", file=sys.stderr)
            os.remove(chrom_info_name_compressed)

    ###### DOWNLOAD REPEATMASKER FILE #################
    if rmsk:
        if verbosity:
            print("Downloading Repeatmasker file..." + "\n", file = sys.stderr)
        rmsk_file=outfolder + "/" + build + "_rmsk.txt"
        rmsk_list=set()
        rmsk_loc="http://hgdownload.cse.ucsc.edu/goldenPath"   + "/" + build  + "/" + "database" + "/"

        urlpath = urlopen(rmsk_loc)
        string = urlpath.read().decode('utf-8')

        pattern = re.compile('\\brmsk.txt.gz\\b')
        filelist = pattern.findall(string)
        for filename in filelist:
            rmsk_list.add(filename)

        pattern = re.compile('chr[0-9][0-9]*_rmsk.txt.gz')
        filelist = pattern.findall(string)
        for filename in filelist:
            rmsk_list.add(filename)

        pattern = re.compile('chr[A-Z]_rmsk.txt.gz')
        filelist = pattern.findall(string)
        for filename in filelist:
            rmsk_list.add(filename)

        if len(rmsk_list) > 1:
            if verbosity:
                print("Multiple Repeatmasker files found, Downloading, Decompressing and combining into a single file..." + "\n", file = sys.stderr)
            with open(rmsk_file,'wb') as outfile:
                for filename in rmsk_list:
                    remotefile=urllib.urlretrieve(rmsk_loc + filename, filename=outfolder +"/" + filename)
                    if verbosity:
                           print("Downloading Compressed Repeatmasker file" + " " + filename + "\n", file=sys.stderr)
                    newfilename=filename.replace(".gz","")
                    decompress(compressed=outfolder + "/" + filename, decompressed=outfolder + "/" + newfilename)
                    with open(outfolder + "/" + newfilename, 'rb') as inrmsk:
                        shutil.copyfileobj(inrmsk, outfile)
                        if verbosity:
                           print("Adding to Repeatmasker file" + " " + rmsk_file + "\n", file=sys.stderr)
                     #Deletes decompressed repeatmasker file
                    if keep == False:
                        if verbosity:
                           print("Deleting Compressed Repeatmasker file" + " " + filename + "\n", file=sys.stderr)
                           os.remove(outfolder +"/" + filename)
                        if verbosity:
                           print("Deleting Decompressed Repeatmasker file" + " " + newfilename + "\n", file=sys.stderr)
                           os.remove(outfolder + "/" + newfilename)

        elif len(rmsk_list) == 1:
            rmsk_list=list(rmsk_list)
            filename=rmsk_list[0]
            remotefile=urllib.urlretrieve(rmsk_loc + filename, filename=outfolder +"/" + filename)
            if verbosity:
                print("Finished Downloading Repeatmasker file, Decompressing..." + "\n", file = sys.stderr)
            decompress(compressed=outfolder + "/" + filename, decompressed=rmsk_file)
            if keep == False:
                if verbosity:
                   print("Deleting Compressed Repeatmasker file" + "\n", file=sys.stderr)
                   os.remove(outfolder +"/" + filename)
        elif not rmsk_list:
            raise Exception("Was not able to download rmsk file from UCSC" + "\n", file = sys.stderr)

        if verbosity:
            print("Finished with Repeatmasker download step" + "\n", file = sys.stderr)


    ####### DOWNLOAD GENE ANNOTATIONS ##########
    if gene:
        if verbosity:
            print("Downloading RefGene file..." + "\n", file = sys.stderr)

        refGene_loc = "http://hgdownload.cse.ucsc.edu/goldenPath"   + "/" + build  + "/" + "database" + "/"+ "refGene.txt.gz"
        refGene_name = outfolder + "/" + build + "_refGene.txt"
        refGene_name_compressed = refGene_name + ".gz"
        #Downloads Chromosome info file
        urllib.urlretrieve(refGene_loc, filename=refGene_name_compressed)

        if verbosity:
            print("Finished Downloading refGene file, Decompressing..." + "\n", file = sys.stderr)

        #Decompresses chromosome info file
        decompress(compressed = refGene_name_compressed, decompressed = refGene_name)

        if verbosity:
            print("Finished Decompressing refGene file: "  + "\t" +  refGene_name + "\n", file = sys.stderr)
        #Deletes compressed chromosome info file
        if keep == False:
            if verbosity:
               print("Deleting Compressed refGene file" + "\n", file=sys.stderr)
            os.remove(refGene_name_compressed)


        if verbosity:
            print("Converting RefGene file to GTF ..." + "\n", file = sys.stderr)
        refGene_gtf=outfolder + "/" + build + "_refGene.gtf"
        refGene_temp=make_tempfile("refGene",outfolder)
        refGene_temp2=make_tempfile("refGene2",outfolder)
        refGene_temp3=make_tempfile("refGene3",outfolder)

        genePredToGtf_commandlist = ["cut","-f2-",refGene_name,"|","genePredToGtf","file","stdin",refGene_temp]
        genePredToGtf_command = " ".join(genePredToGtf_commandlist)
        sp.check_call(["/bin/sh", "-c", genePredToGtf_command])

        replace_command_list = ["awk","-v", "OFS='\\t'", """'{ gsub("stdin","hg38_refGene",$2); print $0 }'""", refGene_temp, ">", refGene_temp2]
        replace_command = " ".join(replace_command_list)
        sp.check_call(["/bin/sh","-c",replace_command])

        sort_commandlist = ["sort","-k1,1", "-k4,4n", refGene_temp2, ">", refGene_temp3]
        sort_command = " ".join(sort_commandlist)
        sp.check_call(["/bin/sh", "-c", sort_command])

        fix_gtf(refGene_temp3, refGene_gtf)

        os.remove(refGene_temp)
        os.remove(refGene_temp2)
        os.remove(refGene_temp3)


        if verbosity:
            print("Finished converting RefGene file to GTF ..." + "\n", file = sys.stderr)

        if verbosity:
            print("Finished downloading genePredToBed ..." + "\n", file = sys.stderr)
        if verbosity:
            print("Converting RefGene file to Bed ..." + "\n", file = sys.stderr)
        refGene_Bed=outfolder + "/" + build + "_refGene.bed"
        gtf_to_bed(refGene_gtf,refGene_Bed)

        if verbosity:
            print("Finished converting RefGene file to Bed ..." + "\n", file = sys.stderr)

    ####### CREATE STAR INDEX ##########
    if index:
        chrom_folder = outfolder + "/" + build + ".chromFa"
        if not os.path.isdir(chrom_folder):
            raise Exception(str(chrom_folder) +  "not found" + "\n", file = sys.stderr)
        fasta_list=find_files(chrom_folder,".fa",1)
        genome_filepath = " ".join(fasta_list)
        index_name = outfolder + "/" + build + "_STAR"
        make_dir(index_name)
        STAR_build_commandlist = ["STAR","""--runThreadN""", str(pthreads), """--runMode genomeGenerate""","""--genomeFastaFiles""",genome_filepath,"""--genomeDir""",index_name]
        STAR_build_command = " ".join(STAR_build_commandlist)
        if verbosity:
            print("Building STAR index" + "\n", file = sys.stderr)
            print(STAR_build_command,file=sys.stderr)
        sp.check_call(["/bin/sh", "-c", STAR_build_command])
    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)   # print end time
        print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)   # print total time


#########################################################
if __name__ == "__main__":
    main()

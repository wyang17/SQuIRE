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
    sp.check_call(["/bin/bash", "-c", sort_command])

def get_basename(filepath):
        filename = os.path.basename(filepath)
        filebase = os.path.splitext(filename)[0]
        return filebase


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main(**kwargs):

    ######## ARGUMENTS ###########
    #check if already args is provided, i.e. main() is called from the top level script
    args = kwargs.get('args', None)  # if no arguments, the below parser statements will be printed
    if args is None: ## i.e. standalone script called from command line in normal way
        parser = argparse.ArgumentParser(description = "Installs required software")
        parser._optionals.title = "Arguments"
        parser.add_argument("-b","--build_folder", help = "Destination folder for downloaded UCSC file(s) (optional; default='squire_build')", type=str, default="squire_build", metavar = "<folder>")
        parser.add_argument("-s","--software", help = "Install required SQuIRE software and add to PATH - specify 'all' or provide comma-separated list (no spaces) of: STAR,bedtools,samtools,stringtie (optional; default = False)" , type=str, metavar = "<software>", default=False)
        parser.add_argument("-v","--verbosity", help = "Want messages and runtime printed to stderr (optional; default=False)",  action = "store_true", default=False)
        args,extra_args = parser.parse_known_args()

    ###### I/O ############
    outfolder=args.build_folder
    software = args.software
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


    if software:
        avail_software=["STAR","samtools","bedtools","stringtie"]
        if software == "all":
            software = avail_software
        else:
            software = software.split(",")

        for i in software:
            if i not in avail_software:
                warnings.warn(i + " not in available software")

        software_list = [f for f in pkg_resources.resource_listdir('software', '') if any(i in f for i in software)]

        software_folder = pkg_resources.resource_filename('software', '')

        homepath = os.path.expanduser('~')
        bashrcfile = homepath + "/" + ".bashrc"
        bashprofile_file = homepath + "/" + ".bash_profile"
        outfolder = os.path.abspath(outfolder)
        if verbosity:
            print("Looking for software in" +  str(software_folder) +"..." + "\n", file = sys.stderr)
        bashrc = open(bashrcfile,'a')
        bashprofile = open(bashprofile_file,'a')
        for package in software_list:
            package_path = software_folder + "/" + package
            package_file = get_basename(package).replace(".tar","")
            new_package_path = outfolder + "/" + package_file
            if verbosity:
                print("Decompressing " +  package_file +"..." + "\n", file = sys.stderr)
            with tarfile.TarFile.open(package_path, 'r') as tarredgzippedFile:
                tarredgzippedFile.extractall(new_package_path)
            if verbosity:
                print("Adding permissions for " +  package_file + "and adding to PATH..." + "\n", file = sys.stderr)
            if "bedtools" in package:
                make_folder = new_package_path + "/" "bedtools2"
                make_command_list = ["make"]
                make_command = "".join(make_command_list)
                sp.check_call(["/bin/bash", "-c", make_command],cwd = make_folder)
                new_package_path = make_folder + "/" + "bin"
            elif "samtools" in package:
                make_folder = new_package_path + "/" "samtools-1.1"
                make_command_list = ["make"]
                make_command = "".join(make_command_list)
                sp.check_call(["/bin/bash", "-c", make_command],cwd = make_folder)
                new_package_path = make_folder
            elif "STAR" in package:
                new_package_path = new_package_path + "/STAR-2.5.3a/bin/Linux_x86_64"
            elif "stringtie" in package:
                new_package_path = new_package_path + "/stringtie-1.3.3b.Linux_x86_64"
            newline = """export PATH='""" + new_package_path + """':$PATH"""
            bashrc.writelines(newline + "\n")
            bashprofile.writelines(newline + "\n")

        bashrc.close()
        bashprofile.close()

    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)   # print end time
        print('it took: ' + str(endTime-startTime) + "\n", file = sys.stderr)   # print total time


#########################################################
if __name__ == "__main__":
    main()

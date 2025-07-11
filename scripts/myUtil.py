#!/usr/bin/python
#		Subroutine getallFiles
#		Subroutine packgz
#		Subroutine unpackgz


import os
import sys
#import subprocess
import gzip
import shutil
import random

def packgz(path):
# gzip file from path and return packed file name
    file = path+'.gz'
    with open(path, 'rb') as src, gzip.open(file, 'wb') as dst:
        dst.writelines(src)
    return file



def unpackgz(path):
# gunzip file from path and return unpacked file name
    file = path[:-3]
    with gzip.open(path, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file
    
def command(string):
    os.system(string)
    
    return
    
    
def unlink(path):
    os.unlink(path)
    return
    
def removeExtension(path):
    return os.path.splitext(path)[0] 

def getExtension(path):
    return os.path.splitext(path)[-1]
     
def getExtension2(path):
    name = getFileName(path)
    if os.path.splitext(path)[-1]:
        return "."+name.split(".",1)[-1]
    else:
        return

def removeExtension2(path):
    name = getFileName(path)
    if os.path.splitext(path)[0]:
        return name.split(".",1)[0]
    else:
        return

def getPath(Path):
    return os.path.split(Path)[0]+"/"

def getAllFiles(directory, ending = 0):
#get all files of a directory with all subdirectories and return a list
    list = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            file = os.path.join(path, name)
            if ending == 0:
                list.append(file)
            elif file.endswith(ending):
                list.append(file)
            
    return list

def getFileName(Path):
    file_name = os.path.basename(Path)
    return file_name

def dir_path(string):
#check if path is valid dir
    if os.path.isdir(string):
        if string[-1] == "/":
            return string[:-1]
        else:
            return string
    else:
        sys.exit(f"\nERROR: {string} is not a valid directory")
        #raise Exception(f"\nERROR: {string} is not a valid directory")

def file_path(string):
#check if path is valid file
    if os.path.isfile(string):
        return string
    else:
        sys.exit(f"\nERROR: {string} is not a valid file")
        #raise Exception(f"\nERROR: {string} is not a valid directory")

def compareFileLists(directory,ext1=0,ext2=0):
#return a list of all files with extension 1 which have no equivalent with extension 2
    if ext1 and ext2:
        Files1 = getAllFiles(directory,ext1)
        Files2 = getAllFiles(directory,ext2)
        CompareList1 = removeExtFromList(Files1,ext1)
        CompareList2 = removeExtFromList(Files2,ext2)
        Difference = set(CompareList1).difference(set(CompareList2))
        
        listing = addExtToList(list(Difference),ext1)
                
        return listing

    return

def removeExtFromList(listing,ext):
#removes extension from the right end of each element of a list
    index = 0
    for element in listing:
        element = element[:-len(ext)]
        listing[index] = element
        index += 1
    return listing
    

def addExtToList(listing,ext):
#removes extension from the right end of each element of a list
    index = 0
    for element in listing:
        element += ext
        listing[index] = element
        index += 1
    return listing

def getGenomeID(Path):
    #should return genomeID according to standard genome identifiers from common source
    #29.8.22
    File = getFileName(Path)
    File = removeExtension2(File)
    return File

def getReportName(Input):
    file_name = removeExtension(Input)
    hmm_report = file_name+".HmmReport"
    return hmm_report

def duplicate(input_list):
    #returns a list with the doublicate values
    return list(set([x for x in input_list if input_list.count(x) > 1]))

def taxonomy_lineage(array,trennzeichen):
    #04.11.22
    #concats the taxonomy lineage into one string. 
    #Returned string shall be added at the end of any string
    try:
        div = ''.join(trennzeichen)
        string = div.join(array)
        string = string.replace(" ","-")
        return string
    except:
        return 0

def clean_empty_files(directory):

    # Loop over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file is not a directory
        if not os.path.isdir(os.path.join(directory, filename)):
            # Check if the file is empty
            if os.path.getsize(os.path.join(directory, filename)) == 0:
                # Remove the file
                os.remove(os.path.join(directory, filename))
                #print(filename, 'was removed')

def generate_color(seed_int):
    random.seed(seed_int)
    color = '#{:06x}'.format(random.randint(0, 0xFFFFFF))
    return color








def print_header(string,verbose=0):
    if not verbose:
        print("\n"+string)
        print(len(string)*"-")

def print_start():
    print("HMSSS v 1.0.0")
    print(20*"-")
    #print("2022 by Tomohisa Sebastian Tanabe")
    

    



"""
allFiles = getAllFiles("ttest",".faa.gz")
missingGff = compareFileLists("ttest",".faa.gz",".gff.gz")
#for all missing files try to do prodigal

for File in allFiles:
    unpackedFile = unpackgz(File)
    filename = removeExtension(unpackedFile)
    print(filename)

    filename += ".gff.gz"
    print(filename)
    unpackedGffFile= unpackgz(filename)    
    unlink(unpackedFile)
    unlink(unpackedGffFile)
"""





















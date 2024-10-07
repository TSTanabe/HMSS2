#!/usr/bin/python

import os
import sys
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
    # Check if the file is a .gz file
    if not path.endswith('.gz'):
        return path
    
    # Determine the name of the unpacked file
    file = path[:-3]
    
    # Check if the unpacked file already exists
    if os.path.exists(file):
        return file  # Unpacked file already exists, no need to unpack
    
    # Unpack the .gz file
    with gzip.open(path, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return file
    
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


def print_header(string,verbose=0):
    if not verbose:
        print("\n"+string)
        print(len(string)*"-")
        
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
        
        
        
    
#### Next three routines are meant to work together
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



def getGenomeID(path):
    #return genomeID according to standard genome identifiers from common source. DO NOT USE '.' Where they do not belong in the filename!!!
    basename = os.path.basename(path)
    genomeID = basename.split('.')[0]
    return genomeID

def getReportName(path):
    file_name = os.path.splitext(path)[0]
    hmm_report = file_name+".HmmReport"
    return hmm_report


def taxonomy_lineage(array,trennzeichen):
    #04.11.22
    #concats the taxonomy lineage into one string. 
    #Returned string shall be added at the end of any string
    try:
        div = ''.join(trennzeichen)
        string = div.join(array)
        string = string.replace(" ","-")
        return str(string)
    except:
        return "NoTaxonomy"


def get_executable_dir():
    """
    Get the directory of the current executable or script.
    This works whether the script is compiled or run directly as a Python script.
    """
    if getattr(sys, 'frozen', False):
        # If the program is compiled, sys.frozen is True, and sys.executable gives the path to the executable
        return os.path.dirname(sys.executable)
    else:
        # If running as a script, __file__ gives the path to the script
        return os.path.dirname(os.path.abspath(__file__))

def find_executable(executable):
    """
    Find the any executable.
    First check in the system's PATH, then in the local ./bin directory relative to the executable/script.
    Returns the path to the executable.
    """
    # Check if MAFFT is in the system's PATH
    executable_path = shutil.which(f"{executable}")
    if executable_path:
        return executable_path
    
    # If not found, check in the local ./bin directory relative to the executable or script
    executable_dir = get_executable_dir()
    local_executable_path = os.path.join(executable_dir, "bin", f"{executable}")
    if os.path.isfile(local_executable_path) and os.access(local_executable_path, os.X_OK):
        return local_executable_path
    
    raise FileNotFoundError(f"{executable} executable not found in system PATH or local bin directory.")


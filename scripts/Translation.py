#!/usr/bin/python

import os
import re
import multiprocessing

from . import myUtil


def parallel_translation(directory,cores):
    """
    18.5.24
        Args:  
            directory   fasta file containing directory
            
        Uses prodigal to translate all nucleotide fasta files with fna or fna.gz or .fasta ending to faa
        Warning: if the directory path includes parentheses function prodigal is not working
    """
    zipFnaFiles = myUtil.compareFileLists(directory,".fna.gz",".faa.gz") 
    FnaFiles = myUtil.compareFileLists(directory,".fna",".faa")
    fastaFiles = myUtil.getAllFiles(directory,".fasta")
    NucleotideFastaFiles = zipFnaFiles + FnaFiles + fastaFiles
    print(f"Found {len(NucleotideFastaFiles)} assemblies in nucleotide or ambigous format for prodigal")
    
    manager  = multiprocessing.Manager()
    counter = manager.Value('i',0)
    lock = manager.Lock()
    length = len(NucleotideFastaFiles)
    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length ,counter,lock) for fasta in NucleotideFastaFiles]
        pool.map(translate_fasta, args_list)

    return



def translate_fasta(args):
    fasta,length, counter, lock = args

    #unpack if required
    if os.path.splitext(fasta)[-1] == ".gz":
        fasta = myUtil.unpackgz(fasta)

    #Run prodigal
    output = os.path.splitext(fasta)[0]
    faa = output + ".faa"
    
    prodigal = myUtil.find_executable("prodigal")
    
    string = f"{prodigal} -a {faa} -i {fasta} >/dev/null 2>&1"
    try:
        os.system(string)
    except Exception as e:
        with lock:
            print(f"\tWARNING: Could not translate {fasta} - {e}")
        
    with lock:
        counter.value += 1
        print(f"\rProcessing assembly {counter.value} of {length}", end ='',flush=True)        
    return




############################################################################
############### Parallel Transcription #####################################
############################################################################


def parallel_transcription(directory,cores):
    """
    8.10.22
        Args:  
            directory   fasta file containing directory
            
        Transcribe for all faa files gff3 files
        Secure a packed and unpacked version is present
        Unlink unpacked versions afterwards
        Warning: if the directory path includes parentheses function prodigal is not working
    """
            
    

    gzfaaFiles = myUtil.getAllFiles(directory,".faa.gz")
    gzgffFiles = myUtil.getAllFiles(directory,".gff.gz")
    print(f"Found {len(gzfaaFiles)} zipped faa files")
    print(f"Found {len(gzgffFiles)} zipped gff files")

    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker,gzgffFiles)
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker,gzfaaFiles)   
    
    faaFiles = myUtil.getAllFiles(directory,".faa")
    gffFiles = myUtil.getAllFiles(directory,".gff")   
    print(f"Found {len(gffFiles)} gff files")
    print(f"Found {len(faaFiles)} faa files")    
        	
    FaaFiles = myUtil.compareFileLists(directory,".faa",".gff")
    print(f"Found {len(FaaFiles)} protein fasta files without gff")
    
    manager = multiprocessing.Manager()
    counter = manager.Value('i',0)
    lock = manager.Lock()
    length = len(FaaFiles)
    
    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length, counter, lock) for fasta in FaaFiles]
        pool.map(transcripe_fasta, args_list)
    
    print("\nFinished faa and gff file preparation")
    return


def transcripe_fasta(args):
    fasta, length, counter, lock = args
    gff = ""
            
    if check_prodigal_format(fasta):        
        gff = prodigalFaaToGff(fasta)
    
        with lock:
            counter.value += 1
            print(f"\rProcessing file {counter.value} of {length}", end='', flush=True)
  
    return

def check_prodigal_format(File):
    
    with open(File, "r") as reader:
        for line in reader.readlines():
            if line[0] == ">":
                line = line[1:]
                ar = line.split("#")
                if len(ar) == 5:
                    return 1
                else:
                    return 0
    return Gff

def prodigalFaaToGff(filepath):
#Translates faa files from prodigal to normal gff3 formated files. returns the gff3 file name    
    dir_path = os.path.dirname(filepath)
    # Extract the filename with extensions
    filename_with_ext = os.path.basename(filepath)
    
    if filename_with_ext.endswith('.gz'):
        filename_with_ext = os.path.splitext(filename_with_ext)[0]

    # Remove the .faa extension if present
    if filename_with_ext.endswith('.faa'):
        filename_without_ext = os.path.splitext(filename_with_ext)[0]
    else:
        filename_without_ext = filename_with_ext
    
    Gff = dir_path+"/"+filename_without_ext+'.gff'
    
    writer = open(Gff,"w")
    
    with open(filepath, "r") as reader:
        genomeID = myUtil.getGenomeID(filepath)
        for line in reader.readlines():
            if line[0] == ">":
                try:
                    line = line[1:]
                    ar = line.split("#")
                    #print(ar)
                    contig = re.split("\_{1}\d+\W+$",ar[0])
                    #print(contig)
                    strand = '+' if ar[3] == ' 1 ' else '-'
                    writer.write(contig[0]+"\tprodigal\tcds\t"+ar[1]+"\t"+ar[2]+"\t0.0\t"+strand+"\t0\tID=cds-"+ar[0]+";Genome="+genomeID+"\n")
                except Exception as e:
                    print(f"Error: Missformated header\n {line} - {e}")
    writer.close()
    return Gff


def prepare_gff_faa_packing(directory,cores):

    #Prepare the packing status
    gffFiles = myUtil.getAllFiles(directory,".gff")
    faaFiles = myUtil.getAllFiles(directory,".faa")
    
    if len(gffFiles) > 0 or len(faaFiles) > 0:
        print(f"Found {len(gffFiles)} gff files to pack")
        print(f"Found {len(faaFiles)} faa files to pack")
        with multiprocessing.Pool(processes=cores) as pool:
            pool.map(packer,gffFiles)
        with multiprocessing.Pool(processes=cores) as pool:
            pool.map(packer,faaFiles) 
    return


def packer(file):
    myUtil.packgz(file)
def unpacker(file):
    myUtil.unpackgz(file)    

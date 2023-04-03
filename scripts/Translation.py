#!/usr/bin/python
#		Module PrepareGenomicData
#		Subroutine FaaToGff
#		Subroutine Translation

from datetime import datetime
from . import myUtil
import re




def prodigalFaaToGff(File):
#Translates faa files from prodigal to normal gff3 formated files. returns the gff3 file name    
    #GffFile = re.sub("faa","gff",File)
    Gff = File[:-3] + 'gff'
    writer = open(Gff,"w")

    with open(File, "r") as reader:
        for line in reader.readlines():
            if line[0] == ">":
                try:
                    line = line[1:]
                    ar = line.split("#")
                    #print(ar)
                    contig = re.split("\_{1}\d+\W+$",ar[0])
                    #print(contig)
                    strand = '+' if ar[3] == '1' else '-'
                    writer.write(contig[0]+"\tprodigal\tcds\t"+ar[1]+"\t"+ar[2]+"\t0.0\t"+strand+"\t0\tID=cds-"+ar[0]+";"+ar[3]+"\n")
                except:
                    print(f"Error: Missformated header\n {line}")
    return Gff

def check_prodigal_format(File):
    
    with open(File, "r") as reader:
        for line in reader.readlines():
            if line[0] == ">":
                line = line[1:]
                ar = line.split("#")
                #print(ar)
                #print(len(ar))
                contig = re.split("\_{1}\d+\W+$",ar[0])
                #print(contig)
                #writer.write(contig[0]+"\tprodigal\tcds\t"+ar[1]+"\t"+ar[2]+"\t0.0\t+\t0\tID=cds-"+ar[0]+";"+ar[3]+"\n")
                if len(ar) == 5:
                    return 1
                else:
                    return 0
    return Gff


def translation(directory):
    """
    3.9.22
        Args:  
            directory   fasta file containing directory
            
        Uses prodigal to translate all nucleotide fasta files with fna or fna.gz or .fasta ending to faa
        Warning: if the directory path includes parentheses function prodigal is not working
    """
    zipFnaFiles = myUtil.compareFileLists(directory,".fna.gz",".faa.gz") 
    FnaFiles = myUtil.compareFileLists(directory,".fna",".faa")
    fastaFiles = myUtil.getAllFiles(directory,".fasta")
    NucleotideFastaFiles = zipFnaFiles + FnaFiles + fastaFiles
    print(f"Found {len(NucleotideFastaFiles)} assemblies in nucleotide format for prodigal")
    for index,fasta in enumerate(NucleotideFastaFiles):
        #entpacken falls notwendig
        now = datetime.now()
        print(f"{now} Processing assembly {index+1} of {len(NucleotideFastaFiles)}")
        if myUtil.getExtension(fasta) == ".gz":
            fasta = myUtil.unpackgz(fasta)
        else:
            myUtil.packgz(fasta)    
        #prodigal
        output = myUtil.removeExtension(fasta)
        faa = output + ".faa"
        features = output + ".features"
        string = "prodigal -a "+faa+" -f gff -i "+fasta+" -o "+features+" >/dev/null 2>&1"
        try:
            myUtil.command(string)
        except:
            print(f"\tWARNING: Could not translate {fasta}")
        else:
            gff = prodigalFaaToGff(faa)
            myUtil.packgz(gff)
            myUtil.unlink(gff)

            myUtil.packgz(faa)
            myUtil.unlink(faa)

            myUtil.packgz(features)
            myUtil.unlink(features)
    
        #pack und unlink für mehr speicher

        
        myUtil.unlink(fasta)
        
    return

def transcription(directory):
    """
    8.10.22
        Args:  
            directory   fasta file containing directory
            
        Transcribe for all faa files gff3 files
        Secure a packed and unpacked version is present
        Unlink unpacked versions afterwards
        Warning: if the directory path includes parentheses function prodigal is not working
    """
    
    gffFiles = myUtil.getAllFiles(directory,".gff")
    print(f"Found {len(gffFiles)} gff files to pack")
    for gff in gffFiles:
    	myUtil.packgz(gff)
    	myUtil.unlink(gff)
    faaFiles = myUtil.getAllFiles(directory,".faa")
    print(f"Found {len(faaFiles)} faa files to pack")
    for faa in faaFiles:
    	myUtil.packgz(faa)
    	myUtil.unlink(faa)
    	
    zipFaaFiles = myUtil.compareFileLists(directory,".faa.gz",".gff.gz")
    print(f"Found {len(zipFaaFiles)} protein fasta files without gff")
    
    for index,fasta in enumerate(zipFaaFiles):
        print(f"Trying to generate gff for file {index+1} of {len(zipFaaFiles)}",end="\r")        
        if myUtil.getExtension(fasta) == ".gz":
            fasta = myUtil.unpackgz(fasta)
        else:
            myUtil.packgz(fasta)
            
        if check_prodigal_format(fasta):        
            gff = prodigalFaaToGff(fasta)
            
            myUtil.packgz(gff)
            myUtil.unlink(fasta)
            myUtil.unlink(gff)
            
    
    print("Finished file preparation")
    return






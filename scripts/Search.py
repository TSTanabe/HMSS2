#!/usr/bin/python
#		Subroutine HMMsearch
#		Subroutine BlastSearch
#		Subroutine getGenomeID
from . import myUtil
import re



def makeThresholdDict(File,threshold_type=1):
    """
        30.3.23
        Threshold_type 1 => optimized
        Threshold_type 2 => noise cutoff
        Threshold_type 3 => trusted_cutoff
    """
    Thresholds = {}
    with open(File, "r") as reader:
        for line in reader.readlines():
            #print(line)
            #lines = "key1=value1;key2=value2;key3=value3"
            l = line.split("\t")

            try:
                Thresholds[l[0]] = float(l[threshold_type])
            except:
                print("Error in cutoff line: "+line)
        #for k, v in Thresholds.items():
        #   print(k, v)
    return Thresholds

def HMMsearch(Path,HMMLibrary,cores = 1):
    #Path to faa File, HMMLibrary for the HMMs and cores number of cores used by HMMER3
    #29.8.22

    Output = myUtil.getReportName(Path)
    #print("hmmsearch -E 0.0001 --cpu "+str(cores)+" "+HMMLibrary+" "+Path+">"+Output)
    myUtil.command(f'hmmsearch -E 0.0001 --cpu {str(cores)} {HMMLibrary} {Path}>{Output}')
    return Output

def BlastSearch():
    #Kommandozeile für blastp herausfinden
    return

#def getGenomeID(Path):
    #should return genomeID according to standard genome identifiers from common source
    #29.8.22
    #9.10.22 deprecated, moved to my Util
#    File = myUtil.getFileName(Path)
#    File = myUtil.removeExtension2(File)
#    return File




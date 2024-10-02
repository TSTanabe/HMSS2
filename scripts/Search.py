#!/usr/bin/python
import os



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
    return Thresholds

def get_HMMreport_file(path):
    dir_path = os.path.dirname(path)
    basename = os.path.basename(path)
    basename_without_ext = os.path.splitext(basename)[0]  # Remove the file extension
    output = os.path.join(dir_path, basename_without_ext + '.hmmreport')
    return output

def HMMsearch(path,query_db,options,cores = 1):
    #Path to faa File, HMMLibrary for the HMMs and cores number of cores used by HMMER3
    #will not overwrite old results
    
    score = options.thrs_score
    output = get_HMMreport_file(path)
    if not os.path.isfile(output) or options.clean_reports or options.redo_search:
        os.system(f'hmmsearch -T {score} --domT {score} --cpu {str(cores)} --noali {query_db} {path}>{output}')
    return output


def remove_existing_reports(path):
    output = get_HMMreport_file(path)
    if os.path.isfile(output):
        os.remove(output)



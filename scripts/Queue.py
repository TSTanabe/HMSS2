#!/usr/bin/python

import os
from . import Database
from . import myUtil

     

def queue_files(options):
    """
    Args:
        directory   directory for the result files to be stored
        finished    set with genome identifiers already processed
        options     current options object
    Operation:
        collect all zip protein fasta files and unzipped protein fasta files
        collect all corresponding gff files and check for existence of this file
        if both files present
        get the genome identifiers
    """
    print("\nFilling the queue with faa files --", end="\r")
    genomeID_queue = set()
    faa_files = {}
    gff_files = {}

    pairs = find_faa_gff_pairs(options.fasta_file_directory)
    
    for faa_file,gff_file in pairs:
        genomeID = myUtil.getGenomeID(faa_file)
        genomeID_queue.add(genomeID)
        faa_files[genomeID] = faa_file
        gff_files[genomeID] = gff_file
        
    # compare two sets
    find_missing_genomes(genomeID_queue, options.fasta_file_directory)
    
    options.queued_genomes = genomeID_queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    print("Filling the queue with faa files -- ok")
    print(f"Queued {len(options.queued_genomes)} faa/gff pairs")
    return
    

def compare_with_existing_database(options,genomeIDs):
    
    genomeIDs = Database.fetch_genomeIDs_from_proteins(options.database_directory)
    for genomeID in genomeIDs:
        if genomeID in options.faa_files.keys():
            print(f"\tFound assembly {genomeID} in database leaving out {options.faa_files[genomeID]}")
            del options.faa_files[genomeID]
            del options.gff_files[genomeID]
            options.queued_genomes.remove(genomeID)
    
    print(f"Queued {len(options.queued_genomes)} for processing")
    if len(options.queued_genomes) == 0:
        print("There were 0 genomes queued, as all were already present in the local result database")
    
    return
    

def find_faa_gff_pairs(directory):
    """
    Find pairs of files with the same name but different extensions (.faa/.faa.gz and .gff/.gff.gz)
    in the given directory and its subdirectories.

    Args:
        directory (str): The directory to search for file pairs.

    Returns:
        list: A list of tuples, each containing the paths to a paired .faa and .gff file.
    """
    # Dictionary to store files with the same basename
    files_dict = {}

    # Traverse the directory and its subdirectories
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            
            # Check for .faa or .faa.gz files
            if file.endswith('.faa') or file.endswith('.faa.gz'):
                basename = file.replace('.faa', '').replace('.gz', '')
                if basename not in files_dict:
                    files_dict[basename] = {}
                files_dict[basename]['faa'] = file_path
            
            # Check for .gff or .gff.gz files
            elif file.endswith('.gff') or file.endswith('.gff.gz'):
                basename = file.replace('.gff', '').replace('.gz', '')
                if basename not in files_dict:
                    files_dict[basename] = {}
                files_dict[basename]['gff'] = file_path

    # Find and store pairs of .faa and .gff files
    pairs = []
    for basename, file_paths in files_dict.items():
        if 'faa' in file_paths and 'gff' in file_paths:
            pairs.append((file_paths['faa'], file_paths['gff']))
    return pairs



def find_missing_genomes(genomeIDs, faa_file_directory):
    """Find .faa and .faa.gz files in the directory whose genome IDs are not in the provided list."""
    
    def list_faa_files(directory):
        """List all .faa and .faa.gz files in the directory."""
        return [f for f in os.listdir(directory) if f.endswith('.faa') or f.endswith('.faa.gz')]

    def extract_genomeID_from_faa(filename):
        """Extract genome ID from filename using myUtil.get_genomeID."""
        return myUtil.getGenomeID(filename)

    missing_files = []
    all_faa_files = list_faa_files(faa_file_directory)
    
    for faa_file in all_faa_files:
        genomeID = extract_genomeID_from_faa(faa_file)
        if genomeID not in genomeIDs:
            missing_files.append(faa_file)
    
    return missing_files


def prepare_HMMlib(execute_location):

    # Define the target location for HMMlib
    output_file_path = os.path.join(execute_location, "src", "HMMlib")
    
    # Ensure the directory for HMMlib exists
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    
    if not os.path.isfile(execute_location+"/src/HMMlib"):
        print("Preparing HMMlib from source")
        # Library does not exists, then find all .hmm files in the directory and concat to build the library
        hmm_files = []
        for root, _, files in os.walk(execute_location+"/src/HMMs"):
            for file in files:
                if file.endswith(".hmm"):
                    hmm_files.append(os.path.join(root, file))
        
        # Concatenate all .hmm files using the cat command
        if hmm_files:  # Check if there are .hmm files to concatenate
            cat_command = f"cat {' '.join(hmm_files)} > {output_file_path}"
            os.system(cat_command)


#def create_database_readme(options):
    





#!/usr/bin/python

import os
import shutil
import tempfile
from typing import List, Tuple, Set, Optional, Any

from . import Database
from . import myUtil
     
logger = myUtil.logger

     

def queue_files(options) -> None:
    """
    Fills the options object with genome IDs, .faa and .gff file mappings.
    
    Args:
        options: options object with at least .fasta_file_directory, will be filled with:
            .queued_genomes (set[str])
            .faa_files (dict[str, str])
            .gff_files (dict[str, str])
        
    Operation:
        - Collect all zipped/unzipped protein fasta files and corresponding gff files.
        - Queue only if both files present, by genome identifier.
        
    Output Example:
        options.queued_genomes = {'GCF_000001405.39', ...}
        options.faa_files = {'GCF_000001405.39': '/dir/xxx.faa', ...}
        options.gff_files = {'GCF_000001405.39': '/dir/xxx.gff', ...}
    """
    
    logger.info("Filling the queue with faa files to be processed")
    genomeID_queue = set()
    faa_files = {}
    gff_files = {}

    pairs = find_faa_gff_pairs(options.fasta_file_directory)
    
    for faa_file,gff_file in pairs:
        genomeID = myUtil.getGenomeID(faa_file)
        genomeID_queue.add(genomeID)
        faa_files[genomeID] = faa_file
        gff_files[genomeID] = gff_file
        
    # compare two sets (find missing)
    find_missing_genomes(genomeID_queue, options.fasta_file_directory)
    
    options.queued_genomes = genomeID_queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    
    logger.info(f"Queued {len(options.queued_genomes)} faa/gff pairs")
    
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
    

def find_faa_gff_pairs(directory: str) -> List[Tuple[str, str]]:
    """
    Find pairs of files with the same name but different extensions (.faa/.faa.gz and .gff/.gff.gz)
    in the given directory and its subdirectories.

    Args:
        directory (str): The directory to search for file pairs.

    Returns:
        list of tuple: Each containing the paths to a paired .faa and .gff file.

    Output Example:
        [('/path/xxx.faa', '/path/xxx.gff'), ...]
    """
    
    # Dictionary to store files with the same basename
    files_dict = {}

    # Traverse the directory and its subdirectories
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            
            # Check for .faa or .faa.gz files
            if file.endswith('.faa'): #or file.endswith('.faa.gz'):
                basename = file.replace('.faa', '').replace('.gz', '')
                if basename not in files_dict:
                    files_dict[basename] = {}
                files_dict[basename]['faa'] = file_path
            
            # Check for .gff or .gff.gz files
            elif file.endswith('.gff'): # or file.endswith('.gff.gz'):
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



def find_missing_genomes(genomeIDs: Set[str], faa_file_directory: str) -> List[str]:
    """
    Find .faa files in the directory whose genome IDs are not in the provided list.

    Args:
        genomeIDs (set): Set of genome IDs
        faa_file_directory (str): Directory to search

    Returns:
        List of missing .faa file names (not present in genomeIDs)

    Output Example:
        ['GCF_000001405.39.faa', ...]
    """
        
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


def prepare_HMMlib_shell(
    options: Any,
    execute_location: str,
    allowed_prefixes: Optional[Set[str]] = None
) -> None:
    """
    Uses 'cat' to concatenate all .hmm files from subdirectories with allowed prefixes
    inside src/HMMs into a single HMMlib file.

    Parameters
    ----------
    options : Any
        Argument container with script options (unused here, but kept for compatibility).
    execute_location : str
        Path to the base execution directory.
    allowed_prefixes : set or list of str, optional
        Folder name prefixes to include (e.g., {"grp0", "grp1"}).
        If None or empty, all subdirectories are included.
    """
    output_file_path = os.path.join(execute_location, "src", "HMMlib")
    hmm_base_dir = os.path.join(execute_location, "src", "HMMs")
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    hmm_files = []

    for entry in os.listdir(hmm_base_dir):
        entry_path = os.path.join(hmm_base_dir, entry)
        if os.path.isdir(entry_path):
            if allowed_prefixes:
                if not any(entry.startswith(prefix) for prefix in allowed_prefixes):
                    continue
            # Recursively collect .hmm files
            for root, _, files in os.walk(entry_path):
                for file in files:
                    if file.endswith(".hmm"):
                        hmm_files.append(os.path.join(root, file))

    if hmm_files:
        cat_command = f"cat {' '.join(map(str, hmm_files))} > {output_file_path}"
        logger.debug(f"Running: {cat_command}")
        os.system(cat_command)
        logger.debug(f"All HMM files concatenated into {output_file_path}")
    else:
        logger.error(f"No matching HMM files found in directory {hmm_base_dir}.")


    
def concatenate_files_shell(
    search_directory: str,
    allowed_prefix: str,
    allowed_suffix: str,
    output_file_path: str
) -> None:
    """
    Rekursiv alle Files suchen, die mit allowed_prefix beginnen und allowed_suffix enden.
    Diese zusammenführen und als output_file_path speichern.
    """
    matched_files = []
    for root, _, files in os.walk(search_directory):
        for fname in files:
            if fname.startswith(allowed_prefix) and fname.endswith(allowed_suffix):
                matched_files.append(os.path.join(root, fname))

    if matched_files:
        os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
        cat_command = 'cat ' + ' '.join(f'"{f}"' for f in matched_files) + f' > "{output_file_path}"'
        #logger.debug(f"Running: {cat_command}")
        os.system(cat_command)
        logger.debug(f"Concatenated {len(matched_files)} files into {output_file_path}")
    else:
        logger.error(f"No matching files found for concatenation in {search_directory}.")


def format_pattern_files_inplace(filename: str, prefix: str, suffix: str):
    tmpfile = tempfile.NamedTemporaryFile("w", delete=False)
    with open(filename, "r") as fin, tmpfile:
        for i, line in enumerate(fin, 1):
            new_line = f"{prefix}{i}{suffix} {line.rstrip()}"
            new_line = new_line.replace(" ", "\t")
            tmpfile.write(new_line + "\n")
    shutil.move(tmpfile.name, filename)




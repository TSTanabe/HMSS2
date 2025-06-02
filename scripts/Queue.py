#!/usr/bin/python

import sys
import os
from . import Database
from . import myUtil

     

def queue_files(options):
    """
    Populate the queue with paired .faa and .gff files, extracting genome IDs.

    Args:
        options: Object containing script options. Expected attributes:
            - fasta_file_directory (str): Directory to search for .faa/.gff files.
    Operation:
        1. Find matching .faa/.gff pairs.
        2. Extract genome IDs from the .faa filenames.
        3. Track missing genomes (if a pair is incomplete).
        4. Store results in options.queued_genomes, options.faa_files, options.gff_files.
    """
    print("\n[INFO] Filling the queue with faa files --", end="\r")

    genome_id_queue = set()
    faa_files_map = {}
    gff_files_map = {}

    # find_faa_gff_pairs returns an iterable of (faa_path, gff_path) tuples
    pairs = find_faa_gff_pairs(options.fasta_file_directory)

    for faa_path, gff_path in pairs:
        genome_id = myUtil.getGenomeID(faa_path)
        genome_id_queue.add(genome_id)
        faa_files_map[genome_id] = faa_path
        gff_files_map[genome_id] = gff_path

    # Compare sets to report any missing genomes (i.e., unmatched files)
    find_missing_genomes(genome_id_queue, options.fasta_file_directory)

    options.queued_genomes = genome_id_queue
    options.faa_files = faa_files_map
    options.gff_files = gff_files_map

    print("[INFO] Filling the queue with faa files -- ok")
    print(f"[INFO] Queued {len(options.queued_genomes)} faa/gff pairs")
    

def compare_with_existing_database(options, genome_ids):
    """
    Remove already-processed genomes from the queue based on the database contents.

    Args:
        options: Object containing script options. Expected attributes:
            - database_directory (str): Path to the existing database.
            - faa_files (dict): Mapping genome_id → path to .faa file.
            - gff_files (dict): Mapping genome_id → path to .gff file.
            - queued_genomes (set): Set of genome IDs currently queued.
        genome_ids: (Ignored) Placeholder to match signature; actual IDs are fetched from the database.

    Operation:
        1. Fetch all genome IDs already present in the database (based on protein entries).
        2. For each genome ID found in both the database and the current queue:
           a. Print a message indicating it will be skipped.
           b. Remove it from options.faa_files, options.gff_files, and options.queued_genomes.
        3. Print the final count of genomes still queued.
        4. If no genomes remain, print an informational message.
    """
    # Fetch genome IDs present in the database (based on protein entries)
    existing_ids = Database.fetch_genomeIDs_from_proteins(options.database_directory)

    for genome_id in existing_ids:
        if genome_id in options.faa_files:
            print(f"\tFound assembly {genome_id} in database; skipping {options.faa_files[genome_id]}")
            # Remove from FAA/GFF maps and queued set
            options.faa_files.pop(genome_id, None)
            options.gff_files.pop(genome_id, None)
            options.queued_genomes.discard(genome_id)

    remaining = len(options.queued_genomes)
    print(f"Queued {remaining} genome(s) for processing")

    if remaining == 0:
        print("No genomes queued; all were already present in the local result database")

    

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


def prepare_HMMlib(options, execute_location, allowed_prefixes=None):
    """
    Prepares the HMM library by concatenating .hmm files from subdirectories
    inside src/HMMs, optionally filtering by given folder-name prefixes.

    Args:
        options: Argument container with script options (unused here but kept for signature).
        execute_location (str): Path to the base execution directory.
        allowed_prefixes (set or list, optional): If provided, only subdirectories
            whose names start with one of these prefixes will be included. If None or
            empty, all subdirectories are used.
    """
    # Define paths
    hmm_base_dir = os.path.join(execute_location, "src", "HMMs")
    output_file_path = os.path.join(execute_location, "src", "HMMlib")

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    print("[INFO] Preparing HMMlib from source")

    # If allowed_prefixes is provided, ensure it's a set for fast lookup
    if allowed_prefixes:
        allowed_prefixes = set(allowed_prefixes)
    else:
        allowed_prefixes = None

    # Collect all .hmm files
    hmm_files = []

    if not os.path.isdir(hmm_base_dir):
        print(f"[ERROR] HMM base directory does not exist: {hmm_base_dir}")
        return

    for entry in os.listdir(hmm_base_dir):
        entry_path = os.path.join(hmm_base_dir, entry)

        # Only consider directories
        if not os.path.isdir(entry_path):
            continue

        # If prefixes are specified, skip directories that don't start with any prefix
        if allowed_prefixes:
            if not any(entry.startswith(prefix) for prefix in allowed_prefixes):
                continue

        # Recursively find all .hmm files under this subdirectory
        for root, _, files in os.walk(entry_path):
            for fname in files:
                if fname.endswith(".hmm"):
                    hmm_files.append(os.path.join(root, fname))

    if not hmm_files:
        print("[INFO] No HMM files found to concatenate.")
        return

    print(f"[INFO] Found {len(hmm_files)} HMM file(s) to concatenate.")

    # Concatenate them into the output file
    try:
        with open(output_file_path, "wb") as outfile:
            for file_path in hmm_files:
                with open(file_path, "rb") as infile:
                    outfile.write(infile.read())
        print(f"[INFO] Successfully wrote concatenated HMMs to {output_file_path}")
    except OSError as e:
        print(f"[ERROR] Failed to write HMMlib: {e}")

    





#!/usr/bin/python
import os
import traceback
import subprocess

import multiprocessing
from multiprocessing import Manager, Pool

from . import myUtil
from . import Search
from . import Csb_finder
from . import Database
from . import ParseReports
from . import Output


######## MAIN ROUTINE ##########
#If only 1 File is given as concat of genomes
def initial_bulk_search(options):
    #2 challenges: The resources for so many parallel searches are enormous, definitly better to use them on an HPC separatly (some HMMs need over 60 Gb RAM alone overload a single node)
    #gets a bulk faa and a bulk gff
    
    #create two folders inside the search directory
    
    #folder for hmmreports that is stored in globreports folder
    
    #folder for the single gff and faa files if not provided in the first place #TODO insert this into options
    
    #fork of the search HMMsearch with the single HMMs provided in the HMMlib
    
    return
    
    


def initial_bulk_reports(options):
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type)
    
    hmmreport_tsv = get_hmmreport_tsv(options, score_threshold_diction)
    
    genomeIDs_set = collect_genomeIDs(hmmreport_tsv) #returns a set of all genomeIDs
    genomeIDs_set = remove_absent_faa_genomeIDs(genomeIDs_set,options) # remove the genomeIDs that are not present in faa files
    print(f"Found hits in {len(genomeIDs_set)} of the given .faa files")
    Database.insert_database_genomeIDs(options.database_directory, genomeIDs_set) # Insert genomeIDs into database
    
    genomeID_batches = split_genomeIDs_into_batches(list(genomeIDs_set), options.cores-1) # sets of genomeIDs
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)

    process_instances = []
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    with Pool(processes = options.cores, initializer=init_worker, 
              initargs=(hmmreport_tsv, score_threshold_diction, csb_patterns_diction, csb_pattern_names)) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options, counter))
        
        pool.map(process_parallel_bulk_parse_batch, [(data_queue, batch, options, counter) for batch in genomeID_batches])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
    print("Finished parsing reports")
    return    

# Initializer function to set up global variables for each worker
def init_worker(hmmreport_tsv, score_threshold_diction, csb_patterns_diction, csb_pattern_names):
    global global_hmmreport_tsv
    global global_score_threshold_diction
    global global_csb_patterns_diction
    global global_csb_pattern_names
    
    # Assign the arguments to the global variables in worker processes
    global_hmmreport_tsv = hmmreport_tsv
    global_score_threshold_diction = score_threshold_diction
    global_csb_patterns_diction = csb_patterns_diction
    global_csb_pattern_names = csb_pattern_names
    
##################### Routines for tsv file preparation ###################
def get_hmmreport_tsv(options, score_threshold_diction):
    #Gets a directory with separate hmmreport files. These are then prepared and concatenated by HMSS2 to a single light weight tsv file
    
    

    #Prepare the hmmreport as light weight tsv

    hmmreport_tsv = None
    print(f"Processing concatenated hmmreport {options.glob_report}") #print control sequence output    
    if os.path.isdir(options.glob_report): #unconcatened files are provided
        output_filepath = os.path.join(options.glob_report, 'concat.hmmreport')
        if os.path.isfile(output_filepath): #if concat.hmmreport exists use this one to spare time
            print(f"Found existing report {output_filepath}. This will be processed")
            hmmreport_tsv = output_filepath    
        else:
            hmmreport_tsv = prepare_hmmreports_from_dir(options.glob_report,score_threshold_diction,options.cores) #Concat to concat.hmmreport if needed
    else:
        print("Error: Invalid directory for .hmmreport files. Separate reports must have the .hmmreport ending")
    
    return hmmreport_tsv


def prepare_hmmreports_from_dir(hmmreports_dir,Thresholds,cores):
    hmmreport_files = myUtil.getAllFiles(hmmreports_dir,"hmmreport")
    args_list = []
    
    for file_path in hmmreport_files:
        query = get_query_from_line(file_path)
        if not query:
            print(f"Warning: No query found for {file_path}. Skipping the loop")
            continue
        print(f"Adding query {query} into the queue")
        threshold = 10
        if query in Thresholds.keys():
                threshold = Thresholds[query]
        args = (file_path, query, threshold)
        args_list.append(args)
    with Pool(processes=cores) as pool:
        pool.map(write_hmmreport_to_tsv, args_list)
        
    print("Now concatenating the filtered hmmreports")
    
    # Define the output file path
    output_filepath = os.path.join(hmmreports_dir, 'concat.hmmreport')
    os.system(f"cat {hmmreports_dir}/*.filtered_report > {output_filepath}")


    return output_filepath


def extract_query_name(file_path):
#given a hmmreport this returns the name of the Query using grep function
    try:
        # Use grep to find the line containing "Query:"
        command = f"grep 'Query:' {file_path}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error reading file: {result.stderr}")
            return None
        
        # Extract the line containing "Query:"
        line = result.stdout.strip()
        
        # Split the line and extract the query name
        parts = line.split()
        if len(parts) >= 2:
            query_name = parts[1]
            return query_name
        else:
            print("Error: Query name not found in the line.")
            return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def get_query_from_line(file_path, line_number=11):
    try:
        with open(file_path, 'r') as file:
            for current_line_number, line in enumerate(file, start=1):
                if current_line_number == line_number:
                    if line.startswith("Query:"):
                        string = line.strip()
                        parts = string.split()
                        if len(parts) >= 2:
                            query_name = parts[1]
                            return query_name
                        else:
                            return None
                    else:
                        print(f"No 'Query:' found at line {line_number}")
                        return None
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
        
        
        

def write_hmmreport_to_tsv(args):
    file_path,query,threshold = args

    output_filepath = file_path.rsplit('.', 1)[0] + ".filtered_report"
    if os.path.isfile(output_filepath):
        return
    print(f"Processing {query} with threshold {threshold}")    
    hit_dict = parse_hmmreport_hit_table(file_path,query,threshold) #parse the first part of the hmmreport
    hit_dict = parse_hmmreport_domain_table(file_path, hit_dict)

    tsvfile = open(output_filepath, 'w')
    for line in hit_dict.values(): 
        tsvfile.write(line + '\n')
    tsvfile.close()
    hit_dict.clear()

def parse_hmmreport_hit_table(file_path,query,threshold):
    # Initialize an empty dictionary to store hits
    hit_dict = dict()
    file = None
    try:
        file = open(file_path, 'r')
        # Read lines from the file
        lines = file.readlines()
        # Ensure we start from the start_line (1-indexed)
        for i in range(15, len(lines)):
            line = lines[i].strip()
            if line == "":
                break
            try:
                # Split the line at whitespace
                array = line.split()
                # Convert the second element to float and compare with threshold
                if float(array[1]) >= float(threshold):
                    # Print the array for debugging
                    #print(line.split()) # array contains: evalue, score, bias, and sequence identifier
                    # Populate hit_dict with the required information
                    hit_dict[array[8]] = f"{array[8]}\t{query}\t{array[0]}\t{array[1]}\t{array[2]}" # hitID, query, evalue, score, bias
                else:
                    # Explicitly close the file before returning
                    file.close()
                    return hit_dict
            except Exception as e:
                # Skip lines that cause an exception and continue with the next line
                print(f"Error: Skipping line due to error: {e}")
                continue
    except Exception as e:
        print(f"Error: reading file {file_path}: {e}")
    finally:
        if file is not None:
            file.close()
    
    return hit_dict

def parse_hmmreport_domain_table(file_path, gene_dict):
    start_string = "Domain annotation for each sequence:"
    start_reading = False
    current_gene = None

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()

                # Start reading after the start_string is found
                if start_reading:
                    if line.startswith(">> "):
                        current_gene = line[3:]
                        if current_gene not in gene_dict:
                            # Stop processing if the gene is not in the dictionary
                            return gene_dict

                    # Skip lines that are uninteresting
                    if line.startswith("#") or line.startswith("---") or line == "":
                        continue

                    # Extract alifrom and ali to values from the relevant line

                    if current_gene and len(line.split()) > 8:
                        columns = line.split()
                        alifrom = columns[6]
                        alito = columns[7]
                        try:
                            # Check if alifrom and alito can be converted to float
                            alifrom_float = float(alifrom)
                            alito_float = float(alito)
                            # Update the dictionary if conversion is successful
                            if current_gene in gene_dict:
                                gene_dict[current_gene] += f"\t{alifrom}\t{alito}"
                                current_gene = None  # Reset current gene after processing
                        except ValueError:
                            # If conversion fails, remove current_gene from the dictionary
                            if current_gene in gene_dict:
                                del gene_dict[current_gene]
                            print(f"Error: Removed {current_gene} due to invalid domain coordinates")
                            current_gene = None  # Reset current gene after processing

                elif start_string in line:
                    start_reading = True

    except Exception as e:
        print(f"Error: {e}")

    return gene_dict


############################ Collect genomeIDs ###################################

def collect_genomeIDs(report_path, div='___'):
    """
    Collect genome IDs from a TSV file where the hit identifier is in the first column.

    Parameters:
    report_path (str): Path to the TSV file.
    div (str): Divider used when the proteinID includes the genomeID.

    Returns:
    set: A set of unique genome IDs.
    """
    genomeIDs = set()
    
    with open(report_path, 'r') as tsvfile:
        for line in tsvfile:
            line = line.strip()  # Remove leading/trailing whitespace
            if line:  # Ensure the line is not empty
                columns = line.split('\t')
                key = columns[0].split(div)[0]
                genomeIDs.add(key)
    
    return genomeIDs

def remove_absent_faa_genomeIDs(genomeID_set, options):
    
    genomeID_set_copy = genomeID_set.copy()
    
    for genomeID in genomeID_set_copy:
        if genomeID not in options.faa_files:
            genomeID_set.remove(genomeID)  # Entferne genomeID aus dem Set
    
    return genomeID_set
################# Prepare the genomeID sets for parallel processes ################
def split_genomeIDs_into_batches(genomeIDs_list, num_batches):
    """
    Splits a list of genomeIDs into nearly equal-sized batches without using math.ceil().
    
    Args:
        genomeIDs_list (list): The list of genomeIDs to be split.
        num_batches (int): The number of batches to create.
    
    Returns:
        list: A list containing `num_batches` sublists, each containing a portion of the genomeIDs.
    """
    # Calculate the approximate size of each batch
    batch_size = len(genomeIDs_list) // num_batches
    remainder = len(genomeIDs_list) % num_batches

    # Create the batches, distributing the remainder across the first few batches
    batches = []
    start = 0
    for i in range(num_batches):
        # Distribute the remainder across the first 'remainder' batches
        end = start + batch_size + (1 if i < remainder else 0)
        batches.append(genomeIDs_list[start:end])
        start = end

    return batches
    
def process_parallel_bulk_parse_batch(args):
    """
    This function processes a batch of genomeIDs in parallel.
    
    Args:
        args: Tuple containing data_queue, genomeID_batch, and other required parameters.
    """
    data_queue, genomeID_batch, options, counter = args
    
    #Get global variables for read only
    hmmreport_tsv = global_hmmreport_tsv
    score_threshold_diction = global_score_threshold_diction
    csb_patterns_diction = global_csb_patterns_diction
    csb_pattern_names = global_csb_pattern_names
    
    
    # Process each genomeID in the batch
    for genomeID in genomeID_batch:
        process_parallel_bulk_parse((data_queue, genomeID, options, hmmreport_tsv, score_threshold_diction, csb_patterns_diction, csb_pattern_names, counter))

    
def process_parallel_bulk_parse(args_tuple):
    data_queue,genomeID,options,hmmreport_tsv,score_threshold_diction, csb_patterns_diction,csb_pattern_names,counter = args_tuple
    
    #counter.value +=1
    #print(f"Processed {counter.value} genomes by worker", end="\r")#
    protein_dict = dict()
    cluster_dict = dict()
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        
        #Parse the hits
        protein_dict = parse_bulk_HMMreport_genomize(genomeID,hmmreport_tsv,score_threshold_diction,options.thrs_score)
        #Complete the hit information
        ParseReports.parseGFFfile(gff_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,protein_dict)
        
        #Find the syntenic regions and insert to database
        #If search is redone do not calculate csb, but do it separately with redo csb after all
        if not options.redo_search:
            cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
            Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,options.min_completeness)


            if options.synthenic_block_support_detection:
                missing_protein_types, missing_proteins_list_dict = Csb_finder.find_csb_pattern_difference(csb_patterns_diction,csb_pattern_names,cluster_dict,3)
                candidate_proteins_dict = ParseReports.parseHMMreport_below_cutoff_hits(missing_protein_types,hmm_report,score_threshold_diction,options.cut_score)
                Csb_finder.synteny_completion(gff_file,protein_dict,cluster_dict,candidate_proteins_dict,missing_proteins_list_dict,options.nucleotide_range)
        data_queue.put((protein_dict,cluster_dict))
    except Exception as e:
        error_message = f"\nError: occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
        return
    
    return


def process_writer(queue, options, counter):
    # This routine handles the output of the search and writes it into the database
    # It gets input from multiple workers as the database connection to sqlite is unique
    
    protein_batch = {}
    cluster_batch = {}
    batch_size = options.glob_chunks
    batch_counter = 0

    while True:
        tup = queue.get()
        if tup is None:
            # Process any remaining data in the batch
            if protein_batch and cluster_batch:
                submit_batches(protein_batch, cluster_batch, options)
            break
        
        else:
            counter.value += 1
            print(f"Processed {counter.value} genomes ", end="\r")#
        
        protein_dict, cluster_dict = tup

        # Concatenate the data
        protein_batch.update(protein_dict)
        cluster_batch.update(cluster_dict)
        batch_counter += 1
        
        # Print text reports if desired
        if options.individual_reports:
            if protein_dict:  # Check if protein_dict is not empty
                first_protein_key = next(iter(protein_dict))  # Get the first key
                genomeID = protein_dict[first_protein_key].genomeID
                filepath = os.path.join(options.fasta_initial_hit_directory, str(genomeID))
                Output.output_genome_report(filepath, protein_dict, cluster_dict)
        
        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            submit_batches(protein_batch, cluster_batch, options)
            protein_batch = {}
            cluster_batch = {}
            batch_counter = 0

    return

def submit_batches(protein_batch, cluster_batch, options):
    # Submit the batches to the database and write to the file

    # Insert into the database
    Database.insert_database_proteins(options.database_directory, protein_batch)
    Database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')




def parse_bulk_HMMreport_genomize(genomeID,Filepath,Thresholds,cut_score=10):
    """
    1.9.22 
    Required input are a path to a Hmmreport File from HMMER3 and a thresholds dictionary with threshold scores for each HMM
    Returns a list of protein objects for further utilization
    Args
        -Inputfile Report from HMMER3
        -Dictionary Thresholds for HMMs
        
    Return
        -list of Protein objects
    """

    protein_dict = {}

    result = subprocess.run(['grep', genomeID, Filepath], stdout=subprocess.PIPE, text=True)
    
    lines = result.stdout.splitlines()  # Split output into lines

    for line in lines:
        columns = line.split('\t')  # Assuming columns are space-separated
        if columns:
            try:
                key = columns[0]  # Last column as key
                #{hit_id}\t{query_id}\t{e_value}\t{score}\t{bias}\t{hsp_start}\t{hsp_end}\t{description}
                hit_proteinID = columns[0]
                query = columns[1]
                hit_bitscore = int(float(columns[3]))
                hsp_start = int(float(columns[5]))
                hsp_end = int(float(columns[6]))
                
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = ParseReports.Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore,genomeID)
            except Exception as e:
                error_message = f"\nError occurred: {str(e)}"
                traceback_details = traceback.format_exc()
                print(f"\tWARNING: Skipped {Filepath} due to an error - {error_message}")
                print(f"\tTraceback details:\n{traceback_details}")
                continue
                
    return protein_dict

def generate_genomeIDs(protein_dict,div='___'):
    #Div is the divider when proteinID includes the genomeID
    for protein in protein_dict.values():
        proteinID = protein.proteinID 
        genomeID = proteinID.split(div,1)[0]
        protein.genomeID = genomeID
    return




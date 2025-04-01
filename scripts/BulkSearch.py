#!/usr/bin/python
import os
import traceback
import subprocess
import pprint

import multiprocessing
from multiprocessing import Manager, Pool
from functools import partial

from . import myUtil
from . import Search
from . import Csb_finder
from . import Database
from . import ParseReports
from . import Output


######## MAIN ROUTINE ##########
#If only 1 File is given as concat of genomes
global_hmmreport_tsv = None
global_score_threshold_dict = None
global_csb_patterns_dict = None
global_csb_pattern_names = None

def initial_bulk_reports(options):
    global global_hmmreport_tsv
    global global_score_threshold_dict
    global global_csb_patterns_dict
    global global_csb_pattern_names
    global global_deconcat_domains_report_dict

    # Make the noise cutoff files for the synteny completion
    score_threshold_dict = Search.makeThresholdDict(options.score_threshold_file, 3, options.thrs_score)
    hmmreport_tsv = get_hmmreport_tsv(options, score_threshold_dict, "concat.noise_report")
    
    # Remove the hits above trusted cutoff first, because they are always found by the HMM
    score_threshold_dict = Search.makeThresholdDict(options.score_threshold_file, 2, options.thrs_score)
    hmmreport_tsv = remove_trusted_hits(score_threshold_dict,hmmreport_tsv)
    separate_domain_report_dict = deconcat_hmmreport_to_per_domain_noise_cutoff_report(hmmreport_tsv)
    
    #remove the .hmmreport.filtered_report files from options.glob_report
    [os.remove(os.path.join(options.glob_report, f)) for f in os.listdir(options.glob_report) if f.endswith('.hmmreport.filtered_report')]
    
    #Make the concat file for the grep of hits
    score_threshold_dict = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type, options.thrs_score)
    print(f"Cutoffs:\n{score_threshold_dict}")
    hmmreport_tsv = get_hmmreport_tsv(options, score_threshold_dict, "concat.report")
    
    
    global_hmmreport_tsv = hmmreport_tsv
    global_score_threshold_dict = score_threshold_dict
    global_csb_patterns_dict, global_csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    global_deconcat_domains_report_dict = separate_domain_report_dict
    
    ### Collect and insert genome identifiers
    genomeIDs_set = collect_genomeIDs(hmmreport_tsv, options) #returns a set of all genomeIDs
    print(f"Found hits in {len(genomeIDs_set)} of the given .faa files")
    Database.insert_database_genomeIDs(options.database_directory, genomeIDs_set) # Insert genomeIDs into database
    
    genomeID_batches = split_into_batches(list(genomeIDs_set), options.cores-1) # sets of genomeIDs
    
    with Manager() as manager:
        data_queue = manager.Queue()
        counter = manager.Value('i', 0)

        with Pool(processes = options.cores) as pool:
            p_writer = pool.apply_async(process_writer, (data_queue, options, counter))
            batch_args = [(data_queue, batch, options, counter) for batch in genomeID_batches]
            pool.map(process_batch, batch_args)

            # Signal the end of data to the process_writer
            for _ in range(options.cores):
                print("Finished parsing reports")
                data_queue.put(None)
        
            p_writer.get()
    print("Finished database insertion")
    return    

    
##################### Routines for tsv file preparation ###################
def get_hmmreport_tsv(options, score_threshold_dict, output_name="concat.report"):
    if os.path.isdir(options.glob_report):
        output_filepath = os.path.join(options.glob_report, output_name)
        if os.path.isfile(output_filepath):
            print(f"Found existing report {output_filepath}. Using this to spare time.")
            return output_filepath
        else:
            return prepare_hmmreports(options.glob_report, score_threshold_dict, options.cores, output_name)
    else:
        raise ValueError("Invalid directory for .hmmreport files.")

def prepare_hmmreports(directory, thresholds, cores, output_name):
    report_files = myUtil.getAllFiles(directory, "hmmreport")
    args_list = [(file, get_query_from_line(file), thresholds.get(get_query_from_line(file), 10)) for file in report_files]

    with Pool(processes=cores) as pool:
        pool.map(write_filtered_report, args_list)

    output_filepath = os.path.join(directory, output_name)
    subprocess.run(f"cat {directory}/*.filtered_report > {output_filepath}", shell=True, check=True)

    return output_filepath

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
        
        
        

def write_filtered_report(args):
    print("Filtering report")
    file_path, query, threshold = args
    output_filepath = f"{file_path}.filtered_report"
    if os.path.isfile(output_filepath):
        return

    print(f"Processing {query} with threshold {threshold}")
    hit_dict = parse_hit_table(file_path, query, threshold)
    hit_dict = parse_domain_table(file_path, hit_dict)

    with open(output_filepath, 'w') as tsvfile:
        for line in hit_dict.values():
            tsvfile.write(line + '\n')


def parse_hit_table(file_path, query, threshold):
    """
    Parses a hit table file and returns a dictionary of hits.

    Parameters:
        file_path (str): Path to the file containing hit data.
        query (str): Query identifier to include in the output.
        threshold (float): Minimum score threshold for hits to be included.

    Returns:
        dict: A dictionary containing hit data.
    """
    hit_dict = {}
    file = None

    try:
        file = open(file_path, 'r')  # Explicitly open the file
        lines = file.readlines()  # Read all lines at once
        
        # Process lines starting from line 16 (0-indexed 15)
        for line in lines[15:]:
            line = line.strip()
            if not line:  # Break on an empty line
                break

            try:
                # Split the line and parse the necessary fields
                array = line.split()
                score = float(array[1])

                if score >= threshold:
                    hit_id = array[8]
                    hit_dict[hit_id] = f"{hit_id}\t{query}\t{array[0]}\t{array[1]}\t{array[2]}"
            except (ValueError, IndexError) as e:
                # Handle parsing errors gracefully
                print(f"Skipping line due to error: {e}")
                continue

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
    finally:
        if file is not None:  # Ensure the file is closed properly
            file.close()

    return hit_dict








def parse_domain_table(file_path, gene_dict):
    """
    Parses a domain table file and updates the gene dictionary with domain annotation details.

    Parameters:
        file_path (str): Path to the file containing domain annotations.
        gene_dict (dict): Dictionary containing genes to be updated.

    Returns:
        dict: Updated gene dictionary with domain annotation details.
    """
    start_string = "Domain annotation for each sequence:"
    start_reading = False
    current_gene = None

    try:
        file = open(file_path, 'r')  # Explicitly open the file
        for line in file:
            line = line.strip()

            # Start processing lines after detecting the start_string
            if start_reading:
                if line.startswith(">> "):
                    # Extract gene name
                    current_gene = line[3:]
                    if current_gene not in gene_dict:
                        # Skip further processing if gene is not in the dictionary
                        continue

                # Skip comments, separators, and blank lines
                if line.startswith("#") or line.startswith("---") or not line:
                    continue

                # Extract and process domain coordinates
                parts = line.split()
                if current_gene and len(parts) > 8:
                    try:
                        # Convert coordinates to floats
                        alifrom = float(parts[6])
                        alito = float(parts[7])

                        # Update gene dictionary with domain information
                        if current_gene in gene_dict:
                            gene_dict[current_gene] += f"\t{int(alifrom)}\t{int(alito)}"
                            current_gene = None  # Reset after processing
                    except ValueError:
                        # Remove gene entry if domain coordinates are invalid
                        if current_gene in gene_dict:
                            del gene_dict[current_gene]
                        print(f"Error: Removed {current_gene} due to invalid domain coordinates: {line}")
                        current_gene = None

            # Enable processing when the start_string is found
            elif start_string in line:
                start_reading = True

        file.close()  # Explicitly close the file

    except Exception as e:
        print(f"Error: {e}")
    finally:
        if file:
            file.close()

    return gene_dict

def deconcat_hmmreport_to_per_domain_noise_cutoff_report(concat_filepath):
    """
    Deconcatenates a concatenated .hmmreport file into subfiles based on domain type,
    writing output files to the same directory as the input file.
    
    Args:
        concat_filepath (str): Path to the concatenated .hmmreport file.
    
    Returns:
        dict: A dictionary where keys are domain types and values are file paths.
    """
    # Determine the directory of the input file
    output_dir = os.path.dirname(concat_filepath)
    
    # Dictionary to store file handles for each domain
    domain_files = {}

    try:
        # Process the file line by line
        with open(concat_filepath, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) < 2:
                    continue  # Skip malformed lines
                
                domain_type = columns[1]
                if domain_type not in domain_files:
                    # Open a new file for this domain if it doesn't exist
                    domain_filename = os.path.join(output_dir, f"combined_{domain_type}.cleaned_report")
                    domain_files[domain_type] = open(domain_filename, 'a')
                
                # Write the line to the corresponding file
                domain_files[domain_type].write(line)
        
    finally:
        # Ensure all file handles are closed
        for f in domain_files.values():
            f.close()

    # Return the mapping of domain types to file paths
    return {domain: handle.name for domain, handle in domain_files.items()}

def remove_trusted_hits(score_threshold_dict, hmmreport_tsv):
    """
    Entfernt Treffer aus der HMM-Report-Tabelle, deren Score über dem definierten Cutoff-Wert liegt.
    
    Parameters:
    - score_threshold_dict (dict): Ein Dictionary mit {domain: cutoff_score}.
    - hmmreport_tsv (str): Der Dateipfad zur TSV-Datei.
    
    Returns:
    - str: Pfad zur gefilterten Datei.
    """
    output_path = hmmreport_tsv.replace(".noise_report", ".cleaned_report")
    
    with open(hmmreport_tsv, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) < 4:
                continue  # Zeilen ignorieren, die nicht genügend Spalten haben
            
            domain = columns[1]
            try:
                score = float(columns[3])
            except ValueError:
                continue  # Überspringen, falls der Score nicht konvertierbar ist
            
            if score < score_threshold_dict.get(domain, float('inf')):
                outfile.write(line)
    
    return output_path
    
############################ Collect genomeIDs ###################################

def collect_genomeIDs(report_path, options, div='___'):
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
                genomeID = line.split('\t')[0].split(div)[0]
                if genomeID in options.faa_files:
                    genomeIDs.add(genomeID)
    
    return genomeIDs

################# Prepare the genomeID sets for parallel processes ################
def split_into_batches(data_list, num_batches):
    batch_size = len(data_list) // num_batches
    return [data_list[i * batch_size:(i + 1) * batch_size] for i in range(num_batches)]
    
def process_batch(args):
    data_queue, genome_ids, options, counter = args
    for genome_id in genome_ids:
        try:
            process_genome(data_queue, genome_id, options)
        except Exception as e:
            print(f"Error processing {genome_id}: {e}")

def process_genome(data_queue, genome_id, options):

    protein_dict = dict()
    cluster_dict = dict()
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genome_id])
        gff_file = myUtil.unpackgz(options.gff_files[genome_id])

        protein_dict = parse_bulk_HMMreport_genomize(genome_id, global_hmmreport_tsv, global_score_threshold_dict, options.thrs_score)
        ParseReports.parseGFFfile(gff_file, protein_dict)
        ParseReports.getProteinSequence(faa_file, protein_dict)

        cluster_dict = Csb_finder.find_syntenicblocks(genome_id, protein_dict, options.nucleotide_range)
        Csb_finder.name_syntenicblocks(global_csb_patterns_dict,global_csb_pattern_names,cluster_dict,options.min_completeness)

        #synteny completion with provided patterns
        missing_domains = Csb_finder.extract_missing_domains_and_coordinates(cluster_dict) #e.g. structure {'SPIREOTU_00146859_4': {'oxDsrMK': {'missing_domains': {'oxDsrO', 'oxDsrP'}, 'start': 1362, 'end': 7719}}}
        
        synteny_completion_candidate_dict = process_missing_domains(genome_id, missing_domains, gff_file, faa_file, options.nucleotide_range)
        
        #Add the candidate proteins to the submitting protein_dict and cluster_dict
        #Add candidate to protein_dict if key is not already present int he protein dict
        for proteinID, protein in synteny_completion_candidate_dict.items():
            if proteinID not in protein_dict:
                protein_dict[proteinID] = protein
                cluster = cluster_dict[protein.clusterID]
                cluster.add_gene(proteinID, protein.get_domains(), protein.gene_start, protein.gene_end)
        
        data_queue.put((protein_dict, cluster_dict))
    except Exception as e:
        error_message = f"\nError: occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")

def process_missing_domains(genomeID, missing_domains, gff_file, faa_file, nt_range=3500):
    """
    Process missing domains by identifying candidate proteins, parsing their features,
    and checking if they can complete the gene clusters.

    Args:
        missing_domains (dict): Dictionary of missing domains for gene clusters.
        gff_file (str): Path to the GFF file.
        global_deconcat_domains_report_dict (dict): Dictionary mapping domain names to report file paths.

    Returns:
        dict: Updated missing_domains with candidate proteins that complete the clusters.
    """
    # Step 1: Prepare a dictionary for protein information based on missing domains
    candidate_protein_dict = {}
    insert_protein_dict = {}
    for csb_id, csb_data in missing_domains.items():
        for domain, domain_data in csb_data.items():
            if 'missing_domains' not in domain_data:
                continue
            for missing_domain in domain_data['missing_domains']: #For each missing domain get all distant hits from the concatenated report for this domain
                if missing_domain in global_deconcat_domains_report_dict:
                    domain_report_file = global_deconcat_domains_report_dict[missing_domain]
                    candidate_protein_dict = parse_bulk_HMMreport_genomize(genomeID,domain_report_file,{},10,candidate_protein_dict) # search for can. proteins that might fill the missing one

    # Step 2: Add features to proteins
    ParseReports.parseGFFfile(gff_file, candidate_protein_dict)

    
    # Step 3: Filter candidates and assign to clusters
    # only take candidates that are within range of the gencluster with the pattern that misses something
    #e.g. structure {'SPIREOTU_00146859_4': {'oxDsrMK': {'missing_domains': {'oxDsrO', 'oxDsrP'}, 'start': 1362, 'end': 7719}}}
    for csb_id, csb_data in missing_domains.items():
        for keyword_values in csb_data.values():
            cluster_start = keyword_values.get('start')
            cluster_end = keyword_values.get('end')
            missing_types = keyword_values.get('missing_domains')
            for protein_id, protein_info in candidate_protein_dict.items():
                gene_start = protein_info.gene_start
                gene_end = protein_info.gene_end

                # Check if the protein falls within the range or 3500 nt outside of the gene cluster
                if (cluster_start - nt_range <= gene_start <= cluster_end + nt_range) or \
                   (cluster_start - nt_range <= gene_end <= cluster_end + nt_range):
                    for domain in protein_info.get_domain_set():
                        if domain in missing_types:
                            # Update protein dictionary data with this protein
                            protein_info.clusterID = csb_id # assign clusterID to the protein
                            insert_protein_dict[protein_id] = protein_info

    # Step 4: Add sequences to the insertion candidates
    ParseReports.getProteinSequence(faa_file, candidate_protein_dict)
    
    return insert_protein_dict
    
                  
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

    #Submit the remaining and file reports
    if protein_batch or cluster_batch:
        submit_batches(protein_batch, cluster_batch, options)
        if options.individual_reports:
            first_protein_key = next(iter(protein_dict))  # Get the first key
            genomeID = protein_dict[first_protein_key].genomeID
            filepath = os.path.join(options.fasta_initial_hit_directory, str(genomeID))
            Output.output_genome_report(filepath, protein_dict, cluster_dict)
    return

def submit_batches(protein_batch, cluster_batch, options):
    # Submit the batches to the database and write to the file
    print(f"Submitted batch received")
    # Insert into the database
    Database.insert_database_proteins(options.database_directory, protein_batch)
    Database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')
    print(f"Submitted batch written to db and gene cluster file")



def parse_bulk_HMMreport_genomize(genomeID,Filepath,Thresholds,cut_score=10, protein_dict=None):
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
    if protein_dict is None:
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
    

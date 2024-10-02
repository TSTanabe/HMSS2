#!/usr/bin/python
import traceback
import multiprocessing
from multiprocessing import Manager, Pool

from . import Csb_finder
from . import Database
import sys

######## MAIN ROUTINE ##########

#### redo csb finding ##########

def redo_csb_reports(options):
    genomeIDs = Database.fetch_genomeIDs(options.database_directory)
    
    open(options.gene_clusters_file, 'w').close()
    
    # Parallel process all genomes
    manager = Manager()
    data_queue = manager.Queue()
    counterA = manager.Value('i', 0)
    counterB = manager.Value('i', 0)

    csb_patterns_diction, csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)

    with Pool(processes=1) as pool:
        # Process the genomes in parallel
        pool.map(process_parallel_bulk_parse, [(data_queue, genomeID, options, csb_patterns_diction, csb_pattern_names, counterA) for genomeID in genomeIDs])

    # After all processes are done, start the writer process
    process_writer(data_queue, options, counterB)

    print("Finished csb prediction")
    return    

def process_parallel_bulk_parse(args_tuple):
    data_queue, genomeID, options, csb_patterns_diction, csb_pattern_names, counter = args_tuple
    cluster_dict = {}
    protein_dict = {}
    try:

        # Fetch data from the database
        protein_dict = Database.fetch_protein_dict(options.database_directory, genomeID,protein_dict)

        # Find the syntenic regions and insert into the data structure
        cluster_dict = Csb_finder.find_syntenicblocks(genomeID, protein_dict, options.nucleotide_range)
        Csb_finder.name_syntenicblocks(csb_patterns_diction, csb_pattern_names, cluster_dict, options.min_completeness)

        data_queue.put(cluster_dict)

    except Exception as e:
        error_message = f"\nError: occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Skipped {genomeID} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
    counter.value += 1
    print(f"Processed {counter.value} genomes ", end="\r")
    
    return

def process_writer(queue, options, counter):
    """
    This routine handles the output of the search and writes it into the database.
    It gets input from multiple workers as the database connection to sqlite is unique.
    """
    cluster_batch = {}
    batch_size = options.glob_chunks
    batch_counter = 0

    while not queue.empty():
        cluster_dict = queue.get()
        if cluster_dict is None:
            break

        counter.value += 1
        print(f"Inserted {counter.value} genomes ", end="\r")

        cluster_batch.update(cluster_dict)
        batch_counter += 1

        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            submit_batches(cluster_batch, options)
            cluster_batch = {}
            batch_counter = 0

    # Process any remaining data in the batch
    if cluster_batch:
        submit_batches(cluster_batch, options)

    return

def submit_batches(cluster_batch, options):
    # Insert all collected clusters into the database
    
    
    genome_ids_to_delete = {cluster_obj.genomeID for cluster_obj in cluster_batch.values()}
    Database.delete_database_clusters_by_genomeID(options.database_directory, genome_ids_to_delete)
    Database.drop_specified_indices(options.database_directory)
    Database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')

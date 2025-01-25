#!/usr/bin/python
import os
import traceback
import pprint

from . import myUtil
from . import Search
from . import ParseReports
from . import Csb_finder
from . import Database
from . import Output

from collections import defaultdict
from multiprocessing import Manager, Pool, Value

# Global dictionaries to be shared as read-only
score_threshold_diction = {}
csb_patterns_diction = {}
csb_pattern_names = {}

def multi_search_process(options):
    global score_threshold_diction, noise_threshold_diction, csb_patterns_diction, csb_pattern_names
    
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type)
    noise_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, 3)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
   
    process_instances = []
    
    with Pool(processes = options.cores) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options))
        
        # Use an iterator for memory efficiency and avoid creating a large argument list
        tasks = ((data_queue, genomeID, options, counter)
                 for genomeID in options.queued_genomes)
                 
        pool.map(process_parallel_search, tasks)

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
        
    print("\nFinished searching")
    print(f"Saved hits to database {options.database_directory}")
    print(f"Saved individual hit lists to {options.fasta_initial_hit_directory}")
    return    
        
     


def process_parallel_search(args_tuple):
    global score_threshold_diction, noise_threshold_diction, csb_patterns_diction, csb_pattern_names
    
    queue,genomeID,options,counter = args_tuple
    counter.value += 1
    print(f"Searching assembly ({counter.value}/{len(options.queued_genomes)})", end="\r")
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        
        hmm_report = Search.HMMsearch(faa_file,options.library,options,1)
        protein_dict = ParseReports.parseHMMreport(hmm_report,score_threshold_diction,options.thrs_score)
        ParseReports.parseGFFfile(gff_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,protein_dict)

        cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,options.min_completeness)

        missing_domains = Csb_finder.extract_missing_domains_and_coordinates(cluster_dict) #e.g. structure {'SPIREOTU_00146859_4': {'oxDsrMK': {'missing_domains': {'oxDsrO', 'oxDsrP'}, 'start': 1362, 'end': 7719}}}
        synteny_completion_candidate_dict = process_missing_domains_single_genome(genomeID, missing_domains, gff_file, faa_file, hmm_report, options.nucleotide_range)

        for proteinID, protein in synteny_completion_candidate_dict.items():
            if proteinID not in protein_dict:
                protein_dict[proteinID] = protein
                cluster = cluster_dict[protein.clusterID]
                cluster.add_gene(proteinID, protein.get_domains(), protein.gene_start, protein.gene_end)
        
#            if options.synthenic_block_support_detection:
#                missing_protein_types, missing_proteins_list_dict = Csb_finder.find_csb_pattern_difference(csb_patterns_diction,csb_pattern_names,cluster_dict,3)
#                candidate_proteins_dict = ParseReports.parseHMMreport_below_cutoff_hits(missing_protein_types,hmm_report,score_threshold_diction,options.cut_score)
#                Csb_finder.synteny_completion(gff_file,protein_dict,cluster_dict,candidate_proteins_dict,missing_proteins_list_dict,options.nucleotide_range)
        
        queue.put((genomeID,protein_dict,cluster_dict))
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
        return

def process_missing_domains_single_genome(genomeID, missing_domains, gff_file, faa_file, hmm_report, nt_range=3500):
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

    candidate_protein_dict = ParseReports.parseHMMreport(hmm_report,noise_threshold_diction)
                


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

     
def process_writer(queue,options):
    #This routine handles the output of the search and writes it into the database
    #It gets input from multiple workers as the database connection to sqlite is unique
    while True:
        tup = queue.get()
        if tup is None :
            break

        genomeID,protein_dict,cluster_dict = tup
        
        Database.insert_database_genomeID(options.database_directory,genomeID)
        Database.insert_database_protein(options.database_directory,genomeID,protein_dict)
        Database.insert_database_cluster(options.database_directory,genomeID,cluster_dict)
        
        with open(options.gene_clusters_file, "a") as file:
                for clusterID, cluster in cluster_dict.items():
                    domains = cluster.get_domains()
                    file.write(clusterID + '\t' + '\t'.join(domains) + '\n')
        
        if options.individual_reports:
            filepath = os.path.join(options.fasta_initial_hit_directory, str(genomeID))
            Output.output_genome_report(filepath, protein_dict, cluster_dict,{} , genomeID)
        
    return
        

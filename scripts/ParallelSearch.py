#!/usr/bin/python
import os
import traceback

from . import myUtil
from . import Search
from . import ParseReports
from . import Csb_finder
from . import Database
from . import Output

import multiprocessing
from multiprocessing import Process, Manager, Pool

def multi_search_process(options):
    
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
   
    process_instances = []
    
    with Pool(processes = options.cores) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options))
        
        pool.map(process_parallel_search, [(data_queue, genomeID, options, counter, score_threshold_diction, csb_patterns_diction,csb_pattern_names) for genomeID in options.queued_genomes])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
    print("\nFinished searching")
    print(f"Saved hits to database {options.database_directory}")
    print(f"Saved individual hit lists to {options.fasta_initial_hit_directory}")
    return    
        
     


def process_parallel_search(args_tuple):
    
    queue,genomeID,options,counter,score_threshold_diction, csb_patterns_diction,csb_pattern_names = args_tuple
    counter.value += 1
    print(f"Searching assembly ({counter.value}/{len(options.queued_genomes)})", end="\r")
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        hmm_report = Search.HMMsearch(faa_file,options.library,options,1)
        protein_dict = ParseReports.parseHMMreport(hmm_report,score_threshold_diction,options.thrs_score)
            
        #complete the hit information
        ParseReports.parseGFFfile(gff_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,protein_dict)
            
            
            
        #find the syntenic regions and insert to database
        if not options.redo_search:
            cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
            Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,options.min_completeness)


            if options.synthenic_block_support_detection:
                missing_protein_types, missing_proteins_list_dict = Csb_finder.find_csb_pattern_difference(csb_patterns_diction,csb_pattern_names,cluster_dict,3)
                candidate_proteins_dict = ParseReports.parseHMMreport_below_cutoff_hits(missing_protein_types,hmm_report,score_threshold_diction,options.cut_score)
                Csb_finder.synteny_completion(gff_file,protein_dict,cluster_dict,candidate_proteins_dict,missing_proteins_list_dict,options.nucleotide_range)
        
        queue.put((genomeID,protein_dict,cluster_dict))
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
        return

    
    
    
     
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
        

     
     
     
     
     
#If only 1 CPU is available        
def initial_search(options):
    
    query = options.library
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    for index,genomeID in enumerate(options.queued_genomes):
        
        print(f"Processing assembly {index+1} of {len(options.queued_genomes)}",end="\r") #print control sequence output
        #Prepare names
        try:
            #Find hits and parse results
            faa_file = myUtil.unpackgz(options.faa_files[genomeID])
            gff_file = myUtil.unpackgz(options.gff_files[genomeID]) 
            hmm_report = Search.HMMsearch(faa_file,options.library,options,options.cores)
            protein_dict = ParseReports.parseHMMreport(hmm_report,score_threshold_diction,options.cut_score)
            
            #complete the hit information
            ParseReports.parseGFFfile(gff_file,protein_dict)
            ParseReports.getProteinSequence(faa_file,protein_dict)
            
            
            
            #find the syntenic regions and insert to database
            cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
            Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,options.min_completeness)


            if options.synthenic_block_support_detection:
                missing_protein_types, missing_proteins_list_dict = Csb_finder.find_csb_pattern_difference(csb_patterns_diction,csb_pattern_names,cluster_dict,3)
                candidate_proteins_dict = ParseReports.parseHMMreport_below_cutoff_hits(missing_protein_types,hmm_report,score_threshold_diction,options.cut_score)
                Csb_finder.synteny_completion(gff_file,protein_dict,cluster_dict,candidate_proteins_dict,missing_proteins_list_dict,options.nucleotide_range)

            #Insert to database
            Database.insert_database_genomeID(options.database_directory,genomeID)
            Database.insert_database_protein(options.database_directory,genomeID,protein_dict)
            Database.insert_database_cluster(options.database_directory,genomeID,cluster_dict)
            
            #Write for csb finder
            with open(options.gene_clusters_file, "a") as file:
                for clusterID, cluster in cluster_diction.items():
                    domains = cluster.get_domains()
                    file.write(clusterID + '\t' + '\t'.join(domains) + '\n')
            Output.output_genome_report(options.fasta_initial_hit_directory+"/"+genomeID,protein_diction,cluster_diction)
        except Exception as e:
            error_message = f"\nError occurred: {str(e)}"
            print(f"\t WARNING: Skipped {faa_file} due to an error - {error_message}")
            continue

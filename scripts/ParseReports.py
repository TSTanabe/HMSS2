#!/usr/bin/python
import os
import re
import subprocess
import traceback
import itertools
import pprint

from . import Database
from . import Csb_finder
from . import myUtil
from . import Output

import numpy as np

from multiprocessing import Pool, Manager
from typing import Dict, Optional, Any
from collections import defaultdict
from scipy.optimize import linear_sum_assignment

logger = myUtil.logger

class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have a higher score than all overlapped domains. If new domain has a lower score than any previously added domain new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein. Here the additional information of other domains is added if it does not interfere with this assumption
    
    Organizes protein domains and related attributes for a protein. 
    New domains are only added if they do not overlap with higher-scoring existing domains.
    All coordinates, scores, and HMM names are accessible as dash-separated strings.

    Args:
        proteinID (str): Unique identifier.
        HMM (str): Domain name.
        start (int): Start coordinate.
        end (int): End coordinate.
        score (float): Domain bitscore.
        genomeID (str): Genome identifier (optional).
        ident (int): Percent identity (default: 25).
        bsr (float): Blast score ratio (default: 1.0).

    Example:
        p = Protein("Prot1", "HMM_A", 10, 120, 40)
        p.add_domain("HMM_B", 130, 200, 30)
    """


    def __init__(
        self,
        proteinID: str,
        HMM: str,
        start: int = 0,
        end: int = 0,
        score: float = 1,
        genomeID: str = "",
        ident: int = 25,
        bsr: float = 1.0
    ):
        self.proteinID: str = proteinID
        self.genomeID: str = genomeID
        self.protein_sequence: str = ""
        self.gene_contig: str = ""
        self.gene_start: int = 0
        self.gene_end: int = 0
        self.gene_strand: str = "."
        self.gene_locustag: str = ""
        self.clusterID: str = ""
        self.keywords: Dict = {}
        self.domains: Dict[int, Domain] = {}  # start coordinate → Domain object
        self.add_domain(HMM, start, end, score, ident, bsr)
        

    ##### Getter ####
            
    def get_domains(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_HMM())
        return '-'.join(listing)
    
    def get_domains_dict(self):
    #return dict
        return self.domains
            
    def get_domain_listing(self):
    #return list
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key])
        return listing
        
    def get_domain_set(self):
        domains = set()
        for v in self.domains.values():
            domains.add(v.get_HMM())
        return domains
        
    def get_domain_coordinates(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_start()}:{self.domains[key].get_end()}")
        return '-'.join(listing)

    def get_domain_scores(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_score()}")
        return '-'.join(listing)   

    def get_domain_count(self):
        return len(self.domains)
        
    def get_protein_string(self):
        #3.9.22 representation of the whole protein in one line
        a = self.proteinID
        b = self.get_domains()
        c = self.get_domain_scores()
        d = self.gene_contig
        e = self.gene_start
        f = self.gene_end
        #d = self.get_protein_sequence()
        string = f"{a} {b} {c} {d} {e} {f}"
        return string
        
    def get_protein_list(self):
        #3.9.22 representation of the whole protein in one line
        listing = []
        listing.append(self.proteinID)
        listing.append(self.get_domains())
        listing.append(str (self.get_domain_scores()))
        listing.append(str (self.get_domain_coordinates()))
        listing.append(self.gene_contig)
        listing.append(str (self.gene_start))
        listing.append(str (self.gene_end))
        listing.append(self.gene_strand)
        listing.append(self.gene_locustag)
        #d = self.get_protein_sequence()
        #string = f"{a} {b} {c} {d} {e} {f}"
        return listing
            
    def get_sequence(self):
        return str(self.protein_sequence)
        
    def get_length(self):
        return len(self.protein_sequence)
    ##### Setter #####

    def check_domain_overlap(self,new_start, new_end,\
    current_start,current_end):
    #2.9.22

        if (current_start <= new_start and new_start <= current_end)\
        or (current_start <= new_end and new_end <= current_end):
        #start oder endpunkt innerhalb er grenzen
            return 1

        elif (new_start <= current_start and current_end <= new_end)\
        or (current_start <= new_start and new_end <= current_end):
        #start und end innerhalb der grenzen oder alte domäne innerhalb der neuen
            return 1
            
        elif (new_start <= current_start and new_end <= current_start)\
        or (new_start >= current_end and new_end >= current_end):
        #start und end kleiner als self start oder start und end größer als self end dann
        #domäne außerhalb der alten domäne und adden egal welcher score
            return 0
        return None
        

    def add_domain(
        self,
        HMM: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0
    ) -> int:
        """
        Adds a domain to the protein only if it does not overlap
        with a higher-scoring existing domain. If overlap exists with lower-scoring
        domain, that domain is removed.

        Returns:
            int: 1 if domain added, 0 if not added.
        """

        del_domains = [] # start coordinates/keys of domains to be replace
        for domain in self.domains.values():
            if self.check_domain_overlap(start,end,domain.get_start(),domain.get_end()):
                if domain.get_score() < score:
                    del_domains.append(domain.get_start())
                else:
                    return 0


        
        for key in del_domains:
            self.domains.pop(key)
        self.domains.update({start:Domain(HMM,start,end,score,ident,bsr)}) # if loop complete
        
        return 1
        
        
        


class Domain:
#2.9.22
    """
    Stores domain information (HMM name, coordinates, score, identity, bsr).

    Args:
        HMM (str): Domain name.
        start (int): Start coord.
        end (int): End coord.
        score (float): Bitscore.
        ident (int): Percent identity.
        bsr (float): Blast score ratio.
    """
    def __init__(self, HMM: str, start: int, end: int, score: float, ident: int = 1, bsr: float = 1.0):
        self.HMM: str = HMM
        self.start: int = int(start)
        self.end: int = int(end)
        self.score: float = float(score)
        self.identity: int = int(ident)
        self.bsr: float = float(bsr)
    
    def __hash__(self):
        return hash((self.HMM, self.start, self.end, self.score))
    
    def __eq__(self, other):
        if isinstance(other, Domain):
            return self.HMM == other.HMM and self.start == other.start and self.end == other.end and self.score == other.score
        return False
    
    def get_HMM(self):
        return self.HMM
    def get_start(self):
        return self.start
    def get_end(self):
        return self.end
    def get_score(self):
        return self.score



#########################################
########   Parsing subroutines ##########
#########################################


def parseGFFfile(
    filepath: str, protein_dict: Dict[str, Protein]
) -> Dict[str, Protein]:
    """
    3.9.22
    Adds GFF attributes to each Protein object in the dictionary.

    Args:
        filepath (str): Path to GFF3 file.
        protein_dict (dict): {proteinID: Protein object}

    Returns:
        dict: Updated protein_dict.
    """
    locustag_pattern = re.compile(r'locus_tag=(\S*?)(?:[;\s]|$)')
    geneID_pattern = re.compile(r'ID=(cds-)?(\S+?)(?:[;\s]|$)')
    
    grep_pattern = "|".join(protein_dict.keys())
    try:
        grep_process = subprocess.Popen(['grep', '-E', grep_pattern, filepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = grep_process.communicate()

        if stderr:
            logger.error("Grep process:", stderr)
            return protein_dict
        
        for line in stdout.decode('utf-8').split('\n'):
            if not line:
                continue
            gff = line.split("\t")
            match = geneID_pattern.search(gff[-1])
            if not match:
                continue
            match = match.group(2)
            if match in protein_dict:
                # Add protein
                protein = protein_dict[match]
                
                protein.gene_contig = str(gff[0])
                protein.gene_start = int(gff[3])
                protein.gene_end = int(gff[4])
                protein.gene_strand = str(gff[6])
                locustag = getLocustag(locustag_pattern, line)
                protein.gene_locustag = str(locustag)
    except Exception as e:
        logger.warning(f"Error occurred while parsing GFF: {str(e)}")
        return protein_dict
    return protein_dict





def get_protein_sequence(filepath, protein_dict):
    """
    Add protein sequences from a FASTA file to Protein objects in a dictionary.

    Iterates through a FASTA file of amino acid sequences and sets the `.protein_sequence`
    attribute for each Protein object in `protein_dict`, keyed by their protein ID (header).
    Only protein IDs present in `protein_dict` are processed.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file containing amino acid sequences.
    protein_dict : dict
        Dictionary with protein IDs as keys and Protein objects as values.

    Returns
    -------
    dict
        The updated protein_dict with sequences added to the Protein objects.

    Examples
    --------
    >>> protein_dict = {'WP_0001': Protein(...), 'WP_0002': Protein(...)}
    >>> get_protein_sequence('proteins.faa', protein_dict)
    {'WP_0001': <Protein with sequence>, 'WP_0002': <Protein with sequence>}
    """
    reader = None
    try:
        reader = open(filepath, "r")
        sequence = ""
        header = None
        save_sequence = False

        for line in reader:
            line = line.strip()
            if line.startswith(">"):
                if header and save_sequence and sequence:
                    # Save sequence for previous protein
                    protein = protein_dict[header]
                    protein.protein_sequence = sequence
                # Parse header up to first whitespace
                header = line[1:].split()[0]
                if header in protein_dict:
                    save_sequence = True
                    sequence = ""
                else:
                    save_sequence = False
            elif save_sequence:
                sequence += line

        # Handle the last protein
        if header and save_sequence and sequence:
            protein = protein_dict[header]
            protein.protein_sequence = sequence

    except IOError as e:
        logger.error(f"Cannot open {filepath}: {e}")
    finally:
        if reader is not None:
            reader.close()

    return protein_dict


def getLocustag(locustag_pattern: re.Pattern, string: str) -> str:
    """
    Extracts locus_tag from a string using a regex pattern.

    Args:
        locustag_pattern (re.Pattern): Regex pattern for locus_tag.
        string (str): Line to search.

    Returns:
        str: locus_tag or empty string if not found.
    """
    match = locustag_pattern.search(string)
    return match.group(1) if match else ""




###############################################################################################################
###############################################################################################################
###############################################################################################################

# Writer for hit reports

def submit_batches(protein_batch, cluster_batch, options):

    # Insert into the database
    Database.insert_database_proteins(options.database_directory, protein_batch)
    Database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')



def process_writer(queue, options):
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
            batch_counter += 1
            print(f"Processed {batch_counter} genomes ", end="\r")#
        
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
                filepath = os.path.join(options.fasta_initial_hit_directory, str(genomeID)+".hit_table_txt")
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
            if protein_dict:
                first_protein_key = next(iter(protein_dict))  # Get the first key
                genomeID = protein_dict[first_protein_key].genomeID
                filepath = os.path.join(options.fasta_initial_hit_directory, str(genomeID)+".hit_table_txt")

                Output.output_genome_report(filepath, protein_dict, cluster_dict)
    return
    
    
#########################################################################################
################ Processing routines for parsing genome hits ############################
#########################################################################################

def main_parse_summary_hmmreport(options):


    genome_ids = list(options.queued_genomes)
    genomeID_batches = split_into_batches(genome_ids, options.cores - 1)

    # Lade Patterns nur 1x im Hauptprozess
    csb_patterns, csb_names = Csb_finder.make_pattern_dict(options.patterns_file)
    
    # Insert genomeIDs in DB
    Database.insert_database_genomeIDs(options.database_directory, set(genome_ids))

    with Manager() as manager:
        data_queue = manager.Queue()

        with Pool(processes=options.cores) as pool:
            # Start writer
            p_writer = pool.apply_async(process_writer, (data_queue, options))
            args = [
                (
                    data_queue,
                    batch,
                    options.faa_files,
                    options.gff_files,
                    options.glob_trusted_hitreport,
                    options.glob_intermediate_hitreport,
                    options.nucleotide_range,
                    options.min_completeness,
                    csb_patterns,
                    csb_names
                )
                for batch in genomeID_batches
            ]

            # Start readers
            pool.starmap(process_batch, args)

            # End writer
            for _ in range(options.cores):
                data_queue.put(None)

            p_writer.get()
    logger.info("Finished parsing of search results and writing to local database")


def split_into_batches(data_list, num_batches):
    if num_batches <= 0:
        raise ValueError("Number of batches must be > 0.")
    if len(data_list) == 0:
        return [[] for _ in range(num_batches)]
    
    batches = [[] for _ in range(num_batches)]
    for idx, item in enumerate(data_list):
        batches[idx % num_batches].append(item)
    return batches


def process_batch(
    data_queue, genome_ids, faa_files, gff_files, hmmreport_path, intermediate_hmmreport,
    nucleotide_range, min_completeness, pattern_dict, pattern_names
):
    for genome_id in genome_ids:
        try:
            process_genome(
                data_queue, genome_id,
                faa_files[genome_id], gff_files[genome_id], hmmreport_path, intermediate_hmmreport,
                nucleotide_range, min_completeness, pattern_dict, pattern_names
            )
        except Exception as e:
            logger.warn(f"Failed to process genome '{genome_id}' — {str(e)}")
            continue
    
    
        
   
    
def process_genome(
    data_queue, genome_id, faa_path, gff_path, hmmreport_path, intermediate_hmmreport,
    nucleotide_range, min_completeness, pattern_dict, pattern_names
):
    try:
        faa_file = myUtil.unpackgz(faa_path)
        gff_file = myUtil.unpackgz(gff_path)
        
        # Get the intermediate hits
        intermediate_protein_dict = parse_bulk_HMMreport_genomize(genome_id, intermediate_hmmreport)
        parseGFFfile(gff_file, intermediate_protein_dict)
                
        # Get the initial hit information
        protein_dict = parse_bulk_HMMreport_genomize(genome_id, hmmreport_path)
        parseGFFfile(gff_file, protein_dict)
        
        # Recognition of named gene clusters
        combined_protein_dict = {**intermediate_protein_dict,**protein_dict}
        cluster_dict = Csb_finder.find_syntenic_blocks(genome_id, combined_protein_dict, nucleotide_range)
        
        # Name syntenic blocks with known patterns
        cluster_dict = Csb_finder.name_syntenic_blocks(pattern_dict, pattern_names, cluster_dict, min_completeness)
	


        # Enhance the completeness of the dsb by swapping assignment
        enhance_syntenic_block_completeness(cluster_dict, combined_protein_dict, intermediate_protein_dict, intermediate_hmmreport, pattern_dict)
        
        # Remove intermediate hits that are not part of a named gene cluster pattern
        #TODO soll nicht removen was gerade noch geändert wurde
        combined_protein_dict = remove_unassigned_intermediate_proteins(combined_protein_dict, protein_dict, cluster_dict)
            # for each cluster
            # if completeness is  0.5 <= completeness < 1, wobei 0.5 hier durch das min(0.5 oder min_completeness) ersetzt werden sollte
            # falls es zusätzliche gibt dann schau nach
            # ob die proteinIDs in intermediate existieren (trusted sollten nicht verändert werden)
            # hole alle domänen einträge für diese proteinID
            
            # kann man eine oder mehrere in proteine umwandeln die in missing enthalten sind
            # zähle die anzahl der wechsel
            # a b c d e x g y 
            # missing: f, h, additional: x, y,
            # a b c d e x->f g y->h
            
            # den block mit der höchsten unvollständigen completeness nehmen, falls mehrere alle davon nehmen
            # falls exakter match mit completeness 1 und missing 0 und additional 0 dann skip
            # für csb(s) mit höchster completeness und gleicher anzahl an missing und additional oder möglichst gleicher anzahl
            # zählen wieviele austausche gemacht werden können
            
            # kalkuliere wie wievle scorepunkte eine änderung bedeuten würde und ob man in eine der fehlenden domänen ändern könnte
            # jetzt wird es kompliziert in der abfolge
            # ich will wissen wieviele änderungen es braucht um den vorhandenen csb in den pattern csb zu überführen oder anzunähren (maximiere durch mögliche änderungen hin zu mehr completeness)
            # wenn mehrere mit gleich vielen änderungen das maximum erreichen können, dann 
        
        
        # Add the protein sequences
        get_protein_sequence(faa_file, combined_protein_dict)

        data_queue.put((combined_protein_dict, cluster_dict))

    except Exception as e:

        logger.error(f"Error: {genome_id} -> {e}")
        logger.error(traceback.format_exc())


def remove_unassigned_intermediate_proteins(
    combined_protein_dict: Dict[str, Any],
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Remove proteins from combined_protein_dict that are not present in protein_dict
    or in the set of covered_protein_ids from any cluster.

    Parameters
    combined_protein_dict : dict
        All protein objects (combined sources; may contain intermediate-only proteins).
    protein_dict : dict
        All primary protein objects from main search.
    cluster_dict : dict
        Dictionary of cluster objects after naming. Each cluster should have a
        .covered_protein_ids attribute containing recognized protein IDs.

    Returns
    dict
        Filtered combined_protein_dict with only protein IDs that are present in the
        primary protein set or are covered by a recognized pattern in any cluster.

    """
    proteinIDs = set(protein_dict.keys())
    for cluster in cluster_dict.values():
        proteinIDs.update(getattr(cluster, 'covered_protein_ids', set()))

    # Remove any protein in combined_protein_dict that is not in proteinIDs
    for protein_id in list(combined_protein_dict.keys()):
        if protein_id not in proteinIDs:
            del combined_protein_dict[protein_id]

    # Remove unassigned genes and update clusters accordingly
    for clusterID in list(cluster_dict.keys()):  # list() so we can delete inside loop
        cluster = cluster_dict[clusterID]
        # Remove genes from cluster.genes and cluster.types that are not in combined_protein_dict
        if hasattr(cluster, 'genes'):
            # Remove genes not present anymore
            filtered_genes = []
            filtered_types = []
            for gene, typ in zip(cluster.genes, getattr(cluster, "types", [])):
                if gene in combined_protein_dict:
                    filtered_genes.append(gene)
                    filtered_types.append(typ)
            cluster.genes = filtered_genes
            if hasattr(cluster, "types"):
                cluster.types = filtered_types

            # If the cluster now has fewer than two proteins, remove the cluster entirely
            if len(cluster.genes) < 2:
                # Also remove the clusterID from the corresponding proteins
                for gene in cluster.genes:
                    if gene in combined_protein_dict:
                        combined_protein_dict[gene].clusterID = ""
                del cluster_dict[clusterID]

    return combined_protein_dict

    
    
def parse_bulk_HMMreport_genomize(genomeID,Filepath,protein_dict=None):
    """
 
    """
    if protein_dict is None:
        protein_dict = {}

    result = subprocess.run(['grep', genomeID, Filepath], stdout=subprocess.PIPE, text=True)
    
    lines = result.stdout.splitlines()  # Split output into lines

    for line in lines:
        columns = line.split('\t')  # Assuming columns are space-separated
        if columns:
            try:
                key = columns[0]
                genomeID, hit_proteinID = key.split('___',1)
                query = columns[3].split('_')[-1]
                hit_bitscore = int(float(columns[7]))
                hsp_start = int(float(columns[17]))
                hsp_end = int(float(columns[18]))
                
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore,genomeID)
            except Exception as e:
                error_message = f"\nError occurred: {str(e)}"
                traceback_details = traceback.format_exc()
                logger.error(f"Skipped {Filepath} due to an error - {error_message}")
                logger.error(f"Traceback details:\n{traceback_details}")
                continue
                
    return protein_dict





def enhance_syntenic_block_completeness(cluster_dict, combined_protein_dict, intermediate_protein_dict, intermediate_hmmreport, pattern_dict, min_completeness=0.5):

    """
    Enhance syntenic block completeness by swapping additional protein domains with missing ones if possible.
    for each cluster that has no directly matching pattern from the given patterns 
    it is tested if a possible conversion of protein types to alternative ones with lower hitscore could reach 
    a better completion
    """


    # Go through all clusters
    for cluster in cluster_dict.values():
        # Skip clusters that have a perfectly matching pattern
        if any(kw.get_completeness() == 1 and not kw.get_additional_domains() for kw in cluster.get_keywords()):
            continue

        
        genes = cluster.genes
        types = cluster.types
        
        # Dictionary to hold the hmmreport lines
        alternative_protein_type_dict = {}
        
        # Possible executions to increase the completeness for each keyword
        possible_optimized_executable_transitions = {}
        
        for keyword in cluster.get_keywords():
            completeness = keyword.get_completeness()
            missing_domains = keyword.get_missing_domains()
            additional_domains = keyword.get_additional_domains()
            transition_dict = defaultdict(set)
            # Structure: alternative_protein_type_dict[(hit_proteinID, query)] = (hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID)
            # Structure: transition_dict[query] => ((hit_proteinID, difference), (hit_proteinID, difference))
            #print(f"\nNew gene cluster {cluster.clusterID}")
            #print(types)
            #print(f"Processing {keyword.keyword} {completeness}")
            #print(f"Missing {missing_domains}")
            #print(f"Additional {additional_domains}")            
            if completeness >= min(0.5, min_completeness) and completeness < 1 and additional_domains:
                # For every "additional" domain: get proteinID and possible transitions to "missing"
                transitions = []
                for index, current_domain in enumerate(types):
                    #remove after debugging
                    proteinID = genes[index]
                    protein = combined_protein_dict[proteinID]
                    
                    # Get the proteinID of the additional domains
                    if current_domain in additional_domains and genes[index] in intermediate_protein_dict:
                        # Parse domtblout to get all potential domains and their scores                        
                        proteinID = genes[index]
                        protein = combined_protein_dict[proteinID]

                        # Updates the alternative_protein_type_dict and transition_dict
                        find_possible_transitions(
                            proteinID, current_domain, missing_domains, protein, intermediate_hmmreport,
                            alternative_protein_type_dict, transition_dict
                        )

                # Optimize the transitions by minimizing the bitscore changes and number of transitions to reach the missing domains
                pattern_length = len(pattern_dict[keyword.keyword_id])
                

                chosen_transitions, posterior_completeness, total_score_diff = get_optimal_transitions(transition_dict, missing_domains, completeness, pattern_length)
                """
                Example output
                chosen_transitions = [
                    ('P1', 'A', 50),
                    ('P2', 'B', 20),
                    ('P3', 'C', 40)
                ]
                completeness = 1.0
                total_score_diff = 110
                """
                if chosen_transitions:
                    possible_optimized_executable_transitions[(pattern_length, posterior_completeness, total_score_diff)] = chosen_transitions
        # Keyword loop finished
        
        
        # Now from all possible keyword completions find the optimum dictionary key (transitions, bitscore, posterior completeness) => [(proteinID to domain), (proteinID to domain), (proteinID to domain)]
        # Highest pattern length
        # Highest posterior completeness
        # Minimal total score difference
        if possible_optimized_executable_transitions:
            # Sort by: pattern_length (desc), posterior_completeness (desc), total_score_diff (asc)
            best_key = max(
                possible_optimized_executable_transitions.keys(),
                key=lambda x: (x[1], x[0], -x[2])  # pattern_length, posterior_completeness, -score_diff
            )
            best_transitions = possible_optimized_executable_transitions[best_key]

            logger.debug(f"For clusterID {cluster.clusterID} following conversion is done")
            logger.debug(best_transitions)

            # For the best transition alter the proteins domain information
            for proteinID, to_domain in best_transitions:
                if proteinID in combined_protein_dict:
                    protein = combined_protein_dict[proteinID]
                    hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID = alternative_protein_type_dict.get((proteinID, to_domain)) #(hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID)
                    protein.domains.clear()
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                
    return combined_protein_dict










def find_possible_transitions(
    proteinID, 
    current_domain, 
    missing_domains, 
    protein, 
    hmmreport,
    alternative_protein_type_dict,   
    transition_dict                 
):

    """
        Looks for up for a given proteinID if alternative hits are in the present in
        a given hmmreport in domtblout format with 18 columns.
        First the lines with the proteinID are grepped with grep
        then the lines are parsed. If query hmm is also fitting to a missing protein domain
        this alternative is saved with all values and the difference to the initial hitscore
        is calculated.
        Returned are the a dictionary with the values from the hmmreport line
        and a dictionary indicating the proteinID and the possible alternative protein type
    """
    result = subprocess.run(['grep', proteinID, hmmreport], stdout=subprocess.PIPE, text=True)
    lines = result.stdout.splitlines()

    for line in lines:
        columns = line.split('\t')
        if columns and len(columns) > 18:
            try:
                key = columns[0]
                genomeID, hit_proteinID = key.split('___', 1)
                query = columns[3].split('_')[-1]
                hit_bitscore = float(columns[7])
                hsp_start = int(float(columns[17]))
                hsp_end = int(float(columns[18]))

                if query in missing_domains:
                    alternative_protein_type_dict[(hit_proteinID, query)] = (
                        hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID
                    )
                    current_domain_score = next(
                        (domain.get_score() for domain in protein.domains.values() if domain.get_HMM() == current_domain),
                        0
                    )
                    
                    difference = abs(current_domain_score - hit_bitscore)
                    transition_dict[query].add((hit_proteinID, difference))
                    
            except Exception as e:
                logger.warn(f"Skipped line in {hmmreport} due to an error - {str(e)}")
                continue

    return alternative_protein_type_dict, transition_dict





def get_optimal_transitions(transition_dict, missing_domains, initial_completeness=0.0, total_domains=None):
    """
    Calculates the optimal set of transitions to cover as many missing domains as possible,
    each proteinID at most once, minimizing transitions and total score difference.

    transition_dict: {missing_domain: set of (proteinID, score_diff)}
    missing_domains: list or set of missing domains to fulfill

    Returns:
        chosen_transitions: [(proteinID, to_domain, score_diff)]
        completeness: number fulfilled / total
        total_score_diff: sum of chosen score diffs
    """
    # Gather all protein candidates
    proteins = set()
    for domain in missing_domains:
        for protein, diff in transition_dict[domain]:
            proteins.add(protein)
    proteins = list(proteins)
    n_proteins = len(proteins)
    n_domains = len(missing_domains)

    if n_proteins == 0 or n_domains == 0:
        return [], 0.0, 0.0

    # Build cost matrix
    cost_matrix = np.full((n_proteins, n_domains), np.inf)
    for j, domain in enumerate(missing_domains):
        for protein, diff in transition_dict[domain]:
            i = proteins.index(protein)
            cost_matrix[i, j] = diff

    # Hungarian assignment: minimize score difference (and thus minimize transitions)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    chosen_transitions = []
    total_score_diff = 0
    fulfilled_domains = set()
    used_proteins = set()
    missing_domains = list(missing_domains)
    
    for i, j in zip(row_ind, col_ind):
        cost = cost_matrix[i, j]
        if np.isfinite(cost):
            proteinID = proteins[i]
            domain = missing_domains[j]
            chosen_transitions.append((proteinID, domain))
            total_score_diff += cost
            fulfilled_domains.add(domain)
            used_proteins.add(proteinID)
    
    # Completeness for the final keyword
    final_completeness = initial_completeness
    if total_domains:
        final_completeness += len(fulfilled_domains) / total_domains
    else:
        final_completeness += len(fulfilled_domains) / max(len(missing_domains), 1)
    
    return chosen_transitions, final_completeness, total_score_diff












  


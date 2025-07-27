#!/usr/bin/python
import os
import re
import subprocess
import traceback

from . import Database
from . import Csb_finder
from . import myUtil
from . import Output

from multiprocessing import Pool, Manager
from typing import Dict, Optional

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



def process_missing_domains(genomeID, missing_domains_dict, candidate_hit_files_dict, gff_file, faa_file, nt_range=3500):
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
    for domain, csb_data in missing_domains_dict.items():
        if domain in candidate_hit_files_dict:
            grep_process = subprocess.Popen(
                ['grep', genomeID, candidate_hit_files_dict[domain]],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            stdout, stderr = grep_process.communicate()
            if stderr:
                print("Error for fetching candidate hits:", stderr.decode('utf-8'))
                return

            lines = stdout.decode('utf-8').split('\n')
            for line in lines:
                if line.strip():
                    columns = line.strip().split('\t')
                    key = columns[0]
                    genomeID, hit_proteinID = key.split('___',1)
                    query = columns[3]
                    hit_bitscore = int(float(columns[7]))
                    hsp_start = int(float(columns[17]))
                    hsp_end = int(float(columns[18]))
                    candidate_protein_dict[hit_proteinID] = Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore,genomeID)

            
    # Step 2: Add features to proteins
    parseGFFfile(gff_file, candidate_protein_dict)

    
    # Step 3: Filter candidates and assign to clusters
    # only take candidates that are within range of the gencluster with the pattern that misses something
    #e.g. structure {'redDsrD': [('GCF_000266945', 'GCF_000266945_6', 'GCF_000266945_000000000001', 1853338, 1859499)], 'TmcA': [('GCF_000266945', 'GCF_000266945_12', 'GCF_000266945_000000000001', 6138287, 6145798)]}}
    for domain, csb_data in missing_domains_dict.items():
        for data in csb_data:
            genome_id, clusterID, contig, cluster_start, cluster_end = data # for the clusterID we are looking for some protein that lays between these coordinates
            
            for proteinID, protein_info in candidate_protein_dict.items():
                if domain in protein_info.get_domain_set():
                    gene_start = protein_info.gene_start
                    gene_end = protein_info.gene_end

                    # Check if the protein falls within the range or 3500 nt outside of the gene cluster
                    if (cluster_start - nt_range <= gene_start <= cluster_end + nt_range) or \
                       (cluster_start - nt_range <= gene_end <= cluster_end + nt_range):
                        protein_info.clusterID = clusterID
                        insert_protein_dict[proteinID] = protein_info # Update dict for the main routine
            
    # Step 4: Add sequences to the insertion candidates
    get_protein_sequence(faa_file, insert_protein_dict)
    
    return insert_protein_dict

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
    csb_patterns, csb_names = Csb_finder.makePatternDict(options.patterns_file)
    
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
                    options.summary_hmmreport,
                    options.intermediate_hmmreport,
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
            print(f"Warning: Failed to process genome '{genome_id}' — {str(e)}")
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
        cluster_dict = Csb_finder.find_syntenicblocks(genome_id, combined_protein_dict, nucleotide_range)
        
        # Name syntenic blocks with known patterns
        cluster_dict = Csb_finder.name_syntenicblocks(pattern_dict, pattern_names, cluster_dict, min_completeness)

        # Remove intermediate hits that are not part of a named gene cluster pattern
        combined_protein_dict = remove_unassigned_intermediate_proteins(combined_protein_dict, protein_dict, cluster_dict)

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
        # covered proteinID includes the proteinIDs that are part of a recognized pattern
        proteinIDs.update(getattr(cluster, 'covered_protein_ids', set()))

    # Remove any protein in combined_protein_dict that is not in proteinIDs
    for protein_id in list(combined_protein_dict.keys()):
        if protein_id not in proteinIDs:
            del combined_protein_dict[protein_id]

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
                print(f"\tWARNING: Skipped {Filepath} due to an error - {error_message}")
                print(f"\tTraceback details:\n{traceback_details}")
                continue
                
    return protein_dict




































  


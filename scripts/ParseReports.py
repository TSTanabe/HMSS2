#!/usr/bin/python
import os
import traceback
import subprocess
from . import Database
from . import Csb_finder
from . import myUtil
from . import Output
from multiprocessing import Pool, Manager
from typing import Any, Dict, Set, Tuple, Optional, List

import re


class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have a higher score than all overlapped domains. If new domain has a lower score than any previously added domain new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein. Here the additional information of other domains is added if it does not interfere with this assumption
    
    Args:
        protein_ID - unique string as identifier
        HMM - protein type designation
        start - start coordinate of the protein type in the AA sequence
        end - end coordinate of the protein type in the AA sequence
        score - bitscore, propability of the dignated protein type
    """


    def __init__(self,proteinID,HMM,start=0,end=0,score=1,genomeID=""): 
        #Protein attributes
        self.proteinID = proteinID
        self.genomeID = genomeID
        self.protein_sequence = ""
        
        #Gene attributes
        self.gene_contig = ""
        self.gene_start = 0
        self.gene_end = 0
        self.gene_strand = "."
        self.gene_locustag = ""
        
        #Cluster attributes
        self.clusterID = "" #reziprok mit Cluster class
        self.keywords = {} #reziprok mit Cluster class
        
        self.domains = {}   # dictionary start coordinate => Domain object
        self.add_domain(HMM,start,end,score)
        

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
        

    def add_domain(self,HMM,start,end,score):
        """
        2.9.22
        Adds a domain to an existing protein. Only if there is no overlap with a 
        currently existing domain or of the overlapping domain scores higher than
        any existing overlapped domain. Overlapped domains with minor scores are deleted
        
        Args:
            HMM - protein type designation
            start - start coordinate of the protein type in the AA sequence
            end - end coordinate of the protein type in the AA sequence
            score - bitscore, propability of the dignated protein type
        Return: 
            0 - no domain was added
            1 - domain was added
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
        self.domains.update({start:Domain(HMM,start,end,score)}) # if loop complete
        
        return 1
        
        
        


class Domain:
#2.9.22
    def __init__(self,HMM,start,end,score):
        self.HMM = HMM
        self.start = int(start)
        self.end = int(end)
        self.score = int(score)
    
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



#####################################################################
######################   Parsing subroutines ########################
#####################################################################


def main_parse_summary_hmmreport(options: Any) -> None:
    """
    Parses the summary HMM report and populates the database with CSB (collinear syntenic block)
    annotations for each genome in parallel, using a writer process to commit results.

    Steps:
      1. Split queued genomes into batches based on available CPU cores.
      2. Load CSB patterns once in the main process.
      3. Build a dictionary mapping each HMM ID to its '.trusted_hits' file path.
      4. Insert all genome IDs into the database.
      5. Spawn a writer process that listens on a shared queue.
      6. Parallelize reading batches: each worker calls `process_batch(...)` to parse hits
         and enqueue results.
      7. Signal the writer to finish by sending `None` for each worker.
      8. Wait for the writer to complete before exiting.

    Args:
        options: An object with at least the following attributes:
            - queued_genomes (Set[str]): Set of genome IDs to process.
            - cores (int): Number of CPU cores to use.
            - patterns_file (str): Path to the CSB patterns file.
            - cross_check_directory (str): Directory containing '*.trusted_hits' files.
            - database_directory (str): Path to the SQLite database file.
            - faa_files (Dict[str, str]): Mapping genomeID → path to .faa file.
            - gff_files (Dict[str, str]): Mapping genomeID → path to .gff file.
            - summary_hmmreport (str): Path to the 'global_trusted_hits_summary.hmmreport'.
            - intermediate_hmmreport (str): Path to the 'global_intermediate_hits_summary.hmmreport'.
            - nucleotide_range (int): Max nucleotide distance for syntenic genes.
            - min_completeness (float): Minimum completeness fraction for CSB detection.
    """
    # 1. Split queued genomes into batches
    genome_ids: List[str] = list(options.queued_genomes)
    if options.cores < 2:
        raise ValueError("[ERROR] At least 2 cores are required (1 for parsing, 1 for writing).")
    genome_batches: List[List[str]] = split_into_batches(genome_ids, options.cores - 1)

    # 2. Load CSB patterns (only once)
    csb_patterns, csb_names = Csb_finder.makePatternDict(options.patterns_file)

    # 3. Build dict: HMM ID → path to '.trusted_hits' file
    cross_dir: str = options.cross_check_directory
    if not os.path.isdir(cross_dir):
        raise FileNotFoundError(f"[ERROR] Cross-check directory not found: {cross_dir}")

    intermediate_hit_dict: Dict[str, str] = {}
    for filename in os.listdir(cross_dir):
        if filename.endswith(".trusted_hits"):
            hmm_id = filename[:-len(".trusted_hits")]
            intermediate_hit_dict[hmm_id] = os.path.join(cross_dir, filename)

    # 4. Insert genome IDs into the database
    Database.insert_database_genomeIDs(options.database_directory, set(genome_ids))

    # 5. Spawn a writer process using a Manager Queue
    with Manager() as manager:
        data_queue = manager.Queue()

        with Pool(processes=options.cores) as pool:
            # 6a. Start the writer process
            writer_proc = pool.apply_async(process_writer, (data_queue, options))

            # 6b. Prepare arguments for each worker (process_batch)
            args_list = [
                (
                    data_queue,
                    batch,
                    options.faa_files,
                    options.gff_files,
                    options.summary_hmmreport,
                    options.intermediate_hmmreport,
                    options.nucleotide_range,
                    options.min_completeness,
                    intermediate_hit_dict,
                    csb_patterns,
                    csb_names,
                    options.synteny_support
                )
                for batch in genome_batches
            ]

            # 6c. Start reader workers in parallel
            pool.starmap(process_batch, args_list)

            # 7. Signal writer to finish by sending `None` once per worker
            for _ in range(options.cores):
                data_queue.put(None)

            # 8. Wait for the writer to complete
            writer_proc.get()

    print("[INFO] Finished parsing of search results")





def process_batch(
    data_queue: Any,
    genome_ids: List[str],
    faa_files: Dict[str, str],
    gff_files: Dict[str, str],
    hmmreport_path: str,
    intermediate_hmmreport: str,
    nucleotide_range: int,
    min_completeness: float,
    intermediate_hit_dict: Dict[str, str],
    pattern_dict: Dict[str, Any],
    pattern_names: List[str],
    synteny_support: bool
) -> None:
    """
    Processes a batch of genome IDs by extracting CSB data and enqueuing results.

    For each genome_id in `genome_ids`, this function calls `process_genome(...)`,
    catches any exceptions, logs a warning, and continues with the next genome.

    Args:
        data_queue:             Multiprocessing queue for sending parsed results to the writer.
        genome_ids:             List of genome IDs to process in this batch.
        faa_files:              Mapping from genome ID to its .faa file path.
        gff_files:              Mapping from genome ID to its .gff file path.
        hmmreport_path:         Path to the global trusted-hits summary (.hmmreport).
        intermediate_hmmreport: Path to the global intermediate-hits summary (.hmmreport).
        nucleotide_range:       Maximum nucleotide distance for syntenic genes.
        min_completeness:       Minimum completeness fraction for CSB detection.
        intermediate_hit_dict:  Mapping from HMM ID to its '.trusted_hits' file path.
        pattern_dict:           Dictionary of CSB patterns (e.g., from makePatternDict).
        pattern_names:          List of pattern names corresponding to pattern_dict keys.
    """
    for genome_id in genome_ids:
        try:
            process_genome(
                data_queue=data_queue,
                genome_id=genome_id,
                faa_path=faa_files[genome_id],
                gff_path=gff_files[genome_id],
                hmmreport_path=hmmreport_path,
                intermediate_hmmreport=intermediate_hmmreport,
                nucleotide_range=nucleotide_range,
                min_completeness=min_completeness,
                intermediate_hit_dict=intermediate_hit_dict,
                pattern_dict=pattern_dict,
                pattern_names=pattern_names,
                synteny_support=synteny_support
            )
        except Exception as exc:
            print(f"[WARNING] Failed to process genome '{genome_id}': {exc}")
            continue



def process_genome(
    data_queue: Any,
    genome_id: str,
    faa_path: str,
    gff_path: str,
    hmmreport_path: str,
    intermediate_hmmreport: str,
    nucleotide_range: int,
    min_completeness: float,
    intermediate_hit_dict: Dict[str, str],
    pattern_dict: Dict[str, Any],
    pattern_names: Set[str],
    synteny_support: bool
) -> None:
    """
    Processes a single genome: unpacks files, parses HMM hits, finds and names
    syntenic blocks, completes missing domains, and enqueues results for writing.

    Steps:
      1. Unpack .faa and .gff (if gzipped).
      2. Parse intermediate and final HMM reports to collect protein hits.
      3. Identify and name syntenic gene clusters.
      4. Remove intermediate hits not belonging to named clusters.
      5. Fetch protein sequences for combined hits.
      6. Put (protein_dict, cluster_dict) onto data_queue for the writer.

    Args:
        data_queue:            Multiprocessing queue used to send parsed results.
        genome_id:             Genome identifier (string).
        faa_path:              Path to genome’s .faa (possibly .gz).
        gff_path:              Path to genome’s .gff (possibly .gz).
        hmmreport_path:        Path to global trusted-hits summary (.hmmreport).
        intermediate_hmmreport: Path to global intermediate-hits summary (.hmmreport).
        nucleotide_range:      Maximum nucleotide distance for syntenic genes.
        min_completeness:      Minimum completeness fraction for CSB detection.
        intermediate_hit_dict: Mapping from HMM ID to its '.trusted_hits' file path.
        pattern_dict:          Dictionary of CSB patterns (from makePatternDict).
        pattern_names:         Set of pattern names corresponding to pattern_dict keys.

    Returns:
        None. Puts a tuple (combined_protein_dict, cluster_dict) onto data_queue.
    """
    try:
        # 1. Unpack FASTA and GFF (handles .gz if present)
        faa_file = myUtil.unpackgz(faa_path)
        gff_file = myUtil.unpackgz(gff_path)

        # 2. Parse intermediate HMM report for this genome
        intermediate_protein_dict = parse_bulk_hmmreport_genomize(genome_id, intermediate_hmmreport)
        parseGFFfile(gff_file, intermediate_protein_dict)

        # 3. Parse final HMM report for this genome
        protein_dict = parse_bulk_hmmreport_genomize(genome_id, hmmreport_path)
        parseGFFfile(gff_file, protein_dict)

        # 4. Synteny completion
        # Combination of intermediate and trusted hits
        if synteny_support:
            combined_protein_dict = {
                **intermediate_protein_dict,
                **protein_dict
            }
            
            # 4a. Find the syntenic clusters for combined intermediate and trusted
            cluster_dict = Csb_finder.find_syntenicblocks(
                genome_id,
                combined_protein_dict,
                nucleotide_range
            )
            
            # 4b. Name the csb from combined
            Csb_finder.name_syntenicblocks(
                pattern_dict,
                pattern_names,
                cluster_dict,
                min_completeness
            )

            # 4c. Remove intermediate hits not assigned to any named cluster
            remove_unassigned_intermediate_proteins(
                combined_protein_dict,
                protein_dict,
                intermediate_protein_dict,
                cluster_dict
            )
        
        else:
            # Do not combine the dictionaries use only the trusted and refseq hits
            combined_protein_dict = protein_dict
            
            # 4a. Find the syntenic clusters for combined intermediate and trusted
            cluster_dict = Csb_finder.find_syntenicblocks(
                genome_id,
                combined_protein_dict,
                nucleotide_range
            )
            
            # 4b. Name the csb from combined
            Csb_finder.name_syntenicblocks(
                pattern_dict,
                pattern_names,
                cluster_dict,
                min_completeness
            )
        
        
        # 6. Add protein sequences for all combined hits
        getProteinSequence(faa_file, combined_protein_dict)


        # 7. Enqueue results: combined hits + cluster assignments
        data_queue.put((combined_protein_dict, cluster_dict))

    except Exception as exc:
        print(f"[ERROR] Processing genome '{genome_id}': {exc}")
        print(traceback.format_exc())

#########################################
########   Parsing subroutines ##########
#########################################






def parseGFFfile(Filepath, protein_dict):
    """
    3.9.22
    
    Adds the general genomic features to Protein Objects in a dictionary
    
    Args:
        Filepath - GFF3 formatted file
        protein_dict - Dictionary with key proteinID and value Protein Objects
    Return:
        protein_dict (even though possibly not necessary)
    """
    locustag_pattern = re.compile(r'locus_tag=(\S*?)[\n;]')
    geneID_pattern = re.compile(r'ID=(cds-)?(\S+?)(?:[;\s]|$)')
    
    grep_pattern = "|".join(protein_dict.keys())
    
    try:
        grep_process = subprocess.Popen(['grep', '-E', grep_pattern, Filepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = grep_process.communicate()

        if stderr:
            print("Error in grep process:", stderr)
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
        error_message = f"\nError occurred: {str(e)}"
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        return protein_dict
    return protein_dict




def getProteinSequence(Filepath, protein_dict):
    """
    Fügt die Proteinsequenz zu den Proteinobjekten in einem Dictionary hinzu. Die Protein-IDs im Dictionary
    müssen mit den Headern der .faa-Datei übereinstimmen.

    Args:
        Filepath: Pfad zur Datei im FASTA-Format mit Aminosäuresequenzen
        protein_dict: Dictionary mit Protein-ID als Schlüssel und Protein-Objekten als Werte
    Return:
        protein_dict (obwohl möglicherweise nicht notwendig)
    """
    
    reader = None
    try:
        reader = open(Filepath, "r")
        sequence = ""
        header = None
        save_sequence = False  # Indikator, ob die Sequenz gespeichert werden muss
        
        for line in reader:
            line = line.strip()
            
            if line.startswith(">"):
                # Speichere die bisherige Sequenz, falls notwendig
                if header and save_sequence and sequence:
                    protein = protein_dict[header]
                    protein.protein_sequence = sequence

                # Extrahiere die Protein-ID aus dem Header (bis zum ersten Leerzeichen)
                header = line[1:].split()[0]
                
                # Prüfe, ob diese Protein-ID im Dictionary ist
                if header in protein_dict:
                    save_sequence = True  # Nur dann sammeln wir die Sequenz
                    sequence = ""  # Setze die Sequenz zurück für das neue Protein
                else:
                    save_sequence = False  # Ignoriere Sequenzen, die nicht im Dictionary sind

            elif save_sequence:
                # Füge die Sequenzzeile nur hinzu, wenn die ID im Dictionary ist
                sequence += line
        
        # Verarbeite die letzte Sequenz, falls nötig
        if header and save_sequence and sequence:
            protein = protein_dict[header]
            protein.protein_sequence = sequence
    
    except IOError:
        print(f"Error: File {Filepath} could not be opened.")
    
    finally:
        if reader is not None:
            reader.close()  # Datei explizit schließen
    
    return protein_dict

def getLocustag(locustag_pattern,string):
    match = locustag_pattern.search(string)
    if match:
        return match.group(1)
    else:
        return ""


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
    getProteinSequence(faa_file, insert_protein_dict)
    
    return insert_protein_dict

##################################################################################################
############################# Writing to database ################################################
##################################################################################################

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
    print(f"Submitted batch written to db and gene cluster file")



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
            print(f"[INFO] Processed {batch_counter} genomes ", end="\r")#
        
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

def split_into_batches(data_list, num_batches):
    if num_batches <= 0:
        raise ValueError("Number of batches must be > 0.")
    if len(data_list) == 0:
        return [[] for _ in range(num_batches)]
    
    batches = [[] for _ in range(num_batches)]
    for idx, item in enumerate(data_list):
        batches[idx % num_batches].append(item)
    return batches

def parse_bulk_hmmreport_genomize(
    genomeID: str,
    file_path: str,
    protein_dict: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Parses a bulk HMM report file for a specific genome and updates a dictionary
    of protein objects with domain hits.

    For each line in the report:
      1. Skip empty lines or comments (lines starting with '#').
      2. Split the line on tabs and confirm the first column starts with 'genome_id___'.
      3. Extract protein ID (after '___'), query name, bitscore, HSP start, and HSP end.
      4. If the protein ID already exists in protein_dict, call `add_domain(...)`.
         Otherwise, create a new Protein object and add it to protein_dict.

    Args:parse_bulk_HMMreport_genomize
        genome_id:     Genome identifier used as prefix in hit IDs (before '___').
        file_path:     Path to the HMM report file (e.g., `.hmmreport` or `.trusted_hits`).
        protein_dict:  Existing mapping from protein ID to Protein. If None, a new dict is created.

    Returns:
        A dictionary mapping protein IDs to Protein objects with updated domain info.

    """

    if protein_dict is None:
        protein_dict = {}

    result = subprocess.run(['grep', genomeID, file_path], stdout=subprocess.PIPE, text=True)
    
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
    
def remove_unassigned_intermediate_proteins(
    combined_protein_dict: Dict[str, Any],
    protein_dict: Dict[str, Any],
    intermediate_protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Removes proteins that originated from intermediate hits but are not part of any named cluster.

    A protein is kept in combined_protein_dict if it:
      - Appears in protein_dict (final hits), or
      - Appears in any cluster that has assigned keywords.

    Args:
        combined_protein_dict: Mapping from protein ID to protein object for all hits (intermediate + final).
        protein_dict:          Mapping from protein ID to protein object for final HMM hits.
        intermediate_protein_dict:
                               Mapping from protein ID to protein object for intermediate HMM hits.
        cluster_dict:          Mapping from cluster ID to cluster object. Cluster objects must support:
                                 - get_keywords() -> Set[str] or similar
                                 - get_genes() -> Set[str]

    Returns:
        The updated combined_protein_dict (with unassigned intermediate proteins removed).
    """
    # Collect all protein IDs that belong to a cluster with assigned keywords
    proteins_in_named_clusters: Set[str] = set()
    for cluster in cluster_dict.values():
        keywords = cluster.get_keywords()
        if keywords:
            proteins_in_named_clusters.update(cluster.get_genes())

    # Remove intermediate-only proteins that are not in named clusters
    for protein_id in list(intermediate_protein_dict.keys()):
        if protein_id not in protein_dict and protein_id not in proteins_in_named_clusters:
            combined_protein_dict.pop(protein_id, None)

    return combined_protein_dict
    














  


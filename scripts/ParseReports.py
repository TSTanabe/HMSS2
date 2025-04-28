#!/usr/bin/python
import os
import traceback
import subprocess
from . import Database
from . import Csb_finder
from . import myUtil
from . import Output
from multiprocessing import Pool, Manager
#import multiprocessing
#from multiprocessing import Process, Manager, Pool, Semaphore
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

def split_into_batches(data_list, num_batches):
    if num_batches <= 0:
        raise ValueError("Number of batches must be > 0.")
    if len(data_list) == 0:
        return [[] for _ in range(num_batches)]
    
    batches = [[] for _ in range(num_batches)]
    for idx, item in enumerate(data_list):
        batches[idx % num_batches].append(item)
    return batches
    
    
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

        
def remove_unassigned_intermediate_proteins(combined_protein_dict, protein_dict, intermediate_protein_dict, cluster_dict):
    """
    Removes proteins that came from the intermediate hits and are not part of any named cluster.
    
    Parameters:
        protein_dict (dict): All protein objects.
        intermediate_protein_dict (dict): Proteins originally from the intermediate HMM search.
        cluster_dict (dict): Dictionary of cluster objects after naming.
    """
    # IDs von Proteinen, die in benannten CSBs vorkommen
    proteins_in_named_clusters = set()

    for cluster in cluster_dict.values():
        if cluster.get_keywords():  # Nur Cluster mit Keywords berücksichtigen
            proteins_in_named_clusters.update(cluster.get_genes())

    # Jetzt alle Intermediate-Proteine durchgehen
    for protein_id in list(intermediate_protein_dict.keys()):
        if protein_id not in protein_dict and protein_id not in proteins_in_named_clusters:
            del combined_protein_dict[protein_id]

    return combined_protein_dict    
    
def process_genome(
    data_queue, genome_id, faa_path, gff_path, hmmreport_path, intermediate_hmmreport,
    nucleotide_range, min_completeness, intermediate_hit_dict,
    pattern_dict, pattern_names
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
        Csb_finder.name_syntenicblocks(pattern_dict, pattern_names, cluster_dict, min_completeness)

        # Add proteinID with named syntenic blocks to the proteinID
        remove_unassigned_intermediate_proteins(combined_protein_dict, protein_dict, intermediate_protein_dict, cluster_dict)

        # Add the protein sequences
        getProteinSequence(faa_file, combined_protein_dict)

        # Synteny completion mechanism for named gene clusters
        missing = Csb_finder.extract_missing_domain_targets(cluster_dict)
        new_proteins = process_missing_domains(genome_id, missing,intermediate_hit_dict, gff_file, faa_file, nucleotide_range)
        for pid, prot in new_proteins.items():
            if pid not in protein_dict:
                protein_dict[pid] = prot
                cluster_dict[prot.clusterID].add_gene(pid, prot.get_domains(), prot.gene_start, prot.gene_end)

        data_queue.put((combined_protein_dict, cluster_dict))

    except Exception as e:

        print(f"Error: {genome_id} -> {e}")
        print(traceback.format_exc())



def process_batch(
    data_queue, genome_ids, faa_files, gff_files, hmmreport_path, intermediate_hmmreport,
    nucleotide_range, min_completeness, intermediate_hit_dict,
    pattern_dict, pattern_names
):
    for genome_id in genome_ids:
        try:
            process_genome(
                data_queue, genome_id,
                faa_files[genome_id], gff_files[genome_id], hmmreport_path, intermediate_hmmreport,
                nucleotide_range, min_completeness, intermediate_hit_dict,
                pattern_dict, pattern_names
            )
        except Exception as e:
            print(f"Warning: Failed to process genome '{genome_id}' — {str(e)}")
            continue


def main_parse_summary_hmmreport(options):


    genome_ids = list(options.queued_genomes)
    genomeID_batches = split_into_batches(genome_ids, options.cores - 1)

    # Lade Patterns nur 1x im Hauptprozess
    csb_patterns, csb_names = Csb_finder.makePatternDict(options.patterns_file)
    
    # Load intermediate hit file
    intermediate_hit_dict = {os.path.splitext(f)[0]: os.path.join(options.Cross_check_directory, f) for f in os.listdir(options.Cross_check_directory) if f.endswith('.trusted_hits')}

    
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
                    intermediate_hit_dict,
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
    print("Finished parsing of search results")



































  


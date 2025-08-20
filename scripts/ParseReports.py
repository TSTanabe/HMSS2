#!/usr/bin/python
import os
import re
import subprocess
import traceback
import itertools
import pprint

from . import Search
from . import Database
from . import Csb_finder
from . import myUtil
from . import Output

import numpy as np

from multiprocessing import Pool, Manager
from typing import Dict, Optional, Any, List, Set
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
        bsr: float = 1.0,
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
        self.selection_comment: Set[str] = set()  # Trusted cutoff Flag or cooccurrence
        self.alternative_hit: str = ""

    ##### Getter ####

    def get_domains(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_HMM())
        return "-".join(listing)

    def get_domains_dict(self):
        # return dict
        return self.domains

    def get_domain_listing(self):
        # return list
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
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(
                f"{self.domains[key].get_start()}:{self.domains[key].get_end()}"
            )
        return "-".join(listing)

    def get_domain_scores(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_score()}")
        return "-".join(listing)

    def get_domain_count(self):
        return len(self.domains)

    def get_protein_string(self):
        # 3.9.22 representation of the whole protein in one line
        a = self.proteinID
        b = self.get_domains()
        c = self.get_domain_scores()
        d = self.gene_contig
        e = self.gene_start
        f = self.gene_end
        # d = self.get_protein_sequence()
        string = f"{a} {b} {c} {d} {e} {f}"
        return string

    def get_protein_list(self):
        # 3.9.22 representation of the whole protein in one line
        listing = []
        listing.append(self.proteinID)
        listing.append(self.get_domains())
        listing.append(str(self.get_domain_scores()))
        listing.append(str(self.get_domain_coordinates()))
        listing.append(self.gene_contig)
        listing.append(str(self.gene_start))
        listing.append(str(self.gene_end))
        listing.append(self.gene_strand)
        listing.append(self.gene_locustag)
        # d = self.get_protein_sequence()
        # string = f"{a} {b} {c} {d} {e} {f}"
        return listing

    def get_sequence(self):
        return str(self.protein_sequence)

    def get_length(self):
        return len(self.protein_sequence)

    def get_selection_comment_csv(self, sep: str = ",") -> str:
        """
        Return the selection_comment set as a sorted, comma-separated string.
        Sorting ensures stable output for logs/tables.
        """
        return sep.join(sorted(self.selection_comment))

    ##### Setter #####

    def check_domain_overlap(self, new_start, new_end, current_start, current_end):
        # 2.9.22

        if (current_start <= new_start and new_start <= current_end) or (
            current_start <= new_end and new_end <= current_end
        ):
            # start oder endpunkt innerhalb er grenzen
            return 1

        elif (new_start <= current_start and current_end <= new_end) or (
            current_start <= new_start and new_end <= current_end
        ):
            # start und end innerhalb der grenzen oder alte domäne innerhalb der neuen
            return 1

        elif (new_start <= current_start and new_end <= current_start) or (
            new_start >= current_end and new_end >= current_end
        ):
            # start und end kleiner als self start oder start und end größer als self end dann
            # domäne außerhalb der alten domäne und adden egal welcher score
            return 0
        return None

    def remove_selection_comment(self) -> None:
        self.selection_comment = set()

    def add_selection_comment(self, comment: str, sep: str = ",") -> None:
        """
        Add one or multiple comment tokens to this protein.
        - Accepts a single token ("Tc") or a CSV string ("Tc,Coo").
        """
        if not comment:
            return
        tokens = [t.strip() for t in str(comment).split(sep) if t.strip()]
        self.selection_comment.update(tokens)

    def add_domain(
        self,
        HMM: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0,
    ) -> int:
        """
        Adds a domain to the protein only if it does not overlap
        with a higher-scoring existing domain. If overlap exists with lower-scoring
        domain, that domain is removed.

        Returns:
            int: 1 if domain added, 0 if not added.
        """

        del_domains = []  # start coordinates/keys of domains to be replace
        for domain in self.domains.values():
            if self.check_domain_overlap(
                start, end, domain.get_start(), domain.get_end()
            ):
                if domain.get_score() < score:
                    del_domains.append(domain.get_start())
                else:
                    return 0

        for key in del_domains:
            self.domains.pop(key)
        self.domains.update(
            {start: Domain(HMM, start, end, score, ident, bsr)}
        )  # if loop complete

        return 1


class Domain:
    # 2.9.22
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

    def __init__(
        self,
        HMM: str,
        start: int,
        end: int,
        score: float,
        ident: int = 1,
        bsr: float = 1.0,
    ):
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
            return (
                self.HMM == other.HMM
                and self.start == other.start
                and self.end == other.end
                and self.score == other.score
            )
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


def parseGFFfile(filepath: str, protein_dict: Dict[str, Protein]) -> Dict[str, Protein]:
    """
    3.9.22
    Adds GFF attributes to each Protein object in the dictionary.

    Args:
        filepath (str): Path to GFF3 file.
        protein_dict (dict): {proteinID: Protein object}

    Returns:
        dict: Updated protein_dict.
    """
    locustag_pattern = re.compile(r"locus_tag=(\S*?)(?:[;\s]|$)")
    geneID_pattern = re.compile(r"ID=(cds-)?(\S+?)(?:[;\s]|$)")

    grep_pattern = "|".join(protein_dict.keys())
    try:
        grep_process = subprocess.Popen(
            ["grep", "-E", grep_pattern, filepath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = grep_process.communicate()

        if stderr:
            logger.error("Grep process:", stderr)
            return protein_dict

        for line in stdout.decode("utf-8").split("\n"):
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
            file.write(clusterID + "\t" + "\t".join(domains) + "\n")


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
            print(f"Processed {batch_counter} genomes ", end="\r")  #

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
                filepath = os.path.join(
                    options.fasta_initial_hit_directory,
                    str(genomeID) + ".hit_table_txt",
                )
                Output.output_genome_report(filepath, protein_dict, cluster_dict)

        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            submit_batches(protein_batch, cluster_batch, options)
            protein_batch = {}
            cluster_batch = {}
            batch_counter = 0

    # Submit the remaining and file reports
    if protein_batch or cluster_batch:
        submit_batches(protein_batch, cluster_batch, options)
        if options.individual_reports:
            if protein_dict:
                first_protein_key = next(iter(protein_dict))  # Get the first key
                genomeID = protein_dict[first_protein_key].genomeID
                filepath = os.path.join(
                    options.fasta_initial_hit_directory,
                    str(genomeID) + ".hit_table_txt",
                )

                Output.output_genome_report(filepath, protein_dict, cluster_dict)
    logger.info(f"Processed {batch_counter} genomes")
    return


#########################################################################################
################ Processing routines for parsing genome hits ############################
#########################################################################################


def main_parse_summary_hmmreport(options):
    genome_ids = list(options.queued_genomes)
    genomeID_batches = split_into_batches(genome_ids, options.cores - 1)

    # Lade Patterns nur 1x im Hauptprozess
    csb_patterns, csb_names = Csb_finder.make_pattern_dict(options.patterns_file)
    cooccurrence_pattern, cooccurence_names = Csb_finder.make_pattern_dict(
        options.cooccurrence_file
    )
    threshold_dict = Search.make_threshold_dict(
        options.score_threshold_file, 3, options.thrs_score
    )
    exclusion_singletons = parse_exclusion_singletons(options.exclusion_singletons)

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
                    csb_names,
                    cooccurrence_pattern,
                    exclusion_singletons,
                    threshold_dict,
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


def parse_exclusion_singletons(file_path: str) -> set:
    """
    Parse a tab-separated file and collect all words/tokens into a single set,
    regardless of line count or number of tokens per line.

    Args:
        file_path (str): Path to the Exclusion_singletons file.

    Returns:
        Set[str]: All unique words found in the file.
    """
    singleton_set = set()
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            # Split line by tab and extend set
            tokens = line.strip().split("\t")
            singleton_set.update(token for token in tokens if token)
    return singleton_set


def process_batch(
    data_queue,
    genome_ids,
    faa_files,
    gff_files,
    trusted_hmmreport_path,
    intermediate_hmmreport,
    nucleotide_range,
    min_completeness,
    pattern_dict,
    pattern_names,
    cooccurrence_pattern,
    exclusion_singletons,
    threshold_dict,
):
    """
    Runs process_genome for each genome in the batch with all necessary arguments.
    """
    for genome_id in genome_ids:
        try:
            process_genome(
                data_queue=data_queue,
                genome_id=genome_id,
                faa_path=faa_files[genome_id],
                gff_path=gff_files[genome_id],
                trusted_hmmreport_path=trusted_hmmreport_path,
                intermediate_hmmreport=intermediate_hmmreport,
                nucleotide_range=nucleotide_range,
                min_completeness=min_completeness,
                pattern_dict=pattern_dict,
                pattern_names=pattern_names,
                cooccurrence_pattern=cooccurrence_pattern,
                exclusion_singletons=exclusion_singletons,
                threshold_dict=threshold_dict,
            )
        except Exception as e:
            logger.warning(f"Failed to process genome '{genome_id}' — {str(e)}")
            continue


def process_genome(
    data_queue: Any,
    genome_id: str,
    faa_path: str,
    gff_path: str,
    trusted_hmmreport_path: str,
    intermediate_hmmreport: str,
    nucleotide_range: int,
    min_completeness: float,
    pattern_dict: Dict[str, list],
    pattern_names: Dict[str, Any],
    cooccurrence_pattern: Dict[str, list],
    exclusion_singletons: Set[str],
    threshold_dict: Dict[str, float],
) -> None:
    """
    Main pipeline to process one genome:
    - Parses protein and cluster data
    - Annotates clusters
    - Removes singleton proteins with specified domains outside clusters
    - Removes unassigned intermediate hits, keeping trusted proteins
    - Attaches protein sequences
    - Returns (combined_protein_dict, cluster_dict) via queue

    Args:
        data_queue: Multiprocessing queue for result transport.
        genome_id: ID of the genome.
        faa_path, gff_path: Paths to input files (can be .gz).
        trusted_hmmreport_path, intermediate_hmmreport: HMM report files.
        nucleotide_range: Nucleotide window for cluster detection.
        min_completeness: Minimum completeness for cluster pattern assignment.
        pattern_dict: Patterns for cluster annotation.
        pattern_names: Names/labels for patterns.
        exclusion_singletons: Domains for singleton exclusion.
        threshold_dict: Score cutoffs per domain.
    """
    try:
        faa_file = myUtil.unpackgz(faa_path)
        gff_file = myUtil.unpackgz(gff_path)

        # Intermediate protein hits
        intermediate_protein_dict = parse_bulk_HMMreport_genomize(
            genome_id, intermediate_hmmreport
        )
        parseGFFfile(gff_file, intermediate_protein_dict)

        # Primary protein hits
        trusted_protein_dict = parse_bulk_HMMreport_genomize(
            genome_id, trusted_hmmreport_path
        )
        parseGFFfile(gff_file, trusted_protein_dict)

        # Combine protein dictionaries
        combined_protein_dict = {**intermediate_protein_dict, **trusted_protein_dict}

        # Detect and annotate syntenic gene clusters
        cluster_dict = Csb_finder.find_syntenic_blocks(
            genome_id, combined_protein_dict, nucleotide_range
        )
        cluster_dict = Csb_finder.name_syntenic_blocks(
            pattern_dict, pattern_names, cluster_dict, min_completeness
        )

        # Attach the reason for selection to protein objects trusted cutoff/reference sequence
        trusted_cutoff_protein_ids = set(trusted_protein_dict.keys())
        combined_protein_dict = add_selection_comment_to_many_proteins(
            combined_protein_dict, trusted_cutoff_protein_ids, "Tc"
        )

        # Enhance cluster completeness if needed
        # alters the combined_protein_dict
        enhance_syntenic_block_completeness(
            cluster_dict,
            combined_protein_dict,
            intermediate_protein_dict,
            intermediate_hmmreport,
            pattern_dict,
        )

        # Collect trusted protein IDs: those with complete pathways and those in the main protein dict
        singletons_with_complete_pathway_set = enhance_pathway_completeness(
            combined_protein_dict, cooccurrence_pattern, threshold_dict
        )

        trusted_protein_ids = trusted_cutoff_protein_ids.union(
            singletons_with_complete_pathway_set
        )

        # Remove unassigned intermediate proteins, but keep trusted ones
        combined_protein_dict = remove_unassigned_intermediate_proteins(
            combined_protein_dict, trusted_protein_ids, cluster_dict
        )

        # Remove genes that should not occur as singletons
        # alters the combined_protein_dict but ignores singletons that complete pathway
        remove_exclusion_singletons(
            combined_protein_dict,
            cluster_dict,
            exclusion_singletons,
            singletons_with_complete_pathway_set,
        )

        # Attach the reason for selection to protein objects complete pathway
        combined_protein_dict = add_selection_comment_to_many_proteins(
            combined_protein_dict, singletons_with_complete_pathway_set, "Coo"
        )

        # Attach protein sequences
        get_protein_sequence(faa_file, combined_protein_dict)

        data_queue.put((combined_protein_dict, cluster_dict))

    except Exception as e:
        logger.error(f"Error: {genome_id} -> {e}")
        logger.error(traceback.format_exc())


def enhance_pathway_completeness(
    protein_dict: Dict[str, "Protein"],
    pattern_dict: Dict[str, List[str]],
    threshold_dict: Dict[str, float],
) -> Set[str]:
    """
    Identify proteins that contribute to fully complete pathway patterns.

    Logic:
    For each pattern in `pattern_dict`, this function checks:
      1. Is there at least one protein in the genome for *each* required domain
         that scores above or equal to the threshold for that domain?
         - (A pattern is considered "complete" only if this is true for *every* domain.)
      2. If so, gather *all* protein IDs that have any required domain with score >= its threshold
         (i.e., for each domain in the pattern, collect all proteins that fulfill the score criterion).

    If any required domain for a pattern has no matching protein above threshold,
    the pattern is ignored and does not contribute to the final set.

    Inputs:
    protein_dict : Dict[str, Protein]
        Mapping of protein IDs to Protein objects.

    pattern_dict : Dict[str, List[str]]
        Mapping of pattern names to lists of required domain names.
                "PatternB": ["X", "Y"]

    threshold_dict : Dict[str, float]
                "A": 40.0,

    Outputs:
    Set[str]
        Set of all protein IDs that fulfill at least one required domain
        (with score >= threshold) for any *fully complete* pattern.
        Each protein ID is included at most once (set semantics).

    """
    found_protein_ids: Set[str] = set()
    domain_to_protein: Dict[str, List[tuple[str, float]]] = {}

    # Build mapping from domain name to all (protein_id, score) tuples in the genome
    for protein_id, protein in protein_dict.items():
        for domain in protein.get_domain_listing():
            domain_name = domain.get_HMM()
            score = domain.get_score()
            domain_to_protein.setdefault(domain_name, []).append((protein_id, score))

    for required_domains in pattern_dict.values():
        domain_hits: Dict[str, Set[str]] = {}
        all_domains_above = True
        for domain in required_domains:
            threshold = threshold_dict.get(domain, 0)
            hits = {
                protein_id
                for protein_id, score in domain_to_protein.get(domain, [])
                if score >= threshold
            }
            if not hits:
                all_domains_above = False
                break
            domain_hits[domain] = hits

        if all_domains_above:
            # Format readable block
            msg = f"Pattern complete: {', '.join(domain_hits.keys())}\n"
            for domain in domain_hits:
                proteins = ", ".join(sorted(domain_hits[domain]))
                msg += f"  {domain}: {proteins}\n"
            logger.debug(msg.rstrip())
            for hits in domain_hits.values():
                found_protein_ids.update(hits)

    return found_protein_ids


def remove_exclusion_singletons(
    combined_protein_dict: Dict[str, "Protein"],
    cluster_dict: Dict[str, "Cluster"],
    exclusion_singletons: Set[str],
    trusted_protein_ids: Set[str],
) -> None:
    """
    Removes from `combined_protein_dict` all Protein objects that:
      - are NOT present in any cluster in `cluster_dict`
      - AND have at least one domain whose name is in `exclusion_singletons`
      - AND are NOT present in `trusted_protein_ids`

    This function modifies `combined_protein_dict` in place.

    Args:
        combined_protein_dict: Mapping of proteinID to Protein object.
        cluster_dict: Mapping of clusterID to Cluster object.
        exclusion_singletons: Set of domain names for exclusion.
        trusted_protein_ids: Set of proteinIDs to protect from removal.
    """
    proteins_in_clusters: Set[str] = set()
    for cluster in cluster_dict.values():
        proteins_in_clusters.update(cluster.genes)

    to_remove: Set[str] = set()
    for protein_id, protein in combined_protein_dict.items():
        if (
            protein_id not in proteins_in_clusters
            and protein_id not in trusted_protein_ids
        ):
            for domain in protein.get_domain_listing():
                if domain.get_HMM() in exclusion_singletons:
                    to_remove.add(protein_id)
                    break

    for protein_id in to_remove:
        protein_obj = combined_protein_dict[protein_id]
        domains = protein_obj.get_domains()  # Gibt z.B. 'HMM_A-HMM_B-HMM_C' zurück
        logger.debug(f"Removed singleton {protein_id} with domains: {domains}")
        del combined_protein_dict[protein_id]


def remove_unassigned_intermediate_proteins(
    combined_protein_dict: Dict[str, Any], proteinIDs: set, cluster_dict: Dict[str, Any]
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
    for cluster in cluster_dict.values():
        proteinIDs.update(getattr(cluster, "covered_protein_ids", set()))

    # Remove any protein in combined_protein_dict that is not in proteinIDs
    for protein_id in list(combined_protein_dict.keys()):
        if protein_id not in proteinIDs:
            del combined_protein_dict[protein_id]

    # Remove unassigned genes and update clusters accordingly
    for clusterID in list(cluster_dict.keys()):  # list() so we can delete inside loop
        cluster = cluster_dict[clusterID]
        # Remove genes from cluster.genes and cluster.types that are not in combined_protein_dict
        if hasattr(cluster, "genes"):
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


def add_selection_comment_to_many_proteins(
    combined_protein_dict: Dict[str, "Protein"],
    protein_ids_set: Set[str],
    comment: str,
) -> Dict[str, "Protein"]:
    """
    Add `comment` once to `selection_comment` for all proteins whose ID is in `protein_ids_set`.
    - Creates `selection_comment` if missing.
    - Preserves order and avoids duplicates.
    - Skips empty/whitespace comments.
    """
    comment = (comment or "").strip()
    if not comment:
        return combined_protein_dict  # nothing to add

    # Iterate over the IDs we actually want to update
    for protein_id in protein_ids_set:
        protein = combined_protein_dict.get(protein_id)
        protein.add_selection_comment(comment)

    return combined_protein_dict


def parse_bulk_HMMreport_genomize(genomeID, Filepath, protein_dict=None):
    """ """
    if protein_dict is None:
        protein_dict = {}

    result = subprocess.run(
        ["grep", genomeID, Filepath], stdout=subprocess.PIPE, text=True
    )

    lines = result.stdout.splitlines()  # Split output into lines

    for line in lines:
        columns = line.split("\t")  # Assuming columns are space-separated
        if columns:
            try:
                key = columns[0]
                genomeID, hit_proteinID = key.split("___", 1)
                query = columns[3].split("_")[-1]
                hit_bitscore = int(float(columns[7]))
                hsp_start = int(float(columns[17]))
                hsp_end = int(float(columns[18]))

                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    protein.add_domain(query, hsp_start, hsp_end, hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = Protein(
                        hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID
                    )
            except Exception as e:
                error_message = f"\nError occurred: {str(e)}"
                traceback_details = traceback.format_exc()
                logger.error(f"Skipped {Filepath} due to an error - {error_message}")
                logger.error(f"Traceback details:\n{traceback_details}")
                continue

    return protein_dict


def enhance_syntenic_block_completeness(
    cluster_dict,
    combined_protein_dict,
    intermediate_protein_dict,
    intermediate_hmmreport,
    pattern_dict,
    min_completeness=0.5,
):
    """
    Enhance syntenic block completeness by swapping additional protein domains with missing ones if possible.
    for each cluster that has no directly matching pattern from the given patterns
    it is tested if a possible conversion of protein types to alternative ones with lower hitscore could reach
    a better completion
    """

    # Go through all clusters
    for cluster in cluster_dict.values():
        # Skip clusters that have a perfectly matching pattern
        if any(
            kw.get_completeness() == 1 and not kw.get_additional_domains()
            for kw in cluster.get_keywords()
        ):
            continue

        genes = cluster.genes
        types = cluster.types

        # Dictionary to hold the hmmreport lines
        alternative_protein_type_dict = {}

        # Possible executions to increase the completeness for each keyword
        possible_optimized_executable_transitions = {}

        # Filter keyword to longest keyword with minimal equal number of missing and additionals
        best_keywords = select_balanced_best_keywords(
            cluster.get_keywords(), pattern_dict
        )
        if not best_keywords:
            continue

        for keyword in best_keywords:
            completeness = keyword.get_completeness()
            missing_domains = keyword.get_missing_domains()
            additional_domains = keyword.get_additional_domains()
            transition_dict = defaultdict(set)
            # Structure: alternative_protein_type_dict[(hit_proteinID, query)] = (hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID)
            # Structure: transition_dict[query] => ((hit_proteinID, difference), (hit_proteinID, difference))
            # print(f"\nNew gene cluster {cluster.clusterID}")
            # print(types)
            # print(f"Processing {keyword.keyword} {completeness}")
            # print(f"Missing {missing_domains}")
            # print(f"Additional {additional_domains}")

            if (
                completeness >= min(0.5, min_completeness)
                and completeness < 1
                and additional_domains
            ):
                # For every "additional" domain: get proteinID and possible transitions to "missing"
                transitions = []
                for index, current_domain in enumerate(types):
                    # Parse domtblout to get all potential domains and their scores
                    proteinID = genes[index]
                    protein = combined_protein_dict[proteinID]

                    # Get the proteinID of the additional domains
                    if (
                        current_domain in additional_domains
                        and genes[index] in intermediate_protein_dict
                        and not "Tc" in protein.selection_comment
                    ):
                        # Updates the alternative_protein_type_dict and transition_dict
                        find_possible_transitions(
                            proteinID,
                            current_domain,
                            missing_domains,
                            protein,
                            intermediate_hmmreport,
                            alternative_protein_type_dict,
                            transition_dict,
                        )

                # Optimize the transitions by minimizing the bitscore changes and number of transitions to reach the missing domains
                pattern_length = len(pattern_dict[keyword.keyword_id])

                chosen_transitions, posterior_completeness, total_score_diff = (
                    get_optimal_transitions(
                        transition_dict, missing_domains, completeness, pattern_length
                    )
                )
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
                    possible_optimized_executable_transitions[
                        (pattern_length, posterior_completeness, total_score_diff)
                    ] = chosen_transitions
        # Keyword loop finished

        # Now from all possible keyword completions find the optimum dictionary key
        # (transitions, bitscore, posterior completeness) => [(proteinID to domain), (proteinID to domain), (proteinID to domain)]
        if possible_optimized_executable_transitions:
            # Sort by: pattern_length (desc), posterior_completeness (desc), total_score_diff (asc)
            best_key = max(
                possible_optimized_executable_transitions.keys(),
                key=lambda x: (
                    x[1],
                    x[0],
                    -x[2],
                ),  # pattern_length, posterior_completeness, -score_diff
            )
            best_transitions = possible_optimized_executable_transitions[best_key]

            logger.debug(
                f"For clusterID {cluster.clusterID} with pattern length, completeness and score difference {best_key} following conversion is done"
            )
            logger.debug(best_transitions)

            # For the best transition alter the proteins domain information
            for proteinID, to_domain, _ in best_transitions:
                if proteinID in combined_protein_dict:
                    protein = combined_protein_dict[proteinID]

                    # Save the original hit as comment
                    original_domains = protein.get_domains()
                    protein.add_selection_comment("Syc")
                    protein.alternative_hit = original_domains

                    # Add the alternative lower hit for synteny completion
                    hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID = (
                        alternative_protein_type_dict.get((proteinID, to_domain))
                    )  # (hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID)
                    protein.domains.clear()
                    protein.add_domain(query, hsp_start, hsp_end, hit_bitscore)

    return combined_protein_dict


def select_balanced_best_keywords(keywords, pattern_dict):
    """
    Filtere und wähle die besten Keywords aus:
    1. Nur mit ausgeglichenem Verhältnis missing/additional (beide > 0, gleich groß)
    2. Nur mit minimaler Anzahl missing/additional
    3. Nur mit maximaler Pattern-Länge (ggf. mehrere)
    Gibt Liste der besten Keyword-Objekte zurück.
    """
    keyword_infos = []
    for keyword in keywords:
        missing_domains = keyword.get_missing_domains()
        additional_domains = keyword.get_additional_domains()
        pattern_length = len(pattern_dict[keyword.keyword_id])
        keyword_infos.append(
            {
                "keyword": keyword,
                "n_missing": len(missing_domains),
                "n_additional": len(additional_domains),
                "pattern_length": pattern_length,
            }
        )

    # Filter: ausgeglichen und mindestens 1 fehlend
    balanced = [
        info
        for info in keyword_infos
        if info["n_missing"] == info["n_additional"] and info["n_missing"] > 0
    ]
    if not balanced:
        return []

    # Minimal missing/additional
    min_missing = min(info["n_missing"] for info in balanced)
    minimal = [info for info in balanced if info["n_missing"] == min_missing]

    # Längstes Pattern
    max_len = max(info["pattern_length"] for info in minimal)
    best_keywords = [
        info["keyword"] for info in minimal if info["pattern_length"] == max_len
    ]
    return best_keywords


def find_possible_transitions(
    proteinID,
    current_domain,
    missing_domains,
    protein,
    hmmreport,
    alternative_protein_type_dict,
    transition_dict,
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
    result = subprocess.run(
        ["grep", proteinID, hmmreport], stdout=subprocess.PIPE, text=True
    )
    lines = result.stdout.splitlines()

    for line in lines:
        columns = line.split("\t")
        if columns and len(columns) > 18:
            try:
                key = columns[0]
                genomeID, hit_proteinID = key.split("___", 1)
                query = columns[3].split("_")[-1]
                hit_bitscore = float(columns[7])
                hsp_start = int(float(columns[17]))
                hsp_end = int(float(columns[18]))

                if query in missing_domains:
                    alternative_protein_type_dict[(hit_proteinID, query)] = (
                        hit_proteinID,
                        query,
                        hsp_start,
                        hsp_end,
                        hit_bitscore,
                        genomeID,
                    )
                    current_domain_score = next(
                        (
                            domain.get_score()
                            for domain in protein.domains.values()
                            if domain.get_HMM() == current_domain
                        ),
                        0,
                    )

                    difference = abs(current_domain_score - hit_bitscore)
                    transition_dict[query].add((hit_proteinID, difference))

            except Exception as e:
                logger.warn(f"Skipped line in {hmmreport} due to an error - {str(e)}")
                continue

    return alternative_protein_type_dict, transition_dict


def get_optimal_transitions(
    transition_dict,
    missing_domains,
    initial_completeness=0.0,
    total_domains=None,
    logger=None,
):
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
    orig_row_ind, orig_col_ind = solve_assignment(cost_matrix, logger=logger)

    if orig_row_ind is None:
        # Matching impossible, keep the cluster as assigned
        return [], initial_completeness, 0.0

    chosen_transitions = []
    total_score_diff = 0
    fulfilled_domains = set()
    used_proteins = set()
    missing_domains = list(missing_domains)

    for i, j in zip(orig_row_ind, orig_col_ind):
        cost = cost_matrix[i, j]
        if np.isfinite(cost):
            proteinID = proteins[i]
            domain = missing_domains[j]
            chosen_transitions.append((proteinID, domain, cost))
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


def solve_assignment(cost_matrix, logger=None, cluster_id=None):
    """
    Robustly solves an assignment problem using the Hungarian algorithm,
    even if some rows are infeasible (only np.inf).
    If no assignment is possible, returns None.
    Debug-Ausgabe: Gibt infeasible Matrix bei Problemen als logger.debug aus.
    """
    orig_matrix = cost_matrix.copy()
    n_rows, n_cols = cost_matrix.shape

    valid_rows = ~np.all(np.isinf(cost_matrix), axis=1)
    valid_cols = ~np.all(np.isinf(cost_matrix), axis=0)
    reduced_matrix = cost_matrix[np.ix_(valid_rows, valid_cols)]

    def matrix_str(matrix):
        with np.printoptions(
            precision=2, suppress=True, linewidth=120, nanstr="nan", infstr="inf"
        ):
            return "\n" + "\n".join(" ".join(f"{x:7}" for x in row) for row in matrix)

    # Check auf infeasibility
    if (
        reduced_matrix.size == 0
        or np.any(np.all(np.isinf(reduced_matrix), axis=1))
        or np.any(np.all(np.isinf(reduced_matrix), axis=0))
    ):
        if logger:
            logger.debug(
                f"[Assignment] Infeasible cost matrix for cluster {cluster_id or ''}:\n{matrix_str(cost_matrix)}"
            )
            logger.error(
                f"Assignment failed for cluster {cluster_id or ''}: cost matrix infeasible after row/col removal. Gencluster bleibt unverändert."
            )
        return None, None

    try:
        row_ind, col_ind = linear_sum_assignment(reduced_matrix)
        orig_row_ind = np.where(valid_rows)[0][row_ind]
        orig_col_ind = np.where(valid_cols)[0][col_ind]
        return orig_row_ind, orig_col_ind
    except Exception as e:
        if logger:
            logger.debug(
                f"[Assignment] Exception on cost matrix for cluster {cluster_id or ''}:\n{matrix_str(cost_matrix)}"
            )
            logger.error(
                f"Assignment error for cluster {cluster_id or ''}: {str(e)}. Gencluster bleibt unverändert."
            )
        return None, None

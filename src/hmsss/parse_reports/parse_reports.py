#!/usr/bin/python
import os
import re
import subprocess
import traceback

from multiprocessing import Pool, Manager
from typing import Dict, Any, Set

from contextlib import contextmanager
from time import perf_counter

from hmsss.cross_check import search_cross_reference, report_db
from hmsss.parse_reports import pattern_completion_synteny, pattern_completion_pathway, csb_finder, csb_trie_algorithm
from hmsss.parse_reports.csb_finder import Cluster
from hmsss.parse_reports.csb_trie_algorithm import TrieIndex
from hmsss.core.logging import get_logger
from hmsss.db import database

logger = get_logger(__name__)


@contextmanager
def tick(label: str):
    t0 = perf_counter()
    try:
        yield
    finally:
        dt = perf_counter() - t0
        logger.info("[TIMER] %s took %.3f s", label, dt)


class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned
    to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores
    and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is
    checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have
    a higher score than all overlapped domains. If new domain has a lower score than any previously added domain
    new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein.
    Here the additional information of other domains is added if it does not interfere with this assumption

    Organizes protein domains and related attributes for a protein.
    New domains are only added if they do not overlap with higher-scoring existing domains.
    All coordinates, scores, and HMM names are accessible as dash-separated strings.

    Args:
        protein_id (str): Unique identifier.
        hmm (str): Domain name.
        start (int): Start coordinate.
        end (int): End coordinate.
        score (float): Domain bitscore.
        genome_id (str): Genome identifier (optional).
        ident (int): Percent identity (default: 25).
        bsr (float): Blast score ratio (default: 1.0).

    Example:
        p = Protein("Prot1", "HMM_A", 10, 120, 40)
        p.add_domain("HMM_B", 130, 200, 30)
    """

    def __init__(
        self,
        protein_id: str,
        hmm: str,
        start: int = 0,
        end: int = 0,
        score: float = 1,
        genome_id: str = "",
        ident: int = 25,
        bsr: float = 1.0,
    ):
        self.proteinID: str = protein_id
        self.genomeID: str = genome_id
        self.protein_sequence: str = ""
        self.gene_contig: str = ""
        self.gene_start: int = 0
        self.gene_end: int = 0
        self.gene_strand: str = "."
        self.gene_locustag: str = ""
        self.clusterID: str = ""
        self.keywords: Dict = {}
        self.domains: Dict[int, Domain] = {}  # start coordinate → Domain object
        self.deleted_domains: Dict[
            str, Domain
        ] = {}  # Domain objects that have been removed
        self.add_domain(hmm, start, end, score, ident, bsr)
        self.selection_comment: Set[str] = set()  # Trusted cutoff Flag or cooccurrence
        self.alternative_hit: str = ""
        self.valid_hit = True

    ##### Getter ####

    def get_domains(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_domain())
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
            domains.add(v.get_domain())
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

    def get_protein_list(self):
        # 3.9.22 representation of the whole protein in one line
        listing = [self.proteinID, self.get_domains(), str(self.get_domain_scores()),
                   str(self.get_domain_coordinates()), self.gene_contig, str(self.gene_start), str(self.gene_end),
                   self.gene_strand, self.gene_locustag]
        # d = self.get_protein_sequence()
        # string = f"{a} {b} {c} {d} {e} {f}"
        return listing

    def get_sequence(self):
        return str(self.protein_sequence)

    def get_selection_comment_csv(self, sep: str = ",") -> str:
        """
        Return the selection_comment as a sep-separated string.
        Accepts both set[str] (preferred) and str (fallback).
        """
        val = self.selection_comment
        if not val:
            return ""
        if isinstance(val, str):
            return val  # bereits CSV-String
        # erwarteter Fall: Menge/Tokens
        return sep.join(sorted(val))

    ##### Setter #####

    @staticmethod
    def check_domain_overlap(new_start, new_end, current_start, current_end):
        # 2.9.22

        if (current_start <= new_start <= current_end) or (
                current_start <= new_end <= current_end
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
        hmm: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0,
        *,
        force: bool = False,
    ) -> int:
        """
        Adds a domain to the protein.

        Standardverhalten (force=False):
          - Wenn Überlappung mit existierender Domäne vorliegt:
              * Entferne überlappende Domänen mit geringerem Score
              * Brich ab (return 0), wenn eine überlappende Domäne >= Score hat.
          - Ansonsten einfügen (return 1).

        Force-Modus (force=True):
          - Ignoriere die Score-Vergleiche bei Überlappung.
          - Entferne alle überlappenden Domänen und füge die neue ein (return 1).

        Returns:
            1 wenn hinzugefügt, 0 wenn nicht hinzugefügt.
        """
        del_domains = []  # start-Koordinaten der zu entfernenden Domänen

        for domain in self.domains.values():
            if self.check_domain_overlap(
                start, end, domain.get_start(), domain.get_end()
            ):
                if force:
                    # im Force-Modus: jede überlappende Domäne räumen
                    del_domains.append(domain.get_start())
                else:
                    # Standard: nur schwächere Domänen räumen, sonst abbrechen
                    if domain.get_score() < score:
                        del_domains.append(domain.get_start())
                    else:
                        return 0

        # überlappende domänen verschieben in deleted_domains
        for key in del_domains:
            dom = self.domains.pop(key, None)
            if dom is not None:
                key = dom.domain
                self.deleted_domains[key] = dom

        # neue Domäne eintragen (Schlüssel = start)
        self.domains[start] = Domain(hmm, start, end, score, ident, bsr)

        return 1


class Domain:
    # 2.9.22
    """
    Stores domain information (HMM name, coordinates, score, identity, bsr).

    Args:
        domain (str): Domain name.
        start (int): Start coord.
        end (int): End coord.
        score (float): Bitscore.
        ident (int): Percent identity.
        bsr (float): Blast score ratio.
    """

    def __init__(
        self,
        domain: str,
        start: int,
        end: int,
        score: float,
        ident: int = 1,
        bsr: float = 1.0,
    ):
        self.domain: str = domain
        self.start: int = int(start)
        self.end: int = int(end)
        self.score: float = float(score)
        self.identity: int = int(ident)
        self.bsr: float = float(bsr)

    def __hash__(self):
        return hash((self.domain, self.start, self.end, self.score))

    def __eq__(self, other):
        if isinstance(other, Domain):
            return (
                self.domain == other.domain
                and self.start == other.start
                and self.end == other.end
                and self.score == other.score
            )
        return False

    def get_domain(self):
        return self.domain

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_score(self):
        return self.score


#########################################
########   Parsing subroutines ##########
#########################################


def parse_gff_file(filepath: str, protein_dict: Dict[str, Protein]) -> Dict[str, Protein]:
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
    gene_id_pattern = re.compile(r"ID=(cds-)?(\S+?)(?:[;\s]|$)")

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
            match = gene_id_pattern.search(gff[-1])
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
                locustag = get_locustag(locustag_pattern, line)
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

        for line in reader: # type: str
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


def get_locustag(locustag_pattern: re.Pattern, string: str) -> str:
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

def output_genome_report(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, str],
    genomeID: str = "",
    writemode: str = "w",
    taxon_divider: str = ".",
) -> None:
    """
    Writes the main genome hit table (TSV).

    Columns:
      genomeID, proteinID, domains, dom_scores, dom_coords, contig, gene_start, gene_end,
      gene_strand, locustag, clusterID, csb,
      Superkingdom, Phylum, Class, Order, Family, Genus, Species
    """
    taxon_cols = [
        "genomeID",
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    header_cols = [
        "genomeID",
        "proteinID",
        "domains",
        "dom_scores",
        "dom_coords",
        "contig",
        "gene_start",
        "gene_end",
        "gene_strand",
        "locustag",
        "selection_comment",
        "alternative hit",
        "clusterID",
        *taxon_cols,
    ]
    header = "\t".join(header_cols)

    # sortiert nach genomeID, contig, start
    proteinID_list = sorted(
        protein_dict,
        key=lambda x: (
            protein_dict[x].genomeID,
            protein_dict[x].gene_contig,
            protein_dict[x].gene_start,
        ),
    )

    with open(output_filepath, writemode) as writer:
        writer.write(header + "\n")

        for proteinID in proteinID_list:
            protein = protein_dict[proteinID]
            # protein.get_protein_list() erwartete Reihenfolge laut Docstring:
            # [proteinID, get_domains(), get_domain_scores(), get_domain_coordinates(),
            #  gene_contig, gene_start, gene_end, gene_strand, gene_locustag]
            #if not protein.valid_hit:
            #    continue

            pl = protein.get_protein_list()

            # clusterID + csb (ohne keyword/completeness)
            clusterID = protein.clusterID
            csb_val = ""
            if clusterID:
                if clusterID in cluster_dict:
                    # cluster.get_cluster_list(",") liefert i. d. R. [clusterID, keyword, completeness, csb]
                    cl = cluster_dict[clusterID].get_cluster_list(",")
                    # defensiv extrahieren:
                    out_clusterID = cl[0] if len(cl) >= 1 else clusterID
                    csb_val = (
                        cl[3] if len(cl) >= 4 else (cl[-1] if len(cl) >= 2 else "")
                    )
                else:
                    out_clusterID = clusterID
            else:
                out_clusterID = ""

            # Taxonomie in 7 Spalten
            taxon_dict = {} if taxon_dict is None else taxon_dict

            gid = protein.genomeID
            taxon_levels = [""] * 8
            if gid in taxon_dict and taxon_dict[gid]:
                parts = [p.strip() for p in str(taxon_dict[gid]).split(taxon_divider)]
                for i in range(min(8, len(parts))):
                    taxon_levels[i] = parts[i]
            row = [
                protein.genomeID,  # genomeID
                pl[0],  # proteinID
                pl[1],  # domains
                pl[2],  # dom_scores
                pl[3],  # dom_coords
                pl[4],  # contig
                pl[5],  # gene_start
                pl[6],  # gene_end
                pl[7],  # gene_strand
                pl[8],  # locustag
                protein.get_selection_comment_csv(),
                protein.alternative_hit,
                out_clusterID,
                csb_val,
                *taxon_levels,
            ]
            writer.write("\t".join(map(str, row)) + "\n")
    return

###############################################################################################################
###############################################################################################################

# Writer for hit reports


def submit_batches(protein_batch, cluster_batch, options):
    # Insert into the database
    database.insert_database_proteins(options.database_directory, protein_batch)
    database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + "\t" + "\t".join(domains) + "\n")


def process_writer(queue, options):
    # This routine handles the output of the search and writes it into the database
    # It gets input from multiple workers as the database connection to sqlite is unique

    global cluster_dict, protein_dict
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
            logger.debug(f"Processed {batch_counter} genomes ")  #

        protein_dict, cluster_dict = tup

        # Concatenate the data
        protein_batch.update(protein_dict)
        cluster_batch.update(cluster_dict)

        # Print text reports if desired
        if options.individual_reports:
            if protein_dict:  # Check if protein_dict is not empty
                first_protein_key = next(iter(protein_dict))  # Get the first key
                genomeID = protein_dict[first_protein_key].genomeID
                filepath = os.path.join(
                    options.fasta_initial_hit_directory,
                    str(genomeID) + ".hit_table_txt",
                )
                output_genome_report(filepath, protein_dict, cluster_dict, {})

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

                output_genome_report(filepath, protein_dict, cluster_dict, {})
    logger.info(f"Processed {batch_counter} genomes")
    return


#########################################################################################
################ Processing routines for parsing genome hits ############################
#########################################################################################


def main_parse_summary_hmmreport(config):
    genome_ids = list(config.queued_genomes)

    genomeID_batches = split_into_batches(genome_ids, config.cores - 1)

    # Lade Patterns nur 1x im Hauptprozess
    csb_patterns = csb_finder.make_pattern_dict(config.patterns_file)
    cooccurrence_pattern = csb_finder.make_pattern_dict(config.cooccurrence_file)
    threshold_dict = search_cross_reference.make_threshold_dict(
        config.score_threshold_file, 3, config.thrs_score
    )
    exclusion_singletons = parse_exclusion_singletons(config.exclusion_singletons)

    # Make csb naming index table
    only_pattern_dict = {name: patset for name, (patset, _) in csb_patterns.items()}
    index_trie = csb_trie_algorithm.build_trie_index(only_pattern_dict)

    # Insert genomeIDs in DB
    database.insert_database_genome_ids(config.database_directory, set(genome_ids))

    with Manager() as manager:
        data_queue = manager.Queue()

        with Pool(processes=config.cores) as pool:
            # Start writer
            p_writer = pool.apply_async(process_writer, (data_queue, config))
            args = [
                (
                    data_queue,
                    batch,
                    config.faa_files,
                    config.gff_files,
                    config.glob_trusted_hitreport,
                    config.glob_intermediate_hitreport,
                    config.nucleotide_range,
                    config.min_completeness,
                    csb_patterns,
                    cooccurrence_pattern,
                    exclusion_singletons,
                    threshold_dict,
                    index_trie,
                    config.use_synteny_completion,
                    config.use_remove_unassigned_intermediates,
                    config.use_remove_unassigned_singletons,
                )
                for batch in genomeID_batches
            ]

            # Start readers
            pool.starmap(process_batch, args)

            # End writer
            for _ in range(config.cores):
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
    cooccurrence_pattern,
    exclusion_singletons,
    threshold_dict,
    index_trie,
    use_synteny_completion: bool,
    use_remove_unassigned_intermediates: bool,
    use_remove_unassigned_singletons: bool,
):
    """
    Runs process_genome for each genome in the batch with all necessary arguments.
    """
    for genome_id in genome_ids:
        try:
            process_genome(
                data_queue=data_queue,
                genome_id=genome_id,
                faa_file=faa_files[genome_id],
                gff_file=gff_files[genome_id],
                trusted_hmmreport_path=trusted_hmmreport_path,
                intermediate_hmmreport=intermediate_hmmreport,
                nucleotide_range=nucleotide_range,
                min_completeness=min_completeness,
                pattern_dict=pattern_dict,
                cooccurrence_pattern=cooccurrence_pattern,
                exclusion_singletons=exclusion_singletons,
                threshold_dict=threshold_dict,
                index_trie=index_trie,
                use_synteny_completion=use_synteny_completion,
                use_remove_unassigned_intermediates=use_remove_unassigned_intermediates,
                use_remove_unassigned_singletons=use_remove_unassigned_singletons,
            )
        except Exception as e:
            logger.warning(f"Failed to process genome '{genome_id}' — {str(e)}")
            continue


def process_genome(
    data_queue: Any,
    genome_id: str,
    faa_file: str,
    gff_file: str,
    trusted_hmmreport_path: str,
    intermediate_hmmreport: str,
    nucleotide_range: int,
    min_completeness: float,
    pattern_dict: dict[str, tuple[set[str], int]],
    cooccurrence_pattern: dict[str, tuple[set[str], int]],
    exclusion_singletons: set[str],
    threshold_dict: dict[str, float],
    index_trie: TrieIndex,
    use_synteny_completion: bool,
    use_remove_unassigned_intermediates: bool,
    use_remove_unassigned_singletons: bool,
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
        trusted_hmmreport_path:
        intermediate_hmmreport:
        use_synteny_completion:
        use_remove_unassigned_intermediates (bool):
        use_remove_unassigned_singletons:
        index_trie: trie for the csb naming routine, for faster lookup
        cooccurrence_pattern (Dict[str, tuple[set[str],int]]):
        data_queue: Multiprocessing queue for result transport.
        genome_id: ID of the genome.
        faa_path, gff_path: Paths to input files (can be .gz).
        trusted_hmmreport_path, intermediate_hmmreport: HMM report files.
        nucleotide_range: Nucleotide window for cluster detection.
        min_completeness: Minimum completeness for cluster pattern assignment.
        pattern_dict: Patterns for cluster annotation.
        exclusion_singletons: Domains for singleton exclusion.
        threshold_dict: Score cutoffs per domain.

        dieser teil hier ist notorisch langsam für besonders große glob files
        Der intermediate file sollte das gleiche sein wie der hmmreport, nur, dass hier
        auch die Tc mit drinstehen. Tc ist aber deutlich kleiner 300 MB vs 9.7 Gb
        Die 300 MB können aber trotzdem mal sortiert und indexiert werden für schnelleren lookup.
        Außerdem einmal den intermediate file auseinander nehmen und in die hmmreports schreiben .nc_hmmreport

        Diese routinen sind etwa 10 mal langsamer als das lookup:
        syntenic block completion
        find and name sytenic block sind die langsamsten schritte

        scheint auch ein fehler bei der erkennung von SQR SDO und co zu haben. die werden im intermediate nicht aufgeführt, was komisch ist.
    """
    try:
        with tick(f"read concatenated trusted hmmreport {genome_id}"):
            # Intermediate protein hits
            intermediate_protein_dict = parse_bulk_HMMreport_genomize(
                genome_id, intermediate_hmmreport
            )
            parse_gff_file(gff_file, intermediate_protein_dict)

        with tick(f"read concatenated intermediate hmmreport {genome_id}"):
            # Primary protein hits
            trusted_protein_dict = parse_bulk_HMMreport_genomize(
                genome_id, trusted_hmmreport_path
            )
            parse_gff_file(gff_file, trusted_protein_dict)

        # Combine protein dictionaries
        combined_protein_dict = {**intermediate_protein_dict, **trusted_protein_dict}

        with tick(f"Find syntenic blocks {genome_id}"):
            # Detect and annotate syntenic gene clusters
            cluster_dict = csb_finder.find_syntenic_blocks(
                genome_id, combined_protein_dict, nucleotide_range
            )
        with tick(f"Name syntenic blocks {genome_id}"):

            #cluster_dict = csb_finder.name_syntenic_blocks(pattern_dict, cluster_dict, min_completeness)
            cluster_dict = csb_finder.name_syntenic_blocks_trie(cluster_dict, index_trie, min_completeness=min_completeness)

        with tick(f"Add selection comment trusted proteins {genome_id}"):
            # Attach the reason for selection to protein objects trusted cutoff/reference sequence
            trusted_cutoff_protein_ids = set(trusted_protein_dict.keys())
            combined_protein_dict = add_selection_comment_to_many_proteins(
                combined_protein_dict, trusted_cutoff_protein_ids, "Tc"
            )

        # Enhance cluster completeness if needed with synteny correction
        # Optional: alters the combined_protein_dict
        if use_synteny_completion:
            with tick(f"Increase syntenic block completeness {genome_id}"):
                pattern_completion_synteny.enhance_syntenic_block_completeness(cluster_dict, combined_protein_dict,
                                                                           pattern_dict)

        with tick(f"Enhance pathway completeness {genome_id}"):
            # Collect trusted protein IDs: those with complete pathways and those in the main protein dict
            singletons_with_complete_pathway_set = (
                pattern_completion_pathway.enhance_pathway_completeness(
                    combined_protein_dict, cooccurrence_pattern, threshold_dict
                )
            )

        with tick(f"Comment selection criteria {genome_id}"):
            # Attach the reason for selection to protein objects complete pathway
            combined_protein_dict = add_selection_comment_to_many_proteins(
                combined_protein_dict, singletons_with_complete_pathway_set, "Coo"
            )

            trusted_protein_ids = trusted_cutoff_protein_ids.union(
                singletons_with_complete_pathway_set
            )

        if use_remove_unassigned_intermediates:
            with tick(f"Remove unassigned intermediate hits {genome_id}"):
                # Remove unassigned intermediate proteins, but keep trusted ones and hits in named gene clusters
                combined_protein_dict = remove_unassigned_intermediate_proteins(
                    combined_protein_dict, trusted_protein_ids, cluster_dict
                )
        if use_remove_unassigned_singletons:
            with tick(f"Remove unassigned singletons {genome_id}"):
                # Remove genes that should not occur as singletons
                # alters the combined_protein_dict but ignores singletons that complete pathway
                remove_exclusion_singletons(
                    combined_protein_dict,
                    cluster_dict,
                    exclusion_singletons,
                    singletons_with_complete_pathway_set,
                )

        with tick(f"Attach protein sequences {genome_id}"):
            # Attach protein sequences
            get_protein_sequence(faa_file, combined_protein_dict)

        data_queue.put((combined_protein_dict, cluster_dict))

    except Exception as e:
        logger.error(f"Error: {genome_id} -> {e}")
        logger.error(traceback.format_exc())


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
                if domain.get_domain() in exclusion_singletons:
                    to_remove.add(protein_id)
                    break

    for protein_id in to_remove:
        protein_obj = combined_protein_dict[protein_id]
        domains = protein_obj.get_domains()  # Gibt z.B. 'HMM_A-HMM_B-HMM_C' zurück
        logger.debug(f"Removed singleton {protein_id} with domains: {domains}")
        del combined_protein_dict[protein_id]


def remove_unassigned_intermediate_proteins(
    combined_protein_dict: Dict[str, Any], trusted_protein_ids: set, cluster_dict: Dict[str, Any]
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

    # Update the trusted proteinIDs with the proteinIDs
    # that are covered by recognized patterns
    for cluster in cluster_dict.values():
        trusted_protein_ids.update(getattr(cluster, "covered_protein_ids", set()))

    # New 091125 mark up the intermediate hits instead of remove
    for pid, protein in combined_protein_dict.items():
        if pid in trusted_protein_ids:
            protein.valid_hit = True
            protein.add_selection_comment("Sc")
        else: # was not in trusted hits nor in a recognized gene cluster
            protein.valid_hit = False
            protein.add_selection_comment("Nc")

    return combined_protein_dict

    # Remove any protein in combined_protein_dict that is not in proteinIDs
    for protein_id in list(combined_protein_dict.keys()):
        if protein_id not in trusted_protein_ids:
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


def parse_bulk_HMMreport_genomize(genomeID, database, protein_dict=None):
    """
    Parse for a single genomeID the results from raw hmmreport results database
    Args:
        genomeID:
        database:
        protein_dict:

    Returns:

    """
    if protein_dict is None:
        protein_dict = {}

    results = report_db.get_hits_by_genome_from_report_db(database, genomeID)
    ProteinClass = Protein  # lookup once

    for (
        combined_id,
        genome_id,
        protein_id,
        query,
        hit_bitscore,
        hsp_start,
        hsp_end,
    ) in results:
        protein = protein_dict.get(protein_id)
        if protein is None:
            protein = ProteinClass(
                protein_id, query, hsp_start, hsp_end, hit_bitscore, genome_id
            )
            protein_dict[protein_id] = protein
        else:
            protein = protein_dict[protein_id]
            protein.add_domain(query, hsp_start, hsp_end, hit_bitscore)

    return protein_dict

#!/usr/bin/python

import re
import subprocess

from typing import Dict, Any, Set

from hmsss.core.logging import get_logger
logger = get_logger(__name__)

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
        self.valid_hit = False

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
        listing = [
            self.proteinID,
            self.get_domains(),
            str(self.get_domain_scores()),
            str(self.get_domain_coordinates()),
            self.gene_contig,
            str(self.gene_start),
            str(self.gene_end),
            self.gene_strand,
            self.gene_locustag,
        ]
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


def parse_gff_file(
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

    """
    reader = None
    try:
        reader = open(filepath, "r")
        sequence = ""
        header = None
        save_sequence = False

        for line in reader:  # type: str
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
    taxon_divider: str = "\t",
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
        #"valid hit",
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
            tax_levels = [""] * len(taxon_cols)

            if taxon_dict and gid in taxon_dict and isinstance(taxon_dict[gid], dict):
                for i, level in enumerate(taxon_cols):
                    tax_levels[i] = taxon_dict[gid].get(level, "")
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
                #protein.valid_hit,
                protein.get_selection_comment_csv(),
                protein.alternative_hit,
                out_clusterID,
                csb_val,
                *tax_levels,
            ]
            writer.write("\t".join(map(str, row)) + "\n")
    return




#########################################################################################
################ Processing routines for parsing genome hits ############################
#########################################################################################


def remove_unassigned_intermediate_proteins(
    combined_protein_dict: Dict[str, Any],
    trusted_protein_ids: set,
    cluster_dict: Dict[str, Any],
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
        if protein.valid_hit is True:
            continue # Skip proteins that are already recognized
        if pid in trusted_protein_ids:
            protein.valid_hit = True
            protein.add_selection_comment("Sc")
        else:  # was not in trusted hits nor in a recognized gene cluster
            protein.valid_hit = False
            protein.add_selection_comment("Nc")

    return combined_protein_dict
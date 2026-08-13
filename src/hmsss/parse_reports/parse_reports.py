#!/usr/bin/python
from __future__ import annotations

import re
import subprocess
from bisect import bisect_right
from dataclasses import dataclass, field, replace

from typing import Dict, Any, Set, List, Tuple, FrozenSet

from hmsss.core.logging import get_logger

logger = get_logger(__name__)


@dataclass(slots=True)
class Domain:
    domain: str
    start: int
    end: int
    score: float
    selection_comment_list: List[str] = field(default_factory=list)
    identity: int = 1
    bsr: float = 1.0

    # 2.9.22
    def __hash__(self):
        return hash((self.domain, self.start, self.end, self.score))

    def add_selection_comment(self, selection_comment: str):
        self.selection_comment_list.append(selection_comment)
        return

    def get_domain(self):
        return self.domain

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_score(self):
        return self.score


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
        selection_comment: str = "",
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
        self.domains: Set[Domain] = (
            set()
        )  # Domain objects that are included due to high score
        self.low_score_domains: Set[Domain] = (
            set()
        )  # Domain objects that are excluded by score
        self.add_domain(
            hmm, start, end, score, ident, bsr, selection_comment=selection_comment
        )
        self.selection_comment = ""
        self.alternative_hit: str = ""
        self.valid_hit = False

    ##### Getter ####

    def get_domains(self) -> str:
        listing = [dom.domain for dom in sorted(self.domains, key=lambda d: d.start)]
        # collapse identical domains, e.g. SoxX-SoxX -> SoxX
        if listing and len(set(listing)) == 1:
            return listing[0]
        return "-".join(listing)

    def get_domain_listing(self):
        # return list
        return list(self.domains)

    def get_domain_coordinates(self) -> str:
        listing = [
            f"{dom.start}:{dom.end}"
            for dom in sorted(self.domains, key=lambda d: (d.start, d.end))
        ]
        return ";".join(listing)

    def get_domain_scores(self) -> str:
        listing = [
            f"{dom.score}" for dom in sorted(self.domains, key=lambda d: d.start)
        ]
        return ";".join(listing)

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
    def _check_domain_overlap(a: Domain, b: Domain) -> bool:
        """
        Explicit overlap checks for two domains with inclusive coordinates.
        a = existing/current domain
        b = new domain
        """

        a_start, a_end = a.start, a.end
        b_start, b_end = b.start, b.end

        # --- Fall 1: Start oder Ende von b liegt innerhalb von a ---
        if a_start <= b_start <= a_end:
            return True
        if a_start <= b_end <= a_end:
            return True

        # --- Fall 2: eine Domäne liegt vollständig in der anderen ---
        if b_start <= a_start and a_end <= b_end:
            return True
        if a_start <= b_start and b_end <= a_end:
            return True

        # --- Fall 3: klar getrennt (b komplett links oder rechts von a) ---
        if b_end < a_start:
            return False
        if b_start > a_end:
            return False

        # --- Sollte logisch nie erreicht werden ---
        return False

    @staticmethod
    def best_nonoverlapping_domain_set(
        self,
        domains: Set[Domain],
        *,
        inclusive: bool = True,
    ) -> Set[Domain]:
        # intern als Liste arbeiten
        doms = sorted(domains, key=lambda d: (d.end, d.start))
        ends = [d.end for d in doms]

        def compatible_end_value(start: int) -> int:
            return start - 1 if inclusive else start

        p: List[int] = []
        for i, d in enumerate(doms):
            j = (
                bisect_right(
                    ends,
                    compatible_end_value(d.start),
                    0,
                    i,
                )
                - 1
            )
            p.append(j)

        n = len(doms)
        dp: List[Tuple[float, int]] = [(0.0, 0)] * n
        take: List[bool] = [False] * n

        def dom_len(d: Domain) -> int:
            return (d.end - d.start + 1) if inclusive else (d.end - d.start)

        def better(a: Tuple[float, int], b: Tuple[float, int]) -> bool:
            if a[0] != b[0]:
                return a[0] > b[0]
            return a[1] > b[1]

        def effective_score(d: Domain) -> float:
            score = float(d.score)

            if (
                "Bc" in d.selection_comment_list
            ):  # if below minimal cutoff, add penalty to domain score for calculation
                return score * 0.01

            return score

        for i, d in enumerate(doms):
            best_skip = dp[i - 1] if i > 0 else (0.0, 0)

            prev = dp[p[i]] if p[i] >= 0 else (0.0, 0)
            best_take = (
                prev[0] + effective_score(d),
                prev[1] + dom_len(d),
            )

            if better(best_take, best_skip):
                dp[i] = best_take
                take[i] = True
            else:
                dp[i] = best_skip
                take[i] = False

        # Reconstruct
        chosen: set[Domain] = set()

        i = n - 1
        while i >= 0:
            if take[i]:
                chosen.add(doms[i])
                i = p[i]
            else:
                i -= 1

        # ------------------------------------------------------------------
        # Collapse consecutive selected domains of the same type
        # ------------------------------------------------------------------

        chosen_sorted = sorted(chosen, key=lambda d: (d.start, d.end))

        if not chosen_sorted:
            return set()

        collapsed: set[Domain] = set()

        current_group = [chosen_sorted[0]]

        def collapse_group(group: list[Domain]) -> Domain:
            """
            Collapse consecutive domains of the same type into one Domain.

            - domain:    same domain name
            - start:     start of first domain
            - end:       end of last domain
            - score:     score of first domain
            - comments:  unique comments from all merged domains
            - identity:  identity of first domain
            - bsr:       BSR of first domain
            """

            first = group[0]
            last = group[-1]

            # Merge selection comments, preserving their order
            comments = []

            for dom in group:
                for comment in dom.selection_comment_list:
                    if comment not in comments:
                        comments.append(comment)

            return Domain(
                domain=first.domain,
                start=first.start,
                end=last.end,
                score=first.score,
                selection_comment_list=comments,
                identity=first.identity,
                bsr=first.bsr,
            )

        for dom in chosen_sorted[1:]:
            if dom.domain == current_group[-1].domain:
                # same domain type -> add to current group
                current_group.append(dom)

            else:
                # different domain type -> finish previous group
                collapsed.add(collapse_group(current_group))

                current_group = [dom]

        # Finish final group
        collapsed.add(collapse_group(current_group))

        return collapsed

    def add_selection_comment_to_domain(
        self,
        domain_name: str,
        comment: str,
        sep: str = ",",
    ) -> None:
        """
        Add one or multiple selection comment tokens to a specific domain.

        Parameters
        ----------
        domain_name : str
            Name of the domain (HMM) to which the comment should be added.
        comment : str
            Single token ("Tc") or CSV string ("Tc,Coo").
        sep : str
            Separator used in comment string.
        """

        comments = [
            t.strip() for t in str(comment).split(sep) if t.strip()
        ]  # split the comments

        for dom in self.domains | self.low_score_domains:
            if dom.domain != domain_name:
                continue

            # additive, no duplicates
            for comment in comments:
                if comment not in dom.selection_comment_list:
                    dom.selection_comment_list.append(comment)

    def add_domain(
        self,
        hmm: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0,
        *,
        selection_comment: str = "",
        force: bool = False,
    ) -> None:
        """
        Adds a domain to the protein.

        """
        # print(
        #    f"{hmm}\t{start}\t{end}\t{score}\tidentity {ident}\t{bsr}\tselection comment {selection_comment}"
        # )
        added_domain = Domain(
            hmm,
            start,
            end,
            score,
            identity=ident,
            bsr=bsr,
        )
        added_domain.add_selection_comment(selection_comment)
        if force:
            to_remove = set()
            for current_domain in self.domains:
                if self._check_domain_overlap(added_domain, current_domain):
                    to_remove.add(current_domain)

            # Namen der entfernten Domains extrahieren
            removed_names = [dom.domain for dom in to_remove]
            if removed_names:
                if self.alternative_hit:
                    self.alternative_hit += "-" + "-".join(removed_names)
                else:
                    self.alternative_hit = "-".join(removed_names)

            # Domains wirklich entfernen
            self.domains.difference_update(to_remove)

            # neue Domäne hinzufügen
            self.domains.add(added_domain)
        else:
            self.low_score_domains.add(added_domain)

    def define_best_scoring_domains(self) -> None:
        self.domains = self.best_nonoverlapping_domain_set(
            self, self.low_score_domains, inclusive=False
        )

    def define_selection_comment(self) -> None:
        """
        Collect selection comments from all domains (sorted by start coordinate)
        and store them as a single string on the Protein.
        """
        parts: list[str] = []

        for dom in sorted(self.domains, key=lambda d: (d.start, d.end)):
            # dom.selection_comment is e.g. frozenset[str] (or set[str])
            if not dom.selection_comment_list:
                continue

            # deterministic order within a domain
            joined_domain_comment = "-".join(dom.selection_comment_list)
            parts.append(joined_domain_comment)

        if parts:
            # join domains in genomic order
            self.selection_comment = ";".join(parts)


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
        # "valid hit",
        "selection_comment",
        "alternative hit",
        "clusterID",
        *taxon_cols,
    ]
    header = "\t".join(header_cols)

    # sortiert nach genomeID, contig, start
    proteinID_list = sorted(
        [
            proteinID
            for proteinID, protein in protein_dict.items()
            if not genomeID or protein.genomeID == genomeID
        ],
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
            # if not protein.valid_hit:
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
                # protein.valid_hit,
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


def update_protein_validity_by_synteny(
    combined_protein_dict: Dict[str, Any],
    trusted_protein_ids: set,
    cluster_dict: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Mark proteins from combined_protein_dict that are trusted hits
    or in the set of covered_protein_ids from any cluster as valid hits.

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
            continue  # Skip proteins that are already recognized
        if pid in trusted_protein_ids:
            protein.valid_hit = True
            for domain in protein.domains:
                domain.add_selection_comment("Sb")  # Detected syntenic block
        else:  # was not in trusted hits nor in a recognized gene cluster
            protein.valid_hit = False
            for domain in protein.domains:
                domain.add_selection_comment("Nb")  # No syntenic block

    return combined_protein_dict


def remove_invalid_bc_only_proteins(
    protein_dict: dict[str, Protein],
) -> dict[str, Protein]:
    to_remove = []

    allowed_comments = {"Bc", "Nb"}

    for protein_id, protein in protein_dict.items():
        if protein.valid_hit:
            continue

        if not protein.domains:
            continue

        only_bc_like = all(
            set(domain.selection_comment_list).issubset(allowed_comments)
            and "Bc" in set(domain.selection_comment_list)
            for domain in protein.domains
        )

        if only_bc_like:
            to_remove.append(protein_id)

    for protein_id in to_remove:
        del protein_dict[protein_id]

    return protein_dict


def define_best_score_hits_for_protein_dict(
    protein_dict: dict[str, Protein],
):
    for protein in protein_dict.values():
        protein.define_best_scoring_domains()


def define_selection_comments_for_protein_dict(
    protein_dict: dict[str, Protein],
):
    for protein in protein_dict.values():
        protein.define_selection_comment()

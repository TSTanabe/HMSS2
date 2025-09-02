#!/usr/bin/python
from itertools import chain
from typing import Any, Dict, List, Optional, Set, Tuple

from hmsss.algorithms.csb_trie_algorithm import TrieIndex
from hmsss.utils import myUtil
from hmsss.algorithms import csb_trie_algorithm

logger = myUtil.log


class Cluster:
    """
    3.9.22
    Cluster organises genes laying in synteny with a specific order. Each cluster has a unique clusterID derived from the assemblyID but with an added index number.

    Args
        clusterID = unique string

    11.04.23 added the cluster_start and cluster_end lines
    """

    def __init__(self, cluster_id: str, distance: int = 3500) -> None:
        self.clusterID: str = cluster_id
        self.genomeID: str = ""
        self.contig: str = ""
        self.distance: int = distance
        self.genes: List[str] = []
        self.types: List[str] = []
        self.keywords: List[Keyword] = []
        self.keywords_dict: Dict[Any, Keyword] = {}
        self.cluster_start: Optional[int] = None
        self.cluster_end: Optional[int] = None

        # Necessary for the pattern matching
        self.type_to_proteins: Dict[
            str, Set[str]
        ] = {}  # Saves type => set(proteinID). Necessary to get proteinIDs with pattern matching types
        self.covered_protein_ids: Set[str] = (
            set()
        )  # proteinIDs which are part of matching patterns

    def add_gene(
        self,
        protein_id: str,
        types: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> None:
        self.genes.append(protein_id)
        self.types.append(types)

        if self.cluster_start is None or (
            start is not None and start < self.cluster_start
        ):
            self.cluster_start = start
        if self.cluster_end is None or (end is not None and end > self.cluster_end):
            self.cluster_end = end

        if types not in self.type_to_proteins:
            # Define the proteinIDs per protein type
            self.type_to_proteins[types] = set()
        self.type_to_proteins[types].add(protein_id)

    def add_keyword(
        self,
        keyword: str,
        completeness: float = 0.0,
        csb: str = ".",
        missing: Optional[tuple[str]] = None,
        additional_elements: Optional[tuple[str]] = None,
        keyword_id: Any = ".",
    ) -> None:
        missing = missing or []
        keyword_id = keyword if keyword_id == "." else keyword_id
        self.keywords_dict[keyword_id] = Keyword(
            keyword, completeness, csb, missing, additional_elements, keyword_id
        )

    def add_covered_types(self, covered_types: Set[str] | List[str] | tuple[str, ...]) -> Set[str]:
        """
        Given a set/list/tuple of covered protein 'types' (domain names),
        look up their protein IDs via `type_to_proteins` and accumulate them
        into `self.covered_protein_ids`.

        Returns the set of IDs that were added on this call.
        """
        # Accept any iterable of strings; avoid accidental string-iteration
        if isinstance(covered_types, str):
            covered_types = {covered_types}
        else:
            covered_types = set(covered_types)

        newly_added: Set[str] = set()
        for typ in covered_types:
            ids = self.type_to_proteins.get(typ)
            if ids:
                newly_added.update(ids)

        # Accumulate into the cluster-level set
        self.covered_protein_ids.update(newly_added)
        return newly_added

    def get_keywords(self) -> List["Keyword"]:
        return list(self.keywords_dict.values())

    def get_genes(self) -> List[str]:
        return self.genes

    def get_domains(self) -> List[str]:
        return self.types

    def get_domain_set(self) -> Set[str]:
        tmp = "-".join(self.types)
        return set(tmp.split("-"))

    def get_cluster_id(self):
        return self.clusterID

    def get_protein_type_list(self) -> List[str]:
        """
        Returns a list of protein types behind the last '_' in each entry of self.types.
        Example: "grp0_ProteinA" → "ProteinA", "TIGRFAM0000_ProteinB" → "ProteinB".
        """
        return [typus.split("_")[-1] if "_" in typus else typus for typus in self.types]

    def get_protein_type_set(self) -> Set[str]:
        """
        Returns a set of unique protein types behind the last '_' in each entry of self.types.
        Example: "grp0_ProteinA" → "ProteinA", "TIGRFAM0000_ProteinB" → "ProteinB".
        """
        return {typus.split("_")[-1] if "_" in typus else typus for typus in self.types}

    def get_cluster_list(self, separator: str) -> List[str]:
        keywords_string = ""
        completeness_string = ""
        csb_string = ""
        listing = list(self.keywords_dict.values())
        if listing:
            element = listing.pop(0)
            keywords_string += str(element.get_keyword())
            completeness_string += str(element.get_completeness())
            csb_string += str(element.get_csb())
        for element in listing:
            keywords_string += separator + str(element.get_keyword())
            completeness_string += separator + str(element.get_completeness())
            csb_string += separator + str(element.get_csb())
        return [self.get_cluster_id(), keywords_string]


class Keyword:
    """
    3.9.22
    Holds the information about a single keyword, its completeness and if it is a csb.
    """

    def __init__(
        self,
        keyword: str,
        completeness: float = 0.0,
        csb: str = ".",
        missing: Optional[tuple[str]] = None,
        additional: Optional[tuple[str]] = None,
        keyword_id: Any = ".",
        transitions=None,
    ) -> None:
        self.keyword: str = str(keyword)
        self.csb: str = csb  # Collinear to the reference pattern
        self.completeness: float = completeness
        self.missing_domains: tuple[str] = missing if missing is not None else ()
        self.additional_domains: tuple[str] = (
            additional if additional is not None else ()
        )
        self.transitions: dict[str, tuple[str, ...]] = (
            transitions if transitions is not None else {}
        )  # additional_domain => (missing_domain1, missing_domain2, ...)
        self.keyword_id: Any = keyword_id

    def get_keyword(self) -> str:
        return self.keyword

    def get_csb(self) -> str:
        return self.csb

    def get_completeness(self) -> float:
        return self.completeness

    def get_missing_domains(self) -> tuple:
        return self.missing_domains

    def get_additional_domains(self) -> tuple:
        return self.additional_domains

    def get_transitions(self, key: Optional[str] = None):
        """
        Gibt entweder das gesamte transitions-Dict zurück oder,
        wenn `key` angegeben ist, nur das Tupel für diesen Key.
        """

        if key is None:
            return self.transitions
        return self.transitions.get(key, ())

    def set_keyword(self, keyword: str) -> None:
        self.keyword = keyword

    def set_csb(self, csb: str) -> None:
        self.csb = csb

    def set_completeness(self, completeness: float) -> None:
        self.completeness = completeness

    def set_transition(self, transition: dict[str, tuple[str, ...]]) -> None:
        self.transitions = transition


def make_pattern_dict(filepath: str) -> Dict[str, Tuple[Set[str], int]]:
    """
    Reads a file containing patterns and their names.

    Returns:
        pattern: Dict[name, (set_of_patterns, length_of_pattern)]
        pattern_names: Dict[name, original_name]
    """
    pattern: Dict[str, Tuple[Set[str], int]] = {}

    try:
        with open(filepath, "r") as reader:
            for line_num, line in enumerate(reader, start=1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if not parts:
                    continue

                name = parts[0].strip()
                patset = {item.strip() for item in parts[1:] if item.strip()}

                if not patset:
                    continue  # leere Patternzeilen überspringen

                if name in pattern:
                    logger.warning(
                        f"Duplicate pattern name '{name}' found on line {line_num}, skipping."
                    )
                    continue

                pattern[name] = (patset, len(patset))

    except FileNotFoundError:
        logger.error(f"File not found: {filepath}")
        return {}
    except IOError as e:
        logger.error(f"Unable to read file: {filepath}. {e}")
        return {}

    return pattern


def find_syntenic_blocks(
    genome_id: str, protein_dict: Dict[str, Any], distance: int = 3500
) -> Dict[str, Cluster]:
    """
    Drop-in replacement: gleicher Output wie vorher, aber schneller.
    Logik:
      - sortiere (contig, start)
      - starte einen Cluster erst dann, wenn eine echte Synergie (<= distance, gleicher contig) vorliegt
      - füge solange hinzu, bis Bedingung bricht; dann Cluster abschließen
    """
    # Lokale Aliasse minimieren Kosten:
    d = protein_dict
    dist = distance

    # Sortierung nur nach wirklich benötigten Keys:
    # (kein str() auf IDs; greift direkt die Objekte)
    protein_id_list = sorted(d, key=lambda pid: (d[pid].gene_contig, d[pid].gene_start))

    clusters: Dict[str, Cluster] = {}
    cluster_idx = 0
    cur_cluster: Optional[Cluster] = None

    # "Vorheriges" Protein für Adjazenzprüfung
    prev_pid: Optional[str] = None
    prev_prot = None  # type: ignore

    for pid in protein_id_list:
        p = d[pid]
        if prev_pid is None:
            prev_pid, prev_prot = pid, p
            continue

        same_contig = prev_prot.gene_contig == p.gene_contig
        close_enough = same_contig and (p.gene_start - prev_prot.gene_end <= dist)

        if close_enough:
            # Start eines neuen Clusters?
            if cur_cluster is None:
                cluster_idx += 1
                cur_cluster = Cluster(f"{genome_id}_{cluster_idx}", dist)
                cur_cluster.genomeID = genome_id
                cur_cluster.contig = prev_prot.gene_contig
                # Vorheriges (Anker) GEN hinzufügen
                cur_cluster.add_gene(
                    prev_pid,
                    prev_prot.get_domains(),
                    prev_prot.gene_start,
                    prev_prot.gene_end,
                )
                prev_prot.clusterID = cur_cluster.clusterID

            # aktuelles GEN hinzufügen
            cur_cluster.add_gene(pid, p.get_domains(), p.gene_start, p.gene_end)
            p.clusterID = cur_cluster.clusterID
        else:
            # ggf. Cluster abschließen
            if cur_cluster is not None:
                clusters[cur_cluster.clusterID] = cur_cluster
                cur_cluster = None

        # Schiebe Fenster weiter
        prev_pid, prev_prot = pid, p

    # letztes offenes Cluster abschließen
    if cur_cluster is not None:
        clusters[cur_cluster.clusterID] = cur_cluster

    return clusters


def name_syntenic_blocks(
    patterns: Dict[str, Any],
    cluster_id_dict: Dict[str, Cluster],
    min_completeness: float = 0.5,
    collinearity_check: int = 1,
) -> Dict[str, Cluster]:
    """
    3.9.22
    Assigns syntenic pattern keywords to gene clusters based on pattern completeness and (optionally) collinearity.

    For each cluster, all provided patterns are checked. If a minimum completeness (fraction of pattern types found in the cluster)
    is reached, the pattern's keyword is assigned to the cluster along with information about completeness and missing pattern types.
    Optionally, the collinearity (order) of pattern matches in the cluster can be checked.

    Parameters
    patterns : Dict[str, Tuple(set,int)]
        A dictionary mapping pattern IDs to lists of protein types/domains (the pattern).
        Example: { 1: ({'X', 'Y', 'Z'}, 3)} }
    cluster_id_dict : Dict[str, Cluster]
        Dictionary of clusterID to Cluster object. Each Cluster contains genes and protein/domain information.
    min_completeness : float, optional
        Minimum fraction (0–1) of pattern types that must be present in a cluster for it to be annotated with the pattern keyword.
        Default is 0.5.

    Returns
    Dict[str, Cluster]
        The same cluster_id_dict, with each Cluster potentially annotated with one or more pattern keywords (and related info)
        if the patterns are sufficiently complete in the cluster.

    """

    for cluster in cluster_id_dict.values():
        cluster_gene_type_set = cluster.get_protein_type_set()

        for keyword, pattern_tuple in patterns.items():
            pattern_set, length = pattern_tuple
            covered_count = 0

            # Early pruning, count only for completeness
            for gene in pattern_set:
                if gene in cluster_gene_type_set:
                    covered_count += 1

            completeness = covered_count / length
            if completeness < min_completeness:
                continue

            missing_domains = (
                (pattern_set - cluster_gene_type_set)
                if covered_count < length
                else set()
            )
            additional_domains = cluster_gene_type_set.difference(pattern_set)

            # Add the matching proteinIDs for covered types
            covered_ids = cluster.covered_protein_ids
            t2p = cluster.type_to_proteins
            for typ in pattern_set:
                ids = t2p.get(typ)
                if ids:
                    covered_ids.update(ids)

            # Add the keyword to the cluster
            cluster.add_keyword(
                keyword=keyword,
                completeness=completeness,
                csb="0",
                missing=tuple(missing_domains),
                additional_elements=tuple(additional_domains),
                keyword_id=keyword,
            )
            # else:
            #    print("Pattern not recognized")
            #    print(cluster.clusterID)
            #    print(pattern_names[pattern_id])
            #    print(completeness)
            #    print(pattern_set)
            #    print(missing_domains)
    return cluster_id_dict

from typing import Dict, Any, List

def name_syntenic_blocks_trie(
    cluster_dict: Dict[str, Any],
    index: TrieIndex,
    *,
    min_completeness: float = 0.5,
    only_terminal_node_patterns: bool = False,   # True: nur exakte Patterns, die am Endknoten terminieren
    include_partials: bool = True,       # True: zusätzlich auch Teiltreffer (unvollständige Patterns) aufnehmen
) -> dict[str, Any]:
    """
    Benennt Syntenie-Cluster anhand des vorbereiteten TrieIndex.
    - Completeness = k / |P|, wobei k die Tiefe des letzten erreichbaren Trie-Knotens ist.
    - missing = Pattern − Cluster (Domänennamen)
    - additional = Cluster − Pattern (Domänennamen, inkl. unbekannter Domains)

    Effekt: schreibt die Keywords direkt in die Cluster-Objekte via cluster.add_keyword(...).

    Args:
        only_terminal_node_patterns (bool):
    """

    # Über alle Cluster iterieren, die annotiert werden sollen
    for cluster_id, cluster in cluster_dict.items():  # <- alle Cluster einmal durchgehen
        # Alle Domains (als Strings) aus dem Cluster holen
        present_names: set[str] = set(cluster.get_protein_type_set())
        if len(present_names) <= 1:
            continue

        # Cluster-Domains auf IDs der Trie-Ordnung abbilden; zusätzlich Originalnamen behalten
        present_ids_sorted, present_mask, present_names_all, unknown_names = \
            csb_trie_algorithm.map_present_to_sorted_ids(present_names, index)

        # Im Trie so weit wie möglich entlang der vorhandenen Domains absteigen
        best_k: dict[int, int] = sb_trie_matching.iter_candidates_anywhere_start(
            index,
            present_ids_sorted,
            present_mask,

        )

        # Candidate pattern declaration
        if only_terminal_node_patterns and not include_partials:
            # Nur exakte Terminale: |P| == k
            pids_k = {pid: k for pid, k in best_k.items() if index.pattern_meta[pid].length == k}
        else:
            # Return all pattern IDs (exact + incomplete) from the subtree
            pids_k = best_k

        covered_types = set() # needed for cluster object
        # Über alle ausgewählten Pattern-Kandidaten iterieren
        for pid, k in pids_k.items():
            meta = index.pattern_meta[pid]
            if meta.length == 0: # Leere Patterns (sollten praktisch nicht vorkommen) überspringen
                continue

            # Completeness berechnen:
            covered = (meta.mask & present_mask).bit_count()
            completeness = covered / meta.length




            # Pattern-Domänen als NAMEN (nicht IDs)
            pat_names = { index.id_to_domain[i] for i in meta.ids_sorted }

            # Mindestschwelle anwenden – Patterns unterhalb werden nicht hinzugefügt
            if completeness < min_completeness:
                continue
            # Fehlende Domänen = im Pattern, aber nicht im Cluster
            missing = sorted(pat_names - present_names_all)

            # Zusätzliche Domänen = im Cluster, aber nicht im Pattern
            # (enthält bewusst auch unbekannte Domänen, da present_names_all Originalnamen umfasst)
            additional = sorted(present_names_all - pat_names)

            # >>> Direktes Hinzufügen ins Cluster-Objekt (ohne Zwischen-dict) <<<
            covered_types.update(pat_names & present_names_all)

            cluster.add_keyword(
                keyword=meta.name,
                completeness=completeness,
                csb="0",  # falls du hier eine echte CSB-ID hast, setze sie entsprechend
                missing=tuple(missing),
                additional_elements=tuple(additional),
                keyword_id=meta.name,  # oder eine echte ID, falls vorhanden
            )
        # Adds the domain types that are covered by the recognized patterns. Later used to remove low score hits
        # from the cluster not matching a pattern
        cluster.add_covered_types(covered_types)

    return cluster_dict
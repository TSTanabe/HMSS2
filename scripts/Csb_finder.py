#!/usr/bin/python

import re
from multiprocessing import Manager, Pool, Value
from typing import Any, Dict, List, Optional, Set, Tuple

from . import Database
from . import myUtil
logger = myUtil.logger


class Cluster:
    """
    3.9.22
    Cluster organises genes laying in synteny with a specific order. Each cluster has a unique clusterID derived from the assemblyID but with an added index number. 
    
    Args
        clusterID = unique string
        
    11.04.23 added the cluster_start and cluster_end lines   
     
    Organizes genes in synteny with a specific order.
    Each cluster has a unique clusterID derived from the assemblyID plus an index number.
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
        
        self.type_to_proteins: Dict[str, Set[str]] = {} # Saves type => set(proteinID). Necessary to get proteinIDs with pattern matching types
        self.covered_protein_ids: Set[str] = set() # proteinIDs which are part of matching patterns

        
    def add_gene(self, protein_id: str, types: str, start: Optional[int] = None, end: Optional[int] = None) -> None:
        self.genes.append(protein_id)
        self.types.append(types)
        
        if self.cluster_start is None or (start is not None and start < self.cluster_start):
            self.cluster_start = start
        if self.cluster_end is None or (end is not None and end > self.cluster_end):
            self.cluster_end = end

        if types not in self.type_to_proteins:
            # Define the proteinIDs per protein type
            self.type_to_proteins[types] = set()
        self.type_to_proteins[types].add(protein_id)


    def add_keyword(self, keyword: str, completeness: float = 0.0, csb: str = ".", missing: Optional[Set[str]] = None, additional_elements: Optional[Set[str]] = None, keyword_id: Any = ".") -> None:
        missing = missing or []
        keyword_id = keyword if keyword_id == "." else keyword_id
        self.keywords_dict[keyword_id] = Keyword(keyword, completeness, csb, missing, additional_elements, keyword_id)

    def get_keywords(self) -> List["Keyword"]:
        return list(self.keywords_dict.values())
       
    def get_genes(self) -> List[str]:
        return self.genes

    def get_domains(self) -> List[str]:
        return self.types
        
    def get_domain_set(self) -> Set[str]:
        tmp = '-'.join(self.types)
        return set(tmp.split('-'))

    def get_clusterID(self):
        return self.clusterID
        
    def get_protein_type_list(self) -> List[str]:
        """
        Returns a list of protein types behind the last '_' in each entry of self.types.
        Example: "grp0_ProteinA" → "ProteinA", "TGIRFAM0000_ProteinB" → "ProteinB".
        """
        return [typus.split('_')[-1] if '_' in typus else typus for typus in self.types]

    def get_protein_type_set(self) -> Set[str]:
        """
        Returns a set of unique protein types behind the last '_' in each entry of self.types.
        Example: "grp0_ProteinA" → "ProteinA", "TGIRFAM0000_ProteinB" → "ProteinB".
        """
        return {typus.split('_')[-1] if '_' in typus else typus for typus in self.types}


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
        return [self.get_clusterID(), keywords_string]
	
class Keyword:
    """
    3.9.22
    Holds the information about a single keyword, its completeness and if it is a csb.
    """

    def __init__(self, keyword: str, completeness: float = 0.0, csb: str = ".", missing: Optional[Set[str]] = None, additional: Optional[Set[str]] = None, keyword_id: Any = ".") -> None:

        if missing is None:
            missing = set()
        if additional is None:
            additional = set()
        self.keyword: str = str(keyword)
        self.csb: str = csb # Collinear to the reference pattern
        self.completeness: float = completeness
        self.missing_domains: Set[str] = missing
        self.additional_domains: Set[str] = additional
        self.keyword_id: Any = keyword_id

    def get_keyword(self) -> str:
        return self.keyword

    def get_csb(self) -> str:
        return self.csb

    def get_completeness(self) -> float:
        return self.completeness

    def get_missing_domains(self) -> set:
        return self.missing_domains
    
    def get_additional_domains(self) -> set:
        return self.additional_domains

    def set_keyword(self, keyword: str) -> None:
        self.keyword = keyword

    def set_csb(self, csb: str) -> None:
        self.csb = csb

    def set_completeness(self, completeness: float) -> None:
        self.completeness = completeness

def check_order(test_list: List[int]) -> int:
    """
    Checks if a list is ordered ascending or descending.
    """
    if all(test_list[i] <= test_list[i + 1] for i in range(len(test_list)-1)):
        return 1
    elif all(test_list[i] >= test_list[i + 1] for i in range(len(test_list)-1)):
        return 1
    else:
        return 0

def make_pattern_dict(filepath: str) -> Tuple[Dict[int, List[str]], Dict[int, str]]:
    """
    Reads a file containing patterns and their names, returning two dictionaries.
    """
    pattern: Dict[int, List[str]] = {}
    pattern_names: Dict[int, str] = {}
    i = 1

    try:
        with open(filepath, "r") as reader:
            for line_num, line in enumerate(reader, start=1):
                line = line.strip()
                if not line:
                    continue
                try:
                    parts = line.split("\t")
                    name = parts.pop(0)
                    pattern_names[i] = name.strip()
                    patlist = [item.strip() for item in parts if item.strip()]
                    if patlist:
                        pattern[i] = patlist
                        pattern_names[i] = name.strip()
                        i += 1
                except IndexError:
                    logger.warn(f"Skipping unrecognized pattern on line {line_num}: {line}")
    except FileNotFoundError:
        logger.error(f"File not found: {filepath}")
        return {}, {}
    except IOError as e:
        logger.error(f"Unable to read file: {filepath}. {e}")
    return pattern, pattern_names
    
    
def find_syntenic_blocks(
    genome_id: str,
    protein_dict: Dict[str, Any],
    distance: int = 3500
) -> Dict[str, Cluster]:
    """
    3.9.22
    Gets a dictionary with protein objects from one genome. Creates a cluster object and adds proteinID belonging to this genetic cluster forming a syntenic block.
    Order of addition to the syntenic block follows the contig and start order. Two different contigs cannot exist inside a syntenic block. Returns a list of cluster objects
    for the cluster analysis.
    
    Args:
        protein_dict - dictionary with key:proteinID and value:proteinObject
        distance - integer of maximal nucleotides distance between genes to consider in synteny
        genomeID - Assembly identifier for unique cluster id
    Return:
        dictionary of cluster objects
    """
    new_sb = 1
    cluster_id_dict: Dict[str, Cluster] = {}
    protein_id_list = sorted(protein_dict, key=lambda x: (protein_dict[x].gene_contig, protein_dict[x].gene_start))
    cluster_id_number = 1
    cluster = Cluster(f"{genome_id}_{cluster_id_number}", distance)
    cluster.genomeID = genome_id

    for index, elem in enumerate(protein_id_list):
        if index - 1 >= 0:
            prev_el_protein_id = str(protein_id_list[index-1])
            curr_el_protein_id = str(elem)
            prev_protein = protein_dict[prev_el_protein_id]
            curr_protein = protein_dict[curr_el_protein_id]

            if (
                prev_protein.gene_contig == curr_protein.gene_contig and
                curr_protein.gene_start - prev_protein.gene_end <= distance
            ):
                if new_sb:
                    new_sb = 0
                    cluster.add_gene(prev_el_protein_id, prev_protein.get_domains(), prev_protein.gene_start, prev_protein.gene_end)
                    prev_protein.clusterID = cluster.clusterID
                    cluster.add_gene(curr_el_protein_id, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    cluster.contig = prev_protein.gene_contig
                    curr_protein.clusterID = cluster.clusterID
                else:
                    cluster.add_gene(curr_el_protein_id, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    curr_protein.clusterID = cluster.clusterID
            elif new_sb == 0:
                cluster_id_dict[f"{genome_id}_{cluster_id_number}"] = cluster
                cluster_id_number += 1
                cluster = Cluster(f"{genome_id}_{cluster_id_number}", distance)
                cluster.genomeID = genome_id
                new_sb = 1

    if not new_sb:
        cluster_id_dict[f"{genome_id}_{cluster_id_number}"] = cluster
    return cluster_id_dict




def name_syntenic_blocks(
    patterns: Dict[int, List[str]],
    pattern_names: Dict[int, str],
    cluster_id_dict: Dict[str, Cluster],
    min_completeness: float = 0.5,
    collinearity_check: int = 1
) -> Dict[str, Cluster]:
    """
    3.9.22
    Assigns syntenic pattern keywords to gene clusters based on pattern completeness and (optionally) collinearity.

    For each cluster, all provided patterns are checked. If a minimum completeness (fraction of pattern types found in the cluster)
    is reached, the pattern's keyword is assigned to the cluster along with information about completeness and missing pattern types.
    Optionally, the collinearity (order) of pattern matches in the cluster can be checked.

    Parameters
    ----------
    patterns : Dict[int, List[str]]
        A dictionary mapping pattern IDs to lists of protein types/domains (the pattern).
        Example: {1: ['A', 'B', 'C'], 2: ['X', 'Y', 'Z']}
    pattern_names : Dict[int, str]
        A dictionary mapping pattern IDs to pattern names/keywords.
        Example: {1: 'ABC_pattern', 2: 'XYZ_pattern'}
    cluster_id_dict : Dict[str, Cluster]
        Dictionary of clusterID to Cluster object. Each Cluster contains genes and protein/domain information.
    min_completeness : float, optional
        Minimum fraction (0–1) of pattern types that must be present in a cluster for it to be annotated with the pattern keyword.
        Default is 0.5.
    collinearity_check : int, optional
        If 1, collinearity (order of pattern matches in the cluster) is checked and stored with the keyword.
        If 0, collinearity is ignored. Default is 1.

    Returns
    -------
    Dict[str, Cluster]
        The same cluster_id_dict, with each Cluster potentially annotated with one or more pattern keywords (and related info)
        if the patterns are sufficiently complete in the cluster.

    Notes
    -----
    - Each Cluster will have updated keyword assignments in its .keywords_dict.
    - For every pattern that passes the completeness threshold, information about missing pattern types,
      completeness, and (optionally) collinearity is recorded.
    - If you require also the covered proteinIDs, extend Cluster as discussed previously.
    """

    for cluster in cluster_id_dict.values():
        
        protein_type_set = cluster.get_protein_type_set()
        
        for pattern_id, pattern in patterns.items():
            keyword = pattern_names[pattern_id] # Define the name of the pattern
            pattern_set = set(pattern)
            
            additional_elements = protein_type_set.difference(pattern_set)
            
            missing_elements = pattern_set.difference(protein_type_set)
            completeness = (len(pattern_set) - len(missing_elements)) / len(pattern_set)
            
            if min_completeness <= completeness:
                
                covered_types = pattern_set & protein_type_set

                covered_protein_ids = set()
                for typ in covered_types:
                    covered_protein_ids.update(cluster.type_to_proteins.get(typ, set()))

                cluster.covered_protein_ids.update(covered_protein_ids)
                cluster.add_keyword(keyword, completeness, "0", missing_elements, additional_elements, pattern_id)
            #else:
            #    print("Pattern not recognized")
            #    print(cluster.clusterID)
            #    print(pattern_names[pattern_id])
            #    print(completeness)
            #    print(pattern_set)
            #    print(missing_elements)
    return cluster_id_dict
    
    



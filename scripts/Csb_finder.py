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
        self.cluster_id: str = cluster_id
        self.genome_id: str = ""
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


    def add_keyword(self, keyword: str, completeness: float = 0.0, csb: str = ".", missing: Optional[List[str]] = None, keyword_id: Any = ".") -> None:
        missing = missing or []
        keyword_id = keyword if keyword_id == "." else keyword_id
        self.keywords_dict[keyword_id] = Keyword(keyword, completeness, csb, missing, keyword_id)

    def get_keywords(self) -> List["Keyword"]:
        return list(self.keywords_dict.values())
       
    def get_genes(self) -> List[str]:
        return self.genes

    def get_domains(self) -> List[str]:
        return self.types
        
    def get_domain_set(self) -> Set[str]:
        tmp = '-'.join(self.types)
        return set(tmp.split('-'))

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
        return [self.get_cluster_id(), keywords_string]
	
class Keyword:
    """
    3.9.22
    Holds the information about a single keyword, its completeness and if it is a csb.
    """

    def __init__(self, keyword: str, completeness: float = 0.0, csb: str = ".", missing: Optional[List[str]] = None, keyword_id: Any = ".") -> None:
        if missing is None:
            missing = []
        self.keyword: str = str(keyword)
        self.csb: str = csb
        self.completeness: float = completeness
        self.missing_domains: List[str] = missing
        self.keyword_id: Any = keyword_id

    def get_keyword(self) -> str:
        return self.keyword

    def get_csb(self) -> str:
        return self.csb

    def get_completeness(self) -> float:
        return self.completeness

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
                    pattern[i] = [item.strip() for item in parts if item.strip()]
                    i += 1
                except IndexError:
                    logger.error(f"WARNING: Skipping unrecognized pattern on line {line_num}: {line}")
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
    cluster.genome_id = genome_id

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
                    prev_protein.cluster_id = cluster.cluster_id
                    cluster.add_gene(curr_el_protein_id, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    cluster.contig = prev_protein.gene_contig
                    curr_protein.cluster_id = cluster.cluster_id
                else:
                    cluster.add_gene(curr_el_protein_id, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    curr_protein.cluster_id = cluster.cluster_id
            elif new_sb == 0:
                cluster_id_dict[f"{genome_id}_{cluster_id_number}"] = cluster
                cluster_id_number += 1
                cluster = Cluster(f"{genome_id}_{cluster_id_number}", distance)
                cluster.genome_id = genome_id
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
            
            missing_elements = pattern_set.difference(protein_type_set)
            completeness = (len(pattern_set) - len(missing_elements)) / len(pattern_set)
            
            if min_completeness <= completeness:
                
                covered_types = pattern_set & protein_type_set

                covered_protein_ids = set()
                for typ in covered_types:
                    covered_protein_ids.update(cluster.type_to_proteins.get(typ, set()))

                cluster.covered_protein_ids.update(covered_protein_ids)
                cluster.add_keyword(keyword, completeness, "0", list(missing_elements), pattern_id)
    
    return cluster_id_dict
    
    

def parallel_name_syntenic_blocks(options, genomeIDs):
    """
    Parallelized CSB naming process
    """
    global csb_patterns_diction, csb_pattern_names
    
    # Initialize pattern dictionaries
    csb_patterns_diction, csb_pattern_names = makePatternDict(options.patterns_file)

    # Setup multiprocessing resources
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
   
    # Create process pool
    with Pool(processes=options.cores) as pool:
        # Start the writer process asynchronously
        p_writer = pool.apply_async(csb_renaming_writer, (data_queue, options)) 

        # Create tasks for parallel execution
        tasks = ((data_queue, genomeID, options, counter) for genomeID in genomeIDs)
        
        # Execute tasks in parallel
        pool.map(process_parallel_naming, tasks)

        # Signal the end of data to the writer process
        for _ in range(options.cores):
            data_queue.put(None)

        # Wait for the writer process to finish
        p_writer.get()

    print("\nFinished processing all genomes.")

def process_parallel_naming(args_tuple):
    """
    Process a single genome in parallel.
    """
    global csb_patterns_diction, csb_pattern_names
    
    # Unpack the arguments
    queue, genomeID, options, counter = args_tuple

    # Increment counter (shared among processes)
    counter.value += 1
    print(f"Searching assembly {counter.value}", end="\r")

    # Fetch and process cluster dictionary
    cluster_diction = Database.fetch_cluster_dict(options.database_directory, genomeID)
    name_syntenicblocks(csb_patterns_diction, csb_pattern_names, cluster_diction, options.min_completeness)

    # Add the result to the queue
    queue.put((genomeID, cluster_diction))

def csb_renaming_writer(queue, options):
    """
    Writes processed results to the database.
    """
    while True:
        tup = queue.get()
        if tup is None:
            break

        # Unpack and write the result to the database
        genomeID, cluster_dict = tup
        Database.insert_database_clusters(options.database_directory, cluster_dict)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
			
def extract_missing_domain_targets(cluster_dict):
    """
    Creates a domain-centric view of all missing domains across all clusters.
    
    Returns a dictionary:
        {
            domain_name: [
                (genomeID, contig, start, end),
                ...
            ]
        }
    """
    result = {}

    for cluster in cluster_dict.values():

        genomeID = cluster.genomeID
        clusterID = cluster.clusterID
        contig = cluster.contig
        start = cluster.cluster_start
        end = cluster.cluster_end

        for keyword in cluster.get_keywords():
            for missing_domain in keyword.missing_domains:
                if missing_domain not in result:
                    result[missing_domain] = []
                result[missing_domain].append((genomeID, clusterID, contig, start, end))

    return result



def find_csb_pattern_difference(patterns,pattern_names,cluster_dict,min_pattern_length=4):
    #11.04.23
    #part of the synteny supported correction of hits below cutoff
    missing_proteins_list_dict = {} #list of types missing per cluster
    for clusterID,cluster in cluster_dict.items():
        protein_type_list = cluster.get_domain_ends_list() #set of proteins in the cluster
        keywords = cluster.get_keywords() #keywords assigned to the cluster
        missing_proteins = set()
        for keyword in keywords:
            if keyword.get_completeness() < 1: # iterates all incomplete keyword patterns

                name = keyword.keyword
                pattern_id_list = [k for k, v in pattern_names.items() if v == name] #key is the ID, value is the keyword/the name of the pattern. Saves all ID that have the keyword assigned. Can be multiple because different patterns may have the same keyword
                for ID in pattern_id_list:
                    difference = set(patterns[ID]) - set(protein_type_list)
                    missing_proteins.update(difference)
                    
        if missing_proteins: #if not empty there are proteins missing
            missing_proteins_list_dict[cluster.clusterID] = missing_proteins

    
    missing_protein_types = set() #from all clusters these types are missing
    for protein_types in missing_proteins_list_dict.values():
        missing_protein_types.update(protein_types)

    return  missing_protein_types, missing_proteins_list_dict


def synteny_completion(gff3_file,protein_dict,cluster_dict,candidate_protein_dict,missing_proteins_list_dict,difference = 3500):
    #11.04.23
    #part of the synteny supported correction of hits below cutoff
    
    with open(gff3_file,"r") as reader:
        for line in reader.readlines():
            if line.startswith("#"):
                continue
            match = re.search('ID=(cds-){0,1}(\S+?)\W{0,1};',line)
            proteinID = match.group(2) #using the match as proteinID in the redo_csb routine possible?
            if proteinID in candidate_protein_dict.keys():
                
                gff = line.split("\t")
                start = gff[3]
                end = gff[4]
                
                for clusterID,cluster_list in missing_proteins_list_dict.items():
                    
                    if candidate_protein_dict[proteinID][0] in cluster_list: #if proteintype in clusterlist of missing proteins then continue testing
                        cluster = cluster_dict[clusterID]
                        
                        if cluster.cluster_start-difference < int(start) < cluster.cluster_end+difference or cluster.cluster_start-difference < int(end) < cluster.cluster_end+difference:
                            
                        #if yes then the below threshold hit has the correct type and is part of the cluster
                            if proteinID in protein_dict:
                                protein = protein_dict[proteinID]
                                protein.add_domain(candidate_protein_dict[proteinID][0],candidate_protein_dict[proteinID][1],candidate_protein_dict[proteinID][2],candidate_protein_dict[proteinID][3])
                            else:
                                protein_dict[proteinID] = ParseReports.Protein(proteinID,candidate_protein_dict[proteinID][0],candidate_protein_dict[proteinID][1],candidate_protein_dict[proteinID][2],candidate_protein_dict[proteinID][3])
                                protein = protein_dict[proteinID]
                                protein.gene_contig = gff[0]
                                protein.gene_start = gff[3]
                                protein.gene_end = gff[4]
                                protein.gene_strand = gff[6]
                                locustag = ParseReports.getLocustag(line)
                                protein.gene_locustag = locustag
                                protein.clusterID = cluster.clusterID

                                cluster.add_gene(proteinID,candidate_protein_dict[proteinID][0],int(start),int(end)) #this cluster should run again trough the naming routine
                                


#!/usr/bin/python


import os
import heapq
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple, Any

from . import Database
from . import Csb_Mp_Algorithm
from . import myUtil

from scipy.spatial import distance
from sklearn.cluster import AgglomerativeClustering

logger = myUtil.logger

            
# For the clustering of csbs by jaccard and agglomerativeClustering
def csb_prediction(options: Any) -> None:
    """
    Main entry for clustering collinear syntenic blocks (CSBs) using the CSB finder and cluster analysis.

    Args:
        options: Options object. Needs:
            .gene_clusters_file (str)
            .glob_chunks (int)
            .insertions (int)
            .occurence (int)
            .min_csb_size (int)
            .csb_name_prefix (str)
            .csb_name_suffix (str)

    Output:
        options.computed_Instances_dict: Dict with CSB patterns (tuple of genes) as keys.
        options.redundant, options.non_redundant: Paths to redundant/non-redundant cluster files.
        options.redundancy_hash: Hash of cluster redundancy.

    Example Input:
        options.gene_clusters_file = "gene_clusters.txt"
        options.glob_chunks = 10000
    """
    logger.info("Detecting and sorting gene clusters for with csb finder algorithm")

    # Sort the all gene clusters that were detected in the search
    options.gene_clusters_file = sort_by_first_column_and_filter_csb(options.gene_clusters_file, options.glob_chunks, options.min_csb_size, options.max_csb_size, options.max_domain_repeats)
    
    #Finds collinear syntenic blocks with the csb finder algorithm using a printed representation of the clusters.
    options.redundant,options.non_redundant = dereplicate(options.gene_clusters_file) #returns two filepaths, dereplicates identical gene clusters
    
    logger.debug("Initilizing hashes for the Csb match point algorithm")
    options.redundancy_hash = create_redundancy_hash(options.redundant) # value is an integer, number of existing replicates
    gene_clusters = create_gene_cluster_hash(options.non_redundant)
    extend_redundancy_hash(options.non_redundant,options.redundancy_hash) # for all which do not have a redundant gene cluster
    #modified CsbfinderS algorithm
    logger.debug("Initilizing Csb match point algorithm for csb pattern recoginition")
    computed_Instances_dict = Csb_Mp_Algorithm.csb_finderS_matchpoint_algorithm(options.redundancy_hash,gene_clusters,options.insertions,options.occurence) #k insertions und q occurences müssen über die optionen festgelegt werden

    # Combine reverse csbs
    reverse_pairs = find_reverse_pairs(computed_Instances_dict)
    computed_Instances_dict = merge_sets_for_pairs(computed_Instances_dict, reverse_pairs)
    
    # Reduce redundancy in the keys
    options.computed_Instances_dict = csb_collapse_to_longest_pattern(computed_Instances_dict)
    




    
def csb_jaccard(options: Any, jaccard_distance: float) -> Dict[str, Set[str]]:
    """
    Agglomerative clustering of CSBs based on Jaccard similarity.

    Args:
        options: Options object with .computed_Instances_dict, .min_csb_size, etc.
        jaccard_distance (float): Threshold for clustering. 0.2 means at least 80% overlap.

    Returns:
        Dict[str, Set[str]]: Maps CSB cluster names to sets of cluster IDs (gene clusters).

    Example Output:
        {'csb*0*': {'GCID_1', 'GCID_3'}, 'csb*1*': {'GCID_2'}}
    """
    computed_Instances_key_list = csb_Instance_key_list(options.computed_Instances_dict, options.min_csb_size)
    cluster_dict = dict()
    if len(computed_Instances_key_list) > 1:
        matrix = calculate_similarity_matrix_jaccard(computed_Instances_key_list) #matrix of similarity and the corresponding clusterID for each row and column as names
        cluster_dict = hierachy_clustering(matrix, jaccard_distance) # 0.2 means that 80 % have to be the same genes
    elif len(computed_Instances_key_list) == 1:
        cluster_dict[0] = [0]
    else:
        return
    #Sorts the csb keywords to the geneclusters based on the present csb
    csb_gene_cluster_dict, grouped_csb_tuples = csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,options.computed_Instances_dict,options.csb_name_prefix,options.csb_name_suffix)
    write_grouped_csb(options.csb_output_file,grouped_csb_tuples)

    #adds the replicates again for saving
    csb_gene_cluster_dict = replicates(csb_gene_cluster_dict,options.redundancy_hash,options.redundant)
    

    return csb_gene_cluster_dict
    









########################################################################################
################ Secondary subroutines #################################################
################ Subroutines used here #################################################
########################################################################################

def sort_by_first_column_and_filter_csb(
    input_file: str, 
    chunk_size: int = 10000, 
    min_columns: int = None, 
    max_columns: int = None,
    max_repeats: int = None
) -> str:
    """
    Sorts a file by column count and then by the first column's string. Merges using external sort to support large files. Also filters die gene clusters by length and repetitive occurence of genes
    
    Do not copy to filtered list if not:
        - minimal number of unique domains is reached
        - maximal csb size is exceeded
        - any domain is more than max_repeats times repeated
    
    Args:
        input_file (str): File to sort (TSV: first column is genome identifier).
        chunk_size (int): Lines per chunk (default: 10000).
        min_columns (int, optional): Minimum columns (inclusive).
        max_columns (int, optional): Maximum columns (inclusive).
        max_repeats (int, optional): Maximum number of repeats per gene in csb (inclusive)

    Returns:
        str: Path to sorted output file.

    Example:
        input_file = "gene_clusters.txt"
        return "sorted_gene_clusters.txt"
    """
    
    logger.debug(
        f"Sorting csb result file by csb length and alphabet for reproducibility. "
        f"min_columns={min_columns}, max_columns={max_columns}, remove repetitives"
    )
    # Generate the output file path with 'sorted_' prefix
    directory, filename = os.path.split(input_file)
    output_file = os.path.join(directory, f"sorted_{filename}")
    
    temp_files = []

    # Step 1: Read file in chunks, sort each chunk, and write to temporary files
    with open(input_file, 'r') as infile:
        while True:
            lines = []
            for _ in range(chunk_size):
                line = infile.readline()
                if not line:
                    break
                line = line.strip()
                if not line:
                    continue
                columns = line.split()
                
                # The csb are filtered at this point by length and composition criteria
                unique_domains = set(columns[1:])
                n_domains = len(columns) - 1 
                if min_columns is not None and len(unique_domains) < min_columns:
                    continue
                if max_columns is not None and n_domains > max_columns:
                    continue
                if max_repeats is not None and csb_too_many_repeats(columns[1:], max_repeats):
                    continue
                
                lines.append(line)
            if not lines:
                break

            # Sort by number of columns, then by first column string
            lines.sort(key=lambda x: (len(x.split()), x.split()[0]))
            temp_file = f'temp_{len(temp_files)}.txt'
            with open(temp_file, 'w') as f:
                f.write('\n'.join(lines) + '\n')
            temp_files.append(temp_file)


    # Step 2: Merge sorted temporary files
    with open(output_file, 'w') as outfile:
        open_files = [open(temp_file, 'r') for temp_file in temp_files]
        sorted_stream = heapq.merge(*(f for f in open_files), key=lambda x: (len(x.split()), x.split()[0]))
        for line in sorted_stream:
            outfile.write(line)
        
        # Close and remove temporary files
        for f in open_files:
            f.close()
        for temp_file in temp_files:
            os.remove(temp_file)
    return output_file
            
            
            
            
def dereplicate(filepath: str) -> Tuple[str, str]:
    """
    Finds and separates redundant and non-redundant lines in a cluster file.

    Args:
        filepath (str): Path to input cluster file.

    Returns:
        Tuple[str, str]: (redundant_file, non_redundant_file)
    """
    logger.debug("Dereplicate csbs with multiple occurences")
    
    # Dictionary to store lines based on variable columns
    line_dict = {}

    # List to store non-redundant lines
    non_redundant_lines = []

    # Read the file and check for duplicates
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            identifier, *variables = line.split('\t')
            variable_container = tuple(variables)
            if variable_container in line_dict:
                line_dict[variable_container].append(identifier)
            elif reversed(variable_container) in line_dict:
                line_dict[variable_container].append(identifier)
            else:
                line_dict[variable_container] = [identifier]
                non_redundant_lines.append(line)

    # Print the list of redundant identifiers
    redundant_identifiers = []
    for identifiers in line_dict.values():
        if len(identifiers) > 1:
            redundant_identifiers.extend(identifiers)

    # Write redundant identifiers to a file
    redundant_file = os.path.join(os.path.dirname(filepath), 'redundant.txt')
    with open(redundant_file, 'w') as f:
        for identifiers in line_dict.values():
            if len(identifiers) > 1:
                f.write('\t'.join(identifiers)+'\n')
                

    # Write non-redundant lines to a file
    non_redundant_file = os.path.join(os.path.dirname(filepath), 'non-redundant.txt')
    with open(non_redundant_file, 'w') as f:
        f.write('\n'.join(non_redundant_lines))

    return redundant_file, non_redundant_file


def create_redundancy_hash(file_path: str) -> Dict[str, int]:
    """
    Counts the number of columns for each key in a file.

    Args:
        file_path (str): Path to file (TSV: first column is key).

    Returns:
        dict: key (str) -> number of columns (int)
    """
    result_dict = {}  # Initialize an empty dictionary

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                value = len(columns)  # Number of columns is the value
                
                # Store the key and value in the dictionary
                result_dict[key] = value
    return result_dict

def create_gene_cluster_hash(file_path: str) -> Dict[str, List[str]]:
    """
    Maps key to value-list from a cluster file.

    Args:
        file_path (str): TSV: first column is key.

    Returns:
        dict: key (str) -> [value1, value2, ...]
    """
    result_dict = {}  # Initialize an empty dictionary

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                values = columns[1:]  # Values from the other columns as a list

                # Store the key and values in the dictionary
                result_dict[key] = values
                

    return result_dict

def extend_redundancy_hash(filepath: str, redundancy_hash: Dict[str, int]) -> Dict[str, int]:
    """
    Ensures every key from file is present in redundancy_hash.

    Args:
        filepath (str): Path to file (first column = key)
        redundancy_hash (dict): key -> count

    Returns:
        dict: updated redundancy_hash
    """
    with open(filepath, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                if not key in redundancy_hash:
                    redundancy_hash[key] = 1
    return redundancy_hash
    
    
def csb_Instance_key_list(Instance_dict: Dict[Any, Any], threshold: int) -> List[Any]:
    """
    Returns keys of Instance_dict with length >= threshold.

    Args:
        Instance_dict (dict): Dict with tuple keys.
        threshold (int): Minimum key length.

    Returns:
        list: Filtered keys

    Example Output:
        [("geneA", "geneB", "geneC"), ...]
    """
    return [k for k in Instance_dict.keys() if len(k) >= threshold]


def write_grouped_csb(filepath: str, data: Dict[str, List[Any]]) -> None:
    """
    Writes CSB clusters and their member tuples to file.

    Args:
        filepath (str): Output path.
        data (dict): csb_name -> [tuple, tuple, ...]
    """
    # Writes down which csb keyword belongs to which csb in a tab-separated file
    with open(filepath, 'w') as f:
        # Iterate through the dictionary items
        for key, value in data.items():
            # Write the key to the file, followed by a tab
            f.write(f"{key}\t")
            # Iterate through the list of tuples
            for item in value:
                # Write each tuple as a tab-separated string
                f.write(f"{item}\t")
            # Add a newline after each entry
            f.write('\n')            
###########################################################################
# Added 19.10.2024 to specifically cluster the reverse csb

def find_reverse_pairs(my_dict: Dict[Tuple, Any]) -> List[Tuple[Tuple, Tuple]]:
    """
    Finds all key pairs where one is the reverse of the other.

    Args:
        my_dict (dict): dict with tuple keys

    Returns:
        list: [(key1, key2), ...] where key2 == key1[::-1]
    """
    reverse_pairs = []
    seen = set()  # Um Duplikate zu vermeiden
    
    for key in my_dict.keys():
        reversed_key = key[::-1]  # Erzeuge das Revers-Tupel
        
        # Überprüfe, ob das Revers in den Keys ist und das Paar noch nicht überprüft wurde
        if reversed_key in my_dict and reversed_key not in seen:
            reverse_pairs.append((key, reversed_key))
            seen.add(key)  # Füge das Paar zu "gesehenen" hinzu
            seen.add(reversed_key)
    
    return reverse_pairs

def merge_sets_for_pairs(my_dict: Dict[Any, Set[Any]], reverse_pairs: List[Tuple[Any, Any]]) -> Dict[Any, Set[Any]]:
    """
    Merges sets for reverse key pairs and keeps only one key.

    Args:
        my_dict (dict): Original dict (tuple -> set).
        reverse_pairs (list): List of key pairs.

    Returns:
        dict: Updated dict with merged sets.
    """
    # Iteriere über die Reverse-Paare und führe die Sets zusammen
    for key1, key2 in reverse_pairs:
        if key1 in my_dict and key2 in my_dict:
            # Führe die Sets von key1 und key2 zusammen
            my_dict[key1] = my_dict[key1].union(my_dict[key2])
            # Entferne key2 aus dem Dictionary
            del my_dict[key2]
    
    return my_dict

###########################################################################
# Read the gene clusters from the file
def jaccard(set1: Set[Any], set2: Set[Any]) -> float:
    """
    Computes Jaccard similarity between two sets.

    Args:
        set1, set2: Sets of genes.

    Returns:
        float: Jaccard similarity (0..1)

    Example:
        jaccard({'A','B'}, {'B','C'}) -> 0.333...
    """
    if not isinstance(set1, set) or not isinstance(set2,set):
        try:
            set1 = set(set1)
            set2 = set(set2)
        except Exception as e:
            logger.error(f"Error converting to sets for Jaccard: {e}")
            return 0
            
    intersection = len(set1.intersection(set2))
    union = len(set1) + len(set2) - intersection
    jaccard_similarity = intersection / union if union != 0 else 0

    return jaccard_similarity


def calculate_similarity_matrix_jaccard(cluster_columns: List[Set[Any]]) -> Any:
    """
    Calculates pairwise Jaccard similarity matrix for clustering.

    Args:
        cluster_columns (list): Each item is a set (genes per CSB).

    Returns:
        numpy.ndarray: similarity matrix (1 on diagonal, [0,1] elsewhere)
    """
    #cluster_columns needs to be a list of sets, not tuple because jaccard relies on the set datatype
    
    # Convert sets to a 2D binary matrix
    unique_elements = sorted(set().union(*cluster_columns))  # Find all unique elements across the sets
    binary_matrix = [[1 if element in cluster else 0 for element in unique_elements] for cluster in cluster_columns]
    
    # Compute the Jaccard distance (1 - Jaccard similarity)
    jaccard_distances = distance.pdist(binary_matrix, metric='jaccard')
    
    # Convert the distances to a full similarity matrix using squareform
    similarity_matrix = 1 - distance.squareform(jaccard_distances)
    
    # Add 1s to the diagonal to represent self-similarity
    for i in range(len(cluster_columns)):
        similarity_matrix[i, i] = 1.0


    return similarity_matrix

def hierachy_clustering(similarity_matrix: Any, threshold: float) -> Dict[int, List[int]]:
    """
    Performs agglomerative clustering on CSBs.

    Args:
        similarity_matrix: 2D array (numpy.ndarray)
        threshold (float): Dissimilarity cutoff (e.g. 0.2)

    Returns:
        dict: cluster_id -> list of row indices (CSBs)

    Example Output:
        {0: [0,2,5], 1: [1,3]}
    """

    avg_dissim_threshold = threshold

    # Compute the dissimilarity matrix
    dissimilarity_matrix = 1 - similarity_matrix
    #print(dissimilarity_matrix)
    # Perform Agglomerative Clustering with average linkage
    clustering = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage='average', distance_threshold=avg_dissim_threshold)
    
    # Fit the clustering model and obtain the cluster labels
    cluster_labels = clustering.fit_predict(dissimilarity_matrix)

    # Create a dictionary to store the clusters
    clusters = defaultdict(list)
    for i, label in enumerate(cluster_labels):
        clusters[label].append(i)

    #returns a dictionary, key is a number identifying the cluster (eine ganze zahl). value is a list with numbers. these numbers
    #are the cluster labels. this dictionary is pushed forward to the subroutine Clusters.keyword_dictionary
    return clusters

def csb_index_to_gene_clusterID(
    cluster_dict: Dict[int, List[int]],
    computed_Instances_key_list: List[Any],
    computed_Instances_dict: Dict[Any, Set[Any]],
    prefix: str = "csb-",
    suffix: str = "_"
) -> Tuple[Dict[str, List[Any]], Dict[str, List[Any]]]:
    """
    Maps CSB cluster index to gene cluster IDs and lists.

    Args:
        cluster_dict: cluster id -> row indices
        computed_Instances_key_list: List of tuple keys (order = matrix rows)
        computed_Instances_dict: key -> set of gene clusters
        prefix, suffix: CSB cluster label

    Returns:
        dict1: csbID -> [gene cluster IDs]
        dict2: csbID -> [tuples]

    Example Output:
        {'csb-0_': ['GC1','GC2'], ...}, {'csb-1_': [(...), (...)]}
    """
    
    result_dict = dict() # csbID to geneclusters
    result_dict2 = dict() # csbID to csb tuples
    
    for key,indices in cluster_dict.items():
        tuples = [computed_Instances_key_list[i] for i in indices] #list of tuples corresponding to the indices from jaccard clustering
        result_dict2[prefix+str(key)+suffix] = tuples
        for i in tuples:
            for e in computed_Instances_dict[i]: #returns all clusterIDs of the csb 'i'
                if prefix+str(key)+suffix in result_dict.keys():
                    result_dict[prefix+str(key)+suffix].append(e)
                else:
                    result_dict[prefix+str(key)+suffix] = [e]
    return result_dict,result_dict2
    

def replicates(
    csb_gene_cluster_dict: Dict[str, List[Any]],
    redundancy_hash: Dict[str, int],
    filepath_redundant: str
) -> Dict[str, Set[Any]]:
    """
    Expands each csb cluster to include redundant gene clusters.

    Args:
        csb_gene_cluster_dict: csb cluster -> set/list of gene clusters
        redundancy_hash: mapping gene cluster id -> count
        filepath_redundant: File with redundancy info

    Returns:
        dict: csb cluster -> set of gene clusters (incl. redundant)
    """
    
    redundant_dict = dict()
    
    # Read cluster IDs from the redundant file
    with open(filepath_redundant, 'r') as file:
        for line in file:
            clusterIDs = line.strip().split('\t')
            first_clusterID = clusterIDs[0]
            redundant_dict[first_clusterID] = clusterIDs[1:]
    # Iterate through the csb_gene_cluster_dict and merge with redundant_dict
    for key, cluster_ids in csb_gene_cluster_dict.items():
        expanded_clusters = set(cluster_ids)  # Create a copy to avoid modifying the original list
        for cluster_id in cluster_ids:
            if cluster_id in redundant_dict:
                expanded_clusters.update(redundant_dict[cluster_id]) 
        csb_gene_cluster_dict[key] = expanded_clusters #remove doublicates
    return csb_gene_cluster_dict  # Return the updated dictionary


######################################################################
# Tertiary subroutines for filtering the computed_Instances_dict

def csb_remove_repetitives(csb_dict: Dict[Any, Any], min_csb_size: int) -> None:
    """
    Removes CSBs where the pattern (key) is smaller than min_csb_size unique genes.

    Args:
        csb_dict (dict): CSB pattern (tuple) -> value
        min_csb_size (int): Minimum unique gene count in key.
    """
    keys_to_remove = [k for k in csb_dict if len(set(k)) < min_csb_size]
    for k in keys_to_remove:
        del csb_dict[k]


def csb_too_many_repeats(domains: list, max_repeats: int = 3) -> bool:
    """
    Returns True if any domain occurs more than max_repeats times.
    Args:
        domains (list): Domains (z.B. columns[1:])
        max_repeats (int): Allowed repeat count (default: 3)
    Example:
        csb_too_many_repeats(['A','A','A','B'], 3) -> False
        csb_too_many_repeats(['A','A','A','A'], 3) -> True
    """
    domain_counts = Counter(domains)
    return any(count > max_repeats for count in domain_counts.values())


def csb_collapse_to_longest_pattern(input_dict: Dict[Any, Set[Any]]) -> Dict[Any, Set[Any]]:
    """
    For CSBs mapping to same gene sets, keep only the longest pattern key.

    Args:
        input_dict: dict of pattern -> set(gene clusters)

    Returns:
        dict: longest pattern for each value -> set(gene clusters)
    """
    # Step 1: Convert set values to sorted tuples (to make them hashable)
    processed_dict = {key: tuple(sorted(value)) for key, value in input_dict.items()}

    # Step 2: Group keys by their values
    value_to_keys = {}
    for key, value in processed_dict.items():
        if value in value_to_keys:
            value_to_keys[value].append(key)
        else:
            value_to_keys[value] = [key]

    # Step 3: Select the longest key for each value
    filtered_dict = {}
    for value, keys in value_to_keys.items():
        # Find the longest key by tuple length
        longest_key = max(keys, key=len)
        filtered_dict[longest_key] = set(value)  # Convert back to set for output consistency

    return filtered_dict





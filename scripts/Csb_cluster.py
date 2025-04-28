#!/usr/bin/python


import os
import numpy as np
from . import Csb_Mp_Algorithm
from scipy.spatial import distance
from sklearn.cluster import AgglomerativeClustering
from itertools import combinations
from collections import defaultdict

# For the clustering of csbs by jaccard and agglomerativeClustering
def csb_prediction(options):

    if not os.path.isfile(options.gene_clusters_file):
        return    
    
    #Finds collinear syntenic blocks with the csb finder algorithm using a printed representation of the clusters.
    options.redundant,options.non_redundant = dereplicate(options.gene_clusters_file) #returns two filepaths, dereplicates identical gene clusters
    
    options.redundancy_hash = create_redundancy_hash(options.redundant) # value is an integer, number of existing replicates
    gene_clusters = create_gene_cluster_hash(options.non_redundant)
    extend_redundancy_hash(options.non_redundant,options.redundancy_hash) # for all which do not have a redundant gene cluster
    #modified CsbfinderS algorithm
    computed_Instances_dict = Csb_Mp_Algorithm.csb_finderS_matchpoint_algorithm(options.redundancy_hash,gene_clusters,options.insertions,options.occurence) #k insertions und q occurences müssen über die optionen festgelegt werden

    options.computed_Instances_dict = computed_Instances_dict
    

def csb_jaccard(options):
    
    if not options.computed_Instances_dict:
        return {}
    
    #convert the keys in computed_instances_dict into a list   #GCA_030823955_1 was not recognized use as example for debugging
    computed_Instances_key_list = csb_Instance_key_list(options.computed_Instances_dict, options.min_csb_size)
    cluster_dict = dict()
    if len(computed_Instances_key_list) > 1:
        matrix = calculate_similarity_matrix_jaccard(computed_Instances_key_list) #matrix of similarity and the corresponding clusterID for each row and column as names
        cluster_dict = hierachy_clustering(matrix,options.jaccard) # 0.2 means that 80 % have to be the same genes
    elif len(computed_Instances_key_list) == 1:
        cluster_dict[0] = [0]
    else:
        return
    #Sorts the csb keywords to the geneclusters based on the present csb
    csb_gene_cluster_dict, grouped_csb_tuples = csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,options.computed_Instances_dict,options.csb_name_prefix,options.csb_name_suffix)
    #adds the replicates again for saving
    replicates(csb_gene_cluster_dict,options.redundancy_hash,options.redundant)
    
    write_grouped_csb(options.csb_output_file,grouped_csb_tuples)

    return csb_gene_cluster_dict
    








def print_cluster(filepath, cluster_dict):
    with open(filepath, 'a') as writer:
        for cluster in cluster_dict:
            writer.write(cluster.clusterID)  # Write cluster ID (assuming cluster.clusterID is a string)
            domains = '\t'.join(cluster.domains)  # Join the domains list into a tab-separated string
            writer.write("\t" + domains + "\n")  # Write cluster domains followed by a newline

    return filepath


def dereplicate(filepath):
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


def create_redundancy_hash(file_path):
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

def create_gene_cluster_hash(file_path):
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

def extend_redundancy_hash(filepath,redundancy_hash):
    with open(filepath, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                if not key in redundancy_hash:
                    redundancy_hash[key] = 1
    return redundancy_hash
    
    
def csb_Instance_key_list(Instance_dict,threshold):
    listing = list()
    for k in Instance_dict.keys():
        if len(k)>=threshold:
            listing.append(k)

    return listing

def write_grouped_csb(filepath, data):
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
#################################


# Read the gene clusters from the file
def jaccard(set1,set2):
    if not isinstance(set1, set) or not isinstance(set2,set):
        try:
            set1 = set(set1)
            set2 = set(set2)
        except Exception as e:
            print("An error occurred while converting to sets:", e)
            return 0
            
    intersection = len(set1.intersection(set2))
    union = len(set1) + len(set2) - intersection
    jaccard_similarity = intersection / union if union != 0 else 0

    return jaccard_similarity

def calculate_similarity_matrix_jaccard(cluster_columns):
    #cluster_columns needs to be a list of sets, not tuple because jaccard relies on the set datatype
    
    # Determine the number of clusters
    num_clusters = len(cluster_columns) if cluster_columns else 0


    # Calculate Jaccard similarity matrix
    jaccard_generator = (jaccard(row1, row2) for row1, row2 in combinations(cluster_columns, r=2))
    flattened_matrix = np.fromiter(jaccard_generator, dtype=np.float64)

    # since flattened_matrix is the flattened upper triangle of the matrix
    # we need to expand it.
    similarity_matrix = distance.squareform(flattened_matrix)
    similarity_matrix += np.identity(len(cluster_columns))# replacing zeros with ones at the diagonal. 
    

    return similarity_matrix

def hierachy_clustering(similarity_matrix,threshold):

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


def csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,computed_Instances_dict,prefix="csb*",suffix="*"):
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
    

def replicates(csb_gene_cluster_dict, redundancy_hash, filepath_redundant):
    redundant_dict = dict()
    # Read cluster IDs from the redundant file
    with open(filepath_redundant, 'r') as file:
        for line in file:
            clusterIDs = line.strip().split('\t')
            first_clusterID = clusterIDs[0]
            redundant_dict[first_clusterID] = clusterIDs[1:]
    # Iterate through the csb_gene_cluster_dict and merge with redundant_dict
    for key, lst in csb_gene_cluster_dict.items():
        new_list = lst.copy()  # Create a copy to avoid modifying the original list
        for cluster_id in lst:
            if cluster_id in redundant_dict.keys():
                new_list.extend(redundant_dict[cluster_id])
        csb_gene_cluster_dict[key] = new_list
    return csb_gene_cluster_dict  # Return the updated dictionary

        




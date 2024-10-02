#!/usr/bin/python
from collections import defaultdict
import os


import sqlite3
import random

from Bio import SeqIO
from . import ParseReports
from . import myUtil
from operator import add
from itertools import combinations
from math import floor


##################################################################
########## Species tree mapping matrix from diction ##############
##################################################################



def main_binary_dataset(options, directory, taxon_dict, protein_dict, cluster_dict):
    """
    Main function to process species binary datasets for proteins and keywords.

    Args:
        options: Object containing options for fetching proteins and keywords.
        directory: Directory where the output files will be written.
        taxon_dict: Dictionary mapping genomeID to taxonomic lineage.
        protein_dict: Dictionary of protein data.
        cluster_dict: Dictionary of cluster data.
    """
    # Process protein dataset
    protein_presence_dict = fetch_protein_presence(protein_dict)
    protein_binary_matrix_dict, prot_headers = create_presence_absence_matrix(taxon_dict, options.fetch_proteins, protein_presence_dict)
    process_binary_dataset(taxon_dict, protein_binary_matrix_dict, prot_headers, directory, options.dataset_divide_sign, "domains")
    
    # Process keyword dataset
    keywords_presence_dict = fetch_keyword_presence(cluster_dict)
    keywords_binary_matrix_dict, key_headers = create_presence_absence_matrix(taxon_dict, options.fetch_keywords, keywords_presence_dict)
    process_binary_dataset(taxon_dict, keywords_binary_matrix_dict, key_headers, directory, options.dataset_divide_sign, "keywords")

    return
    
def process_binary_dataset(taxon_dict, binary_matrix_dict, headers, directory, divide_sign, suffix=""):
    """
    Helper function to process and write binary datasets and taxon group counts.

    Args:
        taxon_dict: Dictionary mapping genomeID to taxonomic lineage.
        binary_matrix_dict: Dictionary mapping genomeID to binary presence/absence data.
        headers: List of headers for the binary data.
        directory: Directory where the output files will be written.
        dataset_name: Name of the dataset (e.g., "domains", "keywords").
    """
    # Create the iTol binary dataline
    dataline_dict = iTol_lineage_to_binary_dataline(taxon_dict, binary_matrix_dict)
    
    # Write the iTol binary dataset
    iTol_binary_file_dataset(directory, dataline_dict, headers, suffix)
    
    # Count taxon group presence by level
    taxon_group_counts_by_level = count_taxon_group_presence_by_level(taxon_dict, binary_matrix_dict,divide_sign)
    
    # Write the taxon group counts to TSV files
    write_taxon_group_counts_to_tsv(taxon_group_counts_by_level, directory, headers, suffix)


def fetch_protein_presence(protein_dict):
    # Dictionary, das die Domains als Schlüssel enthält und Sets von genomeIDs als Werte
    domain_to_genomeIDs = {}

    # Iteriere über alle Protein-Objekte im protein_dict
    for protein in protein_dict.values():
        # Hole die Domain des aktuellen Proteins (single string)
        domain = protein.get_domains()

        # Hole die genomeID des aktuellen Proteins
        genomeID = protein.genomeID

        # Füge die genomeID dem entsprechenden Set für die Domain hinzu
        if domain not in domain_to_genomeIDs:
            # Initialisiere ein neues Set, wenn die Domain noch nicht im Dictionary existiert
            domain_to_genomeIDs[domain] = set()

        # Füge die genomeID zum Set hinzu
        domain_to_genomeIDs[domain].add(genomeID)

    return domain_to_genomeIDs

def fetch_keyword_presence(cluster_dict):
    # Dictionary, das die Keywords als Schlüssel enthält und Sets von genomeIDs als Werte
    keyword_to_genomeIDs = {}

    # Iteriere über alle Cluster-Objekte im cluster_dict
    for cluster in cluster_dict.values():
        # Hole die Keywords des aktuellen Clusters
        keywords = cluster.get_keywords()

        # Hole die genomeID des aktuellen Clusters
        genomeID = cluster.genomeID
        
        # Füge für jedes Keyword die genomeID dem entsprechenden Set hinzu
        for keyword in keywords:
            key = keyword.get_keyword()
            if key not in keyword_to_genomeIDs:
                # Initialisiere ein neues Set, wenn das Keyword noch nicht im Dictionary existiert
                keyword_to_genomeIDs[key] = set()
            
            # Füge die genomeID zum Set hinzu
            keyword_to_genomeIDs[key].add(genomeID)

    return keyword_to_genomeIDs
    
def create_presence_absence_matrix(taxon_dictionary, order_list, genome_presence_dictionary):
    # Dictionary that will store a tab-separated string for each genomeID
    # The string will represent the presence/absence matrix for each genomeID
    presence_absence_dict = {}

    # Identify the additional elements not covered in order_list
    genome_presence = {key for key in genome_presence_dictionary.keys()}
    additional_elements = sorted(set(genome_presence) - set(order_list))
    #valid_order_list = []
    
    # Iterate over all genomeIDs in the taxon_dictionary
    for genomeID in taxon_dictionary:
        presence_absence_list = []

        # Iterate over the list order to check presence in the specified order
        for element in order_list:
            if element in genome_presence_dictionary:
                #valid_order_list.append(element)
                if genomeID in genome_presence_dictionary[element]:
                    presence_absence_list.append('1')
                else:
                    presence_absence_list.append('0')
                
        # Iterate over the additional elements in a consistent order
        for additional_element in additional_elements:
            if genomeID in genome_presence_dictionary[additional_element]:
                presence_absence_list.append('1')
            else:
                presence_absence_list.append('0')

        # Create the tab-separated string and store it in the dictionary
        presence_absence_dict[genomeID] = presence_absence_list
    header_order = order_list + additional_elements
    
    return presence_absence_dict, header_order

    
    
def iTol_lineage_to_binary_dataline(taxon_dict, dataset_dict):
    # 17.08.24
    # Concatenates the lineage to the dataset line by the genomeID
    # Key is the genomeID, value is the taxonomy + presence/absence
    dataline_dict = {}

    # Iterate over the taxon_dict, assuming taxon_dict and dataset_dict share the same keys
    for genomeID, lineage in taxon_dict.items():
        # Check if genomeID exists in dataset_dict
        if genomeID in dataset_dict:
            # Convert the list in dataset_dict[genomeID] to a tab-separated string
            dataline = '\t'.join(dataset_dict[genomeID])
            # Concatenate the lineage with the dataline and add to dataline_dict
            dataline_dict[genomeID] = f"{lineage}\t{dataline}\n"
        else:
            print(f"Warning: genomeID {genomeID} not found in dataset_dict")
    
    return dataline_dict
    
    
    
def iTol_binary_file_dataset(directory,dataline_dict,headers,name=""):
    """
    04.11.22
        Args:
            dataline_dict    dictionary with limiting genomeIDs+taxonomy => binary dataline lineage
            directory      output file
            headers       column headers        
        Return:
        Output:
            iTol binary dataset as filepath
    
    Writes down the iTol dataset for presence/absence of protein or keyword        
    """
    writer = open(directory+f"/iTol_binary_dataset_{name}.txt","w")    
    
    writer.write("DATASET_BINARY\n")
    writer.write("SEPARATOR TAB\n")    
    writer.write("DATASET_LABEL\tdataset_label\n")
    writer.write("COLOR\t#000088\n")
    #
    #space for optional command prints
    #
    writer.write("LEGEND_TITLE\tlegend_title\n")
    writer.write("LEGEND_SHAPES" + "\t2" * len(headers) + "\n")
    writer.write("LEGEND_COLORS" + "\t#2b2c7c" * len(headers) + "\n")
    writer.write("LEGEND_SHAPE_SCALES\t1\t1\t0.5\n")
    writer.write("HEIGHT_FACTOR\t2\n")
    writer.write("SYMBOL_SPACING\t15\n")
    writer.write("FIELD_SHAPES" + "\t2" * len(headers) + "\n")
    writer.write("FIELD_LABELS\t" + '\t'.join(headers) + "\n")
    writer.write("FIELD_COLORS" + "\t#2b2c7c" * len(headers) + "\n")
    writer.write("DATA\n")
    
    for dataline in dataline_dict.values():
        writer.write(dataline)
        
    writer.close()
    
    return    
    
    
    
########################################################################################    
#################### Presence/absence dataset statistics ###############################
########################################################################################
    
def count_taxon_group_presence_by_level(taxon_dict, presence_absence_dict,divide_sign):
    # Define the taxonomic ranks (adjust this list based on your specific data)
    taxonomic_ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    # Dictionary to store the count of '1's for each taxonomic group, grouped by level
    taxon_group_counts_by_level = defaultdict(lambda: defaultdict(list))

    # Iterate over the genomes and their presence/absence lists
    for genomeID, presence_list in presence_absence_dict.items():
        # Get the full taxonomic lineage for the current genomeID
        
        lineage = taxon_dict[genomeID].split(divide_sign)  # Split the lineage into components
        # Iterate over each taxonomic level in the lineage (e.g., Proteobacteria, Gammaproteobacteria)
        for level, taxon_group in enumerate(lineage[1:]):
            # Map the level to the corresponding taxonomic rank
            taxonomic_rank = taxonomic_ranks[level] if level < len(taxonomic_ranks) else f"Level {level + 1}"
            # Initialize the count list for the taxonomic group at this level if not already present
            if not taxon_group_counts_by_level[taxonomic_rank][taxon_group]:
                taxon_group_counts_by_level[taxonomic_rank][taxon_group] = [0] * len(presence_list)

            # Iterate over the presence/absence list and count the '1's at each index
            for index, presence in enumerate(presence_list):
                if presence == '1':
                    taxon_group_counts_by_level[taxonomic_rank][taxon_group][index] += 1

    return taxon_group_counts_by_level #is a dict with dicts with lists


def write_taxon_group_counts_to_tsv(taxon_group_counts_by_level, output_directory, headers, suffix):
    """
    Writes the taxon group counts, grouped by taxonomic rank, into separate TSV files.
    
    Args:
        taxon_group_counts_by_level: Dictionary containing taxon group counts grouped by taxonomic rank.
        output_directory: Directory where the output TSV files will be written.
        headers: List of headers for the presence/absence columns.
    """
    
    # Iterate over each taxonomic rank (e.g., Phylum, Class, etc.)
    for taxonomic_rank, taxon_group_counts in taxon_group_counts_by_level.items():
        # Create the output file path for this taxonomic rank
        output_file = os.path.join(output_directory, f"Absolute_counts_{suffix}_{taxonomic_rank}.tsv")
        
        # Open the file for writing
        with open(output_file, 'w') as tsv_file:
            # Write the header row (e.g., "Taxon_Group", "Presence_1", "Presence_2", ...)
            tsv_file.write("Taxon\t" + "\t".join(headers) + "\n")
            
            # Iterate over each taxon group and its counts
            for taxon_group, counts in taxon_group_counts.items():
                # Convert the counts to strings and join them with tabs
                counts_str = "\t".join(map(str, counts))
                # Write the taxon group name followed by the counts
                tsv_file.write(f"{taxon_group}\t{counts_str}\n")
        
        print(f"Taxon group counts written to {output_file}")









   
    
    
##################################################################
############## Protein tree mapping iTol datasets ################
##################################################################

    

    
def iTol_taxonomy_range_dataset(directory,taxon_dict):
    """
    14.12.22
        Args:
            taxon_dict    dictionary with limiting genomeIDs => taxonomic lineage
            filepath      output file
            headers       column headers        
        Return:
        Output:
            iTol binary dataset as filepath
            
    """
    writer = open(directory+"/iTol_taxon_range_dataset.txt","w")    
    
    writer.write("DATASET_SYMBOL\n")
    writer.write("SEPARATOR TAB\n")    
    writer.write("DATASET_LABEL\tdataset_label\n")
    writer.write("COLOR\t#000088\n")
    #
    #space for optional command prints
    #
    #writer.write("LEGEND_TITLE\tlegend_title\n")
    #writer.write("LEGEND_SHAPES" + "\t2" * len(headers) + "\n")
    #writer.write("LEGEND_COLORS" + "\t#2b2c7c" * len(headers) + "\n")
    #writer.write("LEGEND_SHAPE_SCALES,1,1,0.5\n")
    writer.write("MAXIMUM_SIZE\t5\n")
    writer.write("DATA\n")

    for genomeID in taxon_dict.keys():
        #dataline = '\t'.join(str(item) for item in dataset_dict[genomeID])
        #dataline_dict[genomeID] = taxon_dict[genomeID] + "\t" + dataline + "\n"
        
        #ID,symbol,size,color,fill,position,label
        string = taxon_dict[genomeID] + "\t2" + "\t10" + "\t#ff0000" + "\t1" + "\t1" + "\n"
        writer.write(string)
        
    writer.close()
    return
        
def load_dataline_combination(filepath):
    """
    18.11.22
        Args:
            filepath      output file
        Return:
            all_combined_headers List
        Output:
            
    """
    all_combined_headers = []
    with open(filepath,"r") as reader:
        for line in reader.readlines():
            line = line.replace("\n","")
            line = line.replace(" ","")
            listing = line.split("\t")
            all_combined_headers.extend(tuple(listing))
    return all_combined_headers
    
def dataline_combinations(dataline_dict,headers,filepath, combine=5, min_combine = 1):
    """
    05.11.22
        Args:
            dataline_dict    dictionary with limiting genomeIDs => binary data list
            filepath      output file
            headers       column headers        
        Return:

        Output:
            
    """
    
    #combine all headers available
    all_combined_headers = []
    if combine > 0:
        for i in range(min_combine,combine):
                combined_headers = list(combinations(headers,i+1)) #combinations creates unique permutations by 2 to all elements
                #print(combined_headers)
                all_combined_headers.extend(combined_headers)
    elif myUtil.file_path(filepath):
        all_combined_headers = load_dataline_combination(filepath)
#    combined_headers2 = list(combinations(headers,2)) #combinations creates unique permutations by 2 elements
#    combined_headers3 = list(combinations(headers,3))
#    combined_headers = [*combined_headers2,*combined_headers3]
    
    for genomeID,dataline in dataline_dict.items():
        res = dict(zip(headers, dataline)) #mapping the headline to the value for the dataline of a single genome
        
        listing = []
        for combination in all_combined_headers:
            summe = 0
            for element in combination: #combination is a tuple
                summe += res[element] 
            
            binary = summe/len(combination) #if all present s/len = 1 otherwise < 1
            binary = floor(binary)
            listing.append(binary)
        
                    
        dataline_dict[genomeID] = [*dataline,*listing]

    for combination in all_combined_headers:
        string = '+'.join(combination)
        headers.append(string)    
        
    return dataline_dict,headers

    
def dataset_statistics(database,directory,headers,dataline_dict):
    """
    4.11.22
        database    database of HMSSS
        directory   output directory
        headers     column headers from iTol dataset
        dataline_dict   binary columns from iTol dataset key genomeID+taxonomy => binary data presence absence
        
        output of absolute and relative values for presence absence of fetched protein inside a taxonomic level
    """
    lineage = ["Superkingdom","Clade","Phylum","Class","Ordnung","Family","Genus","Species"]
    
    
    
    
    for lin in lineage:
        taxon_sum = {} #hold overall sum from this taxonomy level Taxon => count list
        taxon_id_count = {} #hold number of genomes for specific taxonomy group
        genomeID_dict = {} #for a taxon lvl gather dict genomeID => Taxon (iterative per lvl)
        order_of_taxonomy = []
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            query = f"SELECT genomeID,{lin} from Genomes ORDER BY Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species"
            cur.execute(query)
            
            
            for row in cur: #relation of genomeID to each taxonomy group
                order_of_taxonomy.append(row[1])
                genomeID_dict[row[0]] = row[1]
                taxon_id_count[row[1]] = 0
                if not row[1] in taxon_sum:
                    taxon_sum[row[1]] = [0] * len(headers)

        for igenomeID in dataline_dict.keys():
            taxon = genomeID_dict[igenomeID] #name of the phylum, or class, etc.
            summery = taxon_sum[taxon]
            dataline = dataline_dict[igenomeID]
            taxon_sum[taxon] = list( map(add, dataline, summery) )
            taxon_id_count[taxon] += 1
        
        #Write absolute values
        writer = open(directory+"/"+lin+"_absolute_count.csv","w")
        order_of_taxonomy = list(dict.fromkeys(order_of_taxonomy)) #makes each lineage unique
        writer.write("\t\t"+"\t".join(headers)+"\n")
        for taxon in order_of_taxonomy:
                writer.write(str(taxon_id_count[taxon])+"\t")
                writer.write(str(taxon)+"\t")
                writer.write('\t'.join(str(item) for item in taxon_sum[taxon]))
                writer.write("\n")
        writer.close()
        
        #Write relative values
        writer = open(directory+"/"+lin+"_relative_count.csv","w")
        order_of_taxonomy = list(dict.fromkeys(order_of_taxonomy)) #makes each lineage unique
        writer.write("\t\t"+"\t".join(headers)+"\n")
        for taxon in order_of_taxonomy:
                id_count = taxon_id_count[taxon]
                writer.write(str(taxon_id_count[taxon])+"\t")
                writer.write(str(taxon)+"\t")
                try:
                    writer.write('\t'.join(str(int(item)/int(id_count)) for item in taxon_sum[taxon]))
                except:
                    #divison with 0 just skip writing down the value
                    #this can happen if a protein has no hits
                    #print(f"{taxon} has 0 hits for dataset",end="\r")
                    pass
                writer.write("\n")

        writer.close()
        #print("\tFinished dataset statistics")
    
def iTol_range_dataset(directory,database,filepath,trennzeichen=';'):
    """
    11.03.23
        Args:
            directory     output directory
            filepath      input file
            trennzeichen  connector between lineages
            database      database sequences were derived from
        Return:
        Input:
            has to be a protein fasta file fetched from the database
            taxonomy must not be assigned
        Output:
            iTol range dataset with unique color per domain type
            headers are written as in binary dataset and add_taxonomy
            routine to match the possible nodes in iTol
            taxonomy is assigned during the process to the dataset, so it
            will be compatible with the other datasets
    """
    print("Creating range dataset for coloring protein types")
    record_dict = {}
    taxon_dict = {}
    color_dict = {}
    color_iterator = 0
    
    
    writer = open(directory+"/iTol_taxon_range_dataset.txt","w")    
    
    writer.write("DATASET_SYMBOL\n")
    writer.write("SEPARATOR TAB\n")    
    writer.write("DATASET_LABEL\tdataset_label\n")
    writer.write("COLOR\t#000088\n")

    writer.write("MAXIMUM_SIZE\t5\n")
    writer.write("DATA\n")
    
    
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()
    
        for record in SeqIO.parse(filepath, "fasta"):
            dataset_range_line = ""
            genomeID, proteinID = record.id.split("-", maxsplit=1)
            genomeID = myUtil.getGenomeID(genomeID)
            #check genome has multiple proteins in the fasta
            if genomeID not in record_dict:
                record_dict[genomeID] = 1
            else:
                record_dict[genomeID] = record_dict[genomeID] + 1

            #Fetch taxonomy
            #Creating the first column for the dataset
            cur.execute("""SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ?""",(genomeID,))
            row = cur.fetchone()
            lineage = myUtil.taxonomy_lineage(row, trennzeichen)
            dataset_range_line = lineage if record_dict[genomeID] == 1 else f"{lineage}{trennzeichen}{record_dict[genomeID]}"
            
            
            #Fetch protein domains
            #Creating the actual range dataset infos
            cur.execute("""SELECT domain FROM Domains WHERE proteinID = ?""",(proteinID,))
            domain_string = ""
            for row in cur:
                result_str += row[0] + "-"
            domain_string = domain_string[:-1]
            
            color = ""
            if label in color_dict:
                color = color_dict[label]
            else:
                color = myUtil.generate_color(color_iterator)
                color_iterator = color_iterator + 42
                color_dict[label] = color
            dataset_range_line = dataset_range_line + "\t2" + "\t10" + f"\t{color}" + "\t1" + "\t1" + "\n" #add to line the color range
            writer.write(dataset_range_line+"\n")
    
    writer.close()
    
    
    
def iTol_domain_dataset(directory,database,filepath,trennzeichen='_'):
    """
    11.03.23
        Args:
            directory     output directory
            filepath      input file
            trennzeichen  connector between lineages
            database      database sequences were derived from
        Return:
        Input:
            has to be a protein fasta file fetched from the database
            taxonomy must not be assigned
        Output:
            iTol protein domain dataset with unique color per domain type in
            a sorted gene cluster
            headers are written as in binary dataset and add_taxonomy
            routine to match the possible nodes in iTol
            taxonomy is assigned during the process to the dataset, so it
            will be compatible with the other datasets
            
    """
    print("Creating range dataset for coloring protein types")
    record_dict = {}
    taxon_dict = {}
    color_dict = {}
    color_iterator = 42
    
    
    writer = open(directory+"/iTol_taxon_domain_dataset.txt","w")    
    
    writer.write("DATASET_DOMAINS\n")
    writer.write("SEPARATOR TAB\n")    
    writer.write("DATASET_LABEL\tdataset_label\n")
    writer.write("COLOR\t#000088\n")
    writer.write("DATA\n")
    
    
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()
    
        for record in SeqIO.parse(filepath, "fasta"):
            dataset_range_line = ""
            genomeID, proteinID = record.id.split("-", maxsplit=1)
            genomeID = myUtil.getGenomeID(genomeID)
            #check genome has multiple proteins in the fasta
            if genomeID not in record_dict:
                record_dict[genomeID] = 1
            else:
                record_dict[genomeID] = record_dict[genomeID] + 1

            #Fetch taxonomy
            #Creating the first column for the dataset
            cur.execute("""SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ?""",(genomeID,))
            row = cur.fetchone()
            lineage = myUtil.taxonomy_lineage(row, trennzeichen)
            dataset_range_line = lineage if record_dict[genomeID] == 1 else str(lineage)+str(trennzeichen)+str(record_dict[genomeID])
            dataset_range_line = lineage if lineage else record.id
            if record_dict[genomeID] > 1:
                dataset_range_line += '_'+str(record_dict[genomeID])
            
            #Fetch protein domains
            #Creating the actual domain dataset infos
            #1st Column identifier| overall range| SHAPE|START|END|COLOR|LABEL
            #node,1200,RE|100|150|#ff0000|SH2,EL|400|500|#0000ff|SH3,OC|700|900|#00ff00|PH
            cur.execute("""SELECT clusterID from Proteins WHERE proteinID = ?""",(proteinID,))
            row = cur.fetchone()
            if not row:
                continue
            clusterID = row[0]
            
            
            
            protein_dict = dict()
            fusion_protIDs = dict()
            query = "SELECT DISTINCT Proteins.proteinID,Proteins.genomeID,Proteins.clusterID,contig,start,end,strand,sequence,domain,domStart,domEnd,score,dom_count from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE Proteins.clusterID = ?"
            cur.execute(query,[clusterID,])
            #Fill the protein dict
            
            for row in cur:
                # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => contig,
                # 4 => start, 5 => end, 6 => strand, 7 => sequence,
                # 8 => domain, 9 => domStart, 10 => domEnd, 11 => score,
                if row[0] in protein_dict.keys():
                    protein = protein_dict[row[0]]
                    protein.add_domain(row[8],row[9],row[10],row[11])
                else:
                    protein = ParseReports.Protein(row[0],row[8],row[9],row[10],row[11])
                    protein.set_genomeID(row[1])
                    protein.set_clusterID(row[2])
                    protein.set_gene_contig(row[3])
                    protein.set_gene_start(row[4])
                    protein.set_gene_end(row[5])
                    protein.set_gene_strand(row[6])
                    protein_dict[row[0]] = protein

                    
            if not protein_dict:
                continue                        
            proteinID_list = sorted(protein_dict, key=lambda x: protein_dict[x].gene_start) # sort proteins into correct order
            
            start_coordinate = protein_dict[proteinID_list[0]].gene_start
            end_coordinate = protein_dict[proteinID_list[-1]].gene_end
            genecluster_range = end_coordinate - start_coordinate
            dataset_range_line += f"\t{genecluster_range}"
            #1st Column identifier| overall range| SHAPE|START|END|COLOR|LABEL
            #node,1200,RE|100|150|#ff0000|SH2,EL|400|500|#0000ff|SH3,OC|700|900|#00ff00|PH
            
            for current_proteinID in proteinID_list:
                protein = protein_dict[current_proteinID]
                
                start = protein.get_gene_start()- start_coordinate
                end = protein.get_gene_end()- start_coordinate
                label = protein.get_domains()
                color = ""
                if label in color_dict:
                    color = color_dict[label]
                else:
                    color = myUtil.generate_color(color_iterator)
                    color_iterator = color_iterator + 42
                    color_dict[label] = color
                dataset_range_line += f"\tRE|{start}|{end}|{color}|{label}"
            
            writer.write(dataset_range_line+"\n")
    writer.close() 
            
            
            
            
            
            
            
            
            
            
    
    
    


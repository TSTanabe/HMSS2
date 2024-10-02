#!/usr/bin/python
import os
import sys
import csv
import sqlite3

from Bio import SeqIO

from . import myUtil
from . import ParseReports
from . import Csb_finder


#########################################################################
####################### MAIN OUTPUT ROUTINE #############################
#########################################################################

def print_fasta_and_hit_table(directory,options):
    print_command_line_args(directory+"1_fetch_command.txt")
    
    path1 = directory+"1_hit_table.txt" #individual hit table in tsv file
    path2 = directory+"2_Protein"       #protein fasta files
    
    #Limit the genomeIDs and fetch the lineage
    taxon_diction = dict()
    if options.limiter:
        taxon_diction = fetch_limiter_data(options)
    #Fetch protein/cluster diction
    if options.fetch_csbs:
        print(f"Collecting gene clusters containing {options.fetch_csbs}")
        csb_listing = find_csbs_with_proteins(options.csb_output_file, options.fetch_csbs)
        print(csb_listing)
        print(f"Found {len(csb_listing)} types of gene clusters containing one the genes for {options.fetch_csbs}")
        options.fetch_keywords.extend(csb_listing)
    

    print("Collecting protein sequences")
    protein_diction,cluster_diction,taxon_diction = fetch_bulk_data(options.database_directory,\
                                                             options.fetch_genomes,options.fetch_proteins,options.fetch_keywords,\
                                                             taxon_diction,options.min_completeness,options.dataset_divide_sign)
    print("Collecting proteins sequences -- ok\n")
    
    print("Writing fasta output to disc\n")
    output_genome_report(path1,protein_diction,cluster_diction,taxon_diction)    #Output report
    files = output_distinct_fasta_reports(path2,protein_diction,cluster_diction) #Output per domain
    singletons(directory,files)                                                  #Output singleton per genome and doublicates
    
    myUtil.clean_empty_files(directory)
    generate_taxonomy_summary(path1,directory+"1_hit_taxonomy_stats.txt")
    print(f"\nWrote fasta output to :\n{directory}")
    
    return protein_diction,cluster_diction,taxon_diction   
    
#########################################################################
#########################################################################
#########################################################################






def print_command_line_args(output_file):
    """
    Print the command-line arguments to a specified file.

    Args:
        output_file: Path to the file where the arguments will be written.
    """
    try:
        with open(output_file, 'w') as f:
            f.write("Command-line arguments passed to the script:\n")
            for index, arg in enumerate(sys.argv):
                f.write(f"Argument {index}: {arg}\n")
    except Exception as e:
        print(f"Failed to write to file {output_file}: {e}")



def print_file_content(file_path):
    """
    Prints the content of the file along with the file path.

    :param file_path: Path to the file
    
    Needed to print csb patterns and other files to terminal
    """
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        print(f"File Path: {file_path}")
        print("File Content:")
        print(content)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")


def get_protein_ids_by_domains(database, domains, keywords=None):
    """
    Retrieve proteinIDs for given domains, with optional filters for keywords.

    Args:
        database (str): Pathway to the database file.
        domains (list of str): List of domains to filter the proteins by.
        keywords (list of str, optional): List of keywords to filter by (connected by 'OR').

    Returns:
        dict: Dictionary with domains as keys and sets of proteinIDs as values.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Base query
        query = """
        SELECT DISTINCT Domains.domain, Proteins.proteinID
        FROM Proteins
        JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE Domains.domain IN ({})
        """.format(','.join(['?'] * len(domains)))

        # Parameters list
        params = domains

        # Add keyword filter if provided
        if keywords:
            query += """
            AND Proteins.clusterID IN (
                SELECT clusterID
                FROM Keywords
                WHERE 
            """
            keyword_conditions = " OR ".join(["Keywords.keyword = ?" for _ in keywords])
            query += f"{keyword_conditions})"
            params.extend(keywords)

        cur.execute(query, params)
        result = cur.fetchall()

    # Construct the dictionary
    domain_protein_dict = {domain: set() for domain in domains}
    for row in result:
        dom, prot_id = row
        domain_protein_dict[dom].add(prot_id)

    return domain_protein_dict


def get_protein_ids_by_keywords(database, keywords, domain=None):
    """
    Retrieve proteinIDs for a given set of keywords, with an optional filter for a domain.

    Args:
        database (str): Pathway to the database file.
        keywords (list of str): List of keywords to filter by (connected by 'OR').
        domain (str, optional): Domain to further filter the proteins by.

    Returns:
        dict: Dictionary with domains as keys and sets of proteinIDs as values.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Base query
        query = """
        SELECT DISTINCT Domains.domain, Proteins.proteinID
        FROM Proteins
        JOIN Clusters ON Proteins.clusterID = Clusters.clusterID
        JOIN Keywords ON Clusters.clusterID = Keywords.clusterID
        JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE 
        """
        
        # Add keyword filter
        keyword_conditions = " OR ".join(["Keywords.keyword = ?" for _ in keywords])
        query += f"{keyword_conditions}"

        # Parameters list
        params = keywords

        # Add domain filter if provided
        if domain:
            query += " AND Domains.domain = ?"
            params.append(domain)

        cur.execute(query, params)
        result = cur.fetchall()

    # Construct the dictionary
    domain_protein_dict = {}
    for row in result:
        dom, prot_id = row
        if dom not in domain_protein_dict:
            domain_protein_dict[dom] = set()
        domain_protein_dict[dom].add(prot_id)

    return domain_protein_dict

# Example usage of routine get_protein_ids_by_keywords:
# database = "path/to/your/database.db"
# keywords = ["keyword1", "keyword2"]
# domain = "example_domain"  # Optional
# domain_protein_dict = get_protein_ids_by_keywords(database, keywords, domain)
# print(domain_protein_dict)


def fetch_protein_details(database, protein_ids):
    """
    Fetches detailed information for a set of proteinIDs.

    Args:
        database (str): Pathway to the database file.
        protein_ids (set): Set of proteinIDs to fetch details for.

    Returns:
        list: List of tuples with detailed protein information.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Construct query to fetch detailed information for given proteinIDs
        query = """
        SELECT 
            Proteins.proteinID, Proteins.genomeID, Proteins.clusterID, 
            Proteins.contig, Proteins.start, Proteins.end, Proteins.strand, 
            Proteins.sequence, Domains.domain, Domains.domStart, 
            Domains.domEnd, Domains.score
        FROM Proteins
        LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE Proteins.proteinID IN ({})
        """.format(','.join(['?'] * len(protein_ids)))

        cur.execute(query, list(protein_ids))
        rows = cur.fetchall()

    return rows


def write_detail_to_protein_dict(details):
    protein_dict = dict()
    cluster_dict = dict()
    count = 1
    for row in details:
        print(f"\tFetched domains: {count}",end="\r")
        count = count +1 
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
            protein.set_protein_sequence(row[7])
            protein_dict[row[0]] = protein

        if not row[2] is None and not row[2] in cluster_dict:
            cluster = Csb_finder.Cluster(row[2])
            cluster.add_gene(row[0],row[8])
            cluster_dict[row[2]] = cluster

    return protein_dict,cluster_dict

#############################################################################################
#############################################################################################
#############################################################################################

def output_genome_report(filepath,protein_dict,cluster_dict,taxon_dict={},genomeID="",writemode="w"):
    """
    11.9.22
    
    Args:
        Filepath    Outputfile
        protein_dict    proteinID => proteinObj
        cluster_dict    key:clusterID => value:clusterObj
        Writemode   default "w" for overwrite/new file otherwise "a" for append
    Output:
        TSV file with columns
        Protein attributes:
        proteinID get_domains domain_scores domain_coordinates
        gene_contig gene_start gene_end gene_strand gene_locustag
        Cluster attributes:
        keyword completeness csb
        Taxonomy lineage 
        
        Reports in a textfile a abbreviated summary of the protein information with cluster and lineage
        
                listing.append(self.proteinID)
        listing.append(self.get_domains())
        listing.append(str (self.get_domain_scores()))
        listing.append(str (self.get_domain_coordinates()))
        listing.append(self.gene_contig)
        listing.append(str (self.gene_start))
        listing.append(str (self.gene_end))
        listing.append(self.gene_strand)
        listing.append(self.gene_locustag)
    """
    header ='\t'.join(["genomeID","proteinID","domains","dom_scores","dom_coords","contig","gene_start","gene_end","gene_strand","locustag","clusterID","keyword","cluster_completeness","","Taxonomy"]) 
    
    with open(filepath, writemode) as writer:
    
        writer.write(header+"\n")
        proteinID_list = sorted(protein_dict, key=lambda x: \
        (protein_dict[x].genomeID, protein_dict[x].gene_contig, protein_dict[x].gene_start)) 
        
        for proteinID in proteinID_list:
        #Write each line includes all information about one protein
            protein = protein_dict[proteinID]
            proteinlist = protein.get_protein_list()  #list representation of a protein object
            
            out = protein.genomeID + '\t' +'\t'.join(proteinlist)
            
            clusterID = protein.clusterID
            if clusterID:
                
                if clusterID in cluster_dict.keys():
                    cluster = cluster_dict[clusterID]
                    clusterlist = cluster.get_cluster_list(",")
                    out += '\t' + '\t'.join(clusterlist)
                
                else:
                    out += '\t'+ clusterID + '\t\t\t'
            else:
                out += '\t\t\t\t'
            
            genomeID = protein.genomeID
            if genomeID in taxon_dict.keys():
                taxon = taxon_dict[genomeID]
                out += '\t' + taxon
            else:
                out += '\t\t\t\t\t\t\t\t\t\t\t'
            out += '\n'                
            writer.write(out)
    return


































##############################################################
########## Fetch information from database routines ##########
##############################################################

def fetch_limiter_data(options):
    """
    Fetch limiter data from the database based on specified conditions.

    Args:
        database (str): Name of the database to be worked on.
        lineage (str): Limits to specific lineage.
        taxon (str): Limits to specific taxon.
        proteins (list): Limits to specific protein type.
        keywords (list): Limits to specific cluster keyword.
        min_cluster_completeness (float): Minimal completeness for cluster keywords to occur in the output.
        trennzeichen (str): Separator for taxonomy information.

    Returns:
        dict: Taxon dictionary mapping genomeID to taxonomy lineage.
    """
    database = options.database_directory
    lineage = options.dataset_limit_lineage
    taxon = options.dataset_limit_taxon
    proteins = options.dataset_limit_proteins
    keywords = options.dataset_limit_keywords
    trennzeichen = options.dataset_divide_sign
    
    
    taxon_dict = {}

    query = "SELECT DISTINCT Genomes.genomeID, Superkingdom, Phylum, Class, Ordnung, Family, Genus, Species FROM Genomes"
    conditions = []
    params = []

    if lineage and taxon:
        conditions.append(f"{lineage} LIKE ?")
        params.append(f"%{taxon}%")
        print(f"\tLimiting to taxonomy {lineage} like {taxon}")

    if proteins:
        protein_conditions = " OR ".join(["domain LIKE ?"] * len(proteins))
        query += " LEFT JOIN Domains ON Genomes.genomeID = Domains.genomeID"
        conditions.append(f"({protein_conditions})")
        params.extend([f"%{protein}%" for protein in proteins])
        print(f"\tLimiting to proteins {proteins}")

    if keywords:
        keyword_conditions = " OR ".join(["keyword LIKE ?"] * len(keywords))
        query += " LEFT JOIN Keywords ON Genomes.genomeID = Keywords.genomeID"
        conditions.append(f"({keyword_conditions})")
        params.extend([f"%{keyword}%" for keyword in keywords])
        print(f"\tLimiting to keywords {keywords}")

    if conditions:
        query += " WHERE " + " AND ".join(conditions)
    
    print("Collecting taxonomy lineage information")

    with sqlite3.connect(database) as con:
        con.execute("PRAGMA foreign_keys = ON;")
        cur = con.cursor()
        cur.execute(query, params)
        for index, row in enumerate(cur):
            print(f"\tSelecting genomes {index + 1}", end="\r")
            if row[0] not in taxon_dict:
                taxon_dict[row[0]] = myUtil.taxonomy_lineage(row, trennzeichen)
    
    print("\nCollecting taxonomy lineage information -- ok\n")
    
    return taxon_dict


def parse_protein_set(protein_set_str):
    """
    Parse a string representation of a protein set into a Python list.
    
    Args:
        protein_set_str (str): String representation of a protein set.
        
    Returns:
        list: List of proteins in the set.
    """
    protein_set_str = protein_set_str.strip("()")
    if not protein_set_str:
        return []
    return [protein.strip().strip("'") for protein in protein_set_str.split(",")]

def find_csbs_with_proteins(file_path, proteins):
    """
    Find all CSB identifiers that contain all given proteins.

    Args:
        file_path (str): Path to the file containing CSB data which is the "./project/Collinear_syntenic_blocks/Csb_output.txt" file.
        proteins (list): List of proteins to search for.

    Returns:
        list: List of CSB identifiers that contain all given proteins.
    """
    csb_list = []
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # skip empty lines
                parts = line.split("\t")
                csb_id = parts[0]
                protein_sets = [parse_protein_set(ps) for ps in parts[1:]]
                if all(any(protein in protein_set for protein_set in protein_sets) for protein in proteins):
                    csb_list.append(csb_id)

    return sorted(csb_list)
    
    
    
def fetch_bulk_data(database, genomes, proteins, keywords, taxon_dict=dict(), min_cluster_completeness=0, trennzeichen=';'):
    """
    Fetch bulk data from the database based on specified conditions.

    Args:
        database (str): Name of the database to be worked on.
        proteins (list): Limits to specific protein type.
        keywords (list): Limits to specific cluster keyword.
        taxon_dict (dict): Dictionary to store taxonomy information.
        min_cluster_completeness (float): Minimal completeness for cluster keywords to occur in the output.
        trennzeichen (str): Separator for taxonomy information.

    Returns:
        tuple: protein dictionary with key:proteinID => value:protein object for a single genome,
               cluster dictionary with key:clusterID => value:cluster object,
               updated taxon dictionary genomeID => taxonomy lineage.
    """
    
    protein_dict = {}
    cluster_dict = {}
    genomeID_set = set()
    fusion_protIDs = set()

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")
        query, args = generate_fetch_query(genomes, proteins, keywords, taxon_dict) # inner connection is or connection between prot, key and tax is and
        cur.execute(query, args)
        # Fetched only proteins with doms in keywords and genomeID in taxon dict, if taxon dict is present
        for index, row in enumerate(cur):
            print(f"\tFetched domains: {index + 1}", end="\r")#
            proteinID, genomeID, clusterID = row[0], row[1], row[2]
            domain, domStart, domEnd, score = row[8], row[9], row[10], row[11]

            if proteinID in protein_dict:
                protein_dict[proteinID].add_domain(domain, domStart, domEnd, score)
            else:
                protein = ParseReports.Protein(proteinID, domain, domStart, domEnd, score)
                protein.genomeID = genomeID
                protein.clusterID = clusterID
                protein.gene_contig = row[3]
                protein.gene_start = row[4]
                protein.gene_end = row[5]
                protein.gene_strand = row[6]
                protein.protein_sequence = row[7]
                protein_dict[proteinID] = protein
                genomeID_set.add(genomeID)

            if clusterID is not None and clusterID not in cluster_dict.keys(): # only needed to know which clusters are involved
                cluster = Csb_finder.Cluster(clusterID)
                cluster.genomeID = genomeID
                cluster.add_gene(proteinID, domain)
                cluster_dict[clusterID] = cluster
                
            if row[12] >= 2:
                fusion_protIDs.add(proteinID)
        print("")
        
        # Add the additional domains to fused proteins with multiple domains
        if fusion_protIDs:
            # Convert fusion_protIDs to a tuple for use in the SQL IN clause
            placeholders = ','.join(['?'] * len(fusion_protIDs))
            query = f"SELECT DISTINCT proteinID, domain, domStart, domEnd, score FROM Domains WHERE proteinID IN ({placeholders})"
            
            # Execute the query with all proteinIDs at once
            cur.execute(query,tuple(fusion_protIDs))
            
            # Fetch all results at once
            rows = cur.fetchall()
            
            # Process the fetched rows
            for index, row in enumerate(rows):
                print(f"\tFetched fused domains: {index} proteinID {row[0]}", end="\r")
                protein_dict[row[0]].add_domain(row[1], row[2], row[3], row[4])
            print("")


        # Get the keyword information for all the gene clusters of the fetched proteins 
        if cluster_dict:
            # Create a comma-separated list of placeholders for the IN clause
            placeholders = ','.join(['?'] * len(cluster_dict.keys()))
            query = f"SELECT DISTINCT Keywords.clusterID, keyword, completeness, collinearity FROM Keywords WHERE Keywords.clusterID IN ({placeholders})"

            # Execute the query with all clusterIDs at once
            cur.execute(query, tuple(cluster_dict.keys()))

            # Fetch all results at once
            rows = cur.fetchall()

            # Iterate over the rows and process the results
            for index, row in enumerate(rows):
                print(f"\tFetched keywords: {index + 1}     ", end="\r")
                
                clusterID = row[0]
                keyword = row[1]
                completeness = row[2]
                collinearity = row[3]
                
                cluster_dict[clusterID].add_keyword(keyword, completeness,collinearity)

        print("")

        # Fetch the lineage information for all proteins in the protein_dict
        if not taxon_dict:
            taxon_dict = {}

            # Convert genomeID_set to a tuple and create placeholders
            placeholders = ','.join(['?'] * len(genomeID_set))
            query = f"SELECT Genomes.genomeID, Superkingdom, Phylum, Class, Ordnung, Family, Genus, Species FROM Genomes WHERE genomeID IN ({placeholders})"
            
            # Execute the query with all genomeIDs at once
            cur.execute(query, tuple(genomeID_set))
            
            # Fetch all results at once
            rows = cur.fetchall()
            
            # Process the fetched rows
            for index, row in enumerate(rows):
                print(f"\tFetched taxonomy: {index + 1}", end="\r")
                taxon_dict[row[0]] = myUtil.taxonomy_lineage(row, trennzeichen)

            print("")  # Ensure the last print finishes on a new line

    return protein_dict, cluster_dict, taxon_dict
    

def generate_fetch_query(genomes, proteins, keywords, taxon_dict):
    """
    Generate a SQL query for fetching data based on proteins and keywords.

    Args:
        proteins (list): List of proteins.
        keywords (list): List of keywords.
        taxon_dict (dict): Dictionary of taxonomy information.

    Returns:
        tuple: SQL query string and list of arguments for the query.
    """
    query = """
        SELECT DISTINCT 
            Proteins.proteinID, Genomes.genomeID, Proteins.clusterID, Proteins.contig,
            Proteins.start, Proteins.end, Proteins.strand, Proteins.sequence,
            Domains.domain, Domains.domStart, Domains.domEnd, Domains.score, Proteins.dom_count
        FROM Proteins
        LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
        LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
        LEFT JOIN Genomes ON Proteins.genomeID = Genomes.genomeID
    """
    conditions = []
    args = []
    if genomes:
        conditions.append("(" + " OR ".join(["Genomes.genomeID LIKE ?"] * len(genomes)) + ")")
        args.extend([f"%{genome}%" for genome in genomes])
        
    if proteins:
        conditions.append("(" + " OR ".join(["Domains.domain LIKE ?"] * len(proteins)) + ")")
        args.extend([f"%{protein}%" for protein in proteins])

    if keywords:
        conditions.append("(" + " OR ".join(["Keywords.keyword LIKE ?"] * len(keywords)) + ")")
        args.extend([f"%{keyword}%" for keyword in keywords])
    
    if taxon_dict:
        conditions.append("Genomes.genomeID IN (" + ",".join("?" * len(taxon_dict)) + ")") #Only exact matches
        args.extend(taxon_dict.keys())

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    return query, args
    





########## File Output Routines ##########




def output_distinct_fasta_reports(directory,protein_dict,cluster_dict,writemode="w"):
    """
    01.11.22
    
    Args:
        Filepath    Outputfile
        protein_dict    proteinID => proteinObj
        cluster_dict    key:clusterID => value:clusterObj
        Writemode   default "w" for overwrite/new file otherwise "a" for append
    Output:
        Fasta file
        Header format: genomeID + 
                       proteinID get_domains domain_scores domain_coordinates +
                       keyword completeness csb
    21.11.22
        output distinct domains from fusion proteins in extra files marked with fusion-domain_{HMM}.faa
    """
    
#    proteinID_list = sorted(protein_dict, key=lambda x: \
#    (protein_dict[x].gene_contig, protein_dict[x].gene_start)) 
    
    HMM_dict = {}
    fusion_dict = {}
    files = set()
    
    #group proteins in a dict by their domain type  domain => list of protein ids with this type
    for proteinID in protein_dict.keys():
        protein = protein_dict[proteinID]
        HMM = protein.get_domains()
        if not HMM in HMM_dict:
            HMM_dict[HMM] = [protein.proteinID]
        else:
            HMM_dict[HMM].append(protein.proteinID)
        if protein.get_domain_count()>1: #fusion protein detected
            
            fusion_dict[proteinID] = protein
    
    for HMM,proteinID_list in HMM_dict.items():
        filepath = directory + f"_{HMM}.faa"
        files.add(filepath)
        with open(filepath, writemode) as writer:
            
            for proteinID in proteinID_list:
            #Write each line includes all information about one protein
                protein = protein_dict[proteinID]
                genomeID = protein.genomeID
                proteinlist = protein.get_protein_list()  #list representation of a protein object
                sequence = str(protein.protein_sequence)
                sequence.replace('*','')
                
                clusterID = protein.clusterID
                if clusterID in cluster_dict:
                    cluster = cluster_dict[protein.clusterID]
                    clusterlist = cluster.get_cluster_list(",")
                    out = '>' + genomeID + '-' +' '.join(proteinlist[:-5]) + ' ' + ' '.join(clusterlist) + '\n'
                    writer.write(out)
                    writer.write(sequence + '\n')
                else:
                    out = '>' + genomeID + '-' +' '.join(proteinlist[:-5]) + '\n'
                    writer.write(out)
                    writer.write(sequence + '\n')
                


    for protein in fusion_dict.values():
    #Ausgabe der domänen aus fusionsproteinen als einzelne domänen     
        #erkennen welche domänen vorhanden sind
        clusterID = protein.clusterID
        domain_dict = protein.get_domains_dict()
        sequence = str(protein.protein_sequence)
        for domain in domain_dict.values():
            domain.HMM
            domain_sequence = sequence[domain.start:domain.end]         #slice sequence
            filepath = directory + f"_fused_domain_{HMM}.faa"
            writer = open(filepath, "a")         #file finden falls vorhanden sonst neuer file
            
            
            if clusterID in cluster_dict:
                cluster = cluster_dict[protein.clusterID]
                clusterlist = cluster.get_cluster_list(",")
                out = '>' + genomeID + '-' +' '.join(proteinlist[:-5]) + ' ' + ' '.join(clusterlist) + '\n'
                writer.write(out)
                writer.write(domain_sequence + '\n')
            else:
                out = '>' + genomeID + '-' +' '.join(proteinlist[:-5]) + '\n'
                writer.write(out)
                writer.write(domain_sequence + '\n')        

            writer.close()        
        
    
    
    
    return files



        







##############################################################
##########          Fasta File Postprocessing       ##########
##############################################################

def singletons(directory,filepaths):
    """
    01.11.22
        Args:
            directory to store new files
            filepaths with the files to be processed
            
        Routine should filter all singletons, meaning each new file will result in one protein per type per genome
        paralogs are filtered out and should be wrote separately. Output is a faa fasta file
    """
    
    for filepath in filepaths:
        record_dict = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        genome_list = []
        for record in record_dict.keys():
            genomeID = record.split("-")[0]
            genome_list.append(genomeID)

        doublicate = set([x for x in genome_list if genome_list.count(x) > 1])
        
        path = myUtil.removeExtension(filepath)
        name = myUtil.getFileName(path)
        single = open(directory+f"/{name}_singleton.faa","w")
        double = open(directory+f"/{name}_doublicate.faa","w")
        
        for record in record_dict.keys():
            #print(record_dict[record])
            sequence = record_dict[record]
            genomeID = record.split("-")[0]
            if not genomeID in doublicate:
                single.write(">"+sequence.description+"\n")
                single.write(str(sequence.seq)+"\n")
            else:
                double.write(">"+sequence.description+"\n")
                double.write(str(sequence.seq)+"\n")
        single.close()
        double.close()        




def generate_taxonomy_summary(input_file, output_file):
    # Initialize a dictionary to store the results based on the adjusted mapping
    taxonomy_levels = {
        'Kingdom': {},
        'Phylum': {},
        'Class': {},
        'Order': {},
        'Family': {},
        'Genus': {},
        'Species': {}
    }

    # Open the input file and process each line
    with open(input_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            genome_id = columns[0]  # First column is the genome ID
            taxonomy = columns[-1].split(';')  # Last column is the concatenated taxonomy string
            
            # Map the taxonomy levels according to your specification
            if len(taxonomy) > 1:
                if taxonomy[1] not in taxonomy_levels['Kingdom']:
                    taxonomy_levels['Kingdom'][taxonomy[1]] = set()
                taxonomy_levels['Kingdom'][taxonomy[1]].add(genome_id)
            
            if len(taxonomy) > 3:
                if taxonomy[3] not in taxonomy_levels['Phylum']:
                    taxonomy_levels['Phylum'][taxonomy[3]] = set()
                taxonomy_levels['Phylum'][taxonomy[3]].add(genome_id)
            
            if len(taxonomy) > 4:
                if taxonomy[4] not in taxonomy_levels['Class']:
                    taxonomy_levels['Class'][taxonomy[4]] = set()
                taxonomy_levels['Class'][taxonomy[4]].add(genome_id)
            
            if len(taxonomy) > 5:
                if taxonomy[5] not in taxonomy_levels['Order']:
                    taxonomy_levels['Order'][taxonomy[5]] = set()
                taxonomy_levels['Order'][taxonomy[5]].add(genome_id)
            
            if len(taxonomy) > 6:
                if taxonomy[6] not in taxonomy_levels['Family']:
                    taxonomy_levels['Family'][taxonomy[6]] = set()
                taxonomy_levels['Family'][taxonomy[6]].add(genome_id)
            
            if len(taxonomy) > 7:
                if taxonomy[7] not in taxonomy_levels['Genus']:
                    taxonomy_levels['Genus'][taxonomy[7]] = set()
                taxonomy_levels['Genus'][taxonomy[7]].add(genome_id)
            
            if len(taxonomy) > 8:
                if taxonomy[8] not in taxonomy_levels['Species']:
                    taxonomy_levels['Species'][taxonomy[8]] = set()
                taxonomy_levels['Species'][taxonomy[8]].add(genome_id)

    # Write the summary to the output file
    with open(output_file, 'w') as file:
        file.write('Taxonomic Level,Taxon,Count\n')
        for level in taxonomy_levels:
            for taxon, genomes in taxonomy_levels[level].items():
                file.write(f'{level}\t{taxon}\t{len(genomes)}\n')

    print(f"Summary file generated: {output_file}")








##############################################################
##########    Fasta File Postprocessing ON COMMAND    ########
##############################################################


def add_taxonomy(database,filepath,trennzeichen='_'):
    """
    12.03.23
    taxonomy lineage and domain type should be added to a fasta file and replace the
    previous headers
    trennzeichen is the character between the separate informations replacing white characters for
    more readeability
    """
    print("Collecting genome identifier")
    record_dict = {}
    taxon_dict = {}
    writer = open(filepath+".taxon_names","w")    
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
            dataset_range_line = lineage if lineage else record.id
            if record_dict[genomeID] > 1:
                dataset_range_line += '_'+str(record_dict[genomeID])
            

            try:
                dom_type = record.description.split(' ')[1]
                dom_type = ' '+dom_type+' '
            except:
                dom_type = ''
            try:
                #print(f"Insert {genomeID}")
                #print(record_dict[genomeID[:-2]])
                writer.write(">"+dataset_range_line+dom_type+"\n")
                writer.write(str(record.seq) + "\n")
            except:
                print(f"WARNING: No taxonomy found for {genomeID}")

    writer.close()

    
    
    
    
    

def add_genomic_context(database,filepath):
    """
    22.02.23
        Args:
            database for the sequence file 
            filepath with the file to be processed
        This routine takes a fasta file, iterates through the sequences and
        writes down the specific genomic context in which this was found together with
        some genomic features
        
        nicht schön aber funktioniert vielleicht so
    """
    cp_dict = {}
    print("Adding genomic context")
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = "SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({})"

        # Group protein IDs into batches of 1000 for efficient querying
        batch_size = 1000
        protein_ids = [f"'{record.id.split('-', 1)[-1]}'" for record in SeqIO.parse(filepath, "fasta")]

        for i in range(0, len(protein_ids), batch_size):
            batch = ','.join(protein_ids[i:i+batch_size])
            cur.execute(query.format(batch))
            results = cur.fetchall()
            for protein_id, cluster_id in results:
                cp_dict[protein_id] = cluster_id

        
        
        cur.execute("""PRAGMA foreign_keys = ON;""")
        with open(filepath+"_gene_vicinity", "a") as writer:
            writer.write("Index\tProteinID\tDomain(s)\tHit_score\thit_align\tcontig\tstart\tend\tstrand\tSuperkingdom\tClade\tPhylum\tClass\tOrdnung\tFamily\tGenus\tSpecies")
                    
        print(f"\tAssigning genetic environment:")
        for current_proteinID,current_clusterID in cp_dict.items():
            print(f"\tFetched: {current_proteinID} {current_clusterID}",end="\r")
            protein_dict = dict()
            fusion_protIDs = dict()
            query = "SELECT DISTINCT Proteins.proteinID,Proteins.genomeID,Proteins.clusterID,contig,start,end,strand,sequence,domain,domStart,domEnd,score,dom_count from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE Proteins.clusterID = ?"
            cur.execute(query,[current_clusterID,])
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

                if row[12] > 1:
                    fusion_protIDs[row[0]] = 1
                    
            #Also collect domains from fusion proteins
            
            if fusion_protIDs:
                count = 1        
                query = """SELECT DISTINCT proteinID,domain,domStart,domEnd,score FROM Domains WHERE proteinID = ? """
                add_array = []
                proteins_array = []
                args = []
                                    
                for proteinID in fusion_protIDs.keys():
                    arguments = args.copy()
                    arguments.insert(0,proteinID)
                    cur.execute(query,arguments)
                    for row in cur:
                        print(f"\tFetched fused domains: {count}",end="\r")
                        count = count + 1
                        protein = protein_dict[row[0]]
                        protein.add_domain(row[1],row[2],row[3],row[4])
                        
                print("")
            #an diesem punkt haben wir ein vollständiges unsortiertes protein_dict aber keine txonomy information
            # taxonomie info holen
            taxon_dict = dict()
            query = "SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ? "
            genomeID = current_proteinID.split('-')[0]
            cur.execute(query,[genomeID,])
            for row in cur:
                taxon_dict[row[0]] = row[1:]
            
            #create indices for the
            proteinID_list = sorted(protein_dict, key=lambda x: protein_dict[x].gene_start)
            try:
                index = proteinID_list.index(current_proteinID) * -1
            except:
                index = 0    
            with open(filepath+"_gene_vicinity", "a") as writer:
                for proteinID in proteinID_list:
                    #Write each line includes all information about one protein
                    protein = protein_dict[proteinID]
                    proteinlist = protein.get_protein_list()  #list representation of a protein object
            
                    out = str(index) + '\t' +'\t'.join(proteinlist) + '\t'.join(taxon_dict[genomeID]) + '\n'
                    index = index + 1
                    writer.write(out)
            
    
    # Specify input and output file names
    input_file = filepath+"_gene_vicinity"
    output_file = filepath+"_gene_vicinity_sorted"

    # Specify the headers of the columns to sort onSuperkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species 
    sort_headers = ['Superkingdom', 'Phylum', 'Class', 'Ordnung', 'Family', 'Genus', 'Species']

    # Read in the input file and store the rows as a list of dictionaries
    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = [row for row in reader]

    # Sort the rows hierarchically based on the sort columns
    sorted_rows = sorted(rows, key=lambda x: (x[sort_headers[0]], x[sort_headers[1]],x[sort_headers[2]],x[sort_headers[3]],x[sort_headers[4]],x[sort_headers[5]]))

    # Write the sorted rows to the output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow(row)



    
def fetch_database_all(database):
    """
    1.10.22
        Args:
            ONLY FOR DEBUGGING
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" SELECT * FROM Proteins; """)
        con.commit()
        #print("----------Proteins--------------")
        #print(cur.fetchall())
        
        cur = con.cursor()
        cur.execute(""" SELECT * FROM Clusters; """)
        con.commit()
        print("----------Clusters--------------")
        print(cur.fetchall())


        cur = con.cursor()
        cur.execute(""" SELECT * FROM Keywords; """)
        con.commit()
        print("----------Keywords--------------")
        print(cur.fetchall())
        
        cur = con.cursor()
        cur.execute(""" SELECT * FROM Genomes; """)
        con.commit()
        print("----------Genomes--------------")
        print(cur.fetchall())
    return



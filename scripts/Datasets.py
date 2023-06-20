#!/usr/bin/python
#		Module PrepareGenomicData
#		Subroutine FaaToGff
#		Subroutine Translation

import sqlite3
import random

from Bio import SeqIO
from . import ParseReports
from . import myUtil
from operator import add
from itertools import combinations
from math import floor

##############################################################
########## Fetch information from database routines ##########
##############################################################

    #TODO the keys and proteins value has to be controlled by the upper layer in main otherwise there occur errors
def fetch_binary_dataline(database,taxon_dict,proteins,keywords,fusions,protkeys,min_cluster_completeness=0,keys_and_proteins=0):
    """
    04.11.22
        Args:
            database        Name of the database to be worked on
            taxon_dict      dictionary with limiting genomeIDs => taxonomic lineage
            protein list    limits the db fetch to specific protein type
            keyword list    limits the db fetch to specific cluster keyword
            protkeys        empty list
            fusions         empty list
            min_cluster_completeness     minimal completeness for cluster keywords to occur in the output
        Return:
            dataset_dict dictionary with limited genomeIDs => list with binary data
            headers     headline/ name of each list            
            
            In the prediction that taxon dict datasets will only be displayed on somewhat species tree
            dataset_dict may be used for genome statistics
    22.11.22
        extend routine to fetch specific fusion proteins
        and
        proteins combined with keyword; additionally select and count these proteins without this keyword
    12.12.22
        proteins are now searched in keyword clusters if any keywords are provided. If proteins are meant to be
        searched alone please use separate fetch commands for iTol datasets    
    
    Searches for presence/absence matrix. If proteins are given, protein presence. If keyword is given keyword presence.
    If both is given it searches for presence of proteins with this keyword
    22.03.23
        remodel the routine
    """
    dataset_dict = {}
    for key in taxon_dict.keys():
        dataset_dict[key] = []
    
    #define the taxon_dict limiter
    genome_ids = "'" + "', '".join(taxon_dict.keys()) + "'"
    limiter = " AND Proteins.genomeID IN ({genome_ids})"
    limiter = limiter.format(genome_ids=genome_ids)
        
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        query = "SELECT Proteins.genomeID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE domain LIKE ? "
        

        
        keywords_array = []
        keyword_args = []
        keyword_query_string = ""
        for keyword in keywords:
            keywords_array.append(" keyword LIKE ? ")
            keyword_args.append('%'+keyword+'%')

        
        # searching for protein domains in clusters with keyword X
        if keys_and_proteins:
            if keywords_array:
                keyword_query_string = " AND (" + "OR".join(keywords_array) + ")"
                query = query + keyword_query_string    
                print("Searching for entries:" + str(len(proteins)) + " proteins with " + str(len(keywords)) + " keywords")
                # For each protein check presence in gene cluster with keyword
                for protein in proteins: 
                    args = ['%'+protein+'%'] + keyword_args
                    cur.execute(query,args)
                    fetched = {} #genomeIDs with presenct domain in gene cluster

                    print(f"\tDomains 0",end= "\r")
                    for count,row in enumerate(cur):
                        print(f"\tFound {count} domains for protein {protein}",end= "\r")
                        fetched[row[0]] = 1
                        
                    print(f"\n\tSaving presence absence of domains {protein}")    
                    for genomeID,listing in dataset_dict.items():
                        if genomeID in fetched:    
                            listing.append(1)
                        else:
                            listing.append(0)


        #If not keywords and proteins where provided do it separately
        else:
            print("Searching for entries: " + str(len(proteins)) + " proteins and " + str(len(keywords)) + " keywords independently")        
            query = "SELECT Proteins.genomeID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE domain LIKE ? " + limiter

            for protein in proteins:
                cur.execute(query,['%'+protein+'%'])
                fetched = {}
                #print("Found: " + str(len(cur)) + " entries like " + protein)
                print(f"\tDomains 0",end= "\r")
                for count,row in enumerate(cur):
                    print(f"\tFound {count} domains for protein {protein}",end= "\r")
                    fetched[row[0]] = 1
                print(f"\n\tSaving presence absence of domains {protein}")    
                for genomeID,listing in dataset_dict.items():
                    if genomeID in fetched:    
                        listing.append(1)
                    else:
                        listing.append(0)            
                        
            query = "SELECT DISTINCT Proteins.genomeID from Proteins LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE keyword LIKE ?" + limiter
            for keyword in keywords:
            
                cur.execute(query,['%'+keyword+'%'])
                fetched = {}
                for count,row in enumerate(cur):
                    fetched[row[0]] = 1
                    print(f"\tFound {count} clusters for keyword {keyword}",end= "\r")
                print(f"\n\tSaving presence absence of keyword {keyword}")
                for genomeID,listing in dataset_dict.items():
                    if genomeID in fetched:    
                        listing.append(1)
                    else:
                        listing.append(0)
         
        #Dealing with use provided fusion proteins 
        for fusion in fusions:
            #Format: A+B
            #parse fusion
            domains = fusion.split("+")
            protein_domains = len(domains)-1 #number of domains in the protein - 1 because of the > comparison
            query = "SELECT DISTINCT Proteins.proteinID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE dom_count > ? AND domain LIKE ? ;"
            domain_dict_set = {} # key domain => value set of proteinIDs
            for domain in domains:
                proteinID_set = set()
                cur.execute(query,[protein_domains,domain])
                
                for row in cur:
                    proteinID_set.add(row[0]) #first collect all multidomain proteins containing each single domain


                domain_dict_set[domain] = proteinID_set

            #find proteinIDs, which are intersect here
            intersection = set.intersection(domain_dict_set.values()) # third find the intersection between the multidomain proteins with either domain, therefore having both of them
            
            query = "SELECT genomeID from Proteins WHERE proteinID = ?"
            fetched = {}
            for proteinID in intersection:
                cur.execute(query,[proteinID])
                for row in cur:
                    fetched[row[0]] = 1

            for genomeID,listing in dataset_dict.items():
                    if genomeID in fetched:    
                        listing.append(1)
                    else:
                        listing.append(0)
         
         #for proteins with keywords and without keyword
        protkeys_headers = []  
        for protkey in protkeys:
            # will find domain with keyA or keyB or ...
            #format p + k{n}
            protkey_keywords = protkey.split("+")
            domain = protkey_keywords.pop(0)
            args['%'+domain+'%']
            add_array = []
            proteins_array = []
            keywords_array =[]
            
            protkeys_headers.append(protkey)
            protkeys_headers.append(domain + "-".join(protkey_keywords))            
            #parsing the query
            query = "SELECT DISTINCT Proteins.genomeID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE domain LIKE ? "
            for protkey_keyword in protkey_keywords:
                keywords_array.append(" keyword LIKE ? ")
                args.append('%'+keyword+'%')
            if keywords_array:
                add_array.append("(" + "OR".join(keywords_array) + ")")
            query = query + " AND ".join(add_array) + ";"
            
            #executing bulk fetch for protkey
            cur.execute(query,args)
            fetched = {}
            for count,row in enumerate(cur):
                fetched[row[0]] = 1
                print(f"\tKeyword {count}",end= "\r")
            print("\n\tSaving presence absence of domains with keyword")  
            for genomeID,listing in dataset_dict.items():
                if genomeID in fetched:    
                    listing.append(1)
                else:
                    listing.append(0)
                     
            #parsing the query
            query = "SELECT DISTINCT Proteins.genomeID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE domain LIKE ? "
            for protkey_keyword in protkey_keywords:
                keywords_array.append(" NOT keyword LIKE ? ")
                args.append('%'+keyword+'%')
            if keywords_array:
                add_array.append("(" + "OR".join(keywords_array) + ")")
            query = query + " AND ".join(add_array) + ";"
            
            #executing bulk fetch for protkey
            cur.execute(query,args)
            fetched = {}
            for count,row in enumerate(cur):
                fetched[row[0]] = 1
                print(f"\tKeyword {count}",end= "\r")
            print("\n\tSaving presence absence of domains with keyword")  
            for genomeID,listing in dataset_dict.items():
                if genomeID in fetched:    
                    listing.append(1)
                else:
                    listing.append(0) 
            
            
            
            
    headers = [*proteins,*keywords,*fusions,*protkeys_headers]
    return dataset_dict,headers

def iTol_lineage_to_binary_dataline(taxon_dict,dataset_dict):
    #04.11.22
    #concats the lineage to the dataset line by the genomeID
    #key is the genomeID value is the taxonomy + presence absence
    dataline_dict = {}
    for genomeID in taxon_dict.keys():
        dataline = '\t'.join(str(item) for item in dataset_dict[genomeID])
        dataline_dict[genomeID] = taxon_dict[genomeID] + "\t" + dataline + "\n"
        #print(taxon_dict[genomeID])
    
    return dataline_dict    
    
def iTol_binary_file_dataset(directory,dataline_dict,headers):
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
    writer = open(directory+"/iTol_binary_dataset.txt","w")    
    
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
    writer.write("LEGEND_SHAPE_SCALES,1,1,0.5\n")
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
    
def iTol_range_dataset(directory,database,filepath,trennzeichen='_'):
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
            
            
            
            
            
            
            
            
            
            
    
    
    


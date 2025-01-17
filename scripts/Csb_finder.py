#!/usr/bin/python

import re

class Cluster:
    """
    3.9.22
    Cluster organises genes laying in synteny with a specific order. Each cluster has a unique clusterID derived from the assemblyID but with an added index number. 
    
    Args
        clusterID = unique string
        
    11.04.23 added the cluster_start and cluster_end lines    
    """
    def __init__(self,clusterID,distance = 3500):
        self.clusterID = clusterID
        self.genomeID = ""
        self.distance = distance
        self.genes = [] #holds the proteinIDs
        self.types = [] #holds the types, defined as domains of each protein
        self.keywords = []
        self.keywords_dict = dict()
        self.cluster_start = None
        self.cluster_end = None
        
        
    def add_gene(self,proteinID,types,start=None,end=None):
        self.genes.append(proteinID)
        self.types.append(types)
        if self.cluster_start is None:
            self.cluster_start = start
        elif self.cluster_start > start:
            self.cluster_start = start
        if self.cluster_end is None:
            self.cluster_end = end
        elif self.cluster_end < end:
            self.cluster_end = end

        
    def add_keyword(self,keyword,completeness=0,csb="."):
        #self.keywords.append(Keyword(keyword,completeness,csb))
        self.keywords_dict[keyword] = Keyword(keyword,completeness,csb)
    	
    def get_keywords(self):
        listing = []
        for k,v in self.keywords_dict.items():
            listing.append(v)
        
        return listing #returns list of keyword objects

    def get_cluster_start(self):
        return self.cluster_start

    def get_cluster_end(self):
        return self.cluster_end

       
    def get_clusterID(self):
        return self.clusterID
    
    def get_genes(self):
        return self.genes    #genes list
        
    def get_domains(self):
        return self.types    #genes list
        
    def get_domain_set(self):
        tmp = '-'.join(self.types)
        s = set(tmp.split('-'))
        return s    #domain set

    def get_domain_ends_list(self):
        tmp = []
        for typus in self.types:
            parts = typus.split('_')
            if len(parts) > 1:
                # Append the last part after the underscore
                tmp.append(parts[-1])
            else:
                tmp.append(typus)
        return tmp
        
    def get_domain_list(self):
        tmp = '-'.join(self.types)
        s = tmp.split('-')
        return s    #domain list
        
    def get_distance(self):
        return self.distance
        
    def get_keywords_string_obj(self):
        #better formated for database entry
        listing = [self.clusterID]
        for keyword in self.keywords:
            a = keyword.get_keyword()
            b = keyword.get_completeness()
            c = keyword.get_csb()
            listing.append(f"{a} {b} {c}")
        return listing
        
    def get_cluster_list(self,separator):
        """
            Returns a list with keys completeness and csb each in one string
        """
        #not for database but column
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
            keywords_string += separator
            completeness_string += separator
            csb_string += separator
            keywords_string += str(element.get_keyword())
            completeness_string += str(element.get_completeness())
            csb_string += str(element.get_csb())
        #print([keywords_string,completeness_string,csb_string])
        return [self.get_clusterID(),keywords_string] # all keys (compl,csb) in einem string. ,completeness_string,csb_string
	                                       # keys mit i =1 etc.
	
class Keyword:
    """
        3.9.22
        Holds the information about a single keyword its completeness and if it is a csb
    """
    
    def __init__(self,keyword,completeness=0,csb="."):
        self.keyword = str(keyword)
        self.csb = csb
        self.completeness = completeness
        
    def get_keyword(self):
        return self.keyword        
        
    def get_csb(self):
        return self.csb
        
    def get_completeness(self):
        return self.completeness
        
    def set_keyword(self,keyword):
        self.keyword = keyword
        
    def set_csb(self,csb):
        self.csb = csb
        
    def set_completeness(self,completeness):
        self.completeness = completeness
        
def check_order(test_list):
    #3.9.22
    if(all(test_list[i] <= test_list[i + 1] for i in range(len(test_list)-1))):
        return 1
    elif(all(test_list[i] >= test_list[i + 1] for i in range(len(test_list)-1))):
        return 1
    else:
        return 0

def makePatternDict(Filepath):
    #Patterns can have the same name
    pattern = {}        # id => pattern
    pattern_names = {}  # id => name
    i = 1
    with open(Filepath, "r") as reader:
        for line in reader.readlines():
            #print(line)
            #lines = "key1=value1;key2=value2;key3=value3"
            line = line.replace("\n","")
            line = line.replace(" ","")
            if line.strip() == "":
                continue
                
            l = line.split("\t")

            try:
                name = l.pop(0)
                pattern_names[i] = name
                pattern[i] = l
                i = i + 1
            except:
                print(f"\tWARNING: Pattern was not recognized \""+line+"\"")
        #for k, v in pattern.items():
        #   print(k, v)
    return pattern,pattern_names

def find_syntenicblocks(genomeID, protein_dict, distance=3500):
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
    new_sb = 1 #new syntenic block flag
    clusterID_dict = {} #key:clusterID value:clusterObject

    #key=lambda x means anonymous function
    #sort by first contig, then start
    proteinID_list = sorted(protein_dict, key=lambda x: \
    (protein_dict[x].gene_contig, protein_dict[x].gene_start)) 

    clusterID_number = 1
    cluster = Cluster(f"{genomeID}_{clusterID_number}", distance)
    cluster.genomeID = genomeID
    
    for index, elem in enumerate(proteinID_list):
        if index - 1 >= 0:  # Check index bounds
            prev_el_proteinID = str(proteinID_list[index-1])
            curr_el_proteinID = str(elem)
            
            prev_protein = protein_dict[prev_el_proteinID]
            curr_protein = protein_dict[curr_el_proteinID]
            
            if prev_protein.gene_contig == curr_protein.gene_contig and \
               curr_protein.gene_start - prev_protein.gene_end <= distance:
                #print(prev_el_proteinID, curr_el_proteinID)
                if new_sb:
                    # add prev and curr protein to new syntenic block (sb) and sb = false
                    new_sb = 0
                    cluster.add_gene(prev_el_proteinID, prev_protein.get_domains(), prev_protein.gene_start, prev_protein.gene_end)
                    prev_protein.clusterID = cluster.clusterID
                    cluster.add_gene(curr_el_proteinID, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    curr_protein.clusterID = cluster.clusterID
                else:
                    # add current protein to current syntenic block
                    cluster.add_gene(curr_el_proteinID, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                    curr_protein.clusterID = cluster.clusterID
            elif new_sb == 0:
                #print(f"{curr_el_proteinID} was not in range")
                # block was extended but the new element is out of range
                # finalize the current block and add to the list, then reset sb
                clusterID_dict[f"{genomeID}_{clusterID_number}"] = cluster
                clusterID_number += 1
                cluster = Cluster(f"{genomeID}_{clusterID_number}", distance)
                cluster.genomeID = genomeID
                new_sb = 1

    # Add the last cluster to the dictionary if it contains any genes
    if not new_sb:
        clusterID_dict[f"{genomeID}_{clusterID_number}"] = cluster

    #print(clusterID_dict)
    return clusterID_dict



def find_syntenicblocks_bulk_deprecated(protein_dict, distance=3500):
    """
    3.9.22
    Gets a dictionary with protein objects from one genome. Creates a cluster object and adds proteinID belonging to this genetic cluster forming a syntenic block.
    Order of addition to the syntenic block follows the contig and start order. Two different contigs cannot exist inside a syntenic block. Returns a list of cluster objects
    for the cluster analysis.
    
    Args:
        protein_dict - dictionary with key:proteinID and value:proteinObject
        distance - integer of maximal nucleotides distance between genes to consider in synteny
    Return:
        dictionary of cluster objects
    """
    # Sort by genomeID, then contig, and then by start position
    proteinID_list = sorted(protein_dict, key=lambda x: (protein_dict[x].genomeID, protein_dict[x].gene_contig, protein_dict[x].gene_start))
    
    clusterID_number = 1
    new_sb = 1  # new syntenic block flag
    clusterID_dict = {}  # key:clusterID value:clusterObject
    cluster = None
    previous_genomeID = None
    
    for index in range(1, len(proteinID_list)):
        curr_proteinID = proteinID_list[index]
        prev_proteinID = proteinID_list[index - 1]
        
        curr_protein = protein_dict[curr_proteinID]
        prev_protein = protein_dict[prev_proteinID]
        
        genomeID = curr_protein.genomeID  # Retrieve genomeID from the current protein
        clusterID_name = genomeID

        if genomeID != previous_genomeID:
            clusterID_number = 1
            previous_genomeID = genomeID
        
        if cluster is None:
            cluster = Cluster(f"{clusterID_name}_{clusterID_number}", distance)
        
        if prev_protein.gene_contig == curr_protein.gene_contig and curr_protein.gene_start - prev_protein.gene_end <= distance:
            if new_sb:
                # Add previous and current protein to new syntenic block
                new_sb = 0
                cluster.add_gene(prev_proteinID, prev_protein.get_domains())
                prev_protein.clusterID = cluster.clusterID
                cluster.add_gene(curr_proteinID, curr_protein.get_domains())
                curr_protein.clusterID = cluster.clusterID
            else:
                # Add current protein to existing syntenic block
                cluster.add_gene(curr_proteinID, curr_protein.get_domains(), curr_protein.gene_start, curr_protein.gene_end)
                curr_protein.clusterID = cluster.clusterID
        elif new_sb == 0:
            # Close the current syntenic block and start a new one
            new_sb = 1
            cluster.genomeID = genomeID
            clusterID_dict[f"{clusterID_name}_{clusterID_number}"] = cluster
            clusterID_number += 1
            cluster = Cluster(f"{clusterID_name}_{clusterID_number}", distance)
    
    # Handle the last cluster if not added
    if cluster and not new_sb:
        cluster.genomeID = prev_genomeID
        clusterID_dict[f"{clusterID_name}_{clusterID_number}"] = cluster
    return clusterID_dict

def name_syntenicblocks(patterns,pattern_names,clusterID_dict,min_completeness=0.5,collinearity_check=1):
    """
    3.9.22
    Assigns names to syntenic blocks, while ignoring collinearity
    
    Args:
        clusterID_dict - dictionary holding cluster objects
        filepath - path to a textfile containing the named patterns of collinear syntenic blocks
    Returns:
    	clusterID_dict	- dictionary holding the cluster objects 
    	key:clusterID => value:clusterObject
    """
    if not type(patterns) is dict:
        patterns = makePatternDict(patterns)
    for cluster in clusterID_dict.values():
        protein_type_set = cluster.get_domain_ends_list()	#Domänen im cluster geordnet wird nachgeordnet zum set umgewandelt
        #print(protein_type_set)
        for key,pattern in patterns.items():
            keyword = pattern_names[key]
            pattern_set = set(pattern)
            difference = pattern_set.difference(set(protein_type_set))
            completeness = (len(pattern_set)-len(difference))/len(pattern_set)
            if min_completeness <= completeness:
                
                """
                Collinearity check:
                	Alter a copy of the current pattern, so in the next iteration the pattern is
                	still complete. Remove all items from pattern which were not present, check
                	for the rest collinearity by indices. Checks only for complete collinearity
                	of all genes, either on + or - strand. Therefore iterate pattern and get the
                	index of the corresponding pattern element in the syntenic block. Then check
                	for completely ascending or descending index order.
                """
                if collinearity_check:
                    collinearity_pattern = pattern.copy()
                    for item in difference:
                        collinearity_pattern.remove(item)
                    
                    protein_type_list = protein_type_set
                    indices = []
                    for item in collinearity_pattern:
                        index = protein_type_list.index(item)
                        indices.append(index)
                        
                    if check_order(indices):
                        cluster.add_keyword(keyword,completeness,"1")
                        #print(f"Added {keyword}")
                    else:
                        cluster.add_keyword(keyword,completeness,"0")
                        #print(f"Added {keyword}")
                else:
                    cluster.add_keyword(keyword,completeness)
                    #print(f"Added {keyword}")
            #else:
            #    print(f"Not added because {min_completeness} < {completeness} for pattern ")
            #    print(pattern)
            #    print("and gene cluster")
            #    print(protein_type_set)
    return clusterID_dict

   



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
                                


#!/usr/bin/python

import bisect
import os

class InstanceList:
    # Saves the matchpoints and the keys of the gene clusters to be extended in the next round
    
    def __init__(self):
        self.instance_hash = {}
        self.active_keys = set()

    def add_instance(self,key,start, end, extend_tuple = (-1,-1)):
        # Das Tuple (start, end) in die Liste einfügen, während die Sortierung beibehalten wird
        
        if key in self.instance_hash.keys():
            instance_list = self.instance_hash[key]
            if extend_tuple in instance_list:
                instance_list.remove(extend_tuple)
                bisect.insort(instance_list, (start, end))
            else:
                bisect.insort(instance_list, (start, end))
            return 1
        else:
            instance_list = []
            bisect.insort(instance_list, (start, end))
            self.instance_hash[key] = instance_list
            return 1

    def get_instance_hash(self):
        return self.instance_hash

    def get_instance_list(self,key):
        if key in self.instance_hash.keys():
            return self.instance_hash[key]
        else:
            print(f"\tError: Key {key} not found.")
            return {}

    def set_active_keys(self,keys_set):
        #keys of the gene clusters that are currently extending
        self.active_keys = keys_set.copy()
        
    def get_active_keys(self):
        return self.active_keys

    def copy(self):
        new_instance_list = InstanceList()
        new_instance_list.instance_hash = self.instance_hash.copy()
        new_instance_list.active_keys = self.active_keys.copy()
        return new_instance_list


def create_alphabet(hash_table,redundancy_hash,q=2):
    # learns the number of relevant genes from the gene cluster hashtable
    # returns a set of genes that are at least q times in the list
    
    unique_elements = set()
    element_counts = {}
    
    for key,lst in hash_table.items():
        for item in lst:
            element_counts[item] = element_counts.get(item,0) + redundancy_hash[key] # bei redundz minimierten darf das  hier kein feste wert sein, sondern muss der anzahl der reduzierten gene_cluster entsprechen
            
    # Durch die Listen in der Hash-Tabelle iterieren und Elemente zum Set hinzufügen
    for item, count in element_counts.items():
        if count >= q:
            unique_elements.add(item)

    return unique_elements



def create_matchlists(Seqs_hash, Alphabet):
    # für jeden buchstaben des alphabets, eine liste mit den indices verlinkd zu jedem wort
    """
    Input data like this
    gene_clusters = {   "alpha":['rhd','tusA','dsrE','dra'],
                        "beta":['lip','ghd','rhd','tusA','dsrE','dra'],
                        "gamma":['ggf','rhd','dra','tusA','dsrE','dra','pfu'],
                        "delta":['rhd','tusA','dsrE','dra'],
                        "epsilon":['rhd','tusA','dsrE','dra']}
        Datastructure like this
        {
         'rhd': {'alpha': [0], 'beta': [2], 'gamma': [1], 'delta': [0], 'epsilon': [0]},
         'dra': {'alpha': [3], 'beta': [5], 'gamma': [2, 5], 'delta': [3], 'epsilon': [3]},
         'tusA': {'alpha': [1], 'beta': [3], 'gamma': [3], 'delta': [1], 'epsilon': [1]},
         'dsrE': {'alpha': [2], 'beta': [4], 'gamma': [4], 'delta': [2], 'epsilon': [2]}
        }
    """
    
    matchlists = {}  # Ein Dictionary zur Speicherung der Ergebnisse
    
    # Initialisieren Sie die Matchlists für jedes Zeichen in Σ und jede Sequenz in S
    for character in Alphabet:
        matchlists[character] = {}
        for key,lst in Seqs_hash.items():
            matchlists[character][key] = []

    # Durchlaufen Sie die Sequenzen in S
    for key,lst in Seqs_hash.items():
        for i, item in enumerate(lst):
            # Fügen Sie den Index i zur entsprechenden Matchlist für das Zeichen char hinzu
            if item in Alphabet:
                matchlists[item][key].append(i)

    return matchlists


def create_bigrams(Seqs_hash,Alphabet):
    bigrams_set = set()
    for key,lst in Seqs_hash.items():
        for i in range (len(lst) -1):
            sigma_1 = lst[i]
            sigma_2 = lst[i+1]
            if sigma_1 in Alphabet and sigma_2 in Alphabet:
                bigrams_set.add((sigma_1,sigma_2))
    return bigrams_set #returns a set with 2 tuple as elements


def create_next_match(Seqs_hash, matchlists, bigrams):
    #Seqs ist die referenz y wird zum vergleich herangezogen
    """
    Returns datastructure:
    {'alpha': {(1, 'dsrE'): 2, (2, 'dra'): 3, (0, 'dra'): 3, (0, 'tusA'): 1}, 'beta': {(3, 'dsrE'): 4, (4, 'dra'): 5, (2, 'dra'): 5, (2, 'tusA'): 3}, 'gamma': {(3, 'dsrE'): 4, (4, 'dra'): 5, (1, 'dra'): 2, (1, 'tusA'): 3, (2, 'tusA'): 3}}
    >Das äußere Dictionary hat Schlüssel ('alpha', 'beta', 'gamma') und jeweils einen Wert.
    >Jeder Wert ist ein weiteres Dictionary, das Schlüssel-Wert-Paare enthält.
    >In den inneren Dictionaries sind die Schlüssel Tupel, die aus einem Integer und einem String bestehen, z. B. (1, 'dsrE').
    >Die Werte in den inneren Dictionaries sind Integer, z. B. 2, 3, usw.
    """
    
    next_match = {}  # Ein Dictionary zur Speicherung der Ergebnisse
    
    for key,lst in Seqs_hash.items():
    # Durchlaufen Sie alle zweistelligen Substrings σ1σ2 in S
        next_matchy = {}
        for sigmas in bigrams:
            sigma1 = sigmas[0]
            sigma2 = sigmas[1]
            
            # Überprüfen, ob σ1σ2 ein Substring einer Eingabesequenz ist
            if sigma1 in matchlists and key in matchlists[sigma1]: 
                # string_1 zeichen in der matchlist, string_2 hat enthält
                # das string_1 zeichen. Dann hole mir die matchpoints für 
                # string2 für das string_1 zeichen
                  
                matchlist_sigma1_y = matchlists[sigma1][key]

                # Durchlaufen Sie die Indizes in M atchlists[σ1][y]
                for idx1 in matchlist_sigma1_y:
                    
                    # Durchlaufen Sie M atchlists[σ2][y] bis zum ersten Index j, für den j > idx1 gilt
                    matchlist_sigma2_y = matchlists.get(sigma2, {}).get(key, [])
                    for idx2 in matchlist_sigma2_y:
                        if idx2 > idx1:
                            # Speichern Sie j als NextMatchy(i, σ2)
                            next_matchy[(idx1, sigma2)] = idx2
                            break  # Wir haben den ersten passenden Index gefunden, brechen Sie die Schleife ab
        next_match[key] = next_matchy


    #Frage: Ist es notwendig einen nächsten punkt zu speichern, der weiter weg ist als k (number of insertions) erlaubt?
    return next_match



def initialize_Instances(Pattern, current_Seq_key, Seq_hash, InstanceListP, MatchLists, NextMatch):
    #Initilizes a new InstanceList of the Pattern length is 1. Corresponds to the matchpoints of the gene
    extended_keys = set()
    for key,lst in Seq_hash.items():
        if key != current_Seq_key:
            add_flag = 0
            sigma = Pattern[-1]  # Letztes Zeichen von P
            MatchList = MatchLists[sigma][key]
            #add_flag = initInstances(MatchList,InstanceListP,key)
            for p in MatchList:
                # Hinzufügen einer minimalen k-Instanz (p, p) zu InstanceListP_y
                add_flag = InstanceListP.add_instance(key, p, p)

            if add_flag:
                extended_keys.add(key)

    return extended_keys




def findInstances(Pattern, current_Seq_key, Seq_hash, InstanceListP, MatchLists, NextMatch, expanding_set,k,q):
    extended_keys = set()
    for key,lst in Seq_hash.items():
        if key != current_Seq_key and key in expanding_set:
            add_flag = 0
            sigma = Pattern[-1]  # Letztes Zeichen von P
            MatchList = MatchLists[sigma][key]
            add_flag = extendInstances(Pattern, key, InstanceListP, MatchList, NextMatch,k,q) #extends a csb in the list returns 1 if successfull
            if add_flag:
                extended_keys.add(key)
    return extended_keys

def extendInstances(Pattern, key, InstanceListP, MatchList, NextMatch,k,q):
    add_flag = 0
    instanceListP_y = InstanceListP.get_instance_list(key).copy()
    for index,chain_tuple in enumerate(instanceListP_y):
        s1,e1 = chain_tuple
        if index + 1 < len(instanceListP_y): # find the next csb in the chain or set infinite
            s2,e2 = instanceListP_y[index+1]
        else:
            e2 = float('inf')
            
        next_character = Pattern[-1]
        try:
            next_idx_j = NextMatch[key][(e1,next_character)]

            if next_idx_j <= e2 and ((next_idx_j-s1)-len(Pattern)+1) <= k:
                add_flag = InstanceListP.add_instance(key,s1,next_idx_j,chain_tuple)
             
        except KeyError:
            
            continue
        #'alpha': {(1, 'dsrE'): 2, (2, 'dra'): 3, (0, 'dra'): 3, (0, 'tusA'): 1}

    return add_flag


def expand_redundancy(redundancy_hash,key,gene_cluster_identifier_set):
    count = redundancy_hash[key]
    for identifier in gene_cluster_identifier_set:
        count = count + redundancy_hash[identifier]
    return count






#redundancy_hash = {"alpha":1,"beta":1,"gamma":1,"delta":1,"epsilon":1,"phi":11}

#gene_clusters = {"alpha":['rhd','tusA','dsrE','dra'],
#"beta":['lip','ghd','rhd','tusA','dsrE','dra','fix','chi','snu'],
#"gamma":['ggf','rhd','dra','tusA','dsrE','dra','pfu','fix','chi','snu'],
#"delta":['rhd','tusA','dsrE','dra','ttr','chi','snu'],
#"epsilon":['rhd','tusA','dsrE','dra','ttr'],
#"phi":['yhhp','yhhe','yhhq','yhhr']}
#gene_clusters = {"alpha":['rhd','tusA','dsrE','dra'],"beta":['lip','ghd','rhd','tusA','dsrE','dra']}


def csb_finderS_matchpoint_algorithm(redundancy_hash,gene_clusters,k=1,q=1):
    """
        number of insertions k
        number of occurence q
        csb is defined as occuring pattern of neigbhouring genes in at least q gene clusters with maximal k insertions with minimal length
        
        create_alphabet: get all genes that are minimally occuring q times in sum of all gene clusters
        matchlist: calculate the matchpoints between all gene clusters
        bigrams: creates tuples of neigbouring genes that have an instance in the alphabet 
        Nextmatch: starting with bigrams, where is the next matchpoint. saves this in a complex datastructure
        
        then iterate through all gene clusters, and start with each gene as starting point of a new csb (i:j) mit i == j
        then increase j and search for the extension of the gene cluster
        only the last character of the csb is needed. it will be searched for its next occurence in all other gene clusters. if it is in range of k than added
        else the gene cluster is saved as it is.
        
        if a csb occures twice in then the precomputed Instance List is taken. This List saves the status of all csb
        
        return is a hash with tuples of csb as keys and the gene cluster identifier as values
        
    """

    alphabet = create_alphabet(gene_clusters,redundancy_hash,q) #alphabet is a set

    MatchLists = create_matchlists(gene_clusters, alphabet)

    bigrams = create_bigrams(gene_clusters,alphabet)

    NextMatch = create_next_match(gene_clusters,MatchLists,bigrams)

    computed_Instances_dict = dict()

    for key,lst in gene_clusters.items():
        length = len(lst)
        
        for i in range(0, length):
            active_expanding_keys = set()
            InstanceListP = InstanceList()
            pattern_tuple = 0
            for j in range(i+1,length+1):
                
                 
                pattern = lst[i:j]
                if not pattern[-1] in alphabet:
                    break
                
                pattern_tuple = tuple(pattern)
                if pattern_tuple in computed_Instances_dict.keys():
                    InstanceListP = computed_Instances_dict[pattern_tuple].copy()
                    active_expanding_keys = InstanceListP.get_active_keys()
                       
                elif len(pattern) == 1:
                    
                    initilized = initialize_Instances(pattern, key, gene_clusters, InstanceListP, MatchLists, NextMatch)
                    count = expand_redundancy(redundancy_hash,key,initilized)
                    if count > q:
                        active_expanding_keys = initilized
                else:
                    extended = findInstances(pattern, key, gene_clusters, InstanceListP, MatchLists, NextMatch,active_expanding_keys,k,q)
                    count = expand_redundancy(redundancy_hash,key,extended)
                    if count > q:
                        discontinued = active_expanding_keys - extended # returns the keys of gene clusters that were not extended
                        active_expanding_keys = active_expanding_keys & extended # returns the keys of the gene clusters that were extended and will be tested in the next round
                        
                        
                        
                    else:
                        break
                        
                #Add pattern giving csb itself to the InstanceList        
                InstanceListP.add_instance(key,i,j,(i,j))
                active_expanding_keys.add(key)
                InstanceListP.set_active_keys(active_expanding_keys)
                #save_computed_Instances pattern => current instance list
                computed_Instances_dict[pattern_tuple] = InstanceListP.copy() 
            
                
            #print("\n\nInstance Listing for pattern \t",InstanceListP.get_instance_hash())


    dictionary = {}
    for key,it in computed_Instances_dict.items(): #for the return, make the dict csb_pattern => gene cluster IDs
        #print(key," ",it.active_keys)
        dictionary[key] = it.active_keys
    return dictionary



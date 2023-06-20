#!/usr/bin/python
#		Module PrepareGenomicData
#		Subroutine FaaToGff
#		Subroutine Translation
import sqlite3
import csv
import re

from Bio import SeqIO
from itertools import combinations

from . import myUtil
from . import ParseReports
from . import Csb_finder

"""
Two dictionaries one with protein objects one with cluster objects
Bei proteins sollte die übergeordnete Einheit Genome lauten, bei cluster sollte die übergeordnete einheit keyword lauten
"""

##############################################################
########## Fetch information from database routines ##########
##############################################################

    
    
def fetch_bulk_data(database,proteins,keywords,taxon_dict=1,min_cluster_completeness=0):
    """
    15.10.22
        Args:
           database     Name of the database to be worked on
           genomeID     List of genomeIDs to be retrieved 
           taxon        limits to specific taxon
           protein list      limits to specific protein type
           keyword list      limits to specific cluster keyword
           min_cluster_completeness     minimal completeness for cluster keywords to occur in the output
        Return:
            protein dictionary with key:proteinID => value:protein object for a single genome
            
            
    01.11.22
        Routine accepts now lists for parameters protein and keyword
        taxon information dict 
    09.11.22
        This routine does not work properly for big datasets because JOIN is far too slow
        therefore use fetch protein with keyword or fetch keyword with protein subroutine
        #SOLVED by indexing foreign keys
    10.11.22
        Streamlining this routine to higher performance out of the experiments with the
        JOIN workaround
    21.11.22
        Addition: For two domain proteins the unselected domain shall be selected too, in order
        to discriminate between the fusion versions and single domain proteins. maybe an option
        in output would be good to output these domains without any fusion domain
    22.03.23
        Included Limit to genomeID in taxon_dict during the fetch
    """
    query = "SELECT DISTINCT Proteins.proteinID,Proteins.genomeID,Proteins.clusterID,contig,start,end,strand,sequence,domain,domStart,domEnd,score,dom_count from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID"
    add_array = []
    proteins_array = []
    keywords_array =[]
    args = []
    
    if proteins:
        for protein in proteins:
            proteins_array.append(" domain LIKE ? ")
            args.append('%'+protein+'%')
        if proteins_array:
            add_array.append("(" + "OR".join(proteins_array) + ")")    
    
    if keywords:
        query += " LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID"
        for keyword in keywords:
            keywords_array.append(" keyword LIKE ? ")
            args.append('%'+keyword+'%')
        if keywords_array:
            add_array.append("(" + "OR".join(keywords_array) + ")")
    
    # Limit to genomeIDs in taxon_dict
    genome_ids = "'" + "', '".join(taxon_dict.keys()) + "'"
    limiter = " AND Proteins.genomeID IN ({genome_ids})"
    limiter = limiter.format(genome_ids=genome_ids)

    
    if len(args) > 0:
        query += " WHERE "
    query = query + " AND ".join(add_array) + " " + limiter+";"
    
    print(f"\tFetching data matching description: {args}\n")
    protein_dict = {}
    cluster_dict = {}
    fusion_protIDs = {}
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(query,args)
        #Fill the protein dict
        count = 1
        for row in cur:
            print(f"\tFetched domains: {count}",end="\r")
            count = count +1 
            # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => contig,
            # 4 => start, 5 => end, 6 => strand, 7 => sequence,
            # 8 => domain, 9 => domStart, 10 => domEnd, 11 => score,
            if row[1] in taxon_dict.keys():
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
            
                if row[12] > 1:
                    fusion_protIDs[row[0]] = 1
                
        #Also collect domains from fusion proteins
        print("")
        if fusion_protIDs:
            count = 1        
            query = """SELECT DISTINCT proteinID,domain,domStart,domEnd,score FROM Domains WHERE proteinID = ? """
            add_array = []
            proteins_array = []
            args = []
            #for protein in proteins:
            #    proteins_array.append(" NOT domain LIKE ? ")
            #    args.append('%'+protein+'%')
            #query += "(" + "AND".join(proteins_array) + ")"
                
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
        count = 1
        query = """SELECT DISTINCT Keywords.clusterID,keyword,completeness,collinearity FROM Keywords WHERE Keywords.clusterID = ?"""
        for clusterID,clusterObj in cluster_dict.items():
            print(f"\tFetched keywords: {count}",end="\r")
            count = count +1
            #    0 => clusterID, 1 => keyword, 2 => completeness, 3 => collinearity, 4 => proteinID
            
            cur.execute(query,[clusterID])
                
            for row in cur:
                if min_cluster_completeness < float(row[2]):
                    clusterObj.add_keyword(row[1],row[2],row[3])              
                
    con.commit()
    con.close()    
    return protein_dict,cluster_dict


def fetch_limiter_data(database,lineage,taxon,proteins,keywords,min_cluster_completeness=0.75,trennzeichen="_"):
    """
    03.11.22
        Args:
           database     Name of the database to be worked on
           taxon        limits to specific taxon
           protein list      limits to specific protein type
           keyword list      limits to specific cluster keyword
           min_cluster_completeness     minimal completeness for cluster keywords to occur in the output
        Return:
            taxon dictionary genomeID => taxonomy lineage
    22.02.23
        copy of the fetch_taxonomy_data() routine but this time with speed up for processes
    11.03.23
        Moved from Datasets to here, because this routine fetches all genomes the search is limited to
        therefore it is rather an output routine than dataset generation routine      
         
    """

    query = "SELECT Genomes.genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species from Genomes"
    if lineage and taxon:
        query += " WHERE "
        query += f" {lineage[0]} LIKE ? "

    taxon_dict = {}
    print(query)
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")

        if taxon:
            cur.execute(query,taxon)
        else:
            cur.execute(query)
        print("\tSelecting taxonomy information")
        for index,row in enumerate(cur):
            print(f"\tSelecting genomes {index+1}",end = "\r")
            if not row[0] in taxon_dict:
                taxon_dict[row[0]] = myUtil.taxonomy_lineage(row,trennzeichen)

    
    # Execute block if limiting domains are given
    # selects relevant genome IDs and removes all others from the taxon dictionary
    if proteins:
        query = "SELECT DISTINCT Proteins.genomeID from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE "
        args = []
        proteins_array = []
        limiter_dict = dict()
        for protein in proteins:
            proteins_array.append(" domain LIKE ? ")
            args.append('%'+protein+'%')
        query += ("(" + "OR".join(proteins_array) + ")")    
        with sqlite3.connect(database) as con:
            cur=con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute(query,args)
            for index,row in enumerate(cur):
                print(f"\tSelecting genomes with given domain(s) {index+1}",end = "\r")
                limiter_dict[row[0]] = 1
        #remove the taxon_dict keys which are also in the limiter dict
        for key in limiter_dict.keys():
            if key in taxon_dict:
                del taxon_dict[key]
    

    # Execute block if limiting keywords are given
    # selects relevant genome IDs and removes all others from the taxon dictionary
    if keywords:
        query = "SELECT DISTINCT Proteins.genomeID from Proteins LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID WHERE "
        args = []
        add_array = []
        keywords_array =[]
        limiter_dict = dict()
        for keyword in keywords:
            keywords_array.append(" keyword LIKE ? ")
            args.append('%'+keyword+'%')
        query += ("(" + "OR".join(keywords_array) + ")")
        with sqlite3.connect(database) as con:
            cur=con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute(query,args)
            for index,row in enumerate(cur):
                print(f"\tSelecting genomes with given keyword(s) {index+1}",end = "\r")
                limiter_dict[row[0]] = 1
        #remove the taxon_dict keys which are also in the limiter dict
        for key in limiter_dict.keys():
            if key in taxon_dict:
                del taxon_dict[key]
                       
    con.commit()
    con.close()
    
    return taxon_dict
    


def fetch_HMM_groupe_proteinIDs(database,HMM):
    #01.11.22
    proteinIDs = set()
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT DISTINCT proteinID from Proteins JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE domain rlike = ? """,[HMM])
        
        for row in cur:
            proteinIDs.add(row[0])
    return proteinIDs



########## File Output Routines ##########
def output_genome_report(filepath,protein_dict,cluster_dict,taxon_dict=dict(),writemode="w"):
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
    """
    with open(filepath, writemode) as writer:
    
    
        proteinID_list = sorted(protein_dict, key=lambda x: \
        (protein_dict[x].genomeID, protein_dict[x].gene_contig, protein_dict[x].gene_start)) 
        
        for proteinID in proteinID_list:
        #Write each line includes all information about one protein
            protein = protein_dict[proteinID]
            proteinlist = protein.get_protein_list()  #list representation of a protein object
            
            out = protein.genomeID + '\t' +'\t'.join(proteinlist)
            
            clusterID = protein.get_clusterID()
            if clusterID:
                
                if clusterID in cluster_dict:
                    cluster = cluster_dict[clusterID]
                    clusterlist = cluster.get_cluster_list(",")
                    out += '\t' + '\t'.join(clusterlist)
                
                else:
                    out += '\t'+ clusterID + '\t\t\t'
            else:
                out += '\t\t\t\t'
            
            genomeID = protein.get_genomeID()
            if genomeID in taxon_dict:
                taxon = taxon_dict[genomeID]
                out += '\t' + taxon
            else:
                out += '\t\t\t\t\t\t\t\t\t\t\t'
            out += '\n'                
            writer.write(out)
    return



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
    
    for proteinID in protein_dict.keys():
        protein = protein_dict[proteinID]
        HMM = protein.get_domains()
        if not HMM in HMM_dict:
            HMM_dict[HMM] = [protein.proteinID]
        else:
            HMM_dict[HMM].append(protein.proteinID)
        if protein.get_domain_count()>1:
            #fusion protein detected
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
                sequence = str(protein.get_protein_sequence())
                sequence.replace('*','')
                
                clusterID = protein.get_clusterID()
                if clusterID in cluster_dict:
                    cluster = cluster_dict[protein.get_clusterID()]
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
        clusterID = protein.get_clusterID()
        domain_dict = protein.get_domains_dict()
        sequence = str(protein.get_protein_sequence())
        for domain in domain_dict.values():
            domain.HMM
            domain_sequence = sequence[domain.start:domain.end]         #slice sequence
            filepath = directory + f"_fused_domain_{HMM}.faa"
            writer = open(filepath, "a")         #file finden falls vorhanden sonst neuer file
            
            
            if clusterID in cluster_dict:
                cluster = cluster_dict[protein.get_clusterID()]
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

def output_combined_keyword_fasta_reports(directory,filepaths,keywords,writemode="w"):
    """
    07.12.22
    Args:
        Directory   Write new files to this directory path
        Filepath    Input fasta files
        keywords    Searched keywords for the combination
        Writemode   default "w" for overwrite/new file otherwise "a" for append
    Output:
        Fasta file
    
    Combine the cluster keywords and create file output for each keyword combination
    Output includes the sequences which are assigned not only to one keyword but to the combined keywords 
    """



    #get all possible combination of keywords
    all_keyword_combinations = []
    for keyword in keywords:
        all_keyword_combinations.append((f"{keyword}",)) #append by a tuple holding a keyword

    if len(keywords) > 1:    
        print("combine keywords")
        for i in range(1,len(keywords)):
	        combined_keys = list(combinations(keywords,i+1)) #combinations creates unique permutations by 2 to all elements
	        all_keyword_combinations.extend(combined_keys)
    print(all_keyword_combinations)    
    
    
    file_list = list()
    for filepath in filepaths:
        path = myUtil.removeExtension(filepath)
        name = myUtil.getFileName(path)
        

        for combine in all_keyword_combinations:

            combined = open(directory+f"/{name}_keywords_{combine}.faa","w")
            file_list.append(directory+f"/{name}_keywords_{combine}.faa")
            for record in SeqIO.parse(filepath, "fasta"):
                description_list = record.description.split(" ")
                if len(description_list)>7:
                    #falls alle aus combine präsent sind aufschreiben
                    keyword_list = description_list[-3].split(",")
                    flag = True
                    for element in combine:
                        element_regex = re.compile(element)
                        if element in keyword_list:
                            flag = True
                        elif any((match := element_regex.match(item)) for item in keyword_list):
                            flag = True
                        else:
                            flag = False 
                            break
                    if flag:
                        SeqIO.write(record,combined,"fasta")
                    
            combined.close()

    return file_list


        







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

##############################################################
##########    Fasta File Postprocessing ON COMMAND    ########
##############################################################

def merge_fasta(directory,filepaths):
    #02.11.22
    #merge to fasta files without doublicates
    outerset = set()
    for filepath in filepaths:
        innerset = set()
        for record in SeqIO.parse(filepath, "fasta"):
            innerset.add(record.id)
        outerset = outerset.union(innerset)

    concat = open(directory,"w")
    for filepath in filepaths:
        for record in SeqIO.parse(filepath, "fasta"):
            if record.id in outerset:
                SeqIO.write(record,concat,"fasta")
                outerset.remove(record.id)
    
    concat.close()
    return directory

def filter_length_fasta(filepath,minimal,maximal):
    #22.11.22
    output = filepath+"_filter_"+str(minimal)+str(maximal)
    filter_writer = open(output,"w")

    for record in SeqIO.parse(filepath, "fasta"):
        sequence_length = len(str(record.seq))
        if sequence_length > minimal and sequence_length < maximal:
            SeqIO.write(record,filter_writer,"fasta")
    
    
    filter_writer.close()    

    return output
    
def concat_alignments(directory,filepaths,add_gaps=0):
    """
    02.11.22
        Args:
            filepaths with the files to be processed
            directory for the concatenated seqs file
        Alignment files should be concated if header is the same            
    """
    concat_dict = {}
    identifier_dict = {}
    alignment_lengths = {}
    #/home/tomohisa/BioprojectMagenta/Program/Compile_v1_1_0/dist/results/Archaea_gtdb/2022-11-02 10:00:30.2846870_0_['LipS1', 'LipS2', 'LplA']_[]_LipS1.faa
    #parse fasta file
    
    #get all genomeIDs and seq lengths
   
    for filepath in filepaths:
        for record in SeqIO.parse(filepath, "fasta"):
            genomeID = record.id.split("-")[0]
            concat_dict[genomeID] = ""
            
            try:
                dom_type = record.description.split(' ')[1]
                dom_type = dom_type.split('_')[-1]
            except:
                dom_type = ''

            if genomeID in identifier_dict.keys():
                description = identifier_dict[genomeID]
                identifier_dict[genomeID] = description + '_' + dom_type
            else:
                identifier_dict[genomeID] = record.id + '_' + dom_type
                
                
            alignment_lengths[filepath] = len(str(record.seq))
            #print(len(str(record.seq)))

    for filepath in filepaths:
        present_set = set()
        for record in SeqIO.parse(filepath, "fasta"):
            genomeID = record.id.split("-")[0]
            if genomeID in concat_dict:
                concat_dict[genomeID] += record.seq
                present_set.add(genomeID)
            
        #absent_genomes = present_set.difference(set(concat_dict.keys()))
        absent_genomes = set(concat_dict.keys()).difference(present_set)
        length = alignment_lengths[filepath]
        
        for absent in absent_genomes:
            if add_gaps:
                print(f"WARNING: Missing sequence for {absent} in {filepath}.\nConcating gaps instead")
                concat_dict[absent] += length*"-"
            else:
                print(f"WARNING: Missing sequence for {absent} in {filepath}.\nRemoving sequence")
                concat_dict.pop(absent)
    concat = open(directory,"w")
    
    for k,v in concat_dict.items():
        h = identifier_dict[k]
        concat.write(">"+h+"\n")
        concat.write(str(v)+"\n") 
        
    concat.close()

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
            dataset_range_line = lineage if record_dict[genomeID] == 1 and lineage else record.id


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



def help_routine_select_random_marker_set(database,filepath):

    print("Collecting genome identifier")
    level_dict = {}
    genomeID_dict = {}
    
        
    print(f"Fetching taxonomic data from {database}")
    #For performance reasons this is written here without a specific subroutine but execute shall be the same select as in Datasets
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        #cur.execute("""SELECT DISTINCT Phylum FROM Genomes WHERE NOT Phylum LIKE "-" AND NOT Phylum LIKE "_" OR NOT Phylum LIKE "1" OR NOT Phylum LIKE "0";""")
        phyla = ["Abyssubacteria","Acidobacteriota","Actinobacteriota","Aerophobota","Aquificota","Armatimonadota","Atribacterota","Aureabacteria","Bacteroidota","Bdellovibrionota","Bipolaricaulota","Caldisericota","Calditrichota","Calescibacterota","Campylobacterota","Chlamydiota","Chloroflexota","Chrysiogenota","Cloacimonadota","Coprothermobacterota","Cyanobacteria","Deferribacterota","Deinococcota","Delongbacteria","Dependentiae","Desantisbacteria","Desulfobacterota","Dictyoglomota","Dormibacterota","Edwardsbacteria","Eisenbacteria","Elusimicrobiota","Eremiobacterota","Fermentibacterota","Fibrobacterota","Firestonebacteria","Firmicutes","Fusobacteriota","Gemmatimonadota","Goldbacteria","Hydrogenedentota","Krumholzibacteriota","Latescibacterota","Lindowbacteria","Margulisbacteria","Marinisomatota","Mcinerneyibacteriota","Methylomirabilota","Moduliflexota","Muirbacteria","Myxococcota","Nitrospinota","Nitrospirota","Omnitrophota","Patescibacteria","Planctomycetota","Poribacteria","Proteobacteria","Ratteibacteria","Riflebacteria","Schekmanbacteria","Spirochaetota","Sumerlaeota","Synergistota","Tectomicrobia","Thermodesulfobiota","Thermosulfidibacterota","Thermotogota","Verrucomicrobiota","WOR-3","Wallbacteria","Zixibacteria"]
        
        for row in phyla:
            level_dict[row] = 1
                
        
        for level in level_dict.keys():
            cur.execute("""SELECT genomeID from Genomes WHERE Phylum = ? LIMIT 5""",[level])        
            for row in cur:
                genomeID_dict[row[0]] = 1
                        
    writer = open(filepath+".limited_phyla","w")
    for record in SeqIO.parse(filepath, "fasta"):
        genomeID = record.id.split("_Bacteria")[0]
        #genomeID = myUtil.getGenomeID(genomeID)
        if genomeID in genomeID_dict.keys():
            writer.write(">"+genomeID+"\n")
            writer.write(str(record.seq) + "\n")
        
    
    
    
    writer.close()
    
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




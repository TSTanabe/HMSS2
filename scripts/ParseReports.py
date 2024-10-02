#!/usr/bin/python
from Bio import SearchIO
from Bio import SeqIO
import traceback
import subprocess
from . import Database
#import multiprocessing
#from multiprocessing import Process, Manager, Pool, Semaphore
import re


class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have a higher score than all overlapped domains. If new domain has a lower score than any previously added domain new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein. Here the additional information of other domains is added if it does not interfere with this assumption
    
    Args:
        protein_ID - unique string as identifier
        HMM - protein type designation
        start - start coordinate of the protein type in the AA sequence
        end - end coordinate of the protein type in the AA sequence
        score - bitscore, propability of the dignated protein type
    """


    def __init__(self,proteinID,HMM,start=0,end=0,score=1,genomeID=""): 
        #Protein attributes
        self.proteinID = proteinID
        self.genomeID = genomeID
        self.protein_sequence = ""
        
        #Gene attributes
        self.gene_contig = ""
        self.gene_start = 0
        self.gene_end = 0
        self.gene_strand = "."
        self.gene_locustag = ""
        
        #Cluster attributes
        self.clusterID = "" #reziprok mit Cluster class
        self.keywords = {} #reziprok mit Cluster class
        
        self.domains = {}   # dictionary start coordinate => Domain object
        self.add_domain(HMM,start,end,score)
        

    ##### Getter ####
            
    def get_domains(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_HMM())
        return '-'.join(listing)
    
    def get_domains_dict(self):
    #return dict
        return self.domains
            
    def get_domain_listing(self):
    #return list
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key])
        return listing
                        
    def get_domain_coordinates(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_start()}:{self.domains[key].get_end()}")
        return '-'.join(listing)

    def get_domain_scores(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_score()}")
        return '-'.join(listing)   

    def get_domain_count(self):
        return len(self.domains)
        
    def get_protein_string(self):
        #3.9.22 representation of the whole protein in one line
        a = self.proteinID
        b = self.get_domains()
        c = self.get_domain_scores()
        d = self.gene_contig
        e = self.gene_start
        f = self.gene_end
        #d = self.get_protein_sequence()
        string = f"{a} {b} {c} {d} {e} {f}"
        return string
        
    def get_protein_list(self):
        #3.9.22 representation of the whole protein in one line
        listing = []
        listing.append(self.proteinID)
        listing.append(self.get_domains())
        listing.append(str (self.get_domain_scores()))
        listing.append(str (self.get_domain_coordinates()))
        listing.append(self.gene_contig)
        listing.append(str (self.gene_start))
        listing.append(str (self.gene_end))
        listing.append(self.gene_strand)
        listing.append(self.gene_locustag)
        #d = self.get_protein_sequence()
        #string = f"{a} {b} {c} {d} {e} {f}"
        return listing
            
    def get_sequence(self):
        return str(self.protein_sequence)
        
    def get_length(self):
        return len(self.protein_sequence)
    ##### Setter #####

    def check_domain_overlap(self,new_start, new_end,\
    current_start,current_end):
    #2.9.22

        if (current_start <= new_start and new_start <= current_end)\
        or (current_start <= new_end and new_end <= current_end):
        #start oder endpunkt innerhalb er grenzen
            return 1

        elif (new_start <= current_start and current_end <= new_end)\
        or (current_start <= new_start and new_end <= current_end):
        #start und end innerhalb der grenzen oder alte domäne innerhalb der neuen
            return 1
            
        elif (new_start <= current_start and new_end <= current_start)\
        or (new_start >= current_end and new_end >= current_end):
        #start und end kleiner als self start oder start und end größer als self end dann
        #domäne außerhalb der alten domäne und adden egal welcher score
            return 0
        return None
        

    def add_domain(self,HMM,start,end,score):
        """
        2.9.22
        Adds a domain to an existing protein. Only if there is no overlap with a 
        currently existing domain or of the overlapping domain scores higher than
        any existing overlapped domain. Overlapped domains with minor scores are deleted
        
        Args:
            HMM - protein type designation
            start - start coordinate of the protein type in the AA sequence
            end - end coordinate of the protein type in the AA sequence
            score - bitscore, propability of the dignated protein type
        Return: 
            0 - no domain was added
            1 - domain was added
        """

        del_domains = [] # start coordinates/keys of domains to be replace
        for domain in self.domains.values():
            if self.check_domain_overlap(start,end,domain.get_start(),domain.get_end()):
                if domain.get_score() < score:
                    del_domains.append(domain.get_start())
                else:
                    return 0


        
        for key in del_domains:
            self.domains.pop(key)
        self.domains.update({start:Domain(HMM,start,end,score)}) # if loop complete
        
        return 1
        
        
        


class Domain:
#2.9.22
    def __init__(self,HMM,start,end,score):
        self.HMM = HMM
        self.start = int(start)
        self.end = int(end)
        self.score = int(score)
    
    def __hash__(self):
        return hash((self.HMM, self.start, self.end, self.score))
    
    def __eq__(self, other):
        if isinstance(other, Domain):
            return self.HMM == other.HMM and self.start == other.start and self.end == other.end and self.score == other.score
        return False
    
    def get_HMM(self):
        return self.HMM
    def get_start(self):
        return self.start
    def get_end(self):
        return self.end
    def get_score(self):
        return self.score
    
#   Parsing subroutines
#------------------------------------------------------------


def parseHMMreport(Filepath, Thresholds, cut_score=10):
    protein_dict = {}
    
    try:
        protein_dict = parseHMMreport_hmmer3_format(Filepath, Thresholds, cut_score)
    except Exception as e:  # Catch all exceptions
        print(f"Error occurred while processing the file: {Filepath}")
        traceback.print_exc()  # This prints the complete traceback of the error
        

    return protein_dict
    
    
def parseHMMreport_hmmer3_format(Filepath,Thresholds,cut_score=10):
    """
    1.9.22 
    Required input are a path to a Hmmreport File from HMMER3 and a thresholds dictionary with threshold scores for each HMM
    Returns a list of protein objects for further utilization
    Args
        -Inputfile Report from HMMER3
        -Dictionary Thresholds for HMMs
        
    Return
        -list of Protein objects
    """

    
    protein_dict = {}


    for hmmer_qresult in SearchIO.parse(Filepath,"hmmer3-text"):
        query = hmmer_qresult.id    # Name of HMM without description
        hit_proteinID = ""
        hit_bitscore = 0
        hit_bias = 0
        hit_evalue = 1
        threshold = cut_score
        if query in Thresholds:
            threshold = Thresholds[query] # specific threshold
        
        
        for hit in hmmer_qresult:
            if threshold<hit.bitscore:
                hit_proteinID = hit.id
                hit_bitscore = hit.bitscore
                hit_bias = hit.bias
                hit_evalue = hit.evalue	
                hsp_bitscore = 0
                hsp_start = 0
                hsp_end = 0
                for hsp in hit:
                    #take highest scoring domain as coordinates
                    if hsp_bitscore < hsp.bitscore:
                        hsp_start = hsp.hit_start
                        hsp_end = hsp.hit_end
                        hsp_bitscore = hsp.bitscore
                        
                #print (f"HMM: {query} ProteinID:{hit_proteinID} Hitscore:{hit_bitscore} Bias:{hit_bias} HspScore:{hsp_bitscore} Start:{hsp_start} End:{hsp_end}") Debugging Zeile
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore)
    return protein_dict



def parseHMMreport_below_cutoff_hits(protein_types,Filepath,Thresholds,cut_score=10):
    """
    11.04.2023
    Routine shall find hits that are below the  threshold but that are still significant
    candidate hits are returned in a dictionary with proteinID as key and list as value
    """
    candidate_dict = {}
    for hmmer_qresults in SearchIO.parse(Filepath,"hmmer3-text"):
        query = hmmer_qresults.id
        if query in protein_types:
            hit_proteinID = ""
            hit_bitscore = 0
            hit_bias = 0
            hit_evalue = 1
            threshold = 10
            if query in Thresholds:
                threshold = Thresholds[query] #Upper limit
                cutoff = threshold * 0.1 # Lower limit
                
                for hit in hmmer_qresults:
                    if hit.bitscore<threshold and hit.bitscore > cutoff:
                        hit_proteinID = hit.id
                        hit_bitscore = hit.bitscore
                        hit_bias = hit.bias
                        hit_evalue = hit.evalue
                        hsp_bitscore = 0
                        hsp_start = 0
                        hsp_end = 0
                        for hsp in hit:#take highest scoring domain as coordinates
                            if hsp_bitscore < hsp.bitscore:
                                hsp_start = hsp.hit_start
                                hsp_end = hsp.hit_end
                                hsp_bitscore = hsp.bitscore
                        candidate_dict[hit_proteinID] = [query,hsp_start,hsp_end,hsp_bitscore]
    
    return candidate_dict


def parseGFFfile(Filepath, protein_dict):
    """
    3.9.22
    
    Adds the general genomic features to Protein Objects in a dictionary
    
    Args:
        Filepath - GFF3 formatted file
        protein_dict - Dictionary with key proteinID and value Protein Objects
    Return:
        protein_dict (even though possibly not necessary)
    """
    locustag_pattern = re.compile(r'locus_tag=(\S*?)[\n;]')
    geneID_pattern = re.compile(r'ID=(cds-)?(\S+?)[;\s]')
    
    grep_pattern = "|".join(protein_dict.keys())
    
    try:
        grep_process = subprocess.Popen(['grep', '-E', grep_pattern, Filepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = grep_process.communicate()

        if stderr:
            print("Error in grep process:", stderr)
            return protein_dict
        
        for line in stdout.decode('utf-8').split('\n'):
            if not line:
                continue
            gff = line.split("\t")
            match = geneID_pattern.search(gff[-1])
            if not match:
                continue
            match = match.group(2)
            if match in protein_dict:
                # Add protein
                protein = protein_dict[match]
                
                protein.gene_contig = str(gff[0])
                protein.gene_start = int(gff[3])
                protein.gene_end = int(gff[4])
                protein.gene_strand = str(gff[6])
                locustag = getLocustag(locustag_pattern, line)
                protein.gene_locustag = str(locustag)
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        return protein_dict
    return protein_dict




def getProteinSequence(Filepath, protein_dict):
    """
    3.9.22
    Adds the protein Sequence to the Protein Objects in a dictionary. ProteinIDs of the dictionary have
    to match the header of the .faa file
    
    Args:
        Filepath - fasta formatted amino acid sequence containing file
        protein_dict - Dictionary with key proteinID and value Protein Objects
    Return:
        protein_dict (even though possibly not necessary)    
    """
    reader = None
    try:
        reader = open(Filepath, "r")
        for record in SeqIO.parse(reader, "fasta"):
            if record.id in protein_dict:
                protein = protein_dict[record.id]
                protein.protein_sequence = str(record.seq)
    finally:
        reader.close()
    return protein_dict

def grapProteinSequence(filepath, protein_dict):
    """
    Adds the protein sequence to the Protein Objects in a dictionary.
    ProteinIDs of the dictionary have to match the header of the .faa file.
    
    Args:
        filepath - fasta formatted amino acid sequence containing file
        protein_dict - Dictionary with key proteinID and value Protein Objects
    Return:
        protein_dict (even though possibly not necessary)    
    """
    for protein_id in protein_dict.keys():
        # Use awk to get the sequence until the next '>'
        awk_command = f"awk '/>{protein_id}/ {{flag=1; next}} /^>/ {{flag=0}} flag {{print}}' {filepath}"
        result = subprocess.run(awk_command, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            output = result.stdout
            # Remove all whitespace characters
            sequence = ''.join(output.split())
            protein_dict[protein_id].protein_sequence = sequence
        else:
            print(f"Protein ID {protein_id} not found in {filepath}")
    
    return protein_dict
        

def getLocustag(locustag_pattern,string):
    match = locustag_pattern.search(string)
    if match:
        return match.group(1)
    else:
        return ""









#print("MAKE THRESHOLDS")
#thrs = makeThresholdDict('ttest/Thresholds.txt')
#print("PARSE HMM REPORTS")
#diction = parseHMMreport("ptest/GCA_000006985.1_ASM698v1_genomic.HmmReport",thrs);
#print("PARSE GFF FILE")
#parseGFFfile("ptest/GCA_000006985.1_ASM698v1_genomic.gff",diction)
#print("PARSE AA SEQUENCE")
#getProteinSequence("ptest/GCA_000006985.1_ASM698v1_genomic.faa",diction)

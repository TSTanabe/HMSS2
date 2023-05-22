#!/usr/bin/python




import os
import sys
import argparse
from datetime import datetime

from . import Translation
from . import Search
from . import ParseReports
from . import Csb_finder
from . import AssemblyStatistics
from . import Database
from . import Datasets
from . import myUtil
from . import Output


#get location of script or executable
#for the output module report 

if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
print(__location__)



class HMSSS:
    """
    Objects of this class include the parameters for the program given by the user or by default
    """
    def __init__(self):
        #Parameters
        self.fasta_file_directory = None #bug from prodigal if parentheses are in path
        self.score_threshold_file_directory = __location__+"/src/Thresholds"
        self.pattern_file_directory = __location__+"/src/Patterns"
        self.taxdump = __location__+"/src/taxdump"
        self.library = __location__+"/src/HMMlib"
        self.ziplibrary = __location__+"compressed_HMMlib.zip"
        self.cores = 1
        self.nucleotide_range = 3500
        self.min_completeness = 0.5
        self.threshold_type = 1
        
        self.database_directory = None                       #Directory for database results
        self.result_files_directory = __location__+"/results" #Directory for all results files
        self.logfile = "logfile"
        self.genomize = None
        
        self.synthenic_block_support_detection = False
        
        self.taxonomyfile = None
        self.assembly_stat_file_directory = ""
        self.assembly_stat_mode = "" #Where to retrieve taxonomy information
        self.db_get_genomeIDs = 0   #if database is provided, this will be extended and taxonomy will be updated

        #Control structures
        self.redo_taxonomy = False
        self.redo_search = False
        self.redo_csb = False
        self.redo_csb_naming = False
        self.stat_genomes = False
        
        #Functional dictionaries
        self.finished_genomes = {}
        self.queued_genomes = {}
        self.faa_files = {}
        self.gff_files = {}
        
        self.parameter_required = ("-f","-fasta","-l","-library","-t","-thresholds","-p","--patterns","-c","--cores","-nt","-nucleotide_range","-db","--database","-r","--results","-as","--assembly_stat","-gtdb","-img")
        #Fetching data
        self.fetch = False
        self.fetch_lineage = 0
        self.fetch_taxon = 0
        self.fetch_proteins = []
        self.fetch_keywords = []
        self.keywords_connector = 'OR'
        self.fetch_all = 0
        
        #Alignment and fasta file processing
        self.align_directory = None
        self.merge_fasta = False            # 1 means file list 2: directory
        self.filter_fasta = False           
        self.filter_limits = []
        self.concat_alignment = False       # 1 means file list 2: directory
        self.add_aln_taxonomy = False
        self.process = False
        self.add_genomic_context = False
        
        #iTol dataset creation
        self.dataset_limit_lineage = 0
        self.dataset_limit_taxon = 0
        self.dataset_limit_proteins = []
        self.dataset_limit_keywords = []
        self.dataset_limit_min_cluster_completeness = 0
        self.dataset_divide_sign = '_'
        
        self.dataset_fetch_proteins = []
        self.dataset_fetch_keywords = []
        self.dataset_fetch_fusion = []
        self.dataset_fetch_protkeys = []
        self.dataset_min_cluster_completeness = 0
        self.dataset_combine_binary = 3
        self.dataset_combine_file = 0
        self.create_type_range_dataset = 0 #range dataset for iTol from a protein fasta file
        self.create_gene_cluster_dataset = 0 #gene cluster dataset for iTol from a protein fasta file
        self.dataset = False
        
        self.project_name = "project"
        self.index_db = False

def parse_arguments(arguments):
    """
       15.10.22 ongoing
       Argument parser organizing the different groups of arguments the user can give.
       This includes Search and directory and workflow operators. Also operators for the
       output of data in sequence fasta format, and iTol dataset format
       Also this includes some useful operators for processing the default output files
       regarding taxonomy and sequence sorting/merging and concatenation
    """
    options = HMSSS()
    
    #check if HMMlib has to be unpacked
    if not os.path.isfile(__location__+"/src/HMMlib"):
        if os.path.isfile(__location__+"/src/HMMlib.gz"):
            myUtil.unpackgz(__location__+"/src/HMMlib.gz")
    
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=72,width =200)
    parser = argparse.ArgumentParser(formatter_class=formatter, description = "HMSS2 version 1.0.5 \nSyntax: HMSSS [OPTIONS]",epilog = "Please cite: Tanabe TS, Dahl C. HMS-S-S: A tool for the identification of sulphur metabolism-related genes and analysis of operon structures in genome and metagenome assemblies. Mol Ecol Resour. 2022;22(7):2758-2774. doi:10.1111/1755-0998.13642")
    parser.add_argument('-n','-name',nargs=1, type=str, default=["project"], metavar='<string>', help='Name new project')
    parser.add_argument('-index_db', action='store_true', help='Create index table for database')
    # Resources
    resources = parser.add_argument_group("Search resources")
    resources.add_argument('-f','-fasta',nargs=1, type=myUtil.dir_path, metavar='<directory>', help='Directory to be searched')
    resources.add_argument('-l','-library',nargs=1, type=myUtil.file_path, default=[__location__+"/src/HMMlib"], metavar='<filepath>', help='Filepath to HMM library')
    resources.add_argument('-t','-thresholds',nargs=1, type=myUtil.file_path, default=[__location__+"/src/Thresholds"], metavar='<filepath>', help='Filepath to threshold file')
    resources.add_argument('-p','-patterns',nargs=1, type=myUtil.file_path, default=[__location__+"/src/Patterns"], metavar='<filepath>', help='Filepath to patterns file')
    # Parameters
    parameters = parser.add_argument_group("Search parameters")
    parameters.add_argument('-c','-cores',nargs=1, type=int, default = [2], metavar='<int>', help='Number of cores used by HMMsearch')    
    parameters.add_argument('-nt','-nucleotides',nargs=1, type=int, default = [3500], metavar='<int>', help='Max. number of nucleotides between synthenic genes')    
    parameters.add_argument('-mc','-min_completeness',nargs=1, type=int, default = [0.5], metavar='<int>', help='Min. fraction of existing pattern name or fetch csb')
    parameters.add_argument('-sbs', action='store_false', help='Use synthenic block support for detection')     
    parameters.add_argument('-noise_cut', action='store_true', help='Use noise cutoff scores as threshold')
    parameters.add_argument('-trusted_cut', action='store_true', help='Use trusted cutoff scores as threshold')
    
    # Result space
    results = parser.add_argument_group("Result directory")
    results.add_argument('-r','-results',nargs=1, type=myUtil.dir_path, default=[__location__+"/results"], metavar='<directory>', help='Directory to project')    
    results.add_argument('-db','-database',nargs=1, type=myUtil.file_path, metavar='<filepath>', help='Filepath to existing database')    
    results.add_argument('-gtdb',nargs=1, type=myUtil.file_path, metavar='<filepath>', help='GTDB metadata filepath')
    #results.add_argument('-img',nargs=1, type=myUtil.file_path, metavar='<filepath>', help='IMG metadata filepath')
    results.add_argument('-cutax',nargs=1, type=myUtil.file_path, metavar='<filepath>', help='Custom taxonomy metadata filepath')
    
    #Flow regulators
    flow = parser.add_argument_group("Work step regulation")
    flow.add_argument('-redo_csb', action='store_true', help='Redo the collinear synthenic block prediction')
    flow.add_argument('-redo_csb_naming', action='store_true', help='Redo pattern matching for collinear synthenic blocks')
    flow.add_argument('-redo_search', action='store_true', help='Do not exclude genomes already stored in the database')
    flow.add_argument('-redo_taxonomy', action='store_true', help='Redo the taxonomy')
    
    #Limiter for genomes to account to
    limiter = parser.add_argument_group("Limit output to genomes with conditions")
    limiter.add_argument('-dll','-dlimit_lineage',nargs=1, type=str, metavar='<string>', choices = ['Superkingdom','Phylum','Class','Ordnung','Family','Genus','Species'],\
                            help='Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]')
    limiter.add_argument('-dlt','-dfetch_taxon',nargs=1, type=str, metavar='<string>', help='Taxonomic name e.g. Proteobacteria. Requires -dll')
    limiter.add_argument('-dlp','-dlimit_proteins',nargs='+', type=str, help='Limit fetch to genomes with <protein>')
    limiter.add_argument('-dlk','-dlimit_keyword',nargs='+', type=str, help='Limit fetch to genomes with <keyword>')
    limiter.add_argument('-dlc','-dlmin_cluster_compl',nargs=1, default = [0], type=int, help='Minimal cluster completeness for -dlk keywords. Default 0.75')

    limiter.add_argument('-dtd','-dtaxon_separator',nargs='+', default = '_', type=str, help='Separator for taxonomy information')
    
    #Fetch operator
    operators = parser.add_argument_group("Output sequences with these conditions")
    operators.add_argument('-fl','-fetch_taxonrank',nargs=1, type=str, metavar='<string>', choices = ['Superkingdom','Phylum','Class','Ordnung','Family','Genus','Species'],\
                            help='Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]')
    operators.add_argument('-ft','-fetch_taxon',nargs=1, type=str, metavar='<string>', help='Taxonomic name e.g. Proteobacteria. Requires -fl')
    operators.add_argument('-fd','-fetch_domains',nargs='+', type=str, help='Select only proteins with this domain')
    operators.add_argument('-fk','-fetch_keywords',nargs='+', type=str, help='Select only proteins in gene cluster with this keyword')
    operators.add_argument('-kc','-key_connector',nargs=1, type=str, default = ['OR'], choices = ['AND','OR'], help='Select cluster with keywords connected by AND or OR')
    operators.add_argument('-fetch_all', action='store_true', help='Select everything in DB without limitation')
    operators.add_argument('-stat_genomes', action='store_true', help='Print taxonomy information from database')
    
    
    dataset = parser.add_argument_group("Dataset generation")
    dataset.add_argument('-dfd','-dfetch_domains',nargs='+', type=str, help='Presence/absence matrix of domains')
    dataset.add_argument('-dfk','-dfetch_keywords',nargs='+', type=str, help='presence/absence matrix of keywords')
    dataset.add_argument('-dff','-dfetch_fusions',nargs='+', type=str, help='Presence/absence matrix of multi domain proteins (Syntax A+B C+D)')
    dataset.add_argument('-dfpk','-dfetch_protkeys',nargs='+', type=str, help='presence/absence matrix of proteins with keywords (Syntax dom+key+key )')

    dataset.add_argument('-dmc','-dmin_cluster_compl',nargs='+', default = [0], type=int, help='Minimal cluster completeness for -dfk keywords. Default 0')
    dataset.add_argument('-dcb','-dcombined_binary',nargs=1, default = [3], type=int, help='Calculate n levels of binary permutations')
    dataset.add_argument('-dcf','-dcombined_file',nargs=1, default = [0], type=myUtil.file_path, metavar='<file>', help='Tab seperated file with specified combinations')
    
    #these two define the genomes to be searched, if none provided than all in database will be in iTol dataset and set statistics

    #Alignment and sequence file processing filter_fasta
    process = parser.add_argument_group("Alignment and sequence file processing")
    process.add_argument('-merge_fasta', nargs='+', type=str, metavar='<file> or <directory>', help='Merges two or more sequence files with extention .faa without doublicate')
    process.add_argument('-filter_fasta', nargs=1, type=str, metavar='<file>', help='Filter fasta file by sequence length')
    process.add_argument('-filter_limits', nargs=2, type=int, metavar='<int>', help='Filters sequences by length upper and lower limit')
    process.add_argument('-concat_alignment', nargs='+', type=str, metavar='<file> or <directory>', help='Concatenates alignment files with extention .fasta_aln')
    process.add_argument('-add_taxonomy', nargs=1, type=str, metavar='<file> or <directory>', help='Adds taxonomy to alignment files in <dir>, requires -db with taxonomy')
    process.add_argument('-add_genomic_context', nargs=1, type=myUtil.file_path, metavar='<file>', help='Adds genomic context to sequences from fasta file, requires -db with taxonomy')
    process.add_argument('-create_type_range_dataset', nargs=1, type=myUtil.file_path, metavar='<file>', help='Create protein type range dataset from sequences fasta file, requires -db with taxonomy')
    process.add_argument('-create_gene_cluster_dataset', nargs=1, type=myUtil.file_path, metavar='<file>', help='Create gene cluster dataset from sequences fasta file, requires -db with taxonomy')
    
    #process.add_argument('-dir',nargs=1, type=myUtil.dir_path, metavar='<directory>', help='Directory with fasta and alignment files')

    

    
    namespace = parser.parse_args()
    
    #Assign Namespace variables to options
    if namespace.f:
        options.fasta_file_directory = namespace.f[0]
        options.assembly_stat_file_directory = namespace.f[0]
    options.library = namespace.l[0]
    options.thresholds = namespace.t[0]
    options.patterns = namespace.p[0]
    
    options.cores = namespace.c[0]
    options.nucleotide_range = namespace.nt[0]
    options.min_completeness = namespace.mc[0]
    
    options.result_files_directory = namespace.r[0]
    if namespace.db:
        options.database_directory = namespace.db[0] # Specify database directory
        options.result_files_directory = os.path.dirname(namespace.db[0]) # Set the project path to the directory of the database

    #if namespace.gtdb and namespace.img:
    #    print("WARNING: Only GTDB or IMG can be used to assign taxonomy. Assigning GTDB taxonomy...")
    
    #if namespace.img:
    #    options.assembly_stat_file_directory = namespace.img[0]
    #    options.assembly_stat_mode = "IMG"
    if namespace.gtdb:
        options.assembly_stat_file_directory = namespace.gtdb[0]
        options.assembly_stat_mode = "GTDB"
    if namespace.cutax:
        options.assembly_stat_file_directory = namespace.cutax[0]
        options.assembly_stat_mode = "custom"
        #options.assembly_stat_mode = "GTDB"
        
    if namespace.noise_cut and namespace.trusted_cut:
        print("WARNING: Choose either trusted or noise cutoff.")
        sys.exit()
    
    if namespace.sbs:
        options.synthenic_block_support_detection = True
    
    if namespace.noise_cut:
        options.threshold_type = 3
    if namespace.trusted_cut:
        options.threshold_type = 2

    
    if namespace.redo_search: #redo search for genomes present in database and all other given assemblies
        options.redo_search = True
    if namespace.redo_csb: #redo csb and csb naming
        options.redo_csb = True    
    if namespace.redo_csb_naming: #redo naming of csb only
        options.redo_csb_naming = True
    if namespace.redo_taxonomy:
        options.redo_taxonomy = True
        if namespace.f:
            options.assembly_stat_file_directory = namespace.f[0]
        if options.database_directory == None:
            sys.exit("\nERROR: Database is required to redo taxonomic assginment.")
            raise Exception(f"\nERROR: No writing rights.")
        else:
            options.taxonomyfile = os.path.dirname(namespace.db[0]) + "/taxonomy"
            
            
    if namespace.stat_genomes: #print database statistics
        options.stat_genomes = True
    if namespace.n:
        options.project_name = namespace.n[0]
    if namespace.index_db:
        options.index_db = True
    
    #Limiter options        do not start fetch by itself    
    if namespace.dll and namespace.dlt:
        options.dataset_limit_lineage = namespace.dll
        options.dataset_limit_taxon = namespace.dlt
    if namespace.dlp:
	    options.dataset_limit_proteins = namespace.dlp
    if namespace.dlk:
        options.dataset_limit_keywords = namespace.dlk
    if namespace.dlc:
        options.dataset_limit_min_cluster_completeness = namespace.dlc
    if namespace.dmc:
        options.dataset_min_cluster_completeness = namespace.dmc
    if namespace.dtd:
        options.dataset_divide_sign = namespace.dtd 
    
    #Sequence retrieve      do start fetch
    if namespace.fl and namespace.ft:
        options.dataset_limit_lineage = namespace.fl[0]
        options.dataset_limit_taxon = namespace.ft[0]
        options.fetch = True
    if namespace.fd:
        options.fetch_proteins = namespace.fd #class list
        options.fetch = True
    if namespace.fk:
        options.fetch_keywords = namespace.fk #class list
        options.fetch = True
    if namespace.kc:
        options.keywords_connector = namespace.kc[0] #AND or OR
    if namespace.fetch_all:
        options.fetch_all = 1
        options.fetch = True
    
    #Dataset generation without sequence output     do start by itself
    if namespace.dfd:
        options.dataset_fetch_proteins = namespace.dfd
        options.dataset = True
    if namespace.dfk:
        options.dataset_fetch_keywords = namespace.dfk
        options.dataset = True   
    if namespace.dff:
        options.dataset_fetch_fusion = namespace.dff
        options.dataset = True
    if namespace.dfpk:
        options.dataset_fetch_protkeys = namespace.dfpk
        options.dataset = True   

    if namespace.dcb:
        options.dataset_combine_binary = namespace.dcb[0]
    if namespace.dcf:
        options.dataset_combine_file = namespace.dcf
       
        

    if namespace.merge_fasta:        
        if os.path.isdir(namespace.merge_fasta[0]):
            options.align_directory = namespace.merge_fasta[0]
            options.merge_fasta = 2
            options.process = True
            
        else:
            for path in namespace.merge_fasta:
                myUtil.file_path(path)
            options.align_directory = namespace.merge_fasta
            options.merge_fasta = 1
            options.process = True
    
    if namespace.concat_alignment:
        if os.path.isdir(namespace.concat_alignment[0]):
            options.align_directory = namespace.concat_alignment[0]
            options.concat_alignment = 2
            options.process = True
            
        else:
            for path in namespace.concat_alignment:
                myUtil.file_path(path)
            options.align_directory = namespace.concat_alignment
            options.concat_alignment = 1
            options.process = True
        
        
    if namespace.filter_fasta:
        options.filter_fasta = namespace.filter_fasta[0]
        options.process = True
    if namespace.filter_limits:
        options.filter_limits.append(namespace.filter_limits[0])
        options.filter_limits.append(namespace.filter_limits[1])
    
    if namespace.add_taxonomy:            
        options.add_aln_taxonomy = namespace.add_taxonomy[0]
        options.process = True    
    if namespace.add_genomic_context:
        options.add_genomic_context = namespace.add_genomic_context[0]
        options.process = True    
    if namespace.create_type_range_dataset:
        options.create_type_range_dataset = namespace.create_type_range_dataset[0]
        options.process = True    
    if namespace.create_gene_cluster_dataset:
        options.create_gene_cluster_dataset = namespace.create_gene_cluster_dataset[0]
        options.process = True            
    return options
    

#############
####   Subroutines for preparation of the result folder
#############
def create_project(directory,projectname="project"):
    now = datetime.now()
    timestamp = datetime.timestamp(now)
    directory = directory + "/" + str(timestamp) + "_" + projectname
    try:
        os.mkdir(directory)
    except:
        sys.exit("\nERROR: Could not create result directory. No writing rights.")
        raise Exception(f"\nERROR: No writing rights.")
    return directory
    
def create_logfile(directory):
    logfile = directory +"/logfile"
    #print(logfile)
    with open (logfile,"w") as writer:
        writer.write("Date\tgenomeID\tLibrary\tdatabase\tfile\n")
    
    return logfile

def create_genomize_dir(directory):
    os.mkdir(directory+"/per_genome")
    return  directory+"/per_genome"

def read_logfile(filepath):
    used_files = set()
    #f = open(filepath,"r")
    #print(f.readlines())
    with open(filepath,"r") as reader:
        line = reader.readline()  # skip header
        for line in reader.readlines():
            line = line.replace("\n","")
            array = line.split("\t")
            #used_files[array[1]] = 1    # genomeID into array
            used_files.add(array[1])
    return used_files

def write_logfile_entry(genomeID,options):
    array = [str(datetime.now()),genomeID,myUtil.getFileName(options.library),options.database_directory,options.faa_files[genomeID]]
    string = "\t".join(array)+"\n"
    with open (options.logfile,"a") as writer:
        writer.write(string)
    writer.close()
    return

def create_taxonomyfile(directory):
    taxon = directory +"/taxonomy"
    #print(logfile)
    with open (taxon,"w") as writer:
                       #filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain
                        writer.write("genomeID\tsuperkingdom\tclade\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\ttaxid\tbiosample\tbioproject\tgenbank\trefseq\tcompleteness\tcontamination\ttypestrain\n")
    return taxon

  
    
def prepare_result_space(options,project="project"):
    """
        6.10.22
        Args:
            directory   directory for the result files to be stored
        Return:
            path        directory path for the result files
            continue    boolean if old project should be continues
        Options taken are
        Flowchart for this routine
            Checks if user input directory is already present
             if not 
                try to create this one
             if yes 
                check if database is already present
                 if no
                    create new one
                 if yes
                    continue writing on this existing database
            when database is created follow with the other report files        
            
    """
    myUtil.print_header("\nPreparing space for the results")
    directory = options.result_files_directory
    database = options.database_directory
    logfile = 0
    
    if os.path.isdir(directory):   #check if directory is existing
        #if exists check for database, without database create new one
        if options.database_directory == None and os.path.isfile(directory+"/database.db"):
            options.database_directory = directory+"/database.db" #old DB found
            print(f"Database: Extending database {options.database_directory}")
            options.db_get_genomeIDs = 1                          #get genomeIDs from database without taxonomy
        elif not options.database_directory == None:
            print(f"Database: Extending database {options.database_directory}")
            options.db_get_genomeIDs = 1
        elif not os.path.isfile(directory+"/database.db"): #dir does not contain any DB, so it is a result folder and new projects are added here

            directory = create_project(directory,project) #create new project in the folder
            
            options.database_directory = directory+"/database.db" #no old DB found creating new one
            Database.create_database(options.database_directory)
            print(f"Database: New database {options.database_directory}")
            
            
        #check for HMMlibrary    
        #check for logfile and continue if present

        if os.path.isfile(directory + "/logfile"):
            print(f"Project: Continue project {directory}")
            options.logfile = directory + "/logfile"
        
        else:
            options.logfile = create_logfile(directory)
        
        if options.taxonomyfile == None and os.path.isfile(directory+"/taxonomy"):
            options.taxonomyfile = directory+"/taxonomy"
        elif not os.path.isfile(directory+"/taxonomy"):
            print("Writing taxonomy to: "+directory+"/taxonomy")
            options.taxonomyfile = create_taxonomyfile(directory)
        
        try:
            options.genomize = create_genomize_dir(directory)
        except:
            options.genomize = directory+"/per_genome"
          
            
    else:                           #if not existing create one
        try:
            os.mkdir(directory)
        except:
            sys.exit("\nERROR: Could not create result directory. Missing writing rights.")
            raise Exception(f"\nERROR: No writing rights.")
    
        directory = create_project(directory,project)
        print(f"Project: New project {directory}")
        options.genomize = create_genomize_dir(directory)
        options.logfile = create_logfile(directory)
        options.taxonomyfile = create_taxonomyfile(directory)
        if options.database_directory == None:
            options.database_directory = directory+"/database.db" #no old DB found creating new one
            print(f"Database: New database {options.database_directory}")
            Database.create_database(options.database_directory)
        else:
            print(f"Created new project but using database {options.database_directory}")
            print(f"There was no logfile found. Searching all given files in {options.fasta_file_directory}")
    
      
def queue_files(options):
    """
        9.10.22
        Args:
            directory   directory for the result files to be stored
            finished    set with genome identifiers already processed
            options     current options object
        Operation:
            collect all zip protein fasta files and unzipped protein fasta files
            collect all corresponding gff files and check for existance of this file
            if both files present
            get the genome identifiers
    """
    print("\nFilling the queue with faa files --")
    queue = set()
    faa_files = {}
    gff_files = {}

    #For all zipped .faa files find corresponding zipped or unzipped gff. Fill hashes with zipped format only
    all_zip_FaaFiles = myUtil.getAllFiles(options.fasta_file_directory,".faa.gz")
    for zip_FaaFile in all_zip_FaaFiles:
        faa_file = myUtil.removeExtension(zip_FaaFile)
        file_name = myUtil.removeExtension(faa_file)
        #print(file_name+".gff.gz")
        if os.path.isfile(file_name+".gff.gz"):
            gff_file = file_name+".gff.gz"
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID) #GenomeIDs that are queued
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        elif os.path.isfile(file_name+".gff"):
            gff_file = myUtil.packgz(file_name+".gff")
            myUtil.unlink(file_name+".gff")
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        else:
            print(f"\tWARNING: No corresponding .gff file found for {zip_FaaFile}")
    
    #For all unzipped .faa files find corresponding zipped or unzipped gff. Fill hashes with zipped format only
    all_FaaFiles = myUtil.getAllFiles(options.fasta_file_directory,".faa")
    for FaaFile in all_FaaFiles:
        zip_FaaFile = myUtil.packgz(FaaFile)
        myUtil.unlink(FaaFile)
        faa_file = myUtil.removeExtension(zip_FaaFile)
        file_name = myUtil.removeExtension(faa_file)
        if os.path.isfile(file_name+".gff.gz"):
            gff_file = file_name+".gff.gz"
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        elif os.path.isfile(file_name+".gff"):
            gff_file = myUtil.packgz(file_name+".gff")
            myUtil.unlink(file_name+".gff")
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        else:
            print(f"\tWARNING: No corresponding .gff file found for {zip_FaaFile}")

    #compare two sets
    options.queued_genomes = queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    
    
    #compare with database leaving out all genomes already searched previously
    if options.redo_search: #redo the search for all
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        options.finished_genomes = {}
        for genomeID in genomeIDs:
            options.finished_genomes[genomeID] = 1 #requires for extending the protein dict in csb assignment with all old results
    else: #leave out all genomeIDs in the database
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        for genomeID in genomeIDs:
            if genomeID in faa_files.keys():
                print(f"\tFound assembly {genomeID} in database leaving out {myUtil.getFileName(faa_files[genomeID])}")
                del faa_files[genomeID]
                del gff_files[genomeID]
                options.queued_genomes.remove(genomeID)
    
    print(f"Queued {len(options.queued_genomes)} for processing")
    if len(options.queued_genomes) == 0:
        print("There were 0 genomes queued. Use -redo_search option if genomes are already present in database")
    print("\nFilling the queue with faa files -- ok")
    

    
def translate_and_transcripe_fasta(directory):
    """
        9.10.22
        Args:
            directory   fasta file directory
        Return:
            nothing
        Process:
            takes all fna files without gff (either .gz or not) and uses prodigal
            for translation to faa and gff
            
            takes all faa without gff (either .gz or not) and transcripes
            to gff if it is a prodigal derived file            
    """
    myUtil.print_header("\nProkaryotic gene recognition and translation via prodigal")
    Translation.translation(directory)
    Translation.transcription(directory)
  
            
def search_and_csb(options):
    """
        9.10.22
        Args:
            options object
        Return:
            nothing
            
        Reads threshold file and stores in a hash
        Reads pattern files and stores in a hash
    """
    myUtil.print_header("\nPerforming protein annotation and collinear synthenic block prediction")
    # prepare thresholds and patterns
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file_directory,options.threshold_type)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.pattern_file_directory)
    
    for index,genomeID in enumerate(options.queued_genomes):
        now = datetime.now()
        print(f"{now} Processing assembly {index+1} of {len(options.queued_genomes)}") #print control sequence output
        #Prepare names
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        file_name = myUtil.removeExtension(faa_file)
        hmm_report = myUtil.getReportName(faa_file)
        
        # Search
        Search.HMMsearch(faa_file,options.library,options.cores)
        
        # Parse results & write to database
        protein_diction = ParseReports.parseHMMreport(hmm_report,score_threshold_diction)
        ParseReports.parseGFFfile(gff_file,protein_diction)
        
        # Update the protein dict with previous results
        if genomeID in options.finished_genomes.keys():
            protein_diction.update(Database.fetch_protein_dict(options.database_directory,genomeID)) #dict1.update(dict2)
            
        # Synteny and cluster naming    
        cluster_diction = Csb_finder.find_syntenicblocks(genomeID,protein_diction,options.nucleotide_range)
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_diction,options.min_completeness)
        
        # Increase detection sensitivity by non-homologous synteny detection criterum
        if options.synthenic_block_support_detection:
            missing_protein_types, missing_proteins_list_dict = Csb_finder.find_csb_pattern_difference(csb_patterns_diction,csb_pattern_names,cluster_diction,3)
            candidate_proteins_dict = ParseReports.parseHMMreport_below_cutoff_hits(missing_protein_types,hmm_report,score_threshold_diction)
            Csb_finder.synteny_completion(gff_file,protein_diction,cluster_diction,candidate_proteins_dict,missing_proteins_list_dict,options.nucleotide_range)
        
        
        ParseReports.getProteinSequence(faa_file,protein_diction)
        Database.insert_database_genomeID(options.database_directory,genomeID)
        Database.insert_database_protein(options.database_directory,genomeID,protein_diction)
        Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)

        """11.4.23 deprecated
        if genomeID in options.finished_genomes.keys(): #present database is extended, protein_dict might be incomplete at the start
            Database.extend_database_protein(options.database_directory,genomeID,protein_diction)

            protein_diction = Database.fetch_protein_dict(options.database_directory,genomeID) #get all proteins of the altered genome and find/name csb
            cluster_diction = Csb_finder.find_syntenicblocks(genomeID,protein_diction,options.nucleotide_range)
            Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_diction,options.min_completeness)
            Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
        
        else: #new genome was searched, protein_dict likely complete
            Database.insert_database_genomeID(options.database_directory,genomeID)
            Database.insert_database_protein(options.database_directory,genomeID,protein_diction)

            cluster_diction = Csb_finder.find_syntenicblocks(genomeID,protein_diction,options.nucleotide_range)
            Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_diction,options.min_completeness)
            Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
        """
        
        # Write output files
        
        Output.output_genome_report(options.genomize+"/"+genomeID,protein_diction,cluster_diction)
        
        write_logfile_entry(genomeID,options)

        myUtil.unlink(faa_file)
        myUtil.unlink(gff_file)
    print("\nFinished")
    
def csb(options):
    """
        18.10.22
        Args:
            options object
        Return:
            nothing
    """
    myUtil.print_header(f"\nRedo collinear syntenic block assignment")
    print("Reading pattern")
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.pattern_file_directory)
    print("Collecting genomeIDs")
    genomeIDs = Database.fetch_genomeIDs(options.database_directory)
    lim = str(len(genomeIDs))
    for index,genomeID in enumerate(genomeIDs):
        print(f"Genome number {index+1} of {lim}", end="\r")
        protein_diction = Database.fetch_protein_dict(options.database_directory,genomeID)
        #print("\t Looking for syntenic blocks")
        cluster_diction = Csb_finder.find_syntenicblocks(genomeID,protein_diction,options.nucleotide_range)
        #print("\t Naming syntenic blocks")
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_diction,options.min_completeness)
        #print("\t Insert into database")
        Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
    print("Finished csb naming\n")

def csb_names(options):
    """
        09.11.22
        Args:
            options object
        Return:
            nothing
    """
    
    myUtil.print_header(f"\nRedo collinear syntenic block assignment")
    print("Reading pattern")
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.pattern_file_directory)
    print("Collecting genomeIDs")
    genomeIDs = Database.fetch_genomeIDs(options.database_directory)
    lim = str(len(genomeIDs))
    print(f"Reiterate pattern naming for {lim} genomes")
    for index,genomeID in enumerate(genomeIDs):
        print(f"Genome number {index+1} of {lim}")
        cluster_diction = Database.fetch_cluster_dict(options.database_directory,genomeID)
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_diction,options.min_completeness)
        Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
    print("\nFinished csb naming")    
    
def collect_taxonomy_information(options):
    """
     12.10.22 no comment
     19.2.23
      can either be invoked by default when a new database is created and searched
      or as a redo with existing database
      
      converts first the input format into a HMSSS compatible format and then adds it 
      to the DB
     22.5.23
        NCBI deactivated, trouble with bioperl script in 
    """  
    if not options.assembly_stat_mode:
        return
        
    myUtil.print_header(f"\nTaxonomy assignment via {options.assembly_stat_mode} taxonomy")
    print("Assigning assembly statistics --",end="\r")
    if options.redo_taxonomy:
        #get all genomeIDs regardless of pre-existing taxonomy
        options.queued_genomes = Database.fetch_genomeIDs(options.database_directory)
    elif options.db_get_genomeIDs:
        #get all genomIDs were any taxonomy level information is missing and all queued genomes for assignment
        DBset = Database.fetch_genomeID_without_phylogeny(options.database_directory)
        options.queued_genomes = set.union(options.queued_genomes,DBset)

    
    #database to insert into the information => options.database_directory
    #genomeIDs to search taxonomic information for => options.queued_genomes
    #directory with the assembly metadata => options.assembly_stat_file_directory
    
    #For NCBI it has to be checked if directory was given as it is default and with -redo_taxonomy
    # it could also mean that only the DB genomeIDs should be redone
    #if options.assembly_stat_mode == "NCBI" and os.path.isdir(options.assembly_stat_file_directory):
    #    AssemblyStatistics.taxdump_perl_script(__location__+"/bin/Assemblystats",options.assembly_stat_file_directory,\
    #                                           options.taxonomyfile,options.taxdump)
    #    print("Assigning assembly statistics -- ok")  
    if options.assembly_stat_mode == "GTDB":
        AssemblyStatistics.parse_gtdb_metadata(options.assembly_stat_file_directory,\
                                               options.queued_genomes,options.taxonomyfile)
        print("Assigning assembly statistics -- ok")  
    #elif options.assembly_stat_mode = "IMG": %TODO
        #AssemblyStatistics.IMG_perl_script(options.assembly_stat_file_directory,options.taxonomyfile)
    elif options.assembly_stat_mode == "custom": # custom set means direct entry of the file 
        AssemblyStatistics.parse_custom_metadata(options.assembly_stat_file_directory,\
                                               options.queued_genomes,options.taxonomyfile)
        print("Assigning assembly statistics -- ok")
    
      
    print("Writing taxonomy to file --", end="\r")
    AssemblyStatistics.insert_database_assembly_statistics(options.database_directory,options.taxonomyfile,options.queued_genomes)
    print("Writing taxonomy to file -- ok")

def output_operator(options):
    """
        17.10.22
        07.11.22
        10.11.22
        Args:
            options object
        Return:
            nothing
        Output:
            File    fasta formatted file
            File    metadata file
    

    """

    myUtil.print_header(f"\nOutput from database")
    
    #Set directory
    date = datetime.now()
    database = options.database_directory
    path = myUtil.getPath(options.database_directory)
    directory = path+str(date)+"_dataset/"
    os.mkdir(directory)
    #path1 = directory+str(options.fetch_lineage)+"_"+str(options.fetch_taxon)+"_"+str(options.fetch_proteins)+"_"+str(options.fetch_keywords)+".txt"
    #path2 = directory+str(options.fetch_lineage)+"_"+str(options.fetch_taxon)+"_"+str(options.fetch_proteins)+"_"+str(options.fetch_keywords)
    path1 = directory+"1_hit_table.txt"
    path2 = directory+"2_Protein"
    
    
    #Limit the genomeIDs and fetch the lineage
    print("Collecting taxonomy lineage information")
    taxon_diction = Output.fetch_limiter_data(database,options.dataset_limit_lineage,\
                                              options.dataset_limit_taxon,\
                                              options.dataset_limit_proteins,\
                                              options.dataset_limit_keywords,\
                                              options.dataset_limit_min_cluster_completeness,\
                                              options.dataset_divide_sign)
    print("\nCollecting taxonomy lineage information -- ok\n")
    #Fetch protein/cluster diction                      
    print("Collecting protein sequences")                        
    protein_diction,cluster_diction = Output.fetch_bulk_data(options.database_directory,\
                                                             options.fetch_proteins,options.fetch_keywords,\
                                                             taxon_diction,options.min_completeness)
    print("Collecting proteins sequences -- ok\n")
    
    
    print("Writing output to disc\n")
    Output.output_genome_report(path1,protein_diction,cluster_diction,taxon_diction)    #Output report
    files = Output.output_distinct_fasta_reports(path2,protein_diction,cluster_diction) #Output per domain
    Output.singletons(directory,files)                                                  #Output singleton per genome and doublicates
    
    files = Output.output_combined_keyword_fasta_reports(directory,files,options.fetch_keywords)
    Output.singletons(directory,files)
    myUtil.clean_empty_files(directory)
    protein_diction.clear()
    cluster_diction.clear()
    print(f"\nWrote output to :\n{path1} \nand\n{path2}")
    
    
    #Fetch dataset dictionary and headers       
    print("\n\nCollecting presence/absence of sequences per genome\n")
    dataset_dict,headers = Datasets.fetch_binary_dataline(database,taxon_diction,\
                                                          options.fetch_proteins,options.fetch_keywords,\
                                                          list(),list(),\
                                                          options.dataset_min_cluster_completeness)
    print("Collecting presence/absence of sequences per genome -- ok")
    print("Writing datasets on disk space")
    #Output iTol dataset and Matrices
    iTol_data_dict = Datasets.iTol_lineage_to_binary_dataline(taxon_diction,dataset_dict)
    Datasets.iTol_binary_file_dataset(directory,iTol_data_dict,headers)
    combined_dataset_dict,combined_headers = Datasets.dataline_combinations(dataset_dict,headers,options.dataset_combine_file,options.dataset_combine_binary)
    Datasets.dataset_statistics(database,directory,combined_headers,combined_dataset_dict)
    print("Writing datasets on disk space -- ok")
    print("Finished")
    
    
def dataset_generation(options):
    """
        06.11.22
        Args:
            options object with directory and mode
                -database
                -taxonomy/protein/keyword limiter
                -matrix proteins/keywords
                -minimal cluster completeness
                -level of permutations
        Return:
            nothing
        Output:
            File    for each taxonomy level presence absence matrix of proteins and keywords
            File    iTol dataset of the 
            
        
    """
    date = datetime.now()
    database = options.database_directory
    path = myUtil.getPath(options.database_directory)
    directory = path+str(date)+"_dataset/"
    os.mkdir(directory)


    print("Collecting taxonomy lineage information")
    taxon_dict = Output.fetch_limiter_data(database,options.dataset_limit_lineage,\
                                              options.dataset_limit_taxon,\
                                              options.dataset_limit_proteins,\
                                              options.dataset_limit_keywords,\
                                              options.dataset_limit_min_cluster_completeness,\
                                              options.dataset_divide_sign) #return hold taxonomy data
    print("\nCollecting taxonomy lineage information -- ok\n")
    #Datasets.iTol_taxonomy_range_dataset(directory,taxon_dict)     Commented out because range datasets have far too many possibilites to just be made here by incidence.
    print("Collecting presence/absence of domains and keywords per genome")                                              
    dataset_dict,headers = Datasets.fetch_binary_dataline(database,taxon_dict,\
                                                          options.dataset_fetch_proteins,options.dataset_fetch_keywords,\
                                                          options.dataset_fetch_fusion,options.dataset_fetch_protkeys,\
                                                          options.dataset_min_cluster_completeness) #return genomeID => binary line tab separated
    print("Collecting presence/absence of domains and keywords per genome -- ok")
    print("Writing datasets on disk space")
    iTol_data_dict = Datasets.iTol_lineage_to_binary_dataline(taxon_dict,dataset_dict)
    Datasets.iTol_binary_file_dataset(directory,iTol_data_dict,headers) #iTol binary dataset output
    combined_dataset_dict,combined_headers = Datasets.dataline_combinations(dataset_dict,headers,options.dataset_combine_file,options.dataset_combine_binary) #combined presence absence
    Datasets.dataset_statistics(database,directory,combined_headers,combined_dataset_dict)
    print("Writing datasets on disk space -- ok")
    print("Finished")

    
    

    
def output_statistics(options):    
    Database.fetch_genome_statistic(options.database_directory)
    Database.fetch_bulk_statistic(options.database_directory)
    
    
def process_operator(options):
    """
        03.11.22
        Args:
            options object with directory and mode
            directory   is essential
        Return:
            nothing
        Output:
            File    fasta formatted file
            
        This process shall merge fasta files, concat alignments, and add taxonomic information to fasta/alignment files    
    """
    myUtil.print_header(f"\nProcessing sequence files")
    #if not options.align_directory and not options.add_aln_taxonomy and not options.filter_fasta and not options.add_genomic_context:
    #    print("WARNING: Missing alignment file directory, please use -dir argument")
    #    return
    
#Merge fasta files
    if options.merge_fasta:
        if options.merge_fasta == 2:
            fasta_files = myUtil.getAllFiles(options.align_directory,".faa")
            output = options.align_directory+"/concat.fasta"
        elif options.merge_fasta == 1:
            fasta_files = options.align_directory
            path = myUtil.getPath(options.align_directory[0])
            output = path+"/concat.faa"
            
        Output.merge_fasta(output,fasta_files)
        print(f"\nWrote merged fasta files without doublicates to:\n {output}")
#Concat alignment files    
    if options.concat_alignment:
        if options.concat_alignment == 2:
             fasta_files = myUtil.getAllFiles(options.align_directory,"fasta_aln")
             output = options.align_directory+"/concat.fasta_aln"
        elif options.concat_alignment == 1:
            fasta_files = options.align_directory
            path = myUtil.getPath(options.align_directory[0])
            output = path+"/concat.faa"
            
        Output.concat_alignments(output,fasta_files)
        print(f"\nWrote concatenated alignments to:\n {output}")    

#Filter fasta files by length    
    if options.filter_fasta:
        if os.path.isfile(options.filter_fasta): 
            Output.filter_length_fasta(options.filter_fasta,options.filter_limits[0],options.filter_limits[1])
    
#Add taxonomy information    
    if options.add_aln_taxonomy:
        if not options.database_directory:
            print("WARNING: Missing database to assign taxonomy, please use -db argument")
            return  
        if os.path.isdir(options.add_aln_taxonomy): 
        #if directory was provided
        #name all sequences and concat all files
            fasta_files = myUtil.getAllFiles(options.add_aln_taxonomy,".faa")
            for fasta in fasta_files:
                Output.add_taxonomy(options.database_directory,fasta)

            print("Concatenating named fasta files")
            fasta_files = myUtil.getAllFiles(options.add_aln_taxonomy,".taxon_names")    
            output = options.add_aln_taxonomy+"/concat.fasta_aln"
            Output.concat_alignments(output,fasta_files)
            print(f"\nWrote concatenated alignments to:\n {output}")
            
        elif os.path.isfile(options.add_aln_taxonomy):
        #if single file was provided
        #name all sequences
            Output.add_taxonomy(options.database_directory,options.add_aln_taxonomy)
            
        else:
            print(f"ERROR: {options.add_aln_taxonomy} is not a file or directory")

#Get a textfile with genomic context based on provided sequence fasta file    
    if options.add_genomic_context:
        if not options.database_directory:
            print("WARNING: Missing database to assign taxonomy, please use -db argument")
            return 
        Output.add_genomic_context(options.database_directory,options.add_genomic_context)

#Get a iTol dataset file with gene cluster dataset    
    if options.create_gene_cluster_dataset:
        if not options.database_directory:
            print("WARNING: Missing database to assign taxonomy, please use -db argument")
            return 
        directory = os.path.dirname(options.create_gene_cluster_dataset)
        Datasets.iTol_domain_dataset(directory,options.database_directory,options.create_gene_cluster_dataset,options.dataset_divide_sign)

#Get a iTol dataset file with range data per protein type
    if options.create_type_range_dataset:
        if not options.database_directory:
            print("WARNING: Missing database to assign taxonomy, please use -db argument")
            return 
        directory = os.path.dirname(options.create_type_range_dataset)
        Dataset.iTol_range_dataset(directory,options.database_directory,options.create_type_range_dataset,options.dataset_divide_sign)
    
    
    


def main(args=None):
    #1
    options = parse_arguments(args)
    if options.redo_csb:
        #8
        csb(options)
    elif options.redo_csb_naming:
        #8
        csb_names(options)
    elif options.redo_taxonomy:
        #7
        collect_taxonomy_information(options)
    elif options.fasta_file_directory:
        #2
        translate_and_transcripe_fasta(options.fasta_file_directory)
        #3
        prepare_result_space(options,options.project_name)
        #4
        queue_files(options)    
        #5 & 6
        search_and_csb(options)
        #7
        collect_taxonomy_information(options)
    
    
    
    
    
    
    
    if options.index_db:
        #12
        Database.index_database(options.database_directory)
    if options.fetch:
        #14
        output_operator(options)
    if options.stat_genomes:
        #9
        output_statistics(options)
       
    if options.process:
        #10 special because not in pipeline
        # file/alignment concat utilities
        process_operator(options)
        
    if options.dataset:
        #11 light version of #14
        #on big datasets faster when separately invoked
        dataset_generation(options) 
    

    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args) #calls the main method of __main__


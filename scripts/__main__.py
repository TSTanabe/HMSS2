#!/usr/bin/python
import os
import sys
import argparse
from datetime import datetime

from . import Translation
from . import Csb_finder
from . import Csb_Mp_Algorithm
from . import Csb_redo
from . import Csb_cluster
from . import Database
from . import Datasets
from . import myUtil
from . import Output
from . import Project
from . import Queue
from . import ParallelSearch
from . import BulkSearch
from . import Processing


# - sollen unterschiedlicher feld informationen trennen
# _ sollen namentrennungen sein, bzw indices
#get location of script or executable
#for the output module report 
#TODO fc option sollte auch die patterns adden, die über den pattern file eingetragen werden
#TODO print a database description file from database (otherwise it is too complex for a short task) should include the 
#TODO syntenycompletion is not properly working


if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#print(__location__)



class HMSSS:
    """
    Objects of this class include the parameters for the program given by the user or by default
    """
    def __init__(self):
        #Parameters
        self.execute_location = __location__
        
        #queue functional dictionaries
        self.finished_genomes = {}
        self.queued_genomes = {}
        self.faa_files = {}
        self.gff_files = {}
        
        #csb prediction dereplication
        self.redundant = 0
        self.non_redundant = 0
        self.redundancy_hash = dict()
        self.computed_Instances_dict = 0 #for json of computed csb
        
        #Limiter
        self.limiter = False
        
        #Fetching data
        self.fetch = False
        
        #Alignment and fasta file processing
        self.process = False
        
        self.project_name = "project"
        self.index_db = False
        
        self.csb_name_prefix = "csb-" #prefix of clusterIDs determined by csb finder algorithm
        self.csb_name_suffix = "." #suffix of clusterIDs determined by csb finder algorithm
        
        self.genomeID_divider = '___' #dividing sign between genomeID and proteinID, first part will be taken as genomeID
        
        self.synthenic_block_support_detection = False #TODO dummy option before implementation
        
def parse_arguments(arguments):
    """
       Argument parser organizing the different groups of arguments the user can give.
       This includes Search and directory and workflow operators. Also operators for the
       output of data in sequence fasta format, and iTol dataset format
       Also this includes some useful operators for processing the default output files
       regarding taxonomy and sequence sorting/merging and concatenation
    """
    options = HMSSS()
    
    #Define the Argpase headers and footers
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=72,width =200)
    description = ""
    epilog = "Please cite: Tanabe TS, Dahl C. HMSS2: An advanced tool for the analysis of sulphur metabolism, including organosulphur compound transformation, in genome and metagenome assemblies. Mol Ecol Resour. 2023;23(8):1930-1945. doi:10.1111/1755-0998.13848"
    prog = "HMSS2 version 1.0.7"
    usage = "HMSSS [OPTIONS] (*) dependent on the -db option"
    
    parser = argparse.ArgumentParser(formatter_class=formatter, description = description, epilog = epilog, usage = usage)
    parser.add_argument('-f', dest= 'fasta_file_directory', type=myUtil.dir_path, default = None, metavar='<directory>', help='Directory to be searched')
    parser.add_argument('-glob_report', dest='glob_report', type=myUtil.dir_path, metavar='<directory>', help = 'Directory with hmmreports. Each report with one HMM queried against the concatenated genomes.')
    parser.add_argument('-r', dest='result_files_directory', type=myUtil.dir_path, metavar='<directory>', default = __location__+"/results", help='Directory for the result files')
    parser.add_argument('-n', dest='name', type=str, default="project", metavar='<string>', help='Name new project')
    parser.add_argument('-s', dest='stage', type=int, default = 0, choices= [0,1,2,3,4,5,6,7,8,9],help='Start at stage')
    
    # Informations
    information = parser.add_argument_group("Information on resources")
    #information.add_argument('-index_db', action='store_true',help='Create index table for database') #moved the index database routine into the output routine as default
    information.add_argument('-stat_keywords', action='store_true', help='Print patterns for keyword naming')
    information.add_argument('-stat_csb', action='store_true', help='Print automatically found csbs')
    information.add_argument('-stat_genomes', action='store_true', help='Print taxonomy information from database')
    
      
    
 
    # Resources
    resources = parser.add_argument_group("Search resources")
    resources.add_argument('-l', dest= 'library', type=myUtil.file_path, default=__location__+"/src/HMMlib", metavar='<filepath>', help='Filepath to HMM library')
    resources.add_argument('-c', dest= 'cores' , type=int, default = 4, metavar='<int>', help='Allocated CPU cores. Default: 4')
    resources.add_argument('-t', dest= 'score_threshold_file', type=myUtil.file_path, default=__location__+"/src/Thresholds", metavar='<filepath>', help='Filepath to threshold file')
    resources.add_argument('-p', dest= 'patterns_file' , type=myUtil.file_path, default=__location__+"/src/Patterns", metavar='<filepath>', help='Filepath to patterns file')
    resources.add_argument('-db',dest='database_directory', type=myUtil.file_path, metavar='<filepath>', help='Filepath to existing database')
    resources.add_argument('-tax', dest='taxonomy_file', type=myUtil.file_path, metavar='<filepath>', help ='Filepath to taxonomy tsv file')
    resources.add_argument('-clean', dest='clean_reports', action='store_true', help = 'Clean up any pre-existing HMMreport file')
    resources.add_argument('-no_reports', dest='individual_reports', action='store_false', help = 'Do not write individual files per genome')
    resources.add_argument('-glob_search', dest='glob_search', action='store_true', help = 'Input reports from search of concatenated genomes')
    resources.add_argument('-glob_chunks', dest='glob_chunks', type=int, default=5000, metavar='<int>', help='Chunk size for parsing results from glob before entering into database')
    
    # Parameters
    parameters = parser.add_argument_group("Search parameters")
    parameters.add_argument('-mc', dest= 'min_completeness', type=int, default = 0.5, metavar='<int>', help='Minimal fraction of predefined csb to be recognized')
    #parameters.add_argument('-sbs', dest= 'synthenic_block_support_detection', action='store_true', help='Use synthenic block support for detection') TODO this has to be added again but only for predefined patterns
    parameters.add_argument('-cut_type', dest='threshold_type', type=int, default = 1, metavar='<int>', choices = [1,2,3], help='Choice of cutoff: 1 optimized; 2 trusted; 3 noise; Default: 1')
    parameters.add_argument('-cut_score', dest='thrs_score', type=int, default = 10, metavar='<int>', help='Default cutoff score. Default: 10')
    
    # Result space
    csb = parser.add_argument_group("Collinear syntenic block prediction")
    csb.add_argument('-nt', dest= 'nucleotide_range', type=int, default = 3500, metavar='<int>', help='Max. nucleotide distance to be considered synthenic genes. Default: 3500')    
    csb.add_argument('-insertions', dest='insertions', type=int,default = 1, metavar='<int>', help='Max. insertions in a csb. Default: 1')
    csb.add_argument('-occurence', dest='occurence', type=int,default = 2, metavar='<int>', help='Min. number occurences to be recognized as csb. Default: 2')
    csb.add_argument('-min_csb_size', dest='min_csb_size', type=int,default = 4, metavar='<int>', help='Min. number of genes in a csb before recognized. Default: 4')
    csb.add_argument('-jaccard', dest='jaccard', type=float,default = 0.2, metavar='<float>', help='Acceptable dissimilarity in jaccard clustering [0-1]. Default: 0.2') #. 0.2 means that 80 percent have to be the same genes
    
    #Flow regulators
    flow = parser.add_argument_group("Work step regulation")
    flow.add_argument('-redo_csb', dest= 'redo_csb', action='store_true', help='Redo the collinear synthenic block prediction *. Default:False')
    flow.add_argument('-redo_csb_naming', dest= 'redo_csb_naming', action='store_true', help='Redo pattern matching for collinear synthenic blocks *. Default: False')
    flow.add_argument('-redo_search', dest= 'redo_search', action='store_true', help='Do not exclude genomes already stored in the database for the hmmsearch *. Default: False')
    flow.add_argument('-redo_taxonomy', dest = 'redo_taxonomy', action='store_true', help='Redo the taxonomy assignment*. Default: False')
    
    
    
    ########################################### Fetch & Process operators ########################
    #Limiter for genomes to account to
    limiter = parser.add_argument_group("Limit output to genomes with conditions *")
    limiter.add_argument('-dll',dest='dataset_limit_lineage', type=str, default = None, metavar='<string>', choices = ['Superkingdom','Phylum','Class','Ordnung','Family','Genus','Species'],\
                            help='Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]')
    limiter.add_argument('-dlt',dest='dataset_limit_taxon', type=str, default = None, metavar='<string>', help='Taxonomic name e.g. Proteobacteria. Requires -dll')
    limiter.add_argument('-dlp',dest='dataset_limit_proteins', type=str, default = 0, metavar='<list>', help='Limit fetch to genomes with <protein>')
    limiter.add_argument('-dlk',dest='dataset_limit_keywords', type=str, default = 0, metavar='<list>', help='Limit fetch to genomes with <keyword>')

    limiter.add_argument('-dtd',dest='dataset_divide_sign', default = '.', type=str, metavar='<string>', help='Separator for taxonomy information. The characters ";" ":" and "," cause strange behavior') #TODimplement this in fetch_bulk and all other dataset routines
    
    #Fetch operator
    operators = parser.add_argument_group("Output sequences with these conditions *")
    operators.add_argument('-fl',dest='dataset_limit_lineage', type=str, metavar='<string>', choices = ['Superkingdom','Phylum','Class','Ordnung','Family','Genus','Species'],\
                            help='Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]')
    operators.add_argument('-ft',dest='dataset_limit_taxon', type=str, metavar='<string>', help='Taxonomic name e.g. Proteobacteria. Requires -fl')
    operators.add_argument('-fg', nargs='+', dest='fetch_genomes', type=str, default=[],metavar='<list>', help='Select only genomes with these identifiers (whitespace separated)')
    operators.add_argument('-fd', nargs='+', dest='fetch_proteins', type=str, default=[],metavar='<list>', help='Select only proteins with these domains (whitespace separated)')
    operators.add_argument('-fc', nargs='+', dest='fetch_csbs', type=str, default=[],metavar='<list>', help='Select only csb containing the given proteins (whitespace separated)')
    operators.add_argument('-fk', nargs='+', dest='fetch_keywords', type=str, default=[],metavar='<list>', help='Select only proteins in gene cluster with this keyword (whitespace separated)')
    operators.add_argument('-kc', dest='keywords_connector', type=str, default = 'OR', choices = ['AND','OR'], help='Select cluster with keywords connected by AND or OR')
  
    #Alignment and sequence file processing filter_fasta
    process = parser.add_argument_group("Alignment and sequence file processing")
    process.add_argument('-merge_fasta',dest='merge_fasta', type=myUtil.dir_path, metavar='<directory>', help='Merges two or more sequence files with extention .faa without doublicates')
    process.add_argument('-filter_fasta', dest='filter_fasta', type=str, metavar='<file> <int> <int>', help='Filter fasta file by length upper and lower limit')
    process.add_argument('-concat_alignment',dest='concat_alignment', type=myUtil.dir_path, metavar='<directory>', help='Concatenates alignment files with extention .fasta_aln')
    process.add_argument('-add_taxonomy',dest='add_taxonomy', type=str, metavar='<file> or <directory>', help='Adds taxonomy to alignment files in <dir>, requires -db with taxonomy')
    process.add_argument('-add_genomic_context',dest='add_genomic_context', type=myUtil.file_path, metavar='<file>', help='Adds genomic context to sequences from fasta file, requires -db with taxonomy')
    process.add_argument('-create_type_range_dataset',dest='create_type_range_dataset', type=myUtil.file_path, metavar='<file>', help='Create protein type range dataset from sequences fasta file, requires -db with taxonomy')
    process.add_argument('-create_gene_cluster_dataset',dest='create_gene_cluster_dataset', type=myUtil.file_path, metavar='<file>', help='Create gene cluster dataset from sequences fasta file, requires -db with taxonomy')
    process.add_argument('-aln_gaps', dest='gaps',action='store_true',help='When concating alignments add gaps for missing sequences')

    
    
    #### Parse the arguments

    if len(sys.argv) == 1: # Check if no arguments were provided and exit
        parser.print_help()
        sys.exit(1)
        
        
    options = HMSSS()
    parser.parse_args(namespace=options)
    
    
    # Get default values dynamically from the parser
    default_values = {action.dest: action.default for action in parser._actions if action.dest != 'help'}
    # Filter out default values only for the limiters fetch operators and process operators
    relevant_groups = ['dataset_limit_lineage', 'dataset_limit_taxon', 'dataset_limit_proteins', 'dataset_limit_keywords', 'dataset_limit_min_cluster_completeness', 'dataset_divide_sign',
                       'fetch_proteins', 'fetch_genomes', 'fetch_csbs', 'fetch_keywords', 'keywords_connector', 'stat_genomes', 'stat_csb',
                       'merge_fasta', 'filter_fasta',
                       'concat_alignment', 'add_taxonomy', 'add_genomic_context', 'create_type_range_dataset', 'create_gene_cluster_dataset', 'gaps']
    default_values = {key: value for key, value in default_values.items() if key in relevant_groups}
    process_args = Project.any_process_args_provided(options, default_values) #returns bool
    
    
    
    
    if options.redo_taxonomy:
        if options.taxonomy_file and options.database_directory:
            options.stage = 5 # direct the stage to taxonomy addition
        else:
            sys.exit("Please use the -db argument to provide a valid database and -tax argument to provide a taxonomy table")
    if options.redo_csb_naming:
    	options.stage = 4
    if options.redo_csb:
        options.stage = 3
    
 #####Processes are at stage 100
    if process_args:
        options.stage = 100
        
        
        if options.fetch_proteins or options.fetch_keywords or options.fetch_csbs or options.fetch_genomes:
            options.fetch = True
            if not options.database_directory:
                sys.exit("Please use the -db argument to provide a valid database")

            if options.dataset_limit_proteins or options.dataset_limit_keywords or options.dataset_limit_lineage or options.dataset_limit_taxon:
                options.limiter = True
        
        if options.filter_fasta or options.concat_alignment or options.merge_fasta:
            options.process = True
        
        if options.add_taxonomy or options.add_genomic_context or options.create_type_range_dataset or options.create_gene_cluster_dataset:
            options.process = True
            if not options.database_directory:
                    sys.exit("Please use the -db argument to provide a valid database")


    
    return options
    

#############
####   Subroutines for preparation of the result folder
#############
     
    
  
            
def initial_search(options):
    
    Queue.queue_files(options) #looks for .faa.gz and .gff.gz file pairs and queues them
    
    #Create database if not exist
    if os.path.isfile(options.database_directory):
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        Queue.compare_with_existing_database(options,genomeIDs)
        Database.drop_specified_indices(options.database_directory)
    else:
        print(f"Saving results to database: {options.database_directory}")
        Database.create_database(options.database_directory)
    
    #Start the search or bulk parse
    if options.glob_report:
        BulkSearch.initial_bulk_reports(options)
    
    elif options.cores <= 2:
        ParallelSearch.initial_search(options)
    
    elif options.cores >= 3:
        ParallelSearch.multi_search_process(options)


def csb_finder(options):
    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(options)
    Database.delete_keywords_from_csb(options.database_directory, options) #remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) #assigns the names of the keywords to the clusters






def redo_block_assignment(options):
    """
        This routine has to be modified to run more efficiently and print the results out to the gene_cluster_file for further processing
        Parallelization is possible with a detached writing process
    """
    
    Database.light_index_database(options.database_directory) #make it light: for proteinIDs and genomeIDs only, not cluster
    Csb_redo.redo_csb_reports(options)
    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(options)
    Database.delete_keywords_from_csb(options.database_directory, options) #remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) #assigns the names of the keywords to the clusters


def redo_csb_names(options):
    """
        09.11.22
        Args:
            options object
        Return:
            nothing
    """
    
    myUtil.print_header(f"\nRedo collinear syntenic block assignment")
    print("Reading pattern")
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
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
    myUtil.print_header(f"\nTaxonomy assignment")
      
    print("Writing taxonomy to file --", end="\r")
    Database.insert_taxonomy_data(options.database_directory, options.taxonomy_file)
    print("Writing taxonomy to file -- ok")

def output_operator(options):
    """
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
    directory = os.path.dirname(options.database_directory)+"/"+str(date)+"_dataset/"
    os.mkdir(directory) #save results in the same folder as the database

    #primary output routine for fasta files
    protein_dict, cluster_dict, taxon_dict = Output.print_fasta_and_hit_table(directory,options)
    
    #dataset generation
    Datasets.main_binary_dataset(options, directory, taxon_dict, protein_dict, cluster_dict)
    return

    
def output_statistics(options):    
    Database.fetch_genome_statistic(options.database_directory)
    
    
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
    
    #Merge fasta files
    if options.merge_fasta:
        Processing.merge_fasta(options)
        
    #Concat alignment files    
    if options.concat_alignment:
        Processing.concat_alignments(options)
        

#Filter fasta files by length    
    if options.filter_fasta:
        Processing.filter_length_fasta(options.filter_fasta[0],options.filter_fasta[1],options.filter_fasta[2])
    
#Add taxonomy information    
    if options.add_taxonomy:
        Processing.taxonomy_comprehension(options)

    #Get a textfile with genomic context based on provided sequence fasta file    
    if options.add_genomic_context:
        Processing.add_genomic_context(options.database_directory,options.add_genomic_context)

    #Get a iTol dataset file with gene cluster dataset    
    if options.create_gene_cluster_dataset:
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
    
    Queue.prepare_HMMlib(__location__)

    options = parse_arguments(args)
    if options.stat_keywords:
        Output.print_file_content(options.patterns_file)
        sys.exit()

            
    #print(f"Start at stage {options.stage}")    
    myUtil.print_header("\nInitilizing result file directory")
    
    Project.prepare_result_space(options,options.name)
    
    
    if options.stat_csb:
        Output.print_file_content(options.csb_output_file)
        sys.exit()
        
    if options.stage <= 1 and not options.glob_report:
        #ignored if bulk is used because nobody should want to translate a glob via prodigal
        myUtil.print_header("\nProkaryotic gene recognition and translation via prodigal")
        Translation.parallel_translation(options.fasta_file_directory, options.cores)      #fine
        Translation.parallel_transcription(options.fasta_file_directory, options.cores)    #fine
    
    
    if options.stage <= 2 or options.redo_search:
        myUtil.print_header("\nSearching for homologoues sequences")
        initial_search(options)
        options.stage = 2
        
    if options.redo_csb or options.redo_search:
        myUtil.print_header(f"\nRedo collinear syntenic block assignment")
        redo_block_assignment(options)
        options.stage = 5

    if options.stage <= 3:
        myUtil.print_header("\nSearching for collinear syntenic blocks")
        csb_finder(options)
   
    if options.redo_csb_naming:
        redo_csb_names(options)

    if options.stage <= 5:
        collect_taxonomy_information(options)
            


    
    
    
# These routines modify the existing data and output    
        
    if options.fetch:
        #14
        Database.index_database(options.database_directory)
        output_operator(options)
        
    if options.stat_genomes:
        #9
        output_statistics(options)
    

    
    if options.process:
        # file/alignment concat utilities
        process_operator(options)
    

    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args) #calls the main method of __main__


#!/usr/bin/python
import os
import sys
import argparse
from datetime import datetime

from . import Csb_finder
from . import Csb_cluster
from . import Database
from . import Datasets
from . import Output
from . import myUtil
from . import ParseReports
from . import Processing
from . import Project
from . import Search
from . import Translation
from . import Queue


# - sollen unterschiedlicher feld informationen trennen
# _ sollen namentrennungen sein, bzw indices
# get location of script or executable
# for the output module report
# TODO print a database description file from database (otherwise it is too complex for a short task) should include the
# TODO make the synteny completion, transitions and removal optional
logger = myUtil.logger

if getattr(sys, "frozen", False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )
# print(__location__)


class HMSSS:
    """
    Objects of this class include the parameters for the program given by the user or by default
    """

    def __init__(self):
        # Parameters
        self.execute_location = __location__

        # queue functional dictionaries
        self.finished_genomes = {}
        self.queued_genomes = {}
        self.faa_files = {}
        self.gff_files = {}

        # csb prediction dereplication
        self.redundant = 0
        self.non_redundant = 0
        self.redundancy_hash = dict()

        # Limiter
        self.limiter = False

        # Fetching data
        self.fetch = False

        # Alignment and fasta file processing
        self.process = False

        self.project_name = "project"
        self.index_db = False

        self.csb_name_prefix = (
            "csb-"  # prefix of clusterIDs determined by csb finder algorithm
        )
        self.csb_name_suffix = (
            "_"  # suffix of clusterIDs determined by csb finder algorithm
        )

        self.genomeID_divider = "___"  # dividing sign between genomeID and proteinID, first part will be taken as genomeID


def parse_arguments(arguments: list):
    """
    Argument parser organizing the different groups of arguments the user can give.
    This includes Search and directory and workflow operators. Also operators for the
    output of data in sequence fasta format, and iTol dataset format
    Also this includes some useful operators for processing the default output files
    regarding taxonomy and sequence sorting/merging and concatenation
    """

    # Check if --help-all is in CLI args
    show_all = "--help-all" in arguments
    clean_args = [a for a in arguments if a != "--help-all"]

    def add_all_groups(parser, show_advanced: bool):
        """
        Add all argument groups to the parser.
        """
        # Basic inputs
        inputdef = parser.add_argument_group("Input definition")
        inputdef.add_argument(
            "-f",
            dest="fasta_file_directory",
            type=myUtil.dir_path,
            default=None,
            metavar="<directory>",
            help="Directory to be searched",
        )
        inputdef.add_argument(
            "-t",
            dest="score_threshold_file",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to tab separated threshold file with optimized, trusted and noise cutoff"
            if show_advanced
            else argparse.SUPPRESS,
        )
        inputdef.add_argument(
            "-l",
            dest="library",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to a custom HMM library for the hmmsearch"
            if show_advanced
            else argparse.SUPPRESS,
        )
        inputdef.add_argument(
            "-r",
            dest="result_files_directory",
            type=myUtil.dir_path,
            metavar="<directory>",
            default=__location__ + "/results",
            help="Directory for the result files",
        )
        inputdef.add_argument(
            "-db",
            dest="database_directory",
            type=myUtil.file_path,
            metavar="<filepath>",
            help="Filepath to existing sqlite database",
        )
        inputdef.add_argument(
            "-c",
            dest="cores",
            type=int,
            default=4,
            metavar="<int>",
            help="Allocated CPU cores" if show_advanced else argparse.SUPPRESS,
        )
        inputdef.add_argument(
            "-glob_report",
            dest="glob_report",
            type=myUtil.file_path,
            metavar="<filepath>",
            help="Filepath to glob hmmreport. Each report with one HMM queried against the concatenated genomes."
            if show_advanced
            else argparse.SUPPRESS,
        )
        parser.add_argument(
            "-v",
            "--verbose",
            type=int,
            default=1,
            choices=[0, 1, 2],
            help="Set logging level: 0=WARNING, 1=INFO, 2=DEBUG",
        )
        parser.add_argument(
            "--help-all",
            action="store_true",
            help="Show all available options, including advanced parameters, and exit.",
        )

        # Search options
        parameters = parser.add_argument_group("Search parameters")
        parameters.add_argument(
            "-cut_type",
            dest="threshold_type",
            type=int,
            default=1,
            metavar="<int>",
            choices=[1, 2, 3],
            help="Choice of cutoff: 1 optimized; 2 trusted; 3 noise"
            if show_advanced
            else argparse.SUPPRESS,
        )
        parameters.add_argument(
            "-cut_score",
            dest="thrs_score",
            type=int,
            default=50,
            metavar="<int>",
            help="Global minimal score cutoff" if show_advanced else argparse.SUPPRESS,
        )
        parameters.add_argument(
            "-taxonomy_info",
            dest="taxonomy_file",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to tab separated taxonomy file",
        )
        parameters.add_argument(
            "-refseq_ident",
            dest="refseq_identity",
            type=int,
            default=90,
            metavar="<int>",
            help="Minimal percent identity to reference sequence set"
            if show_advanced
            else argparse.SUPPRESS,
        )
        parameters.add_argument(
            "-n",
            dest="name",
            type=str,
            default="project",
            metavar="<string>",
            help="Name new project" if show_advanced else argparse.SUPPRESS,
        )
        parameters.add_argument(
            "-s",
            dest="stage",
            type=int,
            default=0,
            choices=[0, 1, 2, 3, 4, 5],
            help="Start at step" if show_advanced else argparse.SUPPRESS,
        )
        parameters.add_argument(
            "-x",
            dest="exit",
            type=int,
            default=10,
            choices=[0, 1, 2, 3, 4, 5],
            help="Exit at step" if show_advanced else argparse.SUPPRESS,
        )

        # Search library resources
        resources = parser.add_argument_group("Search library resources")

        resources.add_argument(
            "-hmms",
            nargs="+",
            dest="HMM_sets",
            type=str,
            default=[],
            metavar="<list>",
            help="Limit to HMM sets (whitespace separated)"
            if show_advanced
            else argparse.SUPPRESS,
        )
        resources.add_argument(
            "-clean",
            dest="clean_reports",
            action="store_true",
            help="Clean up any pre-existing HMMreport file"
            if show_advanced
            else argparse.SUPPRESS,
        )
        resources.add_argument(
            "-no_reports",
            dest="individual_reports",
            action="store_false",
            help="Do not write individual files per genome"
            if show_advanced
            else argparse.SUPPRESS,
        )
        resources.add_argument(
            "-max_seq_per_genome",
            dest="max_seqs_per_genome",
            type=int,
            default=4,
            help="Max. number of sequences per protein per genome for cross check via Diamond"
            if show_advanced
            else argparse.SUPPRESS,
        )
        resources.add_argument(
            "-no_cross_check",
            dest="bool_cross_check",
            action="store_false",
            help="No cross check with reference sequences via Diamond"
            if show_advanced
            else argparse.SUPPRESS,
        )
        resources.add_argument(
            "-optimized_cutoff_cross_check",
            dest="optimized_cutoff_cross_check",
            action="store_true",
            help="Use optimized cutoff instead of cross check with Diamond"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Synteny options
        synteny = parser.add_argument_group("Synteny options")
        synteny.add_argument(
            "-p",
            dest="patterns_file",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to patterns file" if show_advanced else argparse.SUPPRESS,
        )
        synteny.add_argument(
            "-cooccurrence",
            dest="cooccurrence_file",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to co-occurrence file"
            if show_advanced
            else argparse.SUPPRESS,
        )
        synteny.add_argument(
            "-exclude_singletons",
            dest="exclusion_singletons",
            type=myUtil.file_path,
            default=None,
            metavar="<filepath>",
            help="Filepath to tab separated file for singletons that are excluded"
            if show_advanced
            else argparse.SUPPRESS,
        )
        synteny.add_argument(
            "-mc",
            dest="min_completeness",
            type=float,
            default=0.5,
            metavar="<float>",
            help="Minimal fraction of predefined csb to be recognized"
            if show_advanced
            else argparse.SUPPRESS,
        )
        synteny.add_argument(
            "-chunks",
            dest="glob_chunks",
            type=int,
            default=5000,
            metavar="<int>",
            help="Chunk size for parsing results from glob before entering into database"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Information on resources
        information = parser.add_argument_group("Information on resources")
        information.add_argument(
            "-stat_keywords",
            action="store_true",
            help="Print patterns for keyword naming"
            if show_advanced
            else argparse.SUPPRESS,
        )
        information.add_argument(
            "-stat_csb",
            action="store_true",
            help="Print automatically found csbs"
            if show_advanced
            else argparse.SUPPRESS,
        )
        information.add_argument(
            "-stat_genomes",
            action="store_true",
            help="Print taxonomy information from database"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Collinear syntenic block prediction
        csb = parser.add_argument_group("Collinear syntenic block prediction")
        csb.add_argument(
            "-nt",
            dest="nucleotide_range",
            type=int,
            default=3500,
            metavar="<int>",
            help="Max. nucleotide distance to be considered synthenic genes"
            if show_advanced
            else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-insertions",
            dest="insertions",
            type=int,
            default=1,
            metavar="<int>",
            help="Max. insertions in a csb" if show_advanced else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-occurence",
            dest="occurence",
            type=int,
            default=2,
            metavar="<int>",
            help="Min. number occurences to be recognized as csb"
            if show_advanced
            else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-min_csb_size",
            dest="min_csb_size",
            type=int,
            default=4,
            metavar="<int>",
            help="Min. number of genes in a csb before recognized"
            if show_advanced
            else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-max_csb_size",
            dest="max_csb_size",
            type=int,
            default=50,
            metavar="<int>",
            help="Max. number of genes in a csb before recognized"
            if show_advanced
            else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-max_gene_repeats",
            dest="max_domain_repeats",
            type=int,
            default=4,
            metavar="<int>",
            help="Maximum number of repeated genes in a csb."
            if show_advanced
            else argparse.SUPPRESS,
        )
        csb.add_argument(
            "-jaccard",
            dest="jaccard",
            type=float,
            default=0.0,
            metavar="<float>",
            help="Acceptable dissimilarity in jaccard clustering [0.0-1.0]"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Work step regulation
        flow = parser.add_argument_group("Work step regulation")
        flow.add_argument(
            "-redo_taxonomy",
            dest="redo_taxonomy",
            action="store_true",
            help="Redo the taxonomy assignment*. Default: False"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Limiter for genomes to account to
        limiter = parser.add_argument_group("Limit output to genomes with conditions *")
        limiter.add_argument(
            "-dll",
            dest="dataset_limit_lineage",
            type=str,
            default=None,
            metavar="<string>",
            choices=[
                "Superkingdom",
                "Phylum",
                "Class",
                "Ordnung",
                "Family",
                "Genus",
                "Species",
            ],
            help="Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]"
            if show_advanced
            else argparse.SUPPRESS,
        )
        limiter.add_argument(
            "-dlt",
            dest="dataset_limit_taxon",
            type=str,
            default=None,
            metavar="<string>",
            help="Taxonomic name e.g. Proteobacteria. Requires -dll"
            if show_advanced
            else argparse.SUPPRESS,
        )
        limiter.add_argument(
            "-dlp",
            dest="dataset_limit_proteins",
            type=str,
            default=0,
            metavar="<list>",
            help="Limit fetch to genomes with <protein>"
            if show_advanced
            else argparse.SUPPRESS,
        )
        limiter.add_argument(
            "-dlk",
            dest="dataset_limit_keywords",
            type=str,
            default=0,
            metavar="<list>",
            help="Limit fetch to genomes with <keyword>"
            if show_advanced
            else argparse.SUPPRESS,
        )
        limiter.add_argument(
            "-dtd",
            dest="dataset_divide_sign",
            default=".",
            type=str,
            metavar="<string>",
            help='Separator for taxonomy information. The characters ";" ":" and "," cause strange behavior'
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Output sequences with these conditions
        operators = parser.add_argument_group(
            "Output sequences with these conditions *"
        )
        operators.add_argument(
            "-fl",
            dest="dataset_limit_lineage",
            type=str,
            metavar="<string>",
            choices=[
                "Superkingdom",
                "Phylum",
                "Class",
                "Ordnung",
                "Family",
                "Genus",
                "Species",
            ],
            help="Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]",
        )
        operators.add_argument(
            "-ft",
            dest="dataset_limit_taxon",
            type=str,
            metavar="<string>",
            help="Taxonomic name e.g. Proteobacteria. Requires -fl"
            if show_advanced
            else argparse.SUPPRESS,
        )
        operators.add_argument(
            "-fg",
            nargs="+",
            dest="fetch_genomes",
            type=str,
            default=[],
            metavar="<list>",
            help="Select only genomes with these identifiers (whitespace separated)"
            if show_advanced
            else argparse.SUPPRESS,
        )
        operators.add_argument(
            "-fd",
            nargs="+",
            dest="fetch_proteins",
            type=str,
            default=[],
            metavar="<list>",
            help="Select only proteins with these domains (whitespace separated)",
        )
        operators.add_argument(
            "-fc",
            nargs="+",
            dest="fetch_csbs",
            type=str,
            default=[],
            metavar="<list>",
            help="Select only csb encoding the given proteins (whitespace separated)",
        )
        operators.add_argument(
            "-fk",
            nargs="+",
            dest="fetch_keywords",
            type=str,
            default=[],
            metavar="<list>",
            help="Select only proteins in gene cluster with this keyword (whitespace separated)"
            if show_advanced
            else argparse.SUPPRESS,
        )
        operators.add_argument(
            "-kc",
            dest="keywords_connector",
            type=str,
            default="OR",
            choices=["AND", "OR"],
            help="Select cluster with keywords connected by AND or OR"
            if show_advanced
            else argparse.SUPPRESS,
        )

        # Alignment and sequence file processing
        process = parser.add_argument_group("Alignment and sequence file processing")
        process.add_argument(
            "-merge_fasta",
            dest="merge_fasta",
            type=myUtil.dir_path,
            metavar="<directory>",
            help="Merges two or more sequence files with extention .faa without doublicates"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-filter_fasta",
            dest="filter_fasta",
            type=str,
            metavar="<file> <int> <int>",
            help="Filter fasta file by length upper and lower limit"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-concat_alignment",
            dest="concat_alignment",
            type=myUtil.dir_path,
            metavar="<directory>",
            help="Concatenates alignment files with extention .fasta_aln"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-add_taxonomy_to_alignment",
            dest="add_taxonomy",
            type=str,
            metavar="<file> or <directory>",
            help="Adds taxonomy to alignment files in <dir>, requires -db with taxonomy"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-add_genomic_context",
            dest="add_genomic_context",
            type=myUtil.file_path,
            metavar="<file>",
            help="Adds genomic context to sequences from fasta file, requires -db with taxonomy"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-create_type_range_dataset",
            dest="create_type_range_dataset",
            type=myUtil.file_path,
            metavar="<file>",
            help="Create protein type range dataset from sequences fasta file, requires -db with taxonomy"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-create_gene_cluster_dataset",
            dest="create_gene_cluster_dataset",
            type=myUtil.file_path,
            metavar="<file>",
            help="Create gene cluster dataset from sequences fasta file, requires -db with taxonomy"
            if show_advanced
            else argparse.SUPPRESS,
        )
        process.add_argument(
            "-aln_gaps",
            dest="gaps",
            action="store_true",
            help="When concating alignments add gaps for missing sequences"
            if show_advanced
            else argparse.SUPPRESS,
        )

    # ---- Build parser (show_advanced = True <-> --help-all, else False) ----
    formatter = lambda prog: argparse.HelpFormatter(
        prog, max_help_position=96, width=300
    )

    parser = argparse.ArgumentParser(
        description="HMSS2: Sulfur metabolism annotation",
        epilog="Please cite: Tanabe TS, Dahl C. HMSS2: An advanced tool for the analysis of sulphur metabolism, including organosulphur compound transformation, in genome and metagenome assemblies. Mol Ecol Resour. 2023;23(8):1930-1945. doi:10.1111/1755-0998.13848",  # Formatter makes this a one-liner
        usage="HMSSS.py -f <genomes_dir> [options]\n       HMSSS.py --help-all",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_all_groups(parser, show_advanced=show_all)

    # --help-all triggers full help and exits
    if show_all:
        parser.print_help()
        sys.exit(0)

    # Print help if no arguments were provided
    if len(arguments) == 0:
        parser.print_help()
        sys.exit("Please provide arguments. For help use -h or --help_all")

    options = HMSSS()
    parser.parse_args(namespace=options)

    # Set the locations of src directory and major src files
    options.location = __location__

    options.reference_seq_dir = __location__ + "/src/RefSeqs"

    if options.score_threshold_file is None:
        options.score_threshold_file = __location__ + "/src/Thresholds"

    if options.library is None:
        options.library = __location__ + "/src/HMMlib"

    if options.patterns_file is None:
        options.patterns_file = __location__ + "/src/Patterns"

    if options.cooccurrence_file is None:
        options.cooccurrence_file = __location__ + "/src/Cooccurrence"

    if options.exclusion_singletons is None:
        options.exclusion_singletons = __location__ + "/src/Exclusion_singletons"

    # Check if results dir is default location
    if options.result_files_directory == __location__ + "/results":
        # default location is never an existing project
        options.new_project = True
    else:
        options.new_project = False

    # Get default values dynamically from the parser
    default_values = {
        action.dest: action.default
        for action in parser._actions
        if action.dest != "help"
    }
    # Filter out default values only for the limiters fetch operators and process operators
    relevant_groups = [
        "dataset_limit_lineage",
        "dataset_limit_taxon",
        "dataset_limit_proteins",
        "dataset_limit_keywords",
        "dataset_limit_min_cluster_completeness",
        "dataset_divide_sign",
        "fetch_proteins",
        "fetch_genomes",
        "fetch_csbs",
        "fetch_keywords",
        "keywords_connector",
        "stat_genomes",
        "stat_csb",
        "merge_fasta",
        "filter_fasta",
        "concat_alignment",
        "add_taxonomy",
        "add_genomic_context",
        "create_type_range_dataset",
        "create_gene_cluster_dataset",
        "gaps",
    ]
    default_values = {
        key: value for key, value in default_values.items() if key in relevant_groups
    }
    process_args = Project.any_process_args_provided(
        options, default_values
    )  # returns bool

    if (
        options.taxonomy_file
        and options.database_directory
        and not options.fasta_file_directory
    ):
        options.stage = 100  # direct the stage to taxonomy addition
        logger.info("Redo taxonomy assignment with database and taxonomy tsv")
    #####Processes are at stage 100
    if process_args:
        options.stage = 101

        if (
            options.fetch_proteins
            or options.fetch_keywords
            or options.fetch_csbs
            or options.fetch_genomes
        ):
            options.fetch = True
            if not options.database_directory:
                sys.exit("Please use the -db argument to provide a valid database")

            if (
                options.dataset_limit_proteins
                or options.dataset_limit_keywords
                or options.dataset_limit_lineage
                or options.dataset_limit_taxon
            ):
                options.limiter = True

        if options.filter_fasta or options.concat_alignment or options.merge_fasta:
            options.process = True

        if (
            options.add_taxonomy
            or options.add_genomic_context
            or options.create_type_range_dataset
            or options.create_gene_cluster_dataset
        ):
            options.process = True
            if not options.database_directory:
                sys.exit("Please use the -db argument to provide a valid database")

    return options


#############
####   Subroutines for preparation of the result folder
#############
def ressource_preparation(options: HMSSS) -> None:
    """
    Prepares the cutoffs, patterns, and HMMs
    For HMMs and cutoffs these are concatenated
    For patterns, whitespaces are replaced with tabs and
    numerical pattern names are given
    """
    if options.HMM_sets:
        # specific HMM sets only
        allowed_words = (
            options.HMM_sets
            if isinstance(options.HMM_sets, list)
            else options.HMM_sets.split()
        )
        Queue.concatenate_selected_hmms(
            __location__ + "/src", allowed_words, "", ".hmm", options.library
        )

    if not os.path.isfile(options.library):
        Queue.concatenate_files_shell(
            __location__ + "/src", "grp", ".hmm", options.library
        )

    if not os.path.isfile(options.score_threshold_file):
        Queue.concatenate_files_shell(
            __location__ + "/src", "cutoffs", ".txt", options.score_threshold_file
        )

    if not os.path.isfile(options.cooccurrence_file):
        Queue.concatenate_files_shell(
            __location__ + "/src", "cooccurrence", ".txt", options.cooccurrence_file
        )
        Queue.format_pattern_files_inplace(
            options.cooccurrence_file, "cpb-", options.csb_name_suffix
        )  #  cob- is the co-occurring protein block

    if not os.path.isfile(options.patterns_file):
        Queue.concatenate_files_shell(
            __location__ + "/src", "patterns", ".txt", options.patterns_file
        )
        Queue.format_pattern_files_inplace(
            options.patterns_file, "dsb-", options.csb_name_suffix
        )  # dsb- is the defined collinear syntenic block

    return


def fasta_preparation(options: HMSSS) -> None:
    """
    Prepares/queues genome fasta files for analysis, supports concatenation and deconcatenation.

    Args:
        options (Options): Main pipeline options

    Input Example:
        options.glob_faa: str = "all.faa"
        options.fasta_file_directory: str = "./genomes"

    Output:
        - Queues or deconcatenates files for search
        - Generates glob files if needed
    """

    # Unpacks and translates fasta
    Translation.parallel_translation(options.fasta_file_directory, options.cores)
    Translation.parallel_transcription(options.fasta_file_directory, options.cores)

    return


def initial_search(options: object) -> None:
    """Performs the initial HMM search workflow and organizes search hits.

    This function builds a protein hit database, performs HMM searches (unless results already exist),
    filters hits by thresholds or crosschecks against reference sequences, and writes summary reports.

    Args:
        options (object): Configuration object containing:
            - database_directory (str): Path to the output database file.
            - glob_report (str): Path to the global HMM search result file.
            - cores (int): Number of CPU cores to use.
            - bool_cross_check (bool): Whether to use reference crosschecking.
            - Cross_check_directory (str): Folder to hold intermediate .faa and .tsv files.
            - glob_trusted_hitreport (str): Output path for trusted hit reports.
            - result_files_directory (str): Directory for result files.
            - glob_intermediate_hitreport (str): Output path for intermediate summary.

    Returns:
        None

    Example:
        >>> initial_search(options)
    """

    # Writes a database with the protein hits and their gene clusters for later use.

    if not os.path.isfile(options.database_directory):
        Database.create_database(options.database_directory)

    # Create a glob report or use the provided glob report
    if options.glob_report and os.path.isfile(options.glob_report):
        logger.info(
            "Glob report %s was provided, hmmsearch will be omitted",
            options.glob_report,
        )
    else:
        logger.info("Running unified hmmsearch")
        genome_to_hmmreport_dict = Search.unified_search(
            options, int(options.cores / 2)
        )
        logger.info("Concatenating hmmsearch results")
        options.glob_report = Search.concatenate_hmmreports_cat(
            genome_to_hmmreport_dict, options.glob_report
        )

    # Now the glob_report exists with names for proteins with genomeID___proteinID
    logger.info("Filtering hits into trusted, noise, and intermediate hit category")
    Search.filter_trusted_and_noise_hits(
        options, options.glob_report, options.cores
    )  # Sort by TC, above NC and below NC into Cross check directory


def reference_sequence_check(options: object) -> None:
    ### Test hits NC < score < TC with reference sequences
    if options.bool_cross_check and not os.path.isfile(options.glob_trusted_hitreport):
        # Create the .faa fasta files from the intermediate hit lists
        # Promote hits related to ref seqs to trusted_hit list
        logger.info(
            "Initializing testing hits NC < score < TC with reference sequences"
        )
        Search.generate_faa_per_hitfile_parallel(
            options, options.Cross_check_directory, options.cores
        )
        Search.promote_crosschecked_hits(options.Cross_check_directory, options.cores)

        # Perform inclusion of hits above chosen cutoff for the hits without refseqs
        logger.info("Testing hits NC < score < TC with reference sequences")
        refseq_unavailable_list = Search.cross_check_candidates_with_reference_seqs(
            options
        )  # Cross check the with diamond against reference sequences
        Search.promote_by_cutoff(
            options,
            options.Cross_check_directory,
            options.cores,
            refseq_unavailable_list,
        )

    elif options.optimized_cutoff_cross_check:
        # Perform inclusion of hits above chosen cutoff in options
        logger.info(
            "Using optimized cutoff instead of trusted cutoff or cross referencing"
        )
        Search.promote_by_cutoff(
            options, options.Cross_check_directory, options.cores, "all"
        )

    ### Summarize the trusted hits
    logger.info("Summarizing trusted and intermediate hit reports")
    options.glob_trusted_hitreport = Search.summarize_trusted_hits(
        options.result_files_directory,
        options.Cross_check_directory,
        "global_trusted_hits_summary.hmmreport",
        ".trusted_hits",
    )
    options.glob_intermediate_hitreport = Search.summarize_trusted_hits(
        options.result_files_directory,
        options.Cross_check_directory,
        "global_intermediate_hits_summary.hmmreport",
        "intermediate_hits",
    )
    return


def parse_reports_to_database(options: object) -> None:
    logger.info("Parsing summary report into database")
    ParseReports.main_parse_summary_hmmreport(options)

    return


def csb_finder(options: HMSSS) -> None:
    """
    Runs the CSB finder algorithm and updates the database with syntenic block clusters.

    This function uses clustering logic to predict CSB (Conserved Syntenic Blocks),
    updates the keyword database with cluster results, and removes outdated entries.

    Args:
        options (object): Configuration object containing:
            - database_directory (str): Path to the SQLite database.
            - csb_name_prefix (str): Prefix for identifying CSB keyword entries.
            - csb_name_suffix (str): Suffix for identifying CSB keyword entries.

    Input Example:
        options.database_directory = "./results/my_db.sqlite"

    Output:
        Updates the SQLite database with predicted CSB cluster keyword assignments.

    Example:
        >>> csb_finder(options)
    """
    logger.info("Running collinear syntenic block pattern prediction")
    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(
        options, 0.0
    )  # 0.0 does merge csb but create the csb cluster dict

    Database.index_database(options.database_directory)
    Database.delete_keywords_from_csb(
        options.database_directory, options
    )  # remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(
        options.database_directory, csb_gene_cluster_dict
    )  # assigns the names of the keywords to the clusters

    return


def collect_taxonomy_information(options: object) -> None:
    """Collects taxonomy information and stores it in the database.

    This function checks for the presence of a tab separated taxonomy file
    and inserts its contents into the database if available.

    Args:
        options (object): Configuration object containing:
            - taxonomy_file (str): Path to the taxonomy file.
            - database_directory (str): Path to the SQLite database.

    Input Example:
        options.taxonomy_file = "./data/taxonomy.tsv"
        options.database_directory = "./results/my_db.sqlite"

    Output:
        Inserts taxonomy data into the database if the file exists.

    Example:
        >>> collect_taxonomy_information(options)
    """
    if not options.taxonomy_file is None and os.path.isfile(options.taxonomy_file):
        logger.info("Writing taxonomy assignments to database")
        Database.insert_taxonomy_data(options.database_directory, options.taxonomy_file)
    else:
        logger.warning("Taxonomy file was not provided")


def output_operator(options: HMSSS) -> None:
    """Handles the final output of results including FASTA and metadata files.

    This function creates a results directory based on the current timestamp,
    exports FASTA and hit tables, and generates the dataset binary.

    Args:
        options (object): Configuration object containing:
            - database_directory (str): Path to the SQLite database.

    Output:
        Creates:
            - FASTA formatted file with protein sequences.
            - Metadata file with annotation and clustering results.

    Example:
        >>> output_operator(options)
    """

    # Set directory
    date_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # YYYY-MM-DD_HH-MM-SS
    directory = os.path.join(
        os.path.dirname(options.database_directory), f"{date_str}_dataset/"
    )
    os.mkdir(directory)  # save results in the same folder as the database
    logger.info("Created output directory: %s", directory)

    # primary output routine for fasta files
    Output.print_command_line_args(os.path.join(directory, "1_fetch_command.txt"))

    protein_dict, cluster_dict, taxon_dict = Output.fetch_fasta_and_hit_data(options)
    Output.print_fasta_and_hit_outputs(
        directory, protein_dict, cluster_dict, taxon_dict
    )
    # dataset generation
    Datasets.main_binary_dataset(
        options, directory, protein_dict, cluster_dict, taxon_dict
    )
    logger.info("Generated binary dataset")

    return


def output_statistics(options: HMSSS) -> None:
    Database.fetch_genome_statistic(options.database_directory)


def process_operator(options: HMSSS) -> None:
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

    # Merge fasta files
    if options.merge_fasta:
        Processing.merge_fasta(options)

    # Concat alignment files
    if options.concat_alignment:
        Processing.concat_alignments(options)

    # Filter fasta files by length
    if options.filter_fasta:
        Processing.filter_length_fasta(
            options.filter_fasta[0], options.filter_fasta[1], options.filter_fasta[2]
        )

    # Add taxonomy information
    if options.add_taxonomy:
        Processing.taxonomy_comprehension(options)

    # Get a textfile with genomic context based on provided sequence fasta file
    if options.add_genomic_context:
        Processing.add_genomic_context(
            options.database_directory, options.add_genomic_context
        )

    # Get a iTol dataset file with gene cluster dataset
    if options.create_gene_cluster_dataset:
        directory = os.path.dirname(options.create_gene_cluster_dataset)
        Datasets.iTol_domain_dataset(
            directory,
            options.database_directory,
            options.create_gene_cluster_dataset,
            options.dataset_divide_sign,
        )

    # Get a iTol dataset file with range data per protein type
    if options.create_type_range_dataset:
        if not options.database_directory:
            print(
                "WARNING: Missing database to assign taxonomy, please use -db argument"
            )
            return
        directory = os.path.dirname(options.create_type_range_dataset)
        Dataset.iTol_range_dataset(
            directory,
            options.database_directory,
            options.create_type_range_dataset,
            options.dataset_divide_sign,
        )


def main(args=None):
    myUtil.print_header("\nInitilizing result file directory")
    options = parse_arguments(args)

    # Initilize the logger
    # 1
    if options.stage < 100:
        # Prepare results directory and new project
        Project.prepare_result_space(options)
        log_file = os.path.join(options.result_files_directory, "execution_logfile.txt")
        myUtil.setup_logging(getattr(options, "verbose", 0), log_file)

        # Set up library, pattern and cutoff files
        ressource_preparation(options)

    if options.stage <= 1 and options.exit >= 1:
        # ignored if bulk is used because nobody should want to translate a glob via prodigal
        myUtil.print_header(
            "\nProkaryotic gene recognition and translation via prodigal"
        )
        fasta_preparation(options)

    if options.stage <= 2 and options.exit >= 2:
        # Queue the .faa/.gff file pairs
        Queue.queue_files(options)

        myUtil.print_header("\nSearching for homologoues sequences")
        initial_search(options)
        options.stage = 2

    if options.stage <= 3 and options.exit >= 3:
        myUtil.print_header("\nCross check with reference sequences")
        reference_sequence_check(options)
        options.stage = 3

    if options.stage <= 4 and options.exit >= 4:
        myUtil.print_header("\nParse trusted hits and recognized gene clusters")
        parse_reports_to_database(options)
        options.stage = 4

    if options.stage <= 5 and options.exit >= 5:
        myUtil.print_header("\nSearching for collinear syntenic blocks")
        csb_finder(options)
        options.stage = 5

    if options.stage <= 6 and options.exit >= 6:
        # Add taxonomy to existing database
        myUtil.print_header(f"\nAssigning taxonomy information")
        collect_taxonomy_information(options)

    ######### Output routines from main
    if options.stage > 99:
        log_file = os.path.join(options.result_files_directory, "execution_logfile.txt")
        myUtil.setup_logging(getattr(options, "verbose", 0), log_file)

    if options.stage == 100:
        # Add taxonomy to existing database
        myUtil.print_header(f"\nAssigning taxonomy information")
        collect_taxonomy_information(options)

    # These routines modify the existing data and output

    if options.fetch:
        # 14
        myUtil.print_header(f"\nOutput from database")
        Project.prepare_minimal_output_context(options)
        Database.index_database(options.database_directory)
        output_operator(options)

    if options.stat_genomes:
        # 9
        output_statistics(options)

    if options.process:
        # file/alignment concat utilities
        process_operator(options)

    if options.stat_keywords:
        Output.print_file_content(options.patterns_file)
        sys.exit()
    if options.stat_csb:  # TODO move to argument parser for direct execution
        Output.print_file_content(options.csb_output_file)
        sys.exit()


if __name__ == "__main__":
    args = sys.argv[1:]

    main(args)  # calls the main method of __main__

# src/hmsss/cli/parse.py
from __future__ import annotations

import argparse
import os
import sys
from typing import List

from hmsss.utils.paths import DATA_ROOT
from hmsss.core.options import Hmsss as Options
from hmsss.db import project as project


# ---------------------------------------------------------------------------
# Hilfsroutinen (ersetzen die bisher in myUtil verwendeten argparse-Validatoren)
# ---------------------------------------------------------------------------


def dir_path(p: str) -> str:
    """Verlangt ein existierendes Verzeichnis; gibt den absoluten Pfad zurück."""
    ap = os.path.abspath(p)
    if not os.path.isdir(ap):
        raise argparse.ArgumentTypeError(f"directory does not exist: {p}")
    return ap


def file_path(p: str) -> str:
    """Verlangt existierenden Pfad (Datei ODER Verzeichnis); gibt den absoluten Pfad zurück."""
    ap = os.path.abspath(p)
    if not os.path.exists(ap):
        raise argparse.ArgumentTypeError(f"path does not exist: {p}")
    return ap


def path_str(p: str) -> str:
    """Nur Normalisierung: absoluter Pfad; Existenz wird NICHT geprüft (für -r/-db/Output-Ziele)."""
    return os.path.abspath(p)


def _project_root_from_this_file() -> str:
    """Geht von src/hmsss/cli/parse.py drei Ebenen nach oben → HMSS2/"""
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.abspath(os.path.join(here, "..", "..", ".."))


def _data_dir() -> str:
    return str(DATA_ROOT)


def _list_from_csv_or_repeat(values: List[str]) -> List[str]:
    """Erlaubt -hmms A B C oder -hmms A,B,C."""
    out: List[str] = []
    for v in values or []:
        out.extend([p for p in v.split(",") if p])
    return out


# ---------------------------------------------------------------------------
# Hauptfunktion: parse_arguments
# ---------------------------------------------------------------------------


def parse_arguments(arguments: List[str]) -> Options:
    """
    Baut den Argumentparser, liest Args in ein Options(HMSSS)-Objekt ein,
    setzt Ressourcen-Defaults (HMSS2/data) und erkennt Modi/Stages.
    """
    show_all = "--help-all" in (arguments or [])
    clean_args = [a for a in (arguments or []) if a != "--help-all"]

    # ---- Argument-Gruppen genau wie in deiner aktuellen main ----
    parser = argparse.ArgumentParser(
        description="HMSS2: Sulfur metabolism annotation",
        epilog=(
            "Please cite: Tanabe TS, Dahl C. HMSS2: An advanced tool for the analysis of "
            "sulphur metabolism, including organosulphur compound transformation, in genome "
            "and metagenome assemblies. Mol Ecol Resour. 2023;23(8):1930-1945. doi:10.1111/1755-0998.13848"
        ),
        usage="hmsss -f <genomes_dir> [options]\n       hmsss --help-all",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input definition
    inputdef = parser.add_argument_group("Input definition")
    inputdef.add_argument(
        "-f",
        dest="fasta_file_directory",
        type=dir_path,
        default=None,
        metavar="<directory>",
        help="Directory to be searched",
    )
    inputdef.add_argument(
        "-t",
        dest="score_threshold_file",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Filepath to tab separated threshold file with optimized, trusted and noise cutoff"
        if show_all
        else argparse.SUPPRESS,
    )
    inputdef.add_argument(
        "-l",
        dest="library",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Filepath to a custom HMM library for the hmmsearch"
        if show_all
        else argparse.SUPPRESS,
    )
    inputdef.add_argument(
        "-r",
        dest="result_files_directory",
        type=path_str,
        metavar="<directory>",
        default=os.path.join(_project_root_from_this_file(), "results"),
        help="Directory for the result files",
    )
    inputdef.add_argument(
        "-db",
        dest="database_directory",
        type=path_str,
        metavar="<filepath>",
        help="Filepath to sqlite database (created if missing)",
    )
    inputdef.add_argument(
        "-c",
        dest="cores",
        type=int,
        default=4,
        metavar="<int>",
        help="Allocated CPU cores" if show_all else argparse.SUPPRESS,
    )
    inputdef.add_argument(
        "-glob_report",
        dest="glob_report",
        type=file_path,
        metavar="<filepath>",
        help="Filepath to glob hmmreport. Each report with one HMM queried against the concatenated genomes."
        if show_all
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

    # Search parameters
    parameters = parser.add_argument_group("Search parameters")
    parameters.add_argument(
        "-cut_type",
        dest="threshold_type",
        type=int,
        default=1,
        metavar="<int>",
        choices=[1, 2, 3],
        help="Choice of cutoff: 1 optimized; 2 trusted; 3 noise"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-cut_score",
        dest="thrs_score",
        type=int,
        default=50,
        metavar="<int>",
        help="Global minimal score cutoff" if show_all else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-taxonomy_info",
        dest="taxonomy_file",
        type=file_path,
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
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-n",
        dest="name",
        type=str,
        default="project",
        metavar="<string>",
        help="Name new project" if show_all else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-s",
        dest="stage",
        type=int,
        default=0,
        choices=[0, 1, 2, 3, 4, 5],
        help="Start at step" if show_all else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-x",
        dest="exit",
        type=int,
        default=10,
        choices=[0, 1, 2, 3, 4, 5],
        help="Exit at step" if show_all else argparse.SUPPRESS,
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
        help="Limit to HMM sets (whitespace or CSV separated)"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-clean",
        dest="clean_reports",
        action="store_true",
        help="Clean up any pre-existing HMMreport file"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-no_reports",
        dest="individual_reports",
        action="store_false",
        help="Do not write individual files per genome"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-max_seq_per_genome",
        dest="max_seqs_per_genome",
        type=int,
        default=4,
        help="Max. number of sequences per protein per genome for cross check via Diamond"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-no_cross_check",
        dest="bool_cross_check",
        action="store_false",
        help="No cross check with reference sequences via Diamond"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-optimized_cutoff_cross_check",
        dest="optimized_cutoff_cross_check",
        action="store_true",
        help="Use optimized cutoff instead of cross check with Diamond"
        if show_all
        else argparse.SUPPRESS,
    )

    # Synteny options
    synteny = parser.add_argument_group("Synteny options")
    synteny.add_argument(
        "-p",
        dest="patterns_file",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Filepath to patterns file" if show_all else argparse.SUPPRESS,
    )
    synteny.add_argument(
        "-cooccurrence",
        dest="cooccurrence_file",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Filepath to co-occurrence file" if show_all else argparse.SUPPRESS,
    )
    synteny.add_argument(
        "-exclude_singletons",
        dest="exclusion_singletons",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Filepath to tab separated file for singletons that are excluded"
        if show_all
        else argparse.SUPPRESS,
    )
    synteny.add_argument(
        "-mc",
        dest="min_completeness",
        type=float,
        default=0.5,
        metavar="<float>",
        help="Minimal fraction of predefined csb to be recognized"
        if show_all
        else argparse.SUPPRESS,
    )
    synteny.add_argument(
        "-chunks",
        dest="glob_chunks",
        type=int,
        default=5000,
        metavar="<int>",
        help="Chunk size for parsing results from glob before entering into database"
        if show_all
        else argparse.SUPPRESS,
    )

    # Information on resources
    information = parser.add_argument_group("Information on resources")
    information.add_argument(
        "-stat_keywords",
        action="store_true",
        help="Print patterns for keyword naming" if show_all else argparse.SUPPRESS,
    )
    information.add_argument(
        "-stat_csb",
        action="store_true",
        help="Print automatically found csbs" if show_all else argparse.SUPPRESS,
    )
    information.add_argument(
        "-stat_genomes",
        action="store_true",
        help="Print taxonomy information from database"
        if show_all
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
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-insertions",
        dest="insertions",
        type=int,
        default=1,
        metavar="<int>",
        help="Max. insertions in a csb" if show_all else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-occurence",
        dest="occurence",
        type=int,
        default=1,
        metavar="<int>",
        help="Min. number occurences to be recognized as csb"
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-min_csb_size",
        dest="min_csb_size",
        type=int,
        default=4,
        metavar="<int>",
        help="Min. number of genes in a csb before recognized"
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-max_csb_size",
        dest="max_csb_size",
        type=int,
        default=50,
        metavar="<int>",
        help="Max. number of genes in a csb before recognized"
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-max_gene_repeats",
        dest="max_domain_repeats",
        type=int,
        default=4,
        metavar="<int>",
        help="Maximum number of repeated genes in a csb."
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-jaccard",
        dest="jaccard",
        type=float,
        default=0.0,
        metavar="<float>",
        help="Acceptable dissimilarity in jaccard clustering [0.0-1.0]"
        if show_all
        else argparse.SUPPRESS,
    )

    # Work step regulation
    flow = parser.add_argument_group("Work step regulation")
    flow.add_argument(
        "-redo_taxonomy",
        dest="redo_taxonomy",
        action="store_true",
        help="Redo the taxonomy assignment*" if show_all else argparse.SUPPRESS,
    )

    # Limiter (dataset conditions)
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
        if show_all
        else argparse.SUPPRESS,
    )
    limiter.add_argument(
        "-dlt",
        dest="dataset_limit_taxon",
        type=str,
        default=None,
        metavar="<string>",
        help="Taxonomic name e.g. Proteobacteria. Requires -dll"
        if show_all
        else argparse.SUPPRESS,
    )
    limiter.add_argument(
        "-dlp",
        dest="dataset_limit_proteins",
        type=str,
        default="0",
        metavar="<list>",
        help="Limit fetch to genomes with <protein>" if show_all else argparse.SUPPRESS,
    )
    limiter.add_argument(
        "-dlk",
        dest="dataset_limit_keywords",
        type=str,
        default="0",
        metavar="<list>",
        help="Limit fetch to genomes with <keyword>" if show_all else argparse.SUPPRESS,
    )
    limiter.add_argument(
        "-dtd",
        dest="dataset_divide_sign",
        default=".",
        type=str,
        metavar="<string>",
        help='Separator for taxonomy information. The characters ";" ":" and "," cause strange behavior'
        if show_all
        else argparse.SUPPRESS,
    )

    # Output operators / fetch
    operators = parser.add_argument_group("Output sequences with these conditions *")
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
        if show_all
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
        if show_all
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
        if show_all
        else argparse.SUPPRESS,
    )
    operators.add_argument(
        "-kc",
        dest="keywords_connector",
        type=str,
        default="OR",
        choices=["AND", "OR"],
        help="Select cluster with keywords connected by AND or OR"
        if show_all
        else argparse.SUPPRESS,
    )

    # Alignment and sequence file processing
    process = parser.add_argument_group("Alignment and sequence file processing")
    process.add_argument(
        "-merge_fasta",
        dest="merge_fasta",
        type=dir_path,
        metavar="<directory>",
        help="Merges two or more sequence files with extension .faa without duplicates"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-filter_fasta",
        dest="filter_fasta",
        nargs=3,
        metavar=("FILE", "MIN", "MAX"),
        help="Filter FASTA by length MIN..MAX; write to FILE"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-concat_alignment",
        dest="concat_alignment",
        type=dir_path,
        metavar="<directory>",
        help="Concatenates alignment files with extension .fasta_aln"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-add_taxonomy_to_alignment",
        dest="add_taxonomy",
        type=file_path,
        metavar="<file> or <directory>",
        help="Adds taxonomy to alignment files in <dir>, requires -db with taxonomy"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-add_genomic_context",
        dest="add_genomic_context",
        type=file_path,
        metavar="<file>",
        help="Adds genomic context to sequences from fasta file, requires -db with taxonomy"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-create_type_range_dataset",
        dest="create_type_range_dataset",
        type=path_str,
        metavar="<file>",
        help="Create protein type range dataset from sequences fasta file, requires -db with taxonomy"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-create_gene_cluster_dataset",
        dest="create_gene_cluster_dataset",
        type=path_str,
        metavar="<file>",
        help="Create gene cluster dataset from sequences fasta file, requires -db with taxonomy"
        if show_all
        else argparse.SUPPRESS,
    )
    process.add_argument(
        "-aln_gaps",
        dest="gaps",
        action="store_true",
        help="When concatenating alignments add gaps for missing sequences"
        if show_all
        else argparse.SUPPRESS,
    )

    # ---- Hilfe? Keine Args? ----
    if show_all:
        parser.print_help()
        sys.exit(0)
    if len(clean_args) == 0:
        parser.print_help()
        sys.exit("Please provide arguments. For help use -h or --help-all")

    # ---- In Options(HMSSS)-Objekt einlesen ----
    options = Options()
    parser.parse_args(clean_args, namespace=options)

    # ---- Ressourcen-Defaults (HMSS2/data/...) ----
    data = _data_dir()
    options.location = _project_root_from_this_file()
    options.reference_seq_dir = os.path.join(data, "RefSeqs")

    if options.score_threshold_file is None:
        options.score_threshold_file = os.path.join(data, "Thresholds")
    if options.library is None:
        options.library = os.path.join(data, "HMMlib")
    if options.patterns_file is None:
        options.patterns_file = os.path.join(data, "Patterns")
    if options.cooccurrence_file is None:
        options.cooccurrence_file = os.path.join(data, "Cooccurrence")
    if options.exclusion_singletons is None:
        options.exclusion_singletons = os.path.join(data, "Exclusion_singletons")

    # result_files_directory absolut machen
    if not os.path.isabs(options.result_files_directory):
        options.result_files_directory = os.path.abspath(options.result_files_directory)

    # Default-Resultpfad = neues Projekt?
    default_results = os.path.join(_project_root_from_this_file(), "results")
    options.new_project = options.result_files_directory == default_results

    # ---- Betriebsmodi/Stages erkennen (wie in deiner main) ----
    # Parser-Defaults einsammeln (nur relevante Felder)
    default_values = {
        action.dest: action.default
        for action in parser._actions
        if getattr(action, "dest", None) and action.dest != "help"
    }
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
    default_values = {k: v for k, v in default_values.items() if k in relevant_groups}

    process_args = project.any_process_args_provided(
        options, default_values
    )  # returns bool

    # Spezial-Stage 100 (nur Taxonomy)
    if (
        options.taxonomy_file
        and options.database_directory
        and not options.fasta_file_directory
    ):
        options.stage = 100

    # Stage 101: Processing/Fetch/Operators
    if process_args:
        options.stage = 101

        # Fetch?
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

        # Processing?
        if (
            options.filter_fasta
            or options.concat_alignment
            or options.merge_fasta
            or options.add_taxonomy
            or options.add_genomic_context
            or options.create_type_range_dataset
            or options.create_gene_cluster_dataset
        ):
            options.process = True
            # Für diese hier braucht man DB
            if (
                options.add_taxonomy
                or options.add_genomic_context
                or options.create_type_range_dataset
                or options.create_gene_cluster_dataset
            ):
                if not options.database_directory:
                    sys.exit("Please use the -db argument to provide a valid database")

    # Normalisierung einiger Felder
    if options.HMM_sets:
        options.HMM_sets = _list_from_csv_or_repeat(options.HMM_sets)

    if options.filter_fasta:
        # ['file', '100', '200'] -> ['file', 100, 200]
        f, lo, hi = options.filter_fasta
        options.filter_fasta = [str(f), int(lo), int(hi)]

    return options

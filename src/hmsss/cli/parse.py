# src/hmsss/cli/parse.py
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Sequence, Tuple, Any

from hmsss.core.config import (
    Config, PathsCfg,
    CliInput, CliSearchParams, CliResources, CliSynteny, CliInfo,
    CliCsb, CliFlow, CliLimiter, CliOperators, CliProcess,
)
from hmsss.cli import paths as paths
#from hmsss.core.options import Hmsss as Options
from hmsss.db import project as project

# ---------------------------------------------------------------------------
# Hilfsroutinen (ersetzen die bisher in myUtil verwendeten argparse-Validatoren)
# ---------------------------------------------------------------------------


def dir_path(p: str) -> str:
    """Verlangt ein existierendes Verzeichnis; gibt den absoluten Pfad zurück. Prüft die existenz des Pfads."""
    ap = os.path.abspath(p)
    if not os.path.isdir(ap):
        raise argparse.ArgumentTypeError(f"directory does not exist: {p}")
    return ap


def file_path(p: str) -> str:
    """Verlangt existierenden Pfad (Datei ODER Verzeichnis); gibt den absoluten Pfad zurück. Prüft die existenz des Pfads."""
    ap = os.path.abspath(p)
    if not os.path.exists(ap):
        raise argparse.ArgumentTypeError(f"path does not exist: {p}")
    return ap


def path_str(p: str) -> str:
    """Nur Normalisierung: absoluter Pfad; Existenz wird NICHT geprüft (für -r/-db/Output-Ziele)."""
    return os.path.abspath(p)

def _list_from_csv_or_repeat(values: List[str]) -> List[str]:
    """Erlaubt -hmms A B C oder -hmms A,B,C."""
    out: List[str] = []
    for v in values or []:
        out.extend([p for p in v.split(",") if p])
    return out


# ---------------------------------------------------------------------------
# Hauptfunktion: parse_arguments
# ---------------------------------------------------------------------------


def parse_arguments(*,show_all: bool = False) -> argparse.ArgumentParser:
    """
    Baut den Argumentparser für das Programm
    """

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
        type=dir_path,
        metavar="<directory>",
        default=None,
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
        help="Overwrite pre-existing hmmsearch report file"
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


    return parser

def parse_cli(argv: Sequence[str] | None = None) -> argparse.Namespace:

    argv = list(argv) if argv is not None else sys.argv[1:]
    show_all = "--help-all" in argv

    parser = parse_arguments(show_all=show_all)  # baut NUR den Parser

    if show_all:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args(argv)

def _s(ns: Any, name: str) -> str | None:
    """Hilfsfunktion: CLI-Wert als String normalisieren (type= greift nicht auf Defaults)."""
    v = getattr(ns, name, None)
    return os.fspath(v) if v is not None else None

def _paths_cfg_from_paths_module() -> PathsCfg:
    """ Erzeugt ein config objekt mit den konstanten paths"""
    d = paths.as_dict(str_paths=True)
    return PathsCfg(
        root=d["ROOT_DIR"], bin=d["BIN_DIR"], data=d["DATA_DIR"],
        hmms=d["HMMS_DIR"], refseq=d["REFSEQ_DIR"],
        results=d["RESULTS_DIR"], package=d["PACKAGE_DIR"],
    )

def build_config_from_namespace(ns) -> Config:
    """
    Übersetzt argparse.Namespace → Config (ohne Legacy).
    Setzt fehlende Pfad-Defaults aus hmsss.cli.paths.
    """
    paths_cfg = _paths_cfg_from_paths_module()

    # ---------- CLI-Blöcke ----------
    cli_input = CliInput(
        fasta_file_directory=_s(ns, "fasta_file_directory"),
        score_threshold_file=_s(ns, "score_threshold_file") or os.path.join(paths_cfg.data, "Thresholds"),
        library=_s(ns, "library") or paths_cfg.hmms,
        result_files_directory=_s(ns, "result_files_directory") or paths_cfg.results,
        database_directory=_s(ns, "database_directory"),
        cores=int(getattr(ns, "cores", 4)),
        glob_report=_s(ns, "glob_report"),
        verbose=int(getattr(ns, "verbose", 1)),
    )

    cli_params = CliSearchParams(
        threshold_type=int(getattr(ns, "threshold_type", 1)),
        thrs_score=float(getattr(ns, "thrs_score", 50)),
        taxonomy_file=_s(ns, "taxonomy_file"),
        refseq_identity=int(getattr(ns, "refseq_identity", 90)),
        name=getattr(ns, "name", "project"),
        stage=int(getattr(ns, "stage", 0)),
        exit=int(getattr(ns, "exit", 10)),
    )

    cli_resources = CliResources(
        HMM_sets=list(getattr(ns, "HMM_sets", [])),
        clean_reports=bool(getattr(ns, "clean_reports", False)),
        individual_reports=bool(getattr(ns, "individual_reports", True)),
        max_seqs_per_genome=int(getattr(ns, "max_seqs_per_genome", 4)),
        bool_cross_check=bool(getattr(ns, "bool_cross_check", True)),
        optimized_cutoff_cross_check=bool(getattr(ns, "optimized_cutoff_cross_check", False)),
    )

    cli_synteny = CliSynteny(
        patterns_file=_s(ns, "patterns_file") or os.path.join(paths_cfg.data, "Patterns"),
        cooccurrence_file=_s(ns, "cooccurrence_file") or os.path.join(paths_cfg.data, "Cooccurrence"),
        exclusion_singletons=_s(ns, "exclusion_singletons") or os.path.join(paths_cfg.data, "Exclusion_singletons"),
        min_completeness=float(getattr(ns, "min_completeness", 0.5)),
        glob_chunks=int(getattr(ns, "glob_chunks", 5000)),
    )

    cli_info = CliInfo(
        stat_keywords=bool(getattr(ns, "stat_keywords", False)),
        stat_csb=bool(getattr(ns, "stat_csb", False)),
        stat_genomes=bool(getattr(ns, "stat_genomes", False)),
    )

    cli_csb = CliCsb(
        nucleotide_range=int(getattr(ns, "nucleotide_range", 3500)),
        insertions=int(getattr(ns, "insertions", 1)),
        occurence=int(getattr(ns, "occurence", 1)),
        min_csb_size=int(getattr(ns, "min_csb_size", 4)),
        max_csb_size=int(getattr(ns, "max_csb_size", 50)),
        max_domain_repeats=int(getattr(ns, "max_domain_repeats", 4)),
        jaccard=float(getattr(ns, "jaccard", 0.0)),
    )

    cli_flow = CliFlow(redo_taxonomy=bool(getattr(ns, "redo_taxonomy", False)))

    cli_limiter = CliLimiter(
        dataset_limit_lineage=getattr(ns, "dataset_limit_lineage", None),
        dataset_limit_taxon=getattr(ns, "dataset_limit_taxon", None),
        dataset_limit_proteins=getattr(ns, "dataset_limit_proteins", "0"),
        dataset_limit_keywords=getattr(ns, "dataset_limit_keywords", "0"),
        dataset_divide_sign=getattr(ns, "dataset_divide_sign", "."),
    )

    cli_ops = CliOperators(
        fetch_genomes=list(getattr(ns, "fetch_genomes", [])),
        fetch_proteins=list(getattr(ns, "fetch_proteins", [])),
        fetch_csbs=list(getattr(ns, "fetch_csbs", [])),
        fetch_keywords=list(getattr(ns, "fetch_keywords", [])),
        keywords_connector=getattr(ns, "keywords_connector", "OR"),
    )

    cli_process = CliProcess(
        merge_fasta=_s(ns, "merge_fasta"),
        filter_fasta=getattr(ns, "filter_fasta", None),
        concat_alignment=_s(ns, "concat_alignment"),
        add_taxonomy=_s(ns, "add_taxonomy"),
        add_genomic_context=_s(ns, "add_genomic_context"),
        create_type_range_dataset=_s(ns, "create_type_range_dataset"),
        create_gene_cluster_dataset=_s(ns, "create_gene_cluster_dataset"),
        gaps=bool(getattr(ns, "gaps", False)),
    )

    cfg = Config(
        paths=paths_cfg,
        cli_input=cli_input,
        cli_params=cli_params,
        cli_resources=cli_resources,
        cli_synteny=cli_synteny,
        cli_info=cli_info,
        cli_csb=cli_csb,
        cli_flow=cli_flow,
        cli_limiter=cli_limiter,
        cli_ops=cli_ops,
        cli_process=cli_process,
    )
    cfg.validate()
    return cfg

def _needs_stage_100(ns: argparse.Namespace) -> bool:
    """ Prüfe ob ein prozess oder fetch argument oder taxonomie addition vorgebracht wurde """
    # (1) Fetch aus Datenbank angefordert?
    fetch_requested = any([
        bool(getattr(ns, "fetch_genomes", [])),
        bool(getattr(ns, "fetch_proteins", [])),
        bool(getattr(ns, "fetch_csbs", [])),
        bool(getattr(ns, "fetch_keywords", [])),
        bool(getattr(ns, "dataset_limit_lineage", None)),
        bool(getattr(ns, "dataset_limit_taxon", None)),
        getattr(ns, "dataset_limit_proteins", "0") not in (None, "0"),
        getattr(ns, "dataset_limit_keywords", "0") not in (None, "0"),
    ])

    # (2) Änderungen an FASTA/Alignment?
    processing_requested = any([
        bool(getattr(ns, "merge_fasta", None)),
        bool(getattr(ns, "filter_fasta", None)),              # list mit ["FILE","MIN","MAX"]
        bool(getattr(ns, "concat_alignment", None)),
        bool(getattr(ns, "add_taxonomy", None)),              # -add_taxonomy_to_alignment
        bool(getattr(ns, "add_genomic_context", None)),
        bool(getattr(ns, "create_type_range_dataset", None)),
        bool(getattr(ns, "create_gene_cluster_dataset", None)),
    ])

    # (3) redo taxonomy?
    redo_tax = bool(getattr(ns, "redo_taxonomy", False))

    return fetch_requested or processing_requested or redo_tax

def _apply_runtime_defaults(ns: argparse.Namespace) -> argparse.Namespace:
    """
    Füllt fehlende Pfad-Defaults aus hmsss.cli.paths und legt die Stage robust fest.
    (post-parse, damit type= Validierungen nicht auf „phantom defaults“ laufen)
    """
    # 1) Pfad-Defaults: nur setzen, wenn None/leer
    def _set_default(attr: str, value: str) -> None:
        v = getattr(ns, attr, None)
        if v is None or (isinstance(v, str) and v.strip() == ""):
            setattr(ns, attr, value)

    _set_default("score_threshold_file", str(paths.SRC_FILE_THRESHOLDS))
    _set_default("library",               str(paths.SRC_FILE_HMM_LIBRARY))
    _set_default("patterns_file",         str(paths.SRC_FILE_PATTERNS))
    _set_default("cooccurrence_file",     str(paths.SRC_FILE_COOCCURRENCE))
    _set_default("exclusion_singletons",  str(paths.SRC_FILE_EXCLUSION_SINGLETONS))
    _set_default("result_files_directory",str(paths.RESULTS_DIR))

    # Stage normalisieren oder auf 100 forcieren
    raw_stage = getattr(ns, "stage", None)
    try:
        normalized = int(raw_stage) if raw_stage is not None else 0
    except Exception:
        normalized = 0
    # clamp 0..5
    if normalized < 0: normalized = 0
    if normalized > 5: normalized = 5

    if _needs_stage_100(ns):
        setattr(ns, "stage", 100)
    else:
        setattr(ns, "stage", normalized)

    return ns

def parse_to_config(argv: list[str] | None = None) -> Config:
    """
    Komfort-Funktion: parst argv und liefert direkt eine fertige Config.
    """
    ns = parse_cli(argv)  # deine bestehende Routine
    ns = _apply_runtime_defaults(ns)

    return build_config_from_namespace(ns)

"""
Hier werden die Argumente von argparse bzw die default werte dieser Argumente addiert
Danach werden die default Argumente der Laufzeit hinzugefügt. Das betrifft Dateien, die
unter Umständen noch nicht vorhanden sind und daher nicht im Argparse auftauchen können
"""
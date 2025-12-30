# src/hmsss/cli/parse.py
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Sequence, Any

from hmsss.cli.config import (
    Config,
    PathsCfg,
    CliInput,
    CliSearchParams,
    CliResources,
    CliSynteny,
    CliInfo,
    CliCsb,
    CliFlow,
    CliOperators,
    CliReadMapping,
)
from hmsss.cli import paths as paths


"""
Argument parsing and configuration assembly for HMSS2/HMSSS.

This module builds the command-line interface (CLI), parses arguments,
and converts them into a structured `Config` object used by the pipeline.
It also injects runtime defaults that depend on the inferred project paths
(`hmsss.cli.paths`).

Typical usage:
    ns = parse_cli(sys.argv[1:])
    ns = _apply_runtime_defaults(ns)
    cfg = build_config_from_namespace(ns)
"""


# ---------------------------------------------------------------------------
# Utility routines
# ---------------------------------------------------------------------------


def dir_path(p: str) -> str:
    """Validate that `p` is an existing directory and return its absolute path.

    Args:
        p: Directory path provided by the user.

    Returns:
        Absolute path to the directory.

    Raises:
        argparse.ArgumentTypeError: If the directory does not exist.
    """
    ap = os.path.abspath(p)
    if not os.path.isdir(ap):
        raise argparse.ArgumentTypeError(f"directory does not exist: {p}")
    return ap


def file_path(p: str) -> str:
    """Validate that `p` is an existing filesystem path (file or directory).

    Args:
        p: File or directory path provided by the user.

    Returns:
        Absolute path to the file or directory.

    Raises:
        argparse.ArgumentTypeError: If the path does not exist.
    """
    ap = os.path.abspath(p)
    if not os.path.exists(ap):
        raise argparse.ArgumentTypeError(f"path does not exist: {p}")
    return ap


def path_str(p: str) -> str:
    """Validate that `p` is an existing filesystem path (file or directory).

    Args:
        p: File or directory path provided by the user.

    Returns:
        Absolute path to the file or directory.

    Raises:
        argparse.ArgumentTypeError: If the path does not exist.
    """
    return os.path.abspath(p)


def _list_from_csv_or_repeat(values: List[str]) -> List[str]:
    """Split CLI values by comma and flatten whitespace-separated lists.

    Supports both `-hmms A B C` and `-hmms A,B,C`.

    Args:
        values: Raw list of strings provided by argparse.

    Returns:
        A flattened list of items without empty entries.
    """
    out: List[str] = []
    for v in values or []:
        out.extend([p for p in v.split(",") if p])
    return out


# ---------------------------------------------------------------------------
# Main function: parse_arguments
# ---------------------------------------------------------------------------


def parse_arguments(*, show_all: bool = False) -> argparse.ArgumentParser:
    """Construct the top-level argument parser for HMSS2/HMSSS.

    The parser is organized into semantic groups (input, search, resources,
    synteny, CSB prediction, flow/limiters, operators, processing). If
    `show_all` is True, advanced/rarely used options are included in `--help`.

    Args:
        show_all: Include advanced options in help output.

    Returns:
        A fully configured `argparse.ArgumentParser` (no parsing yet).
    """

    never_show = False  # Always hide these arguments, because currently unused

    # ---- Argument groups ----
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
        "-r",
        dest="result_files_directory",
        type=dir_path,
        metavar="<directory>",
        default=None,
        help="Directory for the result files",
    )
    inputdef.add_argument(
        "-c",
        dest="cores",
        type=int,
        default=4,
        metavar="<int>",
        help="Allocated CPU cores" if show_all else argparse.SUPPRESS,
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
        "-taxonomy",
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
        default=["DHPS", "DMS", "Dsr", "SQ", "Sulfonates"],
        metavar="<list>",
        help="Limit to HMM sets (whitespace or CSV separated)"
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
        help="Max. number of initial hits per protein per genome forwarded to cross check via Diamond"
        if show_all
        else argparse.SUPPRESS,
    )
    resources.add_argument(
        "-diamond_speed",
        dest="diamond_speed_mode",
        type=str,
        choices=[
            "faster",
            "fast",
            "mid-sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ],
        default="faster",
        help="DIAMOND blastp speed mode" if show_all else argparse.SUPPRESS,
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
        if never_show
        else argparse.SUPPRESS,
    )

    # Synteny options
    # synteny = parser.add_argument_group("Synteny options")

    # Information on resources
    information = parser.add_argument_group("Information on resources")
    information.add_argument(
        "-stat_keywords",
        action="store_true",
        help="Print patterns for keyword naming" if never_show else argparse.SUPPRESS,
    )
    information.add_argument(
        "-stat_csb",
        action="store_true",
        help="Print automatically found csbs" if never_show else argparse.SUPPRESS,
    )
    information.add_argument(
        "-stat_genomes",
        action="store_true",
        help="Print taxonomy information from database"
        if show_all
        else argparse.SUPPRESS,
    )

    # Collinear syntenic block prediction
    csb = parser.add_argument_group("Collinear syntenic block detection")
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
        "-mc",
        dest="min_completeness",
        type=float,
        default=0.51,
        metavar="<float>",
        help="Minimal fraction of predefined csb to be recognized"
        if show_all
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-insertions",
        dest="insertions",
        type=int,
        default=1,
        metavar="<int>",
        help="Max. insertions in a csb" if never_show else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-occurrence",
        dest="occurence",
        type=int,
        default=1,
        metavar="<int>",
        help="Min. number occurrences to be recognized as csb"
        if never_show
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-min_csb_size",
        dest="min_csb_size",
        type=int,
        default=4,
        metavar="<int>",
        help="Min. number of genes in a csb before recognized"
        if never_show
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-max_csb_size",
        dest="max_csb_size",
        type=int,
        default=50,
        metavar="<int>",
        help="Max. number of genes in a csb before recognized"
        if never_show
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-max_gene_repeats",
        dest="max_domain_repeats",
        type=int,
        default=4,
        metavar="<int>",
        help="Maximum number of repeated genes in a csb."
        if never_show
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-jaccard",
        dest="jaccard",
        type=float,
        default=0.0,
        metavar="<float>",
        help="Acceptable dissimilarity in jaccard clustering [0.0-1.0]"
        if never_show
        else argparse.SUPPRESS,
    )
    csb.add_argument(
        "-chunks",
        dest="glob_chunks",
        type=int,
        default=5000,
        metavar="<int>",
        help="Chunk size for parsing results from glob before entering into database"
        if never_show
        else argparse.SUPPRESS,
    )

    # Read mapping workflow algorithm
    readmap = parser.add_argument_group("Read mapping integration")

    readmap.add_argument(
        "--read-mapping",
        dest="use_read_mapping",
        action="store_true",
        help="Enable read-mapping analysis for FASTQ files located in the -f directory.",
    )

    readmap.add_argument(
        "--rm-threads",
        dest="rm_threads",
        type=int,
        default=14,
        metavar="<int>",
        help="Number of threads to use (default: 14).",
    )

    readmap.add_argument(
        "--rm-evalue",
        dest="rm_evalue",
        type=float,
        default=1e-5,
        metavar="<float>",
        help="E-value threshold for homology search (default: 1e-5).",
    )

    readmap.add_argument(
        "--rm-placements-cutoff",
        dest="rm_placements_cutoff",
        type=float,
        default=0.75,
        metavar="<float>",
        help="Placement cutoff for phylogenetic placement (default: 0.75).",
    )

    readmap.add_argument(
        "--rm-resolve-placements",
        dest="rm_resolve_placements",
        action="store_true",
        help="Resolve ambiguous phylogenetic placements (default: False).",
    )

    readmap.add_argument(
        "--rm-min-orf-length",
        dest="rm_min_orf_length",
        type=int,
        default=96,
        metavar="<int>",
        help="Minimum ORF length (default: 96).",
    )

    readmap.add_argument(
        "--rm-restrict-read-length",
        dest="rm_restrict_read_length",
        type=int,
        default=None,
        metavar="<int>",
        help="Maximum read length (default: None).",
    )

    readmap.add_argument(
        "--rm-translation-table",
        dest="rm_translation_table",
        type=int,
        default=11,
        metavar="<int>",
        help="NCBI translation table to use (default: 11).",
    )

    # Work step regulation
    flow = parser.add_argument_group("Work step regulation")
    flow.add_argument(
        "-no_synteny_completion",
        dest="use_synteny_completion",
        action="store_false",
        help="Use syntenic block completeness enhancement",
    )
    flow.add_argument(
        "-no_remove_intermediate",
        dest="use_remove_unassigned_intermediates",
        action="store_false",
        help="Remove intermediate hits without genetic context",
    )
    flow.add_argument(
        "-no_remove_exclusion_singletons",
        dest="use_remove_exclusion_singletons",
        action="store_false",
        help="Remove genes that not occur as singletons",
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
        help="Select only csb encoding the given proteins (whitespace separated). The ':' without whitespace will be interpreted as logical OR",
    )
    operators.add_argument(
        "-fnd",
        nargs="+",
        dest="exclude_domains",
        type=str,
        default=[],
        metavar="<list>",
        help="Select gene cluster without these proteins (whitespace separated). The ':' without whitespace will be interpreted as logical OR"
        if show_all
        else argparse.SUPPRESS,
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
    operators.add_argument(
        "-include_noise_cut_hits",
        dest="use_valid_hits",
        action="store_false",
        help="Include distant homologs that are considered as noise",
    )
    operators.add_argument(
        "-fasta",
        dest="print_fasta",
        action="store_true",
        help="Print protein sequence fasta files for retrieved hits"
        if show_all
        else argparse.SUPPRESS,
    )
    operators.add_argument(
        "--print-graphs",
        dest="print_graphs",
        action="store_true",
        help="Print graphs for the selected output",
    )
    operators.add_argument(
        "--graph-tax-levels",
        dest="graph_tax_levels",
        nargs="+",
        type=str,
        metavar="<taxlevel>",
        choices=[
            "Superkingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ],
        default=["Phylum"],
        help=(
            "Taxonomic levels to summarize in graphs. Multiple selection is possible. "
            "[Superkingdom, Phylum, Class, Order, Family, Genus, Species]"
        ),
    )

    return parser


def parse_cli(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments and return the argparse namespace.

    If `--help-all` is present, help is printed and the process exits.

    Args:
        argv: Argument vector (defaults to `sys.argv[1:]`).

    Returns:
        Namespace with validated and typed CLI options.
    """
    argv = list(argv) if argv is not None else sys.argv[1:]
    show_all = "--help-all" in argv

    parser = parse_arguments(show_all=show_all)  # baut NUR den Parser

    if show_all:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args(argv)


def _s(ns: Any, name: str) -> str | None:
    """Fetch an attribute from `namespace` and coerce it to a filesystem string.

    Args:
        ns: Parsed argparse namespace.
        name: Attribute name to fetch.

    Returns:
        String path if present, otherwise `None`.
    """
    v = getattr(ns, name, None)
    return os.fspath(v) if v is not None else None


def _paths_cfg_from_paths_module() -> PathsCfg:
    """Build a `PathsCfg` from `hmsss.cli.paths` constants.

    Returns:
        A `PathsCfg` instance that centralizes project directories resolved at import time.
    """
    d = paths.as_dict(str_paths=True)
    return PathsCfg(
        root=d["ROOT_DIR"],
        bin=d["BIN_DIR"],
        data=d["DATA_DIR"],
        hmms=d["HMMS_DIR"],
        refseq=d["REFSEQ_DIR"],
        results=d["RESULTS_DIR"],
        package=d["PACKAGE_DIR"],
    )


def build_config_from_namespace(ns) -> Config:
    """Assemble the high-level `Config` object from parsed arguments.

    This function maps argparse fields into structured dataclasses grouped
    by concern (input, search params, resources, synteny, info, CSB, flow,
    limiters, operators, processing). Path defaults are injected from
    `hmsss.cli.paths` when missing.

    Args:
        ns: Argparse namespace produced by `parse_cli` / `_apply_runtime_defaults`.

    Returns:
        A validated `Config` instance ready for the pipeline.

    Raises:
        ValueError: If semantic constraints are violated.
        FileNotFoundError: If required project directories are missing during validation.
    """
    paths_cfg = _paths_cfg_from_paths_module()

    # ---------- CLI-Blöcke ----------
    cli_input = CliInput(
        fasta_file_directory=_s(ns, "fasta_file_directory"),
        score_threshold_file=_s(ns, "score_threshold_file")
        or os.path.join(paths_cfg.data, "Thresholds"),
        library=_s(ns, "library") or paths_cfg.hmms,
        result_files_directory=_s(ns, "result_files_directory") or paths_cfg.results,
        cores=int(getattr(ns, "cores", 4)),
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
        individual_reports=bool(getattr(ns, "individual_reports", True)),
        max_seqs_per_genome=int(getattr(ns, "max_seqs_per_genome", 4)),
        diamond_speed_mode=str(getattr(ns, "diamond_speed_mode", "faster")),
        bool_cross_check=bool(getattr(ns, "bool_cross_check", True)),
        optimized_cutoff_cross_check=bool(
            getattr(ns, "optimized_cutoff_cross_check", False)
        ),
    )

    cli_synteny = CliSynteny(
        patterns_file=_s(ns, "patterns_file")
        or os.path.join(paths_cfg.data, "Patterns"),
        cooccurrence_file=_s(ns, "cooccurrence_file")
        or os.path.join(paths_cfg.data, "Cooccurrence"),
        exclusion_singletons=_s(ns, "exclusion_singletons")
        or os.path.join(paths_cfg.data, "Exclusion_singletons"),
        min_completeness=float(getattr(ns, "min_completeness", 0.5)),
        glob_chunks=int(getattr(ns, "glob_chunks", 5000)),
    )

    cli_info = CliInfo(
        stat_keywords=bool(getattr(ns, "stat_keywords", False)),
        stat_csb=bool(getattr(ns, "stat_csb", False)),
        stat_genomes=bool(getattr(ns, "stat_genomes", False)),
        metabolic_information=_s(ns, "metabolism_information"),
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

    cli_readmap = CliReadMapping(
        use_read_mapping=bool(getattr(ns, "use_read_mapping", False)),
        threads=ns.rm_threads,
        evalue=ns.rm_evalue,
        placements_cutoff=ns.rm_placements_cutoff,
        resolve_placements=bool(getattr(ns, "rm_resolve_placements", False)),
        min_orf_length=ns.rm_min_orf_length,
        restrict_read_length=ns.rm_restrict_read_length,
        translation_table=ns.rm_translation_table,
    )

    cli_flow = CliFlow(
        redo_taxonomy=bool(getattr(ns, "redo_taxonomy", False)),
        use_synteny_completion=bool(getattr(ns, "use_synteny_completion", True)),
        use_remove_unassigned_intermediates=bool(
            getattr(ns, "use_remove_unassigned_intermediates", True)
        ),
        use_remove_exclusion_singletons=bool(
            getattr(ns, "use_remove_exclusion_singletons", True)
        ),
    )

    cli_ops = CliOperators(
        fetch_genomes=list(getattr(ns, "fetch_genomes", [])),
        fetch_proteins=list(getattr(ns, "fetch_proteins", [])),
        fetch_csbs=list(getattr(ns, "fetch_csbs", [])),
        fetch_not_csb_with_these_domains=list(getattr(ns, "exclude_domains", [])),
        fetch_keywords=list(getattr(ns, "fetch_keywords", [])),
        keywords_connector=getattr(ns, "keywords_connector", "OR"),
        print_fasta=getattr(ns, "print_fasta", False),
        print_graphs=bool(getattr(ns, "print_graphs", False)),
        use_valid_hits=bool(getattr(ns, "use_valid_hits", True)),
        graph_tax_levels=list(getattr(ns, "graph_tax_levels", ["Phylum"])),
    )

    cfg = Config(
        paths=paths_cfg,
        cli_input=cli_input,
        cli_params=cli_params,
        cli_resources=cli_resources,
        cli_synteny=cli_synteny,
        cli_info=cli_info,
        cli_csb=cli_csb,
        cli_readmap=cli_readmap,
        cli_flow=cli_flow,
        cli_ops=cli_ops,
    )
    cfg.validate()
    return cfg


def _needs_stage_50(ns: argparse.Namespace) -> bool:
    """Determine whether stage 50 (read mapping redo) must be forced.

    Returns:
        True if read mapping was requested.
    """
    return bool(getattr(ns, "use_read_mapping", False))


def _needs_stage_100(ns: argparse.Namespace) -> bool:
    """Force stage 100 when user requests taxonomy add/update only.

    Rule:
      - taxonomy_file is provided AND
      - result_files_directory is provided AND
      - fasta_file_directory is NOT provided
    """
    taxonomy_file = getattr(ns, "taxonomy_file", None)
    results_dir = getattr(ns, "result_files_directory", None)
    fasta_dir = getattr(ns, "fasta_file_directory", None)

    return bool(taxonomy_file) and bool(results_dir) and not bool(fasta_dir)


def _needs_stage_101(ns: argparse.Namespace) -> bool:
    """Determine whether stage 101 is required for dataset/output/processing.

    Stage 101 is used for database fetch operations and FASTA/alignment utilities.

    Returns:
        True if any fetch/output/processing option has been requested.
    """
    # Fetch from database
    fetch_requested = any(
        [
            bool(getattr(ns, "fetch_genomes", [])),
            bool(getattr(ns, "fetch_proteins", [])),
            bool(getattr(ns, "fetch_csbs", [])),
            bool(getattr(ns, "fetch_keywords", [])),
            bool(getattr(ns, "dataset_limit_lineage", None)),
            bool(getattr(ns, "dataset_limit_taxon", None)),
            getattr(ns, "dataset_limit_proteins", "0") not in (None, "0"),
            getattr(ns, "dataset_limit_keywords", "0") not in (None, "0"),
        ]
    )

    # Modify fasta/alignment files
    processing_requested = any(
        [
            bool(getattr(ns, "merge_fasta", None)),
            bool(getattr(ns, "filter_fasta", None)),  # list mit ["FILE","MIN","MAX"]
            bool(getattr(ns, "concat_alignment", None)),
            bool(getattr(ns, "add_taxonomy", None)),  # -add_taxonomy_to_alignment
            bool(getattr(ns, "add_genomic_context", None)),
            bool(getattr(ns, "create_type_range_dataset", None)),
            bool(getattr(ns, "create_gene_cluster_dataset", None)),
        ]
    )

    return fetch_requested or processing_requested


def _apply_runtime_defaults(ns: argparse.Namespace) -> argparse.Namespace:
    """Inject runtime defaults for paths and normalize the `stage`.

    This runs *after* parsing to avoid triggering argparse validators on
    non-existent files that are meant to be created later. It also forces
    `stage=100` or `stage=101` when the requested actions require it.

    Args:
        ns: Parsed argparse namespace.

    Returns:
        The mutated namespace with defaults applied and stage normalized.
    """

    # 1) default paths for undefined argparse paths
    def _set_default(attr: str, value: str) -> None:
        v = getattr(ns, attr, None)
        if v is None or (isinstance(v, str) and v.strip() == ""):
            setattr(ns, attr, value)

    _set_default("score_threshold_file", str(paths.SRC_FILE_THRESHOLDS))
    _set_default("library", str(paths.SRC_FILE_HMM_LIBRARY))
    _set_default("patterns_file", str(paths.SRC_FILE_PATTERNS))
    _set_default("cooccurrence_file", str(paths.SRC_FILE_COOCCURRENCE))
    _set_default("exclusion_singletons", str(paths.SRC_FILE_EXCLUSION_SINGLETONS))
    _set_default("metabolism_information", str(paths.SRC_FILE_METABOLISM_INFORMATION))
    _set_default("result_files_directory", str(paths.RESULTS_DIR))

    # Stage normalization
    raw_stage = getattr(ns, "stage", None)
    try:
        normalized_stage = int(raw_stage) if raw_stage is not None else 0
    except Exception:
        normalized_stage = 0
    # clamp 0..5
    if normalized_stage < 0:
        normalized_stage = 0
    if normalized_stage > 5:
        normalized_stage = 5

    if _needs_stage_100(ns):
        # Required for redo taxonomy command
        setattr(ns, "stage", 100)
    elif _needs_stage_101(ns):
        # Required for all dataset, output and processing commands
        setattr(ns, "stage", 101)
    elif _needs_stage_50(ns):
        # Required for read-mapping workflow
        setattr(ns, "stage", 50)
    else:
        setattr(ns, "stage", normalized_stage)

    return ns


def parse_to_config(argv: list[str] | None = None) -> Config:
    """Convenience wrapper: parse arguments, apply defaults, and build `Config`.

    Args:
        argv: Optional CLI arguments; defaults to `sys.argv[1:]`.

    Returns:
        A fully validated `Config` object.
    """
    namespace = parse_cli(argv)
    namespace = _apply_runtime_defaults(namespace)
    config = build_config_from_namespace(namespace)

    return config


"""
Hier werden die Argumente von argparse bzw die default werte dieser Argumente addiert
Danach werden die default Argumente der Laufzeit hinzugefügt. Das betrifft Dateien, die
unter Umständen noch nicht vorhanden sind und daher nicht im Argparse auftauchen können
"""

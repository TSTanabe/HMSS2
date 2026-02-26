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


class FetchHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.print_fetch_operator_help()
        parser.exit(0)

    def print_fetch_operator_help(self) -> None:
        """
        Print a concise explanation of how -fd and -fc are interpreted,
        including AND/OR logic, brackets, optional domains, and combinations.
        """

        text = """
    FETCH OPERATORS: -fd and -fc
    ===========================

    -fd (fetch domains / proteins)
    ------------------------------
    Selects proteins by domain name and restricts results by genome-level co-occurrence.

    Core meaning:
      -fd A B means: select genomes where A and B co-occur in the SAME genome.
      Output then contains the matching proteins/domains from those genomes (subject to other filters).

    Rules:
      - Multiple values are combined using logical AND at the genome level (co-occurrence in one genome).
      - Use ':' (without whitespace) to express OR within a single slot.
      - A trailing ':' makes a domain optional (can be omitted for that slot).
      - Square brackets [...] define alternative requirement groups (logical OR between groups).

    Examples:
      -fd A B
          -> genomes containing A AND B (both must be present in the same genome)

      -fd A:B C
          -> genomes containing (A OR B) AND C

      -fd A:B:C D
          -> genomes containing (A OR B OR C) AND D

      -fd A:
          -> genomes containing A OR nothing (A is optional)

      -fd A: B
          -> genomes containing (A optional) AND B

      -fd [A B] [C D]
          -> (genomes with A AND B) OR (genomes with C AND D)


    -fc (fetch gene clusters / CSBs)
    --------------------------------
    Selects genomes by gene clusters (CSBs) that encode the specified domains.

    Core meaning:
      -fc A B means: select genomes that have at least one gene cluster containing A and B together.

    Rules:
      - Within a cluster definition, whitespace-separated tokens are combined using logical AND.
      - ':' (without whitespace) expresses OR within a single token (alternative domains for that slot).
      - A trailing ':' makes that token optional (can be omitted).
      - Square brackets [...] define alternative cluster definitions (logical OR between groups).

    Examples:
      -fc A B
          -> genomes with a cluster containing A AND B

      -fc A:B C
          -> genomes with a cluster containing (A OR B) AND C

      -fc A:B:C D
          -> genomes with a cluster containing (A OR B OR C) AND D

      -fc A:
          -> genomes with a cluster containing A OR no constraint for that slot

      -fc [A B] [C D]
          -> genomes with a cluster containing (A AND B) OR (C AND D)


    -fc and -fd combined
    -------------------
    When used together:

      -fc defines which genomes are selected based on gene clusters (hard genome filter).
      -fd then fetches/adds proteins by domain, but ONLY within the genomes selected by -fc.
      -fd cannot introduce new genomes beyond the -fc selection.

    In short:
      -fc defines the genome set;
      -fd extends the protein output within that set.
    
    
    -fnd (exclude domains)
    ----------------------
    Exclude results based on the presence of domains.
    
    Core meaning:
      -fnd A B means: exclude any genome, gene cluster, or protein where A or B is present.
      Exclusion is triggered as soon as any specified domain occurs.
    
    Rules:
      - All specified domains are treated as independent exclusion criteria.
      - Multiple values are always combined using logical OR.
      - If any one domain matches, the result is excluded.
    
    Examples:
      -fnd A
          -> exclude results containing A
    
      -fnd A B
          -> exclude results containing A or B
    
    Interaction with other options:
      - Exclusion is applied after genome or cluster selection.
      - Results matching -fnd are removed even if they satisfy -fd or -fc.

    ==============
    Important note
    --------------
    The fetch operators (-fd, -fc, -fnd) are evaluated only when the -r option
    is provided.
    
    The -r option must point to an existing result directory from a previous run.
    This directory is used as the input source for all fetch operations.
    
    If -r is not specified or does not point to a valid result directory:
      - fetch conditions are ignored
      - no fetch-based summaries, FASTA files, or graphs are generated

    """

        print(text.strip())


class FetchOutputHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.print_fetch_output_help()
        parser.exit(0)

    def print_fetch_output_help(self) -> None:
        text = """
    OUTPUT FILES
    ============

    The fetch operation produces tabular summary reports, FASTA sequence files,
    and graphical summaries. All outputs are derived from the same final,
    filtered result set after applying -fd, -fc and -fnd.

    SUMMARY TABLES (TSV)
    --------------------

    summary_hit_table.txt
      Main per-protein hit table (genomeID/proteinID, domains, coordinates, cluster, taxonomy).

    summary_gene_taxonomy.txt
      Protein-to-taxonomy mapping (proteinID -> full lineage string).

    summary_unique_lineages.txt
      Unique taxonomy lineages with genome counts.

    summary_hit_taxonomy_counts.txt
      Taxonomy-level presence summary for all fetched proteins.

    summary_requested_hit_taxonomy_counts.txt
      Same as above, restricted to proteins explicitly requested via -fd.

    summary_metabolic_annotations.txt
      Metabolic/functional annotations per protein/domain plus taxonomy.

    summary_genecluster_overview_table.txt
      Overview of gene cluster compositions associated with requested domains.

    summary_strain_variability_by_species.txt
      Species-level completeness/variability for requested domains.

    FASTA OUTPUT
    ------------

    Protein_sequences/
      <DOMAIN>.faa
        All proteins containing <DOMAIN>.
      multi_domain_<DOMAIN>.faa
        Domain sequences extracted from fusion proteins.
      <NAME>_ortho.faa / <NAME>_paralog.faa
        Split by single-copy vs multi-copy per genome.

    GRAPHICAL OUTPUT (JPG)
    ----------------------

    summary_hit_taxonomy_counts.jpg
      Presence/absence and co-occurrence visualization for all fetched proteins.

    summary_requested_hit_taxonomy_counts.jpg
      Presence/absence visualization for requested proteins only.

    summary_gene_taxonomy.jpg
      Protein-to-taxonomy assignments.

    summary_unique_lineages.jpg
      Unique lineage overview with genome counts.

    summary_metabolic_annotations.jpg
      Metabolic/functional annotation overview across taxa.

    summary_genecluster_overview_table.jpg
      Gene cluster composition overview.

    COMMAND-LINE RECORD
    -------------------

    command_line_args.txt
      Exact command-line arguments used for the run (index + value).
    """
        print(text.strip())


class ReadMappingHelpAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        self.print_read_mapping_help()
        parser.exit(0)

    def print_read_mapping_help(self) -> None:
        # Defaults aus deinem Setup
        print("""
        READ MAPPING — AVAILABLE GPKG PACKAGES & RAM REQUIREMENTS
        =========================================================

        This list shows all available graftM packages (*.gpkg) for read mapping
        together with their estimated peak RAM usage (GB).

        RAM values are based on benchmark runs with:
          - threads = 8
          - reads   = 100
          - metric  = peak RSS (GB)

        Package sets are derived from the parent directory:
          gpkg/v10_gpkg_set_<SET>/.../<NAME>.gpkg


        SET: Apr_Qmo
        ------------
        AprM        : 10.47 GB
        QmoA        : 38.83 GB
        QmoB        : 49.98 GB
        QmoC        : 14.64 GB
        oxAprAI     : 24.14 GB
        oxAprAII    : 16.28 GB
        oxAprBI     : 4.53 GB
        oxAprBII    : 7.13 GB
        oxSat       : 21.58 GB
        qHdrB       : 6.71 GB
        qHdrC       : 11.57 GB
        redAprA     : 21.93 GB
        redAprB     : 2.19 GB
        redSat      : 16.59 GB


        SET: Asr_Mcc_Phs_Ttr
        -------------------
        AsrA        : 7.91 GB
        AsrB        : 6.48 GB
        AsrC        : 7.41 GB
        MccA        : 2.93 GB
        MccB        : 1.33 GB
        MccC        : 1.19 GB
        MccD        : 2.11 GB
        PhsA        : 29.35 GB
        PhsB        : 7.16 GB
        PhsC        : 12.16 GB
        TtrA        : 12.29 GB
        TtrB        : 3.64 GB
        TtrC        : 4.29 GB


        SET: CS_Aryl
        ------------
        AtsA        : 36.46 GB
        AtsB        : 32.47 GB
        CosH        : 61.94 GB
        Cs2H        : 1.24 GB
        DszA        : 5.77 GB
        DszB        : 1.38 GB
        DszC        : 4.97 GB
        DszD        : 0.25 GB
        SncA        : 0.58 GB
        SncB        : 0.66 GB
        SncC        : 0.97 GB
        SsuD        : 1.11 GB
        SsuE        : 0.66 GB
        TcdH        : 0.58 GB
        TcdS        : 0.54 GB
        TcdT        : 0.29 GB


        SET: DHPS_Taurine_Isethionate
        ----------------------------
        AdhE        : 15.01 GB
        ComC        : 49.28 GB
        ComD        : 2.70 GB
        ComE        : 3.13 GB
        HpfD        : 5.16 GB
        HpfG        : 27.43 GB
        HpfH        : 4.03 GB
        HpfX        : 6.42 GB
        HpfY        : 32.31 GB
        HpfZ        : 24.64 GB
        HpsG        : 26.94 GB
        HpsH        : 9.78 GB
        HpsN        : 46.88 GB
        HpsO        : 30.49 GB
        HpsP        : 27.35 GB
        IseJ        : 59.44 GB
        IsfD        : 11.50 GB
        SarD        : 5.08 GB
        SauS        : 10.17 GB
        SauT        : 15.21 GB
        SlcC        : 37.33 GB
        SlcD        : 21.90 GB
        SuyA        : 5.41 GB
        SuyB        : 29.91 GB
        TauX        : 2.25 GB
        TauY        : 7.98 GB
        Toa         : 162.91 GB
        Tpa         : 86.53 GB


        SET: DMS
        --------
        AcuI        : 17.89 GB
        AcuK        : 200.63 GB
        AcuN        : 294.56 GB
        DddA        : 73.16 GB
        DddC        : 108.97 GB
        DddD        : 4.19 GB
        DddL        : 0.25 GB
        DddP        : 4.20 GB
        DddQ        : 0.33 GB
        DddT        : 30.86 GB
        DddW        : 0.44 GB
        DddY        : 0.32 GB
        DdhA        : 2.46 GB
        DdhB        : 0.64 GB
        DdhC        : 0.29 GB
        DdhD        : 0.29 GB
        DmdA        : 5.11 GB
        DmdB        : 112.81 GB
        DmdD        : 5.24 GB
        DmoA        : 59.96 GB
        DmsA        : 115.08 GB
        DmsB        : 20.31 GB
        DmsC        : 2.60 GB
        DmsD        : 2.31 GB
        DorA        : 0.29 GB
        DorC        : 1.06 GB
        DorD        : 0.28 GB
        DsoA        : 0.27 GB
        DsoB        : 7.86 GB
        DsoC        : 1.56 GB
        DsoD        : 17.85 GB
        DsoE        : 0.99 GB
        DsoF        : 15.65 GB
        MarB        : 1.00 GB
        MarD        : 1.17 GB
        MarH        : 1.02 GB
        MarK        : 1.59 GB
        MddA        : 5.11 GB
        MddH        : 2.98 GB
        MsmA        : 3.87 GB
        MsmB        : 1.37 GB
        MsmC        : 0.70 GB
        MsmD        : 3.80 GB
        MsuC        : 58.78 GB
        MsuD        : 29.14 GB
        MsuE        : 5.20 GB
        Mtox        : 10.76 GB
        SnfG        : 16.97 GB


        SET: Sor_Soe
        ------------
        SoeA        : 103.17 GB
        SoeB        : 27.15 GB
        SoeC        : 35.04 GB
        SorA        : 5.48 GB
        SorB        : 1.46 GB


        SET: SQR
        --------
        CstA        : 66.29 GB
        CstB        : 2.46 GB
        SQRI        : 52.44 GB
        SQRII       : 57.44 GB
        SQRIII     : 176.69 GB
        SQRIV       : 28.31 GB
        SQRV        : 539.77 GB
        SQRVI       : 35.51 GB


        SET: SQ_SQDG
        ------------
        SftD        : 0.42 GB
        SftI        : 26.44 GB
        SftT        : 3.58 GB
        SftX        : 1.19 GB
        Sgdh        : 3.67 GB
        SmoB        : 10.29 GB
        SmoC        : 13.10 GB
        SqdA        : 0.22 GB
        SqdB        : 16.14 GB
        SqdC        : 4.12 GB
        SqdX        : 12.05 GB
        Sqald       : 1.35 GB
        Sqdh        : 1.79 GB
        SqgA        : 27.61 GB
        SqiA        : 1.03 GB
        SqiK        : 1.23 GB
        Sql         : 1.79 GB
        SqoD        : 1.23 GB
        SqvB        : 1.69 GB
        SqwD        : 1.51 GB
        SqwF        : 2.76 GB
        SqwG        : 4.58 GB
        SqwH        : 5.13 GB
        SqwI        : 1.67 GB
        SqwK        : 1.00 GB
        SqwL        : 0.84 GB
        YihQ        : 38.88 GB
        YihR        : 2.52 GB
        YihS        : 1.89 GB
        YihT        : 1.70 GB
        YihU        : 0.77 GB
        YihV        : 1.72 GB


        SET: sHdr_Dsr_Sox
        -----------------
        DoxA        : 0.29 GB
        DoxD        : 0.29 GB
        DsrE        : 10.34 GB
        DsrE3A      : 21.29 GB
        DsrE3B      : 0.95 GB
        DsrE3C      : 4.92 GB
        DsrF        : 9.80 GB
        DsrH        : 5.90 GB
        DsrL        : 34.79 GB
        DsrR        : 2.06 GB
        DsrS        : 1.60 GB
        EMO         : 43.40 GB
        FccA        : 3.75 GB
        FccB        : 10.52 GB
        LipS1       : 16.61 GB
        LipS2       : 16.91 GB
        LipT        : 12.86 GB
        LbpA1       : 1.07 GB
        LbpA2       : 7.61 GB
        Rhd442      : 0.90 GB
        SOR         : 1.47 GB
        SoxA        : 16.65 GB
        SoxB        : 47.30 GB
        SoxC        : 28.96 GB
        SoxD        : 16.64 GB
        SoxE        : 3.96 GB
        SoxF        : 14.71 GB
        SoxG        : 5.73 GB
        SoxH        : 8.91 GB
        SoxR        : 3.80 GB
        SoxS        : 4.40 GB
        SoxT1       : 16.11 GB
        SoxT2       : 15.13 GB
        SoxV        : 9.16 GB
        SoxW        : 7.06 GB
        SoxX        : 7.12 GB
        SoxY        : 6.79 GB
        SoxZ        : 7.09 GB
        oxDsrA      : 20.74 GB
        oxDsrB      : 16.19 GB
        oxDsrC      : 10.55 GB
        oxDsrJ      : 2.90 GB
        oxDsrK      : 30.78 GB
        oxDsrM      : 12.78 GB
        oxDsrN      : 22.31 GB
        oxDsrO      : 10.95 GB
        oxDsrP      : 18.26 GB
        redDsrA     : 20.40 GB
        redDsrB     : 19.39 GB
        redDsrC     : 4.73 GB
        redDsrD     : 2.17 GB
        redDsrE     : 2.16 GB
        redDsrF     : 1.86 GB
        redDsrH     : 1.97 GB
        redDsrJ     : 6.80 GB
        redDsrK     : 29.03 GB
        redDsrM     : 17.91 GB
        redDsrN     : 17.01 GB
        redDsrO     : 13.17 GB
        redDsrP     : 24.67 GB
        sHdrA       : 6.90 GB
        sHdrB1      : 7.08 GB
        sHdrB2      : 4.65 GB
        sHdrB3      : 3.20 GB
        sHdrC1      : 4.66 GB
        sHdrC2      : 3.76 GB
        sHdrH       : 5.13 GB
        sHdrI       : 1.53 GB
        sHdrT       : 9.40 GB
        sLplAB      : 26.01 GB
        """)


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
    inputdef = parser.add_argument_group("Input data")
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
        dest="verbose",
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
    parameters = parser.add_argument_group("Annotation parameters")
    parameters.add_argument(
        "-n",
        dest="name",
        type=str,
        default="project",
        metavar="<string>",
        help="Name new project" if show_all else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--hmm-sets",
        nargs="+",
        dest="HMM_sets",
        type=str,
        default=sorted(["DHPS", "DMS", "Dsr", "SQ", "Aryl"]),
        choices=sorted(["DHPS", "DMS", "Dsr", "SQ", "Aryl"]),
        metavar="",
        help="Limit to HMM sets (whitespace separated)"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--hmm-packages",
        nargs="+",
        dest="HMM_packages",
        type=str,
        default=sorted(["v8"]),
        choices=sorted(["v7", "v8", "chen", "disco", "hmss2"]),
        metavar="",
        help="Limit to specific HMM packages (whitespace separated)"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--cut-type",
        dest="threshold_type",
        type=int,
        default=1,
        metavar="<int>",
        choices=[1, 2, 3],
        help="Choice of cutoff: 1 optimized; 2 trusted; 3 noise"
        if never_show
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--cut-score",
        dest="thrs_score",
        type=int,
        default=50,
        metavar="<int>",
        help="Global minimal score cutoff" if show_all else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--threshold-factor",
        dest="threshold_factor",
        type=float,
        default=1.0,
        metavar="<float>",
        help=(
            "Multiplicative factor applied to score thresholds "
            "(>= 0.0, default: 1.0)"
            if show_all
            else argparse.SUPPRESS
        ),
    )

    parameters.add_argument(
        "--refseq-ident",
        dest="refseq_identity",
        type=int,
        default=90,
        metavar="<int>",
        help="Minimal percent identity to reference sequence set cross check"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--blast-cross-check",
        dest="bool_cross_check",
        action="store_false",
        help="Use Diamond blastp cross check for hit selection"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--taxonomy",
        dest="taxonomy_file",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help="Add taxonomy from this tab separated taxonomy file",
    )
    parameters.add_argument(
        "-s",
        dest="stage",
        type=int,
        default=0,
        choices=[0, 1, 2, 3, 4, 5],
        help="Start at step" if never_show else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "-x",
        dest="exit",
        type=int,
        default=10,
        choices=[0, 1, 2, 3, 4, 5],
        help="Exit at step" if never_show else argparse.SUPPRESS,
    )

    parameters.add_argument(
        "--max-seq-per-genome",
        dest="max_seqs_per_genome",
        type=int,
        default=4,
        help="Max. number of paralogs per genome forwarded to cross check via Diamond"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--diamond-speed",
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
    parameters.add_argument(
        "--disable-reports",
        dest="disable_individual_reports",
        action="store_true",
        help="Disable individual reports, only bulk database"
        if show_all
        else argparse.SUPPRESS,
    )
    parameters.add_argument(
        "--optimized-cutoff-cross-check",
        dest="optimized_cutoff_cross_check",
        action="store_true",
        help="Use optimized cutoff for hit selection"
        if never_show
        else argparse.SUPPRESS,
    )

    # Synteny options
    # synteny = parser.add_argument_group("Synteny options")

    # Information on resources
    information = parser.add_argument_group("Information on resources")
    information.add_argument(
        "--stat-keywords",
        action="store_true",
        help="Print patterns for keyword naming" if never_show else argparse.SUPPRESS,
    )
    information.add_argument(
        "--stat-csb",
        action="store_true",
        help="Print automatically found csbs" if never_show else argparse.SUPPRESS,
    )
    information.add_argument(
        "--stat-genomes",
        action="store_true",
        help="Print taxonomy information from database"
        if never_show
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
        help="Max. nucleotide distance to be considered syntenic genes"
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
        "--disable-synteny-completion",
        dest="disable_synteny_completion",
        action="store_true",
        help="Disable syntenic block supported annotation"
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
    readmap = parser.add_argument_group("Read mapping (optional)")

    readmap.add_argument(
        "--read-mapping",
        dest="use_read_mapping",
        action="store_true",
        help="Enable read-mapping analysis for files located in the -f directory",
    )

    readmap.add_argument(
        "--help-read-mapping",
        action=ReadMappingHelpAction,
        nargs=0,
        help="Explain read-mapping workflow and list available GPKG packages (set + estimated RAM) and exit",
    )

    readmap.add_argument(
        "--gpkg-sets",
        nargs="+",
        dest="gpkg_sets",
        type=str,
        choices=sorted(
            [
                "Apr",
                "Qmo",
                "Asr",
                "Mcc",
                "Phs",
                "Ttr",
                "CS",
                "Aryl",
                "DHPS",
                "Taurine",
                "Isethionat",
                "DMS",
                "sHdr",
                "Dsr",
                "Sox",
                "Shy",
                "Sud",
                "Sor",
                "Soe",
                "SQ",
                "SQDG",
                "SQR",
            ]
        ),
        metavar="<SET>",
        default=[],
        help=(
            "Select GPKG package sets (whitespace separated)"
            if show_all
            else argparse.SUPPRESS
        ),
    )

    readmap.add_argument(
        "--gpkg-packs",
        nargs="+",
        dest="gpkg_packs",
        type=str,
        default=[],
        metavar="",
        help=(
            "Select specific GPKG packages (whitespace separated)"
            if show_all
            else argparse.SUPPRESS
        ),
    )

    readmap.add_argument(
        "--rm-threads",
        dest="rm_threads",
        type=int,
        default=14,
        metavar="<int>",
        help="Number of threads to use for read mapping"
        if never_show
        else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-evalue",
        dest="rm_evalue",
        type=float,
        default=1e-5,
        metavar="<float>",
        help="E-value threshold for homology search" if show_all else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-placements-cutoff",
        dest="rm_placements_cutoff",
        type=float,
        default=0.75,
        metavar="<float>",
        help="Placement cutoff for phylogenetic placement."
        if show_all
        else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-resolve-placements",
        dest="rm_resolve_placements",
        action="store_true",
        help="Resolve ambiguous phylogenetic placements"
        if show_all
        else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-min-orf-length",
        dest="rm_min_orf_length",
        type=int,
        default=96,
        metavar="<int>",
        help="Minimum ORF length" if show_all else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-restrict-read-length",
        dest="rm_restrict_read_length",
        type=int,
        default=None,
        metavar="<int>",
        help="Maximum read length" if show_all else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--rm-translation-table",
        dest="rm_translation_table",
        type=int,
        default=11,
        metavar="<int>",
        help="NCBI translation table to use" if show_all else argparse.SUPPRESS,
    )
    readmap.add_argument(
        "--ram-limit-min",
        dest="rm_ram_limit_min",
        type=float,
        default=0.0,
        metavar="<float>",
        help=("Use gpkg with at least minimum estimated RAM")
        if show_all
        else argparse.SUPPRESS,
    )

    readmap.add_argument(
        "--ram-limit-max",
        dest="rm_ram_limit_max",
        type=float,
        default=16.0,
        metavar="<float>",
        help=("Use gpkg with at less than estimated RAM. 0 means unlimited")
        if show_all
        else argparse.SUPPRESS,
    )

    # --- Interleaved FASTQ ---
    readmap.add_argument(
        "--interleaved",
        action="store_true",
        default=False,
        help=(
            "Treat input FASTQ files as interleaved reads "
            "(mutually exclusive with reverse pairs)."
        )
        if show_all
        else argparse.SUPPRESS,
    )

    # Work step regulation
    flow = parser.add_argument_group("Work step regulation")
    flow.add_argument(
        "-no_remove_intermediate",
        dest="use_remove_unassigned_intermediates",
        action="store_false",
        help="Remove intermediate hits without genetic context"
        if never_show
        else argparse.SUPPRESS,
    )
    flow.add_argument(
        "-no_remove_exclusion_singletons",
        dest="use_remove_exclusion_singletons",
        action="store_false",
        help="Remove genes that not occur as singletons"
        if never_show
        else argparse.SUPPRESS,
    )

    # Output operators / fetch
    operators = parser.add_argument_group("Output filtering and export")
    operators.add_argument(
        "--help-fetch",
        action=FetchHelpAction,
        nargs=0,
        help="Explain how -fd and -fc are interpreted (AND/OR logic, brackets, combinations) and exit",
    )
    operators.add_argument(
        "--help-output",
        action=FetchOutputHelpAction,
        nargs=0,
        help="Explain fetch output files (summary tables, FASTA outputs, graphs) and exit",
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
        help="Taxonomy level [Superkingdom,Phylum,Class,Ordnung,Family,Genus,Species]"
        if show_all
        else argparse.SUPPRESS,
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
        metavar="",
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
        metavar="",
        help="Select only proteins with these domains (whitespace separated)",
    )
    operators.add_argument(
        "-fc",
        nargs="+",
        dest="fetch_csbs",
        type=str,
        default=[],
        metavar="",
        help="Select gene clusters encoding the given proteins (whitespace separated). The ':' without whitespace will be interpreted as logical OR",
    )
    operators.add_argument(
        "--allow-fd-adds-genomes",
        dest="fd_can_add_genomes",
        action="store_true",
        help="When -fc and -fd are combined, allow -fd to add genomes beyond the -fc selection.",
    )
    operators.add_argument(
        "-fnd",
        nargs="+",
        dest="exclude_domains",
        type=str,
        default=[],
        metavar="",
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
        if never_show
        else argparse.SUPPRESS,
    )
    operators.add_argument(
        "-kc",
        dest="keywords_connector",
        type=str,
        default="OR",
        choices=["AND", "OR"],
        help="Select cluster with keywords connected by AND or OR"
        if never_show
        else argparse.SUPPRESS,
    )
    operators.add_argument(
        "--disable-filters",
        dest="use_non_valid_hits",
        action="store_true",
        help="Fetch without noise hit filters",
    )
    operators.add_argument(
        "--fasta",
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
        help="Print graphs for the selected output" if show_all else argparse.SUPPRESS,
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
        )
        if show_all
        else argparse.SUPPRESS,
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
        gpkg=d["GPKG_DIR"],
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
        threshold_factor=float(getattr(ns, "threshold_factor", 1.0)),
        taxonomy_file=_s(ns, "taxonomy_file"),
        refseq_identity=int(getattr(ns, "refseq_identity", 90)),
        name=getattr(ns, "name", "project"),
        stage=int(getattr(ns, "stage", 0)),
        exit=int(getattr(ns, "exit", 10)),
    )

    cli_resources = CliResources(
        HMM_sets=list(getattr(ns, "HMM_sets", [])),
        HMM_packages=list(getattr(ns, "HMM_packages", [])),
        disable_individual_reports=bool(
            getattr(ns, "disable_individual_reports", False)
        ),
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
        gpkg_sets=list(getattr(ns, "gpkg_sets", [])),
        gpkg_packs=list(getattr(ns, "gpkg_packs", [])),
        threads=ns.rm_threads,
        evalue=ns.rm_evalue,
        placements_cutoff=ns.rm_placements_cutoff,
        resolve_placements=bool(getattr(ns, "rm_resolve_placements", False)),
        min_orf_length=ns.rm_min_orf_length,
        restrict_read_length=ns.rm_restrict_read_length,
        translation_table=ns.rm_translation_table,
        ram_limit_max=ns.rm_ram_limit_max,
        ram_limit_min=ns.rm_ram_limit_min,
        interleaved=bool(getattr(ns, "interleaved", False)),
        ram_profile_file=ns.ram_profile_file or paths_cfg.gpkg,
    )

    cli_flow = CliFlow(
        redo_taxonomy=bool(getattr(ns, "redo_taxonomy", False)),
        disable_synteny_completion=bool(getattr(ns, "use_synteny_completion", True)),
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
        fd_can_add_genomes=ns.fd_can_add_genomes,
        keywords_connector=getattr(ns, "keywords_connector", "OR"),
        print_fasta=getattr(ns, "print_fasta", False),
        print_graphs=bool(getattr(ns, "print_graphs", False)),
        use_non_valid_hits=bool(getattr(ns, "use_non_valid_hits", True)),
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
    _set_default("ram_profile_file", str(paths.SRC_FILE_GPKG_RAM_INFORMATION))

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

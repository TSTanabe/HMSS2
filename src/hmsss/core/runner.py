# src/hmsss/core/runner.py
from __future__ import annotations

import os

from hmsss.db import project as project
from hmsss.db.database import index_database

from hmsss.stages.ressource_prep import ressource_preparation
from hmsss.stages.fasta_preparation import fasta_preparation
from hmsss.stages.initial_search import initial_search
from hmsss.stages.cross_check import reference_sequence_check
from hmsss.stages.parse_reports import parse_reports_to_database
from hmsss.stages.csb import csb_finder
from hmsss.stages.taxonomy import collect_taxonomy_information
from hmsss.stages.output_dataset import output_operator, output_statistics
from hmsss.stages.process_seqfiles import process_operator

from hmsss.io.queue import queue_protein_annotation_inputs
from hmsss.io.output import print_file_content

from hmsss.core.logging import setup_logging, print_header, get_logger

logger = get_logger(__name__)

"""
Pipeline runner for HMSSS.

Defines `run_pipeline`, which orchestrates all pipeline stages
(resource preparation, FASTA preprocessing, hmmsearch, cross-check,
report parsing, CSB detection, taxonomy, and post-processing).

Stages:
    0: Resource preparation (always executed if stage < 100)
    1: FASTA/Prodigal preparation
    2: Hmmsearch (including queueing)
    3: Cross-check / cutoff promotion
    4: Parse reports into database
    5: Collinear syntenic block (CSB) detection
    6: Taxonomy information collection

Special stages:
    100: Taxonomy-only mode
    101: Database fetch, output, and processing operators
"""


def run_pipeline(config) -> None:
    """Execute the full HMSSS pipeline according to configuration.

    The pipeline proceeds from `config.stage` to `config.exit`, executing
    the appropriate modules. Special stages 100 and 101 are handled for
    taxonomy-only and fetch/processing modes, respectively.

    Args:
        config: A validated `Config` object containing all CLI parameters,
            paths, and runtime state.

    Side Effects:
        - Creates result directories and initializes logging.
        - Populates the SQLite database with parsed results.
        - Produces output datasets, statistics, and processed files.

    Raises:
        RuntimeError: If database is missing for fetch operations.
    """

    # Ergebnisraum anlegen (Unterordner, Projektpfade, etc.)
    project.prepare_result_space(config)

    # Logging initialisieren
    log_file = os.path.join(config.result_files_directory, "execution_logfile.txt")
    setup_logging(getattr(config, "verbose", 1), log_file)

    if config.stage < 6:
        print_header("Initializing resources")
        ressource_preparation(config)

    # --- Stage 1: FASTA/Prodigal ---
    if config.stage <= 1 <= config.exit:
        print_header("Prokaryotic gene recognition and translation (prodigal)")
        fasta_preparation(config)

    # --- Stage 2: Hmmsearch ---
    if config.stage <= 2 <= config.exit:
        print_header("Queueing input files")
        queue_protein_annotation_inputs(config)

        print_header("Searching for homologous sequences (hmmsearch)")
        initial_search(config)
        config.stage = 2

    # --- Stage 3: Cross-Check / Cutoffs ---
    if config.stage <= 3 <= config.exit:
        print_header("Cross check with reference sequences / cutoff optimization")
        reference_sequence_check(config)
        config.stage = 3

    # --- Stage 4: Reports -> DB ---
    if config.stage <= 4 <= config.exit:
        print_header("Parse trusted hits and recognized gene clusters into database")
        queue_protein_annotation_inputs(config)
        parse_reports_to_database(config)
        config.stage = 4

    # --- Stage 5: CSB Finder ---
    if config.stage <= 5 <= config.exit:
        print_header("Searching for collinear syntenic blocks (CSB)")
        csb_finder(config)
        config.stage = 5

    # --- Stage 6: Taxonomy ---
    if config.stage <= 6 <= config.exit:
        print_header("Assigning taxonomy information")
        collect_taxonomy_information(config)

    if config.stage == 100:
        print_header("Assigning taxonomy information (stage 100)")
        collect_taxonomy_information(config)

    # --- Output-/Stats-/Processing-Operatoren ---
    if config.stage == 101:
        print_header("Output from database (fetch)")
        if not config.database_directory:
            logger.error(
                f"Database not found in given project {config.result_files_directory}. Please use a valid project directory or use the -db argument to provide a valid database for fetch operations."
            )
        index_database(config.database_directory)
        output_operator(config)

    if getattr(config, "stat_genomes", False):
        output_statistics(config)

    if getattr(config, "process", False):
        process_operator(config)

    if getattr(config, "stat_keywords", False):
        print_file_content(config.patterns_file)

    if getattr(config, "stat_csb", False):
        print_file_content(config.csb_output_file)


__all__ = ["run_pipeline"]

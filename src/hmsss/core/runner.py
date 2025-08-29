# src/hmsss/core/runner.py
from __future__ import annotations

import os


from hmsss.db import project as project

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
from hmsss.db.database import index_database
from hmsss.io.output import print_file_content
from hmsss.core.logging import setup_logging, print_header, get_logger
logger = get_logger(__name__)

def run_pipeline(config) -> None:
    """
    Orchestriert die komplette HMSSS-Pipeline anhand von 'stage' / 'exit'
    und den Mode-Flags (fetch / limiter / process).

    Erwartete Stages:
      0: Ressourcen vorbereiten (immer, wenn < 100)
      1: FASTA/Prodigal-Vorbereitung
      2: Hmmsearch (inkl. Queueing)
      3: Cross-Check / Cutoff-Promotion
      4: Reports parsen -> DB
      5: CSB Finder
      6: Taxonomy-Informationen sammeln

    Spezialstages:
    100: Nur Taxonomie (wenn -taxonomy_info & -db & kein -f)
    101: Operatoren/Fetch/Processing (wenn any_process_args_provided True)
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
            logger.error(f"Database not found in given project {config.result_files_directory}. Please use a valid project directory or use the -db argument to provide a valid database for fetch operations.")
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

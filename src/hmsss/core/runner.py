# src/hmsss/core/runner.py
from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from hmsss.core.logging import setup_logging, print_header
from hmsss.db import project as project

# Stages (bitte sicherstellen, dass du diese Module/Symbole erzeugst)
from hmsss.stages.ressource_prep import ressource_preparation
from hmsss.stages.fasta_preparation import fasta_preparation
from hmsss.stages.initial_search import initial_search
from hmsss.stages.cross_check import reference_sequence_check
from hmsss.stages.parse_reports import parse_reports_to_database
from hmsss.stages.csb import csb_finder
from hmsss.stages.taxonomy import collect_taxonomy_information
from hmsss.stages.output_dataset import output_operator, output_statistics
from hmsss.stages.process_seqfiles import process_operator

# Sonstige Aufrufer
from hmsss.io.queue import queue_files
from hmsss.db.database import index_database
from hmsss.io.output import print_file_content
from hmsss.io.queue import decompress_gz_recursively

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss


def _ensure_results_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def run_pipeline(options: Hmsss) -> None:
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
    _ensure_results_dir(options.result_files_directory)

    # --- Vorbereitende Schritte & Ressourcen ---
    if options.stage < 100:
        # Ergebnisraum anlegen (Unterordner, Projektpfade, etc.)
        project.prepare_result_space(options)

        # Logging initialisieren (falls nicht schon in __main__ geschehen)
        log_file = os.path.join(options.result_files_directory, "execution_logfile.txt")
        setup_logging(getattr(options, "verbose", 1), log_file)

        print_header("\nInitializing resources")
        ressource_preparation(options)
        decompress_gz_recursively(options.fasta_file_directory)

    # --- Stage 1: FASTA/Prodigal ---
    if options.stage <= 1 <= options.exit:
        print_header("\nProkaryotic gene recognition and translation (prodigal)")
        fasta_preparation(options)

    # --- Stage 2: Hmmsearch ---
    if options.stage <= 2 <= options.exit:
        print_header("\nQueueing input files")
        queue_files(options)

        print_header("\nSearching for homologous sequences (hmmsearch)")
        initial_search(options)
        options.stage = 2

    # --- Stage 3: Cross-Check / Cutoffs ---
    if options.stage <= 3 <= options.exit:
        print_header("\nCross check with reference sequences / cutoff optimization")
        reference_sequence_check(options)
        options.stage = 3

    # --- Stage 4: Reports -> DB ---
    if options.stage <= 4 <= options.exit:
        print_header("\nParse trusted hits and recognized gene clusters into database")
        parse_reports_to_database(options)
        options.stage = 4

    # --- Stage 5: CSB Finder ---
    if options.stage <= 5 <= options.exit:
        print_header("\nSearching for collinear syntenic blocks (CSB)")
        csb_finder(options)
        options.stage = 5

    # --- Stage 6: Taxonomy ---
    if options.stage <= 6 <= options.exit:
        print_header("\nAssigning taxonomy information")
        collect_taxonomy_information(options)

    # --- Spezialfälle >= 100 ---
    if options.stage > 99:
        # eigenes Logging sicherstellen (falls 100/101 direkt aufgerufen werden)
        log_file = os.path.join(options.result_files_directory, "execution_logfile.txt")
        setup_logging(getattr(options, "verbose", 1), log_file)

    if options.stage == 100:
        print_header("\nAssigning taxonomy information (stage 100)")
        collect_taxonomy_information(options)

    # --- Output-/Stats-/Processing-Operatoren ---
    if getattr(options, "fetch", False):
        print_header("\nOutput from database (fetch)")
        if not options.database_directory:
            raise SystemExit(
                "Please use the -db argument to provide a valid database for fetch operations."
            )
        index_database(options.database_directory)
        output_operator(options)

    if getattr(options, "stat_genomes", False):
        output_statistics(options)

    if getattr(options, "process", False):
        process_operator(options)

    if getattr(options, "stat_keywords", False):
        print_file_content(options.patterns_file)

    if getattr(options, "stat_csb", False):
        print_file_content(options.csb_output_file)


__all__ = ["run_pipeline"]

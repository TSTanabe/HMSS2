# src/hmsss/core/runner.py
from __future__ import annotations

import os
import sys
from pathlib import Path

from hmsss.core import project as project
from hmsss.core import ressource_prep
from hmsss.core.logging import get_logger, print_header, setup_logging
from hmsss.cross_check import cross_check
from hmsss.db import database
from hmsss.fasta_preparation import fasta_preparation
from hmsss.graft import initial_read_mapping
from hmsss.io import (
    db_fetch_genome_reports,
    metabolic_pathway_sankey,
    output_dataset,
    print_command_args,
)
from hmsss.read_mapping import initial_read_detection, prepare_libraries
from hmsss.search import initial_search
from hmsss.stages import (
    process_seqfiles,  # import process_operator
    taxonomy,  # import collect_taxonomy_information
)

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

    if config.stage < 6 or 40 >= config.stage >= 30:
        print_header("Initializing resources")
        ressource_prep.ressource_preparation(config)

    # --- Stage 1: FASTA/Prodigal ---
    if config.stage <= 1 <= config.exit:
        print_header("Prokaryotic gene recognition and translation (prodigal)")
        fasta_preparation.fasta_preparation(config)

    # --- Stage 2: Hmmsearch ---
    # 2 und teile von 4 zusammen gelegt
    if config.stage <= 2 <= config.exit:
        print_header("Searching for homologous sequences (hmmsearch)")
        initial_search.initial_search(config)
        config.stage = 2

    # --- Stage 3: Cross-Check / Cutoffs ---
    if config.stage <= 3 <= config.exit:
        print_header("Cross check with reference sequences / cutoff optimization")
        if config.bool_cross_check:
            cross_check.reference_sequence_check(config)
        database.rebuild_presence_tables(config.database_directory)
        database.index_database(config.database_directory)
        config.stage = 3

    # --- Stage 6: Taxonomy ---
    if config.stage <= 6 <= config.exit:
        print_header("Assigning taxonomy information")
        taxonomy.collect_taxonomy_information(config)

    if config.stage <= 7 <= config.exit:
        print_header("Writing individual hit reports and metabolic pathway report")
        db_fetch_genome_reports.write_individual_genome_reports(config)
        print_header("Creating metabolic pathway plots")
        metabolic_pathway_sankey.create_metabolic_pathway_plots(config)

    # The 30 steps are for read mapping with the custom implementation from here

    if config.stage == 30:
        print_header("Library preparation for metagenome search")
        result = prepare_libraries.collect_placement_search_hmms(config)
        prepare_libraries.create_placement_search_hmm_library(config, result)
        config.stage = 31
    print(config.stage)
    if config.stage <= 31 <= config.exit:
        print_header("Initial read detection in metagenomes")
        initial_read_detection.initial_read_detection(
            config,
            min_orf_length=96,
            translation_table=11,
            batch_size=10000,
            evalue=1e-5,
        )
        config.stage = 32

    if config.stage <= 32 <= config.exit:
        print_header("Phylogenetic placement of detected reads")
        cross_check.reference_sequence_check(config)
        # TODO Es ist nochunklar ob das platzieren so funktioniert gemessen am ressourcen verbrauch, code auch noch ungetestet
        # die weitere folge ist das alle hits noch nicht als valid gelten, bzw könnte man auch alle einfach als valid annehmen,
        # da die platzierung ja filtert
        config.stage = 33

    if config.stage <= 33 <= config.exit:
        print_header("Writing summary reports")
        config.stage = 34

    # The 50 steps are implementation of GraftM python module
    if config.stage == 50:
        print_header("Mapping reads with GraftM")
        ressource_prep.ressource_preparation(config)
        initial_read_mapping.initial_read_mapping(config)
        # TODO Ungetestet ist bisher das parsing des outputs. Das muss noch erstellt werden
        #
        sys.exit()

    if config.stage == 100:
        print_header("Assigning taxonomy information")
        taxonomy.collect_taxonomy_information(config)

    # --- Output-/Stats-/Processing-Operatoren ---
    if config.stage == 101:
        print_header("Output from database (fetch)")
        if not config.database_directory:
            logger.error(
                f"Database not found in given project {config.result_files_directory}. Please use a valid project directory or use the -db argument to provide a valid database for fetch operations."
            )
        database.rebuild_presence_tables(config.database_directory)
        output_dataset.output_operator(config)

    if getattr(config, "stat_genomes", False):
        output_dataset.output_statistics(config)

    if getattr(config, "process", False):
        process_seqfiles.process_operator(config)

    if getattr(config, "stat_keywords", False):
        print_command_args.print_file_content(config.patterns_file)

    if getattr(config, "stat_csb", False):
        print_command_args.print_file_content(config.csb_output_file)


__all__ = ["run_pipeline"]

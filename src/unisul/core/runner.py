# src/unisul/core/runner.py
from __future__ import annotations

import os
import sys

from unisul.core import project as project
from unisul.core import ressource_prep
from unisul.core.logging import get_logger, print_header, setup_logging
from unisul.data_management import require_pipeline_data
from unisul.db import database
from unisul.fasta_preparation import fasta_preparation
from unisul.graft import initial_read_mapping
from unisul.io import (
    db_fetch_genome_reports,
    metabolic_pathway_sankey,
    output_dataset,
    print_command_args,
)
from unisul.placement import (
    initial_read_detection,
    mmseqs_reference_assignment,
    phylogenetic_placement,
    prepare_library,
)
from unisul.search import initial_search
from unisul.stages import (
    process_seqfiles,  # import process_operator
    taxonomy,  # import collect_taxonomy_information
)

logger = get_logger(__name__)

"""
Pipeline runner for UNISUL.

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
    """Execute the full UNISUL pipeline according to configuration.

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

    require_pipeline_data(config)

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
        config.exit = 10

    # --- Stage 2: Hmmsearch ---
    if config.stage <= 2 <= config.exit:
        print_header("Searching for homologous sequences")
        initial_search.initial_search(config)
        config.stage = 3

    # --- Stage 3: MMseqs2 competitive reference assignment ---
    if config.stage <= 3 <= config.exit:
        if config.disable_mmseqs_reference_assignment:
            print_header("Comparison with Reference Sequences")
            logger.info(
                "MMseqs2 reference classification disabled; "
                "intermediate-hit FASTA files remain unfiltered"
            )
        else:
            print_header("Comparison with Reference Sequences")
            mmseqs_reference_assignment.mmseqs_placement_package_reference_assignment(
                config=config
            )
        if config.placements:
            print_header("Phylogenetic placement of genome hits")
            phylogenetic_placement.run_placement_stage(config)

        config.stage = 4

    # --- Stage 4: Database indexing ---
    if config.stage <= 4 <= config.exit:
        print_header("Database index generation")

        # if config.bool_cross_check:
        #    cross_check.reference_sequence_check(config)

        database.rebuild_presence_tables(config.database_directory)
        database.index_database(config.database_directory)

        config.stage = 4

    # --- Stage 6: Taxonomy ---
    if config.stage <= 6 <= config.exit:
        print_header("Assigning taxonomy information")
        # taxonomy.collect_taxonomy_information(config)
        database.insert_taxonomy_data(config.database_directory, config.taxonomy_file)

    if config.stage <= 7 <= config.exit:
        print_header("Writing individual hit reports and metabolic pathway report")
        db_fetch_genome_reports.write_individual_genome_reports(config)
        print_header("Creating metabolic pathway plots")
        metabolic_pathway_sankey.create_metabolic_pathway_plots(config)

    # ============================================================
    # METAGENOME ORF / PLACEMENT WORKFLOW
    # ============================================================

    # Stage 31: Build HMM search library
    if config.stage <= 31 <= config.exit:
        print_header("Library preparation for metagenome search")
        placement_packages = prepare_library.collect_placement_search_hmms(config)
        prepare_library.create_placement_search_hmm_library(
            config,
            placement_packages,
        )
        config.stage = 31
        config.exit = 40

    # Stage 32: Detect metagenomic ORFs and HMM candidates
    if config.stage <= 32 <= config.exit:
        print_header("Initial read detection in metagenomes")
        initial_read_detection.initial_read_detection(
            config,
            min_orf_length=config.min_orf_length,
            translation_table=config.translation_table,
            batch_size=config.batch_size,
            evalue=config.evalue,
        )
        config.stage = 32

    # Stage 33: MMseqs reference classification / Outgroup filtering
    if config.stage <= 33 <= config.exit:
        print_header("MMseqs2 reference classification and Outgroup filtering")
        if config.disable_mmseqs_reference_assignment:
            logger.info(
                "MMseqs2 reference classification disabled; "
                "intermediate-hit FASTA files remain unfiltered"
            )
        else:
            mmseqs_reference_assignment.mmseqs_placement_package_reference_assignment(
                config=config
            )

        config.stage = 33

    # Stage 34: EPA-ng phylogenetic placement
    if config.stage <= 34 <= config.exit:
        print_header("Phylogenetic placement of detected reads")
        phylogenetic_placement.run_placement_stage(config)
        config.stage = 34

    # Stage 35: Standard per-genome/metagenome reports
    if config.stage <= 35 <= config.exit:
        print_header("Writing metagenome hit reports")
        db_fetch_genome_reports.write_individual_genome_reports(config)
        config.stage = 35

    # The 50 steps are implementation of GraftM python module
    if config.stage == 50:
        print_header("Mapping reads with GraftM")
        ressource_prep.ressource_preparation(config)
        initial_read_mapping.initial_read_mapping(config)
        # Output bleibt der graftM output
        #
        sys.exit()

    if config.stage == 100:
        print_header("Assigning taxonomy information")
        taxonomy.collect_taxonomy_information(config)
        database.insert_taxonomy_data(config.database_directory, config.taxonomy_file)

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

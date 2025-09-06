from __future__ import annotations

import os

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import search_hmmer as search_hmmer
from hmsss.algorithms import search_cross_reference as search_cross_reference
from hmsss.io.queue import queue_protein_annotation_inputs
from hmsss.db.database import create_database



log = get_logger(__name__)

"""
Initial search stage (hmmsearch) for HMSSS.

This stage ensures the project database exists, runs HMMER searches per genome
(when no global report is present), concatenates individual reports into a
global report, and partitions hits into trusted / intermediate / noise
according to score thresholds.
"""

def initial_search(config) -> None:
    """Run the initial hmmsearch stage and prepare global reports.

    Steps:
      1) Create the project SQLite database if missing.
      2) If a global report already exists, reuse it.
      3) Otherwise, run hmmsearch for each genome and collect per-genome reports.
      4) Concatenate per-genome reports into a global report.
      5) Split global hits into trusted / intermediate / noise via thresholds.

    Args:
        options: Pipeline options/config with fields such as
            `database_directory`, `glob_report`, `cores`, `hmmreport_files`,
            and `score_threshold_file`.

    Side Effects:
        - Creates database file and per-genome `.hmmreport` if needed.
        - Writes a concatenated global report.
        - Writes categorized hit files into the cross-check directory.

    Raises:
        FileNotFoundError: If expected input files are missing during concatenation.
    """

    print_header("Initial search (hmmsearch)", logger=log)
    if not os.path.isfile(config.database_directory):
        create_database(config.database_directory)

    # Falls bereits glob_report existiert, überspringen wir die Suche
    if config.glob_report and os.path.isfile(config.glob_report):
        log.info("Using existing global hmmreport: %s", config.glob_report)
    else:
        log.info("Running hmmsearch for each input genome")
        search_hmmer.consecutive_hmm_search(config, int(config.cores / 2))

        log.info("Concatenating hmmsearch results for cross reference check")
        config.glob_report = search_cross_reference.concatenate_hmmreports_cat_xargs(
            config.hmmreport_files, config.glob_report
        )

    # Globale Auswertung in trusted/noise/intermediate
    log.info("Filtering hits into trusted, noise, and intermediate categories")
    search_cross_reference.filter_trusted_and_noise_hits(
        config, config.glob_report, config.cores
    )

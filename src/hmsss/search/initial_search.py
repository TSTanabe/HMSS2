from __future__ import annotations

import os

from hmsss.core.logging import get_logger, print_header
from hmsss.core import queue
from hmsss.search import search_pyhmmer
from hmsss.db.database import create_database


logger = get_logger(__name__)

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
        config: Pipeline options/config with fields such as
            `database_directory`, `glob_report`, `cores`, `hmmreport_files`,
            and `score_threshold_file`.

    Side Effects:
        - Creates database file and per-genome `.hmmreport` if needed.
        - Writes a concatenated global report.
        - Writes categorized hit files into the cross-check directory.

    Raises:
        FileNotFoundError: If expected input files are missing during concatenation.
    """

    print_header("Initial search (hmmsearch)", logger=logger)
    if os.path.isfile(config.database_directory):

        queue.queue_protein_annotation_inputs(config)
        removed = queue.remove_genomes_already_in_db_from_queue(config)
        logger.info(
            f"Removed {removed} genomes from queue because they are already in the database."
        )
    else:
        create_database(config.database_directory)
        queue.queue_protein_annotation_inputs(config)

    if not config.queued_genomes:
        logger.info("No genomes left to process after DB-filtering. Done.")
        return
    logger.info("Running pyhmmer search + parsing + cluster detection + DB write")
    search_pyhmmer.consecutive_hmm_search(config, processes=int(config.cores / 2))

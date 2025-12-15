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
        If database is present genomeIDs in the database will be ignored.
        Uses pyhmmer to reduce IO operations and processes everything except the cross-check
        without writing extra external files
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
    search_pyhmmer.consecutive_hmm_search(config, processes=int(config.cores))

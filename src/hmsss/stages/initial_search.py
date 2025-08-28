from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import search_hmmer as search_hmmer
from hmsss.algorithms import search_cross_reference as search_cross_reference
from hmsss.io.queue import queue_protein_annotation_inputs
from hmsss.db.database import create_database

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def initial_search(options: Hmsss) -> None:
    """
    Aus __main__.py:
    - ggf. queueing
    - hmmsearch pro Genome
    - globale Report-Datei zusammenführen
    - trusted/noise/intermediate anhand Cutoffs vorsortieren
    """
    print_header("Initial search (hmmsearch)", logger=log)
    if not os.path.isfile(options.database_directory):
        create_database(options.database_directory)

    # Falls bereits glob_report existiert, überspringen wir die Suche
    if options.glob_report and os.path.isfile(options.glob_report):
        log.info("Using existing global hmmreport: %s", options.glob_report)
    else:


        log.info("Running hmmsearch for each input genome")
        search_hmmer.consecutive_hmm_search(options, int(options.cores / 2))

        log.info("Collecting hmmreports for each genome")
        queue_protein_annotation_inputs(options)

        log.info("Concatenating hmmsearch results for cross reference check")
        options.glob_report = search_cross_reference.concatenate_hmmreports_cat_xargs(
            options.hmmreport_files, options.glob_report
        )


    # Globale Auswertung in trusted/noise/intermediate
    log.info("Filtering hits into trusted, noise, and intermediate categories")
    search_cross_reference.filter_trusted_and_noise_hits(options, options.glob_report, options.cores)

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.io import queue as queue
from hmsss.algorithms import search as search

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

    # Falls bereits glob_report existiert, überspringen wir die Suche
    if options.glob_report and os.path.isfile(options.glob_report):
        log.info("Using existing global hmmreport: %s", options.glob_report)
    else:
        log.info("Queueing input files")
        queue.queue_files(options)

        log.info("Running hmmsearch for each input genome")
        genome_to_hmmreport = search.unified_search(options, int(options.cores / 2))

        log.info("Concatenating hmmsearch results")
        options.glob_report = search.concatenate_hmmreports_cat_xargs(
            genome_to_hmmreport, options.glob_report
        )

    # Globale Auswertung in trusted/noise/intermediate
    log.info("Filtering hits into trusted, noise, and intermediate categories")
    search.filter_trusted_and_noise_hits(options, options.glob_report, options.cores)

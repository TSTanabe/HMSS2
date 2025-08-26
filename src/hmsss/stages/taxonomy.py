from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.db import database as database

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def collect_taxonomy_information(options: Hmsss) -> None:
    """
    Aus __main__.py: Taxonomie-Datei (TSV) in DB einspielen, falls vorhanden.
    """
    print_header("Collect taxonomy information", logger=log)

    if options.taxonomy_file and os.path.isfile(options.taxonomy_file):
        log.info("Writing taxonomy assignments to database: %s", options.taxonomy_file)
        database.insert_taxonomy_data(options.database_directory, options.taxonomy_file)
    else:
        log.warning("Taxonomy file was not provided / not found")

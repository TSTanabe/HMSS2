from __future__ import annotations

import os

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.db import database as database


log = get_logger(__name__)


def collect_taxonomy_information(config: Config) -> None:
    """Insert taxonomy assignments into the database.

    If `options.taxonomy_file` exists, its content is written into the
    project database. Otherwise, a warning is logged.

    Args:
        config: Pipeline configuration with `.taxonomy_file`
            and `.database_directory`.
    """
    print_header("Collect taxonomy information", logger=log)

    if config.taxonomy_file and os.path.isfile(config.taxonomy_file):
        log.info("Writing taxonomy assignments to database: %s", config.taxonomy_file)
        database.insert_taxonomy_data(config.database_directory, config.taxonomy_file)
    else:
        log.warning("Taxonomy file was not provided / not found")

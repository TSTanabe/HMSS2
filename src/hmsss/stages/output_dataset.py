from __future__ import annotations

import os
from datetime import datetime
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.io import output as output
from hmsss.io import datasets as datasets

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def output_operator(options: Hmsss) -> None:
    """
    Aus __main__.py:
    - Ergebnisordner (zeitgestempelt) anlegen
    - CLI-Args protokollieren
    - FASTA & Tabellen exportieren
    - Binär-Dataset erzeugen
    """
    print_header("Output operator (fetch/export)", logger=log)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    directory = os.path.join(options.result_files_directory, f"fetch_{ts}")
    os.makedirs(directory, exist_ok=True)

    output.print_command_line_args(os.path.join(directory, "1_fetch_command.txt"))

    protein_dict, cluster_dict, taxon_dict = output.fetch_fasta_and_hit_data(options)
    output.print_fasta_and_hit_outputs(
        directory, protein_dict, cluster_dict, taxon_dict
    )

    datasets.main_binary_dataset(
        options, directory, protein_dict, cluster_dict, taxon_dict
    )
    log.info("Generated binary dataset → %s", directory)


def output_statistics(options: Hmsss) -> None:
    from hmsss.db import database as database

    database.fetch_genome_statistic(options.database_directory)

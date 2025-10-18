from __future__ import annotations

import os
from datetime import datetime
from typing import Tuple, Dict, Any

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.io import print_reports, db_fetch_taxonomy, output as output, db_fetch_protein
from hmsss.db import database as database

logger = get_logger(__name__)

"""
Stage: Output operator.

Fetch/export pipeline results to a timestamped directory:
- Write CLI arguments
- Export FASTA and hit tables
- Generate binary dataset
- Optionally print statistics
"""

def output_operator(config: Config) -> None:
    """Run output operators to fetch/export results.

    Creates a timestamped subdirectory under `options.result_files_directory`
    and writes:
      - Command-line arguments
      - FASTA and hit outputs
      - Binary dataset

    Args:
        config: Configuration object with paths and DB reference.

    Side Effects:
        Creates new output directory with multiple files.
    """
    print_header("Output operator (fetch/export)", logger=logger)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    directory = os.path.join(config.result_files_directory, f"fetch_{ts}")
    os.makedirs(directory, exist_ok=True)

    output.print_command_line_args(os.path.join(directory, "logged_fetch_command.txt"))

    protein_dict, cluster_dict, taxon_dict = fetch_fasta_and_hit_data(config)

    requests = config.fetch_csbs + config.fetch_proteins # These are all domains from -fd and -fc
    print_reports.print_hit_reports(
        directory, protein_dict, cluster_dict, taxon_dict, requests
    )

    if config.print_fasta:
        print_reports.print_fasta_files(directory, protein_dict, cluster_dict)

    #datasets.main_binary_dataset(
    #    config, directory, protein_dict, cluster_dict, taxon_dict
    #)
    logger.info("Generated binary dataset → %s", directory)


def output_statistics(config: Config) -> None:
    """Fetch and print genome statistics from the database.

    Args:
        config: Configuration with `.database_directory`.
    """


    database.fetch_genome_statistic(config.database_directory)

def fetch_fasta_and_hit_data(
    config: Config,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Main function is to acquire the data for the print routines

    Args:
        config (object): Configuration object with output and filter parameters.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
            - protein_dict: Maps proteinID to proteinObj.
            - cluster_dict: Maps clusterID to clusterObj.
            - limiter_dict: Maps genomeID to taxonomy string.

    Example:
        >>> prot, clust, tax = fetch_fasta_and_hit_data(config)
    """
    fetch_from_gene_cluster = False

    # Limit the genomeIDs and fetch the lineage
    limiter_dict = dict()
    if config.dataset_limit_lineage:
        limiter_dict = db_fetch_taxonomy.fetch_limiter_data_keys_only(config)

    # Fetch keywords defined by the routines from the fc command
    required_proteins = set()
    if config.fetch_csbs:
        logger.info(f"Collecting gene clusters containing {config.fetch_csbs}")
        fetch_from_gene_cluster = True
        required_proteins = config.fetch_csbs

    # Fetch the protein domains from anywhere in the genome from fd command
    elif config.fetch_proteins:
        logger.info(f"Collecting proteins containing {config.fetch_proteins}")
        fetch_from_gene_cluster = False
        required_proteins = config.fetch_proteins

    # Collect the data defined by the keywords and the limiter dictionary
    logger.info("Collecting hits from local database")
    protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
        config.database_directory,
        required_proteins,
        limiter_dict,
        fetch_from_gene_cluster
    )
    return protein_dict, cluster_dict, taxon_dict
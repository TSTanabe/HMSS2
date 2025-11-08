from __future__ import annotations

import os
from datetime import datetime
from typing import Dict

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.io import print_reports, print_command_args
from hmsss.io import db_fetch_data_general
from hmsss.io import db_fetch_context
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

def _load_domain_annotations(tsv_path: str) -> Dict[str, Dict[str, str]]:
    """
    Liest Domain-Annotationen aus einer TSV:
      Spalten (Header, Tab-getrennt):
        domain    reaction    protein_description    system    metabolism
    Gibt ein Dict: domain -> { 'reaction', 'protein_description', 'system', 'metabolism' }
    """
    ann: Dict[str, Dict[str, str]] = {}
    with open(tsv_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        # Spaltenindizes robust bestimmen
        idx = {name: header.index(name) for name in
               ["domain", "reaction", "protein_description", "system", "metabolism"]}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            d = parts[idx["domain"]].strip()
            if not d:
                continue
            ann[d] = {
                "reaction": parts[idx["reaction"]].strip() if idx["reaction"] < len(parts) else "",
                "protein_description": parts[idx["protein_description"]].strip() if idx["protein_description"] < len(parts) else "",
                "system": parts[idx["system"]].strip() if idx["system"] < len(parts) else "",
                "metabolism": parts[idx["metabolism"]].strip() if idx["metabolism"] < len(parts) else "",
            }
    return ann

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

    print_command_args.print_command_line_args(os.path.join(directory, "logged_fetch_command.txt"))

    # The fetch function collects the protein_dict and taxon_dict, cluster_dict is empty
    protein_dict, cluster_dict, taxon_dict = db_fetch_data_general.fetch_fasta_and_hit_data(config)

    #metabolic_dict = _load_domain_annotations(config.metabolic_information)
    metabolic_dict= {}

    cluster_context_dict = db_fetch_context.fetch_cluster_context_for_proteins(
        database=config.database_directory,
        base_protein_dict=protein_dict,
        fetch_from_gene_clusters=bool(config.fetch_csbs),
    )

    # Combine CSB and protein requests, flattening any ':'-separated tokens
    requests = [
        token.strip()
        for group in (config.fetch_csbs or []) + (config.fetch_proteins or [])
        for token in group.split(":")
        if token.strip()
    ]

    print_reports.print_hit_reports(
        directory, protein_dict, cluster_dict, taxon_dict, metabolic_dict, cluster_context_dict, requests
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
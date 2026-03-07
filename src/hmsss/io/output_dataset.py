from __future__ import annotations

import os
from datetime import datetime
from typing import Dict

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.io import print_reports, print_command_args, print_graphs, print_reads
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
    Metabolic informatic / Metabolism information
    Liest Domain-Annotationen aus einer TSV:
      Spalten (Header, Tab-getrennt):
        domain    reaction    protein_description    system    metabolism
    Gibt ein Dict: domain -> { 'reaction', 'protein_description', 'system', 'metabolism' }
    """
    ann: Dict[str, Dict[str, str]] = {}
    with open(tsv_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        # Spaltenindizes robust bestimmen
        idx = {
            name: header.index(name)
            for name in [
                "domain",
                "reaction",
                "protein_description",
                "system",
                "metabolism",
            ]
        }
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            d = parts[idx["domain"]].strip()
            if not d:
                continue
            ann[d] = {
                "reaction": parts[idx["reaction"]].strip()
                if idx["reaction"] < len(parts)
                else "",
                "protein_description": parts[idx["protein_description"]].strip()
                if idx["protein_description"] < len(parts)
                else "",
                "system": parts[idx["system"]].strip()
                if idx["system"] < len(parts)
                else "",
                "metabolism": parts[idx["metabolism"]].strip()
                if idx["metabolism"] < len(parts)
                else "",
            }
    return ann


def _make_request_list(fetch_csbs, fetch_proteins):
    """
    Build an ordered list of request tokens from fetch_csbs and fetch_proteins.
    - Removes '[' and ']'
    - Splits groups on ':'
    - Keeps first occurrence of each token (stable order)
    """

    # 1) Tokens extrahieren
    raw_requests = [
        token.strip().replace("[", "").replace("]", "")
        for group in (fetch_csbs or []) + (fetch_proteins or [])
        for token in group.split(":")
        if token.strip()
    ]

    # 2) Duplikate entfernen (erstes Auftreten gewinnt)
    seen = set()
    requests = []
    for r in raw_requests:
        if r not in seen:
            seen.add(r)
            requests.append(r)

    return requests


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

    print_command_args.print_command_line_args(
        os.path.join(directory, "logged_fetch_command.txt")
    )
    if config.fetch_reads or config.fetch_metagenomes:
        run_read_fetch_output(config, directory)
    else:
        run_protein_fetch_output(config, directory)


def run_protein_fetch_output(config: Config, directory: str) -> None:
    # The fetch function collects the protein_dict and taxon_dict, cluster_dict is empty
    protein_dict, cluster_dict, taxon_dict = (
        db_fetch_data_general.fetch_fasta_and_hit_data(config)
    )

    # metabolic_dict = _load_domain_annotations(config.metabolic_information)
    metabolic_dict = {}

    # cluster_context_dict = db_fetch_context.fetch_cluster_context_for_proteins(
    #    database=config.database_directory,
    #    base_protein_dict=protein_dict,
    #    fetch_from_gene_clusters=bool(config.fetch_csbs),
    # )
    cluster_context_dict = {}
    # remove the brackets and double points from argument parsing and replicates
    requests = _make_request_list(config.fetch_csbs, config.fetch_proteins)

    print_reports.print_hit_reports(
        directory=directory,
        protein_dict=protein_dict,
        cluster_dict=cluster_dict,
        taxon_dict=taxon_dict,
        metabolic_dict=metabolic_dict,
        context_dict=cluster_context_dict,
        fetch_proteins=requests,
    )

    if config.print_fasta:
        print_reports.print_fasta_files(directory, protein_dict, cluster_dict)

    if config.print_graphs:
        taxonomy_levels: list[str] = config.graph_tax_levels
        print_graphs.print_hit_graphs(
            directory=directory,
            protein_dict=protein_dict,
            cluster_dict=cluster_dict,
            taxon_dict=taxon_dict,
            metabolic_dict=metabolic_dict,
            context_dict=cluster_context_dict,
            fetch_proteins=requests,
            levels=taxonomy_levels,
        )
    # options. graph taxonomy levels for the taxonomy levels
    # datasets.main_binary_dataset(
    #    config, directory, protein_dict, cluster_dict, taxon_dict
    # )
    #

    # logger.info("Generated binary dataset → %s", directory)


def output_statistics(config: Config) -> None:
    """Fetch and print genome statistics from the database.

    Args:
        config: Configuration with `.database_directory`.
    """

    database.fetch_genome_statistic(config.database_directory)


def run_read_fetch_output(config: Config, directory: str) -> None:
    """
    Fetch and print read/metagenome based output files.

    This routine is the dedicated output path for read data. It retrieves
    read placements from the database and writes read-based summary tables.

    Parameters
    ----------
    config : Config
        Global configuration object.
    directory : str
        Output directory for the fetch results.
    """
    logger.info("Running read/metagenome fetch output")

    read_dict, metagenome_dict, lineage_dict, gpkg_lengths = (
        db_fetch_data_general.fetch_read_and_hit_data(config)
    )

    logger.info(
        "Read fetch returned %s placements across %s metagenomes and %s lineages",
        len(read_dict),
        len(metagenome_dict),
        len(lineage_dict),
    )

    if not read_dict:
        logger.warning("No read placements matched the requested filters.")
        return

    print_reads.print_read_hit_reports(
        directory=directory,
        read_dict=read_dict,
        metagenome_dict=metagenome_dict,
        lineage_dict=lineage_dict,
        gpkg_length_dict=gpkg_lengths,
    )

    if config.print_fasta:
        print_reads.output_read_fastas(
            directory=directory,
            read_dict=read_dict,
        )

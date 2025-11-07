from __future__ import annotations

import os
from itertools import product
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

    output.print_command_line_args(os.path.join(directory, "logged_fetch_command.txt"))

    protein_dict, cluster_dict, taxon_dict = fetch_fasta_and_hit_data(config)

    #metabolic_dict = _load_domain_annotations(config.metabolic_information)
    metabolic_dict= {}

    requests = config.fetch_csbs + config.fetch_proteins # These are all domains from -fd and -fc
    print_reports.print_hit_reports(
        directory, protein_dict, cluster_dict, taxon_dict, metabolic_dict, requests
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

def expand_required_proteins(raw: list[str]) -> list[list[str]]:
    """
    Beispiel:
      ['a', 'b|c'] → [['a', 'b'], ['a', 'c']]
      ['a', 'b|c', 'd|e|f'] → alle 6 Kombinationen.
    """
    # Jede Position in der Eingabe wird zu einer Liste von Alternativen
    option_groups = [[p.strip() for p in token.split(":") if p.strip()] for token in raw]
    # Alle Kombinationen bilden
    combinations = [list(combo) for combo in product(*option_groups)]
    return combinations


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
    sum_protein_dict = {}
    sum_cluster_dict = {}
    sum_taxon_dict = {}

    # Limit the genomeIDs and fetch the lineage
    limiter_dict = dict()
    if config.dataset_limit_lineage:
        limiter_dict = db_fetch_taxonomy.fetch_limiter_data_keys_only(config)
    if config.fetch_genomes:
        # leere Taxonomie als Platzhalter reicht; Schlüssel sind entscheidend
        for gid in config.fetch_genomes:
            limiter_dict.setdefault(gid, {})
    # Fetch keywords defined by the routines from the fc command
    excluded_domains = config.fetch_not_csb_with_these_domains

    if config.fetch_csbs:
        logger.info(f"Collecting gene clusters containing {config.fetch_csbs}")
        fetch_from_gene_cluster = True
        raw_required_proteins = config.fetch_csbs

    # Fetch the protein domains from anywhere in the genome from fd command
    elif config.fetch_proteins:
        logger.info(f"Collecting proteins from genomes containing {config.fetch_proteins}")
        fetch_from_gene_cluster = False
        raw_required_proteins = config.fetch_proteins
    else:
        raw_required_proteins = []

    # Alle Kombinationen erzeugen
    required_combinations = expand_required_proteins(raw_required_proteins)
    logger.info(f"Collecting {len(required_combinations)} combinations of required proteins")

    for required_proteins in required_combinations:
        # Collect the data defined by the keywords and the limiter dictionary
        logger.info("Collecting hits from local database")
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=required_proteins,
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=fetch_from_gene_cluster,
            excluded_domains=excluded_domains
        )
        # Dictionaries zusammenführen (spätere Werte überschreiben ggf.)
        sum_protein_dict.update(protein_dict)
        sum_cluster_dict.update(cluster_dict)
        sum_taxon_dict.update(taxon_dict)
    return sum_protein_dict, sum_cluster_dict, sum_taxon_dict
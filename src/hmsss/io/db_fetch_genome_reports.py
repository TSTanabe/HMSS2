# src/hmsss/stages/individual_reports.py
from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from typing import List, Dict

from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_protein, db_fetch_taxonomy
from hmsss.parse_reports import parse_reports

logger = get_logger(__name__)


def _db_uri_ro_immutable(db_file: str) -> str:
    abs_db = os.path.abspath(db_file)
    return f"file:{abs_db}?mode=ro&immutable=1"


def _list_all_genome_ids(db_file: str) -> List[str]:
    db_path = _db_uri_ro_immutable(db_file)
    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()
        cur.execute("SELECT genomeID FROM Genomes;")
        return [r["genomeID"] for r in cur.fetchall()]


def write_individual_genome_reports(config) -> None:
    """
    Stage 7:
    - nimmt ALLE genomeIDs aus der DB (Genomes-Tabelle)
    - erzeugt pro genomeID eine Datei: <result>/genome_reports/<genomeID>.tsv
    - Inhalt = parse_reports.output_genome_report(...)
    """

    out_dir = Path(config.fasta_initial_hit_directory)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_ids = _list_all_genome_ids(config.database_directory)
    logger.info("Writing %d per-genome reports to %s", len(genome_ids), out_dir)

    # Taxonomie einmal holen (schnell) und pro Genom nur sub-setten
    db_uri = _db_uri_ro_immutable(config.database_directory)
    all_taxon: Dict[str, str] = db_fetch_taxonomy.fetch_taxonomy_dict(db_path=db_uri)

    excluded = getattr(config, "fetch_not_csb_with_these_domains", None)
    use_valid = getattr(config, "use_valid_hits", True)

    for i, gid in enumerate(genome_ids, start=1):
        # limiter_dict: genau dieses Genom
        limiter_dict = {gid: {}}
        logger.info(f"Printing report for {gid}")
        protein_dict, cluster_dict, _taxon_dict_unused = (
            db_fetch_protein.fetch_bulk_data(
                database=config.database_directory,
                syntenic_domains=None,  # => "alles" für dieses Genom
                limiter_dict=limiter_dict,  # => nur dieses Genom
                fetch_from_gene_clusters=False,
                excluded_domains=excluded,
                use_valid_hits=use_valid,
            )
        )

        out_file = out_dir / f"{gid}.tsv"
        parse_reports.output_genome_report(
            output_filepath=str(out_file),
            protein_dict=protein_dict,
            cluster_dict=cluster_dict,
            taxon_dict={gid: all_taxon.get(gid, "")},
            genomeID=gid,
            writemode="w",
            taxon_divider="\t",
        )

        if i % 250 == 0:
            logger.info(
                "Individual report printing progress: %d / %d", i, len(genome_ids)
            )

    logger.info("Finished: %d reports", len(genome_ids))

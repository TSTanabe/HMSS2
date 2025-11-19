"""
Fetch helpers for database exports (used by output/operator code).

Funktionen in diesem Modul:
- expand_required_proteins: unterstützt OR-Gruppen in -fd / -fc Argumenten,
  z.B. "A|B" oder "A+B" → Alternativen.
- fetch_fasta_and_hit_data: kapselt die komplette Logik zum Holen von
  protein_dict, cluster_dict, taxon_dict aus der SQLite-DB inkl.:
    * Taxonomie-Limiter
    * fetch_genomes
    * CSB vs. freie Proteinsuche
    * Alle Kombinationen von OR-Gruppen
    * Zusammenführen der Ergebnisse über alle Kombinationen
"""

import sys
import re
from itertools import product
from typing import Any, Dict, List, Tuple

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_taxonomy, db_fetch_protein

logger = get_logger(__name__)


def _split_or_token(token: str) -> List[str]:
    """
    Zerlegt einen Token anhand von ':' in Alternativen.

    Semantik:
      - 'A:B'  -> ['A', 'B']
      - 'YihQ:' -> ['YihQ', '']   # '' bedeutet: diese Position kann auch entfallen
    """
    token = token.strip()
    if not token:
        return []

    if ":" in token:
        parts = [p.strip() for p in token.split(":")]

        alts: List[str] = []
        for i, p in enumerate(parts):
            if p:
                alts.append(p)
            else:
                # Leerer letzter Teil (trailing ':') => optional
                # Beispiel: 'YihQ:' -> ['YihQ', '']
                if i == len(parts) - 1:
                    alts.append("")
                # Leere mittlere Teile (z.B. 'A::B') ignorieren wir
        return alts

    return [token]


def expand_required_proteins(raw: List[str]) -> List[List[str]]:
    """
    Expandiert eine Liste von Tokens mit OR-Gruppen zu allen Kombinationen.

    Idee:
      Jeder Eingabetoken kann eine Liste von Alternativen definieren.
      Aus dem kartesischen Produkt aller Alternativ-Listen entstehen
      vollständige Anforderungen, die jeweils an fetch_bulk_data übergeben werden.

    Beispiele:
        ["A", "B"] ->
            [["A", "B"]]

        ["A", "B|C"] ->
            [["A", "B"],
             ["A", "C"]]

        ["A", "B|C", "D|E|F"] ->
            [["A", "B", "D"],
             ["A", "B", "E"],
             ["A", "B", "F"],
             ["A", "C", "D"],
             ["A", "C", "E"],
             ["A", "C", "F"]]
    """
    if not raw:
        return []

    argument = " ".join(raw)
    groups = [g.strip() for g in re.split(r"[\[\]]", argument) if g.strip()]

    all_combos: List[List[str]] = []

    for group in groups:
        tokens = group.split(" ")

        option_groups: List[List[str]] = []

        for token in tokens:
            alts = _split_or_token(token)
            # Leere Tokens ignorieren
            if not alts:
                continue
            option_groups.append(alts)

        if not option_groups:
            return []
        # Kartesisches Produkt innerhalb der Gruppe
        for combo in product(*option_groups):
            # '' bedeutet "optional weglassen"
            filtered = [x for x in combo if x != ""]
            if filtered:  # skip completely empty list
                all_combos.append(filtered)

    return all_combos


def _build_limiter_dict(config: Config) -> Dict[str, Any]:
    """
    Erzeugt das limiter_dict:
    - optional eingeschränkt über Taxonomie (dataset_limit_lineage/-taxon)
    - erweitert um explizit angegebene fetch_genomes
    """
    limiter_dict: Dict[str, Any] = {}

    # Taxonomie-Limiter (liefert nur Keys; Inhalte sind hier egal)
    if config.dataset_limit_lineage:
        limiter_dict = db_fetch_taxonomy.fetch_limiter_data_keys_only(config)

    # Explizite Genomliste ergänzt / überschreibt die Keys
    if config.fetch_genomes:
        for gid in config.fetch_genomes:
            limiter_dict.setdefault(gid, {})

    return limiter_dict


def fetch_fasta_and_hit_data(
    config: Config,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    Zentrale Fetch-Routine für Output-Operatoren.

    Unterstützt:
      - -fc (fetch_csbs): Gene-Cluster-Modus
      - -fd (fetch_proteins): freie Proteinsuche im Genom
      - OR-Gruppen mit | oder + (z.B. 'A|B', 'C+D')
      - Kombination aller Alternativen
      - Zusammenführen der Ergebnisse aus allen Kombinationen

    Rückgabe:
        sum_protein_dict: proteinID -> proteinObj
        sum_cluster_dict: clusterID -> clusterObj
        sum_taxon_dict:   genomeID  -> taxonomy / metadata
    """
    # Limiter vorbereiten
    limiter_dict = _build_limiter_dict(config)
    excluded_domains = config.fetch_not_csb_with_these_domains

    raw_required: List[str] = []
    additional_proteins: List[str] = []

    # Quelle bestimmen: CSB oder Proteindomänen
    if config.fetch_csbs and config.fetch_proteins:
        fetch_from_gene_cluster = True
        raw_required = config.fetch_csbs
        required_combinations = expand_required_proteins(raw_required)
        logger.info(f"Collecting gene clusters containing: {required_combinations}")

        raw_required = config.fetch_proteins
        additional_proteins = expand_required_proteins(raw_required)
        logger.info(
            f"Adding proteins to genomes with these gene clusters: {additional_proteins}"
        )

    elif config.fetch_csbs:
        fetch_from_gene_cluster = True
        raw_required = config.fetch_csbs
        required_combinations = expand_required_proteins(raw_required)
        logger.info(f"Collecting gene clusters containing: {raw_required}")
    elif config.fetch_proteins:
        fetch_from_gene_cluster = False
        raw_required = config.fetch_proteins
        required_combinations = expand_required_proteins(raw_required)
        logger.info(f"Collecting proteins containing: {raw_required}")
    else:
        # Keine Angabe von domains, daher alles für die gewünschten Genome
        logger.info(f"Fetching all hits for genomes {limiter_dict.keys()}")
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=raw_required,
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=False,
            excluded_domains=excluded_domains,
            use_valid_hits=config.use_valid_hits,
        )
        return protein_dict, cluster_dict, taxon_dict

    # Sammel-Container über alle Kombinationen
    sum_protein_dict: Dict[str, Any] = {}
    sum_cluster_dict: Dict[str, Any] = {}
    sum_taxon_dict: Dict[str, Any] = {}

    # Jede Kombination sequenziell abfragen und zusammenführen
    for combo in required_combinations:
        logger.debug(f"Fetching combination: {combo}")
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=combo,
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=fetch_from_gene_cluster,
            excluded_domains=excluded_domains,
            use_valid_hits=config.use_valid_hits,
        )

        # Merge-Strategie:
        # - spätere Treffer überschreiben frühere bei gleichen Keys
        #   (vermeidet Duplikate, einfaches Verhalten)
        # - falls nötig, könnte man das später zu Aggregation anpassen
        if protein_dict:
            sum_protein_dict.update(protein_dict)
        if cluster_dict:
            sum_cluster_dict.update(cluster_dict)
        if taxon_dict:
            sum_taxon_dict.update(taxon_dict)

    # Addition von einzelnen proteinen
    for combo in additional_proteins:
        logger.debug(f"Fetching combination: {combo}")
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=combo,
            limiter_dict=sum_taxon_dict,
            fetch_from_gene_clusters=False,
            excluded_domains=excluded_domains,
            use_valid_hits=config.use_valid_hits,
        )

        # Merge-Strategie:
        # - spätere Treffer überschreiben frühere bei gleichen Keys
        #   (vermeidet Duplikate, einfaches Verhalten)
        # - falls nötig, könnte man das später zu Aggregation anpassen
        if protein_dict:
            sum_protein_dict.update(protein_dict)
        if cluster_dict:
            sum_cluster_dict.update(cluster_dict)
        if taxon_dict:
            sum_taxon_dict.update(taxon_dict)
    logger.info(
        "Fetch summary: %d proteins, %d taxa",
        len(sum_protein_dict),
        len(sum_taxon_dict),
    )

    return sum_protein_dict, sum_cluster_dict, sum_taxon_dict

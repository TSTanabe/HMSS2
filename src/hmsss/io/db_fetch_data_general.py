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

import re
import sqlite3
import sys
from itertools import product
from typing import Any, Dict, List, Set, Tuple

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read
from hmsss.io import db_fetch_protein, db_fetch_read, db_fetch_taxonomy
from hmsss.parse_reports import parse_reports

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


def _expand_optional_group_tokens(group_text: str) -> list[str]:
    """
    Für 'A B C' erzeugt:
      ['A: B C', 'A B: C', 'A B C:']
    (genau ein Token optional pro Variante)
    """
    toks = [t for t in group_text.split() if t]
    if len(toks) <= 1:
        # bei 0/1 Token macht die "genau eins optional" Semantik kaum Sinn;
        # wir behandeln es als "Token optional"
        return [f"{toks[0]}:"] if toks else []
    out = []
    for i in range(len(toks)):
        vt = toks.copy()
        vt[i] = vt[i] + ":"
        out.append(" ".join(vt))
    return out


def expand_required_proteins(raw: list[str]) -> list[list[str]]:
    if not raw:
        return []

    argument = " ".join(raw).strip()
    if not argument:
        return []

    # 1) Argument in "Chunks" zerlegen, wobei [...](:?) als Einheit behandelt wird
    #    Wir bauen eine Liste von Gruppenstrings, die später wie bisher verarbeitet werden.
    group_strings: list[str] = []

    pos = 0
    for m in re.finditer(r"\[(.*?)\](:?)", argument):
        # Text vor der Klammer: als "normale" Gruppe(n) behandeln
        prefix = argument[pos : m.start()].strip()
        if prefix:
            # prefix kann selbst mehrere bracket-freie "Gruppen" enthalten.
            # simplest: als eine Gruppe weiterreichen
            group_strings.append(prefix)

        inner = (m.group(1) or "").strip()
        has_colon = m.group(2) == ":"

        if inner:
            if has_colon:
                # 2) Makro-Expansion: genau ein Token optional pro Variante
                for variant in _expand_optional_group_tokens(inner):
                    group_strings.append(variant)
            else:
                group_strings.append(inner)

        pos = m.end()

    # Rest nach letztem Match
    suffix = argument[pos:].strip()
    if suffix:
        group_strings.append(suffix)

    # 3) Bestehende Kombinatorik: jede group_string wie bisher expandieren
    all_combos: list[list[str]] = []
    for group in group_strings:
        tokens = group.split()
        option_groups: list[list[str]] = []
        for token in tokens:
            alts = _split_or_token(token)  # bleibt unverändert
            if not alts:
                continue
            option_groups.append(alts)

        if not option_groups:
            continue

        for combo in product(*option_groups):
            filtered = [x for x in combo if x != ""]
            if filtered:
                all_combos.append(filtered)

    return all_combos


def _limiter_value_is_set(value) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip() not in ("", "0")
    return bool(value)


def _limiter_requested(config: Config) -> bool:
    return bool(
        config.fetch_genomes
        or config.dataset_limit_lineage
        or config.dataset_limit_taxon
        or _limiter_value_is_set(config.dataset_limit_proteins)
        or _limiter_value_is_set(config.dataset_limit_keywords)
    )


def _build_limiter_dict(config: Config) -> Dict[str, Any]:
    limiter_dict: Dict[str, Any] = {}

    dataset_limiter_requested = bool(
        config.dataset_limit_lineage
        or config.dataset_limit_taxon
        or _limiter_value_is_set(config.dataset_limit_proteins)
        or _limiter_value_is_set(config.dataset_limit_keywords)
    )

    if dataset_limiter_requested:
        limiter_dict = db_fetch_taxonomy.fetch_limiter_data_keys_only(config)

    # Explicit genome IDs are added to the limiter.
    if config.fetch_genomes:
        for gid in config.fetch_genomes:
            limiter_dict.setdefault(gid, {})

    return limiter_dict


def fetch_fasta_and_hit_data(
    config: Config,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    limiter_dict = _build_limiter_dict(config)

    if _limiter_requested(config) and not limiter_dict:
        logger.info("Genome limiter matched no genomes.")
        return {}, {}, {}, {}

    excluded_domains = config.fetch_not_csb_with_these_domains

    raw_required: List[str] = []
    additional_proteins: List[List[str]] = []
    required_combinations: List[List[str]] = []

    if config.fetch_csbs and config.fetch_proteins:
        fetch_from_gene_cluster = True
        required_combinations = expand_required_proteins(config.fetch_csbs)
        additional_proteins = expand_required_proteins(config.fetch_proteins)

        logger.info("Collecting gene clusters containing: %s", required_combinations)
        logger.info(
            "Adding proteins to genomes with these gene clusters: %s",
            additional_proteins,
        )

    elif config.fetch_csbs:
        fetch_from_gene_cluster = True
        required_combinations = expand_required_proteins(config.fetch_csbs)
        logger.info("Collecting gene clusters containing: %s", required_combinations)

    elif config.fetch_proteins:
        fetch_from_gene_cluster = False
        required_combinations = expand_required_proteins(config.fetch_proteins)
        logger.info("Collecting proteins containing: %s", required_combinations)

    else:
        logger.info("Fetching all hits for selected genomes.")

        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=[],
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=False,
            excluded_domains=excluded_domains,
            use_non_valid_hits=config.use_non_valid_hits,
        )

        return protein_dict, cluster_dict, taxon_dict, {}

    sum_protein_dict, sum_cluster_dict, sum_taxon_dict, sum_combo_to_genomes_dict = (
        db_fetch_protein.fetch_bulk_data_for_combinations(
            database=config.database_directory,
            combinations=required_combinations,
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=fetch_from_gene_cluster,
            excluded_domains=excluded_domains,
            use_non_valid_hits=config.use_non_valid_hits,
        )
    )

    if additional_proteins:
        if not config.fd_can_add_genomes and not sum_taxon_dict:
            logger.info(
                "Skipping additional -fd proteins because the preceding -fc search matched no genomes."
            )
        else:
            fd_limiter = {} if config.fd_can_add_genomes else sum_taxon_dict

            add_proteins, add_clusters, add_taxa, _ = (
                db_fetch_protein.fetch_bulk_data_for_combinations(
                    database=config.database_directory,
                    combinations=additional_proteins,
                    limiter_dict=fd_limiter,
                    fetch_from_gene_clusters=False,
                    excluded_domains=excluded_domains,
                    use_non_valid_hits=config.use_non_valid_hits,
                )
            )

            sum_protein_dict.update(add_proteins)
            sum_cluster_dict.update(add_clusters)
            sum_taxon_dict.update(add_taxa)

    logger.info(
        "Fetch summary: %d proteins, %d taxa, %d required combinations",
        len(sum_protein_dict),
        len(sum_taxon_dict),
        len(sum_combo_to_genomes_dict),
    )

    return (
        sum_protein_dict,
        sum_cluster_dict,
        sum_taxon_dict,
        sum_combo_to_genomes_dict,
    )


#
# Fetch routines for read mapping of metagenomes
#


def fetch_read_and_hit_data(
    config: Config,
) -> Tuple[
    Dict[Tuple[str, str, str], Read],
    Dict[str, Dict[str, Any]],
    Dict[str, Dict[str, Any]],
    Dict[str, int],
]:
    """
    Central general fetch routine for read/metagenome output mode.

    Returns
    -------
    read_dict
        Keyed by (readID, gpkg_name, metagenomeID)
    metagenome_dict
        Keyed by metagenomeID
    lineage_dict
        Keyed by lineageID
    """

    read_dict, metagenome_dict, lineage_dict = db_fetch_read.fetch_bulk_read_data(
        database=config.database_directory,
        domain_types=config.fetch_reads,  # oder eigener fetch_read_types Operator
        metagenome_ids=config.fetch_metagenomes,
    )

    gpkg_length_dict = db_fetch_read.fetch_gpkg_lengths(
        database=config.database_directory,
        domain_types=config.fetch_reads if config.fetch_reads else None,
    )

    return read_dict, metagenome_dict, lineage_dict, gpkg_length_dict

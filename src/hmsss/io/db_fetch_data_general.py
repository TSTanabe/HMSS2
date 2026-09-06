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
from typing import Any, Dict, List, Tuple, Set

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_taxonomy, db_fetch_protein, db_fetch_read

from hmsss.graft.read_models import Read

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
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
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
            use_non_valid_hits=config.use_non_valid_hits,
        )
        return protein_dict, cluster_dict, taxon_dict, {}

    # Sammel-Container über alle Kombinationen
    sum_protein_dict: Dict[str, Any] = {}
    sum_cluster_dict: Dict[str, Any] = {}
    sum_taxon_dict: Dict[str, Any] = {}
    sum_combo_to_genomes_dict: Dict[str, Any] = {}

    # Jede Kombination sequenziell abfragen und zusammenführen
    for combo in required_combinations:
        logger.debug(f"Fetching combination: {combo}")
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=combo,
            limiter_dict=limiter_dict,
            fetch_from_gene_clusters=fetch_from_gene_cluster,
            excluded_domains=excluded_domains,
            use_non_valid_hits=config.use_non_valid_hits,
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
            combo_key = tuple(combo)
            sum_combo_to_genomes_dict.setdefault(combo_key, set()).update(
                taxon_dict.keys()
            )

    # Addition von einzelnen proteinen
    for combo in additional_proteins:
        logger.debug(f"Fetching combination: {combo}")
        fd_limiter = sum_taxon_dict if not config.fd_can_add_genomes else {}
        protein_dict, cluster_dict, taxon_dict = db_fetch_protein.fetch_bulk_data(
            database=config.database_directory,
            syntenic_domains=combo,
            limiter_dict=fd_limiter,
            fetch_from_gene_clusters=False,
            excluded_domains=excluded_domains,
            use_non_valid_hits=config.use_non_valid_hits,
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

        # Für die strain variability muss hier noch die gesamtheit der genomeIDs gespeichert werden
        # combo => genomeIDs

    logger.info(
        "Fetch summary: %d proteins, %d taxa",
        len(sum_protein_dict),
        len(sum_taxon_dict),
    )

    return sum_protein_dict, sum_cluster_dict, sum_taxon_dict, sum_combo_to_genomes_dict


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

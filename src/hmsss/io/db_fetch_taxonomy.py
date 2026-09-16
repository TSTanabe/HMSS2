from __future__ import annotations

import os
import sqlite3
from typing import Any, Dict, Iterable, Optional

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


LINEAGE_COLUMNS = {
    "domain": "Superkingdom",
    "superkingdom": "Superkingdom",
    "clade": "Clade",
    "phylum": "Phylum",
    "class": "Class",
    "order": "Ordnung",
    "ordnung": "Ordnung",
    "family": "Family",
    "genus": "Genus",
    "species": "Species",
    "strain": "Strain",
}


def _db_uri(database: str) -> str:
    if database.startswith("file:"):
        return database
    return f"file:{os.path.abspath(database)}?mode=ro&immutable=1"


def _normalise_taxon(value: Any, na_value: str = "NA") -> str:
    if value is None:
        return na_value
    value = str(value).strip()
    return na_value if not value or value.upper() in {"NULL", "N/A", "NA"} else value


def _taxonomy_from_row(row: sqlite3.Row, na_value: str = "NA") -> Dict[str, str]:
    rec = {
        "Superkingdom": _normalise_taxon(row["Superkingdom"], na_value),
        "Phylum": _normalise_taxon(row["Phylum"], na_value),
        "Class": _normalise_taxon(row["Class"], na_value),
        "Order": _normalise_taxon(row["Order"], na_value),
        "Family": _normalise_taxon(row["Family"], na_value),
        "Genus": _normalise_taxon(row["Genus"], na_value),
        "Species": _normalise_taxon(row["Species"], na_value),
    }

    rec["DeepestLevel"] = na_value
    rec["DeepestValue"] = na_value
    for rank in (
        "Species",
        "Genus",
        "Family",
        "Order",
        "Class",
        "Phylum",
        "Superkingdom",
    ):
        if rec[rank] != na_value:
            rec["DeepestLevel"], rec["DeepestValue"] = rank, rec[rank]
            break

    return rec


def _lineage_column(lineage: Optional[str]) -> Optional[str]:
    if not lineage:
        return None

    column = LINEAGE_COLUMNS.get(str(lineage).strip().lower())
    if column is None:
        raise ValueError(
            f"Unsupported taxonomy rank '{lineage}'. "
            f"Allowed values: {', '.join(sorted(LINEAGE_COLUMNS))}"
        )
    return column


def _prepare_limiter_temp_tables(
    cur: sqlite3.Cursor,
    domains: Iterable[str] | None,
    keywords: Iterable[str] | None,
) -> tuple[bool, bool]:
    """
    Resolve external domain names to domain_pk once.

    Returns:
        (domain_filter_requested, keyword_filter_requested)
    """
    domain_names = sorted({str(d).strip() for d in (domains or []) if str(d).strip()})
    keyword_names = sorted({str(k).strip() for k in (keywords or []) if str(k).strip()})

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_limit_domain_names (
            domain TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)
    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_limit_domains (
            domain_pk INTEGER PRIMARY KEY
        )
    """)
    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_limit_keywords (
            keyword TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)

    cur.execute("DELETE FROM tmp_limit_domain_names")
    cur.execute("DELETE FROM tmp_limit_domains")
    cur.execute("DELETE FROM tmp_limit_keywords")

    if domain_names:
        cur.executemany(
            "INSERT INTO tmp_limit_domain_names(domain) VALUES (?)",
            ((d,) for d in domain_names),
        )
        cur.execute("""
            INSERT INTO tmp_limit_domains(domain_pk)
            SELECT dt.domain_pk
            FROM DomainTypes dt
            JOIN tmp_limit_domain_names r ON r.domain = dt.domain
        """)

    if keyword_names:
        cur.executemany(
            "INSERT INTO tmp_limit_keywords(keyword) VALUES (?)",
            ((k,) for k in keyword_names),
        )

    return bool(domain_names), bool(keyword_names)


def _build_limiter_query(
    lineage: Optional[str],
    taxon: Optional[str],
    domain_filter: bool,
    keyword_filter: bool,
    taxonomy: bool,
) -> tuple[str, list[Any]]:
    """Build the common limiter query used by both public limiter functions."""

    where: list[str] = []
    params: list[Any] = []

    lineage_column = _lineage_column(lineage)
    if lineage_column and taxon:
        where.append(f'g."{lineage_column}" LIKE ?')
        params.append(f"%{taxon}%")

    # Compatibility with previous behavior:
    # at least ONE requested domain must occur in the genome.
    if domain_filter:
        where.append("""
            EXISTS (
                SELECT 1
                FROM GenomeDomains gd
                JOIN tmp_limit_domains r ON r.domain_pk = gd.domain_pk
                WHERE gd.genome_pk = g.genome_pk
            )
        """)

    # At least ONE requested keyword must occur in a cluster of the genome.
    if keyword_filter:
        where.append("""
            EXISTS (
                SELECT 1
                FROM Clusters c
                JOIN Keywords k ON k.cluster_pk = c.cluster_pk
                JOIN tmp_limit_keywords r ON r.keyword = k.keyword
                WHERE c.genome_pk = g.genome_pk
            )
        """)

    if taxonomy:
        select = """
            SELECT g.genomeID,
                   g.Superkingdom, g.Phylum, g.Class, g.Ordnung AS "Order",
                   g.Family, g.Genus, g.Species
            FROM Genomes g
        """
    else:
        select = "SELECT g.genomeID FROM Genomes g"

    if where:
        select += "\nWHERE " + " AND ".join(where)

    return select, params


def fetch_limiter_data(config: Config) -> Dict[str, Dict[str, str]]:
    """
    Return taxonomy records for genomes passing the configured dataset limiter.

    Domain and keyword filters preserve the previous ANY-of semantics:
    at least one requested domain / keyword must be present.
    """
    lineage = config.dataset_limit_lineage
    taxon = config.dataset_limit_taxon
    domains = list(config.dataset_limit_proteins or [])
    keywords = list(config.dataset_limit_keywords or [])

    taxon_dict: Dict[str, Dict[str, str]] = {}

    with sqlite3.connect(_db_uri(config.database_directory), uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -131072")  # ~128 MiB

        domain_filter, keyword_filter = _prepare_limiter_temp_tables(
            cur, domains, keywords
        )
        sql, params = _build_limiter_query(
            lineage=lineage,
            taxon=taxon,
            domain_filter=domain_filter,
            keyword_filter=keyword_filter,
            taxonomy=True,
        )

        cur.execute(sql, params)
        for row in cur:
            taxon_dict[row["genomeID"]] = _taxonomy_from_row(row)

    logger.info("Dataset limiter retained %d genomes.", len(taxon_dict))
    return taxon_dict


def fetch_taxonomy_dict(
    db_path: str,
    genome_ids: Optional[Iterable[str]] = None,
    *,
    na_value: str = "NA",
) -> Dict[str, Dict[str, str]]:
    """
    Fetch taxonomy.

    genome_ids=None:
        return taxonomy for all genomes (backward compatible)

    genome_ids supplied:
        return taxonomy only for those genomes, using TEMP tables to avoid
        large IN (...) expressions.
    """
    taxon_dict: Dict[str, Dict[str, str]] = {}

    with sqlite3.connect(_db_uri(db_path), uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -131072")

        join = ""

        if genome_ids is not None:
            wanted = sorted({str(g).strip() for g in genome_ids if str(g).strip()})
            if not wanted:
                return {}

            cur.execute("""
                CREATE TEMP TABLE IF NOT EXISTS tmp_tax_genome_ids (
                    genomeID TEXT PRIMARY KEY
                ) WITHOUT ROWID
            """)
            cur.execute("""
                CREATE TEMP TABLE IF NOT EXISTS tmp_tax_genomes (
                    genome_pk INTEGER PRIMARY KEY
                )
            """)

            cur.execute("DELETE FROM tmp_tax_genome_ids")
            cur.execute("DELETE FROM tmp_tax_genomes")
            cur.executemany(
                "INSERT INTO tmp_tax_genome_ids(genomeID) VALUES (?)",
                ((g,) for g in wanted),
            )

            # Text IDs are resolved only once.
            cur.execute("""
                INSERT INTO tmp_tax_genomes(genome_pk)
                SELECT g.genome_pk
                FROM Genomes g
                JOIN tmp_tax_genome_ids t ON t.genomeID = g.genomeID
            """)

            join = "JOIN tmp_tax_genomes t ON t.genome_pk = g.genome_pk"

        cur.execute(f"""
            SELECT g.genomeID,
                   g.Superkingdom, g.Phylum, g.Class, g.Ordnung AS "Order",
                   g.Family, g.Genus, g.Species
            FROM Genomes g
            {join}
        """)

        for row in cur:
            taxon_dict[row["genomeID"]] = _taxonomy_from_row(row, na_value)

    logger.debug("Fetched taxonomy for %d genomes.", len(taxon_dict))
    return taxon_dict


def fetch_limiter_data_keys_only(config: Config) -> Dict[str, Dict[str, str]]:
    """
    Return only genomeIDs passing the configured dataset limiter.

    Values intentionally remain empty dicts to avoid loading taxonomy when
    callers only need a set-like genomeID limiter.
    """
    lineage = config.dataset_limit_lineage
    taxon = config.dataset_limit_taxon
    domains = list(config.dataset_limit_proteins or [])
    keywords = list(config.dataset_limit_keywords or [])

    result: Dict[str, Dict[str, str]] = {}

    with sqlite3.connect(_db_uri(config.database_directory), uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -131072")

        domain_filter, keyword_filter = _prepare_limiter_temp_tables(
            cur, domains, keywords
        )
        sql, params = _build_limiter_query(
            lineage=lineage,
            taxon=taxon,
            domain_filter=domain_filter,
            keyword_filter=keyword_filter,
            taxonomy=False,
        )

        cur.execute(sql, params)
        for row in cur:
            result[row["genomeID"]] = {}

    logger.info("Dataset limiter retained %d genome IDs.", len(result))
    return result

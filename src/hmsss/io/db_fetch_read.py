#!/usr/bin/python
from __future__ import annotations

import os
import sqlite3
from typing import Any, Dict, Iterable, Optional, Tuple

from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


def _clean_values(values: Optional[Iterable[str]]) -> list[str]:
    """
    Normalize iterable input:
    - None -> []
    - strips whitespace
    - removes empty strings
    - de-duplicates while preserving stable sorted output
    """
    if not values:
        return []
    return sorted({str(v).strip() for v in values if str(v).strip()})


def _prepare_required_domains_temp(
    cur: sqlite3.Cursor,
    required_domains: Optional[Iterable[str]],
    table_name: str = "tmp_req_domains",
) -> int:
    """
    Create and fill a TEMP table with requested domain/GPKG names.

    Returns
    -------
    int
        Number of distinct inserted values, or 0 if no filter was requested.
    """
    doms = _clean_values(required_domains)

    cur.execute(
        f"""CREATE TEMP TABLE IF NOT EXISTS {table_name} (
            domain_type TEXT PRIMARY KEY
        ) WITHOUT ROWID"""
    )
    cur.execute(f"DELETE FROM {table_name};")

    if not doms:
        return 0

    cur.executemany(
        f"INSERT OR IGNORE INTO {table_name}(domain_type) VALUES (?)",
        ((d,) for d in doms),
    )
    return cur.rowcount or len(doms)


def _prepare_required_metagenomes_temp(
    cur: sqlite3.Cursor,
    metagenome_ids: Optional[Iterable[str]],
    table_name: str = "tmp_req_metagenomes",
) -> int:
    """
    Create and fill a TEMP table with requested metagenome IDs.

    Returns
    -------
    int
        Number of distinct inserted values, or 0 if no filter was requested.
    """
    mids = _clean_values(metagenome_ids)

    cur.execute(
        f"""CREATE TEMP TABLE IF NOT EXISTS {table_name} (
            metagenomeID TEXT PRIMARY KEY
        ) WITHOUT ROWID"""
    )
    cur.execute(f"DELETE FROM {table_name};")

    if not mids:
        return 0

    cur.executemany(
        f"INSERT OR IGNORE INTO {table_name}(metagenomeID) VALUES (?)",
        ((m,) for m in mids),
    )
    return cur.rowcount or len(mids)


def _fetch_metagenome_metadata(
    cur: sqlite3.Cursor,
    metagenome_ids: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Fetch metagenome metadata.

    If metagenome_ids is None/empty, returns all metagenomes currently present
    in the DB. In normal fetch flow this helper is called with a restricted set.
    """
    mids = _clean_values(metagenome_ids)
    out: Dict[str, Dict[str, Any]] = {}

    if mids:
        cur.execute("""
            CREATE TEMP TABLE IF NOT EXISTS tmp_meta_lookup (
                metagenomeID TEXT PRIMARY KEY
            ) WITHOUT ROWID
        """)
        cur.execute("DELETE FROM tmp_meta_lookup;")
        cur.executemany(
            "INSERT OR IGNORE INTO tmp_meta_lookup(metagenomeID) VALUES (?)",
            ((m,) for m in mids),
        )

        sql = """
        SELECT
            m.metagenomeID,
            m.genomeID,
            m.forward_reads,
            m.reverse_reads,
            m.prokaryotic_fraction
        FROM Metagenomes m
        JOIN tmp_meta_lookup t
          ON t.metagenomeID = m.metagenomeID
        """
        cur.execute(sql)
    else:
        cur.execute(
            """
            SELECT
                m.metagenomeID,
                m.genomeID,
                m.forward_reads,
                m.reverse_reads,
                m.prokaryotic_fraction
            FROM Metagenomes m
            """
        )

    for row in cur:
        out[row["metagenomeID"]] = {
            "metagenomeID": row["metagenomeID"],
            "genomeID": row["genomeID"],
            "forward_reads": row["forward_reads"],
            "reverse_reads": row["reverse_reads"],
            "prokaryotic_fraction": row["prokaryotic_fraction"],
        }

    return out


def _fetch_lineage_metadata(
    cur: sqlite3.Cursor,
    lineage_ids: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Fetch lineage metadata for the supplied lineageIDs.
    """
    lids = _clean_values(lineage_ids)
    out: Dict[str, Dict[str, Any]] = {}

    if not lids:
        return out

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_lineage_lookup (
            lineageID TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_lineage_lookup;")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_lineage_lookup(lineageID) VALUES (?)",
        ((lid,) for lid in lids),
    )

    cur.execute(
        """
        SELECT
            l.lineageID,
            l.root,
            l.kingdom,
            l.phylum,
            l.class,
            l."order" AS tax_order,
            l.family,
            l.genus,
            l.species,
            l.raw_lineage
        FROM Lineage l
        JOIN tmp_lineage_lookup t
          ON t.lineageID = l.lineageID
        """
    )

    for row in cur:
        out[row["lineageID"]] = {
            "lineageID": row["lineageID"],
            "root": row["root"],
            "kingdom": row["kingdom"],
            "phylum": row["phylum"],
            "class": row["class"],
            "order": row["tax_order"],
            "family": row["family"],
            "genus": row["genus"],
            "species": row["species"],
            "raw_lineage": row["raw_lineage"],
        }

    return out


def _hydrate_read_metadata(
    read_dict: Dict[Tuple[str, str, str], Read],
    metagenome_dict: Dict[str, Dict[str, Any]],
    lineage_dict: Dict[str, Dict[str, Any]],
) -> None:
    """Attach metagenome and lineage metadata to Read objects."""

    for read in read_dict.values():
        meta = metagenome_dict.get(read.metagenomeID)
        if meta:
            read.genomeID = meta.get("genomeID") or ""

        lineage = lineage_dict.get(read.lineageID)
        if lineage:
            read.lineage = {
                "root": lineage.get("root") or "",
                "k": lineage.get("kingdom") or "",
                "p": lineage.get("phylum") or "",
                "c": lineage.get("class") or "",
                "o": lineage.get("order") or "",
                "f": lineage.get("family") or "",
                "g": lineage.get("genus") or "",
                "s": lineage.get("species") or "",
            }
        else:
            read.lineage = {}


def generate_fetch_query(
    *, use_domain_filter: bool, use_metagenome_filter: bool
) -> str:
    """Build a narrow Placement query. Metagenome and lineage metadata are hydrated afterwards."""

    join_domains = (
        "JOIN tmp_req_domains rd ON rd.domain_type = p.domain_type"
        if use_domain_filter
        else ""
    )
    join_metas = (
        "JOIN tmp_req_metagenomes rm ON rm.metagenomeID = p.metagenomeID"
        if use_metagenome_filter
        else ""
    )

    return f"""
        SELECT p.domain_type, p.readID, p.metagenomeID, p.proteinID, p.lineageID,
               p.dom_start, p.dom_end, p.coverage, p.sequence, p.alignment
        FROM Placement p
        {join_domains}
        {join_metas}
    """


def build_reads_from_query(
    cur: sqlite3.Cursor,
    sql: str,
) -> Tuple[Dict[Tuple[str, str, str], Read], set[str], set[str]]:
    """Build Read objects and collect referenced metagenome and lineage IDs."""

    read_dict: Dict[Tuple[str, str, str], Read] = {}
    metagenome_ids: set[str] = set()
    lineage_ids: set[str] = set()

    cur.execute(sql)

    for row in cur:
        mid = row["metagenomeID"]
        lid = row["lineageID"]
        key = (row["readID"], row["domain_type"], mid)

        read = Read(
            readID=row["readID"],
            alignment=row["alignment"] or "",
            sequence=row["sequence"] or "",
            gpkg_name=row["domain_type"],
            metagenomeID=mid,
            genomeID="",
        )

        read.start = row["dom_start"] if row["dom_start"] is not None else 0
        read.end = row["dom_end"] if row["dom_end"] is not None else 1
        read.coverage = row["coverage"] if row["coverage"] is not None else 1.0
        read.lineageID = lid or ""
        read.lineage = {}

        read_dict[key] = read
        metagenome_ids.add(mid)

        if lid:
            lineage_ids.add(lid)

    return read_dict, metagenome_ids, lineage_ids


def fetch_gpkg_lengths(
    database: str,
    domain_types: Optional[Iterable[str]] = None,
) -> Dict[str, int]:
    """
    Fetch gpkg / protein lengths from the GpkgLengths table.

    Parameters
    ----------
    database : str
        Path to SQLite database.
    domain_types : iterable[str] | None
        Restrict the fetch to these domain_type / gpkg names.

    Returns
    -------
    dict
        Mapping {domain_type: protein_length}
    """
    abs_db = os.path.abspath(database)
    db_path = f"file:{abs_db}?mode=ro&immutable=1"

    result: Dict[str, int] = {}
    req_domains = _clean_values(domain_types)

    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -131072")

        if req_domains:
            cur.execute("""
                CREATE TEMP TABLE IF NOT EXISTS tmp_req_gpkg_lengths (
                    domain_type TEXT PRIMARY KEY
                ) WITHOUT ROWID
            """)
            cur.execute("DELETE FROM tmp_req_gpkg_lengths;")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_req_gpkg_lengths(domain_type) VALUES (?)",
                ((d,) for d in req_domains),
            )

            cur.execute(
                """
                SELECT
                    g.domain_type,
                    g.protein_length
                FROM GpkgLengths g
                JOIN tmp_req_gpkg_lengths t
                  ON t.domain_type = g.domain_type
                """
            )
        else:
            cur.execute(
                """
                SELECT
                    g.domain_type,
                    g.protein_length
                FROM GpkgLengths g
                """
            )

        for row in cur:
            result[row["domain_type"]] = int(row["protein_length"])

    return result


def fetch_bulk_read_data(
    database: str,
    domain_types: Optional[Iterable[str]] = None,
    metagenome_ids: Optional[Iterable[str]] = None,
) -> Tuple[
    Dict[Tuple[str, str, str], Read],
    Dict[str, Dict[str, Any]],
    Dict[str, Dict[str, Any]],
]:
    """Fetch placements first and hydrate only metadata actually referenced by those placements."""

    abs_db = os.path.abspath(database)
    db_path = f"file:{abs_db}?mode=ro&immutable=1"

    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -262144")

        n_domains = _prepare_required_domains_temp(cur, domain_types)
        n_metas = _prepare_required_metagenomes_temp(cur, metagenome_ids)

        sql = generate_fetch_query(
            use_domain_filter=n_domains > 0,
            use_metagenome_filter=n_metas > 0,
        )

        logger.info(
            "Fetching reads | domain types=%s | metagenomes=%s",
            n_domains if n_domains else "ALL",
            n_metas if n_metas else "ALL",
        )

        read_dict, used_metagenomes, used_lineages = build_reads_from_query(cur, sql)

        metagenome_dict = _fetch_metagenome_metadata(cur, used_metagenomes)
        lineage_dict = _fetch_lineage_metadata(cur, used_lineages)
        _hydrate_read_metadata(read_dict, metagenome_dict, lineage_dict)

    logger.info(
        "Fetched %d placements from %d metagenomes with %d lineages.",
        len(read_dict),
        len(metagenome_dict),
        len(lineage_dict),
    )

    return read_dict, metagenome_dict, lineage_dict

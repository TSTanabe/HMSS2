from __future__ import annotations

import os
import sqlite3
from typing import Any, Dict, List, Optional, Set, Tuple, Iterable

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger

logger = get_logger(__name__)



"""
Here the routines for the protein 
"""


def fetch_limiter_data(config: Config) -> Dict[str, Dict[str,str]]:
    db = config.database_directory
    lineage = config.dataset_limit_lineage  # z.B. 'Genus'
    taxon = config.dataset_limit_taxon
    domains = list(config.dataset_limit_proteins or [])
    keywords = list(config.dataset_limit_keywords or [])
    sep = config.dataset_divide_sign

    taxon_dict: Dict[str, Dict[str,str]] = {}

    # Baue WHERE schmal und nutze EXISTS – keine LEFT JOINs, kein DISTINCT
    where = []
    params: list[Any] = []

    if lineage and taxon:
        # wenn möglich LIKE 'Taxon%' statt '%Taxon%' (Index-freundlicher)
        where.append(f'g.{lineage} LIKE ?')
        params.append(f'%{taxon}%')

    # Proteindomains in derselben GenomeID?
    if domains:
        # EXISTS: Genomes -> Proteins -> Domains
        q_marks = ','.join(['?'] * len(domains))
        where.append(f"""EXISTS (
            SELECT 1
            FROM Proteins p
            JOIN Domains d ON d.proteinID = p.proteinID
            WHERE p.genomeID = g.genomeID
              AND d.domain IN ({q_marks})
        )""")
        params.extend(domains)

    # Keywords in Clustern derselben GenomeID?
    if keywords:
        q_marks = ','.join(['?'] * len(keywords))
        where.append(f"""EXISTS (
            SELECT 1
            FROM Clusters c
            JOIN Keywords k ON k.clusterID = c.clusterID
            WHERE c.genomeID = g.genomeID
              AND k.keyword IN ({q_marks})
        )""")
        params.extend(keywords)

    sql = (
        'SELECT g.genomeID, g.Superkingdom, g.Phylum, g.Class, g."Order", '
        '       g.Family, g.Genus, g.Species '
        'FROM Genomes g ' + (('WHERE ' + ' AND '.join(where)) if where else '')
    )

    # Read-only & immutable öffnet schneller/sicherer auf NFS
    con = sqlite3.connect(f"file:{db}?mode=ro&immutable=1", uri=True)
    con.execute("PRAGMA query_only = ON;")          # verbietet Schreiboperationen
    con.execute("PRAGMA journal_mode = OFF;")       # kein Journal nötig (read-only)
    con.execute("PRAGMA synchronous = OFF;")        # keine Syncs (nur lesend)
    try:
        cur = con.execute(sql, params)
        for row in cur:
            gid = row[0]
            # Rohwerte -> Strings, Leer/None => 'NA'
            raw = {
                "Superkingdom": row[1],
                "Phylum":       row[2],
                "Class":        row[3],
                "Order":        row[4],
                "Family":       row[5],
                "Genus":        row[6],
                "Species":      row[7],
            }
            norm: dict[str, str] = {
                k: (str(v).strip() if (v is not None and str(v).strip() != "") else "NA")
                for k, v in raw.items()
            }

            # tiefste nicht-NA Ebene bestimmen
            depth_order = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Superkingdom"]
            deepest_level = next((lvl for lvl in depth_order if norm.get(lvl, "NA") != "NA"), "NA")
            deepest_value = norm.get(deepest_level, "NA") if deepest_level != "NA" else "NA"
            norm["DeepestLevel"] = deepest_level
            norm["DeepestValue"] = deepest_value


            taxon_dict[gid] = norm

    finally:
        con.close()

    return taxon_dict

def fetch_taxonomy_dict(
    db_path: str,
    *,
    na_value: str = "NA",
) -> Dict[str, Dict[str, str]]:

    def _norm(v, na_value="NA") -> str:
        # Debug hilft mehr mit repr:
        # print(f">>{v!r}<<")
        if v is None:
            return na_value
        s = str(v).strip()
        if s == "" or s.upper() in {"NULL", "N/A", "NA"}:
            return na_value
        return s

    # Vorinitialisieren: jede angefragte ID ist enthalten (alles NA)
    taxon_dict: Dict[str, Dict[str, str]] = {}

    con = sqlite3.connect(db_path, uri=True)
    try:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # Temp schneller, aber kein query_only setzen, weil wir TEMP schreiben
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA temp_store = MEMORY;")


        cur.execute("""
            SELECT
              g.genomeID     AS genomeID,
              g.Superkingdom AS Superkingdom,
              g.Phylum       AS Phylum,
              g.Class        AS Class,
              g.Ordnung      AS "Order",
              g.Family       AS Family,
              g.Genus        AS Genus,
              g.Species      AS Species
            FROM Genomes g
        """)

        for row in cur:
            gid = row["genomeID"]
            rec = {
                "Superkingdom": _norm(row["Superkingdom"]),
                "Phylum":       _norm(row["Phylum"]),
                "Class":        _norm(row["Class"]),
                "Order":        _norm(row["Order"]),
                "Family":       _norm(row["Family"]),
                "Genus":        _norm(row["Genus"]),
                "Species":      _norm(row["Species"]),
            }
            # tiefste nicht-NA Ebene bestimmen (beim Schreiben)
            for lvl in ("Species","Genus","Family","Order","Class","Phylum","Superkingdom"):
                if rec[lvl] != na_value:
                    rec["DeepestLevel"] = lvl
                    rec["DeepestValue"] = rec[lvl]
                    break
                else:
                    rec["DeepestLevel"] = na_value
                    rec["DeepestValue"] = na_value

            taxon_dict[gid] = rec
    finally:
        con.close()

    return taxon_dict


def fetch_limiter_data_keys_only(config: Config) -> Dict[str, Dict[str, str]]:
    """
    Liefert ein taxon_dict, das NUR die genomeIDs als Keys enthält.
    Die Values sind absichtlich leere Dicts {}, um Zeit zu sparen.
    Filtering wie in fetch_limiter_data (lineage/taxon, domains, keywords),
    aber es werden KEINE Taxonomie-Spalten geladen.

    Returns:
        Dict[str, Dict[str, str]]: { genomeID: {} , ... }
    """
    db = config.database_directory
    lineage  = config.dataset_limit_lineage
    taxon    = config.dataset_limit_taxon
    domains  = list(config.dataset_limit_proteins or [])
    keywords = list(config.dataset_limit_keywords or [])

    taxon_dict: Dict[str, Dict[str, str]] = {}

    # WHERE-Bedingungen wie gehabt
    where = []
    params: list[Any] = []

    if lineage and taxon:
        where.append(f'g.{lineage} LIKE ?')
        params.append(f'%{taxon}%')

    if domains:
        q_marks = ','.join(['?'] * len(domains))
        where.append(f"""EXISTS (
            SELECT 1
            FROM Proteins p
            JOIN Domains d ON d.proteinID = p.proteinID
            WHERE p.genomeID = g.genomeID
              AND d.domain IN ({q_marks})
        )""")
        params.extend(domains)

    if keywords:
        q_marks = ','.join(['?'] * len(keywords))
        where.append(f"""EXISTS (
            SELECT 1
            FROM Clusters c
            JOIN Keywords k ON k.clusterID = c.clusterID
            WHERE c.genomeID = g.genomeID
              AND k.keyword IN ({q_marks})
        )""")
        params.extend(keywords)

    sql = (
        'SELECT g.genomeID '
        'FROM Genomes g ' + (('WHERE ' + ' AND '.join(where)) if where else '')
    )

    con = sqlite3.connect(f"file:{db}?mode=ro&immutable=1", uri=True)
    con.execute("PRAGMA query_only = ON;")
    con.execute("PRAGMA journal_mode = OFF;")
    con.execute("PRAGMA synchronous = OFF;")
    try:
        cur = con.execute(sql, params)
        for (gid,) in cur:
            taxon_dict[gid] = {}  # Werte absichtlich leer
    finally:
        con.close()

    return taxon_dict
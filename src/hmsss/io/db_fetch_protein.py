#!/usr/bin/python
import os
import sqlite3
from typing import Any, Dict, List, Optional, Set, Tuple, Iterable

from hmsss.io import db_fetch_taxonomy
from hmsss.parse_reports import parse_reports
from hmsss.utils import myUtil
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def fetch_bulk_data(
    database: str,
    syntenic_domains: Optional[List[str]],
    limiter_dict: Optional[Dict[str, str]] = None,
    fetch_from_gene_clusters: bool = False,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, dict[str, str]]]:
    """
    Fetch bulk data from the database based on specified conditions, using batching
    to avoid SQLite's variable limit.

    This version expects the SELECT produced by `generate_fetch_query(...)` to provide
    stable, unique column aliases. Specifically, the following aliases are used here:

      proteinID, genomeID, clusterID,
      contig, gene_start, gene_end, gene_strand, protein_sequence,
      domain, domStart, domEnd, score,
      dom_count, comment

    Implementation notes:
      - Uses sqlite3.Row for name-based access to row fields (avoids index errors).
      - Keeps your existing flow: build Protein objects on-the-fly, collect Cluster
        stubs (one per clusterID), then enrich clusters with Keywords and genomes
        with taxonomy in batched queries.
      - `min_cluster_completeness` is available for optional filtering after keyword
        hydration (left unchanged here to preserve current behavior).
    """
    protein_dict: Dict[str, Any] = {}
    cluster_dict: Dict[str, Any] = {}
    genome_id_set: Set[str] = set()
    fusion_prot_ids: Set[str] = set()
    if limiter_dict is None:
        limiter_dict = {}

    abs_db = os.path.abspath(database)
    db_path = f"file:{abs_db}?mode=ro&immutable=1"
    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # Pragmas
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA cache_size = 100000;")   # ~100k Pages (~100k * 1.5–2 KB je nach build)
        cur.execute("PRAGMA synchronous = OFF;")

        _prepare_required_domains_temp(cur, syntenic_domains)
        n = _prepare_limiter_genomes_temp(cur, limiter_dict)

        if fetch_from_gene_clusters:
            sql, args = generate_fetch_query_covering_domains(
                set(syntenic_domains), use_limiter=(n>1)
            )
        else:
            sql, args = generate_fetch_query_domains_anywhere( set(syntenic_domains), use_limiter=(n>1))

        #print(sql)
        #print(args)

        cur.execute(sql, args)

        for index, row in enumerate(cur):
            protein_id = row["proteinID"]
            genome_id = row["genomeID"]
            cluster_id = row["clusterID"]

            domain = row["domain"]
            dom_start = row["domStart"]
            dom_end = row["domEnd"]
            score = row["score"]

            # Create-or-extend Protein object
            if protein_id in protein_dict:
                protein_dict[protein_id].add_domain(domain, dom_start, dom_end, score)
            else:
                p = parse_reports.Protein(protein_id, domain, dom_start, dom_end, score)
                p.genomeID = genome_id
                p.clusterID = cluster_id
                p.gene_contig = row["contig"]
                p.gene_start = row["gene_start"]
                p.gene_end = row["gene_end"]
                p.gene_strand = row["gene_strand"]
                p.protein_sequence = row["protein_sequence"]
                p.selection_comment = row["comment"]
                p.alternative_hit = row["alternative_hit"]
                protein_dict[protein_id] = p
                genome_id_set.add(genome_id)

            # Track proteins that appear to be fused (>=2 domains)
            if row["dom_count"] >= 2:
                fusion_prot_ids.add(protein_id)

        logger.info(f"Fetched {len(protein_dict)} proteins.")

    # 2) Add all domains for fused proteins (batched)
    if fusion_prot_ids:
        hydrate_fused_protein_domains(db_path, fusion_prot_ids, protein_dict)

    # 4) Get taxonomy
    taxon_dict = db_fetch_taxonomy.fetch_taxonomy_dict(
        db_path=db_path,
    )

    return protein_dict, cluster_dict, taxon_dict

################
####   Generate fetch query
################

def _batched(iterable: Iterable[str], n: int):
    it = iter(iterable)
    while True:
        chunk = list([x for _, x in zip(range(n), it)])
        if not chunk:
            return
        yield chunk

def _prepare_required_domains_temp(cur: sqlite3.Cursor, required_domains: "Iterable[str]") -> int:
    """
    Erstellt/leer die TEMP-Tabelle tmp_req_domains(domain TEXT PRIMARY KEY)
    und befüllt sie in einem Rutsch (executemany).
    Returns: Anzahl eingefügter unterschiedlicher Domains.
    """
    doms = {d for d in required_domains if d}  # entdoppeln + leere raus
    if not doms:
        raise ValueError("required_domains must be a non-empty iterable of strings.")

    cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_req_domains (domain TEXT PRIMARY KEY);")
    cur.execute("DELETE FROM tmp_req_domains;")
    cur.executemany("INSERT OR IGNORE INTO tmp_req_domains(domain) VALUES (?)", ((d,) for d in doms))
    return cur.rowcount or len(doms)

def _prepare_limiter_genomes_temp(
    cur: sqlite3.Cursor, taxon_dict: Optional[Dict[str, Any]]
) -> int:
    """
    Erstellt TEMP-Tabelle tmp_req_genomes(genomeID TEXT PRIMARY KEY)
    und befüllt sie, falls taxon_dict nicht leer ist.

    Returns:
        Anzahl eingefügter GenomeIDs, oder 0 falls taxon_dict leer war.
    """
    if not taxon_dict:
        # Nichts zu tun → keine Einschränkung
        return 0

    gids = {str(g).strip() for g in taxon_dict.keys() if str(g).strip()}
    if not gids:
        return 0

    cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_req_genomes (genomeID TEXT PRIMARY KEY);")
    cur.execute("DELETE FROM tmp_req_genomes;")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_req_genomes(genomeID) VALUES (?)",
        ((g,) for g in gids)
    )
    return cur.rowcount or len(gids)

def generate_fetch_query_covering_domains(
    required_domains: Iterable[str], use_limiter: bool = True
) -> Tuple[str, List[Any]]:
    """
    Liefert ein SELECT, das ALLE Proteine (mit Domains) aus genau den Clustern zurückgibt,
    die das komplette Set der geforderten Domains enthalten.

    Wenn use_limiter=True, wird zusätzlich tmp_req_genomes benutzt
    (muss vorher befüllt sein). Ist taxon_dict leer, übergibt man use_limiter=False.

    Rückgabe: (SQL, [])
    """
    # Basis-CTEs: req + optional lim
    sql = """
    WITH req AS (
        SELECT domain FROM tmp_req_domains
    )
    """
    if use_limiter:
        sql += """,
    lim AS (
        SELECT genomeID FROM tmp_req_genomes
    )
    """

    sql += """,
    clusters_covering AS (
        SELECT p.clusterID AS clusterID
        FROM Proteins p
        {join_limiter}
        JOIN Domains d ON d.proteinID = p.proteinID
        JOIN req r     ON r.domain    = d.domain
        GROUP BY p.clusterID
        HAVING COUNT(DISTINCT r.domain) = (SELECT COUNT(*) FROM req)
    )
    SELECT
        p.proteinID        AS proteinID,
        p.genomeID         AS genomeID,
        p.clusterID        AS clusterID,
        p.contig           AS contig,
        p.start            AS gene_start,
        p.end              AS gene_end,
        p.strand           AS gene_strand,
        p.sequence         AS protein_sequence,
        d.domain           AS domain,
        d.domStart         AS domStart,
        d.domEnd           AS domEnd,
        d.score            AS score,
        p.dom_count        AS dom_count,
        p.comment          AS comment,
        p.alternative_hit  AS alternative_hit
    FROM Proteins p
        {join_limiter2}
    JOIN clusters_covering c ON c.clusterID = p.clusterID
    JOIN Domains d           ON d.proteinID = p.proteinID
    """

    # Platzhalter fürs optionale JOIN ersetzen
    join_txt = "JOIN lim lg ON lg.genomeID = p.genomeID" if use_limiter else ""
    sql = sql.format(join_limiter=join_txt, join_limiter2=join_txt)

    return sql, []

def generate_fetch_query_domains_anywhere(
    required_domains: Iterable[str],
    use_limiter: bool = True,
) -> Tuple[str, List[Any]]:
    """
    Baut ein SELECT, das *alle* Domain-Hits (egal wo lokalisiert, unabhängig von Clustern)
    für ein gegebenes Set von Domains zurückliefert.

    Limiter-Logik:
      - Wenn use_limiter=True, wird die TEMP-Tabelle tmp_req_genomes genutzt
        (muss zuvor via prepare_limiter_genomes_temp(...) befüllt sein).
      - Wenn use_limiter=False (oder Limiter leer), erfolgt *keine* GenomeID-Einschränkung.

    VORAUSSETZUNGEN (analog zur ersten Routine):
      - prepare_required_domains_temp(cur, required_domains) wurde ausgeführt
        und befüllt die TEMP-Tabelle tmp_req_domains(domain).
      - Optional: prepare_limiter_genomes_temp(cur, taxon_dict) für den Genome-Limiter.

    Rückgabe-Spalten (identisch zu deiner bestehenden generate_fetch_query):
      proteinID, genomeID, clusterID, contig, gene_start, gene_end, gene_strand,
      protein_sequence, domain, domStart, domEnd, score, dom_count, comment, alternative_hit

    Returns:
      (sql, args) — args ist leer, da wir mit TEMP-Tabellen arbeiten (keine 999-Placeholder-Probleme).
    """
    sql = """
    WITH req AS (
        SELECT domain FROM tmp_req_domains
    )
    """
    if use_limiter:
        sql += """,
    lim AS (
        SELECT genomeID FROM tmp_req_genomes
    )
    """

    sql += """
    SELECT
        p.proteinID        AS proteinID,
        p.genomeID         AS genomeID,
        p.clusterID        AS clusterID,
        p.contig           AS contig,
        p.start            AS gene_start,
        p.end              AS gene_end,
        p.strand           AS gene_strand,
        p.sequence         AS protein_sequence,
        d.domain           AS domain,
        d.domStart         AS domStart,
        d.domEnd           AS domEnd,
        d.score            AS score,
        p.dom_count        AS dom_count,
        p.comment          AS comment,
        p.alternative_hit  AS alternative_hit
    FROM Domains d
    JOIN req r      ON r.domain   = d.domain
    JOIN Proteins p ON p.proteinID = d.proteinID
    {join_limiter}
    """

    join_limiter = "JOIN lim lg ON lg.genomeID = p.genomeID" if use_limiter else ""
    sql = sql.format(join_limiter=join_limiter)

    return sql, []

#

def build_proteins_from_rows(
    rows: Iterable[sqlite3.Row],
    protein_dict: Dict[str, "parse_reports.Protein"],
    genome_id_set: Set[str],
) -> Tuple[int, int, Set[str]]:
    """
    Baut/erweitert Protein-Objekte aus SQL-Zeilen.

    Erwartete Spalten je Row:
      proteinID, genomeID, clusterID, contig, gene_start, gene_end, gene_strand,
      protein_sequence, comment, alternative_hit, domain, domStart, domEnd, score, dom_count

    Args:
        rows: Iterable von sqlite3.Row (z.B. direkter Cursor nach cur.execute(...))
        protein_dict: Ziel-Dict proteinID -> Protein-Objekt (wird in-place erweitert)
        genome_id_set: Set der gesehenen genomeIDs (wird in-place erweitert)

    Returns:
        (rows_processed, proteins_new, fusion_prot_ids)
    """
    log = logger
    rows_processed = 0
    proteins_new = 0
    fusion_prot_ids: Set[str] = set()

    for rows_processed, row in enumerate(rows, start=1):
        protein_id = row["proteinID"]
        genome_id  = row["genomeID"]
        cluster_id = row["clusterID"]

        domain     = row["domain"]
        dom_start  = row["domStart"]
        dom_end    = row["domEnd"]
        score      = row["score"]

        p = protein_dict.get(protein_id)
        if p is not None:
            p.add_domain(domain, dom_start, dom_end, score)
        else:
            p = parse_reports.Protein(protein_id, domain, dom_start, dom_end, score)
            p.genomeID          = genome_id
            p.clusterID         = cluster_id
            p.gene_contig       = row["contig"]
            p.gene_start        = row["gene_start"]
            p.gene_end          = row["gene_end"]
            p.gene_strand       = row["gene_strand"]
            p.protein_sequence  = row["protein_sequence"]
            p.selection_comment = row["comment"]
            p.alternative_hit   = row["alternative_hit"]
            protein_dict[protein_id] = p
            genome_id_set.add(genome_id)
            proteins_new += 1

        dom_count = row["dom_count"]
        if dom_count is not None and int(dom_count) >= 2:
            fusion_prot_ids.add(protein_id)

        if log and (rows_processed % 10000 == 0):
            log.info(f"Processed rows: {rows_processed} | Proteins: {len(protein_dict)}")

    if log:
        log.info(f"Fetched {len(protein_dict)} proteins (processed {rows_processed} rows).")
        log.info(f"Found {len(fusion_prot_ids)} fused proteins.")

    return rows_processed, proteins_new, fusion_prot_ids


# Fused protein fetch

def hydrate_fused_protein_domains(
    db_path: str,
    fusion_prot_ids: "set[str] | list[str]",
    protein_dict: dict,
) -> int:
    """
    Fügt für alle Proteine in `fusion_prot_ids` sämtliche Domains hinzu – effizient und 999-sicher.
    Implementierung:
      1) TEMP-Tabelle tmp_fused_ids(proteinID TEXT PRIMARY KEY) befüllen (executemany).
      2) Ein einziges SELECT JOIN Domains d ON d.proteinID = tmp_fused_ids.proteinID.
      3) Iteration über Resultate: protein_dict[pid].add_domain(...)

    Returns:
        Anzahl der hinzugefügten Domains (Zeilen aus Domains).
    """

    if not fusion_prot_ids:
        return 0

    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # TEMP-Tabelle anlegen & leeren
        cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_fused_ids (proteinID TEXT PRIMARY KEY);")
        cur.execute("DELETE FROM tmp_fused_ids;")

        # Atomare Bulk-Inserts ohne manuelles BEGIN/COMMIT:
        # Der 'with con:' Kontext oben sorgt für Transaktion pro Block.
        cur.executemany(
            "INSERT OR IGNORE INTO tmp_fused_ids(proteinID) VALUES (?)",
            ((pid,) for pid in fusion_prot_ids)
        )

        # Domains in einem Rutsch joinen
        cur.execute("""
                    SELECT d.proteinID AS proteinID,
                           d.domain AS domain,
                d.domStart  AS domStart,
                d.domEnd    AS domEnd,
                d.score     AS score
                    FROM Domains d
                        JOIN tmp_fused_ids t
                    ON t.proteinID = d.proteinID
                    """)

        added = 0
        for i, r in enumerate(cur):
            pid = r["proteinID"]
            protein = protein_dict.get(pid)
            if protein is not None:
                protein.add_domain(r["domain"], r["domStart"], r["domEnd"], r["score"])
                added += 1
                if (i + 1) % 25000 == 0:
                    logger.debug(f"Added fused domains rows: {i + 1}")
            # Wenn protein fehlt, silently skip (oder optional warnen)

        logger.info(f"Added {added} fused-domain rows for {len(fusion_prot_ids)} proteins (single-pass).")
    return added



def fetch_taxonomy_dict(
    db_path: str,
    genome_ids: Iterable[str],
    trennzeichen: str,
    existing: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Holt Taxonomie-Infos für die gegebenen genomeIDs in EINEM Query, 999-sicher.
    - Nutzt eine TEMP-Tabelle für die IDs.
    - Verwendet myUtil.taxonomy_lineage(row, trennzeichen) für die Formatierung.
    - 'existing' kann ein bereits teilweise gefülltes taxon_dict sein.

    Returns:
        Dict[str, str]: genomeID -> taxonomy_lineage
    """
    taxon_dict: Dict[str, str] = dict(existing or {})
    # IDs, die noch fehlen
    wanted = [gid for gid in set(genome_ids) if gid not in taxon_dict]
    if not wanted:
        return taxon_dict
    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_tax_fetch (id TEXT PRIMARY KEY);")
        cur.execute("DELETE FROM tmp_tax_fetch;")

        cur.executemany(
            "INSERT OR IGNORE INTO tmp_tax_fetch(id) VALUES (?)",
            ((g,) for g in wanted)
        )

        cur.execute("""
                    SELECT genomeID     AS genomeID,
                           Superkingdom AS Superkingdom,
                           Phylum       AS Phylum,
                           Class        AS Class,
                           Ordnung      AS Ordnung,
                           Family       AS Family,
                           Genus        AS Genus,
                           Species      AS Species
                    FROM Genomes
                             JOIN tmp_tax_fetch t ON t.id = Genomes.genomeID
                    """)

        added = 0
        for i, r in enumerate(cur):
            gid = r["genomeID"]
            taxon_dict[gid] = myUtil.taxonomy_lineage(r, trennzeichen)
            added += 1
            if (i + 1) % 10000 == 0:
                logger.debug(f"fetch_taxonomy_dict: {i + 1} Zeilen verarbeitet.")


        logger.info(f"fetch_taxonomy_dict: {added} Taxonomie-Zeilen für {len(wanted)} genomeIDs hinzugefügt.")
    return taxon_dict
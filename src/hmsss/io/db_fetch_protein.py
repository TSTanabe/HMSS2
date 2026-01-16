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
        excluded_domains: Optional[List[str]] = None,
        use_non_valid_hits: bool = False,
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
        cur.execute(
            "PRAGMA cache_size = 100000;"
        )  # ~100k Pages (~100k * 1.5–2 KB je nach build)
        cur.execute("PRAGMA synchronous = OFF;")
        # excluded_domains = ['sHdrB2']

        _prepare_required_domains_temp(cur, syntenic_domains)
        _prepare_excluded_domains_temp(cur, excluded_domains)
        n = _prepare_limiter_genomes_temp(cur, limiter_dict)

        # print(syntenic_domains)
        # print(excluded_domains)
        # Delete row if sqliteDB downwards compatibility is not an issue anymore
        valid_hit_column_available = has_column(cur, "Proteins", "valid_hit")
        if fetch_from_gene_clusters:
            logger.info(f"Searching for syntenic {syntenic_domains}")
            sql, args = generate_fetch_query_covering_domains(
                set(syntenic_domains),
                use_limiter=(n > 0),
                use_exclusions=True,
                use_non_valid_hits=use_non_valid_hits,  # Default True
                valid_hit_column_available=valid_hit_column_available,
            )  # Fetches all domains that are in a syntenic gene cluster, but not csb including the exclusion
        else:
            logger.info(f"Searching for co-occuring {syntenic_domains}")
            sql, args = generate_fetch_query_domains_anywhere_excluding_clusters(
                use_limiter=(n > 0),  # -fg wirklich anwenden
                use_exclusions=True,
                require_all_domains_in_same_genome=bool(
                    syntenic_domains
                ),
                # nur fordern, wenn explizite Domains übergeben wurden. Kann leer sein, wenn komplettes genom gefordert
                use_non_valid_hits=use_non_valid_hits,  # Default True
                valid_hit_column_available=valid_hit_column_available,
            )

        fusion_prot_ids: Set[str] = set()
        build_proteins_from_query(
            cur=cur,
            sql=sql,
            args=args,
            protein_dict=protein_dict,
            genome_id_set=genome_id_set,
            fusion_prot_ids=fusion_prot_ids,
        )
        logger.info(f"Fetched {len(protein_dict)} proteins.")

    # 2) Add all domains for fused proteins (batched)
    if fusion_prot_ids:
        hydrate_fused_protein_domains(db_path, fusion_prot_ids, protein_dict)

    # 4) Get taxonomy
    taxon_dict = db_fetch_taxonomy.fetch_taxonomy_dict(
        db_path=db_path,
    )
    # 5) Filter taxonomy entries to only include genomes present in protein_dict
    genome_ids_in_proteins = {
        p.genomeID for p in protein_dict.values() if getattr(p, "genomeID", None)
    }
    before = len(taxon_dict)
    taxon_dict = {
        gid: rec for gid, rec in taxon_dict.items() if gid in genome_ids_in_proteins
    }
    after = len(taxon_dict)
    logger.info(f"Hits were present in {after} genome lineages of {before}.")

    parse_reports.define_best_score_hits_for_protein_dict(protein_dict)
    return protein_dict, cluster_dict, taxon_dict


################
####   Generate fetch query
################


def _prepare_required_domains_temp(
        cur: sqlite3.Cursor, required_domains: "Iterable[str]"
) -> int:
    """
    Legt die TEMP-Tabelle tmp_req_domains(domain TEXT PRIMARY KEY) an und befüllt sie.

    Verhalten:
    - Wenn required_domains leer oder None → es werden ALLE Domains aus der DB geladen.
    - Sonst → nur die übergebenen Domains (entdoppelt, ohne leere Strings).

    Returns: Anzahl eingefügter unterschiedlicher Domains.
    """

    # TEMP-Tabelle anlegen & leeren
    cur.execute(
        "CREATE TEMP TABLE IF NOT EXISTS tmp_req_domains (domain TEXT PRIMARY KEY);"
    )
    cur.execute("DELETE FROM tmp_req_domains;")

    # Domains normalisieren
    doms = {d for d in (required_domains or []) if d}

    if not doms:
        # → keine Vorgabe: alle Domains holen
        cur.execute(
            "INSERT OR IGNORE INTO tmp_req_domains(domain) SELECT DISTINCT domain FROM Domains;"
        )
        return cur.rowcount

    # → gewählte Domains einfügen
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_req_domains(domain) VALUES (?)", ((d,) for d in doms)
    )
    return cur.rowcount or len(doms)


def _prepare_excluded_domains_temp(
        cur: sqlite3.Cursor, excluded_domains: "Iterable[str] | None"
) -> int:
    doms = {d for d in (excluded_domains or []) if d}
    cur.execute(
        "CREATE TEMP TABLE IF NOT EXISTS tmp_excl_domains (domain TEXT PRIMARY KEY);"
    )
    cur.execute("DELETE FROM tmp_excl_domains;")
    if not doms:
        return 0
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_excl_domains(domain) VALUES (?)",
        ((d,) for d in doms),
    )
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

    cur.execute(
        "CREATE TEMP TABLE IF NOT EXISTS tmp_req_genomes (genomeID TEXT PRIMARY KEY);"
    )
    cur.execute("DELETE FROM tmp_req_genomes;")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_req_genomes(genomeID) VALUES (?)",
        ((g,) for g in gids),
    )
    return cur.rowcount or len(gids)


def has_column(cur: sqlite3.Cursor, table_name: str, column_name: str) -> bool:
    """
    Prüft, ob eine bestimmte Spalte in einer Tabelle existiert.
    Nutzt den bestehenden Cursor/Connection-Kontext.
    """
    try:
        cur.execute(f"PRAGMA table_info({table_name});")
        return any(row[1] == column_name for row in cur.fetchall())
    except sqlite3.Error as e:
        logger.warning(
            f"Fehler bei Prüfung der Spalte '{column_name}' in Tabelle '{table_name}': {e}"
        )
        return False


def generate_fetch_query_covering_domains(
        required_domains: Iterable[str],
        use_limiter: bool = True,
        use_exclusions: bool = True,
        use_non_valid_hits: bool = True,
        valid_hit_column_available: bool = False,
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

    # CTE: alle Cluster, die das komplette erforderliche Domain-Set abdecken
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
    """

    # CTE: alle Cluster, die mindestens eine ausgeschlossene Domäne enthalten
    if use_exclusions:
        sql += """,
    clusters_excluded AS (
        SELECT p.clusterID AS clusterID
        FROM Proteins p
        {join_limiter2}
        JOIN Domains d   ON d.proteinID = p.proteinID
        JOIN tmp_excl_domains e ON e.domain = d.domain
        GROUP BY p.clusterID
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
    FROM Proteins p
        {join_limiter3}
    JOIN clusters_covering c ON c.clusterID = p.clusterID
    {left_join_excl}
    JOIN Domains d           ON d.proteinID = p.proteinID
    {where_clause}
    """

    join_txt = "JOIN lim lg ON lg.genomeID = p.genomeID" if use_limiter else ""
    left_join_excl = (
        "LEFT JOIN clusters_excluded x ON x.clusterID = p.clusterID"
        if use_exclusions
        else ""
    )

    # flexible WHERE-Klausel
    where_parts = []
    if use_exclusions:
        where_parts.append("x.clusterID IS NULL")
    if not use_non_valid_hits and valid_hit_column_available:
        where_parts.append("p.valid_hit = 1")
    where_clause = f"WHERE {' AND '.join(where_parts)}" if where_parts else ""

    sql = sql.format(
        join_limiter=join_txt,
        join_limiter2=join_txt,
        join_limiter3=join_txt,
        left_join_excl=left_join_excl,
        where_clause=where_clause,
    )
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


# Testing routine
def generate_fetch_query_domains_anywhere_excluding_clusters(
        use_limiter: bool = False,
        use_exclusions: bool = True,
        require_all_domains_in_same_genome: bool = True,
        use_non_valid_hits: bool = True,
        valid_hit_column_available: bool = False,
) -> tuple[str, list]:
    """
    Selektiert alle Domain-Hits aus tmp_req_domains, schließt aber Proteine aus
    Clustern aus, in denen irgendeine Domäne aus tmp_excl_domains vorkommt.
    Optional: nur Genomes zulassen, die *alle* gewünschten Domänen enthalten.

    Erwartete TEMP-Tabellen:
      - tmp_req_domains(domain TEXT)          (Pflicht)
      - tmp_excl_domains(domain TEXT)         (wenn use_exclusions=True)
      - tmp_req_genomes(genomeID TEXT)        (wenn use_limiter=True)
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

    # 1) Genomes finden, die *alle* gewünschten Domänen haben (mind. je 1 Hit)
    #    -> zählt DISTINCT req.domains pro genomeID
    if require_all_domains_in_same_genome:
        sql += """,
    req_count AS (
        SELECT COUNT(*) AS n_req FROM req
    ),
    genomes_ok AS (
        SELECT p.genomeID
        FROM Proteins p
        {join_limiter0}
        JOIN Domains d ON d.proteinID = p.proteinID
        JOIN req r     ON r.domain    = d.domain
        GROUP BY p.genomeID
        HAVING COUNT(DISTINCT r.domain) = (SELECT n_req FROM req_count)
    )
    """.format(
            join_limiter0=(
                "JOIN lim lg0 ON lg0.genomeID = p.genomeID" if use_limiter else ""
            )
        )

    # 2) Cluster ausschließen, die irgendeine Exklusionsdomäne enthalten
    if use_exclusions:
        sql += """,
    clusters_excluded AS (
        SELECT p.clusterID AS clusterID
        FROM Proteins p
        {join_limiter2}
        JOIN Domains d   ON d.proteinID = p.proteinID
        JOIN tmp_excl_domains e ON e.domain = d.domain
        WHERE p.clusterID IS NOT NULL
        GROUP BY p.clusterID
    )
    """.format(
            join_limiter2=(
                "JOIN lim lg2 ON lg2.genomeID = p.genomeID" if use_limiter else ""
            )
        )

    # 3) Finale Auswahl
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
    JOIN req r      ON r.domain    = d.domain
    JOIN Proteins p ON p.proteinID = d.proteinID
    {join_limiter3}
    {join_genomes_ok}
    {left_join_excl}
    {where_clause}
    ORDER BY p.genomeID, COALESCE(p.clusterID, -1), p.contig, p.start, d.domStart
    """

    join_limiter3 = "JOIN lim lg3 ON lg3.genomeID = p.genomeID" if use_limiter else ""
    join_genomes_ok = (
        "JOIN genomes_ok gok ON gok.genomeID = p.genomeID"
        if require_all_domains_in_same_genome
        else ""
    )
    left_join_excl = (
        "LEFT JOIN clusters_excluded x ON x.clusterID = p.clusterID"
        if use_exclusions
        else ""
    )

    where_parts = []
    if use_exclusions:
        where_parts.append("x.clusterID IS NULL")
    if not use_non_valid_hits and valid_hit_column_available:
        where_parts.append("p.valid_hit = 1")
    where_clause = f"WHERE {' AND '.join(where_parts)}" if where_parts else ""

    sql = sql.format(
        join_limiter3=join_limiter3,
        join_genomes_ok=join_genomes_ok,
        left_join_excl=left_join_excl,
        where_clause=where_clause,
    )
    return sql, []


#


def build_proteins_from_query(
        cur: sqlite3.Cursor,
        sql: str,
        args: tuple | list | None,
        protein_dict: Dict[str, parse_reports.Protein],
        genome_id_set: Set[str],
        fusion_prot_ids: Set[str] | None = None,
) -> None:
    """
    Execute SQL and build Protein objects directly from the cursor iterator.

    This avoids cur.fetchall() and processes rows as they stream from SQLite.
    """
    cur.execute(sql, args)

    for row in cur:  # streamed, no full result loaded
        protein_id = row["proteinID"]
        genome_id = row["genomeID"]
        cluster_id = row["clusterID"]

        domain = row["domain"]
        dom_start = row["domStart"]
        dom_end = row["domEnd"]
        score = row["score"]

        existing = protein_dict.get(protein_id)
        if existing is not None:
            existing.add_domain(domain, dom_start, dom_end, score)
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

        if fusion_prot_ids is not None:
            dom_count = row["dom_count"]
            if dom_count is not None and int(dom_count) >= 2:
                fusion_prot_ids.add(protein_id)


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
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_fused_ids (proteinID TEXT PRIMARY KEY);"
        )
        cur.execute("DELETE FROM tmp_fused_ids;")

        # Atomare Bulk-Inserts ohne manuelles BEGIN/COMMIT:
        # Der 'with con:' Kontext oben sorgt für Transaktion pro Block.
        cur.executemany(
            "INSERT OR IGNORE INTO tmp_fused_ids(proteinID) VALUES (?)",
            ((pid,) for pid in fusion_prot_ids),
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

        logger.info(
            f"Added {added} fused-domain rows for {len(fusion_prot_ids)} proteins (single-pass)."
        )
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

        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_tax_fetch (id TEXT PRIMARY KEY);"
        )
        cur.execute("DELETE FROM tmp_tax_fetch;")

        cur.executemany(
            "INSERT OR IGNORE INTO tmp_tax_fetch(id) VALUES (?)", ((g,) for g in wanted)
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

        logger.info(
            f"fetch_taxonomy_dict: {added} Taxonomie-Zeilen für {len(wanted)} genomeIDs hinzugefügt."
        )
    return taxon_dict

#!/usr/bin/python
from __future__ import annotations

import os
import sqlite3
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_taxonomy
from hmsss.parse_reports import parse_reports

logger = get_logger(__name__)


def _db_uri(database: str) -> str:
    """Return a read-only immutable SQLite URI."""
    if database.startswith("file:"):
        return database
    return f"file:{os.path.abspath(database)}?mode=ro&immutable=1"


# ---------------------------------------------------------------------------
# TEMP-table preparation
# ---------------------------------------------------------------------------


def _prepare_required_domains_temp(
    cur: sqlite3.Cursor,
    required_domains: Optional[Iterable[str]],
) -> int:
    """
    Store requested external domain names and resolve them once to domain_pk.

    If required_domains is empty/None, all DomainTypes are inserted. This
    preserves the previous behaviour used when all hits of selected genomes
    should be fetched.

    Unknown requested domain names remain in the table with domain_pk = NULL.
    This is intentional: an ALL-of query containing an unknown domain must
    return no candidate rather than silently ignoring that domain.
    """
    domains = sorted(
        {str(d).strip() for d in (required_domains or []) if str(d).strip()}
    )

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_req_domains (
            domain TEXT PRIMARY KEY,
            domain_pk INTEGER
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_req_domains")

    if not domains:
        cur.execute("""
            INSERT INTO tmp_req_domains(domain, domain_pk)
            SELECT domain, domain_pk FROM DomainTypes
        """)
        return cur.rowcount

    cur.executemany(
        "INSERT INTO tmp_req_domains(domain) VALUES (?)", ((d,) for d in domains)
    )
    cur.execute("""
        UPDATE tmp_req_domains
        SET domain_pk = (
            SELECT dt.domain_pk
            FROM DomainTypes dt
            WHERE dt.domain = tmp_req_domains.domain
        )
    """)
    return len(domains)


def _prepare_excluded_domains_temp(
    cur: sqlite3.Cursor,
    excluded_domains: Optional[Iterable[str]],
) -> int:
    """Store excluded external domain names and resolve them to domain_pk."""
    domains = sorted(
        {str(d).strip() for d in (excluded_domains or []) if str(d).strip()}
    )

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_excl_domains (
            domain TEXT PRIMARY KEY,
            domain_pk INTEGER
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_excl_domains")

    if not domains:
        return 0

    cur.executemany(
        "INSERT INTO tmp_excl_domains(domain) VALUES (?)", ((d,) for d in domains)
    )
    cur.execute("""
        UPDATE tmp_excl_domains
        SET domain_pk = (
            SELECT dt.domain_pk
            FROM DomainTypes dt
            WHERE dt.domain = tmp_excl_domains.domain
        )
    """)
    return len(domains)


def _prepare_limiter_genomes_temp(
    cur: sqlite3.Cursor,
    limiter_dict: Optional[Dict[str, Any]],
) -> int:
    """
    Store external genomeIDs and resolve them once to genome_pk.

    A non-empty limiter containing unknown genomeIDs still counts as an active
    limiter; unresolved rows simply have genome_pk = NULL and cannot match.
    """
    if not limiter_dict:
        return 0

    genome_ids = sorted({str(g).strip() for g in limiter_dict if str(g).strip()})
    if not genome_ids:
        return 0

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_req_genomes (
            genomeID TEXT PRIMARY KEY,
            genome_pk INTEGER
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_req_genomes")
    cur.executemany(
        "INSERT INTO tmp_req_genomes(genomeID) VALUES (?)", ((g,) for g in genome_ids)
    )
    cur.execute("""
        UPDATE tmp_req_genomes
        SET genome_pk = (
            SELECT g.genome_pk
            FROM Genomes g
            WHERE g.genomeID = tmp_req_genomes.genomeID
        )
    """)
    return len(genome_ids)


def has_column(cur: sqlite3.Cursor, table_name: str, column_name: str) -> bool:
    """Compatibility helper retained for callers/tests using the old API."""
    try:
        cur.execute(f"PRAGMA table_info({table_name})")
        return any(row[1] == column_name for row in cur.fetchall())
    except sqlite3.Error as exc:
        logger.warning(
            "Could not inspect column %s.%s: %s", table_name, column_name, exc
        )
        return False


# ---------------------------------------------------------------------------
# Query generators for single-combination fetches
# ---------------------------------------------------------------------------


def generate_fetch_query_covering_domains(
    required_domains: Iterable[str],
    use_limiter: bool = True,
    use_exclusions: bool = True,
    use_non_valid_hits: bool = True,
    valid_hit_column_available: bool = True,
) -> Tuple[str, List[Any]]:
    """
    Return all protein/domain rows from clusters containing every requested domain.

    Candidate selection uses ClusterDomains; full annotations are then read from
    Proteins/Domains. Excluded domains exclude a cluster irrespective of
    valid_hit, matching the previous -fc behaviour.
    """
    del required_domains, valid_hit_column_available

    limiter_join = (
        """
        JOIN Clusters lc ON lc.cluster_pk = cd.cluster_pk
        JOIN tmp_req_genomes lg ON lg.genome_pk = lc.genome_pk
    """
        if use_limiter
        else ""
    )

    presence_valid = "" if use_non_valid_hits else "AND cd.valid_present = 1"
    protein_valid = "" if use_non_valid_hits else "AND p.valid_hit = 1"

    exclusion_cte = ""
    exclusion_join = ""
    exclusion_where = ""
    if use_exclusions:
        exclusion_cte = """,
        clusters_excluded AS (
            SELECT DISTINCT xd.cluster_pk
            FROM ClusterDomains xd
            JOIN tmp_excl_domains e ON e.domain_pk = xd.domain_pk
            WHERE e.domain_pk IS NOT NULL
        )
        """
        exclusion_join = "LEFT JOIN clusters_excluded x ON x.cluster_pk = p.cluster_pk"
        exclusion_where = "AND x.cluster_pk IS NULL"

    sql = f"""
        WITH req AS (
            SELECT domain_pk FROM tmp_req_domains
        ),
        req_count AS (
            SELECT COUNT(*) AS n_req FROM tmp_req_domains
        ),
        clusters_covering AS (
            SELECT cd.cluster_pk
            FROM ClusterDomains cd
            {limiter_join}
            JOIN req r ON r.domain_pk = cd.domain_pk
            WHERE 1=1 {presence_valid}
            GROUP BY cd.cluster_pk
            HAVING COUNT(*) = (SELECT n_req FROM req_count)
        )
        {exclusion_cte}

        SELECT p.proteinID AS proteinID, g.genomeID AS genomeID, c.clusterID AS clusterID,
               p.contig AS contig, p.start AS gene_start, p.end AS gene_end,
               p.strand AS gene_strand, dt.domain AS domain, d.domStart AS domStart,
               d.domEnd AS domEnd, d.score AS score, p.dom_count AS dom_count,
               p.comment AS comment, p.alternative_hit AS alternative_hit
        FROM clusters_covering cc
        JOIN Proteins p ON p.cluster_pk = cc.cluster_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        JOIN Domains d ON d.protein_pk = p.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        {exclusion_join}
        WHERE 1=1 {protein_valid} {exclusion_where}
    """
    return sql, []


def generate_fetch_query_domains_anywhere(
    required_domains: Iterable[str],
    use_limiter: bool = True,
) -> Tuple[str, List[Any]]:
    """Return requested domain hits anywhere in the selected genomes."""
    del required_domains

    limiter_join = (
        "JOIN tmp_req_genomes lg ON lg.genome_pk = p.genome_pk" if use_limiter else ""
    )

    sql = f"""
        SELECT p.proteinID AS proteinID, g.genomeID AS genomeID, c.clusterID AS clusterID,
               p.contig AS contig, p.start AS gene_start, p.end AS gene_end,
               p.strand AS gene_strand, dt.domain AS domain, d.domStart AS domStart,
               d.domEnd AS domEnd, d.score AS score, p.dom_count AS dom_count,
               p.comment AS comment, p.alternative_hit AS alternative_hit
        FROM Domains d
        JOIN tmp_req_domains r ON r.domain_pk = d.domain_pk
        JOIN Proteins p ON p.protein_pk = d.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        LEFT JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        {limiter_join}
    """
    return sql, []


def generate_fetch_query_domains_anywhere_excluding_clusters(
    use_limiter: bool = False,
    use_exclusions: bool = True,
    require_all_domains_in_same_genome: bool = True,
    use_non_valid_hits: bool = True,
    valid_hit_column_available: bool = True,
) -> tuple[str, list]:
    """
    Return requested domain hits from genomes satisfying the complete requested set.

    Genome candidate selection uses GenomeDomains. Cluster exclusions preserve the
    previous -fd semantics: when valid-only mode is active, a cluster is excluded
    only when the excluded domain is present on at least one valid protein.
    """
    del valid_hit_column_available

    presence_valid = "" if use_non_valid_hits else "AND gd.valid_present = 1"
    exclusion_valid = "" if use_non_valid_hits else "AND cd.valid_present = 1"
    protein_valid = "" if use_non_valid_hits else "AND p.valid_hit = 1"

    limiter_candidate = (
        "JOIN tmp_req_genomes lg ON lg.genome_pk = gd.genome_pk" if use_limiter else ""
    )
    limiter_final = (
        "JOIN tmp_req_genomes lg2 ON lg2.genome_pk = p.genome_pk" if use_limiter else ""
    )

    genome_cte = ""
    genome_join = ""
    if require_all_domains_in_same_genome:
        genome_cte = f""",
        req_count AS (
            SELECT COUNT(*) AS n_req FROM tmp_req_domains
        ),
        genomes_ok AS (
            SELECT gd.genome_pk
            FROM GenomeDomains gd
            {limiter_candidate}
            JOIN req r ON r.domain_pk = gd.domain_pk
            WHERE 1=1 {presence_valid}
            GROUP BY gd.genome_pk
            HAVING COUNT(*) = (SELECT n_req FROM req_count)
        )
        """
        genome_join = "JOIN genomes_ok gok ON gok.genome_pk = p.genome_pk"

    exclusion_cte = ""
    exclusion_join = ""
    exclusion_where = ""
    if use_exclusions:
        exclusion_cte = f""",
        clusters_excluded AS (
            SELECT DISTINCT cd.cluster_pk
            FROM ClusterDomains cd
            JOIN tmp_excl_domains e ON e.domain_pk = cd.domain_pk
            WHERE e.domain_pk IS NOT NULL {exclusion_valid}
        )
        """
        exclusion_join = "LEFT JOIN clusters_excluded x ON x.cluster_pk = p.cluster_pk"
        exclusion_where = "AND x.cluster_pk IS NULL"

    sql = f"""
        WITH req AS (
            SELECT domain_pk FROM tmp_req_domains
        )
        {genome_cte}
        {exclusion_cte}

        SELECT p.proteinID AS proteinID, g.genomeID AS genomeID, c.clusterID AS clusterID,
               p.contig AS contig, p.start AS gene_start, p.end AS gene_end,
               p.strand AS gene_strand, dt.domain AS domain, d.domStart AS domStart,
               d.domEnd AS domEnd, d.score AS score, p.dom_count AS dom_count,
               p.comment AS comment, p.alternative_hit AS alternative_hit
        FROM Domains d
        JOIN req r ON r.domain_pk = d.domain_pk
        JOIN Proteins p ON p.protein_pk = d.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        LEFT JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        {limiter_final}
        {genome_join}
        {exclusion_join}
        WHERE 1=1 {protein_valid} {exclusion_where}
        ORDER BY p.genome_pk, p.cluster_pk, p.contig, p.start, d.domStart
    """
    return sql, []


# ---------------------------------------------------------------------------
# Protein-object hydration
# ---------------------------------------------------------------------------


def build_proteins_from_query(
    cur: sqlite3.Cursor,
    sql: str,
    args: tuple | list | None,
    protein_dict: Dict[str, parse_reports.Protein],
    genome_id_set: Set[str],
    fusion_prot_ids: Optional[Set[str]] = None,
) -> None:
    """Execute a query and stream its rows into Protein objects."""
    cur.execute(sql, args or [])

    for row in cur:
        protein_id = row["proteinID"]
        existing = protein_dict.get(protein_id)

        if existing is not None:
            existing.add_domain(
                row["domain"], row["domStart"], row["domEnd"], row["score"]
            )
        else:
            protein = parse_reports.Protein(
                protein_id, row["domain"], row["domStart"], row["domEnd"], row["score"]
            )
            protein.genomeID = row["genomeID"]
            protein.clusterID = row["clusterID"] or ""
            protein.gene_contig = row["contig"] or ""
            protein.gene_start = row["gene_start"] or 0
            protein.gene_end = row["gene_end"] or 0
            protein.gene_strand = row["gene_strand"] or "."
            protein.protein_sequence = ""
            protein.selection_comment = row["comment"] or ""
            protein.alternative_hit = row["alternative_hit"] or ""
            protein_dict[protein_id] = protein

            if row["genomeID"]:
                genome_id_set.add(row["genomeID"])

        if (
            fusion_prot_ids is not None
            and row["dom_count"] is not None
            and int(row["dom_count"]) >= 2
        ):
            fusion_prot_ids.add(protein_id)


def hydrate_protein_sequences(cur: sqlite3.Cursor, protein_dict: Dict[str, Any]) -> int:
    """Load each selected protein sequence once, after hit selection is complete."""
    if not protein_dict:
        return 0

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_sequence_protein_ids (
            proteinID TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_sequence_protein_ids")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_sequence_protein_ids(proteinID) VALUES (?)",
        ((pid,) for pid in protein_dict),
    )

    cur.execute("""
        SELECT p.proteinID, s.sequence
        FROM tmp_sequence_protein_ids t
        JOIN Proteins p ON p.proteinID = t.proteinID
        JOIN ProteinSequences s ON s.protein_pk = p.protein_pk
    """)

    added = 0
    for row in cur:
        protein = protein_dict.get(row["proteinID"])
        if protein is not None:
            protein.protein_sequence = row["sequence"] or ""
            added += 1

    logger.debug("Hydrated sequences for %d proteins.", added)
    return added


def hydrate_fused_protein_domains_from_cursor(
    cur: sqlite3.Cursor,
    fusion_prot_ids: Iterable[str],
    protein_dict: Dict[str, Any],
) -> int:
    """
    Add all domains for selected multi-domain proteins using integer joins after
    resolving external proteinIDs once.
    """
    fusion_ids = {str(pid) for pid in fusion_prot_ids if pid}
    if not fusion_ids:
        return 0

    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_fused_protein_ids (
            proteinID TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_fused_protein_ids")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_fused_protein_ids(proteinID) VALUES (?)",
        ((pid,) for pid in fusion_ids),
    )

    cur.execute("""
        SELECT p.proteinID, dt.domain, d.domStart, d.domEnd, d.score
        FROM tmp_fused_protein_ids t
        JOIN Proteins p ON p.proteinID = t.proteinID
        JOIN Domains d ON d.protein_pk = p.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
    """)

    added = 0
    for row in cur:
        protein = protein_dict.get(row["proteinID"])
        if protein is not None:
            protein.add_domain(
                row["domain"], row["domStart"], row["domEnd"], row["score"]
            )
            added += 1

    logger.debug(
        "Hydrated %d fused-domain rows for %d proteins.", added, len(fusion_ids)
    )
    return added


def hydrate_fused_protein_domains(
    db_path: str,
    fusion_prot_ids: Iterable[str],
    protein_dict: Dict[str, Any],
) -> int:
    """Backward-compatible wrapper opening its own read-only connection."""
    if not fusion_prot_ids:
        return 0

    with sqlite3.connect(_db_uri(db_path), uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()
        cur.execute("PRAGMA temp_store = MEMORY")
        return hydrate_fused_protein_domains_from_cursor(
            cur, fusion_prot_ids, protein_dict
        )


# ---------------------------------------------------------------------------
# Standard single-combination fetch
# ---------------------------------------------------------------------------


def fetch_bulk_data(
    database: str,
    syntenic_domains: Optional[List[str]],
    limiter_dict: Optional[Dict[str, Any]] = None,
    fetch_from_gene_clusters: bool = False,
    excluded_domains: Optional[List[str]] = None,
    use_non_valid_hits: bool = False,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, dict[str, str]]]:
    """Fetch proteins for one requested domain combination."""
    protein_dict: Dict[str, Any] = {}
    cluster_dict: Dict[str, Any] = {}
    genome_id_set: Set[str] = set()
    fusion_prot_ids: Set[str] = set()

    db_path = _db_uri(database)

    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -262144")

        _prepare_required_domains_temp(cur, syntenic_domains)
        n_excluded = _prepare_excluded_domains_temp(cur, excluded_domains)
        n_limiter = _prepare_limiter_genomes_temp(cur, limiter_dict)

        if fetch_from_gene_clusters:
            logger.info("Searching for syntenic %s", syntenic_domains)
            sql, args = generate_fetch_query_covering_domains(
                syntenic_domains or [],
                use_limiter=n_limiter > 0,
                use_exclusions=n_excluded > 0,
                use_non_valid_hits=use_non_valid_hits,
            )
        else:
            logger.info("Searching for co-occurring %s", syntenic_domains)
            sql, args = generate_fetch_query_domains_anywhere_excluding_clusters(
                use_limiter=n_limiter > 0,
                use_exclusions=n_excluded > 0,
                require_all_domains_in_same_genome=bool(syntenic_domains),
                use_non_valid_hits=use_non_valid_hits,
            )

        build_proteins_from_query(
            cur, sql, args, protein_dict, genome_id_set, fusion_prot_ids
        )

        if fusion_prot_ids and not fetch_from_gene_clusters:
            hydrate_fused_protein_domains_from_cursor(
                cur, fusion_prot_ids, protein_dict
            )

        hydrate_protein_sequences(cur, protein_dict)

    taxon_dict = db_fetch_taxonomy.fetch_taxonomy_dict(
        db_path, genome_ids=genome_id_set
    )
    parse_reports.define_best_score_hits_for_protein_dict(protein_dict)

    logger.info(
        "Fetched %d proteins from %d genomes.", len(protein_dict), len(genome_id_set)
    )
    return protein_dict, cluster_dict, taxon_dict


# ---------------------------------------------------------------------------
# Multi-combination candidate collection
# ---------------------------------------------------------------------------


def _prepare_candidate_tables(cur: sqlite3.Cursor) -> None:
    """
    Prepare TEMP tables used while evaluating OR/optional combinations.

    For -fd we store candidate (genome_pk, domain_pk) pairs rather than one
    global domain union. This preserves the exact union of the former
    combination-by-combination fetches even for disjoint combinations.
    """
    ddl = (
        "CREATE TEMP TABLE IF NOT EXISTS tmp_combo_genomes (genome_pk INTEGER PRIMARY KEY)",
        "CREATE TEMP TABLE IF NOT EXISTS tmp_combo_output_genomes (genome_pk INTEGER PRIMARY KEY)",
        "CREATE TEMP TABLE IF NOT EXISTS tmp_combo_clusters (cluster_pk INTEGER PRIMARY KEY)",
        """CREATE TEMP TABLE IF NOT EXISTS tmp_candidate_genome_domains (
               genome_pk INTEGER NOT NULL,
               domain_pk INTEGER NOT NULL,
               PRIMARY KEY (genome_pk, domain_pk)
           ) WITHOUT ROWID""",
        "CREATE TEMP TABLE IF NOT EXISTS tmp_candidate_clusters (cluster_pk INTEGER PRIMARY KEY)",
    )
    for sql in ddl:
        cur.execute(sql)

    for table in (
        "tmp_combo_genomes",
        "tmp_combo_output_genomes",
        "tmp_combo_clusters",
        "tmp_candidate_genome_domains",
        "tmp_candidate_clusters",
    ):
        cur.execute(f"DELETE FROM {table}")


def _collect_genome_candidates_for_combo(
    cur: sqlite3.Cursor,
    combo: Iterable[str],
    *,
    use_limiter: bool,
    use_exclusions: bool,
    use_non_valid_hits: bool,
) -> Set[str]:
    """
    Evaluate one -fd combination without hydrating Protein objects.

    The ALL-of candidate check happens on GenomeDomains. Candidate domain pairs
    are retained only for domains belonging to combinations that the genome
    actually satisfies.
    """
    _prepare_required_domains_temp(cur, combo)
    cur.execute("DELETE FROM tmp_combo_genomes")
    cur.execute("DELETE FROM tmp_combo_output_genomes")

    limiter_join = (
        "JOIN tmp_req_genomes lim ON lim.genome_pk = gd.genome_pk"
        if use_limiter
        else ""
    )
    presence_valid = "" if use_non_valid_hits else "AND gd.valid_present = 1"

    cur.execute(f"""
        INSERT OR IGNORE INTO tmp_combo_genomes(genome_pk)
        SELECT gd.genome_pk
        FROM GenomeDomains gd
        {limiter_join}
        JOIN tmp_req_domains r ON r.domain_pk = gd.domain_pk
        WHERE 1=1 {presence_valid}
        GROUP BY gd.genome_pk
        HAVING COUNT(*) = (SELECT COUNT(*) FROM tmp_req_domains)
    """)

    if not use_exclusions:
        cur.execute("""
            INSERT OR IGNORE INTO tmp_combo_output_genomes(genome_pk)
            SELECT genome_pk FROM tmp_combo_genomes
        """)
        cur.execute("""
            INSERT OR IGNORE INTO tmp_candidate_genome_domains(genome_pk, domain_pk)
            SELECT cg.genome_pk, r.domain_pk
            FROM tmp_combo_genomes cg
            CROSS JOIN tmp_req_domains r
            WHERE r.domain_pk IS NOT NULL
        """)
    else:
        protein_valid = "" if use_non_valid_hits else "AND p.valid_hit = 1"
        exclusion_valid = "" if use_non_valid_hits else "AND xd.valid_present = 1"

        eligible_rows = f"""
            FROM tmp_combo_genomes cg
            JOIN Proteins p ON p.genome_pk = cg.genome_pk
            JOIN Domains d ON d.protein_pk = p.protein_pk
            JOIN tmp_req_domains r ON r.domain_pk = d.domain_pk
            WHERE r.domain_pk IS NOT NULL
              {protein_valid}
              AND (
                  p.cluster_pk IS NULL
                  OR NOT EXISTS (
                      SELECT 1
                      FROM ClusterDomains xd
                      JOIN tmp_excl_domains e ON e.domain_pk = xd.domain_pk
                      WHERE xd.cluster_pk = p.cluster_pk {exclusion_valid}
                  )
              )
        """

        cur.execute(f"""
            INSERT OR IGNORE INTO tmp_combo_output_genomes(genome_pk)
            SELECT DISTINCT cg.genome_pk
            {eligible_rows}
        """)
        cur.execute(f"""
            INSERT OR IGNORE INTO tmp_candidate_genome_domains(genome_pk, domain_pk)
            SELECT DISTINCT cg.genome_pk, d.domain_pk
            {eligible_rows}
        """)

    cur.execute("""
        SELECT g.genomeID
        FROM tmp_combo_output_genomes t
        JOIN Genomes g ON g.genome_pk = t.genome_pk
    """)
    return {row[0] for row in cur}


def _collect_cluster_candidates_for_combo(
    cur: sqlite3.Cursor,
    combo: Iterable[str],
    *,
    use_limiter: bool,
    use_exclusions: bool,
    use_non_valid_hits: bool,
) -> Set[str]:
    """Evaluate one -fc combination using ClusterDomains only."""
    _prepare_required_domains_temp(cur, combo)
    cur.execute("DELETE FROM tmp_combo_clusters")

    limiter_join = (
        """
        JOIN Clusters lc ON lc.cluster_pk = cd.cluster_pk
        JOIN tmp_req_genomes lim ON lim.genome_pk = lc.genome_pk
    """
        if use_limiter
        else ""
    )

    presence_valid = "" if use_non_valid_hits else "AND cd.valid_present = 1"

    exclusion = ""
    if use_exclusions:
        exclusion = """
            AND NOT EXISTS (
                SELECT 1
                FROM ClusterDomains xd
                JOIN tmp_excl_domains e ON e.domain_pk = xd.domain_pk
                WHERE xd.cluster_pk = cd.cluster_pk
            )
        """

    cur.execute(f"""
        INSERT OR IGNORE INTO tmp_combo_clusters(cluster_pk)
        SELECT cd.cluster_pk
        FROM ClusterDomains cd
        {limiter_join}
        JOIN tmp_req_domains r ON r.domain_pk = cd.domain_pk
        WHERE 1=1 {presence_valid} {exclusion}
        GROUP BY cd.cluster_pk
        HAVING COUNT(*) = (SELECT COUNT(*) FROM tmp_req_domains)
    """)

    cur.execute("""
        INSERT OR IGNORE INTO tmp_candidate_clusters(cluster_pk)
        SELECT cluster_pk FROM tmp_combo_clusters
    """)

    cur.execute("""
        SELECT DISTINCT g.genomeID
        FROM tmp_combo_clusters t
        JOIN Clusters c ON c.cluster_pk = t.cluster_pk
        JOIN Genomes g ON g.genome_pk = c.genome_pk
    """)
    return {row[0] for row in cur}


def _fetch_candidate_clusters(
    cur: sqlite3.Cursor,
    *,
    protein_dict: Dict[str, Any],
    genome_id_set: Set[str],
    use_non_valid_hits: bool,
) -> None:
    """Hydrate all protein/domain rows belonging to the union of candidate clusters."""
    valid_where = "" if use_non_valid_hits else "WHERE p.valid_hit = 1"

    sql = f"""
        SELECT p.proteinID AS proteinID, g.genomeID AS genomeID, c.clusterID AS clusterID,
               p.contig AS contig, p.start AS gene_start, p.end AS gene_end,
               p.strand AS gene_strand, dt.domain AS domain, d.domStart AS domStart,
               d.domEnd AS domEnd, d.score AS score, p.dom_count AS dom_count,
               p.comment AS comment, p.alternative_hit AS alternative_hit
        FROM tmp_candidate_clusters x
        JOIN Proteins p ON p.cluster_pk = x.cluster_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        JOIN Domains d ON d.protein_pk = p.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        {valid_where}
        ORDER BY p.genome_pk, p.cluster_pk, p.contig, p.start, d.domStart
    """
    build_proteins_from_query(
        cur, sql, [], protein_dict, genome_id_set, fusion_prot_ids=None
    )


def _fetch_candidate_genomes(
    cur: sqlite3.Cursor,
    *,
    protein_dict: Dict[str, Any],
    genome_id_set: Set[str],
    fusion_prot_ids: Set[str],
    use_exclusions: bool,
    use_non_valid_hits: bool,
) -> None:
    """
    Hydrate the exact union of requested domain hits for qualifying genomes.

    tmp_candidate_genome_domains prevents over-fetching domains that belong to
    another disjoint OR combination not satisfied by the same genome.
    """
    protein_valid = "" if use_non_valid_hits else "AND p.valid_hit = 1"
    exclusion = ""

    if use_exclusions:
        exclusion_valid = "" if use_non_valid_hits else "AND xd.valid_present = 1"
        exclusion = f"""
            AND (
                p.cluster_pk IS NULL
                OR NOT EXISTS (
                    SELECT 1
                    FROM ClusterDomains xd
                    JOIN tmp_excl_domains e ON e.domain_pk = xd.domain_pk
                    WHERE xd.cluster_pk = p.cluster_pk {exclusion_valid}
                )
            )
        """

    sql = f"""
        SELECT p.proteinID AS proteinID, g.genomeID AS genomeID, c.clusterID AS clusterID,
               p.contig AS contig, p.start AS gene_start, p.end AS gene_end,
               p.strand AS gene_strand, dt.domain AS domain, d.domStart AS domStart,
               d.domEnd AS domEnd, d.score AS score, p.dom_count AS dom_count,
               p.comment AS comment, p.alternative_hit AS alternative_hit
        FROM tmp_candidate_genome_domains cg
        JOIN Proteins p ON p.genome_pk = cg.genome_pk
        JOIN Domains d ON d.protein_pk = p.protein_pk AND d.domain_pk = cg.domain_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        LEFT JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        WHERE 1=1 {protein_valid} {exclusion}
        ORDER BY p.genome_pk, p.cluster_pk, p.contig, p.start, d.domStart
    """
    build_proteins_from_query(
        cur, sql, [], protein_dict, genome_id_set, fusion_prot_ids
    )


def fetch_bulk_data_for_combinations(
    database: str,
    combinations: Iterable[Iterable[str]],
    limiter_dict: Optional[Dict[str, Any]] = None,
    fetch_from_gene_clusters: bool = False,
    excluded_domains: Optional[Iterable[str]] = None,
    use_non_valid_hits: bool = False,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[tuple[str, ...], set[str]],
]:
    """
    Evaluate many OR/optional combinations cheaply, then hydrate the union once.

    Returns the same external structures used by the former loop of repeated
    fetch_bulk_data() calls, including combo -> genomeID mapping.
    """
    combos = [list(combo) for combo in combinations if combo]
    if not combos:
        return {}, {}, {}, {}

    protein_dict: Dict[str, Any] = {}
    cluster_dict: Dict[str, Any] = {}
    genome_id_set: Set[str] = set()
    fusion_prot_ids: Set[str] = set()
    combo_to_genomes: Dict[tuple[str, ...], set[str]] = {}

    db_path = _db_uri(database)

    with sqlite3.connect(db_path, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -262144")

        n_limiter = _prepare_limiter_genomes_temp(cur, limiter_dict)
        n_excluded = _prepare_excluded_domains_temp(cur, excluded_domains)
        _prepare_candidate_tables(cur)

        for i, combo in enumerate(combos, start=1):
            logger.info(
                "Searching combination %d/%d: %s",
                i,
                len(combos),
                ", ".join(combo),
            )

            if fetch_from_gene_clusters:
                genomes = _collect_cluster_candidates_for_combo(
                    cur,
                    combo,
                    use_limiter=n_limiter > 0,
                    use_exclusions=n_excluded > 0,
                    use_non_valid_hits=use_non_valid_hits,
                )
            else:
                genomes = _collect_genome_candidates_for_combo(
                    cur,
                    combo,
                    use_limiter=n_limiter > 0,
                    use_exclusions=n_excluded > 0,
                    use_non_valid_hits=use_non_valid_hits,
                )

            combo_to_genomes[tuple(combo)] = genomes

        if fetch_from_gene_clusters:
            _fetch_candidate_clusters(
                cur,
                protein_dict=protein_dict,
                genome_id_set=genome_id_set,
                use_non_valid_hits=use_non_valid_hits,
            )
        else:
            _fetch_candidate_genomes(
                cur,
                protein_dict=protein_dict,
                genome_id_set=genome_id_set,
                fusion_prot_ids=fusion_prot_ids,
                use_exclusions=n_excluded > 0,
                use_non_valid_hits=use_non_valid_hits,
            )

            if fusion_prot_ids:
                hydrate_fused_protein_domains_from_cursor(
                    cur, fusion_prot_ids, protein_dict
                )

        hydrate_protein_sequences(cur, protein_dict)

    taxon_dict = db_fetch_taxonomy.fetch_taxonomy_dict(
        db_path, genome_ids=genome_id_set
    )
    parse_reports.define_best_score_hits_for_protein_dict(protein_dict)

    logger.info(
        "Fetched %d proteins from %d genomes across %d combinations.",
        len(protein_dict),
        len(genome_id_set),
        len(combos),
    )
    return protein_dict, cluster_dict, taxon_dict, combo_to_genomes

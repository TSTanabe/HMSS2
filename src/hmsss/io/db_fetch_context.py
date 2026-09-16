import os
import sqlite3
from typing import Any, Dict, Set, Tuple

from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_protein

logger = get_logger(__name__)


def _collect_cluster_ids_from_proteins(protein_dict: Dict[str, Any]) -> Set[str]:
    """Collect unique non-empty external clusterIDs from Protein objects."""
    return {
        str(p.clusterID) for p in protein_dict.values() if getattr(p, "clusterID", None)
    }


def _prepare_cluster_context_query(
    cur: sqlite3.Cursor, cluster_ids: Set[str]
) -> Tuple[str, Tuple]:
    """
    Resolve external clusterIDs once to internal cluster_pk values and return
    the query for retrieving the complete protein/domain context.
    """

    # External identifiers enter the DB interface here.
    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_context_cluster_ids (
            clusterID TEXT PRIMARY KEY
        ) WITHOUT ROWID
    """)
    cur.execute("DELETE FROM tmp_context_cluster_ids")
    cur.executemany(
        "INSERT OR IGNORE INTO tmp_context_cluster_ids(clusterID) VALUES (?)",
        ((cid,) for cid in cluster_ids),
    )

    # Resolve TEXT clusterIDs once. All subsequent large joins use INTEGER keys.
    cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS tmp_context_clusters (
            cluster_pk INTEGER PRIMARY KEY
        )
    """)
    cur.execute("DELETE FROM tmp_context_clusters")
    cur.execute("""
        INSERT OR IGNORE INTO tmp_context_clusters(cluster_pk)
        SELECT c.cluster_pk
        FROM Clusters c
        JOIN tmp_context_cluster_ids t ON t.clusterID = c.clusterID
    """)

    sql = """
        SELECT
            p.proteinID AS proteinID,
            g.genomeID AS genomeID,
            c.clusterID AS clusterID,
            p.contig AS contig,
            p.start AS gene_start,
            p.end AS gene_end,
            p.strand AS gene_strand,
            dt.domain AS domain,
            d.domStart AS domStart,
            d.domEnd AS domEnd,
            d.score AS score,
            p.dom_count AS dom_count,
            p.comment AS comment,
            p.alternative_hit AS alternative_hit
        FROM tmp_context_clusters t
        JOIN Proteins p ON p.cluster_pk = t.cluster_pk
        JOIN Genomes g ON g.genome_pk = p.genome_pk
        JOIN Clusters c ON c.cluster_pk = p.cluster_pk
        JOIN Domains d ON d.protein_pk = p.protein_pk
        JOIN DomainTypes dt ON dt.domain_pk = d.domain_pk
        ORDER BY g.genomeID, c.clusterID, p.contig, p.start, d.domStart
    """

    return sql, ()


def fetch_cluster_context_for_proteins(
    database: str,
    base_protein_dict: Dict[str, Any],
    fetch_from_gene_clusters: bool,
) -> Dict[str, Any]:
    """
    Retrieve complete cluster context.

    - -fc: base fetch already contains complete clusters, so return unchanged.
    - -fd: collect clusterIDs from initial hits and fetch all proteins/domains
      belonging to these clusters.
    """

    if fetch_from_gene_clusters:
        logger.info("Cluster context: -fc result already contains complete clusters.")
        return base_protein_dict

    if not base_protein_dict:
        logger.info("Cluster context: base protein dictionary is empty.")
        return {}

    cluster_ids = _collect_cluster_ids_from_proteins(base_protein_dict)
    if not cluster_ids:
        logger.info("Cluster context: no clusterIDs present in initial proteins.")
        return {}

    logger.info(
        "Cluster context: fetching complete context for %d clusters.", len(cluster_ids)
    )

    abs_db = os.path.abspath(database)
    db_uri = f"file:{abs_db}?mode=ro&immutable=1"

    context_protein_dict: Dict[str, Any] = {}
    genome_id_set: Set[str] = set()

    with sqlite3.connect(db_uri, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")
        cur.execute("PRAGMA cache_size = -262144")  # ~256 MiB

        sql, args = _prepare_cluster_context_query(cur, cluster_ids)

        db_fetch_protein.build_proteins_from_query(
            cur=cur,
            sql=sql,
            args=args,
            protein_dict=context_protein_dict,
            genome_id_set=genome_id_set,
            fusion_prot_ids=None,
        )

        db_fetch_protein.hydrate_protein_sequences(cur, context_protein_dict)

    logger.info(
        "Cluster context: fetched %d proteins across %d requested clusters from %d genomes.",
        len(context_protein_dict),
        len(cluster_ids),
        len(genome_id_set),
    )

    return context_protein_dict

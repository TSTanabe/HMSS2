import os
import sqlite3
from typing import Dict, Any, Set, Tuple

from hmsss.core.logging import get_logger
from hmsss.io import db_fetch_protein

logger = get_logger(__name__)


def _collect_cluster_ids_from_proteins(
    protein_dict: Dict[str, Any],
) -> Set[str]:
    """
    Collect all unique, non-empty clusterIDs from a protein dictionary.

    Args:
        protein_dict: Mapping proteinID → Protein object

    Returns:
        Set of clusterID strings
    """
    cluster_ids: Set[str] = set()
    for p in protein_dict.values():
        cid = getattr(p, "clusterID", None)
        if cid:
            cluster_ids.add(str(cid))
    return cluster_ids


def _prepare_cluster_context_query(
    cur: sqlite3.Cursor,
    cluster_ids: Set[str],
) -> Tuple[str, Tuple]:
    """
    Prepare the SQL query to fetch full cluster context.

    - Creates and populates a TEMP table with the given clusterIDs.
    - Returns the SELECT statement and its args tuple.

    This keeps SQL construction and temp-table handling in one place.
    """
    # TEMP table for relevant clusterIDs
    cur.execute(
        """
        CREATE TEMP TABLE IF NOT EXISTS tmp_context_clusterIDs (
            clusterID TEXT PRIMARY KEY
        );
        """
    )
    cur.execute("DELETE FROM tmp_context_clusterIDs;")

    cur.executemany(
        "INSERT OR IGNORE INTO tmp_context_clusterIDs(clusterID) VALUES (?);",
        ((cid,) for cid in cluster_ids),
    )

    sql = """
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
        JOIN tmp_context_clusterIDs t
          ON t.clusterID = p.clusterID
        LEFT JOIN Domains d
          ON d.proteinID = p.proteinID
        ORDER BY
            p.genomeID,
            COALESCE(p.clusterID, ''),
            p.contig,
            p.start,
            d.domStart
    """

    # no positional args needed here, but we keep the interface generic
    return sql, ()


def fetch_cluster_context_for_proteins(
    database: str,
    base_protein_dict: Dict[str, Any],
    fetch_from_gene_clusters: bool,
) -> Dict[str, Any]:
    """
    Retrieve the full gene-cluster context for proteins.

    - If fetch_from_gene_clusters=True (-fc mode):
        Clusters are already complete -> return base_protein_dict unchanged.
    - If fetch_from_gene_clusters=False (-fd mode):
        1. Collect clusterIDs from base_protein_dict.
        2. Build a cluster-context query via _prepare_cluster_context_query.
        3. Stream rows into Protein objects using build_proteins_from_query.
    """
    if fetch_from_gene_clusters:
        logger.info(
            "Cluster context: fetch_from_gene_clusters=True, using existing clusters."
        )
        return base_protein_dict

    if not base_protein_dict:
        logger.info(
            "Cluster context: base_protein_dict is empty; skipping context fetch."
        )
        return {}

    cluster_ids = _collect_cluster_ids_from_proteins(base_protein_dict)
    if not cluster_ids:
        logger.info(
            "Cluster context: no clusterIDs found in base_protein_dict; "
            "no cluster context to fetch."
        )
        return {}

    logger.info(
        "Cluster context: fetching full clusters for %d clusterIDs.",
        len(cluster_ids),
    )

    abs_db = os.path.abspath(database)
    db_uri = f"file:{abs_db}?mode=ro&immutable=1"

    context_protein_dict: Dict[str, Any] = {}
    genome_id_set: Set[str] = set()
    fusion_prot_ids: Set[str] = set()

    with sqlite3.connect(db_uri, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # tuned for read-heavy workloads
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA temp_store = MEMORY;")

        sql, args = _prepare_cluster_context_query(cur, cluster_ids)

        db_fetch_protein.build_proteins_from_query(
            cur=cur,
            sql=sql,
            args=args,
            protein_dict=context_protein_dict,
            genome_id_set=genome_id_set,
            fusion_prot_ids=fusion_prot_ids,
        )

    if fusion_prot_ids:
        db_fetch_protein.hydrate_fused_protein_domains(
            db_path=abs_db,
            fusion_prot_ids=fusion_prot_ids,
            protein_dict=context_protein_dict,
        )

    logger.info(
        "Cluster context: final context dict has %d proteins across %d clusters "
        "(base hits: %d).",
        len(context_protein_dict),
        len(cluster_ids),
        len(base_protein_dict),
    )

    return context_protein_dict

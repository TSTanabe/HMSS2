from __future__ import annotations

from hmsss.cli.config import Config
from hmsss.db import database as database
from hmsss.algorithms import csb_cluster as csb_cluster
from hmsss.core.logging import get_logger, print_header

logger = get_logger(__name__)

"""
Stage: Collinear syntenic block (CSB) detection.

Predicts CSBs, clusters them by Jaccard similarity, and updates the database
with new CSB keywords, replacing old ones.
"""

def csb_finder(config: Config) -> None:
    """Run CSB prediction and update the database.

    Steps:
      - Predict collinear syntenic blocks via `csb_cluster.csb_prediction`.
      - Cluster CSBs using Jaccard distance.
      - Re-index database.
      - Delete old CSB keywords and insert new ones.

    Args:
        config: Configuration with database path and clustering parameters.
    """
    csb_gene_cluster_dict = {}
    print_header("CSB finder", logger=logger)

    logger.info("Running collinear syntenic block pattern prediction")
    try:
        csb_instances = csb_cluster.csb_prediction(config)

        csb_gene_cluster_dict = csb_cluster.csb_jaccard(
            config, computed_instances_dict=csb_instances, jaccard_distance=0.0
        )  # 0.0: nur Clusterdict bilden
    except Exception as err:
        logger.error(f"CSB finder failed: \n {err}", logger=logger)

    database.index_database(config.database_directory)
    database.delete_keywords_from_csb(
        config.database_directory
    )  # alte Schlüssel entfernen
    database.update_keywords(
        config.database_directory, csb_gene_cluster_dict
    )  # neue Schlüssel eintragen

from __future__ import annotations

from hmsss.core.config import Config
from hmsss.db import database as database
from hmsss.algorithms import csb_cluster as csb_cluster
from hmsss.core.logging import get_logger, print_header

log = get_logger(__name__)


def csb_finder(config: Config) -> None:
    """
    Aus __main__.py:
    - CSB-Vorhersage & Jaccard-Clustering
    - DB-Index aktualisieren
    - alte CSB-Keywords löschen, neue einspielen
    """
    print_header("CSB finder", logger=log)

    log.info("Running collinear syntenic block pattern prediction")
    csb_instances = csb_cluster.csb_prediction(config)

    csb_gene_cluster_dict = csb_cluster.csb_jaccard(
        config, computed_instances_dict=csb_instances, jaccard_distance= 0.0
    )  # 0.0: nur Clusterdict bilden

    database.index_database(config.database_directory)
    database.delete_keywords_from_csb(
        config.database_directory, config
    )  # alte Schlüssel entfernen
    database.update_keywords(
        config.database_directory, csb_gene_cluster_dict
    )  # neue Schlüssel eintragen

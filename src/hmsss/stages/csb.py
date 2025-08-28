from __future__ import annotations

from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.db import database as database

from hmsss.algorithms import csb_cluster as csb_cluster  # bevorzugt, wenn vorhanden

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def csb_finder(options: Hmsss) -> None:
    """
    Aus __main__.py:
    - CSB-Vorhersage & Jaccard-Clustering
    - DB-Index aktualisieren
    - alte CSB-Keywords löschen, neue einspielen
    """
    print_header("CSB finder", logger=log)

    log.info("Running collinear syntenic block pattern prediction")
    csb_instances = csb_cluster.csb_prediction(options)

    csb_gene_cluster_dict = csb_cluster.csb_jaccard(
        options, computed_Instances_dict=csb_instances, jaccard_distance= 0.0
    )  # 0.0: nur Clusterdict bilden

    database.index_database(options.database_directory)
    database.delete_keywords_from_csb(
        options.database_directory, options
    )  # alte Schlüssel entfernen
    database.update_keywords(
        options.database_directory, csb_gene_cluster_dict
    )  # neue Schlüssel eintragen

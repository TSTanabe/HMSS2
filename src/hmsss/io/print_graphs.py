#!/usr/bin/python
import os
from typing import Dict, Any, Optional, Set, Tuple, List

from hmsss.graphics import graph_presence_absence
from hmsss.graphics import graph_occurence_network
from hmsss.graphics import graph_gene_cluster
from hmsss.graphics import graph_strain_variability
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def print_hit_graphs(
    directory: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    metabolic_dict: Dict[str, Any],
    context_dict: Dict[str, Any] | None,
    fetch_proteins: List[str],
    levels: List[str],
) -> None:
    """
    Main output routine: creates hit tables, taxonomy summaries, and protein FASTA files.

    Args:
        metabolic_dict:
        fetch_proteins: input proteins from the CLI
        directory (str): Output directory path.
        protein_dict: Mapping proteinID -> proteinObj.
        cluster_dict: Mapping clusterID -> clusterObj.
        taxon_dict (Dict[str, str]): Mapping genomeID -> taxonomy string.
    """

    # Printing results
    logger.info("Printing genome information output files")
    # Metadata for taxonomy and hits
    hit_report = os.path.join(
        directory, "summary_hit_table.jpg"
    )  # individual hit table in tsv file
    gene_taxonomy = os.path.join(directory, "summary_gene_taxonomy.jpg")
    unique_file = os.path.join(directory, "summary_unique_lineages.jpg")
    taxonomy_summary = os.path.join(directory, "summary_hit_taxonomy_counts.svg")
    taxonomy_summary2 = os.path.join(
        directory, "summary_requested_hit_taxonomy_counts.svg"
    )
    metabolic_annotation = os.path.join(directory, "summary_metabolic_annotations.jpg")
    cluster_overview_report = os.path.join(
        directory, "summary_genecluster_overview_table.jpg"
    )
    path_strain_variability_tsv = os.path.join(
        directory, "summary_strain_variability_by_taxonomy.txt"
    )
    path_strain_variability_plot = os.path.join(
        directory, "summary_strain_variability_by_taxonomy"
    )

    # Plots the network for presence absence
    graph_occurence_network.plot_taxonomy_cooccurrence_network(
        taxonomy_summary,
        protein_dict,
        taxon_dict,
        allowed_types=None,  # oder Liste von Domain-Namen
        allowed_levels={"Phylum"},
        preferred_order=None,
        min_cooccurrence=3,  # kannst du anpassen
    )

    # Plot includes all proteins that were fetched extended the requested ones
    graph_presence_absence.plot_taxonomy_summary_bubbles(
        taxonomy_summary,
        protein_dict,
        taxon_dict,
        allowed_levels=set(levels),
        preferred_order=fetch_proteins,
    )

    # Plots only for the requested proteins presence absence matrix
    graph_presence_absence.plot_taxonomy_summary_bubbles(
        taxonomy_summary2,
        protein_dict,
        taxon_dict,
        allowed_types=fetch_proteins,
        allowed_levels=set(levels),
        preferred_order=fetch_proteins,
    )

    graph_strain_variability.plot_taxonomy_stacked_bars(
        input_tsv=path_strain_variability_tsv,
        outdir=path_strain_variability_plot,
        top_n=15,
        min_category_fraction=0.05,
    )

    # RAM usage is too high.
    # graph_gene_cluster.plot_gene_cluster_summary(hit_report, protein_dict, taxon_dict,allowed_types=fetch_proteins,allowed_levels=set(levels))

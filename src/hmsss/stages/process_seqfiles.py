from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import processing as processing
from hmsss.io import datasets as datasets

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def process_operator(options: Hmsss) -> None:
    """
    Aus __main__.py:
      - merge_fasta
      - concat_alignment
      - filter_fasta (min/max)
      - add_taxonomy_to_alignment
      - add_genomic_context
      - create_gene_cluster_dataset
      - create_type_range_dataset
    """
    print_header("Processing operators", logger=log)

    # Merge fasta files
    if options.merge_fasta:
        processing.merge_fasta(options)

    # Concat alignment files
    if options.concat_alignment:
        processing.concat_alignments(options)

    # Filter fasta files by length
    if options.filter_fasta:
        target, lo, hi = options.filter_fasta
        processing.filter_length_fasta(target, lo, hi)

    # Add taxonomy information to alignment(s)
    if options.add_taxonomy:
        processing.taxonomy_comprehension(options)

    # Add genomic context: schreibt Textdatei zu bereitgestellten Sequenzen
    if options.add_genomic_context:
        processing.add_genomic_context(
            options.database_directory, options.add_genomic_context
        )

    # iTol domain (Gene-Cluster) dataset
    if options.create_gene_cluster_dataset:
        directory = os.path.dirname(options.create_gene_cluster_dataset)
        datasets.iTol_domain_dataset(
            directory,
            options.database_directory,
            options.create_gene_cluster_dataset,
            options.dataset_divide_sign,
        )

    # iTol range dataset (per protein type)
    if options.create_type_range_dataset:
        if not options.database_directory:
            log.warning("Missing database to assign taxonomy; use -db argument")
            return
        directory = os.path.dirname(options.create_type_range_dataset)
        datasets.iTol_range_dataset(
            directory,
            options.database_directory,
            options.create_type_range_dataset,
            options.dataset_divide_sign,
        )

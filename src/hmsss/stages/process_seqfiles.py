from __future__ import annotations

import os

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import processing as processing
from hmsss.io import datasets as datasets



log = get_logger(__name__)

"""
Stage: Processing operators.

Applies post-processing operators to sequence and alignment files:
- Merge FASTA
- Concatenate alignments
- Filter FASTA by length
- Add taxonomy to alignments
- Add genomic context
- Create iTol datasets (gene cluster or type range)
"""

def process_operator(config: Config) -> None:
    """Apply processing operators according to CLI options.

    The following operations are conditionally executed:
      - `merge_fasta`: Merge FAA files into one without duplicates.
      - `concat_alignment`: Concatenate alignment files into one.
      - `filter_fasta`: Filter sequences by minimum/maximum length.
      - `add_taxonomy`: Add taxonomy information to alignments.
      - `add_genomic_context`: Write genomic context text file for given sequences.
      - `create_gene_cluster_dataset`: Build iTol domain dataset.
      - `create_type_range_dataset`: Build iTol range dataset.

    Args:
        options: Configuration object with operator flags and parameters.

    Side Effects:
        Creates/updates files and directories as required by the operators.
    """
    print_header("Processing operators", logger=log)

    # Merge fasta files
    if config.merge_fasta:
        processing.merge_fasta(config)

    # Concat alignment files
    if config.concat_alignment:
        processing.concat_alignments(config)

    # Filter fasta files by length
    if config.filter_fasta:
        target, lo, hi = config.filter_fasta
        processing.filter_length_fasta(target, lo, hi)

    # Add taxonomy information to alignment(s)
    if config.add_taxonomy:
        processing.taxonomy_comprehension(config)

    # Add genomic context: schreibt Textdatei zu bereitgestellten Sequenzen
    if config.add_genomic_context:
        processing.add_genomic_context(
            config.database_directory, config.add_genomic_context
        )

    # iTol domain (Gene-Cluster) dataset
    if config.create_gene_cluster_dataset:
        directory = os.path.dirname(config.create_gene_cluster_dataset)
        datasets.iTol_domain_dataset(
            directory,
            config.database_directory,
            config.create_gene_cluster_dataset,
            config.dataset_divide_sign,
        )

    # iTol range dataset (per protein type)
    if config.create_type_range_dataset:
        if not config.database_directory:
            log.warning("Missing database to assign taxonomy; use -db argument")
            return
        directory = os.path.dirname(config.create_type_range_dataset)
        datasets.iTol_range_dataset(
            directory,
            config.database_directory,
            config.create_type_range_dataset,
            config.dataset_divide_sign,
        )

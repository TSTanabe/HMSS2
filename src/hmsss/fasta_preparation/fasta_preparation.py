from __future__ import annotations

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.fasta_preparation import translation, pyrodigal_translation
from hmsss.core.queue import queue_fna_inputs, queue_faa_without_gff

logger = get_logger(__name__)

"""
Stage: FASTA preparation.

Translates `.fna` files to `.faa`/`.gff` pairs using Prodigal
and transcribes additional FAA files if no GFF is present.
"""


def fasta_preparation(config: Config) -> None:
    """Prepare FASTA files for the pipeline.

    - Collect `.fna` files and translate them to FAA/GFF using Prodigal.
    - Collect FAA files lacking GFF and transcribe them in parallel.

    Args:
        config: Configuration with `.cores` and input directories.

    Side Effects:
        Creates `.faa` and `.gff` files as required.
    """
    print_header("FASTA preparation (gene calling / translation)", logger=logger)
    fna_files: dict[str, str] = queue_fna_inputs(config)
    pyrodigal_translation.parallel_pyrodigal_translation(fna_files, config.cores)
    # translation.parallel_translation(fna_files, config.cores) Needs a lot more context switches

    faa_files_without_gff: dict[str, str] = queue_faa_without_gff(config)
    translation.parallel_transcription(faa_files_without_gff, config.cores)

    return

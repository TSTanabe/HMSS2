from __future__ import annotations

from hmsss.core.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms.translation import parallel_translation, parallel_transcription
from hmsss.io.queue import queue_fna_inputs, queue_faa_without_gff

log = get_logger(__name__)


def fasta_preparation(config: Config) -> None:
    """
    Translate and transcripe the .fna files to faa/gff pairs
    """
    print_header("FASTA preparation (gene calling / translation)", logger=log)

    fna_files: dict[str, str] = queue_fna_inputs(config)
    parallel_translation(fna_files, config.cores)

    faa_files_without_gff: dict[str, str] = queue_faa_without_gff(config)
    parallel_transcription(faa_files_without_gff, config.cores)

    return

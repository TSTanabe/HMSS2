from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms.translation import parallel_translation, parallel_transcription
from hmsss.io.queue import queue_fna_inputs, queue_faa_without_gff

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def fasta_preparation(options: Hmsss) -> None:
    """
    Translate and transcripe the .fna files to faa/gff pairs
    """
    print_header("FASTA preparation (gene calling / translation)", logger=log)

    fna_files: dict[str, str] = queue_fna_inputs(options)
    parallel_translation(fna_files, options.cores)

    faa_files_without_gff: dict[str, str] = queue_faa_without_gff(options)
    parallel_transcription(faa_files_without_gff, options.cores)

    return

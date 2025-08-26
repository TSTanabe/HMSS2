from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms.translation import parallel_translation, parallel_transcription

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def fasta_preparation(options: Hmsss) -> None:
    """
    Translate and transcripe the .fna files to faa/gff pairs
    """
    print_header("FASTA preparation (gene calling / translation)", logger=log)

    in_dir = Path(options.fasta_file_directory or "")
    if not in_dir.is_dir():
        raise SystemExit(f"Input FASTA directory does not exist: {in_dir}")
    log.info("Input genomes directory: %s", in_dir)

    parallel_translation(options.fasta_file_directory, options.cores)
    parallel_transcription(options.fasta_file_directory, options.cores)

    return

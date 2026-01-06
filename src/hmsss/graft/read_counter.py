import gzip
import os
import multiprocessing as mp
from typing import Any

from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def safe_read_count(path: str | None) -> int:
    """Return read count or 0 if path is None, missing, or empty."""
    if path is None:
        return 0
    if not isinstance(path, str):
        logger.error("Read were not counted, missformated path: {path}")
        return 0
    if not os.path.isfile(path):
        logger.error("Read were not counted, path is not a file: {path}")
        return 0
    if os.path.getsize(path) == 0:
        logger.error("Read were not counted, file is empty: {path}")
        return 0
    try:
        return count_reads_fasta_or_fastq(path)
    except Exception:
        return 0


def count_reads_fasta_or_fastq(path: str) -> int:
    """
    Count reads/sequences in a FASTA/FASTQ file.
    Supports gzipped (.gz) and uncompressed files.

    Rules:
      - FASTQ: count records as total_lines // 4
      - FASTA (and other fasta-like): count headers starting with '>'

    Format detection:
      - If first non-empty line starts with '@' -> FASTQ
      - If first non-empty line starts with '>' -> FASTA
      - Else: treat as FASTA-like and count '>' headers (raises if none found)
    """
    if not isinstance(path, str):
        raise TypeError("path must be a string")
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as f:
        # read first non-empty line to detect format
        first = ""
        while True:
            first = f.readline()
            if first == "":  # EOF
                return 0
            if first.strip():
                break

        f.seek(0)

        # FASTQ: count lines // 4 (valid FASTQ has 4 lines per record)
        if first.startswith("@"):
            n_lines = sum(1 for _ in f)
            return n_lines // 4

        # FASTA: count headers
        if first.startswith(">"):
            return sum(1 for line in f if line.startswith(">"))

        # Fallback: count '>' anyway; error if none found (unknown format)
        n_headers = sum(1 for line in f if line.startswith(">"))
        if n_headers == 0:
            raise ValueError(
                f"Unknown format (neither FASTQ '@' nor FASTA '>'): {path}"
            )
        return n_headers


def guess_extension(path: str) -> str:
    """
    Return the matched extension (including multi-part like '.fastq.gz').
    Raises ValueError if unknown.
    """
    # sort longest first so '.fastq.gz' matches before '.gz'
    exts = [
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".faa.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".fna",
        ".faa",
        ".gz",
    ]
    for ext in exts:
        if path.endswith(ext):
            return ext
    raise ValueError(f"Unable to guess file format of sequence file: {path}")


def read_basename(read_file: str) -> str:
    """
    Return filename without recognized sequencing extension.
    Example: '/x/y/sample_1.fq.gz' -> 'sample_1'
    """
    base = os.path.basename(read_file)
    ext = guess_extension(read_file)
    return base[: -len(ext)]


def metagenome_counts_from_task(task) -> dict[str, str | int | Any]:
    """
    Compute forward/reverse read counts for a task and build the dict required by
    insert_database_metagenomes().

    Returns
    -------
    (forward_reads, reverse_reads, metagenome_dict)

    metagenome_dict format:
      {metagenomeID: (genomeID, forward_reads, reverse_reads)}
    """
    forward_reads = safe_read_count(task.forward)
    reverse_reads = safe_read_count(task.reverse)

    task.metagenome_id = read_basename(task.forward)
    metagenomeID = task.metagenome_id  # set this in task creation, or use files["base"]
    genomeID = task.genome_id

    return {
        "metagenomeID": metagenomeID,
        "genomeID": genomeID,
        "forward_reads": forward_reads,
        "reverse_reads": reverse_reads,
    }


def collect_metagenome_counts_parallel(
        tasks: list[Any],
        *,
        processes: int | None = None,
        chunksize: int = 1,
) -> tuple[dict[str, tuple[str, int, int]], set[str]]:
    """
    Run metagenome_counts_from_task(task) in parallel.

    chunksize defines how many task a single worker gets
    processes defines how many workers exist

    Returns
    -------
    meta_dict : dict
        {metagenomeID: (genomeID, forward_reads, reverse_reads)}
    genome_ids : set
        {genomeID, ...}

    """
    meta_dict: dict[str, tuple[str, int, int]] = {}
    genome_ids: set[str] = set()

    ctx = mp.get_context("spawn")  # safer on many HPC setups
    with ctx.Pool(processes=processes) as pool:
        for res in pool.imap_unordered(
                metagenome_counts_from_task, tasks, chunksize=chunksize
        ):
            metagenomeID = res["metagenomeID"]
            genomeID = res["genomeID"]
            fwd = int(res["forward_reads"])
            rev = int(res["reverse_reads"])

            meta_dict[metagenomeID] = (genomeID, fwd, rev)
            genome_ids.add(genomeID)

    return meta_dict, genome_ids

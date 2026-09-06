import gzip
import os
import multiprocessing as mp
import subprocess
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
        logger.debug("Counting read for path: {path}")
        return count_fast_external(path)
    except Exception:
        return 0


def count_fast_external(path: str) -> int:
    if path.endswith(".fastq.gz") or path.endswith(".fq.gz"):
        cmd = f"zcat {path} | wc -l"
        lines = int(subprocess.check_output(cmd, shell=True))
        # print("Counting reads with zcat")
        return lines // 4
    elif path.endswith(".fasta.gz") or path.endswith(".fa.gz"):
        cmd = f"zgrep -c '^>' {path}"
        # print("Counting reads with zgrep")
        return int(subprocess.check_output(cmd, shell=True))
    elif path.endswith(".fasta") or path.endswith(".fa"):
        cmd = f"grep -c '^>' {path}"
        # print("Counting reads with grep")
        return int(subprocess.check_output(cmd, shell=True, text=True).strip())
    elif path.endswith(".fastq") or path.endswith(".fq"):
        cmd = f"wc -l < {path}"
        # print("Counting reads with wc")
        lines = int(subprocess.check_output(cmd, shell=True, text=True).strip())
        return lines // 4

    else:
        print("Not counting because error")
        raise ValueError(path)


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

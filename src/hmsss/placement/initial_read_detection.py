# src/hmsss/search/initial_read_detection.py

from __future__ import annotations

import logging
import multiprocessing as mp
import shlex
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pyhmmer

from hmsss.cli.config import INTERMEDIATE_HITS_SUFFIX
from hmsss.core import queue
from hmsss.db import database
from hmsss.graft import read_counter

logger = logging.getLogger(__name__)


def _iter_fasta_stream(handle) -> Iterator[Tuple[str, str]]:
    """
    Read FASTA records directly from a binary stream.

    Designed for OrfM stdout.

    Yields:
        (orf_id, amino_acid_sequence)

    Only the current FASTA record is buffered.
    """

    sequence_id = None
    sequence_parts = []

    for raw_line in handle:
        line = raw_line.decode("ascii").strip()

        if not line:
            continue

        if line.startswith(">"):
            if sequence_id is not None:
                yield sequence_id, "".join(sequence_parts)

            sequence_id = line[1:].split()[0]
            sequence_parts = []

            if not sequence_id:
                raise ValueError("Empty FASTA sequence identifier")

        else:
            if sequence_id is None:
                raise ValueError("FASTA sequence found before first header")

            sequence_parts.append(line)

    if sequence_id is not None:
        yield sequence_id, "".join(sequence_parts)


def _write_intermediate_orf_faa(
    detections: list[tuple],
    metagenome_id: str,
    output_dir: str | Path,
) -> int:
    """
    Write detected ORFs into protein-type-specific intermediate FASTA files.

    FASTA identifier:
        proteinID___dom_START_END

    where:
        proteinID = metagenomeID-orfmID

    Example:
        >SRR123456-read42___dom_15_128

    The proteinID itself remains unchanged in the database.
    dom_start/dom_end identify the corresponding Domains row.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(list)

    for (
        orfm_id,
        protein_type,
        sequence,
        _score,
        dom_start,
        dom_end,
    ) in detections:
        protein_id = f"{metagenome_id}-{orfm_id}"

        header = f"{protein_id}___dom_{int(dom_start)}_{int(dom_end)}"

        grouped[protein_type].append((header, sequence))

    written = 0

    for protein_type, sequences in grouped.items():
        output_file = output_dir / f"{protein_type}{INTERMEDIATE_HITS_SUFFIX}"

        with output_file.open(
            "a",
            encoding="utf-8",
        ) as handle:
            for header, sequence in sequences:
                handle.write(f">{header}\n{sequence}\n")
                written += 1

    logger.debug(
        "Wrote %d intermediate ORFs for %d protein types",
        written,
        len(grouped),
    )

    return written


_G_HMM_LIBRARY = None
_G_HMMSCAN_CPUS = 1


def _allocate_hmmscan_cpus(
    total_cores: int,
    target_cpus_per_worker: int = 6,
) -> list[int]:
    """
    Reserve one CPU for OrfM/main and distribute the remaining CPUs
    over hmmscan workers, targeting about 6 CPUs per worker.
    """
    total_cores = max(1, int(total_cores))

    if total_cores <= 2:
        return [1]

    hmmscan_cores = total_cores - 1

    worker_count = max(
        1,
        round(hmmscan_cores / target_cpus_per_worker),
    )
    worker_count = min(worker_count, hmmscan_cores)

    base = hmmscan_cores // worker_count
    remainder = hmmscan_cores % worker_count

    return [base + (1 if i < remainder else 0) for i in range(worker_count)]


def _init_hmmscan_worker(
    library_path: str,
    cpu_queue,
) -> None:
    """Load the HMM library once per worker."""
    global _G_HMM_LIBRARY, _G_HMMSCAN_CPUS

    _G_HMMSCAN_CPUS = max(1, int(cpu_queue.get()))

    with pyhmmer.plan7.HMMFile(library_path) as hmm_file:
        _G_HMM_LIBRARY = list(hmm_file)


def _hmmscan_orf_batch_worker(
    batch_records: list[tuple[str, str]],
    evalue_cutoff: float,
) -> list[tuple]:
    """
    Scan one ORF batch in a worker process.

    Workers only perform hmmscan and return detections.
    Database and intermediate FASTA writes stay in the main process.
    """
    if not batch_records:
        return []

    if _G_HMM_LIBRARY is None:
        raise RuntimeError("HMM worker library was not initialized")

    alphabet = pyhmmer.easel.Alphabet.amino()
    batch = pyhmmer.easel.DigitalSequenceBlock(alphabet)

    for orfm_id, amino_acids in batch_records:
        batch.append(
            pyhmmer.easel.TextSequence(
                name=orfm_id.encode("utf-8"),
                sequence=amino_acids,
            ).digitize(alphabet)
        )

    results = pyhmmer.hmmer.hmmscan(
        batch,
        _G_HMM_LIBRARY,
        cpus=_G_HMMSCAN_CPUS,
        E=evalue_cutoff,
        domE=evalue_cutoff,
    )

    detections = []

    for sequence, hits in zip(batch, results):
        if not hits:
            continue

        best_hit = hits[0]
        best_domain = best_hit.best_domain

        detections.append(
            (
                sequence.name.decode(),
                best_hit.name.decode().split("_")[-1],
                sequence.textize().sequence,
                float(best_hit.score),
                int(best_domain.env_from),
                int(best_domain.env_to),
            )
        )

    return detections


def _store_hmmscan_batch(
    detections: list[tuple],
    database_path: str,
    intermediate_dir: str | Path,
    metagenome_id: str,
) -> Tuple[int, int]:
    """Write one completed worker batch in the main process."""
    if not detections:
        return 0, 0

    inserted = database.insert_metagenome_orf_batch(
        database=database_path,
        metagenome_id=metagenome_id,
        detections=detections,
    )

    _write_intermediate_orf_faa(
        detections=detections,
        metagenome_id=metagenome_id,
        output_dir=intermediate_dir,
    )

    return len(detections), inserted


def _build_orfm_command(
    fastq_file: str | Path,
    min_orf_length: int = 96,
    translation_table: Optional[int] = 11,
    orfm_extra_args: Optional[List[str]] = None,
) -> list[str]:
    """
    Build the OrfM command and log the complete command.

    Returns
    -------
    list[str]
        Command suitable for subprocess.Popen().
    """

    orfm_executable = shutil.which("orfm")

    if orfm_executable is None:
        raise RuntimeError("OrfM executable not found in PATH")

    command = [
        orfm_executable,
        "-m",
        str(min_orf_length),
    ]

    if translation_table is not None:
        command.extend(["-c", str(translation_table)])

    if orfm_extra_args:
        command.extend(str(arg) for arg in orfm_extra_args)

    command.append(str(fastq_file))

    logger.debug("OrfM command: %s", shlex.join(command))

    return command


def _count_orfm_orfs(
    command: list[str],
    stderr_path: str | Path,
    metagenome_id: str,
) -> int:
    """
    Run OrfM once and count the number of ORFs it emits.

    Only FASTA header lines are counted; amino-acid sequences are discarded.
    This gives the exact number of ORFs that the second OrfM pass will emit
    with the same command-line parameters.
    """

    stderr_path = Path(stderr_path)
    total_orfs = 0

    logger.debug(
        "[%s] First OrfM pass: counting total ORFs",
        metagenome_id,
    )

    with stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            bufsize=1024 * 1024,
        )

        try:
            if process.stdout is None:
                raise RuntimeError("Failed to capture OrfM stdout")

            for raw_line in process.stdout:
                if raw_line.startswith(b">"):
                    total_orfs += 1

                    if total_orfs % 1_000_000 == 0:
                        logger.debug(
                            "[%s] First OrfM pass: %d ORFs counted",
                            metagenome_id,
                            total_orfs,
                        )

            process.stdout.close()
            returncode = process.wait()

            if returncode != 0:
                raise RuntimeError(
                    f"First OrfM pass failed for {metagenome_id} "
                    f"with exit code {returncode}. "
                    f"See {stderr_path}"
                )

        except BaseException:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise

        finally:
            if process.stdout is not None:
                process.stdout.close()

    logger.info(
        "[%s] OrfM prediction found %d total ORFs",
        metagenome_id,
        total_orfs,
    )

    return total_orfs


# ============================================================
# 5. PROCESS ONE FASTQ WITH ORFM + PYHMMER
# ============================================================


def search_fastq_with_orfm(
    fastq_file: str | Path,
    metagenome_id: str,
    database_path: str,
    intermediate_dir: str | Path,
    executor: ProcessPoolExecutor,
    hmmscan_workers: int,
    output_dir: str | Path,
    batch_size: int = 10000,
    min_orf_length: int = 96,
    translation_table: int = 11,
    evalue: float = 1e-5,
    orfm_extra_args: Optional[List[str]] = None,
) -> dict:
    """
    Translate one FASTQ with OrfM and stream ORF batches to parallel hmmscan workers.

    OrfM runs continuously in the main process. Worker processes perform only
    hmmscan. Completed detections are written to SQLite and intermediate FASTA
    files by the main process.
    """

    fastq_file = Path(fastq_file)
    output_dir = Path(output_dir)

    if not fastq_file.is_file():
        raise FileNotFoundError(fastq_file)

    # --------------------------------------------------------
    # Count input reads before OrfM
    # --------------------------------------------------------

    total_reads = read_counter.safe_read_count(str(fastq_file))

    logger.info(
        "[%s] Input reads: %d",
        metagenome_id,
        total_reads,
    )

    database.insert_database_genome_ids(database_path, {metagenome_id})
    database.insert_database_metagenomes(
        database=database_path,
        metagenome_dict={
            metagenome_id: (
                metagenome_id,
                total_reads,
                None,
            )
        },
    )

    # --------------------------------------------------------
    # Prepare OrfM command
    # --------------------------------------------------------

    command = _build_orfm_command(
        fastq_file=fastq_file,
        min_orf_length=min_orf_length,
        translation_table=translation_table,
        orfm_extra_args=orfm_extra_args,
    )

    total_orfs = 0
    total_orfs_expected = 0
    detected_orfs = 0
    stored_proteins = 0
    batch_number = 0
    failed_batches = 0
    processed_orfs = 0

    # --------------------------------------------------------
    # First OrfM pass: determine exact number of ORFs
    # --------------------------------------------------------

    count_stderr_path = output_dir / "orfm_count.stderr.log"
    stderr_path = output_dir / "orfm.stderr.log"

    total_orfs_expected = _count_orfm_orfs(
        command=command,
        stderr_path=count_stderr_path,
        metagenome_id=metagenome_id,
    )

    logger.debug(
        "[%s] Second OrfM pass: streaming %d expected ORFs to PyHMMER",
        metagenome_id,
        total_orfs_expected,
    )

    # Keep only a small number of outstanding batches in RAM.
    max_in_flight = max(2, hmmscan_workers * 2)
    in_flight = set()
    future_info = {}
    batch_records: list[tuple[str, str]] = []

    def collect_finished(done_futures) -> None:
        nonlocal detected_orfs, stored_proteins, failed_batches, processed_orfs

        for future in done_futures:
            batch_id, batch_orfs = future_info.pop(future)

            try:
                detections = future.result()
            except Exception:
                failed_batches += 1
                logger.exception(
                    "[%s] hmmscan batch %d failed; skipping batch",
                    metagenome_id,
                    batch_id,
                )
                continue

            detected, stored = _store_hmmscan_batch(
                detections=detections,
                database_path=database_path,
                intermediate_dir=intermediate_dir,
                metagenome_id=metagenome_id,
            )

            processed_orfs += batch_orfs
            detected_orfs += detected
            stored_proteins += stored

            percent = (
                100.0 * processed_orfs / total_orfs_expected
                if total_orfs_expected
                else 0.0
            )

            logger.info(
                "[%s] %d/%d ORFs (%.2f%%) | detected batch: %d | detected total: %d",
                metagenome_id,
                processed_orfs,
                total_orfs_expected,
                percent,
                detected,
                detected_orfs,
            )

    def submit_batch(records: list[tuple[str, str]]) -> None:
        nonlocal batch_number, in_flight

        if not records:
            return

        batch_number += 1
        batch_id = batch_number

        future = executor.submit(
            _hmmscan_orf_batch_worker,
            records,
            evalue,
        )

        in_flight.add(future)
        future_info[future] = (
            batch_id,
            len(records),
        )

        # Backpressure: do not let ORF batches accumulate indefinitely.
        if len(in_flight) >= max_in_flight:
            done, pending = wait(
                in_flight,
                return_when=FIRST_COMPLETED,
            )
            in_flight = pending
            collect_finished(done)

    # --------------------------------------------------------
    # Second OrfM pass: stream ORFs into parallel hmmscan
    # --------------------------------------------------------

    with stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            bufsize=1024 * 1024,
        )

        try:
            if process.stdout is None:
                raise RuntimeError("Failed to capture OrfM stdout")

            for orfm_id, amino_acids in _iter_fasta_stream(process.stdout):
                if not amino_acids:
                    continue

                batch_records.append((orfm_id, amino_acids))
                total_orfs += 1

                if len(batch_records) >= batch_size:
                    submit_batch(batch_records)
                    batch_records = []

            # Process final incomplete batch.
            if batch_records:
                submit_batch(batch_records)
                batch_records = []

            process.stdout.close()
            returncode = process.wait()

            if returncode != 0:
                raise RuntimeError(
                    f"OrfM failed for {metagenome_id} "
                    f"with exit code {returncode}. "
                    f"See {stderr_path}"
                )

            # Collect all outstanding hmmscan batches.
            while in_flight:
                done, pending = wait(
                    in_flight,
                    return_when=FIRST_COMPLETED,
                )
                in_flight = pending
                collect_finished(done)

        except BaseException:
            if process.poll() is None:
                process.terminate()

                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()

            for future in in_flight:
                future.cancel()

            raise

        finally:
            if process.stdout is not None:
                process.stdout.close()

    logger.debug(
        "[%s] Completed: %d ORFs translated, %d detected, %d inserted/updated",
        metagenome_id,
        total_orfs,
        detected_orfs,
        stored_proteins,
    )

    if failed_batches:
        logger.warning(
            "[%s] %d hmmscan batches failed",
            metagenome_id,
            failed_batches,
        )

    return {
        "metagenome_id": metagenome_id,
        "fastq_file": str(fastq_file),
        "total_reads": total_reads,
        "total_orfs_expected": total_orfs_expected,
        "total_orfs": total_orfs,
        "processed_orfs": processed_orfs,
        "detected_orfs": detected_orfs,
        "stored_proteins": stored_proteins,
        "batches": batch_number,
        "failed_batches": failed_batches,
    }


# ============================================================
# 6. STAGE 30 PIPELINE CONTROLLER
# ============================================================


def initial_read_detection(
    config,
    *,
    min_orf_length: int = 96,
    translation_table: int = 11,
    batch_size: int = 10000,
    evalue: float = 1e-5,
    orfm_extra_args: Optional[List[str]] = None,
) -> dict:
    """
    Stage 30: initial metagenomic read detection.

    Workflow:
        1. Collect compatible search HMMs.
        2. Create one combined search library.
        3. Load optimized HMM profiles once.
        4. Collect all FASTQ files from config.
        5. Process each FASTQ with OrfM.
        6. Stream translated ORFs into PyHMMER hmmscan.
        7. Save best HMM detection for every ORF.
        8. Write directly into the HMSS3 database.

    Parameters
    ----------
    config:
        HMSS3 configuration object.

    min_orf_length:
        Minimum ORF length passed to OrfM.

    translation_table:
        NCBI translation table passed to OrfM.

    batch_size:
        Number of translated ORFs per PyHMMER batch.

    evalue:
        Maximum HMMER E-value for accepted hits.

    orfm_extra_args:
        Additional OrfM command-line arguments.

    Returns
    -------
    dict:
        {
            metagenome_id: {
                "total_orfs": ...,
                "detected_orfs": ...,
                "stored_proteins": ...,
                ...
            }
        }
    """

    results = {}
    queue.queue_read_mapping_fastq_inputs(config)
    database.create_database(config.database_directory)

    cpu_allocations = _allocate_hmmscan_cpus(
        total_cores=int(config.cores),
        target_cpus_per_worker=6,
    )
    worker_count = len(cpu_allocations)

    logger.info(
        "Read detection: 1 OrfM/main process + %d hmmscan workers (%s CPUs)",
        worker_count,
        ", ".join(str(cpus) for cpus in cpu_allocations),
    )

    ctx = mp.get_context()
    cpu_queue = ctx.Queue()

    for worker_cpus in cpu_allocations:
        cpu_queue.put(worker_cpus)

    try:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=ctx,
            initializer=_init_hmmscan_worker,
            initargs=(str(config.library), cpu_queue),
        ) as executor:
            for metagenome_id, fastq_file in config.fastq_files.items():
                try:
                    result = search_fastq_with_orfm(
                        fastq_file=fastq_file,
                        metagenome_id=metagenome_id,
                        database_path=config.database_directory,
                        intermediate_dir=config.cross_check_directory,
                        executor=executor,
                        hmmscan_workers=worker_count,
                        output_dir=config.result_files_directory,
                        batch_size=batch_size,
                        min_orf_length=min_orf_length,
                        translation_table=translation_table,
                        evalue=evalue,
                        orfm_extra_args=orfm_extra_args,
                    )

                    results[metagenome_id] = result

                except Exception:
                    logger.exception(
                        "[%s] Read detection failed; continuing with next metagenome",
                        metagenome_id,
                    )

    finally:
        cpu_queue.close()
        cpu_queue.join_thread()

    return results

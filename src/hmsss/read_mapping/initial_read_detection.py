# src/hmsss/search/initial_read_detection.py

from __future__ import annotations

import gzip
import logging
import os
import shlex
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pyhmmer

from hmsss.core import queue
from hmsss.db import database
from hmsss.graft import read_counter

logger = logging.getLogger(__name__)


def _iter_fasta_stream(handle) -> Iterator[Tuple[str, str]]:
    """
    Read FASTA records from a binary file handle.

    Yields
    ------
    tuple[str, str]
        (sequence_id, amino_acid_sequence)

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


def _prepare_orfm_input(
    fastq_file: Path,
    temp_dir: Path,
    metagenome_id: str,
    progress_every: int = 1_000_000,
) -> tuple[Path, int]:
    """
    Prepare an OrfM-readable FASTQ input and determine its read count.

    Plain FASTQ files are used directly. Gzipped FASTQ files are decompressed
    once into the run-local temporary directory while reads are counted. This
    avoids feeding OrfM through stdin and therefore removes the Python -> OrfM
    pipe from the processing stage.

    Returns
    -------
    tuple[Path, int]
        (path passed to OrfM, total_reads)
    """

    name = fastq_file.name.lower()

    if not name.endswith((".fastq.gz", ".fq.gz")):
        total_reads = read_counter.safe_read_count(str(fastq_file))
        return fastq_file, total_reads

    safe_id = metagenome_id.replace(os.sep, "_")
    uncompressed_fastq = temp_dir / f"{safe_id}.fastq"
    total_reads = 0

    logger.info(
        "[%s] Decompressing FASTQ to temporary storage: %s",
        metagenome_id,
        uncompressed_fastq,
    )

    with gzip.open(fastq_file, "rb") as source, uncompressed_fastq.open("wb") as target:
        while True:
            header = source.readline()

            if not header:
                break

            sequence = source.readline()
            plus = source.readline()
            quality = source.readline()

            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record in {fastq_file}")

            target.write(header)
            target.write(sequence)
            target.write(plus)
            target.write(quality)

            total_reads += 1

            if total_reads % progress_every == 0:
                logger.info(
                    "[%s] Decompressed %d reads",
                    metagenome_id,
                    total_reads,
                )

    logger.info(
        "[%s] Temporary FASTQ ready: %d reads | %.2f GiB",
        metagenome_id,
        total_reads,
        uncompressed_fastq.stat().st_size / (1024**3),
    )

    return uncompressed_fastq, total_reads


def _write_intermediate_orf_faa(
    detections: list[tuple],
    metagenome_id: str,
    output_dir: str | Path,
) -> int:
    """
    Write detected ORFs into protein-type-specific FASTA files.

    Each detection contains:
        (
            orfm_id,
            protein_type,
            sequence,
            bitscore,
            dom_start,
            dom_end
        )

    FASTA identifier:
        metagenomeID-orfmID

    Files are appended across batches and metagenomes.

    Returns
    -------
    int
        Number of sequences written.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(list)

    for orfm_id, protein_type, sequence, *_ in detections:
        protein_id = f"{metagenome_id}-{orfm_id}"
        grouped[protein_type].append((protein_id, sequence))

    written = 0

    for protein_type, sequences in grouped.items():
        output_file = output_dir / f"{protein_type}.intermediate_hits.faa"

        with output_file.open("a", encoding="utf-8") as handle:
            for protein_id, sequence in sequences:
                handle.write(f">{protein_id}\n{sequence}\n")
                written += 1

    logger.debug(
        "Wrote %d intermediate ORFs for %d protein types",
        written,
        len(grouped),
    )

    return written


def _hmmscan_orf_batch(
    batch,
    library,
    database_path: str,
    intermediate_dir: str | Path,
    metagenome_id: str,
    cores: int = 20,
    evalue_cutoff: float = 1e-5,
) -> Tuple[int, int]:
    """
    Scan one batch of translated ORFs using PyHMMER hmmscan.

    Every OrfM ORF is processed independently.

    Returns
    -------
    tuple[int, int]
        (detected_orfs, inserted_or_updated_proteins)
    """

    if not batch:
        return 0, 0

    detections = []

    results = pyhmmer.hmmer.hmmscan(
        batch,
        library,
        cpus=cores,
        E=evalue_cutoff,
        domE=evalue_cutoff,
    )

    for sequence, hits in zip(batch, results):
        if not hits:
            continue

        best_hit = hits[0]

        hmm_id = best_hit.name.decode()
        protein_type = hmm_id.split("_")[-1]

        orfm_id = sequence.name.decode()
        amino_acids = sequence.textize().sequence

        best_domain = best_hit.best_domain

        detections.append(
            (
                orfm_id,
                protein_type,
                amino_acids,
                float(best_hit.score),
                int(best_domain.env_from),
                int(best_domain.env_to),
            )
        )

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
    """Build the OrfM command for direct file input."""

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

    # Direct file input: no stdin pipe.
    command.append(str(fastq_file))

    logger.info("OrfM command: %s", shlex.join(command))

    return command


# ============================================================
# 5. PROCESS ONE FASTQ WITH ORFM + TEMPORARY FAA + PYHMMER
# ============================================================


def search_fastq_with_orfm(
    fastq_file: str | Path,
    metagenome_id: str,
    database_path: str,
    intermediate_dir: str | Path,
    profiles: str | Path,
    output_dir: str | Path,
    cores: int = 20,
    batch_size: int = 10000,
    min_orf_length: int = 96,
    translation_table: int = 11,
    evalue: float = 1e-5,
    orfm_extra_args: Optional[List[str]] = None,
) -> dict:
    """
    Translate one FASTQ with OrfM into a temporary FAA and then scan the
    resulting ORFs batchwise with PyHMMER.

    The two CPU-heavy phases are deliberately separated:

        FASTQ -> OrfM -> $TMPDIR/*.orfm.faa
        $TMPDIR/*.orfm.faa -> PyHMMER -> database/intermediate hits

    This removes the OrfM stdout -> Python pipe and prevents OrfM from being
    repeatedly blocked while PyHMMER scans a batch.
    """

    fastq_file = Path(fastq_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not fastq_file.is_file():
        raise FileNotFoundError(fastq_file)

    if not shutil.which("orfm"):
        raise RuntimeError("OrfM executable not found in PATH")

    # A singleton set is required here; set(metagenome_id) would split the
    # identifier into individual characters.
    database.insert_database_genome_ids(database_path, {metagenome_id})

    alphabet = pyhmmer.easel.Alphabet.amino()

    total_reads = 0
    total_orfs = 0
    detected_orfs = 0
    stored_proteins = 0
    batch_number = 0

    safe_id = metagenome_id.replace(os.sep, "_")

    tmp_root = Path(output_dir) / "tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)

    stderr_path = output_dir / f"{safe_id}.orfm.stderr.log"

    logger.info(
        "[%s] Temporary working directory root: %s",
        metagenome_id,
        tmp_root.resolve(),
    )

    with tempfile.TemporaryDirectory(
        prefix=f"{safe_id}_",
        dir=tmp_root,
    ) as tmp_dir:
        tmp_dir = Path(tmp_dir)

        orf_faa = tmp_dir / f"{safe_id}.orfm.faa"

        # weitere temporäre Dateien hier

        # ----------------------------------------------------
        # 1. Prepare input for OrfM
        # ----------------------------------------------------
        orfm_input, total_reads = _prepare_orfm_input(
            fastq_file=fastq_file,
            temp_dir=temp_dir,
            metagenome_id=metagenome_id,
        )

        logger.info(
            "[%s] Input reads: %d",
            metagenome_id,
            total_reads,
        )

        # ----------------------------------------------------
        # 2. Run OrfM and write ALL ORFs to local temporary FAA
        # ----------------------------------------------------
        orf_faa = temp_dir / f"{safe_id}.orfm.faa"

        command = _build_orfm_command(
            fastq_file=orfm_input,
            min_orf_length=min_orf_length,
            translation_table=translation_table,
            orfm_extra_args=orfm_extra_args,
        )

        logger.info(
            "[%s] Running OrfM -> %s",
            metagenome_id,
            orf_faa,
        )

        with (
            orf_faa.open("wb") as orf_handle,
            stderr_path.open("wb") as stderr_handle,
        ):
            completed = subprocess.run(
                command,
                stdout=orf_handle,
                stderr=stderr_handle,
                check=False,
            )

        if completed.returncode != 0:
            raise RuntimeError(
                f"OrfM failed for {metagenome_id} "
                f"with exit code {completed.returncode}. "
                f"See {stderr_path}"
            )

        if not orf_faa.is_file() or orf_faa.stat().st_size == 0:
            logger.warning(
                "[%s] OrfM produced no ORFs",
                metagenome_id,
            )
            return {
                "metagenome_id": metagenome_id,
                "fastq_file": str(fastq_file),
                "total_reads": total_reads,
                "total_orfs": 0,
                "detected_orfs": 0,
                "stored_proteins": 0,
                "batches": 0,
            }

        logger.info(
            "[%s] OrfM finished | temporary FAA: %.2f GiB",
            metagenome_id,
            orf_faa.stat().st_size / (1024**3),
        )

        # ----------------------------------------------------
        # 3. Read temporary FAA and scan ORFs with PyHMMER
        # ----------------------------------------------------
        batch = pyhmmer.easel.DigitalSequenceBlock(alphabet)

        with orf_faa.open("rb") as handle:
            for orfm_id, amino_acids in _iter_fasta_stream(handle):
                if not amino_acids:
                    continue

                digital_sequence = pyhmmer.easel.TextSequence(
                    name=orfm_id.encode("utf-8"),
                    sequence=amino_acids,
                ).digitize(alphabet)

                batch.append(digital_sequence)
                total_orfs += 1

                if len(batch) < batch_size:
                    continue

                batch_number += 1

                detected, stored = _hmmscan_orf_batch(
                    batch=batch,
                    library=profiles,
                    database_path=database_path,
                    intermediate_dir=intermediate_dir,
                    metagenome_id=metagenome_id,
                    cores=cores,
                    evalue_cutoff=evalue,
                )

                detected_orfs += detected
                stored_proteins += stored

                logger.info(
                    "[%s] PyHMMER batch %d | ORFs processed: %d | "
                    "detected: %d | stored/updated: %d",
                    metagenome_id,
                    batch_number,
                    total_orfs,
                    detected_orfs,
                    stored_proteins,
                )

                batch = pyhmmer.easel.DigitalSequenceBlock(alphabet)

        # ----------------------------------------------------
        # 4. Process final incomplete batch
        # ----------------------------------------------------
        if batch:
            batch_number += 1

            detected, stored = _hmmscan_orf_batch(
                batch=batch,
                library=profiles,
                database_path=database_path,
                intermediate_dir=intermediate_dir,
                metagenome_id=metagenome_id,
                cores=cores,
                evalue_cutoff=evalue,
            )

            detected_orfs += detected
            stored_proteins += stored

        del batch

        # TemporaryDirectory removes orf_faa and, for gzipped input, the
        # decompressed temporary FASTQ automatically after this block.

    logger.info(
        "[%s] Completed: %d reads | %d ORFs translated | %d detected | "
        "%d inserted/updated",
        metagenome_id,
        total_reads,
        total_orfs,
        detected_orfs,
        stored_proteins,
    )

    return {
        "metagenome_id": metagenome_id,
        "fastq_file": str(fastq_file),
        "total_reads": total_reads,
        "total_orfs": total_orfs,
        "detected_orfs": detected_orfs,
        "stored_proteins": stored_proteins,
        "batches": batch_number,
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
    Stage 30: initial metagenomic read detection using a temporary OrfM FAA.

    Workflow:
        1. Collect FASTQ files.
        2. Load the HMM library once.
        3. For each metagenome, prepare an OrfM-readable FASTQ in $TMPDIR if
           decompression is necessary.
        4. Run OrfM to completion and write all ORFs to $TMPDIR.
        5. Read the temporary ORF FAA batchwise with PyHMMER.
        6. Save best HMM detections to the database and intermediate FASTAs.
        7. Remove temporary files automatically.
    """

    results = {}
    queue.queue_read_mapping_fastq_inputs(config)
    database.create_database(config.database_directory)

    # Load HMM library once.
    with pyhmmer.plan7.HMMFile(config.library) as hmm_file:
        profiles = list(hmm_file)

    logger.info(
        "Loaded %d HMM profiles from %s",
        len(profiles),
        config.library,
    )

    for metagenome_id, fastq_file in config.fastq_files.items():
        result = search_fastq_with_orfm(
            fastq_file=fastq_file,
            metagenome_id=metagenome_id,
            database_path=config.database_directory,
            intermediate_dir=config.cross_check_directory,
            profiles=profiles,
            output_dir=config.result_files_directory,
            cores=int(config.cores),
            batch_size=batch_size,
            min_orf_length=min_orf_length,
            translation_table=translation_table,
            evalue=evalue,
            orfm_extra_args=orfm_extra_args,
        )

        results[metagenome_id] = result

    return results

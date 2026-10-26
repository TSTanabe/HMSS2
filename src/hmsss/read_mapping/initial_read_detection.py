# src/hmsss/search/initial_read_detection.py

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pyhmmer

from hmsss.cli import config
from hmsss.core import queue
from hmsss.db import database

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

    # Group sequences by protein type
    grouped = defaultdict(list)

    for orfm_id, protein_type, sequence, *_ in detections:
        protein_id = f"{metagenome_id}-{orfm_id}"
        grouped[protein_type].append((protein_id, sequence))

    written = 0

    # Write one FASTA file per protein type
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

    For each ORF:
        1. Identify significant HMM hits.
        2. Select the hit with the highest bit score.
        3. Extract its protein type.
        4. Preserve the original amino acid sequence.
        5. Write the detection into HMSS3.

    Returns:
        (detected_orfs, inserted_or_updated_proteins)
    """

    if not batch:
        return 0, 0

    detections = []

    # --------------------------------------------------------
    # Run hmmscan
    # --------------------------------------------------------

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
    written = _write_intermediate_orf_faa(
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

    logger.info("OrfM command: %s", shlex.join(command))

    return command


# ============================================================
# 5. PROCESS ONE FASTQ WITH ORFM + PYHMMER
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
    Translate a FASTQ file with OrfM and stream the amino acid
    sequences directly into PyHMMER hmmscan.

    No intermediate complete FAA file is generated.
    Only one ORF batch is held in memory.

    Detections are inserted directly into:
        Genomes
        Proteins
        ProteinSequences
        DomainTypes
        Domains

    Returns:
        Statistics for the processed metagenome.
    """

    fastq_file = Path(fastq_file)
    output_dir = Path(output_dir)

    if not fastq_file.is_file():
        raise FileNotFoundError(fastq_file)

    database.insert_database_genome_ids(database_path, set(metagenome_id))

    # --------------------------------------------------------
    # Find OrfM executable
    # --------------------------------------------------------

    orfm_executable = shutil.which("orfm")

    if orfm_executable is None:
        raise RuntimeError("OrfM executable not found in PATH")

    # --------------------------------------------------------
    # Prepare OrfM command
    # --------------------------------------------------------
    logger.debug(f"Detection of ORFs with orfM for {fastq_file}")
    command = _build_orfm_command(
        fastq_file=fastq_file,
        min_orf_length=min_orf_length,
        translation_table=translation_table,
        orfm_extra_args=orfm_extra_args,
    )

    alphabet = pyhmmer.easel.Alphabet.amino()

    total_orfs = 0
    detected_orfs = 0
    stored_proteins = 0
    batch_number = 0

    # --------------------------------------------------------
    # Start OrfM
    # --------------------------------------------------------

    stderr_path = output_dir / "orfm.stderr.log"

    # Write stderr directly to disk to avoid filling
    # an unconsumed stderr pipe.

    with stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            bufsize=1024 * 1024,
        )

        try:
            batch = pyhmmer.easel.DigitalSequenceBlock(alphabet)

            # ------------------------------------------------
            # Consume OrfM output
            # ------------------------------------------------

            for orfm_id, amino_acids in _iter_fasta_stream(process.stdout):
                if not amino_acids:
                    continue

                # Convert ORF to a PyHMMER digital sequence
                digital_sequence = pyhmmer.easel.TextSequence(
                    name=orfm_id.encode("utf-8"),
                    sequence=amino_acids,
                ).digitize(alphabet)

                batch.append(digital_sequence)

                total_orfs += 1

                # Wait until the batch is full
                if len(batch) < batch_size:
                    continue

                # --------------------------------------------
                # Scan complete batch
                # --------------------------------------------

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
                    "[%s] Batch %d | ORFs: %d | detected: %d | stored/updated: %d",
                    metagenome_id,
                    batch_number,
                    total_orfs,
                    detected_orfs,
                    stored_proteins,
                )

                # Release current batch
                batch = pyhmmer.easel.DigitalSequenceBlock(alphabet)

            # ------------------------------------------------
            # Process final incomplete batch
            # ------------------------------------------------

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

            # ------------------------------------------------
            # Wait for OrfM
            # ------------------------------------------------

            process.stdout.close()

            returncode = process.wait()

            if returncode != 0:
                raise RuntimeError(
                    f"OrfM failed for {metagenome_id} "
                    f"with exit code {returncode}. "
                    f"See {stderr_path}"
                )

        except BaseException:
            # Stop OrfM if PyHMMER or database writing fails.
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

    # --------------------------------------------------------
    # Return statistics
    # --------------------------------------------------------

    logger.info(
        "[%s] Completed: %d ORFs translated, %d detected, %d inserted/updated",
        metagenome_id,
        total_orfs,
        detected_orfs,
        stored_proteins,
    )

    return {
        "metagenome_id": metagenome_id,
        "fastq_file": str(fastq_file),
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
    for metagenome_id, fastq_file in config.fastq_files.items():
        result = search_fastq_with_orfm(
            fastq_file=fastq_file,
            metagenome_id=metagenome_id,
            database_path=config.database_directory,
            intermediate_dir=config.cross_check_directory,
            profiles=config.library,
            output_dir=config.results_directory,
            cores=int(config.cores),
            batch_size=batch_size,
            min_orf_length=min_orf_length,
            translation_table=translation_table,
            evalue=evalue,
            orfm_extra_args=orfm_extra_args,
        )

        results[metagenome_id] = result

    return results

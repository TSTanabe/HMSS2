#!/usr/bin/python
from __future__ import annotations

import os
import sqlite3
import subprocess
from multiprocessing import Manager, Pool, Process
from pathlib import Path
from queue import Queue
from typing import List, Set

from hmsss.core.logging import get_logger
from hmsss.utils.myUtil import find_executable

logger = get_logger(__name__)

"""
Cross-reference utilities for HMM hits.

Provides routines to:
- Compare intermediate hits with reference sequences
- Promote validated candidates with close hit with a reference sequence
"""

# Global shared variables
current_counter = None
counter_lock = None


##########################################################################################################################################################
#################### Cross check hits with reference sequences and add the hmmreport lines to the trusted cutoff hmmreport ###############################
##########################################################################################################################################################


def find_file_by_name_partial(root_dir: str, filename: str) -> list[str]:
    """
    Durchsucht ein Verzeichnis rekursiv nach Dateien, bei denen der Teil nach dem
    ersten Unterstrich (_) mit dem im filename angegebenen Teil übereinstimmt.

    Beispiel:
        filename = "abc_target.faa" -> sucht nach Dateien, bei denen der Teil
        nach dem ersten "_" gleich "target.faa" ist.

    Args:
        root_dir (str): Startverzeichnis für die Suche.
        filename (str): Referenzdateiname (z. B. 'abc_target.faa').

    Returns:
        list[str]: Liste der vollständigen Pfade zu gefundenen Dateien.
    """
    matches = []
    # Referenzsuffix = Teil nach erstem Unterstrich
    ref_suffix = filename.split("_", 1)[1] if "_" in filename else filename

    for dirpath, _, files in os.walk(root_dir):
        for f in files:
            file_suffix = f.split("_", 1)[1] if "_" in f else f
            if file_suffix == ref_suffix:
                matches.append(os.path.join(dirpath, f))

    return matches


def cross_check_candidates_with_reference_seqs(config, intermediate_dict) -> List[str]:
    """Cross-check intermediate FASTA files against reference sequences using DIAMOND."""

    logger.info("Cross-checking intermediate hits with reference sequences")

    refseq_dir = Path(config.paths.refseq)
    data_dir = Path(config.paths.data)
    cross_check_dir = Path(config.cross_check_directory)
    cross_check_dir.mkdir(parents=True, exist_ok=True)

    diamond = find_executable("diamond")
    failed = []

    for hmm_id, files in intermediate_dict.items():
        # Support single paths and lists of paths
        if isinstance(files, (str, Path)):
            files = [files]

        hmm_type = hmm_id.split("_")[-1]
        db_path = refseq_dir / f"{hmm_id}.dmnd"

        for inter_file in files:
            inter_file = Path(inter_file)

            if not inter_file.is_file():
                logger.warning("Intermediate FASTA missing: %s", inter_file)
                failed.append(hmm_id)
                continue

            output_file = (
                cross_check_dir
                / f"{inter_file.name.removesuffix('.intermediate_hits_faa')}.crosschecked.tsv"
            )

            if output_file.is_file():
                logger.debug("Cross-check already present: %s", output_file)
                continue

            try:
                # Create DIAMOND database if necessary
                if not db_path.is_file():
                    matches = find_file_by_name_partial(str(data_dir), f"{hmm_id}.faa")

                    if not matches:
                        logger.warning(
                            "No reference sequences available for %s", hmm_id
                        )
                        failed.append(hmm_id)
                        continue

                    faa_path = Path(matches[0])

                    subprocess.run(
                        [
                            diamond,
                            "makedb",
                            "--in",
                            str(faa_path),
                            "-d",
                            str(db_path),
                            "--quiet",
                        ],
                        check=True,
                    )

                # Run DIAMOND
                cmd = [
                    diamond,
                    "blastp",
                    "--query",
                    str(inter_file),
                    "--db",
                    str(db_path),
                    "--out",
                    str(output_file),
                    "--outfmt",
                    "6",
                    "--max-target-seqs",
                    "1",
                    "--id",
                    str(config.refseq_identity),
                    "--threads",
                    str(config.cores),
                    "--quiet",
                    f"--{config.diamond_speed_mode}",
                ]

                logger.debug("Running DIAMOND for %s", hmm_id)
                subprocess.run(cmd, check=True)

                if output_file.stat().st_size == 0:
                    output_file.unlink()
                    logger.debug("No reference matches for %s", hmm_id)

            except (OSError, subprocess.CalledProcessError) as e:
                logger.warning("DIAMOND cross-check failed for %s: %s", hmm_id, e)
                output_file.unlink(missing_ok=True)
                failed.append(hmm_id)

    return list(dict.fromkeys(failed))


#############################################################################
#################### Promote hits with refseq hit ###########################
#############################################################################
def set_valid_hit_true_for_protein_ids(
    database: str,
    protein_ids,
) -> int:
    """
    Set Proteins.valid_hit = 1 for given proteinIDs using a temporary table.
    Scales to arbitrarily many IDs without placeholder limits.
    """
    ids = sorted(set(protein_ids))
    if not ids:
        return 0

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # 1) Temp table anlegen (pro Connection sichtbar)
        cur.execute("""
            CREATE TEMP TABLE IF NOT EXISTS tmp_promote_ids (
                proteinID TEXT PRIMARY KEY
            )
        """)
        cur.execute("DELETE FROM tmp_promote_ids")

        # 2) IDs effizient einfügen
        cur.executemany(
            "INSERT OR IGNORE INTO tmp_promote_ids (proteinID) VALUES (?)",
            ((pid,) for pid in ids),
        )

        # 3) Update via EXISTS / JOIN
        cur.execute("""
            UPDATE Proteins
            SET valid_hit = 1
            WHERE proteinID IN (
                SELECT proteinID FROM tmp_promote_ids
            )
        """)

        con.commit()

    return len(ids)


def p_crosscheck_writer(queue: Queue, database_path: str, *, batch_size: int = 50000):
    """
    Asynchronous writer process.
    Receives sets of proteinIDs and updates SQLite (valid_hit = 1).
    """
    buffer = set()
    total_received = 0
    total_updated = 0

    while True:
        item = queue.get()
        if item is None:  # sentinel
            break

        # item: set[str]
        buffer.update(item)
        total_received += len(item)

        if len(buffer) >= batch_size:
            total_updated += set_valid_hit_true_for_protein_ids(database_path, buffer)
            buffer.clear()

    # final flush
    if buffer:
        total_updated += set_valid_hit_true_for_protein_ids(database_path, buffer)

    logger.info(
        "Crosscheck writer finished: received=%d ids, updated=%d unique ids",
        total_received,
        total_updated,
    )


def _crosscheck_worker(hmm_id: str, crosscheck_dir: str, queue) -> None:
    """
    Worker:
      - reads {hmm_id}.crosschecked.tsv
      - extracts promoted protein IDs (first column)
      - sends them to writer via queue

    NOTE: IDs must match DB proteinID format: genomeID-proteinID
    """
    crosscheck_path = os.path.join(crosscheck_dir, f"{hmm_id}.crosschecked.tsv")
    if not os.path.exists(crosscheck_path):
        logger.warning("Crosscheck file missing: %s", crosscheck_path)
        return

    promoted: Set[str] = set()
    with open(crosscheck_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            promoted.add(line.split("\t", 1)[0].strip())

    if promoted:
        # send as a set to reduce queue traffic
        queue.put(promoted)

    logger.info("%s: queued %d promoted IDs", hmm_id, len(promoted))


def promote_crosschecked_hits_to_db(
    *,
    crosscheck_dir: str,
    database_path: str,
    processes: int = 4,
    writer_batch_size: int = 50000,
):
    """
    Run crosscheck promotion with:
      - N worker processes
      - 1 async writer process (SQLite only there)
    """
    hmm_ids = [
        f.replace(".crosschecked.tsv", "")
        for f in os.listdir(crosscheck_dir)
        if f.endswith(".crosschecked.tsv")
    ]

    if not hmm_ids:
        logger.info("No crosschecked files found in %s", crosscheck_dir)
        return

    with Manager() as manager:
        queue = manager.Queue()

        # --- start writer ---
        writer = Process(
            target=p_crosscheck_writer,
            args=(queue, database_path),
            kwargs={"batch_size": writer_batch_size},
        )
        writer.start()

        # --- start workers ---
        with Pool(processes=processes) as pool:
            pool.starmap(
                _crosscheck_worker,
                [(hmm_id, crosscheck_dir, queue) for hmm_id in hmm_ids],
            )

        # --- stop writer ---
        queue.put(None)
        writer.join()

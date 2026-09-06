#!/usr/bin/python
from __future__ import annotations

import os
import sqlite3
import glob
import subprocess

from multiprocessing import Pool, Manager, Process
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


def cross_check_candidates_with_reference_seqs(config) -> List[str]:
    """Validate candidates via DIAMOND against reference sequences.

    For each `{hmm}.intermediate_hits_faa`, build or reuse a DIAMOND DB from
    reference FASTA (`{hmm}.faa` or a fallback by suffix), run `blastp`, and
    write `{hmm}.crosschecked.tsv` with the top hit (if any). Empty results
    are removed. Returns a list of HMM IDs where reference sequences were
    unavailable or the run failed.

    Args:
        config: Configuration with `paths.refseq`, `refseq_identity`, `cores`,
            and `cross_check_directory`.

    Returns:
        List of HMM IDs that could not be cross-checked.
    """
    logger.info("Cross check hit sequences with reference sequences")

    refseq_dir = config.paths.refseq
    data_dir = config.paths.data
    refseq_unavailable_list = []

    cross_check_dir = config.cross_check_directory
    intermediate_files = glob.glob(
        os.path.join(cross_check_dir, "*.intermediate_hits_faa")
    )
    diamond = find_executable("diamond")

    for inter_file in intermediate_files:
        logger.debug(
            f"Checking reference sequences for candidates sequences in {inter_file}"
        )
        hmm_id = os.path.splitext(os.path.basename(inter_file))[0].replace(
            ".intermediate_hits_faa", ""
        )
        hmm_type = hmm_id.split("_")[-1]

        db_base = os.path.splitext(os.path.join(refseq_dir, hmm_id))[0]
        db_path = db_base + ".dmnd"
        output_file = os.path.join(cross_check_dir, f"{hmm_id}.crosschecked.tsv")

        if os.path.isfile(output_file):
            logger.debug(f"Results already present for {hmm_id}")
            continue
        try:
            # 1. Try exact match for .dmnd
            if not os.path.isfile(db_path):
                # 2. Try exact match for .faa

                matches = find_file_by_name_partial(data_dir, f"{hmm_id}.faa")
                if not matches:
                    continue

                exact_faa = matches[0] if len(matches) > 1 else matches[0]

                if os.path.isfile(exact_faa):
                    faa_path = exact_faa
                    logger.debug(f"Found exact match: {exact_faa}")
                else:
                    # 3. Fallback: any file ending with {hmm_type}.faa
                    logger.warning(
                        f"Skipping {hmm_id}: Reference sequence file not found."
                    )
                    refseq_unavailable_list.append(hmm_id)
                    continue

                logger.debug(f"Creating Diamond DB from {faa_path} for {hmm_id}")
                subprocess.run(
                    [diamond, "makedb", "--in", faa_path, "-d", db_path, "--quiet"],
                    check=True,
                )

            # Run DIAMOND
            cmd = [
                diamond,
                "blastp",
                "--query",
                inter_file,
                "--db",
                db_path,
                "--out",
                output_file,
                "--outfmt",
                "6",
                "--max-target-seqs",
                "1",
                "--id",
                str(config.refseq_identity),
                "--threads",
                str(config.cores),
                "--quiet",
                "--"
                + str(
                    config.diamond_speed_mode
                ),  # Fastest mode for DIAMOND, only suitable for > 80 % identity
            ]
            logger.debug(f"Verifying {hmm_id} hits with reference sequences")
            subprocess.run(cmd)

            if os.path.getsize(output_file) == 0:
                os.remove(output_file)

        except Exception as e:
            logger.error(f"Failed to compare with diamond {hmm_id}\nError: {e}")
            refseq_unavailable_list.append(hmm_id)
            continue

    return refseq_unavailable_list


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

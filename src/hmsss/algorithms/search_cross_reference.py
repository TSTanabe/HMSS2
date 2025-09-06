#!/usr/bin/python
from __future__ import annotations

import os
import sys
import glob
import subprocess
import shlex
import shutil
import tempfile

from multiprocessing import Pool
from typing import Dict, List, Optional

from hmsss.cli.config import Config

from hmsss.core.logging import get_logger
from hmsss.utils.myUtil import find_executable

logger = get_logger(__name__)

"""
Cross-reference utilities for HMM hits.

Provides routines to:
- Concatenate many `.hmmreport` files into a single global report
  (via `xargs -0 cat` with a Python fallback).
- Parse score cutoff tables into dictionaries.
- Partition hits into trusted/intermediate/noise categories in parallel.
- Extract candidate protein FASTA per HMM and cross-check with reference
  sequences (DIAMOND).
- Promote validated candidates and summarize trusted hits.
"""

# Global shared variables
current_counter = None
counter_lock = None


########################################################################################################
#################### Filter the glob report to trusted hits and potential hits #########################
########################################################################################################


def concatenate_hmmreports_cat_xargs(
    report_paths: dict[str, str],
    output_path: str = "global_report.cat_hmmreport",
) -> str:
    """Concatenate multiple `.hmmreport` files into a single global file.

    Uses a NULL-delimited list with `xargs -0 cat` to avoid argv limits;
    falls back to a Python block-copy if the shell method fails.

    Args:
        report_paths: Mapping genome ID → hmmreport path.
        output_path: Destination file path.

    Returns:
        The `output_path` string.

    Raises:
        FileNotFoundError: If no valid input files are present.
        RuntimeError: If the shell concatenation fails unexpectedly.
    """
    logger.info(f"Concatenate raw hit reports to {output_path}")

    # Filter existing files
    valid_paths = [p for p in report_paths.values() if os.path.isfile(p)]
    if not valid_paths:
        logger.error("No valid hmmreport files found.")
        raise FileNotFoundError("No valid hmmreport files found.")

    # Create/overwrite output
    open_mode = "wb"
    with open(output_path, open_mode) as _:
        pass

    # Write a NULL-delimited file list to a temp file
    try:
        with tempfile.NamedTemporaryFile("wb", delete=False) as tf:
            tmp_list = tf.name
            for p in valid_paths:
                tf.write(p.encode("utf-8") + b"\x00")

        cmd = (
            f"xargs -0 -a {shlex.quote(tmp_list)} cat -- >> {shlex.quote(output_path)}"
        )

        completed = subprocess.run(cmd, shell=True)
        if completed.returncode != 0:
            raise RuntimeError(f"xargs/cat failed with code {completed.returncode}")

    except Exception as e:
        logger.warning(f"xargs+cat failed ({e}); falling back to Python copy.")
        # Fast block-wise fallback (no argv limits, nearly as fast as cat)
        with open(output_path, "wb") as outfp:
            for p in valid_paths:
                with open(p, "rb") as infp:
                    shutil.copyfileobj(infp, outfp, length=64 * 1024)
    finally:
        try:
            os.unlink(tmp_list)
        except Exception as e:
            logger.warning(f"Temporary concatenation file as not unlinked {e}")
            pass

        return output_path


def make_threshold_dict(
    file_path: str, threshold_type: int = 1, default_score: float = 50.0
) -> Dict[str, float]:
    """Parse a tab-separated cutoff table into a {hmm_id: threshold} dict.

    Columns are interpreted as different threshold types; values like `-inf`
    are mapped to a very high sentinel cutoff to mark them as unreachable.

    Args:
        file_path: Path to thresholds TSV file.
        threshold_type: Column (0-based) to read for the cutoff.
        default_score: Fallback score if the row lacks the chosen column.

    Returns:
        Mapping from HMM ID to its score threshold.
    """
    thresholds: Dict[str, float] = {}
    with open(file_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            parts = line.strip().split("\t")
            key = parts[0] if parts else None
            score = default_score

            try:
                if len(parts) == 2:
                    # Only one score present, use it
                    score = float(parts[1])
                elif len(parts) > threshold_type:
                    score = float(parts[threshold_type])
                    if parts[threshold_type] == "-inf":
                        score = float(
                            5000
                        )  # Hardcode never reachable score cutoff equal to infinite
            except (ValueError, IndexError) as e:
                logger.warning(
                    f"[Line {line_number}] Problem parsing: {line.strip()} — {e}"
                )
                continue

            if key:
                thresholds[key] = score

    return thresholds


def filter_trusted_and_noise_hits(
    config, glob_report: str, processes: int = 4
) -> str:
    """Split the global report into trusted and intermediate hits per HMM.

    Builds trusted/noise cutoff dicts from `options.score_threshold_file`,
    then processes each HMM ID in parallel to generate:
      - `{hmm}.trusted_hits`
      - `{hmm}.intermediate_hits`

    Args:
        config: Configuration with threshold table path and output directories.
        glob_report: Path to the concatenated global report.
        processes: Number of worker processes.

    Returns:
        Absolute path to the directory containing the filtered hit files.
    """
    trusted_dict = make_threshold_dict(
        config.score_threshold_file, 2, config.thrs_score
    )
    noise_dict = make_threshold_dict(
        config.score_threshold_file, 3, config.thrs_score
    )

    output_dir = config.cross_check_directory
    os.makedirs(output_dir, exist_ok=True)

    args = [
        (
            hmm_id,
            glob_report,
            trusted_dict[hmm_id],
            noise_dict.get(hmm_id, config.thrs_score),
            output_dir,
        )
        for hmm_id in trusted_dict
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(process_single_hmm, args)

    return os.path.abspath(output_dir)


def process_single_hmm(
    hmm_id: str,
    glob_report: str,
    trusted_cutoff: float,
    noise_cutoff: float,
    output_dir: str,
) -> None:
    """Extract hits for one HMM and split them by cutoff.

    Writes `{hmm_id}.trusted_hits` and, if applicable,
    `{hmm_id}.intermediate_hits`.

    Args:
        hmm_id: HMM profile identifier.
        glob_report: Path to the global report (`cat_hmmreport`).
        trusted_cutoff: Score threshold for trusted hits.
        noise_cutoff: Lower bound for intermediate hits.
        output_dir: Destination directory for outputs.
    """

    # Define output list for trusted and intermediate hits
    trusted_path = os.path.join(output_dir, f"{hmm_id}.trusted_hits")
    intermediate_path = os.path.join(output_dir, f"{hmm_id}.intermediate_hits")

    # Check if the trusted and the intermediate hit files are already present and skip existing files
    if os.path.isfile(intermediate_path):
        return
    if os.path.isfile(trusted_path):
        return

    # trusted hits + collect candidates
    candidates = {}
    with open(glob_report, "r") as infile, open(trusted_path, "w") as trusted_out:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")

            target = parts[0]
            hit_hmm = parts[3]

            if hit_hmm != hmm_id:
                continue

            try:
                if int(parts[10]) > 1:
                    # Domain Bit score is in column 13
                    score = float(parts[13])
                else:
                    # Full Bit score is in column 7 (0-based index)
                    score = float(parts[7])

            except ValueError:
                continue

            if score >= trusted_cutoff:
                logger.debug(
                    f"Above trusted cutoff hit {score} >= {trusted_cutoff} {parts[0]} {parts[3]}"
                )
                trusted_out.write(line)
            elif score > noise_cutoff:
                logger.debug(
                    f"Intermediate above noise cutoff hit {score} >= {trusted_cutoff} {parts[0]} {parts[3]}"
                )
                candidates[target] = score

    # remove empty trusted files
    if os.path.getsize(trusted_path) == 0:
        os.remove(trusted_path)
    # and skip if no candidates found
    if not candidates:
        return

    # Write down remaining candidates
    with open(glob_report, "r") as infile, open(intermediate_path, "w") as interm_out:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")

            target = parts[0]
            hit_hmm = parts[3]

            if hit_hmm == hmm_id and target in candidates:
                interm_out.write(line)

    if os.path.getsize(intermediate_path) == 0:
        os.remove(intermediate_path)
    return


def process_hitfile(
    hitfile_path: str,
    intermediate_hit_dir: str,
    faa_files: Dict[str, str],
    max_per_genome: int = 10,
) -> None:
    """Write per-HMM candidate sequences as FASTA from genome FAA files.

    Limits the number of sequences per genome to `max_per_genome`.
    Produces `{hmm_id}.intermediate_hits_faa`.

    Args:
        hitfile_path: Path to `{hmm_id}.intermediate_hits`.
        intermediate_hit_dir: Directory where outputs are written.
        faa_files: Mapping genome ID → FAA path.
        max_per_genome: Maximum sequences to extract per genome.
    """

    hmm_id = os.path.basename(hitfile_path).replace(".intermediate_hits", "")
    output_fasta = os.path.join(intermediate_hit_dir, f"{hmm_id}.intermediate_hits_faa")

    # Skip existing outputs
    if os.path.isfile(output_fasta):
        logger.debug(f"The above noise hit faa file already existed for {hmm_id}")
        return

    genome_hits = {}
    with open(hitfile_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            full_id = parts[0]
            if "___" not in full_id:
                continue
            genome_id, protein_id = full_id.split("___", 1)
            genome_hits.setdefault(genome_id, set()).add(protein_id)

    with open(output_fasta, "w") as out:
        for genome_id, protein_ids in genome_hits.items():
            faa_path = faa_files.get(genome_id)

            if not faa_path or not os.path.isfile(faa_path):
                logger.warning(f"FASTA not found for {genome_id} {faa_path}")

                continue

            written = 0  # ← pro genom zurücksetzen
            with open(faa_path, "r") as faa:
                write = False

                for line in faa:
                    if line.startswith(">"):
                        if max_per_genome is not None and written >= max_per_genome:
                            break  # Stop further writing for this genome
                        header_id = line[1:].split()[0]
                        write = header_id in protein_ids
                        if write:
                            out.write(f">{genome_id}___{header_id}\n")
                            written += 1
                    elif write:
                        out.write(line)


def generate_faa_per_hitfile_parallel(
    config: Config, intermediate_hit_dir: str, processes: int = 4
) -> None:
    """Extract candidate FASTA for all `.intermediate_hits` files in parallel.

    Args:
        config: Configuration with `faa_files` mapping.
        intermediate_hit_dir: Directory containing `.intermediate_hits`.
        processes: Number of parallel workers.
    """

    output_dir = intermediate_hit_dir  # same dir for output
    faa_files = config.faa_files  # dict: genome_id → path

    max_seqs_per_genome = config.max_seqs_per_genome

    hitfiles = [
        os.path.join(intermediate_hit_dir, f)
        for f in os.listdir(intermediate_hit_dir)
        if f.endswith(".intermediate_hits")
    ]

    args = [
        (hitfile, output_dir, faa_files, max_seqs_per_genome) for hitfile in hitfiles
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(process_hitfile, args)


##########################################################################################################################################################
#################### Cross check hits with reference sequences and add the hmmreport lines to the trusted cutoff hmmreport ###############################
##########################################################################################################################################################


def find_refseq_file(base_dir: str, filename: str) -> Optional[str]:
    """Find a file named `filename` under `base_dir` recursively.

    Args:
        base_dir: Search root directory.
        filename: Target filename.

    Returns:
        Full path if found, otherwise `None`.
    """

    for root, _, files in os.walk(base_dir):
        if filename in files:
            return os.path.join(root, filename)
    return None


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
                exact_faa = os.path.join(refseq_dir, f"{hmm_id}.faa")

                if os.path.isfile(exact_faa):
                    faa_path = exact_faa
                    logger.debug(f"Found exact match: {exact_faa}")
                else:
                    # 3. Fallback: any file ending with {hmm_type}.faa
                    logger.debug(
                        f"Exact match for {hmm_id} was not found. Now searching for {hmm_type}.faa"
                    )
                    pattern = os.path.join(refseq_dir, f"**/*{hmm_type}.faa")
                    backup_faa_files = glob.glob(pattern, recursive=True)

                    if backup_faa_files:
                        faa_path = backup_faa_files[0]
                        logger.debug(f"Using backup match: {faa_path}")
                    else:
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
            ]
            logger.info(f"Verifying {hmm_id} hits with reference sequences")
            subprocess.run(cmd)

            if os.path.getsize(output_file) == 0:
                os.remove(output_file)

        except Exception as e:
            logger.error(f"Failed to compare with diamond {hmm_id}\nError: {e}")
            refseq_unavailable_list.append(hmm_id)
            continue

    return refseq_unavailable_list


def process_crosscheck(hmm_id: str, crosscheck_dir: str) -> None:
    """Promote cross-validated candidates to trusted hits.

    Appends lines from `{hmm}.intermediate_hits` whose IDs occur in
    `{hmm}.crosschecked.tsv` to `{hmm}.trusted_hits`.

    Args:
        hmm_id: HMM profile identifier.
        crosscheck_dir: Directory with crosscheck and hit files.
    """
    crosscheck_path = os.path.join(crosscheck_dir, f"{hmm_id}.crosschecked.tsv")
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    if not os.path.exists(crosscheck_path):
        logger.error(f"Crosscheck file missing: '{crosscheck_path}'")
    if not os.path.exists(intermediate_path):
        logger.error(f"Intermediate file missing: '{intermediate_path}'")  #

    # Lade IDs aus crosscheck
    valid_hits = set()
    with open(crosscheck_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            valid_hits.add(line.strip().split("\t")[0])

    if not valid_hits:
        return

    promoted_count = 0

    # Promote direkt
    with open(intermediate_path, "r") as interm, open(trusted_path, "a") as trusted:
        for line in interm:
            if line.startswith("#") or not line.strip():
                continue
            if line.strip().split("\t")[0] in valid_hits:
                trusted.write(line)
                promoted_count += 1

    logger.info(
        f"{hmm_id}: promoted {promoted_count} hits via comparison with reference sequences"
    )


def promote_crosschecked_hits(crosscheck_dir: str, processes: int = 4) -> None:
    """Parallelize promotion of cross-validated candidates for all HMMs.

    Args:
        crosscheck_dir: Directory containing `.crosschecked.tsv` files.
        processes: Number of workers.
    """

    hmm_ids = [
        f.replace(".crosschecked.tsv", "")
        for f in os.listdir(crosscheck_dir)
        if f.endswith(".crosschecked.tsv")
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            process_crosscheck, [(hmm_id, crosscheck_dir) for hmm_id in hmm_ids]
        )


def summarize_trusted_hits(
    directory: str, crosscheck_dir: str, name: str, suffix: str = ".trusted_hits"
) -> str:
    """Concatenate all trusted hit files into a single summary.

    Args:
        directory: Output directory for the summary file.
        crosscheck_dir: Directory where `{hmm}.trusted_hits` live.
        name: Summary filename.
        suffix: Filename suffix to match (default: `.trusted_hits`).

    Returns:
        Path to the summary file.

    Raises:
        SystemExit: If no trusted hit files are found.
    """

    summary_path = os.path.join(directory, name)

    if os.path.isfile(summary_path) and os.path.getsize(summary_path) != 0:
        return summary_path

    trusted_files = [
        os.path.join(crosscheck_dir, f)
        for f in os.listdir(crosscheck_dir)
        if f.endswith(suffix)
    ]

    if not trusted_files:
        logger.error("No trusted hit files found for summary.")
        logger.info("There were no hits found in any genome. Closing the search")
        sys.exit()

    exit_code = os.system(f"cat {' '.join(trusted_files)} > {summary_path}")
    if exit_code != 0:
        logger.error("Could not summarize hit files")

    return summary_path


####################################################################################################
#################### Promote by cutoff when cross check is not available ###########################
####################################################################################################


def promote_by_cutoff(
    config: Config,
    directory: str,
    processes: int = 4,
    hmm_ids: Optional[List[str]] = None,
) -> None:
    """Parallel promotion of intermediate hits to trusted hits based on thresholds.

    Args:
        config (Config): Configuration with threshold file and type.
        directory (str): Directory with .intermediate_hits files.
        processes (int, optional): Number of worker processes.
        hmm_ids (List[str] or None): List of HMM IDs, or 'all' for all present.

    Returns:
        None:
        None

    Example:
        >>> promote_by_cutoff(config, '/tmp/xcheck', 4, hmm_ids=['PF00001', 'PF00002'])
        >>> promote_by_cutoff(config, '/tmp/xcheck', 4, hmm_ids=["all"])
    """
    # If no hmm identifier were defined use all that are in
    if hmm_ids == ["all"]:
        hmm_ids = [
            f.replace(".intermediate_hits", "")
            for f in os.listdir(directory)
            if f.endswith(".intermediate_hits")
        ]

    optimized_dict = make_threshold_dict(
        config.score_threshold_file, config.threshold_type, config.thrs_score
    )

    with Pool(processes=processes) as pool:
        pool.starmap(
            process_optimized_cutoff,
            [(hmm_id, directory, optimized_dict) for hmm_id in hmm_ids],
        )


def process_optimized_cutoff(
    hmm_id: str, crosscheck_dir: str, optimized_dict: Dict[str, float]
) -> None:
    """Promotes intermediate hits to trusted hits based on a score threshold.

    Args:
        hmm_id (str): HMM ID to process.
        crosscheck_dir (str): Directory with hit files.
        optimized_dict (Dict[str, float]): Dict of {hmm_id: threshold_score}.

    Returns:
        None

    Example:
        >>> process_optimized_cutoff('PF00001', '/tmp/xcheck', {'PF00001': 42.0})
    """

    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    threshold_score = optimized_dict.get(hmm_id, 50)

    if not os.path.isfile(intermediate_path):
        logger.warning(f"Intermediate or trusted hit file missing for {hmm_id}")
        return

    promoted_count = 0

    # Promote direkt
    with open(intermediate_path, "r") as interm, open(trusted_path, "a") as trusted:
        for line in interm:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            score = float(parts[7])
            if score >= threshold_score:
                trusted.write(line)
                promoted_count += 1

    logger.info(
        f"{hmm_id}: promoted {promoted_count} hits due to the given threshold {threshold_score}"
    )

#!/usr/bin/python
import gzip
import os

import shutil
import tempfile
import glob

from pathlib import Path
from typing import List, Tuple, Set, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from hmsss.cli.parse import file_path
from hmsss.utils.myUtil import get_all_files
from hmsss.utils import myUtil

log = myUtil.log


def queue_files(options) -> None:
    """
    Fills the options object with genome IDs, .faa and .gff file mappings.

    Args:
        options: options object with at least .fasta_file_directory, will be filled with:
            .queued_genomes (set[str])
            .faa_files (dict[str, str])
            .gff_files (dict[str, str])

    Operation:
        - Collect all zipped/unzipped protein fasta files and corresponding gff files.
        - Queue only if both files present, by genome identifier.

    Output Example:
        options.queued_genomes = {'GCF_000001405.39', ...}
        options.faa_files = {'GCF_000001405.39': '/dir/xxx.faa', ...}
        options.gff_files = {'GCF_000001405.39': '/dir/xxx.gff', ...}

    The hmmreport scheme: output = os.path.splitext(inputpath)[0] + ".hmmreport"
    """

    log.info("Filling the queue with faa files to be processed")
    genome_id_queue = set()
    fna_files = {}
    faa_files = {}
    gff_files = {}
    hmmreport_files = {}

    fna_gz_files = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".fna.gz"
    )
    faa_gz_files = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".faa.gz"
    )
    gff_gz_files = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".gff.gz"
    )

    fna_files = get_genome_id_files_dict(options.fasta_file_directory, extension=".fna")
    faa_files = get_genome_id_files_dict(options.fasta_file_directory, extension=".faa")
    gff_files = get_genome_id_files_dict(options.fasta_file_directory, extension=".gff")

    # Selektiere für Dekomprimierung
    decompress_targets = set()

    faa_all_genomes = set(faa_files.keys()) | set(faa_gz_files.keys())
    for gid, gz_path in fna_gz_files.items():
        if gid not in faa_all_genomes:
            decompress_targets.add(gz_path)

    # Regel 2: alle faa.gz und gff.gz entpacken, sofern noch keine entpackte Datei existiert
    for gid, gz_path in faa_gz_files.items():
        if gid not in faa_files:
            decompress_targets.add(gz_path)

    for gid, gz_path in gff_gz_files.items():
        if gid not in gff_files:
            decompress_targets.add(gz_path)

    log.info(f"Planned to decompress {len(decompress_targets)} file(s).")

    _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    fna_files: dict[str, str] = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".fna"
    )
    faa_files: dict[str, str] = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".faa"
    )
    gff_files: dict[str, str] = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".gff"
    )
    hmmreport_files: dict[str, str] = get_genome_id_files_dict(
        options.fasta_file_directory, extension=".hmmreport"
    )

    # Entferne die fna files die bereits einen faa file haben
    for gid in list(fna_files.keys()):
        if gid in faa_files:
            del fna_files[gid]

    # Definiere die faa/gff paare
    common_ids = set(gff_files.keys()) & set(faa_files.keys())
    for gid in list(faa_files.keys()):
        if gid not in common_ids:
            del faa_files[gid]
    for gid in list(gff_files.keys()):
        if gid not in common_ids:
            del gff_files[gid]
    for gid in list(hmmreport_files.keys()):
        if gid not in common_ids:
            del hmmreport_files[gid]

    # Queue mit den validen Genome-IDs füllen
    genome_id_queue = common_ids

    options.queued_genomes = genome_id_queue
    options.fna_files = fna_files
    options.faa_files = faa_files
    options.gff_files = gff_files
    options.hmmreport_files = hmmreport_files
    log.info(f"Found {len(fna_files)} fna files for transcription.")
    log.info(f"Queued {len(options.queued_genomes)} faa/gff pairs")
    log.info(f"Found {len(hmmreport_files)} existing hmmreports for faa/gff pairs.")
    return


def get_all_files_with_extension(directory: str, extension: str) -> Set[str]:
    """
    Recursively find all files with a given extension inside a directory.

    Args:
        directory (str): Root directory to search.
        extension (str): File extension to look for (with or without leading dot).

    Returns:
        Set[str]: Absolute file paths of matching files.
    """
    root = Path(directory)
    if not root.is_dir():
        raise ValueError(f"Not a directory: {directory}")

    ext = extension if extension.startswith(".") else f".{extension}"

    result: Set[str] = set()
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname.endswith(ext):
                result.add(str(Path(dirpath) / fname))

    return result


def get_genome_id_files_dict(directory: str, extension: str) -> Dict[str, str]:
    """
    Returns a dictionary mapping genome identifiers to file paths.
        Args: directory (str): Root directory to search.
        extension (str): File extension to filter by.
    """

    genome_id_files = {}
    file_path_set = get_all_files_with_extension(directory, extension)

    for file_path in file_path_set:
        genome_id = myUtil.get_genome_id(file_path)
        genome_id_files[genome_id] = file_path

    return genome_id_files


def unpackgz(path: str) -> str:
    """
    Decompresses a .gz file if not already extracted.
    """
    if not path.endswith(".gz"):
        return path
    file = path[:-3]
    if os.path.exists(file):
        return file
    with gzip.open(path, "rb") as f_in:
        with open(file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)  # type: ignore[arg-type]
    return file


# --- kleine Helper-Routine: entpackt nur .gz, sonst no-op ---
def _decompress_gz_only(path: str) -> str:
    """
    Entpackt NUR .gz-Dateien (no-op für andere Pfade).
    Gibt den Pfad zur entpackten Datei zurück.
    """
    try:
        return unpackgz(path)
    except Exception as e:
        log.error(f"Failed to decompress '{path}': {e}")
        return path


def _parallel_decompress(paths, max_workers: Optional[int] = None) -> None:
    """
    Entpackt eine Menge von .gz-Dateien parallel. No-op für Nicht-.gz.
    Keine Rückgabewerte
    """
    if not paths:
        return
    max_workers = max_workers or 4
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(_decompress_gz_only, p): p for p in paths}
        for fut in as_completed(futures):
            _ = fut.result()  # Fehlerlogging passiert in _decompress_gz_only


def concatenate_selected_hmms(
    src_dir: str,
    allowed_words: List[str],
    prefix: str,
    suffix: str,
    output_library: str,
) -> None:
    """
    Concatenate .hmm files from subdirectories where at least one word in the dir name (split by '_')
    is present in allowed_words list.

    Args:
        src_dir (str): Parent directory to search (e.g., __location__ + "/data").
        allowed_words (List[str]): List of allowed words (from whitespace-separated user input).
        prefix (str): File prefix filter (e.g., 'grp').
        suffix (str): File suffix filter (e.g., '.hmm').
        output_library (str): Output concatenated library file path.

    Returns:
        None:
    """

    files_to_concatenate = []

    for subdir, dirs, files in os.walk(src_dir):
        subdir_name = os.path.basename(subdir)
        subdir_words = set(subdir_name.split("_"))
        if subdir_words & set(allowed_words):
            # If intersection is non-empty, at least one word matches
            matched_files = glob.glob(os.path.join(subdir, f"{prefix}*{suffix}"))
            files_to_concatenate.extend(matched_files)

    # Concatenate files
    with open(output_library, "w") as outfile:
        for fname in files_to_concatenate:
            with open(fname) as infile:
                shutil.copyfileobj(infile, outfile)


def concatenate_files_shell(
    search_directory: str,
    allowed_prefix: str,
    allowed_suffix: str,
    output_file_path: str,
) -> None:
    """
    Rekursiv alle Files suchen, die mit allowed_prefix beginnen und allowed_suffix enden.
    Diese zusammenführen und als output_file_path speichern.
    """
    matched_files = []
    for root, _, files in os.walk(search_directory):
        for fname in files:
            if fname.startswith(allowed_prefix) and fname.endswith(allowed_suffix):
                matched_files.append(os.path.join(root, fname))

    if matched_files:
        os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
        cat_command = (
            "cat "
            + " ".join(f'"{f}"' for f in matched_files)
            + f' > "{output_file_path}"'
        )
        # logger.debug(f"Running: {cat_command}")
        os.system(cat_command)
        log.debug(f"Concatenated {len(matched_files)} files into {output_file_path}")
    else:
        log.error(f"No matching files found for concatenation in {search_directory}.")


def format_pattern_files_inplace(filename: str, prefix: str, suffix: str):
    tmpfile = tempfile.NamedTemporaryFile("w", delete=False)
    with open(filename, "r") as fin, tmpfile:
        for i, line in enumerate(fin, 1):
            new_line = f"{prefix}{i}{suffix} {line.rstrip()}"
            new_line = new_line.replace(" ", "\t")
            tmpfile.write(new_line + "\n")
    shutil.move(tmpfile.name, filename)

#!/usr/bin/python
import gzip
import os

import shutil
import tempfile
import glob

from pathlib import Path
from typing import List, Set, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from hmsss.core.logging import get_logger
from hmsss.utils import myUtil

log = get_logger(__name__)


def queue_fna_inputs(config) -> dict[str, str]:
    """
    Sammelt ausschließlich FNA-Inputs für die spätere Translation.
    - Entpackt .fna.gz NUR dann, wenn für die GenomeID KEIN .faa/.faa.gz existiert.
    - Entfernt alle FNA-Einträge, für die bereits ein FAA existiert (gz oder ungezipped).

    Setzt auf `options`:
      .fna_files (dict[genome_id -> .fna-Pfad])
    """
    root = config.fasta_file_directory

    # Aktuelle Lage erfassen
    fna_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna.gz")
    fna_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna")
    faa_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")

    faa_all_genomes: Set[str] = set(faa_files) | set(faa_gz_files)

    # .fna.gz nur entpacken, wenn (noch) kein FAA vorhanden ist
    decompress_targets: Set[str] = set()
    for gid, gz_path in fna_gz_files.items():
        if gid not in faa_all_genomes:
            # Ungezippte .fna fehlt oder wir wollen sicherstellen, dass sie da ist
            if gid not in fna_files:
                decompress_targets.add(gz_path)

    if decompress_targets:
        log.info(f"Planned to decompress {len(decompress_targets)} .fna file(s).")
        _parallel_decompress(decompress_targets, getattr(config, "cores", 4))

    # Nach evtl. Entpacken erneut einlesen
    fna_files = get_genome_id_files_dict(root, extension=".fna")

    # ALLE FNA entfernen, wenn FAA bereits existiert
    for gid in list(fna_files.keys()):
        if gid in faa_all_genomes:
            del fna_files[gid]

    config.fna_files = fna_files
    log.info(f"Found {len(fna_files)} fna files for translation.")
    return fna_files


def queue_protein_annotation_inputs(options) -> None:
    """
    Behandelt FAA/GFF/HMMREPORT vollständig getrennt vom FNA-Teil.
    - Entpackt .faa.gz/.gff.gz sofern ungezippte Pendants fehlen.
    - Bildet gültige Paare (GenomeIDs mit BOTH: .faa UND .gff).
    - Filtert .hmmreport auf diese Paare.

    Setzt auf `options`:
      .queued_genomes (set[str])
      .faa_files (dict[genome_id -> .faa-Pfad])
      .gff_files (dict[genome_id -> .gff-Pfad])
      .hmmreport_files (dict[genome_id -> .hmmreport-Pfad])
    """
    root = options.fasta_file_directory

    # Aktuelle Lage erfassen
    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")
    gff_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff.gz")

    faa_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    gff_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff")

    # Entpacken planen: .faa.gz / .gff.gz nur wenn das ungezippte Pendant fehlt
    decompress_targets: Set[str] = set()
    for gid, gz_path in faa_gz_files.items():
        if gid not in faa_files:
            decompress_targets.add(gz_path)
    for gid, gz_path in gff_gz_files.items():
        if gid not in gff_files:
            decompress_targets.add(gz_path)

    if decompress_targets:
        log.info(f"[ANN] Planned to decompress {len(decompress_targets)} file(s).")
        _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    # Nach evtl. Entpacken erneut einlesen
    faa_files = get_genome_id_files_dict(root, extension=".faa")
    gff_files = get_genome_id_files_dict(root, extension=".gff")
    hmmreport_files: Dict[str, str] = get_genome_id_files_dict(
        root, extension=".hmmreport"
    )

    # Nur GenomeIDs behalten, die FAA UND GFF haben
    common_ids: Set[str] = set(faa_files) & set(gff_files)

    # Dictionaries auf common_ids beschränken
    faa_files = {gid: path for gid, path in faa_files.items() if gid in common_ids}
    gff_files = {gid: path for gid, path in gff_files.items() if gid in common_ids}
    hmmreport_files = {
        gid: path for gid, path in hmmreport_files.items() if gid in common_ids
    }

    options.queued_genomes = common_ids
    options.faa_files = faa_files
    options.gff_files = gff_files
    options.hmmreport_files = hmmreport_files

    log.info(f"Queued {len(common_ids)} faa/gff pairs.")
    log.info(f"Found {len(hmmreport_files)} existing hmmreports for faa/gff pairs.")


def queue_faa_without_gff(options) -> dict[str, str]:
    """
    Sammelt alle FAA-Files, für die KEIN GFF existiert (weder .gff noch .gff.gz).
    - Entpackt .faa.gz vor Aufnahme in die Queue (falls .faa fehlt).

    Setzt auf `options`:
      .missing_gff_genomes : set[str]
      .faa_missing_gff     : dict[str, str]  # genomeID -> Pfad zu entpacktem .faa
    """
    root = options.fasta_file_directory

    # Ist-Zustand erfassen
    faa_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")
    gff_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff")
    gff_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff.gz")

    # GenomeIDs mit vorhandenen GFFs (gezipped oder ungezipped)
    genomes_with_any_gff: Set[str] = set(gff_files) | set(gff_gz_files)

    # Alle GenomeIDs, für die es FAA (gezipped/ungezipped) gibt
    genomes_with_any_faa: Set[str] = set(faa_files) | set(faa_gz_files)

    # Ziel: FAA ohne GFF
    missing_gff_genomes: Set[str] = genomes_with_any_faa - genomes_with_any_gff

    # .faa.gz für diese Ziel-Genome entpacken, falls .faa fehlt
    decompress_targets: Set[str] = set()
    for gid in missing_gff_genomes:
        if gid not in faa_files and gid in faa_gz_files:
            decompress_targets.add(faa_gz_files[gid])

    if decompress_targets:
        log.info(f"Planned to decompress {len(decompress_targets)} file(s).")
        _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    # Nach Entpacken FAA erneut einlesen
    faa_files = get_genome_id_files_dict(root, extension=".faa")

    # Queue-Dict: nur GenomeIDs, die (jetzt) ein ungezippetes .faa haben und weiterhin kein GFF
    faa_missing_gff: Dict[str, str] = {
        gid: faa_files[gid]
        for gid in missing_gff_genomes
        if gid in faa_files  # sicherstellen, dass ungezippte FAA existiert
    }

    log.info(f"Queued {len(faa_missing_gff)} faa files without gff for transcription.")
    return faa_missing_gff


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

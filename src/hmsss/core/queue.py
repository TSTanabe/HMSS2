#!/usr/bin/python
import gzip
import os
import re

import shutil
import tempfile
import glob

from pathlib import Path
from typing import List, Set, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from hmsss.core.logging import get_logger
from hmsss.utils import myUtil

log = get_logger(__name__)

"""
Queue and preprocessing utilities for genome input files.

This module handles:
- Collection and filtering of FNA/FAA/GFF input files by genome ID.
- Parallel decompression of `.gz` archives when needed.
- Construction of dictionaries mapping genome IDs to file paths.
- Concatenation of HMM libraries, cutoff files, co-occurrence, and patterns.
- Formatting of pattern/co-occurrence files with systematic prefixes.
"""


def queue_fna_inputs(config) -> dict[str, str]:
    """Collect FNA inputs for translation.

    - Decompress `.fna.gz` files only if no corresponding `.faa/.faa.gz` exists.
    - Remove all FNA entries if a corresponding FAA already exists.

    Args:
        config: Configuration object with `.fasta_file_directory` and `.cores`.

    Returns:
        Mapping genome ID → `.fna` path. Also stored in `config.fna_files`.

    Side Effects:
        May decompress `.fna.gz` files in parallel.
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


def queue_protein_annotation_inputs(config) -> None:
    """Collect FAA/GFF/HMMREPORT inputs, independent of FNA.

    - Decompress `.faa.gz` and `.gff.gz` if uncompressed counterparts are missing.
    - Keep only genome IDs that have both FAA and GFF.
    - Restrict HMMREPORT files to these genome IDs.

    Args:
        config: Object with `.fasta_file_directory` and `.cores`.

    Side Effects:
        Sets attributes on `config`:
          - `queued_genomes` (set of genome IDs)
          - `faa_files` (dict genome ID → FAA path)
          - `gff_files` (dict genome ID → GFF path)
          - `hmmreport_files` (dict genome ID → HMMREPORT path)
    """
    root = config.fasta_file_directory

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
        log.info(f"Planned to decompress {len(decompress_targets)} .gz file(s).")
        _parallel_decompress(decompress_targets, getattr(config, "cores", None))
    else:
        log.info(f"All files are already decompressed.")
    # Collect the faa and gff files to the dictionary
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

    config.queued_genomes = common_ids
    config.faa_files = faa_files
    config.gff_files = gff_files
    config.hmmreport_files = hmmreport_files

    log.info(f"Queued {len(common_ids)} faa/gff pairs.")
    log.info(f"Found {len(hmmreport_files)} existing hmmreports for faa/gff pairs.")


def queue_faa_without_gff(config) -> dict[str, str]:
    """Collect FAA files for which no GFF exists.

    - Includes both `.gff` and `.gff.gz` in the check.
    - Decompress `.faa.gz` for target genomes if `.faa` is missing.

    Args:
        config: Object with `.fasta_file_directory` and `.cores`.

    Returns:
        Mapping genome ID → `.faa` path for genomes without GFF.

    Side Effects:
        Sets on `options`:
          - `missing_gff_genomes` (set of genome IDs)
          - `faa_missing_gff` (dict genome ID → FAA path)
    """
    root = config.fasta_file_directory

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
        _parallel_decompress(decompress_targets, getattr(config, "cores", None))

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

def queue_graftm_fastq_inputs(config) -> dict[str, str]:
    """
    Sucht im FASTA/FASTQ-Inputverzeichnis und gibt ein Dictionary zurück:

        { genome_id : full_path }

    Falls genomeID sowohl in .fastq als auch in .fastq.gz vorkommt:
        -> .fastq.gz gewinnt
    """

    root = config.fasta_file_directory

    # Sammeln aller FASTQ-Dateien
    fastq: Dict[str, str] = get_genome_id_files_dict(root, extension=".fastq")
    fastq_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".fastq.gz")

    # Merge:
    # zuerst .fastq, dann .fastq.gz drüber -> gz überschreibt fastq bei Konflikt
    merged: Dict[str, str] = {}
    merged.update(fastq)
    merged.update(fastq_gz)

    # Speichern in Config
    config.graftm_fastq_files = merged

    return merged


def get_all_files_with_extension(directory: str, extension: str) -> Set[str]:
    """Recursively find all files with a given extension.

    Args:
        directory: Root directory to search.
        extension: File extension to filter by (with or without leading dot).

    Returns:
        Set of absolute file paths.
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
    """Map genome IDs to file paths for a given extension.

    Args:
        directory: Root directory to search.
        extension: File extension to filter by.

    Returns:
        Dictionary genome ID → file path.
    """

    genome_id_files = {}
    file_path_set = get_all_files_with_extension(directory, extension)

    for file_path in file_path_set:
        genome_id = myUtil.get_genome_id(file_path)
        genome_id_files[genome_id] = file_path

    return genome_id_files


def unpackgz(path: str) -> str:
    """Decompress a `.gz` file if not already extracted.

    Args:
        path: Path to a `.gz` file or plain file.

    Returns:
        Path to the decompressed file. If input is not `.gz`, returns unchanged.
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
    """Decompress `.gz` files only (no-op for others).

    Args:
        path: File path.

    Returns:
        Path to decompressed file or unchanged path if not `.gz`.
    """
    try:
        return unpackgz(path)
    except Exception as e:
        log.error(f"Failed to decompress '{path}': {e}")
        return path


def _parallel_decompress(paths, max_workers: Optional[int] = None) -> None:
    """Decompress a set of `.gz` files in parallel.

    Args:
        paths: Iterable of file paths.
        max_workers: Number of worker processes (default: 4).

    Returns:
        None. Decompressed files are written in place.
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
    """Concatenate HMM files filtered by directory name.

    - Searches subdirectories under `src_dir`.
    - Includes only subdirectories where the name shares a token with `allowed_words`.
    - Concatenates all files matching prefix+suffix pattern.

    Args:
        src_dir: Parent directory to search.
        allowed_words: Tokens (from user input) to match against directory names.
        prefix: Required file prefix (e.g., "grp").
        suffix: Required file suffix (e.g., ".hmm").
        output_library: Path to write the concatenated library.
    """

    files_to_concatenate = []

    for subdir, dirs, files in os.walk(src_dir):
        subdir_name = os.path.basename(subdir)
        subdir_words = set(subdir_name.split("_"))
        if subdir_words & set(allowed_words):
            # If intersection is non-empty, at least one word matches
            # Finds everything below the subdir
            matched_files = glob.glob(os.path.join(subdir, "**", f"{prefix}*{suffix}"))
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
    """Concatenate files matching prefix and suffix into a single output file.

    - Searches recursively in `search_directory`.
    - Concatenates all files that start with `allowed_prefix` and end with `allowed_suffix`.

    Args:
        search_directory: Directory to search.
        allowed_prefix: Required filename prefix.
        allowed_suffix: Required filename suffix.
        output_file_path: Path to write the concatenated output.
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
        log.error(
            f"No matching files found for concatenation in {search_directory} with prefix '{allowed_prefix}' and suffix '{allowed_suffix}'."
        )


def concatenate_hmms_from_selected_metabolism_packages(
    src_dir: str,
    allowed_words: list[str],
    output_library: str,
) -> None:
    """
    Concatenate .hmm Dateien NUR aus den Paketen direkt unter `src_dir`,
    deren Paketname (tokenisiert) mindestens ein Token mit `allowed_words` teilt.
    Innerhalb eines Pakets werden ausschließlich Unterordner berücksichtigt,
    deren Name "HMMs" (case-insensitive) enthält. Darunter wird rekursiv nach *.hmm gesucht.

    Beispielstruktur:
        data/
          v8_Sulfur/
            v8_HMMs_sqr_dsr/
              grpX.hmm, ...
          v8_Methan/
            meta3_meta4/

    allowed_words = ["Sulfur", "Nitrogen"]  -> nimmt nur v8_Sulfur und v8_Nitrogen
    """

    def _tokenize(name: str) -> set[str]:
        """Tokenisiere Ordner-/Set-Namen robust (Unterstriche, Bindestriche, Ziffernblöcke)."""
        return {t.lower() for t in re.findall(r"[A-Za-z0-9]+", name)}

    base = Path(src_dir)
    allowed_tokens = set().union(*(_tokenize(w) for w in allowed_words))
    files_to_concat: list[str] = []

    if not base.is_dir():
        raise ValueError(f"Not a directory: {src_dir}")

    # Nur Top-Level-Pakete (direkte Unterordner von src_dir)
    for pkg in sorted(p for p in base.iterdir() if p.is_dir()):
        pkg_tokens = _tokenize(pkg.name)
        if not (pkg_tokens & allowed_tokens):
            continue  # Paket nicht ausgewählt

        # Nur Ordner innerhalb des Pakets, deren Name "HMMs" enthält
        hmm_dirs: list[Path] = []
        for root, dirs, _files in os.walk(pkg):
            for d in dirs:
                if "hmms" in d.lower():
                    hmm_dirs.append(Path(root) / d)

        # Rekursiv *.hmm aus allen passenden HMM-Ordnern einsammeln
        for hmm_dir in hmm_dirs:
            for path in hmm_dir.rglob("*.hmm"):
                if path.is_file():
                    files_to_concat.append(str(path))

    # deterministische Reihenfolge
    files_to_concat = sorted(set(files_to_concat))

    if not files_to_concat:
        log.warning(
            "No HMM files found for selected packages in '%s' with hmm_sets=%s",
            src_dir,
            allowed_words,
        )
        # absichtlich KEIN Schreiben -> ressource_preparation fällt dann auf den
        # globalen Fallback zurück, falls konfiguriert.
        return

    os.makedirs(os.path.dirname(output_library) or ".", exist_ok=True)
    with open(output_library, "w") as out:
        for fp in files_to_concat:
            with open(fp, "r") as fin:
                shutil.copyfileobj(fin, out)
    log.info("Concatenated %d HMMs into %s", len(files_to_concat), output_library)


def format_pattern_files_inplace(filename: str, prefix: str, suffix: str):
    """Rewrite a pattern/co-occurrence file with systematic prefixes.

    - Adds line numbers as IDs with given prefix/suffix.
    - Replaces spaces with tabs.

    Args:
        filename: Path to the file to reformat.
        prefix: Prefix string for IDs (e.g., "dsb-").
        suffix: Suffix string for IDs (e.g., "_").
    """
    tmpfile = tempfile.NamedTemporaryFile("w", delete=False)
    with open(filename, "r") as fin, tmpfile:
        for i, line in enumerate(fin, 1):
            new_line = f"{prefix}{i}{suffix} {line.rstrip()}"
            new_line = new_line.replace(" ", "\t")
            tmpfile.write(new_line + "\n")
    shutil.move(tmpfile.name, filename)

#!/usr/bin/python
import gzip
import os
import re

import shutil
import tempfile
import glob

from pathlib import Path
from typing import List, Set, Dict, Optional

from hmsss.core.logging import get_logger
from hmsss.db import database
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
    """Collect FNA inputs for translation (NO decompression).

    Rules:
      - If both .fna and .fna.gz exist for the same genomeID -> .fna wins
      - If any FAA exists (.faa or .faa.gz) for the same genomeID -> remove FNA from queue

    Returns:
        Mapping genome ID → chosen FNA path (.fna preferred over .fna.gz).
        Also stored in config.fna_files.
    """
    root = config.fasta_file_directory

    # collect files (no decompression)
    fna_plain: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna")
    fna_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna.gz")

    faa_plain: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")

    # genomes that already have protein FASTA -> skip FNA translation
    faa_all_genomes: Set[str] = set(faa_plain) | set(faa_gz)

    # merge: gz first, then plain (plain wins)
    fna_files: Dict[str, str] = {}
    fna_files.update(fna_gz)
    fna_files.update(fna_plain)

    # drop genomes where FAA already exists
    for gid in list(fna_files.keys()):
        if gid in faa_all_genomes:
            del fna_files[gid]

    config.fna_files = fna_files
    log.info(f"Found {len(fna_files)} fna files for translation.")
    return fna_files


def queue_protein_annotation_inputs(config) -> None:
    """
    Collect FAA/GFF inputs for annotation (gz + plain), WITHOUT decompression.

    Rule:
      - If both plain and gz exist for the same genomeID → plain WINS
      - Genome is queued only if BOTH FAA and GFF exist (gz or plain)

    Side effects:
      - sets config.queued_genomes
      - sets config.faa_files
      - sets config.gff_files
      - sets config.hmmreport_files (legacy)
    """
    root = config.fasta_file_directory

    # --- collect files ---
    faa_plain = get_genome_id_files_dict(root, extension=".faa")
    faa_gz = get_genome_id_files_dict(root, extension=".faa.gz")

    gff_plain = get_genome_id_files_dict(root, extension=".gff")
    gff_gz = get_genome_id_files_dict(root, extension=".gff.gz")

    # --- merge: gz first, then plain (plain wins) ---
    faa_files = {}
    faa_files.update(faa_gz)
    faa_files.update(faa_plain)

    gff_files = {}
    gff_files.update(gff_gz)
    gff_files.update(gff_plain)

    # legacy / optional
    hmmreport_files = get_genome_id_files_dict(root, extension=".hmmreport")

    # --- keep only genomeIDs that have BOTH faa and gff ---
    common_ids = set(faa_files) & set(gff_files)

    faa_files = {gid: faa_files[gid] for gid in common_ids}
    gff_files = {gid: gff_files[gid] for gid in common_ids}
    hmmreport_files = {
        gid: path for gid, path in hmmreport_files.items() if gid in common_ids
    }

    config.queued_genomes = list(common_ids)
    config.faa_files = faa_files
    config.gff_files = gff_files
    config.hmmreport_files = hmmreport_files

    log.info(f"Queued {len(common_ids)} genomes with FAA and GFF.")


def queue_faa_without_gff(config) -> dict[str, str]:
    """Collect FAA files for which no GFF exists (NO decompression).

    Rules:
      - If .faa and .faa.gz both exist -> .faa wins
      - If only .faa.gz exists -> use .faa.gz
      - If any GFF exists (.gff or .gff.gz) -> skip genome

    Returns:
        Mapping genome ID → FAA path (.faa preferred over .faa.gz)
    """
    root = config.fasta_file_directory

    # Ist-Zustand erfassen
    faa_plain: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")

    gff_plain: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff")
    gff_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff.gz")

    # GenomeIDs mit vorhandenen GFFs
    genomes_with_any_gff: Set[str] = set(gff_plain) | set(gff_gz)

    # GenomeIDs mit FAA (egal ob gezippt)
    genomes_with_any_faa: Set[str] = set(faa_plain) | set(faa_gz)

    # Ziel: FAA ohne GFF
    missing_gff_genomes: Set[str] = genomes_with_any_faa - genomes_with_any_gff

    # Merge FAA: gz zuerst, dann plain (plain gewinnt)
    faa_missing_gff: Dict[str, str] = {}
    for gid in missing_gff_genomes:
        if gid in faa_gz:
            faa_missing_gff[gid] = faa_gz[gid]
        if gid in faa_plain:
            faa_missing_gff[gid] = faa_plain[gid]

    log.info(f"Queued {len(faa_missing_gff)} faa files without gff for transcription.")
    return faa_missing_gff


def queue_read_mapping_fastq_inputs(config) -> dict[str, str]:
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
    config.fastq_files = merged

    return merged


def queue_read_mapping_faa_inputs(config) -> dict[str, str]:
    """
    Search the FASTA input directory for protein sequence files and return:

        { genome_id : full_path }

    If both .faa and .faa.gz exist for the same genome_id:
        -> .faa.gz wins
    """

    root = config.fasta_file_directory

    faa: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")

    merged: Dict[str, str] = {}
    merged.update(faa)
    merged.update(faa_gz)

    config.faa_files = merged
    return merged


def queue_read_mapping_fna_inputs(config) -> dict[str, str]:
    """
    Search the FASTA input directory for nucleotide genome files and return:

        { genome_id : full_path }

    If both .fna and .fna.gz exist for the same genome_id:
        -> .fna.gz wins
    """

    root = config.fasta_file_directory

    fna: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna")
    fna_gz: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna.gz")

    merged: Dict[str, str] = {}
    merged.update(fna)
    merged.update(fna_gz)

    config.fna_files = merged
    return merged


def remove_genomes_already_in_db_from_queue(config) -> int:
    """
    Entfernt Genome aus config.queued_genomes, die bereits in der DB existieren.
    Passt außerdem die zugehörigen Dicts (faa_files, gff_files, hmmreport_files) an.

    Returns:
        Anzahl der entfernten Genome.
    """
    if not getattr(config, "database_directory", None):
        return 0

    # falls DB noch nicht existiert -> nichts entfernen
    if not os.path.isfile(config.database_directory):
        return 0

    existing: Set[str] = database.fetch_genome_ids(
        config.database_directory
    )  # :contentReference[oaicite:3]{index=3}

    queued = list(getattr(config, "queued_genomes", []))
    if not queued:
        return 0

    # welche sollen raus?
    to_remove = {gid for gid in queued if gid in existing}
    if not to_remove:
        return 0

    # queued_genomes filtern (Reihenfolge beibehalten)
    config.queued_genomes = [gid for gid in queued if gid not in to_remove]

    # zugehörige Mappings filtern, falls vorhanden
    for attr in ("faa_files", "gff_files", "hmmreport_files"):
        d = getattr(config, attr, None)
        if isinstance(d, dict) and d:
            for gid in list(to_remove):
                d.pop(gid, None)

    return len(to_remove)


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


def collect_gpkg_from_selected_metabolism_packages(
    src_dir: str,
    allowed_words: list[str],
) -> dict[str, str]:
    """
    Return a FLAT dictionary mapping:

        { gpkg_basename : full_path_to_gpkg }

    Only metabolism packages (directories directly under src_dir) whose
    name shares at least one token with allowed_words are inspected.
    Inside each selected metabolism package, subdirectories whose name
    contains 'Gpkg' (case-insensitive) are searched recursively for *.gpkg.

    Example return:
        {
            "dsr_core": "/.../v8_Gpkg_Dsr_sHdr_Sox/dsr_core.gpkg",
            "sox_cluster": "/.../v8_Gpkg_SulfurOxidation/sox_cluster.gpkg",
        }
    """

    def _tokenize(name: str) -> set[str]:
        return {t.lower() for t in re.findall(r"[A-Za-z0-9]+", name)}

    base = Path(src_dir)
    if not base.is_dir():
        raise ValueError(f"Not a directory: {src_dir}")

    allowed_tokens = set().union(*(_tokenize(w) for w in allowed_words))

    gpkg_map: dict[str, str] = {}

    # iterate top-level metabolism packages
    for pkg in sorted(p for p in base.iterdir() if p.is_dir()):
        pkg_tokens = _tokenize(pkg.name)
        if not (pkg_tokens & allowed_tokens):
            continue

        # find directories containing "Gpkg" anywhere in their name
        gpkg_dirs: list[Path] = []
        for root, dirs, _files in os.walk(pkg):
            for d in dirs:
                if "Gpkg" in d.lower() or "gpkg" in d.lower():
                    gpkg_dirs.append(Path(root) / d)

        # now collect *.gpkg files from all gpkg dirs
        for gdir in gpkg_dirs:
            for path in gdir.rglob("*.gpkg"):
                if path.is_dir():
                    key = path.stem  # filename without ".gpkg"
                    value = str(path)
                    gpkg_map[key] = value

    return gpkg_map


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

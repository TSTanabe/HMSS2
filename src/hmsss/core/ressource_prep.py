from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

from hmsss.core.logging import get_logger, print_header
from hmsss.cli.paths import DATA_DIR, GPKG_DIR
from hmsss.core import queue as queue

log = get_logger(__name__)

"""
Resource preparation stage for HMSSS.

This stage sets up result subdirectories, ensures availability of HMM
libraries, thresholds, patterns, and co-occurrence files, and initializes
the cross-check directory. If resources are missing, they are concatenated
from the bundled `DATA_DIR`. Raises `SystemExit` if required files are not found.
"""


def _ensure_dir(p: str | os.PathLike) -> str:
    """Create a directory if it does not exist.

    Args:
        p: Path to the directory.

    Returns:
        Path as string.
    """
    Path(p).mkdir(parents=True, exist_ok=True)
    return str(p)


def _require_path_exists(p: str, desc: str) -> None:
    """Ensure that a required path exists.

    Args:
        p: Path to check.
        desc: Human-readable description of the resource.

    Raises:
        SystemExit: If the path does not exist.
    """
    if not Path(p).exists():
        log.error("%s not found: %s", desc, p)
        raise SystemExit(f"Missing required resource: {desc} -> {p}")


def concatenate_hmms_from_selected_metabolism_packages(
    src_dir: str,
    allowed_words: list[str],
    output_library: str,
    allowed_versions: set[str] | None = None,
) -> None:
    """
    ... wie zuvor ...

    Zusätzlich:
    - `allowed_versions` filtert nach Paket-Version = Prefix bis zum ersten "_" im
      Top-Level-Ordnernamen unter src_dir (z.B. "v8_Sulfur" -> "v8").
    - Wenn allowed_versions None oder leer: kein Versionsfilter.
    """

    def _tokenize(name: str) -> set[str]:
        return {t.lower() for t in re.findall(r"[A-Za-z0-9]+", name)}

    def _pkg_version(pkg_name: str) -> str | None:
        # "v8_Sulfur" -> "v8"; wenn kein "_" vorhanden, keine Version ableitbar
        if "_" not in pkg_name:
            return None
        return pkg_name.split("_", 1)[0]

    base = Path(src_dir)
    allowed_tokens = set().union(*(_tokenize(w) for w in allowed_words))
    files_to_concat: list[str] = []

    if not base.is_dir():
        raise ValueError(f"Not a directory: {src_dir}")

    # normalize versions once
    allowed_versions_norm: set[str] | None = None
    if allowed_versions:
        allowed_versions_norm = {
            v.strip().lower() for v in allowed_versions if v.strip()
        }

    # Nur Top-Level-Pakete (direkte Unterordner von src_dir)
    for pkg in sorted(p for p in base.iterdir() if p.is_dir()):
        # NEW: Versionsfilter auf Paketebene
        if allowed_versions_norm is not None:
            v = _pkg_version(pkg.name)
            if v is None or v.lower() not in allowed_versions_norm:
                continue

        # bisheriger Filter: Paketname muss Tokens mit allowed_words teilen
        pkg_tokens = _tokenize(pkg.name)
        if not (pkg_tokens & allowed_tokens):
            continue

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
            "No HMM files found for selected packages in '%s' with hmm_sets=%s and packages=%s",
            src_dir,
            allowed_words,
            sorted(allowed_versions_norm)
            if allowed_versions_norm is not None
            else None,
        )
        return

    os.makedirs(os.path.dirname(output_library) or ".", exist_ok=True)
    with open(output_library, "w") as out:
        for fp in files_to_concat:
            with open(fp, "r") as fin:
                shutil.copyfileobj(fin, out)
    log.info("Concatenated %d HMMs into %s", len(files_to_concat), output_library)


def ressource_preparation(config) -> None:
    """Prepare the result space and all required resources.

    Steps performed:
      - Create `reports` and `cross_check` subdirectories under results.
      - Build the HMM library and cutoff files, restricted to selected sets
        if specified.
      - Concatenate default pattern and co-occurrence files if missing.
      - Normalize pattern/co-occurrence files via `format_pattern_files_inplace`.
      - Verify presence of HMM library, score thresholds, patterns,
        co-occurrence, exclusion singletons, and reference sequences.
      - Update `config` with paths to cross-check and CSB output files.

    Args:
        config: Configuration object with attributes for paths and resources.

    Side Effects:
        Creates directories and writes combined resource files as needed.

    Raises:
        SystemExit: If required resources are missing after preparation.
    """
    print_header("Preparing result space and resources", logger=log)

    # Ergebnis-Unterordner
    reports_dir = _ensure_dir(Path(config.result_files_directory) / "reports")
    cross_dir = _ensure_dir(Path(config.result_files_directory) / "cross_check")
    config.paths.refseq = _ensure_dir(config.paths.refseq)

    # ---- HMM-Sets (optional eingeschränkt) ----
    if config.hmm_sets:
        log.info(f"For the library collecting HMMs with token {config.hmm_sets}")
        packages = (
            config.hmm_packages
            if isinstance(config.hmm_packages, list)
            else config.hmm_packages.split()
        )

        allowed = (
            config.hmm_sets
            if isinstance(config.hmm_sets, list)
            else config.hmm_sets.split()
        )

        # baut aus DATA_ROOT/<grp>/*.hmm eine Library
        concatenate_hmms_from_selected_metabolism_packages(
            str(DATA_DIR), allowed, config.library, packages
        )

    # ---- Library, Cutoffs, Cooccurrence, Patterns, Metabolism ggf. zusammenführen ----
    if not os.path.isfile(config.library):
        queue.concatenate_files_shell(str(DATA_DIR), "grp", ".hmm", config.library)

    if not os.path.isfile(config.score_threshold_file):
        queue.concatenate_files_shell(
            str(DATA_DIR), "cutoffs", ".txt", config.score_threshold_file
        )

    if not os.path.isfile(config.cooccurrence_file):
        queue.concatenate_files_shell(
            str(DATA_DIR), "cooccurrence", ".txt", config.cooccurrence_file
        )
        queue.format_pattern_files_inplace(
            config.cooccurrence_file, "cpb-", "_"
        )  # co-occurring protein blocks

    if not os.path.isfile(config.patterns_file):
        queue.concatenate_files_shell(
            str(DATA_DIR), "patterns", ".txt", config.patterns_file
        )
        queue.format_pattern_files_inplace(
            config.patterns_file, "dsb-", "_"
        )  # defined syntenic blocks

    if not os.path.isfile(config.exclusion_singletons):
        queue.concatenate_files_shell(
            str(DATA_DIR), "exclusion_singletons", ".txt", config.exclusion_singletons
        )

    if not os.path.isfile(config.metabolic_information):
        queue.concatenate_files_shell(
            str(DATA_DIR), "metabolic_information", ".txt", config.metabolic_information
        )
    if not os.path.isfile(config.rm_ram_profile):
        queue.concatenate_files_shell(
            str(GPKG_DIR), "ram_profile", ".txt", config.rm_ram_profile
        )
    # ---- Existenz der Ressourcen sicherstellen ----
    _require_path_exists(config.library, "HMM library")
    _require_path_exists(config.score_threshold_file, "Score thresholds")
    _require_path_exists(config.patterns_file, "Patterns")
    _require_path_exists(config.cooccurrence_file, "Cooccurrence")
    _require_path_exists(config.exclusion_singletons, "Exclusion_singletons")
    _require_path_exists(config.paths.refseq, "Reference sequences")
    # _require_path_exists(config.metabolic_information, "Metabolism information")

    # ---- Ableitungen & Artefakte in options hinterlegen ----
    config.cross_check_directory = cross_dir

    config.csb_output_file = str(
        Path(config.result_files_directory) / f"{config.name}_csb.tsv"
    )

    log.debug("Result dir: %s", config.result_files_directory)
    log.debug("Reports dir: %s", reports_dir)
    log.debug("Cross-check dir: %s", config.cross_check_directory)

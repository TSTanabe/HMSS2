from __future__ import annotations

import os
from pathlib import Path

from hmsss.core.logging import get_logger, print_header
from hmsss.cli.paths import DATA_DIR
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

    # ---- HMM-Sets (optional eingeschränkt) ----
    if config.hmm_sets:
        allowed = (
            config.hmm_sets
            if isinstance(config.hmm_sets, list)
            else config.hmm_sets.split()
        )
        # baut aus DATA_ROOT/<grp>/*.hmm eine Library
        queue.concatenate_selected_hmms(
            str(DATA_DIR), allowed, "", ".hmm", config.library
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

    # ---- Existenz der Ressourcen sicherstellen ----
    _require_path_exists(config.library, "HMM library")
    _require_path_exists(config.score_threshold_file, "Score thresholds")
    _require_path_exists(config.patterns_file, "Patterns")
    _require_path_exists(config.cooccurrence_file, "Cooccurrence")
    _require_path_exists(config.exclusion_singletons, "Exclusion_singletons")
    _require_path_exists(config.paths.refseq, "Reference sequences")
    _require_path_exists(config.metabolic_information, "Metabolism information")

    # ---- Ableitungen & Artefakte in options hinterlegen ----
    config.cross_check_directory = cross_dir

    config.csb_output_file = str(
        Path(config.result_files_directory) / f"{config.name}_csb.tsv"
    )

    log.debug("Result dir: %s", config.result_files_directory)
    log.debug("Reports dir: %s", reports_dir)
    log.debug("Cross-check dir: %s", config.cross_check_directory)

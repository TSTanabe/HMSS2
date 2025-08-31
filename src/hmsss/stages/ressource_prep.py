from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.cli.paths import DATA_DIR

# IO-Helfer
from hmsss.io import queue as queue

log = get_logger(__name__)


def _ensure_dir(p: str | os.PathLike) -> str:
    Path(p).mkdir(parents=True, exist_ok=True)
    return str(p)


def _require_path_exists(p: str, desc: str) -> None:
    if not Path(p).exists():
        log.error("%s not found: %s", desc, p)
        raise SystemExit(f"Missing required resource: {desc} -> {p}")


def ressource_preparation(config) -> None:
    """
    Aus __main__.py extrahiert:
    - HMM-Library & Cutoffs zusammenbauen (ggf. nur selektierte Sets)
    - Patterns / Cooccurrence normalisieren
    - Cross-check / Reports-Ordner anlegen
    - Abgeleitete Pfade am Options-Objekt setzen
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

    # ---- Library, Cutoffs, Cooccurrence, Patterns ggf. zusammenführen ----
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
            config.cooccurrence_file, "cpb-", config.csb_name_suffix
        )  # co-occurring protein blocks

    if not os.path.isfile(config.patterns_file):
        queue.concatenate_files_shell(
            str(DATA_DIR), "patterns", ".txt", config.patterns_file
        )
        queue.format_pattern_files_inplace(
            config.patterns_file, "dsb-", config.csb_name_suffix
        )  # defined syntenic blocks

    # ---- Existenz der Ressourcen sicherstellen ----
    _require_path_exists(config.library, "HMM library")
    _require_path_exists(config.score_threshold_file, "Score thresholds")
    _require_path_exists(config.patterns_file, "Patterns")
    _require_path_exists(config.cooccurrence_file, "Cooccurrence")
    _require_path_exists(config.exclusion_singletons, "Exclusion_singletons")
    _require_path_exists(config.paths.refseq, "Reference sequences")

    # ---- Ableitungen & Artefakte in options hinterlegen ----
    config.cross_check_directory = cross_dir

    config.csb_output_file = str(
        Path(config.result_files_directory) / f"{config.name}_csb.tsv"
    )

    log.debug("Result dir: %s", config.result_files_directory)
    log.debug("Reports dir: %s", reports_dir)
    log.debug("Cross-check dir: %s", config.cross_check_directory)

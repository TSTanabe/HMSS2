from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.utils.paths import DATA_ROOT

# IO-Helfer
from hmsss.io import queue as queue

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def _ensure_dir(p: str | os.PathLike) -> str:
    Path(p).mkdir(parents=True, exist_ok=True)
    return str(p)


def _require_path_exists(p: str, desc: str) -> None:
    if not Path(p).exists():
        log.error("%s not found: %s", desc, p)
        raise SystemExit(f"Missing required resource: {desc} -> {p}")


def ressource_preparation(options: Hmsss) -> None:
    """
    Aus __main__.py extrahiert:
    - HMM-Library & Cutoffs zusammenbauen (ggf. nur selektierte Sets)
    - Patterns / Cooccurrence normalisieren
    - Cross-check / Reports-Ordner anlegen
    - Abgeleitete Pfade am Options-Objekt setzen
    """
    print_header("Preparing result space and resources", logger=log)

    # Ergebnis-Unterordner
    reports_dir = _ensure_dir(Path(options.result_files_directory) / "reports")
    cross_dir = _ensure_dir(Path(options.result_files_directory) / "cross_check")

    # ---- HMM-Sets (optional eingeschränkt) ----
    if options.HMM_sets:
        allowed = (
            options.HMM_sets
            if isinstance(options.HMM_sets, list)
            else options.HMM_sets.split()
        )
        # baut aus DATA_ROOT/<grp>/*.hmm eine Library
        queue.concatenate_selected_hmms(
            str(DATA_ROOT), allowed, "", ".hmm", options.library
        )

    # ---- Library, Cutoffs, Cooccurrence, Patterns ggf. zusammenführen ----
    if not os.path.isfile(options.library):
        queue.concatenate_files_shell(str(DATA_ROOT), "grp", ".hmm", options.library)

    if not os.path.isfile(options.score_threshold_file):
        queue.concatenate_files_shell(
            str(DATA_ROOT), "cutoffs", ".txt", options.score_threshold_file
        )

    if not os.path.isfile(options.cooccurrence_file):
        queue.concatenate_files_shell(
            str(DATA_ROOT), "cooccurrence", ".txt", options.cooccurrence_file
        )
        queue.format_pattern_files_inplace(
            options.cooccurrence_file, "cpb-", options.csb_name_suffix
        )  # co-occurring protein blocks

    if not os.path.isfile(options.patterns_file):
        queue.concatenate_files_shell(
            str(DATA_ROOT), "patterns", ".txt", options.patterns_file
        )
        queue.format_pattern_files_inplace(
            options.patterns_file, "dsb-", options.csb_name_suffix
        )  # defined syntenic blocks

    # ---- Existenz der Ressourcen sicherstellen ----
    _require_path_exists(options.library, "HMM library")
    _require_path_exists(options.score_threshold_file, "Score thresholds")
    _require_path_exists(options.patterns_file, "Patterns")
    _require_path_exists(options.cooccurrence_file, "Cooccurrence")
    _require_path_exists(options.exclusion_singletons, "Exclusion_singletons")
    _require_path_exists(options.reference_seq_dir, "Reference sequences")

    # ---- Ableitungen & Artefakte in options hinterlegen ----
    options.Cross_check_directory = cross_dir
    options.glob_trusted_hitreport = str(Path(reports_dir) / "trusted.hmmreport")
    options.glob_intermediate_hitreport = str(
        Path(reports_dir) / "intermediate.hmmreport"
    )
    options.csb_output_file = str(
        Path(options.result_files_directory) / f"{options.name}_csb.tsv"
    )

    log.debug("Result dir: %s", options.result_files_directory)
    log.debug("Reports dir: %s", reports_dir)
    log.debug("Cross-check dir: %s", options.Cross_check_directory)

# HMSS2/src/hmsss/cli/paths.py
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict

"""
Project path detection and constants for HMSS2/HMSSS.

This module infers the project root (both for source and frozen builds),
derives standardized subdirectories (bin, data, results, package),
and exposes helpers to refresh/rewrite paths for tests.

The constants are initialized on import; `refresh_paths()` can be used to
recompute them for a different base (e.g., temporary test roots).
"""

__all__ = [
    "ROOT_DIR",
    "BIN_DIR",
    "DATA_DIR",
    "GPKG_DIR",
    "HMMS_DIR",
    "PLACEMENT_DIR",
    "REFSEQ_DIR",
    "RESULTS_DIR",
    "PACKAGE_DIR",
    "SRC_FILE_RESOURCE_METADATA",
    "SRC_FILE_HMM_LIBRARY",
    "SRC_FILE_COOCCURRENCE",
    "SRC_FILE_THRESHOLDS",
    "SRC_FILE_PATTERNS",
    "SRC_FILE_PRECOMPUTED_PATHWAYS",
    "SRC_FILE_PLAUSIBILITY_MODELS",
    "SRC_FILE_EXCLUSION_SINGLETONS",
    "SRC_FILE_METABOLISM_INFORMATION",
    "SRC_FILE_GPKG_RAM_INFORMATION",
    "refresh_paths",
    "as_dict",
]


def _detect_root_from_source(here: Path) -> Path:
    """Ascend from `here` to the project root in source layout.

    Expected layout: HMSS2/src/hmsss/cli/paths.py

    Args:
        here: Absolute path to this file.

    Returns:
        Path to the detected project root.
    """
    # paths.py -> cli -> hmsss -> src -> HMSS2
    try:
        return here.parents[3]  # 0:cli, 1:hmsss, 2:src, 3:HMSS2
    except IndexError:
        # Fallback: find the directory that is the location of the src/hmsss directory
        for candidate in here.parents:
            if (candidate / "src" / "hmsss").is_dir():
                return candidate
        return here.parent  # letzter Fallback


def _detect_root() -> Path:
    """Detect the project root for both frozen and source distributions.

    - Frozen builds: use the directory of the executable.
    - Source builds: resolve from `__file__` via `_detect_root_from_source()`.

    Returns:
        Absolute `Path` to the project root.
    """
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    here = Path(__file__).resolve()
    return _detect_root_from_source(here)


def _make_paths(root: Path) -> Dict[str, Path]:
    """Construct the full set of canonical project paths from `root`.

    Args:
        root: Project root directory.

    Returns:
        Mapping of path labels to absolute `Path` objects, including DATA_DIR,
        HMMS_DIR, REFSEQ_DIR, RESULTS_DIR, PACKAGE_DIR, and source resource files.
    """
    data = root / "data"

    return {
        "ROOT_DIR": root,
        "BIN_DIR": root / "bin",
        "DATA_DIR": data,
        "GPKG_DIR": data / "gpkg",
        "PLACEMENT_DIR": data / "placements",
        # Consistent with parse.py-Defaults: HMMlib
        "HMMS_DIR": data / "HMMs",
        "REFSEQ_DIR": data / "RefSeqs",
        "RESULTS_DIR": root / "results",
        "PACKAGE_DIR": root / "src" / "hmsss",
        # Standard resources
        "SRC_FILE_RESOURCE_METADATA": data / "Resource_module_metadata",
        "SRC_FILE_HMM_LIBRARY": data / "HMMlib",
        "SRC_FILE_COOCCURRENCE": data / "Cooccurrence",
        "SRC_FILE_THRESHOLDS": data / "Thresholds",
        "SRC_FILE_PATTERNS": data / "Patterns",
        "SRC_FILE_PRECOMPUTED_PATHWAYS": data / "Precomputed_pathways",
        "SRC_FILE_PLAUSIBILITY_MODELS": data / "Plausibility_models.jsonl",
        "SRC_FILE_EXCLUSION_SINGLETONS": data / "Exclusion_singletons",
        "SRC_FILE_METABOLISM_INFORMATION": data / "Metabolism_information",
        # gpkg resources
        "SRC_FILE_GPKG_RAM_INFORMATION": data / "gpkg" / "RAM_profile.txt",
    }


# Legacy module-level constants are retained for existing import sites.
# All values originate from _make_paths(), avoiding duplicate assignments.
_paths: Dict[str, Path] = {}


def refresh_paths(base: str | Path | None = None) -> None:
    """Refresh canonical paths (including all exported SRC_FILE_* paths)."""
    root = Path(base).resolve() if base is not None else _detect_root()
    values = _make_paths(root)
    _paths.clear()
    _paths.update(values)
    globals().update(values)


def as_dict(
    str_paths: bool = True, include_sources: bool = False
) -> Dict[str, str | Path]:
    """Return canonical paths. Defaults to the non-resource directory keys."""
    values = {
        key: value
        for key, value in _paths.items()
        if include_sources or not key.startswith("SRC_FILE_")
    }
    if str_paths:
        return {key: str(value) for key, value in values.items()}
    return values


refresh_paths()

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
    "HMMS_DIR",
    "REFSEQ_DIR",
    "RESULTS_DIR",
    "PACKAGE_DIR",
    "SRC_FILE_HMM_LIBRARY",
    "SRC_FILE_COOCCURRENCE",
    "SRC_FILE_THRESHOLDS",
    "SRC_FILE_PATTERNS",
    "SRC_FILE_EXCLUSION_SINGLETONS",
    "SRC_FILE_METABOLISM_INFORMATION",
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
        # Consistent with parse.py-Defaults: HMMlib
        "HMMS_DIR": data / "HMMs",
        "REFSEQ_DIR": data / "RefSeqs",
        "RESULTS_DIR": root / "results",
        "PACKAGE_DIR": root / "src" / "hmsss",
        # Standard resources
        "SRC_FILE_HMM_LIBRARY": data / "HMMlib",
        "SRC_FILE_COOCCURRENCE": data / "Cooccurrence",
        "SRC_FILE_THRESHOLDS": data / "Thresholds",
        "SRC_FILE_PATTERNS": data / "Patterns",
        "SRC_FILE_EXCLUSION_SINGLETONS": data / "Exclusion_singletons",
        "SRC_FILE_METABOLISM_INFORMATION": data / "Metabolism_information",
    }


# Initialisierung bei Modulimport
_paths = _make_paths(_detect_root())

# Public constants
ROOT_DIR: Path = _paths["ROOT_DIR"]
BIN_DIR: Path = _paths["BIN_DIR"]
DATA_DIR: Path = _paths["DATA_DIR"]
HMMS_DIR: Path = _paths["HMMS_DIR"]
REFSEQ_DIR: Path = _paths["REFSEQ_DIR"]
RESULTS_DIR: Path = _paths["RESULTS_DIR"]
PACKAGE_DIR: Path = _paths["PACKAGE_DIR"]

SRC_FILE_HMM_LIBRARY: Path = _paths["SRC_FILE_HMM_LIBRARY"]
SRC_FILE_COOCCURRENCE: Path = _paths["SRC_FILE_COOCCURRENCE"]
SRC_FILE_THRESHOLDS: Path = _paths["SRC_FILE_THRESHOLDS"]
SRC_FILE_PATTERNS: Path = _paths["SRC_FILE_PATTERNS"]
SRC_FILE_EXCLUSION_SINGLETONS: Path = _paths["SRC_FILE_EXCLUSION_SINGLETONS"]
SRC_FILE_METABOLISM_INFORMATION: Path = _paths["SRC_FILE_METABOLISM_INFORMATION"]

def refresh_paths(base: str | Path | None = None) -> None:
    """Reinitialize exported path constants (useful for tests).

    Args:
        base: If provided, force the project root to this directory; otherwise
              auto-detect as in `_detect_root()`.

    Side Effects:
        Updates module-level constants such as ROOT_DIR, DATA_DIR, RESULTS_DIR, etc.
    """
    global \
        ROOT_DIR, \
        BIN_DIR, \
        DATA_DIR, \
        HMMS_DIR, \
        REFSEQ_DIR, \
        RESULTS_DIR, \
        PACKAGE_DIR, \
        SRC_FILE_HMM_LIBRARY, \
        SRC_FILE_COOCCURRENCE, \
        SRC_FILE_THRESHOLDS, \
        SRC_FILE_PATTERNS, \
        SRC_FILE_EXCLUSION_SINGLETONS, \
        _paths

    root = Path(base).resolve() if base else _detect_root()
    _paths = _make_paths(root)

    ROOT_DIR = _paths["ROOT_DIR"]
    BIN_DIR = _paths["BIN_DIR"]
    DATA_DIR = _paths["DATA_DIR"]
    HMMS_DIR = _paths["HMMS_DIR"]
    REFSEQ_DIR = _paths["REFSEQ_DIR"]
    RESULTS_DIR = _paths["RESULTS_DIR"]
    PACKAGE_DIR = _paths["PACKAGE_DIR"]

    SRC_FILE_HMM_LIBRARY = _paths["SRC_FILE_HMM_LIBRARY"]
    SRC_FILE_COOCCURRENCE = _paths["SRC_FILE_COOCCURRENCE"]
    SRC_FILE_THRESHOLDS = _paths["SRC_FILE_THRESHOLDS"]
    SRC_FILE_PATTERNS = _paths["SRC_FILE_PATTERNS"]
    SRC_FILE_EXCLUSION_SINGLETONS = _paths["SRC_FILE_EXCLUSION_SINGLETONS"]
    SRC_FILE_METABOLISM_INFORMATION = _paths["SRC_FILE_METABOLISM_INFORMATION"]


def as_dict(
    str_paths: bool = True, include_sources: bool = False
) -> Dict[str, str | Path]:
    """Return all core paths as a dictionary (for logs/debugging).

    Args:
        str_paths: If True, convert values to strings.
        include_sources: If True, also include the SRC_FILE_* entries.

    Returns:
        Dictionary of path names to values.
    """
    d: Dict[str, Path] = {
        "ROOT_DIR": ROOT_DIR,
        "BIN_DIR": BIN_DIR,
        "DATA_DIR": DATA_DIR,
        "HMMS_DIR": HMMS_DIR,
        "REFSEQ_DIR": REFSEQ_DIR,
        "RESULTS_DIR": RESULTS_DIR,
        "PACKAGE_DIR": PACKAGE_DIR,
    }
    if include_sources:
        d.update(
            {
                "SRC_FILE_HMM_LIBRARY": SRC_FILE_HMM_LIBRARY,
                "SRC_FILE_COOCCURRENCE": SRC_FILE_COOCCURRENCE,
                "SRC_FILE_THRESHOLDS": SRC_FILE_THRESHOLDS,
                "SRC_FILE_PATTERNS": SRC_FILE_PATTERNS,
                "SRC_FILE_EXCLUSION_SINGLETONS": SRC_FILE_EXCLUSION_SINGLETONS,
                "SRC_FILE_METABOLISM_INFORMATION": SRC_FILE_METABOLISM_INFORMATION
            }
        )
    if str_paths:
        return {k: str(v) for k, v in d.items()}
    return d

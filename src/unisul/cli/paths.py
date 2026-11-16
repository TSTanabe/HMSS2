from __future__ import annotations

import sys
from pathlib import Path

"""
Project path detection and constants for UNISUL.

This module detects the project root for source and frozen builds, derives
canonical project/resource paths and exposes helpers to refresh these paths,
for example when using temporary roots during tests.
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
    """
    Detect the project root from the source layout.

    Expected layout:
        <project-root>/src/unisul/cli/paths.py
    """
    try:
        return here.parents[3]
    except IndexError:
        for candidate in here.parents:
            if (candidate / "src" / "unisul").is_dir():
                return candidate

        return here.parent


def _detect_root() -> Path:
    """
    Detect the UNISUL project root.

    Frozen build:
        directory containing the executable

    Source build:
        project root inferred from this module location
    """
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent

    return _detect_root_from_source(Path(__file__).resolve())


def _make_paths(root: Path) -> dict[str, Path]:
    """Construct all canonical UNISUL paths from a project root."""
    data = root / "data"

    return {
        "ROOT_DIR": root,
        "BIN_DIR": root / "bin",
        "DATA_DIR": data,
        "GPKG_DIR": data / "gpkg",
        "HMMS_DIR": data / "HMMs",
        "PLACEMENT_DIR": data / "placements",
        "REFSEQ_DIR": data / "RefSeqs",
        "RESULTS_DIR": Path.cwd().resolve() / "results",
        "PACKAGE_DIR": root / "src" / "unisul",
        "SRC_FILE_RESOURCE_METADATA": data / "Resource_module_metadata",
        "SRC_FILE_HMM_LIBRARY": data / "HMMlib",
        "SRC_FILE_COOCCURRENCE": data / "Cooccurrence",
        "SRC_FILE_THRESHOLDS": data / "Thresholds",
        "SRC_FILE_PATTERNS": data / "Patterns",
        "SRC_FILE_PRECOMPUTED_PATHWAYS": data / "Precomputed_pathways",
        "SRC_FILE_PLAUSIBILITY_MODELS": data / "Plausibility_models.jsonl",
        "SRC_FILE_EXCLUSION_SINGLETONS": data / "Exclusion_singletons",
        "SRC_FILE_METABOLISM_INFORMATION": data / "Metabolism_information",
        "SRC_FILE_GPKG_RAM_INFORMATION": data / "gpkg" / "RAM_profile.txt",
    }


# ---------------------------------------------------------------------------
# Module-level path variables
# ---------------------------------------------------------------------------

ROOT_DIR: Path
BIN_DIR: Path
DATA_DIR: Path
GPKG_DIR: Path
HMMS_DIR: Path
PLACEMENT_DIR: Path
REFSEQ_DIR: Path
RESULTS_DIR: Path
PACKAGE_DIR: Path

SRC_FILE_RESOURCE_METADATA: Path
SRC_FILE_HMM_LIBRARY: Path
SRC_FILE_COOCCURRENCE: Path
SRC_FILE_THRESHOLDS: Path
SRC_FILE_PATTERNS: Path
SRC_FILE_PRECOMPUTED_PATHWAYS: Path
SRC_FILE_PLAUSIBILITY_MODELS: Path
SRC_FILE_EXCLUSION_SINGLETONS: Path
SRC_FILE_METABOLISM_INFORMATION: Path
SRC_FILE_GPKG_RAM_INFORMATION: Path


_paths: dict[str, Path] = {}


def refresh_paths(base: str | Path | None = None) -> None:
    """
    Recompute all canonical paths.

    Args:
        base:
            Optional project root. If omitted, the root is detected
            automatically from the executable or source layout.
    """
    global ROOT_DIR
    global BIN_DIR
    global DATA_DIR
    global GPKG_DIR
    global HMMS_DIR
    global PLACEMENT_DIR
    global REFSEQ_DIR
    global RESULTS_DIR
    global PACKAGE_DIR
    global SRC_FILE_RESOURCE_METADATA
    global SRC_FILE_HMM_LIBRARY
    global SRC_FILE_COOCCURRENCE
    global SRC_FILE_THRESHOLDS
    global SRC_FILE_PATTERNS
    global SRC_FILE_PRECOMPUTED_PATHWAYS
    global SRC_FILE_PLAUSIBILITY_MODELS
    global SRC_FILE_EXCLUSION_SINGLETONS
    global SRC_FILE_METABOLISM_INFORMATION
    global SRC_FILE_GPKG_RAM_INFORMATION

    root = Path(base).resolve() if base is not None else _detect_root()
    values = _make_paths(root)

    _paths.clear()
    _paths.update(values)

    ROOT_DIR = values["ROOT_DIR"]
    BIN_DIR = values["BIN_DIR"]
    DATA_DIR = values["DATA_DIR"]
    GPKG_DIR = values["GPKG_DIR"]
    HMMS_DIR = values["HMMS_DIR"]
    PLACEMENT_DIR = values["PLACEMENT_DIR"]
    REFSEQ_DIR = values["REFSEQ_DIR"]
    RESULTS_DIR = values["RESULTS_DIR"]
    PACKAGE_DIR = values["PACKAGE_DIR"]

    SRC_FILE_RESOURCE_METADATA = values["SRC_FILE_RESOURCE_METADATA"]
    SRC_FILE_HMM_LIBRARY = values["SRC_FILE_HMM_LIBRARY"]
    SRC_FILE_COOCCURRENCE = values["SRC_FILE_COOCCURRENCE"]
    SRC_FILE_THRESHOLDS = values["SRC_FILE_THRESHOLDS"]
    SRC_FILE_PATTERNS = values["SRC_FILE_PATTERNS"]
    SRC_FILE_PRECOMPUTED_PATHWAYS = values["SRC_FILE_PRECOMPUTED_PATHWAYS"]
    SRC_FILE_PLAUSIBILITY_MODELS = values["SRC_FILE_PLAUSIBILITY_MODELS"]
    SRC_FILE_EXCLUSION_SINGLETONS = values["SRC_FILE_EXCLUSION_SINGLETONS"]
    SRC_FILE_METABOLISM_INFORMATION = values["SRC_FILE_METABOLISM_INFORMATION"]
    SRC_FILE_GPKG_RAM_INFORMATION = values["SRC_FILE_GPKG_RAM_INFORMATION"]


def as_dict(
    str_paths: bool = True,
    include_sources: bool = False,
) -> dict[str, str | Path]:
    """
    Return the currently configured paths.

    Args:
        str_paths:
            Convert Path objects to strings.

        include_sources:
            Include SRC_FILE_* resource paths. By default only directory-like
            project paths are returned.
    """
    values = {
        key: value
        for key, value in _paths.items()
        if include_sources or not key.startswith("SRC_FILE_")
    }

    if str_paths:
        return {key: str(value) for key, value in values.items()}

    return dict(values)


# Initialize paths on import.
refresh_paths()

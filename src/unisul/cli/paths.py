from __future__ import annotations

import os
import sys
from pathlib import Path

from unisul.data_management import resolve_data_root

ROOT_DIR: Path
BIN_DIR: Path
DATA_ROOT: Path
DATA_DIR: Path
GPKG_DIR: Path
HMMS_DIR: Path
PLACEMENT_DIR: Path
REFSEQ_DIR: Path
RESULTS_DIR: Path
PACKAGE_DIR: Path
TAXONOMY_FILE: Path
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

__all__ = [
    "ROOT_DIR",
    "BIN_DIR",
    "DATA_ROOT",
    "DATA_DIR",
    "GPKG_DIR",
    "HMMS_DIR",
    "PLACEMENT_DIR",
    "REFSEQ_DIR",
    "RESULTS_DIR",
    "PACKAGE_DIR",
    "TAXONOMY_FILE",
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


def _package_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def _detect_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    package_dir = _package_dir()
    if package_dir.parent.name == "src":
        return package_dir.parent.parent
    return package_dir.parent


def _make_paths(
    root: Path, data_root: Path, *, legacy_flat_data: bool = False
) -> dict[str, Path]:
    """Build canonical paths.

    ``data_root`` always points directly to the complete external UNISUL data
    directory. ``legacy_flat_data`` is kept only for compatibility with tests
    or callers using ``refresh_paths(base=...)``.
    """
    data = data_root
    return {
        "ROOT_DIR": root,
        "BIN_DIR": Path(
            os.environ.get("UNISUL_BIN_DIR")
            or os.environ.get("UNISUL_BIN")
            or root / "bin"
        ).expanduser(),
        "DATA_ROOT": data_root,
        "DATA_DIR": data,
        "GPKG_DIR": data / "gpkg",
        "HMMS_DIR": data / "HMMs",
        "PLACEMENT_DIR": data / "placements",
        "REFSEQ_DIR": data / "RefSeqs",
        "RESULTS_DIR": Path.cwd().resolve() / "results",
        "PACKAGE_DIR": (
            root / "src" / "unisul"
            if legacy_flat_data or getattr(sys, "frozen", False)
            else _package_dir()
        ),
        "TAXONOMY_FILE": data / "RefSeqs" / "combined_globdb_taxonomy.tsv",
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


_paths: dict[str, Path] = {}


def refresh_paths(
    base: str | Path | None = None,
    *,
    data_dir: str | Path | None = None,
    preserve_root: bool = False,
) -> None:
    """Refresh path constants after CLI options are parsed.

    ``base`` is retained for compatibility with callers and tests that model
    the previous repository layout. Runtime data selection uses ``data_dir``.
    """
    if base is not None and data_dir is not None:
        raise ValueError("Specify either base or data_dir, not both")
    root = (
        Path(base).resolve()
        if base is not None
        else _paths["ROOT_DIR"].resolve()
        if preserve_root and "ROOT_DIR" in _paths
        else _detect_root()
    )
    if base is not None:
        selected_data = root / "data"
    else:
        selected_data = resolve_data_root(data_dir)
    values = _make_paths(root, selected_data, legacy_flat_data=base is not None)

    _paths.clear()
    _paths.update(values)
    globals().update(values)


def as_dict(
    str_paths: bool = True,
    include_sources: bool = False,
) -> dict[str, str | Path]:
    values = {
        key: value
        for key, value in _paths.items()
        if key not in {"DATA_ROOT", "TAXONOMY_FILE"}
        if include_sources or not key.startswith("SRC_FILE_")
    }
    if str_paths:
        return {key: str(value) for key, value in values.items()}
    return dict(values)


refresh_paths()

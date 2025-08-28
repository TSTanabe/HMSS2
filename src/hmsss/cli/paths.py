# HMSS2/src/hmsss/cli/paths.py
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict

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
    "refresh_paths",
    "as_dict",
]


def _detect_root_from_source(here: Path) -> Path:
    """
    Läuft vom aktuellen Dateipfad (…/src/hmsss/cli/paths.py) nach oben zum Projektwurzelordner HMSS2.
    Erwartete Struktur: HMSS2/src/hmsss/cli/paths.py
    """
    # paths.py -> cli -> hmsss -> src -> HMSS2
    try:
        return here.parents[3]  # 0:cli, 1:hmsss, 2:src, 3:HMSS2
    except IndexError:
        # Fallback: Suche aufwärts nach einem Verzeichnis, das src/hmsss enthält
        for candidate in here.parents:
            if (candidate / "src" / "hmsss").is_dir():
                return candidate
        return here.parent  # letzter Fallback


def _detect_root() -> Path:
    """
    Ermittelt den Projektwurzelordner HMSS2.

    - Kompilierte Version: sys.frozen == True → Executable liegt direkt unter HMSS2/
      → Wurzel = Parent des Executables.
    - Source-Version: Ableitung über __file__ → .../HMSS2/src/hmsss/cli/paths.py
    """
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    here = Path(__file__).resolve()
    return _detect_root_from_source(here)


def _make_paths(root: Path) -> Dict[str, Path]:
    """
    Baut das Pfad-Set relativ zur festen Projektstruktur.
    """
    data = root / "data"
    return {
        "ROOT_DIR": root,
        "BIN_DIR": root / "bin",
        "DATA_DIR": data,
        # Konsistent zu parse.py-Defaults: HMMlib
        "HMMS_DIR": data / "HMMs",
        "REFSEQ_DIR": data / "RefSeqs",
        "RESULTS_DIR": root / "results",
        "PACKAGE_DIR": root / "src" / "hmsss",
        # Quell-Standardressourcen
        "SRC_FILE_HMM_LIBRARY": data / "HMMlib",
        "SRC_FILE_COOCCURRENCE": data / "Cooccurrence",
        "SRC_FILE_THRESHOLDS": data / "Thresholds",
        "SRC_FILE_PATTERNS": data / "Patterns",
        "SRC_FILE_EXCLUSION_SINGLETONS": data / "Exclusion_singletons",
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


def refresh_paths(base: str | Path | None = None) -> None:
    """
    Reinitialisiert die Pfade (z. B. für Tests).
    - base=None  → Autodetektion (wie oben)
    - base=Pfad  → Erzwingt HMSS2-Wurzel = base
    """
    global ROOT_DIR, BIN_DIR, DATA_DIR, HMMS_DIR, REFSEQ_DIR, RESULTS_DIR, PACKAGE_DIR, \
           SRC_FILE_HMM_LIBRARY, SRC_FILE_COOCCURRENCE, SRC_FILE_THRESHOLDS, \
           SRC_FILE_PATTERNS, SRC_FILE_EXCLUSION_SINGLETONS, _paths

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


def as_dict(str_paths: bool = True, include_sources: bool = False) -> Dict[str, str | Path]:
    """
    Gibt alle Pfade als Dict zurück (für Debug/Logs).
    - str_paths=True → Werte als Strings
    - include_sources=True → zusätzlich die SRC_FILE_* Konstanten aufnehmen
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
            }
        )
    if str_paths:
        return {k: str(v) for k, v in d.items()}
    return d

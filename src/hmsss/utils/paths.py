# src/hmsss/utils/paths.py
from __future__ import annotations
import os
from importlib.resources import files
from pathlib import Path

# Paket-Wurzel (…/src/hmsss)
PACKAGE_ROOT: Path = Path(__file__).resolve().parents[1]
# Projekt-Wurzel (…/HMSS2)
PROJECT_ROOT: Path = PACKAGE_ROOT.parent.parent
# Daten-Wurzel: ENV > Fallback auf HMSS2/data
DATA_ROOT: Path = Path(os.environ.get("HMSSS_DATA_DIR", PROJECT_ROOT / "data"))
# Binary-Wurzel: ENV > Fallback auf HMSS2/bin
BIN_ROOT: Path = Path(os.environ.get("HMSSS_BIN_DIR", PROJECT_ROOT / "bin"))


def res_path(rel: str) -> str:
    """Pfad zu paketierten Ressourcen (falls ihr später welche mitliefert)."""
    return str(files("hmsss.resources") / rel)


def data_path(*parts: str) -> str:
    """Pfad innerhalb des externen Datenordners HMSS2/data (oder ENV-Override)."""
    return str(DATA_ROOT.joinpath(*parts))


__all__ = [
    "PACKAGE_ROOT",
    "PROJECT_ROOT",
    "DATA_ROOT",
    "BIN_ROOT",
    "res_path",
    "data_path",
]

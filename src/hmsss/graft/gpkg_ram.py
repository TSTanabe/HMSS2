from __future__ import annotations

import re
from pathlib import Path
from typing import Dict

from hmsss.cli.paths import SRC_FILE_GPKG_RAM_INFORMATION
from hmsss.core.logging import get_logger

logger = get_logger(__name__)

_RESULT_RE = re.compile(
    r"""
    ^\[RESULT\]                          # prefix
    \s+gpkg=(?P<gpkg>\S+)                # gpkg=XYZ
    .*?peak_rss_gb=(?P<peak>[0-9]+(?:\.[0-9]+)?)  # peak_rss_gb=12.34
    (?:\s+exit_code=(?P<exit>\-?\d+))?   # optional exit_code=0
    """,
    re.VERBOSE,
)


def parse_ram_profiles_file(path: Path) -> Dict[str, float]:
    """
    Parse a ram_profiles file written by the gpkg RAM profiler script.

    Expected lines include:
      [RESULT] gpkg=sqdg.gpkg threads=14 reads=1000 peak_rss_gb=31.84 exit_code=0 wall_s=240.7

    Returns:
      { "sqdg.gpkg": 31.84, ... }
    """
    if not path.exists():
        return {}

    out: Dict[str, float] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("[RESULT]"):
                continue

            m = _RESULT_RE.match(line)
            if not m:
                continue

            gpkg = m.group("gpkg")
            peak = float(m.group("peak"))

            # Optional: nur erfolgreiche Läufe berücksichtigen
            exit_code = m.group("exit")
            if exit_code is not None and exit_code != "0":
                continue

            out[gpkg] = peak

    return out


def collect_gpkg_ram(
    gpkg_packages: Dict[str, str],
) -> Dict[str, float]:
    """
    Collect peak RAM usage per GPKG.

    Returns:
        { gpkg_name : peak_rss_gb }
    """
    profiles_by_id = parse_ram_profiles_file(SRC_FILE_GPKG_RAM_INFORMATION)

    # helper indices from profiler output
    by_basename = profiles_by_id  # "sqdg.gpkg" -> peak
    by_stem = {Path(k).stem: v for k, v in profiles_by_id.items()}  # "sqdg" -> peak

    result: Dict[str, float] = {}

    for gpkg_name, gpkg_dir in gpkg_packages.items():
        p = Path(gpkg_dir)
        basename = p.name  # e.g. "sqdg.gpkg"
        stem = p.stem  # e.g. "sqdg"

        peak = None

        # bevorzugte Reihenfolge:
        # 1) explizit gpkg_name (falls Profiler das so geschrieben hat)
        if gpkg_name in profiles_by_id:
            peak = profiles_by_id[gpkg_name]

        # 2) Verzeichnisname mit .gpkg
        elif basename in by_basename:
            peak = by_basename[basename]

        # 3) Stem (ohne .gpkg)
        elif stem in by_stem:
            peak = by_stem[stem]

        if peak is not None:
            result[gpkg_name] = float(peak)
        else:
            logger.debug(f"No RAM profile found for gpkg {gpkg_name}")

    logger.debug(f"RAM profile per package {result}")
    return result

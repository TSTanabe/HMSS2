from __future__ import annotations

import os
import re
import sys
import subprocess
from pathlib import Path
from typing import Dict

from hmsss.cli.config import Config
from hmsss.cli.paths import GPKG_DIR
from hmsss.core import queue as queue
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def collect_gpkg_from_selected_metabolism_packages(
        src_dir: str,
        allowed_words: list[str],
) -> dict[str, str]:
    """
    Return a FLAT dictionary mapping:

        { gpkg_basename : full_path_to_gpkg }

    Only metabolism packages (directories directly under src_dir) whose
    name shares at least one token with allowed_words are inspected.
    Inside each selected metabolism package, subdirectories whose name
    contains 'Gpkg' (case-insensitive) are searched recursively for *.gpkg.

    Example return:
        {
            "dsr_core": "/.../v8_Gpkg_Dsr_sHdr_Sox/dsr_core.gpkg",
            "sox_cluster": "/.../v8_Gpkg_SulfurOxidation/sox_cluster.gpkg",
        }
    """

    def _tokenize(name: str) -> set[str]:
        return {t.lower() for t in re.findall(r"[A-Za-z0-9]+", name)}

    base = Path(src_dir)
    if not base.is_dir():
        raise ValueError(f"Not a directory: {src_dir}")

    allowed_tokens = set().union(*(_tokenize(w) for w in allowed_words))

    gpkg_map: dict[str, str] = {}

    # iterate top-level metabolism packages
    for pkg in sorted(p for p in base.iterdir() if p.is_dir()):
        pkg_tokens = _tokenize(pkg.name)
        if not (pkg_tokens & allowed_tokens):
            continue

        # find directories containing "Gpkg" anywhere in their name
        gpkg_dirs: list[Path] = []
        for root, dirs, _files in os.walk(pkg):
            for d in dirs:
                if "Gpkg" in d.lower() or "gpkg" in d.lower():
                    gpkg_dirs.append(Path(root) / d)

        # now collect *.gpkg files from all gpkg dirs
        for gdir in gpkg_dirs:
            for path in gdir.rglob("*.gpkg"):
                if path.is_dir():
                    key = path.stem  # filename without ".gpkg"
                    value = str(path)
                    gpkg_map[key] = value

    return gpkg_map


def collect_gpkg_by_filename_tokens(
        src_dir: str,
        allowed_words: list[str],
) -> dict[str, str]:
    """
    Return a FLAT dictionary mapping:

        { gpkg_basename : full_path_to_gpkg }

    All *.gpkg files under src_dir are considered (recursive).
    Only those GPKGs whose *filename* shares at least one token
    with allowed_words are kept.

    Tokenization:
        - case-insensitive
        - alphanumeric tokens (A–Z, a–z, 0–9)

    Example:
        allowed_words = ["dsr", "sox"]
        matches:
            dsr_core.gpkg
            sox_cluster_v2.gpkg
    """

    def _tokenize(name: str) -> set[str]:
        return {t.lower() for t in re.findall(r"[A-Za-z0-9]+", name)}

    base = Path(src_dir)
    if not base.is_dir():
        raise ValueError(f"Not a directory: {src_dir}")

    allowed_tokens: set[str] = set().union(
        *(_tokenize(w) for w in allowed_words)
    )

    gpkg_map: dict[str, str] = {}

    for path in base.rglob("*.gpkg"):
        if not path.is_file():
            continue

        filename_tokens = _tokenize(path.stem)

        if not (filename_tokens & allowed_tokens):
            continue

        key = path.stem
        gpkg_map[key] = str(path)

    return gpkg_map


def prepare_gpkg_packages(config: Config) -> None | dict[str, str]:
    if config.gpkg_sets:
        logger.info(f"Collecting graftM packages for {config.gpkg_sets}")
        allowed = (
            config.gpkg_sets
            if isinstance(config.gpkg_sets, list)
            else config.gpkg_sets.split()
        )
        return collect_gpkg_from_selected_metabolism_packages(
            str(GPKG_DIR), allowed
        )
    else:
        logger.warning("Define a gpkg set for the read mapping")
        sys.exit()


def prepare_gpkg_packs(config: Config) -> None | dict[str, str]:
    if config.gpkg_packs:
        logger.info(f"Collecting graftM packages for {config.gpkg_packs}")
        allowed = (
            config.gpkg_packs
            if isinstance(config.gpkg_packs, list)
            else config.gpkg_packs.split()
        )
        return collect_gpkg_from_selected_metabolism_packages(
            str(GPKG_DIR), allowed
        )
    else:
        logger.warning("Define a gpkg set for the read mapping")
        sys.exit()


def _ensure_dmnd(
        *,
        faa_path: str,
        dmnd_path: str,
        name: str,
        threads: int = 4,
) -> None:
    """
    Ensure a DIAMOND database exists for the given FASTA.
    """
    if os.path.isfile(dmnd_path):
        logger.debug(f"For gpkg {name} database exists exists: {dmnd_path}")
        return

    if not os.path.isfile(faa_path):
        raise FileNotFoundError(f"Missing FASTA for DIAMOND DB: {faa_path}")

    prefix = dmnd_path[:-5] if dmnd_path.endswith(".dmnd") else dmnd_path

    logger.debug(f"Initilizing {name} package: {dmnd_path}")
    cmd = [
        "diamond",
        "makedb",
        "--in",
        faa_path,
        "-d",
        prefix,
        "--threads",
        str(threads),
    ]
    try:
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except Exception as e:
        logger.error(
            f"Initialization of database for gpkg package {name} failed with exception\n{e}"
        )


def initialize_gpkg_packages(
        package_dict: Dict[str, str],
        *,
        threads: int = 4,
) -> None:
    """
    Ensure that each gpkg contains refseq and decoy DIAMOND databases.

    Parameters
    ----------
    package_dict
        Mapping: gpkg_name -> gpkg_path
    threads
        Threads for diamond makedb
    """
    logger.debug("Initializing GraftM package DIAMOND databases")

    for gpkg_name, gpkg_path in package_dict.items():
        logger.debug(f"Initializing package: {gpkg_name}")

        try:
            ref_faa = os.path.join(gpkg_path, "refseq_database.faa")
            ref_dmnd = os.path.join(gpkg_path, "refseq_database.dmnd")
            decoy_faa = os.path.join(gpkg_path, "decoy_database.faa")
            decoy_dmnd = os.path.join(gpkg_path, "decoy_database.dmnd")

            # --- Refseq DB ---
            if not os.path.isfile(ref_faa):
                logger.warning(f"{gpkg_name}.gpkg Missing refseq FASTA {ref_faa}")
            else:
                _ensure_dmnd(
                    faa_path=str(ref_faa),
                    dmnd_path=str(ref_dmnd),
                    name=f"{gpkg_name}:ref",
                    threads=threads,
                )

            # --- Decoy DB ---
            if not os.path.isfile(decoy_faa):
                logger.warning(f"{gpkg_name}.gpkg Missing decoy FASTA {decoy_faa}")
            else:
                _ensure_dmnd(
                    faa_path=str(decoy_faa),
                    dmnd_path=str(decoy_dmnd),
                    name=f"{gpkg_name}:decoy",
                    threads=threads,
                )

        except Exception as e:
            # Catch ALL package-level errors, continue with next package
            logger.exception(
                f"{gpkg_name} Failed to build databases, skipping package: {e}"
            )
            continue

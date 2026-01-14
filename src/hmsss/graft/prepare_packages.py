from __future__ import annotations

import os
import sys
import subprocess
from typing import Dict

from hmsss.cli.config import Config
from hmsss.cli.paths import GPKG_DIR
from hmsss.core import queue as queue
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def prepare_gpkg_packages(config: Config) -> None | dict[str, str]:
    if config.hmm_sets:
        logger.info(f"Collecting graftM packages for {config.hmm_sets}")
        allowed = (
            config.hmm_sets
            if isinstance(config.hmm_sets, list)
            else config.hmm_sets.split()
        )
        return queue.collect_gpkg_from_selected_metabolism_packages(
            str(GPKG_DIR), allowed
        )
    else:
        logger.warning("Define a library set for the read mapping")
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

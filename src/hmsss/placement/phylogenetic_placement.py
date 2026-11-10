from __future__ import annotations

import logging
import shutil
from pathlib import Path

from hmsss.placement import epang_placement

logger = logging.getLogger("hmg.epa")


def run_placement_stage(config) -> None:
    """
    Run phylogenetic placement as an independent pipeline stage.

    Intermediate-hit FASTA files are collected from the cross-check directory.
    config.placement_domains restricts the domains that are submitted.

    Empty placement_domains means: use all available domains.
    """

    if config.disable_placement:
        logger.info("Phylogenetic placement disabled")
        return

    if not shutil.which("epa-ng"):
        logger.warning("EPA-ng executable not found; skipping placement stage")
        return

    intermediate_files = sorted(
        Path(config.cross_check_directory).glob("*.intermediate_hits.faa")
    )

    if not intermediate_files:
        logger.info("No intermediate-hit FASTA files found for phylogenetic placement")
        return

    suffix = ".intermediate_hits.faa"

    intermediate_dict = {}

    for path in intermediate_files:
        name = path.name

        if not name.endswith(suffix):
            continue

        domain = name[: -len(suffix)]

        intermediate_dict[domain] = path

    if not intermediate_dict:
        logger.info("No intermediate domains available for placement")
        return

    selected_domains = set(config.placement_domains or [])

    # Empty set = all domains
    if selected_domains:
        placement_dict = {
            domain: path
            for domain, path in intermediate_dict.items()
            if domain in selected_domains
        }

        missing = selected_domains - set(intermediate_dict)

        if missing:
            logger.warning(
                "Requested placement domains without intermediate hits: %s",
                ", ".join(sorted(missing)),
            )
    else:
        placement_dict = intermediate_dict.copy()

    if not placement_dict:
        logger.info("No domains selected for placement")
        return

    logger.info(
        "Running phylogenetic placement for %d of %d intermediate domains",
        len(placement_dict),
        len(intermediate_dict),
    )

    epang_placement.epa_ng_placement(
        config=config,
        intermediate_files=placement_dict,
    )

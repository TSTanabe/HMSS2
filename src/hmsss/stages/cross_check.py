from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import search_cross_reference as search


log = get_logger(__name__)


def reference_sequence_check(config: Config) -> None:
    """
    Aus __main__.py:
    - Cross-check NC<score<TC mit RefSeq (Diamond)
    - Cross-checked Hits promoten
    - Rest per optimized cutoff promoten
    - Alternative: nur optimized cutoff verwenden

    Args:
        config (Config): 
    """
    print_header("Reference sequence cross-check & cutoff promotion", logger=log)

    if config.bool_cross_check and not os.path.isfile(config.glob_trusted_hitreport):
        log.info("Generating .faa files for intermediate hit lists")
        search.generate_faa_per_hitfile_parallel(
            config, config.cross_check_directory, config.cores
        )

        log.info("Cross-check candidates with reference sequences")
        refseq_unavailable_list = search.cross_check_candidates_with_reference_seqs(
            config
        )

        log.info("Promoting cross-checked hits to trusted list")
        search.promote_crosschecked_hits(config.cross_check_directory, config.cores)

        log.info("Promoting remaining candidates by optimized cutoff")
        search.promote_by_cutoff(
            config,
            config.cross_check_directory,
            config.cores,
            refseq_unavailable_list,
        )

    elif config.optimized_cutoff_cross_check:
        log.info("Using optimized cutoff instead of trusted/cross-check")
        # neuere Signatur akzeptiert Liste; fallback auf String
        try:
            search.promote_by_cutoff(
                config, config.cross_check_directory, config.cores, ["all"]
            )
        except TypeError:
            search.promote_by_cutoff(
                config, config.cross_check_directory, config.cores, ["all"]
            )

    log.info("Summarizing trusted and intermediate hit reports")
    config.glob_trusted_hitreport = search.summarize_trusted_hits(
        config.result_files_directory,
        config.cross_check_directory,
        "global_trusted_hits_summary.hmmreport",
        ".trusted_hits",
    )
    config.glob_intermediate_hitreport = search.summarize_trusted_hits(
        config.result_files_directory,
        config.cross_check_directory,
        "global_intermediate_hits_summary.hmmreport",
        "intermediate_hits",
    )

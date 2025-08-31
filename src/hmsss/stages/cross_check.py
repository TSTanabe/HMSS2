from __future__ import annotations

import os

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import search_cross_reference as search_cross_reference
from hmsss.db.report_db import load_hmmreports_to_sqlite
from hmsss.utils import myUtil

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
        search_cross_reference.generate_faa_per_hitfile_parallel(
            config, config.cross_check_directory, config.cores
        )

        log.info("Cross-check candidates with reference sequences")
        refseq_unavailable_list = (
            search_cross_reference.cross_check_candidates_with_reference_seqs(config)
        )

        log.info("Promoting cross-checked hits to trusted list")
        search_cross_reference.promote_crosschecked_hits(
            config.cross_check_directory, config.cores
        )

        log.info("Promoting remaining candidates by optimized cutoff")
        search_cross_reference.promote_by_cutoff(
            config,
            config.cross_check_directory,
            config.cores,
            refseq_unavailable_list,
        )

    elif config.optimized_cutoff_cross_check:
        log.info("Using optimized cutoff instead of trusted/cross-check")
        # neuere Signatur akzeptiert Liste; fallback auf String
        try:
            search_cross_reference.promote_by_cutoff(
                config, config.cross_check_directory, config.cores, ["all"]
            )
        except TypeError:
            search_cross_reference.promote_by_cutoff(
                config, config.cross_check_directory, config.cores, ["all"]
            )

    log.info("Summarizing trusted and intermediate hit reports")

    trusted_files = myUtil.get_all_files(config.cross_check_directory, ".trusted_hits")
    load_hmmreports_to_sqlite(config.glob_trusted_hitreport, trusted_files)
    intermediate_files = myUtil.get_all_files(
        config.cross_check_directory, ".intermediate_hits"
    )
    load_hmmreports_to_sqlite(config.glob_intermediate_hitreport, intermediate_files)

    return
    config.glob_trusted_hitreport = search_cross_reference.summarize_trusted_hits(
        config.result_files_directory,
        config.cross_check_directory,
        "global_trusted_hits_summary.hmmreport",
        ".trusted_hits",
    )
    config.glob_intermediate_hitreport = search_cross_reference.summarize_trusted_hits(
        config.result_files_directory,
        config.cross_check_directory,
        "global_intermediate_hits_summary.hmmreport",
        "intermediate_hits",
    )

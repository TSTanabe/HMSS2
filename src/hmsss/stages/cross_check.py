from __future__ import annotations

import os
from typing import TYPE_CHECKING

from hmsss.core.logging import get_logger, print_header
from hmsss.algorithms import search_cross_reference as search

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def reference_sequence_check(options: Hmsss) -> None:
    """
    Aus __main__.py:
    - Cross-check NC<score<TC mit RefSeq (Diamond)
    - Cross-checked Hits promoten
    - Rest per optimized cutoff promoten
    - Alternative: nur optimized cutoff verwenden
    """
    print_header("Reference sequence cross-check & cutoff promotion", logger=log)

    if options.bool_cross_check and not os.path.isfile(options.glob_trusted_hitreport):
        log.info("Generating .faa files for intermediate hit lists")
        search.generate_faa_per_hitfile_parallel(
            options, options.cross_check_directory, options.cores
        )

        log.info("Cross-check candidates with reference sequences")
        refseq_unavailable_list = search.cross_check_candidates_with_reference_seqs(
            options
        )

        log.info("Promoting cross-checked hits to trusted list")
        search.promote_crosschecked_hits(options.cross_check_directory, options.cores)

        log.info("Promoting remaining candidates by optimized cutoff")
        search.promote_by_cutoff(
            options,
            options.cross_check_directory,
            options.cores,
            refseq_unavailable_list,
        )

    elif options.optimized_cutoff_cross_check:
        log.info("Using optimized cutoff instead of trusted/cross-check")
        # neuere Signatur akzeptiert Liste; fallback auf String
        try:
            search.promote_by_cutoff(
                options, options.cross_check_directory, options.cores, ["all"]
            )
        except TypeError:
            search.promote_by_cutoff(
                options, options.cross_check_directory, options.cores, ["all"]
            )

    log.info("Summarizing trusted and intermediate hit reports")
    options.glob_trusted_hitreport = search.summarize_trusted_hits(
        options.result_files_directory,
        options.cross_check_directory,
        "global_trusted_hits_summary.hmmreport",
        ".trusted_hits",
    )
    options.glob_intermediate_hitreport = search.summarize_trusted_hits(
        options.result_files_directory,
        options.cross_check_directory,
        "global_intermediate_hits_summary.hmmreport",
        "intermediate_hits",
    )

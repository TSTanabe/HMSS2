from __future__ import annotations

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.cross_check import search_cross_reference as search_cross_reference

log = get_logger(__name__)

"""
Stage: Reference sequence cross-check.

Validates intermediate hits by DIAMOND against reference sequences and
promotes them to trusted. Remaining hits can be promoted by optimized cutoff.
Finally, trusted and intermediate hit reports are summarized and stored.
"""


def reference_sequence_check(config: Config) -> None:
    """Perform reference sequence cross-check and cutoff promotion.

    Workflow:
      - Generate FAA files for intermediate hits.
      - Cross-check with reference sequences using DIAMOND.
      - Promote validated candidates to trusted hits.
      - Promote remaining candidates by optimized cutoff.
      - Alternatively: promote only by optimized cutoff.

    Args:
        config: Pipeline configuration with cutoff settings and directories.

    Side Effects:
        Creates FAA, crosscheck, trusted, and intermediate hit files;
        updates SQLite with summarized reports.
    """
    print_header("Reference sequence cross-check & cutoff promotion", logger=log)

    log.info("Cross-check candidates with reference sequences")
    search_cross_reference.cross_check_candidates_with_reference_seqs(config)

    log.info("Promoting cross-checked hits to trusted list")
    search_cross_reference.promote_crosschecked_hits_to_db(
        crosscheck_dir=config.cross_check_directory, database_path=config.database_directory, processes= config.cores
    )

    return

from __future__ import annotations

from typing import TYPE_CHECKING
from hmsss.core.logging import get_logger, print_header
from hmsss.io import parse_reports as parse_reports

if TYPE_CHECKING:
    from hmsss.core.options import Hmsss

log = get_logger(__name__)


def parse_reports_to_database(options: Hmsss) -> None:
    """
    Aus __main__.py: Summary-HMMreport in die DB schreiben.
    """
    print_header("Parse reports → database", logger=log)
    parse_reports.main_parse_summary_hmmreport(options)

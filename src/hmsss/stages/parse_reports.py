from __future__ import annotations

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger, print_header
from hmsss.io import parse_reports as parse_reports

log = get_logger(__name__)

"""
Stage: Parse reports.

Parses summary HMM reports and writes them into the SQLite database.
"""
def parse_reports_to_database(config: Config) -> None:
    """Parse summary hmmreport and insert into the database.

    Args:
        options: Configuration with `.database_directory`
            and report file paths.
    """
    print_header("Parse reports to database", logger=log)
    parse_reports.main_parse_summary_hmmreport(config)

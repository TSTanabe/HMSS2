# src/hmsss/__main__.py
from __future__ import annotations

import os
import sys

from hmsss.cli.parse import parse_arguments
from hmsss.core.runner import run_pipeline
from hmsss.core.logging import setup_logging


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    opts = parse_arguments(argv)

    # zentrales Log
    log_file = os.path.join(opts.result_files_directory, "execution_logfile.txt")
    setup_logging(getattr(opts, "verbose", 1), log_file)

    run_pipeline(opts)


if __name__ == "__main__":
    main()

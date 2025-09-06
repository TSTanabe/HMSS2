# src/hmsss/__main__.py
from __future__ import annotations

import sys

from hmsss.cli.parse import parse_to_config
from hmsss.core.runner import run_pipeline


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    config = parse_to_config(argv)

    # This function organizes the whole pipeline
    run_pipeline(config)


if __name__ == "__main__":
    main()
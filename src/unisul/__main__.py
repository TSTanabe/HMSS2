from __future__ import annotations

import sys

from unisul.data_management import DataManagementError, data_cli


def _split_global_data_dir(argv: list[str]) -> tuple[list[str], str | None]:
    """Remove global --data-dir from argv so it can be used by data subcommands."""
    remaining: list[str] = []
    data_dir: str | None = None
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--data-dir":
            if i + 1 >= len(argv):
                raise DataManagementError("--data-dir requires a path")
            data_dir = argv[i + 1]
            i += 2
        elif arg.startswith("--data-dir="):
            data_dir = arg.split("=", 1)[1]
            i += 1
        else:
            remaining.append(arg)
            i += 1
    return remaining, data_dir


def main(argv: list[str] | None = None) -> None:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        stripped, data_dir = _split_global_data_dir(arguments)
        if stripped and stripped[0] == "data":
            data_cli(stripped[1:], cli_path=data_dir)
            return

        from unisul.cli.parse import parse_to_config
        from unisul.core.runner import run_pipeline

        run_pipeline(parse_to_config(arguments))
    except DataManagementError as exc:
        raise SystemExit(f"unisul: error: {exc}") from exc

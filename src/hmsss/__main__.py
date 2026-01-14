# src/hmsss/__main__.py
from __future__ import annotations

import sys

from hmsss.cli.parse import parse_to_config
from hmsss.core.runner import run_pipeline
from hmsss.core.load_data import ensure_data_dir_with_marker
from hmsss.core.load_data import ensure_gpkg_dir_with_marker

DATA_ZIP_URL = "https://zenodo.org/records/XXXX/files/hmss2_data_v1.zip?download=1"
GPKG_ZIP_URL = "https://.../HMSS2_gpkg.zip"

# optional, wenn du Hashes hast
DATA_SHA256 = None
GPKG_SHA256 = None


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    config = parse_to_config(argv)
    """
    ensure_data_dir_with_marker(
        data_zip_url=DATA_ZIP_URL,
        expected_sha256=DATA_SHA256,
    )

    # GPKG data ca 8 GB. Load when required
    if config.use_read_mapping:
        #ensure_gpkg_dir_with_marker(
            gpkg_zip_url=GPKG_ZIP_URL,
            expected_sha256=GPKG_SHA256,
        )
    """
    # This function organizes the whole pipeline
    run_pipeline(config)


if __name__ == "__main__":
    main()

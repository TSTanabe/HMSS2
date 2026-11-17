import glob
import os
import shutil
import sys
from pathlib import Path
from typing import List, Sequence, Union

from unisul.cli.config import INTERMEDIATE_HITS_SUFFIX, Config
from unisul.cli.paths import BIN_DIR  # nutzt UNISUL_BIN_DIR
from unisul.core.logging import get_logger

log = get_logger(__name__)


def get_all_files(directory: str, ending: Union[str, int] = 0) -> List[str]:
    """
    Recursively retrieves files from a directory.

    Args:
        directory (str): Root directory to search.
        ending (str or int): Optional file suffix to filter by.

    Returns:
        List[str]: List of matching file paths.
    """
    matched = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            file = os.path.join(path, name)
            if ending == 0 or file.endswith(ending):
                matched.append(file)
    return matched


def collect_intermediate_hits(config: Config) -> dict[str, str]:
    """
    Collect all intermediate HMM-hit FASTA files.

    Expected filenames:
        <domain_type>{INTERMEDIATE_HITS_SUFFIX}

    Returns
    -------
    dict[str, str]
        Mapping:
            domain_type -> intermediate FASTA path

        Example:
            {
                "MorD": "/.../MorD.intermediate_hits_faa",
                "AtrD": "/.../AtrD.intermediate_hits_faa",
            }
    """
    intermediate_files = glob.glob(
        os.path.join(
            config.cross_check_directory,
            f"*{INTERMEDIATE_HITS_SUFFIX}",
        )
    )

    intermediate_dict = {
        os.path.basename(path).removesuffix(INTERMEDIATE_HITS_SUFFIX): path
        for path in intermediate_files
    }

    log.debug(
        "Collected %d intermediate hit FASTA files from %s",
        len(intermediate_dict),
        config.cross_check_directory,
    )

    return intermediate_dict


def get_genome_id(path: str) -> str:
    """Extracts genome ID from filename (before first dot)."""
    return os.path.basename(path).split(".")[0]


def taxonomy_lineage(array: List[str], trennzeichen: str) -> str:
    """
    Joins taxonomy names with separator, replacing spaces with dashes.

    Returns:
        str: A single string like 'Bacteria-Firmicutes-Bacilli'.
    """
    if not array or not all(isinstance(x, str) for x in array):
        return "NoTaxonomy"

    return trennzeichen.join(array).replace(" ", "-")


def get_executable_dir() -> str:
    """Returns directory of script (or executable in case of PyInstaller)."""
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def find_executable(executable: str | Sequence[str]) -> str:
    """
    Locate an executable in PATH, then the configured UNISUL bin directory.
    Raises FileNotFoundError if not found.
    """
    candidates = (executable,) if isinstance(executable, str) else tuple(executable)
    for candidate in candidates:
        path = shutil.which(candidate)
        if path:
            log.debug("Found executable in PATH: %s", path)
            return path

        bin_candidate = Path(BIN_DIR) / candidate
        if bin_candidate.is_file() and os.access(str(bin_candidate), os.X_OK):
            log.debug("Found executable in BIN_ROOT: %s", bin_candidate)
            return str(bin_candidate)

    log.error(
        "Executable not found: %s. Searched PATH and BIN_ROOT=%s",
        ", ".join(candidates),
        BIN_DIR,
    )
    raise FileNotFoundError(f"{', '.join(candidates)} executable not found.")

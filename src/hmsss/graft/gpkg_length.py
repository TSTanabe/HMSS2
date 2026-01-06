import os
from typing import Dict, Iterable
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def _iter_fasta_lengths(faa_path: str) -> Iterable[int]:
    """Iteriere über Sequenzlängen einer FASTA-Datei (ungapped)."""
    length = 0
    with open(faa_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if length > 0:
                    yield length
                length = 0
            else:
                length += len(line.strip())
        if length > 0:
            yield length


def _median(values: list[int]) -> int:
    if not values:
        raise ValueError("Cannot compute median of empty list")
    values.sort()
    n = len(values)
    mid = n // 2
    if n % 2 == 1:
        return values[mid]
    else:
        return (values[mid - 1] + values[mid]) // 2


def collect_gpkg_reference_median_lengths(
    gpkg_packages: Dict[str, str],
) -> Dict[str, int]:
    """
    For each GPKG, read <gpkg_name>.faa and compute median sequence length.

    Returns:
        { gpkg_name : median_length }
    """
    logger.info("Initializing median sequence length calculation for .gpkg packages")
    result: Dict[str, int] = {}

    for gpkg_name, gpkg_dir in gpkg_packages.items():
        faa_path = os.path.join(gpkg_dir, f"{gpkg_name}.faa")
        faa_path2 = os.path.join(gpkg_dir, f"refseq_database.faa")
        try:
            if os.path.isfile(faa_path) and os.path.getsize(faa_path) > 0:
                lengths = list(_iter_fasta_lengths(faa_path))
                result[gpkg_name] = _median(lengths)
            elif os.path.isfile(faa_path2) and os.path.getsize(faa_path2) > 0:
                lengths = list(_iter_fasta_lengths(faa_path2))
                result[gpkg_name] = _median(lengths)
            else:
                logger.warning(f"{gpkg_name}.gpkg Reference FASTA missing {faa_path}")
                continue

        except Exception as e:
            logger.error(f"Median length determination for gpkg packages failed {e}")
    logger.debug(f"Lengths per package {result}")
    return result

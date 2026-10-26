import csv
import shutil
from collections import defaultdict
from pathlib import Path

from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def collect_placement_search_hmms(config) -> dict:
    """
    Match search HMMs to protein types supported by placement packages.

    Search HMM names are derived from filenames by removing everything
    before and including the first underscore.

    Only HMMs located under DisCo, DiSuCy, HMSS2 or Basis are considered.

    Returns:
        {
            package_path: {
                protein_type: [search_hmm_path, ...]
            }
        }

    No files are modified or concatenated.
    """

    placement_root = Path(config.paths.placement).resolve()
    data_root = Path(config.paths.data).resolve()

    allowed_directories = {"basis", "disco", "disucy", "hmss2"}

    # ---------------------------------------------------------
    # 1. Collect protein types from placement packages
    # ---------------------------------------------------------

    package_types = {}

    for package in sorted(placement_root.iterdir()):
        if not package.is_dir():
            continue

        annotations = package / "annotations.tsv"
        manifest = package / "manifest.json"

        if not annotations.is_file() or not manifest.is_file():
            continue

        with annotations.open(encoding="utf-8") as handle:
            protein_types = {
                row[1].strip()
                for row in csv.reader(handle, delimiter="\t")
                if len(row) >= 2 and row[1].strip()
            }

        if protein_types:
            package_types[package] = protein_types

    all_types = set().union(*package_types.values()) if package_types else set()

    logger.info(
        "Found %d placement packages supporting %d protein types",
        len(package_types),
        len(all_types),
    )

    # ---------------------------------------------------------
    # 2. Index available search HMMs by filename
    # ---------------------------------------------------------

    hmm_index = defaultdict(set)

    for hmm_file in data_root.rglob("*.hmm"):
        # Consider only selected resource directories
        relative_parts = {
            part.lower() for part in hmm_file.relative_to(data_root).parts[:-1]
        }

        if not relative_parts.intersection(allowed_directories):
            continue

        # Remove prefix before the first underscore
        filename = hmm_file.stem

        if "_" not in filename:
            continue

        hmm_name = filename.split("_", 1)[1]

        # Match against protein types supported by placement packages
        if hmm_name in all_types:
            hmm_index[hmm_name].add(hmm_file.resolve())

    logger.info(
        "Found search HMMs for %d of %d placement protein types",
        len(hmm_index),
        len(all_types),
    )

    # ---------------------------------------------------------
    # 3. Associate search HMMs with placement packages
    # ---------------------------------------------------------

    result = {}

    for package, protein_types in package_types.items():
        matched = {
            protein_type: sorted(hmm_index[protein_type])
            for protein_type in sorted(protein_types)
            if hmm_index.get(protein_type)
        }

        if not matched:
            logger.warning(
                "No matching search HMMs for placement package %s",
                package.name,
            )
            continue

        missing = protein_types - matched.keys()

        if missing:
            logger.warning(
                "Placement package %s: no search HMMs for %s",
                package.name,
                ", ".join(sorted(missing)),
            )

        result[package] = matched

        logger.info(
            "Placement %-25s | %d protein types | %d search HMMs",
            package.name,
            len(matched),
            sum(len(files) for files in matched.values()),
        )

    # ---------------------------------------------------------
    # 4. Summary
    # ---------------------------------------------------------

    unique_types = (
        set().union(*(types.keys() for types in result.values())) if result else set()
    )

    unique_hmms = {
        hmm for types in result.values() for files in types.values() for hmm in files
    }

    logger.info(
        "Stage 30 resources: %d placement packages, "
        "%d protein types, %d unique search HMM files",
        len(result),
        len(unique_types),
        len(unique_hmms),
    )

    return result


def create_placement_search_hmm_library(
    config,
    result: dict,
) -> Path:
    """
    Concatenate all unique search HMMs collected from placement packages
    into the HMM library specified by config.library.

    Parameters
    ----------
    config:
        HMSS3 configuration object.

    result:
        Output of collect_placement_search_hmms(config):

        {
            package_path: {
                protein_type: [hmm_path_1, hmm_path_2, ...]
            }
        }

    overwrite:
        If True, overwrite an existing HMM library.
        If False, raise FileExistsError if the library exists.

    Returns
    -------
    Path:
        Path to the generated combined HMM library.
    """

    library_path = Path(config.library)

    hmm_files = {
        Path(hmm_file).resolve()
        for protein_hmms in result.values()
        for files in protein_hmms.values()
        for hmm_file in files
    }

    with library_path.open("wb") as output:
        for hmm_file in hmm_files:
            with hmm_file.open("rb") as source:
                shutil.copyfileobj(source, output)

            # Ensure separation between HMM files.
            output.write(b"\n")

    return library_path

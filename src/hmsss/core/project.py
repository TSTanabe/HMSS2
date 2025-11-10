#!/usr/bin/python
from __future__ import annotations
import csv
import os

from datetime import datetime
from dataclasses import is_dataclass, asdict
from pathlib import Path
from typing import Any, Mapping, Iterable, Tuple, Optional

from hmsss.core.logging import get_logger

logger = get_logger(__name__)

"""
Project setup and configuration utilities for HMSSS.

Handles creation of result directories, discovery of existing databases,
and export of configuration parameters to flat TSV files. Provides helpers
to flatten nested dataclass/Namespace structures for logging or saving.
"""


def prepare_result_space(config, project: str = "project") -> None:
    """
    Creates and sets up the results directory for a project, including all needed subdirectories.
    - If no custom `-r` result directory is given, creates a new timestamped
      subdirectory under the default results path.
    - If a result directory is given, checks for an existing database; if
      found, reuses it; otherwise creates a new project directory there.
    - Initializes subdirectories (Hit_list, Sequences, Collinear_syntenic_blocks,
      Filtered_hits).
    - Updates multiple attributes of `config` with new file and directory paths.
    - Writes a TSV summary of the configuration parameters.

    Args:
        config (Options): Has .location, .result_files_directory, etc.
        project (str): Project name suffix (default: 'project')

    Output:
        Updates multiple attributes of options (see docstring).

    Example Input:
        options.location = '/home/user/myproject'
        options.result_files_directory = '/home/user/myproject/results'
    """

    # Use standard result directory if -r was not used
    if config.paths.results == config.cli_result_dir_in:
        # Make a new project in default directory
        config.result_files_directory = create_project(config.paths.results, project)
        new_project = True

    else:
        # -r was given search for existing project
        tuple_database_dir_database_path = find_database_in_directory(
            config.cli_result_dir_in
        )
        if tuple_database_dir_database_path:
            # Existing database was found, first occurrence of .db is used
            database_dir, database_path = tuple_database_dir_database_path
            config.result_files_directory = str(database_dir)

        else:
            # No existing was found, save new project to given location
            config.result_files_directory = create_project(
                config.cli_result_dir_in, project
            )
            new_project = True

    # The project directory is now saved in config.result_files_directory
    # Setup the project folders
    config.database_directory = config.result_files_directory + "/database.db"
    config.fasta_initial_hit_directory = config.result_files_directory + "/Hit_list"
    config.fasta_output_directory = config.result_files_directory + "/Sequences"
    config.csb_directory = config.result_files_directory + "/Collinear_syntenic_blocks"
    config.cross_check_directory = config.result_files_directory + "/Filtered_hits"

    # Output files
    config.glob_report = os.path.join(
        config.result_files_directory, "global_report.cat_hmmreport"
    )
    config.glob_trusted_hitreport = (
        config.result_files_directory + "/global_trusted_hits_summary.db"
    )
    config.glob_intermediate_hitreport = (
        config.result_files_directory + "/global_intermediate_hits_summary.db"
    )
    config.csb_output_file = config.csb_directory + "/Csb_output.txt"
    config.gene_clusters_file = config.csb_directory + "/All_gene_clusters.txt"

    # Create required directories
    for path in [
        config.fasta_initial_hit_directory,
        config.fasta_output_directory,
        config.cross_check_directory,
        config.csb_directory,
    ]:
        if not os.path.exists(path):
            os.mkdir(path)

    # Write down setting
    write_config_to_tsv(config, config.result_files_directory)
    # 5. If using a pre-existing project, force pipeline to start at stage 3
    # if new_project and config.stage < 3:
    #    logger.warning("Existing project directory detected. Setting start stage to 4.")
    #    config.stage = 4

    return


def create_project(directory, projectname="project") -> str:
    """Create a new project directory with a timestamped name.

    Example: `2025-04-16_14-53-21_project`

    Args:
        directory: Parent directory where the project folder is created.
        projectname: Project suffix (default: "project").

    Returns:
        Path of the created directory as string.
    """
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")  # z. B. "2025-04-16_14-53-21"
    directory = os.path.join(directory, f"{timestamp}_{projectname}")

    try:
        os.mkdir(directory)
    except Exception:
        logger.error("Creation of project directory failed, no writing rights")

    return directory


def find_database_in_directory(
    directory: str | Path,
    db_name: str = "database.db",
) -> Optional[Tuple[Path, Path]]:
    """Search for a database file under the given directory.

    - First checks if `db_name` exists directly under `directory`.
    - If not found, searches recursively and returns the first hit.

    Args:
        directory: Directory to search in.
        db_name: Database filename to search for (default: "database.db").

    Returns:
        Tuple of (database directory, database file path) if found, else None.
    """
    root = Path(directory).expanduser().resolve()

    # 1) Erwartete Position: direkt in `directory`
    candidate = root / db_name
    if candidate.is_file():
        return candidate.parent, candidate

    logger.warning(
        "No database file at expected location: %s. Searching recursively under %s.",
        candidate,
        root,
    )

    # 2) Rekursiv suchen – erster Treffer
    for db_path in root.rglob(db_name):
        logger.info("Found database at %s", db_path)
        return db_path.parent, db_path

    logger.warning("No database named %r found under %s.", db_name, root)
    return None


def any_process_args_provided(args, default_values: dict) -> bool:
    """Check if any process-related CLI arguments differ from their defaults.

    Args:
        args: Argument object (e.g. argparse.Namespace).
        default_values: Mapping of argument name -> default value.

    Returns:
        True if any argument is not equal to its default, otherwise False.
    """
    for arg, default in default_values.items():
        if getattr(args, arg) != default:
            # print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False


# Routinen für den config file
def write_config_to_tsv(
    config: Any,
    output_directory: str | os.PathLike,
    filename: str = "parameters_summary.tsv",
) -> str:
    """Write a (nested) Config object into a flat TSV file.

    - Dataclasses are recursively converted using `asdict()`.
    - Nested structures are flattened into dot-path keys.
    - Lists/Tuples/Sets are serialized as comma-separated strings.
    - Path objects are converted to strings.

    Args:
        config: Configuration object (possibly nested dataclasses).
        output_directory: Directory where the TSV file is saved.
        filename: Output filename (default: "parameters_summary.tsv").

    Returns:
        Path to the written TSV file as string.
    """
    os.makedirs(output_directory, exist_ok=True)
    output_path = os.path.join(output_directory, filename)

    flat = _flatten_config(config)

    with open(output_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Parameter", "Value"])
        for k in sorted(flat):
            w.writerow([k, flat[k]])

    logger.info("Parameters saved: %s", output_path)
    return output_path


def _flatten_config(obj: Any, prefix: str = "") -> dict[str, str]:
    """Recursively flatten a configuration object into dot-path keys.

    Handles dataclasses, argparse.Namespace, mappings, iterables, and scalars.
    Paths are converted to strings. Non-scalar iterables fall back to `repr()`.

    Args:
        obj: Object to flatten.
        prefix: Dot-path prefix used for recursion.

    Returns:
        Mapping from dot-path keys to stringified values.
    """
    # Dataclass → Dict (rekursiv)
    if is_dataclass(obj):
        obj = asdict(obj)

    # argparse.Namespace oder normales Objekt mit __dict__
    if hasattr(obj, "__dict__") and not isinstance(obj, (Path, str, bytes)):
        obj = vars(obj)

    # Mapping (dict-ähnlich)
    if isinstance(obj, Mapping):
        out: dict[str, str] = {}
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.update(_flatten_config(v, key))
        return out

    # Iterable (Liste/Tuple/Set), aber keine Strings/Bytes
    if isinstance(obj, Iterable) and not isinstance(obj, (str, bytes, Path)):
        # Versuche, skalare Elemente kommasepariert zu schreiben
        try:
            items = list(obj)
            if all(_is_scalar(x) for x in items):
                return {prefix: ",".join(_to_str(x) for x in items)}
        except TypeError:
            pass  # Nicht indexierbar → repr unten
        # Fallback: repr() der Struktur
        return {prefix: repr(obj)}

    # Skalar oder Path
    return {prefix: _to_str(obj)}


def _is_scalar(x: Any) -> bool:
    """Return True if `x` is a scalar value (str, int, float, bool, None, Path)."""
    return isinstance(x, (str, bytes, int, float, bool, type(None), Path))


def _to_str(x: Any) -> str:
    """Convert supported objects to strings.

    - Path → string
    - None → empty string
    - Otherwise → str(x)
    """
    if isinstance(x, Path):
        return str(x)
    if x is None:
        return ""
    return str(x)

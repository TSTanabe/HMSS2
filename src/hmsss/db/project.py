#!/usr/bin/python

import csv
import os
import sys
from datetime import datetime
from dataclasses import is_dataclass, asdict
from pathlib import Path
from typing import Any, Mapping, Iterable

from hmsss.core.logging import get_logger
logger = get_logger(__name__)


def prepare_result_space(config, project: str = "project") -> None:
    """
    Creates and sets up the results directory for a project, including all needed subdirectories.

    Args:
        config (Options): Has .location, .result_files_directory, etc.
        project (str): Project name suffix (default: 'project')

    Output:
        Updates multiple attributes of options (see docstring).

    Example Input:
        options.location = '/home/user/myproject'
        options.result_files_directory = '/home/user/myproject/results'
    """

    now = datetime.now()
    timestamp = str(datetime.timestamp(now))

    # Use standard result directory if -r was not used
    if config.paths.results == config.result_dir_in:
        # Make a new project in default directory
        config.result_files_directory = create_project(config.paths.results, project)
        write_config_to_tsv(config, config.result_files_directory)

    # User-defined directory: check for project
    elif not is_existing_project_directory(config):
        if not os.path.isdir(config.result_files_directory):
            try:
                os.mkdir(config.result_files_directory)
                logger.info(f"Created results dir: {config.result_files_directory}")
            except Exception as e:
                logger.error(
                    f"No writing rights for directory {config.result_files_directory}\n {e}"
                )
                sys.exit(1)

        config.result_files_directory = create_project(
            config.result_files_directory, project
        )
        config.new_project = True
        write_config_to_tsv(config, config.result_files_directory)

        # Project structure, applies in both cases
        # Directories
    config.database_directory = config.result_files_directory + "/database.db"
    config.fasta_initial_hit_directory = config.result_files_directory + "/Hit_list"
    config.fasta_output_directory = config.result_files_directory + "/Sequences"
    config.csb_directory = (
            config.result_files_directory + "/Collinear_syntenic_blocks"
    )
    config.cross_check_directory = config.result_files_directory + "/Filtered_hits"

    # Output files
    if config.glob_report is None or not os.path.isfile(config.glob_report):
        config.glob_report = os.path.join(
            config.result_files_directory, "global_report.cat_hmmreport"
        )

    config.glob_trusted_hitreport = (
            config.result_files_directory + "/global_trusted_hits_summary.hmmreport"
    )
    config.glob_intermediate_hitreport = (
            config.result_files_directory + "/global_intermediate_hits_summary.hmmreport"
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

    # 5. If using a pre-existing project, force pipeline to start at stage 3
    # if not options.new_project and options.stage < 3:
    #    logger.warning("Existing project directory detected. Setting start stage to 4.")
    #    options.stage = 4

    return


def create_project(directory, projectname="project"):
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")  # z. B. "2025-04-16_14-53-21"
    directory = os.path.join(directory, f"{timestamp}_{projectname}")

    try:
        os.mkdir(directory)
    except Exception:
        logger.error("Creation of project directory failed, no writing rights")

    return directory


def is_existing_project_directory(config) -> bool:
    """
    Checks if a result_files_directory contains a valid project structure.
    If a database is found, adjusts options accordingly.

    Args:
        config (Options): Has .result_files_directory etc.

    Returns:
        bool: True if project folder found or created, else False.
    """

    # Database may be provided, then consider this as a result directory
    if config.database_directory and os.path.isfile(config.database_directory):
        config.result_files_directory = os.path.dirname(config.database_directory)

    try:
        database_path = os.path.join(config.result_files_directory, "database.db")

        if not os.path.isfile(database_path):
            logger.warning(
                "No database file found in expected directory. Searching recursively..."
            )

            # Recursive search for database.db
            found_db = None
            for root, dirs, files in os.walk(config.result_files_directory):
                if "database.db" in files:
                    found_db = os.path.join(root, "database.db")
                    break

            if found_db:
                logger.info(f"Found database at {found_db}")
                config.result_files_directory = os.path.dirname(found_db)
            else:
                logger.warning(
                    "No database found at all. Existing project was not found."
                )
                return False

        # Now ensure required directories exist
        required_dirs = [
            "Sequences",
            "Filtered_hits",
            "Hit_list",
            "Collinear_syntenic_blocks",
        ]

        for subdir in required_dirs:
            path = os.path.join(config.result_files_directory, subdir)
            if not os.path.isdir(path):
                logger.info(f"Missing directory {subdir} created")
                os.makedirs(path, exist_ok=True)

        self_query_path = os.path.join(config.result_files_directory, "self_blast.faa")
        if os.path.isfile(self_query_path):
            config.self_query = self_query_path
        else:
            logger.warning(
                "No internal query file (self_blast.faa) found. Will need to create it later."
            )

        # Find blast table if not already defined
        if config.glob_report is None:
            for file_name in os.listdir(config.result_files_directory):
                if file_name.startswith("global_report.cat_hmmreport"):
                    config.glob_report = os.path.join(
                        config.result_files_directory, file_name
                    )
                    break

    except Exception as e:
        logger.error(
            f"An error occurred while checking for existing project folder: {e}"
        )
        sys.exit(1)

    return True


def any_process_args_provided(args, default_values: dict) -> bool:
    """
    Checks if any CLI/process arguments are different from their default values.

    Args:
        args: The argument object (e.g., argparse.Namespace)
        default_values (dict): Mapping of argname -> default

    Returns:
        bool: True if any arg has a value different from default, else False.
    """
    for arg, default in default_values.items():
        if getattr(args, arg) != default:
            # print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False



def write_config_to_tsv(
    config: Any,
    output_directory: str | os.PathLike,
    filename: str = "parameters_summary.tsv",
) -> str:
    """
    Schreibt ein (verschachteltes) Config-Objekt flach als TSV.
    - Dataclasses werden rekursiv mit asdict() aufgelöst.
    - Verschachtelte Strukturen werden mit Dot-Pfaden als Keys ausgegeben.
    - Listen/Tuples/Sets werden kommasepariert.
    - Path-Objekte werden zu Strings konvertiert.

    Returns:
        Pfad zur geschriebenen TSV-Datei als String.
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
    """
    Rekursive Flachlegung nach 'dot path'-Schlüssel → String-Werte.
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
    return isinstance(x, (str, bytes, int, float, bool, type(None), Path))

def _to_str(x: Any) -> str:
    if isinstance(x, Path):
        return str(x)
    if x is None:
        return ""
    return str(x)

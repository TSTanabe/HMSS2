#!/usr/bin/python

import os
from datetime import datetime
from typing import Union, Any, Dict




def prepare_result_space(options: Any, project: str = "project") -> None:
    """
    Prepares the result directory structure for a new or existing project.

    If the given result directory already contains a project, this function will
    set up the paths to the existing subdirectories and files. Otherwise, it will
    create a new project directory and initialize all required subdirectories and files.

    Args:
        options: An object with at least these attributes:
            - result_files_directory (str): Path to the top‐level results directory.
        project: Name of the project (used as folder name if creating new).

    Raises:
        OSError: If the function cannot create the specified directories due to
                 insufficient permissions or other filesystem errors.
    """
    # Check if this is already a project folder
    if is_project_folder(options):
        # Existing project: simply assign subdirectory and file paths
        base_dir = options.result_files_directory
        options.database_directory = os.path.join(base_dir, "database.db")
        options.fasta_initial_hit_directory = os.path.join(base_dir, "Hit_list")
        options.fasta_output_directory = os.path.join(base_dir, "Sequences")
        options.csb_directory = os.path.join(base_dir, "Collinear_syntenic_blocks")
        options.cross_check_directory = os.path.join(base_dir, "Filtered_hits")

        options.divergent_output_file = os.path.join(base_dir, "div_output_file.faa")
        options.csb_output_file = os.path.join(options.csb_directory, "Csb_output.txt")
        options.gene_clusters_file = os.path.join(
            options.csb_directory, "All_gene_clusters.txt"
        )
    else:
        # New project: create the top‐level directory if it doesn't exist
        top_dir = options.result_files_directory
        if not os.path.isdir(top_dir):
            try:
                os.mkdir(top_dir)
            except OSError as exc:
                raise OSError(f"[ERROR] Cannot create result directory '{top_dir}': {exc}")

        # Create (or rename) the project directory under the specified top‐level path
        project_dir = create_project(top_dir, project)
        options.result_files_directory = project_dir

        # Assign paths for new subdirectories
        options.database_directory = os.path.join(project_dir, "database.db")

        options.fasta_initial_hit_directory = os.path.join(project_dir, "Hit_list")
        os.mkdir(options.fasta_initial_hit_directory)

        options.fasta_output_directory = os.path.join(project_dir, "Sequences")
        os.mkdir(options.fasta_output_directory)

        options.csb_directory = os.path.join(project_dir, "Collinear_syntenic_blocks")
        os.mkdir(options.csb_directory)

        options.cross_check_directory = os.path.join(project_dir, "Filtered_hits")
        os.mkdir(options.cross_check_directory)

        # Assign paths for output files
        options.divergent_output_file = os.path.join(project_dir, "div_output_file.faa")
        options.csb_output_file = os.path.join(options.csb_directory, "Csb_output.txt")
        options.gene_clusters_file = os.path.join(
            options.csb_directory, "All_gene_clusters.txt"
        )

        


def create_project(directory: str, projectname: str = "project") -> str:
    """
    Creates a new project folder inside the given directory. The folder name
    will be prefixed by a timestamp in the format YYYY-MM-DD_HH-MM-SS.

    Args:
        directory: Path to the parent directory where the project folder should be created.
        projectname: Name of the project (appended after the timestamp). Defaults to "project".

    Returns:
        The full path to the newly created project directory.

    Raises:
        OSError: If the directory cannot be created due to permission issues or other filesystem errors.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    project_dir = os.path.join(directory, f"{timestamp}_{projectname}")

    try:
        os.mkdir(project_dir)
    except OSError as exc:
        raise OSError(f"ERROR: Cannot create project directory '{project_dir}': {exc}") from exc

    return project_dir
    

def is_project_folder(options: Any) -> bool:
    """
    Checks whether the given options.result_files_directory represents an existing
    project folder (i.e., contains the expected subdirectories and database file).

    If options.database_directory points to an existing file, its parent directory
    is assumed to be the project folder.

    Args:
        options: Object with at least these attributes:
            - database_directory (str or None): Path to an existing database file.
            - result_files_directory (str): Path to the directory to check.

    Returns:
        True if the directory is a valid project folder (contains 'database.db',
        'Sequences', 'Hit_list', and 'Collinear_syntenic_blocks'); otherwise False.

    Raises:
        OSError: If there is an unexpected filesystem error while checking.
    """
    # If a database file is provided and exists, use its parent directory.
    if getattr(options, "database_directory", None) and os.path.isfile(options.database_directory):
        options.result_files_directory = os.path.dirname(options.database_directory)

    base_dir = options.result_files_directory

    try:
        # Check for 'database.db'
        db_path = os.path.join(base_dir, "database.db")
        if not os.path.isfile(db_path):
            print("No database file found. Creating new project folder.")
            return False

        # Check for 'Sequences' directory
        sequences_dir = os.path.join(base_dir, "Sequences")
        if not os.path.isdir(sequences_dir):
            print("No 'Sequences' directory found. Creating new project folder.")
            return False

        # Check for 'Hit_list' directory
        hit_list_dir = os.path.join(base_dir, "Hit_list")
        if not os.path.isdir(hit_list_dir):
            print("No 'Hit_list' directory found. Creating new project folder.")
            return False

        # Check for 'Collinear_syntenic_blocks' directory
        csb_dir = os.path.join(base_dir, "Collinear_syntenic_blocks")
        if not os.path.isdir(csb_dir):
            print("No 'Collinear_syntenic_blocks' directory found. Creating new project folder.")
            return False

    except OSError as exc:
        raise OSError(f"Error checking project folder '{base_dir}': {exc}") from exc

    return True
    
    
def any_process_args_provided(args: Any, default_values: Dict[str, Any]) -> bool:
    """
    Checks if any argument value in `args` differs from its default in `default_values`.

    Args:
        args: Namespace or object with attributes corresponding to argument names.
        default_values: Dictionary mapping argument names (str) to their default values.

    Returns:
        True if at least one argument in `args` does not equal its default value;
        otherwise False.
    """
    for arg_name, default_val in default_values.items():
        if getattr(args, arg_name) != default_val:
            return True

    return False
    

   
     







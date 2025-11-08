#!/usr/bin/python
import sys
from hmsss.core.logging import get_logger

logger = get_logger(__name__)

def print_command_line_args(output_file: str) -> None:
    """Writes the command-line arguments to a file.

    Args:
        output_file (str): Path to the file to write arguments to.

    """
    try:
        with open(output_file, "w") as f:
            f.write("Command-line arguments passed to the script:\n")
            for index, arg in enumerate(sys.argv):
                f.write(f"Argument {index}: {arg}\n")
    except Exception as e:
        logger.error(f"Failed to write to file {output_file}: {e}")


def print_file_content(file_path: str) -> None:
    """Prints the content of a file along with its path.

    Args:
        file_path (str): Path to the file.

    Example:
        >>> print_file_content("csb_patterns.txt")
    """
    try:
        with open(file_path, "r") as file:
            content = file.read()
        logger.info(f"File Path: {file_path}")
        logger.info("File Content:")
        logger.info(content)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
    except Exception as e:
        logger.error(f"An error occurred while reading the file: {e}")

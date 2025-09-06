#!/usr/bin/python
import sys
from hmsss.core.logging import get_logger

logger = get_logger(__name__)

def print_command_line_args(output_file: str) -> None:
    """Writes the command-line arguments to a file.

    Args:
        output_file (str): Path to the file to write arguments to.

    Example:
        >>> print_command_line_args("args.txt")
    """
    try:
        with open(output_file, "w") as f:
            f.write("Command-line arguments passed to the script:\n")
            for index, arg in enumerate(sys.argv):
                f.write(f"Argument {index}: {arg}\n")
    except Exception as e:
        logger.error(f"Failed to write to file {output_file}: {e}")
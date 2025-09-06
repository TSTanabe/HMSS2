# src/hmsss/core/logging.py
from __future__ import annotations

import logging
import os
import sys
from typing import Optional

"""
Logging setup utilities for HMSSS.

Provides a unified package logger (`hmsss`) with console and optional file
handlers. Also offers convenience functions to get child loggers and to
print formatted section headers into the logs.
"""


_LEVELS = {
    0: logging.WARNING,  # quiet
    1: logging.INFO,  # default
    2: logging.DEBUG,  # verbose
}


def setup_logging(
    verbosity: int = 1, logfile: Optional[str] = None, *, force: bool = False
) -> logging.Logger:
    """Initialize the package logger `hmsss`.

    Sets up a console handler and, optionally, a file handler. Existing
    handlers can be replaced if `force=True`.

    Args:
        verbosity: Logging level (0 = WARNING, 1 = INFO, 2 = DEBUG).
        logfile: Path to a log file. If given, messages are always logged
            at DEBUG level. The directory is created if it does not exist.
        force: If True, remove existing handlers and rebuild the logger.

    Returns:
        Configured `logging.Logger` for the package root.
    """
    logger = logging.getLogger("hmsss")

    if force:
        for h in list(logger.handlers):
            logger.removeHandler(h)

    if logger.handlers and not force:
        # Logger existiert schon → nur Level nachziehen und zurück
        level = _LEVELS.get(verbosity, logging.INFO)
        logger.setLevel(level)
        for h in logger.handlers:
            if isinstance(h, logging.StreamHandler):
                h.setLevel(level)
        return logger

    level = _LEVELS.get(verbosity, logging.INFO)
    logger.setLevel(level)

    # Set formats for console and logging file
    console_format = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_format = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(filename)s:%(lineno)d | %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(console_format)
    logger.addHandler(ch)

    # File (optional, immer DEBUG, inkl. Ordner anlegen)
    if logfile:
        os.makedirs(os.path.dirname(logfile), exist_ok=True)
        fh = logging.FileHandler(logfile, mode="a", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(file_format)
        logger.addHandler(fh)

    logger.propagate = False  # keine Weitergabe an root
    return logger


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a child logger of the package logger.

    Args:
        name: Optional sub-name, e.g. `"stages.initial_search"`. If None,
            returns the root logger `"hmsss"`.

    Returns:
        Logger instance under `"hmsss"` namespace.
    """
    base = "hmsss" if not name else f"hmsss.{name}"
    return logging.getLogger(base)


def print_header(
    text: str, *, char: str = "=", logger: Optional[logging.Logger] = None
) -> None:
    """Log a formatted header line.

    Useful to visually separate major steps in the logfile.

    Args:
        text: Header text to print.
        char: Character used for framing the header (default: "=").
        logger: Logger to use. Defaults to a child logger of `"hmsss"`.
    """
    log = logger or get_logger(__name__)
    text: str = 5 * char + f" {text} " + 5 * char
    log.info(text)


__all__ = ["setup_logging", "get_logger", "print_header"]

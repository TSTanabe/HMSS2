# src/hmsss/core/logging.py
from __future__ import annotations

import logging
import os
import sys
from typing import Optional

_LEVELS = {
    0: logging.WARNING,  # quiet
    1: logging.INFO,  # default
    2: logging.DEBUG,  # verbose
}


def setup_logging(
    verbosity: int = 1, logfile: Optional[str] = None, *, force: bool = False
) -> logging.Logger:
    """
    Initialisiert den Paket-Logger 'hmsss' mit Console- und optionalem File-Handler.
    - verbosity: 0=WARNING, 1=INFO, 2=DEBUG
    - logfile: Pfad zur Logdatei (wird im DEBUG-Level beschrieben)
    - force: vorhandene Handler entfernen & neu aufbauen
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
    """
    Hole einen Child-Logger unterhalb 'hmsss' (z. B. 'hmsss.stages.initial_search').

    Returns:
        logging.Logger:
    """
    base = "hmsss" if not name else f"hmsss.{name}"
    return logging.getLogger(base)


def print_header(
    text: str, *, char: str = "=", logger: Optional[logging.Logger] = None
) -> None:
    """
    Write headers into logfile
    """
    log = logger or get_logger(__name__)
    text: str = 5 * char + f" {text} " + 5 * char
    log.info(text)



__all__ = ["setup_logging", "get_logger", "print_header"]

#!/usr/bin/python

import os

from multiprocessing import Pool, Value, Lock

from hmsss.cli.config import Config
from hmsss.core.logging import get_logger
from hmsss.utils.myUtil import get_genome_id

log = get_logger(__name__)

"""
HMMER-based search utilities.

Implements parallel execution of HMMER `hmmsearch` across genomes,
progress tracking via shared counters, and conversion of `domtblout`
into `.hmmreport` files with genome-prefixed target IDs for global uniqueness.
"""

# Global shared variables
current_counter = None
counter_lock = None


########################################################################################################
#################### HMMsearch routines for individual file input ######################################
########################################################################################################


def consecutive_hmm_search(config: Config, processes: int = 4) -> dict[str, str]:
    """
    Executes a parallelized HMM search across multiple genome protein files and processes the results.

    This routine performs the following steps for each genome in `options.queued_genomes`:
    1. Runs an HMMER search (HMMsearch) against the specified HMM library using the genome's .faa protein file.
    2. Prefixes each domain hit ID in the resulting domtblout file with the genome-specific filename, making all hits uniquely traceable across genomes.
    3. Tracks and prints the progress of processing using a shared counter across worker processes.

    Parameters:
    - options: A configuration object containing:
        - queued_genomes (list of str): Genome IDs to be processed.
        - faa_files (dict): Mapping from genome IDs to paths of .faa files.
        - library (str): Path to the HMM database.
        - thrs_score (float): Minimum score threshold for accepting hits.
        - clean_reports (bool): Whether to delete old report files before search.
        - redo_search (bool): Whether to force re-running existing searches.
    - processes (int): Number of parallel worker processes to use.

    Returns:
    - dict: Mapping from genome ID to the corresponding '.hmmreport' file path, where domain hit IDs have been prefixed.
    """
    total = len(config.queued_genomes)
    return_hmmreports: dict[str, str] = {}
    # Prepare arguments for Pool processing
    args = []
    for genomeID in config.queued_genomes:
        faa_path = config.faa_files[genomeID]
        hmmreport_path = config.hmmreport_files.get(genomeID)

        # Fill with hmmreport filename
        if hmmreport_path:
            return_hmmreports[genomeID] = hmmreport_path
        else:
            hmmreport_path = os.path.splitext(faa_path)[0] + ".hmmreport"

        # append if report is missing oder overwrite is allowed
        if not os.path.isfile(hmmreport_path) or config.clean_reports:
            args.append(
                (
                    genomeID,
                    faa_path,
                    config.library,
                    hmmreport_path,  # future report
                    config.thrs_score,
                    config.clean_reports,
                    total,
                )
            )

    # Shared Counter und Lock erstellen
    counter = Value("i", 0)  # 'i' = integer
    lock = Lock()

    with Pool(
        processes=processes, initializer=init_globals, initargs=(counter, lock)
    ) as pool:
        results = pool.starmap(run_search, args)

        for genomeID, hmmreport_path in results:
            return_hmmreports[genomeID] = hmmreport_path

    return return_hmmreports


def init_globals(counter: Value, lock: Lock) -> None:
    global current_counter
    global counter_lock
    current_counter = counter
    counter_lock = lock


def run_search(
    genomeID: str,
    faa_file: str,
    query_db: str,
    hmmreport_path: str,
    score: float,
    clean_reports: bool,
    total: int,
) -> str:
    """Runs hmmsearch for a single protein FASTA file and creates a prefixed report.

    This function calls HMMER's hmmsearch command on the given protein FASTA file,
    prefixes each hit with the file's basename, and writes the results to a new .hmmreport file.

    Args:
        hmmreport_path (str): Path to the HMM report file
        faa_file (str): Path to the input protein FASTA file.
        query_db (str): Path to the HMM profile database.
        score (float): The minimum bit score threshold for reporting hits.
        clean_reports (bool): If True, existing reports will be overwritten.
        total (int): Total number of genomes/files to process (for progress logging).

    Returns:
        str: Path to the resulting .hmmreport file.

    Example:
        >>> run_search("/tmp/A.faa", "/tmp/db.hmm", 42.0, True, 12)
        '/tmp/A.hmmreport'
    """

    with counter_lock:
        counter = current_counter.value + 1
        current_counter.value = counter
        print(f"Processing file {counter} of {total}", end="\r")

    domtblout_path = hmm_search(faa_file, query_db, score, clean_reports, 2)
    hmmreport = prefix_domtblout_hits(domtblout_path, hmmreport_path, separator="___")

    return genomeID, hmmreport


def hmm_search(
    faa_path: str,
    query_db: str,
    score: float,
    clean_reports: bool = False,
    cores: int = 1,
) -> str:
    """Executes HMMER hmmsearch and returns path to the domtblout file.

    Args:
        faa_path (str): Path to input protein FASTA file.
        query_db (str): Path to HMM profile database.
        score (float): Minimum bit score threshold for reporting hits.
        clean_reports (bool, optional): If True, overwrite existing output files.
        cores (int, optional): Number of CPU cores to use.

    Returns:
        str: Path to the generated domtblout file.

    Example:
        >>> hmm_search("/tmp/A.faa", "/tmp/db.hmm", 42.0)
        '/tmp/A.domtblout'
    """

    output = os.path.splitext(faa_path)[0] + ".domtblout"
    os.system(
        f"hmmsearch -T {score} --domT {score} --cpu {str(cores)} --noali --domtblout {output} {query_db} {faa_path} > /dev/null 2>&1"
    )
    return output


def prefix_domtblout_hits(
    domtblout_path: str, hmmreport_path: str, separator: str = "___"
) -> str:
    """Prefixes each domain hit ID in a domtblout file with the file's basename and writes to a new file.

    Args:
        domtblout_path (str): Path to the domtblout file to process.
        hmmreport_path (str): Path to the hypothetical hmmreport file
        separator (str, optional): String used between basename and original hit ID.

    Returns:
        str: Path to the new .hmmreport file with prefixed IDs.

    Example:
        >>> prefix_domtblout_hits('/tmp/A.domtblout')
        '/tmp/A.hmmreport'
    """

    # Neuen Prefix vorbereiten (basename ohne Endung)
    basename = os.path.splitext(os.path.basename(domtblout_path))[0]
    basename = get_genome_id(basename)
    if os.path.isfile(domtblout_path):
        with open(domtblout_path, "r") as infile, open(hmmreport_path, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue

                parts = line.strip().split()
                if len(parts) < 5:
                    continue

                parts[0] = f"{basename}{separator}{parts[0]}"
                outfile.write("\t".join(parts) + "\n")
        os.remove(domtblout_path)

    return hmmreport_path

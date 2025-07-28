#!/usr/bin/python

import os
import sys
import glob
import subprocess

from multiprocessing import Pool,Value, Lock
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional


from . import myUtil

logger = myUtil.logger

# Global shared variables
current_counter = None
counter_lock = None


########################################################################################################
#################### HMMsearch routines for individual file input ######################################
########################################################################################################

def unified_search(options: Any, processes: int = 4) -> Dict[str, str]:
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
    total = len(options.queued_genomes)

    args = [
        (
            options.faa_files[genomeID],
            options.library,
            options.thrs_score, # static minimal cutoff score
            options.clean_reports,
            total
        )
        for genomeID in options.queued_genomes
    ]

    # Shared Counter und Lock erstellen
    counter = Value('i', 0)  # 'i' = integer
    lock = Lock()

    with Pool(processes=processes, initializer=init_globals, initargs=(counter, lock)) as pool:
        results = pool.starmap(run_search, args)

    return dict(zip(options.queued_genomes, results))
    
    
    
def init_globals(counter: Value, lock: Lock) -> None:
    global current_counter
    global counter_lock
    current_counter = counter
    counter_lock = lock



def run_search(
    faa_file: str,
    query_db: str,
    score: float,
    clean_reports: bool,
    total: int
) -> str:
    """Runs hmmsearch for a single protein FASTA file and creates a prefixed report.

    This function calls HMMER's hmmsearch command on the given protein FASTA file,
    prefixes each hit with the file's basename, and writes the results to a new .hmmreport file.

    Args:
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
    hmmreport = prefix_domtblout_hits(domtblout_path, separator="___", suffix=".hmmreport")
    
    return hmmreport
    
    
def hmm_search(
    path: str,
    query_db: str,
    score: float,
    clean_reports: bool = False,
    cores: int = 1
) -> str:
    """Executes HMMER hmmsearch and returns path to the domtblout file.

    Args:
        path (str): Path to input protein FASTA file.
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
    
    output = os.path.splitext(path)[0] + '.domtblout'
    hmmreport = os.path.splitext(path)[0] + '.hmmreport'
    if not os.path.isfile(hmmreport) or clean_reports:
        os.system(f'hmmsearch -T {score} --domT {score} --cpu {str(cores)} --noali --domtblout {output} {query_db} {path} > /dev/null 2>&1')
    return output


def prefix_domtblout_hits(
    domtblout_path: str,
    separator: str = "___",
    suffix: str = ".hmmreport"
) -> str:
    """Prefixes each domain hit ID in a domtblout file with the file's basename and writes to a new file.

    Args:
        domtblout_path (str): Path to the domtblout file to process.
        separator (str, optional): String used between basename and original hit ID.
        suffix (str, optional): Suffix for the output file.

    Returns:
        str: Path to the new .hmmreport file with prefixed IDs.

    Example:
        >>> prefix_domtblout_hits('/tmp/A.domtblout')
        '/tmp/A.hmmreport'
    """

    # Neuen Prefix vorbereiten (basename ohne Endung)
    basename = os.path.splitext(os.path.basename(domtblout_path))[0]

    output_path = os.path.splitext(domtblout_path)[0] + suffix
    if os.path.isfile(domtblout_path):
        with open(domtblout_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    continue

                parts = line.strip().split()
                if len(parts) < 5:
                    continue

                parts[0] = f"{basename}{separator}{parts[0]}"
                outfile.write('\t'.join(parts) + '\n')
        os.remove(domtblout_path)
        
    return output_path


def concatenate_hmmreports_cat(
    report_paths: Dict[str, str],
    output_path: str = "global_report.cat_hmmreport"
) -> str:
    """Concatenates multiple .hmmreport files into a single global file.

    Args:
        report_paths (Dict[str, str]): Mapping of genome ID to .hmmreport file path.
        output_path (str, optional): Path for the concatenated output file.

    Returns:
        str: Path to the combined global report file.

    Example:
        >>> concatenate_hmmreports_cat({'g1': '/tmp/g1.hmmreport', 'g2': '/tmp/g2.hmmreport'}, '/tmp/all.cat_hmmreport')
        '/tmp/all.cat_hmmreport'
    """
    
    logger.info(f"Concatenate hit reports to {output_path}")
    
    # Filtere nur existierende Dateien
    valid_paths = [path for path in report_paths.values() if os.path.isfile(path)]
    if not valid_paths:
        logger.error(f"No valid hmmreport files found in {report_paths}")
        raise FileNotFoundError(f"No valid hmmreport files found in directory {report_paths}.")

    # Führe das cat-Kommando aus
    cmd = ["cat"] + valid_paths
    with open(output_path, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile)

    return output_path
    
########################################################################################################
#################### Filter the glob report to trusted hits and potential hits #########################
########################################################################################################

def filter_trusted_and_noise_hits(
    options: object,
    glob_report: str,
    processes: int = 4
) -> str:
    """Filters hits in a global report into trusted/intermediate categories and writes summary files.

    Uses score thresholds defined in options.score_threshold_file.

    Args:
        options (object): Configuration object containing:
            - score_threshold_file (str): Path to tab-separated file with thresholds.
            - threshold_type (int): Which column to use for score thresholds.
            - thrs_score (float): Default fallback threshold.
            - Cross_check_directory (str): Directory for output.
        glob_report (str): Path to the concatenated global report file.
        processes (int, optional): Number of worker processes.

    Returns:
        str: Absolute path to the output directory containing filtered hit files.

    Example:
        >>> filter_trusted_and_noise_hits(options, '/tmp/global.cat_hmmreport', 2)
        '/tmp/xcheck'
    """
    trusted_dict = make_threshold_dict(
        options.score_threshold_file, options.threshold_type, options.thrs_score
    )
    noise_dict = make_threshold_dict(
        options.score_threshold_file, 3, options.thrs_score
    )
    
    output_dir = options.Cross_check_directory
    os.makedirs(output_dir, exist_ok=True)

    args = [
        (
            hmm_id,
            glob_report,
            trusted_dict[hmm_id],
            noise_dict.get(hmm_id, options.thrs_score),
            output_dir
        )
        for hmm_id in trusted_dict
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(process_single_hmm, args)

    return os.path.abspath(output_dir)
    
def make_threshold_dict(
    file_path: str,
    threshold_type: int = 1,
    default_score: float = 50.0
) -> Dict[str, float]:
    """Parses a tab-separated score threshold file into a dict.

    Args:
        file_path (str): Path to thresholds TSV file.
        threshold_type (int, optional): Which column (0-based) to use for score.
        default_score (float, optional): Default fallback score if missing.

    Returns:
        Dict[str, float]: Mapping from HMM ID to score threshold.

    Example:
        >>> make_threshold_dict('/tmp/thresholds.tsv', 1)
        {'PF00001': 42.0, 'PF00002': 50.0}
    """
    thresholds: Dict[str, float] = {}
    with open(file_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            parts = line.strip().split("\t")
            key = parts[0] if parts else None
            score = default_score

            try:
                if len(parts) == 2:
                    # Only one score present, use it
                    score = float(parts[1])
                elif len(parts) > threshold_type:
                    score = float(parts[threshold_type])
            except (ValueError, IndexError) as e:
                logger.warning(f"[Line {line_number}] Problem parsing: {line.strip()} — {e}")
                continue

            if key:
                thresholds[key] = score

    return thresholds


def process_single_hmm(
    hmm_id: str,
    glob_report: str,
    trusted_cutoff: float,
    noise_cutoff: float,
    output_dir: str
) -> None:
    """Processes hits for one HMM ID, separating trusted and intermediate hits.

    Writes two files: {hmm_id}.trusted_hits and {hmm_id}.intermediate_hits.

    Args:
        hmm_id (str): HMM profile ID to extract.
        glob_report (str): Path to the global cat_hmmreport file.
        trusted_cutoff (float): Score threshold for trusted hits.
        noise_cutoff (float): Lower bound for intermediate hits.
        output_dir (str): Directory to write result files.

    Returns:
        None

    Example:
        >>> process_single_hmm('PF00001', '/tmp/global.cat_hmmreport', 42.0, 21.0, '/tmp/xcheck')
    """
    
    # Define output list for trusted and intermediate hits
    trusted_path = os.path.join(output_dir, f"{hmm_id}.trusted_hits")
    intermediate_path = os.path.join(output_dir, f"{hmm_id}.intermediate_hits")
    
    # Check if the trusted and the intermediate hit files are already present and skip existing files
    if os.path.isfile(intermediate_path):
        return
    if os.path.isfile(trusted_path):
        return
    
    
    # trusted hits + collect candidates
    candidates = {}
    with open(glob_report, 'r') as infile, open(trusted_path, 'w') as trusted_out:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')

            target = parts[0]
            hit_hmm = parts[3]

            if hit_hmm != hmm_id:
                continue

            try:
                score = float(parts[7])  # Bit-Score
                if int(parts[10]) > 1:
                # Domain Bit score is in column 13
                    score = float(parts[13])
                else:
                # Full Bit score is in column 7 (0-based index)
                    score = float(parts[7])

            except ValueError:
                continue


            if score >= trusted_cutoff:
                trusted_out.write(line)
            elif score > noise_cutoff:
                candidates[target] = score
    
    # remove empty trusted files 
    if os.path.getsize(trusted_path) == 0:
        os.remove(trusted_path)
    # and skip if no candidates found
    if not candidates:
        return

    # Filter out candidates with better hits in other HMMs
    with open(glob_report, 'r') as infile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')

            target = parts[0]
            hit_hmm = parts[3]
            try:
                score = float(parts[2])
            except ValueError:
                continue

            if target in candidates and hit_hmm != hmm_id and score > candidates[target]:
                del candidates[target]

    # If not candidates left leave the routine
    if not candidates:
        return

    # Write down remaining candidates
    with open(glob_report, 'r') as infile, open(intermediate_path, 'w') as interm_out:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')

            target = parts[0]
            hit_hmm = parts[3]

            if hit_hmm == hmm_id and target in candidates:
                interm_out.write(line)
    
    if os.path.getsize(intermediate_path) == 0:
        os.remove(intermediate_path)
    return






def extract_fasta_per_intermediate_hitfile(
    options: object,
    intermediate_hit_dir: str
) -> None:
    """Extracts protein sequences for every .intermediate_hits file and writes .intermediate_hit_faa files.

    Args:
        options (object): Configuration object containing:
            - faa_files (Dict[str, str]): Mapping from genome ID to FASTA path.
        intermediate_hit_dir (str): Directory containing .intermediate_hits files.

    Returns:
        None

    Example:
        >>> extract_fasta_per_intermediate_hitfile(options, '/tmp/xcheck')
    """
    
    for file in os.listdir(intermediate_hit_dir):
        if not file.endswith(".intermediate_hits"):
            continue

        hmm_id = file.replace(".intermediate_hits", "")
        hitfile_path = os.path.join(intermediate_hit_dir, file)
        output_fasta = os.path.join(intermediate_hit_dir, f"{hmm_id}.intermediate_hit_faa")

        # IDs sammeln: genomeID → set(proteinIDs)
        genome_hits = {}
        with open(hitfile_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                full_id = parts[0]
                if "___" not in full_id:
                    continue
                genome_id, protein_id = full_id.split("___", 1)
                genome_hits.setdefault(genome_id, set()).add(protein_id)

        # Write sequential files per genome
        with open(output_fasta, 'w') as out:
            for genome_id, protein_ids in genome_hits.items():
                faa_path = options.faa_files.get(genome_id)
                if not faa_path or not os.path.isfile(faa_path):
                    logger.warning(f".faa file not found for genome: {genome_id}")
                    continue

                with open(faa_path, 'r') as faa:
                    write = False
                    header_id = None

                    for line in faa:
                        if line.startswith(">"):
                            header_id = line[1:].split()[0]
                            write = header_id in protein_ids
                            if write:
                                out.write(f">{genome_id}___{header_id}\n")
                        elif write:
                            out.write(line)

        print(f"[✓] {hmm_id} → {output_fasta}")



def process_hitfile(
    hitfile_path: str,
    intermediate_hit_dir: str,
    faa_files: Dict[str, str]
) -> None:
    """Extracts relevant FASTA sequences for a single intermediate hitfile.

    Args:
        hitfile_path (str): Path to the .intermediate_hits file.
        intermediate_hit_dir (str): Directory to write the output .fasta file.
        faa_files (Dict[str, str]): Mapping from genome ID to FASTA file path.

    Returns:
        None

    Example:
        >>> process_hitfile('/tmp/xcheck/PF00001.intermediate_hits', '/tmp/xcheck', {'g1': '/tmp/g1.faa'})
    """
    
    hmm_id = os.path.basename(hitfile_path).replace(".intermediate_hits", "")
    output_fasta = os.path.join(intermediate_hit_dir, f"{hmm_id}.intermediate_hits_faa")

    genome_hits = {}
    with open(hitfile_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            full_id = parts[0]
            if "___" not in full_id:
                continue
            genome_id, protein_id = full_id.split("___", 1)
            genome_hits.setdefault(genome_id, set()).add(protein_id)

    with open(output_fasta, 'w') as out:
        for genome_id, protein_ids in genome_hits.items():
            faa_path = faa_files.get(genome_id)
            if not faa_path or not os.path.isfile(faa_path):
                logger.warning(f"FASTA not found for {genome_id}")
                continue

            with open(faa_path, 'r') as faa:
                write = False
                header_id = None

                for line in faa:
                    if line.startswith(">"):
                        header_id = line[1:].split()[0]
                        write = header_id in protein_ids
                        if write:
                            out.write(f">{genome_id}___{header_id}\n")
                    elif write:
                        out.write(line)


def generate_faa_per_hitfile_parallel(
    options: object,
    intermediate_hit_dir: str,
    processes: int = 4
) -> None:
    """Parallel extraction of FASTA for all .intermediate_hits files in a directory.

    Args:
        options (object): Configuration object containing:
            - faa_files (Dict[str, str]): Mapping from genome ID to FASTA path.
        intermediate_hit_dir (str): Directory containing .intermediate_hits files.
        processes (int, optional): Number of parallel workers.

    Returns:
        None

    Example:
        >>> generate_faa_per_hitfile_parallel(options, '/tmp/xcheck', 2)
    """
    
    output_dir = intermediate_hit_dir  # same dir for output
    faa_files = options.faa_files      # dict: genome_id → path

    hitfiles = [
        os.path.join(intermediate_hit_dir, f)
        for f in os.listdir(intermediate_hit_dir)
        if f.endswith(".intermediate_hits")
    ]

    args = [(hitfile, output_dir, faa_files) for hitfile in hitfiles]

    with Pool(processes=processes) as pool:
        pool.starmap(process_hitfile, args)

##########################################################################################################################################################
#################### Cross check hits with reference sequences and add the hmmreport lines to the trusted cutoff hmmreport ###############################
##########################################################################################################################################################

def find_refseq_file(base_dir: str, filename: str) -> Optional[str]:
    """Recursively search a base directory for a file with a specific name.

    Args:
        base_dir (str): Directory to search in.
        filename (str): Name of the file to find.

    Returns:
        Optional[str]: Full path to the found file, or None if not found.

    Example:
        >>> find_refseq_file('/data', 'ref.fa')
        '/data/refs/ref.fa'
    """
    
    for root, _, files in os.walk(base_dir):
        if filename in files:
            return os.path.join(root, filename)
    return None


def cross_check_candidates_with_reference_seqs(options) -> List[str]:
    """Cross-checks candidate hit sequences with reference sequences using DIAMOND.

    For each .intermediate_hits_faa file in options.Cross_check_directory, 
    the function searches for a reference database in the RefSeqs directory. 
    If no .dmnd database is present, it tries to create one from a .faa file.
    Then, it runs a DIAMOND search to compare the candidates to the reference sequences.

    Args:
        options (object): Configuration object containing:
            - execute_location (str): Root execution directory.
            - Cross_check_directory (str): Directory containing .intermediate_hits_faa files.
            - refseq_identity (float): Percent identity cutoff for DIAMOND.
            - cores (int): Number of threads for DIAMOND.

    Returns:
        List[str]: List of HMM IDs where no reference sequence or db could be found.

    Example:
        >>> missing = cross_check_candidates_with_reference_seqs(options)
    """
    logger.info("Cross check hit sequences with reference sequences")

    refseq_dir = os.path.join(options.execute_location, "src", "RefSeqs")
    refseq_unavailable_list = []
    
    cross_check_dir = options.Cross_check_directory # directoy with the intermediate hit fasta faa files
    intermediate_files = glob.glob(os.path.join(cross_check_dir, "*.intermediate_hits_faa"))

    # Iterate the intermediate faa files
    for inter_file in intermediate_files:
        logger.debug(f"Checking reference sequences for candidates sequences in {inter_file}")
        hmm_id = os.path.splitext(os.path.basename(inter_file))[0].replace(".intermediate_hits_faa", "")
        db_file = f"{hmm_id}.dmnd"
        db_path = find_file_in_prefixed_subdirs(refseq_dir, db_file, dir_prefix="") #dir_prefix is for version control, possibly uneccessary

        # Prüfen ob .dmnd existiert, sonst erstellen
        if not os.path.isfile(db_path):
            faa_file = f"{hmm_id}.faa"
            faa_path = find_file_in_prefixed_subdirs(refseq_dir, faa_file, dir_prefix="")

            if os.path.isfile(faa_path):
                logger.debug(f"Creating Diamond DB from {faa_path} because {db_path} was not found")

                try:
                    subprocess.run(["diamond", "makedb", "--in", faa_path, "-d", db_base], check=True)
                    db_path = db_base + ".dmnd"
                except subprocess.CalledProcessError:
                    logger.error(f"Failed to create Diamond database for {faa_path}")
                    refseq_unavailable_list.append(hmm_id)
                    continue
            else:
                logger.warning(f"Skipping {hmm_id}: Reference sequence file not found.")
                refseq_unavailable_list.append(hmm_id)
                continue

        output_file = os.path.join(cross_check_dir, f"{hmm_id}.crosschecked.tsv")
        diamond = myUtil.find_executable("diamond")
        cmd = [
            diamond, "blastp",
            "--query", inter_file,
            "--db", db_path,
            "--out", output_file,
            "--outfmt", "6",
            "--max-target-seqs", "1",
            "--id", str(options.refseq_identity),
            "--threads", str(options.cores),
            "--quiet"
        ]

        logger.info(f"Verifying {hmm_id} hits with reference sequences")
        result = subprocess.run(cmd)

        if result.returncode != 0:
            logger.error(f"DIAMOND search failed for {hmm_id}")
            continue
        if os.path.getsize(output_file) == 0:
            os.remove(output_file)

    return refseq_unavailable_list
    
    
def find_file_in_prefixed_subdirs(base_dir: str, filename: str, dir_prefix: str) -> str:
    """Recursively searches for a file in subdirectories with a specific prefix.

    Args:
        base_dir (str): Root directory to search.
        filename (str): Filename to find.
        dir_prefix (str): Subdirectory name prefix to restrict the search (empty string for all).

    Returns:
        str: Full path to the found file, or empty string if not found.

    Example:
        >>> find_file_in_prefixed_subdirs('/data/RefSeqs', 'PF00001.dmnd', '')
        '/data/RefSeqs/PF00001/PF00001.dmnd'
    """
    for root, dirs, files in os.walk(base_dir):
        # Nur Verzeichnisse mit dem gewünschten Prefix betreten
        if not os.path.basename(root).startswith(dir_prefix):
            continue

        if filename in files:
            return os.path.join(root, filename)

    return ""  # nicht gefunden


def process_crosscheck(hmm_id: str, crosscheck_dir: str) -> None:
    """Promotes intermediate hits to trusted hits if they are validated by crosschecking.

    Args:
        hmm_id (str): HMM ID to process.
        crosscheck_dir (str): Directory containing .crosschecked.tsv, .intermediate_hits, .trusted_hits files.

    Returns:
        None

    Example:
        >>> process_crosscheck('PF00001', '/tmp/xcheck')
    """
    crosscheck_path = os.path.join(crosscheck_dir, f"{hmm_id}.crosschecked.tsv")
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    if not os.path.exists(crosscheck_path):
        logger.error(f"Crosscheck file missing: '{crosscheck_path}'")
    if not os.path.exists(intermediate_path):
        logger.error(f"Intermediate file missing: '{intermediate_path}'")#
        
    # Lade IDs aus crosscheck
    valid_hits = set()
    with open(crosscheck_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            valid_hits.add(line.strip().split('\t')[0])

    if not valid_hits:
        return

    promoted_count = 0

    # Promote direkt
    with open(intermediate_path, 'r') as interm, open(trusted_path, 'a') as trusted:
        for line in interm:
            if line.startswith("#") or not line.strip():
                continue
            if line.strip().split('\t')[0] in valid_hits:
                trusted.write(line)
                promoted_count += 1

    logger.info(f"{hmm_id}: promoted {promoted_count} hits via comparison with reference sequences")


def promote_crosschecked_hits(crosscheck_dir: str, processes: int = 4) -> None:
    """Parallelizes promoting of hits from .intermediate_hits to .trusted_hits
    using the .crosschecked.tsv files in a directory.

    Args:
        crosscheck_dir (str): Directory with crosscheck and hit files.
        processes (int): Number of worker processes.

    Returns:
        None

    Example:
        >>> promote_crosschecked_hits('/tmp/xcheck', 4)
    """
    
    hmm_ids = [
        f.replace(".crosschecked.tsv", "")
        for f in os.listdir(crosscheck_dir)
        if f.endswith(".crosschecked.tsv")
    ]

    with Pool(processes=processes) as pool:
        pool.starmap(
            process_crosscheck,
            [(hmm_id, crosscheck_dir) for hmm_id in hmm_ids]
        )
    
    
def summarize_trusted_hits(directory: str, crosscheck_dir: str, name: str, suffix: str = ".trusted_hits") -> str:
    """Creates a summary file by concatenating all .trusted_hits files in a directory.

    Args:
        directory (str): Directory to write the summary file to.
        crosscheck_dir (str): Directory with .trusted_hits files.
        name (str): Name of the summary file.
        suffix (str, optional): File suffix to match (default '.trusted_hits').

    Returns:
        str: Path to the summary file.

    Example:
        >>> summarize_trusted_hits('/tmp/results', '/tmp/xcheck', 'summary.txt')
        '/tmp/results/summary.txt'
    """
    
    summary_path = os.path.join(directory, name)
    
    if os.path.isfile(summary_path) and os.path.getsize(summary_path) != 0:
    	return summary_path
    
    trusted_files = [
        os.path.join(crosscheck_dir, f)
        for f in os.listdir(crosscheck_dir)
        if f.endswith(suffix)
    ]

    if not trusted_files:
        logger.error("No trusted hit files found for summary.")
        logger.info("There were no hits found in any genome. Closing the search")
        sys.exit()
        
    exit_code = os.system(f"cat {' '.join(trusted_files)} > {summary_path}")
    
    return summary_path
    
    

    
####################################################################################################
#################### Promote by cutoff when cross check is not available ###########################
####################################################################################################

def promote_by_cutoff(
    options: object,
    directory: str,
    processes: int = 4,
    hmm_ids: Optional[List[str]] = None
) -> None:
    """Parallel promotion of intermediate hits to trusted hits based on thresholds.

    Args:
        options (object): Configuration with threshold file and type.
        directory (str): Directory with .intermediate_hits files.
        processes (int, optional): Number of worker processes.
        hmm_ids (List[str] or None): List of HMM IDs, or 'all' for all present.

    Returns:
        None

    Example:
        >>> promote_by_cutoff(options, '/tmp/xcheck', 4, hmm_ids=['PF00001', 'PF00002'])
        >>> promote_by_cutoff(options, '/tmp/xcheck', 4, hmm_ids="all")
    """
    # If no hmm identifier were defined use all that are in
    if hmm_ids=="all":
        hmm_ids = [
            f.replace(".intermediate_hits", "")
            for f in os.listdir(directory)
            if f.endswith(".intermediate_hits")
        ]

    optimized_dict = make_threshold_dict(
        options.score_threshold_file, options.threshold_type, options.thrs_score
    )
    
    with Pool(processes=processes) as pool:
        pool.starmap(
            process_optimized_cutoff,
            [(hmm_id, directory, optimized_dict) for hmm_id in hmm_ids]
        )




def process_optimized_cutoff(hmm_id: str, crosscheck_dir: str, optimized_dict: Dict[str, float]) -> None:
    """Promotes intermediate hits to trusted hits based on a score threshold.

    Args:
        hmm_id (str): HMM ID to process.
        crosscheck_dir (str): Directory with hit files.
        optimized_dict (Dict[str, float]): Dict of {hmm_id: threshold_score}.

    Returns:
        None

    Example:
        >>> process_optimized_cutoff('PF00001', '/tmp/xcheck', {'PF00001': 42.0})
    """
    
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")
    
    threshold_score = optimized_dict.get(hmm_id, 50)
    
    if not os.path.isfile(intermediate_path):
        logger.warning(f"Intermediate or trusted hit file missing for {hmm_id}")
        return

    promoted_count = 0

    # Promote direkt
    with open(intermediate_path, 'r') as interm, open(trusted_path, 'a') as trusted:
        for line in interm:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            score = float(parts[7])
            if score >= threshold_score:
                trusted.write(line)
                promoted_count += 1

    logger.info(f"{hmm_id}: promoted {promoted_count} hits due to the given threshold {threshold_score}")

    




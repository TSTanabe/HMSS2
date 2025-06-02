#!/usr/bin/python
import os
import subprocess
from multiprocessing import Pool,Value, Lock
from collections import defaultdict

from typing import Any, Dict, List, Optional, Tuple, Set


from . import myUtil
import glob

# Global shared variables
current_counter = None
counter_lock = None


########################################################################################################
#################### HMMsearch routines for individual file input ######################################
########################################################################################################

def unified_search(options: Any, processes: int = 4) -> Dict[str, str]:
    """
    Executes a parallelized HMM search across multiple genome protein files and processes the results.

    This function performs the following steps for each genome in options.queued_genomes:
      1. Runs an HMMER search (run_search) against the specified HMM library using the genome's .faa protein file.
      2. Prefixes each domain hit ID in the resulting report with the genome-specific filename,
         making all hits uniquely traceable across genomes.
      3. Tracks and prints the progress of processing using a shared counter across worker processes.

    Args:
        options: Configuration object with at least these attributes:
            - queued_genomes (List[str]): Genome IDs to be processed.
            - faa_files (Dict[str, str]): Mapping from genome ID to path of .faa file.
            - library (str): Path to the HMM database file.
            - thrs_score (float): Minimum score threshold for accepting hits.
            - clean_reports (bool): Whether to delete old report files before search.
            - redo_search (bool): Whether to force re-running existing searches.
        processes: Number of parallel worker processes to use (default: 4).

    Returns:
        Dict[str, str]: Mapping from genome ID to the corresponding '.hmmreport' file path,
                        where domain hit IDs have been prefixed.
    """
    total = len(options.queued_genomes)

    args = [
        (
            options.faa_files[genomeID],
            options.library,
            options.thrs_score,
            options.clean_reports,
            options.redo_search,
            total
        )
        for genomeID in options.queued_genomes
    ]

    # Shared Counter und Lock erstellen
    counter = Value('i', 0)  # 'i' = integer
    lock = Lock()

    with Pool(processes=processes, initializer=init_globals, initargs=(counter, lock)) as pool:
        results = pool.starmap(run_search, args)
    
    print("\n")
    return dict(zip(options.queued_genomes, results))
    
    
    
def init_globals(counter, lock):
    global current_counter
    global counter_lock
    current_counter = counter
    counter_lock = lock



def run_search(
    faa_file: str,
    query_db: str,
    score: float,
    clean_reports: bool,
    redo_search: bool,
    total: int
) -> str:
    """
    Executes an HMM search for a single genome and prefixes each domain hit ID.

    This function updates a shared counter to report progress, then:
      1. Runs HMMsearch on the given .faa file against the query database.
      2. Prefixes each domain hit in the resulting domtblout file.
      3. Returns the path to the final .hmmreport file.

    Args:
        faa_file: Path to the genome's .faa protein file.
        query_db: Path to the HMM database file.
        score: Minimum score threshold for accepting hits.
        clean_reports: If True, delete any existing reports before searching.
        redo_search: If True, force re-running the search even if a report exists.
        total: Total number of genomes being processed (for progress reporting).

    Returns:
        Path to the generated .hmmreport file.

    Raises:
        RuntimeError: If HMMsearch or prefixing fails.
    """
    # Update and display progress
    with counter_lock:
        current = current_counter.value + 1
        current_counter.value = current
        print(f"[INFO] Processing file {current} of {total}", end="\r", flush=True)

    # Perform the HMM search
    domtblout_path = hmm_search(
        faa_file,
        query_db,
        score,
        clean_reports,
        redo_search,
        2  # number of threads within each process
    )
    if not domtblout_path or not isinstance(domtblout_path, str):
        raise RuntimeError(f"[ERROR] HMMsearch failed for file: {faa_file}")

    # Prefix domain hit IDs and produce final .hmmreport
    hmmreport_path = prefix_domtblout_hits(
        domtblout_path,
        separator="___",
        suffix=".hmmreport"
    )
    if not hmmreport_path or not isinstance(hmmreport_path, str):
        raise RuntimeError(f"[ERROR] Prefixing domtblout hits failed for file: {domtblout_path}")

    return hmmreport_path

    
    
def hmm_search(
    fasta_path: str,
    query_db: str,
    threshold: float,
    clean_reports: bool = False,
    redo_search: bool = False,
    cores: int = 1
) -> Optional[str]:

    output = os.path.splitext(fasta_path)[0] + '.domtblout'
    hmmreport = os.path.splitext(fasta_path)[0] + '.hmmreport'
    if not os.path.isfile(hmmreport) or clean_reports or redo_search:
        os.system(
            f"hmmsearch -T {threshold} --domT {threshold} --cpu {cores} "
            f"--noali --domtblout {output} {query_db} {fasta_path} > /dev/null 2>&1"
        )
    return output


def prefix_domtblout_hits(
    domtblout_path: str,
    separator: str = "___",
    suffix: str = ".hmmreport"
) -> Optional[str]:
    """
    Reads a domtblout file, prefixes each hit ID (column 1) with the file's basename,
    and writes the modified lines to a new file with the given suffix. Removes the
    original domtblout file afterward.

    Args:
        domtblout_path: Path to the original domtblout file.
        separator: String inserted between the basename and the original hit ID.
        suffix:   Suffix for the new output file (replacing the original extension).

    Returns:
        The path to the newly created .hmmreport file if successful; otherwise None.

    Raises:
        OSError: If reading or writing files fails.
    """
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


def concatenate_hmmreports(
    report_paths: Dict[str, str],
    output_path: str = "global_report.cat_hmmreport"
) -> str:
    """
    Concatenates multiple .hmmreport files into a single global report using the
    UNIX 'cat' command.

    Args:
        report_paths: A mapping from genome ID (or other key) to the path of its
                      .hmmreport file.
        output_path:  Path to the output concatenated file. Defaults to
                      "global_report.cat_hmmreport".

    Returns:
        The path to the concatenated global report (output_path).

    Raises:
        FileNotFoundError: If none of the provided paths exist as files.
        RuntimeError:    If the 'cat' command fails.
    """
    print("[INFO] Concatenating HMM report files")

    # Collect only the paths that actually exist on disk
    valid_paths = [p for p in report_paths.values() if os.path.isfile(p)]
    if not valid_paths:
        raise FileNotFoundError(
            "[ERROR] No valid .hmmreport files found among the provided paths."
        )

    # Execute the 'cat' command to merge them
    try:
        with open(output_path, "w") as outfile:
            subprocess.run(["cat"] + valid_paths, stdout=outfile, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"[ERROR] Failed to concatenate reports: {exc}") from exc

    return output_path
    
########################################################################################################
#################### Filter the glob report to trusted hits and potential hits #########################
########################################################################################################
    
    
def make_threshold_dict(
    file_path: str,
    threshold_type: int = 1,
    default_score: float = 50.0
) -> Dict[str, float]:
    """
    Builds a dictionary of thresholds from a tab-separated file.

    Each line in the file should have at least a key (first column). 
    Additional numeric columns represent different score types.

    Rules:
      - If a line contains exactly two columns, use the second column as the score.
      - If a line has more columns than 'threshold_type', use the column at that index.
      - If no valid score is found, use the default_score.

    Args:
        file_path: Path to the tab-separated thresholds file.
        threshold_type: Index of the score column to use when multiple scores are present.
        default_score: Score to use if parsing fails or no score columns exist.

    Returns:
        A dictionary mapping each key string to its chosen float score.
    """
    thresholds = {}
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
                print(f"[ERROR] Cannot parse [Line {line_number}] {line.strip()} — {e}")
                continue

            if key:
                thresholds[key] = score

    return thresholds

def filter_trusted_and_noise_hits(options: Any, processes: int = 4) -> str:
    """
    Filters HMM hits from a global report into 'trusted' and 'noise' categories
    based on threshold scores, in parallel across multiple HMM IDs.

    This function:
      1. Reads trusted and noise thresholds from a tab-separated threshold file.
      2. Ensures the Cross_check_directory exists.
      3. For each HMM ID in the trusted threshold dictionary, calls process_single_hmm
         in parallel, passing the trusted and noise thresholds along with the global report.
      4. Returns the absolute path of the Cross_check_directory.

    Args:
        options: An object with at least the following attributes:
            - glob_report (str): Path to the concatenated HMM report file.
            - score_threshold_file (str): Path to the tab-separated threshold file.
            - threshold_type (int): Column index for trusted thresholds.
            - thrs_score (float): Default score if threshold lookup fails.
            - Cross_check_directory (str): Directory where filtered outputs will be written.
        processes: Number of parallel worker processes to use (default: 4).

    Returns:
        The absolute path to the Cross_check_directory (str).

    Raises:
        FileNotFoundError: If the glob_report or threshold file cannot be found.
        OSError: If the Cross_check_directory cannot be created.
        KeyError: If a trusted HMM ID is missing from the threshold dictionaries.
    """
    # Validate that the global report exists
    glob_report_path: str = options.glob_report
    if not os.path.isfile(glob_report_path):
        raise FileNotFoundError(f"Global report not found: {glob_report_path}")

    # Build threshold dictionaries
    trusted_dict: Dict[str, float] = make_threshold_dict(
        options.score_threshold_file,
        options.threshold_type, # Here the 'trusted' score is defined by default
        options.thrs_score
    )
    noise_dict: Dict[str, float] = make_threshold_dict(
        options.score_threshold_file,
        threshold_type=3,
        default_score=options.thrs_score
    )

    # Ensure the Cross_check_directory exists
    output_dir: str = options.cross_check_directory
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Cannot create directory '{output_dir}': {exc}") from exc

    # Prepare argument tuples for each HMM ID
    args_list = []
    for hmm_id, trusted_score in trusted_dict.items():
        noise_score = noise_dict.get(hmm_id, options.thrs_score)
        args_list.append(
            (hmm_id, glob_report_path, trusted_score, noise_score, output_dir)
        )

    # Run process_single_hmm() in parallel
    with Pool(processes=processes) as pool:
        pool.starmap(process_single_hmm, args_list)

    return os.path.abspath(output_dir)
    
    
def process_single_hmm(
    hmm_id: str,
    glob_report: str,
    trusted_cutoff: float,
    noise_cutoff: float,
    output_dir: str
) -> None:
    """
    Processes a single HMM ID by filtering hits from a global report into
    trusted and intermediate categories.

    Steps:
      1. Create output paths for trusted and intermediate hit files.
      2. If either file already exists, skip processing.
      3. Read the global report and write lines with scores >= trusted_cutoff
         to the trusted file. Collect lines with scores > noise_cutoff into
         a candidates dict.
      4. If the trusted file ends up empty, delete it.
      5. If no candidates remain, exit.
      6. Re-scan the report to remove any candidate that has a better hit in
         another HMM.
      7. If no candidates remain after filtering, exit.
      8. Write the remaining candidate lines (for this HMM) to the intermediate
         file. If it ends up empty, delete it.

    Args:
        hmm_id:          Identifier of the HMM to process.
        glob_report:     Path to the concatenated HMM report file.
        trusted_cutoff:  Score threshold to classify a hit as trusted.
        noise_cutoff:    Score threshold to classify a hit as intermediate candidate.
        output_dir:      Directory where output files will be written.

    Returns:
        None
    """
    trusted_path = os.path.join(output_dir, f"{hmm_id}.trusted_hits")
    intermediate_path = os.path.join(output_dir, f"{hmm_id}.intermediate_hits")

    # Skip if outputs already exist
    if os.path.isfile(trusted_path) or os.path.isfile(intermediate_path):
        return

    candidates: Dict[str, float] = {}

    # Step 1: Identify trusted hits and initial candidates
    try:
        with open(glob_report, "r") as infile, open(trusted_path, "w") as trusted_out:
            for raw_line in infile:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.strip().split("\t")
                if len(parts) < 8:
                    continue

                target = parts[0]
                hit_hmm = parts[3]

                if hit_hmm != hmm_id:
                    continue

                try:
                    # Bit score is in column 7 (0-based index)
                    score = float(parts[7])
                except ValueError:
                    continue

                if score >= trusted_cutoff:
                    trusted_out.write(raw_line)
                elif score > noise_cutoff:
                    candidates[target] = score
    except OSError:
        # If file operations fail, bail out
        return

    # Remove trusted file if empty
    if os.path.isfile(trusted_path) and os.path.getsize(trusted_path) == 0:
        os.remove(trusted_path)

    # If no candidates, nothing more to do
    if not candidates:
        return

    # Step 2: Filter out candidates with better hits in other HMMs
    try:
        with open(glob_report, "r") as infile:
            for raw_line in infile:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.strip().split("\t")
                if len(parts) < 3:
                    continue

                target = parts[0]
                hit_hmm = parts[3]

                try:
                    score = float(parts[2])
                except ValueError:
                    continue

                # If the same target appears with a different HMM and a higher score,
                # remove it from candidates
                if target in candidates and hit_hmm != hmm_id and score > candidates[target]:
                    candidates.pop(target, None)
    except OSError:
        return

    # If no candidates remain after filtering, exit
    if not candidates:
        return

    # Step 3: Write remaining candidates to intermediate file
    try:
        with open(glob_report, "r") as infile, open(intermediate_path, "w") as interm_out:
            for raw_line in infile:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.strip().split("\t")
                if len(parts) < 4:
                    continue

                target = parts[0]
                hit_hmm = parts[3]

                if hit_hmm == hmm_id and target in candidates:
                    interm_out.write(raw_line)
    except OSError:
        return

    # Remove intermediate file if empty
    if os.path.isfile(intermediate_path) and os.path.getsize(intermediate_path) == 0:
        os.remove(intermediate_path)

    return


def process_hitfile(
    hitfile_path: str,
    intermediate_hit_dir: str,
    faa_files: Dict[str, str]
) -> None:
    """
    Processes a single '.intermediate_hits' file by extracting matching protein sequences
    from corresponding .faa files and writing them into a new FASTA file.

    Args:
        hitfile_path: Path to the '.intermediate_hits' file.
        intermediate_hit_dir: Directory where the output FASTA should be written.
        faa_files: Mapping from genome ID (str) to its .faa file path (str).

    Operation:
      1. Derive the HMM ID from the hitfile name.
      2. Read the hitfile, collecting for each genome ID the set of protein IDs.
      3. For each genome ID, open the corresponding .faa file and write only those
         sequences whose headers match the collected protein IDs.
      4. Prefix sequence headers in the output FASTA with 'genomeID___proteinID'.
    """
    # Derive HMM ID and output FASTA path
    base = os.path.basename(hitfile_path)
    suffix = ".intermediate_hits"
    if base.endswith(suffix):
        hmm_id = base[:-len(suffix)]
    else:
        hmm_id = base
    output_fasta = os.path.join(intermediate_hit_dir, f"{hmm_id}.intermediate_hits_faa")

    # Step 1: Parse hitfile to collect protein IDs per genome
    genome_hits: Dict[str, Set[str]] = {}
    try:
        with open(hitfile_path, "r") as hit_in:
            for raw_line in hit_in:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.strip().split("\t")
                full_id = parts[0]
                if "___" not in full_id:
                    continue

                genome_id, protein_id = full_id.split("___", 1)
                genome_hits.setdefault(genome_id, set()).add(protein_id)
    except OSError as exc:
        print(f"[ERROR] Cannot read hitfile '{hitfile_path}': {exc}")
        return

    # Step 2: Write matching sequences to the output FASTA
    try:
        with open(output_fasta, "w") as fasta_out:
            for genome_id, protein_ids in genome_hits.items():
                faa_path = faa_files.get(genome_id)
                if not faa_path or not os.path.isfile(faa_path):
                    print(f"[ERROR] .faa file not found for genome '{genome_id}'")
                    continue

                try:
                    with open(faa_path, "r") as faa_in:
                        write_sequence = False
                        for raw_line in faa_in:
                            if raw_line.startswith(">"):
                                header = raw_line[1:].split()[0]
                                write_sequence = header in protein_ids
                                if write_sequence:
                                    fasta_out.write(f">{genome_id}___{header}\n")
                            elif write_sequence:
                                fasta_out.write(raw_line)
                except OSError as exc:
                    print(f"[ERROR] Cannot read FASTA '{faa_path}': {exc}")
                    continue
    except OSError as exc:
        print(f"[ERROR] Cannot write output FASTA '{output_fasta}': {exc}")
        return


def generate_faa_per_hitfile_parallel(
    options: Any,
    intermediate_hit_dir: str,
    processes: int = 4
) -> None:
    """
    For each '.intermediate_hits' file in the given directory, run `process_hitfile` in parallel
    to extract matching protein sequences into separate FASTA files.

    Args:
        options: Object with at least the attribute `faa_files` (Dict[str, str]), 
                 mapping genome IDs to their .faa file paths.
        intermediate_hit_dir: Directory containing '.intermediate_hits' files. 
                              Also serves as the output directory for FASTA files.
        processes: Number of parallel worker processes to use (default: 4).

    Raises:
        OSError: If `intermediate_hit_dir` does not exist or is not a directory.
    """
    # Ensure the intermediate_hit_dir exists
    if not os.path.isdir(intermediate_hit_dir):
        raise OSError(f"Directory not found: {intermediate_hit_dir}")

    # Make a list of all '.intermediate_hits' files in that directory
    hitfiles: List[str] = [
        os.path.join(intermediate_hit_dir, filename)
        for filename in os.listdir(intermediate_hit_dir)
        if filename.endswith(".intermediate_hits")
    ]

    # Prepare arguments for each worker: (hitfile_path, output_dir, faa_files_dict)
    args: List[Tuple[str, str, dict]] = [
        (hitfile_path, intermediate_hit_dir, options.faa_files)
        for hitfile_path in hitfiles
    ]

    # Execute process_hitfile in parallel
    with Pool(processes=processes) as pool:
        pool.starmap(process_hitfile, args)

##########################################################################################################################################################
#################### Cross check hits with reference sequences and add the hmmreport lines to the trusted cutoff hmmreport ###############################
##########################################################################################################################################################

def cross_check_candidates_with_reference_seqs(options: Any) -> List[str]:
    """
    Cross-checks intermediate hit FASTA files against reference sequence databases using DIAMOND.

    For each '<hmm_id>.intermediate_hits_faa' in the Cross_check_directory:
      1. Attempts to locate an existing Diamond database '<hmm_id>.dmnd' under 'src/RefSeqs'.
      2. If not found, tries to locate '<hmm_id>.faa' under 'src/RefSeqs' and build a Diamond DB.
      3. Runs 'diamond blastp' using the intermediate FASTA as query against the Diamond DB.
      4. Writes the BLASTP output to '<hmm_id>.crosschecked.tsv' in Cross_check_directory.
      5. If the output file is empty, deletes it.

    Args:
        options: An object with at least these attributes:
            - execute_location (str): Base directory for 'src/RefSeqs'.
            - Cross_check_directory (str): Directory containing intermediate FASTA files.
            - refseq_identity (float): Minimum percent identity for Diamond search.
            - cores (int): Number of threads for Diamond.
    
    Returns:
        A list of HMM IDs (strings) for which no valid reference sequence DB or FASTA was found.
    """
    print("\n\n[INFO] Cross-checking hit sequences with reference sequences\n")

    refseq_root = os.path.join(options.execute_location, "src", "RefSeqs")
    unavailable_hmms: List[str] = []
    cross_check_dir = options.cross_check_directory

    # Collect all intermediate_hits_faa files
    intermediate_patterns = os.path.join(cross_check_dir, "*.intermediate_hits_faa")
    intermediate_files = glob.glob(intermediate_patterns)

    for inter_faa_path in intermediate_files:
        basename = os.path.basename(inter_faa_path)
        hmm_id = basename.replace(".intermediate_hits_faa", "")
        diamond_db_basename = os.path.join(refseq_root, hmm_id)
        diamond_db_path = f"{diamond_db_basename}.dmnd"

        # If Diamond DB doesn't exist, try to build it from a .faa under RefSeqs
        if not os.path.isfile(diamond_db_path):
            query_faa_basename = os.path.join(refseq_root, f"{hmm_id}.faa")
            if os.path.isfile(query_faa_basename):
                print(f"[INFO] Creating Diamond DB from '{query_faa_basename}' for HMM '{hmm_id}'")
                try:
                    subprocess.run(
                        ["diamond", "makedb", "--in", query_faa_basename, "-d", diamond_db_basename],
                        check=True
                    )
                    diamond_db_path = f"{diamond_db_basename}.dmnd"
                except subprocess.CalledProcessError:
                    print(f"[ERROR] Failed to create Diamond DB for '{query_faa_basename}'")
                    unavailable_hmms.append(hmm_id)
                    continue
            else:
                print(f"[WARN] Skipping HMM ID '{hmm_id}': No '.dmnd' or '.faa' found under '{refseq_root}'")
                unavailable_hmms.append(hmm_id)
                continue

        # Prepare output TSV path
        output_tsv = os.path.join(cross_check_dir, f"{hmm_id}.crosschecked.tsv")
        diamond_exe = myUtil.find_executable("diamond")
        cmd = [
            diamond_exe,
            "blastp",
            "--query", inter_faa_path,
            "--db", diamond_db_path,
            "--out", output_tsv,
            "--outfmt", "6",
            "--max-target-seqs", "1",
            "--id", str(options.refseq_identity),
            "--threads", str(options.cores),
            "--quiet"
        ]

        print(f"[INFO] Verifying hits for HMM ID '{hmm_id}' against reference DB")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"[ERROR] DIAMOND search failed for '{hmm_id}'")
            continue

        # If output exists but is empty, remove it
        if os.path.isfile(output_tsv) and os.path.getsize(output_tsv) == 0:
            os.remove(output_tsv)

    # returns a list of hmms where no .dmnd for cross check was available
    return unavailable_hmms
    
    
def promote_crosschecked_hits(crosscheck_dir: str, processes: int = 4) -> None:
    """
    Promotes hits from '.intermediate_hits' to '.trusted_hits' based on the
    results in '.crosschecked.tsv' files found in the specified directory.

    This function:
      1. Scans the directory for all '*.crosschecked.tsv' files.
      2. Extracts the HMM IDs from those filenames.
      3. Calls `process_crosscheck(hmm_id, crosscheck_dir)` in parallel
         for each HMM ID.

    Args:
        crosscheck_dir: Path to the directory containing '.crosschecked.tsv' files.
        processes: Number of parallel worker processes to use (default: 4).

    Raises:
        OSError: If `crosscheck_dir` does not exist or cannot be listed.
    """
    print(f"\n\n[INFO] Analyzing cross reference hits between trusted and noise cutoff")
    
    if not os.path.isdir(crosscheck_dir):
        raise OSError(f"Directory not found: {crosscheck_dir}")

    # Collect HMM IDs by stripping the '.crosschecked.tsv' suffix
    crosschecked_files: List[str] = [
        filename
        for filename in os.listdir(crosscheck_dir)
        if filename.endswith(".crosschecked.tsv")
    ]

    suffix = ".crosschecked.tsv"
    hmm_ids: List[str] = [
        filename[:-len(suffix)] if filename.endswith(suffix) else filename
        for filename in crosschecked_files
    ]


    # Prepare argument tuples for each HMM ID
    args = [(hmm_id, crosscheck_dir) for hmm_id in hmm_ids]

    # Run process_crosscheck in parallel
    with Pool(processes=processes) as pool:
        pool.starmap(process_crosscheck, args)

def process_crosscheck(hmm_id: str, crosscheck_dir: str) -> None:
    """
    Promotes intermediate hits to trusted hits based on cross-checked results.

    This function:
      1. Reads the '<hmm_id>.crosschecked.tsv' file to collect valid hit IDs.
      2. Reads the '<hmm_id>.intermediate_hits' file and appends any lines whose
         ID appears in the valid set to '<hmm_id>.trusted_hits'.
      3. Prints the number of promoted hits.

    Args:
        hmm_id:         Identifier of the HMM.
        crosscheck_dir: Directory containing the crosscheck and intermediate files.

    Output:
        - '<hmm_id>.trusted_hits' wird um die gültigen Zeilen ergänzt.
    """
    crosscheck_path = os.path.join(crosscheck_dir, f"{hmm_id}.crosschecked.tsv")
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    if not os.path.isfile(crosscheck_path):
        print(f"[ERROR] Crosscheck file missing: '{crosscheck_path}'")
        return

    if not os.path.isfile(intermediate_path):
        print(f"[ERROR] Intermediate file missing: '{intermediate_path}'")
        return

    # 1. Load valid hit IDs from the crosscheck file
    valid_hits: Set[str] = set()
    try:
        with open(crosscheck_path, "r") as cross_in:
            for raw_line in cross_in:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue
                hit_id = raw_line.strip().split("\t")[0]
                valid_hits.add(hit_id)
    except OSError as exc:
        print(f"[ERROR] Cannot read '{crosscheck_path}': {exc}")
        return

    if not valid_hits:
        # Keine gültigen Hits → nichts zu promoten
        return

    promoted_count = 0

    # 2. Append valid hits from the intermediate file into the trusted file
    try:
        with open(intermediate_path, "r") as interm_in, open(trusted_path, "a") as trusted_out:
            for raw_line in interm_in:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue
                hit_id = raw_line.strip().split("\t")[0]
                if hit_id in valid_hits:
                    trusted_out.write(raw_line)
                    promoted_count += 1
    except OSError as exc:
        print(f"[ERROR] Error processing files for HMM '{hmm_id}': {exc}")
        return

    print(f"[INFO] {promoted_count} additional hits found for {hmm_id} based on reference comparison")

    
    
def summarize_trusted_hits(
    directory: str,
    crosscheck_dir: str,
    name: str,
    suffix: str = ".trusted_hits"
) -> str:
    """
    Concatenates all trusted-hit files in `crosscheck_dir` into a single summary
    file named `name` under `directory`. If the summary already exists and is
    non-empty, it is returned immediately. Otherwise, the function collects all
    files ending with `suffix`, concatenates their contents (in alphabetical order),
    and writes them to the summary file.

    Args:
        directory:     Directory where the summary file should be written.
        crosscheck_dir: Directory containing individual trusted-hit files.
        name:          Filename for the summary (e.g., "global_trusted_hits_summary.hmmreport").
        suffix:        File suffix used to identify trusted-hit files (default: ".trusted_hits").

    Returns:
        The absolute path to the summary file.

    """
    summary_path = os.path.join(directory, name)
    if os.path.isfile(summary_path) and os.path.getsize(summary_path) != 0:
    	return summary_path
    
    trusted_files = [
        os.path.join(crosscheck_dir, f)
        for f in os.listdir(crosscheck_dir)
        if f.endswith(suffix)
    ]

    exit_code = os.system(f"cat {' '.join(trusted_files)} > {summary_path}")
    
    return summary_path
    
    

    
####################################################################################################
#################### Promote by cutoff when cross check is not available ###########################
####################################################################################################

def process_optimized_cutoff(
    hmm_id: str,
    crosscheck_dir: str,
    optimized_dict: Dict[str, float]
) -> None:
    """
    Promotes hits from '.intermediate_hits' to '.trusted_hits' based on an
    optimized score cutoff when no cross-check is available.

    Steps:
      1. Construct paths for intermediate and trusted hit files.
      2. Determine the threshold score for this HMM ID from optimized_dict,
         defaulting to 50 if not present.
      3. If the intermediate file does not exist, warn and return.
      4. Read the intermediate file line by line:
         - Skip empty lines or comments.
         - Parse the bit score from column 8 (0-based index 7).
         - If score ≥ threshold, append the line to trusted file.
      5. Print how many hits were promoted.

    Args:
        hmm_id:         Identifier of the HMM.
        crosscheck_dir: Directory containing '.intermediate_hits' files.
        optimized_dict: Mapping from HMM ID to its optimized score cutoff.
    """
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    threshold_score = optimized_dict.get(hmm_id, 50.0)

    if not os.path.isfile(intermediate_path):
        print(f"[WARN] Intermediate hit file not found for HMM ID '{hmm_id}'")
        return

    promoted_count = 0

    try:
        with open(intermediate_path, "r") as interm_in, open(trusted_path, "a") as trusted_out:
            for raw_line in interm_in:
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.strip().split("\t")
                if len(parts) < 8:
                    # Not enough columns to extract score
                    continue

                try:
                    score = float(parts[7])
                except ValueError:
                    # If score cannot be parsed, skip the line
                    continue

                if score >= threshold_score:
                    trusted_out.write(raw_line)
                    promoted_count += 1
    except OSError as exc:
        print(f"[ERROR] Failed to process '{intermediate_path}' or write '{trusted_path}': {exc}")
        return

    print(f"[INFO] {hmm_id}: Promoted {promoted_count} hits (threshold = {threshold_score})")

    
def promote_by_cutoff(options, directory, processes: int = 4, hmm_ids: Optional[List[str]] = None):

    """
    Parallelisiert das Promoten von Hits aus .intermediate_hits zu .trusted_hits
    anhand der crosschecked.tsv-Dateien im gegebenen Verzeichnis.
    """
    
    print("\n\n[INFO] Selecting proteins that were not cross checked by optimized threshold")
    
    # If no hmm identifier were defined use all that are in
    if hmm_ids=="all":
        hmm_ids = [
            f.replace(".intermediate_hits", "")
            for f in os.listdir(directory)
            if f.endswith(".intermediate_hits")
        ]

    optimized_dict = make_threshold_dict(
        options.score_threshold_file, 1, options.thrs_score
    )
    
    with Pool(processes=processes) as pool:
        pool.starmap(
            process_optimized_cutoff,
            [(hmm_id, directory, optimized_dict) for hmm_id in hmm_ids]
        )

































#!/usr/bin/python
import os
import subprocess
from multiprocessing import Pool,Value, Lock
from collections import defaultdict



from . import myUtil
import glob

# Global shared variables
current_counter = None
counter_lock = None


########################################################################################################
#################### HMMsearch routines for individual file input ######################################
########################################################################################################

def unified_search(options, processes=4):
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

    return dict(zip(options.queued_genomes, results))
    
    
    
def init_globals(counter, lock):
    global current_counter
    global counter_lock
    current_counter = counter
    counter_lock = lock



def run_search(faa_file, query_db, score, clean_reports, redo_search, total):

    with counter_lock:
        counter = current_counter.value + 1
        current_counter.value = counter
        print(f"Processing file {counter} of {total}", end="\r")
        
    domtblout_path = HMMsearch(faa_file, query_db, score, clean_reports, redo_search, 2)
    hmmreport = prefix_domtblout_hits(domtblout_path, separator="___", suffix=".hmmreport")
    return hmmreport
    
    
def HMMsearch(path,query_db,score,clean_reports=False, redo_search=False, cores = 1):

    output = os.path.splitext(path)[0] + '.domtblout'
    hmmreport = os.path.splitext(path)[0] + '.hmmreport'
    if not os.path.isfile(hmmreport) or clean_reports or redo_search:
        os.system(f'hmmsearch -T {score} --domT {score} --cpu {str(cores)} --noali --domtblout {output} {query_db} {path} > /dev/null 2>&1')
    return output


def prefix_domtblout_hits(domtblout_path, separator="___", suffix=".hmmreport"):
    """
    Schreibt eine neue domtblout-Datei, in der Spalte 1 (Hit-ID) mit dem Dateinamen prefixiert wird.
    Nutzt shared counter nur für Fortschrittsanzeige.
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


def concatenate_hmmreports_cat(report_paths, output_path="global_report.cat_hmmreport"):
    """
    Verwendet das UNIX 'cat'-Kommando, um .hmmreport-Dateien zu einer globalen Datei zusammenzuführen.

    Parameters:
        report_paths (dict): Liste von Pfaden zu .hmmreport-Dateien
        output_path (str): Pfad zur Ausgabedatei
    
    Returns:
        str: Pfad zur globalen Datei
    """
    
    print("Concatenate hit reports")
    
    # Filtere nur existierende Dateien
    valid_paths = [path for path in report_paths.values() if os.path.isfile(path)]
    if not valid_paths:
        raise FileNotFoundError(f"No valid hmmreport files found in directory {report_paths}.")

    # Führe das cat-Kommando aus
    cmd = ["cat"] + valid_paths
    with open(output_path, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile)

    return output_path
    
########################################################################################################
#################### Filter the glob report to trusted hits and potential hits #########################
########################################################################################################
    
    
def make_threshold_dict(file_path, threshold_type=1, default_score=50.0):
    """
    Builds a dictionary of thresholds from a tab-separated file.
    
    Rules:
        - If only one score is present, use that.
        - If multiple scores are present, use the one matching threshold_type.
        - If no score is present, use default_score.
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
                print(f"[Line {line_number}] Problem parsing: {line.strip()} — {e}")
                continue

            if key:
                thresholds[key] = score

    return thresholds


def process_single_hmm(hmm_id, glob_report, trusted_cutoff, noise_cutoff, output_dir):

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
            except ValueError:
                continue


            if score >= trusted_cutoff:
                trusted_out.write(line)
            elif score > noise_cutoff:
                candidates[target] = score
    
    #Remove empty files and skip if no candidates found
    if os.path.getsize(trusted_path) == 0:
        os.remove(trusted_path)
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


def filter_trusted_and_noise_hits(options, processes=4):
    glob_report = options.glob_report
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



def extract_fasta_per_intermediate_hitfile(options, intermediate_hit_dir):
    """
    Für jede .intermediate_hits Datei wird ein separates FASTA-File mit den
    passenden Proteinsequenzen erstellt – ohne Umbrechen oder Biopython.
    
    Die Ausgabe wird zeilenweise verarbeitet und nur passende Sequenzen werden
    direkt geschrieben. Minimaler Speicherverbrauch.
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

        # Sequentielles Schreiben pro Genome
        with open(output_fasta, 'w') as out:
            for genome_id, protein_ids in genome_hits.items():
                faa_path = options.faa_files.get(genome_id)
                if not faa_path or not os.path.isfile(faa_path):
                    print(f"[!] .faa file not found for genome: {genome_id}")
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



def process_hitfile(hitfile_path, intermediate_hit_dir, faa_files):
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
                print(f"Warning: FASTA not found for {genome_id}")
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


def generate_faa_per_hitfile_parallel(options, intermediate_hit_dir, processes=4):
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

def find_refseq_file(base_dir, filename):
    """
    Sucht rekursiv in base_dir nach genau einer Datei mit dem gegebenen Dateinamen.
    Gibt den vollständigen Pfad zurück, wenn gefunden, sonst None.
    """
    for root, _, files in os.walk(base_dir):
        if filename in files:
            return os.path.join(root, filename)
    return None


def cross_check_candidates_with_reference_seqs(options):
    print("Cross check hit sequences with reference sequences")    

    refseq_dir = os.path.join(options.execute_location, "src", "RefSeqs")
    refseq_unavailable_list = []
    
    cross_check_dir = options.Cross_check_directory # directoy with the intermediate hit fasta faa files

    intermediate_files = glob.glob(os.path.join(cross_check_dir, "*.intermediate_hits_faa"))

    # Iterate the intermediate faa files
    for inter_file in intermediate_files:
        hmm_id = os.path.splitext(os.path.basename(inter_file))[0].replace(".intermediate_hits_faa", "")
        db_file = f"{hmm_id}.dmnd"
        db_path = find_file_in_prefixed_subdirs(refseq_dir, db_file, dir_prefix="") #dir_prefix is for version control, possibly uneccessary

        # Prüfen ob .dmnd existiert, sonst erstellen
        if not os.path.isfile(db_path):
            faa_file = f"{hmm_id}.faa"
            faa_path = find_file_in_prefixed_subdirs(refseq_dir, faa_file, dir_prefix="")

            if os.path.isfile(faa_path):
                print(f"Info: Creating Diamond DB from {faa_path} because {db_path} was not found")

                try:
                    subprocess.run(["diamond", "makedb", "--in", faa_path, "-d", db_base], check=True)
                    db_path = db_base + ".dmnd"
                except subprocess.CalledProcessError:
                    print(f"Error: Failed to create Diamond database for {faa_path}")
                    refseq_unavailable_list.append(hmm_id)
                    continue
            else:
                print(f"Warning: Skipping {hmm_id}: Neither .dmnd nor .faa file found recursively.")
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

        print(f"Verifying {hmm_id} hits with reference sequences")
        result = subprocess.run(cmd)

        if result.returncode != 0:
            print(f"Error: DIAMOND search failed for {hmm_id}")
            continue
        if os.path.getsize(output_file) == 0:
            os.remove(output_file)

    return refseq_unavailable_list
    
    
def find_file_in_prefixed_subdirs(base_dir, filename, dir_prefix):
    for root, dirs, files in os.walk(base_dir):
        # Nur Verzeichnisse mit dem gewünschten Prefix betreten
        if not os.path.basename(root).startswith(dir_prefix):
            continue

        if filename in files:
            return os.path.join(root, filename)

    return ""  # nicht gefunden


def process_crosscheck(hmm_id, crosscheck_dir):
    crosscheck_path = os.path.join(crosscheck_dir, f"{hmm_id}.crosschecked.tsv")
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")

    if not os.path.exists(crosscheck_path):
        print(f"ERROR: Crosscheck file missing: '{crosscheck_path}'")
    if not os.path.exists(intermediate_path):
        print(f"ERROR: Intermediate file missing: '{intermediate_path}'")

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

    print(f"{hmm_id}: promoted {promoted_count} hits via comparison with reference sequences")


def promote_crosschecked_hits(crosscheck_dir, processes=4):
    """
    Parallelisiert das Promoten von Hits aus .intermediate_hits zu .trusted_hits
    anhand der crosschecked.tsv-Dateien im gegebenen Verzeichnis.
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
    
    
def summarize_trusted_hits(directory, crosscheck_dir, name ,suffix=".trusted_hits"):
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

def process_optimized_cutoff(hmm_id, crosscheck_dir, optimized_dict):
    intermediate_path = os.path.join(crosscheck_dir, f"{hmm_id}.intermediate_hits")
    trusted_path = os.path.join(crosscheck_dir, f"{hmm_id}.trusted_hits")
    
    threshold_score = optimized_dict.get(hmm_id, 50)
    
    if not os.path.isfile(intermediate_path):
        print(f"Warning: Intermediate or trusted hit file missing for {hmm_id}")
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

    print(f"{hmm_id}: promoted {promoted_count} hits due to the given threshold {threshold_score}")

    
def promote_by_cutoff(options, directory, processes=4, hmm_ids=[]):
    """
    Parallelisiert das Promoten von Hits aus .intermediate_hits zu .trusted_hits
    anhand der crosschecked.tsv-Dateien im gegebenen Verzeichnis.
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

































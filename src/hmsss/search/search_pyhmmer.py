from __future__ import annotations
import time
from multiprocessing.queues import SimpleQueue

import pyhmmer

from multiprocessing import Pool
from typing import Dict, List, Union, Optional, Any

from hmsss.cli.config import Config
from hmsss.cross_check import generate_cross_check_fasta

from hmsss.fasta_preparation.materialize import materialize_pair_gz_next_to_input

from hmsss.parse_reports import parse_reports, pattern_completion_synteny, pattern_completion_pathway
from hmsss.parse_reports import csb_trie_algorithm
from hmsss.parse_reports import csb_finder
from hmsss.db import database

# --- pyhmmer imports (wenn du wirklich pyhmmer nutzt) ---

from hmsss.core.logging import get_logger
from hmsss.parse_reports.csb_trie_algorithm import TrieIndex
from hmsss.parse_reports.parse_reports import Protein

logger = get_logger(__name__)



def make_threshold_dict(
    file_path: str,
    threshold_type: int = 1,
    default_score: float = 50.0,
) -> Dict[str, Union[float, Dict[str, float]]]:
    """
    Parse a tab-separated cutoff table.

    Column meaning (0-based after HMM ID):
        1: optimized cutoff
        2: trusted cutoff
        3: noise cutoff

    threshold_type:
        0 -> return ALL thresholds as dict:
             {hmm_id: {"optimized": x, "trusted": y, "noise": z}}
        1 -> optimized only   (backward compatible)
        2 -> trusted only
        3 -> noise only
    """
    thresholds: Dict[str, Union[float, Dict[str, float]]] = {}

    def _parse(val: str) -> float:
        if val == "-inf":
            return 5000.0  # sentinel: unreachable
        return float(val)
    with open(file_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue

            hmm_id = parts[0]



            try:
                # normalize columns
                optimized = default_score
                trusted   = default_score
                noise     = default_score

                if len(parts) > 1:
                    optimized = _parse(parts[1])
                if len(parts) > 2:
                    trusted = _parse(parts[2])
                if len(parts) > 3:
                    noise = _parse(parts[3])

                if threshold_type == 0:
                    thresholds[hmm_id] = {
                        "optimized": optimized,
                        "trusted": trusted,
                        "noise": noise,
                    }
                elif threshold_type == 1:
                    thresholds[hmm_id] = optimized
                elif threshold_type == 2:
                    thresholds[hmm_id] = trusted
                elif threshold_type == 3:
                    thresholds[hmm_id] = noise
                else:
                    raise ValueError(f"Unknown threshold_type={threshold_type}")

            except (ValueError, IndexError) as e:
                logger.warning(
                    f"[Line {line_number}] Problem parsing: {line.strip()} — {e}"
                )
                continue

    return thresholds


# -----------------------------
# Helpers: batching
# -----------------------------
def split_into_batches(items: List[str], n_batches: int) -> List[List[str]]:
    if n_batches <= 0:
        raise ValueError("n_batches must be > 0")
    batches = [[] for _ in range(n_batches)]
    for i, item in enumerate(items):
        batches[i % n_batches].append(item)
    return batches


# -----------------------------
# Globals inside workers (loaded once per worker)
# -----------------------------
_G_HMMS: list = []
_G_THRESH: dict = {}
_G_CSB_PATTERNS = None
_G_COOCCURRENCE: Optional[Dict[str, tuple[set, int]]] = None
_G_INDEX_TRIE: TrieIndex = None
_G_CONFIG_LIGHT = None  # falls du einzelne Config-Flags brauchst
_G_OPT_THRESH: dict[str, float] = {}

_G_FAA_FILES: Optional[Dict[str, str]] = None
_G_GFF_FILES: Optional[Dict[str, str]] = None

_G_QUEUE: Optional[SimpleQueue] = None
_G_NUCLEOTIDE_RANGE: Optional[int] = None
_G_MIN_COMPLETENESS: Optional[float] = None
_G_USE_SYNTENY_COMPLETION: Optional[bool] = None

def _build_optimized_suffix_thresholds(
    threshold_dict: dict,
) -> dict[str, float]:
    """
    Build optimized-cutoff dict keyed by suffix (after first '_').
    If multiple HMMs map to the same suffix, keep the LOWER optimized cutoff.
    """
    opt_dict: dict[str, float] = {}

    for hmm_id, thr in threshold_dict.items():
        # optimized cutoff extrahieren (robust)
        opt = thr.get("optimized")

        parts = hmm_id.split("_", 1)
        suffix = parts[1] if len(parts) == 2 else hmm_id

        if suffix in opt_dict:
            opt_dict[suffix] = min(opt_dict[suffix], opt)
        else:
            opt_dict[suffix] = opt

    return opt_dict


def _init_worker(
    hmm_path: str,
    threshold_dict: Dict[str, float],
    config_light: dict,
    faa_files: Dict[str, str],
    gff_files: Dict[str, str],
    nucleotide_range: int,
    min_completeness: float,
    use_synteny_completion: bool,
):
    """
    Lädt schwere/konstante Daten einmal pro Worker.
    - HMMlib wird genau einmal in RAM gebracht
    - threshold_dict wird im Worker verfügbar gemacht
    - CSB naming index einmal bauen
    """
    global _G_HMMS, _G_THRESH,_G_OPT_THRESH, _G_CSB_PATTERNS, _G_COOCCURRENCE, _G_INDEX_TRIE, _G_CONFIG_LIGHT, _G_FAA_FILES, _G_GFF_FILES, _G_NUCLEOTIDE_RANGE, _G_MIN_COMPLETENESS, _G_USE_SYNTENY_COMPLETION

    _G_CONFIG_LIGHT = config_light
    _G_THRESH = threshold_dict
    _G_OPT_THRESH = _build_optimized_suffix_thresholds(threshold_dict)

    _G_FAA_FILES = faa_files
    _G_GFF_FILES = gff_files

    _G_NUCLEOTIDE_RANGE = nucleotide_range
    _G_MIN_COMPLETENESS = min_completeness
    _G_USE_SYNTENY_COMPLETION = use_synteny_completion

    # 1) HMMlib laden (einmal pro Worker)
    with pyhmmer.plan7.HMMFile(hmm_path) as hf:
        _G_HMMS = list(hf)

    # 2) optional: Cutoffs aus threshold_dict in HMMs setzen (wenn du trusted/noise pro HMM willst)
    #    Hier beispielhaft: trusted Cutoff als "trusted" in hmm.cutoffs eintragen.
    #    (Noise brauchst du ggf. separat; je nachdem wie du klassifizierst.)
    for hmm in _G_HMMS:
        hmm_id = hmm.name.decode() if isinstance(hmm.name, (bytes, bytearray)) else str(hmm.name)
        thr = _G_THRESH.get(hmm_id)
        if not thr:
            noise = 10
            trusted = 1000
        else:
            noise = thr.get("noise")
            trusted = thr.get("trusted")

        if noise is not None:
            # (seq, dom) beide auf denselben Wert, falls du keinen getrennten dom-noise hast
            hmm.cutoffs.noise = (float(noise), float(noise))
        else:
            logger.warning("There was not cutoff for {hmm_id} defined")

        if trusted is not None:
            hmm.cutoffs.trusted = (float(trusted), float(trusted))
        else:
            logger.warning("There was not cutoff for {hmm_id} defined")

    # 3) CSB Patterns + Trie einmal laden/bauen (wie in parse_reports.main_parse_summary_hmmreport) :contentReference[oaicite:4]{index=4}
    patterns_file = config_light["patterns_file"]
    cooc_file = config_light["cooccurrence_file"]

    _G_CSB_PATTERNS = csb_finder.make_pattern_dict(patterns_file)
    _G_COOCCURRENCE = csb_finder.make_pattern_dict(cooc_file)

    only_pattern_dict = {name: patset for name, (patset, _) in _G_CSB_PATTERNS.items()}

    # Trie Index für pattern hierarchical tree
    _G_INDEX_TRIE = csb_trie_algorithm.build_trie_index(only_pattern_dict)


#
# parse pyhmmer hits to protein dict
#

def add_pyhmmer_hits_to_protein_dict(
    *,
    genome_id: str,
    tophits_iter,
) -> dict[str, Protein]:
    """
    Füllt protein_dict mit Domains aus pyhmmer hmmsearch.
    - verwirft dom.score < noise (pro HMM)
    - setzt Tc für dom.score >= trusted (optional, wenn du das gleich markieren willst)
    - alle >= noise bleiben als Domains im Protein (wie im report-parser) :contentReference[oaicite:2]{index=2}
    """
    protein_dict = {}

    protein_object = Protein

    for tophits in tophits_iter:
        hmm_id = (
            tophits.query.name.decode()
            if isinstance(tophits.query.name, (bytes, bytearray))
            else str(tophits.query.name)
        )
        parts = hmm_id.split("_", 1)
        hmm_name = parts[1] if len(parts) == 2 else hmm_id

        thr = _G_THRESH.get(hmm_id)
        trusted = float(thr["trusted"]) if not thr is None else 1000

        for hit in tophits:
            prot_id = (
                hit.name.decode()
                if isinstance(hit.name, (bytes, bytearray))
                else str(hit.name)
            )
            score = hit.score

            start,end = 0,1
            for dom in hit.domains:
                # Koordinaten: HMMER/pyhmmer nutzt i.d.R. 1-based inkl. Endpunkt
                start = min(start, int(dom.env_from))
                end = max(end, int(dom.env_to))

            protein = protein_dict.get(prot_id)
            if protein is None:
                protein = protein_object(prot_id, hmm_name, start, end, int(score), genome_id)
                protein_dict[prot_id] = protein
            else:
                protein.add_domain(hmm_name, start, end, int(score))
            if score >= trusted:
                protein.add_selection_comment("Tc")
                protein.valid_hit = True

    return protein_dict


# -----------------------------
# Worker: verarbeitet ein Batch
# -----------------------------
def _process_genome1(
    genome_id: str,
) -> tuple[dict[str, Protein], dict[str, Any]] | None:
    """
    Pro Batch: pro Genom
      1) pyhmmer hmmsearch (HMMlib liegt global im Worker)
      2) parse gff, attach coords/strand/locustag
      3) gene cluster finden + benennen
      4) queue.put((protein_dict, cluster_dict))
    """

    try:
        faa_in = _G_FAA_FILES[genome_id]
        gff_in = _G_GFF_FILES[genome_id]

        with materialize_pair_gz_next_to_input(faa_in, gff_in) as (faa_file, gff_file):

            # load genome with all sequences into RAM
            abc = pyhmmer.easel.Alphabet.amino()
            with pyhmmer.easel.SequenceFile(faa_file, "fasta", digital=True, alphabet=abc) as sf:
                seqs = sf.read_block()

            tophits_iter = pyhmmer.hmmer.hmmsearch(_G_HMMS, seqs, cpus=2, bit_cutoffs="noise")

            # -- 2) combined protein hits, includes intermdiates, valid hits = >trusted cutoff ---
            protein_dict = add_pyhmmer_hits_to_protein_dict(
                genome_id=genome_id,
                tophits_iter=tophits_iter,
            )

            # --- 3) GFF parsing + attach info to all proteins ---
            parse_reports.parse_gff_file(gff_file, protein_dict)

            # --- 4) Gene clusters finden + benennen ---
            cluster_dict = csb_finder.find_syntenic_blocks(genome_id, protein_dict, _G_NUCLEOTIDE_RANGE)
            cluster_dict = csb_finder.name_syntenic_blocks_trie(cluster_dict, _G_INDEX_TRIE,
                                                                min_completeness=_G_MIN_COMPLETENESS)

            # --- 5) Pattern completion mechanism ---
            pattern_completion_synteny.enhance_syntenic_block_completeness(
                cluster_dict, protein_dict, _G_CSB_PATTERNS
            )

            # --- 6) Co-occurrence patterns added to valid hits
            if _G_USE_SYNTENY_COMPLETION:
                pattern_completion_pathway.enhance_pathway_completeness(
                        protein_dict, _G_COOCCURRENCE, _G_OPT_THRESH
                    )

            # --- 7) all named gene clusters are marked as valid hits --
            parse_reports.remove_unassigned_intermediate_proteins(protein_dict, set() , cluster_dict) # labelled alles was in gencluster liegt oder unter trusted fällt

            # --- 8) add sequences to proteins
            parse_reports.get_protein_sequence(faa_file, protein_dict)

            # --- 9) an writer prozess geben ---
            #_G_QUEUE.put((protein_dict, cluster_dict))
            return protein_dict, cluster_dict

    except Exception as e:
        logger.exception(f"Error processing {genome_id}: In search_pyhmmer _process_genome1:\n {e}")
        return {}, {}


def _process_genome(genome_id: str) -> tuple[dict[str, Protein], dict[str, Any]] | None:
    t_total0 = time.perf_counter()
    try:
        faa_in = _G_FAA_FILES[genome_id]
        gff_in = _G_GFF_FILES[genome_id]

        with materialize_pair_gz_next_to_input(faa_in, gff_in) as (faa_file, gff_file):

            # 1) load sequences
            t0 = time.perf_counter()
            with pyhmmer.easel.SequenceFile(faa_file, "fasta", digital=True) as sf:
                seqs = sf.read_block()
            t_load = time.perf_counter() - t0

            # 2) hmmsearch (NUR das)
            t0 = time.perf_counter()
            tophits_iter = pyhmmer.hmmer.hmmsearch(_G_HMMS, seqs, cpus=2, bit_cutoffs="noise")
            # wichtig: Iterator materialisieren, sonst misst du später unabsichtlich hmmsearch mit
            tophits_list = list(tophits_iter)
            t_hmmsearch = time.perf_counter() - t0

            # 3) parsing (NUR das)
            t0 = time.perf_counter()
            protein_dict = add_pyhmmer_hits_to_protein_dict(
                genome_id=genome_id,
                tophits_iter=tophits_list,
                protein_dict=None,
            )
            t_parse_hits = time.perf_counter() - t0

            # 4) Rest wie gehabt
            t0 = time.perf_counter()
            parse_reports.parse_gff_file(gff_file, protein_dict)
            t_gff = time.perf_counter() - t0

            t0 = time.perf_counter()
            cluster_dict = csb_finder.find_syntenic_blocks(genome_id, protein_dict, _G_NUCLEOTIDE_RANGE)
            cluster_dict = csb_finder.name_syntenic_blocks_trie(
                cluster_dict, _G_INDEX_TRIE, min_completeness=_G_MIN_COMPLETENESS
            )
            t_csb = time.perf_counter() - t0

            t0 = time.perf_counter()
            pattern_completion_synteny.enhance_syntenic_block_completeness(
                cluster_dict, protein_dict, _G_CSB_PATTERNS
            )
            t_syn = time.perf_counter() - t0

            t_path = 0.0
            if _G_USE_SYNTENY_COMPLETION:
                t0 = time.perf_counter()
                pattern_completion_pathway.enhance_pathway_completeness(
                    protein_dict, _G_COOCCURRENCE, _G_OPT_THRESH
                )
                t_path = time.perf_counter() - t0

            t0 = time.perf_counter()
            parse_reports.remove_unassigned_intermediate_proteins(protein_dict, set(), cluster_dict)
            t_rm = time.perf_counter() - t0

            t0 = time.perf_counter()
            parse_reports.get_protein_sequence(faa_file, protein_dict)
            t_seq = time.perf_counter() - t0

        t_total = time.perf_counter() - t_total0

        # Beispiel: 1% sampling log (sonst zu viel)
        if 0 == 0:
            logger.info(
                f"[timing {genome_id}] total={t_total:.3f}s "
                f"load={t_load:.3f} hmmsearch={t_hmmsearch:.3f} parse_hits={t_parse_hits:.3f} "
                f"gff={t_gff:.3f} csb={t_csb:.3f} syn={t_syn:.3f} path={t_path:.3f} "
                f"rm={t_rm:.3f} seq={t_seq:.3f}"
            )

        return protein_dict, cluster_dict

    except Exception as e:
        logger.exception(f"Error processing {genome_id}: {e}")
        return {},{}

# -----------------------------
# Writer: SQLite insert (du hast das bereits)
# -----------------------------
def process_writer(queue, config):
    # This routine handles the output of the search and writes it into the database
    # It gets input from multiple workers as the database connection to sqlite is unique
    total_genomes = len(config.queued_genomes)
    genomes_done = 0
    log_step = max(1, total_genomes // 100)

    protein_batch = {}
    cluster_batch = {}
    batch_size = config.glob_chunks
    batch_counter = 0

    while True:
        tup = queue.get()
        if tup is None:
            break

        else:
            batch_counter += 1
            logger.debug(f"Processed {batch_counter} genomes ")  #

        protein_dict, cluster_dict = tup

        # Concatenate the data
        protein_batch.update(protein_dict)
        cluster_batch.update(cluster_dict)

        genomes_done += 1
        if (genomes_done % log_step == 0) or (genomes_done == total_genomes):
            pct = (genomes_done * 100) // max(1, total_genomes)
            logger.info(f"[Genome progress] {genomes_done}/{total_genomes} ({pct}%) genomes processed")

        # Print text reports if desired
        #if config.individual_reports:
        #    if protein_dict:  # Check if protein_dict is not empty
        #        first_protein_key = next(iter(protein_dict))  # Get the first key
        #        genome_id = protein_dict[first_protein_key].genomeID
        #        filepath = os.path.join(
        #            config.fasta_initial_hit_directory,
        #            str(genome_id) + ".hit_table_txt",
        #        )
        #        parse_reports.output_genome_report(filepath, protein_dict, cluster_dict, {})

        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            database.insert_database_proteins(config.database_directory, protein_batch)
            database.insert_database_clusters(config.database_directory, cluster_batch)

            # write the intermediate hits for cross check
            generate_cross_check_fasta.write_intermediate_hits_faa(protein_batch, config.cross_check_directory, suffix=".intermediate_hits_faa", max_open_files=32)

            protein_batch.clear()
            cluster_batch.clear()
            batch_counter = 0


    # Submit the remaining and file reports
    if protein_batch or cluster_batch:
        database.insert_database_proteins(config.database_directory, protein_batch)
        database.insert_database_clusters(config.database_directory, cluster_batch)
        generate_cross_check_fasta.write_intermediate_hits_faa(protein_batch, config.cross_check_directory,
                                                               suffix=".intermediate_hits_faa", max_open_files=32)
        logger.info(f"Processed {batch_counter} genomes")
    return




# -----------------------------
# pyhmmer consecutive hmm search
# -----------------------------
def consecutive_hmm_search(config: Config, processes: int = 4) -> None:
    """
    Neue Pipeline-Variante:
      - bildet genomeID/Faa/Gff batches
      - initialisiert Worker mit (HMMlib, threshold_dict, CSB patterns)
      - Worker machen search + gff + cluster
      - Writer-Prozess schreibt in sqlite

      - alle proteine die entweder >TC oder IN names gene cluster oder Pattern completion
      werden als valid hits geführt
      - alle non valid hits werden in cross check dir geschrieben als intermediate faas

      Es bleibt danach auszuführen der cross check mit den Referenz sequenzen
    """

    # Insert genomeIDs into database Genomes table
    genome_ids = list(config.queued_genomes)
    database.insert_database_genome_ids(config.database_directory, set(genome_ids))

    threshold_dict = make_threshold_dict(config.score_threshold_file, 0, config.thrs_score)

    # Für CSB pattern building im worker initializer: nur “leichte” config-Infos
    config_light = {
        "patterns_file": config.patterns_file,
        "cooccurrence_file": config.cooccurrence_file,
    }

    n_genomes: int = len(genome_ids)
    hmm_cpus = 2
    worker_processes = max(1, (processes - 1) // hmm_cpus)
    worker_processes = min(worker_processes, n_genomes)
    #chunksize: int = 1 if n_genomes < 5 else 5
    chunksize = 1

    # Writer process that gets the dictionaries
    #ctx = multiprocessing.get_context()  # nutzt default start method (fork auf Linux)
    #q = ctx.SimpleQueue()
    #writer = Process(target=process_writer, args=(q, config), daemon=True)
    #writer.start()


    # Batch buffers im Main
    protein_batch: dict[str, Protein] = {}
    cluster_batch: dict = {}
    batch_size: int = config.glob_chunks
    batch_counter = 0

    # Fortschritt
    genomes_done = 0
    log_step = max(1, n_genomes // 100)

    with Pool(
        processes=worker_processes,
        initializer=_init_worker,
        initargs=(config.library, threshold_dict, config_light, config.faa_files, config.gff_files, config.nucleotide_range, config.min_completeness, config.use_synteny_completion), # arguments for the init worker
    ) as pool:
        for protein_dict, cluster_dict in pool.imap_unordered(_process_genome1, genome_ids, chunksize=chunksize):
            genomes_done += 1
            if (genomes_done % log_step == 0) or (genomes_done == n_genomes):
                pct = (genomes_done * 100) // max(1, n_genomes)
                logger.info(f"[Genome progress] {genomes_done}/{n_genomes} ({pct}%) genomes processed")

            # Batch sammeln
            protein_batch.update(protein_dict)
            cluster_batch.update(cluster_dict)
            batch_counter += 1

            # Flush
            if batch_counter >= batch_size:
                database.insert_database_proteins(config.database_directory, protein_batch)
                database.insert_database_clusters(config.database_directory, cluster_batch)
                generate_cross_check_fasta.write_intermediate_hits_faa(
                    protein_batch,
                    config.cross_check_directory,
                    suffix=".intermediate_hits_faa",
                    max_open_files=32,
                )
                protein_batch.clear()
                cluster_batch.clear()
                batch_counter = 0

            # Rest flush
        if protein_batch or cluster_batch:
            database.insert_database_proteins(config.database_directory, protein_batch)
            database.insert_database_clusters(config.database_directory, cluster_batch)
            generate_cross_check_fasta.write_intermediate_hits_faa(
                protein_batch,
                config.cross_check_directory,
                suffix=".intermediate_hits_faa",
                max_open_files=32,
            )


    logger.info("Finished pyhmmer search + parsing + DB write.")



from __future__ import annotations

import os
import math
import sys
from dataclasses import dataclass
from types import SimpleNamespace
from collections import defaultdict
from typing import Iterable, Tuple, Dict

from hmsss.core import queue
from hmsss.core.logging import get_logger
from hmsss.graft import prepare_packages, gpkg_length, gpkg_ram
from hmsss.utils import myUtil

logger = get_logger(__name__)


@dataclass
class GraftMTask:
    # takes all variable arguments for the read mapping task
    gpkg: str
    gpkg_name: str
    outdir: str
    threads: int = 5
    evalue: str = "1e-5"
    length: int = 1
    min_coverage: float = 0.3
    diamond_db = None
    decoy_db = None

    # Memory limits für das mapping
    mem_est_gb: float | None = None  # for scheduling / tokens
    mem_cap_gb: float | None = None  # hard per-process limit

    metagenome_id: str = ""
    genome_id: str = ""

    forward: str | None = None
    reverse: str | None = None  # nur zusammen mit forward
    interleaved: bool | None = None  # exklusiv (statt forward/reverse)

    def __post_init__(self):
        # For a given gpkg in the task define the paths to diamond and decoy database
        gpkg_dir = self.gpkg

        diamond_candidate = os.path.join(gpkg_dir, "refseq_database.dmnd")
        diamond_db = diamond_candidate if os.path.isfile(diamond_candidate) else None

        decoy_candidate = os.path.join(gpkg_dir, "decoy_database.dmnd")
        decoy_db = decoy_candidate if os.path.isfile(decoy_candidate) else None

        object.__setattr__(self, "diamond_db", diamond_db)
        object.__setattr__(self, "decoy_db", decoy_db)

    def apply_ram_profile(
            self,
            peak_gb: float | None,
            *,
            est_factor: float = 1.02,
            cap_factor: float = 1.1,
            fallback_est_gb: float = 16.0,
            fallback_cap_gb: float = 24.0,
            rounding: str = "ceil",
    ) -> None:
        """
        Set mem_est_gb and mem_cap_gb based on a measured peak RAM (GB).
        If peak_gb is None, use conservative fallback values.
        """
        if peak_gb is None:
            self.mem_est_gb = float(fallback_est_gb)
            self.mem_cap_gb = float(fallback_cap_gb)
            return

        est = float(peak_gb) * float(est_factor)
        cap = float(peak_gb) * float(cap_factor)

        if rounding == "ceil":
            self.mem_est_gb = float(math.ceil(est))
            self.mem_cap_gb = float(math.ceil(cap))
        elif rounding == "round":
            self.mem_est_gb = float(round(est))
            self.mem_cap_gb = float(round(cap))
        else:
            # no rounding
            self.mem_est_gb = est
            self.mem_cap_gb = cap


def merge_read_mapping_inputs(
        *,
        fna_files: Dict[str, str],
        faa_files: Dict[str, str],
        fastq_files: Dict[str, str],
) -> Dict[str, str]:
    """
    Merge FNA, FAA and FASTQ input dictionaries into a single mapping:

        { genome_id : path }

    Conflict priority:
        FAA  >  FNA  >  FASTQ
    """
    merged: Dict[str, str] = {}

    # lowest priority → add first
    merged.update(fastq_files)

    # medium priority → overwrite fastq if key exists
    merged.update(fna_files)

    # highest priority → overwrite fna + fastq
    merged.update(faa_files)

    return merged


def automatic_forward_reverse_file_detection(
        files: Dict[str, str],
        forward_extension: str | None = None,
        reverse_extension: str | None = None,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    # immer absolute Pfade

    files = [os.path.abspath(p) for p in files.values()]

    # ------------------------------------------------------------
    # 1) Endungen automatisch bestimmen
    # ------------------------------------------------------------
    if forward_extension is None or reverse_extension is None:
        suffix_counter = defaultdict(int)

        for path in files:
            name = os.path.basename(path)
            for token in ("_R1", "_R2", "_1", "_2"):
                if token in name:
                    suffix = name[name.index(token):]
                    suffix_counter[suffix] += 1

        if len(suffix_counter) == 0:
            # kein Muster erkannt → alles forward
            forward_extension = None
            reverse_extension = None

        elif len(suffix_counter) == 1:
            forward_extension = next(iter(suffix_counter))
            reverse_extension = None

        else:
            sorted_suffixes = sorted(
                suffix_counter.items(), key=lambda x: x[1], reverse=True
            )
            s1, s2 = sorted_suffixes[:2]

            if "R1" in s1[0] or s1[0].endswith("_1.fastq.gz"):
                forward_extension, reverse_extension = s1[0], s2[0]
            else:
                forward_extension, reverse_extension = s2[0], s1[0]

    # Generate forward and reverse dictionary
    forward: Dict[str, str] = {}
    reverse: Dict[str, str] = {}

    for path in files:
        name = os.path.basename(path)

        # Reverse found
        if reverse_extension and name.endswith(reverse_extension):
            base = myUtil.get_genome_id(name[: -len(reverse_extension)])
            reverse[base] = path

        # Forward found
        elif forward_extension and name.endswith(forward_extension):
            base = myUtil.get_genome_id(name[: -len(forward_extension)])
            forward[base] = path

        # Fallback: unmatched fastq files into forward reads, invalids are in later in try catch
        else:
            base = os.path.splitext(name)[0]
            forward[base] = path

    return forward, reverse


def build_graft_args(task: GraftMTask) -> SimpleNamespace:
    gpkg = str(task.gpkg)
    outdir = str(task.outdir)

    # --- GraftM erwartet forward/reverse oft als LISTEN (oder None) ---
    if task.interleaved:
        if task.forward is None:
            raise ValueError(
                "Interleaved mode: forward must be paths to interleaved FASTQ/FASTA files."
            )
        if task.reverse is not None:
            raise ValueError("Interleaved mode: reverse must not be set.")
        forward = None
        reverse = None
        interleaved = [str(task.forward)]  # GraftM CLI: --interleaved nargs='+'
        file = interleaved[0]

    else:
        # Non-interleaved: forward required, reverse optional
        if task.forward is None:
            raise ValueError("Non-interleaved mode: forward must be set.")
        forward = [str(task.forward)]
        reverse = [str(task.reverse)] if task.reverse is not None else None
        interleaved = None
        file = forward[0]

    # Define output directory. GraftM needs a non-existing one, to not overwrite results
    file = os.path.basename(file).split(".")[0]
    output = os.path.join(outdir, task.gpkg_name + "_" + file)

    return SimpleNamespace(
        # Dispatch
        subparser_name="graft",
        # Inputs
        graftm_package=gpkg,
        forward=forward,
        reverse=reverse,
        interleaved=interleaved,
        # running options
        input_sequence_type=None,
        # Output
        output_directory=output,
        force=False,
        verbosity=2,
        log=False,
        # Pipeline controls
        threads=int(task.threads),
        evalue=str(task.evalue),
        search_only=False,
        search_and_align_only=False,
        # Merge-Reads Verhalten deterministisch halten (optional)
        merge_reads=(task.reverse is not None and not task.interleaved),
        no_merge_reads=not (task.reverse is not None and not task.interleaved),
        # Search / assignment
        search_method="hmmsearch+diamond",
        assignment_method="pplacer",
        placements_cutoff=0.75,
        search_diamond_file=task.diamond_db,
        diamond_performance_parameters="",
        decoy_database=task.decoy_db,
        expand_search_contigs=None,
        maximum_range=None,
        filter_minimum=None,
        # ORF / translation
        min_orf_length=90,
        restrict_read_length=None,
        translation_table=11,
        euk_check=False,
        euk_hmm_file=None,
        # Wird später (durch HouseKeeping) aus dem gpkg gesetzt
        search_hmm_files=[],
        search_hmm_list_file=None,
        aln_hmm_file=None,
        reference_package=None,
        resolve_placements=False,
        max_samples_for_krona=100,
    )


def guess_extension(path: str) -> str:
    """
    Return the matched extension (including multi-part like '.fastq.gz').
    Raises ValueError if unknown.
    """
    # sort longest first so '.fastq.gz' matches before '.gz'
    exts = [
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".faa.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".fna",
        ".faa",
        ".gz",
    ]
    for ext in exts:
        if path.endswith(ext):
            return ext
    raise ValueError(f"Unable to guess file format of sequence file: {path}")


def read_basename(read_file: str) -> str:
    """
    Return filename without recognized sequencing extension.
    Example: '/x/y/sample_1.fq.gz' -> 'sample_1'
    """
    base = os.path.basename(read_file)
    ext = guess_extension(read_file)
    return base[: -len(ext)]


def create_task_list(
        gpkg_packages: dict[str, str],
        forward_dict: dict[str, str],
        reverse_dict: dict[str, str],
        output_directory: str,
        length_dict: dict[str, int],
        gpkg_ram_profile: dict[str, float],
        *,
        threads: int = 5,
        evalue: str = "1e-5",
        interleaved: bool = False,
) -> list[GraftMTask]:
    tasks: list[GraftMTask] = []

    for key, forward_path in forward_dict.items():
        reverse_path = reverse_dict.get(key)  # None falls nicht vorhanden
        metagenome_id = read_basename(forward_path)
        genome_id = myUtil.get_genome_id(forward_path)
        for gpkg_name, gpkg in gpkg_packages.items():
            length = length_dict.get(gpkg_name)
            peak = gpkg_ram_profile.get(gpkg_name, 4.0)  # default 4 GB RAM
            task = GraftMTask(
                gpkg_name=gpkg_name,
                gpkg=gpkg,
                length=length,
                metagenome_id=metagenome_id,
                genome_id=genome_id,
                forward=forward_path,
                reverse=reverse_path,
                interleaved=interleaved,
                outdir=output_directory,
                threads=threads,
                evalue=evalue,
            )
            task.apply_ram_profile(peak)

            tasks.append(task)

    return tasks


def filter_gpkg_by_ram(
        *,
        gpkg_packages: dict[str, str],
        gpkg_ram_profile: dict[str, float],
        ram_limit_gb: float,
        cap_factor: float = 1.1,
        fallback_cap_gb: float = 16.0,
) -> dict[str, str]:
    """
    Remove GPKGs that would exceed the RAM limit even in isolation.
    """
    kept: dict[str, str] = {}

    for gpkg_name, gpkg_path in gpkg_packages.items():
        peak = gpkg_ram_profile.get(gpkg_name)

        if peak is None:
            cap = fallback_cap_gb
        else:
            cap = math.ceil(float(peak) * cap_factor)

        if 0 < ram_limit_gb < cap:  # limit 0 means no limit
            logger.warning(
                f"[GPKG filter] Removing gpkg '{gpkg_name}': "
                f"Estimated needed RAM={cap} GB > RAM limit={ram_limit_gb} GB"
            )
            continue

        kept[gpkg_name] = gpkg_path

    return kept


def initialize_task_list(config):
    # _check_dependencies()  # Check if graftM dependencies are present

    combined_inputs = merge_read_mapping_inputs(
        fna_files=config.fna_files,
        faa_files=config.faa_files,
        fastq_files=config.fastq_files,
    )
    forward_dict, reverse_dict = automatic_forward_reverse_file_detection(
        files=combined_inputs, forward_extension=None, reverse_extension=None
    )
    logger.info(
        f"Forward read files: {len(forward_dict)} and reverse read files: {len(reverse_dict)}"
    )

    # Select sets and/or specific gpkgs from the library
    gpkg_packages = prepare_packages.prepare_gpkg_packages(config)
    gpkg_specifics = prepare_packages.prepare_gpkg_packs(config)
    gpkg_packages = gpkg_packages | gpkg_specifics
    logger.info(f"Initialized {len(gpkg_packages)} gpkg packages")
    logger.debug(f"{gpkg_packages}")
    if not len(gpkg_packages):
        logger.warning("Define a gpkg set for the read mapping")
        sys.exit()

    prepare_packages.initialize_gpkg_packages(gpkg_packages, threads=4)
    gpkg_length_dict = gpkg_length.collect_gpkg_reference_median_lengths(gpkg_packages)
    gpkg_ram_dict = gpkg_ram.collect_gpkg_ram(gpkg_packages)

    gpkg_packages = filter_gpkg_by_ram(
        gpkg_packages=gpkg_packages,
        gpkg_ram_profile=gpkg_ram_dict,
        ram_limit_gb=float(config.rm_ram_limit_max),
        # Ein min limit wäre auch sinnvoll wenn man gewisse gpkgs eineln laufen lassen will
    )

    # create task list can also define the cpu threads and the minimal e value
    task_list = create_task_list(
        gpkg_packages=gpkg_packages,
        forward_dict=forward_dict,
        reverse_dict=reverse_dict,
        output_directory=config.fasta_output_directory,
        length_dict=gpkg_length_dict,
        gpkg_ram_profile=gpkg_ram_dict,
        interleaved=config.rm_interleaved,
        threads=20,
    )

    return task_list

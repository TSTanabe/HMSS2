# src/hmsss/graft/adapter_core.py
# This is the adaptor between the main program to
# the implementation of graftM. GraftM code was adjusted. Now only
# read identification and placement are done and minimal output
# No other functionality including the package create is included
# Packages have to be provided by the user or the /data


from __future__ import annotations

import os
from types import SimpleNamespace
from typing import Dict

from hmsss.cli.config import Config
from hmsss.graft.adaptor_execute import Run
from hmsss.core.logging import get_logger
from hmsss.cli.paths import DATA_DIR
from hmsss.core import queue as queue

log = get_logger(__name__)


def prepare_gpkg_packages(config: Config) -> None | dict[str, str]:
    if config.hmm_sets:
        log.info(f"For the library collecting HMMs with token {config.hmm_sets}")
        allowed = (
            config.hmm_sets
            if isinstance(config.hmm_sets, list)
            else config.hmm_sets.split()
        )
        return queue.collect_gpkg_from_selected_metabolism_packages(
            str(DATA_DIR), allowed
        )
    return


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


def run_read_mapping(
    *,
    # --- Pflichtargumente ---
    forward: str,
    graftm_package: str,
    decoy_database: str = "",  # located in the package
    search_diamond_file: str = "",
    # --- HMSSS-Metainfo ---
    package_name: str | None = None,
    genome_id: str = "",
    # --- optional, aber nützlich ---
    reverse: list[str] | None = None,
    interleaved: list[str] | None = None,
    threads: int = 14,
    evalue: float = 1e-5,
    placements_cutoff: float = 0.75,
    resolve_placements: bool = False,
    min_orf_length: int = 96,
    restrict_read_length: int | None = None,
    translation_table: int = 11,
    euk_check: bool = False,
    input_sequence_type: str
    | None = None,  # z.B. UnpackRawReads.NUCLEOTIDE_SEQUENCE_TYPE
    output_directory: str = "results",
    force: bool = False,
    verbosity: int = 4,
    log: str | None = None,
) -> None:
    """
    HMSSS-Wrapper für GraftM 'graft' mit folgenden Festlegungen:

    - assignment_method: pplacer
    - search_method: hmmsearch+diamond
    - immer gpkg + Diamond + Decoy-DB
    - Krona deaktiviert
    """

    # Fixe Modi
    search_method = "hmmsearch+diamond"
    assignment_method = "pplacer"
    max_samples_for_krona = 0  # Krona effektiv aus

    # Parse the arguments into a namespace, needed by graftM Run class
    args = SimpleNamespace(
        subparser_name="graft",
        # I/O
        forward=forward,
        reverse=reverse,
        interleaved=interleaved,
        graftm_package=graftm_package,
        package_name=package_name,
        genome_id=genome_id,
        # Pipeline-Steuerung
        threads=threads,
        input_sequence_type=input_sequence_type,
        filter_minimum=None,
        evalue=evalue,
        search_and_align_only=False,
        search_only=False,
        euk_check=euk_check,
        search_method=search_method,
        decoy_database=decoy_database,
        maximum_range=None,
        expand_search_contigs=None,
        search_hmm_files=None,
        search_hmm_list_file=None,
        search_diamond_file=search_diamond_file,
        aln_hmm_file=None,
        assignment_method=assignment_method,
        placements_cutoff=placements_cutoff,
        resolve_placements=resolve_placements,
        no_merge_reads=False,
        diamond_performance_parameters="",
        euk_hmm_file=None,
        min_orf_length=min_orf_length,
        restrict_read_length=restrict_read_length,
        translation_table=translation_table,
        verbosity=verbosity,
        log=log,
        output_directory=output_directory,
        force=force,
        max_samples_for_krona=max_samples_for_krona,
    )

    Run(args).main()

    return


def read_mapping(config: Config):
    """

    Args:
        config:

    Returns:

    """
    # Define packages
    gpkg_dict = prepare_gpkg_packages(config)
    if not gpkg_dict:
        log.warning(
            f"No read_mapping_packages with .gpkg extension configured in any set: {config.hmm_sets}."
        )
        return

    # Define input files
    fna_files: Dict[str, str] = getattr(config, "fna_files", {}) or {}
    faa_files: Dict[str, str] = getattr(config, "faa_files", {}) or {}
    fastq_files: Dict[str, str] = getattr(config, "fastq_files", {}) or {}
    input_files = merge_read_mapping_inputs(
        fna_files=fna_files,
        faa_files=faa_files,
        fastq_files=fastq_files,
    )

    # Parameters from config, default if none
    rm_threads = getattr(config, "rm_threads", None) or 14
    rm_evalue = getattr(config, "rm_evalue", None) or 1e-5
    rm_placements_cutoff = getattr(config, "rm_placements_cutoff", None) or 0.75
    rm_resolve_placements = bool(getattr(config, "rm_resolve_placements", False))
    rm_min_orf_length = getattr(config, "rm_min_orf_length", None) or 96
    rm_restrict_read_length = getattr(config, "rm_restrict_read_length", None)
    rm_translation_table = getattr(config, "rm_translation_table", None) or 11

    # Iterate input files, graftM of each input against all
    for genome_id, input_path in input_files.items():
        for pkg_name, gpkg_dirpath in gpkg_dict.items():
            # config.result_files_directory
            # config.fasta_initial_hit_directory
            diamond_database = os.path.join(gpkg_dirpath, "Q")
            decoy_database = os.path.join(gpkg_dirpath, "7_bait_db.dmnd")
            try:
                run_read_mapping(
                    forward=input_path,
                    graftm_package=gpkg_dirpath,
                    decoy_database=decoy_database,
                    search_diamond_file=diamond_database,
                    package_name=pkg_name,
                    genome_id=genome_id,  # Paketname explizit übergeben
                    # Pipeline-Parameter
                    threads=rm_threads,
                    evalue=rm_evalue,
                    placements_cutoff=rm_placements_cutoff,
                    resolve_placements=rm_resolve_placements,
                    min_orf_length=rm_min_orf_length,
                    restrict_read_length=rm_restrict_read_length,
                    translation_table=rm_translation_table,
                    # IO / Logging
                    output_directory=config.fasta_initial_hit_directory,
                    force=False,
                    verbosity=getattr(config, "rm_verbosity", 4),
                    log=None,
                )

            except Exception:
                # ALLE Fehler abfangen, aber vollständig loggen
                log.exception(
                    "Read mapping failed for genome_id=%s with package=%s "
                    "(input=%s, gpkg=%s). Skipping this combination.",
                    genome_id,
                    pkg_name,
                    input_path,
                    gpkg_dirpath,
                )
                # continue → einfach zur nächsten Kombination
                continue

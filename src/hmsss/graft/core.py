# src/hmsss/graft/core.py

from __future__ import annotations
from types import SimpleNamespace
from hmsss.graft.execute import Run

def run_graft(
    *,
    # --- Pflichtargumente ---
    forward: list[str],
    graftm_package: str,
    decoy_database: str,

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
    input_sequence_type: str | None = None,  # z.B. UnpackRawReads.NUCLEOTIDE_SEQUENCE_TYPE

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
    search_method = 'hmmsearch+diamond'
    assignment_method = 'pplacer'
    max_samples_for_krona = 0  # Krona effektiv aus

    args = SimpleNamespace(
        subparser_name="graft",

        # I/O
        forward=forward,
        reverse=reverse,
        interleaved=interleaved,
        graftm_package=graftm_package,

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
        search_diamond_file=None,  # Diamond-DB kommt aus gpkg
        aln_hmm_file=None,         # wird von HouseKeeping.set_attributes gesetzt

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

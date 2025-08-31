# tests/test_config_build_from_cli.py
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

# Module importieren
from hmsss.cli import paths as p
from hmsss.cli import parse as cli_parse
from hmsss.cli.config import (
    Config,
    PathsCfg,
    CliInput,
    CliSearchParams,
    CliResources,
    CliSynteny,
    CliInfo,
    CliCsb,
    CliFlow,
    CliLimiter,
    CliOperators,
    CliProcess,
)


def _mk_fake_tree(root: Path) -> None:
    """Erzeuge die erwartete HMSS2-Struktur unter 'root'."""
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)
    # ein minimaler Eingabeordner für -f
    (root / "inputs").mkdir(parents=True, exist_ok=True)


def _paths_cfg_from_module() -> PathsCfg:
    """Baue einen PathsCfg aus dem aktuellen paths-Modul."""
    d = p.as_dict(str_paths=True)
    return PathsCfg(
        root=d["ROOT_DIR"],
        bin=d["BIN_DIR"],
        data=d["DATA_DIR"],
        hmms=d["HMMS_DIR"],
        refseq=d["REFSEQ_DIR"],
        results=d["RESULTS_DIR"],
        package=d["PACKAGE_DIR"],
    )


def _fill_config_from_namespace(ns) -> Config:
    """
    Minimaler Builder direkt im Test:
    Überträgt CLI-Werte in die Config-Blöcke und setzt sinnvolle Defaults aus paths.
    (Wenn du später einen offiziellen Builder hast, kannst du den hier verwenden.)
    """
    paths_cfg = _paths_cfg_from_module()

    # -------- CLI-Blöcke befüllen (nur die wichtigsten Felder; Rest bleibt per Default) --------
    cli_input = CliInput(
        fasta_file_directory=getattr(ns, "fasta_file_directory", None),
        score_threshold_file=getattr(ns, "score_threshold_file", None)
        or str(Path(paths_cfg.data) / "Thresholds"),
        library=getattr(ns, "library", None) or paths_cfg.hmms,
        result_files_directory=getattr(ns, "result_files_directory", None)
        or paths_cfg.results,
        database_directory=getattr(
            ns, "database_directory", None
        ),  # kann später via project gesetzt werden
        cores=int(getattr(ns, "cores", 4)),
        glob_report=getattr(ns, "glob_report", None),
        verbose=int(getattr(ns, "verbose", 1)),
    )

    cli_params = CliSearchParams(
        threshold_type=int(getattr(ns, "threshold_type", 1)),
        thrs_score=float(getattr(ns, "thrs_score", 50)),
        taxonomy_file=getattr(ns, "taxonomy_file", None),
        refseq_identity=int(getattr(ns, "refseq_identity", 90)),
        name=getattr(ns, "name", "project"),
        stage=int(getattr(ns, "stage", 0)),
        exit=int(getattr(ns, "exit", 10)),
    )

    cli_resources = CliResources(
        HMM_sets=list(getattr(ns, "HMM_sets", [])),
        clean_reports=bool(getattr(ns, "clean_reports", False)),
        individual_reports=bool(getattr(ns, "individual_reports", True)),
        max_seqs_per_genome=int(getattr(ns, "max_seqs_per_genome", 4)),
        bool_cross_check=bool(getattr(ns, "bool_cross_check", True)),
        optimized_cutoff_cross_check=bool(
            getattr(ns, "optimized_cutoff_cross_check", False)
        ),
    )

    cli_synteny = CliSynteny(
        patterns_file=getattr(ns, "patterns_file", None)
        or str(Path(paths_cfg.data) / "Patterns"),
        cooccurrence_file=getattr(ns, "cooccurrence_file", None)
        or str(Path(paths_cfg.data) / "Cooccurrence"),
        exclusion_singletons=getattr(ns, "exclusion_singletons", None)
        or str(Path(paths_cfg.data) / "Exclusion_singletons"),
        min_completeness=float(getattr(ns, "min_completeness", 0.5)),
        glob_chunks=int(getattr(ns, "glob_chunks", 5000)),
    )

    cli_info = CliInfo(
        stat_keywords=bool(getattr(ns, "stat_keywords", False)),
        stat_csb=bool(getattr(ns, "stat_csb", False)),
        stat_genomes=bool(getattr(ns, "stat_genomes", False)),
    )

    cli_csb = CliCsb(
        nucleotide_range=int(getattr(ns, "nucleotide_range", 3500)),
        insertions=int(getattr(ns, "insertions", 1)),
        occurence=int(getattr(ns, "occurence", 1)),
        min_csb_size=int(getattr(ns, "min_csb_size", 4)),
        max_csb_size=int(getattr(ns, "max_csb_size", 50)),
        max_domain_repeats=int(getattr(ns, "max_domain_repeats", 4)),
        jaccard=float(getattr(ns, "jaccard", 0.0)),
    )

    cli_flow = CliFlow(redo_taxonomy=bool(getattr(ns, "redo_taxonomy", False)))

    cli_limiter = CliLimiter(
        dataset_limit_lineage=getattr(ns, "dataset_limit_lineage", None),
        dataset_limit_taxon=getattr(ns, "dataset_limit_taxon", None),
        dataset_limit_proteins=getattr(ns, "dataset_limit_proteins", "0"),
        dataset_limit_keywords=getattr(ns, "dataset_limit_keywords", "0"),
        dataset_divide_sign=getattr(ns, "dataset_divide_sign", "."),
    )

    cli_ops = CliOperators(
        fetch_genomes=list(getattr(ns, "fetch_genomes", [])),
        fetch_proteins=list(getattr(ns, "fetch_proteins", [])),
        fetch_csbs=list(getattr(ns, "fetch_csbs", [])),
        fetch_keywords=list(getattr(ns, "fetch_keywords", [])),
        keywords_connector=getattr(ns, "keywords_connector", "OR"),
    )

    cli_process = CliProcess(
        merge_fasta=getattr(ns, "merge_fasta", None),
        filter_fasta=getattr(ns, "filter_fasta", None),
        concat_alignment=getattr(ns, "concat_alignment", None),
        add_taxonomy=getattr(ns, "add_taxonomy", None),
        add_genomic_context=getattr(ns, "add_genomic_context", None),
        create_type_range_dataset=getattr(ns, "create_type_range_dataset", None),
        create_gene_cluster_dataset=getattr(ns, "create_gene_cluster_dataset", None),
        gaps=bool(getattr(ns, "gaps", False)),
    )

    cfg = Config(
        paths=paths_cfg,
        cli_input=cli_input,
        cli_params=cli_params,
        cli_resources=cli_resources,
        cli_synteny=cli_synteny,
        cli_info=cli_info,
        cli_csb=cli_csb,
        cli_flow=cli_flow,
        cli_limiter=cli_limiter,
        cli_ops=cli_ops,
        cli_process=cli_process,
    )
    # Optional: Basiskonsistenz prüfen
    cfg.validate()
    return cfg


def _pretty_dump(cfg: Config) -> str:
    """Human-readable Flattening der wichtigsten Felder für den Printout."""
    rows: Dict[str, Any] = {
        # Pfade
        "paths.root": cfg.paths.root,
        "paths.results": cfg.paths.results,
        "paths.hmms": cfg.paths.hmms,
        # Input/Runtime
        "cli_input.fasta_file_directory": cfg.cli_input.fasta_file_directory,
        "cli_input.result_files_directory": cfg.cli_input.result_files_directory,
        "cli_input.library": cfg.cli_input.library,
        "cli_input.cores": cfg.cli_input.cores,
        # Search params
        "cli_params.threshold_type": cfg.cli_params.threshold_type,
        "cli_params.thrs_score": cfg.cli_params.thrs_score,
        "cli_params.name": cfg.cli_params.name,
        "cli_params.stage": cfg.cli_params.stage,
        # Synteny
        "cli_synteny.min_completeness": cfg.cli_synteny.min_completeness,
        # CSB
        "cli_csb.nucleotide_range": cfg.cli_csb.nucleotide_range,
        # Operators
        "cli_ops.keywords_connector": cfg.cli_ops.keywords_connector,
    }
    width = max(len(k) for k in rows)
    return "\n".join(f"{k:<{width}} : {v}" for k, v in rows.items())


def test_build_config_with_defaults_and_print(tmp_path, capsys):
    """
    Nutzt paths.py (refresh_paths) und parse.py (parse_cli) und baut config.Config.
    Prüft einige Defaults und gibt am Ende eine gut lesbare Übersicht aus.
    """
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)

    # Pfade auf unseren Fake-Baum setzen
    p.refresh_paths(root)

    # Minimal-CLI: nur -f (genomes). -r/-l lassen wir weg → Defaults greifen aus paths.
    argv = ["-f", str(root / "inputs"), "-c", "8", "-n", "unittest"]
    ns = cli_parse.parse_cli(argv)

    cfg = _fill_config_from_namespace(ns)

    # Assertions auf ein paar Kernwerte
    assert cfg.cli_input.fasta_file_directory == str(root / "inputs")
    assert cfg.cli_input.cores == 8
    # library kommt aus paths.hmms (Default)
    assert cfg.cli_input.library == str(root / "data" / "HMMlib")
    # results kommt aus paths.results (Default)
    assert cfg.cli_input.result_files_directory == str(root / "results")
    # Name übernommen
    assert cfg.cli_params.name == "unittest"

    # Human-readable Printout
    print("\nCONFIG OVERVIEW\n" + _pretty_dump(cfg))
    out = capsys.readouterr().out
    assert "CONFIG OVERVIEW" in out
    assert "cli_input.cores" in out and "8" in out


def test_build_config_with_custom_values_and_print(tmp_path, capsys):
    """
    Setzt benutzerdefinierte -r und weitere Optionen und prüft Mapping + Printout.
    """
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)
    p.refresh_paths(root)

    custom_results = root / "custom_results"
    custom_results.mkdir(parents=True, exist_ok=True)

    argv = [
        "-f",
        str(root / "inputs"),
        "-r",
        str(custom_results),
        "-c",
        "12",
        "-cut_type",
        "2",
        "-cut_score",
        "42",
        "-jaccard",
        "0.25",
        "-kc",
        "AND",
        "-n",
        "projB",
        "-s",
        "3",
    ]
    ns = cli_parse.parse_cli(argv)
    cfg = _fill_config_from_namespace(ns)

    # Assertions auf custom Werte
    assert cfg.cli_input.result_files_directory == str(custom_results)
    assert cfg.cli_input.cores == 12
    assert cfg.cli_params.threshold_type == 2
    assert cfg.cli_params.thrs_score == 42.0
    assert cfg.cli_csb.jaccard == 0.25
    assert cfg.cli_ops.keywords_connector == "AND"
    assert cfg.cli_params.name == "projB"
    assert cfg.cli_params.stage == 3

    print("\nCONFIG OVERVIEW (CUSTOM)\n" + _pretty_dump(cfg))
    out = capsys.readouterr().out
    assert "CONFIG OVERVIEW (CUSTOM)" in out
    assert "cli_params.thrs_score" in out and "42.0" in out

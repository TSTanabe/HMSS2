# tests/test_show_options.py
from __future__ import annotations

import sys
import importlib
from pathlib import Path


def test_show_parsed_options_uses_project_data(tmp_path, monkeypatch):
    """
    Zeigt und prüft, wie das Options-Objekt nach argparse aussieht,
    OHNE HMSSS_DATA_DIR (Fallback = PROJECT_ROOT/data).
    Ausgabe via print(...) – sichtbar mit: pytest -s tests/test_show_options.py
    """
    # --- src/ in den Pfad aufnehmen ---
    ROOT = Path(__file__).resolve().parents[1]
    SRC = ROOT / "src"
    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))

    # --- sicherstellen: KEINE ENV-Override ---
    monkeypatch.delenv("HMSSS_DATA_DIR", raising=False)

    # --- Module reloaden, damit der Fallback greift ---
    import hmsss.utils.paths as paths

    importlib.reload(paths)  # aktualisiert PROJECT_ROOT/DATA_ROOT je nach ENV
    import hmsss.cli.parse as parse

    importlib.reload(parse)

    # Projekt- & Datenwurzel nach Parser-/Paths-Logik
    project_root = paths.PROJECT_ROOT
    data_root = project_root / "data"

    # Diese Verzeichnisse MÜSSEN im Repo vorhanden sein:
    must_exist = [
        data_root / "RefSeqs",
        data_root / "Thresholds",
        data_root / "HMMlib",
        data_root / "Patterns",
        data_root / "Cooccurrence",
        data_root / "Exclusion_singletons",
    ]
    for p in must_exist:
        print(f"expect exists: {p}")
        assert p.exists(), (
            f"Fehlt: {p} (lege {data_root}/... an oder setze HMSSS_DATA_DIR)"
        )

    # --- Minimal-Argumente: -f (muss existieren), -db (darf neu sein) ---
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    db_file = tmp_path / "db.sqlite"

    opts = parse.parse_arguments(["-f", str(genomes_dir), "-db", str(db_file)])

    # --- Alles ausgeben, damit du es siehst ---
    print(f"PROJECT_ROOT: {project_root}")
    print(f"result_files_directory: {opts.result_files_directory}")
    print(f"reference_seq_dir:      {opts.reference_seq_dir}")
    print(f"score_threshold_file:   {opts.score_threshold_file}")
    print(f"library:                {opts.library}")
    print(f"patterns_file:          {opts.patterns_file}")
    print(f"cooccurrence_file:      {opts.cooccurrence_file}")
    print(f"exclusion_singletons:   {opts.exclusion_singletons}")

    # --- Prüfen, dass die Defaults auf PROJECT_ROOT/data zeigen ---
    assert Path(opts.reference_seq_dir).resolve() == (data_root / "RefSeqs").resolve()
    assert (
        Path(opts.score_threshold_file).resolve()
        == (data_root / "Thresholds").resolve()
    )
    assert Path(opts.library).resolve() == (data_root / "HMMlib").resolve()
    assert Path(opts.patterns_file).resolve() == (data_root / "Patterns").resolve()
    assert (
        Path(opts.cooccurrence_file).resolve() == (data_root / "Cooccurrence").resolve()
    )
    assert (
        Path(opts.exclusion_singletons).resolve()
        == (data_root / "Exclusion_singletons").resolve()
    )

    # kleine Sanity-Checks
    assert opts.fasta_file_directory == str(genomes_dir)
    assert opts.database_directory == str(db_file)

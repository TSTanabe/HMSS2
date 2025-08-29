# tests/test_parse_arguments.py
"""
Ausführlich kommentierter pytest für hmsss/cli/parse.py (HMSS2, Branch hades_diamond_v9).

Dieser Test:
- importiert hmsss.cli.parse direkt aus dem Repo (ohne Installation),
- emuliert die HMSS2-Projektstruktur im Temp-Verzeichnis,
- patched _project_root_from_this_file() und DATA_ROOT, damit Defaults auf die Teststruktur zeigen,
- testet unterschiedliche Argument-Kombinationen unabhängig voneinander.

Voraussetzungen:
- Projektlayout wie im Repo (src/hmsss/...) ist vorhanden.
- pytest wird aus dem Repo-Root ausgeführt (wo auch /src und /tests liegen).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Tuple, List

import pytest


# ------------------------------------------------------------
# Hilfs-Funktionen / -Fixtures
# ------------------------------------------------------------

def _add_src_to_syspath():
    """
    Sorgt dafür, dass 'from hmsss.cli import parse' funktioniert, ohne Paket zu installieren.
    Erwartung: Diese Testdatei liegt unter <REPO>/tests/, also ist <REPO>/src die Import-Root.
    """
    repo_root = Path(__file__).resolve().parents[1]
    src = repo_root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


@pytest.fixture(scope="session")
def parse_module():
    """
    Importiert das zu testende Modul hmsss.cli.parse einmal pro Test-Session.
    """
    _add_src_to_syspath()
    from hmsss.cli import parse  # type: ignore
    return parse


@pytest.fixture
def fake_project(tmp_path: Path) -> Dict[str, Path]:
    """
    Baut eine minimale HMSS2-Projektstruktur im temporären Verzeichnis.

    Struktur (nur das, was parse.py als Defaults nutzt):
    HMSS2/
      ├─ data/
      │   ├─ HMMlib/
      │   ├─ RefSeqs/
      │   ├─ Thresholds/
      │   ├─ Patterns/
      │   ├─ Cooccurrence/
      │   └─ Exclusion_singletons/
      ├─ results/
      └─ genomes/                (simuliertes Eingabeverzeichnis für -f)

    Zusätzlich werden ein paar Dummy-Dateien erzeugt, falls bestimmte Argumente
    (z. B. -taxonomy_info) einen existierenden Pfad verlangen.
    """
    root = tmp_path / "HMSS2"
    data = root / "data"
    results = root / "results"
    genomes = root / "genomes"

    # Verzeichnisse
    for p in [
        root,
        data / "HMMlib",
        data / "RefSeqs",
        data / "Thresholds",
        data / "Patterns",
        data / "Cooccurrence",
        data / "Exclusion_singletons",
        results,
        genomes,
    ]:
        p.mkdir(parents=True, exist_ok=True)

    # Beispiel-Datei für -taxonomy_info
    (tmp_path / "taxonomy.tsv").write_text("genomeID\tSuperkingdom\n", encoding="utf-8")

    return {
        "root": root,
        "data": data,
        "results": results,
        "genomes": genomes,
        "taxonomy_file": tmp_path / "taxonomy.tsv",
    }


@pytest.fixture
def patched_defaults(parse_module, fake_project, monkeypatch):
    """
    Robust:
    - Falls parse._project_root_from_this_file existiert -> patchen.
    - Andernfalls hmsss.cli.paths.refresh_paths(fake_root) aufrufen.
    - DATA_ROOT in parse (falls vorhanden) und hmsss.utils.paths patchen.
    """
    root = fake_project["root"]
    data = fake_project["data"]

    # 1) Versuche die alte Helper-Funktion in parse zu patchen (falls vorhanden)
    monkeypatch.setattr(parse_module, "_project_root_from_this_file", lambda: str(root), raising=False)

    # 2) Neue Pfad-API verwenden, falls vorhanden
    try:
        from hmsss.cli import paths as paths_mod  # neue Struktur
    except Exception:
        paths_mod = None

    if paths_mod is not None:
        # Projektwurzel & abgeleitete Konstanten (RESULTS, DATA, …) neu setzen
        paths_mod.refresh_paths(root)

    # 3) DATA_ROOT patchen (je nachdem, wo parse es herbezieht)
    #    a) direkt im parse-Modul (nur wenn dort vorhanden)
    monkeypatch.setattr(parse_module, "DATA_ROOT", data, raising=False)

    #    b) im alten utils.paths-Modul (falls parse "from hmsss.utils.paths import DATA_ROOT" verwendet)
    try:
        from hmsss.utils import paths as utils_paths_mod
        monkeypatch.setattr(utils_paths_mod, "DATA_ROOT", data, raising=False)
    except Exception:
        pass

    return {
        "root": root,
        "data": data,
        "results": fake_project["results"],
        "genomes": fake_project["genomes"],
    }


import sys
import inspect

def _call_parse(parse_module, args):
    """
    Ruft die Parser-Funktion robust auf:
    - parse_to_config(argv) bevorzugen, wenn vorhanden (nimmt argv als Liste)
    - sonst parse_arguments(argv) wenn die Signatur Parameter hat
    - sonst parse_arguments() (ohne Parameter) und sys.argv temporär setzen
    """
    # 1) parse_to_config(argv) vorhanden?
    if hasattr(parse_module, "parse_to_config"):
        return parse_module.parse_to_config(args)

    # 2) parse_arguments – Signatur prüfen
    if not hasattr(parse_module, "parse_arguments"):
        raise AttributeError("Neither parse_to_config nor parse_arguments found in hmsss.cli.parse")

    fn = parse_module.parse_arguments
    sig = inspect.signature(fn)
    if len(sig.parameters) == 0:
        old_argv = sys.argv[:]
        try:
            sys.argv = ["hmsss", *args]
            return fn()  # keine Args
        finally:
            sys.argv = old_argv
    else:
        return fn(args)  # nimmt argv-Liste


# ------------------------------------------------------------
# Tests
# ------------------------------------------------------------

def test_minimal_args_defaults(parse_module, patched_defaults, fake_project):
    """
    Minimalfall: nur -f <genomes_dir>. Erwartung:
    - Pfade werden absolut und aus Defaults abgeleitet.
    - 'new_project' ist True (Standard-Results-Pfad).
    - Bibliothek/Thresholds/Patterns/Cooccurrence/Exclusion werden aus DATA_ROOT ergänzt.
    """
    genomes = str(fake_project["genomes"])

    opts = _call_parse(parse_module, ["-f", genomes])

    # Eingabe-Verzeichnis
    assert Path(opts.fasta_file_directory) == Path(genomes)

    # Projektwurzel & Defaults
    assert Path(opts.location) == patched_defaults["root"]
    assert Path(opts.result_files_directory) == patched_defaults["results"]
    assert opts.new_project is True  # default results → "neues Projekt"

    data = patched_defaults["data"]
    assert Path(opts.reference_seq_dir) == data / "RefSeqs"
    assert Path(opts.library) == data / "HMMlib"
    assert Path(opts.score_threshold_file) == data / "Thresholds"
    assert Path(opts.patterns_file) == data / "Patterns"
    assert Path(opts.cooccurrence_file) == data / "Cooccurrence"
    assert Path(opts.exclusion_singletons) == data / "Exclusion_singletons"

    # Standard-Verbosity
    assert opts.verbose == 1  # aus parse.py Defaults


def test_custom_results_dir_disables_new_project(parse_module, patched_defaults, fake_project, tmp_path):
    """
    Wenn -r ein benutzerdefiniertes Results-Verzeichnis ist, dann ist new_project=False,
    und der Pfad wird absolut normalisiert.
    """
    genomes = str(fake_project["genomes"])
    user_results = tmp_path / "my_results"

    opts = _call_parse(parse_module, ["-f", genomes, "-r", str(user_results)])

    assert Path(opts.result_files_directory) == user_results.resolve()
    # Unterschied zum Default-Results → kein "neues Projekt"
    assert opts.new_project is False


def test_stage_101_when_process_or_fetch_args(parse_module, patched_defaults, fake_project, tmp_path):
    """
    Sobald "process/fetch"-relevante Argumente gesetzt werden, soll Stage=101 erkannt werden.
    Beispiel 1: Fetch-Operatoren (-fd) → erfordern -db.
    Beispiel 2: Processing-Operatoren (-filter_fasta) → erfordern ggf. -db, wenn add_taxonomy/... genutzt wird.
    """
    genomes = str(fake_project["genomes"])

    # --- Beispiel 1: Fetch (requires -db present, sonst sys.exit in parse) ---
    db_path = tmp_path / "database.db"  # muss nicht existieren, -db wird nicht auf Existenz geprüft
    opts_fetch = _call_parse(
        parse_module,
        ["-f", genomes, "-db", str(db_path), "-fd", "DsrA", "DsrB"]
    )
    assert opts_fetch.stage == 101
    assert opts_fetch.fetch is True
    # limiter bleibt False (keine -dll/-dlt/-dlp/-dlk)
    assert getattr(opts_fetch, "limiter", False) is False

    # --- Beispiel 2: Processing (hier: -filter_fasta) ---
    # -filter_fasta: drei Werte [FILE, MIN, MAX] → werden zu [str, int, int] gecastet
    out_file = tmp_path / "out.faa"
    opts_proc = _call_parse(
        parse_module,
        ["-f", genomes, "-filter_fasta", str(out_file), "100", "250"]
    )
    assert opts_proc.stage == 101
    assert opts_proc.process is True
    # Casting prüfen
    assert isinstance(opts_proc.filter_fasta[0], str)
    assert opts_proc.filter_fasta[1:] == [100, 250]


def test_stage_100_taxonomy_only_mode(parse_module, patched_defaults, fake_project, tmp_path):
    """
    Spezial-Stage 100: Wenn -taxonomy_info und -db gesetzt sind, aber KEIN -f,
    soll parse.py Stage=100 setzen (nur Taxonomy-Verarbeitung).
    """
    tax_file = str(fake_project["taxonomy_file"])
    db_path = tmp_path / "db.sqlite"

    opts = _call_parse(parse_module, ["-taxonomy_info", tax_file, "-db", str(db_path)])
    assert opts.stage == 100
    # In diesem Modus sollte kein fasta_file_directory gesetzt sein
    assert getattr(opts, "fasta_file_directory", None) in (None,)


def test_hmms_csv_and_whitespace_list_parsing(parse_module, patched_defaults, fake_project):
    """
    -hmms akzeptiert CSV ODER Whitespace. parse.py normalisiert das via _list_from_csv_or_repeat().
    """
    genomes = str(fake_project["genomes"])

    # CSV-Variante
    opts_csv = _call_parse(parse_module, ["-f", genomes, "-hmms", "A,B,C"])
    assert opts_csv.HMM_sets == ["A", "B", "C"]

    # Whitespace-Variante
    opts_ws = _call_parse(parse_module, ["-f", genomes, "-hmms", "A", "B", "C"])
    assert opts_ws.HMM_sets == ["A", "B", "C"]


def test_result_dir_is_absolute_even_if_relative_given(parse_module, patched_defaults, fake_project, monkeypatch):
    """
    Wenn -r relativ übergeben wird, wandelt parse.py den Pfad in einen absoluten um.
    """
    genomes = str(fake_project["genomes"])
    # simuliere ein relatives Ziel wie "rel_results"
    rel = "rel_results"

    # Arbeitsverzeichnis für den Test deterministisch setzen,
    # damit resolve() reproduzierbar ist.
    cwd = Path.cwd()
    try:
        os.chdir(str(fake_project["root"]))  # beliebiger stabiler Ort
        opts = _call_parse(parse_module, ["-f", genomes, "-r", rel])
        assert Path(opts.result_files_directory).is_absolute()
        # passt grob zur erwarteten Auflösung
        assert Path(opts.result_files_directory).name == "rel_results"
    finally:
        os.chdir(str(cwd))

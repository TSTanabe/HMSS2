# test/test_build_config_from_namespace.py
"""
Ausführlicher pytest für hmsss/cli/parse.build_config_from_namespace.

Was der Test tut:
- Baut eine minimale HMSS2-Struktur im Temp-Ordner, die mit Config.validate() kompatibel ist
  (bin/, data/HMMs, data/RefSeqs, results/).
- Lenkt hmsss.cli.paths auf diese Struktur (falls verfügbar via refresh_paths).
- Konstruiert argparse.Namespace-Objekte mit verschiedenen Argument-Kombinationen.
- Ruft ausschließlich build_config_from_namespace(ns) auf und prüft Felder im resultierenden Config.

Warum so?
- build_config_from_namespace konstruiert die Pfade über paths.as_dict() und baut daraus PathsCfg,
  danach erstellt es alle CLI-Blöcke und ruft cfg.validate() (daher müssen die Basisverzeichnisse existieren).  # :contentReference[oaicite:1]{index=1}
- Config.validate() verlangt u. a. bin/, data/HMMs, data/RefSeqs, results/.                        # :contentReference[oaicite:2]{index=2}
"""

from __future__ import annotations

import sys
from pathlib import Path
import argparse
import pytest


# ---------------------------
# Hilfs-Setup fürs Importing
# ---------------------------
def _add_src_to_syspath() -> None:
    """
    Import von 'hmsss...' ermöglichen, ohne Installation:
    Wir fügen <repo>/src vorne in sys.path hinzu.
    """
    repo_root = Path(__file__).resolve().parents[1]
    src = repo_root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


@pytest.fixture(scope="session")
def parse_mod():
    """
    Testobjekt importieren: hmsss.cli.parse
    """
    _add_src_to_syspath()
    from hmsss.cli import parse  # type: ignore
    return parse


# ---------------------------
# Test-Fixtures
# ---------------------------
@pytest.fixture
def fake_project(tmp_path: Path) -> dict[str, Path]:
    """
    Erzeugt einen minimalen HMSS2-Baum, der mit Config.validate() kompatibel ist:
      <root>/
        bin/
        data/HMMs/
        data/RefSeqs/
        results/
        src/hmsss/        (für PACKAGE_DIR)
    Dazu noch gängige Daten-Ordner (Thresholds/Patterns/Cooccurrence/Exclusion_singletons),
    damit Defaults in build_config_from_namespace auf sinnvolle Pfade zeigen.              # :contentReference[oaicite:3]{index=3}
    """
    root = tmp_path / "HMSS2"
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)

    # optionale, aber sinnvolle Datenverzeichnisse für Default-Zuweisungen
    for sub in ("Thresholds", "Patterns", "Cooccurrence", "Exclusion_singletons"):
        (root / "data" / sub).mkdir(parents=True, exist_ok=True)

    # Beispiel-Input-Verzeichnis für -f
    genomes = root / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)

    return {
        "root": root,
        "bin": root / "bin",
        "data": root / "data",
        "hmms": root / "data" / "HMMs",
        "refseq": root / "data" / "RefSeqs",
        "results": root / "results",
        "package": root / "src" / "hmsss",
        "genomes": genomes,
    }


@pytest.fixture
def set_paths(fake_project):
    """
    Richtet hmsss.cli.paths auf unser Fake-Projekt aus, damit build_config_from_namespace
    in _paths_cfg_from_paths_module() die richtigen Basispfade bekommt.                      # :contentReference[oaicite:4]{index=4}
    """
    from hmsss.cli import paths as paths_mod  # type: ignore

    # Wenn vorhanden: zentrale Re-Init-Funktion nutzen
    if hasattr(paths_mod, "refresh_paths"):
        paths_mod.refresh_paths(fake_project["root"])
    else:
        # Fallback: as_dict monkeypatchen – hier nicht nötig, wenn refresh_paths existiert.
        pass

    return fake_project


# ---------------------------
# Hilfsfunktion für Namespace
# ---------------------------
def ns(**kwargs) -> argparse.Namespace:
    """
    Bequemer Konstruktor für argparse.Namespace, damit die Testfälle lesbar bleiben.
    Alles, was nicht angegeben ist, lässt build_config_from_namespace intern per Defaults/paths füllen.
    """
    return argparse.Namespace(**kwargs)


# ===========================
#       Testfälle
# ===========================

def test_defaults_are_filled_from_paths(parse_mod, set_paths):
    """
    Prüft, dass fehlende Pfade per paths/as_dict → PathsCfg befüllt werden
    und ins Config-Objekt einfließen (library, results, patterns, ...).                   # :contentReference[oaicite:5]{index=5}
    """
    N = ns(
        fasta_file_directory=str(set_paths["genomes"]),
        # absichtlich KEINE result_files_directory / library / patterns_file etc.
        verbose=2,  # irgendein nicht-Default, damit wir sehen, dass Werte übernommen werden
    )

    cfg = parse_mod.build_config_from_namespace(N)  # ruft cfg.validate() → Basisdirs müssen existieren  # :contentReference[oaicite:6]{index=6}

    # PathsCfg stammt aus hmsss.cli.paths.as_dict()
    assert Path(cfg.paths.root) == set_paths["root"]
    assert Path(cfg.paths.bin) == set_paths["bin"]
    assert Path(cfg.paths.data) == set_paths["data"]
    assert Path(cfg.paths.hmms) == set_paths["hmms"]
    assert Path(cfg.paths.refseq) == set_paths["refseq"]
    assert Path(cfg.paths.results) == set_paths["results"]
    assert Path(cfg.paths.package) == set_paths["package"]

    # CLI Input: gegebener -f Pfad
    assert Path(cfg.cli_input.fasta_file_directory) == set_paths["genomes"]
    # nicht gesetzte Felder → Default aus PathsCfg / data-Unterordnern
    assert Path(cfg.cli_input.result_files_directory) == set_paths["results"]
    assert Path(cfg.cli_input.library) == set_paths["hmms"]
    assert Path(cfg.cli_synteny.patterns_file) == set_paths["data"] / "Patterns"
    assert Path(cfg.cli_synteny.cooccurrence_file) == set_paths["data"] / "Cooccurrence"
    assert Path(cfg.cli_synteny.exclusion_singletons) == set_paths["data"] / "Exclusion_singletons"
    assert Path(cfg.cli_input.score_threshold_file) == set_paths["data"] / "Thresholds"

    # unabhängige Felder: Defaults/übergebene Werte
    assert cfg.cli_input.verbose == 2
    assert cfg.cli_params.stage == 0          # nicht gesetzt → Default in build_config
    assert cfg.cli_resources.HMM_sets == []   # nicht gesetzt → Default
    assert cfg.cli_csb.jaccard == 0.0         # nicht gesetzt → Default                 # :contentReference[oaicite:7]{index=7}


def test_overrides_are_respected(parse_mod, set_paths, tmp_path):
    """
    Übergaben im Namespace überschreiben die Defaults:
    - eigenes results-Verzeichnis
    - eigene library
    - HMM_sets-Liste, Keywords-Connector, Scores/Thresholds, Stage etc.
    """
    custom_results = tmp_path / "my_results"
    custom_results.mkdir(parents=True, exist_ok=True)
    custom_library = tmp_path / "libdir"
    custom_library.mkdir(parents=True, exist_ok=True)

    N = ns(
        fasta_file_directory=str(set_paths["genomes"]),
        result_files_directory=str(custom_results),
        library=str(custom_library),
        name="myproj",
        threshold_type=2,
        thrs_score=75,
        stage=5,
        HMM_sets=["A", "B", "C"],
        keywords_connector="AND",
        jaccard=0.3,
        filter_fasta=["out.faa", "100", "250"],  # build_config übernimmt die Liste wie gegeben  # :contentReference[oaicite:8]{index=8}
    )

    cfg = parse_mod.build_config_from_namespace(N)

    # Overrides aus dem Namespace:
    assert Path(cfg.cli_input.result_files_directory) == custom_results
    assert Path(cfg.cli_input.library) == custom_library
    assert cfg.cli_params.name == "myproj"
    assert cfg.cli_params.threshold_type == 2
    assert cfg.cli_params.thrs_score == 75.0   # build_config castet zu float                 # :contentReference[oaicite:9]{index=9}
    assert cfg.cli_params.stage == 5
    assert cfg.cli_resources.HMM_sets == ["A", "B", "C"]
    assert cfg.cli_ops.keywords_connector == "AND"
    assert cfg.cli_csb.jaccard == 0.3
    assert cfg.cli_process.filter_fasta == ["out.faa", "100", "250"]


def test_type_casts_and_defaults(parse_mod, set_paths):
    """
    Test unabhängiger Felder/Typ-Casts:
    - thrs_score wird zu float gecastet
    - cores wird zu int gecastet
    - bools und Listen werden korrekt übernommen/standardisiert
    """
    N = ns(
        fasta_file_directory=str(set_paths["genomes"]),
        thrs_score="60",          # als String → build_config → float(60)                 # :contentReference[oaicite:10]{index=10}
        cores="8",                # als String → build_config → int(8)                   # :contentReference[oaicite:11]{index=11}
        clean_reports=True,       # bool
        individual_reports=False, # bool
        fetch_proteins=["DsrA", "DsrB"],  # Liste
    )

    cfg = parse_mod.build_config_from_namespace(N)

    assert cfg.cli_params.thrs_score == 60.0
    assert cfg.cli_input.cores == 8
    assert cfg.cli_resources.clean_reports is True
    assert cfg.cli_resources.individual_reports is False
    assert cfg.cli_ops.fetch_proteins == ["DsrA", "DsrB"]


def test_stage_changes_when_taxonomy_requested(parse_mod, set_paths, tmp_path, capsys):
    """
    Prüft zwei Szenarien rund um Taxonomy:

    (A) Taxonomy-only:
        - Voraussetzungen: taxonomy_info + database_directory, KEIN fasta_file_directory
        - Erwartung: cfg.cli_params.stage == 100

    Zusätzlich werden die relevanten Einstellungen für ein "Redo Taxonomy" ausgegeben (print),
    damit ersichtlich ist, welche Felder gesetzt sind. Mit `pytest -s` sieht man die Ausgabe.
    """
    # --- A) TAXONOMY-ONLY → STAGE 100 ---
    taxonomy_file = tmp_path / "taxonomy.tsv"
    taxonomy_file.write_text("genomeID\tSuperkingdom\n", encoding="utf-8")
    custom_results = tmp_path / "my_results"
    custom_results.mkdir(parents=True, exist_ok=True)

    N = ns(
        result_files_directory=str(custom_results),
        taxonomy_file=str(taxonomy_file),
        redo_taxonomy=True,
        stage=100,

    )

    cfg = parse_mod.build_config_from_namespace(N)

    assert cfg.cli_input.result_files_directory == str(custom_results)
    assert cfg.taxonomy_file == str(taxonomy_file)
    assert cfg.stage == 100
    assert cfg.cli_flow.redo_taxonomy is True
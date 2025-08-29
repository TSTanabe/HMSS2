# test/test_stage_taxonomy.py
"""
Tests für die Stage-Logik rund um Taxonomy:
- Taxonomy-only Modus (stage == 100), wenn nur Taxonomie gefordert wird.
- Redo Taxonomy als Processing-Schritt (stage == 101), wenn add_taxonomy aktiviert ist.

Die Routine gibt zusätzlich die relevanten Einstellungen für ein "Redo Taxonomy" aus.
"""

from __future__ import annotations

import sys
from pathlib import Path
import argparse
import pytest


# ------------------------------------------------------------
# Kleine Helpers/Fixtures (minimal, damit der Test unabhängig läuft)
# ------------------------------------------------------------
def _add_src_to_syspath() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    src = repo_root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


@pytest.fixture(scope="session")
def parse_mod():
    _add_src_to_syspath()
    from hmsss.cli import parse  # type: ignore
    return parse


@pytest.fixture
def fake_project(tmp_path: Path) -> dict[str, Path]:
    """
    Baut eine minimale Struktur, die mit Config.validate() kompatibel ist:
      <root>/
        bin/
        data/HMMs/
        data/RefSeqs/
        results/
        genomes/
        src/hmsss/
      plus: data/Thresholds, data/Patterns, data/Cooccurrence, data/Exclusion_singletons
    """
    root = tmp_path / "HMSS2"
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)
    for sub in ("Thresholds", "Patterns", "Cooccurrence", "Exclusion_singletons"):
        (root / "data" / sub).mkdir(parents=True, exist_ok=True)

    genomes = root / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)

    taxonomy_file = tmp_path / "taxonomy.tsv"
    taxonomy_file.write_text("genomeID\tSuperkingdom\n", encoding="utf-8")

    return {
        "root": root,
        "bin": root / "bin",
        "data": root / "data",
        "hmms": root / "data" / "HMMs",
        "refseq": root / "data" / "RefSeqs",
        "results": root / "results",
        "package": root / "src" / "hmsss",
        "genomes": genomes,
        "taxonomy": taxonomy_file,
    }


@pytest.fixture
def set_paths(fake_project):
    """
    Lenkt hmsss.cli.paths (falls vorhanden) auf unser Fake-Projekt,
    damit build_config_from_namespace() die korrekten Basis-Pfade bekommt.
    """
    try:
        from hmsss.cli import paths as paths_mod  # type: ignore
    except Exception:
        paths_mod = None
    if paths_mod is not None and hasattr(paths_mod, "refresh_paths"):
        paths_mod.refresh_paths(fake_project["root"])
    return fake_project


def ns(**kwargs) -> argparse.Namespace:
    """Bequemer Konstruktor für argparse.Namespace."""
    return argparse.Namespace(**kwargs)


# ------------------------------------------------------------
# Der eigentliche Test
# ------------------------------------------------------------
def test_stage_changes_when_taxonomy_requested(parse_mod, set_paths, tmp_path, capsys):
    """
    Prüft zwei unabhängige Szenarien:

    (A) Taxonomy-only:
        - Voraussetzungen: taxonomy_info + database_directory, KEIN fasta_file_directory
        - Erwartung: cfg.cli_params.stage == 100

    (B) Redo Taxonomy als Processing-Schritt:
        - Voraussetzungen: add_taxonomy=True, database_directory, (optional) fasta_file_directory
        - Erwartung: cfg.cli_params.stage == 101
        - Zusätzlich: es werden die relevanten Einstellungen "ausgedruckt"
          (für Nachvollziehbarkeit der benötigten Config-Felder).
    """
    # --- A) TAXONOMY-ONLY → STAGE 100 ---
    N_tax_only = ns(
        taxonomy_info=str(set_paths["taxonomy"]),
        database_directory=str(tmp_path / "db_tax_only.sqlite"),
        # kein fasta_file_directory → taxonomy-only Mode
    )
    cfg_tax_only = parse_mod.build_config_from_namespace(N_tax_only)
    assert cfg_tax_only.cli_params.stage == 100

    # --- B) REDO TAXONOMY (Processing) → STAGE 101 ---
    # Wir setzen add_taxonomy=True, geben eine DB an und (zur Sicherheit) ein -f-Verzeichnis.
    N_redo = ns(
        add_taxonomy=True,
        database_directory=str(tmp_path / "db_redo.sqlite"),
        fasta_file_directory=str(set_paths["genomes"]),
        # optional: weitere Prozess-Flags könnten hier ergänzt werden
    )
    cfg_redo = parse_mod.build_config_from_namespace(N_redo)
    assert cfg_redo.cli_params.stage == 101
    assert cfg_redo.cli_process.add_taxonomy is True

    # --- Relevante Einstellungen für "Redo Taxonomy" ausgeben ---
    print("\n=== Redo Taxonomy: relevante Einstellungen ===")
    print(f"stage: {cfg_redo.cli_params.stage}")
    print(f"add_taxonomy: {cfg_redo.cli_process.add_taxonomy}")
    print(f"database_directory: {cfg_redo.cli_input.database_directory}")
    # Falls du beim Redo zusätzlich eine Taxonomie-Tabelle/Datei einliest:
    print(f"taxonomy_info: {getattr(cfg_redo.cli_input, 'taxonomy_info', None)}")
    # Falls die Pipeline auch FASTA-Verzeichnis braucht (abhängig von deiner Logik):
    print(f"fasta_file_directory: {cfg_redo.cli_input.fasta_file_directory}")

    # Optional: die Prints validieren (so erscheinen sie auch ohne -s bei Fehlern im Output)
    captured = capsys.readouterr().out
    assert "Redo Taxonomy: relevante Einstellungen" in captured
    assert "add_taxonomy: True" in captured
    assert f"database_directory: {cfg_redo.cli_input.database_directory}" in captured

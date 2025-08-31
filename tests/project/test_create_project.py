# tests/cli/test_prepare_result_space.py
from __future__ import annotations

import sys
from pathlib import Path
import argparse
import pytest

# ganz oben
from pathlib import Path
from typing import (
    Any,
    Mapping,
    Iterable,
    Optional,
    Tuple,
)  # <- Optional, Tuple ergänzen


# ---------------- Import/Fixtures wie in deinen anderen Tests ----------------
def _add_src_to_syspath() -> None:
    repo_root = Path(__file__).resolve().parents[2]  # .../tests/cli -> repo-root
    src = repo_root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


@pytest.fixture(scope="session")
def parse_mod():
    _add_src_to_syspath()
    from hmsss.cli import parse  # type: ignore

    return parse


@pytest.fixture(scope="session")
def project_mod():
    _add_src_to_syspath()
    from hmsss.db import project  # type: ignore

    return project


@pytest.fixture
def fake_project(tmp_path: Path) -> dict[str, Path]:
    """
    Minimaler HMSS2-Basisbaum, der von deiner Config.validate() akzeptiert wird.
    """
    root = tmp_path / "HMSS2"
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)
    genomes = root / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    # optionale Ordner, die parse/build_config ggf. als Defaults setzt:
    for sub in ("Thresholds", "Patterns", "Cooccurrence", "Exclusion_singletons"):
        (root / "data" / sub).mkdir(parents=True, exist_ok=True)
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
    Lenkt hmsss.cli.paths auf den Fake-Baum.
    """
    try:
        from hmsss.cli import paths as paths_mod  # type: ignore
    except Exception:
        paths_mod = None
    if paths_mod is not None and hasattr(paths_mod, "refresh_paths"):
        paths_mod.refresh_paths(fake_project["root"])
    return fake_project


def ns(**kwargs) -> argparse.Namespace:
    return argparse.Namespace(**kwargs)


# ------------------------------- Hilfen --------------------------------------
def _db_path_from_cfg(cfg):
    """
    Hole die Datenbankpfad-Angabe aus Config.
    Priorität: cfg.project.database_directory -> cfg.cli_input.database_directory
    """
    db = None
    if getattr(cfg, "project", None) is not None:
        db = getattr(cfg.project, "database_directory", None)
    if not db:
        db = getattr(cfg.cli_input, "database_directory", None)
    return Path(db) if db else None


def _results_dir_from_cfg(cfg):
    """
    Hole das Ergebnis-Wurzelverzeichnis aus Config (das, worauf -r zeigt).
    Priorität: cfg.project.result_files_directory -> cfg.cli_input.result_files_directory
    """
    res = None
    if getattr(cfg, "project", None) is not None:
        res = getattr(cfg.project, "result_files_directory", None)
    if not res:
        res = getattr(cfg.cli_input, "result_files_directory", None)
    return Path(res) if res else None


def _mk_existing_project_tree(base: Path) -> dict[str, Path]:
    """
    Erzeuge eine Beispielstruktur eines bestehenden Projekts unter `base`.
    """
    base.mkdir(parents=True, exist_ok=True)
    db = base / "database.db"  # passe Name ggf. an (z.B. hmsss.sqlite)
    db.write_text("", encoding="utf-8")
    csb_dir = base / "CSB"
    csb_dir.mkdir(parents=True, exist_ok=True)
    (base / "global_report.hmmreport").write_text("", encoding="utf-8")
    (base / "global_trusted_hits.hmmreport").write_text("", encoding="utf-8")
    (base / "global_intermediate_hits.hmmreport").write_text("", encoding="utf-8")
    (csb_dir / "CSB_output.tsv").write_text("", encoding="utf-8")
    (csb_dir / "gene_clusters.tsv").write_text("", encoding="utf-8")
    return {
        "db": db,
        "csb_dir": csb_dir,
        "glob_report": base / "global_report.hmmreport",
        "glob_trusted": base / "global_trusted_hits.hmmreport",
        "glob_intermediate": base / "global_intermediate_hits.hmmreport",
        "csb_output": csb_dir / "CSB_output.tsv",
        "gene_clusters": csb_dir / "gene_clusters.tsv",
    }


# ================================ TESTS ======================================


def test_prepare_default_results_dir(parse_mod, project_mod, set_paths):
    """
    -r wird NICHT gesetzt → parse_to_config nutzt das Default-Results-Verzeichnis.
    Erwartung:
      * prepare_result_space richtet den Raum ein
      * results-Dir existiert
      * ggf. vorgesehene Subverzeichnisse werden angelegt
      * DB-Pfad liegt (später) unterhalb dieses results-Dirs
    """
    cfg = parse_mod.parse_to_config(["-f", str(set_paths["genomes"])])
    # Vorbedingung: results ist Default
    default_results = set_paths["results"]
    assert _results_dir_from_cfg(cfg) == default_results

    # Aktion
    project_mod.prepare_result_space(cfg)

    # Prüfungen
    res_dir = _results_dir_from_cfg(cfg)
    assert res_dir == default_results
    assert res_dir.is_dir()

    db = _db_path_from_cfg(cfg)
    # Beim „neuen“ Projekt existiert die DB evtl. noch nicht — aber
    # der vorgesehene Speicherort muss im results liegen:
    if db:
        assert str(db).startswith(str(res_dir))


def test_prepare_empty_custom_results_dir(parse_mod, project_mod, set_paths, tmp_path):
    """
    -r zeigt auf ein LEERES benutzerdefiniertes Verzeichnis.
    Erwartung:
      * prepare_result_space richtet den Raum ein
      * DB-Ziel (falls gesetzt) liegt unterhalb von -r
    """
    custom_res = tmp_path / "my_results"
    custom_res.mkdir(parents=True, exist_ok=True)
    cfg = parse_mod.parse_to_config(
        ["-f", str(set_paths["genomes"]), "-r", str(custom_res)]
    )

    project_mod.prepare_result_space(cfg)

    res_dir = _results_dir_from_cfg(cfg)
    assert res_dir == custom_res
    assert res_dir.is_dir()

    db = _db_path_from_cfg(cfg)
    if db:
        assert str(db).startswith(str(res_dir))


def test_prepare_results_with_db_in_subdirectory(
    parse_mod, project_mod, set_paths, tmp_path
):
    """
    -r zeigt auf ein Verzeichnis, in dessen *Unterordner* bereits eine Datenbank liegt.
    Erwartung:
      * prepare_result_space erkennt die DB unterhalb von -r
      * DB-Pfad im Config zeigt auf die gefundene Datei (unterhalb -r)
    """
    root = tmp_path / "existing_parent"
    root.mkdir(parents=True, exist_ok=True)
    sub = root / "run1"
    sub.mkdir(parents=True, exist_ok=True)
    dbfile = sub / "database.db"  # ggf. anpassen, falls deine DB anders heißt
    dbfile.write_text("", encoding="utf-8")

    cfg = parse_mod.parse_to_config(["-f", str(set_paths["genomes"]), "-r", str(root)])

    project_mod.prepare_result_space(cfg)

    res_dir = _results_dir_from_cfg(cfg)
    assert res_dir == root

    db = _db_path_from_cfg(cfg)
    assert db is not None
    assert db.exists()
    # die DB liegt in einem Unterordner von -r:
    assert str(db).startswith(str(root))
    assert db.name in ("database.db", "hmsss.sqlite", "hmsss.db")  # zur Not anpassen


def test_prepare_results_points_to_existing_project(
    parse_mod, project_mod, set_paths, tmp_path
):
    """
    -r zeigt direkt auf einen vorhandenen Projektordner (mit DB + Standard-Files).
    Erwartung:
      * prepare_result_space übernimmt diesen Ordner samt DB/File-Pfaden.
    """
    existing = tmp_path / "existing_project"
    files = _mk_existing_project_tree(existing)

    cfg = parse_mod.parse_to_config(
        ["-f", str(set_paths["genomes"]), "-r", str(existing)]
    )

    project_mod.prepare_result_space(cfg)

    res_dir = _results_dir_from_cfg(cfg)
    assert res_dir == existing

    db = _db_path_from_cfg(cfg)
    assert db == files["db"]

    # Wenn dein Config-Objekt ein `project`-Feld mit weiteren Pfaden trägt, prüfe es:
    if getattr(cfg, "project", None):
        assert Path(cfg.project.csb_directory) == files["csb_dir"]
        assert Path(cfg.project.csb_output_file) == files["csb_output"]
        assert Path(cfg.project.gene_clusters_file) == files["gene_clusters"]
        assert Path(cfg.project.glob_report) == files["glob_report"]
        assert Path(cfg.project.glob_trusted_hitreport) == files["glob_trusted"]
        assert (
            Path(cfg.project.glob_intermediate_hitreport) == files["glob_intermediate"]
        )

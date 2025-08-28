# tests/test_cli_paths.py
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

from hmsss.cli import paths as paths_mod


def _mk_fake_tree(root: Path) -> None:
    """
    Erzeuge die erwartete Ordnerstruktur unter root:
    bin, data/HMMs, data/RefSeqs, results, src/hmsss
    """
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)


def format_paths_human_readable(d: dict[str, str | Path]) -> str:
    """
    Gut lesbares Printout der Konstanten.
    """
    # Als Strings ausgeben, Keys auf gleiche Breite ausrichten
    items = [(k, str(v)) for k, v in d.items()]
    width = max(len(k) for k, _ in items) if items else 0
    lines = [f"{k:<{width}} : {v}" for k, v in items]
    return "\n".join(lines)


def test_refresh_paths_with_base_sets_expected_constants(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)

    # explizit setzen (keine Autodetektion), damit Test deterministisch ist
    paths_mod.refresh_paths(root)

    assert paths_mod.ROOT_DIR == root
    assert paths_mod.BIN_DIR == root / "bin"
    assert paths_mod.DATA_DIR == root / "data"
    assert paths_mod.HMMS_DIR == root / "data" / "HMMs"
    assert paths_mod.REFSEQ_DIR == root / "data" / "RefSeqs"
    assert paths_mod.RESULTS_DIR == root / "results"
    assert paths_mod.PACKAGE_DIR == root / "src" / "hmsss"

    # as_dict → Strings
    d_str = paths_mod.as_dict(str_paths=True)
    assert isinstance(d_str["ROOT_DIR"], str)
    assert d_str["ROOT_DIR"].endswith("HMSS2")

    # as_dict → Paths
    d_path = paths_mod.as_dict(str_paths=False)
    assert isinstance(d_path["ROOT_DIR"], Path)
    assert d_path["BIN_DIR"] == root / "bin"


def test_human_readable_printout(tmp_path: Path, capsys):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)
    paths_mod.refresh_paths(root)

    pretty = format_paths_human_readable(paths_mod.as_dict(str_paths=False))
    print(pretty)  # „gut lesbarer“ Output

    out = capsys.readouterr().out
    # Minimalprüfungen auf Inhalt/Format
    assert "ROOT_DIR" in out and str(root) in out
    assert "HMMS_DIR" in out and str(root / "data" / "HMMs") in out
    # Jede Zeile soll genau ein „ : “ enthalten
    assert all(" : " in line for line in out.strip().splitlines())


def test_as_dict_returns_expected_keys(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)
    paths_mod.refresh_paths(root)

    d = paths_mod.as_dict(str_paths=False)
    expected = {
        "ROOT_DIR",
        "BIN_DIR",
        "DATA_DIR",
        "HMMS_DIR",
        "REFSEQ_DIR",
        "RESULTS_DIR",
        "PACKAGE_DIR",
    }
    assert set(d.keys()) == expected


def test_compiled_binary_root_detection(monkeypatch, tmp_path: Path):
    """
    Simuliere: kompilierte Version — Executable liegt direkt unter HMSS2/.
    Erwartung: ROOT_DIR == Parent des Executables, restliche Struktur relativ dazu.
    """
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)

    # Fake-Binary anlegen (Name egal; Annahme: liegt direkt in HMSS2/)
    fake_exe = root / "HMSS2"
    fake_exe.write_bytes(b"\x00")  # Platzhalter

    # sys-Flags für „frozen“ Umgebung setzen
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "executable", str(fake_exe), raising=False)

    # Modul neu laden, damit _detect_root() erneut läuft
    reloaded = importlib.reload(paths_mod)

    # Assertions gegen das re-geladene Modul
    assert reloaded.ROOT_DIR == root
    assert reloaded.BIN_DIR == root / "bin"
    assert reloaded.DATA_DIR == root / "data"
    assert reloaded.HMMS_DIR == root / "data" / "HMMs"
    assert reloaded.REFSEQ_DIR == root / "data" / "RefSeqs"
    assert reloaded.RESULTS_DIR == root / "results"
    assert reloaded.PACKAGE_DIR == root / "src" / "hmsss"

    # Zur Sicherheit: Danach wieder auf einen definierten Zustand zurücksetzen,
    # damit andere Tests nicht vom Reload betroffen sind.
    reloaded.refresh_paths(root)

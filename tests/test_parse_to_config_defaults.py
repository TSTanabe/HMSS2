# tests/test_parse_to_config_defaults.py
from __future__ import annotations

from pathlib import Path

from hmsss.cli import paths as p
from hmsss.cli.parse import parse_to_config


def _mk_fake_tree(root: Path) -> None:
    """Erzeuge minimalen HMSS2-Baum (nur Verzeichnisse; Dateien sind für Defaults nicht nötig)."""
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)
    (root / "inputs").mkdir(parents=True, exist_ok=True)


def test_parse_to_config_sets_all_expected_defaults(tmp_path: Path):
    # Arrange: Projektstruktur simulieren & Pfade auf tmp setzen
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)
    p.refresh_paths(root)

    # Act: nur minimal erforderliche CLI-Args
    argv = ["-f", str(root / "inputs")]
    cfg = parse_to_config(argv)

    # Assert: Basis-Pfade (aus paths.py) sind korrekt übernommen
    assert cfg.paths.root    == str(root)
    assert cfg.paths.bin     == str(root / "bin")
    assert cfg.paths.data    == str(root / "data")
    assert cfg.paths.hmms    == str(root / "data" / "HMMs")
    assert cfg.paths.refseq  == str(root / "data" / "RefSeqs")
    assert cfg.paths.results == str(root / "results")
    assert cfg.paths.package == str(root / "src" / "hmsss")

    # Assert: CLI-Defaults nach _apply_runtime_defaults (weil nicht übergeben)
    assert cfg.cli_input.score_threshold_file == str(root / "data" / "Thresholds")
    assert cfg.cli_input.library              == str(root / "data" / "HMMlib")
    assert cfg.cli_synteny.patterns_file      == str(root / "data" / "Patterns")
    assert cfg.cli_synteny.cooccurrence_file  == str(root / "data" / "Cooccurrence")
    assert cfg.cli_synteny.exclusion_singletons == str(root / "data" / "Exclusion_singletons")
    assert cfg.cli_input.result_files_directory == str(root / "results")

    # Stage: ohne Prozess-/Fetch-/Redo-Argumente kein Override → 0 (Default)
    assert cfg.cli_params.stage == 0

    # Bonus: ein paar sinnvolle Defaults vom Parser
    assert cfg.cli_input.verbose == 1
    assert cfg.cli_input.cores >= 1  # Standard 4 oder was dein Parser setzt

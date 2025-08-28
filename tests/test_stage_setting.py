# tests/test_stage_setting.py
from __future__ import annotations

from pathlib import Path
import pytest

from hmsss.cli import paths as p
from hmsss.cli.parse import parse_to_config


def _mk_fake_tree(root: Path) -> None:
    """Erzeuge minimalen HMSS2-Baum plus Default-Dateien, damit file_path/dir_path-Checks nicht scheitern."""
    (root / "bin").mkdir(parents=True, exist_ok=True)
    (root / "data" / "HMMs").mkdir(parents=True, exist_ok=True)
    (root / "data" / "RefSeqs").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "src" / "hmsss").mkdir(parents=True, exist_ok=True)
    (root / "inputs").mkdir(parents=True, exist_ok=True)
    # Defaults, falls dein Parser sie als file_path prüft
    (root / "data" / "Thresholds").write_text("", encoding="utf-8")
    (root / "data" / "Patterns").write_text("", encoding="utf-8")
    (root / "data" / "Cooccurrence").write_text("", encoding="utf-8")
    (root / "data" / "Exclusion_singletons").write_text("", encoding="utf-8")


def _base_argv(root: Path) -> list[str]:
    """Minimales argv für alle Tests."""
    return ["-f", str(root / "inputs")]


def _cfg_for(root: Path, extra: list[str]) -> "Config":
    argv = _base_argv(root) + extra
    return parse_to_config(argv)


# -------------------------
# Kontrolltest ohne Prozess
# -------------------------
def test_stage_remains_when_no_processing(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root)
    p.refresh_paths(root)

    cfg = _cfg_for(root, ["-s", "3"])
    assert cfg.cli_params.stage == 3  # nichts triggert → kein 100er-Override


# -------------------------------------------
# Einzeltests: jedes Prozess-Argument für sich
# -------------------------------------------
def test_stage_100_on_merge_fasta(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    merge_dir = root / "merge_in"
    merge_dir.mkdir(parents=True, exist_ok=True)

    cfg = _cfg_for(root, ["-merge_fasta", str(merge_dir)])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_filter_fasta(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    out_file = root / "filtered.faa"
    cfg = _cfg_for(root, ["-filter_fasta", str(out_file), "10", "1000"])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_concat_alignment(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    aln_dir = root / "aln_dir"
    aln_dir.mkdir(parents=True, exist_ok=True)

    cfg = _cfg_for(root, ["-concat_alignment", str(aln_dir)])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_add_taxonomy_to_alignment(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    tax_file = root / "taxonomy.tsv"
    tax_file.write_text("id\ttax\n", encoding="utf-8")

    cfg = _cfg_for(root, ["-add_taxonomy_to_alignment", str(tax_file)])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_add_genomic_context(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    ctx_file = root / "context.faa"
    ctx_file.write_text(">seq\nM\n", encoding="utf-8")

    cfg = _cfg_for(root, ["-add_genomic_context", str(ctx_file)])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_create_type_range_dataset(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    out_ds = root / "type_range.tsv"
    cfg = _cfg_for(root, ["-create_type_range_dataset", str(out_ds)])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_create_gene_cluster_dataset(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    out_gc = root / "gene_clusters.tsv"
    cfg = _cfg_for(root, ["-create_gene_cluster_dataset", str(out_gc)])
    assert cfg.cli_params.stage == 100


# ------------------------------------
# Weitere Trigger laut deiner Vorgabe:
# - Fetch aus DB
# - redo_taxonomy
# ------------------------------------
def test_stage_100_on_db_fetch_keywords(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    cfg = _cfg_for(root, ["-fk", "sqr", "ddh"])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_db_fetch_lineage_taxon(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    cfg = _cfg_for(root, ["-fl", "Phylum", "-ft", "Proteobacteria"])
    assert cfg.cli_params.stage == 100


def test_stage_100_on_redo_taxonomy(tmp_path: Path):
    root = tmp_path / "HMSS2"
    _mk_fake_tree(root); p.refresh_paths(root)

    cfg = _cfg_for(root, ["-redo_taxonomy"])
    assert cfg.cli_params.stage == 100

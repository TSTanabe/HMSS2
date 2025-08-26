# tests/test_genomeID_files_dict.py
from pathlib import Path
import pytest
import hmsss.io.queue as file_finder  # dein Modul mit den Funktionen


def _ensure_myutil(monkeypatch):
    """Sorgt dafür, dass myUtil.get_genome_id wie gewünscht funktioniert."""
    import os

    def gid_fn(path: str) -> str:
        return os.path.basename(path).split(".")[0]

    FakeUtil = type("FakeUtil", (), {"get_genome_id": staticmethod(gid_fn)})
    monkeypatch.setattr(file_finder, "myUtil", FakeUtil, raising=True)


def test_mapping_recursive_and_extension(monkeypatch, tmp_path: Path):
    _ensure_myutil(monkeypatch)

    sub = tmp_path / "subdir"
    sub.mkdir()

    f1 = tmp_path / "ABC123_genomic.fna"
    f2 = sub / "XYZ789_something.fna"
    f1.write_text("ATGC")
    f2.write_text("TGCA")

    (tmp_path / "ignore.fa").write_text("AAAA")
    (sub / "ignore.fna.gz").write_text("BBBB")

    mapping = file_finder.get_genome_id_files_dict(str(tmp_path), "fna")

    assert "ABC123_genomic" in mapping
    assert "XYZ789_something" in mapping
    assert all(p.endswith(".fna") for p in mapping.values())


def test_genome_id_parsing_with_dots(monkeypatch, tmp_path: Path):
    _ensure_myutil(monkeypatch)

    f1 = tmp_path / "GCA_963695115.1_MFD09964.bin.1.60_genomic.fna"
    f2 = tmp_path / "GCF_000001405.40_GRCh38.p14_genomic.fna"
    f1.write_text("N")
    f2.write_text("N")

    mapping = file_finder.get_genome_id_files_dict(str(tmp_path), ".fna")

    assert "GCA_963695115" in mapping
    assert "GCF_000001405" in mapping


def test_empty_dir_returns_empty_dict(monkeypatch, tmp_path: Path):
    _ensure_myutil(monkeypatch)
    mapping = file_finder.get_genome_id_files_dict(str(tmp_path), ".fna")
    assert mapping == {}


def test_find_faa_gz(monkeypatch, tmp_path: Path):
    """Prüft, ob auch Endungen mit zwei Punkten (.faa.gz) gefunden werden."""
    _ensure_myutil(monkeypatch)

    f1 = tmp_path / "AAA111_predicted.faa.gz"
    f2 = tmp_path / "BBB222_predicted.faa.gz"
    f1.write_text("PROTSEQ1")
    f2.write_text("PROTSEQ2")

    f3 = tmp_path / "CCC333_predicted.faa"  # darf NICHT gefunden werden
    f3.write_text("PROTSEQ3")

    result = file_finder.get_all_files_with_extension(str(tmp_path), ".faa.gz")

    assert str(f1) in result
    assert str(f2) in result
    assert all(p.endswith(".faa.gz") for p in result)
    assert str(f3) not in result

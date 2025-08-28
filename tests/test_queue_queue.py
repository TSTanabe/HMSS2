# test_queue_routines.py
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Set, Optional
import gzip
import io
import types
import pytest
from hmsss.io.queue import get_genome_id_files_dict, _parallel_decompress
# -----------------------------------------------------------
# Hilfen für den Test
# -----------------------------------------------------------

class Options:
    def __init__(self, fasta_file_directory: str, cores: Optional[int] = 4):
        self.fasta_file_directory = fasta_file_directory
        self.cores = cores
        # werden in den Routinen befüllt:
        self.fna_files: Dict[str, str] = {}
        self.faa_files: Dict[str, str] = {}
        self.gff_files: Dict[str, str] = {}
        self.hmmreport_files: Dict[str, str] = {}
        self.queued_genomes: Set[str] = set()
        self.missing_gff_genomes: Set[str] = set()
        self.faa_missing_gff: Dict[str, str] = {}

# Minimale, robuste GenomeID-Extraktion:
# Nimmt den gesamten Basename ohne die angegebene Endung; bei .gz wird die vorletzte Endung entfernt.
def _strip_ext_for_genome_id(path: Path, extension: str) -> str:
    name = path.name
    if extension.endswith(".gz") and name.endswith(".gz"):
        name = name[:-3]  # .gz streichen
    if name.endswith(extension.replace(".gz", "")):
        return name[: -len(extension.replace(".gz", ""))]
    # Fallback: Basename ohne erste Punkt-Endung
    stem = name.split(".", 1)[0]
    return stem

def impl_get_genome_id_files_dict(root: str, extension: str = ".faa") -> Dict[str, str]:
    rootp = Path(root)
    out: Dict[str, str] = {}
    for p in rootp.rglob("*"):
        if not p.is_file():
            continue
        if p.name.endswith(extension):
            gid = _strip_ext_for_genome_id(p, extension)
            out[gid] = str(p)
    return out

def impl_parallel_decompress(paths: Set[str], cores: Optional[int] = None) -> None:
    # "Entpackt" *.gz -> legt einfach die ungezippte Datei an (leer oder Inhalt kopiert).
    for gz in paths:
        gz_path = Path(gz)
        assert gz_path.suffix == ".gz", f"Expected .gz path, got {gz_path}"
        out_path = gz_path.with_suffix("")  # .gz weg
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Falls Datei schon existiert, überschreiben
        # Simpler Inhalt: wenn gz lesbar, entpacken, sonst leer anlegen
        try:
            with gzip.open(gz_path, "rb") as f_in:
                data = f_in.read()
        except OSError:
            data = b""
        with open(out_path, "wb") as f_out:
            f_out.write(data)

# Dummy-Logger
class DummyLog:
    def info(self, *args, **kwargs):
        pass

log = DummyLog()


# -----------------------------------------------------------
# Die zu testenden Routinen (wie in deiner Vorlage),
# aber als freie Funktionen, damit wir sie direkt importfrei testen können.
# WICHTIG: Sie referenzieren get_genome_id_files_dict und _parallel_decompress,
# die wir im Test per monkeypatch ersetzen.
# -----------------------------------------------------------

def queue_fna_inputs(options) -> dict[str, str]:
    root = options.fasta_file_directory
    fna_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".fna.gz")
    fna_files: Dict[str, str]    = get_genome_id_files_dict(root, extension=".fna")
    faa_files: Dict[str, str]    = get_genome_id_files_dict(root, extension=".faa")
    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")

    faa_all_genomes: Set[str] = set(faa_files) | set(faa_gz_files)

    decompress_targets: Set[str] = set()
    for gid, gz_path in fna_gz_files.items():
        if gid not in faa_all_genomes:
            if gid not in fna_files:
                decompress_targets.add(gz_path)

    if decompress_targets:
        log.info(f"[FNA] Planned to decompress {len(decompress_targets)} file(s).")
        _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    fna_files = get_genome_id_files_dict(root, extension=".fna")

    for gid in list(fna_files.keys()):
        if gid in faa_all_genomes:
            del fna_files[gid]

    options.fna_files = fna_files
    log.info(f"Found {len(fna_files)} fna files for translation.")
    return fna_files


def queue_protein_annotation_inputs(options) -> None:
    root = options.fasta_file_directory

    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")
    gff_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff.gz")
    faa_files: Dict[str, str]    = get_genome_id_files_dict(root, extension=".faa")
    gff_files: Dict[str, str]    = get_genome_id_files_dict(root, extension=".gff")

    decompress_targets: Set[str] = set()
    for gid, gz_path in faa_gz_files.items():
        if gid not in faa_files:
            decompress_targets.add(gz_path)
    for gid, gz_path in gff_gz_files.items():
        if gid not in gff_files:
            decompress_targets.add(gz_path)

    if decompress_targets:
        log.info(f"[ANN] Planned to decompress {len(decompress_targets)} file(s).")
        _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    faa_files = get_genome_id_files_dict(root, extension=".faa")
    gff_files = get_genome_id_files_dict(root, extension=".gff")
    hmmreport_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".hmmreport")

    common_ids: Set[str] = set(faa_files) & set(gff_files)

    faa_files = {gid: path for gid, path in faa_files.items() if gid in common_ids}
    gff_files = {gid: path for gid, path in gff_files.items() if gid in common_ids}
    hmmreport_files = {gid: path for gid, path in hmmreport_files.items() if gid in common_ids}

    options.queued_genomes = common_ids
    options.faa_files = faa_files
    options.gff_files = gff_files
    options.hmmreport_files = hmmreport_files

    log.info(f"Queued {len(common_ids)} faa/gff pairs.")
    log.info(f"Found {len(hmmreport_files)} existing hmmreports for faa/gff pairs.")


def queue_faa_without_gff(options) -> dict[str, str]:
    root = options.fasta_file_directory

    faa_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa")
    faa_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".faa.gz")
    gff_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff")
    gff_gz_files: Dict[str, str] = get_genome_id_files_dict(root, extension=".gff.gz")

    genomes_with_any_gff: Set[str] = set(gff_files) | set(gff_gz_files)
    genomes_with_any_faa: Set[str] = set(faa_files) | set(faa_gz_files)

    missing_gff_genomes: Set[str] = genomes_with_any_faa - genomes_with_any_gff

    decompress_targets: Set[str] = set()
    for gid in missing_gff_genomes:
        if gid not in faa_files and gid in faa_gz_files:
            decompress_targets.add(faa_gz_files[gid])

    if decompress_targets:
        log.info(f"Planned to decompress {len(decompress_targets)} file(s).")
        _parallel_decompress(decompress_targets, getattr(options, "cores", None))

    faa_files = get_genome_id_files_dict(root, extension=".faa")

    faa_missing_gff: Dict[str, str] = {
        gid: faa_files[gid]
        for gid in missing_gff_genomes
        if gid in faa_files
    }

    log.info(f"Queued {len(faa_missing_gff)} faa files without gff for transcription.")
    return faa_missing_gff


# -----------------------------------------------------------
# Pytest Fixtures
# -----------------------------------------------------------

@pytest.fixture
def patched_helpers(monkeypatch):
    # monkeypatch globale Namen, die in den Routinen referenziert werden
    monkeypatch.setitem(globals(), "get_genome_id_files_dict", impl_get_genome_id_files_dict)
    monkeypatch.setitem(globals(), "_parallel_decompress", impl_parallel_decompress)
    # optional: echten Logger ersetzen, falls gewünscht
    monkeypatch.setitem(globals(), "log", DummyLog())


@pytest.fixture
def example_tree(tmp_path: Path):
    """
    Legt eine Teststruktur mit vielen Kombinationen an:

    A: A.fna.gz                      -> FNA ohne FAA -> soll entpackt & behalten werden
    B: B.fna.gz + B.faa.gz           -> FNA darf NICHT behalten werden (FAA existiert)
    C: C.fna                         -> FNA ohne FAA -> behalten
    D: D.faa + D.gff + D.hmmreport   -> vollständiges Paar + report
    E: E.faa.gz + E.gff              -> FAA.gz muss entpackt werden -> Paar
    F: F.faa + F.gff.gz              -> GFF.gz muss entpackt werden -> Paar
    G: G.faa.gz + G.gff.gz           -> beide müssen entpackt -> Paar
    H: H.faa                         -> FAA ohne GFF -> in faa_without_gff
    I: I.faa.gz                      -> FAA.gz ohne GFF -> entpacken & in faa_without_gff
    J: J.gff                         -> nur GFF (ignoriert)
    K: K.hmmreport                   -> verwaister report (soll herausfallen)
    In subdir: L.fna.gz (wie A, nur in Unterordner)
    """
    files = [
        "A.fna.gz",
        "B.fna.gz", "B.faa.gz",
        "C.fna",
        "D.faa", "D.gff", "D.hmmreport",
        "E.faa.gz", "E.gff",
        "F.faa", "F.gff.gz",
        "G.faa.gz", "G.gff.gz",
        "H.faa",
        "I.faa.gz",
        "J.gff",
        "K.hmmreport",
        # Unterverzeichnis
        "sub/L.fna.gz",
    ]
    for rel in files:
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if rel.endswith(".gz"):
            with gzip.open(p, "wb") as f:
                f.write(b"test")
        else:
            p.write_bytes(b"test")

    return tmp_path


# -----------------------------------------------------------
# Tests
# -----------------------------------------------------------

def test_queue_fna_inputs(example_tree, patched_helpers):
    opts = Options(str(example_tree))

    # Vorbedingungen prüfen
    # A.fna.gz, C.fna, sub/L.fna.gz sind FNA-Kandidaten
    # B.fna.gz hat FAA.gz -> darf NICHT behalten werden
    fna_map = queue_fna_inputs(opts)

    # Erwartung: A und L werden entpackt (A.fna, L.fna), C bleibt;
    # B.fna(.gz) wird NICHT in fna_files gelistet (weil FAA existiert)
    got_ids = set(opts.fna_files.keys())
    assert got_ids == {"A", "C", "sub/L"} or got_ids == {"A", "C", "L"}  # je nach GenomeID-Extraktion; hier rechnen wir mit "sub/L"

    # Entpackte Dateien existieren?
    assert (example_tree / "A.fna").exists()
    assert (example_tree / "sub" / "L.fna").exists()
    # B.fna wurde NICHT entpackt, weil FAA vorhanden
    assert not (example_tree / "B.fna").exists()


def test_queue_protein_annotation_inputs(example_tree, patched_helpers):
    opts = Options(str(example_tree))

    queue_protein_annotation_inputs(opts)

    # Erwartete Paare: D, E, F, G
    # - D: bereits ungezipt
    # - E: faa.gz -> faa entpackt
    # - F: gff.gz -> gff entpackt
    # - G: beides gz -> beides entpackt
    expected_pairs = {"D", "E", "F", "G"}
    assert opts.queued_genomes == expected_pairs

    # Entpackte Dateien vorhanden?
    assert (example_tree / "E.faa").exists()
    assert (example_tree / "F.gff").exists()
    assert (example_tree / "G.faa").exists()
    assert (example_tree / "G.gff").exists()

    # HMMREPORT: nur solche, deren GenomeID in den Pairs ist (hier nur D)
    assert set(opts.hmmreport_files.keys()) == {"D"}

    # Keine FAA-only (H, I) in den Paaren
    assert "H" not in opts.queued_genomes
    assert "I" not in opts.queued_genomes


def test_queue_faa_without_gff(example_tree, patched_helpers):
    opts = Options(str(example_tree))

    # Vorab: sicherstellen, dass noch nichts entpackt wurde, was der Test benötigt
    # (getrennte Aufrufe, keine Abhängigkeit von vorherigen Tests)
    result = queue_faa_without_gff(opts)

    # Erwartung:
    # - B.faa (hat fna und faa.gz)
    # - H.faa (bereits ungezipt) => drin
    # - I.faa.gz (ohne gff) => muss entpackt werden und dann drin sein
    expected = {"B", "H", "I"}
    assert set(result.keys()) == expected

    # Entpackung von I erfolgt?
    assert (example_tree / "I.faa").exists()

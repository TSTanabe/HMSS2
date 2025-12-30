"""
pyrodigal_translation_mp.py

Parallel translation stage for genome assemblies:
- Inputs: dict[genome_id -> fna_path (.fna/.fa/.fasta optionally .gz)]
- Outputs (written by workers next to the input file):
    <basename>.faa.gz
    <basename>.gff.gz

Best-of-both-worlds writing:
- Stream records directly into a TEMP gzip file (low RAM, no giant IPC strings)
- On success: os.replace(temp, final) (atomic commit)
- On failure: temp files are removed; final files remain untouched

Design:
- Workers do everything (materialize .gz -> train -> find_genes -> write files -> gzip).
- Main process only schedules work and logs progress (like search_pyhmmer).
- Worker processes are capped (default 16) to avoid hammering shared filesystems.

Notes:
- Uses pyrodigal in genome mode (meta=False) with explicit training.
- Uses materialize_single_next_to_input for .gz inputs.
"""

from __future__ import annotations

import gzip
import os
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

from hmsss.core.logging import get_logger
from hmsss.fasta_preparation import materialize

logger = get_logger(__name__)

# -------------------------
# Worker globals (per process)
# -------------------------
_G_GENE_FINDER = None
_G_MIN_TRAIN_BP = 20000
_G_SKIP_IF_UPTODATE = True


def _init_pyrodigal_worker(min_train_bp: int, skip_if_uptodate: bool):
    global _G_GENE_FINDER, _G_MIN_TRAIN_BP, _G_SKIP_IF_UPTODATE
    import pyrodigal

    _G_GENE_FINDER = pyrodigal.GeneFinder(meta=False)
    _G_MIN_TRAIN_BP = min_train_bp
    _G_SKIP_IF_UPTODATE = skip_if_uptodate


# -------------------------
# Helpers
# -------------------------
def _iter_fasta_records(fasta_path: str) -> Iterator[tuple[str, str]]:
    """Minimal FASTA reader; yields (contig_id, nucleotide_sequence)."""
    header: Optional[str] = None
    seq_chunks: list[str] = []
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def _output_paths_next_to_input(fna_path: str) -> tuple[Path, Path]:
    """
    Determine output paths next to input.
    /x/genome.fna      -> /x/genome.faa.gz, /x/genome.gff.gz
    /x/genome.fna.gz   -> /x/genome.faa.gz, /x/genome.gff.gz
    /x/genome.fa.gz    -> /x/genome.faa.gz, /x/genome.gff.gz
    """
    p = Path(fna_path)
    name = p.name[:-3] if p.name.endswith(".gz") else p.name  # strip .gz
    stem = Path(name).stem  # strip last extension (.fna/.fa/.fasta)
    base = p.with_name(stem)
    return Path(str(base) + ".faa.gz"), Path(str(base) + ".gff.gz")


def _is_up_to_date(fna_path: str, faa_gz: Path, gff_gz: Path) -> bool:
    """
    Skip if outputs exist and are newer than input.
    For .gz input, compare against the .gz mtime (fine for pipeline usage).
    """
    try:
        if not faa_gz.exists() or not gff_gz.exists():
            return False
        in_mtime = Path(fna_path).stat().st_mtime
        return faa_gz.stat().st_mtime >= in_mtime and gff_gz.stat().st_mtime >= in_mtime
    except OSError:
        return False


def _select_training_seqs(
    contigs: list[tuple[str, str]], min_train_bp: int
) -> list[str]:
    """Pick longest contigs until reaching min_train_bp."""
    contigs_sorted = sorted(contigs, key=lambda x: len(x[1]), reverse=True)
    train: list[str] = []
    total = 0
    for _cid, seq in contigs_sorted:
        if not seq:
            continue
        train.append(seq)
        total += len(seq)
        if total >= min_train_bp:
            break
    return train


@dataclass(frozen=True)
class _TmpOutputs:
    faa_tmp: Path
    gff_tmp: Path


def _make_tmp_paths_for(final_faa_gz: Path, final_gff_gz: Path) -> _TmpOutputs:
    """
    Create unique temp paths in the same directory as final outputs.
    Same directory is required for atomic os.replace on most filesystems.
    """
    out_dir = final_faa_gz.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Unique filenames (mkstemp gives us a real file path; we close the fd and re-open with gzip)
    faa_fd, faa_tmp = tempfile.mkstemp(
        prefix=final_faa_gz.name + ".tmp_", dir=str(out_dir)
    )
    gff_fd, gff_tmp = tempfile.mkstemp(
        prefix=final_gff_gz.name + ".tmp_", dir=str(out_dir)
    )
    os.close(faa_fd)
    os.close(gff_fd)

    return _TmpOutputs(faa_tmp=Path(faa_tmp), gff_tmp=Path(gff_tmp))


def _cleanup_tmp(paths: _TmpOutputs) -> None:
    for p in (paths.faa_tmp, paths.gff_tmp):
        try:
            if p.exists():
                p.unlink()
        except OSError:
            pass


def _commit_tmp(paths: _TmpOutputs, final_faa_gz: Path, final_gff_gz: Path) -> None:
    """
    Atomically replace final outputs with the completed temp files.
    """
    os.replace(paths.faa_tmp, final_faa_gz)
    os.replace(paths.gff_tmp, final_gff_gz)


def _write_wrapped_fasta_record(fh, header: str, seq: str, width: int = 60) -> None:
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


# -------------------------
# Worker status (small IPC)
# -------------------------
@dataclass(frozen=True)
class TranslationStatus:
    gid: str
    fna_path: str
    faa_gz: str
    gff_gz: str
    status: str  # "ok" | "skipped" | "error"
    message: str = ""


def pyrodigal_translate_and_write_worker(
    item: tuple[str, str],
    *,
    min_train_bp: int = 20000,
    translation_table: int = 11,
    force_nonsd: bool = False,
    start_weight: float = 4.35,
    skip_if_uptodate: bool = True,
) -> TranslationStatus:
    """
    Worker:
      - materialize .fna.gz if needed
      - train (meta=False requires train)
      - find genes for each contig
      - stream-write FAA+GFF into temp gzip files
      - atomic commit (os.replace) to final .faa.gz/.gff.gz paths
    Returns a small status object for main-process progress/logging.
    """
    gid, fna_path = item

    global _G_GENE_FINDER
    if _G_GENE_FINDER is None:
        raise RuntimeError(
            "Worker not initialized: GeneFinder is None (missing initializer?)."
        )

    final_faa_gz, final_gff_gz = _output_paths_next_to_input(fna_path)

    if skip_if_uptodate and _is_up_to_date(fna_path, final_faa_gz, final_gff_gz):
        return TranslationStatus(
            gid=gid,
            fna_path=fna_path,
            faa_gz=str(final_faa_gz),
            gff_gz=str(final_gff_gz),
            status="skipped",
            message="up-to-date",
        )

    tmp_paths = _make_tmp_paths_for(final_faa_gz, final_gff_gz)

    try:
        with materialize.materialize_single_next_to_input(fna_path) as plain_fna:
            contigs = [(cid, dna) for cid, dna in _iter_fasta_records(plain_fna) if dna]
            if not contigs:
                # Still write minimal valid outputs
                with gzip.open(tmp_paths.faa_tmp, "wt") as faa_fh:
                    pass
                with gzip.open(tmp_paths.gff_tmp, "wt") as gff_fh:
                    gff_fh.write("##gff-version 3\n")

                _commit_tmp(tmp_paths, final_faa_gz, final_gff_gz)
                return TranslationStatus(
                    gid=gid,
                    fna_path=fna_path,
                    faa_gz=str(final_faa_gz),
                    gff_gz=str(final_gff_gz),
                    status="ok",
                    message="no contigs; wrote empty outputs",
                )

            train_seqs = _select_training_seqs(contigs, min_train_bp=min_train_bp)
            if not train_seqs:
                raise RuntimeError("could not select training contigs")

            # meta=False -> must train
            _G_GENE_FINDER.train(
                train_seqs[0],
                *train_seqs[1:],
                translation_table=translation_table,
                force_nonsd=force_nonsd,
                start_weight=start_weight,
            )

            # Stream outputs into temp gzip files (low RAM)
            with (
                gzip.open(tmp_paths.faa_tmp, "wt") as faa_fh,
                gzip.open(tmp_paths.gff_tmp, "wt") as gff_fh,
            ):
                gff_fh.write("##gff-version 3\n")
                gene_counter = 1
                for contig_id, dna in contigs:
                    genes = _G_GENE_FINDER.find_genes(dna)

                    for i, gene in enumerate(genes, start=1):
                        aa = gene.translate()
                        if not aa:
                            continue

                        prot_id = f"{gid}___{gene_counter}"
                        gene_counter += 1

                        # FAA
                        _write_wrapped_fasta_record(faa_fh, prot_id, aa, width=60)

                        # GFF3 (1-based inclusive coordinates)
                        start = int(gene.begin) + 1
                        end = int(gene.end)
                        strand = "+" if gene.strand == 1 else "-"
                        attrs = f"ID={prot_id};locus_tag={prot_id}"
                        gff_fh.write(
                            f"{contig_id}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n"
                        )

        # Atomic commit only after successful completion
        _commit_tmp(tmp_paths, final_faa_gz, final_gff_gz)

        return TranslationStatus(
            gid=gid,
            fna_path=fna_path,
            faa_gz=str(final_faa_gz),
            gff_gz=str(final_gff_gz),
            status="ok",
        )

    except Exception as e:
        _cleanup_tmp(tmp_paths)
        return TranslationStatus(
            gid=gid,
            fna_path=fna_path,
            faa_gz=str(final_faa_gz),
            gff_gz=str(final_gff_gz),
            status="error",
            message=str(e),
        )


# -------------------------
# Public API: parallel stage
# -------------------------
def parallel_pyrodigal_translation(
    fna_files: Dict[str, str],
    cores: int,
    *,
    max_workers: int = 16,
    min_train_bp: int = 20000,
    chunksize: int = 1,
    skip_if_uptodate: bool = True,
) -> None:
    """
    Parallel translation stage.

    Requirements satisfied:
      - outputs (.faa.gz/.gff.gz) written next to input .fna
      - input .fna.gz materialized temporarily via materialize_single_next_to_input
      - limit worker processes to 16
      - workers write outputs themselves (no large IPC payloads)
      - main process logs writer progress like search_pyhmmer

    CPU assignment:
      - pyrodigal/prodigal is effectively single-core per genome in this setup
      - so: 1 CPU per worker; cap workers at min(cores, 16, n_genomes)
    """
    if not fna_files:
        logger.info("No FNA inputs queued for pyrodigal translation.")
        return

    items = list(fna_files.items())
    n_genomes = len(items)

    worker_processes = min(max_workers, max(1, cores), n_genomes)

    logger.info(
        "Starting pyrodigal translation: %d genomes, parallel=%d (cap=%d, cores=%d)",
        n_genomes,
        worker_processes,
        max_workers,
        cores,
    )

    genomes_done = 0
    log_step = max(1, n_genomes // 100)

    errors = 0
    skipped = 0

    with Pool(
        processes=worker_processes,
        initializer=_init_pyrodigal_worker,
        initargs=(min_train_bp, skip_if_uptodate),
    ) as pool:
        for status in pool.imap_unordered(
            pyrodigal_translate_and_write_worker, items, chunksize=chunksize
        ):
            genomes_done += 1

            if status.status == "error":
                errors += 1
                logger.error(
                    "[pyrodigal] %s failed for %s: %s",
                    status.gid,
                    status.fna_path,
                    status.message,
                )
            elif status.status == "skipped":
                skipped += 1

            if (genomes_done % log_step == 0) or (genomes_done == n_genomes):
                pct = (genomes_done * 100) // max(1, n_genomes)
                logger.info(
                    f"[ORF prediction] {genomes_done}/{n_genomes} ({pct}%) genomes processed"
                )

    logger.info(
        "Finished pyrodigal translation stage: %d total, %d skipped, %d errors.",
        n_genomes,
        skipped,
        errors,
    )
    if errors:
        raise RuntimeError(
            f"pyrodigal translation: {errors} genomes failed (see logs)."
        )

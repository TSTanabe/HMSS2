#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import random
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterator, Tuple, List, Optional


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Streaming FASTA iterator.
    Yields: (header_without_>, sequence_string)
    """
    header: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def reservoir_sample_fasta(
    path: Path, k: int, seed: int | None = None
) -> List[Tuple[str, str]]:
    """
    Reservoir sampling over FASTA records: O(k) memory, 1 pass.
    Returns up to k (header, seq) tuples.
    """
    rng = random.Random(seed)
    reservoir: List[Tuple[str, str]] = []
    for i, rec in enumerate(iter_fasta(path), start=1):
        if len(reservoir) < k:
            reservoir.append(rec)
        else:
            j = rng.randint(1, i)
            if j <= k:
                reservoir[j - 1] = rec
    return reservoir


def make_reads_from_refseq(
    ref_faa: Path,
    out_faa: Path,
    *,
    n_reads: int = 1000,
    min_len: int = 30,
    max_len: int = 60,
    seed: int | None = None,
    reservoir_k: int = 5000,
) -> int:
    """
    Generates n_reads random fragments (length min_len..max_len) from sequences in ref_faa.
    If insufficient distinct sequences exist, samples with replacement automatically.
    """
    rng = random.Random(seed)

    # Sample a manageable pool of reference sequences
    pool = reservoir_sample_fasta(ref_faa, k=reservoir_k, seed=seed)

    if not pool:
        raise RuntimeError(f"No FASTA records found in: {ref_faa}")

    # Filter for sequences that can produce required fragment lengths
    viable = [(h, s) for (h, s) in pool if len(s) >= min_len]
    if not viable:
        raise RuntimeError(
            f"None of the sampled sequences in {ref_faa} are >= {min_len} aa. "
            f"Try increasing reservoir_k or check the input."
        )

    out_faa.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    attempts = 0
    max_attempts = n_reads * 50  # safety

    with out_faa.open("w", encoding="utf-8") as out:
        while written < n_reads and attempts < max_attempts:
            attempts += 1
            header, seq = rng.choice(viable)

            L = rng.randint(min_len, max_len)
            if len(seq) < L:
                # If seq shorter than drawn L but still >= min_len, clip L down
                L = rng.randint(min_len, min(max_len, len(seq)))

            if len(seq) < min_len:
                continue

            start = rng.randint(0, len(seq) - L)
            frag = seq[start : start + L]

            # deterministic-ish id: read_000001|src=...
            rid = f"read_{written + 1:06d}|src={header.split()[0]}|pos={start}|len={L}"
            out.write(f">{rid}\n{frag}\n")
            written += 1

    if written < n_reads:
        raise RuntimeError(
            f"Could only generate {written}/{n_reads} reads from {ref_faa}. "
            f"Check min_len/max_len and reference sequence lengths."
        )

    return written


def kill_process_group(p: subprocess.Popen) -> None:
    """Kill entire process group (graftM + children)"""
    try:
        os.killpg(p.pid, signal.SIGTERM)
    except Exception:
        pass
    time.sleep(1.0)
    try:
        os.killpg(p.pid, signal.SIGKILL)
    except Exception:
        pass


def peak_rss_of_tree_bytes(ps_proc) -> int:
    """
    Sum RSS of a process + its children recursively.
    Requires psutil.
    """
    total = 0
    try:
        total += ps_proc.memory_info().rss
        for child in ps_proc.children(recursive=True):
            try:
                total += child.memory_info().rss
            except Exception:
                continue
    except Exception:
        return 0
    return total


def run_graftm_and_measure_peak_rss(
    *,
    gpkg_dir: Path,
    reads_faa: Path,
    out_dir: Path,
    threads: int = 8,
    poll_s: float = 0.2,
) -> tuple[int, float]:
    """
    Runs 'graftM graft' and measures peak RSS (GB) of graftM process tree.
    """
    try:
        import psutil  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            "psutil is required for peak RSS measurement. Install via: pip install psutil"
        ) from e

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "graftM",
        "graft",
        "--graftm_package",
        str(gpkg_dir),
        "--forward",
        str(reads_faa),
        "--output_directory",
        str(out_dir),
        "--threads",
        str(threads),
        "--force",
    ]

    # Start in its own process group so we can cleanly kill everything
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        preexec_fn=os.setsid,
    )

    ps = psutil.Process(p.pid)
    peak = 0

    try:
        while p.poll() is None:
            rss = peak_rss_of_tree_bytes(ps)
            if rss > peak:
                peak = rss
            time.sleep(poll_s)

        stdout, stderr = p.communicate()

        # Echo graftM logs (optional but useful)
        if stdout:
            sys.stdout.write(stdout)
        if stderr:
            sys.stderr.write(stderr)

        peak_gb = peak / 1024**3
        return p.returncode, peak_gb

    except KeyboardInterrupt:
        kill_process_group(p)
        raise
    except Exception:
        kill_process_group(p)
        raise


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Profile peak RAM for a GraftM gpkg by placing 1000 random AA fragments via graftM graft (pplacer)."
    )
    ap.add_argument("gpkg", help="Path to a *.gpkg directory")
    ap.add_argument("--threads", type=int, default=8, help="Threads for graftM")
    ap.add_argument(
        "--n-reads", type=int, default=1000, help="Number of reads to generate"
    )
    ap.add_argument("--min-len", type=int, default=30, help="Min AA length per read")
    ap.add_argument("--max-len", type=int, default=60, help="Max AA length per read")
    ap.add_argument("--seed", type=int, default=1, help="Random seed")
    ap.add_argument(
        "--reservoir-k",
        type=int,
        default=5000,
        help="Reservoir size for refseq sampling",
    )
    ap.add_argument(
        "--keep", action="store_true", help="Keep temp directory (for debugging)"
    )
    args = ap.parse_args(argv)

    gpkg_dir = Path(args.gpkg).resolve()
    if not gpkg_dir.is_dir():
        ap.error(f"Not a directory: {gpkg_dir}")

    ref_faa = gpkg_dir / "refseq_database.faa"
    if not ref_faa.is_file():
        ap.error(f"Missing {ref_faa} (expected inside gpkg directory)")

    tmp_root = Path(tempfile.mkdtemp(prefix=f"gpkg_ram_profile_{gpkg_dir.stem}_"))
    try:
        reads_faa = tmp_root / "reads_1000.faa"
        out_dir = tmp_root / "graftm_out"

        n = make_reads_from_refseq(
            ref_faa,
            reads_faa,
            n_reads=args.n_reads,
            min_len=args.min_len,
            max_len=args.max_len,
            seed=args.seed,
            reservoir_k=args.reservoir_k,
        )
        print(f"[INFO] Generated {n} reads: {reads_faa}")

        rc, peak_gb = run_graftm_and_measure_peak_rss(
            gpkg_dir=gpkg_dir,
            reads_faa=reads_faa,
            out_dir=out_dir,
            threads=args.threads,
        )

        pkg_name = gpkg_dir.name  # z.B. "sqdg.gpkg" oder "dsr_core.gpkg"
        print(
            f"[RESULT] gpkg={pkg_name} threads={args.threads} reads={args.n_reads} "
            f"peak_rss_gb={peak_gb:.2f} exit_code={rc}"
        )

        return rc

    finally:
        if args.keep:
            print(f"[DEBUG] Keeping temp dir: {tmp_root}")
        else:
            # best-effort cleanup
            try:
                import shutil

                shutil.rmtree(tmp_root, ignore_errors=True)
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main())

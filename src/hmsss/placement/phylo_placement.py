#!/usr/bin/env python3
"""Align protein queries to an existing reference MSA and place using EPA-ng.

Requires: hmmalign (HMMER 3) and epa-ng on PATH. No Python dependencies.
The seed MSA must be the one used to build the supplied HMM.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
from pathlib import Path

LOG = logging.getLogger("hmg.epa")


def fasta_records(path: Path):
    """Stream FASTA records as (first-word ID, sequence)."""
    name, parts = None, []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Empty FASTA identifier in {path}")
                parts = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)


def load_alignment(path: Path) -> dict[str, str]:
    records = {}
    length = None
    for name, seq in fasta_records(path):
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier: {name}")
        if not seq:
            raise ValueError(f"Empty sequence: {name}")
        if length is None:
            length = len(seq)
        if len(seq) != length:
            raise ValueError(f"Unaligned FASTA: inconsistent lengths in {path}")
        records[name] = seq
    if not records:
        raise ValueError(f"No sequences in {path}")
    return records


def write_query_alignment(seed_msa: Path, combined_msa: Path, output: Path) -> int:
    """Remove hmmalign insertion columns and verify original MSA is unchanged."""
    seed = load_alignment(seed_msa)
    combined = load_alignment(combined_msa)
    seed_ids = set(seed)
    if not seed_ids.issubset(combined):
        raise ValueError(
            f"hmmalign output missing reference IDs: {sorted(seed_ids - set(combined))[:5]}"
        )
    queries = {k: v for k, v in combined.items() if k not in seed_ids}
    if not queries:
        raise ValueError("hmmalign produced no aligned queries")
    aln_length = len(next(iter(combined.values())))
    keep = [
        i
        for i in range(aln_length)
        if any(combined[rid][i] not in "-." for rid in seed)
    ]
    expected = len(next(iter(seed.values())))
    if len(keep) != expected:
        raise ValueError(
            "Cannot reconstruct original seed alignment columns from hmmalign output. "
            "Check that the HMM was built from exactly this seed alignment."
        )
    for rid, original in seed.items():
        reconstructed = (
            "".join(combined[rid][i] for i in keep).upper().replace(".", "-")
        )
        if reconstructed != original.upper().replace(".", "-"):
            raise ValueError(
                f"Seed alignment does not match HMM map at reference {rid}"
            )
    with output.open("w", encoding="utf-8") as handle:
        for qid, seq in queries.items():
            aligned = "".join(seq[i] for i in keep).upper().replace(".", "-")
            if not aligned.strip("-"):
                raise ValueError(f"Query has no aligned residues: {qid}")
            handle.write(f">{qid}\n{aligned}\n")
    return len(queries)


def run(
    tree: Path,
    seed: Path,
    hmm: Path,
    queries: Path,
    model: str,
    output_dir: Path,
    threads: int,
    overwrite: bool = False,
) -> Path:
    for path in (tree, seed, hmm, queries):
        if not path.is_file():
            raise FileNotFoundError(path)
    if threads < 1:
        raise ValueError("--threads must be >= 1")
    for executable in ("hmmalign", "epa-ng"):
        if not shutil.which(executable):
            raise RuntimeError(f"{executable} is not available on PATH")
    seed_records = load_alignment(seed)
    query_ids = set()
    n_queries = 0
    for qid, _ in fasta_records(queries):
        if qid in seed_records or qid in query_ids:
            raise ValueError(f"Query ID duplicated or overlaps references: {qid}")
        query_ids.add(qid)
        n_queries += 1
    if not n_queries:
        raise ValueError("No query sequences")
    output_dir.mkdir(parents=True, exist_ok=True)
    combined = output_dir / "combined_alignment.fasta"
    aligned_queries = output_dir / "aligned_queries.fasta"
    jplace = output_dir / "epa_result.jplace"
    if jplace.exists() and not overwrite:
        raise FileExistsError(f"{jplace} already exists (use --overwrite to replace)")
    if overwrite and jplace.exists():
        jplace.unlink()

    hmm_cmd = [
        "hmmalign",
        "--amino",
        "--outformat",
        "afa",
        "--mapali",
        str(seed),
        str(hmm),
        str(queries),
    ]
    LOG.info("Aligning queries: %s", " ".join(hmm_cmd))
    with combined.open("w", encoding="utf-8") as handle:
        subprocess.run(hmm_cmd, stdout=handle, check=True)
    aligned_count = write_query_alignment(seed, combined, aligned_queries)
    LOG.info("Aligned %d query sequences", aligned_count)
    cmd = [
        "epa-ng",
        "--tree",
        str(tree),
        "--ref-msa",
        str(seed),
        "--query",
        str(aligned_queries),
        "--model",
        model,
        "--threads",
        str(threads),
        "--outdir",
        str(output_dir),
    ]
    LOG.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    if not jplace.is_file() or jplace.stat().st_size == 0:
        raise RuntimeError("EPA-ng did not create a nonempty epa_result.jplace")
    with jplace.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not all(key in data for key in ("tree", "fields", "placements")):
        raise ValueError("Output JPLACE is missing required fields")
    LOG.info("Produced %s (%d placement groups)", jplace, len(data["placements"]))
    return jplace


def main():
    # Runs the epa-ng placement and returns the output json file
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--tree", required=True, type=Path, help="Reference tree (Newick)"
    )
    parser.add_argument(
        "--seed-alignment",
        required=True,
        type=Path,
        help="Reference aligned FASTA used to construct HMM",
    )
    parser.add_argument(
        "--hmm", required=True, type=Path, help="Profile HMM built from reference MSA"
    )
    parser.add_argument(
        "--queries",
        required=True,
        type=Path,
        help="Unaligned amino-acid queries in FASTA",
    )
    parser.add_argument(
        "--model",
        required=True,
        help="EPA-ng model descriptor or RAxML-NG bestModel file",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    output = run(
        args.tree.resolve(),
        args.seed_alignment.resolve(),
        args.hmm.resolve(),
        args.queries.resolve(),
        args.model,
        args.output_dir.resolve(),
        args.threads,
        args.overwrite,
    )
    print(output)


if __name__ == "__main__":
    main()

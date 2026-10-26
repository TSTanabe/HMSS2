#!/usr/bin/env python3
"""Align protein queries to an existing reference MSA and place using EPA-ng.

Requires: hmmalign (HMMER 3) and epa-ng on PATH. No Python dependencies.
The seed MSA must be the one used to build the supplied HMM.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import shutil
import subprocess
from pathlib import Path

import pyhmmer

from hmsss.db import database

logger = logging.getLogger("hmg.epa")


def find_matching_placement_packages(placement_root, intermediate_files):
    """Yield reference packages and matching intermediate FASTA files."""

    for package_dir in sorted(Path(placement_root).iterdir()):
        annotations = package_dir / "annotations.tsv"

        if not package_dir.is_dir() or not annotations.is_file():
            continue
        try:
            with annotations.open() as f:
                protein_types = {
                    row[1] for row in csv.reader(f, delimiter="\t") if len(row) >= 2
                }
        except (OSError, UnicodeError, csv.Error) as e:
            logger.warning("Skipping package %s: %s", package_dir.name, e)
            continue

        matched_types = protein_types & intermediate_files.keys()
        if not matched_types:
            continue

        faa_files = []
        for protein_type in sorted(matched_types):
            files = intermediate_files[protein_type]
            faa_files.extend([files] if isinstance(files, (str, Path)) else files)

        yield package_dir, list(dict.fromkeys(map(Path, faa_files)))


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


def merge_placement_fastas(faa_files, output_file):
    """Merge valid FASTA files, skipping faulty files and duplicate IDs."""

    output_file = Path(output_file)
    seen_ids = set()

    with output_file.open("w") as out:
        for faa_file in faa_files:
            try:
                faa_file = Path(faa_file)

                for seq_id, sequence in fasta_records(faa_file):
                    if seq_id in seen_ids:
                        continue

                    seen_ids.add(seq_id)
                    out.write(f">{seq_id}\n{sequence}\n")

            except (OSError, ValueError, UnicodeError) as e:
                logger.warning("Skipping faulty FASTA file %s: %s", faa_file, e)
                continue

    return output_file


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
            if not aligned.strip("-."):
                logger.warning("Skipping query without aligned residues: %s", qid)
                continue
            handle.write(f">{qid}\n{aligned}\n")
    return len(queries)


def load_placement_package(package_dir: str | Path) -> dict[str, Path]:
    """Load and validate EPA-ng reference files using SHA-256 checksums."""

    package_dir = Path(package_dir).resolve()
    manifest_file = package_dir / "manifest.json"

    with manifest_file.open() as f:
        manifest = json.load(f)

    required = {
        "tree": "tree",
        "seed": "trimmed_alignment",
        "hmm": "hmm",
    }

    package_files = {}

    for argument, role in required.items():
        filename = manifest["roles"][role]
        path = (package_dir / filename).resolve()

        if not path.is_relative_to(package_dir):
            raise ValueError(f"Invalid package path: {filename}")

        if not path.is_file():
            raise FileNotFoundError(path)

        # Verify SHA-256 checksum
        expected_sha256 = manifest["files"][filename]["sha256"]

        digest = hashlib.sha256()
        with path.open("rb") as f:
            while chunk := f.read(1024 * 1024):
                digest.update(chunk)

        actual_sha256 = digest.hexdigest()

        if actual_sha256.lower() != expected_sha256.lower():
            raise ValueError(
                f"SHA-256 checksum mismatch for {path}\n"
                f"Expected: {expected_sha256}\n"
                f"Actual:   {actual_sha256}"
            )

        package_files[argument] = path

    logger.info(
        "Reference package '%s': all SHA-256 checksums verified",
        package_dir.name,
    )

    return package_files


def _fasta_records_from_handle(handle):
    """Yield FASTA records from a binary or text file handle."""

    name, parts = None, []

    for line in handle:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        line = line.strip()

        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)

            name = line[1:].split()[0]
            parts = []

        elif line:
            parts.append(line)

    if name is not None:
        yield name, "".join(parts)


def align_queries_pyhmmer(
    hmm: Path,
    seed: Path,
    queries: Path,
    output: Path,
    batch_size: int = 256,
    threads: int = 1,
) -> int:
    """Align protein queries to an HMM in batches using PyHMMER."""

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.plan7.HMMFile(str(hmm)) as hmm_file:
        model = hmm_file.read()

    if model is None:
        raise ValueError(f"No HMM found: {hmm}")

    # HMM and reference alignment must have identical column counts
    with pyhmmer.easel.MSAFile(str(seed), format="afa") as msa_file:
        reference = msa_file.read()

    aligned_count = 0

    with (
        pyhmmer.easel.SequenceFile(
            str(queries), digital=True, alphabet=alphabet
        ) as seq_file,
        output.open("w") as out,
    ):
        while True:
            batch = seq_file.read_block(sequences=batch_size)

            if not batch:
                break

            # Align queries to the HMM
            msa = pyhmmer.hmmer.hmmalign(
                model,
                batch,
                trim=False,
                all_consensus_cols=True,
            )

            # A2M marks insertion residues with lowercase letters
            buffer = io.BytesIO()
            msa.write(buffer, "a2m")
            buffer.seek(0)

            for query_id, sequence in _fasta_records_from_handle(buffer):
                # Retain only HMM match-state columns
                aligned = "".join(c for c in sequence if c.isupper() or c == "-")

                if len(aligned) != model.M or not aligned.strip("-"):
                    logger.warning("Skipping invalid alignment: %s", query_id)
                    continue

                out.write(f">{query_id}\n{aligned}\n")
                aligned_count += 1

    logger.info("PyHMMER aligned %d query sequences", aligned_count)
    return aligned_count


def run_epa_ng(
    tree: Path,
    seed: Path,
    hmm: Path,
    queries: Path,
    model: str,
    output_dir: Path,
    threads: int,
) -> Path | None:
    """Align queries with PyHMMER and perform phylogenetic placement."""

    if not shutil.which("epa-ng"):
        raise RuntimeError("EPA-ng executable not found")

    output_dir.mkdir(parents=True, exist_ok=True)

    aligned_queries = output_dir / "aligned_queries.fasta"
    jplace = output_dir / "epa_result.jplace"

    # Remove previous GAPPA results
    for existing_jplace in output_dir.glob("*.jplace"):
        existing_jplace.unlink()

    # 1. Align queries with PyHMMER
    aligned_count = align_queries_pyhmmer(
        hmm=hmm,
        seed=seed,
        queries=queries,
        output=aligned_queries,
        batch_size=2056,
        threads=4,
    )

    if not aligned_count:
        logger.warning("No alignable queries for EPA-ng")
        return None

    # 2. EPA-ng placement
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

    logger.info("Running EPA-ng with %d queries", aligned_count)
    subprocess.run(cmd, check=True)

    # 3. Validate JPLACE output
    if not jplace.is_file() or jplace.stat().st_size == 0:
        raise RuntimeError(f"EPA-ng produced no valid output: {jplace}")

    with jplace.open() as f:
        data = json.load(f)

    if not {"tree", "fields", "placements"} <= data.keys():
        raise ValueError(f"Invalid JPLACE format: {jplace}")

    logger.info(
        "EPA-ng completed: %s (%d placement groups)",
        jplace,
        len(data["placements"]),
    )

    return jplace


def assign_protein_clades(
    jplace: Path, annotations: Path, output_dir: Path, threshold: float = 0.75
) -> Path:
    """Assign placements to protein-type clades using GAPPA."""

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gappa",
        "edit",
        "extract",
        "--jplace-path",
        str(jplace.resolve()),
        "--clade-list-file",
        str(annotations.resolve()),
        "--threshold",
        str(threshold),
    ]

    logger.info("Assigning protein clades: %s", jplace)
    subprocess.run(cmd, cwd=output_dir, check=True)

    return output_dir


def parse_function_assignments(
    function_dir: str | Path,
    database_path: str,
) -> int:
    """Parse GAPPA clade assignments and update domains or selection comments."""

    assignments = {}

    for file in sorted(Path(function_dir).glob("*.jplace")):
        domain_type = file.stem

        with file.open() as f:
            data = json.load(f)

        for placement in data["placements"]:
            names = placement.get("n", [])
            names += [entry[0] for entry in placement.get("nm", [])]

            for name in names:
                domain = database.parse_domain_header(name)

                if domain in assignments and assignments[domain] != domain_type:
                    logger.warning(
                        f"Conflicting assignments for {domain}: "
                        f"{assignments[domain]} vs {domain_type}"
                    )
                    continue
                assignments[domain] = domain_type

    grouped = {}
    for domain, domain_type in assignments.items():
        grouped.setdefault(domain_type, set()).add(domain)

    updated = 0

    for domain_type, domain_hits in grouped.items():
        if domain_type in {"basal", "uncertain"}:
            database.update_selection_comments(
                database=database_path,
                domain_hits=domain_hits,
                comment=f"{domain_type}_placement",
            )
        else:
            updated += database.update_domains(
                database=database_path,
                domain_hits=domain_hits,
                new_domain=domain_type,
            )

    logger.debug("Updated %d domain assignments", updated)
    return updated


def assign_taxonomy(
    jplace: Path,
    taxonomy_file: Path,
    output_dir: Path,
    threads: int = 1,
) -> Path:
    """Assign taxonomy to EPA-ng placements using GAPPA."""

    if not jplace.is_file():
        raise FileNotFoundError(jplace)
    if not taxonomy_file.is_file():
        raise FileNotFoundError(taxonomy_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gappa",
        "examine",
        "assign",
        "--jplace-path",
        str(jplace.resolve()),
        "--taxon-file",
        str(taxonomy_file.resolve()),
        "--per-query-results",
        "--best-hit",
        "--threads",
        str(threads),
        "--out-dir",
        str(output_dir.resolve()),
    ]

    logger.info("Assigning taxonomy: %s", jplace)
    subprocess.run(cmd, check=True)

    return output_dir


def parse_taxonomy_assignments(taxonomy_dir: str | Path) -> dict[str, str]:
    """Parse GAPPA per-query taxonomy assignments."""

    domain_lineages = {}

    for file in Path(taxonomy_dir).glob("*per_query.tsv"):
        with file.open(newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                domain_lineages[row["name"]] = row["taxopath"]

    return domain_lineages


def remove_processed_files(intermediate_files, processed_files):
    """Remove successfully processed FASTA files from intermediate_files."""

    processed_files = set(map(Path, processed_files))

    for protein_type, files in list(intermediate_files.items()):
        if isinstance(files, (str, Path)):
            if Path(files) in processed_files:
                del intermediate_files[protein_type]
        else:
            remaining = [file for file in files if Path(file) not in processed_files]

            if remaining:
                intermediate_files[protein_type] = remaining
            else:
                del intermediate_files[protein_type]


def epa_ng_placement(config, intermediate_files):
    """
    Group intermediate hits by reference package.

    Each package is processed only once, even if multiple
    protein types match.
    """

    packages = find_matching_placement_packages(
        config.paths.placement, intermediate_files
    )

    for package_dir, faa_files in packages:
        try:
            output_dir = (
                Path(config.phylogenetic_placement_directory) / package_dir.name
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            queries = merge_placement_fastas(
                faa_files, output_dir / "placement_candidates.faa"
            )
            package = load_placement_package(package_dir)

            logger.info("Running phylogenetic placement for %s", package_dir.name)

            # EPA-ng placement
            jplace = run_epa_ng(
                tree=package["tree"],
                seed=package["seed"],
                hmm=package["hmm"],
                queries=queries,
                model="LG",
                output_dir=output_dir,
                threads=config.cores,
            )
            if jplace is None:
                logger.warning(
                    "Skipping package %s: no valid placements", package_dir.name
                )
                continue

            # Functional classification
            function_dir = assign_protein_clades(
                jplace=jplace,
                annotations=package_dir / "annotations.tsv",
                output_dir=output_dir / "function",
                threshold=0.75,
            )

            parse_function_assignments(
                function_dir=function_dir,
                database_path=config.database_directory,
            )

            # Taxonomic classification (metagenomes only)
            taxonomy_dir = assign_taxonomy(
                jplace=jplace,
                taxonomy_file=package_dir / "reference_taxonomy.tsv",
                output_dir=output_dir / "taxonomy",
                threads=config.cores,
            )

            domain_lineages = parse_taxonomy_assignments(taxonomy_dir)

            database.insert_domain_placements(
                database=config.database_directory,
                domain_lineages=domain_lineages,
            )

            # Remove successfully processed files from the dictionary
            remove_processed_files(intermediate_files, faa_files)

        except Exception as e:
            logger.exception("Skipping placement package %s: %s", package_dir.name, e)
            continue

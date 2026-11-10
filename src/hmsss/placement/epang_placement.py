#!/usr/bin/env python3
"""Phylogenetic placement of protein-domain candidates using EPA-ng and GAPPA.

Workflow:
    query proteins
        -> PyHMMER alignment against package HMM
        -> EPA-ng placement
        -> GAPPA functional/domain assignment
        -> GAPPA taxonomic assignment

Functional and taxonomic GAPPA workflows are independent.
Failures affect only the current operation/package and do not terminate
the complete pipeline.

Required config attributes:
    config.paths.placement
    config.placement_results_directory
    config.database_directory
    config.cores

    config.placements
        Optional package restriction. Empty list = all matching packages.

    config.placement_threshold
        GAPPA functional clade threshold.

    config.placement_lwr
        EPA-ng dynamic heuristic cumulative LWR threshold.
        Recommended default: 0.99999.

    config.taxonomy_lwr
        Minimum accumulated LWR (aLWR) for taxonomic assignment.
        Recommended default: 0.5.
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
from typing import Dict, Optional, Union

import pyhmmer

from hmsss.db import database

logger = logging.getLogger("hmg.epa")


# ======================================================================
# Placement package discovery
# ======================================================================


def find_matching_placement_packages(placement_root, intermediate_files):
    """Yield placement packages matching domains present in intermediate files."""

    placement_root = Path(placement_root)

    if not placement_root.is_dir():
        logger.error("Placement package directory does not exist: %s", placement_root)
        return

    try:
        package_dirs = sorted(placement_root.iterdir())
    except OSError as exc:
        logger.error(
            "Cannot read placement package directory %s: %s", placement_root, exc
        )
        return

    for package_dir in package_dirs:
        annotations = package_dir / "annotations.tsv"

        if not package_dir.is_dir() or not annotations.is_file():
            continue

        try:
            with annotations.open(encoding="utf-8") as handle:
                protein_types = {
                    row[1]
                    for row in csv.reader(handle, delimiter="\t")
                    if len(row) >= 2
                }
        except (OSError, UnicodeError, csv.Error) as exc:
            logger.warning("Skipping placement package %s: %s", package_dir.name, exc)
            continue

        matched_types = protein_types & intermediate_files.keys()

        if not matched_types:
            continue

        faa_files = []

        for protein_type in sorted(matched_types):
            files = intermediate_files[protein_type]
            faa_files.extend([files] if isinstance(files, (str, Path)) else files)

        yield package_dir, list(dict.fromkeys(map(Path, faa_files)))


# ======================================================================
# FASTA handling
# ======================================================================


def fasta_records(path: Path):
    """Stream FASTA records as (first-word identifier, sequence)."""

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
                continue

            if name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")

            parts.append(line)

    if name is not None:
        yield name, "".join(parts)


def merge_placement_fastas(faa_files, output_file):
    """Merge FASTA files while preserving identifiers and removing duplicates."""

    output_file = Path(output_file)
    seen_ids = set()

    with output_file.open("w", encoding="utf-8") as out:
        for faa_file in faa_files:
            faa_file = Path(faa_file)

            try:
                for seq_id, sequence in fasta_records(faa_file):
                    if seq_id in seen_ids:
                        continue

                    if not sequence:
                        logger.warning(
                            "Skipping empty sequence %s in %s", seq_id, faa_file
                        )
                        continue

                    seen_ids.add(seq_id)
                    out.write(f">{seq_id}\n{sequence}\n")

            except (OSError, ValueError, UnicodeError) as exc:
                logger.warning("Skipping faulty FASTA file %s: %s", faa_file, exc)

    return output_file


def _fasta_records_from_handle(handle):
    """Yield FASTA records from a binary or text handle."""

    name, parts = None, []

    for line in handle:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)

            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if name is not None:
        yield name, "".join(parts)


# ======================================================================
# Placement package loading
# ======================================================================


def load_placement_package(
    package_dir: Union[str, Path],
) -> Dict[str, Path]:
    """Load and validate EPA-ng reference files using SHA-256 checksums."""

    package_dir = Path(package_dir).resolve()
    manifest_file = package_dir / "manifest.json"

    with manifest_file.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

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

        expected_sha256 = manifest["files"][filename]["sha256"]

        digest = hashlib.sha256()

        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)

        actual_sha256 = digest.hexdigest()

        if actual_sha256.lower() != expected_sha256.lower():
            raise ValueError(
                f"SHA-256 checksum mismatch for {path}\n"
                f"Expected: {expected_sha256}\n"
                f"Actual:   {actual_sha256}"
            )

        package_files[argument] = path

    logger.debug(
        "Reference package '%s': SHA-256 checksums verified",
        package_dir.name,
    )

    return package_files


# ======================================================================
# PyHMMER alignment
# ======================================================================


def _sequence_name(sequence) -> str:
    """Return a PyHMMER sequence identifier as str."""

    name = sequence.name
    return name.decode("utf-8") if isinstance(name, bytes) else str(name)


def _align_block(model, block, output_handle):
    """Align one PyHMMER sequence block and retain HMM match-state columns."""

    block_ids = {_sequence_name(sequence) for sequence in block}

    try:
        msa = pyhmmer.hmmer.hmmalign(
            model,
            block,
            trim=False,
            all_consensus_cols=True,
        )
    except Exception as exc:
        logger.warning(
            "hmmalign failed for batch of %d queries: %s",
            len(block_ids),
            exc,
        )
        return set(), block_ids

    buffer = io.BytesIO()
    msa.write(buffer, "a2m")
    buffer.seek(0)

    valid_ids = set()
    failed_ids = set()

    for query_id, sequence in _fasta_records_from_handle(buffer):
        # A2M:
        # uppercase = HMM match state
        # lowercase = insertion state
        # "-"       = gap
        aligned = "".join(char for char in sequence if char.isupper() or char == "-")

        if len(aligned) != model.M:
            logger.warning(
                "Invalid masked alignment length for %s: %d != %d",
                query_id,
                len(aligned),
                model.M,
            )
            failed_ids.add(query_id)
            continue

        if not aligned.strip("-"):
            failed_ids.add(query_id)
            continue

        output_handle.write(f">{query_id}\n{aligned}\n")
        valid_ids.add(query_id)

    failed_ids.update(block_ids - valid_ids - failed_ids)

    return valid_ids, failed_ids


def align_queries_pyhmmer(
    hmm: Path,
    queries: Path,
    output: Path,
    batch_size: int = 2056,
) -> int:
    """
    Align queries to the package HMM.

    Failed batch members are retried individually. Queries without
    residues in HMM match-state columns are excluded from EPA-ng.
    """

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.plan7.HMMFile(str(hmm)) as hmm_file:
        model = hmm_file.read()

    if model is None:
        raise ValueError(f"No HMM found: {hmm}")

    output.parent.mkdir(parents=True, exist_ok=True)

    aligned_ids = set()
    failed_ids = set()

    # First pass: process queries in batches.
    with (
        pyhmmer.easel.SequenceFile(
            str(queries),
            digital=True,
            alphabet=alphabet,
        ) as seq_file,
        output.open("w", encoding="utf-8") as out,
    ):
        while True:
            block = seq_file.read_block(sequences=batch_size)

            if not block:
                break

            valid, failed = _align_block(model, block, out)
            aligned_ids.update(valid)
            failed_ids.update(failed)

    failed_ids -= aligned_ids
    rescued_ids = set()

    # Retry failed batch members individually.
    if failed_ids:
        rescue_input = output.with_name("alignment_rescue_queries.faa")

        with rescue_input.open("w", encoding="utf-8") as handle:
            for query_id, sequence in fasta_records(queries):
                if query_id in failed_ids and sequence:
                    handle.write(f">{query_id}\n{sequence}\n")

        if rescue_input.stat().st_size:
            with (
                pyhmmer.easel.SequenceFile(
                    str(rescue_input),
                    digital=True,
                    alphabet=alphabet,
                ) as seq_file,
                output.open("a", encoding="utf-8") as out,
            ):
                while True:
                    block = seq_file.read_block(sequences=1)

                    if not block:
                        break

                    valid, _ = _align_block(model, block, out)
                    rescued_ids.update(valid)

        rescue_input.unlink(missing_ok=True)

    aligned_ids.update(rescued_ids)
    unresolved_ids = failed_ids - rescued_ids

    unplaceable_file = output.with_name("unplaceable_queries.faa")

    if unresolved_ids:
        with unplaceable_file.open("w", encoding="utf-8") as handle:
            for query_id, sequence in fasta_records(queries):
                if query_id in unresolved_ids:
                    handle.write(f">{query_id}\n{sequence}\n")

        logger.warning(
            "%d queries could not be aligned to HMM match-state columns "
            "and were written to %s",
            len(unresolved_ids),
            unplaceable_file,
        )
    else:
        unplaceable_file.unlink(missing_ok=True)

    logger.info(
        "PyHMMER alignment: %d usable, %d rescued, %d unplaceable",
        len(aligned_ids),
        len(rescued_ids),
        len(unresolved_ids),
    )

    return len(aligned_ids)


# ======================================================================
# EPA-ng
# ======================================================================


def run_epa_ng(
    tree: Path,
    seed: Path,
    hmm: Path,
    queries: Path,
    model: str,
    output_dir: Path,
    threads: int,
    placement_lwr: float = 0.99999,
) -> Optional[Path]:
    """
    Align queries and run EPA-ng.

    placement_lwr controls EPA-ng's dynamic preplacement heuristic.
    It is the cumulative LWR that candidate branches must cover.
    """

    if not 0.0 < placement_lwr <= 1.0:
        logger.error(
            "Invalid EPA-ng placement LWR %.6f; expected 0 < LWR <= 1",
            placement_lwr,
        )
        return None

    if not shutil.which("epa-ng"):
        logger.warning("EPA-ng executable not found")
        return None

    output_dir.mkdir(parents=True, exist_ok=True)

    aligned_queries = output_dir / "aligned_queries.fasta"
    jplace = output_dir / "epa_result.jplace"

    jplace.unlink(missing_ok=True)

    try:
        aligned_count = align_queries_pyhmmer(
            hmm=hmm,
            queries=queries,
            output=aligned_queries,
        )
    except Exception:
        logger.exception("Query alignment failed for EPA-ng package")
        return None

    if not aligned_count:
        logger.warning("No alignable queries available for EPA-ng")
        return None

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
        "--dyn-heur",
        str(placement_lwr),
        "--outdir",
        str(output_dir),
    ]

    logger.info(
        "Running EPA-ng with %d queries; dynamic LWR threshold %.5f",
        aligned_count,
        placement_lwr,
    )

    try:
        result = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        logger.error("Could not execute EPA-ng: %s", exc)
        return None

    if result.returncode != 0:
        logger.error(
            "EPA-ng failed with exit code %d\n%s",
            result.returncode,
            result.stderr.strip(),
        )
        return None

    if not jplace.is_file() or jplace.stat().st_size == 0:
        logger.error("EPA-ng produced no valid JPLACE file: %s", jplace)
        return None

    try:
        with jplace.open(encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        logger.error("Could not parse EPA-ng output %s: %s", jplace, exc)
        return None

    if not {"tree", "fields", "placements"} <= data.keys():
        logger.error("Invalid JPLACE structure: %s", jplace)
        return None

    if not data["placements"]:
        logger.warning("EPA-ng returned an empty placement set: %s", jplace)
        return None

    logger.info(
        "EPA-ng completed: %d placement groups",
        len(data["placements"]),
    )

    return jplace


# ======================================================================
# GAPPA functional/domain assignment
# ======================================================================


def assign_protein_clades(
    jplace: Path,
    annotations: Path,
    output_dir: Path,
    threshold: float = 0.75,
) -> Optional[Path]:
    """
    Assign placements to functional protein clades.

    If GAPPA fails, progressively less strict thresholds are tried.
    """

    if not shutil.which("gappa"):
        logger.warning("GAPPA executable not found; functional assignment skipped")
        return None

    if not jplace.is_file():
        logger.warning(
            "JPLACE file missing for functional assignment: %s",
            jplace,
        )
        return None

    if not annotations.is_file():
        logger.warning(
            "Functional annotation file missing: %s",
            annotations,
        )
        return None

    output_dir.mkdir(parents=True, exist_ok=True)

    thresholds = [
        threshold,
        threshold * 0.8,
        threshold * 0.6,
        threshold * 0.4,
        threshold * 0.2,
        0.0,
    ]
    thresholds = list(dict.fromkeys(max(0.0, value) for value in thresholds))

    for current_threshold in thresholds:
        for old_file in output_dir.glob("*.jplace"):
            old_file.unlink()

        cmd = [
            "gappa",
            "edit",
            "extract",
            "--jplace-path",
            str(jplace.resolve()),
            "--clade-list-file",
            str(annotations.resolve()),
            "--threshold",
            str(current_threshold),
        ]

        logger.info(
            "Running GAPPA functional assignment with threshold %.3f",
            current_threshold,
        )

        try:
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except OSError as exc:
            logger.error(
                "Could not execute GAPPA functional assignment: %s",
                exc,
            )
            return None

        if result.returncode != 0:
            logger.warning(
                "GAPPA functional assignment failed at threshold %.3f: %s",
                current_threshold,
                result.stderr.strip(),
            )
            continue

        if not list(output_dir.glob("*.jplace")):
            logger.warning(
                "GAPPA functional assignment at threshold %.3f "
                "produced no JPLACE files",
                current_threshold,
            )
            continue

        if current_threshold != threshold:
            logger.warning(
                "Functional placement required relaxed threshold %.3f -> %.3f",
                threshold,
                current_threshold,
            )

        return output_dir

    logger.error("GAPPA functional assignment failed at all tested thresholds")

    return None


def parse_function_assignments(
    function_dir: Union[str, Path],
    database_path: str,
) -> int:
    """Store functional GAPPA assignments in the existing database."""

    function_dir = Path(function_dir)

    if not function_dir.is_dir():
        return 0

    assignments = {}

    for file in sorted(function_dir.glob("*.jplace")):
        domain_type = file.stem

        try:
            with file.open(encoding="utf-8") as handle:
                data = json.load(handle)
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning(
                "Cannot parse functional GAPPA output %s: %s",
                file,
                exc,
            )
            continue

        for placement in data.get("placements", []):
            names = list(placement.get("n", []))
            names.extend(entry[0] for entry in placement.get("nm", []) if entry)

            for name in names:
                try:
                    domain_hit = database.parse_domain_header(name)
                except ValueError as exc:
                    logger.warning(
                        "Cannot map placement query '%s' to database hit: %s",
                        name,
                        exc,
                    )
                    continue

                previous = assignments.get(domain_hit)

                if previous is not None and previous != domain_type:
                    logger.warning(
                        "Conflicting placement assignments for %s: "
                        "%s versus %s; keeping first",
                        domain_hit,
                        previous,
                        domain_type,
                    )
                    continue

                assignments[domain_hit] = domain_type

    if not assignments:
        logger.info("No functional placement assignments found")
        return 0

    grouped = {}

    for domain_hit, domain_type in assignments.items():
        grouped.setdefault(domain_type, set()).add(domain_hit)

    updated = 0

    for domain_type, domain_hits in grouped.items():
        if domain_type.lower() in {"basal", "uncertain"}:
            continue

        updated += database.update_domains(
            database=database_path,
            domain_hits=domain_hits,
            new_domain=domain_type,
        )

        database.update_selection_comments(
            database=database_path,
            domain_hits=domain_hits,
            comment="pp",
        )

    logger.info(
        "Placement changed %d functional domain assignments",
        updated,
    )

    return updated


# ======================================================================
# GAPPA taxonomic assignment
# ======================================================================


def assign_taxonomy(
    jplace: Path,
    taxonomy_file: Path,
    output_dir: Path,
    threads: int = 1,
) -> Optional[Path]:
    """
    Calculate per-query taxonomic LWR profiles using GAPPA.

    --best-hit is deliberately NOT used. The complete per-query profile is
    retained so the deepest sufficiently supported taxonomic path can later
    be selected using accumulated LWR.
    """

    if not shutil.which("gappa"):
        logger.warning("GAPPA executable not found; taxonomic assignment skipped")
        return None

    if not jplace.is_file():
        logger.warning(
            "JPLACE file missing for taxonomic assignment: %s",
            jplace,
        )
        return None

    if not taxonomy_file.is_file():
        logger.warning(
            "Reference taxonomy file missing: %s",
            taxonomy_file,
        )
        return None

    output_dir.mkdir(parents=True, exist_ok=True)

    for old_file in output_dir.glob("*per_query.tsv"):
        old_file.unlink()

    cmd = [
        "gappa",
        "examine",
        "assign",
        "--jplace-path",
        str(jplace.resolve()),
        "--taxon-file",
        str(taxonomy_file.resolve()),
        "--per-query-results",
        "--threads",
        str(threads),
        "--out-dir",
        str(output_dir.resolve()),
    ]

    logger.info("Running GAPPA taxonomic assignment")

    try:
        result = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        logger.error(
            "Could not execute GAPPA taxonomic assignment: %s",
            exc,
        )
        return None

    if result.returncode != 0:
        logger.error(
            "GAPPA taxonomic assignment failed: %s",
            result.stderr.strip(),
        )
        return None

    if not list(output_dir.glob("*per_query.tsv")):
        logger.error("GAPPA taxonomic assignment produced no per-query result file")
        return None

    return output_dir


def parse_taxonomy_assignments(
    taxonomy_dir: Union[str, Path],
    minimum_lwr: float = 0.5,
) -> Dict[str, str]:
    """
    Select the deepest sufficiently supported taxonomic assignment.

    GAPPA per-query output contains:

        name
        LWR
        fract
        aLWR
        afract
        taxopath

    For each query, the deepest taxopath with:

        aLWR >= minimum_lwr

    is selected.

    Using accumulated LWR rather than a single-placement best hit allows
    uncertain placements to be retained at a higher, better-supported
    taxonomic rank instead of being discarded completely or assigned
    spuriously to a weakly supported species.
    """

    if not 0.0 <= minimum_lwr <= 1.0:
        raise ValueError(f"taxonomy_lwr must be between 0 and 1, got {minimum_lwr}")

    candidates = {}

    for file in Path(taxonomy_dir).glob("*per_query.tsv"):
        try:
            with file.open(
                newline="",
                encoding="utf-8",
            ) as handle:
                reader = csv.DictReader(
                    handle,
                    delimiter="\t",
                )

                required = {
                    "name",
                    "aLWR",
                    "taxopath",
                }

                if not reader.fieldnames or not required <= set(reader.fieldnames):
                    logger.warning(
                        "Invalid GAPPA taxonomy table: %s",
                        file,
                    )
                    continue

                for row in reader:
                    name = (row.get("name") or "").strip()
                    taxopath = (row.get("taxopath") or "").strip()
                    alwr_raw = (row.get("aLWR") or "").strip()

                    if not name or not taxopath or not alwr_raw:
                        continue

                    try:
                        alwr = float(alwr_raw)
                    except ValueError:
                        continue

                    if alwr < minimum_lwr:
                        continue

                    # Number of taxonomic levels in this path.
                    depth = sum(1 for part in taxopath.split(";") if part.strip())

                    previous = candidates.get(name)

                    # Prefer the deepest supported taxonomic path.
                    # If depth is equal, prefer the one with higher aLWR.
                    if (
                        previous is None
                        or depth > previous[0]
                        or (depth == previous[0] and alwr > previous[1])
                    ):
                        candidates[name] = (
                            depth,
                            alwr,
                            taxopath,
                        )

        except (OSError, UnicodeError, csv.Error) as exc:
            logger.warning(
                "Could not parse GAPPA taxonomy file %s: %s",
                file,
                exc,
            )

    assignments = {name: values[2] for name, values in candidates.items()}

    logger.info(
        "Taxonomic placement: %d queries assigned at minimum aLWR %.3f",
        len(assignments),
        minimum_lwr,
    )

    return assignments


# ======================================================================
# Complete placement workflow
# ======================================================================


def epa_ng_placement(config, intermediate_files) -> None:
    """
    Run phylogenetic placement for all matching packages.

    Errors are isolated:
        package failure -> skip package
        functional GAPPA failure -> still run taxonomy
        taxonomy failure -> retain functional result
        database failure -> skip only affected result
    """

    processed_packages = 0
    failed_packages = 0

    try:
        packages = find_matching_placement_packages(
            config.paths.placement,
            intermediate_files,
        )

        if packages is None:
            return

        for package_dir, faa_files in packages:
            try:
                # Empty list = all matching packages.
                if config.placements and package_dir.name not in config.placements:
                    continue

                output_dir = Path(config.placement_results_directory) / package_dir.name

                output_dir.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                # --------------------------------------------------
                # Load package
                # --------------------------------------------------

                try:
                    package = load_placement_package(package_dir)
                except Exception:
                    logger.exception(
                        "Skipping package %s: package loading failed",
                        package_dir.name,
                    )
                    failed_packages += 1
                    continue

                # --------------------------------------------------
                # Prepare queries
                # --------------------------------------------------

                try:
                    queries = merge_placement_fastas(
                        faa_files,
                        output_dir / "placement_candidates.faa",
                    )
                except Exception:
                    logger.exception(
                        "Skipping package %s: query preparation failed",
                        package_dir.name,
                    )
                    failed_packages += 1
                    continue

                if not queries.is_file() or queries.stat().st_size == 0:
                    logger.warning(
                        "Skipping package %s: no placement candidates",
                        package_dir.name,
                    )
                    continue

                logger.info(
                    "Running placement package %s",
                    package_dir.name,
                )

                # --------------------------------------------------
                # EPA-ng
                # --------------------------------------------------

                try:
                    jplace = run_epa_ng(
                        tree=package["tree"],
                        seed=package["seed"],
                        hmm=package["hmm"],
                        queries=queries,
                        model="LG",
                        output_dir=output_dir,
                        threads=config.cores,
                        placement_lwr=config.placement_lwr,
                    )
                except Exception:
                    logger.exception(
                        "Skipping package %s: EPA-ng failed unexpectedly",
                        package_dir.name,
                    )
                    failed_packages += 1
                    continue

                if jplace is None:
                    logger.warning(
                        "Skipping package %s: no usable EPA-ng placement",
                        package_dir.name,
                    )
                    failed_packages += 1
                    continue

                # ==================================================
                # GAPPA workflow 1:
                # functional/domain assignment
                # ==================================================

                try:
                    function_dir = assign_protein_clades(
                        jplace=jplace,
                        annotations=package_dir / "annotations.tsv",
                        output_dir=output_dir / "function",
                        threshold=config.placement_threshold,
                    )

                    if function_dir is not None:
                        try:
                            parse_function_assignments(
                                function_dir=function_dir,
                                database_path=config.database_directory,
                            )
                        except Exception:
                            logger.exception(
                                "Functional assignments for package %s "
                                "could not be stored; skipping them",
                                package_dir.name,
                            )

                except Exception:
                    logger.exception(
                        "Functional GAPPA failed for package %s; "
                        "continuing with taxonomy",
                        package_dir.name,
                    )

                # ==================================================
                # GAPPA workflow 2:
                # taxonomic assignment
                # ==================================================

                try:
                    taxonomy_dir = assign_taxonomy(
                        jplace=jplace,
                        taxonomy_file=(package_dir / "reference_taxonomy.tsv"),
                        output_dir=output_dir / "taxonomy",
                        threads=config.cores,
                    )

                    if taxonomy_dir is not None:
                        domain_lineages = parse_taxonomy_assignments(
                            taxonomy_dir,
                            minimum_lwr=config.taxonomy_lwr,
                        )

                        if domain_lineages:
                            database.insert_domain_placements(
                                database=config.database_directory,
                                domain_lineages=domain_lineages,
                            )

                except Exception:
                    logger.exception(
                        "Taxonomic GAPPA failed for package %s; "
                        "continuing with next package",
                        package_dir.name,
                    )

                processed_packages += 1

            except Exception:
                logger.exception(
                    "Unexpected error while processing package %s; skipping package",
                    package_dir.name,
                )
                failed_packages += 1
                continue

    except Exception:
        logger.exception(
            "Unexpected error in placement stage; "
            "remaining placement processing skipped"
        )

    logger.info(
        "Placement stage finished: %d packages processed, %d packages failed",
        processed_packages,
        failed_packages,
    )

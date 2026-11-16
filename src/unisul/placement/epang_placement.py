#!/usr/bin/env python3
"""EPA-ng placement followed by GAPPA functional and rank-aware taxonomy assignment.

Workflow
--------
1. Match intermediate FASTA files to selected placement packages.
2. Run EPA-ng for ALL matching packages first.
3. Persist one package-specific JPLACE file in:
       config.placement_results_directory/<package>/<package>.epa_result.jplace
4. Only after placement has finished for all packages:
       a. run GAPPA functional/clade assignment and update domain types;
       b. run GAPPA taxonomy assignment without --best-hit;
       c. apply rank-aware accumulated-support classification using ``afract``;
       d. write taxonomy to DomainPlacements with method="epa_ng".

The implementation is intentionally disk-oriented for large metagenomic datasets:
- intermediate FASTA files are streamed;
- aligned queries are processed in bounded PyHMMER batches;
- package HMM and reference MSA share the same alignment coordinate system;
- every query uses the canonical header ``proteinID___dom_START_END``;
- no global query-ID or taxonomy dictionaries are retained;
- GAPPA JPLACE and TSV outputs are consumed directly in bounded batches;
- no temporary SQLite databases are created.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import logging
import re
import subprocess
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

import pyhmmer

from unisul.db import database

logger = logging.getLogger("hmg.epa")

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
DB_BATCH_SIZE = 5000
ALIGN_BATCH_SIZE = 2048


# ============================================================================
# Basic FASTA / package handling
# ============================================================================


def fasta_records(path: Union[str, Path]) -> Iterator[Tuple[str, str]]:
    """Stream FASTA records as ``(identifier, sequence)``."""
    name = None
    parts: List[str] = []

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
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


def find_matching_placement_packages(
    placement_root: Union[str, Path],
    intermediate_files,
):
    """Yield placement packages that share protein types with intermediate files."""
    root = Path(placement_root)
    if not root.is_dir():
        logger.warning("Placement package directory not found: %s", root)
        return

    intermediate_types = set(intermediate_files)

    try:
        package_dirs = sorted(root.iterdir())
    except OSError as exc:
        logger.error("Cannot read placement package directory %s: %s", root, exc)
        return

    for package_dir in package_dirs:
        if not package_dir.is_dir():
            continue

        annotations = package_dir / "annotations.tsv"
        if not annotations.is_file():
            continue

        try:
            protein_types = set()
            with annotations.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.reader(handle, delimiter="\t")
                for row in reader:
                    if len(row) >= 2 and row[1] and row[1] != "protein_type":
                        protein_types.add(row[1].strip())
        except Exception as exc:
            logger.exception(
                "Skipping placement package %s: cannot read annotations.tsv: %s",
                package_dir.name,
                exc,
            )
            continue

        matched = sorted(protein_types & intermediate_types)
        if not matched:
            continue

        faa_files: List[Path] = []
        try:
            for protein_type in matched:
                value = intermediate_files[protein_type]
                if isinstance(value, (str, Path)):
                    faa_files.append(Path(value))
                else:
                    faa_files.extend(Path(path) for path in value)
        except Exception as exc:
            logger.exception(
                "Skipping placement package %s: invalid intermediate file mapping: %s",
                package_dir.name,
                exc,
            )
            continue

        yield package_dir, faa_files


def load_placement_package(package_dir: Union[str, Path]) -> Dict[str, Path]:
    """Resolve the three files required by EPA-ng from the package manifest."""
    package_dir = Path(package_dir)
    manifest = json.loads((package_dir / "manifest.json").read_text(encoding="utf-8"))
    roles = manifest["roles"]

    files = {
        "tree": package_dir / roles["tree"],
        "seed": package_dir / roles["trimmed_alignment"],
        "hmm": package_dir / roles["hmm"],
    }

    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing placement-package file(s): " + ", ".join(missing)
        )

    return files


def merge_placement_fastas(faa_files: Iterable[Path], output_file: Path) -> int:
    """Stream-concatenate intermediate FASTAs and return the query count.

    Intermediate read detection assigns each ORF to one best HMM type, so query
    identifiers are expected to be unique across these files. No global ID set
    is retained in RAM.
    """
    count = 0
    with output_file.open("w", encoding="utf-8") as out:
        for faa_file in faa_files:
            for seq_id, sequence in fasta_records(faa_file):
                if not sequence:
                    continue
                out.write(f">{seq_id}\n{sequence}\n")
                count += 1
    return count


def _reference_alignment_width(path: Path) -> int:
    """Return the column count of the reference alignment."""
    width = None

    for _, sequence in fasta_records(path):
        if width is None:
            width = len(sequence)
        elif len(sequence) != width:
            raise ValueError(f"Inconsistent reference alignment width in {path}")

    if width is None:
        raise ValueError(f"Empty reference alignment: {path}")

    return width


# ============================================================================
# Query alignment
# ============================================================================


def _sequence_name(sequence) -> str:
    name = sequence.name
    return name.decode("utf-8") if isinstance(name, bytes) else str(name)


def _fasta_records_from_handle(handle):
    name = None
    parts = []
    for raw in handle:
        line = raw.decode("utf-8").strip() if isinstance(raw, bytes) else raw.strip()
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


def _align_block(
    model,
    block,
    output_handle,
) -> Tuple[int, set]:
    """Align one bounded query block and retain HMM match-state columns."""
    block_ids = {_sequence_name(sequence) for sequence in block}

    try:
        msa = pyhmmer.hmmer.hmmalign(
            model,
            block,
            trim=False,
            all_consensus_cols=True,
        )
    except Exception:
        return 0, block_ids

    buffer = io.BytesIO()
    msa.write(buffer, "a2m")
    buffer.seek(0)

    written = 0
    returned = set()
    failed = set()

    for query_id, sequence in _fasta_records_from_handle(buffer):
        returned.add(query_id)

        aligned = "".join(char for char in sequence if char.isupper() or char == "-")

        if len(aligned) != model.M or not aligned.strip("-"):
            failed.add(query_id)
            continue

        output_handle.write(f">{query_id}\n{aligned}\n")
        written += 1

    failed.update(block_ids - returned)
    return written, failed


def align_queries_pyhmmer(
    hmm: Path,
    seed: Path,
    queries: Path,
    output: Path,
    unplaceable: Path,
    batch_size: int = ALIGN_BATCH_SIZE,
) -> Tuple[int, int]:
    """Align queries with PyHMMER in the package reference coordinate system.

    HMSS3 placement packages guarantee that the HMM has one match state per
    column of ``trimmed.fasta_aln``. The HMM can still be trained only from
    ingroup sequences; the coordinate system remains identical to the full
    Ingroup+Outgroup reference alignment used by EPA-ng.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.plan7.HMMFile(str(hmm)) as hmm_file:
        model = hmm_file.read()

    if model is None:
        raise ValueError(f"No HMM found in {hmm}")

    reference_width = _reference_alignment_width(seed)

    if model.M != reference_width:
        raise ValueError(
            "Placement package coordinate mismatch: "
            f"HMM match states={model.M}, "
            f"reference alignment width={reference_width}"
        )

    aligned_count = 0
    unplaceable_count = 0

    with (
        pyhmmer.easel.SequenceFile(
            str(queries),
            digital=True,
            alphabet=alphabet,
        ) as seq_file,
        output.open("w", encoding="utf-8") as aligned_out,
        unplaceable.open("w", encoding="utf-8") as failed_out,
    ):
        while True:
            block = seq_file.read_block(sequences=batch_size)
            if not block:
                break

            written, failed_ids = _align_block(
                model,
                block,
                aligned_out,
            )
            aligned_count += written

            if not failed_ids:
                continue

            for sequence in block:
                query_id = _sequence_name(sequence)
                if query_id not in failed_ids:
                    continue

                rescue = pyhmmer.easel.DigitalSequenceBlock(alphabet)
                rescue.append(sequence)

                rescued, still_failed = _align_block(
                    model,
                    rescue,
                    aligned_out,
                )

                if rescued:
                    aligned_count += rescued
                elif query_id in still_failed:
                    failed_out.write(f">{query_id}\n{sequence.textize().sequence}\n")
                    unplaceable_count += 1

    if unplaceable_count == 0:
        unplaceable.unlink(missing_ok=True)

    return aligned_count, unplaceable_count


# ============================================================================
# EPA-ng
# ============================================================================


def run_epa_ng(
    package_name: str,
    tree: Path,
    seed: Path,
    hmm: Path,
    queries: Path,
    output_dir: Path,
    threads: int,
    placement_lwr: float,
) -> Optional[Path]:
    """Run one EPA-ng package and persist a package-specific JPLACE file."""
    aligned_queries = output_dir / f"{package_name}.aligned_queries.fasta"
    unplaceable = output_dir / f"{package_name}.unplaceable_queries.faa"
    default_jplace = output_dir / "epa_result.jplace"
    final_jplace = output_dir / f"{package_name}.epa_result.jplace"

    default_jplace.unlink(missing_ok=True)

    aligned_count, failed_count = align_queries_pyhmmer(
        hmm=hmm,
        seed=seed,
        queries=queries,
        output=aligned_queries,
        unplaceable=unplaceable,
    )

    if aligned_count == 0:
        logger.warning("%s: no alignable queries", package_name)
        aligned_queries.unlink(missing_ok=True)
        return None

    logger.info(
        "%s: EPA-ng placement of %d queries (%d unplaceable)",
        package_name,
        aligned_count,
        failed_count,
    )

    cmd = [
        "epa-ng",
        "--tree",
        str(tree),
        "--ref-msa",
        str(seed),
        "--query",
        str(aligned_queries),
        "--model",
        "LG",
        "--threads",
        str(threads),
        "--dyn-heur",
        str(placement_lwr),
        "--outdir",
        str(output_dir),
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    aligned_queries.unlink(missing_ok=True)

    if result.returncode != 0:
        logger.error("%s: EPA-ng failed: %s", package_name, result.stderr.strip())
        default_jplace.unlink(missing_ok=True)
        return None

    if not default_jplace.is_file() or default_jplace.stat().st_size == 0:
        logger.error("%s: EPA-ng produced no JPLACE output", package_name)
        return None

    # One explicit filename per package. Package result directories remain
    # unchanged, but JPLACE files cannot collide with generic EPA-ng names.
    default_jplace.replace(final_jplace)
    return final_jplace


# ============================================================================
# Canonical domain headers
# ============================================================================


def _parse_domain_hits(query_ids: Sequence[str]) -> set:
    """Parse canonical ``proteinID___dom_START_END`` query identifiers."""
    hits = set()

    for query_id in query_ids:
        try:
            hits.add(database.parse_domain_header(query_id))
        except ValueError:
            logger.warning("Ignoring invalid placement query header: %s", query_id)

    return hits


def _domain_lineages(
    assignments: Sequence[Tuple[str, str]],
) -> Dict[str, str]:
    """Return canonical GAPPA query headers mapped to assigned taxonomy."""
    result = {}

    for query_id, taxopath in assignments:
        if not taxopath:
            continue

        try:
            database.parse_domain_header(query_id)
        except ValueError:
            logger.warning("Ignoring invalid taxonomy query header: %s", query_id)
            continue

        result[query_id] = taxopath

    return result


# ============================================================================
# GAPPA functional assignment
# ============================================================================


def _annotation_types_and_clade_file(
    annotations: Path,
    output_file: Path,
) -> set:
    """Write GAPPA's headerless clade list and return the small type set."""
    protein_types = set()

    with (
        annotations.open(newline="", encoding="utf-8-sig") as source,
        output_file.open("w", encoding="utf-8") as target,
    ):
        reader = csv.reader(source, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            sequence_id = row[0].strip()
            protein_type = row[1].strip()
            if not sequence_id or not protein_type:
                continue
            if sequence_id == "sequence_id" and protein_type == "protein_type":
                continue
            target.write(f"{sequence_id}\t{protein_type}\n")
            protein_types.add(protein_type)

    return protein_types


def assign_protein_clades(
    jplace: Path,
    annotations: Path,
    output_dir: Path,
    threshold: float,
) -> Tuple[Optional[Path], set]:
    """Run one fixed-threshold GAPPA functional assignment."""
    output_dir.mkdir(parents=True, exist_ok=True)

    for old in output_dir.glob("*.jplace*"):
        old.unlink()

    clade_file = output_dir / "gappa_clades.tsv"
    expected_types = _annotation_types_and_clade_file(annotations, clade_file)

    result = subprocess.run(
        [
            "gappa",
            "edit",
            "extract",
            "--jplace-path",
            str(jplace.resolve()),
            "--clade-list-file",
            str(clade_file.resolve()),
            "--threshold",
            str(threshold),
        ],
        cwd=output_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        logger.error("GAPPA functional assignment failed: %s", result.stderr.strip())
        return None, expected_types

    return output_dir, expected_types


def _protein_type_from_filename(path: Path, expected_types: set) -> Optional[str]:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".jplace"):
        name = name[:-7]

    lowered = name.lower()
    if re.search(r"(^|[._-])(basal|uncertain)($|[._-])", lowered):
        return None

    if name in expected_types:
        return name

    matches = [
        protein_type
        for protein_type in expected_types
        if re.search(
            r"(^|[._-])" + re.escape(protein_type) + r"($|[._-])",
            name,
        )
    ]

    return matches[0] if len(matches) == 1 else None


def _iter_jplace_query_ids(path: Path) -> Iterator[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        data = json.load(handle)

    for placement in data.get("placements", []):
        for entry in placement.get("n", []):
            if isinstance(entry, str):
                yield entry
            elif isinstance(entry, (list, tuple)) and entry:
                yield entry[0]

        for entry in placement.get("nm", []):
            if isinstance(entry, (list, tuple)) and entry:
                yield entry[0]


def store_function_assignments(
    function_dir: Path,
    expected_types: set,
    database_path: str,
    batch_size: int = DB_BATCH_SIZE,
) -> int:
    """Read GAPPA clade JPLACE files directly and update matching domain hits."""
    paths = sorted(function_dir.glob("*.jplace"))
    paths.extend(sorted(function_dir.glob("*.jplace.gz")))

    updated = 0

    for path in paths:
        protein_type = _protein_type_from_filename(path, expected_types)
        if protein_type is None:
            continue

        query_batch: List[str] = []

        for query_id in _iter_jplace_query_ids(path):
            query_batch.append(query_id)

            if len(query_batch) >= batch_size:
                hits = _parse_domain_hits(query_batch)
                if hits:
                    updated += database.update_domains(
                        database=database_path,
                        domain_hits=hits,
                        new_domain=protein_type,
                    )
                    database.update_selection_comments(
                        database=database_path,
                        domain_hits=hits,
                        comment="pp",
                    )
                query_batch.clear()

        if query_batch:
            hits = _parse_domain_hits(query_batch)
            if hits:
                updated += database.update_domains(
                    database=database_path,
                    domain_hits=hits,
                    new_domain=protein_type,
                )
                database.update_selection_comments(
                    database=database_path,
                    domain_hits=hits,
                    comment="pp",
                )

    return updated


# ============================================================================
# GAPPA rank-aware taxonomy
# ============================================================================


def assign_taxonomy(
    jplace: Path,
    taxonomy_file: Path,
    output_dir: Path,
    threads: int,
    consensus_threshold: Optional[float] = None,
) -> Optional[Path]:
    """Generate complete per-query GAPPA taxonomy profiles (no --best-hit)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    for old in output_dir.glob("*per_query*.tsv"):
        old.unlink()

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

    if consensus_threshold is not None:
        cmd.extend(["--consensus-thresh", str(consensus_threshold)])

    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        logger.error("GAPPA taxonomy assignment failed: %s", result.stderr.strip())
        return None

    return output_dir


def _split_taxopath(taxopath: str) -> Tuple[str, ...]:
    return tuple(token.strip() for token in taxopath.split(";") if token.strip())


def rank_aware_classification(
    candidates: Sequence[Tuple[int, str, float]],
    support_threshold: float,
) -> str:
    """Exact rank-wise logic used by the validation benchmark.

    At every successive rank, choose the best-supported taxopath consistent
    with the already accepted parent path. Support is GAPPA ``afract``.
    Classification stops at the first rank below ``support_threshold``.
    """
    accepted: Tuple[str, ...] = ()

    for depth in range(1, len(RANKS) + 1):
        best_path = None
        best_support = -1.0

        for row_depth, taxopath, support in candidates:
            if row_depth != depth:
                continue
            tokens = _split_taxopath(taxopath)
            if accepted and tokens[: len(accepted)] != accepted:
                continue
            if support > best_support or (
                support == best_support and taxopath > (best_path or "")
            ):
                best_path = taxopath
                best_support = support

        if best_path is None or best_support < support_threshold:
            break

        accepted = _split_taxopath(best_path)

    return ";".join(accepted)


def _iter_rank_aware_taxonomy(
    path: Path,
    support_threshold: float,
) -> Iterator[Tuple[str, str]]:
    """Stream rank-aware assignments from one GAPPA per-query TSV.

    GAPPA writes the rows belonging to a query together. Only the current
    query's taxonomic candidates are retained in memory.
    """
    current_query: Optional[str] = None
    candidates: List[Tuple[int, str, float]] = []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            query_id = (row.get("name") or row.get("query") or "").strip()
            taxopath = (row.get("taxopath") or "").strip()
            support_raw = (row.get("afract") or "").strip()

            if not query_id or not taxopath or not support_raw:
                continue

            try:
                support = float(support_raw)
            except ValueError:
                continue

            depth = len(_split_taxopath(taxopath))
            if depth == 0:
                continue

            if current_query is not None and query_id != current_query:
                taxonomy = rank_aware_classification(
                    candidates,
                    support_threshold=support_threshold,
                )
                if taxonomy:
                    yield current_query, taxonomy
                candidates = []

            current_query = query_id
            candidates.append((depth, taxopath, support))

    if current_query is not None:
        taxonomy = rank_aware_classification(
            candidates,
            support_threshold=support_threshold,
        )
        if taxonomy:
            yield current_query, taxonomy


def store_rank_aware_taxonomy(
    taxonomy_dir: Path,
    database_path: str,
    support_threshold: float,
    batch_size: int = DB_BATCH_SIZE,
) -> int:
    """Read GAPPA TSV files directly and store rank-aware domain taxonomy."""
    taxonomy_files = sorted(taxonomy_dir.glob("*per_query*.tsv"))

    if not taxonomy_files:
        logger.warning("No GAPPA per-query taxonomy table found in %s", taxonomy_dir)
        return 0

    total_written = 0
    assignment_batch: List[Tuple[str, str]] = []

    for path in taxonomy_files:
        for query_id, taxonomy in _iter_rank_aware_taxonomy(
            path,
            support_threshold=support_threshold,
        ):
            assignment_batch.append((query_id, taxonomy))

            if len(assignment_batch) >= batch_size:
                domain_lineages = _domain_lineages(assignment_batch)

                if domain_lineages:
                    total_written += database.insert_domain_placements(
                        database=database_path,
                        domain_lineages=domain_lineages,
                        method="epa_ng",
                    )

                assignment_batch.clear()

    if assignment_batch:
        domain_lineages = _domain_lineages(assignment_batch)

        if domain_lineages:
            total_written += database.insert_domain_placements(
                database=database_path,
                domain_lineages=domain_lineages,
                method="epa_ng",
            )

    return total_written


# ============================================================================
# Complete two-phase workflow
# ============================================================================


def epa_ng_placement(config, intermediate_files) -> None:
    """Place all packages first, then run GAPPA and update the database."""
    placement_lwr = float(getattr(config, "placement_lwr", 0.99999))
    function_threshold = float(getattr(config, "placement_threshold", 0.75))
    taxonomy_lwr = float(getattr(config, "taxonomy_lwr", 0.80))
    consensus_threshold = getattr(config, "taxonomy_consensus_threshold", None)

    if not 0.0 < placement_lwr <= 1.0:
        raise ValueError("placement_lwr must satisfy 0 < value <= 1")
    if not 0.0 <= function_threshold <= 1.0:
        raise ValueError("placement_threshold must be between 0 and 1")
    if not 0.0 <= taxonomy_lwr <= 1.0:
        raise ValueError("taxonomy_lwr must be between 0 and 1")

    placement_results = []

    # ------------------------------------------------------------------
    # PHASE 1: all EPA-ng placements
    # ------------------------------------------------------------------
    for package_dir, faa_files in find_matching_placement_packages(
        config.paths.placement,
        intermediate_files,
    ):
        output_dir = Path(config.placement_results_directory) / package_dir.name
        output_dir.mkdir(parents=True, exist_ok=True)
        query_file = output_dir / f"{package_dir.name}.placement_candidates.faa"
        final_jplace = output_dir / f"{package_dir.name}.epa_result.jplace"

        # Resume support: a completed EPA-ng placement is reused instead of
        # being recalculated. It is still passed to the later GAPPA phase.
        if final_jplace.is_file() and final_jplace.stat().st_size > 0:
            logger.info(
                "%s: existing EPA-ng placement found; skipping EPA-ng",
                package_dir.name,
            )
            placement_results.append((package_dir, output_dir, final_jplace))
            continue

        try:
            query_count = merge_placement_fastas(faa_files, query_file)
            if query_count == 0:
                query_file.unlink(missing_ok=True)
                continue

            package = load_placement_package(package_dir)

            logger.info(
                "%s: starting EPA-ng placement for %d intermediate sequences",
                package_dir.name,
                query_count,
            )

            jplace = run_epa_ng(
                package_name=package_dir.name,
                tree=package["tree"],
                seed=package["seed"],
                hmm=package["hmm"],
                queries=query_file,
                output_dir=output_dir,
                threads=int(config.cores),
                placement_lwr=placement_lwr,
            )

            if jplace is not None:
                placement_results.append((package_dir, output_dir, jplace))

        except Exception:
            logger.exception(
                "%s: placement failed; continuing with remaining packages",
                package_dir.name,
            )
        finally:
            # Source intermediate files remain untouched; the merged copy is
            # unnecessary once EPA-ng has finished.
            query_file.unlink(missing_ok=True)

    logger.info(
        "EPA-ng phase completed: %d package JPLACE files ready for GAPPA",
        len(placement_results),
    )

    # ------------------------------------------------------------------
    # PHASE 2: GAPPA functional + rank-aware taxonomy
    # ------------------------------------------------------------------
    for package_dir, output_dir, jplace in placement_results:
        logger.info("%s: starting GAPPA assignments", package_dir.name)

        # Functional/domain assignment.
        try:
            function_dir, expected_types = assign_protein_clades(
                jplace=jplace,
                annotations=package_dir / "annotations.tsv",
                output_dir=output_dir / "function",
                threshold=function_threshold,
            )

            if function_dir is not None:
                updated = store_function_assignments(
                    function_dir=function_dir,
                    expected_types=expected_types,
                    database_path=config.database_directory,
                )
                logger.info(
                    "%s: GAPPA updated %d functional domain assignments",
                    package_dir.name,
                    updated,
                )
        except Exception:
            logger.exception(
                "%s: functional GAPPA failed; taxonomy will still be attempted",
                package_dir.name,
            )

        # Rank-aware taxonomy assignment.
        try:
            taxonomy_dir = assign_taxonomy(
                jplace=jplace,
                taxonomy_file=package_dir / "reference_taxonomy.tsv",
                output_dir=output_dir / "taxonomy",
                threads=int(config.cores),
                consensus_threshold=consensus_threshold,
            )

            if taxonomy_dir is not None:
                stored = store_rank_aware_taxonomy(
                    taxonomy_dir=taxonomy_dir,
                    database_path=config.database_directory,
                    support_threshold=taxonomy_lwr,
                )
                logger.info(
                    "%s: stored %d rank-aware EPA-ng taxonomy assignments "
                    "(minimum afract %.3f)",
                    package_dir.name,
                    stored,
                    taxonomy_lwr,
                )
        except Exception:
            logger.exception(
                "%s: taxonomic GAPPA failed; continuing with remaining packages",
                package_dir.name,
            )

    logger.info(
        "EPA-ng/GAPPA stage finished: %d packages placed and evaluated",
        len(placement_results),
    )

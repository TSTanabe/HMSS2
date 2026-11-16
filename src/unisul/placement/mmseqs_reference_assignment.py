#!/usr/bin/env python3
"""MMseqs2-based competitive reference and taxonomy assignment for UNISUL.

Design goals
------------
- Reuse placement-package files directly.
- Never copy ingroup.faa or outgroup.faa.
- Never create combined/reference FASTA files.
- Never create intermediate.faa or accepted_ingroup.faa.
- Preserve original reference identifiers throughout.
- Build one persistent MMseqs DB for ingroup and, if present, one for outgroup.
- Add taxonomy directly to the ingroup MMseqs DB and reuse that DB for both
  homology search and taxonomic assignment.
- Keep query/search MMseqs databases in a temporary directory and remove them
  automatically after each package.
- Keep only one combined assignments.tsv plus taxonomy/summary outputs in the result directory.\n- Stream MMseqs hit TSVs so all queries are evaluated without retaining all raw hits in RAM.\n- Remove per-package MMseqs search/taxonomy temporary databases immediately after use.
- Reuse persistent reference caches after interrupted/failed runs and remove the cache root only after a fully successful run.

Expected placement package
--------------------------
<package>/
    ingroup.faa
    outgroup.faa              # optional
    annotations.tsv           # reference sequence ID -> protein type
    reference_taxonomy.tsv    # reference sequence/genome ID -> GTDB lineage

Intermediate hit files
----------------------
<domain_type>{INTERMEDIATE_HITS_SUFFIX}

Intermediate sequence IDs are expected to follow:
    <proteinID>___dom_<start>_<end>

Requires
--------
- mmseqs on PATH
- unisul.db.database.insert_domain_placements(..., method="mmseqs")
- DomainPlacements primary key must include method if EPA-ng and MMseqs
  assignments shall coexist for the same domain hit.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

from unisul.cli.config import INTERMEDIATE_HITS_SUFFIX
from unisul.core.logging import get_logger
from unisul.cross_check.search_cross_reference import (
    set_valid_hit_true_for_protein_ids,
)
from unisul.db import database
from unisul.placement.epang_placement import find_matching_placement_packages
from unisul.utils import myUtil

logger = get_logger(__name__)


@dataclass
class HitCandidate:
    target: str
    protein_type: Optional[str]
    pident: float
    evalue: float
    bitscore: float
    qcov: float
    tcov: float


@dataclass
class QueryHitSummary:
    """Decision-relevant MMseqs hits for one query.

    Every query is still evaluated. We simply avoid retaining all raw hits.
    """

    best_ingroup: Optional[HitCandidate] = None
    best_outgroup: Optional[HitCandidate] = None
    best_by_type: Optional[Dict[str, HitCandidate]] = None

    def __post_init__(self) -> None:
        if self.best_by_type is None:
            self.best_by_type = {}


@dataclass
class QueryAssignment:
    query: str
    original_hmm_type: str
    status: str
    ingroup_assignment: bool
    valid_hit: bool
    final_type: str
    function_status: str
    best_ingroup_reference: Optional[str] = None
    best_ingroup_type: Optional[str] = None
    best_ingroup_bitscore: Optional[float] = None
    best_ingroup_identity: Optional[float] = None
    best_ingroup_qcov: Optional[float] = None
    best_ingroup_tcov: Optional[float] = None
    best_outgroup_reference: Optional[str] = None
    best_outgroup_bitscore: Optional[float] = None
    ingroup_outgroup_margin: Optional[float] = None
    second_type: Optional[str] = None
    second_type_bitscore: Optional[float] = None
    functional_margin: Optional[float] = None
    taxonomy: Optional[str] = None
    taxid: Optional[int] = None
    tax_rank: Optional[str] = None
    tax_name: Optional[str] = None


_GTDB_RANKS = {
    "d": "superkingdom",
    "k": "superkingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def _run(cmd: List[str], *, cwd: Optional[Path] = None) -> None:
    """Run an external command without writing directly to the terminal.

    stdout/stderr are captured so normal MMseqs progress output stays hidden.
    Successful command output is only emitted at DEBUG level; failed command
    output is logged at ERROR level before the exception is propagated.
    """
    command = [str(x) for x in cmd]
    logger.debug("Running: %s", " ".join(command))

    result = subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    if result.stdout:
        if result.returncode == 0:
            logger.debug("%s", result.stdout.rstrip())
        else:
            logger.error(
                "Command failed (%d): %s\n%s",
                result.returncode,
                " ".join(command),
                result.stdout.rstrip(),
            )

    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            command,
            output=result.stdout,
        )


def _safe_fraction(value: object) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError("Expected fraction in [0, 1], got %s" % value)
    return value


def _strip_header(name: str) -> str:
    return name.removeprefix(">").split()[0].strip()


def genome_id_from_sequence_id(sequence_id: str) -> str:
    sequence_id = _strip_header(sequence_id)
    return sequence_id.split("-", 1)[0].split("___", 1)[0]


def load_annotations(path: Union[str, Path]) -> Dict[str, str]:
    """Load reference sequence ID -> protein type."""
    mapping: Dict[str, str] = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2:
                continue
            seq_id, protein_type = row[0].strip(), row[1].strip()
            if not seq_id or not protein_type:
                continue
            if seq_id.lower() in {"sequence_id", "sequence", "protein_id", "proteinid"}:
                continue
            previous = mapping.get(seq_id)
            if previous is not None and previous != protein_type:
                logger.warning(
                    "Conflicting annotation for %s: %s / %s; keeping first",
                    seq_id,
                    previous,
                    protein_type,
                )
                continue
            mapping[seq_id] = protein_type
    return mapping


def load_taxonomy(path: Union[str, Path]) -> Dict[str, str]:
    """Load reference sequence/genome ID -> GTDB lineage from first two columns."""
    mapping: Dict[str, str] = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2:
                continue
            key, lineage = row[0].strip(), row[1].strip()
            if not key or not lineage:
                continue
            if key.lower() in {
                "sequence_id",
                "sequence",
                "protein_id",
                "proteinid",
                "genome_id",
                "genomeid",
            }:
                continue
            mapping.setdefault(key, lineage)
    return mapping


def resolve_taxonomy(sequence_id: str, taxonomy: Dict[str, str]) -> Optional[str]:
    """Resolve taxonomy by exact reference ID, then by genome ID."""
    lineage = taxonomy.get(sequence_id)
    if lineage:
        return lineage
    return taxonomy.get(genome_id_from_sequence_id(sequence_id))


def fasta_ids(path: Union[str, Path]) -> Iterable[str]:
    """Yield FASTA first-word identifiers without loading sequence strings."""
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                seq_id = raw[1:].split()[0].strip()
                if seq_id:
                    yield seq_id


def load_mmseqs_lookup(db_prefix: Union[str, Path]) -> List[str]:
    """Return original FASTA identifiers from an MMseqs sequence DB .lookup."""
    lookup = Path(str(db_prefix) + ".lookup")
    if not lookup.is_file():
        raise FileNotFoundError("MMseqs lookup missing: %s" % lookup)

    ids: List[str] = []
    with lookup.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1].strip():
                ids.append(fields[1].strip())
    return ids


def _db_marker(db_prefix: Union[str, Path]) -> Path:
    return Path(str(db_prefix) + ".dbtype")


def _remove_db_prefix(db_prefix: Union[str, Path]) -> None:
    """Remove files belonging to one MMseqs DB prefix."""
    db_prefix = Path(db_prefix)
    for path in db_prefix.parent.glob(db_prefix.name + "*"):
        if path.is_file() or path.is_symlink():
            path.unlink(missing_ok=True)
        elif path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


def _db_is_fresh(db_prefix: Union[str, Path], source_files: Iterable[Path]) -> bool:
    marker = _db_marker(db_prefix)
    lookup = Path(str(db_prefix) + ".lookup")
    if not marker.is_file() or not lookup.is_file():
        return False
    db_mtime = min(marker.stat().st_mtime, lookup.stat().st_mtime)
    return all(db_mtime >= path.stat().st_mtime for path in source_files)


def ensure_sequence_db(
    mmseqs: str,
    fasta_files: List[Path],
    db_prefix: Path,
) -> Path:
    """Create/reuse an MMseqs sequence DB directly from existing FASTA files."""
    if not fasta_files:
        raise ValueError("No FASTA files supplied for MMseqs database")
    for path in fasta_files:
        if not path.is_file():
            raise FileNotFoundError(path)

    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    if _db_is_fresh(db_prefix, fasta_files):
        return db_prefix

    _remove_db_prefix(db_prefix)
    _run([mmseqs, "createdb"] + [str(path) for path in fasta_files] + [str(db_prefix)])
    return db_prefix


def _parse_gtdb_lineage(lineage: str) -> List[Tuple[str, str, str]]:
    parsed: List[Tuple[str, str, str]] = []
    for token in lineage.split(";"):
        token = token.strip()
        if not token or "__" not in token:
            continue
        prefix, name = token.split("__", 1)
        prefix, name = prefix.strip(), name.strip()
        if prefix and name:
            parsed.append((prefix[0], name, "%s__%s" % (prefix[0], name)))
    return parsed


def build_local_taxdump(
    ingroup_ids: Iterable[str],
    taxonomy: Dict[str, str],
    taxdump_dir: Union[str, Path],
) -> Tuple[Path, Path, int, int]:
    """Build a small NCBI-style taxdump for the package-local GTDB lineages."""
    taxdump_dir = Path(taxdump_dir)
    taxdump_dir.mkdir(parents=True, exist_ok=True)

    nodes = taxdump_dir / "nodes.dmp"
    names = taxdump_dir / "names.dmp"
    merged = taxdump_dir / "merged.dmp"
    delnodes = taxdump_dir / "delnodes.dmp"
    mapping = taxdump_dir / "sequence_to_taxid.tsv"
    taxid_lineages = taxdump_dir / "taxid_to_lineage.tsv"

    path_to_taxid: Dict[Tuple[str, ...], int] = {(): 1}
    taxid_to_path: Dict[int, Tuple[str, ...]] = {1: ()}
    taxid_to_rank: Dict[int, str] = {1: "no rank"}
    taxid_to_name: Dict[int, str] = {1: "root"}
    sequence_taxid: Dict[str, int] = {}

    next_taxid = 2
    total = 0
    missing: List[str] = []

    for seq_id in ingroup_ids:
        total += 1
        lineage = resolve_taxonomy(seq_id, taxonomy)
        if not lineage:
            missing.append(seq_id)
            continue

        parsed = _parse_gtdb_lineage(lineage)
        if not parsed:
            missing.append(seq_id)
            continue

        current_path: Tuple[str, ...] = ()
        for prefix, name, token in parsed:
            child_path = current_path + (token,)
            if child_path not in path_to_taxid:
                taxid = next_taxid
                next_taxid += 1
                path_to_taxid[child_path] = taxid
                taxid_to_path[taxid] = child_path
                taxid_to_rank[taxid] = _GTDB_RANKS.get(prefix, "no rank")
                taxid_to_name[taxid] = name
            current_path = child_path
        sequence_taxid[seq_id] = path_to_taxid[current_path]

    with (
        nodes.open("w", encoding="utf-8") as node_out,
        names.open("w", encoding="utf-8") as name_out,
    ):
        node_out.write("1\t|\t1\t|\tno rank\t|\t-\t|\n")
        name_out.write("1\t|\troot\t|\t-\t|\tscientific name\t|\n")

        for taxid in sorted(taxid_to_path):
            if taxid == 1:
                continue
            path = taxid_to_path[taxid]
            parent_taxid = path_to_taxid[path[:-1]]
            node_out.write(
                "%d\t|\t%d\t|\t%s\t|\t-\t|\n"
                % (taxid, parent_taxid, taxid_to_rank[taxid])
            )
            name_out.write(
                "%d\t|\t%s\t|\t-\t|\tscientific name\t|\n"
                % (taxid, taxid_to_name[taxid])
            )

    merged.touch()
    delnodes.touch()

    with mapping.open("w", encoding="utf-8") as handle:
        for seq_id in sorted(sequence_taxid):
            handle.write("%s\t%d\n" % (seq_id, sequence_taxid[seq_id]))

    with taxid_lineages.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxid", "lineage"])
        for taxid in sorted(taxid_to_path):
            writer.writerow([taxid, ";".join(taxid_to_path[taxid])])

    if missing:
        logger.debug(
            "%d/%d ingroup references have no usable taxonomy. Examples: %s",
            len(missing),
            total,
            ", ".join(missing[:10]),
        )
    logger.debug(
        "Taxonomy mapping: %d/%d ingroup references mapped to local taxids",
        len(sequence_taxid),
        total,
    )

    if not sequence_taxid:
        raise ValueError("No ingroup reference could be mapped to a taxonomy")

    return mapping, taxid_lineages, len(sequence_taxid), total


def load_taxid_lineages(path: Union[str, Path]) -> Dict[int, str]:
    result: Dict[int, str] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[int(row["taxid"])] = (row.get("lineage") or "").strip()
    return result


def _taxdb_ready(ingroup_db: Path) -> bool:
    mapping = Path(str(ingroup_db) + "_mapping")
    taxonomy_binary = Path(str(ingroup_db) + "_taxonomy")
    nodes = Path(str(ingroup_db) + "_nodes.dmp")
    return mapping.is_file() and (taxonomy_binary.is_file() or nodes.is_file())


def _remove_taxdb_files(ingroup_db: Path) -> None:
    for suffix in (
        "_mapping",
        "_taxonomy",
        "_nodes.dmp",
        "_names.dmp",
        "_merged.dmp",
        "_delnodes.dmp",
    ):
        Path(str(ingroup_db) + suffix).unlink(missing_ok=True)


def ensure_taxonomy_db(
    mmseqs: str,
    ingroup_db: Path,
    ingroup_ids: List[str],
    taxonomy: Dict[str, str],
    taxonomy_source: Path,
    cache_dir: Path,
    threads: int,
) -> Tuple[Dict[int, str], int, int]:
    """Attach/reuse package-local GTDB taxonomy on the existing ingroup DB."""
    taxdump_dir = cache_dir / "taxonomy"
    mapping = taxdump_dir / "sequence_to_taxid.tsv"
    taxid_lineages = taxdump_dir / "taxid_to_lineage.tsv"
    nodes = taxdump_dir / "nodes.dmp"
    names = taxdump_dir / "names.dmp"

    lookup = Path(str(ingroup_db) + ".lookup")
    source_mtime = max(taxonomy_source.stat().st_mtime, lookup.stat().st_mtime)
    taxdump_fresh = (
        mapping.is_file()
        and taxid_lineages.is_file()
        and nodes.is_file()
        and names.is_file()
        and min(
            mapping.stat().st_mtime,
            taxid_lineages.stat().st_mtime,
            nodes.stat().st_mtime,
            names.stat().st_mtime,
        )
        >= source_mtime
    )

    if taxdump_fresh:
        mapped = 0
        with mapping.open(encoding="utf-8") as handle:
            mapped = sum(1 for line in handle if line.strip())
        total = len(ingroup_ids)
    else:
        mapping, taxid_lineages, mapped, total = build_local_taxdump(
            ingroup_ids=ingroup_ids,
            taxonomy=taxonomy,
            taxdump_dir=taxdump_dir,
        )

    tax_source_mtime = max(
        mapping.stat().st_mtime,
        nodes.stat().st_mtime,
        names.stat().st_mtime,
    )
    taxdb_files = [
        Path(str(ingroup_db) + "_mapping"),
        Path(str(ingroup_db) + "_taxonomy"),
        Path(str(ingroup_db) + "_nodes.dmp"),
    ]
    existing_taxdb = [path for path in taxdb_files if path.is_file()]
    taxdb_fresh = (
        _taxdb_ready(ingroup_db)
        and existing_taxdb
        and min(path.stat().st_mtime for path in existing_taxdb) >= tax_source_mtime
    )

    if not taxdb_fresh:
        _remove_taxdb_files(ingroup_db)
        tmp_tax = cache_dir / "tmp_createtaxdb"
        shutil.rmtree(tmp_tax, ignore_errors=True)
        tmp_tax.mkdir(parents=True, exist_ok=True)
        _run(
            [
                mmseqs,
                "createtaxdb",
                str(ingroup_db),
                str(tmp_tax),
                "--ncbi-tax-dump",
                str(taxdump_dir),
                "--tax-mapping-file",
                str(mapping),
                "--tax-mapping-mode",
                "0",
                "--threads",
                str(threads),
            ]
        )
        shutil.rmtree(tmp_tax, ignore_errors=True)

    return load_taxid_lineages(taxid_lineages), mapped, total


def build_query_type_map(faa_files: Iterable[Path]) -> Dict[str, str]:
    """Map query ID -> original HMM type by scanning only FASTA headers."""
    query_hmm_type: Dict[str, str] = {}

    for faa_file in map(Path, faa_files):
        filename = faa_file.name
        hmm_type = filename
        if filename.endswith(INTERMEDIATE_HITS_SUFFIX):
            hmm_type = filename.removesuffix(INTERMEDIATE_HITS_SUFFIX)
        else:
            hmm_type = faa_file.stem

        for query_id in fasta_ids(faa_file):
            previous = query_hmm_type.get(query_id)
            if previous is not None and previous != hmm_type:
                logger.warning(
                    "Query %s occurs in multiple intermediate types (%s/%s); keeping first",
                    query_id,
                    previous,
                    hmm_type,
                )
                continue
            query_hmm_type.setdefault(query_id, hmm_type)

    return query_hmm_type


def run_search(
    mmseqs: str,
    query_db: Path,
    target_db: Path,
    work_dir: Path,
    label: str,
    threads: int,
    sensitivity: float,
    evalue: float,
    min_identity: float,
    min_qcov: float,
    max_hits: int,
) -> Tuple[Path, Path]:
    """Run one MMseqs search and return its result DB plus TSV.

    Identity, query coverage and E-value thresholds are enforced by MMseqs
    itself so the same alignment DB can be reused directly for LCA.
    """
    result_db = work_dir / (label + "_result")
    tmp_dir = work_dir / ("tmp_" + label)
    output_tsv = work_dir / (label + "_hits.tsv")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _run(
        [
            mmseqs,
            "search",
            str(query_db),
            str(target_db),
            str(result_db),
            str(tmp_dir),
            "--threads",
            str(threads),
            "-s",
            str(sensitivity),
            "-e",
            str(evalue),
            "--min-seq-id",
            str(min_identity),
            "-c",
            str(min_qcov),
            "--cov-mode",
            "2",
            "--max-seqs",
            str(max_hits),
            "--alignment-mode",
            "3",
        ]
    )

    _run(
        [
            mmseqs,
            "convertalis",
            str(query_db),
            str(target_db),
            str(result_db),
            str(output_tsv),
            "--format-output",
            "query,target,pident,evalue,bits,qcov,tcov,alnlen,qlen,tlen",
        ]
    )
    return result_db, output_tsv


def _hit_is_better(candidate: HitCandidate, current: Optional[HitCandidate]) -> bool:
    if current is None:
        return True
    return (
        -candidate.bitscore,
        candidate.evalue,
        -candidate.pident,
        -candidate.qcov,
    ) < (
        -current.bitscore,
        current.evalue,
        -current.pident,
        -current.qcov,
    )


def accumulate_search_hits(
    path: Union[str, Path],
    summaries: Dict[str, QueryHitSummary],
    annotations: Dict[str, str],
    group: str,
) -> None:
    """Stream MMseqs TSV and retain only hits needed for final decisions.

    All queries represented in the search output are processed. For each query
    we retain only:
      - best ingroup hit
      - best outgroup hit
      - best ingroup hit per annotated protein type

    This preserves the original classification logic while avoiding storage of
    every raw alignment hit in RAM.
    """
    is_ingroup = group == "ingroup"

    with Path(path).open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 10:
                if raw.strip():
                    logger.warning("Malformed MMseqs row %d in %s", line_no, path)
                continue

            query, target = fields[0], fields[1]
            hit = HitCandidate(
                target=target,
                protein_type=annotations.get(target) if is_ingroup else None,
                pident=float(fields[2]) / 100.0,
                evalue=float(fields[3]),
                bitscore=float(fields[4]),
                qcov=float(fields[5]),
                tcov=float(fields[6]),
            )

            summary = summaries.get(query)
            if summary is None:
                summary = QueryHitSummary()
                summaries[query] = summary

            if is_ingroup:
                if _hit_is_better(hit, summary.best_ingroup):
                    summary.best_ingroup = hit

                if hit.protein_type:
                    previous = summary.best_by_type.get(hit.protein_type)
                    if previous is None or hit.bitscore > previous.bitscore:
                        summary.best_by_type[hit.protein_type] = hit
            elif _hit_is_better(hit, summary.best_outgroup):
                summary.best_outgroup = hit


def _score_margin(best: float, second: float) -> float:
    return 0.0 if best <= 0 else (best - second) / best


def classify_queries(
    query_hmm_type: Dict[str, str],
    hit_summaries: Dict[str, QueryHitSummary],
    package_has_outgroup: bool,
    ingroup_min_identity: float,
    valid_hit_identity: float,
    domain_reassign_identity: float,
    functional_margin: float,
) -> Dict[str, QueryAssignment]:
    """Classify every intermediate query using compact hit summaries."""
    assignments: Dict[str, QueryAssignment] = {}

    for query, original_type in query_hmm_type.items():
        summary = hit_summaries.get(query)
        best_in = summary.best_ingroup if summary is not None else None

        if best_in is None:
            assignments[query] = QueryAssignment(
                query=query,
                original_hmm_type=original_type,
                status="no_ingroup_hit",
                ingroup_assignment=False,
                valid_hit=False,
                final_type=original_type,
                function_status="not_evaluated",
            )
            continue

        best_out = summary.best_outgroup
        if best_in.pident < ingroup_min_identity:
            ingroup_assignment = False
            status = "rejected_low_ingroup_identity"
        elif (
            package_has_outgroup
            and best_out is not None
            and best_out.bitscore >= best_in.bitscore
        ):
            ingroup_assignment = False
            status = "rejected_outgroup_closer"
        else:
            ingroup_assignment = True
            if best_out is not None:
                status = "accepted_ingroup"
            elif package_has_outgroup:
                status = "accepted_ingroup_no_outgroup_alignment"
            else:
                status = "accepted_ingroup_no_outgroup"

        margin = (
            _score_margin(best_in.bitscore, best_out.bitscore)
            if best_out is not None
            else None
        )
        valid_hit = ingroup_assignment and best_in.pident >= valid_hit_identity

        assignment = QueryAssignment(
            query=query,
            original_hmm_type=original_type,
            status=status,
            ingroup_assignment=ingroup_assignment,
            valid_hit=valid_hit,
            final_type=original_type,
            function_status="pending" if ingroup_assignment else "not_evaluated",
            best_ingroup_reference=best_in.target,
            best_ingroup_type=best_in.protein_type,
            best_ingroup_bitscore=best_in.bitscore,
            best_ingroup_identity=best_in.pident,
            best_ingroup_qcov=best_in.qcov,
            best_ingroup_tcov=best_in.tcov,
            best_outgroup_reference=best_out.target if best_out else None,
            best_outgroup_bitscore=best_out.bitscore if best_out else None,
            ingroup_outgroup_margin=margin,
        )

        if not ingroup_assignment:
            assignments[query] = assignment
            continue

        ranked_types = sorted(
            summary.best_by_type.items(),
            key=lambda item: (-item[1].bitscore, item[1].evalue, -item[1].pident),
        )

        if not ranked_types:
            assignment.function_status = "functional_uncertain"
            assignments[query] = assignment
            continue

        winner_type, winner_hit = ranked_types[0]
        if winner_hit.pident < domain_reassign_identity:
            assignment.function_status = "functional_uncertain"
            assignments[query] = assignment
            continue

        if len(ranked_types) == 1:
            assignment.final_type = winner_type
            assignment.function_status = (
                "same_function_high_identity"
                if winner_type == original_type
                else "single_ingroup_type"
            )
            assignments[query] = assignment
            continue

        second_type, second_hit = ranked_types[1]
        fmargin = _score_margin(winner_hit.bitscore, second_hit.bitscore)
        assignment.second_type = second_type
        assignment.second_type_bitscore = second_hit.bitscore
        assignment.functional_margin = fmargin

        if fmargin >= functional_margin:
            assignment.final_type = winner_type
            assignment.function_status = (
                "same_function_high_identity"
                if winner_type == original_type
                else "function_confident"
            )
        else:
            assignment.function_status = "functional_uncertain"

        assignments[query] = assignment

    return assignments


def filter_intermediate_fastas_by_assignments(
    faa_files: Iterable[Path],
    assignments: Dict[str, QueryAssignment],
    hit_summaries: Dict[str, QueryHitSummary],
) -> Tuple[int, int]:
    """Remove Outgroup-like queries from intermediate FASTA files in one linear pass.

    No additional MMseqs search, database, sort step, or global reject set is
    created. The already existing classification result is reused directly.

    Removed:
      - status == "rejected_outgroup_closer"
      - no Ingroup hit but an Outgroup hit exists

    Kept:
      - accepted Ingroup queries
      - low-Ingroup-identity queries without evidence that Outgroup is closer
      - queries without either Ingroup or Outgroup MMseqs alignment

    The 90% valid-hit rule is unchanged and independent of this FASTA filter.
    """
    kept_total = 0
    removed_total = 0

    for faa_file in map(Path, faa_files):
        tmp_file = faa_file.with_name(faa_file.name + ".tmp_outgroup_filter")
        kept = 0
        removed = 0
        write_record = False

        try:
            with (
                faa_file.open("r", encoding="utf-8") as source,
                tmp_file.open("w", encoding="utf-8") as target,
            ):
                for raw in source:
                    if raw.startswith(">"):
                        query_id = raw[1:].split()[0].strip()
                        assignment = assignments.get(query_id)
                        summary = hit_summaries.get(query_id)

                        outgroup_only = (
                            assignment is not None
                            and assignment.status == "no_ingroup_hit"
                            and summary is not None
                            and summary.best_outgroup is not None
                        )
                        reject = assignment is not None and (
                            assignment.status == "rejected_outgroup_closer"
                            or outgroup_only
                        )

                        write_record = not reject
                        if write_record:
                            kept += 1
                            target.write(raw)
                        else:
                            removed += 1
                        continue

                    if write_record:
                        target.write(raw)

            tmp_file.replace(faa_file)

        except Exception:
            tmp_file.unlink(missing_ok=True)
            raise

        kept_total += kept
        removed_total += removed
        logger.debug(
            "%s: Outgroup filter kept=%d removed=%d",
            faa_file.name,
            kept,
            removed,
        )

    return kept_total, removed_total


def update_database_functional_assignments(
    database_path: str,
    assignments: Dict[str, QueryAssignment],
) -> Tuple[int, int]:
    """Apply high-confidence MMseqs functional evidence to SQLite.

    ``valid_hit`` and domain reassignment are intentionally independent:
    - valid_hit=True requires the configured very-high identity threshold
      (90% by default);
    - a domain/function reassignment requires >=70% identity plus the
      functional ambiguity check performed in ``classify_queries``.
    """
    valid_protein_ids: Set[str] = set()
    domain_updates: Dict[str, Set[Tuple[str, int, int]]] = defaultdict(set)

    for query, assignment in assignments.items():
        try:
            protein_id, start, end = database.parse_domain_header(query)
        except ValueError:
            logger.warning("Invalid UNISUL domain header: %s", query)
            continue

        if assignment.valid_hit:
            valid_protein_ids.add(protein_id)

        if (
            assignment.ingroup_assignment
            and assignment.function_status
            in {
                "function_confident",
                "single_ingroup_type",
            }
            and assignment.final_type != assignment.original_hmm_type
        ):
            domain_updates[assignment.final_type].add((protein_id, start, end))

    promoted = set_valid_hit_true_for_protein_ids(database_path, valid_protein_ids)
    updated_domains = 0
    for protein_type, domain_hits in domain_updates.items():
        updated_domains += database.update_domains(
            database=database_path,
            domain_hits=domain_hits,
            new_domain=protein_type,
        )
    return promoted, updated_domains


def run_mmseqs_taxonomy(
    mmseqs: str,
    query_db: Path,
    ingroup_db: Path,
    work_dir: Path,
    output_tsv: Path,
    threads: int,
    lca_mode: int,
    sensitivity: float,
    evalue: float,
) -> Path:
    """Run the MMseqs taxonomy workflow using true 2bLCA when lca_mode=2.

    The existing query DB and taxonomically annotated ingroup DB are reused
    directly. No FASTA files or reference databases are copied/rebuilt here.
    """
    taxonomy_result = work_dir / "taxonomy_result"
    tmp_dir = work_dir / "tmp_taxonomy"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _run(
        [
            mmseqs,
            "taxonomy",
            str(query_db),
            str(ingroup_db),
            str(taxonomy_result),
            str(tmp_dir),
            "--lca-mode",
            str(lca_mode),
            "--tax-lineage",
            "1",
            "--blacklist",
            "",
            "-s",
            str(sensitivity),
            "-e",
            str(evalue),
            "--threads",
            str(threads),
        ]
    )

    _run(
        [
            mmseqs,
            "createtsv",
            str(query_db),
            str(taxonomy_result),
            str(output_tsv),
        ]
    )
    return output_tsv


def parse_mmseqs_taxonomy(
    taxonomy_tsv: Union[str, Path],
    taxid_to_lineage: Dict[int, str],
    assignments: Dict[str, QueryAssignment],
) -> Tuple[Dict[str, str], Dict[str, Tuple[int, str, str]]]:
    """
    Parse MMseqs taxonomy output and convert assigned local taxids back to
    UNISUL/GTDB lineages.

    Expected MMseqs TSV columns:
        query, taxid, rank, taxname[, lineage]

    The package-local MMseqs taxonomy database uses synthetic taxids generated
    from GTDB lineages. Therefore ``taxid_to_lineage`` is the authoritative
    source for the lineage stored in UNISUL.

    Example:
        taxid 8446
        -> d__Bacteria;p__Desulfobacterota;c__Desulfobacteria;...

    The lineage emitted directly by MMseqs is used only as a fallback and is
    normalized from e.g.

        d_Bacteria;p_Desulfobacterota

    to

        d__Bacteria;p__Desulfobacterota

    Only queries already accepted as ingroup assignments are considered. Queries assigned only to taxid 0 (unclassified) or taxid 1 (artificial root) are not stored as taxonomic assignments.
    """
    domain_lineages: Dict[str, str] = {}
    details: Dict[str, Tuple[int, str, str]] = {}

    valid_prefixes = {"d", "k", "p", "c", "o", "f", "g", "s"}

    with Path(taxonomy_tsv).open(encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                logger.warning(
                    "Malformed MMseqs taxonomy row %d in %s",
                    line_no,
                    taxonomy_tsv,
                )
                continue

            query = fields[0].strip()

            assignment = assignments.get(query)
            if assignment is None or not assignment.ingroup_assignment:
                continue

            try:
                taxid = int(fields[1])
            except ValueError:
                logger.warning(
                    "Invalid MMseqs taxid on line %d: %s",
                    line_no,
                    fields[1],
                )
                continue

            # 0 = unclassified
            # 1 = artificial root of the package-local taxonomy
            if taxid in {0, 1}:
                continue

            rank = fields[2].strip()
            tax_name = fields[3].strip()

            # ---------------------------------------------------------
            # Preferred source:
            # package-local taxid -> original GTDB lineage
            # ---------------------------------------------------------
            lineage = (taxid_to_lineage.get(taxid) or "").strip()

            # ---------------------------------------------------------
            # Fallback:
            # use lineage reported by MMseqs itself.
            #
            # MMseqs may emit:
            #   d_Bacteria;p_Desulfobacterota;...
            #
            # UNISUL expects:
            #   d__Bacteria;p__Desulfobacterota;...
            # ---------------------------------------------------------
            if not lineage and len(fields) >= 5:
                raw_lineage = fields[4].strip()

                if raw_lineage.lower() not in {
                    "",
                    "unclassified",
                    "unknown",
                    "na",
                    "n/a",
                    "none",
                    "-",
                }:
                    lineage_parts = []

                    for token in raw_lineage.split(";"):
                        token = token.strip()

                        if not token:
                            continue

                        # Already correct GTDB/UNISUL format.
                        if "__" in token:
                            prefix, name = token.split("__", 1)
                            prefix = prefix.strip()
                            name = name.strip()

                            if prefix and name and prefix[0] in valid_prefixes:
                                lineage_parts.append(f"{prefix[0]}__{name}")
                            continue

                        # MMseqs representation:
                        # d_Bacteria -> d__Bacteria
                        if (
                            len(token) >= 3
                            and token[0] in valid_prefixes
                            and token[1] == "_"
                        ):
                            name = token[2:].strip()

                            if name:
                                lineage_parts.append(f"{token[0]}__{name}")

                    lineage = ";".join(lineage_parts)

            if not lineage:
                logger.warning(
                    "No usable lineage available for MMseqs taxid %d, query %s",
                    taxid,
                    query,
                )
                continue

            # Final safeguard: only store lineages can actually parse.
            if not any(
                token.strip().startswith(
                    ("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
                )
                for token in lineage.split(";")
            ):
                logger.warning(
                    "Ignoring malformed MMseqs lineage for query %s: %s",
                    query,
                    lineage,
                )
                continue

            domain_lineages[query] = lineage
            details[query] = (
                taxid,
                rank,
                tax_name,
            )

    return domain_lineages, details


def store_mmseqs_taxonomy(
    database_path: str,
    domain_lineages: Dict[str, str],
) -> int:
    """Store complete MMseqs taxonomic lineages separately from EPA-ng."""
    return database.insert_domain_placements(
        database=database_path,
        domain_lineages=domain_lineages,
        method="mmseqs",
    )


def write_assignment_table(
    assignments: Dict[str, QueryAssignment],
    output_file: Union[str, Path],
) -> None:
    columns = [
        "query",
        "original_hmm_type",
        "status",
        "ingroup_assignment",
        "valid_hit",
        "final_type",
        "function_status",
        "best_ingroup_reference",
        "best_ingroup_type",
        "best_ingroup_bitscore",
        "best_ingroup_identity",
        "best_ingroup_qcov",
        "best_ingroup_tcov",
        "best_outgroup_reference",
        "best_outgroup_bitscore",
        "ingroup_outgroup_margin",
        "second_type",
        "second_type_bitscore",
        "functional_margin",
        "taxid",
        "tax_rank",
        "tax_name",
        "taxonomy",
    ]

    with Path(output_file).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for query, a in assignments.items():
            writer.writerow(
                {
                    "query": a.query,
                    "original_hmm_type": a.original_hmm_type,
                    "status": a.status,
                    "ingroup_assignment": int(a.ingroup_assignment),
                    "valid_hit": int(a.valid_hit),
                    "final_type": a.final_type,
                    "function_status": a.function_status,
                    "best_ingroup_reference": a.best_ingroup_reference or "",
                    "best_ingroup_type": a.best_ingroup_type or "",
                    "best_ingroup_bitscore": ""
                    if a.best_ingroup_bitscore is None
                    else a.best_ingroup_bitscore,
                    "best_ingroup_identity": ""
                    if a.best_ingroup_identity is None
                    else a.best_ingroup_identity,
                    "best_ingroup_qcov": ""
                    if a.best_ingroup_qcov is None
                    else a.best_ingroup_qcov,
                    "best_ingroup_tcov": ""
                    if a.best_ingroup_tcov is None
                    else a.best_ingroup_tcov,
                    "best_outgroup_reference": a.best_outgroup_reference or "",
                    "best_outgroup_bitscore": ""
                    if a.best_outgroup_bitscore is None
                    else a.best_outgroup_bitscore,
                    "ingroup_outgroup_margin": ""
                    if a.ingroup_outgroup_margin is None
                    else a.ingroup_outgroup_margin,
                    "second_type": a.second_type or "",
                    "second_type_bitscore": ""
                    if a.second_type_bitscore is None
                    else a.second_type_bitscore,
                    "functional_margin": ""
                    if a.functional_margin is None
                    else a.functional_margin,
                    "taxid": "" if a.taxid is None else a.taxid,
                    "tax_rank": a.tax_rank or "",
                    "tax_name": a.tax_name or "",
                    "taxonomy": a.taxonomy or "",
                }
            )


def process_package(
    config,
    package_dir: Path,
    faa_files: List[Path],
    output_root: Path,
    cache_root: Path,
) -> Dict[str, Union[int, str]]:
    """Process one placement package without copying any FASTA reference data."""
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        raise RuntimeError("MMseqs2 executable 'mmseqs' not found on PATH")

    package_output = output_root / package_dir.name
    package_cache = cache_root / package_dir.name
    package_output.mkdir(parents=True, exist_ok=True)
    package_cache.mkdir(parents=True, exist_ok=True)

    ingroup_faa = package_dir / "ingroup.faa"
    outgroup_faa = package_dir / "outgroup.faa"
    annotations_file = package_dir / "annotations.tsv"
    taxonomy_file = package_dir / "reference_taxonomy.tsv"

    if not ingroup_faa.is_file():
        raise FileNotFoundError("Missing ingroup FASTA: %s" % ingroup_faa)
    if not annotations_file.is_file():
        raise FileNotFoundError("Missing annotations: %s" % annotations_file)
    if not taxonomy_file.is_file():
        raise FileNotFoundError("Missing reference taxonomy: %s" % taxonomy_file)

    sensitivity = float(getattr(config, "mmseqs_sensitivity", 7.5))
    search_evalue = float(getattr(config, "mmseqs_evalue", 1e-5))
    # Search is deliberately more permissive than the biological decisions so
    # short metagenomic fragments are not discarded before classification.
    search_min_identity = _safe_fraction(getattr(config, "mmseqs_min_identity", 0.20))
    search_min_qcov = _safe_fraction(getattr(config, "mmseqs_min_qcov", 0.50))
    max_hits = int(getattr(config, "mmseqs_max_hits", 100))

    # Decision thresholds.
    ingroup_min_identity = _safe_fraction(
        getattr(config, "mmseqs_ingroup_min_identity", 0.30)
    )
    domain_reassign_identity = _safe_fraction(
        getattr(config, "mmseqs_domain_reassign_identity", 0.70)
    )
    valid_hit_identity = _safe_fraction(
        getattr(config, "mmseqs_valid_hit_identity", 0.90)
    )
    functional_margin = _safe_fraction(getattr(config, "mmseqs_function_margin", 0.05))

    # Match the validated benchmark workflow: non-approximated 2bLCA by default.
    taxonomy_lca_mode = int(getattr(config, "mmseqs_lca_mode", 2))
    if taxonomy_lca_mode not in {1, 2, 3, 4}:
        raise ValueError("mmseqs_lca_mode must be 1, 2, 3, or 4")

    annotations = load_annotations(annotations_file)
    taxonomy = load_taxonomy(taxonomy_file)

    ingroup_db = ensure_sequence_db(
        mmseqs=mmseqs,
        fasta_files=[ingroup_faa],
        db_prefix=package_cache / "ingroup_db",
    )
    ingroup_ids = load_mmseqs_lookup(ingroup_db)

    outgroup_db: Optional[Path] = None
    outgroup_ids: List[str] = []
    if outgroup_faa.is_file() and outgroup_faa.stat().st_size > 0:
        outgroup_db = ensure_sequence_db(
            mmseqs=mmseqs,
            fasta_files=[outgroup_faa],
            db_prefix=package_cache / "outgroup_db",
        )
        outgroup_ids = load_mmseqs_lookup(outgroup_db)

    taxid_to_lineage, mapped_tax_refs, total_tax_refs = ensure_taxonomy_db(
        mmseqs=mmseqs,
        ingroup_db=ingroup_db,
        ingroup_ids=ingroup_ids,
        taxonomy=taxonomy,
        taxonomy_source=taxonomy_file,
        cache_dir=package_cache,
        threads=config.cores,
    )

    query_hmm_type = build_query_type_map(faa_files)
    if not query_hmm_type:
        logger.debug("%s: no intermediate queries", package_dir.name)
        return {
            "package": package_dir.name,
            "queries": 0,
            "valid_hits": 0,
            "domain_updates": 0,
            "taxonomy_updates": 0,
        }

    tmp_base = os.environ.get("TMPDIR")
    with tempfile.TemporaryDirectory(
        prefix="unisul_mmseqs_%s_" % package_dir.name,
        dir=tmp_base if tmp_base and Path(tmp_base).is_dir() else None,
    ) as tmp_name:
        work_dir = Path(tmp_name)
        query_db = ensure_sequence_db(
            mmseqs=mmseqs,
            fasta_files=list(map(Path, faa_files)),
            db_prefix=work_dir / "query_db",
        )

        ingroup_result_db, ingroup_tsv = run_search(
            mmseqs=mmseqs,
            query_db=query_db,
            target_db=ingroup_db,
            work_dir=work_dir,
            label="ingroup",
            threads=config.cores,
            sensitivity=sensitivity,
            evalue=search_evalue,
            min_identity=search_min_identity,
            min_qcov=search_min_qcov,
            max_hits=max_hits,
        )
        hit_summaries: Dict[str, QueryHitSummary] = {}
        accumulate_search_hits(
            path=ingroup_tsv,
            summaries=hit_summaries,
            annotations=annotations,
            group="ingroup",
        )
        _remove_db_prefix(ingroup_result_db)
        ingroup_tsv.unlink(missing_ok=True)
        shutil.rmtree(work_dir / "tmp_ingroup", ignore_errors=True)

        if outgroup_db is not None:
            outgroup_result_db, outgroup_tsv = run_search(
                mmseqs=mmseqs,
                query_db=query_db,
                target_db=outgroup_db,
                work_dir=work_dir,
                label="outgroup",
                threads=config.cores,
                sensitivity=sensitivity,
                evalue=search_evalue,
                min_identity=search_min_identity,
                min_qcov=search_min_qcov,
                max_hits=max_hits,
            )
            accumulate_search_hits(
                path=outgroup_tsv,
                summaries=hit_summaries,
                annotations=annotations,
                group="outgroup",
            )
            _remove_db_prefix(outgroup_result_db)
            outgroup_tsv.unlink(missing_ok=True)
            shutil.rmtree(work_dir / "tmp_outgroup", ignore_errors=True)

        assignments = classify_queries(
            query_hmm_type=query_hmm_type,
            hit_summaries=hit_summaries,
            package_has_outgroup=outgroup_db is not None,
            ingroup_min_identity=ingroup_min_identity,
            valid_hit_identity=valid_hit_identity,
            domain_reassign_identity=domain_reassign_identity,
            functional_margin=functional_margin,
        )

        # Reuse the existing competitive Ingroup/Outgroup classification.
        # No second MMseqs run is needed. Outgroup-like sequences are removed
        # from the intermediate FASTAs before downstream taxonomy/placement.
        intermediate_kept, intermediate_removed_outgroup = (
            filter_intermediate_fastas_by_assignments(
                faa_files=faa_files,
                assignments=assignments,
                hit_summaries=hit_summaries,
            )
        )

        promoted, domain_updates = update_database_functional_assignments(
            database_path=config.database_directory,
            assignments=assignments,
        )

        ingroup_assignment_count = sum(
            1 for assignment in assignments.values() if assignment.ingroup_assignment
        )
        valid_hit_count = sum(
            1 for assignment in assignments.values() if assignment.valid_hit
        )

        taxonomy_updates = 0
        taxonomy_tsv = package_output / "taxonomy.tsv"
        if ingroup_assignment_count:
            # Use the validated MMseqs taxonomy workflow with true 2bLCA.
            # query_db and the existing taxonomically annotated ingroup_db are
            # reused directly; no FASTA/reference copies are created.
            run_mmseqs_taxonomy(
                mmseqs=mmseqs,
                query_db=query_db,
                ingroup_db=ingroup_db,
                work_dir=work_dir,
                output_tsv=taxonomy_tsv,
                threads=config.cores,
                lca_mode=taxonomy_lca_mode,
                sensitivity=sensitivity,
                evalue=search_evalue,
            )

            domain_lineages, tax_details = parse_mmseqs_taxonomy(
                taxonomy_tsv=taxonomy_tsv,
                taxid_to_lineage=taxid_to_lineage,
                assignments=assignments,
            )
            _remove_db_prefix(work_dir / "taxonomy_result")
            shutil.rmtree(work_dir / "tmp_taxonomy", ignore_errors=True)

            for query, lineage in domain_lineages.items():
                assignment = assignments.get(query)
                if assignment is None:
                    continue
                assignment.taxonomy = lineage
                detail = tax_details.get(query)
                if detail is not None:
                    assignment.taxid, assignment.tax_rank, assignment.tax_name = detail

            taxonomy_updates = store_mmseqs_taxonomy(
                database_path=config.database_directory,
                domain_lineages=domain_lineages,
            )
        else:
            taxonomy_tsv.unlink(missing_ok=True)

    write_assignment_table(assignments, package_output / "assignments.tsv")

    summary: Dict[str, Union[int, str]] = {
        "package": package_dir.name,
        "queries": len(query_hmm_type),
        "intermediate_hits_kept": intermediate_kept,
        "intermediate_hits_removed_outgroup": intermediate_removed_outgroup,
        "valid_hits": valid_hit_count,
        "ingroup_assignments": ingroup_assignment_count,
        "promoted_proteins": promoted,
        "domain_updates": domain_updates,
        "taxonomy_updates": taxonomy_updates,
        "ingroup_references": len(ingroup_ids),
        "outgroup_references": len(outgroup_ids),
        "taxonomically_mapped_references": mapped_tax_refs,
        "taxonomic_reference_total": total_tax_refs,
        "mmseqs_lca_mode": taxonomy_lca_mode,
    }

    with (package_output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    logger.info(
        "[%s] queries=%d filtered_outgroup=%d valid=%d domain_updates=%d taxonomy_updates=%d",
        package_dir.name,
        len(query_hmm_type),
        intermediate_removed_outgroup,
        valid_hit_count,
        domain_updates,
        taxonomy_updates,
    )
    return summary


def mmseqs_placement_package_reference_assignment(
    config,
) -> List[Dict[str, Union[int, str]]]:
    """Run MMseqs2 assignment for all matching placement packages.

    MMseqs reference assignment is independent of whether EPA-ng placement
    was requested. Matching packages are selected automatically from the
    intermediate-hit domain files.

    Failure of one package does not abort processing of the remaining
    packages.
    """

    if shutil.which("mmseqs") is None:
        logger.warning("MMseqs2 is not available on PATH; skipping MMseqs assignment")
        return []

    intermediate_files = myUtil.collect_intermediate_hits(config)

    if not intermediate_files:
        logger.info("No intermediate hits available for MMseqs assignment")
        return []

    output_root = (
        Path(config.cross_check_directory) / "mmseqs_reference_sequence_assignment"
    )
    cache_root = output_root / "_reference_cache"

    output_root.mkdir(parents=True, exist_ok=True)
    cache_root.mkdir(parents=True, exist_ok=True)

    packages = list(
        find_matching_placement_packages(
            config.paths.placement,
            intermediate_files,
        )
    )

    if not packages:
        logger.info(
            "No matching placement packages found for MMseqs reference assignment"
        )
        return []

    summaries: List[Dict[str, Union[int, str]]] = []
    failed_packages: List[str] = []

    for package_dir, faa_files in packages:
        package_dir = Path(package_dir)

        try:
            summaries.append(
                process_package(
                    config=config,
                    package_dir=package_dir,
                    faa_files=[Path(path) for path in faa_files],
                    output_root=output_root,
                    cache_root=cache_root,
                )
            )

        except Exception as exc:
            failed_packages.append(package_dir.name)

            logger.exception(
                "MMseqs assignment failed for package %s; "
                "continuing with next package: %s",
                package_dir.name,
                exc,
            )
            continue

    if failed_packages:
        logger.warning(
            "Keeping MMseqs reference cache because %d package(s) failed: %s",
            len(failed_packages),
            ", ".join(failed_packages[:20]),
        )
    else:
        shutil.rmtree(
            cache_root,
            ignore_errors=True,
        )
        logger.debug(
            "Removed MMseqs reference cache after successful completion: %s",
            cache_root,
        )

    return summaries

# src/hmsss/stages/individual_reports.py
# Streaming replacement for individual genome report writing.

from __future__ import annotations

import os
import sqlite3
import csv
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, List, Set, TextIO

from hmsss.core.logging import get_logger
from hmsss.parse_reports import parse_reports

logger = get_logger(__name__)

TAXON_COLS = [
    "genomeID",
    "Superkingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]


def _db_uri_ro_immutable(db_file: str) -> str:
    """Return a read-only immutable SQLite URI."""
    abs_db = os.path.abspath(db_file)
    return f"file:{abs_db}?mode=ro&immutable=1"


def _has_column(cur: sqlite3.Cursor, table_name: str, column_name: str) -> bool:
    """Return True if table_name contains column_name."""
    try:
        cur.execute(f"PRAGMA table_info({table_name});")
        return any(row[1] == column_name for row in cur.fetchall())
    except sqlite3.Error as exc:
        logger.warning("Could not inspect column %s.%s: %s", table_name, column_name, exc)
        return False


def _count_total_reports(
        cur: sqlite3.Cursor,
        *,
        write_empty: bool,
        use_non_valid_hits: bool,
        valid_hit_column_available: bool,
) -> int:
    """Return the number of genome report files expected to be written."""
    if write_empty:
        cur.execute("SELECT COUNT(*) AS n FROM Genomes;")
        row = cur.fetchone()
        return int(row["n"] if row is not None else 0)

    valid_where = ""
    if not use_non_valid_hits and valid_hit_column_available:
        valid_where = "WHERE p.valid_hit = 1"

    cur.execute(f"""
        SELECT COUNT(DISTINCT p.genomeID) AS n
        FROM Proteins p
        JOIN Domains d ON d.proteinID = p.proteinID
        {valid_where};
    """)
    row = cur.fetchone()
    return int(row["n"] if row is not None else 0)


def _log_report_progress(
        *,
        written_reports: int,
        total_reports: int,
        n_rows: int,
        force: bool = False,
        log_every: int = 10000,
) -> None:
    """Log report-writing progress including percentage."""
    if total_reports <= 0:
        return

    if not force and written_reports % log_every != 0:
        return

    percent = (written_reports / total_reports) * 100.0
    msg = (
        f"Wrote {written_reports}/{total_reports} genome reports "
        f"({percent:.2f}%); streamed {n_rows} domain rows"
    )
    logger.info(msg)
    # print(msg, flush=True)


def _protein_from_row(row: sqlite3.Row) -> parse_reports.Protein:
    """Create a Protein object from one protein-domain SQL row.

    Important: p.comment is stored as the domain-level selection comment.
    The downstream output routine rebuilds protein.selection_comment from
    Domain.selection_comment_list, so the comment must be passed into the
    Protein constructor and not only assigned to protein.selection_comment.
    """
    protein = parse_reports.Protein(
        row["proteinID"],
        row["domain"],
        row["domStart"],
        row["domEnd"],
        row["score"],
        selection_comment=row["comment"] or "",
    )
    protein.genomeID = row["genomeID"] or ""
    protein.clusterID = row["clusterID"] or ""
    protein.gene_contig = row["contig"] or ""
    protein.gene_start = row["gene_start"] or 0
    protein.gene_end = row["gene_end"] or 0
    protein.gene_strand = row["gene_strand"] or "."
    protein.gene_locustag = row["locustag"] or ""
    protein.protein_sequence = row["protein_sequence"] or ""
    protein.selection_comment = row["comment"] or ""
    protein.alternative_hit = row["alternative_hit"] or ""
    return protein


def _taxon_from_row(row: sqlite3.Row) -> Dict[str, str]:
    """Build a taxonomy record compatible with parse_reports.output_genome_report()."""
    return {
        "genomeID": row["genomeID"] or "",
        "Superkingdom": row["Superkingdom"] or "",
        "Phylum": row["Phylum"] or "",
        "Class": row["Class"] or "",
        "Order": row["Ordnung"] or "",  # DB column is Ordnung; report column is Order.
        "Family": row["Family"] or "",
        "Genus": row["Genus"] or "",
        "Species": row["Species"] or "",
    }


def _empty_taxon_from_genomes_row(row: sqlite3.Row) -> Dict[str, str]:
    """Build a taxonomy record from a Genomes-only row."""
    return {
        "genomeID": row["genomeID"] or "",
        "Superkingdom": row["Superkingdom"] or "",
        "Phylum": row["Phylum"] or "",
        "Class": row["Class"] or "",
        "Order": row["Ordnung"] or "",
        "Family": row["Family"] or "",
        "Genus": row["Genus"] or "",
        "Species": row["Species"] or "",
    }


PATHWAY_INPUT_COL = "input_sulfur_species_or_compound"
PATHWAY_OUTPUT_COL = "output_sulfur_species_or_compound"
PATHWAY_ENZYME_COL = "enzyme_abbreviation"
PATHWAY_PATHWAY_COL = "pathway"


def _split_enzyme_field(value: str) -> Set[str]:
    if value is None:
        return set()
    return {x.strip() for x in str(value).split(";") if x.strip()}


def _normalize_domain_name(domain: str) -> str:
    domain = str(domain).strip()
    if not domain:
        return ""
    return domain.rsplit("_", 1)[-1]


def _protein_domain_set(protein_dict: Dict[str, parse_reports.Protein]) -> Set[str]:
    domains: Set[str] = set()
    for protein in protein_dict.values():
        dom_string = protein.get_domains()
        if not dom_string:
            continue
        for part in str(dom_string).split("-"):
            part = _normalize_domain_name(part)
            if part:
                domains.add(part)
    return domains


def _load_pathway_rows(pathway_file: str | Path | None) -> List[Dict[str, str]]:
    if not pathway_file:
        return []
    path = Path(pathway_file)
    if not path.is_file():
        raise FileNotFoundError(f"Pathway file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(handle, dialect=dialect)
        required = [PATHWAY_INPUT_COL, PATHWAY_OUTPUT_COL, PATHWAY_ENZYME_COL, PATHWAY_PATHWAY_COL]
        missing = [col for col in required if col not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Pathway file {path} is missing required columns: {missing}. Found: {reader.fieldnames}")
        rows: List[Dict[str, str]] = []
        for row in reader:
            enzymes = _split_enzyme_field(row.get(PATHWAY_ENZYME_COL, ""))
            if not enzymes:
                continue
            rows.append({
                PATHWAY_INPUT_COL: row.get(PATHWAY_INPUT_COL, "") or "",
                PATHWAY_OUTPUT_COL: row.get(PATHWAY_OUTPUT_COL, "") or "",
                PATHWAY_ENZYME_COL: ";".join(sorted(enzymes)),
                PATHWAY_PATHWAY_COL: row.get(PATHWAY_PATHWAY_COL, "") or "",
                "_enzyme_set": enzymes,
            })
    return rows


def _collapse_nested_pathways(calls):
    kept = []

    for i, call in enumerate(calls):
        enzymes_i = call["_enzyme_set"]
        remove = False

        for j, other in enumerate(calls):
            if i == j:
                continue

            enzymes_j = other["_enzyme_set"]

            if enzymes_i < enzymes_j:
                remove = True
                break

        if not remove:
            kept.append(call)

    return kept


def _write_pathway_report_header(writer: TextIO) -> None:
    writer.write("genomeID\tspecies\tinput\toutput\tenzymes\n")


def _write_genome_pathways(*, writer, genome_id, species, protein_dict, pathway_rows) -> int:
    if not pathway_rows or not protein_dict:
        return 0

    present = _protein_domain_set(protein_dict)

    calls = []
    for pathway in pathway_rows:
        required = pathway["_enzyme_set"]
        if required <= present:
            calls.append(pathway)

    calls = _collapse_nested_pathways(calls)

    n_written = 0
    for pathway in calls:
        writer.write("\t".join([
            genome_id,
            species or "",
            pathway[PATHWAY_INPUT_COL],
            pathway[PATHWAY_OUTPUT_COL],
            pathway[PATHWAY_ENZYME_COL],
        ]) + "\n")
        n_written += 1

    return n_written


def _write_one_report(
        *,
        out_dir: Path,
        genome_id: str,
        protein_dict: Dict[str, parse_reports.Protein],
        taxon_rec: Optional[Dict[str, str]],
        pathway_rows: Optional[List[Dict[str, str]]] = None,
        pathway_writer: Optional[TextIO] = None,
        write_individual_report: bool = True,
) -> int:
    """
    Finalize one genome-sized protein dict and write one genome TSV report.

    If pathway_rows and pathway_writer are supplied, also append matching
    pathway calls to the global pathway report.

    Returns
    -------
    int
        Number of pathway rows written for this genome.
    """
    if protein_dict:
        parse_reports.define_best_score_hits_for_protein_dict(protein_dict)
        parse_reports.define_selection_comments_for_protein_dict(protein_dict)

    if write_individual_report:
        out_file = out_dir / f"{genome_id}.tsv"
        parse_reports.output_genome_report(
            output_filepath=str(out_file),
            protein_dict=protein_dict,
            cluster_dict={},
            taxon_dict={genome_id: taxon_rec or {}},
            genomeID="",
            writemode="w",
            taxon_divider="\t",
        )

    n_pathway_rows = 0
    if pathway_rows and pathway_writer is not None:
        species = ""
        if taxon_rec:
            species = taxon_rec.get("Species", "") or ""

        n_pathway_rows = _write_genome_pathways(
            writer=pathway_writer,
            genome_id=genome_id,
            species=species,
            protein_dict=protein_dict,
            pathway_rows=pathway_rows,
        )

    return n_pathway_rows


def _stream_hit_rows(
        cur: sqlite3.Cursor,
        *,
        use_non_valid_hits: bool,
        valid_hit_column_available: bool,
) -> Iterable[sqlite3.Row]:
    """
    Stream all protein-domain rows ordered by genome.

    The ORDER BY is intentional: it makes all rows for a genome contiguous,
    so only one genome has to be held in memory at any time. For performance,
    create these indexes once on the database:

        CREATE INDEX IF NOT EXISTS idx_proteins_report_order
        ON Proteins(genomeID, contig, start, proteinID);

        CREATE INDEX IF NOT EXISTS idx_domains_protein_start
        ON Domains(proteinID, domStart);
    """
    valid_where = ""
    if not use_non_valid_hits and valid_hit_column_available:
        valid_where = "WHERE p.valid_hit = 1"

    sql = f"""
        SELECT
            p.proteinID       AS proteinID,
            p.genomeID        AS genomeID,
            p.clusterID       AS clusterID,
            p.contig          AS contig,
            p.start           AS gene_start,
            p.end             AS gene_end,
            p.strand          AS gene_strand,
            p.locustag        AS locustag,
            p.sequence        AS protein_sequence,
            p.comment         AS comment,
            p.alternative_hit AS alternative_hit,
            d.domain          AS domain,
            d.domStart        AS domStart,
            d.domEnd          AS domEnd,
            d.score           AS score,
            g.Superkingdom    AS Superkingdom,
            g.Phylum          AS Phylum,
            g.Class           AS Class,
            g.Ordnung         AS Ordnung,
            g.Family          AS Family,
            g.Genus           AS Genus,
            g.Species         AS Species
        FROM Proteins p
        JOIN Domains d ON d.proteinID = p.proteinID
        LEFT JOIN Genomes g ON g.genomeID = p.genomeID
        {valid_where}
        ORDER BY p.genomeID, p.contig, p.start;
    """

    cur.execute(sql)
    yield from cur


def _write_empty_reports_for_missing_genomes(
        cur: sqlite3.Cursor,
        *,
        out_dir: Path,
        written_genome_ids: set[str],
        total_reports: int,
        n_rows: int,
        log_every: int = 10000,
        write_individual_reports: bool = True,
) -> int:
    """
    Write empty report files for genomes that had no streamed hit rows.

    This preserves the previous behavior of producing one file per genome.
    Disable this through config.write_empty_genome_reports = False if only genomes
    with hits should be written.
    """
    sql = """
        SELECT
            genomeID,
            Superkingdom,
            Phylum,
            Class,
            Ordnung,
            Family,
            Genus,
            Species
        FROM Genomes
        ORDER BY genomeID;
    """
    cur.execute(sql)

    n_empty = 0
    for row in cur:
        gid = row["genomeID"]
        if gid in written_genome_ids:
            continue

        _write_one_report(
            out_dir=out_dir,
            genome_id=gid,
            protein_dict={},
            taxon_rec=_empty_taxon_from_genomes_row(row),
            write_individual_report=write_individual_reports,
        )
        written_genome_ids.add(gid)
        n_empty += 1

        _log_report_progress(
            written_reports=len(written_genome_ids),
            total_reports=total_reports,
            n_rows=n_rows,
            log_every=log_every,
        )

    return n_empty


def write_individual_genome_reports(config) -> None:
    """
    Streaming replacement for individual genome report writing.

    This function avoids chunked fetches entirely:
      - one cursor streams Proteins JOIN Domains rows ordered by genomeID;
      - only the current genome's proteins are kept in memory;
      - the report is written as soon as the genomeID changes;
      - taxonomy is read via a LEFT JOIN instead of loading the whole taxonomy dict;
      - optional empty reports are written in a second lightweight Genomes stream.

    Optional config attributes:
      - use_non_valid_hits: bool, default False
      - write_empty_genome_reports: bool, default True
      - report_stream_log_every: int, default 10000
        Progress is logged as written/total reports plus percentage.
    """
    out_dir = Path(config.fasta_initial_hit_directory)
    write_individual_reports = not bool(getattr(config, "disable_individual_reports", False))

    if write_individual_reports:
        out_dir.mkdir(parents=True, exist_ok=True)

    pathway_file = getattr(config, "pathway_file", None)
    pathway_report_file = getattr(config, "pathway_report_file", None)
    pathway_rows = _load_pathway_rows(pathway_file) if pathway_file else []
    if pathway_rows and not pathway_report_file:
        pathway_report_file = str(out_dir / "genome_pathway_report.tsv")

    db_path = _db_uri_ro_immutable(config.database_directory)
    use_non_valid = bool(getattr(config, "use_non_valid_hits", False))
    write_empty = bool(getattr(config, "write_empty_genome_reports", False))
    log_every = int(getattr(config, "report_stream_log_every", 10000))

    written_genome_ids: set[str] = set()
    current_gid: Optional[str] = None
    current_taxon: Optional[Dict[str, str]] = None
    current_proteins: Dict[str, parse_reports.Protein] = {}

    n_rows = 0
    n_reports_with_hits = 0
    n_pathway_rows = 0

    if write_individual_reports:
        start_msg = f"Streaming individual genome reports to {out_dir}"
    else:
        start_msg = "Individual genome reports disabled; streaming genomes for pathway report only"
    logger.info(start_msg)
    # print(start_msg, flush=True)

    if pathway_rows:
        msg1 = f"Loaded {len(pathway_rows)} precomputed pathway definitions from {pathway_file}"
        msg2 = f"Writing genome pathway report to {pathway_report_file}"
        logger.info(msg1)
        logger.info(msg2)
        # print(msg1, flush=True)
        # print(msg2, flush=True)

    pathway_handle = None
    try:
        if pathway_rows:
            pathway_report_path = Path(pathway_report_file)
            pathway_report_path.parent.mkdir(parents=True, exist_ok=True)
            pathway_handle = pathway_report_path.open("w", encoding="utf-8")
            _write_pathway_report_header(pathway_handle)

        with sqlite3.connect(db_path, uri=True) as con:
            con.row_factory = sqlite3.Row
            cur = con.cursor()

            cur.execute("PRAGMA foreign_keys = ON;")
            cur.execute("PRAGMA temp_store = MEMORY;")
            cur.execute("PRAGMA cache_size = 200000;")

            valid_col = _has_column(cur, "Proteins", "valid_hit")
            total_reports = _count_total_reports(
                cur,
                write_empty=(write_empty and write_individual_reports),
                use_non_valid_hits=use_non_valid,
                valid_hit_column_available=valid_col,
            )
            if write_individual_reports:
                expected_msg = f"Expected genome reports to write: {total_reports}"
            else:
                expected_msg = f"Expected genomes with hits to stream for pathway report: {total_reports}"
            logger.info(expected_msg)
            # print(expected_msg, flush=True)

            for row in _stream_hit_rows(
                    cur,
                    use_non_valid_hits=use_non_valid,
                    valid_hit_column_available=valid_col,
            ):
                n_rows += 1
                gid = row["genomeID"]
                pid = row["proteinID"]

                if current_gid is None:
                    current_gid = gid
                    current_taxon = _taxon_from_row(row)

                if gid != current_gid:
                    n_pathway_rows += _write_one_report(
                        out_dir=out_dir,
                        genome_id=current_gid,
                        protein_dict=current_proteins,
                        taxon_rec=current_taxon,
                        pathway_rows=pathway_rows,
                        pathway_writer=pathway_handle,
                        write_individual_report=write_individual_reports,
                    )
                    written_genome_ids.add(current_gid)
                    n_reports_with_hits += 1

                    _log_report_progress(
                        written_reports=len(written_genome_ids),
                        total_reports=total_reports,
                        n_rows=n_rows,
                        log_every=log_every,
                    )

                    current_gid = gid
                    current_taxon = _taxon_from_row(row)
                    current_proteins = {}

                protein = current_proteins.get(pid)
                if protein is None:
                    current_proteins[pid] = _protein_from_row(row)
                else:
                    protein.add_domain(
                        row["domain"],
                        row["domStart"],
                        row["domEnd"],
                        row["score"],
                        selection_comment=row["comment"] or "",
                    )

            # Flush final genome with hits.
            if current_gid is not None:
                n_pathway_rows += _write_one_report(
                    out_dir=out_dir,
                    genome_id=current_gid,
                    protein_dict=current_proteins,
                    taxon_rec=current_taxon,
                    pathway_rows=pathway_rows,
                    pathway_writer=pathway_handle,
                    write_individual_report=write_individual_reports,
                )
                written_genome_ids.add(current_gid)
                n_reports_with_hits += 1
                _log_report_progress(
                    written_reports=len(written_genome_ids),
                    total_reports=total_reports,
                    n_rows=n_rows,
                    force=True,
                    log_every=log_every,
                )

            n_empty = 0
            if write_empty and write_individual_reports:
                n_empty = _write_empty_reports_for_missing_genomes(
                    cur,
                    out_dir=out_dir,
                    written_genome_ids=written_genome_ids,
                    total_reports=total_reports,
                    n_rows=n_rows,
                    log_every=log_every,
                    write_individual_reports=write_individual_reports,
                )

    finally:
        if pathway_handle is not None:
            pathway_handle.close()

    if 'n_empty' not in locals():
        n_empty = 0

    if write_individual_reports:
        final_msg = (
            f"Finished streaming genome reports: {len(written_genome_ids)}/"
            f"{total_reports if 'total_reports' in locals() else len(written_genome_ids)} written; "
            f"{n_reports_with_hits} with hits, {n_empty} empty, "
            f"{n_rows} domain rows streamed"
        )
    else:
        final_msg = (
            f"Finished streaming genomes for pathway report: {len(written_genome_ids)}/"
            f"{total_reports if 'total_reports' in locals() else len(written_genome_ids)} genomes with hits processed; "
            f"{n_rows} domain rows streamed"
        )
    logger.info(final_msg)
    # print(final_msg, flush=True)

    if pathway_rows:
        pathway_msg = f"Finished genome pathway report: {n_pathway_rows} pathway rows written"
        logger.info(pathway_msg)
        # print(pathway_msg, flush=True)

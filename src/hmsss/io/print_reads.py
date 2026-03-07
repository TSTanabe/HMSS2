#!/usr/bin/python
from __future__ import annotations

import os
import re
from collections import defaultdict
from typing import Any, Dict, Tuple

from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


def print_read_hit_reports(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
        metagenome_dict: Dict[str, Dict[str, Any]],
        lineage_dict: Dict[str, Dict[str, Any]],
        gpkg_length_dict: Dict[str, int],
) -> None:
    """
    Write all read-based report files.

    Parameters
    ----------
    directory : str
        Output directory for fetch reports.
    read_dict : dict
        Mapping (readID, gpkg_name, metagenomeID) -> Read
    metagenome_dict : dict
        Mapping metagenomeID -> metadata dictionary
    lineage_dict : dict
        Mapping lineageID -> lineage metadata dictionary
    """
    os.makedirs(directory, exist_ok=True)

    logger.info(
        "Printing read hit reports to %s | reads=%s | metagenomes=%s | lineages=%s",
        directory,
        len(read_dict),
        len(metagenome_dict),
        len(lineage_dict),
    )

    output_read_hit_table(
        directory=directory,
        read_dict=read_dict,
        metagenome_dict=metagenome_dict,
        lineage_dict=lineage_dict,
        filename="summary_read_hit_table.txt",
    )

    output_read_hit_details(
        directory=directory,
        read_dict=read_dict,
        metagenome_dict=metagenome_dict,
        lineage_dict=lineage_dict,
        filename="summary_read_hits_detailed.txt",
    )

    output_read_lineage_counts(
        directory=directory,
        read_dict=read_dict,
        metagenome_dict=metagenome_dict,
        lineage_dict=lineage_dict,
        filename="summary_read_lineage_counts.txt",
    )

    output_read_taxonomy_level_counts(
        directory=directory,
        read_dict=read_dict,
        metagenome_dict=metagenome_dict,
        lineage_dict=lineage_dict,
        gpkg_length_dict=gpkg_length_dict,
        filename="summary_read_taxonomy_level_counts.txt",
    )


def _safe_str(value: Any, default: str = "NA") -> str:
    """
    Convert a value to string, replacing None / empty values by a default.
    """
    if value is None:
        return default
    s = str(value).strip()
    return s if s else default


def _build_lineage_string(lineage_meta: Dict[str, Any]) -> str:
    """
    Return a printable lineage string.

    Priority:
    1. raw_lineage
    2. reconstructed lineage from rank fields
    3. NA
    """
    if not lineage_meta:
        return "NA"

    raw = lineage_meta.get("raw_lineage")
    if raw:
        raw = str(raw).strip()
        if raw:
            return raw

    parts = []

    root = _safe_str(lineage_meta.get("root"), default="")
    kingdom = _safe_str(lineage_meta.get("kingdom"), default="")
    phylum = _safe_str(lineage_meta.get("phylum"), default="")
    class_ = _safe_str(lineage_meta.get("class"), default="")
    order = _safe_str(lineage_meta.get("order"), default="")
    family = _safe_str(lineage_meta.get("family"), default="")
    genus = _safe_str(lineage_meta.get("genus"), default="")
    species = _safe_str(lineage_meta.get("species"), default="")

    if root:
        parts.append(root)
    if kingdom:
        parts.append(f"k__{kingdom}")
    if phylum:
        parts.append(f"p__{phylum}")
    if class_:
        parts.append(f"c__{class_}")
    if order:
        parts.append(f"o__{order}")
    if family:
        parts.append(f"f__{family}")
    if genus:
        parts.append(f"g__{genus}")
    if species:
        parts.append(f"s__{species}")

    return "; ".join(parts) if parts else "NA"


def _aggregate_reads_by_metagenome_gpkg_lineage(
        read_dict: Dict[Tuple[str, str, str], Read],
) -> Dict[Tuple[str, str, str], int]:
    """
    Aggregate reads by (metagenomeID, gpkg_name, lineageID).

    Returns
    -------
    dict
        Mapping (metagenomeID, gpkg_name, lineageID) -> count
    """
    counts: Dict[Tuple[str, str, str], int] = defaultdict(int)

    for (_, gpkg_name, metagenomeID), read in read_dict.items():
        lineageID = _safe_str(read.lineageID)
        counts[(metagenomeID, gpkg_name, lineageID)] += 1

    return counts


def output_read_hit_table(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
        metagenome_dict: Dict[str, Dict[str, Any]],
        lineage_dict: Dict[str, Dict[str, Any]],
        filename: str = "summary_read_hit_table.txt",
) -> None:
    """
    Write aggregated read hit table.

    Output columns
    --------------
    metagenomeID
    genomeID
    gpkg_name
    read_count
    lineageID
    lineage
    forward_reads
    reverse_reads
    prokaryotic_fraction
    """
    outpath = os.path.join(directory, filename)
    counts = _aggregate_reads_by_metagenome_gpkg_lineage(read_dict)

    def taxonomy_sort_key(key):
        metagenomeID, gpkg_name, lineageID = key
        lineage_meta = lineage_dict.get(lineageID, {})

        return (
            _safe_str(lineage_meta.get("root"), ""),
            _safe_str(lineage_meta.get("kingdom"), ""),
            _safe_str(lineage_meta.get("phylum"), ""),
            _safe_str(lineage_meta.get("class"), ""),
            _safe_str(lineage_meta.get("order"), ""),
            _safe_str(lineage_meta.get("family"), ""),
            _safe_str(lineage_meta.get("genus"), ""),
            _safe_str(lineage_meta.get("species"), ""),
            _safe_str(metagenomeID),
            _safe_str(gpkg_name),
        )

    sorted_keys = sorted(counts, key=taxonomy_sort_key)
    with open(outpath, "w") as out:
        out.write(
            "metagenomeID\tgenomeID\tgpkg_name\tread_count\tlineage\tforward_reads\treverse_reads\tprokaryotic_fraction\n"
        )

        for metagenomeID, gpkg_name, lineageID in sorted_keys:
            meta = metagenome_dict.get(metagenomeID, {})
            genomeID = _safe_str(meta.get("genomeID"))
            forward_reads = _safe_str(meta.get("forward_reads"))
            reverse_reads = _safe_str(meta.get("reverse_reads"))
            prokaryotic_fraction = _safe_str(meta.get("prokaryotic_fraction"))

            lineage_meta = lineage_dict.get(lineageID, {})
            lineage = _build_lineage_string(lineage_meta)

            out.write(
                f"{_safe_str(metagenomeID)}\t"
                f"{genomeID}\t"
                f"{_safe_str(gpkg_name)}\t"
                f"{counts[(metagenomeID, gpkg_name, lineageID)]}\t"
                f"{lineage}\t"
                f"{forward_reads}\t"
                f"{reverse_reads}\t"
                f"{prokaryotic_fraction}\n"
            )

    logger.info("Wrote aggregated read hit table: %s", outpath)


def output_read_hit_details(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
        metagenome_dict: Dict[str, Dict[str, Any]],
        lineage_dict: Dict[str, Dict[str, Any]],
        filename: str = "summary_read_hits_detailed.txt",
) -> None:
    """
    Write detailed per-read placement table.

    Output columns
    --------------
    metagenomeID
    genomeID
    readID
    gpkg_name
    lineageID
    lineage
    coverage
    dom_start
    dom_end
    sequence_length
    alignment_length
    """
    outpath = os.path.join(directory, filename)

    with open(outpath, "w") as out:
        out.write(
            "metagenomeID\tgenomeID\treadID\tgpkg_name\tlineageID\tlineage\tcoverage\tdom_start\tdom_end\tsequence_length\talignment_length\n"
        )

        for (readID, gpkg_name, metagenomeID), read in sorted(read_dict.items()):
            meta = metagenome_dict.get(metagenomeID, {})
            genomeID = _safe_str(meta.get("genomeID"))

            lineageID = _safe_str(read.lineageID)
            lineage_meta = lineage_dict.get(lineageID, {})
            lineage = _build_lineage_string(lineage_meta)

            sequence_length = len(read.sequence) if getattr(read, "sequence", None) else 0
            alignment_length = len(read.alignment) if getattr(read, "alignment", None) else 0

            out.write(
                f"{_safe_str(metagenomeID)}\t"
                f"{genomeID}\t"
                f"{_safe_str(readID)}\t"
                f"{_safe_str(gpkg_name)}\t"
                f"{lineageID}\t"
                f"{lineage}\t"
                f"{getattr(read, 'coverage', 'NA')}\t"
                f"{getattr(read, 'start', 'NA')}\t"
                f"{getattr(read, 'end', 'NA')}\t"
                f"{sequence_length}\t"
                f"{alignment_length}\n"
            )

    logger.info("Wrote detailed read hit table: %s", outpath)


def output_read_lineage_counts(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
        metagenome_dict: Dict[str, Dict[str, Any]],
        lineage_dict: Dict[str, Dict[str, Any]],
        filename: str = "summary_read_lineage_counts.txt",
) -> None:
    """
    Write lineage summary separated by gpkg/domain type.

    Aggregation level
    -----------------
    (lineageID, gpkg_name)

    Output columns
    --------------
    lineageID
    lineage
    gpkg_name
    read_count
    metagenome_count
    """

    outpath = os.path.join(directory, filename)

    lineage_read_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    lineage_metagenomes: Dict[Tuple[str, str], set[str]] = defaultdict(set)

    # Aggregation
    for (_, gpkg_name, metagenomeID), read in read_dict.items():
        lineageID = read.lineageID if read.lineageID else "NA"

        key = (lineageID, gpkg_name)

        lineage_read_counts[key] += 1
        lineage_metagenomes[key].add(metagenomeID)

    # ---------- taxonomy sort key ----------
    def taxonomy_sort_key(item):
        lineageID, gpkg_name = item

        meta = lineage_dict.get(lineageID, {})

        return (
            str(meta.get("root", "")),
            str(meta.get("kingdom", "")),
            str(meta.get("phylum", "")),
            str(meta.get("class", "")),
            str(meta.get("order", "")),
            str(meta.get("family", "")),
            str(meta.get("genus", "")),
            str(meta.get("species", "")),
            gpkg_name,
        )

    sorted_keys = sorted(lineage_read_counts.keys(), key=taxonomy_sort_key)

    with open(outpath, "w") as out:
        out.write(
            "gpkg_name\tread_count\tmetagenome_count\tlineage\n"
        )

        for lineageID, gpkg_name in sorted_keys:
            lineage_meta = lineage_dict.get(lineageID, {})

            # lineage string
            raw = lineage_meta.get("raw_lineage")
            if raw:
                lineage = raw
            else:
                parts = []

                if lineage_meta.get("root"):
                    parts.append(str(lineage_meta["root"]))
                if lineage_meta.get("kingdom"):
                    parts.append(f"k__{lineage_meta['kingdom']}")
                if lineage_meta.get("phylum"):
                    parts.append(f"p__{lineage_meta['phylum']}")
                if lineage_meta.get("class"):
                    parts.append(f"c__{lineage_meta['class']}")
                if lineage_meta.get("order"):
                    parts.append(f"o__{lineage_meta['order']}")
                if lineage_meta.get("family"):
                    parts.append(f"f__{lineage_meta['family']}")
                if lineage_meta.get("genus"):
                    parts.append(f"g__{lineage_meta['genus']}")
                if lineage_meta.get("species"):
                    parts.append(f"s__{lineage_meta['species']}")

                lineage = "; ".join(parts) if parts else "NA"

            out.write(
                f"{gpkg_name}\t"
                f"{lineage_read_counts[(lineageID, gpkg_name)]}\t"
                f"{len(lineage_metagenomes[(lineageID, gpkg_name)])}\t"
                f"{lineage}\n"
            )

    logger.info("Wrote lineage-by-type summary table: %s", outpath)


def output_read_fastas(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
) -> None:
    """
    Write separate FASTA files per protein/gpkg type.

    Output files
    ------------
    <TYPE>.fna
        DNA sequences for reads assigned to this type

    <TYPE>.faa_aln
        aligned protein sequences for reads assigned to this type
    """

    def sanitize_filename(name: str) -> str:
        """
        Make a filesystem-safe filename stem.
        """
        name = str(name).strip()
        name = re.sub(r"[^\w.\-]+", "_", name)
        return name if name else "NA"

    grouped_reads: Dict[str, list[Tuple[str, str, str, Read]]] = defaultdict(list)

    for (readID, gpkg_name, metagenomeID), read in read_dict.items():
        grouped_reads[gpkg_name].append((readID, gpkg_name, metagenomeID, read))

    total_dna_written = 0
    total_protein_written = 0

    for gpkg_name in sorted(grouped_reads):
        safe_name = sanitize_filename(gpkg_name)

        dna_path = os.path.join(directory, f"{safe_name}.fna")
        protein_path = os.path.join(directory, f"{safe_name}.faa_aln")

        dna_written = 0
        protein_written = 0

        with open(dna_path, "w") as dna_out, open(protein_path, "w") as prot_out:
            for readID, gpkg_name, metagenomeID, read in sorted(grouped_reads[gpkg_name]):
                header = f"{readID}|{metagenomeID}|{gpkg_name}"

                if getattr(read, "sequence", None):
                    dna_out.write(f">{header}\n{read.sequence}\n")
                    dna_written += 1

                if getattr(read, "alignment", None):
                    prot_out.write(f">{header}\n{read.alignment}\n")
                    protein_written += 1

        total_dna_written += dna_written
        total_protein_written += protein_written

        logger.info(
            "Wrote FASTA files for %s | DNA: %s -> %s | protein: %s -> %s",
            gpkg_name,
            dna_written,
            dna_path,
            protein_written,
            protein_path,
        )

    logger.info(
        "Finished writing per-type FASTA outputs | DNA total: %s | protein total: %s | types: %s",
        total_dna_written,
        total_protein_written,
        len(grouped_reads),
    )


def _calculate_effective_library_reads(meta: Dict[str, Any]) -> float:
    """
    Calculate effective library size for RPKM.

    If prokaryotic_fraction is given, multiply forward and reverse reads
    individually by that fraction before summing them.
    """
    try:
        forward_reads = float(meta.get("forward_reads") or 0.0)
    except Exception:
        forward_reads = 0.0

    try:
        reverse_reads = float(meta.get("reverse_reads") or 0.0)
    except Exception:
        reverse_reads = 0.0

    try:
        fraction = float(meta.get("prokaryotic_fraction")) if meta.get("prokaryotic_fraction") is not None else 1.0
    except Exception:
        fraction = 1.0

    return (forward_reads * fraction) + (reverse_reads * fraction)


def output_read_taxonomy_level_counts(
        directory: str,
        read_dict: Dict[Tuple[str, str, str], Read],
        metagenome_dict: Dict[str, Dict[str, Any]],
        lineage_dict: Dict[str, Dict[str, Any]],
        gpkg_length_dict: Dict[str, int | float],
        filename: str = "summary_read_taxonomy_level_counts.txt",
) -> None:
    """
    Summarise reads per metagenome, protein type, and taxonomic level,
    and calculate RPKM values.

    Aggregation level
    -----------------
    (metagenomeID, gpkg_name, taxonomic_level, taxon_name)

    Output columns
    --------------
    metagenomeID
    protein_name
    taxonomic_level
    read_count
    taxon_name
    forward_reads
    reverse_reads
    prokaryotic_fraction
    protein_length
    effective_library_reads
    rpkm
    """
    outpath = os.path.join(directory, filename)

    tax_levels = [
        ("kingdom", "Kingdom"),
        ("phylum", "Phylum"),
        ("class", "Class"),
        ("order", "Order"),
        ("family", "Family"),
        ("genus", "Genus"),
        ("species", "Species"),
    ]

    counts: Dict[Tuple[str, str, str, str], int] = defaultdict(int)

    # Count reads per metagenome, protein, taxonomic level, and taxon name
    for (_, gpkg_name, metagenomeID), read in read_dict.items():
        lineageID = _safe_str(getattr(read, "lineageID", None), default="NA")
        lineage_meta = lineage_dict.get(lineageID, {})

        for dict_key, label in tax_levels:
            taxon_name = _safe_str(lineage_meta.get(dict_key), default="")
            if not taxon_name or taxon_name == "NA":
                continue

            counts[(metagenomeID, gpkg_name, label, taxon_name)] += 1

    level_order = {label: i for i, (_, label) in enumerate(tax_levels)}

    sorted_keys = sorted(
        counts.keys(),
        key=lambda x: (
            _safe_str(x[0]).casefold(),  # metagenomeID
            _safe_str(x[1]).casefold(),  # protein_name / gpkg_name
            level_order.get(x[2], 999),  # taxonomic level order
            _safe_str(x[3]).casefold(),  # taxon name
        )
    )

    with open(outpath, "w") as out:
        out.write(
            "metagenomeID\tprotein_name\ttaxonomic_level\tread_count\ttaxon_name\tforward_reads\treverse_reads\tprokaryotic_fraction\tprotein_length\teffective_library_reads\trpkm\n"
        )

        for metagenomeID, gpkg_name, level, taxon_name in sorted_keys:
            meta = metagenome_dict.get(metagenomeID, {})

            forward_reads_raw = meta.get("forward_reads")
            reverse_reads_raw = meta.get("reverse_reads")
            prokaryotic_fraction_raw = meta.get("prokaryotic_fraction")

            forward_reads = _safe_str(forward_reads_raw)
            reverse_reads = _safe_str(reverse_reads_raw)
            prokaryotic_fraction = _safe_str(prokaryotic_fraction_raw)

            effective_library_reads = _calculate_effective_library_reads(meta)

            try:
                protein_length = float(gpkg_length_dict.get(gpkg_name, 0))
            except Exception:
                protein_length = 0.0

            read_count = counts[(metagenomeID, gpkg_name, level, taxon_name)]

            if effective_library_reads > 0 and protein_length > 0:
                rpkm = (float(read_count) * 1_000_000_000.0) / (
                        effective_library_reads * protein_length
                )
            else:
                rpkm = 0.0

            out.write(
                f"{metagenomeID}\t"
                f"{gpkg_name}\t"
                f"{level}\t"
                f"{read_count}\t"
                f"{taxon_name}\t"
                f"{forward_reads}\t"
                f"{reverse_reads}\t"
                f"{prokaryotic_fraction}\t"
                f"{protein_length}\t"
                f"{effective_library_reads}\t"
                f"{rpkm:.10f}\n"
            )

    logger.info(
        "Wrote metagenome/protein/taxonomy level summary table with RPKM: %s",
        outpath,
    )

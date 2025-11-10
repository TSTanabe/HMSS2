#!/usr/bin/env python3
from __future__ import annotations

import bisect
import re
from typing import Any, Dict, List, Tuple

from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def prepare_extended_cluster_ranges(
    cluster_dict: Dict[str, Any],
    flank_extension: int = 3000,
) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Prepare extended genomic ranges (±flank_extension) for each cluster
    and group them by contig.

    Args:
        cluster_dict: Dictionary with Cluster objects (clusterID -> Cluster)
                      Cluster must provide attributes:
                        - contig
                        - cluster_start
                        - cluster_end
        flank_extension: Number of base pairs to extend the cluster boundaries.

    Returns:
        Dict[contig, List[(start, end, clusterID)]]:
            For each contig, a sorted list of extended cluster intervals.
    """
    ranges: Dict[str, List[Tuple[int, int, str]]] = {}

    for cid, cluster in cluster_dict.items():
        contig = getattr(cluster, "contig", None)
        c_start = getattr(cluster, "cluster_start", None)
        c_end = getattr(cluster, "cluster_end", None)

        if contig and c_start is not None and c_end is not None:
            start = max(0, int(c_start) - flank_extension)
            end = int(c_end) + flank_extension
            ranges.setdefault(contig, []).append((start, end, cid))

    # Sort per contig for efficient interval search
    for contig in ranges:
        ranges[contig].sort(key=lambda x: x[0])

    return ranges


def add_context_genes_from_gff_fast(
    genome_id: str,
    gff_file: str,
    combined_protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    flank_extension: int = 3000,
) -> None:
    """
    Stream through a GFF file and add all CDS genes that fall inside
    any extended cluster range (cluster_start - flank_extension to
    cluster_end + flank_extension).

    Behavior:
      - GFF is read line-by-line (low memory footprint).
      - Only CDS features are considered.
      - Any CDS that lies fully inside one of the extended intervals is added
        as a non-valid Protein object (valid_hit = False) if not already present.

    Args:
        genome_id: Genome identifier.
        gff_file: Path to the GFF file.
        combined_protein_dict: Dict[proteinID -> Protein object]; will be extended.
        cluster_dict: Cluster objects for this genome.
        flank_extension: Extension in bp around each cluster.
    """
    # Build extended intervals from cluster_dict
    cluster_ranges = prepare_extended_cluster_ranges(
        cluster_dict,
        flank_extension=flank_extension,
    )
    if not cluster_ranges:
        logger.debug(f"[{genome_id}] No cluster intervals found, skipping context gene addition.")
        return

    # Precompile regex patterns
    id_pattern = re.compile(r"ID=(?:cds-)?(\S+?)(?:[;\s]|$)")
    locus_pattern = re.compile(r"locus_tag=(\S+?)(?:[;\s]|$)")

    # Precompute list of starts per contig for bisect
    range_index: Dict[str, List[int]] = {
        contig: [r[0] for r in ranges] for contig, ranges in cluster_ranges.items()
    }

    try:
        with open(gff_file, "r") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                ftype = parts[2].lower()
                if ftype != "cds":
                    continue

                contig = parts[0]
                if contig not in cluster_ranges:
                    continue

                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue

                attrs = parts[8]

                # Extract protein/gene ID
                m_id = id_pattern.search(attrs)
                if not m_id:
                    continue
                gene_id = m_id.group(1)

                # Skip if already present
                if gene_id in combined_protein_dict:
                    continue

                # Optional: locus_tag
                locus_tag = ""
                if "locus_tag" in attrs:
                    m_lt = locus_pattern.search(attrs)
                    if m_lt:
                        locus_tag = m_lt.group(1)

                ranges = cluster_ranges[contig]
                starts = range_index[contig]

                # Find rightmost interval whose start <= start
                idx = bisect.bisect_right(starts, start) - 1
                if idx < 0:
                    continue

                added = False

                # Scan leftwards from idx to check if any interval fully contains the CDS
                j = idx
                while j >= 0:
                    cstart, cend, cid = ranges[j]
                    if cstart <= start:
                        if end <= cend:
                            _add_nonvalid_protein_entry_fast(
                                genome_id=genome_id,
                                combined_protein_dict=combined_protein_dict,
                                gene_id=gene_id,
                                contig=contig,
                                start=start,
                                end=end,
                                strand=parts[6],
                                locus_tag=locus_tag,
                                comment=f"context_{cid}",
                            )
                            added = True
                            break
                        # If this interval ends before the CDS starts,
                        # earlier intervals cannot contain it.
                        if cend < start:
                            break
                    j -= 1

                # Edge case: check the next interval to the right (rare overlap layout)
                if not added and idx + 1 < len(ranges):
                    cstart, cend, cid = ranges[idx + 1]
                    if cstart <= start and end <= cend:
                        _add_nonvalid_protein_entry_fast(
                            genome_id=genome_id,
                            combined_protein_dict=combined_protein_dict,
                            gene_id=gene_id,
                            contig=contig,
                            start=start,
                            end=end,
                            strand=parts[6],
                            locus_tag=locus_tag,
                            comment=f"context_{cid}",
                        )

    except IOError as e:
        logger.warning(f"[{genome_id}] Error reading GFF '{gff_file}': {e}")


def _resolve_protein_class(combined_protein_dict: Dict[str, Any]):
    """
    Dynamically determine the Protein class:

      - Prefer the class of an existing Protein object in combined_protein_dict.
      - Fallback: try to import Protein from hmsss.parse_reports.parse_reports.

    Avoids hard cyclic imports at module level.
    """
    for obj in combined_protein_dict.values():
        return obj.__class__

    try:
        # Lazy import to avoid circular dependencies
        from hmsss.parse_reports.parse_reports import Protein  # type: ignore
        return Protein
    except Exception as e:
        logger.error(f"Could not resolve Protein class: {e}")
        return None


def _add_nonvalid_protein_entry_fast(
    genome_id: str,
    combined_protein_dict: Dict[str, Any],
    gene_id: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    locus_tag: str,
    comment: str,
) -> None:
    """
    Add a new non-valid Protein object if the gene_id does not yet exist.

    Behavior:
      - domain is set to "flank" (context gene)
      - score = 0.0
      - valid_hit = False (if attribute exists)
      - selection_comment is extended with the provided comment (if supported)
    """
    if gene_id in combined_protein_dict:
        return

    ProteinClass = _resolve_protein_class(combined_protein_dict)
    if ProteinClass is None:
        return

    try:
        # Expected typical signature:
        # Protein(proteinID, domain, domStart, domEnd, score, genomeID)
        p = ProteinClass(gene_id, "flank", 0, 1, 0.0, genome_id)
    except TypeError:
        # Fallback in case of a slightly different constructor
        p = ProteinClass(gene_id, "flank", 0, 1, 0.0)
        if hasattr(p, "genomeID"):
            setattr(p, "genomeID", genome_id)

    if hasattr(p, "gene_contig"):
        p.gene_contig = contig
    if hasattr(p, "gene_start"):
        p.gene_start = start
    if hasattr(p, "gene_end"):
        p.gene_end = end
    if hasattr(p, "gene_strand"):
        p.gene_strand = strand
    if hasattr(p, "gene_locustag"):
        p.gene_locustag = locus_tag
    if hasattr(p, "valid_hit"):
        p.valid_hit = False
    if hasattr(p, "add_selection_comment"):
        p.add_selection_comment(comment)

    combined_protein_dict[gene_id] = p

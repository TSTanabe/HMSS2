from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Dict, Tuple, Iterator


@dataclass
class Read:
    """
    Represents one mapped read with raw sequence + aligned sequence and taxonomy.

    Attributes
    ----------
    readID : str
        Identifier of the read (FASTA header without ">").
    alignment : str
        Aligned sequence (usually contains '-' gaps).
    sequence : str
        Raw / ungapped sequence.
    coverage : float
        Fraction of the read that is aligned (0..1).
    lineageID : str
        String identifier of the lineage (e.g., full lineage string).
    lineage : dict
        Parsed lineage mapping (rank -> taxon), if available.
    """

    readID: str
    alignment: str
    sequence: str
    gpkg_name: str
    metagenomeID: str = ""
    genomeID: str = ""
    base: str = ""
    start: int = 0
    end: int = 1
    coverage: float = 1.0
    lineageID: str = ""
    lineage: Dict[str, str] = field(default_factory=dict)

    def compute_coverage(self, *, hmm_length: int) -> float:
        """
        Compute coverage from current alignment + sequence.
        """

        def _ungapped_len(s: str) -> int:
            return sum(1 for c in s if c not in "-. \n\r\t")

        aligned_len = _ungapped_len(self.alignment or "")
        seq_len = _ungapped_len(self.sequence or "")

        if seq_len <= 0:
            self.coverage = 0.0
            return self.coverage

        denominator = seq_len
        if hmm_length is not None and 0 < hmm_length < seq_len:
            denominator = hmm_length

        cov = aligned_len / denominator
        if cov < 0.0:
            cov = 0.0
        elif cov > 1.0:
            cov = 1.0

        self.coverage = cov
        return self.coverage

    def aligned_span_in_sequence(self):
        """
        Locate the ungapped alignment sequence within the original read sequence.
        Determines only from where to where the alignment goes in the sequence

        The alignment is assumed to be against the HMM (not the read), so we:
          1) remove gaps from the alignment
          2) search the resulting sequence as a substring of `self.sequence`

        Returns
        -------
        (start, end) : tuple[int, int]
            0-based, end-exclusive coordinates in `self.sequence`,
            or None if not found.
        """
        aln = getattr(self, "alignment", None)
        seq = getattr(self, "sequence", None)

        if not aln or not seq:
            self.start = 0
            self.end = 1
            return None

        # build ungapped query from alignment
        query = "".join(c for c in self.alignment if c not in "-. \n\r\t")
        start = self.sequence.find(query[0:6])  # Workaround for less big seqs
        if start == -1:
            self.start = 0
        else:
            self.start = start

        self.end = self.start + len(query)


# ---------------------------------------------------------------------
# Subroutines for parsing GraftM output files
# ---------------------------------------------------------------------


def iter_fasta(fasta_path: str) -> Iterator[Tuple[str, str]]:
    """Yield (read_id, sequence) from a FASTA file (read_id = first token after '>')."""
    rid = None
    chunks: list[str] = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(chunks)
                rid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if rid is not None:
        yield rid, "".join(chunks)


def load_sequences_into_reads(
        reads: Dict[tuple, "Read"], fasta_path: str, gpkg_name: str
) -> None:
    """
    Stream a FASTA and write sequences into `reads` (keyed by readID).

    - If readID exists: set/overwrite `Read.sequence`
    - If readID does not exist: create `Read(readID=..., sequence=..., alignment="")`
    type is the gpkg package name, needed for keys
    """
    for rid, seq in iter_fasta(fasta_path):
        rid_key = (rid, gpkg_name)
        r = reads.get(rid_key)
        if r is None:
            reads[rid_key] = Read(
                readID=rid, sequence=seq, alignment="", gpkg_name=gpkg_name
            )
        else:
            r.sequence = seq


def load_alignments_into_reads(
        reads: Dict[tuple, "Read"], alignment_fasta: str, gpkg_name: str
) -> None:
    """
    Stream an alignment FASTA and write alignments into `reads` (keyed by readID).

    - If readID exists: set/overwrite `Read.alignment`
    - If readID does not exist: create `Read(readID=..., alignment=..., sequence="")`
    """
    for rid, aln in iter_fasta(alignment_fasta):
        rid_key = (rid, gpkg_name)
        r = reads.get(rid_key)
        if r is None:
            reads[rid_key] = Read(
                readID=rid, alignment=aln, sequence="", gpkg_name=gpkg_name
            )
        else:
            r.alignment = aln


def _lineage_to_id(lineage: str, n_hex: int = 12) -> str:
    canonical = ";".join(p.strip() for p in lineage.split(";") if p.strip())
    return hashlib.blake2s(
        canonical.encode("utf-8"), digest_size=n_hex // 2
    ).hexdigest()


def _parse_lineage_to_dict(lineage: str) -> Dict[str, str]:
    """
    Parse a lineage string into a dict mapping rank -> taxon.

    Expected input examples:
      "Root; k__Bacteria; p__Actinomycetota; c__Acidimicrobiia; o__UBA5794"
      "k__Bacteria;p__Firmicutes"
      "Root; Bacteria; Firmicutes"   (fallback mode)

    Returns
    -------
    dict[str, str]
        Example:
        {
          "root": "Root",
          "k": "Bacteria",
          "p": "Actinomycetota",
          "c": "Acidimicrobiia",
          "o": "UBA5794"
        }
    """
    lineage = lineage.strip()
    if not lineage:
        return {}

    parts = [p.strip() for p in lineage.split(";") if p.strip()]
    out: Dict[str, str] = {}

    for part in parts:
        # Case 1: rank-prefixed taxonomy (k__Bacteria)
        if "__" in part:
            rank, taxon = part.split("__", 1)
            rank = rank.strip()
            taxon = taxon.strip()
            if rank and taxon:
                out[rank] = taxon

        # Case 2: Root without prefix
        elif part.lower() == "root":
            out["root"] = "Root"

        # Case 3: fallback (no rank info)
        else:
            # store unranked entries in order
            key = f"lvl{len(out)}"
            out[key] = part

    return out


def load_read_taxonomy_into_reads(
        taxonomy_tsv: str, reads: Dict[tuple, "Read"], gpkg_name: str
) -> None:
    """
    Stream-read taxonomy file and write taxonomy fields directly into existing Read objects.

    Expected format per line:
      <readID><whitespace><lineage string>

    Notes:
    - Uses only the first token as readID; the rest of the line is lineage.
    - If readID is not present in `reads`, the line is ignored.
    """
    with open(taxonomy_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # rid = first token; lineage = rest of line
            rid, lineage = line.split(None, 1)
            lineage = lineage.strip()

            read = reads.get((rid, gpkg_name))  # Tuple key
            if read is None:
                continue

            read.lineageID = _lineage_to_id(lineage, n_hex=12)
            read.lineage = _parse_lineage_to_dict(lineage)


def finalize_reads(
        reads: Dict[tuple, Read],
        *,
        min_coverage: float,
        hmm_length: int,
) -> None:
    """
    Compute coverage for all reads and delete reads with coverage < min_coverage.

    Coverage rule:
      aligned_len = number of non-gap characters in Read.alignment (excludes '-' and '.')
      denom = len(Read.sequence) (ungapped) unless sequence is longer than hmm_length,
              then denom = hmm_length
      coverage = aligned_len / denom  (clamped to [0, 1])

    Mutates `reads` in-place (computes coverage + deletes low-coverage entries).
    """
    if hmm_length <= 0:
        hmm_length = 1
    if not (0.0 <= min_coverage <= 1.0):
        min_coverage = 0.1

    to_delete = []

    for key, r in reads.items():
        r.compute_coverage(hmm_length=hmm_length)
        r.aligned_span_in_sequence()
        if r.coverage < min_coverage:
            to_delete.append(key)

    for key in to_delete:
        del reads[key]


def build_reads_from_outputs(
        *,
        gpkg_name: str,
        taxonomy_csv: str,
        alignment_fasta: str,
        sequence_fasta: str,
        hmm_length: int,
        min_coverage: float,
) -> Dict[tuple, Read]:
    """
    Build Read objects keyed by readID by streaming the three output files.

    Order:
      1) sequences  -> create/fill Read.sequence
      2) alignments -> create/fill Read.alignment
      3) taxonomy   -> fill Read.lineageID + Read.lineage

    Returns
    -------
    dict[str, Read]
        Mapping readID -> Read
    """
    reads: Dict[tuple, Read] = {}

    load_sequences_into_reads(
        reads=reads, fasta_path=sequence_fasta, gpkg_name=gpkg_name
    )
    load_alignments_into_reads(
        reads=reads, alignment_fasta=alignment_fasta, gpkg_name=gpkg_name
    )

    # remove reads below coverage cutoff.
    finalize_reads(reads, min_coverage=min_coverage, hmm_length=hmm_length)

    load_read_taxonomy_into_reads(
        taxonomy_tsv=taxonomy_csv, reads=reads, gpkg_name=gpkg_name
    )

    return reads

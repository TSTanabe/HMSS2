import gzip
import os


def safe_read_count(path: str | None) -> int:
    """Return read count or 0 if path is None, missing, or empty."""
    if path is None:
        return 0
    if not isinstance(path, str):
        return 0
    if not os.path.isfile(path):
        return 0
    if os.path.getsize(path) == 0:
        return 0
    try:
        return count_reads_fasta_or_fastq(path)
    except Exception:
        return 0


def count_reads_fasta_or_fastq(path: str) -> int:
    """
    Count reads/sequences in a FASTA/FASTQ file.
    Supports gzipped (.gz) and uncompressed files.

    Rules:
      - FASTQ: count records as total_lines // 4
      - FASTA (and other fasta-like): count headers starting with '>'

    Format detection:
      - If first non-empty line starts with '@' -> FASTQ
      - If first non-empty line starts with '>' -> FASTA
      - Else: treat as FASTA-like and count '>' headers (raises if none found)
    """
    if not isinstance(path, str):
        raise TypeError("path must be a string")
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as f:
        # read first non-empty line to detect format
        first = ""
        while True:
            first = f.readline()
            if first == "":  # EOF
                return 0
            if first.strip():
                break

        f.seek(0)

        # FASTQ: count lines // 4 (valid FASTQ has 4 lines per record)
        if first.startswith("@"):
            n_lines = sum(1 for _ in f)
            return n_lines // 4

        # FASTA: count headers
        if first.startswith(">"):
            return sum(1 for line in f if line.startswith(">"))

        # Fallback: count '>' anyway; error if none found (unknown format)
        n_headers = sum(1 for line in f if line.startswith(">"))
        if n_headers == 0:
            raise ValueError(
                f"Unknown format (neither FASTQ '@' nor FASTA '>'): {path}"
            )
        return n_headers

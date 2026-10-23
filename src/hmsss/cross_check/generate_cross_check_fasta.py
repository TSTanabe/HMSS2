import os
from collections import OrderedDict
from hmsss.core.logging import get_logger

logger = get_logger(__name__)

class FastaHandleCache:
    """
    LRU cache for FASTA file handles to avoid 'too many open files'.
    """

    def __init__(self, out_dir: str, suffix: str, max_open: int = 32):
        self.out_dir = out_dir
        self.suffix = suffix
        self.max_open = max_open
        self.handles: OrderedDict[str, object] = OrderedDict()
        os.makedirs(out_dir, exist_ok=True)

    def get(self, domain_id: str):
        if domain_id in self.handles:
            fh = self.handles.pop(domain_id)
            self.handles[domain_id] = fh
            return fh

        # LRU eviction
        if len(self.handles) >= self.max_open:
            _, old_fh = self.handles.popitem(last=False)
            old_fh.close()

        path = os.path.join(self.out_dir, f"{domain_id}{self.suffix}")
        fh = open(path, "a")
        self.handles[domain_id] = fh
        return fh

    def close_all(self):
        for fh in self.handles.values():
            try:
                fh.close()
            except Exception:
                pass
        self.handles.clear()


def write_intermediate_hits_faa(
    protein_dict,
    out_dir: str,
    *,
    suffix: str = ".intermediate_hits.faa",
    max_open_files: int = 32,
    add_valid_hits: bool = True,
):
    """Write individual domain sequences to domain-specific FASTA files."""

    cache = FastaHandleCache(
        out_dir=out_dir, suffix=suffix, max_open=max_open_files
    )

    try:
        for protein in protein_dict.values():
            if not add_valid_hits and getattr(protein, "valid_hit", False):
                continue

            seq = protein.get_sequence()
            if not seq:
                continue

            protein_id = f"{protein.genomeID}-{protein.proteinID}"

            for dom in protein.get_domain_listing():
                start, end = int(dom.start), int(dom.end)

                # Coordinates are 1-based and inclusive
                if not (1 <= start <= end <= len(seq)):
                    logger.warning(
                        f"Invalid domain coordinates: {protein_id}, {start}-{end}"
                    )

                header = f"{protein_id}___dom_{start}_{end}"
                domain_seq = seq[start - 1:end]

                fh = cache.get(dom.get_domain())
                fh.write(f">{header}\n{domain_seq}\n")

    finally:
        cache.close_all()

import os
from collections import OrderedDict


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
):
    """
    Write non-valid hits (protein.valid_hit == False) as FASTA files,
    one file per domain/HMM.

    Header format:
        >genomeID-proteinID

    Each protein is written to *all* domain FASTAs for which it has hits.
    """
    cache = FastaHandleCache(
        out_dir=out_dir,
        suffix=suffix,
        max_open=max_open_files,
    )

    try:
        for protein in protein_dict.values():
            # Only intermediate (non-valid) hits
            if getattr(protein, "valid_hit", False):
                continue

            seq = protein.get_sequence()
            if not seq:
                continue

            header = f"{protein.genomeID}-{protein.proteinID}"

            for dom in protein.get_domain_listing():
                domain_id = dom.get_domain()
                fh = cache.get(domain_id)
                fh.write(f">{header}\n{seq}\n")

    finally:
        cache.close_all()

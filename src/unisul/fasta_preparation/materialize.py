from __future__ import annotations

import gzip
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


@contextmanager
def materialize_gz_next_to_input(path: str) -> Iterator[str]:
    """
    Yield a usable *plain* path.

    - If `path` is not gzipped: yield it unchanged.
    - If `path` ends with `.gz`: decompress into a temporary file in the SAME directory
      as the input file, yield that temp path, then delete it afterwards.

    Notes:
    - The temp file is hidden (prefix '.') and unique.
    - Cleanup happens even if an exception is raised in the caller.
    """
    p = Path(path)

    if not p.name.endswith(".gz"):
        yield str(p)
        return

    out_dir = p.parent
    plain_name = p.name[:-3]  # remove ".gz"

    fd, tmp_path = tempfile.mkstemp(
        prefix=f".{plain_name}.tmp_",
        dir=str(out_dir),
    )
    os.close(fd)

    try:
        with gzip.open(str(p), "rb") as fin, open(tmp_path, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)  # 1 MB chunks
        yield tmp_path
    finally:
        try:
            os.remove(tmp_path)
        except FileNotFoundError:
            pass


@contextmanager
def materialize_pair_gz_next_to_input(
    path_a: str, path_b: str
) -> Iterator[tuple[str, str]]:
    """
    Convenience wrapper to materialize two paths in a single 'with' statement.
    """
    with (
        materialize_gz_next_to_input(path_a) as a,
        materialize_gz_next_to_input(path_b) as b,
    ):
        yield a, b


@contextmanager
def materialize_single_next_to_input(path: str) -> Iterator[str]:
    """
    Materialize exactly ONE input file.

    - If `path` is not gzipped: yield unchanged path.
    - If `path` ends with `.gz`: decompress once into temp file next to input,
      yield temp path, then delete afterwards.

    This avoids using the pair-materializer with (path, path) which would
    decompress twice.
    """
    with materialize_gz_next_to_input(path) as p:
        yield p

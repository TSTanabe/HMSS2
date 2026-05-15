#!/usr/bin/env python3

import argparse
import os
import sqlite3
import tempfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

TAX_COLS = [
    "Superkingdom", "Phylum", "Class",
    "Ordnung", "Family", "Genus", "Species"
]


def connect_readonly(database: str) -> sqlite3.Connection:
    uri = f"file:{Path(database).resolve()}?mode=ro"
    con = sqlite3.connect(uri, uri=True, timeout=60)
    con.execute("PRAGMA query_only = TRUE;")
    con.execute("PRAGMA temp_store = MEMORY;")
    con.execute("PRAGMA cache_size = -262144;")
    con.execute("PRAGMA mmap_size = 2147483648;")
    return con


def fetch_all_domains(database: str) -> list[str]:
    print("[INFO] Fetching all valid domain types.")

    with connect_readonly(database) as con:
        df = pd.read_sql_query(
            """
            SELECT DISTINCT d.domain
            FROM Proteins p
            JOIN Domains d ON p.proteinID = d.proteinID
            WHERE p.valid_hit = 1
              AND d.domain IS NOT NULL
              AND d.domain != ''
            ORDER BY d.domain
            """,
            con,
        )

    domains = df["domain"].tolist()
    print(f"[INFO] Found {len(domains):,} domain types.")
    return domains


def fetch_genome_ids(database: str) -> list[str]:
    print("[INFO] Fetching genome IDs with valid hits.")

    with connect_readonly(database) as con:
        df = pd.read_sql_query(
            """
            SELECT DISTINCT genomeID
            FROM Proteins
            WHERE valid_hit = 1
            ORDER BY genomeID
            """,
            con,
        )

    genome_ids = df["genomeID"].tolist()
    print(f"[INFO] Found {len(genome_ids):,} genomes with valid hits.")
    return genome_ids


def chunk_list(values: list[str], chunk_size: int) -> list[list[str]]:
    return [
        values[i:i + chunk_size]
        for i in range(0, len(values), chunk_size)
    ]


def fetch_chunk_and_write(args) -> str:
    (
        database,
        genome_ids,
        domains,
        include_taxonomy,
        tmpdir,
        chunk_index,
    ) = args

    placeholders = ",".join(["?"] * len(genome_ids))

    hit_query = f"""
    SELECT DISTINCT
        p.genomeID,
        d.domain
    FROM Proteins p
    JOIN Domains d
        ON p.proteinID = d.proteinID
    WHERE p.valid_hit = 1
      AND d.domain IS NOT NULL
      AND d.domain != ''
      AND p.genomeID IN ({placeholders})
    """

    with connect_readonly(database) as con:
        hits = pd.read_sql_query(
            hit_query,
            con,
            params=genome_ids,
        )

        pam = pd.DataFrame(
            0,
            index=genome_ids,
            columns=domains,
            dtype="int8",
        )

        if not hits.empty:
            hits = hits.drop_duplicates()
            row_idx = pd.Categorical(
                hits["genomeID"],
                categories=genome_ids,
            ).codes
            col_idx = pd.Categorical(
                hits["domain"],
                categories=domains,
            ).codes

            valid = (row_idx >= 0) & (col_idx >= 0)
            pam.values[row_idx[valid], col_idx[valid]] = 1

        pam.insert(0, "genomeID", pam.index)
        pam = pam.reset_index(drop=True)

        if include_taxonomy:
            tax_query = f"""
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
            WHERE genomeID IN ({placeholders})
            """

            tax = pd.read_sql_query(
                tax_query,
                con,
                params=genome_ids,
            )

            pam = tax.merge(pam, on="genomeID", how="right")

    out = Path(tmpdir) / f"pam_chunk_{chunk_index:06d}.tsv"
    pam.to_csv(out, sep="\t", index=False)

    return str(out)


def concatenate_chunk_files(chunk_files: list[str], output: str) -> None:
    print("[INFO] Concatenating chunk files.")

    with open(output, "w") as out_handle:
        for i, file in enumerate(sorted(chunk_files)):
            with open(file, "r") as in_handle:
                for line_no, line in enumerate(in_handle):
                    if i > 0 and line_no == 0:
                        continue
                    out_handle.write(line)

            if (i + 1) % 10 == 0 or (i + 1) == len(chunk_files):
                print(f"[PROGRESS] Concatenated {i + 1}/{len(chunk_files)} chunks.")


def build_pam_parallel_chunked(
        database: str,
        output: str,
        threads: int,
        chunk_size: int,
        include_taxonomy: bool,
        keep_tmp: bool,
        tmpdir: str | None,
) -> None:
    domains = fetch_all_domains(database)
    genome_ids = fetch_genome_ids(database)

    genome_chunks = chunk_list(genome_ids, chunk_size)

    print(f"[INFO] Number of genome chunks: {len(genome_chunks):,}")
    print(f"[INFO] Threads: {threads}")
    print(f"[INFO] Chunk size: {chunk_size}")

    if tmpdir is None:
        tmp_context = tempfile.TemporaryDirectory(prefix="pam_chunks_")
        tmp_path = tmp_context.name
    else:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
        tmp_context = None
        tmp_path = tmpdir

    print(f"[INFO] Temporary directory: {tmp_path}")

    worker_args = [
        (
            database,
            chunk,
            domains,
            include_taxonomy,
            tmp_path,
            i,
        )
        for i, chunk in enumerate(genome_chunks)
    ]

    chunk_files = []

    print("[INFO] Starting parallel chunk processing.")

    try:
        with Pool(processes=threads) as pool:
            for i, chunk_file in enumerate(
                    pool.imap_unordered(fetch_chunk_and_write, worker_args),
                    1,
            ):
                chunk_files.append(chunk_file)

                print(
                    f"[PROGRESS] Finished chunk {i}/{len(worker_args)} "
                    f"({(i / len(worker_args)) * 100:.1f}%)"
                )

        concatenate_chunk_files(chunk_files, output)

        print(f"[DONE] Wrote PAM to: {output}")

    finally:
        if tmp_context is not None and not keep_tmp:
            tmp_context.cleanup()
        elif keep_tmp:
            print(f"[INFO] Temporary files kept in: {tmp_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Create a chunked parallel presence/absence matrix from SQLite."
    )

    parser.add_argument("--database", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--no-taxonomy", action="store_true")
    parser.add_argument("--tmpdir", default=None)
    parser.add_argument("--keep-tmp", action="store_true")

    args = parser.parse_args()

    build_pam_parallel_chunked(
        database=args.database,
        output=args.output,
        threads=args.threads,
        chunk_size=args.chunk_size,
        include_taxonomy=not args.no_taxonomy,
        keep_tmp=args.keep_tmp,
        tmpdir=args.tmpdir,
    )


if __name__ == "__main__":
    main()

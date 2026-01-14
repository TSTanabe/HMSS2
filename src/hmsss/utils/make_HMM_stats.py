#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sqlite3
import sys


def iter_threshold_rows(threshold_file: Path):
    """
    Yields (domain, threshold_type, cutoff) rows.

    Input format per line (whitespace or TAB separated):
        HMM_NAME  GATHERING  TRUSTED  NOISE

    The domain name is derived from HMM_NAME by removing the prefix up to the first underscore:
        grp0_CosH  ->  CosH
    If no underscore is present, the full HMM_NAME is used as domain name.
    """
    with threshold_file.open("r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()  # robust for mixed whitespace/tabs
            if len(parts) < 4:
                raise ValueError(f"<4 columns at {threshold_file}:{line_no}: {parts!r}")

            hmm_name = parts[0]
            gathering, trusted, noise = parts[1:4]

            domain = hmm_name.split("_", 1)[1] if "_" in hmm_name else hmm_name

            yield (domain, "gathering", float(gathering))
            yield (domain, "trusted", float(trusted))
            yield (domain, "noise", float(noise))


def main(db_path: str, threshold_path: str, out_tsv: str) -> None:
    db_p = Path(db_path)
    thr_p = Path(threshold_path)
    out_p = Path(out_tsv)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    if not db_p.exists():
        raise SystemExit(f"ERROR: DB not found: {db_p}")
    if not thr_p.exists():
        raise SystemExit(f"ERROR: Threshold file not found: {thr_p}")

    # IMPORTANT: open normally (not mode=ro), otherwise TEMP TABLE creation may fail.
    con = sqlite3.connect(str(db_p))
    cur = con.cursor()

    # Performance pragmas (safe: do not modify DB schema; TEMP objects only in-session)
    cur.execute("PRAGMA temp_store = MEMORY;")
    # cache_size: negative => KB. Adjust to your machine if needed.
    cur.execute("PRAGMA cache_size = -200000;")  # ~200 MB cache
    cur.execute("PRAGMA mmap_size = 268435456;")  # 256 MB mmap

    # TEMP threshold lookup table
    cur.execute("DROP TABLE IF EXISTS temp.thresholds;")
    cur.execute(
        """
        CREATE TEMP TABLE thresholds (
            domain TEXT NOT NULL,
            threshold_type TEXT NOT NULL,
            cutoff REAL NOT NULL
        );
        """
    )
    cur.execute("CREATE TEMP INDEX idx_thr_domain ON thresholds(domain);")

    thr_rows = list(iter_threshold_rows(thr_p))
    cur.executemany(
        "INSERT INTO thresholds(domain, threshold_type, cutoff) VALUES (?, ?, ?);",
        thr_rows,
    )

    # Single-pass aggregation + F1
    cur.execute(
        """
        WITH agg AS (
            SELECT
                t.domain AS domain,
                t.threshold_type AS threshold_type,
                t.cutoff AS threshold,

                SUM(CASE WHEN p.valid_hit = 1 AND d.score >= t.cutoff THEN 1 ELSE 0 END) AS TP,
                SUM(CASE WHEN p.valid_hit = 1 AND d.score <  t.cutoff THEN 1 ELSE 0 END) AS FN,
                SUM(CASE WHEN p.valid_hit = 0 AND d.score >= t.cutoff THEN 1 ELSE 0 END) AS FP,
                SUM(CASE WHEN p.valid_hit = 0 AND d.score <  t.cutoff THEN 1 ELSE 0 END) AS TN
            FROM thresholds t
            JOIN Domains  d ON d.domain = t.domain
            JOIN Proteins p ON p.proteinID = d.proteinID
            GROUP BY t.domain, t.threshold_type, t.cutoff
        )
        SELECT
            domain,
            threshold_type,
            threshold,
            TP, FP, FN, TN,
            CASE
                WHEN (2*TP + FP + FN) = 0 THEN 0.0
                ELSE (2.0*TP) / (2.0*TP + FP + FN)
            END AS F1
        FROM agg
        ORDER BY domain, threshold_type;
        """
    )

    header = ["domain", "threshold_type", "threshold", "TP", "FP", "FN", "TN", "F1"]

    with out_p.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(header) + "\n")
        for domain, thr_type, thr, tp, fp, fn, tn, f1 in cur.fetchall():
            out.write(
                "\t".join(
                    [
                        str(domain),
                        str(thr_type),
                        f"{float(thr):.6f}".rstrip("0").rstrip("."),
                        str(int(tp or 0)),
                        str(int(fp or 0)),
                        str(int(fn or 0)),
                        str(int(tn or 0)),
                        f"{float(f1):.6f}".rstrip("0").rstrip("."),
                    ]
                )
                + "\n"
            )

    con.close()
    print(f"[OK] wrote: {out_p}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: make_HMM_stats.py <db.sqlite> <Thresholds> <out.tsv>")
        raise SystemExit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])

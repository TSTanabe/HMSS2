import argparse
import sqlite3
from collections import defaultdict
from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, norm
from concurrent.futures import ProcessPoolExecutor, as_completed


# ---------------------------------------------------------------------
# Threshold parsing (mirrors make_threshold_dict logic)
# ---------------------------------------------------------------------
def load_thresholds(threshold_file: str):
    """
    Thresholds file has NO header.
    Format per line (tab-separated):
      col0: HMM name with prefix (e.g. v10_xxx_DHPS_Taurine_Isethionat)
      col1: optimized cutoff
      col2: trusted cutoff
      col3: noise cutoff

    Protein type is everything after the first underscore in col0.
    If multiple lines map to the same protein_type, keep the MIN cutoff (per cutoff type).
    """
    thresholds = {
        "optimized": {},
        "trusted": {},
        "noise": {},
    }

    with open(threshold_file, "r") as fh:
        for ln, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"Thresholds parse error: expected 4 tab-separated columns at line {ln}, "
                    f"got {len(parts)}: {line!r}"
                )

            hmm_name = parts[0].strip()
            if "_" not in hmm_name:
                # If you *want* to hard-fail instead, replace with raise ValueError(...)
                continue

            protein_type = hmm_name.split("_", 1)[1]

            def _maybe_float(s: str):
                s = s.strip()
                if not s or s.upper() == "NA":
                    return None
                return float(s)

            opt = _maybe_float(parts[1])
            tru = _maybe_float(parts[2])
            noi = _maybe_float(parts[3])

            for cutoff_type, cutoff in (
                ("optimized", opt),
                ("trusted", tru),
                ("noise", noi),
            ):
                if cutoff is None:
                    continue
                prev = thresholds[cutoff_type].get(protein_type)
                if prev is None or cutoff < prev:
                    thresholds[cutoff_type][protein_type] = cutoff

    return thresholds


# ---------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------
# --- Worker: läuft in separatem Prozess ---
def _confmat_worker(
    db_path: str,
    total_proteins: int,
    protein_type: str,
    cutoff_type: str,
    cutoff: float,
):
    # read-only Verbindung (URI) + query_only Schutz
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute("PRAGMA query_only=ON;")

    query = """
    WITH dom_best AS (
        SELECT
            d.proteinID,
            MAX(d.score) AS best_score,
            p.valid_hit AS valid_hit
        FROM Domains d
        JOIN Proteins p ON p.proteinID = d.proteinID
        WHERE d.domain = ?
        GROUP BY d.proteinID
    )
    SELECT
        SUM(CASE WHEN valid_hit = 1 AND best_score >= ? THEN 1 ELSE 0 END) AS TP,
        SUM(CASE WHEN valid_hit = 0 AND best_score >= ? THEN 1 ELSE 0 END) AS FP,
        SUM(CASE WHEN valid_hit = 1 AND best_score <  ? THEN 1 ELSE 0 END) AS FN,
        COUNT(*) AS assigned_hits
    FROM dom_best
    """

    row = cur.execute(query, (protein_type, cutoff, cutoff, cutoff)).fetchone()

    TP = row["TP"] or 0
    FP = row["FP"] or 0
    FN = row["FN"] or 0
    assigned_hits = row["assigned_hits"] or 0
    TN = total_proteins - assigned_hits

    conn.close()

    return (protein_type, cutoff_type, cutoff, TP, FP, FN, TN, assigned_hits)


def compute_confusion_matrices_parallel(
    db_path,
    thresholds,
    total_proteins,
    *,
    max_workers: int = 8,
    chunksize: int = 50,
):
    """
    Parallelisierte Variante.
    Yields tuples:
    (protein_type, cutoff_type, cutoff, TP, FP, FN, TN, assigned_hits)
    """

    # Aufgabenliste bauen (damit wir Fortschritt sauber im Main-Process loggen können)
    tasks = []
    for cutoff_type, cutoff_dict in thresholds.items():
        items = list(cutoff_dict.items())
        print(f"[INFO] processing cutoff type: {cutoff_type} ({len(items)} HMMs)")
        for i, (protein_type, cutoff) in enumerate(items, start=1):
            # Nur Main-Process printet
            print(
                f"[INFO]  {cutoff_type}: HMM {i}/{len(items)} → {protein_type} (cutoff={cutoff})"
            )
            tasks.append((protein_type, cutoff_type, float(cutoff)))

    # Prozesse starten
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        # Optional: batching reduziert Overhead bei vielen tausend Jobs
        # Wir submitten in Blöcken, um nicht zehntausende Futures auf einmal zu halten.
        for start in range(0, len(tasks), chunksize):
            block = tasks[start : start + chunksize]
            futures = [
                ex.submit(
                    _confmat_worker,
                    db_path,
                    total_proteins,
                    protein_type,
                    cutoff_type,
                    cutoff,
                )
                for (protein_type, cutoff_type, cutoff) in block
            ]
            for fut in as_completed(futures):
                yield fut.result()


def compute_confusion_matrices(db_path, thresholds, total_proteins):
    """
    Yields tuples:
    (protein_type, cutoff_type, cutoff, TP, FP, FN, TN, assigned_hits)
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    for cutoff_type, cutoff_dict in thresholds.items():
        print(f"[INFO] processing cutoff type: {cutoff_type} ({len(cutoff_dict)} HMMs)")
        for i, (protein_type, cutoff) in enumerate(cutoff_dict.items(), start=1):
            print(
                f"[INFO]  {cutoff_type}: HMM {i}/{len(cutoff_dict)} → {protein_type} "
                f"(cutoff={cutoff})"
            )

            # Aggregate best score per protein for this domain
            query = """
            WITH dom_best AS (
                SELECT
                    d.proteinID,
                    MAX(d.score) AS best_score,
                    p.valid_hit AS valid_hit
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                WHERE d.domain = ?
                GROUP BY d.proteinID
            )
            SELECT
                SUM(CASE WHEN valid_hit = 1 AND best_score >= ? THEN 1 ELSE 0 END) AS TP,
                SUM(CASE WHEN valid_hit = 0 AND best_score >= ? THEN 1 ELSE 0 END) AS FP,
                SUM(CASE WHEN valid_hit = 1 AND best_score <  ? THEN 1 ELSE 0 END) AS FN,
                COUNT(*) AS assigned_hits
            FROM dom_best
            """

            row = cur.execute(query, (protein_type, cutoff, cutoff, cutoff)).fetchone()

            TP = row["TP"] or 0
            FP = row["FP"] or 0
            FN = row["FN"] or 0
            assigned_hits = row["assigned_hits"] or 0

            TN = total_proteins - assigned_hits

            yield (
                protein_type,
                cutoff_type,
                cutoff,
                TP,
                FP,
                FN,
                TN,
                assigned_hits,
            )

    conn.close()


def _safe_div(num, den):
    return None if den == 0 else (num / den)


def balanced_accuracy(tp, fp, fn, tn):
    tpr = _safe_div(tp, tp + fn)  # sensitivity/recall
    tnr = _safe_div(tn, tn + fp)  # specificity
    if tpr is None or tnr is None:
        return None
    return 0.5 * (tpr + tnr)


def f1_score(tp, fp, fn, tn=None):
    den = 2 * tp + fp + fn
    return None if den == 0 else (2 * tp) / den


def mcc(tp, fp, fn, tn):
    a = tp + fp
    b = tp + fn
    c = tn + fp
    d = tn + fn
    den = a * b * c * d
    if den == 0:
        return None
    return (tp * tn - fp * fn) / math.sqrt(den)


def fmt(x):
    return "NA" if x is None else f"{x:.6f}"


#
#
#


def plot_f1_boxplots(tsv_path: str, out_prefix: str, alpha: float = 0.05):
    """
    Headerless TSV, fixed indices.
    - Boxplot + all points (black), jittered, no double outliers (showfliers=False)
    - exclude rows with TP=FP=TN=0
    - Kruskal-Wallis across trusted/optimized/noise
    - Dunn post-hoc (Holm corrected) for pairwise comparisons
    - Display stats bottom-left in plotting area with significance labels
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import kruskal, rankdata, norm

    # ---- fixed column indices (0-based) ----
    CUTOFF_TYPE_COL = 1
    TP_COL = 3
    FP_COL = 4
    TN_COL = 6
    F1_COL = 9

    # ---- helpers ----
    def sig_stars(p):
        if p is None or np.isnan(p):
            return "NA"
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    def dunn_test_holm(groups_dict):
        """
        Dunn pairwise test with Holm correction.
        Returns list of dicts: {g1,g2,z,p,p_adj,stars,better}
        better indicates which group has higher mean rank.
        """
        names = list(groups_dict.keys())
        values = np.concatenate([groups_dict[k] for k in names])
        labels = np.concatenate([[k] * len(groups_dict[k]) for k in names])

        N = len(values)
        ranks = rankdata(values)

        mean_ranks = {k: ranks[labels == k].mean() for k in names}
        ns = {k: len(groups_dict[k]) for k in names}

        raw = []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                g1, g2 = names[i], names[j]
                num = mean_ranks[g1] - mean_ranks[g2]
                den = np.sqrt((N * (N + 1) / 12.0) * (1.0 / ns[g1] + 1.0 / ns[g2]))
                z = num / den
                p = 2.0 * (1.0 - norm.cdf(abs(z)))
                better = g1 if z > 0 else g2
                raw.append({"g1": g1, "g2": g2, "z": z, "p": p, "better": better})

        # Holm correction
        ps = np.array([r["p"] for r in raw], dtype=float)
        order = np.argsort(ps)
        m = len(ps)
        p_adj = np.empty_like(ps)
        for k, idx in enumerate(order):
            p_adj[idx] = min(ps[idx] * (m - k), 1.0)

        for r, pa in zip(raw, p_adj):
            r["p_adj"] = float(pa)
            r["stars"] = sig_stars(r["p_adj"])

        return raw

    # ---- load headerless TSV ----
    df = pd.read_csv(tsv_path, sep="\t", header=None, dtype=str)
    need = max(CUTOFF_TYPE_COL, TP_COL, FP_COL, TN_COL, F1_COL) + 1
    if df.shape[1] < need:
        raise ValueError(f"TSV has {df.shape[1]} columns, need at least {need}.")

    cutoff_type = df.iloc[:, CUTOFF_TYPE_COL].astype(str).str.strip().str.lower()
    tp = pd.to_numeric(df.iloc[:, TP_COL], errors="coerce")
    fp = pd.to_numeric(df.iloc[:, FP_COL], errors="coerce")
    tn = pd.to_numeric(df.iloc[:, TN_COL], errors="coerce")
    f1 = pd.to_numeric(df.iloc[:, F1_COL], errors="coerce")

    # valid rows
    mask = cutoff_type.notna() & tp.notna() & fp.notna() & tn.notna() & f1.notna()
    cutoff_type = cutoff_type[mask]
    tp, fp, tn, f1 = tp[mask], fp[mask], tn[mask], f1[mask]

    # exclude degenerate rows: TP=FP=TN=0
    mask_nonempty = ~((tp == 0) & (fp == 0) & (tn == 0))
    cutoff_type = cutoff_type[mask_nonempty]
    f1 = f1[mask_nonempty]

    order = ["trusted", "optimized", "noise"]
    groups = {k: f1[cutoff_type == k].to_numpy(dtype=float) for k in order}

    present = sorted(cutoff_type.unique().tolist())
    for k in order:
        if groups[k].size == 0:
            raise ValueError(
                f"No usable F1 values for cutoff_type={k!r} after filtering.\n"
                f"Present cutoff types: {present}"
            )

    # ---- global test ----
    H, p_kw = kruskal(groups["trusted"], groups["optimized"], groups["noise"])

    # ---- post-hoc ----
    dunn = dunn_test_holm(groups)

    # summarize significant pairs
    sig_pairs = [r for r in dunn if r["p_adj"] < alpha]
    if len(sig_pairs) == 0:
        sig_line = "Dunn (Holm): none significant"
    else:
        # show direction + stars
        parts = []
        for r in sig_pairs:
            g1, g2 = r["g1"], r["g2"]
            better = r["better"]
            worse = g2 if better == g1 else g1
            parts.append(f"{better}>{worse} {r['stars']} (p={r['p_adj']:.2g})")
        sig_line = "Dunn (Holm): " + "; ".join(parts)

    # also include all pair results compactly (including ns)
    # trusted vs optimized, trusted vs noise, optimized vs noise
    def _pair_line(g1, g2):
        for r in dunn:
            if (r["g1"] == g1 and r["g2"] == g2) or (r["g1"] == g2 and r["g2"] == g1):
                better = r["better"]
                worse = g2 if better == g1 else g1
                return (
                    f"{g1} vs {g2}: {better}>{worse} {r['stars']} (p={r['p_adj']:.2g})"
                )
        return f"{g1} vs {g2}: NA"

    pair_lines = [
        _pair_line("trusted", "optimized"),
        _pair_line("trusted", "noise"),
        _pair_line("optimized", "noise"),
    ]

    n_tr, n_op, n_no = (
        groups["trusted"].size,
        groups["optimized"].size,
        groups["noise"].size,
    )

    stats_text = (
        f"Kruskal–Wallis: H={H:.3g}, p={p_kw:.2g}\n"
        f"n(trusted,opt,noise)=({n_tr},{n_op},{n_no})\n" + "\n".join(pair_lines)
    )

    # ---- plot ----
    plt.figure()

    data = [groups[k] for k in order]
    try:
        plt.boxplot(data, tick_labels=order, showfliers=False)
    except TypeError:
        # older matplotlib
        plt.boxplot(data, labels=order, showfliers=False)

    # overlay all points in black
    rng = np.random.default_rng(42)
    for i, k in enumerate(order, start=1):
        y = groups[k]
        x = rng.normal(loc=i, scale=0.06, size=y.size)
        plt.scatter(x, y, s=12, alpha=0.7, c="black")

    plt.ylabel("F1 score")
    plt.title("F1 score by cutoff type (boxplot + all points)")
    ax = plt.gca()

    # place stats bottom-left inside axes
    ax.text(
        0.02,
        0.02,
        stats_text,
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
    )

    plt.tight_layout()

    png_path = f"{out_prefix}.f1_boxplot.png"
    pdf_path = f"{out_prefix}.f1_boxplot.pdf"
    plt.savefig(png_path, dpi=200)
    plt.savefig(pdf_path)
    plt.close()

    # Keep return signature stable (2 values), so main doesn't break
    return png_path, pdf_path


def dunn_test(groups: dict, method: str = "holm"):
    """
    Pairwise Dunn test with multiple-testing correction.
    groups: dict {name: np.array}
    Returns list of tuples:
      (group1, group2, z, p_uncorrected, p_corrected)
    """

    names = list(groups.keys())
    values = np.concatenate([groups[k] for k in names])
    labels = np.concatenate([[k] * len(groups[k]) for k in names])

    ranks = rankdata(values)
    N = len(values)

    # Mean rank per group
    mean_ranks = {k: ranks[labels == k].mean() for k in names}
    ns = {k: len(groups[k]) for k in names}

    results = []

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            g1, g2 = names[i], names[j]

            num = mean_ranks[g1] - mean_ranks[g2]
            den = np.sqrt((N * (N + 1) / 12.0) * (1.0 / ns[g1] + 1.0 / ns[g2]))

            z = num / den
            p = 2 * (1 - norm.cdf(abs(z)))  # two-sided
            results.append([g1, g2, z, p])

    # ---- multiple testing correction ----
    ps = np.array([r[3] for r in results])

    if method == "holm":
        order = np.argsort(ps)
        adj = np.empty_like(ps)
        m = len(ps)
        for i, idx in enumerate(order):
            adj[idx] = min((m - i) * ps[idx], 1.0)
    elif method == "bh":
        order = np.argsort(ps)
        adj = np.empty_like(ps)
        m = len(ps)
        for i, idx in enumerate(order):
            adj[idx] = ps[idx] * m / (i + 1)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
    else:
        raise ValueError("method must be 'holm' or 'bh'")

    for r, p_adj in zip(results, adj):
        r.append(p_adj)

    return results


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compute confusion matrices from SQLite DB"
    )
    ap.add_argument("--db", required=True, help="SQLite database file")
    ap.add_argument("--thresholds", required=True, help="Thresholds file")
    ap.add_argument("--out", required=True, help="Output TSV file")
    ap.add_argument(
        "--total-proteins",
        type=int,
        default=430842000,
        help="Total number of proteins (for TN calculation)",
    )
    ap.add_argument(
        "--plot",
        action="store_true",
        help="After writing TSV, create boxplots from metrics.",
    )
    ap.add_argument(
        "--plot-prefix",
        default=None,
        help="Output prefix for plots (default: same as --out without extension).",
    )

    args = ap.parse_args()

    thresholds = load_thresholds(args.thresholds)

    with open(args.out, "w") as out:
        for row in compute_confusion_matrices_parallel(
            args.db, thresholds, args.total_proteins
        ):
            (protein_type, cutoff_type, cutoff, TP, FP, FN, TN, assigned_hits) = row

            ba = balanced_accuracy(TP, FP, FN, TN)
            f1 = f1_score(TP, FP, FN)
            mm = mcc(TP, FP, FN, TN)

            line = "\t".join(
                [
                    protein_type,
                    cutoff_type,
                    str(cutoff),
                    str(TP),
                    str(FP),
                    str(FN),
                    str(TN),
                    str(assigned_hits),
                    fmt(ba),
                    fmt(f1),
                    fmt(mm),
                ]
            )

            out.write(line + "\n")
            if getattr(args, "verbose", 0) >= 1:
                print(line)

    if args.plot:
        if args.plot_prefix is None:
            # default prefix: strip .tsv if present
            p = Path(args.out)
            prefix = str(p.with_suffix("")) if p.suffix else str(p)
        else:
            prefix = args.plot_prefix

        png_path, pdf_path = plot_f1_boxplots(args.out, prefix)
        print(f"[plot] wrote: {png_path}")
        print(f"[plot] wrote: {pdf_path}")


if __name__ == "__main__":
    main()

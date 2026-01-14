#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import math
import matplotlib.pyplot as plt

from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path

DEFAULT_FILENAME = "evaluation_stats.txt"
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


# ------------------------------------------------------------
# Parsing helpers
# ------------------------------------------------------------
def _coerce_value(val: str) -> Any:
    val = val.strip()
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return val


def parse_eval_stats(path: Path) -> Dict[str, Any]:
    """
    Parse a key=value text file (evaluation_stats.txt) into a dictionary.

    - Ignores empty lines and comments starting with '#'
    - Coerces values to int, then float, else keeps as string
    """
    d: Dict[str, Any] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            d[k.strip()] = _coerce_value(v)
    return d


def find_eval_files(root: Path, filename: str) -> List[Path]:
    """
    Recursively search for files named `filename` under `root`.
    Returns a sorted list of paths.
    """
    hits: List[Path] = []
    for dirpath, _, files in os.walk(root):
        if filename in files:
            hits.append(Path(dirpath) / filename)
    return sorted(hits)


# ------------------------------------------------------------
# Table formatting
# ------------------------------------------------------------
def _fmt(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        if x != x:  # NaN
            return "nan"
        return f"{x:.6f}".rstrip("0").rstrip(".")
    return str(x)


def print_table(rows: List[Dict[str, Any]], columns: List[str]) -> None:
    """
    Print a simple fixed-width table to stdout.
    """
    widths = {c: len(c) for c in columns}
    for r in rows:
        for c in columns:
            widths[c] = max(widths[c], len(_fmt(r.get(c))))

    header = "  ".join(c.ljust(widths[c]) for c in columns)
    sep = "  ".join("-" * widths[c] for c in columns)
    print(header)
    print(sep)

    for r in rows:
        print("  ".join(_fmt(r.get(c)).ljust(widths[c]) for c in columns))


# ------------------------------------------------------------
# Core API (callable as a normal function)
# ------------------------------------------------------------
def summarize_eval_stats(
        root: Path,
        *,
        filename: str = DEFAULT_FILENAME,
        sort: str = "overall_MCC",
        desc: bool = False,
        show_path: bool = False,
) -> int | list[dict[str, Any]]:
    """
    Summarize evaluation_stats files under `root` as a table (one row per package).

    Parameters
    ----------
    root:
        Root directory to search recursively.
    filename:
        Stats filename to look for (default: evaluation_stats.txt).
    sort:
        Column name to sort by (default: overall_MCC).
    desc:
        Sort descending if True (default: False).
    show_path:
        Include a relative path column if True (default: False).

    Returns
    -------
    int
        Exit-like status code:
          0 on success,
          1 if no files found,
          2 for invalid root directory.
    """
    root = root.resolve()
    if not root.is_dir():
        print(f"ERROR: not a directory: {root}", file=sys.stderr)
        return 2

    files = find_eval_files(root, filename)
    if not files:
        print(f"No {filename} found under {root}", file=sys.stderr)
        return 1

    rows: List[Dict[str, Any]] = []

    for fp in files:
        d = parse_eval_stats(fp)
        row: Dict[str, Any] = {}

        row["package"] = fp.parent.name
        row["protein_type"] = d.get("protein_type", "")

        if show_path:
            row["path"] = str(fp.parent.relative_to(root))

        # --- Overall confusion ---
        for k in (
                "overall_TP",
                "overall_FP",
                "overall_FN",
                "overall_TN",
                "overall_balanced_accuracy",
                "overall_F1",
                "overall_MCC",
        ):
            row[k] = d.get(k)

        # --- Placement confusion ---
        for k in (
                "placement_TP",
                "placement_FP",
                "placement_FN",
                "placement_TN",
                "placement_balanced_accuracy",
                "placement_F1",
                "placement_MCC",
        ):
            row[k] = d.get(k)

        # --- Fragment & key counts ---
        row["TP_frags"] = d.get("TP_frags")
        row["TN_frags"] = d.get("TN_frags")
        row["n_test_keys"] = d.get("n_test_keys")
        row["fp_test_keys"] = d.get("fp_test_keys")
        row["tp_lowest_rank"] = d.get("tp_lowest_rank")

        # --- Per-rank stats ---
        for r in RANKS:
            row[f"assigned_{r}"] = d.get(f"assigned_{r}")
            row[f"correct_{r}"] = d.get(f"correct_{r}")
            row[f"wrong_{r}"] = d.get(f"wrong_{r}")
            row[f"not_assigned_{r}"] = d.get(f"not_assigned_{r}")

        rows.append(row)

    try:
        # Sorting
        # Note: if the sort column mixes types (str vs float), Python may raise on comparison.
        rows.sort(
            key=lambda r: (r.get(sort) is None, r.get(sort)),
            reverse=desc,
        )
        return rows
    except Exception as e:
        print("")
        return rows


def print_eval_stats(
        rows: list[dict],
        *,
        root: Path | None = None,
) -> int:
    """
    Print evaluation statistics for multiple packages.

    Parameters
    ----------
    rows
        List of dictionaries, one per package (output of summarize_eval_stats).
    root
        Optional root directory that was searched (for summary line).
    """
    if not rows:
        print("No evaluation statistics to display.")
        return 0

    # Column order
    cols = ["package", "protein_type"]

    cols += [
        "overall_TP",
        "overall_FP",
        "overall_FN",
        "overall_TN",
        "overall_balanced_accuracy",
        "overall_F1",
        "overall_MCC",
        "placement_TP",
        "placement_FP",
        "placement_FN",
        "placement_TN",
        "placement_balanced_accuracy",
        "placement_F1",
        "placement_MCC",
        "TP_frags",
        "TN_frags",
        "n_test_keys",
        "fp_test_keys",
        "tp_lowest_rank",
    ]

    for r in RANKS:
        cols += [
            f"assigned_{r}",
            f"correct_{r}",
            f"wrong_{r}",
            f"not_assigned_{r}",
        ]

    print_table(rows, cols)

    if root is not None:
        print(f"\nFound {len(rows)} package(s) under {root}")
    else:
        print(f"\nFound {len(rows)} package(s)")

    return 0


# ------------------------------------------------------------
# Print plots routines
# ------------------------------------------------------------

def plot_overall_metrics_boxplots(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
) -> None:
    """
    Create boxplots for overall balanced accuracy, F1, and MCC
    from an eval stat list (one dict per package).

    Parameters
    ----------
    eval_stat_list
        List of dictionaries with evaluation statistics (one per package).
    out_dir
        Directory where PNG files will be written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = {
        "overall_balanced_accuracy": "Overall balanced accuracy",
        "overall_F1": "Overall F1 score",
        "overall_MCC": "Overall MCC",
    }

    def _collect(metric: str) -> List[float]:
        vals: List[float] = []
        for row in eval_stat_list:
            v = row.get(metric)
            if isinstance(v, (int, float)):
                if isinstance(v, float) and math.isnan(v):
                    continue
                vals.append(float(v))
        return vals

    for key, title in metrics.items():
        values = _collect(key)

        if not values:
            print(f"[WARN] No valid values for {key}, skipping plot.")
            continue

        plt.figure()
        plt.boxplot(values)
        plt.ylabel(title)
        plt.title(title)

        plt.ylim(0.0, 1.0)
        out_path = out_dir / f"{key}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"[OK] Wrote boxplot: {out_path}")


def plot_placement_metrics_boxplots(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
) -> None:
    """
    Create boxplots for placement balanced accuracy, F1, and MCC
    from an eval stat list (one dict per package).

    Parameters
    ----------
    eval_stat_list
        List of dictionaries with evaluation statistics (one per package).
    out_dir
        Directory where PNG files will be written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics = {
        "placement_balanced_accuracy": "Placement balanced accuracy",
        "placement_F1": "Placement F1 score",
        "placement_MCC": "Placement MCC",
    }

    def _collect(metric: str) -> List[float]:
        vals: List[float] = []
        for row in eval_stat_list:
            v = row.get(metric)
            if isinstance(v, (int, float)):
                if isinstance(v, float) and math.isnan(v):
                    continue
                vals.append(float(v))
        return vals

    for key, title in metrics.items():
        values = _collect(key)

        if not values:
            print(f"[WARN] No valid values for {key}, skipping plot.")
            continue

        plt.figure()
        plt.boxplot(values)
        plt.ylabel(title)
        plt.title(title)

        out_path = out_dir / f"{key}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"[OK] Wrote boxplot: {out_path}")


def plot_rank_assignment_fractions_boxplot(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        ranks: List[str] = RANKS,
        filename: str = "rank_assignment_fractions_boxplot.png",
        title: str = "Assigned fraction per rank (assigned / TP_frags)",
) -> None:
    """
    For each taxonomic rank, compute the fraction of assigned reads:
        assigned_<rank> / TP_frags
    across all packages, and show one boxplot per rank.

    Notes / assumptions
    -------------------
    - TP_frags is treated as the total number of possible hits (denominator).
    - assigned_<rank> is treated as the number of assigned hits at that rank (numerator).
    - Rows with missing assigned_<rank> or invalid TP_frags are skipped for that rank.

    Parameters
    ----------
    eval_stat_list
        List of dicts, one per package, containing TP_frags and assigned_<rank>.
    out_dir
        Output directory for the plot.
    ranks
        Rank order to show on the x-axis.
    filename
        Output image filename (PNG).
    title
        Plot title.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    # Collect per-rank fractions across all packages
    data_per_rank: List[List[float]] = []
    n_used_per_rank: List[int] = []

    for r in ranks:
        key = f"assigned_{r}"
        vals: List[float] = []

        for row in eval_stat_list:
            tp = row.get("TP_frags")
            assg = row.get(key)

            if not isinstance(tp, (int, float)) or _is_nan(float(tp)) or float(tp) <= 0:
                continue
            if not isinstance(assg, (int, float)) or _is_nan(float(assg)):
                continue

            frac = float(assg) / float(tp)

            # Defensive: fractions should be in [0, 1], but allow small numerical drift.
            if _is_nan(frac) or frac < 0:
                continue
            vals.append(frac)

        data_per_rank.append(vals)
        n_used_per_rank.append(len(vals))

    if not any(len(v) > 0 for v in data_per_rank):
        print("[WARN] No valid rank assignment fractions found. Skipping plot.")
        return

    # Matplotlib boxplot: one box per rank (empty ranks are allowed but will error if passed)
    # Therefore we filter to ranks that have at least one value, and keep labels aligned.
    filtered_data: List[List[float]] = []
    filtered_labels: List[str] = []
    filtered_counts: List[int] = []

    for r, vals, n in zip(ranks, data_per_rank, n_used_per_rank):
        if len(vals) == 0:
            continue
        filtered_data.append(vals)
        filtered_labels.append(r)
        filtered_counts.append(n)

    if not filtered_data:
        print("[WARN] All ranks empty after filtering. Skipping plot.")
        return

    plt.figure()
    plt.boxplot(filtered_data, labels=filtered_labels)
    plt.ylabel("Assigned fraction (assigned / TP_frags)")
    plt.title(title)

    # Optional: annotate how many packages contributed per rank (no colors; plain text)
    # Place just above x tick labels.
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    y_text = ymin + 0.02 * (ymax - ymin)
    for i, n in enumerate(filtered_counts, start=1):
        ax.text(i, y_text, f"n={n}", ha="center", va="bottom", fontsize=8)

    out_path = out_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote rank assignment boxplot: {out_path}")


def plot_rank_correct_fractions_boxplot(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        ranks: List[str] = RANKS,
        filename: str = "rank_correct_fractions_boxplot.png",
        title: str = "Correct fraction per rank (correct / TP_frags)",
) -> None:
    """
    For each taxonomic rank, compute the fraction of correct reads:
        correct_<rank> / TP_frags
    across all packages, and show one boxplot per rank.

    Notes / assumptions
    -------------------
    - TP_frags is treated as the total number of possible hits (denominator).
    - correct_<rank> is treated as the number of correct hits at that rank (numerator).
    - Rows with missing correct_<rank> or invalid TP_frags are skipped for that rank.

    Parameters
    ----------
    eval_stat_list
        List of dicts, one per package, containing TP_frags and correct_<rank>.
    out_dir
        Output directory for the plot.
    ranks
        Rank order to show on the x-axis.
    filename
        Output image filename (PNG).
    title
        Plot title.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    data_per_rank: List[List[float]] = []
    n_used_per_rank: List[int] = []

    for r in ranks:
        key = f"correct_{r}"
        vals: List[float] = []

        for row in eval_stat_list:
            tp = row.get("TP_frags")
            corr = row.get(key)

            if not isinstance(tp, (int, float)) or _is_nan(float(tp)) or float(tp) <= 0:
                continue
            if not isinstance(corr, (int, float)) or _is_nan(float(corr)):
                continue

            frac = float(corr) / float(tp)

            if _is_nan(frac) or frac < 0:
                continue
            vals.append(frac)

        data_per_rank.append(vals)
        n_used_per_rank.append(len(vals))

    if not any(len(v) > 0 for v in data_per_rank):
        print("[WARN] No valid rank correct fractions found. Skipping plot.")
        return

    filtered_data: List[List[float]] = []
    filtered_labels: List[str] = []
    filtered_counts: List[int] = []

    for r, vals, n in zip(ranks, data_per_rank, n_used_per_rank):
        if len(vals) == 0:
            continue
        filtered_data.append(vals)
        filtered_labels.append(r)
        filtered_counts.append(n)

    if not filtered_data:
        print("[WARN] All ranks empty after filtering. Skipping plot.")
        return

    plt.figure()
    plt.boxplot(filtered_data, labels=filtered_labels)
    plt.ylabel("Correct fraction (correct / TP_frags)")
    plt.title(title)

    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    y_text = ymin + 0.02 * (ymax - ymin)
    for i, n in enumerate(filtered_counts, start=1):
        ax.text(i, y_text, f"n={n}", ha="center", va="bottom", fontsize=8)

    out_path = out_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote rank correct boxplot: {out_path}")


def plot_rank_correct_given_assigned_boxplot(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        ranks: List[str] = RANKS,
        filename: str = "rank_correct_given_assigned_boxplot.png",
        title: str = "Correct among assigned per rank (correct / assigned)",
) -> None:
    """
    For each taxonomic rank, compute:
        correct_<rank> / assigned_<rank>
    across all packages, and show one boxplot per rank.

    Interpretation
    --------------
    This estimates the precision at each rank *conditional on being assigned*:
    "Of the reads that received an assignment at this rank, how many are correct?"

    Parameters
    ----------
    eval_stat_list
        List of dicts, one per package, containing assigned_<rank> and correct_<rank>.
    out_dir
        Output directory for the plot.
    ranks
        Rank order to show on the x-axis.
    filename
        Output image filename (PNG).
    title
        Plot title.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    data_per_rank: List[List[float]] = []
    n_used_per_rank: List[int] = []

    for r in ranks:
        key_assigned = f"assigned_{r}"
        key_correct = f"correct_{r}"
        vals: List[float] = []

        for row in eval_stat_list:
            assg = row.get(key_assigned)
            corr = row.get(key_correct)

            if not isinstance(assg, (int, float)) or _is_nan(float(assg)) or float(assg) <= 0:
                continue
            if not isinstance(corr, (int, float)) or _is_nan(float(corr)):
                continue

            frac = float(corr) / float(assg)

            # Defensive checks
            if _is_nan(frac) or frac < 0:
                continue
            vals.append(frac)

        data_per_rank.append(vals)
        n_used_per_rank.append(len(vals))

    if not any(len(v) > 0 for v in data_per_rank):
        print("[WARN] No valid correct/assigned fractions found. Skipping plot.")
        return

    # Filter empty ranks to avoid matplotlib errors
    filtered_data: List[List[float]] = []
    filtered_labels: List[str] = []
    filtered_counts: List[int] = []

    for r, vals, n in zip(ranks, data_per_rank, n_used_per_rank):
        if len(vals) == 0:
            continue
        filtered_data.append(vals)
        filtered_labels.append(r)
        filtered_counts.append(n)

    if not filtered_data:
        print("[WARN] All ranks empty after filtering. Skipping plot.")
        return

    plt.figure()
    plt.boxplot(filtered_data, labels=filtered_labels)
    plt.ylabel("Correct among assigned (correct / assigned)")
    plt.title(title)

    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    y_text = ymin + 0.02 * (ymax - ymin)
    for i, n in enumerate(filtered_counts, start=1):
        ax.text(i, y_text, f"n={n}", ha="center", va="bottom", fontsize=8)

    out_path = out_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote rank correct-given-assigned boxplot: {out_path}")


def plot_rank_wrong_given_assigned_boxplot(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        ranks: List[str] = RANKS,
        filename: str = "rank_wrong_given_assigned_boxplot.png",
        title: str = "Wrong among assigned per rank (wrong / assigned)",
) -> None:
    """
    For each taxonomic rank, compute:
        wrong_<rank> / assigned_<rank>
    across all packages, and show one boxplot per rank.

    Interpretation
    --------------
    This estimates the false-assignment rate at each rank *conditional on being assigned*:
    "Of the reads that received an assignment at this rank, how many are wrong?"

    Parameters
    ----------
    eval_stat_list
        List of dicts, one per package, containing assigned_<rank> and wrong_<rank>.
    out_dir
        Output directory for the plot.
    ranks
        Rank order to show on the x-axis.
    filename
        Output image filename (PNG).
    title
        Plot title.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    data_per_rank: List[List[float]] = []
    n_used_per_rank: List[int] = []

    for r in ranks:
        key_assigned = f"assigned_{r}"
        key_wrong = f"wrong_{r}"
        vals: List[float] = []

        for row in eval_stat_list:
            assg = row.get(key_assigned)
            wrong = row.get(key_wrong)

            if not isinstance(assg, (int, float)) or _is_nan(float(assg)) or float(assg) <= 0:
                continue
            if not isinstance(wrong, (int, float)) or _is_nan(float(wrong)):
                continue

            frac = float(wrong) / float(assg)

            # Defensive checks
            if _is_nan(frac) or frac < 0:
                continue
            vals.append(frac)

        data_per_rank.append(vals)
        n_used_per_rank.append(len(vals))

    if not any(len(v) > 0 for v in data_per_rank):
        print("[WARN] No valid wrong/assigned fractions found. Skipping plot.")
        return

    # Filter empty ranks (matplotlib cannot boxplot empty lists)
    filtered_data: List[List[float]] = []
    filtered_labels: List[str] = []
    filtered_counts: List[int] = []

    for r, vals, n in zip(ranks, data_per_rank, n_used_per_rank):
        if len(vals) == 0:
            continue
        filtered_data.append(vals)
        filtered_labels.append(r)
        filtered_counts.append(n)

    if not filtered_data:
        print("[WARN] All ranks empty after filtering. Skipping plot.")
        return

    plt.figure()
    plt.boxplot(filtered_data, labels=filtered_labels)
    plt.ylabel("Wrong among assigned (wrong / assigned)")
    plt.title(title)

    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    y_text = ymin + 0.02 * (ymax - ymin)
    for i, n in enumerate(filtered_counts, start=1):
        ax.text(i, y_text, f"n={n}", ha="center", va="bottom", fontsize=8)

    out_path = out_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote rank wrong-given-assigned boxplot: {out_path}")


#
#
#

def write_metric_values_tsv(
        eval_stat_list: List[Dict[str, Any]],
        *,
        metric_key: str,
        out_tsv: Path,
        include_package_cols: bool = True,
) -> Path:
    """
    Write the exact values used for a single-metric boxplot (e.g. overall_MCC)
    to a TSV file.

    Output columns
    --------------
    metric, value, [package, protein_type]
    """
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    lines: List[str] = []
    header = ["metric", "value"]
    if include_package_cols:
        header += ["package", "protein_type"]
    lines.append("\t".join(header))

    n_written = 0
    for row in eval_stat_list:
        v = row.get(metric_key)
        if not isinstance(v, (int, float)) or _is_nan(float(v)):
            continue

        rec = [metric_key, str(float(v))]
        if include_package_cols:
            rec += [str(row.get("package", "")), str(row.get("protein_type", ""))]
        lines.append("\t".join(rec))
        n_written += 1

    out_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[OK] Wrote TSV for {metric_key}: {out_tsv} (n={n_written})")
    return out_tsv


def write_rank_fraction_values_tsv(
        eval_stat_list: List[Dict[str, Any]],
        *,
        numerator_prefix: str,  # e.g. "assigned" / "correct" / "wrong"
        denominator: str,  # either "TP_frags" OR another prefix like "assigned"
        out_tsv: Path,
        ranks: List[str] = RANKS,
) -> Path:
    """
    Write the exact per-rank fractions used for rank boxplots to TSV.

    Fractions per rank:
      numerator = f"{numerator_prefix}_{rank}"
      denominator =
        - if denominator == "TP_frags": row["TP_frags"]
        - else: f"{denominator}_{rank}"  (e.g. correct/assigned)

    Output columns
    --------------
    package, protein_type, rank, numerator_key, denominator_key, numerator, denominator, fraction
    """
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    header = [
        "package",
        "protein_type",
        "rank",
        "numerator_key",
        "denominator_key",
        "numerator",
        "denominator",
        "fraction",
    ]
    lines: List[str] = ["\t".join(header)]

    n_written = 0

    for row in eval_stat_list:
        package = str(row.get("package", ""))
        protein_type = str(row.get("protein_type", ""))

        for r in ranks:
            num_key = f"{numerator_prefix}_{r}"
            num = row.get(num_key)

            # denominator: either global TP_frags or rank-specific prefix
            if denominator == "TP_frags":
                den_key = "TP_frags"
                den = row.get("TP_frags")
            else:
                den_key = f"{denominator}_{r}"
                den = row.get(den_key)

            # validate
            if not isinstance(num, (int, float)) or _is_nan(float(num)):
                continue
            if not isinstance(den, (int, float)) or _is_nan(float(den)) or float(den) <= 0:
                continue

            frac = float(num) / float(den)
            if _is_nan(frac) or frac < 0:
                continue

            lines.append(
                "\t".join(
                    [
                        package,
                        protein_type,
                        r,
                        num_key,
                        den_key,
                        str(float(num)),
                        str(float(den)),
                        str(frac),
                    ]
                )
            )
            n_written += 1

    out_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(
        f"[OK] Wrote TSV for rank fractions ({numerator_prefix}/{denominator}): "
        f"{out_tsv} (n={n_written})"
    )
    return out_tsv


def analyze_eval_stats(eval_stat_list):
    for package in eval_stat_list:
        print(package)
        # the package is a gpkg evaluation stat dict
        # was muss berechnet werden?


#
#
#

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


# ------------------------------------------------------------
# Generic plotters
# ------------------------------------------------------------
def plot_single_metric_boxplot(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        metric_key: str,
        title: str,
        ylabel: str,
        filename: str | None = None,
) -> Path | None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    values: List[float] = []
    for row in eval_stat_list:
        v = row.get(metric_key)
        if isinstance(v, (int, float)) and not _is_nan(float(v)):
            values.append(float(v))

    if not values:
        print(f"[WARN] No valid values for {metric_key}. Skipping plot.")
        return None

    if filename is None:
        filename = f"{metric_key}.png"

    fig, ax = plt.subplots()

    # Boxplot ohne Fliers (keine doppelten Punkte)
    ax.boxplot(
        values,
        showfliers=False,
        flierprops={"marker": ""},  # extra Sicherheit
    )

    # deterministic horizontal jitter
    n = len(values)
    if n == 1:
        xs = [1.0]
    else:
        width = 0.08
        xs = [1.0 + width * (i / (n - 1) - 0.5) for i in range(n)]

    # Punkte (genau einmal), klein, schwarz/weiß, hohl
    ax.scatter(
        xs,
        values,
        s=8,
        facecolors="none",
        edgecolors="black",
        linewidths=0.8,
        alpha=0.9,
    )

    # optional: title (du hattest es auskommentiert)
    # ax.set_title(title)

    ax.set_ylabel(ylabel)
    ax.set_ylim(0.0, 1.0)

    out_path = out_dir / filename
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Wrote plot: {out_path}")
    return out_path


def plot_rank_fraction_boxplots(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
        *,
        numerator_prefix: str,
        denominator: str,
        ranks: List[str] = RANKS,
        title: str | None = None,
        ylabel: str | None = None,
        filename: str | None = None,
) -> Path | None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_nan(x: float) -> bool:
        return isinstance(x, float) and math.isnan(x)

    def _jitter_x(x: float, n: int, width: float = 0.12) -> List[float]:
        # deterministic jitter centered on x
        if n <= 1:
            return [x]
        return [x + width * (i / (n - 1) - 0.5) for i in range(n)]

    data: List[List[float]] = []
    labels: List[str] = []
    counts: List[int] = []

    for r in ranks:
        num_key = f"{numerator_prefix}_{r}"
        den_key = "TP_frags" if denominator == "TP_frags" else f"{denominator}_{r}"

        vals: List[float] = []
        for row in eval_stat_list:
            num = row.get(num_key)
            den = row.get(den_key)

            if not isinstance(num, (int, float)) or _is_nan(float(num)):
                continue
            if not isinstance(den, (int, float)) or _is_nan(float(den)) or float(den) <= 0:
                continue

            frac = float(num) / float(den)
            if _is_nan(frac) or frac < 0:
                continue
            vals.append(frac)

        if vals:
            data.append(vals)
            labels.append(r)
            counts.append(len(vals))

    if not data:
        print(f"[WARN] No valid data for {numerator_prefix}/{denominator}. Skipping plot.")
        return None

    if filename is None:
        filename = f"rank_{numerator_prefix}_over_{denominator}.png"
    if title is None:
        title = f"{numerator_prefix} / {denominator} per rank"
    if ylabel is None:
        ylabel = f"{numerator_prefix} / {denominator}"

    fig, ax = plt.subplots()

    # Boxplots ohne Fliers (damit Punkte nicht doppelt erscheinen)
    ax.boxplot(
        data,
        labels=labels,
        showfliers=False,
        flierprops={"marker": ""},
        medianprops={
            "color": "black",
            "linewidth": 1.2,
        },
    )

    # Punkte: klein, schwarz/weiß, hohl, deterministischer jitter
    for i, vals in enumerate(data, start=1):
        xs = _jitter_x(float(i), len(vals), width=0.18)
        ax.scatter(
            xs,
            vals,
            s=8,
            facecolors="none",
            edgecolors="black",
            linewidths=0.8,
            alpha=0.9,
        )

    # ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0.0, 1.0)

    # n-Annotationen (bei fixem ylim ist das stabil)
    # n-Annotationen oberhalb der Achse (axes coordinates)
    for i, n in enumerate(counts, start=1):
        ax.text(
            (i - 0.5) / len(labels),  # x in axes coords
            1.02,  # leicht über der Achse
            f"n={n}",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=8,
        )

    out_path = out_dir / filename
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Wrote plot: {out_path}")
    return out_path


# ------------------------------------------------------------
# Batch routine: "all desired plots"
# ------------------------------------------------------------
def make_all_evaluation_plots(
        eval_stat_list: List[Dict[str, Any]],
        out_dir: Path,
) -> List[Path]:
    """
    Create all requested plots:
      - Overall: balanced accuracy, F1, MCC (boxplot each)
      - Placement: balanced accuracy, F1, MCC (boxplot each)
      - Rank-based fractions (one boxplot per rank):
          assigned/TP_frags
          correct/TP_frags
          correct/assigned
          wrong/assigned

    Returns
    -------
    List[Path]
        Paths of successfully written plot files.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: List[Path] = []

    # --- Overall metrics ---
    overall_specs: List[Tuple[str, str, str]] = [
        ("overall_balanced_accuracy", "Overall balanced accuracy", "overall balanced accuracy"),
        ("overall_F1", "Overall F1 score", "overall F1"),
        ("overall_MCC", "Overall MCC", "overall MCC"),
    ]
    overall_dir = out_dir / "overall"
    for key, title, ylabel in overall_specs:
        p = plot_single_metric_boxplot(
            eval_stat_list,
            overall_dir,
            metric_key=key,
            title=title,
            ylabel=ylabel,
            filename=f"{key}.png",
        )
        if p is not None:
            written.append(p)

    # --- Placement metrics ---
    placement_specs: List[Tuple[str, str, str]] = [
        ("placement_balanced_accuracy", "Placement balanced accuracy", "placement balanced accuracy"),
        ("placement_F1", "Placement F1 score", "placement F1"),
        ("placement_MCC", "Placement MCC", "placement MCC"),
    ]
    placement_dir = out_dir / "placement"
    for key, title, ylabel in placement_specs:
        p = plot_single_metric_boxplot(
            eval_stat_list,
            placement_dir,
            metric_key=key,
            title=title,
            ylabel=ylabel,
            filename=f"{key}.png",
        )
        if p is not None:
            written.append(p)

    # --- Rank fraction plots ---
    rank_specs: List[Tuple[str, str, str, str, str]] = [
        # numerator, denominator, subdir, title, ylabel
        ("assigned", "TP_frags", "rank", "Assigned fraction per rank (assigned / TP_frags)", "assigned / TP_frags"),
        ("correct", "TP_frags", "rank", "Correct fraction per rank (correct / TP_frags)", "correct / TP_frags"),
        ("correct", "assigned", "rank", "Correct among assigned per rank (correct / assigned)", "correct / assigned"),
        ("wrong", "assigned", "rank", "Wrong among assigned per rank (wrong / assigned)", "wrong / assigned"),
    ]
    rank_dir = out_dir / "rank"
    for num, den, _, title, ylabel in rank_specs:
        p = plot_rank_fraction_boxplots(
            eval_stat_list,
            rank_dir,
            numerator_prefix=num,
            denominator=den,
            title=title,
            ylabel=ylabel,
            filename=f"rank_{num}_over_{den}.png",
        )
        if p is not None:
            written.append(p)

    return written


# ------------------------------------------------------------
# CLI wrapper
# ------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize evaluation_stats.txt files (one row per package)."
    )
    ap.add_argument("root", help="Root directory to search recursively")
    ap.add_argument(
        "--file",
        default=DEFAULT_FILENAME,
        help="Evaluation stats filename (default: evaluation_stats.txt)",
    )
    ap.add_argument(
        "--sort",
        default="overall_MCC",
        help="Column to sort by (default: overall_MCC)",
    )
    ap.add_argument("--desc", action="store_true", help="Sort descending")
    ap.add_argument("--show-path", action="store_true", help="Include relative path column")
    args = ap.parse_args(argv)

    eval_list = summarize_eval_stats(
        root=Path(args.root),
        filename=args.file,
        sort=args.sort,
        desc=args.desc,
        show_path=args.show_path,
    )
    make_all_evaluation_plots(
        eval_stat_list=eval_list,
        out_dir=Path("/home/tomohisa/PycharmProjects/HMSS2/tests/project"),
    )
    print_eval_stats(eval_list)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

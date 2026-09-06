import os
from pathlib import Path
from typing import List

import pandas as pd
import hashlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

UNASSIGNED_LABEL = "Target genes absent"
OTHER_LABEL = "Collapsed categories"


def stable_color_for_label(label: str, palette: list):
    """
    Deterministic color assignment.
    Same label always receives same color.
    """
    digest = hashlib.md5(label.encode("utf-8")).hexdigest()
    idx = int(digest, 16) % len(palette)
    return palette[idx]


def plot_taxonomy_stacked_bars(
    input_tsv: str,
    outdir: str = "taxonomy_barplots",
    top_n: int = 20,
    min_category_fraction: float = 0.01,
    other_threshold_n_categories: int = 10,
) -> List[str]:
    """
    Create stacked horizontal barplots per taxonomic level.

    Parameters
    ----------
    input_tsv : str
        Path to summary_strain_variability_by_taxonomy.txt
    outdir : str
        Output directory for plots
    top_n : int
        Maximum number of taxa shown per taxonomic level
    min_category_fraction : float
        Categories contributing less than this fraction of the displayed total
        are merged into 'Other categories' if the number of categories exceeds
        `other_threshold_n_categories`
    other_threshold_n_categories : int
        Apply 'Other categories' only if the number of category columns is greater
        than this threshold

    Returns
    -------
    List[str]
        List of written plot file paths
    """
    if top_n < 1:
        raise ValueError("top_n must be >= 1")
    if not (0 <= min_category_fraction <= 1):
        raise ValueError("min_category_fraction must be between 0 and 1")
    if other_threshold_n_categories < 0:
        raise ValueError("other_threshold_n_categories must be >= 0")

    Path(outdir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_tsv, sep="\t")

    fixed_cols = [
        "taxonomic_level",
        "taxon_name",
        "genome_count_total_database",
    ]

    missing = [c for c in fixed_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    levels = df["taxonomic_level"].dropna().astype(str).unique().tolist()
    written: List[str] = []

    # Globales Farbmapping: gleiche Kategorie = gleiche Farbe in allen Plots
    all_category_cols = [c for c in df.columns if c not in fixed_cols]

    all_category_cols = sorted(all_category_cols, key=lambda s: s.casefold())

    # Farbpalette ohne Grautöne
    palette = list(cm.get_cmap("tab20b").colors) + list(cm.get_cmap("tab20c").colors)

    # explizit Grautöne entfernen
    filtered_palette = [
        c for c in palette if not (abs(c[0] - c[1]) < 0.08 and abs(c[1] - c[2]) < 0.08)
    ]

    color_map = {
        col: stable_color_for_label(col, filtered_palette) for col in all_category_cols
    }

    # Sonderkategorien fest setzen
    color_map[UNASSIGNED_LABEL] = "#d9d9d9"  # hellgrau
    color_map[OTHER_LABEL] = "#7f7f7f"  # deutlich dunkleres grau

    for level in levels:
        sub = df[df["taxonomic_level"] == level].copy()
        if sub.empty:
            continue

        category_cols = [c for c in sub.columns if c not in fixed_cols]

        sub["genome_count_total_database"] = (
            pd.to_numeric(sub["genome_count_total_database"], errors="coerce")
            .fillna(0)
            .astype(int)
        )

        for col in category_cols:
            sub[col] = pd.to_numeric(sub[col], errors="coerce").fillna(0).astype(int)

        # Summe aller zugeordneten Kategorien pro Taxon
        if category_cols:
            sub["assigned_total"] = sub[category_cols].sum(axis=1)
        else:
            sub["assigned_total"] = 0

        # Top-N nach Zahl der zugeordneten Genome
        sub = (
            sub.sort_values(
                by=["assigned_total", "taxon_name"],
                ascending=[False, True],
            )
            .head(top_n)
            .copy()
        )

        # Unassigned berechnen
        assigned_sum = sub[category_cols].sum(axis=1) if category_cols else 0
        sub[UNASSIGNED_LABEL] = sub["genome_count_total_database"] - assigned_sum
        sub[UNASSIGNED_LABEL] = sub[UNASSIGNED_LABEL].clip(lower=0)

        # "Other categories" nur anwenden, wenn genug Kategorien vorhanden sind
        kept_category_cols: List[str] = []
        small_category_cols: List[str] = []

        use_other = len(category_cols) > other_threshold_n_categories

        category_totals = {col: int(sub[col].sum()) for col in category_cols}

        total_assigned_counts = sum(category_totals.values())

        top_always_keep_n = 5

        top_category_cols = {
            col
            for col, total in sorted(
                category_totals.items(),
                key=lambda x: x[1],
                reverse=True,
            )[:top_always_keep_n]
        }

        if use_other and category_cols and total_assigned_counts > 0:
            for col in category_cols:
                frac = category_totals[col] / total_assigned_counts

                if col in top_category_cols:
                    kept_category_cols.append(col)
                elif frac < min_category_fraction:
                    small_category_cols.append(col)
                else:
                    kept_category_cols.append(col)
        else:
            kept_category_cols = list(category_cols)

        if use_other and small_category_cols:
            sub[OTHER_LABEL] = sub[small_category_cols].sum(axis=1)
            plot_cols = kept_category_cols + [OTHER_LABEL, UNASSIGNED_LABEL]
        else:
            plot_cols = kept_category_cols + [UNASSIGNED_LABEL]

        # Für barh: größte Zuordnung oben -> DataFrame aufsteigend sortieren,
        # weil die letzte Zeile oben erscheint
        sub = sub.sort_values(
            by=["assigned_total", "taxon_name"],
            ascending=[True, True],
        )

        n_rows = len(sub)
        fig_height = max(4, n_rows * 0.45)

        fig, ax = plt.subplots(figsize=(14, fig_height))

        y = range(len(sub))
        left = [0] * len(sub)

        # Plot color palette
        colors = cm.get_cmap("Set2", len(plot_cols))  # z.B. tab20, Set3, viridis

        for i, col in enumerate(plot_cols):
            values = sub[col].tolist()
            ax.barh(
                y,
                values,
                left=left,
                label=col,
                color=color_map[col],
            )
            left = [l + v for l, v in zip(left, values)]

        ax.set_yticks(list(y))
        ax.set_yticklabels(sub["taxon_name"].tolist())
        ax.set_xlabel("Number of genomes")
        ax.set_ylabel(level)
        ax.set_title(f"Taxonomic level: {level}")

        ax.legend(
            title="Category",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=8,
        )

        plt.tight_layout()

        outfile = os.path.join(outdir, f"stacked_barplot_{level}.svg")
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)

        written.append(outfile)

    return written

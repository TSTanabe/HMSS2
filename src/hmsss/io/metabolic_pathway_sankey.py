#!/usr/bin/env python3
"""
HMSS3 metabolic pathway Sankey plotting module.

Input table columns:
    genomeID
    species
    input
    output
    enzymes

Supported modes:
    input_output_enzymes
        Input -> Output -> Enzymes

    species_input_enzymes
        Species -> Input -> Enzymes

    species_output_enzymes
        Species -> Output -> Enzymes

    input_output_enzymes_species
        Input -> Output -> Enzymes -> Species

    all
        Create all four plots.

Color continuity:
    Link colors are assigned by the first layer of each plot and propagated
    through all downstream layers. Downstream links are drawn from full path
    rows to keep values traceable across merged nodes.

Main importable functions:
    create_sankey_input_output_enzymes(...)
    create_sankey_species_input_enzymes(...)
    create_sankey_species_output_enzymes(...)
    create_sankey_input_output_enzymes_species(...)
    create_all_standard_sankeys(...)

Dependencies:
    pandas
    plotly
    kaleido      # required for static svg/pdf/png export

Example CLI:
    python3 hmss3_metabolic_pathway_sankey.py \
        -i metabolic_pathway_report.tsv \
        -o sankey_out \
        --mode all \
        --top-n 40 \
        --top-n-species 40 \
        --split-middle-by-previous
"""

from __future__ import annotations

import argparse
import colorsys
import re
from pathlib import Path
from typing import Sequence

import pandas as pd
import plotly.graph_objects as go

try:
    from hmsss.core.logging import get_logger

    logger = get_logger(__name__)
except Exception:
    import logging

    logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = ["genomeID", "species", "input", "output", "enzymes"]

CONTRAST_HEX = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#6A3D9A",
    "#B15928",
    "#1B9E77",
    "#E7298A",
    "#66A61E",
    "#A6761D",
    "#666666",
]

MODE_TO_COLUMNS = {
    "input_output_enzymes": ("input", "output", "enzymes"),
    "species_input_enzymes": ("species", "input", "enzymes"),
    "species_output_enzymes": ("species", "output", "enzymes"),
    "input_output_enzymes_species": ("input", "output", "enzymes", "species"),
}

MODE_TO_FILENAME = {
    "input_output_enzymes": "sankey_input_output_enzymes",
    "species_input_enzymes": "sankey_species_input_enzymes",
    "species_output_enzymes": "sankey_species_output_enzymes",
    "input_output_enzymes_species": "sankey_input_output_enzymes_species",
}

MODE_TO_TITLE = {
    "input_output_enzymes": "Dominant transformations: input → output → enzymes",
    "species_input_enzymes": "Dominant pathways by species: species → input → enzymes",
    "species_output_enzymes": "Dominant pathways by species: species → output → enzymes",
    "input_output_enzymes_species": "Dominant transformations by species: input → output → enzymes → species",
}

DISPLAY_PREFIX = {
    "species": "Species",
    "input": "Input",
    "output": "Output",
    "enzymes": "Enzymes",
}


def normalize_text(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"\s+", " ", value)
    return value


def compact_label(value: str, max_chars: int) -> str:
    value = normalize_text(value)
    if len(value) <= max_chars:
        return value
    return value[: max_chars - 1] + "…"


def compact_enzymes(value: str, max_chars: int = 55) -> str:
    value = normalize_text(value)
    if len(value) <= max_chars:
        return value

    parts = [p.strip() for p in value.split(";") if p.strip()]
    if len(parts) <= 4:
        return value[: max_chars - 1] + "…"

    return ";".join(parts[:4]) + f";…(+{len(parts) - 4})"


def hex_to_rgb01(hex_color: str) -> tuple[float, float, float]:
    h = hex_color.lstrip("#")
    return (
        int(h[0:2], 16) / 255.0,
        int(h[2:4], 16) / 255.0,
        int(h[4:6], 16) / 255.0,
    )


def rgba_from_hex(hex_color: str, alpha: float, saturation_factor: float = 1.0) -> str:
    r, g, b = hex_to_rgb01(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    s = max(0.0, min(1.0, s * saturation_factor))
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return f"rgba({int(round(r2 * 255))},{int(round(g2 * 255))},{int(round(b2 * 255))},{alpha})"


def read_pathway_report(path: str | Path) -> pd.DataFrame:
    path = Path(path)

    df = pd.read_csv(path, sep="\t")
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=",")

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {missing}. Found columns: {list(df.columns)}"
        )

    for col in REQUIRED_COLUMNS:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df.loc[df["species"].isin(["", "nan", "NaN", "None"]), "species"] = "Unknown species"

    df = df[
        (df["genomeID"] != "")
        & (df["input"] != "")
        & (df["output"] != "")
        & (df["enzymes"] != "")
        ].copy()

    for col in ["species", "input", "output", "enzymes"]:
        df[col] = df[col].map(normalize_text)

    return df


def limit_species(
        df: pd.DataFrame,
        *,
        top_n_species: int = 40,
        other_label: str = "Other species",
) -> pd.DataFrame:
    out = df.copy()

    if top_n_species is None or top_n_species <= 0:
        return out

    top_species = (
        out.groupby("species")["genomeID"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_n_species)
        .index
    )
    keep = set(top_species.tolist())
    out["species"] = out["species"].where(out["species"].isin(keep), other_label)
    return out


def aggregate_paths(
        df: pd.DataFrame,
        *,
        columns: Sequence[str],
        top_n: int | None = 40,
        min_genomes: int = 1,
) -> pd.DataFrame:
    use_cols = list(columns)
    work = df.copy()

    if "enzymes" in use_cols:
        work["enzymes_display"] = work["enzymes"].map(compact_enzymes)
        use_cols = ["enzymes_display" if c == "enzymes" else c for c in use_cols]

    grouped = (
        work.groupby(use_cols, dropna=False)
        .agg(
            n_genomes=("genomeID", "nunique"),
            n_species=("species", "nunique"),
        )
        .reset_index()
    )

    grouped = grouped[grouped["n_genomes"] >= min_genomes].copy()
    grouped = grouped.sort_values(
        ["n_genomes"] + use_cols,
        ascending=[False] + [True] * len(use_cols),
    )

    if top_n is not None and top_n > 0:
        grouped = grouped.head(top_n).copy()

    grouped = grouped.rename(columns={"enzymes_display": "enzymes"})
    return grouped.reset_index(drop=True)


def weighted_barycentric_order(
        *,
        items: list[str],
        edges: pd.DataFrame,
        item_col: str,
        neighbor_col: str,
        neighbor_order: list[str],
        weight_col: str = "n_genomes",
) -> list[str]:
    neighbor_rank = {name: i for i, name in enumerate(neighbor_order)}

    tmp = edges[[item_col, neighbor_col, weight_col]].copy()
    tmp["_neighbor_rank"] = tmp[neighbor_col].map(neighbor_rank)
    tmp = tmp.dropna(subset=["_neighbor_rank"])
    tmp["_weighted_rank"] = tmp["_neighbor_rank"] * tmp[weight_col]

    score = (
        tmp.groupby(item_col, dropna=False)
        .agg(
            weighted_rank_sum=("_weighted_rank", "sum"),
            total=(weight_col, "sum"),
        )
        .reset_index()
    )

    if score.empty:
        return items

    score["_barycenter"] = score["weighted_rank_sum"] / score["total"]
    score = score.sort_values(
        ["_barycenter", "total", item_col],
        ascending=[True, False, True],
    )

    ordered = score[item_col].tolist()
    missing = [x for x in items if x not in set(ordered)]
    return ordered + missing


def optimize_layer_orders(path_df: pd.DataFrame, columns: Sequence[str]) -> list[list[str]]:
    columns = list(columns)

    first_col = columns[0]
    first_order = (
        path_df.groupby(first_col, dropna=False)["n_genomes"]
        .sum()
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    orders = [first_order]

    for left_col, right_col in zip(columns[:-1], columns[1:]):
        all_items = path_df[right_col].drop_duplicates().tolist()
        edges = (
            path_df.groupby([left_col, right_col], dropna=False)["n_genomes"]
            .sum()
            .reset_index()
        )
        right_order = weighted_barycentric_order(
            items=all_items,
            edges=edges,
            item_col=right_col,
            neighbor_col=left_col,
            neighbor_order=orders[-1],
        )
        orders.append(right_order)

    return orders


def make_layer_node_id(layer_index: int, column: str, value: str) -> str:
    return f"L{layer_index}:{column}:{value}"


def make_label(column: str, value: str) -> str:
    if column == "enzymes":
        value = compact_enzymes(value)
    elif column == "species":
        value = compact_label(value, 52)
    else:
        value = compact_label(value, 60)

    return f"{DISPLAY_PREFIX.get(column, column)}: {value}"


def _build_internal_columns_and_dataframe(
        path_df: pd.DataFrame,
        columns: list[str],
        *,
        split_middle_by_previous: bool,
) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    """
    Optionally duplicate the second layer by the previous layer.

    Returns:
        plot_df
        internal_columns
        lane_visible_column
            Maps lane column names to the original visible column name.
    """
    plot_df = path_df.copy()
    internal_columns = columns.copy()
    lane_visible_column: dict[str, str] = {}

    if split_middle_by_previous and len(columns) >= 3:
        prev_col = columns[0]
        middle_col = columns[1]
        lane_col = f"{middle_col}__lane"
        plot_df[lane_col] = plot_df[prev_col] + "|||@@@|||" + plot_df[middle_col]
        internal_columns[1] = lane_col
        lane_visible_column[lane_col] = middle_col

    return plot_df, internal_columns, lane_visible_column


def _visible_value(column: str, value: str, lane_visible_column: dict[str, str]) -> tuple[str, str]:
    """
    Return visible column name and visible value for normal or lane columns.
    """
    if column in lane_visible_column:
        _, visible_value = value.split("|||@@@|||", 1)
        return lane_visible_column[column], visible_value
    return column, value


def _dominant_first_value(
        plot_df: pd.DataFrame,
        *,
        first_col: str,
        column: str,
        value: str,
) -> str:
    sub = plot_df[plot_df[column] == value]
    return (
        sub.groupby(first_col)["n_genomes"]
        .sum()
        .sort_values(ascending=False)
        .index[0]
    )


def build_sankey_figure(
        path_df: pd.DataFrame,
        *,
        columns: Sequence[str],
        title: str,
        split_middle_by_previous: bool = False,
        link_alpha: float = 0.62,
        right_link_alpha: float = 0.42,
        node_pad: int = 12,
        node_thickness: int = 18,
        font_size: int = 12,
        width: int = 1600,
        height: int | None = None,
        compact_height_per_node: int = 24,
        min_height: int = 700,
        margin_left: int = 50,
        margin_right: int = 50,
        margin_top: int = 55,
        margin_bottom: int = 50,
) -> go.Figure:
    """
    Build a Sankey plot for arbitrary ordered columns.

    This function supports 3-layer and 4-layer plots.

    If split_middle_by_previous=True and there are at least three columns,
    the second layer is internally duplicated by the first layer. This makes
    streams look more like visual continuations.
    """
    columns = list(columns)

    if len(columns) not in {3, 4}:
        raise ValueError("This module expects three or four columns per Sankey plot.")

    plot_df, internal_columns, lane_visible_column = _build_internal_columns_and_dataframe(
        path_df,
        columns,
        split_middle_by_previous=split_middle_by_previous,
    )

    orders = optimize_layer_orders(plot_df, internal_columns)

    node_ids: list[str] = []
    labels: list[str] = []
    node_colors: list[str] = []

    first_layer_values = orders[0]
    first_to_hex = {
        value: CONTRAST_HEX[i % len(CONTRAST_HEX)]
        for i, value in enumerate(first_layer_values)
    }

    node_to_dominant_first: dict[str, str] = {}

    for layer_idx, (column, order) in enumerate(zip(internal_columns, orders)):
        for raw_value in order:
            node_id = make_layer_node_id(layer_idx, column, raw_value)
            node_ids.append(node_id)

            visible_col, visible_val = _visible_value(column, raw_value, lane_visible_column)
            labels.append(make_label(visible_col, visible_val))

            if layer_idx == 0:
                dominant = raw_value
            elif column in lane_visible_column:
                dominant, _visible = raw_value.split("|||@@@|||", 1)
            else:
                dominant = _dominant_first_value(
                    plot_df,
                    first_col=internal_columns[0],
                    column=column,
                    value=raw_value,
                )

            node_to_dominant_first[node_id] = dominant

            if layer_idx == 0:
                node_colors.append(rgba_from_hex(first_to_hex[dominant], 0.88, 1.0))
            elif layer_idx == 1:
                node_colors.append(rgba_from_hex(first_to_hex[dominant], 0.34, 0.45))
            elif layer_idx == 2:
                node_colors.append(rgba_from_hex(first_to_hex[dominant], 0.26, 0.35))
            else:
                node_colors.append(rgba_from_hex(first_to_hex[dominant], 0.20, 0.25))

    node_index = {node_id: i for i, node_id in enumerate(node_ids)}

    source: list[int] = []
    target: list[int] = []
    value: list[int] = []
    customdata: list[str] = []
    link_colors: list[str] = []

    # IMPORTANT FOR COLOR CONTINUITY:
    # Draw links from complete path rows, not from collapsed adjacent node pairs.
    # If different first-layer values pass through the same downstream node pair,
    # they remain separate links with their original color instead of being merged
    # and recolored by a dominant node. This makes values traceable through the
    # whole Sankey path.
    rank_maps = [
        {v: i for i, v in enumerate(order)}
        for order in orders
    ]
    plot_df = plot_df.copy()
    for idx, col in enumerate(internal_columns):
        plot_df[f"_rank_{idx}"] = plot_df[col].map(rank_maps[idx])
    plot_df = plot_df.sort_values(
        [f"_rank_{i}" for i in range(len(internal_columns))] + ["n_genomes"],
        ascending=[True] * len(internal_columns) + [False],
    )

    for _, row in plot_df.iterrows():
        first_value = row[internal_columns[0]]
        first_hex = first_to_hex[first_value]

        for layer_idx, (left_col, right_col) in enumerate(zip(internal_columns[:-1], internal_columns[1:])):
            left_value = row[left_col]
            right_value = row[right_col]

            left_node = make_layer_node_id(layer_idx, left_col, left_value)
            right_node = make_layer_node_id(layer_idx + 1, right_col, right_value)

            color_alpha = link_alpha if layer_idx == 0 else right_link_alpha

            source.append(node_index[left_node])
            target.append(node_index[right_node])
            value.append(int(row["n_genomes"]))

            visible_left_col, visible_left_val = _visible_value(left_col, left_value, lane_visible_column)
            visible_right_col, visible_right_val = _visible_value(right_col, right_value, lane_visible_column)

            customdata.append(
                f"{DISPLAY_PREFIX.get(visible_left_col, visible_left_col)}: {visible_left_val}<br>"
                f"{DISPLAY_PREFIX.get(visible_right_col, visible_right_col)}: {visible_right_val}<br>"
                f"{DISPLAY_PREFIX.get(columns[0], columns[0])}: {first_value}<br>"
                f"Genomes: {row['n_genomes']}"
            )

            # Keep exactly the same color family through all downstream links.
            # Later links are only slightly more transparent; they are not
            # recolored by output/enzyme/species nodes.
            link_colors.append(rgba_from_hex(first_hex, color_alpha, 0.88))

    if height is None:
        height = max(min_height, compact_height_per_node * len(labels))

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=node_pad,
                    thickness=node_thickness,
                    line=dict(color="black", width=0.4),
                    label=labels,
                    color=node_colors,
                    hovertemplate="%{label}<extra></extra>",
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_colors,
                    customdata=customdata,
                    hovertemplate="%{customdata}<extra></extra>",
                ),
            )
        ]
    )

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        font=dict(size=font_size),
        margin=dict(
            l=margin_left,
            r=margin_right,
            t=margin_top,
            b=margin_bottom,
        ),
        width=width,
        height=height,
    )

    return fig


def write_figure(fig: go.Figure, output_file: str | Path) -> None:
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    suffix = output_file.suffix.lower()
    if suffix == ".html":
        fig.write_html(str(output_file))
    elif suffix in {".svg", ".pdf", ".png", ".jpeg", ".jpg", ".webp"}:
        fig.write_image(str(output_file))
    else:
        raise ValueError(
            f"Unsupported output extension: {suffix}. "
            "Use .svg, .pdf, .png, .html, .jpeg, .jpg, or .webp."
        )


def create_sankey_for_mode(
        report_file: str | Path,
        output_file: str | Path,
        *,
        mode: str,
        top_n: int | None = 40,
        min_genomes: int = 1,
        top_n_species: int = 40,
        split_middle_by_previous: bool = False,
        title: str | None = None,
        width: int = 1600,
        height: int | None = None,
        node_pad: int = 8,
        node_thickness: int = 22,
        font_size: int = 16,
        compact_height_per_node: int = 24,
        min_height: int = 700,
        link_alpha: float = 0.62,
        right_link_alpha: float = 0.42,
) -> go.Figure:
    if mode not in MODE_TO_COLUMNS:
        raise ValueError(f"Unknown mode {mode!r}. Valid modes: {sorted(MODE_TO_COLUMNS)}")

    columns = MODE_TO_COLUMNS[mode]
    df = read_pathway_report(report_file)

    if "species" in columns:
        df = limit_species(df, top_n_species=top_n_species)

    path_df = aggregate_paths(
        df,
        columns=columns,
        top_n=top_n,
        min_genomes=min_genomes,
    )

    if path_df.empty:
        raise ValueError(f"No rows left for mode {mode!r} after filtering.")

    fig = build_sankey_figure(
        path_df,
        columns=columns,
        title=title or MODE_TO_TITLE[mode],
        split_middle_by_previous=split_middle_by_previous,
        width=width,
        height=height,
        node_pad=node_pad,
        node_thickness=node_thickness,
        font_size=font_size,
        compact_height_per_node=compact_height_per_node,
        min_height=min_height,
        link_alpha=link_alpha,
        right_link_alpha=right_link_alpha,
    )
    write_figure(fig, output_file)
    return fig


def create_sankey_input_output_enzymes(
        report_file: str | Path,
        output_file: str | Path,
        **kwargs,
) -> go.Figure:
    return create_sankey_for_mode(
        report_file,
        output_file,
        mode="input_output_enzymes",
        **kwargs,
    )


def create_sankey_species_input_enzymes(
        report_file: str | Path,
        output_file: str | Path,
        **kwargs,
) -> go.Figure:
    return create_sankey_for_mode(
        report_file,
        output_file,
        mode="species_input_enzymes",
        **kwargs,
    )


def create_sankey_species_output_enzymes(
        report_file: str | Path,
        output_file: str | Path,
        **kwargs,
) -> go.Figure:
    return create_sankey_for_mode(
        report_file,
        output_file,
        mode="species_output_enzymes",
        **kwargs,
    )


def create_sankey_input_output_enzymes_species(
        report_file: str | Path,
        output_file: str | Path,
        **kwargs,
) -> go.Figure:
    return create_sankey_for_mode(
        report_file,
        output_file,
        mode="input_output_enzymes_species",
        **kwargs,
    )


def create_all_standard_sankeys(
        report_file: str | Path,
        output_dir: str | Path,
        *,
        extension: str = "svg",
        top_n: int | None = 40,
        min_genomes: int = 1,
        top_n_species: int = 40,
        split_middle_by_previous: bool = False,
        width: int = 1600,
        height: int | None = None,
        node_pad: int = 8,
        node_thickness: int = 22,
        font_size: int = 16,
        compact_height_per_node: int = 24,
        min_height: int = 700,
        link_alpha: float = 0.62,
        right_link_alpha: float = 0.42,
) -> dict[str, Path]:
    output_dir = Path(output_dir)
    extension = extension.lstrip(".")

    outputs: dict[str, Path] = {}

    for mode in MODE_TO_COLUMNS:
        out_file = output_dir / f"{MODE_TO_FILENAME[mode]}.{extension}"
        create_sankey_for_mode(
            report_file,
            out_file,
            mode=mode,
            top_n=top_n,
            min_genomes=min_genomes,
            top_n_species=top_n_species,
            split_middle_by_previous=split_middle_by_previous,
            width=width,
            height=height,
            node_pad=node_pad,
            node_thickness=node_thickness,
            font_size=font_size,
            compact_height_per_node=compact_height_per_node,
            min_height=min_height,
            link_alpha=link_alpha,
            right_link_alpha=right_link_alpha,
        )
        outputs[mode] = out_file

    return outputs


def create_metabolic_pathway_plots(config) -> None:
    """Create standard metabolic pathway Sankey plots after pathway report generation.

    Expected config attributes:
        pathway_report_file
        metabolic_pathway_directory

    Optional config attributes:
        metabolic_sankey_top_n
        metabolic_sankey_min_genomes
        metabolic_sankey_top_n_species
        metabolic_sankey_split_middle
        metabolic_sankey_width
        metabolic_sankey_height
        metabolic_sankey_node_pad
        metabolic_sankey_node_thickness
        metabolic_sankey_font_size
        metabolic_sankey_height_per_node
        metabolic_sankey_min_height
        metabolic_sankey_link_alpha
        metabolic_sankey_right_link_alpha
    """
    pathway_report_file = getattr(config, "pathway_report_file", None)
    pathway_plot_dir = getattr(config, "metabolic_pathway_directory", None)

    if not pathway_report_file:
        logger.info("No metabolic pathway report file configured; skipping metabolic pathway plots")
        return

    pathway_report_path = Path(pathway_report_file)
    if not pathway_report_path.is_file():
        logger.info("No metabolic pathway report found at %s; skipping metabolic pathway plots", pathway_report_path)
        return

    if pathway_report_path.stat().st_size == 0:
        logger.info("Metabolic pathway report is empty; skipping metabolic pathway plots")
        return

    if pathway_plot_dir is None:
        pathway_plot_dir = pathway_report_path.parent

    pathway_plot_dir = Path(pathway_plot_dir)
    pathway_plot_dir.mkdir(parents=True, exist_ok=True)

    try:
        outputs = create_all_standard_sankeys(
            pathway_report_path,
            pathway_plot_dir,
            extension="svg",
            top_n=int(getattr(config, "metabolic_sankey_top_n", 45)),
            min_genomes=int(getattr(config, "metabolic_sankey_min_genomes", 1)),
            top_n_species=int(getattr(config, "metabolic_sankey_top_n_species", 35)),
            split_middle_by_previous=bool(getattr(config, "metabolic_sankey_split_middle", True)),
            width=int(getattr(config, "metabolic_sankey_width", 1600)),
            height=getattr(config, "metabolic_sankey_height", None),
            node_pad=int(getattr(config, "metabolic_sankey_node_pad", 12)),
            node_thickness=int(getattr(config, "metabolic_sankey_node_thickness", 26)),
            font_size=int(getattr(config, "metabolic_sankey_font_size", 14)),
            compact_height_per_node=int(getattr(config, "metabolic_sankey_height_per_node", 24)),
            min_height=int(getattr(config, "metabolic_sankey_min_height", 700)),
            link_alpha=float(getattr(config, "metabolic_sankey_link_alpha", 0.62)),
            right_link_alpha=float(getattr(config, "metabolic_sankey_right_link_alpha", 0.42)),
        )
    except ImportError as exc:
        logger.warning("Could not create metabolic pathway Sankey plots because a dependency is missing: %s", exc)
        return
    except Exception as exc:
        logger.warning("Could not create metabolic pathway Sankey plots: %s", exc)
        return

    for mode, out_file in outputs.items():
        logger.info("Created metabolic pathway Sankey plot [%s]: %s", mode, out_file)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create HMSS3 metabolic pathway Sankey plots.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="metabolic_pathway_report.tsv")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output file or output directory for --mode all.")
    parser.add_argument("--mode", choices=["all"] + sorted(MODE_TO_COLUMNS), default="all")
    parser.add_argument("--extension", default="svg", help="Used only with --mode all. Default: svg.")
    parser.add_argument("--top-n", type=int, default=45, help="Top N paths. Use 0 for all.")
    parser.add_argument("--min-genomes", type=int, default=1)
    parser.add_argument("--top-n-species", type=int, default=35, help="Top N species. Use 0 for all.")
    parser.add_argument(
        "--split-middle-by-previous",
        action="store_true",
        help=(
            "Duplicate the second layer by the first layer to make flows look like continuations. "
            "Especially useful for input→output→enzyme plots."
        ),
    )
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=None)
    parser.add_argument("--node-pad", type=int, default=12)
    parser.add_argument("--node-thickness", type=int, default=26)
    parser.add_argument("--font-size", type=int, default=14)
    parser.add_argument("--compact-height-per-node", type=int, default=24)
    parser.add_argument("--min-height", type=int, default=700)
    parser.add_argument("--link-alpha", type=float, default=0.62)
    parser.add_argument("--right-link-alpha", type=float, default=0.42)
    parser.add_argument("--title", default=None, help="Only used for single-mode output.")

    args = parser.parse_args()
    top_n = None if args.top_n == 0 else args.top_n

    if args.mode == "all":
        outputs = create_all_standard_sankeys(
            args.input,
            args.output,
            extension=args.extension,
            top_n=top_n,
            min_genomes=args.min_genomes,
            top_n_species=args.top_n_species,
            split_middle_by_previous=args.split_middle_by_previous,
            width=args.width,
            height=args.height,
            node_pad=args.node_pad,
            node_thickness=args.node_thickness,
            font_size=args.font_size,
            compact_height_per_node=args.compact_height_per_node,
            min_height=args.min_height,
            link_alpha=args.link_alpha,
            right_link_alpha=args.right_link_alpha,
        )
        for mode, path in outputs.items():
            print(f"[OK] {mode}: {path}")
    else:
        create_sankey_for_mode(
            args.input,
            args.output,
            mode=args.mode,
            top_n=top_n,
            min_genomes=args.min_genomes,
            top_n_species=args.top_n_species,
            split_middle_by_previous=args.split_middle_by_previous,
            title=args.title,
            width=args.width,
            height=args.height,
            node_pad=args.node_pad,
            node_thickness=args.node_thickness,
            font_size=args.font_size,
            compact_height_per_node=args.compact_height_per_node,
            min_height=args.min_height,
            link_alpha=args.link_alpha,
            right_link_alpha=args.right_link_alpha,
        )
        print(f"[OK] {args.mode}: {args.output}")


if __name__ == "__main__":
    main()

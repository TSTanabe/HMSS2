#!/usr/bin/python
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import Wedge
from typing import Dict, Any, Optional, Set, Tuple, List

import math


def plot_taxonomy_cooccurrence_network(
    output_file: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    allowed_types: Optional[List[str]] = None,
    allowed_levels: Optional[Set[str]] = None,
    preferred_order: Optional[List[str]] = None,
    *,
    min_cooccurrence: int = 3,
    node_size_min: float = 0.05,
    node_size_max: float = 0.12,
) -> None:
    """
    Zeichnet ein Protein-Typ-Co-Occurrence-Netzwerk:

    - Knoten: Protein-Typen (Domains), bestimmt wie in plot_taxonomy_summary_bubbles()
    - Kanten: Co-Occurrence in denselben Genomen, Gewicht = #Genome mit beiden Typen
    - Knotengröße: Anzahl Genome mit diesem Typ
    - Knoten als Pie-Charts: Anteil der Genome pro Phylum (oder anderem Tax-Level)

    Diese Variante:
    - benutzt ein gewichtetes spring_layout (weight='weight')
    - skaliert das Layout global auf eine kompakte Ausdehnung
    - verzichtet auf teure pairwise Anti-Overlap-Iterationen
    """

    # ---- Taxonomie-Level festlegen (für Pies) ----
    default_level = "Phylum"
    tax_levels = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    if allowed_levels is None or len(allowed_levels) == 0:
        level_for_pies = default_level
    else:
        if default_level in allowed_levels:
            level_for_pies = default_level
        else:
            level_for_pies = sorted(allowed_levels)[0]

    def _norm(v: str | None) -> str:
        return "" if v is None else str(v).strip()

    def _protein_type(p: Any) -> str:
        """wie in plot_taxonomy_summary_bubbles"""
        dom_string = ""
        try:
            dom_string = p.get_domains() or ""
        except Exception:
            dom_string = ""
        dom_string = dom_string.strip()
        if dom_string:
            return dom_string

        for dom in p.domains:
            name = getattr(dom, "domain", "") or ""
            name = str(name).strip()
            if name:
                return name

        return "UNK"

    # ---- allowed_types normalisieren ----
    if allowed_types is not None and len(allowed_types) == 0:
        allowed_types = None
    allowed_types_set: Optional[Set[str]] = (
        set(allowed_types) if allowed_types is not None else None
    )

    # ------------------------------------------------------------------
    # 1) Presence/Absence auf Genome x Protein-Typ-Ebene erstellen
    # ------------------------------------------------------------------
    genome_to_types: Dict[str, Set[str]] = {}
    type_to_genomes: Dict[str, Set[str]] = {}

    for prot in protein_dict.values():
        genome_id = getattr(prot, "genomeID", None)
        if not genome_id:
            continue

        ptype = _protein_type(prot) or "UNK"
        if allowed_types_set is not None and ptype not in allowed_types_set:
            continue

        genome_to_types.setdefault(genome_id, set()).add(ptype)
        type_to_genomes.setdefault(ptype, set()).add(genome_id)

    if not type_to_genomes:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No data for given filters", ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_file, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    protein_types = sorted(type_to_genomes.keys(), key=str.casefold)

    # ------------------------------------------------------------------
    # 2) Co-Occurrence zählen: für jeden Genome die Paare seiner Typen
    # ------------------------------------------------------------------
    co_counts: Dict[Tuple[str, str], int] = {}

    for genome, types in genome_to_types.items():
        types = sorted(types)
        n = len(types)
        for i in range(n):
            for j in range(i + 1, n):
                a = types[i]
                b = types[j]
                key = (a, b)
                co_counts[key] = co_counts.get(key, 0) + 1

    # ------------------------------------------------------------------
    # 3) Phylum-Verteilung pro Protein-Typ (für Pies)
    # ------------------------------------------------------------------
    all_phyla: Set[str] = set()
    type_to_phylum_counts: Dict[str, Dict[str, int]] = {}

    for ptype, genomes in type_to_genomes.items():
        p_counts: Dict[str, int] = {}
        for gid in genomes:
            rec = taxon_dict.get(gid, {})
            ph = _norm(rec.get(level_for_pies))
            if not ph:
                ph = "NA"
            p_counts[ph] = p_counts.get(ph, 0) + 1
        type_to_phylum_counts[ptype] = p_counts
        all_phyla.update(p_counts.keys())

    all_phyla = sorted(all_phyla)

    # ------------------------------------------------------------------
    # 4) Netzwerk mit networkx aufbauen
    # ------------------------------------------------------------------
    G = nx.Graph()

    type_counts = {ptype: len(genomes) for ptype, genomes in type_to_genomes.items()}
    for ptype in protein_types:
        G.add_node(ptype, genome_count=type_counts[ptype])

    for (a, b), cnt in co_counts.items():
        if cnt >= min_cooccurrence:
            G.add_edge(a, b, weight=cnt)

    if G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(
            0.5,
            0.5,
            "No nodes after filtering / min_cooccurrence",
            ha="center",
            va="center",
        )
        ax.axis("off")
        fig.savefig(output_file, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    # ------------------------------------------------------------------
    # 5) Layout berechnen (gewichtetes Spring, kompakt skaliert)
    # ------------------------------------------------------------------
    n_nodes = max(1, G.number_of_nodes())
    # kleinere k -> kompakter, weight='weight' -> starke Co-Occurrence zieht zusammen
    k = 0.8 / math.sqrt(n_nodes + 1)

    pos = nx.spring_layout(
        G,
        k=k,
        weight="weight",
        seed=42,
        iterations=150,
        scale=1.0,
        center=(0.0, 0.0),
    )

    # globales Rescaling, damit das Netzwerk im Bereich [-max_extent, max_extent] bleibt
    coords = np.array([pos[n] for n in G.nodes()], dtype=float)
    max_abs = np.abs(coords).max()
    max_extent = 1.2
    if max_abs > 0:
        scale = max_extent / max_abs
        coords *= scale
        for n, (x, y) in zip(G.nodes(), coords):
            pos[n] = np.array([x, y])

    # ------------------------------------------------------------------
    # 6) Radii im Layoutraum für die Pies (klein!)
    # ------------------------------------------------------------------
    max_count = max(type_counts.values()) if type_counts else 1
    xs = [pos[p][0] for p in G.nodes()]
    ys = [pos[p][1] for p in G.nodes()]
    span_x = max(xs) - min(xs) if xs else 1.0
    span_y = max(ys) - min(ys) if ys else 1.0
    span = max(span_x, span_y) or 1.0

    # Basisradius im Layoutraum: bewusst konservativ
    base_radius = (
        span * 0.02 / math.sqrt(n_nodes)
    )  # ggf. 0.015 oder 0.01, falls noch zu groß
    r_min = base_radius * 0.6
    r_max = base_radius * 1.4

    def node_radius(ptype: str) -> float:
        if max_count <= 0:
            return r_min
        frac = type_counts[ptype] / max_count
        return r_min + (r_max - r_min) * frac

    radii = {ptype: node_radius(ptype) for ptype in G.nodes()}

    # ------------------------------------------------------------------
    # 7) Farben für Phyla & Pie-Helpers
    # ------------------------------------------------------------------
    cmap = plt.get_cmap("tab20c")
    phylum_colors = {p: cmap(i % 20) for i, p in enumerate(all_phyla)}

    def draw_pie_node(
        ax, center, fractions, colors, radius, edgecolor="black", linewidth=0.5
    ):
        x, y = center
        start_angle = 0.0
        for frac, color in zip(fractions, colors):
            if frac <= 0:
                continue
            theta1 = start_angle * 360.0
            theta2 = (start_angle + frac) * 360.0
            wedge = Wedge(
                (x, y),
                radius,
                theta1,
                theta2,
                facecolor=color,
                edgecolor=edgecolor,
                linewidth=linewidth,
            )
            ax.add_patch(wedge)
            start_angle += frac

    # ------------------------------------------------------------------
    # 8) Plotten
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 8))

    # Kanten
    edges = G.edges(data=True)
    if edges:
        weights = [d.get("weight", 1) for (_, _, d) in edges]
        max_w = max(weights)
        widths = [0.5 + 2.0 * (w / max_w) for w in weights]
        nx.draw_networkx_edges(
            G, pos, ax=ax, width=widths, alpha=0.5, edge_color="grey"
        )

    # Knoten als Pies
    for ptype in G.nodes():
        cx, cy = pos[ptype]
        p_counts = type_to_phylum_counts.get(ptype, {})
        counts_list = [p_counts.get(p, 0) for p in all_phyla]
        total = sum(counts_list)
        if total > 0:
            fractions = [c / total for c in counts_list]
        else:
            fractions = [0.0 for _ in all_phyla]
        colors = [phylum_colors[p] for p in all_phyla]
        r = radii[ptype]
        draw_pie_node(ax, (cx, cy), fractions, colors, radius=r)

    # Labels auf alle Knoten
    for ptype in G.nodes():
        cx, cy = pos[ptype]
        ax.text(
            cx,
            cy,
            ptype,
            ha="center",
            va="center",
            fontsize=6,
            color="black",
            clip_on=True,  # wichtig: nicht außerhalb der Achsen wachsen lassen
        )

    # Legende für Phyla
    handles = []
    for p in all_phyla:
        handles.append(
            plt.Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                markersize=6,
                markerfacecolor=phylum_colors[p],
                markeredgecolor="black",
                label=p,
            )
        )
    ax.legend(
        handles=handles,
        title=level_for_pies,
        loc="upper right",
        fontsize=8,
        frameon=False,
    )

    ax.set_title(
        "Protein-type co-occurrence network\n"
        f"(node pies = genome {level_for_pies} composition)"
    )
    ax.axis("off")

    # etwas manueller Rand, statt tight_layout + bbox_inches="tight"
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02)

    # kein bbox_inches="tight", geringere DPI → viel weniger RAM-Bedarf
    fig.savefig(output_file, dpi=300)
    plt.close(fig)

#!/usr/bin/python
import matplotlib.pyplot as plt

from typing import Dict, Any, Optional, Set, Tuple, List


def plot_taxonomy_summary_bubbles(
        output_file: str,
        protein_dict: Dict[str, Any],
        taxon_dict: Dict[str, Dict[str, str]],
        allowed_types: Optional[List[str]] = None,
        allowed_levels: Optional[Set[str]] = None,
        preferred_order: Optional[List[str]] = None,
) -> None:
    """
    Erzeugt eine Abbildung einer "Presence/Absence"-ähnlichen Matrix:

    - Zeilen: Taxa (standardmäßig Phyla)
    - Spalte 1: Anzahl Genome im Taxon
    - weitere Spalten: Protein-Typen
    - Zellen: Kreise, deren Größe und Farbe den Anteil der Genome in diesem Taxon
      codieren, die mindestens einen Treffer für diesen Protein-Typ haben.

    Args:
        preferred_order: Defines the order of all or some protein types
        output_file: Pfad für die Bilddatei (z.B. "summary.png").
        protein_dict: {proteinID: Protein-Objekt mit .genomeID, .get_domains(), ...}
        taxon_dict: {genomeID: {tax_level: name, ...}} – gleiche Struktur wie in
                    _output_taxonomy_summary.
        allowed_types: optionale **Liste** von Protein-Typen. Wenn gesetzt und nicht leer,
                       werden NUR diese Typen verwendet – in GENAU dieser Reihenfolge.
        allowed_levels: optionales Set von Taxonomie-Ebenen (z.B. {"Phylum"}).
                        Wenn None oder leer, wird standardmäßig {"Phylum"} genutzt.
    """

    tax_levels = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    def _is_empty_or_na(v: str | None) -> bool:
        if v is None:
            return True
        s = str(v).strip()
        return not s or s.upper() == "NA"

    def _norm(v: str | None) -> str:
        return "" if v is None else str(v).strip()

    def _protein_type(p: Any) -> str:
        """
        Protein-Typ aus Protein-Objekt ableiten (wie in _output_taxonomy_summary):

          1. p.get_domains() wenn nicht leer
          2. erster Domain-Name aus p.get_domains_dict()
          3. sonst 'UNK'
        """
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

    # ------------------------------------------------------------------
    # Filter-Argumente normalisieren
    # ------------------------------------------------------------------
    # allowed_types: None oder leere Liste => kein Filter
    if allowed_types is not None and len(allowed_types) == 0:
        allowed_types = None

    # Für Membership-Check später ein Set bauen, Reihenfolge aber aus der LISTE
    allowed_types_set: Optional[Set[str]] = (
        set(allowed_types) if allowed_types is not None else None
    )

    if allowed_levels is not None and len(allowed_levels) == 0:
        allowed_levels = None

    # Default: nur Phylum, wenn kein Level-Filter explizit gesetzt
    if allowed_levels is None:
        levels_to_use = {"Phylum"}
    else:
        levels_to_use = set(allowed_levels)

    if preferred_order is not None and len(preferred_order) == 0:
        preferred_order = None

    # ------------------------------------------------------------------
    # 1) Counts vorbereiten – PRESENCE pro Genom
    # ------------------------------------------------------------------
    # key = (level, taxon_name)
    # counts[key] = { protein_type: set(genomeID) }
    counts: Dict[Tuple[str, str], Dict[str, Set[str]]] = {}
    # genomes_per_taxon[key] = set(genomeID)
    genomes_per_taxon: Dict[Tuple[str, str], Set[str]] = {}
    all_types: Set[str] = set()

    for prot in protein_dict.values():
        genome_id = getattr(prot, "genomeID", None)
        if not genome_id:
            continue

        rec = taxon_dict.get(genome_id)
        if not isinstance(rec, dict):
            continue

        ptype = _protein_type(prot) or "UNK"

        # Typen-Filter: Membership über Set, Reihenfolge später über Liste
        if allowed_types_set is not None and ptype not in allowed_types_set:
            continue

        all_types.add(ptype)

        # Nur gewünschte Ebenen betrachten
        for level in tax_levels:
            if level not in levels_to_use:
                continue

            tax_name = _norm(rec.get(level))
            if _is_empty_or_na(tax_name):
                continue

            key = (level, tax_name)

            # Set aller Genome im Taxon
            genomes_per_taxon.setdefault(key, set()).add(genome_id)

            # Presence pro Protein-Typ
            level_counts = counts.setdefault(key, {})
            level_counts.setdefault(ptype, set()).add(genome_id)

    # Nix zu plotten?
    if not counts:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No data for given filters", ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_file, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    # ------------------------------------------------------------------
    # 2) Protein-Typen: Reihenfolge
    # ------------------------------------------------------------------
    if allowed_types is not None:
        # expliziter Limiter: nur diese Typen, in genau dieser Reihenfolge (Duplikate entfernen)
        seen_cols: set[str] = set()
        protein_types: List[str] = []
        for t in allowed_types:
            if t in seen_cols:
                continue
            seen_cols.add(t)
            if t in all_types:
                protein_types.append(t)
    else:
        if preferred_order is not None:
            # 1) gewünschte Reihenfolge (requests) vorn
            seen_cols: set[str] = set()
            leading: List[str] = []
            for t in preferred_order:
                if t in seen_cols:
                    continue
                if t in all_types:
                    seen_cols.add(t)
                    leading.append(t)

            # 2) alle restlichen Typen alphabetisch hinten anhängen
            remaining = sorted(
                [t for t in all_types if t not in seen_cols],
                key=lambda s: s.casefold(),
            )
            protein_types = leading + remaining
        else:
            # ganz klassisch: alphabetisch
            protein_types = sorted(all_types, key=lambda s: s.casefold())

    # ------------------------------------------------------------------
    # 3) Zeilen vorbereiten (nur gewählte Levels)
    # ------------------------------------------------------------------
    rows: List[Tuple[Tuple[str, str], int, Dict[str, Set[str]]]] = []
    level_index = {lvl: i for i, lvl in enumerate(tax_levels)}

    for key, level_counts in counts.items():
        level, tax_name = key
        if level not in levels_to_use:
            continue
        genome_count = len(genomes_per_taxon.get(key, set()))
        rows.append((key, genome_count, level_counts))

    # sortiert nach Level-Hierarchie, dann nach übergeordneten Taxa (Lineage),
    # dann Taxon-Name. Damit werden z.B. Archaea/Bacteria sauber getrennt.
    def _row_sort_key(row):
        key, _, _ = row
        level, tax_name = key

        # Level-Priorität bleibt wie vorher
        lvl_rank = level_index.get(level, 999)

        gids = genomes_per_taxon.get(key, set())
        if not gids:
            return (lvl_rank, "", "", "", "", "", "", "", tax_name.casefold())

        # bis zu welchem Rang soll sortiert werden? (inkl. aktuellem level)
        max_i = level_index.get(level, 999)

        def norm(v):
            return "" if v is None else str(v).strip()

        # repräsentative Lineage bestimmen (deterministisch: lexikographisch kleinste)
        lineage_candidates = []
        for gid in gids:
            rec = taxon_dict.get(gid)
            if not isinstance(rec, dict):
                continue
            lineage = [norm(rec.get(lvl)).casefold() for lvl in tax_levels[: max_i + 1]]
            lineage_candidates.append(tuple(lineage))

        if not lineage_candidates:
            return (lvl_rank, "", "", "", "", "", "", "", tax_name.casefold())

        rep = min(lineage_candidates)
        pad = ("",) * (len(tax_levels) - len(rep))

        # Key: (Level-Rank, Superkingdom, Phylum, ..., bis Level, [padding], Taxonname)
        return (lvl_rank,) + rep + pad + (tax_name.casefold(),)

    rows.sort(key=_row_sort_key)

    multiple_levels = len({lvl for (lvl, _), _, _ in rows}) > 1

    # ------------------------------------------------------------------
    # 4) Koordinaten & Bubbles
    # ------------------------------------------------------------------
    n_rows = len(rows)
    n_cols = len(protein_types) + 1  # +1 für "#Genomes"

    x_positions = list(range(n_cols))
    x_type_positions = list(range(1, n_cols))

    fractions: List[float] = []
    xs: List[int] = []
    ys: List[int] = []
    sizes: List[float] = []
    colors: List[float] = []
    genome_counts: List[int] = []

    for row_idx, (key, genome_count, level_counts) in enumerate(rows):
        genome_counts.append(genome_count)

        for col_idx, ptype in enumerate(protein_types, start=1):
            gset_type = level_counts.get(ptype, set())
            count = len(gset_type)
            frac = count / genome_count if genome_count > 0 else 0.0

            # 0 → NICHT zeichnen
            if frac <= 0:
                continue

            fractions.append(frac)
            xs.append(col_idx)
            ys.append(row_idx)
            colors.append(frac)

    if fractions:
        max_frac = max(fractions)
    else:
        max_frac = 0.0

    min_radius = 2.0
    max_radius = 12.0

    for frac in fractions:
        if max_frac <= 0:
            radius = min_radius
        else:
            radius = min_radius + (max_radius - min_radius) * (frac / max_frac)
        sizes.append(radius ** 2)

    fig_width = max(4, 1.0 + 0.6 * n_cols)
    fig_height = max(3, 0.4 * n_rows + 1.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sc = ax.scatter(
        xs,
        ys,
        s=sizes if sizes else 0,
        c=colors if colors else 0,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        edgecolors="k",
        linewidths=0.3,
    )

    y_ticks = list(range(n_rows))
    if multiple_levels:
        y_labels = [f"{lvl}: {name}" for (lvl, name), _, _ in rows]
    else:
        y_labels = [name for (lvl, name), _, _ in rows]

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    x_labels = ["#Genomes", *protein_types]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")

    for row_idx, gc in enumerate(genome_counts):
        ax.text(
            0,
            row_idx,
            str(gc),
            ha="center",
            va="center",
            fontsize=8,
        )

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.invert_yaxis()

    ax.set_xlabel("Protein types")
    ax.set_ylabel("Taxa")

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Fraction of genomes with ≥1 hit")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _row_lineage_key(key: Tuple[str, str]) -> tuple:
    """
    Liefert einen Sortierschlüssel, der übergeordnete Ebenen berücksichtigt.

    key = (level, tax_name)
    Wir nehmen die GenomeIDs der Zeile und lesen deren Taxonomie aus taxon_dict.
    Dann bilden wir eine Lineage-Tuple (Superkingdom, Phylum, Class, ...),
    und wählen als repräsentativen Schlüssel die lexikographisch kleinste Lineage.
    """
    level, tax_name = key
    gids = genomes_per_taxon.get(key, set())
    if not gids:
        return (999, "", "", "", "", "", "", "")

    # bis zu welchem Rang soll sortiert werden?
    # (bei Phylum nur Superkingdom+Phylum, bei Class bis Class, usw.)
    try:
        max_i = level_index.get(level, 999)
    except Exception:
        max_i = 999

    def norm(v):
        return "" if v is None else str(v).strip()

    lineage_candidates = []
    for gid in gids:
        rec = taxon_dict.get(gid)
        if not isinstance(rec, dict):
            continue
        lineage = [norm(rec.get(lvl)) for lvl in tax_levels[: max_i + 1]]
        # case-insensitive sort, aber originalwerte bleiben egal: wir nutzen casefold
        lineage_candidates.append(tuple(x.casefold() for x in lineage))

    if not lineage_candidates:
        return (999, "", "", "", "", "", "", "")

    # Repräsentant: kleinste Lineage (stabil, deterministisch)
    rep = min(lineage_candidates)

    # Primär: Level-Hierarchie, Sekundär: Lineage
    # rep enthält bis max_i; wir padden für einheitliche Tuple-Länge
    pad = ("",) * (len(tax_levels) - len(rep))
    return (level_index.get(level, 999),) + rep + pad

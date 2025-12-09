#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple
from matplotlib.patches import FancyArrowPatch
from matplotlib.gridspec import GridSpec

# ===================== Datenstrukturen =====================

@dataclass
class Gene:
    start: int
    end: int
    strand: int = 1          # +1 oder -1
    color: str = "tab:blue"
    label: str = ""


@dataclass
class ContigSegment:
    name: str
    length: int
    features: List[Gene] = field(default_factory=list)
    track_start: float = 0.0   # in "bp" (Plotkoordinaten)
    track_end: float = 0.0

    def add_gene(self, gene: Gene):
        if not (0 <= gene.start < gene.end <= self.length):
            raise ValueError(
                f"Gene {gene} außerhalb des Contigs '{self.name}' (0, {self.length})"
            )
        self.features.append(gene)

    def transform_coord(self, x: float) -> float:
        """Contig-Koordinate (0..length) -> globale Track-Koordinate (track_start..track_end)."""
        if not (0 <= x <= self.length):
            raise ValueError(
                f"Position {x} außerhalb des Contigs '{self.name}' (0, {self.length})"
            )
        frac = x / self.length
        return self.track_start + frac * (self.track_end - self.track_start)


# ===================== Track =====================

class Track:
    def __init__(
        self,
        name: str,
        contigs: Dict[str, int],
        space: float = 0.02,
        labelsize: float = 12,
        labelmargin: float = 0.01,
        align_label: bool = True,
        gene_colors: Dict[str, str] | None = None,
    ):
        self.name = name
        self.space = space
        self.segments: List[ContigSegment] = [
            ContigSegment(cname, length) for cname, length in contigs.items()
        ]
        self.ax = None
        self.xlim = (0.0, 1.0)

        # Label-Parameter ähnlich GenomeViz FeatureTrack
        self._label: str | None = None
        self._labelsize = labelsize
        self._labelmargin = labelmargin
        self._align_label = align_label

        # Mapping Genlabel -> Farbe
        self._gene_colors: Dict[str, str] = gene_colors or {}

        self._layout_segments()

    # -------- Layout --------

    @property
    def label(self) -> str:
        return self.name if self._label is None else self._label

    @property
    def total_length(self) -> int:
        if not self.segments:
            return 0
        return int(self.segments[-1].track_end)

    def _layout_segments(self):
        """Setzt track_start/track_end der Contigs hintereinander mit Gaps."""
        if not self.segments:
            return

        max_len = max(s.length for s in self.segments)
        if self.space < 1.0:
            gap = max(100, int(max_len * self.space))
        else:
            gap = int(self.space)

        pos = 0
        for seg in self.segments:
            seg.track_start = pos
            seg.track_end = pos + seg.length
            pos = seg.track_end + gap

        # xlim erstmal track-lokal; wird später auf globalen Max-Wert normiert
        self.xlim = (0, self.segments[-1].track_end)

    # -------- API --------

    def add_gene(
        self,
        contig: str,
        start: int,
        end: int,
        strand: int = 1,
        color: str = "tab:blue",
        label: str = "",
    ):
        """Gen auf bestimmtem Contig hinzufügen (Koordinaten relativ zum Contig).

        Farbe:
        - falls `label` im gene_colors-Dict vorkommt -> diese Farbe
        - sonst `color`-Argument (Default: 'tab:blue')
        """
        for seg in self.segments:
            if seg.name == contig:
                # Farbe ggf. aus Mapping überschreiben
                gene_color = self._gene_colors.get(label, color)
                seg.add_gene(Gene(start, end, strand, gene_color, label))
                return
        raise ValueError(f"Unbekannter Contig '{contig}' in Track '{self.name}'")

    def set_ax(self, ax):
        self.ax = ax
        self.ax.axis("off")

    # -------- Zeichnen --------

    def plot(self):
        if self.ax is None:
            raise RuntimeError("Axis für Track nicht gesetzt.")

        ax = self.ax

        # Y-Positionen innerhalb dieses Tracks (0..1)
        baseline_y = 0.35          # Contig-Linie
        contig_label_y = 0.05      # Contig-Label
        gene_y = 0.55              # Gen-Pfeile
        gene_label_y = 0.70        # Gen-Labels
        major_tick = 0.06          # Start/Ende
        minor_tick = 0.03          # 1000-bp
        self.baseline_y = baseline_y  # für externe Referenz (falls nötig)

        # nur innerhalb dieses Panels
        ax.set_xlim(*self.xlim)
        ax.set_ylim(0.0, 1.0)
        ax.axis("off")

        # kombinierter Transform: x in Daten, y in Axes-Koordinaten
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

        # ----- Contigs, Ticks, Contig-Labels -----
        for seg in self.segments:
            # Hauptlinie
            ax.plot(
                [seg.track_start, seg.track_end],
                [baseline_y, baseline_y],
                color="black",
                linewidth=1.0,
                transform=trans,
            )

            # Start/Ende: volle Ticks
            for bp in (0, seg.length):
                x = seg.transform_coord(bp)
                ax.plot(
                    [x, x],
                    [baseline_y - major_tick, baseline_y + major_tick],
                    color="black",
                    linewidth=0.6,
                    transform=trans,
                )

            # alle 1000 bp: halbe Ticks nach unten
            step = 1000
            for bp in range(step, seg.length, step):
                x = seg.transform_coord(bp)
                ax.plot(
                    [x, x],
                    [baseline_y, baseline_y - minor_tick],
                    color="black",
                    linewidth=0.4,
                    transform=trans,
                )

            # Contig-Beschriftung
            mid = (seg.track_start + seg.track_end) / 2.0
            ax.text(
                mid,
                contig_label_y,
                f"{seg.name} (0 - {seg.length} bp)",
                ha="center",
                va="bottom",
                fontsize=6,
                transform=trans,
                clip_on=False,
            )

        # ----- Gene: Pfeile & Labels -----
        for seg in self.segments:
            for g in seg.features:
                x_start = seg.transform_coord(g.start)
                x_end = seg.transform_coord(g.end)
                if x_end == x_start:
                    continue

                if g.strand >= 0:
                    x_tail = x_start
                    x_head = x_end
                else:
                    x_tail = x_end
                    x_head = x_start

                # Pfeil über die volle Gen-Länge
                arrow = FancyArrowPatch(
                    (x_tail, gene_y),
                    (x_head, gene_y),
                    arrowstyle="Simple,head_width=0.8,head_length=0.9,tail_width=0.4",
                    mutation_scale=10,
                    linewidth=0,          # Körper ist gefüllte Fläche, keine Linien nötig
                    facecolor=g.color,
                    edgecolor="none",
                    transform=trans,
                    clip_on=True,
                )
                ax.add_patch(arrow)


                ax.add_patch(arrow)

                # Label
                if g.label:
                    center = (x_start + x_end) / 2.0
                    ax.text(
                        center,
                        gene_label_y,
                        g.label,
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        rotation=45,
                        transform=trans,
                        clip_on=False,
                    )

        # ----- Track-Label (wie GenomeViz FeatureTrack) -----
        if self._align_label:
            x = -self._labelmargin
        else:
            first_seg = self.segments[0]
            first_seg_start_x = first_seg.track_start / (self.xlim[1] or 1)
            x = first_seg_start_x - self._labelmargin

        ax.text(
            x,
            0.5,
            self.label,
            ha="right",
            va="center_baseline",
            fontsize=self._labelsize,
            transform=ax.transAxes,
        )


# ===================== GenomeFigure =====================

class GenomeFigure:
    """
    Einfacher Genome-Plotter:
    - Höhe skaliert mit Anzahl der Tracks (wie GenomeViz)
    - Nutzung von GridSpec (ein Panel pro Track)
    - Export via .save()
    """

    def __init__(
        self,
        fig_width: float = 12.0,
        fig_track_height_cm: float = 2.0,
        bp_per_inch: float = 800.0,
    ):
        """
        Parameters
        ----------
        fig_width : float
            Breite der Abbildung in Zoll.
        fig_track_height_cm : float
            Mindesthöhe pro Track (cm).
        bp_per_inch : float
            Horizontale Skalierung: Basenpaare pro Zoll.
            Kleinere Werte -> breiterer Plot.
        """
        self._fig_width = fig_width
        self._fig_track_height_in = max(0.1, fig_track_height_cm / 2.54)
        self._bp_per_inch = bp_per_inch
        self.tracks: List[Track] = []

    # ---- API ----

    def add_track(
        self,
        name: str,
        contigs: Dict[str, int],
        space: float = 0.02,
        gene_colors: Dict[str, str] | None = None,
    ) -> Track:
        t = Track(name, contigs, space=space, gene_colors=gene_colors)
        self.tracks.append(t)
        return t

    # ---- Layout-Helfer ----

    def _update_xlims(self):
        if not self.tracks:
            return
        max_x = max(t.total_length for t in self.tracks) or 1
        for t in self.tracks:
            t.xlim = (0, max_x)

    @property
    def figsize(self) -> Tuple[float, float]:
        """Figurgröße (wie bei GenomeViz): Breite fix, Höhe = track_height * n_tracks."""
        n = max(1, len(self.tracks))
        max_x = max((t.xlim[1] for t in self.tracks), default=1)
        width_by_bp = max_x / self._bp_per_inch
        width = max(self._fig_width, width_by_bp, 6.0)
        height = n * self._fig_track_height_in
        return (width, height)

    # ---- Plot & Save ----

    def plot(self, dpi: int = 150):
        if not self.tracks:
            raise ValueError("Keine Tracks zum Plotten vorhanden.")

        self._update_xlims()
        fig_w, fig_h = self.figsize
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)

        gs = GridSpec(
            nrows=len(self.tracks),
            ncols=1,
            height_ratios=[1.0] * len(self.tracks),
            figure=fig,
        )
        gs.update(
            left=0.10,   # etwas mehr Platz links für Track-Label
            right=0.99,
            top=0.98,
            bottom=0.02,
            hspace=0.4,  # genug Abstand zwischen Tracks
        )

        axes = []
        for idx, track in enumerate(self.tracks):
            ax = fig.add_subplot(gs[idx])
            track.set_ax(ax)
            track.plot()
            axes.append(ax)

        return fig, axes

    def save(self, filename: str, dpi: int = 300, pad_inches: float = 0.1, **kwargs):
        fig, _ = self.plot(dpi=dpi)
        fig.savefig(filename, dpi=dpi, bbox_inches="tight", pad_inches=pad_inches, **kwargs)
        plt.close(fig)

def plot_gene_cluster_summary(
    output_file: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    allowed_types: Optional[List[str]] = None,
    allowed_levels: Optional[Set[str]] = None,
) -> None:
    """
    Erzeugt eine GenomeFigure mit einem Track pro Genom und zeichnet dort die
    Gene als Pfeile in ihren Contigs.

    - Ein Track = ein Genome (genomeID)
    - Ein Contig-Segment = gene_contig, Länge = max(gene_end) in diesem Contig
    - Ein Gen-Pfeil = (gene_start, gene_end, gene_strand)
    - Label/Farbe = abgeleiteter 'protein type' (Domain-Name)

    allowed_types:
        Liste von Protein-Typen (Domains), die dargestellt werden sollen.
        None oder leere Liste => keine Filterung.

    allowed_levels:
        Wird hier nicht zum Filtern verwendet; könnte später genutzt werden,
        um die Reihenfolge/ Auswahl der Genomes nach Taxonomie zu steuern.
    """

    # ---- Helper: Protein-Typ bestimmen (wie in _output_taxonomy_summary) ----
    def _protein_type(p: Any) -> str:
        # 1) get_domains()
        dom_string = ""
        try:
            dom_string = p.get_domains() or ""
        except Exception:
            dom_string = ""
        dom_string = dom_string.strip()
        if dom_string:
            return dom_string

        # 2) erster Name aus get_domains_dict()
        try:
            dct = p.get_domains_dict()
        except Exception:
            dct = {}

        for dom in dct.values():
            name = getattr(dom, "domain", "") or ""
            name = str(name).strip()
            if name:
                return name

        # 3) Fallback
        return "UNK"

    # ---- 1) Proteine nach Genom und Contig gruppieren ----
    # genomeID -> contig -> Liste von (proteinObj, ptype)
    genome_to_contigs: Dict[str, Dict[str, List[tuple[Any, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )

    # Alle beobachteten Typen sammeln (für Farbzuordnung)
    all_types: Set[str] = set()

    for p in protein_dict.values():
        gid = getattr(p, "genomeID", None)
        contig = getattr(p, "gene_contig", None)
        start = getattr(p, "gene_start", None)
        end = getattr(p, "gene_end", None)

        if not gid or contig is None or start is None or end is None:
            continue

        ptype = _protein_type(p)
        if allowed_types is not None and ptype not in allowed_types:
            continue

        genome_to_contigs[gid][contig].append((p, ptype))
        all_types.add(ptype)

    if not genome_to_contigs:
        # Nichts zu plotten
        return

    # ---- 2) Farb-Mapping nach Protein-Typ ----
    # Reihenfolge stabil halten (damit Farben reproduzierbar sind)
    sorted_types = sorted(all_types, key=lambda s: s.casefold())
    cmap = plt.get_cmap("tab20", max(1, len(sorted_types)))
    gene_colors: Dict[str, Any] = {
        ptype: cmap(i) for i, ptype in enumerate(sorted_types)
    }

    # ---- 3) Genomes sortieren (z.B. nach Taxonomie, dann genomeID) ----
    tax_levels = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    def genome_sort_key(gid: str):
        rec = taxon_dict.get(gid, {}) or {}
        tax_key = tuple(
            (rec.get(lvl) or "").strip().casefold() for lvl in tax_levels
        )
        return (*tax_key, gid)

    genomes_sorted = sorted(genome_to_contigs.keys(), key=genome_sort_key)

    # ---- 4) Figure anlegen und Tracks hinzufügen ----
    gv = GenomeFigure(fig_width=12, fig_track_height_cm=3.0, bp_per_inch=800)

    for gid in genomes_sorted:
        contig_map = genome_to_contigs[gid]

        # Contig-Längen: max(gene_end) in jedem Contig
        contig_lengths: Dict[str, int] = {}
        for cname, plist in contig_map.items():
            max_end = 0
            for p, _ptype in plist:
                end = int(getattr(p, "gene_end", 0) or 0)
                if end > max_end:
                    max_end = end
            if max_end <= 0:
                continue
            contig_lengths[cname] = max_end

        if not contig_lengths:
            continue

        # Track-Label: z.B. Species + genomeID
        rec = taxon_dict.get(gid, {}) or {}
        label_lines: List[str] = [gid]

        if allowed_types is not None:
            for lvl in allowed_levels:
                if lvl == "Species":
                    continue
                name = (rec.get(lvl) or "").strip()
                if name:
                        label_lines.append(name)
        species = (rec.get("Species") or "").strip()
        if species:
            label_lines.append(species)
        track_name = "\n".join(label_lines)
        track = gv.add_track(track_name, contig_lengths, gene_colors=gene_colors)

        # Gene hinzufügen, sortiert nach Contig und Start
        for cname, plist in sorted(
            contig_map.items(),
            key=lambda kv: kv[0],
        ):
            # pro Contig nach Startposition sortieren
            plist_sorted = sorted(
                plist,
                key=lambda t: int(getattr(t[0], "gene_start", 0) or 0),
            )
            for p, ptype in plist_sorted:
                start = int(getattr(p, "gene_start", 0) or 0)
                end = int(getattr(p, "gene_end", 0) or 0)
                strand = getattr(p, "gene_strand", "+")
                if strand == "+":
                    strand = 1
                else:
                    strand = -1
                # Basisschutz
                if end <= start:
                    continue

                track.add_gene(
                    cname,
                    start,
                    end,
                    strand=strand,
                    label=ptype,  # Label = Protein-Typ (Domain-Name)
                )

    # ---- 5) Speichern ----
    gv.save(output_file, dpi=300)
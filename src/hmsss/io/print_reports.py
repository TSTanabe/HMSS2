#!/usr/bin/python
import os
from collections import defaultdict, Counter
from typing import Any, Dict, List, Optional, Set, Tuple, Iterable
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def _get_cluster_info(clusterID, cluster_dict):
    csb_val = ""
    if clusterID:
        if clusterID in cluster_dict:
            # cluster.get_cluster_list(",") liefert i. d. R. [clusterID, keyword, completeness, csb]
            cl = cluster_dict[clusterID].get_cluster_list(",")
            # defensiv extrahieren:
            out_cluster_id = cl[0] if len(cl) >= 1 else clusterID
            csb_val = cl[3] if len(cl) >= 4 else (cl[-1] if len(cl) >= 2 else "")
        else:
            out_cluster_id = clusterID
    else:
        out_cluster_id = ""

    if not csb_val:
        csb_val = ""

    return out_cluster_id, csb_val


def _output_genome_report(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, Any],
    genomeID: str = "",
    writemode: str = "w",
    taxon_divider: str = ".",
) -> None:
    """
    Writes the main genome hit table (TSV).

    Columns:
      genomeID, proteinID, domains, dom_scores, dom_coords, contig, gene_start, gene_end,
      gene_strand, locustag, clusterID, csb,
      Superkingdom, Phylum, Class, Order, Family, Genus, Species
    """
    taxon_cols = [
        "genomeID",
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    header_cols = [
        "genomeID",
        "proteinID",
        "domains",
        "dom_scores",
        "dom_coords",
        "contig",
        "gene_start",
        "gene_end",
        "gene_strand",
        "locustag",
        "selection_comment",
        "alternative hit",
        "clusterID",
        *taxon_cols,
    ]
    header = "\t".join(header_cols)

    # sortiert nach genomeID, contig, start
    proteinID_list = sorted(
        protein_dict,
        key=lambda x: (
            protein_dict[x].genomeID,
            protein_dict[x].gene_contig,
            protein_dict[x].gene_start,
        ),
    )

    with open(output_filepath, writemode) as writer:
        writer.write(header + "\n")

        for proteinID in proteinID_list:
            protein = protein_dict[proteinID]
            # protein.get_protein_list() erwartete Reihenfolge laut Docstring:
            # [proteinID, get_domains(), get_domain_scores(), get_domain_coordinates(),
            #  gene_contig, gene_start, gene_end, gene_strand, gene_locustag]
            pl = protein.get_protein_list()

            # clusterID + csb (ohne keyword/completeness), needed for parse report routine

            cluster_id = protein.clusterID
            out_cluster_id, csb_val = _get_cluster_info(cluster_id, cluster_dict)

            gid = protein.genomeID
            rec = taxon_dict.get(gid, {})  # falls Genome nicht im Dict: leeres Dict

            taxon_levels = [
                rec.get("Superkingdom"),
                rec.get("Phylum"),
                rec.get("Class"),
                rec.get("Order"),
                rec.get("Family"),
                rec.get("Genus"),
                rec.get("Species"),
            ]

            row = [
                protein.genomeID,  # genomeID
                pl[0],  # proteinID
                pl[1],  # domains
                pl[2],  # dom_scores
                pl[3],  # dom_coords
                pl[4],  # contig
                pl[5],  # gene_start
                pl[6],  # gene_end
                pl[7],  # gene_strand
                pl[8],  # locustag
                protein.get_selection_comment_csv(),
                protein.alternative_hit,
                out_cluster_id,
                csb_val,
                *taxon_levels,
            ]
            writer.write("\t".join(map(str, row)) + "\n")
    return


def _output_protein_taxonomy(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    writemode: str = "w",
) -> None:
    """
    Writes a 2-column TSV file:
      proteinID    taxonomy_lineage

    Taxonomy lineage format:
      k__Superkingdom; p__Phylum; c__Class; o__Order; f__Family; g__Genus; s__Species

    Missing or 'NA' taxonomy entries are replaced with the proteinID.
    """

    header = "proteinID\ttaxonomy"

    levels = [
        ("Superkingdom", "k__"),
        ("Phylum", "p__"),
        ("Class", "c__"),
        ("Order", "o__"),
        ("Family", "f__"),
        ("Genus", "g__"),
        ("Species", "s__"),
    ]

    proteinID_list = sorted(
        protein_dict,
        key=lambda x: (
            protein_dict[x].genomeID,
            protein_dict[x].gene_contig,
            protein_dict[x].gene_start,
        ),
    )

    with open(output_filepath, writemode) as writer:
        writer.write(header + "\n")

        for proteinID in proteinID_list:
            protein = protein_dict[proteinID]
            gid = protein.genomeID

            rec = taxon_dict.get(gid, {})

            def get_tax_val(level: str) -> str:
                val = rec.get(level)
                # Replace None, empty string, or 'NA' (case-insensitive) with proteinID
                if not val or str(val).strip().upper() == "NA":
                    return proteinID
                return val

            lineage = "; ".join(
                f"{prefix}{get_tax_val(level)}" for level, prefix in levels
            )

            writer.write(f"{proteinID}\t{lineage}\n")


def _output_taxonomy_summary(
    output_file: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    allowed_types: Optional[List[str]] = None,
    preferred_order: Optional[List[str]] = None,
) -> None:
    """
    Summarise protein-type presence per taxonomic level and taxon.

    Output format (TSV):
        taxonomic_level    taxon_name    genome_count    <prot_type_1> ...

    - taxonomic_level: one of Superkingdom, Phylum, Class, Order, Family, Genus, Species
    - taxon_name: name at this level
    - genome_count: number of distinct genomes that have at least one of the
      considered proteins and are assigned to this taxon at this level.
    - each subsequent column: number of genomes in this taxon that have at least
      one protein of that type (presence/absence per genome).
    - Protein types are derived from `protein_dict` and sorted alphabetically.
    - Rows are sorted hierarchically by level (Superkingdom..Species), then
      alphabetically by taxon_name.
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
        Derive a 'protein type' label from a Protein object.

        Heuristik:
          1. p.get_domains() wenn nicht leer
          2. erster Domain-Name aus p.get_domains_dict()
          3. sonst 'UNK'
        """
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

    # ------------------------------------------------------------------
    # Normalize allowed_types
    # ------------------------------------------------------------------
    if allowed_types is not None:
        # leeres Set bedeutet: "nichts filtern"
        if len(allowed_types) == 0:
            allowed_types = None  # treat as no filtering

    # ------------------------------------------------------------------
    # 1) Counts vorbereiten – jetzt PRESENCE pro Genom
    # ------------------------------------------------------------------
    # (level, taxon_name) -> { protein_type: set(genomeID) }
    counts: Dict[Tuple[str, str], Dict[str, Set[str]]] = {}
    # (level, taxon_name) -> set(genomeID)  (alle Genome in diesem Taxon)
    genomes_per_taxon: Dict[Tuple[str, str], Set[str]] = {}

    all_types: Set[str] = set()

    for prot in protein_dict.values():
        genome_id = getattr(prot, "genomeID", None)
        if not genome_id:
            continue

        rec = taxon_dict.get(genome_id)
        if not isinstance(rec, dict):
            # keine Taxonomie für dieses Genom
            continue

        ptype = _protein_type(prot) or "UNK"

        # Wenn Filtering aktiv: nur erlaubte Protein-Typen sammeln
        if allowed_types is not None and ptype not in allowed_types:
            continue

        all_types.add(ptype)

        for level in tax_levels:
            tax_name = _norm(rec.get(level))
            if _is_empty_or_na(tax_name):
                continue

            key = (level, tax_name)

            # Set aller Genome in diesem Taxon
            gset_tax = genomes_per_taxon.setdefault(key, set())
            gset_tax.add(genome_id)

            # Presence pro Protein-Typ: set(genomeID)
            level_counts = counts.setdefault(key, {})
            gset_type = level_counts.setdefault(ptype, set())
            gset_type.add(genome_id)

    # Wenn gar keine Daten → minimale Datei mit Header
    if not counts:
        with open(output_file, "w", newline="") as out:
            out.write("taxonomic_level\ttaxon_name\tgenome_count\n")
        return

    # ------------------------------------------------------------------
    # 2) Protein-Spalten alphabetisch oder nach vorgabe
    # ------------------------------------------------------------------
    if allowed_types is not None:
        # expliziter Limiter: nur diese Typen, in genau dieser Reihenfolge
        seen_cols: Set[str] = set()
        protein_types: List[str] = []
        for t in allowed_types:
            if t in seen_cols:
                continue
            if t in all_types:
                seen_cols.add(t)
                protein_types.append(t)
    else:
        if preferred_order is not None:
            # zuerst preferred_order (requests), dann Rest alphabetisch
            seen_cols: Set[str] = set()
            leading: List[str] = []
            for t in preferred_order:
                if t in seen_cols:
                    continue
                if t in all_types:
                    seen_cols.add(t)
                    leading.append(t)

            remaining = sorted(
                [t for t in all_types if t not in seen_cols],
                key=lambda s: s.casefold(),
            )
            protein_types = leading + remaining
        else:
            # klassisch: alle beobachteten Typen alphabetisch
            protein_types = sorted(all_types, key=lambda s: s.casefold())

    # ------------------------------------------------------------------
    # 3) Zeilen bauen und sortieren
    # ------------------------------------------------------------------
    rows: List[Tuple[str, str, int, Dict[str, Set[str]]]] = []
    level_index = {lvl: i for i, lvl in enumerate(tax_levels)}

    for (level, tax_name), level_counts in counts.items():
        genome_count = len(genomes_per_taxon.get((level, tax_name), set()))
        rows.append((level, tax_name, genome_count, level_counts))

    # sortiert nach Level-Hierarchie, dann Taxon-Name alphabetisch
    rows.sort(key=lambda r: (level_index.get(r[0], 999), r[1].casefold()))

    # ------------------------------------------------------------------
    # 4) Schreiben
    # ------------------------------------------------------------------
    header = ["taxonomic_level", "taxon_name", "genome_count", *protein_types]
    with open(output_file, "w", newline="") as out:
        out.write("\t".join(header) + "\n")
        for level, tax_name, genome_count, level_counts in rows:
            # pro Protein-Typ: Anzahl Genome mit ≥1 Kopie dieses Typs
            counts_per_type = [
                str(len(level_counts.get(t, set()))) for t in protein_types
            ]
            out.write(
                "\t".join([level, tax_name, str(genome_count), *counts_per_type]) + "\n"
            )


def _output_unique_taxonomy_table(
    output_file: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
) -> str:
    """
    Write a non-redundant taxonomy table from a structured taxonomy dict.

    - Deduplicates by the 7 ranks (Superkingdom..Species).
    - Sorts only by the taxonomy columns.
    - Appends the count of genomes per unique lineage as the last column.

    Returns:
        Path to the written file.
    """
    unique_file = output_file.replace(".txt", "")

    tax_cols = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    def _norm(v: str | None) -> str:
        return "" if v is None else v.strip()

    # Zähle Genome pro eindeutiger 7er-Linie
    counts: Dict[Tuple[str, ...], int] = {}
    hit_genomes = {p.genomeID for p in protein_dict.values()}
    for gid, rec in taxon_dict.items():
        if gid not in hit_genomes:
            continue

        if not isinstance(rec, dict):
            continue
        tail7 = tuple(_norm(rec.get(k)) for k in tax_cols)
        counts[tail7] = counts.get(tail7, 0) + 1

    # Zeilen: 7 Taxonomie-Spalten + Count
    rows: List[Tuple[str, ...]] = [(*taxo, str(n)) for taxo, n in counts.items()]

    # Sortierung: nur nach den 7 Taxonomie-Spalten; leere/NA ans Ende
    def _is_empty_or_na(s: str) -> bool:
        s = (s or "").strip()
        return s == "" or s.upper() == "NA"

    def _sort_key(row: Tuple[str, ...]):
        taxo_cols = row[:7]
        return tuple((_is_empty_or_na(s), (s or "").casefold()) for s in taxo_cols)

    rows_sorted = sorted(rows, key=_sort_key)

    # Schreiben
    header = [*tax_cols, "genome_count"]
    with open(unique_file, "w", newline="") as uf:
        uf.write("\t".join(header) + "\n")
        for row in rows_sorted:
            uf.write("\t".join(map(str, row)) + "\n")

    return unique_file


def _output_strain_variability_by_species(
    directory: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    required_domains: Set[str],
    *,
    unknown_label: str = "Unknown",
) -> str:
    import os
    from collections import defaultdict

    os.makedirs(directory, exist_ok=True)
    outpath = os.path.join(directory, "summary_strain_variability_by_species.txt")

    def _species_label(rec: Dict[str, str] | None) -> str:
        if rec is None:
            return unknown_label
        sp = (rec.get("Species") or "").strip()
        if not sp or sp.upper() == "NA":
            deepest = (rec.get("DeepestValue") or "").strip() or "NA"
            return f"{deepest}_prokaryote"
        return sp

    # 1) Genome -> Species
    genome_to_species: Dict[str, str] = {}
    species_to_genomes: Dict[str, Set[str]] = defaultdict(set)
    for gid, rec in (taxon_dict or {}).items():
        label = _species_label(rec)
        genome_to_species[gid] = label
        species_to_genomes[label].add(gid)

    # 2) Domain-Typen pro Genom aus protein_dict
    genome_to_domains: Dict[str, Set[str]] = defaultdict(set)
    for prot in protein_dict.values():
        gid = getattr(prot, "genomeID", None)
        if not gid:
            continue
        for dom in prot.get_domains_dict().values():
            name = getattr(dom, "domain", None)
            if name:
                genome_to_domains[gid].add(name)

        if gid not in genome_to_species:
            genome_to_species[gid] = unknown_label
            species_to_genomes[unknown_label].add(gid)

    # Hilfsfunktion: Hat irgendein Genom der Spezies mindestens eine Required-Domain?
    def _has_any_required(gids: Set[str]) -> bool:
        if not required_domains:
            return True  # nichts zu filtern, wenn keine Targets vorgegeben sind
        for gid in gids:
            if genome_to_domains.get(gid, set()) & required_domains:
                return True
        return False

    # 3) Zählen je Species (mit Filter "komplett abwesend" überspringen)
    rows: list[Tuple[str, int, int, int, float]] = []
    for species, gids in species_to_genomes.items():
        if not _has_any_required(gids):
            continue  # Spezies komplett ohne Required-Hits -> nicht ausgeben

        total = len(gids)
        with_all = sum(
            1
            for gid in gids
            if required_domains.issubset(genome_to_domains.get(gid, set()))
        )
        missing_any = total - with_all
        pct_complete = (with_all / total * 100.0) if total else 0.0
        rows.append((species, total, with_all, missing_any, pct_complete))

    # 4) Schreiben
    rows.sort(key=lambda r: (r[0] == unknown_label, r[0].casefold()))
    with open(outpath, "w", newline="") as w:
        w.write(
            "Species\tgenomes_total\tgenomes_with_all\tgenomes_missing_any\tpct_complete\tgene_cluster\n"
        )
        for species, total, with_all, missing_any, pct in rows:
            w.write(
                f"{species}\t{total}\t{with_all}\t{missing_any}\t{pct:.1f}\t{required_domains}\n"
            )

    return outpath


def _output_cluster_overview_by_required(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    required_domains: Set[str],
    *,
    writemode: str = "w",
) -> None:
    """
    Eine Zeile pro Gencluster, der mindestens eines der `required_domains` enthält.

    Spalten:
      genomeID, clusterID, cluster_types (space-separiert, nach Position sortiert), occurrences

    Sortierung der Zeilen:
      1) genomeID
      2) clusterID

    'cluster_types' enthält pro Protein den gewählten Typ (Domain-Name).
    'occurrences' zählt, wie oft exakt dieselbe Typen-Sequenz insgesamt vorkommt.
    """
    req = {d for d in (required_domains or set()) if d}

    # Cluster -> Proteine sammeln
    cluster_to_proteins: Dict[str, List[Any]] = defaultdict(list)
    for p in protein_dict.values():
        cid = getattr(p, "clusterID", None)
        if cid:
            cluster_to_proteins[cid].append(p)

    # Enthält Cluster mindestens einen required Typ?
    def cluster_contains_required(proteins: List[Any]) -> bool:
        if not req:
            return True
        for p in proteins:
            for dom in p.get_domains_dict().values():
                if getattr(dom, "domain", None) in req:
                    return True
        return False

    # Reihenfolge innerhalb des Clusters: entlang Genomkoordinate
    def sort_key_protein(p: Any) -> Tuple[str, int, str]:
        return (
            str(getattr(p, "gene_contig", "")),
            int(getattr(p, "gene_start", 0) or 0),
            str(getattr(p, "proteinID", "")),
        )

    # Typ eines Proteins bestimmen (bevorzugt required-Domain, sonst "beste" Domain)
    def protein_type(p: Any) -> str:
        doms = list(getattr(p, "get_domains_dict")().values())
        if not doms:
            return "UNK"
        preferred = [d for d in doms if getattr(d, "domain", None) in req]
        pool = preferred if preferred else doms

        def score_key(d) -> Tuple[float, int, int, str]:
            s = float(getattr(d, "score", 0) or 0.0)
            start = int(getattr(d, "domStart", 0) or 0)
            end = int(getattr(d, "domEnd", 0) or 0)
            length = end - start
            name = str(getattr(d, "domain", "") or "")
            # sortiert aufsteigend → Maxima negativieren, Länge negativieren
            return (-s, -length, start, name)

        best = sorted(pool, key=score_key)[0]
        return str(getattr(best, "domain", "") or "UNK")

    # Zeilen vorbereiten (mit Typen-Liste)
    rows_raw: List[Tuple[str, str, str]] = []
    for cid, plist in cluster_to_proteins.items():
        if not cluster_contains_required(plist):
            continue
        plist_sorted = sorted(plist, key=sort_key_protein)
        genome_id = getattr(plist_sorted[0], "genomeID", "") if plist_sorted else ""
        types_str = " ".join(protein_type(p) for p in plist_sorted)
        rows_raw.append((genome_id, cid, types_str))

    # Häufigkeiten identischer Typen-Sequenzen
    occurrences = Counter(r[2] for r in rows_raw)

    # Sortierung: genomeID, dann clusterID
    rows_raw.sort(key=lambda r: (r[0], r[1]))

    # Schreiben
    header = ["genomeID", "clusterID", "cluster_types", "occurrences"]
    with open(output_filepath, writemode) as writer:
        writer.write("\t".join(header) + "\n")
        for genome_id, cluster_id, types_str in rows_raw:
            writer.write(
                f"{genome_id}\t{cluster_id}\t{types_str}\t{occurrences[types_str]}\n"
            )


def _output_distinct_fasta_reports(
    directory: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    writemode: str = "w",
) -> Set[str]:
    """Writes all protein sequences into distinct FASTA files by domain class and for fusion domains.

    Args:
        directory (str): Output directory, e.g. 'output/'.
        protein_dict (Dict[str, Any]): Maps proteinID to Protein object.
        cluster_dict (Dict[str, Any]): Maps clusterID to Cluster object.
        writemode (str, optional): File mode, "w" for overwrite/new, "a" for append (default: "w").

    Returns:
        Set[str]: Set of all written file paths.

    Example:
        >>> files = _output_distinct_fasta_reports("output/", protein_dict, cluster_dict)
        >>> print(files)
        {'output/_PF00001.faa', 'output/_fused_domain_PF00001.faa'}
    """
    # Proteine sequences
    logger.info(f"Printing protein fasta files to {directory}")
    directory = os.path.join(directory, "Protein_sequences")  # protein fasta files
    os.mkdir(directory)

    proteine_type: Dict[str, list] = {}
    fusion_dict: Dict[str, Any] = {}
    files: Set[str] = set()

    # Group proteins by their domain type
    for proteinID, protein in protein_dict.items():
        domain_name = protein.get_domains()
        proteine_type.setdefault(domain_name, []).append(proteinID)
        if protein.get_domain_count() > 1:  # fusion protein detected
            fusion_dict[proteinID] = protein

    # Output: per domain type
    for domain_name, proteinID_list in proteine_type.items():
        filepath = os.path.join(
            directory, f"{domain_name}.faa"
        )  # it is important to just concatenate
        files.add(filepath)
        with open(filepath, writemode) as writer:
            for proteinID in proteinID_list:
                protein = protein_dict[proteinID]
                genome_id = protein.genomeID
                proteinlist = protein.get_protein_list()
                sequence = str(protein.protein_sequence).replace("*", "")

                cluster_id = protein.clusterID
                if cluster_id in cluster_dict:
                    cluster = cluster_dict[cluster_id]
                    cluster_list = cluster.get_cluster_list(",")
                    out = f">{' '.join(proteinlist[:-5])} {' '.join(cluster_list)}\n"
                    writer.write(out)
                    writer.write(sequence + "\n")
                else:
                    out = f">{' '.join(proteinlist[:-5])}\n"
                    writer.write(out)
                    writer.write(sequence + "\n")

    # Output: fused domains as separate FASTA
    for protein in fusion_dict.values():
        cluster_id = protein.clusterID
        domain_dict = protein.get_domains_dict()
        genome_id = protein.genomeID
        proteinlist = protein.get_protein_list()
        sequence = str(protein.protein_sequence).replace("*", "")
        for domain in domain_dict.values():
            domain_name = domain.domain
            domain_sequence = sequence[domain.start : domain.end]
            filepath = os.path.join(directory, f"multi_domain_{domain_name}.faa")
            files.add(filepath)
            with open(filepath, "a") as writer:
                if cluster_id in cluster_dict:
                    cluster = cluster_dict[cluster_id]
                    cluster_list = cluster.get_cluster_list(",")
                    out = f">{genome_id}-{' '.join(proteinlist[:-5])} {' '.join(cluster_list)}\n"
                    writer.write(out)
                    writer.write(domain_sequence + "\n")
                else:
                    out = f">{genome_id}-{' '.join(proteinlist[:-5])}\n"
                    writer.write(out)
                    writer.write(domain_sequence + "\n")
    # logger.info(f"FASTA output written to: {files}")
    _singletons(directory, files)  # Output singleton per genome and duplicates
    _clean_empty_files(directory)
    return files


def _singletons(directory, filepaths):
    """
    01.11.22
    Args:
        directory: Zielverzeichnis für neue Dateien
        filepaths: Liste der Eingabe-FASTA-Dateien

    Routine filtert alle Singletons:
      - Jede neue Datei erzeugt zwei Ausgaben:
        *_ortho.faa = nur Genome mit genau einem Protein
        *_paralog.faa = Genome mit >=2 Proteinen (Paraloge)
    """
    os.makedirs(directory, exist_ok=True)

    for filepath in filepaths:
        # FASTA parsen: header -> (description, sequence)
        records = {}
        current_header = None
        current_seq = []

        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_header:
                        records[current_header] = (
                            " ".join(current_header.split()),
                            "".join(current_seq),
                        )
                    current_header = line[1:]  # header ohne ">"
                    current_seq = []
                else:
                    current_seq.append(line)
            # letztes record hinzufügen
            if current_header:
                records[current_header] = (
                    " ".join(current_header.split()),
                    "".join(current_seq),
                )

        # GenomeIDs sammeln (angenommen: Teil vor erstem '-' im Header)
        genome_list = [header.split("-")[0] for header in records.keys()]
        duplicate = {gid for gid in genome_list if genome_list.count(gid) > 1}

        path = os.path.splitext(filepath)[0]
        name = os.path.basename(path)

        single_path = os.path.join(directory, f"{name}_ortho.faa")
        double_path = os.path.join(directory, f"{name}_paralog.faa")

        with open(single_path, "w") as single, open(double_path, "w") as double:
            for header, (desc, seq) in records.items():
                genomeID = header.split("-")[0]
                fasta_entry = f">{desc}\n{seq}\n"
                if genomeID not in duplicate:
                    single.write(fasta_entry)
                else:
                    double.write(fasta_entry)


def _output_cluster_context_report(
    output_filepath: str,
    context_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
) -> None:
    """
    Write an extended cluster context table.

    One row per protein in the provided context_dict.

    Columns:
        genomeID, clusterID, proteinID,
        contig, start, end, strand,
        domains,
        Superkingdom, Phylum, Class, Order, Family, Genus, Species
    """
    header_cols = [
        "genomeID",
        "clusterID",
        "proteinID",
        "contig",
        "start",
        "end",
        "strand",
        "domains",
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    # sort proteins for readability
    def sort_key(p):
        return (
            getattr(p, "genomeID", ""),
            getattr(p, "clusterID", ""),
            getattr(p, "gene_contig", ""),
            getattr(p, "gene_start", 0),
        )

    proteins_sorted = sorted(context_dict.values(), key=sort_key)

    with open(output_filepath, "w") as fh:
        fh.write("\t".join(header_cols) + "\n")

        for p in proteins_sorted:
            gid = getattr(p, "genomeID", "")
            rec = taxon_dict.get(gid, {}) or {}

            # domains as comma-separated list
            domains = ""
            doms = getattr(p, "domains", None)
            if doms:
                # domains is typically a list of tuples: (domain, start, end, score)
                try:
                    domains = ",".join(str(d[0]) for d in doms if d)
                except TypeError:
                    print(doms)
            row = [
                str(gid),
                str(getattr(p, "clusterID", "")),
                str(getattr(p, "proteinID", "")),
                str(getattr(p, "gene_contig", "")),
                str(getattr(p, "gene_start", "")),
                str(getattr(p, "gene_end", "")),
                str(getattr(p, "gene_strand", "")),
                domains,
                str(rec.get("Superkingdom", "")),
                str(rec.get("Phylum", "")),
                str(rec.get("Class", "")),
                str(rec.get("Order", "")),
                str(rec.get("Family", "")),
                str(rec.get("Genus", "")),
                str(rec.get("Species", "")),
            ]

            fh.write("\t".join(row) + "\n")

    # logger.info("Wrote cluster context report: %s", output_filepath)


def _clean_empty_files(directory: str) -> None:
    """Removes all empty files in a given directory."""
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            os.remove(file_path)
            logger.debug(f"Removed empty file: {file_path}")


def _output_domain_function_report(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    domain_annotations: Dict[str, Dict[str, str]],
    writemode: str = "w",
) -> None:
    """
    Schreibt pro Protein/Domain Funktionszeilen:
      genomeID, proteinID, domain, reaction, protein_description, system, metabolism,
      Superkingdom, Phylum, Class, Order, Family, Genus, Species

    Multi-Domain-Proteine:
      - Wenn (reaction, protein_description, system, metabolism) für alle gefundenen Domains identisch sind,
        schreibe nur EINE Zeile (domain = alle Domainnamen mit '-').
      - Sonst schreibe eine Zeile pro Domain (nur für Domains, die in domain_annotations vorkommen).
    """
    tax_cols = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    header = "\t".join(
        [
            "genomeID",
            "proteinID",
            "domain",
            "reaction",
            "protein_description",
            "system",
            "metabolism",
            *tax_cols,
        ]
    )

    # sortiere analog zu deinem Hauptreport (genomeID, contig, start)
    proteinID_list = sorted(
        protein_dict,
        key=lambda x: (
            protein_dict[x].genomeID,
            protein_dict[x].gene_contig,
            protein_dict[x].gene_start,
        ),
    )

    with open(output_filepath, writemode) as w:
        w.write(header + "\n")

        for proteinID in proteinID_list:
            p = protein_dict[proteinID]
            gid = p.genomeID
            rec = taxon_dict.get(gid, {}) or {}
            tax_vals = [rec.get(k) for k in tax_cols]

            # Domains extrahieren: bevorzugt String "A-B-C" aus get_domains(), sonst aus get_domains_dict()
            dom_string = ""
            try:
                dom_string = p.get_domains() or ""
            except Exception:
                dom_string = ""

            if dom_string:
                dom_list = [d.strip() for d in dom_string.split("-") if d.strip()]
            else:
                dct = getattr(p, "get_domains_dict")()
                dom_list = [
                    getattr(d, "domain", "").strip()
                    for d in dct.values()
                    if getattr(d, "domain", "").strip()
                ]

            # Nur Domains berücksichtigen, die in den Annotationen vorkommen
            dom_list = [d for d in dom_list if d in domain_annotations]
            if not dom_list:
                continue  # nichts zu schreiben

            # Funktions-Tupel je Domain einsammeln
            func_tuples: List[
                Tuple[str, str, str, str, str]
            ] = []  # (domain, reaction, desc, system, metabolism)
            for d in dom_list:
                info = domain_annotations.get(d, {})
                func_tuples.append(
                    (
                        d,
                        info.get("reaction", ""),
                        info.get("protein_description", ""),
                        info.get("system", ""),
                        info.get("metabolism", ""),
                    )
                )

            # Prüfen, ob alle Funktionswerte (ohne Domain) identisch sind
            unique_payloads = {
                (ft[1], ft[2], ft[3], ft[4])  # reaction, desc, system, metabolism
                for ft in func_tuples
            }

            if len(unique_payloads) == 1:
                # eine Zeile, Domainnamen zusammengeführt
                reaction, desc, system, metab = next(iter(unique_payloads))
                merged_domain = "-".join(dom_list)
                row = [
                    gid,
                    p.proteinID if hasattr(p, "proteinID") else proteinID,
                    merged_domain,
                    reaction,
                    desc,
                    system,
                    metab,
                    *tax_vals,
                ]
                w.write("\t".join(map(str, row)) + "\n")
            else:
                # mehrere Zeilen, je Domain separat
                for d, reaction, desc, system, metab in func_tuples:
                    row = [
                        gid,
                        p.proteinID if hasattr(p, "proteinID") else proteinID,
                        d,
                        reaction,
                        desc,
                        system,
                        metab,
                        *tax_vals,
                    ]
                    w.write("\t".join(map(str, row)) + "\n")


def print_hit_reports(
    directory: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, Dict[str, str]],
    metabolic_dict: Dict[str, Any],
    context_dict: Dict[str, Any] | None,
    fetch_proteins: List[str],
) -> None:
    """
    Main output routine: creates hit tables, taxonomy summaries, and protein FASTA files.

    Args:
        metabolic_dict:
        fetch_proteins: input proteins from the CLI
        directory (str): Output directory path.
        protein_dict: Mapping proteinID -> proteinObj.
        cluster_dict: Mapping clusterID -> clusterObj.
        taxon_dict (Dict[str, str]): Mapping genomeID -> taxonomy string.
    """

    # Printing results
    logger.info("Printing genome information output files")
    # Metadata for taxonomy and hits
    hit_report = os.path.join(
        directory, "summary_hit_table.txt"
    )  # individual hit table in tsv file
    gene_taxonomy = os.path.join(directory, "summary_gene_taxonomy.txt")
    unique_file = os.path.join(directory, "summary_unique_lineages.txt")
    taxonomy_summary = os.path.join(directory, "summary_hit_taxonomy_counts.txt")
    taxonomy_summary2 = os.path.join(
        directory, "summary_requested_hit_taxonomy_counts.txt"
    )
    metabolic_annotation = os.path.join(directory, "summary_metabolic_annotations.txt")
    cluster_overview_report = os.path.join(
        directory, "summary_genecluster_overview_table.txt"
    )

    # Output hit report
    _output_genome_report(hit_report, protein_dict, cluster_dict, taxon_dict)

    # Output gene taxonomy report
    _output_protein_taxonomy(gene_taxonomy, protein_dict, taxon_dict)

    # Output unique taxonomy report
    _output_unique_taxonomy_table(unique_file, protein_dict, taxon_dict)

    # Output domain annotation for metabolism
    _output_domain_function_report(
        output_filepath=metabolic_annotation,
        protein_dict=protein_dict,
        taxon_dict=taxon_dict,
        domain_annotations=metabolic_dict,
        writemode="w",
    )
    # Output taxonomy summary
    _output_taxonomy_summary(taxonomy_summary, protein_dict, taxon_dict)

    _output_taxonomy_summary(
        taxonomy_summary2, protein_dict, taxon_dict, fetch_proteins
    )

    # Output strain variability summary
    _output_strain_variability_by_species(
        directory, protein_dict, taxon_dict, set(fetch_proteins)
    )

    # Output for each protein the genomic context
    _output_cluster_overview_by_required(
        cluster_overview_report, context_dict, set(fetch_proteins)
    )


def print_fasta_files(directory, protein_dict, cluster_dict):
    logger.info("Writing fasta formated output files to disk")
    # Output fasta files for hits
    _output_distinct_fasta_reports(directory, protein_dict, cluster_dict)

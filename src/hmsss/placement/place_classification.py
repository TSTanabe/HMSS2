#!/usr/bin/env python3
"""Classify JPLACE queries to protein types and conservative reference taxonomy.

The reference tree is a PROTEIN tree. Protein types are loaded from sequence annotations, while organism
taxonomy is joined separately using the genome identifier extracted from each sequence ID. Only monophyletic protein types are accepted.
The algorithm is a demonstration, not an exact reimplementation of GraftM.
"""

import argparse
import csv
import json
from pathlib import Path


class Node:
    def __init__(self, name=None, edge=None, children=None):
        self.name, self.edge, self.children = name, edge, children or []


def parse_tree(s):
    i = 0

    def ws():
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def parse():
        nonlocal i
        ws()
        children = []
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                children.append(parse())
                ws()
                if i < len(s) and s[i] == ",":
                    i += 1
                elif i < len(s) and s[i] == ")":
                    i += 1
                    break
                else:
                    raise ValueError("Malformed Newick")
        ws()
        begin = i
        while i < len(s) and s[i] not in ":,);{":
            i += 1
        name = s[begin:i].strip() or None
        if i < len(s) and s[i] == ":":
            i += 1
            while i < len(s) and s[i] not in "{,);":
                i += 1
        edge = None
        if i < len(s) and s[i] == "{":
            i += 1
            start = i
            while i < len(s) and s[i] != "}":
                i += 1
            if i == len(s):
                raise ValueError("Unterminated edge number")
            edge = int(s[start:i])
            i += 1
        return Node(name, edge, children)

    root = parse()
    ws()
    if i >= len(s) or s[i] != ";" or s[i + 1 :].strip():
        raise ValueError("Malformed Newick terminator")
    return root


def common_prefix(paths):
    shared = []
    for rank_items in zip(*paths):
        if len(set(rank_items)) != 1:
            break
        shared.append(rank_items[0])
    return shared


def load_annotations(filename, taxonomy_file):
    """Join protein annotations with genome taxonomy using the genome ID."""

    # 1. Load taxonomy (two columns, optional header)
    taxonomy = {}

    with open(taxonomy_file, newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if not row or all(not value.strip() for value in row):
                continue

            if len(row) < 2:
                raise ValueError(f"Invalid taxonomy row at line {line_number}: {row}")

            genome_id = row[0].strip()
            tax_string = row[1].strip()

            # Skip optional header
            if not taxonomy and genome_id.lower() in {
                "genomeid",
                "genome_id",
                "genome",
            }:
                continue

            if not genome_id or genome_id in taxonomy:
                raise ValueError(f"Empty/duplicate genome ID: {genome_id}")

            # Greengenes/GTDB taxonomy; remove consecutive duplicated ranks.
            raw_tax = [part.strip() for part in tax_string.split(";") if part.strip()]
            tax = []
            for rank in raw_tax:
                if not tax or tax[-1] != rank:
                    tax.append(rank)

            if not tax:
                raise ValueError(f"Missing taxonomy for genome: {genome_id}")

            taxonomy[genome_id] = tax

    # 2. Load protein annotations and join taxonomy per reference ID.
    data = {}
    missing = []
    matched_exact = 0
    matched_prefix = 0

    with open(filename, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = set(reader.fieldnames or [])
        reference_column = next(
            (c for c in ("sequence_id", "reference_id") if c in fields), None
        )
        if reference_column is None or "protein_type" not in fields:
            raise ValueError("Annotations TSV requires sequence_id and protein_type")

        for row in reader:
            ref = row[reference_column].strip()
            if not ref or ref in data:
                raise ValueError(f"Empty/duplicate reference ID: {ref}")
            protein_type = row["protein_type"].strip()
            if not protein_type:
                raise ValueError(f"Missing protein type for reference: {ref}")

            # Genome ID is the segment before the FIRST '___'.
            # Some reference headers add an alias before '-' (e.g.
            # GCA_031264015-GCA_031264015___2278).
            full_genome_id = ref.split("___", 1)[0]
            if full_genome_id in taxonomy:
                genome_id = full_genome_id
                matched_exact += 1
            else:
                short_genome_id = full_genome_id.split("-", 1)[0]
                genome_id = short_genome_id
                if genome_id in taxonomy:
                    matched_prefix += 1

            if genome_id not in taxonomy:
                missing.append((ref, full_genome_id, genome_id))
                # Missing taxonomies must not be silently classified as
                # another genome. Empty path contributes no taxonomic support.
                data[ref] = (protein_type, [])
            else:
                data[ref] = (protein_type, taxonomy[genome_id])

    print(f"Loaded {len(taxonomy):,} genome taxonomies")
    print(f"Loaded {len(data):,} protein annotations")
    print(
        f"Taxonomy matches: exact={matched_exact:,}, prefix={matched_prefix:,}, missing={len(missing):,}"
    )
    if missing:
        import sys

        for ref, full_id, attempted in missing[:10]:
            print(
                f"WARNING: no genome taxonomy for {ref} (tried: {full_id}, {attempted})",
                file=sys.stderr,
            )

    return data


def annotate_edges(root, annotations):
    edges, clade_sets, leaves = {}, [], set()

    def visit(node):
        if node.children:
            descendants = set().union(*(visit(child) for child in node.children))
        else:
            if not node.name or node.name in leaves:
                raise ValueError(f"Duplicate/empty leaf name: {node.name}")
            if node.name not in annotations:
                raise ValueError(f"Missing annotation: {node.name}")
            leaves.add(node.name)
            descendants = {node.name}
        if node.edge is not None:
            if node.edge in edges:
                raise ValueError(f"Duplicate edge ID {node.edge}")
            kinds = {annotations[leaf][0] for leaf in descendants}
            edges[node.edge] = {
                "protein_type": next(iter(kinds)) if len(kinds) == 1 else None,
                "taxonomy": common_prefix(
                    [annotations[leaf][1] for leaf in sorted(descendants)]
                ),
                "descendants": descendants,
            }
        clade_sets.append(descendants)
        return descendants

    visit(root)
    if leaves != set(annotations):
        raise ValueError(f"Tree/annotation mismatch: {set(annotations) ^ leaves}")

    return edges


def classify(jplace, annotations_file, taxonomy_file, threshold=0.75):
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    with open(jplace) as f:
        data = json.load(f)
    edge_info = annotate_edges(
        parse_tree(data["tree"]), load_annotations(annotations_file, taxonomy_file)
    )
    edge_ix = data["fields"].index("edge_num")
    lwr_ix = data["fields"].index("like_weight_ratio")
    result = {}
    for block in data["placements"]:
        weights = []
        for row in block["p"]:
            edge = int(row[edge_ix])
            lwr = float(row[lwr_ix])
            if edge not in edge_info or not 0 <= lwr <= 1:
                raise ValueError(f"Invalid edge or LWR: {edge}, {lwr}")
            weights.append((edge_info[edge], lwr))
        total = sum(w for _, w in weights)
        if total > 1.000001:
            raise ValueError("LWR sum above one")
        # Protein: aggregated LWR on edges entirely within one type.
        support_by_type = {}
        for info, w in weights:
            typ = info["protein_type"]
            if typ is not None:
                support_by_type[typ] = support_by_type.get(typ, 0) + w
        best_type, type_support = (
            max(support_by_type.items(), key=lambda kv: (kv[1], kv[0]))
            if support_by_type
            else (None, 0.0)
        )
        assigned_type = best_type if type_support >= threshold else None

        # Taxonomy: only use the whole placement distribution, not just the
        # winning placement; stop at the first unsupported rank.
        path = []
        rank_support = []
        if assigned_type:
            maxdepth = max((len(info["taxonomy"]) for info, _ in weights), default=0)
            for level in range(maxdepth):
                support = {}
                for info, weight in weights:
                    tax = info["taxonomy"]
                    if (
                        len(tax) > level
                        and tax[:level] == path
                        and tax[level].lower() not in {"unknown", "unclassified"}
                    ):
                        label = tax[level]
                        support[label] = support.get(label, 0) + weight
                if not support:
                    break
                label, w = max(support.items(), key=lambda kv: (kv[1], kv[0]))
                if w < threshold:
                    break
                path.append(label)
                rank_support.append(w)

        names = block.get("n")
        if names is None:
            names = [row[0] for row in block.get("nm", [])]
        if not names:
            raise ValueError("Placement has no query identifiers")
        for name in names:
            if name in result:
                raise ValueError(f"Duplicate query: {name}")
            result[name] = {
                "protein_type": assigned_type or "Unclassified",
                "type_lwr": round(type_support, 6),
                "taxonomy": ";".join(path) or "Unclassified",
                "taxonomic_depth": len(path),
                "last_rank_lwr": round(rank_support[-1], 6) if rank_support else 0.0,
                "stored_lwr": round(total, 6),
            }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--jplace", type=Path, default=Path(__file__).with_name("example.jplace")
    )
    p.add_argument(
        "--annotations",
        type=Path,
        default=Path(__file__).with_name("reference_annotations.tsv"),
    )
    p.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="TSV with genomeID and taxonomy columns",
    )
    p.add_argument("--min-lwr", type=float, default=0.63)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = classify(args.jplace, args.annotations, args.taxonomy, args.min_lwr)
    columns = [
        "read_id",
        "protein_type",
        "type_lwr",
        "taxonomy",
        "taxonomic_depth",
        "last_rank_lwr",
        "stored_lwr",
    ]
    rows = [{"read_id": name, **r} for name, r in sorted(result.items())]
    if args.output:
        with args.output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    print("\t".join(columns))
    for row in rows:
        print("\t".join(map(str, (row[c] for c in columns))))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
pplacer_tax_scampp_graftm.py

A SCAMPP-style wrapper around pplacer that mimics the *GraftM* pplacer() interface:
  - input_path: combined alignment FASTA (reference taxa + query sequences already aligned)
  - output_path + output_file: where to write <output_file>.jplace
  - threads: forwarded to pplacer (-j)

This script runs:
  1) subtree selection per query (Hamming nearest neighbor + subtree of size n)
  2) taxit create (build a temporary refpkg for the subtree)
  3) pplacer on the subtree
  4) remap placements back to the *backbone* tree edges
  5) write a final backbone-referenced .jplace file (as GraftM expects)

Refpkg layout assumption (as you specified):
  - *.tree
  - *.tree_log
  - *.json
  - deduplicated_aligned.fasta
  - *.seqinfo.csv
  - taxonomy.csv

We resolve needed files from the provided --refpkg directory:
  - tree_file      : first "*.tree"
  - ref_alignment  : "deduplicated_aligned.fasta"
  - tree_stats/info: prefer "*.json" containing "tree_stats"/"stats" keys; else fallback to "*.tree_log"; else any "*.json"
  - taxonomy/seqinfo are not required for placement, but may exist in the package.

Executables:
  - pplacer and taxit are expected to be available on PATH (e.g., from a conda environment).
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import treeswift

# -----------------------------------------------------------------------------
# Utils (adapted from the user's uploaded utils.py; corrected/trimmed for use)
# Source reference: /mnt/data/utils.py
# -----------------------------------------------------------------------------
# NOTE: The original uploaded utils.py contains multiple syntactic/semantic issues
# (e.g., duplicate function names, malformed string formatting). The functions
# below are a faithful *functional* adaptation of the intended behavior.
# (You asked to include the needed utils from the uploaded file.)
# -----------------------------------------------------------------------------

BRACKET = {"[": "]", "{": "}", "'": "'", '"': '"'}


def read_fasta_to_dict(path: str) -> Dict[str, str]:
    """Read FASTA into dict[label]=sequence (keeps alignment characters)."""
    result: Dict[str, str] = {}
    label: Optional[str] = None
    seq_chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if label is not None:
                    result[label] = "".join(seq_chunks)
                label = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if label is not None:
            result[label] = "".join(seq_chunks)
    return result


def separate_ref_and_query(aln_dict: Dict[str, str], backbone_leaf_labels: set[str]) -> Tuple[
    Dict[str, str], Dict[str, str]]:
    """Split sequences into ref (labels in backbone) and queries (labels not in backbone)."""
    ref: Dict[str, str] = {}
    query: Dict[str, str] = {}
    for k, v in aln_dict.items():
        if k in backbone_leaf_labels:
            ref[k] = v
        else:
            query[k] = v
    return ref, query


def hamming(seq1: str, seq2: str) -> int:
    return sum(1 for ch1, ch2 in zip(seq1, seq2) if ch1 != ch2)


def set_fragment_indices(x: str) -> Tuple[int, int]:
    """
    Return (si, ei) indices excluding leading/trailing gaps '-' (fragment handling).
    """
    e = len(x)
    si = 0
    ei = e
    # leading
    while si < ei and x[si] == "-":
        si += 1
    # trailing
    while ei > si and x[ei - 1] == "-":
        ei -= 1
    return si, ei


def find_closest_hamming(x: str, ref: Dict[str, str], n: int, fragment_flag: bool) -> List[str]:
    """
    Return list of n closest reference labels by (fragment-aware) Hamming distance.
    Uses a chunked priority queue to avoid computing full Hamming for all refs up front.
    """
    if n <= 0:
        return []

    if fragment_flag:
        si, ei = set_fragment_indices(x)
    else:
        si, ei = 0, len(x)

    # chunk size
    c = 200
    queue: List[Tuple[int, int, int, str]] = []
    counter = 0

    for name, seq in ref.items():
        # initial chunk
        first = hamming(seq[si: si + c], x[si: si + c])
        sites_left = (ei - si) - c
        heapq.heappush(queue, (first, sites_left, counter, name))
        counter += 1

    closest: List[str] = []
    while queue:
        ham_dist, sites_left, _, name = heapq.heappop(queue)
        if sites_left < 0:
            closest.append(name)
            if len(closest) >= n:
                return closest
            continue

        # advance chunk window
        ind = ei - sites_left
        new_ham = hamming(ref[name][ind: ind + c], x[ind: ind + c])
        heapq.heappush(queue, (ham_dist + new_ham, sites_left - c, counter, name))
        counter += 1

    return closest


def subtree_nodes(tree: treeswift.Tree, leaf_y: treeswift.Node, n: int) -> List[str]:
    """Topological BFS-style subtree selection around a seed leaf (n leaves)."""
    queue: List[Tuple[int, int, treeswift.Node]] = [(0, 0, leaf_y.get_parent())]
    leaves = [leaf_y]
    visited = {leaf_y}
    counter = 1

    while len(leaves) < n and queue:
        length, _, node = heapq.heappop(queue)
        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        adjacent = list(node.child_nodes())
        if not node.is_root():
            adjacent.append(node.get_parent())

        for neighbor in adjacent:
            if neighbor not in visited:
                heapq.heappush(queue, (length + 1, counter, neighbor))
                counter += 1

    return [x.get_label() for x in leaves]


def subtree_nodes_with_edge_length(tree: treeswift.Tree, leaf_y: treeswift.Node, n: int) -> List[str]:
    """Edge-length weighted subtree selection around a seed leaf (n leaves)."""
    queue: List[Tuple[float, treeswift.Node]] = [(leaf_y.get_edge_length() or 0.0, leaf_y.get_parent())]
    leaves = [leaf_y]
    visited = {leaf_y}

    while len(leaves) < n and queue:
        length, node = heapq.heappop(queue)
        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        adjacent = list(node.child_nodes())
        if not node.is_root():
            adjacent.append(node.get_parent())

        for neighbor in adjacent:
            if neighbor in visited:
                continue
            if neighbor == node.get_parent():
                heapq.heappush(queue, (length + (node.get_edge_length() or 0.0), neighbor))
            else:
                heapq.heappush(queue, (length + (neighbor.get_edge_length() or 0.0), neighbor))

    return [x.get_label() for x in leaves]


def add_edge_numbers(tree: treeswift.Tree) -> None:
    """
    Add "%%<id>" suffix to each node label (postorder numbering).
    Root also gets a token (matches the intended behavior in the uploaded utils.py).
    """
    counter = 0
    for node in tree.traverse_postorder():
        counter += 1
        label = node.get_label()
        if label is None:
            node.set_label(f"%%{counter}")
        else:
            node.set_label(f"{label}%%{counter}")


def remove_edge_numbers(tree: treeswift.Tree) -> None:
    """Strip trailing '%%<id>' from node labels."""
    for node in tree.traverse_postorder():
        lab = node.get_label()
        if not lab:
            continue
        parts = lab.split("%%", 1)
        node.set_label(parts[0] if parts[0] else None)


def newick_with_edge_tokens(tree: treeswift.Tree) -> str:
    """
    Produce a jplace-style tree string with edge tokens: ':<len>{edge_num}'.
    We use the '%%' numbers attached to child node labels as edge_num.
    """

    def node_to_str(n: treeswift.Node) -> str:
        lab = n.get_label() or ""
        if "%%" in lab:
            label_part, edge_id = lab.split("%%", 1)
        else:
            label_part, edge_id = lab, None

        if n.is_leaf():
            return str(label_part) if label_part else ""
        else:
            children = []
            for c in n.child_nodes():
                c_str = node_to_str(c)
                # edge token comes from child
                c_lab = c.get_label() or ""
                if "%%" in c_lab:
                    _, c_edge_id = c_lab.split("%%", 1)
                else:
                    c_edge_id = None
                el = c.get_edge_length()
                if el is None:
                    children.append(c_str)
                else:
                    children.append(f"{c_str}:{el}{{{int(c_edge_id) if c_edge_id is not None else 0}}}")
            inside = ",".join(children)
            return f"({inside}){label_part}" if label_part else f"({inside})"

    # treeswift doesn't store root edge in Newick meaningfully; jplace wants tokens on edges to children.
    return node_to_str(tree.root) + ";"


def read_tree_newick_edge_tokens(newick_str: str) -> Tuple[treeswift.Tree, Dict[str, treeswift.Node]]:
    """
    Parse a Newick tree that has edge tokens in '{...}' after edge lengths.
    Returns (tree, edge_dict) mapping edge_num(str)->node that is child of that edge.
    """
    edge_dict: Dict[str, treeswift.Node] = {}

    ts = newick_str.strip()
    if ts.startswith("[&R]"):
        ts = ts.split("]", 1)[1].strip()

    # Minimal parser (adapted from uploaded utils.py intent).
    t = treeswift.Tree()
    n = t.root
    i = 0

    while i < len(ts):
        ch = ts[i]

        if ch == ";":
            break
        elif ch == "(":
            c = treeswift.Node()
            n.add_child(c)
            n = c
        elif ch == ")":
            n = n.parent
        elif ch == ",":
            n = n.parent
            c = treeswift.Node()
            n.add_child(c)
            n = c
        elif ch == ":":
            i += 1
            ls = ""
            while i < len(ts) and ts[i] not in ",);{":
                ls += ts[i]
                i += 1
            n.edge_length = float(ls) if ls else None
            i -= 1
        elif ch == "{":
            i += 1
            token = ""
            while i < len(ts) and ts[i] != "}":
                token += ts[i]
                i += 1
            edge_dict[token] = n
        else:
            # label
            label = ""
            bracket = None
            while i < len(ts) and (
                    bracket is not None
                    or ts[i] in BRACKET
                    or ts[i] not in ":,;){"
            ):
                if ts[i] in BRACKET and bracket is None:
                    bracket = ts[i]
                elif bracket is not None and ts[i] == BRACKET[bracket]:
                    bracket = None
                label += ts[i]
                i += 1
            n.label = label if label else None
            i -= 1

        i += 1

    return t, edge_dict


def find_closest_leaf_or_path(x: treeswift.Node, visited: set[treeswift.Node], y: Optional[treeswift.Node] = None):
    """
    If y is None:
      return (closest_leaf, path_nodes_from_x_to_leaf_exclusive_of_x)
    If y is provided:
      return (y, path_nodes_from_x_towards_y) avoiding visited
    Adapted conceptually from uploaded utils.py.
    """
    queue: List[Tuple[float, int, List[treeswift.Node], treeswift.Node]] = []
    cnt = 1
    visited.add(x)

    if x.get_parent() and x.get_parent() not in visited:
        heapq.heappush(queue, (x.get_edge_length() or 0.0, cnt, [x], x.get_parent()))
        cnt += 1

    for child in x.child_nodes():
        if child not in visited:
            heapq.heappush(queue, (child.get_edge_length() or 0.0, cnt, [child], child))
            cnt += 1

    while queue:
        length, _, path, node = heapq.heappop(queue)
        visited.add(node)

        if node.is_leaf():
            if y is None or (node.get_label() == y.get_label()):
                return node, path
            continue

        if node.get_parent() and node.get_parent() not in visited:
            tmp = path.copy()
            tmp.append(node)
            heapq.heappush(queue, (length + (node.get_edge_length() or 0.0), cnt, tmp, node.get_parent()))
            cnt += 1

        for child in node.child_nodes():
            if child not in visited:
                tmp = path.copy()
                tmp.append(node)
                heapq.heappush(queue, (length + (child.get_edge_length() or 0.0), cnt, tmp, child))
                cnt += 1

    return x, [x]


# -----------------------------------------------------------------------------
# Refpkg resolution
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class RefPkgFiles:
    tree_file: str
    ref_alignment: str
    tree_info: str  # tree_stats/info used for `taxit create --tree-stats`
    json_file: Optional[str] = None
    taxonomy_csv: Optional[str] = None
    seqinfo_csv: Optional[str] = None
    tree_log: Optional[str] = None


def resolve_refpkg_files(refpkg_dir: str) -> RefPkgFiles:
    rp = Path(refpkg_dir)
    if not rp.exists() or not rp.is_dir():
        raise FileNotFoundError(f"refpkg directory not found: {refpkg_dir}")

    tree_candidates = sorted(rp.glob("*.tree"))
    if not tree_candidates:
        raise RuntimeError(f"No '*.tree' found in refpkg: {refpkg_dir}")
    tree_file = str(tree_candidates[0])

    ref_alignment = rp / "deduplicated_aligned.fasta"
    if not ref_alignment.exists():
        raise RuntimeError(f"Expected 'deduplicated_aligned.fasta' in refpkg: {refpkg_dir}")

    json_candidates = sorted(rp.glob("*.json"))
    tree_log_candidates = sorted(rp.glob("*.tree_log"))

    taxonomy_csv = rp / "taxonomy.csv"
    seqinfo_csv = next(iter(rp.glob("*.seqinfo.csv")), None)

    # Choose tree_info:
    #  1) any *.json containing keys like "tree_stats"/"stats"
    #  2) else *.tree_log
    #  3) else first *.json
    chosen_json = None
    tree_info = None

    for j in json_candidates:
        try:
            with open(j, "r", encoding="utf-8") as fh:
                obj = json.load(fh)
            # heuristic key check
            if isinstance(obj, dict) and any(k in obj for k in ("tree_stats", "stats", "tree_stats_file", "treeStats")):
                chosen_json = str(j)
                tree_info = str(j)
                break
        except Exception:
            continue

    if tree_info is None and tree_log_candidates:
        tree_info = str(tree_log_candidates[0])

    if tree_info is None and json_candidates:
        chosen_json = str(json_candidates[0])
        tree_info = str(json_candidates[0])

    if tree_info is None:
        raise RuntimeError(f"Could not determine a tree-stats/info file in refpkg: {refpkg_dir}")

    return RefPkgFiles(
        tree_file=tree_file,
        ref_alignment=str(ref_alignment),
        tree_info=tree_info,
        json_file=chosen_json,
        taxonomy_csv=str(taxonomy_csv) if taxonomy_csv.exists() else None,
        seqinfo_csv=str(seqinfo_csv) if seqinfo_csv else None,
        tree_log=str(tree_log_candidates[0]) if tree_log_candidates else None,
    )


# -----------------------------------------------------------------------------
# External tool runners (conda env PATH)
# -----------------------------------------------------------------------------

def run_cmd(cmd: List[str], *, cwd: Optional[str] = None) -> None:
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  cwd: {cwd}\n"
            f"  exit: {proc.returncode}\n"
            f"  stdout:\n{proc.stdout}\n"
            f"  stderr:\n{proc.stderr}\n"
        )


def ensure_on_path(exe: str) -> str:
    path = shutil.which(exe)
    if not path:
        raise RuntimeError(f"Required executable not found on PATH: {exe}. Activate the conda env providing it.")
    return path


# -----------------------------------------------------------------------------
# Main SCAMPP placement logic
# -----------------------------------------------------------------------------

def write_fasta(path: str, records: Dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for k, v in records.items():
            fh.write(f">{k}\n")
            # wrap
            for i in range(0, len(v), 80):
                fh.write(v[i: i + 80] + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run pplacer-tax-SCAMPP with GraftM-like pplacer() inputs and write a GraftM-compatible .jplace."
    )
    # Inputs aligned with the GraftM pplacer() call signature conceptually:
    ap.add_argument("--output_file", required=True,
                    help="Basename (without extension) for final .jplace (like GraftM output_file).")
    ap.add_argument("--output_path", required=True,
                    help="Directory to write the final .jplace (like GraftM output_path).")
    ap.add_argument("--input_path", required=True,
                    help="Combined alignment FASTA containing backbone taxa + query sequences (like GraftM input_path).")
    ap.add_argument("--threads", type=int, default=1, help="Threads passed to pplacer (-j).")

    # Additional required input that GraftM typically has as self.refpkg:
    ap.add_argument("--refpkg", required=True,
                    help="Path to the GraftM refpkg directory (contains *.tree, deduplicated_aligned.fasta, etc.).")

    # SCAMPP parameters (defaults follow the SCAMPP scripts' typical defaults)
    ap.add_argument("--subtreesize", type=int, default=2000, help="Number of leaves in the selected subtree.")
    ap.add_argument(
        "--subtreetype",
        choices=["d", "n", "h"],
        default="d",
        help="Subtree selection type: d=edge-length weighted, n=topological, h=take n nearest by Hamming directly.",
    )
    ap.add_argument("--fragmentflag", action="store_true",
                    help="Treat queries as fragments (ignore leading/trailing gaps).")
    ap.add_argument("--model", default="GTR", help="pplacer model string (passed to pplacer -m).")

    args = ap.parse_args()

    ensure_on_path("pplacer")
    ensure_on_path("taxit")

    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    final_jplace = output_path / f"{args.output_file}.jplace"

    refpkg_files = resolve_refpkg_files(args.refpkg)

    # Load backbone tree
    backbone_tree = treeswift.read_tree_newick(refpkg_files.tree_file)
    # Build backbone leaf-label set and leaf_dict for mapping
    backbone_leaf_labels = set(n.get_label() for n in backbone_tree.traverse_leaves() if n.get_label() is not None)
    leaf_dict = {n.get_label(): n for n in backbone_tree.traverse_leaves() if n.get_label() is not None}

    # Read combined alignment and split into ref vs query (based on labels present in backbone tree)
    aln_dict = read_fasta_to_dict(args.input_path)
    ref_dict, q_dict = separate_ref_and_query(aln_dict, backbone_leaf_labels)

    # Prepare final backbone jplace scaffold
    add_edge_numbers(backbone_tree)
    jplace = {
        "tree": newick_with_edge_tokens(backbone_tree),
        "placements": [],
        "metadata": {
            "invocation": " ".join(sys.argv),
            "refpkg": str(Path(args.refpkg).resolve()),
        },
        "version": 3,
        "fields": ["distal_length", "edge_num", "like_weight_ratio", "likelihood", "pendant_length"],
    }

    # For each query: build subtree, run pplacer, remap placements
    with tempfile.TemporaryDirectory(prefix="tax_scampp_") as tmpd:
        tmpd_p = Path(tmpd)

        # We'll need a backbone tree without edge numbers for subtree extraction
        backbone_plain = treeswift.read_tree_newick(refpkg_files.tree_file)

        for qi, (q_name, q_seq) in enumerate(q_dict.items(), start=1):
            # 1) choose closest refs by (fragment-aware) Hamming
            if args.subtreetype == "h":
                closest = find_closest_hamming(q_seq, ref_dict, args.subtreesize, args.fragmentflag)
                if not closest:
                    continue
                labels = closest
            else:
                # nearest neighbor only
                nearest = find_closest_hamming(q_seq, ref_dict, 1, args.fragmentflag)
                if not nearest:
                    continue
                seed_label = nearest[0]

                # locate seed node in the backbone_plain tree
                seed_node = backbone_plain.label_to_node(selection="leaves").get(seed_label)
                if seed_node is None:
                    # Should not happen if labels match, but stay safe
                    continue

                if args.subtreetype == "n":
                    labels = subtree_nodes(backbone_plain, seed_node, args.subtreesize)
                else:
                    labels = subtree_nodes_with_edge_length(backbone_plain, seed_node, args.subtreesize)

            labels = [lab for lab in labels if lab in ref_dict]
            if not labels:
                continue

            # 2) build temporary subtree tree
            subtree = backbone_plain.extract_tree_with(labels)
            subtree.resolve_polytomies()

            tmp_tree = tmpd_p / f"subtree_{qi}.nwk"
            subtree.write_tree_newick(str(tmp_tree))

            # 3) build tmp alignment: query + subtree refs
            tmp_aln = tmpd_p / f"aln_{qi}.fasta"
            tmp_records = {q_name: q_seq}
            for lab in labels:
                tmp_records[lab] = ref_dict[lab]
            write_fasta(str(tmp_aln), tmp_records)

            # 4) build subtree refpkg with taxit create
            tmp_refpkg = tmpd_p / f"refpkg_{qi}"
            tmp_refpkg.mkdir(parents=True, exist_ok=True)

            # taxit create requires ref alignment and tree stats/info
            taxit_cmd = [
                "taxit",
                "create",
                "-P",
                str(tmp_refpkg),
                "-l",
                f"subtree_{qi}",
                "--aln-fasta",
                refpkg_files.ref_alignment,
                "--tree-file",
                str(tmp_tree),
                "--tree-stats",
                refpkg_files.tree_info,
            ]
            run_cmd(taxit_cmd)

            # 5) run pplacer on tmp alignment
            tmp_jplace = tmpd_p / f"place_{qi}.jplace"
            pplacer_cmd = [
                "pplacer",
                "-m",
                args.model,
                "-c",
                str(tmp_refpkg),
                "-o",
                str(tmp_jplace),
                "-j",
                str(max(1, int(args.threads))),
                str(tmp_aln),
            ]
            run_cmd(pplacer_cmd)

            # 6) remap placements from subtree edges to backbone edges
            with open(tmp_jplace, "r", encoding="utf-8") as fh:
                place_json = json.load(fh)

            added_tree, edge_dict = read_tree_newick_edge_tokens(place_json["tree"])

            # backbone_tree currently has edge numbers in labels; leaf_dict refers to the numbered backbone_tree?
            # leaf_dict was built from backbone_tree *before* edge numbers. Rebuild leaf_dict on numbered tree.
            # The numbering does not change leaf labels *before* '%%', but labels now have '%%'. We want mapping by base label.
            backbone_leaf_map = {}
            for n in backbone_tree.traverse_leaves():
                lab = n.get_label()
                if not lab:
                    continue
                base = lab.split("%%", 1)[0]
                backbone_leaf_map[base] = n

            for placement in place_json.get("placements", []):
                # placement: {"p": [[distal, edge_num, lwr, like, pendant], ...], "n": ["query_id", ...]}
                p_list = placement["p"]

                for p in p_list:
                    distal_len = float(p[0])
                    edge_num = str(p[1])

                    if edge_num not in edge_dict:
                        continue

                    right_n = edge_dict[edge_num]
                    left_n = right_n.get_parent()
                    if left_n is None:
                        continue

                    # find closest leaf on each side (subtree)
                    visited = {left_n}
                    right_leaf, _ = find_closest_leaf_or_path(right_n, visited=set(), y=None)
                    left_leaf, _ = find_closest_leaf_or_path(left_n, visited=visited, y=None)

                    # map those leaves to backbone nodes (by label)
                    rlab = right_leaf.get_label()
                    llab = left_leaf.get_label()
                    if not rlab or not llab:
                        continue

                    # strip any '%%' if present
                    rlab = rlab.split("%%", 1)[0]
                    llab = llab.split("%%", 1)[0]

                    br = backbone_leaf_map.get(rlab)
                    bl = backbone_leaf_map.get(llab)
                    if br is None or bl is None:
                        continue

                    # Find a path from bl to br in the backbone tree, then walk along it to locate target edge by length.
                    # We approximate the SCAMPP remapping logic: move from left leaf towards right leaf, subtracting edge lengths.
                    _, path = find_closest_leaf_or_path(bl, visited={bl}, y=br)

                    # path is a list of nodes encountered; we need edges along this path.
                    remaining = distal_len
                    target_edge_node = None

                    # path list contains nodes; interpret consecutive nodes as edges. Use child edge_length to parent.
                    # For safety, we traverse nodes and consume edge lengths where defined.
                    for node in path:
                        el = node.get_edge_length() or 0.0
                        if remaining <= el:
                            target_edge_node = node
                            break
                        remaining -= el

                    if target_edge_node is None:
                        # fall back: place at the last node on the path
                        target_edge_node = path[-1] if path else bl

                    # Determine backbone edge_num from numbered labels (child carries edge id)
                    tlabel = target_edge_node.get_label() or ""
                    if "%%" in tlabel:
                        _, back_edge = tlabel.split("%%", 1)
                    else:
                        # if no token, cannot remap reliably
                        continue

                    p[0] = float(distal_len)  # keep distal (approx)
                    p[1] = int(back_edge)  # remapped edge_num

                # Add (remapped) placement record to the final placements
                jplace["placements"].append(placement)

    # Write final jplace (GraftM-compatible output location/name)
    with open(final_jplace, "w", encoding="utf-8") as fh:
        json.dump(jplace, fh)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

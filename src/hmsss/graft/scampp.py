# pplacer_tax_scampp_backend.py
from __future__ import annotations

import re
import heapq
import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import warnings

warnings.filterwarnings("ignore")
import treeswift


# -----------------------------------------------------------------------------
# Minimal utils (needed subset) – adapted to be correct and self-contained
# -----------------------------------------------------------------------------

def read_fasta_to_dict(path: str) -> Dict[str, str]:
    """Read FASTA into dict[label]=sequence (keeps alignment chars)."""
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
    """Split into (ref, query) by whether label exists in backbone tree leaf labels."""
    ref: Dict[str, str] = {}
    query: Dict[str, str] = {}
    for k, v in aln_dict.items():
        if k in backbone_leaf_labels:
            ref[k] = v
        else:
            query[k] = v
    return ref, query


def read_and_split_fasta(
        fasta_path: str,
        backbone_leaf_labels: Set[str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Stream FASTA and split records into:
      - ref_dict: headers that exist in backbone_leaf_labels
      - q_dict  : all others

    This avoids materializing a full aln_dict in memory.
    """
    ref_dict: Dict[str, str] = {}
    q_dict: Dict[str, str] = {}

    label: Optional[str] = None
    seq_chunks: List[str] = []

    def flush_record():
        nonlocal label, seq_chunks
        if label is None:
            return
        seq = "".join(seq_chunks)
        if label in backbone_leaf_labels:
            ref_dict[label] = seq
        else:
            q_dict[label] = seq
        label = None
        seq_chunks = []

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):

                flush_record()
                label = line[1:].strip()
            else:
                seq_chunks.append(line)

    flush_record()
    return ref_dict, q_dict


def set_fragment_indices(x: str) -> Tuple[int, int]:
    """Trim leading/trailing '-' for fragment-aware distance."""
    si = 0
    ei = len(x)
    while si < ei and x[si] == "-":
        si += 1
    while ei > si and x[ei - 1] == "-":
        ei -= 1
    return si, ei


def hamming(a: str, b: str) -> int:
    return sum(1 for ca, cb in zip(a, b) if ca != cb)


def find_closest_hamming(x: str, ref: Dict[str, str], n: int, fragment_flag: bool) -> List[str]:
    """
    Return n closest ref labels by (fragment-aware) Hamming distance.
    Chunked evaluation to reduce peak compute.
    """
    if n <= 0 or not ref:
        return []
    si, ei = set_fragment_indices(x) if fragment_flag else (0, len(x))
    c = 200

    heap: List[Tuple[int, int, int, str]] = []
    counter = 0
    for name, seq in ref.items():
        first = hamming(seq[si:si + c], x[si:si + c])
        left = (ei - si) - c
        heapq.heappush(heap, (first, left, counter, name))
        counter += 1
    out: List[str] = []
    while heap:
        dist, left, _, name = heapq.heappop(heap)
        if left < 0:
            out.append(name)
            if len(out) >= n:
                return out
            continue
        ind = ei - left
        dist2 = hamming(ref[name][ind:ind + c], x[ind:ind + c])
        heapq.heappush(heap, (dist + dist2, left - c, counter, name))
        counter += 1

    return out


def subtree_nodes(tree: treeswift.Tree, leaf_y: treeswift.Node, n: int) -> List[str]:
    """Topological (unweighted) nearest-leaves expansion around seed leaf."""
    leaves = [leaf_y]
    visited = {leaf_y}
    heap: List[Tuple[int, int, treeswift.Node]] = []
    counter = 0

    p = leaf_y.get_parent()
    if p is not None:
        heapq.heappush(heap, (0, counter, p))
        counter += 1

    while heap and len(leaves) < n:
        _, _, node = heapq.heappop(heap)
        if node in visited:
            continue
        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        nbrs = list(node.child_nodes())
        if not node.is_root():
            nbrs.append(node.get_parent())

        for nb in nbrs:
            if nb is not None and nb not in visited:
                heapq.heappush(heap, (1, counter, nb))
                counter += 1

    return [x.get_label() for x in leaves if x.get_label() is not None]


def subtree_nodes_with_edge_length(tree: treeswift.Tree, leaf_y: treeswift.Node, n: int) -> List[str]:
    """Edge-length weighted nearest-leaves expansion around seed leaf."""
    leaves = [leaf_y]
    visited = {leaf_y}
    heap: List[Tuple[float, int, treeswift.Node]] = []
    counter = 0

    p = leaf_y.get_parent()
    if p is not None:
        heapq.heappush(heap, (leaf_y.get_edge_length() or 0.0, counter, p))
        counter += 1

    while heap and len(leaves) < n:
        dist, _, node = heapq.heappop(heap)
        if node in visited:
            continue
        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        # neighbors: parent + children
        if not node.is_root():
            par = node.get_parent()
            if par is not None and par not in visited:
                heapq.heappush(heap, (dist + (node.get_edge_length() or 0.0), counter, par))
                counter += 1

        for ch in node.child_nodes():
            if ch is not None and ch not in visited:
                heapq.heappush(heap, (dist + (ch.get_edge_length() or 0.0), counter, ch))
                counter += 1

    return [x.get_label() for x in leaves if x.get_label() is not None]


def add_edge_numbers(tree: treeswift.Tree) -> None:
    """Attach %%<id> to each node label so child-label encodes edge id."""
    counter = 0
    for node in tree.traverse_postorder():
        counter += 1
        lab = node.get_label()
        if lab is None:
            node.set_label(f"%%{counter}")
        else:
            node.set_label(f"{lab}%%{counter}")


def newick_with_edge_tokens(tree: treeswift.Tree) -> str:
    """
    Build jplace-style tree with edge tokens {edge_num} after each child edge length.
    Edge id is derived from child's label suffix '%%<id>'.
    """

    def child_edge_token(child: treeswift.Node) -> str:
        lab = child.get_label() or ""
        if "%%" in lab:
            _, tok = lab.split("%%", 1)
            return tok
        return "0"

    def rec(node: treeswift.Node) -> str:
        label = node.get_label() or ""
        base_label = label.split("%%", 1)[0] if "%%" in label else label

        if node.is_leaf():
            return base_label

        parts = []
        for ch in node.child_nodes():
            ch_str = rec(ch)
            el = ch.get_edge_length()
            if el is None:
                parts.append(ch_str)
            else:
                parts.append(f"{ch_str}:{el}{{{int(child_edge_token(ch))}}}")
        inside = ",".join(parts)
        return f"({inside}){base_label}" if base_label else f"({inside})"

    return rec(tree.root) + ";"


# Minimal parser for jplace tree strings (with {edge} tokens)
BRACKET = {"[": "]", "{": "}", "'": "'", '"': '"'}


def read_tree_newick_edge_tokens(newick_str: str) -> Tuple[treeswift.Tree, Dict[str, treeswift.Node]]:
    """
    Parse Newick containing edge tokens {...} and return (tree, edge_dict)
    where edge_dict[token] = child-node of that edge.
    """
    edge_dict: Dict[str, treeswift.Node] = {}
    s = newick_str.strip()
    if s.startswith("[&R]"):
        s = s.split("]", 1)[1].strip()

    t = treeswift.Tree()
    n = t.root
    i = 0
    while i < len(s):
        ch = s[i]
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
            num = ""
            while i < len(s) and s[i] not in ",);{":
                num += s[i]
                i += 1
            n.edge_length = float(num) if num else None
            i -= 1
        elif ch == "{":
            i += 1
            tok = ""
            while i < len(s) and s[i] != "}":
                tok += s[i]
                i += 1
            edge_dict[tok] = n
        else:
            # label
            label = ""
            bracket = None
            while i < len(s) and (
                    bracket is not None
                    or s[i] in BRACKET
                    or s[i] not in ":,;){"
            ):
                if s[i] in BRACKET and bracket is None:
                    bracket = s[i]
                elif bracket is not None and s[i] == BRACKET[bracket]:
                    bracket = None
                label += s[i]
                i += 1
            n.label = label if label else None
            i -= 1
        i += 1

    return t, edge_dict


def sanitize_newick_for_pplacer_inplace(path: str | Path) -> None:
    """
    In-place removal of Newick/NEXUS comment blocks (e.g. '[&R]')
    to make trees pplacer-compatible.
    """
    path = Path(path)

    text = path.read_text().strip()

    # Remove all [ ... ] comment blocks (non-greedy)
    text = re.sub(r"\[.*?\]", "", text)

    # Clean whitespace
    text = text.strip()

    # Ensure Newick terminator
    if not text.endswith(";"):
        text += ";"

    path.write_text(text + "\n")


def find_closest(x, visited, y=None):
    """ Returns leaf label for closest leaf to the node x through path not travelling through visited.
    If y is populated returns path from x to y not travelling through nodes in visited.

    Parameters
    ----------
    x : dendropy node object
    visited : list containing dendropy node objects
    y : dendropy node object

    Returns
    -------
    If y == None : dendropy node object of closest leaf y to the node x through path not travelling through nodes in visited,
                   list containing dendropy node objects on path to that leaf y from node x
    If y != None : dendropy node object y,
                   list containing dendropy node objects on path from node x to leaf y not travelling through nodes in visited

    """
    queue = []
    cnt = 1
    visited.add(x)

    if x.get_parent() and x.get_parent() not in visited:
        tmp = []
        tmp.append(x)
        heapq.heappush(queue, [x.get_edge_length(), cnt, tmp, x.get_parent()])
        cnt += 1

    for child in x.child_nodes():
        if child and child not in visited:
            tmp = []
            tmp.append(child)
            heapq.heappush(queue, [child.get_edge_length(), cnt, tmp, child])
            cnt += 1

    while len(queue) > 0:
        try:
            [length, _, path, node] = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf():
            if (not y) or node.get_label() == y.get_label():
                return node, path
            else:
                continue

        if node.get_parent() and node.get_parent() not in visited:
            tmp = path.copy()
            tmp.append(node)
            heapq.heappush(queue, [length + node.get_edge_length(), cnt, tmp, node.get_parent()])
            cnt += 1

        for child in node.child_nodes():
            if child and child not in visited:
                tmp = path.copy()
                tmp.append(child)
                heapq.heappush(queue, [length + child.get_edge_length(), cnt, tmp, child])
                cnt += 1

    return x, [x]


def run_cmd(cmd: List[str], *, cwd: Optional[str] = None) -> None:
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  cwd: {cwd}\n"
            f"  exit: {p.returncode}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )


def ensure_on_path(exe: str) -> str:
    path = shutil.which(exe)
    if not path:
        raise RuntimeError(f"Required executable not found on PATH: {exe}. Activate the conda env that provides it.")
    return path


def write_fasta(path: str, records: Dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for k, v in records.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 80):
                fh.write(v[i:i + 80] + "\n")


@dataclass(frozen=True)
class RefPkgFiles:
    tree_file: str
    ref_alignment: str
    tree_stats: str  # --tree-stats input for taxit create


def resolve_refpkg_files(refpkg_dir: str) -> RefPkgFiles:
    """
    Resolve required components from a GraftM refpkg directory, given your file layout:
      - *.tree
      - deduplicated_aligned.fasta
      - *.json / *.tree_log (tree stats)
    """
    rp = Path(refpkg_dir)
    if not rp.exists() or not rp.is_dir():
        raise FileNotFoundError(f"refpkg directory not found: {refpkg_dir}")

    tree_candidates = sorted(rp.glob("*.tree"))
    if not tree_candidates:
        raise RuntimeError(f"No '*.tree' found in refpkg: {refpkg_dir}")
    tree_file = str(tree_candidates[0])

    ref_aln = sorted(rp.glob("*_aligned.fasta"))[0]
    if not ref_aln.exists():
        raise RuntimeError(f"Expected 'deduplicated_aligned.fasta' in refpkg: {refpkg_dir}")

    # Prefer *.tree_log if present; else first *.json; else error.
    # (You can tighten this to your exact tree-stats choice once confirmed.)
    tree_log = sorted(rp.glob("*.tree_log"))
    jsons = sorted(rp.glob("*.json"))

    if tree_log:
        tree_stats = str(tree_log[0])
    elif jsons:
        tree_stats = str(jsons[0])
    else:
        raise RuntimeError(f"No '*.tree_log' or '*.json' found in refpkg: {refpkg_dir}")

    return RefPkgFiles(tree_file=tree_file, ref_alignment=str(ref_aln), tree_stats=tree_stats)


# -----------------------------------------------------------------------------
# Public routine: same *call shape* as GraftM pplacer() + returns full jplace dict
# -----------------------------------------------------------------------------

def pplacer_tax_scampp_like_graftm(
        *,
        output_file: str,
        output_path: str,
        input_path: str,
        threads: int,
        refpkg: str,
        # SCAMPP knobs (defaults match typical tax-SCAMPP usage)
        model: str = "GTR",
        subtreesize: int = 2000,
        subtreetype: str = "d",  # "d" edge-length weighted, "n" topological, "h" take top-n by Hamming directly
        fragmentflag: bool = True,
        tmpfilenbr: int = 0,
) -> Path:
    """
    Run pplacer-tax-SCAMPP style placement but with GraftM-like inputs.

    Inputs match the GraftM pplacer() call:
      - output_file: basename for <output_file>.jplace
      - output_path: directory to write the final jplace
      - input_path : combined alignment (ref + queries)
      - threads    : pplacer -j threads
    Additional:
      - refpkg     : directory containing *.tree, deduplicated_aligned.fasta, etc.

    Returns:
      - The fully assembled final jplace JSON as a Python dict (also written to disk).
    """
    # ensure_on_path("pplacer")
    # ensure_on_path("taxit")
    CANON_FIELDS = ["distal_length", "edge_num", "like_weight_ratio", "likelihood", "pendant_length"]
    outdir = Path(output_path)
    outdir.mkdir(parents=True, exist_ok=True)
    final_jplace_path = outdir / f"{output_file}.jplace"

    ref = resolve_refpkg_files(refpkg)

    # Load backbone tree
    backbone_tree = treeswift.read_tree_newick(ref.tree_file)
    backbone_leaf_labels = {n.get_label() for n in backbone_tree.traverse_leaves() if n.get_label() is not None}

    # Read combined alignment and split into ref vs query by tree labels
    # aln_dict = read_fasta_to_dict(input_path)
    # ref_dict, q_dict = read_and_split_fasta(input_path, backbone_leaf_labels)
    q_dict = read_fasta_to_dict(input_path)  # Queries only
    ref_dict = read_fasta_to_dict(ref.ref_alignment)  # e.g. deduplicated_aligned.fasta

    print("tree leaves (sample):", list(sorted(backbone_leaf_labels))[:5])
    print("ref_dict size:", len(ref_dict), "q_dict size:", len(q_dict))
    print("intersection size:", len(set(ref_dict.keys()) & backbone_leaf_labels))
    # Prepare numbered backbone tree + jplace scaffold
    add_edge_numbers(backbone_tree)
    jplace = {
        # "tree": newick_with_edge_tokens(backbone_tree),
        "placements": [],
        "metadata": {
            "invocation": "pplacer_tax_scampp_like_graftm",
            "refpkg": str(Path(refpkg).resolve()),
            "model": model,
            "subtreesize": subtreesize,
            "subtreetype": subtreetype,
            "fragmentflag": fragmentflag,
            "tmpfilenbr": tmpfilenbr,
        },
        "version": 3,
    }
    jplace["fields"] = CANON_FIELDS

    # Map base leaf labels -> numbered backbone leaf nodes (for remap)
    numbered_leaf_map: Dict[str, treeswift.Node] = {}
    for n in backbone_tree.traverse_leaves():
        lab = n.get_label()
        if not lab:
            continue
        base = lab.split("%%", 1)[0]
        numbered_leaf_map[base] = n
    with tempfile.TemporaryDirectory(prefix=f"tax_scampp_{tmpfilenbr}_") as tmpd:
        tmpd_p = Path(tmpd)

        # Iterate queries: SCAMPP per query
        for qi, (q_name, q_seq) in enumerate(q_dict.items(), start=1):
            # For subtree extraction, use an unmodified copy (no edge-number labels)
            backbone_plain = treeswift.read_tree_newick(ref.tree_file)
            plain_leaf_index = backbone_plain.label_to_node(selection="leaves")

            # 1) choose subtree leaf labels
            if subtreetype == "h":
                labels = find_closest_hamming(q_seq, ref_dict, subtreesize, fragmentflag)
            else:
                nearest = find_closest_hamming(q_seq, ref_dict, 1, fragmentflag)
                if not nearest:
                    continue
                seed_label = nearest[0]
                seed_node = plain_leaf_index.get(seed_label)
                if seed_node is None:
                    continue
                if subtreetype == "n":
                    labels = subtree_nodes(backbone_plain, seed_node, subtreesize)
                else:
                    labels = subtree_nodes_with_edge_length(backbone_plain, seed_node, subtreesize)
            labels = [lab for lab in labels if lab in ref_dict]
            if not labels:
                continue
            # 2) subtree tree
            subtree = backbone_plain.extract_tree_with(labels)
            subtree.resolve_polytomies()
            tmp_tree = tmpd_p / f"subtree_{qi}.nwk"
            subtree.write_tree_newick(str(tmp_tree))
            # 3) tmp alignment: query + subtree refs (already aligned in combined alignment)
            tmp_aln = tmpd_p / f"aln_{qi}.fasta"
            records = {q_name: q_seq}
            for lab in labels:
                records[lab] = ref_dict[lab]
            write_fasta(str(tmp_aln), records)
            # 4) build subtree refpkg with taxit
            tmp_refpkg = tmpd_p / f"refpkg_{qi}"
            # tmp_refpkg.mkdir(parents=True, exist_ok=True)
            sanitize_newick_for_pplacer_inplace(tmp_tree)
            taxit_cmd = [
                "taxit", "create",
                "-P", str(tmp_refpkg),
                "-l", f"subtree_{qi}",
                "--aln-fasta", ref.ref_alignment,
                "--tree-file", str(tmp_tree),
                "--tree-stats", ref.tree_stats,
            ]
            run_cmd(taxit_cmd)

            # 5) pplacer on subtree
            tmp_jplace = tmpd_p / f"place_{qi}.jplace"
            pplacer_cmd = [
                "pplacer",
                "-m", model,
                "-c", str(tmp_refpkg),
                "-o", str(tmp_jplace),
                "-j", str(max(1, int(threads))),
                str(tmp_aln),
            ]
            print("5 Running pplacer")
            run_cmd(pplacer_cmd)

            placements = []
            tmp_output = tmp_jplace
            print("6 load subtree")
            # 6) load subtree placements and remap edge_num onto backbone

            # load the jplace file and find placements in the original backbone tree
            place_file = open(tmp_output, 'r')
            place_json = json.load(place_file)

            if len(place_json["placements"]) > 0:

                added_tree, edge_dict = read_tree_newick_edge_tokens(place_json["tree"])

                tmp_place = place_json["placements"][0]
                for i in range(len(tmp_place["p"])):
                    edge_num = tmp_place["p"][i][1]  # edge number in subtree
                    edge_distal = tmp_place["p"][i][0]  # distal length from parent node

                    # find placement edge according to edge number
                    right_n = edge_dict[str(edge_num)]
                    left_n = right_n.get_parent()

                    # obtain a path from leaf left to leaf right containing placement edge through the subtree
                    left, path_l = find_closest(left_n, {left_n, right_n})
                    right, path_r = find_closest(right_n, {left_n, right_n})

                    # obtain the corresponding path in backbone tree
                    left = plain_leaf_index[left.get_label()]
                    right = plain_leaf_index[right.get_label()]
                    _, path = find_closest(left, {left}, y=right)

                    # find the length of placement along the path from leaf left to leaf right in subtree
                    length = sum([x.get_edge_length() for x in path_l]) + edge_distal

                    # find the target placement edge in backbone tree
                    target_edge = path[-1]
                    for j in range(len(path)):
                        length -= path[j].get_edge_length()
                        if length < 0:
                            target_edge = path[j]
                            break

                    tmp_place["p"][i][0] = 0

                    label = target_edge.get_label()
                    [taxon, target_edge_nbr] = label.split('%%', 1)
                    tmp_place["p"][i][0] = target_edge.get_edge_length() + length
                    tmp_place["p"][i][1] = int(target_edge_nbr)

                # append the placement to the output jplace
                placements.append(tmp_place.copy())

            place_file.close()
        print("14")
        # build jplace file
        jplace["placements"] = placements
        jplace["metadata"]["invocation"] = {
            "routine": "pplacer_tax_scampp_like_graftm",
            "output_file": output_file,
            "output_path": str(Path(output_path).resolve()),
            "input_path": str(Path(input_path).resolve()),
            "threads": int(threads),
            "refpkg": str(Path(refpkg).resolve()),
            "model": model,
            "subtreesize": subtreesize,
            "subtreetype": subtreetype,
            "fragmentflag": fragmentflag,
            "tmpfilenbr": tmpfilenbr,
        }
        jplace["version"] = 3
        jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
                            "likelihood", "pendant_length"]

        # output = open('{}/{}.jplace'.format(output, outFile), 'w')
        final_jplace_path = (Path(output_path) / f"{output_file}.jplace").resolve()
        print("15")
        with open(final_jplace_path, "w", encoding="utf-8") as fh:
            json.dump(jplace, fh, sort_keys=True, indent=2)
            fh.write("\n")
        return final_jplace_path

# tests/test_settrie_completeness.py
from __future__ import annotations
from dataclasses import dataclass
import math

# -----------------------------
# Minimal-Implementierung (inline)
# -----------------------------


class Node:
    __slots__ = ("children", "terminals", "sub_lo", "sub_hi", "max_terminal_depth")

    def __init__(self):
        self.children = {}
        self.terminals = []
        self.sub_lo = -1
        self.sub_hi = -1
        self.max_terminal_depth = 0


@dataclass(frozen=True)
class PatternMeta:
    name: str
    length: int
    mask: int
    ids_sorted: tuple[int, ...]


@dataclass
class TrieIndex:
    root: Node
    TERMS: list[int]
    pattern_meta: list[PatternMeta]
    domain_to_id: dict[str, int]
    id_to_domain: list[str]
    rank: list[int]


def _build_domain_vocab(pattern_dict: dict[str, set[str]]):
    all_domains = sorted({d for pat in pattern_dict.values() for d in pat})
    domain_to_id = {d: i for i, d in enumerate(all_domains)}
    id_to_domain = all_domains
    return domain_to_id, id_to_domain


def _domain_frequency_in_patterns(pattern_dict: dict[str, set[str]]):
    from collections import Counter

    c = Counter()
    for pat in pattern_dict.values():
        c.update(set(pat))
    return c


def _build_rank(domain_to_id: dict[str, int], pattern_dict: dict[str, set[str]]):
    freq = _domain_frequency_in_patterns(pattern_dict)
    # seltene zuerst, bei Gleichstand alphabetisch
    ordered_names = sorted(domain_to_id.keys(), key=lambda d: (freq[d], d))
    name_to_rank = {name: r for r, name in enumerate(ordered_names)}
    rank = [0] * len(domain_to_id)
    for name, did in domain_to_id.items():
        rank[did] = name_to_rank[name]
    return rank


def _encode_mask(ids):
    m = 0
    for i in ids:
        m |= 1 << i
    return m


def _build_pattern_meta(pattern_dict, domain_to_id, rank):
    pattern_meta = []
    for pid, (name, doms) in enumerate(pattern_dict.items()):
        ids = [domain_to_id[d] for d in doms if d in domain_to_id]
        ids_sorted = tuple(sorted(ids, key=lambda x: rank[x]))
        meta = PatternMeta(
            name=name,
            length=len(ids_sorted),
            mask=_encode_mask(ids_sorted),
            ids_sorted=ids_sorted,
        )
        pattern_meta.append(meta)
    return pattern_meta


def _insert_pattern(root: Node, ids_sorted: tuple[int, ...], pid: int):
    n = root
    depth = 0
    for d in ids_sorted:
        n = n.children.setdefault(d, Node())
        depth += 1
    n.terminals.append(pid)
    if depth > n.max_terminal_depth:
        n.max_terminal_depth = depth


def _propagate_max_terminal_depth(root: Node):
    mtd = root.max_terminal_depth
    for child in root.children.values():
        cmtd = _propagate_max_terminal_depth(child)
        if cmtd > mtd:
            mtd = cmtd
    root.max_terminal_depth = mtd
    return mtd


def _dfs_flatten_terminals(root: Node):
    TERMS: list[int] = []

    def dfs(n: Node):
        n.sub_lo = len(TERMS)
        if n.terminals:
            TERMS.extend(n.terminals)
        for d in sorted(n.children.keys()):
            dfs(n.children[d])
        n.sub_hi = len(TERMS)

    dfs(root)
    return TERMS


def build_trie_index(pattern_dict: dict[str, set[str]]) -> TrieIndex:
    if not pattern_dict:
        return TrieIndex(
            root=Node(),
            TERMS=[],
            pattern_meta=[],
            domain_to_id={},
            id_to_domain=[],
            rank=[],
        )
    domain_to_id, id_to_domain = _build_domain_vocab(pattern_dict)
    rank = _build_rank(domain_to_id, pattern_dict)
    pattern_meta = _build_pattern_meta(pattern_dict, domain_to_id, rank)
    root = Node()
    for pid, meta in enumerate(pattern_meta):
        if meta.length == 0:
            continue
        _insert_pattern(root, meta.ids_sorted, pid)
    _propagate_max_terminal_depth(root)
    TERMS = _dfs_flatten_terminals(root)
    return TrieIndex(
        root=root,
        TERMS=TERMS,
        pattern_meta=pattern_meta,
        domain_to_id=domain_to_id,
        id_to_domain=id_to_domain,
        rank=rank,
    )


def map_present_to_sorted_ids(present_domains: set[str], index: TrieIndex):
    """Gibt (ids_sorted, present_mask, present_names) zurück.
    present_names ist die **Original**-Menge (inkl. unbekannter Domänen)."""
    did = index.domain_to_id
    rk = index.rank
    ids = [did[d] for d in present_domains if d in did]
    ids_sorted = sorted(ids, key=lambda x: rk[x])
    present_mask = 0
    for i in ids_sorted:
        present_mask |= 1 << i
    # WICHTIG: Originale Namen-Menge behalten (inkl. unbekannter)
    present_names = set(present_domains)
    return ids_sorted, present_mask, present_names


def descend_last_reachable(root: Node, present_ids_sorted: list[int]):
    node = root
    k = 0
    for d in present_ids_sorted:
        child = node.children.get(d)
        if child is None:
            continue
        node = child
        k += 1
    return node, k


def iter_subtree_pattern_ids(node: Node, TERMS: list[int]):
    if node.sub_lo < 0 or node.sub_hi < 0:
        return ()
    return TERMS[node.sub_lo : node.sub_hi]


# -----------------------------
# Der eigentliche Test
# -----------------------------


def test_trie_completeness_and_missing_additional():
    # 1) Pattern-Korpus
    pattern_dict = {
        "KW_A": {"A", "B", "E"},
        "KW_B": {"A", "D"},
        "KW_C": {"C", "E"},
    }
    index = build_trie_index(pattern_dict)

    # 2) Cluster-Fälle (inkl. unbekannter Domänen X,Y)
    clusters = [
        ("{A,B}", {"A", "B"}),
        ("{A,D}", {"A", "D"}),
        ("{A,B,E}", {"A", "B", "E"}),
        ("{C,E}", {"C", "E"}),
        ("{X,Y}", {"X", "Y"}),  # komplett unbekannt
        ("{A,B,X}", {"A", "B", "X"}),  # teilweise unbekannt
        ("{B,E}", {"B", "E"}),  # Teilmenge von KW_A
        ("{C}", {"C"}),  # Teilmenge von KW_C
    ]

    print("\n=== Ergebnisse ===")
    for label, present in clusters:
        present_ids_sorted, present_mask, present_names = map_present_to_sorted_ids(
            present, index
        )
        node, k = descend_last_reachable(index.root, present_ids_sorted)
        pids = list(iter_subtree_pattern_ids(node, index.TERMS))

        print(f"\nCluster {label} -> Tiefe k={k}")
        if not pids:
            print("  (Keine Kandidaten im Trie-Subtree)")
            continue

        # Hilfsmenge: Pattern-Domänen als Namen
        def pattern_domains_as_names(meta: PatternMeta) -> set[str]:
            return {index.id_to_domain[i] for i in meta.ids_sorted}

        # Ausgabe sortiert nach Pattern-Namen
        for pid in sorted(pids, key=lambda p: index.pattern_meta[p].name):
            meta = index.pattern_meta[pid]
            completeness = k / meta.length if meta.length else 0.0

            pat_names = pattern_domains_as_names(meta)

            # NEU: missing/additional in **Namens-Logik** (inkl. unbekannter present-Namen)
            missing_names = sorted(pat_names - present_names)
            additional_names = sorted(present_names - pat_names)

            print(
                f"  Pattern {meta.name}: "
                f"completeness={completeness:.3f}, "
                f"missing={missing_names}, additional={additional_names}"
            )

    # -----------------
    # Assertions (Kernaussagen)
    # -----------------

    # {A,B} -> KW_A: completeness 2/3, missing ['E'], additional []
    ids_sorted, pmask, pnames = map_present_to_sorted_ids({"A", "B"}, index)
    node, k = descend_last_reachable(index.root, ids_sorted)
    pids = list(iter_subtree_pattern_ids(node, index.TERMS))
    meta = index.pattern_meta[pids[0]]
    assert meta.name == "KW_A"
    assert math.isclose(k / meta.length, 2 / 3, rel_tol=1e-9)
    pat_names = {index.id_to_domain[i] for i in meta.ids_sorted}
    assert sorted(pat_names - pnames) == ["E"]
    assert sorted(pnames - pat_names) == []  # additional leer

    # {A,B,X} -> KW_A: additional enthält 'X'
    ids_sorted, pmask, pnames = map_present_to_sorted_ids({"A", "B", "X"}, index)
    node, k = descend_last_reachable(index.root, ids_sorted)
    pids = list(iter_subtree_pattern_ids(node, index.TERMS))
    meta = index.pattern_meta[pids[0]]
    pat_names = {index.id_to_domain[i] for i in meta.ids_sorted}
    assert "X" in (set(pnames) - pat_names)  # unknown zählt als additional

    # {X,Y} -> alle Patterns, completeness 0.0; additional = ['X','Y'] bei jedem Pattern
    ids_sorted, pmask, pnames = map_present_to_sorted_ids({"X", "Y"}, index)
    node, k = descend_last_reachable(index.root, ids_sorted)
    pids = list(iter_subtree_pattern_ids(node, index.TERMS))
    for pid in pids:
        meta = index.pattern_meta[pid]
        assert math.isclose(k / meta.length, 0.0, rel_tol=1e-9)
        pat_names = {index.id_to_domain[i] for i in meta.ids_sorted}
        assert sorted(pnames - pat_names) == ["X", "Y"]

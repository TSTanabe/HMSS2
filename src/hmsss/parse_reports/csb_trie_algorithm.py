from __future__ import annotations
from dataclasses import dataclass
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional, DefaultDict


# =========================
# Eingabe-Interface
# =========================
# Erwartet wird ein pattern_dict im Format:
#   pattern_dict[name] = list[str] oder set[str]  (Domain-Bezeichner)
#
# Beispiel:
# pattern_dict = {
#   "KW_A": {"A","B","E"},
#   "KW_B": {"A","D"},
#   "KW_C": {"C","E"},
# }

# =========================
# Kern-Datenstrukturen
# =========================


class Node:
    """Leichter Trie-Knoten. Kinder über domain_id (int) -> Node."""

    __slots__ = ("children", "terminals", "sub_lo", "sub_hi", "max_terminal_depth")

    def __init__(self):
        self.children: Dict[
            int, Node
        ] = {}  # Kanten zu Kindknoten, beschriftet mit domain ids
        self.terminals: List[int] = []  # pattern_ids, die in diesem Konten enden
        self.sub_lo: int = -1  # inklusiver Start-Index in TERMS
        self.sub_hi: int = -1  # exklusiver End-Index in TERMS
        self.max_terminal_depth: int = 0  # maximale Pattern-Länge im Teilbaum


@dataclass(frozen=True)
class PatternMeta:
    name: str  # keyword
    length: int  # Länge des pattern
    mask: int  # bitset der domain
    ids_sorted: Tuple[int, ...]  # Domain-IDs in globaler Ordnung


@dataclass
class TrieIndex:
    # Persistierbare, read-only Struktur für die Laufzeit-Abfragen
    root: Node
    TERMS: List[int]  # globale Terminal-Liste (pattern_ids), DFS-Flattened
    pattern_meta: List[PatternMeta]  # Metadaten je pattern_id
    domain_to_id: Dict[str, int]
    id_to_domain: List[str]
    rank: List[int]  # rank[domain_id] -> Ordnungsrang (kleiner = früher)
    nodes_by_domain: Dict[int, List[Node]]  # Reverse-Index für Anywhere-Start


# =========================
# Hilfsfunktionen: Vokabular & Ordnung
# =========================


def _build_domain_vocab(
    pattern_dict: Dict[str, Iterable[str]],
) -> Tuple[Dict[str, int], List[str]]:
    """Vergibt dichte Integer-IDs für Domains in stabiler Reihenfolge (alphabetisch)."""
    all_domains = sorted({d for pat in pattern_dict.values() for d in pat})
    domain_to_id = {d: i for i, d in enumerate(all_domains)}
    id_to_domain = all_domains
    return domain_to_id, id_to_domain


def _domain_frequency_in_patterns(pattern_dict: Dict[str, Iterable[str]]) -> Counter:
    """Häufigkeit je Domain über Patterns (Anzahl Patterns, die die Domain enthalten)."""
    c = Counter()
    for pat in pattern_dict.values():
        seen = set(pat)
        c.update(seen)
    return c


def _build_rank(
    domain_to_id: Dict[str, int], pattern_dict: Dict[str, Iterable[str]]
) -> List[int]:
    """Bestimmt globale Ordnung: seltener -> früher, bei Gleichstand alphabetisch stabil."""
    freq = _domain_frequency_in_patterns(pattern_dict)
    # (häufigkeit, domain_name) -> sortieren, dann auf IDs mappen
    ordered_names = sorted(domain_to_id.keys(), key=lambda d: (freq[d], d))
    name_to_rank = {name: r for r, name in enumerate(ordered_names)}
    # rank als Liste über domain_id
    rank = [0] * len(domain_to_id)
    for name, did in domain_to_id.items():
        rank[did] = name_to_rank[name]
    return rank


# =========================
# Pattern-Metadaten
# =========================


def _encode_mask(ids: Iterable[int]) -> int:
    m = 0
    for i in ids:
        m |= 1 << i
    return m


def _build_pattern_meta(
    pattern_dict: Dict[str, Iterable[str]],
    domain_to_id: Dict[str, int],
    rank: List[int],
) -> Tuple[List[PatternMeta], Dict[str, int]]:
    """Erzeugt PatternMeta pro Pattern und ein name->pid Mapping."""
    pattern_meta: List[PatternMeta] = []
    name_to_pid: Dict[str, int] = {}
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
        name_to_pid[name] = pid
    return pattern_meta, name_to_pid


# =========================
# Trie-Aufbau & Flattening
# =========================


def _insert_pattern(root: Node, ids_sorted: Tuple[int, ...], pid: int) -> None:
    n = root
    depth = 0
    for d in ids_sorted:
        n = n.children.setdefault(d, Node())
        depth += 1
    n.terminals.append(pid)
    if depth > n.max_terminal_depth:
        n.max_terminal_depth = depth


def _propagate_max_terminal_depth(root: Node) -> int:
    """Falls Terminals nur im Blatt erfasst wurden: aggregiere max_terminal_depth nach oben."""
    mtd = root.max_terminal_depth
    for child in root.children.values():
        cmtd = _propagate_max_terminal_depth(child)
        if cmtd > mtd:
            mtd = cmtd
    root.max_terminal_depth = mtd
    return mtd


def _dfs_flatten_terminals(root: Node) -> List[int]:
    """Erstellt die globale TERMS-Liste und setzt an jedem Knoten (sub_lo, sub_hi)."""
    TERMS: List[int] = []

    def dfs(n: Node):
        n.sub_lo = len(TERMS)
        if n.terminals:
            TERMS.extend(n.terminals)
        # Kinder in aufsteigender Domain-ID (oder besser: nach Rang, wenn gewünscht)
        for d in sorted(n.children.keys()):
            dfs(n.children[d])
        n.sub_hi = len(TERMS)

    dfs(root)
    return TERMS


def _collect_nodes_by_domain(root: Node) -> Dict[int, List[Node]]:
    bucket: DefaultDict[int, List[Node]] = defaultdict(list)
    stack = [root]
    while stack:
        n = stack.pop()
        for did, child in n.children.items():
            bucket[did].append(child)
            stack.append(child)
    return dict(bucket)


# =========================
# Öffentliche API
# =========================


def build_trie_index(pattern_dict: Dict[str, Iterable[str]]) -> TrieIndex:
    """
    Haupt-Einstieg: Baut die komplette Index-Struktur (Trie + Metas) aus dem pattern_dict.
    """
    if not pattern_dict:
        # Leere Struktur
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
    pattern_meta, name_to_pid = _build_pattern_meta(pattern_dict, domain_to_id, rank)

    root = Node()
    # Einfügen aller Patterns in den Trie
    for pid, meta in enumerate(pattern_meta):
        if meta.length == 0:
            # leeres Pattern ignorieren
            continue
        _insert_pattern(root, meta.ids_sorted, pid)

    # max_terminal_depth bottom-up aggregieren
    _propagate_max_terminal_depth(root)
    # TERMS erzeugen und sub_lo/sub_hi setzen
    TERMS = _dfs_flatten_terminals(root)
    nodes_by_domain = _collect_nodes_by_domain(root)

    return TrieIndex(
        root=root,
        TERMS=TERMS,
        pattern_meta=pattern_meta,
        domain_to_id=domain_to_id,
        id_to_domain=id_to_domain,
        rank=rank,
        nodes_by_domain=nodes_by_domain,
    )


# =========================
# Laufzeit-Helfer (optional)
# =========================


def map_present_to_sorted_ids(present_domains: set[str], index: TrieIndex):
    """Mappt Cluster-Domainnamen -> sortierte IDs nach globaler Ordnung, gibt auch die Bitmaske zurück."""
    did = index.domain_to_id
    rk = index.rank
    ids = [did[d] for d in present_domains if d in did]
    ids_sorted = sorted(ids, key=lambda x: rk[x])
    present_mask = 0
    for i in ids_sorted:
        present_mask |= 1 << i
    # NEU: Original-Namen (inkl. unbekannter Domains) mit zurückgeben
    present_names = set(present_domains)
    # Optional: unknowns separat, falls du sie explizit loggen willst
    unknown_names = {d for d in present_domains if d not in did}
    return ids_sorted, present_mask, present_names, unknown_names


def descend_last_reachable(
    root: Node, present_ids_sorted: List[int]
) -> Tuple[Node, int]:
    """
    Geht so tief wie möglich entlang der vorhandenen Domains.
    Liefert den letzten erreichbaren Knoten und dessen Tiefe (k = Trefferzahl).
    """
    node = root
    k = 0
    # Wir „springen“ über nicht existierende Kinder einfach hinweg (sie sind im Trie nicht vorhanden).
    for d in present_ids_sorted:
        child = node.children.get(d)
        if child is None:
            continue
        node = child
        k += 1
    return node, k


def iter_subtree_pattern_ids(node: Node, TERMS: List[int]) -> Iterable[int]:
    """Liefert alle pattern_ids im Teilbaum des Knotens (TERMS-Slice)."""
    if node.sub_lo < 0 or node.sub_hi < 0:
        return ()
    return TERMS[node.sub_lo : node.sub_hi]


# + in csb_trie_algorithm.py
def iter_candidates_anywhere_start(
    index: TrieIndex,
    present_ids_sorted: List[int],
    present_mask: int,
) -> Dict[int, int]:
    """
    Anywhere-Start: Für jede Cluster-Domäne an allen passenden Knoten starten.
    Liefert dict: pattern_id -> bestes k (maximale erreichte Tiefe ab Startknoten).
    """
    best_k_by_pid: Dict[int, int] = {}
    ids = present_ids_sorted
    nodes_by_domain = index.nodes_by_domain
    TERMS = index.TERMS
    meta = index.pattern_meta

    for start_id in ids:
        for start_node in nodes_by_domain.get(start_id, ()):
            node = start_node
            k = 1  # Einstieg über start_id zählt als 1
            # restliche Cluster-Domänen in der gleichen globalen Reihenfolge prüfen
            for d in ids:
                if d == start_id:
                    continue
                nxt = node.children.get(d)
                if nxt is not None:
                    node = nxt
                    k += 1

            # PIDs im Subtree + schneller Overlap-Filter
            slice_iter = iter_subtree_pattern_ids(node, TERMS)
            for pid in slice_iter:
                if meta[pid].mask & present_mask:
                    prev = best_k_by_pid.get(pid)
                    if (prev is None) or (k > prev):
                        best_k_by_pid[pid] = k

    return best_k_by_pid


# =========================
# (Optionale) Persistenz
# =========================


def save_trie_index(index: TrieIndex, outdir: Path) -> None:
    """Speichert Index in einfachem, robustem Format (Pickle + NumPy)."""
    import pickle, json, numpy as np

    outdir.mkdir(parents=True, exist_ok=True)

    # Trie als Pickle
    with open(outdir / "trie.pkl", "wb") as f:
        pickle.dump(index.root, f, protocol=pickle.HIGHEST_PROTOCOL)

    # TERMS
    np.save(outdir / "TERMS.npy", np.asarray(index.TERMS, dtype=np.uint32))

    # Vokabular & Ordnung
    with open(outdir / "domain_to_id.json", "w") as f:
        json.dump(index.domain_to_id, f)
    with open(outdir / "id_to_domain.json", "w") as f:
        json.dump(index.id_to_domain, f)
    with open(outdir / "rank.npy", "wb") as f:
        np.save(f, np.asarray(index.rank, dtype=np.uint32))

    # Pattern-Metas (kompakt zerlegt)
    names = [m.name for m in index.pattern_meta]
    lengths = [m.length for m in index.pattern_meta]
    masks = [m.mask for m in index.pattern_meta]
    ids_sorted = [list(m.ids_sorted) for m in index.pattern_meta]

    with open(outdir / "pattern_names.pkl", "wb") as f:
        pickle.dump(names, f, protocol=pickle.HIGHEST_PROTOCOL)
    np.save(outdir / "pattern_lengths.npy", lengths)
    with open(outdir / "pattern_masks.pkl", "wb") as f:
        pickle.dump(masks, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(outdir / "pattern_ids_sorted.pkl", "wb") as f:
        pickle.dump(ids_sorted, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_trie_index(indir: Path) -> TrieIndex:
    """Lädt den gespeicherten Index."""
    import pickle, json, numpy as np

    with open(indir / "trie.pkl", "rb") as f:
        root = pickle.load(f)
    TERMS = np.load(indir / "TERMS.npy", mmap_mode=None).tolist()
    with open(indir / "domain_to_id.json") as f:
        domain_to_id = json.load(f)
    with open(indir / "id_to_domain.json") as f:
        id_to_domain = json.load(f)
    rank = np.load(indir / "rank.npy").tolist()

    with open(indir / "pattern_names.pkl", "rb") as f:
        names = pickle.load(f)
    lengths = np.load(indir / "pattern_lengths.npy").tolist()
    with open(indir / "pattern_masks.pkl", "rb") as f:
        masks = pickle.load(f)
    with open(indir / "pattern_ids_sorted.pkl", "rb") as f:
        ids_sorted_lists = pickle.load(f)

    pattern_meta = [
        PatternMeta(
            name=names[i],
            length=int(lengths[i]),
            mask=int(masks[i]),
            ids_sorted=tuple(ids_sorted_lists[i]),
        )
        for i in range(len(names))
    ]
    return TrieIndex(
        root=root,
        TERMS=TERMS,
        pattern_meta=pattern_meta,
        domain_to_id=domain_to_id,
        id_to_domain=id_to_domain,
        rank=rank,
    )

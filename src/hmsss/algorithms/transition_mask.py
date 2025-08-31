from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, Iterable, Tuple


@dataclass(frozen=True)
class TransitionMasks:
    """
    Datenstruktur für ultraschnellen Coverage-Check per Bitmasken.

    missing_id:  Missing-Typ -> Bit-Index (0..N-1)
    allow_mask:  Additional-Typ -> Bitmaske erlaubter Missing-Typen
    """

    missing_id: Dict[str, int]
    allow_mask: Dict[str, int]


def load_transition_masks_from_tsv(
    path: str,
    *,
    delimiter: str = "\t",
) -> TransitionMasks:
    """
    Lese eine TSV-Datei mit zwei Spalten (additional<TAB>missing) ein
    und baue TransitionMasks für schnelle Abdeckbarkeitsprüfungen.

    Parameter
    ---------
    path : str
        Pfad zur TSV-Datei.
    delimiter : str
        Spaltentrenner (default: Tab).
    skip_header : bool | int
        Wenn True, wird die erste Zeile übersprungen.
        Wenn int>0, werden so viele Zeilen am Anfang übersprungen.
    strip_spaces : bool
        Führende/nachgestellte Leerzeichen pro Feld entfernen.
    ignore_empty : bool
        Leere Felder/Zeilen ignorieren.

    Returns
    -------
    TransitionMasks
        missing_id-Map und allow_mask-Map (Bitmasken).
    """
    pairs: list[Tuple[str, str]] = []

    with open(path, "r", newline="") as f:
        rdr = csv.reader(f, delimiter=delimiter)

        for row in rdr:
            if not row:
                continue

            if len(row) < 2:
                continue

            add_raw, miss_raw = row[0], row[1]
            add_raw = add_raw.strip()
            miss_raw = miss_raw.strip()

            if not add_raw or not miss_raw:
                continue

            # (additional, missing)
            pairs.append((add_raw, miss_raw))

    # Alle Missing-Typen sammeln und indizieren
    missing_types = sorted({m for _, m in pairs})
    missing_id: Dict[str, int] = {m: i for i, m in enumerate(missing_types)}

    # Für jeden Additional-Typ Maske der erlaubten Missing-Typen aufbauen
    allow_mask: Dict[str, int] = {}
    for a, m in pairs:
        bit = 1 << missing_id[m]
        allow_mask[a] = allow_mask.get(a, 0) | bit

    return TransitionMasks(missing_id=missing_id, allow_mask=allow_mask)


def can_cover(
    additional_domains: Iterable[str],
    missing_domains: Iterable[str],
    masks: TransitionMasks,
) -> dict[str, tuple[str, ...]] | bool:
    """
    Prüft, welche additional_domains konkrete missing_domains abdecken können.

    Rückgabe:
      - Dict[additional -> Tuple von missing], wenn es mind. eine gültige Transition gibt.
      - False, wenn gar keine Überführung möglich ist.

    Beispiel dict:
    {
    "extra_domain_A": ("missing1", "missing3"),
    "extra_domain_B": ("missing2",)
    }
    """
    # Indexmap der relevanten Missing-Domains aufbauen
    missing_indices: dict[int, str] = {}
    for m in missing_domains:
        idx = masks.missing_id.get(m)
        if idx is not None:
            missing_indices[idx] = m

    if not missing_indices:
        # keine bekannten Missing-Typen
        return False

    transitions: dict[str, tuple[str, ...]] = {}

    for a in additional_domains:
        allowed_mask = masks.allow_mask.get(a, 0)
        covered: list[str] = []
        for idx, m in missing_indices.items():
            if allowed_mask & (1 << idx):
                covered.append(m)
        if covered:
            transitions[a] = tuple(covered)

    return transitions if transitions else False

import os
from pathlib import Path
import io
import pytest

# Passe den Importpfad ggf. an dein Projekt an
from hmsss.algorithms.search_cross_reference import process_optimized_cutoff


def _write_intermediate(dirpath: Path, hmm_id: str) -> Path:
    """
    Erzeugt eine typische .intermediate_hits-Datei:
    - Kommentarzeile
    - Leere Zeile
    - Drei Treffer mit Scores 49, 50, 60 (Score liegt in Spalte index=7)
    """
    p = dirpath / f"{hmm_id}.intermediate_hits"
    lines = [
        "# domtblout-like header\n",
        "\n",
        "q1\ta\tb\tc\td\te\tf\t49\tX\n",
        "q2\ta\tb\tc\td\te\tf\t50\tX\n",
        "q3\ta\tb\tc\td\te\tf\t60\tX\n",
    ]
    p.write_text("".join(lines), encoding="utf-8")
    return p


def _read_promoted_scores(trusted_path: Path):
    """
    Liest Scores (Spalte 8 / index 7) aus der trusted_hits-Datei,
    ignoriert nicht-tabgetrennte Zeilen (z.B. vorbestehende 'old\n').
    """
    scores = []
    if not trusted_path.exists():
        return scores
    with trusted_path.open("r", encoding="utf-8") as f:
        for line in f:
            if "\t" not in line or not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            scores.append(float(parts[7]))
    return scores


def test_promote_with_specific_threshold(tmp_path: Path):
    cross = tmp_path / "xcheck"
    cross.mkdir()
    hmm_id = "PF00001"

    _write_intermediate(cross, hmm_id)

    # Vorbestehender Inhalt prüfen: sollte im Append-Modus erhalten bleiben
    trusted_path = cross / f"{hmm_id}.trusted_hits"
    trusted_path.write_text("old\n", encoding="utf-8")

    # Threshold 55 → erwartet: Scores 55+ (also 60) und 55 (gleich) → q2=50 NICHT, q3=60 JA, q?=55 (falls vorhanden) JA
    process_optimized_cutoff(hmm_id, str(cross), {hmm_id: 55.0})

    scores = _read_promoted_scores(trusted_path)
    # Es gibt nur 49, 50, 60 → erwartet: nur 60
    assert scores == [60.0], f"Promoted scores mismatch: {scores}"

    # Vorbestandene Zeile blieb erhalten
    raw = trusted_path.read_text(encoding="utf-8").splitlines()
    assert "old" in raw[0]


def test_promote_with_default_threshold_50(tmp_path: Path):
    cross = tmp_path / "xcheck"
    cross.mkdir()
    hmm_id = "PF00002"

    _write_intermediate(cross, hmm_id)

    # Kein Eintrag im Dict → Default 50
    process_optimized_cutoff(hmm_id, str(cross), {})

    trusted_path = cross / f"{hmm_id}.trusted_hits"
    scores = _read_promoted_scores(trusted_path)
    # Erwartet: >=50 → 50 und 60
    assert sorted(scores) == [50.0, 60.0]


def test_all_fallback_threshold_applied(tmp_path: Path):
    """
    Prüft die verbreitete 'all'-Fallback-Konvention:
    Wenn optimized_dict['all'] gesetzt ist und kein spezifischer Eintrag
    für hmm_id vorhanden ist, gilt der 'all'-Wert.

    Da process_optimized_cutoff selbst 'all' nicht interpretiert,
    wird das üblicherweise VOR dem Aufruf aufgelöst. Wir simulieren das hier,
    indem wir das Dict für den Aufruf entsprechend vorbereiten.
    """
    cross = tmp_path / "xcheck"
    cross.mkdir()
    hmm_id = "PF00003"

    _write_intermediate(cross, hmm_id)

    optimized = {"all": 57.0}
    # Üblicher Pre-Processing-Schritt in der Pipeline:
    resolved = {hmm_id: optimized.get(hmm_id, optimized.get("all", 50.0))}

    process_optimized_cutoff(hmm_id, str(cross), resolved)

    trusted_path = cross / f"{hmm_id}.trusted_hits"
    scores = _read_promoted_scores(trusted_path)
    # Threshold 57 → nur 60 wird promotet
    assert scores == [60.0]


def test_missing_intermediate_logs_warning_and_no_trusted(tmp_path: Path, caplog):
    cross = tmp_path / "xcheck"
    cross.mkdir()
    hmm_id = "PF00004"

    # KEINE intermediate-Datei erzeugt
    process_optimized_cutoff(hmm_id, str(cross), {hmm_id: 42.0})

    trusted_path = cross / f"{hmm_id}.trusted_hits"
    assert not trusted_path.exists()

    # Warnung wurde geloggt
    # (je nach Logger-Setup kann das Level variieren; wir prüfen auf Textinhalt)
    assert any(
        "Intermediate or trusted hit file missing" in rec.getMessage()
        for rec in caplog.records
    )

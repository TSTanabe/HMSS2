#!/usr/bin/env python3
import os
import sys
import subprocess
from ete4 import Tree
from typing import Optional, Dict, Any, Set, Iterator, Tuple
import shutil
import random
import csv
from datetime import datetime
import logging
import tempfile

# -------------------------------------------------------------------
# Logger einrichten
# -------------------------------------------------------------------
logger = logging.getLogger("create_test_package")


def setup_logger() -> None:
    """Setzt einen Stream-Logger auf stdout mit DEBUG-Level auf."""
    if logger.handlers:
        # Schon konfiguriert (z.B. wenn als Modul importiert)
        return

    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)


# TESTED on LISC 12.11.2025 for approx. 8000 sequences


def run(cmd: str) -> None:
    logger.info("RUN: %s", cmd)
    subprocess.run(cmd, shell=True, check=True)


# ---------- Helpers für Bait FASTA ----------
def _read_fasta_ids(input_faa: str) -> set[str]:
    ids = set()
    with open(input_faa, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                prot_id = header.rsplit("-", 1)[-1]
                ids.add(prot_id)
    logger.debug(
        "[DEBUG] Parsed %d exclude IDs (first 10: %s)", len(ids), list(ids)[:10]
    )
    return ids


def _grep_blast_lines_with_progress(blast_report: str, protein_type: str) -> list[str]:
    lines, total = [], 0
    logger.info("[INFO] Scanning DIAMOND report for '%s'...", protein_type)
    with open(blast_report, "r") as fh:
        for line in fh:
            total += 1
            if protein_type in line:
                lines.append(line.rstrip("\n"))
            if total % 100000 == 0:
                logger.info(
                    "  processed: %s | matches: %s", f"{total:,}", f"{len(lines):,}"
                )
    logger.info("[DONE] Total lines: %s | matches: %s", f"{total:,}", f"{len(lines):,}")
    return lines


def _collect_subject_ids_from_blast_lines(lines: list[str]) -> set[str]:
    sids = set()
    for ln in lines:
        raw_id = ln.split("\t", 2)[0]  # qseqid (outfmt 6)
        prot_id = raw_id.rsplit("-", 1)[-1]
        sids.add(prot_id)
    logger.debug(
        "[DEBUG] Subject IDs collected: %d (first 10: %s)",
        len(sids),
        list(sids)[:10],
    )
    return sids


def _read_fasta_ids(input_faa: str) -> Set[str]:
    """
    Wie in deinem bestehenden Code: FASTA-IDs lesen und als Suffix nach letztem '-' speichern.
    """
    ids: Set[str] = set()
    with open(input_faa, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                prot_id = header.rsplit("-", 1)[-1]
                ids.add(prot_id)
    logger.debug(
        "[DEBUG] Parsed %d exclude IDs (first 10: %s)", len(ids), list(ids)[:10]
    )
    return ids


def _iter_fasta_records(fp: str):
    """
    Minimaler FASTA-Iterator:
    yield (header_line_with_newline, seq_lines_list, id_suffix_after_last_dash)
    """
    header = None
    seq_lines = []
    current_id = None

    def flush():
        if header is None:
            return None
        return header, seq_lines, current_id

    with open(fp, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                rec = flush()
                if rec is not None:
                    yield rec
                header = line if line.endswith("\n") else line + "\n"
                seq_lines = []
                token = line[1:].strip().split()[0]
                current_id = token.rsplit("-", 1)[-1]
            else:
                seq_lines.append(line)
        rec = flush()
        if rec is not None:
            yield rec


def _ensure_diamond_db(globdb_faa: str, db_dmnd: str, *, threads: int) -> str:
    """
    Stellt sicher, dass die DIAMOND DB existiert und gibt den DB-PREFIX zurück
    (DIAMOND erwartet Prefix ohne '.dmnd' beim --db Parameter).
    """
    db_prefix = db_dmnd[:-5] if db_dmnd.endswith(".dmnd") else db_dmnd
    if os.path.exists(db_dmnd):
        return db_prefix

    cmd = [
        "diamond",
        "makedb",
        "--threads",
        str(threads),
        "--in",
        globdb_faa,
        "-d",
        db_prefix,
    ]
    logger.info("RUN: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return db_prefix


def _diamond_collect_subject_ids(
        *,
        db_prefix: str,
        query_faa: str,
        threads: int,
        evalue: float,
        max_target_seqs: int,
        min_ident: Optional[float],
) -> Set[str]:
    """
    DIAMOND blastp → sammelt nur sseqid (als Suffix nach letztem '-').
    Kein Hit-Parsing außer Spalte 2.
    """
    with tempfile.TemporaryDirectory(prefix="bait_diamond_") as tmpd:
        out_tsv = os.path.join(tmpd, "hits.tsv")

        cmd = [
            "diamond",
            "blastp",
            "--threads",
            str(threads),
            "--db",
            db_prefix,
            "--query",
            query_faa,
            "--out",
            out_tsv,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "--evalue",
            str(evalue),
            "--max-target-seqs",
            str(max_target_seqs),
        ]
        if min_ident is not None:
            cmd.extend(["--id", str(min_ident)])

        logger.info("RUN: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        sids: Set[str] = set()
        with open(out_tsv, "r") as fh:
            for ln in fh:
                if not ln.strip():
                    continue
                # qseqid \t sseqid
                parts = ln.rstrip("\n").split("\t", 2)
                if len(parts) < 2:
                    continue
                sseqid = parts[1]
                sid = sseqid.rsplit("-", 1)[-1]
                sids.add(sid)

        logger.debug(
            "[DEBUG] DIAMOND subject IDs: %d (first 10: %s)", len(sids), list(sids)[:10]
        )
        return sids


def _diamond_collect_subject_ids_with_hi_ident_flag(
        *,
        db_prefix: str,
        query_faa: str,
        threads: int,
        evalue: float,
        max_target_seqs: int,
        hi_ident_threshold: float,
        sensitivity: str = "very-sensitive",
        min_id: float = 25.0,
):
    """
    Single DIAMOND run.
    Liefert:
      candidate_ids    – alle getroffenen subject IDs
      too_similar_ids  – subject IDs mit mind. einem Hit >= hi_ident_threshold
    """
    import tempfile, subprocess

    candidate_ids = set()
    too_similar_ids = set()

    with tempfile.TemporaryDirectory(prefix="diamond_bait_") as tmpd:
        out_tsv = os.path.join(tmpd, "hits.tsv")

        cmd = [
            "diamond",
            "blastp",
            f"--{sensitivity}",
            "--threads",
            str(threads),
            "--db",
            db_prefix,
            "--query",
            query_faa,
            "--out",
            out_tsv,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "--evalue",
            str(evalue),
            "--max-target-seqs",
            str(max_target_seqs),
            "--id",
            str(min_id),
        ]
        subprocess.run(cmd, check=True)

        with open(out_tsv) as fh:
            for ln in fh:
                if not ln.strip():
                    continue
                _, sseqid, pident = ln.rstrip("\n").split("\t")
                sid = sseqid.rsplit("-", 1)[-1]
                candidate_ids.add(sid)
                if float(pident) >= hi_ident_threshold:
                    too_similar_ids.add(sid)

    return candidate_ids, too_similar_ids


def _iter_fasta_simple(fp: str) -> Iterator[Tuple[str, str]]:
    """
    Minimaler FASTA-Iterator: yield (header_without_>, seq_string).
    Nimmt an, dass Sequenzen ggf. über mehrere Zeilen gehen.
    """
    header = None
    seq_chunks = []
    with open(fp, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def _count_fasta_records(fp: str) -> int:
    """Zählt FASTA-Records via '>'-Header."""
    with open(fp, "r") as fh:
        return sum(1 for ln in fh if ln.startswith(">"))


def build_query_faa_for_bait(
        input_faa: str,
        *,
        max_query: int,
        rng: Optional[random.Random] = None,
) -> str:
    """
    Erzeugt eine Query-FASTA für DIAMOND:
    - Wenn input_faa <= max_query Sequenzen enthält: return input_faa
    - Sonst: wählt zufällig max_query Sequenzen und schreibt sie in eine neue FASTA
      im selben Directory wie input_faa. Rückgabe: Pfad zur neuen FASTA.

    Auswahl: uniform über Record-Indizes (0..N-1), ohne Zurücklegen.
    """
    assert os.path.isfile(input_faa), f"Missing FASTA: {input_faa}"
    assert max_query > 0, "max_query must be > 0"

    rng = rng or random

    total = _count_fasta_records(input_faa)
    if total == 0:
        raise ValueError(f"No FASTA records found in: {input_faa}")

    if total <= max_query:
        logger.info(
            "[BAIT][query] Using full query FASTA (%s seqs): %s",
            f"{total:,}",
            input_faa,
        )
        return input_faa

    # Zufällige Indizes bestimmen
    selected = set(rng.sample(range(total), max_query))

    in_dir = os.path.dirname(os.path.abspath(input_faa))
    base = os.path.basename(input_faa)
    stem, ext = os.path.splitext(base)
    out_query = os.path.join(
        in_dir, f"{stem}.query_subsample_{max_query}{ext or '.faa'}"
    )

    # Schreiben (2. Pass)
    written = 0
    seen = 0
    with open(out_query, "w") as out:
        for header, seq in _iter_fasta_simple(input_faa):
            if seen in selected:
                # Header unverändert (ohne '>') wieder schreiben
                out.write(f">{header}\n")
                # Optional: wrap auf 60/80 Zeichen; hier: eine Zeile
                out.write(f"{seq}\n")
                written += 1
                if written >= max_query:
                    break
            seen += 1

    logger.info(
        "[BAIT][query] Subsampled query FASTA written: %s (kept %s of %s)",
        out_query,
        f"{written:,}",
        f"{total:,}",
    )
    return out_query


def _reservoir_sample_fasta_records(
        fasta_path: str,
        *,
        k: int,
        seed: int | None = None,
) -> list[tuple[str, str]]:
    """
    Reservoir sampling über FASTA-Records.
    Return: Liste von (header_line_including_>, seq_with_newlines) mit genau k Records,
            falls mindestens k Records vorhanden sind (sonst weniger).
    1-pass, O(k) RAM.
    """
    rng = random.Random(seed)

    reservoir: list[tuple[str, str]] = []

    def iter_records(fp: str):
        header = None
        seq_lines: list[str] = []
        with open(fp, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_lines)
                    header = line if line.endswith("\n") else line + "\n"
                    seq_lines = []
                else:
                    seq_lines.append(line if line.endswith("\n") else line + "\n")
            if header is not None:
                yield header, "".join(seq_lines)

    for i, (hdr, seq) in enumerate(iter_records(fasta_path), start=1):
        if len(reservoir) < k:
            reservoir.append((hdr, seq))
        else:
            j = rng.randint(1, i)  # 1..i
            if j <= k:
                reservoir[j - 1] = (hdr, seq)

    return reservoir


def make_bait_fasta(
        input_faa: str,
        globdb_faa: str,
        output_bait_faa: str,
        output_dmnd: str,
        *,
        threads: int,
        diamond_db_dmnd: Optional[str] = None,
        evalue: float = 1e-10,
        max_target_seqs: int = 25000,
        hi_ident_threshold: float = 90.0,
        empty_on_no_bait: bool = False,
) -> None:
    assert os.path.isfile(input_faa), f"Missing input_faa: {input_faa}"
    assert os.path.isfile(globdb_faa), f"Missing globdb_faa: {globdb_faa}"
    os.makedirs(os.path.dirname(os.path.abspath(output_bait_faa)) or ".", exist_ok=True)

    query_faa = build_query_faa_for_bait(input_faa=input_faa, max_query=50)

    # Exclude sequences with ids included in input fasta
    exclude_ids = _read_fasta_ids(input_faa)
    max_target_seqs = min(10000, len(exclude_ids) + 5000)  # 500 mehr als

    # DB bestimmen/erzeugen
    if diamond_db_dmnd is None:
        diamond_db_dmnd = globdb_faa + ".dmnd"
    db_prefix = _ensure_diamond_db(globdb_faa, diamond_db_dmnd, threads=threads)

    candidate_ids, too_similar_ids = _diamond_collect_subject_ids_with_hi_ident_flag(
        db_prefix=db_prefix,
        query_faa=query_faa,
        threads=threads,
        evalue=1e-15,  # lockerer als im normalen Lauf
        max_target_seqs=max_target_seqs,  # ggf. höher
        hi_ident_threshold=hi_ident_threshold,  # muss gesetzt sein, wird hier aber nicht zwingend genutzt
        sensitivity="very-sensitive",
        min_id=20.0,  # lockerer als 25 (oder noch niedriger)
    )

    bait_ids = candidate_ids - exclude_ids - too_similar_ids

    logger.info(
        "[BAIT] After exclusion: %s in the baid(removed exclude=%s, too_similar=%s)",
        f"{len(bait_ids):,}",
        f"{len(candidate_ids & exclude_ids):,}",
        f"{len(candidate_ids & too_similar_ids):,}",
    )

    if not bait_ids:
        logger.warning(
            "[BAIT] No bait_ids after filtering → fallback: random %d sequences from globdb",
            1000,
        )

        fallback_k = 1000
        oversample_k = (
            5000  # > fallback_k, damit exclude-filter nicht zu viele rauswirft
        )
        seed = None

        def _norm_id_from_header(hdr_line: str) -> str:
            # hdr_line enthält führendes '>'
            token = hdr_line[1:].strip().split()[0]
            return token.rsplit("-", 1)[-1]  # exakt wie im globdb-scan

        # 1) Ziehe zufällige Records (Reservoir Sampling)
        sampled_records = _reservoir_sample_fasta_records(
            globdb_faa, k=oversample_k, seed=seed
        )

        if not sampled_records:
            logger.warning(
                "[BAIT][fallback-random] Fatal: globdb_faa contains no sequences: %s",
                globdb_faa,
            )
            open(output_bait_faa, "w").close()
            sys.exit(1)
            return

        # 2) Filtere: nicht im exclude, keine Duplikate; cappe auf fallback_k
        kept_records: list[tuple[str, str]] = []
        kept_ids: set[str] = set()

        excluded_hits = 0
        dup_hits = 0

        for hdr, seq in sampled_records:
            cid = _norm_id_from_header(hdr)
            if cid in exclude_ids:
                excluded_hits += 1
                continue
            if cid in kept_ids:
                dup_hits += 1
                continue
            kept_ids.add(cid)
            kept_records.append((hdr, seq))
            if len(kept_records) >= fallback_k:
                break

        if not kept_records:
            logger.warning(
                "[BAIT][fallback-random] Fatal: after filtering exclude_ids, still empty bait FASTA: %s",
                output_bait_faa,
            )
            open(output_bait_faa, "w").close()
            sys.exit(1)
            return

        if len(kept_records) < fallback_k:
            logger.warning(
                "[BAIT][fallback-random] Only %s/%s sequences kept (excluded=%s, dups=%s). "
                "Consider increasing oversample_k.",
                f"{len(kept_records):,}",
                f"{fallback_k:,}",
                f"{excluded_hits:,}",
                f"{dup_hits:,}",
            )

        # 3) Schreibe bait FASTA
        os.makedirs(
            os.path.dirname(os.path.abspath(output_bait_faa)) or ".", exist_ok=True
        )
        with open(output_bait_faa, "w") as out:
            for hdr, seq in kept_records:
                out.write(hdr if hdr.endswith("\n") else hdr + "\n")
                out.write(seq)
                if not seq.endswith("\n"):
                    out.write("\n")

        logger.info(
            "[DONE][fallback-random] Wrote %s random sequences to bait FASTA: %s (excluded=%s, dups=%s, sampled=%s)",
            f"{len(kept_records):,}",
            output_bait_faa,
            f"{excluded_hits:,}",
            f"{dup_hits:,}",
            f"{len(sampled_records):,}",
        )

        # 4) DIAMOND-DB für bait bauen (wie gehabt)
        if os.path.exists(output_dmnd):
            logger.info("→ Skip DIAMOND DB (bait, exists)")
        else:
            dmnd_prefix = output_dmnd[:-5]  # assumes output_dmnd endswith ".dmnd"
            run(
                f"diamond makedb --threads 4 --in '{output_bait_faa}' -d '{dmnd_prefix}'"
            )

        return

    # GlobDB scannen und schreiben
    written_ids: Set[str] = set()
    scanned = 0
    written = 0

    logger.info("[BAIT] Scanning globdb FASTA for bait sequences: %s", globdb_faa)
    with open(globdb_faa, "r") as inp, open(output_bait_faa, "w") as out:
        current_header, current_seq, current_id = None, [], None

        def flush():
            nonlocal written
            if current_header is None:
                return
            if current_id in bait_ids and current_id not in written_ids:
                out.write(current_header)
                out.write("".join(current_seq))
                if not current_seq or not current_seq[-1].endswith("\n"):
                    out.write("\n")
                written_ids.add(current_id)
                written += 1

        for line in inp:
            if line.startswith(">"):
                flush()
                current_header = line if line.endswith("\n") else line + "\n"
                current_seq = []
                token = line[1:].strip().split()[0]
                current_id = token.rsplit("-", 1)[-1]
                scanned += 1
                if scanned % 500000 == 0:
                    logger.info(
                        "  processed %s sequences | written %s",
                        f"{scanned:,}",
                        f"{written:,}",
                    )
            else:
                current_seq.append(line)
        flush()

    # Make the bait dmnd
    if os.path.exists(output_dmnd):
        logger.info("→ Skip DIAMOND DB (bait, exists)")
    else:
        dmnd_prefix = output_dmnd[:-5]
        run(f"diamond makedb --threads 4 --in '{output_bait_faa}' -d '{dmnd_prefix}'")
    logger.info(
        "[DONE] Total scanned: %s | Bait sequences written: %s",
        f"{scanned:,}",
        f"{written:,}",
    )


# ---------- Pipeline-Schritte ----------
def step3_align_mafft(input_faa: str, output_aln: str, threads: int) -> None:
    run(f"mafft --auto --thread {threads} {input_faa} > {output_aln}")


def step4_trimal(input_aln: str, output_trimmed: str) -> None:
    run(f"trimal -in {input_aln} -out {output_trimmed} -gt 0.2")


def step5_build_tree(
        input_alignment: str, output_tree: str, veryfasttree_bin: str
) -> None:
    cmd = f"{veryfasttree_bin} -out {output_tree} {input_alignment}"
    logger.info("RUN: %s", cmd)
    subprocess.run(cmd, shell=True, check=True)


def step6_midpoint_root(input_tree: str, output_rooted_tree: str) -> None:
    t = Tree(input_tree)
    t.set_outgroup(t.get_midpoint_outgroup())
    t.write(outfile=output_rooted_tree)
    logger.info("Midpoint rooted tree written to %s", output_rooted_tree)


# ---------- GraftM-Package bauen ----------
def step9_create_gpkg(
        output_gpkg: str,
        seqs_faa: str,
        trimmed_aln: str,
        rooted_tree: str,
        taxonomy_tsv: str,
        threads: int = 12,
) -> None:
    """
    Erstellt ein GraftM gpkg Paket.

    Erwartet:
      seqs_faa      = unalignierte Referenzsequenzen (FASTA)
      trimmed_aln   = getrimmtes Alignment (FASTA/PHYLIP)
      rooted_tree   = gerooteter Newick-Baum
      taxonomy_tsv  = 2-Spalten-Taxonomie-Datei (ID<TAB>Taxon)
    """
    if os.path.exists(output_gpkg):
        logger.info("→ Skip GraftM gpkg (exists): %s", output_gpkg)
        return

    assert os.path.exists(seqs_faa), f"Missing sequences file: {seqs_faa}"
    assert os.path.exists(trimmed_aln), f"Missing alignment: {trimmed_aln}"
    assert os.path.exists(rooted_tree), f"Missing rooted tree: {rooted_tree}"
    assert os.path.exists(taxonomy_tsv), f"Missing taxonomy: {taxonomy_tsv}"

    cmd = (
        f"graftM create "
        f"--output '{output_gpkg}' "
        f"--sequences '{seqs_faa}' "
        f"--alignment '{trimmed_aln}' "
        f"--rerooted_tree '{rooted_tree}' "
        f"--taxonomy '{taxonomy_tsv}' "
        f"--threads {threads}"
    )
    run(cmd)
    logger.info("✔ Created gpkg: %s", output_gpkg)


# ---------- Step 10 – Bait-Artefakte ins gpkg kopieren ----------
def step10_attach_bait_artifacts(
        out_dir: str,
        protein_type: str,
        bait_faa: str,
        bait_dmnd: str,
        target_dir: str,
) -> Optional[str]:
    """
    Kopiert (falls vorhanden) 6_bait.faa und 7_bait_db.dmnd in den Ordner
    <out_dir>/<protein_type>.gpkg/.
    Legt KEIN gpkg-Verzeichnis an, damit 'graftM create' nicht irrtümlich
    übersprungen wird.

    Idempotent: wiederholte Aufrufe sind ok.

    Returns
    -------
    Optional[str]
        Pfad zur kopierten DIAMOND-Datenbank (.dmnd) im Zielverzeichnis,
        oder None, falls keine Datei kopiert wurde.
    """
    if not os.path.isdir(target_dir):
        logger.info(
            "ℹ gpkg-Verzeichnis noch nicht vorhanden, überspringe Copy: %s", target_dir
        )
        return None

    # 6_bait.faa kopieren
    if os.path.exists(bait_faa):
        shutil.copy2(bait_faa, os.path.join(target_dir, "6_bait.faa"))
        logger.info(" Bait FASTA → %s", target_dir)
    else:
        logger.warning("⚠ Keine bait FASTA gefunden: %s", bait_faa)

    # 7_bait_db.dmnd kopieren
    if os.path.exists(bait_dmnd):
        target_dmnd = os.path.join(target_dir, "7_bait_db.dmnd")
        shutil.copy2(bait_dmnd, target_dmnd)
        logger.info(" DIAMOND DB → %s", target_dmnd)
        return target_dmnd
    else:
        logger.warning("⚠ Keine DIAMOND DB gefunden: %s", bait_dmnd)
        return None


def step_make_mock_fragments(
        input_faa: str,
        out_dir: str,
        output_name: str,
        *,
        min_len: int = 40,
        max_len: int = 60,
        min_step: int = 15,
        max_step: int = 45,
        max_sequences: int = 100000,
) -> int:
    """
    Erzeugt Mock-Fragmente aus Protein-FASTA.
    - Output: out_dir/<output_name>
    - Fragmentnamen: <orig_name>_<1..N> (pro Protein fortlaufend)
    - Return: Gesamtzahl erzeugter Fragmente

    Regeln:
      • L <= 80 aa: genau 5 Fragmente
          1: gesamte Sequenz
          2–5: N- und C-terminale Teilstücke (ohne Padding, nur echte Substrings)
      • L > 80 aa: Sliding-Fenster mit jitter (60–80 aa) über (nahezu) gesamte Länge

    Zusätze:
      • Wenn die Ausgabedatei bereits existiert, wird NICHT neu generiert; es werden
        nur die Fragmente gezählt (Anzahl '>'-Header) und diese Zahl zurückgegeben.
      • Verarbeite höchstens `max_sequences` zufällig ausgewählte Eingabesequenzen.
    """
    assert os.path.isfile(input_faa), f"FASTA not found: {input_faa}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, output_name)

    # Falls schon vorhanden: nur zählen und zurückgeben
    if os.path.isfile(out_path):
        with open(out_path, "r") as fh:
            existing = sum(1 for line in fh if line.startswith(">"))
        logger.info(
            "[Mock] Output existiert bereits → überspringe Generierung. Fragmente: %d",
            existing,
        )
        return existing

    # Hilfsfunktionen
    def _iter_fasta(fp: str):
        name, seq = None, []
        with open(fp, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(seq)
                    name = line[1:].strip().split()[0]
                    seq = []
                else:
                    seq.append(line.strip())
            if name is not None:
                yield name, "".join(seq)

    # Pass 1: Gesamtzahl der Sequenzen zählen
    with open(input_faa, "r") as fh:
        total_sequences = sum(1 for l in fh if l.startswith(">"))
    if total_sequences == 0:
        logger.warning("[Mock] Keine Sequenzen gefunden.")
        return 0

    # Stichprobe bestimmen (Indices 0..N-1)
    k = min(max_sequences, total_sequences)
    if k < total_sequences:
        selected_indices = set(random.sample(range(total_sequences), k))
        sampled = True
    else:
        selected_indices = None
        sampled = False

    # Pass 2: iterieren und nur ausgewählte Sequenzen verarbeiten
    processed = 0  # wie viele Eingabesequenzen wurden angesehen
    taken = 0  # wie viele davon verarbeitet (bei Sampling ≤ k)
    total_frags = 0

    with open(out_path, "w") as out_fh:
        last_pct_printed = -10  # VOR der for-Schleife initialisieren

        for base_name, aa in _iter_fasta(input_faa):
            processed += 1

            if selected_indices is not None and (processed - 1) not in selected_indices:
                pct = (processed / total_sequences) * 100.0
                pct_bucket = int(pct // 10) * 10
                if pct_bucket > 100:
                    pct_bucket = 100

                if pct_bucket >= last_pct_printed + 10:
                    last_pct_printed = pct_bucket
                    sys.stdout.write(
                        f"\r[Mock] Seen {processed}/{total_sequences} ({pct_bucket:3d}%) | kept {taken}/{k}"
                    )
                    sys.stdout.flush()

                continue

            taken += 1
            L = len(aa)
            frag_idx = 0

            if L <= max_len:
                # 1) gesamte Sequenz
                frag_idx += 1
                total_frags += 1
                out_fh.write(f">{base_name}_{frag_idx}\n{aa}\n")

                # 2)–5) N- und C-terminale Substrings (ohne Padding)
                n_len1 = max(1, min(L - 1, int(round(0.50 * L))))
                n_len2 = max(1, min(L - 1, int(round(0.70 * L))))
                c_len1 = n_len1
                c_len2 = n_len2

                # N-terminal
                for nlen in (n_len1, n_len2):
                    frag_idx += 1
                    total_frags += 1
                    out_fh.write(f">{base_name}_{frag_idx}\n{aa[:nlen]}\n")

                # C-terminal
                for clen in (c_len1, c_len2):
                    frag_idx += 1
                    total_frags += 1
                    out_fh.write(f">{base_name}_{frag_idx}\n{aa[L - clen:]}\n")

            else:
                # Sliding-Fenster mit jitter
                pos = 0
                while pos < L:
                    frag_len = random.randint(min_len, max_len)
                    end = pos + frag_len
                    if end > L:
                        end = L
                        pos = max(0, end - frag_len)
                    frag = aa[pos:end]
                    if frag:
                        frag_idx += 1
                        total_frags += 1
                        out_fh.write(f">{base_name}_{frag_idx}\n{frag}\n")
                    step = random.randint(min_step, max_step)
                    pos += step
                    # Abschluss-Guard, falls Ende noch nicht abgedeckt war
                    if pos >= L and end < L:
                        tail_start = max(0, L - random.randint(min_len, max_len))
                        tail = aa[tail_start:L]
                        if tail:
                            frag_idx += 1
                            total_frags += 1
                            out_fh.write(f">{base_name}_{frag_idx}\n{tail}\n")
                        break

            # Fortschritt live (eine Zeile)
            pct = (processed / total_sequences) * 100
            base_info = f"kept {taken}/{k}" if sampled else f"kept {taken}"
            sys.stdout.write(
                f"\r[Mock] Seen {processed}/{total_sequences} ({pct:.0f}%) | {base_info}"
            )
            sys.stdout.flush()

            # Abbruch, wenn Sampling erreicht
            if sampled and taken >= k:
                break

    sys.stdout.write("\n")
    logger.info("✔ Mock-Fragmente geschrieben: %s  (n=%d)", out_path, total_frags)
    if sampled:
        logger.info(
            "[Mock] Stichprobe: %d von %d Sequenzen verarbeitet (max_sequences=%d)",
            k,
            total_sequences,
            max_sequences,
        )
    return total_frags


def step11_graft_with_package(
        gpkg_dir: str,
        input_faa: str,
        search_dmnd: str,
        bait_dmnd: str,
        threads: int = 4,
        result_dir: str = "graft_results",
) -> Optional[str]:
    """
    Führt 'graftM graft' mit dem erzeugten gpkg gegen eine Protein-FASTA aus.

    FAILSAFE:
    - Falls graftM mit einem Fehlercode endet, wird ein ERROR geloggt
      und None zurückgegeben, anstatt das gesamte Script zu beenden.
    """
    if os.path.isdir(result_dir):
        logger.info("[SKIP] Result directory already exists: %s", result_dir)
        return result_dir
    else:
        os.makedirs(result_dir, exist_ok=True)

    cmd = [
        "graftM",
        "graft",
        "--graftm_package",
        gpkg_dir,
        "--forward",
        input_faa,
        "--output_directory",
        result_dir,
        "--threads",
        str(threads),
        "--decoy_database",
        bait_dmnd,
        "--search_diamond_file",
        "Q",
        "--search_method",
        "hmmsearch+diamond",
        "--force",
    ]

    logger.info(
        "[INFO] Running graftM on %s using package %s",
        os.path.basename(input_faa),
        os.path.basename(gpkg_dir),
    )
    logger.info("RUN: %s", " ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(
            "graftM run FAILED for %s (result_dir=%s): %s",
            input_faa,
            result_dir,
            e,
        )
        return None

    logger.info("✔ graftM run abgeschlossen: %s", result_dir)
    return result_dir


def safe_move(src: str, dst: str) -> str:
    """
    Kopiert src nach dst.
    dst kann ein Verzeichnis ODER ein vollständiger Dateipfad sein.
    """
    if os.path.isdir(dst):
        dst_path = os.path.join(dst, os.path.basename(src))
    else:
        dst_path = dst

    if os.path.exists(dst_path):
        logger.warning("[WARN] Überschreibe bestehende Datei: %s", dst_path)
    else:
        logger.info("[COPY] %s → %s", src, dst_path)

    shutil.copy2(src, dst_path)
    return dst_path


# ---------------------------------------------------------------------


def find_read_tax_file(root: str) -> Optional[str]:
    """
    Durchsucht root rekursiv nach einer Datei, die auf 'read_tax.tsv' endet.
    Gibt den vollständigen Pfad zurück oder None, falls keine gefunden wird.
    """
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.endswith("read_tax.tsv"):
                return os.path.join(dirpath, f)
    return None


RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX_MAP = {
    "k": "kingdom",
    "d": "kingdom",  # domain -> kingdom-Ebene
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def parse_taxonomy_string(tax_str: str) -> Dict[str, Optional[str]]:
    """
    Zerlegt einen Taxonomie-String aus GraftM in einzelne Ränge.

    Unterstützt z.B.:
      - 'k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; ...'
      - 'Bacteria; Proteobacteria; Gammaproteobacteria; ...'
    """
    result = {r: None for r in RANK_ORDER}

    if not tax_str:
        return result

    parts = [p.strip() for p in tax_str.split(";") if p.strip()]
    for i, part in enumerate(parts):
        # Fall 1: hat Präfix wie k__, p__, ...
        if "__" in part:
            prefix, name = part.split("__", 1)
            prefix = prefix.strip().lower()
            name = name.strip() or None

            key = None
            if prefix:
                key = RANK_PREFIX_MAP.get(prefix[0])

            if key and name:
                result[key] = name

        else:
            # Fall 2: keine Präfixe, dann nach Position zuordnen
            if i < len(RANK_ORDER):
                rank = RANK_ORDER[i]
                name = part.strip() or None
                if name:
                    result[rank] = name

    return result


def parse_graftm_read_tax(path: str) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Parst eine GraftM read_tax- oder Taxonomie-TSV-Datei.

    Unterstützte Formate:
    - Ohne Header (typische GraftM read_tax.tsv):
        <read_id>\t<tax_string>
    - Mit Header:
        id\ttaxonomy
        <read_id>\t<tax_string>

    Rückgabe:
        { read_id: {kingdom, phylum, ..., species}, ... }
    """
    result: Dict[str, Dict[str, Optional[str]]] = {}

    with open(path, encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        return result

    # Prüfen, ob erste Zeile wie ein Header aussieht
    first_cols = lines[0].split("\t")
    is_header = False
    id_idx = 0
    tax_idx = 1 if len(first_cols) > 1 else None

    lower_first = [c.strip().lower() for c in first_cols]
    if any(
            x in lower_first
            for x in ("id", "read", "read_id", "sequence", "seqid", "seq_id")
    ) or any(
        x in lower_first
        for x in ("taxonomy", "classification", "graftm_taxonomy", "taxon")
    ):
        is_header = True
        for i, col in enumerate(lower_first):
            if any(
                    x in col
                    for x in ("id", "read", "read_id", "sequence", "seqid", "seq_id")
            ):
                id_idx = i
                break
        tax_idx = None
        for i, col in enumerate(lower_first):
            if col in ("taxonomy", "classification", "graftm_taxonomy", "taxon"):
                tax_idx = i
                break
        if tax_idx is None and len(first_cols) > 1:
            tax_idx = 1
        data_lines = lines[1:]
    else:
        data_lines = lines

    if tax_idx is None:
        raise ValueError(
            f"Keine Taxonomie-Spalte in {path} erkannt. Erste Zeile: {lines[0]!r}"
        )

    for ln in data_lines:
        cols = ln.split("\t")
        if len(cols) <= max(id_idx, tax_idx):
            continue
        read_id = cols[id_idx].strip()
        tax_str = cols[tax_idx].strip()
        result[read_id] = parse_taxonomy_string(tax_str)

    return result


RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def evaluate_tax_assignments(
        reference: Dict[str, Dict[str, Any]],
        test: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    n_test_keys = len(test)
    n_correct_lowest_rank = 0

    correct_per_rank = {rank: 0 for rank in RANK_ORDER}
    assigned_per_rank = {rank: 0 for rank in RANK_ORDER}
    wrong_per_rank = {rank: 0 for rank in RANK_ORDER}
    not_assigned_per_rank = {rank: 0 for rank in RANK_ORDER}

    for test_key, test_lineage in test.items():
        ref_key = test_key.rsplit("_", 1)[0]
        ref_lineage = reference.get(ref_key)
        if ref_lineage is None:
            continue

        # ---- deepest-level correctness ----
        lowest_rank, lowest_value = None, None
        for rank in reversed(RANK_ORDER):
            val = test_lineage.get(rank)
            if val is not None:
                lowest_rank = rank
                lowest_value = val
                break

        if lowest_rank is not None:
            ref_val = ref_lineage.get(lowest_rank)
            if ref_val is not None and ref_val == lowest_value:
                n_correct_lowest_rank += 1

        # ---- per-rank stats ----
        for rank in RANK_ORDER:
            test_val = test_lineage.get(rank)
            ref_val = ref_lineage.get(rank)

            if test_val is None:
                not_assigned_per_rank[rank] += 1
            else:
                assigned_per_rank[rank] += 1

            if test_val is not None and ref_val is not None and test_val == ref_val:
                correct_per_rank[rank] += 1

            if test_val is not None and ref_val is not None and test_val != ref_val:
                wrong_per_rank[rank] += 1

    return {
        "n_test_keys": n_test_keys,
        "n_correct_lowest_rank": n_correct_lowest_rank,
        "correct_per_rank": correct_per_rank,
        "assigned_per_rank": assigned_per_rank,
        "wrong_per_rank": wrong_per_rank,
        "not_assigned_per_rank": not_assigned_per_rank,
    }


# ---------------------------------------------------------------------


def write_eval_txt(
        gpkg_dir: str,
        protein_type: str,
        TP_frags: int,
        TN_frags: int,
        eval_stats: Dict[str, Any],
        fp_test_keys: int,
) -> str:
    """
    Schreibt alle Evaluationswerte in ein einfaches maschinenlesbares Textfile.
    Format:
        key = value
    """
    import math

    n_test_keys = eval_stats["n_test_keys"]
    tp_lowest_rank = eval_stats["n_correct_lowest_rank"]

    correct_per_rank = eval_stats["correct_per_rank"]
    assigned_per_rank = eval_stats.get("assigned_per_rank", {})
    wrong_per_rank = eval_stats.get("wrong_per_rank", {})
    not_assigned_per_rank = eval_stats.get("not_assigned_per_rank", {})

    # Overall Detection Confusion
    TP_overall = n_test_keys
    FP_overall = fp_test_keys
    FN_overall = TP_frags - TP_overall
    TN_overall = TN_frags - FP_overall

    # Placement Confusion (tiefste Ebene)
    TP_place = tp_lowest_rank
    FP_place = fp_test_keys
    FN_place = TP_frags - TP_place
    TN_place = TN_frags - FP_place

    def _compute_binary_metrics(TP: int, FP: int, FN: int, TN: int):
        tp = float(TP)
        fp = float(FP)
        fn = float(FN)
        tn = float(TN)

        tpr = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        tnr = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
        if math.isnan(tpr) or math.isnan(tnr):
            bal_acc = float("nan")
        else:
            bal_acc = 0.5 * (tpr + tnr)

        denom_f1 = 2.0 * tp + fp + fn
        if denom_f1 > 0:
            f1 = 2.0 * tp / denom_f1
        else:
            f1 = float("nan")

        denom_mcc = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        if denom_mcc > 0:
            mcc = (tp * tn - fp * fn) / math.sqrt(denom_mcc)
        else:
            mcc = float("nan")

        return bal_acc, f1, mcc

    overall_bal_acc, overall_f1, overall_mcc = _compute_binary_metrics(
        TP_overall, FP_overall, FN_overall, TN_overall
    )
    placement_bal_acc, placement_f1, placement_mcc = _compute_binary_metrics(
        TP_place, FP_place, FN_place, TN_place
    )

    out_path = os.path.join(gpkg_dir, "evaluation_stats.txt")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"protein_type = {protein_type}\n")
        fh.write(f"TP_frags = {TP_frags}\n")
        fh.write(f"TN_frags = {TN_frags}\n\n")

        fh.write("# Fixpunkte aus evaluate_tax_assignments\n")
        fh.write(f"n_test_keys = {n_test_keys}\n")
        fh.write(f"fp_test_keys = {fp_test_keys}\n")
        fh.write(f"tp_lowest_rank = {tp_lowest_rank}\n\n")

        fh.write("# Pro Taxonomieebene: Anzahl korrekt zugeordneter Reads\n")
        for rank, val in correct_per_rank.items():
            fh.write(f"accuracy_{rank} = {val}\n")
        fh.write("\n")

        fh.write("# Pro Taxonomieebene: Zuordnungsstatistiken\n")
        for rank in RANK_ORDER:
            corr = correct_per_rank.get(rank, 0)
            assg = assigned_per_rank.get(rank, 0)
            wrong = wrong_per_rank.get(rank, 0)
            not_assg = not_assigned_per_rank.get(rank, 0)

            fh.write(f"assigned_{rank} = {assg}\n")
            fh.write(f"correct_{rank} = {corr}\n")
            fh.write(f"wrong_{rank} = {wrong}\n")
            fh.write(f"not_assigned_{rank} = {not_assg}\n")
        fh.write("\n")

        fh.write("# Overall Confusion Matrix (Detection)\n")
        fh.write(f"overall_TP = {TP_overall}\n")
        fh.write(f"overall_FP = {FP_overall}\n")
        fh.write(f"overall_FN = {FN_overall}\n")
        fh.write(f"overall_TN = {TN_overall}\n")
        fh.write(f"overall_balanced_accuracy = {overall_bal_acc:.6f}\n")
        fh.write(f"overall_F1 = {overall_f1:.6f}\n")
        fh.write(f"overall_MCC = {overall_mcc:.6f}\n\n")

        fh.write("# Placement Confusion Matrix (tiefste Ebene)\n")
        fh.write(f"placement_TP = {TP_place}\n")
        fh.write(f"placement_FP = {FP_place}\n")
        fh.write(f"placement_FN = {FN_place}\n")
        fh.write(f"placement_TN = {TN_place}\n")
        fh.write(f"placement_balanced_accuracy = {placement_bal_acc:.6f}\n")
        fh.write(f"placement_F1 = {placement_f1:.6f}\n")
        fh.write(f"placement_MCC = {placement_mcc:.6f}\n")

    logger.info("[EVAL] Stats → %s", out_path)
    return out_path


# ---------- Verarbeitung einer einzelnen FASTA-Datei ----------
def process_protein_fasta(
        input_faa: str,
        out_dir: str,
        threads: int,
        veryfasttree_bin: str,
        globdb_faa: str,
        taxonomy_tsv: str,
) -> None:
    protein_type = os.path.splitext(os.path.basename(input_faa))[0]
    os.makedirs(out_dir, exist_ok=True)

    in_faa = input_faa
    out_aln = os.path.join(out_dir, "2_aligned.faa")
    out_trim = os.path.join(out_dir, "3_trimmed.faa")
    out_tree = os.path.join(out_dir, "4_tree.nwk")
    out_root = os.path.join(out_dir, "5_tree_rooted.nwk")
    out_bait = os.path.join(out_dir, "6_bait.faa")
    out_dmnd = os.path.join(out_dir, "7_bait_db.dmnd")
    out_refdmnd = os.path.join(out_dir, "1_refseq.dmnd")
    out_refdmnd_single_char = os.path.join(out_dir, "Q")
    out_gpkg = os.path.join(out_dir, f"{protein_type}.gpkg")

    out_graft_test_positive = os.path.join(out_dir, "11_graft_positive")
    out_graft_test_negative = os.path.join(out_dir, "11_graft_negative")

    logger.info("\n=== %s ===", protein_type)
    logger.info("Input: %s  | Outdir: %s", input_faa, out_dir)

    if os.path.exists(out_aln):
        logger.info("→ Skip MAFFT (exists)")
    else:
        step3_align_mafft(input_faa, out_aln, threads)

    if os.path.exists(out_trim):
        logger.info("→ Skip trimAl (exists)")
    else:
        step4_trimal(out_aln, out_trim)

    if os.path.exists(out_tree):
        logger.info("→ Skip VeryFastTree (exists)")
    else:
        step5_build_tree(out_trim, out_tree, veryfasttree_bin)

    if os.path.exists(out_root):
        logger.info("→ Skip midpoint rooting (exists)")
    else:
        step6_midpoint_root(out_tree, out_root)

    if os.path.exists(out_bait):
        logger.info("→ Skip bait FASTA (exists)")
    else:
        make_bait_fasta(
            input_faa=input_faa,
            globdb_faa=globdb_faa,
            output_bait_faa=out_bait,
            output_dmnd=out_dmnd,
            threads=4,
        )

    # Make the bait dmnd
    if os.path.exists(out_dmnd):
        logger.info("→ Skip DIAMOND DB (bait, exists)")
    else:
        dmnd_prefix = out_dmnd[:-5]
        run(f"diamond makedb --threads 4 --in '{out_bait}' -d '{dmnd_prefix}'")

    # Make the refseq dmnd
    if os.path.exists(out_refdmnd) or os.path.exists(out_refdmnd_single_char):
        logger.info("→ Skip DIAMOND DB (refseq, exists)")
    else:
        dmnd_prefix = out_refdmnd[:-5]
        run(f"diamond makedb --threads 4 --in '{in_faa}' -d '{dmnd_prefix}'")
        os.rename(out_refdmnd, out_refdmnd_single_char)  # das sollte nun Q heißen

    # GPKG Paket bauen
    step9_create_gpkg(
        output_gpkg=out_gpkg,
        seqs_faa=input_faa,
        trimmed_aln=out_trim,
        rooted_tree=out_root,
        taxonomy_tsv=taxonomy_tsv,
        threads=14,
    )

    # Verschiebe die bait.faa und .dmnd in das gpkg paket
    safe_move(out_bait, os.path.join(out_gpkg, "decoy_database.faa"))
    safe_move(
        in_faa, os.path.join(out_gpkg, "refseq_database.faa")
    )  # moves refseq to gpkg package
    # safe_move(out_refdmnd_single_char, out_gpkg)

    cwd = os.getcwd()  # Hier wird ins working directory kopiert, sonst scheitert graftM
    safe_move(out_refdmnd_single_char, cwd)

    gpkg_dmnd = out_dmnd
    gpkg_refdmnd = os.path.join(cwd, "Q")

    # GPKG Paket testen: Baue TP und TN Sets
    TP_name = "8_TP_reads.faa"
    TP_frags = step_make_mock_fragments(
        input_faa=input_faa,
        out_dir=out_dir,
        output_name=TP_name,
        max_sequences=50000,
    )
    TN_name = "9_TN_reads.faa"
    TN_frags = step_make_mock_fragments(
        input_faa=out_bait,
        out_dir=out_dir,
        output_name=TN_name,
        max_sequences=200000,
    )

    # --- graftM-Runs (TP/TN) mit Failover-Logik ---
    tp_input = os.path.join(out_dir, TP_name)
    graft_out_dir_positive = step11_graft_with_package(
        gpkg_dir=out_gpkg,
        input_faa=tp_input,
        search_dmnd=gpkg_refdmnd,
        bait_dmnd=gpkg_dmnd,
        threads=14,
        result_dir=out_graft_test_positive,
    )

    tn_input = os.path.join(out_dir, TN_name)
    graft_out_dir_negative = step11_graft_with_package(
        gpkg_dir=out_gpkg,
        input_faa=tn_input,
        search_dmnd=gpkg_refdmnd,
        bait_dmnd=gpkg_dmnd,
        threads=14,
        result_dir=out_graft_test_negative,
    )

    logger.info("[CONFUSION MATRIX]")

    # Taxonomie einlesen (ohne die geht die rank-basierte Evaluation nicht sinnvoll)
    try:
        taxonomy = parse_graftm_read_tax(taxonomy_tsv)
    except Exception as e:
        logger.error(
            "[EVAL] Konnte Taxonomie-TSV %s nicht parsen: %s. "
            "Überspringe rank-basierte Evaluation für %s.",
            taxonomy_tsv,
            e,
            protein_type,
        )
        logger.info("✔ Fertig (ohne evaluation_stats.txt): %s", protein_type)
        return

    # -----------------------------
    # TP-Ergebnisse (positive Reads)
    # -----------------------------
    positive_read_parsed: Dict[str, Dict[str, Any]] = {}

    if graft_out_dir_positive is None:
        logger.error(
            "[EVAL] Positive graftM-Ausführung fehlgeschlagen – keine TP-Assignments vorhanden."
        )
    else:
        read_tax_path = find_read_tax_file(graft_out_dir_positive)
        if read_tax_path is None:
            logger.error(
                "[EVAL] Keine read_tax.tsv im positiven graftM-Output gefunden (%s). "
                "Keine TP-Assignments.",
                graft_out_dir_positive,
            )
        else:
            try:
                positive_read_parsed = parse_graftm_read_tax(read_tax_path)
            except Exception as e:
                logger.error(
                    "[EVAL] Fehler beim Parsen von %s: %s – keine TP-Assignments.",
                    read_tax_path,
                    e,
                )

    # Falls positive_read_parsed leer ist, wird evaluate_tax_assignments
    # trotzdem aufgerufen – dann sind n_test_keys=0 etc.
    eval_stats = evaluate_tax_assignments(
        reference=taxonomy,
        test=positive_read_parsed,
    )

    # -----------------------------
    # TN-Ergebnisse (negative Reads)
    # -----------------------------
    negative_read_parsed: Dict[str, Dict[str, Any]] = {}

    if graft_out_dir_negative is None:
        logger.error(
            "[EVAL] Negative graftM-Ausführung fehlgeschlagen – keine TN/FP-Assignments vorhanden."
        )
    else:
        neg_read_tax_path = find_read_tax_file(graft_out_dir_negative)
        if neg_read_tax_path is None:
            logger.error(
                "[EVAL] Keine read_tax.tsv im negativen graftM-Output gefunden (%s). "
                "Keine TN/FP-Assignments.",
                graft_out_dir_negative,
            )
        else:
            try:
                negative_read_parsed = parse_graftm_read_tax(neg_read_tax_path)
            except Exception as e:
                logger.error(
                    "[EVAL] Fehler beim Parsen von %s: %s – keine TN/FP-Assignments.",
                    neg_read_tax_path,
                    e,
                )

    # Anzahl der als positiv klassifizierten Reads aus dem negativen Datensatz
    # (= False Positives). Wenn negative_read_parsed leer ist, ist fp_test_keys=0.
    fp_test_keys = len(negative_read_parsed)

    # evaluation_stats.txt immer schreiben – basierend auf den vorhandenen Daten
    write_eval_txt(
        gpkg_dir=out_gpkg,
        protein_type=protein_type,
        TP_frags=TP_frags,
        TN_frags=TN_frags,
        eval_stats=eval_stats,
        fp_test_keys=fp_test_keys,
    )

    logger.info(
        "✔ Fertig (Evaluation geschrieben, TP-reads=%d, TN-reads=%d, FP_reads=%d): %s",
        eval_stats["n_test_keys"],
        TN_frags,
        fp_test_keys,
        protein_type,
    )


def main():
    setup_logger()

    root = sys.argv[1]  # enthält NUR *.faa /
    threads = int(sys.argv[2])
    veryfasttree = sys.argv[3]
    globdb_faa = sys.argv[4]
    taxonomy_tsv = sys.argv[5]

    # Robustheit: Eingaben prüfen
    assert os.path.isdir(root), f"Input dir not found: {root}"
    assert os.path.exists(veryfasttree), f"VeryFastTree bin not found: {veryfasttree}"
    assert os.path.exists(globdb_faa), f"GloBDB FASTA not found: {globdb_faa}"
    assert os.path.exists(taxonomy_tsv), f"Taxonomy TSV not found: {taxonomy_tsv}"

    entries = sorted(os.listdir(root))
    fasta_ext = {".faa"}

    for name in entries:
        fp = os.path.join(root, name)
        if not os.path.isfile(fp):
            continue
        if os.path.splitext(name)[1].lower() not in fasta_ext:
            continue

        protein = os.path.splitext(name)[0]
        out_dir = os.path.join(root, protein)  # separater Ordner pro Protein
        process_protein_fasta(
            fp,
            out_dir,
            threads,
            veryfasttree,
            globdb_faa,
            taxonomy_tsv,
        )


if __name__ == "__main__":
    main()

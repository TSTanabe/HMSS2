#!/usr/bin/env python3
import sqlite3
import csv
import sys
from typing import Iterable, List, Tuple, Iterator, Optional


# --------------------------- Schema & Indizes ---------------------------


def create_schema(con: sqlite3.Connection) -> None:
    cur = con.cursor()
    # Performance-PRAGMAs NUR für den Import (riskanter bei Crash)
    cur.executescript("""
    PRAGMA page_size = 8192;        -- nur bei neuer DB wirksam
    PRAGMA journal_mode = OFF;      -- Max-Speed (nicht crash-sicher!)
    PRAGMA synchronous = OFF;       -- kein fsync -> schneller
    PRAGMA cache_size = -200000;    -- ~200 MB Cache
    PRAGMA temp_store = MEMORY;
    PRAGMA locking_mode = EXCLUSIVE;

    CREATE TABLE IF NOT EXISTS hmm_hits (
      id           INTEGER PRIMARY KEY,   -- Surrogat-PK (ROWID); auto-indiziert
      combined_id  TEXT    NOT NULL,      -- Spalte 0
      genome_id    TEXT    NOT NULL,      -- Prefix vor "___"
      protein_id   TEXT    NOT NULL,      -- Suffix nach "___"
      query        TEXT    NOT NULL,      -- Spalte 3
      score        INTEGER,               -- Spalte 7
      start        INTEGER,               -- Spalte 17
      end          INTEGER                -- Spalte 18
    );
    """)
    con.commit()


def ensure_secondary_indexes(con: sqlite3.Connection) -> None:
    cur = con.cursor()
    cur.executescript("""
    CREATE INDEX IF NOT EXISTS ix_hmm_combined  ON hmm_hits(combined_id);
    CREATE INDEX IF NOT EXISTS ix_hmm_genome    ON hmm_hits(genome_id);
    CREATE INDEX IF NOT EXISTS ix_hmm_protein   ON hmm_hits(protein_id);
    CREATE INDEX IF NOT EXISTS ix_hmm_query     ON hmm_hits(query);
    CREATE INDEX IF NOT EXISTS ix_hmm_genome_query  ON hmm_hits(genome_id, query);
    CREATE INDEX IF NOT EXISTS ix_hmm_protein_query ON hmm_hits(protein_id, query);
    """)
    con.commit()


# --------------------------- Parser & Loader ---------------------------


def _split_ids(combined_id: str) -> Tuple[str, str]:
    # "genome___protein" -> (genome, protein); robust degradiert
    if "___" in combined_id:
        return combined_id.split("___", 1)
    return "", combined_id


def _iter_rows(path: str) -> Iterable[Tuple[str, str, str, str, str, str, str]]:
    # Erwartet: TSV; nutzt Spalten 0,3,7,17,18
    with open(path, "r", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        for cols in rdr:
            if len(cols) <= 18:
                continue
            combined = cols[0]
            query = cols[3]
            c7 = cols[7]
            c17 = cols[17]
            c18 = cols[18]
            genome, protein = _split_ids(combined)
            yield (combined, genome, protein, query, c7, c17, c18)



def get_hits_by_genome_from_report_db(
    db_path: str,
    genome_id: str,
    chunk_size: int = 100_000,
) -> Iterator[Tuple[str, str, str, str, int, int, int]]:
    """
    Rein lesend:
      - öffnet DB mit mode=ro&immutable=1 (keine Locks/Journale/Schreibzugriffe).
      - benötigt einen Index auf (genome_id) für Speed.
      - streamt Ergebnisse in Batches via fetchmany(chunk_size).

    Returns:
        Iterator über (combined_id, genome_id, protein_id, query_suffix, bitscore, hsp_start, hsp_end)
    """
    # Read-only & immutable öffnen (keine Writes, keine Journal/Lock-Updates)
    uri = f"file:{db_path}?mode=ro&immutable=1"
    con = sqlite3.connect(uri, uri=True)
    try:
        cur = con.cursor()

        # Nur benötigte Spalten selektieren; CASTs sind hier optional,
        # falls die Eingabedatei numerische Felder als Text geliefert hat.
        cur.execute(
            """
            SELECT 
                combined_id,
                genome_id,
                protein_id,
                query,
                CAST(score AS INTEGER) AS bitscore,
                CAST(start AS INTEGER) AS hsp_start,
                CAST(end   AS INTEGER) AS hsp_end
            FROM hmm_hits
            WHERE genome_id = ?
            """,
            (genome_id,),
        )

        while True:
            rows = cur.fetchmany(chunk_size)
            if not rows:
                break
            for (
                combined_id,
                g_id,
                protein_id,
                query,
                bitscore,
                hsp_start,
                hsp_end,
            ) in rows:
                # Prefix am ersten/letzten "_" entfernen? (du nutzt bisher das letzte "_")
                query_suffix = query.rsplit("_", 1)[-1] if "_" in query else query
                yield (
                    combined_id,
                    g_id,
                    protein_id,
                    query_suffix,
                    bitscore,
                    hsp_start,
                    hsp_end,
                )
    finally:
        con.close()


def fetch_hits_for_pairs_ro(
    db_path: str,
    pairs: Iterable[Tuple[str, str]],
    attached_db_path: Optional[str] = None,
    *,
    max_pairs: int = 499,
) -> List[sqlite3.Row]:
    """
    Liest aus hmm_hits alle Zeilen zu gegebenen (protein_id, query)-Paaren.
    - Öffnet Datenbank(en) strikt read-only & immutable (kein Schreiben, kein Lock/Journal).
    - Verwendet gechunkte Row-Value IN (VALUES …) wegen 999-Placeholder-Limit.
    - Optional: zweite DB anhängen und beide Ergebnisse via UNION ALL zusammenführen.

    Voraussetzungen für Geschwindigkeit:
      CREATE INDEX IF NOT EXISTS idx_hits_protein_query ON hmm_hits (protein_id, query);

    Args:
        db_path: Pfad zur Haupt-DB.
        pairs: Iterable von (protein_id, query).
        attached_db_path: Optionaler Pfad zu einer zweiten DB (read-only angehängt).
        max_pairs: Max. Anzahl Paare pro Chunk (Default 499 => 998 Platzhalter).

    Returns:
        Liste von sqlite3.Row (Zeilen aus main.hmm_hits und ggf. db2.hmm_hits).
    """
    pair_list = list(pairs)
    if not pair_list:
        return []

    rows: List[sqlite3.Row] = []
    main_uri = f"file:{db_path}?mode=ro&immutable=1"

    with sqlite3.connect(main_uri, uri=True) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        if attached_db_path:
            attach_uri = f"file:{attached_db_path}?mode=ro&immutable=1"
            # URI sicher per Parameter binden
            cur.execute("ATTACH DATABASE ? AS db2", (attach_uri,))

        try:
            for i in range(0, len(pair_list), max_pairs):
                chunk = pair_list[i : i + max_pairs]
                placeholders = ",".join(["(?,?)"] * len(chunk))
                flat_params = [x for pair in chunk for x in pair]

                if attached_db_path:
                    sql = f"""
                        SELECT *
                        FROM hmm_hits
                        WHERE (protein_id, query) IN (VALUES {placeholders})
                        UNION ALL
                        SELECT *
                        FROM db2.hmm_hits
                        WHERE (protein_id, query) IN (VALUES {placeholders})
                    """
                    params = flat_params + flat_params
                else:
                    sql = f"""
                        SELECT *
                        FROM hmm_hits
                        WHERE (protein_id, query) IN (VALUES {placeholders})
                    """
                    params = flat_params

                rows.extend(cur.execute(sql, params).fetchall())
        finally:
            if attached_db_path:
                # Optionales sauberes Abhängen (beim Schließen der Connection auch ok)
                try:
                    cur.execute("DETACH DATABASE db2")
                except sqlite3.Error:
                    pass

    return rows

#!/usr/bin/env python3
import sqlite3
import sys


def show_proteins_without_sequence(db_path: str, limit: int = 20):
    """Zeigt die Zeilen aus der Tabelle Proteins ohne die Sequenzspalte."""
    with sqlite3.connect(db_path) as con:
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # Alle Spaltennamen außer 'sequence' holen
        cur.execute("PRAGMA table_info(Proteins);")
        cols = [row[1] for row in cur.fetchall() if row[1] != "sequence"]
        col_list = ", ".join(cols)

        # Daten abrufen (limit optional)
        cur.execute(f"SELECT {col_list} FROM Proteins LIMIT {limit};")
        rows = cur.fetchall()

        # Ausgabe
        print(f"Found {len(rows)} rows (showing up to {limit}):\n")
        for r in rows:
            print({k: r[k] for k in cols})


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python show_proteins.py <database.db> [limit]")
        sys.exit(1)

    db_file = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    show_proteins_without_sequence(db_file, limit)

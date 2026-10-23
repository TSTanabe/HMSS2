import sqlite3
import sys

db_path = sys.argv[1]

with sqlite3.connect(db_path) as conn:
    conn.row_factory = sqlite3.Row

    tables = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
    ).fetchall()

    for table in tables:
        name = table["name"]
        print(f"\n{'=' * 80}\nTABLE: {name}\n{'=' * 80}")

        for row in conn.execute(f'SELECT * FROM "{name}"'):
            print(dict(row))

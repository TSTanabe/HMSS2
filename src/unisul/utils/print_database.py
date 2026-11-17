
import sqlite3
import sys


def dump_database(database_path):
    conn = sqlite3.connect(database_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    tables = cur.execute(
        """
        SELECT name
        FROM sqlite_master
        WHERE type='table'
          AND name NOT LIKE 'sqlite_%'
        ORDER BY name
        """
    ).fetchall()

    for table_row in tables:
        table = table_row["name"]

        print("\n" + "=" * 80)
        print(f"TABLE: {table}")
        print("=" * 80)

        rows = cur.execute(f'SELECT * FROM "{table}"').fetchall()

        if not rows:
            print("<empty>")
            continue

        columns = rows[0].keys()
        print("\t".join(columns))

        for row in rows:
            print(
                "\t".join("" if row[col] is None else str(row[col]) for col in columns)
            )

    conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} database.sqlite")

    dump_database(sys.argv[1])

import sqlite3
from typing import List

"""
With a given list of protein domains get all keywords that includes all these proteins 
"""

def find_csbs_with_proteins_db(
    database: str,
    proteins: List[str],
    keyword_prefix: str = None,  # optional: only csb or dsb keywords
) -> List[str]:
    """
    Liefert die CSB-Keyword-Namen aus der Tabelle Keywords, für die es mindestens
    einen Cluster gibt, der *für jeden* gewünschten Proteintyp (proteins) ein
    passendes Domain-Vorkommen enthält.

    Matching-Regel für den Proteintyp:
      - exakte Domain-Gleichheit ODER
      - Domain endet auf '_<Proteintyp>' (z. B. 'grp3_SQRI' → 'SQRI').

    Nutzt EXISTS pro Typ -> gute Nutzung der vorhandenen Indizes.
    """
    if not proteins:
        return []

    sql = """\
    SELECT DISTINCT k.keyword
    FROM Keywords k
    WHERE (? = '' OR k.keyword LIKE ?)
    """
    params: List[str] = ["", ""]  # default: kein Prefix-Filter
    if keyword_prefix:
        params = ["x", f"{keyword_prefix}%"]  # aktiviere Prefix-Filter

    # Für jeden geforderten Proteintyp eine EXISTS-Klausel
    for _ in proteins:
        sql += """
        AND EXISTS (
          SELECT 1
          FROM Proteins p
          JOIN Domains d ON d.proteinID = p.proteinID
          WHERE p.clusterID = k.clusterID
            AND (d.domain = ?
                 OR d.domain LIKE '%' || '_' || ?)  -- suffix match
        )"""

    # Parameter anhängen: je Protein 2 Stück (exakt, suffix)
    for typ in proteins:
        params.extend([typ, typ])

    with (
        sqlite3.connect(f"file:{database}?mode=ro&immutable=1", uri=True) as con):
        con.execute("PRAGMA foreign_keys = ON;")
        con.execute("PRAGMA query_only = ON;")  # no writing
        con.execute("PRAGMA journal_mode = OFF;")  # no journal for read only
        con.execute("PRAGMA synchronous = OFF;")  # no sync
        cur = con.cursor()
        cur.execute(sql, params)
        return [row[0] for row in cur.fetchall()]
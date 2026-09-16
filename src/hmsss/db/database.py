#!/usr/bin/python
import os
import sqlite3
import sys
import time
import traceback
from typing import Any, Dict, Iterable, List, Set, Tuple

from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


########## Write output to Database Routines ##########


def create_database(database: str) -> None:
    """
    Create the HMSS3 SQLite database.

    Internal relations use INTEGER primary/foreign keys for efficient joins.

    Human-readable/external identifiers such as genomeID, proteinID,
    clusterID and domain remain UNIQUE and are used at the Python/API level.

    Large protein sequences are stored separately from Proteins.

    Presence tables GenomeDomains and ClusterDomains are used later for
    fast candidate selection during -fd and -fc database fetches.
    """

    logger.info("Creating database %s", database)

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON;")

        # ==========================================================
        # Genomes
        # ==========================================================

        cur.execute("""
            CREATE TABLE Genomes (
                genome_pk       INTEGER PRIMARY KEY,
                genomeID        TEXT NOT NULL UNIQUE,

                Superkingdom    TEXT DEFAULT NULL,
                Clade           TEXT DEFAULT NULL,
                Phylum          TEXT DEFAULT NULL,
                Class           TEXT DEFAULT NULL,
                Ordnung         TEXT DEFAULT NULL,
                Family          TEXT DEFAULT NULL,
                Genus           TEXT DEFAULT NULL,
                Species         TEXT DEFAULT NULL,
                Strain          TEXT DEFAULT NULL,

                TypeStrain      INTEGER DEFAULT NULL,
                Completeness    REAL DEFAULT NULL,
                Contamination   REAL DEFAULT NULL,
                dRep            INTEGER DEFAULT NULL,

                NCBITaxon       INTEGER DEFAULT NULL,
                NCBIProject     INTEGER DEFAULT NULL,
                NCBIBioproject  TEXT DEFAULT NULL,
                NCBIBiosample   TEXT DEFAULT NULL,
                NCBIAssembly    TEXT DEFAULT NULL
            );
        """)

        # ==========================================================
        # Clusters
        # ==========================================================

        cur.execute("""
            CREATE TABLE Clusters (
                cluster_pk  INTEGER PRIMARY KEY,
                clusterID   TEXT NOT NULL UNIQUE,
                genome_pk   INTEGER NOT NULL,

                FOREIGN KEY (genome_pk)
                    REFERENCES Genomes(genome_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE
            );
        """)

        # ==========================================================
        # Proteins
        # ==========================================================

        cur.execute("""
            CREATE TABLE Proteins (
                protein_pk       INTEGER PRIMARY KEY,
                proteinID        TEXT NOT NULL UNIQUE,

                genome_pk        INTEGER NOT NULL,
                cluster_pk       INTEGER DEFAULT NULL,

                locustag         TEXT DEFAULT NULL,
                contig           TEXT DEFAULT NULL,
                start            INTEGER DEFAULT NULL,
                end              INTEGER DEFAULT NULL,
                strand           TEXT DEFAULT NULL,

                comment          TEXT DEFAULT NULL,
                alternative_hit  TEXT DEFAULT NULL,

                dom_count        INTEGER DEFAULT NULL,
                valid_hit        INTEGER NOT NULL DEFAULT 0,

                FOREIGN KEY (genome_pk)
                    REFERENCES Genomes(genome_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                FOREIGN KEY (cluster_pk)
                    REFERENCES Clusters(cluster_pk)
                    ON DELETE SET NULL
                    ON UPDATE CASCADE,

                CHECK (valid_hit IN (0, 1))
            );
        """)

        # ==========================================================
        # Protein sequences
        #
        # Large payload deliberately separated from Proteins.
        # ==========================================================

        cur.execute("""
            CREATE TABLE ProteinSequences (
                protein_pk  INTEGER PRIMARY KEY,
                sequence    TEXT NOT NULL,

                FOREIGN KEY (protein_pk)
                    REFERENCES Proteins(protein_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE
            );
        """)

        # ==========================================================
        # Domain types
        #
        # Domain names are stored once and represented internally by
        # small INTEGER keys.
        # ==========================================================

        cur.execute("""
            CREATE TABLE DomainTypes (
                domain_pk   INTEGER PRIMARY KEY,
                domain      TEXT NOT NULL UNIQUE
            );
        """)

        # ==========================================================
        # Domains
        #
        # No artificial domain-hit ID required.
        # The natural hit identity is the composite primary key.
        # ==========================================================

        cur.execute("""
            CREATE TABLE Domains (
                protein_pk         INTEGER NOT NULL,
                domain_pk          INTEGER NOT NULL,

                domStart           INTEGER NOT NULL,
                domEnd             INTEGER NOT NULL,

                score              REAL DEFAULT NULL,
                blast_score_ratio  REAL DEFAULT NULL,
                identity           REAL DEFAULT NULL,

                PRIMARY KEY (
                    protein_pk,
                    domain_pk,
                    domStart,
                    domEnd
                ),

                FOREIGN KEY (protein_pk)
                    REFERENCES Proteins(protein_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                FOREIGN KEY (domain_pk)
                    REFERENCES DomainTypes(domain_pk)
                    ON DELETE RESTRICT
                    ON UPDATE CASCADE
            ) WITHOUT ROWID;
        """)

        # ==========================================================
        # Keywords
        # ==========================================================

        cur.execute("""
            CREATE TABLE Keywords (
                cluster_pk    INTEGER NOT NULL,
                keyword       TEXT NOT NULL,

                completeness  TEXT DEFAULT NULL,
                collinearity  TEXT DEFAULT NULL,

                PRIMARY KEY (
                    cluster_pk,
                    keyword
                ),

                FOREIGN KEY (cluster_pk)
                    REFERENCES Clusters(cluster_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE
            ) WITHOUT ROWID;
        """)

        # ==========================================================
        # Genome/domain presence table
        #
        # Used for fast -fd candidate selection.
        # ==========================================================

        cur.execute("""
            CREATE TABLE GenomeDomains (
                genome_pk      INTEGER NOT NULL,
                domain_pk      INTEGER NOT NULL,
                valid_present  INTEGER NOT NULL DEFAULT 0,

                PRIMARY KEY (
                    genome_pk,
                    domain_pk
                ),

                FOREIGN KEY (genome_pk)
                    REFERENCES Genomes(genome_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                FOREIGN KEY (domain_pk)
                    REFERENCES DomainTypes(domain_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                CHECK (valid_present IN (0, 1))
            ) WITHOUT ROWID;
        """)

        # ==========================================================
        # Cluster/domain presence table
        #
        # Used for fast -fc candidate selection.
        # ==========================================================

        cur.execute("""
            CREATE TABLE ClusterDomains (
                cluster_pk     INTEGER NOT NULL,
                domain_pk      INTEGER NOT NULL,
                valid_present  INTEGER NOT NULL DEFAULT 0,

                PRIMARY KEY (
                    cluster_pk,
                    domain_pk
                ),

                FOREIGN KEY (cluster_pk)
                    REFERENCES Clusters(cluster_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                FOREIGN KEY (domain_pk)
                    REFERENCES DomainTypes(domain_pk)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE,

                CHECK (valid_present IN (0, 1))
            ) WITHOUT ROWID;
        """)

        # ==========================================================
        # Metagenomes
        #
        # Kept largely compatible for now. The read/placement part
        # can be converted to integer surrogate keys separately later.
        # ==========================================================

        cur.execute("""
            CREATE TABLE Metagenomes (
                metagenomeID         TEXT PRIMARY KEY NOT NULL,
                genomeID             TEXT NOT NULL,

                forward_reads        INTEGER DEFAULT NULL,
                reverse_reads        INTEGER DEFAULT NULL,
                prokaryotic_fraction REAL DEFAULT 1.0,

                FOREIGN KEY (genomeID)
                    REFERENCES Genomes(genomeID)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE
            );
        """)

        # ==========================================================
        # GraftM package lengths
        # ==========================================================

        cur.execute("""
            CREATE TABLE GpkgLengths (
                domain_type     TEXT PRIMARY KEY NOT NULL,
                protein_length  INTEGER NOT NULL
            );
        """)

        # ==========================================================
        # Read lineage
        # ==========================================================

        cur.execute("""
            CREATE TABLE Lineage (
                lineageID   TEXT PRIMARY KEY NOT NULL,

                root        TEXT DEFAULT NULL,
                kingdom     TEXT DEFAULT NULL,
                phylum      TEXT DEFAULT NULL,
                class       TEXT DEFAULT NULL,
                "order"     TEXT DEFAULT NULL,
                family      TEXT DEFAULT NULL,
                genus       TEXT DEFAULT NULL,
                species     TEXT DEFAULT NULL,

                raw_lineage TEXT DEFAULT NULL
            );
        """)

        # ==========================================================
        # Read placements
        #
        # Kept compatible for the first migration phase.
        # ==========================================================

        cur.execute("""
            CREATE TABLE Placement (
                domain_type   TEXT NOT NULL,
                readID        TEXT NOT NULL,
                metagenomeID  TEXT NOT NULL,

                proteinID     TEXT DEFAULT NULL,
                lineageID     TEXT DEFAULT NULL,

                dom_start     INTEGER DEFAULT NULL,
                dom_end       INTEGER DEFAULT NULL,
                coverage      REAL DEFAULT NULL,

                sequence      TEXT DEFAULT NULL,
                alignment     TEXT DEFAULT NULL,

                PRIMARY KEY (
                    metagenomeID,
                    domain_type,
                    readID
                ),

                FOREIGN KEY (proteinID)
                    REFERENCES Proteins(proteinID)
                    ON DELETE SET NULL
                    ON UPDATE CASCADE,

                FOREIGN KEY (lineageID)
                    REFERENCES Lineage(lineageID)
                    ON DELETE SET NULL
                    ON UPDATE CASCADE,

                FOREIGN KEY (metagenomeID)
                    REFERENCES Metagenomes(metagenomeID)
                    ON DELETE CASCADE
                    ON UPDATE CASCADE
            );
        """)

    logger.info("Created database")


def index_database(database: str) -> None:
    """Create indexes for the normalized HMSS3 database and update planner statistics."""

    indexes = [
        (
            "idx_clusters_genome",
            "CREATE INDEX IF NOT EXISTS idx_clusters_genome ON Clusters(genome_pk)",
        ),
        (
            "idx_proteins_genome",
            "CREATE INDEX IF NOT EXISTS idx_proteins_genome ON Proteins(genome_pk)",
        ),
        (
            "idx_proteins_cluster",
            "CREATE INDEX IF NOT EXISTS idx_proteins_cluster ON Proteins(cluster_pk)",
        ),
        (
            "idx_proteins_genome_order",
            "CREATE INDEX IF NOT EXISTS idx_proteins_genome_order ON Proteins(genome_pk, contig, start, protein_pk)",
        ),
        (
            "idx_proteins_cluster_order",
            "CREATE INDEX IF NOT EXISTS idx_proteins_cluster_order ON Proteins(cluster_pk, genome_pk, contig, start, protein_pk)",
        ),
        (
            "idx_proteins_cluster_valid",
            "CREATE INDEX IF NOT EXISTS idx_proteins_cluster_valid ON Proteins(cluster_pk, valid_hit, protein_pk)",
        ),
        # Domains PRIMARY KEY already starts with protein_pk.
        # This index provides the reverse access direction: domain -> proteins.
        (
            "idx_domains_domain_protein",
            "CREATE INDEX IF NOT EXISTS idx_domains_domain_protein ON Domains(domain_pk, protein_pk)",
        ),
        (
            "idx_domains_protein_start",
            "CREATE INDEX IF NOT EXISTS idx_domains_protein_start ON Domains(protein_pk, domStart, domain_pk)",
        ),
        # Keywords PRIMARY KEY already starts with cluster_pk.
        (
            "idx_keywords_keyword_cluster",
            "CREATE INDEX IF NOT EXISTS idx_keywords_keyword_cluster ON Keywords(keyword, cluster_pk)",
        ),
        # Reverse indexes for presence searches.
        (
            "idx_genomedomains_domain",
            "CREATE INDEX IF NOT EXISTS idx_genomedomains_domain ON GenomeDomains(domain_pk, valid_present, genome_pk)",
        ),
        (
            "idx_clusterdomains_domain",
            "CREATE INDEX IF NOT EXISTS idx_clusterdomains_domain ON ClusterDomains(domain_pk, valid_present, cluster_pk)",
        ),
        # Read/metagenome side remains largely unchanged.
        (
            "idx_metagenomes_genome",
            "CREATE INDEX IF NOT EXISTS idx_metagenomes_genome ON Metagenomes(genomeID)",
        ),
        (
            "idx_placement_protein",
            "CREATE INDEX IF NOT EXISTS idx_placement_protein ON Placement(proteinID)",
        ),
        (
            "idx_placement_lineage",
            "CREATE INDEX IF NOT EXISTS idx_placement_lineage ON Placement(lineageID)",
        ),
        (
            "idx_placement_domain_meta_read",
            "CREATE INDEX IF NOT EXISTS idx_placement_domain_meta_read ON Placement(domain_type, metagenomeID, readID)",
        ),
    ]

    logger.info("Indexing database %s", database)

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")

        for name, sql in indexes:
            logger.debug("Creating index %s", name)
            cur.execute(sql)

        logger.info("Updating SQLite query planner statistics")
        cur.execute("ANALYZE")
        cur.execute("PRAGMA optimize")

    logger.info("Finished indexing database %s", database)


def rebuild_presence_tables(database: str) -> None:
    """
    Rebuild GenomeDomains and ClusterDomains from the current Proteins/Domains state.

    valid_present = 1 if at least one corresponding protein hit has valid_hit = 1.
    """

    logger.info("Building domain presence tables")

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute("PRAGMA temp_store = MEMORY")

        cur.execute("DELETE FROM GenomeDomains")
        cur.execute("""
            INSERT INTO GenomeDomains(genome_pk, domain_pk, valid_present)
            SELECT p.genome_pk, d.domain_pk, MAX(CASE WHEN p.valid_hit = 1 THEN 1 ELSE 0 END)
            FROM Proteins p
            JOIN Domains d ON d.protein_pk = p.protein_pk
            GROUP BY p.genome_pk, d.domain_pk
        """)

        genome_rows = cur.rowcount

        cur.execute("DELETE FROM ClusterDomains")
        cur.execute("""
            INSERT INTO ClusterDomains(cluster_pk, domain_pk, valid_present)
            SELECT p.cluster_pk, d.domain_pk, MAX(CASE WHEN p.valid_hit = 1 THEN 1 ELSE 0 END)
            FROM Proteins p
            JOIN Domains d ON d.protein_pk = p.protein_pk
            WHERE p.cluster_pk IS NOT NULL
            GROUP BY p.cluster_pk, d.domain_pk
        """)

        cluster_rows = cur.rowcount

    logger.info(
        "Presence tables rebuilt: %d genome-domain and %d cluster-domain entries",
        genome_rows,
        cluster_rows,
    )


def _lookup_pk_map(
    cur: sqlite3.Cursor,
    values: Iterable[str],
    *,
    temp_table: str,
    source_table: str,
    external_col: str,
    pk_col: str,
    strict: bool = True,
) -> Dict[str, int]:
    """
    Resolve external TEXT identifiers to internal INTEGER primary keys.

    Large identifier sets are passed through a TEMP table instead of IN (...).

    The table/column names passed here are internal constants only and must
    never originate from user input.
    """

    ids = sorted(
        {
            str(value).strip()
            for value in values
            if value is not None and str(value).strip()
        }
    )

    if not ids:
        return {}

    cur.execute(
        f"""
        CREATE TEMP TABLE IF NOT EXISTS {temp_table} (
            external_id TEXT PRIMARY KEY
        ) WITHOUT ROWID
        """
    )

    cur.execute(f"DELETE FROM {temp_table}")

    cur.executemany(
        f"""
        INSERT OR IGNORE INTO {temp_table}(external_id)
        VALUES (?)
        """,
        ((value,) for value in ids),
    )

    cur.execute(
        f"""
        SELECT
            s.{external_col} AS external_id,
            s.{pk_col} AS internal_pk
        FROM {source_table} s
        JOIN {temp_table} t
          ON t.external_id = s.{external_col}
        """
    )

    result = {str(row[0]): int(row[1]) for row in cur}

    if strict and len(result) != len(ids):
        missing = [value for value in ids if value not in result]

        raise KeyError(
            f"Could not resolve {len(missing)} identifiers in "
            f"{source_table}.{external_col}. "
            f"Examples: {missing[:10]}"
        )

    return result


def _ensure_domain_types(
    cur: sqlite3.Cursor,
    domains: Iterable[str],
) -> Dict[str, int]:
    """
    Ensure all requested domain names exist in DomainTypes and return

        domain_name -> domain_pk
    """

    domain_names = sorted(
        {
            str(domain).strip()
            for domain in domains
            if domain is not None and str(domain).strip()
        }
    )

    if not domain_names:
        return {}

    cur.executemany(
        """
        INSERT OR IGNORE INTO DomainTypes(domain)
        VALUES (?)
        """,
        ((domain,) for domain in domain_names),
    )

    return _lookup_pk_map(
        cur,
        domain_names,
        temp_table="tmp_lookup_domain_types",
        source_table="DomainTypes",
        external_col="domain",
        pk_col="domain_pk",
    )


def insert_database_genome_ids(database: str, genome_ids: Set[str]) -> None:
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA journal_mode = OFF;")

        cur.executemany(
            """
            INSERT OR IGNORE INTO Genomes(genomeID)
            VALUES (?)
            """,
            ((genomeID,) for genomeID in genome_ids),
        )


def insert_database_proteins(
    database: str,
    protein_dict: Dict[str, Any],
) -> None:
    """
    Insert proteins, sequences and domains using internal INTEGER keys.

    External identifiers remain unchanged:
        genomeID
        proteinID
        domain

    Internally stored relations use:
        genome_pk
        protein_pk
        domain_pk
    """

    if not protein_dict:
        return

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()

            cur.execute("PRAGMA foreign_keys = ON;")
            cur.execute("PRAGMA synchronous = OFF;")
            cur.execute("PRAGMA journal_mode = OFF;")
            cur.execute("PRAGMA temp_store = MEMORY;")

            protein_list = sorted(
                protein_dict.values(),
                key=lambda x: (
                    x.genomeID or "",
                    x.gene_contig or "",
                    x.gene_start or 0,
                ),
            )

            # ------------------------------------------------------
            # 1. Collect external identifiers
            # ------------------------------------------------------

            genome_ids: Set[str] = set()
            domain_names: Set[str] = set()

            for protein in protein_list:
                if protein.genomeID:
                    genome_ids.add(str(protein.genomeID))

                for domain in protein.domains:
                    if domain.domain:
                        domain_names.add(str(domain.domain))

            # ------------------------------------------------------
            # 2. Resolve genomeID -> genome_pk
            # ------------------------------------------------------

            genome_pk_map = _lookup_pk_map(
                cur,
                genome_ids,
                temp_table="tmp_lookup_protein_genomes",
                source_table="Genomes",
                external_col="genomeID",
                pk_col="genome_pk",
            )

            # ------------------------------------------------------
            # 3. Ensure domain types and obtain domain -> domain_pk
            # ------------------------------------------------------

            domain_pk_map = _ensure_domain_types(
                cur,
                domain_names,
            )

            # ------------------------------------------------------
            # 4. Build protein records
            # ------------------------------------------------------

            protein_records = []

            # Keep sequence/domain records temporarily using the
            # external proteinID. protein_pk does not exist until after
            # Proteins has been inserted.
            sequence_records_pending = []
            domain_records_pending = []

            protein_ids: Set[str] = set()

            for protein in protein_list:
                genome_id = str(protein.genomeID or "").strip()

                if not genome_id:
                    logger.warning(
                        "Skipping protein without genomeID: %s",
                        protein.proteinID,
                    )
                    continue

                genome_pk = genome_pk_map[genome_id]

                # Keep the existing HMSS3 protein identifier convention.
                protein_id = f"{genome_id}-{protein.proteinID}"

                protein_ids.add(protein_id)

                domains = protein.domains

                protein_records.append(
                    (
                        protein_id,
                        genome_pk,
                        protein.gene_locustag,
                        protein.gene_contig,
                        protein.gene_start,
                        protein.gene_end,
                        protein.gene_strand,
                        protein.get_selection_comment_csv(),
                        protein.alternative_hit,
                        len(domains),
                        int(bool(protein.valid_hit)),
                    )
                )

                sequence = protein.get_sequence()

                if sequence is not None:
                    sequence_records_pending.append(
                        (
                            protein_id,
                            sequence,
                        )
                    )

                for domain in domains:
                    domain_name = str(domain.domain)

                    domain_records_pending.append(
                        (
                            protein_id,
                            domain_pk_map[domain_name],
                            domain.start,
                            domain.end,
                            domain.score,
                            domain.identity,
                            domain.bsr,
                        )
                    )

            # ------------------------------------------------------
            # 5. Insert Proteins
            # ------------------------------------------------------

            cur.executemany(
                """
                INSERT INTO Proteins (
                    proteinID,
                    genome_pk,
                    locustag,
                    contig,
                    start,
                    end,
                    strand,
                    comment,
                    alternative_hit,
                    dom_count,
                    valid_hit
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)

                ON CONFLICT(proteinID) DO UPDATE SET

                    locustag = COALESCE(
                        Proteins.locustag,
                        excluded.locustag
                    ),

                    contig = COALESCE(
                        Proteins.contig,
                        excluded.contig
                    ),

                    start = COALESCE(
                        Proteins.start,
                        excluded.start
                    ),

                    end = COALESCE(
                        Proteins.end,
                        excluded.end
                    ),

                    strand = COALESCE(
                        Proteins.strand,
                        excluded.strand
                    ),

                    comment = COALESCE(
                        Proteins.comment,
                        excluded.comment
                    ),

                    alternative_hit = COALESCE(
                        Proteins.alternative_hit,
                        excluded.alternative_hit
                    ),

                    dom_count = COALESCE(
                        COALESCE(Proteins.dom_count, 0),
                        COALESCE(excluded.dom_count, 0)
                    ),

                    valid_hit = COALESCE(
                        Proteins.valid_hit,
                        excluded.valid_hit
                    )
                """,
                protein_records,
            )

            # ------------------------------------------------------
            # 6. Resolve proteinID -> protein_pk
            # ------------------------------------------------------

            protein_pk_map = _lookup_pk_map(
                cur,
                protein_ids,
                temp_table="tmp_lookup_inserted_proteins",
                source_table="Proteins",
                external_col="proteinID",
                pk_col="protein_pk",
            )

            # ------------------------------------------------------
            # 7. Insert sequences
            # ------------------------------------------------------

            sequence_records = [
                (
                    protein_pk_map[protein_id],
                    sequence,
                )
                for protein_id, sequence in sequence_records_pending
            ]

            if sequence_records:
                cur.executemany(
                    """
                    INSERT INTO ProteinSequences (
                        protein_pk,
                        sequence
                    )
                    VALUES (?, ?)

                    ON CONFLICT(protein_pk) DO UPDATE SET
                        sequence = excluded.sequence
                    """,
                    sequence_records,
                )

            # ------------------------------------------------------
            # 8. Insert domain hits
            # ------------------------------------------------------

            domain_records = [
                (
                    protein_pk_map[protein_id],
                    domain_pk,
                    dom_start,
                    dom_end,
                    score,
                    identity,
                    bsr,
                )
                for (
                    protein_id,
                    domain_pk,
                    dom_start,
                    dom_end,
                    score,
                    identity,
                    bsr,
                ) in domain_records_pending
            ]

            if domain_records:
                cur.executemany(
                    """
                    INSERT INTO Domains (
                        protein_pk,
                        domain_pk,
                        domStart,
                        domEnd,
                        score,
                        identity,
                        blast_score_ratio
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?)

                    ON CONFLICT(
                        protein_pk,
                        domain_pk,
                        domStart,
                        domEnd
                    ) DO NOTHING
                    """,
                    domain_records,
                )

            logger.info(
                "Inserted/updated %d proteins, %d domain hits and %d sequences.",
                len(protein_records),
                len(domain_records),
                len(sequence_records),
            )

    except Exception as e:
        logger.warning(
            f"Proteins were not inserted for {genome_id}.\nError: - {str(e)}\nTraceback: {traceback.format_exc()}"
        )

    return


def insert_database_clusters(
    database: str,
    cluster_dict: Dict[str, Any],
) -> None:
    """
    Insert clusters and assign proteins/keywords using internal cluster_pk.
    """

    if not cluster_dict:
        return

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()

            cur.execute("PRAGMA foreign_keys = ON;")
            cur.execute("PRAGMA synchronous = OFF;")
            cur.execute("PRAGMA journal_mode = OFF;")
            cur.execute("PRAGMA temp_store = MEMORY;")

            cluster_records_pending = []
            protein_updates_pending = []
            keyword_records_pending = []

            genome_ids: Set[str] = set()
            cluster_ids: Set[str] = set()
            protein_ids: Set[str] = set()

            # ------------------------------------------------------
            # Collect identifiers
            # ------------------------------------------------------

            for cluster in cluster_dict.values():
                cluster_id = cluster.get_cluster_id()

                if not cluster_id:
                    continue

                genome_id = str(cluster.genomeID)

                cluster_ids.add(cluster_id)
                genome_ids.add(genome_id)

                cluster_records_pending.append(
                    (
                        cluster_id,
                        genome_id,
                    )
                )

                for protein_id in cluster.get_genes():
                    full_protein_id = f"{genome_id}-{protein_id}"
                    protein_ids.add(full_protein_id)
                    protein_updates_pending.append(
                        (
                            cluster_id,
                            full_protein_id,
                        )
                    )

                for keyword in cluster.get_keywords():
                    keyword_records_pending.append(
                        (
                            cluster_id,
                            keyword.get_keyword(),
                            keyword.get_completeness(),
                            keyword.get_csb(),
                        )
                    )

            # ------------------------------------------------------
            # Resolve genome IDs
            # ------------------------------------------------------

            genome_pk_map = _lookup_pk_map(
                cur,
                genome_ids,
                temp_table="tmp_lookup_cluster_genomes",
                source_table="Genomes",
                external_col="genomeID",
                pk_col="genome_pk",
            )

            # ------------------------------------------------------
            # Insert clusters
            # ------------------------------------------------------

            cluster_records = [
                (
                    cluster_id,
                    genome_pk_map[genome_id],
                )
                for cluster_id, genome_id in cluster_records_pending
            ]

            cur.executemany(
                """
                INSERT INTO Clusters (
                    clusterID,
                    genome_pk
                )
                VALUES (?, ?)

                ON CONFLICT(clusterID) DO NOTHING
                """,
                cluster_records,
            )

            # ------------------------------------------------------
            # Resolve cluster_pk
            # ------------------------------------------------------

            cluster_pk_map = _lookup_pk_map(
                cur,
                cluster_ids,
                temp_table="tmp_lookup_inserted_clusters",
                source_table="Clusters",
                external_col="clusterID",
                pk_col="cluster_pk",
            )

            # ------------------------------------------------------
            # Resolve protein_pk
            # ------------------------------------------------------

            protein_pk_map = _lookup_pk_map(
                cur,
                protein_ids,
                temp_table="tmp_lookup_cluster_proteins",
                source_table="Proteins",
                external_col="proteinID",
                pk_col="protein_pk",
                strict=False,
            )

            missing_proteins = protein_ids - set(protein_pk_map)

            if missing_proteins:
                logger.warning(
                    "%d proteins referenced by clusters were not present "
                    "in Proteins. Examples: %s",
                    len(missing_proteins),
                    sorted(missing_proteins)[:10],
                )

            # ------------------------------------------------------
            # Assign cluster_pk to proteins
            # ------------------------------------------------------

            protein_updates = [
                (
                    cluster_pk_map[cluster_id],
                    protein_pk_map[protein_id],
                )
                for cluster_id, protein_id in protein_updates_pending
                if protein_id in protein_pk_map
            ]

            if protein_updates:
                cur.executemany(
                    """
                    UPDATE Proteins
                    SET cluster_pk = ?
                    WHERE protein_pk = ?
                    """,
                    protein_updates,
                )

            # ------------------------------------------------------
            # Insert keywords
            # ------------------------------------------------------

            keyword_records = [
                (
                    cluster_pk_map[cluster_id],
                    keyword,
                    completeness,
                    collinearity,
                )
                for (
                    cluster_id,
                    keyword,
                    completeness,
                    collinearity,
                ) in keyword_records_pending
            ]

            if keyword_records:
                cur.executemany(
                    """
                    INSERT INTO Keywords (
                        cluster_pk,
                        keyword,
                        completeness,
                        collinearity
                    )
                    VALUES (?, ?, ?, ?)

                    ON CONFLICT(cluster_pk, keyword) DO UPDATE SET
                        completeness = excluded.completeness,
                        collinearity = excluded.collinearity
                    """,
                    keyword_records,
                )

            logger.info(
                "Inserted/updated %d clusters, assigned %d proteins "
                "and stored %d keywords.",
                len(cluster_records),
                len(protein_updates),
                len(keyword_records),
            )

    except Exception as e:
        logger.warning(
            f"Due to an error - {str(e)}\nTraceback: {traceback.format_exc()}"
        )

    return


def insert_taxonomy_data(database: str, taxonomy_file: str) -> None:
    """
    Insert taxonomy data from a file into the Genomes table of the database.

    Args:
        database (str): Path to SQLite database.
        taxonomy_file (str): Path to the parsed taxonomy file (tab-separated).
    """
    try:
        # Check if the database file path is valid
        if not os.path.exists(os.path.dirname(database)):
            logger.error(
                f"Directory for database does not exist: {os.path.dirname(database)}"
            )
            return

        # Connect to the SQLite database
        with sqlite3.connect(database) as con:
            cur = con.cursor()

            # Read the taxonomy file and insert data into the Genomes table
            with open(taxonomy_file, "r") as file:
                # Skip the header
                next(file)

                # Read each line in the file
                for line in file:
                    fields = line.strip().split("\t")

                    if len(fields) < 8:
                        fields = parse_taxonomy_line(line, "NA")
                    if len(fields) < 8:
                        logger.warning(f"Line has insufficient columns: {line}")
                        continue

                    # Prepare the data for insertion
                    genome_id = fields[0]
                    superkingdom = fields[1]
                    clade = ""
                    phylum = fields[2]
                    class_ = fields[3]
                    order = fields[4]
                    family = fields[5]
                    genus = fields[6]
                    species = fields[7]

                    # Insert the data into the Genomes table
                    cur.execute(
                        """
                        INSERT INTO Genomes (genomeID, Superkingdom, Clade, Phylum, Class, Ordnung, Family, Genus, Species)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(genomeID) DO UPDATE SET
                        Superkingdom=excluded.Superkingdom,
                        Clade=excluded.Clade,
                        Phylum=excluded.Phylum,
                        Class=excluded.Class,
                        Ordnung=excluded.Ordnung,
                        Family=excluded.Family,
                        Genus=excluded.Genus,
                        Species=excluded.Species
                    """,
                        (
                            genome_id,
                            superkingdom,
                            clade,
                            phylum,
                            class_,
                            order,
                            family,
                            genus,
                            species,
                        ),
                    )

            # Commit the transaction
            con.commit()

        logger.info("Inserted taxonomy data from %s", taxonomy_file)
    except sqlite3.OperationalError as e:
        logger.error(f"SQLite Operational Error: {e}\nDatabase path: {database}")
    except sqlite3.Error as e:
        logger.error(f"SQLite Error: {e}")
    except FileNotFoundError as e:
        logger.error(f"File Not Found Error: {e}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")


def parse_taxonomy_line(line: str, na: str = "") -> List[str]:
    rank_keys = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    parts = line.rstrip("\n").split("\t")
    if len(parts) < 2:
        logger.error(
            f"Line must contain at least two tab-separated fields: <ID> and <taxonomy>\nLine: {parts}"
        )

    genome_id = parts[0].strip()
    tax_str = parts[1].strip()

    raw_tokens = [t.strip() for t in tax_str.split(";") if t.strip()]
    ranks = {k: na for k in rank_keys}

    prefix_to_rank = {
        "d__": "domain",
        "k__": "kingdom",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }

    for token in raw_tokens:
        for pre in prefix_to_rank:
            if token.startswith(pre):
                rank = prefix_to_rank[pre]
                value = token[len(pre) :].strip()
                # Leerzeichen in Unterstrich nur bei species
                if rank == "species":
                    value = value.replace(" ", "_")
                ranks[rank] = value
                break

    return [genome_id] + [ranks[k] for k in rank_keys]


##############################################################
######## Metagenome information to database routines #########
##############################################################


def insert_database_metagenomes(
    database: str,
    metagenome_dict: Dict[str, Tuple[str, int | None, int | None]],
) -> None:
    """
    Insert metagenomes into the Metagenomes table.

    Inputs
    ------
    database : str
        Path to sqlite database.
    metagenome_dict : dict
        {metagenomeID: (genomeID, forward_reads, reverse_reads)}
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""PRAGMA synchronous = OFF;""")
        cur.execute("""PRAGMA journal_mode = OFF;""")

        records = [
            (metagenomeID, genomeID, fwd, rev)
            for metagenomeID, (genomeID, fwd, rev) in metagenome_dict.items()
        ]

        cur.executemany(
            """
            INSERT OR IGNORE INTO Metagenomes
              (metagenomeID, genomeID, forward_reads, reverse_reads)
            VALUES (?, ?, ?, ?)
            """,
            records,
        )
        con.commit()
    con.close()
    return


def insert_database_lineages(database: str, reads: Dict[tuple, "Read"]) -> None:
    """
    Insert lineage information for all Read objects into the Lineage table.

    - Uses lineageID as PRIMARY KEY
    - Existing lineageIDs are ignored (INSERT OR IGNORE)

    Parameters
    ----------
    database : str
        Path to sqlite database.
    reads : dict
        Dict[(...), Read] or Dict[key, Read]; values must be Read objects
        with attributes:
          - lineageID : str
          - lineage   : dict[str, str]
    """
    records = []

    for r in reads.values():
        if not r.lineageID or not r.lineage:
            continue

        records.append(
            (
                r.lineageID,
                r.lineage.get("root", "NA"),
                r.lineage.get("k", "NA"),
                r.lineage.get("p", "NA"),
                r.lineage.get("c", "NA"),
                r.lineage.get("o", "NA"),
                r.lineage.get("f", "NA"),
                r.lineage.get("g", "NA"),
                r.lineage.get("s", "NA"),
            )
        )

    if not records:
        return

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA journal_mode = OFF;")

        cur.executemany(
            """
            INSERT OR IGNORE INTO Lineage
              (lineageID, root, kingdom, phylum, class, "order", family, genus, species)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            records,
        )
        con.commit()
    con.close()
    return


def insert_database_stub_proteins_from_reads(
    database: str,
    reads: Dict[tuple, "Read"],
) -> None:
    pending = {}
    genome_ids = set()

    for read in reads.values():
        protein_id = read.readID
        genome_id = read.genomeID

        if not protein_id or not genome_id:
            continue

        pending[protein_id] = genome_id
        genome_ids.add(genome_id)

    if not pending:
        return

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA journal_mode = OFF;")
        cur.execute("PRAGMA temp_store = MEMORY;")

        genome_pk_map = _lookup_pk_map(
            cur,
            genome_ids,
            temp_table="tmp_lookup_stub_genomes",
            source_table="Genomes",
            external_col="genomeID",
            pk_col="genome_pk",
        )

        records = [
            (
                protein_id,
                genome_pk_map[genome_id],
            )
            for protein_id, genome_id in pending.items()
        ]

        cur.executemany(
            """
            INSERT OR IGNORE INTO Proteins (
                proteinID,
                genome_pk
            )
            VALUES (?, ?)
            """,
            records,
        )


def insert_database_placements(database: str, reads: Dict[tuple, "Read"]) -> None:
    """
    Insert Read objects into Placement.

    Expects each Read to have at least:
      - type (domain_type)
      - readID
      - metagenomeID
      - lineageID (optional)
      - start, end, coverage (optional but recommended)
      - sequence, alignment (optional)

    Uses INSERT OR IGNORE to avoid duplicates for PRIMARY KEY (domain_type, readID).
    """
    records = []

    for r in reads.values():
        if not r.gpkg_name or not r.readID or not r.metagenomeID:
            continue

        records.append(
            (
                r.gpkg_name,  # domain_type
                r.readID,  # readID
                r.metagenomeID,  # metagenomeID
                r.readID,  # proteinID (stub mapping: proteinID == readID)
                r.lineageID or None,
                r.start if r.start is not None else None,
                r.end if r.end is not None else None,
                r.coverage if r.coverage is not None else None,
                r.sequence or None,
                r.alignment or None,
            )
        )

    if not records:
        return

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA journal_mode = OFF;")

        cur.executemany(
            """
            INSERT OR IGNORE INTO Placement
              (domain_type, readID, metagenomeID, proteinID, lineageID, dom_start, dom_end, coverage, sequence, alignment)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            records,
        )
        con.commit()
    con.close()
    return


def insert_database_gpkg_lengths(
    database: str,
    gpkg_length_dict: Dict[str, int],
) -> None:
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA journal_mode = OFF;")

        records = [(k, int(v)) for k, v in gpkg_length_dict.items()]

        cur.executemany(
            """
            INSERT INTO GpkgLengths (domain_type, protein_length)
            VALUES (?, ?)
            ON CONFLICT(domain_type) DO UPDATE SET
                protein_length = excluded.protein_length
            """,
            records,
        )
        con.commit()


##############################################################
########## Alter information from database routines ##########
##############################################################


def update_domain(
    database: str, protein_diction: Dict[str, Any], old_tag: str, new_tag: str
) -> None:
    """Change one domain type to another for the selected proteins only."""

    if not protein_diction or old_tag == new_tag:
        return

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("PRAGMA foreign_keys = ON")
            cur.execute("PRAGMA temp_store = MEMORY")

            old_row = cur.execute(
                "SELECT domain_pk FROM DomainTypes WHERE domain = ?", (old_tag,)
            ).fetchone()

            if old_row is None:
                logger.warning("Domain %s does not exist in DomainTypes", old_tag)
                return

            old_domain_pk = int(old_row[0])
            new_domain_pk = _ensure_domain_types(cur, [new_tag])[new_tag]

            protein_pk_map = _lookup_pk_map(
                cur,
                protein_diction.keys(),
                temp_table="tmp_lookup_update_domain_proteins",
                source_table="Proteins",
                external_col="proteinID",
                pk_col="protein_pk",
                strict=False,
            )

            if not protein_pk_map:
                return

            cur.execute("""
                CREATE TEMP TABLE IF NOT EXISTS tmp_update_domain_proteins (
                    protein_pk INTEGER PRIMARY KEY
                )
            """)
            cur.execute("DELETE FROM tmp_update_domain_proteins")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_update_domain_proteins(protein_pk) VALUES (?)",
                ((pk,) for pk in protein_pk_map.values()),
            )

            # Copy old hits with the new domain_pk. INSERT OR IGNORE also handles
            # the case where the corresponding new-tag hit already exists.
            cur.execute(
                """
                INSERT OR IGNORE INTO Domains
                    (protein_pk, domain_pk, domStart, domEnd, score, blast_score_ratio, identity)
                SELECT d.protein_pk, ?, d.domStart, d.domEnd, d.score, d.blast_score_ratio, d.identity
                FROM Domains d
                JOIN tmp_update_domain_proteins t ON t.protein_pk = d.protein_pk
                WHERE d.domain_pk = ?
            """,
                (new_domain_pk, old_domain_pk),
            )

            # Remove the corresponding old-tag hits.
            cur.execute(
                """
                DELETE FROM Domains
                WHERE domain_pk = ?
                  AND EXISTS (
                      SELECT 1 FROM tmp_update_domain_proteins t
                      WHERE t.protein_pk = Domains.protein_pk
                  )
            """,
                (old_domain_pk,),
            )

            changed = cur.rowcount

        logger.info("Updated %d domain hits from %s to %s", changed, old_tag, new_tag)

    except Exception:
        logger.exception("Error updating domains from %s to %s", old_tag, new_tag)
        raise


##############################################################
########## Fetch information from database routines ##########
##############################################################


def update_keywords(
    database: str, keyword_dict: Dict[str, Set[str]], batch_size: int = 10000
) -> None:
    """Add keywords to clusters while resolving external clusterID to cluster_pk."""

    if not keyword_dict:
        return

    cluster_ids = {cid for cluster_ids in keyword_dict.values() for cid in cluster_ids}

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("PRAGMA foreign_keys = ON")
            cur.execute("PRAGMA temp_store = MEMORY")

            cluster_pk_map = _lookup_pk_map(
                cur,
                cluster_ids,
                temp_table="tmp_lookup_keyword_clusters",
                source_table="Clusters",
                external_col="clusterID",
                pk_col="cluster_pk",
                strict=False,
            )

            missing = cluster_ids - set(cluster_pk_map)
            if missing:
                logger.warning(
                    "%d clusterIDs for keywords were not found. Examples: %s",
                    len(missing),
                    sorted(missing)[:10],
                )

            inserts = [
                (cluster_pk_map[cid], keyword)
                for keyword, cluster_ids_for_keyword in keyword_dict.items()
                for cid in cluster_ids_for_keyword
                if cid in cluster_pk_map
            ]

            for i in range(0, len(inserts), batch_size):
                cur.executemany(
                    "INSERT OR IGNORE INTO Keywords(cluster_pk, keyword) VALUES (?, ?)",
                    inserts[i : i + batch_size],
                )

        logger.info("Updated keywords with %d entries", len(inserts))

    except Exception:
        logger.exception("Error updating keywords")
        raise


def delete_keywords_from_csb(
    database: str, prefix: str = "csb-", suffix: str = "_"
) -> None:
    """
    Remove keywords from the database that match the pattern options.csb_name_prefix + a number + options.csb_name_suffix.

    Args:
        database: Name of the database to be worked on.
        prefix: prefix of the keyword to be deleted
        suffix: suffix of the keyword to be deleted
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Construct the pattern
        pattern = f"{prefix}%{suffix}"

        # SQL query to delete matching keywords
        delete_query = "DELETE FROM Keywords WHERE keyword LIKE ?"

        try:
            cur.execute(delete_query, (pattern,))
            con.commit()
        except sqlite3.Error as e:
            print("[ERROR] SQLite error:", e)
            raise

    return


def fetch_genome_ids(database: str) -> Set[str]:
    """
    Fetches the distinct genome_ids from the Genomes table in the SQLite database.

    Args:
        database: Path to the SQLite database.

    Returns:
        A set of distinct genome_ids.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT genomeID FROM Genomes")

        # Use set comprehension to create the set of genome_ids
        genome_ids = {row[0] for row in cur.fetchall()}

    return genome_ids


def fetch_genome_ids_with_proteins(database: str) -> Set[str]:
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        cur.execute("""
            SELECT DISTINCT g.genomeID
            FROM Proteins p
            JOIN Genomes g
              ON g.genome_pk = p.genome_pk
        """)

        return {row[0] for row in cur}


def clean_database_locks(database_path, wait_seconds=10):
    """
    Ensure the database is unlocked and all pending writes are flushed.

    Args:
        database_path: Path to the SQLite .db file
        wait_seconds: How long to wait if database is busy (default 10 seconds)
    """
    wal_path = database_path + "-wal"
    shm_path = database_path + "-shm"

    print(f"[INFO] Checking database for locks: {database_path}")

    # Step 1: Check if WAL and SHM files exist
    wal_exists = os.path.exists(wal_path)
    shm_exists = os.path.exists(shm_path)

    if wal_exists or shm_exists:
        print(f"[INFO] Found WAL/SHM files. WAL: {wal_exists}, SHM: {shm_exists}")
    else:
        print("[INFO] No WAL/SHM files found. Database seems clean.")

    # Step 2: Try to connect and perform a WAL checkpoint
    print("[INFO] Attempting to checkpoint database")
    success = False
    start_time = time.time()

    while not success and (time.time() - start_time) < wait_seconds:
        try:
            with sqlite3.connect(database_path, timeout=5) as con:
                # Force complete WAL checkpoint
                con.execute("PRAGMA wal_checkpoint(FULL);")
                con.commit()
                success = True
                print("[INFO] Checkpoint successful. Database flushed.")
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                print("[WARN] Database is locked. Waiting a bit...")
                time.sleep(1)
            else:
                print(f"Unexpected SQLite error: {e}")
                raise

    if not success:
        raise RuntimeError(
            f"Could not unlock the database after {wait_seconds} seconds."
        )

    # Step 3: Check again if WAL and SHM still exist
    wal_exists = os.path.exists(wal_path)
    shm_exists = os.path.exists(shm_path)

    if not wal_exists and not shm_exists:
        logger.info("WAL/SHM files cleaned up successfully.")
    else:
        logger.warning(
            "WAL/SHM files still exist. Database might have unclean shutdown earlier."
        )
    logger.info("Database cleaning routine finished.")


def fetch_genome_statistic(database):
    """
    18.10.22
        Args:
           database     Name of the database to be worked on
           filepath     Path to directory, will be extended by taxons
        Return:

        This one should return all Taxonomic groupings found in the database with number of associated genomes
    """

    taxons = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Ordnung",
        "Family",
        "Genus",
    ]  # Species left out because better to make own file for it
    con = sqlite3.connect(database)
    cur = con.cursor()

    for index, taxon in enumerate(taxons):
        query = (
            f"SELECT {taxon},count(*) FROM Genomes GROUP BY {taxon} ORDER BY {taxon}"
        )
        cur.execute(query)
        rows = cur.fetchall()

        try:
            writer = open(database + "_statistics_" + taxon, "w")
            writer.write(f"{taxon}\tcount[#]\n")
            for row in rows:
                print(row)
                writer.write(f"{row[0]}" + "\t" + f"{row[1]}\n")
        except (FileNotFoundError, IsADirectoryError, PermissionError) as e:
            # Pfad existiert nicht / ist ein Verzeichnis / fehlt Schreibrecht
            raise OSError(f"Cannot write statistics from '{database}': {e}") from e
        except OSError as e:
            # Sonstige OS-bezogene I/O-Fehler (IOError alias)
            raise OSError(f"I/O error while writing '{database}': {e}") from e
        else:
            writer.close()

    query = "SELECT * FROM Genomes ORDER BY Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species"
    cur.execute(query)
    rows = cur.fetchall()
    writer = open(database + "_statistics_Species", "w")
    names = list(map(lambda x: x[0], cur.description))
    writer.write("\t".join(names) + "\n")
    for row in rows:
        row = ["None" if v is None else str(v) for v in row]
        writer.write("\t".join(row) + "\n")

    writer.close()
    con.close()

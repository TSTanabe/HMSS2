#!/usr/bin/python
import sqlite3
import os
import sys
import traceback
import time
from typing import List, Dict, Set, Any

from hmsss.core.logging import get_logger

logger = get_logger(__name__)

########## Write output to Database Routines ##########


def create_database(database: str) -> None:
    """
    11.9.22
    Creating a database file for midterm storage of results

    Args:
        database    Pathway to database file
    """
    # create database
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""CREATE TABLE Genomes (
        genomeID    varchar(32)     PRIMARY KEY     NOT NULL,
        Superkingdom      varchar(128)              DEFAULT 'NULL',
        Clade       varchar(128)                    DEFAULT 'NULL',
        Phylum      varchar(128)                    DEFAULT 'NULL',
        Class       varchar(128)                    DEFAULT 'NULL',
        Ordnung     varchar(128)                    DEFAULT 'NULL',
        Family      varchar(128)                    DEFAULT 'NULL',
        Genus       varchar(128)                    DEFAULT 'NULL',
        Species     varchar(128)                    DEFAULT 'NULL',
        Strain      varchar(128)                    DEFAULT 'NULL',
        TypeStrain  tinyint(4)                      DEFAULT NULL,
        Completeness decimal(5,2)                   DEFAULT NULL,
        Contamination decimal(5,2)                  DEFAULT NULL,
        dRep        tinyint(1)                      DEFAULT NULL,
        NCBITaxon   int(11)                         DEFAULT NULL,
        NCBIProject int(11)                         DEFAULT NULL,
        NCBIBioproject varchar(32)                  DEFAULT NULL,
        NCBIBiosample varchar(32)                   DEFAULT NULL,
        NCBIAssembly varchar(32)                    DEFAULT NULL
        );""")

        cur.execute("""CREATE TABLE Clusters (
        clusterID   varchar(32)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(32)                     NOT NULL,
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE
        );""")

        cur.execute("""CREATE TABLE Keywords (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        clusterID   varchar(32)                     NOT NULL,
        keyword     varchar(32)                     NOT NULL,
        completeness    varchar(32)                 DEFAULT NULL,
        collinearity    varchar(32)                 DEFAULT NULL,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE CASCADE ON UPDATE CASCADE
        );""")

        cur.execute("""CREATE TABLE Proteins (
        proteinID   varchar(128)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(64)                     NOT NULL,
        clusterID   varchar(64)                     DEFAULT NULL, 
        locustag    varchar(32)                     DEFAULT NULL,
        contig      varchar(32)                     DEFAULT NULL,
        start       int(11)                         DEFAULT NULL,
        end         int(11)                         DEFAULT NULL,
        strand      varchar(1)                      DEFAULT NULL,
        comment     varchar(12)                     DEFAULT NULL,
        alternative_hit     varchar(64)                     DEFAULT NULL,
        dom_count   smallint(6)                     DEFAULT NULL,
        sequence    varchar(4096)                   DEFAULT NULL,
        UNIQUE(proteinID,genomeID),
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE SET NULL ON UPDATE CASCADE
        );""")

        cur.execute("""CREATE TABLE Domains (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        proteinID   varchar(32)                     NOT NULL,
        domain      varchar(32)                     DEFAULT NULL,
        score       smallint(6)                     DEFAULT NULL,
        blast_score_ratio       smallint(6)                     DEFAULT NULL,
        identity       smallint(6)                     DEFAULT NULL,
        domStart    int(11)                             DEFAULT NULL,
        domEnd      int(11)                             DEFAULT NULL,
        CONSTRAINT fk_proteinID FOREIGN KEY (proteinID) REFERENCES Proteins(proteinID) ON DELETE CASCADE ON UPDATE CASCADE
        );""")

    con.commit()
    con.close()

    return


def index_database(database: str) -> None:
    """
    12.11.22
    Indexes the relevant columns in the database to improve search and join performance.
    This is especially important for large databases (> 3000 genomes).
    """
    logger.info(f"Indexing database {database}")

    # List of indexes to create, each represented as a tuple (index_name, create_statement)
    indexes = [
        (
            "tab_dom_pid_index",
            "CREATE INDEX IF NOT EXISTS tab_dom_pid_index ON Domains(proteinID)",
        ),
        (
            "tab_prot_gid_index",
            "CREATE INDEX IF NOT EXISTS tab_prot_gid_index ON Proteins(genomeID)",
        ),
        (
            "tab_prot_cid_index",
            "CREATE INDEX IF NOT EXISTS tab_prot_cid_index ON Proteins(clusterID)",
        ),
        (
            "tab_key_cid_index",
            "CREATE INDEX IF NOT EXISTS tab_key_cid_index ON Keywords(clusterID)",
        ),
        (
            "tab_clus_gid_index",
            "CREATE INDEX IF NOT EXISTS tab_clus_gid_index ON Clusters(genomeID)",
        ),
        (
            "tab_dom_did_index",
            "CREATE INDEX IF NOT EXISTS tab_dom_did_index ON Domains(domain)",
        ),
        (
            "tab_key_kid_index",
            "CREATE INDEX IF NOT EXISTS tab_key_kid_index ON Keywords(keyword)",
        ),
    ]

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")

            # Start a transaction for batch index creation
            with con:
                for index_name, create_stmt in indexes:
                    cur.execute(create_stmt)

    except sqlite3.Error as e:
        logger.error(f"Error while indexing the database: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during indexing: {e}")
        sys.exit(1)


def insert_database_genome_ids(database: str, genome_ids: Set[str]) -> None:
    """
    22.6.24
        Args:
            Database    Name of the database to be appended
            GenomeIDs   List of genomeIDs to be added to the Genomes table
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""PRAGMA synchronous = OFF;""")
        cur.execute("""PRAGMA journal_mode = OFF;""")
        # Prepare a list of tuples, each containing one genomeID
        genome_id_tuples = [(genomeID,) for genomeID in genome_ids]
        # Use executemany to insert all genomeIDs in a single batch
        cur.executemany(
            """INSERT OR IGNORE INTO Genomes (genomeID) VALUES (?)""", genome_id_tuples
        )
        con.commit()
    con.close()
    return


def insert_database_proteins(database: str, protein_dict: Dict[str, Any]) -> None:
    """
    Inserts for concatenated glob hmm searches. GenomeId must be defined within the protein object
    Args:
        database (str): Path to database file.
        protein_dict (dict): Dictionary of protein objects. protein_id => obj
    """

    try:
        with sqlite3.connect(database) as con:
            # Inside 'with con', if an exception occurs, it will rollback automatically
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""PRAGMA synchronous = OFF;""")
            cur.execute("""PRAGMA journal_mode = OFF;""")

            protein_list = sorted(
                protein_dict.values(), key=lambda x: (x.gene_contig, x.gene_start)
            )
            protein_records = []
            domain_records = []

            for protein in protein_list:
                genome_id = protein.genomeID
                protein_id = f"{genome_id}-{protein.proteinID}"
                domains = protein.domains
                if not genome_id or not protein_id:
                    print(f"[WARN] Undefined error for {protein.proteinID}")
                    continue

                # Prepare the protein record
                protein_record = (
                    protein_id,
                    genome_id,
                    protein.gene_locustag,
                    protein.gene_contig,
                    protein.gene_start,
                    protein.gene_end,
                    protein.gene_strand,
                    protein.get_selection_comment_csv(),
                    protein.alternative_hit,
                    len(domains),
                    protein.get_sequence(),
                )
                protein_records.append(protein_record)

                # Prepare domain records
                for domain in domains.values():
                    domain_record = (
                        protein_id,
                        domain.HMM,
                        domain.start,
                        domain.end,
                        domain.score,
                        domain.identity,
                        domain.bsr,
                    )
                    domain_records.append(domain_record)

            # Batch insert for proteins
            cur.executemany(
                """INSERT OR IGNORE INTO Proteins
                (proteinID, genomeID, locustag, contig, start, end, strand, comment, alternative_hit, dom_count, sequence)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                protein_records,
            )

            # Batch insert for domains
            cur.executemany(
                """INSERT OR IGNORE INTO Domains
                (proteinID, domain, domStart, domEnd, score, identity, blast_score_ratio)
                VALUES (?, ?, ?, ?, ?, ?, ?)""",
                domain_records,
            )

            # No need to manually call con.commit() — 'with' will handle it

    except Exception as e:
        logger.warning(
            f"Proteins were not inserted for {genome_id}.\nError: - {str(e)}\nTraceback: {traceback.format_exc()}"
        )

    return


def insert_database_clusters(database: str, cluster_dict: Dict[str, Any]) -> None:
    """
    Inserts clusters and associated proteins/keywords.

    Args:
        database (str): Path to database file.
        cluster_dict (dict): Cluster objects.
    """

    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""PRAGMA synchronous = OFF;""")
            cur.execute("""PRAGMA journal_mode = OFF;""")

            cluster_inserts = []
            protein_updates = []
            keyword_inserts = []

            for cluster_index, cluster in cluster_dict.items():
                cluster_id = cluster.get_cluster_id()
                if cluster_id is None:
                    continue

                cluster_inserts.append((cluster_id, cluster.genomeID))

                for proteinID in cluster.get_genes():
                    full_protein_id = f"{cluster.genomeID}-{proteinID}"
                    protein_updates.append((cluster_id, full_protein_id))

                for Keyword in cluster.get_keywords():
                    keyword_inserts.append(
                        (
                            cluster_id,
                            Keyword.get_keyword(),
                            Keyword.get_completeness(),
                            Keyword.get_csb(),
                        )
                    )

            # Batch insert clusters
            cur.executemany(
                """INSERT OR IGNORE INTO Clusters (clusterID, genomeID) VALUES (?, ?)""",
                cluster_inserts,
            )

            # Batch update proteins
            cur.executemany(
                """UPDATE Proteins SET clusterID = ? WHERE proteinID = ?""",
                protein_updates,
            )

            # Batch insert or replace keywords
            cur.executemany(
                """INSERT OR REPLACE INTO Keywords (clusterID, keyword, completeness, collinearity) VALUES (?, ?, ?, ?)""",
                keyword_inserts,
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
        raise ValueError(
            "Line must contain at least two tab-separated fields: <ID> and <taxonomy>"
        )

    genome_id = parts[0].strip()
    tax_str = parts[1].strip()

    raw_tokens = [t.strip() for t in tax_str.split(";") if t.strip()]
    ranks = {k: na for k in rank_keys}

    prefix_to_rank = {
        "d__": "domain",
        "k__": "domain",
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
########## Alter information from database routines ##########
##############################################################


def update_domain(
    database: str, protein_diction: Dict[str, Any], old_tag: str, new_tag: str
) -> None:
    """
    18.11.22
        Args:
           database     Name of the database to be worked on
           protein_diction Dictionary of protein IDs to be updated
           old_tag     The old domain tag that needs to be updated
           new_tag     The new domain tag that will replace the old tag
    """
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            query = (
                """UPDATE Domains SET domain = ? WHERE proteinID = ? AND domain = ?"""
            )
            data = [(new_tag, proteinID, old_tag) for proteinID in protein_diction]
            cur.executemany(query, data)
        logger.info(
            f"Updated {len(protein_diction)} domains from {old_tag} to {new_tag}"
        )
    except Exception as e:
        logger.error(f"Error updating domains: {e}")
    return


##############################################################
########## Fetch information from database routines ##########
##############################################################


def update_keywords(
    database: str, keyword_dict: Dict[str, Set[str]], batch_size: int = 400
) -> None:
    """
    Update keywords in the database in batches.

    Args:
        database: Name of the database to be worked on.
        keyword_dict: Dictionary containing the new keywords and associated clusterIDs.
        batch_size: Number of (clusterID, keyword) pairs per batch insert.
    """
    if keyword_dict:
        try:
            with sqlite3.connect(database) as con:
                cur = con.cursor()
                query = """INSERT INTO Keywords (clusterID, keyword) VALUES (?, ?)"""
                inserts = []
                for new_keyword, clusterIDs in keyword_dict.items():
                    for clusterID in clusterIDs:
                        inserts.append((clusterID, new_keyword))
                for i in range(0, len(inserts), batch_size):
                    batch = inserts[i : i + batch_size]
                    cur.executemany(query, batch)
            logger.info(f"Updated keywords with {len(inserts)} entries.")
        except Exception as e:
            logger.error(f"Error updating keywords: {e}")
    return


def delete_keywords_from_csb(database: str, prefix: str = "csb-", suffix: str = "_") -> None:
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

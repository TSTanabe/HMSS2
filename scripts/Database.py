#!/usr/bin/python
import sqlite3
import os
import re
import traceback
from typing import Set, List

from . import myUtil
from . import ParseReports
from . import Csb_finder

"""
This module keeps all routines for database input and working routines on the database
if looking for output routines invoked by the user got to Output.py module
"""

########## Write output to Database Routines ##########

def create_database(database: str) -> None:
    """
    11.9.22
    Creates a new SQLite database file with the required schema for storing
    genome, cluster, keyword, protein, and domain information.

    Args:
        database: Path to the SQLite database file to create.

    Raises:
        sqlite3.Error: If there is an error executing any SQL statement.
        OSError:     If the database file cannot be created due to filesystem issues.
    """
    schema_statements = [
        """
        PRAGMA foreign_keys = ON;
        """,
        """
        CREATE TABLE Genomes (
            genomeID        VARCHAR(32)  PRIMARY KEY NOT NULL,
            Superkingdom    VARCHAR(128) DEFAULT 'NULL',
            Clade           VARCHAR(128) DEFAULT 'NULL',
            Phylum          VARCHAR(128) DEFAULT 'NULL',
            Class           VARCHAR(128) DEFAULT 'NULL',
            Ordnung         VARCHAR(128) DEFAULT 'NULL',
            Family          VARCHAR(128) DEFAULT 'NULL',
            Genus           VARCHAR(128) DEFAULT 'NULL',
            Species         VARCHAR(128) DEFAULT 'NULL',
            Strain          VARCHAR(128) DEFAULT 'NULL',
            TypeStrain      TINYINT(4)   DEFAULT NULL,
            Completeness    DECIMAL(5,2) DEFAULT NULL,
            Contamination   DECIMAL(5,2) DEFAULT NULL,
            dRep            TINYINT(1)   DEFAULT NULL,
            NCBITaxon       INT(11)      DEFAULT NULL,
            NCBIProject     INT(11)      DEFAULT NULL,
            NCBIBioproject  VARCHAR(32)  DEFAULT NULL,
            NCBIBiosample   VARCHAR(32)  DEFAULT NULL,
            NCBIAssembly    VARCHAR(32)  DEFAULT NULL
        );
        """,
        """
        CREATE TABLE Clusters (
            clusterID    VARCHAR(32) NOT NULL PRIMARY KEY,
            genomeID     VARCHAR(32) NOT NULL,
            CONSTRAINT fk_genomeID
                FOREIGN KEY (genomeID)
                REFERENCES Genomes(genomeID)
                ON DELETE CASCADE
                ON UPDATE CASCADE
        );
        """,
        """
        CREATE TABLE Keywords (
            clusterID     VARCHAR(32) NOT NULL,
            keyword       VARCHAR(32) NOT NULL,
            completeness  VARCHAR(32) DEFAULT NULL,
            collinearity  VARCHAR(32) DEFAULT NULL,
            PRIMARY KEY (clusterID, keyword),
            CONSTRAINT fk_clusterID
                FOREIGN KEY (clusterID)
                REFERENCES Clusters(clusterID)
                ON DELETE CASCADE
                ON UPDATE CASCADE
        );
        """,
        """
        CREATE TABLE Proteins (
            proteinID    VARCHAR(128) NOT NULL PRIMARY KEY,
            genomeID     VARCHAR(64)  NOT NULL,
            clusterID    VARCHAR(64)  DEFAULT NULL,
            locustag     VARCHAR(32)  DEFAULT NULL,
            contig       VARCHAR(32)  DEFAULT NULL,
            start        INT(11)      DEFAULT NULL,
            end          INT(11)      DEFAULT NULL,
            strand       VARCHAR(1)   DEFAULT NULL,
            dom_count    SMALLINT(6)  DEFAULT NULL,
            sequence     VARCHAR(4096) DEFAULT NULL,
            UNIQUE(proteinID, genomeID),
            CONSTRAINT fk_genomeID
                FOREIGN KEY (genomeID)
                REFERENCES Genomes(genomeID)
                ON DELETE CASCADE
                ON UPDATE CASCADE,
            CONSTRAINT fk_clusterID
                FOREIGN KEY (clusterID)
                REFERENCES Clusters(clusterID)
                ON DELETE SET NULL
                ON UPDATE CASCADE
        );
        """,
        """
        CREATE TABLE Domains (
            proteinID  VARCHAR(128) NOT NULL,
            domain     VARCHAR(32)  DEFAULT NULL,
            score      SMALLINT(6)  DEFAULT NULL,
            domStart   INT(11)      DEFAULT NULL,
            domEnd     INT(11)      DEFAULT NULL,
            PRIMARY KEY (proteinID, domain, domStart, domEnd),
            CONSTRAINT fk_proteinID
                FOREIGN KEY (proteinID)
                REFERENCES Proteins(proteinID)
                ON DELETE CASCADE
                ON UPDATE CASCADE
        );
        """
    ]

    try:
        with sqlite3.connect(database) as conn:
            cursor = conn.cursor()
            for statement in schema_statements:
                cursor.execute(statement)
    except sqlite3.Error as exc:
        raise sqlite3.Error(f"SQLite error while creating schema: {exc}") from exc
    except OSError as exc:
        raise OSError(f"Filesystem error creating database '{database}': {exc}") from exc

    
        
    return 
    
def index_database(database):
    """
    12.11.22
    Indexes the relevant columns in the database to improve search and join performance.
    This is especially important for large databases (> 3000 genomes).
    """
    print(f"Indexing database {database}")
    
    # List of indexes to create, each represented as a tuple (index_name, create_statement)
    indexes = [
        ("tab_dom_pid_index", "CREATE INDEX IF NOT EXISTS tab_dom_pid_index ON Domains(proteinID)"),
        ("tab_prot_gid_index", "CREATE INDEX IF NOT EXISTS tab_prot_gid_index ON Proteins(genomeID)"),
        ("tab_prot_cid_index", "CREATE INDEX IF NOT EXISTS tab_prot_cid_index ON Proteins(clusterID)"),
        ("tab_key_cid_index", "CREATE INDEX IF NOT EXISTS tab_key_cid_index ON Keywords(clusterID)"),
        ("tab_clus_gid_index", "CREATE INDEX IF NOT EXISTS tab_clus_gid_index ON Clusters(genomeID)"),
        ("tab_dom_did_index", "CREATE INDEX IF NOT EXISTS tab_dom_did_index ON Domains(domain)"),
        ("tab_key_kid_index", "CREATE INDEX IF NOT EXISTS tab_key_kid_index ON Keywords(keyword)")
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
        print(f"An error occurred while indexing the database: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    print("Indexing complete.")




def drop_specified_indices(database: str) -> None:
    """
    Drops predefined indices from the SQLite database if they exist and
    are not tied to UNIQUE or PRIMARY KEY constraints.

    Args:
        database: Path to the SQLite database file.
    """
    indices_to_drop: List[str] = [
        "tab_dom_pid_index",
        "tab_prot_gid_index",
        "tab_prot_cid_index",
        "tab_key_cid_index",
        "tab_clus_gid_index",
        "tab_dom_did_index",
        "tab_key_kid_index",
    ]

    try:
        with sqlite3.connect(database) as conn:
            cursor = conn.cursor()

            # Build query with placeholders for index names
            placeholders = ", ".join("?" for _ in indices_to_drop)
            query = f"""
                SELECT name, sql
                FROM sqlite_master
                WHERE type = 'index'
                  AND name IN ({placeholders});
            """
            cursor.execute(query, indices_to_drop)
            indices = cursor.fetchall()

            for index_name, index_sql in indices:
                # Skip indices tied to UNIQUE or PRIMARY KEY constraints
                if "UNIQUE" in index_sql or "PRIMARY KEY" in index_sql:
                    print(f"Skipping index '{index_name}' (associated with UNIQUE or PRIMARY KEY)")
                else:
                    cursor.execute(f"DROP INDEX IF EXISTS {index_name};")

            print("Finished dropping specified indices.")
    except sqlite3.Error as exc:
        print(f"An error occurred while dropping indices: {exc}")


    

def insert_database_genomeIDs(database, genomeIDs):
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
        genomeID_tuples = [(genomeID,) for genomeID in genomeIDs]
        # Use executemany to insert all genomeIDs in a single batch
        cur.executemany('''INSERT OR IGNORE INTO Genomes (genomeID) VALUES (?)''', genomeID_tuples)
        con.commit()
    con.close()
    return
    
    
        
def insert_database_proteins(database, protein_dict):
    """
        Inserts for concated glob hmmsearches. GenomeId must be defined within the protein object
        Args:
            protein_dict    dictionary with protein objects
    """
    with sqlite3.connect(database) as con:
        try:
            cur = con.cursor()
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""PRAGMA synchronous = OFF;""")
            cur.execute("""PRAGMA journal_mode = OFF;""")
            
            protein_list = sorted(protein_dict.values(), key=lambda x: (x.gene_contig, x.gene_start))
            #protein_list = protein_dict.values()
            protein_records = []
            domain_records = []

            for protein in protein_list: #protein_dict.values()
                genomeID = protein.genomeID
                proteinID = f"{genomeID}-{protein.proteinID}"
                domains = protein.domains
                if not genomeID or not proteinID:
                    print(f"Warning error annotated protein {protein.proteinID}")
                    continue
                # Prepare the protein record
                protein_record = (
                    proteinID, genomeID, protein.gene_locustag, protein.gene_contig,
                    protein.gene_start, protein.gene_end, protein.gene_strand,
                    len(domains), protein.get_sequence()
                )
                protein_records.append(protein_record)
                # Prepare domain records
                for domain in domains.values():
                    domain_record = (
                        proteinID, domain.HMM, domain.start, domain.end, domain.score
                    )
                    domain_records.append(domain_record)
            
            # Batch insert for proteins
            cur.executemany('''INSERT OR IGNORE INTO Proteins
                (proteinID, genomeID, locustag, contig, start, end, strand, dom_count, sequence)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''', protein_records)
            
            # Batch insert for domains
            cur.executemany('''INSERT OR IGNORE INTO Domains
                (proteinID, domain, domStart, domEnd, score)
                VALUES (?, ?, ?, ?, ?)''', domain_records)
            
            con.commit()
        except Exception as e:
                error_message = f"\nError occurred: {str(e)}"
                traceback_details = traceback.format_exc()
                print(f"\tWARNING: Due to an error - {error_message}")
                print(f"\tTraceback details:\n{traceback_details}")
                
    return
    

def insert_database_clusters(database, cluster_dict):
    """
    22.06.2024
        Args:
            database - path to the database file
            cluster_dict - dictionary with cluster objects
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
                clusterID = cluster.get_clusterID()
                if clusterID is None:
                    continue

                cluster_inserts.append((clusterID, cluster.genomeID))
                
                for proteinID in cluster.get_genes():
                    full_proteinID = f"{cluster.genomeID}-{proteinID}"
                    protein_updates.append((clusterID, full_proteinID))

                for Keyword in cluster.get_keywords():
                    keyword_inserts.append((clusterID, Keyword.get_keyword(), Keyword.get_completeness(), Keyword.get_csb()))

            # Batch insert clusters
            cur.executemany("""INSERT OR IGNORE INTO Clusters (clusterID, genomeID) VALUES (?, ?)""", cluster_inserts)

            # Batch update proteins
            cur.executemany("""UPDATE Proteins SET clusterID = ? WHERE proteinID = ?""", protein_updates)

            # Batch insert or replace keywords
            cur.executemany("""INSERT OR REPLACE INTO Keywords (clusterID, keyword, completeness, collinearity) VALUES (?, ?, ?, ?)""", keyword_inserts)
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\tWARNING: Due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
 

    return


def insert_taxonomy_data(database, taxonomy_file):
    """
    Insert taxonomy data from a file into the Genomes table of the database.

    Args:
        database: Path to the SQLite database file.
        taxonomy_file: Path to the parsed taxonomy file (tab-separated).
    """
    try:
        # Check if the database file path is valid
        if not os.path.exists(os.path.dirname(database)):
            raise FileNotFoundError(f"Directory for database does not exist: {os.path.dirname(database)}")

        # Connect to the SQLite database
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            
            # Read the taxonomy file and insert data into the Genomes table
            with open(taxonomy_file, 'r') as file:
                # Skip the header
                next(file)
                
                # Read each line in the file
                for line in file:
                    fields = line.strip().split('\t')
                   
                    if len(fields) < 8:
                        print(f"Skipping line due to insufficient columns: {line}")
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
                    cur.execute('''
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
                    ''', (genome_id, superkingdom, clade, phylum, class_, order, family, genus, species))
            
            # Commit the transaction
            con.commit()
    
    except sqlite3.OperationalError as e:
        print(f"SQLite Operational Error: {e}")
        print(f"Database path: {database}")
        print(f"Ensure the file exists and you have read/write permissions.")
    except sqlite3.Error as e:
        print(f"SQLite Error: {e}")
    except FileNotFoundError as e:
        print(f"File Not Found Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
        
        
##############################################################
########## Alter information from database routines ##########
##############################################################

def update_domain(database, protein_diction, old_tag, new_tag):
    """
    18.11.22
        Args:
           database     Name of the database to be worked on
           protein_diction Dictionary of protein IDs to be updated
           old_tag     The old domain tag that needs to be updated
           new_tag     The new domain tag that will replace the old tag
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = """UPDATE Domains SET domain = ? WHERE proteinID = ? AND domain = ?"""
        
        # Prepare the data for bulk update
        data = [(new_tag, proteinID, old_tag) for proteinID in protein_diction]
        
        # Use executemany for batch processing
        cur.executemany(query, data)
        
        # Commit the transaction (optional, as it is done automatically with 'with' context)
        con.commit()

    return

##############################################################
########## Fetch information from database routines ##########
##############################################################


def fetch_genomeID_list(database):
    """
    10.11.22
        Args:
           database     Name of the database to be worked on
        Return:
           dictionary genomeID => 1
    """
    genomeIDs = []
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""SELECT genomeID FROM Genomes;""")
        
        for count,c in enumerate(cur):
            print(f"\tSelected genomeIDs {count}",end="\r")
            genomeIDs.append(c[0])
        print(f"\tSelected genomeIDs {count} for csb finder")
    return genomeIDs








def fetch_genome_statistic(database):
    """
    18.10.22
        Args:
           database     Name of the database to be worked on
           filepath     Path to directory, will be extended by taxons
        Return:
        
        This one should return all Taxonomic groupings found in the database with number of associated genomes
    """
    
    taxons = ["Superkingdom","Phylum","Class","Ordnung","Family","Genus"] #Species left out because better to make own file for it
    con = sqlite3.connect(database)
    cur=con.cursor()
    
    
    for index,taxon in enumerate(taxons):
        query = f"SELECT {taxon},count(*) FROM Genomes GROUP BY {taxon} ORDER BY {taxon}"
        cur.execute(query)
        rows = cur.fetchall()

        try:
            writer = open(database+"_statistics_"+taxon, "w")
            writer.write(f"{taxon}\tcount[#]\n")
            for row in rows:
                print(row)
                writer.write(f"{row[0]}"+"\t"+f"{row[1]}\n")
        except:
            print("ERROR: Could not open write database statistics")
        else:
            writer.close()
    
    
    query = "SELECT * FROM Genomes ORDER BY Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species"
    cur.execute(query)
    rows = cur.fetchall()
    writer = open(database+"_statistics_Species", "w")
    names = list(map(lambda x: x[0], cur.description))
    writer.write("\t".join(names)+"\n")
    for row in rows:

        row = ['None' if v is None else str(v) for v in row]
        writer.write("\t".join(row)+"\n")
        
    writer.close()    
    con.close()
    
    




def update_keywords(database, keyword_dict):
    """
    Update keywords in the database.

    Args:
        database: Name of the database to be worked on.
        keyword_dict: Dictionary containing the clusterID as key and the new keyword as value.
    """
    if keyword_dict:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            query = """INSERT OR IGNORE INTO Keywords (clusterID, keyword) VALUES (?, ?)"""
            for new_keyword,clusterIDs in keyword_dict.items():
                for clusterID in clusterIDs:
                    cur.execute(query, (clusterID, new_keyword))

    return




def delete_keywords_from_csb(database, options):
    """
    Remove keywords from the database that match the pattern options.csb_name_prefix + a number + options.csb_name_suffix.

    Args:
        database: Name of the database to be worked on.
        options: An object containing csb_name_prefix and csb_name_suffix attributes.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        
        # Construct the pattern
        pattern = f"{options.csb_name_prefix}%{options.csb_name_suffix}"
        
        # SQL query to delete matching keywords
        delete_query = "DELETE FROM Keywords WHERE keyword LIKE ?"
        
        try:
            cur.execute(delete_query, (pattern,))
            con.commit()
        except sqlite3.Error as e:
            print("SQLite error:", e)
            raise

    return










def fetch_genomeIDs_from_proteins(database):
    """
    Fetches the distinct genomeIDs from the Genomes table in the SQLite database.
    
    Args:
        database: Path to the SQLite database.
    
    Returns:
        A set of distinct genomeIDs.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT genomeID FROM Proteins")
        
        # Use set comprehension to create the set of genomeIDs
        genomeIDs = {row[0] for row in cur.fetchall()}
    
    return genomeIDs

            
def fetch_database_all(database):
    """
    1.10.22
        Args:
            ONLY FOR DEBUGGING
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" SELECT Proteins.proteinID,domain,clusterID FROM Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE genomeID =  "GB_GCA_000016605" """)
        con.commit()
        print("----------Proteins--------------")
        print(cur.fetchall())
        # WHERE genomeID = GB_GCA_000016605
        cur = con.cursor()
        cur.execute(""" SELECT * FROM Clusters WHERE genomeID = "GB_GCA_000016605" """)
        con.commit()
        print("----------Clusters--------------")
        print(cur.fetchall())


        cur = con.cursor()
        cur.execute(""" SELECT * FROM Keywords  LEFT JOIN Clusters ON Keywords.clusterID = Clusters.clusterID WHERE genomeID = "GB_GCA_000016605" """)
        con.commit()
        print("----------Keywords--------------")
        print(cur.fetchall())
        
        cur = con.cursor()
        cur.execute(""" SELECT * FROM Genomes LIMIT 10 """)
        con.commit()
        print("----------Genomes--------------")
        #print(cur.fetchall())
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" SELECT DISTINCT Genomes.genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID LEFT JOIN Genomes ON Proteins.genomeID = Genomes.genomeID LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID LIMIT 100 """)
        con.commit()
        print("----------Proteins--------------")
        #print(cur.fetchall())
    return






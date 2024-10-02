#!/usr/bin/python
import sqlite3
import os
import re
import traceback


from . import myUtil
from . import ParseReports
from . import Csb_finder

"""
This module keeps all routines for database input and working routines on the database
if looking for output routines invoked by the user got to Output.py module
"""

########## Write output to Database Routines ##########

def create_database(database):
    """
    11.9.22
    Creating a database file for midterm storage of results
    
    Args:
        database    Pathway to database file
    """
    #create database
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute('''CREATE TABLE Genomes (
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
        );''')
        
        cur.execute('''CREATE TABLE Clusters (
        clusterID   varchar(32)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(32)                     NOT NULL,
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Keywords (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        clusterID   varchar(32)                     NOT NULL,
        keyword     varchar(32)                     NOT NULL,
        completeness    varchar(32)                 DEFAULT NULL,
        collinearity    varchar(32)                 DEFAULT NULL,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Proteins (
        proteinID   varchar(128)     PRIMARY KEY     NOT NULL,
        genomeID    varchar(64)                     NOT NULL,
        clusterID   varchar(64)                     DEFAULT NULL, 
        locustag    varchar(32)                     DEFAULT NULL,
        contig      varchar(32)                     DEFAULT NULL,
        start       int(11)                         DEFAULT NULL,
        end         int(11)                         DEFAULT NULL,
        strand      varchar(1)                      DEFAULT NULL,
        dom_count   smallint(6)                     DEFAULT NULL,
        sequence    varchar(4096)                   DEFAULT NULL,
        UNIQUE(proteinID,genomeID),
        CONSTRAINT fk_genomeID FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE,
        CONSTRAINT fk_clusterID FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE SET NULL ON UPDATE CASCADE
        );''')
        
        cur.execute('''CREATE TABLE Domains (
        ID          integer         PRIMARY KEY     AUTOINCREMENT,
        proteinID   varchar(32)                     NOT NULL,
        domain      varchar(32)                     DEFAULT NULL,
        score       smallint(6)                     DEFAULT NULL,
        domStart    int(11)                             DEFAULT NULL,
        domEnd      int(11)                             DEFAULT NULL,
        CONSTRAINT fk_proteinID FOREIGN KEY (proteinID) REFERENCES Proteins(proteinID) ON DELETE CASCADE ON UPDATE CASCADE
        );''')
        
              
    con.commit()
    con.close()
        #res = cur.execute("Select name from sqlite_master")
        #print(res.fetchall())
    
        
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

def light_index_database(database):
    """
    12.11.22
    Indexes the relevant columns in the database to improve search and join performance.
    This is especially important for large databases (> 3000 genomes).
    """
    
    # List of indexes to create, each represented as a tuple (index_name, create_statement)
    indexes = [
        ("tab_dom_pid_index", "CREATE INDEX IF NOT EXISTS tab_dom_pid_index ON Domains(proteinID)"),
        ("tab_prot_gid_index", "CREATE INDEX IF NOT EXISTS tab_prot_gid_index ON Proteins(genomeID)"),
        ("tab_clus_gid_index", "CREATE INDEX IF NOT EXISTS tab_clus_gid_index ON Clusters(genomeID)")
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

def drop_specified_indices(database):
    """
    Drops specified indices from the SQLite database if they exist and are not associated with UNIQUE or PRIMARY KEY constraints.

    Args:
        database (str): Path to the SQLite database file.
        indices_to_drop (list): List of index names that should be dropped.
    """
    indices_to_drop = [
    "tab_dom_pid_index",
    "tab_prot_gid_index",
    "tab_prot_cid_index",
    "tab_key_cid_index",
    "tab_clus_gid_index",
    "tab_dom_did_index",
    "tab_key_kid_index"
    ]
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            
            # Query to get information about the specified indices
            cur.execute(f"""
            SELECT name, sql 
            FROM sqlite_master 
            WHERE type = 'index' 
            AND name IN ({','.join('?' for _ in indices_to_drop)});
            """, indices_to_drop)
            
            indices = cur.fetchall()
            
            for index_name, index_sql in indices:
                # Check if the index is associated with UNIQUE or PRIMARY KEY constraints
                if "UNIQUE" in index_sql or "PRIMARY KEY" in index_sql:
                    print(f"Skipping index '{index_name}' (associated with UNIQUE or PRIMARY KEY constraint)")
                else:
                    #print(f"Dropping index: {index_name}")
                    cur.execute(f"DROP INDEX IF EXISTS {index_name};")
            
            print("Finished dropping specified indices.")

    except sqlite3.Error as e:
        print(f"An error occurred while dropping indices: {e}")


    
def insert_database_genomeID(database,genomeID):
    """
    11.9.22
        Args: 
            Database    Name of the database to be appended
            GenomeID    genomeID to be added to the Genomes table
            
        
    """
    with sqlite3.connect(database) as con:
        
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(''' INSERT OR IGNORE INTO Genomes (genomeID) VALUES (?) ON CONFLICT(genomeID) DO NOTHING;''',[genomeID])
    con.commit()
    con.close()
    return

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
        # Prepare a list of tuples, each containing one genomeID
        genomeID_tuples = [(genomeID,) for genomeID in genomeIDs]
        # Use executemany to insert all genomeIDs in a single batch
        cur.executemany('''INSERT OR IGNORE INTO Genomes (genomeID) VALUES (?)''', genomeID_tuples)
        con.commit()
    con.close()
    return
    
def insert_database_protein_deprecated(database,genomeID,protein_dict):
    """
    1.10.22
        Args:
            protein_dict    dictionary with protein objects
            
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        
        
        protein_list = sorted(protein_dict, key=lambda x: \
        (protein_dict[x].gene_contig, protein_dict[x].gene_start)) 
            
        for key in protein_list:
            protein = protein_dict.get(key)
            proteinID = protein.proteinID
            proteinID = f"{genomeID}-{proteinID}"   #genomeID added to proteinID to deal with the multispecies information
            domains = protein.get_domains_dict()
            
            
            #print(proteinID, genomeID, protein.get_gene_locustag(), protein.get_gene_contig(), protein.get_gene_start(), protein.get_gene_end(), protein.get_gene_strand(), protein.get_domain_count(), protein.get_sequence())    
            cur.execute('''INSERT OR IGNORE INTO Proteins
            (proteinID,genomeID,locustag,contig,start,end,strand,dom_count,sequence)
            VALUES (?,?,?,?,?,?,?,?,?) ''',\
            (proteinID, genomeID, protein.gene_locustag, protein.gene_contig,\
             protein.gene_start, protein.gene_end, protein.gene_strand,\
             protein.get_domain_count(), protein.get_sequence())\
             )
            
            for domain_index in domains:
                domain = domains.get(domain_index)
                cur.execute(""" INSERT OR IGNORE INTO Domains 
                                (proteinID,domain,domStart,domEnd,score) 
                                VALUES (?,?,?,?,?)""",\
                                (proteinID,domain.get_HMM(),domain.get_start(),\
                                domain.get_end(),domain.get_score())\
                                )
            
    con.commit()
    con.close()
    return
    
def insert_database_protein(database, genomeID, protein_dict):
    """
    1.10.22
    Inserts protein data into the database, including associated domains.
    
    Args:
        database: Path to the SQLite database.
        genomeID: The ID of the genome to which the proteins belong.
        protein_dict: Dictionary containing protein objects.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        
        # Sort the proteins by contig and start position
        protein_list = sorted(protein_dict, key=lambda x: 
                              (protein_dict[x].gene_contig, protein_dict[x].gene_start))
        
        # Lists to store batch insert data
        protein_data = []
        domain_data = []
        
        for key in protein_list:
            protein = protein_dict[key]
            proteinID = f"{genomeID}-{protein.proteinID}"  # Add genomeID to proteinID
            
            # Cache frequently accessed attributes to minimize redundant lookups
            gene_locustag = protein.gene_locustag
            gene_contig = protein.gene_contig
            gene_start = protein.gene_start
            gene_end = protein.gene_end
            gene_strand = protein.gene_strand
            dom_count = protein.get_domain_count()
            sequence = protein.get_sequence()
            domains = protein.get_domains_dict()
            
            # Append protein data to batch list
            protein_data.append((proteinID, genomeID, gene_locustag, gene_contig,
                                 gene_start, gene_end, gene_strand, dom_count, sequence))
            
            # Append domain data to batch list
            for domain_index in domains:
                domain = domains[domain_index]
                domain_data.append((proteinID, domain.get_HMM(), domain.get_start(),
                                    domain.get_end(), domain.get_score()))
        
        # Perform batch inserts for proteins and domains
        cur.executemany('''INSERT OR IGNORE INTO Proteins
                           (proteinID, genomeID, locustag, contig, start, end, strand, dom_count, sequence)
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''', protein_data)
        
        cur.executemany('''INSERT OR IGNORE INTO Domains
                           (proteinID, domain, domStart, domEnd, score)
                           VALUES (?, ?, ?, ?, ?)''', domain_data)
        
        # Commit all changes at once
        con.commit()
        
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
    
def update_protein_sequences(database, protein_dict):
    """
    Updates all protein sequences in the database from the protein objects in the hash.

    Args:
        database (str): Path to the database file.
        protein_dict (dict): Dictionary with protein objects.
    """
    update_data = []

    for proteinID, protein in protein_dict.items():
        sequence = protein.protein_sequence
        genomeID = protein.genomeID
        update_data.append((sequence, f"{genomeID}-{proteinID}"))

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        
        # Batch update protein sequences
        cur.executemany("""
            UPDATE Proteins
            SET sequence = ?
            WHERE proteinID = ?;
        """, update_data)

        con.commit()

    con.close()
    
    
def insert_database_cluster(database, genomeID, cluster_dict):
    """
    Inserts cluster data into the database, including associated proteins and keywords.

    Args:
        database: Path to the SQLite database.
        genomeID: The ID of the genome to which the clusters belong.
        cluster_dict: Dictionary containing cluster objects.
    """
    # TODO: Behavior -> Do not delete old entries on rename, only handle duplicates.
    # TODO: Add keyword selection behavior: if found, do nothing; otherwise, insert.
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")


        # Prepare lists for batch operations
        cluster_data = []
        protein_update_data = []
        keyword_data = []

        for cluster_index, cluster in cluster_dict.items():
            clusterID = cluster.clusterID
            if clusterID is None:
                continue

            # Insert cluster data
            cluster_data.append((clusterID, genomeID))

            # Update protein data with the clusterID
            cluster_proteins = cluster.get_genes()
            protein_update_data.extend([(clusterID, f"{genomeID}-{proteinID}") for proteinID in cluster_proteins])

            # Insert or update keyword data
            for Keyword in cluster.get_keywords():
                keyword_data.append((clusterID, Keyword.get_keyword(), Keyword.get_completeness(), Keyword.get_csb()))

        # Perform batch insert for Clusters
        cur.executemany("""INSERT OR IGNORE INTO Clusters (clusterID, genomeID) VALUES (?, ?)""", cluster_data)

        # Perform batch update for Proteins
        cur.executemany("""UPDATE Proteins SET clusterID = ? WHERE proteinID = ?""", protein_update_data)

        # Perform batch insert or replace for Keywords
        cur.executemany("""INSERT OR REPLACE INTO Keywords (clusterID, keyword, completeness, collinearity) 
                           VALUES (?, ?, ?, ?)""", keyword_data)


        # Commit all changes once after all operations
        con.commit()


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

  
    
        
def extend_database_protein(database,genomeID,protein_dict):
    """
    18.2.23
        Args:
            protein_dict    dictionary with possibly incomplete protein objects may be incomplete due to new protein domains
            
        Routine should extend proteins by domains and add new proteins to an existing proteins table
        Currently also overwrites data already existing. This is possibly more natural since the new HMM might be the only correct one
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        
        for key in protein_dict.keys():
            protein = protein_dict.get(key)
            proteinID = protein.get_proteinID()
            proteinID = f"{genomeID}-{proteinID}"   #genomeID added to proteinID to deal with the multispecies information
            domains = protein.get_domains_dict()
            new_domains_list = protein.get_domain_listing()
            cur.execute('SELECT EXISTS(SELECT 1 FROM Proteins WHERE proteinID = ?)', (proteinID,))
            result = cur.fetchone()[0]
            
            if result == 1:
                #'The value exists in the Proteins table' therefore add domain if possible
                #wiederherstellen des domain table aus datenbank
                #check if insert is possible, if yes insert into db else do nothing, but maybe an alert or protocol output
                cur.execute('SELECT domain,domStart,domEnd,score,ID FROM Domains WHERE proteinID = ?', (proteinID,))
                db_domains = cur.fetchall()
                
                IDs = []
                for db_domain in db_domains:
                    protein.add_domain(db_domain[0],db_domain[1],db_domain[2],db_domain[3])
                    IDs.append(db_domain[4])
                    
                annotated_domains = protein.get_domain_listing()
                if set(annotated_domains).intersection(set(new_domains_list)): # returns the new domains which are compatible
                    
                    domains = protein.get_domains_dict()
                    for domain_index in domains:
                        domain = domains.get(domain_index)
                        cur.execute(""" INSERT OR IGNORE INTO Domains 
                                    (proteinID,domain,domStart,domEnd,score) 
                                    VALUES (?,?,?,?,?)""",\
                                    (proteinID,domain.get_HMM(),domain.get_start(),\
                                    domain.get_end(),domain.get_score())\
                                    )
                    for ID in IDs:
                        cur.execute('DELETE FROM Domains WHERE ID = ?', (ID,))

                
            else:
                #'The value does not exist in the Proteins table' new protein found, add to the existing table
                cur.execute('''INSERT OR IGNORE INTO Proteins
                (proteinID,genomeID,locustag,contig,start,end,strand,dom_count,sequence)
                VALUES (?,?,?,?,?,?,?,?,?) ''',\
                (proteinID, genomeID, protein.gene_locustag, protein.gene_contig,\
                 protein.get_gene_start, protein.gene_end, protein.gene_strand,\
                 protein.get_domain_count(), protein.get_sequence())\
                 )
                
                for domain_index in domains:
                    domain = domains.get(domain_index)
                    cur.execute(""" INSERT OR IGNORE INTO Domains 
                                    (proteinID,domain,domStart,domEnd,score) 
                                    VALUES (?,?,?,?,?)""",\
                                    (proteinID,domain.get_HMM(),domain.get_start(),\
                                    domain.get_end(),domain.get_score())\
                                    )




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

def fetch_genomes_dict(database):
    """
    10.11.22
        Args:
           database     Name of the database to be worked on
        Return:
           dictionary genomeID => proteinID list
    """
    genome_dict = {}
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""SELECT proteinID,genomeID FROM Proteins;""")
        
        for count,c in enumerate(cur):
            print(f"\tSelected entries {count}",end="\r")
            if c[1] in genome_dict.keys():
                genome_dict[c[1]].append(c[0])
            else:
                genome_dict[c[1]] = [c[0]]
        
    return genome_dict

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


def fetch_protein_dict(database,genomeID,protein_dict={}):
    """
    1.10.22
        Args:
           database     Name of the database to be worked on
           genomeID     List of genomeIDs to be retrieved 
        Return:
            protein dictionary with key:proteinID => value:protein object for a single genome
    22.02.23
    	removed the Distinct and order by contig,start,end statement in the query because it is
    	possibly unnecessary        
    """
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT Proteins.proteinID, Proteins.genomeID, Proteins.contig, Domains.domStart, Domains.domEnd, Proteins.strand, Proteins.sequence, Domains.domain, Domains.score from Proteins JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE genomeID = ? ;""",[genomeID])
        for row in cur:
            #print(row)
            # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => locustag, 4 => contig,
            # 5 => start, 6 => end, 7 => strand, 8 => domain_count, 9 => sequence,
            #10 => id, 11 => proteinID, 12 => HMM, 13 => dom_start, 14 => dom_end, 15 => score
            proteinID, genomeID, contig, start, end, strand, sequence, domain, score = row
            if proteinID in protein_dict:
                protein = protein_dict[proteinID]
                protein.add_domain(domain, int(start), int(end), score)

            else:
                protein = ParseReports.Protein(proteinID, domain, int(start), int(end), score)
                protein.genomeID = genomeID
                protein.gene_contig = contig
                protein.gene_start = int(start)
                protein.gene_end = int(end)
                protein.gene_strand = strand
                protein.protein_sequence = sequence
                protein_dict[proteinID] = protein
    con.commit()
    con.close()    
    return protein_dict

def fetch_cluster_dict(database, genomeID):
    """
    Fetches cluster data from the database for a specific genomeID.
    
    Args:
        database (str): Path to the SQLite database file.
        genomeID (str): The genomeID to retrieve data for.
    
    Returns:
        dict: A dictionary with clusterID as keys and Cluster objects as values.
    """
    cluster_dict = {}  # clusterID => cluster obj
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        # Query to fetch cluster and protein data for the given genomeID
        query = """
        SELECT DISTINCT Proteins.clusterID, Proteins.proteinID, Domains.domain, Proteins.start, Proteins.end
        FROM Proteins
        JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE Proteins.genomeID = ?;
        """
        cur.execute(query, (genomeID,))
        
        # Loop through the rows returned by the query
        for count, row in enumerate(cur):
            print(f"\tSelected proteins {count}", end="\r")
            clusterID, proteinID, domain, start, end = row
            
            # Get or create the Cluster object for the current clusterID
            cluster = cluster_dict.get(clusterID)
            if cluster is None:
                cluster = Csb_finder.Cluster(clusterID)
                cluster_dict[clusterID] = cluster
            
            # Add the gene (proteinID, domain, start, end) to the cluster
            cluster.add_gene(proteinID, domain, int(start), int(end))
    
    return cluster_dict


def delete_database_clusters_by_genomeID(database, genome_ids_to_delete):
    """
    Deletes clusters from the database based on the genomeID provided in the cluster_batch.

    Args:
        database (str): Path to the SQLite database file.
        cluster_batch (dict): Dictionary where values are cluster objects with a genomeID attribute.
    """
    try:
        # Collect all unique genomeIDs from the cluster_batch
        
        if not genome_ids_to_delete:
            return
        
        with sqlite3.connect(database) as con:
            cur = con.cursor()

            # Start a single transaction for all deletions
            con.execute("BEGIN TRANSACTION")

            # Prepare the SQL statement to delete clusters by genomeID
            placeholders = ', '.join('?' for _ in genome_ids_to_delete)
            sql_delete = f"DELETE FROM Clusters WHERE genomeID IN ({placeholders})"

            # Execute the delete statement with all genomeIDs
            cur.execute(sql_delete, list(genome_ids_to_delete))

            # Commit the transaction after all deletions are done
            con.commit()

    except sqlite3.Error as e:
        print(f"An error occurred while deleting clusters: {e}")



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
    
    


def delete_database_all_keywords(database):
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("DELETE FROM Keywords WHERE keyword LIKE 'csb*%';")
        con.commit()

    return


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
            query = """INSERT INTO Keywords (clusterID, keyword) VALUES (?, ?)"""
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










def fetch_genomeIDs(database):
    """
    Fetches the distinct genomeIDs from the Genomes table in the SQLite database.
    
    Args:
        database: Path to the SQLite database.
    
    Returns:
        A set of distinct genomeIDs.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT genomeID FROM Genomes")
        
        # Use set comprehension to create the set of genomeIDs
        genomeIDs = {row[0] for row in cur.fetchall()}
    
    return genomeIDs

def delete_database_genomeID(database,genomeID):
    """
    2.10.22
        Args: genomeID to be deleted.
        
        Will delete a complete genome with all constrains on protein and cluster table
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" Delete FROM Genomes WHERE genomeID = ?; """,[genomeID])
        con.commit()
    
    return
    
def delete_database_proteins_genomeID(database,genomeID):
    """
    2.10.22
        Args: genomeID to be deleted.
        
        Will delete a complete genome with all constrains on protein and cluster table
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" Delete FROM Proteins WHERE genomeID = ?; """,[genomeID])
        con.commit()
    
    return  
      
def delete_database_proteinID(database,proteinID):
    """
    2.10.22
        Args: proteinID to be deleted.
        
        Will delete a protein with all constrains on protein and domain table
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" Delete FROM Proteins WHERE proteinID = ?; """,[proteinID])
        con.commit()
    
    return
            
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




#fetch_genome_statistic("results/Database2")
#fetch_bulk_data("results/Database2","Phylum","Proteobacteria","Hdr","hdr2")
#fetch_database_all("/home/tomohisa/BioprojectGTDB/Archaea_gtdb/database.db" )



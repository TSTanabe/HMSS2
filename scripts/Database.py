#!/usr/bin/python
#		Module PrepareGenomicData
#		Subroutine FaaToGff
#		Subroutine Translation

import sqlite3
import re
from . import myUtil
from . import ParseReports
from . import Csb_finder #delete after debugging

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
        Strain      varchar(128)                    DEFAULT NULL,
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
    For finished database index the IDs to accelerate the search and join functions. Otherwise
    big databases > 3000 genomes get very slow to impossible slow
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        try:
            cur.execute("""PRAGMA foreign_keys = ON;""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_dom_pid_index ON Domains(proteinID)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_prot_gid_index ON Proteins(genomeID)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_prot_cid_index ON Proteins(clusterID)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_key_cid_index ON Keywords(clusterID)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_clus_gid_index ON Clusters(genomeID)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_dom_did_index ON Domains(domain)""")
            cur.execute("""CREATE INDEX IF NOT EXISTS tab_key_kid_index ON Keywords(keyword)""")
        except sqlite3.Error as e:
            print("An error occurred:", e.args[0])

    
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

def insert_database_protein(database,genomeID,protein_dict):
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
            proteinID = protein.get_proteinID()
            proteinID = f"{genomeID}-{proteinID}"   #genomeID added to proteinID to deal with the multispecies information
            domains = protein.get_domains_dict()
            
            
            #print(proteinID, genomeID, protein.get_gene_locustag(), protein.get_gene_contig(), protein.get_gene_start(), protein.get_gene_end(), protein.get_gene_strand(), protein.get_domain_count(), protein.get_sequence())    
            cur.execute('''INSERT OR IGNORE INTO Proteins
            (proteinID,genomeID,locustag,contig,start,end,strand,dom_count,sequence)
            VALUES (?,?,?,?,?,?,?,?,?) ''',\
            (proteinID, genomeID, protein.get_gene_locustag(), protein.get_gene_contig(),\
             protein.get_gene_start(), protein.get_gene_end(), protein.get_gene_strand(),\
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
                (proteinID, genomeID, protein.get_gene_locustag(), protein.get_gene_contig(),\
                 protein.get_gene_start(), protein.get_gene_end(), protein.get_gene_strand(),\
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



def insert_database_cluster(database,genomeID,cluster_dict):
    """
    1.10.22
        Args:
          
    """
    #TODO behavior -> bei rename sollten alte nicht gelöscht werden, nur doppelte
    #bei keyword select einbauen. wenn gefunden dann nichts tun, sonst reinschreiben.
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT Keywords.id from Keywords LEFT JOIN Clusters ON Keywords.clusterID = Clusters.clusterID WHERE genomeID = ?""", [genomeID])
        old_ids = cur.fetchall()
        
        for cluster_index in cluster_dict:
            cluster = cluster_dict.get(cluster_index)
            clusterID = cluster.get_clusterID()
            cluster_proteins = cluster.get_genes()
            keywords = cluster.get_keywords()
            if clusterID == None:   
                continue
            
            #print(clusterID)
            #cur.execute(""" SELECT clusterID,genomeID from Clusters WHERE genomeID = ? """, [genomeID])
            #print(cur.fetchall())
            cur.execute(""" INSERT OR IGNORE INTO Clusters (clusterID,genomeID) VALUES (?,?) """, \
            (clusterID,genomeID))
            for proteinID in cluster_proteins:
                proteinID = f"{genomeID}-{proteinID}"
                cur.execute(""" UPDATE Proteins SET clusterID = ? WHERE proteinID = ?;""", [clusterID,proteinID])
            
            
            for keyword in keywords:
                #print("BEFORE")
                cur.execute(""" SELECT id,keyword from Keywords LEFT JOIN Clusters ON Keywords.clusterID = Clusters.clusterID WHERE genomeID = ?""", [genomeID])
                #print(cur.fetchall())
                #print("INSERT")
                #print([clusterID,keyword.get_keyword(),keyword.get_completeness(),keyword.get_csb()])
                cur.execute(""" INSERT OR REPLACE INTO Keywords (clusterID,keyword,completeness,collinearity) VALUES (?,?,?,?) """, \
                [clusterID,keyword.get_keyword(),keyword.get_completeness(),keyword.get_csb()])
                #print("AFTER")
                #cur.execute(""" SELECT id,keyword from Keywords LEFT JOIN Clusters ON Keywords.clusterID = Clusters.clusterID WHERE genomeID = ?""", [genomeID])
                #print(cur.fetchall())
            for index in old_ids:
                cur.execute(""" DELETE FROM Keywords where id = ?""",index)
                
                
    con.commit()
    con.close()
    return

##############################################################
########## Alter information from database routines ##########
##############################################################

def update_domain(database,protein_diction,old_tag,new_tag):
    """
    18.11.22
        Args:
           database     Name of the database to be worked on
    """
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        query = """UPDATE Domains SET domain = ? WHERE proteinID = ? and domain = ? """
        
        for proteinID in protein_diction:
            cur.execute(query,[new_tag,proteinID,old_tag])
        
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
    

def fetch_protein_dict(database,genomeID):
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
    protein_dict = {}
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT * from Proteins JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE genomeID = ? ;""",[genomeID])
        
        for row in cur:
            #print(row)
            # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => locustag, 4 => contig,
            # 5 => start, 6 => end, 7 => strand, 8 => domain_count, 9 => sequence,
            #10 => id, 11 => proteinID, 12 => HMM, 13 => dom_start, 14 => dom_end, 15 => score
            if row[0] in protein_dict:
                protein = protein_dict[row[0]]
                protein.add_domain(row[12],row[13],row[14],row[15])

            else:
                protein = ParseReports.Protein(row[0],row[12],row[13],row[14],row[15])
                protein.set_gene_contig(row[4])
                protein.set_gene_start(row[5])
                protein.set_gene_end(row[6])
                protein.set_gene_strand(row[7])
                protein.set_protein_sequence(row[9])
                protein_dict[row[0]] = protein

    con.commit()
    con.close()    
    return protein_dict

def fetch_cluster_dict(database,genomeID):
    """
    1.10.22
        Args:
           database     Name of the database to be worked on
           genomeID     List of genomeIDs to be retrieved 
        Return:
            protein dictionary with key:proteinID => value:protein object for a single genome
    """
    cluster_dict = {} #clusterID => cluster obj
    proteinIDs = {}   #proteinID => clusterID
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        #query = "SELECT clusterID,proteinID from Proteins WHERE proteinID = ?"
        cur.execute("""SELECT DISTINCT Proteins.clusterID,Proteins.proteinID,domain from Proteins JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE genomeID = ? ;""",[genomeID])
        for count,row in enumerate(cur):
            print(f"\tSelected proteins {count}",end="\r")
            if not row[0] in cluster_dict:
                cluster = Csb_finder.Cluster(row[0])
                cluster.add_gene(row[1],row[2])
                cluster_dict[row[0]] = cluster
            else:
                cluster = cluster_dict[row[0]]
                cluster.add_gene(row[1],row[2])

    con.close()    
    return cluster_dict

def fetch_genome_statistic(database):
    """
    18.10.22
        Args:
           database     Name of the database to be worked on
           filepath     Path to directory, will be extended by taxons
        Return:
        
        This one should return all Taxonomic groupings found in the database with number of associated genomes
    """
    
    taxons = ["Superkingdom","Clade","Phylum","Class","Ordnung","Family","Genus"] #Species left out because better to make own file for it
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
    
    
def fetch_bulk_statistic(database):
    """
    18.10.22
        Args:
           database     Name of the database to be worked on
           filepath     Path to directory, will be extended by taxons
        Return:
        
        This one should return all Taxonomic groupings found in the database with number of associated genomes
    """    
    
    con = sqlite3.connect(database)
    cur=con.cursor()
    query = "SELECT DISTINCT Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes ORDER BY Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species"
    cur.execute(query)
    print("\n\nSuperkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species present in the database")
    for row in cur:
        print(row)
    
    
    query = "SELECT domain,count(Proteins.proteinID) FROM Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID GROUP BY domain ORDER BY domain"
    cur.execute(query)
    print("\n\nDomains present in the database")
    for row in cur:
        print(row)
        
    query = "SELECT keyword,count(clusterID) FROM Keywords GROUP BY keyword ORDER BY keyword"
    cur.execute(query)
    print("\n\nKeywords present in the database")
    for row in cur:
        print(row)
        
        
    con.close()

def fetch_genomeID_without_phylogeny(database):
    """
    14.10.22
        Args:
           database     Name of the database to be worked on
        Return:
            set of genomeIDs without phylogeny
    """
    
    genomeIDs = set()
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT DISTINCT genomeID from Genomes WHERE Phylum = "NULL" OR Class = "NULL" OR Family = "NULL" OR Ordnung = "NULL" OR Genus = "NULL" OR Species = "NULL" """)
        
        for row in cur:
            genomeIDs.add(row[0])
    return genomeIDs

def fetch_genomeIDs(database):
    #18.10.22
    genomeIDs = set()
    with sqlite3.connect(database) as con:
        cur=con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute("""SELECT DISTINCT genomeID from Genomes """)
        
        for row in cur:
            genomeIDs.add(row[0])
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



#!/usr/bin/env python3
import sqlite3
from pathlib import Path
import argparse

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS Genomes (
    genomeID varchar(32) PRIMARY KEY NOT NULL,
    Superkingdom varchar(128) DEFAULT 'NULL',
    Clade varchar(128) DEFAULT 'NULL',
    Phylum varchar(128) DEFAULT 'NULL',
    Class varchar(128) DEFAULT 'NULL',
    Ordnung varchar(128) DEFAULT 'NULL',
    Family varchar(128) DEFAULT 'NULL',
    Genus varchar(128) DEFAULT 'NULL',
    Species varchar(128) DEFAULT 'NULL',
    Strain varchar(128) DEFAULT 'NULL',
    TypeStrain tinyint(4) DEFAULT NULL,
    Completeness decimal(5,2) DEFAULT NULL,
    Contamination decimal(5,2) DEFAULT NULL,
    dRep tinyint(1) DEFAULT NULL,
    NCBITaxon int(11) DEFAULT NULL,
    NCBIProject int(11) DEFAULT NULL,
    NCBIBioproject varchar(32) DEFAULT NULL,
    NCBIBiosample varchar(32) DEFAULT NULL,
    NCBIAssembly varchar(32) DEFAULT NULL
);

CREATE TABLE IF NOT EXISTS Clusters (
    clusterID varchar(32) PRIMARY KEY NOT NULL,
    genomeID varchar(32) NOT NULL,
    FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE IF NOT EXISTS Keywords (
    ID integer PRIMARY KEY AUTOINCREMENT,
    clusterID varchar(32) NOT NULL,
    keyword varchar(32) NOT NULL,
    completeness varchar(32) DEFAULT NULL,
    collinearity varchar(32) DEFAULT NULL,
    FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE IF NOT EXISTS Proteins (
    proteinID varchar(128) PRIMARY KEY NOT NULL,
    genomeID varchar(64) NOT NULL,
    clusterID varchar(64) DEFAULT NULL,
    locustag varchar(32) DEFAULT NULL,
    contig varchar(32) DEFAULT NULL,
    start int(11) DEFAULT NULL,
    end int(11) DEFAULT NULL,
    strand varchar(1) DEFAULT NULL,
    comment varchar(12) DEFAULT NULL,
    alternative_hit varchar(64) DEFAULT NULL,
    dom_count smallint(6) DEFAULT NULL,
    valid_hit tinyint(1) DEFAULT 0,
    sequence varchar(4096) DEFAULT NULL,
    UNIQUE(proteinID, genomeID),
    FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE,
    FOREIGN KEY (clusterID) REFERENCES Clusters(clusterID) ON DELETE SET NULL ON UPDATE CASCADE
);

CREATE TABLE IF NOT EXISTS Domains (
    ID integer PRIMARY KEY AUTOINCREMENT,
    proteinID varchar(32) NOT NULL,
    domain varchar(32) DEFAULT NULL,
    score smallint(6) DEFAULT NULL,
    blast_score_ratio smallint(6) DEFAULT NULL,
    identity smallint(6) DEFAULT NULL,
    domStart int(11) DEFAULT NULL,
    domEnd int(11) DEFAULT NULL,
    FOREIGN KEY (proteinID) REFERENCES Proteins(proteinID) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE IF NOT EXISTS Metagenomes (
    metagenomeID varchar(128) PRIMARY KEY NOT NULL,
    genomeID varchar(32) NOT NULL,
    forward_reads INTEGER DEFAULT NULL,
    reverse_reads INTEGER DEFAULT NULL,
    prokaryotic_fraction REAL DEFAULT 1.0,
    FOREIGN KEY (genomeID) REFERENCES Genomes(genomeID) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE IF NOT EXISTS GpkgLengths (
    domain_type VARCHAR(64) PRIMARY KEY NOT NULL,
    protein_length INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS Lineage (
    lineageID VARCHAR(64) PRIMARY KEY NOT NULL,
    root VARCHAR(64) DEFAULT NULL,
    kingdom VARCHAR(256) DEFAULT NULL,
    phylum VARCHAR(256) DEFAULT NULL,
    class VARCHAR(256) DEFAULT NULL,
    "order" VARCHAR(256) DEFAULT NULL,
    family VARCHAR(256) DEFAULT NULL,
    genus VARCHAR(256) DEFAULT NULL,
    species VARCHAR(256) DEFAULT NULL,
    raw_lineage TEXT DEFAULT NULL
);

CREATE TABLE IF NOT EXISTS Placement (
    domain_type VARCHAR(64) NOT NULL,
    readID VARCHAR(256) NOT NULL,
    metagenomeID VARCHAR(128) NOT NULL,
    proteinID VARCHAR(128) DEFAULT NULL,
    lineageID VARCHAR(64) DEFAULT NULL,
    dom_start INTEGER DEFAULT NULL,
    dom_end INTEGER DEFAULT NULL,
    coverage REAL DEFAULT NULL,
    sequence TEXT DEFAULT NULL,
    alignment TEXT DEFAULT NULL,
    PRIMARY KEY (metagenomeID, domain_type, readID),
    FOREIGN KEY (proteinID) REFERENCES Proteins(proteinID) ON DELETE SET NULL ON UPDATE CASCADE,
    FOREIGN KEY (lineageID) REFERENCES Lineage(lineageID) ON DELETE SET NULL ON UPDATE CASCADE,
    FOREIGN KEY (metagenomeID) REFERENCES Metagenomes(metagenomeID) ON DELETE CASCADE ON UPDATE CASCADE
);
"""


def init_db(con: sqlite3.Connection) -> None:
    con.executescript(SCHEMA_SQL)
    con.commit()


def merge_one_db(src_path: Path, out_con: sqlite3.Connection) -> None:
    print(f"Merging: {src_path}")
    src_con = sqlite3.connect(src_path)
    src_con.row_factory = sqlite3.Row

    try:
        out_con.execute("PRAGMA foreign_keys = OFF;")
        src_con.execute("PRAGMA foreign_keys = OFF;")

        # Genomes
        for r in src_con.execute("SELECT * FROM Genomes"):
            out_con.execute(
                """
                INSERT INTO Genomes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(genomeID) DO UPDATE SET
                    Superkingdom   = COALESCE(Genomes.Superkingdom, excluded.Superkingdom),
                    Clade          = COALESCE(Genomes.Clade, excluded.Clade),
                    Phylum         = COALESCE(Genomes.Phylum, excluded.Phylum),
                    Class          = COALESCE(Genomes.Class, excluded.Class),
                    Ordnung        = COALESCE(Genomes.Ordnung, excluded.Ordnung),
                    Family         = COALESCE(Genomes.Family, excluded.Family),
                    Genus          = COALESCE(Genomes.Genus, excluded.Genus),
                    Species        = COALESCE(Genomes.Species, excluded.Species),
                    Strain         = COALESCE(Genomes.Strain, excluded.Strain),
                    TypeStrain     = COALESCE(Genomes.TypeStrain, excluded.TypeStrain),
                    Completeness   = COALESCE(Genomes.Completeness, excluded.Completeness),
                    Contamination  = COALESCE(Genomes.Contamination, excluded.Contamination),
                    dRep           = COALESCE(Genomes.dRep, excluded.dRep),
                    NCBITaxon      = COALESCE(Genomes.NCBITaxon, excluded.NCBITaxon),
                    NCBIProject    = COALESCE(Genomes.NCBIProject, excluded.NCBIProject),
                    NCBIBioproject = COALESCE(Genomes.NCBIBioproject, excluded.NCBIBioproject),
                    NCBIBiosample  = COALESCE(Genomes.NCBIBiosample, excluded.NCBIBiosample),
                    NCBIAssembly   = COALESCE(Genomes.NCBIAssembly, excluded.NCBIAssembly)
                """,
                tuple(r),
            )

        # Clusters
        for r in src_con.execute("SELECT * FROM Clusters"):
            out_con.execute(
                """
                INSERT INTO Clusters(clusterID, genomeID) VALUES (?, ?)
                ON CONFLICT(clusterID) DO NOTHING
                """,
                tuple(r),
            )

        # Proteins
        for r in src_con.execute("SELECT * FROM Proteins"):
            out_con.execute(
                """
                INSERT INTO Proteins
                (proteinID, genomeID, clusterID, locustag, contig, start, end, strand,
                 comment, alternative_hit, dom_count, valid_hit, sequence)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(proteinID) DO UPDATE SET
                    genomeID        = COALESCE(Proteins.genomeID, excluded.genomeID),
                    clusterID       = COALESCE(Proteins.clusterID, excluded.clusterID),
                    locustag        = COALESCE(Proteins.locustag, excluded.locustag),
                    contig          = COALESCE(Proteins.contig, excluded.contig),
                    start           = COALESCE(Proteins.start, excluded.start),
                    end             = COALESCE(Proteins.end, excluded.end),
                    strand          = COALESCE(Proteins.strand, excluded.strand),
                    comment         = COALESCE(Proteins.comment, excluded.comment),
                    alternative_hit = COALESCE(Proteins.alternative_hit, excluded.alternative_hit),
                    dom_count       = COALESCE(Proteins.dom_count, excluded.dom_count),
                    valid_hit       = COALESCE(Proteins.valid_hit, excluded.valid_hit),
                    sequence        = COALESCE(Proteins.sequence, excluded.sequence)
                """,
                tuple(r),
            )

        # Metagenomes
        for r in src_con.execute("SELECT * FROM Metagenomes"):
            out_con.execute(
                """
                INSERT INTO Metagenomes
                (metagenomeID, genomeID, forward_reads, reverse_reads, prokaryotic_fraction)
                VALUES (?,?,?,?,?)
                ON CONFLICT(metagenomeID) DO UPDATE SET
                    genomeID             = COALESCE(Metagenomes.genomeID, excluded.genomeID),
                    forward_reads        = COALESCE(Metagenomes.forward_reads, excluded.forward_reads),
                    reverse_reads        = COALESCE(Metagenomes.reverse_reads, excluded.reverse_reads),
                    prokaryotic_fraction = COALESCE(Metagenomes.prokaryotic_fraction, excluded.prokaryotic_fraction)
                """,
                tuple(r),
            )

        # GpkgLengths
        for r in src_con.execute("SELECT * FROM GpkgLengths"):
            out_con.execute(
                """
                INSERT INTO GpkgLengths(domain_type, protein_length) VALUES (?, ?)
                ON CONFLICT(domain_type) DO UPDATE SET
                    protein_length = excluded.protein_length
                """,
                tuple(r),
            )

        # Lineage
        for r in src_con.execute(
                'SELECT lineageID, root, kingdom, phylum, class, "order", family, genus, species, raw_lineage FROM Lineage'):
            out_con.execute(
                """
                INSERT INTO Lineage
                (lineageID, root, kingdom, phylum, class, "order", family, genus, species, raw_lineage)
                VALUES (?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(lineageID) DO UPDATE SET
                    root       = COALESCE(Lineage.root, excluded.root),
                    kingdom    = COALESCE(Lineage.kingdom, excluded.kingdom),
                    phylum     = COALESCE(Lineage.phylum, excluded.phylum),
                    class      = COALESCE(Lineage.class, excluded.class),
                    "order"    = COALESCE(Lineage."order", excluded."order"),
                    family     = COALESCE(Lineage.family, excluded.family),
                    genus      = COALESCE(Lineage.genus, excluded.genus),
                    species    = COALESCE(Lineage.species, excluded.species),
                    raw_lineage= COALESCE(Lineage.raw_lineage, excluded.raw_lineage)
                """,
                tuple(r),
            )

        # Keywords: ohne ID, mit Deduplikation
        for r in src_con.execute("SELECT clusterID, keyword, completeness, collinearity FROM Keywords"):
            out_con.execute(
                """
                INSERT INTO Keywords (clusterID, keyword, completeness, collinearity)
                SELECT ?, ?, ?, ?
                WHERE NOT EXISTS (
                    SELECT 1 FROM Keywords
                    WHERE clusterID = ?
                      AND keyword = ?
                      AND IFNULL(completeness, '') = IFNULL(?, '')
                      AND IFNULL(collinearity, '') = IFNULL(?, '')
                )
                """,
                (*tuple(r), r["clusterID"], r["keyword"], r["completeness"], r["collinearity"]),
            )

        # Domains: ohne ID, mit Deduplikation
        for r in src_con.execute(
                "SELECT proteinID, domain, score, blast_score_ratio, identity, domStart, domEnd FROM Domains"):
            out_con.execute(
                """
                INSERT INTO Domains
                (proteinID, domain, score, blast_score_ratio, identity, domStart, domEnd)
                SELECT ?, ?, ?, ?, ?, ?, ?
                WHERE NOT EXISTS (
                    SELECT 1 FROM Domains
                    WHERE proteinID = ?
                      AND IFNULL(domain, '') = IFNULL(?, '')
                      AND IFNULL(score, -999999) = IFNULL(?, -999999)
                      AND IFNULL(blast_score_ratio, -999999) = IFNULL(?, -999999)
                      AND IFNULL(identity, -999999) = IFNULL(?, -999999)
                      AND IFNULL(domStart, -999999) = IFNULL(?, -999999)
                      AND IFNULL(domEnd, -999999) = IFNULL(?, -999999)
                )
                """,
                (
                    r["proteinID"], r["domain"], r["score"], r["blast_score_ratio"],
                    r["identity"], r["domStart"], r["domEnd"],
                    r["proteinID"], r["domain"], r["score"], r["blast_score_ratio"],
                    r["identity"], r["domStart"], r["domEnd"],
                ),
            )

        # Placement
        for r in src_con.execute(
                "SELECT domain_type, readID, metagenomeID, proteinID, lineageID, dom_start, dom_end, coverage, sequence, alignment FROM Placement"
        ):
            out_con.execute(
                """
                INSERT INTO Placement
                (domain_type, readID, metagenomeID, proteinID, lineageID, dom_start, dom_end, coverage, sequence, alignment)
                VALUES (?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(metagenomeID, domain_type, readID) DO UPDATE SET
                    proteinID = COALESCE(Placement.proteinID, excluded.proteinID),
                    lineageID = COALESCE(Placement.lineageID, excluded.lineageID),
                    dom_start = COALESCE(Placement.dom_start, excluded.dom_start),
                    dom_end   = COALESCE(Placement.dom_end, excluded.dom_end),
                    coverage  = COALESCE(Placement.coverage, excluded.coverage),
                    sequence  = COALESCE(Placement.sequence, excluded.sequence),
                    alignment = COALESCE(Placement.alignment, excluded.alignment)
                """,
                tuple(r),
            )

        out_con.commit()

    finally:
        src_con.close()
        out_con.execute("PRAGMA foreign_keys = ON;")


def find_databases(root: Path):
    return sorted(p for p in root.rglob("database.db") if p.is_file())


def main():
    parser = argparse.ArgumentParser(
        description="Merge all recursive database.db SQLite files into one output database.")
    parser.add_argument("input_dir", help="Root directory to search recursively")
    parser.add_argument("-o", "--output", default="merged_database.db", help="Output SQLite database")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_db = Path(args.output).resolve()

    dbs = [p for p in find_databases(input_dir) if p != output_db]

    if not dbs:
        raise SystemExit(f"No database.db files found under: {input_dir}")

    print(f"Found {len(dbs)} database(s).")
    out_con = sqlite3.connect(output_db)

    try:
        init_db(out_con)
        for db in dbs:
            merge_one_db(db, out_con)
        out_con.commit()
    finally:
        out_con.close()

    print(f"Done. Merged database written to: {output_db}")


if __name__ == "__main__":
    main()

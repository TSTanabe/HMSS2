#!/usr/bin/python
import os
import sys
import csv
import sqlite3
from typing import Any, Dict, List, Optional, Set, Tuple, Union, Iterable

from Bio import SeqIO

from . import myUtil
from . import ParseReports
from . import Csb_finder

logger = myUtil.logger

#########################################################################
####################### MAIN OUTPUT ROUTINE #############################
#########################################################################


def fetch_fasta_and_hit_data(
    options: Any,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Main output routine: creates hit tables and writes protein FASTA files.

    Args:
        directory (str): Output directory path.
        options (object): Configuration object with output and filter parameters.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
            - protein_dict: Maps proteinID to proteinObj.
            - cluster_dict: Maps clusterID to clusterObj.
            - taxon_dict: Maps genomeID to taxonomy string.

    Example:
        >>> prot, clust, tax = print_fasta_and_hit_table("/out/", options)
    """

    # Limit the genomeIDs and fetch the lineage
    taxon_dict = dict()
    if options.limiter:
        taxon_dict = fetch_limiter_data(options)
    # Fetch protein/cluster dict
    if options.fetch_csbs:
        logger.info(f"Collecting gene clusters containing {options.fetch_csbs}")
        csb_listing = find_csbs_with_proteins_db(
            options.database_directory, options.fetch_csbs, None
        )
        logger.debug(csb_listing)
        logger.info(
            f"Found {len(csb_listing)} types of gene clusters containing {options.fetch_csbs}"
        )
        if len(csb_listing) == 0:
            sys.exit()
        options.fetch_keywords.extend(csb_listing)

    logger.info("Collecting hits from local database")
    protein_dict, cluster_dict, taxon_dict = fetch_bulk_data(
        options.database_directory,
        options.fetch_genomes,
        options.fetch_proteins,
        options.fetch_keywords,
        taxon_dict,
        options.min_completeness,
        options.dataset_divide_sign,
    )
    return protein_dict, cluster_dict, taxon_dict


def print_fasta_and_hit_outputs(
    directory: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, str],
) -> None:
    """
    Main output routine: creates hit tables, taxonomy summaries, and protein FASTA files.

    Args:
        directory (str): Output directory path.
        protein_dict: Mapping proteinID -> proteinObj.
        cluster_dict: Mapping clusterID -> clusterObj.
        taxon_dict: Mapping genomeID -> taxonomy string.
        taxon_divider: Separator used within taxonomy strings.
    """

    # Printing results
    logger.info("Finished collecting protein sequences")
    logger.info("Writing fasta formated output files to disk")
    logger.info("Printing genome information output files")

    # Metadata for taxonomy and hits
    hit_report = directory + "1_hit_table.txt"  # individual hit table in tsv file
    taxonomy_report = directory + "1_taxonomy_table.txt"
    taxonomy_summary = directory + "1_hit_taxonomy_counts.txt"

    # Output hit report
    output_genome_report(
        hit_report, protein_dict, cluster_dict, taxon_dict
    )
    
    # Output unique taxonomy report
    output_unique_taxonomy_table(
        taxonomy_report, taxon_dict
    )
    
    # Output taxonomy summary
    output_taxonomy_summary(taxonomy_summary, taxon_dict)

    # Proteine sequences
    logger.info(f"Printing protein fasta files to {directory}")

    fasta_file_directory = directory + "2_Protein"  # protein fasta files
    files = output_distinct_fasta_reports(
        fasta_file_directory, protein_dict, cluster_dict
    )  # Output per domain
    singletons(directory, files)  # Output singleton per genome and doublicates

    myUtil.clean_empty_files(directory)


#########################################################################
#########################################################################
#########################################################################


def print_command_line_args(output_file: str) -> None:
    """Writes the command-line arguments to a file.

    Args:
        output_file (str): Path to the file to write arguments to.

    Example:
        >>> print_command_line_args("args.txt")
    """
    try:
        with open(output_file, "w") as f:
            f.write("Command-line arguments passed to the script:\n")
            for index, arg in enumerate(sys.argv):
                f.write(f"Argument {index}: {arg}\n")
    except Exception as e:
        logger.error(f"Failed to write to file {output_file}: {e}")


def print_file_content(file_path: str) -> None:
    """Prints the content of a file along with its path.

    Args:
        file_path (str): Path to the file.

    Example:
        >>> print_file_content("csb_patterns.txt")
    """
    try:
        with open(file_path, "r") as file:
            content = file.read()
        logger.info(f"File Path: {file_path}")
        logger.info("File Content:")
        logger.info(content)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
    except Exception as e:
        logger.error(f"An error occurred while reading the file: {e}")


def get_protein_ids_by_domains(
    database: str, domains: List[str], keywords: Optional[List[str]] = None
) -> Dict[str, Set[str]]:
    """Retrieve proteinIDs for given domains, optionally filtered by keywords.

    Args:
        database (str): Path to the database.
        domains (List[str]): List of domain strings.
        keywords (List[str], optional): Keywords for filtering (OR).

    Returns:
        Dict[str, Set[str]]: Maps domain to set of proteinIDs.

    Example:
        >>> get_protein_ids_by_domains("db.sqlite", ["PF00001"])
    """

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Base query
        query = """
        SELECT DISTINCT Domains.domain, Proteins.proteinID
        FROM Proteins
        JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE Domains.domain IN ({})
        """.format(",".join(["?"] * len(domains)))

        # Parameters list
        params = domains

        # Add keyword filter if provided
        if keywords:
            query += """
            AND Proteins.clusterID IN (
                SELECT clusterID
                FROM Keywords
                WHERE 
            """
            keyword_conditions = " OR ".join(["Keywords.keyword = ?" for _ in keywords])
            query += f"{keyword_conditions})"
            params.extend(keywords)

        cur.execute(query, params)
        result = cur.fetchall()

    # Construct the dictionary
    domain_protein_dict = {domain: set() for domain in domains}
    for row in result:
        dom, prot_id = row
        domain_protein_dict[dom].add(prot_id)

    return domain_protein_dict


def get_protein_ids_by_keywords(
    database: str, keywords: List[str], domain: Optional[str] = None
) -> Dict[str, Set[str]]:
    """Retrieve proteinIDs for a set of keywords, optionally filtered by domain.

    Args:
        database (str): Path to the database.
        keywords (List[str]): List of keyword strings.
        domain (str, optional): Domain string to filter by.

    Returns:
        Dict[str, Set[str]]: Maps domain to set of proteinIDs.

    Example:
        >>> get_protein_ids_by_keywords("db.sqlite", ["motifX"])
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Base query
        query = """
        SELECT DISTINCT Domains.domain, Proteins.proteinID
        FROM Proteins
        JOIN Clusters ON Proteins.clusterID = Clusters.clusterID
        JOIN Keywords ON Clusters.clusterID = Keywords.clusterID
        JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE 
        """

        # Add keyword filter
        keyword_conditions = " OR ".join(["Keywords.keyword = ?" for _ in keywords])
        query += f"{keyword_conditions}"

        # Parameters list
        params = keywords

        # Add domain filter if provided
        if domain:
            query += " AND Domains.domain = ?"
            params.append(domain)

        cur.execute(query, params)
        result = cur.fetchall()

    # Construct the dictionary
    domain_protein_dict = {}
    for row in result:
        dom, prot_id = row
        if dom not in domain_protein_dict:
            domain_protein_dict[dom] = set()
        domain_protein_dict[dom].add(prot_id)

    return domain_protein_dict


# Example usage of routine get_protein_ids_by_keywords:
# database = "path/to/your/database.db"
# keywords = ["keyword1", "keyword2"]
# domain = "example_domain"  # Optional
# domain_protein_dict = get_protein_ids_by_keywords(database, keywords, domain)
# print(domain_protein_dict)


def fetch_protein_details(
    database: str, protein_ids: Set[str]
) -> List[Tuple[Any, ...]]:
    """Fetches detailed information for a set of proteinIDs.

    Args:
        database (str): Path to the database.
        protein_ids (Set[str]): Set of protein IDs.

    Returns:
        List[Tuple]: List of protein details (as tuples).

    Example:
        >>> fetch_protein_details("db.sqlite", {"prot1", "prot2"})
    """

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Construct query to fetch detailed information for given proteinIDs
        query = """
        SELECT 
            Proteins.proteinID, Proteins.genomeID, Proteins.clusterID, 
            Proteins.contig, Proteins.start, Proteins.end, Proteins.strand, 
            Proteins.sequence, Domains.domain, Domains.domStart, 
            Domains.domEnd, Domains.score
        FROM Proteins
        LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
        WHERE Proteins.proteinID IN ({})
        """.format(",".join(["?"] * len(protein_ids)))

        cur.execute(query, list(protein_ids))
        rows = cur.fetchall()

    return rows


def write_detail_to_protein_dict(
    details: List[Tuple[Any, ...]],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Converts SQL protein details to protein and cluster dictionaries.

    Args:
        details (List[Tuple]): Tuples as returned from fetch_protein_details.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any]]: protein_dict, cluster_dict.

    Example:
        >>> prot, clust = write_detail_to_protein_dict(details)
    """

    protein_dict = dict()
    cluster_dict = dict()
    count = 1
    for row in details:
        print(f"\tFetched protein sequences: {count}", end="\r")
        count = count + 1
        # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => contig,
        # 4 => start, 5 => end, 6 => strand, 7 => sequence,
        # 8 => domain, 9 => domStart, 10 => domEnd, 11 => score,
        if row[0] in protein_dict.keys():
            protein = protein_dict[row[0]]
            protein.add_domain(row[8], row[9], row[10], row[11])

        else:
            protein = ParseReports.Protein(row[0], row[8], row[9], row[10], row[11])
            protein.set_genomeID(row[1])
            protein.set_clusterID(row[2])
            protein.set_gene_contig(row[3])
            protein.set_gene_start(row[4])
            protein.set_gene_end(row[5])
            protein.set_gene_strand(row[6])
            protein.set_protein_sequence(row[7])
            protein_dict[row[0]] = protein

        if not row[2] is None and not row[2] in cluster_dict:
            cluster = Csb_finder.Cluster(row[2])
            cluster.add_gene(row[0], row[8])
            cluster_dict[row[2]] = cluster

    return protein_dict, cluster_dict


#############################################################################################
#############################################################################################
#############################################################################################


def output_genome_report(
    output_filepath: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    taxon_dict: Dict[str, str] = {},
    genomeID: str = "",
    writemode: str = "w",
    taxon_divider: str = ".",
) -> None:
    """
    Writes the main genome hit table (TSV).

    Columns:
      genomeID, proteinID, domains, dom_scores, dom_coords, contig, gene_start, gene_end,
      gene_strand, locustag, clusterID, csb,
      Superkingdom, Phylum, Class, Order, Family, Genus, Species
    """
    taxon_cols = [
        "genomeID",
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    header_cols = [
        "genomeID",
        "proteinID",
        "domains",
        "dom_scores",
        "dom_coords",
        "contig",
        "gene_start",
        "gene_end",
        "gene_strand",
        "locustag",
        "selection_comment",
        "clusterID",
        *taxon_cols,
    ]
    header = "\t".join(header_cols)

    # sortiert nach genomeID, contig, start
    proteinID_list = sorted(
        protein_dict,
        key=lambda x: (
            protein_dict[x].genomeID,
            protein_dict[x].gene_contig,
            protein_dict[x].gene_start,
        ),
    )

    with open(output_filepath, writemode) as writer:
        writer.write(header + "\n")

        for proteinID in proteinID_list:
            protein = protein_dict[proteinID]
            # protein.get_protein_list() erwartete Reihenfolge laut Docstring:
            # [proteinID, get_domains(), get_domain_scores(), get_domain_coordinates(),
            #  gene_contig, gene_start, gene_end, gene_strand, gene_locustag]
            pl = protein.get_protein_list()

            # clusterID + csb (ohne keyword/completeness)
            clusterID = protein.clusterID
            csb_val = ""
            if clusterID:
                if clusterID in cluster_dict:
                    # cluster.get_cluster_list(",") liefert i. d. R. [clusterID, keyword, completeness, csb]
                    cl = cluster_dict[clusterID].get_cluster_list(",")
                    # defensiv extrahieren:
                    out_clusterID = cl[0] if len(cl) >= 1 else clusterID
                    csb_val = (
                        cl[3] if len(cl) >= 4 else (cl[-1] if len(cl) >= 2 else "")
                    )
                else:
                    out_clusterID = clusterID
            else:
                out_clusterID = ""

            # Selection comment

            # Taxonomie in 7 Spalten
            gid = protein.genomeID
            taxon_levels = [""] * 8
            if gid in taxon_dict and taxon_dict[gid]:
                parts = [p.strip() for p in str(taxon_dict[gid]).split(taxon_divider)]
                for i in range(min(8, len(parts))):
                    taxon_levels[i] = parts[i]
            row = [
                protein.genomeID,  # genomeID
                pl[0],  # proteinID
                pl[1],  # domains
                pl[2],  # dom_scores
                pl[3],  # dom_coords
                pl[4],  # contig
                pl[5],  # gene_start
                pl[6],  # gene_end
                pl[7],  # gene_strand
                pl[8],  # locustag
                protein.selection_comment,
                out_clusterID,
                *taxon_levels,
            ]
            writer.write("\t".join(map(str, row)) + "\n")
    return


def output_taxonomy_summary(
    output_file: str,
    taxon_dict: Dict[str, str],
    taxon_divider: str = ".",
) -> None:
    """
    Build taxonomy statistics directly from `taxon_dict` (no file I/O for input).
    Parsing mirrors `output_genome_report`: split taxonomy by `taxon_divider`,
    use the last 7 parts as ranks (Superkingdom..Species), pad with "" if fewer.

    Writes:
        - `output_file`: a TSV with per-level counts of distinct genomes per taxon.
          Columns: Taxonomic Level, Taxon, distinct genomes
    """
    tax_cols = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    # level -> {taxon -> set(genomeIDs)}
    taxonomy_levels: Dict[str, Dict[str, Set[str]]] = {level: {} for level in tax_cols}
    # set of unique 7-rank lineages (not used further here, but kept for completeness)
    unique_taxonomies: Set[Tuple[str, str, str, str, str, str, str]] = set()

    # Parse taxonomy directly from the dict
    for genome_id, lineage_str in (taxon_dict or {}).items():
        if not genome_id or not lineage_str:
            continue

        # Split and take the last 7 parts as ranks; if fewer, left-pad with ""
        parts = [p.strip() for p in str(lineage_str).split(taxon_divider)]
        tail7: List[str] = (
            parts[-7:] if len(parts) >= 7 else ([""] * (7 - len(parts)) + parts)
        )
        if len(tail7) != 7:
            tail7 = (tail7 + [""] * 7)[:7]  # harden against odd inputs

        # Update per-level sets and unique lineage set
        for level_name, value in zip(tax_cols, tail7):
            if value:
                taxonomy_levels[level_name].setdefault(value, set()).add(genome_id)
        unique_taxonomies.add(tuple(tail7))

    # Write per-level statistics
    with open(output_file, "w", newline="") as out:
        out.write("Taxonomic Level\tTaxon\tdistinct genomes\n")
        for level in tax_cols:
            # Alphabetical within each level; empty strings (if any) go last
            for taxon in sorted(
                taxonomy_levels[level].keys(), key=lambda s: (s == "", s.casefold())
            ):
                out.write(f"{level}\t{taxon}\t{len(taxonomy_levels[level][taxon])}\n")


def output_unique_taxonomy_table(
    output_file: str, taxon_dict: Dict[str, str], taxon_divider: str = "."
) -> str:
    """
    Write a non-redundant taxonomy table.
    - Deduplicates entries based on the 7 taxonomy ranks (Superkingdom..Species).
    - Sorts only by the taxonomy columns, ignoring the GenomeID.
    - Prints GenomeID in the first column.

    Args:
        output_file: Base file name used to generate the output file path.
        taxon_dict: Mapping from genomeID to taxonomy string.
        taxon_divider: Separator used in taxonomy strings (default: ".").

    Returns:
        Path to the written file.
    """
    # Output file name: add suffix "_unique_taxonomies.txt"
    unique_file = output_file.replace(".txt", "") + "_unique_taxonomies.txt"

    # Map from taxonomy tuple (7 ranks) -> representative genomeID
    representative: Dict[Tuple[str, str, str, str, str, str, str], str] = {}

    for gid, lineage_str in (taxon_dict or {}).items():
        if not lineage_str:
            continue
        # Split taxonomy string into parts
        parts = [p.strip() for p in str(lineage_str).split(taxon_divider)]
        # Take the last 7 parts as the taxonomy ranks, fill with "" if fewer
        tail7 = parts[-7:] if len(parts) >= 7 else ([""] * (7 - len(parts)) + parts)
        if len(tail7) != 7:
            tail7 = (tail7 + [""] * 7)[:7]
        tail_tup = tuple(tail7)
        # Keep the first genomeID encountered for each unique taxonomy
        representative.setdefault(tail_tup, gid)

    # Build rows: GenomeID first, then taxonomy ranks
    rows: List[Tuple[str, ...]] = [
        (rep_gid, *taxo) for taxo, rep_gid in representative.items()
    ]

    # Column headers
    header_cols = [
        "GenomeID",
        "Superkingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]

    # Sorting key: use only taxonomy columns (ignore GenomeID in column 0)
    # Empty fields are sorted last; otherwise sort case-insensitively
    def _sort_key(row: Tuple[str, ...]):
        taxo_cols = row[1:]
        return tuple((s == "", (s or "").casefold()) for s in taxo_cols)

    # Sort rows using taxonomy columns
    rows_sorted = sorted(rows, key=_sort_key)

    # Write file
    with open(unique_file, "w", newline="") as uf:
        uf.write("\t".join(header_cols) + "\n")
        for row in rows_sorted:
            uf.write("\t".join(map(str, row)) + "\n")

    return unique_file


##############################################################
########## Fetch information from database routines ##########
##############################################################


def fetch_limiter_data(options: object) -> Dict[str, str]:
    """Fetch limiter data from the database based on specified conditions.

    Args:
        options (object): Configuration object with attributes:
            - database_directory (str): Path to the SQLite database.
            - dataset_limit_lineage (str): Taxonomy field to restrict (e.g., 'Genus').
            - dataset_limit_taxon (str): Taxon name/value to match.
            - dataset_limit_proteins (list of str): Protein domains to limit.
            - dataset_limit_keywords (list of str): Cluster keywords to limit.
            - dataset_divide_sign (str): Separator for taxonomy lineage.

    Returns:
        Dict[str, str]: Maps genomeID to taxonomy lineage string.

    Example:
        >>> tax_dict = fetch_limiter_data(options)
    """
    database: str = options.database_directory
    lineage: str = options.dataset_limit_lineage
    taxon: str = options.dataset_limit_taxon
    proteins: Any = options.dataset_limit_proteins
    keywords: Any = options.dataset_limit_keywords
    trennzeichen: str = options.dataset_divide_sign

    taxon_dict: Dict[str, str] = {}

    query = (
        "SELECT DISTINCT Genomes.genomeID, Superkingdom, Phylum, Class, Ordnung, "
        "Family, Genus, Species FROM Genomes"
    )
    conditions = []
    params = []

    if lineage and taxon:
        conditions.append(f"{lineage} LIKE ?")
        params.append(f"%{taxon}%")
        logger.info(f"Limiting to taxonomy {lineage} like {taxon}")

    if proteins:
        protein_conditions = " OR ".join(["domain LIKE ?"] * len(proteins))
        query += " LEFT JOIN Domains ON Genomes.genomeID = Domains.genomeID"
        conditions.append(f"({protein_conditions})")
        params.extend([f"%{protein}%" for protein in proteins])
        logger.info(f"Limiting to proteins {proteins}")

    if keywords:
        keyword_conditions = " OR ".join(["keyword LIKE ?"] * len(keywords))
        query += " LEFT JOIN Keywords ON Genomes.genomeID = Keywords.genomeID"
        conditions.append(f"({keyword_conditions})")
        params.extend([f"%{keyword}%" for keyword in keywords])
        logger.info(f"Limiting to keywords {keywords}")

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    with sqlite3.connect(database) as con:
        con.execute("PRAGMA foreign_keys = ON;")
        cur = con.cursor()
        cur.execute(query, params)
        for index, row in enumerate(cur):
            logger.debug(f"Selecting genome {index + 1}")
            if row[0] not in taxon_dict:
                taxon_dict[row[0]] = myUtil.taxonomy_lineage(row, trennzeichen)

    return taxon_dict


def parse_protein_set(protein_set_str: str) -> List[str]:
    """Parses a string representation of a protein set into a Python list.

    Args:
        protein_set_str (str): String representation of a protein set,
            e.g. "('prot1','prot2','prot3')"

    Returns:
        List[str]: List of protein names.

    Example:
        >>> parse_protein_set("('protA','protB','protC')")
        ['protA', 'protB', 'protC']
        >>> parse_protein_set("")
        []
    """
    protein_set_str = protein_set_str.strip("()")
    if not protein_set_str:
        return []
    return [protein.strip().strip("'") for protein in protein_set_str.split(",")]


def find_csbs_with_proteins(file_path: str, proteins: List[str]) -> List[str]:
    """Find all CSB identifiers that contain all given proteins.

    Args:
        file_path (str): Path to the CSB file (e.g. './project/Collinear_syntenic_blocks/Csb_output.txt').
        proteins (List[str]): List of protein names to search for.

    Returns:
        List[str]: List of CSB identifiers that contain all given proteins.

    Example:
        >>> find_csbs_with_proteins("Csb_output.txt", ["protA", "protB"])
        ['CSB_01', 'CSB_15']
    """
    logger.debug(f"Csb file filepath is {file_path}")
    csb_list: List[str] = []

    try:
        with open(file_path, "r") as file:
            for line in file:
                if line.strip():
                    parts = line.rstrip("\n").split("\t")
                    csb_id = parts[0]
                    protein_sets = [parse_protein_set(ps) for ps in parts[1:]]
                    # Each 'protein' must occur in at least one of the sets in this CSB
                    if all(
                        any(protein in protein_set for protein_set in protein_sets)
                        for protein in proteins
                    ):
                        csb_list.append(csb_id)
    except Exception as e:
        logger.error(f"Could not read CSB file '{file_path}': {e}")
    return sorted(csb_list)

def find_csbs_with_proteins_db(
    database: str,
    proteins: List[str],
    keyword_prefix: str = "csb-"  # optional: nur CSB-Keywords berücksichtigen
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

    with sqlite3.connect(database) as con:
        con.execute("PRAGMA foreign_keys = ON;")
        cur = con.cursor()
        cur.execute(sql, params)
        return [row[0] for row in cur.fetchall()]
################### Fetch batch results ##############################


def fetch_bulk_data(
    database: str,
    genomes: Optional[List[str]],
    proteins: Optional[List[str]],
    keywords: Optional[List[str]],
    taxon_dict: Optional[Dict[str, str]] = None,
    min_cluster_completeness: float = 0,
    trennzeichen: str = ";",
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
    """
    Fetch bulk data from the database based on specified conditions, using batching
    to avoid SQLite's variable limit.

    This version expects the SELECT produced by `generate_fetch_query(...)` to provide
    stable, unique column aliases. Specifically, the following aliases are used here:

      proteinID, genomeID, clusterID,
      contig, gene_start, gene_end, gene_strand, protein_sequence,
      domain, domStart, domEnd, score,
      dom_count, comment

    Implementation notes:
      - Uses sqlite3.Row for name-based access to row fields (avoids index errors).
      - Keeps your existing flow: build Protein objects on-the-fly, collect Cluster
        stubs (one per clusterID), then enrich clusters with Keywords and genomes
        with taxonomy in batched queries.
      - `min_cluster_completeness` is available for optional filtering after keyword
        hydration (left unchanged here to preserve current behavior).
    """
    protein_dict: Dict[str, Any] = {}
    cluster_dict: Dict[str, Any] = {}
    genomeID_set: Set[str] = set()
    fusion_protIDs: Set[str] = set()
    if taxon_dict is None:
        taxon_dict = {}

    with sqlite3.connect(database) as con:
        # Enable name-based access: row["column_alias"]
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # Pragmas for speed (same as before)
        cur.execute("PRAGMA foreign_keys = ON;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute("PRAGMA synchronous = OFF;")
        cur.execute("PRAGMA temp_store = MEMORY;")

        # 1) Main streaming SELECT: Proteins JOIN Domains (+ Genomes/Keywords) with clear aliases
        query, args = generate_fetch_query(genomes, proteins, keywords, taxon_dict)
        cur.execute(query, args)

        for index, row in enumerate(cur):
            protein_id  = row["proteinID"]
            genome_id   = row["genomeID"]
            cluster_id  = row["clusterID"]

            domain      = row["domain"]
            dom_start   = row["domStart"]
            dom_end     = row["domEnd"]
            score       = row["score"]

            # Progress / debug
            logger.debug(f"Fetched protein for {genome_id}. Total proteins: {index + 1}")

            # Create-or-extend Protein object
            if protein_id in protein_dict:
                protein_dict[protein_id].add_domain(domain, dom_start, dom_end, score)
            else:
                p = ParseReports.Protein(protein_id, domain, dom_start, dom_end, score)
                p.genomeID          = genome_id
                p.clusterID         = cluster_id
                p.gene_contig       = row["contig"]
                p.gene_start        = row["gene_start"]
                p.gene_end          = row["gene_end"]
                p.gene_strand       = row["gene_strand"]
                p.protein_sequence  = row["protein_sequence"]
                p.selection_comment           = row["comment"]
                protein_dict[protein_id] = p
                genomeID_set.add(genome_id)

            # Create cluster stub on first encounter
            if cluster_id is not None and cluster_id not in cluster_dict:
                cl = Csb_finder.Cluster(cluster_id)
                cl.genomeID = genome_id
                cl.add_gene(protein_id, domain)
                cluster_dict[cluster_id] = cl

            # Track proteins that appear to be fused (>=2 domains)
            if row["dom_count"] >= 2:
                fusion_protIDs.add(protein_id)

        logger.info(f"Fetched {len(protein_dict)} proteins.")

        # 2) Add all domains for fused proteins (batched)
        if fusion_protIDs:
            base_query = (
                "SELECT DISTINCT "
                "  proteinID AS proteinID, "
                "  domain    AS domain, "
                "  domStart  AS domStart, "
                "  domEnd    AS domEnd, "
                "  score     AS score "
                "FROM Domains WHERE proteinID IN ({})"
            )
            rows: List[sqlite3.Row] = []
            for batch in batched(list(fusion_protIDs), 500):
                placeholders = ",".join(["?"] * len(batch))
                q = base_query.format(placeholders)
                cur.execute(q, batch)
                rows.extend(cur.fetchall())
            for i, r in enumerate(rows):
                logger.debug(f"Fetched fused domains: {i} proteinID {r['proteinID']}")
                protein_dict[r["proteinID"]].add_domain(
                    r["domain"], r["domStart"], r["domEnd"], r["score"]
                )

        # 3) Hydrate clusters with keywords (batched)
        if cluster_dict:
            base_query = (
                "SELECT DISTINCT "
                "  clusterID    AS clusterID, "
                "  keyword      AS keyword, "
                "  completeness AS completeness, "
                "  collinearity AS collinearity "
                "FROM Keywords WHERE clusterID IN ({})"
            )
            rows: List[sqlite3.Row] = []
            for batch in batched(list(cluster_dict), 500):
                placeholders = ",".join(["?"] * len(batch))
                q = base_query.format(placeholders)
                cur.execute(q, batch)
                rows.extend(cur.fetchall())
            for i, r in enumerate(rows):
                logger.debug(f"Fetched keywords: {i + 1}")
                cid = r["clusterID"]
                cluster_dict[cid].add_keyword(
                    r["keyword"], r["completeness"], r["collinearity"]
                )

            # Optional (kept disabled to preserve current semantics):
            # If you want to drop clusters that don't meet the completeness threshold:
            # if min_cluster_completeness > 0:
            #     cluster_dict = {
            #         cid: cl for cid, cl in cluster_dict.items()
            #         if any(kw.completeness >= min_cluster_completeness for kw in cl.get_keywords())
            #     }

        # 4) Taxonomy info for any missing genomes (batched)
        if not taxon_dict and genomeID_set:
            base_query = (
                "SELECT "
                "  genomeID     AS genomeID, "
                "  Superkingdom AS Superkingdom, "
                "  Phylum       AS Phylum, "
                "  Class        AS Class, "
                "  Ordnung      AS Ordnung, "
                "  Family       AS Family, "
                "  Genus        AS Genus, "
                "  Species      AS Species "
                "FROM Genomes WHERE genomeID IN ({})"
            )
            rows: List[sqlite3.Row] = []
            for batch in batched(list(genomeID_set), 500):
                placeholders = ",".join(["?"] * len(batch))
                q = base_query.format(placeholders)
                cur.execute(q, batch)
                rows.extend(cur.fetchall())
            for i, r in enumerate(rows):
                logger.debug(f"Fetched taxonomy: {i + 1}")
                # myUtil.taxonomy_lineage expects columns in the order shown above;
                # sqlite3.Row supports both index- and name-based access, so passing r is fine.
                taxon_dict[r["genomeID"]] = myUtil.taxonomy_lineage(r, trennzeichen)

    return protein_dict, cluster_dict, taxon_dict


def batched(iterable: Iterable[Any], n: int = 500) -> Iterable[List[Any]]:
    """Yield successive n-sized batches from an iterable.

    Args:
        iterable (Iterable[Any]): Input items to batch.
        n (int): Batch size.

    Yields:
        List[Any]: Next batch of up to n items.

    Example:
        >>> list(batched([1,2,3,4,5], 2))
        [[1, 2], [3, 4], [5]]
    """
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == n:
            yield batch
            batch = []
    if batch:
        yield batch


def generate_fetch_query(
    genomes: Optional[List[str]],
    proteins: Optional[List[str]],
    keywords: Optional[List[str]],
    taxon_dict: Optional[Dict[str, Any]],
) -> Tuple[str, List[Any]]:
    """
    Build a parameterized SQL SELECT that joins Proteins, Domains, Keywords, and Genomes,
    returning a row-per-(proteinID, domain) with stable column aliases suitable for
    name-based access (sqlite3.Row). This function only constructs the SQL and its bound
    parameters; it does not execute the query.

    Filtering semantics:
      - genomes: partial match (LIKE) on Genomes.genomeID
      - proteins: partial match (LIKE) on Domains.domain
      - keywords: exact match (=) on Keywords.keyword
      - taxon_dict: restrict to the set of genomeIDs present as keys (IN (...))

    Notes on JOINs and NULL-handling:
      - We use LEFT JOIN for Domains and Keywords so we can still fetch Proteins rows even if
        a given protein has no recorded domain hits or keyworded cluster (depending on filters).
      - HOWEVER: when *no* explicit proteins filter is specified, we add "d.domain IS NOT NULL"
        to avoid yielding rows where the domain fields would be NULL (downstream code expects
        valid domain values when constructing Protein objects). If you truly want bare proteins
        without domains, remove that condition.

    Returns:
      (sql, args): sql is the SELECT statement with "?" placeholders; args is the list of
      parameters in correct order. Always use the returned args with cursor.execute(sql, args)
      to ensure correctness and protection against SQL injection.

    Performance considerations:
      - The query is DISTINCT to guard against duplicate rows from JOINs.
      - Ensure indexes exist on:
          Proteins(proteinID), Proteins(genomeID), Proteins(clusterID)
          Domains(proteinID), Domains(domain)
          Keywords(clusterID), Keywords(keyword)
          Genomes(genomeID)
        so that LIKE/IN/EXISTS patterns remain efficient.
    """
    # Base SELECT with explicit, unique aliases for every output column.
    # This prevents ambiguity in sqlite3.Row lookups when different tables have identically named columns.
    query = """
        SELECT DISTINCT
            p.proteinID        AS proteinID,
            g.genomeID         AS genomeID,
            p.clusterID        AS clusterID,
            p.contig           AS contig,
            p.start            AS gene_start,
            p.end              AS gene_end,
            p.strand           AS gene_strand,
            p.sequence         AS protein_sequence,
            d.domain           AS domain,
            d.domStart         AS domStart,
            d.domEnd           AS domEnd,
            d.score            AS score,
            p.dom_count        AS dom_count,
            p.comment          AS comment
        FROM Proteins p
        LEFT JOIN Domains  d ON d.proteinID = p.proteinID
        LEFT JOIN Keywords k ON k.clusterID  = p.clusterID
        LEFT JOIN Genomes  g ON g.genomeID   = p.genomeID
    """

    conditions: List[str] = []
    args: List[Any] = []

    # genomes: partial match (LIKE) over Genomes.genomeID
    # Build (g.genomeID LIKE ? OR g.genomeID LIKE ? OR ...)
    if genomes:
        conditions.append("(" + " OR ".join(["g.genomeID LIKE ?"] * len(genomes)) + ")")
        args.extend([f"%{gid}%" for gid in genomes])

    # proteins: partial match (LIKE) over Domains.domain
    # Build (d.domain LIKE ? OR d.domain LIKE ? OR ...)
    if proteins:
        conditions.append("(" + " OR ".join(["d.domain LIKE ?"] * len(proteins)) + ")")
        args.extend([f"%{prot}%" for prot in proteins])

    # keywords: exact match over Keywords.keyword
    # Build (k.keyword = ? OR k.keyword = ? OR ...)
    if keywords:
        conditions.append("(" + " OR ".join(["k.keyword = ?"] * len(keywords)) + ")")
        args.extend(list(keywords))

    # taxon_dict: hard restriction to specific genomeIDs (IN (...))
    # Note: order of keys is preserved as we only append parameters; SQLite doesn't care about order for IN.
    if taxon_dict:
        conditions.append("g.genomeID IN (" + ",".join("?" * len(taxon_dict)) + ")")
        args.extend(list(taxon_dict.keys()))

    # If no protein-domain filter is provided, ensure we only return rows where a domain exists.
    # This keeps downstream logic simple (it expects to call .add_domain(...) with non-NULL values).
    if not proteins:
        conditions.append("d.domain IS NOT NULL")

    # Final WHERE clause assembly (if any conditions were added).
    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    return query, args



########## File Output Routines ##########


def output_distinct_fasta_reports(
    directory: str,
    protein_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    writemode: str = "w",
) -> Set[str]:
    """Writes all protein sequences into distinct FASTA files by domain class and for fusion domains.

    Args:
        directory (str): Output directory, e.g. 'output/'.
        protein_dict (Dict[str, Any]): Maps proteinID to Protein object.
        cluster_dict (Dict[str, Any]): Maps clusterID to Cluster object.
        writemode (str, optional): File mode, "w" for overwrite/new, "a" for append (default: "w").

    Returns:
        Set[str]: Set of all written file paths.

    Example:
        >>> files = output_distinct_fasta_reports("output/", protein_dict, cluster_dict)
        >>> print(files)
        {'output/_PF00001.faa', 'output/_fused_domain_PF00001.faa'}
    """
    HMM_dict: Dict[str, list] = {}
    fusion_dict: Dict[str, Any] = {}
    files: Set[str] = set()

    # Group proteins by their domain type
    for proteinID, protein in protein_dict.items():
        HMM = protein.get_domains()
        HMM_dict.setdefault(HMM, []).append(proteinID)
        if protein.get_domain_count() > 1:  # fusion protein detected
            fusion_dict[proteinID] = protein

    # Output: per domain class
    for HMM, proteinID_list in HMM_dict.items():
        filepath = directory + f"_{HMM}.faa"
        files.add(filepath)
        with open(filepath, writemode) as writer:
            for proteinID in proteinID_list:
                protein = protein_dict[proteinID]
                genomeID = protein.genomeID
                proteinlist = protein.get_protein_list()
                sequence = str(protein.protein_sequence).replace("*", "")

                clusterID = protein.clusterID
                if clusterID in cluster_dict:
                    cluster = cluster_dict[clusterID]
                    clusterlist = cluster.get_cluster_list(",")
                    out = f">{genomeID}-{' '.join(proteinlist[:-5])} {' '.join(clusterlist)}\n"
                    writer.write(out)
                    writer.write(sequence + "\n")
                else:
                    out = f">{genomeID}-{' '.join(proteinlist[:-5])}\n"
                    writer.write(out)
                    writer.write(sequence + "\n")

    # Output: fused domains as separate FASTA
    for protein in fusion_dict.values():
        clusterID = protein.clusterID
        domain_dict = protein.get_domains_dict()
        genomeID = protein.genomeID
        proteinlist = protein.get_protein_list()
        sequence = str(protein.protein_sequence).replace("*", "")
        for domain in domain_dict.values():
            HMM = domain.HMM
            domain_sequence = sequence[domain.start : domain.end]
            filepath = os.path.join(directory, f"_fused_domain_{HMM}.faa")
            files.add(filepath)
            with open(filepath, "a") as writer:
                if clusterID in cluster_dict:
                    cluster = cluster_dict[clusterID]
                    clusterlist = cluster.get_cluster_list(",")
                    out = f">{genomeID}-{' '.join(proteinlist[:-5])} {' '.join(clusterlist)}\n"
                    writer.write(out)
                    writer.write(domain_sequence + "\n")
                else:
                    out = f">{genomeID}-{' '.join(proteinlist[:-5])}\n"
                    writer.write(out)
                    writer.write(domain_sequence + "\n")
    #logger.info(f"FASTA output written to: {files}")
    return files


##############################################################
##########          Fasta File Postprocessing       ##########
##############################################################


def singletons(directory, filepaths):
    """
    01.11.22
        Args:
            directory to store new files
            filepaths with the files to be processed

        Routine should filter all singletons, meaning each new file will result in one protein per type per genome
        paralogs are filtered out and should be wrote separately. Output is a faa fasta file
    """

    for filepath in filepaths:
        record_dict = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        genome_list = []
        for record in record_dict.keys():
            genomeID = record.split("-")[0]
            genome_list.append(genomeID)

        doublicate = set([x for x in genome_list if genome_list.count(x) > 1])

        path = os.path.splitext(filepath)[0]
        name = os.path.basename(path)
        single = open(directory + f"/{name}_singleton.faa", "w")
        double = open(directory + f"/{name}_doublicate.faa", "w")

        for record in record_dict.keys():
            # print(record_dict[record])
            sequence = record_dict[record]
            genomeID = record.split("-")[0]
            if not genomeID in doublicate:
                single.write(">" + sequence.description + "\n")
                single.write(str(sequence.seq) + "\n")
            else:
                double.write(">" + sequence.description + "\n")
                double.write(str(sequence.seq) + "\n")
        single.close()
        double.close()


##############################################################
##########    Fasta File Postprocessing ON COMMAND    ########
##############################################################


def add_taxonomy(database, filepath, trennzeichen="_"):
    """
    12.03.23
    taxonomy lineage and domain type should be added to a fasta file and replace the
    previous headers
    trennzeichen is the character between the separate informations replacing white characters for
    more readeability
    """
    print("Collecting genome identifier")
    record_dict = {}
    taxon_dict = {}
    writer = open(filepath + ".taxon_names", "w")
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        for record in SeqIO.parse(filepath, "fasta"):
            dataset_range_line = ""
            genomeID, proteinID = record.id.split("-", maxsplit=1)
            genomeID = myUtil.getGenomeID(genomeID)
            # check genome has multiple proteins in the fasta
            if genomeID not in record_dict:
                record_dict[genomeID] = 1
            else:
                record_dict[genomeID] = record_dict[genomeID] + 1

            # Fetch taxonomy
            # Creating the first column for the dataset
            cur.execute(
                """SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ?""",
                (genomeID,),
            )
            row = cur.fetchone()
            lineage = myUtil.taxonomy_lineage(row, trennzeichen)
            dataset_range_line = lineage if lineage else record.id
            if record_dict[genomeID] > 1:
                dataset_range_line += "_" + str(record_dict[genomeID])

            try:
                dom_type = record.description.split(" ")[1]
                dom_type = " " + dom_type + " "
            except:
                dom_type = ""
            try:
                # print(f"Insert {genomeID}")
                # print(record_dict[genomeID[:-2]])
                writer.write(">" + dataset_range_line + dom_type + "\n")
                writer.write(str(record.seq) + "\n")
            except:
                print(f"WARNING: No taxonomy found for {genomeID}")

    writer.close()


def add_genomic_context(database, filepath):
    """
    22.02.23
        Args:
            database for the sequence file
            filepath with the file to be processed
        This routine takes a fasta file, iterates through the sequences and
        writes down the specific genomic context in which this was found together with
        some genomic features

        nicht schön aber funktioniert vielleicht so
    """
    cp_dict = {}
    print("Adding genomic context")
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = "SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({})"

        # Group protein IDs into batches of 1000 for efficient querying
        batch_size = 1000
        protein_ids = [
            f"'{record.id.split('-', 1)[-1]}'"
            for record in SeqIO.parse(filepath, "fasta")
        ]

        for i in range(0, len(protein_ids), batch_size):
            batch = ",".join(protein_ids[i : i + batch_size])
            cur.execute(query.format(batch))
            results = cur.fetchall()
            for protein_id, cluster_id in results:
                cp_dict[protein_id] = cluster_id

        cur.execute("""PRAGMA foreign_keys = ON;""")
        with open(filepath + "_gene_vicinity", "a") as writer:
            writer.write(
                "Index\tProteinID\tDomain(s)\tHit_score\thit_align\tcontig\tstart\tend\tstrand\tSuperkingdom\tClade\tPhylum\tClass\tOrdnung\tFamily\tGenus\tSpecies"
            )

        print(f"\tAssigning genetic environment:")
        for current_proteinID, current_clusterID in cp_dict.items():
            print(f"\tFetched: {current_proteinID} {current_clusterID}", end="\r")
            protein_dict = dict()
            fusion_protIDs = dict()
            query = "SELECT DISTINCT Proteins.proteinID,Proteins.genomeID,Proteins.clusterID,contig,start,end,strand,sequence,domain,domStart,domEnd,score,dom_count from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE Proteins.clusterID = ?"
            cur.execute(
                query,
                [
                    current_clusterID,
                ],
            )
            # Fill the protein dict

            for row in cur:
                # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => contig,
                # 4 => start, 5 => end, 6 => strand, 7 => sequence,
                # 8 => domain, 9 => domStart, 10 => domEnd, 11 => score,
                if row[0] in protein_dict.keys():
                    protein = protein_dict[row[0]]
                    protein.add_domain(row[8], row[9], row[10], row[11])
                else:
                    protein = ParseReports.Protein(
                        row[0], row[8], row[9], row[10], row[11]
                    )
                    protein.set_genomeID(row[1])
                    protein.set_clusterID(row[2])
                    protein.set_gene_contig(row[3])
                    protein.set_gene_start(row[4])
                    protein.set_gene_end(row[5])
                    protein.set_gene_strand(row[6])
                    protein_dict[row[0]] = protein

                if row[12] > 1:
                    fusion_protIDs[row[0]] = 1

            # Also collect domains from fusion proteins

            if fusion_protIDs:
                count = 1
                query = """SELECT DISTINCT proteinID,domain,domStart,domEnd,score FROM Domains WHERE proteinID = ? """
                add_array = []
                proteins_array = []
                args = []

                for proteinID in fusion_protIDs.keys():
                    arguments = args.copy()
                    arguments.insert(0, proteinID)
                    cur.execute(query, arguments)
                    for row in cur:
                        print(f"\tFetched fused domains: {count}", end="\r")
                        count = count + 1
                        protein = protein_dict[row[0]]
                        protein.add_domain(row[1], row[2], row[3], row[4])

                print("")
            # an diesem punkt haben wir ein vollständiges unsortiertes protein_dict aber keine txonomy information
            # taxonomie info holen
            taxon_dict = dict()
            query = "SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ? "
            genomeID = current_proteinID.split("-")[0]
            cur.execute(
                query,
                [
                    genomeID,
                ],
            )
            for row in cur:
                taxon_dict[row[0]] = row[1:]

            # create indices for the
            proteinID_list = sorted(
                protein_dict, key=lambda x: protein_dict[x].gene_start
            )
            try:
                index = proteinID_list.index(current_proteinID) * -1
            except:
                index = 0
            with open(filepath + "_gene_vicinity", "a") as writer:
                for proteinID in proteinID_list:
                    # Write each line includes all information about one protein
                    protein = protein_dict[proteinID]
                    proteinlist = (
                        protein.get_protein_list()
                    )  # list representation of a protein object

                    out = (
                        str(index)
                        + "\t"
                        + "\t".join(proteinlist)
                        + "\t".join(taxon_dict[genomeID])
                        + "\n"
                    )
                    index = index + 1
                    writer.write(out)

    # Specify input and output file names
    input_file = filepath + "_gene_vicinity"
    output_file = filepath + "_gene_vicinity_sorted"

    # Specify the headers of the columns to sort onSuperkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species
    sort_headers = [
        "Superkingdom",
        "Phylum",
        "Class",
        "Ordnung",
        "Family",
        "Genus",
        "Species",
    ]

    # Read in the input file and store the rows as a list of dictionaries
    with open(input_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = [row for row in reader]

    # Sort the rows hierarchically based on the sort columns
    sorted_rows = sorted(
        rows,
        key=lambda x: (
            x[sort_headers[0]],
            x[sort_headers[1]],
            x[sort_headers[2]],
            x[sort_headers[3]],
            x[sort_headers[4]],
            x[sort_headers[5]],
        ),
    )

    # Write the sorted rows to the output file
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow(row)


def fetch_database_all(database):
    """
    1.10.22
        Args:
            ONLY FOR DEBUGGING
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        cur.execute(""" SELECT * FROM Proteins; """)
        con.commit()
        # print("----------Proteins--------------")
        # print(cur.fetchall())

        cur = con.cursor()
        cur.execute(""" SELECT * FROM Clusters; """)
        con.commit()
        print("----------Clusters--------------")
        print(cur.fetchall())

        cur = con.cursor()
        cur.execute(""" SELECT * FROM Keywords; """)
        con.commit()
        print("----------Keywords--------------")
        print(cur.fetchall())

        cur = con.cursor()
        cur.execute(""" SELECT * FROM Genomes; """)
        con.commit()
        print("----------Genomes--------------")
        print(cur.fetchall())
    return

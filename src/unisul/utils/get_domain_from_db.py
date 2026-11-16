#!/usr/bin/env python3

import argparse
import os
import sqlite3
import statistics
from collections import Counter, defaultdict

TAXON_RANKS = ["Family", "Genus", "Species"]


def connect_readonly(db_path):
    uri = f"file:{os.path.abspath(db_path)}?mode=ro"
    con = sqlite3.connect(uri, uri=True)
    con.execute("PRAGMA query_only = TRUE;")
    con.execute("PRAGMA temp_store = MEMORY;")
    return con


def normalize_taxon(value):
    if value is None:
        return "NA"

    value = str(value).strip()

    if value == "" or value.upper() in {"NULL", "NONE", "NA", "N/A"}:
        return "NA"

    return value


def summarize_counts(counts):
    if not counts:
        return 0.0, 0.0, 0

    mean_val = sum(counts) / len(counts)
    median_val = statistics.median(counts)
    max_val = max(counts)

    return mean_val, median_val, max_val


def read_domains_from_file(path):
    domains = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue
            if line.startswith("#"):
                continue

            domains.append(line.split()[0])

    return domains


def write_tsv(rows, out_path, fieldnames):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    with open(out_path, "w") as handle:
        handle.write("\t".join(fieldnames) + "\n")

        for row in rows:
            handle.write(
                "\t".join(str(row.get(field, "")) for field in fieldnames) + "\n"
            )


def fetch_paralog_counts_for_rank(con, domains, taxon_rank):
    """
    Returns per-genome paralog counts for one taxonomy rank.

    Output rows:
        domain
        rank
        taxon
        genomeID
        total_paralogs
        valid_paralogs
        invalid_paralogs
    """

    rows_out = []

    for domain in domains:
        print(f"[INFO] Processing domain={domain} rank={taxon_rank}", flush=True)

        query = f"""
            SELECT
                d.domain AS domain,
                COALESCE(g.{taxon_rank}, 'NA') AS taxon,
                p.genomeID AS genomeID,

                COUNT(DISTINCT d.proteinID) AS total_paralogs,

                COUNT(DISTINCT CASE
                    WHEN p.valid_hit = 1 THEN d.proteinID
                END) AS valid_paralogs,

                COUNT(DISTINCT CASE
                    WHEN p.valid_hit = 0 THEN d.proteinID
                END) AS invalid_paralogs

            FROM Domains d

            JOIN Proteins p
              ON d.proteinID = p.proteinID

            LEFT JOIN Genomes g
              ON p.genomeID = g.genomeID

            WHERE d.domain = ?

            GROUP BY
                d.domain,
                g.{taxon_rank},
                p.genomeID

            ORDER BY
                taxon,
                total_paralogs DESC,
                p.genomeID
        """

        rows = con.execute(query, (domain,)).fetchall()

        for (
            domain_name,
            taxon,
            genome_id,
            total_paralogs,
            valid_paralogs,
            invalid_paralogs,
        ) in rows:
            rows_out.append(
                {
                    "domain": domain_name,
                    "rank": taxon_rank,
                    "taxon": normalize_taxon(taxon),
                    "genomeID": genome_id,
                    "total_paralogs": int(total_paralogs),
                    "valid_paralogs": int(valid_paralogs),
                    "invalid_paralogs": int(invalid_paralogs),
                }
            )

    return rows_out


def summarize_taxon_paralogs(genome_rows):
    """
    Creates one summary row per domain x rank x taxon.

    Includes:
        total/valid/invalid hit sums
        number of genomes
        mean/median/max paralogs per positive genome
        number of genomes with >=2 paralogs
    """

    grouped = defaultdict(list)

    for row in genome_rows:
        key = (row["domain"], row["rank"], row["taxon"])
        grouped[key].append(row)

    summary_rows = []

    for (domain, rank, taxon), rows in sorted(grouped.items()):
        total_counts = [r["total_paralogs"] for r in rows]
        valid_counts = [r["valid_paralogs"] for r in rows]
        invalid_counts = [r["invalid_paralogs"] for r in rows]

        mean_total, median_total, max_total = summarize_counts(total_counts)
        mean_valid, median_valid, max_valid = summarize_counts(valid_counts)
        mean_invalid, median_invalid, max_invalid = summarize_counts(invalid_counts)

        summary_rows.append(
            {
                "domain": domain,
                "rank": rank,
                "taxon": taxon,
                "genomes_with_domain": len(rows),
                "total_hits": sum(total_counts),
                "valid_hits": sum(valid_counts),
                "invalid_hits": sum(invalid_counts),
                "genomes_with_2plus_total_paralogs": sum(
                    1 for x in total_counts if x >= 2
                ),
                "genomes_with_2plus_valid_paralogs": sum(
                    1 for x in valid_counts if x >= 2
                ),
                "genomes_with_2plus_invalid_paralogs": sum(
                    1 for x in invalid_counts if x >= 2
                ),
                "mean_total_paralogs_per_positive_genome": mean_total,
                "median_total_paralogs_per_positive_genome": median_total,
                "max_total_paralogs_per_genome": max_total,
                "mean_valid_paralogs_per_positive_genome": mean_valid,
                "median_valid_paralogs_per_positive_genome": median_valid,
                "max_valid_paralogs_per_genome": max_valid,
                "mean_invalid_paralogs_per_positive_genome": mean_invalid,
                "median_invalid_paralogs_per_positive_genome": median_invalid,
                "max_invalid_paralogs_per_genome": max_invalid,
            }
        )

    return summary_rows


def make_paralog_distribution(genome_rows):
    """
    Creates frequency distribution of paralog numbers.

    Output:
        domain
        rank
        taxon
        paralog_type
        paralogs
        genomes
    """

    counter = Counter()

    for row in genome_rows:
        for paralog_type in ["total", "valid", "invalid"]:
            count_col = f"{paralog_type}_paralogs"

            key = (
                row["domain"],
                row["rank"],
                row["taxon"],
                paralog_type,
                row[count_col],
            )

            counter[key] += 1

    distribution_rows = []

    for (domain, rank, taxon, paralog_type, paralogs), genomes in sorted(
        counter.items()
    ):
        distribution_rows.append(
            {
                "domain": domain,
                "rank": rank,
                "taxon": taxon,
                "paralog_type": paralog_type,
                "paralogs": paralogs,
                "genomes": genomes,
            }
        )

    return distribution_rows


def summarize_domain_level(genome_rows):
    """
    Domain-level summary across all taxa and ranks.
    """

    grouped = defaultdict(list)

    for row in genome_rows:
        key = row["domain"]
        grouped[key].append(row)

    summary_rows = []

    for domain, rows in sorted(grouped.items()):
        total_counts = [r["total_paralogs"] for r in rows]
        valid_counts = [r["valid_paralogs"] for r in rows]
        invalid_counts = [r["invalid_paralogs"] for r in rows]

        mean_total, median_total, max_total = summarize_counts(total_counts)
        mean_valid, median_valid, max_valid = summarize_counts(valid_counts)
        mean_invalid, median_invalid, max_invalid = summarize_counts(invalid_counts)

        summary_rows.append(
            {
                "domain": domain,
                "genome_taxon_rows": len(rows),
                "total_hits": sum(total_counts),
                "valid_hits": sum(valid_counts),
                "invalid_hits": sum(invalid_counts),
                "genome_taxon_rows_with_2plus_total_paralogs": sum(
                    1 for x in total_counts if x >= 2
                ),
                "genome_taxon_rows_with_2plus_valid_paralogs": sum(
                    1 for x in valid_counts if x >= 2
                ),
                "genome_taxon_rows_with_2plus_invalid_paralogs": sum(
                    1 for x in invalid_counts if x >= 2
                ),
                "mean_total_paralogs_per_positive_genome_taxon_row": mean_total,
                "median_total_paralogs_per_positive_genome_taxon_row": median_total,
                "max_total_paralogs_per_genome_taxon_row": max_total,
                "mean_valid_paralogs_per_positive_genome_taxon_row": mean_valid,
                "median_valid_paralogs_per_positive_genome_taxon_row": median_valid,
                "max_valid_paralogs_per_genome_taxon_row": max_valid,
                "mean_invalid_paralogs_per_positive_genome_taxon_row": mean_invalid,
                "median_invalid_paralogs_per_positive_genome_taxon_row": median_invalid,
                "max_invalid_paralogs_per_genome_taxon_row": max_invalid,
            }
        )

    return summary_rows


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create paralog count summaries for each taxon from Family to Species. "
            "Outputs per-genome counts, per-taxon mean/median/max summaries, "
            "and paralog-number frequency distributions."
        )
    )

    parser.add_argument(
        "-db",
        "--database",
        required=True,
        help="Path to HAMSTER SQLite database.",
    )

    parser.add_argument(
        "-d",
        "--domains",
        nargs="+",
        default=None,
        help="Domain names to analyse, e.g. Sgdh Sgald Sladh.",
    )

    parser.add_argument(
        "-df",
        "--domain-file",
        default=None,
        help="Optional file with one domain name per line.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--include-na",
        action="store_true",
        help="Keep NA taxonomy values. Default: remove NA taxa from summaries.",
    )

    args = parser.parse_args()

    domains = []

    if args.domains:
        domains.extend(args.domains)

    if args.domain_file:
        domains.extend(read_domains_from_file(args.domain_file))

    domains = list(dict.fromkeys(domains))

    if not domains:
        raise ValueError("No domains provided. Use --domains or --domain-file.")

    print(f"[INFO] Domains: {', '.join(domains)}", flush=True)
    print(f"[INFO] Taxonomic ranks: {', '.join(TAXON_RANKS)}", flush=True)

    con = connect_readonly(args.database)

    try:
        all_genome_rows = []

        for rank in TAXON_RANKS:
            rank_rows = fetch_paralog_counts_for_rank(
                con=con,
                domains=domains,
                taxon_rank=rank,
            )
            all_genome_rows.extend(rank_rows)

    finally:
        con.close()

    if not args.include_na:
        before = len(all_genome_rows)
        all_genome_rows = [row for row in all_genome_rows if row["taxon"] != "NA"]
        print(
            f"[INFO] Removed NA taxonomy rows: {len(all_genome_rows):,}/{before:,} retained",
            flush=True,
        )

    print(f"[INFO] Per-genome taxon rows: {len(all_genome_rows):,}", flush=True)

    taxon_summary_rows = summarize_taxon_paralogs(all_genome_rows)
    distribution_rows = make_paralog_distribution(all_genome_rows)
    domain_summary_rows = summarize_domain_level(all_genome_rows)

    genome_fields = [
        "domain",
        "rank",
        "taxon",
        "genomeID",
        "total_paralogs",
        "valid_paralogs",
        "invalid_paralogs",
    ]

    taxon_summary_fields = [
        "domain",
        "rank",
        "taxon",
        "genomes_with_domain",
        "total_hits",
        "valid_hits",
        "invalid_hits",
        "genomes_with_2plus_total_paralogs",
        "genomes_with_2plus_valid_paralogs",
        "genomes_with_2plus_invalid_paralogs",
        "mean_total_paralogs_per_positive_genome",
        "median_total_paralogs_per_positive_genome",
        "max_total_paralogs_per_genome",
        "mean_valid_paralogs_per_positive_genome",
        "median_valid_paralogs_per_positive_genome",
        "max_valid_paralogs_per_genome",
        "mean_invalid_paralogs_per_positive_genome",
        "median_invalid_paralogs_per_positive_genome",
        "max_invalid_paralogs_per_genome",
    ]

    distribution_fields = [
        "domain",
        "rank",
        "taxon",
        "paralog_type",
        "paralogs",
        "genomes",
    ]

    domain_summary_fields = [
        "domain",
        "genome_taxon_rows",
        "total_hits",
        "valid_hits",
        "invalid_hits",
        "genome_taxon_rows_with_2plus_total_paralogs",
        "genome_taxon_rows_with_2plus_valid_paralogs",
        "genome_taxon_rows_with_2plus_invalid_paralogs",
        "mean_total_paralogs_per_positive_genome_taxon_row",
        "median_total_paralogs_per_positive_genome_taxon_row",
        "max_total_paralogs_per_genome_taxon_row",
        "mean_valid_paralogs_per_positive_genome_taxon_row",
        "median_valid_paralogs_per_positive_genome_taxon_row",
        "max_valid_paralogs_per_genome_taxon_row",
        "mean_invalid_paralogs_per_positive_genome_taxon_row",
        "median_invalid_paralogs_per_positive_genome_taxon_row",
        "max_invalid_paralogs_per_genome_taxon_row",
    ]

    per_genome_out = os.path.join(
        args.output_dir,
        "domain_family_to_species_genome_paralog_counts.tsv",
    )

    taxon_summary_out = os.path.join(
        args.output_dir,
        "domain_family_to_species_taxon_paralog_summary.tsv",
    )

    distribution_out = os.path.join(
        args.output_dir,
        "domain_family_to_species_paralog_distribution.tsv",
    )

    domain_summary_out = os.path.join(
        args.output_dir,
        "domain_family_to_species_domain_summary.tsv",
    )

    write_tsv(
        all_genome_rows,
        per_genome_out,
        genome_fields,
    )

    write_tsv(
        taxon_summary_rows,
        taxon_summary_out,
        taxon_summary_fields,
    )

    write_tsv(
        distribution_rows,
        distribution_out,
        distribution_fields,
    )

    write_tsv(
        domain_summary_rows,
        domain_summary_out,
        domain_summary_fields,
    )

    print(f"[SAVE] {per_genome_out}", flush=True)
    print(f"[SAVE] {taxon_summary_out}", flush=True)
    print(f"[SAVE] {distribution_out}", flush=True)
    print(f"[SAVE] {domain_summary_out}", flush=True)


if __name__ == "__main__":
    main()

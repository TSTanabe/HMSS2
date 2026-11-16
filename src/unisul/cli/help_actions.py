import argparse
import csv
import re
import shutil
import textwrap
from collections import defaultdict
from pathlib import Path

from unisul.cli import paths


class FetchHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.print_fetch_operator_help()
        parser.exit(0)

    def print_fetch_operator_help(self) -> None:
        """
        Print a concise explanation of how -fd and -fc are interpreted,
        including AND/OR logic, brackets, optional domains, and combinations.
        """

        text = """
    FETCH OPERATORS: -fd and -fc
    ===========================

    -fd (fetch domains / proteins)
    ------------------------------
    Selects proteins by domain name and restricts results by genome-level co-occurrence.

    Core meaning:
      -fd A B means: select genomes where A and B co-occur in the SAME genome.
      Output then contains the matching proteins/domains from those genomes (subject to other filters).

    Rules:
      - Multiple values are combined using logical AND at the genome level (co-occurrence in one genome).
      - Use ':' (without whitespace) to express OR within a single slot.
      - A trailing ':' makes a domain optional (can be omitted for that slot).
      - Square brackets [...] define alternative requirement groups (logical OR between groups).

    Examples:
      -fd A B
          -> genomes containing A AND B (both must be present in the same genome)

      -fd A:B C
          -> genomes containing (A OR B) AND C

      -fd A:B:C D
          -> genomes containing (A OR B OR C) AND D

      -fd A:
          -> genomes containing A OR nothing (A is optional)

      -fd A: B
          -> genomes containing (A optional) AND B

      -fd [A B] [C D]
          -> (genomes with A AND B) OR (genomes with C AND D)


    -fc (fetch gene clusters / CSBs)
    --------------------------------
    Selects genomes by gene clusters (CSBs) that encode the specified domains.

    Core meaning:
      -fc A B means: select genomes that have at least one gene cluster containing A and B together.

    Rules:
      - Within a cluster definition, whitespace-separated tokens are combined using logical AND.
      - ':' (without whitespace) expresses OR within a single token (alternative domains for that slot).
      - A trailing ':' makes that token optional (can be omitted).
      - Square brackets [...] define alternative cluster definitions (logical OR between groups).

    Examples:
      -fc A B
          -> genomes with a cluster containing A AND B

      -fc A:B C
          -> genomes with a cluster containing (A OR B) AND C

      -fc A:B:C D
          -> genomes with a cluster containing (A OR B OR C) AND D

      -fc A:
          -> genomes with a cluster containing A OR no constraint for that slot

      -fc [A B] [C D]
          -> genomes with a cluster containing (A AND B) OR (C AND D)


    -fc and -fd combined
    -------------------
    When used together:

      -fc defines which genomes are selected based on gene clusters (hard genome filter).
      -fd then fetches/adds proteins by domain, but ONLY within the genomes selected by -fc.
      -fd cannot introduce new genomes beyond the -fc selection.

    In short:
      -fc defines the genome set;
      -fd extends the protein output within that set.


    -fnd (exclude domains)
    ----------------------
    Exclude results based on the presence of domains.

    Core meaning:
      -fnd A B means: exclude any genome, gene cluster, or protein where A or B is present.
      Exclusion is triggered as soon as any specified domain occurs.

    Rules:
      - All specified domains are treated as independent exclusion criteria.
      - Multiple values are always combined using logical OR.
      - If any one domain matches, the result is excluded.

    Examples:
      -fnd A
          -> exclude results containing A

      -fnd A B
          -> exclude results containing A or B

    Interaction with other options:
      - Exclusion is applied after genome or cluster selection.
      - Results matching -fnd are removed even if they satisfy -fd or -fc.

    ==============
    Important note
    --------------
    The fetch operators (-fd, -fc, -fnd) are evaluated only when the -r option
    is provided.

    The -r option must point to an existing result directory from a previous run.
    This directory is used as the input source for all fetch operations.

    If -r is not specified or does not point to a valid result directory:
      - fetch conditions are ignored
      - no fetch-based summaries, FASTA files, or graphs are generated

    """

        print(text.strip())


class FetchOutputHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.print_fetch_output_help()
        parser.exit(0)

    def print_fetch_output_help(self) -> None:
        text = """
    OUTPUT FILES
    ============

    The fetch operation produces tabular summary reports, FASTA sequence files,
    and graphical summaries. All outputs are derived from the same final,
    filtered result set after applying -fd, -fc and -fnd.

    SUMMARY TABLES (TSV)
    --------------------

    summary_hit_table.txt
      Main per-protein hit table (genomeID/proteinID, domains, coordinates, cluster, taxonomy).

    summary_gene_taxonomy.txt
      Protein-to-taxonomy mapping (proteinID -> full lineage string).

    summary_unique_lineages.txt
      Unique taxonomy lineages with genome counts.

    summary_hit_taxonomy_counts.txt
      Taxonomy-level presence summary for all fetched proteins.

    summary_requested_hit_taxonomy_counts.txt
      Same as above, restricted to proteins explicitly requested via -fd.

    summary_metabolic_annotations.txt
      Metabolic/functional annotations per protein/domain plus taxonomy.

    summary_genecluster_overview_table.txt
      Overview of gene cluster compositions associated with requested domains.

    summary_strain_variability_by_species.txt
      Species-level completeness/variability for requested domains.

    FASTA OUTPUT
    ------------

    Protein_sequences/
      <DOMAIN>.faa
        All proteins containing <DOMAIN>.
      multi_domain_<DOMAIN>.faa
        Domain sequences extracted from fusion proteins.
      <NAME>_ortho.faa / <NAME>_paralog.faa
        Split by single-copy vs multi-copy per genome.

    GRAPHICAL OUTPUT (JPG)
    ----------------------

    summary_hit_taxonomy_counts.jpg
      Presence/absence and co-occurrence visualization for all fetched proteins.

    summary_requested_hit_taxonomy_counts.jpg
      Presence/absence visualization for requested proteins only.

    summary_gene_taxonomy.jpg
      Protein-to-taxonomy assignments.

    summary_unique_lineages.jpg
      Unique lineage overview with genome counts.

    summary_metabolic_annotations.jpg
      Metabolic/functional annotation overview across taxa.

    summary_genecluster_overview_table.jpg
      Gene cluster composition overview.

    COMMAND-LINE RECORD
    -------------------

    command_line_args.txt
      Exact command-line arguments used for the run (index + value).
    """
        print(text.strip())


class ReadMappingHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.print_read_mapping_help()
        parser.exit(0)

    def print_read_mapping_help(self) -> None:
        # Defaults aus deinem Setup
        print("""
        READ MAPPING — AVAILABLE GPKG PACKAGES & RAM REQUIREMENTS
        =========================================================

        This list shows all available graftM packages (*.gpkg) for read mapping
        together with their estimated peak RAM usage (GB).

        RAM values are based on benchmark runs with:
          - threads = 8
          - reads   = 100
          - metric  = peak RSS (GB)

        Package sets are derived from the parent directory:
          gpkg/v10_gpkg_set_<SET>/.../<NAME>.gpkg


        SET: Apr_Qmo
        ------------
        AprM        : 10.47 GB
        QmoA        : 38.83 GB
        QmoB        : 49.98 GB
        QmoC        : 14.64 GB
        oxAprAI     : 24.14 GB
        oxAprAII    : 16.28 GB
        oxAprBI     : 4.53 GB
        oxAprBII    : 7.13 GB
        oxSat       : 21.58 GB
        qHdrB       : 6.71 GB
        qHdrC       : 11.57 GB
        redAprA     : 21.93 GB
        redAprB     : 2.19 GB
        redSat      : 16.59 GB


        SET: Asr_Mcc_Phs_Ttr
        -------------------
        AsrA        : 7.91 GB
        AsrB        : 6.48 GB
        AsrC        : 7.41 GB
        MccA        : 2.93 GB
        MccB        : 1.33 GB
        MccC        : 1.19 GB
        MccD        : 2.11 GB
        PhsA        : 29.35 GB
        PhsB        : 7.16 GB
        PhsC        : 12.16 GB
        TtrA        : 12.29 GB
        TtrB        : 3.64 GB
        TtrC        : 4.29 GB


        SET: CS_Aryl
        ------------
        AtsA        : 36.46 GB
        AtsB        : 32.47 GB
        CosH        : 61.94 GB
        Cs2H        : 1.24 GB
        DszA        : 5.77 GB
        DszB        : 1.38 GB
        DszC        : 4.97 GB
        DszD        : 0.25 GB
        SncA        : 0.58 GB
        SncB        : 0.66 GB
        SncC        : 0.97 GB
        SsuD        : 1.11 GB
        SsuE        : 0.66 GB
        TcdH        : 0.58 GB
        TcdS        : 0.54 GB
        TcdT        : 0.29 GB


        SET: DHPS_Taurine_Isethionate
        ----------------------------
        AdhE        : 15.01 GB
        ComC        : 49.28 GB
        ComD        : 2.70 GB
        ComE        : 3.13 GB
        HpfD        : 5.16 GB
        HpfG        : 27.43 GB
        HpfH        : 4.03 GB
        HpfX        : 6.42 GB
        HpfY        : 32.31 GB
        HpfZ        : 24.64 GB
        HpsG        : 26.94 GB
        HpsH        : 9.78 GB
        HpsN        : 46.88 GB
        HpsO        : 30.49 GB
        HpsP        : 27.35 GB
        IseJ        : 59.44 GB
        IsfD        : 11.50 GB
        SarD        : 5.08 GB
        SauS        : 10.17 GB
        SauT        : 15.21 GB
        SlcC        : 37.33 GB
        SlcD        : 21.90 GB
        SuyA        : 5.41 GB
        SuyB        : 29.91 GB
        TauX        : 2.25 GB
        TauY        : 7.98 GB
        Toa         : 162.91 GB
        Tpa         : 86.53 GB


        SET: DMS
        --------
        AcuI        : 17.89 GB
        AcuK        : 200.63 GB
        AcuN        : 294.56 GB
        DddA        : 73.16 GB
        DddC        : 108.97 GB
        DddD        : 4.19 GB
        DddL        : 0.25 GB
        DddP        : 4.20 GB
        DddQ        : 0.33 GB
        DddT        : 30.86 GB
        DddW        : 0.44 GB
        DddY        : 0.32 GB
        DdhA        : 2.46 GB
        DdhB        : 0.64 GB
        DdhC        : 0.29 GB
        DdhD        : 0.29 GB
        DmdA        : 5.11 GB
        DmdB        : 112.81 GB
        DmdD        : 5.24 GB
        DmoA        : 59.96 GB
        DmsA        : 115.08 GB
        DmsB        : 20.31 GB
        DmsC        : 2.60 GB
        DmsD        : 2.31 GB
        DorA        : 0.29 GB
        DorC        : 1.06 GB
        DorD        : 0.28 GB
        DsoA        : 0.27 GB
        DsoB        : 7.86 GB
        DsoC        : 1.56 GB
        DsoD        : 17.85 GB
        DsoE        : 0.99 GB
        DsoF        : 15.65 GB
        MarB        : 1.00 GB
        MarD        : 1.17 GB
        MarH        : 1.02 GB
        MarK        : 1.59 GB
        MddA        : 5.11 GB
        MddH        : 2.98 GB
        MsmA        : 3.87 GB
        MsmB        : 1.37 GB
        MsmC        : 0.70 GB
        MsmD        : 3.80 GB
        MsuC        : 58.78 GB
        MsuD        : 29.14 GB
        MsuE        : 5.20 GB
        Mtox        : 10.76 GB
        SnfG        : 16.97 GB


        SET: Sor_Soe
        ------------
        SoeA        : 103.17 GB
        SoeB        : 27.15 GB
        SoeC        : 35.04 GB
        SorA        : 5.48 GB
        SorB        : 1.46 GB


        SET: SQR
        --------
        CstA        : 66.29 GB
        CstB        : 2.46 GB
        SQRI        : 52.44 GB
        SQRII       : 57.44 GB
        SQRIII     : 176.69 GB
        SQRIV       : 28.31 GB
        SQRV        : 539.77 GB
        SQRVI       : 35.51 GB


        SET: SQ_SQDG
        ------------
        SftD        : 0.42 GB
        SftI        : 26.44 GB
        SftT        : 3.58 GB
        SftX        : 1.19 GB
        Sgdh        : 3.67 GB
        SmoB        : 10.29 GB
        SmoC        : 13.10 GB
        SqdA        : 0.22 GB
        SqdB        : 16.14 GB
        SqdC        : 4.12 GB
        SqdX        : 12.05 GB
        Sqald       : 1.35 GB
        Sqdh        : 1.79 GB
        SqgA        : 27.61 GB
        SqiA        : 1.03 GB
        SqiK        : 1.23 GB
        Sql         : 1.79 GB
        SqoD        : 1.23 GB
        SqvB        : 1.69 GB
        SqwD        : 1.51 GB
        SqwF        : 2.76 GB
        SqwG        : 4.58 GB
        SqwH        : 5.13 GB
        SqwI        : 1.67 GB
        SqwK        : 1.00 GB
        SqwL        : 0.84 GB
        YihQ        : 38.88 GB
        YihR        : 2.52 GB
        YihS        : 1.89 GB
        YihT        : 1.70 GB
        YihU        : 0.77 GB
        YihV        : 1.72 GB


        SET: sHdr_Dsr_Sox
        -----------------
        DoxA        : 0.29 GB
        DoxD        : 0.29 GB
        DsrE        : 10.34 GB
        DsrE3A      : 21.29 GB
        DsrE3B      : 0.95 GB
        DsrE3C      : 4.92 GB
        DsrF        : 9.80 GB
        DsrH        : 5.90 GB
        DsrL        : 34.79 GB
        DsrR        : 2.06 GB
        DsrS        : 1.60 GB
        EMO         : 43.40 GB
        FccA        : 3.75 GB
        FccB        : 10.52 GB
        LipS1       : 16.61 GB
        LipS2       : 16.91 GB
        LipT        : 12.86 GB
        LbpA1       : 1.07 GB
        LbpA2       : 7.61 GB
        Rhd442      : 0.90 GB
        SOR         : 1.47 GB
        SoxA        : 16.65 GB
        SoxB        : 47.30 GB
        SoxC        : 28.96 GB
        SoxD        : 16.64 GB
        SoxE        : 3.96 GB
        SoxF        : 14.71 GB
        SoxG        : 5.73 GB
        SoxH        : 8.91 GB
        SoxR        : 3.80 GB
        SoxS        : 4.40 GB
        SoxT1       : 16.11 GB
        SoxT2       : 15.13 GB
        SoxV        : 9.16 GB
        SoxW        : 7.06 GB
        SoxX        : 7.12 GB
        SoxY        : 6.79 GB
        SoxZ        : 7.09 GB
        oxDsrA      : 20.74 GB
        oxDsrB      : 16.19 GB
        oxDsrC      : 10.55 GB
        oxDsrJ      : 2.90 GB
        oxDsrK      : 30.78 GB
        oxDsrM      : 12.78 GB
        oxDsrN      : 22.31 GB
        oxDsrO      : 10.95 GB
        oxDsrP      : 18.26 GB
        redDsrA     : 20.40 GB
        redDsrB     : 19.39 GB
        redDsrC     : 4.73 GB
        redDsrD     : 2.17 GB
        redDsrE     : 2.16 GB
        redDsrF     : 1.86 GB
        redDsrH     : 1.97 GB
        redDsrJ     : 6.80 GB
        redDsrK     : 29.03 GB
        redDsrM     : 17.91 GB
        redDsrN     : 17.01 GB
        redDsrO     : 13.17 GB
        redDsrP     : 24.67 GB
        sHdrA       : 6.90 GB
        sHdrB1      : 7.08 GB
        sHdrB2      : 4.65 GB
        sHdrB3      : 3.20 GB
        sHdrC1      : 4.66 GB
        sHdrC2      : 3.76 GB
        sHdrH       : 5.13 GB
        sHdrI       : 1.53 GB
        sHdrT       : 9.40 GB
        sLplAB      : 26.01 GB
        """)


class PlacementPackageHelpAction(argparse.Action):
    """
    Show all available phylogenetic placement packages and the
    protein types contained in their annotations.tsv files.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        self.print_placement_package_help()
        parser.exit(0)

    @staticmethod
    def _read_package_annotations(placement_root: Path) -> dict[str, list[str]]:
        packages = {}

        if not placement_root.is_dir():
            raise FileNotFoundError(
                f"Placement package directory not found: {placement_root}"
            )

        for package_dir in sorted(
            placement_root.iterdir(),
            key=lambda path: path.name.casefold(),
        ):
            if not package_dir.is_dir():
                continue

            annotations = package_dir / "annotations.tsv"

            if not annotations.is_file():
                continue

            protein_types = set()

            with annotations.open(
                "r",
                encoding="utf-8-sig",
                newline="",
            ) as handle:
                reader = csv.reader(
                    handle,
                    delimiter="\t",
                )

                for row in reader:
                    if len(row) < 2:
                        continue

                    protein_type = row[1].strip()

                    if not protein_type or protein_type == "protein_type":
                        continue

                    protein_types.add(protein_type)

            if protein_types:
                packages[package_dir.name] = sorted(
                    protein_types,
                    key=lambda x: (
                        x.casefold() == "outgroup",
                        x.casefold(),
                    ),
                )

        return packages

    def print_placement_package_help(self) -> None:
        placement_root = paths.PLACEMENT_DIR

        packages = self._read_package_annotations(placement_root)

        if not packages:
            print(
                f"No placement packages with annotations.tsv found in: {placement_root}"
            )
            return

        print()
        print("AVAILABLE PHYLOGENETIC PLACEMENT PACKAGES")
        print("=========================================")
        print()
        print(
            "These placement packages provide the recommended phylogenetic classification\n"
            "workflow for UNISUL reads and protein sequences and are an alternative to the\n"
            "GraftM workflow. The pipeline combines HMM-based candidate detection,\n"
            "competitive MMseqs2 reference filtering, EPA-ng phylogenetic placement and\n"
            "GAPPA-based functional and taxonomic classification.\n"
        )

        print(
            "Genome workflow\n"
            "---------------\n"
            "Genome-derived proteins are first detected by the standard UNISUL HMM workflow.\n"
            "Candidate proteins are written to protein-type-specific intermediate-hit FASTA\n"
            "files and automatically matched to placement packages based on the protein\n"
            "types listed in annotations.tsv.\n"
            "\n"
            "MMseqs2 then compares each candidate against the package ingroup and, where\n"
            "available, outgroup reference sequences. Candidates are retained as ingroup\n"
            "when they show sufficient similarity to the ingroup and are not matched at\n"
            "least as well by the outgroup. High-confidence MMseqs2 matches can additionally\n"
            "promote valid hits and refine functional protein-type assignments. MMseqs2 also\n"
            "provides an independent taxonomic assignment based on the package references.\n"
            "\n"
            "If phylogenetic placement is enabled, the remaining intermediate hits are\n"
            "aligned to the package reference coordinate system, placed into the reference\n"
            "tree with EPA-ng and evaluated with GAPPA for functional and rank-aware\n"
            "taxonomic assignment.\n"
        )

        print(
            "Metagenome workflow\n"
            "-------------------\n"
            "Metagenomic FASTQ reads are first translated into protein ORFs with OrfM.\n"
            "Translated ORFs are streamed in batches through PyHMMER hmmscan against the\n"
            "combined search-HMM library derived from the available placement packages.\n"
            "The best HMM hit provides the initial protein type and matching ORFs are written\n"
            "to protein-type-specific intermediate-hit FASTA files.\n"
            "\n"
            "These candidate ORFs then undergo the same competitive MMseqs2 ingroup/outgroup\n"
            "comparison used for genome hits. Outgroup-like ORFs are removed before\n"
            "placement, while retained ingroup-like ORFs can be promoted to valid hits,\n"
            "functionally reassigned when sufficiently similar to a reference sequence, and\n"
            "taxonomically classified with MMseqs2.\n"
            "\n"
            "The remaining ORFs are then matched to the appropriate placement packages,\n"
            "aligned to the package reference alignment, placed with EPA-ng and subsequently\n"
            "classified with GAPPA. EPA-ng/GAPPA taxonomy is stored separately from MMseqs2\n"
            "taxonomy so that both assignment approaches remain available.\n"
        )

        print(
            "EPA-ng and MMseqs2 can use multiple CPU cores, and placement packages are\n"
            "processed independently, making this workflow well suited for parallel\n"
            "execution on large genome and metagenome datasets.\n"
        )

        print(
            "Each package is listed below together with the protein types contained in its\n"
            "annotations.tsv file."
        )
        print()

        package_width = max(
            len("Placement package"),
            max(len(name) for name in packages),
        )

        print(f"{'Placement package':<{package_width}}  Protein types")
        print(f"{'-' * package_width}  {'-' * 50}")

        for package_name, protein_types in packages.items():
            protein_text = ", ".join(protein_types)

            print(f"{package_name:<{package_width}}  {protein_text}")

        unique_types = {
            protein_type
            for protein_types in packages.values()
            for protein_type in protein_types
        }

        print()
        print(
            f"{len(packages)} placement packages | "
            f"{len(unique_types)} unique protein types"
        )
        print()


class ProteinAnnotationHelpAction(argparse.Action):
    """
    Show all HMMs/modules available for protein annotation.

    Information is read dynamically from Resource_module_metadata.
    """

    TRUE_VALUES = {"true", "1", "yes"}

    def __call__(self, parser, namespace, values, option_string=None):
        self.print_protein_annotation_help()
        parser.exit(0)

    @staticmethod
    def _as_int(value, default=999999):
        try:
            return int(value)
        except (TypeError, ValueError):
            return default

    def _read_metadata(self):
        metadata_file = paths.SRC_FILE_RESOURCE_METADATA

        if not metadata_file.is_file():
            raise FileNotFoundError(
                f"Resource module metadata not found: {metadata_file}"
            )

        # set -> module -> module information
        data = defaultdict(dict)

        with metadata_file.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            required = {
                "package",
                "set",
                "module",
                "module_name",
                "default",
                "hmms",
            }

            if not reader.fieldnames:
                raise ValueError("Resource module metadata has no header.")

            missing = required - set(reader.fieldnames)

            if missing:
                raise ValueError(
                    "Resource module metadata is missing columns: "
                    + ", ".join(sorted(missing))
                )

            for row in reader:
                package = row["package"].strip()
                set_name = row["set"].strip()
                module = row["module"].strip()
                module_name = row["module_name"].strip()

                set_order = self._as_int(row.get("set_order"))
                module_order = self._as_int(row.get("module_order"))

                is_default = row["default"].strip().lower() in self.TRUE_VALUES

                hmms = [
                    hmm.strip() for hmm in re.split(r"[,;]", row["hmms"]) if hmm.strip()
                ]

                if module not in data[set_name]:
                    data[set_name][module] = {
                        "module_name": module_name,
                        "set_order": set_order,
                        "module_order": module_order,
                        "hmms": defaultdict(
                            lambda: {
                                "packages": set(),
                                "defaults": set(),
                            }
                        ),
                    }

                module_entry = data[set_name][module]

                for hmm in hmms:
                    module_entry["hmms"][hmm]["packages"].add(package)

                    if is_default:
                        module_entry["hmms"][hmm]["defaults"].add(package)

        return data

    @staticmethod
    def _print_row(
        hmm,
        module,
        module_name,
        default,
        packages,
        widths,
        terminal_width,
    ):
        hmm_w, module_w, name_w, default_w = widths

        prefix = (
            f"{hmm:<{hmm_w}}  "
            f"{module:<{module_w}}  "
            f"{module_name:<{name_w}}  "
            f"{default:<{default_w}}  "
        )

        continuation = " " * len(prefix)

        package_width = max(
            15,
            terminal_width - len(prefix),
        )

        package_lines = textwrap.wrap(
            packages,
            width=package_width,
            break_long_words=False,
            break_on_hyphens=False,
        ) or [""]

        print(prefix + package_lines[0])

        for line in package_lines[1:]:
            print(continuation + line)

    def print_protein_annotation_help(self):
        data = self._read_metadata()

        if not data:
            print("No protein annotation metadata found.")
            return

        print()
        print("AVAILABLE PROTEIN ANNOTATION MODULES")
        print("====================================")
        print()

        print(
            "UNISUL uses a modular HMM-based annotation workflow for genome and protein\n"
            "sequence analysis. Protein families are organized into functional modules,\n"
            "and each module can be provided by one or more HMM packages. The module\n"
            "definition determines which HMMs are included in the search library.\n"
        )

        print(
            "Standard annotation workflow\n"
            "----------------------------\n"
            "For genome input, protein-coding genes are first predicted and translated.\n"
            "The resulting protein sequences are searched with PyHMMER against the HMMs\n"
            "selected from the annotation modules listed below. Each HMM has associated\n"
            "score thresholds that are used to distinguish high-confidence hits from\n"
            "intermediate and low-confidence candidates.\n"
            "\n"
            "The initial HMM search is deliberately sensitive and uses the noise cutoff.\n"
            "Detected domains are subsequently evaluated using their HMM score, alignment\n"
            "coverage and sequence identity. Hits reaching the trusted cutoff are directly\n"
            "accepted as valid annotations. Lower-scoring candidates can subsequently be\n"
            "supported by gene synteny, predefined collinear gene-cluster patterns and\n"
            "pathway/co-occurrence information. Candidate sequences are also retained in\n"
            "intermediate-hit FASTA files for downstream reference comparison and optional\n"
            "phylogenetic placement.\n"
        )

        print(
            "Annotation modules and HMM selection\n"
            "------------------------------------\n"
            "A module represents a biologically related annotation unit and may contain\n"
            "one or several HMMs. For example, a pathway or protein complex can be defined\n"
            "as one module while containing separate HMMs for its individual subunits.\n"
            "\n"
            "The Resource_module_metadata table defines which HMMs belong to each module,\n"
            "which package provides them, and which package variant is used by default.\n"
            "At startup, UNISUL selects every module/package combination marked as Default\n"
            "and builds the HMM search library from the HMMs associated with those selected\n"
            "module variants.\n"
            "\n"
            "Some modules are available in multiple packages. These package variants can\n"
            "contain alternative HMM implementations for the same biological module, for\n"
            "example models derived from different reference datasets or specialized\n"
            "classification workflows. Only the variants explicitly selected for a module\n"
            "are added to the search library.\n"
        )

        print(
            "Module selection\n"
            "----------------\n"
            "Without additional arguments, all module/package combinations marked as\n"
            "Default are selected automatically.\n"
            "\n"
            "  --add-module MODULE@PACKAGE\n"
            "      Add another package variant in addition to the currently selected\n"
            "      variant(s). The default variant remains active.\n"
            "\n"
            "  --replace-module MODULE@PACKAGE\n"
            "      Remove all currently selected variants of that module and use the\n"
            "      specified package variant instead.\n"
            "\n"
            "Examples:\n"
            "\n"
            "  --add-module redDsr@DiSCo\n"
            "      Adds the DiSCo redDsr HMMs in addition to the default redDsr models.\n"
            "\n"
            "  --replace-module redDsr@DiSCo\n"
            "      Removes the default redDsr package variant and uses the DiSCo variant.\n"
            "\n"
            "Multiple --add-module and --replace-module arguments can be supplied in the\n"
            "same run. The value in the 'Module' column below is the MODULE identifier\n"
            "used for these command-line arguments.\n"
        )

        print(
            "Understanding the table\n"
            "-----------------------\n"
            "HMM\n"
            "    Protein/domain HMM included by this module variant.\n"
            "\n"
            "Module\n"
            "    Logical annotation module used with --add-module and --replace-module.\n"
            "\n"
            "Module name\n"
            "    Human-readable description of the biological module.\n"
            "\n"
            "Default\n"
            "    Package variant selected automatically by the standard workflow.\n"
            "\n"
            "Packages\n"
            "    All packages that provide the corresponding HMM for this module.\n"
        )

        # ---------------------------------------------------------
        # Determine fixed column widths globally
        # ---------------------------------------------------------

        all_hmms = []
        all_modules = []
        all_module_names = []
        all_defaults = []

        for modules in data.values():
            for module, info in modules.items():
                all_modules.append(module)
                all_module_names.append(info["module_name"])

                for hmm, hmm_info in info["hmms"].items():
                    all_hmms.append(hmm)

                    default = (
                        ", ".join(
                            sorted(
                                hmm_info["defaults"],
                                key=str.casefold,
                            )
                        )
                        or "-"
                    )

                    all_defaults.append(default)

        hmm_w = min(
            20,
            max(len("HMM"), *(len(x) for x in all_hmms)),
        )

        module_w = min(
            18,
            max(len("Module"), *(len(x) for x in all_modules)),
        )

        name_w = min(
            42,
            max(
                len("Module name"),
                *(len(x) for x in all_module_names),
            ),
        )

        default_w = min(
            15,
            max(
                len("Default"),
                *(len(x) for x in all_defaults),
            ),
        )

        widths = (
            hmm_w,
            module_w,
            name_w,
            default_w,
        )

        terminal_width = shutil.get_terminal_size(fallback=(140, 24)).columns

        sorted_sets = sorted(
            data.items(),
            key=lambda item: (
                min(module["set_order"] for module in item[1].values()),
                item[0].casefold(),
            ),
        )

        for set_name, modules in sorted_sets:
            print()
            print(f"SET: {set_name}")
            print("=" * min(terminal_width, 120))

            header = (
                f"{'HMM':<{hmm_w}}  "
                f"{'Module':<{module_w}}  "
                f"{'Module name':<{name_w}}  "
                f"{'Default':<{default_w}}  "
                f"Packages"
            )

            print(header)

            print(
                f"{'-' * hmm_w}  "
                f"{'-' * module_w}  "
                f"{'-' * name_w}  "
                f"{'-' * default_w}  "
                f"{'-' * 20}"
            )

            sorted_modules = sorted(
                modules.items(),
                key=lambda item: (
                    item[1]["module_order"],
                    item[1]["module_name"].casefold(),
                    item[0].casefold(),
                ),
            )

            for module, info in sorted_modules:
                first_hmm = True

                for hmm in sorted(
                    info["hmms"],
                    key=str.casefold,
                ):
                    hmm_info = info["hmms"][hmm]

                    packages = ", ".join(
                        sorted(
                            hmm_info["packages"],
                            key=str.casefold,
                        )
                    )

                    default = (
                        ", ".join(
                            sorted(
                                hmm_info["defaults"],
                                key=str.casefold,
                            )
                        )
                        or "-"
                    )

                    if first_hmm:
                        module_display = module
                        name_display = info["module_name"]
                        first_hmm = False
                    else:
                        module_display = ""
                        name_display = ""

                    self._print_row(
                        hmm=hmm,
                        module=module_display,
                        module_name=name_display,
                        default=default,
                        packages=packages,
                        widths=widths,
                        terminal_width=terminal_width,
                    )

            print()

# hmsss/core/config.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Any
from pathlib import Path

"""
Typed configuration model (dataclasses) for HMSS2/HMSSS.

The `Config` aggregate captures:
- CLI parameter groups (input, search, resources, synteny, info, CSB, flow,
  limiters, operators, processing),
- resolved project paths,
- project-specific output fields,
- runtime state that grows during the pipeline.

It also provides convenience properties that proxy into nested dataclasses
(e.g., `config.stage`, `config.patterns_file`) and a `validate()` method
for basic consistency checks.
"""


@dataclass(slots=True)
class PathsCfg:
    """Canonical project directories resolved from `hmsss.cli.paths`.

    Attributes:
        root: Project root directory.
        bin:  Executables/tools directory.
        data: Top-level data directory.
        gpkg: Top-level gpkg directory.
        hmms: Directory containing HMM libraries.
        refseq: Directory with reference sequences.
        results: Default results directory.
        package: Python package root (`src/hmsss`).
    """

    root: str
    bin: str
    data: str
    gpkg: str
    hmms: str
    refseq: str
    results: str
    package: str


# ==========================================
# 1) CLI-Argumentgruppen (thematisch getrennt)
# ==========================================
@dataclass(slots=True)
class CliInput:
    """Input and I/O related CLI options.

    Attributes:
        fasta_file_directory: Directory to search for genomes.
        score_threshold_file: TSV file with optimized/trusted/noise cutoffs.
        library: HMM library path.
        result_files_directory: Output directory for results.
        database_directory: Path to SQLite database (created if missing).
        cores: Number of CPU cores to use.
        glob_report: Path to global hmmreport.
        verbose: Logging level (0=WARNING, 1=INFO, 2=DEBUG).
    """

    fasta_file_directory: Optional[str] = None
    score_threshold_file: Optional[str] = None
    library: Optional[str] = None
    result_files_directory: Optional[str] = None
    database_directory: Optional[str] = None
    cores: int = 4
    glob_report: Optional[str] = None
    verbose: int = 1


@dataclass(slots=True)
class CliSearchParams:
    """Search parameters and global thresholds.

    Attributes:
        threshold_type: 1=optimized, 2=trusted, 3=noise.
        thrs_score: Minimal global score cutoff.
        taxonomy_file: Optional TSV with taxonomy.
        refseq_identity: Minimal % identity vs. reference set.
        name: Project name.
        stage: Pipeline stage to start from.
        exit: Pipeline stage to exit after.
    """

    threshold_type: int = 1  # 1=optimized, 2=trusted, 3=noise
    thrs_score: float = 50.0
    threshold_factor: float = 1.0
    taxonomy_file: Optional[str] = None
    refseq_identity: int = 90
    name: str = "project"
    stage: int = 0
    exit: int = 10


@dataclass(slots=True)
class CliResources:
    """Resources and toggles for the search library and cross-checks.

    Attributes:
        HMM_sets: Optional subset of HMM sets to include.
        clean_reports: Overwrite existing hmmsearch reports.
        disable_individual_reports: Write per-genome reports if True.
        max_seqs_per_genome: Cap sequences per protein per genome (Diamond check).
        bool_cross_check: Enable reference cross-check via Diamond.
        optimized_cutoff_cross_check: Use optimized cutoff instead of Diamond.
    """

    HMM_sets: List[str] = field(default_factory=list)
    clean_reports: bool = False
    disable_individual_reports: bool = True
    max_seqs_per_genome: int = 4
    diamond_speed_mode: str = "fast"
    bool_cross_check: bool = True
    optimized_cutoff_cross_check: bool = False


@dataclass(slots=True)
class CliSynteny:
    """Synteny and pattern options for CSB naming.

    Attributes:
        patterns_file: Path to pattern definitions.
        cooccurrence_file: Path to co-occurrence matrix/file.
        exclusion_singletons: Proteins to exclude when unclustered.
        min_completeness: Minimal fraction of CSB required for recognition.
        glob_chunks: Chunk size for parsing `glob` results before DB insert.
    """

    patterns_file: Optional[str] = None
    cooccurrence_file: Optional[str] = None
    exclusion_singletons: Optional[str] = None
    min_completeness: float = 0.5
    glob_chunks: int = 5000


@dataclass(slots=True)
class CliInfo:
    """Toggles to print auxiliary statistics."""

    stat_keywords: bool = False
    stat_csb: bool = False
    stat_genomes: bool = False
    metabolic_information: Optional[str] = None


@dataclass(slots=True)
class CliCsb:
    """Parameters for collinear syntenic block detection.

    Attributes:
        nucleotide_range: Max nucleotide distance for synteny.
        insertions: Allowed insertions within a CSB.
        occurence: Minimal occurrences to accept a CSB (spelling kept).
        min_csb_size: Minimum number of genes per CSB.
        max_csb_size: Maximum number of genes per CSB.
        max_domain_repeats: Max repeats for a domain within a CSB.
        jaccard: Dissimilarity threshold [0.0–1.0].
    """

    nucleotide_range: int = 3500
    insertions: int = 1
    occurence: int = 1  # Schreibweise wie im Parser beibehalten
    min_csb_size: int = 4
    max_csb_size: int = 50
    max_domain_repeats: int = 4
    jaccard: float = 0.0


@dataclass(slots=True)
class CliReadMapping:
    """CLI options for the hidden read-mapping module."""

    use_read_mapping: bool = False
    gpkg_sets: List[str] = field(default_factory=list)
    gpkg_packs: List[str] = field(default_factory=list)
    threads: int = 5
    evalue: float = 1e-5
    placements_cutoff: float = 0.75
    resolve_placements: bool = False

    min_orf_length: int = 96
    restrict_read_length: Optional[int] = None
    translation_table: int = 11

    ram_limit_min: Optional[float] = None
    ram_limit_max: Optional[float] = None

    interleaved: bool = False

    ram_profile_file: Optional[str] = None


@dataclass(slots=True)
class CliFlow:
    """Global flow modifiers (e.g., recompute taxonomy)."""

    redo_taxonomy: bool = False
    disable_synteny_completion: bool = True
    use_remove_unassigned_intermediates: bool = True
    use_remove_exclusion_singletons: bool = True


@dataclass(slots=True)
class CliLimiter:
    """Dataset filters to constrain DB output.

    Attributes:
        dataset_limit_lineage: Taxonomic rank name.
        dataset_limit_taxon: Specific taxon string.
        dataset_limit_proteins: Protein/domain filter expression.
        dataset_limit_keywords: Keyword filter expression.
        dataset_divide_sign: Separator used in taxonomy strings.
    """

    dataset_limit_lineage: Optional[str] = None
    dataset_limit_taxon: Optional[str] = None
    dataset_limit_proteins: str = "0"
    dataset_limit_keywords: str = "0"
    dataset_divide_sign: str = "."


@dataclass(slots=True)
class CliOperators:
    """Operators controlling which entities are fetched from the DB.

    Attributes:
        fetch_genomes: List of genome identifiers.
        fetch_proteins: Protein domain names.
        fetch_csbs: CSB identifiers.
        fetch_keywords: Cluster naming keywords.
        keywords_connector: Logical connector for keyword filters ("AND"/"OR").
    """

    fetch_genomes: List[str] = field(default_factory=list)
    fetch_proteins: List[str] = field(default_factory=list)
    fetch_csbs: List[str] = field(default_factory=list)
    fetch_keywords: List[str] = field(default_factory=list)
    fetch_not_csb_with_these_domains: List[str] = field(default_factory=list)
    keywords_connector: str = "OR"
    print_fasta: bool = False
    print_graphs: bool = False
    use_non_valid_hits: bool = True
    graph_tax_levels: List[str] = field(default_factory=lambda: ["Phylum"])


# ==================================================
# 2) Internal runtime states
# ==================================================
@dataclass(slots=True)
class RuntimeState:
    """Mutable runtime state filled as the pipeline proceeds.

    Tracks queued/finished genomes, discovered files, and redundancy sets/maps.
    """

    queued_genomes: Set[str] = field(default_factory=set)
    finished_genomes: Set[str] = field(default_factory=set)

    fna_files: Dict[str, str] = field(default_factory=dict)
    faa_files: Dict[str, str] = field(default_factory=dict)
    gff_files: Dict[str, str] = field(default_factory=dict)
    hmmreport_files: Dict[str, str] = field(default_factory=dict)
    fastq_files: Dict[str, str] = field(default_factory=dict)

    redundant: Dict[str, List[str]] = field(default_factory=dict)
    non_redundant: Dict[str, str] = field(default_factory=dict)
    redundancy_hash: Dict[str, str] = field(default_factory=dict)


# ==================================================
# 3) Project specific attributes
# ==================================================
@dataclass(slots=True)
class ProjectFields:
    """Project-derived output locations and filenames.

    Populated by project setup code after parsing (in project.py)
    """

    result_files_directory: Optional[str] = None  # finaler Projekt-Results-Pfad
    fasta_initial_hit_directory: Optional[str] = None
    fasta_output_directory: Optional[str] = None
    cross_check_directory: Optional[str] = None
    csb_directory: Optional[str] = None

    # Dateien
    database_directory: Optional[str] = None
    glob_report: Optional[str] = None
    glob_trusted_hitreport: Optional[str] = None
    glob_intermediate_hitreport: Optional[str] = None
    csb_output_file: Optional[str] = None
    gene_clusters_file: Optional[str] = None


# ======================
# Aggregiertes Objekt
# ======================


def _cfg_get(obj: Any, path: str) -> Any:
    cur = obj
    for name in path.split("."):
        cur = getattr(cur, name)
    return cur


def _cfg_set(obj: Any, path: str, value: Any) -> None:
    parts = path.split(".")
    cur = obj
    for name in parts[:-1]:
        cur = getattr(cur, name)
    setattr(cur, parts[-1], value)


def prop(path: str) -> property:
    """Create a property proxy to a nested dataclass field using a dotted path.

    Example:
        `stage = prop("cli_params.stage")` makes `config.stage` access
        `config.cli_params.stage`.

    Args:
        path: Dotted path to nested field.

    Returns:
        A Python `property` object that gets/sets the nested value.
    """

    def fget(self):
        return _cfg_get(self, path)

    def fset(self, v):
        _cfg_set(self, path, v)

    return property(fget, fset)


@dataclass(slots=True)
class Config:
    """Central configuration aggregate for HMSS2/HMSSS.

    This class holds all CLI groups, runtime state, and project fields.
    It exposes convenience properties (created with `prop`) for frequently
    accessed fields to avoid deep attribute chains in the pipeline code.
    """

    # Pflicht: Pfade müssen zur Initialisierung übergeben werden
    paths: PathsCfg

    # CLI-Gruppen
    cli_input: CliInput = field(default_factory=CliInput)
    cli_params: CliSearchParams = field(default_factory=CliSearchParams)
    cli_resources: CliResources = field(default_factory=CliResources)
    cli_synteny: CliSynteny = field(default_factory=CliSynteny)
    cli_info: CliInfo = field(default_factory=CliInfo)
    cli_csb: CliCsb = field(default_factory=CliCsb)
    cli_flow: CliFlow = field(default_factory=CliFlow)
    cli_limiter: CliLimiter = field(default_factory=CliLimiter)
    cli_ops: CliOperators = field(default_factory=CliOperators)
    cli_readmap: CliReadMapping = field(default_factory=CliReadMapping)

    # Laufzeit-State & Projektfelder
    state: RuntimeState = field(default_factory=RuntimeState)
    project: ProjectFields = field(default_factory=ProjectFields)

    # Häufig genutzte CLI-Parameter (kurzer Zugriff)
    stage = prop("cli_params.stage")  # :contentReference[oaicite:1]{index=1}
    name = prop("cli_params.name")  # :contentReference[oaicite:2]{index=2}
    exit = prop("cli_params.exit")  # :contentReference[oaicite:3]{index=3}
    thrs_score = prop("cli_params.thrs_score")  # :contentReference[oaicite:4]{index=4}
    threshold_factor = prop("cli_params.threshold_factor")
    threshold_type = prop(
        "cli_params.threshold_type"
    )  # :contentReference[oaicite:5]{index=5}
    refseq_identity = prop("cli_params.refseq_identity")
    taxonomy_file = prop("cli_params.taxonomy_file")

    fasta_file_directory = prop(
        "cli_input.fasta_file_directory"
    )  # :contentReference[oaicite:6]{index=6}
    database_in = prop(
        "cli_input.database_directory"
    )  # Eingabe-DB (CLI) :contentReference[oaicite:7]{index=7}
    cli_result_dir_in = prop(
        "cli_input.result_files_directory"
    )  # Eingabe-Results (CLI) :contentReference[oaicite:8]{index=8}
    cores = prop("cli_input.cores")  # :contentReference[oaicite:9]{index=9}
    verbose = prop("cli_input.verbose")  # :contentReference[oaicite:10]{index=10}
    score_threshold_file = prop(
        "cli_input.score_threshold_file"
    )  # :contentReference[oaicite:11]{index=11}
    library = prop("cli_input.library")  # :contentReference[oaicite:12]{index=12}
    glob_report_in = prop(
        "cli_input.glob_report"
    )  # :contentReference[oaicite:13]{index=13}

    # Ressourcen / Operators / Limiter
    hmm_sets = prop("cli_resources.HMM_sets")  # :contentReference[oaicite:14]{index=14}
    clean_reports = prop("cli_resources.clean_reports")
    bool_cross_check = prop("cli_resources.bool_cross_check")
    disable_individual_reports = prop("cli_resources.disable_individual_reports")
    max_seqs_per_genome = prop("cli_resources.max_seqs_per_genome")
    diamond_speed_mode = prop("cli_resources.diamond_speed_mode")
    optimized_cutoff_cross_check = prop("cli_resources.optimized_cutoff_cross_check")

    # Database operations
    keywords_connector = prop(
        "cli_ops.keywords_connector"
    )  # :contentReference[oaicite:15]{index=15}
    fetch_genomes = prop(
        "cli_ops.fetch_genomes"
    )  # :contentReference[oaicite:16]{index=16}
    fetch_proteins = prop(
        "cli_ops.fetch_proteins"
    )  # :contentReference[oaicite:17]{index=17}
    fetch_csbs = prop("cli_ops.fetch_csbs")  # :contentReference[oaicite:18]{index=18}
    fetch_keywords = prop(
        "cli_ops.fetch_keywords"
    )  # :contentReference[oaicite:19]{index=19}
    fetch_not_csb_with_these_domains = prop("cli_ops.fetch_not_csb_with_these_domains")
    print_fasta = prop("cli_ops.print_fasta")
    print_graphs = prop("cli_ops.print_graphs")
    use_non_valid_hits = prop("cli_ops.use_non_valid_hits")
    graph_tax_levels = prop("cli_ops.graph_tax_levels")

    dataset_limit_lineage = prop(
        "cli_limiter.dataset_limit_lineage"
    )  # :contentReference[oaicite:20]{index=20}
    dataset_limit_taxon = prop(
        "cli_limiter.dataset_limit_taxon"
    )  # :contentReference[oaicite:21]{index=21}
    dataset_limit_proteins = prop("cli_limiter.dataset_limit_proteins")
    dataset_limit_keywords = prop("cli_limiter.dataset_limit_keywords")
    dataset_divide_sign = prop(
        "cli_limiter.dataset_divide_sign"
    )  # :contentReference[oaicite:22]{index=22}

    # Synteny / CSB / Process
    patterns_file = prop(
        "cli_synteny.patterns_file"
    )  # :contentReference[oaicite:23]{index=23}
    cooccurrence_file = prop(
        "cli_synteny.cooccurrence_file"
    )  # :contentReference[oaicite:24]{index=24}
    exclusion_singletons = prop(
        "cli_synteny.exclusion_singletons"
    )  # :contentReference[oaicite:25]{index=25}
    min_completeness = prop(
        "cli_synteny.min_completeness"
    )  # :contentReference[oaicite:26]{index=26}
    glob_chunks = prop(
        "cli_synteny.glob_chunks"
    )  # :contentReference[oaicite:27]{index=27}

    nucleotide_range = prop(
        "cli_csb.nucleotide_range"
    )  # :contentReference[oaicite:28]{index=28}
    insertions = prop("cli_csb.insertions")  # :contentReference[oaicite:29]{index=29}
    occurence = prop("cli_csb.occurence")  # :contentReference[oaicite:30]{index=30}
    min_csb_size = prop(
        "cli_csb.min_csb_size"
    )  # :contentReference[oaicite:31]{index=31}
    max_csb_size = prop(
        "cli_csb.max_csb_size"
    )  # :contentReference[oaicite:32]{index=32}
    max_domain_repeats = prop(
        "cli_csb.max_domain_repeats"
    )  # :contentReference[oaicite:33]{index=33}
    jaccard = prop("cli_csb.jaccard")  # :contentReference[oaicite:34]{index=34}

    # Read mapping (shortcuts)
    use_read_mapping = prop("cli_readmap.use_read_mapping")
    gpkg_sets = prop("cli_readmap.gpkg_sets")
    gpkg_packs = prop("cli_readmap.gpkg_packs")
    rm_threads = prop("cli_readmap.threads")
    rm_evalue = prop("cli_readmap.evalue")
    rm_placements_cutoff = prop("cli_readmap.placements_cutoff")
    rm_resolve_placements = prop("cli_readmap.resolve_placements")
    rm_min_orf_length = prop("cli_readmap.min_orf_length")
    rm_restrict_read_length = prop("cli_readmap.restrict_read_length")
    rm_translation_table = prop("cli_readmap.translation_table")
    rm_ram_limit_min = prop("cli_readmap.ram_limit_min")
    rm_ram_limit_max = prop("cli_readmap.ram_limit_max")
    rm_interleaved = prop("cli_readmap.interleaved")
    rm_ram_profile = prop("cli_readmap.ram_profile_file")

    merge_fasta = prop(
        "cli_process.merge_fasta"
    )  # :contentReference[oaicite:35]{index=35}
    filter_fasta = prop(
        "cli_process.filter_fasta"
    )  # :contentReference[oaicite:36]{index=36}
    concat_alignment = prop(
        "cli_process.concat_alignment"
    )  # :contentReference[oaicite:37]{index=37}
    add_taxonomy = prop(
        "cli_process.add_taxonomy"
    )  # :contentReference[oaicite:38]{index=38}
    add_genomic_context = prop(
        "cli_process.add_genomic_context"
    )  # :contentReference[oaicite:39]{index=39}
    create_type_range_dataset = prop(
        "cli_process.create_type_range_dataset"
    )  # :contentReference[oaicite:40]{index=40}
    create_gene_cluster_dataset = prop(
        "cli_process.create_gene_cluster_dataset"
    )  # :contentReference[oaicite:41]{index=41}
    gaps = prop("cli_process.gaps")  # :contentReference[oaicite:42]{index=42}

    # Projekt-Ausgaben (vom Project-Setup befüllt)
    result_files_directory = prop(
        "project.result_files_directory"
    )  # :contentReference[oaicite:43]{index=43}
    database_directory = prop(
        "project.database_directory"
    )  # :contentReference[oaicite:44]{index=44}
    glob_report = prop("project.glob_report")  # :contentReference[oaicite:45]{index=45}
    fasta_initial_hit_directory = prop(
        "project.fasta_initial_hit_directory"
    )  # :contentReference[oaicite:46]{index=46}
    fasta_output_directory = prop(
        "project.fasta_output_directory"
    )  # :contentReference[oaicite:47]{index=47}
    cross_check_directory = prop(
        "project.cross_check_directory"
    )  # :contentReference[oaicite:48]{index=48}
    csb_directory = prop(
        "project.csb_directory"
    )  # :contentReference[oaicite:49]{index=49}
    csb_output_file = prop(
        "project.csb_output_file"
    )  # :contentReference[oaicite:50]{index=50}
    gene_clusters_file = prop(
        "project.gene_clusters_file"
    )  # :contentReference[oaicite:51]{index=51}
    glob_trusted_hitreport = prop("project.glob_trusted_hitreport")
    glob_intermediate_hitreport = prop("project.glob_intermediate_hitreport")

    # Laufzeit und State Ausgaben
    queued_genomes = prop("state.queued_genomes")
    finished_genomes = prop("state.finished_genomes")
    fna_files = prop("state.fna_files")
    faa_files = prop("state.faa_files")
    gff_files = prop("state.gff_files")
    fastq_files = prop("state.fastq_files")
    hmmreport_files = prop("state.hmmreport_files")
    redundant = prop("state.redundant")
    non_redundant = prop("state.non_redundant")
    redundancy_hash = prop("state.redundancy_hash")

    # Metabolic information and statistics
    metabolic_information = prop("cli_info.metabolic_information")

    # Flow control for parsing and redo taxonomy
    redo_taxonomy = prop("cli_flow.redo_taxonomy")
    disable_synteny_completion = prop("cli_flow.disable_synteny_completion")
    use_remove_unassigned_intermediates = prop(
        "cli_flow.use_remove_unassigned_intermediates"
    )
    use_remove_unassigned_singletons = prop("cli_flow.use_remove_exclusion_singletons")

    def validate(self) -> None:
        """Run basic consistency checks and ensure required directories exist.

        Raises:
            ValueError: On invalid parameter ranges (e.g., `jaccard` not in [0,1]).
            FileNotFoundError: If one of the canonical project directories is missing.
        """
        if self.cli_params.threshold_type not in (1, 2, 3):
            raise ValueError(
                "threshold_type must be 1 (optimized), 2 (trusted) or 3 (noise)"
            )
        if self.cli_params.threshold_factor < 0.0:
            raise ValueError("threshold_factor must be >= 0.0")

        if self.cli_ops.keywords_connector not in ("AND", "OR"):
            raise ValueError("keywords_connector must be 'AND' or 'OR'")
        if (
                self.cli_synteny.min_completeness < 0.0
                or self.cli_synteny.min_completeness > 1.0
        ):
            raise ValueError("min_completeness must be within [0.0, 1.0]")
        if self.cli_csb.jaccard < 0.0 or self.cli_csb.jaccard > 1.0:
            raise ValueError("jaccard must be within [0.0, 1.0]")
        if self.cli_input.cores < 1:
            raise ValueError("cores must be >= 1")
        if (self.cli_readmap.gpkg_sets or self.cli_readmap.gpkg_packs) and not self.cli_readmap.use_read_mapping:
            self.stage = 50
            self.cli_readmap.use_read_mapping = True

        # Pfade der Basisstruktur prüfen
        missing: list[str] = []

        def _req_dir(path: str, label: str) -> None:
            if not Path(path).is_dir():
                missing.append(f"{label}: {path}")

        _req_dir(self.paths.root, "ROOT_DIR")
        _req_dir(self.paths.bin, "BIN_DIR")
        _req_dir(self.paths.data, "DATA_DIR")
        # _req_dir(self.paths.hmms, "HMMS_DIR")
        _req_dir(self.paths.refseq, "REFSEQ_DIR")
        _req_dir(self.paths.results, "RESULTS_DIR")

        if missing:
            details = "\n - ".join(missing)
            raise FileNotFoundError("Required directories are missing:\n - " + details)

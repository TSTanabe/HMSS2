# src/hmsss/core/options.py
from __future__ import annotations

from typing import Optional, List, Dict, Any
from hmsss.utils.paths import PACKAGE_ROOT

# Kompatibler Ersatz für das frühere __location__
__location__ = str(PACKAGE_ROOT)


class Hmsss:
    """
    Zentrales Options-/Konfigurationsobjekt für HMSSS.

    - Alle Attribute sind vorab vorhanden (mit passenden Defaults),
      sodass IDE-Inspektionen (PyCharm) keine "unresolved attribute"-Warnungen mehr zeigen.
    - argparse schreibt bei parse_arguments(...) in dieses Objekt (namespace=opts)
      und überschreibt die hier gesetzten Defaults bei Bedarf.
    """

    def __init__(
        self,
        # --- Bisherige HMSSS-Defaults / interne Felder ---
        execute_location: str = __location__,
        project_name: str = "project",
        index_db: bool = False,
        csb_name_prefix: str = "csb-",
        csb_name_suffix: str = "_",
        genome_id_divider: str = "___",
        # Laufzeit-/Queue-Container
        finished_genomes: Optional[Dict[str, Any]] = None,
        queued_genomes: Optional[Dict[str, Any]] = None,
        fna_files: Optional[Dict[str, str]] = None,
        faa_files: Optional[Dict[str, str]] = None,
        gff_files: Optional[Dict[str, str]] = None,
        hmmreport_files: Optional[Dict[str, str]] = None,
        # Redundanz
        redundant: int = 0,
        non_redundant: int = 0,
        redundancy_hash: Optional[Dict[str, Any]] = None,
        # Pipeline-Flags (werden teilweise in parse_arguments gesetzt)
        limiter: bool = False,
        fetch: bool = False,
        process: bool = False,
        # --- argparse: Input definition ---
        fasta_file_directory: Optional[str] = None,  # -f
        score_threshold_file: Optional[str] = None,  # -t
        library: Optional[str] = None,  # -l
        result_files_directory: str = __location__ + "/results",  # -r
        database_directory: Optional[str] = None,  # -db
        cores: int = 4,  # -c
        glob_report: Optional[str] = None,  # -glob_report
        verbose: int = 1,  # -v/--verbose (0,1,2)
        # --- argparse: Search parameters ---
        threshold_type: int = 1,  # -cut_type (1,2,3)
        thrs_score: int = 50,  # -cut_score
        taxonomy_file: Optional[str] = None,  # -taxonomy_info
        refseq_identity: int = 90,  # -refseq_ident
        name: str = "project",  # -n
        start_stage: int = 0,  # -s
        exit_stage: int = 10,  # -x
        # --- argparse: Search library resources ---
        hmm_sets: Optional[List[str]] = None,  # -hmms (list)
        clean_reports: bool = False,  # -clean (store_true)
        individual_reports: bool = True,  # -no_reports (store_false)
        max_seqs_per_genome: int = 4,  # -max_seq_per_genome
        bool_cross_check: bool = True,  # -no_cross_check (store_false)
        optimized_cutoff_cross_check: bool = False,  # -optimized_cutoff_cross_check
        # --- argparse: Synteny options ---
        patterns_file: Optional[str] = None,  # -p
        cooccurrence_file: Optional[str] = None,  # -cooccurrence
        exclusion_singletons: Optional[str] = None,  # -exclude_singletons
        min_completeness: float = 0.5,  # -mc
        glob_chunks: int = 5000,  # -chunks
        # --- argparse: Information flags ---
        stat_keywords: bool = False,  # -stat_keywords
        stat_csb: bool = False,  # -stat_csb
        stat_genomes: bool = False,  # -stat_genomes
        # --- argparse: CSB prediction ---
        nucleotide_range: int = 3500,  # -nt
        insertions: int = 1,  # -insertions
        occurence: int = 1,  # -occurence (sic)
        min_csb_size: int = 4,  # -min_csb_size
        max_csb_size: int = 50,  # -max_csb_size
        max_domain_repeats: int = 4,  # -max_gene_repeats
        jaccard: float = 0.0,  # -jaccard
        # --- argparse: Work step regulation ---
        redo_taxonomy: bool = False,  # -redo_taxonomy
        # --- argparse: Limiter / Operators ---
        dataset_limit_lineage: Optional[str] = None,  # -dll / -fl
        dataset_limit_taxon: Optional[str] = None,  # -dlt / -ft
        dataset_limit_proteins: str = "0",  # -dlp (Parser-Default "0")
        dataset_limit_keywords: str = "0",  # -dlk (Parser-Default "0")
        dataset_divide_sign: str = ".",  # -dtd
        # (in Code teils referenziert, im Parser evtl. nicht mehr vorhanden)
        dataset_limit_min_cluster_completeness: Optional[float] = None,
        fetch_genomes: Optional[List[str]] = None,  # -fg
        fetch_proteins: Optional[List[str]] = None,  # -fd
        fetch_csbs: Optional[List[str]] = None,  # -fc
        fetch_keywords: Optional[List[str]] = None,  # -fk
        keywords_connector: str = "OR",  # -kc ("AND"|"OR")
        # --- argparse: Alignment / sequence processing ---
        merge_fasta: Optional[str] = None,  # -merge_fasta (dir)
        # Hinweis: In deinem process_operator wird filter_fasta wie eine Sequenz verwendet.
        # Wenn du sicher gehen willst, kannst du hier List[str] annehmen und im Parser entsprechend setzen.
        filter_fasta: Optional[Any] = None,  # -filter_fasta (file,int,int)
        concat_alignment: Optional[str] = None,  # -concat_alignment
        add_taxonomy: Optional[str] = None,  # -add_taxonomy
        add_genomic_context: Optional[str] = None,  # -add_genomic_context
        create_type_range_dataset: Optional[str] = None,  # -create_type_range_dataset
        create_gene_cluster_dataset: Optional[
            str
        ] = None,  # -create_gene_cluster_dataset
        gaps: bool = False,  # -aln_gaps
        # --- Weitere intern genutzte Felder (nicht direkt aus argparse) ---
        location: Optional[str] = None,
        reference_seq_dir: Optional[str] = None,
        new_project: Optional[bool] = None,
        # Platz für zukünftige Felder
        **kwargs,
    ):
        # --- feste Basiswerte / bisherige HMSSS-Felder ---
        self.execute_location = execute_location
        self.project_name = project_name
        self.index_db = index_db
        self.csb_name_prefix = csb_name_prefix
        self.csb_name_suffix = csb_name_suffix
        self.genomeID_divider = genome_id_divider

        # Laufzeit-/Queue-Container
        self.finished_genomes = {} if finished_genomes is None else finished_genomes
        self.queued_genomes = {} if queued_genomes is None else queued_genomes
        self.fna_files = {} if faa_files is None else faa_files
        self.faa_files = {} if faa_files is None else faa_files
        self.gff_files = {} if gff_files is None else gff_files
        self.hmmreport_files = {} if hmmreport_files is None else hmmreport_files

        # Redundanz
        self.redundant = redundant
        self.non_redundant = non_redundant
        self.redundancy_hash = {} if redundancy_hash is None else redundancy_hash

        # Pipeline-Flags
        self.limiter = limiter
        self.fetch = fetch
        self.process = process

        # --- argparse: Input definition ---
        self.fasta_file_directory = fasta_file_directory
        self.score_threshold_file = score_threshold_file
        self.library = library
        self.result_files_directory = result_files_directory
        self.database_directory = database_directory
        self.cores = cores
        self.glob_report = glob_report
        self.verbose = verbose

        # --- argparse: Search parameters ---
        self.threshold_type = threshold_type
        self.thrs_score = thrs_score
        self.taxonomy_file = taxonomy_file
        self.refseq_identity = refseq_identity
        self.name = name
        self.stage = start_stage
        self.exit = exit_stage

        # --- argparse: Search library resources ---
        self.HMM_sets = [] if hmm_sets is None else hmm_sets
        self.clean_reports = clean_reports
        self.individual_reports = individual_reports
        self.max_seqs_per_genome = max_seqs_per_genome
        self.bool_cross_check = bool_cross_check
        self.optimized_cutoff_cross_check = optimized_cutoff_cross_check

        # --- argparse: Synteny options ---
        self.patterns_file = patterns_file
        self.cooccurrence_file = cooccurrence_file
        self.exclusion_singletons = exclusion_singletons
        self.min_completeness = min_completeness
        self.glob_chunks = glob_chunks

        # --- argparse: Information flags ---
        self.stat_keywords = stat_keywords
        self.stat_csb = stat_csb
        self.stat_genomes = stat_genomes

        # --- argparse: CSB prediction ---
        self.nucleotide_range = nucleotide_range
        self.insertions = insertions
        self.occurence = occurence
        self.min_csb_size = min_csb_size
        self.max_csb_size = max_csb_size
        self.max_domain_repeats = max_domain_repeats
        self.jaccard = jaccard

        # --- argparse: Work step regulation ---
        self.redo_taxonomy = redo_taxonomy

        # --- argparse: Limiter / Operators ---
        self.dataset_limit_lineage = dataset_limit_lineage
        self.dataset_limit_taxon = dataset_limit_taxon
        self.dataset_limit_proteins = dataset_limit_proteins
        self.dataset_limit_keywords = dataset_limit_keywords
        self.dataset_divide_sign = dataset_divide_sign
        self.dataset_limit_min_cluster_completeness = (
            dataset_limit_min_cluster_completeness
        )

        self.fetch_genomes = [] if fetch_genomes is None else fetch_genomes
        self.fetch_proteins = [] if fetch_proteins is None else fetch_proteins
        self.fetch_csbs = [] if fetch_csbs is None else fetch_csbs
        self.fetch_keywords = [] if fetch_keywords is None else fetch_keywords
        self.keywords_connector = keywords_connector

        # --- argparse: Processing ---
        self.merge_fasta = merge_fasta
        self.filter_fasta = filter_fasta
        self.concat_alignment = concat_alignment
        self.add_taxonomy = add_taxonomy
        self.add_genomic_context = add_genomic_context
        self.create_type_range_dataset = create_type_range_dataset
        self.create_gene_cluster_dataset = create_gene_cluster_dataset
        self.gaps = gaps

        # --- weitere intern genutzte Felder (werden später oft befüllt) ---
        # Pfade/Orte
        self.location = __location__ if location is None else location
        self.reference_seq_dir = (
            (__location__ + "/src/RefSeqs")
            if reference_seq_dir is None
            else reference_seq_dir
        )
        self.new_project = new_project  # wird in parse_arguments auf True/False gesetzt

        # Dateien/Verzeichnisse, die im Lauf erzeugt werden
        self.Cross_check_directory: Optional[str] = getattr(
            self, "Cross_check_directory", None
        )
        self.glob_trusted_hitreport: Optional[str] = getattr(
            self, "glob_trusted_hitreport", None
        )
        self.glob_intermediate_hitreport: Optional[str] = getattr(
            self, "glob_intermediate_hitreport", None
        )
        self.csb_output_file: Optional[str] = getattr(self, "csb_output_file", None)

        # Beliebige zusätzliche Schlüssel zulassen (Vorwärtskompatibilität)
        for k, v in kwargs.items():
            setattr(self, k, v)


__all__ = ["Hmsss", "__location__"]

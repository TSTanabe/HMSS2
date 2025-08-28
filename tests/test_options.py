import types


def test_hmsss_defaults():
    from hmsss.core.options import Hmsss, __location__

    o = Hmsss()

    # Basale Defaults
    assert o.execute_location == __location__
    assert o.project_name == "project"
    assert o.index_db is False
    assert o.csb_name_prefix == "csb-"
    assert o.csb_name_suffix == "_"
    assert o.genomeID_divider == "___"

    # Container sind vorhanden (keine shared mutable defaults)
    assert isinstance(o.finished_genomes, dict) and o.finished_genomes == {}
    assert isinstance(o.queued_genomes, dict) and o.queued_genomes == {}
    assert isinstance(o.faa_files, dict) and o.faa_files == {}
    assert isinstance(o.gff_files, dict) and o.gff_files == {}
    assert isinstance(o.hmmreport_files, dict) and o.hmmreport_files == {}
    assert isinstance(o.redundancy_hash, dict) and o.redundancy_hash == {}

    # Pipeline-Flags
    assert o.fetch is False
    assert o.process is False
    assert o.limiter is False

    # Wichtige Parser-Defaults
    assert o.result_files_directory == __location__ + "/results"
    assert o.database_directory is None
    assert o.cores == 4
    assert o.verbose == 1
    assert o.threshold_type == 1
    assert o.thrs_score == 50
    assert o.refseq_identity == 90
    assert o.name == "project"
    assert o.stage == 0 and o.exit == 10

    # Search-Library-Parameter
    assert o.hmm_sets == []
    assert o.clean_reports is False
    assert o.individual_reports is True
    assert o.max_seqs_per_genome == 4
    assert o.bool_cross_check is True
    assert o.optimized_cutoff_cross_check is False

    # Synteny / Chunks
    assert o.patterns_file is None
    assert o.cooccurrence_file is None
    assert o.exclusion_singletons is None
    assert o.min_completeness == 0.5
    assert o.glob_chunks == 5000

    # Info-Flags
    assert o.stat_keywords is False
    assert o.stat_csb is False
    assert o.stat_genomes is False

    # CSB prediction
    assert o.nucleotide_range == 3500
    assert o.insertions == 1
    assert o.occurence == 1
    assert o.min_csb_size == 4
    assert o.max_csb_size == 50
    assert o.max_domain_repeats == 4
    assert o.jaccard == 0.0

    # Limiter/Operators
    assert o.dataset_limit_lineage is None
    assert o.dataset_limit_taxon is None
    assert o.dataset_limit_proteins == "0"
    assert o.dataset_limit_keywords == "0"
    assert o.dataset_divide_sign == "."
    assert hasattr(o, "dataset_limit_min_cluster_completeness")

    assert o.fetch_genomes == []
    assert o.fetch_proteins == []
    assert o.fetch_csbs == []
    assert o.fetch_keywords == []
    assert o.keywords_connector == "OR"

    # Processing
    assert o.merge_fasta is None
    assert o.filter_fasta is None
    assert o.concat_alignment is None
    assert o.add_taxonomy is None
    assert o.add_genomic_context is None
    assert o.create_type_range_dataset is None
    assert o.create_gene_cluster_dataset is None
    assert o.gaps is False

    # Interne Felder / dynamische Felder
    assert isinstance(o.location, str) and o.location == __location__
    assert isinstance(o.reference_seq_dir, str)
    assert hasattr(o, "new_project")

    # dynamisch gesetzte spätere Pfade existieren als Attribute
    assert hasattr(o, "Cross_check_directory") and o.cross_check_directory is None
    assert hasattr(o, "glob_trusted_hitreport") and o.glob_trusted_hitreport is None
    assert (
        hasattr(o, "glob_intermediate_hitreport")
        and o.glob_intermediate_hitreport is None
    )
    assert hasattr(o, "csb_output_file")


def test_hmsss_overrides_and_types(tmp_path):
    from hmsss.core.options import Hmsss

    out_dir = tmp_path / "results"
    db_file = out_dir / "db.sqlite"
    o = Hmsss(
        result_files_directory=str(out_dir),
        database_directory=str(db_file),
        cores=8,
        fetch=True,
        process=True,
        limiter=True,
        HMM_sets=["setA", "setB"],
        fetch_genomes=["g1"],
        fetch_keywords=["k1", "k2"],
        dataset_limit_min_cluster_completeness=0.7,
    )

    assert o.result_files_directory == str(out_dir)
    assert o.database_directory == str(db_file)
    assert o.cores == 8
    assert o.fetch is True and o.process is True and o.limiter is True
    assert o.hmm_sets == ["setA", "setB"]
    assert o.fetch_genomes == ["g1"]
    assert o.fetch_keywords == ["k1", "k2"]
    assert o.dataset_limit_min_cluster_completeness == 0.7


def test_mutable_defaults_are_not_shared():
    from hmsss.core.options import Hmsss

    a = Hmsss()
    b = Hmsss()
    a.hmm_sets.append("X")
    a.fetch_genomes.append("A")
    a.redundancy_hash["k"] = 1

    assert b.hmm_sets == []  # nicht geteilt
    assert b.fetch_genomes == []  # nicht geteilt
    assert "k" not in b.redundancy_hash  # nicht geteilt


def test_kwargs_passthrough():
    from hmsss.core.options import Hmsss

    o = Hmsss(custom_field_foo=123, another="bar")
    assert getattr(o, "custom_field_foo") == 123
    assert getattr(o, "another") == "bar"

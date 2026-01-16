#!/usr/bin/env python3

import os
import tempfile
import shutil
import time

from graftm.sequence_search_results import SequenceSearchResult
from graftm.graftm_output_paths import GraftMFiles
from graftm.search_table import SearchTableWriter
from graftm.sequence_searcher import SequenceSearcher
from graftm.hmmsearcher import NoInputSequencesException
from graftm.housekeeping import HouseKeeping
from graftm.summarise import Stats_And_Summary
from graftm.pplacer import Pplacer

# from graftm.update import Update
from graftm.unpack_sequences import UnpackRawReads
from graftm.graftm_package import GraftMPackage
from graftm.expand_searcher import ExpandSearcher
from graftm.diamond import Diamond
from graftm.getaxnseq import Getaxnseq
from graftm.sequence_io import SequenceIO
from graftm.timeit import Timer
from graftm.clusterer import Clusterer

# from graftm.decorator import Decorator
# from graftm.external_program_suite import ExternalProgramSuite
# from graftm.archive import Archive
from graftm.decoy_filter import DecoyFilter
from hmsss.core.logging import get_logger

logging = get_logger(__name__)
T = Timer()


class UnrecognisedSuffixError(Exception):
    pass


class Run:
    PIPELINE_AA = "P"
    PIPELINE_NT = "D"

    _MIN_VERBOSITY_FOR_ART = 3  # with 2 then, only errors are printed

    PPLACER_TAXONOMIC_ASSIGNMENT = "pplacer"
    DIAMOND_TAXONOMIC_ASSIGNMENT = "diamond"

    MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES = 95
    MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES = 30

    DEFAULT_MAX_SAMPLES_FOR_KRONA = 100

    NO_ORFS_EXITSTATUS = 128

    def __init__(self, args):
        self.args = args
        self.setattributes(self.args)

    def setattributes(self, args):
        self.hk = HouseKeeping()
        self.s = Stats_And_Summary()

        self.hk.set_attributes(self.args)
        self.hk.set_euk_hmm(self.args)
        if args.euk_check:
            self.args.search_hmm_files.append(self.args.euk_hmm_file)

        self.ss = SequenceSearcher(
            self.args.search_hmm_files,
            (None if self.args.search_only else self.args.aln_hmm_file),
        )
        self.sequence_pair_list = self.hk.parameter_checks(args)
        if hasattr(args, "reference_package"):
            self.p = Pplacer(self.args.reference_package)

    def output_filepaths(self, base_list):
        """
        summarise - write summary information to file, including otu table, biom
                    file, krona plot, and timing information

        Parameters
        ----------
        base_list : array
            list of each of the files processed by graftm, with the path and
            and suffixed removed
        trusted_placements : dict
            dictionary of placements with entry as the key, a taxonomy string
            as the value
        reverse_pipe : bool
            True = run reverse pipe, False = run normal pipeline
        """

        # Summary steps.
        placements_list = []
        filepaths = []  # filepaths to taxonomy, alignment and sequence files
        for base in base_list:
            # First assign the hash that contains all of the trusted placements
            # to a variable to it can be passed to otu_builder, to be written
            # to a file. :)
            taxonomy_file = GraftMFiles(
                base, self.args.output_directory, False
            ).read_tax_output_path(base)

            alignment_file = GraftMFiles(
                base, self.args.output_directory, False
            ).aligned_fasta_output_path(base)

            sequence_file = GraftMFiles(
                base, self.args.output_directory, False
            ).fa_output_path(base)

            filepaths.append(
                {
                    "base": base,
                    "taxonomy": taxonomy_file,
                    "alignment": alignment_file,
                    "sequences": sequence_file,
                }
            )

        return filepaths

    def summarise(self, base_list, trusted_placements, reverse_pipe):
        """
        summarise - write summary information to file, including otu table, biom
                    file, krona plot, and timing information

        Parameters
        ----------
        base_list : array
            list of each of the files processed by graftm, with the path and
            and suffixed removed
        trusted_placements : dict
            dictionary of placements with entry as the key, a taxonomy string
            as the value
        reverse_pipe : bool
            True = run reverse pipe, False = run normal pipeline
        times : array
            list of the recorded times for each step in the pipeline in the
            format: [search_step_time, alignment_step_time, placement_step_time]
        hit_read_count_list : array
            list containing sublists, one for each file run through the GraftM
            pipeline, each two entries, the first being the number of putative
            eukaryotic reads (when searching 16S), the second being the number
            of hits aligned and placed in the tree.
        max_samples_for_krona: int
            If the number of files processed is greater than this number, then
            do not generate a krona diagram.
        Returns
        -------
        """

        # Summary steps.
        placements_list = []
        for base in base_list:
            # First assign the hash that contains all of the trusted placements
            # to a variable to it can be passed to otu_builder, to be written
            # to a file. :)
            placements = trusted_placements[base]
            self.s.readTax(
                placements,
                GraftMFiles(
                    base, self.args.output_directory, False
                ).read_tax_output_path(base),
            )
            placements_list.append(placements)

        # Generate coverage table
        # logging.info('Building coverage table for %s' % base)
        # self.s.coverage_of_hmm(self.args.aln_hmm_file,
        #                         self.gmf.summary_table_output_path(base),
        #                         self.gmf.coverage_table_path(base),
        #                         summary_dict[base]['read_length'])

        with open(self.gmf.combined_summary_table_output_path(), "w") as f:
            self.s.write_tabular_otu_table(base_list, placements_list, f)

        # logging.info('Writing biom file')
        # with biom_open(self.gmf.combined_biom_output_path(), 'w') as f:
        #    biom_successful = self.s.write_biom(base_list, placements_list, f)
        # if not biom_successful:
        #    os.remove(self.gmf.combined_biom_output_path())

        #
        # --- Skip Krona Output ---
        #
        # logging.info('Building summary krona plot')
        # if len(base_list) > max_samples_for_krona:
        #    logging.warning(
        #        "Skipping creation of Krona diagram since there are too many input files. The maximum can be overridden using --max_samples_for_krona")
        # else:
        #    self.s.write_krona_plot(base_list, placements_list, self.gmf.krona_output_path())

        # Basic statistics
        # placed_reads = [len(trusted_placements[base]) for base in base_list]
        # self.s.build_basic_statistics(times, hit_read_count_list, placed_reads, \
        #                              base_list, self.gmf.basic_stats_path())

        # Delete unnecessary files
        for base in base_list:
            directions = ["forward", "reverse"]
            if reverse_pipe:
                for i in range(0, 2):
                    self.gmf = GraftMFiles(
                        base, self.args.output_directory, directions[i]
                    )
                    self.hk.delete(
                        [
                            self.gmf.for_aln_path(base),
                            self.gmf.rev_aln_path(base),
                            self.gmf.conv_output_rev_path(base),
                            self.gmf.conv_output_for_path(base),
                            self.gmf.euk_free_path(base),
                            self.gmf.euk_contam_path(base),
                            self.gmf.readnames_output_path(base),
                            self.gmf.sto_output_path(base),
                            self.gmf.orf_titles_output_path(base),
                            self.gmf.orf_output_path(base),
                            self.gmf.output_for_path(base),
                            self.gmf.output_rev_path(base),
                        ]
                    )
            else:
                self.gmf = GraftMFiles(base, self.args.output_directory, False)
                self.hk.delete(
                    [
                        self.gmf.for_aln_path(base),
                        self.gmf.rev_aln_path(base),
                        self.gmf.conv_output_rev_path(base),
                        self.gmf.conv_output_for_path(base),
                        self.gmf.euk_free_path(base),
                        self.gmf.euk_contam_path(base),
                        self.gmf.readnames_output_path(base),
                        self.gmf.sto_output_path(base),
                        self.gmf.orf_titles_output_path(base),
                        self.gmf.orf_output_path(base),
                        self.gmf.output_for_path(base),
                        self.gmf.output_rev_path(base),
                    ]
                )

        # logging.info('Done, thanks for using graftM!\n')

    def graft(self):
        # The Graft pipeline:
        # Searches for reads using hmmer, and places them in phylogenetic
        # trees to derive a community structure.
        if self.args.graftm_package:
            gpkg = GraftMPackage.acquire(self.args.graftm_package)
        else:
            gpkg = None

        REVERSE_PIPE = True if self.args.reverse else False
        INTERLEAVED = True if self.args.interleaved else False
        base_list = []
        seqs_list = []
        search_results = []
        hit_read_count_list = []
        db_search_results = []

        if gpkg:
            maximum_range = gpkg.maximum_range()

            if self.args.search_diamond_file:
                self.args.search_method = self.hk.DIAMOND_SEARCH_METHOD
                diamond_db = (
                    self.args.search_diamond_file
                )  # this var is a string. Remove the [0] to fix this bug
            else:
                diamond_db = gpkg.diamond_database_path()
                if self.args.search_method == self.hk.DIAMOND_SEARCH_METHOD:
                    if not diamond_db:
                        logging.error(
                            "%s search method selected, but no diamond database specified. \
                        Please either provide a gpkg to the --graftm_package flag, or a diamond \
                        database to the --search_diamond_file flag."
                            % self.args.search_method
                        )
                        raise Exception()
        else:
            logging.error(
                "%s search method selected, but no gpkg provided "
                % self.args.search_method
            )

        if self.args.assignment_method == Run.DIAMOND_TAXONOMIC_ASSIGNMENT:
            if self.args.reverse:
                logging.warning(
                    "--reverse reads specified with --assignment_method diamond. Reverse reads will be ignored."
                )
                self.args.reverse = None

        # If merge reads is specified, check that there are reverse reads to merge with
        if self.args.merge_reads and not hasattr(self.args, "reverse"):
            raise Exception(
                f"Reverse reads are missing for file {hasattr(self.args, 'forward')}."
            )

        # Set the output directory if not specified and create that directory
        logging.debug("Creating working directory: %s" % self.args.output_directory)
        self.hk.make_working_directory(self.args.output_directory, self.args.force)

        # Set pipeline and evalue by checking HMM format
        if self.args.search_only:
            if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
                hmm_type, hmm_tc = self.hk.setpipe(self.args.search_hmm_files[0])
                logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))
        else:
            hmm_type, hmm_tc = self.hk.setpipe(self.args.aln_hmm_file)
            logging.debug("HMM type: %s Trusted Cutoff: %s" % (hmm_type, hmm_tc))

        if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
            setattr(self.args, "type", hmm_type)
            if hmm_tc:
                setattr(self.args, "evalue", "--cut_tc")
        else:
            setattr(self.args, "type", self.PIPELINE_AA)

        if self.args.filter_minimum is not None:
            filter_minimum = self.args.filter_minimum
        else:
            if self.args.type == self.PIPELINE_NT:
                filter_minimum = Run.MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES
            else:
                filter_minimum = Run.MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES

        # Generate expand_search database if required
        if self.args.expand_search_contigs:
            if self.args.graftm_package:
                pkg = GraftMPackage.acquire(self.args.graftm_package)
            else:
                pkg = None
            boots = ExpandSearcher(
                search_hmm_files=self.args.search_hmm_files,
                maximum_range=self.args.maximum_range,
                threads=self.args.threads,
                evalue=self.args.evalue,
                min_orf_length=self.args.min_orf_length,
                graftm_package=pkg,
            )

            # this is a hack, it should really use GraftMFiles but that class isn't currently flexible enough
            new_database = (
                os.path.join(self.args.output_directory, "expand_search.hmm")
                if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD
                else os.path.join(self.args.output_directory, "expand_search")
            )

            if boots.generate_expand_search_database_from_contigs(
                self.args.expand_search_contigs, new_database, self.args.search_method
            ):
                if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
                    self.ss.search_hmm.append(new_database)
                else:
                    diamond_db = new_database

        first_search_method = self.args.search_method
        if self.args.decoy_database:
            decoy_filter = DecoyFilter(
                Diamond(diamond_db, threads=self.args.threads),
                Diamond(self.args.decoy_database, threads=self.args.threads),
            )
            doing_decoy_search = True
        elif self.args.search_method == self.hk.HMMSEARCH_AND_DIAMOND_SEARCH_METHOD:
            decoy_filter = DecoyFilter(Diamond(diamond_db, threads=self.args.threads))
            doing_decoy_search = True
            first_search_method = self.hk.HMMSEARCH_SEARCH_METHOD
        else:
            doing_decoy_search = False

        # For each pair (or single file passed to GraftM)
        logging.debug("Working with %i file(s)" % len(self.sequence_pair_list))
        for pair in self.sequence_pair_list:
            # Guess the sequence file type, if not already specified to GraftM
            unpack = UnpackRawReads(pair[0], self.args.input_sequence_type, INTERLEAVED)

            # Set the basename, and make an entry to the summary table.
            base = unpack.basename()
            pair_direction = ["forward", "reverse"]
            logging.info("Working on %s" % base)

            # Make the working base subdirectory
            self.hk.make_working_directory(
                os.path.join(self.args.output_directory, base), self.args.force
            )

            # for each of the paired end read files
            for read_file in pair:
                unpack = UnpackRawReads(
                    read_file, self.args.input_sequence_type, INTERLEAVED
                )
                if read_file is None:
                    # placeholder for interleaved (second file is None)
                    continue

                if not os.path.isfile(read_file):  # Check file exists
                    logging.info("%s does not exist! Skipping this file.." % read_file)
                    continue

                # Set the output file_name
                if len(pair) == 2:
                    direction = (
                        "interleaved" if pair[1] is None else pair_direction.pop(0)
                    )
                    logging.info("Working on %s reads" % direction)
                    self.gmf = GraftMFiles(base, self.args.output_directory, direction)
                    self.hk.make_working_directory(
                        os.path.join(self.args.output_directory, base, direction),
                        self.args.force,
                    )
                else:
                    direction = False
                    self.gmf = GraftMFiles(base, self.args.output_directory, direction)

                t0 = time.perf_counter()
                if self.args.type == self.PIPELINE_AA:
                    logging.debug("Running protein pipeline")
                    try:
                        search_time, (result, complement_information) = (
                            self.ss.aa_db_search(
                                self.gmf,
                                base,
                                unpack,
                                first_search_method,
                                maximum_range,
                                self.args.threads,
                                self.args.evalue,
                                self.args.min_orf_length,
                                self.args.restrict_read_length,
                                self.args.translation_table,
                                diamond_db,
                                self.args.diamond_performance_parameters,
                            )
                        )
                    except NoInputSequencesException as e:
                        logging.error(
                            "No sufficiently long open reading frames were found, indicating"
                            " either the input sequences are too short or the min orf length"
                            " cutoff is too high. Cannot continue sorry. Alternatively, there"
                            " is something amiss with the installation of OrfM. The specific"
                            " command that failed was: %s" % e.command
                        )
                        exit(Run.NO_ORFS_EXITSTATUS)
                    dt = time.perf_counter() - t0
                    logging.debug(f"[TIME] aa_db_search took {dt:.3f} seconds")

                # Or the DNA pipeline
                elif self.args.type == self.PIPELINE_NT:
                    logging.debug("Running nucleotide pipeline")
                    search_time, (result, complement_information) = (
                        self.ss.nt_db_search(
                            self.gmf,
                            base,
                            unpack,
                            self.args.euk_check,
                            self.args.search_method,
                            maximum_range,
                            self.args.threads,
                            self.args.evalue,
                        )
                    )

                reads_detected = True
                if not result.hit_fasta() or os.path.getsize(result.hit_fasta()) == 0:
                    logging.info("No reads found in %s" % base)
                    reads_detected = False

                if self.args.search_only:
                    db_search_results.append(result)
                    base_list.append(base)
                    continue

                t0 = time.perf_counter()
                # Filter out decoys if specified
                if reads_detected and doing_decoy_search:
                    with tempfile.NamedTemporaryFile(
                        prefix="graftm_decoy", suffix=".fa"
                    ) as f:
                        tmpname = f.name
                    any_remaining = decoy_filter.filter(result.hit_fasta(), tmpname)
                    if any_remaining:
                        shutil.move(tmpname, result.hit_fasta())
                    else:
                        # No hits remain after decoy filtering.
                        os.remove(result.hit_fasta())
                        continue
                dt = time.perf_counter() - t0
                logging.debug(f"[TIME] decoy filtering took {dt:.3f} seconds")

                t0 = time.perf_counter()
                if self.args.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
                    logging.info("aligning reads to reference package database")
                    hit_aligned_reads = self.gmf.aligned_fasta_output_path(base)

                    if reads_detected:
                        aln_time, aln_result = self.ss.align(
                            result.hit_fasta(),
                            hit_aligned_reads,
                            complement_information,
                            self.args.type,
                            filter_minimum,
                        )
                    else:
                        aln_time = "n/a"
                    if not os.path.exists(
                        hit_aligned_reads
                    ):  # If all were filtered out, or there just was none..
                        with open(hit_aligned_reads, "w") as f:
                            pass  # just touch the file, nothing else
                    seqs_list.append(hit_aligned_reads)
                dt = time.perf_counter() - t0
                logging.debug(f"[TIME] pplacer preparation took {dt:.3f} seconds")
                db_search_results.append(result)
                base_list.append(base)
                search_results.append(result.search_result)
                hit_read_count_list.append(result.hit_count)

        # Write summary table
        srchtw = SearchTableWriter()
        srchtw.build_search_otu_table(
            [x.search_objects for x in db_search_results],
            base_list,
            self.gmf.search_otu_table(),
        )

        if self.args.search_only:
            logging.info("Stopping before alignment and taxonomic assignment phase\n")
            exit(0)

        if (
            self.args.merge_reads
        ):  # not run when diamond is the assignment mode- enforced by argparse grokking
            logging.debug("Running merge reads output")
            if self.args.interleaved:
                fwd_seqs = seqs_list
                rev_seqs = []
            else:
                base_list = base_list[0::2]
                fwd_seqs = seqs_list[0::2]
                rev_seqs = seqs_list[1::2]
            merged_output = [
                GraftMFiles(
                    base, self.args.output_directory, False
                ).aligned_fasta_output_path(base)
                for base in base_list
            ]
            logging.debug("merged reads to %s", merged_output)
            self.ss.merge_forev_aln(fwd_seqs, rev_seqs, merged_output)
            seqs_list = merged_output
            REVERSE_PIPE = False

        elif REVERSE_PIPE:
            base_list = base_list[0::2]

        t0 = time.perf_counter()
        # Leave the pipeline if search only was specified
        if self.args.search_and_align_only:
            logging.info("Stopping before taxonomic assignment phase\n")
            exit(0)
        elif not any(base_list):
            logging.info(
                "No hits in any of the provided files. Cannot continue with no reads to assign taxonomy to.\n"
            )
            exit(0)
        self.gmf = GraftMFiles("", self.args.output_directory, False)

        if self.args.assignment_method == Run.PPLACER_TAXONOMIC_ASSIGNMENT:
            clusterer = Clusterer()
            # Classification steps
            seqs_list = clusterer.cluster(seqs_list, REVERSE_PIPE)
            logging.info("Placing reads into phylogenetic tree")
            taxonomic_assignment_time, assignments = self.p.place(
                REVERSE_PIPE,
                seqs_list,
                self.args.resolve_placements,
                self.gmf,
                self.args,
                result.slash_endings,
                gpkg.taxtastic_taxonomy_path(),
                clusterer,
            )
            assignments = clusterer.uncluster_annotations(assignments, REVERSE_PIPE)

        elif self.args.assignment_method == Run.DIAMOND_TAXONOMIC_ASSIGNMENT:
            logging.info("Assigning taxonomy with diamond")
            taxonomic_assignment_time, assignments = self._assign_taxonomy_with_diamond(
                base_list,
                db_search_results,
                gpkg,
                self.gmf,
                self.args.diamond_performance_parameters,
            )
            aln_time = "n/a"
        else:
            raise Exception(
                "Unexpected assignment method encountered: %s"
                % self.args.placement_method
            )

        dt = time.perf_counter() - t0
        logging.debug(f"[TIME] pplacer assignment phase took {dt:.3f} seconds")
        # Prepare read mapping und alignments for return
        read_tax_paths = {}
        alignment_paths = {}

        for base in base_list:
            gmf = GraftMFiles(base, self.args.output_directory, False)
            read_tax_paths[base] = gmf.read_tax_output_path(base)

            # alignments (abhängig von direction)
            aln_files = []

            if REVERSE_PIPE:
                for direction in ("forward", "reverse"):
                    gmf_dir = GraftMFiles(base, self.args.output_directory, direction)
                    aln_files.append(gmf_dir.aligned_fasta_output_path(base))
            else:
                aln_files.append(gmf.aligned_fasta_output_path(base))

            alignment_paths[base] = aln_files

        # Remove temp files and create read mapping summary
        self.summarise(
            base_list=base_list,
            trusted_placements=assignments,
            reverse_pipe=REVERSE_PIPE,
        )
        filepaths = self.output_filepaths(base_list=base_list)
        return filepaths  # filepaths to taxonomy, alignment and sequence files. Each field has a dict for the

    @T.timeit
    def _assign_taxonomy_with_diamond(
        self,
        base_list,
        db_search_results,
        graftm_package,
        graftm_files,
        diamond_performance_parameters,
    ):
        """Run diamond to assign taxonomy

        Parameters
        ----------
        base_list: list of str
            list of sequence block names
        db_search_results: list of DBSearchResult
            the result of running hmmsearches
        graftm_package: GraftMPackage object
            Diamond is run against this database
        graftm_files: GraftMFiles object
            Result files are written here
        diamond_performance_parameters : str
            extra args for DIAMOND

        Returns
        -------
        list of
        1. time taken for assignment
        2. assignments i.e. dict of base_list entry to dict of read names to
            to taxonomies, or None if there was no hit detected.
        """
        runner = Diamond(
            graftm_package.diamond_database_path(), self.args.threads, self.args.evalue
        )
        taxonomy_definition = Getaxnseq().read_taxtastic_taxonomy_and_seqinfo(
            open(graftm_package.taxtastic_taxonomy_path()),
            open(graftm_package.taxtastic_seqinfo_path()),
        )
        results = {}

        # For each of the search results,
        for i, search_result in enumerate(db_search_results):
            if search_result.hit_fasta() is None:
                sequence_id_to_taxonomy = {}
            else:
                sequence_id_to_hit = {}
                # Run diamond
                logging.debug("Running diamond on %s" % search_result.hit_fasta())
                diamond_result = runner.run(
                    search_result.hit_fasta(),
                    UnpackRawReads.PROTEIN_SEQUENCE_TYPE,
                    daa_file_basename=graftm_files.diamond_assignment_output_basename(
                        base_list[i]
                    ),
                    extra_args=diamond_performance_parameters,
                )
                for res in diamond_result.each(
                    [
                        SequenceSearchResult.QUERY_ID_FIELD,
                        SequenceSearchResult.HIT_ID_FIELD,
                    ]
                ):
                    if res[0] in sequence_id_to_hit:
                        # do not accept duplicates
                        if sequence_id_to_hit[res[0]] != res[1]:
                            raise Exception(
                                "Diamond unexpectedly gave two hits for a single query sequence for %s"
                                % res[0]
                            )
                    else:
                        sequence_id_to_hit[res[0]] = res[1]

                # Extract taxonomy of the best hit, and add in the no hits
                sequence_id_to_taxonomy = {}
                for seqio in SequenceIO().read_fasta_file(search_result.hit_fasta()):
                    name = seqio.name
                    if name in sequence_id_to_hit:
                        # Add Root; to be in line with pplacer assignment method
                        sequence_id_to_taxonomy[name] = ["Root"] + taxonomy_definition[
                            sequence_id_to_hit[name]
                        ]
                    else:
                        # picked up in the initial search (by hmmsearch, say), but diamond misses it
                        sequence_id_to_taxonomy[name] = ["Root"]

            results[base_list[i]] = sequence_id_to_taxonomy
        return results

    def main(self):
        if self.args.subparser_name == "graft":
            return self.graft()
        else:
            raise Exception(
                "Unexpected graftM subparser name %s" % self.args.subparser_name
            )

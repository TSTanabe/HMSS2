#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
#
# graftM - A pipeline for gene centric analyses of metagenome datasets
#
###############################################################################

__author__ = "Joel Boyd, Ben Woodcroft"
__copyright__ = "Copyright 2014"
__credits__ = ["Joel Boyd", "Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd, Ben Woodcroft"
__email__ = "joel.boyd near uq.net.au, b.woodcroft near uq.edu.au"
__status__ = "Development"

# This is the 'run' module from the graftM package, adjusted to am
# minimal functionality. Only read detection and placement for a single
# metagnome is possible and output was limited to the placement of reads
# Other outputs and functionalities of the package, except for the graftM graft
# were removed

import os
from hmsss.core.logging import get_logger, print_header
logger = get_logger(__name__)

import tempfile
import shutil

from hmsss.graft.housekeeping import HouseKeeping
from hmsss.graft.summarise import Stats_And_Summary
from hmsss.graft.graftm_package import GraftMPackage
from hmsss.graft.sequence_search_results import SequenceSearchResult
from hmsss.graft.graftm_output_paths import GraftMFiles
from hmsss.graft.search_table import SearchTableWriter
from hmsss.graft.sequence_searcher import SequenceSearcher
from hmsss.graft.hmmsearcher import NoInputSequencesException


from hmsss.graft.pplacer import Pplacer
from hmsss.graft.unpack_sequences import UnpackRawReads
from hmsss.graft.expand_searcher import ExpandSearcher
from hmsss.graft.diamond import Diamond
from hmsss.graft.getaxnseq import Getaxnseq # package not required?
from hmsss.graft.sequence_io import SequenceIO
from hmsss.graft.clusterer import Clusterer
from hmsss.graft.external_program_suite import ExternalProgramSuite
from hmsss.graft.decoy_filter import DecoyFilter


class UnrecognisedSuffixError(Exception):
    pass


class Run:
    """
    Minimale GraftM-Run-Klasse, reduziert auf:

    - Subcommand 'graft'
    - Placement immer mit pplacer
    - Search immer hmmsearch+diamond
    - Nutzung eines GraftM-Pakets (gpkg) + Diamond-DB aus gpkg
    - Decoy-Datenbank (für Decoy-Filter)
    - Keine Krona-Ausgabe
    """

    PIPELINE_AA = "P"
    PIPELINE_NT = "D"

    _MIN_VERBOSITY_FOR_ART = 3  # ab welcher verbosity das ASCII-Art gedruckt wird

    PPLACER_TAXONOMIC_ASSIGNMENT = 'pplacer'
    DIAMOND_TAXONOMIC_ASSIGNMENT = 'diamond'  # wird hier nicht aktiv benutzt, aber behalten für Kompatibilität

    MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES = 95
    MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES = 30

    NO_ORFS_EXITSTATUS = 128

    def __init__(self, args):
        self.args = args
        self.setattributes(self.args)

    def setattributes(self, args):
        """
        Initialisierung für die 'graft'-Pipeline.
        Andere Subcommands (create, update, tree, archive, expand_search) wurden entfernt.
        """

        self.hk = HouseKeeping()
        self.s = Stats_And_Summary()

        if args.subparser_name != 'graft':
            raise Exception(f"This class only supports subparser_name='graft', got: {args.subparser_name!r}")

        # Deine Vorgaben erzwingen:
        # - immer hmmsearch+diamond
        # - immer pplacer
        # - Pflicht: gpkg + decoy_database
        if not getattr(args, 'graftm_package', None):
            raise Exception("graftm_package is required for this Run implementation.")

        if not getattr(args, 'decoy_database', None):
            raise Exception("decoy_database is required for this Run implementation.")

        # Logging/externes Tool-Setup (ohne ktImportText, da keine Krona-Ausgabe)
        commands = ExternalProgramSuite([
            'orfm', 'nhmmer', 'hmmsearch',
            'mfqe', 'pplacer',
            'diamond',
        ])

        # Graft-spezifische Attribute setzen
        self.hk.set_attributes(self.args)
        self.hk.set_euk_hmm(self.args)
        if args.euk_check:
            self.args.search_hmm_files.append(self.args.euk_hmm_file)

        # SequenceSearcher-Objekt
        self.ss = SequenceSearcher(
            self.args.search_hmm_files,
            (None if self.args.search_only else self.args.aln_hmm_file),
        )

        # Eingabedateien prüfen und Sequenz-Paare ermitteln
        self.sequence_pair_list = self.hk.parameter_checks(args)

        if len(self.sequence_pair_list) != 1:
            raise Exception(
                f"This Run implementation expects exactly one metagenome/sample, "
                f"but got {len(self.sequence_pair_list)}."
            )

        # Pplacer vorbereiten, falls Referenzpaket aus gpkg bekannt
        if hasattr(args, 'reference_package'):
            self.p = Pplacer(self.args.reference_package)

    def summarise(self, base_list, trusted_placements, reverse_pipe, times,
                  hit_read_count_list):
        """
        Zusammenfassung für genau EIN Metagenom/Sample.

        Outputs im self.args.output_directory:
        - combined_count_table.txt  (Taxon → Counts)
        - read_tax.tsv              (Read → Taxonomie, optional mit Sample)
        """

        placements_list = []

        # 1) Gemeinsame read_tax.tsv im Output-Directory schreiben
        read_tax_path = os.path.join(self.args.output_directory, "read_tax.tsv")
        logger.info("Writing read→taxonomy assignments to %s", read_tax_path)

        with open(read_tax_path, "w") as out:
            # Header (kannst du weglassen, aber sehr praktisch für Downstream)
            out.write("sample\tread_id\ttaxonomy\n")

            for base in base_list:
                placements = trusted_placements[base]
                placements_list.append(placements)

                # Annahme: placements ist dict: read_id -> Sequenz von Taxon-Strings
                for read_id, tax in placements.items():
                    out.write(
                        "{}\t{}\t{}\n".format(
                            base,
                            read_id,
                            "; ".join(tax),
                        )
                    )

        # 2) combined_count_table.txt wie gehabt (eine Spalte = dein Sample)
        logger.info('Writing summary table')
        with open(self.gmf.combined_summary_table_output_path(), 'w') as f:
            self.s.write_tabular_otu_table(base_list, placements_list, f)

        # Delete unnecessary files
        logger.info('Cleaning up')
        for base in base_list:
            directions = ['forward', 'reverse']
            if reverse_pipe:
                for i in range(0, 2):
                    self.gmf = GraftMFiles(base, self.args.output_directory, directions[i])
                    self.hk.delete([
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
                    ])
            else:
                self.gmf = GraftMFiles(base, self.args.output_directory, False)
                self.hk.delete([
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
                ])

        logger.info('Done, thanks for using graftM!\n')

    def main(self):
        """
        Graft-Pipeline:

        - Reads suchen (hmmsearch + diamond)
        - ggf. expand_search_contigs
        - Decoy-Filtering
        - Alignment der Hits
        - pplacer-Placement
        - Summary & Cleanup
        """

        gpkg = GraftMPackage.acquire(self.args.graftm_package)

        REVERSE_PIPE = bool(self.args.reverse)
        INTERLEAVED = bool(self.args.interleaved)
        base_list = []
        seqs_list = []
        search_results = []
        hit_read_count_list = []
        db_search_results = []

        # Diamond-Datenbank aus gpkg
        maximum_range = gpkg.maximum_range()
        diamond_db = self.args.search_diamond_file # Prefer arugment diamond db
        if not diamond_db:
            # Fallback auf Diamond db aus gpkg
            diamond_db = gpkg.diamond_database_path()
            if not diamond_db:
                logger.error(
                    "No diamond database found in GraftM package %s, but this Run "
                    "implementation requires a diamond DB.",
                    self.args.graftm_package,
                )
                raise Exception("Missing diamond database in GraftM package")

        if self.args.assignment_method == self.DIAMOND_TAXONOMIC_ASSIGNMENT:
            # Sollte durch unsere setattributes() nie eintreten, aber sicherheitshalber:
            if self.args.reverse:
                logger.warning(
                    "--reverse reads specified with --assignment_method diamond. "
                    "Reverse reads will be ignored."
                )
                self.args.reverse = None

        if self.args.merge_reads and not hasattr(self.args, 'reverse'):
            # If merge reads is specified, check that there are reverse reads to merge with
            raise Exception("Programming error: merge_reads True but no reverse reads present")

        # Output-Directory anlegen
        logger.debug('Creating working directory: %s', self.args.output_directory)
        self.hk.make_working_directory(self.args.output_directory, self.args.force)

        # HMM-Typ bestimmen
        if self.args.search_only:
            if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
                hmm_type, hmm_tc = self.hk.setpipe(self.args.search_hmm_files[0])
                logger.debug("HMM type: %s Trusted Cutoff: %s", hmm_type, hmm_tc)
        else:
            hmm_type, hmm_tc = self.hk.setpipe(self.args.aln_hmm_file)
            logger.debug("HMM type: %s Trusted Cutoff: %s", hmm_type, hmm_tc)

        if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
            setattr(self.args, 'type', hmm_type)
            if hmm_tc:
                setattr(self.args, 'evalue', '--cut_tc')
        else:
            setattr(self.args, 'type', self.PIPELINE_AA)

        # Filter-Minimum setzen
        if self.args.filter_minimum is not None:
            filter_minimum = self.args.filter_minimum
        else:
            if self.args.type == self.PIPELINE_NT:
                filter_minimum = self.MIN_ALIGNED_FILTER_FOR_NUCLEOTIDE_PACKAGES
            else:
                filter_minimum = self.MIN_ALIGNED_FILTER_FOR_AMINO_ACID_PACKAGES

        # Expand_search (optional)
        if self.args.expand_search_contigs:
            pkg = gpkg
            boots = ExpandSearcher(
                search_hmm_files=self.args.search_hmm_files,
                maximum_range=self.args.maximum_range,
                threads=self.args.threads,
                evalue=self.args.evalue,
                min_orf_length=self.args.min_orf_length,
                graftm_package=pkg,
            )

            new_database = (
                os.path.join(self.args.output_directory, "expand_search.hmm")
                if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD
                else os.path.join(self.args.output_directory, "expand_search")
            )

            if boots.generate_expand_search_database_from_contigs(
                self.args.expand_search_contigs,
                new_database,
                self.args.search_method,
            ):
                if self.args.search_method == self.hk.HMMSEARCH_SEARCH_METHOD:
                    self.ss.search_hmm.append(new_database)
                else:
                    diamond_db = new_database

        # Decoy-Filter: für deinen Use Case immer aktiv (decoy_database Pflicht)
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

        logger.debug('Working with %i file(s)', len(self.sequence_pair_list))

        for pair in self.sequence_pair_list:
            # Dateityp raten falls nötig
            unpack = UnpackRawReads(
                pair[0],
                self.args.input_sequence_type,
                INTERLEAVED,
            )
            base = unpack.basename()
            pair_direction = ['forward', 'reverse']
            logger.info("Working on %s", base)

            # Subdirectory für Base
            self.hk.make_working_directory(
                os.path.join(self.args.output_directory, base),
                self.args.force,
            )

            # Durch Reads iterieren (forward/reverse oder interleaved)
            for read_file in pair:
                unpack = UnpackRawReads(
                    read_file,
                    self.args.input_sequence_type,
                    INTERLEAVED,
                )
                if read_file is None:
                    continue

                if not os.path.isfile(read_file):
                    logger.info('%s does not exist! Skipping this file..', read_file)
                    continue

                if len(pair) == 2:
                    direction = 'interleaved' if pair[1] is None else pair_direction.pop(0)
                    logger.info("Working on %s reads", direction)
                    self.gmf = GraftMFiles(base, self.args.output_directory, direction)
                    self.hk.make_working_directory(
                        os.path.join(self.args.output_directory, base, direction),
                        self.args.force,
                    )
                else:
                    direction = False
                    self.gmf = GraftMFiles(base, self.args.output_directory, direction)

                # Protein-Pipeline
                if self.args.type == self.PIPELINE_AA:
                    logger.debug("Running protein pipeline")
                    try:
                        search_time, (result, complement_information) = self.ss.aa_db_search(
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
                    except NoInputSequencesException as e:
                        logger.error(
                            "No sufficiently long open reading frames were found. "
                            "Cannot continue. Command that failed was: %s",
                            e.command,
                        )
                        exit(self.NO_ORFS_EXITSTATUS)

                # DNA-Pipeline
                elif self.args.type == self.PIPELINE_NT:
                    logger.debug("Running nucleotide pipeline")
                    search_time, (result, complement_information) = self.ss.nt_db_search(
                        self.gmf,
                        base,
                        unpack,
                        self.args.euk_check,
                        self.args.search_method,
                        maximum_range,
                        self.args.threads,
                        self.args.evalue,
                    )
                else:
                    raise Exception(f"Unexpected pipeline type: {self.args.type!r}")

                reads_detected = True
                if not result.hit_fasta() or os.path.getsize(result.hit_fasta()) == 0:
                    logger.info('No reads found in %s', base)
                    reads_detected = False

                if self.args.search_only:
                    db_search_results.append(result)
                    base_list.append(base)
                    continue

                # Decoy-Filter anwenden
                if reads_detected and doing_decoy_search:
                    with tempfile.NamedTemporaryFile(prefix="graftm_decoy", suffix='.fa') as f:
                        tmpname = f.name
                    any_remaining = decoy_filter.filter(
                        result.hit_fasta(),
                        tmpname,
                    )
                    if any_remaining:
                        shutil.move(tmpname, result.hit_fasta())
                    else:
                        os.remove(result.hit_fasta())
                        continue

                # pplacer-Placement vorbereiten
                if self.args.assignment_method == self.PPLACER_TAXONOMIC_ASSIGNMENT:
                    logger.info('aligning reads to reference package database')
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
                        aln_time = 'n/a'

                    if not os.path.exists(hit_aligned_reads):
                        with open(hit_aligned_reads, 'w'):
                            pass  # Datei „anfassen“
                    seqs_list.append(hit_aligned_reads)

                db_search_results.append(result)
                base_list.append(base)
                search_results.append(result.search_result)
                hit_read_count_list.append(result.hit_count)

        # Search-OTU-Tabelle
        srchtw = SearchTableWriter()
        srchtw.build_search_otu_table(
            [x.search_objects for x in db_search_results],
            base_list,
            self.gmf.search_otu_table(),
        )

        if self.args.search_only:
            logger.info('Stopping before alignment and taxonomic assignment phase\n')
            exit(0)

        # Merge-Reads falls gewünscht
        if self.args.merge_reads:
            logger.debug("Running merge reads output")
            if self.args.interleaved:
                fwd_seqs = seqs_list
                rev_seqs = []
            else:
                base_list = base_list[0::2]
                fwd_seqs = seqs_list[0::2]
                rev_seqs = seqs_list[1::2]
            merged_output = [
                GraftMFiles(base, self.args.output_directory, False).aligned_fasta_output_path(base)
                for base in base_list
            ]
            logger.debug("merged reads to %s", merged_output)
            self.ss.merge_forev_aln(fwd_seqs, rev_seqs, merged_output)
            seqs_list = merged_output
            REVERSE_PIPE = False

        elif REVERSE_PIPE:
            base_list = base_list[0::2]

        if self.args.search_and_align_only:
            logger.info('Stopping before taxonomic assignment phase\n')
            exit(0)
        elif not any(base_list):
            logger.error(
                'No hits in any of the provided files. Cannot continue with no reads to assign taxonomy to.\n'
            )
            exit(0)

        self.gmf = GraftMFiles('', self.args.output_directory, False)

        # pplacer-Assignment (immer, gemäß deiner Vorgabe)
        clusterer = Clusterer()
        seqs_list = clusterer.cluster(seqs_list, REVERSE_PIPE)
        logger.info("Placing reads into phylogenetic tree")
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

        # Zusammenfassung (ohne Krona)
        self.summarise(
            base_list,
            assignments,
            REVERSE_PIPE,
            [search_time, aln_time, taxonomic_assignment_time],
            hit_read_count_list,
        )
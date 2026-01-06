#!/usr/bin/env python3

from graftm.external_program_suite import ExternalProgramSuite

from hmsss.db import database
from hmsss.graft import read_counter, graft_mp
from hmsss.graft import generate_task
from hmsss.db.database import create_database
from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def _check_dependencies():
    commands = ExternalProgramSuite(
        [
            "orfm",
            "nhmmer",
            "hmmsearch",
            "mfqe",
            "pplacer",
            "diamond",
        ]
    )


def initial_read_mapping(config):
    # _check_dependencies()  # Check if graftM dependencies are present
    create_database(config.database_directory)
    task_list = generate_task.initialize_task_list(config)

    # TODO Generate the database for the mapped data
    logger.info("Counting reads for all (meta-)genomes")
    meta_dict, genome_id_set = read_counter.collect_metagenome_counts_parallel(
        task_list, processes=4, chunksize=4
    )
    # 1) GenomeIDs sicherstellen (FK-Voraussetzung)
    database.insert_database_genome_ids(
        config.database_directory,
        genome_id_set,
    )

    # 2) Metagenomes einfügen
    database.insert_database_metagenomes(
        config.database_directory,
        meta_dict,
    )

    logger.info(f"Executing read assignments for {len(task_list)} tasks")
    # Perform the graft for all tasks in the list
    graft_mp.graft_mp(task_list, config.glob_chunks, config)

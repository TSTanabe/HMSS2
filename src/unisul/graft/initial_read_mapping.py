from typing import Any

from graftm.external_program_suite import ExternalProgramSuite

from unisul.core import queue
from unisul.core.logging import get_logger
from unisul.db import database
from unisul.db.database import create_database
from unisul.graft import generate_task, graft_mp, read_counter

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
    logger.info(f"Initializing queue")
    queue.queue_read_mapping_faa_inputs(config)  # declares config.faa_files
    queue.queue_read_mapping_fna_inputs(config)  # declares config.fna_files
    queue.queue_read_mapping_fastq_inputs(config)  # declares config.fastq_files

    create_database(config.database_directory)
    task_list, gpkg_length_dict = generate_task.initialize_task_list(config)

    # Insert gpkg lengths into database
    database.insert_database_gpkg_lengths(
        config.database_directory,
        gpkg_length_dict,
    )
    logger.info("Counting reads for all (meta-)genomes")

    # deduplicate by metagenome_id so each fastq is counted once
    unique_tasks_by_meta: dict[str, Any] = {}
    for t in task_list:
        # keep first occurrence; all have same forward/reverse for same metagenome_id
        unique_tasks_by_meta.setdefault(t.metagenome_id, t)

    dedup_tasks = list(unique_tasks_by_meta.values())

    meta_dict, genome_id_set = read_counter.collect_metagenome_counts_parallel(
        dedup_tasks, processes=4, chunksize=4
    )  # TODO hier gibt es noch probleme beim zählen. das hängt sich irgendwie auf. Möglicherweise das zählen aussetzen

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

    # Perform the graft for all tasks in the list
    if config.rm_ram_limit_max:
        logger.info(
            f"Executing read assignments for {len(task_list)} tasks with {config.rm_ram_limit_max} GB RAM"
        )
        graft_mp.graft_mp_tokenized_executor(task_list, config.glob_chunks, config)
    else:
        logger.info(f"Executing read assignments for {len(task_list)} tasks")
        graft_mp.graft_mp(task_list, config.glob_chunks, config)

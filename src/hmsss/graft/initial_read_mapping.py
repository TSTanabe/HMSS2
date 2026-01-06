#!/usr/bin/env python3
import sys
import os
import traceback
from multiprocessing import get_context
from graftm.external_program_suite import ExternalProgramSuite

from hmsss.graft import graft_runner, read_counter, gpkg_length, prepare_packages
from hmsss.core import queue
from hmsss.graft import generate_task

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


def _run_graft_task(task):
    """
    Worker: Namespace bauen -> Run(args).main()
    Catch errors without stopping the multiprocessing pool
    """
    # Execute the graftM read mapping
    try:
        args = generate_task.build_graft_args(task)
        logger.debug(args)
        forward_read_number = read_counter.safe_read_count(task.forward)
        reverse_read_number = read_counter.safe_read_count(task.reverse)
        hmm_length = task.length
        logger.info(hmm_length, forward_read_number, reverse_read_number)

        graft_runner.Run(args).main()

    except SystemExit as e:
        # graftM verwendet exit() an mehreren Stellen
        return {"ok": False, "task": task, "error": f"SystemExit({e.code})"}

    except Exception:
        return {"ok": False, "task": task, "error": traceback.format_exc()}

    # parse the graftM read mapping results
    # TODO hier noch die Werte berechnen lassen die wir später haben wollen, wie TPM, RPKM, FPKM usw.
    # args.output_directory # Directory with the detected reads
    # Für jedes Taxonomie level einzeln berechnen und die Gesamtheit.
    path = os.path.join(args.output_directory, "combined_count_table.txt")
    if os.path.isfile(path):
        print("")

    return {"ok": True, "task": task, "error": None}


def initial_read_mapping(config):
    # _check_dependencies()  # Check if graftM dependencies are present
    queue.queue_read_mapping_faa_inputs(config)  # declares config.faa_files
    queue.queue_read_mapping_fna_inputs(config)  # declares config.fna_files
    queue.queue_read_mapping_fastq_inputs(config)  # declares config.fastq_files
    combined_inputs = generate_task.merge_read_mapping_inputs(
        fna_files=config.fna_files,
        faa_files=config.faa_files,
        fastq_files=config.fastq_files,
    )
    forward_dict, reverse_dict = generate_task.automatic_forward_reverse_file_detection(
        files=combined_inputs, forward_extension=None, reverse_extension=None
    )
    gpkg_packages = prepare_packages.prepare_gpkg_packages(config)
    prepare_packages.initialize_gpkg_packages(gpkg_packages, threads=4)
    gpkg_length_dict = gpkg_length.collect_gpkg_reference_median_lengths(gpkg_packages)
    # create task list can also define the cpu threads and the minimal e value
    task_list = generate_task.create_task_list(
        gpkg_packages=gpkg_packages,
        forward_dict=forward_dict,
        reverse_dict=reverse_dict,
        output_directory=config.fasta_output_directory,
        length_dict=gpkg_length_dict,
    )

    # TODO Generate the database for the mapped data

    # for each task make a process that generates the argument space for the task and runs the
    ctx = get_context("spawn")
    ok, fail = [], []

    with ctx.Pool(processes=config.cores - 1, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(_run_graft_task, task_list, chunksize=1):
            if res["ok"]:
                ok.append(res)
            else:
                fail.append(res)
                # kurz loggen und weiter
                t = res["task"]
                print(
                    "FAIL:",
                    getattr(t, "gpkg", "?"),
                    getattr(t, "forward", "?"),
                    str(res["error"]).splitlines()[-1],
                )

    # in return summarize the results for db input

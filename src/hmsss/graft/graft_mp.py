#!/usr/bin/env python3
import os
import traceback
from multiprocessing import get_context

from hmsss.cli.config import Config
from hmsss.db import database
from hmsss.graft import generate_task, graft_runner, read_models
from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


def _run_graft_task(task):
    """
    Worker: Namespace bauen -> Run(args).main()
    Catch errors without stopping the multiprocessing pool
    """
    # Execute the graftM read mapping
    try:
        args = generate_task.build_graft_args(task)
        # logger.debug(args)
        # forward_read_number = read_counter.safe_read_count(task.forward)
        # reverse_read_number = read_counter.safe_read_count(task.reverse)
        hmm_length = task.length
        min_coverage = task.min_coverage
        # logger.info(hmm_length, forward_read_number, reverse_read_number)

        # Make the graftM read assignment
        filepaths = graft_runner.Run(args).main()

        # Parse result files to read dictionary
        files = filepaths[0]
        taxonomy_csv = files.get("taxonomy")
        sequence_fasta = files.get("sequences")
        alignment_fasta = files.get("alignment")

        if (
                not taxonomy_csv
                or not sequence_fasta
                or not alignment_fasta
                or not os.path.isfile(alignment_fasta)
                or not os.path.isfile(taxonomy_csv)
                or not os.path.isfile(sequence_fasta)
        ):
            raise FileNotFoundError(
                f"Missing or empty input file(s): "
                f"taxonomy={taxonomy_csv}, "
                f"sequence={sequence_fasta}, "
                f"alignment={alignment_fasta}"
            )

        read_dict = read_models.build_reads_from_outputs(
            gpkg_name=task.gpkg_name,
            taxonomy_csv=taxonomy_csv,
            alignment_fasta=alignment_fasta,
            sequence_fasta=sequence_fasta,
            hmm_length=hmm_length,
            min_coverage=min_coverage,
        )
        _dump_read_batch(read_dict)
        for read in read_dict.values():
            read.metagenomeID = task.metagenomeID
            read.genomeID = task.genome_id

        return {
            "ok": True,
            "task": task,
            "error": "",
            "base": files.get("base"),
            "read_dict": read_dict,
        }

    except SystemExit as e:
        # graftM verwendet exit() an mehreren Stellen
        return {
            "ok": False,
            "task": task,
            "error": f"SystemExit({e.code})",
            "error_code": False,
            "base": "",
            "read_dict": {},
        }

    except Exception:
        return {
            "ok": False,
            "task": task,
            "error": traceback.format_exc(),
            "error_code": True,
            "base": "",
            "read_dict": {},
        }


def _flush_read_batch_to_db(database_path: str, read_batch: dict) -> None:
    # 1) Lineages (INSERT OR IGNORE)
    database.insert_database_lineages(database_path, read_batch)

    # 2) Stub proteins (INSERT OR IGNORE) – nur nötig, wenn Placement.proteinID nicht NULL sein soll
    database.insert_database_stub_proteins_from_reads(database_path, read_batch)

    # 3) Placements (INSERT OR IGNORE)
    database.insert_database_placements(database_path, read_batch)


def _dump_read_batch(read_batch: dict, *, max_reads: int | None = None) -> None:
    print(f"[DEBUG] read_batch contains {len(read_batch)} reads")

    for i, (rid, r) in enumerate(read_batch.items(), start=1):
        print(f"\n--- READ {i} ---")
        print(f"dict key: {rid}")

        # Alle Attribute des Read-Objekts anzeigen
        if hasattr(r, "__dict__"):
            for k, v in r.__dict__.items():
                print(f"  {k}: {v}")
        else:
            print("  [no __dict__] repr:", repr(r))

        if max_reads is not None and i >= max_reads:
            print(f"\n[DEBUG] stopped after {max_reads} reads")
            break


def graft_mp(task_list: list, batch_size: int, config: Config) -> None:
    read_batch: dict[str, Read] = {}
    batch_counter: int = 0

    # Fortschritt
    read_mappings_done: int = 0
    n_tasks = len(task_list)
    log_step = max(1, n_tasks // 100)  # ~5% Schritte
    # for each task make a process that generates the argument space for the task and runs the
    ctx = get_context("spawn")

    with ctx.Pool(processes=config.cores - 1, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(_run_graft_task, task_list, chunksize=1):
            read_mappings_done += 1
            if not res["ok"]:
                continue

            # Progress logger
            if (read_mappings_done % log_step == 0) or (read_mappings_done == n_tasks):
                pct = (read_mappings_done * 100) // max(1, n_tasks)
                logger.info(
                    f"[Read-mapping progress] {read_mappings_done}/{n_tasks} ({pct}%) tasks processed"
                )

            # Buffer read dicts
            res_dict = res["read_dict"]
            read_batch.update(res_dict)
            batch_counter += 1

            print("Batch batch_counter:", batch_counter)
            print("Length", len(read_batch))
            _dump_read_batch(read_batch)

            # Insert into database
            if batch_counter >= batch_size:
                _flush_read_batch_to_db(
                    database_path=config.database_directory, read_batch=read_batch
                )
                read_batch.clear()
                batch_counter = 0

        # Insert the remaining reads
        if read_batch:
            _flush_read_batch_to_db(
                database_path=config.database_directory, read_batch=read_batch
            )
            read_batch.clear()
    # in return summarize the results for db input

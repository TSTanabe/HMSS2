#!/usr/bin/env python3
import os
import resource
import time
import traceback
from multiprocessing import get_context

from hmsss.cli.config import Config
from hmsss.db import database
from hmsss.graft import generate_task, graft_runner, read_models
from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


def _set_worker_mem_limit_gb(limit_gb: float) -> None:
    """
    Set a hard per-process address-space limit (RLIMIT_AS).
    Works on Linux/Unix. Subprocesses typically inherit this limit.
    """
    if limit_gb <= 0:
        return

    bytes_limit = int(limit_gb * 1024 ** 3)

    # Hard+soft limit
    resource.setrlimit(resource.RLIMIT_AS, (bytes_limit, bytes_limit))


def _run_graft_task(task):
    """
    Worker: Namespace bauen -> Run(args).main()
    Catch errors without stopping the multiprocessing pool
    """
    # Execute the graftM read mapping
    try:
        # Set RAM limit in GB
        mem_cap_gb = getattr(task, "mem_cap_gb", None)
        if mem_cap_gb is not None:
            _set_worker_mem_limit_gb(float(mem_cap_gb))

        args = generate_task.build_graft_args(task)
        # logger.debug(args)
        # Commented because is done separately in the main process
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

        for read in read_dict.values():
            read.metagenomeID = task.metagenome_id
            read.genomeID = task.genome_id
        # _dump_read_batch(read_dict)
        return {
            "ok": True,
            "task": task,
            "error": "",
            "base": files.get("base"),
            "read_dict": read_dict,
        }

    except SystemExit as e:
        # graftM verwendet exit() an mehreren Stellen
        print("SYSTEM EXIT ERROR")
        return {
            "ok": False,
            "task": task,
            "error": f"SystemExit({e.code})",
            "error_code": False,
            "base": "",
            "read_dict": {},
        }

    except MemoryError as e:
        print("MEMORY EXIT ERROR")
        print("EXCEPTION:", e)
        print(traceback.format_exc())
        return {
            "ok": False,
            "task": task,
            "error": traceback.format_exc(),
            "error_code": True,
            "base": "",
            "read_dict": {},
        }

    except Exception as e:
        print("EXCEPTION:", e)
        print(traceback.format_exc())
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

            # print("Batch batch_counter:", batch_counter)
            # print("Length", len(read_batch))
            # _dump_read_batch(read_batch)
            # print("----- original res dict -----")
            # _dump_read_batch(res_dict)

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


def _task_mem_est_gb(task, default_gb: float = 4.0) -> float:
    v = getattr(task, "mem_est_gb", None)
    try:
        return float(v) if v is not None else float(default_gb)
    except Exception:
        return float(default_gb)


import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import get_context

from concurrent.futures import wait, FIRST_COMPLETED


def _collect_done_futures(
        *,
        future_to_tokens: dict,
        # future_to_task: dict,
        timeout: float | None,
        available_tokens_gb: float,
        handle_result_fn,
) -> float:
    """
    Collect finished futures, free reserved tokens, and process results.

    Args:
        timeout:
            - None  => block until at least one future finishes (if any exist)
            - 0.0   => non-blocking poll
            - >0.0  => wait up to that many seconds
    Returns:
        Updated available_tokens_gb
    """
    if not future_to_tokens:
        return available_tokens_gb

    done, _ = wait(
        future_to_tokens.keys(),
        timeout=timeout,
        return_when=FIRST_COMPLETED,
    )

    for fut in done:
        reserved = future_to_tokens.pop(fut)
        # future_to_task.pop(fut, None)
        available_tokens_gb += reserved

        try:
            res = fut.result()
        except Exception as e:
            print(e)
            res = {"ok": False}

        handle_result_fn(res)

    return available_tokens_gb


def graft_mp_tokenized_executor(task_list: list, batch_size: int, config: Config) -> None:
    """
    Token-aware multiprocessing using ProcessPoolExecutor + wait(FIRST_COMPLETED).

    Properties:
      - Skip tasks whose mem_est_gb > ram_budget_gb (warn)
      - Start tasks only when enough RAM tokens are available
      - Free tokens as soon as ANY task finishes (FIRST_COMPLETED)
      - Same DB flush logic as graft_mp / graft_mp_tokenized
    """
    read_batch: dict[str, Read] = {}
    batch_counter: int = 0

    read_mappings_done: int = 0
    n_tasks = len(task_list)
    log_step = max(1, n_tasks // 100)

    # ---- RAM token budget (GB) ----
    total_tokens_gb = float(getattr(config, "ram_budget_gb", 8.0))
    available_tokens_gb = total_tokens_gb

    # Optional: largest-first gegen Fragmentierung
    task_list = sorted(task_list, key=lambda t: _task_mem_est_gb(t), reverse=True)

    # spawn context
    ctx = get_context("spawn")
    max_workers = max(1, int(config.cores) - 1)

    # in-flight: future -> reserved_tokens_gb
    future_to_tokens: dict = {}

    def _handle_result(res: dict) -> None:
        nonlocal read_mappings_done, batch_counter, read_batch

        read_mappings_done += 1

        # Progress logger
        if (read_mappings_done % log_step == 0) or (read_mappings_done == n_tasks):
            pct = (read_mappings_done * 100) // max(1, n_tasks)
            logger.info(f"[Read-mapping progress] {read_mappings_done}/{n_tasks} ({pct}%) tasks processed")

        if not res.get("ok"):
            return

        res_dict = res["read_dict"]
        read_batch.update(res_dict)
        batch_counter += 1

        if batch_counter >= batch_size:
            _flush_read_batch_to_db(database_path=config.database_directory, read_batch=read_batch)
            read_batch.clear()
            batch_counter = 0

    with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as ex:
        idx = 0

        while idx < len(task_list) or future_to_tokens:
            # 1) So viele Tasks starten wie möglich (Tokens & worker budget)
            started_any = False
            while idx < len(task_list):
                task = task_list[idx]
                need = _task_mem_est_gb(task)

                # Policy: zu groß für den Node => verwerfen
                if need > total_tokens_gb:
                    logger.warning(
                        f"[RAM tokens] Skipping task {getattr(task, 'gpkg_name', '?')} "
                        f"(mem_est_gb={need:.1f} > ram_budget_gb={total_tokens_gb:.1f})"
                    )
                    idx += 1
                    # zählt als "processed", damit Progress/Ende konsistent bleibt
                    read_mappings_done += 1
                    continue

                # nicht genug Tokens frei => jetzt nicht starten
                if need > available_tokens_gb:
                    break

                # Workerzahl nicht überschreiten: executor blockiert nicht,
                # aber wir begrenzen in-flight, um Token-Accounting stabil zu halten
                if len(future_to_tokens) >= max_workers:
                    break

                available_tokens_gb -= need
                fut = ex.submit(_run_graft_task, task)
                future_to_tokens[fut] = need
                idx += 1
                started_any = True

            # 2) Wenn nichts startbar ist: auf mindestens EIN fertiges Future warten
            if future_to_tokens and (not started_any):
                available_tokens_gb = _collect_done_futures(
                    future_to_tokens=future_to_tokens,
                    timeout=None,  # <- blockierend bis FIRST_COMPLETED
                    available_tokens_gb=available_tokens_gb,
                    handle_result_fn=_handle_result,
                )
                continue

            # 3) Wenn wir gestartet haben, können wir ebenfalls fertige Futures einsammeln (ohne zu blockieren)
            available_tokens_gb = _collect_done_futures(
                future_to_tokens=future_to_tokens,
                timeout=0.0,  # <- non-blocking poll
                available_tokens_gb=available_tokens_gb,
                handle_result_fn=_handle_result,
            )

            # kleine Pause gegen Busy-wait
            if future_to_tokens:
                time.sleep(0.02)

    # remaining reads flush
    if read_batch:
        _flush_read_batch_to_db(database_path=config.database_directory, read_batch=read_batch)
        read_batch.clear()

#!/usr/bin/env python3
import os

# --- hard cap for nested threading libraries (critical under spawn + ProcessPoolExecutor) ---
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import resource
import traceback
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import get_context

from typing import Callable, Any
from hmsss.cli.config import Config
from hmsss.db import database
from hmsss.graft import generate_task, graft_runner, read_models
from hmsss.core.logging import get_logger
from hmsss.graft.read_models import Read

logger = get_logger(__name__)


def _worker_init_thread_limits():
    import os
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"


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
        # Set RAM limit in GB, RSS limiter thorughs exceptions when additional libraries are loaded
        # mem_cap_gb = getattr(task, "mem_cap_gb", None)
        # if mem_cap_gb is not None:
        #    _set_worker_mem_limit_gb(float(mem_cap_gb))

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
        non_decoy_sequences = files.get("non_decoy_sequences")

        if (
                not taxonomy_csv
                or not sequence_fasta
                or not alignment_fasta
                or not os.path.isfile(alignment_fasta)
                or not os.path.isfile(taxonomy_csv)
                # or not os.path.isfile(sequence_fasta) does not work if fw & rv are provided
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
        _dump_read_batch(read_dict)
        return {
            "ok": True,
            "task": task,
            "error": "",
            "base": files.get("base"),
            "read_dict": read_dict,
        }

    except SystemExit as e:
        # graftM verwendet exit() an mehreren Stellen
        print("SYSTEM EXIT ERROR:", e)
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


def _start_bestfit_tasks(
        *,
        ex: object,
        pending: list,
        future_to_tokens: dict,
        available_tokens_gb: float,
        total_tokens_gb: float,
        max_workers: int,
        k_scan: int = 100,
        handle_result_fn: Callable[[dict], None],
) -> tuple[float, bool]:
    """
    Start tasks using Best-Fit from Top-K pending tasks.

    Returns:
      (updated_available_tokens_gb, did_progress)

    did_progress is True if we either:
      - started at least one task, or
      - skipped (discarded) at least one oversize task.
    """
    did_progress = False

    while pending and (len(future_to_tokens) < max_workers):
        best_i = None
        best_need = 0.0

        scan_n = min(k_scan, len(pending))

        # Scan Top-K for best-fit
        for i in range(scan_n):
            task = pending[i]
            need = _task_mem_est_gb(task)

            # Policy: oversize -> discard with warning
            if need > total_tokens_gb:
                logger.warning(
                    f"[RAM tokens] Skipping task {getattr(task, 'gpkg_name', '?')} "
                    f"(ESTIMATED NEEDED RAM={need:.1f} GB > RAM LIMIT={total_tokens_gb:.1f}) GB"
                )
                pending.pop(i)
                handle_result_fn({"ok": False})  # counts as processed consistently
                did_progress = True
                best_i = None
                best_need = 0.0
                break  # restart scanning (indices shifted)

            # Best fit: largest that still fits
            if need <= available_tokens_gb and need > best_need:
                best_need = need
                best_i = i

        # If we removed an oversize task, restart the outer while and try again
        if best_i is None and best_need == 0.0 and did_progress:
            continue

        # No task fits into remaining tokens right now -> stop starting
        if best_i is None:
            break

        # Start the selected task
        task = pending.pop(best_i)
        logger.info(
            f"Starting task {task.gpkg_name} Need: {need:.1f} GB; Available: {available_tokens_gb:.1f} GB; Total tokens: {total_tokens_gb:.1f} GB"
        )

        available_tokens_gb -= best_need

        fut = ex.submit(_run_graft_task, task)
        fut._hmss_gpkg = getattr(task, "gpkg_name", "?")
        fut._hmss_mem_gb = best_need

        future_to_tokens[fut] = best_need
        did_progress = True

    return available_tokens_gb, did_progress


import time
from typing import Callable


def _start_bestfit_tasks_debug(
        *,
        ex: object,
        pending: list,
        future_to_tokens: dict,
        available_tokens_gb: float,
        total_tokens_gb: float,
        max_workers: int,
        k_scan: int = 100,
        handle_result_fn: Callable[[dict], None],
        debug_scan: int = 12,  # how many candidates to print per tick
        debug_level: str = "INFO",  # "INFO" or "DEBUG"
) -> tuple[float, bool]:
    """
    Debug-instrumented variant of _start_bestfit_tasks().
    Prints scheduler state and decision reasons (worker slots vs token fit).
    """
    did_progress = False

    def _log(msg: str, *args):
        # Use logger if you want; prints are simplest for HPC logs
        # Swap to logger.info/debug if preferred.
        print(msg % args if args else msg, flush=True)

    # Print entry status once per call
    _log(
        "[SCHED] enter: pending=%d running=%d max_workers=%d avail=%.1f total=%.1f k_scan=%d",
        len(pending),
        len(future_to_tokens),
        max_workers,
        available_tokens_gb,
        total_tokens_gb,
        k_scan,
    )

    # If no worker slots, we can immediately explain why nothing starts
    if len(future_to_tokens) >= max_workers:
        _log(
            "[SCHED] no worker slot: running=%d >= max_workers=%d (tokens avail=%.1f)",
            len(future_to_tokens),
            max_workers,
            available_tokens_gb,
        )
        return available_tokens_gb, False

    while pending and (len(future_to_tokens) < max_workers):
        t0 = time.time()

        best_i = None
        best_need = 0.0

        scan_n = min(k_scan, len(pending))

        # Quick peek at the top-N tasks and their estimated RAM
        peek_n = min(debug_scan, scan_n)
        peek = []
        for j in range(peek_n):
            tj = pending[j]
            nj = _task_mem_est_gb(tj)
            peek.append((getattr(tj, "gpkg_name", "?"), nj))
        _log(
            "[SCHED] scan peek top-%d/%d (avail=%.1f): %s",
            peek_n,
            scan_n,
            available_tokens_gb,
            ", ".join([f"{n}:{gb:.1f}" for n, gb in peek]) if peek else "(none)",
        )

        removed_oversize = False

        # Scan Top-K for best-fit
        for i in range(scan_n):
            task = pending[i]
            need_i = _task_mem_est_gb(task)
            name = getattr(task, "gpkg_name", "?")

            # Policy: oversize -> discard with warning
            if need_i > total_tokens_gb:
                _log(
                    "[SCHED] OVERSIZE -> drop %s (need=%.1f > total=%.1f). pending before=%d",
                    name,
                    need_i,
                    total_tokens_gb,
                    len(pending),
                )
                pending.pop(i)
                handle_result_fn({"ok": False})
                did_progress = True
                removed_oversize = True
                break  # indices shifted; restart scan

            # Best fit: largest that still fits
            if need_i <= available_tokens_gb and need_i > best_need:
                best_need = need_i
                best_i = i

        # If we removed an oversize task, restart the loop
        if removed_oversize:
            _log(
                "[SCHED] restart scan after oversize removal (pending now=%d)",
                len(pending),
            )
            continue

        # No task fits into remaining tokens right now -> stop starting
        if best_i is None:
            # Distinguish between "tokens too low" and "no candidates" (rare)
            # Estimate smallest need in scanned window to explain quickly.
            min_need = None
            for i in range(scan_n):
                ni = _task_mem_est_gb(pending[i])
                min_need = ni if (min_need is None or ni < min_need) else min_need

            _log(
                "[SCHED] no fit: avail=%.1f, scanned=%d, min_need_in_scan=%s; running=%d/%d",
                available_tokens_gb,
                scan_n,
                f"{min_need:.1f}" if min_need is not None else "n/a",
                len(future_to_tokens),
                max_workers,
            )
            break

        # Start the selected task
        task = pending.pop(best_i)
        name = getattr(task, "gpkg_name", "?")

        _log(
            "[SCHED] START %s need=%.1f avail_before=%.1f total=%.1f running=%d/%d pending_left=%d",
            name,
            best_need,
            available_tokens_gb,
            total_tokens_gb,
            len(future_to_tokens),
            max_workers,
            len(pending),
        )

        available_tokens_gb -= best_need

        fut = ex.submit(_run_graft_task, task)
        fut._hmss_gpkg = name
        fut._hmss_mem_gb = best_need
        future_to_tokens[fut] = best_need

        did_progress = True

        _log(
            "[SCHED] submitted %s; avail_after=%.1f running_now=%d/%d (tick=%.3fs)",
            name,
            available_tokens_gb,
            len(future_to_tokens),
            max_workers,
            time.time() - t0,
        )

    _log(
        "[SCHED] exit: did_progress=%s pending=%d running=%d avail=%.1f",
        did_progress,
        len(pending),
        len(future_to_tokens),
        available_tokens_gb,
    )
    return available_tokens_gb, did_progress


def graft_mp_tokenized_executor(
        task_list: list, batch_size: int, config: "Config"
) -> None:
    read_batch: dict[str, "Read"] = {}
    batch_counter: int = 0

    processed: int = 0
    n_tasks = len(task_list)
    log_step = max(1, n_tasks // 100)

    total_tokens_gb = float(getattr(config, "rm_ram_limit_max", 16.0))
    available_tokens_gb = total_tokens_gb

    pending = sorted(task_list, key=lambda t: _task_mem_est_gb(t), reverse=True)

    ctx = get_context("spawn")
    max_workers = max(1, (int(config.cores) - 1) // 5)

    future_to_tokens: dict = {}

    def _handle_result(res: dict) -> None:
        nonlocal processed, batch_counter, read_batch

        processed += 1
        if (processed % log_step == 0) or (processed == n_tasks):
            pct = (processed * 100) // max(1, n_tasks)
            logger.info(
                f"[Read-mapping progress] {processed}/{n_tasks} ({pct}%) tasks processed"
            )

        if not res.get("ok"):
            return

        read_batch.update(res["read_dict"])
        batch_counter += 1

        if batch_counter >= batch_size:
            _flush_read_batch_to_db(
                database_path=config.database_directory, read_batch=read_batch
            )
            read_batch.clear()
            batch_counter = 0

    with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx, initializer=_worker_init_thread_limits()) as ex:
        while pending or future_to_tokens:
            # 1) Refill: start as many as possible
            available_tokens_gb, _ = _start_bestfit_tasks_debug(
                ex=ex,
                pending=pending,
                future_to_tokens=future_to_tokens,
                available_tokens_gb=available_tokens_gb,
                total_tokens_gb=total_tokens_gb,
                max_workers=max_workers,
                k_scan=100,
                handle_result_fn=_handle_result,
            )

            # 2) If nothing is running, we're done (or only oversize were skipped)
            if not future_to_tokens:
                break

            # 3) Always block until at least one finishes
            done, _ = wait(future_to_tokens.keys(), return_when=FIRST_COMPLETED)

            for fut in done:
                reserved = future_to_tokens.pop(fut)
                available_tokens_gb += reserved
                logger.debug(
                    "Reserved tokens returned: %.1f GB | gpkg=%s",
                    reserved,
                    getattr(fut, "_hmss_gpkg", "?"),
                )
                try:
                    res = fut.result()
                except Exception as e:
                    print(e)
                    res = {"ok": False}
                _handle_result(res)

    if read_batch:
        _flush_read_batch_to_db(
            database_path=config.database_directory, read_batch=read_batch
        )
        read_batch.clear()

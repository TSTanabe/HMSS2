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

import os
import time
import csv
import threading
from dataclasses import dataclass, asdict
from pathlib import Path

import psutil


@dataclass
class Sample:
    t_s: float
    cpu_percent: float
    rss_bytes: int
    vms_bytes: int
    read_bytes: int
    write_bytes: int
    read_count: int
    write_count: int
    num_fds: int | None = None


class ProcSampler:
    """
    Periodisches Sampling von CPU/RAM/IO des *aktuellen* Prozesses.
    Schreibt nach Ende optional CSV.
    """

    def __init__(self, *, interval_s: float = 0.2, out_csv: str | None = None):
        self.interval_s = float(interval_s)
        self.out_csv = out_csv
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.samples: list[Sample] = []
        self._proc = psutil.Process(os.getpid())
        self._t0 = None

    def __enter__(self):
        self._t0 = time.perf_counter()

        # cpu_percent braucht eine "Baseline"-Messung
        self._proc.cpu_percent(interval=None)

        self._thread = threading.Thread(target=self._run, name="proc-sampler", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2.0)

        if self.out_csv:
            self._write_csv(self.out_csv)

        # Exceptions NICHT schlucken – das machst du außen in deinem try/except
        return False

    def _run(self):
        while not self._stop.is_set():
            self._take_sample()
            time.sleep(self.interval_s)

    def _take_sample(self):
        t_s = time.perf_counter() - self._t0

        mem = self._proc.memory_info()
        io = self._proc.io_counters()  # read/write bytes + counts (plattformabhängig zuverlässig unter Linux)
        try:
            num_fds = self._proc.num_fds()
        except Exception:
            num_fds = None

        s = Sample(
            t_s=t_s,
            cpu_percent=self._proc.cpu_percent(interval=None),
            rss_bytes=mem.rss,
            vms_bytes=getattr(mem, "vms", 0),
            read_bytes=getattr(io, "read_bytes", 0),
            write_bytes=getattr(io, "write_bytes", 0),
            read_count=getattr(io, "read_count", 0),
            write_count=getattr(io, "write_count", 0),
            num_fds=num_fds,
        )
        self.samples.append(s)

    def _write_csv(self, out_csv: str):
        outp = Path(out_csv)
        outp.parent.mkdir(parents=True, exist_ok=True)
        with outp.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(self.samples[0]).keys()) if self.samples else [])
            if self.samples:
                w.writeheader()
                for s in self.samples:
                    w.writerow(asdict(s))


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

        prof_dir = Path(task.outdir) / "profile"
        prof_csv = prof_dir / f"graft_run_pid{os.getpid()}.csv"
        t0 = time.perf_counter()
        with ProcSampler(interval_s=0.2, out_csv=str(prof_csv)) as sampler:
            graft_runner.Run(args).main()
        wall = time.perf_counter() - t0
        print("Profile")
        print(
            f"graftM done | wall={wall:.2f}s | samples={len(sampler.samples)} | profile={prof_csv}",
            flush=True
        )

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

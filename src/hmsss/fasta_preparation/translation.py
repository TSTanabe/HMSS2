#!/usr/bin/python

import os
import re
import multiprocessing
import gzip
import shutil
import tempfile
from pathlib import Path

from typing import Tuple, Set

from hmsss.core.logging import get_logger
from hmsss.utils import myUtil

logger = get_logger(__name__)


def parallel_translation(fna_files: dict[str, str], cores: int) -> None:
    """
    Uses prodigal to translate all nucleotide fasta files in a directory to protein fasta (.faa).

    Args:
        fna_files (dict[str,str]): Dictionary of fasta files to translate.:
        cores (int): Number of CPU cores to use (multiprocessing).

    Output:
        Translated .faa files are written next to source files.

    Example:
        parallel_translation('./genomes', 8)

    Warning:
        If the filepaths include parentheses prodigal is not working
    """

    logger.info(
        f"Found {len(fna_files)} assemblies in nucleotide or ambiguous format for prodigal"
    )

    manager = multiprocessing.Manager()
    counter = manager.Value("i", 0)
    lock = manager.Lock()
    length = len(fna_files)
    prodigal = myUtil.find_executable("prodigal")

    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [
            (fna_fasta, length, counter, lock, prodigal)
            for fna_fasta in fna_files.values()
        ]
        pool.map(translate_fasta, args_list)
    logger.info(f"Processing assembly {counter.value} of {length}")
    return


def translate_fasta(
    args: Tuple[str, int, multiprocessing.Value, multiprocessing.Lock, str],
) -> None:
    """
    Runs prodigal for a single fasta file.
    If fasta is .gz -> decompress to a temporary file first, then translate, then delete temp.
    """
    fasta, length, counter, lock, prodigal = args

    fasta_p = Path(fasta)

    # output should be next to the source file, but WITHOUT the .gz suffix in the basename
    # e.g. genome.fna.gz -> genome.faa
    stem = fasta_p.name
    if stem.endswith(".gz"):
        stem = stem[:-3]  # remove .gz
    output_base = fasta_p.with_name(Path(stem).stem)  # remove .fna or other last suffix
    faa = str(output_base) + ".faa"

    tmp_path = None
    input_for_prodigal = str(fasta_p)

    try:
        if fasta_p.name.endswith(".gz"):
            # decompress to temp file (same filesystem preferred; /tmp is fine too)
            # suffix helps prodigal / debugging
            with tempfile.NamedTemporaryFile(
                mode="wb", suffix=".fna", delete=False
            ) as tmp:
                tmp_path = tmp.name
                with gzip.open(str(fasta_p), "rb") as f_in:
                    shutil.copyfileobj(f_in, tmp)

            input_for_prodigal = tmp_path

        cmd = f'{prodigal} -a "{faa}" -i "{input_for_prodigal}" >/dev/null 2>&1'
        os.system(cmd)

    except Exception as e:
        logger.warning(f"Could not translate {fasta} - {e}")
        return

    finally:
        if tmp_path is not None:
            try:
                os.unlink(tmp_path)
            except Exception:
                pass

    with lock:
        counter.value += 1
        print(
            f"Prodigal processing assembly {counter.value} of {length}",
            end="\r",
            flush=True,
        )
    return


############################################################################
############### Parallel Transcription #####################################
############################################################################


def parallel_transcription(faa_files: dict[str, str], cores: int) -> None:
    """
    8.10.22
        Args:
            cores: Number of cores
            faa_files: Dictionary of fasta files to translate.

        Transcribe for all faa files gff3 files
        Secure a packed and unpacked version is present
        Unlink unpacked versions afterward
        Warning: if the directory path includes parentheses function prodigal is not working
        For all .faa files in a directory, transcribe GFF3 files using prodigal header format.

    """

    manager = multiprocessing.Manager()
    counter = manager.Value("i", 0)
    lock = manager.Lock()
    length = len(faa_files)

    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [
            (faa_fasta, length, counter, lock) for faa_fasta in faa_files.values()
        ]
        pool.map(transcripe_fasta, args_list)

    logger.info(f"Generated corresponding gff files")
    return


def _output_prefix_from_original(fasta_p: Path) -> str:
    """
    genome.faa      -> genome
    genome.faa.gz   -> genome
    """
    name = fasta_p.name
    if name.endswith(".gz"):
        name = name[:-3]
    return str(fasta_p.with_name(Path(name).stem))


def transcripe_fasta(
    args: Tuple[str, int, multiprocessing.Value, multiprocessing.Lock],
) -> None:
    fasta, length, counter, lock = args
    fasta_p = Path(fasta)

    tmp_path: str | None = None
    input_for_processing = str(fasta_p)

    # Output-GFF soll neben dem ORIGINAL liegen:
    # genome.faa      -> genome.gff
    # genome.faa.gz   -> genome.gff
    out_base = fasta_p
    if out_base.name.endswith(".gz"):
        out_base = out_base.with_name(out_base.name[:-3])  # strip .gz
    if out_base.name.endswith(".faa"):
        out_gff = str(out_base.with_suffix(".gff"))
    else:
        out_gff = str(out_base) + ".gff"

    try:
        if fasta_p.name.endswith(".gz"):
            with tempfile.NamedTemporaryFile(
                mode="wb", suffix=".faa", delete=False
            ) as tmp:
                tmp_path = tmp.name
                with gzip.open(str(fasta_p), "rb") as f_in:
                    shutil.copyfileobj(f_in, tmp)
            input_for_processing = tmp_path

        if check_prodigal_format(input_for_processing):
            prodigal_faa_to_gff(input_for_processing, output_gff=out_gff)

            with lock:
                counter.value += 1
                print(
                    f"Processing file {counter.value} of {length}",
                    end="\r",
                    flush=True,
                )

    finally:
        if tmp_path is not None:
            try:
                os.unlink(tmp_path)
            except Exception:
                pass


def check_prodigal_format(gff_file: str) -> int:
    """
    Checks if a fasta file contains prodigal header (5 fields).

    Args:
        gff_file (str): Path to fasta file.

    Returns:
        int: 1 if prodigal format, 0 otherwise.
    """

    with open(gff_file, "r") as reader:
        for line in reader.readlines():
            if line[0] == ">":
                line = line[1:]
                ar = line.split("#")
                if len(ar) == 5:
                    return 1
                else:
                    return 0
    return 0


def prodigal_faa_to_gff(filepath: str, *, output_gff: str | None = None) -> str:
    """
    Translates prodigal .faa to GFF3 format. Returns the GFF3 filename.
    """
    if output_gff is None:
        dir_path = os.path.dirname(filepath)
        filename_with_ext = os.path.basename(filepath)

        if filename_with_ext.endswith(".gz"):
            filename_with_ext = os.path.splitext(filename_with_ext)[0]

        if filename_with_ext.endswith(".faa"):
            filename_without_ext = os.path.splitext(filename_with_ext)[0]
        else:
            filename_without_ext = filename_with_ext

        gff = os.path.join(dir_path, filename_without_ext + ".gff")
    else:
        gff = output_gff

    with open(gff, "w") as writer:
        with open(filepath, "r") as reader:
            genome_id = myUtil.get_genome_id(filepath)
            for line in reader.readlines():
                if line and line[0] == ">":
                    try:
                        line = line[1:]
                        ar = line.split("#")
                        contig = re.split(r"_\d+\W+$", ar[0])
                        strand = "+" if ar[3] == " 1 " else "-"
                        writer.write(
                            contig[0]
                            + "\tprodigal\tcds\t"
                            + ar[1]
                            + "\t"
                            + ar[2]
                            + "\t0.0\t"
                            + strand
                            + "\t0\tID=cds-"
                            + ar[0]
                            + ";Genome="
                            + genome_id
                            + "\n"
                        )
                    except Exception as e:
                        logger.error(f"Missformated header\n {line} - {e}")

    return gff


def compare_file_lists(directory: str, ext1: str, ext2: str) -> Set[str]:
    """
    Liefert alle Dateien mit Endung `ext1`, für die KEIN Gegenstück mit Endung `ext2` existiert.
    Unterstützt mehrstufige Endungen (z.B. ".fna.gz", ".faa.gz").

    Beispiel:
      dir:  ["genome1.fna.gz", "genome1.faa", "genome2.fna.gz"]
      Aufruf: compare_file_lists(dir, ".fna.gz", ".faa")
      Rückgabe: {"<pfad>/genome2.fna.gz"}  # genome1 hat Gegenstück, genome2 nicht
    """
    p = Path(directory)

    # mapping: Normalisierter_Basisname -> Path
    files1 = {}
    for f in p.glob(f"*{ext1}"):
        if f.is_file():
            base = _strip_suffix(f.name, ext1)  # z.B. "genome1" aus "genome1.fna.gz"
            files1[base] = f

    files2 = {}
    for f in p.glob(f"*{ext2}"):
        if f.is_file():
            base = _strip_suffix(f.name, ext2)  # z.B. "genome1" aus "genome1.faa"
            files2[base] = f

    # diejenigen ext1-Dateien, die keinen passenden ext2-Basisnamen haben
    missing = {str(files1[base]) for base in files1.keys() - files2.keys()}
    return missing


def _strip_suffix(filename: str, suffix: str) -> str:
    """
    Entfernt den exakten Suffix `suffix` nur dann, wenn er am Ende steht.
    Beispiel: ("a.b.c.fna.gz", ".fna.gz") -> "a.b.c"
             ("a.b.c.fna", ".fna.gz")     -> unverändert
    """
    return (
        filename[: -len(suffix)] if suffix and filename.endswith(suffix) else filename
    )

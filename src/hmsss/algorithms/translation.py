#!/usr/bin/python

import os
import re
import multiprocessing
from pathlib import Path
from typing import Tuple, Set

from src.hmsss.utils import myUtil

logger = myUtil.log


def parallel_translation(fna_files: dict[str,str], cores: int) -> None:
    """
    Uses prodigal to translate all nucleotide fasta files in a directory to protein fasta (.faa).

    Args:
        fna_files (dict[str,str]): Dictionary of fasta files to translate.:
        directory (str): Path to the directory with fasta/fna files.
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
            (fasta, length, counter, lock, prodigal) for fasta in fna_files
        ]
        pool.map(translate_fasta, args_list)
    logger.info(f"Processing assembly {counter.value} of {length}")
    return


def translate_fasta(
    args: Tuple[str, int, multiprocessing.Value, multiprocessing.Lock, str],
) -> None:
    """
    Runs prodigal for a single fasta file.

    Args:
        args: Tuple of (fasta path, total length, shared counter, shared lock, prodigal_executable)
    """
    fasta, length, counter, lock, prodigal = args

    # Run prodigal
    output = os.path.splitext(fasta)[0]
    faa = output + ".faa"

    string = f"{prodigal} -a {faa} -i {fasta} >/dev/null 2>&1"
    try:
        os.system(string)
    except Exception as e:
        logger.warning(f"Could not translate {fasta} - {e}")
        return

    with lock:
        counter.value += 1
        print(
            f"[INFO] Processing assembly {counter.value} of {length}",
            end="\r",
            flush=True,
        )
    return


############################################################################
############### Parallel Transcription #####################################
############################################################################


def parallel_transcription(faa_files: dict[str,str], cores: int) -> None:
    """
    8.10.22
        Args:
            directory   fasta file containing directory

        Transcribe for all faa files gff3 files
        Secure a packed and unpacked version is present
        Unlink unpacked versions afterward
        Warning: if the directory path includes parentheses function prodigal is not working
    For all .faa files in a directory, transcribe GFF3 files using prodigal header format.

    Args:
        directory (str): Directory with .faa and .gff files.
        cores (int): Number of CPU cores.

    Output:
        GFF files are generated in-place.
    """

    manager = multiprocessing.Manager()
    counter = manager.Value("i", 0)
    lock = manager.Lock()
    length = len(faa_files)

    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length, counter, lock) for fasta in faa_files]
        pool.map(transcripe_fasta, args_list)

    logger.info(f"Generated corresponding gff files")
    return


def transcripe_fasta(
    args: Tuple[str, int, multiprocessing.Value, multiprocessing.Lock],
) -> None:
    """
    Converts prodigal .faa to GFF3 format if prodigal header detected.

    Args:
        args: Tuple of (fasta path, total length, shared counter, shared lock)
    """
    fasta, length, counter, lock = args

    if check_prodigal_format(fasta):
        prodigal_faa_to_gff(fasta)

        with lock:
            counter.value += 1
            print(
                f"[INFO] Processing file {counter.value} of {length}",
                end="",
                flush=True,
            )

    return


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


def prodigal_faa_to_gff(filepath: str) -> str:
    """
    Translates prodigal .faa to GFF3 format. Returns the GFF3 filename.

    Args:
        filepath (str): Path to prodigal .faa file.

    Returns:
        str: Path to generated .gff file.
    """
    # Translates faa files from prodigal to normal gff3 formated files. returns the gff3 file name
    dir_path = os.path.dirname(filepath)
    filename_with_ext = os.path.basename(filepath)

    if filename_with_ext.endswith(".gz"):
        filename_with_ext = os.path.splitext(filename_with_ext)[0]

    # Remove the .faa extension if present
    if filename_with_ext.endswith(".faa"):
        filename_without_ext = os.path.splitext(filename_with_ext)[0]
    else:
        filename_without_ext = filename_with_ext

    gff = dir_path + "/" + filename_without_ext + ".gff"

    writer = open(gff, "w")

    with open(filepath, "r") as reader:
        genome_id = myUtil.get_genome_id(filepath)
        for line in reader.readlines():
            if line[0] == ">":
                try:
                    line = line[1:]
                    ar = line.split("#")
                    # print(ar)
                    contig = re.split(r"_{1}\d+\W+$", ar[0])
                    # print(contig)
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
    writer.close()
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

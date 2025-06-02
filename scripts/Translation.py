#!/usr/bin/python

import os
import re
import multiprocessing
from typing import Any, List, Tuple
from . import myUtil


def parallel_translation(directory,cores):
    """
    18.5.24
        Args:  
            directory   fasta file containing directory
            
        Uses prodigal to translate all nucleotide fasta files with fna or fna.gz or .fasta ending to faa
        Warning: if the directory path includes parentheses function prodigal is not working
    """
    zipFnaFiles = myUtil.compareFileLists(directory,".fna.gz",".faa.gz") 
    FnaFiles = myUtil.compareFileLists(directory,".fna",".faa")
    fastaFiles = myUtil.getAllFiles(directory,".fasta")
    NucleotideFastaFiles = zipFnaFiles + FnaFiles + fastaFiles
    print(f"[INFO] Found {len(NucleotideFastaFiles)} assemblies in nucleotide or ambigous format for prodigal")
    
    manager  = multiprocessing.Manager()
    counter = manager.Value('i',0)
    lock = manager.Lock()
    length = len(NucleotideFastaFiles)
    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length ,counter,lock) for fasta in NucleotideFastaFiles]
        pool.map(translate_fasta, args_list)

    return



def translate_fasta(args):
    fasta,length, counter, lock = args

    #unpack if required
    if os.path.splitext(fasta)[-1] == ".gz":
        fasta = myUtil.unpackgz(fasta)

    #Run prodigal
    output = os.path.splitext(fasta)[0]
    faa = output + ".faa"
    
    prodigal = myUtil.find_executable("prodigal")
    
    string = f"{prodigal} -a {faa} -i {fasta} >/dev/null 2>&1"
    try:
        os.system(string)
    except Exception as e:
        with lock:
            print(f"[WARN] Cannot translate {fasta} - {e}")
        
    with lock:
        counter.value += 1
        print(f"\r[INFO] Processing assembly {counter.value} of {length}", end ='',flush=True)        
    return




############################################################################
############### Parallel Transcription #####################################
############################################################################


def parallel_transcription(directory: str, cores: int) -> None:
    """
    Transcribe all .faa files to .gff files in parallel, ensuring both
    compressed and uncompressed versions are handled. Removes uncompressed
    files after processing.

    Warning: If the directory path contains parentheses, Prodigal may fail.

    Args:
        directory: Path to the directory containing FASTA/GFF files.
        cores: Number of parallel processes to use.
    """
    # Find all compressed FASTA and GFF files
    gz_faa_files: List[str] = myUtil.getAllFiles(directory, ".faa.gz")
    gz_gff_files: List[str] = myUtil.getAllFiles(directory, ".gff.gz")
    print(f"[INFO] Found {len(gz_faa_files)} zipped .faa files")
    print(f"[INFO] Found {len(gz_gff_files)} zipped .gff files")

    # Unpack all .gff.gz files in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker, gz_gff_files)

    # Unpack all .faa.gz files in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker, gz_faa_files)

    # Find all uncompressed .faa and .gff files
    faa_files: List[str] = myUtil.getAllFiles(directory, ".faa")
    gff_files: List[str] = myUtil.getAllFiles(directory, ".gff")
    print(f"[INFO] Found {len(gff_files)} .gff files")
    print(f"[INFO] Found {len(faa_files)} .faa files")

    # Identify .faa files without matching .gff
    unmatched_faa: List[str] = myUtil.compareFileLists(directory, ".faa", ".gff")
    print(f"[INFO] Found {len(unmatched_faa)} protein FASTA files without GFF")

    # Prepare for parallel transcription
    manager: Any = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()
    total = len(unmatched_faa)

    # Create argument tuples for each FASTA file
    args_list: List[Tuple[str, int, Any, Any]] = [
        (fasta_path, total, counter, lock)
        for fasta_path in unmatched_faa
    ]

    # Transcribe all unmatched .faa files in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(_transcribe_fasta, args_list)

    print("\n[INFO] Finished FASTA and GFF file preparation")
    return


def _transcribe_fasta(args):
    fasta, length, counter, lock = args
    gff = ""
            
    if _check_prodigal_format(fasta):        
        gff = _prodigal_faa_to_gff(fasta)
    
        with lock:
            counter.value += 1
            print(f"\rProcessing file {counter.value} of {length}", end='', flush=True)
  
    return

def _check_prodigal_format(filepath: str) -> bool:
    """
    Checks whether a given FASTA file uses the expected Prodigal header format.

    A valid Prodigal header (after the leading '>') should split into exactly
    five parts when using '#' as a delimiter:
        >contig#start#end#strand_indicator#other_info

    Args:
        filepath: Path to the FASTA file to check.

    Returns:
        True if the first header line splits into exactly five segments;
        False otherwise (including if no header is found).

    Raises:
        OSError: If the file cannot be opened for reading.
    """
    try:
        with open(filepath, "r") as infile:
            for raw_line in infile:
                if not raw_line.startswith(">"):
                    continue

                # Remove leading '>' and trailing newline/whitespace
                header = raw_line[1:].strip()
                parts = header.split("#")

                return len(parts) == 5
    except OSError as exc:
        raise OSError(f"Cannot open file '{filepath}': {exc}") from exc

    # No header lines were found
    return False


def _prodigal_faa_to_gff(filepath: str) -> str:
    """
    Translate a Prodigal-generated .faa file into GFF3 format.

    Parses each header line (starting with '>') assuming a Prodigal-specific format:
      ><contig>#<start>#<end>#<strand_indicator>#...
    Extracts contig, start, end, strand, and genome ID, then writes a GFF3 entry
    for each CDS.

    Args:
        filepath: Path to the Prodigal .faa file (optionally gzipped).

    Returns:
        The path to the generated .gff file.

    Raises:
        OSError: If the input file cannot be read or the output file cannot be written.
    """
    directory = os.path.dirname(filepath)
    filename = os.path.basename(filepath)

    # Remove ".gz" extension if present
    if filename.endswith(".gz"):
        filename = os.path.splitext(filename)[0]

    # Remove ".faa" extension if present
    if filename.endswith(".faa"):
        filename_without_ext = os.path.splitext(filename)[0]
    else:
        filename_without_ext = filename

    output_gff = os.path.join(directory, f"{filename_without_ext}.gff")
    genome_id = myUtil.getGenomeID(filepath)

    # Compile a regex to strip "_<digits><non-word>" at end of contig name
    contig_pattern = re.compile(r"_\d+\W+$")

    try:
        with open(filepath, "r") as infile, open(output_gff, "w") as outfile:
            for line in infile:
                if not line.startswith(">"):
                    continue

                header = line[1:].strip()
                try:
                    parts = header.split("#")
                    # parts example: [contig_raw, start, end, strand_indicator, ...]
                    contig_raw = parts[0]
                    contig = contig_pattern.split(contig_raw)[0]
                    start = parts[1]
                    end = parts[2]
                    strand = "+" if parts[3].strip() == "1" else "-"
                    gene_id = parts[0]

                    gff_line = (
                        f"{contig}\tprodigal\tcds\t{start}\t{end}\t0.0\t"
                        f"{strand}\t0\tID=cds-{gene_id};Genome={genome_id}\n"
                    )
                    outfile.write(gff_line)
                except IndexError as exc:
                    print(f"Error: Malformed header '{header}' - {exc}")
    except OSError as exc:
        raise OSError(f"Cannot process file '{filepath}': {exc}") from exc

    return output_gff


def unpacker(file):
    myUtil.unpackgz(file)    

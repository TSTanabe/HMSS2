import os
import pickle
import sys
import gzip
import shutil
import random
import logging
from datetime import datetime
from typing import Optional, Union, List, Dict, Set, Any
from pathlib import Path

logger = logging.getLogger(__name__)

def setup_logging(verbose_level: int = 0, log_file: Optional[str] = None) -> None:
    """
    Sets up logging:
      - Console: Level as specified (0=WARNING, 1=INFO, 2=DEBUG)
      - File (if provided): DEBUG level, full detail
      - Both: YYYY-MM-DD HH:MM:SS timestamp (to the second)
      - All logs go to both handlers (if file specified)

    Args:
        verbose_level (int): 0=WARNING, 1=INFO, 2=DEBUG
        log_file (str, optional): File path to write logs
    """
    # Determine user console log level
    level = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG
    }.get(verbose_level, logging.DEBUG)

    # Remove all handlers if present (for repeated calls in notebooks/scripts)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Console Handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_format = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_handler.setFormatter(console_format)
    handlers = [console_handler]

    # File Handler (if requested)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  # Always maximum detail
        file_format = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(processName)s | %(message)s | %(filename)s:%(lineno)d",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        file_handler.setFormatter(file_format)
        handlers.append(file_handler)

    # Root logger: accept everything, handlers do the filtering
    logging.basicConfig(
        level=logging.DEBUG,
        handlers=handlers
    )


def packgz(path: str) -> str:
    """
    Compress a file using gzip.

    Args:
        path (str): Path to input file.

    Returns:
        str: Path to compressed .gz file.
    """
    file = path + '.gz'
    with open(path, 'rb') as src, gzip.open(file, 'wb') as dst:
        dst.writelines(src)
    return file

def unpackgz(path: str) -> str:
    """
    Decompresses a .gz file if not already extracted.

    Args:
        path (str): Path to .gz file.

    Returns:
        str: Path to decompressed file.
    """
    if not path.endswith('.gz'):
        return path
    file = path[:-3]
    if os.path.exists(file):
        return file
    with gzip.open(path, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file

def dir_path(string: str) -> str:
    """Validates that a string points to a valid directory."""
    if os.path.isdir(string):
        return string.rstrip('/')
    sys.exit(f"\nERROR: {string} is not a valid directory")

def file_path(string: str) -> str:
    """Validates that a string points to a valid file."""
    if os.path.isfile(string):
        return string
    sys.exit(f"\nERROR: {string} is not a valid file")

def auto_float(value):
    if value == 'auto':
        return value
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("%r is not a valid float or 'auto'" % value)



def print_header(string: str, verbose: int = 0) -> None:
    """Prints a timestamped header if verbose is set."""
    if not verbose:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("\n" + string)
        print(current_time)
        print((len(current_time) + 3 + len(string)) * "-")

def clean_empty_files(directory: str) -> None:
    """Removes all empty files in a given directory."""
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            os.remove(file_path)
            logger.debug(f"Removed empty file: {file_path}")

def generate_color(seed_int: int) -> str:
    """Generates a consistent random color hex string from a seed integer."""
    random.seed(seed_int)
    return '#{:06x}'.format(random.randint(0, 0xFFFFFF))




def compare_file_lists(directory: str, ext1: str, ext2: str) -> Set[str]:
    """
    Compares files in a directory: finds files with ext1 that do not have a corresponding file with ext2.
    
    Returns:
        Set[str]: Set of filepaths (with ext1) that lack a matching file with ext2.
    """
    directory = Path(directory)
    files1 = {f.stem: f for f in directory.glob(f'*{ext1}') if f.is_file()}
    files2 = {f.stem: f for f in directory.glob(f'*{ext2}') if f.is_file()}
    
    # Set comprehension to collect missing file paths
    missing = {str(files1[stem]) for stem in files1 if stem not in files2}
    return missing



def getAllFiles(directory: str, ending: Union[str, int] = 0) -> List[str]:
    """
    Recursively retrieves files from a directory.

    Args:
        directory (str): Root directory to search.
        ending (str or int): Optional file suffix to filter by.

    Returns:
        List[str]: List of matching file paths.
    """
    matched = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            file = os.path.join(path, name)
            if ending == 0 or file.endswith(ending):
                matched.append(file)
    return matched

def compareFileLists(directory: str, ext1: Union[str, int] = 0, ext2: Union[str, int] = 0) -> List[str]:
    """
    Compares file lists by extension and returns those missing a match.

    Example:
        ext1 = '.faa', ext2 = '.txt'
        => returns ['genome1.faa', 'genome2.faa'] if .txt files are missing.
    """
    if ext1 and ext2:
        files1 = getAllFiles(directory, ext1)
        files2 = getAllFiles(directory, ext2)
        compare1 = removeExtFromList(files1, ext1)
        compare2 = removeExtFromList(files2, ext2)
        difference = set(compare1).difference(set(compare2))
        return addExtToList(list(difference), ext1)
    return []

def removeExtFromList(listing: List[str], ext: str) -> List[str]:
    """Removes extension and returns only basenames without directory."""
    return [os.path.splitext(os.path.basename(element))[0] for element in listing]

def addExtToList(listing: List[str], ext: str) -> List[str]:
    """Appends the specified extension to each string in the list."""
    return [element + ext for element in listing]

def getGenomeID(path: str) -> str:
    """Extracts genome ID from filename (before first dot)."""
    return os.path.basename(path).split('.')[0]

def getReportName(path: str) -> str:
    """Appends '.HmmReport' extension to a given path (removes original extension)."""
    return os.path.splitext(path)[0] + ".HmmReport"

def taxonomy_lineage(array: List[str], trennzeichen: str) -> str:
    """
    Joins taxonomy names with separator, replacing spaces with dashes.

    Returns:
        str: A single string like 'Bacteria-Firmicutes-Bacilli'.
    """
    try:
        string = trennzeichen.join(array).replace(" ", "-")
        return str(string)
    except Exception:
        return "NoTaxonomy"


def get_executable_dir() -> str:
    """Returns directory of script (or executable in case of PyInstaller)."""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def find_executable(executable: str) -> str:
    """
    Attempts to locate an executable in PATH or local ./bin directory.

    Raises:
        FileNotFoundError if the binary is not found.
    """
    path = shutil.which(executable)
    if path:
        logger.debug(f"Found executable in PATH: {path}")
        return path
    local_path = os.path.join(get_executable_dir(), "bin", executable)
    if os.path.isfile(local_path) and os.access(local_path, os.X_OK):
        logger.debug(f"Found local executable: {local_path}")
        return local_path
    logger.error(f"Executable not found: {executable}")
    raise FileNotFoundError(f"{executable} executable not found.")


def remove_directory(directory_path: str) -> None:
    """
    Deletes a directory and all its contents recursively.

    Logs a message whether the directory was deleted or not found.
    """
    if os.path.exists(directory_path):
        for root, dirs, files in os.walk(directory_path, topdown=False):
            for file in files:
                os.remove(os.path.join(root, file))
            for dir in dirs:
                os.rmdir(os.path.join(root, dir))
        os.rmdir(directory_path)
        logger.debug(f"Removed directory and contents: {directory_path}")
    else:
        logger.debug(f"Directory not found: {directory_path}")



def save_cache(options: Any, name: str, data: object, redirect: Optional[str] = None, overwrite: bool = False) -> None:
    """
    Saves a Python object using pickle into a defined cache directory.

    Args:
        options: Options object with a .result_files_directory attribute.
        name (str): Filename for the pickle file.
        data (object): Any serializable Python object.
        redirect (str, optional): Alternate output directory.
        overwrite (bool): Whether to overwrite an existing file.
    """
    cache_dir = redirect if redirect else os.path.join(options.result_files_directory, "pkl_cache")
    os.makedirs(cache_dir, exist_ok=True)
    file_path = os.path.join(cache_dir, name)
    if os.path.exists(file_path) and not overwrite:
        logger.debug(f"Cache file exists and overwrite disabled: {file_path}")
        return
    with open(file_path, "wb") as f:
        pickle.dump(data, f)
    logger.debug(f"Saved cache file: {file_path}")



def load_cache(options: Any, name: str, file_path: Optional[str] = None) -> Optional[object]:
    """
    Loads a Python object from a pickle file in the cache directory.

    Args:
        options: Options object with .result_files_directory.
        name (str): Name of the cache file.
        file_path (str, optional): Full path override.

    Returns:
        The loaded Python object or None if not found or on error.
    """
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache")
    file_path = file_path or os.path.join(cache_dir, name)
    if os.path.exists(file_path):
        try:
            logger.debug(f"Loading cache file: {file_path}")
            with open(file_path, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Failed to load cache {file_path}: {e}")
            return None
    logger.debug(f"Precomputed cache file not found: {file_path}")
    return None



def merge_grouped_refseq_dicts_simple(grouped_3_dict: Dict[str, Set[str]], grouped_4_dict: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    """
    Merges two protein dictionaries.

    Input example:
        grouped_3_dict = {"domainA": {"prot1", "prot2"}}
        grouped_4_dict = {"domainA": {"prot3"}}

    Output:
        {"domainA": {"prot1", "prot2", "prot3"}}
    """
    merged = {}
    all_domains = set(grouped_3_dict.keys()) | set(grouped_4_dict.keys())
    for domain in all_domains:
        merged[domain] = grouped_3_dict.get(domain, set()) | grouped_4_dict.get(domain, set())
    return merged



def merge_score_limits(dict1: Dict[str, Dict[str, float]], dict2: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """
    Merges two dictionaries of score limits.

    For each domain, computes the minimum lower_limit and maximum upper_limit.

    Input:
        dict1 = {"domainA": {"lower_limit": 20.0, "upper_limit": 50.0}}
        dict2 = {"domainA": {"lower_limit": 15.0, "upper_limit": 55.0}}

    Output:
        {"domainA": {"lower_limit": 15.0, "upper_limit": 55.0}}
    """
    merged = {}
    all_keys = set(dict1.keys()).union(dict2.keys())
    for key in all_keys:
        val1 = dict1.get(key, {})
        val2 = dict2.get(key, {})
        lower = min(val1.get("lower_limit", float('inf')), val2.get("lower_limit", float('inf')))
        upper = max(val1.get("upper_limit", float('-inf')), val2.get("upper_limit", float('-inf')))
        merged[key] = {"lower_limit": lower, "upper_limit": upper}
    return merged


def move_HMMs(input_folder,output_folder,file_extension):
    logger.info(f"Saving HMMs in the directory: {output_folder}")
    for datafile in os.listdir(input_folder):
        if datafile.endswith(file_extension):
            source = os.path.join(input_folder, datafile)
            target = os.path.join(output_folder, datafile)
            shutil.move(source,target)

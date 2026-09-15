#!/usr/bin/env python3

"""
Metadata-driven HMM library construction for HMSS3.

The metadata table determines:
    - available packages
    - available modules
    - default module/package combinations
    - HMMs belonging to each module/package combination

The CLI/config determines:
    - additional modules
    - replacement modules
    - metadata file
    - resource data directory
    - output library

No package names, HMM prefixes or biological module names are hard-coded here.

HMM selection uses filenames only.

Examples:
    grp3_redDsrA.hmm       -> redDsrA
    grpC_redDsrA.hmm       -> redDsrA
    grpCFirmi_redDsrA.hmm  -> redDsrA
    redDsrA.hmm            -> redDsrA

Everything before the first "_" in the HMM filename is ignored.

Package directories are identified by the part before the first "_":

    v8_Dsr_sHdr_Sox        -> v8
    DiSuCy_Dsr_sHdr        -> DiSuCy
    HMSS2_Dsr_DMS          -> HMSS2
"""

import csv
import hashlib
import os
import tempfile
from collections import defaultdict
from pathlib import Path

from hmsss.core.logging import get_logger

log = get_logger(__name__)


def load_metadata(path):
    """
    Read the module metadata TSV.

    Required columns:
        package
        module
        default
        hmms

    'default' must be exactly 'true' or 'false' (case-insensitive).

    All additional columns remain in the row dictionaries and can be used
    elsewhere in HMSS3 without this module having to know their meaning.
    """

    rows = []

    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Resource metadata has no header: {path}")

        required = {"package", "module", "default", "hmms"}
        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "Missing required resource metadata columns: "
                + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in row.items()}

            default = row["default"].lower()

            if default not in {"true", "false"}:
                raise ValueError(
                    f"Invalid default value in metadata line {line_number}: "
                    f"{row['default']!r}"
                )

            hmms = tuple(
                hmm.strip()
                for hmm in row["hmms"].replace(";", ",").split(",")
                if hmm.strip()
            )

            if not row["package"] or not row["module"] or not hmms:
                raise ValueError(
                    f"Incomplete resource definition in metadata line {line_number}"
                )

            row["default"] = default == "true"
            row["hmms"] = hmms
            row["_line"] = line_number

            rows.append(row)

    return rows


def _normalize_module_arguments(values):
    """
    Support both:

        ["redDsr@DiSCo", "Mcc@HMSS2"]

    and:

        ["redDsr@DiSCo,Mcc@HMSS2"]
    """

    return [
        item.strip()
        for value in (values or [])
        for item in value.split(",")
        if item.strip()
    ]


def _build_module_indices(metadata):
    """
    Build small lookup dictionaries once.

    by_selector:
        ("redDsr", "DiSCo") -> metadata row

    by_module:
        "redDsr" -> [all package variants]
    """

    by_selector = {}
    by_module = defaultdict(list)

    for row in metadata:
        by_selector[(row["module"], row["package"])] = row
        by_module[row["module"]].append(row)

    return by_selector, by_module


def _resolve_module(selector, by_selector, by_module):
    """
    Resolve:

        redDsr@DiSCo

    A package-less selector is also accepted when the module exists in only
    one package.
    """

    if "@" in selector:
        module, package = (part.strip() for part in selector.rsplit("@", 1))

        try:
            return by_selector[(module, package)]
        except KeyError:
            variants = ", ".join(
                f"{row['module']}@{row['package']}" for row in by_module.get(module, [])
            )

            if variants:
                raise ValueError(
                    f"Unknown module selection {selector!r}. "
                    f"Available variants: {variants}"
                ) from None

            raise ValueError(f"Unknown module: {module!r}") from None

    module = selector.strip()
    variants = by_module.get(module, [])

    if not variants:
        raise ValueError(f"Unknown module: {module!r}")

    if len(variants) != 1:
        choices = ", ".join(f"{row['module']}@{row['package']}" for row in variants)

        raise ValueError(
            f"Module {module!r} exists in multiple packages. "
            f"Specify one explicitly: {choices}"
        )

    return variants[0]


def select_modules(metadata, add_modules=None, replace_modules=None):
    """
    Select final module/package combinations.

    1. Start with every row where default=true.

    2. --add-module redDsr@DiSCo
       adds this variant without removing redDsr@v8.

    3. --replace-module redDsr@DiSCo
       removes every currently selected redDsr variant and selects DiSCo.

    Returns:
        selected metadata rows in original metadata order
        bool indicating whether the resulting selection equals the default
    """

    by_selector, by_module = _build_module_indices(metadata)

    default_keys = {
        (row["module"], row["package"]) for row in metadata if row["default"]
    }

    selected_keys = set(default_keys)

    for selector in _normalize_module_arguments(add_modules):
        row = _resolve_module(selector, by_selector, by_module)
        selected_keys.add((row["module"], row["package"]))

    for selector in _normalize_module_arguments(replace_modules):
        row = _resolve_module(selector, by_selector, by_module)

        selected_keys = {key for key in selected_keys if key[0] != row["module"]}

        selected_keys.add((row["module"], row["package"]))

    selected = [
        row for row in metadata if (row["module"], row["package"]) in selected_keys
    ]

    return selected, selected_keys == default_keys


def _find_package_directories(data_dir):
    """
    Index immediate subdirectories of DATA_DIR.

    Package = directory name before the first "_".

        v8_Dsr_Sox    -> v8
        DiSCo_Dsr     -> DiSCo

    Multiple directories may belong to the same package.
    """

    directories = defaultdict(list)

    for path in Path(data_dir).iterdir():
        if path.is_dir():
            package = path.name.split("_", 1)[0].lower()
            directories[package].append(path)

    return directories


def _hmm_name_from_filename(filename):
    """
    Convert concrete HMM filename to metadata HMM name.

        grp3_redDsrA.hmm      -> redDsrA
        grpC_redDsrA.hmm      -> redDsrA
        redDsrA.hmm           -> redDsrA
    """

    stem = Path(filename).stem
    return stem.split("_", 1)[1] if "_" in stem else stem


def find_selected_hmm_files(data_dir, selected_modules):
    """
    Translate selected metadata modules directly into source HMM files.

    No general HMM catalog is constructed.

    Only packages required by the selected modules are scanned and only HMMs
    requested by those modules are retained.

    Return:
        selected_files:
            list of dictionaries:
                package
                hmm
                path

        file_modules:
            Path -> set of module selectors using this HMM
    """

    required = defaultdict(set)

    for row in selected_modules:
        required[row["package"]].update(row["hmms"])

    package_dirs = _find_package_directories(data_dir)
    found = defaultdict(list)

    # Scan every required package exactly once.
    for package, required_hmms in required.items():
        directories = package_dirs.get(package.lower(), [])

        if not directories:
            raise ValueError(
                f"No resource directory found for package {package!r}. "
                f"Expected a top-level directory beginning with '{package}_'."
            )

        for package_dir in directories:
            for path in package_dir.rglob("*.hmm"):
                if not path.is_file():
                    continue

                hmm = _hmm_name_from_filename(path.name)

                if hmm in required_hmms:
                    found[(package, hmm)].append(path)

    missing = [
        (package, hmm)
        for package, hmms in required.items()
        for hmm in hmms
        if (package, hmm) not in found
    ]

    if missing:
        details = "\n".join(f"  {package}: {hmm}" for package, hmm in sorted(missing))

        raise ValueError(
            "The following HMMs requested by the resource metadata "
            f"were not found:\n{details}"
        )

    # Restore metadata-defined ordering.
    selected_files = []
    file_modules = defaultdict(set)
    seen_files = set()

    for row in selected_modules:
        selector = f"{row['module']}@{row['package']}"

        for hmm in row["hmms"]:
            for path in sorted(found[(row["package"], hmm)], key=str):
                file_modules[path].add(selector)

                if path in seen_files:
                    continue

                seen_files.add(path)

                selected_files.append(
                    {
                        "package": row["package"],
                        "hmm": hmm,
                        "path": path,
                    }
                )

    return selected_files, file_modules


def _hash_file(path, chunk_size=1024 * 1024):
    """Calculate SHA256 of one existing file."""

    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)

    return digest.digest()


def _hash_concatenated_sources(selected_files, chunk_size=1024 * 1024):
    """
    Calculate the SHA256 that the final concatenated library would have.

    No temporary library is created.

    Important:
        Source files are hashed exactly in the same order and without adding
        extra bytes between files. This mirrors write_library().
    """

    digest = hashlib.sha256()

    for entry in selected_files:
        with open(entry["path"], "rb") as handle:
            while chunk := handle.read(chunk_size):
                digest.update(chunk)

    return digest.digest()


def library_matches_selection(library, selected_files):
    """
    Check whether the existing library is byte-identical to the library that
    would be produced from selected_files.

    Fast path:
        compare total source size to existing library size first.

    Only if sizes match are SHA256 hashes calculated.

    Therefore no output or temporary file is created when the library is
    already current.
    """

    library = Path(library)

    if not library.is_file():
        return False

    expected_size = sum(entry["path"].stat().st_size for entry in selected_files)

    if library.stat().st_size != expected_size:
        return False

    return _hash_file(library) == _hash_concatenated_sources(selected_files)


def write_library(library, selected_files):
    """
    Atomically build the HMM library.

    A temporary file is only created when the current library is known to
    differ from the selected source HMMs.
    """

    library = Path(library)
    library.parent.mkdir(parents=True, exist_ok=True)

    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=library.parent,
            prefix=f".{library.name}.",
            delete=False,
        ) as output:
            temp_path = Path(output.name)

            for entry in selected_files:
                with open(entry["path"], "rb") as source:
                    while chunk := source.read(1024 * 1024):
                        output.write(chunk)

        os.replace(temp_path, library)

    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)

        raise


def build_library(config):
    """
    Build or reuse config.library based entirely on resource metadata and config.

    Required config attributes:
        config.resource_metadata
        config.resource_data_dir
        config.library

    Optional CLI/config attributes:
        config.add_modules
        config.replace_modules

    The function:

        1. reads metadata
        2. determines the module selection
        3. determines whether the selection equals metadata defaults
        4. finds the required HMM files
        5. compares the expected library to the existing config.library
        6. writes a new library only when it would actually differ

    Useful state is written back to config:
        config.library_is_default
        config.library_rebuilt
        config.selected_hmm_modules
        config.selected_hmm_files
    """

    metadata = load_metadata(config.resource_metadata)

    selected_modules, is_default = select_modules(
        metadata,
        add_modules=getattr(config, "add_modules", None),
        replace_modules=getattr(config, "replace_modules", None),
    )

    selected_files, _ = find_selected_hmm_files(
        config.paths.data,
        selected_modules,
    )

    selectors = [f"{row['module']}@{row['package']}" for row in selected_modules]

    config.library_is_default = is_default
    config.selected_hmm_modules = selectors
    config.selected_hmm_files = [str(entry["path"]) for entry in selected_files]

    if library_matches_selection(config.library, selected_files):
        config.library_rebuilt = False

        log.info(
            "Existing HMM library is identical to the selected resources; reusing %s",
            config.library,
        )

    else:
        write_library(config.library, selected_files)
        config.library_rebuilt = True

        log.info(
            "Built HMM library with %d HMM files: %s",
            len(selected_files),
            config.library,
        )

    if is_default:
        log.info("HMM library selection corresponds to metadata defaults")
    else:
        log.info(
            "Using custom HMM module selection: %s",
            ", ".join(selectors),
        )

    return {
        "is_default": is_default,
        "rebuilt": config.library_rebuilt,
        "modules": selectors,
        "hmm_count": len(selected_files),
    }

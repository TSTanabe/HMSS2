#!/usr/bin/env python3

"""
python3 check_resource_consistency.py --hmms data --cutoffs data/ --patterns data/Patterns --cooccurrence data/Cooccurrence --metadata data/Resource_module_metadata
Check consistency of HMSS3 resource files.

Two different HMM identities are deliberately used:

1. Library / metadata selection
   ----------------------------
   Uses the PACKAGE DIRECTORY and HMM FILE NAME.

   Package:
       DiSuCy_Dsr_Sox/...  -> DiSuCy

   HMM:
       grpC1_FccA.hmm      -> FccA
       grpC1_FccA1.hmm     -> FccA1
       grpCFirmi_FccA.hmm  -> FccA
       FccA.hmm            -> FccA

   Matching of the HMM part is exact and case-sensitive because this is
   also how library_generation.py selects HMM files.

   Therefore:

       metadata FccA != file FccA1.hmm
       metadata FccA != file fcca.hmm

   Both situations are reported.


2. Cutoffs / patterns / co-occurrence
   ----------------------------------
   These resources refer to the INTERNAL HMMER NAME.

   Example:

       file: grpC1_FccA1.hmm

       content:
           NAME  grpC_FccA

       cutoff:
           grpC_FccA ...

   The cutoff is compared to grpC_FccA, not to the filename.


The checker performs:

    physical HMM files <-> metadata
    internal HMM names <-> cutoffs
    internal HMM names <-> patterns
    internal HMM names <-> co-occurrence

Detailed individual errors are written to a TSV file.
Terminal output is aggregated.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.:+/\-]*$")
QUOTED_RE = re.compile(r"""['"]([^'"]+)['"]""")


# =============================================================================
# Reporter
# =============================================================================


class Reporter:
    """
    Detailed report:
        every individual occurrence is written to TSV.

    Terminal:
        repeated identical problems are aggregated.
    """

    def __init__(self, report_file=None):
        self.errors = 0
        self.warnings = 0
        self.summary = {}

        self.handle = None
        self.writer = None

        if report_file:
            self.handle = open(report_file, "w", encoding="utf-8", newline="")
            self.writer = csv.writer(self.handle, delimiter="\t", lineterminator="\n")

            self.writer.writerow(
                [
                    "severity",
                    "category",
                    "package",
                    "file",
                    "line",
                    "name",
                    "message",
                ]
            )

    def _record(
        self, severity, category, package="", filename="", line="", name="", message=""
    ):
        if self.writer:
            self.writer.writerow(
                [
                    severity,
                    category,
                    package,
                    str(filename),
                    line,
                    name,
                    message,
                ]
            )

        key = (severity, category, package, name, message)

        if key not in self.summary:
            self.summary[key] = {
                "count": 0,
                "locations": set(),
            }

        self.summary[key]["count"] += 1

        if filename:
            location = str(filename)
            if line:
                location += f":{line}"

            self.summary[key]["locations"].add(location)

    def error(self, category, filename="", line="", name="", message="", package=""):
        self.errors += 1
        self._record("ERROR", category, package, filename, line, name, message)

    def warning(self, category, filename="", line="", name="", message="", package=""):
        self.warnings += 1
        self._record("WARNING", category, package, filename, line, name, message)

    def print_problems(self):
        if not self.summary:
            return

        print()
        print("=" * 80)
        print("DETECTED PROBLEMS")
        print("=" * 80)

        grouped = defaultdict(list)

        for (severity, category, package, name, message), info in self.summary.items():
            grouped[(severity, category)].append((package, name, message, info))

        severity_order = {"ERROR": 0, "WARNING": 1}

        for severity, category in sorted(
            grouped,
            key=lambda value: (severity_order.get(value[0], 99), value[1]),
        ):
            print()
            print(f"[{severity}] {category}")

            entries = sorted(
                grouped[(severity, category)],
                key=lambda value: (value[0], value[1], value[2]),
            )

            for package, name, message, info in entries:
                package_label = package or "UNKNOWN"
                name_label = name or "<unnamed>"

                if info["count"] == 1:
                    print(f"  - [{package_label}] {name_label}")
                else:
                    print(
                        f"  - [{package_label}] {name_label} ({info['count']} occurrences)"
                    )

                # For these errors the concrete location is especially useful.
                if category in {
                    "duplicate_hmm",
                    "hmm_file_not_selectable",
                    "hmm_file_without_package",
                }:
                    for location in sorted(info["locations"]):
                        print(f"      -> {location}")

                # Naming mismatch hints are useful directly in the terminal.
                if category in {
                    "hmm_file_not_selectable",
                    "metadata_hmm_without_file",
                }:
                    print(f"      {message}")

        print()

    def unique_error_count(self):
        return sum(1 for key in self.summary if key[0] == "ERROR")

    def unique_warning_count(self):
        return sum(1 for key in self.summary if key[0] == "WARNING")

    def close(self):
        if self.handle:
            self.handle.close()


# =============================================================================
# HMM file discovery
# =============================================================================


def _iter_hmms_in_package(package_dir):
    """Recursively yield *.hmm files below one package directory."""

    for root, dirs, files in os.walk(package_dir):
        dirs[:] = [
            directory
            for directory in dirs
            if directory != "gpkg"
            and directory != "results"
            and not directory.endswith(".gpkg")
        ]

        for filename in files:
            # Keep this exact because library_generation.py uses *.hmm.
            if filename.endswith(".hmm"):
                yield Path(root) / filename


def iter_hmm_files(paths):
    """
    Yield (path, package).

    Expected normal use:

        --hmms /path/to/data

    with:

        data/
            v8_.../
            HMSS2_.../
            DiSuCy_.../
            DiSCo_.../

    Package = part of the package-directory name before the first "_".

    A package directory may also be supplied directly.
    """

    for raw_path in paths:
        path = Path(raw_path)

        if not path.exists():
            raise FileNotFoundError(f"Path does not exist: {path}")

        # Individual HMM file.
        if path.is_file():
            if path.name.endswith(".hmm"):
                parent = path.parent.name
                package = parent.split("_", 1)[0] if "_" in parent else None
                yield path, package

            continue

        # If the supplied directory itself looks like a package directory,
        # treat it as one package.
        if "_" in path.name:
            package = path.name.split("_", 1)[0]

            for hmm_path in _iter_hmms_in_package(path):
                yield hmm_path, package

            continue

        # Otherwise assume this is the resource root and direct children
        # are package directories.
        for package_dir in sorted(p for p in path.iterdir() if p.is_dir()):
            if package_dir.name in {"gpkg", "results"} or package_dir.name.endswith(
                ".gpkg"
            ):
                continue

            package = package_dir.name.split("_", 1)[0]

            for hmm_path in _iter_hmms_in_package(package_dir):
                yield hmm_path, package


# =============================================================================
# HMM filename identity used by library generation
# =============================================================================


def hmm_name_from_filename(path):
    """
    Derive the metadata HMM name exactly as library_generation.py does.

    Examples:

        grp3_redDsrA.hmm       -> redDsrA
        grpC1_FccA.hmm         -> FccA
        grpC1_FccA1.hmm        -> FccA1
        grpCFirmi_redDsrK.hmm  -> redDsrK
        FccA.hmm               -> FccA

    Everything before the first "_" is ignored.

    Important:
        Case and numerical suffixes after "_" are retained exactly.
    """

    stem = path.stem
    return stem.split("_", 1)[1] if "_" in stem else stem


# =============================================================================
# Read physical HMM inventory + internal HMMER NAMEs
# =============================================================================


def read_hmm_inventory(paths, reporter):
    """
    Read every source HMM once.

    Returns
    -------
    files:
        list of dictionaries:
            path
            package
            file_hmm
            internal_names

    internal_names:
        set of all HMMER NAME values

    internal_name_packages:
        internal HMM NAME -> packages containing that NAME

    n_files:
        number of physical HMM files

    File name and internal NAME are deliberately kept separate.
    """

    files = []
    internal_names = set()
    name_locations = defaultdict(list)
    internal_name_packages = defaultdict(set)

    for path, package in iter_hmm_files(paths):
        names_in_file = []

        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for lineno, line in enumerate(handle, 1):
                if not line.startswith("NAME"):
                    continue

                fields = line.split()

                if len(fields) < 2:
                    reporter.error(
                        "malformed_hmm_name",
                        path,
                        lineno,
                        "",
                        "Malformed HMM NAME line",
                        package=package or "UNKNOWN",
                    )
                    continue

                name = fields[1].strip()

                names_in_file.append(name)
                internal_names.add(name)
                name_locations[name].append((path, lineno))

                if package:
                    internal_name_packages[name].add(package)

        if not names_in_file:
            reporter.error(
                "invalid_hmm_file",
                path,
                "",
                path.name,
                "No HMMER NAME field found",
                package=package or "UNKNOWN",
            )

        files.append(
            {
                "path": path,
                "package": package,
                "file_hmm": hmm_name_from_filename(path),
                "internal_names": tuple(names_in_file),
            }
        )

    if not files:
        raise RuntimeError("No .hmm files found")

    # Report ALL locations of duplicated internal HMMER NAMEs.
    for name, locations in sorted(name_locations.items()):
        if len(locations) <= 1:
            continue

        packages = {
            package
            for path, _ in locations
            for package in [
                next(
                    (entry["package"] for entry in files if entry["path"] == path),
                    None,
                )
            ]
            if package
        }

        package_label = ",".join(sorted(packages)) if packages else "UNKNOWN"

        for path, lineno in locations:
            reporter.error(
                "duplicate_hmm",
                path,
                lineno,
                name,
                "Internal HMMER NAME occurs in more than one HMM location",
                package=package_label,
            )

    return files, internal_names, internal_name_packages, len(files)


def package_for_internal_name(name, internal_name_packages):
    """Return package label for an internal HMM name if it is known."""

    packages = internal_name_packages.get(name)

    if not packages:
        return "UNKNOWN"

    return ",".join(sorted(packages))


# =============================================================================
# Named text resource discovery
# =============================================================================


def iter_named_resource_files(paths, prefix, suffix=".txt"):
    """
    Yield files matching:

        <prefix>*<suffix>

    Examples:
        cutoffs*.txt
        patterns*.txt
        cooccurrence*.txt
    """

    prefix_lower = prefix.lower()
    suffix_lower = suffix.lower()

    for raw_path in paths:
        path = Path(raw_path)

        if not path.exists():
            raise FileNotFoundError(f"Path does not exist: {path}")

        if path.is_file():
            name_lower = path.name.lower()

            if name_lower.startswith(prefix_lower) and name_lower.endswith(
                suffix_lower
            ):
                yield path

            continue

        for root, dirs, files in os.walk(path):
            dirs[:] = [
                directory
                for directory in dirs
                if directory != "gpkg"
                and directory != "results"
                and not directory.endswith(".gpkg")
            ]

            for filename in files:
                filename_lower = filename.lower()

                if filename_lower.startswith(prefix_lower) and filename_lower.endswith(
                    suffix_lower
                ):
                    yield Path(root) / filename


# =============================================================================
# Cutoff parsing
# =============================================================================


def read_cutoff_names(paths, reporter, internal_name_packages):
    """
    Read cutoffs*.txt files.

    The first whitespace-delimited column is interpreted as the INTERNAL
    HMMER NAME.

    This deliberately does NOT use HMM filenames.
    """

    names = set()
    locations = defaultdict(list)
    n_files = 0

    for path in iter_named_resource_files(paths, prefix="cutoffs", suffix=".txt"):
        n_files += 1

        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for lineno, line in enumerate(handle, 1):
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                fields = line.split()

                if not fields:
                    continue

                name = fields[0]

                if name.lower() in {"hmm", "hmm_name", "domain", "name", "model"}:
                    continue

                package = package_for_internal_name(name, internal_name_packages)

                if not NAME_RE.match(name):
                    reporter.warning(
                        "invalid_cutoff_line",
                        path,
                        lineno,
                        name,
                        "First column does not look like an HMM name",
                        package=package,
                    )
                    continue

                names.add(name)
                locations[name].append((path, lineno))

    if n_files == 0:
        raise RuntimeError("No cutoffs*.txt files found")

    for name, occurrences in sorted(locations.items()):
        if len(occurrences) <= 1:
            continue

        package = package_for_internal_name(name, internal_name_packages)

        for path, lineno in occurrences:
            reporter.error(
                "duplicate_cutoff",
                path,
                lineno,
                name,
                "Cutoff occurs more than once",
                package=package,
            )

    return names, n_files


# =============================================================================
# Internal HMM NAME <-> cutoff consistency
# =============================================================================


def compare_hmms_cutoffs(hmm_names, cutoff_names, internal_name_packages, reporter):
    """
    Compare INTERNAL HMMER NAMEs to cutoff names.

    This is independent of metadata/library filename selection.
    """

    for name in sorted(hmm_names - cutoff_names):
        reporter.error(
            "hmm_without_cutoff",
            "",
            "",
            name,
            "Internal HMMER NAME has no cutoff",
            package=package_for_internal_name(name, internal_name_packages),
        )

    for name in sorted(cutoff_names - hmm_names):
        reporter.error(
            "cutoff_without_hmm",
            "",
            "",
            name,
            "Cutoff has no corresponding internal HMMER NAME",
            package="UNKNOWN",
        )


# =============================================================================
# Metadata
# =============================================================================


def read_metadata(path, reporter):
    """
    Read all potentially selectable HMMs from metadata.

    Required columns:
        package
        module
        default
        hmms

    Every metadata row is considered, independent of default=true/false.

    Metadata key:
        (package.casefold(), exact HMM name)

    Package matching is case-insensitive.
    HMM matching is intentionally case-sensitive.
    """

    metadata_keys = set()
    metadata_info = {}
    package_names = {}
    module_package_rows = defaultdict(list)

    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise RuntimeError("Metadata has no header")

        required = {"package", "module", "default", "hmms"}
        missing = required - set(reader.fieldnames)

        if missing:
            raise RuntimeError(
                "Missing metadata columns: " + ", ".join(sorted(missing))
            )

        for lineno, row in enumerate(reader, 2):
            package = (row.get("package") or "").strip()
            module = (row.get("module") or "").strip()
            default = (row.get("default") or "").strip().lower()
            hmm_field = (row.get("hmms") or "").strip()

            if not package:
                reporter.error(
                    "empty_package",
                    path,
                    lineno,
                    "",
                    "Package is empty",
                    package="UNKNOWN",
                )
                continue

            if not module:
                reporter.error(
                    "empty_module",
                    path,
                    lineno,
                    "",
                    "Module name is empty",
                    package=package,
                )

            if default not in {"true", "false"}:
                reporter.error(
                    "invalid_default",
                    path,
                    lineno,
                    module,
                    "default must be true or false",
                    package=package,
                )

            package_key = package.casefold()
            package_names.setdefault(package_key, package)

            module_package_rows[(package_key, module)].append(lineno)

            for hmm in hmm_field.replace(";", ",").split(","):
                hmm = hmm.strip()

                if not hmm:
                    continue

                key = (package_key, hmm)
                metadata_keys.add(key)

                info = metadata_info.setdefault(
                    key,
                    {
                        "package": package,
                        "hmm": hmm,
                        "modules": set(),
                        "lines": [],
                    },
                )

                info["modules"].add(module)
                info["lines"].append(lineno)

    # Same package/module combination is ambiguous for CLI module selection.
    for (package_key, module), lines in sorted(module_package_rows.items()):
        if len(lines) <= 1:
            continue

        package = package_names.get(package_key, package_key)

        reporter.error(
            "duplicate_module_package",
            path,
            lines[0],
            module,
            f"Module/package combination occurs in metadata lines {','.join(map(str, lines))}",
            package=package,
        )

    return metadata_keys, metadata_info, package_names


# =============================================================================
# Similar-name diagnostics
# =============================================================================


def find_similar_names(name, candidates):
    """
    Find likely naming inconsistencies.

    These are diagnostic hints ONLY and never count as valid matches.

    Examples:
        FccA  vs fcca   -> capitalization differs
        FccA  vs FccA1  -> numerical suffix differs

    Exact matching remains required.
    """

    similar = []

    for candidate in sorted(candidates):
        if candidate == name:
            continue

        if candidate.casefold() == name.casefold():
            similar.append(f"{candidate} [capitalization differs]")
            continue

        stripped_name = re.sub(r"\d+$", "", name)
        stripped_candidate = re.sub(r"\d+$", "", candidate)

        if stripped_name.casefold() == stripped_candidate.casefold():
            similar.append(f"{candidate} [numerical suffix/capitalization differs]")

    return similar


# =============================================================================
# Physical HMM files <-> metadata consistency
# =============================================================================


def compare_hmm_files_metadata(hmm_files, metadata_keys, metadata_info, reporter):
    """
    Bidirectional library-selectability test.

    A) Every physical *.hmm file must be represented in metadata.

       File identity:
           package from package directory
           HMM name from filename after first "_"

    B) Every metadata (package, HMM) must have at least one physical *.hmm file.

    Matching:
        package -> case-insensitive
        HMM     -> exact and case-sensitive

    Returns
    -------
    selectable_paths:
        physical HMM files that can actually be selected via metadata

    selectable_internal_names:
        INTERNAL HMMER NAMEs contained in metadata-selectable files.
        This is used for pattern/co-occurrence consistency.
    """

    filesystem = defaultdict(list)
    filesystem_package_names = {}
    file_by_path = {}

    for entry in hmm_files:
        path = entry["path"]
        package = entry["package"]
        file_by_path[path] = entry

        if not package:
            reporter.error(
                "hmm_file_without_package",
                path,
                "",
                entry["file_hmm"],
                "Could not determine package from package directory",
                package="UNKNOWN",
            )
            continue

        package_key = package.casefold()
        filesystem_package_names.setdefault(package_key, package)

        filesystem[(package_key, entry["file_hmm"])].append(entry)

    metadata_names_by_package = defaultdict(set)
    filesystem_names_by_package = defaultdict(set)

    for package_key, hmm in metadata_keys:
        metadata_names_by_package[package_key].add(hmm)

    for package_key, hmm in filesystem:
        filesystem_names_by_package[package_key].add(hmm)

    selectable_paths = set()

    # -------------------------------------------------------------------------
    # A) Every physical HMM file must be represented by metadata
    # -------------------------------------------------------------------------

    for (package_key, file_hmm), entries in sorted(
        filesystem.items(),
        key=lambda item: (item[0][0], item[0][1]),
    ):
        key = (package_key, file_hmm)

        if key in metadata_keys:
            selectable_paths.update(entry["path"] for entry in entries)
            continue

        package = filesystem_package_names.get(package_key, package_key)

        similar = find_similar_names(
            file_hmm,
            metadata_names_by_package.get(package_key, set()),
        )

        hint = ""
        if similar:
            hint = " Similar metadata name(s): " + ", ".join(similar) + "."

        for entry in entries:
            reporter.error(
                "hmm_file_not_selectable",
                entry["path"],
                "",
                file_hmm,
                (f"No exact metadata entry exists for {package}:{file_hmm}.{hint}"),
                package=package,
            )

    # -------------------------------------------------------------------------
    # B) Every metadata HMM must have at least one physical HMM file
    # -------------------------------------------------------------------------

    for package_key, metadata_hmm in sorted(metadata_keys):
        if (package_key, metadata_hmm) in filesystem:
            continue

        info = metadata_info[(package_key, metadata_hmm)]
        package = info["package"]

        similar = find_similar_names(
            metadata_hmm,
            filesystem_names_by_package.get(package_key, set()),
        )

        hint = ""
        if similar:
            hint = " Similar physical HMM name(s): " + ", ".join(similar) + "."

        modules = ",".join(sorted(info["modules"]))
        lines = ",".join(str(line) for line in info["lines"])

        reporter.error(
            "metadata_hmm_without_file",
            "",
            info["lines"][0] if info["lines"] else "",
            metadata_hmm,
            (
                f"Metadata entry {package}:{metadata_hmm} has no exact *.hmm file "
                f"in the corresponding package. Modules: {modules}; metadata lines: {lines}."
                f"{hint}"
            ),
            package=package,
        )

    # Internal names belonging to files which the library can actually select.
    selectable_internal_names = {
        internal_name
        for path in selectable_paths
        for internal_name in file_by_path[path]["internal_names"]
    }

    return selectable_paths, selectable_internal_names


# =============================================================================
# Pattern / co-occurrence parsing
# =============================================================================


def extract_names_from_field(field):
    """
    Extract potential INTERNAL HMMER NAMEs.

    Examples:

        grp3_SoxA

        ('grp3_SoxA', 'grp3_SoxB')

        ['grp3_SoxA', 'grp3_SoxB']

        grp3_SoxA,grp3_SoxB
    """

    field = field.strip()

    if not field:
        return

    quoted = QUOTED_RE.findall(field)

    if quoted:
        for value in quoted:
            value = value.strip()

            if NAME_RE.match(value):
                yield value

        return

    if NAME_RE.match(field):
        yield field
        return

    translation = str.maketrans(
        {
            "(": " ",
            ")": " ",
            "[": " ",
            "]": " ",
            "{": " ",
            "}": " ",
            ",": " ",
            ";": " ",
        }
    )

    cleaned = field.translate(translation)

    for token in cleaned.split():
        token = token.strip("'\"")

        if NAME_RE.match(token):
            yield token


def check_reference_files(
    paths,
    resource_prefix,
    label,
    hmm_names,
    selectable_hmms,
    internal_name_packages,
    reporter,
):
    """
    Check pattern/co-occurrence resources against INTERNAL HMMER NAMEs.

    First TAB-separated field is interpreted as block identifier and ignored.

    selectable_hmms contains internal HMMER NAMEs whose physical files can
    actually be selected through metadata.
    """

    n_files = 0
    n_references = 0

    for path in iter_named_resource_files(paths, prefix=resource_prefix, suffix=".txt"):
        n_files += 1

        with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
            for lineno, line in enumerate(handle, 1):
                line = line.rstrip("\n")

                if not line.strip() or line.lstrip().startswith("#"):
                    continue

                fields = line.split("\t")

                if len(fields) < 2:
                    reporter.warning(
                        f"{label}_format",
                        path,
                        lineno,
                        "",
                        "Line has fewer than two tab-separated fields",
                        package="UNKNOWN",
                    )
                    continue

                # First field = block ID.
                for field in fields[1:]:
                    for name in extract_names_from_field(field):
                        # Preserve the previous conservative behaviour:
                        # unknown bare words are ignored, while grp* values
                        # are interpreted as explicit HMM references.
                        if not name.startswith("grp") and name not in hmm_names:
                            continue

                        n_references += 1
                        package = package_for_internal_name(
                            name, internal_name_packages
                        )

                        if name not in hmm_names:
                            reporter.error(
                                f"{label}_hmm_not_found",
                                path,
                                lineno,
                                name,
                                "Referenced internal HMMER NAME does not exist",
                                package=package,
                            )
                            continue

                        if name not in selectable_hmms:
                            reporter.error(
                                f"{label}_hmm_not_selectable",
                                path,
                                lineno,
                                name,
                                (
                                    "Referenced internal HMMER NAME exists, but its "
                                    "physical HMM file is not selectable through metadata"
                                ),
                                package=package,
                            )

    return n_references, n_files


# =============================================================================
# CLI
# =============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Check consistency between HMSS3 HMM files, internal HMMER names, "
            "cutoffs, module metadata, patterns and co-occurrence resources."
        )
    )

    parser.add_argument(
        "--hmms",
        nargs="+",
        required=True,
        metavar="PATH",
        help=(
            "Resource root or package directories containing *.hmm files. "
            "Package names are derived from package-directory names."
        ),
    )

    parser.add_argument(
        "--cutoffs",
        nargs="+",
        required=True,
        metavar="PATH",
        help="Cutoff directories/files. Only cutoffs*.txt files are read.",
    )

    parser.add_argument(
        "--metadata",
        required=True,
        metavar="TSV",
        help="Resource module metadata TSV.",
    )

    parser.add_argument(
        "--patterns",
        nargs="*",
        default=[],
        metavar="PATH",
        help="Pattern directories/files. Only patterns*.txt files are read.",
    )

    parser.add_argument(
        "--cooccurrence",
        nargs="*",
        default=[],
        metavar="PATH",
        help="Co-occurrence directories/files. Only cooccurrence*.txt files are read.",
    )

    parser.add_argument(
        "--report",
        default="resource_consistency_report.tsv",
        metavar="TSV",
        help="Detailed TSV report retaining every individual occurrence.",
    )

    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================


def main():
    args = parse_args()
    reporter = Reporter(args.report)

    try:
        # ---------------------------------------------------------------------
        # Physical HMMs + internal HMMER NAMEs
        # ---------------------------------------------------------------------

        print("Reading HMM files...")

        (
            hmm_files,
            hmm_names,
            internal_name_packages,
            n_hmm_files,
        ) = read_hmm_inventory(
            args.hmms,
            reporter,
        )

        print(f"  physical HMM files: {n_hmm_files}")
        print(f"  unique internal HMMER NAMEs: {len(hmm_names)}")

        # ---------------------------------------------------------------------
        # Cutoffs
        # ---------------------------------------------------------------------

        print()
        print("Reading cutoffs...")

        cutoff_names, n_cutoff_files = read_cutoff_names(
            args.cutoffs,
            reporter,
            internal_name_packages,
        )

        print(f"  cutoff files: {n_cutoff_files}")
        print(f"  unique cutoff names: {len(cutoff_names)}")

        print()
        print("Checking internal HMMER NAME <-> cutoff...")

        compare_hmms_cutoffs(
            hmm_names,
            cutoff_names,
            internal_name_packages,
            reporter,
        )

        # ---------------------------------------------------------------------
        # Metadata
        # ---------------------------------------------------------------------

        print()
        print("Reading metadata...")

        (
            metadata_keys,
            metadata_info,
            metadata_package_names,
        ) = read_metadata(
            args.metadata,
            reporter,
        )

        print(f"  metadata package/HMM keys: {len(metadata_keys)}")

        print()
        print("Checking physical HMM files <-> metadata...")

        (
            selectable_paths,
            selectable_internal_hmms,
        ) = compare_hmm_files_metadata(
            hmm_files,
            metadata_keys,
            metadata_info,
            reporter,
        )

        print(
            f"  selectable physical HMM files: {len(selectable_paths)} / {n_hmm_files}"
        )
        print(
            f"  internal HMMER NAMEs in selectable files: "
            f"{len(selectable_internal_hmms)} / {len(hmm_names)}"
        )

        # ---------------------------------------------------------------------
        # Patterns
        # ---------------------------------------------------------------------

        pattern_refs = 0
        n_pattern_files = 0

        if args.patterns:
            print()
            print("Checking patterns...")

            pattern_refs, n_pattern_files = check_reference_files(
                args.patterns,
                resource_prefix="patterns",
                label="pattern",
                hmm_names=hmm_names,
                selectable_hmms=selectable_internal_hmms,
                internal_name_packages=internal_name_packages,
                reporter=reporter,
            )

            print(f"  pattern files: {n_pattern_files}")
            print(f"  HMM references: {pattern_refs}")

        # ---------------------------------------------------------------------
        # Co-occurrence
        # ---------------------------------------------------------------------

        cooccurrence_refs = 0
        n_cooccurrence_files = 0

        if args.cooccurrence:
            print()
            print("Checking co-occurrence...")

            cooccurrence_refs, n_cooccurrence_files = check_reference_files(
                args.cooccurrence,
                resource_prefix="cooccurrence",
                label="cooccurrence",
                hmm_names=hmm_names,
                selectable_hmms=selectable_internal_hmms,
                internal_name_packages=internal_name_packages,
                reporter=reporter,
            )

            print(f"  co-occurrence files: {n_cooccurrence_files}")
            print(f"  HMM references: {cooccurrence_refs}")

        # ---------------------------------------------------------------------
        # Problems
        # ---------------------------------------------------------------------

        reporter.print_problems()

        # ---------------------------------------------------------------------
        # Summary
        # ---------------------------------------------------------------------

        print()
        print("=" * 80)
        print("RESOURCE CONSISTENCY SUMMARY")
        print("=" * 80)

        print(f"Physical HMM files                  : {n_hmm_files}")
        print(f"Unique internal HMMER NAMEs         : {len(hmm_names)}")
        print(f"Cutoff files                        : {n_cutoff_files}")
        print(f"Unique cutoff names                 : {len(cutoff_names)}")
        print(f"Metadata package/HMM keys           : {len(metadata_keys)}")
        print(f"Metadata-selectable HMM files       : {len(selectable_paths)}")
        print(f"Selectable internal HMMER NAMEs     : {len(selectable_internal_hmms)}")
        print(f"Pattern files                       : {n_pattern_files}")
        print(f"Pattern references                  : {pattern_refs}")
        print(f"Co-occurrence files                 : {n_cooccurrence_files}")
        print(f"Co-occurrence references            : {cooccurrence_refs}")

        print("-" * 80)

        print(f"Error occurrences                   : {reporter.errors}")
        print(f"Unique errors                       : {reporter.unique_error_count()}")
        print(f"Warning occurrences                 : {reporter.warnings}")
        print(
            f"Unique warnings                     : {reporter.unique_warning_count()}"
        )

        print("=" * 80)
        print(f"Detailed report: {args.report}")

        if reporter.errors == 0:
            print()
            print("[OK] Resources are consistent.")
            return 0

        print()
        print("[FAILED] Resource inconsistencies detected.")
        return 1

    except Exception as error:
        print(f"[FATAL] {error}", file=sys.stderr)
        return 2

    finally:
        reporter.close()


if __name__ == "__main__":
    sys.exit(main())

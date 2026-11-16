#!/usr/bin/env python3
"""Build phylogenetic placement packages from a recursive tree of .gpkg folders.

Edit PROTEIN_GROUPS below. Each group becomes one output package.
For each group, the first protein's .gpkg supplies decoy_database.faa.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
import sys
from pathlib import Path

LOGGER = logging.getLogger("gpkg_package_builder")

# ============================================================
# EDIT THIS SECTION: proteins that should share one tree
# ============================================================
# Keys: new placement package directory names.
# Values: protein names, ordered! The FIRST protein provides the outgroup.
# These are illustrative names; replace or extend with your actual groups.
PROTEIN_GROUPS: dict[str, list[str]] = {
    "YihQ": ["YihQ"],
}
PROTEIN_GROUPS2: dict[str, list[str]] = {
    # "MarD_MorD_NflD": ["MarD", "MorD", "NflD-like"],
    # "AsrC_DsrAB": ["AsrC", "oxDsrA", "oxDsrB", "redDsrA", "redDsrB"],
    # "DsrM": ["QmoC", "redDsrM", "oxDsrM"],
    # "DsrK": ["redDsrK", "oxDsrK"],
    # "DsrJ": ["redDsrJ", "oxDsrJ"],
    "DsrO": ["redDsrO", "oxDsrO", "TtrB", "SoeB", "PhsB", "DmsB", "DdhB"],
    "DsrP": ["redDsrP", "oxDsrP", "TtrC", "SoeC", "PhsC", "DmsC", "DdhC"],
    "SoeA": ["TtrA", "SoeA", "PhsA", "DmsA", "DdhA"],
    "AprA": ["redAprA", "oxAprAI", "oxAprAII"],
    "AprB": ["redAprB", "oxAprBI", "oxAprBII"],
    "QmoA": ["QmoA"],
    "QmoB": ["QmoB", "sHdrB3"],
    "HdrA": ["sHdrA", "dHdrA"],
    "HdrB": ["qHdrB", "sHdrB1", "sHdrB2"],
    "HdrC": ["dHdrC", "qHdrC", "sHdrC1", "sHdrC2"],
    "Lbp": ["LbpA1", "LbpA2"],
    "LipS": ["LipS1", "LipS2"],
    "Sat": ["redSat", "oxSat"],
    "MccA": ["MccA"],
    "AsrA": ["ShyB", "AsrA"],
    "DsrC": ["redDsrC", "oxDsrC"],
    "DsrE": ["DsrE", "DsrE3A", "DsrE3B", "DsrE3C"],
    "SoxA": ["TsdA", "SoxA"],
    "SoxB": ["SoxB"],
    "SoxC": ["SoxC", "SorA"],
    "SoxD": ["SoxD"],
    "SoxX": ["SoxX"],
    "SoxY": ["SoxY"],
    "DoxD": ["DoxD"],
    "DoxA": ["DoxA"],
    "DddW": ["DddW", "DddK"],
    "DddD": ["DddD"],
    "DddA": ["DddA"],
    "DddX": ["DddX"],
    "DddL": ["DddL"],
    "DddC": ["DddC"],
    "DddP": ["DddP"],
    "DddQ": ["DddQ"],
    "DddT": ["DddT"],
    "DddX": ["DddX"],
    "DddX": ["DddY"],
    "Ssu": ["SnfG", "DmoA", "SsuD", "MsuD"],
    "MddA": ["MddA"],
    "MddH": ["MddH"],
    "DdhB": ["DdhB"],
    "DdhC": ["DdhC"],
    "DmdA": ["DmdA"],
    "DmdB": ["DmdB"],
    "DmdA": ["DmdA"],
    "DorA": ["DorA"],
    "DorC": ["DorC"],
    "DorD": ["DorD"],
    "DsoA": ["DsoA"],
    "DsoB": ["DsoB"],
    "DsoC": ["DsoC"],
    "DsoD": ["DsoD"],
    "YihQ": ["YihQ"],
    "SftI": ["SftI"],
    "YihS": ["YihS"],
    "YihT": ["YihT"],
    "YihU": ["YihU"],
    "YihV": ["YihV"],
    "SqwK": ["SqwK"],
    "SqwL": ["SqwK"],
    "SqwK": ["SqwL"],
    "SqwI": ["SqwI"],
    "SqwG": ["SqwG"],
    "SqwH": ["SqwH"],
    "SqwD": ["SqwD"],
    "SqwF": ["SqwF"],
    "SqvB": ["SqvB"],
    "SqoD": ["SqoD"],
    "SmoC": ["SmoC"],
    "SqiK": ["SqiK"],
    "SqiA": ["SqiA"],
    "SqgA": ["SqgA"],
    "Sqdh": ["Sqdh"],
    "Sgdh": ["Sgdh"],
    "Sgald": ["Sgald"],
    "Sgl": ["Sgl"],
    "SftD": ["SftD"],
    "SftT": ["SftT"],
    "Sgl": ["Sgl"],
    "SuyA": ["SuyA"],
    "SuyB": ["SuyB"],
    "HpsG": ["HpsG"],
    "HpsH": ["HpsH"],
    "HpfG": ["HpfG"],
    "HpfH": ["HpfH"],
    "HpfD": ["HpfD"],
    "HpsN": ["HpsN"],
    "HpsP": ["HpsP"],
    "HpfX": ["HpfX"],
    "HpfY": ["HpfY"],
    "HpfZ": ["HpfZ"],
    "IseJ": ["IseJ"],
    "IsfD": ["IsfD"],
    "SarD": ["SarD"],
    "SauS": ["SauS"],
    "SauT": ["SauT"],
    "SlcC": ["SlcC"],
    "SlcD": ["SlcD"],
    "TauX": ["TauX"],
    "TauY": ["TauY"],
    "Toa": ["Toa"],
    "Tpa": ["Tpa"],
    "ComC": ["ComC"],
    "ComD": ["ComD"],
    "ComE": ["ComE"],
    # "Other_family": ["ProteinA", "ProteinB"],
}

OUTGROUP_LABEL = "Outgroup"
GENERATOR_FILENAME = "generate_ppgk_alignment_workflows.py"


def index_gpkg_directories(root: Path) -> dict[str, list[Path]]:
    """Index folders with exact names <protein>.gpkg (case sensitive)."""
    index: dict[str, list[Path]] = {}
    for folder in sorted(root.rglob("*.gpkg")):
        if folder.is_dir() and folder.name.endswith(".gpkg"):
            protein = folder.name[:-5]
            index.setdefault(protein, []).append(folder)
    return index


def resolve_group(
    name: str,
    proteins: list[str],
    index: dict[str, list[Path]],
) -> tuple[list[tuple[Path, str]], Path]:
    """Resolve every protein to <protein>.gpkg/<protein>.faa; get first decoy."""
    if not proteins:
        raise ValueError(f"Group {name!r} is empty")
    if len(set(proteins)) != len(proteins):
        raise ValueError(f"Group {name!r} contains duplicate protein names")

    seeds: list[tuple[Path, str]] = []
    first_folder: Path | None = None
    for protein in proteins:
        folders = index.get(protein, [])
        if not folders:
            raise FileNotFoundError(f"{name}: no folder named {protein}.gpkg")
        if len(folders) > 1:
            paths = "\n  ".join(str(p) for p in folders)
            raise ValueError(
                f"{name}: ambiguous {protein}.gpkg; found {len(folders)} folders:\n  {paths}"
            )

        folder = folders[0]
        if first_folder is None:
            first_folder = folder
        # Search for seed FASTA (priority: protein name, then refseq_database)
        seed_candidates = [
            folder / f"{protein}.faa",
            folder / "refseq_database.faa",
        ]

        seed = next(
            (path for path in seed_candidates if path.is_file()),
            None,
        )

        if seed is None:
            raise FileNotFoundError(
                f"{name}: No seed FASTA found for {protein} in {folder}"
            )

        seeds.append((seed.resolve(), protein))

    assert first_folder is not None
    outgroup = first_folder / "decoy_database.faa"
    if not outgroup.is_file():
        raise FileNotFoundError(
            f"{name}: missing outgroup {outgroup} (from first protein {proteins[0]!r})"
        )
    return seeds, outgroup.resolve()


def make_command(
    generator: Path,
    seeds: list[tuple[Path, str]],
    outgroup: Path,
    output: Path,
    taxonomy: Path,
    workflow: str,
    threads: int,
    package_version: str,
) -> list[str]:
    command = [
        sys.executable,
        str(generator),
        "--taxonomy",
        str(taxonomy),
        "--output",
        str(output),
        "--alignment-workflow",
        workflow,
        "--threads",
        str(threads),
        "--package-version",
        package_version,
    ]
    for fasta, protein in seeds:
        command.extend(["--seed", str(fasta), protein])
    command.extend(["--outgroup", str(outgroup), OUTGROUP_LABEL])
    return command


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build placement packages from hard-coded protein groups and .gpkg folders."
    )
    parser.add_argument(
        "--gpkg-dir",
        required=True,
        type=Path,
        help="Root directory recursively containing <protein>.gpkg folders",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output root: each group gets its own subdirectory",
    )
    parser.add_argument(
        "--taxonomy",
        required=True,
        type=Path,
        help="Genome taxonomy TSV passed to placement-package generator",
    )
    parser.add_argument(
        "--alignment-workflow",
        choices=("A", "B", "D"),
        required=True,
        help="A: add then trim; B: add --keeplength then trim; "
        "D: trim ingroup then add --keeplength",
    )
    parser.add_argument(
        "--group-index",
        type=int,
        nargs="+",
        metavar="N",
        help="Run only the selected zero-based PROTEIN_GROUPS indices (e.g. 0 3 5); default: all groups.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--package-version", default="1.0.0")
    parser.add_argument(
        "--generator",
        type=Path,
        default=Path(__file__).resolve().with_name(GENERATOR_FILENAME),
        help="Path to generate_ppgk_alignment_workflows.py",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve inputs and show commands without executing",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Try subsequent groups after errors; still return nonzero",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow using an already existing output group folder",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    root = args.gpkg_dir.resolve()
    output_root = args.output_dir.resolve()
    taxonomy = args.taxonomy.resolve()
    generator = args.generator.resolve()
    if not root.is_dir():
        parser.error(f"GPKG root does not exist: {root}")
    if not taxonomy.is_file():
        parser.error(f"Taxonomy TSV does not exist: {taxonomy}")
    if not generator.is_file():
        parser.error(f"Generator script does not exist: {generator}")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if not PROTEIN_GROUPS:
        parser.error("PROTEIN_GROUPS is empty; configure it at the top of this script")

    # Python dictionaries preserve insertion order. Indices are zero-based.
    groups = list(PROTEIN_GROUPS.items())
    if args.group_index is not None:
        invalid = [i for i in args.group_index if not 0 <= i < len(groups)]
        if invalid:
            parser.error(
                f"Invalid group index(es): {invalid}; valid range is 0..{len(groups) - 1}"
            )
        selected_indices = list(dict.fromkeys(args.group_index))
    else:
        selected_indices = list(range(len(groups)))

    LOGGER.info("Selected %d of %d groups", len(selected_indices), len(groups))
    for idx in selected_indices:
        LOGGER.info("Group %d: %s", idx, groups[idx][0])

    index = index_gpkg_directories(root)
    LOGGER.info("Indexed %d distinct protein packages from %s", len(index), root)
    failed = 0
    succeeded = 0

    for idx in selected_indices:
        group_name, proteins = groups[idx]
        try:
            # Prevent traversal or directory collisions through editable group names.
            if group_name in {"", ".", ".."} or Path(group_name).name != group_name:
                raise ValueError(f"Invalid output group name: {group_name!r}")
            seeds, outgroup = resolve_group(group_name, proteins, index)
            output = output_root / group_name
            if output.exists() and not args.overwrite:
                raise FileExistsError(
                    f"Output directory already exists: {output}; "
                    "use --overwrite to reuse it"
                )
            command = make_command(
                generator,
                seeds,
                outgroup,
                output,
                taxonomy,
                args.alignment_workflow,
                args.threads,
                args.package_version,
            )
            LOGGER.info("[%s] Proteins: %s", group_name, ", ".join(proteins))
            LOGGER.info("[%s] Outgroup: %s", group_name, outgroup)
            LOGGER.info("[%s] Command: %s", group_name, shlex.join(command))
            if not args.dry_run:
                output_root.mkdir(parents=True, exist_ok=True)
                subprocess.run(command, check=True)
            succeeded += 1
        except (OSError, ValueError, subprocess.CalledProcessError) as exc:
            failed += 1
            LOGGER.error("[%s] Failed: %s", group_name, exc)
            if not args.continue_on_error:
                break

    LOGGER.info("Finished: %d successful/planned; %d failed", succeeded, failed)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

import argparse
import colorsys
import csv
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

TRIMAL = "/home/tomohisa/PycharmProjects/HMSS2/bin/trimal"
VERYFASTTREE = "/home/tomohisa/PycharmProjects/HMSS2/bin/VeryFastTree-avx2"


def load_genome_taxonomy(taxonomy_file: str | Path) -> dict[str, str]:
    """Load non-empty genome taxonomic lineages from a TSV."""
    taxonomy = {}
    with Path(taxonomy_file).open(newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 2 or not row[0].strip():
                continue
            genome, lineage = row[0].strip(), row[1].strip()
            if genome.lower() in {"genome_id", "genomeid", "genome"}:
                continue
            if not lineage or lineage.lower() in {
                "unknown",
                "unclassified",
                "na",
                "n/a",
            }:
                continue
            if genome in taxonomy:
                raise ValueError(f"Duplicate genome ID in taxonomy: {genome}")
            taxonomy[genome] = lineage
    return taxonomy


def combine_seeds(
    seed_files: list[str | Path],
    output_dir: str | Path,
    taxonomy: dict[str, str],
    protein_types: dict[Path, str] | None = None,
    allow_missing_taxonomy: bool = False,
) -> tuple[Path, Path]:
    """Keep only sequences with taxonomy; first occurrence wins on duplicates.

    protein_types maps resolved FASTA paths to explicit annotations; when absent,
    derive the protein type from the filename (legacy behavior).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    combined = output_dir / "combined.faa"
    annotations = output_dir / "annotations.tsv"
    seen_sequences = set()
    seen_ids = set()
    retained = 0
    missing_taxonomy = 0
    duplicates = 0
    missing_genome_ids = set()

    with combined.open("w") as fasta_out, annotations.open("w", newline="") as tsv_out:
        writer = csv.writer(tsv_out, delimiter="\t")
        writer.writerow(["sequence_id", "protein_type"])

        for seed_file in seed_files:
            seed_file = Path(seed_file)
            protein_type = (
                protein_types[seed_file.resolve()]
                if protein_types is not None and seed_file.resolve() in protein_types
                else seed_file.stem.split("_", 1)[-1]
            )

            with seed_file.open() as f:
                sequence_id = None
                sequence = []

                def process_record():
                    nonlocal retained, missing_taxonomy, duplicates
                    if sequence_id is None:
                        return
                    seq = "".join(sequence).upper()
                    if not seq:
                        return
                    genome_id = sequence_id.split("-", 1)[0].split("___", 1)[0]

                    if genome_id not in taxonomy:
                        missing_taxonomy += 1

                        # Collect up to 5 genome IDs without taxonomy
                        if len(missing_genome_ids) < 5:
                            missing_genome_ids.add(genome_id)

                        if not allow_missing_taxonomy:
                            return
                    if seq in seen_sequences or sequence_id in seen_ids:
                        duplicates += 1
                        return
                    seen_sequences.add(seq)
                    seen_ids.add(sequence_id)
                    fasta_out.write(f">{sequence_id}\n{seq}\n")
                    writer.writerow([sequence_id, protein_type])
                    retained += 1

                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        process_record()
                        sequence_id = line[1:].split()[0]
                        sequence = []
                    elif line:
                        sequence.append(line)
                process_record()

    print(
        f"Seeds: retained={retained:,}, no_taxonomy={missing_taxonomy:,}, duplicates={duplicates:,}"
    )
    if retained < (1 if allow_missing_taxonomy else 3):
        raise ValueError(
            "Fewer than 3 taxonomically annotated unique sequences retained"
        )
    if missing_genome_ids:
        print(
            "First 5 genome IDs without taxonomy:",
            ", ".join(sorted(missing_genome_ids)),
        )
    return combined, annotations


def filter_outgroup_against_ingroup(
    ingroup_fasta: str | Path,
    outgroup_fasta: str | Path,
    outgroup_annotations: str | Path,
    filtered_fasta: str | Path,
    filtered_annotations: str | Path,
) -> tuple[int, int]:
    """Exclude outgroup records sharing a FASTA ID or AA sequence with ingroup.

    Also exclude internal outgroup duplicates as a safeguard. Keep the outgroup
    annotation TSV perfectly synchronized with the filtered FASTA.
    """

    def records(path):
        name = None
        parts = []
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(parts).upper()
                    name = line[1:].split()[0]
                    parts = []
                elif name is None:
                    raise ValueError(f"FASTA sequence without header: {path}")
                else:
                    parts.append(line)
            if name is not None:
                yield name, "".join(parts).upper()

    seen_ids = set()
    seen_sequences = set()
    for name, sequence in records(ingroup_fasta):
        seen_ids.add(name)
        seen_sequences.add(sequence)

    with Path(outgroup_annotations).open(newline="") as handle:
        labels = {
            row["sequence_id"]: row["protein_type"]
            for row in csv.DictReader(handle, delimiter="\t")
        }

    kept = 0
    removed = 0
    with (
        Path(filtered_fasta).open("w", encoding="utf-8") as fasta_out,
        Path(filtered_annotations).open("w", newline="", encoding="utf-8") as ann_out,
    ):
        writer = csv.writer(ann_out, delimiter="\t")
        writer.writerow(["sequence_id", "protein_type"])
        for name, sequence in records(outgroup_fasta):
            if name in seen_ids or sequence in seen_sequences:
                removed += 1
                continue
            if name not in labels:
                raise ValueError(f"Missing outgroup annotation for {name}")
            fasta_out.write(f">{name}\n{sequence}\n")
            writer.writerow([name, labels[name]])
            seen_ids.add(name)
            seen_sequences.add(sequence)
            kept += 1

    if kept == 0:
        raise ValueError(
            "No outgroup sequences remain after excluding ingroup duplicates"
        )
    print(f"Outgroup vs ingroup: retained={kept:,}, duplicates_removed={removed:,}")
    return kept, removed


def align_and_trim(
    combined_fasta: str | Path,
    output_dir: str | Path,
    threads: int = 8,
    max_gap_fraction: float = 0.95,
    outgroup_fasta: str | Path | None = None,
    workflow: str = "B",
) -> tuple[Path, Path]:
    """Build the reference MSA according to a selectable outgroup workflow.

    A: Align ingroup -> MAFFT --add outgroup -> trim combined alignment.
    B: Align ingroup -> MAFFT --add --keeplength -> trim combined alignment.
    D: Align ingroup -> trim ingroup -> MAFFT --add --keeplength.

    With no outgroup, simply align and trim the ingroup in all workflows.
    Returns (alignment.fasta_aln, trimmed.fasta_aln). The latter is always
    the FINAL reference MSA used by hmmbuild and VeryFastTree.
    """
    combined_fasta = Path(combined_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    if workflow not in {"A", "B", "D"}:
        raise ValueError(f"Unknown alignment workflow: {workflow}")
    if not combined_fasta.is_file():
        raise FileNotFoundError(combined_fasta)
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if not 0 <= max_gap_fraction <= 1:
        raise ValueError("max_gap_fraction must be between 0 and 1")
    if outgroup_fasta is not None:
        outgroup_fasta = Path(outgroup_fasta).resolve()
        if not outgroup_fasta.is_file():
            raise FileNotFoundError(outgroup_fasta)

    output_dir.mkdir(parents=True, exist_ok=True)
    alignment = output_dir / "alignment.fasta_aln"
    trimmed = output_dir / "trimmed.fasta_aln"
    ingroup_alignment = (
        output_dir / "ingroup_alignment.fasta_aln"
        if outgroup_fasta is not None
        else alignment
    )

    print(
        f"Alignment workflow: {workflow}"
        + (" (no outgroup)" if outgroup_fasta is None else "")
    )

    # Stage 1: MAFFT never sees the outgroup during the ingroup MSA.
    with ingroup_alignment.open("w") as handle:
        subprocess.run(
            ["mafft", "--auto", "--thread", str(threads), str(combined_fasta)],
            stdout=handle,
            check=True,
        )

    def run_trimal(source: Path, destination: Path) -> None:
        subprocess.run(
            [
                TRIMAL,
                "-in",
                str(source),
                "-out",
                str(destination),
                "-fasta",
                "-gt",
                str(1 - max_gap_fraction),
            ],
            check=True,
        )

    if outgroup_fasta is None:
        run_trimal(ingroup_alignment, trimmed)
        return alignment, trimmed

    if workflow == "D":
        # Column selection is based entirely on the ingroup.
        ingroup_trimmed = output_dir / "ingroup_trimmed.fasta_aln"
        run_trimal(ingroup_alignment, ingroup_trimmed)
        mafft_reference = ingroup_trimmed
    else:
        mafft_reference = ingroup_alignment

    # Stage 2: add the divergent outgroup without realigning the ingroup.
    command = ["mafft", "--thread", str(threads), "--add", str(outgroup_fasta)]
    if workflow in {"B", "D"}:
        command.append("--keeplength")
    command.append(str(mafft_reference))
    with alignment.open("w") as handle:
        subprocess.run(command, stdout=handle, check=True)

    if workflow == "D":
        # The final MAFFT output is already column-filtered; do NOT re-trim it.
        # Preserve the standard trimmed.fasta_aln downstream interface.
        shutil.copyfile(alignment, trimmed)
    else:
        run_trimal(alignment, trimmed)
    return alignment, trimmed


def build_hmm(
    trimmed_alignment: str | Path,
    output_dir: str | Path,
    threads: int = 8,
) -> Path:
    """
    Build a protein profile HMM from a trimmed FASTA alignment.

    Returns
    -------
    Path
        Path to the generated HMM file.
    """

    trimmed_alignment = Path(trimmed_alignment).resolve()
    output_dir = Path(output_dir).resolve()

    if not trimmed_alignment.is_file():
        raise FileNotFoundError(trimmed_alignment)

    output_dir.mkdir(parents=True, exist_ok=True)

    hmm_file = output_dir / "alignment.hmm"

    subprocess.run(
        [
            "hmmbuild",
            "--amino",
            "--informat",
            "afa",
            "--symfrac",
            "0.0",
            "--cpu",
            str(threads),
            str(hmm_file),
            str(trimmed_alignment),
        ],
        check=True,
    )

    return hmm_file


def build_tree(
    trimmed_alignment: str | Path,
    output_dir: str | Path,
    threads: int = 8,
) -> Path:
    """
    Build a phylogenetic protein tree using VeryFastTree.

    Model: LG + Gamma approximation

    Returns
    -------
    Path
        Path to the generated Newick tree.
    """

    trimmed_alignment = Path(trimmed_alignment).resolve()
    output_dir = Path(output_dir).resolve()

    if not trimmed_alignment.is_file():
        raise FileNotFoundError(trimmed_alignment)

    if threads < 1:
        raise ValueError("threads must be >= 1")

    executable = VERYFASTTREE

    if executable is None:
        raise RuntimeError("VeryFastTree not found in PATH")

    output_dir.mkdir(parents=True, exist_ok=True)

    tree_file = output_dir / "tree.nwk"

    command = [
        executable,
        "-lg",
        "-gamma",
        "-double-precision",
        "-pseudo",
        "1.0",
        "-threads",
        str(threads),
        str(trimmed_alignment),
    ]

    with tree_file.open("w") as output:
        subprocess.run(
            command,
            stdout=output,
            check=True,
        )

    return tree_file


def create_itol_protein_types(
    annotations_tsv: str | Path,
    output_file: str | Path,
) -> Path:
    """
    Create an iTOL color strip dataset for protein types.

    Input TSV columns:
        sequence_id, protein_type
    """

    annotations_tsv = Path(annotations_tsv)
    output_file = Path(output_file)

    with annotations_tsv.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    protein_types = sorted({row["protein_type"] for row in rows})

    # Generate distinct colors
    colors = {}

    for i, protein in enumerate(protein_types):
        hue = i / len(protein_types)
        rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.9)

        colors[protein] = "#{:02X}{:02X}{:02X}".format(*(round(c * 255) for c in rgb))

    with output_file.open("w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATASET_LABEL\tProtein types\n")
        out.write("COLOR\t#000000\n")
        out.write("STRIP_WIDTH\t25\n")
        out.write("MARGIN\t5\n")
        out.write("SHOW_INTERNAL\t0\n")

        out.write("LEGEND_TITLE\tProtein types\n")
        out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(protein_types)) + "\n")
        out.write(
            "LEGEND_COLORS\t" + "\t".join(colors[p] for p in protein_types) + "\n"
        )
        out.write("LEGEND_LABELS\t" + "\t".join(protein_types) + "\n")

        out.write("DATA\n")

        for row in rows:
            protein = row["protein_type"]

            out.write(f"{row['sequence_id']}\t{colors[protein]}\t{protein}\n")

    return output_file


def create_gappa_taxonomy(
    annotations_tsv: str | Path,
    taxonomy: dict[str, str],
    output_file: str | Path,
) -> Path:
    """Write headerless leaf-ID <TAB> taxonomy for gappa examine assign."""
    matched = 0
    with (
        Path(annotations_tsv).open(newline="") as source,
        Path(output_file).open("w") as out,
    ):
        for row in csv.DictReader(source, delimiter="\t"):
            sequence_id = row["sequence_id"]
            genome_id = sequence_id.split("-", 1)[0].split("___", 1)[0]
            out.write(f"{sequence_id}\t{taxonomy.get(genome_id, 'unknown')}\n")
            matched += 1
    print(f"GAPPA taxonomy: matched={matched:,}")
    return Path(output_file)


def make_tree_bifurcating(tree_file: str | Path) -> Path:
    """Resolve Newick polytomies using plain Python (no ETE3/Biopython).

    All original leaf labels and branch lengths are retained. Artificial
    branches are assigned length 0.0; this does not establish a biological root.
    """
    tree_file = Path(tree_file)
    newick = tree_file.read_text().strip()
    if not newick.endswith(";"):
        raise ValueError("Invalid Newick: missing terminating semicolon")

    # Parse using a stack, avoiding recursion limits for large trees.
    # Nodes are either leaf strings or [children, suffix_after_closing_parenthesis].
    stack = []
    root = None
    i = 0
    polytomies = 0
    while i < len(newick):
        char = newick[i]
        if char == "(":
            stack.append([])
            i += 1
        elif char == ")":
            if not stack or len(stack[-1]) < 2:
                raise ValueError("Invalid Newick: malformed internal node")
            children = stack.pop()
            if len(children) > 2:
                polytomies += 1
            i += 1
            start = i
            while i < len(newick) and newick[i] not in ",);":
                i += 1
            node = (children, newick[start:i])
            if stack:
                stack[-1].append(node)
            elif root is None:
                root = node
            else:
                raise ValueError("Invalid Newick: multiple roots")
        elif char == ",":
            i += 1
        elif char == ";":
            i += 1
            if i != len(newick):
                raise ValueError("Unexpected content after Newick terminator")
        else:
            start = i
            while i < len(newick) and newick[i] not in "(),;":
                i += 1
            leaf = newick[start:i]
            if not leaf.strip() or not stack:
                raise ValueError("Invalid Newick: unexpected leaf")
            stack[-1].append(leaf)

    if stack or root is None:
        raise ValueError("Invalid Newick: unbalanced parentheses")

    # Iterative postorder serialisation; original existing edges are unmodified.
    converted = {}
    pending = [(root, False)]
    while pending:
        node, ready = pending.pop()
        if isinstance(node, str):
            continue
        children, suffix = node
        if not ready:
            pending.append((node, True))
            pending.extend(
                (child, False) for child in children if not isinstance(child, str)
            )
        else:
            text_children = [
                child if isinstance(child, str) else converted[id(child)]
                for child in children
            ]
            while len(text_children) > 2:
                first, second, *remaining = text_children
                text_children = [f"({first},{second}):0.0", *remaining]
            converted[id(node)] = f"({','.join(text_children)}){suffix}"

    output_tree = tree_file.with_name("tree_bifurcating.nwk")
    output_tree.write_text(converted[id(root)] + ";\n")
    print(f"Resolved {polytomies} polytomous nodes: {output_tree}")
    return output_tree


## Von hier an ist es die ausführung der routinen


def sha256_file(path: str | Path) -> str:
    """Hash a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tool_version(command: list[str]) -> str | None:
    """Best-effort tool-version detection; never blocks package creation."""
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        lines = (result.stdout + "\n" + result.stderr).splitlines()
        for line in lines:
            if re.search(r"\b(version|v\d+[.]\d+|MAFFT v|HMMER)\b", line, re.I):
                return line.strip()[:250]
        return next((line.strip()[:250] for line in lines if line.strip()), None)
    except (OSError, subprocess.TimeoutExpired):
        return None


def create_manifest(
    output_dir: str | Path,
    package_version: str,
    seed_files: list[Path],
    taxonomy_file: Path,
    threads: int,
    outgroup_files: list[dict] | None = None,
    workflow: str = "B",
) -> Path:
    """Write manifest last; checksums cover every generated package file."""
    output_dir = Path(output_dir).resolve()
    file_roles = {
        "sequences": "combined.faa",
        "protein_annotations": "annotations.tsv",
        "taxonomy": "reference_taxonomy.tsv",
        "itol_annotations": "itol_protein_types.txt",
        "alignment": "alignment.fasta_aln",
        "trimmed_alignment": "trimmed.fasta_aln",
        "hmm": "alignment.hmm",
        "raw_tree": "tree.nwk",
        "tree": "tree_bifurcating.nwk",
    }
    # Include any additional generated files, while excluding the manifest itself.
    files = {}
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "manifest.json":
            files[path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
    for role, name in file_roles.items():
        if name not in files:
            raise FileNotFoundError(f"Missing required package file ({role}): {name}")

    manifest = {
        "schema_version": "1.0",
        "package_id": output_dir.name,
        "package_version": package_version,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files": files,
        "roles": file_roles,
        "parameters": {
            "threads": threads,
            "max_gap_fraction": 0.95,
            "outgroup_alignment_method": "mafft --add" if outgroup_files else None,
            "alignment_workflow": workflow,
            "outgroup_keeplength": bool(outgroup_files) and workflow in {"B", "D"},
            "trimming_scope": "ingroup"
            if not outgroup_files or workflow == "D"
            else "combined",
        },
        "outgroups": outgroup_files or [],
        "tools": {
            "mafft": tool_version(["mafft", "--version"]),
            "trimal": tool_version([TRIMAL, "-version"]),
            "hmmbuild": tool_version(["hmmbuild", "-h"]),
            "veryfasttree": tool_version([VERYFASTTREE, "-help"]),
        },
        "sources": {
            "seed_fastas": [
                {"path": str(p.resolve()), "sha256": sha256_file(p)} for p in seed_files
            ],
            "genome_taxonomy": {
                "path": str(taxonomy_file.resolve()),
                "sha256": sha256_file(taxonomy_file),
            },
        },
    }
    target = output_dir / "manifest.json"
    temporary = output_dir / "manifest.json.tmp"
    try:
        temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def main():
    parser = argparse.ArgumentParser(
        description="Build a phylogenetic reference package from protein seed FASTAs."
    )

    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=False,
        type=Path,
        help="Seed FASTAs with protein types inferred from filenames (legacy mode).",
    )

    parser.add_argument(
        "--seed",
        nargs=2,
        action="append",
        metavar=("FASTA", "PROTEIN_TYPE"),
        default=[],
        help=(
            "Explicit FASTA and protein-type pairing. Repeat as needed, e.g. "
            "--seed /data/MarD.faa MarD. "
            "Can be combined with -i; explicitly named seeds override the filename rule."
        ),
    )

    parser.add_argument(
        "--outgroup",
        nargs=2,
        action="append",
        metavar=("FASTA", "LABEL"),
        default=[],
        help="Outgroup FASTA and its annotation label; repeat as needed. "
        "Added after the ingroup MSA with MAFFT --add.",
    )
    parser.add_argument(
        "--alignment-workflow",
        choices=("A", "B", "D"),
        default="B",
        help=(
            "A: add outgroup then trim all (no --keeplength); "
            "B: add --keeplength then trim all (default); "
            "D: trim ingroup first then add --keeplength (no second trimming)."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Number of threads (default: 8).",
    )

    parser.add_argument(
        "--taxonomy",
        required=True,
        type=Path,
        help="Genome taxonomy TSV (genome_id, semicolon-separated lineage).",
    )

    parser.add_argument(
        "--package-version",
        default="1.0.0",
        help="Reference package version (default: 1.0.0).",
    )

    args = parser.parse_args()

    if not args.taxonomy.is_file():
        parser.error(f"Taxonomy file not found: {args.taxonomy}")

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    if not args.input and not args.seed:
        parser.error(
            "Provide at least one seed via -i/--input or --seed FASTA PROTEIN_TYPE"
        )

    # Keep legacy inputs first, followed by explicitly annotated seeds.
    # A file is accepted only once, and conflicting type assignments are rejected.
    seed_files = []
    protein_types = {}
    seen_paths = set()

    for seed_file in args.input or []:
        path = seed_file.resolve()
        if not path.is_file():
            parser.error(f"Seed file not found: {path}")
        if path not in seen_paths:
            seed_files.append(path)
            seen_paths.add(path)

    for raw_path, protein_type in args.seed:
        path = Path(raw_path).resolve()
        protein_type = protein_type.strip()
        if not path.is_file():
            parser.error(f"Seed file not found: {path}")
        if not protein_type:
            parser.error(f"Empty protein type for {path}")
        existing = protein_types.get(path)
        if existing is not None and existing != protein_type:
            parser.error(
                f"Conflicting protein types for {path}: {existing} vs {protein_type}"
            )
        protein_types[path] = protein_type
        if path not in seen_paths:
            seed_files.append(path)
            seen_paths.add(path)

    outgroup_files = []
    outgroup_types = {}
    for raw_path, label in args.outgroup:
        path = Path(raw_path).resolve()
        label = label.strip()
        if not path.is_file():
            parser.error(f"Outgroup FASTA not found: {path}")
        if not label:
            parser.error(f"Empty outgroup label for {path}")
        if path in seen_paths:
            parser.error(f"File appears as both seed and outgroup: {path}")
        if path in outgroup_types and outgroup_types[path] != label:
            parser.error(f"Conflicting outgroup labels for {path}")
        if path not in outgroup_types:
            outgroup_files.append(path)
            outgroup_types[path] = label

    for executable in (TRIMAL, VERYFASTTREE):
        if not Path(executable).is_file():
            parser.error(f"Executable not found: {executable}")

    output_dir = args.output.resolve()

    # 1. Load taxonomy and discard unannotated seeds before any downstream analysis
    taxonomy = load_genome_taxonomy(args.taxonomy)
    combined, annotations = combine_seeds(
        seed_files=seed_files,
        output_dir=output_dir,
        taxonomy=taxonomy,
        protein_types=protein_types,
    )

    outgroup_fasta = None
    if outgroup_files:
        outgroup_dir = output_dir / "_outgroup_inputs"
        outgroup_combined, outgroup_annotations = combine_seeds(
            seed_files=outgroup_files,
            output_dir=outgroup_dir,
            taxonomy=taxonomy,
            protein_types=outgroup_types,
            allow_missing_taxonomy=True,
        )
        # Preserve unmodified ingroup FASTA for the initial MAFFT alignment.
        ingroup_fasta = output_dir / "ingroup.faa"
        ingroup_fasta.write_bytes(combined.read_bytes())

        # Filter the outgroup against the ingroup BEFORE MAFFT --add.
        # Synchronize all downstream FASTA and annotation files with the filter.
        outgroup_fasta = output_dir / "outgroup.faa"
        filtered_outgroup_annotations = output_dir / "outgroup_annotations.tsv"
        filter_outgroup_against_ingroup(
            ingroup_fasta=ingroup_fasta,
            outgroup_fasta=outgroup_combined,
            outgroup_annotations=outgroup_annotations,
            filtered_fasta=outgroup_fasta,
            filtered_annotations=filtered_outgroup_annotations,
        )

        # Append only retained outgroup annotations.
        with (
            annotations.open("a") as target,
            filtered_outgroup_annotations.open() as source,
        ):
            next(source)
            for line in source:
                target.write(line)

        # combined.faa must match the final reference alignment exactly.
        with combined.open("ab") as target, outgroup_fasta.open("rb") as source:
            target.write(b"\n")
            for block in iter(lambda: source.read(1024 * 1024), b""):
                target.write(block)
    else:
        ingroup_fasta = combined

    # 2. Create iTOL annotation
    create_itol_protein_types(
        annotations_tsv=annotations,
        output_file=output_dir / "itol_protein_types.txt",
    )

    # 3. GAPPA taxonomy for each reference leaf
    create_gappa_taxonomy(
        annotations_tsv=annotations,
        taxonomy=taxonomy,
        output_file=output_dir / "reference_taxonomy.tsv",
    )

    # 4. MAFFT alignment and trimAl
    alignment, trimmed = align_and_trim(
        combined_fasta=ingroup_fasta,
        output_dir=output_dir,
        threads=args.threads,
        max_gap_fraction=0.95,
        outgroup_fasta=outgroup_fasta,
        workflow=args.alignment_workflow,
    )

    # 5. Build HMM
    build_hmm(
        trimmed_alignment=trimmed,
        output_dir=output_dir,
        threads=args.threads,
    )

    # 6. Build phylogenetic tree
    tree_file = build_tree(
        trimmed_alignment=trimmed,
        output_dir=output_dir,
        threads=args.threads,
    )

    # 7. Produce fully bifurcating reference tree for EPA-ng
    make_tree_bifurcating(tree_file)

    # 8. Record checksums, input provenance and tool versions
    manifest_file = create_manifest(
        output_dir=output_dir,
        package_version=args.package_version,
        seed_files=seed_files + outgroup_files,
        taxonomy_file=args.taxonomy,
        threads=args.threads,
        outgroup_files=[
            {"path": str(p), "protein_type": outgroup_types[p]} for p in outgroup_files
        ],
        workflow=args.alignment_workflow,
    )
    print(f"Manifest created: {manifest_file}")

    print(f"Reference package created: {output_dir}")


if __name__ == "__main__":
    main()

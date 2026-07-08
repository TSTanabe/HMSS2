#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

INPUT_COL = "input_sulfur_species_or_compound"
OUTPUT_COL = "output_sulfur_species_or_compound"
ENZYME_COL = "enzyme_abbreviation"
PATHWAY_COL = "pathway"


def split_compounds(value):
    if value is None:
        return []
    value = str(value).strip()
    if not value:
        return []
    parts = re.split(r"\s*;\s*|\s*/\s*", value)
    return [p.strip() for p in parts if p.strip()]


def normalize_compound(value):
    value = str(value).strip().lower()
    value = value.replace("–", "-").replace("—", "-")
    value = re.sub(r"\s+", " ", value)
    return value


def split_enzymes(value):
    if value is None:
        return []
    return [x.strip() for x in str(value).strip().split(";") if x.strip()]


def join_unique(values):
    seen = set()
    out = []
    for v in values:
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return ";".join(out)


def load_reactions(path):
    with path.open("r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(f, dialect=dialect)

        required = [INPUT_COL, OUTPUT_COL, ENZYME_COL, PATHWAY_COL]
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {missing}\nFound: {reader.fieldnames}")

        reactions = []
        for idx, row in enumerate(reader, start=1):
            in_raw = row.get(INPUT_COL, "").strip()
            out_raw = row.get(OUTPUT_COL, "").strip()
            enzyme_raw = row.get(ENZYME_COL, "").strip()
            pathway_raw = row.get(PATHWAY_COL, "").strip()

            in_parts = split_compounds(in_raw)
            out_parts = split_compounds(out_raw)

            reactions.append({
                "id": idx,
                "input_raw": in_raw,
                "output_raw": out_raw,
                "enzyme_raw": enzyme_raw,
                "pathway_raw": pathway_raw,
                "input_norm": {normalize_compound(x) for x in in_parts},
                "output_norm": {normalize_compound(x) for x in out_parts},
                "enzymes": split_enzymes(enzyme_raw),
            })

    return reactions


def is_direct_reverse(a, b):
    return (
        bool(a["input_norm"])
        and bool(a["output_norm"])
        and a["input_norm"] == b["output_norm"]
        and a["output_norm"] == b["input_norm"]
    )


def creates_cycle(upstream_inputs, next_outputs):
    return bool(set(upstream_inputs) & set(next_outputs))


def combine_two_step(reactions):
    by_input = defaultdict(list)
    for r in reactions:
        for compound in r["input_norm"]:
            by_input[compound].append(r)

    combos = []
    seen = set()

    for a in reactions:
        for intermediate in a["output_norm"]:
            for b in by_input.get(intermediate, []):
                if a["id"] == b["id"]:
                    continue
                if is_direct_reverse(a, b):
                    continue
                if creates_cycle(a["input_norm"], b["output_norm"]):
                    continue

                key = (a["id"], b["id"], intermediate)
                if key in seen:
                    continue
                seen.add(key)

                pathway = a["pathway_raw"]
                if b["pathway_raw"] and b["pathway_raw"] != pathway:
                    pathway = f"{pathway} + {b['pathway_raw']}" if pathway else b["pathway_raw"]

                combos.append({
                    "input": a["input_raw"],
                    "output": b["output_raw"],
                    "enzymes": join_unique(a["enzymes"] + b["enzymes"]),
                    "pathway": pathway,
                    "reaction_ids": f"{a['id']};{b['id']}",
                    "intermediate": intermediate,
                })

    return combos


def combine_paths(reactions, max_steps):
    by_input = defaultdict(list)
    by_id = {r["id"]: r for r in reactions}

    for r in reactions:
        for c in r["input_norm"]:
            by_input[c].append(r)

    active = []
    for r in reactions:
        if r["input_norm"] and r["output_norm"]:
            active.append({
                "ids": [r["id"]],
                "input_seen": set(r["input_norm"]),
                "outputs": set(r["output_norm"]),
            })

    all_paths = []
    seen_final = set()

    for _length in range(2, max_steps + 1):
        new_active = []
        seen_paths_this_round = set()

        for path in active:
            last = by_id[path["ids"][-1]]

            for intermediate in path["outputs"]:
                for nxt in by_input.get(intermediate, []):
                    if nxt["id"] in path["ids"]:
                        continue
                    if is_direct_reverse(last, nxt):
                        continue
                    if creates_cycle(path["input_seen"], nxt["output_norm"]):
                        continue

                    new_ids = path["ids"] + [nxt["id"]]
                    path_key = tuple(new_ids)
                    if path_key in seen_paths_this_round:
                        continue
                    seen_paths_this_round.add(path_key)

                    new_active.append({
                        "ids": new_ids,
                        "input_seen": set(path["input_seen"]) | set(nxt["input_norm"]),
                        "outputs": set(nxt["output_norm"]),
                    })

                    first = by_id[new_ids[0]]
                    last_r = by_id[new_ids[-1]]
                    enzymes = []
                    pathway_parts = []

                    for rid in new_ids:
                        r = by_id[rid]
                        enzymes.extend(r["enzymes"])
                        if r["pathway_raw"] and r["pathway_raw"] not in pathway_parts:
                            pathway_parts.append(r["pathway_raw"])

                    final_key = (tuple(new_ids), first["input_raw"], last_r["output_raw"])
                    if final_key in seen_final:
                        continue
                    seen_final.add(final_key)

                    all_paths.append({
                        "input": first["input_raw"],
                        "output": last_r["output_raw"],
                        "enzymes": join_unique(enzymes),
                        "pathway": " + ".join(pathway_parts),
                        "reaction_ids": ";".join(map(str, new_ids)),
                    })

        active = new_active

    return all_paths


def write_reaction_format_tsv(rows, out_path):
    fieldnames = [
        "pathway_id",
        INPUT_COL,
        OUTPUT_COL,
        ENZYME_COL,
        PATHWAY_COL,
    ]

    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for idx, row in enumerate(rows, start=1):
            writer.writerow({
                "pathway_id": f"{idx:06d}",
                INPUT_COL: row["input"],
                OUTPUT_COL: row["output"],
                ENZYME_COL: row["enzymes"],
                PATHWAY_COL: "",
            })


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine sulfur reaction rows into pathway candidates by matching "
            "output compounds to input compounds. Output is written in a simple "
            "reaction-table format with pathway_id, input, output, enzymes, pathway."
        )
    )
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--print-all", action="store_true")
    parser.add_argument("--max-print", type=int, default=50)
    args = parser.parse_args()

    reactions = load_reactions(args.input)
    usable = [r for r in reactions if r["input_norm"] and r["output_norm"] and r["enzymes"]]

    if args.max_steps == 2:
        combinations = combine_two_step(usable)
    else:
        combinations = combine_paths(usable, args.max_steps)

    write_reaction_format_tsv(combinations, args.output)

    print("=== Reaction pathway precomputation ===")
    print(f"Input reactions total       : {len(reactions)}")
    print(f"Usable reactions            : {len(usable)}")
    print(f"Skipped incomplete reactions: {len(reactions) - len(usable)}")
    print(f"Maximum pathway length      : {args.max_steps}")
    print(f"Combinations found          : {len(combinations)}")
    print(f"[OK] wrote TSV             : {args.output}")
    print()

    n = len(combinations) if args.print_all else min(len(combinations), args.max_print)
    print(f"=== Combinations shown: {n}/{len(combinations)} ===")

    for i, c in enumerate(combinations[:n], start=1):
        print()
        print(f"[pathway_{i:06d}]")
        print(f"    input   : {c['input']}")
        print(f"    output  : {c['output']}")
        print(f"    enzymes : {c['enzymes']}")
        print(f"    pathway : {c['pathway']}")
        print(f"    source reaction IDs: {c.get('reaction_ids', '')}")
        if c.get("intermediate"):
            print(f"    matched via        : {c['intermediate']}")

    if len(combinations) > n:
        print()
        print(f"... {len(combinations) - n} further combinations not shown. Use --print-all to print all.")


if __name__ == "__main__":
    main()

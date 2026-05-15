#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import random
import sqlite3
import statistics
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, List, Tuple


def resolve_next_to_script(filename: str) -> str:
    return str(Path(__file__).resolve().parent / filename)


def chunked(values: List[Any], size: int):
    for i in range(0, len(values), size):
        yield values[i:i + size]


def readonly_connection(database: str) -> sqlite3.Connection:
    db_path = os.path.abspath(database)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, timeout=60)
    con.execute("PRAGMA query_only = TRUE;")
    con.execute("PRAGMA cache_size = 100000;")
    con.execute("PRAGMA temp_store = MEMORY;")
    con.execute("PRAGMA mmap_size = 2147483648;")
    return con


def read_domain_list(path: str) -> List[str]:
    domains = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                domains.append(line.split()[0])
    return domains


def read_excluded_features(path: str | None) -> set[str]:
    if path is None:
        return set()

    excluded = set()

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            excluded.add(line.split()[0])

    return excluded


def fetch_positive_genomes_for_domain(database: str, target_domain: str) -> List[str]:
    query = """
        SELECT DISTINCT p.genomeID
        FROM Proteins p
        JOIN Domains d ON p.proteinID = d.proteinID
        WHERE p.valid_hit = 1
          AND d.domain = ?
          AND p.genomeID IS NOT NULL
          AND p.genomeID != ''
        ORDER BY p.genomeID
    """

    with readonly_connection(database) as con:
        cur = con.cursor()
        cur.execute(query, (target_domain,))
        return [row[0] for row in cur.fetchall()]


def fetch_decoy_genomes_for_domain(
        database: str,
        target_domain: str,
        n_decoys: int,
        seed: int,
        exclude_genomes: set[str] | None = None,
) -> List[str]:
    exclude_genomes = exclude_genomes or set()

    query = """
        SELECT DISTINCT p.genomeID
        FROM Proteins p
        WHERE p.genomeID IS NOT NULL
          AND p.genomeID != ''
          AND p.genomeID NOT IN (
              SELECT DISTINCT p2.genomeID
              FROM Proteins p2
              JOIN Domains d2 ON p2.proteinID = d2.proteinID
              WHERE p2.valid_hit = 1
                AND d2.domain = ?
          )
        ORDER BY p.genomeID
    """

    with readonly_connection(database) as con:
        cur = con.cursor()
        cur.execute(query, (target_domain,))
        candidates = [
            row[0]
            for row in cur.fetchall()
            if row[0] not in exclude_genomes
        ]

    rng = random.Random(seed)

    if len(candidates) > n_decoys:
        return sorted(rng.sample(candidates, n_decoys))

    return candidates


def fetch_domains_for_genome_chunk(
        args: Tuple[str, List[str]],
) -> Dict[str, Dict[str, List[str]]]:
    database, genome_chunk = args

    if not genome_chunk:
        return {}

    placeholders = ",".join("?" for _ in genome_chunk)

    query = f"""
        SELECT
            p.genomeID,
            d.domain,
            p.proteinID
        FROM Proteins p
        JOIN Domains d ON p.proteinID = d.proteinID
        WHERE p.valid_hit = 1
          AND p.genomeID IN ({placeholders})
          AND d.domain IS NOT NULL
          AND d.domain != ''
    """

    out = defaultdict(lambda: defaultdict(list))

    with readonly_connection(database) as con:
        cur = con.cursor()
        cur.execute(query, genome_chunk)

        for genome_id, domain, protein_id in cur:
            out[genome_id][domain].append(protein_id)

    return {gid: dict(dommap) for gid, dommap in out.items()}


def fetch_pam_for_genomes(
        database: str,
        genome_ids: List[str],
        chunk_size: int,
        cores: int,
) -> Dict[str, Dict[str, List[str]]]:
    if not genome_ids:
        return {}

    genome_chunks = list(chunked(genome_ids, chunk_size))
    worker_args = [(database, chunk) for chunk in genome_chunks]

    pam = defaultdict(lambda: defaultdict(list))

    with Pool(processes=cores) as pool:
        for i, partial in enumerate(
                pool.imap_unordered(fetch_domains_for_genome_chunk, worker_args),
                1,
        ):
            for genome_id, domain_map in partial.items():
                for domain, protein_ids in domain_map.items():
                    pam[genome_id][domain].extend(protein_ids)

            if i == 1 or i % 10 == 0 or i == len(genome_chunks):
                print(
                    f"    [FETCH] chunk {i}/{len(genome_chunks)} "
                    f"({(i / len(genome_chunks)) * 100:.1f}%)"
                )

    return {gid: dict(dommap) for gid, dommap in pam.items()}


def split_positive_genomes(
        genome_ids: List[str],
        seed: int,
) -> Tuple[List[str], List[str], List[str]]:
    ids = list(genome_ids)
    rng = random.Random(seed)
    rng.shuffle(ids)

    n = len(ids)
    n_train = int(n * 0.70)
    n_calibration = int(n * 0.10)

    train = ids[:n_train]
    calibration = ids[n_train:n_train + n_calibration]
    test = ids[n_train + n_calibration:]

    return train, calibration, test


def train_positive_support_model_from_pam(
        pam: Dict[str, Dict[str, List[str]]],
        target_domain: str,
        alpha: float,
        min_feature_genomes: int,
        max_model_features: int | None,
        excluded_features: set[str],
) -> Dict[str, Any]:
    positive_genomes = [
        genome_id
        for genome_id, dommap in pam.items()
        if target_domain in dommap and dommap[target_domain]
    ]

    if not positive_genomes:
        raise ValueError(f"No positive training genomes for {target_domain}")

    n_train = len(positive_genomes)
    feature_counts: Dict[str, int] = {}

    for genome_id in positive_genomes:
        dommap = pam[genome_id]

        for domain, protein_ids in dommap.items():
            if domain == target_domain:
                continue

            if protein_ids:
                feature_counts[domain] = feature_counts.get(domain, 0) + 1

    initial_features = sorted(
        feature
        for feature, count in feature_counts.items()
        if (
                count >= min_feature_genomes
                and feature not in excluded_features
                and feature != target_domain
        )
    )

    feature_probs = {}
    feature_strengths = {}

    for feature in initial_features:
        count = feature_counts[feature]
        p = (count + alpha) / (n_train + 2.0 * alpha)
        p = min(max(p, 1e-12), 1.0 - 1e-12)

        feature_probs[feature] = p
        feature_strengths[feature] = -math.log(1.0 - p)

    ranked_features = sorted(
        initial_features,
        key=lambda f: (
            feature_strengths.get(f, 0.0),
            feature_probs.get(f, 0.0),
            feature_counts.get(f, 0),
        ),
        reverse=True,
    )

    if max_model_features is not None and max_model_features > 0:
        features = ranked_features[:max_model_features]
    else:
        features = ranked_features

    feature_probs = {f: feature_probs[f] for f in features}
    feature_strengths = {f: feature_strengths[f] for f in features}

    total_strength = sum(feature_strengths.values())

    if total_strength > 0:
        feature_weights = {
            feature: strength / total_strength
            for feature, strength in feature_strengths.items()
        }
    else:
        feature_weights = {feature: 0.0 for feature in features}

    return {
        "model_type": "positive_only_bernoulli_support",
        "score_definition": "weighted_feature_sum",
        "target_domain": target_domain,
        "alpha": alpha,
        "min_feature_genomes": min_feature_genomes,
        "max_model_features": max_model_features,
        "n_train_genomes": n_train,
        "n_initial_features": len(initial_features),
        "n_features": len(features),
        "features": features,
        "feature_counts": {
            feature: int(feature_counts[feature])
            for feature in features
        },
        "feature_probs": feature_probs,
        "feature_strengths": feature_strengths,
        "feature_weights": feature_weights,
    }


def score_genome(dommap: Dict[str, List[str]], model: Dict[str, Any]) -> float:
    return float(
        sum(
            weight
            for feature, weight in model["feature_weights"].items()
            if dommap.get(feature)
        )
    )


def score_pam(
        pam: Dict[str, Dict[str, List[str]]],
        genome_ids: List[str],
        model: Dict[str, Any],
) -> Dict[str, float]:
    return {
        genome_id: score_genome(pam.get(genome_id, {}), model)
        for genome_id in genome_ids
    }


def confusion_at_threshold(
        positive_scores: List[float],
        decoy_scores: List[float],
        threshold: float,
) -> Dict[str, float]:
    tp = sum(1 for s in positive_scores if s >= threshold)
    fn = sum(1 for s in positive_scores if s < threshold)

    fp = sum(1 for s in decoy_scores if s >= threshold)
    tn = sum(1 for s in decoy_scores if s < threshold)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    return {
        "threshold": float(threshold),
        "TP": int(tp),
        "FP": int(fp),
        "TN": int(tn),
        "FN": int(fn),
        "precision": float(precision),
        "recall": float(recall),
        "F1": float(f1),
    }


def find_best_f1_threshold(
        calibration_scores: List[float],
        decoy_scores: List[float],
) -> Dict[str, float]:
    candidate_thresholds = sorted(
        set(calibration_scores + decoy_scores),
        reverse=True,
    )

    if not candidate_thresholds:
        return {
            "threshold": 1.0,
            "TP": 0,
            "FP": 0,
            "TN": 0,
            "FN": 0,
            "precision": 0.0,
            "recall": 0.0,
            "F1": 0.0,
        }

    best = None

    for threshold in candidate_thresholds:
        cm = confusion_at_threshold(
            positive_scores=calibration_scores,
            decoy_scores=decoy_scores,
            threshold=threshold,
        )

        if best is None:
            best = cm
            continue

        if (
                cm["F1"] > best["F1"]
                or (
                cm["F1"] == best["F1"]
                and cm["FP"] < best["FP"]
        )
                or (
                cm["F1"] == best["F1"]
                and cm["FP"] == best["FP"]
                and cm["threshold"] > best["threshold"]
        )
        ):
            best = cm

    return best


def safe_median(values: List[float]) -> float:
    return float(statistics.median(values)) if values else 0.0


def safe_min(values: List[float]) -> float:
    return float(min(values)) if values else 0.0


def safe_max(values: List[float]) -> float:
    return float(max(values)) if values else 0.0


def append_jsonl(path: str, record: Dict[str, Any]) -> None:
    with open(path, "a") as handle:
        handle.write(json.dumps(record) + "\n")


def write_evaluation_header(path: str) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "target_domain",
            "status",
            "total_positive_genomes",
            "n_train",
            "n_calibration_positive",
            "n_calibration_decoy",
            "n_test_positive",
            "n_test_decoy",
            "n_initial_features",
            "n_features",
            "max_model_features",
            "threshold",
            "calibration_precision",
            "calibration_recall",
            "calibration_F1",
            "calibration_TP",
            "calibration_FP",
            "calibration_TN",
            "calibration_FN",
            "test_precision",
            "test_recall",
            "test_F1",
            "test_TP",
            "test_FP",
            "test_TN",
            "test_FN",
            "runtime_seconds",
            "error",
        ])


def append_evaluation_row(path: str, row: Dict[str, Any]) -> None:
    with open(path, "a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            row.get("target_domain", ""),
            row.get("status", ""),
            row.get("total_positive_genomes", ""),
            row.get("n_train", ""),
            row.get("n_calibration_positive", ""),
            row.get("n_calibration_decoy", ""),
            row.get("n_test_positive", ""),
            row.get("n_test_decoy", ""),
            row.get("n_initial_features", ""),
            row.get("n_features", ""),
            row.get("max_model_features", ""),
            row.get("threshold", ""),
            row.get("calibration_precision", ""),
            row.get("calibration_recall", ""),
            row.get("calibration_F1", ""),
            row.get("calibration_TP", ""),
            row.get("calibration_FP", ""),
            row.get("calibration_TN", ""),
            row.get("calibration_FN", ""),
            row.get("test_precision", ""),
            row.get("test_recall", ""),
            row.get("test_F1", ""),
            row.get("test_TP", ""),
            row.get("test_FP", ""),
            row.get("test_TN", ""),
            row.get("test_FN", ""),
            row.get("runtime_seconds", ""),
            row.get("error", ""),
        ])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Train positive-only Bernoulli support models. "
            "Threshold is optimized on 10% calibration positives vs equal-sized decoys. "
            "Final performance is evaluated on 20% test positives vs equal-sized independent decoys."
        )
    )

    parser.add_argument("--database", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument(
        "--excluded-features",
        default=None,
        help=(
            "Optional file containing domains that should never be "
            "used as model features."
        ),
    )
    parser.add_argument("--models-jsonl", default="positive_support_models.jsonl")
    parser.add_argument("--evaluation-report", default="positive_support_test_evaluation.tsv")

    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=900)
    parser.add_argument("--random-seed", type=int, default=42)

    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--min-feature-genomes", type=int, default=2)
    parser.add_argument("--min-positives", type=int, default=30)

    parser.add_argument("--max-model-features", type=int, default=20)
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    if args.max_model_features <= 0:
        args.max_model_features = None

    args.models_jsonl = resolve_next_to_script(args.models_jsonl)
    args.evaluation_report = resolve_next_to_script(args.evaluation_report)

    for path in [args.models_jsonl, args.evaluation_report]:
        if os.path.exists(path):
            if args.overwrite:
                os.remove(path)
            else:
                raise FileExistsError(f"Output exists: {path}. Use --overwrite.")

    write_evaluation_header(args.evaluation_report)

    domains = read_domain_list(args.domains)
    excluded_features = read_excluded_features(args.excluded_features)

    trained = 0
    skipped = 0
    start_all = time.time()

    for i, domain in enumerate(domains, 1):
        start = time.time()

        print(f"\n[DOMAIN {i}/{len(domains)}] {domain}")

        try:
            positives = fetch_positive_genomes_for_domain(args.database, domain)
            total_pos = len(positives)

            print(f"    [INFO] Total positive genomes: {total_pos:,}")

            if total_pos < args.min_positives:
                skipped += 1

                model = {
                    "target_domain": domain,
                    "status": "skipped_too_few_positives",
                    "total_positive_genomes": total_pos,
                    "min_positives": args.min_positives,
                }

                append_jsonl(args.models_jsonl, model)

                append_evaluation_row(args.evaluation_report, {
                    "target_domain": domain,
                    "status": "skipped_too_few_positives",
                    "total_positive_genomes": total_pos,
                    "runtime_seconds": round(time.time() - start, 3),
                })

                continue

            train_ids, calibration_ids, test_ids = split_positive_genomes(
                positives,
                seed=args.random_seed + i,
            )

            print(
                f"    [SPLIT] train={len(train_ids):,}, "
                f"calibration_positive={len(calibration_ids):,}, "
                f"test_positive={len(test_ids):,}"
            )

            print("    [INFO] Fetching training PAM")
            train_pam = fetch_pam_for_genomes(
                args.database,
                train_ids,
                args.chunk_size,
                args.cores,
            )

            model = train_positive_support_model_from_pam(
                pam=train_pam,
                target_domain=domain,
                alpha=args.alpha,
                min_feature_genomes=args.min_feature_genomes,
                max_model_features=args.max_model_features,
                excluded_features=excluded_features,
            )
            calibration_decoy_ids = fetch_decoy_genomes_for_domain(
                database=args.database,
                target_domain=domain,
                n_decoys=len(calibration_ids),
                seed=args.random_seed + i + 10_000,
            )

            test_decoy_ids = fetch_decoy_genomes_for_domain(
                database=args.database,
                target_domain=domain,
                n_decoys=len(test_ids),
                seed=args.random_seed + i + 20_000,
                exclude_genomes=set(calibration_decoy_ids),
            )

            print(
                f"    [DECOYS] calibration_decoy={len(calibration_decoy_ids):,}, "
                f"test_decoy={len(test_decoy_ids):,}"
            )

            calibration_pam = fetch_pam_for_genomes(
                args.database,
                calibration_ids,
                args.chunk_size,
                args.cores,
            )

            calibration_decoy_pam = fetch_pam_for_genomes(
                args.database,
                calibration_decoy_ids,
                args.chunk_size,
                args.cores,
            )

            test_pam = fetch_pam_for_genomes(
                args.database,
                test_ids,
                args.chunk_size,
                args.cores,
            )

            test_decoy_pam = fetch_pam_for_genomes(
                args.database,
                test_decoy_ids,
                args.chunk_size,
                args.cores,
            )

            calibration_scores = list(
                score_pam(calibration_pam, calibration_ids, model).values()
            )

            calibration_decoy_scores = list(
                score_pam(calibration_decoy_pam, calibration_decoy_ids, model).values()
            )

            test_scores = list(
                score_pam(test_pam, test_ids, model).values()
            )

            test_decoy_scores = list(
                score_pam(test_decoy_pam, test_decoy_ids, model).values()
            )

            calibration_result = find_best_f1_threshold(
                calibration_scores=calibration_scores,
                decoy_scores=calibration_decoy_scores,
            )

            threshold = calibration_result["threshold"]

            test_result = confusion_at_threshold(
                positive_scores=test_scores,
                decoy_scores=test_decoy_scores,
                threshold=threshold,
            )

            model["status"] = "trained"

            model["split"] = {
                "total_positive_genomes": total_pos,
                "n_train": len(train_ids),
                "n_calibration_positive": len(calibration_ids),
                "n_calibration_decoy": len(calibration_decoy_ids),
                "n_test_positive": len(test_ids),
                "n_test_decoy": len(test_decoy_ids),
                "random_seed": args.random_seed + i,
            }

            model["classification"] = {
                "score_definition": model["score_definition"],
                "threshold_selection": "max_F1_on_calibration_positives_vs_equal_sized_decoys",
                "plausibility_threshold": threshold,
                "calibration": calibration_result,
                "test": test_result,
            }

            append_jsonl(args.models_jsonl, model)

            runtime = time.time() - start

            append_evaluation_row(args.evaluation_report, {
                "target_domain": domain,
                "status": "trained",
                "total_positive_genomes": total_pos,
                "n_train": len(train_ids),
                "n_calibration_positive": len(calibration_ids),
                "n_calibration_decoy": len(calibration_decoy_ids),
                "n_test_positive": len(test_ids),
                "n_test_decoy": len(test_decoy_ids),
                "n_initial_features": model["n_initial_features"],
                "n_features": model["n_features"],
                "max_model_features": model["max_model_features"],
                "threshold": threshold,
                "calibration_precision": calibration_result["precision"],
                "calibration_recall": calibration_result["recall"],
                "calibration_F1": calibration_result["F1"],
                "calibration_TP": calibration_result["TP"],
                "calibration_FP": calibration_result["FP"],
                "calibration_TN": calibration_result["TN"],
                "calibration_FN": calibration_result["FN"],
                "test_precision": test_result["precision"],
                "test_recall": test_result["recall"],
                "test_F1": test_result["F1"],
                "test_TP": test_result["TP"],
                "test_FP": test_result["FP"],
                "test_TN": test_result["TN"],
                "test_FN": test_result["FN"],
                "runtime_seconds": round(runtime, 3),
            })

            trained += 1

            print(
                f"    [CALIBRATION] threshold={threshold:.4f} | "
                f"F1={calibration_result['F1']:.3f} | "
                f"P={calibration_result['precision']:.3f} | "
                f"R={calibration_result['recall']:.3f} | "
                f"TP={calibration_result['TP']} FP={calibration_result['FP']} "
                f"TN={calibration_result['TN']} FN={calibration_result['FN']}"
            )

            print(
                f"    [TEST] F1={test_result['F1']:.3f} | "
                f"P={test_result['precision']:.3f} | "
                f"R={test_result['recall']:.3f} | "
                f"TP={test_result['TP']} FP={test_result['FP']} "
                f"TN={test_result['TN']} FN={test_result['FN']}"
            )

        except Exception as exc:
            skipped += 1
            runtime = time.time() - start

            append_jsonl(args.models_jsonl, {
                "target_domain": domain,
                "status": "error",
                "error": str(exc),
            })

            append_evaluation_row(args.evaluation_report, {
                "target_domain": domain,
                "status": "error",
                "runtime_seconds": round(runtime, 3),
                "error": str(exc),
            })

            print(f"    [ERROR] {domain}: {exc}")

        print(
            f"    [PROGRESS] trained={trained:,}, skipped={skipped:,}, "
            f"processed={i:,}/{len(domains):,}"
        )

    print("\n========== SUMMARY ==========")
    print(f"Domains processed: {len(domains):,}")
    print(f"Models trained:    {trained:,}")
    print(f"Models skipped:    {skipped:,}")
    print(f"Runtime:           {time.time() - start_all:.1f}s")
    print(f"Models JSONL:      {args.models_jsonl}")
    print(f"Evaluation report: {args.evaluation_report}")
    print("=============================")


if __name__ == "__main__":
    main()

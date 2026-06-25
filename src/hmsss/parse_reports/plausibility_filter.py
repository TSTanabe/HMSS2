from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable


def load_plausibility_models(
        model_jsonl_path: str | Path,
) -> dict[str, dict[str, Any]]:
    """
    Load positive-only plausibility models from a JSONL file.

    Returns
    -------
    dict
        {target_domain: model}
    """

    model_jsonl_path = Path(model_jsonl_path)

    models: dict[str, dict[str, Any]] = {}

    with open(model_jsonl_path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            record = json.loads(line)

            if record.get("status") not in {None, "trained"}:
                continue

            target_domain = record.get("target_domain")

            if not target_domain:
                continue

            if "feature_weights" not in record:
                continue

            if "classification" not in record:
                continue

            models[target_domain] = record

    return models


def get_model_threshold(model: dict[str, Any]) -> float | None:
    """
    Accepts both currently used model formats:
    - model["classification"]["threshold"]
    - model["classification"]["plausibility_threshold"]
    """

    classification = model.get("classification", {})

    if "threshold" in classification:
        return float(classification["threshold"])

    if "plausibility_threshold" in classification:
        return float(classification["plausibility_threshold"])

    return None


def collect_present_domains_from_protein_dict(
        protein_dict: dict[str, Any],
) -> set[str]:
    """
    Collect all domains present in the current protein_dict.

    Domains are collected from protein.domains, not from low_score_domains.
    This means the plausibility background reflects the currently selected
    best/non-overlapping domains.
    """

    present_domains: set[str] = set()

    for protein in protein_dict.values():
        for domain in getattr(protein, "domains", set()):
            domain_name = getattr(domain, "domain", None)

            if domain_name:
                present_domains.add(str(domain_name))

    return present_domains


def score_present_domains_against_model(
        present_domains: set[str],
        model: dict[str, Any],
) -> float:
    """
    Score a set of present domains against one positive-only support model.

    This matches the training-time scoring logic:
    weighted sum over model features that are present in the genome.
    """

    target_domain = model.get("target_domain")

    feature_weights = model.get("feature_weights", {})

    score = 0.0

    for feature, weight in feature_weights.items():
        if feature == target_domain:
            continue

        if feature in present_domains:
            score += float(weight)

    return float(score)


def iter_protein_domain_names(protein: Any) -> Iterable[str]:
    """
    Yield domain names from protein.domains.
    """

    for domain in getattr(protein, "domains", set()):
        domain_name = getattr(domain, "domain", None)

        if domain_name:
            yield str(domain_name)


def add_plausibility_comment_to_matching_domains(
        protein: Any,
        target_domain: str,
        comment: str = "Pp",
) -> None:
    """
    Add plausibility comment to all domains of this protein matching target_domain.
    """

    for domain in getattr(protein, "domains", set()):
        domain_name = getattr(domain, "domain", None)

        if domain_name != target_domain:
            continue

        if comment not in domain.selection_comment_list:
            domain.add_selection_comment(comment)


def apply_plausibility_thresholds_to_protein_dict(
        protein_dict: dict[str, Any],
        plausibility_models: dict[str, dict[str, Any]],
        *,
        comment: str = "Pt",
) -> dict[str, Any]:
    """
    Rescue only the best-scoring non-valid candidate per target domain.

    Rules:
    1. Only non-valid proteins are candidates.
    2. If any protein with the same target domain is already valid, skip that domain.
       This prevents plausibility-based validation of additional paralogs.
    3. For each remaining target domain, only the candidate protein/domain with the
       highest domain score is tested.
    4. A candidate is rescued only if a model exists and score >= threshold.
    """

    if not protein_dict or not plausibility_models:
        return protein_dict

    present_domains = collect_present_domains_from_protein_dict(protein_dict)

    if not present_domains:
        return protein_dict

    # Domains that already have at least one valid protein in this genome.
    already_valid_domains: set[str] = set()

    for protein in protein_dict.values():
        if not getattr(protein, "valid_hit", False):
            continue

        for domain_name in iter_protein_domain_names(protein):
            already_valid_domains.add(domain_name)

    # Best non-valid candidate per domain:
    # target_domain -> (protein, domain_score)
    best_candidate_by_domain: dict[str, tuple[Any, float]] = {}

    for protein in protein_dict.values():
        if getattr(protein, "valid_hit", False):
            continue

        for domain in getattr(protein, "domains", set()):
            target_domain = getattr(domain, "domain", None)

            if not target_domain:
                continue

            target_domain = str(target_domain)

            # Skip domains that already have a valid hit in this genome.
            # This prevents plausibility rescue of paralogs.
            if target_domain in already_valid_domains:
                continue

            # Skip if no plausibility model exists.
            if target_domain not in plausibility_models:
                continue

            score = float(getattr(domain, "score", 0.0))

            previous = best_candidate_by_domain.get(target_domain)

            if previous is None or score > previous[1]:
                best_candidate_by_domain[target_domain] = (protein, score)

    # Cache model decisions per domain.
    domain_decision_cache: dict[str, tuple[bool, float, float]] = {}

    for target_domain, (protein, domain_score) in best_candidate_by_domain.items():
        model = plausibility_models.get(target_domain)

        if model is None:
            continue

        threshold = get_model_threshold(model)

        if threshold is None:
            continue

        if target_domain not in domain_decision_cache:
            plausibility_score = score_present_domains_against_model(
                present_domains=present_domains,
                model=model,
            )

            is_plausible = plausibility_score >= threshold

            domain_decision_cache[target_domain] = (
                is_plausible,
                plausibility_score,
                threshold,
            )

        is_plausible, plausibility_score, threshold = domain_decision_cache[target_domain]

        if not is_plausible:
            continue

        protein.valid_hit = True

        add_plausibility_comment_to_matching_domains(
            protein=protein,
            target_domain=target_domain,
            comment=comment,
        )

        setattr(protein, "plausibility_rescued", True)

    return protein_dict

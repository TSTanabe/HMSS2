#!/usr/bin/python
from typing import Dict, List, Set

from hmsss.core.logging import get_logger

logger = get_logger(__name__)


def enhance_pathway_completeness(
    protein_dict: Dict[str, "Protein"],
    pattern_dict: Dict[str, tuple[set, int]],
    threshold_dict: Dict[str, float],
    *,
    mark_valid: bool = True,
    selection_comment: str = "Coo",
) -> Set[str]:
    """
    Identify proteins that contribute to fully complete pathway patterns.

    Logic:
    For each pattern in `pattern_dict`, this function checks:
      1. Is there at least one protein in the genome for *each* required domain
         that scores above or equal to the threshold for that domain?
         - (A pattern is considered "complete" only if this is true for *every* domain.)
      2. If so, gather *all* protein IDs that have any required domain with score >= its threshold
         (i.e., for each domain in the pattern, collect all proteins that fulfill the score criterion).

    If any required domain for a pattern has no matching protein above threshold,
    the pattern is ignored and does not contribute to the final set.

    Inputs:
    protein_dict : Dict[str, Protein]
        Mapping of protein IDs to Protein objects.

    pattern_dict : Dict[str, tuple[set, int]]
        Mapping of pattern names to lists of required domain names.
                "PatternB": ["X", "Y"] int is the pattern length

    threshold_dict : Dict[str, float]
                "A": 40.0,

    Outputs:
    Set[str]
        Set of all protein IDs that fulfill at least one required domain
        (with score >= threshold) for any *fully complete* pattern.
        Each protein ID is included at most once (set semantics).

    """
    found_protein_ids: Set[str] = set()
    domain_to_protein: Dict[str, List[tuple[str, float]]] = {}

    # Build mapping from domain name to all (protein_id, score) tuples in the genome
    for protein_id, protein in protein_dict.items():
        for domain in protein.get_domain_listing():
            domain_name = domain.get_domain()
            score = domain.get_score()
            domain_to_protein.setdefault(domain_name, []).append((protein_id, score))

    for required_domains, pattern_length in pattern_dict.values():
        domain_hits: Dict[str, Set[str]] = {}
        all_domains_present = True

        for domain in required_domains:
            threshold = threshold_dict.get(domain, 0)

            pairs = domain_to_protein.get(domain, [])
            if not pairs:
                all_domains_present = False
                break

            # hits above/equal threshold
            hits_above = {pid for pid, score in pairs if score >= threshold}

            if hits_above:
                domain_hits[domain] = hits_above
            else:
                # Fallback: pick the best highest-scoring hit to "complete" the pattern
                best_pid, best_score = max(pairs, key=lambda t: t[1])
                domain_hits[domain] = {best_pid}
                # print(f"Fallback: best_pid: {best_pid}, best_score: {best_score}")

        if not all_domains_present:
            continue

        # debug
        msg = f"Pattern complete (threshold+fallback-lowest): {', '.join(domain_hits.keys())}\n"
        for d, pids in domain_hits.items():
            msg += f"  {d}: {', '.join(sorted(pids))}\n"
        logger.debug(msg.rstrip())

        for domain_name, hits in domain_hits.items():
            for hit in hits:
                protein = protein_dict.get(hit)
                if protein is None:
                    continue
                protein.valid_hit = True
                protein.add_selection_comment(selection_comment)

    return found_protein_ids

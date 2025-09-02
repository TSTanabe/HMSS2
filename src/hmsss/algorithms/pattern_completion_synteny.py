#!/usr/bin/python

from collections import defaultdict
from typing import List, Iterable, Dict, Tuple, Any, Set, Optional

import numpy as np
from numpy._typing import NDArray
from scipy.optimize import linear_sum_assignment

from hmsss.algorithms.csb_finder import Keyword
from hmsss.core.logging import get_logger


logger = get_logger(__name__)


def _select_balanced_best_keywords(
    keywords: Iterable[Keyword],
    pattern_dict: dict[str, tuple[List[str], int]],
    min_length: int,
    *,
    strict_equal: bool = False,
) -> List[Keyword]:
    """
    Wähle die besten Keywords in einem Pass:
      - Balance: n_missing > 0 und (== n_additional, falls strict_equal) sonst <= n_additional.
      - Minimales n_missing.
      - Bei Gleichstand: maximale pattern_length.
      - Nur Keywords, deren pattern_length > min_length.

    Args:
        keywords: Iterable von Keyword-Objekten.
        pattern_dict: Map von keyword_id -> (Pattern-Set, Pattern-Länge).
        min_length: Untergrenze: nur Patterns mit Länge > min_length werden berücksichtigt.
        strict_equal: Wenn True, ist Balance n_missing == n_additional (statt <=).

    Returns:
        Liste der besten Keyword-Objekte (kann mehrere bei Gleichstand enthalten).
    """
    best_fitting_keywords: List["Keyword"] = []
    best_missing: int = 10**9
    best_len: int = -1

    for kw in keywords:
        pat = pattern_dict.get(kw.keyword_id)
        if not pat:
            continue

        if not kw.completeness > 0.5:
            continue

        _, pattern_length = pat
        if pattern_length <= min_length:
            continue

        n_missing = len(kw.get_missing_domains())
        if n_missing == 0:
            continue

        n_additional = len(kw.get_additional_domains())
        balanced = (
            (n_missing == n_additional) if strict_equal else (n_missing <= n_additional)
        )
        if not balanced:
            continue

        if n_missing < best_missing:
            best_missing = n_missing
            best_len = pattern_length
            best_fitting_keywords = [kw]
        elif n_missing == best_missing:
            if pattern_length > best_len:
                best_len = pattern_length
                best_fitting_keywords = [kw]
            elif pattern_length == best_len:
                best_fitting_keywords.append(kw)

    return best_fitting_keywords


def _find_best_keywords(cluster_dict: Dict[str, Any], pattern_dict: Dict[str, Any]):
    best_keywords_dict: Dict[str, Any] = {}
    for cluster_id, cluster in cluster_dict.items():
        keywords = cluster.get_keywords()
        max_length = 0
        # get the longest complete pattern and return length
        for keyword in keywords:
            if keyword.get_completeness() == 1:
                pattern_set, pattern_length = pattern_dict[keyword.keyword_id]
                if pattern_length > max_length:
                    max_length = pattern_length

        best_keywords = _select_balanced_best_keywords(
            keywords, pattern_dict, max_length
        )
        best_keywords_dict[cluster_id] = best_keywords
    return best_keywords_dict


def _find_possible_transitions(
    protein_id: str,
    protein: object,
    current_domain: str,
    missing_domains: Set[str],
    transition_dict: Dict[str, Set[Tuple[str, float]]],
) -> Dict[str, Set[Tuple[str, float]]]:
    """
    Find possible transitions from a current domain to any missing domains
    using previously deleted domains of a given protein.

    For the given `protein_id`, this function:
      1. Determines the score of the current domain.
      2. Iterates over the protein's `deleted_domains`.
      3. If a deleted domain is in the set of missing domains,
         it computes the absolute score difference to the current domain
         and records a possible transition.

    Args:
        protein_id (str): Identifier of the protein being evaluated.
        protein (Protein): Protein object containing active and deleted domains.
        current_domain (str): Domain name currently assigned to the protein.
        missing_domains (Set[str]): Domains expected but currently missing.
        transition_dict (Dict[str, Set[Tuple[str, float]]]):
            Mapping of missing_domain -> set of (protein_id, score_difference).

    Returns:
        Dict[str, Set[Tuple[str, float]]]:
            Updated transition_dict including any newly identified transitions.
    """
    current_score = 0.0
    for domain in protein.domains.values():
        if domain.get_domain() == current_domain:
            current_score = domain.get_score()

    for domain_name, domain_obj in protein.deleted_domains.items():
        if domain_name in missing_domains:
            score = domain_obj.get_score()
            difference = abs(current_score - score)
            transition_dict[domain_name].add((protein_id, difference))

    return transition_dict


def _solve_assignment(
    cost_matrix: NDArray[np.float_],
    cluster_id: Optional[str] = None,
) -> Tuple[Optional[NDArray[np.int_]], Optional[NDArray[np.int_]]]:
    """
    Solve a (possibly sparse) assignment problem using the Hungarian algorithm.

    The function first removes rows/columns that are entirely infeasible (all np.inf),
    solves the reduced problem, and then maps the assignment indices back to the
    original matrix coordinates. If the reduced problem is still infeasible or an
    error occurs, (None, None) is returned.

    Args:
        cost_matrix: 2D array of shape (n_rows, n_cols) with non-negative costs.
                     Use np.inf to indicate impossible assignments.
        cluster_id:  Optional identifier used to enrich log messages.

    Returns:
        (row_indices, col_indices):
            - row_indices: indices into the original rows (or None if infeasible)
            - col_indices: indices into the original columns (or None if infeasible)
    """


    # Mask rows/cols that contain at least one feasible entry
    valid_rows = ~np.all(np.isinf(cost_matrix), axis=1)
    valid_cols = ~np.all(np.isinf(cost_matrix), axis=0)
    reduced = cost_matrix[np.ix_(valid_rows, valid_cols)]

    def _matrix_str(matrix: NDArray[np.float_]) -> str:
        with np.printoptions(precision=2, suppress=True, linewidth=120, nanstr="nan", infstr="inf"):
            return "\n" + "\n".join(" ".join(f"{x:7}" for x in row) for row in matrix)

    # Infeasible if empty or any remaining row/col is still all inf
    if (
        reduced.size == 0
        or np.any(np.all(np.isinf(reduced), axis=1))
        or np.any(np.all(np.isinf(reduced), axis=0))
    ):
        logger.debug(
            f"[Assignment] Infeasible cost matrix for cluster {cluster_id or ''}:{_matrix_str(cost_matrix)}"
        )
        logger.error(
            f"Assignment failed for cluster {cluster_id or ''}: infeasible after row/col removal."
        )
        return None, None

    try:
        r_red, c_red = linear_sum_assignment(reduced)
        # Map reduced indices back to original coordinates
        row_idx = np.where(valid_rows)[0][r_red].astype(np.int_)
        col_idx = np.where(valid_cols)[0][c_red].astype(np.int_)
        return row_idx, col_idx
    except Exception as exc:
        logger.debug(
            f"[Assignment] Exception on cost matrix for cluster {cluster_id or ''}:{_matrix_str(cost_matrix)}"
        )
        logger.error(f"Assignment error for cluster {cluster_id or ''}: {exc}.")
        return None, None


def _get_optimal_transitions(
    transition_dict,
    missing_domains,
    initial_completeness=0.0,
    total_domains=None,
):
    """
    Calculates the optimal set of transitions to cover as many missing domains as possible,
    each proteinID at most once, minimizing transitions and total score difference.

    transition_dict: {missing_domain: set of (proteinID, score_diff)}
    missing_domains: list or set of missing domains to fulfill

    Returns:
        chosen_transitions: [(proteinID, to_domain, score_diff)]
        completeness: number fulfilled / total
        total_score_diff: sum of chosen score diffs
    """
    # Gather all protein candidates
    proteins = set()
    for domain in missing_domains:
        for protein, diff in transition_dict[domain]:
            proteins.add(protein)
    proteins = list(proteins)
    n_proteins = len(proteins)
    n_domains = len(missing_domains)

    if n_proteins == 0 or n_domains == 0:
        return [], 0.0, 0.0

    # Build cost matrix
    cost_matrix = np.full((n_proteins, n_domains), np.inf)
    for j, domain in enumerate(missing_domains):
        for protein, diff in transition_dict[domain]:
            i = proteins.index(protein)
            cost_matrix[i, j] = diff

    # Hungarian assignment: minimize score difference (and thus minimize transitions)
    orig_row_ind, orig_col_ind = _solve_assignment(cost_matrix)

    if orig_row_ind is None:
        # Matching impossible, keep the cluster as assigned
        return [], initial_completeness, 0.0

    chosen_transitions = []
    total_score_diff = 0
    fulfilled_domains = set()
    used_proteins = set()
    missing_domains = list(missing_domains)

    for i, j in zip(orig_row_ind, orig_col_ind):
        cost = cost_matrix[i, j]
        if np.isfinite(cost):
            protein_id = proteins[i]
            domain = missing_domains[j]
            chosen_transitions.append((protein_id, domain, cost))
            total_score_diff += cost
            fulfilled_domains.add(domain)
            used_proteins.add(protein_id)

    # Completeness for the final keyword
    final_completeness = initial_completeness
    if total_domains:
        final_completeness += len(fulfilled_domains) / total_domains
    else:
        final_completeness += len(fulfilled_domains) / max(len(missing_domains), 1)

    return chosen_transitions, final_completeness, total_score_diff


def get_all_optimized_transitions(
    best_keywords_dict: Dict[str, Any],
    cluster_dict: Dict[str, Any],
    protein_dict: Dict[str, Any],
    pattern_dict: Dict[str, Any],
):
    transitions_per_cluster = {}
    for cluster_id, best_keywords in best_keywords_dict.items():
        possible_transitions = {}
        for keyword in best_keywords:
            completeness = keyword.get_completeness()
            missing_domains = keyword.get_missing_domains()
            additional_domains = keyword.get_additional_domains()
            transition_dict = defaultdict(set)
            # Structure: alternative_protein_type_dict[(hit_proteinID, query)] = (hit_proteinID, query, hsp_start, hsp_end, hit_bitscore, genomeID)
            # Structure: transition_dict[query] => ((hit_proteinID, difference), (hit_proteinID, difference))
            # print(f"\nNew gene cluster {cluster.clusterID}")
            # print(types)
            # print(f"Processing {keyword.keyword} {completeness}")
            # print(f"Missing {missing_domains}")
            # print(f"Additional {additional_domains}")

            # For every "additional" domain: get proteinID and possible transitions to "missing"
            cluster = cluster_dict[cluster_id]
            genes = cluster.get_genes()
            domains = cluster.get_domains()
            for current_domain, protein_id in zip(domains, genes):
                if (
                    current_domain in additional_domains
                    # and protein_id in intermediate_protein_dict
                    # and "Tc" not in protein.selection_comment
                ):
                    protein = protein_dict[protein_id]
                    _find_possible_transitions(
                        protein_id,
                        protein,
                        current_domain,
                        missing_domains,
                        transition_dict,
                    )

            # Optimize the transitions by minimizing the bitscore changes and number of transitions to reach the missing domains
            pattern, pattern_length = pattern_dict[keyword.keyword_id]

            chosen_transitions, posterior_completeness, total_score_diff = (
                _get_optimal_transitions(
                    transition_dict, missing_domains, completeness, pattern_length
                )
            )
            """
            Example output
            chosen_transitions = [
                ('P1', 'A', 50),
                ('P2', 'B', 20),
                ('P3', 'C', 40)
            ]
            completeness = 1.0
            total_score_diff = 110
            """

            if chosen_transitions:
                possible_transitions[
                    (pattern_length, posterior_completeness, total_score_diff)
                ] = chosen_transitions

        if possible_transitions:
            transitions_per_cluster[cluster_id] = (
                possible_transitions
            )

    return transitions_per_cluster


def execute_pattern_completion(transition_dict, protein_dict):
    for clusterID, optimal_transitions_per_keyword in transition_dict.items():
        # Sort by: pattern_length (desc), posterior_completeness (desc), total_score_diff (asc)

        # Skip empty transitions
        if not optimal_transitions_per_keyword:
            continue

        best_key = max(
            optimal_transitions_per_keyword.keys(),
            key=lambda x: (
                x[1],
                x[0],
                -x[2],
            ),  # pattern_length, posterior_completeness, -score_diff
        )
        best_transitions = optimal_transitions_per_keyword[best_key]

        logger.debug(
            f"For clusterID {clusterID} with pattern length, completeness and score difference {best_key} following conversion is done"
        )
        logger.debug(best_transitions)

        # For the best transition alter the proteins domain information
        for proteinID, to_domain, _ in best_transitions:
            if proteinID in protein_dict:
                protein = protein_dict[proteinID]

                # Save the original hit as comment
                original_domains = protein.get_domains()
                protein.add_selection_comment("Syc")
                protein.alternative_hit = original_domains

                # Get alternative domain
                new_domain = protein.deleted_domains[to_domain]

                protein.add_domain(
                    new_domain.domain,
                    new_domain.start,
                    new_domain.end,
                    new_domain.score,
                    force=True,
                )


def enhance_syntenic_block_completeness(
    cluster_dict,
    combined_protein_dict,
    pattern_dict,
):
    """
    Enhance syntenic block completeness by swapping additional protein domains with missing ones if possible.
    for each cluster that has no directly matching pattern from the given patterns
    it is tested if a possible conversion of protein types to alternative ones with lower hitscore could reach
    a better completion

    This can overwrite hit above trusted cutoff/reference sequence hits if the completion is better
    These transitions are displayed in the database and the hit report output

    Keywords that are used to find a possible transition have either same number of additionals to missing
    or more additionals than missing.
    The transition with the highest pattern lenght, highest final completion and lowest lost score is
    chosen for the final transition selection
    """
    best_keywords_dict = _find_best_keywords(cluster_dict, pattern_dict)

    # Optional task filter out keywords that do not comply with allowed and possible transitions

    optimized_executable_transitions_dict = get_all_optimized_transitions(
        best_keywords_dict, cluster_dict, combined_protein_dict, pattern_dict
    )

    execute_pattern_completion(
        optimized_executable_transitions_dict, combined_protein_dict
    )

    return combined_protein_dict

#!/usr/bin/python
import os
import re
import subprocess
import traceback
from collections import defaultdict
from multiprocessing import Pool, Manager
from typing import Dict, Any, List, Set

from contextlib import contextmanager
from time import perf_counter

import numpy as np
from scipy.optimize import linear_sum_assignment

from hmsss.algorithms import csb_finder
from hmsss.algorithms import search_cross_reference
from hmsss.core.logging import get_logger
from hmsss.db import database
from hmsss.io import output

logger = get_logger(__name__)


def enhance_pathway_completeness(
    protein_dict: Dict[str, "Protein"],
    pattern_dict: Dict[str, tuple[set, int]],
    threshold_dict: Dict[str, float],
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

    pattern_dict : Dict[str, List[str]]
        Mapping of pattern names to lists of required domain names.
                "PatternB": ["X", "Y"]

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
        all_domains_above = True
        for domain in required_domains:
            threshold = threshold_dict.get(domain, 0)
            hits = {
                protein_id
                for protein_id, score in domain_to_protein.get(domain, [])
                if score >= threshold
            }
            if not hits:
                all_domains_above = False
                break
            domain_hits[domain] = hits

        if all_domains_above:
            # Format readable block
            msg = f"Pattern complete: {', '.join(domain_hits.keys())}\n"
            for domain in domain_hits:
                proteins = ", ".join(sorted(domain_hits[domain]))
                msg += f"  {domain}: {proteins}\n"
            logger.debug(msg.rstrip())
            for hits in domain_hits.values():
                found_protein_ids.update(hits)

    return found_protein_ids

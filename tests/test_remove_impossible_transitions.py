# test_remove_impossible_transitions.py
from typing import Any, Dict, List, Tuple
import pytest

# import the function under test (FIX THIS PATH to your module!)
from hmsss.algorithms.pattern_completion_synteny import _set_possible_transitions

# bring in the mask type from your transition_mask module
from hmsss.algorithms.transition_mask import TransitionMasks
from hmsss.algorithms import transition_mask


class KeywordStub:
    """Minimal stand-in for your keyword objects."""

    def __init__(self, kid: str, additional: Tuple[str, ...], missing: Tuple[str, ...]):
        self.keyword_id = kid
        self.additional_domains = additional
        self.missing_domains = missing

    def __repr__(self) -> str:
        return f"KW({self.keyword_id}, add={self.additional_domains}, miss={self.missing_domains})"


def test_remove_impossible_transitions_filters(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    # --- Input keywords grouped by cluster ---
    keyword_dict: Dict[str, List[Any]] = {
        "clusterA": [
            KeywordStub("kw1", ("A", "B"), ("X",)),  # we will mark this as possible
            KeywordStub("kw2", ("C",), ("Y",)),  # impossible
        ],
        "clusterB": [
            KeywordStub("kw3", ("B",), ("X", "Z")),  # possible
            KeywordStub("kw4", ("D",), ("Q",)),  # impossible
        ],
    }

    # --- Build a minimal mask (content not used by our mock) ---
    masks = TransitionMasks(
        missing_id={"X": 0, "Y": 1, "Z": 2, "Q": 3}, allow_mask={"A": 1, "B": 2}
    )

    # --- Decide which (additional, missing) tuple pairs should be considered "possible" ---
    # We key by the exact tuples that KeywordStub carries.
    POSSIBLE = {
        (("A", "B"), ("X",)),  # kw1 (clusterA) -> possible
        (("B",), ("X", "Z")),  # kw3 (clusterB) -> possible
        # all others -> False (impossible)
    }

    def mock_can_cover(additional, missing, mask) -> bool:
        # Ensure tuple semantics for mapping (your routine already passes tuples)
        return (tuple(additional), tuple(missing)) in POSSIBLE

    # Monkeypatch transition_mask.can_cover used by the routine
    monkeypatch.setattr(transition_mask, "can_cover", mock_can_cover)

    # --- Run the routine ---
    out = _set_possible_transitions(keyword_dict, masks)

    # --- Print comparison (as requested) ---
    print("INPUT:")
    for cid, kws in keyword_dict.items():
        print(f"  {cid}: {[kw.keyword_id for kw in kws]}")

    print("\nOUTPUT (filtered):")
    for cid, kws in out.items():
        print(f"  {cid}: {[kw.keyword_id for kw in kws]}")

    # Flush prints to test output
    captured = capsys.readouterr()
    print(captured.out)  # so it's visible when you run pytest -s

    # --- Assertions ---
    # Only kw1 should remain in clusterA; only kw3 in clusterB
    assert set(out.keys()) == {"clusterA", "clusterB"}
    assert [kw.keyword_id for kw in out["clusterA"]] == ["kw1"]
    assert [kw.keyword_id for kw in out["clusterB"]] == ["kw3"]

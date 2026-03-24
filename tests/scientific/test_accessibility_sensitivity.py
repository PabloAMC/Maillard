"""
tests/scientific/test_accessibility_sensitivity.py

P1 gate tests: verify that changing denaturation state (0.3 → 0.8) on pea_iso
produces a measurable change in predictions, and that AccessibilityState profiles
are correctly assigned and surfaced.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_correction import classify_accessibility_state, AccessibilityState
from src.usability_reports import DomainOfValidityChecker
from src.pipeline import MaillardPipeline
from src.conditions import ReactionConditions


# ---------------------------------------------------------------------------
# Unit-level: AccessibilityState profile mappings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("denaturation,expected_profile,expect_warning", [
    (0.10, "protein_embedded", True),
    (0.15, "protein_embedded", True),
    (0.30, "peptide_bound", True),
    (0.55, "peptide_bound", True),
    (0.70, "partially_opened", False),
    (0.88, "partially_opened", False),
    (0.95, "free_like", False),
])
def test_accessibility_state_profiles(denaturation, expected_profile, expect_warning):
    state = classify_accessibility_state("pea_iso", denaturation)
    assert state.profile == expected_profile, (
        f"denaturation={denaturation}: expected profile={expected_profile!r}, got {state.profile!r}"
    )
    assert state.accessibility_warning == expect_warning, (
        f"denaturation={denaturation}: expected warning={expect_warning}, got {state.accessibility_warning}"
    )


def test_free_protein_type_always_free_like():
    state = classify_accessibility_state("free", 0.1)
    assert state.profile == "free_like"
    assert not state.accessibility_warning


def test_accessibility_warning_triggers_domain_warning():
    """Embedded state on pea_iso should produce a CAUTION DomainWarning."""
    checker = DomainOfValidityChecker("meaty")
    warnings = checker.check(
        precursor_names=["ribose", "cysteine"],
        protein_type="pea_iso",
        temp_c=100.0,
        ph=6.5,
        aw=0.95,
        matrix_explainability={
            "accessibility_profile": "protein_embedded",
            "accessibility_warning": True,
            "accessibility_dominant_source": "estimated_from_conditions",
        },
    )
    accessibility_warnings = [w for w in warnings if w.category == "ACCESSIBILITY"]
    assert accessibility_warnings, "Expected an ACCESSIBILITY DomainWarning for an embedded pea_iso state"


# ---------------------------------------------------------------------------
# Integration: denaturation sensitivity produces measurable ppb difference
# ---------------------------------------------------------------------------

def _run_pea_predictions(denaturation_state: float) -> dict:
    """Run the pea_iso pipeline with a given denaturation state."""
    conditions = ReactionConditions(
        pH=6.5,
        temperature_celsius=100.0,
        water_activity=0.95,
        protein_type="pea_iso",
    )
    formulation = {
        "name": f"pea_iso_denat_{denaturation_state}",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine"],
        "lipids": [],
        "additives": [],
        "molar_ratios": {},
        "ph": 6.5,
        "temp": 100.0,
        "aw": 0.95,
        "time_minutes": 45.0,
        "protein_type": "pea_iso",
        "denaturation_state": denaturation_state,
        "catalyst": None,
    }
    designer = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    result = designer.evaluate_single(formulation, conditions)
    return result


def test_denaturation_rerankings_pea_isolate():
    """
    P1 exit gate: denaturation 0.3 → 0.8 on pea_iso must produce either:
    (a) a change in the top-ranked observable compound, OR
    (b) >20% change in the top compound's observable_ppb.

    This proves that AccessibilityState actually changes predictions,
    not just labels them.
    """
    result_low = _run_pea_predictions(denaturation_state=0.3)
    result_high = _run_pea_predictions(denaturation_state=0.8)

    # Get top ppb for each run
    def top_compound_and_ppb(result):
        ranked = sorted(
            result.projection_metadata.items(),
            key=lambda item: float(item[1].get("observable_ppb", 0.0)),
            reverse=True,
        )
        if not ranked:
            return None, 0.0
        canon, meta = ranked[0]
        return str(meta.get("compound", canon)), float(meta.get("observable_ppb", 0.0))

    top_compound_low, ppb_low = top_compound_and_ppb(result_low)
    top_compound_high, ppb_high = top_compound_and_ppb(result_high)

    # Check accessibility profiles differ
    matrix_low = result_low.matrix_explainability or {}
    matrix_high = result_high.matrix_explainability or {}
    profile_low = matrix_low.get("accessibility_profile", "unknown")
    profile_high = matrix_high.get("accessibility_profile", "unknown")

    assert profile_low != profile_high, (
        f"Expected different accessibility profiles for denat=0.3 vs 0.8, "
        f"but both are: {profile_low!r}"
    )

    # Check top compound re-ranks OR ppb changes significantly
    top_differs = top_compound_low != top_compound_high
    if not top_differs:
        # Must have measurable ppb change
        reference = max(ppb_low, ppb_high)
        if reference > 0.0:
            relative_change = abs(ppb_high - ppb_low) / reference
            assert relative_change > 0.20, (
                f"Top compound is the same ({top_compound_low}) and ppb change is only "
                f"{relative_change:.1%} (<20%). Accessibility modeling appears to have no effect. "
                f"denat=0.3 → ppb={ppb_low:.4f}, denat=0.8 → ppb={ppb_high:.4f}"
            )

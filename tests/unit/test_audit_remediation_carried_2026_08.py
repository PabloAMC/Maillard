"""Regression tests for the three carried audit-remediation fixes of 2026-08-27.

See tasks/audit_remediation.md:

A. Campaign leaderboard roll-up printed pH 0.00 / Temp 0.0 (and Protein "free")
   for every row because it read the conditions off the per-run formulation
   dict, which carries none.
B. The reaction-network diagram's colour map was keyed on the pre-Wave-G1
   reaction-family names, so real chemistry drew in the grey default colour.
C. Four genuinely DOI-less sources kept their identifiers in `doi` /
   `source_doi` fields. They now carry a typed `identifier` +
   `identifier_scheme` pair, and the evidence classification must not move —
   above all the Liu hold-out bundle's `external_validation_only`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from src.benchmark_validation import (  # noqa: E402
    _matrix_external_data_status,
    _matrix_source_origin,
    matrix_source_anchor,
)
from src.input_normalization import resolve_condition_float, resolve_condition_value  # noqa: E402
from src.reporting import build_campaign_leaderboard  # noqa: E402


# --------------------------------------------------------------------------
# A. Campaign leaderboard condition resolution
# --------------------------------------------------------------------------


class _StubResult:
    """Minimal stand-in for FormulationResult as the leaderboard reads it."""

    def __init__(self, name: str, target_score: float) -> None:
        self.name = name
        self.target_score = target_score
        self.off_flavour_risk = 1.0
        self.safety_score = 2.0
        self.confidence_metadata = {"tier": "low", "prediction_mode": "directional_only"}


@pytest.mark.parametrize(
    "shared",
    [
        # the spelling data/campaigns/*.yml uses
        {"ph": 5.5, "temp": 105.0, "protein_type": "pea_iso"},
        # the spelling a scientist types / Formulation dataclass names
        {"pH": 5.5, "temperature": 105.0, "protein": "pea_iso"},
        # the benchmark-payload spelling
        {"ph": 5.5, "temp_C": 105.0, "protein_type": "pea_iso"},
    ],
)
def test_leaderboard_resolves_every_condition_key_spelling(shared):
    # Grid formulations carry NO process conditions — that is the whole bug.
    conditions_list = [{"name": "Cysteine Enrichment (Basic)", "sugars": ["ribose"]}]
    rows = build_campaign_leaderboard(
        results=[_StubResult("Cysteine Enrichment (Basic)", 4.4)],
        conditions_list=conditions_list,
        shared_conditions=shared,
    )
    assert rows[0]["ph"] == pytest.approx(5.5)
    assert rows[0]["temp"] == pytest.approx(105.0)
    assert rows[0]["protein_type"] == "pea_iso"


def test_leaderboard_prefers_per_run_override_then_shared_then_effective():
    results = [_StubResult("overridden", 2.0), _StubResult("inherits", 1.0)]
    conditions_list = [
        {"name": "overridden", "ph": 8.0, "temp": 180.0},
        {"name": "inherits"},
    ]
    rows = {
        row["name"]: row
        for row in build_campaign_leaderboard(
            results=results,
            conditions_list=conditions_list,
            shared_conditions={"ph": 5.5},
            effective_conditions={"ph": 6.0, "temp": 150.0, "protein_type": "free"},
        )
    }
    # per-run override wins
    assert rows["overridden"]["ph"] == pytest.approx(8.0)
    assert rows["overridden"]["temp"] == pytest.approx(180.0)
    # shared conditions beat the global fallback, which supplies what is missing
    assert rows["inherits"]["ph"] == pytest.approx(5.5)
    assert rows["inherits"]["temp"] == pytest.approx(150.0)


def test_leaderboard_reports_unknown_condition_as_none_not_zero():
    # The original defect: an unresolvable pH printed as 0.00, which is a
    # chemistry claim. Unknown must stay unknown.
    rows = build_campaign_leaderboard(
        results=[_StubResult("bare", 1.0)],
        conditions_list=[{"name": "bare"}],
    )
    assert rows[0]["ph"] is None
    assert rows[0]["temp"] is None
    assert rows[0]["ph"] != 0.0


def test_leaderboard_is_ranked_by_target_score():
    rows = build_campaign_leaderboard(
        results=[_StubResult("low", 0.1), _StubResult("high", 4.4)],
        conditions_list=[{"name": "low"}, {"name": "high"}],
        shared_conditions={"ph": 6.0, "temp": 120.0},
    )
    assert [row["name"] for row in rows] == ["high", "low"]


def test_condition_resolver_handles_formulation_objects():
    from src.formulation import Formulation

    form = Formulation(name="typed", ph=7.2, temperature=133.0)
    assert resolve_condition_float("ph", (form,)) == pytest.approx(7.2)
    assert resolve_condition_float("temp", (form,)) == pytest.approx(133.0)
    assert resolve_condition_value("protein_type", (form,)) == "free"


def test_condition_resolver_ignores_uncoercible_values():
    assert resolve_condition_float("ph", ({"ph": "not-a-number"},), default=None) is None


# --------------------------------------------------------------------------
# B. Reaction-network colour map
# --------------------------------------------------------------------------


def test_every_curated_reaction_family_has_an_edge_colour():
    import generate_reaction_network as network

    missing = network.uncoloured_reaction_families()
    assert not missing, (
        "reaction-family colour map is stale; no colour for "
        f"{sorted(missing)} — they would draw in {network.DEFAULT_EDGE_COLOR}"
    )


def test_colour_map_carries_no_family_the_engine_stopped_emitting():
    import generate_reaction_network as network

    stale = network.stale_colour_map_keys()
    assert not stale, f"colour map keys no curated pathway emits: {sorted(stale)}"


def test_colour_map_colours_are_distinct():
    import generate_reaction_network as network

    colours = list(network.FAMILY_EDGE_COLORS.values())
    assert len(set(colours)) == len(colours)
    assert network.DEFAULT_EDGE_COLOR not in colours


# --------------------------------------------------------------------------
# C. Typed identifiers for DOI-less sources
# --------------------------------------------------------------------------

HOLDOUT_BUNDLE = (
    ROOT
    / "data"
    / "benchmarks"
    / "external_validation"
    / "external_validation_liu_2023_ppi_offnote_baseline.json"
)
INTAKE_REGISTRY = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

# id -> expected identifier scheme, for the four sources that have no DOI.
DOI_LESS_REGISTRY_SOURCES = {
    "huang_2022_thiamine_metal_catalysis": "url",
    "uspto_ptacts_2023_yeast_extract_anchor": "patent",
    "cantre_2007_corned_beef_furosine": "journal_locator",
}


def test_doi_less_registry_sources_carry_a_typed_identifier_pair():
    payload = json.loads(INTAKE_REGISTRY.read_text(encoding="utf-8"))
    by_id = {row["id"]: row for row in payload["eligible_references"]}
    for source_id, scheme in DOI_LESS_REGISTRY_SOURCES.items():
        row = by_id[source_id]
        assert "doi" not in row, f"{source_id} still stores a non-DOI in a doi field"
        assert row["identifier"], f"{source_id} has no identifier"
        assert row["identifier_scheme"] == scheme
        assert row["identifier_note"], f"{source_id} must record why it has no DOI"


def test_holdout_bundle_carries_typed_identifier_and_keeps_its_evidence_class():
    bench = json.loads(HOLDOUT_BUNDLE.read_text(encoding="utf-8"))
    assert "source_doi" not in bench
    assert bench["identifier_scheme"] == "citation"
    assert bench["identifier"]
    # The invariant that matters: retyping the identifier must not touch the
    # hold-out classification.
    assert bench["evidence_class"] == "external_validation_only"
    assert str(_matrix_external_data_status(bench)) == "external_validation_only"


def test_typed_identifier_counts_as_an_external_anchor():
    # A DOI-less external source must not be downgraded to "unspecified origin"
    # just because its identifier moved out of a doi-named field.
    typed = {
        "benchmark_id": "typed_identifier_probe",
        "identifier": "https://example.invalid/thesis/1",
        "identifier_scheme": "url",
        "measured_volatiles": {"Hexanal": {"conc_ppb": 100.0}},
    }
    assert matrix_source_anchor(typed) == "https://example.invalid/thesis/1"
    assert _matrix_source_origin(typed) == "external_literature"
    assert str(_matrix_external_data_status(typed)) == "external_quantitative"

    # ... and a payload with neither is still unspecified.
    bare = {"benchmark_id": "bare", "measured_volatiles": {"Hexanal": {"conc_ppb": 1.0}}}
    assert matrix_source_anchor(bare) == ""
    assert _matrix_source_origin(bare) == "unspecified"


def test_source_doi_still_wins_when_present():
    bench = {"source_doi": "10.1016/j.foodchem.2022.134998", "identifier": "fallback"}
    assert matrix_source_anchor(bench) == "10.1016/j.foodchem.2022.134998"


def test_citation_gate_baseline_is_empty():
    import importlib

    gate = importlib.import_module("ci.citation_gate")
    assert gate.WAIVERS == (), (
        "the citation-gate waiver baseline ratchets down only; re-adding a waiver "
        "needs a dated justification in tasks/audit_remediation.md"
    )
    blocking, waived, stale, _dois, files_scanned = gate.run_offline()
    assert files_scanned > 0
    assert not blocking, [str(v) for v in blocking]
    assert not waived
    assert not stale

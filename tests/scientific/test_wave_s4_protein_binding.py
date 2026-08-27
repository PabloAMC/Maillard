"""Wave S4 — measured protein-binding physics as a matrix observability mode.

2026-08-27.  See tasks/audit_remediation.md, "Wave S4".

WHAT THESE TESTS DO AND DO NOT PIN.

They pin (a) that the shipped prediction path is UNCHANGED unless a caller explicitly
selects another mode, (b) the arithmetic of the binding model, (c) the provenance
discipline of `data/lit/binding_constants.yml`, and (d) that the no-double-count guard
actually fires.

They deliberately do NOT pin the mode-vs-mode SCORES.  Pinning a number this wave just
produced is the circularity Rounds 1-3 of this audit removed, and Wave U set the
precedent explicitly.  The scores live in
`results/validation/matrix_binding_mode_comparison.{json,md}`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from src import protein_binding as pb
from src.benchmark_validation import _run_matrix_only_benchmark_prediction

ROOT = Path(__file__).resolve().parents[2]
PEA_BENCH = ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json"
LIU_BENCH = (
    ROOT
    / "data"
    / "benchmarks"
    / "external_validation"
    / "external_validation_liu_2023_ppi_offnote_baseline.json"
)


def _bench(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------------------
# (1) The shipped path is inert unless a mode is explicitly selected.
# --------------------------------------------------------------------------------------

def test_default_mode_is_the_fitted_incumbent():
    assert pb.observability_mode() == pb.MODE_FITTED
    assert pb.binding_mode_active() is False


def test_shipped_pea_prediction_is_unchanged_by_this_wave():
    """The wave must not move a shipped number. These are the pre-Wave-S4 values.

    1125.278 ppb hexanal is `oxidation load x base marker yield x 4.31725`, the Wave O
    refit constant. If this moves, Wave S4 has leaked into the default path.
    """
    result = _run_matrix_only_benchmark_prediction(_bench(PEA_BENCH))
    assert result["predicted_ppb"]["Hexanal"] == pytest.approx(1125.278, rel=1e-4)
    assert result["predicted_ppb"]["2-Pentylfuran"] == pytest.approx(638.267, rel=1e-4)
    assert result["projection_metadata"]["Hexanal"]["total_observable_factor"] == pytest.approx(
        4.31725, rel=1e-6
    )
    assert result["projection_metadata"]["Hexanal"]["observability_mode"] == pb.MODE_FITTED


def test_mode_context_manager_restores_the_default():
    with pb.use_observability_mode(pb.MODE_BINDING):
        assert pb.observability_mode() == pb.MODE_BINDING
    assert pb.observability_mode() == pb.MODE_FITTED
    with pytest.raises(ValueError):
        with pb.use_observability_mode("not_a_mode"):
            pass


# --------------------------------------------------------------------------------------
# (2) The model's arithmetic, in explicit units.
# --------------------------------------------------------------------------------------

def test_free_fraction_is_the_single_site_langmuir_in_declared_units():
    """f_free = 1 / (1 + a_p * Pow * c_p), with a_p in L per gram of protein.

    Computed here from the YAML's own numbers rather than a hard-coded expectation, so
    the test checks the FORM and the unit handling; if a constant is ever edited the
    expectation moves with it and the arithmetic is still verified.
    """
    ap = pb.hydrophobic_interaction_parameter("pea_iso", "ketone")
    pow_entry = pb.octanol_water_partition("hexanal")
    assert ap is not None and pow_entry is not None
    assert ap["units"] == "L_per_g_protein"

    c_p = 100.0
    expected = 1.0 / (1.0 + float(ap["value"]) * float(pow_entry["pow"]) * c_p)
    result = pb.free_fraction_from_ap(
        "hexanal", protein_source="pea_iso", protein_concentration_g_per_l=c_p
    )
    assert result.f_free == pytest.approx(expected, rel=1e-12)
    assert result.k_eff_l_per_g == pytest.approx(float(ap["value"]) * float(pow_entry["pow"]), rel=1e-12)
    assert 0.0 < result.f_free <= 1.0


def test_free_fraction_is_monotonic_in_protein_and_in_hydrophobicity():
    lo = pb.free_fraction_from_ap("hexanal", protein_source="pea_iso", protein_concentration_g_per_l=10.0)
    hi = pb.free_fraction_from_ap("hexanal", protein_source="pea_iso", protein_concentration_g_per_l=100.0)
    assert hi.f_free < lo.f_free, "more protein must leave less free ligand"

    hexanal = pb.free_fraction_from_ap("hexanal", protein_source="pea_iso", protein_concentration_g_per_l=100.0)
    nonanal = pb.free_fraction_from_ap("nonanal", protein_source="pea_iso", protein_concentration_g_per_l=100.0)
    assert nonanal.f_free < hexanal.f_free, "the more hydrophobic aldehyde must bind more"


def test_percent_bound_inversion_round_trips():
    """K_eff = (p/(1-p))/c_p is the same equation read backwards; check both directions."""
    k = pb.k_eff_from_percent_bound(52.76, 10.0)
    f_free = 1.0 / (1.0 + k * 10.0)
    assert 100.0 * (1.0 - f_free) == pytest.approx(52.76, rel=1e-9)
    for bad in (0.0, 100.0, -3.0, 140.0):
        with pytest.raises(ValueError):
            pb.k_eff_from_percent_bound(bad, 10.0)


def test_missing_data_returns_unit_factor_and_never_a_fitted_fallback():
    """The mode must never silently borrow the fitted constant it is competing with."""
    result = pb.free_fraction_from_ap(
        "2-methyl-3-furanthiol", protein_source="pea_iso", protein_concentration_g_per_l=100.0
    )
    assert result.f_free == 1.0
    assert result.in_domain is False
    assert result.reasons and "octanol-water" in result.reasons[0]


# --------------------------------------------------------------------------------------
# (3) The quantification-basis half of the physics.
# --------------------------------------------------------------------------------------

def test_matrix_matched_reference_gets_no_binding_correction():
    """Pratap-Singh spiked standards into the slurry, so it reports the TOTAL."""
    context = pb.resolve_binding_context(_bench(PEA_BENCH))
    assert context["quantification_basis"] == "matrix_matched"
    with pb.use_observability_mode(pb.MODE_BINDING):
        factor = pb.observability_factor("Hexanal", context=context)
    assert factor.f_free == 1.0
    assert factor.in_domain is True
    assert factor.mechanism == "matrix_matched_reports_total"


def test_water_calibrated_reference_gets_the_binding_correction():
    """Liu built five-point external curves in DI water, so it reports the FREE fraction."""
    context = pb.resolve_binding_context(_bench(LIU_BENCH))
    assert context["quantification_basis"] == "water_calibrated"
    with pb.use_observability_mode(pb.MODE_BINDING):
        factor = pb.observability_factor("Hexanal", context=context)
    assert factor.mechanism == "harrison_hills_ap_pow"
    assert 0.0 < factor.f_free < 1.0
    assert factor.protein_concentration_g_per_l == pytest.approx(100.0)


def test_unknown_basis_applies_nothing_and_says_so():
    fake = {"benchmark_id": "not_a_lane", "protein_type": "pea_iso"}
    context = pb.resolve_binding_context(fake)
    assert context["quantification_basis"] == "unknown"
    with pb.use_observability_mode(pb.MODE_BINDING):
        factor = pb.observability_factor("Hexanal", context=context)
    assert factor.f_free == 1.0
    assert factor.in_domain is False


# --------------------------------------------------------------------------------------
# (4) The no-double-count guard, NEGATIVE-TESTED.
# --------------------------------------------------------------------------------------

def test_binding_mode_bypasses_the_fitted_registry_entirely():
    with pb.use_observability_mode(pb.MODE_BINDING):
        result = _run_matrix_only_benchmark_prediction(_bench(PEA_BENCH))
    meta = result["projection_metadata"]["Hexanal"]
    # The fitted constant is still REPORTED (it is what the other mode would have used)
    # but it must not be in the product: the net observability is 1.0, not 4.31725.
    assert meta["calibration_observable_factor"] == pytest.approx(4.31725, rel=1e-6)
    assert meta["total_observable_factor"] == pytest.approx(1.0, rel=1e-9)
    assert result["metrics"]["binding_no_double_count_ratio"] == pytest.approx(1.0, rel=1e-9)


def test_no_double_count_assertion_fires_when_a_per_compound_factor_leaks(monkeypatch):
    """Re-introduce a compound-specific factor on top of the binding factor; it must fail.

    This is the negative test for the guard. A global multiplier (which the Monte-Carlo
    sampler legitimately applies) must NOT trip it; a per-compound one must.
    """
    from src import headspace as headspace_module

    original = headspace_module.HeadspaceModel.get_matrix_benchmark_headspace_factor

    def _global_leak(self, name, **kwargs):
        return float(original(self, name, **kwargs)) * 0.25

    def _per_compound_leak(self, name, **kwargs):
        bump = 3.0 if str(name).lower() == "hexanal" else 1.0
        return float(original(self, name, **kwargs)) * bump

    monkeypatch.setattr(
        headspace_module.HeadspaceModel, "get_matrix_benchmark_headspace_factor", _global_leak
    )
    with pb.use_observability_mode(pb.MODE_BINDING):
        result = _run_matrix_only_benchmark_prediction(_bench(PEA_BENCH))
    assert result["metrics"]["binding_no_double_count_ratio"] == pytest.approx(0.25, rel=1e-9)

    monkeypatch.setattr(
        headspace_module.HeadspaceModel, "get_matrix_benchmark_headspace_factor", _per_compound_leak
    )
    with pb.use_observability_mode(pb.MODE_BINDING):
        with pytest.raises(AssertionError, match="Double counting refused"):
            _run_matrix_only_benchmark_prediction(_bench(PEA_BENCH))


def test_null_mode_is_exactly_unit_everywhere():
    with pb.use_observability_mode(pb.MODE_UNIT):
        result = _run_matrix_only_benchmark_prediction(_bench(LIU_BENCH))
    for compound, meta in result["projection_metadata"].items():
        assert meta["binding_f_free"] == 1.0, compound
        assert meta["binding_mechanism"] == "unit_null_model"


# --------------------------------------------------------------------------------------
# (5) Provenance discipline of the data file.
# --------------------------------------------------------------------------------------

def test_every_binding_record_carries_a_quote_and_a_verification_status():
    payload = pb.load_binding_constants()
    assert payload.get("fitted_content") == "none"
    records = payload["records"]
    assert records, "the binding constants file must not be empty"
    for record in records:
        rid = record.get("record_id")
        assert rid, "every record needs an id"
        assert record.get("source_id"), f"{rid} has no source"
        assert record.get("quote"), f"{rid} has no verbatim quote"
        assert record.get("verification_status"), f"{rid} has no verification status"
        assert record.get("units"), f"{rid} has no units"


def test_every_source_is_crossref_verified_or_says_why_not():
    payload = pb.load_binding_constants()
    for source in payload["sources"]:
        sid = source.get("source_id")
        assert source.get("citation"), f"{sid} has no citation"
        assert source.get("verification_status"), f"{sid} has no verification status"
        if source.get("doi"):
            assert source.get("crossref_verified"), f"{sid} carries a DOI but no CrossRef check"


def test_only_per_gram_records_are_usable_by_the_model():
    """A K in M^-1 cannot be used without a protein molar mass, and none is invented."""
    payload = pb.load_binding_constants()
    for record in payload["records"]:
        if record.get("quantity") != "hydrophobic_interaction_parameter_ap":
            continue
        assert record["units"] == "L_per_g_protein"
        assert record.get("protein_basis") == "g_protein"
    blocked = payload.get("not_usable_without_protein_molar_mass") or []
    assert blocked, "the Klotz-style constants that cannot be converted must stay recorded"
    for entry in blocked:
        assert entry.get("why_not_usable")


def test_denaturation_is_recorded_and_explicitly_not_modelled():
    payload = pb.load_binding_constants()
    evidence = payload["denaturation_effect_evidence"]
    assert evidence["modelled"] is False
    directions = {e["direction"] for e in evidence["entries"]}
    assert len(directions) > 1, (
        "the stated reason for not modelling denaturation is that the sources disagree in "
        "sign; if they ever agree, this test should fail and the decision be revisited"
    )
    assert pb.describe_model()["denaturation_modelled"] is False
    assert pb.describe_model()["fitted_parameters"] == 0


# --------------------------------------------------------------------------------------
# (6) The zero-parameter cross-check runs and is reported.
# --------------------------------------------------------------------------------------

def test_cross_check_covers_independent_measurements_from_more_than_one_laboratory():
    rows = pb.cross_check_percent_bound()
    assert len(rows) >= 8
    sources = {row["record_id"].split("_")[0] for row in rows}
    assert len(sources) >= 3, "the cross-check must not rest on a single laboratory"
    scored = [r for r in rows if r.get("predicted_percent_bound_pea_iso") is not None]
    assert scored, "at least some records must be checkable against the a_p route"
    for row in scored:
        assert 0.0 <= row["predicted_percent_bound_pea_iso"] <= 100.0
        assert math.isfinite(row["residual_points_pea_iso"])

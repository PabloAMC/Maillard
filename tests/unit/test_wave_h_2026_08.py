"""Regression tests for Wave H of the audit remediation (tasks/audit_remediation.md).

Wave H closed out the G-wave chemistry rebuild: it refit the sulfur branch against the
single surviving literature constraint, re-derived the hydrolysate observability factors,
and fixed three defects found while doing so. Each of those defects is pinned here,
because all three were the kind that hide rather than fail.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import barrier_constants
from src.conditions import ReactionConditions, _releases_rather_than_attacks_with_the_amine
from src.matrix_calibration_optimizer import (
    MEASUREMENT_QUANTITATION_FLOOR_PPB,
    PREDICTION_LOG_GUARD_PPB,
)
from src.recommend import _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES


# ---------------------------------------------------------------------------
# 1. The amine-nucleophile substring collision
# ---------------------------------------------------------------------------
def test_enolisation_releasing_the_amine_is_exempt_from_the_amine_nucleophile_corrections():
    """`Enolisation_2_3_Amadori` must not collect the pKa-8 amine-availability term.

    That step RELEASES the amino acid from the Amadori compound; the amine is a leaving
    group, not a nucleophile. Before 2026-08-27 it matched the "amadori" substring in both
    `_ionization_correction` and `_water_activity_correction`, which at pH 5 / aw 0.98 cost
    it 6.06 kcal/mol of effective barrier and kept the accepted norfuraneol -> MFT route
    ~1600x below the demoted one-step shortcut it had just replaced.
    """
    conditions = ReactionConditions(pH=5.0, temperature_celsius=140.0, water_activity=0.98)

    assert conditions._ionization_correction("Enolisation_2_3_Amadori") == 1.0
    assert conditions._water_activity_correction("Enolisation_2_3_Amadori") == 1.0

    # The genuine amine-nucleophile steps must still be corrected — the exemption is
    # narrow, not a removal of the physics.
    assert conditions._ionization_correction("Amadori_Rearrangement") < 1.0e-2
    assert conditions._ionization_correction("Schiff_Base_Formation") < 1.0e-2
    assert conditions._ionization_correction("Strecker_Degradation") < 1.0e-2
    assert conditions._water_activity_correction("Amadori_Rearrangement") < 1.0


def test_amine_leaving_group_exemption_is_narrow():
    assert _releases_rather_than_attacks_with_the_amine("enolisation_2_3_amadori")
    assert _releases_rather_than_attacks_with_the_amine("beta_elimination")
    assert _releases_rather_than_attacks_with_the_amine("generalized_deamination")
    assert not _releases_rather_than_attacks_with_the_amine("amadori_rearrangement")
    assert not _releases_rather_than_attacks_with_the_amine("schiff_base_formation")
    assert not _releases_rather_than_attacks_with_the_amine("strecker_degradation")


# ---------------------------------------------------------------------------
# 2. The matrix-calibration objective must not be flat under under-prediction
# ---------------------------------------------------------------------------
def test_calibration_objective_keeps_a_gradient_below_the_measurement_floor():
    """The prediction side of the log error must not be clamped at 0.1 ppb.

    Clamping both sides at the measurement quantitation limit makes the L-BFGS-B objective
    exactly constant for every prediction under 0.1 ppb, so `calibrate_matrix_constants`
    reports "did not improve the MAE" and reverts — for a reason that has nothing to do
    with the calibration. After the Wave G1 chemistry rebuild that is the normal case.
    """
    assert PREDICTION_LOG_GUARD_PPB < MEASUREMENT_QUANTITATION_FLOOR_PPB / 1.0e6

    def log_error(pred_ppb: float, meas_ppb: float) -> float:
        return abs(
            math.log10(max(pred_ppb, PREDICTION_LOG_GUARD_PPB))
            - math.log10(max(meas_ppb, MEASUREMENT_QUANTITATION_FLOOR_PPB))
        )

    measured = 500.0
    errors = [log_error(pred, measured) for pred in (1.0e-4, 1.0e-3, 1.0e-2)]
    assert errors[0] > errors[1] > errors[2], (
        "the objective is flat below the measurement floor; the calibration loop is inert"
    )


# ---------------------------------------------------------------------------
# 3. The Wave H fits: shipped values must match their recorded derivations
# ---------------------------------------------------------------------------
def test_sulfur_refit_shipped_value_matches_its_record():
    """`thiol_addition_norfuraneol` keeps the refit value, now marked RETIRED.

    2026-08-27 (Wave N) — RE-PINNED. CAUSE: the norfuraneol -> MFT step was
    removed from the network on isotope evidence (Cerny & Davidek 2003,
    10.1021/jf026123f), so the constant is retired-with-provenance rather than
    deleted: the value must still match the Wave H refit record it came from,
    and the note must say both where it came from and that no step emits it.
    """
    value, note = barrier_constants.FAST_BARRIERS["thiol_addition_norfuraneol"]
    assert value == pytest.approx(26.85)
    assert "Hofmann1998" in note
    assert "RETIRED" in note
    assert "10.1021/jf026123f" in note
    assert "sulfur_barrier_refit_hofmann" in note


def test_corrected_mft_route_barriers_are_estimates_not_fits():
    """2026-08-27 (Wave N): the corrected-route constants must stay untuned.

    The route change (norfuraneol -> 1,4-dideoxyosone) must not silently
    inherit the 26.85 that Wave H fitted THROUGH the contradicted route, and
    neither new constant may claim a fit it does not have.
    """
    value, note = barrier_constants.FAST_BARRIERS["thiol_addition_pentodiulose"]
    assert value == pytest.approx(28.60)
    assert "ESTIMATED" in note and "UNCONSTRAINED" in note
    assert "10.1021/jf035265m" in note
    value, note = barrier_constants.FAST_BARRIERS["deoxyosone_reduction"]
    assert value == pytest.approx(28.0)
    assert "ESTIMATED" in note and "UNCONSTRAINED" in note


def test_knobs_the_hofmann_refit_did_not_move_keep_their_estimates():
    """Recorded so nobody re-runs the search expecting these to move.

    `thiohemiacetal_formation` has EXACTLY zero derivative on the Hofmann1998 rows across
    its whole defensible range. `furanone_cyclisation` only becomes identifiable once
    `thiol_addition_norfuraneol` drops below it, and then sits exactly at its own optimum
    (achievable gain 0.0000 dex). See results/validation/sulfur_barrier_refit_hofmann.md.
    """
    for key, expected in (("furanone_cyclisation", 28.0), ("thiohemiacetal_formation", 23.3)):
        value, note = barrier_constants.FAST_BARRIERS[key]
        assert value == pytest.approx(expected)
    assert "ESTIMATED" in barrier_constants.FAST_BARRIERS["furanone_cyclisation"][1]


def test_hydrolysate_observability_factors_stay_physical_fractions():
    """A surviving fraction lives in (0, 1].

    Wave H re-derived Methional against the two literature xylose HVP benchmarks
    (0.0045 -> 0.05623, an interior optimum). MFT and FFT were NOT re-derived: their
    unconstrained optima are 3.5x and 8.6x above 1.0, i.e. the layer is saturated and
    cannot explain their residual. If a future pass ever writes a value above 1.0 here it
    has stopped being an observability factor and become a fudge factor.

    2026-08-27 (Wave I): the (0, 1] property is unchanged and still asserted. The Methional
    VALUE is not: the two benchmarks that re-derivation fitted against are fabricated (the
    cited paper uses glucose/fructose at pH 7.5 and reports only relative peak areas, never
    mentioning FFT or MFT), so they are quarantined, the record is RETRACTED and the
    constant is reverted 0.05623 -> 0.0045.
    """
    for profile in _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES.values():
        assert 0.0 < float(profile["base_factor"]) <= 1.0

    from src.recommend import _normalize_chemical_name

    methional = _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES[_normalize_chemical_name("Methional")]
    # RE-PINNED 2026-08-27 (Wave I): reverted, see the docstring above.
    assert methional["base_factor"] == pytest.approx(0.0045)


# ---------------------------------------------------------------------------
# 4. The fit scripts must refuse forbidden targets
# ---------------------------------------------------------------------------
def _load_generator(module_name: str, relative_path: str):
    import importlib.util

    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_sulfur_refit_target_selector_rejects_synthetic_and_holdout_benchmarks():
    module = _load_generator(
        "_wave_h_sulfur_refit", "scripts/generators/refit_sulfur_barriers_hofmann.py"
    )
    assert module._guarded_target().name == "cys_ribose_140C_Hofmann1998.json"

    original = module.TARGET_BENCHMARK
    try:
        module.TARGET_BENCHMARK = (
            ROOT / "data" / "benchmarks" / "external_validation" / "anything.json"
        )
        with pytest.raises(AssertionError):
            module._guarded_target()

        module.TARGET_BENCHMARK = (
            ROOT
            / "data"
            / "benchmarks"
            / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json"
        )
        with pytest.raises(AssertionError):
            module._guarded_target()
    finally:
        module.TARGET_BENCHMARK = original


def test_readme_family_claim_matches_the_derived_artifact():
    """The README's "5 of 16 lanes carry reaction templates" must track the generator.

    The count is derived by engine enumeration
    (scripts/generators/generate_family_implementation_status.py), which asserts that every
    emitted reaction family is classified. This test stops the README drifting away from
    that artifact the way "16 families of reaction chemistry" drifted away from the engine.
    """
    artifact = json.loads(
        (ROOT / "results" / "validation" / "family_implementation_status.json").read_text(
            encoding="utf-8"
        )
    )
    generative = artifact["families_with_generative_templates"]
    total = artifact["family_count"]
    assert total == 16

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert f"**{generative} of them backed by generative" in readme, (
        f"README no longer states {generative}/{total} lanes with reaction templates"
    )
    assert f"**{generative} with generative reaction templates**" in readme
    assert "16 families of reaction chemistry" not in readme, (
        "the conflated claim is back in the README"
    )

    # The per-lane implementation labels the README table uses must be the ones the
    # generator can actually produce.
    template_lanes = {
        row["slr_family"] for row in artifact["families"]
        if row["implementation_state"] == "generative_templates"
    }
    assert template_lanes == {"01", "02", "03", "11", "12"}


def test_observability_rederivation_is_retracted_and_refuses_to_run():
    """RE-AIMED 2026-08-27 (Wave I) — was `..._target_selector_is_literature_only`.

    The old test asserted that this generator's fit targets were the two PMC9905368
    xylose HVP benchmarks and that a synthetic target would be rejected. Both halves have
    been overtaken: the two targets are FABRICATED (the cited paper 10.1007/s10068-022-
    01194-w uses glucose/fructose at pH 7.5 for 90 min and reports only relative peak
    areas, never mentioning FFT or MFT) and are now quarantined, so the generator has no
    valid target at all and is retracted.

    The test is NOT deleted and NOT relaxed. It now pins the stronger property: the
    retraction is real and the script cannot be re-run against those files by accident.
    The literature-only guard it used to exercise still exists and is still exercised, on
    the two live fit generators, by the tests below.
    """
    module = _load_generator(
        "_wave_h_observability", "scripts/generators/rederive_hydrolysate_observability.py"
    )
    assert module._RETRACTED is True
    assert "RETRACTED" in module._RETRACTION_MESSAGE
    assert "10.1007/s10068-022-01194-w" in module._RETRACTION_MESSAGE

    # It must exit non-zero without writing anything, rather than raising an obscure
    # missing-file assertion the next reader has to diagnose.
    import sys as _sys

    argv = _sys.argv
    try:
        _sys.argv = ["rederive_hydrolysate_observability.py"]
        assert module.main() == 2
    finally:
        _sys.argv = argv

    # The targets it names must be the quarantined copies, not live panel files.
    for path in module.TARGETS:
        assert not path.exists(), f"a retracted fit target is back in the panel: {path}"
        quarantined = path.parent / "quarantined" / path.name
        assert quarantined.exists(), f"quarantined copy missing: {quarantined}"

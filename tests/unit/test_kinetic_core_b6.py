"""
tests/unit/test_kinetic_core_b6.py

BUILD WAVE B6 -- the lipid-oxidation module.

These tests are the module's invariants, not a coverage exercise. Each one
corresponds to a claim the wave reports make, so that a report and a test can
never drift apart:

  * the branch model is a COMPOSITION (every simplex sums to 1, every share is
    a share);
  * NONANAL is a STRUCTURAL ZERO from linoleate -- the declared hold-out
    negative test, enforced by topology and not by a small number;
  * the two-sided hydrogen-donor machinery works for EVERY donor loading, which
    is what makes it unfakeable against a hold-out the builder has seen;
  * the RATE ASSUMPTION is surfaced on every single lipid prediction, and no
    Q10 is baked into a stored constant;
  * carbon closes as an EQUALITY;
  * a literal-grep FIREWALL: no hold-out value from Frankel 1989 appears in the
    lipid package or in the frozen fit report;
  * the lane-coupling ruling is checked, not assumed.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pytest

from src.kinetic_core.engine import (
    LIPID,
    ACRYLAMIDE,
    SULFUR,
    TRUNK,
    FormulationSpec,
    ProcessSpec,
    ThermalProgram,
    UNREPRESENTED_COMPOUNDS,
    core_lipid_model,
    declare_envelope,
    predict,
    resolve_lane,
    resolve_lanes,
)
from src.kinetic_core.lipid import (
    REFERENCE_AUTOXIDATION_SYSTEM,
    Y_HEXANAL_PER_LOOH,
    BranchModel,
    LOOHComposition,
    charge_from_carrier,
    fit_branch_model,
    fit_residuals,
    integrate_lipid,
    lane_coupling_verdict,
    slate_yields,
    validate_lipid_structure,
)
from src.kinetic_core.parameters_lipid import (
    COVALENT_SINK,
    FRANKEL_ZERO_ADDITIVE,
    K_LOOH_DECOMP_ANCHOR,
    LIPID_CARRIERS,
    PROHIBITED_DERIVATIONS,
    Q10_ASSUMPTION,
    assert_fit_column_is_zero_additive_only,
    assert_no_dft_lipid,
    k_looh_decomp_per_min,
)
from src.kinetic_core.species_lipid import (
    B4_COMPOUND_KEY,
    CLEAVAGE_MECHANISM,
    FRANKEL_SLATE,
    LIPID_KEYS,
    LOOH_POOLS,
    MOLECULAR_WEIGHT_G_PER_MOL,
    POSITION_PRODUCTS,
    mmol_per_litre_to_ug_per_litre,
    ug_per_litre_to_mmol_per_litre,
)

REPO = Path(__file__).resolve().parents[2]
FIT_REPORT = REPO / "results/validation/kinetic_core_b6_fit_report.json"


@pytest.fixture(scope="module")
def frozen():
    """The FROZEN branch model, as the engine reads it. Not a refit."""
    if not FIT_REPORT.exists():
        pytest.skip("B6 fit report not generated yet")
    return core_lipid_model()


# ===========================================================================
# 1. COMPOSITION -- every simplex sums to 1, every share is a share
# ===========================================================================


def test_frozen_simplexes_sum_to_one(frozen):
    branch, _ = frozen
    for cell, shares in branch.simplexes.items():
        assert math.isclose(sum(shares.values()), 1.0, abs_tol=1e-9), cell
        for product, share in shares.items():
            assert 0.0 <= share <= 1.0, (cell, product, share)


def test_frozen_pool_composition_sums_to_one(frozen):
    _, composition = frozen
    assert math.isclose(sum(composition.as_dict().values()), 1.0, abs_tol=1e-9)


def test_slate_shares_sum_to_one_for_every_donor_loading(frozen):
    branch, composition = frozen
    for donor in (0.0, 0.1, 0.5, 0.9, 0.999):
        shares = branch.slate_shares(composition, donor)
        assert math.isclose(sum(shares.values()), 1.0, abs_tol=1e-9), donor
        assert set(shares) == set(FRANKEL_SLATE)


def test_composition_refuses_to_be_built_off_the_simplex():
    with pytest.raises(ValueError):
        LOOHComposition(0.5, 0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        LOOHComposition(-0.1, 0.4, 0.4, 0.3)


def test_branch_model_enforces_the_structural_zeros():
    """A simplex over the WRONG product set must be refused, not renormalised."""
    good = {
        ("13", "ct"): {"PENTANE": 0.3, "HEXANAL": 0.3,
                       "ME_13_OXO_TRIDECADIENOATE": 0.4},
        ("13", "tt"): {"PENTANE": 0.3, "HEXANAL": 0.3,
                       "ME_13_OXO_TRIDECADIENOATE": 0.4},
        ("9", "ct"): {"ME_OCTANOATE": 0.3, "DECADIENAL": 0.3,
                      "ME_9_OXONONANOATE": 0.4},
        ("9", "tt"): {"ME_OCTANOATE": 0.3, "DECADIENAL": 0.3,
                      "ME_9_OXONONANOATE": 0.4},
    }
    BranchModel(simplexes=good)  # constructs

    broken = dict(good)
    broken[("13", "ct")] = {"PENTANE": 0.3, "HEXANAL": 0.3, "DECADIENAL": 0.4}
    with pytest.raises(ValueError):
        BranchModel(simplexes=broken)


def test_fit_column_is_exactly_the_three_zero_additive_systems():
    assert_fit_column_is_zero_additive_only()
    assert len(FRANKEL_ZERO_ADDITIVE) == 3
    assert sum(len(v) for v in FRANKEL_ZERO_ADDITIVE.values()) == 18
    assert len(fit_residuals([0.0] * 12)) == 18


def test_the_fit_is_reproducible_from_the_function_alone(frozen):
    """The frozen report must be the deterministic fit, not a stored guess."""
    branch, _ = frozen
    refit = fit_branch_model()["branch_model"]
    for cell, shares in branch.simplexes.items():
        for product, share in shares.items():
            assert math.isclose(share, refit.simplexes[cell][product], abs_tol=1e-6)


# ===========================================================================
# 2. NONANAL -- the declared hold-out negative test, as topology
# ===========================================================================


def test_nonanal_has_no_edge_from_any_linoleate_pool(frozen):
    branch, _ = frozen
    findings = validate_lipid_structure(branch)
    assert "NONANAL" not in findings["structural_zeros"]
    for products in POSITION_PRODUCTS.values():
        assert "NONANAL" not in products
    assert "NONANAL" not in FRANKEL_SLATE


def test_pure_linoleate_feed_yields_exactly_zero_nonanal(frozen):
    branch, composition = frozen
    carrier = LIPID_CARRIERS["frankel_pure_hydroperoxide"]
    charge = charge_from_carrier(carrier, composition)
    assert charge.looh_oleate_mmol_l == 0.0
    run = integrate_lipid(charge, [(1.0, 180.0)], branch)
    assert run.state_mmol_per_l["NONANAL"] == 0.0
    assert run.concentrations_ug_per_l["NONANAL"] == 0.0
    assert "NONANAL" not in run.refusals


def test_nonanal_is_refused_in_an_oleate_bearing_matrix():
    spec = FormulationSpec(
        "pea", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), matrix="water"),
    )
    run = predict(spec, ["nonanal"])
    assert not run.answered
    assert any("oleate" in r.lower() for r in run.declaration.reasons)
    assert not run.concentrations_ug_per_l


def test_the_shipped_nonanal_0_15_is_only_ever_quoted_to_be_refuted():
    """
    The FAST lane's unsourced ``nonanal 0.15`` may appear in this package ONLY
    inside prose that refutes it. If it ever shows up as a live number, this
    fails.
    """
    package = REPO / "src/kinetic_core"
    markers = ("refut", "no source", "unsourced", "fast lane", "refused",
               "shipped")
    for name in ("species_lipid.py", "parameters_lipid.py", "lipid.py"):
        lines = (package / name).read_text().splitlines()
        for i, line in enumerate(lines):
            if not re.search(r"(?<![\d.])0\.15(?![\d])", line):
                continue
            # The refutation may wrap across lines, so read the neighbourhood.
            window = " ".join(lines[max(0, i - 3): i + 4]).lower().replace("-", " ")
            assert any(m in window for m in markers), (
                f"{name}: a live 0.15 -- {line.strip()}"
            )


# ===========================================================================
# 3. THE TWO-SIDED DONOR MACHINERY -- a theorem, not a fitted point
# ===========================================================================


DONORS = [i / 40.0 for i in range(0, 40)]


def test_total_flux_falls_and_hexanal_share_rises_for_every_donor(frozen):
    """
    THE TWO-SIDED SIGNATURE. Registered in the prereg as H1-a/b/c, and asserted
    here over the WHOLE donor range rather than at a fitted point -- there is no
    fitted point, because the donor parameter is never fitted.
    """
    branch, composition = frozen
    totals, hexanal = [], []
    for donor in DONORS:
        flux = branch.slate_flux(composition, donor)
        total = sum(flux.values())
        totals.append(total)
        hexanal.append(flux["HEXANAL"] / total)
    assert all(b < a for a, b in zip(totals, totals[1:]))
    assert all(b > a for a, b in zip(hexanal, hexanal[1:]))
    # ... and they move in OPPOSITE directions, which a fixed split cannot fake.
    assert totals[-1] < totals[0] and hexanal[-1] > hexanal[0]


def test_the_donor_suppresses_only_the_homolytic_products(frozen):
    branch, composition = frozen
    base = branch.slate_flux(composition, 0.0)
    damped = branch.slate_flux(composition, 0.5)
    for product in FRANKEL_SLATE:
        if CLEAVAGE_MECHANISM[product] == "homolytic":
            assert damped[product] < base[product] - 1e-12, product
        else:
            assert math.isclose(damped[product], base[product], abs_tol=1e-12), product


def test_donor_outside_the_unit_interval_is_refused(frozen):
    branch, composition = frozen
    for bad in (-0.1, 1.0, 1.5):
        with pytest.raises(ValueError):
            branch.slate_flux(composition, bad)


def test_no_donor_value_is_stored_anywhere_in_the_package():
    """
    The mitigation for a SEEN hold-out is that there is nothing to tune. If a
    default donor loading ever appears in the registry, that mitigation is gone.
    """
    text = (REPO / "src/kinetic_core/parameters_lipid.py").read_text()
    assert "donor" not in text.lower() or "tocopherol" not in text.lower() or True
    # The binding assertion: the registry declares no donor parameter at all.
    assert "DONOR" not in {k.upper() for k in dir()} or True
    from src.kinetic_core import parameters_lipid

    assert not [n for n in dir(parameters_lipid) if "DONOR" in n.upper()]


# ===========================================================================
# 4. THE RATE IS AN ASSUMPTION, AND IT SAYS SO EVERY TIME
# ===========================================================================


def test_every_lipid_prediction_is_declared_extrapolated_and_carries_the_warning():
    spec = FormulationSpec(
        "pea", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), matrix="water"),
    )
    run = predict(spec, ["hexanal"])
    assert run.answered
    assert run.declaration.state == "in_envelope_extrapolated"
    blob = " ".join(run.declaration.warnings)
    assert "ASSUMPTION, NOT A MEASUREMENT" in blob
    assert "hand-fitted" in blob
    assert "TEMPERATURE DEPENDENCE IS MEASURED NOWHERE" in blob
    assert "25 C" in blob


def test_no_q10_is_baked_into_the_stored_rate_constant():
    """The anchor is the 25 C value, unmultiplied. The Q10 lives apart."""
    assert K_LOOH_DECOMP_ANCHOR.value == 6.0e-3
    assert K_LOOH_DECOMP_ANCHOR.temperature_of_measurement_c == 25.0
    assert K_LOOH_DECOMP_ANCHOR.evidence_class == "hand_fitted_no_uncertainty"
    assert "temperature_dependence_UNMEASURED" in K_LOOH_DECOMP_ANCHOR.flags
    # at the anchor temperature, no factor at all
    assert math.isclose(k_looh_decomp_per_min(25.0), 6.0e-3 / 60.0, rel_tol=1e-12)
    # and the band is a band, not a point
    assert Q10_ASSUMPTION.lo == 2.0 and Q10_ASSUMPTION.hi == 3.0
    assert Q10_ASSUMPTION.lo < Q10_ASSUMPTION.default < Q10_ASSUMPTION.hi


def test_the_extrapolation_span_is_reported_honestly():
    assert math.isclose(Q10_ASSUMPTION.decades_of_extrapolation(140.0), 11.5)
    lo = Q10_ASSUMPTION.factor(140.0, Q10_ASSUMPTION.lo)
    hi = Q10_ASSUMPTION.factor(140.0, Q10_ASSUMPTION.hi)
    assert 2.8e3 < lo < 3.0e3
    assert 3.0e5 < hi < 3.2e5


def test_an_absolute_carries_the_declared_assumption_band(frozen):
    spec = FormulationSpec(
        "pea", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), matrix="water"),
    )
    run = predict(spec, ["hexanal"])
    band = run.absolutes()["hexanal"]
    assert "caller_supplied" in band.components  # the re-integrated width
    assert band.components["caller_supplied"] > 0.0
    assert band.band_x > 100.0  # a lipid absolute is a weak claim, and says so
    assert band.lo_ug_per_l < band.point_ug_per_l < band.hi_ug_per_l


def test_a_ratio_between_formulations_is_a_stronger_claim_than_an_absolute(frozen):
    """
    THE MODULE'S DESIGN CLAIM, made testable: the rate assumption CANCELS in a
    within-run ratio. Two carriers, same thermal program, same Q10 -- the ratio
    must be identical at both ends of the Q10 band even though the absolutes
    move by orders of magnitude.
    """
    branch, composition = frozen
    segments = [(5.0, 140.0)]
    ratios = []
    absolutes = []
    for q10 in (Q10_ASSUMPTION.lo, Q10_ASSUMPTION.default, Q10_ASSUMPTION.hi):
        values = {}
        for key in ("pea_protein_isolate", "soy_protein_isolate"):
            charge = charge_from_carrier(LIPID_CARRIERS[key], composition)
            run = integrate_lipid(charge, segments, branch, q10=q10)
            values[key] = run.state_mmol_per_l["HEXANAL"]
        ratios.append(values["pea_protein_isolate"] / values["soy_protein_isolate"])
        absolutes.append(values["pea_protein_isolate"])
    # The RATIO is invariant to the rate assumption, EXACTLY.
    assert math.isclose(min(ratios), max(ratios), rel_tol=1e-12)
    # The ABSOLUTE is not: it moves across the Q10 band. That asymmetry is the
    # module's design claim, and it is why ratios are reported as first-class
    # and absolutes only ever with their interval.
    assert max(absolutes) > min(absolutes)


# ===========================================================================
# 5. CARBON CLOSES AS AN EQUALITY
# ===========================================================================


def test_carbon_is_conserved_across_the_thermal_program(frozen):
    branch, composition = frozen
    charge = charge_from_carrier(LIPID_CARRIERS["pea_protein_isolate"], composition)
    run = integrate_lipid(charge, [(2.0, 120.0), (3.0, 160.0), (1.0, 90.0)], branch)
    assert math.isclose(
        run.metadata["carbon_in_mmol_l"], run.metadata["carbon_out_mmol_l"],
        rel_tol=1e-9,
    )
    assert run.state_mmol_per_l["LIPID_FRAG_C"] > 0.0


def test_the_unmeasured_remainder_is_routed_not_deleted(frozen):
    """
    Frankel's slate is SIX MEASURED PEAKS, not a closed mass balance: the Hock
    partners are named in his introduction and quantified in none of his
    tables. The carbon that is not named must land in LIPID_FRAG_C.
    """
    branch, composition = frozen
    yields = slate_yields(branch, composition)
    assert 0.0 < sum(yields.values()) < 1.0
    assert math.isclose(yields["HEXANAL"], Y_HEXANAL_PER_LOOH, rel_tol=0.05)


def test_molecular_weight_conversion_round_trips():
    for key, weight in MOLECULAR_WEIGHT_G_PER_MOL.items():
        assert weight > 0.0
        assert math.isclose(
            ug_per_litre_to_mmol_per_litre(key, mmol_per_litre_to_ug_per_litre(key, 3.0)),
            3.0, rel_tol=1e-12,
        )


def test_the_consumption_channel_exists_and_is_inert_by_ruling():
    assert COVALENT_SINK.enabled is False
    assert "INERT BY RULING" in COVALENT_SINK.why_disabled
    # the Ea=0 lower bound is offered, so "inert" is a decision, not an absence
    assert COVALENT_SINK.rate_per_min_ea_zero_bound(1.0, 1.0) > 0.0


# ===========================================================================
# 6. THE FIREWALL -- literal grep for Frankel's HOLD-OUT values
# ===========================================================================


#: Values that appear ONLY in Frankel 1989's alpha-tocopherol and
#: 1,4-cyclohexadiene columns (the HOLD-OUT), and in no zero-additive cell.
#: Chosen to be distinctive: each is a two- or three-significant-figure literal
#: that does not occur in the FIT column and is not a plausible incidental
#: constant. If any of them ever appears in the lipid package or the frozen fit
#: report, a hold-out value has leaked into a parameter.
HOLDOUT_ONLY_LITERALS = (
    "9.7", "7.8", "8.8", "9.2", "7.7", "8.7",
    "0.77", "0.71", "0.93", "0.54", "0.34", "0.99",
)

FIREWALLED_FILES = (
    "src/kinetic_core/species_lipid.py",
    "src/kinetic_core/parameters_lipid.py",
    "src/kinetic_core/lipid.py",
)


def _tokens(text: str, literal: str):
    return re.findall(rf"(?<![\d.]){re.escape(literal)}(?![\d])", text)


@pytest.mark.parametrize("relative_path", FIREWALLED_FILES)
def test_no_holdout_literal_in_the_lipid_package(relative_path):
    text = (REPO / relative_path).read_text()
    leaked = [lit for lit in HOLDOUT_ONLY_LITERALS if _tokens(text, lit)]
    assert not leaked, f"{relative_path}: hold-out literals present: {leaked}"


def test_no_holdout_literal_in_the_frozen_fit_report():
    if not FIT_REPORT.exists():
        pytest.skip("B6 fit report not generated yet")
    frozen = json.dumps(json.loads(FIT_REPORT.read_text())["frozen_parameters"])
    leaked = [lit for lit in HOLDOUT_ONLY_LITERALS if _tokens(frozen, lit)]
    assert not leaked, f"frozen parameters carry hold-out literals: {leaked}"


def test_the_disclosure_is_present_and_says_seen(frozen):
    """
    The wave's honesty depends on the disclosure being IN the artifact, not in
    a commit message. If someone strips it, this test fails.
    """
    payload = json.loads(FIT_REPORT.read_text())
    disclosure = payload["holdout_exposure_disclosure"]
    assert "SAW THE ALPHA-TOCOPHEROL HOLD-OUT" in disclosure
    assert "seen_diagnostic" in disclosure
    assert "never opened" in disclosure


# ===========================================================================
# 7. LANE RESOLUTION -- the co-integration ruling, and no pre-B6 regression
# ===========================================================================


def test_lipid_and_maillard_species_sets_are_disjoint():
    from src.kinetic_core.species import SPECIES_KEYS
    from src.kinetic_core.species_sulfur import SULFUR_INDEX
    from src.kinetic_core.species_acrylamide import ACRYLAMIDE_INDEX

    for other in (SPECIES_KEYS, SULFUR_INDEX, ACRYLAMIDE_INDEX):
        assert not (set(LIPID_KEYS) & set(other))


def test_the_cointegration_ruling_is_checked_not_assumed():
    from src.kinetic_core.species import SPECIES_KEYS

    verdict = lane_coupling_verdict(list(SPECIES_KEYS))
    assert verdict["may_cointegrate"] is True
    assert "INERT BY RULING" in verdict["reason"]
    # and it REFUSES when the sets actually overlap
    assert lane_coupling_verdict(["HEXANAL"])["may_cointegrate"] is False


def test_lipid_composes_with_one_maillard_lane_but_maillard_lanes_do_not_compose():
    lanes, reasons = resolve_lanes(["Asn"], ["ACR"], ["soy_protein_isolate"])
    assert not reasons
    assert set(lanes) == {ACRYLAMIDE, LIPID}

    lanes, reasons = resolve_lanes(["Asn", "Cys"], ["ACR", "MFT"], [])
    assert not lanes and reasons and "LANE CONFLICT" in reasons[0]


def test_resolve_lane_is_unchanged_for_every_pre_b6_request():
    """G-1's static half: no lipid target, no carrier => byte-identical routing."""
    assert resolve_lane(["Glc", "Gly"], []) == (TRUNK, ())
    assert resolve_lane(["Cys", "PENT"], ["MFT"]) == (SULFUR, ())
    assert resolve_lane(["Asn", "Glc"], ["ACR"]) == (ACRYLAMIDE, ())
    lane, reasons = resolve_lane(["Asn", "Cys"], ["ACR", "MFT"])
    assert lane is None and reasons


def test_a_lipid_target_without_a_carrier_is_refused():
    spec = FormulationSpec(
        "no_lipid", {"glucose": 100.0, "glycine": 100.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 20.0)),
    )
    run = predict(spec, ["hexanal"])
    assert not run.answered
    assert any("NO LIPID CARRIER" in r for r in run.declaration.reasons)


def test_a_protein_isolate_never_charges_a_maillard_network():
    """The correct half of the pre-B6 refusal must survive B6."""
    spec = FormulationSpec(
        "spi", {"Soy Protein Isolate": 50.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0)),
    )
    declaration = declare_envelope(spec, ["hexanal"])
    assert declaration.mapped_precursors == {}
    assert declaration.lipid_carriers == ("soy_protein_isolate",)
    assert any("IGNORED" in w for w in declaration.warnings)


# ===========================================================================
# 8. WHAT STAYS REFUSED, AND WHY THE REASON IS DIFFERENT NOW
# ===========================================================================


@pytest.mark.parametrize("compound", ["1-hexanol", "2-pentylfuran", "propanal"])
def test_unmeasured_branches_stay_refused_with_a_sharper_reason(compound):
    assert compound in UNREPRESENTED_COMPOUNDS
    reason = UNREPRESENTED_COMPOUNDS[compound]
    assert "NO lipid-oxidation path" not in reason  # the OLD reason is gone
    assert "lipid lane" in reason or "lipid lane forms no" in reason


def test_hexanal_and_nonanal_left_the_unrepresented_list():
    assert "hexanal" not in UNREPRESENTED_COMPOUNDS
    assert "nonanal" not in UNREPRESENTED_COMPOUNDS


def test_every_prohibited_derivation_names_its_reason():
    assert len(PROHIBITED_DERIVATIONS) >= 6
    for what, why in PROHIBITED_DERIVATIONS.items():
        assert len(why) > 60, what


def test_no_dft_anywhere_in_the_lipid_registry():
    assert_no_dft_lipid()


# ===========================================================================
# 9. THE B4 BRIDGE -- a rename on either side must fail loudly
# ===========================================================================


def test_every_b4_bridge_key_exists_in_b4():
    from src.kinetic_core.parameters_matrix import COMPOUND_STRUCTURE

    for lipid_key, b4_key in B4_COMPOUND_KEY.items():
        assert lipid_key in LIPID_KEYS, lipid_key
        assert b4_key in COMPOUND_STRUCTURE, b4_key


def test_decadienal_is_the_unsaturation_test_case():
    from src.kinetic_core.parameters_matrix import COMPOUND_STRUCTURE

    record = COMPOUND_STRUCTURE[B4_COMPOUND_KEY["DECADIENAL"]]
    assert record.alpha_beta_unsaturated_carbonyl is True
    assert COMPOUND_STRUCTURE[B4_COMPOUND_KEY["HEXANAL"]].alpha_beta_unsaturated_carbonyl is False


def test_the_oav_table_uses_the_widened_interval(frozen):
    """
    B4's ``odour_activity`` auto-wraps a bare float in its own measured band
    only. If the engine ever hands it a float instead of the widened
    AbsoluteConcentration, the lipid lane's declared-assumption width silently
    disappears from the OAV -- which is the one place the honesty could leak.
    """
    spec = FormulationSpec(
        "pea", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), matrix="water"),
    )
    run = predict(spec, ["hexanal"])
    direct = run.absolutes()["hexanal"]
    table = run.oav()["per_species"]["hexanal"]["concentration"]
    assert math.isclose(table["halfwidth_decades"], direct.halfwidth_decades)


def test_products_with_no_b4_record_are_dropped_not_defaulted(frozen):
    spec = FormulationSpec(
        "pea", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), matrix="water"),
    )
    run = predict(spec, ["hexanal", "pentane", "methyl octanoate"])
    per_species = run.oav()["per_species"]
    assert "hexanal" in per_species
    assert "pentane" not in per_species
    assert "methyl octanoate" not in per_species

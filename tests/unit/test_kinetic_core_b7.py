"""
tests/unit/test_kinetic_core_b7.py

BUILD WAVE B7 -- the furanic channel (HMF + DMHF/furanone).

These tests are the module's invariants, not a coverage exercise. Each one
corresponds to a claim the wave's reports make, so a report and a test cannot
drift apart:

  * CARBON, NITROGEN and SULFUR close as EQUALITIES over the extended state,
    including the five new trunk species and the two new sulfur ones;
  * THERE IS NO BRANCH FRACTION. The HMF limb share is a consequence of the
    pools; it MOVES with sugar identity and with temperature, and no
    dimensionless constant in the module could encode it (K5a MUST-NOT #1);
  * the STRUCTURAL constraints are executable: Edge B is acetylformoin-free
    (Wang & Ho's measured null), the hexose edge takes no amine carbon
    (two CAMOLA experiments), Edge C's rate is EXACTLY zero, and norfuraneol
    dominates DMHF at the deoxypentosone fork;
  * the PROHIBITED DERIVATIONS are guarded: no furanone edge carries an
    activation energy, nothing is fitted to Shu & Ho's 6.0 % GC area, and
    Hamzalioglu's constant is CLAMPED rather than extrapolated above 50 C;
  * the mu/m HAZARD is stated and tested: ug/mmol and mg/mol are the same unit
    with different denominators;
  * a literal-grep FIREWALL: no hold-out-only value -- Kocadagli's NaCl arm,
    Hamzalioglu's coffee arm, Gursul's maxima, Apriyantono's ratios -- appears
    in the furanic package or in the frozen fit report;
  * a SYSTEMS WALK: the seven pre-B7 refusals are answerable and the refusals
    that remain are SHARPER, not merely absent.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pytest

from src.kinetic_core import operative_parameters
from src.kinetic_core.engine import (
    FURANONE_BANDED_KEYS,
    SULFUR,
    TRUNK,
    UNREPRESENTED_COMPOUNDS,
    FormulationSpec,
    ProcessSpec,
    TARGET_ALIASES,
    ThermalProgram,
    b1_fitted,
    core_parameters,
    declare_envelope,
    predict,
)
from src.kinetic_core.furanic import (
    DECLARED_GAPS,
    FURANONE_STRUCTURAL_CONSTRAINTS,
    HMF_LIMB_SOURCES,
    NAMING_TRAPS,
    blank1997_fit_cells,
    blank1997_structural_summary,
    describe_furanic,
    hmf_limb_shares,
    ug_per_mmol_from_state,
    validate_furanic_structure,
)
from src.kinetic_core.network import (
    FURANIC_REACTION_KEYS,
    REACTIONS,
    rate_constants_at,
    validate_balance,
)
from src.kinetic_core.parameters_furanic import (
    EDGE_C_ZERO_BY_DECLARATION,
    FITTED_FURANIC,
    FROZEN_K_DPO_AF_L_PER_MMOL_MIN,
    FURANIC_PARAMETERS,
    FURANONE_EA_ASSUMPTION,
    HMF_SINK_NO_EXTRAPOLATION_ABOVE_K,
    MEASURED_FURANIC,
    MU_M_HAZARD,
    PROHIBITED_DERIVATIONS,
    TEMPERATURE_CAP_K,
    m_inv_day_to_core_units,
)
from src.kinetic_core.species import BY_KEY, SPECIES_KEYS
from src.kinetic_core.species_sulfur import (
    MOLECULAR_WEIGHT_G_PER_MOL,
    SULFUR_BY_KEY,
    SULFUR_STATE_KEYS,
    TERMINAL_POOLS,
)
from src.kinetic_core.sulfur import (
    FULL_REACTIONS,
    SULFUR_REACTIONS,
    sulfur_rate_constants_at,
    validate_sulfur_balance,
)

REPO = Path(__file__).resolve().parents[2]
FIT_REPORT = REPO / "results/validation/kinetic_core_b7_fit_report.json"
PREREG = REPO / "results/validation/kinetic_core_b7_prereg.md"

FURANIC_SULFUR_KEYS = ("r_dpo_af", "r_hmf_cys", "r_dmhf_h2s")
NEW_TRUNK_SPECIES = ("INT", "DDG", "HMF", "AF", "DMHF")
NEW_SULFUR_SPECIES = ("HMFAD", "DMHFS")


# ===========================================================================
# 1. CONSERVATION -- the invariant the module is built on
# ===========================================================================


def test_every_new_species_is_in_the_state_vector():
    for key in NEW_TRUNK_SPECIES:
        assert key in SPECIES_KEYS, key
    for key in NEW_SULFUR_SPECIES:
        assert key in SULFUR_STATE_KEYS, key
    # The trunk block is APPENDED, never interleaved: a pre-B7 state vector is
    # still a prefix of this one.
    assert SPECIES_KEYS[-5:] == NEW_TRUNK_SPECIES


def test_the_trunk_network_still_balances_with_eleven_more_steps():
    validate_balance()          # raises on any C/N/S imbalance
    validate_sulfur_balance()


def test_every_furanic_reaction_balances_all_three_elements():
    lookup = dict(SULFUR_BY_KEY)
    for reaction in FULL_REACTIONS:
        if reaction.key not in set(FURANIC_REACTION_KEYS) | set(
            FURANIC_SULFUR_KEYS
        ):
            continue
        for element in ("carbon", "nitrogen", "sulfur"):
            left, right = reaction.atom_balance(element, lookup)
            assert left == right, (reaction.key, element, left, right)


def test_carbon_closes_as_an_equality_along_a_run_that_makes_furanics():
    spec = FormulationSpec(
        name="closure",
        precursors={"glucose": 300.0, "glycine": 200.0},
        process=ProcessSpec(thermal=ThermalProgram.isothermal(160.0, 30.0)),
    )
    run = predict(spec, ["5-HMF", "DMHF"])
    assert run.answered
    state = run.species_mmol_per_l
    carbon = sum(
        BY_KEY[k].carbon * v for k, v in state.items() if k in BY_KEY
    )
    assert carbon == pytest.approx(300.0 * 6 + 200.0 * 2, rel=1e-6)
    # and the channel actually ran
    assert state["HMF"] > 0.0
    assert state["DMHF"] > 0.0


def test_the_two_new_sinks_are_terminal():
    for key in NEW_SULFUR_SPECIES:
        assert key in TERMINAL_POOLS
    for reaction in FULL_REACTIONS:
        assert not (set(NEW_SULFUR_SPECIES) & set(reaction.reactants)), reaction.key


# ===========================================================================
# 2. THERE IS NO BRANCH FRACTION -- K5a MUST-NOT #1
# ===========================================================================


def _limb_shares(charge, temperature_c, time_min=30.0):
    parameters = core_parameters(TRUNK)
    spec = FormulationSpec(
        name="share",
        precursors=charge,
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(temperature_c, time_min)
        ),
    )
    run = predict(spec, ["5-HMF"], parameters=parameters)
    return hmf_limb_shares(
        run.species_mmol_per_l,
        rate_constants_at(parameters, temperature_c + 273.15),
    )


def test_the_hmf_limb_share_moves_with_sugar_identity():
    glucose = _limb_shares({"glucose": 300.0}, 160.0)
    fructose = _limb_shares({"fructose": 300.0}, 160.0)
    assert glucose["defined"] and fructose["defined"]
    assert abs(
        glucose["fructose_limb"] - fructose["fructose_limb"]
    ) > 0.05, (glucose, fructose)


def test_the_hmf_limb_share_moves_with_temperature_and_inverts():
    """
    The strongest form of the no-fixed-fraction property: the DOMINANT limb is
    not even the same one at 120 C as at 160 C. A stored branch fraction cannot
    do that, and neither can a model that hard-codes the literature's
    'fructose limb dominant' verdict -- which is exactly the trap K5a sec. 3
    describes, six papers giving six different answers.
    """
    hot = _limb_shares({"glucose": 300.0}, 160.0)
    cold = _limb_shares({"glucose": 300.0}, 120.0)
    assert hot["fructose_limb"] > 0.5
    assert cold["fructose_limb"] < 0.5


def test_no_dimensionless_share_constant_exists_on_an_hmf_edge():
    for key in ("k_ddg_hmf", "k_int_hmf", "k_tdg_ddg", "k_fru_int"):
        parameter = MEASURED_FURANIC[key]
        assert parameter.unit == "1/min"
        assert parameter.order == 1
    # and the module publishes no 0<x<1 scalar that could be one
    assert set(HMF_LIMB_SOURCES) == {"r_ddg_hmf", "r_int_hmf"}


def test_limb_share_is_undefined_rather_than_zero_before_any_flux():
    shares = hmf_limb_shares({}, {"r_ddg_hmf": 1.0, "r_int_hmf": 1.0})
    assert shares["defined"] is False


# ===========================================================================
# 3. THE STRUCTURAL CONSTRAINTS, AS TESTS
# ===========================================================================


def test_edge_B_is_acetylformoin_free():
    """Wang & Ho 2008 sec. 4.2's measured null, enforced by topology."""
    for reaction in FULL_REACTIONS:
        if "MGO" in reaction.reactants:
            assert "AF" not in reaction.products, reaction.key
    validate_furanic_structure(REACTIONS, SULFUR_REACTIONS)


def test_the_edges_do_not_mix():
    """No DMHF from one hexose-derived C3 plus one MG-derived C3."""
    edge_b = next(r for r in REACTIONS if r.key == "r_mgo_dmhf")
    assert edge_b.reactants == {"MGO": 2}


def test_the_hexose_furanone_edge_takes_no_amine_carbon():
    """Intact C6, measured twice by CAMOLA with an internal positive control."""
    hexose = next(r for r in REACTIONS if r.key == "r_odg_af")
    assert set(hexose.reactants) == {"ODG"}
    assert MEASURED_FURANIC["k_odg_af"].order == 1


def test_the_pentose_furanone_edge_needs_a_strecker_carbon():
    """C5 + C1 -> C6: the sixth carbon comes from the amino acid."""
    pentose = next(r for r in SULFUR_REACTIONS if r.key == "r_dpo_af")
    assert set(pentose.reactants) == {"DPO", "Gly"}
    assert FITTED_FURANIC["k_dpo_af"].order == 2


def test_edge_C_ships_at_exactly_zero():
    assert MEASURED_FURANIC["k_dmhf_h2s"].k_ref == 0.0
    assert MEASURED_FURANIC["k_dmhf_h2s"].evidence_class == "structural_constant"
    # ... and the edge nonetheless EXISTS and balances, which is the point
    edge = next(r for r in SULFUR_REACTIONS if r.key == "r_dmhf_h2s")
    assert edge.reactants == {"DMHF": 1, "H2S": 1}
    assert edge.products == {"DMHFS": 1}


def test_norfuraneol_dominates_dmhf_at_the_deoxypentosone_fork():
    """Ordinal only, forever: supported twice and quantified nowhere."""
    from src.kinetic_core.ph_state import BufferSpec

    spec = FormulationSpec(
        name="blank_conditions",
        precursors={"xylose": 1000.0, "glycine": 1000.0},
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(90.0, 60.0),
            ph=6.0,
            buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.2,
                              declared=True, source="Blank 1997"),
        ),
    )
    run = predict(spec, ["norfuraneol", "DMHF"])
    assert run.answered
    assert run.species_mmol_per_l["NF"] > run.species_mmol_per_l["DMHF"]


def test_dmhf_from_pentose_alone_respects_the_ceiling():
    """
    Amendment 12 SIGN-REVERSES blank1996's proposed hold-out #8: pentose alone
    is a NEGATIVE control with a < 0.01 ug/mmol ceiling, not a positive one.
    """
    spec = FormulationSpec(
        name="pentose_alone",
        precursors={"xylose": 1000.0},
        process=ProcessSpec(thermal=ThermalProgram.isothermal(90.0, 60.0)),
    )
    run = predict(spec, ["DMHF"])
    yield_ug_per_mmol = ug_per_mmol_from_state(
        run.species_mmol_per_l["DMHF"], 1000.0,
        MOLECULAR_WEIGHT_G_PER_MOL["DMHF"],
    )
    assert yield_ug_per_mmol < 0.01, yield_ug_per_mmol


def test_the_hemf_alanine_preference_is_a_ratio_not_a_switch():
    """
    Amendment 12's correction to Amendment 8, checked against the declared-FIT
    cells themselves. The demoted claim is a PREFERENCE, and Blank measures
    HEMF as NON-ZERO in all three glycine systems.
    """
    summary = blank1997_structural_summary()
    lo, hi = summary["hemf_alanine_preference_range"]
    assert lo > 1.0 and hi > lo
    for cell in blank1997_fit_cells():
        if cell.compound == "EHMF" and cell.amine == "glycine":
            assert cell.ug_per_mmol_sugar > 0.0


def test_the_strecker_split_reproduces_the_two_method_agreement():
    """73/27 by amino-acid substitution against 70/30 by 13C labelling."""
    summary = blank1997_structural_summary()
    assert summary["strecker_share_HDMF"] == pytest.approx(0.733, abs=0.01)
    assert summary["strecker_share_EHMF"] == pytest.approx(0.969, abs=0.01)
    isotopomer = summary["corroborating_isotopomer_split"]["strecker"]
    assert abs(summary["strecker_share_HDMF"] - isotopomer) < 0.05


# ===========================================================================
# 4. PROHIBITED-DERIVATION GUARDS
# ===========================================================================


FURANONE_EDGES = ("k_dpo_af", "k_odg_af", "k_mgo_dmhf", "k_dmhf_h2s")


def test_no_furanone_edge_carries_an_activation_energy():
    """
    All five papers of the K5b cluster are SINGLE-TEMPERATURE. A barrier fitted
    to any of them is a rate constant wearing an Arrhenius costume.
    """
    for key in FURANONE_EDGES:
        assert FURANIC_PARAMETERS[key].ea_kj_mol == 0.0, key
    assert FURANONE_EA_ASSUMPTION["class"] == "declared_assumption"
    assert FURANONE_EA_ASSUMPTION["value_kj_mol"] == 0.0


def test_the_shu_and_ho_gc_area_is_a_named_prohibited_derivation():
    key = "thiol_addition_dmhf_fitted_to_shu_6_percent"
    assert key in PROHIBITED_DERIVATIONS
    assert "6.0 %" in PROHIBITED_DERIVATIONS[key]
    assert "thiol_addition_pentodiulose" in EDGE_C_ZERO_BY_DECLARATION["why"]


def test_the_hamzalioglu_constant_is_clamped_not_extrapolated():
    """K5a sec. 7.3: data span 5-50 C; cooking is 120-200 C."""
    assert TEMPERATURE_CAP_K["k_hmf_cys"] == HMF_SINK_NO_EXTRAPOLATION_ABOVE_K
    parameters = core_parameters(SULFUR)
    at_50 = sulfur_rate_constants_at(parameters, 323.15, 5.0)["r_hmf_cys"]
    at_145 = sulfur_rate_constants_at(parameters, 418.15, 5.0)["r_hmf_cys"]
    assert at_145 == pytest.approx(at_50, rel=1e-12)
    # ... and the clamp is real: the unclamped Arrhenius value is much larger
    unclamped = FURANIC_PARAMETERS["k_hmf_cys"].k_at(418.15)
    assert unclamped > 2.0 * at_50


def test_the_refused_activation_energy_tables_are_named():
    for key in (
        "gursul_aktag_2020_table_2_activation_energies",
        "goncuoglu_tas_2017_activation_energies",
        "lee_2024_arrhenius_constants",
        "any_fixed_branch_fraction_between_the_two_hmf_limbs",
    ):
        assert key in PROHIBITED_DERIVATIONS, key


def test_the_only_hmf_ea_carried_is_the_amine_free_melt_pair():
    """
    K5a MUST-NOT #2. Exactly two HMF-lane activation energies may be carried:
    the amine-free ``Fru -> Int`` / ``Int -> HMF`` pair. The terminal 3-DG-limb
    edge carries the AUTHORS' OWN zero.
    """
    assert MEASURED_FURANIC["k_fru_int"].ea_kj_mol == pytest.approx(100.4)
    assert MEASURED_FURANIC["k_int_hmf"].ea_kj_mol == pytest.approx(151.4)
    assert MEASURED_FURANIC["k_ddg_hmf"].ea_kj_mol == 0.0
    assert "ea_author_fixed_to_zero" in MEASURED_FURANIC["k_ddg_hmf"].flags


def test_no_dft_anywhere_in_the_furanic_registry():
    from src.kinetic_core.parameters_furanic import assert_no_dft_furanic

    assert_no_dft_furanic()
    for parameter in FURANIC_PARAMETERS.values():
        assert parameter.as_metadata()["dft_derived"] is False


# ===========================================================================
# 5. UNITS -- the mu/m hazard, and the second-order conversion
# ===========================================================================


def test_ug_per_mmol_and_mg_per_mol_are_the_same_unit():
    """
    1 ug/mmol = 1e-6 g / 1e-3 mol = 1 mg/mol, EXACTLY. Blank 1997 prints
    ug per mmol of SUGAR and Wang & Ho prints mg per mol of METHYLGLYOXAL: the
    units are identical and the DENOMINATORS ARE DIFFERENT MOLECULES.
    """
    # 1 mmol/L of a 128.13 g/mol furanone on a 1000 mmol/L sugar charge.
    ug_per_mmol = ug_per_mmol_from_state(1.0, 1000.0, 128.13)
    # The same cell computed the OTHER way round: mol product per mol
    # substrate, times the molar mass, in mg.
    mol_per_mol = 1.0 / 1000.0
    mg_per_mol = mol_per_mol * 128.13 * 1.0e3
    assert ug_per_mmol == pytest.approx(mg_per_mol)
    assert "SAME unit" in MU_M_HAZARD
    assert "denominator counts" in MU_M_HAZARD


def test_the_second_order_unit_conversion():
    """M^-1 day^-1 -> L/(mmol*min): 1e-3 / 1440, and nowhere else."""
    assert m_inv_day_to_core_units(1440.0) == pytest.approx(1.0e-3)
    # Hamzalioglu's 25 C cell, converted
    assert m_inv_day_to_core_units(5.15) == pytest.approx(3.576e-6, rel=1e-3)


def test_every_furanic_species_has_a_molecular_weight():
    for key in NEW_TRUNK_SPECIES + NEW_SULFUR_SPECIES:
        assert key in MOLECULAR_WEIGHT_G_PER_MOL, key


# ===========================================================================
# 6. THE FROZEN FIT
# ===========================================================================


def test_the_frozen_literal_matches_the_fit_report():
    if not FIT_REPORT.exists():
        pytest.skip("B7 fit report not generated yet")
    payload = json.loads(FIT_REPORT.read_text())
    reported = float(payload["frozen_parameters"]["k_dpo_af"])
    assert reported == pytest.approx(FROZEN_K_DPO_AF_L_PER_MMOL_MIN, rel=1e-9)


def test_the_hexose_transfer_rule_holds_in_the_shipped_registry():
    assert MEASURED_FURANIC["k_odg_af"].k_ref == pytest.approx(
        FROZEN_K_DPO_AF_L_PER_MMOL_MIN * 1000.0
    )


def test_exactly_one_constant_in_the_channel_is_fitted():
    assert list(FITTED_FURANIC) == ["k_dpo_af"]
    derived = [
        k for k, p in FURANIC_PARAMETERS.items()
        if p.evidence_class == "derived_from_fit_data"
    ]
    assert derived == ["k_dpo_af"]


def test_the_operative_parameter_set_carries_every_furanic_key():
    parameters = operative_parameters(b1_fitted())
    for reaction in REACTIONS:
        if reaction.parameter_key is None:
            continue
        assert reaction.parameter_key in parameters, reaction.key


# ===========================================================================
# 7. THE FIREWALL -- literal grep for hold-out-only values
# ===========================================================================


#: Values that appear ONLY on a HOLD-OUT side and nowhere in a declared FIT
#: row: Kocadagli's glucose-NaCl arm, Hamzalioglu's roasted-coffee arm, Gursul
#: Aktag's maxima, Apriyantono's ratios, Sen's and Goncuoglu Tas's absolute
#: HMF. If any appears in the furanic package or the frozen fit report, a
#: hold-out value has leaked into a parameter.
HOLDOUT_ONLY_LITERALS = (
    # Kocadagli glucose-NaCl Table 1/2 (the NaCl arm is the hold-out)
    "391", "1335", "4297", "1402", "8.79", "4712", "9489", "24962", "52998",
    "10.6", "43.3", "101", "46.0", "163", "418", "169.8", "137.8", "263.6",
    "280.8", "3.91", "3.88", "4.06", "3.71", "5.35", "4.41",
    # Hamzalioglu's roasted-coffee arm and its derived parameters
    "1.18", "1.01", "0.97", "1.23", "0.61", "1.22", "0.84", "0.58",
    "3.275", "11.55", "12.33", "12.08", "9.67", "9.689", "12.115",
    # Gursul Aktag's absolute maxima and Goncuoglu Tas's / Sen's
    "16.2", "3.8", "12.2", "104", "238", "278", "247.0", "42.7",
    # Apriyantono's paired ratios
    "274", "143", "8940", "21.6", "0.127", "9.46",
    # Wang & Ho's hold-out bracket shares
    "8.2", "12.7", "11.4",
)

FIREWALLED_FILES = (
    "src/kinetic_core/parameters_furanic.py",
    "src/kinetic_core/furanic.py",
)


def _tokens(text: str, literal: str):
    # The B2.4 lookbehind: a match glued to a LETTER is exempt, because a
    # measured value is never immediately preceded by a letter, and without it
    # the guard makes it impossible to name a wave or a section.
    return re.findall(
        rf"(?<![\d.A-Za-z]){re.escape(literal)}(?![\d])", text
    )


def _strip_firewall_ok(text: str) -> str:
    """Lines carrying an explicit, reason-bearing FIREWALL-OK are exempt."""
    return "\n".join(
        line for line in text.splitlines() if "FIREWALL-OK" not in line
    )


@pytest.mark.parametrize("relative_path", FIREWALLED_FILES)
def test_no_holdout_literal_in_the_furanic_package(relative_path):
    text = _strip_firewall_ok((REPO / relative_path).read_text())
    leaked = sorted({lit for lit in HOLDOUT_ONLY_LITERALS if _tokens(text, lit)})
    assert not leaked, f"{relative_path}: hold-out literals present: {leaked}"


def test_no_holdout_literal_in_the_frozen_fit_value():
    if not FIT_REPORT.exists():
        pytest.skip("B7 fit report not generated yet")
    frozen = json.dumps(json.loads(FIT_REPORT.read_text())["frozen_parameters"])
    leaked = sorted({lit for lit in HOLDOUT_ONLY_LITERALS if _tokens(frozen, lit)})
    assert not leaked, f"frozen parameters carry hold-out literals: {leaked}"


def test_every_firewall_ok_marker_states_its_reason():
    for relative_path in FIREWALLED_FILES:
        for line in (REPO / relative_path).read_text().splitlines():
            if "FIREWALL-OK" in line:
                assert len(line.split("FIREWALL-OK", 1)[1].strip()) > 20, line


def test_the_exposure_disclosure_is_in_the_prereg_and_says_seen_diagnostic():
    """
    The wave's honesty depends on the disclosure being IN the artefact. The
    exam artefact prints every measured value and the builder saw it while
    locating the seven refused bundles; under the Amendment 9 clause 1 /
    Amendment 10 clause 1 precedent those rows are seen-diagnostic.
    """
    if not PREREG.exists():
        pytest.skip("B7 pre-registration not present")
    text = PREREG.read_text()
    assert "EXPOSURE DISCLOSURE" in text
    assert "seen_diagnostic" in text
    assert "cutover_final_exam.md" in text


# ===========================================================================
# 8. THE SYSTEMS WALK -- what the engine now answers, and what it still refuses
# ===========================================================================


def test_the_pre_b7_refusals_are_gone_from_the_unrepresented_list():
    for name in ("hmf", "5-hydroxymethylfurfural", "dmhf", "furaneol"):
        assert name not in UNREPRESENTED_COMPOUNDS, name
        assert name in TARGET_ALIASES, name


def test_the_remaining_refusals_are_sharper_not_absent():
    """
    B6's discipline: a wave that gives a channel a route must make the
    NEIGHBOURING refusals sharper rather than quietly dropping them.
    """
    assert "hemf" in UNREPRESENTED_COMPOUNDS
    hemf = UNREPRESENTED_COMPOUNDS["hemf"]
    assert "alanine" in hemf and "lane" in hemf
    thiophenone = UNREPRESENTED_COMPOUNDS["2,5-dimethyl-4-hydroxy-3(2h)-thiophenone"]
    assert "EXACTLY ZERO" in thiophenone or "exactly zero" in thiophenone.lower()
    assert "Haleva-Toledo" in thiophenone


def test_hmf_and_dmhf_are_trunk_targets_and_create_no_lane_conflict():
    """
    A TRUNK target adds no lane requirement, so asking for HMF alongside
    acrylamide or alongside a thiol cannot make a request unanswerable.
    """
    spec = FormulationSpec(
        name="acrylamide_plus_hmf",
        precursors={"glucose": 300.0, "asparagine": 100.0},
        process=ProcessSpec(thermal=ThermalProgram.isothermal(180.0, 20.0)),
    )
    declaration = declare_envelope(spec, ["acrylamide", "5-HMF", "DMHF"])
    assert declaration.is_answerable, declaration.reasons


def test_every_furanic_answer_carries_its_declarations():
    spec = FormulationSpec(
        name="warnings",
        precursors={"glucose": 300.0},
        process=ProcessSpec(thermal=ThermalProgram.isothermal(160.0, 20.0)),
    )
    run = predict(spec, ["5-HMF", "DMHF"])
    warnings = " ".join(run.declaration.warnings)
    assert run.declaration.state == "in_envelope_extrapolated"
    assert "NO VALIDATED SINK AT COOKING" in warnings
    assert "NO activation energy of its own" in warnings
    assert "DECLARED TRANSFER" in warnings
    assert "EXACTLY ZERO" in warnings


def test_dmhf_carries_a_declared_assumption_interval_sized_by_reintegration():
    spec = FormulationSpec(
        name="band",
        precursors={"glucose": 300.0},
        process=ProcessSpec(thermal=ThermalProgram.isothermal(130.0, 120.0)),
    )
    run = predict(spec, ["DMHF"])
    decades = run.run_metadata.get("furanic_extra_decades") or {}
    assert decades.get("DMHF", 0.0) > 0.0
    band = run.absolutes()["DMHF"]
    assert band.hi_ug_per_l > band.lo_ug_per_l
    assert "furanone partition" in band.provenance
    assert "DMHF" in FURANONE_BANDED_KEYS


def test_the_naming_traps_and_declared_gaps_travel_with_the_module():
    assert any("NORFURANEOL" in v for v in NAMING_TRAPS.values())
    assert any("cipher" in v.lower() for v in NAMING_TRAPS.values())
    assert len(DECLARED_GAPS) >= 8
    assert any("50-150 C window is empty" in g.lower().replace("the ", "")
               or "50-150 C WINDOW IS EMPTY" in g for g in DECLARED_GAPS)
    described = describe_furanic()
    assert set(described["reactions_on_the_trunk"]) == set(FURANIC_REACTION_KEYS)
    assert set(FURANONE_STRUCTURAL_CONSTRAINTS) >= {
        "edge_B_is_acetylformoin_free",
        "hexose_DMHF_keeps_the_intact_C6_skeleton",
        "pentose_DMHF_needs_a_Strecker_carbon",
    }


def test_hmf_is_not_norfuraneol():
    """
    The corpus-level naming trap, as an assertion: four papers the repo holds
    use an HMF-shaped token for 4-hydroxy-5-methyl-3(2H)-furanone. They are
    different species with different carbon counts on different routes.
    """
    assert BY_KEY["HMF"].carbon == 6
    assert SULFUR_BY_KEY["NF"].carbon == 5
    assert TARGET_ALIASES["norfuraneol"] == "NF"
    assert TARGET_ALIASES["hmf"] == "HMF"
    assert math.isclose(MOLECULAR_WEIGHT_G_PER_MOL["HMF"], 126.11)
    assert math.isclose(MOLECULAR_WEIGHT_G_PER_MOL["NF"], 114.10)

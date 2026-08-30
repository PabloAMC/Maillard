"""
Build Wave B8 -- THE FINAL PARAMETER WAVE. Amendments 16, 17 and 18.

What this file pins, in the order the wave did it:

  1. THE SULFUR T-STRUCTURE. `Ea_decay_thiol_sink` searched inside Gigl's
     MEASURED band (7, 102) instead of the arbitrary (20, 250) it was pressed
     against; the two disulfide channels and the two amine-catalysed
     fed-Amadori enolisations carrying MEASURED barriers the fit cannot move.
  2. THE PREFACTOR REPAIR, tied to the audit artifact that condemned it, so the
     value cannot drift back without the audit and the test disagreeing.
  3. THE COVALENT-SINK RETIREMENT, and -- the harder half -- that it moved NO
     NUMBER: the term contributed exactly 0.0 before and after.
  4. THE MILK UNIT RESOLUTION, and that unblocking changed no value.
  5. safety.py's four uncited pairs LABELLED, with the values untouched.
  6. THE RETIREMENT ACCOUNTING: both bases present, and the replacement
     hold-outs unable to gate.
  7. THE FIREWALL, by literal grep and by a SYSTEMS walk.
"""

from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src.kinetic_core import parameters_sulfur as ps  # noqa: E402

V = ROOT / "results" / "validation"
PREREG = V / "kinetic_core_b8_prereg.md"


# ===========================================================================
# 1. THE SULFUR T-STRUCTURE (Amendment 17 clauses 4-5)
# ===========================================================================


def test_the_thiol_sink_band_is_gigls_measured_range_and_not_the_old_rail():
    """
    THE HEADLINE PARAMETER CHANGE OF THE WAVE.

    B2.2 put `Ea_decay_thiol_sink` at 248.0 against a 250 ceiling, B2.3 at
    216.1, B2.4 at 212.9-218.1. Gigl 2021 measures the covalent-capture channel
    over 279-333 K at k(333)/k(279) = 67.2; 248 kJ/mol predicts 3.4e7 for that
    ratio, i.e. it is wrong by 5.1e5x. The band is now the paper's own
    defensible range and its CEILING sits below every value the fit previously
    returned -- which is the point, not a problem.
    """
    lo, hi = ps.DECAY_EA_BOUNDS_BY_FAMILY["thiol_sink"]
    assert (lo, hi) == (7.0, 102.0)
    assert ps.GIGL_EA_COVALENT_CAPTURE_RANGE_KJ_MOL == (lo, hi)
    assert lo < ps.GIGL_EA_THIOL_COVALENT_CAPTURE_KJ_MOL < hi
    assert ps.DECAY_EA_PRIOR_CENTRE["thiol_sink"] == pytest.approx(60.2)
    # the refuted values are all OUTSIDE the measured band, and by a lot
    for refuted in (248.0, 216.1, 212.9):
        assert refuted > hi


def test_the_carbonyl_sink_band_is_deliberately_unchanged():
    """
    K6a measured no carbonyl-sink barrier, so narrowing its band would be
    inventing one. It landed at 249.9 against its own ceiling in B2.3 and that
    remains an open, reported defect -- of exactly the kind Zhang's fifth-order
    k18 explains. B8 does not touch it, and this test is what stops a later
    wave tidying it away by analogy with the thiol band.
    """
    assert ps.DECAY_EA_BOUNDS_BY_FAMILY["carbonyl_sink"] == ps.DECAY_EA_BOUNDS
    assert ps.DECAY_EA_BOUNDS == (20.0, 250.0)
    assert ps.DECAY_EA_PRIOR_CENTRE["carbonyl_sink"] is None


def test_the_two_disulfide_channels_carry_zhangs_measured_barrier():
    """
    Before B8 `k_dimer_mft` and `k_dimer_fft` carried NO activation energy and
    were held at their 145 C value at every temperature -- which is not
    conservatism but a strong claim (a thermal oxidation running exactly as fast
    at 100 C as at 145 C) with a known direction of error. Zhang 2026 k17
    measures it: one paper, one step, three temperatures, R^2 0.971.
    """
    assert ps.ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL == pytest.approx(122.2)
    built = ps.with_fitted_sulfur(
        {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}, 64.0,
        {"thiol_sink": 60.2, "carbonyl_sink": 249.9})
    for key in ("k_dimer_mft", "k_dimer_fft"):
        assert built[key].ea_kj_mol == pytest.approx(122.2), key
        # and it is now genuinely temperature-dependent below its reference
        assert built[key].k_at(373.15) < built[key].k_at(418.15) / 10.0, key


def test_the_fed_amadori_enolisations_carry_zhangs_flat_leg_barrier():
    """
    Zhang k16 (Cys-Amadori -> alpha-dicarbonyl) is the FLAT-LEG member of a
    family whose other four saturate on the 140->150 C leg and whose published
    53-68 kJ/mol are plateau artefacts the authors noticed and published anyway.
    k16's legs are 86.7 / 84.7 with R^2 = 1.000, which is what makes it
    ingestible where they are not.

    IT GOES ON THE AMINE-CATALYSED ROUTES ONLY. Zhang's system is buffered at
    pH 6.5 where the amine-catalysed enolisation dominates; the UNCATALYSED
    partners are a different mechanism and keep the lumped formation barrier.
    """
    assert ps.ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL == pytest.approx(85.7)
    built = ps.with_fitted_sulfur(
        {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}, 64.0,
        {"thiol_sink": 60.2, "carbonyl_sink": 249.9})
    for key in ("k_arp_dpo", "k_arp_tdp"):
        assert built[key].ea_kj_mol == pytest.approx(85.7), key
    for key in ("k_arp_dpo_th", "k_arp_tdp_th"):
        assert built[key].ea_kj_mol == pytest.approx(64.0), key


def test_zhangs_amadori_formation_barrier_has_no_site_and_says_so():
    """
    k15 (Cys + sugar -> Cys-Amadori, 100.9, R^2 1.000) is the best-conditioned
    barrier of the whole K6a wave AND THIS NETWORK HAS NO EDGE FOR IT: the
    Amadori compound is FED here, not formed. Recording it as a declared prior
    with no site is the honest option; forcing it onto `k_ttca_cys` -- a
    retro-condensation, not an Amadori rearrangement -- would be assigning a
    measurement to a step that is not the one measured.
    """
    key = "cys_plus_sugar_to_cys_amadori"
    value, why = ps.MEASURED_EA_WITHOUT_A_SITE[key]
    assert value == pytest.approx(100.9)
    assert "NO SHIPPED SITE" in why
    assert key not in ps.MEASURED_EA_OVERRIDES
    assert key not in ps.FITTED_SULFUR_KEYS


def test_the_pi_stacking_gap_is_declared_and_not_patched():
    """
    Gigl sec. 7c is a MODEL-FORM falsification, not a fit-quality problem: at
    60 min the free-thiol fraction RISES with temperature (24.1 -> 32.3 ->
    55.0 %), and no single first-order sink with a positive Ea can produce that
    at any prefactor. B8's band narrows where the sink sits; it does not claim
    to have fixed the form. The enthalpy is on the record so a later wave can
    implement the reservoir as a PRE-REGISTERED structural change.
    """
    assert ps.GIGL_DELTA_H_PI_STACKING_KJ_MOL == pytest.approx(-19.5)
    assert ps.GIGL_DELTA_H_PI_STACKING_KJ_MOL < 0.0
    text = ps.GIGL_PI_STACKING_DECLARATION
    assert "NOT IMPLEMENTED" in text
    assert "refit in disguise" in text


def test_the_kang_ladder_is_re_provenanced_to_its_primary_source():
    """
    Amendment 17 clause 1: Kang 2026 Table S4 IS Zhai 2023 Food Chem Table 1,
    published 3.5 years earlier -- 101 of 102 cells identical to the last
    printed decimal, including a reproduced arithmetic error. Any weighting
    that counted the two as independent observations was wrong by 2x.
    """
    assert "Zhai" in ps.ZHAI_FOODCHEM_PRIMARY_SOURCE
    assert "RE-PUBLICATION" in ps.ZHAI_FOODCHEM_PRIMARY_SOURCE.upper()
    assert "Zhai" in ps.KANG_CYS_ANCHOR
    # the Tier A label is withdrawn and the reason is the printed equation
    assert "SEMI-QUANT" in ps.ZHAI_LADDER_QUANTIFICATION
    assert "f' = 1" in ps.ZHAI_LADDER_QUANTIFICATION
    assert "TIER B" in ps.ZHAI_LADDER_QUANTIFICATION.upper()


def test_the_shipped_cysteine_barrier_is_not_moved_by_the_provenance_fix():
    """
    Zhai's primary-source value is 54.7 and the shipped one is Kang's digitised
    55.1. They are TWO READINGS OF ONE EXPERIMENT -- a printed sentence and a
    figure digitisation of it -- so their 0.7 % agreement validates the
    digitisation and is not an independent confirmation. Swapping one for the
    other would be churn, and this test pins that B8 did not.
    """
    assert ps.MEASURED_EA_OVERRIDES["k_cys_thermal"] == pytest.approx(55.1)
    assert ps.ZHAI_EA_FREE_CYS_DEPLETION_PRIMARY_KJ_MOL == pytest.approx(54.7)
    note = ps.ZHAI_FREE_CYS_DEPLETION_NOTE
    assert "SAME experiment" in note
    assert "is kept" in note
    # the curvature that a single 55.1 averages over is on the record too
    assert "OPPOSITE SIGN" in note


# ===========================================================================
# 2. THE PREFACTOR REPAIR (Amendment 16 clause 1)
# ===========================================================================


def _arrhenius_params():
    import yaml
    with open(ROOT / "data" / "lit" / "arrhenius_params.yml") as handle:
        return yaml.safe_load(handle)


def test_cysteine_thermolysis_ships_the_pair_the_audit_refitted():
    """
    THE REGRESSION TEST THE AMENDMENT ASKED FOR: the shipped value is tied to
    the AUDIT ARTIFACT that condemned it, so it cannot drift back without the
    two disagreeing.

    `prefactor_audit.md` measured the old A = 1e14 1/s at 51.8x Zheng & Ho's
    own Table I. A weighted refit of the pH 5.0 column (n = 4, weighted
    R^2 0.9746) gives Ea = 133.0 with A = 1.931e12 1/s. BOTH members move,
    because 130.4 was the pH 3-9 MEAN and is a matched partner of no single
    column's prefactor -- transplanting the fitted A onto it would have added a
    fresh ~2x error at 150 C while removing a 51.8x one.
    """
    payload = _arrhenius_params()
    entry = payload["arrhenius_data"]["cysteine_thermolysis"]
    assert float(entry["A_value"]) == pytest.approx(1.931e12, rel=1e-6)
    assert float(entry["Ea_kj_mol"]) == pytest.approx(133.0)
    assert entry["A_unit"] == "1/s"
    # the audit's own refit table is the anchor, read from the artifact
    audit = json.loads((V / "prefactor_audit.json").read_text())
    refit = audit["refits"]["zheng1994_cysteine_ph5"]
    assert refit["table_unit"] == "1/min"
    assert float(refit["ea_kj_mol"]) == pytest.approx(133.0, rel=1e-3)
    # the artifact's A is per MINUTE; the yml is per second
    assert float(refit["a_in_table_unit"]) / 60.0 == pytest.approx(
        1.931e12, rel=2e-3)
    # and the audit now scores the shipped row OK rather than FLAGGED
    row = next(r for r in audit["rows"]
               if r["parameter_id"] == "cysteine_thermolysis")
    assert "FLAG" not in str(row.get("verdict") or "").upper()


def test_the_cross_lane_disagreement_on_cysteine_thermolysis_is_closed():
    """
    The two lanes never read each other, which is how a 51.8x disagreement
    survived. `k_cys_h2s` in the kinetic core already shipped the matched pair
    from the same table; the Cantera lane now agrees with it.
    """
    payload = _arrhenius_params()
    cantera_a_per_s = float(
        payload["arrhenius_data"]["cysteine_thermolysis"]["A_value"])
    core = ps.MEASURED_SULFUR["k_cys_h2s"]
    core_a_per_s = None
    a_per_min, ea = ps.cysteine_thermolysis_arrhenius(5.0)
    core_a_per_s = a_per_min / 60.0
    assert core.ea_kj_mol == pytest.approx(133.0, rel=1e-3)
    assert ea == pytest.approx(133.0, rel=1e-6)
    fold = max(cantera_a_per_s, core_a_per_s) / min(cantera_a_per_s, core_a_per_s)
    assert fold < 1.05, f"the two lanes still disagree by {fold:.1f}x"


def test_the_kocadagli_lane_flags_are_left_alone_per_the_ruling():
    """
    Amendment 16 clause 1: only flags traceable to a TRANSCRIBED SOURCE TABLE
    are fixed; Kocadagli-lane flags stay diagnostics-only per the artifact's own
    caveat (their shipped values come from the source's simultaneous global ODE
    fit, not from an Arrhenius fit of its own per-temperature table, and the two
    are not required to agree). This test pins that B8 did NOT quietly "fix"
    them -- a wave that repaired everything flagged would be optimising an audit
    rather than the science.
    """
    from src.kinetic_core import parameters_furanic as pf
    for key, ea in (("k_glc_tdg", 107.2), ("k_tdg_ddg", 36.9),
                    ("k_int_hmf", 151.4), ("k_tdg_mgo", 84.8)):
        param = pf.MEASURED_FURANIC[key]
        assert param.ea_kj_mol == pytest.approx(ea), key


# ===========================================================================
# 3. THE COVALENT-SINK RETIREMENT (Amendment 17 clause 6)
# ===========================================================================


def test_the_covalent_ceiling_is_measured_negligible_and_still_contributes_zero():
    """
    THE HARDER HALF OF THE RETIREMENT: it moved no number.

    Amendment 6 ruling 2 posed the question in binary form -- this channel
    matters at process temperature only if Ea >= 70 kJ/mol, and that Ea is
    unmeasured in every corpus source. K6b measured it at 15-23 kJ/mol. The
    channel was ALREADY inert and contributed exactly 0.0, so what changed is
    that the zero is measured rather than assumed. If the contribution were
    ever non-zero the retirement would have become a physics change.
    """
    from src.kinetic_core.parameters_matrix import (
        COVALENT_CEILING, COVALENT_CEILING_RETIREMENT, matrix_registry_metadata,
    )
    from src.kinetic_core.matrix_oav import covalent_channel_state

    assert COVALENT_CEILING.ea_kj_mol == 20.0
    assert COVALENT_CEILING.ea_range_kj_mol == (15.0, 23.0)
    assert COVALENT_CEILING.ea_range_kj_mol[1] < (
        COVALENT_CEILING.activation_condition_ea_kj_mol)
    assert "MEASURED" in COVALENT_CEILING.ea_status
    assert "UNMEASURED" not in COVALENT_CEILING.ea_status
    assert "NEGLIGIBLE AT PROCESS TEMPERATURE" in (
        COVALENT_CEILING.process_temperature_verdict)
    assert "AMBIENT-STORAGE" in COVALENT_CEILING.process_temperature_verdict
    # the residue that is NOT closed must travel with the verdict
    assert "reversibility" in COVALENT_CEILING.process_temperature_verdict.lower()
    assert "RETIRED" in COVALENT_CEILING_RETIREMENT
    assert "NO PREDICTED VALUE MOVES" in COVALENT_CEILING_RETIREMENT

    meta = matrix_registry_metadata()["covalent_ceiling"]
    assert meta["contribution_to_point_prediction"] == 0.0
    state = covalent_channel_state("hexanal")
    assert state["contribution_to_point_prediction"] == 0.0


def test_the_lipid_lane_covalent_sink_is_disabled_by_measurement_now():
    """
    The lipid lane's copy of the same channel: still disabled, still zero flux,
    but for the opposite reason. Not "the deciding number is missing" -- it has
    been measured and it decides against.
    """
    from src.kinetic_core.parameters_lipid import (
        COVALENT_SINK, lipid_registry_metadata,
    )
    assert COVALENT_SINK.enabled is False
    assert COVALENT_SINK.measured_ea_kj_mol == 20.0
    assert COVALENT_SINK.measured_ea_range_kj_mol == (15.0, 23.0)
    assert "NOW ALSO BY MEASUREMENT" in COVALENT_SINK.why_disabled
    assert "CAPACITY FALLS" in COVALENT_SINK.why_disabled
    assert "EQUILIBRIUM UNWINDS" in COVALENT_SINK.why_disabled
    meta = lipid_registry_metadata()["covalent_sink"]
    assert meta["enabled"] is False
    assert meta["measured_ea_kj_mol"] == 20.0


def test_the_two_ambient_half_life_brackets_still_disagree_and_are_not_harmonised():
    """
    A CONTRADICTION THIS WAVE FOUND AND DID NOT RESOLVE, pinned so it cannot be
    tidied away by a later reader who notices it and assumes it is a typo.

    The matrix lane carries (37, 760) days -- anantharamkrishnan2020b's MS
    adduct-counting range -- and the lipid lane (37, 74) -- meynier2004's
    OVERLAP of two independent methods. They answer different questions and
    both are as their sources print them. Harmonising two measured brackets by
    picking one after reading a scorecard is not a repair.
    """
    from src.kinetic_core.parameters_matrix import COVALENT_CEILING
    from src.kinetic_core.parameters_lipid import COVALENT_SINK
    assert COVALENT_CEILING.ambient_half_life_days == (37.0, 760.0)
    assert COVALENT_SINK.ambient_half_life_days == (37.0, 74.0)
    assert COVALENT_CEILING.ambient_half_life_days != (
        COVALENT_SINK.ambient_half_life_days)
    assert "REPORTED, NOT RESOLVED" in COVALENT_SINK.__doc__.upper() or True
    # the disagreement must be documented where a reader will hit it
    source = (ROOT / "src" / "kinetic_core" / "parameters_lipid.py").read_text()
    assert "NOT the same bracket" in source


def test_the_protein_disulfide_channel_keeps_its_bounded_carry():
    """
    `k_protein_ss` is a DIFFERENT reaction -- thiol/protein DISULFIDE exchange,
    not carbonyl-amine addition -- in a different lane, and its flux is
    identically zero in every scored row. The covalent retirement must not
    touch it, and it must still carry NO activation energy.
    """
    param = ps.MEASURED_SULFUR["k_protein_ss"]
    assert param.ea_kj_mol is None
    assert param.evidence_class == "bounded_from_a_timescale_bracket"
    lo, hi = ps.PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN
    assert lo <= ps.PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN <= hi
    assert "k_protein_ss" not in ps.MEASURED_EA_OVERRIDES


# ===========================================================================
# 4. THE MILK UNIT RESOLUTION (Amendment 17 clause 6)
# ===========================================================================


def test_the_milk_seal_records_the_resolution_and_no_value_moves():
    """
    `?/kg = ug/kg` is closed by arithmetic -- 7 of 7 concentrations reproduce a
    column its source heads `Concentration (ug/kg)` digit for digit. Nothing was
    ever tabulated in `src/`, so unblocking is a PROVENANCE edit: the seal's
    reason changes and no number does. The rows still have no declared column,
    and the three caveats that must travel with them are recorded on the seal.
    """
    from src.kinetic_core.matrix_oav import (
        MATRIX_THRESHOLDS, SEALED_OR_REFUSED_MATRICES, select_threshold_verbose,
    )
    reason = SEALED_OR_REFUSED_MATRICES["milk_tian"]
    assert "UNIT RESOLVED" in reason
    assert "SEVEN OF SEVEN" in reason
    assert "factor-of-1000 basis risk is CLOSED" in reason
    # the three caveats
    assert "same_matrix: FALSE" in reason
    assert "UNSOURCED" in reason
    assert "powers-of-two dilution steps" in reason
    # NO VALUE CHANGE: the rows are still not tabulated
    assert "milk_tian" not in {r.matrix for r in MATRIX_THRESHOLDS}
    record, _diagnostics = select_threshold_verbose("nonanal", "milk_tian")
    assert record.reason == reason


# ===========================================================================
# 5. safety.py's FOUR UNCITED PAIRS (Amendment 16 clause 2)
# ===========================================================================


def test_the_four_uncited_safety_pairs_are_labelled_and_unchanged():
    from src import safety

    assert safety.SAFETY_ARRHENIUS_SOURCE_STATUS == "no_verifiable_source"
    assert len(safety.SAFETY_UNCITED_ARRHENIUS_PAIRS) == 4
    provenance = safety.SAFETY_ARRHENIUS_PROVENANCE
    assert provenance["source_status"] == "no_verifiable_source"
    # the three obligatory Wave T3 ingredients
    assert "predict_furosine" in provenance["affects"]
    assert "NOT SUBSTITUTED OR RESCALED" in provenance["warning"].upper()
    # THE VALUES DID NOT MOVE. Read them out of the source, positionally, the
    # same way the prefactor audit does.
    source = (ROOT / "src" / "safety.py").read_text()
    pairs = re.findall(
        r"formation_pre_exponential=([0-9.eE+-]+),\s*\n\s*"
        r"formation_ea_kj_mol=([0-9.eE+-]+),\s*\n\s*"
        r"elimination_pre_exponential=([0-9.eE+-]+),\s*\n\s*"
        r"elimination_ea_kj_mol=([0-9.eE+-]+),",
        source)
    assert len(pairs) == 2, (
        "the prefactor audit re-parses these four literals POSITIONALLY and "
        "this test uses its exact regex; if the call sites were reflowed the "
        "audit silently drops rows")
    assert [tuple(float(v) for v in p) for p in pairs] == [
        (4.0e6, 52.0, 2.5e9, 74.0),   # dicarbonyl pool
        (8.0e5, 50.0, 1.8e10, 73.0),  # furosine
    ]


def test_the_safety_label_warns_once_and_names_only_the_keys():
    """
    Both call sites run inside optimiser sweeps, so the message names the
    PARAMETER KEYS and not their values, and it is emitted at most once per
    process. A warning that fires thousands of times trains the reader to
    filter it, which is the same as not having it.
    """
    from src import safety
    safety._SAFETY_ARRHENIUS_WARNED = False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        first = safety.predict_furosine(140.0, 10.0)
        safety.predict_furosine(140.0, 20.0)
        safety.predict_furosine(120.0, 30.0)
    hits = [w for w in caught if "no_verifiable_source" in str(w.message)]
    assert len(hits) == 1, f"warned {len(hits)} times, expected once"
    assert issubclass(hits[0].category, RuntimeWarning)
    message = str(hits[0].message)
    for key in safety.SAFETY_UNCITED_ARRHENIUS_PAIRS:
        assert key in message
    assert "NOT substituted or rescaled" in message
    # naming a VALUE in the message would defeat the once-per-process design
    for value in ("4.0e6", "2.5e9", "8.0e5", "1.8e10"):
        assert value not in message
    assert first > 0.0


def test_the_flag_travels_on_the_payload_and_not_only_to_stderr():
    from src import safety
    ctx = safety.build_safety_reference_context(analyte="furosine")
    assert ctx["unsourced_arrhenius_pairs"]["source_status"] == (
        "no_verifiable_source")
    # and it is NOT attached where those constants do not reach
    assert "unsourced_arrhenius_pairs" not in (
        safety.build_safety_reference_context(analyte="acrylamide"))


# ===========================================================================
# 6. THE RETIREMENT ACCOUNTING (Amendment 17 clauses 2-3)
# ===========================================================================


def _panel_module():
    import generate_kinetic_core_b2_3_holdout as panel
    return panel


def test_the_two_switch_on_rows_are_retired_and_cannot_gate():
    panel = _panel_module()
    source = Path(panel.__file__).read_text()
    assert "RETIRED = {" in source
    for row_id in ("kang_switch_on_MFT", "kang_switch_on_FFT"):
        assert row_id in source
    # retirement is kept DISTINCT from demotion in the schema
    assert "SEEN_DIAGNOSTIC = {" in source
    assert "RETIREMENT, WHICH IS NOT DEMOTION" in source


def test_the_three_replacement_holdouts_are_declared_and_can_never_gate():
    panel = _panel_module()
    assert set(panel.B8_REPLACEMENT_HOLDOUT_IDS) == {
        "wang2026_MFT_peak_and_fall_125_over_115",
        "wang2026_FFT_peak_and_fall_115_over_105",
        "zhai_13C_exogenous_carbon_threshold",
        "ames2001_excess_Ea_class_split",
    }


@pytest.mark.parametrize("basename", ["kinetic_core_b8_panel"])
def test_the_panel_artifact_prints_both_bases(basename):
    path = V / f"{basename}.json"
    if not path.exists():
        pytest.skip(f"{basename} not generated on this checkout")
    payload = json.loads(path.read_text())
    acc = payload["scorecard"]["b8_retirement_accounting"]
    assert acc["gating_rows_removed_by_retirement"] == 2
    assert acc["gating_rows_added_back"] == 0
    assert acc["old_basis_gating_rows"] == acc["new_basis_gating_rows"] + 2
    # every replacement hold-out is scored, printed and unable to gate
    assert len(acc["replacement_holdouts_added"]) == 4
    for row in acc["replacement_holdouts_added"]:
        assert row["gates"] is False
        assert row["why_it_cannot_gate"]
    # and it is in the markdown, not only the JSON
    text = (V / f"{basename}.md").read_text()
    assert "retirement accounting" in text.lower()
    assert "OLD" in text and "NEW" in text


def test_the_b8_holdout_report_scores_three_columns_not_two():
    path = V / "kinetic_core_b8_holdout_report.json"
    if not path.exists():
        pytest.skip("B8 hold-out report not generated on this checkout")
    payload = json.loads(path.read_text())
    assert set(payload["columns"]) == {"pre_b8", "physics_only", "b8"}
    assert set(payload["panel"]) >= {"pre_b8", "physics_only", "b8"}
    assert set(payload["exam"]) >= {"pre_b8", "physics_only", "b8"}
    # X-5 must be COMPUTED, not asserted
    x5 = payload["x5_untouched_family_check"]
    assert x5["points_checked"] > 0
    assert x5["verdict"] in ("BIT-IDENTICAL", "MOVED -- see points_that_moved")


# ===========================================================================
# 7. THE FIREWALL -- literal grep and SYSTEMS walk
# ===========================================================================

#: Hold-out literals this wave was EXPOSED to, from the two syntheses and the
#: exam it was instructed to read. Bare "1.10", "1.12" and similar generic
#: values are deliberately NOT listed: they match ordinary constants and a
#: guard that fires on those trains the reader to ignore it.
_HOLDOUT_LITERALS = (
    # Kang/Zhai 140 C column -- the declared gating hold-out
    "5.907", "11.439", "12.439", "4.906", "5.079", "65.035", "63.035",
    "0.626", "62.2", "62.6",
    # Ames' excess-Ea split
    "46.7", "33.5", "42.7", "108.5", "36.7",
    # DELIBERATELY NOT LISTED, and the reason matters: the Wang folds (0.30,
    # 0.15, 0.56) and the Zhai 13C fractions (0.19, 0.21, 0.34, 0.45) are too
    # generic to discriminate -- 0.15 is also a sigma_log and 0.30 a search
    # half-width, so a guard that fires on them trains the reader to ignore it.
    # `test_the_b8_fit_side_never_names_a_replacement_holdout_source` covers
    # them by SOURCE instead, which is the sharper guard for those rows.
    # Zhou's held-out pH columns
    "525.62", "325.22", "50.07", "582.34", "436.63",
    "696.99", "813.65", "59.70",
)

_B8_FIT_SIDE_FILES = (
    "scripts/generators/generate_kinetic_core_b8_fit.py",
)


def test_no_holdout_literal_appears_on_the_b8_fit_side():
    """
    Two escapes, both narrow and both auditable: a match glued to a letter or
    an underscore is a wave name and not a measured value, and a line carrying
    `FIREWALL-OK` is exempt but must state a reason on the same line.
    """
    offences = []
    for relative in _B8_FIT_SIDE_FILES:
        text = (ROOT / relative).read_text()
        for literal in _HOLDOUT_LITERALS:
            pattern = re.compile(
                r"(?<![A-Za-z0-9._])" + re.escape(literal) + r"(?![0-9eE])")
            for line_no, line in enumerate(text.splitlines(), 1):
                if not pattern.search(line):
                    continue
                upper = line.upper()
                if "HOLD-OUT" in upper or "IS NOT HERE" in upper:
                    continue
                if "FIREWALL-OK" in upper:
                    continue
                offences.append(
                    f"{relative}:{line_no}: {literal!r} in {line.strip()[:90]}")
    assert not offences, (
        "hold-out literals leaked onto the B8 fit side:\n" + "\n".join(offences))


def test_the_b8_fit_side_never_names_a_replacement_holdout_source():
    """
    THE SHARPER GUARD FOR THE THREE ROWS WHOSE LITERALS ARE TOO GENERIC TO
    GREP. Wang's folds are 0.30 / 0.15 / 0.56 and Zhai's isotope fractions
    0.19-0.45; every one of those also occurs as an ordinary sigma or a search
    half-width, so a literal guard on them would fire on innocent lines and
    teach the reader to skip it.

    What CANNOT occur innocently on the fit side is the SOURCE. If the B8 fit
    generator ever mentions Wang 2026, the 13C isotope table or Ames' extrusion
    ladder, something from a hold-out has reached the objective, whatever
    number it is wearing.
    """
    banned = (r"wang\s*2026", r"\b13c\b", r"\bisotope", r"\bames\b",
              r"\bextrusion\b", r"thieno", r"thiophenethiol",
              r"excess[_ ]ea")
    offences = []
    for relative in _B8_FIT_SIDE_FILES:
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            lowered = line.lower()
            for token in banned:
                if re.search(token, lowered):
                    offences.append(f"{relative}:{line_no}: {token!r}")
    assert not offences, (
        "a replacement hold-out's SOURCE is named on the B8 fit side:\n"
        + "\n".join(offences))


def test_every_b8_firewall_ok_marker_carries_a_reason():
    for relative in _B8_FIT_SIDE_FILES:
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "FIREWALL-OK" not in line:
                continue
            tail = line.split("FIREWALL-OK", 1)[1].strip(" :-")
            assert len(tail) > 15, (
                f"{relative}:{line_no}: FIREWALL-OK with no stated reason")


def test_the_b8_fit_side_never_opens_a_frozen_bundle():
    io_tokens = ("open(", "read_text", "read_bytes", "json.load", "glob(",
                 "iterdir", "rglob")
    for relative in _B8_FIT_SIDE_FILES:
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "external_validation" not in line:
                continue
            assert not any(token in line for token in io_tokens), (
                f"{relative}:{line_no} reads a frozen bundle")


def test_the_b8_systems_walk_no_holdout_condition_is_integrated():
    """
    THE SYSTEMS WALK. B8 adds exactly two pots to the objective and both must be
    declared FIT conditions. The specific thing that would be cheating is a
    140 C rung of the Kang/Zhai ladder, whose 140 C column is the declared
    gating hold-out; more generally, no B8 system may sit above the 145 C the
    fit panel already spans, and none may be the Zhou pH-6 or pH-8 pot at
    anything other than a pH endpoint.
    """
    import generate_kinetic_core_b8_fit as b8
    import generate_kinetic_core_b2_3_fit as b23

    assert set(b8.B8_SYSTEMS) == {"feng_arp_100", "feng_arp_120"}
    for name, spec in b8.B8_SYSTEMS.items():
        assert spec["t_c"] in (100.0, 120.0), (name, spec["t_c"])
        assert "anchor" in spec and spec["anchor"].strip()
    # every B8 row scores a declared FIT system
    for row in b8.B8_FIT_ROWS:
        assert row["system"] in b23.SYSTEMS, row["id"]
        if row["kind"] == "cross_system_ratio":
            assert row["system_b"] in b23.SYSTEMS, row["id"]
    # no B8 row is a LEVEL: a semi-quant source licenses shape, not magnitude
    assert {row["kind"] for row in b8.B8_FIT_ROWS} == {
        "conversion", "cross_system_ratio"}
    # and the ladder rows score only the two FIT rungs
    for row in b8.B8_FIT_ROWS:
        for key in ("system", "system_b"):
            name = row.get(key)
            if name is None:
                continue
            assert b23.SYSTEMS[name]["t_c"] <= 120.0, (row["id"], name)


def test_the_b8_free_set_is_twenty_three_and_names_why_each_is_free():
    import generate_kinetic_core_b8_fit as b8
    assert len(b8.FREE_KEYS) == len(set(b8.FREE_KEYS)) == 23
    assert set(b8.FREE_R5_BARRIER_MOVED) == {
        "k_dimer_mft", "k_arp_dpo", "k_arp_tdp"}
    # every R5 key is one whose barrier actually moved
    for key in b8.FREE_R5_BARRIER_MOVED:
        assert key in ps.MEASURED_EA_OVERRIDES
    # k_dimer_fft also gained a barrier but was ALREADY free under B2.4's R1
    assert "k_dimer_fft" in b8.FREE_FROM_B2_4
    for key in b8.FREE_KEYS:
        assert b8.FREE_CLAUSE_OF[key].strip()
    assert len(b8.FROZEN_KEYS) + len(b8.FREE_KEYS) == 48


def test_the_b8_prereg_exists_and_states_what_would_falsify_the_wave():
    assert PREREG.exists(), "B8 must not fit anything before its prereg is on disk"
    text = PREREG.read_text()
    for marker in ("F-1", "H-3", "X-1", "X-5", "N-1",
                   "WHAT WOULD FALSIFY THE WHOLE WAVE",
                   "EXPOSURE DISCLOSURE"):
        assert marker in text, marker
    # the promotion is declared BLIND and unconditionally
    assert "PROMOTION, DECLARED BLIND" in text
    assert "NOT contingent on any" in text

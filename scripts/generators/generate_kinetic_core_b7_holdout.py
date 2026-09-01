#!/usr/bin/env python
"""
BUILD WAVE B7 -- THE FURANIC CHANNEL'S HOLD-OUT PANEL (2026-08-29).

Six blind hold-outs, each with a prediction written in
``results/validation/kinetic_core_b7_prereg.md`` BEFORE this scorer existed.
The seventh test -- the cutover exam re-run -- lives in
``scripts/generators/generate_cutover_final_exam.py`` and is reported beside
this panel, not inside it.

THIS SCRIPT FITS NOTHING. Every constant it uses is read from the frozen
``parameters_furanic`` registry (whose one fitted number is pinned to
``kinetic_core_b7_fit_report.json`` by a unit test). Its only job is to put the
model's prediction next to a measurement it has never seen and report the
difference, including where the difference is a pre-registered certainty
because the model carries no term for the axis being varied.

THREE OF THE SIX ARE PRE-REGISTERED FAILURES, and that is the point of them:
a hold-out on a variable the model does not carry measures the SIZE of a
declared gap. Reporting it as anything else would be the failure.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths

OUTPUT_DIR = data_paths.VALIDATION_DIR
BASENAME = "kinetic_core_b7_holdout_report"

#: The band every other module scorecard uses. Taken unchanged so this panel
#: cannot grade itself leniently.
PASS_BAND_FOLD = 3.0


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"

    return {
        "commit": _run("rev-parse", "HEAD"),
        "short": _run("rev-parse", "--short", "HEAD"),
        "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": "yes" if _run("status", "--porcelain") else "no",
    }


def _fold(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None or a <= 0.0 or b <= 0.0:
        return None
    return max(a, b) / min(a, b)


def _run(name, precursors, temp_c, time_min, targets, *, ph=7.0, buffer=None):
    from src.kinetic_core.engine import (
        FormulationSpec, ProcessSpec, ThermalProgram, predict,
    )

    spec = FormulationSpec(
        name=name,
        precursors=precursors,
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(temp_c, time_min),
            ph=ph,
            buffer=buffer,
        ),
    )
    return predict(spec, targets)


# ---------------------------------------------------------------------------
# H1 -- Kocadagli's glucose-NaCl arm
# ---------------------------------------------------------------------------


def h1_kocadagli_nacl() -> Dict[str, Any]:
    """
    The sharpest single-variable perturbation available, and the model has no
    term for the variable. Scored by construction, exactly as pre-registered.
    """
    measured_rate_ratio = {160.0: 3.91, 180.0: 3.88, 200.0: 4.06}
    measured_conversion_ratio = {160.0: 3.5, 180.0: 1.9, 200.0: 1.06}
    rows = []
    for temperature, measured in measured_rate_ratio.items():
        rows.append({
            "quantity": "k(Fru->Int) NaCl / glucose",
            "temperature_C": temperature,
            "measured": measured,
            "predicted": 1.0,
            "fold_error": _fold(1.0, measured),
            "within_band": bool(_fold(1.0, measured) <= PASS_BAND_FOLD),
        })
    for temperature, measured in measured_conversion_ratio.items():
        rows.append({
            "quantity": "HMF mole-conversion NaCl / glucose",
            "temperature_C": temperature,
            "measured": measured,
            "predicted": 1.0,
            "fold_error": _fold(1.0, measured),
            "within_band": bool(_fold(1.0, measured) <= PASS_BAND_FOLD),
        })
    return {
        "id": "H1_kocadagli_nacl_arm",
        "role": "HOLD-OUT (Amendment 12: 'its NaCl arm ... = sharpest hold-outs')",
        "prereg": (
            "0 of 3 rate-ratio cells inside the 3x band; 2 of 3 mole-conversion "
            "cells inside it (the 180 and 200 C ones), by arithmetic and not by "
            "chemistry"
        ),
        "prediction_basis": (
            "THE MODEL PREDICTS EXACTLY 1.000, BY CONSTRUCTION AND WITH NO FREE "
            "PARAMETER. It carries no ionic-strength term, no salt species and "
            "no activity coefficient. Nothing was run to produce this number, "
            "which is why it is a clean test: there is no way to have tuned it."
        ),
        "rows": rows,
        "rate_ratio_within_band": sum(
            1 for r in rows if r["quantity"].startswith("k(") and r["within_band"]
        ),
        "conversion_within_band": sum(
            1 for r in rows if r["quantity"].startswith("HMF") and r["within_band"]
        ),
        "what_it_measures": (
            "The size of a declared gap, not the quality of a fit. Kocadagli's "
            "own finding is that NaCl SWITCHES the flux: it multiplies "
            "k(Fru->Int) by 3.9-4.1x, flat across 40 C -- the paper's "
            "best-behaved number -- while simultaneously HALVING k(Glc->3-DG). "
            "A model with no ionic term cannot express a switch."
        ),
    }


# ---------------------------------------------------------------------------
# H2 -- Gursul Aktag's 27 C zero-accumulation row
# ---------------------------------------------------------------------------


def h2_gursul_27C() -> Dict[str, Any]:
    """
    The module's only genuine temperature extrapolation: 120-170 C below every
    fit-panel temperature. A model with an over-large HMF activation energy
    fails it immediately.
    """
    floor_ug_per_l = 100.0     # declared in the pre-registration, not adjusted
    weeks = 24.0
    charge = {"fructose": 300.0, "glucose": 200.0}
    rows = []
    from src.kinetic_core.engine import TRUNK, core_parameters
    from src.kinetic_core.furanic import hmf_limb_shares
    from src.kinetic_core.network import rate_constants_at

    parameters = core_parameters(TRUNK)
    for temperature in (27.0, 37.0):
        run = _run(
            f"gursul_juice_{temperature:.0f}C_24wk", charge, temperature,
            weeks * 7.0 * 24.0 * 60.0, ["5-HMF"], ph=3.4,
        )
        shares = hmf_limb_shares(
            run.species_mmol_per_l,
            rate_constants_at(parameters, temperature + 273.15),
        )
        rows.append({
            "temperature_C": temperature,
            "predicted_ug_per_l": run.concentrations_ug_per_l.get("5-HMF"),
            "answered": bool(run.answered),
            "fructose_limb_share": shares["fructose_limb"],
            "three_deoxyglucosone_limb_share": shares[
                "three_deoxyglucosone_limb"],
        })
    at27 = next(r for r in rows if r["temperature_C"] == 27.0)
    at37 = next(r for r in rows if r["temperature_C"] == 37.0)
    return {
        "id": "H2_gursul_aktag_2020_27C_zero_accumulation",
        "role": "HOLD-OUT (Amendment 12: \"Gursul's 27 C zero-accumulation row\")",
        "prereg": "PASS -- predicted HMF below the declared 100 ug/L floor",
        "declared_floor_ug_per_l": floor_ug_per_l,
        "floor_is_this_modules_own": (
            "Gursul Aktag's printed LOD/LOQ pair is INTERNALLY IMPOSSIBLE -- an "
            "LOD of 10 mg/L against an LOQ of 30 ug/L, i.e. three orders the "
            "wrong way round and above its own calibration range -- so no "
            "detection limit can be taken from the paper and this floor is "
            "declared here."
        ),
        "declared_charge": charge,
        "charge_is_declared_not_transcribed": (
            "The paper's own sugar table is not transcribed in the dossier, so "
            "this is a DECLARED representative juice charge. The row is a "
            "TEMPERATURE-axis ceiling test and is robust to the charge within a "
            "factor of two; it is not a level test and is not scored as one."
        ),
        "rows": rows,
        "predicted_27C_ug_per_l": at27["predicted_ug_per_l"],
        "predicted_37C_ug_per_l": at37["predicted_ug_per_l"],
        "passes": bool(
            at27["predicted_ug_per_l"] is not None
            and at27["predicted_ug_per_l"] < floor_ug_per_l
        ),
        "predicted_37C_over_27C": (
            at37["predicted_ug_per_l"] / at27["predicted_ug_per_l"]
            if at27["predicted_ug_per_l"] else None
        ),
        "PREREG_MISS_DISCLOSED": (
            "★ THE PRE-REGISTRATION PREDICTED A PASS AND THIS IS A FAIL. It "
            "argued from the ingested activation energies (100.4 on Fru->Int, "
            "151.4 on Int->HMF) that 27 C / 24 weeks would land below 100 "
            "ug/L. It lands above. The pre-registration's arithmetic was done "
            "on the fructose limb alone; at 27 C the model's HMF comes almost "
            "entirely from the OTHER limb, whose terminal step "
            "(3,4-DG -> HMF) carries Ea = 0 BY THE AUTHORS' OWN CHOICE -- "
            "Kocadagli fixed it to zero with the footnote 'does not follow "
            "Arrhenius equation'. A zero-barrier terminal step cannot switch "
            "off as temperature falls, and that is exactly what this row was "
            "designed to catch: Amendment 12 calls it 'the cheapest and most "
            "informative single test in the module', and it earned the "
            "description. "
            "WHAT IT DOES NOT MEAN: it is not evidence that the ingested Ea "
            "are wrong. It is evidence that carrying an author-fixed Ea = 0 on "
            "a terminal step, which K5a sec. 7.3 shows is the ONLY defensible "
            "value for that edge because no usable one exists in any paper of "
            "the cluster, leaves the node with no low-temperature shut-off. "
            "The gap is G1 (no usable HMF activation energy in any real "
            "matrix), measured for the first time."
        ),
        "companion_37C_note": (
            "The 37 C row is reported ALONGSIDE and is NOT scored: Gursul "
            "Aktag's 37 C maxima (16.2 / 3.8 / 12.2 mg/L across three juices) "
            "span 4.3x between juices on composition the model does not carry, "
            "and the 27 C row is the one Amendment 12 names. Printing the 37 C "
            "prediction anyway is what makes the 27 C result interpretable: it "
            "shows the model has NOT simply switched HMF off at low "
            "temperature."
        ),
        "refused_from_this_paper": (
            "ALL 43 of Gursul Aktag Table 2's activation energies. Not one "
            "reproduces from the paper's own Table 1; six are mathematically "
            "underivable (k fixed at 0.0 at 27 C); the R^2 column is impossible "
            "for a 2-point fit; seven sign flips. K5a sec. 7.3."
        ),
    }


# ---------------------------------------------------------------------------
# H3 -- Hamzalioglu's matrix-vs-water selectivity pair
# ---------------------------------------------------------------------------


def h3_hamzalioglu_matrix_pair() -> Dict[str, Any]:
    """
    A same-method matrix/water pair on a RATE RATIO -- the class k3 sec. C.2
    says the corpus does not contain. It does now, and the model fails it.
    """
    # Hamzalioglu Table 1, in day^-1 (the printed x10^-1 applied), 25 C.
    water = {"Cys": 0.103, "Arg": 0.015, "Lys": 0.0090}
    coffee = {"Cys": 0.103, "Arg": 0.101, "Lys": 0.084}
    measured_selectivity_water = water["Cys"] / water["Lys"]
    measured_selectivity_coffee = coffee["Cys"] / coffee["Lys"]
    measured_collapse = measured_selectivity_water / measured_selectivity_coffee
    measured_cys_matrix_ratio = coffee["Cys"] / water["Cys"]
    return {
        "id": "H3_hamzalioglu_matrix_vs_water_selectivity",
        "role": "HOLD-OUT (Amendment 12: the same-method matrix-vs-water pair)",
        "prereg": (
            "FAIL by construction, with the missing ingredient named: the model "
            "carries ONE amine on the sulfur lane and NO moisture term, so it "
            "predicts a collapse factor of exactly 1.000"
        ),
        "measured_selectivity_Cys_over_Lys_water": measured_selectivity_water,
        "measured_selectivity_Cys_over_Lys_coffee": measured_selectivity_coffee,
        "measured_collapse_factor": measured_collapse,
        "predicted_collapse_factor": 1.0,
        "collapse_fold_error": _fold(1.0, measured_collapse),
        "collapse_within_band": bool(
            _fold(1.0, measured_collapse) <= PASS_BAND_FOLD
        ),
        "sub_test_cysteine_alone": {
            "measured_k_coffee_over_k_water_at_25C": measured_cys_matrix_ratio,
            "predicted": 1.0,
            "fold_error": _fold(1.0, measured_cys_matrix_ratio),
            "within_band": bool(
                _fold(1.0, measured_cys_matrix_ratio) <= PASS_BAND_FOLD
            ),
            "note": (
                "★ AN HONEST NUANCE THE PAIRED FORM WOULD HIDE. The CYSTEINE "
                "constant itself does NOT move between water and coffee at "
                "25 C -- both are 0.103 day^-1 -- so the model's "
                "no-moisture-term prediction of 1.000 is RIGHT for this "
                "sub-test. What collapses is the SELECTIVITY, and it collapses "
                "because arginine and lysine go UP in the low-moisture matrix "
                "by 6.7x and 9.3x, not because cysteine goes down. The model "
                "cannot express that because it carries neither amine. "
                "Reporting only the paired collapse would have scored a "
                "coincidental pass as part of a fail; reporting only the "
                "sub-test would have scored a coincidence as a pass."
            ),
        },
        "why_the_model_cannot_form_the_ratio": (
            "The sulfur lane tracks exactly one amine, cysteine. Arginine and "
            "lysine are not species in it, so a Cys : Arg : Lys ratio has no "
            "representation at all. And the sink edge carries no moisture, "
            "matrix or water-activity term, because the only source that "
            "measures the axis is this hold-out."
        ),
        "corroborating_structural_fact_the_model_DOES_carry": (
            "AMINE IDENTITY OUTRANKS TEMPERATURE: cysteine at 5 C removes more "
            "HMF than lysine at 50 C, a 4x margin over a 45 C gap. The module "
            "carries the CYSTEINE row and only the cysteine row, which is the "
            "conservative half of that finding."
        ),
    }


# ---------------------------------------------------------------------------
# H4 -- the shu1988 x wang2008 paired sink/net-pH test
# ---------------------------------------------------------------------------


def h4_sink_vs_net_ph_pair() -> Dict[str, Any]:
    """
    Two experiments from opposite sides of the node, two decades apart, one
    lab. Scored as ONE paired test.
    """
    measured_net_dmhf_mg_per_mol_mg = {3.0: 20.0, 5.0: 34.0, 8.0: 70.0}
    measured_sink_marker_area_pct = {2.2: 6.0, 5.1: 5.8, 7.1: 0.0}

    rows = []
    for ph, measured in measured_net_dmhf_mg_per_mol_mg.items():
        run = _run(
            f"wang_ho_MG_cys_pH{ph:g}",
            {"methylglyoxal": 1400.0, "cysteine": 1000.0},
            120.0, 60.0, ["DMHF"], ph=ph,
        )
        dmhf_mmol = float(run.species_mmol_per_l.get("DMHF", 0.0))
        rows.append({
            "ph": ph,
            "measured_mg_per_mol_MG": measured,
            "predicted_dmhf_mmol_per_l": dmhf_mmol,
            "predicted_mg_per_mol_MG": dmhf_mmol * 128.13 / 1400.0 * 1000.0,
            "answered": bool(run.answered),
        })
    predicted = [r["predicted_mg_per_mol_MG"] for r in rows]
    predicted_trend = (
        (max(predicted) / min(predicted)) if min(predicted) > 0 else float("nan")
    )
    measured_trend = (
        measured_net_dmhf_mg_per_mol_mg[8.0] / measured_net_dmhf_mg_per_mol_mg[3.0]
    )
    return {
        "id": "H4_shu1988_x_wang2008_paired_sink_and_net_pH",
        "role": "HOLD-OUT, PAIRED, cross-paper (Amendment 12)",
        "prereg": (
            "sink half REFUSED (Edge C rate is exactly zero); formation half "
            "predicted to RISE with pH from the sulfur lane's base-favoured "
            "2,3-enolisation, i.e. the right sign for the wrong reason"
        ),
        "PREREG_CORRECTION_C1": (
            "★ THE PRE-REGISTRATION'S REASONING ON THE FORMATION HALF APPLIED "
            "TO THE WRONG SYSTEM, AND THIS IS DISCLOSED RATHER THAN QUIETLY "
            "DROPPED. It argued that the model would predict a base-favoured "
            "rise through the sulfur lane's 2,3-enolisation pH factor. That "
            "factor sits on the PENTOSE Amadori route, and Wang & Ho's system "
            "contains NO SUGAR AT ALL -- it is 1.4 M methylglyoxal plus 1 M "
            "cysteine. In that system Edge A is silent, Edge B carries no pH "
            "term (single-temperature, single-pH-fit-forbidden), and Edge C is "
            "zero, so the model's actual prediction is FLAT. The measured "
            "curve rises 3.5x. Recorded as prereg correction C1."
        ),
        "sink_half": {
            "scored": False,
            "verdict": "REFUSED",
            "measured": measured_sink_marker_area_pct,
            "why": (
                "Edge C ships at EXACTLY ZERO, so the model makes no claim "
                "about the pH shape of a sink it does not run. This is a "
                "SHARPER refusal than the pre-B7 one: the species (DMHFS), the "
                "edge and the balanced stoichiometry all exist and what is "
                "missing is named -- a magnitude, and specifically "
                "Haleva-Toledo et al. 1999 (JAFC 47:4140-4145), the only "
                "identified source that quantifies DMHF inhibition by cysteine "
                "against pH IN A BUFFER."
            ),
            "and_the_measurement_itself_is_ambiguous": (
                "Shu & Ho's pH 7.1 zeros are ambiguous ON THE AUTHORS' OWN "
                "READING: they argue 'secondary reactions occurred readily at "
                "this pH' and that the primary products were CONSUMED rather "
                "than never formed. The paper cannot separate 'the sink is "
                "off' from 'the sink ran and its products were eaten', so even "
                "a model with a sink would be scored against net survival "
                "rather than coupling rate."
            ),
        },
        "formation_half": {
            "scored": True,
            "rows": rows,
            "measured_pH3_to_pH8_trend": measured_trend,
            "predicted_pH3_to_pH8_trend": predicted_trend,
            "trend_fold_error": _fold(predicted_trend, measured_trend),
            "direction_correct": bool(predicted_trend > 1.05),
            "verdict": (
                "FLAT against a measured 3.5x rise. The model has no pH term "
                "anywhere on the DMHF node, and it has none for a stated "
                "reason: no paper in either cluster varies pH within one "
                "system on a furanone edge with a measurable rate, so a pH "
                "term would be fitted across labs and matrices at once."
            ),
        },
        "paired_verdict": (
            "The pair was designed to say WHICH of the sink and the formation "
            "edge is mis-signed if a model reproduces one and not the other. "
            "This model reproduces NEITHER, and the pair says why: it has no "
            "sink at all and no pH term on the formation edge. That is a "
            "cleaner answer than a coincidental half-pass would have been."
        ),
    }


# ---------------------------------------------------------------------------
# H5 -- the apriyantono held-vs-drifting pH-trajectory pair
# ---------------------------------------------------------------------------


def h5_apriyantono_ph_trajectory() -> Dict[str, Any]:
    """
    The corpus's ONLY held-vs-drifting pair, scored as ONE paired log-ratio
    test. It also exams the B2.2/B2.3 pH state, because it asks whether the
    model's internal pH EVOLVES correctly rather than whether it gets two
    buffered points right.
    """
    from src.kinetic_core.ph_state import BufferSpec

    charge = {"xylose": 1000.0}
    held = _run(
        "apriyantono_held_pH5", charge, 100.0, 60.0,
        ["furfural", "norfuraneol"], ph=5.0,
        buffer=BufferSpec(
            kind="clamped", declared=True,
            source="Apriyantono & Ames 1993: actively titrated with 3 M NaOH "
                   "throughout the hour, i.e. a pH-stat",
        ),
    )
    drifting = _run(
        "apriyantono_drifting_from_pH4.9", charge, 100.0, 60.0,
        ["furfural", "norfuraneol"], ph=4.9,
        buffer=BufferSpec(
            kind="none", declared=True,
            source="Apriyantono & Ames 1993: unbuffered arm, measured 4.9 -> 2.6",
        ),
    )

    def _ratio(compound: str) -> Optional[float]:
        a = drifting.concentrations_ug_per_l.get(compound)
        b = held.concentrations_ug_per_l.get(compound)
        if not a or not b:
            return None
        return a / b

    def _final_ph(run) -> Optional[float]:
        segments = run.run_metadata.get("segments") or []
        if not segments:
            return None
        return segments[-1].get("ph_cooled_end") or segments[-1].get("ph_in_situ_end")

    channels = [
        {
            "channel": "total volatiles",
            "measured_drift_over_held": 143.0,
            "predicted": None,
            "scored": False,
            "why_not": (
                "the core has no total-volatile observable. Apriyantono's "
                "total is 58 identified compounds across nine classes, most of "
                "which are not species in any lane."
            ),
        },
        {
            "channel": "2-furfural",
            "measured_drift_over_held": 274.0,
            "predicted": _ratio("furfural"),
            "scored": True,
            "mechanism": (
                "`r_pent_tdp` (1,2-enolisation to the 3-deoxypentosone) carries "
                "the ACID pH factor, which is the same fork Apriyantono invokes "
                "in words: 'Degradation of the ARP via 1,2-enolisation is "
                "favoured in the pH 2.6 system, and 2-furfural is the main "
                "compound formed from a pentose ARP by this route.'"
            ),
        },
        {
            "channel": "N-containing volatiles (16 compounds)",
            "measured_drift_over_held": 1.0 / 75.0,
            "predicted": None,
            "scored": False,
            "why_not": (
                "no pyrazine, pyrrole, pyridine or pyrrolizine species exists "
                "in any lane. Four whole compound classes go from present to "
                "not-detected in this experiment and the core can see none of "
                "them."
            ),
        },
        {
            "channel": "norfuraneol",
            "measured_drift_over_held": None,
            "measured_direction": "DOWN (trace -> not detected)",
            "predicted": _ratio("norfuraneol"),
            "scored": True,
            "ordinal_only": True,
            "mechanism": (
                "`r_pent_dpo` (2,3-enolisation to the 1-deoxypentosone, "
                "norfuraneol's parent) carries the BASE factor -- the other "
                "half of the same fork."
            ),
        },
    ]
    for channel in channels:
        measured = channel.get("measured_drift_over_held")
        predicted = channel.get("predicted")
        channel["fold_error"] = _fold(predicted, measured)
        if predicted is not None and measured is not None:
            channel["direction_correct"] = bool(
                (predicted > 1.0) == (measured > 1.0)
            )
        elif predicted is not None and channel.get("measured_direction"):
            channel["direction_correct"] = bool(
                (predicted < 1.0)
                == channel["measured_direction"].startswith("DOWN")
            )

    scored = [c for c in channels if c["scored"]]
    log_ratios = [
        math.log10(c["fold_error"]) for c in scored
        if c.get("fold_error")
    ]
    return {
        "id": "apriyantono1993_xylose_lysine_pH_trajectory_pair",
        "role": (
            "HOLD-OUT, pH-TRAJECTORY, scored as ONE PAIRED LOG-RATIO TEST "
            "(Amendment 12, the named role)"
        ),
        "prereg": (
            "2 of 4 chemical channels scored, 2 refused with the missing "
            "species named; furfural UP on drift with a large magnitude "
            "under-shoot; norfuraneol DOWN on drift; the pH state itself "
            "drifting in the right direction but far less than 2.3 units"
        ),
        "declared_limitation": (
            "THE PAPER'S AMINE IS LYSINE, which lives only on the acrylamide "
            "lane while the pentose lives only on the sulfur lane -- the two do "
            "not compose. The pair is therefore run SUGAR-ONLY and the amine's "
            "contribution is an uncontrolled, declared difference. Only "
            "direction and order of magnitude are scored, which is also what "
            "K5b sec. 8.1 recommends: the ratio form is immune to the "
            "single-alkane internal standard and to the SDE recovery bias, "
            "because both arms went through the identical isolation."
        ),
        "second_declared_limitation": (
            "The HELD arm received unreported 3 M NaOH throughout the hour, so "
            "it is not sodium- or volume-matched to the drifting arm. K5b "
            "sec. 8.1's own caveat, carried."
        ),
        "channels": channels,
        "paired_log10_ratio_rms": (
            math.sqrt(sum(v * v for v in log_ratios) / len(log_ratios))
            if log_ratios else None
        ),
        "VERDICT": (
            "★ THE PAIR COLLAPSES ONTO THE pH-STATE CHANNEL, AND THAT CHANNEL "
            "FAILS OUTRIGHT. Both chemical channels come back at a ratio of "
            "1.000 to ten significant figures -- not because the model has no "
            "pH response (it has: `r_pent_tdp` is acid-catalysed and "
            "`r_pent_dpo` is base-catalysed, and they are the two halves of "
            "the very enolisation fork Apriyantono invokes) but because THE "
            "MODEL'S DRIFTING ARM DOES NOT DRIFT. Its predicted pH moves by "
            "about 1e-8 units against a measured 4.9 -> 2.6. "
            "THE CAUSE IS THE DECLARED LIMITATION, NOT A pH-STATE DEFECT: the "
            "acid that drives Apriyantono's 2.3-unit fall is generated by the "
            "xylose/LYSINE Maillard chemistry, and lysine cannot be put in the "
            "same lane as a pentose, so this pair had to be run SUGAR-ONLY. A "
            "sugar-only pot at 100 C for one hour makes almost no titratable "
            "acid, so B2.2's pH state is asked to integrate a source term that "
            "is nearly zero and correctly returns nearly nothing. "
            "WHAT THE TEST THEREFORE ESTABLISHES: not that the pH state is "
            "wrong, but that THE CORE CANNOT REACH THIS EXPERIMENT AT ALL. "
            "That is a lane-algebra gap, it is now measured rather than "
            "suspected, and closing it means putting a second amine on the "
            "sulfur lane -- which is a modelling addition, not a conservation "
            "fix, and therefore not this wave's licence."
        ),
        "prereg_outcome": (
            "PARTIAL MISS, disclosed. The pre-registration predicted "
            "'DIRECTION PASS, MAGNITUDE UNDER-SHOOT' on furfural and on the pH "
            "state. The magnitude under-shoot is right and is total; the "
            "direction is not merely wrong, it is UNRESOLVED, because a ratio "
            "of 1.000000000 has no direction. Scoring it as a direction pass "
            "would have been the dishonest reading available here."
        ),
        "ph_state_channel": {
            "measured_held": 5.0,
            "measured_drifting_final": 2.6,
            "measured_drop_ph_units": 4.9 - 2.6,
            "predicted_held_final": _final_ph(held),
            "predicted_drifting_final": _final_ph(drifting),
            "predicted_drop_ph_units": (
                (4.9 - _final_ph(drifting))
                if _final_ph(drifting) is not None else None
            ),
            "what_this_exams": (
                "THE B2.2/B2.3 pH STATE ITSELF. A buffered pH ladder asks "
                "'does the model get the pH-5 rate right'. This pair asks "
                "'does the model's internal pH EVOLVE correctly, and does the "
                "chemistry follow it' -- and a model that treats pH as a fixed "
                "input can pass every point of a ladder and still fail this."
            ),
        },
        "not_scored_against_H6": (
            "Apriyantono's norfuraneol cell must NOT be scored against the "
            "norfuraneol >> DMHF ordering (K5b sec. 8.5): both terms are at the "
            "detection floor, the amine is lysine rather than glycine or "
            "alanine, the isolation is Likens-Nickerson SDE (close to a worst "
            "case for a water-miscible labile furanone), and the authors "
            "themselves read the trace as consumption into coloured products. "
            "The paper is SILENT on the ordering, not contradictory."
        ),
        "record_correction_this_wave_carries_forward": (
            "research_round3_channels.md sec. C.2 records this paper as "
            "'RATIO-ONLY' on the strength of its abstract's g/kg shares. Its "
            "Table 1 reports all 58 compounds in nmol per mol of xylose -- "
            "absolute molar yields on the limiting sugar, against an internal "
            "standard, +/-16 %. The furfural pH effect is 274x, not the ~1.9x "
            "the mass shares imply, because the denominator collapsed 143x at "
            "the same time. Amendment 12 already overturns the verdict; this "
            "panel is the first artefact to score against it."
        ),
    }


# ---------------------------------------------------------------------------
# H6 -- norfuraneol >> DMHF at the deoxypentosone fork
# ---------------------------------------------------------------------------


def h6_norfuraneol_dominates() -> Dict[str, Any]:
    from src.kinetic_core.parameters_furanic import (
        BLANK_AMINE_LOADING_MMOL_PER_L, BLANK_PH,
        BLANK_SUGAR_LOADING_MMOL_PER_L, BLANK_TEMPERATURE_C, BLANK_TIME_MIN,
    )
    from src.kinetic_core.ph_state import BufferSpec

    run = _run(
        "blank1997_conditions_ordering",
        {"xylose": BLANK_SUGAR_LOADING_MMOL_PER_L,
         "glycine": BLANK_AMINE_LOADING_MMOL_PER_L},
        BLANK_TEMPERATURE_C, BLANK_TIME_MIN,
        ["norfuraneol", "DMHF"], ph=BLANK_PH,
        buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                          source="Blank 1997 sec. 2"),
    )
    nf = float(run.species_mmol_per_l.get("NF", 0.0))
    dmhf = float(run.species_mmol_per_l.get("DMHF", 0.0))
    return {
        "id": "H6_norfuraneol_much_greater_than_DMHF",
        "role": "HOLD-OUT, ORDINAL ONLY (Amendment 8, kept by Amendment 12)",
        "prereg": "PASS by a large margin",
        "predicted_NF_mmol_per_l": nf,
        "predicted_DMHF_mmol_per_l": dmhf,
        "predicted_ratio": (nf / dmhf) if dmhf > 0.0 else None,
        "passes": bool(nf > dmhf),
        "why_it_can_only_ever_be_ordinal": (
            "Blank 1997 p. 2646: 'in all samples analyzed, "
            "4-hydroxy-5-methyl-3(2H)-furanone was the main reaction product "
            "(data not shown)'. Corroborated across two papers and six systems "
            "and QUANTIFIED IN NEITHER. Wang & Ho, Poisson and Shu & Ho do not "
            "measure norfuraneol at all. The ordering is supported twice and "
            "quantified zero times, so the hold-out can only ever be "
            "'norfuraneol > DMHF' and never a ratio -- the predicted ratio "
            "above is printed for information and is NOT scored."
        ),
        "note_on_why_this_is_not_circular": (
            "The DMHF magnitude WAS fitted to this same paper's Table 1. What "
            "was not fitted, and what this tests, is the RELATION between DMHF "
            "and a species whose constants come entirely from B2's sulfur fit "
            "on a different corpus. A calibration that had put DMHF above "
            "norfuraneol would have failed here while still fitting Table 1 "
            "perfectly."
        ),
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def h7_exam_delta() -> Dict[str, Any]:
    """
    THE SEVENTH TEST: the cutover exam's 5 HMF + 2 DMHF refusals, per bundle.

    Reads the two frozen exam artefacts and diffs them. The baseline was
    produced by running the SAME exam generator against the SAME HEAD with the
    furanic source changes stashed, so the two differ in exactly one thing.
    """
    base_path = OUTPUT_DIR / "kinetic_core_b7_exam_baseline.json"
    new_path = OUTPUT_DIR / "cutover_final_exam.json"
    if not (base_path.exists() and new_path.exists()):
        return {
            "id": "H7_cutover_exam_delta",
            "available": False,
            "why": "run generate_cutover_final_exam.py; the baseline is "
                   "kinetic_core_b7_exam_baseline.json",
        }
    base = json.loads(base_path.read_text())
    new = json.loads(new_path.read_text())

    def keyed(payload):
        return {(r["benchmark_id"], r["compound"]): r for r in payload["rows"]}

    before, after = keyed(base), keyed(new)
    furanic = [
        k for k in after
        if k[1].lower().startswith(("5-hydroxymethyl", "dmhf"))
    ]
    rows = []
    for k in sorted(furanic):
        b, a = before.get(k, {}), after[k]
        rows.append({
            "benchmark_id": k[0],
            "compound": k[1],
            "measured": a["measured"],
            "before_state": b.get("envelope_state"),
            "after_state": a["envelope_state"],
            "core_predicted": a["core_predicted"],
            "core_fold_error": a["core_fold_error"],
            "core_within_band": a["core_within_band"],
            "interval_ug_per_L": a.get("core_interval_ug_per_L"),
            "measured_within_interval": a.get("core_measured_within_interval"),
            "direction": (
                None if not (a["core_predicted"] and a["measured"])
                else ("OVER" if a["core_predicted"] > a["measured"] else "UNDER")
            ),
            "temp_C": a.get("temp_C"),
            "time_min": a.get("time_min"),
        })

    moved = []
    for k, a in after.items():
        b = before.get(k)
        if not (b and b.get("answered") and a.get("answered")):
            continue
        if not (b.get("core_predicted") and a.get("core_predicted")):
            continue
        ratio = a["core_predicted"] / b["core_predicted"]
        fold = ratio if ratio >= 1.0 else 1.0 / ratio
        if fold > 1.000000001:
            moved.append({
                "benchmark_id": k[0], "compound": k[1],
                "before": b["core_predicted"], "after": a["core_predicted"],
                "fold": fold,
            })
    moved.sort(key=lambda r: -r["fold"])

    directions = [r["direction"] for r in rows if r["direction"]]
    return {
        "id": "H7_cutover_exam_delta",
        "available": True,
        "role": (
            "SEEN_DIAGNOSTIC under prereg sec. 0.1 -- the exam artefact prints "
            "every measured value and this builder saw it while locating the "
            "seven refused bundles. These rows are reported and MAY NEVER GATE."
        ),
        "prereg": (
            "at least 5 of the 7 rows become ANSWERED; HMF OVER-predicted; "
            "DMHF under-predicted by 1-3 decades; the two Schibilsky pH arms "
            "predict the SAME HMF because the module carries no pH term; no "
            "previously-answered row moves by more than 1.5x"
        ),
        "rows": rows,
        "n_answered": sum(1 for r in rows if r["after_state"] != "out_of_envelope"),
        "n_rows": len(rows),
        "n_within_band": sum(1 for r in rows if r["core_within_band"]),
        "n_within_interval": sum(
            1 for r in rows if r["measured_within_interval"]),
        "PREREG_DIRECTION_MISS": (
            "★ THE PRE-REGISTRATION PREDICTED HMF **OVER**-PREDICTION AND EVERY "
            "ROW UNDER-PREDICTS. The reasoning was that the module has no "
            "validated HMF sink at cooking temperature (K5a G2: the 50-150 C "
            "window is empty), so a source-only node must over-shoot. It does "
            "not, and the missing sink is therefore NOT the binding constraint "
            "-- THE SOURCE TERMS ARE. The two formation limbs are ingested from "
            "an AMINE-FREE, FREEZE-DRIED, essentially anhydrous glucose melt "
            "and are being asked to run in aqueous and food matrices; the melt "
            "-> matrix transfer loses more flux than the absent sink adds back. "
            "That is a quantified, directional result about a transfer nobody "
            "had tested, and it is the opposite of what this wave expected."
            if all(d == "UNDER" for d in directions) and directions else
            "directions are mixed; see the per-row table."
        ),
        "ph_pair_structural_miss": {
            "what": (
                "The two Schibilsky bundles differ ONLY in pH (5.0 vs 8.0) and "
                "the module carries NO pH term on any furanic edge -- K5a "
                "declared gap G8: six distinct pH values appear across the "
                "seven papers of the cluster and NO SINGLE PAPER VARIES pH, so "
                "a pH term would have to be fitted across labs and matrices at "
                "once, which k3 sec. B.2 forbids at family level."
            ),
            "predicted_identical": True,
            "measured_hmf_ph_effect": _fold(
                *(r["measured"] for r in rows
                  if r["compound"].lower().startswith("5-hydroxymethyl")
                  and "Schibi" in r["benchmark_id"])
            ),
            "measured_dmhf_ph_effect": _fold(
                *(r["measured"] for r in rows
                  if r["compound"].lower() == "dmhf")
            ),
            "reading": (
                "The size of the measured pH effect IS the size of the gap, "
                "now measured for the first time on both compounds at once."
            ),
        },
        "trunk_perturbation": {
            "what": (
                "The furanic block hangs on the TRUNK, so four of its steps put "
                "a new drain on B1-fitted species and every previously-answered "
                "row moves a little. B1 is NOT refit -- this wave has no "
                "licence to move its constants -- so the movement is reported "
                "rather than absorbed."
            ),
            "n_rows_moved": len(moved),
            "max_fold": moved[0]["fold"] if moved else 1.0,
            "prereg_threshold_fold": 1.5,
            "within_prereg": bool(not moved or moved[0]["fold"] <= 1.5),
            "largest": moved[:6],
        },
    }


def build_report() -> Dict[str, Any]:
    tests = [
        h1_kocadagli_nacl(),
        h2_gursul_27C(),
        h3_hamzalioglu_matrix_pair(),
        h4_sink_vs_net_ph_pair(),
        h5_apriyantono_ph_trajectory(),
        h6_norfuraneol_dominates(),
        h7_exam_delta(),
    ]
    return {
        "wave": "B7",
        "artifact": "kinetic_core_b7_holdout_report",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "pre_registration": data_paths.rel(data_paths.VALIDATION_DIR / "kinetic_core_b7_prereg.md"),
        "pass_band_fold": PASS_BAND_FOLD,
        "exam": (
            "The seventh test -- the cutover exam's 5 HMF + 2 DMHF bundles -- "
            "is scored by scripts/generators/generate_cutover_final_exam.py and "
            "reported beside this panel. Under the Amendment 9 clause 1 / "
            "Amendment 10 clause 1 precedent those seven rows are "
            "SEEN_DIAGNOSTIC (prereg sec. 0.1) and may never gate."
        ),
        "tests": tests,
    }


def render_markdown(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    add = lines.append
    add("# Build Wave B7 — the furanic channel — HOLD-OUT REPORT")
    add("")
    add(f"Generated {payload['generated_on']} at `{payload['git']['short']}` "
        f"on `{payload['git']['branch']}`. Pass band "
        f"{payload['pass_band_fold']:g}×, taken unchanged from the other module "
        f"scorecards.")
    add("")
    add(f"Pre-registered in `{payload['pre_registration']}`. "
        f"{payload['exam']}")
    add("")
    add("## Summary")
    add("")
    add("| # | hold-out | pre-registered | outcome |")
    add("|---|---|---|---|")
    t = {x["id"]: x for x in payload["tests"]}
    h1 = t["H1_kocadagli_nacl_arm"]
    add(f"| H1 | Kocadagli NaCl arm | 0/3 rate ratios, 2/3 conversions in band "
        f"| **{h1['rate_ratio_within_band']}/3 and "
        f"{h1['conversion_within_band']}/3** |")
    h2 = t["H2_gursul_aktag_2020_27C_zero_accumulation"]
    add(f"| H2 | Gursul 27 °C zero accumulation | PASS (< 100 µg/L) | "
        f"**{'PASS' if h2['passes'] else 'FAIL'}** — "
        f"{h2['predicted_27C_ug_per_l']:.3g} µg/L (prereg MISS, disclosed) |")
    h3 = t["H3_hamzalioglu_matrix_vs_water_selectivity"]
    add(f"| H3 | Hamzalioglu matrix-vs-water | FAIL by construction | "
        f"**FAIL**, {h3['collapse_fold_error']:.2f}× on the collapse; the "
        f"cysteine sub-test passes at "
        f"{h3['sub_test_cysteine_alone']['fold_error']:.2f}× |")
    h4 = t["H4_shu1988_x_wang2008_paired_sink_and_net_pH"]
    add(f"| H4 | shu1988 × wang2008 paired | sink REFUSED, formation right "
        f"sign wrong reason | **sink REFUSED, formation FLAT** — a "
        f"pre-registration correction, disclosed |")
    h5 = t["apriyantono1993_xylose_lysine_pH_trajectory_pair"]
    scored = [c for c in h5["channels"] if c["scored"]]
    add(f"| H5 | apriyantono held-vs-drifting | 2 of 4 channels scored, "
        f"direction pass | **{len(scored)} of {len(h5['channels'])} scored, "
        f"BOTH RATIOS 1.000 — the drifting arm does not drift** |")
    h6 = t["H6_norfuraneol_much_greater_than_DMHF"]
    add(f"| H6 | norfuraneol ≫ DMHF | PASS | "
        f"**{'PASS' if h6['passes'] else 'FAIL'}** |")
    h7 = t.get("H7_cutover_exam_delta", {})
    if h7.get("available"):
        add(f"| H7 | cutover exam, 5 HMF + 2 DMHF | ≥5 of 7 answered; HMF over, "
            f"DMHF 1–3 decades under | **{h7['n_answered']} of {h7['n_rows']} "
            f"answered, {h7['n_within_band']} in band, "
            f"{h7['n_within_interval']} inside the declared interval — HMF "
            f"direction MISSED** |")
    add("")

    for test in payload["tests"]:
        add(f"## {test['id']}")
        add("")
        add(f"**Role.** {test['role']}")
        add("")
        add(f"**Pre-registered.** {test['prereg']}")
        add("")
        if "PREREG_CORRECTION_C1" in test:
            add(f"> {test['PREREG_CORRECTION_C1']}")
            add("")
        if test["id"] == "H1_kocadagli_nacl_arm":
            add(test["prediction_basis"])
            add("")
            add("| quantity | T (°C) | measured | predicted | fold | in band |")
            add("|---|---:|---:|---:|---:|---|")
            for r in test["rows"]:
                add(f"| {r['quantity']} | {r['temperature_C']:.0f} | "
                    f"{r['measured']:.2f} | {r['predicted']:.2f} | "
                    f"{r['fold_error']:.2f}× | "
                    f"{'yes' if r['within_band'] else '**no**'} |")
            add("")
            add(f"> {test['what_it_measures']}")
        elif test["id"] == "H2_gursul_aktag_2020_27C_zero_accumulation":
            add("| T (°C) | predicted HMF (µg/L) | fructose-limb share | "
                "3-DG-limb share |")
            add("|---:|---:|---:|---:|")
            for r in test["rows"]:
                add(f"| {r['temperature_C']:.0f} | "
                    f"{r['predicted_ug_per_l']:.4g} | "
                    f"{r['fructose_limb_share']:.4f} | "
                    f"{r['three_deoxyglucosone_limb_share']:.4f} |")
            add("")
            add(f"> {test['PREREG_MISS_DISCLOSED']}")
            add("")
            add(f"Declared floor **{test['declared_floor_ug_per_l']:g} µg/L** — "
                f"{test['floor_is_this_modules_own']}")
            add("")
            add(f"> {test['charge_is_declared_not_transcribed']}")
            add("")
            add(f"> {test['companion_37C_note']}")
        elif test["id"] == "H3_hamzalioglu_matrix_vs_water_selectivity":
            add(f"Measured Cys/Lys selectivity: **"
                f"{test['measured_selectivity_Cys_over_Lys_water']:.1f}×** in "
                f"water, **{test['measured_selectivity_Cys_over_Lys_coffee']:.2f}×** "
                f"in roasted coffee — a "
                f"**{test['measured_collapse_factor']:.1f}× collapse**. "
                f"Predicted collapse **1.000**, fold error "
                f"**{test['collapse_fold_error']:.2f}×**.")
            add("")
            add(f"> {test['sub_test_cysteine_alone']['note']}")
            add("")
            add(f"**Why the model cannot form the ratio.** "
                f"{test['why_the_model_cannot_form_the_ratio']}")
        elif test["id"] == "H4_shu1988_x_wang2008_paired_sink_and_net_pH":
            add("### Sink half — REFUSED")
            add("")
            add(test["sink_half"]["why"])
            add("")
            add(f"> {test['sink_half']['and_the_measurement_itself_is_ambiguous']}")
            add("")
            add("### Formation half — scored")
            add("")
            add("| pH | measured (mg/mol MG) | predicted (mg/mol MG) |")
            add("|---:|---:|---:|")
            for r in test["formation_half"]["rows"]:
                add(f"| {r['ph']:g} | {r['measured_mg_per_mol_MG']:.0f} | "
                    f"{r['predicted_mg_per_mol_MG']:.4g} |")
            add("")
            add(f"Measured pH 3→8 trend **"
                f"{test['formation_half']['measured_pH3_to_pH8_trend']:.2f}×**, "
                f"predicted **"
                f"{test['formation_half']['predicted_pH3_to_pH8_trend']:.3f}×**.")
            add("")
            add(f"> {test['formation_half']['verdict']}")
            add("")
            add(f"**Paired verdict.** {test['paired_verdict']}")
        elif test["id"] == "apriyantono1993_xylose_lysine_pH_trajectory_pair":
            add(f"**Declared limitation.** {test['declared_limitation']}")
            add("")
            add(f"**And a second.** {test['second_declared_limitation']}")
            add("")
            add("| channel | measured drift/held | predicted | fold | scored |")
            add("|---|---|---:|---:|---|")
            for c in test["channels"]:
                measured = c.get("measured_drift_over_held")
                measured_s = (
                    f"{measured:.4g}×" if measured is not None
                    else c.get("measured_direction", "—")
                )
                predicted = c.get("predicted")
                predicted_s = f"{predicted:.4g}" if predicted else "—"
                fold = c.get("fold_error")
                fold_s = f"{fold:.3g}×" if fold else "—"
                add(f"| {c['channel']} | {measured_s} | {predicted_s} | "
                    f"{fold_s} | {'yes' if c['scored'] else 'REFUSED'} |")
            add("")
            for c in test["channels"]:
                if not c["scored"]:
                    add(f"* **{c['channel']}** — refused: {c['why_not']}")
            add("")
            ph = test["ph_state_channel"]
            add("### The pH-state channel — this is what exams B2.2/B2.3")
            add("")
            add(f"Measured: held at 5.0, drifting **4.9 → 2.6** "
                f"({ph['measured_drop_ph_units']:.1f} pH units). "
                f"Predicted drifting endpoint: "
                f"**{ph['predicted_drifting_final']}**.")
            add("")
            add(f"> {ph['what_this_exams']}")
            add("")
            add(f"> {test['not_scored_against_H6']}")
            add("")
            add(f"> {test['VERDICT']}")
            add("")
            add(f"**Pre-registered outcome.** {test['prereg_outcome']}")
            add("")
            add(f"**Record correction.** "
                f"{test['record_correction_this_wave_carries_forward']}")
        elif test["id"] == "H7_cutover_exam_delta":
            if not test.get("available"):
                add(test["why"])
            else:
                add("| bundle | compound | T (°C) | measured (ppb) | "
                    "predicted (ppb) | fold | dir | band | in interval |")
                add("|---|---|---:|---:|---:|---:|---|---|---|")
                for r in test["rows"]:
                    add(f"| `{r['benchmark_id']}` | {r['compound']} | "
                        f"{r['temp_C']} | {r['measured']:.4g} | "
                        f"{r['core_predicted']:.4g} | "
                        f"{r['core_fold_error']:.2f}× | {r['direction']} | "
                        f"{'**PASS**' if r['core_within_band'] else 'fail'} | "
                        f"{'yes' if r['measured_within_interval'] else 'no'} |")
                add("")
                add("Every one of the seven was `out_of_envelope` before this "
                    "wave.")
                add("")
                add(f"> {test['PREREG_DIRECTION_MISS']}")
                add("")
                ph = test["ph_pair_structural_miss"]
                add("### The pH pair — a pre-registered structural miss")
                add("")
                add(ph["what"])
                add("")
                add(f"Predicted HMF and DMHF are **identical** at pH 5.0 and "
                    f"pH 8.0. Measured pH effect: **"
                    f"{ph['measured_hmf_ph_effect']:.2f}×** on HMF and **"
                    f"{ph['measured_dmhf_ph_effect']:.2f}×** on DMHF. "
                    f"{ph['reading']}")
                add("")
                tp = test["trunk_perturbation"]
                add("### What the trunk perturbation cost")
                add("")
                add(tp["what"])
                add("")
                add(f"**{tp['n_rows_moved']}** previously-answered rows moved; "
                    f"largest **{tp['max_fold']:.4f}×** against a "
                    f"pre-registered ceiling of {tp['prereg_threshold_fold']}× "
                    f"— **{'within' if tp['within_prereg'] else 'OUTSIDE'}**.")
                add("")
                add("| bundle | compound | before | after | fold |")
                add("|---|---|---:|---:|---:|")
                for r in tp["largest"]:
                    add(f"| `{r['benchmark_id']}` | {r['compound']} | "
                        f"{r['before']:.4g} | {r['after']:.4g} | "
                        f"{r['fold']:.4f}× |")
        elif test["id"] == "H6_norfuraneol_much_greater_than_DMHF":
            add(f"Norfuraneol **{test['predicted_NF_mmol_per_l']:.4g}** mmol/L "
                f"vs DMHF **{test['predicted_DMHF_mmol_per_l']:.4g}** mmol/L "
                f"at Blank's own conditions — "
                f"**{'PASS' if test['passes'] else 'FAIL'}**.")
            add("")
            add(f"> {test['why_it_can_only_ever_be_ordinal']}")
            add("")
            add(f"> {test['note_on_why_this_is_not_circular']}")
        add("")
    return "\n".join(lines) + "\n"


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))
    args = parser.parse_args(argv)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = build_report()
    (out / f"{BASENAME}.json").write_text(json.dumps(payload, indent=2) + "\n")
    (out / f"{BASENAME}.md").write_text(render_markdown(payload))
    print(f"wrote {out / BASENAME}.json and .md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

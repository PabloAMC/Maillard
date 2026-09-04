#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b3_holdout.py

THE PRE-REGISTERED HOLD-OUT SCORING OF BUILD WAVE B3 (acrylamide / safety).

Reads the FROZEN parameter set from
`results/validation/kinetic_core_b3_fit_report.json`, changes NOTHING, predicts
every declared Module 3 hold-out row, and writes
`results/validation/kinetic_core_b3_holdout_report.{json,md}`.

===========================================================================
THE RULE THIS SCRIPT OBEYS
===========================================================================
No parameter is touched here. There is no optimiser in this file, no bounds and
no `least_squares` import. It integrates the FROZEN network at conditions the
fit never saw and prints the fold error, row by row, pass or fail, with no
averaging away of failures and no per-row weighting that could hide one.

===========================================================================
PRE-REGISTERED FAILURES AND UNSCOREABLES -- stated BEFORE the numbers
===========================================================================
Written down here so that they cannot later be presented as discoveries.

  1. **De Vleeschouwer II's GLUTAMINE system WILL FAIL, in a known direction.**
     Every competition channel in this module can only REMOVE precursor or
     REMOVE acrylamide, so a competitor can never raise the yield above the
     control's. Glutamine is measured to RAISE it. The model will predict a
     ratio at or below 1 against a measured promotion. No term was added,
     deliberately: the promotion's TEMPERATURE SHAPE is the B5.5 sign-crossing
     and its a_w-0.92 half is this very hold-out, so a term fitted to the
     Claeys half would be a term built toward the hold-out. See
     `parameters_acrylamide.DELIBERATE_UNDERFITS`.

  2. **De Vleeschouwer I's SUCROSE system is UNSCOREABLE.** The module has no
     sucrose species, because the only measurement of sucrose hydrolysis in the
     corpus sits in the hold-out column. Inventing one in order to be able to
     score the hold-out would defeat the hold-out. Scored `unscoreable`, not
     `fail`, because the two mean different things.

  3. **De Vleeschouwer I's FRUCTOSE system is scored DIAGNOSTIC, not gating,
     and its prediction is reported as an INTERVAL.** Two independent reasons,
     both stated in the declaration itself. (a) The source's own 95 % HPD on
     the fructose formation constant SPANS ZERO, so the target is not a
     determination -- the declaration's own words are that the model "should
     reproduce [it] as *low confidence*, not as a fitted value". (b) The model
     has NO fructose-specific parameter: fructose reaches acrylamide only by
     isomerising to glucose through Martins' measured constants, which were
     measured at 80-120 C and are being evaluated at 160-200 C. The prediction
     therefore comes with a Monte-Carlo interval over the published parameter
     uncertainties, and the interval is the answer.

  4. **Knol 2010's Table 2 step constants are UNSCOREABLE**, and this is a
     corpus fact rather than a model limitation: the seven steps are not
     transcribed anywhere in the extraction dossiers. Only Knol 2010's
     end-to-end acrylamide YIELD is available, and that is scored.

  5. **Knol 2009's real-food band is UNSCOREABLE.** No precursor loading for
     its matrix is transcribed in the corpus, so predicting into the band would
     require inventing an asparagine and glucose content. The source's author
     refused to transfer his own model to real food (inventory sec. B8.5),
     which is exactly what this row was held out to test -- and the honest
     result of that test is that the transfer cannot even be attempted from the
     evidence the repository holds.

  6. **The one row that could go either way is Knol 2010's ~3 % molar yield.**
     It is GATING. The model's yield scales with the sugar concentration
     through the bimolecular initiation, so predicting a 0.1 M system's yield
     from constants identified on a 0.01 M one is a genuine test of the
     concentration handling that Z1 sec. E diagnosed as the shipped model's
     central defect.

===========================================================================
A DISCLOSURE ABOUT FIREWALL HYGIENE
===========================================================================
The build brief directed this wave to read
`data/lit/extraction_dossiers/k3_final_parameter_inventory.md` sec. A.2 and
`k1_kinetic_parameters.md` sec. 2, and BOTH print the declared HOLD-OUT columns
in the same tables as the FIT ones: De Vleeschouwer 2009 Part I's fructose and
sucrose columns, and Part II's glutamine column. Inventory sec. B5.5
additionally prints the glutamine promotion percentages. Those values were
therefore SEEN before the fit was run. They were not used: they appear in no
fit row, no bound, no initialisation and in no file under `src/kinetic_core/`,
and `tests/unit/test_kinetic_core_b3.py` enforces that mechanically with a
literal grep over executable code.

The frozen `mp_holdout_*` bundles under `data/benchmarks/external_validation/`
were NOT opened by this wave at any point.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import operative_acrylamide_parameters  # noqa: E402
from src.kinetic_core.acrylamide import (  # noqa: E402
    apparent_activation_energy,
    apparent_lumped_constants,
)

CELSIUS = 273.15
FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b3_fit_report.json"

B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

#: Pass bands, declared BEFORE the numbers are seen, and identical to B2's.
#: A LEVEL row passes within 3x; a RATIO row within 2x. A DIAGNOSTIC row is
#: scored and reported but does not gate.
PASS_BAND_LEVEL = 3.0
PASS_BAND_RATIO = 2.0

#: Monte-Carlo draws for the interval on the low-confidence rows.
N_DRAWS = 200
MC_SEED = 20260828

DV_LOADING = 3000.0
KNOL_LOADING = 100.0
WINDOW = 45.0


def load_frozen() -> Dict[str, Any]:
    if not FIT_REPORT.exists():
        raise SystemExit(
            f"{FIT_REPORT} not found. Run generate_kinetic_core_b3_fit.py "
            f"first; the hold-out scorer never fits anything itself."
        )
    payload = json.loads(FIT_REPORT.read_text())
    frozen = payload["frozen_parameters"]
    parameters = operative_acrylamide_parameters(
        B1_FITTED, frozen["log10_k_ref_at_160C"], frozen["fitted_Ea_kJ_mol"]
    )
    return {"parameters": parameters, "fit_payload": payload}


def observe(parameters, initial, temperature_c, aw, minutes=WINDOW, n_points=61):
    return apparent_lumped_constants(
        parameters, float(temperature_c) + CELSIUS, initial, float(minutes),
        n_points=n_points, water_activity=aw, rtol=1e-8, atol=1e-14,
    )


def score(observed: float, predicted: float, band: float) -> Dict[str, Any]:
    if not np.isfinite(predicted) or predicted <= 0 or observed <= 0:
        return {"observed": observed, "predicted": float(predicted),
                "fold_error": float("nan"), "pass_band": band, "pass": False}
    fold = predicted / observed
    return {
        "observed": observed,
        "predicted": float(predicted),
        "fold_error": fold,
        "pass_band": band,
        "pass": bool((1.0 / band) <= fold <= band),
    }


# ---------------------------------------------------------------------------
# The interval machinery -- how "low confidence" is expressed as a number
# ---------------------------------------------------------------------------
# Only the PUBLISHED uncertainties of the MEASURED parameters are propagated:
# the seven acrylamide constants of `MEASURED_ACRYLAMIDE` and the two Martins
# isomerisation barriers that carry fructose into the glucose lane. The FITTED
# parameters' Gauss-Newton intervals are deliberately NOT propagated, because
# they are conditional on the model being right and would make the band look
# more informative than it is. The band is therefore a LOWER bound on the true
# uncertainty, and saying so is the point.

_PERTURBED = (
    "k_asn_glc", "k_int1_acr", "k_asn_asp", "k_acr_cys", "k_cys_sink",
    "k_cys_glc", "k_asp_sink", "k_glc_fru", "k_fru_glc",
)


def draw_parameters(parameters, rng):
    """One Monte-Carlo draw over the PUBLISHED parameter uncertainties."""
    out = dict(parameters)
    for key in _PERTURBED:
        parameter = out.get(key)
        if parameter is None:
            continue
        k_ref = float(parameter.k_ref)
        ea = float(parameter.ea_kj_mol)
        k_ci = getattr(parameter, "k_ref_ci95", None)
        ea_ci = getattr(parameter, "ea_ci95_kj_mol", None)
        if k_ci:
            # A normal draw truncated at a positive floor: several of these
            # intervals reach zero or below, which is a statement about the
            # measurement and must not become a negative rate constant.
            k_ref = max(1e-12 * k_ref, rng.normal(k_ref, float(k_ci) / 1.96))
        if ea_ci:
            ea = rng.normal(ea, float(ea_ci) / 1.96)
        out[key] = replace(parameter, k_ref=k_ref, ea_kj_mol=ea)
    return out


def interval(function, parameters, n_draws=N_DRAWS, seed=MC_SEED):
    """Median and 95 % band of ``function(parameters)`` over the draws."""
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_draws):
        try:
            value = float(function(draw_parameters(parameters, rng)))
        except Exception:
            value = float("nan")
        if np.isfinite(value):
            values.append(value)
    if not values:
        return {"median": float("nan"), "low": float("nan"),
                "high": float("nan"), "n_draws": 0}
    array = np.array(values)
    return {
        "median": float(np.median(array)),
        "low": float(np.percentile(array, 2.5)),
        "high": float(np.percentile(array, 97.5)),
        "spread_factor": float(
            np.percentile(array, 97.5) / max(np.percentile(array, 2.5), 1e-30)
        ),
        "n_draws": int(array.size),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Score Build Wave B3's (acrylamide / safety) pre-registered "
            "hold-out rows against the frozen B3 parameters, fitting nothing; "
            "writes "
            "results/validation/kinetic_core_b3_holdout_report.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(data_paths.VALIDATION_DIR),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    frozen = load_frozen()
    p = frozen["parameters"]
    rows: List[Dict[str, Any]] = []

    def add(row_id, group, gating, anchor, result, comment="", extra=None):
        entry = {"id": row_id, "group": group, "gating": gating,
                 "source_anchor": anchor, "comment": comment}
        entry.update(result)
        if extra:
            entry.update(extra)
        rows.append(entry)

    # =====================================================================
    # 1. DE VLEESCHOUWER I, THE FRUCTOSE SYSTEM -- DIAGNOSTIC, WITH AN
    #    INTERVAL, BECAUSE THE SOURCE'S OWN HPD SPANS ZERO
    # =====================================================================
    fructose_initial = {"Asn": DV_LOADING, "Fru": DV_LOADING}
    glucose_initial = {"Asn": DV_LOADING, "Glc": DV_LOADING}

    def _ratio(params):
        f = observe(params, fructose_initial, 160.0, 0.92)["peak_acrylamide_ppb"]
        g = observe(params, glucose_initial, 160.0, 0.92)["peak_acrylamide_ppb"]
        return f / g if g > 0 else float("nan")

    ratio_point = _ratio(p)
    ratio_band = interval(_ratio, p)
    # The target is the ratio of the HELD-OUT fructose formation constant to
    # the FIT glucose one, 7.40e-3 / 3.57e-3.
    target_ratio = 7.40e-3 / 3.57e-3
    add(
        "dv1_fructose_over_glucose_acrylamide",
        "De Vleeschouwer 2009 I -- the sugar-transfer test",
        "DIAGNOSTIC",
        "De Vleeschouwer 2009 I Table 3, fructose k_Fref 7.40 +/- 9.48 e-3 "
        "against the FIT glucose 3.57 +/- 1.38 e-3 min^-1",
        score(target_ratio, ratio_point, PASS_BAND_RATIO),
        comment=(
            "LOW CONFIDENCE, and that is the required answer rather than a "
            "hedge. The source's own 95 % HPD on the fructose constant is "
            "+/- 9.48 against an estimate of 7.40, i.e. IT SPANS ZERO: the "
            "target is not a determination and a model that matched it "
            "precisely would be matching noise. The declaration says so "
            "explicitly ('fructose's HPDs span zero, which the model should "
            "reproduce as LOW CONFIDENCE, not as a fitted value'). The model "
            "carries NO fructose-specific parameter -- fructose reaches "
            "acrylamide only by isomerising to glucose through Martins' "
            "measured constants, evaluated 40-80 C above the window they were "
            "measured in. The interval below, not the point value, is the "
            "prediction."
        ),
        extra={
            "confidence": "LOW",
            "prediction_interval_95": ratio_band,
            "target_own_interval_spans_zero": True,
            "fitted_parameters_for_this_sugar": 0,
        },
    )

    def _fructose_level(params):
        return observe(params, fructose_initial, 160.0, 0.92)["peak_acrylamide_ppb"]

    add(
        "dv1_fructose_peak_acrylamide_ppb",
        "De Vleeschouwer 2009 I -- the sugar-transfer test",
        "DIAGNOSTIC -- REPORTED, NOT SCORED",
        "no absolute acrylamide concentration is printed for the De "
        "Vleeschouwer systems; only rate constants",
        {"observed": None, "predicted": _fructose_level(p),
         "fold_error": float("nan"), "pass_band": None, "pass": None},
        comment=(
            "There is nothing to score this against: the source publishes "
            "constants, not concentrations. It is reported so that the "
            "interval width on an absolute level is visible next to the "
            "interval on the ratio."
        ),
        extra={"confidence": "LOW",
               "prediction_interval_95": interval(_fructose_level, p)},
    )

    # =====================================================================
    # 2. DE VLEESCHOUWER I, THE SUCROSE SYSTEM -- UNSCOREABLE
    # =====================================================================
    add(
        "dv1_sucrose_system",
        "De Vleeschouwer 2009 I -- the sugar-transfer test",
        "UNSCOREABLE",
        "De Vleeschouwer 2009 I Table 3, sucrose column",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        comment=(
            "The module has NO SUCROSE SPECIES. Sucrose enters the acrylamide "
            "lane only by hydrolysing to glucose and fructose, and the only "
            "measurement of that hydrolysis rate in the corpus is in this same "
            "held-out column. Adding a sucrose species would therefore have "
            "meant either inventing its hydrolysis rate or reading the "
            "hold-out. Reported UNSCOREABLE rather than failed, because the "
            "two are different findings: this one says the module's scope is "
            "narrower than the declaration's, not that its chemistry is wrong."
        ),
    )

    # =====================================================================
    # 3. DE VLEESCHOUWER II, THE GLUTAMINE SYSTEM -- PRE-REGISTERED FAIL
    # =====================================================================
    gln_initial = {"Asn": DV_LOADING, "Glc": DV_LOADING, "Gln": DV_LOADING}
    for temperature_c, promotion in ((120.0, 2.67), (200.0, 1.20)):
        gln = observe(p, gln_initial, temperature_c, 0.92)["peak_acrylamide_ppb"]
        control = observe(p, glucose_initial, temperature_c, 0.92)["peak_acrylamide_ppb"]
        add(
            f"dv2_glutamine_promotion_{int(temperature_c)}C",
            "De Vleeschouwer 2009 II -- the glutamine system (B5.5)",
            "GATING",
            f"De Vleeschouwer 2009 II, glutamine promotion of acrylamide at "
            f"a_w 0.92, {int(temperature_c)} C",
            score(promotion, gln / control if control > 0 else float("nan"),
                  PASS_BAND_RATIO),
            comment=(
                "PRE-REGISTERED FAILURE, direction and cause both known before "
                "the fit. Competition in this module can only SUPPRESS: a "
                "competitor either eats the shared glucose or scavenges the "
                "acrylamide, so the predicted ratio cannot exceed 1. Glutamine "
                "is measured to PROMOTE. The model was deliberately not given "
                "a promotion mechanism, because the promotion's temperature "
                "SHAPE (growing with T in liquid, shrinking with T at a_w "
                "0.92) is the B5.5 sign-crossing whose a_w-0.92 half is THIS "
                "ROW. A term fitted to the Claeys half would have been a term "
                "built toward this hold-out. The failure localises a missing "
                "MECHANISM, which is what the declaration says this row is for."
            ),
        )

    gln_kf = observe(p, gln_initial, 160.0, 0.92)["k_F_app_per_min"]
    control_kf = observe(p, glucose_initial, 160.0, 0.92)["k_F_app_per_min"]
    add(
        "dv2_glutamine_kF_ratio_160C",
        "De Vleeschouwer 2009 II -- the glutamine system (B5.5)",
        "GATING",
        "De Vleeschouwer 2009 II Table 3, glutamine k_Fref 8.05e-3 against the "
        "control's 3.57e-3 min^-1 (a 2.25x promotion)",
        score(8.05e-3 / 3.57e-3, gln_kf / control_kf if control_kf > 0 else float("nan"),
              PASS_BAND_RATIO),
        comment=(
            "PRE-REGISTERED FAILURE, the same one, on the formation constant "
            "rather than on the yield. The model's ratio can approach 1 and "
            "cannot exceed it, because suppression is the only thing a "
            "competitor in this module can do."
        ),
    )

    # =====================================================================
    # 4. KNOL 2010 -- THE CROSS-LAB YIELD, THE ONE GENUINELY OPEN ROW
    # =====================================================================
    knol = observe(p, {"Asn": KNOL_LOADING, "Glc": KNOL_LOADING}, 160.0, 1.0)
    add(
        "knol2010_molar_yield_on_asparagine",
        "Knol 2010 -- the cross-lab trunk test",
        "GATING",
        "Knol 2010: acrylamide yield ~3 % of the initial asparagine in an "
        "aqueous 0.1 M asparagine/glucose system",
        score(0.03, knol["molar_yield_on_asparagine"], PASS_BAND_LEVEL),
        comment=(
            "THE ROW TO READ. It is a THIRD lab, and -- more importantly -- it "
            "is at a TEN TIMES higher precursor loading than the Claeys system "
            "the fitted partition was identified on. Because initiation here "
            "is genuinely bimolecular, the predicted molar yield RISES with "
            "loading; a model with a fixed yield fraction (which is what the "
            "shipped lane had) cannot move on this axis at all. Z1 sec. E "
            "diagnosed exactly that saturation as the central defect, and this "
            "row is the first out-of-sample test of the fix."
        ),
        extra={"predicted_peak_ppb": knol["peak_acrylamide_ppb"],
               "predicted_t_max_min": knol["time_of_peak_min"],
               "loading_mmol_L": KNOL_LOADING,
               "loading_is_an_inference_from": "Knol 2010 prose, 'at 0.1 M'"},
    )

    knol_low = observe(p, {"Asn": KNOL_LOADING, "Glc": KNOL_LOADING}, 120.0, 1.0)
    knol_high = observe(p, {"Asn": KNOL_LOADING, "Glc": KNOL_LOADING}, 200.0, 1.0)
    ea_formation = apparent_activation_energy(
        knol_low["k_F_app_per_min"], 120.0 + CELSIUS,
        knol_high["k_F_app_per_min"], 200.0 + CELSIUS,
    )
    add(
        "knol2010_formation_barrier_ceiling",
        "Knol 2010 -- the cross-lab trunk test",
        "DIAGNOSTIC (one-sided)",
        "Knol 2010's largest published activation energy is 93 +/- 12 kJ/mol",
        {"observed": 93.0, "predicted": float(ea_formation),
         "fold_error": float("nan"), "pass_band": None,
         "pass": bool(ea_formation <= 93.0 + 12.0)},
        comment=(
            "A ONE-SIDED CEILING, not a level. What is known about Knol 2010 "
            "from the corpus is that no barrier in it exceeds 93 +/- 12; the "
            "seven step constants themselves are not transcribed. The model's "
            "apparent formation barrier is compared against that ceiling. "
            "Note the tension this row inherits: Claeys measures 168.25 +/- "
            "3.80 for the same transformation and Claeys is a FIT row, so the "
            "model is being pulled in two directions by two datasets that "
            "disagree by 75 kJ/mol."
        ),
    )

    add(
        "knol2010_table2_seven_steps",
        "Knol 2010 -- the cross-lab trunk test",
        "UNSCOREABLE",
        "Knol 2010 Table 2, 7 steps x 5 temperatures + Ea + HPD",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        comment=(
            "NOT TRANSCRIBED ANYWHERE IN THE CORPUS. The extraction dossiers "
            "carry only Knol 2010's three organic-acid and isomerisation "
            "barriers, which orchestrator decision 2 assigns to Module 4 as "
            "FIT. This is a CORPUS gap, not a model limitation, and the "
            "declaration's strongest acrylamide hold-out is therefore only "
            "partly available: what can be scored is the end-to-end yield row "
            "above."
        ),
    )

    # =====================================================================
    # 5. KNOL 2009 REAL-FOOD BAND -- UNSCOREABLE
    # =====================================================================
    add(
        "knol2009_real_food_band",
        "Knol 2009 -- real food",
        "UNSCOREABLE",
        "Knol 2009, 9.3e3 - 2.6e4 ug/kg dm in real food",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        comment=(
            "No asparagine or glucose content for Knol's matrix is transcribed "
            "in the corpus, so a prediction into this band would begin by "
            "inventing the precursor loading it is most sensitive to. The "
            "source's own author refused to transfer his model to real food -- "
            "'the parameters are only applicable for specific experimental "
            "conditions such as time-temperature profile of frying, potato "
            "genotype, slice thickness and initial concentration of "
            "precursors' (inventory sec. B8.5). Reporting UNSCOREABLE is the "
            "honest outcome of a test the evidence cannot support, and it is "
            "NOT a pass."
        ),
    )

    # =====================================================================
    # Scorecard
    # =====================================================================
    gating = [r for r in rows if r["gating"] == "GATING"]
    gating_passed = [r for r in gating if r["pass"]]
    unscoreable = [r for r in rows if r["gating"] == "UNSCOREABLE"]
    diagnostic = [r for r in rows if r["gating"].startswith("DIAGNOSTIC")]

    payload: Dict[str, Any] = {
        "wave": "B3 -- acrylamide / safety, HOLD-OUT",
        "generated_by": "scripts/generators/generate_kinetic_core_b3_holdout.py",
        "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md D.5 (Module 3)",
        "frozen_from": str(FIT_REPORT.relative_to(ROOT)),
        "frozen_parameters": frozen["fit_payload"]["frozen_parameters"],
        "no_parameter_was_changed": True,
        "pass_bands": {"level": PASS_BAND_LEVEL, "ratio": PASS_BAND_RATIO},
        "monte_carlo": {
            "n_draws": N_DRAWS,
            "seed": MC_SEED,
            "what_is_propagated": (
                "the PUBLISHED uncertainties of the seven measured acrylamide "
                "constants and of the two Martins isomerisation barriers that "
                "carry fructose into the glucose lane"
            ),
            "what_is_NOT_propagated": (
                "the fitted parameters' Gauss-Newton intervals, because they "
                "are conditional on the model being right. Every interval "
                "printed here is therefore a LOWER BOUND on the true "
                "uncertainty."
            ),
        },
        "scorecard": {
            "gating_rows": len(gating),
            "gating_passed": len(gating_passed),
            "diagnostic_rows": len(diagnostic),
            "unscoreable_rows": len(unscoreable),
            "pre_registered_failures": [
                r["id"] for r in rows
                if "PRE-REGISTERED FAILURE" in r.get("comment", "")
            ],
        },
        "rows": rows,
    }

    out_json = Path(args.output_dir) / "kinetic_core_b3_holdout_report.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, default=str))

    lines: List[str] = []
    a = lines.append
    a("# Kinetic core, Build Wave B3 -- the acrylamide hold-out scorecard")
    a("")
    a("Frozen parameters from `results/validation/kinetic_core_b3_fit_report.json`. "
      "**Nothing was fitted in this script and there is no optimiser in it.**")
    a("")
    a(f"- gating rows: **{len(gating_passed)} / {len(gating)} pass** "
      f"(level band {PASS_BAND_LEVEL}x, ratio band {PASS_BAND_RATIO}x)")
    a(f"- diagnostic rows: {len(diagnostic)}")
    a(f"- unscoreable rows: {len(unscoreable)} -- and each one names WHY, "
      f"because 'unscoreable' and 'pass' are not the same thing")
    a("")
    pre_registered = payload["scorecard"]["pre_registered_failures"]
    open_rows = [
        r for r in gating if r["id"] not in pre_registered
    ]
    a("## What the scorecard actually says")
    a("")
    a(f"**{len(pre_registered)} of the {len(gating)} gating rows were "
      f"pre-registered as failures BEFORE the fit was run**, with the "
      f"mechanism named: this module has no glutamine promotion route and was "
      f"deliberately not given one, because the promotion's temperature shape "
      f"is the hold-out itself. Those rows failing is the expected outcome and "
      f"is not new information; what they buy is a localisation -- the missing "
      f"thing is a MECHANISM, not a mis-set constant.")
    a("")
    if open_rows:
        plural = "rows are" if len(open_rows) > 1 else "row is"
        a(f"**The genuinely open {plural} "
          f"`{'`, `'.join(r['id'] for r in open_rows)}`.** Result:")
        a("")
        for r in open_rows:
            fold = r.get("fold_error", float("nan"))
            verdict = "PASS" if r["pass"] else "FAIL"
            a(f"- `{r['id']}`: predicted {r['predicted']:.4g} against "
              f"{r['observed']:.4g}, a factor of "
              f"{(1.0 / fold if np.isfinite(fold) and fold < 1 else fold):.2f} "
              f"{'LOW' if np.isfinite(fold) and fold < 1 else 'HIGH'}, "
              f"against a {r['pass_band']}x band -- **{verdict}**.")
        a("")
    a("Read the three UNSCOREABLE rows as what they are. Two of them "
      "(Knol 2010's step table, Knol 2009's real-food band) are CORPUS gaps: "
      "the numbers are not transcribed anywhere in the extraction dossiers, so "
      "the declaration's strongest acrylamide hold-out is only partly "
      "available. The third (sucrose) is a scope limit this module chose "
      "deliberately in order not to read the hold-out. None of the three is a "
      "pass.")
    a("")
    a("## Row by row")
    a("")
    a("| row | gating | observed | predicted | fold | result |")
    a("|---|---|---:|---:|---:|---|")
    for r in rows:
        observed = "--" if r["observed"] is None else f"{r['observed']:.4g}"
        predicted = "--" if r["predicted"] is None else f"{r['predicted']:.4g}"
        fold = (
            f"{r['fold_error']:.2f}x" if np.isfinite(r.get("fold_error", np.nan))
            else "--"
        )
        if r["pass"] is None:
            result = "n/a"
        elif r["pass"]:
            result = "**PASS**"
        else:
            result = "FAIL"
        a(f"| `{r['id']}` | {r['gating']} | {observed} | {predicted} | "
          f"{fold} | {result} |")
    a("")
    a("## Intervals on the low-confidence rows")
    a("")
    for r in rows:
        band = r.get("prediction_interval_95")
        if not band:
            continue
        a(f"- `{r['id']}`: median **{band['median']:.4g}**, 95 % band "
          f"**{band['low']:.4g} - {band['high']:.4g}** "
          f"(a factor of {band.get('spread_factor', float('nan')):.1f}), "
          f"from {band['n_draws']} draws.")
    a("")
    a("Read the two bands together, and note that they are NOT the same kind "
      "of statement. The RATIO band is narrow because most of the parameter "
      "uncertainty cancels between the two sugars -- the two systems share "
      "every constant except the isomerisation that carries fructose into the "
      "glucose lane. The LEVEL band is wider. NEITHER contains the structural "
      "uncertainty, which here is the dominant term: the model has no "
      "fructose-specific chemistry at all, and the isomerisation constants "
      "that do all the work were measured 40-80 C below where they are being "
      "evaluated. A narrow band around a structurally wrong route is a "
      "precise wrong answer, and that is the honest reading of this row.")
    a("")
    a("## What each row means")
    a("")
    for r in rows:
        a(f"### `{r['id']}` -- {r['gating']}")
        a("")
        a(f"Source: {r['source_anchor']}")
        a("")
        a(r["comment"])
        a("")

    out_md = Path(args.output_dir) / "kinetic_core_b3_holdout_report.md"
    out_md.write_text("\n".join(lines))
    print(f"wrote {out_json}")
    print(f"wrote {out_md}")
    print(f"gating: {len(gating_passed)}/{len(gating)} pass; "
          f"{len(unscoreable)} unscoreable; {len(diagnostic)} diagnostic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

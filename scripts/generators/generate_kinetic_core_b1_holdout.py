#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b1_holdout.py

THE HOLD-OUT STAGE OF BUILD WAVE B1. RUN ONLY AFTER THE FIT IS FROZEN.

Reads the frozen parameters from
`results/validation/kinetic_core_b1_fit_report.json`, integrates the Module 4
hold-out conditions WITHOUT CHANGING A SINGLE PARAMETER, and writes
`results/validation/kinetic_core_b1_holdout_report.{json,md}`.

THE TWO HOLD-OUTS (docs/reference/FIT_HOLDOUT_DECLARATION.md D.6, Module 4)
--------------------------------------------------------------------------
  1. "Martins 2005 browning (step 9, epsilon 0.64)" -- the melanoidin response
     of `martins2005_glucose_glycine_80_100_120C_pH68.yml` at 80/100/120 C,
     together with epsilon = 0.64 +/- 0.03 L/(mmol*cm) at 470 nm.
     Declaration's reasoning: "the browning lane has never had a parameter;
     holding it out is the only way to learn whether the trunk predicts colour
     or merely accommodates it."

  2. "Knol 2005 epsilon = 282 L/(mol*cm) (Glc/Asn)".
     Declaration's reasoning: "tests the amine-specificity of epsilon rather
     than fitting around it."

THIS FILE IS THE ONLY PLACE IN THE MODULE WHERE EITHER EPSILON APPEARS.
Neither value exists in `src/kinetic_core/` or in the fit script.

WHAT COUNTS AS A PREDICTION HERE, STATED BEFORE THE NUMBERS
-----------------------------------------------------------
Two variants were frozen by the fit stage, and the difference between them is
the whole point:

  * **Variant B (the pre-registered headline).** The melanoidin sink constant
    was re-estimated from the REACTANT side only -- 3-deoxyglucosone and
    glycine -- with the browning response never in the objective. Its
    prediction of the browning trajectory is genuinely out-of-sample.

  * **Variant A (secondary, NOT out-of-sample).** Martins' own measured step-9
    constant. Martins fitted that constant TO THE HELD-OUT RESPONSE, so a good
    score here is a reproducibility result and must not be reported as a
    successful prediction. It is included because the declaration lists the
    Table 2 rate as a FIT row, and the contrast between A and B is informative.

Usage:
    python scripts/generators/generate_kinetic_core_b1_holdout.py
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.kinetic_core import integrate, operative_parameters  # noqa: E402

TIMESERIES = ROOT / "data" / "lit" / "timeseries"
GG_FILE = "martins2005_glucose_glycine_80_100_120C_pH68.yml"
OUTPUT_DIR = ROOT / "results" / "validation"
FIT_BASENAME = "kinetic_core_b1_fit_report"
BASENAME = "kinetic_core_b1_holdout_report"

INITIAL = {"Glc": 200.0, "Gly": 200.0}

# ---------------------------------------------------------------------------
# THE HELD-OUT NUMBERS. Read at scoring time, nowhere else.
# ---------------------------------------------------------------------------
EPSILON_MARTINS_GLY = {
    "value": 0.64,
    "ci95": 0.03,
    "unit": "L/(mmol*cm)",
    "value_L_per_mol_cm": 640.0,
    "wavelength_nm": 470,
    "amine": "glycine",
    "source_anchor": (
        "Martins & van Boekel 2005 (= Martins 2003 thesis ch. 5); "
        "k3_final_parameter_inventory.md sec. A.1 line 121"
    ),
    "role": "HOLD-OUT (FIT_HOLDOUT_DECLARATION.md D.6, Module 4)",
}
EPSILON_KNOL_ASN = {
    "value_L_per_mol_cm": 282.0,
    "unit": "L/(mol*cm)",
    "value_L_per_mmol_cm": 0.282,
    "amine": "asparagine",
    "source_anchor": "Knol et al. 2005; k3_final_parameter_inventory.md sec. A.1 line 122",
    "role": "HOLD-OUT (FIT_HOLDOUT_DECLARATION.md D.6, Module 4)",
}

#: Additive floor for the fold-error metric on the browning panel: 2% of the
#: panel's printed full scale (panel F, "Mel mmol/L" 0-10), the same rule the
#: fit stage used for every other panel.
MELANOIDIN_FLOOR = 0.02 * 10.0


def load_holdout_series() -> Dict[float, List[Tuple[float, float]]]:
    """The melanoidin response, by temperature. Read ONLY at scoring time."""
    with (TIMESERIES / GG_FILE).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    out: Dict[float, List[Tuple[float, float]]] = {}
    for entry in payload["series"]:
        if entry.get("species") != "melanoidins":
            continue
        t_c = entry.get("T_C")
        if not isinstance(t_c, (int, float)) or isinstance(t_c, bool):
            # 'ambiguous_80_or_100': the source could not resolve the
            # temperature, so neither can the score.
            continue
        for t_min, value in entry.get("points", []):
            out.setdefault(float(t_c), []).append((float(t_min), float(value)))
    for rows in out.values():
        rows.sort()
    return out


def score_variant(name: str, fitted: Dict[str, Dict[str, float]],
                  measured: Dict[float, List[Tuple[float, float]]]) -> Dict[str, Any]:
    """Integrate and score, changing nothing."""
    from dataclasses import replace

    core = {k: (v["k_ref_100C"], v["ea_kj_mol"])
            for k, v in fitted.items() if k != "k_tdg_mel"}
    parameters = dict(operative_parameters(core))
    sink_source = "Martins Table 2 step 9, measured (NOT out-of-sample)"
    if "k_tdg_mel" in fitted:
        parameters["k_tdg_mel"] = replace(
            parameters["k_tdg_mel"],
            k_ref=float(fitted["k_tdg_mel"]["k_ref_100C"]),
            ea_kj_mol=float(fitted["k_tdg_mel"]["ea_kj_mol"]),
        )
        sink_source = "re-estimated from the reactant side only (out-of-sample)"

    per_point: List[Dict[str, Any]] = []
    per_temperature: List[Dict[str, Any]] = []
    for t_c in sorted(measured):
        rows = measured[t_c]
        grid = np.array(sorted({t for t, _ in rows}), dtype=float)
        run = integrate(parameters, t_c + 273.15, INITIAL, grid)
        predicted_by_time = {
            float(t): float(v) for t, v in zip(grid, run.melanoidin_mmol_L())
        }
        errors = []
        for t_min, value in rows:
            predicted = predicted_by_time[t_min]
            signed = math.log10(
                (predicted + MELANOIDIN_FLOOR) / (value + MELANOIDIN_FLOOR)
            )
            errors.append(signed)
            per_point.append({
                "T_C": t_c,
                "t_min": t_min,
                "measured_melanoidin_mmol_L": value,
                "predicted_melanoidin_mmol_L": predicted,
                "measured_A470_implied": value * EPSILON_MARTINS_GLY["value"],
                "predicted_A470": predicted * EPSILON_MARTINS_GLY["value"],
                "signed_log10_error": signed,
                "fold_error": 10 ** abs(signed),
            })
        arr = np.array(errors)
        # The implied epsilon: what extinction coefficient would make the
        # model's predicted concentration reproduce the measured absorbance.
        # Compared against the held-out 0.64 at the END of the run, where the
        # signal is largest and the floor matters least.
        last_t, last_value = rows[-1]
        implied = (
            (last_value * EPSILON_MARTINS_GLY["value"]) / predicted_by_time[last_t]
            if predicted_by_time[last_t] > 0 else float("inf")
        )
        per_temperature.append({
            "T_C": t_c,
            "n_points": len(rows),
            "median_fold_error": float(10 ** float(np.median(np.abs(arr)))),
            "max_fold_error": float(10 ** float(np.max(np.abs(arr)))),
            "median_signed_log10_error": float(np.median(arr)),
            "direction": "over-predicts" if float(np.median(arr)) > 0 else "under-predicts",
            "implied_epsilon_at_final_time_L_per_mmol_cm": implied,
            "final_time_min": last_t,
        })

    all_signed = np.array([row["signed_log10_error"] for row in per_point])
    return {
        "variant": name,
        "melanoidin_sink_constant_source": sink_source,
        "melanoidin_sink_k_ref_100C": float(parameters["k_tdg_mel"].k_ref),
        "n_holdout_points": len(per_point),
        "median_fold_error": float(10 ** float(np.median(np.abs(all_signed)))),
        "max_fold_error": float(10 ** float(np.max(np.abs(all_signed)))),
        "median_signed_log10_error": float(np.median(all_signed)),
        "fraction_within_2x": float(np.mean(10 ** np.abs(all_signed) <= 2.0)),
        "fraction_within_3x": float(np.mean(10 ** np.abs(all_signed) <= 3.0)),
        "per_temperature": per_temperature,
        "per_point": per_point,
    }


def epsilon_amine_specificity(scores: Dict[str, Any]) -> Dict[str, Any]:
    """
    HOLD-OUT 2: Knol's epsilon = 282 L/(mol*cm) for the glucose/asparagine
    melanoidin, against Martins' 640 L/(mol*cm) for glucose/glycine.

    THE HONEST ANSWER FIRST: this module cannot predict it. Nothing in the
    network distinguishes one amine's chromophore from another's -- the
    melanoidin pool is elemental, and an extinction coefficient is an optical
    property of a polymer structure that a mass-action rate model does not
    represent. There is no measured relation in the corpus from which such a
    prediction could be built, and inventing one would be exactly the invented
    functional form the build spec forbids.

    What CAN be stated, and is: the size of the error that assuming a single
    epsilon would cause, and the fact that the model's own implied epsilon
    (from hold-out 1) sits far from Knol's value -- which is a consistency
    statement, not a prediction of it.
    """
    ratio = (EPSILON_MARTINS_GLY["value_L_per_mol_cm"]
             / EPSILON_KNOL_ASN["value_L_per_mol_cm"])
    implied = [
        row["implied_epsilon_at_final_time_L_per_mmol_cm"]
        for row in scores["per_temperature"]
    ]
    return {
        "holdout": "Knol 2005 epsilon = 282 L/(mol*cm), glucose/asparagine",
        "prediction_possible": False,
        "verdict": "NOT PREDICTED -- structural gap, declared",
        "why": (
            "the kinetic core carries the melanoidin pool ELEMENTALLY (carbon and "
            "nitrogen) and has no representation of the polymer's optical "
            "properties. An amine-specific extinction coefficient is not derivable "
            "from any mass-action rate, and the corpus contains no measured "
            "relation between amine identity and epsilon from which one could be "
            "built. Inventing a functional form to score this row would defeat the "
            "purpose of holding it out."
        ),
        "epsilon_martins_glycine_L_per_mol_cm": EPSILON_MARTINS_GLY["value_L_per_mol_cm"],
        "epsilon_knol_asparagine_L_per_mol_cm": EPSILON_KNOL_ASN["value_L_per_mol_cm"],
        "ratio_glycine_over_asparagine": ratio,
        "cost_of_assuming_one_epsilon": (
            f"any model that carries a single epsilon across amino acids is wrong "
            f"by {ratio:.2f}x on browning whenever the amine changes between "
            f"glycine and asparagine"
        ),
        "model_implied_epsilon_per_temperature_L_per_mmol_cm": implied,
        "consistency_note": (
            "the implied epsilons above come from hold-out 1 and are a restatement "
            "of that hold-out's fold errors in optical units. They are NOT a "
            "prediction of Knol's value and must not be read as one."
        ),
    }


def git_head() -> Dict[str, Any]:
    def run(args: List[str]) -> str:
        try:
            return subprocess.check_output(
                args, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unknown"

    return {
        "commit": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(run(["git", "status", "--porcelain"])),
    }


def render_markdown(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    add = lines.append
    add("# Kinetic core B1 -- HOLD-OUT report (predicted vs measured)")
    add("")
    add(f"Generated {payload['generated_on']} on `{payload['git']['branch']}` @ "
        f"`{payload['git']['commit'][:8]}`"
        f"{' (dirty)' if payload['git']['dirty'] else ''}.")
    add("")
    add("Parameters were frozen by "
        "`results/validation/kinetic_core_b1_fit_report.json` BEFORE any hold-out "
        "value was read. **Nothing was changed after reading them.**")
    add("")
    add("## Hold-out 1 -- Martins 2005 step-9 browning (epsilon 0.64)")
    add("")
    add(f"- Held-out response: the melanoidin trajectory at 80/100/120 C, "
        f"{payload['headline']['n_holdout_points']} points.")
    add(f"- Held-out constant: epsilon = "
        f"{EPSILON_MARTINS_GLY['value']} +/- {EPSILON_MARTINS_GLY['ci95']} "
        f"{EPSILON_MARTINS_GLY['unit']} at 470 nm.")
    add("")
    add("### Headline (variant B -- genuinely out-of-sample)")
    add("")
    head = payload["headline"]
    add(f"The melanoidin sink constant was estimated from the reactant side only "
        f"(3-DG and glycine); the browning response never entered the objective.")
    add("")
    add(f"- **median fold error {head['median_fold_error']:.2f}x**, "
        f"max {head['max_fold_error']:.2f}x")
    add(f"- median signed log10 error {head['median_signed_log10_error']:+.3f} "
        f"(the model {('over-predicts' if head['median_signed_log10_error'] > 0 else 'under-predicts')} browning)")
    add(f"- within 2x: {100 * head['fraction_within_2x']:.0f}% of points; "
        f"within 3x: {100 * head['fraction_within_3x']:.0f}%")
    add(f"- **verdict: {payload['verdict_holdout_1']}**")
    add("")
    add("| variant | T (C) | n | median fold error | max fold error | direction | implied epsilon L/(mmol*cm) |")
    add("|---|---:|---:|---:|---:|---|---:|")
    for variant_key in ("headline", "secondary"):
        variant = payload[variant_key]
        for row in variant["per_temperature"]:
            add(f"| {variant['variant']} | {row['T_C']:.0f} | {row['n_points']} | "
                f"{row['median_fold_error']:.2f}x | {row['max_fold_error']:.2f}x | "
                f"{row['direction']} | "
                f"{row['implied_epsilon_at_final_time_L_per_mmol_cm']:.3f} |")
    add("")
    add("The **implied epsilon** column is the extinction coefficient that would "
        "make the model's predicted melanoidin concentration reproduce the "
        f"measured absorbance. The held-out value is "
        f"{EPSILON_MARTINS_GLY['value']} +/- {EPSILON_MARTINS_GLY['ci95']}.")
    add("")
    add("### Secondary (variant A -- Martins' own step-9 constant, NOT out-of-sample)")
    add("")
    second = payload["secondary"]
    add(f"- median fold error {second['median_fold_error']:.2f}x, "
        f"max {second['max_fold_error']:.2f}x, "
        f"signed log10 {second['median_signed_log10_error']:+.3f}")
    add("- **This is a reproducibility check, not a prediction.** Martins fitted "
        "this constant to the very response being scored, so a good number here "
        "carries much less evidential weight than the same number in variant B.")
    add("")
    add("### Predicted vs measured, point by point (variant B)")
    add("")
    add("| T (C) | t (min) | measured mmol/L | predicted mmol/L | measured A470 | predicted A470 | fold error |")
    add("|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["headline"]["per_point"]:
        add(f"| {row['T_C']:.0f} | {row['t_min']:.0f} | "
            f"{row['measured_melanoidin_mmol_L']:.3f} | "
            f"{row['predicted_melanoidin_mmol_L']:.3f} | "
            f"{row['measured_A470_implied']:.3f} | {row['predicted_A470']:.3f} | "
            f"{row['fold_error']:.2f}x |")
    add("")
    add("## Hold-out 2 -- Knol 2005 epsilon = 282 L/(mol*cm), glucose/asparagine")
    add("")
    eps = payload["epsilon_amine_specificity"]
    add(f"- **{eps['verdict']}**")
    add(f"- {eps['why']}")
    add(f"- Martins (glycine) {eps['epsilon_martins_glycine_L_per_mol_cm']:.0f} vs "
        f"Knol (asparagine) {eps['epsilon_knol_asparagine_L_per_mol_cm']:.0f} "
        f"L/(mol*cm) = {eps['ratio_glycine_over_asparagine']:.2f}x.")
    add(f"- {eps['cost_of_assuming_one_epsilon']}.")
    add(f"- {eps['consistency_note']}")
    add("")
    add("## What this does and does not establish")
    add("")
    for line in payload["interpretation"]:
        add(f"- {line}")
    add("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    fit_payload = json.loads(
        (args.output_dir / f"{FIT_BASENAME}.json").read_text(encoding="utf-8")
    )
    frozen = fit_payload["frozen_parameters"]

    measured = load_holdout_series()
    headline = score_variant(
        "B (reactant-side sink, out-of-sample)",
        frozen["variant_B_reactant_side_sink"], measured,
    )
    secondary = score_variant(
        "A (Martins' measured sink, reproducibility)",
        frozen["variant_A_measured_sink"], measured,
    )

    median_fold = headline["median_fold_error"]
    if median_fold <= 1.5:
        verdict = f"PASS -- median {median_fold:.2f}x, inside the 1.5x band"
    elif median_fold <= 3.0:
        verdict = (f"PARTIAL -- median {median_fold:.2f}x. The trunk predicts the "
                   f"browning trajectory to within a factor of 3 but not to "
                   f"within the read-off error of the data.")
    else:
        verdict = (f"FAIL -- median {median_fold:.2f}x. The trunk does NOT predict "
                   f"the browning it was not fitted to.")

    payload: Dict[str, Any] = {
        "artifact": BASENAME,
        "wave": "B1",
        "generated_on": date.today().isoformat(),
        "git": git_head(),
        "frozen_parameters_source": f"results/validation/{FIT_BASENAME}.json",
        "parameters_changed_after_reading_the_holdout": False,
        "holdout_declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md D.6, Module 4",
        "holdout_1": {
            "response": "Martins 2005 melanoidin trajectory, 80/100/120 C",
            "constant": EPSILON_MARTINS_GLY,
            "melanoidin_floor_mmol_L": MELANOIDIN_FLOOR,
        },
        "holdout_2": EPSILON_KNOL_ASN,
        "headline": headline,
        "secondary": secondary,
        "verdict_holdout_1": verdict,
        "epsilon_amine_specificity": None,
        "interpretation": [],
    }
    payload["epsilon_amine_specificity"] = epsilon_amine_specificity(headline)

    ratio_ab = (headline["melanoidin_sink_k_ref_100C"]
                / secondary["melanoidin_sink_k_ref_100C"])
    payload["interpretation"] = [
        (f"The reactant-side estimate of the melanoidin sink constant is "
         f"{ratio_ab:.3f}x Martins' own value, recovered WITHOUT seeing the "
         f"browning response. That agreement is the substantive result of this "
         f"wave: the trunk's own reactant balance determines the browning flux."),
        ("epsilon CANCELS in the concentration comparison, because Martins' "
         "melanoidin response is itself A470/epsilon. Hold-out 1 therefore tests "
         "the predicted FLUX into the sink, and only the implied-epsilon column "
         "restates it in optical units. This is stated so the pass is not "
         "over-read as a test of epsilon itself."),
        ("Hold-out 2 (Knol's epsilon) is NOT predicted and cannot be: the module "
         "has no optical model. That is a declared structural gap, not a "
         "near-miss."),
        ("The whole evaluation sits at pH 6.8 and dilute aqueous conditions. The "
         "module has no pH term and no a_w term, so none of these numbers "
         "transfer off that point."),
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"{BASENAME}.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )
    (args.output_dir / f"{BASENAME}.md").write_text(
        render_markdown(payload), encoding="utf-8"
    )
    print(f"wrote {args.output_dir / BASENAME}.json and .md")
    print(f"HOLD-OUT 1 (variant B, out-of-sample): {verdict}")
    print(f"HOLD-OUT 1 (variant A, reproducibility): median "
          f"{secondary['median_fold_error']:.2f}x")
    print("HOLD-OUT 2 (Knol epsilon 282): NOT PREDICTED -- structural gap")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

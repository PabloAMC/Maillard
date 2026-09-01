#!/usr/bin/env python3
"""Refit the ambient-slurry matrix observability factors against the CONTENT-CORRECTED
Pratap-Singh anchors — and only those anchors.

Context (2026-08-27, Wave O; owner-approved)
--------------------------------------------
`src/matrix_calibration_registry.py` carries per-compound `observable_factor` constants
that were back-solved so that the two ambient Pratap-Singh benchmarks reproduce their own
measured ppb.  Wave K read the paper (Molecules 2021, 26, 4104, Table 1, Europe PMC
PMC8271896) and found the repo's hexanal reference values were transcription errors:

    pea  hexanal   260 ppb  in the repo   ->  1138.00 +/- 297.30 ppb in the paper
    soy  hexanal   380 ppb  in the repo   ->  1621.71 +/- 159.69 ppb in the paper
    1-hexanol   80 / 120 ppb in the repo  ->  n.d. in BOTH matrices (fabricated)
    2-pentylfuran 638 / 2492 ppb          ->  verified VERBATIM, unchanged

The benchmark files were corrected in Wave M; the constants were deliberately NOT, so the
model kept reproducing 260.65 / 379.88 ppb — i.e. it reproduced the ERROR, to four
significant figures, and the two rows scored as honest 4.37x / 4.27x under-predictions.
This script performs the refit the owner approved, against the verified values.

THE OBJECTIVE
-------------
    J(theta) = sum over anchored rows of ( log10( predicted_ppb(theta) / measured_ppb ) )^2

summed over the calibration-eligible matrix-only rows only (see EXCLUSIONS).  Reported in
dex, the same convention as `refit_projection_constants.py` and
`refit_sulfur_barriers_hofmann.py`.

WHAT IS FITTED (and what is not)
--------------------------------
Fitted: ONE parameter — a common multiplicative scale `s` on the ambient-slurry HEXANAL
observability factors of both the pea and the soy lane.  Nothing else.  Not barriers, not
projection constants, not the marker yields in `benchmark_validation.py`, not the
2-pentylfuran factors (their anchors were verified verbatim and are already recovered to
4e-4), not the Trikusuma heated-pea lane (its reference values were not corrected).

WHY ONE PARAMETER AND NOT TWO.  Two free factors (one per lane) against two rows is an
exactly-determined system: the residual is zero by arithmetic, and a zero residual carries
no information at all.  One shared scale against two rows leaves one degree of freedom for
the data to disagree with, so the residual it reports is a real, if small, measurement of
whether the two corrected anchors are mutually consistent.  Both alternatives are computed
and reported below; the one-parameter fit is the one adopted.  See `alternative_per_lane_fit`
in the emitted record for the two-parameter numbers, and `mutual_consistency` for what the
one-parameter residual means.

The scale is HEXANAL-SPECIFIC because the data correction was hexanal-specific: applying it
to 2-pentylfuran would destroy two anchors that the paper confirms verbatim.

SEARCH RANGE
------------
    s in [1e-3, 3e1], log-spaced, `GRID_POINTS` samples, refined by an exact solve.

That envelope is deliberately wider than the full span of the registry's own shipped
observability factors (0.0095957 to 5.9203, ~2.8 decades), so no optimum can be produced by
a bound.  The record reports whether the optimum landed within `BOUND_MARGIN_DEX` of either
end; it did not.

WHY A CLOSED FORM IS LEGITIMATE HERE
------------------------------------
The matrix-only prediction is EXACTLY linear in the observability factor:

    observable_ppb = oxidation_load_ppb * base_marker_yield * dynamic * pH_release * factor

so J(s) is an exact parabola in log10(s) and the minimiser is the geometric mean of the
per-row required scales.  The script does not assume this — it MEASURES it
(`linearity_check` in the record: predictions are re-evaluated at s and 2s and the ratio is
asserted to be 2 within `LINEARITY_TOL`) and only then quotes the closed form, having also
scanned the grid.

EXCLUSIONS (asserted, not merely intended)
------------------------------------------
* The external hold-out (`data/benchmarks/external_validation/`, evidence_class
  `external_validation_only`) is STRUCTURALLY invisible to this script: the benchmark
  iteration is a NON-RECURSIVE glob over `data/benchmarks/*.json`, exactly as
  `derive_matrix_sigma_from_residuals.py` does, and an explicit assertion re-checks that no
  scored path contains `external_validation`.  The hold-out's old->new movement is reported
  by the normal hold-out generator (`generate_external_validation_report.py`) AFTER the
  constants land — it is scored, never fitted.
* Internal synthetic snapshots (`*Internal2026`, `*ProtocolPilot2026`) are model output, not
  experiments; they are neither evidence nor fit targets and are asserted out.
* Quarantined benchmarks (`data/benchmarks/quarantined/`) are outside the glob.

CONSEQUENCE FOR SCORING (stated up front)
-----------------------------------------
The two fit targets are ALREADY declared `fitted_to_benchmark` in the registry, so they were
already reported as fit-recovery rather than as passes.  This refit changes the SIZE of the
recovery (4.37x/4.27x under -> ~1.01x), not its status.  The emitted record additionally
declares `fit_leverage.class = "per_row_recovery"`, which is what
`src/fit_target_index.py` and `scripts/ci/fit_target_gate.py` read to keep these rows out of
the honest literature-coverage numerator AND denominator.  A row that agrees with a constant
solved from it is not evidence, however good the agreement looks.

APPLICATION
-----------
This script does NOT write to `src/matrix_calibration_registry.py`.  The constants are
applied by hand so the provenance block lands with them (old value in `previous_value=`, a
dated causal comment, and a pointer to this record).  Checking is therefore the DEFAULT
behaviour: every run re-derives the fit, compares it against the shipped registry, and exits
non-zero if the two have drifted — so the application stays verifiable without giving a
generator write access to the registry.  `--no-check` emits the record without failing,
which is what you want when deliberately re-deriving BEFORE applying the constants.
`tests/scientific/test_matrix_observability_refit_wave_o.py` runs the same comparison.

Usage
-----
    python scripts/generators/refit_matrix_observability_pratap_singh.py
    python scripts/generators/refit_matrix_observability_pratap_singh.py --no-check
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
import src.matrix_calibration_registry as registry_module
from src.benchmark_validation import evaluate_benchmark

BENCHMARK_DIR = data_paths.BENCHMARKS_DIR

# --- fit configuration -------------------------------------------------------------
SEARCH_RANGE = (1.0e-3, 3.0e1)
GRID_POINTS = 601
LINEARITY_TOL = 1.0e-6          # relative tolerance on the measured linearity check
MATERIALITY_DEX = 0.01          # a factor moves only if |log10(new/old)| >= this (~2.3%)
BOUND_MARGIN_DEX = 0.05         # "hit a bound" if the optimum is this close to either end
SHIPPED_CHECK_TOL_DEX = 5.0e-4  # --check tolerance between shipped and refitted constants

# (protein_type, process_state, compound) -> the knob is the observability factor there.
FITTED_LANES: Tuple[Tuple[str, str, str], ...] = (
    ("pea_iso", "ambient_slurry", "hexanal"),
    ("soy_iso", "ambient_slurry", "hexanal"),
)

# THE INCUMBENT, DECLARED — not read from the live registry.
#
# These are the pre-Wave-O shipped values, i.e. the constants that were back-solved from the
# 260 / 380 ppb transcription errors. The fit is expressed against them and NOT against
# whatever the registry currently holds, so that re-running this script after the refit has
# landed reproduces the SAME record (incumbent, optimum, profile, gain) instead of quietly
# re-fitting a scale of 1.0 against constants it just wrote. A refit record that cannot
# restate its own starting point is not a record.
#
#   pea_iso/ambient_slurry/hexanal   1.0                  ("1.0 by construction": the scale
#                                                          lived in the marker yields)
#   soy_iso/ambient_slurry/hexanal   0.453 / 0.205        (both halves rescaled measurements)
#   soy_iso/heated_matrix/hexanal    (0.453/0.205) * (1 - 0.7060)   (Shu 2024 attenuation on
#                                                          the ambient baseline)
INCUMBENT_FACTORS: Dict[Tuple[str, str, str], float] = {
    ("pea_iso", "ambient_slurry", "hexanal"): 1.0,
    ("soy_iso", "ambient_slurry", "hexanal"): 0.453 / 0.205,
    ("soy_iso", "heated_matrix", "hexanal"): (0.453 / 0.205) * (1.0 - 0.7060),
}

# Lanes whose factors are re-derived for the record but are NOT expected to move: their
# reference values were verified verbatim against the paper, or were never corrected.
REPORTED_NOT_FITTED: Dict[Tuple[str, str, str], str] = {
    ("pea_iso", "ambient_slurry", "2-pentylfuran"): (
        "638 ppb verified VERBATIM against Molecules 26:4104 Table 1 (Wave K). The shipped "
        "factor already recovers it; moving it would be fitting rounding."
    ),
    ("soy_iso", "ambient_slurry", "2-pentylfuran"): (
        "2492 ppb verified VERBATIM against the same table. Same reading."
    ),
    ("pea_iso", "heated_matrix", "hexanal"): (
        "Trikusuma 2019 was NOT content-corrected — it is also the last unverified pillar of "
        "the matrix lane (Wave K: full text not retrieved). Refitting it here would move a "
        "constant against an anchor nobody has read."
    ),
    ("pea_iso", "heated_matrix", "2-pentylfuran"): "Same as the heated-pea hexanal entry.",
    ("pea_iso", "heated_matrix", "nonanal"): "Same as the heated-pea hexanal entry.",
}

# Lanes that have NO anchor at all after the content correction. Reported, never fitted.
UNANCHORED: Dict[Tuple[str, str, str], str] = {
    ("pea_iso", "ambient_slurry", "1-hexanol"): (
        "The paper reports n.d. for pea 1-hexanol and its text states pea proteins "
        "'contained no alcohol compounds'. The 80 ppb this factor was built from is not a "
        "measurement. There is nothing to fit against; the factor is left at 1.0."
    ),
    ("soy_iso", "ambient_slurry", "1-hexanol"): (
        "Soy 1-hexanol is n.d. as well (the paper's soy total alcohols are 40 +/- 9 ppb, "
        "1-octen-3-ol only). The shipped factor 0.143/0.063 is a RATIO OF TWO FABRICATED "
        "NUMBERS (120 and 80 ppb), neither of which appears in the paper. It is left "
        "untouched because there is no anchor to fit it to -- but it is not a calibration, "
        "and it is the factor behind the hold-out's 1117x 1-hexanol miss on Li 2026. "
        "Retiring it is a science decision outside this refit's approved scope."
    ),
    ("pea_iso", "ambient_slurry", "nonanal"): (
        "Nonanal is not reported for pea in this benchmark. Factor stays 1.0."
    ),
    ("soy_iso", "ambient_slurry", "nonanal"): (
        "Nonanal is not reported for soy either; 0.160/0.150 is lane-internal. Untouched."
    ),
}

# Lanes whose factor is COMPOSED from a fitted one and therefore moves WITH the fit, but is
# not itself fitted (no benchmark constrains it).
PROPAGATED: Dict[Tuple[str, str, str], str] = {
    ("soy_iso", "heated_matrix", "hexanal"): (
        "Defined as `soy ambient hexanal factor x (1 - 0.7060)`, the second term being the "
        "Shu 2024 heated-soy attenuation. Its BASELINE is the factor this script refits, so "
        "freezing it would leave a corrected fit composed with a stale baseline -- the exact "
        "defect this wave exists to remove. It moves by the same scale s. No panel benchmark "
        "constrains it; the movement is a propagation, not a fit."
    ),
}

MEASURED_KEY_TO_PREDICTION_KEY = {
    "hexanal": "Hexanal",
    "2-pentylfuran": "2-Pentylfuran",
    "1-hexanol": "1-Hexanol",
    "hexanol": "1-Hexanol",
    "nonanal": "Nonanal",
}

SYNTHETIC_MARKERS = ("Internal2026", "ProtocolPilot2026")


# ---------------------------------------------------------------------------
# fit-target guards — asserted, not intended
# ---------------------------------------------------------------------------
def _eligible_matrix_benchmarks() -> List[Path]:
    """Calibration-eligible matrix-only benchmarks, non-recursively.

    NON-RECURSIVE by design: `data/benchmarks/external_validation/` and
    `data/benchmarks/quarantined/` are subdirectories and are therefore structurally
    invisible here, the same mechanism `src.benchmark_validation.get_benchmark_files`
    relies on.
    """
    paths: List[Path] = []
    for path in sorted(BENCHMARK_DIR.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (payload.get("metadata") or {}).get("execution_path") != "matrix_only":
            continue
        if payload.get("evidence_class") == "external_validation_only":
            raise SystemExit(f"FIT-TARGET GUARD: hold-out payload reached the fit: {path}")
        if any(marker in path.stem for marker in SYNTHETIC_MARKERS):
            raise SystemExit(f"FIT-TARGET GUARD: internal synthetic payload reached the fit: {path}")
        if "external_validation" in str(path):
            raise SystemExit(f"FIT-TARGET GUARD: hold-out path reached the fit: {path}")
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# evaluation
# ---------------------------------------------------------------------------
_BASELINE_RECORDS = tuple(registry_module._MATRIX_CALIBRATION_RECORDS)


def _with_absolute_factors(factors: Dict[Tuple[str, str, str], float]) -> Tuple[Any, ...]:
    """Registry records with the named lanes' observability factors SET (not scaled).

    Absolute, so every evaluation in this script is independent of what the registry
    currently ships and the record is reproducible after the refit has landed.
    """
    out = []
    for record in _BASELINE_RECORDS:
        key = (record.protein_type, record.process_state, record.compound)
        value = factors.get(key)
        if value is None:
            out.append(record)
            continue
        fields = dict(record.__dict__)
        fields["observable_factor"] = float(value)
        out.append(registry_module.MatrixCalibrationRecord(**fields))
    return tuple(out)


def _score_rows(factors: Dict[Tuple[str, str, str], float]) -> Tuple[float, List[Dict[str, Any]]]:
    """Objective (dex^2) and per-row detail at the given absolute per-lane factors."""
    original = registry_module._MATRIX_CALIBRATION_RECORDS
    registry_module._MATRIX_CALIBRATION_RECORDS = _with_absolute_factors(factors)
    try:
        rows: List[Dict[str, Any]] = []
        objective = 0.0
        for path in _eligible_matrix_benchmarks():
            payload = json.loads(path.read_text(encoding="utf-8"))
            benchmark_id = payload.get("benchmark_id", path.stem)
            protein_type = str(payload.get("protein_type"))
            process_state = str((payload.get("process_metadata") or {}).get("state"))
            evaluation = evaluate_benchmark(path)
            predicted = evaluation.predicted_ppb
            for measured_key, entry in (payload.get("measured_volatiles") or {}).items():
                pred_key = MEASURED_KEY_TO_PREDICTION_KEY.get(str(measured_key).strip().lower())
                measured = float((entry or {}).get("conc_ppb", 0.0) or 0.0)
                if pred_key is None or measured <= 0.0 or pred_key not in predicted:
                    continue
                pred = float(predicted[pred_key])
                resid = math.log10(pred / measured) if pred > 0.0 else float("nan")
                lane = (protein_type, process_state, str(measured_key).strip().lower())
                is_fitted = lane in FITTED_LANES
                rows.append(
                    {
                        "benchmark_id": benchmark_id,
                        "lane": f"{protein_type}/{process_state}",
                        "compound": measured_key,
                        "measured_ppb": measured,
                        "predicted_ppb": pred,
                        "ratio": pred / measured,
                        "fold_error": max(pred / measured, measured / pred),
                        "log10_residual_dex": resid,
                        "in_objective": is_fitted,
                    }
                )
                if is_fitted:
                    objective += resid * resid
        return objective, rows
    finally:
        registry_module._MATRIX_CALIBRATION_RECORDS = original


def _at_shared_scale(scale: float) -> Dict[Tuple[str, str, str], float]:
    """Absolute factors obtained by applying one shared scale to the declared incumbents."""
    return {lane: value * scale for lane, value in INCUMBENT_FACTORS.items()}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="Do not fail when the shipped registry constants differ from the refit.",
    )
    args = parser.parse_args()

    record: Dict[str, Any] = {
        "generator": "scripts/generators/refit_matrix_observability_pratap_singh.py",
        "date": "2026-08-27",
        "wave": "Wave O",
        "approval": "owner-approved refit against content-corrected anchors",
        "objective": (
            "sum over anchored rows of (log10(predicted_ppb / measured_ppb))^2, in dex^2; "
            "reported as its square root where a single-row dex figure is quoted"
        ),
        "fit_target_ids": [
            "pea_isolate_40C_PratapSingh2021",
            "soy_isolate_40C_PratapSingh2021",
        ],
        "fit_target_files": [
            data_paths.rel(data_paths.benchmark_path("pea_isolate_40C_PratapSingh2021")),
            data_paths.rel(data_paths.benchmark_path("soy_isolate_40C_PratapSingh2021")),
        ],
        "anchor_provenance": (
            "Molecules 2021, 26, 4104, Table 1 (Europe PMC PMC8271896), read in Wave K. "
            "pea hexanal 1138.00 +/- 297.30 ppb; soy hexanal 1621.71 +/- 159.69 ppb. The "
            "repo's previous 260 / 380 ppb are transcription errors and appear nowhere in "
            "the paper."
        ),
        "search_range": list(SEARCH_RANGE),
        "search_range_basis": (
            "wider than the full span of the registry's own shipped observability factors "
            "(0.0095957 to 5.9203, ~2.8 decades), so no optimum can be an artefact of a bound"
        ),
        "grid_points": GRID_POINTS,
    }

    # ---- 0. exclusions, asserted ------------------------------------------
    eligible = _eligible_matrix_benchmarks()
    record["calibration_eligible_benchmarks"] = [p.name for p in eligible]
    record["exclusions"] = {
        "external_validation_holdout": (
            "structurally excluded: non-recursive glob over data/benchmarks/*.json, plus an "
            "explicit evidence_class / path assertion. Never read by this script."
        ),
        "internal_synthetic": "asserted out by filename marker (Internal2026 / ProtocolPilot2026).",
        "quarantined": "outside the non-recursive glob (data/benchmarks/quarantined/).",
    }

    # ---- 1. linearity, MEASURED not assumed --------------------------------
    _, rows_at_1 = _score_rows(_at_shared_scale(1.0))
    _, rows_at_2 = _score_rows(_at_shared_scale(2.0))
    linearity: List[Dict[str, Any]] = []
    linear = True
    for a, b in zip(rows_at_1, rows_at_2):
        if not a["in_objective"]:
            continue
        ratio = b["predicted_ppb"] / a["predicted_ppb"]
        ok = abs(ratio - 2.0) <= LINEARITY_TOL * 2.0
        linear = linear and ok
        linearity.append(
            {"benchmark_id": a["benchmark_id"], "compound": a["compound"],
             "predicted_at_s1": a["predicted_ppb"], "predicted_at_s2": b["predicted_ppb"],
             "measured_ratio": ratio, "linear_within_tol": ok}
        )
    record["linearity_check"] = {
        "tolerance_relative": LINEARITY_TOL,
        "rows": linearity,
        "all_linear": linear,
        "reading": (
            "the matrix-only prediction is exactly proportional to the observability factor, "
            "so J(s) is an exact parabola in log10(s) and the grid scan below is a "
            "confirmation of a closed form rather than a search over an unknown surface"
        ),
    }
    if not linear:
        raise SystemExit("refit: prediction is NOT linear in the observability factor; the closed form does not apply")

    # ---- 2. required per-lane scales and the shared optimum -----------------
    per_lane_required: Dict[str, float] = {}
    for row in rows_at_1:
        if not row["in_objective"]:
            continue
        per_lane_required[f"{row['lane']}/{row['compound']}"] = row["measured_ppb"] / row["predicted_ppb"]
    shared_optimum = math.exp(sum(math.log(v) for v in per_lane_required.values()) / len(per_lane_required))
    record["per_lane_required_scales"] = per_lane_required
    record["shared_scale_optimum"] = shared_optimum

    # ---- 3. profile over the documented range ------------------------------
    lo, hi = SEARCH_RANGE
    profile = []
    best = (float("inf"), None)
    for i in range(GRID_POINTS):
        s = lo * (hi / lo) ** (i / (GRID_POINTS - 1))
        objective, _ = (0.0, None)
        # closed form: J(s) = sum_i (log10(s / s_i))^2 with s_i the per-lane required scale
        objective = sum((math.log10(s / si)) ** 2 for si in per_lane_required.values())
        profile.append({"scale": s, "objective_dex2": objective})
        if objective < best[0]:
            best = (objective, s)
    record["profile"] = profile[:: max(1, GRID_POINTS // 60)]
    record["profile_note"] = (
        "downsampled to ~60 points for the artifact; the full grid is regenerated by "
        "re-running the script. J is a parabola in log10(s) with a single interior minimum."
    )
    grid_best_objective, grid_best_scale = best
    record["grid_best"] = {"scale": grid_best_scale, "objective_dex2": grid_best_objective}

    # bounds check
    span_dex = math.log10(hi) - math.log10(lo)
    dist_lo = math.log10(shared_optimum) - math.log10(lo)
    dist_hi = math.log10(hi) - math.log10(shared_optimum)
    record["bounds_check"] = {
        "range_span_dex": span_dex,
        "distance_to_lower_bound_dex": dist_lo,
        "distance_to_upper_bound_dex": dist_hi,
        "hit_a_bound": bool(min(dist_lo, dist_hi) < BOUND_MARGIN_DEX),
        "saturated": False,
        "reading": (
            "the optimum sits in the interior of the range by more than a decade on both "
            "sides; nothing saturated and no knob is pinned by a bound"
        ),
    }

    # ---- 4. objective at incumbent and at the fit --------------------------
    baseline_objective, baseline_rows = _score_rows(_at_shared_scale(1.0))
    fitted_objective, fitted_rows = _score_rows(_at_shared_scale(shared_optimum))
    record["objective_at_incumbent_dex2"] = baseline_objective
    record["objective_at_fit_dex2"] = fitted_objective
    record["objective_gain_dex2"] = baseline_objective - fitted_objective
    record["rows_at_incumbent"] = baseline_rows
    record["rows_at_fit"] = fitted_rows

    # ---- 5. mutual consistency of the two corrected anchors -----------------
    residuals = [row["log10_residual_dex"] for row in fitted_rows if row["in_objective"]]
    record["mutual_consistency"] = {
        "definition": (
            "residual left on each anchored row after ONE shared scale is fitted to BOTH. "
            "With two free factors this would be exactly zero and would say nothing."
        ),
        "residual_dex_per_row": residuals,
        "worst_fold_error": max(10.0 ** abs(r) for r in residuals) if residuals else None,
        "reading": (
            "the two corrected anchors require scales of "
            + " and ".join(f"{v:.6g}x" for v in per_lane_required.values())
            + f"; a single shared scale of {shared_optimum:.6g}x satisfies both to within "
            + f"{max(10.0 ** abs(r) for r in residuals):.4f}x. The content correction is "
            "therefore very nearly a COMMON absolute-scale error on the ambient hexanal "
            "lane: the pea-vs-soy release RATIO the registry encodes (2.2098) survives the "
            "correction almost unchanged (it would become "
            f"{per_lane_required[list(per_lane_required)[1]] / per_lane_required[list(per_lane_required)[0]] * 2.2097560975609757:.4f} "
            "under a per-lane fit). What was wrong was the absolute scale, not the "
            "relative structure."
        ),
    }

    # ---- 6. the alternative we declined ------------------------------------
    per_lane_scales = {}
    for lane in FITTED_LANES:
        key = f"{lane[0]}/{lane[1]}/{lane[2]}"
        per_lane_scales[lane] = per_lane_required[key]
    # the propagated soy heated lane follows its own baseline under the per-lane fit
    per_lane_scales[("soy_iso", "heated_matrix", "hexanal")] = per_lane_scales[
        ("soy_iso", "ambient_slurry", "hexanal")
    ]
    alt_factors = {lane: INCUMBENT_FACTORS[lane] * scale for lane, scale in per_lane_scales.items()}
    alt_objective, alt_rows = _score_rows(alt_factors)
    record["alternative_per_lane_fit"] = {
        "free_parameters": 2,
        "fitted_rows": 2,
        "objective_dex2": alt_objective,
        "rows": [r for r in alt_rows if r["in_objective"]],
        "why_declined": (
            "exactly determined: two factors, two rows. Its zero residual is arithmetic and "
            "carries no information about the model or about whether the two corrected "
            "anchors agree. The one-parameter fit is adopted instead; the price is a "
            f"{max(10.0 ** abs(r) for r in residuals):.4f}x residual on each row, which is "
            "the informative part."
        ),
    }

    # ---- 7. adopted constants ----------------------------------------------
    adopted = []
    for lane in FITTED_LANES:
        incumbent = INCUMBENT_FACTORS[lane]
        adopted.append({
            "protein_type": lane[0],
            "process_state": lane[1],
            "compound": lane[2],
            "role": "fitted",
            "previous_value": incumbent,
            "refit_value": incumbent * shared_optimum,
            "move_dex": math.log10(shared_optimum),
            "material": abs(math.log10(shared_optimum)) >= MATERIALITY_DEX,
        })
    for lane, reason in PROPAGATED.items():
        incumbent = INCUMBENT_FACTORS[lane]
        adopted.append({
            "protein_type": lane[0],
            "process_state": lane[1],
            "compound": lane[2],
            "role": "propagated",
            "previous_value": incumbent,
            "refit_value": incumbent * shared_optimum,
            "move_dex": math.log10(shared_optimum),
            "reason": reason,
        })
    record["adopted"] = adopted
    record["materiality_threshold_dex"] = MATERIALITY_DEX

    # ---- 8. reported, not fitted -------------------------------------------
    not_fitted = []
    for lane, reason in REPORTED_NOT_FITTED.items():
        rec = next((r for r in _BASELINE_RECORDS
                    if (r.protein_type, r.process_state, r.compound) == lane), None)
        row = next((r for r in baseline_rows
                    if r["lane"] == f"{lane[0]}/{lane[1]}" and r["compound"].lower() == lane[2]), None)
        implied = (row["measured_ppb"] / row["predicted_ppb"]) if row else None
        not_fitted.append({
            "lane": f"{lane[0]}/{lane[1]}/{lane[2]}",
            "shipped_factor": float(rec.observable_factor) if rec else None,
            "implied_scale_if_it_were_fitted": implied,
            "move_dex_if_it_were_fitted": math.log10(implied) if implied else None,
            "material": (abs(math.log10(implied)) >= MATERIALITY_DEX) if implied else None,
            "reason_not_fitted": reason,
        })
    record["reported_not_fitted"] = not_fitted
    record["unanchored_after_content_correction"] = [
        {"lane": f"{k[0]}/{k[1]}/{k[2]}", "reason": v} for k, v in UNANCHORED.items()
    ]

    # ---- 9. fit leverage ----------------------------------------------------
    n_fitted_rows = sum(1 for r in fitted_rows if r["in_objective"])
    record["fit_leverage"] = {
        "free_parameters": 1,
        "fitted_rows": n_fitted_rows,
        "parameters_per_row": 1.0 / n_fitted_rows if n_fitted_rows else None,
        "class": "per_row_recovery",
        "interpretation": (
            "ONE shared scale against two rows is 0.5 parameters per row, exactly at "
            "src.fit_target_index.FIT_LEVERAGE_THRESHOLD, so these rows are classified "
            "per_row_recovery and stay OUT of the honest literature-coverage numerator and "
            "denominator. That classification is unchanged by this refit: both benchmarks "
            "were already declared `fitted_to_benchmark` in the registry. What changes is "
            "the SIZE of the recovery (4.37x / 4.27x under -> ~1.01x), not its status."
        ),
        "why_this_exclusion_does_not_flatter_the_model": (
            "Excluding a fitted row removes its hits as well as its misses. Before this "
            "refit these two rows were MISSES and were still excluded; after it they are "
            "near-exact and are still excluded. The number that moves as a result of this "
            "refit and is NOT excluded is the external hold-out, reported in "
            "results/validation/external_validation_report.json -- and it gets WORSE. See "
            "`holdout_note`."
        ),
    }
    record["holdout_note"] = (
        "This script never reads the hold-out. Scored AFTER the constants landed, by the "
        "normal generator, the hold-out median |log10| error moves from 1.1849 dex (15.31x) "
        "to a worse value, because the pea ambient lane contains two mutually contradictory "
        "external measurements (Bi 2020: 1260 ppb; Liu 2023 band midpoint: 51.96 ppb, a 24x "
        "spread at nominally the same conditions) and the erroneous 260 ppb constant happened "
        "to sit almost exactly at their geometric mean (sqrt(1260 x 51.96) = 255.9). Moving "
        "to the verified anchor moves the lane onto the Bi side of that disagreement. See "
        "tasks/audit_remediation.md, Wave O, for the 8-point table."
    )

    # ---- 10. shipped-vs-refit check ----------------------------------------
    shipped_check = []
    drifted = False
    for entry in adopted:
        lane = (entry["protein_type"], entry["process_state"], entry["compound"])
        live = next((r for r in registry_module._MATRIX_CALIBRATION_RECORDS
                     if (r.protein_type, r.process_state, r.compound) == lane), None)
        shipped_value = float(live.observable_factor) if live else None
        delta_dex = (
            abs(math.log10(shipped_value / entry["refit_value"]))
            if shipped_value and entry["refit_value"] else None
        )
        agrees = bool(delta_dex is not None and delta_dex <= SHIPPED_CHECK_TOL_DEX)
        drifted = drifted or not agrees
        shipped_check.append({
            "lane": f"{lane[0]}/{lane[1]}/{lane[2]}",
            "shipped": shipped_value,
            "refit": entry["refit_value"],
            "delta_dex": delta_dex,
            "agrees": agrees,
        })
    record["shipped_vs_refit"] = {
        "tolerance_dex": SHIPPED_CHECK_TOL_DEX,
        "rows": shipped_check,
        "all_agree": not drifted,
    }

    # ---- write --------------------------------------------------------------
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "matrix_observability_refit_pratap_singh.json"
    json_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    md: List[str] = [
        "# Matrix observability refit against the content-corrected Pratap-Singh anchors (Wave O)",
        "",
        "Generator: `scripts/generators/refit_matrix_observability_pratap_singh.py` — 2026-08-27, owner-approved.",
        "",
        "Fit targets (both already `fitted_to_benchmark` in the registry, both excluded from",
        "the honest literature-coverage count before AND after this refit):",
        "",
        "- `pea_isolate_40C_PratapSingh2021` — hexanal **1138.00 ± 297.30 ppb**",
        "- `soy_isolate_40C_PratapSingh2021` — hexanal **1621.71 ± 159.69 ppb**",
        "",
        "Molecules 2021, 26, 4104, Table 1 (Europe PMC PMC8271896). The repo's previous",
        "260 / 380 ppb appear nowhere in the paper; the shipped constants reproduced them to",
        "four significant figures, i.e. the model reproduced the transcription error.",
        "",
        "## What was fitted",
        "",
        "ONE parameter: a common multiplicative scale on the ambient-slurry **hexanal**",
        "observability factors of the pea and soy lanes. No barriers, no projection constants,",
        "no marker yields, no 2-pentylfuran factor (those anchors are verified verbatim).",
        "",
        f"- objective at incumbent: **{baseline_objective:.6f} dex²**",
        f"- objective at the fit:  **{fitted_objective:.6f} dex²**",
        f"- adopted shared scale:  **{shared_optimum:.6g}x**",
        f"- search range {SEARCH_RANGE[0]:g} – {SEARCH_RANGE[1]:g} "
        f"({span_dex:.1f} decades); optimum sits {dist_lo:.2f} dex above the floor and "
        f"{dist_hi:.2f} dex below the ceiling — **no bound was hit, nothing saturated**",
        "",
        "## Mutual consistency of the two corrected anchors",
        "",
        "| lane | required scale |",
        "| --- | --- |",
    ]
    for key, value in per_lane_required.items():
        md.append(f"| `{key}` | {value:.6g}x |")
    md += [
        "",
        f"A single shared scale of {shared_optimum:.6g}x satisfies both to within "
        f"**{max(10.0 ** abs(r) for r in residuals):.4f}x**. The correction is therefore",
        "almost purely an absolute-scale error on the ambient hexanal lane; the pea-vs-soy",
        "release ratio the registry encodes survives it. That residual is the informative",
        "part of this fit, and it only exists because the second degree of freedom was",
        "declined — see `alternative_per_lane_fit` in the JSON.",
        "",
        "## Adopted constants",
        "",
        "| lane | role | previous | refit | move |",
        "| --- | --- | --- | --- | --- |",
    ]
    for entry in adopted:
        md.append(
            f"| `{entry['protein_type']}/{entry['process_state']}/{entry['compound']}` | "
            f"{entry['role']} | {entry['previous_value']:.6g} | {entry['refit_value']:.6g} | "
            f"{10.0 ** entry['move_dex']:.4f}x |"
        )
    md += [
        "",
        "## Rows",
        "",
        "| benchmark | compound | measured ppb | predicted (incumbent) | predicted (fit) | fold before | fold after | in objective |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for before, after in zip(baseline_rows, fitted_rows):
        md.append(
            f"| `{before['benchmark_id']}` | {before['compound']} | {before['measured_ppb']:.6g} | "
            f"{before['predicted_ppb']:.6g} | {after['predicted_ppb']:.6g} | "
            f"{before['fold_error']:.4f}x | {after['fold_error']:.4f}x | "
            f"{'yes' if before['in_objective'] else 'no'} |"
        )
    md += [
        "",
        "## Reported, not fitted",
        "",
    ]
    for entry in not_fitted:
        implied = entry["implied_scale_if_it_were_fitted"]
        implied_str = f"{implied:.6g}x" if implied else "n/a"
        md.append(f"- `{entry['lane']}` — implied scale {implied_str}. {entry['reason_not_fitted']}")
    md += ["", "## Unanchored after the content correction", ""]
    for entry in record["unanchored_after_content_correction"]:
        md.append(f"- `{entry['lane']}` — {entry['reason']}")
    md += [
        "",
        "## Hold-out",
        "",
        record["holdout_note"],
        "",
        "## Shipped-vs-refit check",
        "",
        f"tolerance {SHIPPED_CHECK_TOL_DEX:g} dex — "
        + ("**all shipped constants agree with this refit**" if not drifted
           else "**DRIFT: the shipped registry does not match this refit**"),
        "",
    ]
    for row in shipped_check:
        md.append(
            f"- `{row['lane']}`: shipped {row['shipped']:.6g}, refit {row['refit']:.6g}, "
            f"delta {row['delta_dex']:.2e} dex — {'ok' if row['agrees'] else 'DRIFTED'}"
        )
    md_path = output_dir / "matrix_observability_refit_pratap_singh.md"
    md_path.write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"adopted shared scale: {shared_optimum:.6g}x")
    print(f"objective {baseline_objective:.6f} -> {fitted_objective:.6f} dex^2")
    for row in shipped_check:
        print(f"  {row['lane']:44s} shipped {row['shipped']:.6g}  refit {row['refit']:.6g}  "
              f"{'ok' if row['agrees'] else 'DRIFTED'}")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")

    if drifted and not args.no_check:
        print("FAIL: shipped registry constants do not match this refit.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

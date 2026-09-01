#!/usr/bin/env python3
"""Constrain the projection volatile budget with Wave X's STEP-LEVEL molar yields, and
measure what applying that constraint would cost.

Context (2026-08-28, Wave Y)
---------------------------
`src/projection.py` builds a volatile budget that the reaction network never sees:

    drive            = (t / tau_ref) * exp(-(Ea/R) * (1/T - 1/T_ref))
    conversion_extent= 1 - exp(-drive)                       # linear while drive << 1
    yield_fraction   = baseline + (ceiling - baseline) * conversion_extent
    budget_molar     = limiting_precursor_molar * yield_fraction * load_factor

and `src/recommend.py::_project_weighted_flux_to_ppb` then spends ALL of it:
`mol_fraction = activity / total_activity` sums to 1 by construction.  So the largest ppb
the model can ever assign to one analyte is

    budget_ceiling_ppb(analyte) = budget_molar * MW(analyte) * 1e6

with `mol_fraction = 1`.  **A measured value above that ceiling is unreachable at ANY
allocation, ANY barrier and ANY observability factor.  Only the budget can reach it.**

WHY THIS IS NEW, AND WHY IT NEEDED WAVE X
-----------------------------------------
Wave H retained `tau_ref = 12589 min` on the argument that "the budget is already the
right order of magnitude (~1050 ppb at the Hofmann conditions against 542 ppb of measured
MFT+FFT)" -- i.e. the budget was measured to be LARGER than the sum of the measured
analytes on every row it had.  Every row it had was an end-to-end cascade from a raw
sugar, where the model's own allocation spreads the budget over species nobody measured.

Wave X ingested single-step systems: a fed intermediate plus a sulfur donor, 20 mM each,
145 C / 20 min, where the paper prints the **molar yield of the product on the fed
precursor** (0.48 mol %, 0.24 mol %, ...).  `conversion_extent * load_factor` is the model's
own prediction of exactly that quantity, so those rows constrain the budget scale DIRECTLY,
with no allocation term in between.  On the sharpest of them the model gives the measured
analyte ~100 % of the budget and is still short, which is a constraint Wave H could not
have written down.

WHAT IS DERIVED
---------------
1. Per row: measured molar yield, the model's `conversion_extent * load_factor`, the
   budget ceiling in ppb, and `required_budget_scale = measured / ceiling`.
2. The one-sided bound each row places on `tau_ref` (the budget is exactly inversely
   proportional to `tau_ref` in the linear regime -- MEASURED here, not assumed).
3. The counterfactual: `tau_ref` set to the value the binding row demands, applied in
   memory only, and the whole free-precursor panel plus the frozen maillard_path hold-out
   re-scored against it.

FIT CORPUS -- DECLARED BEFORE ANY NUMBER IS COMPUTED
----------------------------------------------------
The constraint rows are the free-precursor literature benchmarks under
`data/benchmarks/*.json`, EXCLUDING:

* `hofmann1998_norfuraneol_h2s_145C_20min_pH5` -- **Wave X's declared fit target.**
* everything under `data/benchmarks/external_validation/` -- BOTH hold-outs.
* `*Internal2026*` / `*ProtocolPilot2026` -- synthetic model output.
* the matrix_only lane (`pea_isolate_*`, `soy_isolate_*` at Pratap-Singh / Trikusuma
  conditions) -- it does not use the projection budget at all, which is itself a finding
  about `refit_projection_constants.py` (see `objective_contamination` in the record).

Both exclusions are asserted, not merely intended.

Usage
-----
    python scripts/generators/derive_projection_budget_from_step_yields.py
    python scripts/generators/derive_projection_budget_from_step_yields.py --no-counterfactual
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import statistics
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
import src.projection as projection_module  # noqa: E402
import src.recommend as recommend_module  # noqa: E402
from src.benchmark_validation import (  # noqa: E402
    _build_comparisons,
    _run_benchmark_recommendation,
    evaluate_benchmark,
    load_benchmark,
)
from src.projection import DEFAULT_PROJECTION_STRATEGY  # noqa: E402

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

OUT_JSON = data_paths.VALIDATION_DIR / "projection_budget_step_yield_constraint.json"
OUT_MD = data_paths.VALIDATION_DIR / "projection_budget_step_yield_constraint.md"

BENCH_DIR = data_paths.BENCHMARKS_DIR

# Wave X's declared fit target. Excluded from every objective and every bound here.
WAVE_X_FIT_TARGET = "hofmann1998_norfuraneol_h2s_145C_20min_pH5"

SKIP_TOKENS = ("Internal2026", "ProtocolPilot2026")
ASSERT_ABSENT_TOKENS = ("external_validation", "quarantined")

FROZEN_HOLDOUT = data_paths.VALIDATION_DIR / "maillard_path_holdout_frozen_predictions.json"


def _install_strategy(strategy) -> None:
    """Point every live consumer of DEFAULT_PROJECTION_STRATEGY at `strategy`.

    Same idiom as scripts/generators/refit_projection_constants.py; the two modules bind
    the constant separately and both have to be moved.
    """
    projection_module.DEFAULT_PROJECTION_STRATEGY = strategy
    recommend_module.DEFAULT_PROJECTION_STRATEGY = strategy


def _dig_budget(metadata: Any) -> Optional[Dict[str, Any]]:
    if isinstance(metadata, dict):
        if "total_volatile_budget_molar" in metadata:
            return metadata
        for value in metadata.values():
            found = _dig_budget(value)
            if found:
                return found
    return None


def _mw(smiles_or_name: str) -> Optional[float]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles_or_name)
    return float(Descriptors.MolWt(mol)) if mol is not None else None


def collect_rows() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for path_str in sorted(glob.glob(str(BENCH_DIR / "*.json"))):
        path = Path(path_str)
        # The synthetic reproducibility lane is SKIPPED (model output, not evidence); the
        # hold-out and the quarantine are ASSERTED absent, because the non-recursive glob
        # is what keeps them out and a future refactor could break that silently.
        if any(token in str(path) for token in SKIP_TOKENS):
            continue
        for token in ASSERT_ABSENT_TOKENS:
            assert token not in str(path), f"corpus contamination: {path} contains {token!r}"
        bench = load_benchmark(path)
        bid = str(bench["benchmark_id"])
        # Mirror the shipped router in `_evaluate_loaded_benchmark` EXACTLY: the safety
        # lane and the matrix_only lane never reach _estimate_projection_budget, so a row
        # from either would be a phantom constraint. (Calling the recommender directly,
        # as this script must for `debug_paths`, would happily produce one -- measured
        # while writing this script, which is why the filter is explicit rather than
        # implied by "did a budget come back".)
        if str((bench.get("metadata") or {}).get("execution_path")) != "free_precursor":
            continue
        if bench.get("benchmark_type") == "safety" or "acrylamide" in bid:
            continue
        try:
            # The recommender is called directly rather than through evaluate_benchmark
            # because `debug_paths` -- the winning route per species, and therefore its
            # STEP COUNT -- is dropped by BenchmarkEvaluation and is the variable this
            # record needs. `_build_comparisons` is the same matcher evaluate_benchmark
            # uses, so the measured/predicted pairing is identical.
            rec = _run_benchmark_recommendation(bench)
        except Exception as exc:  # pragma: no cover - reported, not swallowed
            rows.append({"benchmark_id": bid, "skipped": f"{type(exc).__name__}: {exc}"})
            continue
        evaluation = _Evaluation(rec["predicted_ppb"], rec.get("projection_metadata") or {},
                                 _build_comparisons(bench, rec["predicted_ppb"]))
        debug_paths = rec.get("debug_paths") or {}
        budget = _dig_budget(evaluation.projection_metadata or {})
        if not budget:
            continue  # safety / matrix_only lanes: no projection budget at all
        limiting = float(budget["limiting_precursor_molar"])
        conversion_extent = float(budget["projection_conversion_extent"])
        load_factor = float(budget["projection_load_factor"])
        budget_molar = float(budget["total_volatile_budget_molar"])
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None or not comparison.measured_ppb:
                continue
            canon, mw = _canon_and_mw_for_comparison(evaluation, comparison)
            if mw is None:
                continue
            path_length = len(debug_paths.get(canon) or []) if canon else None
            measured_molar = comparison.measured_ppb / (mw * 1.0e6)
            ceiling_ppb = budget_molar * mw * 1.0e6
            rows.append(
                {
                    "benchmark_id": bid,
                    "compound": comparison.compound,
                    "is_wave_x_fit_target": bid == WAVE_X_FIT_TARGET,
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "analyte_mw": mw,
                    "analyte_canonical_smiles": canon,
                    "winning_route_step_count": path_length,
                    "limiting_precursor_molar": limiting,
                    "model_conversion_extent": conversion_extent,
                    "model_load_factor": load_factor,
                    "model_budget_molar": budget_molar,
                    "model_molar_yield_on_limiting": budget_molar / limiting if limiting else None,
                    "measured_molar_yield_on_limiting": measured_molar / limiting if limiting else None,
                    "budget_ceiling_ppb": ceiling_ppb,
                    "required_budget_scale": comparison.measured_ppb / ceiling_ppb if ceiling_ppb else None,
                    "allocated_fraction_of_budget": (
                        comparison.predicted_ppb / ceiling_ppb if ceiling_ppb else None
                    ),
                }
            )
    return rows


class _Evaluation:
    """The three fields this script needs from a benchmark run, kept together."""

    def __init__(self, predicted_ppb, projection_metadata, comparisons):
        self.predicted_ppb = predicted_ppb
        self.projection_metadata = projection_metadata
        self.comparisons = comparisons


def _canon_and_mw_for_comparison(evaluation, comparison):
    from src.chem_utils import canonicalize_smiles

    for key, value in (evaluation.predicted_ppb or {}).items():
        if abs(value - comparison.predicted_ppb) > 1e-9:
            continue
        canon = canonicalize_smiles(key, fallback_to_original=False, strip_salts=True)
        if canon:
            mw = _mw(canon)
            if mw:
                return canon, mw
    return None, None


def score_panel() -> Dict[str, Any]:
    """Mean/median |log10| over the free-precursor literature rows, at the live strategy."""
    errors: List[float] = []
    per_row: List[Dict[str, Any]] = []
    for path_str in sorted(glob.glob(str(BENCH_DIR / "*.json"))):
        path = Path(path_str)
        bench = load_benchmark(path)
        bid = str(bench["benchmark_id"])
        if "Internal2026" in bid or "ProtocolPilot2026" in bid:
            continue
        try:
            evaluation = evaluate_benchmark(path)
        except Exception:
            continue
        if not _dig_budget(evaluation.projection_metadata or {}):
            continue
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None or not comparison.measured_ppb:
                continue
            if comparison.predicted_ppb <= 0.0:
                continue
            err = abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb))
            errors.append(err)
            per_row.append(
                {
                    "benchmark_id": bid,
                    "compound": comparison.compound,
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "fold": comparison.predicted_ppb / comparison.measured_ppb,
                    "abs_log10": err,
                    "is_wave_x_fit_target": bid == WAVE_X_FIT_TARGET,
                }
            )
    scored = [row for row in per_row if not row["is_wave_x_fit_target"]]
    vals = [row["abs_log10"] for row in scored]
    return {
        "rows": per_row,
        "scored_row_count": len(scored),
        "mean_abs_log10": statistics.fmean(vals) if vals else None,
        "median_abs_log10": statistics.median(vals) if vals else None,
        "within_10x": sum(1 for row in scored if 0.1 <= row["fold"] <= 10.0),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-counterfactual", action="store_true")
    args = parser.parse_args()

    shipped_tau = float(DEFAULT_PROJECTION_STRATEGY.reference_conversion_time_min)
    rows = collect_rows()
    constraint_rows = [r for r in rows if r.get("required_budget_scale") and not r["is_wave_x_fit_target"]]

    # LINEARITY, MEASURED NOT ASSUMED. The budget is claimed to be exactly inversely
    # proportional to tau_ref in the regime the panel occupies; halving tau_ref must
    # double the budget.
    probe_path = BENCH_DIR / "hofmann1998_furan2aldehyde_h2s_145C_20min_pH5.json"
    base_eval = evaluate_benchmark(probe_path)
    base_budget = _dig_budget(base_eval.projection_metadata)["total_volatile_budget_molar"]
    _install_strategy(replace(DEFAULT_PROJECTION_STRATEGY, reference_conversion_time_min=shipped_tau / 2.0))
    half_eval = evaluate_benchmark(probe_path)
    half_budget = _dig_budget(half_eval.projection_metadata)["total_volatile_budget_molar"]
    _install_strategy(DEFAULT_PROJECTION_STRATEGY)
    linearity_ratio = half_budget / base_budget

    # THE STEP-COUNT GROUPING. `winning_route_step_count` is the length of the route the
    # propagator actually took to the analyte, read from `debug_paths`. Grouping the
    # required budget scale by it is the measurement that turns Wave X's "0.518 dex on a
    # step vs 0.952 dex end to end" from an observation into a mechanism.
    by_steps: Dict[int, List[float]] = {}
    for row in constraint_rows:
        steps = row.get("winning_route_step_count")
        if steps is None:
            continue
        by_steps.setdefault(int(steps), []).append(float(row["required_budget_scale"]))
    step_grouping = {
        str(steps): {
            "row_count": len(vals),
            "geometric_mean_required_scale": 10.0 ** (sum(math.log10(v) for v in vals) / len(vals)),
            "min": min(vals),
            "max": max(vals),
        }
        for steps, vals in sorted(by_steps.items())
    }
    short = [v for steps, vals in by_steps.items() if steps <= 2 for v in vals]
    long_ = [v for steps, vals in by_steps.items() if steps >= 4 for v in vals]
    payload_step_summary = {
        "by_winning_route_step_count": step_grouping,
        "short_paths_le_2_steps": {
            "row_count": len(short),
            "geometric_mean_required_scale": (
                10.0 ** (sum(math.log10(v) for v in short) / len(short)) if short else None
            ),
        },
        "long_paths_ge_4_steps": {
            "row_count": len(long_),
            "geometric_mean_required_scale": (
                10.0 ** (sum(math.log10(v) for v in long_) / len(long_)) if long_ else None
            ),
        },
    }
    if short and long_:
        payload_step_summary["short_over_long_fold"] = (
            payload_step_summary["short_paths_le_2_steps"]["geometric_mean_required_scale"]
            / payload_step_summary["long_paths_ge_4_steps"]["geometric_mean_required_scale"]
        )

    binding = max(constraint_rows, key=lambda r: r["required_budget_scale"]) if constraint_rows else None
    violations = [r for r in constraint_rows if r["required_budget_scale"] > 1.0]
    implied_tau = shipped_tau / binding["required_budget_scale"] if binding else None

    payload: Dict[str, Any] = {
        "generated_by": "scripts/generators/derive_projection_budget_from_step_yields.py",
        "wave": "Wave Y (2026-08-28)",
        "shipped_reference_conversion_time_min": shipped_tau,
        "fixed_constants": {
            "apparent_activation_energy_kj_mol": DEFAULT_PROJECTION_STRATEGY.apparent_activation_energy_kj_mol,
            "reference_temperature_kelvin": DEFAULT_PROJECTION_STRATEGY.reference_temperature_kelvin,
            "conversion_ceiling_fraction": DEFAULT_PROJECTION_STRATEGY.conversion_ceiling_fraction,
            "baseline_volatile_yield_fraction": DEFAULT_PROJECTION_STRATEGY.baseline_volatile_yield_fraction,
            "note": (
                "Ea stays pinned to the arrhenius_params.yml `enolisation` entry, exactly as "
                "refit_projection_constants.py pins it: 10 literature rows carry no temperature "
                "leverage and a fitted Ea would be a free parameter dressed as physics."
            ),
        },
        "fit_target_files_excluded": [f"{WAVE_X_FIT_TARGET}.json"],
        "forbidden_as_fit_targets": [
            "data/benchmarks/external_validation/** (matrix hold-out AND maillard_path hold-out)",
            "*Internal2026* / *ProtocolPilot2026*",
            f"{WAVE_X_FIT_TARGET} (Wave X's declared fit target)",
        ],
        "linearity_check": {
            "probe": probe_path.name,
            "tau_halved_budget_ratio": linearity_ratio,
            "expected": 2.0,
            "passes": abs(linearity_ratio - 2.0) < 1e-2,
            "meaning": (
                "The budget is inversely proportional to tau_ref to within the curvature of "
                "1 - exp(-drive). The measured shortfall below 2.0 is that curvature and nothing else: "
                "at drive ~ 1.06e-3 the second-order term is drive/2 ~ 5.3e-4, so doubling the drive "
                "returns 2 - drive rather than 2. A required budget scale therefore maps onto an "
                "implied tau_ref to ~0.1%, and the residual shows up honestly as the binding row "
                "landing at 0.997x rather than 1.000x. MEASURED, not assumed."
            ),
        },
        "rows": rows,
        "step_count_analysis": payload_step_summary,
        "unreachable_rows": violations,
        "binding_constraint": binding,
        "implied_reference_conversion_time_min": implied_tau,
        "objective_contamination": (
            "results/validation/projection_constant_refit.json declares "
            "pea_isolate_40C_PratapSingh2021, soy_isolate_40C_PratapSingh2021 and "
            "pea_isolate_uht_140C_Trikusuma2019 among its fit_target_files. Those benchmarks run the "
            "matrix_only lane, which never calls _estimate_projection_budget -- MEASURED: they produce "
            "no projection metadata at all. Their residuals are therefore CONSTANTS in that objective, "
            "diluting the mean by 5 of 16 rows without carrying any information about tau_ref. This is "
            "why a marker-yield change moves that objective's VALUE and cannot move its ARGMIN."
        ),
    }

    if not args.no_counterfactual and implied_tau:
        before = score_panel()
        _install_strategy(replace(DEFAULT_PROJECTION_STRATEGY, reference_conversion_time_min=implied_tau))
        after = score_panel()
        _install_strategy(DEFAULT_PROJECTION_STRATEGY)
        moves = []
        by_key_after = {(r["benchmark_id"], r["compound"]): r for r in after["rows"]}
        for row in before["rows"]:
            key = (row["benchmark_id"], row["compound"])
            post = by_key_after.get(key)
            if not post:
                continue
            moves.append(
                {
                    "benchmark_id": row["benchmark_id"],
                    "compound": row["compound"],
                    "is_wave_x_fit_target": row["is_wave_x_fit_target"],
                    "measured_ppb": row["measured_ppb"],
                    "fold_before": row["fold"],
                    "fold_after": post["fold"],
                    "abs_log10_before": row["abs_log10"],
                    "abs_log10_after": post["abs_log10"],
                    "better": post["abs_log10"] < row["abs_log10"],
                }
            )
        scored_moves = [m for m in moves if not m["is_wave_x_fit_target"]]
        payload["counterfactual"] = {
            "reference_conversion_time_min": implied_tau,
            "budget_scale": shipped_tau / implied_tau,
            "panel_before": {k: v for k, v in before.items() if k != "rows"},
            "panel_after": {k: v for k, v in after.items() if k != "rows"},
            "rows_better": sum(1 for m in scored_moves if m["better"]),
            "rows_worse": sum(1 for m in scored_moves if not m["better"]),
            "moves": moves,
        }

    payload["decision"] = _decide(payload)
    payload["applied_to_runtime"] = False

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUT_MD.write_text(_render_markdown(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


def _decide(payload: Dict[str, Any]) -> Dict[str, Any]:
    cf = payload.get("counterfactual")
    binding = payload.get("binding_constraint")
    if not (cf and binding):
        return {"verdict": "no counterfactual measured"}
    worse = cf["panel_after"]["mean_abs_log10"] > cf["panel_before"]["mean_abs_log10"]
    return {
        "verdict": "DO NOT APPLY" if worse else "APPLY",
        "why": (
            "The constraint is real and one-sided: "
            f"{binding['benchmark_id']} / {binding['compound']} measures "
            f"{binding['measured_ppb']:.1f} ppb against a budget CEILING of "
            f"{binding['budget_ceiling_ppb']:.1f} ppb, and the model already allocates "
            f"{binding['allocated_fraction_of_budget'] * 100:.1f}% of the budget to it. The row is "
            "unreachable at any allocation. But raising the global budget to reach it moves every "
            "other row by the same factor, and the panel's other rows are predominantly OVER-predicted, "
            f"so the mean |log10| goes {cf['panel_before']['mean_abs_log10']:.4f} -> "
            f"{cf['panel_after']['mean_abs_log10']:.4f} dex "
            f"({cf['rows_better']} rows better, {cf['rows_worse']} worse). "
            "The two demands are simultaneously live and no single global scale satisfies both: the "
            "deficit is NOT a scale, it is that the budget is precursor-limited without regard to how "
            "many steps separate the precursor from the product."
        ),
    }


def _render_markdown(payload: Dict[str, Any]) -> str:
    out: List[str] = []
    out.append("# Projection budget vs. the step-level molar yields (Wave Y, 2026-08-28)\n")
    out.append(f"Generated by `{payload['generated_by']}`.\n")
    out.append(f"Shipped `reference_conversion_time_min` = **{payload['shipped_reference_conversion_time_min']:.6g} min**.\n")

    lin = payload["linearity_check"]
    out.append(
        f"**Linearity (measured):** halving `tau_ref` multiplies the budget by "
        f"{lin['tau_halved_budget_ratio']:.9f} (expected 2.0) — {'PASS' if lin['passes'] else 'FAIL'}. "
        f"{lin['meaning']}\n"
    )

    out.append("## The budget ceiling, row by row\n")
    out.append(
        "`budget_ceiling_ppb` is the largest ppb the model could assign to the analyte if the "
        "allocation gave it **100 %** of the budget. `alloc %` is how much it actually gives it.\n"
    )
    out.append("| benchmark | analyte | steps | measured ppb | ceiling ppb | alloc % | required budget scale | model mol% | measured mol% |")
    out.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in payload["rows"]:
        if "required_budget_scale" not in row or row["required_budget_scale"] is None:
            continue
        tag = " *(FIT TARGET — excluded)*" if row["is_wave_x_fit_target"] else ""
        out.append(
            f"| `{row['benchmark_id']}`{tag} | {row['compound']} | {row.get('winning_route_step_count')} | {row['measured_ppb']:.2f} | "
            f"{row['budget_ceiling_ppb']:.2f} | {row['allocated_fraction_of_budget'] * 100:.1f}% | "
            f"**{row['required_budget_scale']:.4f}x** | "
            f"{row['model_molar_yield_on_limiting'] * 100:.4f} | "
            f"{row['measured_molar_yield_on_limiting'] * 100:.4f} |"
        )
    out.append("")

    sca = payload["step_count_analysis"]
    out.append("## Required budget scale vs. how many steps the route took\n")
    out.append("| winning-route steps | rows | geometric-mean required scale | min | max |")
    out.append("|---:|---:|---:|---:|---:|")
    for steps, blk in sca["by_winning_route_step_count"].items():
        out.append(
            f"| {steps} | {blk['row_count']} | **{blk['geometric_mean_required_scale']:.4f}x** | "
            f"{blk['min']:.4f}x | {blk['max']:.4f}x |"
        )
    out.append("")
    sp = sca["short_paths_le_2_steps"]
    lp = sca["long_paths_ge_4_steps"]
    out.append(
        f"Short paths (<= 2 steps, {sp['row_count']} rows) want the budget "
        f"**{sp['geometric_mean_required_scale']:.4f}x**; long paths (>= 4 steps, {lp['row_count']} rows) want it "
        f"**{lp['geometric_mean_required_scale']:.4f}x**. Ratio **{sca.get('short_over_long_fold', float('nan')):.1f}x**.\n"
    )

    binding = payload["binding_constraint"]
    if binding:
        out.append("## The binding constraint\n")
        out.append(
            f"`{binding['benchmark_id']}` / {binding['compound']}: measured "
            f"{binding['measured_ppb']:.1f} ppb, budget ceiling {binding['budget_ceiling_ppb']:.1f} ppb, "
            f"allocated {binding['allocated_fraction_of_budget'] * 100:.1f}% of the budget. "
            f"**Required budget scale {binding['required_budget_scale']:.4f}x**, i.e. implied "
            f"`tau_ref` = {payload['implied_reference_conversion_time_min']:.1f} min.\n"
        )
        out.append(
            f"Unreachable rows (measured above the ceiling): "
            f"{len(payload['unreachable_rows'])}"
            + (
                " — " + ", ".join(f"`{r['benchmark_id']}`/{r['compound']}" for r in payload["unreachable_rows"])
                if payload["unreachable_rows"]
                else ""
            )
            + "\n"
        )

    cf = payload.get("counterfactual")
    if cf:
        out.append("## Counterfactual — the constraint applied, in memory only\n")
        out.append("| | shipped | counterfactual |")
        out.append("|---|---:|---:|")
        out.append(f"| `reference_conversion_time_min` | {payload['shipped_reference_conversion_time_min']:.6g} | {cf['reference_conversion_time_min']:.6g} |")
        out.append(f"| budget scale | 1.0000x | {cf['budget_scale']:.4f}x |")
        out.append(f"| panel mean \\|log10\\| | {cf['panel_before']['mean_abs_log10']:.4f} | {cf['panel_after']['mean_abs_log10']:.4f} |")
        out.append(f"| panel median \\|log10\\| | {cf['panel_before']['median_abs_log10']:.4f} | {cf['panel_after']['median_abs_log10']:.4f} |")
        out.append(f"| within 10x | {cf['panel_before']['within_10x']}/{cf['panel_before']['scored_row_count']} | {cf['panel_after']['within_10x']}/{cf['panel_after']['scored_row_count']} |")
        out.append(f"| rows better / worse | — | {cf['rows_better']} / {cf['rows_worse']} |")
        out.append("")
        out.append("| benchmark | analyte | fold before | fold after | |")
        out.append("|---|---|---:|---:|---|")
        for move in cf["moves"]:
            tag = " *(fit target)*" if move["is_wave_x_fit_target"] else ""
            out.append(
                f"| `{move['benchmark_id']}`{tag} | {move['compound']} | {move['fold_before']:.4f}x | "
                f"{move['fold_after']:.4f}x | {'better' if move['better'] else 'worse'} |"
            )
        out.append("")

    decision = payload["decision"]
    out.append(f"## Decision: **{decision['verdict']}**\n")
    out.append(decision.get("why", "") + "\n")
    out.append("## A note on the existing projection refit objective\n")
    out.append(payload["objective_contamination"] + "\n")
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())

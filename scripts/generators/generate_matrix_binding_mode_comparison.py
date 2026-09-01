#!/usr/bin/env python
"""Score the matrix observability lane under BOTH observability modes and compare.

2026-08-27 (Wave S4).  See tasks/audit_remediation.md, "Wave S4".

WHAT THIS ANSWERS.  The matrix lane's observability factors are FITTED -- several
back-solved from the benchmark they are then scored against, two of them
(`1-hexanol`, `0.143 / 0.063`) back-solved from values that appear in no paper.  Wave S4
built the alternative: `src/protein_binding.py`, which computes the same factor from
MEASURED binding constants with zero fitted parameters.  This generator scores both.

THE MODES
  fitted_factors                 -- the incumbent, unchanged.  Default everywhere.
  unit_observability             -- the NULL MODEL: factor 1.0 on every lane.  It is here
                                    because without it a binding-mode win is
                                    unattributable: the binding mode has usable binding
                                    data on only a minority of hold-out rows and reduces
                                    to exactly this on the rest, so part of any
                                    improvement is "the fitted factor was worse than
                                    nothing", not "the binding physics is right".
                                    Scoring all three separates those two claims.
  binding_physics                -- f_free = 1/(1 + a_p*Pow*c_p), applied ONLY where the
                                    lane is an aqueous dispersion with a declared protein
                                    loading AND the reference is water-calibrated.  A
                                    matrix-matched reference measures the TOTAL, so its
                                    correct factor is 1.0 and the mode returns 1.0.

A FOURTH MODE EXISTS AND IS NOT SCORED HERE.  `binding_physics_out_of_domain` would apply
the aqueous Langmuir to the low-water-activity lanes too, so that the size of what the
in-domain rule declines to do would be visible.  It is currently a NO-OP and therefore
not reported: every non-aqueous lane in this repository (Bi roasted pea, Li 2026 HME) also
has NO declared protein loading and NO established quantification basis, so the mode has
nothing to compute with.  What blocks it is missing DATA, not the domain rule.

NOTHING HERE IS FITTED.  This script only evaluates; it writes no constant back anywhere.
The hold-out bundles are read, never modified, and no mode was selected by looking at the
score it produces.

THE OUTPUT FILENAME IS LOAD-BEARING.  `scripts/ci/holdout_guard.py` check 4 treats every
`results/validation/*refit*.json` and `*rederivation*.json` as a FIT RECORD and fails if
one names a hold-out benchmark.  This artifact names every hold-out benchmark and is not a
fit record, so it must never be renamed to contain either substring.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.benchmark_validation import _run_matrix_only_benchmark_prediction  # noqa: E402
from src.external_validation import build_external_validation_report, get_holdout_benchmark_files  # noqa: E402
from src.protein_binding import (  # noqa: E402
    MODE_BINDING,
    MODE_FITTED,
    MODE_UNIT,
    cross_check_percent_bound,
    describe_model,
    free_fraction_from_ap,
    load_binding_constants,
    matrix_protein_loading,
    use_observability_mode,
)

OUT_JSON = data_paths.VALIDATION_DIR / "matrix_binding_mode_comparison.json"
OUT_MD = data_paths.VALIDATION_DIR / "matrix_binding_mode_comparison.md"

MODES = (MODE_FITTED, MODE_UNIT, MODE_BINDING)

# The two rows whose observability factor is composed entirely of fabricated values
# (Wave T3): `0.143 / 0.063` = 120 ppb / 80 ppb, neither of which appears in
# Pratap-Singh et al. (Molecules 2021, 26, 4104).
FABRICATED_FACTOR_ROWS = (
    ("pea_isolate_40C_PratapSingh2021", "1-Hexanol"),
    ("soy_isolate_40C_PratapSingh2021", "1-Hexanol"),
    ("external_validation_li_2026_spi_wg_hme_control", "1-Hexanol"),
)

IN_PANEL_FILES = (
    data_paths.benchmark_path("pea_isolate_40C_PratapSingh2021"),
    data_paths.benchmark_path("soy_isolate_40C_PratapSingh2021"),
    data_paths.benchmark_path("pea_isolate_uht_140C_Trikusuma2019"),
)


def _git_head() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:  # pragma: no cover - informational only
        return "unknown"


def _fold(predicted: Optional[float], measured: Optional[float]) -> Optional[float]:
    if not predicted or not measured or predicted <= 0 or measured <= 0:
        return None
    ratio = predicted / measured
    return max(ratio, 1.0 / ratio)


def _measured_map(bench: Dict[str, Any]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for name, row in (bench.get("measured_volatiles") or {}).items():
        value = row.get("conc_ppb") if isinstance(row, dict) else row
        if value is None:
            continue
        try:
            out[str(name).strip().lower()] = float(value)
        except (TypeError, ValueError):
            continue
    return out


def deterministic_rows(files: List[Path], mode: str) -> List[Dict[str, Any]]:
    """Deterministic (non-Monte-Carlo) predictions, so the modes differ only by mode."""
    rows: List[Dict[str, Any]] = []
    with use_observability_mode(mode):
        for path in files:
            bench = json.loads(Path(path).read_text(encoding="utf-8"))
            result = _run_matrix_only_benchmark_prediction(bench)
            measured = _measured_map(bench)
            for compound, predicted in result["predicted_ppb"].items():
                meta = result["projection_metadata"][compound]
                key = compound.strip().lower()
                rows.append(
                    {
                        "benchmark_id": bench.get("benchmark_id"),
                        "compound": compound,
                        "measured_ppb": measured.get(key),
                        "predicted_ppb": float(predicted),
                        "net_observability_factor": meta.get("total_observable_factor"),
                        "fitted_registry_factor": meta.get("calibration_observable_factor"),
                        "binding_f_free": meta.get("binding_f_free"),
                        "binding_in_domain": meta.get("binding_in_domain"),
                        "binding_mechanism": meta.get("binding_mechanism"),
                        "binding_reasons": meta.get("binding_reasons"),
                        "fold_error": _fold(float(predicted), measured.get(key)),
                    }
                )
    return rows


def _summarise(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    folds = [r["fold_error"] for r in rows if r.get("fold_error")]
    logs = [abs(math.log10(f)) for f in folds]
    return {
        "scored_rows": len(folds),
        "median_fold_error": statistics.median(folds) if folds else None,
        "median_abs_log10_error": statistics.median(logs) if logs else None,
        "within_10x": sum(1 for f in folds if f <= 10.0),
        "worst_fold_error": max(folds) if folds else None,
    }


def holdout_report(mode: str, n_samples: int, seed: int) -> Dict[str, Any]:
    with use_observability_mode(mode):
        payload = build_external_validation_report(n_samples=n_samples, seed=seed)
    rows: List[Dict[str, Any]] = []
    for bench in payload.get("benchmarks", []):
        for compound in bench.get("compounds", []):
            rows.append(
                {
                    "benchmark_id": bench.get("benchmark_id"),
                    "compound": compound.get("compound"),
                    "measured_ppb": compound.get("measured_ppb"),
                    "predicted_p50": compound.get("predicted_p50"),
                    "fold_error": compound.get("fold_error"),
                    "abs_log10_error": compound.get("abs_log10_error"),
                    "inside_ci": compound.get("inside_ci"),
                    "value_provenance": compound.get("value_provenance"),
                }
            )
    summary = payload.get("summary", {})
    return {
        "rows": rows,
        "median_fold_error": summary.get("median_accuracy_fold"),
        "median_abs_log10_error": summary.get("median_abs_log10_error"),
        "ci_coverage_hits": sum(1 for r in rows if r.get("inside_ci")),
        "ci_coverage_total": len(rows),
    }


def pratap_vs_liu_test() -> Dict[str, Any]:
    """Does the binding model explain the Pratap-Singh vs Liu hexanal gap?

    Pratap-Singh spiked standards INTO the slurry (matrix-matched) -> reports the TOTAL.
    Liu built five-point external curves in DI WATER -> reports the FREE-equivalent, i.e.
    total x f_free.  So the binding model makes a falsifiable, zero-parameter prediction:
    a water-calibrated measurement of the SAME material must read 1/f_free times LOWER
    than a matrix-matched one.  This computes that factor and checks it against what the
    two sources actually report.
    """
    liu = matrix_protein_loading("external_validation_liu_2023_ppi_offnote_baseline") or {}
    pratap = matrix_protein_loading("pea_isolate_40C_PratapSingh2021") or {}
    c_p_liu = float(liu.get("protein_concentration_g_per_L") or 0.0)
    f = free_fraction_from_ap(
        "hexanal",
        protein_source="pea_iso",
        protein_concentration_g_per_l=c_p_liu,
        protein_basis=str(liu.get("protein_basis") or "g_protein"),
    )
    under_read = 1.0 / f.f_free if f.f_free else None

    # The measured numbers, both verified in earlier waves and quoted in their bundles.
    pratap_total = 1138.00           # ug/L of the 1:7 w/v slurry, Molecules 2021 Table 1
    liu_band = (2445.0, 52454.0)     # ug/L of the 10%-solids slurry, thesis Table 2.7
    liu_geo_mid = math.sqrt(liu_band[0] * liu_band[1])
    liu_corrected = tuple(v * under_read for v in liu_band) if under_read else None

    return {
        "question": (
            "Does protein binding explain why Liu's water-calibrated pea hexanal numbers "
            "differ from Pratap-Singh's matrix-matched one?"
        ),
        "binding_prediction": {
            "liu_protein_loading_g_per_l": c_p_liu,
            "liu_protein_basis": liu.get("protein_basis"),
            "pratap_protein_loading_g_per_l": pratap.get("protein_concentration_g_per_L"),
            "f_free_hexanal_pea_at_liu_loading": f.f_free,
            "k_eff_l_per_g": f.k_eff_l_per_g,
            "predicted_water_calibrated_under_read_fold": under_read,
            "direction": (
                "Liu's reported values must be LOWER than a matrix-matched measurement of "
                "the same material by this factor."
            ),
        },
        "observed": {
            "pratap_singh_matrix_matched_ppb": pratap_total,
            "liu_water_calibrated_band_ppb": list(liu_band),
            "liu_geometric_mid_ppb": liu_geo_mid,
            "liu_lowest_lot_over_pratap_fold": liu_band[0] / pratap_total,
            "liu_geo_mid_over_pratap_fold": liu_geo_mid / pratap_total,
        },
        "verdict_inputs": {
            "liu_band_corrected_to_total_ppb": list(liu_corrected) if liu_corrected else None,
            "liu_geo_mid_corrected_ppb": (liu_geo_mid * under_read) if under_read else None,
            "gap_before_correction_fold": liu_geo_mid / pratap_total,
            "gap_after_correction_fold": (liu_geo_mid * under_read / pratap_total) if under_read else None,
        },
    }


def build() -> Dict[str, Any]:
    holdout_files = [Path(p) for p in get_holdout_benchmark_files()]
    in_panel_files = list(IN_PANEL_FILES)

    payload: Dict[str, Any] = {
        "generated_by": "scripts/generators/generate_matrix_binding_mode_comparison.py",
        "wave": "S4",
        "git_head": _git_head(),
        "fitted_parameters_in_this_wave": 0,
        "model": describe_model(),
        "binding_constants_file": data_paths.rel(data_paths.BINDING_CONSTANTS),
        "binding_constants_sources": [
            s.get("source_id") for s in (load_binding_constants().get("sources") or [])
        ],
        "modes": list(MODES),
        "deterministic": {},
        "holdout_monte_carlo": {},
        "fabricated_factor_rows": {},
        "cross_check_percent_bound": cross_check_percent_bound(),
        "pratap_vs_liu_binding_test": pratap_vs_liu_test(),
    }

    for mode in MODES:
        holdout_rows = deterministic_rows(holdout_files, mode)
        panel_rows = deterministic_rows(in_panel_files, mode)
        payload["deterministic"][mode] = {
            "holdout": {"rows": holdout_rows, "summary": _summarise(holdout_rows)},
            "in_panel": {"rows": panel_rows, "summary": _summarise(panel_rows)},
        }

    # ATTRIBUTION. A binding-mode win over the fitted mode is not automatically evidence
    # for the binding physics: on any lane where the binding data does not apply, the
    # binding mode IS the null model. This counts which is which so the headline cannot be
    # read as stronger than it is.
    binding_rows = payload["deterministic"][MODE_BINDING]["holdout"]["rows"]
    applied = [r for r in binding_rows if r.get("binding_mechanism") == "harrison_hills_ap_pow"]
    payload["attribution"] = {
        "holdout_rows_total": len(binding_rows),
        "holdout_rows_where_binding_actually_applied": len(applied),
        "holdout_rows_where_binding_reduces_to_null_model": len(binding_rows) - len(applied),
        "rows_where_binding_applied": [
            {"benchmark_id": r["benchmark_id"], "compound": r["compound"], "f_free": r["binding_f_free"]}
            for r in applied
        ],
        "reading": (
            "Any difference between `unit_observability` and `binding_physics` is the "
            "binding physics. Any difference between `fitted_factors` and "
            "`unit_observability` is the fitted factors being better or worse than no "
            "factor at all. Report the two separately."
        ),
    }

    for mode in MODES:
        rows = {
            (r["benchmark_id"], r["compound"]): r
            for r in payload["deterministic"][mode]["holdout"]["rows"]
            + payload["deterministic"][mode]["in_panel"]["rows"]
        }
        payload["fabricated_factor_rows"][mode] = [
            {
                "benchmark_id": bid,
                "compound": compound,
                **{
                    k: rows[(bid, compound)].get(k)
                    for k in (
                        "measured_ppb",
                        "predicted_ppb",
                        "net_observability_factor",
                        "fitted_registry_factor",
                        "binding_f_free",
                        "binding_in_domain",
                        "fold_error",
                    )
                },
            }
            for bid, compound in FABRICATED_FACTOR_ROWS
            if (bid, compound) in rows
        ]

    return payload


def render_markdown(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    lines.append("# Matrix observability: fitted factors vs measured binding physics")
    lines.append("")
    lines.append(
        "<!-- Auto-generated by scripts/generators/generate_matrix_binding_mode_comparison.py. "
        "Manual edits will be overwritten. -->"
    )
    lines.append("")
    lines.append(f"Wave S4, git HEAD `{payload.get('git_head')}`. **Fitted parameters in the binding mode: 0.**")
    lines.append("")
    lines.append(
        "The binding mode computes the matrix observability factor as "
        "`f_free = 1 / (1 + a_p * Pow * c_p)` from constants that were measured in other "
        "laboratories and transcribed with verbatim quotes into "
        "`data/lit/binding_constants.yml`. The fitted mode is the shipped incumbent."
    )
    lines.append("")

    lines.append("## Zero-parameter cross-check of the binding model itself")
    lines.append("")
    lines.append(
        "Before comparing modes: does the binding model reproduce percent-bound "
        "measurements it was not built from? `a_p` was fitted by Snel et al. (2023) on "
        "2-alkanones in demineralised water by APCI-MS. The rows below come from other "
        "laboratories and other methods."
    )
    lines.append("")
    lines.append("| record | compound | c_p g/L | measured % bound | predicted (pea a_p) | residual (pts) |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
    for row in payload["cross_check_percent_bound"]:
        pred = row.get("predicted_percent_bound_pea_iso")
        res = row.get("residual_points_pea_iso")
        lines.append(
            f"| `{row['record_id']}` | {row['compound']} | {row['protein_concentration_g_per_l_used']:.2f} | "
            f"{row['measured_percent_bound']:.2f} | "
            f"{('%.2f' % pred) if pred is not None else 'n/a'} | "
            f"{('%+.2f' % res) if res is not None else 'n/a'} |"
        )
    lines.append("")

    lines.append("## Hold-out, mode vs mode (Monte Carlo, the shipped scoring path)")
    lines.append("")
    mc = payload.get("holdout_monte_carlo", {})
    if mc:
        header = "| point | measured |" + "".join(f" p50 {m} | fold {m} |" for m in MODES)
        lines.append(header)
        lines.append("| --- | ---: |" + "".join(" ---: | ---: |" for _ in MODES))
        keys = [(r["benchmark_id"], r["compound"]) for r in mc[MODES[0]]["rows"]]
        indexed = {m: {(r["benchmark_id"], r["compound"]): r for r in mc[m]["rows"]} for m in MODES}
        for key in keys:
            base = indexed[MODES[0]][key]
            cells = ""
            for m in MODES:
                row = indexed[m].get(key, {})
                p50 = row.get("predicted_p50")
                fold = row.get("fold_error")
                cells += f" {('%.4g' % p50) if p50 else 'n/a'} | {('%.4g' % fold) if fold else 'n/a'} |"
            lines.append(f"| {key[0]} / {key[1]} | {base.get('measured_ppb')} |{cells}")
        lines.append("")
        for m in MODES:
            lines.append(
                f"* **{m}** — median fold error **{mc[m]['median_fold_error']}**, "
                f"median |log10| {mc[m]['median_abs_log10_error']}, "
                f"CI coverage {mc[m]['ci_coverage_hits']}/{mc[m]['ci_coverage_total']}"
            )
        lines.append("")

    lines.append("## Attribution — what the mode difference is actually measuring")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps(payload.get("attribution", {}), indent=2))
    lines.append("```")
    lines.append("")

    lines.append("## Deterministic scores (no Monte Carlo sampling, so the modes differ only by mode)")
    lines.append("")
    lines.append("| surface | mode | scored rows | median fold | within 10x | worst fold |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")
    for surface in ("holdout", "in_panel"):
        for mode in MODES:
            s = payload["deterministic"][mode][surface]["summary"]
            lines.append(
                f"| {surface} | {mode} | {s['scored_rows']} | "
                f"{('%.4g' % s['median_fold_error']) if s['median_fold_error'] else 'n/a'} | "
                f"{s['within_10x']} | "
                f"{('%.4g' % s['worst_fold_error']) if s['worst_fold_error'] else 'n/a'} |"
            )
    lines.append("")

    lines.append("## The rows whose fitted factor is composed of fabricated values")
    lines.append("")
    lines.append(
        "`0.143 / 0.063` is a ratio of 120 ppb and 80 ppb, neither of which appears in "
        "Pratap-Singh et al. (Molecules 2021, 26, 4104); Wave T3 traced both to this "
        "repository's own abstract-reconstructed brief."
    )
    lines.append("")
    lines.append("| benchmark | compound | mode | observability factor | predicted ppb | measured ppb | fold |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: |")
    for mode in MODES:
        for row in payload["fabricated_factor_rows"][mode]:
            lines.append(
                f"| {row['benchmark_id']} | {row['compound']} | {mode} | "
                f"{('%.5g' % row['net_observability_factor']) if row.get('net_observability_factor') else 'n/a'} | "
                f"{('%.5g' % row['predicted_ppb'])} | "
                f"{row.get('measured_ppb') if row.get('measured_ppb') else 'n/a'} | "
                f"{('%.5g' % row['fold_error']) if row.get('fold_error') else 'n/a'} |"
            )
    lines.append("")

    lines.append("## The Pratap-Singh vs Liu test")
    lines.append("")
    test = payload["pratap_vs_liu_binding_test"]
    lines.append("```json")
    lines.append(json.dumps(test, indent=2))
    lines.append("```")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--skip-monte-carlo", action="store_true")
    args = parser.parse_args()

    payload = build()
    if not args.skip_monte_carlo:
        for mode in MODES:
            payload["holdout_monte_carlo"][mode] = holdout_report(mode, args.n_samples, args.seed)
        payload["monte_carlo_settings"] = {"n_samples": args.n_samples, "seed": args.seed}

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(render_markdown(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

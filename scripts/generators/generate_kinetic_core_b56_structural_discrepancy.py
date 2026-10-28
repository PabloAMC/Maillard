#!/usr/bin/env python
"""
WAVE B56 -- the per-lane STRUCTURAL-DISCREPANCY term, estimated from IN-SAMPLE fit residuals only.

Pre-registered in ``results/validation/kinetic_core_b55_b56_prereg.md`` sec. 3, before any coverage
number under it existed. Writes ``results/validation/kinetic_core_b56_structural_discrepancy.{json,md}``,
which ``src/kinetic_core/uncertainty.py`` reads when the B56 switch is on.

THE ESTIMATOR (a declared model-error ASSUMPTION, not a measurement). For each lane L, over the fits the
shipped engine reads for that lane, on the rows those fits read:

    sigma_L^2 = sum_f c_f * sum_{i in f} r_i^2 / sum_f n_used_f,     c_f = n_f / (n_f - p_f)

r_i = log10(predicted / observed) at the frozen optimum of a TWO-SIDED level / yield / ratio / rate /
share row (pH-unit rows, kJ/mol rows and one-sided ceiling / floor rows are excluded); n_f the fit's
total objective rows, p_f its free coordinates. c_f is the degrees-of-freedom correction (residuals at
an optimum are shrunk by the p_f parameters spent on them: E[sum e^2] = sigma^2 (n - p)); it is the
parameter-uncertainty accounting at the fit rows, so the parameter share is not subtracted twice.

NO PANEL ROW IS READ. Every residual comes from a fit report's own row table; the panel's bundles,
measured values and the envelope artifact are never opened here.

Usage:
    python scripts/generators/generate_kinetic_core_b56_structural_discrepancy.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_b56_structural_discrepancy.json"
PREREG = "results/validation/kinetic_core_b55_b56_prereg.md"

#: Row kinds whose residual is not a two-sided log10 fold (sulfur row table vocabulary).
_EXCLUDED_SULFUR_KINDS = frozenset({"ph_endpoint", "ceiling", "peak_fraction_ceiling", "floor"})
#: Normal-consistency constant: sigma = 1.4826 * median |r| for a normal residual.
MAD_TO_SIGMA = 1.4826


def _read(name: str) -> Dict[str, Any]:
    return json.loads((V / name).read_text(encoding="utf-8"))


def sulfur_block() -> Dict[str, Any]:
    """The shipped sulfur wave's Laplace residuals (weighted) times each row's sigma_log."""
    import generate_kinetic_core_b2_3_fit as B23  # the row table the objective was built from
    import generate_kinetic_core_b8_fit  # noqa: F401  (installs B8's four rows into B23)

    from src.kinetic_core import engine

    lap_path = V / engine._B2_FIT_REPORT.name.replace("_fit_report.json", "_laplace_covariance.json")
    lap = json.loads(lap_path.read_text(encoding="utf-8"))
    rows = {row["id"]: row for row in B23.FIT_ROWS}
    used: List[Tuple[str, float]] = []
    for row_id, weighted in zip(lap["row_ids"], lap["residuals_at_optimum"], strict=True):
        row = rows[row_id]
        if row.get("kind") in _EXCLUDED_SULFUR_KINDS or row.get("sigma_log") is None:
            continue
        used.append((row_id, float(weighted) * float(row["sigma_log"])))
    return {
        "fit": data_paths.rel(lap_path),
        "n_rows_in_objective": int(lap["n_rows"]),
        "n_free": int(lap["n_free"]),
        "residuals_log10": dict(used),
        "how": "residuals_at_optimum x the row's sigma_log (B2.3/B8 row table); pH, ceiling and floor rows excluded",
    }


def trunk_blocks() -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    b1 = _read("kinetic_core_b1_fit_report.json")
    va = b1["variant_A"]
    groups = {}
    for g in va["per_response_fit"]:
        r = MAD_TO_SIGMA * math.log10(float(g["median_fold_error"]))
        groups[f"{g['species']}@{g['T_C']:g}C"] = {"n": int(g["n"]), "robust_sigma_log10": r}
    out.append({
        "fit": "results/validation/kinetic_core_b1_fit_report.json (variant_A, the variant every lane inherits)",
        "n_rows_in_objective": int(va["n_observations"]),
        "n_free": int(va["n_free_parameters"]),
        "grouped_robust": groups,
        "how": ("the report prints no row residuals: per (species, T) group, sum r^2 ~ n_g x "
                "(1.4826 x log10 median fold error)^2 (a declared robust approximation)"),
    })
    b7 = _read("kinetic_core_b7_fit_report.json")
    out.append({
        "fit": "results/validation/kinetic_core_b7_fit_report.json",
        "n_rows_in_objective": len(b7["fit"]["residual_log10_fold"]),
        "n_free": 1,
        "residuals_log10": {k: float(v) for k, v in b7["fit"]["residual_log10_fold"].items()},
        "how": "fit.residual_log10_fold",
    })
    b18 = _read("kinetic_core_b18_fit_report.json")
    out.append({
        "fit": "results/validation/kinetic_core_b18_fit_report.json",
        "n_rows_in_objective": int(b18["objective"]["n_rows"]),
        "n_free": int(b18["objective"]["n_free_parameters"]),
        "residuals_log10": {k: float(v) for k, v in b18["residual_by_row_dex"].items()},
        "how": "residual_by_row_dex",
    })
    b41 = _read("kinetic_core_b41_fit_report.json")
    rows41 = [r for r in b41["rows"] if r.get("decisive") and r.get("after_dex") is not None]
    out.append({
        "fit": "results/validation/kinetic_core_b41_fit_report.json",
        "n_rows_in_objective": int(b41["objective"]["n_rows"]),
        "n_free": int(b41["objective"]["n_free"]),
        "residuals_log10": {r["id"]: float(r["after_dex"]) for r in rows41},
        "how": "rows[].after_dex over the decisive rows",
    })
    return out


def acrylamide_block() -> Dict[str, Any]:
    b3 = _read("kinetic_core_b3_fit_report.json")
    used = {
        r["id"]: math.log10(float(r["fold_error"]))
        for r in b3["rows"]
        if r["kind"] in ("rate", "ppb") and r.get("fold_error") not in (None,) and float(r["fold_error"]) > 0
    }
    return {
        "fit": "results/validation/kinetic_core_b3_fit_report.json",
        "n_rows_in_objective": int(b3["objective"]["n_rows"]),
        "n_free": int(b3["objective"]["n_free_parameters"]),
        "residuals_log10": used,
        "how": "rows[] of kind rate / ppb: log10(fold_error); Ea rows (kJ/mol) and the ceiling row excluded",
    }


def lipid_block() -> Dict[str, Any]:
    b6 = _read("kinetic_core_b6_fit_report.json")
    used = {
        f"{r['system']}:{r['product']}": math.log10(
            float(r["predicted_relative_percent"]) / float(r["measured_relative_percent"]))
        for r in b6["fit_table"]
    }
    return {
        "fit": "results/validation/kinetic_core_b6_fit_report.json",
        "n_rows_in_objective": int(b6["fit_quality"]["n_fit_values"]),
        "n_free": int(b6["fit_quality"]["n_free_parameters"]),
        "residuals_log10": used,
        "how": "fit_table: log10(predicted / measured relative %)",
    }


def lane_sigma(blocks: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Pool the blocks with each fit's dof correction; record every term."""
    weighted_ss = 0.0
    n_used = 0
    raw_ss = 0.0
    for b in blocks:
        n, p = int(b["n_rows_in_objective"]), int(b["n_free"])
        if n - p <= 0:
            raise ValueError(f"{b['fit']}: n {n} <= p {p}; no residual degrees of freedom")
        c = n / (n - p)
        if "grouped_robust" in b:
            ss = sum(g["n"] * g["robust_sigma_log10"] ** 2 for g in b["grouped_robust"].values())
            k = sum(g["n"] for g in b["grouped_robust"].values())
        else:
            ss = sum(v * v for v in b["residuals_log10"].values())
            k = len(b["residuals_log10"])
        b["dof_correction"] = c
        b["n_used"] = k
        b["rms_log10_uncorrected"] = math.sqrt(ss / k) if k else None
        weighted_ss += c * ss
        raw_ss += ss
        n_used += k
    return {
        "sigma_log10": math.sqrt(weighted_ss / n_used) if n_used else 0.0,
        "rms_log10_uncorrected": math.sqrt(raw_ss / n_used) if n_used else 0.0,
        "n_used": n_used,
        "fits": blocks,
    }


def build() -> Dict[str, Any]:
    lanes = {
        "sulfur": lane_sigma([sulfur_block()]),
        "trunk": lane_sigma(trunk_blocks()),
        "acrylamide": lane_sigma([acrylamide_block()]),
        "lipid": lane_sigma([lipid_block()]),
    }
    inputs = sorted({ROOT / b["fit"].split(" ")[0] for lane in lanes.values() for b in lane["fits"]})
    return {
        "artifact": "kinetic_core_b56_structural_discrepancy",
        "wave": "B56",
        "prereg": PREREG,
        "provenance": provenance.provenance_block(
            "kinetic_core_b56_structural_discrepancy",
            generated_by="scripts/generators/generate_kinetic_core_b56_structural_discrepancy.py",
            inputs=tuple(inputs),
        ),
        "estimator": (
            "sigma_L^2 = sum_f c_f sum_i r_i^2 / sum_f n_used_f, c_f = n_f/(n_f - p_f); r_i the log10 "
            "residual at the frozen optimum of a two-sided row the fit read. A DECLARED MODEL-ERROR "
            "ASSUMPTION estimated in sample only, not a measurement; the parameter share is the dof "
            "correction and is not subtracted again."
        ),
        "sigma_log10_by_lane": {k: v["sigma_log10"] for k, v in lanes.items()},
        "lanes": lanes,
    }


def render(payload: Dict[str, Any]) -> str:
    out = [
        "# B56 structural-discrepancy term (in-sample fit residuals only)",
        "",
        f"Pre-registration: `{payload['prereg']}`. {payload['estimator']}",
        "",
        "| lane | sigma (log10) | RMS uncorrected | rows used | fits |",
        "|---|---:|---:|---:|---|",
    ]
    for lane, v in payload["lanes"].items():
        fits = "; ".join(
            f"{Path(b['fit'].split(' ')[0]).name} (n={b['n_rows_in_objective']}, p={b['n_free']}, used {b['n_used']}, "
            f"c={b['dof_correction']:.2f})"
            for b in v["fits"]
        )
        out.append(f"| {lane} | {v['sigma_log10']:.3f} | {v['rms_log10_uncorrected']:.3f} | {v['n_used']} | {fits} |")
    out.append("")
    return "\n".join(out)


def main() -> int:
    payload = build()
    json_path, md_path = artifact_io.write_artifact(payload, OUT, render=render, sort_keys=True)
    print(f"wrote {json_path}\nwrote {md_path}")
    for lane, s in payload["sigma_log10_by_lane"].items():
        print(f"  {lane:10s} sigma = {s:.3f} dex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

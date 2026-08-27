"""Derive the uncalibrated-tier matrix_headspace ln-sigma from in-panel residuals.

Context (audit 2026-08-26, tasks/audit_remediation.md): the external hold-out uses
DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS["matrix_headspace"], documented in
src/uncertainty_propagation.py as a structural-ignorance prior "NOT fitted to the
external hold-out values" — but the residual-based sizing promised in tasks/todo.md was
never done, so the value then shipped (ln-sigma 2.0) had no derivation. This script
computes that derivation. Its result (2.86) was adopted as the shipped value on
2026-08-26; the script now reads the live constant rather than restating a literal, so
the artifact can never drift from what the runtime actually uses.

Method — leave-lane-out transfer error, in-panel data only:
  In-panel matrix lanes are calibrated (their registry factors / marker yields absorb
  their own measurements), so their calibrated residuals are ~0 by construction and
  useless for sizing uncertainty. What the uncalibrated tier must describe is the error
  of predicting a lane the model has NOT calibrated. We emulate that in-panel: predict
  every non-reference matrix lane using ONLY the reference-lane calibration (the pea
  ambient marker yields, headspace factor 1.0) plus the raw oxidation kinetics, and
  score the ln residual against that lane's measured values.

  The external hold-out is untouched: bundles under data/benchmarks/external_validation/
  are structurally invisible to this script (non-recursive iteration over the panel
  directory), preserving the exclusion contract.

Caveats printed into the artifact:
  - n is small (one residual per lane x compound); the chi-square interval on sigma is
    reported alongside the point estimate.
  - RMS is taken around zero (bias included): the uncalibrated tier cannot assume a
    bias correction it does not have.

2026-08-27 (Wave O) — WHY THE REFIT DID NOT MOVE THIS ARTIFACT, and why that is correct.
  Wave O refitted the ambient-slurry hexanal observability factors against the
  content-verified Pratap-Singh anchors (pea 1.0 -> 4.31725, soy 2.2097561 -> 9.54007;
  results/validation/matrix_observability_refit_pratap_singh.json). Re-running this script
  afterwards produced a BIT-IDENTICAL artifact. That is structural, not luck:
  `_uncalibrated_prediction_ppb` below multiplies the oxidation load by the base MARKER
  YIELD only — it never reads `observable_factor` at all, because the whole point of the
  uncalibrated tier is to describe the error of a lane for which NO observability
  calibration exists. So this derivation is invariant to every constant in
  `src/matrix_calibration_registry.py` by construction, and no refit of those constants can
  ever be used to justify a narrower prior here.
  Two consequences worth stating out loud:
    * The shipped 2.86 remains inside the 90% CI [2.186, 6.796] around rms_ln_sigma 3.2520,
      exactly as it was after Wave M. Wave O supplies no new evidence about it either way.
    * The header's phrase "the reference-lane calibration (the pea ambient marker yields,
      headspace factor 1.0)" is now literally true only of the yields: the pea ambient lane
      does carry a non-unit hexanal observability factor since Wave O. The emulation is
      unchanged and still correct — "factor 1.0" here means "no lane-specific observability
      correction is available", which is the condition the uncalibrated tier describes.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    MATRIX_BENCHMARK_BASE_MARKER_YIELDS,
    MATRIX_BENCHMARK_PROFILES,
)
from src.lipid_oxidation import (
    hydroperoxide_pool_key_for_marker,
    predict_hexanal_generation,
)
from src.uncertainty_propagation import DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS

REFERENCE_LANE = ("pea_iso", "ambient_slurry")

_MEASURED_KEY_TO_YIELD_KEY = {
    "hexanal": "Hexanal",
    "2-pentylfuran": "2-Pentylfuran",
    "hexanol": "1-Hexanol",
    "1-hexanol": "1-Hexanol",
    "nonanal": "Nonanal",
}

# chi-square 5%/95% quantiles for n degrees of freedom, for the CI on sigma
_CHI2_BOUNDS = {6: (1.635, 12.592), 5: (1.145, 11.070), 4: (0.711, 9.488), 3: (0.352, 7.815)}


def _uncalibrated_prediction_ppb(protein_type: str, temp_c: float, time_min: float, yield_key: str) -> float:
    profile = MATRIX_BENCHMARK_PROFILES[protein_type]["lipid_profile"]
    oxidation = predict_hexanal_generation(
        profile, temp_C=float(temp_c), time_min=float(time_min), oxygen_availability=1.0
    )
    # Wave P item 4 (2026-08-27): the uncalibrated load is PER MARKER now, because
    # nonanal is cleaved from the OLEATE pool and the rest from the LINOLEATE pool.
    # This is NOT a no-op: one of the five scored residuals is the Trikusuma heated-pea
    # Nonanal row, whose uncalibrated prediction falls 3238.93 -> 1425.13 ppb (exactly
    # the pea oleic/linoleic ratio, 0.4400x). Measured effect on the derived prior:
    # rms_ln_sigma 3.2520 -> 3.0166, bias_fold 3.9065 -> 3.3149, 90% CI
    # [2.1856, 6.7958] -> [2.0274, 6.3038]. The SHIPPED 2.86 was NOT moved and is
    # still inside the CI -- it moved closer to the derivation, which is a consequence
    # of a substrate correction and not a reason to narrow the prior.
    load_ppb = float(oxidation[hydroperoxide_pool_key_for_marker(yield_key)]) * 1000.0
    return load_ppb * float(MATRIX_BENCHMARK_BASE_MARKER_YIELDS[yield_key])


def build_payload() -> dict:
    residual_rows = []
    for bench_file in sorted((ROOT / "data" / "benchmarks").glob("*.json")):
        bench = json.loads(bench_file.read_text())
        if (bench.get("metadata") or {}).get("execution_path") != "matrix_only":
            continue
        protein_type = str(bench.get("protein_type"))
        process_state = str((bench.get("process_metadata") or {}).get("state"))
        if (protein_type, process_state) == REFERENCE_LANE:
            continue  # the reference lane defines the yields; its residual is 0 by construction
        conditions = bench.get("conditions", {})
        for measured_key, entry in (bench.get("measured_volatiles") or {}).items():
            yield_key = _MEASURED_KEY_TO_YIELD_KEY.get(measured_key.strip().lower())
            measured = float((entry or {}).get("conc_ppb", 0.0) or 0.0)
            if yield_key is None or measured <= 0.0:
                continue
            predicted = _uncalibrated_prediction_ppb(
                protein_type, conditions["temp_C"], conditions["time_min"], yield_key
            )
            ln_resid = math.log(predicted / measured)
            residual_rows.append(
                {
                    "benchmark_id": bench.get("benchmark_id", bench_file.stem),
                    "lane": f"{protein_type}/{process_state}",
                    "compound": yield_key,
                    "measured_ppb": measured,
                    "uncalibrated_predicted_ppb": predicted,
                    "ln_residual": ln_resid,
                    "fold_error": math.exp(abs(ln_resid)),
                }
            )

    n = len(residual_rows)
    ln_residuals = [row["ln_residual"] for row in residual_rows]
    mean_ln = sum(ln_residuals) / n if n else None
    rms_ln = math.sqrt(sum(r * r for r in ln_residuals) / n) if n else None
    centered_sd = (
        math.sqrt(sum((r - mean_ln) ** 2 for r in ln_residuals) / (n - 1)) if n > 1 else None
    )
    lo_q, hi_q = _CHI2_BOUNDS.get(n, (None, None))
    ci = (
        {
            "sigma_low_90": math.sqrt(n * rms_ln * rms_ln / hi_q),
            "sigma_high_90": math.sqrt(n * rms_ln * rms_ln / lo_q),
        }
        if rms_ln is not None and lo_q
        else {}
    )

    # Read the LIVE shipped sigma instead of restating a literal. It was raised
    # 2.0 -> 2.86 on 2026-08-26 (owner-approved, residual-derived point estimate)
    # and this artifact silently kept asserting 2.0 until 2026-08-27.
    shipped = float(DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS["matrix_headspace"])
    return {
        "schema_version": "1.0",
        "description": (
            "Leave-lane-out residual derivation of the uncalibrated-tier matrix_headspace "
            "ln-sigma, using in-panel matrix anchors only (hold-out untouched)."
        ),
        "method": {
            "reference_lane": "/".join(REFERENCE_LANE),
            "estimator": "RMS of ln(pred/meas) around zero (bias included)",
            "holdout_exclusion": "data/benchmarks/external_validation/ never read by this script",
        },
        "residuals": residual_rows,
        "summary": {
            "n_residuals": n,
            "mean_ln_residual": mean_ln,
            "bias_fold": math.exp(mean_ln) if mean_ln is not None else None,
            "rms_ln_sigma": rms_ln,
            "centered_sd_ln": centered_sd,
            "sigma_90pct_ci": ci,
            "shipped_uncalibrated_sigma": shipped,
            "shipped_within_ci": (
                bool(ci and ci["sigma_low_90"] <= shipped <= ci["sigma_high_90"]) if ci else None
            ),
        },
    }


def render_markdown(payload: dict) -> str:
    s = payload["summary"]
    lines = [
        "<!-- Auto-generated by scripts/generators/derive_matrix_sigma_from_residuals.py. Manual edits will be overwritten. -->",
        "",
        "# Uncalibrated-Tier Matrix Sigma: Residual Derivation",
        "",
        "_Leave-lane-out transfer error computed from in-panel matrix anchors only; the",
        "external hold-out is structurally excluded. This is the residual-based sizing",
        f"promised in tasks/todo.md for the `uncalibrated` prior tier "
        f"(shipped ln-sigma {s['shipped_uncalibrated_sigma']:.2f})._",
        "",
        "| Lane | Compound | Measured ppb | Uncalibrated pred. ppb | ln residual | fold |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["residuals"]:
        lines.append(
            f"| {row['lane']} | {row['compound']} | {row['measured_ppb']:.6g} | "
            f"{row['uncalibrated_predicted_ppb']:.6g} | {row['ln_residual']:+.3f} | {row['fold_error']:.2f}x |"
        )
    ci = s.get("sigma_90pct_ci") or {}
    lines += [
        "",
        f"**n = {s['n_residuals']}** residuals · mean ln residual {s['mean_ln_residual']:+.3f} "
        f"(bias {s['bias_fold']:.2f}x) · **RMS ln-sigma = {s['rms_ln_sigma']:.3f}** "
        f"(centered SD {s['centered_sd_ln']:.3f})",
        "",
        (
            f"90% chi-square interval on sigma: [{ci.get('sigma_low_90', float('nan')):.2f}, "
            f"{ci.get('sigma_high_90', float('nan')):.2f}] — shipped value "f"{s['shipped_uncalibrated_sigma']:.2f} is "
            f"{'INSIDE' if s['shipped_within_ci'] else 'OUTSIDE'} this interval."
            if ci
            else "Too few residuals for a chi-square interval."
        ),
        "",
        "> **Caveats:** n is very small and the lanes are not independent draws (two share a",
        "> paper; one shares a matrix). RMS around zero includes systematic bias — the",
        "> uncalibrated tier has no bias correction to apply, so its sigma must cover bias.",
        "> This derivation sizes transfer to *moderate* new lanes; extreme processing",
        "> (roasting, HME) sits outside even this envelope, as the hold-out shows.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    payload = build_payload()
    out_dir = ROOT / "results" / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "matrix_sigma_residual_derivation.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    (out_dir / "matrix_sigma_residual_derivation.md").write_text(render_markdown(payload))
    s = payload["summary"]
    print(
        f"n={s['n_residuals']} rms_ln_sigma={s['rms_ln_sigma']:.3f} "
        f"bias_fold={s['bias_fold']:.2f} shipped={s['shipped_uncalibrated_sigma']:.2f} "f"within_ci={s['shipped_within_ci']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

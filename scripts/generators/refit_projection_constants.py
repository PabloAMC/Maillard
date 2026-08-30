#!/usr/bin/env python3
"""Refit the free constants of the volatile-projection budget.

Context (2026-08-27 audit remediation, Part 1)
---------------------------------------------
`src/projection.py` v1 scaled the volatile budget by a severity sigmoid centred on
110 C with width 18 C. That sigmoid saturates by construction: it reaches 0.966 at
170 C and 0.988 at 190 C, so the budget grew only ~1.11x from 150 C to 190 C while
the model's own Arrhenius drive grows roughly 20x over the same span. Dividing a
nearly-constant budget by a still-growing softmax denominator manufactured a
spurious furfural maximum near 145-150 C and capped all high-temperature chemistry.

v2 replaces the budget's thermal dependence with a first-order conversion extent
under an apparent Arrhenius rate:

    drive(T, t) = (t / tau_ref) * exp(-(Ea_app/R) * (1/T - 1/T_ref))
    yield       = baseline + (ceiling - baseline) * (1 - exp(-drive))

with `ceiling` pinned at 1.0 (mass conservation on the limiting precursor pool),
`Ea_app` pinned at 120 kJ/mol (literature; see the comment in src/projection.py),
`T_ref` = 423.15 K. That leaves TWO free constants, refit here:

    * baseline_volatile_yield_fraction
    * reference_conversion_time_min  (tau_ref)

Fit targets
-----------
LITERATURE-SOURCED benchmark rows ONLY. Row eligibility uses the repo's own
classifier, `src.uncertainty_propagation._benchmark_signal_origin`, and keeps only
`external_literature`. Synthetic rows (Internal2026 / ProtocolPilot2026) and
everything under `data/benchmarks/external_validation/` are FORBIDDEN as fit
targets -- the hold-out guard and the evidence_class conventions are audit-critical
and this script must never quietly launder them into the calibration set. Both
exclusions are asserted below, not merely intended.

Objective
---------
    J = mean over eligible matched rows of |log10(predicted_ppb / measured_ppb)|

Rows that fail to match, or that predict/measure zero, are not silently dropped
(that would let the optimiser win by zeroing predictions): they carry a fixed
penalty of `UNSCORED_ROW_PENALTY` dex so the count of scored rows stays part of
the objective.

Usage
-----
    python scripts/generators/refit_projection_constants.py            # full refit
    python scripts/generators/refit_projection_constants.py --quick    # coarse only
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import projection as projection_module
from src import recommend as recommend_module
from src.benchmark_validation import evaluate_benchmark, get_benchmark_files
from src.projection import DEFAULT_PROJECTION_STRATEGY
from src.uncertainty_propagation import _benchmark_signal_origin

UNSCORED_ROW_PENALTY = 3.0  # dex; ~the worst genuine residual in the panel
BASELINE_MOVE_THRESHOLD_DEX = 0.01  # minimum objective gain to move off the incumbent


def _install_strategy(strategy) -> None:
    """Point every live consumer of DEFAULT_PROJECTION_STRATEGY at `strategy`."""
    projection_module.DEFAULT_PROJECTION_STRATEGY = strategy
    recommend_module.DEFAULT_PROJECTION_STRATEGY = strategy


def _literature_benchmark_files() -> List[Path]:
    files: List[Path] = []
    for path in get_benchmark_files():
        # Hard guard: the hold-out set must never become a fit target.
        assert "external_validation" not in path.parts, (
            f"external_validation benchmark reached the fit-target selector: {path}"
        )
        if _benchmark_signal_origin(path) != "external_literature":
            continue
        files.append(path)
    assert files, "no literature-sourced benchmarks found"
    for path in files:
        # Hard guard: the two synthetic lanes must never become fit targets.
        assert "Internal2026" not in path.name and "ProtocolPilot2026" not in path.name, (
            f"synthetic benchmark reached the fit-target selector: {path}"
        )
    return files


def _score(files: List[Path]) -> Tuple[float, List[Dict[str, Any]]]:
    rows: List[Dict[str, Any]] = []
    for path in files:
        evaluation = evaluate_benchmark(path)
        if not evaluation.supported:
            continue
        for comparison in evaluation.comparisons:
            scored = (
                comparison.matched_name is not None
                and comparison.measured_ppb > 0.0
                and comparison.predicted_ppb > 0.0
            )
            error = (
                abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb))
                if scored
                else UNSCORED_ROW_PENALTY
            )
            rows.append(
                {
                    "benchmark_id": evaluation.benchmark_id,
                    "compound": comparison.compound,
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "scored": scored,
                    "abs_log10_error": error,
                }
            )
    if not rows:
        return float("inf"), rows
    return sum(row["abs_log10_error"] for row in rows) / len(rows), rows


def _evaluate_point(files, *, baseline: float, tau_ref: float, ea_kj: float) -> Tuple[float, List[Dict[str, Any]]]:
    strategy = replace(
        DEFAULT_PROJECTION_STRATEGY,
        baseline_volatile_yield_fraction=baseline,
        reference_conversion_time_min=tau_ref,
        apparent_activation_energy_kj_mol=ea_kj,
    )
    _install_strategy(strategy)
    try:
        return _score(files)
    finally:
        _install_strategy(DEFAULT_PROJECTION_STRATEGY)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--quick", action="store_true", help="coarse scan only")
    args = parser.parse_args()

    files = _literature_benchmark_files()
    print(f"Fit targets: {len(files)} literature-sourced benchmarks")
    for path in files:
        print(f"  - {path.name}")

    ea_kj = DEFAULT_PROJECTION_STRATEGY.apparent_activation_energy_kj_mol
    baseline0 = DEFAULT_PROJECTION_STRATEGY.baseline_volatile_yield_fraction

    record: Dict[str, Any] = {
        "generated_by": "scripts/generators/refit_projection_constants.py",
        "objective": "mean |log10(predicted_ppb / measured_ppb)| over literature-sourced benchmark rows",
        "unscored_row_penalty_dex": UNSCORED_ROW_PENALTY,
        "fit_target_files": [path.name for path in files],
        "forbidden_as_fit_targets": [
            "data/benchmarks/external_validation/** (hold-out)",
            "*Internal2026*  (synthetic reproducibility lane)",
            "*ProtocolPilot2026*  (synthetic reproducibility lane)",
        ],
        # 2026-08-27 (Wave I). Disclosed contamination of THIS record's own history.
        "contamination_note": {
            "date": "2026-08-27",
            "wave": "Wave I (cold-start red-team remediation)",
            "removed_fit_targets": [
                "spi_hvp_xylose_120C_PMC9905368.json",
                "wheat_gluten_hvp_xylose_120C_PMC9905368.json",
            ],
            "why": (
                "Both files were fabricated. The cited paper 10.1007/s10068-022-01194-w "
                "reacts protein hydrolysates with GLUCOSE and FRUCTOSE at pH 7.5 for 90 min "
                "and reports only RELATIVE GC-MS peak areas; it never mentions "
                "2-furfurylthiol or 2-methyl-3-furanthiol and reports no absolute "
                "concentration for any analyte. They are now quarantined and this glob no "
                "longer sees them."
            ),
            "affected_prior_versions": (
                "Every projection_constant_refit record written before 2026-08-27 included "
                "6 rows from those two files (3 analytes each) in the fit objective, i.e. "
                "6 of 15 literature rows. The fitted constant of record, "
                "reference_conversion_time_min = 1.2589e4, was therefore fitted against a "
                "population containing fabricated data."
            ),
            "what_survives": (
                "The QUALITATIVE conclusions of the Wave H projection analysis survive and "
                "are unaffected: (a) tau_ref is a single GLOBAL SCALE, so it cannot fix an "
                "ALLOCATION deficit; (b) the deliberate decision NOT to apply the re-fitted "
                "optimum stands, and stands for a reason (Hofmann MFT vs Resconi furfural "
                "trade-off) that does not involve these files; (c) the Boltzmann "
                "de-duplication finding is a structural result, not a fit. What does NOT "
                "survive is the numeric value of the objective in earlier records and any "
                "reading of the fitted tau_ref as literature-anchored across six of its rows. "
                "The re-run that produced the current record excludes the quarantined files."
            ),
        },
        "fixed_constants": {
            "apparent_activation_energy_kj_mol": ea_kj,
            "reference_temperature_kelvin": DEFAULT_PROJECTION_STRATEGY.reference_temperature_kelvin,
            "conversion_ceiling_fraction": DEFAULT_PROJECTION_STRATEGY.conversion_ceiling_fraction,
        },
        "free_constants": ["baseline_volatile_yield_fraction", "reference_conversion_time_min"],
        # 2026-08-27 (Wave I). Required by scripts/ci/fit_target_gate.py: a fit record must
        # state how much freedom it had, so a reader can tell "one global scale across the
        # whole panel" from "one knob per row". Filled in below once the row count is known.
        "fit_leverage": None,
        "deliberately_not_refit": {
            "recommend.py depth_bias_strength (0.85 offset, 1.0 slope)": "unconstrained legacy fit",
            "recommend.py direct_sulfur_bonus (0.8 coefficient)": "unconstrained legacy fit",
            "rationale": (
                "These two allocation heuristics are free constants of the same layer and "
                "were fitted in an era whose targets (Mottram1994, Farmer1999) are now "
                "quarantined, so they are re-anchoring candidates. They are deliberately "
                "NOT refit here: direct_sulfur_bonus in particular is the exact knob that "
                "would absorb the MFT residuals this panel is supposed to expose, and "
                "tuning it now would launder the Boltzmann finding instead of reporting it. "
                "Flagged for a dedicated workstream."
            ),
        },
        # Measured on this same objective against the pre-retune tree
        # (severity-sigmoid budget + double Boltzmann), 2026-08-27. Recorded here so
        # the cost of the retune stays visible without needing a git checkout.
        "v1_reference": {
            "description": (
                "v1: severity-sigmoid volatile budget, Boltzmann selectivity applied twice "
                "(effective temperature T/2.19)"
            ),
            "objective_all_literature_rows": 0.1469,
            "objective_projection_sensitive_rows_only": 0.3089,
            "objective_projection_sensitive_rows_excluding_bolton": 0.1186,
            "note": (
                "The v1 objective is BETTER. Essentially all of the gap is the "
                "de-duplication of the Boltzmann factor, not the budget retune: the "
                "budget retune is a single global scale that the tau_ref fit absorbs, "
                "whereas removing the unphysical T/2.19 sharpening changes how the "
                "budget is split between competing species. In other words, ~0.27 dex "
                "of the panel's previous agreement on multi-analyte free-precursor "
                "benchmarks was being carried by a selectivity term evaluated at less "
                "than half the physical temperature. That is a finding, not a "
                "regression to be tuned away."
            ),
        },
    }

    # ---- Stage 1: coarse 1-D scan over tau_ref (decades) -----------------
    coarse: List[Dict[str, float]] = []
    best = (float("inf"), None)
    for exponent in range(1, 12):
        tau = 10.0 ** exponent
        value, _ = _evaluate_point(files, baseline=baseline0, tau_ref=tau, ea_kj=ea_kj)
        coarse.append({"log10_tau_ref_min": float(exponent), "objective": value})
        print(f"  stage1 tau_ref=1e{exponent:<2d}  J={value:.4f}")
        if value < best[0]:
            best = (value, exponent)
    record["stage1_coarse_tau_scan"] = coarse
    best_exponent = best[1]

    # ---- Stage 2: fine 1-D scan around the coarse optimum ----------------
    fine: List[Dict[str, float]] = []
    best_tau_log = float(best_exponent)
    if not args.quick:
        candidates = [best_exponent - 1.0 + 0.1 * step for step in range(0, 21)]
        best = (float("inf"), best_tau_log)
        for log_tau in candidates:
            tau = 10.0 ** log_tau
            value, _ = _evaluate_point(files, baseline=baseline0, tau_ref=tau, ea_kj=ea_kj)
            fine.append({"log10_tau_ref_min": round(log_tau, 3), "objective": value})
            if value < best[0]:
                best = (value, log_tau)
        best_tau_log = best[1]
        print(f"  stage2 optimum log10(tau_ref) = {best_tau_log:.2f}  J={best[0]:.4f}")
    record["stage2_fine_tau_scan"] = fine

    tau_star = 10.0 ** best_tau_log

    # ---- Stage 3: baseline identifiability at the tau optimum ------------
    # The incumbent is the reference point and wins ties: the baseline is a floor
    # term that only matters where the Arrhenius drive is negligible, so it is
    # essentially unidentified by this panel. Moving it on a tie would be noise.
    baseline_scan: List[Dict[str, float]] = []
    best_baseline = baseline0
    best_baseline_obj, _ = _evaluate_point(files, baseline=baseline0, tau_ref=tau_star, ea_kj=ea_kj)
    for exponent in (-9, -8, -7, -6, -5, -4):
        baseline = 10.0 ** exponent
        value, _ = _evaluate_point(files, baseline=baseline, tau_ref=tau_star, ea_kj=ea_kj)
        baseline_scan.append({"log10_baseline": float(exponent), "objective": value})
        print(f"  stage3 baseline=1e{exponent:<3d} J={value:.4f}")
        # Only move off the incumbent for a materially better objective. Anything
        # under a hundredth of a dex on a 12-benchmark panel is noise on a flat
        # direction, and chasing it would be overfitting.
        if value < best_baseline_obj - BASELINE_MOVE_THRESHOLD_DEX:
            best_baseline_obj = value
            best_baseline = baseline
    record["stage3_baseline_scan"] = baseline_scan

    # ---- Stage 4: Ea sensitivity (reported, NOT fitted) ------------------
    ea_scan: List[Dict[str, float]] = []
    if not args.quick:
        for ea in (97.0, 110.0, 120.0, 138.0, 160.0):
            value, _ = _evaluate_point(files, baseline=best_baseline, tau_ref=tau_star, ea_kj=ea)
            ea_scan.append({"apparent_activation_energy_kj_mol": ea, "objective": value})
            print(f"  stage4 Ea={ea:5.1f} kJ/mol  J={value:.4f}")
    record["stage4_ea_sensitivity_not_fitted"] = ea_scan

    # ---- Final: install the fit and record per-row residuals -------------
    final_obj, final_rows = _evaluate_point(
        files, baseline=best_baseline, tau_ref=tau_star, ea_kj=ea_kj
    )
    # Which rows does the projection layer actually move? Half the panel runs on the
    # matrix-only / safety lanes and never touches the volatile budget, so the
    # all-rows objective understates how much leverage the fit target really has.
    _, probe_rows = _evaluate_point(
        files, baseline=best_baseline, tau_ref=tau_star * 1.05, ea_kj=ea_kj
    )
    probe_by_key = {(row["benchmark_id"], row["compound"]): row for row in probe_rows}
    for row in final_rows:
        probe = probe_by_key.get((row["benchmark_id"], row["compound"]))
        row["projection_sensitive"] = bool(
            probe is not None
            and abs(probe["predicted_ppb"] - row["predicted_ppb"])
            > 1.0e-9 * max(abs(row["predicted_ppb"]), 1.0)
        )
    sensitive = [row for row in final_rows if row["projection_sensitive"]]

    # ---- Shipped-vs-fitted decision record -------------------------------
    # Added 2026-08-27 (Wave H). This script has never written to src/projection.py;
    # a human applies the fit. After the Wave G1 chemistry rebuild the two are no
    # longer the same value, so the script now states what is shipped, what the fit
    # wants, and what the gap actually buys — otherwise the artifact reads as if the
    # runtime were tracking it.
    shipped_tau = DEFAULT_PROJECTION_STRATEGY.reference_conversion_time_min
    shipped_baseline = DEFAULT_PROJECTION_STRATEGY.baseline_volatile_yield_fraction
    shipped_obj, shipped_rows = _evaluate_point(
        files, baseline=shipped_baseline, tau_ref=shipped_tau, ea_kj=ea_kj
    )
    record["shipped_constants"] = {
        "baseline_volatile_yield_fraction": shipped_baseline,
        "reference_conversion_time_min": shipped_tau,
        "objective": shipped_obj,
        "rows": shipped_rows,
    }
    record["applied_to_runtime"] = False
    record["fitted_constants"] = {
        "baseline_volatile_yield_fraction": best_baseline,
        "reference_conversion_time_min": tau_star,
    }
    record["objective_at_fit"] = final_obj
    record["objective_projection_sensitive_rows_only"] = (
        sum(row["abs_log10_error"] for row in sensitive) / len(sensitive) if sensitive else None
    )
    record["projection_sensitive_rows"] = len(sensitive)
    record["scored_rows"] = sum(1 for row in final_rows if row["scored"])
    record["total_rows"] = len(final_rows)
    record["rows"] = final_rows

    # 2026-08-27 (Wave I). Fit leverage — see src/fit_target_index.py.
    _scored = record["scored_rows"]
    _free = len(record["free_constants"])
    record["fit_leverage"] = {
        "free_parameters": _free,
        "fitted_rows": _scored,
        "parameters_per_row": (_free / _scored) if _scored else None,
        "class": "global_low_leverage" if _scored and (_free / _scored) < 0.5 else "per_row_recovery",
        "interpretation": (
            f"{_free} GLOBAL constants fitted across {_scored} rows. Two constants cannot "
            "make that many rows agree, so agreement on any individual row still carries "
            "information and these rows STAY in the literature-coverage denominator. They "
            "are annotated `fit_target_of` in prediction_uncertainty.json so a reader knows "
            "they are not strictly out-of-sample. Contrast the matrix observability "
            "factors, which are one free factor PER ROW and therefore reproduce their "
            "targets exactly -- those are excluded from the count entirely."
        ),
        "honest_caveat": (
            "This does mean there is NO literature benchmark in the panel that is fully "
            "out-of-sample with respect to the projection scale. The only genuinely "
            "out-of-sample evidence in the repo is the external_validation/ hold-out set."
        ),
    }

    # 2026-08-27 (Wave M). The two residuals this decision turns on used to be
    # HARD-CODED in the prose below, so the artifact kept quoting Wave H numbers
    # after both the Wave N route correction and the Wave K/M content corrections
    # had moved them. They are now read out of the rows this run just computed.
    def _dex(rows: List[Dict[str, Any]], bench: str, compound_sub: str) -> float | None:
        for row in rows:
            if row["benchmark_id"] == bench and compound_sub.lower() in row["compound"].lower():
                return float(row["abs_log10_error"])
        return None

    _mft_shipped = _dex(shipped_rows, "cys_ribose_140C_Hofmann1998", "2-methyl-3-furanthiol")
    _mft_fitted = _dex(final_rows, "cys_ribose_140C_Hofmann1998", "2-methyl-3-furanthiol")
    _furfural_shipped = _dex(shipped_rows, "resconi_2023_pbma_beef_identity_benchmark", "furfural")
    _furfural_fitted = _dex(final_rows, "resconi_2023_pbma_beef_identity_benchmark", "furfural")
    _budget_scale = shipped_tau / tau_star if tau_star else float("nan")

    # 2026-08-27 (Wave P) — SECOND PASS OF THE WAVE M FIX. Wave M derived the two
    # RESIDUALS this prose quotes but left the ALLOCATION arithmetic (~1126 ppb,
    # 542 ppb, ~51 %, ~2800 ppb) and the MFT/FFT fold errors hard-coded, so they went
    # stale again as soon as the Wave P chemistry and the barrier refit moved them.
    # They are now computed from the live prediction, every time. This is the pattern
    # flagged as [P] in Wave M: a generator that writes narrative next to numbers must
    # derive ALL of them or it publishes fiction with a timestamp on it.
    def _hofmann_allocation() -> Dict[str, Any]:
        from src.benchmark_validation import evaluate_benchmark as _eval

        bench = ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json"
        evaluation = _eval(bench)
        # `predicted_ppb` carries each species under several aliases (SMILES, lower-case
        # label, display label); select one canonical key per species so the sum is not
        # multiplied by the alias count.
        canonical = [
            "furfural",
            "2-methyl-3-furanthiol",
            "2-furfurylthiol",
            "bis(2-methyl-3-furyl) disulfide",
            "2,5-dimethylpyrazine",
        ]
        per_species = {
            name: float(evaluation.predicted_ppb.get(name, 0.0)) for name in canonical
        }
        total = sum(per_species.values())
        measured = {
            c.compound: float(c.measured_ppb) for c in evaluation.comparisons
        }
        measured_total = sum(measured.values())
        ratios = {
            c.compound: (float(c.predicted_ppb) / float(c.measured_ppb))
            for c in evaluation.comparisons
            if c.measured_ppb
        }
        return {
            "per_species_ppb": per_species,
            "tracked_total_ppb": total,
            "measured_total_ppb": measured_total,
            "furfural_share": (per_species["furfural"] / total) if total else float("nan"),
            "ratios": ratios,
        }

    _alloc = _hofmann_allocation()

    def _fold(name: str) -> str:
        ratio = _alloc["ratios"].get(name)
        if ratio is None or ratio <= 0.0:
            return "n/a"
        return (
            f"{ratio:.2f}x over" if ratio >= 1.0 else f"{1.0 / ratio:.2f}x under"
        )

    def _f(value: float | None, fmt: str) -> str:
        return format(value, fmt) if value is not None else "n/a"

    record["wave_h_decision"] = {
        "date": "2026-08-27",
        "decision": "FIT NOT APPLIED — incumbent reference_conversion_time_min retained",
        "gap": (
            f"shipped tau_ref = {shipped_tau:.6g} min (objective {shipped_obj:.4f} dex); "
            f"refit optimum = {tau_star:.6g} min (objective {final_obj:.4f} dex), i.e. the "
            f"fit wants the global volatile budget scaled by {_budget_scale:.2f}x. The "
            "optimum has moved repeatedly as the chemistry and the reference data were "
            "corrected (Wave G1 chemistry rebuild, Wave I quarantines, Wave N MFT route "
            "correction, Wave K/M benchmark content corrections); it is a scale that "
            "absorbs whatever the rest of the model gets wrong, which is exactly why it "
            "is reported and not applied."
        ),
        "why_not_applied": (
            "tau_ref is a SINGLE GLOBAL SCALE on the volatile budget. Moving it does not "
            "improve any mechanism; it trades one lane against another, and this panel's "
            "mean-|log10| objective takes that trade happily. At the refit optimum the "
            f"Hofmann1998 MFT residual collapses to {_f(_mft_fitted, '.4f')} dex (from "
            f"{_f(_mft_shipped, '.4f')} dex shipped) while the Resconi 2023 furfural "
            f"residual grows to {_f(_furfural_fitted, '.4f')} dex (from "
            f"{_f(_furfural_shipped, '.4f')} dex shipped), i.e. from "
            f"{_f(10.0 ** _furfural_shipped if _furfural_shipped is not None else None, '.1f')}x "
            f"to {_f(10.0 ** _furfural_fitted if _furfural_fitted is not None else None, '.1f')}x "
            "OVER. The arithmetic is decisive: at the shipped value the total predicted "
            f"pool over the five species this benchmark's system tracks is "
            f"{_alloc['tracked_total_ppb']:.1f} ppb against a measured MFT+FFT of "
            f"{_alloc['measured_total_ppb']:.0f} ppb (derived live, "
            + ", ".join(
                f"{name} {value:.1f}"
                for name, value in _alloc["per_species_ppb"].items()
            )
            + "), so the budget is already the right order of magnitude and the sulfur "
            "deficit is an ALLOCATION deficit — furfural, unmeasured in that benchmark, "
            f"takes {_alloc['furfural_share'] * 100.0:.0f}% of the pool. See "
            "results/validation/sulfur_barrier_refit_pentodiulose.md (the Wave P refit of "
            "`thiol_addition_pentodiulose` against this same benchmark; it SUPERSEDES "
            "sulfur_barrier_refit_hofmann.md, which is stale because it profiles "
            "`thiol_addition_norfuraneol`, a family the isotope evidence retired). "
            f"Raising the global budget {_budget_scale:.2f}x to close the current MFT gap "
            f"({_fold('2-methyl-3-furanthiol')}) means supplying "
            f"~{_alloc['tracked_total_ppb'] * _budget_scale:.0f} ppb of volatiles at these "
            f"conditions to explain {_alloc['measured_total_ppb']:.0f} ppb of measured "
            "product, sending the rest to species nobody measured, and pushing FFT from "
            f"{_fold('2-furfurylthiol')} further over. That is not a calibration, it is "
            "a way of making the allocation finding invisible. The incumbent is therefore "
            "retained and the gap is reported."
        ),
        "what_would_justify_applying_it": (
            "A refit of the ALLOCATION (the depth-bias and direct-sulfur heuristics, both "
            "unconstrained legacy fits from the quarantined-target era), done first and "
            "with its own literature constraints, after which a global scale would be "
            "fitting a scale rather than absorbing a selectivity error."
        ),
    }

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "projection_constant_refit.json"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    lines = [
        "# Projection constant refit (v2 Arrhenius budget)",
        "",
        f"Objective: `{record['objective']}`",
        "",
        f"Fit targets: {len(files)} literature-sourced benchmarks "
        f"({record['scored_rows']}/{record['total_rows']} rows scored).",
        "",
        "Synthetic (Internal2026 / ProtocolPilot2026) rows and the "
        "`external_validation/` hold-out are excluded by assertion, not convention.",
        "",
        "## Shipped vs fitted (2026-08-27, Wave H)",
        "",
        "**This fit is NOT applied to the runtime.** The script has never written to "
        "`src/projection.py`; after the Wave G1 chemistry rebuild the shipped value and the "
        "refit optimum are no longer the same number, so both are stated here.",
        "",
        "| Constant | Shipped (runtime) | Refit optimum |",
        "| --- | --- | --- |",
        f"| baseline_volatile_yield_fraction | {shipped_baseline:.3e} | {best_baseline:.3e} |",
        f"| reference_conversion_time_min | {shipped_tau:.4g} | {tau_star:.4g} |",
        f"| objective (all literature rows) | {shipped_obj:.4f} dex | {final_obj:.4f} dex |",
        "",
        f"**Decision: {record['wave_h_decision']['decision']}.**",
        "",
        record["wave_h_decision"]["why_not_applied"],
        "",
        "What would justify applying it: "
        + record["wave_h_decision"]["what_would_justify_applying_it"],
        "",
        "## Fitted constants",
        "",
        "| Constant | Value |",
        "| --- | --- |",
        f"| baseline_volatile_yield_fraction | {best_baseline:.3e} |",
        f"| reference_conversion_time_min | {tau_star:.4g} |",
        "",
        "## Fixed (not fitted)",
        "",
        "| Constant | Value | Basis |",
        "| --- | --- | --- |",
        f"| apparent_activation_energy_kj_mol | {ea_kj:.1f} | data/lit/arrhenius_params.yml `enolisation` |",
        f"| reference_temperature_kelvin | {DEFAULT_PROJECTION_STRATEGY.reference_temperature_kelvin:.2f} | 150 C panel anchor |",
        f"| conversion_ceiling_fraction | {DEFAULT_PROJECTION_STRATEGY.conversion_ceiling_fraction:.1f} | mass conservation |",
        "",
        f"Objective at fit: **{final_obj:.4f} dex** over all "
        f"{record['total_rows']} literature rows; "
        f"**{record['objective_projection_sensitive_rows_only']:.4f} dex** over the "
        f"{len(sensitive)} rows the projection layer actually moves (the rest run on the "
        "matrix-only and safety lanes, which never touch the volatile budget).",
        "",
        "## Versus v1 (severity sigmoid + double Boltzmann)",
        "",
        "| Objective | v1 | v2 (this fit) |",
        "| --- | --- | --- |",
        f"| All literature rows | {record['v1_reference']['objective_all_literature_rows']:.4f} | {final_obj:.4f} |",
        f"| Projection-sensitive rows | {record['v1_reference']['objective_projection_sensitive_rows_only']:.4f} | "
        f"{record['objective_projection_sensitive_rows_only']:.4f} |",
        "",
        record["v1_reference"]["note"],
        "",
        "## Per-row residuals at the fit",
        "",
        "| Benchmark | Compound | Measured ppb | Predicted ppb | \\|log10 ratio\\| | Projection-sensitive |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in sorted(final_rows, key=lambda item: -item["abs_log10_error"]):
        lines.append(
            f"| {row['benchmark_id']} | {row['compound']} | {row['measured_ppb']:.4g} | "
            f"{row['predicted_ppb']:.4g} | {row['abs_log10_error']:.3f} | "
            f"{'yes' if row['projection_sensitive'] else 'no'} |"
        )
    lines.append("")
    md_path = output_dir / "projection_constant_refit.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print()
    print(f"Fitted baseline_volatile_yield_fraction = {best_baseline:.6e}")
    print(f"Fitted reference_conversion_time_min    = {tau_star:.6g}")
    print(f"Objective at fit                        = {final_obj:.4f} dex")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

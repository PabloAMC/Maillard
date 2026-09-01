#!/usr/bin/env python3
"""Re-derive the per-lane hydrolysate sulfur `upstream_observability_factor`s.

RETRACTED 2026-08-27 (Wave I) — THIS SCRIPT REFUSES TO RUN
----------------------------------------------------------
Both of its fit targets are fabricated benchmarks and are now quarantined.  The cited
paper `10.1007/s10068-022-01194-w` reacts protein hydrolysates with GLUCOSE and FRUCTOSE
at pH 7.5 for 90 min and reports only RELATIVE GC-MS PEAK AREAS; it never mentions
2-furfurylthiol or 2-methyl-3-furanthiol and reports no absolute concentration for any
analyte.  The `conc_ppb` values this script fitted against therefore have no possible
source.  The one constant it moved (Methional `base_factor` 0.0045 -> 0.05623) has been
REVERTED in `src/recommend.py`, and the record it wrote
(`results/validation/hydrolysate_observability_rederivation.{json,md}`) is marked
RETRACTED.  The body below is preserved as the reproducible forensic record of what the
Wave H fit did — not as a live procedure.

Context (2026-08-27, Wave H) — historical
------------------------------------------
`src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` holds one
`base_factor` per compound — the fraction of the modelled volatile that is taken to
survive to the headspace in a hydrolysed-vegetable-protein (HVP) matrix.  The
forensics record (tasks/audit_remediation.md, "Re-anchor the WHOLE sulfur branch")
established what they are: they were fit to the two xylose HVP benchmarks
(PMC9905368) WITH THE BARRIERS FROZEN.  They constrain the PRODUCT of the lane, not
any barrier, so when the chemistry underneath them changes they must be re-derived,
not compensated for by bending a barrier.

Wave G1 changed the chemistry underneath them (the fabricated radical chain was
removed and MFT moved onto the real norfuraneol route), and Wave H refit
`thiol_addition_norfuraneol` against Hofmann1998.  This script re-derives the
factors against the new engine exactly the way they were originally derived:
solve for the factor that best reproduces the measured HVP concentrations.

The one thing that has NOT changed is the physics of the quantity: an observability
factor is a surviving FRACTION.  It lives in (0, 1].  `_resolve_upstream_observability_factor`
enforces that with `max(1.0e-4, min(1.0, ...))`.  So the re-derivation is a
CONSTRAINED one, and a solution that wants to sit outside the constraint is not a
fit — it is a statement that this layer cannot explain the residual.

Fit targets
-----------
The two literature xylose HVP benchmarks ONLY:

    spi_hvp_xylose_120C_PMC9905368
    wheat_gluten_hvp_xylose_120C_PMC9905368

Synthetic lanes (Internal2026 / ProtocolPilot2026) and the `external_validation/`
hold-out are forbidden and are excluded by assertion.  That has a consequence worth
stating: `bis(2-methyl-3-furyl) disulfide` appears in NO literature benchmark of this
lane — its only comparators are the synthetic snapshots — so its base_factor is NOT
re-derivable here and is reported as such rather than quietly fitted to a synthetic
target.

Decision rules — identical to scripts/generators/refit_sulfur_barriers_hofmann.py
---------------------------------------------------------------------------------
1. ADMISSIBILITY.  If the unconstrained optimum lies outside (0, 1], the parameter
   is reported SATURATED and the incumbent is kept.  A boundary value is not a fit,
   and pinning it at 1.0 would additionally collapse the modelled soy-vs-wheat
   observability ranking against the same clamp.
2. MATERIALITY.  A factor moves only for a gain of at least MIN_OBJECTIVE_GAIN_DEX.
3. CONSERVATIVE EDGE.  Among admissible values within INDIFFERENCE_DEX of the
   profile minimum, adopt the SMALLEST factor — the least lane-favourable point the
   data cannot distinguish from the optimum.

Usage
-----
    python scripts/generators/rederive_hydrolysate_observability.py
    python scripts/generators/rederive_hydrolysate_observability.py --apply
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src import recommend as recommend_module
from src.benchmark_validation import evaluate_benchmark
from src.recommend import _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES, _normalize_chemical_name
from src.uncertainty_propagation import _benchmark_signal_origin

TARGETS = (
    data_paths.benchmark_path("spi_hvp_xylose_120C_PMC9905368"),
    data_paths.benchmark_path("wheat_gluten_hvp_xylose_120C_PMC9905368"),
)
RECOMMEND_SOURCE = ROOT / "src" / "recommend.py"

# The compounds carrying a base_factor, in the spelling used in the source table.
PROFILE_COMPOUNDS = (
    "Methional",
    "2-Furfurylthiol",
    "2-Methyl-3-furanthiol",
    "bis(2-methyl-3-furyl) disulfide",
)

FACTOR_FLOOR = 1.0e-4     # the clamp `_resolve_upstream_observability_factor` applies
FACTOR_CEILING = 1.0      # a surviving fraction cannot exceed unity
GRID_POINTS = 241         # log-spaced over [FACTOR_FLOOR, FACTOR_CEILING]
MIN_OBJECTIVE_GAIN_DEX = 0.02
INDIFFERENCE_DEX = 0.01


def _guarded_targets() -> List[Path]:
    paths: List[Path] = []
    for path in TARGETS:
        assert path.exists(), f"fit target missing: {path}"
        assert "external_validation" not in path.parts, (
            f"hold-out benchmark reached the fit-target selector: {path}"
        )
        assert "Internal2026" not in path.name and "ProtocolPilot2026" not in path.name, (
            f"synthetic benchmark reached the fit-target selector: {path}"
        )
        origin = _benchmark_signal_origin(path)
        assert origin == "external_literature", (
            f"fit target must be literature-sourced, got signal origin {origin!r}"
        )
        paths.append(path)
    return paths


def _row_profile_key(comparison) -> str:
    """Which observability profile a benchmark row belongs to.

    Benchmark files label these rows "2-Methyl-3-furanthiol (MFT)" etc., while the
    runtime resolves the factor from the SPECIES name ("2-methyl-3-furanthiol").
    Match on the matched species name first, then on the benchmark label with any
    trailing parenthetical abbreviation stripped.
    """
    candidates = [comparison.matched_name, comparison.compound]
    for candidate in candidates:
        if not candidate:
            continue
        normalized = _normalize_chemical_name(candidate)
        if normalized in _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES:
            return normalized
        stripped = _normalize_chemical_name(re.sub(r"\s*\([^)]*\)\s*$", "", str(candidate)))
        if stripped in _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES:
            return stripped
    return _normalize_chemical_name(comparison.compound)


def _rows(paths: List[Path]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for path in paths:
        evaluation = evaluate_benchmark(path)
        if not evaluation.supported:
            continue
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            if comparison.measured_ppb <= 0.0 or comparison.predicted_ppb <= 0.0:
                continue
            rows.append(
                {
                    "benchmark_id": evaluation.benchmark_id,
                    "compound": comparison.compound,
                    "normalized": _row_profile_key(comparison),
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "ratio": comparison.predicted_ppb / comparison.measured_ppb,
                    "abs_log10_error": abs(
                        math.log10(comparison.predicted_ppb / comparison.measured_ppb)
                    ),
                }
            )
    return rows


def _objective_for(compound: str, paths: List[Path]) -> Tuple[float, List[Dict[str, Any]]]:
    normalized = _normalize_chemical_name(compound)
    rows = [row for row in _rows(paths) if row["normalized"] == normalized]
    if not rows:
        return float("nan"), rows
    return sum(row["abs_log10_error"] for row in rows) / len(rows), rows


def _evaluate_with(compound: str, factor: float, paths: List[Path]):
    key = _normalize_chemical_name(compound)
    profile = _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES[key]
    saved = float(profile["base_factor"])
    try:
        profile["base_factor"] = float(factor)
        return _objective_for(compound, paths)
    finally:
        profile["base_factor"] = saved


def _log_grid() -> List[float]:
    low, high = math.log10(FACTOR_FLOOR), math.log10(FACTOR_CEILING)
    step = (high - low) / (GRID_POINTS - 1)
    return [10.0 ** (low + index * step) for index in range(GRID_POINTS)]


def _rewrite_base_factor(text: str, compound: str, value: float) -> str:
    pattern = re.compile(
        r'(_normalize_chemical_name\("%s"\)\s*:\s*\{"base_factor":\s*)([0-9.eE+-]+)'
        % re.escape(compound)
    )
    match = pattern.search(text)
    if match is None:
        raise SystemExit(f"could not locate base_factor for {compound!r}")
    return text[: match.start()] + f"{match.group(1)}{value:.4g}" + text[match.end():]


# RETRACTED 2026-08-27 (Wave I). Both fit targets are fabricated benchmarks and are now
# quarantined: the cited paper 10.1007/s10068-022-01194-w uses glucose/fructose at pH 7.5
# for 90 min and reports only RELATIVE PEAK AREAS, never mentions FFT or MFT, and cannot be
# the source of any absolute ppb value in those files. This script therefore has no valid
# target and refuses to run. It is kept (not deleted) as the reproducible record of what the
# Wave H fit did, and so that the retraction is discoverable from the code as well as from
# results/validation/hydrolysate_observability_rederivation.{json,md}.
# To revive it, a REAL literature HVP benchmark must exist and be named here.
_RETRACTED = True
_RETRACTION_MESSAGE = (
    "rederive_hydrolysate_observability.py is RETRACTED (2026-08-27, Wave I).\n"
    "Its only fit targets -- spi_hvp_xylose_120C_PMC9905368 and "
    "wheat_gluten_hvp_xylose_120C_PMC9905368 -- are fabricated and quarantined:\n"
    "  the cited paper 10.1007/s10068-022-01194-w reacts hydrolysates with glucose/fructose\n"
    "  at pH 7.5 for 90 min and reports only RELATIVE GC-MS peak areas. It never mentions\n"
    "  2-furfurylthiol or 2-methyl-3-furanthiol and reports no absolute concentrations.\n"
    "The Methional base_factor it fitted (0.0045 -> 0.05623) has been REVERTED in\n"
    "src/recommend.py. Do not re-run this script against those files.\n"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    if _RETRACTED:
        sys.stderr.write(_RETRACTION_MESSAGE)
        return 2

    paths = _guarded_targets()
    grid = _log_grid()

    record: Dict[str, Any] = {
        "generated_by": "scripts/generators/rederive_hydrolysate_observability.py",
        "fit_targets": [path.name for path in paths],
        "objective": (
            "mean |log10(predicted_ppb / measured_ppb)| over the matched rows of the two "
            "literature xylose HVP benchmarks, per compound"
        ),
        "constraint": (
            f"base_factor is a surviving fraction and lives in "
            f"[{FACTOR_FLOOR}, {FACTOR_CEILING}] — the same clamp "
            "`_resolve_upstream_observability_factor` applies at runtime"
        ),
        "forbidden_as_fit_targets": [
            "data/benchmarks/external_validation/** (hold-out; excluded by assertion)",
            "*Internal2026*  (synthetic; excluded by assertion)",
            "*ProtocolPilot2026*  (synthetic; excluded by assertion)",
        ],
        "decision_rules": {
            "min_objective_gain_dex": MIN_OBJECTIVE_GAIN_DEX,
            "indifference_dex": INDIFFERENCE_DEX,
            "conservative_edge": (
                "among admissible values within indifference_dex of the profile minimum, "
                "adopt the SMALLEST factor"
            ),
            "admissibility": (
                "an unconstrained optimum outside (0, 1] is reported SATURATED and the "
                "incumbent is kept — a boundary value is not a fit"
            ),
        },
        # What this re-derivation ACTUALLY changed when it was first run, 2026-08-27. Kept
        # here because the script is idempotent: re-running it after the value has been
        # applied reports "no further move" (the right convergence check, and what the
        # committed artifact shows), which would otherwise erase the record of the move.
        "applied_history": [
            {
                "date": "2026-08-27",
                "compound": "Methional",
                "from": 0.0045,
                "to": 0.05623,
                "unconstrained_optimum": 0.06391,
                "objective_before_dex": 1.1560,
                "objective_after_dex": 0.0619,
                "residuals_before": "16.2x under (SPI) / 12.4x under (wheat gluten)",
                "residuals_after": "1.30x under (SPI) / 1.01x (wheat gluten)",
                "basis": (
                    "conservative edge of the 0.05623-0.07356 indifference band; interior "
                    "and identified, unlike the two thiols whose optima sit above the "
                    "physical maximum of 1.0"
                ),
            }
        ],
        "compounds": {},
    }

    adopted: Dict[str, float] = {}
    for compound in PROFILE_COMPOUNDS:
        key = _normalize_chemical_name(compound)
        incumbent = float(_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES[key]["base_factor"])
        baseline_objective, baseline_rows = _objective_for(compound, paths)
        entry: Dict[str, Any] = {
            "incumbent_base_factor": incumbent,
            "literature_rows": baseline_rows,
            "incumbent_objective": baseline_objective,
        }

        if not baseline_rows:
            entry["decision"] = (
                "NOT DERIVABLE — this compound has no row in any LITERATURE benchmark of "
                "this lane; its only comparators are the synthetic Internal2026 / "
                "ProtocolPilot2026 snapshots, which are forbidden as fit targets. The "
                "incumbent is kept as an unconstrained legacy estimate."
            )
            record["compounds"][compound] = entry
            print(f"\n{compound}: {entry['decision'].splitlines()[0]}")
            continue

        # Unconstrained optimum, read off the residuals: the factor scales the
        # prediction linearly until the clamp bites, so the per-row optimum is
        # incumbent / ratio and the joint optimum is their geometric mean.
        unconstrained = incumbent * 10.0 ** (
            sum(-math.log10(row["ratio"]) for row in baseline_rows) / len(baseline_rows)
        )
        entry["unconstrained_optimum"] = unconstrained
        entry["unconstrained_optimum_admissible"] = (
            FACTOR_FLOOR <= unconstrained <= FACTOR_CEILING
        )

        profile = []
        for factor in grid:
            objective, _rows_at = _evaluate_with(compound, factor, paths)
            profile.append({"base_factor": factor, "objective": objective})
        best = min(point["objective"] for point in profile)
        entry["profile_min_objective"] = best
        entry["achievable_gain_dex"] = baseline_objective - best
        entry["profile"] = profile

        if not entry["unconstrained_optimum_admissible"]:
            entry["decision"] = (
                f"SATURATED — the unconstrained optimum is {unconstrained:.3g}, "
                f"{unconstrained / FACTOR_CEILING:.1f}x ABOVE the physical maximum of "
                f"{FACTOR_CEILING:.1f} for a surviving fraction. This layer cannot explain "
                f"the residual (it can only suppress, and the model already under-predicts "
                f"this lane by {10 ** baseline_objective:.1f}x). Incumbent kept."
            )
        elif entry["achievable_gain_dex"] < MIN_OBJECTIVE_GAIN_DEX:
            entry["decision"] = (
                f"IMMATERIAL — best achievable gain {entry['achievable_gain_dex']:.4f} dex "
                f"< {MIN_OBJECTIVE_GAIN_DEX} dex; incumbent kept"
            )
        else:
            band = [
                point["base_factor"]
                for point in profile
                if point["objective"] <= best + INDIFFERENCE_DEX
            ]
            chosen = min(band)
            adopted[compound] = chosen
            entry["indifference_band"] = [min(band), max(band)]
            entry["adopted_base_factor"] = chosen
            entry["decision"] = (
                f"RE-DERIVED {incumbent:.4g} -> {chosen:.4g} (unconstrained optimum "
                f"{unconstrained:.4g}, admissible; conservative edge of the "
                f"{min(band):.4g}-{max(band):.4g} indifference band; gain "
                f"{entry['achievable_gain_dex']:.4f} dex)"
            )
        record["compounds"][compound] = entry
        print(f"\n{compound}: {entry['decision']}")

    record["adopted"] = adopted

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "hydrolysate_observability_rederivation.json"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    lines = [
        "# Hydrolysate sulfur observability factors — re-derivation (Wave H)",
        "",
        "`src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` holds one "
        "`base_factor` per compound: the fraction of the modelled volatile taken to survive "
        "to the headspace in an HVP matrix. These constrain the PRODUCT of the lane, not a "
        "barrier, so after Wave G1's chemistry change and Wave H's Hofmann-only barrier "
        "refit they are re-derived here — the same way they were originally derived, "
        "against the same two literature benchmarks — rather than being compensated for by "
        "bending a barrier.",
        "",
        f"Fit targets: {', '.join(path.name for path in paths)} (literature only; the "
        "synthetic lanes and the hold-out are excluded by assertion).",
        "",
        f"Constraint: {record['constraint']}.",
        "",
        "## What this re-derivation changed (2026-08-27)",
        "",
        "| Compound | From | To | Objective | Residuals |",
        "| --- | --- | --- | --- | --- |",
    ] + [
        f"| {item['compound']} | {item['from']:.4g} | {item['to']:.4g} | "
        f"{item['objective_before_dex']:.4f} -> {item['objective_after_dex']:.4f} dex | "
        f"{item['residuals_before']} -> {item['residuals_after']} |"
        for item in record["applied_history"]
    ] + [
        "",
        "This script is idempotent: re-run after the value is applied and it reports no "
        "further move, which is the convergence check. The table below is that re-run, so "
        "its residuals are the CURRENT ones.",
        "",
        "| Compound | Incumbent | Unconstrained optimum | Decision |",
        "| --- | --- | --- | --- |",
    ]
    for compound, entry in record["compounds"].items():
        optimum = entry.get("unconstrained_optimum")
        lines.append(
            f"| {compound} | {entry['incumbent_base_factor']:.4g} | "
            f"{('%.4g' % optimum) if optimum is not None else 'n/a'} | "
            f"{entry['decision']} |"
        )
    lines += [
        "",
        "## Per-row residuals at the incumbent",
        "",
        "| Benchmark | Compound | Measured ppb | Predicted ppb | Fold error |",
        "| --- | --- | --- | --- | --- |",
    ]
    for entry in record["compounds"].values():
        for row in entry.get("literature_rows", []):
            fold = 1.0 / row["ratio"] if row["ratio"] else float("nan")
            direction = "under" if fold >= 1.0 else "over"
            lines.append(
                f"| {row['benchmark_id']} | {row['compound']} | {row['measured_ppb']:.4g} | "
                f"{row['predicted_ppb']:.4g} | {max(fold, 1.0 / fold):.1f}x {direction} |"
            )
    lines.append("")
    md_path = output_dir / "hydrolysate_observability_rederivation.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {json_path}")
    print(f"Wrote {md_path}")

    if args.apply and adopted:
        text = RECOMMEND_SOURCE.read_text(encoding="utf-8")
        for compound, value in adopted.items():
            text = _rewrite_base_factor(text, compound, value)
        RECOMMEND_SOURCE.write_text(text, encoding="utf-8")
        print("Applied: " + ", ".join(
            f"{compound} -> {value:.4g}" for compound, value in adopted.items()
        ))
    elif args.apply:
        print("No factor moved; src/recommend.py left untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

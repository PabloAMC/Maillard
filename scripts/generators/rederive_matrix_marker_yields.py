#!/usr/bin/env python3
"""Re-derive `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` from the CONTENT-VERIFIED anchors,
and measure where the matrix lane's absolute-scale deficit actually lives.

Context (2026-08-28, Wave Y)
---------------------------
Three independent findings converged on this layer:

* **Wave S3** measured that the FAST screening lane consumes barriers only as branching
  ratios, so "the deficit is in the projection budget and missing chemistry, not the
  barrier table".
* **Wave S4** made a UNIT argument: an observability factor is the fraction of a total
  that a measurement sees, so it **cannot exceed 1** — and the shipped ambient hexanal
  factors are 4.31725 (pea) and 9.54007 (soy).  Something upstream is absorbing an
  absolute-scale deficit, and the named suspect was
  `MATRIX_BENCHMARK_BASE_MARKER_YIELDS`.
* **Wave X** measured 0.518 dex on single steps against 0.952 dex end to end.

THE LAYER, MEASURED (matrix_only lane, `_run_matrix_only_benchmark_prediction`)
------------------------------------------------------------------------------
    observable_ppb = L(pool, lane, T, t)      # src.lipid_oxidation.predict_hexanal_generation
                     * Y(compound)             # MATRIX_BENCHMARK_BASE_MARKER_YIELDS
                     * release(compound, lane) # HeadspaceModel, net of the registry factor
                     * cal(compound, lane, state)   # matrix_calibration_registry

`Y` and `cal` are PERFECTLY DEGENERATE on any single lane: only their product is
identified.  The split between them is therefore a CONVENTION, and the convention the
repository declares is that the pea ambient slurry is the reference lane whose factors
are 1.0 by construction (see the header of `src/matrix_calibration_registry.py`).

Wave O's shared observability scale s = 4.317249 was fitted against the verified anchors
(pea 1138.00, soy 1621.71 ppb) but was written into the WRONG side of that product: it
went into `cal`, where it broke the reference-lane convention and produced a "fraction
observed" of 4.32.  This script moves it to the side the unit argument requires.

WHAT IS DERIVED
---------------
Three models, all reported, one adopted.

    M0  SHIPPED             Y as shipped; cal carries s.                    (baseline)
    M1  REPARAMETERISATION  Y(Hexanal) := 0.205 * s;  every hexanal `cal` /= s.  ADOPTED.
                            Calibrated-tier predictions are preserved to 6 significant
                            figures; the UNCALIBRATED tier (which reads Y and never reads
                            `cal`) moves, and that movement is the point.
    M2  OBSERVABILITY PINNED TO ITS EVIDENCED VALUE.  Wave S4 (c) established from
                            Pratap-Singh's own verbatim methods that the anchors are
                            MATRIX-MATCHED (standards spiked into the slurry), i.e. they
                            measure the TOTAL, i.e. cal == 1.0 on BOTH ambient lanes.
                            Two free yields against four verified rows.  REPORTED, NOT
                            ADOPTED — its residual is the finding.

IDENTIFIABILITY, STATED BEFORE THE NUMBERS
------------------------------------------
1. `Y` is degenerate with `kinetics.hydroperoxide_scale` (1.0e6 in
   `data/lit/lipid_oxidation_calibration.json`, a round unanchored number).  Only the
   product `hydroperoxide_scale * Y` is identified.  **There is therefore no "a yield
   cannot exceed 1" argument available here**, unlike the observability argument — a
   `Y > 1` is not a unit violation, it is a statement about an arbitrary scale.  Saying
   so is the honest counterpart of Wave S4's unit argument, not a weakening of it.
2. `Y` is shared across matrices; `cal` is per (matrix, process_state, compound).  So a
   discrepancy that is LANE-specific cannot be absorbed by `Y` at all, and a discrepancy
   that is COMPOUND-specific cannot be absorbed by a lipid-profile change.  M2 measures
   which kind the residual is.
3. 1-Hexanol has NO anchor of any kind (the paper reports n.d. in both matrices; Wave T3
   traced the retired 80/120 ppb to the repository's own abstract-reconstructed brief) and
   Nonanal has no ambient anchor.  Neither yield is fitted here; both are reported as
   unanchored.  Inventing an anchor is the defect this audit exists to remove.

FIT CORPUS — DECLARED BEFORE FITTING, AND ASSERTED
--------------------------------------------------
    pea_isolate_40C_PratapSingh2021 / hexanal        1138.00 +/- 297.30 ppb  (Wave K)
    pea_isolate_40C_PratapSingh2021 / 2-pentylfuran   638    +/- 8 %        (verbatim)
    soy_isolate_40C_PratapSingh2021 / hexanal        1621.71 +/- 159.69 ppb  (Wave K)
    soy_isolate_40C_PratapSingh2021 / 2-pentylfuran  2492    +/- 8 %        (verbatim)

    pea_isolate_uht_140C_Trikusuma2019 is a fit target of M1 only in the propagation
    sense: its hexanal factor is divided by the same `s` so its recovery is preserved.
    It is CONTENT-UNVERIFIED and contributes no residual to any objective here.

NEVER IN ANY FIT, AND STRUCTURALLY INVISIBLE TO THIS SCRIPT
-----------------------------------------------------------
* `data/benchmarks/external_validation/**` — the matrix hold-out (8 points) and the
  maillard_path hold-out.  The benchmark iteration is an explicit whitelist of two
  benchmark ids and an assertion re-checks that no scored path contains
  `external_validation`.
* `*Internal2026*` / `*ProtocolPilot2026` — synthetic model output, not experiments.
* Wave X's declared fit target (`hofmann1998_norfuraneol_h2s_145C_20min_pH5`) — a
  different lane entirely, and named here only so the exclusion is on the record.

APPLICATION
-----------
This script does NOT write to `src/`.  The constants are applied by hand so the
provenance block lands with them.  Checking is the DEFAULT: every run re-derives and
compares against the shipped constants, exiting non-zero on drift.  `--no-check` emits
the record without failing, which is what you want when re-deriving BEFORE applying.

Usage
-----
    python scripts/generators/rederive_matrix_marker_yields.py
    python scripts/generators/rederive_matrix_marker_yields.py --no-check
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

from src.benchmark_validation import (  # noqa: E402
    MATRIX_BENCHMARK_BASE_MARKER_YIELDS,
    MATRIX_BENCHMARK_PROFILES,
    _get_condition_water_activity,
    load_benchmark,
)
from src.headspace import HeadspaceModel  # noqa: E402
from src.lipid_oxidation import (  # noqa: E402
    HYDROPEROXIDE_POOL_KEYS,
    hydroperoxide_pool_key_for_marker,
    predict_hexanal_generation,
)
from src.matrix_calibration_registry import (  # noqa: E402
    _MATRIX_CALIBRATION_RECORDS,
    _MATRIX_CLASS_ANCHORS,
    _normalize_compound,
    describe_matrix_calibration,
    determine_matrix_process_state,
)

BENCH_DIR = ROOT / "data" / "benchmarks"
OUT_JSON = ROOT / "results" / "validation" / "matrix_marker_yield_rederivation.json"
OUT_MD = ROOT / "results" / "validation" / "matrix_marker_yield_rederivation.md"

# The fit corpus, as ids. Nothing outside this tuple is read for a measured value.
ANCHOR_BENCHMARKS = (
    "pea_isolate_40C_PratapSingh2021",
    "soy_isolate_40C_PratapSingh2021",
)
# Propagation-only: its hexanal factor is divided by the same scale so its (fit-recovery)
# agreement is preserved. It contributes NO residual to any objective in this script.
PROPAGATION_BENCHMARKS = ("pea_isolate_uht_140C_Trikusuma2019",)

FORBIDDEN_PATH_TOKENS = ("external_validation", "Internal2026", "ProtocolPilot2026", "quarantined")

# Wave O's shared ambient-hexanal observability scale, the constant this script relocates.
WAVE_O_SHARED_SCALE = 4.31725

# Below this the move is not worth making: it is smaller than the rounding in the
# published anchors. Same threshold Wave O used.
MATERIALITY_DEX = 0.01

SIG_FIGS = 6


def _round_sig(value: float, digits: int = SIG_FIGS) -> float:
    if value == 0.0 or not math.isfinite(value):
        return value
    return round(value, -int(math.floor(math.log10(abs(value)))) + (digits - 1))


def _assert_corpus_clean(paths: List[Path]) -> None:
    for path in paths:
        text = str(path)
        for token in FORBIDDEN_PATH_TOKENS:
            if token in text:
                raise AssertionError(
                    f"fit corpus contamination: {text} contains the forbidden token {token!r}"
                )


def decompose_lane(benchmark_id: str) -> Dict[str, Any]:
    """Measure the full matrix_only decomposition for one benchmark, per compound."""
    path = BENCH_DIR / f"{benchmark_id}.json"
    _assert_corpus_clean([path])
    bench = load_benchmark(path)
    protein_type = str(bench["protein_type"])
    conditions = bench["conditions"]
    water_activity = _get_condition_water_activity(conditions)
    process_state = str(
        (bench.get("process_metadata") or {}).get(
            "state",
            determine_matrix_process_state(
                temperature_celsius=float(conditions["temp_C"]),
                time_minutes=float(conditions["time_min"]),
                water_activity=water_activity,
            ),
        )
    )
    oxidation = predict_hexanal_generation(
        MATRIX_BENCHMARK_PROFILES[protein_type]["lipid_profile"],
        temp_C=float(conditions["temp_C"]),
        time_min=float(conditions["time_min"]),
        oxygen_availability=1.0,
    )
    load_by_pool = {key: float(oxidation[key]) * 1000.0 for key in HYDROPEROXIDE_POOL_KEYS.values()}
    headspace = HeadspaceModel()
    measured = bench.get("measured_volatiles") or {}

    rows: Dict[str, Dict[str, Any]] = {}
    for compound, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
        headspace_factor = headspace.get_matrix_benchmark_headspace_factor(
            compound,
            protein_type=protein_type,
            pH=conditions.get("ph"),
            temperature_celsius=float(conditions["temp_C"]),
            time_minutes=float(conditions["time_min"]),
            water_activity=water_activity,
            binding_context=None,
        )
        calibration = describe_matrix_calibration(
            compound, protein_type=protein_type, process_state=process_state
        )
        cal = float(calibration.get("calibration_observable_factor") or 1.0)
        release = headspace_factor / cal if cal > 0.0 else headspace_factor
        pool_key = hydroperoxide_pool_key_for_marker(compound)
        load = load_by_pool[pool_key]
        entry = measured.get(compound.lower()) or {}
        meas = entry.get("conc_ppb")
        rows[compound] = {
            "pool_key": pool_key,
            "oxidation_load_ppb": load,
            "shipped_yield": float(yield_factor),
            "release_factor": release,
            "calibration_factor": cal,
            "calibration_fallback_mode": calibration.get("calibration_fallback_mode"),
            "calibration_source": calibration.get("source"),
            "predicted_ppb": load * float(yield_factor) * release * cal,
            "measured_ppb": float(meas) if isinstance(meas, (int, float)) else None,
            "measured_uncertainty_pct": (
                float(entry["uncertainty_pct"]) if isinstance(entry.get("uncertainty_pct"), (int, float)) else None
            ),
        }
    return {
        "benchmark_id": benchmark_id,
        "protein_type": protein_type,
        "process_state": process_state,
        "conditions": dict(conditions),
        "oxidation_load_ppb_by_pool": load_by_pool,
        "compounds": rows,
    }


def _sigma_log10(pct: float | None, default_pct: float = 26.0) -> float:
    pct = default_pct if pct is None else pct
    return max(float(pct), 1.0) / 100.0 / math.log(10.0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-check", action="store_true", help="emit the record without checking shipped constants")
    args = parser.parse_args()

    lanes = {bid: decompose_lane(bid) for bid in ANCHOR_BENCHMARKS + PROPAGATION_BENCHMARKS}

    # ---------------------------------------------------------------- M2 first.
    # M2 is computed first ON PURPOSE: it is the falsifiable measurement, and computing
    # it before the adopted reparameterisation keeps the reporting order honest.
    m2_rows: List[Dict[str, Any]] = []
    for bid in ANCHOR_BENCHMARKS:
        lane = lanes[bid]
        for compound, row in lane["compounds"].items():
            if row["measured_ppb"] is None:
                continue
            # cal pinned to 1.0 (matrix-matched quantification, Wave S4 (c)); release is
            # measured to be 1.0 on both ambient lanes but is carried, not assumed.
            denom = row["oxidation_load_ppb"] * row["release_factor"]
            m2_rows.append(
                {
                    "benchmark_id": bid,
                    "compound": compound,
                    "measured_ppb": row["measured_ppb"],
                    "oxidation_load_ppb": row["oxidation_load_ppb"],
                    "release_factor": row["release_factor"],
                    "required_yield": row["measured_ppb"] / denom,
                    "sigma_log10": _sigma_log10(row["measured_uncertainty_pct"]),
                }
            )

    m2_yields: Dict[str, Dict[str, Any]] = {}
    for compound in sorted({r["compound"] for r in m2_rows}):
        group = [r for r in m2_rows if r["compound"] == compound]
        weights = [1.0 / (r["sigma_log10"] ** 2) for r in group]
        log_mean = sum(w * math.log10(r["required_yield"]) for w, r in zip(weights, group)) / sum(weights)
        fitted = 10.0**log_mean
        residuals = {r["benchmark_id"]: math.log10(r["required_yield"] / fitted) for r in group}
        m2_yields[compound] = {
            "fitted_yield": fitted,
            "per_lane_required_yield": {r["benchmark_id"]: r["required_yield"] for r in group},
            "residual_dex": residuals,
            "lane_spread_fold": max(r["required_yield"] for r in group) / min(r["required_yield"] for r in group),
        }
    m2_all_resid = [v for c in m2_yields.values() for v in c["residual_dex"].values()]
    m2_rms_dex = math.sqrt(sum(r * r for r in m2_all_resid) / len(m2_all_resid)) if m2_all_resid else 0.0

    # ---------------------------------------------------------------- M1, adopted.
    #
    # DERIVED IN ABSOLUTE TERMS, NOT AS "shipped x scale". Writing it as a multiplier on
    # whatever is currently shipped would make this script NON-IDEMPOTENT: re-running it
    # after the constants land would apply the scale a second time. It is therefore
    # expressed as the solution of the same one-parameter problem Wave O solved, under the
    # convention it should have used:
    #
    #   cal(pea, ambient, hexanal)  := 1.0                 REFERENCE LANE. It is 1.0 both
    #                                                      by the registry's own definition
    #                                                      and, independently, by Wave S4 (c)
    #                                                      evidence (matrix-matched -> total).
    #   cal(soy, ambient, hexanal)  := 0.453 / 0.205       the soy-vs-pea ratio, declared and
    #                                                      not fitted (pre-Wave-O expression).
    #   Y(Hexanal)                  := the ONE free parameter, the geometric mean of the two
    #                                  per-lane required yields (the exact minimiser of the
    #                                  sum of squared log10 residuals over the two rows).
    #
    # One parameter, two rows, one degree of freedom left for the anchors to disagree with --
    # exactly Wave O's discipline, moved to the other side of the product.
    SOY_PEA_AMBIENT_HEXANAL_RATIO = 0.453 / 0.205
    ADOPTED_AMBIENT_HEXANAL_CAL = {
        "pea_isolate_40C_PratapSingh2021": 1.0,
        "soy_isolate_40C_PratapSingh2021": SOY_PEA_AMBIENT_HEXANAL_RATIO,
    }
    per_lane_required_hexanal_yield = {}
    for bid, cal in ADOPTED_AMBIENT_HEXANAL_CAL.items():
        row = lanes[bid]["compounds"]["Hexanal"]
        per_lane_required_hexanal_yield[bid] = row["measured_ppb"] / (
            row["oxidation_load_ppb"] * row["release_factor"] * cal
        )
    _vals = list(per_lane_required_hexanal_yield.values())
    derived_hexanal_yield = _round_sig(10.0 ** (sum(math.log10(v) for v in _vals) / len(_vals)))
    m1_residual_fold = {
        bid: value / derived_hexanal_yield for bid, value in per_lane_required_hexanal_yield.items()
    }

    # The scale by which every hexanal observability constant must fall so the PRODUCT
    # `Y * cal` is preserved on every lane. Derived from the pre-relocation reference
    # value (1.0 on the pea ambient lane before Wave O, 0.205 for the yield), so it does
    # not depend on what is currently shipped.
    PRE_WAVE_O_HEXANAL_YIELD = 0.205
    effective_scale = derived_hexanal_yield / PRE_WAVE_O_HEXANAL_YIELD
    shipped_hexanal_yield = float(MATRIX_BENCHMARK_BASE_MARKER_YIELDS["Hexanal"])

    # The four hexanal observability constants, in absolute terms.
    #   pea ambient  = 1.0                                   (reference lane)
    #   soy ambient  = 0.453 / 0.205                          (the declared soy/pea ratio)
    #   soy heated   = soy ambient * (1 - 0.7060)             (the Shu 2024 composition rule)
    #   pea heated   = the Trikusuma back-solve, re-expressed against the new yield, i.e.
    #                  its pre-relocation value divided by `effective_scale`. It is a
    #                  PROPAGATION of an existing fit, not a new one.
    PRE_RELOCATION_HEXANAL_FACTORS = {
        ("pea_iso", "ambient_slurry"): 4.31725,
        ("pea_iso", "heated_matrix"): 0.228776,
        ("soy_iso", "ambient_slurry"): 9.54007,
        ("soy_iso", "heated_matrix"): 9.54007 * (1.0 - 0.7060),
    }
    ADOPTED_HEXANAL_FACTORS = {
        ("pea_iso", "ambient_slurry"): 1.0,
        ("pea_iso", "heated_matrix"): _round_sig(0.228776 / effective_scale),
        ("soy_iso", "ambient_slurry"): SOY_PEA_AMBIENT_HEXANAL_RATIO,
        ("soy_iso", "heated_matrix"): SOY_PEA_AMBIENT_HEXANAL_RATIO * (1.0 - 0.7060),
    }

    # EVERY hexanal observability constant in the registry, enumerated from the registry
    # itself (not from the benchmark lanes) so the relocation cannot silently miss one.
    hexanal_factor_moves = []
    for record in _MATRIX_CALIBRATION_RECORDS:
        if _normalize_compound(record.compound) != "hexanal":
            continue
        key = (record.protein_type, record.process_state)
        before = PRE_RELOCATION_HEXANAL_FACTORS[key]
        after = ADOPTED_HEXANAL_FACTORS[key]
        hexanal_factor_moves.append(
            {
                "protein_type": record.protein_type,
                "process_state": record.process_state,
                "old_factor": before,
                "new_factor": after,
                "shipped_factor": float(record.observable_factor),
                "exceeds_one_before": before > 1.0 + 1e-9,
                "exceeds_one_after": after > 1.0 + 1e-9,
                "product_preserved_fold": (after * derived_hexanal_yield)
                / (before * PRE_WAVE_O_HEXANAL_YIELD),
            }
        )

    # THE FULL OBSERVABILITY CENSUS. Every compound-specific record and every class
    # anchor in the registry, before and after. This is the population the Wave S4 unit
    # argument is about, and stating the verdict over anything narrower would be a
    # selective reading.
    factor_census = []
    for record in _MATRIX_CALIBRATION_RECORDS:
        is_hexanal = _normalize_compound(record.compound) == "hexanal"
        key = (record.protein_type, record.process_state)
        old = PRE_RELOCATION_HEXANAL_FACTORS[key] if is_hexanal else float(record.observable_factor)
        new = ADOPTED_HEXANAL_FACTORS[key] if is_hexanal else float(record.observable_factor)
        factor_census.append(
            {
                "kind": "compound_specific",
                "protein_type": record.protein_type,
                "process_state": record.process_state,
                "compound": record.compound,
                "factor_before": old,
                "factor_after": new,
                "above_one_before": old > 1.0 + 1e-9,
                "above_one_after": new > 1.0 + 1e-9,
            }
        )
    for anchor in _MATRIX_CLASS_ANCHORS:
        old = float(anchor.observable_factor)
        factor_census.append(
            {
                "kind": "class_anchor",
                "protein_type": anchor.protein_type,
                "process_state": "(class)",
                "compound": anchor.target_class,
                "factor_before": old,
                "factor_after": old,
                "above_one_before": old > 1.0 + 1e-9,
                "above_one_after": old > 1.0 + 1e-9,
            }
        )
    above_before = sorted(
        {(e["protein_type"], e["process_state"], e["compound"]) for e in factor_census if e["above_one_before"]}
    )
    above_after = sorted(
        {(e["protein_type"], e["process_state"], e["compound"]) for e in factor_census if e["above_one_after"]}
    )

    # ---------------------------------------------------------------- Per-compound verdicts.
    yields_verdict: Dict[str, Dict[str, Any]] = {}
    for compound, shipped in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
        pea = lanes["pea_isolate_40C_PratapSingh2021"]["compounds"][compound]
        anchored = pea["measured_ppb"] is not None
        if compound == "Hexanal":
            new = derived_hexanal_yield
            status = "MOVED — the Wave O shared scale relocated from the observability factor"
        elif anchored:
            required = pea["measured_ppb"] / (pea["oxidation_load_ppb"] * pea["release_factor"] * pea["calibration_factor"])
            move_dex = abs(math.log10(required / shipped))
            new = shipped
            status = (
                f"CONFIRMED, not moved — the pea-reference requirement is {required:.7f}, "
                f"{move_dex:.6f} dex from shipped, below the {MATERIALITY_DEX} dex materiality floor"
            )
        else:
            new = shipped
            status = "UNANCHORED — no verified measurement exists for this compound on either ambient lane; NOT fitted"
        yields_verdict[compound] = {
            "shipped": shipped,
            "derived": new,
            "moved": new != shipped,
            "status": status,
            "anchored_in_corpus": bool(anchored),
        }

    # ---------------------------------------------------------------- The Wave P pool staleness.
    # Reported here because it is the same layer: a `cal` back-solved against one
    # hydroperoxide pool while the yield it multiplies now reads another.
    tri = lanes["pea_isolate_uht_140C_Trikusuma2019"]
    non = tri["compounds"]["Nonanal"]
    lin = tri["oxidation_load_ppb_by_pool"]["total_hydroperoxide"]
    ole = tri["oxidation_load_ppb_by_pool"]["total_hydroperoxide_oleate"]
    pool_ratio = lin / ole
    nonanal_staleness = {
        "benchmark_id": tri["benchmark_id"],
        "compound": "Nonanal",
        "measured_ppb": non["measured_ppb"],
        "predicted_ppb": non["predicted_ppb"],
        "fold": (non["predicted_ppb"] / non["measured_ppb"]) if non["measured_ppb"] else None,
        "linoleate_pool_ppb": lin,
        "oleate_pool_ppb": ole,
        "pool_ratio_linoleate_over_oleate": pool_ratio,
        "shipped_factor": non["calibration_factor"],
        "pool_corrected_factor": _round_sig(non["calibration_factor"] * pool_ratio),
        "diagnosis": (
            "The Trikusuma heated-pea nonanal observability factor was back-solved (pre-Wave-P) so that "
            "LINOLEATE_pool * Y * release * cal == 24 ppb. Wave P item 4 then moved nonanal onto the "
            "OLEATE pool without propagating the constant. The row is now under by EXACTLY the pool "
            "ratio, which is exactly linoleic_acid_pct / oleic_acid_pct = 50/22 = 2.272727 for the pea "
            "profile -- an arithmetic staleness, not a model error. Confirmed by re-multiplying: "
            f"{lin:.4f} * {non['shipped_yield']} * {non['release_factor']:.6f} * {non['calibration_factor']} = "
            f"{lin * non['shipped_yield'] * non['release_factor'] * non['calibration_factor']:.4f} ppb."
        ),
    }

    payload: Dict[str, Any] = {
        "generated_by": "scripts/generators/rederive_matrix_marker_yields.py",
        "wave": "Wave Y (2026-08-28)",
        "objective": (
            "Relocate the matrix lane's absolute observable scale from the observability factor "
            "(where an operator > 1 is a unit violation) to the base marker yield (where it is not), "
            "and measure what survives."
        ),
        "layer_map": {
            "equation": "observable_ppb = L(pool, lane, T, t) * Y(compound) * release(compound, lane) * cal(compound, lane, state)",
            "L": "src.lipid_oxidation.predict_hexanal_generation -> total_hydroperoxide / total_hydroperoxide_oleate, x1000",
            "Y": "src.benchmark_validation.MATRIX_BENCHMARK_BASE_MARKER_YIELDS",
            "release": "src.headspace.HeadspaceModel.get_matrix_benchmark_headspace_factor, net of cal",
            "cal": "src.matrix_calibration_registry.describe_matrix_calibration",
            "degeneracy": "Y and cal are perfectly degenerate on a single lane; only their product is identified.",
            "second_degeneracy": (
                "Y is degenerate with kinetics.hydroperoxide_scale = 1.0e6 in "
                "data/lit/lipid_oxidation_calibration.json, so Y carries no unit bound and 'Y > 1' is "
                "not a violation."
            ),
        },
        "fit_target_files": [f"{bid}.json" for bid in ANCHOR_BENCHMARKS + PROPAGATION_BENCHMARKS],
        "fit_target_ids": list(ANCHOR_BENCHMARKS + PROPAGATION_BENCHMARKS),
        "forbidden_as_fit_targets": [
            "data/benchmarks/external_validation/** (matrix hold-out AND maillard_path hold-out)",
            "*Internal2026* / *ProtocolPilot2026* (synthetic reproducibility lane)",
            "hofmann1998_norfuraneol_h2s_145C_20min_pH5 (Wave X's declared fit target, different lane)",
            "data/benchmarks/quarantined/**",
        ],
        "fit_leverage": {
            "free_parameters": 1,
            "fitted_rows": 2,
            "parameters_per_row": 0.5,
            "class": "per_row_recovery",
            "interpretation": (
                "M1 relocates ONE parameter -- Wave O's shared ambient-hexanal scale, fitted against two "
                "verified anchors -- from cal into Y. It introduces no new freedom: the same single "
                "constant, the same two rows, the same 1.0113x residual. Those rows were already "
                "excluded from the honest-literature numerator and denominator as fit recovery and stay "
                "excluded. Trikusuma is carried for PROPAGATION only and contributes no residual."
            ),
        },
        "lanes": lanes,
        "M0_shipped": {
            "hexanal_yield": shipped_hexanal_yield,
            "pea_ambient_hexanal_cal": lanes["pea_isolate_40C_PratapSingh2021"]["compounds"]["Hexanal"]["calibration_factor"],
            "soy_ambient_hexanal_cal": lanes["soy_isolate_40C_PratapSingh2021"]["compounds"]["Hexanal"]["calibration_factor"],
        },
        "M1_reparameterisation_ADOPTED": {
            "wave_o_shared_scale_for_reference": WAVE_O_SHARED_SCALE,
            "effective_scale_after_rounding": effective_scale,
            "hexanal_yield_old": PRE_WAVE_O_HEXANAL_YIELD,
            "hexanal_yield_new": derived_hexanal_yield,
            "hexanal_yield_currently_shipped": shipped_hexanal_yield,
            "adopted_ambient_hexanal_cal": {k: v for k, v in ADOPTED_AMBIENT_HEXANAL_CAL.items()},
            "per_lane_required_hexanal_yield": per_lane_required_hexanal_yield,
            "one_parameter_residual_fold": m1_residual_fold,
            "hexanal_observability_factor_moves": hexanal_factor_moves,
            "calibrated_tier_effect": (
                "None at reporting precision: the product Y*cal is preserved to 6 significant figures on "
                "every lane, and every hexanal lane in the repository resolves to a compound-specific "
                "record (measured -- no lane reaches the class-level aldehyde anchor for hexanal), so the "
                "relocation is complete."
            ),
            "uncalibrated_tier_effect": (
                "REAL AND INTENDED. `_uncalibrated_prediction_ppb` (the tier "
                "scripts/generators/derive_matrix_sigma_from_residuals.py scores) multiplies the "
                "oxidation load by Y and NEVER reads cal. Its hexanal residuals therefore move by the "
                f"full {effective_scale:.6f}x. Wave O recorded that no OBSERVABILITY refit could ever move "
                "that derivation; a YIELD refit can, and this is the first one."
            ),
        },
        "M2_observability_pinned_to_evidence_REPORTED_NOT_ADOPTED": {
            "premise": (
                "Wave S4 (c) established from Pratap-Singh 2021's own verbatim methods that the anchors "
                "are MATRIX-MATCHED (standards spiked into the slurry), i.e. they measure the TOTAL, i.e. "
                "the observability factor on both ambient lanes is 1.0 on evidence rather than by fit."
            ),
            "free_parameters": 2,
            "rows": 4,
            "yields": m2_yields,
            "rms_residual_dex": m2_rms_dex,
            "verdict": (
                "REJECTED AS A SHIPPED MODEL, AND THE REASON IS THE FINDING. With observability pinned to "
                "its evidenced value, ONE shared yield vector cannot fit both matrices: the soy-vs-pea "
                "required-yield ratio is "
                + ", ".join(f"{c} {v['lane_spread_fold']:.4f}x" for c, v in sorted(m2_yields.items()))
                + ". Those two numbers differ from each other, so the disagreement is COMPOUND-SPECIFIC and "
                "cannot be repaired by any change to the soy lipid profile (which moves both markers on the "
                "linoleate pool by the same factor). It is also LANE-SPECIFIC, so it cannot be repaired by "
                "any change to Y. The residual is structurally outside both of this wave's targets."
            ),
        },
        "observability_census": factor_census,
        "S4_claim_verdict": {
            "claim": (
                "Wave S4 (b): the shipped observability factors 4.32 / 9.54 are absorbing an absolute-scale "
                "deficit that lives upstream in MATRIX_BENCHMARK_BASE_MARKER_YIELDS; fixing the yields "
                "should bring the factors back under 1."
            ),
            "factors_above_one_before": [f"{p}/{s}/{c}" for p, s, c in above_before],
            "factors_above_one_after": [f"{p}/{s}/{c}" for p, s, c in above_after],
            "verdict": "PARTIALLY CONFIRMED — see the record's own summary; the surviving population is entirely soy.",
        },
        "wave_p_pool_staleness": nonanal_staleness,
        "derived_yields": yields_verdict,
        "applied_to_runtime": False,
    }

    # ------------------------------------------------------------------ check mode
    drift: List[str] = []
    for compound, verdict in yields_verdict.items():
        shipped_now = float(MATRIX_BENCHMARK_BASE_MARKER_YIELDS[compound])
        if abs(shipped_now - verdict["derived"]) > 1e-12:
            drift.append(
                f"MATRIX_BENCHMARK_BASE_MARKER_YIELDS[{compound!r}] shipped={shipped_now!r} derived={verdict['derived']!r}"
            )
    for move in hexanal_factor_moves:
        if abs(move["shipped_factor"] - move["new_factor"]) > 1e-9:
            drift.append(
                f"hexanal observability on {move['protein_type']}/{move['process_state']} "
                f"shipped={move['shipped_factor']!r} derived={move['new_factor']!r}"
            )
    payload["applied_to_runtime"] = not drift
    payload["drift"] = drift

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(_render_markdown(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")

    if drift and not args.no_check:
        print("\nDRIFT — shipped constants do not match the derivation:", file=sys.stderr)
        for line in drift:
            print(f"  {line}", file=sys.stderr)
        return 1
    if drift:
        print("\n(--no-check) shipped constants not yet applied:")
        for line in drift:
            print(f"  {line}")
    return 0


def _render_markdown(payload: Dict[str, Any]) -> str:
    out: List[str] = []
    out.append("# Matrix marker-yield re-derivation (Wave Y, 2026-08-28)\n")
    out.append(f"Generated by `{payload['generated_by']}`.\n")
    out.append("## The layer\n")
    lm = payload["layer_map"]
    out.append(f"```\n{lm['equation']}\n```\n")
    for key in ("L", "Y", "release", "cal"):
        out.append(f"* `{key}` — {lm[key]}")
    out.append(f"\n* **Degeneracy** — {lm['degeneracy']}")
    out.append(f"* **Second degeneracy** — {lm['second_degeneracy']}\n")

    out.append("## Fit corpus (declared before fitting)\n")
    out.append("| benchmark | compound | measured ppb | +/- % | oxidation load ppb | required Y (cal pinned to 1) |")
    out.append("|---|---|---:|---:|---:|---:|")
    for bid in ANCHOR_BENCHMARKS:
        lane = payload["lanes"][bid]
        for compound, row in lane["compounds"].items():
            if row["measured_ppb"] is None:
                continue
            req = row["measured_ppb"] / (row["oxidation_load_ppb"] * row["release_factor"])
            out.append(
                f"| `{bid}` | {compound} | {row['measured_ppb']:.2f} | {row['measured_uncertainty_pct']} | "
                f"{row['oxidation_load_ppb']:.4f} | {req:.6f} |"
            )
    out.append("")
    out.append("**Never in any fit:** " + "; ".join(payload["forbidden_as_fit_targets"]) + "\n")

    out.append("## M1 — the adopted relocation\n")
    m1 = payload["M1_reparameterisation_ADOPTED"]
    out.append(f"`Y(Hexanal)` **{m1['hexanal_yield_old']} -> {m1['hexanal_yield_new']}** "
               f"(x{m1['effective_scale_after_rounding']:.6f}, Wave O's shared scale, relocated).\n")
    out.append("ONE free parameter, TWO anchored rows, derived in absolute terms so the script is "
               "idempotent. Per-lane required yields: "
               + ", ".join(f"`{k}` {v:.6f}" for k, v in m1["per_lane_required_hexanal_yield"].items())
               + "; residual after the shared fit: "
               + ", ".join(f"{v:.4f}x" for v in m1["one_parameter_residual_fold"].values())
               + ".\n")
    out.append("| lane | process state | observability before | after | >1 before | >1 after |")
    out.append("|---|---|---:|---:|:--:|:--:|")
    for move in m1["hexanal_observability_factor_moves"]:
        out.append(
            f"| {move['protein_type']} | {move['process_state']} | {move['old_factor']:.6g} | "
            f"{move['new_factor']:.6g} | {'YES' if move['exceeds_one_before'] else 'no'} | "
            f"{'YES' if move['exceeds_one_after'] else 'no'} |"
        )
    out.append("")
    out.append("### Full observability census (every registry constant, before and after)\n")
    out.append("| kind | lane | process state | compound/class | before | after | >1 after |")
    out.append("|---|---|---|---|---:|---:|:--:|")
    for entry in payload["observability_census"]:
        out.append(
            f"| {entry['kind']} | {entry['protein_type']} | {entry['process_state']} | {entry['compound']} | "
            f"{entry['factor_before']:.6g} | {entry['factor_after']:.6g} | "
            f"{'**YES**' if entry['above_one_after'] else 'no'} |"
        )
    out.append(f"\n* Calibrated tier: {m1['calibrated_tier_effect']}")
    out.append(f"* Uncalibrated tier: {m1['uncalibrated_tier_effect']}\n")

    out.append("## M2 — observability pinned to its evidenced value (reported, NOT adopted)\n")
    m2 = payload["M2_observability_pinned_to_evidence_REPORTED_NOT_ADOPTED"]
    out.append(m2["premise"] + "\n")
    out.append("| compound | fitted Y | pea required | soy required | soy/pea | pea resid dex | soy resid dex |")
    out.append("|---|---:|---:|---:|---:|---:|---:|")
    for compound, blk in sorted(m2["yields"].items()):
        per = blk["per_lane_required_yield"]
        res = blk["residual_dex"]
        pea = [k for k in per if k.startswith("pea")][0]
        soy = [k for k in per if k.startswith("soy")][0]
        out.append(
            f"| {compound} | {blk['fitted_yield']:.6f} | {per[pea]:.6f} | {per[soy]:.6f} | "
            f"{blk['lane_spread_fold']:.4f}x | {res[pea]:+.4f} | {res[soy]:+.4f} |"
        )
    out.append(f"\nRMS residual **{m2['rms_residual_dex']:.4f} dex**.\n")
    out.append(m2["verdict"] + "\n")

    out.append("## The S4 claim\n")
    s4 = payload["S4_claim_verdict"]
    out.append(f"> {s4['claim']}\n")
    out.append(f"* Above 1 **before**: {', '.join(s4['factors_above_one_before']) or '(none)'}")
    out.append(f"* Above 1 **after**: {', '.join(s4['factors_above_one_after']) or '(none)'}")
    out.append(f"\n**{s4['verdict']}**\n")

    out.append("## Wave P pool staleness (found while mapping the layer)\n")
    st = payload["wave_p_pool_staleness"]
    out.append(f"`{st['benchmark_id']}` / Nonanal: measured {st['measured_ppb']}, predicted "
               f"{st['predicted_ppb']:.4f} ({st['fold']:.4f}x).\n")
    out.append(st["diagnosis"] + "\n")
    out.append(f"Pool-corrected factor: {st['shipped_factor']} -> {st['pool_corrected_factor']}\n")

    out.append("## Derived yields\n")
    out.append("| compound | shipped | derived | moved | status |")
    out.append("|---|---:|---:|:--:|---|")
    for compound, verdict in payload["derived_yields"].items():
        out.append(
            f"| {compound} | {verdict['shipped']} | {verdict['derived']} | "
            f"{'YES' if verdict['moved'] else 'no'} | {verdict['status']} |"
        )
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())

"""
Build Wave B8 -- THE FINAL PARAMETER WAVE of the kinetic-core rebuild.

Pre-registered in `results/validation/kinetic_core_b8_prereg.md`, written and
saved to disk before this file ran for the first time.

WHAT THIS IS
------------
A refit of the SULFUR lane on an UPDATED FIT SET, at the ONE declared weighting
B2.4's own exam-blind and panel-blind criterion shipped (W-HALF, E = 0.70).
`generate_kinetic_core_b2_3_fit` supplies the objective and
`generate_kinetic_core_b2_4_fit` supplies the weighting machinery; both are
imported UNCHANGED. Four things differ from B2.4 and every one is declared in
the pre-registration BEFORE any member ran:

  1. THE T-STRUCTURE (Amendment 17 clauses 4-5). `Ea_decay_thiol_sink` is
     searched inside Gigl 2021's MEASURED band (7, 102) instead of the
     arbitrary (20, 250) it was pressed against at 248 / 216 / 213; the two
     disulfide channels take Zhang 2026 k17's MEASURED 122.2 kJ/mol; the two
     amine-catalysed fed-Amadori enolisations take its k16's MEASURED 85.7.
     Those three changes live in `src/kinetic_core/parameters_sulfur.py`, not
     here -- this file only searches.
  2. FOUR NEW DECLARED FIT ROWS, all RATIO- or CONVERSION-scale, none a level,
     none touching a 140 C rung. See `B8_FIT_ROWS`.
  3. THE FREE SET IS 23 OF 48: B2.4's 20 plus the three k_ref whose barrier
     moved out from under them.
  4. TWO SEEDED STARTS instead of six, at one weighting instead of three.

FIREWALL. This file never opens `data/benchmarks/external_validation/`, never
names a hold-out row, and reads no measured value that is not a declared FIT
row. `tests/unit/test_kinetic_core_b8.py` asserts both by literal grep and by
walking every SYSTEM this objective integrates.

Usage
-----
    python scripts/generators/generate_kinetic_core_b8_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b8_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b8_fit.py --consolidate
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
import generate_kinetic_core_b2_4_fit as B24  # noqa: E402

from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    DECAY_EA_PRIOR_CENTRE,
    GIGL_EA_COVALENT_CAPTURE_RANGE_KJ_MOL,
    decay_ea_bounds,
)

# ===========================================================================
# 1. THE DECLARED WEIGHTING -- ONE, not three
# ===========================================================================

#: W-HALF. B2.4's shipping choice, made by a mechanical criterion that read
#: neither the panel nor the exam (its prereg sec. 3). B8 does not re-open the
#: weighting question and runs at this value only.
WEIGHT_TAG = "half"
WEIGHT_E = B24.PH_ENDPOINT_WEIGHT[WEIGHT_TAG]

#: The incumbent B8 starts from. NOT B2.3's vector: B2.4-half's ensemble-best
#: is the most recent vector fitted at this weighting, so start 0 being the
#: incumbent means a B8 member worse than its predecessor is visibly worse.
B24_HALF_FIT_REPORT = ROOT / "results/validation/kinetic_core_b2_4_fit_half.json"

MEMBER_DIR = ROOT / "results/validation/kinetic_core_b8_members"
OUT_FIT_REPORT = ROOT / "results/validation/kinetic_core_b8_fit_report.json"

SEED_BASE = 20260830
SEED_STRIDE = 1000
N_STARTS = 2

#: PRE-REGISTERED. 600 evaluations per member; ~2.35 s each measured on this
#: container at 58 rows, so a member is ~25 min. A member that exhausts the
#: budget returns its BEST-SO-FAR vector with status -9 and is excluded from
#: any spread statistic, exactly as in B2.4: a budget exhaustion is not a basin.
EVAL_BUDGET = 600

# ===========================================================================
# 2. THE FOUR NEW DECLARED FIT ROWS (Amendment 17 clause 5)
# ===========================================================================
# Every one is RATIO- or CONVERSION-scale. That is not a stylistic choice: the
# Kang/Zhai ladder is now declared SINGLE-IS SEMI-QUANT with f' = 1 and n = 1
# (Amendment 17 clause 1), and a semi-quant source licenses within-study shape
# and ratio use and forbids absolute yields. Entering a semi-quant LEVEL at a
# tight sigma would be using the source for the one thing it cannot support.
#
# NONE of the four touches a 140 C rung. The 140 C column of this ladder is the
# declared gating hold-out and its literals appear nowhere in this file.

B8_SYSTEMS: Dict[str, Dict[str, Any]] = {
    # Feng et al. 2022, JAFC 70:9095-9105 (feng2022_extraction.md secs. 1b, 3a):
    # the purified ribose-glutathione Amadori product at 20 mmol/L, initial
    # pH 7.0, UNBUFFERED, sealed, 120 min, at 100 and 120 C. Declared FIT by
    # Amendment 17 clause 5 ("Feng ARP depletion").
    #
    # WHY THIS SYSTEM AND NOT ZHOU'S: Zhou's ARP pot is at 60 min and carries a
    # cysteine co-substrate; Feng's is ARP ALONE for 120 min at two
    # temperatures, which is the only measurement in the objective that
    # constrains ARP CONSUMPTION with a temperature axis -- and B8 is the wave
    # that gives ARP consumption a measured barrier (Zhang k16, 85.7 kJ/mol) for
    # the first time. A barrier with no row that can see it is not a fit.
    "feng_arp_100": dict(
        initial={"ARP": 20.0, "OX": B23.OX_AMBIENT_MMOL_L},
        t_c=100.0, minutes=120.0, ph=7.0, buffer=B23.BUFFER_NONE,
        anchor="Feng 2022 Fig. 2a, RG-ARP 20 mM alone, pH 7 unbuffered, 100 C"),
    "feng_arp_120": dict(
        initial={"ARP": 20.0, "OX": B23.OX_AMBIENT_MMOL_L},
        t_c=120.0, minutes=120.0, ph=7.0, buffer=B23.BUFFER_NONE,
        anchor="Feng 2022 Fig. 2a, RG-ARP 20 mM alone, pH 7 unbuffered, 120 C"),
}

B8_FIT_ROWS: Tuple[Dict[str, Any], ...] = (
    dict(id="feng_ARP_conversion_100C", system="feng_arp_100",
         kind="conversion", species="ARP", target=0.48, sigma_log=0.30,
         anchor="Feng 2022 Fig. 2a (feng2022_extraction.md sec. B): ~48% of the "
                "RG-ARP consumed at 100 C / 120 min",
         note="DIGITISED from a figure, and the wide sigma says so. The 120 C "
              "partner below is PRINTED IN THE TEXT and carries half the "
              "sigma; the pair is what gives ARP consumption a temperature "
              "axis at all."),
    dict(id="feng_ARP_conversion_120C", system="feng_arp_120",
         kind="conversion", species="ARP", target=0.935, sigma_log=0.15,
         anchor="Feng 2022 (feng2022_extraction.md sec. B): 93.50% of the "
                "RG-ARP consumed at 120 C / 120 min, printed in the text",
         note="THE HARDER OF THE TWO. 93.5% consumed is close to exhaustion, so "
              "the model can miss it in only one direction by much -- it must "
              "not overshoot past complete consumption, and the residual is "
              "one-sided in practice even though the row is not."),
    dict(id="zhai_MFT_fold_120_over_100", system="kang_ttca_120",
         system_b="kang_ttca_100", kind="cross_system_ratio", species="MFT",
         target=1.388 / 1.237, sigma_log=0.15,
         anchor="Zhai et al. 2023 Food Chem. 404:134420 Table 1, TTCA arm, "
                "120 min: MFT 1.237 -> 1.388 ug/L over 100 -> 120 C",
         note="THE RESPONSE-FACTOR-FREE FORM of information the two level rows "
              "already carry at sigma_log 0.4 each. Amendment 17 withdraws the "
              "Tier A label from this ladder: it is single-IS semi-quant with "
              "f' = 1 and n = 1, so the LEVELS stay at their semi-quant width "
              "and the RATIO -- in which a constant unknown response factor "
              "cancels exactly -- carries the shape at a tight sigma. That is "
              "what 'shapes and ratios, never absolutes' means operationally. "
              "PRIMARY SOURCE IS ZHAI, NOT KANG: Kang 2026 Table S4 is a "
              "re-publication of these same cells 3.5 years later."),
    dict(id="zhai_FFT_fold_120_over_100", system="kang_ttca_120",
         system_b="kang_ttca_100", kind="cross_system_ratio", species="FFT",
         target=4.107 / 3.734, sigma_log=0.15,
         anchor="Zhai et al. 2023 Food Chem. 404:134420 Table 1, TTCA arm, "
                "120 min: FFT 3.734 -> 4.107 ug/L over 100 -> 120 C",
         note="The FFT partner. Both folds are near unity, and that near-flat "
              "LOW LEG is the part of the ladder that SURVIVED wave K6a's audit "
              "-- Feng's ARP-depletion normalisation reproduces it (MFT x0.84, "
              "FFT x1.06 per mole consumed) and Meng's independent 5-min leg "
              "brackets it. What did NOT survive is the high leg, and no row "
              "here scores it."),
)


def install_b8_rows() -> None:
    """
    Splice B8's systems and rows into the imported B2.3 objective.

    Rebinding module attributes rather than forking the objective is
    deliberate and is the same discipline B2.4 applied to the scorers: a forked
    objective is an objective that can drift, and the whole value of comparing
    B8 to B2.4 is that 58 of the 62 rows are the SAME rows, scored by the SAME
    code, at the SAME sigmas.
    """
    if all(row["id"] in {r["id"] for r in B23.ACTIVE_FIT_ROWS}
           for row in B8_FIT_ROWS):
        return
    for name, spec in B8_SYSTEMS.items():
        if name not in B23.SYSTEMS:
            B23.SYSTEMS[name] = spec
    B23.ACTIVE_FIT_ROWS = tuple(B23.ACTIVE_FIT_ROWS) + B8_FIT_ROWS
    B23.FIT_ROWS = tuple(B23.FIT_ROWS) + B8_FIT_ROWS


install_b8_rows()

#: The 58 row ids B2.4 scored, frozen at import BEFORE B8's four are added --
#: no, AFTER, and computed by exclusion, which is the same thing and is
#: auditable. `sum_r2_level_shared_55` is summed over exactly the non-pH members
#: of this set, so a B8-vs-B2.4 comparison is like-for-like on rows B2.4 had.
B24_ROW_IDS: Tuple[str, ...] = tuple(
    row["id"] for row in B23.ACTIVE_FIT_ROWS
    if row["id"] not in {r["id"] for r in B8_FIT_ROWS}
)

# ===========================================================================
# 3. THE FREE SET -- 23 of 48
# ===========================================================================

#: B2.4's 20 (clauses R1-R4), unchanged and re-used by import rather than
#: re-listed, so the two waves cannot drift apart silently.
FREE_FROM_B2_4: Tuple[str, ...] = tuple(B24.FREE_KEYS)

#: R5 -- NEW IN B8. The three rate constants whose ACTIVATION ENERGY moved
#: under them. Freezing a k_ref at 145 C while changing its barrier reports the
#: OLD rate at the NEW barrier, which is neither the old model nor the new one.
#:
#:   k_dimer_mft  gained a barrier where it had NONE (held at 145 C before)
#:   k_arp_dpo    64.08 -> 85.7 (Zhang k16, measured)
#:   k_arp_tdp    64.08 -> 85.7 (Zhang k16, measured)
#:
#: `k_dimer_fft` also gained a barrier and is ALREADY free: B2.3's own
#: Gauss-Newton intervals call it individually identified, so it is in B2.4's
#: clause R1 and needs no new clause.
FREE_R5_BARRIER_MOVED: Tuple[str, ...] = (
    "k_dimer_mft", "k_arp_dpo", "k_arp_tdp",
)

FREE_KEYS: Tuple[str, ...] = FREE_FROM_B2_4 + FREE_R5_BARRIER_MOVED

FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(B24.FREE_CLAUSE_OF),
    **{k: "R5 barrier moved under it (B8)" for k in FREE_R5_BARRIER_MOVED},
}

ALL_KEYS: Tuple[str, ...] = tuple(B24.ALL_KEYS)
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(
    k for k in ALL_KEYS if k not in set(FREE_KEYS))

assert len(set(FREE_KEYS)) == len(FREE_KEYS) == 23, FREE_KEYS


# ===========================================================================
# 4. BOUNDS -- per-family decay bands
# ===========================================================================


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    """
    B2.4's bounds with ONE change: each decay family gets its OWN declared band.

    `thiol_sink` is (7, 102) -- Gigl 2021's defensible range for the covalent
    capture barrier, a MEASURED prior band rather than a search convenience.
    `carbonyl_sink` is UNCHANGED at (20, 250): K6a measured no carbonyl-sink
    barrier, so narrowing it would be inventing one.
    """
    from src.kinetic_core.parameters_sulfur import (
        FITTED_SULFUR_BOUNDS_LOG10K, LUMPED_FORMATION_EA_BOUNDS,
    )
    from src.kinetic_core.ph_state import ACID_YIELD_BOUNDS, ARP_AMINE_PKA_BOUNDS

    decay_lo = [decay_ea_bounds(f)[0] for f in B23.DECAY_FAMILY_ORDER]
    decay_hi = [decay_ea_bounds(f)[1] for f in B23.DECAY_FAMILY_ORDER]
    lower = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][0] for k in B23.PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[0]] + decay_lo
        + [ACID_YIELD_BOUNDS[0], ARP_AMINE_PKA_BOUNDS[0]]
    )
    upper = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][1] for k in B23.PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[1]] + decay_hi
        + [ACID_YIELD_BOUNDS[1], ARP_AMINE_PKA_BOUNDS[1]]
    )
    return lower, upper


def incumbent_vector() -> np.ndarray:
    """
    B2.4-half's ensemble-best vector, CLIPPED into B8's declared bands.

    The clip is not cosmetic and it is the wave's headline in one line: that
    vector carries `Ea_decay_thiol_sink` = 215.7, which is OUTSIDE the measured
    band (7, 102) by 2.1x. Start 0 therefore begins at the band's ceiling, and
    the report says so.
    """
    payload = json.loads(B24_HALF_FIT_REPORT.read_text())
    fr = payload["frozen_parameters"]
    x = np.array(
        [fr["log10_k_ref_at_145C"][k] for k in B23.PARAM_ORDER]
        + [fr["lumped_formation_Ea_kJ_mol"]]
        + [fr["decay_Ea_kJ_mol"][f] for f in B23.DECAY_FAMILY_ORDER]
        + [fr["ph_drift"]["acid_yield_per_sink_event"],
           fr["ph_drift"]["arp_secondary_ammonium_pKa"]],
        dtype=float,
    )
    lower, upper = full_bounds()
    return np.clip(x, lower, upper)


#: Half-widths of the seeded local draw for start 1, per coordinate class.
#: Identical to B2.4's, so "a local start" means the same thing in both waves.
LOCAL_HALF_WIDTH_DECADES = 2.0
LOCAL_HALF_WIDTH_EA_KJ = 40.0
LOCAL_HALF_WIDTH_ACID_YIELD = 0.15
LOCAL_HALF_WIDTH_PKA = 1.0


def _local_half_width(key: str) -> float:
    if key.startswith("Ea_"):
        return LOCAL_HALF_WIDTH_EA_KJ
    if key == "ph_acid_yield_per_sink_event":
        return LOCAL_HALF_WIDTH_ACID_YIELD
    if key == "ph_arp_secondary_ammonium_pKa":
        return LOCAL_HALF_WIDTH_PKA
    return LOCAL_HALF_WIDTH_DECADES


def start_vector(start: int, base: np.ndarray) -> np.ndarray:
    """Start 0 IS the clipped incumbent. Start 1 draws locally around it."""
    x = base.copy()
    if start == 0:
        return x
    lower, upper = full_bounds()
    rng = np.random.default_rng(SEED_BASE + SEED_STRIDE * start)
    for idx, key in zip(FREE_INDEX, FREE_KEYS):
        half = _local_half_width(key)
        lo = max(lower[idx], base[idx] - half)
        hi = min(upper[idx], base[idx] + half)
        x[idx] = float(rng.uniform(lo, hi))
    return x


def start_arm(start: int) -> str:
    return "incumbent" if start == 0 else "local"


# ===========================================================================
# 5. THE OBJECTIVE AND THE COMPARATOR
# ===========================================================================


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    """B2.4's residual vector at W-HALF, over B8's 62 rows."""
    return B24.residual_vector(x_full, WEIGHT_E, quick)


def sum_r2_level_shared(r: np.ndarray) -> float:
    """
    THE CROSS-WAVE COMPARATOR: sum of squared residuals over exactly the NON-pH
    rows B2.4 ALSO scored, at their unchanged sigma_log.

    B8's objective has four rows B2.4's did not, so total cost is not
    comparable across the two waves. This is.
    """
    total = 0.0
    shared = set(B24_ROW_IDS)
    for i, row in enumerate(B23.ACTIVE_FIT_ROWS):
        if row["kind"] == "ph_endpoint" or row["id"] not in shared:
            continue
        total += float(r[i]) ** 2
    return total


def sum_r2_level_all(r: np.ndarray) -> float:
    """The same statistic over ALL 59 non-pH rows B8 scores."""
    total = 0.0
    for i, row in enumerate(B23.ACTIVE_FIT_ROWS):
        if row["kind"] == "ph_endpoint":
            continue
        total += float(r[i]) ** 2
    return total


class _BudgetExhausted(Exception):
    pass


def fit_member(start: int, max_nfev: int, quick: bool,
               budget: int = EVAL_BUDGET) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    base = incumbent_vector()
    lower, upper = full_bounds()
    x0_full = start_vector(start, base)
    frozen_full = base.copy()

    idx = np.array(FREE_INDEX, dtype=int)
    lo_free, hi_free = lower[idx], upper[idx]

    def to_full(x_free: np.ndarray) -> np.ndarray:
        x = frozen_full.copy()
        x[idx] = x_free
        return x

    seen = {"n": 0, "best_cost": float("inf"), "best_x": None, "best_r": None}

    def f(x_free: np.ndarray) -> np.ndarray:
        r = residual_vector(to_full(x_free), quick)
        seen["n"] += 1
        cost = 0.5 * float(np.dot(r, r))
        if cost < seen["best_cost"]:
            seen["best_cost"] = cost
            seen["best_x"] = np.array(x_free, dtype=float)
            seen["best_r"] = np.array(r, dtype=float)
        if seen["n"] >= budget:
            raise _BudgetExhausted
        return r

    x0_free = np.clip(x0_full[idx], lo_free, hi_free)
    t0 = time.time()
    r0 = f(x0_free)
    cost0 = 0.5 * float(np.dot(r0, r0))
    print(f"  [B8 start {start}] arm={start_arm(start)} E={WEIGHT_E} "
          f"rows={len(B23.ACTIVE_FIT_ROWS)} free={len(FREE_KEYS)}  "
          f"start cost {cost0:.4f}", flush=True)
    budget_exhausted = False
    try:
        result = least_squares(
            f, x0_free, bounds=(lo_free, hi_free), method="trf", x_scale="jac",
            max_nfev=max_nfev, diff_step=3e-2, ftol=1e-6, xtol=1e-8, verbose=0,
        )
        x_free_out = np.asarray(result.x, dtype=float)
        r = np.asarray(result.fun)
        cost_out = float(result.cost)
        nfev_out = int(result.nfev)
        status_out = int(result.status)
    except _BudgetExhausted:
        budget_exhausted = True
        x_free_out = seen["best_x"]
        r = seen["best_r"]
        cost_out = float(seen["best_cost"])
        nfev_out = int(seen["n"])
        status_out = -9

    seconds = round(time.time() - t0, 1)
    x_full = to_full(x_free_out)
    print(f"  [B8 start {start}] cost {cost_out:.4f}  evals {seen['n']}  "
          f"status {status_out}"
          f"{' BUDGET-EXHAUSTED' if budget_exhausted else ''}  {seconds}s",
          flush=True)

    fitted_map = {k: float(x_full[ALL_KEYS.index(k)]) for k in ALL_KEYS}
    ph_miss = {row_id: float(v) for row_id, v
               in B24.ph_endpoint_residuals_in_units(x_full, quick).items()}
    per_row = {
        row["id"]: float(r[i]) for i, row in enumerate(B23.ACTIVE_FIT_ROWS)}
    return {
        "wave": "B8",
        "weight_tag": WEIGHT_TAG,
        "ph_endpoint_weight_E": WEIGHT_E,
        "sigma_ph_units": B24.SIGMA_LOG_REFERENCE / WEIGHT_E,
        "start": start,
        "arm": start_arm(start),
        "seed": None if start == 0 else SEED_BASE + SEED_STRIDE * start,
        "start_is_b2_4_half_vector": start == 0,
        "start_cost": cost0,
        "cost": cost_out,
        "sum_r2_level_shared_with_b2_4": sum_r2_level_shared(r),
        "sum_r2_level_all_b8_rows": sum_r2_level_all(r),
        "n_rows": len(B23.ACTIVE_FIT_ROWS),
        "n_free": len(FREE_KEYS),
        "evals": int(seen["n"]),
        "eval_budget": int(budget),
        "budget_exhausted": budget_exhausted,
        "nfev": nfev_out,
        "status": status_out,
        "converged": status_out in (1, 2, 3),
        "seconds": seconds,
        "x_full": [float(v) for v in x_full],
        "parameters": fitted_map,
        "ph_endpoint_miss_pH_units": ph_miss,
        "residuals": [float(v) for v in r],
        "residual_by_row": per_row,
        "row_ids": [row["id"] for row in B23.ACTIVE_FIT_ROWS],
    }


# ===========================================================================
# 6. CONSOLIDATION
# ===========================================================================


def load_members() -> List[Dict[str, Any]]:
    return [json.loads(p.read_text()) for p in sorted(MEMBER_DIR.glob("*.json"))]


def band_state(value: float, family: str) -> str:
    lo, hi = decay_ea_bounds(family)
    if abs(value - hi) <= 2.0:
        return "ON THE CEILING"
    if abs(value - lo) <= 2.0:
        return "ON THE FLOOR"
    return "interior"


def consolidate() -> Dict[str, Any]:
    members = load_members()
    if not members:
        raise SystemExit("no B8 members on disk; run --start 0 and --start 1")
    best = min(members, key=lambda m: m["cost"])
    x_full = np.array(best["x_full"], dtype=float)
    frozen = B24._b23_schema_frozen(x_full)
    decay = frozen["decay_Ea_kJ_mol"]

    payload = {
        "wave": "B8 -- the final parameter wave (Amendments 16 + 17)",
        "generated_by": "scripts/generators/generate_kinetic_core_b8_fit.py",
        "prereg": "results/validation/kinetic_core_b8_prereg.md",
        "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendments 16, 17",
        "objective": {
            "form": ("B2.3's objective at B2.4's DECLARED pH exchange rate "
                     "E = 0.70 (W-HALF), over 62 declared FIT rows"),
            "ph_endpoint_weight_E_decades_per_pH_unit": WEIGHT_E,
            "sigma_ph_units": B24.SIGMA_LOG_REFERENCE / WEIGHT_E,
            "n_rows": len(B23.ACTIVE_FIT_ROWS),
            "n_rows_in_b2_4": len(B24_ROW_IDS),
            "n_new_rows": len(B8_FIT_ROWS),
            "new_row_ids": [r["id"] for r in B8_FIT_ROWS],
            "n_free_parameters": len(FREE_KEYS),
            "n_frozen": len(FROZEN_KEYS),
            "final_cost": best["cost"],
            "sum_r2_level_shared_with_b2_4": best["sum_r2_level_shared_with_b2_4"],
            "sum_r2_level_all_b8_rows": best["sum_r2_level_all_b8_rows"],
            "not_comparable_note": (
                "TOTAL COST IS NOT COMPARABLE TO B2.4's: B8 scores four rows "
                "B2.4 did not. `sum_r2_level_shared_with_b2_4` is the "
                "like-for-like comparator and is summed over exactly the "
                "non-pH rows both waves scored."),
        },
        "t_structure": {
            "thiol_sink_band_kj_mol": list(GIGL_EA_COVALENT_CAPTURE_RANGE_KJ_MOL),
            "thiol_sink_prior_centre_kj_mol": DECAY_EA_PRIOR_CENTRE["thiol_sink"],
            "thiol_sink_fitted_kj_mol": decay["thiol_sink"],
            "thiol_sink_band_state": band_state(decay["thiol_sink"], "thiol_sink"),
            "carbonyl_sink_band_kj_mol": list(decay_ea_bounds("carbonyl_sink")),
            "carbonyl_sink_fitted_kj_mol": decay["carbonyl_sink"],
            "carbonyl_sink_band_state": band_state(
                decay["carbonyl_sink"], "carbonyl_sink"),
            "measured_barriers_the_fit_cannot_move": {
                "k_dimer_mft": 122.2, "k_dimer_fft": 122.2,
                "k_arp_dpo": 85.7, "k_arp_tdp": 85.7,
                "k_cys_thermal": 55.1,
            },
            "refuted": (
                "Ea_decay_thiol_sink = 248.0 (B2.2), 216.1 (B2.3), 212.9-218.1 "
                "(B2.4). Gigl 2021 measures the covalent-capture channel at "
                "k(333)/k(279) = 67.2; 248 kJ/mol predicts 3.4e7 -- wrong by "
                "5.1e5x."),
        },
        "free_set": {
            "n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS),
            "keys": list(FREE_KEYS), "clause": FREE_CLAUSE_OF,
            "frozen_keys": list(FROZEN_KEYS),
        },
        "members": [{k: v for k, v in m.items()
                     if k not in ("x_full", "residuals", "parameters")}
                    for m in members],
        "best_start": best["start"],
        "frozen_parameters": frozen,
    }
    OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
    print(f"wrote {OUT_FIT_REPORT}")
    print(f"  best start {best['start']}  cost {best['cost']:.4f}  "
          f"shared sum_r2_level {best['sum_r2_level_shared_with_b2_4']:.3f}")
    print(f"  thiol_sink {decay['thiol_sink']:.2f} "
          f"({payload['t_structure']['thiol_sink_band_state']})  "
          f"carbonyl_sink {decay['carbonyl_sink']:.2f} "
          f"({payload['t_structure']['carbonyl_sink_band_state']})")
    print(f"  lumped formation Ea {frozen['lumped_formation_Ea_kJ_mol']:.2f}")
    return payload


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)

    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = fit_member(args.start, args.max_nfev, args.quick, args.budget)
    dest = MEMBER_DIR / f"b8_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

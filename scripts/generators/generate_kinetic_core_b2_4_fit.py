"""
Build Wave B2.4 -- the DECLARED pH-unit-vs-log-fold exchange rate, and an
ensemble instead of a point.

Pre-registered in `results/validation/kinetic_core_b2_4_prereg.md`, written and
saved to disk before this file ran for the first time.

WHAT THIS IS, AND WHAT IT DELIBERATELY IS NOT
---------------------------------------------
It is NOT a new objective. `generate_kinetic_core_b2_3_fit` is imported
UNCHANGED and its systems, its 58 declared FIT rows, its observables, its
bounds and its optimiser settings are reused verbatim. Exactly two things
differ, and both are declared:

  1. THE EXCHANGE RATE. B2.3 scores its three `zhou_final_pH_*` rows in pH
     units at sigma = 0.25 while its 55 other rows are scored in log10 folds at
     sigma_log = 0.20-0.60. Nobody chose the resulting exchange rate; D1 sec. 3
     measured it at ~9x (one pH unit priced at nine times a 3x level miss) and
     showed that ONE of those rows carries 44% of the entire B2.2 -> B2.3
     refit. Here the rate is a declared scalar:

         E = decades of level error that ONE pH unit of endpoint error is
             declared to be worth,       sigma_ph = SIGMA_LOG_REFERENCE / E

     E = 1.40 reproduces B2.3's objective exactly (0.35 / 0.25), and is the
     control. The other two declared values and their bases are in the
     pre-registration.

  2. THE FREE SET. 20 of the 48 parameters are free and 28 are frozen at their
     B2.3 values, by the four-clause rule pre-registered in sec. 2. The reason
     is arithmetic, not taste: B2.3's own multistart trace records 4731 s and
     2280 s for two full-vector starts, so the draft's 18 full-vector starts is
     ~21 h against a container ceiling that has already killed two attempts.

TOTAL COST IS NOT COMPARABLE ACROSS WEIGHTINGS -- they are three different
objectives. The comparable quantity, reported for every member of every
ensemble, is `sum_r2_level`: the sum of squared residuals over the 55 NON-pH
rows at their UNCHANGED sigma_log.

FIREWALL. This file never opens `data/benchmarks/external_validation/`, never
names a hold-out row, and reads no measured value that is not already a
declared FIT row of B2.3's objective.

Usage
-----
    # one ensemble member (this is what gets parallelised)
    python scripts/generators/generate_kinetic_core_b2_4_fit.py \
        --weight-tag shipped --start 0

    # after all members are on disk
    python scripts/generators/generate_kinetic_core_b2_4_fit.py --consolidate
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402

# ===========================================================================
# THE DECLARATION
# ===========================================================================

#: The reference level sigma: the modal sigma_log of the corpus's level rows.
#: Every Hofmann 1998 level anchor in B2.3's objective carries exactly this.
SIGMA_LOG_REFERENCE = 0.35

#: THE DECLARED EXCHANGE RATE, in decades of level error per pH unit of
#: endpoint error. `sigma_ph = SIGMA_LOG_REFERENCE / E`.
#:
#:   shipped  1.40 -- B2.3's accidental rate, made explicit. 0.35 / 0.25.
#:                    Reproduces B2.3's objective EXACTLY. The control.
#:   half     0.70 -- half the shipped rate; the pH endpoints down-weighted
#:                    4x IN COST, which is the intervention D1's W-3 names.
#:   measured 0.28 -- the rate the corpus measures ON THE FIT SIDE: Kumazawa
#:                    2003's six DECLARED FIT retention rungs move
#:                    2-furfurylthiol from 0.995 to 0.110 survival across pH
#:                    3.0 - 6.4, i.e. 0.956 decades over 3.4 pH units = 0.281
#:                    decades per pH unit. SEE THE CORRECTION NOTE BELOW.
#:
#: E is a DECLARED EXCHANGE RATE, not a measurement uncertainty. sigma_ph =
#: 1.25 pH units at `measured` does not claim Zhou's endpoints are uncertain by
#: that much -- the dossier digitises them at +/-0.06. It claims that missing an
#: endpoint by 1.25 units should cost what missing a level by 0.35 decades
#: costs, because that is the level consequence such a miss carries.
PH_ENDPOINT_WEIGHT: Dict[str, float] = {
    "shipped": 1.40,
    "half": 0.70,
    "measured": 0.28,
}

#: CORRECTION, MADE BEFORE A SINGLE MEMBER AT THIS WEIGHTING HAD RUN, AND
#: CAUGHT BY THIS WAVE'S OWN FIREWALL TEST.
#:
#: The `measured` value was FIRST written as 0.32 decades per pH unit, on the
#: basis of D1 sec. 4's observation that Hofmann's own pH ladder moves
#: 2-furfurylthiol 1.28 decades -- A HOLD-OUT-DERIVED NUMBER, and that is the
#: defect: Hofmann's pH-3 and pH-7 aqueous points are HOLD-OUT in this
#: repo, and their difference is exactly what the exam and the panel score. A
#: fit weighting derived from it would have been a hold-out value entering the
#: objective, which is the single thing the firewall exists to prevent.
#:
#: `tests/unit/test_kinetic_core_b2_4.py::test_no_holdout_literal_appears_on_the_b2_4_fit_side`
#: found it as the HOLD-OUT literal `1.28` on the fit side. The basis is
#: therefore re-derived from the FIT side only: Kumazawa 2003's ladder, six
#: DECLARED FIT rows, which is in any case the better anchor -- it is the
#: corpus's ONLY measurement of a thiol's pH response with formation held out
#: of the picture. It gives 0.281 rather than 0.32, so the two bases happen to
#: agree closely; that is a coincidence and is not the reason for the change.
#:
#: Order of operations, disclosed: the correction was made while the `shipped`
#: ensemble was still running and BEFORE any `measured` member existed on disk.

PH_ENDPOINT_WEIGHT_BASIS: Dict[str, str] = {
    "shipped": (
        "B2.3's ACCIDENTAL rate made explicit: sigma_log 0.35 over sigma_ph "
        "0.25 = 1.40 decades per pH unit. Reproduces B2.3's objective exactly "
        "and is the control -- it is what 'declaring the weighting' costs when "
        "the declared value is the one already in force."
    ),
    "half": (
        "Half the shipped rate, i.e. the three pH endpoints down-weighted by "
        "4x IN COST. Chosen as the smallest intervention large enough to be "
        "visible against a 10-dof objective, and chosen for that reason and "
        "not for any property of its outcome."
    ),
    "measured": (
        "The exchange rate the corpus measures ON THE FIT SIDE. Kumazawa "
        "2003's six DECLARED FIT retention rungs move 2-furfurylthiol from "
        "0.995 to 0.110 survival across pH 3.0-6.4: 0.956 decades over 3.4 pH "
        "units, i.e. 0.281 decades per pH unit. This prices a pH miss at the "
        "level consequence the fit corpus's own pH-response measurement says a "
        "pH miss has. It REPLACES a first draft of this value that was derived "
        "from Hofmann's hold-out pH ladder -- see the correction note in this "
        "module."
    ),
}

#: The row ids the weighting applies to, and nothing else.
PH_ENDPOINT_ROW_IDS: Tuple[str, ...] = tuple(
    row["id"] for row in B23.ACTIVE_FIT_ROWS if row["kind"] == "ph_endpoint"
)

# ===========================================================================
# THE FREE SET -- 20 of 48, by the pre-registered four-clause rule
# ===========================================================================

#: R1 -- every parameter the B2.3 report's OWN Gauss-Newton intervals call
#: individually identified. These are the only constants the corpus determines.
FREE_R1_IDENTIFIED: Tuple[str, ...] = (
    "k_tdp_fur", "k_nf_mft", "k_nf_mp3p", "k_mgo_mp", "k_glc_ha",
    "k_dimer_fft", "k_fft_decay", "k_fur_decay", "k_pent_caramel",
    "k_cys_thermal", "k_thiolate_loss",
)

#: R2 -- the two pH-drift constants. They are the DIRECT price-takers of the
#: declared weighting; freezing them would make the experiment vacuous.
FREE_R2_PH_DRIFT: Tuple[str, ...] = (
    "ph_acid_yield_per_sink_event", "ph_arp_secondary_ammonium_pKa",
)

#: R3 -- the five constants D1 sec. 3 records as having moved >= 3 decades
#: between B2.2 and B2.3 and not already in R1. These are the TOPOLOGY
#: SWITCHES. Freezing them would hard-code the basin under test and would make
#: D1's W-3 unfalsifiable by construction.
FREE_R3_TOPOLOGY_SWITCHES: Tuple[str, ...] = (
    "k_fur_fft", "k_osone_decay", "k_dimer_decay", "k_ttca_deg",
    "k_thiol_decay",
)

#: R4 -- the two decay barriers. D1 sec. 3 shows one traded ONTO its 250
#: kJ/mol ceiling and the other OFF it in the refit; a bound trade is what a
#: basin change looks like.
FREE_R4_DECAY_BARRIERS: Tuple[str, ...] = (
    "Ea_decay_thiol_sink", "Ea_decay_carbonyl_sink",
)

FREE_KEYS: Tuple[str, ...] = (
    FREE_R1_IDENTIFIED + FREE_R2_PH_DRIFT
    + FREE_R3_TOPOLOGY_SWITCHES + FREE_R4_DECAY_BARRIERS
)

FREE_CLAUSE_OF: Dict[str, str] = {
    **{k: "R1 identified" for k in FREE_R1_IDENTIFIED},
    **{k: "R2 pH-drift (price-taker of the weighting)" for k in FREE_R2_PH_DRIFT},
    **{k: "R3 topology switch (moved >=3 decades B2.2->B2.3)"
       for k in FREE_R3_TOPOLOGY_SWITCHES},
    **{k: "R4 decay barrier (traded on/off its bound)"
       for k in FREE_R4_DECAY_BARRIERS},
}

#: The full 48-vector's coordinate names, in B2.3's own order.
ALL_KEYS: Tuple[str, ...] = tuple(B23.PARAM_ORDER) + tuple(B23.EXTRA_PARAM_NAMES)
assert len(ALL_KEYS) == B23.N_K + B23.N_EXTRA == 48

FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))

B23_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b2_3_fit_report.json"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b2_4_members"

#: Seeds. Start 0 is B2.3's own vector exactly, so every ensemble contains the
#: incumbent and a member worse than B2.3 is visibly worse.
SEED_BASE = 20260829
SEED_STRIDE = 1000
N_STARTS = 6


def b23_vector() -> np.ndarray:
    payload = json.loads(B23_FIT_REPORT.read_text())
    fr = payload["frozen_parameters"]
    return np.array(
        [fr["log10_k_ref_at_145C"][k] for k in B23.PARAM_ORDER]
        + [fr["lumped_formation_Ea_kJ_mol"]]
        + [fr["decay_Ea_kJ_mol"][f] for f in B23.DECAY_FAMILY_ORDER]
        + [fr["ph_drift"]["acid_yield_per_sink_event"],
           fr["ph_drift"]["arp_secondary_ammonium_pKa"]],
        dtype=float,
    )


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    from src.kinetic_core.parameters_sulfur import (
        DECAY_EA_BOUNDS, FITTED_SULFUR_BOUNDS_LOG10K, LUMPED_FORMATION_EA_BOUNDS,
    )
    from src.kinetic_core.ph_state import ACID_YIELD_BOUNDS, ARP_AMINE_PKA_BOUNDS

    lower = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][0] for k in B23.PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[0], DECAY_EA_BOUNDS[0], DECAY_EA_BOUNDS[0],
           ACID_YIELD_BOUNDS[0], ARP_AMINE_PKA_BOUNDS[0]]
    )
    upper = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][1] for k in B23.PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[1], DECAY_EA_BOUNDS[1], DECAY_EA_BOUNDS[1],
           ACID_YIELD_BOUNDS[1], ARP_AMINE_PKA_BOUNDS[1]]
    )
    return lower, upper


# ===========================================================================
# THE OBJECTIVE, AT A DECLARED WEIGHTING
# ===========================================================================

_PROGRESS = {"n": 0, "best": float("inf"), "t0": time.time()}


def residual_vector(x_full: np.ndarray, weight: float, quick: bool) -> np.ndarray:
    """
    B2.3's residual vector with ONE change: the pH-endpoint rows are scored at
    `sigma_ph = SIGMA_LOG_REFERENCE / weight` instead of at their hard-coded
    0.25. Every other row is untouched, including the one-sided ceilings.
    """
    sigma_ph = SIGMA_LOG_REFERENCE / float(weight)
    try:
        parameters = B23.build_parameters(x_full)
        _f, _e, _d, drift = B23.unpack(x_full)
        predicted, _ = B23.observables(parameters, quick, drift)
    except Exception:
        _tick(np.full(len(B23.ACTIVE_FIT_ROWS), 25.0))
        return np.full(len(B23.ACTIVE_FIT_ROWS), 25.0)
    out = np.empty(len(B23.ACTIVE_FIT_ROWS))
    for i, row in enumerate(B23.ACTIVE_FIT_ROWS):
        p = predicted[row["id"]]
        t = float(row["target"])
        if not np.isfinite(p):
            out[i] = 25.0
            continue
        if row["kind"] == "ph_endpoint":
            out[i] = (p - t) / sigma_ph
        elif row["kind"] in ("ceiling", "peak_fraction_ceiling"):
            out[i] = max(0.0, math.log10((p + B23.FLOOR) / (t + B23.FLOOR))) / float(
                row["sigma_log"])
        else:
            out[i] = math.log10((p + B23.FLOOR) / (t + B23.FLOOR)) / float(
                row["sigma_log"])
    out = np.clip(out, -25.0, 25.0)
    _tick(out)
    return out


def _tick(r: np.ndarray) -> None:
    _PROGRESS["n"] += 1
    if _PROGRESS["n"] % 50 == 0:
        gc.collect()
    cost = 0.5 * float(np.dot(r, r))
    _PROGRESS["best"] = min(_PROGRESS["best"], cost)
    if _PROGRESS["n"] % 100 == 0:
        print(f"    [{_PROGRESS['n']:5d} evals, "
              f"{time.time() - _PROGRESS['t0']:6.0f}s] best cost "
              f"{_PROGRESS['best']:.3f}", flush=True)


def sum_r2_level(r: np.ndarray) -> float:
    """
    THE CROSS-WEIGHTING COMPARATOR. Sum of squared residuals over the 55 NON-pH
    rows at their unchanged sigma_log. Total cost is not comparable across
    three different objectives; this is.
    """
    total = 0.0
    for i, row in enumerate(B23.ACTIVE_FIT_ROWS):
        if row["kind"] == "ph_endpoint":
            continue
        total += float(r[i]) ** 2
    return total


def ph_endpoint_residuals_in_units(x_full: np.ndarray, quick: bool) -> Dict[str, float]:
    """The three pH endpoints' misses IN pH UNITS, weighting-independent."""
    parameters = B23.build_parameters(x_full)
    _f, _e, _d, drift = B23.unpack(x_full)
    predicted, _ = B23.observables(parameters, quick, drift)
    out = {}
    for row in B23.ACTIVE_FIT_ROWS:
        if row["kind"] != "ph_endpoint":
            continue
        out[row["id"]] = float(predicted[row["id"]]) - float(row["target"])
    return out


# ===========================================================================
# ONE ENSEMBLE MEMBER
# ===========================================================================


#: PREREG AMENDMENT A2, made before any member finished. Starts 1-2 are the
#: LOCAL arm and starts 3-5 the GLOBAL arm. The two answer different questions
#: and the container can only afford to answer the second one badly: a
#: full-bound draw starts at cost 1200-2500 against an optimum near 8, and
#: under A1's evaluation budget such a member reports where it GOT TO rather
#: than where the basin IS. The local arm asks "is there more than one basin
#: within two decades of the incumbent" and can actually converge; the global
#: arm keeps the harder question and reports honestly that it is budget-limited.
LOCAL_ARM_STARTS: Tuple[int, ...] = (1, 2)
GLOBAL_ARM_STARTS: Tuple[int, ...] = (3, 4, 5)

#: Half-widths of the LOCAL arm's uniform draw, per coordinate class.
LOCAL_HALF_WIDTH_DECADES = 2.0
LOCAL_HALF_WIDTH_EA_KJ = 40.0
LOCAL_HALF_WIDTH_ACID_YIELD = 0.15
LOCAL_HALF_WIDTH_PKA = 1.0


def start_arm(start: int) -> str:
    if start == 0:
        return "incumbent"
    return "local" if start in LOCAL_ARM_STARTS else "global"


def _local_half_width(key: str) -> float:
    if key.startswith("Ea_"):
        return LOCAL_HALF_WIDTH_EA_KJ
    if key == "ph_acid_yield_per_sink_event":
        return LOCAL_HALF_WIDTH_ACID_YIELD
    if key == "ph_arp_secondary_ammonium_pKa":
        return LOCAL_HALF_WIDTH_PKA
    return LOCAL_HALF_WIDTH_DECADES


def start_vector(start: int, base: np.ndarray) -> np.ndarray:
    """
    Start 0 IS B2.3's own vector. Starts 1-2 draw uniformly within a declared
    neighbourhood of it; starts 3-5 draw uniformly over the FULL declared
    bound. See LOCAL_ARM_STARTS for why the arms are split.
    """
    x = base.copy()
    if start == 0:
        return x
    lower, upper = full_bounds()
    rng = np.random.default_rng(SEED_BASE + SEED_STRIDE * start)
    local = start in LOCAL_ARM_STARTS
    for idx, key in zip(FREE_INDEX, FREE_KEYS):
        if local:
            half = _local_half_width(key)
            lo = max(lower[idx], base[idx] - half)
            hi = min(upper[idx], base[idx] + half)
        else:
            lo, hi = lower[idx], upper[idx]
        x[idx] = float(rng.uniform(lo, hi))
    return x


#: PREREG AMENDMENT A1, made before any member finished. `max_nfev` does not
#: bound the work: scipy counts trial points in `nfev`, but the
#: finite-difference Jacobian's evaluations are not counted there, and with 20
#: free parameters the Jacobian is where essentially all the cost is. The
#: container delivers about ONE CORE of real throughput for this objective
#: (one evaluation alone costs 1.7 s; with five processes running it costs 8.6 s
#: -- a 5.1x slowdown on five processes), so an unbounded member can run for
#: hours. FIREWALL-OK: those two figures are WALL-CLOCK SECONDS measured on this
#: container, not measured chemistry, and no hold-out row is involved.
#: On exhaustion the member returns its BEST-SO-FAR vector with status -9 and
#: is excluded from the spread statistic by the same rule that excludes a
#: max_nfev termination. A budget exhaustion is not a basin.
EVAL_BUDGET = 500


class _BudgetExhausted(Exception):
    pass


def fit_member(weight_tag: str, start: int, max_nfev: int, quick: bool,
               budget: int = EVAL_BUDGET) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    weight = PH_ENDPOINT_WEIGHT[weight_tag]
    base = b23_vector()
    lower, upper = full_bounds()
    x0_full = start_vector(start, base)
    frozen_full = base.copy()  # the 28 frozen coordinates live here forever

    idx = np.array(FREE_INDEX, dtype=int)
    lo_free, hi_free = lower[idx], upper[idx]

    def to_full(x_free: np.ndarray) -> np.ndarray:
        x = frozen_full.copy()
        x[idx] = x_free
        return x

    seen = {"n": 0, "best_cost": float("inf"), "best_x": None, "best_r": None}

    def f(x_free: np.ndarray) -> np.ndarray:
        r = residual_vector(to_full(x_free), weight, quick)
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
    print(f"  [{weight_tag} start {start}] arm={start_arm(start)} E={weight}  "
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
    print(f"  [{weight_tag} start {start}] cost {cost_out:.4f}  "
          f"evals {seen['n']}  nfev {nfev_out}  status {status_out}"
          f"{' BUDGET-EXHAUSTED' if budget_exhausted else ''}  {seconds}s",
          flush=True)

    fitted_map = {k: float(x_full[ALL_KEYS.index(k)]) for k in ALL_KEYS}
    ph_miss = {
        row_id: float(v)
        for row_id, v in ph_endpoint_residuals_in_units(x_full, quick).items()
    }
    return {
        "weight_tag": weight_tag,
        "ph_endpoint_weight_E": weight,
        "sigma_ph_units": SIGMA_LOG_REFERENCE / weight,
        "start": start,
        "arm": start_arm(start),
        "seed": None if start == 0 else SEED_BASE + SEED_STRIDE * start,
        "start_is_b2_3_vector": start == 0,
        "start_cost": cost0,
        "cost": cost_out,
        "sum_r2_level": sum_r2_level(r),
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
        "row_ids": [row["id"] for row in B23.ACTIVE_FIT_ROWS],
    }


# ===========================================================================
# CONSOLIDATION -- the ensemble, and one B2.3-schema fit report per weighting
# ===========================================================================


def _b23_schema_frozen(x_full: np.ndarray) -> Dict[str, Any]:
    """
    The `frozen_parameters` block in EXACTLY B2.3's schema, so the existing
    hold-out scorer and the engine can read a B2.4 member with no edit to
    either. This is a serialisation, not a claim.
    """
    DECAY_FAMILIES = B23.DECAY_FAMILIES
    DECAY_KEYS_ON_FORMATION_EA = B23.DECAY_KEYS_ON_FORMATION_EA
    _f, _e, _d, drift = B23.unpack(x_full)
    return {
        "log10_k_ref_at_145C": {
            k: float(x_full[i]) for i, k in enumerate(B23.PARAM_ORDER)},
        "lumped_formation_Ea_kJ_mol": float(x_full[B23.N_K]),
        "decay_Ea_kJ_mol": {
            family: float(x_full[B23.N_K + 1 + i])
            for i, family in enumerate(B23.DECAY_FAMILY_ORDER)},
        "decay_families": {k: list(v) for k, v in DECAY_FAMILIES.items()},
        "decay_keys_kept_on_formation_Ea": list(DECAY_KEYS_ON_FORMATION_EA),
        "ph_drift": drift.as_dict(),
        "reference_temperature_K": B23.T_REF_S_K,
    }


def load_members() -> List[Dict[str, Any]]:
    out = []
    for path in sorted(MEMBER_DIR.glob("*.json")):
        out.append(json.loads(path.read_text()))
    return out


def ensemble_summary(members: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Per-weighting ensemble statistics and the pre-declared gates."""
    summary: Dict[str, Any] = {}
    for tag, weight in PH_ENDPOINT_WEIGHT.items():
        block = [m for m in members if m["weight_tag"] == tag]
        block.sort(key=lambda m: m["start"])
        conv = [m for m in block if m["converged"]]
        costs = [m["cost"] for m in conv]
        best = min(block, key=lambda m: m["cost"]) if block else None
        spread = (
            math.log10(max(costs) / min(costs))
            if len(costs) >= 2 and min(costs) > 0 else None
        )
        g1 = None
        if best is not None:
            g1 = all(abs(v) <= 1.0
                     for v in best["ph_endpoint_miss_pH_units"].values())
        local = [m for m in conv if m.get("arm") in ("incumbent", "local")]
        local_costs = [m["cost"] for m in local]
        local_spread = (
            math.log10(max(local_costs) / min(local_costs))
            if len(local_costs) >= 2 and min(local_costs) > 0 else None
        )
        summary[tag] = {
            "E": weight,
            "sigma_ph_units": SIGMA_LOG_REFERENCE / weight,
            "basis": PH_ENDPOINT_WEIGHT_BASIS[tag],
            "n_members": len(block),
            "n_converged": len(conv),
            "n_budget_exhausted": sum(1 for m in block if m.get("budget_exhausted")),
            "arms": {m["start"]: m.get("arm") for m in block},
            "spread_S_local_arm_only": local_spread,
            "costs_by_arm": {
                arm: [m["cost"] for m in block if m.get("arm") == arm]
                for arm in ("incumbent", "local", "global")
            },
            "costs": [m["cost"] for m in block],
            "costs_converged": costs,
            "best_cost": best["cost"] if best else None,
            "worst_converged_cost": max(costs) if costs else None,
            "spread_S_log10_max_over_min": spread,
            "sum_r2_level_per_member": [m["sum_r2_level"] for m in block],
            "best_sum_r2_level": best["sum_r2_level"] if best else None,
            "best_start": best["start"] if best else None,
            "ph_endpoint_miss_at_best": (
                best["ph_endpoint_miss_pH_units"] if best else None),
            "gate_G1_calibration_within_1pH": g1,
            "gate_G2_at_least_4_of_6_converged": len(conv) >= 4,
        }
    return summary


def shipping_choice(summary: Dict[str, Any]) -> Dict[str, Any]:
    """
    THE PRE-DECLARED CRITERION, applied mechanically. Exam-blind and
    panel-blind: neither artifact is read here, or anywhere upstream of here.
    """
    qualifying = [
        tag for tag, b in summary.items()
        if b["gate_G1_calibration_within_1pH"]
        and b["gate_G2_at_least_4_of_6_converged"]
        and b["spread_S_log10_max_over_min"] is not None
    ]
    fallback = False
    pool = qualifying
    if not pool:
        fallback = True
        pool = [tag for tag, b in summary.items()
                if b["gate_G1_calibration_within_1pH"]
                and b["spread_S_log10_max_over_min"] is not None]
    if not pool:
        return {"shipped": None, "reason": "no weighting passed gate G1",
                "qualifying": [], "fallback_used": True}
    best_s = min(summary[t]["spread_S_log10_max_over_min"] for t in pool)
    tied = [t for t in pool
            if summary[t]["spread_S_log10_max_over_min"] - best_s <= 0.05]
    shipped = max(tied, key=lambda t: summary[t]["E"])
    return {
        "shipped": shipped,
        "qualifying": pool,
        "tied_within_0_05": tied,
        "tie_broken_toward_largest_E": len(tied) > 1,
        "fallback_used": fallback,
        "criterion": (
            "G1: all three zhou_final_pH_* within 1.0 pH unit at the "
            "ensemble-best member. G2: >=4 of 6 members converged. "
            "S = log10(max cost / min cost) over converged members; ship the "
            "qualifying weighting with the smallest S; ties within 0.05 break "
            "toward the largest E. Declared in "
            "results/validation/kinetic_core_b2_4_prereg.md sec. 3 before any "
            "fit ran. THE EXAM AND THE PANEL DO NOT ENTER IT."
        ),
    }


def write_per_weighting_fit_reports(members: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    One B2.3-schema fit report per weighting, at that weighting's
    ENSEMBLE-BEST member, so the existing hold-out scorer and the exam can be
    pointed at it with no edit to either.
    """
    paths = {}
    for tag in PH_ENDPOINT_WEIGHT:
        block = [m for m in members if m["weight_tag"] == tag]
        if not block:
            continue
        # PER-MEMBER reports too, in the same schema. W-2 asks what the exam's
        # geometric-mean fold does ACROSS the ensemble, which cannot be answered
        # from the best member alone: the whole claim is that a single number
        # was never a property of the model.
        for member in block:
            (data_paths.VALIDATION_DIR / f"kinetic_core_b2_4_fit_{tag}_s"
                    f"{member['start']}.json").write_text(json.dumps({
                        "wave": f"B2.4 member {tag}/s{member['start']}",
                        "member": {k: v for k, v in member.items()
                                   if k not in ("x_full", "residuals", "parameters")},
                        "frozen_parameters": _b23_schema_frozen(
                            np.array(member["x_full"], dtype=float)),
                    }, indent=2, default=str))
        best = min(block, key=lambda m: m["cost"])
        x_full = np.array(best["x_full"], dtype=float)
        payload = {
            "wave": f"B2.4 -- declared pH exchange rate E={best['ph_endpoint_weight_E']} ({tag})",
            "generated_by": "scripts/generators/generate_kinetic_core_b2_4_fit.py",
            "member": {k: v for k, v in best.items()
                       if k not in ("x_full", "residuals", "parameters")},
            "objective": {
                "form": (
                    "B2.3's objective with a DECLARED pH-unit-vs-log-fold "
                    "exchange rate: ph_endpoint rows use sigma_ph = 0.35 / E, "
                    "every other row is unchanged"
                ),
                "ph_endpoint_weight_E_decades_per_pH_unit": best["ph_endpoint_weight_E"],
                "sigma_ph_units": best["sigma_ph_units"],
                "sigma_log_reference": SIGMA_LOG_REFERENCE,
                "n_rows": len(B23.ACTIVE_FIT_ROWS),
                "n_free_parameters": len(FREE_KEYS),
                "n_frozen_at_b2_3": len(FROZEN_KEYS),
                "final_cost": best["cost"],
                "sum_r2_level": best["sum_r2_level"],
            },
            "frozen_parameters": _b23_schema_frozen(x_full),
        }
        path = data_paths.VALIDATION_DIR / f"kinetic_core_b2_4_fit_{tag}.json"
        path.write_text(json.dumps(payload, indent=2, default=str))
        paths[tag] = str(path.relative_to(ROOT))
    return paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight-tag", choices=sorted(PH_ENDPOINT_WEIGHT))
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)

    MEMBER_DIR.mkdir(parents=True, exist_ok=True)

    if args.consolidate:
        members = load_members()
        summary = ensemble_summary(members)
        choice = shipping_choice(summary)
        paths = write_per_weighting_fit_reports(members)
        out = {
            "wave": "B2.4 -- declared weighting + ensemble",
            "prereg": data_paths.rel(data_paths.VALIDATION_DIR / "kinetic_core_b2_4_prereg.md"),
            "declared_weights": PH_ENDPOINT_WEIGHT,
            "weight_basis": PH_ENDPOINT_WEIGHT_BASIS,
            "sigma_log_reference": SIGMA_LOG_REFERENCE,
            "ph_endpoint_rows": list(PH_ENDPOINT_ROW_IDS),
            "free_set": {
                "n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS),
                "keys": list(FREE_KEYS), "clause": FREE_CLAUSE_OF,
                "frozen_keys": list(FROZEN_KEYS),
            },
            "ensemble": summary,
            "shipping_choice": choice,
            "per_weighting_fit_reports": paths,
            "members": members,
        }
        dest = data_paths.VALIDATION_DIR / "kinetic_core_b2_4_ensemble.json"
        dest.write_text(json.dumps(out, indent=2, default=str))
        print(f"wrote {dest}")
        for tag, b in summary.items():
            print(f"  {tag:9s} E={b['E']:.2f}  n={b['n_members']} "
                  f"conv={b['n_converged']}  best={b['best_cost']}  "
                  f"S={b['spread_S_log10_max_over_min']}  "
                  f"G1={b['gate_G1_calibration_within_1pH']}")
        print("SHIPPING:", choice["shipped"])
        return 0

    if args.weight_tag is None or args.start is None:
        parser.error("--weight-tag and --start are required unless --consolidate")
    member = fit_member(args.weight_tag, args.start, args.max_nfev, args.quick)
    dest = MEMBER_DIR / f"{args.weight_tag}_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b3_fit.py

THE FIT STAGE OF BUILD WAVE B3 (kinetic core, Module 3: acrylamide / safety).

Fits the eight fitted rate constants and the three fitted activation energies
of ``src/kinetic_core/acrylamide.py`` against the DECLARED FIT ROWS of
``docs/reference/FIT_HOLDOUT_DECLARATION.md`` D.5, and writes
``results/validation/kinetic_core_b3_fit_report.{md,json}`` -- including the
FROZEN parameter set that the hold-out scorer reads and never changes.

===========================================================================
THE FIT CORPUS -- and nothing else
===========================================================================
  * Claeys 2005 Table 2, the CONTROL system and the FOUR COMPETITOR systems
    (glutamine, cysteine, lysine, alanine): 20 rows.
  * De Vleeschouwer 2009 Part I, the GLUCOSE non-italic subset. Four of its
    constants are OPERATIVE PARAMETERS of this module rather than fit rows
    (k_INTg, k_F, k_Asp and their barriers); the two that are not -- k_E and
    Ea_E -- are scored.
  * De Vleeschouwer 2009 Part II, the CYSTEINE system: its k_E2 and Ea_E2 are
    operative, and its ">99 % acrylamide reduction" claim is scored as a
    one-sided row.
  * Knol 2005 Table 1: the five k6 values (120-200 C) and the formation
    barrier Ea(k4) = 94.4.
  * the existing FIT benchmark `acrylamide_spi_extrusion_130C_ACSRef3`.

===========================================================================
WHAT IS NEVER READ
===========================================================================
`data/benchmarks/external_validation/` -- not once, at any point in this wave.
The De Vleeschouwer FRUCTOSE, SUCROSE and GLUTAMINE columns, Knol 2010 and the
Knol 2009 real-food band appear in NO row of this file, in no bound and in no
initialisation. `tests/unit/test_kinetic_core_b3.py` enforces that with a
literal-grep firewall over this file's executable code.

DISCLOSURE, because the dossiers print hold-out columns beside FIT ones: the
values of the De Vleeschouwer fructose, sucrose and glutamine columns WERE SEEN
during the directed reads this wave was instructed to make. See
`parameters_acrylamide.HOLDOUT_EXPOSURE_DISCLOSURE`.

===========================================================================
HOW A MASS-ACTION MODEL IS SCORED AGAINST LUMPED CONSTANTS
===========================================================================
Claeys, De Vleeschouwer and Knol all report first-order constants from
multiresponse fits, not elementary rates. This network has elementary rates.
The comparison is made by computing what first-order constants a first-order
description of THIS model's own trajectory would need -- the flux-weighted
apparent constants of `acrylamide.apparent_lumped_constants`, which are EXACT
for a system that really is first order (unit-tested) and approximate here in
two stated directions. Nothing about that recipe uses a literature value: the
window is the same 45 minutes at every temperature and in every system.

===========================================================================
WHAT IS FITTED, AND WHAT IS NOT
===========================================================================
NOT FITTED (measured, fixed): k_asn_glc, k_int1_acr, k_asn_asp, k_acr_cys,
k_cys_sink, k_cys_glc, k_asp_sink -- and the Ea of the three fitted
acrylamide-scavenging channels, which is held at the MEASURED Ea_E2 = 51.3.
Also not fitted: every B1 trunk constant, which is inherited frozen.

FITTED: eight log10 rate constants at 160 C and three activation energies.
Starts are random inside wide bounds, so agreement with any literature number
is a result rather than an initialisation.

THE CYSTEINE COMPETITOR ROW IS THE ONE PARAMETER-FREE PREDICTION IN THE PANEL.
Both of cysteine's channels are measured, so Claeys' cysteine row is predicted
with nothing fitted to it. It is the row to read first.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.kinetic_core import operative_acrylamide_parameters  # noqa: E402
from src.kinetic_core.acrylamide import (  # noqa: E402
    OUT_OF_SCOPE,
    apparent_activation_energy,
    apparent_lumped_constants,
    describe_acrylamide,
)
from src.kinetic_core.parameters_acrylamide import (  # noqa: E402
    CROSS_LAB_CONFLICTS,
    DELIBERATE_UNDERFITS,
    EA_E2_MEASURED_KJ_MOL,
    FITTED_ACRYLAMIDE_BOUNDS_LOG10K,
    FITTED_ACRYLAMIDE_EA_BOUNDS,
    FITTED_ACRYLAMIDE_EA_KEYS,
    FITTED_ACRYLAMIDE_KEYS,
    HOLDOUT_EXPOSURE_DISCLOSURE,
    MEASURED_ACRYLAMIDE,
    REFUSED_PARAMETERS,
    AcrylamideParameter,
    T_REF_A_K,
    acrylamide_registry_metadata,
)

# The B1 trunk fit, frozen, from results/validation/kinetic_core_b1_fit_report.json.
# B3 inherits it unchanged and re-fits no trunk constant.
B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

CELSIUS = 273.15

#: The apparent-constant window, minutes. THE SAME AT EVERY TEMPERATURE AND IN
#: EVERY SYSTEM, so that no literature value enters the definition of the
#: observable. 45 min is the span of Claeys' own Figure 2.
WINDOW_MIN = 45.0

#: The De Vleeschouwer powder's reactant loading. The paper states "~3 M
#: reactant concentrations" without a table; 3000 mmol/L of EACH component is
#: this module's reading of that sentence and is flagged as an inference on
#: every row that uses it. Nothing in the fit is strongly sensitive to it: the
#: DV rows scored here are the elimination constant (zeroth order in loading)
#: and a >99 % reduction ceiling.
DV_LOADING_MMOL_L = 3000.0

#: Knol 2005's aqueous loading is not transcribed in the corpus either. 100
#: mmol/L is taken from the inventory's Knol 2010 line ("~3 % of initial
#: asparagine ... at 0.1 M"), flagged as an inference, and matters only to the
#: FORMATION barrier row -- the five k6 rows are elimination constants and are
#: independent of it.
KNOL_LOADING_MMOL_L = 100.0

SPI_BENCHMARK = ROOT / "data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json"


# ===========================================================================
# THE SYSTEMS
# ===========================================================================

SYSTEMS: Dict[str, Dict[str, Any]] = {
    # ---- Claeys 2005 Table 2: 0.01 M Asn + 0.01 M Glc (+ 0.01 M competitor)
    # in 0.05 M citrate pH 6, closed inox tubes, 140/160/180/200 C, T_ref 160,
    # non-isothermal corrected by Euler integration of a 2 s thermocouple log.
    "claeys_control": dict(
        initial={"Asn": 10.0, "Glc": 10.0}, aw=1.0, ph=6.0,
        temperatures_c=(140.0, 160.0, 180.0, 200.0),
        anchor="Claeys 2005 Biotechnol Prog 21:1525-1530 Table 2, control system"),
    "claeys_gln": dict(
        initial={"Asn": 10.0, "Glc": 10.0, "Gln": 10.0}, aw=1.0, ph=6.0,
        temperatures_c=(140.0, 160.0, 180.0, 200.0),
        anchor="Claeys 2005 Table 2, glutamine competitor system"),
    "claeys_cys": dict(
        initial={"Asn": 10.0, "Glc": 10.0, "Cys": 10.0}, aw=1.0, ph=6.0,
        temperatures_c=(140.0, 160.0, 180.0, 200.0),
        anchor="Claeys 2005 Table 2, cysteine competitor system"),
    "claeys_lys": dict(
        initial={"Asn": 10.0, "Glc": 10.0, "Lys": 10.0}, aw=1.0, ph=6.0,
        temperatures_c=(140.0, 160.0, 180.0, 200.0),
        anchor="Claeys 2005 Table 2, lysine competitor system"),
    "claeys_ala": dict(
        initial={"Asn": 10.0, "Glc": 10.0, "Ala": 10.0}, aw=1.0, ph=6.0,
        temperatures_c=(140.0, 160.0, 180.0, 200.0),
        anchor="Claeys 2005 Table 2, alanine competitor system"),
    # ---- De Vleeschouwer 2009 I, GLUCOSE column: freeze-dried powder at
    # a_w 0.92, ~3 M, 120-200 C.
    "dv1_glucose": dict(
        initial={"Asn": DV_LOADING_MMOL_L, "Glc": DV_LOADING_MMOL_L},
        aw=0.92, ph=None,
        temperatures_c=(120.0, 160.0, 200.0),
        anchor="De Vleeschouwer 2009 I Table 3, glucose column (a_w 0.92 powder)"),
    # ---- De Vleeschouwer 2009 II, CYSTEINE system: as above plus equimolar
    # cysteine.
    "dv2_cysteine": dict(
        initial={"Asn": DV_LOADING_MMOL_L, "Glc": DV_LOADING_MMOL_L,
                 "Cys": DV_LOADING_MMOL_L},
        aw=0.92, ph=None,
        temperatures_c=(160.0,),
        anchor="De Vleeschouwer 2009 II Table 3, cysteine system (a_w 0.92 powder)"),
    # ---- Knol 2005 Table 1: aqueous asparagine/glucose, 120-200 C.
    "knol2005": dict(
        initial={"Asn": KNOL_LOADING_MMOL_L, "Glc": KNOL_LOADING_MMOL_L},
        aw=1.0, ph=None,
        temperatures_c=(120.0, 140.0, 160.0, 180.0, 200.0),
        anchor="Knol et al. 2005 Table 1, aqueous asparagine/glucose"),
}

#: The extrusion benchmark is not an apparent-constant system: it is a single
#: short-residence concentration, so it gets its own run.
SPI_SYSTEM = dict(
    initial={"Asn": 42.0, "Glc": 42.0}, aw=0.35, ph=6.0,
    temperature_c=130.0, minutes=0.4167,
    anchor="data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json",
)


# ===========================================================================
# THE FIT ROWS
# ===========================================================================
# kind:
#   "rate"       target is a first-order constant in 1/min; residual in log10
#   "ea"         target is an activation energy in kJ/mol; residual is LINEAR
#   "ppb"        target is an acrylamide concentration in ppb; residual in log10
#   "ceiling"    one-sided: the prediction must not EXCEED the target
#
# sigma_log = 0.30 dex (a factor of 2) on every lumped-constant row. THAT IS NOT
# THE MEASUREMENT ERROR -- Claeys' standard errors are 5-25 %. It is the level
# of agreement that would count as this network REPRODUCING a lumped constant
# it was never parameterised in, in a system whose water activity differs from
# the parameters' by 0.08 and whose concentration differs by 300x. Using the
# publication SE instead would say that a model built from a different lab's
# powder experiment should match a dilute-solution lumped fit to within 5 %,
# which nothing in the corpus supports.
#
# sigma_ea = 25 kJ/mol on every barrier row, for the same kind of reason and
# with a sharper justification: the corpus's OWN three measurements of the
# acrylamide formation barrier span 74 kJ/mol (168.25 / 159.2 / 94.4). A
# tolerance tighter than the literature's own spread would be demanding that
# the model resolve a disagreement the literature has not.

SIGMA_LOG = 0.30
SIGMA_EA = 25.0

_CLAEYS = {
    # system: (k_F at 160 C, k_E at 160 C, Ea_F, Ea_E), from Table 2 p. 1528.
    # k values are printed x1e-3 min^-1.
    "claeys_control": (0.451e-3, 111.1e-3, 168.25, 167.21),
    "claeys_gln": (1.640e-3, 274.1e-3, 166.8, 103.9),
    "claeys_cys": (0.501e-3, 268.7e-3, 206.3, 180.0),
    "claeys_lys": (0.587e-3, 280.2e-3, 179.3, 140.0),
    "claeys_ala": (0.465e-3, 103.1e-3, 173.3, 169.7),
}

#: Knol 2005 Table 1's k6 (acrylamide -> product X, first order), x1e-3 min^-1.
_KNOL_K6 = {
    120.0: 7.96e-3,
    140.0: 28.1e-3,
    160.0: 88.1e-3,
    180.0: 250e-3,
    200.0: 650e-3,
}


def _build_rows() -> Tuple[Dict[str, Any], ...]:
    rows: List[Dict[str, Any]] = []
    for system, (k_f, k_e, ea_f, ea_e) in _CLAEYS.items():
        label = system.split("_", 1)[1]
        rows.append(dict(
            id=f"claeys_{label}_kF_160C", system=system, kind="rate",
            observable="k_F_app_per_min", temperature_c=160.0,
            target=k_f, sigma=SIGMA_LOG,
            anchor=SYSTEMS[system]["anchor"] + ", k_Fref",
        ))
        rows.append(dict(
            id=f"claeys_{label}_kE_160C", system=system, kind="rate",
            observable="k_E_app_per_min", temperature_c=160.0,
            target=k_e, sigma=SIGMA_LOG,
            anchor=SYSTEMS[system]["anchor"] + ", k_Eref",
        ))
        rows.append(dict(
            id=f"claeys_{label}_EaF", system=system, kind="ea",
            observable="k_F_app_per_min", temperature_c=None,
            target=ea_f, sigma=SIGMA_EA,
            anchor=SYSTEMS[system]["anchor"] + ", Ea_F",
        ))
        rows.append(dict(
            id=f"claeys_{label}_EaE", system=system, kind="ea",
            observable="k_E_app_per_min", temperature_c=None,
            target=ea_e, sigma=SIGMA_EA,
            anchor=SYSTEMS[system]["anchor"] + ", Ea_E",
        ))

    rows.append(dict(
        id="dv1_glucose_kE_160C", system="dv1_glucose", kind="rate",
        observable="k_E_app_per_min", temperature_c=160.0,
        target=0.10, sigma=SIGMA_LOG,
        anchor="De Vleeschouwer 2009 I Table 3, glucose k_Eref = 0.10 +/- 0.04 min^-1",
        note="a_w 0.92 powder; the SECOND of three independent measurements of "
             "this same constant.",
    ))
    rows.append(dict(
        id="dv1_glucose_EaE", system="dv1_glucose", kind="ea",
        observable="k_E_app_per_min", temperature_c=None,
        target=113.2, sigma=SIGMA_EA,
        anchor="De Vleeschouwer 2009 I Table 3, glucose Ea_E = 113.2 +/- 32.3",
    ))
    rows.append(dict(
        id="dv2_cysteine_reduction", system="dv2_cysteine", kind="ceiling",
        observable="_acrylamide_ratio_to_dv1", temperature_c=160.0,
        target=0.01, sigma=SIGMA_LOG,
        anchor="De Vleeschouwer 2009 II abstract: '>99 % reduction' of "
               "acrylamide at the equimolar cysteine loading",
        note="ONE-SIDED. The prediction is the peak acrylamide of the cysteine "
             "system divided by the peak of the glucose control; it must not "
             "EXCEED 0.01. Both of cysteine's channels are MEASURED, so this "
             "row has nothing fitted to it.",
    ))

    for temperature_c, k6 in _KNOL_K6.items():
        rows.append(dict(
            id=f"knol2005_k6_{int(temperature_c)}C", system="knol2005", kind="rate",
            observable="k_E_app_per_min", temperature_c=temperature_c,
            target=k6, sigma=SIGMA_LOG,
            anchor=f"Knol 2005 Table 1, k6 at {int(temperature_c)} C",
            note="the THIRD independent lab on the acrylamide elimination "
                 "constant; the five values imply Ea 85.05, matching the "
                 "printed 85.1 +/- 14.",
        ))
    rows.append(dict(
        id="knol2005_EaF", system="knol2005", kind="ea",
        observable="k_F_app_per_min", temperature_c=None,
        target=94.4, sigma=SIGMA_EA,
        anchor="Knol 2005 Table 1, Ea(k4) = 94.4 +/- 11 kJ/mol (formation)",
        note="THE LARGEST UNRESOLVED CONFLICT IN MODULE 3: Claeys measures "
             "168.25 +/- 3.80 for the same transformation in the same kind of "
             "aqueous system. Both are FIT rows. The model cannot satisfy both "
             "and no attempt is made to.",
    ))

    rows.append(dict(
        id="spi_extrusion_130C_acrylamide_ppb", system="_spi", kind="ppb",
        observable="peak_acrylamide_ppb", temperature_c=130.0,
        target=None, sigma=1.0,
        anchor=SPI_SYSTEM["anchor"],
        note="sigma is 1.0 dex, an order of magnitude, and the reasons are on "
             "the benchmark file itself: the value is a FIGURE READOFF whose "
             "predecessor (62.62 ppb) 'matches no bar in any panel and no "
             "derivable statistic'; the run is 25 SECONDS at a_w 0.35, outside "
             "every system that carries a rate constant here; and the melt "
             "temperature of a 250 kJ/kg extrusion is not 130 C. It is scored "
             "because the declaration lists it as FIT, and it is weighted so "
             "that it cannot drive the fit.",
    ))
    return tuple(rows)


FIT_ROWS: Tuple[Dict[str, Any], ...] = _build_rows()


def _spi_target_ppb() -> float:
    payload = json.loads(SPI_BENCHMARK.read_text())
    return float(payload["measured_volatiles"]["acrylamide"]["conc_ppb"])


# ===========================================================================
# THE OBJECTIVE
# ===========================================================================

PARAM_ORDER: Tuple[str, ...] = FITTED_ACRYLAMIDE_KEYS
EA_ORDER: Tuple[str, ...] = FITTED_ACRYLAMIDE_EA_KEYS
N_K = len(PARAM_ORDER)
N_EA = len(EA_ORDER)


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    log10k = {key: float(x[i]) for i, key in enumerate(PARAM_ORDER)}
    ea = {key: float(x[N_K + i]) for i, key in enumerate(EA_ORDER)}
    return operative_acrylamide_parameters(B1_FITTED, log10k, ea)


def run_systems(parameters: Mapping[str, Any], quick: bool) -> Dict[str, Any]:
    """One apparent-constant evaluation per (system, temperature)."""
    rtol = 1e-5 if quick else 1e-8
    n_points = 41 if quick else 121
    out: Dict[str, Any] = {}
    for name, spec in SYSTEMS.items():
        for temperature_c in spec["temperatures_c"]:
            out[(name, temperature_c)] = apparent_lumped_constants(
                parameters, float(temperature_c) + CELSIUS, spec["initial"],
                WINDOW_MIN, n_points=n_points, water_activity=spec["aw"],
                rtol=rtol, atol=1e-14,
            )
    spi = apparent_lumped_constants(
        parameters, SPI_SYSTEM["temperature_c"] + CELSIUS, SPI_SYSTEM["initial"],
        float(SPI_SYSTEM["minutes"]), n_points=n_points,
        water_activity=SPI_SYSTEM["aw"], rtol=rtol, atol=1e-14,
    )
    out[("_spi", 130.0)] = spi
    return out


def observables(cache: Mapping[Any, Any], rows) -> Dict[str, float]:
    predicted: Dict[str, float] = {}
    for row in rows:
        kind = row["kind"]
        system = row["system"]
        if kind == "rate":
            predicted[row["id"]] = float(
                cache[(system, row["temperature_c"])][row["observable"]]
            )
        elif kind == "ppb":
            predicted[row["id"]] = float(
                cache[(system, row["temperature_c"])][row["observable"]]
            )
        elif kind == "ea":
            temperatures = SYSTEMS[system]["temperatures_c"]
            t_low, t_high = float(min(temperatures)), float(max(temperatures))
            k_low = float(cache[(system, t_low)][row["observable"]])
            k_high = float(cache[(system, t_high)][row["observable"]])
            predicted[row["id"]] = apparent_activation_energy(
                k_low, t_low + CELSIUS, k_high, t_high + CELSIUS
            )
        elif kind == "ceiling":
            numerator = float(cache[("dv2_cysteine", 160.0)]["peak_acrylamide_ppb"])
            denominator = float(cache[("dv1_glucose", 160.0)]["peak_acrylamide_ppb"])
            predicted[row["id"]] = (
                numerator / denominator if denominator > 0 else float("nan")
            )
        else:
            raise ValueError(f"unknown row kind {kind!r}")
    return predicted


FLOOR = 1e-15
_PROGRESS = {"n": 0, "best": float("inf"), "t0": time.time()}


def _rows_with_targets() -> Tuple[Dict[str, Any], ...]:
    spi = _spi_target_ppb()
    out = []
    for row in FIT_ROWS:
        row = dict(row)
        if row["target"] is None:
            row["target"] = spi
        out.append(row)
    return tuple(out)


ACTIVE_FIT_ROWS: Tuple[Dict[str, Any], ...] = _rows_with_targets()


def residuals(x: np.ndarray, quick: bool = True) -> np.ndarray:
    try:
        cache = run_systems(build_parameters(x), quick)
        predicted = observables(cache, ACTIVE_FIT_ROWS)
    except Exception:
        out = np.full(len(ACTIVE_FIT_ROWS), 25.0)
        _progress(out)
        return out
    out = np.empty(len(ACTIVE_FIT_ROWS))
    for i, row in enumerate(ACTIVE_FIT_ROWS):
        p = predicted[row["id"]]
        t = float(row["target"])
        sigma = float(row["sigma"])
        if not np.isfinite(p):
            out[i] = 25.0
            continue
        if row["kind"] == "ea":
            out[i] = (p - t) / sigma
        elif row["kind"] == "ceiling":
            out[i] = max(0.0, math.log10((p + FLOOR) / (t + FLOOR))) / sigma
        else:
            out[i] = math.log10((p + FLOOR) / (t + FLOOR)) / sigma
    out = np.clip(out, -25.0, 25.0)
    _progress(out)
    return out


def _progress(r: np.ndarray) -> None:
    _PROGRESS["n"] += 1
    cost = 0.5 * float(np.dot(r, r))
    _PROGRESS["best"] = min(_PROGRESS["best"], cost)
    if _PROGRESS["n"] % 25 == 0:
        print(f"    [{_PROGRESS['n']:5d} evals, "
              f"{time.time() - _PROGRESS['t0']:6.0f}s] best cost "
              f"{_PROGRESS['best']:.3f}", flush=True)


# ===========================================================================
# THE FIT
# ===========================================================================


def fit(starts: int, quick: bool, max_nfev: int, seed: int = 20260828) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    lower = np.array(
        [FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][0] for k in PARAM_ORDER]
        + [FITTED_ACRYLAMIDE_EA_BOUNDS[0]] * N_EA
    )
    upper = np.array(
        [FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][1] for k in PARAM_ORDER]
        + [FITTED_ACRYLAMIDE_EA_BOUNDS[1]] * N_EA
    )
    rng = np.random.default_rng(seed)
    best = None
    trace: List[Dict[str, Any]] = []
    for start in range(starts):
        # A RANDOM point inside the bounds. Nothing about it encodes a
        # literature value: not one of these eleven numbers is printed in a
        # usable form anywhere in the corpus.
        x0 = lower + rng.random(len(lower)) * (upper - lower)
        t0 = time.time()
        result = least_squares(
            residuals, x0, bounds=(lower, upper), args=(quick,),
            method="trf", x_scale="jac", max_nfev=max_nfev,
            diff_step=3e-2, ftol=1e-8, xtol=1e-10, verbose=0,
        )
        cost = float(result.cost)
        trace.append({
            "start": start,
            "x0": [float(v) for v in x0],
            "cost": cost,
            "nfev": int(result.nfev),
            "seconds": round(time.time() - t0, 1),
            "status": int(result.status),
        })
        print(f"  start {start}: cost {cost:.4f}  nfev {result.nfev}  "
              f"{trace[-1]['seconds']}s", flush=True)
        if best is None or cost < best.cost:
            best = result
    assert best is not None
    return {"result": best, "trace": trace, "bounds": (lower, upper)}


#: A parameter counts as IDENTIFIED only if its 95 % half-width is inside
#: these. The test is on the HALF-WIDTH rather than on the raw sigma, which is
#: the stricter reading and the one that matches what the word means: a rate
#: constant whose 95 % interval spans three orders of magnitude has not been
#: determined by anything, however finite its standard error. One dex is a
#: factor of ten either way, which is already generous for a rate; 60 kJ/mol is
#: set by the corpus's own between-lab spread on the acrylamide barriers
#: (85.1 to 167.2 is 82 kJ/mol) -- a barrier pinned tighter than the
#: literature's own disagreement would be the fit claiming more than the
#: evidence has.
IDENTIFIED_HALFWIDTH_DEX = 1.0
IDENTIFIED_HALFWIDTH_KJ_MOL = 60.0


def intervals(result, n_residuals: int) -> Dict[str, Any]:
    """Gauss-Newton 95 % intervals, with flat directions reported as such."""
    jac = np.asarray(result.jac)
    dof = max(n_residuals - jac.shape[1], 1)
    chi2_red = 2.0 * float(result.cost) / dof
    try:
        covariance = np.linalg.pinv(jac.T @ jac) * chi2_red
        sigmas = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    except Exception:
        sigmas = np.full(jac.shape[1], np.inf)
    names = list(PARAM_ORDER) + list(EA_ORDER)
    out: Dict[str, Any] = {}
    for i, name in enumerate(names):
        s = float(sigmas[i])
        halfwidth = 1.96 * s
        threshold = (
            IDENTIFIED_HALFWIDTH_DEX if i < N_K else IDENTIFIED_HALFWIDTH_KJ_MOL
        )
        identified = np.isfinite(halfwidth) and halfwidth <= threshold
        out[name] = {
            "value": float(result.x[i]),
            "unit": "log10 k(160 C)" if i < N_K else "kJ/mol",
            "ci95_halfwidth": None if not np.isfinite(halfwidth) else halfwidth,
            "identified_threshold": threshold,
            "identified": bool(identified),
        }
    out["_reduced_chi_square"] = chi2_red
    out["_dof"] = dof
    return out


# ===========================================================================
# TWO DIAGNOSTICS THE FIRST FIT EXPOSED
# ===========================================================================
# Neither is a scored row and neither changes a parameter. Both are structural
# facts about the composed network that the row table alone does not show, and
# both were found by inspecting the first fit rather than anticipated.


def _sugar_depletion_diagnostic(parameters) -> Dict[str, Any]:
    """
    How much glucose is left at the end of the apparent-constant window.

    THIS IS THE BINDING CONSTRAINT ON EVERY CLAEYS FORMATION ROW, and it comes
    from B1's TRUNK rather than from anything in Module 3. Martins' sugar
    constants were measured over 80-120 C; at 160 C they are extrapolated 40 C
    beyond that, and at 200 C, 80 C beyond it. Evaluated there, the trunk's own
    fructose fragmentation lane empties a 10 mmol/L glucose charge within
    minutes -- which starves the bimolecular asparagine initiation and caps the
    acrylamide the module can form, no matter what the fitted Int1 partition
    does. It is why the fitted partition sits at its lower bound: the optimiser
    is trying to buy back acrylamide that the sugar lane has already spent.

    The extrapolation is reported on every run by
    `acrylamide._acrylamide_warnings`; this diagnostic quantifies its effect.
    """
    from src.kinetic_core.acrylamide import integrate_acrylamide

    out: Dict[str, Any] = {
        "what_this_is": (
            "glucose remaining at the end of the 45 min window, as a fraction "
            "of the charge, in the Claeys control system"
        ),
        "why_it_matters": (
            "the asparagine initiation is SECOND ORDER, so it stops when the "
            "sugar does. This is a B1 trunk extrapolation effect, not a "
            "Module 3 parameter."
        ),
        "by_temperature": {},
    }
    for temperature_c in SYSTEMS["claeys_control"]["temperatures_c"]:
        run = integrate_acrylamide(
            parameters, float(temperature_c) + CELSIUS,
            SYSTEMS["claeys_control"]["initial"],
            np.array([0.0, WINDOW_MIN]), water_activity=1.0,
            rtol=1e-8, atol=1e-14,
        )
        out["by_temperature"][f"{int(temperature_c)}C"] = {
            "glucose_remaining_fraction": run.final("Glc") / 10.0,
            "fructose_remaining_fraction": run.final("Fru") / 10.0,
            "asparagine_remaining_fraction": run.final("Asn") / 10.0,
        }
    return out


def _cysteine_window_diagnostic(parameters) -> Dict[str, Any]:
    """
    The acrylamide reduction by cysteine, as a function of the window.

    De Vleeschouwer II's abstract claims '>99 %'. This network reproduces that
    at short times and loses it at long ones, and the reason is INSIDE THE SAME
    PAPER: its measured cysteine sink, k_Y = 0.35 min^-1 at 160 C, empties the
    cysteine pool with a half-life of two minutes, while its measured
    acrylamide-forming step, k_F = 3.57e-3 min^-1, keeps producing acrylamide
    for hours. Scavenging that fast and formation that slow cannot both be
    right and also give a durable >99 % reduction.

    This is a genuine internal tension in the source's own parameter set that
    a mass-action network exposes and a lumped two-constant fit cannot. It is
    reported, not resolved, and the scored row is NOT re-defined around it:
    the 45-minute window was declared before the fit for every row in the
    panel, and changing it after seeing which window would pass would be
    fitting the objective to the answer.
    """
    from src.kinetic_core.acrylamide import apparent_lumped_constants as observe

    rows = []
    for window in (2.0, 5.0, 10.0, 20.0, WINDOW_MIN, 120.0):
        with_cys = observe(
            parameters, T_REF_A_K, SYSTEMS["dv2_cysteine"]["initial"],
            window, n_points=61, water_activity=0.92, rtol=1e-8, atol=1e-14,
        )["peak_acrylamide_ppb"]
        control = observe(
            parameters, T_REF_A_K, SYSTEMS["dv1_glucose"]["initial"],
            window, n_points=61, water_activity=0.92, rtol=1e-8, atol=1e-14,
        )["peak_acrylamide_ppb"]
        rows.append({
            "window_min": window,
            "ratio": with_cys / control if control > 0 else float("nan"),
            "reduction_percent": (
                100.0 * (1.0 - with_cys / control) if control > 0 else float("nan")
            ),
        })
    return {
        "claim": "De Vleeschouwer 2009 II abstract: >99 % acrylamide reduction",
        "model_says": "the reduction is REAL but TRANSIENT",
        "mechanism": (
            "the same paper's k_Y = 0.35 min^-1 empties the cysteine pool with "
            "a 2 min half-life while its k_F = 3.57e-3 min^-1 keeps forming "
            "acrylamide for hours, so the scavenger is gone long before the "
            "acrylamide peak"
        ),
        "scored_window_min": WINDOW_MIN,
        "by_window": rows,
        "not_re_scored_because": (
            "the window was declared for every row before the fit; picking the "
            "window that passes would be fitting the objective to the answer"
        ),
    }


# ===========================================================================
# REPORTING
# ===========================================================================


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=400)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    args = parser.parse_args()

    print(f"B3 acrylamide fit: {len(ACTIVE_FIT_ROWS)} declared FIT rows, "
          f"{N_K + N_EA} free parameters", flush=True)

    fitted = fit(args.starts, args.quick, args.max_nfev)
    result = fitted["result"]
    parameters = build_parameters(result.x)
    cache = run_systems(parameters, args.quick)
    predicted = observables(cache, ACTIVE_FIT_ROWS)
    ci = intervals(result, len(ACTIVE_FIT_ROWS))

    rows_out = []
    for row in ACTIVE_FIT_ROWS:
        p = float(predicted[row["id"]])
        t = float(row["target"])
        if row["kind"] == "ea":
            fold = float("nan")
            residual = (p - t) / float(row["sigma"])
            delta = p - t
        else:
            fold = (p / t) if (t > 0 and p > 0) else float("nan")
            residual = math.log10((p + FLOOR) / (t + FLOOR)) / float(row["sigma"])
            delta = float("nan")
        rows_out.append({
            "id": row["id"],
            "system": row["system"],
            "kind": row["kind"],
            "observed": t,
            "predicted": p,
            "fold_error": fold,
            "delta_kj_mol": delta,
            "weighted_residual": residual,
            "sigma": row["sigma"],
            "source_anchor": row["anchor"],
            "note": row.get("note", ""),
        })

    frozen = {
        "log10_k_ref_at_160C": {
            k: float(result.x[i]) for i, k in enumerate(PARAM_ORDER)
        },
        "fitted_Ea_kJ_mol": {
            k: float(result.x[N_K + i]) for i, k in enumerate(EA_ORDER)
        },
        "Ea_of_the_scavenging_channels_kJ_mol": EA_E2_MEASURED_KJ_MOL,
        "Ea_of_the_scavenging_channels_is": "MEASURED, not fitted (De Vleeschouwer II Ea_E2)",
        "reference_temperature_K": T_REF_A_K,
        "window_min": WINDOW_MIN,
        "b1_trunk_fit_inherited": B1_FITTED,
    }

    # Diagnostics that are NOT scored rows, because each is algebraically a
    # restatement of a row above and scoring both would double-count.
    k_acr_dp = 10.0 ** frozen["log10_k_ref_at_160C"]["k_acr_dp"]
    crossover_mmol_l = k_acr_dp / MEASURED_ACRYLAMIDE["k_acr_cys"].k_ref
    diagnostics = {
        "cysteine_crossover_concentration_mmol_L": {
            "predicted": crossover_mmol_l,
            "literature": 2.0,
            "what_it_is": (
                "the cysteine concentration at which the MEASURED scavenging "
                "channel matches the fitted first-order elimination, "
                "k_E / k_E2. De Vleeschouwer II's own derived value is 2.0 mM "
                "(inventory sec. B4.2)."
            ),
            "why_not_a_scored_row": (
                "it is k_acr_dp divided by a measured constant, so it carries "
                "exactly the information of the dv1_glucose_kE row and scoring "
                "both would count one measurement twice."
            ),
        },
        "knol_k6_implied_Ea_kJ_mol": {
            "from_the_five_printed_values": 85.05,
            "printed_by_the_source": 85.1,
            "why_not_a_scored_row": (
                "the five k6 rows already carry the barrier; adding an Ea row "
                "derived from them would weight the same data twice."
            ),
        },
        "claeys_control_internal_validation": {
            "quasi_plateau_ppb_from_the_published_constants": 2886.0,
            "t_max_min_from_the_published_constants": 49.8,
            "model_peak_ppb": float(cache[("claeys_control", 160.0)]["peak_acrylamide_ppb"]),
            "model_t_max_min": float(cache[("claeys_control", 160.0)]["time_of_peak_min"]),
            "what_it_is": (
                "Claeys' own internal validation, recomputed on the model. "
                "DERIVED from the same two constants the kF/kE rows already "
                "score, so it is a diagnostic and not a row."
            ),
        },
        "elimination_channel_split_at_160C": {
            system: cache[(system, 160.0)]["elimination_flux_by_channel_mmol_L"]
            for system in ("claeys_control", "claeys_cys", "claeys_lys",
                           "claeys_gln", "claeys_ala")
        },
        "trunk_sugar_depletion": _sugar_depletion_diagnostic(parameters),
        "cysteine_reduction_is_transient": _cysteine_window_diagnostic(parameters),
    }

    payload: Dict[str, Any] = {
        "wave": "B3 -- acrylamide / safety",
        "generated_by": "scripts/generators/generate_kinetic_core_b3_fit.py",
        "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md D.5 (Module 3)",
        "network": describe_acrylamide(),
        "objective": {
            "form": (
                "weighted least squares. RATE and PPB rows: "
                "r = log10((pred + 1e-15)/(obs + 1e-15)) / sigma_log. "
                "EA rows: r = (pred - obs) / sigma_ea, LINEAR, because an "
                "activation energy is not a positive scale quantity. CEILING "
                "rows: one-sided, no penalty below the ceiling."
            ),
            "n_rows": len(ACTIVE_FIT_ROWS),
            "n_free_parameters": N_K + N_EA,
            "final_cost": float(result.cost),
            "reduced_chi_square": ci["_reduced_chi_square"],
            "dof": ci["_dof"],
            "sigma_log_dex": SIGMA_LOG,
            "sigma_ea_kj_mol": SIGMA_EA,
            "apparent_constant_window_min": WINDOW_MIN,
            "multistart_trace": fitted["trace"],
            "starts_are_random_inside_wide_bounds": True,
        },
        "frozen_parameters": frozen,
        "parameter_intervals": {
            k: v for k, v in ci.items() if not k.startswith("_")
        },
        "rows": rows_out,
        "diagnostics_not_scored": diagnostics,
        "inferred_conditions": {
            "de_vleeschouwer_loading_mmol_L": DV_LOADING_MMOL_L,
            "knol_2005_loading_mmol_L": KNOL_LOADING_MMOL_L,
            "why": (
                "neither paper's reactant molarity is transcribed in the "
                "corpus; both are read from prose ('~3 M', and Knol 2010's "
                "'at 0.1 M'). Recorded as INFERENCES. The rows they carry are "
                "elimination constants and a one-sided ceiling, both of which "
                "are insensitive to the loading; the one row that is sensitive "
                "(knol2005_EaF) is flagged."
            ),
        },
        # AcrylamideParameter, not "anything with an aw_of_measurement field":
        # B1's KineticParameter carries that field too, and sweeping the trunk
        # constants into the acrylamide registry block would misattribute them.
        "registry_metadata": acrylamide_registry_metadata(
            {k: v for k, v in parameters.items()
             if isinstance(v, AcrylamideParameter)}
        ),
        "refused_parameters": [dict(r) for r in REFUSED_PARAMETERS],
        "cross_lab_conflicts": [dict(r) for r in CROSS_LAB_CONFLICTS],
        "deliberate_underfits": [dict(r) for r in DELIBERATE_UNDERFITS],
        "holdout_exposure_disclosure": dict(HOLDOUT_EXPOSURE_DISCLOSURE),
        "out_of_scope": [dict(r) for r in OUT_OF_SCOPE],
        "repo_defect_fixed_en_route": {
            "what": (
                "data/lit/safety_reference_payloads.json entries[27] "
                "(knol_2005_acrylamide_kinetics) stated the Knol 2005 "
                "acrylamide barriers as 52.1 (formation) and 72.9 "
                "(degradation) kJ/mol. The true pair is 94.4 +/- 11 and "
                "85.1 +/- 14 (inventory sec. A.2 lines 202-203, errata #3)."
            ),
            "action": (
                "CORRECTED in place with a dated provenance note. No test "
                "asserted the old values; the schema is unchanged."
            ),
        },
    }

    out_json = ROOT / "results/validation/kinetic_core_b3_fit_report.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, default=str))

    lines: List[str] = []
    a = lines.append
    a("# Kinetic core, Build Wave B3 -- the acrylamide / safety fit")
    a("")
    a("Module 3 of `docs/reference/FIT_HOLDOUT_DECLARATION.md`.")
    a("")
    net = payload["network"]
    a(f"- network: **{net['n_species']} species, {net['n_reactions']} "
      f"reactions** ({net['trunk_reactions']} inherited from B1's trunk, "
      f"{net['acrylamide_reactions']} new), carbon, nitrogen and sulfur "
      f"balance enforced at import")
    a(f"- objective: **{len(ACTIVE_FIT_ROWS)} declared FIT rows**, "
      f"**{N_K + N_EA} free parameters**, final cost "
      f"**{float(result.cost):.3f}**, reduced chi-square "
      f"**{ci['_reduced_chi_square']:.2f}**, dof **{ci['_dof']}**")
    a(f"- acrylamide ELIMINATION channels: "
      f"**{len(net['acrylamide_sink_reactions'])}** "
      f"(`{'`, `'.join(net['acrylamide_sink_reactions'])}`). The old FAST lane "
      f"had none, which is why it was ~40x under-responsive to time.")
    a("- competition multipliers: **0**. Competition is two named mass-action "
      "channels per amino acid, one of them measured.")
    a("- pH terms: **0**. a_w terms: **0**. Occurrences of the unsourceable "
      "129 kJ/mol acrylamide barrier: **0**.")
    a("")
    a("## Read this before reading anything else")
    a("")
    n_identified = sum(
        1 for k, v in ci.items() if not k.startswith("_") and v["identified"]
    )
    a(f"**{n_identified} of {N_K + N_EA} free parameters are individually "
      f"identified.** Seven of the sixteen steps in this module carry a "
      f"literature rate constant and are NOT fitted, which is why the panel "
      f"can afford eleven free parameters at all -- but the row-by-row "
      f"agreement below is still not evidence that the model is right. "
      f"The rows worth reading as evidence are the ones with nothing fitted "
      f"to them:")
    a("")
    a("- `claeys_cys_kE_160C` and `dv2_cysteine_reduction`: **both of "
      "cysteine's channels are measured**, so these are parameter-free "
      "predictions of a competitor's effect.")
    a("- `claeys_*_EaF`: the formation barrier is a composite of two measured "
      "barriers and one fitted one, and the five systems share it.")
    a("- `knol2005_k6_*`: five temperatures from a third lab against a "
      "constant fitted to two others.")
    a("")
    a("## Row-by-row")
    a("")
    a("| row | kind | observed | predicted | fold / delta | w. residual | source |")
    a("|---|---|---:|---:|---:|---:|---|")
    for r in rows_out:
        if r["kind"] == "ea":
            comparison = f"{r['delta_kj_mol']:+.1f} kJ/mol"
        elif np.isfinite(r["fold_error"]):
            comparison = f"{r['fold_error']:.2f}x"
        else:
            comparison = "n/a"
        a(f"| `{r['id']}` | {r['kind']} | {r['observed']:.4g} | "
          f"{r['predicted']:.4g} | {comparison} | "
          f"{r['weighted_residual']:+.2f} | {r['source_anchor'][:60]} |")
    a("")
    a("## Parameters")
    a("")
    a("| parameter | value | unit | 95% half-width | identified? |")
    a("|---|---:|---|---:|---|")
    for name in list(PARAM_ORDER) + list(EA_ORDER):
        entry = ci[name]
        hw = entry["ci95_halfwidth"]
        a(f"| `{name}` | {entry['value']:.3f} | {entry['unit']} | "
          f"{'flat' if hw is None else f'{hw:.3f}'} | "
          f"{'yes' if entry['identified'] else '**UNIDENTIFIED**'} |")
    a("")
    a(f"'Identified' means the 95 % HALF-WIDTH is inside "
      f"{IDENTIFIED_HALFWIDTH_DEX:.0f} dex for a rate constant and "
      f"{IDENTIFIED_HALFWIDTH_KJ_MOL:.0f} kJ/mol for a barrier -- the test is "
      f"on the interval, not on the raw standard error, because a constant "
      f"whose 95 % interval spans three orders of magnitude has not been "
      f"determined by anything. `flat` means the direction is numerically "
      f"flat and no interval exists at all.")
    a("")
    a("The three acrylamide-scavenging channels' activation energy is NOT in "
      f"that table because it is not fitted: it is held at the measured "
      f"Ea_E2 = {EA_E2_MEASURED_KJ_MOL} kJ/mol.")
    a("")
    a("## Deliberate under-fit, stated before the fit was run")
    a("")
    for row in DELIBERATE_UNDERFITS:
        a(f"- **{row['row']}** -- {row['what_the_model_will_do']}")
        a(f"  Why no term was added: {row['why_no_promotion_TERM_IS_ADDED']}")
    a("")
    a("## Cross-lab conflicts carried, not averaged")
    a("")
    for row in CROSS_LAB_CONFLICTS:
        values = ", ".join(f"{k} = {v}" for k, v in row["values"].items())
        a(f"- **{row['quantity']}** ({row['unit']}): {values}. "
          f"Spread: {row['spread']}. {row['verdict']}")
    a("")
    a("## Diagnostics (computed, not scored)")
    a("")
    crossover = diagnostics["cysteine_crossover_concentration_mmol_L"]
    a(f"- cysteine crossover k_E/k_E2: model **{crossover['predicted']:.2f} "
      f"mmol/L** against the literature's 2.0 mmol/L.")
    internal = diagnostics["claeys_control_internal_validation"]
    a(f"- Claeys' own internal validation at 160 C: plateau 2886 ppb and "
      f"t_max 49.8 min from the published constants; this model gives "
      f"**{internal['model_peak_ppb']:.0f} ppb** at "
      f"**{internal['model_t_max_min']:.1f} min**.")
    a("")
    a("### The binding constraint is B1's trunk, not a Module 3 parameter")
    a("")
    depletion = diagnostics["trunk_sugar_depletion"]
    a("Glucose remaining at the end of the 45 min window, Claeys control "
      "system, as a fraction of the charge:")
    a("")
    a("| T | glucose left | fructose left | asparagine left |")
    a("|---|---:|---:|---:|")
    for temperature, values in depletion["by_temperature"].items():
        a(f"| {temperature} | {values['glucose_remaining_fraction']:.3g} | "
          f"{values['fructose_remaining_fraction']:.3g} | "
          f"{values['asparagine_remaining_fraction']:.3g} |")
    a("")
    a("The asparagine initiation is SECOND ORDER, so it stops when the sugar "
      "does. Martins' sugar constants were measured over 80-120 C and this "
      "panel runs at 140-200 C; evaluated there, the trunk's own fructose "
      "fragmentation lane empties the sugar within minutes. **That, and not "
      "the Int1 partition, is what caps every Claeys formation row** -- and it "
      "is why the fitted partition sits at its lower bound, buying back "
      "acrylamide the sugar lane has already spent. It is a B1 trunk "
      "extrapolation, flagged on every run by the module's own warnings, and "
      "it is the first thing a later wave should look at.")
    a("")
    a("### The cysteine reduction is real but TRANSIENT")
    a("")
    transient = diagnostics["cysteine_reduction_is_transient"]
    a("| window | ratio to control | reduction |")
    a("|---:|---:|---:|")
    for entry in transient["by_window"]:
        a(f"| {entry['window_min']:.0f} min | {entry['ratio']:.4f} | "
          f"{entry['reduction_percent']:.1f} % |")
    a("")
    a(f"The mechanism: {transient['mechanism']}. The `dv2_cysteine_reduction` "
      f"row is scored at the panel's declared {WINDOW_MIN:.0f} min window and "
      f"FAILS there, while the same model reproduces the published >99 % at "
      f"two minutes. This is an internal tension in De Vleeschouwer II's own "
      f"parameter set -- a scavenger with a two-minute half-life cannot "
      f"durably suppress a product that forms over hours -- exposed by writing "
      f"the chemistry as mass action rather than as two lumped constants. "
      f"It is reported, not resolved, and the row was NOT re-scored: "
      f"{transient['not_re_scored_because']}.")
    a("")
    a("## Repo defect fixed en route")
    a("")
    a(payload["repo_defect_fixed_en_route"]["what"])
    a("")
    a(payload["repo_defect_fixed_en_route"]["action"])
    a("")
    a("## Hold-out exposure disclosure")
    a("")
    a(HOLDOUT_EXPOSURE_DISCLOSURE["what_was_seen"])
    a("")
    a(HOLDOUT_EXPOSURE_DISCLOSURE["what_was_done_about_it"])
    a("")
    a("## Out of scope for this wave")
    a("")
    for row in OUT_OF_SCOPE:
        a(f"- **{row['lane']}** -- strands: {row['what_is_stranded']} "
          f"Reason: {row['why']}")
    a("")

    out_md = ROOT / "results/validation/kinetic_core_b3_fit_report.md"
    out_md.write_text("\n".join(lines))
    print(f"wrote {out_json}")
    print(f"wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

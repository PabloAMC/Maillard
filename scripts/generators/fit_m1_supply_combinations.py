#!/usr/bin/env python
"""
fit_m1_supply_combinations.py -- plan step M1: thiol supply, R0 (the shipped structure) against C1 (the TTCA /
cysteine-Amadori route), built exactly to its pre-registration.

Pre-registration (the specification; every rule below is its, and a rule that differs from it is a bug here):
results/validation/kinetic_core_m1_supply_loss_combinations_prereg.md (v3.3, frozen at 68247583). Case
docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md; conflict CF-107; M1.0 = EX-060.

The two rivals (prereg sec. 2):
  R0  the shipped switches; B60's 16 identified coordinates + k_mft_decay (17), on B60's 52 rows + the prior.
  C1  ``sulfur.cys_amadori_route`` on (TTCA formation by detailed balance, r_pent_dpo / r_pent_tdp retired); R0's 17 +
      k_ttca_cys + the route's ea_form (19), on B60's 52 rows + the route block (route_submodel.yml less Zhang 2025b,
      the five Zhai 140 C pots and the Kang duplicate) + the prior.
Both fit under the encounter-limit ceiling on k_dimer_* (``fit.dimer_encounter_ceiling``) and the weak k_mft_decay
prior (``fit.mft_decay_prior``), on ``ph.coupling = "inline"`` (fixed-point fallback if V8 fails). Every pot runs under
a wall-clock alarm; a failed pot run is counted in the record (``failed_pot_runs``) and such a record is re-run.

Modes (run in this order; sec. 10):
    --checks [--only V1,V4]          V1, V2, V4 - V8 (numerics: predictions against predictions, never measurements)
    --probe                          the feasibility screen (sec. 3); needs tau_ABS(all) (sec. 6)
    --probe --confirm                the stacking guard's re-check at C1's start-0 fit (only if the rule would kill)
    --fit {R0,C1} --start N --freeze-commit HASH      one start (resumes from its checkpoint)
    --laplace {R0,C1} --freeze-commit HASH            comparison rule, C1's noise continuation, Laplace, profiles, G1
    --score --freeze-commit HASH [--tau-abs X]        the decisive comparison (sec. 7) under every reading, and sec. 8
    --all --freeze-commit HASH       R0 and C1 in parallel lanes (MAILLARD_WORKERS=4 each), sec. 10's order; the
                                     free-memory check first
Writes results/validation/kinetic_core_m1_supply_combinations_{checks,fit_<R>_s<N>,carried_<R>,side_effects}.json,
kinetic_core_m1_feasibility_probe.{json,md} and kinetic_core_m1_decisive_comparison.{json,md}.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import math
import os
import re
import resource
import signal
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# B60's objective (B47's 52 rows, B51's pH state, B60's quick setting). Its import is grandfathered by the import-time
# switch gate and sets ``ph.b51`` to its shipped value only; this file sets no switch at import.
import generate_kinetic_core_b60_fit as B60  # noqa: E402
import build_submodel_probe as probe  # noqa: E402
import m1_preflight_checks as M10  # noqa: E402  (G1's pot runner; no switch set at import)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import config, fit_rows, fitting, observations, parallel, sulfur  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402

B23, B8 = B60.B23, B60.B8
CELSIUS = 273.15
V = data_paths.VALIDATION_DIR
GENERATED_BY = "scripts/generators/fit_m1_supply_combinations.py"
STEM = "kinetic_core_m1_supply_combinations"
PREREG = V / "kinetic_core_m1_supply_loss_combinations_prereg.md"
CHECKS_OUT = V / f"{STEM}_checks.json"
PROBE_OUT = V / "kinetic_core_m1_feasibility_probe.json"
COMPARISON_OUT = V / "kinetic_core_m1_decisive_comparison.json"
SIDE_OUT = V / f"{STEM}_side_effects.json"
B60_REPORT = V / "kinetic_core_b60_fit_report.json"
B60_START1 = V / "kinetic_core_b60_members" / "b60_s1.json"
PREFLIGHT = V / "kinetic_core_m1_preflight_checks.json"
#: EX-059 (H0)'s scored record: tau_ABS(all) is ``D2.tau_ABS_all_used`` there (the two-set maximum the scorer uses)
TAU_FILE = V / "kinetic_core_thiol_lolo_nested_tau.json"

# ----------------------------------------------------------------------------------------------------------------
# The rivals' coordinates (prereg secs. 2.1, 2.5, 4.2)
# ----------------------------------------------------------------------------------------------------------------
#: B60's Laplace-identified coordinates, free in both rivals, B60's bounds (except the dimer's upper bound)
FREE_B60: Tuple[str, ...] = (
    "k_tdp_fur", "k_nf_mft", "k_nf_mp3p", "k_mgo_mp", "k_fft_decay", "k_fur_decay", "k_pent_caramel", "k_cys_thermal",
    "k_thiolate_loss", "ph_arp_secondary_ammonium_pKa", "k_osone_decay", "k_ttca_deg", "Ea_decay_thiol_sink",
    "k_dimer_mft", "k_arp_dpo", "k_arp_tdp")
#: B60's unidentified free coordinates, held at B60 in both rivals
HELD_B60_UNIDENTIFIED: Tuple[str, ...] = (
    "k_glc_ha", "k_dimer_fft", "ph_acid_yield_per_sink_event", "k_fur_fft", "k_dimer_decay", "k_thiol_decay",
    "Ea_decay_carbonyl_sink")
PRIOR_KEY = "k_mft_decay"
#: C1's route field fitted outside B23's 48-vector (a CysAmadoriRoute field, banded in config._BANDS)
ROUTE_FIELD_KEYS = ("ea_form_kj_mol",)
#: the probe's route fields (sec. 3): K on TtcaFormation, the rest on CysAmadoriRoute
PROBE_ROUTE_FIELDS = ("k_l_per_mmol", "k_form_ref_per_min", "ea_decay_kj_mol")
K_TTCA_CYS_HALF_WIDTH_DEX = 3.0          # sec. 2.2: +-3 decades about the start (stage-1 convention)
EA_FORM_BAND = (40.0, 135.0)             # sec. 2.2: sugar_trunk class band
EA_DECAY_PROBE_BAND = (40.0, 135.0)
K_TTCA_PROBE_BAND = (0.045, 0.075)
K_FORM_PROBE_HALF_DEX = 0.44


@dataclass(frozen=True)
class Rival:
    name: str
    keys: Tuple[str, ...]
    route: bool          # the route switch on (and the route block in the objective)


RIVALS: Dict[str, Rival] = {
    "R0": Rival("R0", FREE_B60 + (PRIOR_KEY,), False),
    "C1": Rival("C1", FREE_B60 + (PRIOR_KEY, "k_ttca_cys") + ROUTE_FIELD_KEYS, True),
}
assert len(RIVALS["R0"].keys) == 17 and len(RIVALS["C1"].keys) == 19
assert set(FREE_B60) | set(HELD_B60_UNIDENTIFIED) == set(B60.FREE_KEYS)

# ----------------------------------------------------------------------------------------------------------------
# Optimiser, convergence, budgets (sec. 4.3; sec. 7 profiles; sec. 4.4 numerics)
# ----------------------------------------------------------------------------------------------------------------
BUDGET = 1200
MAX_SEGMENTS = 2                 # one continuation of a status-0 start
COMPARISON_NATS = 1.92           # Delta chi^2 3.84 / 2
NOISE_CONTINUATION = 300
PROFILE_BUDGET = 300
PROFILE_OFFSETS = {"log10k": (0.25, 0.5, 1.0), "Ea": (1.0, 3.0, 10.0)}
MC_DRAWS = 64
MC_SEED = 20261104
PROBE_DRAWS = 256
PROBE_SEED = 20261104
PROBE_LOCAL_BUDGET = 200
POT_SECONDS = 120                # the house wall-clock limit per pot (build_submodel_probe.run_pot_limited)
ATOL = 1e-14                     # the engine's atol on these lanes (B23.run_systems, build_submodel_probe)
B60_RECORD_COST = 27.13196188533709   # kinetic_core_b60_fit_report.json, start 0 (prereg: 27.132)
V1_REL_TOL = 1e-3                # fit_thiol_lolo.REPRODUCTION_TOL, the house reproduction tolerance
V2_TOL_DEX = 1e-6
V7_DECISION_FRACTION = 0.1
V8_STEPS = {"log10k": (0.025, 0.05), "Ea": (0.25, 0.5)}
V8_COORDINATES = ("k_mft_decay", "Ea_decay_thiol_sink", "k_dimer_mft", "k_ttca_cys", "ea_form_kj_mol")
V5_LOG10_K_DIMER = 4.0
V5_MFT_CAP_TRIGGER = 0.1         # mmol/L: above it, the cap becomes 10x the panel's maximum
FIT_LANE_WORKERS = 4


@dataclass(frozen=True)
class Numerics:
    """One numerical setting: the pH coupling, the integrator's rtol and (fixed point only) the node settings."""

    name: str
    coupling: str
    rtol: float
    ph_nodes: int
    ph_iterations: int
    quick: bool          # B23.observables' flux-budget grid (21 points quick, 101 careful)


def fit_inline() -> Numerics:
    return Numerics("fit_inline_rtol1e-5", "inline", 1e-5, B23.PH_NODES_QUICK, B23.PH_ITERS_QUICK, True)


def fit_fixed_point() -> Numerics:
    """B60's quick setting (17 nodes, rtol 1e-5): V1's reproduction and the fallback if V8 fails."""
    choice = B60.QUICK_SETTING
    return Numerics("fit_fixed_point_quick", "fixed_point_nodes", float(choice["rtol"]), int(choice["ph_nodes"]),
                    B23.PH_ITERS_QUICK, True)


def score_inline() -> Numerics:
    return Numerics("score_inline_rtol1e-8", "inline", 1e-8, B23.PH_NODES_CAREFUL, B23.PH_ITERS_CAREFUL, False)


def score_fixed_point_81() -> Numerics:
    return Numerics("score_fixed_point_81_nodes_rtol1e-8", "fixed_point_nodes", 1e-8, 81, B23.PH_ITERS_CAREFUL, False)


def shipped_fixed_point() -> Numerics:
    """R0-shipped (a): the shipped coupling (node fixed point) at the scoring rtol, B23's careful node setting."""
    return Numerics("shipped_fixed_point_rtol1e-8", "fixed_point_nodes", 1e-8, B23.PH_NODES_CAREFUL,
                    B23.PH_ITERS_CAREFUL, False)


# ----------------------------------------------------------------------------------------------------------------
# The decisive data (sec. 5), their noise (sec. 6), the floors (sec. 5) and the readings (sec. 6, v3.3)
# ----------------------------------------------------------------------------------------------------------------
D1_TABLE, D1_POT = "thiol_submodel", "schieberle2000_rib_cys_ph5_100C"
D1_SPECIES = ("MFT", "FFT")
D1_TIMES = (30.0, 360.0, 720.0)            # the 60 min point is not used (as in v1)
D1_RATIOS = ((360.0, 30.0), (720.0, 360.0))   # 720/30 is their sum, reported only
#: sec. 6: Hofmann & Schieberle 1998 Table 1 footnote b, triplicates within 10 %, as the sd of a printed value
SIGMA_P_DEX = math.log10(1.10)
#: sec. 6 sensitivity (b): the s042 series noise (M1.0), a misfit property
S042_NOISE_DEX = {"MFT": 0.80, "FFT": 1.16}
#: sec. 5: the lowest printed point per series, ug/L (4.5 / 2.0 ug per 100 mL); alternative the table's "< 0.1 ug"
D1_FLOOR_UG_L = {"MFT": 45.0, "FFT": 20.0}
D1_FLOOR_ALT_UG_L = {"MFT": 1.0, "FFT": 1.0}
RES_TABLE, RES_POT = "route_submodel", "zhang2025b_blank_fig2gh"
RES_SPECIES = ("Cys", "PENT")
#: sec. 5: Fig. 2g - h at 50 min; alternative the shared ELSD curve's zero crossing
RES_FLOOR_MMOL_L = {"Cys": 0.80, "PENT": 1.22}
RES_FLOOR_ALT_MMOL_L = {"Cys": 0.22, "PENT": 0.22}
CAO_POTS = ("cao2017_blank_fig4", "cao2017_blank_fig5")
#: sec. 6: the Hofmann autoclave's come-up, 125 K / 20 min
HOFMANN_RATE_K_PER_MIN = 125.0 / 20.0
START_C = 20.0
NEWTON_TAUS = {"ramp_only_newtonian_tau4": 4.0, "ramp_only_newtonian_tau6.7": 20.0 / 3.0}
READINGS = ("isothermal", "ramp_only_linear", "ramp_plus_hold", "ramp_only_newtonian_tau4",
            "ramp_only_newtonian_tau6.7", "ramp_only_linear_with_cooling")
#: v3.3: the reservoir pot is isothermal under every reading; the stand-in ramp is a sensitivity outside agreement
RESERVOIR_STAND_IN = "reservoir_stand_in_ramp_20_to_120C_16min"

#: sec. 5 / 8: the reservoir record's structural zeros of R0 (no step forms TTCA or CAMA from cysteine and pentose)
STRUCTURAL_ZERO_POTS = {"cao2017_blank_fig5": ("TTCA", "CAMA"),
                        **{f"zhai2019_{t}c": ("TTCA", "CAMA") for t in (70, 90, 110, 130)},
                        "zhang2025b_blank": ("TTCA", "CAMA")}

# ----------------------------------------------------------------------------------------------------------------
# The route block (sec. 4.1)
# ----------------------------------------------------------------------------------------------------------------
ROUTE_OUT_POTS = ("zhang2025b_blank", "zhang2025b_blank_fig2gh",          # the reservoir block and its record
                  "zhai2021_ttca_ph55_140c_fig3a", "zhai2023foodchem_ttca_ph7_140c_fig3a",
                  "zhai2023foodchem_ttcaxyl_ph7_140c_fig3a", "zhai2023jafc_ttca_ph7_140c_fig5a",
                  "zhai2023jafc_ttcacys_ph7_140c_fig5a")                    # the declared Zhai >= 140 C hold-out
#: sec. 4.1, the duplicates check: B60's ``kang_free_cys_yield_ceiling`` reads Kang 2026 SI Fig. S3's 120 C curve,
#: and route_submodel.yml's duplicate_note says this pot's TTCA-alone 120 / 140 C curves ARE Fig. S3 (lab_ids.yml:
#: kang2026_SI copies zhai2023foodchem). B60's row stays; the route block drops the pot. Kang's cysteine-conversion
#: rows (Fig. S4, the free-cysteine control) have no counterpart pot in the route table.
ROUTE_DUPLICATE_POTS = ("zhai2023foodchem_ttca_ph7_120c_fig3a",)
ROUTE_DUPLICATE_EVIDENCE = "Kang 2026 SI Fig. S3"


def route_series(row) -> str:
    """sec. 4.1: one extra-noise series per data set x species (the ledger entry is the data set)."""
    return f"{row.ledger_id}|{row.species}"


def route_model(noise: Optional[Mapping[str, float]] = None) -> fitting.MeasurementModel:
    """Absolute rows (every route row is 'calibrated': no free scale); sigma^2 = the row's own error^2 + s_series^2."""
    return fitting.MeasurementModel(series_of=route_series, sigma_of=lambda r: r.sigma, method_of=lambda r: r.method,
                                    extra_noise=dict(noise or {}))


_ROUTE_CACHE: Dict[str, Any] = {}


def route_block() -> Tuple[Tuple[Any, ...], Dict[str, Dict[str, Any]], Dict[str, dict]]:
    """(rows, jobs by pot, pots) of the route block: route_submodel.yml less the declared pots, every pot runnable."""
    if "rows" not in _ROUTE_CACHE:
        pots, rows = load_table(RES_TABLE)
        out = set(ROUTE_OUT_POTS) | set(ROUTE_DUPLICATE_POTS)
        for p in ROUTE_DUPLICATE_POTS:   # the duplicates check is read off the table, not asserted from memory
            assert ROUTE_DUPLICATE_EVIDENCE in str(pots[p].get("duplicate_note", "")), p
        kept = tuple(r for r in rows if r.pot not in out)
        jobs = {j["pot"]: j for j in probe.plan(pots, kept) if j["pot"] in {r.pot for r in kept}}
        skipped = [p for p, j in jobs.items() if "skip" in j]
        if skipped:
            raise SystemExit(f"route pots the lane cannot run: {skipped}")
        _ROUTE_CACHE.update(rows=kept, jobs=jobs, pots=pots)
    return _ROUTE_CACHE["rows"], _ROUTE_CACHE["jobs"], _ROUTE_CACHE["pots"]


_TABLES: Dict[str, Any] = {}


def load_table(table: str):
    """fit_rows.load, once per process (the tables are read-only at run time)."""
    if table not in _TABLES:
        _TABLES[table] = fit_rows.load(table)
    return _TABLES[table]


def table_job(table: str, pot: str, times: Sequence[float]) -> Dict[str, Any]:
    pots, rows = load_table(table)
    (job,) = [j for j in probe.plan({pot: pots[pot]}, rows) if j["pot"] == pot]
    if "skip" in job:
        raise SystemExit(f"{table}:{pot}: {job['skip']}")
    job["times"] = sorted({0.0, *map(float, times)})
    return job


# ----------------------------------------------------------------------------------------------------------------
# Pure rules (unit-tested): floors, the encounter ceiling, the prior, the comparison rule, the decision rule
# ----------------------------------------------------------------------------------------------------------------

def floored(value: float, floor: float) -> Tuple[float, bool]:
    """sec. 5 (re-review B1a): a predicted level below the assay's floor is scored AT the floor."""
    v = float(value)
    if not math.isfinite(v):
        return v, False          # a failed run is not a floored prediction
    return (float(floor), True) if v < floor else (v, False)


def objective_temperatures_c(route: bool) -> List[float]:
    temps = [float(s["t_c"]) for s in B23.SYSTEMS.values()]
    if route:
        temps += [float(j["t_c"]) for j in route_block()[1].values()]
    return temps


def dimer_ceiling(ceiling: ps.DimerEncounterCeiling, route: bool) -> float:
    """sec. 2.5: the minimum over the objective's pot temperatures of the encounter ceiling at 145 C."""
    return ceiling.log10_ceiling(objective_temperatures_c(route))


def mft_cap_from_panel_max(panel_max_mmol_l: float) -> float:
    """sec. 2.5: [MFT] capped at 1 mmol/L unless the panel's maximum exceeds 0.1 mmol/L; then 10x that maximum."""
    return 1.0 if panel_max_mmol_l <= V5_MFT_CAP_TRIGGER else 10.0 * float(panel_max_mmol_l)


def block_se(d: Sequence[float]) -> float:
    """sec. 7: SE_b = sqrt(n / (n - 1) x sum (d_u - mean d)^2); with n = 2 it is |d1 - d2|."""
    d = np.asarray(d, dtype=float)
    n = len(d)
    if n < 2:
        return float("inf")
    return float(math.sqrt(n / (n - 1.0) * float(np.sum((d - d.mean()) ** 2))))


def unit_side(signs: Sequence[float]) -> str:
    """v3.3 tie band: a unit favours a rival only if its d_u has the same sign under the linearised and plug-in
    densities and under the floor and tau x0.5 / x2 sensitivities; otherwise it is a tie."""
    s = [float(x) for x in signs]
    if s and all(x > 0.0 for x in s):
        return "C1"
    if s and all(x < 0.0 for x in s):
        return "R0"
    return "tie"


C1_FAVOURED, R0_FAVOURED, NOT_SEPARATED = "C1 favoured", "R0 favoured", "not separated / conflict"


def block_outcome(reservoir_sides: Sequence[str], d1_sides: Sequence[str]) -> str:
    """sec. 7 rules 2 - 4 with v3.3's tie band. With two units per block the one-SE rule is a sign test: C1 favoured
    iff both reservoir units favour C1 AND not both D1 units favour R0 (a D1 tie counts as not worse); R0 favoured iff
    both reservoir units favour R0; otherwise not separated (a reservoir tie gives not separated)."""
    if reservoir_sides and all(s == "C1" for s in reservoir_sides):
        d1_worse = bool(d1_sides) and all(s == "R0" for s in d1_sides)
        return NOT_SEPARATED if d1_worse else C1_FAVOURED
    if reservoir_sides and all(s == "R0" for s in reservoir_sides):
        return R0_FAVOURED
    return NOT_SEPARATED


def verdict_across_readings(outcomes: Mapping[str, str]) -> Dict[str, Any]:
    """sec. 7 rule 5: the verdict counts only if it is the same under every reading; otherwise not separated, naming
    the readings whose outcome differs from the isothermal reading's."""
    values = set(outcomes.values())
    if len(values) == 1:
        return {"verdict": next(iter(values)), "agree": True, "disagreeing_readings": []}
    iso = outcomes.get("isothermal")
    return {"verdict": NOT_SEPARATED, "agree": False,
            "disagreeing_readings": sorted(r for r, o in outcomes.items() if o != iso)}


def comparison_rule(starts: Sequence[Mapping[str, Any]], predictions: Optional[Mapping[int, np.ndarray]] = None,
                    tolerance: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """sec. 4.3 (re-review M2): compared iff (a) >= 2 starts end with status > 0 and costs within 1.92 nats of the
    lowest, and (b) their decisive predictions agree within max(sigma_decision, Laplace predictive sd) at every
    point (without reading y). (a) and not (b): compared only if the per-start verdicts agree (decided at scoring)."""
    converged = [s for s in starts if int(s["status"]) > 0]
    out: Dict[str, Any] = {"converged_starts": [int(s["start"]) for s in converged]}
    if not converged:
        out.update(a=False, b=None, status="not_compared", best_start=None, qualifying_starts=[])
        return out
    best = min(converged, key=lambda s: float(s["cost"]))
    lowest = float(best["cost"])
    qualifying = [s for s in converged if float(s["cost"]) - lowest <= COMPARISON_NATS]
    a = len(qualifying) >= 2
    out.update(best_start=int(best["start"]), lowest_cost=lowest, qualifying_starts=[int(s["start"]) for s in qualifying],
               cost_gaps_nats={int(s["start"]): float(s["cost"]) - lowest for s in converged}, a=a)
    if not a:
        out.update(b=None, status="not_compared")
        return out
    if predictions is None or tolerance is None:
        out.update(b=None, status="pending_predictions")
        return out
    ref = np.asarray(predictions[int(best["start"])], dtype=float)
    tol = np.asarray(tolerance, dtype=float)
    gaps = {int(s["start"]): float(np.max(np.abs(np.asarray(predictions[int(s["start"])]) - ref) / tol))
            for s in qualifying}
    b = all(g <= 1.0 for g in gaps.values())
    out.update(b=b, max_gap_over_tolerance=gaps, status="compared" if b else "compared_if_verdicts_agree")
    return out


def prior_log10_rate(x_full: np.ndarray, prior: ps.MftDecayPrior) -> float:
    p = B23.build_parameters(np.asarray(x_full, dtype=float))[PRIOR_KEY]
    return prior.log10_rate(p.k_ref, p.ea_kj_mol, p.t_ref_k)


def lambda_log10_gradient() -> Dict[str, float]:
    """d log10 lambda_mft(100 C) / d(coordinate): 1 on log10 k_mft_decay(145 C), -(1/T100 - 1/T145)/(R ln10) per
    kJ/mol on Ea_decay_thiol_sink (the engine's Arrhenius form, t_ref 145 C)."""
    t, t_ref = 100.0 + CELSIUS, 145.0 + CELSIUS
    return {PRIOR_KEY: 1.0, "Ea_decay_thiol_sink": -(1.0 / t - 1.0 / t_ref) / (ps.R_KJ * math.log(10.0))}


# ----------------------------------------------------------------------------------------------------------------
# Points and configurations
# ----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FitSide:
    ceiling: ps.DimerEncounterCeiling
    prior: ps.MftDecayPrior

    @classmethod
    def off(cls) -> "FitSide":
        return cls(ps.DimerEncounterCeiling(), ps.MftDecayPrior())

    @classmethod
    def on(cls, mft_cap: float) -> "FitSide":
        return cls(ps.DimerEncounterCeiling(enabled=True, mft_cap_mmol_l=float(mft_cap)), ps.MftDecayPrior(enabled=True))


@dataclass(frozen=True)
class Point:
    """Everything a configuration of one rival at one vector needs: B23's 48-vector, the route switch and its fields."""

    x_full: Tuple[float, ...]
    route: bool
    route_fields: Tuple[Tuple[str, float], ...] = ()

    def x(self) -> np.ndarray:
        return np.asarray(self.x_full, dtype=float)


def shipped_value(name: str) -> Any:
    return config.KEY_BY_NAME[name].default


def config_of(point: Point, numerics: Numerics, fit_side: FitSide) -> config.Config:
    """The Config of one evaluation, built from the SHIPPED declarations (not this process's state): the coupling, the
    route on (C1) with its fields, and the fit-side keys. B59a, B61 - B63, B7a and every other switch stay shipped."""
    changes: Dict[str, Any] = {"ph.coupling": numerics.coupling,
                               "fit.dimer_encounter_ceiling": fit_side.ceiling, "fit.mft_decay_prior": fit_side.prior}
    if point.route:
        fields = dict(point.route_fields)
        k = fields.pop("k_l_per_mmol", None)
        changes["sulfur.cys_amadori_route"] = dataclasses.replace(shipped_value("sulfur.cys_amadori_route"),
                                                                  enabled=True, **fields)
        if k is not None:
            changes["sulfur.ttca_formation"] = dataclasses.replace(shipped_value("sulfur.ttca_formation"),
                                                                   k_l_per_mmol=float(k))
    elif point.route_fields:
        raise ValueError("route fields on a point with the route off")
    return config.Config.shipped().replace(changes)


def b60_vector() -> np.ndarray:
    lower, upper = B60.full_bounds()
    return np.clip(frozen_vector(json.loads(B60_REPORT.read_text(encoding="utf-8"))), lower, upper)


def b60_point(route: bool = False) -> Point:
    fields = (("ea_form_kj_mol", float(shipped_value("sulfur.cys_amadori_route").ea_form_kj_mol)),) if route else ()
    return Point(tuple(map(float, b60_vector())), route, fields)


def kind_of(key: str) -> str:
    if key.startswith(("Ea_", "ea_")):
        return "Ea"
    if key == "ph_arp_secondary_ammonium_pKa":
        return "pka"
    if key == "ph_acid_yield_per_sink_event":
        return "acid_yield"
    return "log10k"


def all_index(key: str) -> int:
    return B60.ALL_KEYS.index(key)


# ----------------------------------------------------------------------------------------------------------------
# Running pots under a wall-clock alarm, failures counted where the results arrive
# ----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class PotFailure:
    name: str
    reason: str
    timeout: bool


class PotFailed(RuntimeError):
    pass


class _Timeout(RuntimeError):
    pass


def with_alarm(seconds: int, fn: Callable[[], Any]) -> Any:
    """``fn()`` under a Unix wall-clock alarm (works in forked workers, as build_submodel_probe.run_pot_limited)."""
    def _alarm(signum, frame):
        raise _Timeout(f"exceeded {seconds} s")

    previous = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def _integrate_system(item) -> Any:
    """One B60 system (B23.run_systems' integration, every argument the same) under the alarm."""
    name, spec, parameters, drift, numerics = item
    minutes = float(spec["minutes"])
    points = int(spec.get("grid_points", 0))
    grid = np.linspace(0.0, minutes, points) if points >= 3 else np.array([0.0, minutes])
    try:
        return with_alarm(POT_SECONDS, lambda: B23.integrate_sulfur(
            parameters, float(spec["t_c"]) + CELSIUS, spec["initial"], grid, ph=float(spec["ph"]),
            buffer_spec=spec.get("buffer"), ph_drift=drift, ph_nodes=numerics.ph_nodes,
            ph_iterations=numerics.ph_iterations, rtol=numerics.rtol, atol=ATOL))
    except _Timeout as exc:
        return PotFailure(f"b60:{name}", str(exc), True)
    except Exception as exc:  # noqa: BLE001 - counted; the 52-row block then scores as a failed evaluation
        return PotFailure(f"b60:{name}", f"{type(exc).__name__}: {exc}", False)


def _run_pot_safe(job: Dict[str, Any]) -> Dict[str, Any]:
    """A fit-table pot (build_submodel_probe's path) under the alarm; a failure is returned, not raised."""
    try:
        return with_alarm(POT_SECONDS, lambda: probe.run_pot(job))
    except _Timeout as exc:
        return {"failed": str(exc), "timeout": True}
    except Exception as exc:  # noqa: BLE001
        return {"failed": f"{type(exc).__name__}: {exc}", "timeout": False}


@contextlib.contextmanager
def b60_runner(numerics: Numerics, failures: List[Dict[str, Any]]):
    """B23's ``run_systems`` / flux budget routed through the alarm and the given numerics, restored on exit. A failed
    system raises inside B23's objective (which then scores the 52 rows as a failed evaluation, as in B60) after the
    failure is recorded HERE, in the process the results came back to."""
    saved = (B23.run_systems, B23.sulfur_flux_budget, B23.branch_shares)
    original_budget, original_shares = saved[1], saved[2]

    def run_systems(parameters, quick, drift=None):
        names = list(B23.SYSTEMS)
        items = [(n, B23.SYSTEMS[n], parameters, drift, numerics) for n in names]
        runs = parallel.map_ordered(_integrate_system, items, config=config.Config.from_globals())
        bad = [r for r in runs if isinstance(r, PotFailure)]
        if bad:
            failures.extend({"pot": f.name, "reason": f.reason, "timeout": f.timeout} for f in bad)
            raise PotFailed(", ".join(f.name for f in bad))
        return {n: {"run": r, "spec": B23.SYSTEMS[n]} for n, r in zip(names, runs)}

    def flux_budget(*args, **kwargs):
        try:
            return with_alarm(POT_SECONDS, lambda: original_budget(*args, **kwargs))
        except _Timeout as exc:
            return PotFailure("b60:flux_budget", str(exc), True)
        except Exception as exc:  # noqa: BLE001
            return PotFailure("b60:flux_budget", f"{type(exc).__name__}: {exc}", False)

    def branch_shares(fluxes):
        if isinstance(fluxes, PotFailure):
            failures.append({"pot": fluxes.name, "reason": fluxes.reason, "timeout": fluxes.timeout})
            raise PotFailed(fluxes.name)
        return original_shares(fluxes)

    B23.run_systems, B23.sulfur_flux_budget, B23.branch_shares = run_systems, flux_budget, branch_shares
    try:
        yield
    finally:
        B23.run_systems, B23.sulfur_flux_budget, B23.branch_shares = saved


def table_tasks(point: Point, jobs: Sequence[Dict[str, Any]], numerics: Numerics) -> List[Dict[str, Any]]:
    x = point.x()
    parameters, drift = B23.build_parameters(x), B23.unpack(x)[3]
    return [dict(j, parameters=parameters, drift=drift, rtol=numerics.rtol, atol=ATOL, ph_nodes=numerics.ph_nodes,
                 ph_iterations=numerics.ph_iterations) for j in jobs]


def run_table_pots(point: Point, jobs: Sequence[Dict[str, Any]], numerics: Numerics, cfg: config.Config,
                   failures: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    tasks = table_tasks(point, jobs, numerics)
    results = parallel.map_ordered(_run_pot_safe, tasks, config=cfg)
    for t, r in zip(tasks, results):
        if "failed" in r:
            failures.append({"pot": f"{t.get('pot')}", "reason": r["failed"], "timeout": bool(r.get("timeout"))})
    return tasks, results


# ----------------------------------------------------------------------------------------------------------------
# The objective (sec. 4.1): B60's 52 rows [+ the route block] + the prior
# ----------------------------------------------------------------------------------------------------------------

def b60_row_ids() -> List[str]:
    return [r["id"] for r in B23.ACTIVE_FIT_ROWS]


def route_raw(point: Point, numerics: Numerics, cfg: config.Config, failures: List[Dict[str, Any]]) -> np.ndarray:
    rows, jobs, _ = route_block()
    pots = sorted({r.pot for r in rows})
    tasks, results = run_table_pots(point, [jobs[p] for p in pots], numerics, cfg, failures)
    res, task = dict(zip(pots, results)), dict(zip(pots, tasks))
    out = []
    for r in rows:
        got = res[r.pot]
        if "failed" in got:
            out.append(float("nan"))
            continue
        predicted = probe.predict(r, task[r.pot], got)
        out.append(fit_rows.residual(r, predicted) if predicted > 0 else float("nan"))
    return np.asarray(out, dtype=float)


@dataclass(frozen=True)
class EvalSpec:
    numerics: Numerics
    fit_side: FitSide
    route_block: bool
    noise: Tuple[Tuple[str, float], ...] = ()


def evaluate(spec: EvalSpec, point: Point) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
    """The objective's residual vector at one point and the pot failures met computing it (in this process)."""
    cfg = config_of(point, spec.numerics, spec.fit_side)
    failures: List[Dict[str, Any]] = []
    with config.applied(cfg), b60_runner(spec.numerics, failures):
        parts = [np.asarray(B60.residual_vector(point.x(), spec.numerics.quick), dtype=float)]
        if spec.route_block:
            rows = route_block()[0]
            parts.append(route_model(dict(spec.noise)).whiten(route_raw(point, spec.numerics, cfg, failures), rows))
        if spec.fit_side.prior.enabled:
            parts.append(np.array([spec.fit_side.prior.residual(prior_log10_rate(point.x(), spec.fit_side.prior))]))
    return np.concatenate(parts), failures


#: the pot failures of this process's evaluations (and the ones workers sent back with their results)
FAILURES: List[Dict[str, Any]] = []


class Objective:
    """One rival's objective over its free coordinates (some may be FIXED: a profile point)."""

    def __init__(self, rival: Rival, numerics: Numerics, fit_side: FitSide, noise: Optional[Mapping[str, float]] = None,
                 base: Optional[Point] = None, fixed: Optional[Mapping[str, float]] = None):
        self.rival, self.numerics, self.fit_side = rival, numerics, fit_side
        self.fixed = dict(fixed or {})
        self.keys = tuple(k for k in rival.keys if k not in self.fixed)
        self.kinds = tuple(kind_of(k) for k in self.keys)
        self.base = base or b60_point(rival.route)
        self.spec = EvalSpec(numerics, fit_side, rival.route, tuple(sorted((noise or {}).items())))
        lo, hi = self.bounds(rival.keys)
        self.lower = np.array([lo[rival.keys.index(k)] for k in self.keys])
        self.upper = np.array([hi[rival.keys.index(k)] for k in self.keys])
        n_route = len(route_block()[0]) if rival.route else 0
        self.row_labels = tuple([f"b60:{i}" for i in b60_row_ids()]
                                + ([f"route:{r.pot}|{r.species}|{r.t_min:g}" for r in route_block()[0]]
                                   if rival.route else [])
                                + (["prior:lambda_mft_100C"] if fit_side.prior.enabled else []))
        self.n_b60, self.n_route = len(B23.ACTIVE_FIT_ROWS), n_route

    def bounds(self, keys: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
        lo_full, hi_full = B60.full_bounds()
        b60 = b60_vector()
        lo, hi = [], []
        for k in keys:
            if k in ROUTE_FIELD_KEYS:
                a, b = EA_FORM_BAND
            elif k == "k_ttca_cys":
                s = float(b60[all_index(k)])
                a, b = s - K_TTCA_CYS_HALF_WIDTH_DEX, s + K_TTCA_CYS_HALF_WIDTH_DEX
            else:
                a, b = float(lo_full[all_index(k)]), float(hi_full[all_index(k)])
                if k in ("k_dimer_mft", "k_dimer_fft") and self.fit_side.ceiling.enabled:
                    b = dimer_ceiling(self.fit_side.ceiling, self.rival.route)
            lo.append(a)
            hi.append(b)
        return np.array(lo), np.array(hi)

    def point(self, theta: Sequence[float]) -> Point:
        values = dict(self.fixed)
        values.update(zip(self.keys, map(float, theta)))
        x = list(self.base.x_full)
        fields = dict(self.base.route_fields)
        for k, v in values.items():
            if k in ROUTE_FIELD_KEYS or k in PROBE_ROUTE_FIELDS:
                fields[k] = v
            else:
                x[all_index(k)] = v
        return Point(tuple(x), self.rival.route, tuple(sorted(fields.items())) if self.rival.route else ())

    def theta_of(self, point: Point) -> np.ndarray:
        fields = dict(point.route_fields)
        return np.array([fields[k] if k in ROUTE_FIELD_KEYS else point.x_full[all_index(k)] for k in self.keys])

    def residuals(self, theta, rows=None) -> np.ndarray:
        r, failed = evaluate(self.spec, self.point(theta))
        FAILURES.extend(failed)
        return r

    def residuals_many(self, points, rows=None) -> List[np.ndarray]:
        spec, pts = self.spec, [self.point(t) for t in points]
        results = parallel.map_ordered(_evaluate_job, [(spec, p) for p in pts])
        out = []
        for r, failed in results:
            FAILURES.extend(failed)
            out.append(r)
        return out

    def problem(self) -> fitting.Problem:
        return fitting.Problem(keys=self.keys, kinds=self.kinds, lower=self.lower, upper=self.upper,
                               rows=self.row_labels, residuals=self.residuals, evaluate_many=self.residuals_many)

    def blocks(self, r: np.ndarray) -> Dict[str, float]:
        r = np.asarray(r, dtype=float)
        a, b = self.n_b60, self.n_b60 + self.n_route
        out = {"b60_52_rows": 0.5 * float(r[:a] @ r[:a]), "route_block": 0.5 * float(r[a:b] @ r[a:b]),
               "prior": 0.5 * float(r[b:] @ r[b:])}
        total = sum(out.values())
        return {**out, "total": total, "share": {k: (v / total if total else None) for k, v in out.items()}}


def _evaluate_job(item):
    spec, point = item
    return evaluate(spec, point)


# ----------------------------------------------------------------------------------------------------------------
# Decisive predictions (sec. 5 / 6): log10 levels at the decisive points, per reading
# ----------------------------------------------------------------------------------------------------------------

def _ramp_end(t_target_c: float) -> float:
    return (t_target_c - START_C) / HOFMANN_RATE_K_PER_MIN


def _newton_profile(t_target_c: float, tau: float, minutes: float) -> Tuple[Tuple[float, float], ...]:
    ts = sorted({*np.linspace(0.0, 10.0 * tau, 81), float(minutes)})
    return tuple((float(t), float(t_target_c - (t_target_c - START_C) * math.exp(-t / tau))) for t in ts)


def d1_runs(reading: str) -> List[Dict[str, Any]]:
    """sec. 6: the Schieberle pot's runs under one reading: [(profile, times integrated, {nominal t: integrated t})]."""
    t_c, end = 100.0, max(D1_TIMES)
    ramp = _ramp_end(t_c)
    if reading == "isothermal":
        return [{"profile": None, "observe": {t: t for t in D1_TIMES}}]
    if reading == "ramp_only_linear":
        return [{"profile": ((0.0, START_C), (ramp, t_c), (end, t_c)), "observe": {t: t for t in D1_TIMES}}]
    if reading == "ramp_plus_hold":
        return [{"profile": ((0.0, START_C), (ramp, t_c), (ramp + end, t_c)),
                 "observe": {t: ramp + t for t in D1_TIMES}}]
    if reading in NEWTON_TAUS:
        return [{"profile": _newton_profile(t_c, NEWTON_TAUS[reading], end), "observe": {t: t for t in D1_TIMES}}]
    if reading == "ramp_only_linear_with_cooling":
        # reading 2 plus a cool-down of the ramp's duration after each sampling time: one run per sampling time
        return [{"profile": ((0.0, START_C), (ramp, t_c), (t, t_c), (t + ramp, START_C)), "observe": {t: t + ramp}}
                for t in D1_TIMES]
    raise ValueError(reading)


def reservoir_runs(stand_in: bool = False) -> List[Dict[str, Any]]:
    times = sorted({r.t_min for r in decisive_rows()["reservoir"]})
    if not stand_in:
        return [{"profile": None, "observe": {t: t for t in times}}]
    ramp, end = _ramp_end(120.0), max(times)
    return [{"profile": ((0.0, START_C), (ramp, 120.0), (end, 120.0)), "observe": {t: t for t in times}}]


_DECISIVE: Dict[str, Any] = {}


def decisive_rows() -> Dict[str, Any]:
    if not _DECISIVE:
        _, rows = load_table(D1_TABLE)
        d1 = {(r.species, float(r.t_min)): r for r in rows if r.pot == D1_POT and r.species in D1_SPECIES}
        _, rrows = load_table(RES_TABLE)
        res = tuple(r for r in rrows if r.pot == RES_POT and r.species in RES_SPECIES)
        assert sorted(Counter_(r.species for r in res).values()) == [8, 8], "the reservoir block: 8 points per species"
        _DECISIVE.update(d1=d1, reservoir=res)
    return _DECISIVE


def Counter_(it):
    out: Dict[Any, int] = defaultdict(int)
    for x in it:
        out[x] += 1
    return out


def point_labels() -> List[str]:
    """The decisive LEVEL points (log10, unfloored): D1 MFT / FFT at 30, 360, 720 min; Zhang Cys / PENT at 8 times."""
    labels = [f"D1|{s}|{t:g}" for s in D1_SPECIES for t in D1_TIMES]
    labels += [f"RES|{r.species}|{r.t_min:g}" for r in decisive_rows()["reservoir"]]
    return labels


def decisive_jobs(reading: str, stand_in: bool = False, label_ph_held: bool = False
                  ) -> List[Tuple[str, Dict[str, Any], Dict[str, Any]]]:
    """(block, job, run) per integration of one reading."""
    out = []
    for run in d1_runs(reading):
        job = table_job(D1_TABLE, D1_POT, run["observe"].values())
        if run["profile"] is not None:
            job["temperature_profile_c"] = run["profile"]
        out.append(("D1", job, run))
    for run in reservoir_runs(stand_in):
        job = table_job(RES_TABLE, RES_POT, run["observe"].values())
        if run["profile"] is not None:
            job["temperature_profile_c"] = run["profile"]
        if label_ph_held:   # sec. 7 sensitivity: the label pH held fixed instead of the B51 balance's drift
            job["buffer"] = "pH-stat (M1 sensitivity: the label pH held) phosphate 0.20 M"
        out.append(("RES", job, run))
    return out


def decisive_levels(point: Point, reading: str, numerics: Numerics, fit_side: Optional[FitSide] = None,
                    stand_in: bool = False, label_ph_held: bool = False) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
    """log10 of the predicted level at every decisive point (unfloored; mmol/L for the reservoir, ug/L for D1)."""
    cfg = config_of(point, numerics, fit_side or FitSide.off())
    failures: List[Dict[str, Any]] = []
    plan_ = decisive_jobs(reading, stand_in, label_ph_held)
    with config.applied(cfg):
        tasks, results = run_table_pots(point, [j for _, j, _ in plan_], numerics, cfg, failures)
    level: Dict[str, float] = {}
    for (block, _, run), task, res in zip(plan_, tasks, results):
        if "failed" in res:
            continue
        times = task["times"]
        for nominal, t in run["observe"].items():
            state = res["states"][times.index(float(t))]
            species = D1_SPECIES if block == "D1" else RES_SPECIES
            for s in species:
                unit = "ug_L" if block == "D1" else "mmol_L"
                v = observations.amount_in(unit, s, state[s])
                level[f"{block}|{s}|{nominal:g}"] = math.log10(v) if v > 0 else -30.0
    out = np.array([level.get(lbl, float("nan")) for lbl in point_labels()])
    return out, failures


def decisive_many(points: Sequence[Point], reading: str, numerics: Numerics, **kw) -> List[np.ndarray]:
    """Many points in parallel (the point is the grain); failures recorded in FAILURES."""
    results = parallel.map_ordered(_decisive_job, [(p, reading, numerics, kw) for p in points])
    out = []
    for levels, failed in results:
        FAILURES.extend(failed)
        out.append(levels)
    return out


def _decisive_job(item):
    p, reading, numerics, kw = item
    return decisive_levels(p, reading, numerics, **kw)


# ----------------------------------------------------------------------------------------------------------------
# The measurement model of the decisive units (sec. 6) and the densities (sec. 7)
# ----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Scoring:
    """One scoring setting: the D1 sigma scale, the floors and tau."""

    tau: float
    d1_sigma: str = "decision"       # "decision" | "read_error" | "s042"
    floors: str = "primary"          # "primary" | "alternative"


def d1_point_variance(species: str, t: float, scale: str) -> float:
    read = float(decisive_rows()["d1"][(species, t)].sigma)
    if scale == "decision":
        return SIGMA_P_DEX ** 2 + read ** 2
    if scale == "read_error":
        return read ** 2
    if scale == "s042":
        return S042_NOISE_DEX[species] ** 2 + read ** 2
    raise ValueError(scale)


def units(levels: np.ndarray, jac: Optional[np.ndarray], sc: Scoring) -> Dict[str, Dict[str, Any]]:
    """The four decisive units from the level points: y (data), f (floored prediction), J (zero rows where floored),
    Sigma (sec. 6). D1 per species: the bivariate log ratio 360/30, 720/360, sharing the 360 min point; the reservoir
    per species: 8 log levels with the integrated level offset N(0, tau^2)."""
    labels = point_labels()
    idx = {lbl: i for i, lbl in enumerate(labels)}
    rows = decisive_rows()
    out: Dict[str, Dict[str, Any]] = {}
    for s in D1_SPECIES:
        floor = (D1_FLOOR_UG_L if sc.floors == "primary" else D1_FLOOR_ALT_UG_L)[s]
        lv, jr, n_floored = {}, {}, 0
        for t in D1_TIMES:
            i = idx[f"D1|{s}|{t:g}"]
            v, fl = floored(10.0 ** levels[i] if np.isfinite(levels[i]) else float("nan"), floor)
            lv[t] = math.log10(v) if np.isfinite(v) else float("nan")
            n_floored += int(fl)
            jr[t] = np.zeros(jac.shape[1]) if (jac is not None and fl) else (jac[i] if jac is not None else None)
        meas = {t: math.log10(rows["d1"][(s, t)].value) for t in D1_TIMES}
        y = np.array([meas[a] - meas[b] for a, b in D1_RATIOS])
        f = np.array([lv[a] - lv[b] for a, b in D1_RATIOS])
        v = {t: d1_point_variance(s, t, sc.d1_sigma) for t in D1_TIMES}
        cov = np.array([[v[360.0] + v[30.0], -v[360.0]], [-v[360.0], v[720.0] + v[360.0]]])
        j = None if jac is None else np.vstack([jr[a] - jr[b] for a, b in D1_RATIOS])
        out[f"D1_{s}"] = {"block": "D1", "y": y, "f": f, "J": j, "Sigma": cov, "n_floored": n_floored,
                          "labels": [f"{s} {int(a)}/{int(b)}" for a, b in D1_RATIOS]}
    for s in RES_SPECIES:
        floor = (RES_FLOOR_MMOL_L if sc.floors == "primary" else RES_FLOOR_ALT_MMOL_L)[s]
        rs = [r for r in rows["reservoir"] if r.species == s]
        f, jr, n_floored = [], [], 0
        for r in rs:
            i = idx[f"RES|{s}|{r.t_min:g}"]
            v, fl = floored(10.0 ** levels[i] if np.isfinite(levels[i]) else float("nan"), floor)
            f.append(math.log10(v) if np.isfinite(v) else float("nan"))
            n_floored += int(fl)
            if jac is not None:
                jr.append(np.zeros(jac.shape[1]) if fl else jac[i])
        y = np.array([math.log10(r.value) for r in rs])
        sig = np.array([float(r.sigma) for r in rs])     # the row's SD (dex, where printed) and read error, quadrature
        cov = np.diag(sig ** 2) + sc.tau ** 2 * np.ones((len(rs), len(rs)))
        out[f"RES_{s}"] = {"block": "reservoir", "y": y, "f": np.array(f), "J": None if jac is None else np.vstack(jr),
                           "Sigma": cov, "n_floored": n_floored, "labels": [f"{s} {r.t_min:g} min" for r in rs]}
    return out


def log_normal(y: np.ndarray, f: np.ndarray, cov: np.ndarray) -> float:
    r = np.asarray(y, dtype=float) - np.asarray(f, dtype=float)
    if not np.all(np.isfinite(r)):
        return float("-inf")
    sign, logdet = np.linalg.slogdet(cov)
    return float(-0.5 * (r @ np.linalg.solve(cov, r) + logdet + len(r) * math.log(2.0 * math.pi)))


def unit_lppd(u: Mapping[str, Any], sigma_theta: Optional[np.ndarray]) -> Dict[str, float]:
    plug = log_normal(u["y"], u["f"], u["Sigma"])
    if sigma_theta is None or u["J"] is None:
        return {"plug_in": plug, "linearised": plug}
    extra = u["J"] @ sigma_theta @ u["J"].T
    return {"plug_in": plug, "linearised": log_normal(u["y"], u["f"], u["Sigma"] + extra)}


def mc_lppd(y: np.ndarray, fs: Sequence[np.ndarray], cov: np.ndarray) -> Dict[str, float]:
    """log mean_s N(y; f_s, Sigma) and its MC SE by the delta method (log units)."""
    lp = np.array([log_normal(y, f, cov) for f in fs])
    if not np.any(np.isfinite(lp)):
        return {"mc": float("-inf"), "mc_se": float("nan")}
    m = float(np.max(lp[np.isfinite(lp)]))
    w = np.exp(np.where(np.isfinite(lp), lp - m, -np.inf))
    mean = float(np.mean(w))
    se = float(np.std(w, ddof=1) / (math.sqrt(len(w)) * mean)) if len(w) > 1 and mean > 0 else float("nan")
    return {"mc": m + math.log(mean), "mc_se": se}


# ----------------------------------------------------------------------------------------------------------------
# Profiles of bound-active coordinates (sec. 7) and the identified covariance
# ----------------------------------------------------------------------------------------------------------------

def profile_density(bound: float, side: str, values: Sequence[float], costs: Sequence[float], c_min: float,
                    far: float) -> Dict[str, Any]:
    """The one-sided profile density exp(-(c(v) - c_min)), linear in log density between the profiled points (the
    bound itself at c_min), truncated at the bound; beyond the last point the last segment's slope is extended while
    it falls (declared), to the coordinate's other bound ``far``. Returns grid, density (normalised), mean, variance."""
    sgn = -1.0 if side == "upper" else 1.0
    pts = sorted(zip([0.0] + [abs(v - bound) for v in values], [0.0] + [-(c - c_min) for c in costs]))
    d, ld = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    span = abs(far - bound)
    grid = np.linspace(0.0, span, 4001)
    lg = np.interp(grid, d, ld)
    slope = (ld[-1] - ld[-2]) / (d[-1] - d[-2]) if len(d) >= 2 and d[-1] > d[-2] else 0.0
    beyond = grid > d[-1]
    flagged = False
    if np.any(beyond):
        if slope < 0:
            lg[beyond] = ld[-1] + slope * (grid[beyond] - d[-1])
        else:
            lg[beyond] = -np.inf
            flagged = True
    w = np.exp(lg - np.max(lg))
    w /= np.trapezoid(w, grid)
    coord = bound + sgn * grid
    mean = float(np.trapezoid(w * coord, grid))
    var = float(np.trapezoid(w * (coord - mean) ** 2, grid))
    return {"grid_coordinate": coord.tolist(), "density": w.tolist(), "mean": mean, "variance": var,
            "truncated_beyond_last_point": flagged, "last_slope_log_density_per_unit": float(slope)}


def sample_profile(prof: Mapping[str, Any], rng: np.random.Generator, n: int) -> np.ndarray:
    x, w = np.asarray(prof["grid_coordinate"]), np.asarray(prof["density"])
    order = np.argsort(x)
    x, w = x[order], w[order]
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(x))])
    cdf /= cdf[-1]
    return np.interp(rng.uniform(size=n), cdf, x)


def reflect(x: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    x = np.array(x, dtype=float)
    for _ in range(8):
        x = np.where(x < lo, 2 * lo - x, x)
        x = np.where(x > hi, 2 * hi - x, x)
    return np.clip(x, lo, hi)


# ----------------------------------------------------------------------------------------------------------------
# Records, freeze checks, refusals
# ----------------------------------------------------------------------------------------------------------------

def _write(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(provenance.dumps_compact(_jsonable(payload), default=float) + "\n", encoding="utf-8")


def _read(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}



def failure_accepted(rec: Mapping[str, Any]) -> Tuple[bool, str]:
    """Amendment A1 (2026-11-04, before any decisive number): a fit record with failed pot runs is ACCEPTED only if
    (a) every failure happened inside the optimiser's segments, so the re-evaluation at the optimum was clean;
    (b) the cost recomputed from the optimum's residuals equals the record's cost (relative 1e-9);
    (c) a same-settings re-run reproduced it exactly (kept as ``<stem>_failed_run1.json``: the same cost to 1e-12
        relative, the same evaluations, the same failure count and the same failed pots), i.e. the failure is a
        deterministic exploratory point, not machine state. Otherwise the record must be re-run (sec. 4.3)."""
    n = int(rec.get("failed_pot_runs", 0))
    if n == 0:
        return True, "no failed pot runs"
    seg = sum(int(h.get("failed_pot_runs_segment", 0)) for h in rec.get("history", []))
    if seg != n:
        return False, f"{n - seg} failure(s) at the optimum's own evaluation"
    res = np.array(list((rec.get("residual_by_row") or {}).values()), dtype=float)
    cost = float(rec.get("cost", float("nan")))
    if not res.size or not np.isfinite(cost) or abs(0.5 * float(res @ res) - cost) > 1e-9 * abs(cost):
        return False, "the optimum's residuals do not reproduce its cost"
    twin = data_paths.VALIDATION_DIR / f"{rec.get('artifact')}_failed_run1.json"
    if not twin.exists():
        return False, "no same-settings re-run to show the failure is deterministic"
    t = json.loads(twin.read_text(encoding="utf-8"))
    pots = lambda r: sorted(f.get("pot") for f in (r.get("failures_sample") or []))  # noqa: E731
    same = (abs(float(t.get("cost", float("nan"))) - cost) <= 1e-12 * abs(cost)
            and int(t.get("evaluations", -1)) == int(rec.get("evaluations", -2))
            and int(t.get("failed_pot_runs", -1)) == n and pots(t) == pots(rec))
    if not same:
        return False, "the same-settings re-run did not reproduce the record"
    return True, f"{n} deterministic exploratory failure(s) {pots(rec)}: optimum clean, reproduced (amendment A1)"

def fit_path(rival: str, start: int) -> Path:
    return V / f"{STEM}_fit_{rival}_s{start}.json"


def carried_path(rival: str) -> Path:
    return V / f"{STEM}_carried_{rival}.json"


def _provenance(name: str, extra: Sequence[Path] = ()) -> Dict[str, Any]:
    return provenance.provenance_block(
        name, generated_by=GENERATED_BY,
        inputs=[PREREG, B60_REPORT, B60_START1, fit_rows.FIT_ROWS_DIR / f"{RES_TABLE}.yml",
                fit_rows.FIT_ROWS_DIR / f"{D1_TABLE}.yml", ROOT / "src/kinetic_core/sulfur.py",
                ROOT / "src/kinetic_core/parameters_sulfur.py", ROOT / "src/kinetic_core/fitting.py",
                ROOT / "src/kinetic_core/config.py", ROOT / "scripts/generators/build_submodel_probe.py",
                Path(__file__)] + list(extra))


def check_freeze(commit: str, root: Path = ROOT, prereg: Path = PREREG) -> str:
    """--fit / --score: the freeze commit exists, is an ancestor of HEAD, and the prereg is unchanged since it
    (committed and working tree alike). Returns the full hash."""
    if not commit:
        raise SystemExit("--freeze-commit is required: the prereg is frozen first (prereg sec. 10)")
    full = subprocess.run(["git", "rev-parse", "--verify", f"{commit}^{{commit}}"], cwd=str(root),
                          capture_output=True, text=True)
    if full.returncode != 0:
        raise SystemExit(f"--freeze-commit {commit!r} is not a commit")
    full_hash = full.stdout.strip()
    if subprocess.run(["git", "merge-base", "--is-ancestor", full_hash, "HEAD"], cwd=str(root)).returncode != 0:
        raise SystemExit(f"--freeze-commit {full_hash} is not an ancestor of HEAD")
    rel = str(Path(prereg).resolve().relative_to(Path(root).resolve()))
    changed = subprocess.run(["git", "diff", "--name-only", full_hash, "--", rel], cwd=str(root),
                             capture_output=True, text=True).stdout.split()
    if changed:
        raise SystemExit(f"the prereg changed since the freeze commit {full_hash[:10]}: {changed}")
    return full_hash


def read_tau(cli_value: Optional[float], path: Path = TAU_FILE) -> Optional[float]:
    """H0's tau_ABS(all) (EX-059): --tau-abs, else ``D2.tau_ABS_all_used`` of the scored nested-tau record."""
    if cli_value is not None:
        return float(cli_value)
    rec = _read(path)
    v = (rec.get("D2") or {}).get("tau_ABS_all_used")
    return None if v is None or not math.isfinite(float(v)) else float(v)


def score_refusals(freeze_commit: str, tau: Optional[float], records: Mapping[str, Sequence[Mapping[str, Any]]],
                   carried: Mapping[str, Mapping[str, Any]], freeze_check: Callable[[str], str] = check_freeze
                   ) -> List[str]:
    """Every reason --score refuses (empty: it may run). The prereg unchanged since the freeze; both rivals' fit
    records with status > 0; no record with failed_pot_runs > 0; tau_ABS(all) supplied (sec. 6: never tau = 0)."""
    problems = []
    try:
        freeze_check(freeze_commit)
    except SystemExit as exc:
        problems.append(str(exc))
    for name in ("R0", "C1"):
        recs = records.get(name) or []
        if not recs:
            problems.append(f"{name}: no fit record")
        elif not any(int(r.get("status", 0)) > 0 for r in recs):
            problems.append(f"{name}: no fit record with status > 0")
        for r in recs:
            ok, why = failure_accepted(r)
            if not ok:
                problems.append(f"{name} start {r.get('start')}: failed_pot_runs {r['failed_pot_runs']} > 0, {why} (re-run it)")
        c = carried.get(name)
        if not c:
            problems.append(f"{name}: no carried record (run --laplace {name})")
        elif int(c.get("failed_pot_runs", 0)) > 0:
            problems.append(f"{name}: the carried record has failed_pot_runs {c['failed_pot_runs']} > 0")
    if tau is None:
        problems.append("no tau_ABS(all): pass --tau-abs or put EX-059's scored record on disk (sec. 6: the scoring "
                        "waits; it never runs at tau = 0)")
    elif not (tau > 0.0 and math.isfinite(tau)):
        problems.append(f"tau_ABS(all) = {tau!r}: must be finite and > 0")
    return problems


def fit_refusals(freeze_commit: str, checks: Mapping[str, Any], freeze_check: Callable[[str], str] = check_freeze
                 ) -> List[str]:
    problems = []
    try:
        freeze_check(freeze_commit)
    except SystemExit as exc:
        problems.append(str(exc))
    c = checks.get("checks") or {}
    for name in ("V1", "V2", "V5", "V6", "V7", "V8"):
        if name not in c:
            problems.append(f"{name} has not been run (--checks first)")
    if "V7" in c and not c["V7"].get("passed"):
        problems.append("V7 failed: a build error, investigated before any fit (prereg sec. 4.4)")
    if "V5" in c and not c["V5"].get("integrated_without_failure"):
        problems.append("V5: an integration at the ceiling failed (a build error to fix before any fit)")
    if "V1" in c and not c["V1"].get("passed"):
        problems.append("V1: B60's cost does not reproduce")
    return problems


def numerics_for_fits(checks: Mapping[str, Any]) -> Numerics:
    """sec. 4.4: if V8 fails, every fit runs on fixed_point_nodes at the quick setting; the scoring stays inline."""
    v8 = (checks.get("checks") or {}).get("V8") or {}
    return fit_inline() if v8.get("passed") else fit_fixed_point()


def fit_side_from_checks(checks: Mapping[str, Any]) -> FitSide:
    v5 = (checks.get("checks") or {}).get("V5") or {}
    return FitSide.on(float(v5.get("mft_cap_mmol_l", 1.0)))


# ----------------------------------------------------------------------------------------------------------------
# Starts (sec. 4.3)
# ----------------------------------------------------------------------------------------------------------------

def perturbed(theta0: np.ndarray, keys: Sequence[str], lower: np.ndarray, upper: np.ndarray, start: int
              ) -> Tuple[np.ndarray, int]:
    """B8's perturbation protocol (generate_kinetic_core_b8_fit.start_vector) on M1's free coordinates: a seeded
    uniform draw within B8's half-widths per coordinate class (2 decades, 40 kJ/mol, 0.15, 1 pH), clipped to the
    bounds. Seed B8.SEED_BASE + B8.SEED_STRIDE x start."""
    seed = B8.SEED_BASE + B8.SEED_STRIDE * start
    rng = np.random.default_rng(seed)
    x = np.array(theta0, dtype=float)
    for i, k in enumerate(keys):
        kind = kind_of(k)
        half = {"Ea": B8.LOCAL_HALF_WIDTH_EA_KJ, "acid_yield": B8.LOCAL_HALF_WIDTH_ACID_YIELD,
                "pka": B8.LOCAL_HALF_WIDTH_PKA}.get(kind, B8.LOCAL_HALF_WIDTH_DECADES)
        x[i] = float(rng.uniform(max(lower[i], theta0[i] - half), min(upper[i], theta0[i] + half)))
    return x, seed


def lowest_converged(rival: str) -> Optional[Dict[str, Any]]:
    recs = [r for r in (_read(fit_path(rival, s)) for s in range(3)) if r.get("complete")]
    ok = [r for r in recs if int(r.get("status", 0)) > 0]
    return min(ok, key=lambda r: float(r["cost"])) if ok else None


def start_theta(obj: Objective, rival: str, start: int) -> Tuple[np.ndarray, Dict[str, Any]]:
    """sec. 4.3's table. R0: s0 B60's optimum; s1 B8's perturbation of s0; s2 B60's start-1 endpoint (its free
    coordinates; the held ones stay at B60's best, as 'held at B60's values' says). C1: s0 B60's vector with the new
    coordinates at their sec. 2 starts; s1 the perturbation; s2 R0's carried optimum (its lowest converged start,
    which is the carried optimum when R0 is compared and its declared fallback when it is not) + the new coordinates."""
    theta0 = obj.theta_of(b60_point(obj.rival.route))
    if start == 0:
        return theta0, {"origin": "B60's optimum" + (", new coordinates at their sec. 2 starts" if obj.rival.route else "")}
    if start == 1:
        x, seed = perturbed(theta0, obj.keys, obj.lower, obj.upper, 1)
        return x, {"origin": "B8's perturbation of start 0", "seed": seed}
    if start == 2 and rival == "R0":
        member = json.loads(B60_START1.read_text(encoding="utf-8"))
        x1 = np.array(member["x_full"], dtype=float)
        th = np.array([x1[all_index(k)] if k in B60.FREE_KEYS else theta0[i] for i, k in enumerate(obj.keys)])
        return np.clip(th, obj.lower, obj.upper), {"origin": f"B60's start-1 endpoint ({data_paths.rel(B60_START1)})"}
    if start == 2 and rival == "C1":
        best = lowest_converged("R0")
        if best is None:
            raise SystemExit("C1's start 2 waits for R0's carried optimum: no converged R0 start yet")
        th = np.array([best["theta"][k] if k in best["theta"] else theta0[i] for i, k in enumerate(obj.keys)])
        return np.clip(th, obj.lower, obj.upper), {"origin": f"R0's carried optimum (start {best['start']})"}
    raise SystemExit(f"start {start}: the prereg declares starts 0 - 2")


def start0_noise(fit_side: FitSide, numerics: Numerics) -> Dict[str, float]:
    """sec. 4.1: s_series per data set x species by fitting.estimate_series_noise at C1's start 0 (held for its starts)."""
    obj = Objective(RIVALS["C1"], numerics, fit_side)
    point = obj.point(obj.theta_of(b60_point(True)))
    failures: List[Dict[str, Any]] = []
    cfg = config_of(point, numerics, fit_side)
    with config.applied(cfg):
        raw = route_raw(point, numerics, cfg, failures)
    FAILURES.extend(failures)
    rows = route_block()[0]
    ok = [i for i, v in enumerate(raw) if np.isfinite(v)]
    return fitting.estimate_series_noise(route_model(), raw[ok], tuple(rows[i] for i in ok), len(obj.keys))


# ----------------------------------------------------------------------------------------------------------------
# --fit
# ----------------------------------------------------------------------------------------------------------------

def run_fit(rival_name: str, start: int, freeze_commit: str, budget: int = BUDGET) -> Dict[str, Any]:
    checks = _read(CHECKS_OUT)
    problems = fit_refusals(freeze_commit, checks)
    if problems:
        raise SystemExit("--fit refuses:\n  " + "\n  ".join(problems))
    commit = check_freeze(freeze_commit)
    rival = RIVALS[rival_name]
    numerics, fit_side = numerics_for_fits(checks), fit_side_from_checks(checks)
    out = fit_path(rival_name, start)
    state = _read(out)
    if int(state.get("failed_pot_runs", 0)) > 0 and failure_accepted(state)[0]:
        print(f"{out.name}: {failure_accepted(state)[1]}; not re-run", flush=True)
        return state
    if int(state.get("failed_pot_runs", 0)) > 0:   # sec. 4.3: re-run with the same settings; the failed one kept, labelled
        k = 1
        while out.with_name(f"{out.stem}_failed_run{k}.json").exists():
            k += 1
        label = out.with_name(f"{out.stem}_failed_run{k}.json")
        state["label"] = f"FAILED POT RUNS ({state['failed_pot_runs']}): kept, not used; re-run as {out.name}"
        _write(label, state)
        print(f"{out.name}: failed_pot_runs {state['failed_pot_runs']} -> kept as {label.name}; re-running", flush=True)
        state = {}
    if state.get("complete"):
        print(f"{out.name}: complete already (status {state['status']}, cost {state['cost']:.4f})")
        return state
    noise = dict(state.get("series_noise") or {})
    if rival.route and not noise:
        noise = start0_noise(fit_side, numerics)
    obj = Objective(rival, numerics, fit_side, noise)
    theta0, origin = start_theta(obj, rival_name, start)
    failed_before = int(state.get("failed_pot_runs", 0))
    del FAILURES[:]
    state.update({
        "artifact": out.stem, "prereg": data_paths.rel(PREREG), "freeze_commit": commit, "rival": rival_name,
        "start": start, "start_origin": origin, "keys": list(obj.keys), "kinds": list(obj.kinds),
        "lower": dict(zip(obj.keys, obj.lower.tolist())), "upper": dict(zip(obj.keys, obj.upper.tolist())),
        "start_theta": dict(zip(obj.keys, theta0.tolist())), "numerics": dataclasses.asdict(numerics),
        "dimer_ceiling_log10_145C": dimer_ceiling(fit_side.ceiling, rival.route),
        "mft_cap_mmol_l": fit_side.ceiling.mft_cap_mmol_l, "series_noise": noise, "budget_per_segment": budget,
        "config_digest_at_start": config_of(obj.point(theta0), numerics, fit_side).digest(),
        "pot_seconds_limit": POT_SECONDS, "workers": parallel.workers()})
    history = state.setdefault("history", [])
    theta = np.array([history[-1]["theta"][k] for k in obj.keys]) if history else theta0
    while not (history and (int(history[-1]["status"]) > 0 or len(history) >= MAX_SEGMENTS)):
        n0 = len(FAILURES)
        f_ = fitting.fit(obj.problem(), theta, budget=budget)
        seg = f_.as_dict(obj.keys)
        seg["failed_pot_runs_segment"] = len(FAILURES) - n0
        history.append(seg)
        theta = f_.theta
        state["failed_pot_runs"] = failed_before + len(FAILURES)
        state["failures_sample"] = (state.get("failures_sample") or []) + FAILURES[n0:n0 + 10]
        _write(out, state)   # checkpoint after every segment
        print(f"{rival_name} s{start} segment {len(history)}: status {f_.status}, cost {f_.cost:.4f}, "
              f"{f_.evaluations} ev, {f_.seconds} s, failed pots {seg['failed_pot_runs_segment']}", flush=True)
    last = history[-1]
    r = obj.problem().evaluate(theta)
    state.update({"status": int(last["status"]), "cost": float(last["cost"]), "theta": dict(last["theta"]),
                  "converged": int(last["status"]) > 0,
                  "not_converged_label": None if int(last["status"]) > 0 else
                  "not converged (status 0 after its continuation): not used (sec. 4.3)",
                  "evaluations": sum(int(h["evaluations"]) for h in history),
                  "seconds": sum(float(h["seconds"]) for h in history),
                  "cost_by_block": obj.blocks(r), "residual_by_row": dict(zip(obj.row_labels, map(float, r))),
                  "failed_pot_runs": failed_before + len(FAILURES), "provenance": _provenance(out.stem),
                  "complete": True})
    _write(out, state)
    print(f"{out.name}: status {state['status']}, cost {state['cost']:.4f}, failed pot runs {state['failed_pot_runs']}")
    return state


# ----------------------------------------------------------------------------------------------------------------
# --laplace: comparison rule, C1's noise continuation, Laplace, profiles, G1 (sec. 4.3, sec. 7)
# ----------------------------------------------------------------------------------------------------------------

def decisive_tolerance(levels: np.ndarray, jac: Optional[np.ndarray], sigma_theta: Optional[np.ndarray]
                       ) -> Tuple[np.ndarray, np.ndarray]:
    """The decisive POINTS (floored D1 ratios and reservoir levels) and their agreement tolerance max(sigma_decision,
    Laplace predictive sd); sigma_decision is the per-ratio decision sd for D1 and the row's sigma (no tau: a shared
    level) for the reservoir."""
    us = units(levels, jac, Scoring(tau=0.0))
    pts, tol = [], []
    for name, u in us.items():
        pts += list(u["f"])
        sd_dec = np.sqrt(np.diag(u["Sigma"]))
        sd_pred = (np.sqrt(np.clip(np.diag(u["J"] @ sigma_theta @ u["J"].T), 0.0, None))
                   if (sigma_theta is not None and u["J"] is not None) else np.zeros(len(u["f"])))
        tol += list(np.maximum(sd_dec, sd_pred))
    return np.array(pts), np.array(tol)


def central_jacobian(fn: Callable[[List[np.ndarray]], List[np.ndarray]], theta: np.ndarray, kinds: Sequence[str],
                     lower: np.ndarray, upper: np.ndarray, coords: Sequence[int]) -> np.ndarray:
    """Central differences of a vector function at ``theta`` on ``coords`` with the module's steps (one-sided at a
    bound), all points evaluated in one batch."""
    pts, plan = [], []
    for j in coords:
        h = fitting._step(kinds[j])
        up, down = theta[j] + h <= upper[j], theta[j] - h >= lower[j]
        xp, xm = theta.copy(), theta.copy()
        xp[j] += h
        xm[j] -= h
        plan.append((j, h, up, down))
        pts += [xp if up else theta.copy(), xm if down else theta.copy()]
    vals = fn(pts)
    cols = []
    for k, (j, h, up, down) in enumerate(plan):
        vp, vm = vals[2 * k], vals[2 * k + 1]
        width = (h if up else 0.0) + (h if down else 0.0)
        cols.append((vp - vm) / width if width > 0 else np.zeros_like(vp))
    return np.column_stack(cols) if cols else np.zeros((len(point_labels()), 0))


def identified_parts(lap: Mapping[str, Any], keys: Sequence[str], profiles: Mapping[str, Any]
                     ) -> Tuple[List[int], np.ndarray]:
    """sec. 7: the coordinates that enter the linearised predictive and Sigma_theta over them: identified ones with
    the Laplace covariance; bound-active ones with their profile variance, uncorrelated (declared approximation);
    unidentified held (left out)."""
    ident = [i for i, k in enumerate(keys) if lap["identified"][k] and not lap["at_bound"][k]]
    bound = [i for i, k in enumerate(keys) if lap["at_bound"][k] and k in profiles and "variance" in profiles[k]]
    cov = np.asarray(lap["covariance"], dtype=float)
    idx = ident + bound
    out = np.zeros((len(idx), len(idx)))
    out[:len(ident), :len(ident)] = cov[np.ix_(ident, ident)]
    for n, i in enumerate(bound):
        out[len(ident) + n, len(ident) + n] = float(profiles[keys[i]]["variance"])
    return idx, out


def run_profiles(obj: Objective, theta: np.ndarray, lap: Mapping[str, Any], c_min: float, record: Dict[str, Any],
                 save: Callable[[], None]) -> Dict[str, Any]:
    profiles = record.setdefault("profiles", {})
    for i, k in enumerate(obj.keys):
        if not lap["at_bound"][k]:
            continue
        kind = obj.kinds[i]
        if kind not in PROFILE_OFFSETS:
            profiles[k] = {"note": f"on a bound; the prereg declares no profile for a {kind} coordinate: held at its "
                                   "value in the linearised predictive (no variance), flagged"}
            save()
            continue
        side = "upper" if (obj.upper[i] - theta[i]) <= (theta[i] - obj.lower[i]) else "lower"
        bound = float(obj.upper[i] if side == "upper" else obj.lower[i])
        prof = profiles.setdefault(k, {"side": side, "bound": bound, "points": []})
        done = {round(p["offset"], 6) for p in prof["points"]}
        for off in PROFILE_OFFSETS[kind]:
            if round(off, 6) in done:
                continue
            v = bound - off if side == "upper" else bound + off
            sub = Objective(obj.rival, obj.numerics, obj.fit_side, dict(obj.spec.noise), base=obj.base,
                            fixed={**obj.fixed, k: v})
            th = np.array([theta[obj.keys.index(kk)] for kk in sub.keys])
            n0 = len(FAILURES)
            f_ = fitting.fit(sub.problem(), th, budget=PROFILE_BUDGET)
            segs = [f_.as_dict(sub.keys)]
            if f_.status == 0:
                f_ = fitting.fit(sub.problem(), f_.theta, budget=PROFILE_BUDGET)
                segs.append(f_.as_dict(sub.keys))
            prof["points"].append({"offset": off, "value": v, "cost": float(f_.cost), "status": int(f_.status),
                                   "segments": segs, "failed_pot_runs": len(FAILURES) - n0,
                                   "flag": None if f_.status > 0 else
                                   "not converged: its cost is an upper bound, its density a lower bound"})
            save()
        pts = sorted(prof["points"], key=lambda p: p["offset"])
        far = float(obj.lower[i] if side == "upper" else obj.upper[i])
        prof.update(profile_density(bound, side, [p["value"] for p in pts], [p["cost"] for p in pts], c_min, far))
        save()
    return profiles


def g1_at(point: Point, numerics: Numerics) -> Dict[str, Any]:
    """M1.0's G1 (B27's gate) at an optimum: oxidant use and the phi = 1 rise in the three ambient pots."""
    cfg = config_of(point, numerics, FitSide.off())
    x = point.x()
    params, drift = B23.build_parameters(x), B23.unpack(x)[3]
    with config.applied(cfg):
        pots = {name: M10.g1_pot(params, drift, name) for name in M10.B27SR.AMBIENT_POTS}
    m10 = _read(PREFLIGHT).get("G1", {}).get("pots", {})
    return {"pots": {n: {k: p[k] for k in ("consumed_fraction", "rise_at_phi_1_b27_sum", "rise_at_phi_1_r_nf_mp3p_only",
                                           "dimer_flux_mmol_l")} for n, p in pots.items()},
            "m1_0_b60": {n: {"consumed_fraction": v["B60"]["consumed_fraction"],
                             "rise_at_phi_1_b27_sum": v["B60"]["rise_at_phi_1_b27_sum"]} for n, v in m10.items()}}


def laplace_with_cov(obj: Objective, theta: np.ndarray) -> Dict[str, Any]:
    return fitting.laplace(obj.problem(), theta, covariance=True)


def decisive_jac_at(obj: Objective, theta: np.ndarray, coords: Sequence[int], reading: str, numerics: Numerics
                    ) -> Tuple[np.ndarray, np.ndarray]:
    lv = decisive_many([obj.point(theta)], reading, numerics)[0]
    jac = central_jacobian(lambda pts: decisive_many([obj.point(p) for p in pts], reading, numerics), theta,
                           obj.kinds, obj.lower, obj.upper, coords)
    return lv, jac


def run_laplace(rival_name: str, freeze_commit: str) -> Dict[str, Any]:
    checks = _read(CHECKS_OUT)
    problems = fit_refusals(freeze_commit, checks)
    if problems:
        raise SystemExit("--laplace refuses:\n  " + "\n  ".join(problems))
    commit = check_freeze(freeze_commit)
    rival = RIVALS[rival_name]
    numerics, fit_side, scoring = numerics_for_fits(checks), fit_side_from_checks(checks), score_inline()
    recs = [_read(fit_path(rival_name, s)) for s in range(3)]
    missing = [s for s, r in enumerate(recs) if not r.get("complete")]
    if missing:
        raise SystemExit(f"{rival_name}: starts {missing} have not ended")
    bad = [r["start"] for r in recs if not failure_accepted(r)[0]]
    for r in recs:
        if int(r.get("failed_pot_runs", 0)) > 0:
            print(f"{rival_name} s{r['start']}: {failure_accepted(r)[1]}", flush=True)
    if bad:
        raise SystemExit(f"{rival_name}: starts {bad} have failed_pot_runs > 0: re-run them with --fit first")
    out = carried_path(rival_name)
    state = _read(out)
    del FAILURES[:]

    def save():
        state["failed_pot_runs"] = int(state.get("failed_pot_runs_before", 0)) + len(FAILURES)
        _write(out, state)

    state.setdefault("failed_pot_runs_before", int(state.get("failed_pot_runs", 0)))
    state.update({"artifact": out.stem, "rival": rival_name, "freeze_commit": commit, "prereg": data_paths.rel(PREREG),
                  "numerics_fits": dataclasses.asdict(numerics), "numerics_scoring": dataclasses.asdict(scoring)})
    noise0 = dict(recs[0].get("series_noise") or {})
    obj = Objective(rival, numerics, fit_side, noise0)
    rule = comparison_rule([{"start": r["start"], "status": r["status"], "cost": r["cost"]} for r in recs])
    state["comparison_a"] = rule
    if rule["best_start"] is None:      # no converged start
        if rival_name == "C1":
            state.update(status="not_compared", note="no C1 start converged: no verdict and no switch changes (sec. 4.3)")
            save()
            return state
        state.update(status="fallback_R0_shipped_b",
                     note="no R0 start converged: the verdict uses R0-shipped (b), flagged 'different dimer bound and "
                          "k_mft_decay treatment' (review B4)", carried_theta=None)
        save()
        return state
    theta_of = {int(r["start"]): np.array([r["theta"][k] for k in obj.keys]) for r in recs}
    best = rule["best_start"]
    if "laplace_best" not in state:
        state["laplace_best"] = laplace_with_cov(obj, theta_of[best])
        save()
    lap_b = state["laplace_best"]
    if rule["a"] and "comparison_b" not in state:
        coords = [i for i, k in enumerate(obj.keys) if lap_b["identified"][k] and not lap_b["at_bound"][k]]
        lv, jac = decisive_jac_at(obj, theta_of[best], coords, "isothermal", scoring)
        cov = np.asarray(lap_b["covariance"])[np.ix_(coords, coords)]
        _, tol = decisive_tolerance(lv, jac, cov)
        preds = {}
        for s in rule["qualifying_starts"]:
            lvs = lv if s == best else decisive_many([obj.point(theta_of[s])], "isothermal", scoring)[0]
            preds[s] = decisive_tolerance(lvs, None, None)[0]
        state["comparison_b"] = comparison_rule([{"start": r["start"], "status": r["status"], "cost": r["cost"]}
                                                 for r in recs], preds, tol)
        state["comparison_b"]["predictions"] = {str(k): v.tolist() for k, v in preds.items()}
        state["comparison_b"]["tolerance"] = tol.tolist()
        save()
    status = state["comparison_b"]["status"] if rule["a"] else "not_compared"
    if status == "not_compared":
        status = ("fallback_lowest_converged_start" if rival_name == "R0" else "not_compared")
    state["status"] = status
    if status == "not_compared":     # C1 not compared: no verdict
        state["note"] = ("C1 not compared, flagged: every start's cost, vector and decisive predictions are listed; no "
                         "verdict and no switch changes (sec. 4.3)")
        state["starts"] = [{k: r[k] for k in ("start", "status", "cost", "theta")} for r in recs]
        save()
        return state
    carried, noise = theta_of[best], noise0
    if rival_name == "C1" and "noise_continuation" not in state:
        # sec. 4.3 (review M11; re-review MINOR 5): s_series re-estimated at the carried optimum, the fit continued
        failures: List[Dict[str, Any]] = []
        point = obj.point(carried)
        cfg = config_of(point, numerics, fit_side)
        with config.applied(cfg):
            raw = route_raw(point, numerics, cfg, failures)
        FAILURES.extend(failures)
        rows = route_block()[0]
        ok = [i for i, v in enumerate(raw) if np.isfinite(v)]
        noise1 = fitting.estimate_series_noise(route_model(), raw[ok], tuple(rows[i] for i in ok), len(obj.keys))
        obj1 = Objective(rival, numerics, fit_side, noise1)
        segs, th = [], carried
        for _ in range(2):
            f_ = fitting.fit(obj1.problem(), th, budget=NOISE_CONTINUATION)
            segs.append(f_.as_dict(obj1.keys))
            th = f_.theta
            if f_.status > 0:
                break
        ok_cont = segs[-1]["status"] > 0
        state["noise_continuation"] = {
            "noise_start": noise0, "noise_reestimated": noise1, "segments": segs, "status": segs[-1]["status"],
            "used": ok_cont, "flag": None if ok_cont else
            "status 0 twice: the pre-re-estimation optimum is carried, flagged, with both noise sets recorded"}
        save()
    if rival_name == "C1" and state["noise_continuation"]["used"]:
        nc = state["noise_continuation"]
        noise = nc["noise_reestimated"]
        carried = np.array([nc["segments"][-1]["theta"][k] for k in obj.keys])
    obj_c = Objective(rival, numerics, fit_side, noise)
    state.update(carried_theta=dict(zip(obj.keys, carried.tolist())), carried_noise=noise,
                 carried_from_start=best)
    if "laplace" not in state:
        state["laplace"] = laplace_with_cov(obj_c, carried)
        save()
    lap = state["laplace"]
    r = obj_c.problem().evaluate(carried)
    c_min = 0.5 * float(r @ r)
    state["carried_cost"] = c_min
    state["cost_by_block"] = obj_c.blocks(r)
    save()
    run_profiles(obj_c, carried, lap, c_min, state, save)
    if status == "compared_if_verdicts_agree":
        per = state.setdefault("laplace_per_qualifying_start", {})
        for s in state["comparison_b"]["qualifying_starts"]:
            if str(s) not in per:
                per[str(s)] = {"theta": dict(zip(obj.keys, theta_of[s].tolist())),
                               "laplace": laplace_with_cov(obj, theta_of[s])}
                save()
    if "G1" not in state:
        state["G1"] = g1_at(obj_c.point(carried), scoring)
        save()
    i_ea = obj.keys.index("Ea_decay_thiol_sink")
    grad = lambda_log10_gradient()
    gi = np.array([grad.get(k, 0.0) for k in obj.keys])
    cov = np.asarray(lap["covariance"])
    state["lambda_mft_100C"] = {"log10": prior_log10_rate(obj_c.point(carried).x(), fit_side.prior),
                                "laplace_sd_dex": float(math.sqrt(max(gi @ cov @ gi, 0.0)))}
    state["gigl_ceiling_binding"] = bool(lap["at_bound"]["Ea_decay_thiol_sink"]
                                         and carried[i_ea] >= obj.upper[i_ea] - 1e-3 * (obj.upper[i_ea] - obj.lower[i_ea]))
    state["provenance"] = _provenance(out.stem, [fit_path(rival_name, s) for s in range(3)])
    state["complete"] = True
    save()
    print(f"{out.name}: {status}; carried cost {c_min:.4f}; on bound "
          f"{[k for k, v in lap['at_bound'].items() if v]}; failed pot runs {state['failed_pot_runs']}")
    return state


# ----------------------------------------------------------------------------------------------------------------
# --score: the decisive comparison (sec. 7) and the side effects (sec. 8)
# ----------------------------------------------------------------------------------------------------------------

@dataclass
class Form:
    name: str
    point: Point
    numerics: Numerics
    obj: Optional[Objective] = None
    theta: Optional[np.ndarray] = None
    lap: Optional[Dict[str, Any]] = None
    profiles: Dict[str, Any] = field(default_factory=dict)
    flag: Optional[str] = None


def forms_for_score(checks: Mapping[str, Any]) -> Tuple[Dict[str, Form], Dict[str, Any]]:
    numerics, fit_side = numerics_for_fits(checks), fit_side_from_checks(checks)
    forms = {"R0_shipped_a": Form("R0_shipped_a", b60_point(), shipped_fixed_point(), flag="printed only"),
             "R0_shipped_b": Form("R0_shipped_b", b60_point(), score_inline())}
    info = {}
    for name in ("R0", "C1"):
        c = _read(carried_path(name))
        info[name] = {k: c.get(k) for k in ("status", "note", "carried_from_start", "gigl_ceiling_binding",
                                             "lambda_mft_100C", "cost_by_block")}
        if c.get("status") == "fallback_R0_shipped_b":
            forms["R0"] = dataclasses.replace(forms["R0_shipped_b"], name="R0",
                                              flag="R0-shipped (b): different dimer bound and k_mft_decay treatment")
            continue
        if c.get("status") == "not_compared":
            continue
        obj = Objective(RIVALS[name], numerics, fit_side, c.get("carried_noise") or {})
        th = np.array([c["carried_theta"][k] for k in obj.keys])
        forms[name] = Form(name, obj.point(th), score_inline(), obj, th, c["laplace"], c.get("profiles") or {},
                           flag="lowest converged start, flagged" if c["status"] == "fallback_lowest_converged_start"
                           else None)
    return forms, info


def form_levels(form: Form, readings: Sequence[str], **kw) -> Dict[str, Dict[str, Any]]:
    """Per reading: the raw log10 levels and (for a fitted form) their Jacobian on the coordinates that enter the
    linearised predictive."""
    out = {}
    coords, sig = [], None
    if form.lap is not None:
        coords, sig = identified_parts(form.lap, form.obj.keys, form.profiles)
    for reading in readings:
        if form.obj is None:
            lv = decisive_many([form.point], reading, form.numerics, **kw)[0]
            out[reading] = {"levels": lv, "jac": None, "sigma_theta": None}
            continue
        lv, jac = decisive_jac_at(form.obj, form.theta, coords, reading, form.numerics) if coords else (
            decisive_many([form.point], reading, form.numerics, **kw)[0], None)
        out[reading] = {"levels": lv, "jac": jac, "sigma_theta": sig}
    return out


def mc_draws(form: Form, rng: np.random.Generator, wide: bool = False) -> List[np.ndarray]:
    """sec. 7: S draws. Identified coordinates from the Laplace (reflected at bounds), bound-active from the profile,
    unidentified held (wide: drawn uniform across their bounds -- the prior-predictive sensitivity)."""
    keys = form.obj.keys
    lap = form.lap
    ident = [i for i, k in enumerate(keys) if lap["identified"][k] and not lap["at_bound"][k]]
    cov = np.asarray(lap["covariance"])[np.ix_(ident, ident)]
    draws = np.tile(form.theta, (MC_DRAWS, 1))
    if ident:
        z = rng.multivariate_normal(form.theta[ident], cov, size=MC_DRAWS, check_valid="ignore")
        draws[:, ident] = reflect(z, form.obj.lower[ident], form.obj.upper[ident])
    for i, k in enumerate(keys):
        if lap["at_bound"][k] and "density" in form.profiles.get(k, {}):
            draws[:, i] = sample_profile(form.profiles[k], rng, MC_DRAWS)
        elif wide and not lap["identified"][k] and not lap["at_bound"][k]:
            draws[:, i] = rng.uniform(form.obj.lower[i], form.obj.upper[i], size=MC_DRAWS)
    return list(draws)


def score_units(lv: Dict[str, Any], sc: Scoring, mc_levels: Optional[Sequence[np.ndarray]] = None
                ) -> Dict[str, Dict[str, Any]]:
    us = units(lv["levels"], lv["jac"], sc)
    out = {}
    mcu = [units(m, None, sc) for m in mc_levels] if mc_levels else None
    for name, u in us.items():
        rec = {**unit_lppd(u, lv["sigma_theta"]), "f": u["f"].tolist(), "y": u["y"].tolist(),
               "n_floored": u["n_floored"], "block": u["block"], "labels": u["labels"]}
        if mcu:
            rec.update(mc_lppd(u["y"], [m[name]["f"] for m in mcu], u["Sigma"]))
        out[name] = rec
    return out


def decide(per_form: Mapping[str, Mapping[str, Any]], c1: str = "C1", r0: str = "R0") -> Dict[str, Any]:
    """d_u = lppd(C1) - lppd(R0) per unit, the tie band, the blocks, the outcome: one reading."""
    prim = "primary"
    u_names = list(per_form[c1][prim].keys())
    rec: Dict[str, Any] = {"units": {}, "blocks": {}}
    for u in u_names:
        signs = {"linearised": per_form[c1][prim][u]["linearised"] - per_form[r0][prim][u]["linearised"],
                 "plug_in": per_form[c1][prim][u]["plug_in"] - per_form[r0][prim][u]["plug_in"]}
        for sens in ("floors_alternative", "tau_x0.5", "tau_x2"):
            signs[sens] = per_form[c1][sens][u]["linearised"] - per_form[r0][sens][u]["linearised"]
        rec["units"][u] = {"d": signs, "side": unit_side(list(signs.values())),
                           "block": per_form[c1][prim][u]["block"]}
    for block in ("reservoir", "D1"):
        d = [v["d"]["linearised"] for v in rec["units"].values() if v["block"] == block]
        rec["blocks"][block] = {"Delta_L": float(np.sum(d)), "SE": block_se(d),
                                "sides": [v["side"] for v in rec["units"].values() if v["block"] == block]}
    rec["outcome"] = block_outcome(rec["blocks"]["reservoir"]["sides"], rec["blocks"]["D1"]["sides"])
    return rec


CF063_NOTE = ("The reservoir block is one pot from one laboratory (Zhang 2025b, BTBU) whose Cys and Xyl print one "
              "identical ELSD calibration curve (CF-063): 'both units' is close to one datum, so any supply credit "
              "rests on Zhang 2025b alone.")
CREDIT_NOTE = ("Supply is credited on the reservoir block only; a D1 gain of C1 is 'consistent with supply, not "
               "attributed' (C2 / R2 not run, sec. 5.1).")


def verdict_text(verdict: str, gigl: bool, extra: Sequence[str] = ()) -> str:
    parts = [f"Verdict: {verdict}."]
    if verdict == NOT_SEPARATED:
        parts.append("CF-107 stays open, no switch changes, and the result is not evidence against C1.")
    if verdict == R0_FAVOURED:
        parts.append("No switch changes and CF-107 stays open; the loss law is not decided.")
    if verdict == C1_FAVOURED:
        parts.append("C1's structure ships only if C1 was compared, the side-effect record is complete with a "
                     "candidate compensating error registered for each loss, no review objection is unanswered, and "
                     "the coupling flip ships first or is declared part of the ship (sec. 7 rule 7).")
    if gigl:
        parts.append("D1 is conditional on the Gigl ceiling (Ea_decay_thiol_sink binds at 102 kJ/mol in a rival).")
    parts += [CREDIT_NOTE, CF063_NOTE, *extra]
    return " ".join(parts)


def run_score(freeze_commit: str, tau_cli: Optional[float], side_effects: bool = True) -> Dict[str, Any]:
    checks = _read(CHECKS_OUT)
    tau = read_tau(tau_cli)
    records = {n: [r for r in (_read(fit_path(n, s)) for s in range(3)) if r] for n in ("R0", "C1")}
    carried = {n: _read(carried_path(n)) for n in ("R0", "C1")}
    problems = score_refusals(freeze_commit, tau, records, carried)
    if problems:
        raise SystemExit("--score refuses:\n  " + "\n  ".join(problems))
    commit = check_freeze(freeze_commit)
    forms, info = forms_for_score(checks)
    del FAILURES[:]
    result: Dict[str, Any] = {"artifact": COMPARISON_OUT.stem, "prereg": data_paths.rel(PREREG), "freeze_commit": commit,
                              "tau_ABS_all": tau, "tau_source": "--tau-abs" if tau_cli is not None else data_paths.rel(TAU_FILE),
                              "rivals": info, "forms": {k: {"flag": f.flag, "numerics": dataclasses.asdict(f.numerics)}
                                                        for k, f in forms.items()}}
    if "C1" not in forms:
        result["verdict"] = {"verdict": "no verdict (C1 not compared)", "text": "No verdict and no switch changes."}
        artifact_io.write_artifact(result, COMPARISON_OUT, render=render_comparison)
        return result
    settings = {"primary": Scoring(tau), "floors_alternative": Scoring(tau, floors="alternative"),
                "tau_x0.5": Scoring(0.5 * tau), "tau_x2": Scoring(2.0 * tau),
                "sigma_read_error": Scoring(tau, d1_sigma="read_error"), "sigma_s042": Scoring(tau, d1_sigma="s042")}
    rng = np.random.default_rng(MC_SEED)
    levels, scored = {}, {}
    for name, form in forms.items():
        levels[name] = form_levels(form, READINGS)
        draws = mc_draws(form, rng) if form.lap is not None else None
        wide = mc_draws(form, rng, wide=True) if form.lap is not None else None
        scored[name] = {}
        for reading in READINGS:
            mc_lv = decisive_many([form.obj.point(d) for d in draws], reading, form.numerics) if draws else None
            wide_lv = decisive_many([form.obj.point(d) for d in wide], reading, form.numerics) if wide else None
            scored[name][reading] = {s: score_units(levels[name][reading], sc, mc_lv if s == "primary" else None)
                                     for s, sc in settings.items()}
            if wide_lv:
                scored[name][reading]["primary_unidentified_wide_mc"] = score_units(
                    levels[name][reading], settings["primary"], wide_lv)
    # the non-linear flag: MC and linearised differ by more than 2 MC SE
    for name in scored:
        for reading in READINGS:
            for u, rec in scored[name][reading]["primary"].items():
                if "mc" in rec and math.isfinite(rec.get("mc_se", float("nan"))):
                    rec["non_linear_predictive"] = bool(abs(rec["mc"] - rec["linearised"]) > 2.0 * rec["mc_se"])
    per_reading = {r: decide({n: scored[n][r] for n in ("C1", "R0")}) for r in READINGS}
    across = verdict_across_readings({r: v["outcome"] for r, v in per_reading.items()})
    gigl = bool(info["R0"].get("gigl_ceiling_binding") or info["C1"].get("gigl_ceiling_binding"))
    verdict = dict(across)
    extra = []
    for name in ("R0", "C1"):
        st = carried[name].get("status")
        if st == "compared_if_verdicts_agree":
            extra.append(f"{name}'s qualifying starts disagree in their decisive predictions; the verdict is computed "
                         "at each and counts only if they agree (sec. 4.3).")
        if forms[name].flag and name in ("R0",):
            extra.append(f"R0: {forms[name].flag}.")
    per_start = per_start_verdicts(forms, carried, settings, rng) if any(
        carried[n].get("status") == "compared_if_verdicts_agree" for n in ("R0", "C1")) else {}
    if per_start and len({v["verdict"] for v in per_start.values()}) > 1:
        verdict = {"verdict": NOT_SEPARATED, "agree": False, "per_start_disagree": True,
                   "disagreeing_readings": across.get("disagreeing_readings", [])}
        extra.append("The per-start verdicts disagree: the rival is not compared, flagged.")
    verdict["text"] = verdict_text(verdict["verdict"], gigl, extra)
    result.update({"readings": list(READINGS), "per_reading": per_reading, "verdict": verdict,
                   "per_start_verdicts": per_start, "scored": _jsonable(scored),
                   "levels": {n: {r: {"levels": v["levels"].tolist(),
                                      "jac": None if v["jac"] is None else v["jac"].tolist()} for r, v in lv.items()}
                              for n, lv in levels.items()},
                   "point_labels": point_labels()})
    _write(COMPARISON_OUT, result)    # checkpoint before the sensitivities
    result["sensitivities"] = sensitivities(forms, settings, gigl)
    result["reservoir_record"] = reservoir_record(forms, tau)
    if not ((checks.get("checks") or {}).get("V8") or {}).get("passed"):
        result["v7_at_carried_optima"] = {n: v7_diff(forms[n].point) for n in ("R0", "C1") if n in forms}
    result["failed_pot_runs"] = len(FAILURES)
    result["provenance"] = _provenance(COMPARISON_OUT.stem, [carried_path("R0"), carried_path("C1"), CHECKS_OUT])
    artifact_io.write_artifact(_jsonable(result), COMPARISON_OUT, render=render_comparison)
    print(verdict["text"])
    if side_effects:
        run_side_effects(forms, commit)
    return result


def per_start_verdicts(forms, carried, settings, rng) -> Dict[str, Any]:
    """sec. 4.3: when (a) holds and (b) fails, the verdict at each qualifying start (the other rival at its carried
    optimum), under every reading."""
    out = {}
    for name in ("R0", "C1"):
        if carried[name].get("status") != "compared_if_verdicts_agree":
            continue
        for s, rec in (carried[name].get("laplace_per_qualifying_start") or {}).items():
            base = forms[name]
            th = np.array([rec["theta"][k] for k in base.obj.keys])
            f = Form(name, base.obj.point(th), base.numerics, base.obj, th, rec["laplace"], {})
            sc_f = {r: {k: score_units(lv, sc) for k, sc in settings.items()}
                    for r, lv in form_levels(f, READINGS).items()}
            other = "C1" if name == "R0" else "R0"
            sc_o = {r: {k: score_units(lv, sc) for k, sc in settings.items()}
                    for r, lv in form_levels(forms[other], READINGS).items()}
            outcomes = {r: decide({name: sc_f[r], other: sc_o[r]})["outcome"] for r in READINGS}
            out[f"{name}_start_{s}"] = verdict_across_readings(outcomes)
    return out


def sensitivities(forms: Mapping[str, Form], settings: Mapping[str, Scoring], gigl: bool) -> Dict[str, Any]:
    """sec. 7, every one computed once at the carried optimum, no refit, never used to choose."""
    out: Dict[str, Any] = {}
    prim = settings["primary"]
    for name in ("R0", "C1"):
        f = forms[name]
        rec: Dict[str, Any] = {}
        lv = decisive_many([f.point], "isothermal", f.numerics, stand_in=True)[0]
        rec[RESERVOIR_STAND_IN] = score_units({"levels": lv, "jac": None, "sigma_theta": None}, prim)
        lv = decisive_many([f.point], "isothermal", f.numerics, label_ph_held=True)[0]
        rec["reservoir_label_ph_held"] = score_units({"levels": lv, "jac": None, "sigma_theta": None}, prim)
        for label, change in (("K_0.045", {"k_l_per_mmol": 0.045}), ("K_0.075", {"k_l_per_mmol": 0.075}),
                              ("ea_decay_40", {"ea_decay_kj_mol": 40.0}), ("ea_decay_135", {"ea_decay_kj_mol": 135.0})):
            if not f.point.route:
                rec[label] = "inert: the route is off in this rival"
                continue
            p = Point(f.point.x_full, True, tuple(sorted({**dict(f.point.route_fields), **change}.items())))
            lv = decisive_many([p], "isothermal", f.numerics)[0]
            rec[label] = {"decisive": score_units({"levels": lv, "jac": None, "sigma_theta": None}, prim),
                          "b60_145C_rows": b60_rows_at(p, f.numerics, only_145=True)}
        if gigl:
            for ea in (112.0, 122.0):
                x = list(f.point.x_full)
                x[all_index("Ea_decay_thiol_sink")] = ea
                p = Point(tuple(x), f.point.route, f.point.route_fields)
                rec[f"gigl_relaxed_Ea_{ea:g}"] = {
                    r: score_units({"levels": decisive_many([p], r, f.numerics)[0], "jac": None, "sigma_theta": None},
                                   prim) for r in READINGS}
        out[name] = rec
    return out


def b60_rows_at(point: Point, numerics: Numerics, only_145: bool = False) -> Dict[str, float]:
    r, failed = evaluate(EvalSpec(numerics, FitSide.off(), False), point)
    FAILURES.extend(failed)
    out = {}
    for i, row in enumerate(B23.ACTIVE_FIT_ROWS):
        if only_145 and float(B23.SYSTEMS[row["system"]]["t_c"]) != 145.0:
            continue
        out[row["id"]] = float(r[i])
    return out


def reservoir_record(forms: Mapping[str, Form], tau: float) -> Dict[str, Any]:
    """sec. 5: R0's structural zeros (reported, never summed), C1's log densities on those rows (in-sample for Cao and
    Zhai, out-of-sample for Zhang 2025b, seen at start), and both rivals' on Cao Fig. 4 (in-sample for C1). Densities:
    per paper x species series, the row sigma with the integrated N(0, tau^2) level offset (declared)."""
    out: Dict[str, Any] = {"structural_zero_label": "R0 cannot represent the measured reservoir (EX-060)"}
    pots_needed = sorted(set(STRUCTURAL_ZERO_POTS) | {"cao2017_blank_fig4"})
    _, rows = load_table(RES_TABLE)
    for name in ("R0", "C1"):
        f = forms[name]
        jobs = [table_job(RES_TABLE, p, sorted({r.t_min for r in rows if r.pot == p})) for p in pots_needed]
        cfg = config_of(f.point, f.numerics, FitSide.off())
        failures: List[Dict[str, Any]] = []
        with config.applied(cfg):
            tasks, results = run_table_pots(f.point, jobs, f.numerics, cfg, failures)
            res = dict(zip(pots_needed, zip(tasks, results)))
            rec = {}
            for pot in pots_needed:
                task, got = res[pot]
                species = sorted({r.species for r in rows if r.pot == pot})
                for s in species:
                    rs = [r for r in rows if r.pot == pot and r.species == s]
                    if name == "R0" and s in STRUCTURAL_ZERO_POTS.get(pot, ()):
                        rec[f"{pot}|{s}"] = {"structural_zero": True, "n": len(rs)}
                        continue
                    if "failed" in got:
                        rec[f"{pot}|{s}"] = {"failed": got["failed"]}
                        continue
                    pred = np.array([probe.predict(r, task, got) for r in rs])
                    y = np.array([math.log10(r.value) for r in rs])
                    f_ = np.log10(np.maximum(pred, 1e-30))
                    sig = np.array([r.sigma for r in rs])
                    cov = np.diag(sig ** 2) + tau ** 2 * np.ones((len(rs), len(rs)))
                    rec[f"{pot}|{s}"] = {"log_density": log_normal(y, f_, cov), "n": len(rs),
                                         "rms_dex": float(np.sqrt(np.mean((f_ - y) ** 2)))}
        FAILURES.extend(failures)
        out[name] = rec
    return out


# ----------------------------------------------------------------------------------------------------------------
# Side effects (sec. 8)
# ----------------------------------------------------------------------------------------------------------------

E10_E12_TEST_CONTRASTS = ("F04#26", "F28#8", "F28#20", "F28#25", "F28#26")
DISULFIDE_TABLES = {"thiol_submodel": ("zhou2022_axarp_cys_ph7_120C", "wang2026b_xyl_cys_120c_ph55"),
                    "thiol_zhou2023_ph_series": ("zhou2023_arp_cys_120c_ph6", "zhou2023_arp_cys_120c_ph8")}
HOLDOUT_THIOL_LEDGER_PREFIXES = ("yiltirak2026", "hofmann1998_t1_")
HOLDOUT_THIOL_POTS = ("hofmann1998_cys_ribose_ph3", "hofmann1998_cys_ribose_ph7")
TAGS = {
    "de Roos 2005 Fig. 4 (heldout_deroos2005_fig4.yml, E10 TEST)": "spent in C1's rows (a declared constant source)",
    "BTBU (Cao 2017, Zhao 2017b fitted; Zhang 2025b chosen on)": "lost laboratory independence: F28#25 (E12), F23#5, "
        "F23#7, F24#4, F25#23 (E6 / E11), F27#5 (E18), heldout_zhao2017b_tab1.yml (E7)",
    "DFA (Schieberle 2000 chosen on)": "lost laboratory independence: F27#15 - F27#17 (E18, Grosch 1992)",
    "Reading (Whitfield 1999)": "pre-existing, R0's objective (B60's three Whitfield MFT rows); M1 adds none",
    "zhou_pH7_dimer_over_MFT, zhang_fig1_cys_dimer_over_MFT": "declared hold-out, in-sample since B60",
}
COMPENSATING_ERRORS = (
    ("Kang fed-TTCA rows (kang_100C_MFT, kang_120C_MFT, FFT, MFT_over_FFT)", "the route retires r_pent_dpo",
     "the instant cysteine-catalysed enolisation standing in for the route; k_nf_mft ~x21 high (Cerny pot D, M2)"),
    ("zhai_MFT_fold_120_over_100, zhai_FFT_fold_120_over_100", "the supply's barrier moves to ea_form",
     "one barrier for every route (E11, M2); the hexose transfer of 100.9 kJ/mol"),
    ("the 145 C fed + H2S rows (fed_ribose_h2s_*, fed_tdp_h2s_FFT)", "k_mft_decay moves under the prior",
     "the frozen fast sink (CF-107); the H2S headspace bias (CF-065); the isothermal reading of the Hofmann pots (CR-03)"),
    ("free MFT in oxidant-bearing pots (Zhang 2024 cystine arm; Zhou 2023)", "the dimer ceiling lets k_dimer_mft rise",
     "the declared ambient OX (1.0 mM) and the irreversible dimer (the C3 follow-up)"),
    ("cerny_ternary_thiamine_share, fed_thiamine_MFT; F27#15", "cysteine sheltered in TTCA",
     "the thiamine route (k_thi_hmp, k_hmp_mft) fitted under instant supply"),
    ("zhou_pH7_FUR_arp_alone, hofmann_ribose_FUR_insitu", "the pentose's fate changes",
     "k_pent_caramel and k_fur_decay fitted against an instant pentose sink"),
    ("the Zhou 2023 final-pH rows", "TTCA formation and CAMA titration move the charge balance",
     "the B51 balance's acid allocation (CF-084, M1b)"),
    ("Kumazawa FFT retention", "none expected", "-- (a loss here is a build error until shown otherwise)"),
    ("Yiltirak 100 - 130 C levels", "the sink's T-dependence moves", "laboratory and work-up offsets (CF-068, R4)"),
    ("any row, by up to the numerics row's size", "inline vs fixed point", "none: numerics"),
)


def table_raw(point: Point, numerics: Numerics, table: str, pot_ids: Sequence[str], species: Optional[Sequence[str]] = None
              ) -> Dict[str, Dict[str, Any]]:
    """Raw log10 residuals (predicted / measured) of a fit table's rows on some pots: a side-effect record."""
    pots, rows = load_table(table)
    rs = [r for r in rows if r.pot in pot_ids and (species is None or r.species in species)]
    want = sorted({r.pot for r in rs})
    jobs = [table_job_from(pots, rows, table, p) for p in want]
    jobs = [j for j in jobs if j is not None]
    cfg = config_of(point, numerics, FitSide.off())
    failures: List[Dict[str, Any]] = []
    out: Dict[str, Dict[str, Any]] = {}
    with config.applied(cfg):
        tasks, results = run_table_pots(point, jobs, numerics, cfg, failures)
        by = {t["pot"]: (t, res) for t, res in zip(tasks, results)}
        for r in rs:
            key = f"{r.pot}|{r.species}|{r.t_min:g}"
            if r.pot not in by:
                out[key] = {"note": "pot not runnable"}
                continue
            t, res = by[r.pot]
            if "failed" in res:
                out[key] = {"failed": res["failed"]}
                continue
            p = probe.predict(r, t, res)
            out[key] = {"predicted": p, "measured": r.value, "unit": r.unit, "sigma": r.sigma,
                        "raw_residual": fit_rows.residual(r, p) if p > 0 or not r.log_scale else None,
                        "label_ph": t["label_ph"]}
    FAILURES.extend(failures)
    return out


def table_job_from(pots, rows, table, pot) -> Optional[Dict[str, Any]]:
    jobs = [j for j in probe.plan({pot: pots[pot]}, rows) if j["pot"] == pot]
    return None if not jobs or "skip" in jobs[0] else jobs[0]


def disulfide_ratios(raw: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    """MFTD / MFT and FFTD / FFT at matching times (predicted vs measured, log10), noise the rows' errors in
    quadrature: per source, never summed."""
    out = {}
    for key, rec in raw.items():
        pot, s, t = key.split("|")
        if s not in ("MFTD", "FFTD") or "predicted" not in rec:
            continue
        base = raw.get(f"{pot}|{s[:-1]}|{t}")
        if not base or "predicted" not in base:
            continue
        out[f"{pot}|{s}/{s[:-1]}|{t}"] = {
            "log10_predicted": math.log10(max(rec["predicted"], 1e-30) / max(base["predicted"], 1e-30)),
            "log10_measured": math.log10(rec["measured"] / base["measured"]),
            "sigma_dex": math.hypot(rec["sigma"], base["sigma"])}
    return out


def run_side_effects(forms: Mapping[str, Form], commit: str) -> Dict[str, Any]:
    import score_contrast_benchmark as SCB

    state = _read(SIDE_OUT)
    state.update({"artifact": SIDE_OUT.stem, "prereg": data_paths.rel(PREREG), "freeze_commit": commit,
                  "never_used_to_choose": True, "tags": TAGS,
                  "candidate_compensating_errors": [dict(zip(("expected_loss", "why", "candidate"), c))
                                                    for c in COMPENSATING_ERRORS]})
    names = [n for n in ("R0_shipped_a", "R0_shipped_b", "R0", "C1") if n in forms]
    rows, _, _ = route_block()
    route_pots = sorted({r.pot for r in rows})
    for n in names:
        f = forms[n]
        rec = state.setdefault("forms", {}).setdefault(n, {})
        if "b60_rows" not in rec:
            rec["b60_rows"] = b60_rows_at(f.point, f.numerics)
            _write(SIDE_OUT, state)
        if "route_rows" not in rec:
            rec["route_rows"] = table_raw(f.point, f.numerics, RES_TABLE, route_pots)
            if not f.point.route:   # sec. 5: a structural zero is reported, never summed
                for k, v in rec["route_rows"].items():
                    if k.split("|")[1] in ("TTCA", "CAMA"):
                        v["structural_zero"] = True
                        v["raw_residual"] = None
            by_ph: Dict[str, List[float]] = defaultdict(list)
            for v in rec["route_rows"].values():
                if v.get("raw_residual") is not None:
                    by_ph[f"{v['label_ph']:g}"].append(v["raw_residual"])
            rec["route_residuals_by_ph"] = {k: {"n": len(v), "mean": float(np.mean(v)),
                                                "rms": float(np.sqrt(np.mean(np.square(v))))}
                                            for k, v in sorted(by_ph.items())}
            _write(SIDE_OUT, state)
        if "declared_holdouts" not in rec:
            pots, trows = load_table("thiol_submodel")
            hp = sorted({r.pot for r in trows if r.ledger_id.startswith(HOLDOUT_THIOL_LEDGER_PREFIXES)
                         or r.pot in HOLDOUT_THIOL_POTS})
            rec["declared_holdouts"] = {"thiol_submodel": table_raw(f.point, f.numerics, "thiol_submodel", hp),
                                        "zhou2023_ph6_ph8": table_raw(f.point, f.numerics, "thiol_zhou2023_ph_series",
                                                                      DISULFIDE_TABLES["thiol_zhou2023_ph_series"]),
                                        "route_zhai_140C": table_raw(f.point, f.numerics, RES_TABLE,
                                                                     [p for p in ROUTE_OUT_POTS if "140c" in p])}
            _write(SIDE_OUT, state)
        if "disulfide_record" not in rec:
            raw = {}
            for table, pots_ in DISULFIDE_TABLES.items():
                raw.update(table_raw(f.point, f.numerics, table, pots_, ("MFT", "MFTD", "FFT", "FFTD")))
            rec["disulfide_record"] = {
                "ratios": disulfide_ratios(raw),
                "b60_rows": {k: v for k, v in rec["b60_rows"].items()
                             if k in ("zhou_pH7_dimer_over_MFT", "zhang_fig1_cys_dimer_over_MFT",
                                      "zhang_fig1_gcys_dimer_over_MFT")},
                "zhou2022_note": "not independent of the in-sample Zhou 2023 pH 7 pot it replicates (review M3)",
                "not_built": "Whitfield 1999's two disulfide arms and Zhang 2024 Fig. 1's free-MFT cystine / cysteine "
                             "ratio: no keyed fit-table rows on disk for them (B60's systems carry no such row); "
                             "registered for the C3 follow-up"}
            _write(SIDE_OUT, state)
        if n != "R0_shipped_a" and "contrast_benchmark" not in rec:
            cfg = config_of(f.point, f.numerics, FitSide.off())
            x = f.point.x()
            tag = f"m1_{n.lower()}"
            payload = SCB.score_named(tag, cfg, B23.build_parameters(x), B23.unpack(x)[3], f"M1 {n}", GENERATED_BY,
                                      inputs=[carried_path("R0"), carried_path("C1")])
            rec["contrast_benchmark"] = {"artifact": f"kinetic_core_contrast_benchmark_score_{tag}",
                                         "readings": payload.get("readings")}
            _write(SIDE_OUT, state)
    for n in ("R0", "C1"):
        if n in forms and "G1" not in state.setdefault("forms", {}).get(n, {}):
            state["forms"][n]["G1"] = _read(carried_path(n)).get("G1")
    if "R0_shipped_a" in state.get("forms", {}) and "R0_shipped_b" in state["forms"]:
        a, b = state["forms"]["R0_shipped_a"], state["forms"]["R0_shipped_b"]
        state["numerics_row"] = {
            "b60_rows_inline_minus_fixed_point": {k: b["b60_rows"][k] - a["b60_rows"][k] for k in a["b60_rows"]},
            "route_rows_max_abs_dex": _max_abs_diff(a["route_rows"], b["route_rows"]),
            "contrast_benchmark_before_record": "kinetic_core_contrast_benchmark_score_v1 (B60, fixed point; frozen)"}
    items = _read(V / "kinetic_core_contrast_benchmark_items.json")
    atlas = {str(c.get("contrast", "")).split("@")[0]: c for c in (items.get("q1_atlas") or {}).get("items", [])}
    state["e10_e12_test_contrasts"] = {cid: atlas.get(cid, "not in the benchmark's atlas record: the scorer cannot "
                                                           "answer it") for cid in E10_E12_TEST_CONTRASTS}
    state["failed_pot_runs"] = len(FAILURES)
    state["provenance"] = _provenance(SIDE_OUT.stem, [COMPARISON_OUT])
    state["complete"] = True
    _write(SIDE_OUT, state)
    return state


def _max_abs_diff(a: Mapping[str, Any], b: Mapping[str, Any]) -> Optional[float]:
    d = [abs(a[k]["raw_residual"] - b[k]["raw_residual"]) for k in a
         if k in b and a[k].get("raw_residual") is not None and b[k].get("raw_residual") is not None]
    return max(d) if d else None


# ----------------------------------------------------------------------------------------------------------------
# --checks: V1, V2, V4 - V8 (sec. 4.4)
# ----------------------------------------------------------------------------------------------------------------

def machine_state() -> Dict[str, Any]:
    load = os.getloadavg()
    return {"load_average_1_5_15": [round(x, 2) for x in load], "cpu_count": os.cpu_count(),
            "free_memory_bytes": free_memory_bytes(), "swap": swap_usage(), "workers": parallel.workers(),
            "time": time.strftime("%Y-%m-%d %H:%M:%S")}


def free_memory_bytes() -> Optional[int]:
    """vm_stat: (free + inactive + speculative + purgeable) pages x page size; None off macOS."""
    try:
        txt = subprocess.run(["vm_stat"], capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    m = re.search(r"page size of (\d+) bytes", txt)
    page = int(m.group(1)) if m else 4096
    pages = 0
    for name in ("Pages free", "Pages inactive", "Pages speculative", "Pages purgeable"):
        mm = re.search(rf"{name}:\s+(\d+)", txt)
        pages += int(mm.group(1)) if mm else 0
    return pages * page


def swap_usage() -> Optional[str]:
    try:
        return subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def v1() -> Dict[str, Any]:
    t0 = time.time()
    r, failed = evaluate(EvalSpec(fit_fixed_point(), FitSide.off(), False), b60_point())
    cost = 0.5 * float(r @ r)
    rel = abs(cost - B60_RECORD_COST) / B60_RECORD_COST
    return {"what": "the route key off: B60's 52-row cost on the node fixed point at B60's quick setting",
            "cost": cost, "record_cost": B60_RECORD_COST, "relative_difference": rel, "tolerance": V1_REL_TOL,
            "passed": bool(rel <= V1_REL_TOL and not failed), "failed_pot_runs": len(failed),
            "seconds": round(time.time() - t0, 1),
            "tests": "every existing test: run separately (pytest tests/unit); recorded in the build report"}


def v2() -> Dict[str, Any]:
    """C1's switch on at B60's vector: the Cao residuals through M1's path (B23.build_parameters of the vector, the
    route block's jobs) against a reference recomputed on main's own probe path (build_submodel_probe.run_pot with
    the shipped core parameters, no override) with M1's switch set (route on, inline; B59a and B61 off)."""
    t0 = time.time()
    point = b60_point(True)
    jobs = [table_job(RES_TABLE, p, sorted({r.t_min for r in load_table(RES_TABLE)[1] if r.pot == p}))
            for p in CAO_POTS]
    _, rows = load_table(RES_TABLE)
    out: Dict[str, Any] = {"what": v2.__doc__.split("\n")[0], "pots": list(CAO_POTS)}
    for label, numerics in (("score_rtol", score_inline()), ("fit_rtol", fit_inline())):
        cfg = config_of(point, numerics, FitSide.off())
        failures: List[Dict[str, Any]] = []
        with config.applied(cfg):
            tasks, results = run_table_pots(point, jobs, numerics, cfg, failures)
            ref = [probe.run_pot(dict(j)) for j in jobs]   # main's path: shipped core parameters, rtol 1e-8
        diffs, resid = [], {}
        for job, task, got, rr in zip(jobs, tasks, results, ref):
            for r in rows:
                if r.pot != job["pot"]:
                    continue
                a = fit_rows.residual(r, probe.predict(r, task, got))
                b = fit_rows.residual(r, probe.predict(r, job, rr))
                resid[f"{r.pot}|{r.species}|{r.t_min:g}"] = {"m1": a, "reference": b}
                diffs.append(abs(a - b))
        out[label] = {"max_abs_diff_dex": max(diffs), "failed_pot_runs": len(failures), "residuals": resid}
    out["tolerance_dex_at_score_rtol"] = V2_TOL_DEX
    out["passed"] = bool(out["score_rtol"]["max_abs_diff_dex"] <= V2_TOL_DEX and not out["score_rtol"]["failed_pot_runs"])
    out["seconds"] = round(time.time() - t0, 1)
    return out


def v4_measure(rival: str) -> Dict[str, Any]:
    """One objective evaluation timed on inline at fit tolerance (run in its own process for a clean peak RSS)."""
    fit_side = FitSide.on(1.0)
    noise = start0_noise(fit_side, fit_inline()) if RIVALS[rival].route else {}
    obj = Objective(RIVALS[rival], fit_inline(), fit_side, noise)
    th = obj.theta_of(b60_point(RIVALS[rival].route))
    obj.residuals(th)      # warm (imports, caches)
    t0 = time.time()
    r = obj.residuals(th)
    sec = time.time() - t0
    t1 = time.time()
    obj.residuals_many([th] * max(1, parallel.workers()))
    sec_many = time.time() - t1
    self_ = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    scale = 1 if sys.platform == "darwin" else 1024
    return {"rival": rival, "seconds_one_evaluation": round(sec, 2),
            "seconds_workers_points_in_parallel": round(sec_many, 2), "cost": 0.5 * float(r @ r),
            "peak_rss_parent_bytes": int(self_ * scale), "peak_rss_largest_child_bytes": int(child * scale),
            "workers": parallel.workers(), "failed_pot_runs": len(FAILURES)}


def v4() -> Dict[str, Any]:
    out = {"machine": machine_state()}
    for rival in ("R0", "C1"):
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--v4-measure", rival], cwd=str(ROOT),
                              capture_output=True, text=True)
        line = [ln for ln in proc.stdout.splitlines() if ln.startswith("V4MEASURE ")]
        if proc.returncode != 0 or not line:
            out[rival] = {"error": proc.stderr[-2000:]}
            continue
        m = json.loads(line[-1][len("V4MEASURE "):])
        w = max(1, int(m["workers"]))
        m["peak_rss_estimate_one_fit_at_4_workers_bytes"] = (m["peak_rss_parent_bytes"]
                                                              + FIT_LANE_WORKERS * m["peak_rss_largest_child_bytes"])
        m["measured_at_workers"] = w
        out[rival] = m
    out["peak_rss_one_fit_bytes"] = max((out[r].get("peak_rss_estimate_one_fit_at_4_workers_bytes", 0)
                                         for r in ("R0", "C1")), default=0)
    out["note"] = ("the estimate at 4 workers is parent + 4 x the largest child (a forked pool); measured at "
                   f"MAILLARD_WORKERS={parallel.workers()} under the machine state above. Re-run "
                   "'--checks --only V4' at MAILLARD_WORKERS=4 before --all (prereg sec. 12: with the other fit running)")
    out["passed"] = all("error" not in out[r] for r in ("R0", "C1"))
    return out


def panel_max_mft(point: Point) -> Dict[str, Any]:
    """The objective's pots' largest MFT over time at a point (B60's systems and the route block, dense grids)."""
    numerics = score_inline()
    cfg = config_of(point, numerics, FitSide.off())
    x = point.x()
    params, drift = B23.build_parameters(x), B23.unpack(x)[3]
    best = {"mmol_l": 0.0, "pot": None}
    with config.applied(cfg):
        items = []
        for name, spec in B23.SYSTEMS.items():
            minutes = float(spec["minutes"])
            dense = dict(spec, grid_points=0)
            items.append((name, dense, params, drift, numerics, minutes))
        runs = parallel.map_ordered(_dense_system, items, config=cfg)
        for (name, *_), m in zip(items, runs):
            if m is not None and m > best["mmol_l"]:
                best = {"mmol_l": float(m), "pot": f"b60:{name}"}
        rows, jobs, _ = route_block()
        dense_jobs = []
        for p, j in jobs.items():
            j = dict(j)
            j["times"] = sorted({*j["times"], *np.linspace(0.0, max(j["times"]), 121).tolist()})
            dense_jobs.append(j)
        failures: List[Dict[str, Any]] = []
        tasks, results = run_table_pots(point, dense_jobs, numerics, cfg, failures)
        for t, res in zip(tasks, results):
            if "failed" in res:
                continue
            m = max(s["MFT"] for s in res["states"])
            if m > best["mmol_l"]:
                best = {"mmol_l": float(m), "pot": f"route:{t['pot']}"}
    return best


def _dense_system(item) -> Optional[float]:
    name, spec, params, drift, numerics, minutes = item
    grid = np.unique(np.concatenate([[0.0], np.geomspace(minutes * 1e-6, minutes, 121), np.linspace(0.0, minutes, 121)]))
    try:
        run = with_alarm(POT_SECONDS, lambda: B23.integrate_sulfur(
            params, float(spec["t_c"]) + CELSIUS, spec["initial"], grid, ph=float(spec["ph"]),
            buffer_spec=spec.get("buffer"), ph_drift=drift, rtol=numerics.rtol, atol=ATOL))
        return float(np.max(run.series("MFT")))
    except Exception:  # noqa: BLE001
        return None


def v5(noise: Mapping[str, float]) -> Dict[str, Any]:
    t0 = time.time()
    pm = panel_max_mft(b60_point(False))
    cap = mft_cap_from_panel_max(pm["mmol_l"])
    fs = FitSide.on(cap)
    out: Dict[str, Any] = {"panel_max_mft_at_B60": pm, "mft_cap_mmol_l": cap,
                           "ceiling_log10_145C_by_rival": {r: dimer_ceiling(fs.ceiling, RIVALS[r].route) for r in RIVALS},
                           "ceiling_at_145C_alone": fs.ceiling.log10_ceiling_at(145.0), "evaluations": {}}
    obj = Objective(RIVALS["C1"], fit_inline(), fs, noise)
    th0 = obj.theta_of(b60_point(True))
    i = obj.keys.index("k_dimer_mft")
    ok = True
    for label, value in (("log10_k_dimer_mft_4", V5_LOG10_K_DIMER), ("at_the_ceiling", float(obj.upper[i]))):
        th = th0.copy()
        th[i] = value
        t1 = time.time()
        r, failed = evaluate(obj.spec, obj.point(th))
        sec = time.time() - t1
        good = not failed and bool(np.all(np.isfinite(r))) and not np.all(r[:obj.n_b60] == 25.0)
        ok &= good
        out["evaluations"][label] = {"log10_k_dimer_mft_145C": value, "seconds": round(sec, 2), "failed": failed,
                                     "cost": 0.5 * float(r @ r), "integrated": good}
    out["integrated_without_failure"] = bool(ok)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def v6(cap: float) -> Dict[str, Any]:
    fs = FitSide.on(cap)
    obj = Objective(RIVALS["R0"], fit_inline(), fs)
    th = obj.theta_of(b60_point(False))
    r, failed = evaluate(obj.spec, obj.point(th))
    return {"what": "R0 at B60's vector with the ceiling and the prior, inline, fit tolerance: R0's start-0 reference",
            "cost": 0.5 * float(r @ r), "cost_by_block": obj.blocks(r), "prior_residual": float(r[-1]),
            "prior_term_nats": 0.5 * float(r[-1]) ** 2,
            "log10_lambda_mft_100C_at_B60": prior_log10_rate(b60_point().x(), fs.prior),
            "prior_centre": fs.prior.centre_log10_per_min, "failed_pot_runs": len(failed)}


def b60_predicted(point: Point, numerics: Numerics) -> Dict[str, float]:
    cfg = config_of(point, numerics, FitSide.off())
    failures: List[Dict[str, Any]] = []
    with config.applied(cfg), b60_runner(numerics, failures):
        x = point.x()
        predicted, _ = B23.observables(B23.build_parameters(x), numerics.quick, B23.unpack(x)[3])
    return predicted


def v7_diff(point: Point) -> Dict[str, Any]:
    """Inline (scoring rtol) against fixed_point_nodes at 81 nodes, rtol 1e-8: decisive pots, Cao 2017 Figs. 4 - 5
    and B60's 52 rows."""
    a_lv = decisive_levels(point, "isothermal", score_inline())[0]
    b_lv = decisive_levels(point, "isothermal", score_fixed_point_81())[0]
    pa, _ = decisive_tolerance(a_lv, None, None)
    pb, _ = decisive_tolerance(b_lv, None, None)
    us = units(a_lv, None, Scoring(tau=0.0))
    tol = []
    for u in us.values():   # 0.1 x the per-ratio decision sd (D1), 0.1 x the row sigma (reservoir; tau 0 here)
        tol += list(V7_DECISION_FRACTION * np.sqrt(np.diag(u["Sigma"])))
    labels = [lbl for u in us.values() for lbl in u["labels"]]
    lab_all = point_labels()
    for s_ in D1_SPECIES:   # the D1 level points too (~0.004 dex each)
        for t in D1_TIMES:
            i = lab_all.index(f"D1|{s_}|{t:g}")
            pa = np.append(pa, a_lv[i])
            pb = np.append(pb, b_lv[i])
            tol.append(V7_DECISION_FRACTION * math.sqrt(d1_point_variance(s_, t, "decision")))
            labels.append(f"{s_} level {t:g} min")
    diff = np.abs(pa - pb)
    tol = np.array(tol)
    cao = {}
    _, rows = load_table(RES_TABLE)
    for numerics in (score_inline(), score_fixed_point_81()):
        cao[numerics.name] = table_raw(point, numerics, RES_TABLE, CAO_POTS)
    cao_diff = {k: abs(cao[score_inline().name][k]["raw_residual"] - cao[score_fixed_point_81().name][k]["raw_residual"])
                for k in cao[score_inline().name] if cao[score_inline().name][k].get("raw_residual") is not None}
    cao_sigma = {f"{r.pot}|{r.species}|{r.t_min:g}": r.sigma for r in rows if r.pot in CAO_POTS}
    pa52, pb52 = b60_predicted(point, score_inline()), b60_predicted(point, score_fixed_point_81())
    kinds = {r["id"]: r["kind"] for r in B23.ACTIVE_FIT_ROWS}
    rows52 = {k: (pa52[k] - pb52[k]) if kinds[k] == "ph_endpoint" else
              (math.log10(max(pa52[k], 1e-30) / max(pb52[k], 1e-30)) if pa52[k] > 0 and pb52[k] > 0 else None)
              for k in pa52}
    return {"decisive_points": {lbl: {"diff_dex": float(d), "tolerance_dex": float(t)}
                                for lbl, d, t in zip(labels, diff, tol)},
            "decisive_max_diff_over_tolerance": float(np.max(diff / tol)),
            "passed": bool(np.all(diff <= tol)),
            "cao_figs_4_5_abs_diff_dex": cao_diff,
            "cao_max_diff_over_0.1_row_sigma": max(cao_diff[k] / (0.1 * cao_sigma[k]) for k in cao_diff) if cao_diff else None,
            "b60_52_rows_inline_minus_fixed_point": rows52,
            "b60_52_rows_max_abs": max(abs(v) for v in rows52.values() if v is not None)}


def v7() -> Dict[str, Any]:
    t0 = time.time()
    out = {"what": "inline (rtol 1e-8) vs fixed_point_nodes at 81 nodes, rtol 1e-8; tolerance 0.1 sigma_decision per "
                   "decisive point (D1 ratios; reservoir row sigma); Cao and the 52 rows printed, not gated",
           "at_B60": v7_diff(b60_point(False)), "at_C1_start": v7_diff(b60_point(True))}
    out["passed"] = bool(out["at_B60"]["passed"] and out["at_C1_start"]["passed"])
    out["seconds"] = round(time.time() - t0, 1)
    return out


def v8(noise: Mapping[str, float], cap: float) -> Dict[str, Any]:
    """Smoothness on five coordinates at C1's start: central differences at two step sizes, inline (fit tolerance)
    against the fixed point (quick setting). Fails if an inline difference changes sign between the step sizes where
    the fixed point's does not, or differs from the fixed point's by more than a factor 2."""
    t0 = time.time()
    fs = FitSide.on(cap)
    out: Dict[str, Any] = {"coordinates": {}}
    passed = True
    for numerics in (fit_inline(), fit_fixed_point()):
        obj = Objective(RIVALS["C1"], numerics, fs, noise)
        th0 = obj.theta_of(b60_point(True))
        pts, plan_ = [th0], []
        for k in V8_COORDINATES:
            i = obj.keys.index(k)
            for h in V8_STEPS[obj.kinds[i]]:
                for sgn in (+1, -1):
                    th = th0.copy()
                    th[i] += sgn * h
                    pts.append(th)
                    plan_.append((k, h, sgn))
        rs = obj.residuals_many(pts)
        costs = [0.5 * float(r @ r) for r in rs]
        c0 = costs[0]
        by = defaultdict(dict)
        for (k, h, sgn), c in zip(plan_, costs[1:]):
            by[k][(h, sgn)] = c
        for k in V8_COORDINATES:
            steps = sorted({h for h, _ in by[k]})
            fd = {h: (by[k][(h, 1)] - by[k][(h, -1)]) / (2 * h) for h in steps}
            out["coordinates"].setdefault(k, {})[numerics.name] = {
                "centre_cost": c0, "costs": {f"{sgn * h:+g}": c for (h, sgn), c in by[k].items()},
                "central_difference": {f"{h:g}": v for h, v in fd.items()}}
    for k, rec in out["coordinates"].items():
        a, b = rec[fit_inline().name]["central_difference"], rec[fit_fixed_point().name]["central_difference"]
        hs = sorted(a, key=float)
        sign_flip_inline = np.sign(a[hs[0]]) != np.sign(a[hs[1]])
        sign_flip_fixed = np.sign(b[hs[0]]) != np.sign(b[hs[1]])
        ratios = {h: (a[h] / b[h] if b[h] != 0 else (1.0 if a[h] == 0 else float("inf"))) for h in hs}
        bad_ratio = any(not (0.5 <= r <= 2.0) for r in ratios.values())
        fail = bool((sign_flip_inline and not sign_flip_fixed) or bad_ratio)
        rec.update({"inline_over_fixed": ratios, "sign_change_inline": bool(sign_flip_inline),
                    "sign_change_fixed_point": bool(sign_flip_fixed), "failed": fail})
        passed &= not fail
    out["passed"] = bool(passed)
    out["consequence"] = ("fits stay on inline" if passed else
                          "every fit runs on fixed_point_nodes at the quick setting (B60's numerics); the scoring "
                          "stays as declared; V7's diff is printed at the carried optima (v3.3)")
    out["seconds"] = round(time.time() - t0, 1)
    return out


def row_checks() -> Dict[str, Any]:
    """sec. 4.1's three checks on the rows: duplicates, declared hold-outs inside B60's objective, block sizes."""
    rows, jobs, pots = route_block()
    return {"route_block_rows": len(rows), "route_block_pots": len(jobs),
            "prereg_states": "296 rows in 34 pots before the duplicates check",
            "duplicates_dropped": {p: pots[p].get("duplicate_note") for p in ROUTE_DUPLICATE_POTS},
            "declared_holdouts_in_b60_objective": [i for i in ("zhou_pH7_dimer_over_MFT", "zhang_fig1_cys_dimer_over_MFT")
                                                    if i in b60_row_ids()],
            "b60_rows": len(B23.ACTIVE_FIT_ROWS),
            "route_rows_by_label_ph": {f"{k:g}": v for k, v in sorted(Counter_(float(jobs[r.pot]["label_ph"])
                                                                                for r in rows).items())}}


V_CHECKS = ("ROWS", "V1", "V2", "V4", "V5", "V6", "V7", "V8")


def run_checks(only: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    state = _read(CHECKS_OUT)
    state.update({"artifact": CHECKS_OUT.stem, "prereg": data_paths.rel(PREREG),
                  "what": "numerics checks (sec. 4.4): predictions against predictions, never against measurements"})
    checks = state.setdefault("checks", {})
    want = list(only) if only else [c for c in V_CHECKS if c not in checks]

    def save():
        state["provenance"] = _provenance(CHECKS_OUT.stem)
        _write(CHECKS_OUT, state)

    noise = state.get("c1_start0_series_noise")
    for name in want:
        print(f"== {name} ({time.strftime('%H:%M:%S')}; load {os.getloadavg()[0]:.1f})", flush=True)
        del FAILURES[:]
        if name in ("V5", "V8", "V4") and noise is None:
            noise = start0_noise(FitSide.on(1.0), fit_inline())
            state["c1_start0_series_noise"] = noise
        cap = float((checks.get("V5") or {}).get("mft_cap_mmol_l", 1.0))
        rec = {"ROWS": row_checks, "V1": v1, "V2": v2, "V4": v4, "V5": lambda: v5(noise), "V6": lambda: v6(cap),
               "V7": v7, "V8": lambda: v8(noise, cap)}[name]()
        rec["machine"] = machine_state()
        checks[name] = rec
        save()
        print(f"   {name}: " + json.dumps({k: v for k, v in rec.items() if k in (
            "passed", "cost", "relative_difference", "max_abs_diff_dex", "integrated_without_failure", "mft_cap_mmol_l",
            "prior_term_nats", "seconds", "peak_rss_one_fit_bytes", "route_block_rows")}, default=str), flush=True)
    return state


# ----------------------------------------------------------------------------------------------------------------
# --probe (sec. 3)
# ----------------------------------------------------------------------------------------------------------------

PROBE_KEYS = ("k_ttca_cys", "ea_form_kj_mol", "k_l_per_mmol", "k_form_ref_per_min", "ea_decay_kj_mol",
              "Ea_decay_thiol_sink", "lambda_mft_100C")


def lambda_max_log10(ea_kj_mol: float) -> float:
    """sec. 3 (re-review MINOR 3): the fit's reachable maximum of log10 lambda(100 C) at a draw's Ea: k(145 C) <= 10^0.5."""
    t, t_ref = 100.0 + CELSIUS, 145.0 + CELSIUS
    return float(ps.FITTED_SULFUR_BOUNDS_LOG10K[PRIOR_KEY][1]) - ea_kj_mol / (ps.R_KJ * math.log(10.0)) * (1.0 / t - 1.0 / t_ref)


def probe_point(u: Sequence[float]) -> Tuple[Point, Dict[str, float]]:
    """A unit-cube point -> C1 at that draw (every other constant at B60)."""
    base = b60_point(True)
    x = list(base.x_full)
    s = float(x[all_index("k_ttca_cys")])
    vals = {"k_ttca_cys": s - 3.0 + 6.0 * u[0], "ea_form_kj_mol": EA_FORM_BAND[0] + u[1] * (EA_FORM_BAND[1] - EA_FORM_BAND[0]),
            "k_l_per_mmol": K_TTCA_PROBE_BAND[0] + u[2] * (K_TTCA_PROBE_BAND[1] - K_TTCA_PROBE_BAND[0]),
            "k_form_ref_per_min": float(shipped_value("sulfur.cys_amadori_route").k_form_ref_per_min)
            * 10.0 ** (-K_FORM_PROBE_HALF_DEX + 2 * K_FORM_PROBE_HALF_DEX * u[3]),
            "ea_decay_kj_mol": EA_DECAY_PROBE_BAND[0] + u[4] * (EA_DECAY_PROBE_BAND[1] - EA_DECAY_PROBE_BAND[0]),
            "Ea_decay_thiol_sink": 7.0 + u[5] * (102.0 - 7.0)}
    lo = ps.MftDecayPrior().centre_log10_per_min - 2.0
    hi = lambda_max_log10(vals["Ea_decay_thiol_sink"])
    vals["lambda_mft_100C"] = lo + u[6] * (hi - lo)
    x[all_index("k_ttca_cys")] = vals["k_ttca_cys"]
    x[all_index("Ea_decay_thiol_sink")] = vals["Ea_decay_thiol_sink"]
    t, t_ref = 100.0 + CELSIUS, 145.0 + CELSIUS
    x[all_index(PRIOR_KEY)] = vals["lambda_mft_100C"] + vals["Ea_decay_thiol_sink"] / (ps.R_KJ * math.log(10.0)) * (
        1.0 / t - 1.0 / t_ref)
    fields = {"ea_form_kj_mol": vals["ea_form_kj_mol"], "k_l_per_mmol": vals["k_l_per_mmol"],
              "k_form_ref_per_min": vals["k_form_ref_per_min"], "ea_decay_kj_mol": vals["ea_decay_kj_mol"]}
    return Point(tuple(x), True, tuple(sorted(fields.items()))), vals


def probe_value(levels: np.ndarray, tau: float) -> Tuple[float, Dict[str, float]]:
    """sec. 3: the plug-in log density on the decisive units (D1 and the reservoir), floors applied."""
    us = units(levels, None, Scoring(tau))
    per = {n: log_normal(u["y"], u["f"], u["Sigma"]) for n, u in us.items()}
    return float(sum(per.values())), per


def probe_whitened(levels: np.ndarray, tau: float) -> np.ndarray:
    us = units(levels, None, Scoring(tau))
    out = []
    for u in us.values():
        r = u["f"] - u["y"]
        if not np.all(np.isfinite(r)):
            return np.full(sum(len(v["y"]) for v in us.values()), 25.0)
        out.append(np.linalg.solve(np.linalg.cholesky(u["Sigma"]), r))
    return np.concatenate(out)


def s1_s2(point: Point, numerics: Numerics) -> Dict[str, Any]:
    """v2's S1 (the Schieberle rise) and S2 (the Cao reservoir) at the best draw: descriptions, no pass mark."""
    lv = decisive_many([point], "isothermal", numerics)[0]
    us = units(lv, None, Scoring(0.0))
    s1 = {n: {"predicted_dex": u["f"].tolist(), "measured_dex": u["y"].tolist(), "labels": u["labels"]}
          for n, u in us.items() if u["block"] == "D1"}
    cao = table_raw(point, numerics, RES_TABLE, CAO_POTS)
    return {"S1_schieberle_rise": s1, "S2_cao_reservoir_raw_residual_dex": {k: v.get("raw_residual") for k, v in cao.items()}}


def run_probe(tau_cli: Optional[float], confirm: bool) -> Dict[str, Any]:
    tau = read_tau(tau_cli)
    if tau is None:
        raise SystemExit("--probe waits for H0's tau_ABS(all) (prereg secs. 6, 10): pass --tau-abs or score EX-059 first")
    numerics = fit_inline()
    state = _read(PROBE_OUT)
    del FAILURES[:]
    if confirm:
        if not state.get("kill_would_apply"):
            raise SystemExit("--confirm: the probe's rule did not kill C1; no confirmation fit is needed")
        rec = _read(fit_path("C1", 0))
        if not rec.get("complete") or int(rec.get("failed_pot_runs", 0)) > 0:
            raise SystemExit("--confirm needs C1's start-0 fit, complete and with no failed pot runs")
        obj = Objective(RIVALS["C1"], numerics, FitSide.off(), rec.get("series_noise") or {})
        p = obj.point(np.array([rec["theta"][k] for k in obj.keys]))
        v, per = probe_value(decisive_many([p], "isothermal", numerics)[0], tau)
        state["confirmation"] = {"c1_start0_value": v, "per_unit": per, "r0_shipped_b": state["r0_shipped_b"]["value"],
                                 "killed": bool(v <= state["r0_shipped_b"]["value"])}
        _write(PROBE_OUT, state)
        artifact_io.write_artifact(state, PROBE_OUT, render=render_probe)
        return state
    from scipy.stats import qmc

    sob = qmc.Sobol(d=len(PROBE_KEYS), scramble=True, seed=PROBE_SEED)
    u = sob.random(PROBE_DRAWS)
    t0 = time.time()
    pts = [probe_point(row) for row in u]
    levels = decisive_many([p for p, _ in pts], "isothermal", numerics)
    draws = []
    for (p, vals), lv in zip(pts, levels):
        v, per = probe_value(lv, tau)
        draws.append({"values": vals, "value": v, "per_unit": per})
    lv0 = decisive_many([b60_point(False)], "isothermal", numerics)[0]
    v0, per0 = probe_value(lv0, tau)
    best_i = int(np.argmax([d["value"] for d in draws]))
    # the local search from the best draw, on the same log likelihood, inside the same bands
    from scipy.optimize import least_squares  # noqa: F401  (fitting.fit wraps it)

    u0 = u[best_i]

    def resid(theta, rows=None):
        lv = decisive_many([probe_point(np.clip(theta, 0.0, 1.0))[0]], "isothermal", numerics)[0]
        return probe_whitened(lv, tau)

    n_res = len(probe_whitened(levels[best_i], tau))
    prob = fitting.Problem(keys=tuple(f"u_{k}" for k in PROBE_KEYS), kinds=("acid_yield",) * len(PROBE_KEYS),
                           lower=np.zeros(len(PROBE_KEYS)), upper=np.ones(len(PROBE_KEYS)),
                           rows=tuple(range(n_res)), residuals=resid)
    f_ = fitting.fit(prob, u0, budget=PROBE_LOCAL_BUDGET)
    p_end, vals_end = probe_point(np.clip(f_.theta, 0.0, 1.0))
    v_end, per_end = probe_value(decisive_many([p_end], "isothermal", numerics)[0], tau)
    best_value = max(draws[best_i]["value"], v_end)
    kill = bool(best_value <= v0)
    best_point = probe_point(u0)[0] if draws[best_i]["value"] >= v_end else p_end
    state = {"artifact": PROBE_OUT.stem, "prereg": data_paths.rel(PREREG), "seed": PROBE_SEED, "draws_n": PROBE_DRAWS,
             "tau_ABS_all": tau, "numerics": dataclasses.asdict(numerics), "reading": "isothermal",
             "bands": {"k_ttca_cys": "B60 +- 3 dex", "ea_form_kj_mol": EA_FORM_BAND, "k_l_per_mmol": K_TTCA_PROBE_BAND,
                       "k_form_ref_per_min": f"3.29e-3 x 10^+-{K_FORM_PROBE_HALF_DEX}", "ea_decay_kj_mol": EA_DECAY_PROBE_BAND,
                       "Ea_decay_thiol_sink": (7.0, 102.0),
                       "lambda_mft_100C": "log10, prior centre - 2 dex up to the reachable maximum at the draw's Ea"},
             "r0_shipped_b": {"value": v0, "per_unit": per0},
             "best_draw": {"index": best_i, **draws[best_i]},
             "local_search": {"start_u": u0.tolist(), "end_values": vals_end, "value": v_end, "per_unit": per_end,
                              "status": f_.status, "evaluations": f_.evaluations, "budget": PROBE_LOCAL_BUDGET},
             "best_value": best_value, "kill_would_apply": kill,
             "rule": "C1 is killed only if its best value (best draw or local search end) is not higher than "
                     "R0-shipped (b)'s on the decisive units (sec. 3)",
             "consequence": ("the stacking guard: C1 is fitted once (--fit C1 --start 0), then --probe --confirm "
                             "re-applies the rule at that optimum" if kill else
                             "C1 passes the screen (screening, never evidence in the verdict)"),
             "descriptions_at_best": s1_s2(best_point, numerics),
             "draws": draws, "failed_pot_runs": len(FAILURES), "seconds": round(time.time() - t0, 1),
             "provenance": _provenance(PROBE_OUT.stem)}
    artifact_io.write_artifact(_jsonable(state), PROBE_OUT, render=render_probe)
    print(f"probe: best {best_value:.2f} vs R0-shipped (b) {v0:.2f} -> {'KILL pending confirmation' if kill else 'pass'}")
    return state


# ----------------------------------------------------------------------------------------------------------------
# --all (sec. 10)
# ----------------------------------------------------------------------------------------------------------------

def memory_plan(peak_one_fit: int, free: Optional[int]) -> str:
    """'parallel' if 2 x one fit's peak fits in free memory, 'sequential' if one does, 'refuse' otherwise."""
    if not peak_one_fit or free is None:
        return "refuse"
    if 2 * peak_one_fit <= free:
        return "parallel"
    return "sequential" if peak_one_fit <= free else "refuse"


def run_all(freeze_commit: str, log_dir: str, tau_cli: Optional[float] = None) -> int:
    commit = check_freeze(freeze_commit)
    checks = _read(CHECKS_OUT)
    problems = fit_refusals(freeze_commit, checks)
    if problems:
        raise SystemExit("--all refuses:\n  " + "\n  ".join(problems))
    probe_rec = _read(PROBE_OUT)
    if not probe_rec:
        raise SystemExit("--all: run --probe first (sec. 10 step 3)")
    confirm_first = bool(probe_rec.get("kill_would_apply") and not probe_rec.get("confirmation"))
    if (probe_rec.get("confirmation") or {}).get("killed"):
        raise SystemExit("--all: the probe's rule killed C1 at its confirmation fit (sec. 3); no M1 fits are run")
    v4_ = checks["checks"].get("V4") or {}
    peak = int(v4_.get("peak_rss_one_fit_bytes") or 0)
    free = free_memory_bytes()
    plan_ = memory_plan(peak, free)
    print(f"free memory {free} bytes; one fit's peak (V4, at 4 workers) {peak} bytes; {swap_usage()} -> {plan_}",
          flush=True)
    if plan_ == "refuse":
        raise SystemExit("--all refuses: one fit's peak memory does not fit in free memory (or V4 not run)")
    env = {**os.environ, "MAILLARD_WORKERS": str(FIT_LANE_WORKERS), "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    me = str(Path(__file__).resolve())
    tau_args = ["--tau-abs", repr(float(tau_cli))] if tau_cli is not None else []

    def cmd(*args):
        return ["nice", "-n", "10", sys.executable, me, *args, "--freeze-commit", commit]

    # sec. 10 step 4: one R0 fit and one C1 fit at a time; C1's start 2 waits for R0's carried optimum (its lowest
    # converged start), so it follows R0's starts in R0's lane while R0's Laplace and profiles run in C1's lane
    phases = [{"A": [cmd("--fit", "R0", "--start", str(s)) for s in range(3)],
               "B": [cmd("--fit", "C1", "--start", str(s)) for s in range(2)]},
              {"A": [cmd("--fit", "C1", "--start", "2")], "B": [cmd("--laplace", "R0")]},
              {"A": [cmd("--laplace", "C1")]}]
    if plan_ == "sequential":
        phases = [{"seq": [c for ph in phases for lane in ph.values() for c in lane]}]
    if confirm_first:
        # sec. 3, the stacking guard: C1 is fitted once (start 0) and the kill stands only if the rule fails there
        print("the probe's rule would kill C1: its start 0 is fitted first, then --probe --confirm", flush=True)
        if run_lanes({"confirm": [cmd("--fit", "C1", "--start", "0"),
                                  ["nice", "-n", "10", sys.executable, me, "--probe", "--confirm", *tau_args]]}, env, log_dir):
            return 1
        if (_read(PROBE_OUT).get("confirmation") or {}).get("killed"):
            print("C1 killed at its confirmation fit (sec. 3): no further M1 fits")
            return 0
    for ph in phases:
        failed = run_lanes(ph, env, log_dir)
        if failed:
            print(f"lanes that exited with an error: {failed}; stopping (sec. 10's order)")
            return 1
    tau = read_tau(tau_cli)
    if tau is None:
        print("fits and Laplace done; --score waits for H0's tau_ABS(all) (sec. 6)")
        return 0
    return subprocess.run(cmd("--score", *tau_args), cwd=str(ROOT), env=env).returncode


def run_lanes(lanes: Mapping[str, Sequence[Sequence[str]]], env: Mapping[str, str], log_dir: str) -> List[str]:
    """Each lane runs its commands in order; lanes run side by side; nothing is ever paused (no SIGSTOP)."""
    queue = {k: list(v) for k, v in lanes.items()}
    running: Dict[str, subprocess.Popen] = {}
    failed: List[str] = []
    while queue or running:
        for lane in list(queue):
            if lane not in running and queue[lane]:
                c = queue[lane].pop(0)
                log = open(Path(log_dir) / f"m1_{lane}.log", "a") if log_dir else None
                print(f"== {lane}: {' '.join(c[3:])} ({time.strftime('%H:%M:%S')})", flush=True)
                running[lane] = subprocess.Popen(c, cwd=str(ROOT), env=dict(env), stdout=log,
                                                 stderr=subprocess.STDOUT if log else None)
            if lane in queue and not queue[lane] and lane not in running:
                del queue[lane]
        time.sleep(5)
        for lane, proc in list(running.items()):
            if proc.poll() is not None:
                del running[lane]
                print(f"-- {lane} exited {proc.returncode} ({time.strftime('%H:%M:%S')})", flush=True)
                if proc.returncode != 0:
                    failed.append(lane)
                    queue.pop(lane, None)
    return failed


# ----------------------------------------------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------------------------------------------

def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def render_probe(p: Mapping[str, Any]) -> str:
    lines = ["# M1 feasibility probe (prereg sec. 3)", "",
             f"*Generated by `{GENERATED_BY} --probe`; seed {p.get('seed')}; {p.get('draws_n')} Sobol draws over C1's "
             f"bands, every other constant at B60; plug-in density on the decisive units (floors applied), isothermal, "
             f"inline, fit tolerance; tau_ABS(all) {p.get('tau_ABS_all')}. Screening, never evidence in the verdict.*", "",
             f"- R0-shipped (b): {p['r0_shipped_b']['value']:.3f}" if "r0_shipped_b" in p else "",
             f"- best draw: {p['best_draw']['value']:.3f}; local search end: {p['local_search']['value']:.3f}"
             if "best_draw" in p else "",
             f"- rule: {p.get('rule')}", f"- **kill would apply: {p.get('kill_would_apply')}**; {p.get('consequence')}", ""]
    if p.get("confirmation"):
        lines += [f"- confirmation at C1's start-0 fit: {p['confirmation']}", ""]
    return "\n".join(lines) + "\n"


def render_comparison(p: Mapping[str, Any]) -> str:
    v = p.get("verdict") or {}
    lines = ["# M1: the decisive comparison, R0 against C1 (prereg sec. 7)", "",
             f"*Generated by `{GENERATED_BY} --score` under `{p.get('prereg')}` (freeze `{str(p.get('freeze_commit'))[:10]}`); "
             f"tau_ABS(all) {p.get('tau_ABS_all')} ({p.get('tau_source')}). d_u = lppd(C1) - lppd(R0).*", "",
             f"**{v.get('text', v.get('verdict'))}**", ""]
    for reading, rec in (p.get("per_reading") or {}).items():
        lines += [f"## Reading `{reading}`: {rec['outcome']}", "",
                  "| unit | block | d linearised | d plug-in | d floors alt | d tau x0.5 | d tau x2 | side |",
                  "|---|---|---|---|---|---|---|---|"]
        for u, ur in rec["units"].items():
            d = ur["d"]
            lines.append(f"| {u} | {ur['block']} | {d['linearised']:+.2f} | {d['plug_in']:+.2f} | "
                         f"{d['floors_alternative']:+.2f} | {d['tau_x0.5']:+.2f} | {d['tau_x2']:+.2f} | {ur['side']} |")
        for b, br in rec["blocks"].items():
            lines.append(f"\n- {b}: Delta L {br['Delta_L']:+.2f}, SE {br['SE']:.2f}, sides {br['sides']}")
        lines.append("")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--checks", action="store_true", help="V1, V2, V4 - V8 and the row checks (sec. 4.4)")
    mode.add_argument("--probe", action="store_true", help="the feasibility screen (sec. 3)")
    mode.add_argument("--fit", choices=sorted(RIVALS))
    mode.add_argument("--laplace", choices=sorted(RIVALS))
    mode.add_argument("--score", action="store_true")
    mode.add_argument("--all", action="store_true")
    mode.add_argument("--v4-measure", dest="v4_measure", choices=sorted(RIVALS), help=argparse.SUPPRESS)
    ap.add_argument("--start", type=int)
    ap.add_argument("--freeze-commit", dest="freeze_commit", default="")
    ap.add_argument("--tau-abs", dest="tau_abs", type=float, default=None)
    ap.add_argument("--only", default="", help="--checks: a comma-separated subset (re-runs those)")
    ap.add_argument("--confirm", action="store_true", help="--probe: the stacking guard's re-check")
    ap.add_argument("--budget", type=int, default=BUDGET)
    ap.add_argument("--no-side-effects", dest="side_effects", action="store_false")
    ap.add_argument("--log-dir", dest="log_dir", default="")
    args = ap.parse_args(argv)
    if config.Config.from_globals() != config.Config.shipped():
        raise SystemExit(f"the process does not start from the shipped Config: {sorted(config.Config.from_globals().diff())}")
    if args.v4_measure:
        print("V4MEASURE " + json.dumps(v4_measure(args.v4_measure)), flush=True)
        return 0
    if args.checks:
        run_checks([s.strip() for s in args.only.split(",") if s.strip()] or None)
        return 0
    if args.probe:
        run_probe(args.tau_abs, args.confirm)
        return 0
    if args.fit:
        if args.start is None:
            ap.error("--fit needs --start")
        if not args.freeze_commit:
            raise SystemExit("--fit refuses without --freeze-commit")
        run_fit(args.fit, args.start, args.freeze_commit, args.budget)
        return 0
    if args.laplace:
        run_laplace(args.laplace, args.freeze_commit)
        return 0
    if args.score:
        run_score(args.freeze_commit, args.tau_abs, args.side_effects)
        return 0
    if args.all:
        return run_all(args.freeze_commit, args.log_dir, args.tau_abs)
    return 0


if __name__ == "__main__":
    sys.exit(main())

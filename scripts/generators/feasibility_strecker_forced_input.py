#!/usr/bin/env python
"""
feasibility_strecker_forced_input.py -- can the two-route forced-input Strecker layer reproduce its signatures anywhere?

Roadmap step 2 (STATUS_AND_PLAN sec. 5), under "mechanisms before constants": a FEASIBILITY check, not a fit. The
layer adopted by the design review (EX-042; results/validation/kinetic_core_strecker_forced_input_prereg.md, review
section), driven by the MEASURED courses of the glucose + Phe pot (no engine):

  dicarbonyl route   J_D = k_D [D](t) [Phe](t),  D in {MGO, GO, TDG, G}        (Hofmann 2000 Fig. 3, forced)
                     PA += (1 - phi_D) J_D;  PAA += (phi_D + omega_class(D) * O2) J_D          (arm A': C2/C3, C6)
  Amadori route      PA += [(1 - psi) k_ox O2 + (1 - psi0) k_a0] [ARP-Phe](t);  PAA += the rest
                                                                               (hofmann2000b Table 3, forced)
  aldehyde loss      d[PA]/dt -= k_loss [PA]

The fed pots fix k_D, phi_D, omega (hofmann2000 Table 4) and k_ox, psi, k_a0, psi0 (hofmann2000b Table 2) by EXACT
inversion (calibrations, no fit credit) given the drawn structural constants (CONFIG_KEYS with bands). Each draw is
scored on scale-free signatures; a bounded violation SEARCH (differential evolution on max_k dev_k / tol_k inside the
bands -- no level enters, so it is not a fit) asks whether ANY point passes. Chance baselines (uninformed bands,
scrambled forcing, flat forcing), negative controls (one route only; a slope-1 free-amine law), an oracle (synthetic
data from one draw + realistic noise: does the truth pass its own thresholds?) and a model-recovery test are reported.

Module layout follows docs/guides/A2_MODULE_CONTRACT.md sec. 3: SPECIES, REACTIONS (the engine's Species / Reaction
types), CONFIG_KEYS (bands, basis), rate_constants(params, ctx), OBSERVABLES, forced inputs {species: (t_min, value)}
interpolated piecewise-linearly; no module-level switches (routes and interpolation are config keys).

    python scripts/generators/feasibility_strecker_forced_input.py --draws 4000 --out-dir <dir>
Writes <out-dir>/kinetic_core_strecker_forced_input_feasibility[_tag].{json,md} (default out-dir results/validation).

Rows are read literature-only (data_access; status 'ready'; every pot's source verified by two readers): PHE, PHEAL,
PHEAC, PHEARP are provisional keys that fit_rows.load refuses by design. Nothing in src/ or data/ is changed.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_access, data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402
from src.kinetic_core.network import Reaction  # noqa: E402
from src.kinetic_core.species import BY_KEY as TRUNK_BY_KEY, Species  # noqa: E402
from src.kinetic_core.species_sulfur import SULFUR_SPECIES  # noqa: E402

TABLES = ("strecker_hofmann2000", "strecker_hofmann2000_fig2", "strecker_hofmann2000b")
DONORS = ("MGO", "GO", "TDG", "G")
DONOR_CLASS = {"MGO": "C23", "GO": "C23", "TDG": "C6", "G": "C6"}
OMEGA_CALIBRATOR = {"C23": "MGO", "C6": "TDG"}          # arm A': GO and G are the oxygen hold-outs
DT = 0.5                                                 # min
T_END = 500.0
D0 = AA0 = 100.0                                         # mmol/L, every fed pot (1 mmol + 1 mmol in 10 mL)
GRID = np.arange(0.0, T_END + DT / 2, DT)

# ------------------------------------------------------------------------------------------- the module as data

#: New keys (none collides with an engine key: checked in the unit test). Engine keys reused: MGO, GO, TDG, G, FRAG_C,
#: FRAG_N. Oxygen is not a species: O2 is a categorical config value (closed air vial = 1, argon = 0).
SPECIES: Tuple[Species, ...] = (
    Species("PHE", "L-phenylalanine", 9, 1, "reactant", True, "provisional key (strecker_hofmann2000.yml ruling 8)"),
    Species("PHEAL", "phenylacetaldehyde", 8, 0, "product", True, "provisional key"),
    Species("PHEAC", "phenylacetic acid", 8, 0, "product", True, "provisional key; stable (declared)"),
    Species("PHEARP", "N-(1-deoxy-D-fructos-1-yl)-L-phenylalanine (ARP-Phe)", 15, 1, "intermediate", True,
            "provisional key (strecker_hofmann2000b.yml ruling S); FORCED in the glucose pot"),
)
_POOLS = {"FRAG_C": TRUNK_BY_KEY["FRAG_C"], "FRAG_N": next(s for s in SULFUR_SPECIES if s.key == "FRAG_N")}
LOOKUP: Dict[str, Species] = {**{d: TRUNK_BY_KEY[d] for d in DONORS}, **_POOLS, **{s.key: s for s in SPECIES}}


def _strecker(d: str, exit_: str, tag: str, note: str) -> Reaction:
    """D + Phe -> PA (or PAA) + CO2 + the donor's aminocarbonyl: the CO2 and the aminocarbonyl go to FRAG_C / FRAG_N."""
    n = TRUNK_BY_KEY[d].carbon
    return Reaction(f"r_strecker_{d.lower()}_{tag}", {d: 1, "PHE": 1}, {exit_: 1, "FRAG_C": n + 1, "FRAG_N": 1},
                    None, note)


REACTIONS: Tuple[Reaction, ...] = tuple(
    [_strecker(d, "PHEAL", "aldehyde", "dicarbonyl Strecker step, aldehyde exit: (1 - phi_D) k_D") for d in DONORS]
    + [_strecker(d, "PHEAC", "acid", "parallel acid exit (NOT via PA: Hofmann 2000 Table 5): phi_D k_D")
       for d in DONORS]
    + [_strecker(d, "PHEAC", "o2_acid", "oxygen's ADDITIONAL acid channel: omega_class * O2 * k_D") for d in DONORS]
    + [Reaction(f"r_donor_sink_{d.lower()}", {d: 2}, {"FRAG_C": 2 * TRUNK_BY_KEY[d].carbon}, None,
                "fed-pot donor sink, amine-free second order (review EX-042): kappa_D [D]^2; the amine taken by the "
                "donor's other reactions is the declared nu (AA = AA0 - nu (D0 - D)); not active in forced pots")
       for d in DONORS]
    + [Reaction("r_arp_ox_aldehyde", {"PHEARP": 1}, {"PHEAL": 1, "FRAG_C": 7, "FRAG_N": 1}, None,
                "Amadori direct oxidative route (hofmann2000b Fig. 4): (1 - psi) k_ox O2"),
       Reaction("r_arp_ox_acid", {"PHEARP": 1}, {"PHEAC": 1, "FRAG_C": 7, "FRAG_N": 1}, None, "psi k_ox O2"),
       Reaction("r_arp_anox_aldehyde", {"PHEARP": 1}, {"PHEAL": 1, "FRAG_C": 7, "FRAG_N": 1}, None,
                "the fed ARP pot's argon yield, lumped (its own dicarbonyls + released Phe): (1 - psi0) k_a0"),
       Reaction("r_arp_anox_acid", {"PHEARP": 1}, {"PHEAC": 1, "FRAG_C": 7, "FRAG_N": 1}, None, "psi0 k_a0"),
       Reaction("r_pheal_loss", {"PHEAL": 1}, {"FRAG_C": 8}, None, "aldehyde loss, first order (aldol / adducts)")])


@dataclass(frozen=True)
class ConfigKey:
    """docs/guides/A2_MODULE_CONTRACT.md sec. 2 (local copy until A3 ships the type)."""
    name: str
    default: Any
    band: Optional[tuple]
    phase: str            # "build" | "rate" | "solver"
    basis: str


CONFIG_KEYS: Tuple[ConfigKey, ...] = (
    # structural (build): the routes and the forcing interpolation
    ConfigKey("strecker.route.dicarbonyl", True, (False, True), "build", "the dicarbonyl Strecker route"),
    ConfigKey("strecker.route.amadori", True, (False, True), "build", "hofmann2000b direct oxidative ARP route"),
    ConfigKey("strecker.route.o2_acid", True, (False, True), "build", "omega: oxygen's extra acid channel (arm A')"),
    ConfigKey("strecker.forcing.interp", "linear", ("linear", "early_step"), "build",
              "prereg sec. 4: piecewise-linear from [D](0) = 0; 'early_step' holds the first measured value back to "
              "t = 0+ (upper bracket on the 0 - 10 min supply; post hoc sensitivity)"),
    # drawn constants (rate): bands and basis
    ConfigKey("strecker.log10_c_mgo", -1.15, (-2.5, 0.5), "rate",
              "kappa_MGO * D0 (/min): MGO lifetime 0.3 - 300 min; review: kappa D0 ~0.07 /min from Table 3"),
    ConfigKey("strecker.nu", 0.5, (0.1, 1.0), "rate", "amine per donor consumed (review: nu = 1 contradicts Phe)"),
    ConfigKey("strecker.log10_k_loss", -3.4, (-5.0, -2.0), "rate",
              "PA loss /min (prereg fit band; Balagiannis 3-MB comparator ~3.6e-4 at 98 C)"),
    ConfigKey("strecker.off_kappa_go", 0.0, (-0.5, 0.4), "rate",
              "log10 kappa_GO / kappa_MGO (prereg +-0.5 dex; Phe mass-balance ceiling +0.4); floor 7e-5 L/(mmol min)"),
    ConfigKey("strecker.off_kappa_g", 0.0, (-0.5, 0.4), "rate", "log10 kappa_G / kappa_MGO, as GO; floor 1.7e-4"),
    ConfigKey("strecker.off_kappa_tdg", 0.0, (-0.5, 0.4), "rate",
              "log10 kappa_TDG / 8.12e-4 (Martins 2005 k_tdg_mel, 3-DG + Gly, 100 C); floor 7e-5"),
    ConfigKey("strecker.log10_k_arp_decay", -1.7, (-2.2, -1.15), "rate",
              "ARP-Phe decay /min, 100 C: >= 0.0062 (glucose pot's 30 -> 300 min fall, a lower bound) to ~2x Martins' "
              "DFG exits (0.034)"),
    ConfigKey("strecker.arp0", 100.0, (90.0, 100.0), "rate", "fed ARP-Phe mmol/L (purity 'about 90 %')"),
    ConfigKey("strecker.o2", 1.0, (0.0, 1.0), "rate", "glucose pots' O2 relative to the closed air vial (review)"),
    ConfigKey("strecker.cal_mgo", 0.0, (-0.30, 0.30), "rate", "MGO forcing calibration, log10 (+-2 x 0.15 dex)"),
    ConfigKey("strecker.cal_go", 0.0, (-0.40, 0.40), "rate", "GO forcing calibration (+-2 x 0.20 dex)"),
    ConfigKey("strecker.cal_tdg", 0.0, (-0.20, 0.20), "rate", "3-DG forcing calibration (HPLC-DAD; +-2 x 0.10 dex)"),
    ConfigKey("strecker.cal_g", 0.0, (-0.20, 0.20), "rate", "glucosone forcing calibration (+-2 x 0.10 dex)"),
    ConfigKey("strecker.cal_arp", 0.0, (-0.20, 0.20), "rate", "ARP-Phe course calibration (+-2 x 0.10 dex)"),
)
DRAWN = tuple(c for c in CONFIG_KEYS if c.phase == "rate")
KAPPA_FLOOR = {"GO": 7e-5, "G": 1.7e-4, "TDG": 7e-5}
KAPPA_TDG_MARTINS = 8.12e-4
#: The prereg's uninformed fit bands (sec. 2 / 6.2), for the chance baseline C1.
UNINFORMED = {"log10_k_d": (-7.0, -3.0), "phi": (0.0, 0.9), "omega": (0.0, 5.0),
              "log10_k_ox": (-6.0, -2.0), "log10_k_a0": (-7.0, -3.0), "psi": (0.0, 0.9)}

ARMS: Dict[str, Dict[str, Any]] = {
    "DC_ARP": {"strecker.route.dicarbonyl": True, "strecker.route.amadori": True, "strecker.route.o2_acid": True},
    "DC": {"strecker.route.dicarbonyl": True, "strecker.route.amadori": False, "strecker.route.o2_acid": True},
    "ARP": {"strecker.route.dicarbonyl": False, "strecker.route.amadori": True, "strecker.route.o2_acid": True},
    "DC_ARP_NOOMEGA": {"strecker.route.dicarbonyl": True, "strecker.route.amadori": True,
                       "strecker.route.o2_acid": False},
}


def config(arm: str = "DC_ARP", **over: Any) -> Dict[str, Any]:
    cfg = {c.name: c.default for c in CONFIG_KEYS}
    cfg.update(ARMS[arm])
    cfg.update(over)
    return cfg


def rate_constants(params: Mapping[str, Any], ctx: Mapping[str, Any]) -> Dict[str, Any]:
    """Per reaction key, the rate constant (arrays over draws allowed). ``params`` holds the calibrated constants
    (k, phi per donor; omega per class; k_ox, psi, k_a0, psi0; k_loss; kappa per donor); ``ctx`` the pot (o2) and the
    config (routes)."""
    cfg, o2 = ctx["cfg"], ctx["o2"]
    dc, arp, om = (cfg["strecker.route.dicarbonyl"], cfg["strecker.route.amadori"], cfg["strecker.route.o2_acid"])
    k: Dict[str, Any] = {}
    for d in DONORS:
        kd = params["k"][d] if dc else 0.0 * params["k"][d]
        k[f"r_strecker_{d.lower()}_aldehyde"] = (1.0 - params["phi"][d]) * kd
        k[f"r_strecker_{d.lower()}_acid"] = params["phi"][d] * kd
        k[f"r_strecker_{d.lower()}_o2_acid"] = (params["omega"][DONOR_CLASS[d]] * o2 if om else 0.0) * kd
        k[f"r_donor_sink_{d.lower()}"] = params["kappa"][d]
    a = 1.0 if arp else 0.0
    k["r_arp_ox_aldehyde"] = a * (1.0 - params["psi"]) * params["k_ox"] * o2
    k["r_arp_ox_acid"] = a * params["psi"] * params["k_ox"] * o2
    k["r_arp_anox_aldehyde"] = a * (1.0 - params["psi0"]) * params["k_a0"]
    k["r_arp_anox_acid"] = a * params["psi0"] * params["k_a0"]
    k["r_pheal_loss"] = params["k_loss"]
    return k


OBSERVABLES: Dict[str, Callable[[Mapping[str, Any]], Any]] = {
    "PA": lambda c: c["PHEAL"],
    "PAA": lambda c: c["PHEAC"],
    "PA+PAA": lambda c: c["PHEAL"] + c["PHEAC"],
    "PA/PAA": lambda c: c["PHEAL"] / c["PHEAC"],
}

#: v1 tolerances (prereg sec. 6.1). T1b (in-sample) added after the smoke run: Table 3's PA tail, k_loss's identifier
#: (prereg sec. 2). G6b in dex of the crossover time (x2).
TOL = {"T1": 0.10, "T1b": 0.06, "G3": 0.15, "G4": 0.15, "G5": 0.15, "G6a": 0.15, "G6b": math.log10(2.0)}
JOINT = ("T1", "T1b", "G3", "G4", "G5", "G6a", "G6b")
TOL_OXYGEN = 0.20
T1_TIMES = (10.0, 30.0, 60.0, 120.0)
G5_TIMES = (10.0, 30.0, 60.0, 240.0, 500.0)


# ----------------------------------------------------------------------------------------------------------- data

def _read_table(name: str) -> Tuple[Dict[str, dict], List[dict]]:
    """Literature-only read: status 'ready' and every pot verified by two readers (the loader's own gates), without the
    engine-species check that refuses the provisional Strecker keys."""
    table = data_access.load_yaml(fit_rows.FIT_ROWS_DIR / f"{name}.yml")
    if table.get("status") != "ready":
        raise PermissionError(f"{name}: status {table.get('status')!r}, not 'ready'")
    verified = fit_rows._verified(name)
    pots = table.get("pots") or {}
    for pot_id, pot in pots.items():
        if not verified.get(pot["ledger_id"], False):
            raise PermissionError(f"{name}: pot {pot_id}'s source {pot['ledger_id']} is not verified by two readers")
    rows = [dict(r, sigma=fit_rows._sigma(r)) for r in table.get("rows") or []]
    return pots, rows


@dataclass(frozen=True)
class Data:
    """Every observation the check uses (mmol/L), by role. CF-074: hofmann2000b Table 2 Glc/Phe air is the Fig. 2 run
    reported again -- it is never scored (diagnostic glc_2000b only); CF-070: the Phe course is hofmann2000b's."""
    t3: Dict[str, Dict[float, float]]                  # hofmann2000 Table 3 MGO fed pot: PHEAL / PHEAC by time
    fed: Dict[str, Dict[str, Tuple[float, float]]]     # donor -> atmosphere -> (PA, PAA) at 60 min (Table 4)
    fed_ph: Dict[str, Dict[float, float]]              # donor -> pH -> PA + PAA (Figs. 4-5; negative control only)
    arp_fed: Dict[str, Tuple[float, float]]            # 'argon' / 'air' -> (PA, PAA), 120 min, 100 C (2000b Table 2)
    glc_t2: Tuple[float, float]                        # hofmann2000 Table 2, 30 min (PA, PAA)
    glc_ph7: Tuple[float, float]                       # Figs. 4-5 pH 7, 60 min (declared conditions)
    fig2: Dict[str, Dict[float, float]]                # Fig. 2: PHEAL / PHEAC by time
    glc_2000b: Dict[str, Tuple[float, float]]          # hofmann2000b Table 2 glucose pot, 120 min (diagnostic)
    forced_inputs: Dict[str, Tuple[np.ndarray, np.ndarray]]   # {species: (t_min, value)} incl. t = 0


def load_data() -> Data:
    rows: Dict[Tuple[str, str, float], float] = {}
    for name in TABLES:
        _, rs = _read_table(name)
        for r in rs:
            rows[(r["pot"], r["species"], float(r["t_min"]))] = float(r["value"])

    def v(pot, sp, t):
        return rows[(pot, sp, float(t))]

    t3pot = "hofmann2000_t3_mgo_phe_98c_closed"
    t3 = {sp: {t: v(t3pot, sp, t) for t in (10, 30, 60, 120, 240, 300, 500)} for sp in ("PHEAL", "PHEAC")}
    fed = {d: {atm: (v(f"hofmann2000_t4_{d.lower()}_phe_60min_{atm}", "PHEAL", 60),
                     v(f"hofmann2000_t4_{d.lower()}_phe_60min_{atm}", "PHEAC", 60)) for atm in ("argon", "air")}
           for d in DONORS}
    fed_ph = {d: {ph: sum(v(f"hofmann2000_ph_{d.lower()}_phe_ph{ph}", sp, 60) for sp in ("PHEAL", "PHEAC"))
                  for ph in (3, 7)} for d in DONORS}
    arp_fed = {atm: (v(f"hofmann2000b_t2_arp_phe_100c_120min_{atm}", "PHEAL", 120),
                     v(f"hofmann2000b_t2_arp_phe_100c_120min_{atm}", "PHEAC", 120)) for atm in ("argon", "air")}
    glc_2000b = {atm: (v(f"hofmann2000b_t2_glc_phe_100c_120min_{atm}", "PHEAL", 120),
                       v(f"hofmann2000b_t2_glc_phe_100c_120min_{atm}", "PHEAC", 120)) for atm in ("argon", "air")}
    f2pot = "hofmann2000_fig2_glc_phe_98c"
    fig2 = {sp: {t: v(f2pot, sp, t) for t in (10, 30, 60, 120, 240, 500)} for sp in ("PHEAL", "PHEAC")}
    f3pot = "hofmann2000_fig3_glc_phe_98c"
    forced_inputs = {}
    for d in DONORS:
        ts = [0.0] + [float(t) for t in (10, 30, 60, 120, 240, 300, 500) if (f3pot, d, float(t)) in rows]
        vs = [0.0] + [v(f3pot, d, t) for t in ts[1:]]
        if d == "G":   # glucosone's 300 / 500-min markers sit on the axis: zero (keyed header, left out as rows)
            ts += [300.0, 500.0]
            vs += [0.0, 0.0]
        forced_inputs[d] = (np.array(ts), np.array(vs))
    t3b = "hofmann2000b_t3_glc_phe_reflux"
    for sp, start in (("PHE", AA0), ("PHEARP", 0.0)):
        ts = [0.0, 10.0, 30.0, 60.0, 120.0, 300.0]
        vs = [start] + [v(t3b, sp, t) for t in ts[1:]]
        slope = math.log(vs[-1] / vs[-2]) / (ts[-1] - ts[-2])   # declared: log-linear 120 -> 300 extended to 500
        forced_inputs[sp] = (np.array(ts + [500.0]), np.array(vs + [vs[-1] * math.exp(slope * 200.0)]))
    return Data(t3=t3, fed=fed, fed_ph=fed_ph, arp_fed=arp_fed,
                glc_t2=(v("hofmann2000_t2_glc_phe_reflux_30min", "PHEAL", 30),
                        v("hofmann2000_t2_glc_phe_reflux_30min", "PHEAC", 30)),
                glc_ph7=(v("hofmann2000_ph_glc_phe_ph7", "PHEAL", 60), v("hofmann2000_ph_glc_phe_ph7", "PHEAC", 60)),
                fig2=fig2, glc_2000b=glc_2000b, forced_inputs=forced_inputs)


# ------------------------------------------------------------------------------------------------------ numerics

def idx(t: float) -> int:
    return int(round(t / DT))


def cum_trapz(f: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid along the last axis (time) on GRID."""
    out = np.zeros_like(f)
    out[..., 1:] = np.cumsum(0.5 * DT * (f[..., 1:] + f[..., :-1]), axis=-1)
    return out


def cum_with_loss(f: np.ndarray, k_loss) -> np.ndarray:
    """L(t) = int_0^t exp(-k (t - s)) f(s) ds = exp(-k t) int_0^t exp(k s) f(s) ds, per draw (rows); k <= 0.01 /min
    on 500 min keeps exp(k t) <= e^5."""
    k = np.asarray(k_loss, dtype=float).reshape(-1, 1) if f.ndim == 2 else float(k_loss)
    return np.exp(-k * GRID) * cum_trapz(f * np.exp(k * GRID))


def fed_donor_courses(kappa: np.ndarray, nu: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fed pot, amine-free second-order donor sink: D = D0 / (1 + kappa D0 t); AA = AA0 - nu (D0 - D) (>= 0)."""
    d = D0 / (1.0 + kappa[:, None] * D0 * GRID[None, :])
    aa = np.maximum(AA0 - nu[:, None] * (D0 - d), 0.0)
    return d, aa


def interpolate(forced_inputs: Mapping[str, Tuple[np.ndarray, np.ndarray]], interp: str = "linear"
                ) -> Dict[str, np.ndarray]:
    """{species: course on GRID}, piecewise-linear in concentration (A2 contract sec. 3.2)."""
    out = {}
    for key, (ts, vs) in forced_inputs.items():
        if interp == "early_step" and key != "PHE":
            ts, vs = np.concatenate([[0.0, DT], ts[1:]]), np.concatenate([[0.0, vs[1]], vs[1:]])
        out[key] = np.interp(GRID, ts, vs)
    return out


# --------------------------------------------------------------------------------------------------------- draws

def _band(c: ConfigKey, o2_band: Tuple[float, float]) -> Tuple[float, float]:
    return o2_band if c.name == "strecker.o2" else c.band


def draw(n: int, rng: np.random.Generator, o2_band: Tuple[float, float] = (0.0, 1.0)) -> Dict[str, np.ndarray]:
    th = {}
    for c in DRAWN:
        lo, hi = _band(c, o2_band)
        th[c.name.split(".", 1)[1]] = lo + (hi - lo) * rng.random(n)
    return th


def from_unit(u: np.ndarray, o2_band: Tuple[float, float] = (0.0, 1.0)) -> Dict[str, np.ndarray]:
    """(n_keys, n) unit-cube coordinates -> theta (the search's parametrisation)."""
    th = {}
    for i, c in enumerate(DRAWN):
        lo, hi = _band(c, o2_band)
        th[c.name.split(".", 1)[1]] = lo + (hi - lo) * u[i]
    return th


def kappas(th: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    k_mgo = 10.0 ** th["log10_c_mgo"] / D0
    return {"MGO": k_mgo,
            "GO": np.maximum(k_mgo * 10.0 ** th["off_kappa_go"], KAPPA_FLOOR["GO"]),
            "G": np.maximum(k_mgo * 10.0 ** th["off_kappa_g"], KAPPA_FLOOR["G"]),
            "TDG": np.maximum(KAPPA_TDG_MARTINS * 10.0 ** th["off_kappa_tdg"], KAPPA_FLOOR["TDG"])}


def calibrate(th: Mapping[str, np.ndarray], data: Data) -> Dict[str, Any]:
    """Exact inversion of the fed pots for the drawn structure -> the module's params: k_D, phi_D (Table 4 argon),
    omega per class (Table 4 air of the class calibrator), k_ox, psi, k_a0, psi0 (hofmann2000b Table 2 ARP pots);
    plus the Table 3 courses for T1 / T1b."""
    n = len(th["nu"])
    k_loss = 10.0 ** th["log10_k_loss"]
    i60 = idx(60.0)
    kap = kappas(th)
    cal: Dict[str, Any] = {"k": {}, "phi": {}, "omega": {}, "kappa": kap, "k_loss": k_loss}
    courses = {}
    for d in DONORS:
        dd, aa = fed_donor_courses(kap[d], th["nu"])
        f = dd * aa
        integ, lossy = cum_trapz(f), cum_with_loss(f, k_loss)
        courses[d] = (integ, lossy)
        pa, paa = data.fed[d]["argon"]
        a, b = pa / lossy[:, i60], paa / integ[:, i60]
        cal["k"][d] = a + b
        cal["phi"][d] = b / (a + b)
    for cls, d in OMEGA_CALIBRATOR.items():
        cal["omega"][cls] = data.fed[d]["air"][1] / (cal["k"][d] * courses[d][0][:, i60]) - cal["phi"][d]
    integ, lossy = courses["MGO"]
    k, phi, om = cal["k"]["MGO"][:, None], cal["phi"]["MGO"][:, None], cal["omega"]["C23"][:, None]
    cal["t3_pa"] = (1 - phi) * k * lossy                 # Table 3: closed vial, air (O2 = 1)
    cal["t3_total"] = cal["t3_pa"] + (phi + om) * k * integ
    kd = 10.0 ** th["log10_k_arp_decay"]
    arp_int = th["arp0"] * (1.0 - np.exp(-120.0 * kd)) / kd
    pa_ar, paa_ar = data.arp_fed["argon"]
    pa_air, paa_air = data.arp_fed["air"]
    cal["k_a0"], cal["psi0"] = (pa_ar + paa_ar) / arp_int, np.full(n, paa_ar / (pa_ar + paa_ar))
    ex_pa, ex_paa = pa_air - pa_ar, paa_air - paa_ar
    cal["k_ox"], cal["psi"] = (ex_pa + ex_paa) / arp_int, np.full(n, ex_paa / (ex_pa + ex_paa))
    return cal


def uninformed_calibration(th: Mapping[str, np.ndarray], data: Data, rng: np.random.Generator) -> Dict[str, Any]:
    """Chance baseline C1 (prereg sec. 6.2): the calibrated constants replaced by uniform draws in their fit bands."""
    n = len(th["nu"])
    u = lambda lo_hi: lo_hi[0] + (lo_hi[1] - lo_hi[0]) * rng.random(n)  # noqa: E731
    cal = calibrate(th, data)
    cal["k"] = {d: 10.0 ** u(UNINFORMED["log10_k_d"]) for d in DONORS}
    cal["phi"] = {d: u(UNINFORMED["phi"]) for d in DONORS}
    cal["omega"] = {c: u(UNINFORMED["omega"]) for c in ("C23", "C6")}
    dd, aa = fed_donor_courses(cal["kappa"]["MGO"], th["nu"])
    f = dd * aa
    k, phi, om = cal["k"]["MGO"][:, None], cal["phi"]["MGO"][:, None], cal["omega"]["C23"][:, None]
    cal["t3_pa"] = (1 - phi) * k * cum_with_loss(f, cal["k_loss"])
    cal["t3_total"] = cal["t3_pa"] + (phi + om) * k * cum_trapz(f)
    cal["k_ox"], cal["k_a0"] = 10.0 ** u(UNINFORMED["log10_k_ox"]), 10.0 ** u(UNINFORMED["log10_k_a0"])
    cal["psi"], cal["psi0"] = u(UNINFORMED["psi"]), u(UNINFORMED["psi"])
    return cal


# --------------------------------------------------------------------------------------------------- the layer

def glucose_pot(th: Mapping[str, np.ndarray], cal: Mapping[str, Any], courses: Mapping[str, np.ndarray],
                cfg: Mapping[str, Any], o2: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """PA, PAA (n x T) in the forced glucose + Phe pot from the module's rate constants, and the routes' parts."""
    n = len(th["nu"])
    k = rate_constants(cal, {"cfg": cfg, "o2": th["o2"] if o2 is None else o2})
    src = {"PA_dc": np.zeros((n, GRID.size)), "PAA_dc": np.zeros((n, GRID.size))}
    for d in DONORS:
        enc = (10.0 ** th[f"cal_{d.lower()}"])[:, None] * (courses[d] * courses["PHE"])[None, :]
        src["PA_dc"] += np.asarray(k[f"r_strecker_{d.lower()}_aldehyde"])[:, None] * enc
        src["PAA_dc"] += np.asarray(k[f"r_strecker_{d.lower()}_acid"] + k[f"r_strecker_{d.lower()}_o2_acid"]
                                    )[:, None] * enc
    arp = (10.0 ** th["cal_arp"])[:, None] * courses["PHEARP"][None, :]
    src["PA_arp"] = np.asarray(k["r_arp_ox_aldehyde"] + k["r_arp_anox_aldehyde"])[:, None] * arp
    src["PAA_arp"] = np.asarray(k["r_arp_ox_acid"] + k["r_arp_anox_acid"])[:, None] * arp
    out = {"PA_dc": cum_with_loss(src["PA_dc"], k["r_pheal_loss"]),
           "PA_arp": cum_with_loss(src["PA_arp"], k["r_pheal_loss"]),
           "PAA_dc": cum_trapz(src["PAA_dc"]), "PAA_arp": cum_trapz(src["PAA_arp"])}
    out["PA"], out["PAA"] = out["PA_dc"] + out["PA_arp"], out["PAA_dc"] + out["PAA_arp"]
    return out


def _crossover(pa: np.ndarray, paa: np.ndarray) -> np.ndarray:
    """First time PAA >= PA (after the first step); inf if never by 500 min."""
    hit = paa[:, 2:] >= pa[:, 2:]
    t = GRID[2:][np.argmax(hit, axis=1)]
    return np.where(hit.any(axis=1), t, np.inf)


def measured_signatures(data: Data) -> Dict[str, Any]:
    t3 = {t: data.t3["PHEAL"][t] + data.t3["PHEAC"][t] for t in data.t3["PHEAL"]}
    f2 = {t: data.fig2["PHEAL"][t] + data.fig2["PHEAC"][t] for t in data.fig2["PHEAL"]}
    lr = {t: math.log10(data.fig2["PHEAL"][t] / data.fig2["PHEAC"][t]) for t in data.fig2["PHEAL"]}
    times = sorted(lr)
    cross = math.inf
    for a, b in zip(times, times[1:]):
        if lr[a] > 0 >= lr[b]:
            cross = a + (b - a) * lr[a] / (lr[a] - lr[b])
            break
    o2_ratio = {d: data.fed[d]["air"][1] / data.fed[d]["argon"][1] for d in DONORS}
    return {"T1": [math.log10(t3[t] / t3[500.0]) for t in T1_TIMES],
            "T1b": math.log10(data.t3["PHEAL"][500.0] / data.t3["PHEAL"][240.0]),
            "G3": math.log10(data.glc_t2[0] / data.glc_t2[1]),
            "G4": math.log10(data.glc_ph7[0] / data.glc_ph7[1]),
            "G5": [math.log10(f2[t] / f2[120.0]) for t in G5_TIMES],
            "G6a": math.log10(data.fig2["PHEAC"][500.0] / data.fig2["PHEAL"][500.0]),
            "G6b": cross,
            "O2": {d: math.log10(r) for d, r in o2_ratio.items()},
            "O3": sorted(DONORS, key=lambda d: -o2_ratio[d]),
            "null_a_G3": math.log10(float(np.mean([data.fed[d]["air"][0] / data.fed[d]["air"][1] for d in DONORS]))),
            "G1_total_30": sum(data.glc_t2), "G2_total_60": sum(data.glc_ph7), "fig2_total_120": f2[120.0],
            "A1_pa_air_argon": data.glc_2000b["air"][0] / data.glc_2000b["argon"][0],
            "A2_paa_air_argon": data.glc_2000b["air"][1] / data.glc_2000b["argon"][1]}


def score(th: Mapping[str, np.ndarray], data: Data, cfg: Mapping[str, Any], *, cal: Optional[Mapping] = None,
          courses: Optional[Mapping[str, np.ndarray]] = None, meas: Optional[Mapping] = None
          ) -> Dict[str, np.ndarray]:
    """Per-draw deviations |log10 model - log10 measured| (dex) of every scored signature, plus diagnostics."""
    cal = cal or calibrate(th, data)
    meas = meas or measured_signatures(data)
    courses = courses or interpolate(data.forced_inputs, cfg["strecker.forcing.interp"])
    g = glucose_pot(th, cal, courses, cfg)
    pa, paa = np.maximum(g["PA"], 1e-300), np.maximum(g["PAA"], 1e-300)
    tot = pa + paa
    t3, t3pa = cal["t3_total"], np.maximum(cal["t3_pa"], 1e-300)
    out: Dict[str, np.ndarray] = {}
    out["T1"] = np.max([np.abs(np.log10(t3[:, idx(t)] / t3[:, idx(500.0)]) - m) for t, m in zip(T1_TIMES, meas["T1"])],
                       axis=0)
    out["T1b"] = np.abs(np.log10(t3pa[:, idx(500.0)] / t3pa[:, idx(240.0)]) - meas["T1b"])
    out["G3"] = np.abs(np.log10(pa[:, idx(30)] / paa[:, idx(30)]) - meas["G3"])
    out["G4"] = np.abs(np.log10(pa[:, idx(60)] / paa[:, idx(60)]) - meas["G4"])
    g5 = [np.log10(tot[:, idx(t)] / tot[:, idx(120)]) - m for t, m in zip(G5_TIMES, meas["G5"])]
    out["G5"] = np.max(np.abs(g5), axis=0)
    out["G5_no10"] = np.max(np.abs(g5[1:]), axis=0)
    out["G5_at10"], out["G5_at500"] = g5[0], g5[-1]
    out["G6a"] = np.abs(np.log10(paa[:, idx(500)] / pa[:, idx(500)]) - meas["G6a"])
    cross = _crossover(pa, paa)
    out["crossover_min"] = cross
    with np.errstate(invalid="ignore", divide="ignore"):
        if math.isfinite(meas["G6b"]):
            out["G6b"] = np.where(np.isfinite(cross), np.abs(np.log10(np.maximum(cross, DT) / meas["G6b"])), np.inf)
        else:   # no crossover measured by 500 min (a synthetic truth): a model crossover is judged against 500 min
            out["G6b"] = np.where(np.isfinite(cross), np.abs(np.log10(np.maximum(cross, DT) / T_END)), 0.0)
    out["G6b_search"] = np.where(np.isfinite(out["G6b"]), out["G6b"],   # finite stand-in when PAA never overtakes
                                 1.0 + np.log10(pa[:, idx(500)] / paa[:, idx(500)]).clip(min=0.0))
    # diagnostics (levels and route shares are NOT scale-free and NOT scored)
    out["R30"] = meas["G1_total_30"] / tot[:, idx(30)]
    out["R60"] = meas["G2_total_60"] / tot[:, idx(60)]
    out["R120_fig2"] = meas["fig2_total_120"] / tot[:, idx(120)]
    out["arp_share_PA_30"] = g["PA_arp"][:, idx(30)] / pa[:, idx(30)]
    out["arp_share_PA_120"] = g["PA_arp"][:, idx(120)] / pa[:, idx(120)]
    out["late_gain_500_240_model"] = np.log10(tot[:, idx(500)] / tot[:, idx(240)])
    g0 = glucose_pot(th, cal, courses, cfg, o2=np.zeros(len(th["nu"])))   # same forcing under argon (CF-074)
    out["A1_pa_air_argon"] = pa[:, idx(120)] / np.maximum(g0["PA"][:, idx(120)], 1e-300)
    out["A2_paa_air_argon"] = paa[:, idx(120)] / np.maximum(g0["PAA"][:, idx(120)], 1e-300)
    ratios = {d: np.maximum(1.0 + cal["omega"][DONOR_CLASS[d]] / cal["phi"][d], 1e-6) for d in DONORS}
    for d in ("GO", "G"):                                 # arm A' hold-outs: 1 + omega_class / phi_D
        out[f"O2_{d}"] = np.abs(np.log10(ratios[d]) - meas["O2"][d])
    order = np.argsort(-np.vstack([ratios[d] for d in DONORS]), axis=0)
    want = np.array([DONORS.index(d) for d in meas["O3"]])[:, None]
    out["O3_exact"] = np.all(order == want, axis=0)
    out["O3_class_order"] = (np.minimum(ratios["MGO"], ratios["GO"]) > np.maximum(ratios["TDG"], ratios["G"]))
    out["G3_beats_null_a"] = out["G3"] < abs(meas["null_a_G3"] - meas["G3"])
    return out


def violation(sig: Mapping[str, np.ndarray], tol: Mapping[str, float], keys: Sequence[str] = JOINT) -> np.ndarray:
    return np.max([(sig["G6b_search"] if k == "G6b" else sig[k]) / tol[k] for k in keys], axis=0)


# --------------------------------------------------------------------------------------- sampling and search

def _chance_courses(data: Data, rng: np.random.Generator, kind: str, interp: str) -> Dict[str, np.ndarray]:
    base = interpolate(data.forced_inputs, interp)
    fc = dict(base)
    if kind == "scrambled":
        perm = list(DONORS)
        while perm == list(DONORS):
            perm = list(rng.permutation(DONORS))
        for d, src in zip(DONORS, perm):
            fc[d] = base[src]
        ts, vs = data.forced_inputs["PHEARP"]
        inner = vs[1:-1].copy()
        while np.allclose(inner, vs[1:-1]):
            inner = rng.permutation(vs[1:-1])
        fc["PHEARP"] = np.interp(GRID, ts, np.concatenate([[0.0], inner, [vs[-1]]]))
    elif kind == "flat":
        for key in DONORS + ("PHEARP",):
            fc[key] = np.full(GRID.size, float(np.mean(base[key])))
    return fc


def _sl(d: Mapping[str, Any], sl) -> Dict[str, Any]:
    return {k: (_sl(v, sl) if isinstance(v, dict) else v[sl]) for k, v in d.items()}


def run_arm(cfg: Mapping[str, Any], n: int, seed: int, data: Data, tol: Mapping[str, float], *,
            chance: Optional[str] = None, o2_band=(0.0, 1.0), chunk: int = 1000) -> Dict[str, Any]:
    """n uniform draws of one structure (or a chance baseline on it): reach per signature, joint passes."""
    rng = np.random.default_rng(seed)
    meas = measured_signatures(data)
    sigs: Dict[str, List[np.ndarray]] = {}
    thetas: List[Dict[str, np.ndarray]] = []
    for start in range(0, n, chunk):
        m = min(chunk, n - start)
        th = draw(m, rng, o2_band)
        cal = uninformed_calibration(th, data, rng) if chance == "uninformed" else calibrate(th, data)
        if chance in ("scrambled", "flat"):   # one scrambled forcing per block of 50 draws
            parts = []
            for s0 in range(0, m, 50):
                sl = slice(s0, min(s0 + 50, m))
                fc = _chance_courses(data, rng, chance, cfg["strecker.forcing.interp"])
                parts.append(score(_sl(th, sl), data, cfg, cal=_sl(cal, sl), courses=fc, meas=meas))
            s = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        else:
            s = score(th, data, cfg, cal=cal, meas=meas)
        for k, v in s.items():
            sigs.setdefault(k, []).append(v)
        thetas.append(th)
    sig = {k: np.concatenate(v) for k, v in sigs.items()}
    theta = {k: np.concatenate([t[k] for t in thetas]) for k in thetas[0]}
    p = {k: sig[k] <= tol[k] for k in JOINT}
    p["joint"] = np.all([p[k] for k in JOINT], axis=0)
    norm = violation(sig, tol)
    return {"n": n, "seed": seed, "chance": chance, "sig": sig, "theta": theta, "pass": p,
            "reach": {k: int(p[k].sum()) for k in JOINT + ("joint",)}, "best_index": int(np.argmin(norm)),
            "best_norm": float(np.min(norm)), "joint_idx": np.flatnonzero(p["joint"])}


SEARCH_REPORT = JOINT + ("G5_no10", "G5_at10", "G5_at500", "R30", "R120_fig2", "arp_share_PA_30", "arp_share_PA_120",
                         "crossover_min", "A1_pa_air_argon", "A2_paa_air_argon", "O2_GO", "O2_G")


def search(cfg: Mapping[str, Any], data: Data, tol: Mapping[str, float], seed: int, *,
           courses: Optional[Mapping[str, np.ndarray]] = None, keys: Sequence[str] = JOINT, o2_band=(0.0, 1.0),
           maxiter: int = 60, popsize: int = 12) -> Dict[str, Any]:
    """Bounded violation search: min over the bands of max_k dev_k / tol_k (differential evolution, vectorised).
    A minimum <= 1 is a constant set inside the bands that passes every signature in ``keys``."""
    from scipy.optimize import differential_evolution
    meas = measured_signatures(data)
    courses = courses or interpolate(data.forced_inputs, cfg["strecker.forcing.interp"])

    def f(u):
        u = u[:, None] if u.ndim == 1 else u
        return violation(score(from_unit(u, o2_band), data, cfg, courses=courses, meas=meas), tol, keys)

    res = differential_evolution(f, [(0.0, 1.0)] * len(DRAWN), seed=seed, maxiter=maxiter, popsize=popsize,
                                 tol=1e-8, polish=False, vectorized=True, updating="deferred", init="sobol")
    th = from_unit(np.asarray(res.x)[:, None], o2_band)
    s = score(th, data, cfg, courses=courses, meas=meas)
    return {"min_violation": float(np.atleast_1d(res.fun)[0]), "passes": bool(np.atleast_1d(res.fun)[0] <= 1.0),
            "theta": {k: float(v[0]) for k, v in th.items()},
            "dev": {k: float(s[k][0]) for k in SEARCH_REPORT}, "nfev": int(res.nfev)}


# --------------------------------------------------------------------------------------- negative control 1

def free_amine_control(data: Data, n: int, seed: int) -> Dict[str, Any]:
    """Fed totals pH 3 / pH 7 (Figs. 4-5, four donors, declared conditions). The layer's total-amine law predicts 1;
    a slope-1 free-amine law with pKa(Phe alpha-NH3+, 98 C) 7.6 +- 0.6 predicts f(3) / f(7). Tolerance 0.30 dex."""
    rng = np.random.default_rng(seed)
    pka = 7.0 + 1.2 * rng.random(n)
    f = lambda ph: 1.0 / (1.0 + 10.0 ** (pka - ph))  # noqa: E731
    pred_free = np.log10(f(3.0) / f(7.0))
    out = {}
    for d in DONORS:
        m = math.log10(data.fed_ph[d][3] / data.fed_ph[d][7])
        out[d] = {"measured_dex": round(m, 3), "total_amine_dev": round(abs(m), 3),
                  "free_amine_dev_min": round(float(np.min(np.abs(pred_free - m))), 2),
                  "free_amine_dev_max": round(float(np.max(np.abs(pred_free - m))), 2)}
    return out


# ------------------------------------------------------------------------------------------------------- oracle

NOISE_REALISTIC = {"PA": 0.02, "PAA": 0.06, "pot_fed": 0.02, "pot_glc": 0.10, "forcing_point": 0.05,
                   "cal": {"MGO": 0.15, "GO": 0.20, "TDG": 0.10, "G": 0.10, "PHEARP": 0.10, "PHE": 0.02}}
NOISE_STRESS = {"PA": 0.15, "PAA": 0.15, "pot_fed": 0.10, "pot_glc": 0.10, "forcing_point": 0.15,
                "cal": {"MGO": 0.15, "GO": 0.20, "TDG": 0.10, "G": 0.10, "PHEARP": 0.10, "PHE": 0.02}}


def truth_rows(th1: Mapping[str, np.ndarray], data: Data, cfg: Mapping[str, Any]) -> Data:
    """Noise-free synthetic observations from one draw (n = 1): the fed pots from its calibrated constants, the
    glucose pots from the layer under the central forcing."""
    cal = calibrate(th1, data)
    i60 = idx(60.0)
    fed = {}
    for d in DONORS:
        dd, aa = fed_donor_courses(cal["kappa"][d], th1["nu"])
        integ, lossy = cum_trapz(dd * aa)[0, i60], cum_with_loss(dd * aa, cal["k_loss"])[0, i60]
        k, phi, om = cal["k"][d][0], cal["phi"][d][0], cal["omega"][DONOR_CLASS[d]][0]
        pa = (1 - phi) * k * lossy
        fed[d] = {"argon": (pa, phi * k * integ), "air": (pa, (phi + om) * k * integ)}
    t3pa, t3tot = cal["t3_pa"][0], cal["t3_total"][0]
    t3 = {"PHEAL": {t: float(t3pa[idx(t)]) for t in data.t3["PHEAL"]},
          "PHEAC": {t: float(t3tot[idx(t)] - t3pa[idx(t)]) for t in data.t3["PHEAC"]}}
    kd = 10.0 ** th1["log10_k_arp_decay"][0]
    arp_int = th1["arp0"][0] * (1 - math.exp(-120 * kd)) / kd
    ka0, psi0, kox, psi = cal["k_a0"][0], cal["psi0"][0], cal["k_ox"][0], cal["psi"][0]
    ar = (ka0 * arp_int * (1 - psi0), ka0 * arp_int * psi0)
    arp_fed = {"argon": ar, "air": (ar[0] + kox * arp_int * (1 - psi), ar[1] + kox * arp_int * psi)}
    courses = interpolate(data.forced_inputs, cfg["strecker.forcing.interp"])
    g = glucose_pot(th1, cal, courses, cfg)
    pa, paa = g["PA"][0], g["PAA"][0]
    g0 = glucose_pot(th1, cal, courses, cfg, o2=np.zeros(1))
    return replace(data, t3=t3, fed=fed, arp_fed=arp_fed, glc_t2=(float(pa[idx(30)]), float(paa[idx(30)])),
                   glc_ph7=(float(pa[idx(60)]), float(paa[idx(60)])),
                   fig2={"PHEAL": {t: float(pa[idx(t)]) for t in data.fig2["PHEAL"]},
                         "PHEAC": {t: float(paa[idx(t)]) for t in data.fig2["PHEAC"]}},
                   glc_2000b={"argon": (float(g0["PA"][0, idx(120)]), float(g0["PAA"][0, idx(120)])),
                              "air": (float(pa[idx(120)]), float(paa[idx(120)]))})


def noisy(truth: Data, rng: np.random.Generator, noise: Mapping[str, Any]) -> Data:
    """Observed = truth x log-normal noise: a per-pot run factor (common to PA and PAA of a pot; fed vs glucose pots)
    and a per-measurement term (PA / PAA floors); the forcing gets a per-series calibration and a per-point term."""
    def lp(x, s):
        return x * 10.0 ** (s * rng.standard_normal())

    def pot(pa_paa, kind):
        f = 10.0 ** (noise[kind] * rng.standard_normal())
        return (lp(pa_paa[0], noise["PA"]) * f, lp(pa_paa[1], noise["PAA"]) * f)

    def course(c, kind):
        f = 10.0 ** (noise[kind] * rng.standard_normal())
        return {"PHEAL": {t: lp(v, noise["PA"]) * f for t, v in c["PHEAL"].items()},
                "PHEAC": {t: lp(v, noise["PAA"]) * f for t, v in c["PHEAC"].items()}}

    forced_inputs = {}
    for key, (ts, vs) in truth.forced_inputs.items():
        new = vs * 10.0 ** (noise["cal"][key] * rng.standard_normal())
        if key != "PHE":
            new = new * 10.0 ** (noise["forcing_point"] * rng.standard_normal(vs.size))
        else:
            new[0] = vs[0]
        forced_inputs[key] = (ts, new)
    return replace(truth, fed={d: {a: pot(truth.fed[d][a], "pot_fed") for a in ("argon", "air")} for d in DONORS},
                   t3=course(truth.t3, "pot_fed"), fig2=course(truth.fig2, "pot_glc"),
                   arp_fed={a: pot(truth.arp_fed[a], "pot_fed") for a in ("argon", "air")},
                   glc_t2=pot(truth.glc_t2, "pot_glc"), glc_ph7=pot(truth.glc_ph7, "pot_glc"),
                   glc_2000b={a: pot(truth.glc_2000b[a], "pot_glc") for a in ("argon", "air")},
                   forced_inputs=forced_inputs)


def oracle(th1: Mapping[str, np.ndarray], data: Data, cfg: Mapping[str, Any], reps: int, seed: int,
           noise: Mapping[str, Any]) -> Dict[str, Any]:
    """Does the TRUE draw pass its own thresholds on synthetic data with noise? The truth is re-calibrated from the
    noisy fed rows and driven by the noisy forcing; its deviations give a noise-calibrated tolerance per signature."""
    truth = truth_rows(th1, data, cfg)
    rng = np.random.default_rng(seed)
    keys = JOINT + ("O2_GO", "O2_G")
    devs: Dict[str, List[float]] = {k: [] for k in keys}
    o3, o3c = [], []
    for _ in range(reps):
        s = score(th1, noisy(truth, rng, noise), cfg)
        for k in keys:
            devs[k].append(float(s[k][0]))
        o3.append(bool(s["O3_exact"][0]))
        o3c.append(bool(s["O3_class_order"][0]))
    arr = {k: np.array(v) for k, v in devs.items()}
    tol_all = dict(TOL, O2_GO=TOL_OXYGEN, O2_G=TOL_OXYGEN)
    q95 = {k: float(np.quantile(np.where(np.isfinite(arr[k]), arr[k], 10.0), 0.95)) for k in keys}   # inf -> 10 dex
    cal_tol = {k: max(TOL[k], q95[k]) for k in JOINT}
    return {"reps": reps, "truth_pass_rate": {k: float(np.mean(arr[k] <= tol_all[k])) for k in keys},
            "truth_joint_rate": float(np.mean(np.all([arr[k] <= TOL[k] for k in JOINT], axis=0))),
            "q95": q95, "calibrated_tol": cal_tol,
            "truth_joint_rate_calibrated": float(np.mean(np.all([arr[k] <= cal_tol[k] for k in JOINT], axis=0))),
            "O3_exact_rate": float(np.mean(o3)), "O3_class_order_rate": float(np.mean(o3c)), "truth": truth}


# -------------------------------------------------------------------------------------------------- reporting

def _summ(x: np.ndarray) -> str:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return "n/a"
    q = np.quantile(x, [0.1, 0.5, 0.9])
    return f"{q[1]:.3g} [{q[0]:.3g}, {q[2]:.3g}]"


def _verdict_shares(r30: np.ndarray) -> Dict[str, float]:
    n = max(r30.size, 1)
    return {"R in [1/2, 2]": round(float(np.sum((r30 >= 0.5) & (r30 <= 2)) / n), 3),
            "R in (2, 10]": round(float(np.sum((r30 > 2) & (r30 <= 10)) / n), 3),
            "R < 1/2": round(float(np.sum(r30 < 0.5) / n), 3), "R > 10": round(float(np.sum(r30 > 10) / n), 3)}


def summarise(res: Dict[str, Any]) -> Dict[str, Any]:
    sig, p = res["sig"], res["pass"]
    j = res["joint_idx"]
    return {"draws": res["n"], "seed": res["seed"], "chance": res["chance"], "reach": res["reach"],
            "best_norm_violation_sampled": round(res["best_norm"], 3),
            "best_dev_per_signature": {k: round(float(np.min(sig[k])), 3) for k in JOINT},
            "pairs": {f"{a}+{b}": int(np.sum(p[a] & p[b])) for i, a in enumerate(JOINT) for b in JOINT[i + 1:]},
            "all_draws": {"R30": _summ(sig["R30"]), "R120_fig2": _summ(sig["R120_fig2"]),
                          "G5_at10": _summ(sig["G5_at10"]), "G5_at500": _summ(sig["G5_at500"]),
                          "late_gain_500_240": _summ(sig["late_gain_500_240_model"]),
                          "arp_share_PA_120": _summ(sig["arp_share_PA_120"]),
                          "verdict_shares_R30": _verdict_shares(sig["R30"])},
            "joint_passes": {"n": int(j.size), "R30": _summ(sig["R30"][j]) if j.size else "n/a",
                             "verdict_shares_R30": _verdict_shares(sig["R30"][j]) if j.size else {}},
            "O2_GO_dev": _summ(sig["O2_GO"]), "O2_G_dev": _summ(sig["O2_G"]),
            "O3_exact_rate": float(np.mean(sig["O3_exact"])),
            "O3_class_order_rate": float(np.mean(sig["O3_class_order"])),
            "G3_beats_null_a_rate": float(np.mean(sig["G3_beats_null_a"]))}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--draws", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--oracle-reps", type=int, default=300)
    parser.add_argument("--recovery-reps", type=int, default=6)
    parser.add_argument("--search-maxiter", type=int, default=60)
    parser.add_argument("--scrambled-searches", type=int, default=5)
    parser.add_argument("--interp", choices=("linear", "early_step"), default="linear")
    parser.add_argument("--o2-high", type=float, default=1.0, help="upper edge of the glucose pots' O2 band")
    parser.add_argument("--tag", default="")
    parser.add_argument("--out-dir", type=Path, default=data_paths.VALIDATION_DIR)
    parser.add_argument("--blind-fig1", action="store_true", help="score the recorded DC_ARP optima on hofmann2000b "
                        "Fig. 1, the blind test of prereg sec. 12 (no search, no new constants)")
    parser.add_argument("--optima-dir", type=Path, default=data_paths.VALIDATION_DIR,
                        help="where the EX-049 feasibility records (primary + early_step) are read from")
    args = parser.parse_args(argv)
    if args.blind_fig1:
        return blind_main(args)
    data = load_data()
    meas = measured_signatures(data)
    o2_band = (0.0, args.o2_high)
    cfgs = {a: config(a, **{"strecker.forcing.interp": args.interp}) for a in ARMS}
    sk = dict(o2_band=o2_band, maxiter=args.search_maxiter)

    arms = {a: run_arm(cfgs[a], args.draws, args.seed + i, data, TOL, o2_band=o2_band) for i, a in enumerate(ARMS)}
    searches = {a: search(cfgs[a], data, TOL, args.seed + 50 + i, **sk) for i, a in enumerate(ARMS)}
    tol_no10 = dict(TOL, G5_no10=TOL["G5"])
    keys_no10 = tuple("G5_no10" if k == "G5" else k for k in JOINT)
    searches_no10 = {a: search(cfgs[a], data, tol_no10, args.seed + 60 + i, keys=keys_no10, **sk)
                     for i, a in enumerate(ARMS)}
    chance = {c: run_arm(cfgs["DC_ARP"], args.draws if c == "uninformed" else min(args.draws, 2000),
                         args.seed + 100 + i, data, TOL, chance=c, o2_band=o2_band)
              for i, c in enumerate(("uninformed", "scrambled", "flat"))}
    rng = np.random.default_rng(args.seed + 3)
    chance_search = {"scrambled": [search(cfgs["DC_ARP"], data, TOL, args.seed + 200 + i,
                                          courses=_chance_courses(data, rng, "scrambled", args.interp), **sk)
                                   ["min_violation"] for i in range(args.scrambled_searches)],
                     "flat": [search(cfgs["DC_ARP"], data, TOL, args.seed + 300,
                                     courses=_chance_courses(data, rng, "flat", args.interp), **sk)["min_violation"]]}

    def truth_theta(a):   # the search optimum: a joint pass if one exists, else the least-violating point
        return {k: np.array([v]) for k, v in searches[a]["theta"].items()}

    oracles = {}
    for a in ("DC_ARP", "DC"):
        for label, noise in (("realistic", NOISE_REALISTIC), ("stress_0.15", NOISE_STRESS)):
            o = oracle(truth_theta(a), data, cfgs[a], args.oracle_reps, args.seed + 7, noise)
            o["truth_violation_on_real_data"] = searches[a]["min_violation"]
            oracles[f"{a}|{label}"] = o
    cal_tol = oracles["DC_ARP|realistic"]["calibrated_tol"]
    arms_cal = {a: {"reach": {k: int((r["sig"][k] <= cal_tol[k]).sum()) for k in JOINT},
                    "joint": int(np.all([r["sig"][k] <= cal_tol[k] for k in JOINT], axis=0).sum())}
                for a, r in arms.items()}
    chance_cal = {c: {"rate": {k: round(float((r["sig"][k] <= cal_tol[k]).mean()), 3) for k in JOINT},
                      "joint_rate": float(np.all([r["sig"][k] <= cal_tol[k] for k in JOINT], axis=0).mean())}
                  for c, r in chance.items()}
    searches_cal = {a: search(cfgs[a], data, cal_tol, args.seed + 70 + i, **sk) for i, a in enumerate(ARMS)}
    searches_cal_no10 = {a: search(cfgs[a], data, dict(cal_tol, G5_no10=cal_tol["G5"]), args.seed + 80 + i,
                                   keys=keys_no10, **sk) for i, a in enumerate(ARMS)}

    # model recovery (positive control of the SEARCH): synthetic data from each truth, the search of every arm
    recovery: Dict[str, Dict[str, List[float]]] = {}
    rrng = np.random.default_rng(args.seed + 11)
    for ta in ("DC_ARP", "DC"):
        truth = truth_rows(truth_theta(ta), data, cfgs[ta])
        recovery[ta] = {a: [] for a in ARMS}
        for rep in range(args.recovery_reps):
            obs = noisy(truth, rrng, NOISE_REALISTIC)
            for i, a in enumerate(ARMS):
                recovery[ta][a].append(round(search(cfgs[a], obs, cal_tol, args.seed + 1000 + 10 * rep + i,
                                                    **sk)["min_violation"], 3))

    stem = args.out_dir / ("kinetic_core_strecker_forced_input_feasibility" + (f"_{args.tag}" if args.tag else ""))
    payload = {
        "artifact": stem.name,
        "provenance": provenance.provenance_block(
            stem.name, generated_by="scripts/generators/feasibility_strecker_forced_input.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in TABLES]),
        "design": "results/validation/kinetic_core_strecker_forced_input_prereg.md (review section, EX-042)",
        "single_laboratory": "every row is DFA Garching (Hofmann / Schieberle): no leave-one-lab-out is possible",
        "interp": args.interp, "o2_band": list(o2_band), "tolerances": TOL, "tol_oxygen": TOL_OXYGEN,
        "config_keys": [[c.name, c.default, c.band, c.phase, c.basis] for c in CONFIG_KEYS],
        "measured": meas,
        "arms_sampled": {a: summarise(r) for a, r in arms.items()},
        "arms_search": searches, "arms_search_G5_without_10min": searches_no10,
        "calibrated_tol": cal_tol, "arms_sampled_at_calibrated_tol": arms_cal,
        "arms_search_at_calibrated_tol": searches_cal, "arms_search_at_calibrated_tol_G5_without_10min":
            searches_cal_no10,
        "chance_sampled": {c: summarise(r) for c, r in chance.items()}, "chance_sampled_at_calibrated_tol": chance_cal,
        "chance_search_min_violation": chance_search,
        "oracle": {k: {kk: vv for kk, vv in o.items() if kk != "truth"} for k, o in oracles.items()},
        "model_recovery_min_violation": recovery,
        "negative_control_free_amine": free_amine_control(data, 2000, args.seed + 5),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=float) + "\n")
    stem.with_suffix(".md").write_text(render_md(payload, args) + "\n")
    print(json.dumps({a: (r["reach"]["joint"], round(searches[a]["min_violation"], 3)) for a, r in arms.items()}))
    print(f"wrote {stem.with_suffix('.md')}")
    return 0


def render_md(p: Dict[str, Any], args) -> str:
    L = ["# Strecker forced-input layer: feasibility (mechanisms before constants)", "",
         f"*Generated by `scripts/generators/feasibility_strecker_forced_input.py --draws {args.draws} --seed {args.seed}"
         f" --interp {args.interp}`. Design: {p['design']}. SINGLE LABORATORY: {p['single_laboratory']}.*", "",
         "## Reach per signature, uniform draws (prereg tolerances)", "",
         "| structure | " + " | ".join(JOINT) + " | joint | search min violation |", "|---" * (len(JOINT) + 3) + "|"]
    for a, s in p["arms_sampled"].items():
        L.append(f"| {a} | " + " | ".join(str(s["reach"][k]) for k in JOINT) +
                 f" | {s['reach']['joint']} / {s['draws']} | {p['arms_search'][a]['min_violation']:.2f} |")
    for c, s in p["chance_sampled"].items():
        L.append(f"| chance:{c} | " + " | ".join(str(s["reach"][k]) for k in JOINT) +
                 f" | {s['reach']['joint']} / {s['draws']} | "
                 f"{', '.join(f'{x:.2f}' for x in p['chance_search_min_violation'].get(c, [])) or '-'} |")
    L += ["", "## Search optimum per arm (deviations, dex)", "", "| arm | " + " | ".join(JOINT) +
          " | G5 at 10 | G5 at 500 | R30 | ARP share PA 120 |", "|---" * (len(JOINT) + 5) + "|"]
    for a, s in p["arms_search"].items():
        d = s["dev"]
        L.append(f"| {a} | " + " | ".join(f"{d[k]:.3f}" for k in JOINT) +
                 f" | {d['G5_at10']:+.3f} | {d['G5_at500']:+.3f} | {d['R30']:.2f} | {d['arp_share_PA_120']:.2f} |")
    L += ["", "G5 without its 10-min point (post hoc): " + json.dumps(
        {a: round(s["min_violation"], 3) for a, s in p["arms_search_G5_without_10min"].items()}), "",
          "## Calibrated tolerances (oracle, DC_ARP truth, realistic noise; max(prereg, q95))", "",
          json.dumps({k: round(v, 3) for k, v in p["calibrated_tol"].items()}), "",
          "search min violation at calibrated tol: " + json.dumps(
              {a: round(s["min_violation"], 3) for a, s in p["arms_search_at_calibrated_tol"].items()}),
          "; without the 10-min point: " + json.dumps(
              {a: round(s["min_violation"], 3) for a, s in p["arms_search_at_calibrated_tol_G5_without_10min"].items()}),
          "", "## Oracle", ""]
    for k, o in p["oracle"].items():
        L.append(f"- {k}: truth joint {o['truth_joint_rate']:.2f} (prereg tol), {o['truth_joint_rate_calibrated']:.2f}"
                 f" (calibrated); per signature {json.dumps({kk: round(v, 2) for kk, v in o['truth_pass_rate'].items()})};"
                 f" q95 {json.dumps({kk: round(v, 3) for kk, v in o['q95'].items()})}; O3 exact {o['O3_exact_rate']:.2f},"
                 f" class order {o['O3_class_order_rate']:.2f}")
    L += ["", "## Model recovery (search min violation per synthetic set, calibrated tolerances)", "",
          json.dumps(p["model_recovery_min_violation"]), "", "## Negative control: slope-1 free-amine law", "",
          json.dumps(p["negative_control_free_amine"])]
    return "\n".join(L)


# ------------------------------------------------------------------------------- blind test: hofmann2000b Fig. 1

#: prereg sec. 12 (declared before Fig. 1 was read): frozen calibrated tolerances; rule 6 drops a blind point whose
#: adjudicated read error exceeds half its tolerance; O2 at its band ceiling (rule 5); the forcing unchanged.
BLIND_TABLE = "strecker_hofmann2000b_fig1"
BLIND_POTS = {"glc": "hofmann2000b_fig1_glc_phe_reflux", "arp": "hofmann2000b_fig1_arp_phe_reflux"}
BLIND_TOL = {"G5p": 0.154, "G6a_p": 0.160, "G6b_p": math.log10(2.0), "G7": 0.154}
BLIND_TIMES = (10.0, 30.0, 60.0, 240.0, 300.0, 500.0)
BLIND_REF_T = 120.0
BLIND_O2 = 1.0
EARLY_T = 30.0
FEASIBILITY_RECORDS = ("kinetic_core_strecker_forced_input_feasibility.json",
                       "kinetic_core_strecker_forced_input_feasibility_early_step.json")


def load_fig1() -> Dict[str, Dict[str, Dict[float, Tuple[float, float]]]]:
    """{'glc' | 'arp': {PHEAL | PHEAC: {t: (mmol/L, read error dex)}}} from the keyed, two-reader table."""
    _, rows = _read_table(BLIND_TABLE)
    back = {v: k for k, v in BLIND_POTS.items()}
    out: Dict[str, Dict[str, Dict[float, Tuple[float, float]]]] = {"glc": {}, "arp": {}}
    for r in rows:
        out[back[r["pot"]]].setdefault(r["species"], {})[float(r["t_min"])] = (float(r["value"]), float(r["sigma"]))
    return out


def _abs_err(value: float, dex: float) -> float:
    return value * dex / math.log10(math.e)


def fig1_measured(fig1: Mapping[str, Any]) -> Dict[str, Any]:
    """The measured blind signatures and each point's read error (dex), with rule 6 applied."""
    g, a = fig1["glc"], fig1["arp"]
    lg = math.log10(math.e)

    def total(t):
        pa, paa = g["PHEAL"][t], g["PHEAC"][t]
        v = pa[0] + paa[0]
        return v, lg * math.hypot(_abs_err(*pa), _abs_err(*paa)) / v

    p120, e120 = total(BLIND_REF_T)
    a120, ea120 = a["PHEAL"][BLIND_REF_T]
    g5, g7 = {}, {}
    for t in BLIND_TIMES:
        v, e = total(t)
        err = math.hypot(e, e120)
        g5[t] = {"measured": math.log10(v / p120), "read_error": err, "dropped": err > BLIND_TOL["G5p"] / 2}
        av, ae = a["PHEAL"][t]
        err = math.hypot(ae, ea120)
        g7[t] = {"measured": math.log10(av / a120), "read_error": err, "dropped": err > BLIND_TOL["G7"] / 2}
    lr = {t: math.log10(g["PHEAL"][t][0] / g["PHEAC"][t][0]) for t in sorted(g["PHEAL"])}
    times = sorted(lr)
    cross = math.inf
    for t0, t1 in zip(times, times[1:]):
        if lr[t0] > 0 >= lr[t1]:
            cross = t0 + (t1 - t0) * lr[t0] / (lr[t0] - lr[t1])
            break
    return {"G5p": g5, "G7": g7, "G6a_p": math.log10(g["PHEAC"][500.0][0] / g["PHEAL"][500.0][0]), "G6b_p": cross,
            "arp_pa_gt_paa": {t: a["PHEAL"][t][0] > a["PHEAC"][t][0] for t in sorted(a["PHEAL"])}}


def blind_predict(theta: Mapping[str, float], data: Data, interp: str = "linear") -> Dict[str, np.ndarray]:
    """One recorded optimum (O2 at its ceiling, rule 5): the forced Glc/Phe pot and the fed ARP-Phe pot on GRID.
    The fed ARP-Phe decays first order at the optimum's k_decay from its arp0; cal_arp (the FORCED course's
    calibration) does not apply to a fed charge."""
    th = {k: np.array([float(v)]) for k, v in theta.items()}
    th["o2"] = np.array([BLIND_O2])
    cfg = config("DC_ARP", **{"strecker.forcing.interp": interp})
    cal = calibrate(th, data)
    g = glucose_pot(th, cal, interpolate(data.forced_inputs, interp), cfg)
    k = rate_constants(cal, {"cfg": cfg, "o2": th["o2"]})
    arp = th["arp0"][0] * np.exp(-(10.0 ** th["log10_k_arp_decay"][0]) * GRID)[None, :]
    pa = cum_with_loss(np.asarray(k["r_arp_ox_aldehyde"] + k["r_arp_anox_aldehyde"])[:, None] * arp, k["r_pheal_loss"])
    paa = cum_trapz(np.asarray(k["r_arp_ox_acid"] + k["r_arp_anox_acid"])[:, None] * arp)
    return {"glc_PA": g["PA"][0], "glc_PAA": g["PAA"][0], "arp_PA": pa[0], "arp_PAA": paa[0]}


def blind_score(pred: Mapping[str, np.ndarray], meas: Mapping[str, Any]) -> Dict[str, Any]:
    """Signed deviations (model - measured, dex) at every blind point; in-sample G6a' / G6b'; the PA > PAA check."""
    tot = pred["glc_PA"] + pred["glc_PAA"]
    g5 = {t: float(np.log10(tot[idx(t)] / tot[idx(BLIND_REF_T)]) - m["measured"]) for t, m in meas["G5p"].items()}
    a = pred["arp_PA"]
    g7 = {t: float(np.log10(a[idx(t)] / a[idx(BLIND_REF_T)]) - m["measured"]) for t, m in meas["G7"].items()}
    cross = float(_crossover(pred["glc_PA"][None, :], pred["glc_PAA"][None, :])[0])
    g6b = abs(math.log10(max(cross, DT) / meas["G6b_p"])) if math.isfinite(cross) and math.isfinite(meas["G6b_p"]) \
        else math.inf
    sign = all(pred["arp_PA"][idx(t)] > pred["arp_PAA"][idx(t)] for t in meas["arp_pa_gt_paa"])
    kept = [("G5p", t) for t, m in meas["G5p"].items() if not m["dropped"]] + \
           [("G7", t) for t, m in meas["G7"].items() if not m["dropped"]]
    devs = {"G5p": g5, "G7": g7}
    blind_ok = all(abs(devs[s][t]) <= BLIND_TOL[s] for s, t in kept)
    early = [t for t, m in meas["G5p"].items() if t <= EARLY_T and not m["dropped"]]
    return {"G5p": g5, "G7": g7, "G6a_p": float(abs(np.log10(pred["glc_PAA"][idx(500)] / pred["glc_PA"][idx(500)])
                                                    - meas["G6a_p"])),
            "G6b_p": g6b, "crossover_min": cross, "G7_pa_gt_paa": sign, "blind_points_within_tol": blind_ok,
            "worst_blind": max((abs(devs[s][t]) / BLIND_TOL[s], f"{s}@{t:g}") for s, t in kept),
            "early_deficit": bool(early) and all(g5[t] < -BLIND_TOL["G5p"] for t in early)}


def blind_outcome(scores: Sequence[Mapping[str, Any]]) -> str:
    """Prereg sec. 12, made computable before scoring: PASS needs >= 3 of 4 optima within tolerance on every blind
    point; with five recorded optima, every 4-subset must meet it (at most one failure). EARLY DEFICIT REPRODUCED: every
    optimum under-predicts the retained Glc/Phe points at t <= 30 min beyond tolerance. Anything else: INCONCLUSIVE."""
    if len(scores) >= 4 and all(sum(1 for s in sub if s["blind_points_within_tol"]) >= 3
                                for sub in itertools.combinations(scores, 4)):
        return "PASS"
    if scores and all(s["early_deficit"] for s in scores):
        return "EARLY DEFICIT REPRODUCED"
    return "INCONCLUSIVE"


def recorded_optima(optima_dir: Path) -> List[Dict[str, Any]]:
    """The shape-feasible DC_ARP search optima of EX-049 (min violation <= 1), exactly as recorded."""
    out = []
    for name in FEASIBILITY_RECORDS:
        rec = json.loads((optima_dir / name).read_text())
        for key, block in rec.items():
            if key.startswith("arms_search") and isinstance(block, dict) and "DC_ARP" in block:
                s = block["DC_ARP"]
                if s["min_violation"] <= 1.0:
                    out.append({"record": name, "search": key, "interp_found": rec["interp"],
                                "min_violation": s["min_violation"], "theta": s["theta"]})
    return out


def blind_main(args) -> int:
    data = load_data()
    meas = fig1_measured(load_fig1())
    optima = recorded_optima(args.optima_dir)
    result: Dict[str, Any] = {}
    for interp in ("linear", "early_step"):
        scores = [dict(blind_score(blind_predict(o["theta"], data, interp), meas), **{k: o[k] for k in
                  ("record", "search", "interp_found", "min_violation")}) for o in optima]
        result[interp] = {"scores": scores, "outcome": blind_outcome(scores)}
    stem = args.out_dir / "kinetic_core_strecker_forced_input_feasibility_blind_fig1"
    payload = {
        "artifact": stem.name,
        "provenance": provenance.provenance_block(
            stem.name, generated_by="scripts/generators/feasibility_strecker_forced_input.py --blind-fig1",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{BLIND_TABLE}.yml"] + [args.optima_dir / n for n in FEASIBILITY_RECORDS]),
        "design": "results/validation/kinetic_core_strecker_forced_input_prereg.md sec. 12 (declared before Fig. 1 was "
                  "read); single laboratory",
        "tolerances": BLIND_TOL, "o2": BLIND_O2, "measured": meas, "optima": optima,
        "primary": result["linear"], "secondary_early_step_reported_only": result["early_step"],
        "outcome": result["linear"]["outcome"],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=float) + "\n")
    L = ["# Strecker forced-input layer: the blind test on hofmann2000b Fig. 1", "",
         "*Generated by `scripts/generators/feasibility_strecker_forced_input.py --blind-fig1`. Design: prereg sec. 12 "
         "(declared before Fig. 1 was read). SINGLE LABORATORY. O2 = 1 (rule 5); forcing unchanged; no new constants.*", "",
         f"**Outcome (primary, linear forcing): {payload['outcome']}**; early-step (reported only): "
         f"{result['early_step']['outcome']}.", "",
         "Blind points dropped by rule 6 (read error > half the tolerance): " + json.dumps(
             {s: [f"{t:g} min ({m['read_error']:.3f} dex)" for t, m in meas[s].items() if m["dropped"]]
              for s in ("G5p", "G7")}), ""]
    for interp in ("linear", "early_step"):
        L += [f"## {interp} forcing", "", "| optimum | G5' 10 / 30 / 60 / 240 / 300 / 500 | G7 10 / 30 / 60 / 240 / 300 / "
              "500 | G6a' | G6b' (cross, min) | blind ok | early deficit |", "|---|---|---|---|---|---|---|"]
        for s in result[interp]["scores"]:
            L.append(f"| {s['record'].split('feasibility')[-1] or '_'}{s['search']} | "
                     + " / ".join(f"{s['G5p'][t]:+.2f}" for t in BLIND_TIMES) + " | "
                     + " / ".join(f"{s['G7'][t]:+.2f}" for t in BLIND_TIMES)
                     + f" | {s['G6a_p']:.3f} | {s['G6b_p']:.3f} ({s['crossover_min']:.0f}) | "
                     f"{s['blind_points_within_tol']} | {s['early_deficit']} |")
        L.append("")
    stem.with_suffix(".md").write_text("\n".join(L) + "\n")
    print(f"blind Fig. 1: primary {payload['outcome']}; early-step {result['early_step']['outcome']}; "
          f"wrote {stem.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

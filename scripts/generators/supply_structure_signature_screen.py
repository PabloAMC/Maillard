#!/usr/bin/env python
"""
supply_structure_signature_screen.py -- the thiol supply structures screened on three mechanism signatures, plug-in,
no fits (2026-11-05; the M1 + M2 prereg v2, sec. 0B; owner rule "mechanisms before constants").

Each candidate structure for the cysteine / pentose reservoir (TTCA <-> Cys + PENT, TTCA -> CAMA) is run through the
engine's own route pots (``fit_m1_supply_combinations.run_table_pots``, inline pH, scoring rtol) with every constant
at M1's carried C1 optimum except the searched ones. A bounded least-squares search over the reservoir coordinates
(inside their bands) looks for ANY point that meets each signature's declared screen tolerance:

  (a) Cao 2017 blank, 90 C: K_app = [TTCA] / ([Cys][Xyl]) inside 45 - 75 L/mol and flat (range <= 0.10 dex) over
      0.5 - 5 h, while the measured pH drifts 4.67 -> 4.15 (cao2017_extraction.md, arithmetic 1).
  (b) Zhai 2021 Fig. 4b, TTCA alone at 100 C: the ln-fraction lost by 80 min at pH 8 at least half the measured
      multiple of that at pH 5.5 (measured 9.6x; zhai2021_extraction.md sec. 2 and 5).
  (c) Zhang 2025b blank, 120 C: TTCA within 20 % of 17.7 / 18.2 mmol/L at 30 / 50 min, CAMA within x2 of 1.61 and free
      Cys / Xyl within x2 of 0.80 / 1.22 mmol/L at 50 min (the TTCA plateau with free reactants near 1 mmol/L is the
      ~19 mmol/L 1 : 1 missing pool, zhang2025b_extraction.md sec. 3.4; CF-063, CF-076).

Structures are implemented by wrapping ``sulfur.sulfur_rate_constants_at`` in this process (inherited by the forked
workers); nothing in src/ is changed. No band is set from Zhang 2025b. Also written: two arithmetic bounds from the
measured data (the pH slope of Cao's K_app; the 1 : 1 drain ceiling from Cao's pool closure) and the engine's pH in
Cao's pot under C1 on the COOLED scale Cao measured (``observations.ph_cooled``, as build_submodel_probe scores it),
against the keyed readings. The engine's in-situ (90 C) pH is not comparable with Cao's readings (vials cooled under
tap water before the meter) and is recorded only as such.

  --run STRUCTURE --work-dir DIR   one structure's search (one part file; ~15 min at 3 workers)
  --report --work-dir DIR          assemble results/validation/kinetic_core_supply_structure_signature_screen.{json,md}

The record is pinned by the commit that wrote it (a pre-design feasibility record, like the feasibility_* records):
BLINDING -- the session that ran it has seen engine predictions on Zhang 2025b's blank pot, so that pot is development
data for any later supply design, whose decisive reservoir data must be others.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ARTIFACT = "results/validation/kinetic_core_supply_structure_signature_screen"
GENERATED_BY = "scripts/generators/supply_structure_signature_screen.py"
CARRIED_C1 = "results/validation/kinetic_core_m1_supply_combinations_carried_C1.json"
ROUTE_TABLE = "data/lit/fit_rows/route_submodel.yml"

R_KJ = 8.314462618e-3
EA_LUMPED = 64.0822869313293          # B60's Ea_lumped_formation (kJ/mol)
EA_DRAIN_START = 64.08                # the drain barrier's search start (B60's lumped value, rounded as first run)
PH_REF_IN_SITU = 3.88                 # the engine's in-situ pH in Cao's pot at 90 C (the S2 laws' reference)
T_REF_DRAIN_K = 363.15                # 90 C: the drain's band is set at Cao's temperature

# ---- measured inputs (keyed / dossier values; never engine numbers) --------------------------------------------
CAO_K_APP = (58.6, 64.6, 60.2, 59.0, 58.5, 61.1)          # L/mol at 0.5 / 1 / 2 / 3 / 4 / 5 h (dossier arithmetic 1)
CAO_PH = (4.67, 4.56, 4.49, 4.33, 4.26, 4.15)              # measured after cooling (dossier table; the slope's x)
#: the keyed cooled readings (ph_acid_submodel.yml, pot cao2017_blank, Fig. 4c, vector): t = 0 and 0.5 - 5 h
CAO_PH_KEYED = {0: 5.50, 30: 4.679, 60: 4.559, 120: 4.482, 180: 4.324, 240: 4.253, 300: 4.148}
CAO_K_READ_REL = 0.07                                      # +-7 % read-error propagation on K_app (dossier)
ZHAI_FRAC80 = {"5.5": 8.19 / 10.0, "8": 1.47 / 10.0}       # TTCA left at 80 min of 10 mmol/L (route_submodel.yml)
ZHANG = {"TTCA30": 17.7, "TTCA50": 18.2, "CAMA50": 1.61, "Cys50": 0.8016, "PENT50": 1.22}

# ---- the screen's declared tolerances -----------------------------------------------------------------------------
CAO_K_BAND = (45.0, 75.0)       # Cao's plausible K_app range (dossier)
CAO_FLAT_DEX = 0.10             # measured +-4 % spread (0.043 dex) + two read errors (2 x 0.029 dex)
ZHAI_RATIO_FRACTION = 0.5       # predicted loss multiple at least half the measured one
ZHANG_TTCA_REL = 0.2
ZHANG_FACTOR = 2.0

POTS: Dict[str, Tuple[float, ...]] = {
    "cao2017_blank_fig5": (30, 60, 120, 180, 240, 300),
    "zhai2021_ttca_ph55_100c_fig4b": (20, 40, 60, 80, 100, 120),
    "zhai2021_ttca_ph7_100c_fig4b": (20, 40, 60, 80, 100, 120),
    "zhai2021_ttca_ph8_100c_fig4b": (20, 40, 60, 80),
    "zhang2025b_blank": (5, 10, 15, 20, 25, 30, 40, 50),
}
#: keyed pot -> (the run pot that reproduces it, signature)
ROW_SOURCES = {"cao2017_blank_fig4": ("cao2017_blank_fig5", "a"), "cao2017_blank_fig5": ("cao2017_blank_fig5", "a"),
               "zhai2021_ttca_ph55_100c_fig4b": ("zhai2021_ttca_ph55_100c_fig4b", "b"),
               "zhai2021_ttca_ph7_100c_fig4b": ("zhai2021_ttca_ph7_100c_fig4b", "b"),
               "zhai2021_ttca_ph8_100c_fig4b": ("zhai2021_ttca_ph8_100c_fig4b", "b"),
               "zhang2025b_blank": ("zhang2025b_blank", "c"), "zhang2025b_blank_fig2gh": ("zhang2025b_blank", "c")}
SPECIES = ("Cys", "PENT", "TTCA", "CAMA")
SEARCHED_BASE = ("k_ttca_cys", "ea_form_kj_mol", "k_pent_caramel", "k_cys_thermal", "k_ttca_deg")
SEARCH_BUDGET = {"joint": 150, "a": 60, "b": 60, "c": 60}
DRAIN_CEILING_LOG10 = math.log10(3e-3)

#: name -> (role, extra coordinates {name: (start, lo, hi)}, description, band sources)
STRUCTURES: Dict[str, Dict[str, Any]] = {
    "S1": {"role": "candidate", "extra": {},
           "what": "M1's C1: K 0.060 L/mmol (Cao 2017), pH- and T-flat; the ring opening k_ttca_cys searched",
           "bands": "k_ttca_cys B60 +-3 dex (stage-1 convention); ea_form 40 - 135 (sugar_trunk class); the rest B60's bounds"},
    "S2pow_rel": {"role": "candidate", "extra": {"m": (0.5, 0.0, 1.0)},
                  "what": "K ~ 10^(m pH) about the engine's in-situ pH in Cao's pot (3.88), the slope carried on the ring opening",
                  "bands": "m 0 - 1 (one proton; v1's declared band)"},
    "S2k3": {"role": "candidate",
             "extra": {"pKa_T": (5.4, 5.2, 6.1), "log10_drain90": (-3.0, -6.0, DRAIN_CEILING_LOG10),
                       "Ea_drain": (EA_DRAIN_START, 40.0, 135.0)},
             "what": "CF-076's Kallen-shaped K(pH) ~ (1 + 10^(pKa_T - pH)) 10^pH on the ring opening, plus S3's drain",
             "bands": "pKa_T 5.2 - 6.1 (Radomski 1989, as CF-076); drain as S3"},
    "S3": {"role": "candidate",
           "extra": {"K": (0.060, 0.045, 0.075), "log10_drain90": (-3.0, -6.0, DRAIN_CEILING_LOG10),
                     "Ea_drain": (EA_DRAIN_START, 40.0, 135.0)},
           "what": "K pH-flat inside Cao's band plus an equal first-order drain on free Cys and free PENT (1 : 1 at equal "
                   "free levels; added to r_cys_thermal and r_pent_caramel)",
           "bands": "K 0.045 - 0.075 (Cao 2017); drain at 90 C <= 3e-3 /min (Cao's pool closure, see bounds); "
                    "barrier 40 - 135 (sugar_trunk class)"},
    "S3X": {"role": "diagnostic",
            "extra": {"K": (0.060, 0.045, 0.075), "log10_drain90": (-3.0, -6.0, DRAIN_CEILING_LOG10),
                      "Ea_drain": (EA_DRAIN_START, 40.0, 135.0), "exit_n": (0.5, 0.0, 1.0)},
            "what": "S3 plus a base-favoured TTCA exit: r_cama_form and r_ttca_deg x 10^(n (pH - 5.5))",
            "bands": "n 0 - 1 (declared; diagnostic only)"},
    "S2k3X": {"role": "diagnostic",
              "extra": {"pKa_T": (5.4, 5.2, 6.1), "log10_drain90": (-3.0, -6.0, DRAIN_CEILING_LOG10),
                        "Ea_drain": (EA_DRAIN_START, 40.0, 135.0), "exit_n": (0.5, 0.0, 1.0)},
              "what": "the Kallen-shaped K + drain + the base-favoured exit",
              "bands": "as S2k3 and S3X (diagnostic only)"},
}

FEASIBILITY_KEYS = ("a_level", "a_flat", "b_ph_order", "c_ttca", "c_cama", "c_free")


# ================================================================================================================
# Pure functions (unit-tested)
# ================================================================================================================

def cao_k_slope(k_app: Sequence[float] = CAO_K_APP, ph: Sequence[float] = CAO_PH,
                read_rel: float = CAO_K_READ_REL) -> Dict[str, float]:
    """The slope of log10 K_app on the measured pH and its 90 % interval (t, n - 2 degrees of freedom), from the
    residual scatter and the read error in quadrature: the bound on m in K ~ 10^(m pH) over Cao's drift."""
    from scipy.stats import t as student_t

    x, y = np.asarray(ph, float), np.log10(np.asarray(k_app, float))
    a = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(a, y, rcond=None)
    r = y - a @ coef
    dof = len(x) - 2
    var_resid = float((r @ r) / dof * np.linalg.inv(a.T @ a)[0, 0])
    sxx = float(((x - x.mean()) ** 2).sum())
    var_read = (math.log10(1 + read_rel)) ** 2 / sxx
    t = float(student_t.ppf(0.95, dof))
    half = t * math.sqrt(var_resid + var_read)
    return {"slope": float(coef[0]), "se_residual": math.sqrt(var_resid), "se_read": math.sqrt(var_read),
            "ci90_lo": float(coef[0]) - half, "ci90_hi": float(coef[0]) + half, "dof": dof}


def kallen_shift_dex(pka_t: float, ph_from: float, ph_to: float) -> float:
    """log10 K(ph_to) - log10 K(ph_from) under K ~ (1 + 10^(pKa_T - pH)) 10^pH (TTCA ring N protonated below pKa_T,
    the cysteine amine fully protonated): flat below pKa_T, slope 1 above."""
    g = lambda p: math.log10(1 + 10 ** (pka_t - p)) + p  # noqa: E731
    return g(ph_to) - g(ph_from)


def drain_ceiling_per_min(pool_tolerance_mm: float = 60.0, free_mm: float = 65.0, minutes: float = 300.0) -> float:
    """The largest first-order drain on a free reactant that Cao 2017's cysteine pool (Cys + TTCA + CAMA) allows: the
    pool closes within its over-closure (6 - 17 % of the 400 mmol/L charge; ~15 % = 60 mmol/L) over 5 h with free Cys at
    ~65 mmol/L (dossier arithmetic 2). 60 / (65 x 300) = 3.1e-3 /min at 90 C."""
    return pool_tolerance_mm / (free_mm * minutes)


def zhai_loss_ratio(frac8: float, frac55: float) -> float:
    """ln(fraction left at pH 8) / ln(fraction left at pH 5.5) at the same time: how many times more TTCA is lost.
    A fraction of 0 counts as 1e-12 left; a pH 5.5 fraction of 1 or more as 1 - 1e-12 (nothing lost)."""
    if not frac8 or not frac55:
        return float("nan")
    return math.log(max(frac8, 1e-12)) / math.log(min(frac55, 1 - 1e-12))


def feasibility(sig: Mapping[str, float]) -> Dict[str, bool]:
    """The declared screen tolerances of the three signatures (module docstring)."""
    lr_meas = zhai_loss_ratio(ZHAI_FRAC80["8"], ZHAI_FRAC80["5.5"])
    lr = sig.get("zhai_lossratio_8_over_55", float("nan"))
    zh = lambda k, ref: sig.get(k, 0.0) / ref  # noqa: E731
    return {
        "a_level": CAO_K_BAND[0] <= sig.get("cao_K_mean", 0.0) <= CAO_K_BAND[1],
        "a_flat": sig.get("cao_K_range_dex", 9.0) <= CAO_FLAT_DEX,
        "b_ph_order": bool(lr == lr and lr >= ZHAI_RATIO_FRACTION * lr_meas),
        "c_ttca": abs(zh("zh_TTCA50", ZHANG["TTCA50"]) - 1) <= ZHANG_TTCA_REL
                  and abs(zh("zh_TTCA30", ZHANG["TTCA30"]) - 1) <= ZHANG_TTCA_REL,
        "c_cama": 1 / ZHANG_FACTOR <= zh("zh_CAMA50", ZHANG["CAMA50"]) <= ZHANG_FACTOR,
        "c_free": (1 / ZHANG_FACTOR <= zh("zh_Cys50", ZHANG["Cys50"]) <= ZHANG_FACTOR
                   and 1 / ZHANG_FACTOR <= zh("zh_PENT50", ZHANG["PENT50"]) <= ZHANG_FACTOR),
    }


def signature_verdicts(flags: Mapping[str, bool]) -> Dict[str, bool]:
    """Per signature: (a) both parts, (b), (c) all three parts -- at ONE point."""
    return {"a": bool(flags.get("a_level") and flags.get("a_flat")), "b": bool(flags.get("b_ph_order")),
            "c": bool(flags.get("c_ttca") and flags.get("c_cama") and flags.get("c_free"))}


def law_shift_dex(mod: Mapping[str, Any], ph: float) -> float:
    """log10 K(ph) - log10 K(ref) for the structure's pH law (0 when it has none)."""
    law = mod.get("law")
    if law is None:
        return 0.0
    if law[0] == "pow":
        return float(law[1]) * (ph - mod["ph_ref"])
    return kallen_shift_dex(float(law[1]), mod["ph_ref"], ph)


def modify_rates(k: Dict[str, float], mod: Mapping[str, Any], temperature_k: float, ph: float) -> Dict[str, float]:
    """The structure's changes to the engine's rate constants at (T, in-situ pH). The order of operations is fixed
    (law, K scale, exit, drain) and each touches its own keys."""
    if not mod:
        return k
    g = law_shift_dex(mod, ph)
    if g:
        k["r_ttca_cys"] *= 10 ** (-g)          # the slope carried on the ring opening: K = k_f / k_r rises with g
    if "K_scale" in mod:
        k["r_ttca_form"] *= mod["K_scale"]
    if "exit_n" in mod:
        f = 10 ** (mod["exit_n"] * (ph - 5.5))
        k["r_cama_form"] *= f
        k["r_ttca_deg"] *= f
    if "drain90" in mod:
        kd = mod["drain90"] * math.exp(-mod.get("Ea_drain", EA_LUMPED) / R_KJ * (1.0 / temperature_k - 1.0 / T_REF_DRAIN_K))
        k["r_cys_thermal"] += kd
        k["r_pent_caramel"] += kd
    return k


def mod_of(structure: str, v: Mapping[str, float]) -> Dict[str, Any]:
    if structure == "S1":
        return {}
    if structure == "S2pow_rel":
        return {"law": ("pow", v["m"]), "ph_ref": PH_REF_IN_SITU}
    m: Dict[str, Any] = {}
    if "pKa_T" in v:
        m.update(law=("spec", v["pKa_T"]), ph_ref=PH_REF_IN_SITU)
    if "K" in v:
        m["K_scale"] = v["K"] / 0.060
    if "log10_drain90" in v:
        m.update(drain90=10 ** v["log10_drain90"], Ea_drain=v["Ea_drain"])
    if "exit_n" in v:
        m["exit_n"] = v["exit_n"]
    return m


# ================================================================================================================
# The engine part
# ================================================================================================================

class Engine:
    """M1's C1 at its carried optimum, the route pots, and the rate-constant wrapper."""

    def __init__(self):
        from scripts.generators import fit_m1_supply_combinations as M
        from src.kinetic_core import config, sulfur

        self.M, self.config, self.sulfur = M, config, sulfur
        self.mod: Dict[str, Any] = {}
        self.failed_pot_runs = 0      # machine-state failures (memory, timeouts) make a part unusable (M1's rule)
        original = sulfur.sulfur_rate_constants_at

        def wrapped(parameters, temperature_k, ph, *a, **kw):
            return modify_rates(dict(original(parameters, temperature_k, ph, *a, **kw)), self.mod,
                                float(temperature_k), float(ph))

        sulfur.sulfur_rate_constants_at = wrapped
        carried = json.loads((ROOT / CARRIED_C1).read_text(encoding="utf-8"))
        self.theta0 = dict(carried["carried_theta"])
        self.num = M.score_inline()
        self.obj = M.Objective(M.RIVALS["C1"], self.num, M.FitSide.off())
        self.jobs = [M.table_job("route_submodel", p, t) for p, t in POTS.items()]
        lo, hi = self.obj.bounds(tuple(self.obj.keys))
        self.bounds = {k: (float(lo[self.obj.keys.index(k)]), float(hi[self.obj.keys.index(k)])) for k in SEARCHED_BASE}

    def cao_ph(self) -> Dict[str, Any]:
        """C1 at its carried optimum on Cao's blank pH pot (ph_acid_submodel.yml, cao2017_blank): the cooled pH the
        probe scores (observations.ph_cooled via build_submodel_probe.predict) and, for the record only, the in-situ
        pH at 90 C."""
        from types import SimpleNamespace

        self.mod.clear()
        times = sorted(CAO_PH_KEYED)
        job = self.M.table_job("ph_acid_submodel", "cao2017_blank", times)
        out: Dict[str, Any] = {"times_min": times, "measured_cooled": [CAO_PH_KEYED[t] for t in times]}
        points = {"C1_carried": self.obj.point([self.theta0[k] for k in self.obj.keys]),   # the screen's held point
                  "route_on_B60": self.M.b60_point(route=True)}                          # B60's constants, route on
        for label, pt in points.items():
            cfg = self.M.config_of(pt, self.num, self.M.FitSide.off())
            fails: List[Dict[str, Any]] = []
            with self.config.applied(cfg):
                (task,), (res,) = self.M.run_table_pots(pt, [job], self.num, cfg, fails)
                cooled = {t: float(self.M.probe.predict(SimpleNamespace(species="pH", t_min=float(t)), task, res))
                          for t in times}
            in_situ = dict(zip(task["times"], map(float, res.get("ph_in_situ", []))))
            out[label] = {"engine_cooled": [cooled[t] for t in times],
                          "residual_cooled": [cooled[t] - CAO_PH_KEYED[t] for t in times],
                          "engine_in_situ_90C_not_comparable": [in_situ.get(float(t)) for t in times],
                          "failed_pot_runs": fails}
        return out

    def run(self, theta: Mapping[str, float], mod: Mapping[str, Any]) -> Dict[str, Any]:
        self.mod.clear()
        self.mod.update(mod)
        pt = self.obj.point([theta[k] for k in self.obj.keys])
        cfg = self.M.config_of(pt, self.num, self.M.FitSide.off())
        fails: List[Dict[str, Any]] = []
        with self.config.applied(cfg):
            tasks, res = self.M.run_table_pots(pt, self.jobs, self.num, cfg, fails)
        self.failed_pot_runs += len(fails)
        out: Dict[str, Any] = {}
        for (pot, times), task, r in zip(POTS.items(), tasks, res):
            if "failed" in r:
                out[pot] = None
                continue
            tt = task["times"]
            out[pot] = {float(t): {s: float(r["states"][tt.index(float(t))][s]) for s in SPECIES} for t in times}
            out[pot]["_ph"] = [float(p) for p in r.get("ph_in_situ", [])]
        return out


def signatures(o: Mapping[str, Any]) -> Dict[str, Any]:
    s: Dict[str, Any] = {}
    c = o.get("cao2017_blank_fig5")
    if c:
        times = POTS["cao2017_blank_fig5"]
        k = np.array([1000 * c[float(t)]["TTCA"] / max(c[float(t)]["Cys"] * c[float(t)]["PENT"], 1e-300) for t in times])
        s["cao_K_mean"] = float(np.exp(np.mean(np.log(k))))
        s["cao_K_range_dex"] = float(np.log10(k.max() / k.min()))
        s.update(cao_Cys300=c[300.0]["Cys"], cao_PENT300=c[300.0]["PENT"], cao_TTCA30=c[30.0]["TTCA"],
                 cao_TTCA300=c[300.0]["TTCA"], cao_CAMA300=c[300.0]["CAMA"])
    for ph, pot in (("5.5", "zhai2021_ttca_ph55_100c_fig4b"), ("7", "zhai2021_ttca_ph7_100c_fig4b"),
                    ("8", "zhai2021_ttca_ph8_100c_fig4b")):
        if o.get(pot):
            s[f"zhai_frac80_ph{ph}"] = o[pot][80.0]["TTCA"] / 10.0
    if "zhai_frac80_ph8" in s and "zhai_frac80_ph5.5" in s:
        s["zhai_lossratio_8_over_55"] = zhai_loss_ratio(s["zhai_frac80_ph8"], s["zhai_frac80_ph5.5"])
    z = o.get("zhang2025b_blank")
    if z:
        s.update(zh_TTCA30=z[30.0]["TTCA"], zh_TTCA50=z[50.0]["TTCA"], zh_CAMA50=z[50.0]["CAMA"],
                 zh_Cys50=z[50.0]["Cys"], zh_PENT50=z[50.0]["PENT"])
    return s


def measured_rows() -> List[Tuple[str, str, str, float, float, float]]:
    """(signature, run pot, species, t, log10 value, sigma dex) from the keyed route table: the read error and the
    printed SD in quadrature, floor 0.02 dex."""
    import yaml

    rows = yaml.safe_load((ROOT / ROUTE_TABLE).read_text(encoding="utf-8"))["rows"]
    out = []
    for r in rows:
        if r["pot"] in ROW_SOURCES and r["species"] in SPECIES:
            pot, sig = ROW_SOURCES[r["pot"]]
            if r["t_min"] not in POTS[pot]:
                continue
            sd = r.get("sd")
            s = math.sqrt(r.get("read_error", 0.0) ** 2 + ((sd / r["value"] / math.log(10)) ** 2 if sd else 0.0))
            out.append((sig, pot, r["species"], float(r["t_min"]), math.log10(r["value"]), max(s, 0.02)))
    return out


def run_structure(structure: str, work_dir: Path) -> Path:
    from scipy.optimize import least_squares

    spec = STRUCTURES[structure]
    eng = Engine()
    data = measured_rows()
    names = list(SEARCHED_BASE) + list(spec["extra"])
    x0 = np.array([eng.theta0[k] for k in SEARCHED_BASE] + [v[0] for v in spec["extra"].values()])
    lo = np.array([eng.bounds[k][0] for k in SEARCHED_BASE] + [v[1] for v in spec["extra"].values()])
    hi = np.array([eng.bounds[k][1] for k in SEARCHED_BASE] + [v[2] for v in spec["extra"].values()])
    log: List[Dict[str, Any]] = []
    count = {"n": 0}

    class Budget(Exception):
        pass

    def evaluate(x):
        v = dict(zip(names, map(float, x)))
        theta = dict(eng.theta0, **{k: v[k] for k in SEARCHED_BASE})
        o = eng.run(theta, mod_of(structure, v))
        s = signatures(o)
        f = feasibility(s)
        res: Dict[str, List[float]] = {}
        for sig, pot, sp, t, y, sd in data:
            p = o.get(pot)
            val = p[t][sp] if p else 0.0
            res.setdefault(sig, []).append((math.log10(max(val, 1e-6)) - y) / sd)
        log.append({"x": v, "sig": s, "feasible": f, "chi2": {k: float(np.sum(np.square(r))) for k, r in res.items()},
                    "n": {k: len(r) for k, r in res.items()}})
        count["n"] += 1
        return res

    def search(which: str) -> None:
        start = count["n"]

        def fun(x):
            if count["n"] - start >= SEARCH_BUDGET[which]:
                raise Budget
            r = evaluate(x)
            keys = ("a", "b", "c") if which == "joint" else (which,)
            return np.concatenate([np.asarray(r[k]) for k in keys])

        try:
            least_squares(fun, x0.copy(), bounds=(lo, hi), x_scale=np.maximum(hi - lo, 1e-9) / 4, diff_step=1e-2)
        except Budget:
            pass

    t0 = time.time()
    for which in ("joint", "a", "b", "c"):
        search(which)
    part = {"structure": structure, "names": names, "lo": lo.tolist(), "hi": hi.tolist(), "evaluations": len(log),
            "failed_pot_runs": eng.failed_pot_runs, "wall_seconds": round(time.time() - t0, 1), "log": log}
    work_dir.mkdir(parents=True, exist_ok=True)
    path = work_dir / f"part_{structure}.json"
    path.write_text(json.dumps(part, default=float), encoding="utf-8")
    print(structure, "evaluations", len(log), "failed pot runs", eng.failed_pot_runs, "s", part["wall_seconds"],
          flush=True)
    return path


# ================================================================================================================
# The report
# ================================================================================================================

def summarise(part: Mapping[str, Any]) -> Dict[str, Any]:
    log = part["log"]
    any_flag = {k: any(e["feasible"].get(k) for e in log) for k in FEASIBILITY_KEYS}
    any_sig = {k: any(signature_verdicts(e["feasible"])[k] for e in log) for k in ("a", "b", "c")}
    joint = any(all(signature_verdicts(e["feasible"]).values()) for e in log)
    best = {}
    for sig in ("a", "b", "c"):
        i = min(range(len(log)), key=lambda j: log[j]["chi2"][sig] / log[j]["n"][sig])
        best[sig] = {"chi2_per_point": log[i]["chi2"][sig] / log[i]["n"][sig], "x": log[i]["x"],
                     "sig": {k: v for k, v in log[i]["sig"].items() if k != "cao_ph_in_situ"}}  # older parts carried it
    lr = [e["sig"].get("zhai_lossratio_8_over_55") for e in log]
    lr = [v for v in lr if v is not None and v == v]
    c_ok = [e for e in log if e["feasible"].get("c_ttca") and e["feasible"].get("c_cama")]
    return {"evaluations": part["evaluations"], "failed_pot_runs": part["failed_pot_runs"],
            "wall_seconds": part["wall_seconds"],
            "any_point_flags": any_flag, "any_point_signature": any_sig, "joint_at_one_point": joint,
            "max_zhai_loss_ratio": max(lr) if lr else None,
            "min_free_cys50_with_ttca_and_cama_ok": min((e["sig"]["zh_Cys50"] for e in c_ok), default=None),
            "min_free_pent50_with_ttca_and_cama_ok": min((e["sig"]["zh_PENT50"] for e in c_ok), default=None),
            "best_by_signature": best}


def build_report(work_dir: Path) -> Dict[str, Any]:
    from src import provenance

    results = {}
    for name in STRUCTURES:
        p = work_dir / f"part_{name}.json"
        if not p.exists():
            raise SystemExit(f"missing {p}: run --run {name} first")
        part = json.loads(p.read_text(encoding="utf-8"))
        if part.get("failed_pot_runs") is None or part["failed_pot_runs"] > 0:
            raise SystemExit(f"{p}: failed_pot_runs = {part.get('failed_pot_runs')}: re-run --run {name} on a quiet "
                             "machine (M1's rule: a record with failed pot runs is not used)")
        results[name] = summarise(part)
    slope = cao_k_slope()
    record = {
        "provenance": provenance.provenance_block(f"{ARTIFACT}.json", generated_by=GENERATED_BY,
                                                  inputs=[ROOT / CARRIED_C1, ROOT / ROUTE_TABLE]),
        "purpose": "pre-design feasibility of the thiol supply structures on three mechanism signatures (M1 + M2 "
                   "prereg v2 sec. 0B); plug-in engine runs, no fits; decides only which structures may enter a supply design",
        "blinding": "The session that ran this screen has seen engine predictions on Zhang 2025b's blank pot (free Cys, "
                    "free Xyl, TTCA, CAMA at 5 - 50 min, for every screened structure). That pot is development data for "
                    "any later supply design; the decisive reservoir data of that design must be others.",
        "signatures": {
            "a": {"what": "Cao 2017 blank, 90 C, 400 + 400 mmol/L, 0.2 M phosphate pH 5.5: K_app = [TTCA]/([Cys][Xyl]) "
                          "at 0.5 - 5 h", "measured_K_app_L_per_mol": list(CAO_K_APP), "measured_pH": list(CAO_PH),
                  "source": "cao2017_extraction.md arithmetic 1 (digitised Figs. 4 - 5); route_submodel.yml",
                  "tolerance": f"mean in {CAO_K_BAND} L/mol and range <= {CAO_FLAT_DEX} dex"},
            "b": {"what": "Zhai 2021 Fig. 4b, TTCA 10 mmol/L alone, 100 C, unbuffered: fraction left at 80 min, pH 8 vs 5.5",
                  "measured_fraction_left_80min": ZHAI_FRAC80,
                  "measured_loss_ratio": zhai_loss_ratio(ZHAI_FRAC80["8"], ZHAI_FRAC80["5.5"]),
                  "source": "zhai2021_extraction.md secs. 2, 5; route_submodel.yml zhai2021_*_fig4b",
                  "tolerance": f"predicted loss ratio >= {ZHAI_RATIO_FRACTION} x measured"},
            "c": {"what": "Zhang 2025b blank, 120 C, 40 + 40 mmol/L, 0.2 M phosphate pH 6.5", "measured": ZHANG,
                  "source": "zhang2025b_extraction.md secs. 3.3 - 3.4; route_submodel.yml zhang2025b_blank(_fig2gh); "
                            "CF-063, CF-076",
                  "tolerance": f"TTCA 30 / 50 min within {ZHANG_TTCA_REL:.0%}; CAMA 50 and free Cys / Xyl 50 within "
                               f"x{ZHANG_FACTOR:g} -- all at one point"},
        },
        "bounds": {
            "cao_K_pH_slope": dict(slope, note="log10 K_app regressed on the measured pH (after cooling), residual and "
                                                "+-7 % read error in quadrature, t with 4 dof: m <= ci90_hi"),
            "kallen_shift_over_cao_drift_dex": {str(p): kallen_shift_dex(p, CAO_PH[-1], CAO_PH[0]) for p in (5.2, 5.4, 6.1, 6.24)},
            "s3_drain_ceiling_per_min_90C": drain_ceiling_per_min(),
            "s3_drain_note": "Cao 2017's cysteine pool (Cys + TTCA + CAMA) closes within its 6 - 17 % over-closure over 5 h "
                             "(dossier arithmetic 2): ~60 mmol/L / (65 mmol/L free Cys x 300 min); the band uses 3e-3 /min",
        },
        "cao_pH": dict(Engine().cao_ph(),
                       note="C1_carried is M1's carried C1 optimum (the point the screen holds); route_on_B60 is B60's "
                            "constants with the route on (the setting of the 2026-11-05 pot-input audit). On the cooled "
                            "scale Cao measured (tap-water cooling, then the meter) both start at the measured 5.50 and, "
                            "from 30 min, drift by about the measured amount while reading below it (the audit traced its "
                            "~0.2 to cysteine carboxyl carried with its amine declared destroyed, ph_state.AMINE_FATE_BASIS; "
                            "CF-084 explanation 1); the drift-size match rests on compensating errors (TTCA at 5 h far "
                            "below the measured 236 mmol/L). The in-situ (90 C) trajectory is not comparable with Cao's "
                            "readings and is listed only as such. The screen's pH laws act on the engine's in-situ pH, so "
                            "the slope bound on K comes from the measured drift, not from the engine."),
        "structures": {n: dict(role=STRUCTURES[n]["role"], what=STRUCTURES[n]["what"], bands=STRUCTURES[n]["bands"],
                               **results[n]) for n in STRUCTURES},
        "settings": {"pots": {k: list(v) for k, v in POTS.items()}, "searched_base": list(SEARCHED_BASE),
                     "search_budget": SEARCH_BUDGET, "held": "every other constant at M1's carried C1 optimum",
                     "numerics": "inline pH, scoring rtol (fit_m1_supply_combinations.score_inline)"},
    }
    return record


def render_md(rec: Mapping[str, Any]) -> str:
    sl = rec["bounds"]["cao_K_pH_slope"]
    st = rec["structures"]
    yes = lambda b: "**feasible**" if b else "no"  # noqa: E731
    lines = [
        "# Supply structures screened on three mechanism signatures (pre-design, no fits)",
        "",
        f"*Generated by `{GENERATED_BY}` from `{ARTIFACT}.json`; do not edit by hand. Pre-registration context: "
        "[kinetic_core_m1m2_supply_routes_refit_prereg.md](kinetic_core_m1m2_supply_routes_refit_prereg.md) sec. 0B "
        "(v2, PARKED). Pinned by the commit that wrote it.*",
        "",
        f"**Blinding.** {rec['blinding']}",
        "",
        "**Verdict: no candidate supply structure is feasible on all three signatures.** (b) needs a TTCA exit that is "
        "faster at high pH, which no candidate has; (c) is reached by no structure inside bands set from data other than "
        "Zhang 2025b.",
        "",
        "## The signatures",
        "",
        "| | data | measured | tolerance (declared for the screen) |",
        "|---|---|---|---|",
    ]
    for k, s in rec["signatures"].items():
        meas = s.get("measured_K_app_L_per_mol") or s.get("measured_loss_ratio") or s.get("measured")
        if isinstance(meas, float):
            meas = f"loss ratio {meas:.2f}"
        lines.append(f"| ({k}) | {s['what']} ({s['source']}) | {meas} | {s['tolerance']} |")
    lines += [
        "",
        "## Bounds from the measured data (arithmetic)",
        "",
        f"- **K's pH slope over Cao's drift:** log10 K_app on pH, slope {sl['slope']:+.3f}, 90 % interval "
        f"{sl['ci90_lo']:+.3f} to {sl['ci90_hi']:+.3f} (se {sl['se_residual']:.3f} residual, {sl['se_read']:.3f} read). "
        f"So m <= {sl['ci90_hi']:.2f} in K ~ 10^(m pH).",
        "- **The Kallen-shaped law over the same drift** (log10 K change, 4.15 -> 4.67): "
        + "; ".join(f"pKa_T {p}: {v:.3f}" for p, v in rec["bounds"]["kallen_shift_over_cao_drift_dex"].items()) + ".",
        f"- **S3's drain ceiling at 90 C:** {rec['bounds']['s3_drain_ceiling_per_min_90C']:.2e} /min "
        f"({rec['bounds']['s3_drain_note']}).",
        "- **Cao's pH on the cooled scale** (`observations.ph_cooled`, as `build_submodel_probe` scores it; minutes: "
        "measured / engine): "
        + " ".join(f"{label} -- " + "; ".join(
            f"{t}: {m:.3f} / {e:.3f}" for t, m, e in zip(rec["cao_pH"]["times_min"], rec["cao_pH"]["measured_cooled"],
                                                       rec["cao_pH"][label]["engine_cooled"])) + "."
            for label in ("C1_carried", "route_on_B60"))
        + f" {rec['cao_pH']['note']}",
        "",
        "## Results (any evaluated point inside the bands)",
        "",
        "| structure | role | (a) | (b) | (c) | all three at one point | max (b) loss ratio | free Cys / Xyl at 50 min when TTCA and CAMA pass (mmol/L) | evaluations |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for n, r in st.items():
        sig = r["any_point_signature"]
        fc, fp = r["min_free_cys50_with_ttca_and_cama_ok"], r["min_free_pent50_with_ttca_and_cama_ok"]
        free = f"{fc:.2f} / {fp:.2f}" if fc is not None else "never passes"
        lines.append(f"| {n} | {r['role']} | {yes(sig['a'])} | {yes(sig['b'])} | {yes(sig['c'])} | "
                     f"{yes(r['joint_at_one_point'])} | {r['max_zhai_loss_ratio']:.2f} | {free} | {r['evaluations']} |")
    lines += ["", "Structures and bands:", ""]
    for n, r in st.items():
        lines.append(f"- **{n}** ({r['role']}): {r['what']}. Bands: {r['bands']}.")
    lines += ["", "## Best point per signature (chi^2 per point against the keyed rows)", "",
              "| structure | (a) | (b) | (c) |", "|---|---|---|---|"]
    for n, r in st.items():
        b = r["best_by_signature"]
        lines.append(f"| {n} | {b['a']['chi2_per_point']:.2f} | {b['b']['chi2_per_point']:.2f} | {b['c']['chi2_per_point']:.2f} |")
    lines += [
        "",
        "## Caveats",
        "",
        "- About 40 other constants are held at M1's carried C1 optimum (fitted under C1); the searches are local and "
        "budgeted; the tolerances are a screen's. (c)'s failure is structural for a pH-flat K in Cao's band: an 18 mmol/L "
        "TTCA plateau with free reactants near 1 mmol/L needs K_app >= ~19 L/mmol at equilibrium, or a release so slow "
        "that formation (k_f = K k_r) is equally slow.",
        "- In the engine the S2 laws act on the in-situ pH, which is not the scale Cao measured; (a)'s flatness in the "
        "engine is therefore weak evidence about m, and the slope bound comes from the measured (cooled) drift (above).",
        "- S2pow_rel's closest approach to (c) (best chi^2 per point above) is at m = 1, which the measured slope bound "
        "excludes; a K rising with pH also moves (b) the wrong way (TTCA more stable at high pH).",
        "- The diagnostics are not candidates: their exit law is declared, not measured.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--run", choices=sorted(STRUCTURES))
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--work-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.run:
        run_structure(args.run, args.work_dir)
    if args.report:
        from src import provenance

        rec = build_report(args.work_dir)
        (ROOT / f"{ARTIFACT}.json").write_text(provenance.dumps_compact(rec) + "\n", encoding="utf-8")
        (ROOT / f"{ARTIFACT}.md").write_text(render_md(rec), encoding="utf-8")
        print("wrote", f"{ARTIFACT}.json", f"{ARTIFACT}.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

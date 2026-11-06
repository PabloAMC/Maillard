#!/usr/bin/env python
"""
fit_m1b_ph_balance.py -- plan step M1b PART (i): the pot's charge-balance pH. B0 (the shipped balance, its two
in-sample constants integrated), B1 (amine survival), B2 (acid lag), B3 (acid supply), B13 (B1 + B3), on CF-084's
controls, built to its pre-registration.

Pre-registration (the specification; a rule here that differs from it is a bug here):
results/validation/kinetic_core_m1b_ph_balance_pentose_steps_prereg.md, v3.6 DRAFT; part (i) freezes after its screen
and power check (secs. 3.5 - 3.6). The engine side is src/kinetic_core/ph_balance_rivals.py (every switch off by
default); the measurement model's pieces are src/kinetic_core/predictive_densities.py (M1.1). Structural template:
scripts/generators/fit_m1_supply_combinations.py (modes, freeze checks, run-once, failed-pot refusal, pot alarms).

Modes (in this order; sec. 8):
    --checks                         V1 - V7, V-density and the default's bit-identity against the base commit
    --screen                         the mechanism screen (sec. 3.5): >= 600 Sobol draws per rival, coverage rule
    --power                          the synthetic power check (sec. 7) on the screen's draws: no new engine run
    --score --freeze-commit HASH     the decisive comparison (secs. 3.7, 5): RUN ONCE, after the freeze only
    --fit RIVAL --start N --freeze-commit HASH   the stacking guard's refit of a favoured rival (sec. 3.8)
Every pot runs under a wall-clock alarm; every failed pot run is counted, and a part with a failed pot run is
refused (its record says so and nothing downstream reads it). Writes results/validation/
kinetic_core_m1b_part_i_{checks,screen,power,score}.{json,md} and kinetic_core_m1b_part_i_fit_<R>_s<N>.json.
Memory (the owner's rule, 2026-11-05): engine-heavy modes check free memory first (``--min-free-percent``, default
40) and refuse to start below it; run with MAILLARD_WORKERS <= 2, under nice.

DECLARED READINGS (where the prereg is not explicit; each is printed in the records, none is a band or a rule):
  * Signatures (sec. 3.5) are read on the decisive (forced) C-Hou run, with the unforced run printed beside them:
    a level or shape is reproduced when it lies inside its unit's two-sided 90 % band (Z_90 x sqrt(noise^2 +
    delta-pH^2)); sigma-Mey when every pot's censored probability P(|d + e| <= 0.2) exceeds 0.05; sigma-Zhou when the
    endpoint rise from label 6 to 8 is below B0-pinned's re-printed value (V6).
  * The power check's "best draw" maximises the decisive units' likelihood (C-Mey (Cys) + C-Hou + C-Zhou-ARP) at the
    measured data; the tie band of each unit is computed once per pair and truth, at the truth's expected data; the
    rule's required-agreement forms in the power check are those the screen's draws carry (omega 0.5x / 2x, r 0.02,
    the unforced C-Hou, B0-pinned); the engine-changing forms (log-uniform yield, delta_T 0, the added-sink branch,
    Sohn's phi_Cys, C1) need new draws and are scored only in --score (sec. 7: "no new engine run").
  * The guards' Holm step uses p_g = Phi(d_g / sd_rep,g) (one-sided: the guard favours b when Holm rejects and its
    reading is b).
  * C-Zhou: the keyed Fig. 2 values and their printed replicate SDs (ph_acid_submodel.yml; sec. 3.3's rule "the
    printed replicate SD where one is readable"), not the unit table's "read 0.06".
  * Meynier's glycine pots are run for every rival: B0's integrated lumped yield acidifies them through
    r_pent_thermal_dpo -> r_osone_decay_dpo, so the prereg's "~0 drift" for B0 / B1 / B2 is printed, not assumed.
  * Hou's ARP, pH and glycine (Fig. 3, hou2026_extraction.md) are not keyed in data/lit/fit_rows; they are carried
    here as the dossier prints them (``HOU``), flagged.
  * Numerics: LSODA, with ONE declared BDF re-integration of a pot whose LSODA run fails or times out
    (``integrate_with_fallback``; LSODA stalled deterministically at an ordinary B0 draw on Meynier Cys pH 5.5).
    Every fallback is listed in the record; a pot failing both is a failed pot run. V-solver measures the agreement.
  * C-Mey's replicate-drawn edge (0.2 + N(0, 2 r^2)) is floored at 1e-3 (a negative edge is a statement no drift
    can meet).
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# B60's objective and vector (its import is grandfathered: it sets ph.b51 to its shipped value only), and M1's
# runner helpers (no switch set at import).
import generate_kinetic_core_b60_fit as B60  # noqa: E402
import fit_m1_supply_combinations as M1  # noqa: E402
import build_submodel_probe as probe  # noqa: E402
from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import config, fit_rows, fitting, observations, parallel, ph_state, sulfur  # noqa: E402
from src.kinetic_core import ph_balance_rivals as pbr  # noqa: E402
from src.kinetic_core import predictive_densities as pdens  # noqa: E402
from src.kinetic_core.ea_class_prior import CARBOHYDRATE, class_values, empirical_quantile  # noqa: E402
from src.kinetic_core.ph_state import LABEL_TEMPERATURE_K, BufferSpec, PhDrift  # noqa: E402

B23 = B60.B23
CELSIUS = 273.15
V = data_paths.VALIDATION_DIR
GENERATED_BY = "scripts/generators/fit_m1b_ph_balance.py"
STEM = "kinetic_core_m1b_part_i"
PREREG = V / "kinetic_core_m1b_ph_balance_pentose_steps_prereg.md"
OUT = {name: V / f"{STEM}_{name}.json" for name in ("checks", "screen", "power", "score")}
M1_CHECKS = V / "kinetic_core_m1_supply_combinations_checks.json"

# ----------------------------------------------------------------------------------------------------------------
# The prereg's numbers (sec. 3.1 bands, sec. 3.3 error model, sec. 3.5 screen, sec. 7 power)
# ----------------------------------------------------------------------------------------------------------------
PKA_BAND = (7.07, 8.18)                       # B60's 7.07 to Mossine 1994's DFG 8.18
YIELD_BAND = (0.0, 1.0)                       # ACID_YIELD_BOUNDS, primary
YIELD_LOG10_BAND = (-4.0, 0.0)                # required agreement: log-uniform on [1e-4, 1]
PHI_BAND = (0.0, 1.0)
PHI_CYS_SOHN_BAND = (0.76, 1.0)               # v3.5, required agreement (Sohn & Ho 1995)
F_RET_BAND = (0.0, 1.0)
LOG10_K_LAG_BAND = (math.log10(0.0097), 1.0)  # Zhai 2023's 0.0097 /min to the numerical 10 /min
DELTA_T_DEX = 0.30                            # the declared transfer term (sec. 3.1.1)
Y_LOWER_PRINTED = {"y5": 0.20, "y6": 0.56, "y_hi": 0.80}   # Davidek 2006a Fig. 3, lower bounds
Y5_UPPER = 0.56                               # monotonicity: y5 <= the pH 6 value
LOG10_K_ISO_BAND = (math.log10(0.025), math.log10(0.24))
EA_ISO_KJ_MOL = 76.0                          # ea_dpo_aa_kj_mol, Martins 2005's hexose analogue
DAVIDEK_LABELS = (5.0, 6.0, 6.5)              # 25 C labels of ~0.16 M phosphate
DAVIDEK_PHOSPHATE_M = 0.16
DAVIDEK_T_C = 120.0
OMEGAS = {"primary": pdens.OMEGA_DEX, "omega_half": 0.5 * pdens.OMEGA_DEX, "omega_double": 2.0 * pdens.OMEGA_DEX,
          "omega_reset_0p390": pdens.OMEGA_DEX_RESET}
METER_R = 0.05                                # decided; sensitivity 0.02
METER_R_SENSITIVITY = 0.02
MEYNIER_BOUND = 0.2
EDGE_FLOOR = 1e-3                             # declared reading: the replicate-drawn edge is floored here
SCREEN_DRAWS = 640                            # >= 600 (sec. 3.5)
SCREEN_SEED = 20261105
COVERAGE_MIN_PASSING = 10
POWER_R = 200
TIE_R = 4000
POWER_SEED = 20261106
ADEQUACY_P = 0.05
POT_SECONDS = 120
ATOL = 1e-14
SCORE_RTOL = 1e-8                             # sec. 3.5: scoring rtol 1e-8, inline
MIN_FREE_PERCENT = 30.0          # the coordinator's gate after the owner's VM stopped (2026-11-05)
MAX_WORKERS = 4                  # the coordinator's cap (M3's power check runs alongside at <= 4)
PAIRS: Tuple[Tuple[str, str], ...] = (("B1", "B0"), ("B2", "B0"), ("B3", "B0"), ("B13", "B0"), ("B13", "B3"))
RIVAL_NAMES = ("B0", "B1", "B2", "B3", "B13")

# ----------------------------------------------------------------------------------------------------------------
# The data
# ----------------------------------------------------------------------------------------------------------------
#: Hou 2026 Fig. 3 (hou2026_extraction.md, the table under "Fig. 3, model pots, 120 C [F]"): GX-ARP 5 mM alone,
#: water, initial pH 7.0 (1 M NaOH), 120 C. NOT KEYED in data/lit/fit_rows (flagged); read errors +-0.1 mM ARP,
#: +-0.05 pH, +-2 mg/L glycine (dossier header).
HOU = {"times_min": (20.0, 40.0, 60.0, 80.0, 100.0, 120.0),
       "ph_cooled": (5.05, 4.70, 4.30, 3.97, 3.72, 3.55),
       "arp_mmol_l": (3.20, 2.18, 1.27, 1.10, 0.75, 0.68),
       "gly_mg_l": (54.0, 68.0, 81.0, 83.5, 86.5, 81.5),
       "arp0_mmol_l": 5.0, "label_ph": 7.0, "t_c": 120.0,
       "read_ph": 0.05, "read_arp_mmol_l": 0.1, "read_gly_mg_l": 2.0, "gly_g_mol": 75.07,
       "source": "data/lit/extraction_dossiers/hou2026_extraction.md, Fig. 3 table (not keyed)"}
MEYNIER_LABELS = (4.5, 5.0, 5.5, 6.0, 6.5)
MEYNIER_TABLE = "heldout_meynier1995_tab1_2_4"
ZHOU_TABLE = "ph_acid_submodel"
ZHOU_LABELS = (6.0, 7.0, 8.0)
#: B51's buffer of Meynier & Mottram 1995 (ph_state.PYROPHOSPHATE): probe.buffer_of would read "pyrophosphate" as
#: orthophosphate, so the spec is written out
PYRO = BufferSpec(kind="pyrophosphate", pyrophosphate_mol_l=0.2, declared=True,
                  source="Meynier & Mottram 1995 Methods: 0.2 M sodium pyrophosphate (heldout_meynier1995_tab1_2_4.yml)")
WATER = BufferSpec(kind="none", declared=True, source="water, pH set with NaOH (Hou 2026; Zhou 2023)")
OX0 = 1.0   # B23's OX_AMBIENT_MMOL_L, as B60's Zhou systems (reproduces CF-084's numbers; V6)


@dataclass(frozen=True)
class Pot:
    name: str
    unit: str                  # the unit it serves
    initial: Tuple[Tuple[str, float], ...]
    t_c: float
    label_ph: float
    times: Tuple[float, ...]   # from 0
    buffer: BufferSpec
    forced: bool = False

    def initial_map(self) -> Dict[str, float]:
        return dict(self.initial)


def meynier_conditions() -> Dict[str, Any]:
    """Meynier's glycine pots as keyed (their CONDITIONS only: the table's rows are TEST data and never read here);
    the cysteine pots are the same design with cysteine for glycine (Table 3, not keyed; CF-084's control)."""
    from src import data_access

    table = data_access.load_yaml(fit_rows.FIT_ROWS_DIR / f"{MEYNIER_TABLE}.yml")
    return {k: v for k, v in table["pots"].items() if k.startswith("meynier1995_gly_ribose_pH")}


def control_pots() -> List[Pot]:
    """The ~18 pots of one draw (sec. 13: 5 Meynier Cys, Hou forced + unforced, 6 Zhou; + the 5 glycine guards)."""
    gly = meynier_conditions()
    out: List[Pot] = []
    for lab in MEYNIER_LABELS:
        key = f"meynier1995_gly_ribose_pH{str(lab).replace('.', 'p')}"
        p = gly[key]
        amine = float(p["initial_mmol_L"]["Gly"])
        pent = float(p["initial_mmol_L"]["PENT"])
        t_c, minutes = float(p["t_c"]), float(p["minutes"])
        assert float(p["label_ph"]) == lab
        out.append(Pot(f"mey_cys_pH{lab:g}", "C-Mey-Cys", (("Cys", amine), ("PENT", pent), ("OX", OX0)), t_c, lab,
                       (0.0, minutes), PYRO))
        out.append(Pot(f"mey_gly_pH{lab:g}", "C-Mey-Gly", (("Gly", amine), ("PENT", pent), ("OX", OX0)), t_c, lab,
                       (0.0, minutes), PYRO))
    hou_t = (0.0,) + tuple(HOU["times_min"])
    for forced in (True, False):
        out.append(Pot("hou_forced" if forced else "hou_unforced", "C-Hou" if forced else "C-Hou-unforced",
                       (("ARP", HOU["arp0_mmol_l"]), ("OX", OX0)), HOU["t_c"], HOU["label_ph"], hou_t, WATER, forced))
    for lab in ZHOU_LABELS:
        out.append(Pot(f"zhou_arp_pH{lab:g}", "C-Zhou-ARP", (("ARP", 20.0), ("OX", OX0)), 120.0, lab, (0.0, 60.0),
                       WATER))
        out.append(Pot(f"zhou_cys_pH{lab:g}", "C-Zhou-Cys", (("ARP", 20.0), ("Cys", 20.0), ("OX", OX0)), 120.0, lab,
                       (0.0, 60.0), WATER))
    return out


def zhou_measured() -> Dict[str, Dict[str, float]]:
    """{pot name: {value, sd}} from the keyed Fig. 2 rows (ph_acid_submodel.yml)."""
    _pots, rows = fit_rows.load(ZHOU_TABLE)
    out = {}
    for r in rows:
        if r.species != "pH" or not r.pot.startswith("zhou2023_arp"):
            continue
        lab = float(r.pot.rsplit("_ph", 1)[1])
        name = f"zhou_{'cys' if '_cys_' in r.pot else 'arp'}_pH{lab:g}"
        sd = r.raw.get("sd")
        out[name] = {"value": float(r.value),
                     "sd": float(sd) if sd else max(float(r.raw.get("read_error") or 0.0), METER_R),
                     "sd_basis": "printed replicate SD (sec. 3.3)" if sd else "max(read, r)", "row": r.pot}
    return out


def hou_forcing() -> pbr.ArpForcing:
    return pbr.ArpForcing(enabled=True, times_min=(0.0,) + tuple(HOU["times_min"]),
                          arp_mmol_l=(HOU["arp0_mmol_l"],) + tuple(HOU["arp_mmol_l"]))


def hou_gly_ratio() -> Tuple[np.ndarray, np.ndarray]:
    """C-Hou-Gly: glycine released / ARP consumed (mol / mol) at the six printed times and its propagated read sd."""
    gly = np.array(HOU["gly_mg_l"]) / HOU["gly_g_mol"]
    used = HOU["arp0_mmol_l"] - np.array(HOU["arp_mmol_l"])
    ratio = gly / used
    sd = ratio * np.sqrt((HOU["read_gly_mg_l"] / np.array(HOU["gly_mg_l"])) ** 2
                         + (HOU["read_arp_mmol_l"] / used) ** 2)
    return ratio, sd


# ----------------------------------------------------------------------------------------------------------------
# B60's vector, the rivals' coordinates and their bands
# ----------------------------------------------------------------------------------------------------------------
_CACHE: Dict[str, Any] = {}


def b60() -> Tuple[np.ndarray, Dict[str, Any], PhDrift]:
    if "b60" not in _CACHE:
        x = M1.b60_vector()
        _CACHE["b60"] = (x, B23.build_parameters(x), B23.unpack(x)[3])
    return _CACHE["b60"]


def b0_pinned_drift() -> PhDrift:
    return b60()[2]


@dataclass(frozen=True)
class Coord:
    name: str
    lo: float
    hi: float
    scale: str = "linear"      # "linear" | "log10" (lo / hi are log10) | "ea_class" (the carbohydrate class)

    def value(self, u: float) -> float:
        if self.scale == "ea_class":
            return empirical_quantile(class_values()[CARBOHYDRATE], float(u))
        v = self.lo + float(u) * (self.hi - self.lo)
        return 10.0 ** v if self.scale == "log10" else v


#: the engine-changing required-agreement forms (sec. 5.1) and which rivals carry them
FORMS = ("primary", "yield_loguniform", "delta_t_0", "added_sink", "sohn_phi_cys", "supply_c1")
FORM_RIVALS = {"primary": RIVAL_NAMES, "yield_loguniform": ("B0", "B1", "B2"), "delta_t_0": ("B3", "B13"),
               "added_sink": ("B3", "B13"), "sohn_phi_cys": ("B1", "B13"), "supply_c1": RIVAL_NAMES}


def coordinates(rival: str, form: str = "primary") -> Tuple[Coord, ...]:
    """The rival's integrated coordinates under a form (sec. 3.1)."""
    pka = Coord("arp_amine_pka", *PKA_BAND)
    y = (Coord("acid_yield", *YIELD_LOG10_BAND, "log10") if form == "yield_loguniform"
         else Coord("acid_yield", *YIELD_BAND))
    phi_cys = Coord("phi_cys", *(PHI_CYS_SOHN_BAND if form == "sohn_phi_cys" else PHI_BAND))
    phi = (phi_cys, Coord("phi_arp", *PHI_BAND))
    widen = 1.0 if form == "delta_t_0" else 10.0 ** (-DELTA_T_DEX)
    supply = (Coord("y5", Y_LOWER_PRINTED["y5"] * widen, Y5_UPPER), Coord("y6", Y_LOWER_PRINTED["y6"] * widen, 1.0),
              Coord("y_hi", Y_LOWER_PRINTED["y_hi"] * widen, 1.0), Coord("k_iso", *LOG10_K_ISO_BAND, "log10"))
    if rival == "B0":
        return (y, pka)
    if rival == "B1":
        return (y, pka) + phi
    if rival == "B2":
        return (y, pka, Coord("f_ret", *F_RET_BAND), Coord("k_lag", *LOG10_K_LAG_BAND, "log10"),
                Coord("ea_lag", 0.0, 1.0, "ea_class"))
    if rival == "B3":
        return (pka,) + supply
    if rival == "B13":
        return (pka,) + supply + phi
    raise KeyError(rival)


def order_ok(values: Mapping[str, float]) -> bool:
    return "y5" not in values or values["y5"] <= values["y6"] <= values["y_hi"]


def draw_coordinates(rival: str, form: str, n: int, seed: int, skip: int = 0) -> List[Dict[str, float]]:
    """n Sobol draws of the rival's coordinates (scrambled, seeded), B3's order y5 <= y6 <= y_hi enforced by
    rejection (sec. 3.1: 'enforced per draw'); ``skip`` continues the same sequence (the coverage rule's doubling)."""
    from scipy.stats import qmc

    coords = coordinates(rival, form)
    sob = qmc.Sobol(len(coords), scramble=True, seed=seed)
    out: List[Dict[str, float]] = []
    taken = 0
    while len(out) < n + skip:
        for u in sob.random(256):
            v = {c.name: c.value(x) for c, x in zip(coords, u)}
            if order_ok(v):
                taken += 1
                out.append(v)
                if len(out) >= n + skip:
                    break
    return out[skip:]


def davidek_breakpoints() -> Dict[str, Any]:
    """V5 / sec. 3.1.1: Davidek's labels (25 C, ~0.16 M phosphate) as in-situ pH at 120 C by the balance's T-shift."""
    buf = BufferSpec(kind="phosphate", phosphate_mol_l=DAVIDEK_PHOSPHATE_M, declared=True,
                     source="Davidek 2006a: ~0.16 M phosphate (dossier)")
    drift = b0_pinned_drift()
    out = []
    for lab in DAVIDEK_LABELS:
        sid = ph_state.strong_ion_difference(lab, {}, drift, buf)
        hot = ph_state.solve_ph(DAVIDEK_T_C + CELSIUS, sid, ph_state.titratable_inventory({}, drift, buf))
        out.append({"label_ph": lab, "label_scale": "label_25C", "label_T_C": 25.0, "in_situ_ph": float(hot),
                    "in_situ_scale": "in_situ", "in_situ_T_C": DAVIDEK_T_C})
    return {"breakpoints": out, "in_situ": tuple(r["in_situ_ph"] for r in out)}


def rival_changes(rival: str, form: str, values: Mapping[str, float], forced: bool,
                  breakpoints: Sequence[float]) -> Tuple[Dict[str, Any], PhDrift]:
    """The Config changes and the balance constants of one draw of one rival on one pot."""
    changes: Dict[str, Any] = {"ph.coupling": "inline"}
    if rival == "B0_pinned":
        drift = b0_pinned_drift()
    else:
        yld = 0.0 if rival in ("B3", "B13") else float(values["acid_yield"])
        drift = PhDrift(acid_yield=yld, arp_amine_pka=float(values["arp_amine_pka"]))
    if rival in ("B1", "B13"):
        changes["ph.amine_survival"] = pbr.AmineSurvival(True, float(values["phi_cys"]), float(values["phi_arp"]))
    if rival == "B2":
        changes["ph.acid_lag"] = pbr.AcidLag(True, float(values["f_ret"]), float(values["k_lag"]),
                                             float(values["ea_lag"]))
    if rival in ("B3", "B13"):
        changes["ph.acid_supply"] = pbr.AcidSupply(
            True, float(values["y5"]), float(values["y6"]), float(values["y_hi"]), tuple(map(float, breakpoints)),
            float(values["k_iso"]), EA_ISO_KJ_MOL, "added_sink" if form == "added_sink" else "hold_total")
    if forced:
        changes["ph.arp_forcing"] = hou_forcing()
    if form == "supply_c1":
        changes["sulfur.cys_amadori_route"] = dataclasses.replace(
            M1.shipped_value("sulfur.cys_amadori_route"), enabled=True)
    return changes, drift


# ----------------------------------------------------------------------------------------------------------------
# One pot of one draw
# ----------------------------------------------------------------------------------------------------------------
PRODUCED_ACIDS = ("CBX", "ACID", "AA", "FA")
PRODUCED_BASES = ("NH4", "RET_AA")


def state_at(run, i: int) -> Dict[str, float]:
    """The network's state at grid index i plus a rival's pools (the charge balance reads both)."""
    st = {k: float(run.series(k)[i]) for k in sulfur.SULFUR_STATE_KEYS}
    pools = (run.metadata.get("ph_balance_rivals") or {}).get("pools") or {}
    st.update(pbr.balance_pools({k: v[i] for k, v in pools.items()}))
    return st


def _cooled_compiled(state: Mapping[str, float], sid: float, drift: PhDrift, spec: BufferSpec, guess: float) -> float:
    inv = ph_state.titratable_inventory(state, drift, spec)
    if ph_state_kernel_compiled():
        return float(ph_state.solve_ph_and_ionic_strength_compiled(LABEL_TEMPERATURE_K, sid, inv, guess)[0])
    return float(ph_state.solve_ph(LABEL_TEMPERATURE_K, sid, inv))


def ph_state_kernel_compiled() -> bool:
    from src.kinetic_core import ph_kernel

    return bool(ph_kernel.COMPILED)


def delta_ph(state: Mapping[str, float], sid: float, drift: PhDrift, spec: BufferSpec, central: float
             ) -> Dict[str, float]:
    """sec. 3.3: delta-pH by two exact re-solves per omega, every produced titratable pool scaled independently
    (``predictive_densities.ph_species_discrepancy``). The pools' present amounts: each is made by the network
    (never charged) and non-negative; FA and AA are fragmented by the trunk's steps, so their present amount is what
    the network made and kept."""
    produced = {k: float(max(state.get(k, 0.0), 0.0)) for k in PRODUCED_ACIDS + PRODUCED_BASES}

    def resolve(amounts):
        return _cooled_compiled({**state, **amounts}, sid, drift, spec, central)

    out = {}
    for name, om in OMEGAS.items():
        r = pdens.ph_species_discrepancy(resolve, produced, PRODUCED_ACIDS, omega_dex=om, scale="cooled")
        out[name] = float(r["delta_ph"])
        if name == "primary":
            out["corners_bracket_central"] = bool(r["corners_bracket_central"])
    return out


def acid_to_measured(state: Mapping[str, float], sid: float, drift: PhDrift, spec: BufferSpec, target: float
                     ) -> Optional[float]:
    """CF-084's remedy: the strong-acid equivalents (mmol/L) the rival's state needs to read ``target`` cooled."""
    from scipy.optimize import brentq

    def f(x):
        return float(observations.ph_cooled(state, sid - x, drift, spec)) - float(target)

    try:
        return float(brentq(f, -100.0, 100.0, xtol=1e-6))
    except ValueError:
        return None


def _pot_task(item) -> Dict[str, Any]:
    """One pot under one draw's Config (installed by the parallel runner): cooled pH at the unit's times with its
    delta-pH at every omega, and the extras the units print. A failure is returned, never raised. The acid
    equivalents to the measured pH (C-Hou's printed inversion) are computed only when the tag asks (best draws)."""
    pot_name, drift_values, tag = item
    want_acid_eq = len(tag) > 3 and bool(tag[3])
    pot = POTS_BY_NAME[pot_name]
    drift = PhDrift(**drift_values)
    t0 = time.time()
    fallback = None
    try:
        run = integrate_with_fallback(pot, drift)
        fallback = run.metadata.get("m1b_solver_fallback")
    except M1._Timeout as exc:
        return {"tag": tag, "pot": pot_name, "failed": str(exc), "timeout": True}
    except Exception as exc:  # noqa: BLE001 - counted by the caller
        return {"tag": tag, "pot": pot_name, "failed": f"{type(exc).__name__}: {exc}", "timeout": False}
    sid = float(run.metadata["strong_ion_difference_mmol_L"])
    res: Dict[str, Any] = {"tag": tag, "pot": pot_name, "ph": [], "delta": {k: [] for k in OMEGAS}, "brackets": [],
                           "arp": [], "ret_aa": [], "aa_plus_laa": [], "acid_eq": [], "own_acid": []}
    pools = (run.metadata.get("ph_balance_rivals") or {}).get("pools") or {}
    for i, t in enumerate(pot.times):
        if t == 0.0:
            continue
        st = state_at(run, i)
        ph = float(observations.ph_cooled(st, sid, drift, pot.buffer))
        res["ph"].append(ph)
        d = delta_ph(st, sid, drift, pot.buffer, ph)
        for k in OMEGAS:
            res["delta"][k].append(d[k])
        res["brackets"].append(d["corners_bracket_central"])
        res["arp"].append(st.get("ARP", 0.0))
        res["ret_aa"].append(float(pools["RET_AA"][i]) if "RET_AA" in pools else 0.0)
        res["aa_plus_laa"].append(st.get("AA", 0.0) + (float(pools["LAA"][i]) if "LAA" in pools else 0.0))
        if pot.unit == "C-Hou":
            j = HOU["times_min"].index(t)
            res["acid_eq"].append(acid_to_measured(st, sid, drift, pot.buffer, HOU["ph_cooled"][j])
                                  if want_acid_eq else None)
            res["own_acid"].append(st.get("CBX", 0.0) + drift.acid_yield * st.get("ACID", 0.0) + st.get("AA", 0.0)
                                   + st.get("FA", 0.0))
    diag = (run.metadata.get("ph_balance_rivals") or {}).get("diagnostics") or {}
    res["diag"] = {k: diag.get(k) for k in ("forced_arp_max_rel_error", "forced_loss_never_negative",
                                             "forced_unresolvable_calls", "dpo_total_rel_dev_max",
                                             "added_sink_cap_reached", "supply_share_min", "supply_share_max")}
    res["ph_fixed_point_converged"] = run.metadata.get("ph_fixed_point_converged")
    res["solver_fallback"] = fallback
    res["seconds"] = round(time.time() - t0, 3)
    return res


#: Declared numerics (2026-11-05, before any screen record): LSODA (the lane's default, M1's) stalls
#: deterministically at ordinary draws (B0, yield 0.395, pKa 7.63, Meynier Cys pH 5.5: > 600 s), where BDF solves
#: the same pot in ~3 s and Radau agrees with it to 3e-9 pH. A pot whose LSODA run fails or exceeds its alarm is
#: re-integrated ONCE with BDF under the same alarm; the fallback is recorded per pot and counted per part; a pot
#: that fails both is a failed pot run (the part refuses). --checks' V-solver measures LSODA vs BDF agreement.
FALLBACK_METHOD = "BDF"


def integrate_with_fallback(pot: Pot, drift: PhDrift, method: str = "LSODA"):
    def run(m):
        return M1.with_alarm(POT_SECONDS, lambda: sulfur.integrate_sulfur(
            b60()[1], pot.t_c + CELSIUS, pot.initial_map(), np.array(pot.times, dtype=float), ph=pot.label_ph,
            buffer_spec=pot.buffer, ph_drift=drift, rtol=SCORE_RTOL, atol=ATOL, method=m))

    try:
        return run(method)
    except Exception as first:  # noqa: BLE001 - one declared retry; a second failure propagates
        out = run(FALLBACK_METHOD)
        out.metadata["m1b_solver_fallback"] = {"from": method, "to": FALLBACK_METHOD,
                                               "reason": f"{type(first).__name__}: {first}"[:200]}
        return out


POTS: List[Pot] = []
POTS_BY_NAME: Dict[str, Pot] = {}


def _init_pots() -> None:
    if not POTS:
        POTS.extend(control_pots())
        POTS_BY_NAME.update({p.name: p for p in POTS})


#: the rivals whose pairs read the C-Mey (Gly) guard (sec. 3.7: B3 and B13, against B0 and B3); B1 and B2 never do
GLY_GUARD_RIVALS = ("B0", "B0_pinned", "B3", "B13")


def pots_for(form: str, rival: Optional[str] = None) -> List[Pot]:
    """C1 changes nothing in a cysteine-free pot (sec. 3.7 item 3): the supply form runs the cysteine pots only.
    The glycine guard pots run only for the rivals a pair with that guard reads (``GLY_GUARD_RIVALS``)."""
    _init_pots()
    pots = list(POTS)
    if form == "supply_c1":
        pots = [p for p in pots if "Cys" in p.initial_map()]
    if rival is not None and rival not in GLY_GUARD_RIVALS:
        pots = [p for p in pots if p.unit != "C-Mey-Gly"]
    return pots


def build_tasks(rival: str, form: str, draws: Sequence[Mapping[str, float]], breakpoints: Sequence[float],
                pots: Sequence[Pot], acid_eq: bool = False) -> List[parallel.Configured]:
    base = config.Config.shipped()
    jobs = []
    for s, values in enumerate(draws):
        for pot in pots:
            changes, drift = rival_changes(rival, form, values, pot.forced, breakpoints)
            cfg = base.replace(changes)
            jobs.append(parallel.Configured((pot.name, dataclasses.asdict(drift), (rival, form, s, acid_eq)), cfg))
    return jobs


def run_jobs(jobs: Sequence[parallel.Configured], failures: List[Dict[str, Any]], chunk: int = 360,
             n_workers: Optional[int] = None) -> List[Dict[str, Any]]:
    """Run in chunks (memory: a fork per chunk), counting failures where the results arrive."""
    out: List[Dict[str, Any]] = []
    for i in range(0, len(jobs), chunk):
        part = parallel.map_ordered(_pot_task, list(jobs[i:i + chunk]), n_workers=n_workers)
        for r in part:
            if "failed" in r:
                failures.append({"pot": r["pot"], "tag": list(r["tag"]), "reason": r["failed"],
                                 "timeout": bool(r.get("timeout"))})
        out.extend(part)
    return out


# ----------------------------------------------------------------------------------------------------------------
# An ensemble: one rival's predictions on every unit, one row per draw
# ----------------------------------------------------------------------------------------------------------------
@dataclass
class Ensemble:
    rival: str
    form: str
    draws: List[Dict[str, float]]
    f: Dict[str, np.ndarray] = field(default_factory=dict)          # unit -> (S, m) cooled pH (Mey: the drift)
    delta: Dict[str, Dict[str, np.ndarray]] = field(default_factory=dict)   # unit -> omega form -> (S, m)
    extra: Dict[str, np.ndarray] = field(default_factory=dict)
    failed: int = 0
    fallbacks: List[Dict[str, Any]] = field(default_factory=list)   # pots solved by the declared BDF fallback

    @property
    def size(self) -> int:
        return len(self.draws)

    def subset(self, index: Sequence[int]) -> "Ensemble":
        ix = np.asarray(index, dtype=int)
        return Ensemble(self.rival, self.form, [self.draws[i] for i in ix],
                        {u: a[ix] for u, a in self.f.items()},
                        {u: {k: a[ix] for k, a in d.items()} for u, d in self.delta.items()},
                        {u: a[ix] for u, a in self.extra.items()}, self.failed, list(self.fallbacks))


UNIT_POTS = {"C-Mey-Cys": [f"mey_cys_pH{x:g}" for x in MEYNIER_LABELS],
             "C-Mey-Gly": [f"mey_gly_pH{x:g}" for x in MEYNIER_LABELS],
             "C-Hou": ["hou_forced"], "C-Hou-unforced": ["hou_unforced"],
             "C-Zhou-ARP": [f"zhou_arp_pH{x:g}" for x in ZHOU_LABELS],
             "C-Zhou-Cys": [f"zhou_cys_pH{x:g}" for x in ZHOU_LABELS]}


def assemble(rival: str, form: str, draws: List[Dict[str, float]], results: Sequence[Mapping[str, Any]]
             ) -> Ensemble:
    by = {(r["tag"][2], r["pot"]): r for r in results}
    ens = Ensemble(rival, form, list(draws))
    ens.failed = sum(1 for r in results if "failed" in r)
    ens.fallbacks = [{"pot": r["pot"], "coordinates": draws[r["tag"][2]], **r["solver_fallback"]} for r in results
                     if r.get("solver_fallback")]
    for unit, pots in UNIT_POTS.items():
        if not all((0, p) in by for p in pots):
            continue
        rows, deltas = [], {k: [] for k in OMEGAS}
        for s in range(len(draws)):
            vals, dl = [], {k: [] for k in OMEGAS}
            for p in pots:
                r = by.get((s, p))
                if r is None or "failed" in r:
                    vals.extend([np.nan] * (len(POTS_BY_NAME[p].times) - 1))
                    for k in OMEGAS:
                        dl[k].extend([np.nan] * (len(POTS_BY_NAME[p].times) - 1))
                    continue
                if unit.startswith("C-Mey"):
                    vals.extend([v - POTS_BY_NAME[p].label_ph for v in r["ph"]])
                else:
                    vals.extend(r["ph"])
                for k in OMEGAS:
                    dl[k].extend(r["delta"][k])
            rows.append(vals)
            for k in OMEGAS:
                deltas[k].append(dl[k])
        ens.f[unit] = np.array(rows, dtype=float)
        ens.delta[unit] = {k: np.array(v, dtype=float) for k, v in deltas.items()}
    for key, pot in (("hou_arp_forced", "hou_forced"), ("hou_arp_unforced", "hou_unforced"),
                     ("hou_ret_forced", "hou_forced")):
        if (0, pot) in by:
            field_name = "arp" if "arp" in key else "ret_aa"
            ens.extra[key] = np.array([by[(s, pot)].get(field_name, [np.nan] * 6) if (s, pot) in by else
                                       [np.nan] * 6 for s in range(len(draws))], dtype=float)
    if (0, "hou_forced") in by:
        ens.extra["hou_acid_eq"] = np.array([[np.nan if v is None else v for v in by[(s, "hou_forced")].get(
            "acid_eq", [None] * 6)] for s in range(len(draws))], dtype=float)
        ens.extra["hou_own_acid"] = np.array([by[(s, "hou_forced")].get("own_acid", [np.nan] * 6)
                                              for s in range(len(draws))], dtype=float)
    mey = [p for p in UNIT_POTS["C-Mey-Cys"] if (0, p) in by]
    if mey:
        ens.extra["mey_cys_aa_plus_laa_60"] = np.array([[by[(s, p)].get("aa_plus_laa", [np.nan])[-1] for p in mey]
                                                        for s in range(len(draws))], dtype=float)
    return ens


def ensemble_to_json(ens: Ensemble) -> Dict[str, Any]:
    return {"rival": ens.rival, "form": ens.form, "draws": ens.draws, "failed_pot_runs": ens.failed,
            "f": {u: a.tolist() for u, a in ens.f.items()},
            "delta": {u: {k: a.tolist() for k, a in d.items()} for u, d in ens.delta.items()},
            "extra": {u: a.tolist() for u, a in ens.extra.items()}, "solver_fallbacks": ens.fallbacks}


def ensemble_from_json(d: Mapping[str, Any]) -> Ensemble:
    return Ensemble(d["rival"], d["form"], list(d["draws"]),
                    {u: np.array(a, dtype=float) for u, a in d["f"].items()},
                    {u: {k: np.array(a, dtype=float) for k, a in dd.items()} for u, dd in d["delta"].items()},
                    {u: np.array(a, dtype=float) for u, a in d["extra"].items()}, int(d.get("failed_pot_runs", 0)),
                    list(d.get("solver_fallbacks", [])))


# ----------------------------------------------------------------------------------------------------------------
# The measurement model (sec. 3.3) and the predictive density per unit (sec. 5.1)
# ----------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class UnitData:
    """A unit's data: pH values with measurement sds, or Meynier's censored statement (``within`` per pot)."""

    unit: str
    y: Tuple[float, ...] = ()
    sd: Tuple[float, ...] = ()
    within: Tuple[bool, ...] = ()
    censored: bool = False


def measured_units(r_meter: float = METER_R) -> Dict[str, UnitData]:
    z = zhou_measured()
    hou_sd = tuple(max(HOU["read_ph"], r_meter) for _ in HOU["times_min"])
    out = {"C-Mey-Cys": UnitData("C-Mey-Cys", within=(True,) * 5, censored=True),
           "C-Mey-Gly": UnitData("C-Mey-Gly", within=(True,) * 5, censored=True),
           "C-Hou": UnitData("C-Hou", tuple(HOU["ph_cooled"]), hou_sd),
           "C-Hou-unforced": UnitData("C-Hou-unforced", tuple(HOU["ph_cooled"]), hou_sd)}
    for unit, kind in (("C-Zhou-ARP", "arp"), ("C-Zhou-Cys", "cys")):
        names = [f"zhou_{kind}_pH{x:g}" for x in ZHOU_LABELS]
        out[unit] = UnitData(unit, tuple(z[n]["value"] for n in names), tuple(z[n]["sd"] for n in names))
    return out


def unit_loglik(data: UnitData, f: np.ndarray, delta: np.ndarray, r_meter: float = METER_R,
                y: Optional[np.ndarray] = None, within: Optional[np.ndarray] = None,
                edge: Optional[np.ndarray] = None) -> np.ndarray:
    """log L_u per draw (S,), or per (replicate, draw) when y / within / edge carry a leading replicate axis.

    pH units: independent Gaussians, variance noise^2 + delta-pH^2 (V5: both sides 'cooled', checked by the caller).
    C-Mey: censored, P(|d + e| <= edge) with e ~ N(0, 2 r^2 + delta-pH^2) for a pot read within, 1 - that otherwise."""
    f = np.asarray(f, dtype=float)
    dl = np.asarray(delta, dtype=float)
    if data.censored:
        s = np.sqrt(2.0 * r_meter ** 2 + dl ** 2)                      # (S, m)
        e = np.full(f.shape[1], MEYNIER_BOUND) if edge is None else np.asarray(edge, dtype=float)
        w = np.asarray(data.within if within is None else within, dtype=bool)
        if e.ndim == 1 and w.ndim == 1:
            hi, lo = (e - f) / s, (-e - f) / s
            return np.sum(np.where(w[None, :], _log_inside(lo, hi), _log_outside(lo, hi)), axis=1)
        e2 = e if e.ndim == 2 else np.broadcast_to(e, (w.shape[0], f.shape[1]))
        w2 = w if w.ndim == 2 else np.broadcast_to(w, e2.shape)
        hi = (e2[:, None, :] - f[None]) / s[None]
        lo = (-e2[:, None, :] - f[None]) / s[None]
        return np.sum(np.where(w2[:, None, :], _log_inside(lo, hi), _log_outside(lo, hi)), axis=2)
    sd = np.asarray(data.sd, dtype=float)
    var = sd[None, :] ** 2 + dl ** 2
    yy = np.asarray(data.y if y is None else y, dtype=float)
    if yy.ndim == 1:
        r = yy[None, :] - f
        return -0.5 * np.sum(r ** 2 / var + np.log(2.0 * math.pi * var), axis=1)
    r = yy[:, None, :] - f[None]
    return -0.5 * np.sum(r ** 2 / var[None] + np.log(2.0 * math.pi * var)[None], axis=2)


def _log_inside(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """log(Phi(hi) - Phi(lo)) for lo < hi, without cancellation in either tail."""
    from scipy.special import log_ndtr

    with np.errstate(divide="ignore", invalid="ignore"):
        right = log_ndtr(-lo) + np.log1p(-np.exp(log_ndtr(-hi) - log_ndtr(-lo)))      # both in the upper tail
        left = log_ndtr(hi) + np.log1p(-np.exp(log_ndtr(lo) - log_ndtr(hi)))          # both in the lower tail
        mid = np.log(np.clip(np.exp(log_ndtr(hi)) - np.exp(log_ndtr(lo)), 1e-300, None))
    return np.where(lo > 0.0, right, np.where(hi < 0.0, left, mid))


def _log_outside(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """log(Phi(lo) + 1 - Phi(hi))."""
    from scipy.special import log_ndtr

    return np.logaddexp(log_ndtr(lo), log_ndtr(-hi))


def _logmeanexp(a: np.ndarray, axis: int = -1) -> np.ndarray:
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    return np.squeeze(m, axis=axis) + np.log(np.mean(np.exp(a - m), axis=axis))


@dataclass
class UnitScore:
    lppd: float
    accepted: bool
    reason: str
    k_hat: Optional[float]
    ess: Optional[float]
    mean: np.ndarray
    p_adequacy: float
    discrepancy_scale: float


def score_unit(data: UnitData, f: np.ndarray, delta: np.ndarray, r_meter: float = METER_R,
               y: Optional[np.ndarray] = None, within: Optional[np.ndarray] = None,
               precision_ess: float = 0.0) -> UnitScore:
    """lppd_u (sec. 5.1) of one rival's draw ensemble, its PSIS diagnostics, and the adequacy p (sec. 5.4)."""
    ll = unit_loglik(data, f, delta, r_meter, y=y, within=within)
    ok = np.isfinite(ll)
    if not ok.all():
        raise pdens.MeasurementModelRefusal("a non-finite likelihood: a failed pot run reached the scorer")
    if ll.size == 1:
        est = pdens.DensityEstimate(float(ll[0]), "plug_in", True)
    else:
        est = pdens.importance_lppd(ll, np.zeros_like(ll), precision_ess=precision_ess)
    w = np.exp(ll - np.max(ll))
    w = w / w.sum()
    mean = w @ f
    scale = float(w @ np.mean(delta, axis=1))
    if data.censored:
        # the censored likelihood's PPC: the posterior-predictive probability of the observed outcome pattern
        p = float(np.sum(w * np.exp(ll)))
    else:
        yy = np.asarray(data.y if y is None else y, dtype=float)
        var = np.asarray(data.sd, float) ** 2 + delta ** 2
        cov = np.diag(w @ var) + (f - mean).T @ ((f - mean) * w[:, None])
        p = pdens.posterior_predictive_p(yy, mean, cov)
    return UnitScore(float(est.lppd), bool(est.accepted), est.reason, est.k_hat, est.ess, mean, p, scale)


def replicate_data(data: UnitData, rng: np.random.Generator, n: int, r_meter: float = METER_R
                   ) -> Dict[str, np.ndarray]:
    """sec. 5.2: R redraws of y_u from its MEASUREMENT noise only; for C-Mey the bound's edge redrawn as 0.2 +
    N(0, 2 r^2) per replicate (the outcome kept)."""
    if data.censored:
        m = len(data.within)
        # a redrawn edge below zero (2.8 sd at r = 0.05) is a statement no drift can meet: floored at EDGE_FLOOR
        edge = MEYNIER_BOUND + rng.normal(0.0, math.sqrt(2.0) * r_meter, size=(n, m))
        return {"edge": np.maximum(edge, EDGE_FLOOR),
                "within": np.broadcast_to(np.asarray(data.within, bool), (n, m))}
    return {"y": np.asarray(data.y, float)[None, :] + rng.normal(size=(n, len(data.y))) * np.asarray(data.sd)[None]}


def tie_band(data: UnitData, ens_a: Tuple[np.ndarray, np.ndarray], ens_b: Tuple[np.ndarray, np.ndarray],
             rng: np.random.Generator, r_meter: float = METER_R, n: int = TIE_R) -> pdens.TieBand:
    """b_u = z_0.95 sd_rep(d_u) at fixed predictive distributions (the two ensembles)."""
    rep = replicate_data(data, rng, n, r_meter)
    d = []
    for i in range(0, n, 500):
        sl = slice(i, min(i + 500, n))
        kw = ({"within": rep["within"][sl], "edge": rep["edge"][sl]} if data.censored else {"y": rep["y"][sl]})
        la = _logmeanexp(unit_loglik(data, *ens_a, r_meter, **kw), axis=1)
        lb = _logmeanexp(unit_loglik(data, *ens_b, r_meter, **kw), axis=1)
        d.append(la - lb)
    dd = np.concatenate(d)
    if not np.all(np.isfinite(dd)):
        raise pdens.MeasurementModelRefusal("a replicate gave a non-finite d_u")
    sd = float(np.std(dd, ddof=1))
    return pdens.TieBand(sd, pdens.Z_90 * sd, int(dd.size))


def _censored_loglik_edge(data, f, delta, r_meter, edge):
    return unit_loglik(data, f, delta, r_meter, within=np.asarray(data.within, bool)[None, :],
                       edge=np.asarray(edge, float)[None, :])[0]


# ----------------------------------------------------------------------------------------------------------------
# The decision rule (sec. 3.7 item 1, sec. 5.3 - 5.5)
# ----------------------------------------------------------------------------------------------------------------
DECISIVE = ("C-Mey-Cys", "C-Hou", "C-Zhou-ARP")


def guards_of(a: str) -> Tuple[str, ...]:
    return ("C-Zhou-Cys", "C-Mey-Gly") if a in ("B3", "B13") else ("C-Zhou-Cys",)


@dataclass(frozen=True)
class RuleForm:
    """One required-agreement form of the rule that the draws already carry (likelihood-only, the unforced C-Hou,
    B0-pinned as b)."""

    name: str
    omega: str = "primary"
    r_meter: float = METER_R
    hou_unit: str = "C-Hou"
    b_pinned: bool = False


LIKELIHOOD_FORMS = (RuleForm("primary"), RuleForm("omega_half", omega="omega_half"),
                    RuleForm("omega_double", omega="omega_double"), RuleForm("r_0p02", r_meter=METER_R_SENSITIVITY),
                    RuleForm("hou_unforced", hou_unit="C-Hou-unforced"))


def _arrays(ens: Ensemble, unit: str, omega: str) -> Tuple[np.ndarray, np.ndarray]:
    return ens.f[unit], ens.delta[unit][omega]


def compare(a: str, ens_a: Ensemble, ens_b: Ensemble, data: Mapping[str, UnitData], form: RuleForm,
            rng: np.random.Generator, overrides: Optional[Mapping[str, Mapping[str, Any]]] = None,
            bands: Optional[Mapping[str, pdens.TieBand]] = None, require_resolved: bool = True) -> Dict[str, Any]:
    """The rule for 'a favoured over b' under one form: every unit's d_u, tie band, adequacy and reading; the guards'
    Holm step. ``overrides`` replaces a unit's data (synthetic y / within); ``bands`` reuses precomputed tie bands.
    ``require_resolved`` False keeps an unaccepted PSIS estimate (the power check's diagnostic reading)."""
    units = {}
    for u in DECISIVE + guards_of(a):
        unit = form.hou_unit if u == "C-Hou" else u
        dat = data[unit]
        ov = dict((overrides or {}).get(u) or {})
        fa, fb = _arrays(ens_a, unit, form.omega), _arrays(ens_b, unit, form.omega)
        band = (bands or {}).get(f"{form.name}|{u}") or tie_band(dat, fa, fb, rng, form.r_meter)
        prec = pdens.tie_band_precision_ess(band.sd_rep)
        sa = score_unit(dat, *fa, form.r_meter, precision_ess=prec, **ov)
        sb = score_unit(dat, *fb, form.r_meter, precision_ess=prec, **ov)
        d = sa.lppd - sb.lppd
        resolved = sa.accepted and sb.accepted
        if resolved or not require_resolved:
            reading = pdens.unit_reading(d, band.band, sa.p_adequacy, sb.p_adequacy, sa.mean, sb.mean,
                                         max(sa.discrepancy_scale, sb.discrepancy_scale), ADEQUACY_P)
        else:
            reading = "tie"
        units[u] = {"unit_data": unit, "d": d, "band": band.band, "sd_rep": band.sd_rep, "reading": reading,
                    "resolved": resolved, "reasons": [x for x in (sa.reason, sb.reason) if x],
                    "p_a": sa.p_adequacy, "p_b": sb.p_adequacy, "k_hat": [sa.k_hat, sb.k_hat], "ess": [sa.ess, sb.ess]}
    from scipy.stats import norm

    guard_p = {g: float(norm.cdf(units[g]["d"] / units[g]["sd_rep"])) if units[g]["sd_rep"] > 0 else 0.5
               for g in guards_of(a)}
    holm = pdens.holm(guard_p)
    guard_favours_b = {g: bool(holm[g]["rejected"] and units[g]["reading"] == "b") for g in guard_p}
    favoured = (units["C-Mey-Cys"]["reading"] == "a" and units["C-Hou"]["reading"] == "a"
                and units["C-Zhou-ARP"]["reading"] != "b" and not any(guard_favours_b.values()))
    both_inadequate = any(units[u]["reading"] == "both inadequate" for u in DECISIVE)
    return {"form": form.name, "units": units, "guards_holm": holm, "guard_favours_b": guard_favours_b,
            "favoured": bool(favoured), "both_inadequate": bool(both_inadequate)}


def pinned_ensemble(b0_pinned: Ensemble) -> Ensemble:
    return b0_pinned


def rule_over_forms(a: str, b: str, ens: Mapping[str, Ensemble], data: Mapping[str, UnitData],
                    rng: np.random.Generator, forms: Sequence[RuleForm] = LIKELIHOOD_FORMS,
                    overrides=None, bands=None, require_resolved: bool = True) -> Dict[str, Any]:
    """'a favoured over b' iff favoured under every form (sec. 5.3), plus B0-pinned as b when b is B0 (N2)."""
    per = {}
    fs = list(forms) + ([RuleForm("vs_b0_pinned", b_pinned=True)] if b == "B0" and "B0_pinned" in ens else [])
    for form in fs:
        eb = ens["B0_pinned"] if form.b_pinned else ens[b]
        per[form.name] = compare(a, ens[a], eb, data, form, rng, overrides, bands, require_resolved)
    favoured = all(v["favoured"] for v in per.values())
    return {"a": a, "b": b, "favoured": favoured, "per_form": per,
            "both_inadequate_primary": per["primary"]["both_inadequate"],
            "breaking_forms": [k for k, v in per.items() if not v["favoured"]]}


# ----------------------------------------------------------------------------------------------------------------
# The screen's signatures (sec. 3.5)
# ----------------------------------------------------------------------------------------------------------------
def signatures(ens: Ensemble, b0_spread: float, hou_unit: str = "C-Hou") -> Dict[str, np.ndarray]:
    """Per draw: sigma-Hou-level, sigma-Hou-shape, sigma-Mey, sigma-Zhou (declared readings: module docstring)."""
    data = measured_units()
    f, dl = ens.f[hou_unit], ens.delta[hou_unit]["primary"]
    sd = np.asarray(data[hou_unit].sd)
    level = np.abs(f[:, 0] - HOU["ph_cooled"][0]) <= pdens.Z_90 * np.sqrt(sd[0] ** 2 + dl[:, 0] ** 2)
    i60, i120 = HOU["times_min"].index(60.0), HOU["times_min"].index(120.0)
    fall_pred = f[:, i60] - f[:, i120]
    fall_meas = HOU["ph_cooled"][i60] - HOU["ph_cooled"][i120]
    shape = np.abs(fall_pred - fall_meas) <= pdens.Z_90 * np.sqrt(sd[i60] ** 2 + sd[i120] ** 2 + dl[:, i60] ** 2
                                                                   + dl[:, i120] ** 2)
    fm, dm = ens.f["C-Mey-Cys"], ens.delta["C-Mey-Cys"]["primary"]
    s = np.sqrt(2.0 * METER_R ** 2 + dm ** 2)
    from scipy.stats import norm

    p_in = norm.cdf((MEYNIER_BOUND - fm) / s) - norm.cdf((-MEYNIER_BOUND - fm) / s)
    mey = np.all(p_in > ADEQUACY_P, axis=1)
    fz = ens.f["C-Zhou-ARP"]
    zhou = (fz[:, -1] - fz[:, 0]) < b0_spread
    return {"sigma_hou_level": level, "sigma_hou_shape": shape, "sigma_mey": mey, "sigma_zhou": zhou,
            "all": level & shape & mey & zhou, "hou_level_ph20": f[:, 0], "hou_fall_60_120": fall_pred,
            "zhou_spread_6_8": fz[:, -1] - fz[:, 0], "mey_min_p_within": np.min(p_in, axis=1)}


def best_draw(ens: Ensemble, data: Mapping[str, UnitData]) -> int:
    ll = sum(unit_loglik(data[u], ens.f[u], ens.delta[u]["primary"]) for u in DECISIVE)
    return int(np.argmax(ll))


# ----------------------------------------------------------------------------------------------------------------
# Records
# ----------------------------------------------------------------------------------------------------------------
def _jsonable(x):
    return M1._jsonable(x)


def write_record(name: str, payload: Mapping[str, Any], md: str) -> None:
    path = OUT[name]
    path.write_text(provenance.dumps_compact(_jsonable(payload), default=float) + "\n", encoding="utf-8")
    path.with_suffix(".md").write_text(md, encoding="utf-8")
    print(f"wrote {data_paths.rel(path)} and .md", flush=True)


def read_record(name: str) -> Dict[str, Any]:
    p = OUT[name]
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _provenance(name: str) -> Dict[str, Any]:
    return provenance.provenance_block(
        f"{STEM}_{name}", generated_by=GENERATED_BY, wave="M1b part (i)",
        inputs=[PREREG, M1.B60_REPORT, fit_rows.FIT_ROWS_DIR / f"{ZHOU_TABLE}.yml",
                fit_rows.FIT_ROWS_DIR / f"{MEYNIER_TABLE}.yml", ROOT / "src/kinetic_core/ph_balance_rivals.py",
                ROOT / "src/kinetic_core/ph_state.py", ROOT / "src/kinetic_core/sulfur.py",
                ROOT / "src/kinetic_core/predictive_densities.py", ROOT / "src/kinetic_core/config.py",
                Path(__file__)])


def memory_gate(min_free_percent: float) -> Dict[str, Any]:
    """The owner's rule (2026-11-05): an engine-heavy run starts only at >= ``min_free_percent`` free memory
    (macOS memory_pressure); swap use is recorded. Refuses (SystemExit) below it."""
    try:
        txt = subprocess.run(["memory_pressure"], capture_output=True, text=True, check=True).stdout
        import re

        m = re.search(r"System-wide memory free percentage:\s*(\d+)%", txt)
        free = float(m.group(1)) if m else None
    except (OSError, subprocess.CalledProcessError):
        free = None
    state = {"free_percent": free, "swap": M1.swap_usage(), "load": [round(x, 2) for x in os.getloadavg()],
             "workers": parallel.workers(), "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    if free is not None and free < float(min_free_percent):
        raise SystemExit(f"memory gate: {free:.0f} % free < {min_free_percent:.0f} % (swap {state['swap']}); wait")
    return state


# ----------------------------------------------------------------------------------------------------------------
# --checks
# ----------------------------------------------------------------------------------------------------------------
FINGERPRINT_CODE = r"""
import hashlib, json, sys
sys.path[:0] = [sys.argv[1], sys.argv[1] + "/scripts/generators"]
import numpy as np
import generate_kinetic_core_b60_fit as B60
from generate_kinetic_core_b8_laplace import frozen_vector
from src.kinetic_core import config, sulfur
from src.kinetic_core.ph_state import BufferSpec
B23 = B60.B23
lo, hi = B60.full_bounds()
x = np.clip(frozen_vector(json.loads(open(sys.argv[2]).read())), lo, hi)
params, drift = B23.build_parameters(x), B23.unpack(x)[3]
pots = json.loads(sys.argv[3])
h = hashlib.sha256()
cfg = config.Config.shipped().replace({"ph.coupling": "inline"})
out = {}
with config.applied(cfg):
    for p in pots:
        buf = BufferSpec(**p["buffer"])
        r = sulfur.integrate_sulfur(params, p["t_c"] + 273.15, dict(p["initial"]), np.array(p["times"], float),
                                    ph=p["label"], buffer_spec=buf, ph_drift=drift, rtol=1e-8, atol=1e-14)
        a = np.ascontiguousarray(r.concentrations, dtype=float).tobytes() + np.ascontiguousarray(r.ph_series).tobytes()
        out[p["name"]] = hashlib.sha256(a).hexdigest()
print(json.dumps(out))
"""


def _buffer_dict(b: BufferSpec) -> Dict[str, Any]:
    return {f.name: getattr(b, f.name) for f in dataclasses.fields(b)}


def bit_identity(base: Optional[str] = None) -> Dict[str, Any]:
    """The shipped default is bit-identical: every control pot's run (state matrix and in-situ pH series) at B60
    pinned on inline, hashed with THIS tree's src and with the base commit's (git archive), compared."""
    _init_pots()
    base = base or subprocess.run(["git", "merge-base", "HEAD", "origin/main"], cwd=str(ROOT), capture_output=True,
                                  text=True, check=True).stdout.strip()
    pots = [{"name": p.name, "initial": list(p.initial), "t_c": p.t_c, "label": p.label_ph, "times": list(p.times),
             "buffer": _buffer_dict(p.buffer)} for p in POTS if not p.forced]
    env = dict(os.environ, MAILLARD_DATA_ROOT=str(data_paths.DATA_ROOT),
               MAILLARD_RESULTS_ROOT=str(data_paths.RESULTS_ROOT), MAILLARD_WORKERS="1")
    with tempfile.TemporaryDirectory() as tmp:
        arch = subprocess.run(["git", "archive", base, "src", "scripts/generators"], cwd=str(ROOT),
                              capture_output=True, check=True).stdout
        subprocess.run(["tar", "-x", "-C", tmp], input=arch, check=True)
        out = {}
        for label, where in (("base", tmp), ("this_tree", str(ROOT))):
            r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", FINGERPRINT_CODE, where,
                                str(M1.B60_REPORT), json.dumps(pots)], capture_output=True, text=True, env=env,
                               cwd=where)
            if r.returncode != 0:
                raise RuntimeError(f"fingerprint on {label} failed: {r.stderr[-2000:]}")
            out[label] = json.loads(r.stdout.strip().splitlines()[-1])
    same = {k: out["base"][k] == out["this_tree"][k] for k in out["base"]}
    return {"what": "every rival key off: each control pot's state matrix and in-situ pH series (B60 pinned, inline, "
                    "rtol 1e-8), sha256 with the base commit's src and with this tree's",
            "base_commit": base, "pots": same, "passed": all(same.values()), "hashes": out}


def v1_cost() -> Dict[str, Any]:
    """V1: every new key off reproduces B60's cost on inline (M1's V6 reference: R0 at B60's vector with the ceiling
    and the prior, inline, fit tolerance)."""
    m1 = json.loads(M1_CHECKS.read_text(encoding="utf-8"))
    ref = float(m1["checks"]["V6"]["cost"])
    numerics, fit_side = M1.numerics_for_fits(m1), M1.fit_side_from_checks(m1)
    t0 = time.time()
    r, failed = M1.evaluate(M1.EvalSpec(numerics, fit_side, False), M1.b60_point())
    cost = 0.5 * float(r @ r)
    rel = abs(cost - ref) / abs(ref)
    return {"what": V1_WHAT, "cost": cost, "m1_v6_reference": ref, "relative_difference": rel,
            "bit_identical": cost == ref, "tolerance": M1.V1_REL_TOL, "passed": rel <= M1.V1_REL_TOL and not failed,
            "failed_pot_runs": len(failed), "seconds": round(time.time() - t0, 1)}


V1_WHAT = "every new key off: R0 at B60's vector, inline, fit tolerance, ceiling + prior (M1's V6 reference)"


def _run_single(rival: str, values: Mapping[str, float], pot: Pot, breakpoints, form: str = "primary",
                extra_changes: Optional[Mapping[str, Any]] = None):
    """One pot run in this process (the checks)."""
    changes, drift = rival_changes(rival, form, values, pot.forced, breakpoints)
    changes.update(extra_changes or {})
    with config.applied(config.Config.shipped().replace(changes)):
        run = sulfur.integrate_sulfur(b60()[1], pot.t_c + CELSIUS, pot.initial_map(), np.array(pot.times, float),
                                      ph=pot.label_ph, buffer_spec=pot.buffer, ph_drift=drift, rtol=SCORE_RTOL,
                                      atol=ATOL)
    sid = float(run.metadata["strong_ion_difference_mmol_L"])
    cooled = [float(observations.ph_cooled(state_at(run, i), sid, drift, pot.buffer)) for i in range(len(pot.times))]
    return run, np.array(cooled), drift, sid


def corner_checks(breakpoints) -> Dict[str, Any]:
    """V2, V3: each rival at its B0 corner reproduces B0-pinned on every control to 1e-6 pH (B2 at f_ret 0 and the
    lag's fast end; B1 at phi 0; B3 at y 0 with the lumped yield at B60's value)."""
    _init_pots()
    d0 = b0_pinned_drift()
    pinned = {"acid_yield": d0.acid_yield, "arp_amine_pka": d0.arp_amine_pka}
    corners = {"V3_B1_phi0": ("B1", {**pinned, "phi_cys": 0.0, "phi_arp": 0.0}),
               "V2_B2_fret0_klag1e6": ("B2", {**pinned, "f_ret": 0.0, "k_lag": 1e6, "ea_lag": 103.7}),
               "V2_B2_fret0_klag10_band_end": ("B2", {**pinned, "f_ret": 0.0, "k_lag": 10.0, "ea_lag": 103.7})}
    out: Dict[str, Any] = {}
    refs = {p.name: _run_single("B0_pinned", {}, p, breakpoints)[1] for p in POTS}
    for name, (rival, vals) in corners.items():
        diffs = {p.name: float(np.max(np.abs(_run_single(rival, vals, p, breakpoints)[1] - refs[p.name])))
                 for p in POTS}
        out[name] = {"max_abs_diff_ph": max(diffs.values()), "per_pot": diffs,
                     "passed": max(diffs.values()) <= 1e-6 if "band_end" not in name else None}
    # B3 at y = 0 with the lumped yield at B60's value: the supply key on at zero shares, the yield kept
    diffs = {}
    for p in POTS:
        changes, _ = rival_changes("B0_pinned", "primary", {}, p.forced, breakpoints)
        changes["ph.acid_supply"] = pbr.AcidSupply(True, 0.0, 0.0, 0.0, tuple(breakpoints), 0.1, EA_ISO_KJ_MOL)
        with config.applied(config.Config.shipped().replace(changes)):
            run = sulfur.integrate_sulfur(b60()[1], p.t_c + CELSIUS, p.initial_map(), np.array(p.times, float),
                                          ph=p.label_ph, buffer_spec=p.buffer, ph_drift=d0, rtol=SCORE_RTOL, atol=ATOL)
        sid = float(run.metadata["strong_ion_difference_mmol_L"])
        got = np.array([observations.ph_cooled(state_at(run, i), sid, d0, p.buffer) for i in range(len(p.times))])
        diffs[p.name] = float(np.max(np.abs(got - refs[p.name])))
    out["V3_B3_y0_yield_B60"] = {"max_abs_diff_ph": max(diffs.values()), "per_pot": diffs,
                                  "passed": max(diffs.values()) <= 1e-6}
    out["passed"] = all(v["passed"] for k, v in out.items() if isinstance(v, dict) and v.get("passed") is not None)
    out["note"] = ("V2's prereg form is 'k_lag at its fast end'; the band's numerical end is 10 /min (lifetime "
                   "0.1 min), which leaves a latent carboxyl of ~(ARP loss rate / k_lag) at every time: the 1e-6 "
                   "test is run at k_lag 1e6 /min and the 10 /min distance printed beside it")
    return out


def sid_invariance(breakpoints) -> Dict[str, Any]:
    """V4: the conserved SID closes the charge balance at every node: the in-situ pH re-solved from each state (with
    the rival's pools) at the run's SID equals the run's own in-situ pH."""
    _init_pots()
    mid = {"B0": {"acid_yield": 0.5, "arp_amine_pka": 7.6},
           "B1": {"acid_yield": 0.5, "arp_amine_pka": 7.6, "phi_cys": 0.5, "phi_arp": 0.5},
           "B2": {"acid_yield": 0.5, "arp_amine_pka": 7.6, "f_ret": 0.5, "k_lag": 0.3, "ea_lag": 103.7},
           "B3": {"arp_amine_pka": 7.6, "y5": 0.3, "y6": 0.6, "y_hi": 0.8, "k_iso": 0.08},
           "B13": {"arp_amine_pka": 7.6, "y5": 0.3, "y6": 0.6, "y_hi": 0.8, "k_iso": 0.08, "phi_cys": 0.5,
                   "phi_arp": 0.5}}
    worst: Dict[str, float] = {}
    for rival, vals in mid.items():
        w = 0.0
        for p in (POTS_BY_NAME["mey_cys_pH5.5"], POTS_BY_NAME["hou_forced"], POTS_BY_NAME["zhou_cys_pH7"]):
            changes, drift = rival_changes(rival, "primary", vals, p.forced, breakpoints)
            with config.applied(config.Config.shipped().replace(changes)):
                run = sulfur.integrate_sulfur(b60()[1], p.t_c + CELSIUS, p.initial_map(), np.array(p.times, float),
                                              ph=p.label_ph, buffer_spec=p.buffer, ph_drift=drift, rtol=SCORE_RTOL,
                                              atol=ATOL)
                sid = float(run.metadata["strong_ion_difference_mmol_L"])
                for i in range(len(p.times)):
                    inv = ph_state.titratable_inventory(state_at(run, i), drift, p.buffer)
                    hot = ph_state.solve_ph(p.t_c + CELSIUS, sid, inv)
                    w = max(w, abs(float(hot) - float(run.ph_series[i])))
        worst[rival] = w
    return {"max_abs_insitu_ph_resolve_diff": worst, "passed": max(worst.values()) <= 1e-6,
            "what": "the in-situ pH re-solved at the run's single conserved SID from every grid state equals the "
                    "run's in-situ pH (a mid-band draw of each rival; Meynier Cys pH 5.5, Hou forced, Zhou Cys pH 7)"}


def v6_reprint(breakpoints) -> Dict[str, Any]:
    """V6: CF-084's engine numbers re-printed at B60, cooled, inline (B0-pinned), before the screen; sec. 3.2's
    expected signs checked against them."""
    _init_pots()
    data = measured_units()
    rows: Dict[str, Any] = {}
    for p in POTS:
        _run, cooled, _d, _s = _run_single("B0_pinned", {}, p, breakpoints)
        rows[p.name] = {"label_ph": p.label_ph, "times": list(p.times), "cooled_ph": cooled.tolist(),
                        "scale": "cooled", "T_C_reading": 25.0}
    mey = {p: rows[p]["cooled_ph"][-1] - rows[p]["label_ph"] for p in UNIT_POTS["C-Mey-Cys"]}
    mey_gly = {p: rows[p]["cooled_ph"][-1] - rows[p]["label_ph"] for p in UNIT_POTS["C-Mey-Gly"]}
    hou = np.array(rows["hou_unforced"]["cooled_ph"][1:])
    hou_f = np.array(rows["hou_forced"]["cooled_ph"][1:])
    res = hou - np.array(HOU["ph_cooled"])
    res_f = hou_f - np.array(HOU["ph_cooled"])
    zhou_a = [rows[f"zhou_arp_pH{x:g}"]["cooled_ph"][-1] for x in ZHOU_LABELS]
    zhou_c = [rows[f"zhou_cys_pH{x:g}"]["cooled_ph"][-1] for x in ZHOU_LABELS]
    pka_rms = {}
    p = POTS_BY_NAME["hou_unforced"]
    for pka in (7.07, 7.6, 8.18):
        _r, c, _d, _s = _run_single("B0", {"acid_yield": b0_pinned_drift().acid_yield, "arp_amine_pka": pka}, p,
                                    breakpoints)
        pka_rms[f"{pka:g}"] = float(np.sqrt(np.mean((c[1:] - np.array(HOU["ph_cooled"])) ** 2)))
    cao = cao_reprint()
    signs = {
        "meynier_cys_too_acid_beyond_0p2": {"expected": "drift < -0.2 at three labels (CF-084)",
                                            "got": {k: round(v, 3) for k, v in mey.items()},
                                            "n_below_minus_0p2": int(sum(v < -0.2 for v in mey.values()))},
        "hou_too_basic_every_time": {"expected": "pred - meas > 0 at all six times", "unforced": res.tolist(),
                                     "forced": res_f.tolist(), "unforced_holds": bool(np.all(res > 0)),
                                     "forced_holds": bool(np.all(res_f > 0))},
        "hou_level_fall_0_20": {"expected": "B0 0.65 (measured 1.95)", "unforced": 7.0 - float(hou[0]),
                                "forced": 7.0 - float(hou_f[0])},
        "zhou_spread_6_8": {"expected": "B0 3.90 (measured 1.00)", "got": zhou_a[-1] - zhou_a[0],
                            "measured": data["C-Zhou-ARP"].y[-1] - data["C-Zhou-ARP"].y[0]},
    }
    differ = []
    if signs["meynier_cys_too_acid_beyond_0p2"]["n_below_minus_0p2"] < 1:
        differ.append("Meynier: B0 is not too acid at B60")
    if not signs["hou_too_basic_every_time"]["unforced_holds"]:
        differ.append("Hou (unforced): B0 is not too basic at every time at B60")
    if not signs["hou_too_basic_every_time"]["forced_holds"]:
        differ.append("Hou (forced, the decisive form): B0 is not too basic at every time at B60")
    if not (signs["zhou_spread_6_8"]["got"] > signs["zhou_spread_6_8"]["measured"]):
        differ.append("Zhou: B0's endpoint spread is not larger than the measured one")
    return {"what": "CF-084's numbers re-printed at B60 (B0-pinned: yield %.3g, pKa %.4f), cooled, inline, rtol 1e-8"
                    % (b0_pinned_drift().acid_yield, b0_pinned_drift().arp_amine_pka),
            "meynier_cys_drift": mey, "meynier_gly_drift": mey_gly, "hou_unforced_cooled": hou.tolist(),
            "hou_forced_cooled": hou_f.tolist(), "hou_residuals_unforced": res.tolist(),
            "hou_residuals_forced": res_f.tolist(), "hou_rms_by_arp_pka_unforced": pka_rms,
            "zhou_arp_alone_final": zhou_a, "zhou_arp_cys_final": zhou_c,
            "zhou_measured": {u: list(data[u].y) for u in ("C-Zhou-ARP", "C-Zhou-Cys")}, "cao": cao,
            "sec_3_2_expectations": signs, "signs_that_differ": differ, "rows": rows,
            "b0_pinned_zhou_spread": zhou_a[-1] - zhou_a[0]}


def cao_reprint() -> Dict[str, Any]:
    """Cao 2017's two arms (ph_acid_submodel.yml), cooled, under R0 and C1 at B60 pinned (sec. 3.9's descriptive
    cross-check; printed, never scored)."""
    pots, rows = fit_rows.load(ZHOU_TABLE)
    out = {}
    x, params, drift = b60()
    for pot in sorted(k for k in pots if k.startswith("cao2017")):
        prs = [r for r in rows if r.pot == pot and r.species == "pH"]
        times = sorted({0.0, *(float(r.t_min) for r in prs)})
        job = M1.table_job(ZHOU_TABLE, pot, times)
        job.update(parameters=params, drift=drift, rtol=SCORE_RTOL, atol=ATOL)
        for structure in ("R0", "C1"):
            changes: Dict[str, Any] = {"ph.coupling": "inline"}
            if structure == "C1":
                changes["sulfur.cys_amadori_route"] = dataclasses.replace(
                    M1.shipped_value("sulfur.cys_amadori_route"), enabled=True)
            with config.applied(config.Config.shipped().replace(changes)):
                res = probe.run_pot(dict(job))
            pred = {float(t): float(observations.ph_cooled(st, res["sid"], drift, probe.buffer_of(
                {"buffer": job["buffer"]}, job["label_ph"]))) for t, st in zip(job["times"], res["states"])}
            out[f"{pot}|{structure}"] = [{"t_min": float(r.t_min), "measured": float(r.value), "predicted":
                                          pred[float(r.t_min)], "scale": "cooled"} for r in prs]
    return out


def v7_forcing(breakpoints) -> Dict[str, Any]:
    """V7: with forcing on, the rival's ARP reproduces Hou's interpolant to 1e-6 relative at every node; with forcing
    off the run is bit-identical to the unforced rival; the forced loss is never negative; B3's branch rule holds
    DPO's total loss to 1e-9 relative at every right-hand-side call."""
    _init_pots()
    pf, pu = POTS_BY_NAME["hou_forced"], POTS_BY_NAME["hou_unforced"]
    vals = {"B0": {"acid_yield": 0.3, "arp_amine_pka": 7.6},
            "B3": {"arp_amine_pka": 7.6, "y5": 0.3, "y6": 0.6, "y_hi": 0.9, "k_iso": 0.08}}
    out: Dict[str, Any] = {}
    for rival, v in vals.items():
        run_f, *_ = _run_single(rival, v, pf, breakpoints)
        diag = run_f.metadata["ph_balance_rivals"]["diagnostics"]
        run_u1, c1, *_ = _run_single(rival, v, pu, breakpoints)
        run_u2, c2, *_ = _run_single(rival, v, pu, breakpoints,
                                     extra_changes={"ph.arp_forcing": dataclasses.replace(hou_forcing(),
                                                                                          enabled=False)})
        same = bool(np.array_equal(run_u1.concentrations, run_u2.concentrations) and np.array_equal(c1, c2))
        out[rival] = {"forced_arp_max_rel_error": diag["forced_arp_max_rel_error"],
                      "forced_loss_never_negative": diag["forced_loss_never_negative"],
                      "forced_unresolvable_calls": diag["forced_unresolvable_calls"],
                      "forcing_off_bit_identical_to_unforced": same,
                      "forced_arp": run_f.series("ARP").tolist(), "measured_arp": [HOU["arp0_mmol_l"]] + list(
                          HOU["arp_mmol_l"])}
        if rival == "B3":
            out[rival]["dpo_total_rel_dev_max"] = diag["dpo_total_rel_dev_max"]
            run_m, *_ = _run_single(rival, v, POTS_BY_NAME["mey_cys_pH6.5"], breakpoints)
            out[rival]["dpo_total_rel_dev_max_meynier_pH6.5"] = (
                run_m.metadata["ph_balance_rivals"]["diagnostics"]["dpo_total_rel_dev_max"])
    ok = all((r["forced_arp_max_rel_error"] or 1.0) <= 1e-6 and r["forced_loss_never_negative"]
             and r["forcing_off_bit_identical_to_unforced"] and r["forced_unresolvable_calls"] == 0
             for r in out.values())
    ok = ok and out["B3"]["dpo_total_rel_dev_max"] <= 1e-9 and out["B3"]["dpo_total_rel_dev_max_meynier_pH6.5"] <= 1e-9
    out["passed"] = bool(ok)
    return out


V_SOLVER_TOL_PH = 1e-4
#: the draw at which LSODA stalled in the first screen attempt (2026-11-05; the attempt was stopped and discarded)
LSODA_STALL_DRAW = ("B0", {"acid_yield": 0.3952675089240074, "arp_amine_pka": 7.62728870975785}, "mey_cys_pH5.5")


def v_solver(breakpoints) -> Dict[str, Any]:
    """V-solver: the declared BDF fallback agrees with LSODA (cooled pH) on every control pot at B0-pinned and at a
    mid-band B13 draw, and with Radau at the draw where LSODA stalled."""
    _init_pots()
    out: Dict[str, Any] = {"per_pot": {}}
    worst = 0.0
    mid13 = {"arp_amine_pka": 7.6, "y5": 0.3, "y6": 0.6, "y_hi": 0.8, "k_iso": 0.08, "phi_cys": 0.5, "phi_arp": 0.5}
    for rival, vals in (("B0_pinned", {}), ("B13", mid13)):
        for p in POTS:
            changes, drift = rival_changes(rival, "primary", vals, p.forced, breakpoints)
            got = {}
            for m in ("LSODA", FALLBACK_METHOD):
                with config.applied(config.Config.shipped().replace(changes)):
                    run = integrate_with_fallback(p, drift, method=m)
                sid = float(run.metadata["strong_ion_difference_mmol_L"])
                got[m] = np.array([observations.ph_cooled(state_at(run, i), sid, drift, p.buffer)
                                   for i in range(len(p.times))])
            d = float(np.max(np.abs(got["LSODA"] - got[FALLBACK_METHOD])))
            out["per_pot"][f"{rival}|{p.name}"] = d
            worst = max(worst, d)
    rival, vals, pot = LSODA_STALL_DRAW
    p = POTS_BY_NAME[pot]
    changes, drift = rival_changes(rival, "primary", vals, False, breakpoints)
    finals = {}
    for m in (FALLBACK_METHOD, "Radau"):
        with config.applied(config.Config.shipped().replace(changes)):
            run = M1.with_alarm(600, lambda: sulfur.integrate_sulfur(
                b60()[1], p.t_c + CELSIUS, p.initial_map(), np.array(p.times, float), ph=p.label_ph,
                buffer_spec=p.buffer, ph_drift=drift, rtol=SCORE_RTOL, atol=ATOL, method=m))
        sid = float(run.metadata["strong_ion_difference_mmol_L"])
        finals[m] = float(observations.ph_cooled(state_at(run, len(p.times) - 1), sid, drift, p.buffer))
    stall = abs(finals[FALLBACK_METHOD] - finals["Radau"])
    out.update({"max_abs_lsoda_vs_bdf_ph": worst, "stall_draw": {"rival": rival, "coordinates": vals, "pot": pot,
                "cooled_ph": finals, "bdf_vs_radau_ph": stall},
                "tolerance_ph": V_SOLVER_TOL_PH, "passed": bool(worst <= V_SOLVER_TOL_PH and stall <= V_SOLVER_TOL_PH),
                "what": "LSODA vs the declared BDF fallback on every control pot (B0-pinned, a mid-band B13 draw); "
                        "BDF vs Radau at LSODA's stall draw"})
    return out


def v_density() -> Dict[str, Any]:
    """V-density: the estimator reproduces M1's known case (R0's RES_Cys, alternative floor, ~ -746): M1.1's tests on
    its engine-scan fixture, run here."""
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests/unit/test_predictive_densities.py", "-k", "exact_density_of_r0_cys"],
                       cwd=str(ROOT), capture_output=True, text=True)
    tail = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-300:]
    return {"what": "tests/unit/test_predictive_densities.py -k exact_density_of_r0_cys (M1.1's -746 by quadrature from an engine scan; "
                    "PSIS with a width-scaled proposal agrees)", "returncode": r.returncode, "summary": tail,
            "passed": r.returncode == 0}


def v5_scales(breakpoints_record) -> Dict[str, Any]:
    """V5: every pH comparison declares both sides' scale and temperature; a mixed-scale comparison is a build
    error (``predictive_densities.require_ph_scale``)."""
    comparisons = []
    for unit, pots in UNIT_POTS.items():
        for p in pots:
            pot = POTS_BY_NAME[p]
            comparisons.append({"unit": unit, "pot": p, "predicted": {"scale": "cooled", "T_C": 25.0,
                                                                       "reaction_T_C": pot.t_c},
                                "measured": {"scale": "cooled", "T_C": 25.0}})
            pdens.require_ph_scale("cooled", "cooled")
    try:
        pdens.ph_difference(pdens.PhReading(4.0, "in_situ", 120.0), pdens.PhReading(4.0, "cooled", 25.0))
        guard = False
    except pdens.MeasurementModelRefusal:
        guard = True
    return {"comparisons": comparisons, "davidek_breakpoints": breakpoints_record, "mixed_scale_refused": guard,
            "passed": guard}


def run_checks(base: Optional[str]) -> Dict[str, Any]:
    _init_pots()
    t0 = time.time()
    bp = davidek_breakpoints()
    br = bp["in_situ"]
    checks: Dict[str, Any] = {}
    for name, fn in (("BIT", lambda: bit_identity(base)), ("V1", v1_cost), ("V2_V3", lambda: corner_checks(br)),
                     ("V4", lambda: sid_invariance(br)), ("V5", lambda: v5_scales(bp)),
                     ("V6", lambda: v6_reprint(br)), ("V7", lambda: v7_forcing(br)), ("V_solver", lambda: v_solver(br)), ("V_density", v_density)):
        s = time.time()
        print(f"--checks {name} ...", flush=True)
        checks[name] = fn()
        checks[name]["seconds"] = round(time.time() - s, 1)
        print(f"--checks {name}: passed {checks[name].get('passed')} ({checks[name]['seconds']} s)", flush=True)
    rec = {"artifact": f"{STEM}_checks", "prereg": data_paths.rel(PREREG), "checks": checks,
           "davidek_breakpoints_in_situ": br, "machine": M1.machine_state(), "seconds": round(time.time() - t0, 1),
           "passed": all(v.get("passed") is not False for v in checks.values()),
           "failed_pot_runs": 0, "provenance": _provenance("checks")}
    write_record("checks", rec, render_checks(rec))
    return rec


# ----------------------------------------------------------------------------------------------------------------
# --screen
# ----------------------------------------------------------------------------------------------------------------
SCREEN_CHUNK = 32                 # draws per checkpointed chunk
SCREEN_PARTIAL = V / f"{STEM}_screen_partial.json"
MEMORY_WAIT_SECONDS = 600         # the owner's rule: check at most every 10 min
SWAP_CLIMB_BYTES = 256 * 2 ** 20


def _swap_used_bytes() -> Optional[float]:
    import re

    s = M1.swap_usage() or ""
    m = re.search(r"used = ([\d.]+)M", s)
    return float(m.group(1)) * 2 ** 20 if m else None


def memory_wait(min_free: float, events: List[Dict[str, Any]], last_swap: List[Optional[float]],
                workers_max: int) -> int:
    """Before every engine-heavy chunk (the owner's rule, 2026-11-05, held strictly): a chunk starts only at >=
    ``min_free`` % free memory (memory_pressure) with swap not climbing; otherwise wait, re-checking every 10 min.
    Every decision is recorded. Returns the worker count (``workers_max`` <= 2)."""
    while True:
        try:
            st = memory_gate(0.0)
        except SystemExit:
            st = {}
        free, swap = st.get("free_percent"), _swap_used_bytes()
        climbing = (swap is not None and last_swap[0] is not None and swap > last_swap[0] + SWAP_CLIMB_BYTES)
        last_swap[0] = swap if swap is not None else last_swap[0]
        if not climbing and free is not None and free >= min_free:
            events.append({**st, "decision": f"run, {workers_max} workers"})
            return workers_max
        why = "free memory unknown" if free is None else (f"{free:.0f} % free < {min_free:.0f} %" if free < min_free
                                                           else "swap climbing")
        events.append({**st, "decision": f"wait {MEMORY_WAIT_SECONDS} s ({why})"})
        print(f"memory gate: {why} ({st.get('swap')}); waiting {MEMORY_WAIT_SECONDS} s", flush=True)
        time.sleep(MEMORY_WAIT_SECONDS)


def concat(a: Optional[Ensemble], b: Ensemble) -> Ensemble:
    if a is None:
        return b
    return Ensemble(a.rival, a.form, a.draws + b.draws, {u: np.vstack([a.f[u], b.f[u]]) for u in a.f},
                    {u: {k: np.vstack([a.delta[u][k], b.delta[u][k]]) for k in a.delta[u]} for u in a.delta},
                    {u: np.vstack([a.extra[u], b.extra[u]]) for u in a.extra}, a.failed + b.failed, a.fallbacks + b.fallbacks)


def screen_draws(rival: str, n: int) -> List[Dict[str, float]]:
    return draw_coordinates(rival, "primary", n, SCREEN_SEED + RIVAL_NAMES.index(rival))


def screen_chunks(rival: str, n_total: int, breakpoints, partial: Dict[str, Any], failures: List[Dict[str, Any]],
                  min_free: float, workers_max: int, last_swap: List[Optional[float]]) -> Ensemble:
    """The rival's first ``n_total`` Sobol draws, in checkpointed chunks (a stop costs one chunk)."""
    draws = screen_draws(rival, n_total)
    rec = partial.setdefault(rival, {"chunks": [], "memory": [], "seconds": 0.0})
    done = sum(len(c["draws"]) for c in rec["chunks"])
    while done < n_total:
        chunk = draws[done:done + SCREEN_CHUNK]
        n_w = memory_wait(min_free, rec["memory"], last_swap, workers_max)
        t0 = time.time()
        res = run_jobs(build_tasks(rival, "primary", chunk, breakpoints, pots_for("primary", rival)), failures,
                       n_workers=n_w)
        ens = assemble(rival, "primary", chunk, res)
        rec["chunks"].append(ensemble_to_json(ens))
        rec["seconds"] += time.time() - t0
        rec["failures"] = failures[-50:]
        done += len(chunk)
        SCREEN_PARTIAL.write_text(json.dumps(_jsonable(partial)) + "\n", encoding="utf-8")
        print(f"screen {rival}: {done}/{n_total} draws, chunk {time.time() - t0:.0f} s at {n_w} worker(s), "
              f"failed so far {len(failures)}", flush=True)
    out: Optional[Ensemble] = None
    for c in rec["chunks"]:
        out = concat(out, ensemble_from_json(c))
    assert out is not None
    return out.subset(range(n_total)) if out.size > n_total else out


def best_draw_extras(rival: str, ens: Ensemble, s: int, breakpoints, failures) -> Dict[str, Any]:
    """C-Hou's printed inversion at the best draw: the acid equivalents the rival's titration needs to reach Hou's
    measured cooled pH, beside the rival's own produced acid (CF-084's remedy; printed, never scored)."""
    pot = POTS_BY_NAME["hou_forced"]
    res = run_jobs(build_tasks(rival, "primary", [ens.draws[s]], breakpoints, [pot], acid_eq=True), failures)[0]
    if "failed" in res:
        return {"failed": res["failed"]}
    used = HOU["arp0_mmol_l"] - np.array(HOU["arp_mmol_l"])
    eq = [None if v is None else float(v) for v in res["acid_eq"]]
    return {"acid_eq_to_measured_mmol_l": eq, "own_acid_mmol_l": res["own_acid"],
            "measured_acid_per_arp_consumed": [None if e is None else float((o + e) / u)
                                               for o, e, u in zip(res["own_acid"], eq, used)],
            "own_acid_per_arp_consumed": [float(o / u) for o, u in zip(res["own_acid"], used)]}


def coverage(ens: Ensemble, passed: np.ndarray) -> Dict[str, Any]:
    coords = [c.name for c in coordinates(ens.rival)]
    pts = np.array([[d[c] for c in coords] for d, ok in zip(ens.draws, passed) if ok])
    return {"fraction_passing": float(np.mean(passed)), "n_passing": int(passed.sum()),
            "bounding_box": ({c: [float(pts[:, i].min()), float(pts[:, i].max())] for i, c in enumerate(coords)}
                             if pts.size else None)}


def run_screen(min_free: float) -> Dict[str, Any]:
    checks = read_record("checks")
    if not checks.get("passed"):
        raise SystemExit("--screen refuses: --checks has not passed (V1 - V7, BIT, V-density first)")
    if read_record("screen").get("complete"):
        raise SystemExit("--screen has run (sec. 8: a screen is run once); its record stands")
    _init_pots()
    br = tuple(checks["davidek_breakpoints_in_situ"])
    b0_spread = float(checks["checks"]["V6"]["b0_pinned_zhou_spread"])
    workers_max = min(MAX_WORKERS, parallel.workers(1))
    gate = memory_gate(0.0)
    partial = json.loads(SCREEN_PARTIAL.read_text(encoding="utf-8")) if SCREEN_PARTIAL.exists() else {}
    failures: List[Dict[str, Any]] = list(partial.get("_failures", []))
    last_swap: List[Optional[float]] = [None]
    t0 = time.time()
    ens_json: Dict[str, Any] = {}
    per: Dict[str, Any] = {}
    data = measured_units()
    # B0-pinned, the single point (N2)
    pin_results = run_jobs(build_tasks("B0_pinned", "primary", [{}], br, pots_for("primary")), failures)
    pinned = assemble("B0_pinned", "primary", [{}], pin_results)
    ens_json["B0_pinned"] = ensemble_to_json(pinned)
    for rival in RIVAL_NAMES:
        ens = screen_chunks(rival, SCREEN_DRAWS, br, partial, failures, min_free, workers_max, last_swap)
        partial["_failures"] = failures
        sig = signatures(ens, b0_spread)
        doubled = False
        if int(sig["all"].sum()) < COVERAGE_MIN_PASSING:
            ens = screen_chunks(rival, 2 * SCREEN_DRAWS, br, partial, failures, min_free, workers_max, last_swap)
            sig = signatures(ens, b0_spread)
            doubled = True
        secs = partial[rival]["seconds"]
        memory_events = partial[rival]["memory"]
        sig_u = signatures(ens, b0_spread, hou_unit="C-Hou-unforced")
        best = best_draw(ens, data)
        inversion = best_draw_extras(rival, ens, best, br, failures)
        per[rival] = {
            "draws": ens.size, "doubled_for_coverage": doubled, "seconds": round(secs, 1),
            "failed_pot_runs": ens.failed, "solver_fallbacks": ens.fallbacks, "memory": memory_events, "hou_inversion_at_best_draw": inversion,
            "pass_fraction": {k: float(np.mean(sig[k])) for k in ("sigma_hou_level", "sigma_hou_shape", "sigma_mey",
                                                                   "sigma_zhou", "all")},
            "pass_fraction_unforced_hou": {k: float(np.mean(sig_u[k])) for k in ("sigma_hou_level", "sigma_hou_shape",
                                                                                  "all")},
            "would_be_killed": bool(not sig["all"].any()),
            "kill_note": ("a would-be kill is confirmed at one fit (start 0, sec. 3.8), not run here"
                          if not sig["all"].any() else None),
            "coverage": coverage(ens, sig["all"]),
            "best_draw": {"index": best, "coordinates": ens.draws[best],
                          "hou_forced": ens.f["C-Hou"][best].tolist(), "hou_unforced": ens.f["C-Hou-unforced"][best]
                          .tolist(), "mey_cys_drift": ens.f["C-Mey-Cys"][best].tolist(),
                          "mey_gly_drift": ens.f["C-Mey-Gly"][best].tolist() if "C-Mey-Gly" in ens.f else None,
                          "zhou_arp": ens.f["C-Zhou-ARP"][best].tolist(), "zhou_cys": ens.f["C-Zhou-Cys"][best].tolist(),
                          "hou_arp_unforced": ens.extra["hou_arp_unforced"][best].tolist(),
                          "hou_ret_aa_forced": ens.extra["hou_ret_forced"][best].tolist()},
            "signature_summaries": {k: _quantiles(sig[k]) for k in ("hou_level_ph20", "hou_fall_60_120",
                                                                    "zhou_spread_6_8", "mey_min_p_within")},
            "mey_gly_drift_quantiles": (_quantiles(np.min(ens.f["C-Mey-Gly"], axis=1)) if "C-Mey-Gly" in ens.f
                                        else "not run: no pair of this rival reads the C-Mey (Gly) guard"),
            "mey_cys_aa_plus_laa_60_quantiles": _quantiles(np.max(ens.extra["mey_cys_aa_plus_laa_60"], axis=1)),
            "delta_ph_primary_quantiles": {u: _quantiles(np.mean(ens.delta[u]["primary"], axis=1)) for u in ens.delta},
        }
        ens_json[rival] = ensemble_to_json(ens)
        print(f"screen {rival}: {ens.size} draws, {per[rival]['pass_fraction']}, failed {ens.failed}, "
              f"{secs:.0f} s", flush=True)
    refused = bool(failures)
    if refused:
        print(f"screen: {len(failures)} failed pot runs -> the record is REFUSED (nothing downstream reads it)",
              flush=True)
    rec = {"artifact": f"{STEM}_screen", "prereg": data_paths.rel(PREREG), "complete": True,
           "refused_for_failed_pot_runs": refused, "failed_pot_runs": len(failures), "failures_sample": failures[:20],
           "memory_at_start": gate, "per_rival": per, "b0_pinned": ensemble_to_json(pinned),
           "b0_pinned_zhou_spread": b0_spread, "davidek_breakpoints_in_situ": br,
           "declared_readings": __doc__.split("DECLARED READINGS")[1].split('"""')[0].strip(),
           "demotions_implied": [], "seconds": round(time.time() - t0, 1), "machine": M1.machine_state(),
           "provenance": _provenance("screen")}
    import gzip

    with gzip.open(DRAWS_FILE, "wt", encoding="utf-8") as fh:
        fh.write(json.dumps(_jsonable(ens_json)) + "\n")
    rec["draws_file"] = data_paths.rel(DRAWS_FILE)
    write_record("screen", rec, render_screen(rec))
    SCREEN_PARTIAL.unlink(missing_ok=True)
    return rec


def _quantiles(x: np.ndarray) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return {}
    q = np.quantile(x, [0.0, 0.05, 0.5, 0.95, 1.0])
    return dict(zip(("min", "q05", "median", "q95", "max"), map(float, q)))


#: the screen's draw ensembles, gzipped (coordinator, 2026-11-06: frozen evidence of a non-reproducible 6.4 h run)
DRAWS_FILE = V / f"{STEM}_screen_draws.json.gz"


def load_screen_ensembles(path: Optional[Path] = None) -> Dict[str, Ensemble]:
    import gzip

    p = Path(path) if path is not None else DRAWS_FILE
    if not p.exists():
        raise SystemExit("no screen draws on disk (run --screen first)")
    raw = json.loads(gzip.open(p, "rt", encoding="utf-8").read())
    return {k: ensemble_from_json(v) for k, v in raw.items()}


# ----------------------------------------------------------------------------------------------------------------
# --power (sec. 7)
# ----------------------------------------------------------------------------------------------------------------
def synthetic(ens: Ensemble, s: int, rng: np.random.Generator, units: Iterable[str]) -> Dict[str, Dict[str, Any]]:
    """One synthetic data set from the truth draw s: the full Sigma_u (measurement + delta-pH at the primary omega).
    pH units: y = f + e; C-Mey: the drift observed with e ~ N(0, 2 r^2 + delta^2), read within the 0.2 edge or not.
    The unforced C-Hou is scored on the forced form's synthetic data (one data set)."""
    data = measured_units()
    out = {}
    for u in units:
        f, dl = ens.f[u][s], ens.delta[u]["primary"][s]
        if data[u].censored:
            d_obs = f + rng.normal(size=f.size) * np.sqrt(2.0 * METER_R ** 2 + dl ** 2)
            out[u] = {"within": np.abs(d_obs) <= MEYNIER_BOUND}
        else:
            out[u] = {"y": f + rng.normal(size=f.size) * np.sqrt(np.asarray(data[u].sd) ** 2 + dl ** 2)}
    if "C-Hou" in out:
        out["C-Hou-unforced"] = out["C-Hou"]
    return out


def _with(data: Mapping[str, UnitData], syn: Mapping[str, Mapping[str, Any]]) -> Dict[str, UnitData]:
    out = dict(data)
    for u, v in syn.items():
        if "y" in v:
            out[u] = dataclasses.replace(data[u], y=tuple(map(float, v["y"])))
        else:
            out[u] = dataclasses.replace(data[u], within=tuple(map(bool, v["within"])))
    return out


def truth_bands(a: str, ens: Mapping[str, Ensemble], truth: Ensemble, s: int, rng: np.random.Generator,
                b: str) -> Dict[str, pdens.TieBand]:
    """Tie bands once per pair and truth, at the truth's expected data (declared reading)."""
    data = measured_units()
    mean = {u: {"y": truth.f[u][s]} for u in truth.f if not data[u].censored}
    mean.update({u: {"within": np.abs(truth.f[u][s]) <= MEYNIER_BOUND} for u in truth.f if data[u].censored})
    dd = _with(data, mean)
    bands = {}
    forms = list(LIKELIHOOD_FORMS) + ([RuleForm("vs_b0_pinned", b_pinned=True)] if b == "B0" else [])
    for form in forms:
        eb = ens["B0_pinned"] if form.b_pinned else ens[b]
        for u in DECISIVE + guards_of(a):
            unit = form.hou_unit if u == "C-Hou" else u
            bands[f"{form.name}|{u}"] = tie_band(dd[unit], _arrays(ens[a], unit, form.omega),
                                                 _arrays(eb, unit, form.omega), rng, form.r_meter, n=TIE_R)
    return bands


def wilson(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def power_pair(a: str, b: str, ens: Mapping[str, Ensemble], rng: np.random.Generator, r: int = POWER_R
               ) -> Dict[str, Any]:
    data = measured_units()
    out: Dict[str, Any] = {"a": a, "b": b}
    units = DECISIVE + guards_of(a)
    for truth_name, key in ((a, "power"), (b, "false_credit")):
        truth = ens[truth_name]
        s = best_draw(truth, data)
        bands = truth_bands(a, ens, truth, s, rng, b)
        credited = credited_unresolved_tie = both = 0
        breaking: Dict[str, int] = {}
        for _ in range(r):
            syn = synthetic(truth, s, rng, units)
            dd = _with(data, syn)
            res = rule_over_forms(a, b, ens, dd, rng, bands=bands, require_resolved=False)
            credited += int(res["favoured"])
            both += int(res["both_inadequate_primary"])
            for f_ in res["breaking_forms"]:
                breaking[f_] = breaking.get(f_, 0) + 1
            if res["favoured"]:
                strict = rule_over_forms(a, b, ens, dd, rng, bands=bands, require_resolved=True)
                credited_unresolved_tie += int(strict["favoured"])
        out[key] = {"truth": truth_name, "truth_draw": s, "truth_coordinates": truth.draws[s], "n": r,
                    "credited": credited, "rate": credited / r, "wilson_95": wilson(credited, r),
                    "rate_unresolved_density_as_tie": credited_unresolved_tie / r,
                    "wilson_95_unresolved_as_tie": wilson(credited_unresolved_tie, r),
                    "p_both_inadequate": both / r, "forms_breaking_credit": breaking}
    out["decidable"] = bool(out["power"]["rate"] > out["false_credit"]["rate"])
    out["low_power"] = bool(out["power"]["rate"] < 0.5)
    out["verdict"] = ("undecidable on disk (power does not exceed false credit): demoted"
                      if not out["decidable"] else ("decidable, low power" if out["low_power"] else "decidable"))
    return out


def run_power() -> Dict[str, Any]:
    screen = read_record("screen")
    if not screen.get("complete"):
        raise SystemExit("--power refuses: the screen has not run")
    if screen.get("refused_for_failed_pot_runs"):
        raise SystemExit(f"--power refuses: the screen had {screen.get('failed_pot_runs')} failed pot runs")
    if read_record("power").get("complete"):
        raise SystemExit("--power has run; its record stands")
    ens = load_screen_ensembles()
    rng = np.random.default_rng(POWER_SEED)
    t0 = time.time()
    pairs = {}
    for a, b in PAIRS:
        s = time.time()
        pairs[f"{a}_vs_{b}"] = power_pair(a, b, ens, rng)
        pairs[f"{a}_vs_{b}"]["seconds"] = round(time.time() - s, 1)
        p = pairs[f"{a}_vs_{b}"]
        print(f"power {a} vs {b}: power {p['power']['rate']:.3f} false credit {p['false_credit']['rate']:.3f} "
              f"-> {p['verdict']}", flush=True)
    rec = {"artifact": f"{STEM}_power", "prereg": data_paths.rel(PREREG), "complete": True, "R": POWER_R,
           "tie_replicates": TIE_R, "pairs": pairs,
           "demoted": [k for k, v in pairs.items() if not v["decidable"]],
           "low_power": [k for k, v in pairs.items() if v["decidable"] and v["low_power"]],
           "forms_in_power_check": [f.name for f in LIKELIHOOD_FORMS] + ["vs_b0_pinned (b = B0)"],
           "forms_not_in_power_check": ["yield_loguniform", "delta_t_0", "added_sink", "sohn_phi_cys", "supply_c1"],
           "note": ("sec. 7: truth = the screen's best draw; both rivals scored on the screen's draw ensembles (the "
                    "truth draw included: no refit, overstates power); the full rule over the forms the draws carry. "
                    "'rate' reads an unaccepted PSIS estimate as estimated; 'rate_unresolved_density_as_tie' reads it "
                    "as a tie (the score's rule, which there uses unit-tuned proposals)."),
           "seconds": round(time.time() - t0, 1), "machine": M1.machine_state(), "failed_pot_runs": 0,
           "provenance": _provenance("power")}
    write_record("power", rec, render_power(rec))
    return rec


# ----------------------------------------------------------------------------------------------------------------
# --score (built and tested; RUN ONCE after the freeze; never run in this build)
# ----------------------------------------------------------------------------------------------------------------
SCORE_S0, SCORE_S_MAX = 256, 4096


def score_refusals(freeze_commit: str, freeze_check: Callable[[str], str] = None) -> List[str]:
    """Every reason --score refuses: the prereg unchanged since a freeze commit, the checks passed, the screen and the
    power check complete without failed pot runs, and the score not yet run (run once)."""
    problems = []
    check = freeze_check or (lambda c: M1.check_freeze(c, prereg=PREREG))
    try:
        check(freeze_commit)
    except SystemExit as exc:
        problems.append(str(exc))
    if not read_record("checks").get("passed"):
        problems.append("--checks has not passed")
    s = read_record("screen")
    if not s.get("complete"):
        problems.append("the screen has not run")
    elif s.get("refused_for_failed_pot_runs"):
        problems.append(f"the screen had {s.get('failed_pot_runs')} failed pot runs")
    if not read_record("power").get("complete"):
        problems.append("the power check has not run")
    if read_record("score"):
        problems.append("the score has run (run once): its record stands")
    return problems


def score_ensemble(rival: str, form: str, breakpoints, s: int, failures: List[Dict[str, Any]],
                   seed: int) -> Ensemble:
    draws = draw_coordinates(rival, form, s, seed)
    results = run_jobs(build_tasks(rival, form, draws, breakpoints, pots_for(form, rival)), failures)
    return assemble(rival, form, draws, results)


def run_score(freeze_commit: str, min_free: float) -> Dict[str, Any]:
    """secs. 3.7 and 5: every rival's lppd on every unit under every form, the rule per pair with required agreement
    over every form (likelihood-only forms reweight the same draws; engine-changing forms get their own draws), the
    adequacy gate, the guards' Holm step; S doubles from 256 to 4096 until every unit's density is accepted, and a
    unit still unresolved is a tie.

    KNOWN GAPS AGAINST SEC. 5.1 (to close before the freeze, or to amend in the prereg): the proposal is the
    rival's prior bands (Sobol), not the unit-tuned Student-t at the mode of likelihood x posterior; the
    Gauss-Hermite fallback on the <= 2 dominant coordinates is not built (an unresolved unit reads as a tie); the
    Laplace / profile parts of sec. 5.1 apply to refitted coordinates only (sec. 3.8), which --fit produces."""
    problems = score_refusals(freeze_commit)
    if problems:
        raise SystemExit("--score refuses:\n  " + "\n  ".join(problems))
    commit = M1.check_freeze(freeze_commit, prereg=PREREG)
    gate = memory_gate(min_free)
    checks = read_record("checks")
    br = tuple(checks["davidek_breakpoints_in_situ"])
    failures: List[Dict[str, Any]] = []
    rng = np.random.default_rng(POWER_SEED + 1)
    data = measured_units()
    ens: Dict[str, Dict[str, Ensemble]] = {f: {} for f in FORMS}
    pin = run_jobs(build_tasks("B0_pinned", "primary", [{}], br, pots_for("primary")), failures)
    pinned = assemble("B0_pinned", "primary", [{}], pin)
    for form in FORMS:
        for rival in FORM_RIVALS[form]:
            s = SCORE_S0
            while True:
                e = score_ensemble(rival, form, br, s, failures, POWER_SEED + 100 * FORMS.index(form)
                                   + RIVAL_NAMES.index(rival))
                if failures:
                    break
                unresolved = [u for u in e.f if not pdens.importance_lppd(
                    unit_loglik(data[u], e.f[u], e.delta[u]["primary"]), np.zeros(e.size)).accepted]
                if not unresolved or s >= SCORE_S_MAX:
                    break
                s *= 2
            ens[form][rival] = e
    if failures:
        rec = {"artifact": f"{STEM}_score", "refused_for_failed_pot_runs": True, "failed_pot_runs": len(failures),
               "failures_sample": failures[:20], "freeze_commit": commit, "provenance": _provenance("score")}
        write_record("score", rec, "# M1b part (i) score: REFUSED (failed pot runs)\n")
        return rec
    verdicts = {}
    for a, b in PAIRS:
        primary = {"B0_pinned": pinned, **ens["primary"]}
        res = rule_over_forms(a, b, primary, data, rng)
        engine_forms = {}
        for form in FORMS[1:]:
            if a in ens[form] or b in ens[form]:
                # a form's own draws for the units it runs (C1: the cysteine pots only); each unit's lppd is a
                # marginal over its own draws, so the units the form leaves unchanged keep the primary draws
                eset = dict(primary)
                for k, e in ens[form].items():
                    eset[k] = Ensemble(k, form, e.draws, {**primary[k].f, **e.f}, {**primary[k].delta, **e.delta},
                                       {**primary[k].extra, **e.extra}, e.failed)
                engine_forms[form] = rule_over_forms(a, b, eset, data, rng, forms=(RuleForm("primary"),))
        favoured = res["favoured"] and all(v["favoured"] for v in engine_forms.values())
        verdicts[f"{a}_vs_{b}"] = {"favoured": favoured, "likelihood_forms": res, "engine_forms": engine_forms}
    rec = {"artifact": f"{STEM}_score", "prereg": data_paths.rel(PREREG), "freeze_commit": commit,
           "verdicts": verdicts, "memory_at_start": gate, "failed_pot_runs": 0, "machine": M1.machine_state(),
           "provenance": _provenance("score")}
    write_record("score", rec, f"# M1b part (i): the decisive comparison\n\n{json.dumps(_jsonable({k: v['favoured'] for k, v in verdicts.items()}))}\n")
    return rec


# ----------------------------------------------------------------------------------------------------------------
# --fit: the stacking guard's refit of a favoured rival (sec. 3.8; built and tested; not run in this build)
# ----------------------------------------------------------------------------------------------------------------
FIT_COORDS = {"B1": ("phi_cys", "phi_arp"), "B2": ("f_ret", "k_lag"), "B3": ("y5", "y6", "y_hi", "k_iso"),
              "B13": ("y5", "y6", "y_hi", "k_iso", "phi_cys", "phi_arp")}


def fit_refusals(rival: str, freeze_commit: str, freeze_check: Callable[[str], str] = None) -> List[str]:
    problems = []
    check = freeze_check or (lambda c: M1.check_freeze(c, prereg=PREREG))
    try:
        check(freeze_commit)
    except SystemExit as exc:
        problems.append(str(exc))
    if rival not in FIT_COORDS:
        problems.append(f"{rival}: no stacking refit (B0 is the incumbent)")
    score = read_record("score")
    if not score:
        problems.append("the score has not run: the stacking guard refits favoured rivals only")
    elif not any(v.get("favoured") and k.startswith(f"{rival}_vs_") for k, v in (score.get("verdicts") or {}).items()):
        problems.append(f"{rival} is not favoured in the score: no refit (sec. 3.8)")
    return problems


class RivalObjective(M1.Objective):
    """M1's R0 objective (B60's 52 rows + the k_mft_decay prior, M1's settings) with a rival's switches on and its
    coordinates appended to R0's 17 (sec. 3.8). The B3 order is kept by clamping y6 >= y5 and y_hi >= y6."""

    def __init__(self, rival: str, numerics, fit_side, breakpoints):
        super().__init__(M1.RIVALS["R0"], numerics, fit_side)
        self.rival_name, self.breakpoints = rival, tuple(breakpoints)
        self.extra = FIT_COORDS[rival]
        bands = {c.name: (c.lo, c.hi) for c in coordinates(rival)}
        self.keys = self.keys + self.extra
        self.kinds = self.kinds + tuple("acid_yield" for _ in self.extra)
        self.lower = np.concatenate([self.lower, [bands[k][0] for k in self.extra]])
        self.upper = np.concatenate([self.upper, [bands[k][1] for k in self.extra]])

    def split(self, theta):
        n = len(self.keys) - len(self.extra)
        values = dict(zip(self.extra, map(float, theta[n:])))
        if "y5" in values:
            values["y6"] = max(values["y6"], values["y5"])
            values["y_hi"] = max(values["y_hi"], values["y6"])
        for k in ("k_iso", "k_lag"):
            if k in values:
                values[k] = 10.0 ** values[k]
        return theta[:n], values

    def rival_config(self, values) -> Dict[str, Any]:
        point_values = {"arp_amine_pka": 7.6, "acid_yield": 0.0, "ea_lag": 103.7, **values}
        changes, _ = rival_changes(self.rival_name, "primary", point_values, False, self.breakpoints)
        changes.pop("ph.coupling", None)
        return changes

    def residuals(self, theta, rows=None):
        base, values = self.split(np.asarray(theta, float))
        point = self.point_r0(base)
        cfg = M1.config_of(point, self.numerics, self.fit_side).replace(self.rival_config(values))
        failures: List[Dict[str, Any]] = []
        with config.applied(cfg), M1.b60_runner(self.numerics, failures):
            parts = [np.asarray(B60.residual_vector(point.x(), self.numerics.quick), dtype=float)]
            if self.fit_side.prior.enabled:
                parts.append(np.array([self.fit_side.prior.residual(M1.prior_log10_rate(point.x(),
                                                                                         self.fit_side.prior))]))
        M1.FAILURES.extend(failures)
        return np.concatenate(parts)

    def residuals_many(self, points, rows=None):
        return [self.residuals(t) for t in points]

    def point_r0(self, base_theta):
        keys, self.keys = self.keys, self.keys[:len(base_theta)]
        try:
            return M1.Objective.point(self, base_theta)
        finally:
            self.keys = keys


def run_fit(rival: str, start: int, freeze_commit: str, budget: int = M1.BUDGET) -> Dict[str, Any]:
    problems = fit_refusals(rival, freeze_commit)
    if problems:
        raise SystemExit("--fit refuses:\n  " + "\n  ".join(problems))
    commit = M1.check_freeze(freeze_commit, prereg=PREREG)
    m1 = json.loads(M1_CHECKS.read_text(encoding="utf-8"))
    numerics, fit_side = M1.numerics_for_fits(m1), M1.fit_side_from_checks(m1)
    br = tuple(read_record("checks")["davidek_breakpoints_in_situ"])
    obj = RivalObjective(rival, numerics, fit_side, br)
    out = V / f"{STEM}_fit_{rival}_s{start}.json"
    state = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {}
    if state.get("complete"):
        print(f"{out.name}: complete already")
        return state
    if int(state.get("failed_pot_runs", 0)) > 0:
        raise SystemExit(f"{out.name}: failed_pot_runs {state['failed_pot_runs']}: move it aside and re-run")
    carried = M1.lowest_converged("R0")
    base = (obj.theta_of(M1.b60_point()) if carried is None else
            np.array([carried["theta"].get(k, v) for k, v in zip(obj.keys[:17], obj.theta_of(M1.b60_point()))]))
    centre = np.array([0.5 * (lo + hi) for lo, hi in zip(obj.lower[17:], obj.upper[17:])])
    theta0 = np.concatenate([base, centre])
    if start == 1:
        theta0, _seed = M1.perturbed(theta0, obj.keys, obj.lower, obj.upper, 1)
    elif start == 2:
        theta0 = np.concatenate([obj.theta_of(M1.b60_point()), obj.lower[17:]])
    del M1.FAILURES[:]
    f_ = fitting.fit(obj.problem(), np.clip(theta0, obj.lower, obj.upper), budget=budget)
    state = {"artifact": out.stem, "rival": rival, "start": start, "freeze_commit": commit, "keys": list(obj.keys),
             "status": int(f_.status), "cost": float(f_.cost), "theta": dict(zip(obj.keys, map(float, f_.theta))),
             "failed_pot_runs": len(M1.FAILURES), "complete": True, "provenance": _provenance(f"fit_{rival}_s{start}")}
    out.write_text(provenance.dumps_compact(_jsonable(state), default=float) + "\n", encoding="utf-8")
    return state


# ----------------------------------------------------------------------------------------------------------------
# Markdown
# ----------------------------------------------------------------------------------------------------------------
def _fmt(x, nd=3):
    if x is None:
        return "--"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if abs(x) >= 1e-3 or x == 0 else f"{x:.2e}"
    return str(x)


def render_checks(rec: Mapping[str, Any]) -> str:
    c = rec["checks"]
    lines = ["# M1b part (i): verification before any score (V1 - V7, BIT, V-density)", "",
             f"*Generated by `{GENERATED_BY} --checks`; prereg `{rec['prereg']}` sec. 3.4. Overall: "
             f"**{'PASSED' if rec['passed'] else 'FAILED'}**.*", "",
             "| check | passed | what |", "|---|---|---|"]
    for k, v in c.items():
        lines.append(f"| {k} | {v.get('passed')} | {str(v.get('what', '')).replace('|', '/')[:160]} |")
    v6 = c.get("V6", {})
    if v6:
        lines += ["", "## V6: CF-084's numbers at B60 (B0-pinned, cooled 25 C, inline)", "",
                  "| control | B0 at B60 | measured | CF-084 (B51) |", "|---|---|---|---|"]
        for p, d in v6["meynier_cys_drift"].items():
            lines.append(f"| Meynier Cys {p.split('pH')[1]} drift | {d:+.3f} | <= 0.2 | -0.79 / -0.65 / -0.39 (4.5 / 5.5 / 6.5) |")
        for t, pu, pf, m in zip(HOU["times_min"], v6["hou_unforced_cooled"], v6["hou_forced_cooled"], HOU["ph_cooled"]):
            lines.append(f"| Hou {t:g} min (unforced / forced) | {pu:.3f} / {pf:.3f} | {m:.2f} | -- |")
        for lab, za, zc, ma, mc in zip(ZHOU_LABELS, v6["zhou_arp_alone_final"], v6["zhou_arp_cys_final"],
                                       v6["zhou_measured"]["C-Zhou-ARP"], v6["zhou_measured"]["C-Zhou-Cys"]):
            lines.append(f"| Zhou label {lab:g}: ARP alone / + Cys | {za:.3f} / {zc:.3f} | {ma:.3f} / {mc:.3f} | -- |")
        lines += ["", f"Hou RMS by ARP pKa (unforced): {v6['hou_rms_by_arp_pka_unforced']}.",
                  f"Signs that differ from sec. 3.2: {v6['signs_that_differ'] or 'none'}."]
        for k, rows in v6["cao"].items():
            lines.append(f"Cao {k}: " + "; ".join(f"{r['t_min']:g} min {r['predicted']:.3f} vs {r['measured']:.3f}"
                                                  for r in rows))
    v7 = c.get("V7", {})
    if v7:
        lines += ["", f"V7: forced ARP max relative error B0 {_fmt(v7['B0']['forced_arp_max_rel_error'])}, B3 "
                      f"{_fmt(v7['B3']['forced_arp_max_rel_error'])}; DPO held total (B3) "
                      f"{_fmt(v7['B3']['dpo_total_rel_dev_max'])}; forcing off bit-identical: "
                      f"{v7['B0']['forcing_off_bit_identical_to_unforced']} / {v7['B3']['forcing_off_bit_identical_to_unforced']}."]
    b = c.get("BIT", {})
    if b:
        lines += ["", f"BIT: base {b['base_commit'][:10]}, every control pot identical: {b['passed']}."]
    v5 = c.get("V5", {})
    if v5:
        lines += ["", "V5: Davidek 2006a's breakpoints, label (25 C, 0.16 M phosphate) -> in situ (120 C): "
                  + "; ".join(f"{r['label_ph']:g} -> {r['in_situ_ph']:.3f}" for r in v5["davidek_breakpoints"]
                              ["breakpoints"]) + ". Every scored comparison is cooled (25 C) on both sides."]
    return "\n".join(lines) + "\n"


def render_screen(rec: Mapping[str, Any]) -> str:
    lines = ["# M1b part (i): the mechanism screen (sec. 3.5)", "",
             f"*Generated by `{GENERATED_BY} --screen`; B60, inline, rtol 1e-8; >= {SCREEN_DRAWS} scrambled Sobol "
             f"draws per rival over its integrated coordinates (both balance constants included). Failed pot runs: "
             f"{rec['failed_pot_runs']}{' -- REFUSED' if rec['refused_for_failed_pot_runs'] else ''}.*", "",
             "| rival | draws | sigma-Hou-level | sigma-Hou-shape | sigma-Mey | sigma-Zhou | all four | would-be kill |",
             "|---|---|---|---|---|---|---|---|"]
    for r, v in rec["per_rival"].items():
        pf = v["pass_fraction"]
        lines.append(f"| {r} | {v['draws']} | {pf['sigma_hou_level']:.3f} | {pf['sigma_hou_shape']:.3f} | "
                     f"{pf['sigma_mey']:.3f} | {pf['sigma_zhou']:.3f} | {pf['all']:.3f} | {v['would_be_killed']} |")
    lines += ["", "Unforced C-Hou (printed beside the decisive forced run):", ""]
    for r, v in rec["per_rival"].items():
        lines.append(f"- {r}: {v['pass_fraction_unforced_hou']}")
    lines += ["", "Coverage (sec. 3.5: reported where fewer than 10 draws pass):", ""]
    for r, v in rec["per_rival"].items():
        lines.append(f"- {r}: {v['coverage']['n_passing']} pass ({v['coverage']['fraction_passing']:.3f}); box "
                     f"{v['coverage']['bounding_box']}")
    lines += ["", "Best draw per rival (the power check's truth):", ""]
    for r, v in rec["per_rival"].items():
        bd = v["best_draw"]
        lines.append(f"- {r}: {json.dumps({k: round(x, 4) for k, x in bd['coordinates'].items()})}; Hou forced "
                     f"{[round(x, 2) for x in bd['hou_forced']]}; Mey Cys drift {[round(x, 2) for x in bd['mey_cys_drift']]}; "
                     f"Zhou ARP {[round(x, 2) for x in bd['zhou_arp']]}")
    lines += ["", "Meynier glycine pots, most negative drift per draw (quantiles): "
              + "; ".join(f"{r} {v['mey_gly_drift_quantiles']}" for r, v in rec["per_rival"].items()),
              "", "Declared readings: see the record's `declared_readings`.", ""]
    return "\n".join(lines) + "\n"


def render_power(rec: Mapping[str, Any]) -> str:
    lines = ["# M1b part (i): the synthetic power check (sec. 7)", "",
             f"*Generated by `{GENERATED_BY} --power`; R = {rec['R']} synthetic data sets per truth; tie bands from "
             f"{rec['tie_replicates']} measurement-noise replicates. {rec['note']}*", "",
             "| pair | power (95 % CI) | false credit (95 % CI) | P(both inadequate) | power, unresolved = tie | "
             "verdict |", "|---|---|---|---|---|---|"]
    for k, v in rec["pairs"].items():
        p, fc = v["power"], v["false_credit"]
        lines.append(f"| {k} | {p['rate']:.3f} ({p['wilson_95'][0]:.3f} - {p['wilson_95'][1]:.3f}) | "
                     f"{fc['rate']:.3f} ({fc['wilson_95'][0]:.3f} - {fc['wilson_95'][1]:.3f}) | "
                     f"{p['p_both_inadequate']:.3f} / {fc['p_both_inadequate']:.3f} | "
                     f"{p['rate_unresolved_density_as_tie']:.3f} | {v['verdict']} |")
    lines += ["", f"Demoted to 'undecidable on disk': {rec['demoted'] or 'none'}. Low power: {rec['low_power'] or 'none'}.",
              f"Forms in the power check: {rec['forms_in_power_check']}; not in it (need new draws, scored in "
              f"--score): {rec['forms_not_in_power_check']}.", ""]
    for k, v in rec["pairs"].items():
        lines.append(f"- {k}: forms breaking credit under truth a {v['power']['forms_breaking_credit']}; under "
                     f"truth b {v['false_credit']['forms_breaking_credit']}")
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--checks", action="store_true")
    mode.add_argument("--screen", action="store_true")
    mode.add_argument("--power", action="store_true")
    mode.add_argument("--score", action="store_true")
    mode.add_argument("--fit", metavar="RIVAL")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--freeze-commit", default="")
    ap.add_argument("--base", default=None, help="--checks: the base commit of the bit-identity check")
    ap.add_argument("--min-free-percent", type=float, default=MIN_FREE_PERCENT)
    args = ap.parse_args(argv)
    if args.checks:
        run_checks(args.base)
    elif args.screen:
        run_screen(args.min_free_percent)
    elif args.power:
        run_power()
    elif args.score:
        run_score(args.freeze_commit, args.min_free_percent)
    else:
        run_fit(args.fit, args.start, args.freeze_commit)
    return 0


if __name__ == "__main__":
    sys.exit(main())

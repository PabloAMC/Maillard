"""
src/m3_strecker_layer.py -- M3 PHASE 1: the standalone forced-input Strecker layer (one shared ARP-Phe decay, two
routes), its rivals, forms, units, noise model and nulls.

Design: results/validation/kinetic_core_m3_strecker_layer_prereg.md (v2 DRAFT, 8dda42e5). Section numbers below are
that file's. This module is NEW code: it imports nothing from the shipped lanes except B18's frozen constants and the
pyrazine pH factor for R0-pinned (V4), and it extends EX-049's forced layer
(scripts/generators/feasibility_strecker_forced_input.py, frozen, not edited). Route A's calibration is EX-049's
arithmetic, re-implemented here and checked bit-exact against it (V1, tests/unit/test_m3_strecker_layer.py).

WHAT IS HERE
============
- sec. 3.1 route A (EX-049): exact inversion of hofmann2000 Tables 3 - 4 per draw (``calibrate_route_a``).
- sec. 3.2 the ARP-Phe exits (e1, e3, eF, eG, eFox, eB), shared by the fed pots and the glucose pot.
- sec. 3.3 oxygen: closed air vial (unlimited, primary), finite headspace budgets {5, 20, 100 mL} x nu_O2 {0.5, 1},
  argon residual o_Ar {0, 0.1}, one shared O2_reflux, Cu at 5 mM by exact inversion in [1, 25], Cu at 25 mM saturated
  or linear.
- sec. 3.4 the in-situ pH: an acid load against the phosphate buffer per pot (``buffer_ph``; V-acid), and the
  k_d(t) = k_d,7 10^(s_pH dpH) law (R4 everywhere; "pH drift on" in the trapped pot for the other rivals).
- sec. 4 the rivals R0-pinned, R1, R2 (widened: FOX, Voigt, trap efficiency from k0_DG1), R3, R4; the coordinates and
  their bands; sigma_ox and c_Cu5 solved per draw (non-circular feasibility).
- sec. 6 the units' value bases, measurement noise (analyser CV band, sigma_run between points, ARP isomers), omega,
  the shared-reference covariance, censoring; the nulls.

Every array function is vectorised over draws (axis 0). Times are minutes, concentrations mmol/L, logs dex unless a
name says ln. Nothing here reads hofmann2000b Figs. 2 - 3: the trapped pot runs on synthetic inputs until the
V-blind check passes (sec. 7; scripts/generators/fit_m3_strecker_layer.py).

DECLARED CHOICES WHERE THE PREREG IS SILENT (each is reported in the build record; none sets a band)
====================================================================================================
D-a  Past HPO4^2- exhaustion the pH is floored at 4.76 (acetic acid's pKa: the acids formed buffer there); flagged.
D-b  The ARP-made donors' sink is kappa_D [D][Phe] (sec. 3.1); the FED-DONOR calibration keeps EX-049's amine-free
     kappa_D [D]^2 (V1). The same kappa_D (L/(mmol min)) serves both.
D-c  A donor's Strecker share s_X = k_X / kappa_X. s_min / s_max are the draw's calibrated glucosone / MGO shares
     (sec. 4.3's "lowest / highest measured DFA share"); s_FOX log-uniform on [s_G, s_MGO]; s_DG1 on
     [s_G, s_MGO 10^0.82] (form CF-073: 10^1.24).
D-d  phi_FOX is drawn on the measured phi_D range [0.14, 0.30] like phi_DG1; FOX carries omega_C23, 1-DG omega_C6;
     kappa_FOX is set like kappa_DG1 (kappa_MGO 10^off, floor 7e-5).
D-e  At 85 C every rate (exits, sinks, k0_DG1, Strecker) carries the exits' Arrhenius factor (one Ea); the
     independent-Ea form applies a second Ea to the O2 exits only.
D-f  In the trapped pot acid comes from osone lost to non-trap channels plus the fragment, FOX and route-B exits:
     the dynamic form of sec. 3.4's (1 - eps_X) y_acid.
D-g  Finite O2 in the trapped pot uses the same headspace per 10 mL of solution as the closed vials.
D-h  Phe liberated (U3) is the cumulative Phe released by the exits (sec. 3.2), not net of route-A consumption.
D-i  U3's read error is the prereg's literal 0.011 umol/mmol (1.1e-5 mmol/mmol) plus 0.03 dex.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from src import data_access
from src.kinetic_core import fit_rows

# ======================================================================================================== constants

DT = 0.5                                   # min, the forced layer's grid (EX-049)
T_END = 500.0
GRID = np.arange(0.0, T_END + DT / 2, DT)
D0 = AA0 = 100.0                           # mmol/L, every fed pot (1 mmol in 10 mL)
DONORS = ("MGO", "GO", "TDG", "G")
DONOR_CLASS = {"MGO": "C23", "GO": "C23", "TDG": "C6", "G": "C6"}
OMEGA_CALIBRATOR = {"C23": "MGO", "C6": "TDG"}
KAPPA_FLOOR = {"GO": 7e-5, "G": 1.7e-4, "TDG": 7e-5, "DG1": 7e-5, "FOX": 7e-5}
KAPPA_TDG_MARTINS = 8.12e-4
TABLES = ("strecker_hofmann2000", "strecker_hofmann2000_fig2", "strecker_hofmann2000b", "strecker_hofmann2000b_fig1")
#: the keyed trapped-pot file: NEVER read by this module (sec. 7); the generator's V-blind check guards it
KEYED_TRAPPED_TABLE = "strecker_hofmann2000b_fig2_fig3"

R_GAS = 8.314462618e-3                     # kJ/(mol K)
T100, T85 = 373.15, 358.15
# sec. 3.4: phosphate; the acid titrates HPO4^2- (pKa2' 6.79: HPO4 fraction 0.62 at pH 7.0 in 0.5 M)
PKA2 = 6.79
PH_LABEL = 7.0
PH_SCALE = "label_25C"
PH_ACID_FLOOR = 4.76                       # D-a
PHOSPHATE_M = 0.5
PHOSPHATE_UNBUFFERED_M = 0.0005            # the Methods' "0.5 mmol/L" reading (form)
# sec. 3.3: O2 per mL of headspace air at 25 C, and the dissolved O2 of 10 mL
O2_UMOL_PER_ML_AIR = 0.2095 / 24.465 * 1000.0
O2_DISSOLVED_UMOL_10ML = 2.6
SOLUTION_ML = 10.0
HEADSPACES_ML = (5.0, 20.0, 100.0)
NU_O2 = (0.5, 1.0)

# sec. 4: bands
LOG10_KD_BAND = (-1.95, -0.88)
LOG10_KD_FLOOR = -2.30
S_PH_BAND = (0.20, 0.54)
Y_ACID_BAND = (0.4, 1.1)
EA_BAND = (65.0, 85.0)
LOG10_K_TR_BAND = (-1.0, 1.0)              # k_tr [0.1, 10] /min, log-uniform (numerical ends)
LOG10_K0_DG1_BAND = (math.log10(0.07), 0.0)
PHI_DG1_BAND = (0.14, 0.30)
GRANVOGL_DEX = 0.82
CF073_DEX = 1.24
VOIGT_BAND = (1.5, 1.6)
C_CU_BOUND = (1.0, 25.0)
ORDER_BAND = (1.0, 2.0)
ORDER_REF = 100.0                          # mM: rate = k_d (ARP / 100 mM)^(n - 1) ARP
CV_ARP_BAND = (0.02, 0.10)
F_ISO_BAND = (0.0, 0.3)

# sec. 6: error model
OMEGA_DEX = 0.419
SIGMA_RUN_DEX = 0.075
S_REP_DEX = 0.025
ROUND_UMOL_MMOL_SD = 0.05 / math.sqrt(3.0)          # +-0.05 umol/mmol uniform: sd 0.029
U3_READ_ABS = 0.011e-3                              # D-i, mmol/mmol
U3_READ_DEX = 0.03
U3_CENSOR = 0.004                                   # mmol/mmol, the synthetic design's axis resolution (sec. 12.1)
U3_CAL_SD = {"DG1": 0.20, "TDG": 0.10, "G": 0.10, "PHE": 0.05}
U3_TIMES = (20.0, 40.0, 80.0, 160.0, 320.0)
U3_SPECIES = ("DG1", "TDG", "G", "PHE")
U3_ARMS = ("argon", "air_cu25")
TRAP_SCALE = 5.0                                    # 0.5 M ARP-Phe against the 0.1 M pots
TRAP_T_END = 320.0

RIVALS = ("R0", "R1", "R2", "R3", "R4")


# ============================================================================================================ forms

@dataclass(frozen=True)
class Form:
    """One setting of every unsourced choice and declared transfer (secs. 3 - 6, 9.1). The primary is Form().
    Layer-changing fields change the draws' predictions; likelihood-only fields only reweight (sec. 9.1)."""
    # layer-changing
    o_ar: float = 0.0
    headspace_ml: Optional[float] = None          # None = closed vial, unlimited O2 (primary)
    nu_o2: float = 1.0
    cu_law: str = "saturated"                     # trapped pot at 25 mM: "saturated" | "linear"
    t0: float = 0.0
    cf070: bool = False                           # constant Phe in the glucose pot
    early_step: bool = False
    independent_ea: bool = False
    cf073: bool = False
    voigt: bool = False
    ph_drift_trapped: bool = False                # non-R4 rivals: the pH term in the trapped pot
    ph_drift_fed: bool = False                    # diagnostic only: R1's drift in the 0.1 M pots
    unbuffered_trapped: bool = False
    kd_floor: bool = False
    isomers: bool = False
    per_pot_o2: bool = False                      # R3 verdicts: the glucose pot's O2 separate from the fed reflux
    kappa_ceiling: float = 0.4
    # likelihood-only
    omega_scale: float = 1.0
    rho: float = 0.0
    censor_scale: float = 1.0
    cv_arp: Optional[float] = None                # None: integrated over the band (a coordinate); else fixed
    dg1_cal_sd: float = 0.20
    run_between_points: bool = False

    def label(self) -> str:
        base = Form()
        diff = {k: v for k, v in self.__dict__.items() if v != getattr(base, k)}
        return "primary" if not diff else ",".join(f"{k}={v}" for k, v in sorted(diff.items()))


PRIMARY = Form()
LIKELIHOOD_ONLY_FIELDS = ("omega_scale", "rho", "censor_scale", "cv_arp", "dg1_cal_sd", "run_between_points")


def layer_part(form: Form) -> Form:
    """The layer-changing part of a form: two forms with equal layer parts share one draw set (reweighting)."""
    return replace(form, **{k: getattr(PRIMARY, k) for k in LIKELIHOOD_ONLY_FIELDS})


# =================================================================================================== coordinates

@dataclass(frozen=True)
class Coord:
    name: str
    lo: float
    hi: float
    basis: str
    rivals: Tuple[str, ...] = ("R1", "R2", "R3", "R4")


_ALL = ("R1", "R2", "R3", "R4")
COORDS_BASE: Tuple[Coord, ...] = (
    # route A (EX-049's bands, forced prereg sec. 2; EX-042 review)
    Coord("log10_c_mgo", -2.5, 0.5, "kappa_MGO D0 /min (EX-049)"),
    Coord("nu", 0.1, 1.0, "amine per donor consumed (EX-049)"),
    Coord("log10_k_loss", -5.0, -2.0, "PA loss /min (EX-049)"),
    Coord("off_kappa_go", -0.5, 0.4, "log10 kappa_GO / kappa_MGO (EX-049; +0.5 form)"),
    Coord("off_kappa_g", -0.5, 0.4, "log10 kappa_G / kappa_MGO (EX-049; +0.5 form)"),
    Coord("off_kappa_tdg", -0.5, 0.4, "log10 kappa_TDG / Martins 8.12e-4 (EX-049; +0.5 form)"),
    Coord("cal_mgo", -0.30, 0.30, "MGO forcing calibration dex (EX-049)", ("R0",) + _ALL),
    Coord("cal_go", -0.40, 0.40, "GO forcing calibration dex (EX-049)", ("R0",) + _ALL),
    Coord("cal_tdg", -0.20, 0.20, "3-DG forcing calibration dex (EX-049)"),
    Coord("cal_g", -0.20, 0.20, "glucosone forcing calibration dex (EX-049)"),
    Coord("cal_arp", -0.20, 0.20, "ARP-Phe course calibration dex (EX-049)"),
    # sec. 4.1
    Coord("o2_reflux", 0.0, 1.0, "one O2 for every reflux pot (sec. 3.3)"),
    Coord("arp0", 90.0, 100.0, "fed ARP-Phe mM (purity 'about 90 %')"),
    Coord("imp_phe", 0.0, 1.0, "free-Phe fraction of the ~10 % impurity (definition)"),
    Coord("u_a1", 0.0, 1.0, "anoxic shares a1, a3, aF ~ Dirichlet(1,1,1): stick-breaking 1"),
    Coord("u_a2", 0.0, 1.0, "anoxic shares: stick-breaking 2"),
    Coord("u_g1", 0.0, 1.0, "O2 split g_G, g_F, g_B ~ Dirichlet(1,1,1) (R2: g_G, g_F ~ Dirichlet(1,1))"),
    Coord("u_g2", 0.0, 1.0, "O2 split: stick-breaking 2 (unused by R2)", ("R1", "R3", "R4")),
    Coord("r_f", 0.0, 1.0, "Phe-releasing fraction of the anoxic fragment exit (definition)"),
    Coord("m_frag", 0.0, 1.0, "MGO fraction of the anoxic fragments (definition)"),
    Coord("log10_kd7", LOG10_KD_BAND[0], LOG10_KD_BAND[1], "k_d at label pH 7.0, air, 100 C (sec. 4.2)"),
    Coord("s_ph", S_PH_BAND[0], S_PH_BAND[1], "dec / pH (Deng 100 C; Martins drift-corrected)"),
    Coord("y_acid", Y_ACID_BAND[0], Y_ACID_BAND[1], "mol acid / mol ARP decayed (sec. 2 item 11)"),
    Coord("ea", EA_BAND[0], EA_BAND[1], "kJ/mol, the exits (Deng 2021)"),
    Coord("log10_k_tr", LOG10_K_TR_BAND[0], LOG10_K_TR_BAND[1], "trap rate /min at 0.6 M (numerical ends)"),
    Coord("u_s_fox", 0.0, 1.0, "s_FOX log-uniform on the measured DFA range (D-c)"),
    Coord("u_s_dg1", 0.0, 1.0, "s_DG1 log-uniform on [s_G, s_MGO 10^0.82] (sec. 4.3)"),
    Coord("phi_dg1", PHI_DG1_BAND[0], PHI_DG1_BAND[1], "1-DG acid share (measured phi_D range)"),
    Coord("phi_fox", PHI_DG1_BAND[0], PHI_DG1_BAND[1], "FOX acid share (D-d)"),
    Coord("off_kappa_dg1", -0.5, 0.4, "log10 kappa_DG1 / kappa_MGO (as the unmeasured donors)"),
    Coord("off_kappa_fox", -0.5, 0.4, "log10 kappa_FOX / kappa_MGO (D-d)"),
    Coord("log10_k0_dg1", LOG10_K0_DG1_BAND[0], LOG10_K0_DG1_BAND[1], "1-DG first-order loss /min (Davidek 2006a)"),
    Coord("psi", 0.0, 1.0, "route B acid share (definition)", ("R1", "R3", "R4")),
    Coord("cv_arp", CV_ARP_BAND[0], CV_ARP_BAND[1], "Table 3 analyser CV, integrated (sec. 6)"),
)
COORDS_FORM: Tuple[Coord, ...] = (
    Coord("ea_ox", EA_BAND[0], EA_BAND[1], "O2 exits' Ea (independent-Ea form)"),
    Coord("voigt_f", VOIGT_BAND[0], VOIGT_BAND[1], "Voigt 2010 O2 factor on s_DG1 (form)"),
    Coord("f_iso", F_ISO_BAND[0], F_ISO_BAND[1], "isomer share of the measured ARP at 300 min (form)"),
    Coord("o2_glc", 0.0, 1.0, "the glucose pot's own O2 (per-pot form, R3 verdicts)"),
)
COORD_R3 = Coord("log10_kd7_glc", LOG10_KD_BAND[0], LOG10_KD_BAND[1], "R3: the glucose pot's k_d,7", ("R3",))
COORD_R4 = Coord("n_order", ORDER_BAND[0], ORDER_BAND[1], "R4: decay order (Martins / Deng 1; Wei 2022 2)", ("R4",))


def coords(rival: str, form: Form = PRIMARY) -> Tuple[Coord, ...]:
    """The integrated coordinates of ``rival`` under ``form``, in a fixed order (the unit cube's axes)."""
    out: List[Coord] = []
    for c in COORDS_BASE + (COORD_R3, COORD_R4):
        if rival not in c.rivals:
            continue
        if c.name == "log10_kd7" and form.kd_floor:
            c = replace(c, lo=LOG10_KD_FLOOR)
        if c.name == "log10_kd7_glc" and form.kd_floor:
            c = replace(c, lo=LOG10_KD_FLOOR)
        if c.name.startswith("off_kappa_") and form.kappa_ceiling != 0.4:
            c = replace(c, hi=form.kappa_ceiling)
        if c.name == "cv_arp" and form.cv_arp is not None:
            continue
        out.append(c)
    if rival != "R0":
        extra = {"ea_ox": form.independent_ea, "voigt_f": form.voigt, "f_iso": form.isomers,
                 "o2_glc": form.per_pot_o2 and rival == "R3"}
        out += [c for c in COORDS_FORM if extra[c.name]]
    return tuple(out)


def from_unit(u: np.ndarray, rival: str, form: Form = PRIMARY) -> Dict[str, np.ndarray]:
    """(n, d) unit-cube points -> {coordinate: (n,) values}. Uniform on every band (log10 coordinates are bands of
    the log)."""
    cs = coords(rival, form)
    u = np.atleast_2d(np.asarray(u, dtype=float))
    if u.shape[1] != len(cs):
        raise ValueError(f"{rival}: expected {len(cs)} coordinates, got {u.shape[1]}")
    th = {c.name: c.lo + (c.hi - c.lo) * u[:, i] for i, c in enumerate(cs)}
    if form.cv_arp is not None:
        th["cv_arp"] = np.full(u.shape[0], float(form.cv_arp))
    return th


def sobol_unit(n: int, d: int, seed: int) -> np.ndarray:
    from scipy.stats import qmc
    m = int(math.ceil(math.log2(max(n, 2))))
    return qmc.Sobol(d, scramble=True, seed=seed).random_base2(m)[:n]


def stick3(u1: np.ndarray, u2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dirichlet(1, 1, 1) by stick-breaking: x1 ~ Beta(1, 2) by inversion, x2 = (1 - x1) u2."""
    x1 = 1.0 - np.sqrt(1.0 - u1)
    x2 = (1.0 - x1) * u2
    return x1, x2, 1.0 - x1 - x2


# ============================================================================================================ data

def _read_table(name: str) -> List[dict]:
    """EX-049's literature-only read: status 'ready', every pot verified by two readers."""
    if name == KEYED_TRAPPED_TABLE:
        raise PermissionError(f"{name} is read only by the scorer after V-blind passes (prereg sec. 7)")
    table = data_access.load_yaml(fit_rows.FIT_ROWS_DIR / f"{name}.yml")
    if table.get("status") != "ready":
        raise PermissionError(f"{name}: status {table.get('status')!r}, not 'ready'")
    verified = fit_rows._verified(name)
    for pot_id, pot in (table.get("pots") or {}).items():
        if not verified.get(pot["ledger_id"], False):
            raise PermissionError(f"{name}: pot {pot_id}'s source {pot['ledger_id']} is not verified by two readers")
    return [dict(r, sigma=fit_rows._sigma(r)) for r in table.get("rows") or []]


@dataclass(frozen=True)
class Data:
    """Every observation phase 1 uses (mmol/L, read errors dex), by role."""
    t3: Dict[str, Dict[float, float]]                        # hofmann2000 Table 3 MGO pot (G-T3)
    fed: Dict[str, Dict[str, Tuple[float, float]]]           # Table 4 donors, argon / air (PA, PAA) at 60 min
    fed_cu_mgo: Tuple[float, float]                          # Table 4 MGO air + Cu (null only)
    arp_t2: Dict[str, Tuple[float, float]]                   # hofmann2000b Table 2 fed ARP: argon / air / air_cu5
    arp_t2_err: Dict[str, Tuple[float, float]]
    fig1_arp: Dict[str, Dict[float, Tuple[float, float]]]    # Fig. 1 circles {PHEAL|PHEAC: {t: (value, err)}}
    fig1_glc: Dict[str, Dict[float, Tuple[float, float]]]    # Fig. 1 squares (U5b)
    fig2: Dict[str, Dict[float, Tuple[float, float]]]        # hofmann2000 Fig. 2 (U5a)
    arp_glc: Dict[float, float]                              # hofmann2000b Table 3 ARP-Phe (U-ARPglc)
    forced_inputs: Dict[str, Tuple[np.ndarray, np.ndarray]]


def load_data() -> Data:
    rows: Dict[Tuple[str, str, float], dict] = {}
    for name in TABLES:
        for r in _read_table(name):
            rows[(r["pot"], r["species"], float(r["t_min"]))] = r

    def v(pot, sp, t):
        return float(rows[(pot, sp, float(t))]["value"])

    def ve(pot, sp, t):
        r = rows[(pot, sp, float(t))]
        return float(r["value"]), float(r["sigma"])

    t3pot = "hofmann2000_t3_mgo_phe_98c_closed"
    t3 = {sp: {t: v(t3pot, sp, t) for t in (10, 30, 60, 120, 240, 300, 500)} for sp in ("PHEAL", "PHEAC")}
    fed = {d: {atm: (v(f"hofmann2000_t4_{d.lower()}_phe_60min_{atm}", "PHEAL", 60),
                     v(f"hofmann2000_t4_{d.lower()}_phe_60min_{atm}", "PHEAC", 60)) for atm in ("argon", "air")}
           for d in DONORS}
    fed_cu = (v("hofmann2000_t4_mgo_phe_60min_air_cu", "PHEAL", 60), v("hofmann2000_t4_mgo_phe_60min_air_cu",
                                                                       "PHEAC", 60))
    t2 = {}
    t2e = {}
    for atm in ("argon", "air", "air_cu"):
        pot = f"hofmann2000b_t2_arp_phe_100c_120min_{atm}"
        key = "air_cu5" if atm == "air_cu" else atm
        (a, ea), (b, eb) = ve(pot, "PHEAL", 120), ve(pot, "PHEAC", 120)
        t2[key], t2e[key] = (a, b), (ea, eb)

    def course(pot, times):
        return {sp: {float(t): ve(pot, sp, t) for t in times} for sp in ("PHEAL", "PHEAC")}

    f1 = (10, 30, 60, 120, 240, 300, 500)
    fig1_arp = course("hofmann2000b_fig1_arp_phe_reflux", f1)
    fig1_glc = course("hofmann2000b_fig1_glc_phe_reflux", f1)
    fig2 = course("hofmann2000_fig2_glc_phe_98c", (10, 30, 60, 120, 240, 500))
    t3b = "hofmann2000b_t3_glc_phe_reflux"
    arp_glc = {float(t): v(t3b, "PHEARP", t) for t in (10, 30, 60, 120, 300)}
    # the forcing: EX-049's construction, verbatim (V1)
    f3pot = "hofmann2000_fig3_glc_phe_98c"
    forced = {}
    for d in DONORS:
        ts = [0.0] + [float(t) for t in (10, 30, 60, 120, 240, 300, 500) if (f3pot, d, float(t)) in rows]
        vs = [0.0] + [v(f3pot, d, t) for t in ts[1:]]
        if d == "G":
            ts += [300.0, 500.0]
            vs += [0.0, 0.0]
        forced[d] = (np.array(ts), np.array(vs))
    for sp, start in (("PHE", AA0), ("PHEARP", 0.0)):
        ts = [0.0, 10.0, 30.0, 60.0, 120.0, 300.0]
        vs = [start] + [v(t3b, sp, t) for t in ts[1:]]
        slope = math.log(vs[-1] / vs[-2]) / (ts[-1] - ts[-2])
        forced[sp] = (np.array(ts + [500.0]), np.array(vs + [vs[-1] * math.exp(slope * 200.0)]))
    return Data(t3=t3, fed=fed, fed_cu_mgo=fed_cu, arp_t2=t2, arp_t2_err=t2e, fig1_arp=fig1_arp, fig1_glc=fig1_glc,
                fig2=fig2, arp_glc=arp_glc, forced_inputs=forced)


# ===================================================================================== route A (EX-049, V1 exact)

def idx(t: float) -> int:
    return int(round(t / DT))


def cum_trapz(f: np.ndarray) -> np.ndarray:
    out = np.zeros_like(f)
    out[..., 1:] = np.cumsum(0.5 * DT * (f[..., 1:] + f[..., :-1]), axis=-1)
    return out


def cum_with_loss(f: np.ndarray, k_loss) -> np.ndarray:
    k = np.asarray(k_loss, dtype=float).reshape(-1, 1) if f.ndim == 2 else float(k_loss)
    return np.exp(-k * GRID) * cum_trapz(f * np.exp(k * GRID))


def fed_donor_courses(kappa: np.ndarray, nu: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    d = D0 / (1.0 + kappa[:, None] * D0 * GRID[None, :])
    aa = np.maximum(AA0 - nu[:, None] * (D0 - d), 0.0)
    return d, aa


def interpolate(forced: Mapping[str, Tuple[np.ndarray, np.ndarray]], interp: str = "linear") -> Dict[str, np.ndarray]:
    out = {}
    for key, (ts, vs) in forced.items():
        if interp == "early_step" and key != "PHE":
            ts, vs = np.concatenate([[0.0, DT], ts[1:]]), np.concatenate([[0.0, vs[1]], vs[1:]])
        out[key] = np.interp(GRID, ts, vs)
    return out


def kappas(th: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    k_mgo = 10.0 ** th["log10_c_mgo"] / D0
    out = {"MGO": k_mgo,
           "GO": np.maximum(k_mgo * 10.0 ** th["off_kappa_go"], KAPPA_FLOOR["GO"]),
           "G": np.maximum(k_mgo * 10.0 ** th["off_kappa_g"], KAPPA_FLOOR["G"]),
           "TDG": np.maximum(KAPPA_TDG_MARTINS * 10.0 ** th["off_kappa_tdg"], KAPPA_FLOOR["TDG"])}
    if "off_kappa_dg1" in th:
        out["DG1"] = np.maximum(k_mgo * 10.0 ** th["off_kappa_dg1"], KAPPA_FLOOR["DG1"])
        out["FOX"] = np.maximum(k_mgo * 10.0 ** th["off_kappa_fox"], KAPPA_FLOOR["FOX"])
    return out


def calibrate_route_a(th: Mapping[str, np.ndarray], data: Data) -> Dict[str, Any]:
    """EX-049's exact inversion (``calibrate`` there, route-A part): k_D, phi_D from Table 4 argon; omega per class
    from the class calibrator's air row; the Table 3 courses for G-T3 (T1, T1b). Bit-exact against EX-049 (V1)."""
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
    cal["t3_pa"] = (1 - phi) * k * lossy
    cal["t3_total"] = cal["t3_pa"] + (phi + om) * k * integ
    return cal


# ================================================================================================ pH (sec. 3.4)

def buffer_ph(acid_m, phosphate_m: float = PHOSPHATE_M, ph0: float = PH_LABEL, pka2: float = PKA2,
              floor: float = PH_ACID_FLOOR):
    """The in-situ pH on the 25 C label scale (rule 13) after ``acid_m`` mol/L of (fully dissociated) acid has
    titrated HPO4^2- -> H2PO4^- in a phosphate buffer at label pH ``ph0`` (Henderson-Hasselbalch at pKa2'). Past
    HPO4^2- exhaustion the pH sits at ``floor`` (D-a). Vectorised; returns (pH, exhausted mask)."""
    a = np.asarray(acid_m, dtype=float)
    r0 = 10.0 ** (ph0 - pka2)
    hpo4 = phosphate_m * r0 / (1.0 + r0) - a
    h2po4 = phosphate_m / (1.0 + r0) + a
    exhausted = hpo4 <= 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        ph = pka2 + np.log10(np.where(exhausted, 1.0, hpo4) / h2po4)
    ph = np.where(exhausted, floor, np.maximum(ph, floor))
    return ph, exhausted | (ph <= floor)


def ph_factor(ph, s_ph) -> np.ndarray:
    """k_d(pH) / k_d,7 = 10^(s_pH (pH - 7.0)), both on the label scale."""
    return 10.0 ** (np.asarray(s_ph) * (np.asarray(ph) - PH_LABEL))


def v_acid_table() -> List[Dict[str, Any]]:
    """Sec. 3.4's acid-load table, recomputed by ``buffer_ph`` (V-acid). Each row: the prereg's printed pH drop(s)
    and the layer's. The pots' decayed amounts are the prereg's (81 % / 100 % of 0.1 M; 50 - 55 % of 0.5 M; Phe
    consumed 49 / 60 mM as the sugar proxy)."""
    rows = []

    def drop(p, acid):
        ph, ex = buffer_ph(acid, p)
        return float(PH_LABEL - ph), bool(ex)

    spec = [
        ("fed ARP 0.1 M, 120 min, y 0.4 / 1.1", 0.5, (0.4 * 0.081, 1.1 * 0.081), (0.12, 0.31)),
        ("fed ARP 0.1 M, 500 min, y 0.4 / 1.1", 0.5, (0.4 * 0.100, 1.1 * 0.100), (0.14, 0.39)),
        ("trapped 0.5 M, 320 min, y 0.4 at 50 %", 0.5, (0.4 * 0.25,), (0.35,)),
        ("trapped 0.5 M, 320 min, y 1.0 at 50 / 55 %", 0.5, (1.0 * 0.25, 1.0 * 0.275), (1.1, 1.3)),
        ("glucose pot, 120 min, y 0.4 / 1.1 (Phe consumed 49 mM)", 0.5, (0.4 * 0.049, 1.1 * 0.049), (0.07, 0.19)),
        ("glucose pot, 300 min, y 0.4 / 1.1 (Phe consumed 60 mM)", 0.5, (0.4 * 0.060, 1.1 * 0.060), (0.09, 0.23)),
    ]
    # the prereg's standalone script (m3_acid_load.py, run through the wrapper) printed these drops to 2 decimals;
    # the sec. 3.4 table rounded its 1.34 (y 1.0, 55 %) to "1.3". V-acid compares with the script's value where it
    # printed one, else with the table's.
    script = {("trapped 0.5 M, 320 min, y 1.0 at 50 / 55 %", 1): 1.34,
              ("trapped 0.5 M, 320 min, y 1.0 at 50 / 55 %", 0): 1.08,
              ("fed ARP 0.1 M, 120 min, y 0.4 / 1.1", 1): 0.31, ("fed ARP 0.1 M, 500 min, y 0.4 / 1.1", 1): 0.39,
              ("glucose pot, 120 min, y 0.4 / 1.1 (Phe consumed 49 mM)", 1): 0.19,
              ("glucose pot, 300 min, y 0.4 / 1.1 (Phe consumed 60 mM)", 1): 0.23}
    for name, p, acids, printed in spec:
        got = [drop(p, a) for a in acids]
        ref = [script.get((name, i), pr) for i, pr in enumerate(printed)]
        rows.append({"pot": name, "phosphate_M": p, "acid_mM": [round(1000 * a, 2) for a in acids],
                     "printed_drop": list(printed), "script_drop": [script.get((name, i)) for i in range(len(acids))],
                     "layer_drop": [round(g[0], 4) for g in got],
                     "abs_diff": [round(abs(g[0] - r), 4) for g, r in zip(got, ref)],
                     "abs_diff_vs_printed": [round(abs(g[0] - pr), 4) for g, pr in zip(got, printed)],
                     "exhausted": [g[1] for g in got]})
    # y 1.1 at 55 %: "HPO4 nearly exhausted (pH <= ~5)"
    ph, ex = buffer_ph(1.1 * 0.275, 0.5)
    rows.append({"pot": "trapped 0.5 M, y 1.1 at 55 %", "phosphate_M": 0.5, "acid_mM": [302.5],
                 "printed_drop": ["pH <= ~5"], "layer_ph": round(float(ph), 3), "pass": bool(float(ph) <= 5.05),
                 "exhausted": [bool(ex)]})
    # Martins check: 6.8 mM acids, 0.1 M phosphate at pH 6.8: predicted drop 0.12 (measured 0.19)
    ph, _ = buffer_ph(0.0068, 0.1, ph0=6.8)
    rows.append({"pot": "Martins fed DFG 10 mM, 0.1 M phosphate pH 6.8, 6.8 mM acids", "phosphate_M": 0.1,
                 "acid_mM": [6.8], "printed_drop": [0.12], "layer_drop": [round(6.8 - float(ph), 4)],
                 "abs_diff": [round(abs(6.8 - float(ph) - 0.12), 4)], "exhausted": [False]})
    return rows


# ========================================================================================= derived parameters

def o2_capacity_mM(headspace_ml: float) -> float:
    """sec. 3.3: headspace air plus the dissolved O2, per litre of the 10 mL solution (umol/mL = mM)."""
    return (headspace_ml * O2_UMOL_PER_ML_AIR + O2_DISSOLVED_UMOL_10ML) / SOLUTION_ML


def arrhenius(ea_kj: np.ndarray, t_k: float) -> np.ndarray:
    return np.exp(-np.asarray(ea_kj) / R_GAS * (1.0 / t_k - 1.0 / T100))


@dataclass
class Params:
    """A draw's derived constants (arrays over draws), before sigma_ox / c_Cu5 are solved."""
    rival: str
    form: Form
    th: Dict[str, np.ndarray]
    cal: Dict[str, Any]
    n: int
    k_d7: np.ndarray
    k_d7_glc: np.ndarray
    a: Tuple[np.ndarray, np.ndarray, np.ndarray]           # a1, a3, aF
    g: Tuple[np.ndarray, np.ndarray, np.ndarray]           # g_G, g_F, g_B
    donor: Dict[str, Dict[str, np.ndarray]]                # X -> {k, phi, omega, kappa, k0}
    n_order: np.ndarray
    ph_on_fed: bool
    ph_on_trap: bool
    ph_on_glc: bool
    sigma_ox: Optional[np.ndarray] = None
    c_cu5: Optional[np.ndarray] = None
    feasible: Optional[np.ndarray] = None
    infeasible_reason: Dict[str, np.ndarray] = field(default_factory=dict)


def derive(th: Mapping[str, np.ndarray], rival: str, form: Form, data: Data) -> Params:
    th = {k: np.asarray(v, dtype=float) for k, v in th.items()}
    n = len(th["nu"])
    cal = calibrate_route_a(th, data)
    k_d7 = 10.0 ** th["log10_kd7"]
    k_d7_glc = 10.0 ** th["log10_kd7_glc"] if rival == "R3" else k_d7
    a = stick3(th["u_a1"], th["u_a2"])
    if rival == "R2":
        g_g = th["u_g1"]
        g = (g_g, 1.0 - g_g, np.zeros(n))
    else:
        g = stick3(th["u_g1"], th["u_g2"])
    kap = cal["kappa"]
    s_g, s_mgo = cal["k"]["G"] / kap["G"], cal["k"]["MGO"] / kap["MGO"]
    lo, hi = np.log10(np.minimum(s_g, s_mgo)), np.log10(np.maximum(s_g, s_mgo))
    s_fox = 10.0 ** (lo + th["u_s_fox"] * (hi - lo))
    hi_dg1 = hi + (CF073_DEX if form.cf073 else GRANVOGL_DEX)
    s_dg1 = 10.0 ** (lo + th["u_s_dg1"] * (hi_dg1 - lo))
    donor = {}
    for d in ("MGO", "TDG", "G"):
        donor[d] = {"k": cal["k"][d], "phi": cal["phi"][d], "omega": cal["omega"][DONOR_CLASS[d]],
                    "kappa": kap[d], "k0": np.zeros(n)}
    donor["DG1"] = {"k": s_dg1 * kap["DG1"], "phi": th["phi_dg1"], "omega": cal["omega"]["C6"], "kappa": kap["DG1"],
                    "k0": 10.0 ** th["log10_k0_dg1"]}
    donor["FOX"] = {"k": s_fox * kap["FOX"], "phi": th["phi_fox"], "omega": cal["omega"]["C23"], "kappa": kap["FOX"],
                    "k0": np.zeros(n)}
    n_order = th["n_order"] if rival == "R4" else np.ones(n)
    r4 = rival == "R4"
    return Params(rival=rival, form=form, th=th, cal=cal, n=n, k_d7=k_d7, k_d7_glc=k_d7_glc, a=a, g=g, donor=donor,
                  n_order=n_order, ph_on_fed=r4 or form.ph_drift_fed, ph_on_trap=r4 or form.ph_drift_trapped,
                  ph_on_glc=r4)


# ============================================================================================== the ARP pots

@dataclass(frozen=True)
class Pot:
    name: str
    atmosphere: str            # "argon" | "air" | "air_cu5" | "air_cu25" | "reflux"
    t_end: float
    out_times: Tuple[float, ...]
    scale: float = 1.0         # 1: 0.1 M ARP; 5: the trapped 0.5 M pot
    temp_k: float = T100
    trap: bool = False
    phosphate_m: float = PHOSPHATE_M


POT_T2 = {atm: Pot(f"t2_{atm}", atm, 120.0, (120.0,)) for atm in ("argon", "air", "air_cu5")}
POT_FIG1 = Pot("fig1_arp_reflux", "reflux", 500.0, (10.0, 30.0, 60.0, 120.0, 240.0, 300.0, 500.0))
POT_TRAP = {atm: Pot(f"trap_{atm}", atm, TRAP_T_END, U3_TIMES, scale=TRAP_SCALE, temp_k=T85, trap=True)
            for atm in U3_ARMS}
ARP_SPECIES = ("DG1", "TDG", "G", "FOX", "MGO")


def _o2_state(p: Params, pot: Pot, sigma_ox: np.ndarray, c_cu5: Optional[np.ndarray]):
    """(constant O2 level or None if budgeted, Cu factor, budget capacity mM or None)."""
    f = p.form
    atm = pot.atmosphere
    if atm == "argon":
        return np.full(p.n, f.o_ar), np.ones(p.n), None
    if atm == "reflux":
        return p.th["o2_reflux"], np.ones(p.n), None
    if atm == "air":
        cu = np.ones(p.n)
    elif atm == "air_cu5":
        cu = c_cu5
    else:
        cu = c_cu5 if f.cu_law == "saturated" else 1.0 + 5.0 * (c_cu5 - 1.0)
    if f.headspace_ml is None:
        return np.ones(p.n), cu, None
    return None, cu, o2_capacity_mM(f.headspace_ml)  # D-g: the same headspace per 10 mL


def solve_arp_pot(p: Params, pot: Pot, sigma_ox: np.ndarray, c_cu5: Optional[np.ndarray] = None,
                  ledger: bool = False) -> Dict[str, np.ndarray]:
    """The fed (or trapped) ARP-Phe pot: ARP decays through the exits of sec. 3.2 at the pot's O2, Cu, pH and
    temperature; released and impurity Phe feed route A on every ARP-made donor; route B makes PA / PAA directly.
    Exponential updates per step (stable for any rate on DT); every amount moved is booked (V2). Returns the
    observables at ``pot.out_times - t0`` (columns) and, with ``ledger``, the balance pools at the end."""
    n = p.n
    th, f = p.th, p.form
    t0 = f.t0 if not pot.trap else 0.0
    steps = int(round(pot.t_end / DT))
    out_idx = {int(round(max(t - t0, 0.0) / DT)): j for j, t in enumerate(pot.out_times)}
    fT = arrhenius(th["ea"], pot.temp_k)
    fTox = arrhenius(th["ea_ox"], pot.temp_k) if f.independent_ea else fT
    k_anox7 = (1.0 - sigma_ox) * p.k_d7
    k_ox7 = sigma_ox * p.k_d7
    a1, a3, aF = p.a
    gG, gF, gB = p.g
    psi = th.get("psi", np.zeros(n))
    o_const, cu, cap = _o2_state(p, pot, sigma_ox, c_cu5)
    budget = cap is not None
    ph_on = p.ph_on_trap if pot.trap else p.ph_on_fed
    phos = PHOSPHATE_UNBUFFERED_M if (pot.trap and f.unbuffered_trapped) else pot.phosphate_m
    s_ph, y_acid = th["s_ph"], th["y_acid"]
    voigt = th["voigt_f"] if (f.voigt and pot.atmosphere != "argon") else np.ones(n)
    k_tr = 10.0 ** th["log10_k_tr"] if pot.trap else np.zeros(n)
    k_loss = p.cal["k_loss"]
    arp0 = th["arp0"] * pot.scale
    arp = arp0.copy()
    phe = th["imp_phe"] * (100.0 - th["arp0"]) * pot.scale
    phe_start = phe.copy()
    X = {x: np.zeros(n) for x in ARP_SPECIES}
    TR = {x: np.zeros(n) for x in ARP_SPECIES}
    pa = np.zeros(n)
    paa = np.zeros(n)
    acid = np.zeros(n)
    o2 = np.full(n, cap) if budget else None
    book = {k: np.zeros(n) for k in ("pa_lost", "phe_sink", "phe_bound", "phe_lib", "sugar_frag", "sugar_routeB",
                                     "donor_sink", "o2_used", "o2_clipped", "phe_limited")}
    dk = {x: p.donor[x] for x in ARP_SPECIES}
    out = {key: np.zeros((n, len(pot.out_times))) for key in
           ("PA", "PAA", "ARP", "PHE_LIB", "TR_DG1", "TR_TDG", "TR_G", "PH", "O2")}

    def record(j):
        out["PA"][:, j], out["PAA"][:, j], out["ARP"][:, j] = pa, paa, arp
        out["PHE_LIB"][:, j] = book["phe_lib"]
        for x in ("DG1", "TDG", "G"):
            out[f"TR_{x}"][:, j] = TR[x]
        out["PH"][:, j] = buffer_ph(acid / 1000.0, phos)[0] if ph_on else PH_LABEL
        out["O2"][:, j] = (o2 / cap) if budget else o_const

    if 0 in out_idx:
        record(out_idx[0])
    for s in range(1, steps + 1):
        # rate factors at the step start
        fac = fT.copy()
        if ph_on:
            ph, _ = buffer_ph(acid / 1000.0, phos)
            fac = fac * ph_factor(ph, s_ph)
        fo = (arp / (ORDER_REF * 1.0)) ** (p.n_order - 1.0) if p.rival == "R4" else 1.0
        o = (o2 / cap) if budget else o_const
        k_an = k_anox7 * fac * fo
        k_oxr = k_ox7 * (fac / fT) * fTox * fo * o * cu
        e = {"e1": a1 * k_an, "e3": a3 * k_an, "eF": aF * k_an, "eG": gG * k_oxr, "eFox": gF * k_oxr,
             "eB": gB * k_oxr}
        K = sum(e.values())
        dec = arp * (1.0 - np.exp(-K * DT))
        with np.errstate(divide="ignore", invalid="ignore"):
            amt = {k: np.where(K > 0, dec * v / K, 0.0) for k, v in e.items()}
        if budget:   # O2 consumed by the oxidative exits this step may not exceed what is left
            need = p.form.nu_o2 * (amt["eG"] + amt["eFox"] + amt["eB"])
            with np.errstate(divide="ignore", invalid="ignore"):
                scale = np.where(need > o2, o2 / np.where(need > 0, need, 1.0), 1.0)
            for k in ("eG", "eFox", "eB"):
                amt[k] = amt[k] * scale
            dec = sum(amt.values())
        arp = arp - dec
        released = amt["e1"] + amt["e3"] + amt["eG"] + amt["eFox"] + th["r_f"] * amt["eF"]
        book["phe_bound"] += (1.0 - th["r_f"]) * amt["eF"]
        book["phe_lib"] += released
        book["sugar_frag"] += (1.0 - th["m_frag"]) * amt["eF"]
        book["sugar_routeB"] += amt["eB"]
        prod = {"DG1": amt["e1"], "TDG": amt["e3"], "G": amt["eG"], "FOX": amt["eFox"],
                "MGO": th["m_frag"] * amt["eF"]}
        # route B
        b_pa, b_paa = (1.0 - psi) * amt["eB"], psi * amt["eB"]
        # route A on the ARP-made donors (Phe at the step start)
        s_pa, s_paa = np.zeros(n), np.zeros(n)
        phe_use = np.zeros(n)
        o2_events = np.zeros(n)
        acid_src = (amt["eF"] + amt["eFox"] + amt["eB"]) if pot.trap else dec.copy()
        part = {}
        demand = np.zeros(n)
        for x in ARP_SPECIES:
            q = dk[x]
            om = q["omega"] * o
            kstr = q["k"] * voigt if x == "DG1" else q["k"]
            L_str = kstr * (1.0 + om) * phe * fac
            L_sink = q["kappa"] * phe * fac
            L_k0 = q["k0"] * fac
            L = L_str + L_sink + L_k0 + k_tr
            ex = np.exp(-L * DT)
            with np.errstate(divide="ignore", invalid="ignore"):
                grow = np.where(L > 0, (prod[x] / DT) / np.where(L > 0, L, 1.0) * (1.0 - ex), prod[x])
                w = {kk: np.where(L > 0, vv / np.where(L > 0, L, 1.0), 0.0) for kk, vv in
                     (("str", L_str), ("sink", L_sink), ("k0", L_k0), ("tr", k_tr))}
            x_new = X[x] * ex + grow
            lost = X[x] + prod[x] - x_new
            st, sk = lost * w["str"], lost * w["sink"]
            demand += st + th["nu"] * sk
            part[x] = (q, om, x_new, lost, w, st, sk)
        # Phe-limited step (exact booking): the Phe-dependent losses are scaled so that no more Phe is used than
        # the pot holds; what is not consumed stays in the donor pool
        avail = phe + released
        with np.errstate(divide="ignore", invalid="ignore"):
            lim = np.where(demand > avail, avail / np.where(demand > 0, demand, 1.0), 1.0)
        book["phe_limited"] += (lim < 1.0)
        for x in ARP_SPECIES:
            q, om, x_new, lost, w, st, sk = part[x]
            st_l, sk_l = st * lim, sk * lim
            x_new = x_new + (st - st_l) + (sk - sk_l)
            str_al = st_l * (1.0 - q["phi"]) / (1.0 + om)
            s_pa += str_al
            s_paa += st_l - str_al
            o2_events += st_l * om / (1.0 + om)
            sink = sk_l + lost * w["k0"]
            book["donor_sink"] += sink
            phe_use += st_l + th["nu"] * sk_l
            book["phe_sink"] += th["nu"] * sk_l
            TR[x] += lost * w["tr"]
            if pot.trap and x in ("DG1", "TDG", "G"):   # D-f: untrapped osone yields acid
                acid_src += sink
            X[x] = x_new
        phe = np.maximum(phe + released - phe_use, 0.0)
        # PA with loss
        src = s_pa + b_pa
        pa_new = pa * np.exp(-k_loss * DT) + src
        book["pa_lost"] += pa + src - pa_new
        pa = pa_new
        paa = paa + s_paa + b_paa
        if budget:
            used = p.form.nu_o2 * (amt["eG"] + amt["eFox"] + amt["eB"] + o2_events)
            o2_new = o2 - used
            book["o2_clipped"] += np.maximum(-o2_new, 0.0)
            book["o2_used"] += np.minimum(used, o2)
            o2 = np.maximum(o2_new, 0.0)
        acid = acid + y_acid * acid_src
        if s in out_idx:
            record(out_idx[s])
    if ledger:
        out["ledger"] = {"arp0": arp0, "phe0": phe_start, "arp": arp, "phe": phe, "pa": pa, "paa": paa,
                         "donors": sum(X.values()), "trapped": sum(TR.values()), "cap": cap, "o2": o2,
                         **{k: v for k, v in book.items()}}
    out["phe_limited_steps"] = book["phe_limited"]
    return out


def illinois(f, a: np.ndarray, b: np.ndarray, fa: np.ndarray, fb: np.ndarray, iters: int) -> np.ndarray:
    """Vectorised regula falsi with the Illinois modification on brackets [a, b] (fa, fb of opposite sign, or a
    bracket that is not one: then the returned point is meaningless and the caller masks it)."""
    a, b, fa, fb = (np.array(v, dtype=float) for v in (a, b, fa, fb))
    side = np.zeros(a.shape, int)
    x = 0.5 * (a + b)
    for _ in range(iters):
        den = fb - fa
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(den != 0, (a * fb - b * fa) / np.where(den != 0, den, 1.0), 0.5 * (a + b))
        x = np.clip(x, np.minimum(a, b), np.maximum(a, b))
        fx = f(x)
        same_a = np.sign(fx) == np.sign(fa)
        a = np.where(same_a, x, a)
        fa = np.where(same_a, fx, fa)
        fb = np.where(same_a & (side == 1), fb / 2.0, fb)
        b = np.where(same_a, b, x)
        fa = np.where(~same_a & (side == -1), fa / 2.0, fa)
        fb = np.where(same_a, fb, fx)
        side = np.where(same_a, 1, -1)
    return x


def _excess_air(p: Params, sig: np.ndarray, data: Data) -> np.ndarray:
    ar = solve_arp_pot(p, POT_T2["argon"], sig)
    air = solve_arp_pot(p, POT_T2["air"], sig)
    return (air["PA"][:, 0] + air["PAA"][:, 0]) - (ar["PA"][:, 0] + ar["PAA"][:, 0])


def solve_sigma_ox(p: Params, data: Data, grid_pts: int = 9, iters: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """sec. 4.4: sigma_ox in [0, 1) such that the fed ARP's predicted air excess total (Table 2, air minus argon,
    PA + PAA at 120 min) equals the measured one. The first sign change on a grid is refined by the Illinois
    method. No root: infeasible."""
    meas = sum(data.arp_t2["air"]) - sum(data.arp_t2["argon"])
    grid = np.linspace(0.0, 0.999, grid_pts)
    vals = np.stack([_excess_air(p, np.full(p.n, s), data) - meas for s in grid], axis=1)
    j_hit = np.full(p.n, -1)
    for j in range(grid_pts - 1):
        hit = (j_hit < 0) & ((vals[:, j] == 0) | (np.sign(vals[:, j]) * np.sign(vals[:, j + 1]) < 0))
        j_hit = np.where(hit, j, j_hit)
    found = j_hit >= 0
    jj = np.maximum(j_hit, 0)
    a, b = grid[jj], grid[jj + 1]
    exact = vals[np.arange(p.n), jj] == 0
    fa_v = vals[np.arange(p.n), jj]
    fb_v = vals[np.arange(p.n), jj + 1]
    root = illinois(lambda x: _excess_air(p, x, data) - meas, a, b, fa_v, fb_v, iters)
    root = np.where(exact, grid[jj], root)
    return np.where(found, root, np.nan), found


def solve_c_cu5(p: Params, data: Data, sigma_ox: np.ndarray, iters: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """sec. 3.3: c_Cu5 in [1, 25] such that the fed ARP's air + Cu 5 mM total (Table 2) is matched (Illinois in
    log c); outside the bound: infeasible."""
    meas = sum(data.arp_t2["air_cu5"])
    sig = np.where(np.isfinite(sigma_ox), sigma_ox, 0.0)

    def tot(c):
        r = solve_arp_pot(p, POT_T2["air_cu5"], sig, c)
        return r["PA"][:, 0] + r["PAA"][:, 0] - meas

    lo, hi = np.zeros(p.n), np.full(p.n, math.log(C_CU_BOUND[1]))
    f_lo, f_hi = tot(np.exp(lo)), tot(np.exp(hi))
    ok = (f_lo <= 0) & (f_hi >= 0)
    root = illinois(lambda x: tot(np.exp(x)), lo, hi, f_lo, f_hi, iters)
    return np.where(ok, np.exp(root), np.nan), ok


def calibrate(p: Params, data: Data) -> Params:
    """sigma_ox, then c_Cu5; the draw's feasibility (sec. 4.4) and why it failed."""
    sig, ok_s = solve_sigma_ox(p, data)
    cu, ok_c = solve_c_cu5(p, data, sig)
    p.sigma_ox, p.c_cu5 = sig, cu
    route_a_ok = np.all([np.isfinite(p.cal["k"][d]) & (p.cal["k"][d] > 0) & (p.cal["phi"][d] >= 0)
                         for d in DONORS], axis=0)
    p.infeasible_reason = {"no_sigma_ox_root": ~ok_s, "c_cu5_out_of_bound": ok_s & ~ok_c,
                           "route_a_calibration": ~route_a_ok}
    p.feasible = ok_s & ok_c & route_a_ok
    return p


# ============================================================================================== the glucose pot

def glucose_pot(p: Params, data: Data, courses: Optional[Mapping[str, np.ndarray]] = None) -> Dict[str, np.ndarray]:
    """The forced glucose + Phe pot (U5a, U5b, U-ARPglc): route A on the forced donors (EX-049's arithmetic), on 1-DG
    and FOX made from the forced ARP-Phe, and route B on the forced ARP-Phe, at k_d,glc (R3's own; R4's order and pH
    terms) and the reflux O2. Returns PA, PAA (n x GRID) and the in-pot exit rate k_glc,eff(t)."""
    f, th, n = p.form, p.th, p.n
    courses = courses or interpolate(data.forced_inputs, "early_step" if f.early_step else "linear")
    phe = np.full(GRID.size, AA0) if f.cf070 else courses["PHE"]
    arp = (10.0 ** th["cal_arp"])[:, None] * courses["PHEARP"][None, :]
    o2 = th["o2_glc"] if (f.per_pot_o2 and p.rival == "R3") else th["o2_reflux"]
    sig = np.where(np.isfinite(p.sigma_ox), p.sigma_ox, 0.0)
    k_an = ((1.0 - sig) * p.k_d7_glc)[:, None]
    k_ox = (sig * p.k_d7_glc * o2)[:, None]
    fac = np.ones((n, GRID.size))
    acid = np.zeros((n, GRID.size))
    if p.ph_on_glc:
        if f.cf070:
            acid = th["y_acid"][:, None] * cum_trapz((k_an + k_ox) * arp)
        else:
            acid = th["y_acid"][:, None] * (AA0 - phe)[None, :]
        ph, _ = buffer_ph(acid / 1000.0, PHOSPHATE_M)
        fac = fac * ph_factor(ph, th["s_ph"][:, None])
    if p.rival == "R4":
        fac = fac * (np.maximum(arp, 1e-300) / ORDER_REF) ** (p.n_order[:, None] - 1.0)
    k_eff = (k_an + k_ox) * fac
    a1, a3, aF = p.a
    gG, gF, gB = p.g
    # route A on the forced donors: EX-049 verbatim
    cal = p.cal
    src_pa = np.zeros((n, GRID.size))
    src_paa = np.zeros((n, GRID.size))
    for d in DONORS:
        k, phi, om = cal["k"][d], cal["phi"][d], cal["omega"][DONOR_CLASS[d]]
        enc = (10.0 ** th[f"cal_{d.lower()}"])[:, None] * (courses[d] * phe)[None, :]
        src_pa += ((1.0 - phi) * k)[:, None] * enc
        src_paa += (phi * k + om * o2 * k)[:, None] * enc
    out = {"PA_dc": cum_with_loss(src_pa, cal["k_loss"]), "PAA_dc": cum_trapz(src_paa)}
    # 1-DG and FOX from the forced ARP-Phe (exponential steps on the grid)
    prod = {"DG1": a1[:, None] * k_an * fac * arp, "FOX": gF[:, None] * k_ox * fac * arp}
    s_pa_x = np.zeros((n, GRID.size))
    s_paa_x = np.zeros((n, GRID.size))
    voigt = th["voigt_f"] if f.voigt else np.ones(n)
    for x in ("DG1", "FOX"):
        q = p.donor[x]
        kstr = q["k"] * (voigt if x == "DG1" else 1.0)
        om = q["omega"] * o2
        xs = np.zeros(n)
        for s in range(1, GRID.size):
            ph_s = phe[s - 1]
            L = kstr * (1.0 + om) * ph_s + q["kappa"] * ph_s + q["k0"]
            L = L * fac[:, s - 1]
            ex = np.exp(-L * DT)
            pr = prod[x][:, s - 1]
            xs_new = xs * ex + np.where(L > 0, pr / np.where(L > 0, L, 1.0) * (1.0 - ex), pr * DT)
            lost = xs + pr * DT - xs_new
            st = lost * np.where(L > 0, kstr * (1.0 + om) * ph_s * fac[:, s - 1] / np.where(L > 0, L, 1.0), 0.0)
            al = st * (1.0 - q["phi"]) / (1.0 + om)
            s_pa_x[:, s] = al / DT
            s_paa_x[:, s] = (st - al) / DT
            xs = xs_new
    # route B
    rb = gB[:, None] * k_ox * fac * arp
    psi = th.get("psi", np.zeros(n))[:, None]
    out["PA_arp"] = cum_with_loss((1.0 - psi) * rb, cal["k_loss"]) + _cum_rect_with_loss(s_pa_x, cal["k_loss"])
    out["PAA_arp"] = cum_trapz(psi * rb) + np.cumsum(s_paa_x * DT, axis=1)
    out["PA"] = out["PA_dc"] + out["PA_arp"]
    out["PAA"] = out["PAA_dc"] + out["PAA_arp"]
    out["k_eff"] = k_eff
    out["acid"] = acid
    return out


def _cum_rect_with_loss(rate: np.ndarray, k_loss) -> np.ndarray:
    """Amount from a per-step source (rate x DT added at each step) with first-order loss."""
    k = np.asarray(k_loss, dtype=float).reshape(-1, 1)
    out = np.zeros_like(rate)
    acc = np.zeros(rate.shape[0])
    ex = np.exp(-k[:, 0] * DT)
    for s in range(1, rate.shape[1]):
        acc = acc * ex + rate[:, s] * DT
        out[:, s] = acc
    return out


# ================================================================================================ R0-pinned (V4)

def r0_pinned_constants() -> Dict[str, float]:
    """B18's SHIPPED k_go_ak / k_mgo_ak (L/(mmol min), 100 C) and the engine's pyrazine pH factor at label 7.0, read
    from the code (V4)."""
    from src.kinetic_core.parameters_pyrazine import FROZEN_B18, pyrazine_ph_factor
    return {"log10_k_go_ak": float(FROZEN_B18["log10_k_go_ak_100C"]),
            "log10_k_mgo_ak": float(FROZEN_B18["log10_k_mgo_ak_100C"]),
            "ph_factor_7": float(pyrazine_ph_factor(PH_LABEL))}


def r0_glucose_pa(th: Mapping[str, np.ndarray], data: Data, form: Form = PRIMARY) -> np.ndarray:
    """R0-pinned: the glucose pot's PA from route A on GO and MGO only, at B18's shipped constants and pH factor,
    Phe charged as glycine (no ARP exits, no route B, no acid, no PA loss). n x GRID."""
    c = r0_pinned_constants()
    courses = interpolate(data.forced_inputs, "early_step" if form.early_step else "linear")
    phe = np.full(GRID.size, AA0) if form.cf070 else courses["PHE"]
    rate = np.zeros((len(th["cal_go"]), GRID.size))
    for d, key in (("GO", "log10_k_go_ak"), ("MGO", "log10_k_mgo_ak")):
        k = 10.0 ** c[key] * c["ph_factor_7"]
        rate += k * (10.0 ** th[f"cal_{d.lower()}"])[:, None] * (courses[d] * phe)[None, :]
    return cum_trapz(rate)


# ===================================================================================================== evaluate

@dataclass
class Evaluation:
    """One rival's draws under one layer form: the coordinates, feasibility and every unit's predictions."""
    rival: str
    form: Form
    theta: Dict[str, np.ndarray]
    feasible: np.ndarray
    reasons: Dict[str, np.ndarray]
    pred: Dict[str, np.ndarray]          # unit -> (n, values) predicted log10 points (pre-design)
    diag: Dict[str, np.ndarray]
    failed: int = 0


def evaluate(u: np.ndarray, rival: str, form: Form, data: Data, trapped: bool = True) -> Evaluation:
    """Every unit's predicted points for draws ``u`` (n, d) of ``rival`` under ``form``. Points are log10 values on
    each unit's basis BEFORE the design (ratios are formed by the unit's design matrix)."""
    th = from_unit(u, rival, form)
    n = u.shape[0]
    if rival == "R0":
        pa = r0_glucose_pa(th, data, form)
        pred = {"U5a_PA_shape": np.log10(np.maximum(pa[:, [idx(t) for t in FIG2_T]], 1e-300)),
                "U5b_PA_shape": np.log10(np.maximum(pa[:, [idx(t) for t in FIG1_T]], 1e-300))}
        return Evaluation(rival, form, th, np.ones(n, bool), {}, pred, {"R0_PA_120": pa[:, idx(120)]})
    p = calibrate(derive(th, rival, form, data), data)
    sig = np.where(np.isfinite(p.sigma_ox), p.sigma_ox, 0.0)
    cu = np.where(np.isfinite(p.c_cu5), p.c_cu5, 1.0)
    lg = lambda x: np.log10(np.maximum(x, 1e-300))  # noqa: E731
    pred: Dict[str, np.ndarray] = {}
    diag: Dict[str, np.ndarray] = {"sigma_ox": p.sigma_ox, "c_cu5": p.c_cu5}
    failed = 0
    t2 = {atm: solve_arp_pot(p, POT_T2[atm], sig, cu) for atm in ("argon", "air", "air_cu5")}
    f1 = solve_arp_pot(p, POT_FIG1, sig, cu)
    diag["phe_limited_steps_fed"] = sum(r["phe_limited_steps"] for r in list(t2.values()) + [f1])
    # U4: log PA, log PAA for Table 2 (argon, air, air_cu5) then Fig. 1 circles at 7 times
    u4 = [lg(t2[a]["PA"][:, 0]) for a in ("argon", "air", "air_cu5")] + \
         [lg(t2[a]["PAA"][:, 0]) for a in ("argon", "air", "air_cu5")]
    pred["U4"] = np.column_stack(u4 + [lg(f1["PA"]), lg(f1["PAA"])])
    pred["U1"] = lg(f1["PA"])                                   # 10, 30, 60, 120, 240, 300, 500
    diag["t2_total"] = np.column_stack([t2[a]["PA"][:, 0] + t2[a]["PAA"][:, 0] for a in ("argon", "air", "air_cu5")])
    diag["fig1_PA"] = f1["PA"]
    diag["fig1_PH_500"] = f1["PH"][:, -1]
    g = glucose_pot(p, data)
    pred["U5_PA"] = lg(g["PA"][:, [idx(t) for t in GLC_T]])
    pred["U5_PAA"] = lg(g["PAA"][:, [idx(t) for t in GLC_T]])
    diag["glc_ARP_share_PA_120"] = g["PA_arp"][:, idx(120)] / np.maximum(g["PA"][:, idx(120)], 1e-300)
    diag["glc_crossover"] = _crossover(g["PA"], g["PAA"])
    diag["glc_acid_300"] = g["acid"][:, idx(300)]
    # U-ARPglc: integrated exit rate over each Table 3 interval; the isomer form's correction
    ke = g["k_eff"]
    cumk = cum_trapz(ke)
    pred["UARP_c"] = np.column_stack([cumk[:, idx(b)] - cumk[:, idx(a)] for a, b in ARP_INTERVALS])
    diag["k_glc_eff_45"] = ke[:, idx(45)]
    diag["k_eff_1min"] = ke[:, ::2]
    pred["T1"] = np.column_stack([np.log10(p.cal["t3_total"][:, idx(t)] / p.cal["t3_total"][:, idx(500)])
                                  for t in (10.0, 30.0, 60.0, 120.0)])
    pred["T1b"] = np.log10(np.maximum(p.cal["t3_pa"][:, idx(500)], 1e-300) /
                           np.maximum(p.cal["t3_pa"][:, idx(240)], 1e-300))[:, None]
    if trapped:
        tr = {arm: solve_arp_pot(p, POT_TRAP[arm], sig, cu) for arm in U3_ARMS}
        diag["phe_limited_steps_trap"] = sum(r["phe_limited_steps"] for r in tr.values())
        cols = []
        for arm in U3_ARMS:
            a0 = (p.th["arp0"] * TRAP_SCALE)[:, None]
            for sp in U3_SPECIES:
                key = "PHE_LIB" if sp == "PHE" else f"TR_{sp}"
                cols.append(tr[arm][key] / a0)
            diag[f"trap_{arm}_ph_320"] = tr[arm]["PH"][:, -1]
            diag[f"trap_{arm}_arp_left_320"] = tr[arm]["ARP"][:, -1] / a0[:, 0]
        pred["U3_lin"] = np.concatenate(cols, axis=1)            # linear yields, arm x species x time
    # a FAILED run: a feasible draw whose predictions are not finite (counted; the reports refuse on any)
    bad = np.zeros(n, bool)
    for k, v in pred.items():
        if k == "U3_lin":
            bad |= ~np.all(np.isfinite(v), axis=1)
        else:
            bad |= ~np.all(np.isfinite(v) & (v > -299.0), axis=1)
    failed = int(np.sum(bad & p.feasible))
    diag["failed_mask"] = bad & p.feasible
    return Evaluation(rival, form, th, p.feasible, p.infeasible_reason, pred, diag, failed)


def _crossover(pa: np.ndarray, paa: np.ndarray) -> np.ndarray:
    hit = paa[:, 2:] >= pa[:, 2:]
    t = GRID[2:][np.argmax(hit, axis=1)]
    return np.where(hit.any(axis=1), t, np.inf)


FIG1_T = (10.0, 30.0, 60.0, 120.0, 240.0, 300.0, 500.0)
FIG2_T = (10.0, 30.0, 60.0, 120.0, 240.0, 500.0)
GLC_T = (10.0, 30.0, 60.0, 120.0, 240.0, 300.0, 500.0)
ARP_INTERVALS = ((10.0, 30.0), (30.0, 60.0), (60.0, 120.0), (120.0, 300.0))


# ======================================================================================================== units

@dataclass(frozen=True)
class Unit:
    """A unit's data on its value basis: y (dex, after the design), the points' measurement sd (dex), the design
    mapping points -> values, the run offset (tau-like, absolute bases only), the censoring and the role."""
    name: str
    value_basis: str
    y: np.ndarray
    point_sd: np.ndarray
    design: np.ndarray
    run_offset_dex: float = 0.0
    run_groups: Optional[np.ndarray] = None             # per point, the run offset's group (U3: per arm)
    cal_groups: Optional[np.ndarray] = None             # per point, a shared calibration group (U3: per species)
    cal_sd: Optional[np.ndarray] = None
    censored: Optional[np.ndarray] = None
    log_floor: Optional[np.ndarray] = None
    role: str = "claim"
    c: float = 1.0                                      # omega multiplier for adequacy's |dmu| (2 for a ratio)


def _ratio_design(n_points: int, ref: int, keep: Sequence[int]) -> np.ndarray:
    a = np.zeros((len(keep), n_points))
    for i, j in enumerate(keep):
        a[i, j], a[i, ref] = 1.0, -1.0
    return a


def _fig_point_sd(v_err: Tuple[float, float], umol_per_mmol: float) -> float:
    """A keyed figure point: read error (dex), rounding +-0.05 umol/mmol (sd 0.029), s_rep."""
    rnd = ROUND_UMOL_MMOL_SD / max(umol_per_mmol, 1e-12) / math.log(10.0)
    return math.sqrt(v_err[1] ** 2 + rnd ** 2 + S_REP_DEX ** 2)


def build_units(data: Data, u3_synthetic: Optional[np.ndarray] = None) -> Dict[str, Unit]:
    """The phase-1 units on real data (U3 only from a SYNTHETIC vector of linear yields: sec. 7)."""
    lg = math.log10
    units: Dict[str, Unit] = {}
    # U1e / U1l: log PA(t) / PA(120), Fig. 1 circles
    pa = [data.fig1_arp["PHEAL"][t] for t in FIG1_T]
    sd = np.array([_fig_point_sd(x, x[0] * 10.0) for x in pa])
    ylog = np.array([lg(x[0]) for x in pa])
    ref = FIG1_T.index(120.0)
    for name, keep, role in (("U1e", (0, 1, 2), "claim"), ("U1l", (4, 5, 6), "guard")):
        a = _ratio_design(len(FIG1_T), ref, keep)
        units[name] = Unit(name, "ratio", a @ ylog, sd, a, role=role, c=2.0)
    # U4: log PA/PAA, Table 2 (3) + Fig. 1 circles (7); points ordered as evaluate()'s pred["U4"]
    t2 = data.arp_t2
    pts = [lg(t2[a][0]) for a in ("argon", "air", "air_cu5")] + [lg(t2[a][1]) for a in ("argon", "air", "air_cu5")]
    psd = [math.hypot(data.arp_t2_err[a][0], S_REP_DEX) for a in ("argon", "air", "air_cu5")] + \
          [math.hypot(data.arp_t2_err[a][1], S_REP_DEX) for a in ("argon", "air", "air_cu5")]
    pts += [lg(data.fig1_arp["PHEAL"][t][0]) for t in FIG1_T] + [lg(data.fig1_arp["PHEAC"][t][0]) for t in FIG1_T]
    psd += [_fig_point_sd(data.fig1_arp["PHEAL"][t], data.fig1_arp["PHEAL"][t][0] * 10) for t in FIG1_T]
    psd += [_fig_point_sd(data.fig1_arp["PHEAC"][t], data.fig1_arp["PHEAC"][t][0] * 10) for t in FIG1_T]
    a = np.zeros((10, 20))
    for i in range(3):
        a[i, i], a[i, 3 + i] = 1.0, -1.0
    for i in range(7):
        a[3 + i, 6 + i], a[3 + i, 13 + i] = 1.0, -1.0
    units["U4"] = Unit("U4", "ratio", a @ np.array(pts), np.array(psd), a, role="claim", c=2.0)
    # U5a / U5b: absolute log levels PA and PAA, sigma_run as the run's shared offset
    for name, course, times in (("U5a", data.fig2, FIG2_T), ("U5b", data.fig1_glc, FIG1_T)):
        cols = [GLC_T.index(t) for t in times]
        pts = [lg(course["PHEAL"][t][0]) for t in times] + [lg(course["PHEAC"][t][0]) for t in times]
        psd = [_fig_point_sd(course[s][t], course[s][t][0] * 10) for s in ("PHEAL", "PHEAC") for t in times]
        a = np.eye(len(pts))
        units[name] = Unit(name, "absolute", np.array(pts), np.array(psd), a, run_offset_dex=SIGMA_RUN_DEX,
                           role="claim")
        # the PA shape (R0-pinned's comparison)
        ref = times.index(120.0)
        keep = tuple(i for i in range(len(times)) if i != ref)
        ad = _ratio_design(len(times), ref, keep)
        units[f"{name}_PA_shape"] = Unit(f"{name}_PA_shape", "ratio", ad @ np.array(pts[:len(times)]),
                                         np.array(psd[:len(times)]), ad, role="r0_shape", c=2.0)
        del cols
    # G-T3: Table 3 shape (in-sample precondition)
    t3 = {t: data.t3["PHEAL"][t] + data.t3["PHEAC"][t] for t in data.t3["PHEAL"]}
    units["GT3"] = Unit("GT3", "ratio", np.array([lg(t3[t] / t3[500.0]) for t in (10.0, 30.0, 60.0, 120.0)]),
                        np.full(4, 0.05), np.eye(4), role="precondition")
    if u3_synthetic is not None:
        units["U3"] = u3_unit(u3_synthetic)
    return units


def u3_layout() -> List[Tuple[str, str, float]]:
    return [(arm, sp, t) for arm in U3_ARMS for sp in U3_SPECIES for t in U3_TIMES]


def u3_unit(y_lin: np.ndarray, censor: float = U3_CENSOR) -> Unit:
    """U3 from a vector of linear yields (mol / mol ARP charged) in ``u3_layout`` order: absolute log yields,
    sigma_run one offset per arm, a calibration offset per species, left-censoring at the axis resolution."""
    lay = u3_layout()
    y_lin = np.asarray(y_lin, dtype=float)
    cens = y_lin < censor
    y = np.log10(np.where(cens, censor, y_lin))
    sd = np.sqrt(U3_READ_DEX ** 2 + (U3_READ_ABS / np.maximum(y_lin, censor) / math.log(10.0)) ** 2)
    return Unit("U3", "absolute", y, sd, np.eye(len(lay)), run_offset_dex=SIGMA_RUN_DEX,
                run_groups=np.array([U3_ARMS.index(a) for a, _, _ in lay]),
                cal_groups=np.array([U3_SPECIES.index(s) for _, s, _ in lay]),
                cal_sd=np.array([U3_CAL_SD[s] for s in U3_SPECIES]), censored=cens,
                log_floor=np.full(len(lay), math.log10(censor)), role="claim")


def unit_covariance(unit: Unit, form: Form = PRIMARY, measurement_only: bool = False) -> np.ndarray:
    """Sigma_u in dex^2 on the unit's values: A [diag(sd^2 (+ between-points sigma_run^2)) + omega^2 R] A^T, plus the
    run offset (absolute bases only, per run group) and the shared calibration offsets. omega's within-pot
    correlation rho is a form. ``measurement_only`` drops omega (the tie band's redraws, sec. 9.2)."""
    from src.kinetic_core.predictive_densities import MeasurementModelRefusal
    a = unit.design
    m = a.shape[1]
    sd2 = unit.point_sd ** 2
    if form.run_between_points:
        sd2 = sd2 + SIGMA_RUN_DEX ** 2
    cov_p = np.diag(sd2)
    if not measurement_only:
        om = OMEGA_DEX * form.omega_scale
        r = (1.0 - form.rho) * np.eye(m) + form.rho * np.ones((m, m))
        cov_p = cov_p + om ** 2 * r
    if unit.run_offset_dex:
        if unit.value_basis != "absolute":
            raise MeasurementModelRefusal(f"{unit.name}: a run offset on a {unit.value_basis!r} basis cancels")
        grp = unit.run_groups if unit.run_groups is not None else np.zeros(m, int)
        cov_p = cov_p + unit.run_offset_dex ** 2 * (grp[:, None] == grp[None, :])
    if unit.cal_groups is not None:
        sdc = unit.cal_sd.copy()
        if unit.name == "U3":
            sdc[U3_SPECIES.index("DG1")] = form.dg1_cal_sd
        cg = unit.cal_groups
        cov_p = cov_p + (sdc[cg][:, None] * sdc[cg][None, :]) * (cg[:, None] == cg[None, :])
    cov = a @ cov_p @ a.T
    return 0.5 * (cov + cov.T)


def unit_loglik(unit: Unit, f_points: np.ndarray, form: Form = PRIMARY, cov: Optional[np.ndarray] = None
                ) -> np.ndarray:
    """log L_u per draw (nats, dex basis) for predicted points ``f_points`` (n, m). Censored values enter as the
    probability below the floor conditional on the observed ones; with several censored values the conditional
    marginals are multiplied (declared approximation; exact for one)."""
    from scipy.stats import norm
    f = np.atleast_2d(f_points) @ unit.design.T
    cov = unit_covariance(unit, form) if cov is None else cov
    n = f.shape[0]
    bad = ~np.all(np.isfinite(f), axis=1)
    f = np.where(np.isfinite(f), f, 0.0)
    c = unit.censored if unit.censored is not None else np.zeros(len(unit.y), bool)
    c = c.copy()
    floor = unit.log_floor + math.log10(form.censor_scale) if unit.log_floor is not None else None
    o = ~c
    out = np.zeros(n)
    if o.any():
        s_oo = cov[np.ix_(o, o)]
        ch = np.linalg.cholesky(s_oo)
        r = (unit.y[o][None, :] - f[:, o])
        z = np.linalg.solve(ch, r.T)
        out += -0.5 * np.sum(z ** 2, axis=0) - np.sum(np.log(np.diag(ch))) - 0.5 * o.sum() * math.log(2 * math.pi)
    if c.any():
        if o.any():
            s_co = cov[np.ix_(c, o)]
            k = np.linalg.solve(cov[np.ix_(o, o)], s_co.T).T
            mu = f[:, c] + (unit.y[o][None, :] - f[:, o]) @ k.T
            s_cc = cov[np.ix_(c, c)] - k @ s_co.T
        else:
            mu, s_cc = f[:, c], cov[np.ix_(c, c)]
        sdv = np.sqrt(np.diag(s_cc))
        out += np.sum(norm.logcdf((floor[c][None, :] - mu) / sdv[None, :]), axis=1)
    out[bad] = -np.inf
    return out


def arp_interval_measured(data: Data, f_iso: Optional[np.ndarray] = None) -> np.ndarray:
    """Table 3's ln(ARP(t1) / ARP(t2)) per interval (cal_ARP cancels); with the isomers form the measured value is
    ARP + I, I inert, formed from ARP (proportional to the integral of the forced course) and reaching f_iso of the
    measured value at 300 min. Returns (n or 1, intervals)."""
    ts = sorted(data.arp_glc)
    m = np.array([data.arp_glc[t] for t in ts])
    if f_iso is None:
        true = m[None, :]
    else:
        cur = interpolate(data.forced_inputs)["PHEARP"]
        ci = cum_trapz(cur)
        frac = np.array([ci[idx(t)] / ci[idx(300.0)] for t in ts])
        iso = np.asarray(f_iso)[:, None] * m[ts.index(300.0)] * frac[None, :]
        true = np.maximum(m[None, :] - iso, 1e-9)
    return np.column_stack([np.log(true[:, ts.index(a)] / true[:, ts.index(b)]) for a, b in ARP_INTERVALS])


def uarp_loglik(data: Data, c_int: np.ndarray, cv: np.ndarray, form: Form = PRIMARY,
                f_iso: Optional[np.ndarray] = None, measured: Optional[np.ndarray] = None) -> np.ndarray:
    """U-ARPglc (sec. 6): formation >= 0 implies ln(ARP1 / ARP2) <= integral of k_glc,eff over the interval. The
    observed ln ratio carries rounding (uniform +-0.5 umol on each count), the analyser CV (a coordinate on [2, 10] %)
    and, as a form, sigma_run between points; the likelihood of a rival is sum log Phi((c - y) / s)."""
    from scipy.stats import norm
    ts = sorted(data.arp_glc)
    umol = np.array([data.arp_glc[t] * 10.0 for t in ts])                     # umol in the pot
    y = arp_interval_measured(data, f_iso) if measured is None else measured
    cv = np.asarray(cv, dtype=float)[:, None]
    rnd = (0.5 / math.sqrt(3.0)) / umol                                       # ln-scale rounding sd per point
    var_pt = cv ** 2 + rnd[None, :] ** 2
    if form.run_between_points:
        var_pt = var_pt + (SIGMA_RUN_DEX * math.log(10.0)) ** 2
    s = np.column_stack([np.sqrt(var_pt[:, ts.index(a)] + var_pt[:, ts.index(b)]) for a, b in ARP_INTERVALS])
    return np.sum(norm.logcdf((c_int - y) / s), axis=1)


# ========================================================================================================= nulls

def nulls(data: Data) -> Dict[str, Any]:
    """sec. 4.4's chemistry-free nulls (points on each unit's pre-design layout)."""
    lg = math.log10
    t4 = data.fed
    mean_ratio = {atm: float(np.mean([t4[d][atm][0] / t4[d][atm][1] for d in DONORS])) for atm in ("argon", "air")}
    mean_ratio["air_cu5"] = data.fed_cu_mgo[0] / data.fed_cu_mgo[1]
    # U4 null: a constant PA/PAA per atmosphere (Fig. 1's reflux pot as air), placed on the PA points with PAA = 0
    u4 = np.zeros(20)
    for i, atm in enumerate(("argon", "air", "air_cu5")):
        u4[i] = lg(mean_ratio[atm])
    u4[6:13] = lg(mean_ratio["air"])
    # shape null: the fed MGO Table 3 PA shape at Fig. 1's times
    t3 = data.t3["PHEAL"]
    u1 = np.array([lg(t3[t]) for t in FIG1_T])
    # level null: the fed donors' mean yield per unit of donor exposure x the integrated forcing
    courses = interpolate(data.forced_inputs)
    dose = cum_trapz(sum(courses[d] for d in DONORS))
    y_pa = float(np.mean([t4[d]["air"][0] for d in DONORS])) / (D0 * 60.0)
    y_paa = float(np.mean([t4[d]["air"][1] for d in DONORS])) / (D0 * 60.0)
    u5 = {"PA": np.array([lg(y_pa * dose[idx(t)]) for t in GLC_T]),
          "PAA": np.array([lg(y_paa * dose[idx(t)]) for t in GLC_T])}
    return {"U4": u4, "U1": u1, "U5": u5, "mean_ratio": mean_ratio}


def u3_null_points(y_lin: np.ndarray) -> np.ndarray:
    """U3's atmosphere-free null: each species' yield equal in both arms (the arms' geometric mean per time); Phe
    liberated linear in time through the arm mean (least squares through the origin on the linear values)."""
    lay = u3_layout()
    y = np.maximum(np.asarray(y_lin, float), U3_CENSOR)
    out = np.zeros(len(lay))
    for sp in U3_SPECIES:
        for t in U3_TIMES:
            i = [k for k, (a, s, tt) in enumerate(lay) if s == sp and tt == t]
            out[i] = np.mean(np.log10(y[i]))
        if sp == "PHE":
            idxs = [k for k, (a, s, tt) in enumerate(lay) if s == sp]
            ts = np.array([lay[k][2] for k in idxs])
            mean_t = {t: 10.0 ** np.mean(np.log10(y[[k for k in idxs if lay[k][2] == t]])) for t in U3_TIMES}
            slope = sum(mean_t[t] * t for t in U3_TIMES) / sum(t * t for t in U3_TIMES)
            out[idxs] = np.log10(slope * ts)
    return out


def points_for(unit: str, ev: Evaluation) -> np.ndarray:
    """A unit's predicted POINTS (n, m) on its pre-design layout, from an Evaluation."""
    pr = ev.pred
    if unit in ("U1e", "U1l"):
        return pr["U1"]
    if unit == "U4":
        return pr["U4"]
    if unit == "GT3":
        return pr["T1"]
    if unit == "U3":
        return np.log10(np.maximum(pr["U3_lin"], 1e-300))
    if unit in ("U5a", "U5b"):
        times = FIG2_T if unit == "U5a" else FIG1_T
        cols = [GLC_T.index(t) for t in times]
        return np.concatenate([pr["U5_PA"][:, cols], pr["U5_PAA"][:, cols]], axis=1)
    if unit in ("U5a_PA_shape", "U5b_PA_shape"):
        if unit in pr:
            return pr[unit]
        times = FIG2_T if unit.startswith("U5a") else FIG1_T
        return pr["U5_PA"][:, [GLC_T.index(t) for t in times]]
    raise KeyError(unit)


def null_points(unit: str, data: Data, u3_y_lin: Optional[np.ndarray] = None) -> np.ndarray:
    """The null's points (1, m) on a unit's layout (sec. 4.4)."""
    nl = nulls(data)
    if unit in ("U1e", "U1l"):
        return nl["U1"][None, :]
    if unit == "U4":
        return nl["U4"][None, :]
    if unit in ("U5a", "U5b"):
        times = FIG2_T if unit == "U5a" else FIG1_T
        cols = [GLC_T.index(t) for t in times]
        return np.concatenate([nl["U5"]["PA"][cols], nl["U5"]["PAA"][cols]])[None, :]
    if unit in ("U5a_PA_shape", "U5b_PA_shape"):
        times = FIG2_T if unit.startswith("U5a") else FIG1_T
        t3 = data.t3["PHEAL"]
        return np.array([[math.log10(t3[t]) for t in times]])
    if unit == "U3":
        return u3_null_points(u3_y_lin)[None, :]
    raise KeyError(unit)

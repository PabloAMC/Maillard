"""
src/m6_reduced_power.py -- case M6 (pyrazines): the falsification arithmetic and the reduced fed-pot mini-model whose
SYNTHETIC data size the designed decider under the adopted sec. 7 item 5 / M1.1 judge.

The case is docs/reference/cases/M6_pyrazines.md (secs. 4 - 5). The generator is
scripts/generators/m6_reduced_power_check.py. Nothing here runs the engine, fits an engine constant, reads a decisive,
reserved or TEST value, or writes anywhere. Printed values used by job (a) are the spent / seen numbers the case quotes
from their dossiers (cited at each constant).

UNITS: mol/L (M) and minutes throughout; log values in dex (log10).

THE MINI-MODEL (one fed pot: amino acid A + methylglyoxal M, pH held, T held; M is the MEASURED forcing of the design,
held at its charge in the reduced model: a declared limit)
- Strecker flux F = kS [A][M] (one event = one aminoketone AK + one Strecker aldehyde + CO2).
- Part C (the aminoketone fate):
  * C0 (R0 / R2: incumbent, condensation instant): AK -> 1/2 DMP at once; no AK pool; DMP rate F / 2.
  * C1 (R1, the claim, as it survives job (a)): a free AK pool, d AK/dt = F - k2 AK^2 - ks AK, ks = kp [M] (the
    pyrrole exit, pseudo-first-order in MGO); DMP rate k2 / 2 AK^2; pyrrole rate ks AK. k2 at Elfsson 1999's class.
  * C3 (R3, the ring bypasses free AK): DMP from two dicarbonyl - amino-acid imines, rate kI ([A][M])^2 (imine in its
    dilute pre-equilibrium); no free AK.
- Part S (the amine budget): amine loss rate (kS + kad) [A][M]; S0 (R0, R3 as written): kad = 0; S1 (R1 / R2, the
  claim's pH-flat, buffer-flat sink): kad / kS in a declared band; S5 (MF-25 (b), CF-095 explanation 2: the sink is a
  phosphate / mixture effect): kad = 0 without phosphate, the S1 band with it.

THE DESIGNED UNITS (all within one pot, so a run-level term cancels: sec. 7 item 5 lets within-run values credit)
- "onset": log10 DMP(10 min) / DMP(40 min) (one pot; a lag in C1 only).
- "mgo_step": log10 of the DMP rate over 45 - 60 min / 20 - 35 min after MGO is raised x3 at 40 min in the same pot.
- "dilution_step": the same ratio after the pot is diluted x4 at 40 min (C3: x 1/256; C0 and C1: x 1/16).
- "ak_trap": log10 (trapped AK at 40 min x the trap's recovery + the trap blank) / the trap's reporting limit.
- "budget_plain", "budget_phosphate": log10 (Strecker events / amine consumed) at 60 min in an unbuffered (pH-stat)
  pot and in a 0.4 M phosphate pot.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from src.kinetic_core import predictive_densities as pdn

L10 = math.log10
R_GAS = 8.314462618

# ---------------------------------------------------------------------------------------------------------------------
# Printed values (spent / seen; dossiers named)
# ---------------------------------------------------------------------------------------------------------------------

#: Zhou 2024 Table 2 (zhou2024 dossier): zero-order rates (uM/min) and R2 of c = kt + b over t = 0, 30, 60, 90, 120 min,
#: Ala + GO -> pyrazine and Ala + MGO -> 2,5-DMP, 20 + 20 mM, pH 8 (initial), 100 / 110 / 120 C. B18 fit rows: SPENT.
ZHOU_TIMES = np.array([0.0, 30.0, 60.0, 90.0, 120.0])
ZHOU_UNITS = {
    "PZ_100C": {"rate_uM_min": 0.0279, "r2": 0.9789, "T_C": 100.0},
    "PZ_110C": {"rate_uM_min": 0.0791, "r2": 0.9918, "T_C": 110.0},
    "PZ_120C": {"rate_uM_min": 0.1507, "r2": 0.9634, "T_C": 120.0},
    "DMP_100C": {"rate_uM_min": 0.0035, "r2": 0.9891, "T_C": 100.0},
    "DMP_110C": {"rate_uM_min": 0.0100, "r2": 0.9946, "T_C": 110.0},
    "DMP_120C": {"rate_uM_min": 0.0230, "r2": 0.9853, "T_C": 120.0},
}
ZHOU_CONC_M = 0.020

#: Elfsson 1999 (elfsson1999 dossier): 5-aminolevulinic acid, -d[ALA]/dt = k2 [ALA]^2; 50 C: k2 5.53 M^-1 h^-1 at
#: pH 7.42 (SE 0.58), 0.069 at pH 4.81; none at pH 2.35; Ea 43.7 kJ/mol (slope CI -> 39.4 - 48.0) at pH 7.53, 37 - 85 C.
ELFSSON_K2_50C_PH742 = 5.53 / 60.0          # M^-1 min^-1
ELFSSON_PH_SLOPE = L10(5.53 / 0.069) / (7.42 - 4.81)   # 0.729 dex per pH unit
ELFSSON_EA = (39.4e3, 43.7e3, 48.0e3)      # J/mol (low, printed, high)

#: Hofmann 2000 Table 3 (hofmann2000 dossier): MGO + Phe 0.1 + 0.1 M, 0.5 M phosphate pH 7.0, 98 C, air:
#: (min, PA umol, PAA umol, Phe degraded %) of 1000 umol Phe. DEV; the Phe column carries an integrity flag
#: (it coincides with Hofmann 2000b's glucose / Phe column).
HOFMANN_T3 = ((10, 6.2, 10.1, 9), (30, 9.5, 17.2, 18), (60, 10.6, 20.6, 36), (120, 12.3, 23.1, 49),
              (240, 13.5, 24.8, 58), (300, 13.3, 25.8, 61), (500, 12.2, 27.6, 69))

#: Piloty 1979 Table 4 (piloty1979 dossier): MGO 312 mM + 13 amino acids 16 mM each, 0.4 M phosphate, 100 C, 3 h:
#: Ala loss % at pH 4.52 / 6.97 / 9.50 (amino-acid analyser; irreversible loss). SPENT (CF-095, EX-058).
PILOTY_ALA_LOSS = (0.32, 0.27, 0.35)
PILOTY_MGO_M, PILOTY_T_MIN = 0.312, 180.0

#: H0 (EX-059) terms, the only measured within- and between-laboratory terms in the programme (thiol lane). The pyrazine
#: lane has none (src/evidence_tiers.py: "this record does not yet carry" a pooled reaction-class spread). Used here
#: ONLY as "the largest noise any listed source supports" for a falsification (sec. 7 item 5), never to credit.
OMEGA_DEX = pdn.OMEGA_DEX
TAU_ABS_DEX = pdn.TAU_ABS_DEX

#: B18's shipped MGO Strecker constant (CF-095: 3e-8 L/(mmol min) = 3e-5 M^-1 min^-1) and EX-021's gap to Hofmann's fed
#: MGO + Phe rate (~2.3 dex above it).
B18_KS_MGO = 3.0e-5
EX021_GAP_DEX = 2.3


def elfsson_k2(T_C: float, pH: float, ea: float = ELFSSON_EA[1]) -> float:
    """Elfsson's ALA consumption constant (M^-1 min^-1) carried to T and pH by its own Ea and pH slope (an analogue)."""
    T = T_C + 273.15
    k50 = ELFSSON_K2_50C_PH742 * math.exp(-ea / R_GAS * (1.0 / T - 1.0 / 323.15))
    return k50 * 10.0 ** (ELFSSON_PH_SLOPE * (pH - 7.42))


# ---------------------------------------------------------------------------------------------------------------------
# Job (a): the falsification arithmetic
# ---------------------------------------------------------------------------------------------------------------------

def lag_shape(t: np.ndarray, kind: str, x: float) -> np.ndarray:
    """Normalised product curve P(t) (any scale) for a free-AK pool fed at a constant flux.

    kind "exit": AK relaxes at the first-order exit rate x (1/min) and the condensation is minor: P ~ int (1 - e^{-x s})^2.
    kind "cond": no exit, condensation only, relaxation time x (min): AK = A_ss tanh(s / (2 x)) ... P = int tanh^2."""
    t = np.asarray(t, float)
    if kind == "exit":
        k = x
        return t - 2.0 * (1.0 - np.exp(-k * t)) / k + (1.0 - np.exp(-2.0 * k * t)) / (2.0 * k)
    if kind == "cond":
        tau = x
        # d a/ds = 1 - a^2 in units of the steady state and tau' = 2 tau: a = tanh(s / tau'); int tanh^2 = s - tau' tanh
        tp = 2.0 * tau
        return t - tp * np.tanh(t / tp)
    if kind == "linear":
        return t.copy()
    raise ValueError(kind)


def r2_linear(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """R2 of the least-squares line with intercept through (t, y); y may be (R x n)."""
    y = np.atleast_2d(y)
    tc = t - t.mean()
    yc = y - y.mean(axis=1, keepdims=True)
    sxy = yc @ tc
    syy = np.sum(yc ** 2, axis=1)
    return (sxy ** 2) / (np.sum(tc ** 2) * np.where(syy > 0, syy, np.nan))


def p_r2_at_least(shape: np.ndarray, r2_obs: float, rng: np.random.Generator, reps: int = 20000,
                  sigmas: Sequence[float] = tuple(np.geomspace(0.003, 1.5, 40))) -> Dict[str, float]:
    """max over the noise SD (additive, as a fraction of the final value; and proportional) of P(R2 >= r2_obs) when the
    truth is ``shape`` at ZHOU_TIMES: "the largest noise any listed source supports" taken to its least favourable
    value for the falsification. Returns the max p and the sigma that attains it."""
    s = shape / shape[-1]
    best = (0.0, None, None)
    z = rng.standard_normal((reps, s.size))
    for form in ("additive", "proportional"):
        for sg in sigmas:
            y = s[None, :] + sg * z * (1.0 if form == "additive" else s[None, :])
            p = float(np.mean(r2_linear(ZHOU_TIMES, y) >= r2_obs))
            if p > best[0]:
                best = (p, float(sg), form)
    return {"p_max": best[0], "sigma_at_max": best[1], "noise_form": best[2],
            "r2_noiseless": float(r2_linear(ZHOU_TIMES, s[None, :])[0])}


def onset_falsification(rng: np.random.Generator, reps: int = 20000) -> Dict[str, Any]:
    """FX-C1: Zhou 2024's six zero-order R2 against a free-AK pool.

    (i) R1 with Elfsson's class and NO first-order AK exit: relaxation tau = 1 / (2 sqrt(F k2)) with F = 2 x the DMP
    (or pyrazine) rate (yield 1 at steady state); every corner of Elfsson's band (Ea 39.4 - 48.0, the pH slope to pH 8).
    (ii) the smallest first-order exit rate ks a free-AK pool needs for the six R2 (Holm), and what it implies."""
    out: Dict[str, Any] = {"units": {}, "no_exit": {}, "exit_grid": {}}
    k2_corners = {f"Ea{int(ea/1e3)}": elfsson_k2(100.0, 8.0, ea) for ea in ELFSSON_EA}
    for u, d in ZHOU_UNITS.items():
        F = 2.0 * d["rate_uM_min"] * 1e-6
        k2s = [elfsson_k2(d["T_C"], 8.0, ea) for ea in ELFSSON_EA]
        k2 = max(k2s)                                   # the fastest corner is the most favourable for R1
        tau = 1.0 / (2.0 * math.sqrt(F * k2))
        res = p_r2_at_least(lag_shape(ZHOU_TIMES, "cond", tau), d["r2"], rng, reps)
        # the k2 needed for a lag-free onset: tau <= 10 min
        k2_need = 1.0 / (4.0 * 10.0 ** 2 * F)
        out["no_exit"][u] = {"F_M_min": F, "k2_fastest_corner": k2, "tau_min": tau, "tau_h": tau / 60.0, **res,
                             "k2_needed_for_tau_10min": k2_need, "gap_dex": L10(k2_need / k2)}
    p_no_exit = {u: v["p_max"] for u, v in out["no_exit"].items()}
    out["no_exit_holm"] = pdn.holm(p_no_exit)
    out["no_exit_falsified"] = any(v["p_holm"] < 0.05 for v in out["no_exit_holm"].values())
    out["k2_corners_100C_pH8"] = k2_corners
    # (ii) exit grid
    grid = np.geomspace(0.005, 2.0, 60)
    ks_min = None
    for ks in grid:
        ps = {u: p_r2_at_least(lag_shape(ZHOU_TIMES, "exit", ks), d["r2"], rng, reps // 4)["p_max"]
              for u, d in ZHOU_UNITS.items()}
        h = pdn.holm(ps)
        ok = all(v["p_holm"] >= 0.05 for v in h.values())
        out["exit_grid"][f"{ks:.4g}"] = {"p": ps, "survives": ok}
        if ok and ks_min is None:
            ks_min = float(ks)
    out["ks_min_per_min"] = ks_min
    if ks_min is not None:
        rate = ZHOU_UNITS["DMP_100C"]["rate_uM_min"] * 1e-6
        k2 = elfsson_k2(100.0, 8.0, ELFSSON_EA[2])
        ak = math.sqrt(2.0 * rate / k2)                 # DMP rate = k2 / 2 AK^2
        F = ks_min * ak                                 # exit-dominated: F ~ ks AK
        kS = F / ZHOU_CONC_M ** 2
        out["implied_at_ks_min"] = {
            "AK_steady_uM": ak * 1e6, "strecker_flux_uM_min": F * 1e6,
            "kS_M_min": kS, "kS_over_B18_dex": L10(kS / B18_KS_MGO),
            "kS_over_hofmann_implied_dex": L10(kS / (B18_KS_MGO * 10 ** EX021_GAP_DEX)),
            "yield_DMP_per_strecker_event": 2.0 * rate / F,
            "pyrrole_uM_at_120min": F * 120.0 * 1e6,
            "note": "DMP at Elfsson's fastest corner at 100 C, pH 8; a slower k2 needs a larger flux"}
    return out


def strecker_only_amine_exit() -> Dict[str, Any]:
    """FX-S1: the incumbent's 'the Strecker step is the only amine exit' (consumed = Strecker events; constant-free).

    (i) within one pot, Hofmann 2000 Table 3: log10 (PA + PAA) / Phe consumed; (ii) across laboratories, Piloty's Ala loss
    against the Strecker constant Zhou's DMP rate implies under R0's yield 1. Each with the break-even noise: the SD
    at which the one-sided p reaches 0.05, against the largest listed terms (omega per value, tau_ABS per laboratory)."""
    from scipy.stats import norm

    rows = []
    for t, pa, paa, deg in HOFMANN_T3:
        ratio = (pa + paa) / (10.0 * deg)
        rows.append({"t_min": t, "products_over_consumed": ratio, "log10": L10(ratio)})
    d10 = -rows[0]["log10"]
    sd_omega = math.sqrt(2.0) * OMEGA_DEX            # a ratio of two assays, omega on each value
    hof = {"rows": rows, "deficit_10min_dex": d10, "sd_two_values_omega": sd_omega,
           "p_one_sided_at_omega": float(1.0 - norm.cdf(d10 / sd_omega)),
           "break_even_sd_dex": d10 / norm.ppf(0.95),
           "integrity_flag": "the Phe column coincides with Hofmann 2000b's glucose / Phe column (dossier caveat)"}
    # (ii) Piloty vs Zhou
    rate = ZHOU_UNITS["DMP_100C"]["rate_uM_min"] * 1e-6
    kS = 2.0 * rate / ZHOU_CONC_M ** 2               # R0: two Strecker events per DMP
    exposure_max = kS * PILOTY_MGO_M * PILOTY_T_MIN  # no MGO loss: an upper bound on R0's Ala loss
    out_p = {}
    for pH, x in zip((4.52, 6.97, 9.50), PILOTY_ALA_LOSS):
        need = -math.log(1.0 - x)
        gap = L10(need / exposure_max)
        out_p[f"pH{pH}"] = {"loss": x, "required_exposure": need, "r0_max_exposure": exposure_max, "gap_dex": gap}
    gap7 = out_p["pH6.97"]["gap_dex"]
    buffer_max_dex = L10(1.0 + 40.0 * 0.4)            # M5's kappa_P upper edge 40 /M at 0.4 M (a sensitivity)
    sd_cross = math.sqrt(2.0 * TAU_ABS_DEX ** 2 + 2.0 * OMEGA_DEX ** 2)
    pil = {"kS_R0_from_zhou_M_min": kS, "rows": out_p, "sd_cross_lab_tau_abs": sd_cross,
           "buffer_factor_max_dex": buffer_max_dex,
           "p_one_sided_no_buffer": float(1.0 - norm.cdf(gap7 / sd_cross)),
           "p_one_sided_buffer_max": float(1.0 - norm.cdf((gap7 - buffer_max_dex) / sd_cross)),
           "break_even_between_lab_sd_no_buffer": math.sqrt(max((gap7 / norm.ppf(0.95)) ** 2 / 2.0 - OMEGA_DEX ** 2, 0)),
           "break_even_between_lab_sd_buffer_max":
               math.sqrt(max(((gap7 - buffer_max_dex) / norm.ppf(0.95)) ** 2 / 2.0 - OMEGA_DEX ** 2, 0))}
    p = {"hofmann_within_pot": hof["p_one_sided_at_omega"], "piloty_vs_zhou": pil["p_one_sided_no_buffer"]}
    return {"hofmann": hof, "piloty_vs_zhou": pil, "holm": pdn.holm(p),
            "falsified": any(v["p_holm"] < 0.05 for v in pdn.holm(p).values())}


#: Kanzler 2026 (kanzler2026 dossier, Figs. 1F, 2B; first reader; SEEN by CK-46 and the M4 sweep): MGO 100 mM +
#: Ala 100 mM, unbuffered pH 5.0 (initial), 100 C, triplicate pots: MGO remaining (% of C0) at 0 / 30 / 60 / 120 min
#: 100 / 31.8 / 14.0 / 4.3; Ala remaining at 120 min 100 % (SD 1, read +- 3); MGO alone second order, k2 1.85e-4
#: mM^-1 min^-1.
KANZLER_T = np.array([0.0, 30.0, 60.0, 120.0])
KANZLER_MGO_ALA = np.array([100.0, 31.8, 14.0, 4.3])
KANZLER_ALA_LEFT, KANZLER_ALA_SD, KANZLER_READ = 100.0, 1.0, 3.0
KANZLER_K2_ALONE = 1.85e-4


def kanzler_one_to_one() -> Dict[str, Any]:
    """FX-S2: 'the amine-dependent MGO loss is a 1:1 amine-consuming adduct' (the form that would let one sink carry
    both Piloty's amine loss and Kanzler's x2.6 - 60 MGO acceleration). Constant-free: the MGO lost beyond the most the
    MGO-alone (second-order) path can take on the with-Ala MGO course must take one Ala each. Within one study, run SD
    of triplicate pots (n - 1 = 2 df), the figure read error as a floor."""
    from scipy.stats import t as student_t

    m = KANZLER_MGO_ALA
    integ_m2 = float(np.sum(0.5 * (m[1:] ** 2 + m[:-1] ** 2) * np.diff(KANZLER_T)))
    alone_max = KANZLER_K2_ALONE * integ_m2
    amine_dep = (m[0] - m[-1]) - alone_max
    pred_left = 100.0 - amine_dep
    sd = max(KANZLER_ALA_SD, KANZLER_READ)
    tstat = (KANZLER_ALA_LEFT - pred_left) / sd
    p = float(2.0 * student_t.sf(abs(tstat), 2))
    return {"mgo_lost_mM": float(m[0] - m[-1]), "alone_path_max_mM": alone_max, "amine_dependent_min_mM": amine_dep,
            "ala_left_predicted_max_pct": pred_left, "ala_left_observed_pct": KANZLER_ALA_LEFT, "sd_used": sd,
            "t_2df": tstat, "p_two_sided": p, "falsified": p < 0.05,
            "consistent_check_piloty_sized_sink": _piloty_sink_at_kanzler()}


def _piloty_sink_at_kanzler() -> Dict[str, float]:
    """The Piloty-sized pH-flat sink (k >= -ln(1 - 0.27) / (312 mM x 180 min)) on Kanzler's with-Ala MGO course."""
    k = -math.log(1.0 - PILOTY_ALA_LOSS[1]) / (PILOTY_MGO_M * 1e3 * PILOTY_T_MIN)       # mM^-1 min^-1, a lower bound
    m = KANZLER_MGO_ALA
    integ = float(np.sum(0.5 * (m[1:] + m[:-1]) * np.diff(KANZLER_T)))
    return {"k_lower_bound_mM_min": k, "mgo_exposure_mM_min": integ,
            "ala_consumed_pct_at_lower_bound": 100.0 * (1.0 - math.exp(-k * integ))}


def adams_ph_budget() -> Dict[str, Any]:
    """D2 (Adams 2008 Table 3, 130 C, 30 min, 1:1; SPENT by EX-052): what R1's single free-AK pool needs to hold the
    pH 4 (acetate) vs pH 7 (phosphate) contrast under Elfsson's pH law (k2 x 1/80 - 1/154 at pH 4)."""
    dmp7, pyr7, dmp4, pyr4 = 0.02, 0.020, 0.042, 1.6
    F7, F4 = 2 * dmp7 + pyr7, 2 * dmp4 + pyr4
    lo, hi = 10 ** (ELFSSON_PH_SLOPE * (7.42 - 4.81)), 10 ** (ELFSSON_PH_SLOPE * 3.0)
    ak_ratio = (math.sqrt(dmp4 / dmp7 * lo), math.sqrt(dmp4 / dmp7 * hi))
    return {"AK_flux_ratio_pH4_over_pH7": F4 / F7, "AK_pool_ratio": ak_ratio,
            "kp_ratio_needed": ((pyr4 / pyr7) / ak_ratio[1], (pyr4 / pyr7) / ak_ratio[0]),
            "reading": "no contradiction: R1 needs the Strecker flux x28 at pH 4 (acetate) over pH 7 (phosphate), the "
                       "buffer confounded with pH; a prediction, not a falsification"}


# ---------------------------------------------------------------------------------------------------------------------
# The mini-model and the designed units
# ---------------------------------------------------------------------------------------------------------------------

UNITS_C = ("onset", "mgo_step", "dilution_step", "ak_trap")
UNITS_S = ("budget_plain", "budget_phosphate")
ONSET_T = np.array([10.0, 20.0])
POT_C = 0.020                     # M, Ala and MGO at 1:1
POT_PH = 7.0
AK_REPORTING_LIMIT = 0.5e-6       # M: the trap's reporting limit (a design requirement, declared)
HEATUP_MIN = 1.0                  # the vessel's measured heat-up (a design requirement): all rivals alike


@dataclass(frozen=True)
class Bands:
    """Declared prior bands of the reduced model (log10 unless stated); every band end is reported."""
    log_kS_C0: Tuple[float, float] = (L10(B18_KS_MGO), L10(B18_KS_MGO) + EX021_GAP_DEX + 0.6)
    #: R1's Strecker constant: job (a) puts it within ~0.3 dex of Hofmann 2000's implied rate (6e-3 M^-1 min^-1);
    #: declared +- 0.8 dex around it
    log_kS_C1: Tuple[float, float] = (L10(B18_KS_MGO) + EX021_GAP_DEX - 0.8, L10(B18_KS_MGO) + EX021_GAP_DEX + 0.8)
    #: R1's first-order AK exit at the pot's MGO (1/min): lower end = job (a)'s ks_min (0.023 /min at Zhou's 20 mM
    #: MGO; the designed pot has the same MGO); upper end 1 /min, a declared, unsourced width (no kp is measured)
    log_ks: Tuple[float, float] = (L10(0.023), L10(1.0))
    #: the ALA -> aminoacetone transfer of Elfsson's k2 (dex): held at 0 in the primary run (in the exit-dominated
    #: regime job (a) leaves, k2 moves no designed unit at first order); -1 / +1 dex are the declared sensitivity
    log_k2_dex_around_elfsson: Tuple[float, float] = (0.0, 0.0)
    log_kI: Tuple[float, float] = (-1.0, 3.0)              # M^-3 min^-1 (R3): scale only
    log_kad_over_kS: Tuple[float, float] = (0.0, 1.5)      # S1 / S5: 1 - 30x (Hofmann 0.06 - 0.18, Piloty)


def _rk4_ak(F: np.ndarray, k2: np.ndarray, ks: np.ndarray, t_end: float, dt: float, ak0: Optional[np.ndarray] = None
            ) -> Tuple[np.ndarray, np.ndarray]:
    """Integrate d AK/dt = F - k2 AK^2 - ks AK and d P/dt = k2/2 AK^2 from 0 to t_end (vectorised over draws).
    Returns (AK(t_end), P(t_end))."""
    ak = np.zeros_like(F) if ak0 is None else ak0.copy()
    P = np.zeros_like(F)
    t_end = np.asarray(t_end, float)
    n = int(math.ceil(float(np.max(t_end)) / dt))
    h = t_end / n

    def f(a):
        return F - k2 * a * a - ks * a, 0.5 * k2 * a * a
    for _ in range(n):
        a1, p1 = f(ak)
        a2, p2 = f(ak + 0.5 * h * a1)
        a3, p3 = f(ak + 0.5 * h * a2)
        a4, p4 = f(ak + h * a3)
        ak = np.maximum(ak + h / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4), 0.0)
        P = P + h / 6.0 * (p1 + 2 * p2 + 2 * p3 + p4)
    return ak, P


def _defensive_uniform(lo: float, hi: float, n: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Draws for a coordinate whose prior is uniform on [lo, hi], from the defensive proposal q = 0.9 U + 0.05
    Beta(0.3, 1) + 0.05 Beta(1, 0.3) on the scaled coordinate (a posterior pinned at a band edge, as a wrong rival's
    is, keeps an effective sample); returns (x, log p(x) - log q(x))."""
    from scipy.stats import beta
    comp = rng.choice(3, size=n, p=(0.9, 0.05, 0.05))
    z = np.where(comp == 0, rng.uniform(0, 1, n), np.where(comp == 1, rng.beta(0.3, 1.0, n), rng.beta(1.0, 0.3, n)))
    z = np.clip(z, 1e-12, 1 - 1e-12)
    q = 0.9 + 0.05 * beta.pdf(z, 0.3, 1.0) + 0.05 * beta.pdf(z, 1.0, 0.3)
    return lo + (hi - lo) * z, -np.log(q)


def _defensive_normal(sd: float, n: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """N(0, sd) prior, proposal 0.9 N(0, sd) + 0.1 N(0, 4 sd); returns (x, log p - log q)."""
    wide = rng.random(n) < 0.1
    x = np.where(wide, rng.normal(0, 4 * sd, n), rng.normal(0, sd, n))
    lp = -0.5 * (x / sd) ** 2 - math.log(sd)
    lq = np.logaddexp(math.log(0.9) - 0.5 * (x / sd) ** 2 - math.log(sd),
                      math.log(0.1) - 0.5 * (x / (4 * sd)) ** 2 - math.log(4 * sd))
    return x, lp - lq


def draw(rival: str, n: int, rng: np.random.Generator, bands: Bands = Bands(), defensive: bool = True
         ) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Parameter draws of a rival: its declared bands (uniform in log10) and the measured design nuisances (the same
    for every rival: the step doses +- 0.01 dex, the heat-up 0.8 - 1.2 min, the trap blank 0.2 - 1 x its reporting
    limit, the Strecker-event readout's recovery 0.9 - 1.0). With ``defensive`` the draws come from the defensive
    proposal and the second return value is log p / q (the judge's ``log_prior_ratio``); else from the prior (zeros)."""
    th: Dict[str, np.ndarray] = {}
    lr = np.zeros(n)

    def u(name, b):
        nonlocal lr
        if b[0] == b[1]:
            th[name] = np.full(n, float(b[0]))
            return
        if defensive:
            x, w = _defensive_uniform(b[0], b[1], n, rng)
            lr = lr + w
        else:
            x = rng.uniform(b[0], b[1], n)
        th[name] = x

    def g(name, sd):
        nonlocal lr
        if defensive:
            x, w = _defensive_normal(sd, n, rng)
            lr = lr + w
        else:
            x = rng.normal(0, sd, n)
        th[name] = x
    if rival == "C0":
        u("log_kS", bands.log_kS_C0)
    elif rival == "C1":
        u("log_kS", bands.log_kS_C1)
        u("log_ks", bands.log_ks)
        u("dlog_k2", bands.log_k2_dex_around_elfsson)
    elif rival == "C3":
        u("log_kI", bands.log_kI)
        u("log_kS", bands.log_kS_C0)
    elif rival in ("S0", "S1", "S5"):
        u("log_kad_over_kS", bands.log_kad_over_kS)
    else:
        raise ValueError(rival)
    if rival.startswith("C"):
        g("dlog_m", 0.01)
        g("dlog_c", 0.01)
        u("heat", (0.8, 1.2))
        u("log_blank", (L10(0.2), 0.0))
        u("log_trap_recovery", (-1.0, 0.0))     # the trap's aminoacetone recovery, uncalibrated: 0.1 - 1
    else:
        u("log_recovery", (L10(0.9), 0.0))
    return th, lr


def _dmp_course(th: Mapping[str, np.ndarray], rival: str, c: float, m_mult_after: float, c_mult_after: float
                ) -> Dict[str, np.ndarray]:
    """DMP (M) at the onset times and 40 min, the AK pool at 40 min, and the DMP formed over 20 - 35 and 45 - 60 min
    when, at 40 min, MGO is raised x m_mult_after (Ala unchanged) or the whole pot is diluted x c_mult_after."""
    A = M = c
    th0 = th["heat"] * HEATUP_MIN
    m_mult_after = m_mult_after * 10 ** th["dlog_m"] if m_mult_after != 1.0 else 1.0
    c_mult_after = c_mult_after * 10 ** th["dlog_c"] if c_mult_after != 1.0 else 1.0
    if rival == "C0":
        kS = 10 ** th["log_kS"]
        rate = 0.5 * kS * A * M
        on = {t: rate * np.maximum(t - th0, 0.0) for t in ONSET_T}
        d40 = rate * (40 - th0)
        before = rate * 15.0
        A2, M2 = A * c_mult_after, M * m_mult_after * c_mult_after
        after = 0.5 * kS * A2 * M2 * 15.0
        return {"onset": {t: on[t] for t in ONSET_T}, "d40": d40, "ak40": np.zeros_like(kS), "before": before,
                "after": after}
    if rival == "C3":
        kI = 10 ** th["log_kI"]
        rate = kI * (A * M) ** 2
        on = {t: rate * np.maximum(t - th0, 0.0) for t in ONSET_T}
        A2, M2 = A * c_mult_after, M * m_mult_after * c_mult_after
        return {"onset": on, "d40": rate * (40 - th0), "ak40": np.zeros_like(kI), "before": rate * 15.0,
                "after": kI * (A2 * M2) ** 2 * 15.0}
    if rival == "C1":
        kS = 10 ** th["log_kS"]
        ks = 10 ** th["log_ks"]                         # at the pot's MGO (pseudo-first-order in M)
        k2 = elfsson_k2(100.0, POT_PH) * 10 ** th["dlog_k2"]
        F = kS * A * M
        on = {}
        ak = np.zeros_like(kS)
        P = np.zeros_like(kS)
        tprev = th0
        for t in ONSET_T:
            ak, dP = _rk4_ak(F, k2, ks, t - tprev, 0.25, ak)
            P = P + dP
            on[t] = P.copy()
            tprev = t
        ak20 = ak
        ak, dP = _rk4_ak(F, k2, ks, 15.0, 0.25, ak20)      # 20 - 35
        before = dP
        P35 = P + dP
        ak40, dP = _rk4_ak(F, k2, ks, 5.0, 0.25, ak)
        d40 = P35 + dP
        A2, M2 = A * c_mult_after, M * m_mult_after * c_mult_after
        F2 = kS * A2 * M2
        ks2 = ks * (M2 / M)
        ak45, _ = _rk4_ak(F2, k2, ks2, 5.0, 0.25, ak40 * c_mult_after)
        _, after = _rk4_ak(F2, k2, ks2, 15.0, 0.25, ak45)
        return {"onset": on, "d40": d40, "ak40": ak40, "before": before, "after": after}
    raise ValueError(rival)


def predict_c(th: Mapping[str, np.ndarray], rival: str) -> Dict[str, np.ndarray]:
    """The four part-C units (S x values each, dex)."""
    base = _dmp_course(th, rival, POT_C, 1.0, 1.0)
    mstep = _dmp_course(th, rival, POT_C, 3.0, 1.0)
    dstep = _dmp_course(th, rival, POT_C, 1.0, 0.25)
    tiny = 1e-300
    onset = np.log10(np.maximum(base["onset"][ONSET_T[0]], tiny) / base["d40"])[:, None]
    ak = np.log10((base["ak40"] * 10 ** th["log_trap_recovery"] + AK_REPORTING_LIMIT * 10 ** th["log_blank"])
                  / AK_REPORTING_LIMIT)
    return {"onset": onset,
            "mgo_step": np.log10(np.maximum(mstep["after"], tiny) / mstep["before"])[:, None],
            "dilution_step": np.log10(np.maximum(dstep["after"], tiny) / dstep["before"])[:, None],
            "ak_trap": ak[:, None]}


def predict_s(th: Mapping[str, np.ndarray], rival: str) -> Dict[str, np.ndarray]:
    """The two part-S units: log10 (Strecker events / amine consumed) without and with 0.4 M phosphate (early-time
    ratio of the two rates; dex)."""
    r = 10 ** th["log_kad_over_kS"]
    rec = th["log_recovery"]
    sink = rec - np.log10(1.0 + r)
    zero = rec
    if rival == "S0":
        return {"budget_plain": zero[:, None], "budget_phosphate": zero[:, None]}
    if rival == "S1":
        return {"budget_plain": sink[:, None], "budget_phosphate": sink[:, None]}
    if rival == "S5":
        return {"budget_plain": zero[:, None], "budget_phosphate": sink[:, None]}
    raise ValueError(rival)


# ---------------------------------------------------------------------------------------------------------------------
# The declared measurement model (sec. 7 item 5, as implemented in predictive_densities)
# ---------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class NoiseSpec:
    """Per-point within-run CV (dex) of one aliquot's reading, the number of replicate analyses per aliquot, and the
    instrument floor (dex). The within-run SD is the printed replicate SD of the design (n - 1 df per point, pooled
    over the unit's points); its estimate here is the declared design target."""
    cv_dex: float = 0.08
    n_rep: int = 3
    floor_dex: float = 0.01


def noise_models(units: Sequence[str], spec: NoiseSpec = NoiseSpec()) -> Dict[str, pdn.UnitNoiseModel]:
    """One UnitNoiseModel per designed unit: within-run level (a run term cancels), a printed within-run replicate SD
    of the mean of n_rep analyses (df = points x (n_rep - 1)), the instrument floor. A ratio of two readings carries
    sqrt(2) the per-reading SD; a rate ratio (four readings) carries 2x."""
    n_values = {"onset": 1, "mgo_step": 1, "dilution_step": 1, "ak_trap": 1, "budget_plain": 1, "budget_phosphate": 1}
    readings = {"onset": 2.0, "mgo_step": 4.0, "dilution_step": 4.0, "ak_trap": 2.0, "budget_plain": 2.0,
                "budget_phosphate": 2.0}
    out = {}
    for u in units:
        obs = f"m6_{u}"
        sd = spec.cv_dex * math.sqrt(readings[u] / spec.n_rep)
        df = float(n_values[u] * readings[u] * (spec.n_rep - 1))
        out[u] = pdn.UnitNoiseModel(u, obs, "within_run",
                                    (pdn.NoiseTerm("within_run_sd", sd, obs, df=df, source="design replicate SD"),
                                     pdn.NoiseTerm("instrument_floor", spec.floor_dex, obs)))
    return out


def misspecified_extra_sd(units: Sequence[str], spec: NoiseSpec, factor: float) -> Dict[str, float]:
    """The false-credit generator's extra per-point SD (beyond the scoring noise): the same-observable within-run
    reproducibility the design does not yet have, taken as ``factor`` x the declared per-point SD (the pyrazine lane has
    no within-laboratory term; a run offset cancels in every unit by construction)."""
    m = noise_models(units, spec)
    return {u: factor * float(m[u].terms[0].sd_array[0]) for u in units}


# ---------------------------------------------------------------------------------------------------------------------
# Judge plumbing
# ---------------------------------------------------------------------------------------------------------------------

def rival_draws(rival: str, units: Sequence[str], n: int, rng: np.random.Generator, bands: Bands = Bands(),
                defensive: bool = True) -> pdn.RivalDraws:
    """A rival's draws as the judge takes them: scoring ensembles from the defensive proposal (with log p / q), truths
    (``defensive=False``) from the prior."""
    th, lr = draw(rival, n, rng, bands, defensive)
    if not defensive:
        th_p = predict_c(th, rival) if rival.startswith("C") else predict_s(th, rival)
        return pdn.RivalDraws({u: th_p[u] for u in units})
    p = predict_c(th, rival) if rival.startswith("C") else predict_s(th, rival)
    return pdn.RivalDraws({u: p[u] for u in units}, log_prior_ratio=lr)


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi

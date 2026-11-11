"""
src/kinetic_core/ph_balance_rivals.py -- M1b PART (i): THE CHARGE-BALANCE RIVALS, EVERY ONE OFF BY DEFAULT.

Pre-registration: results/validation/kinetic_core_m1b_ph_balance_pentose_steps_prereg.md (v3.6 DRAFT), sec. 3.1 -
3.1.1 (the rivals), sec. 3.3 (C-Hou's measured ARP forcing) and sec. 3.4 (V2, V3, V7). Conflict CF-084.

WHAT THIS IS
============
The shipped charge balance (``ph_state``, B51) carries the carboxyl of every released amino acid into CBX and declares
its amine destroyed (``ph_state.AMINE_FATE_BASIS``); the acid of the lumped deoxyosone sink is ``ACID`` x a fitted
yield. CF-084 finds it too acid in Meynier's buffered cysteine pots and too basic in Hou's unbuffered ARP pot. Part (i)
scores four rival balances against it (B0 = the shipped one with its two in-sample constants integrated):

* **B1, amine survival** (``AmineSurvival``): a fraction of the amine the ledger destroys survives as ammonium (pKa
  9.25), ``phi_cys`` for cysteine thermolysis and its Strecker routes, ``phi_arp`` for the Amadori enolisations and the
  cysteine trap of 3-deoxypentosone (sec. 3.1's split, ``CYSTEINE_AMINE_STEPS`` / ``AMADORI_AMINE_STEPS``).
* **B2, acid lag** (``AcidLag``): the Amadori compound's loss returns its amino acid with fraction ``f_ret`` (as the
  zwitterion), and the carboxyl B0 carries at that loss goes to a latent, uncharged pool LCX released as CBX first
  order at ``k_lag`` (``AMADORI_LOSS_STEPS``).
* **B3, acid supply** (``AcidSupply``): a share y(pH) of the 1-deoxypentosone's (DPO's) loss goes to a latent,
  uncharged 2,4-diulose pool LAA, released as acetic acid (``AA``, one carboxyl) first order at ``k_iso``; y(pH) is
  linear in the IN-SITU pH between Davidek 2006a's three breakpoints (converted to in situ by the caller, V5) and flat
  outside. The branch rule is declared in two forms, both required agreement (sec. 3.1.1): ``"hold_total"`` (DPO's
  total loss is held, every other DPO sink scaled by 1 - y) and ``"added_sink"`` (the branch added on top at
  y / (1 - y) x the other sinks). The lumped ``ACID`` yield is the caller's to set to 0 (``PhDrift``).
* **B13** = B1 and B3 together.
* **C-Hou's forcing** (``ArpForcing``): the measured ARP(t) as an input. ARP's loss is forced to -dARP_meas/dt (a
  monotone PCHIP interpolant of the printed points) plus its gross production, shared over the rival's own ARP
  branches in their predicted proportions at each state; everything downstream is free (sec. 3.3).

NOTHING HERE MOVES A SHIPPED NUMBER. With every switch off ``sulfur.integrate_sulfur`` never enters this module (its
inline branch reads ``any_active()``), so the shipped engine is the same code, bit for bit. A rival's run is a
separate integration of the network PLUS five bookkeeping pools (``TRACKERS``); the network's species, stoichiometry
and ledger are unchanged, and the pools are reported in the run's metadata, never in its state matrix. Rivals run only
on the inline dynamic pH (``ph.coupling = "inline"``): any other mode is refused.

THE POOLS (none is a network species; none carries atoms the element audit counts)
* ``NH4_CYS``, ``NH4_ARP``: the surviving ammonium of B1, summed into the balance as ``NH4`` (``ph_state.AMMONIUM``).
* ``RET_AA``: B2's returned amino acid, titrated as glycine (``ph_state.GLYCINE``: the zwitterion between its two
  pKa, hence no net centre in the window these pots span). DECLARED READING: the prereg writes "as its zwitterion, no
  net centre"; titrating it as the zwitterion it is gives that in the window and is exact outside it.
* ``LCX``: B2's latent carboxyl (uncharged). ``LAA``: B3's latent 2,4-diulose (uncharged; its two acetic carbons
  stand outside the network's element audit while latent, so a B3 run's carbon closure drifts by 2 x LAA: reported).

DECLARED READINGS OF THE PREREG (flagged in the build record, not invented silently)
* B1's step sets follow sec. 3.1's wording; ``r_ttca_deg`` (TTCA's ring nitrogen) and ``r_hmf_cys`` (the HMF adduct)
  are named in neither set and stay at B0; ``r_cys_actz`` keeps its own basis (the thiazole ring). The CAMA
  enolisations already return cysteine, so no amine of theirs is destroyed to survive.
* B2's set is every Amadori loss that carries a carboxyl into CBX at B0 (the four enolisations and the lumped
  ``r_arp_decay``): "the carboxyl equivalent B0 carries at ARP loss".
* B2 with B1 (v1's B12) and B2 with B3 are not declared rivals and are refused.
* The added-sink form is undefined at y = 1 (an infinite rate); y is capped at ``ADDED_SINK_Y_CAP`` there and the
  run's metadata says whether the cap was reached.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

R_KJ = 8.314462618e-3
#: the temperature the prereg anchors k_lag and k_iso at (120 C; sec. 3.1: "log10 k_lag (120 C)", k_iso "at 120 C")
ANCHOR_T_K = 393.15

#: sec. 3.1, B1: cysteine thermolysis and its Strecker routes (phi_cys)
CYSTEINE_AMINE_STEPS: Tuple[str, ...] = (
    "r_cys_h2s", "r_cys_thermal", "ch_cys_ox", "r_dpo_ddp", "r_dpo_ddp_h2s", "r_mgo_cys_strecker", "r_mgo_cys_h2s")
#: sec. 3.1, B1: the Amadori / CAMA enolisations and r_tdp_cys_trap (phi_arp)
AMADORI_AMINE_STEPS: Tuple[str, ...] = (
    "r_arp_dpo", "r_arp_tdp", "r_arp_dpo_thermal", "r_arp_tdp_thermal", "r_arp_decay", "r_tdp_cys_trap")
#: sec. 3.1, B2: the Amadori losses whose carboxyl B0 carries into CBX
AMADORI_LOSS_STEPS: Tuple[str, ...] = (
    "r_arp_dpo", "r_arp_tdp", "r_arp_dpo_thermal", "r_arp_tdp_thermal", "r_arp_decay")
#: the bookkeeping pools appended to a rival's integration, in this order
TRACKERS: Tuple[str, ...] = ("NH4_CYS", "NH4_ARP", "RET_AA", "LCX", "LAA")
BRANCH_FORMS: Tuple[str, ...] = ("hold_total", "added_sink")
#: the added-sink form's rate is y / (1 - y) x the other sinks: undefined at y = 1
ADDED_SINK_Y_CAP = 1.0 - 1e-6
#: 2026-11-07 (the S3 application of M1b part (i); ``ph_balance_refusal``): an inline rival integration that needs more
#: than this many right-hand-side plus Jacobian evaluations is stopped (``RivalBudgetExceeded``), so a refusal for want
#: of a rival never depends on the machine's speed. Basis: the rival points that finish on the slowest recorded pots
#: need at most ~38 000 (0.5 M phosphate, pH 7.5, 145 C; 2026-11-07 probe); Chen 2000 / 2002's pH 5 ribose + cysteine
#: pots at 180 C passed 500 000 without finishing. Off for no run that finishes inside it: it changes no number.
#: 2026-11-11: the budget is per attempt; an LSODA run that exceeds it is retried once on BDF (``integrate_inline``).
RIVAL_RHS_BUDGET = 200_000


class RivalBudgetExceeded(RuntimeError):
    """An inline rival integration ran past ``RIVAL_RHS_BUDGET`` evaluations."""


def _check_fraction(name: str, value: float) -> None:
    if not (0.0 <= float(value) <= 1.0):
        raise ValueError(f"{name} = {value!r}: a fraction in [0, 1]")


@dataclass(frozen=True)
class AmineSurvival:
    """B1: the surviving fraction of the amine the ledger destroys, per step family (sec. 3.1).

    ``retained_steps`` (A2, prereg sec. 21.2): cysteine steps whose amine survives in FULL (fraction 1, whatever
    ``phi_cys`` is). It is the chemically right encoding of a step whose own chemistry keeps the amine: ``r_cys_thermal``
    (its lump's largest member is cystine, which keeps both centres) and ``ch_cys_ox`` (autoxidation to cystine). The
    shipped ledger declares both amines destroyed (``sulfur.CENTRE_LEDGER``); flipping that default moves shipped
    outputs, so the right encoding is opt-in here and the flip is a separate case."""

    enabled: bool = False
    phi_cys: float = 0.0
    phi_arp: float = 0.0
    retained_steps: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _check_fraction("phi_cys", self.phi_cys)
        _check_fraction("phi_arp", self.phi_arp)
        bad = [k for k in self.retained_steps if k not in CYSTEINE_AMINE_STEPS]
        if bad:
            raise ValueError(f"retained_steps {bad}: only cysteine amine steps {CYSTEINE_AMINE_STEPS}")


@dataclass(frozen=True)
class AcidLag:
    """B2: the returned amino-acid fraction and the latent carboxyl's release (k at 120 C, its barrier)."""

    enabled: bool = False
    f_ret: float = 0.0
    k_lag_per_min_120c: float = 10.0
    ea_lag_kj_mol: float = 103.7

    def __post_init__(self) -> None:
        _check_fraction("f_ret", self.f_ret)
        if not (float(self.k_lag_per_min_120c) > 0.0 and math.isfinite(float(self.k_lag_per_min_120c))):
            raise ValueError("k_lag_per_min_120c must be positive and finite")
        if not (float(self.ea_lag_kj_mol) >= 0.0):
            raise ValueError("ea_lag_kj_mol must be >= 0")


@dataclass(frozen=True)
class AcidSupply:
    """B3: the acetic share of DPO's loss at three in-situ breakpoints, the latent step, the branch form."""

    enabled: bool = False
    y5: float = 0.0
    y6: float = 0.0
    y_hi: float = 0.0
    #: Davidek 2006a's labels pH 5 / 6 / 6.5 (~0.16 M phosphate, 25 C) as IN-SITU pH (V5 prints the conversion)
    ph_in_situ_breakpoints: Tuple[float, ...] = ()
    k_iso_per_min_120c: float = 0.08
    ea_iso_kj_mol: float = 76.0
    branch_form: str = "hold_total"

    def __post_init__(self) -> None:
        for name in ("y5", "y6", "y_hi"):
            _check_fraction(name, getattr(self, name))
        if self.branch_form not in BRANCH_FORMS:
            raise ValueError(f"branch_form {self.branch_form!r}: one of {BRANCH_FORMS}")
        if not (float(self.k_iso_per_min_120c) > 0.0):
            raise ValueError("k_iso_per_min_120c must be positive")
        if self.enabled:
            if not (float(self.y5) <= float(self.y6) <= float(self.y_hi)):
                raise ValueError(f"the order y5 <= y6 <= y_hi is enforced (sec. 3.1): {self.y5, self.y6, self.y_hi}")
            bp = tuple(float(v) for v in self.ph_in_situ_breakpoints)
            if len(bp) != 3 or not (bp[0] < bp[1] < bp[2]):
                raise ValueError("three increasing in-situ breakpoints are required when B3 is enabled")

    def share(self, ph_in_situ: float) -> float:
        """y(pH): linear in the in-situ pH between the breakpoints, flat outside."""
        return float(np.interp(float(ph_in_situ), self.ph_in_situ_breakpoints, (self.y5, self.y6, self.y_hi)))


@dataclass(frozen=True)
class ArpForcing:
    """C-Hou: the measured ARP(t) imposed as an input (monotone PCHIP through the printed points, t = 0 included)."""

    enabled: bool = False
    times_min: Tuple[float, ...] = ()
    arp_mmol_l: Tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if self.enabled:
            t, a = np.asarray(self.times_min, float), np.asarray(self.arp_mmol_l, float)
            if t.size < 2 or t.shape != a.shape or t[0] != 0.0 or np.any(np.diff(t) <= 0):
                raise ValueError("forcing needs >= 2 points from t = 0, strictly increasing times")
            if np.any(np.diff(a) > 0) or np.any(a <= 0):
                raise ValueError("the forced ARP must be positive and non-increasing (a monotone loss)")


#: THE SWITCHES (config keys ph.amine_survival, ph.acid_lag, ph.acid_supply, ph.arp_forcing). All off.
AMINE_SURVIVAL = AmineSurvival()
ACID_LAG = AcidLag()
ACID_SUPPLY = AcidSupply()
ARP_FORCING = ArpForcing()


def any_active() -> bool:
    return bool(AMINE_SURVIVAL.enabled or ACID_LAG.enabled or ACID_SUPPLY.enabled or ARP_FORCING.enabled)


def check_combination() -> None:
    """Only the declared rivals: B1, B2, B3, B13 (= B1 + B3); the forcing goes with any. B2 + B1 (v1's B12) and
    B2 + B3 are not declared (sec. 3.1.1: 'B12 is not built')."""
    if ACID_LAG.enabled and (AMINE_SURVIVAL.enabled or ACID_SUPPLY.enabled):
        raise ValueError("B2 combined with B1 or B3 is not a declared rival (B12 is not built; no B23)")


def arrhenius_from_120c(k_120c: float, ea_kj_mol: float, temperature_k: float) -> float:
    """k(T) from its 120 C value and barrier (k(120 C) held)."""
    return float(k_120c) * math.exp(-float(ea_kj_mol) / R_KJ * (1.0 / float(temperature_k) - 1.0 / ANCHOR_T_K))


def balance_pools(trackers: Mapping[str, float]) -> Dict[str, float]:
    """The pools the charge balance reads (``ph_state.titratable_inventory``'s ``NH4`` and ``RET_AA``)."""
    return {"NH4": float(max(trackers.get("NH4_CYS", 0.0), 0.0)) + float(max(trackers.get("NH4_ARP", 0.0), 0.0)),
            "RET_AA": float(max(trackers.get("RET_AA", 0.0), 0.0))}


# ----------------------------------------------------------------------------------------------------------------
# The augmented system: the network's stoichiometry, plus the pools and three extra fluxes
# ----------------------------------------------------------------------------------------------------------------
#: the three extra fluxes, in order: B3's branch into LAA, B2's latent-carboxyl release, B3's acetic release
EXTRA_FLUXES: Tuple[str, ...] = ("b3_branch_to_laa", "b2_lcx_release", "b3_laa_release")


def _indices(keys: Sequence[str], all_keys: Sequence[str]) -> np.ndarray:
    pos = {k: i for i, k in enumerate(all_keys)}
    missing = [k for k in keys if k not in pos]
    if missing:
        raise KeyError(f"steps not in the network: {missing}")
    return np.array([pos[k] for k in keys], dtype=int)


def augmented_stoichiometry() -> Tuple[np.ndarray, Dict[str, Any]]:
    """S_aug ((n + 5) x (n_rx + 3)) for the switches in force: the network's matrix, B1's ammonium rows, B2's
    redirection of the Amadori losses' CBX, and the three extra fluxes' columns. Returns it with the index sets."""
    from . import sulfur as _s

    n, n_rx = _s.N_SULFUR_STATE, len(_s.FULL_REACTIONS)
    keys = _s.FULL_REACTION_KEYS
    s = np.zeros((n + len(TRACKERS), n_rx + len(EXTRA_FLUXES)))
    s[:n, :n_rx] = _s.SULFUR_STOICHIOMETRY
    t = {k: n + i for i, k in enumerate(TRACKERS)}
    idx = _s.SULFUR_INDEX
    if AMINE_SURVIVAL.enabled:
        for k, j in zip(CYSTEINE_AMINE_STEPS, _indices(CYSTEINE_AMINE_STEPS, keys)):
            s[t["NH4_CYS"], j] += 1.0 if k in AMINE_SURVIVAL.retained_steps else float(AMINE_SURVIVAL.phi_cys)
        for j in _indices(AMADORI_AMINE_STEPS, keys):
            s[t["NH4_ARP"], j] += float(AMINE_SURVIVAL.phi_arp)
    if ACID_LAG.enabled:
        for j in _indices(AMADORI_LOSS_STEPS, keys):
            carried = s[idx["CBX"], j]
            s[idx["CBX"], j] = 0.0
            s[t["LCX"], j] += (1.0 - float(ACID_LAG.f_ret)) * carried
            s[t["RET_AA"], j] += float(ACID_LAG.f_ret) * carried
    c = n_rx
    # B3's branch: DPO -> LAA + FRAG_C 3 (the triose partner; the two acetic carbons latent in LAA)
    s[idx["DPO"], c] -= 1.0
    s[idx["FRAG_C"], c] += 3.0
    s[t["LAA"], c] += 1.0
    # B2's release: LCX -> CBX
    s[t["LCX"], c + 1] -= 1.0
    s[idx["CBX"], c + 1] += 1.0
    # B3's release: LAA -> AA
    s[t["LAA"], c + 2] -= 1.0
    s[idx["AA"], c + 2] += 1.0
    dpo_row = _s.SULFUR_STOICHIOMETRY[idx["DPO"]]
    arp_row = _s.SULFUR_STOICHIOMETRY[idx["ARP"]]
    sets = {"n": n, "n_rx": n_rx, "tracker_index": t,
            "dpo_consumers": np.flatnonzero(dpo_row < 0.0), "dpo_coefficients": -dpo_row[dpo_row < 0.0],
            "arp_consumers": np.flatnonzero(arp_row < 0.0), "arp_coefficients": -arp_row[arp_row < 0.0],
            "arp_producers": np.flatnonzero(arp_row > 0.0), "arp_production": arp_row[arp_row > 0.0]}
    return s, sets


def _forcing_interpolant():
    from scipy.interpolate import PchipInterpolator

    f = PchipInterpolator(np.asarray(ARP_FORCING.times_min, float), np.asarray(ARP_FORCING.arp_mmol_l, float),
                          extrapolate=False)
    return f, f.derivative()


def integrate_inline(parameters, y0, eval_times, duration, temperature_of, sid, drift, spec, method, rtol, atol):
    """A rival's inline-coupled integration (``sulfur._integrate_inline_ph`` with the pools and fluxes above).

    Returns ``(solution, ph_hot, ionic, ph_cool, k_rows, trackers, diagnostics)``: ``solution.y`` holds the
    network's state only (as the shipped path's), ``trackers`` the pools at ``eval_times``, ``diagnostics`` the
    per-call checks V7 reads (the forced scale's sign, B3's held total, the added-sink cap)."""
    from scipy.integrate import solve_ivp

    from . import ph_kernel
    from . import sulfur as _s
    from .ph_state import (LABEL_TEMPERATURE_K, solve_ph_and_ionic_strength, solve_ph_and_ionic_strength_compiled,
                           solve_ph_and_ionic_strength_warm, titratable_inventory)
    from . import ph_state as _ph

    check_combination()
    warm_solve = solve_ph_and_ionic_strength_compiled if ph_kernel.COMPILED else solve_ph_and_ionic_strength_warm
    s_aug, sets = augmented_stoichiometry()
    n, n_rx, tix = sets["n"], sets["n_rx"], sets["tracker_index"]
    keys = _s._TITRATABLE_KEYS + (_s._TRUNK_SOLUTE_KEYS if _ph.TRUNK_SOLUTES_IN_BALANCE else ())
    index = [_s.SULFUR_INDEX[k] for k in keys]
    use_label = bool(_ph.B51.label_scale_tables)
    memo = {"ph": None, "i": 0.0, "cool": None}
    diag = {"rhs_calls": 0, "forced_scale_min": float("inf"), "forced_scale_max": 0.0, "forced_unresolvable_calls": 0,
            "dpo_total_rel_dev_max": 0.0, "added_sink_cap_reached": False, "supply_share_min": float("inf"),
            "supply_share_max": float("-inf")}
    supply, lag, forcing = ACID_SUPPLY, ACID_LAG, ARP_FORCING
    hold = supply.branch_form == "hold_total"
    arp_meas = d_arp_meas = None
    if forcing.enabled:
        arp_meas, d_arp_meas = _forcing_interpolant()
        t_last = float(forcing.times_min[-1])

    def pools(y):
        return balance_pools({k: y[tix[k]] for k in ("NH4_CYS", "NH4_ARP", "RET_AA")})

    def state_ph(y, temperature):
        conc = {k: float(max(y[i], 0.0)) for k, i in zip(keys, index)}
        conc.update(pools(y))
        inventory = titratable_inventory(conc, drift, spec)
        if memo["ph"] is None:
            ph, ionic = solve_ph_and_ionic_strength(temperature, sid, inventory)
        else:
            ph, ionic = warm_solve(temperature, sid, inventory, memo["ph"], memo["i"])
        cool = None
        if use_label:
            guess = memo["cool"] if memo["cool"] is not None else ph
            cool = warm_solve(LABEL_TEMPERATURE_K, sid, inventory, guess, ionic)[0]
        memo.update(ph=ph, i=ionic, cool=cool)
        return ph, ionic, cool

    def evaluate(t, y, want_jac):
        temperature = temperature_of(t)
        ph, ionic, cool = state_ph(y, temperature)
        kv = _s.rate_vector(_s.sulfur_rate_constants_at(parameters, temperature, float(ph), ph_label=cool,
                                                         ionic_strength_m=float(ionic)))
        y_main = y[:n]
        r = _s._rates_vectorised(y_main, kv)
        scale = np.ones(n_rx)
        if forcing.enabled:
            ac, aw = sets["arp_consumers"], sets["arp_coefficients"]
            loss = float(aw @ r[ac])
            production = float(sets["arp_production"] @ r[sets["arp_producers"]]) if sets["arp_producers"].size else 0.0
            slope = float(d_arp_meas(min(float(t), t_last))) if float(t) <= t_last else 0.0
            target = max(-slope, 0.0) + production
            if loss > 0.0:
                s_f = target / loss
            else:
                s_f = 0.0
                if target > 0.0:
                    diag["forced_unresolvable_calls"] += 1
            scale[ac] *= s_f
            diag["forced_scale_min"] = min(diag["forced_scale_min"], s_f)
            diag["forced_scale_max"] = max(diag["forced_scale_max"], s_f)
        branch, coef = 0.0, 0.0
        dc, dw = sets["dpo_consumers"], sets["dpo_coefficients"]
        if supply.enabled:
            yv = supply.share(ph)
            diag["supply_share_min"] = min(diag["supply_share_min"], yv)
            diag["supply_share_max"] = max(diag["supply_share_max"], yv)
            loss = float(dw @ r[dc])
            if hold:
                scale[dc] *= (1.0 - yv)
                coef = yv
                held = float(dw @ (r[dc] * scale[dc])) + coef * loss
                if loss > 0.0:
                    diag["dpo_total_rel_dev_max"] = max(diag["dpo_total_rel_dev_max"], abs(held - loss) / loss)
            else:
                if yv > ADDED_SINK_Y_CAP:
                    yv, diag["added_sink_cap_reached"] = ADDED_SINK_Y_CAP, True
                coef = yv / (1.0 - yv)
            branch = coef * loss
        k_lag = arrhenius_from_120c(lag.k_lag_per_min_120c, lag.ea_lag_kj_mol, temperature) if lag.enabled else 0.0
        k_iso = (arrhenius_from_120c(supply.k_iso_per_min_120c, supply.ea_iso_kj_mol, temperature)
                 if supply.enabled else 0.0)
        r_aug = np.concatenate([r * scale, [branch, k_lag * max(y[tix["LCX"]], 0.0),
                                            k_iso * max(y[tix["LAA"]], 0.0)]])
        if not want_jac:
            return s_aug @ r_aug
        jr = _s._rate_jacobian(y_main, kv)
        j_aug = np.zeros((n_rx + len(EXTRA_FLUXES), n + len(TRACKERS)))
        j_aug[:n_rx, :n] = jr * scale[:, None]
        if supply.enabled:
            j_aug[n_rx, :n] = coef * (dw @ jr[dc])
        j_aug[n_rx + 1, tix["LCX"]] = k_lag
        j_aug[n_rx + 2, tix["LAA"]] = k_iso
        return s_aug @ j_aug

    def spend():
        diag["evaluations"] = diag.get("evaluations", 0) + 1
        if diag["evaluations"] > RIVAL_RHS_BUDGET:
            raise RivalBudgetExceeded(f"the carried rival's inline integration ({diag['method']}) needed more than "
                                      f"RIVAL_RHS_BUDGET = {RIVAL_RHS_BUDGET} right-hand-side and Jacobian evaluations")

    def rhs(t, y):
        diag["rhs_calls"] += 1
        spend()
        return evaluate(t, y, False)

    def jac(t, y):
        spend()
        return evaluate(t, y, True)

    y_aug0 = np.concatenate([np.asarray(y0, float), np.zeros(len(TRACKERS))])
    # 2026-11-11 (platform robustness; SUG-01): an LSODA run that exceeds the budget or fails is re-integrated ONCE with
    # BDF under a fresh budget, as fit_m1b_ph_balance.integrate_with_fallback does for the screen's pots; the fallback
    # is recorded in the diagnostics. With the pH re-solved inside the right-hand side, ~1 % of LSODA runs of a pH 5 pot
    # lock in the first ~1e-6 min into order-1 Adams steps of ~1e-13 min that never switch to BDF (no Jacobian call;
    # the Jacobian's eigenvalues there are > -1e4 /min, so the cap is not the problem's stiffness). Which run locks is
    # chaotic in the last bits: a 1e-8 relative change of a precursor, or the platform's arithmetic, moves it (SUG-01's
    # glucose arm, B1_d345, locks on ubuntu x86-64 and not on macOS arm64; 8 of 920 runs of SUG-01's two arms under
    # 3e-8 - 6e-7 perturbations lock on macOS, 0 of 920 on BDF). A run that finishes on LSODA is untouched, bit for bit;
    # BDF agrees with LSODA to < 4e-6 dex on every rival output of 294 S3 recipes.
    fresh = dict(diag)
    attempts = [method] + (["BDF"] if method == "LSODA" else [])
    fallback = None
    for attempt, m in enumerate(attempts):
        diag.clear()
        diag.update(fresh, method=m)
        memo.update(ph=None, i=0.0, cool=None)
        kwargs = {"jac": jac} if m in ("BDF", "Radau", "LSODA") else {}
        last = attempt == len(attempts) - 1
        try:
            sol = solve_ivp(rhs, (0.0, duration), y_aug0, t_eval=eval_times, method=m, rtol=rtol, atol=atol, **kwargs)
        except RivalBudgetExceeded as exc:
            if last:
                raise
            fallback = {"from": m, "to": attempts[attempt + 1], "reason": str(exc)[:200]}
            continue
        if sol.success:
            break
        if last:
            raise RuntimeError(f"sulfur-network integration (M1b rival) failed: {sol.message}")
        fallback = {"from": m, "to": attempts[attempt + 1], "reason": str(sol.message)[:200]}
    diag["solver_fallback"] = fallback
    states = np.clip(sol.y.T, 0.0, None)
    memo.update(ph=None, i=0.0, cool=None)
    ph_hot, ionic, ph_cool, k_rows = [], [], [], []
    for t, y in zip(eval_times, states):
        temperature = temperature_of(float(t))
        p, i, c = state_ph(y, temperature)
        ph_hot.append(p), ionic.append(i), ph_cool.append(c)
        k_rows.append(_s.rate_vector(_s.sulfur_rate_constants_at(parameters, temperature, float(p), ph_label=c,
                                                                  ionic_strength_m=float(i))))
    trackers = {k: states[:, tix[k]].copy() for k in TRACKERS}
    if forcing.enabled:
        inside = np.asarray(eval_times, float) <= t_last
        meas = np.full(len(eval_times), np.nan)
        meas[inside] = arp_meas(np.asarray(eval_times, float)[inside])
        got = states[:, _s.SULFUR_INDEX["ARP"]]
        rel = np.abs(got[inside] - meas[inside]) / meas[inside]
        diag["forced_arp_max_rel_error"] = float(np.max(rel)) if rel.size else None
        diag["forced_arp_measured_interpolant"] = [None if not math.isfinite(v) else float(v) for v in meas]
        diag["forced_loss_never_negative"] = bool(diag["forced_scale_min"] >= 0.0)
    for k in ("forced_scale_min", "supply_share_min", "supply_share_max"):
        if not math.isfinite(diag[k]):
            diag[k] = None
    main = SimpleNamespace(y=sol.y[:n], t=sol.t, success=sol.success, message=sol.message)
    return main, np.array(ph_hot), np.array(ionic), ph_cool, np.array(k_rows), trackers, diag


def switches_record() -> Dict[str, Any]:
    """The switches in force, for a run's metadata."""
    from dataclasses import asdict

    return {"amine_survival": asdict(AMINE_SURVIVAL), "acid_lag": asdict(ACID_LAG),
            "acid_supply": asdict(ACID_SUPPLY),
            "arp_forcing": {"enabled": ARP_FORCING.enabled, "n_points": len(ARP_FORCING.times_min)}}


__all__ = ["RIVAL_RHS_BUDGET", "RivalBudgetExceeded", "ACID_LAG", "ACID_SUPPLY", "AMADORI_AMINE_STEPS", "AMADORI_LOSS_STEPS", "AMINE_SURVIVAL", "ARP_FORCING",
           "AcidLag", "AcidSupply", "AmineSurvival", "ArpForcing", "BRANCH_FORMS", "CYSTEINE_AMINE_STEPS", "TRACKERS",
           "any_active", "arrhenius_from_120c", "augmented_stoichiometry", "balance_pools", "check_combination",
           "integrate_inline", "switches_record"]

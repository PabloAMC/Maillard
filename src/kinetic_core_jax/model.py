"""
src/kinetic_core_jax/model.py -- the DYNAMIC half of the port, in JAX (float64): the charge balance, the rate
constants, the right-hand side and the diffrax integration of one isothermal pot.

THE pH COUPLING (the hard part). The engine's inline mode solves the charge balance at every right-hand-side call:
the in-situ (pH, I) at the reaction temperature, with I the Davies fixed point, and the cooled 25 C label pH that the
measured pH tables read (B51). Two ways to put that into a differentiable ODE were open:

* pH as an algebraic STATE (an index-1 DAE). diffrax has no mass-matrix / DAE solver, so this would mean a hand-made
  DAE integrator or a stiff relaxation (an artificial d pH/dt = -F/eps), which changes the model.
* the root solved INSIDE the right-hand side, differentiated by the implicit-function theorem. Chosen: it is exactly
  the engine's physics (the same algebraic constraint at every call, so the trajectories are comparable term by
  term), it keeps every diffrax ODE solver usable, and its derivative is exact: for F(u; c) = 0 with u = (pH, I) and c
  the titratable inventory, du = -(dF/du)^-1 (dF/dc) dc, a 2 x 2 solve (``solve_ph``'s custom JVP). Reverse mode
  transposes that linear rule, so adjoints work too. The forward solve is a bisection (the residual is monotone in
  pH: unconditionally safe, no warm start needed in a pure function) followed by Newton on the coupled (pH, I).

Everything else mirrors ``src/kinetic_core/sulfur.sulfur_rate_constants_at`` and ``ph_state`` line by line; the unit
tests compare both against the engine on grids.
"""
from __future__ import annotations

from typing import Any, Dict, Sequence

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from .lane import Lane  # noqa: E402

LN10 = float(np.log(10.0))
BISECTION_STEPS = 56      # [-1, 15] to 2e-16 pH (the SID's and the tests' full bisection)
NEWTON_STEPS = 6
COARSE_BISECTION_STEPS = 14   # (-1, 15) to 1e-3 pH before Newton
#: unroll the root-find's fixed-count loops into straight-line code (read at trace time). Measured 2026-11-04: no gain
#: inside a solve (1.22 vs 1.18 s), so off.
UNROLL = False


# --------------------------------------------------------------------------- pH physics (ph_state, B51 all on)
def pk_at(pk25, dh, dcp, t_k, r_kj, t0=298.15):
    """ph_state.pk_at with the heat-capacity term (B51)."""
    out = pk25 + dh / (r_kj * LN10) * (1.0 / t_k - 1.0 / t0)
    return out - (1e-3 * dcp / (r_kj * LN10)) * (jnp.log(t_k / t0) + t0 / t_k - 1.0)


def davies_a(t_k, c):
    t_c = t_k - 273.15
    eps = 87.740 - 0.40008 * t_c + 9.398e-4 * t_c ** 2 - 1.410e-6 * t_c ** 3
    rho = jnp.interp(t_c, jnp.asarray(c["rho_table_c"]), jnp.asarray(c["rho_table"]))
    return 1.82483e6 * jnp.sqrt(rho) / (eps * t_k) ** 1.5


def davies_term(ionic_m, t_k, c):
    safe = jnp.maximum(ionic_m, 1e-300)
    s = jnp.sqrt(safe)
    return jnp.where(ionic_m > 0.0, davies_a(t_k, c) * (s / (1.0 + s) - 0.3 * safe), 0.0)


def pkw_at(t_k, c):
    """Linear in 1/T on the tabulated pKw, clamped at the table's ends (ph_state.pkw_at)."""
    x = 1.0 / jnp.asarray(c["pkw_table_k"])[::-1]
    return jnp.interp(1.0 / t_k, x, jnp.asarray(c["pkw_table"])[::-1])


def temperature_terms(t_k, lane: Lane) -> Dict[str, Any]:
    """What the charge balance needs that depends on T alone -- every group's pK at I = 0 (van 't Hoff + dCp),
    pKw, Davies's A -- computed once per right-hand-side call, not in every iteration of the root-find."""
    c = lane.const
    pk0 = pk_at(jnp.asarray(lane.group_pk25), jnp.asarray(lane.group_dh), jnp.asarray(lane.group_dcp), t_k, c["R_KJ"])
    return {"pk0": pk0, "pkw": pkw_at(t_k, c), "a": davies_a(t_k, c)}


def _davies(ionic_m, a):
    safe = jnp.maximum(ionic_m, 1e-300)
    s = jnp.sqrt(safe)
    return jnp.where(ionic_m > 0.0, a * (s / (1.0 + s) - 0.3 * safe), 0.0)


def _group_pk(ionic_m, tt, lane: Lane):
    """Every group's mixed pK at (T, I) (TitratableGroup.pkas_at: pK' = pK + A f(I) (2 z_HA - 1)), and the Davies
    term: constant in pH, so a bisection at fixed I computes them once."""
    af = _davies(ionic_m, tt["a"])
    zmax = jnp.asarray(lane.group_zmax)[:, None]
    return tt["pk0"] + af * (2.0 * (zmax - jnp.arange(4.0)) - 1.0), af


def _group_terms(ph, pk, lane: Lane):
    """Per group: mean charge and SUM_k x_k z_k^2 (TitratableGroup.average_charge / ionic_strength_contribution),
    from log10 weights carried relative to their maximum (no overflow)."""
    log_w = jnp.concatenate([jnp.zeros((pk.shape[0], 1)), jnp.cumsum(ph - pk, axis=1)], axis=1)
    w = 10.0 ** (log_w - jnp.max(log_w, axis=1, keepdims=True))
    w = w / jnp.sum(w, axis=1, keepdims=True)
    i = jnp.arange(5.0)
    z = jnp.asarray(lane.group_zmax)
    return z - jnp.sum(i * w, axis=1), jnp.sum(w * (z[:, None] - i) ** 2, axis=1)


def _residual_given(ph, pk, af, pkw, sid, conc, lane):
    charge, ionic = _group_terms(ph, pk, lane)
    inv_g = 10.0 ** af
    h = 1000.0 * 10.0 ** (-ph) * inv_g
    oh = 1000.0 * 10.0 ** (ph - pkw) * inv_g
    res = sid + h - oh + jnp.sum(conc * charge)
    i_new = 0.5 * (jnp.abs(sid) + h + oh + jnp.sum(conc * ionic)) / 1000.0
    return res, i_new


def _residual_tt(ph, ionic_m, tt, sid, conc, lane):
    pk, af = _group_pk(ionic_m, tt, lane)
    return _residual_given(ph, pk, af, tt["pkw"], sid, conc, lane)


def charge_residual(ph, ionic_m, t_k, sid, conc, lane: Lane):
    """SID + [H+] - [OH-] + SUM c z, mmol/L (ph_state._charge_residual), and I(pH, I) in mol/L
    (ph_state.ionic_strength_of)."""
    return _residual_tt(ph, ionic_m, temperature_terms(t_k, lane), sid, conc, lane)


def _equations(u, tt, sid, conc, lane):
    res, i_new = _residual_tt(u[0], u[1], tt, sid, conc, lane)
    return jnp.stack([res, u[1] - i_new])


def _bisect(tt, sid, conc, ionic_m, lane, steps=BISECTION_STEPS):
    """Bisection on the monotone residual over (-1, 15) at fixed I (ph_state._bisect_ph; at an edge it returns the
    edge, as the engine does)."""
    pk, af = _group_pk(ionic_m, tt, lane)

    def body(_, bounds):
        lo, hi = bounds
        mid = 0.5 * (lo + hi)
        up = _residual_given(mid, pk, af, tt["pkw"], sid, conc, lane)[0] > 0.0
        return jnp.where(up, mid, lo), jnp.where(up, hi, mid)

    lo, hi = jax.lax.fori_loop(0, steps, body, (jnp.asarray(-1.0), jnp.asarray(15.0)), unroll=UNROLL)
    return 0.5 * (lo + hi)


def _solve_forward(tt, sid, conc, lane):
    """A coarse bisection at I = 0 (to 1e-3 pH), one ionic-strength sweep, then Newton on the coupled (pH, I), each
    pH step limited to one unit: the same root as the engine's bisection / Brent + I fixed point."""
    ph0 = _bisect(tt, sid, conc, 0.0, lane, steps=COARSE_BISECTION_STEPS)
    u = jnp.stack([ph0, _residual_tt(ph0, 0.0, tt, sid, conc, lane)[1]])

    def newton(_, u):
        f = _equations(u, tt, sid, conc, lane)
        jac = jax.jacfwd(lambda uu: _equations(uu, tt, sid, conc, lane))(u)
        step = jnp.linalg.solve(jac, f)
        d_ph = jnp.clip(step[0], -1.0, 1.0)
        return jnp.stack([u[0] - d_ph, jnp.maximum(u[1] - step[1], 0.1 * u[1])])

    return jax.lax.fori_loop(0, NEWTON_STEPS, newton, u, unroll=UNROLL)


def make_ph_solver(lane: Lane):
    """(pH, I) = solve(T, SID, inventory): forward by bisection + Newton, derivative by the implicit function."""

    @jax.custom_jvp
    def solve_tt(tt, sid, conc):
        return _solve_forward(jax.lax.stop_gradient(tt), jax.lax.stop_gradient(sid),
                              jax.lax.stop_gradient(conc), lane)

    @solve_tt.defjvp
    def solve_jvp(primals, tangents):
        tt, sid, conc = primals
        u = solve_tt(tt, sid, conc)
        jac_u = jax.jacfwd(lambda uu: _equations(uu, tt, sid, conc, lane))(u)
        _, d_f = jax.jvp(lambda t, s, cc: _equations(u, t, s, cc, lane), (tt, sid, conc), tangents)
        return u, -jnp.linalg.solve(jac_u, d_f)

    def solve(t_k, sid, conc):
        return solve_tt(temperature_terms(t_k, lane), sid, conc)

    return solve



def strong_ion_difference(label_ph, conc0, lane: Lane, sweeps: int = 60):
    """ph_state.strong_ion_difference: the SID that balances the declared label pH at 25 C, with I solved with it."""
    t_l = lane.const["label_t_k"]
    sid = -charge_residual(label_ph, 0.0, t_l, 0.0, conc0, lane)[0]

    def body(_, state):
        sid, ionic = state
        fresh = charge_residual(label_ph, ionic, t_l, sid, conc0, lane)[1]
        return -charge_residual(label_ph, fresh, t_l, 0.0, conc0, lane)[0], fresh

    return jax.lax.fori_loop(0, sweeps, body, (sid, jnp.asarray(0.0)))[0]


def inventory(y, lane: Lane, buffer=None):
    """The titratable inventory (mmol/L per group) of a state: the lane's linear map plus the buffer's own groups
    (``buffer``: a pot's, as an argument, so one compiled function serves pots with different buffers)."""
    offset = jnp.asarray(lane.group_offset) if buffer is None else buffer
    return jnp.asarray(lane.group_weights) @ jnp.maximum(y, 0.0) + offset


# --------------------------------------------------------------------------- rate constants
def _speciation_pk(name, t_k, ionic_m, c):
    """parameters_sulfur.pka_at under B51 (Suleimenov H2S, GKL cysteine pair), minus Davies."""
    if name == "h2s_1":
        log_k1 = (782.43945 + 0.361261 * t_k - 1.6722e-4 * t_k * t_k - 20565.7315 / t_k - 142.741722 * jnp.log(t_k))
        pk = -log_k1
    else:
        pk25, dh, dcp = c["pk_thiol"] if name == "thiol" else c["pk_amine"]
        pk = pk_at(pk25, dh, dcp, t_k, c["R_KJ"], c["pka_ref_t_k"])
    return pk - davies_term(ionic_m, t_k, c)


def ph_factors(ph, t_k, ionic_m, c):
    """parameters_sulfur.ph_factor for every kind, as a vector indexed by lane.KIND_CODES."""
    p5 = c["ph_norm"]
    amine = lambda x, i: 1.0 / (1.0 + 10.0 ** (_speciation_pk("cysteine_amine", t_k, i, c) - x))  # noqa: E731
    neutral = lambda x, i: 1.0 / (1.0 + 10.0 ** (x - _speciation_pk("h2s_1", t_k, i, c)))  # noqa: E731
    thiolate = lambda x, i: 1.0 / (1.0 + 10.0 ** (_speciation_pk("thiol", t_k, i, c) - x))  # noqa: E731
    return jnp.stack([
        jnp.asarray(1.0),
        10.0 ** (p5 - ph),
        amine(ph, ionic_m) / amine(p5, 0.0),
        neutral(ph, ionic_m) / neutral(p5, 0.0),
        (1.0 - neutral(ph, ionic_m)) / (1.0 - neutral(p5, 0.0)),
        thiolate(ph, ionic_m) / thiolate(p5, 0.0),
        10.0 ** (ph - p5),
    ])


def cysteine_thermolysis_k(ph_label, t_k, c):
    pts, la, ea = c["zheng_ph"], c["zheng_log10_a"], c["zheng_ea"]
    x = jnp.minimum(jnp.maximum(ph_label, pts[1]), pts[-1])
    lo = x <= pts[2]
    w = jnp.where(lo, (x - pts[1]) / (pts[2] - pts[1]), (x - pts[2]) / (pts[3] - pts[2]))
    log_a = jnp.where(lo, la[1] + w * (la[2] - la[1]), la[2] + w * (la[3] - la[2]))
    e = jnp.where(lo, ea[1] + w * (ea[2] - ea[1]), ea[2] + w * (ea[3] - ea[2]))
    below = ph_label < c["zheng_boundary"]
    log_a, e = jnp.where(below, la[0], log_a), jnp.where(below, ea[0], e)
    return 10.0 ** log_a * jnp.exp(-e / (c["R_KJ"] * t_k))


def glutathione_thermolysis_k(ph_label, t_k, c):
    (l3, e3), (l5, e5), (l7, e7), (l9, e9) = c["gsh_pairs"]
    x = jnp.minimum(ph_label, 9.0)
    lo = x <= 7.0
    w = jnp.where(lo, (x - 5.0) / 2.0, (x - 7.0) / 2.0)
    ln_a = jnp.where(lo, l5 + w * (l7 - l5), l7 + w * (l9 - l7))
    e = jnp.where(lo, e5 + w * (e7 - e5), e7 + w * (e9 - e7))
    below = ph_label < c["zheng_boundary"]
    ln_a, e = jnp.where(below, l3, ln_a), jnp.where(below, e3, e)
    return jnp.exp(ln_a) * jnp.exp(-e / (c["R_KJ"] * t_k))


def thiamine_h2s_share(ph_label, c):
    (p0, f0), (p1, f1) = c["thi_share_points"]
    mid = 10.0 ** (np.log10(f0) + (ph_label - p0) / (p1 - p0) * (np.log10(f1) - np.log10(f0)))
    return jnp.where(ph_label < c["thi_share_floor"], 0.0, jnp.where(ph_label >= p1, f1, mid))


def cama_decay_at_ref(ph, c):
    phs, logs = jnp.asarray(c["cama_decay_ph"]), jnp.asarray(c["cama_decay_log10k"])
    slope = (logs[-1] - logs[-2]) / (phs[-1] - phs[-2])
    return 10.0 ** jnp.where(ph > phs[-1], logs[-1] + slope * (ph - phs[-1]), jnp.interp(ph, phs, logs))


def _arrhenius_100c(k_100, ea, t_k, c):
    return k_100 * jnp.exp(-ea / c["R_KJ"] * (1.0 / t_k - 1.0 / 373.15))


def coordinate_values(theta, lane: Lane) -> Dict[str, Any]:
    """The named coordinates from the parameter vector (or their fixed values when not free)."""
    out = dict(lane.const["fixed"])
    for i, key in enumerate(lane.coordinate_keys):
        out[key] = theta[i]
    return out


def rate_constants(theta, t_k, ph, ph_label, ionic_m, lane: Lane):
    """Every reaction's constant (sulfur.sulfur_rate_constants_at, R3En configuration), in FULL_REACTIONS order."""
    c, ix = lane.const, lane.index
    t_eff = jnp.minimum(t_k, jnp.asarray(lane.cap_k))
    arr = jnp.asarray(lane.k_ref) * jnp.exp(-jnp.asarray(lane.ea_over_r) * (1.0 / t_eff - jnp.asarray(lane.inv_t_ref)))
    arr = jnp.where(jnp.asarray(lane.ea_over_r) > 0.0, arr, jnp.asarray(lane.k_ref))
    k = arr * 10.0 ** (jnp.asarray(lane.offset_matrix) @ theta) * ph_factors(ph, t_k, ionic_m, c)[lane.kind]
    v = coordinate_values(theta, lane)
    s = {}
    s["r_cys_h2s"] = cysteine_thermolysis_k(ph_label, t_k, c)
    s["r_gsh_h2s"] = glutathione_thermolysis_k(ph_label, t_k, c)
    share = thiamine_h2s_share(ph_label, c)
    s["r_thi_h2s"] = jnp.where(share > 0.0, share / (1.0 - share), 0.0) * (k[ix["r_thi_hmp"]] + k[ix["r_thi_mesh"]])
    s["r_ttca_form"] = c["ttca_k"] * k[ix["r_ttca_cys"]]
    s["r_cama_form"] = c["cama_k_form_ref"] * jnp.exp(-c["cama_ea_form"] / c["R_KJ"]
                                                      * (1.0 / t_k - 1.0 / c["cama_t_form_ref_k"]))
    decay = cama_decay_at_ref(ph, c) * jnp.exp(-c["cama_ea_decay"] / c["R_KJ"]
                                               * (1.0 / t_k - 1.0 / c["cama_t_decay_ref_k"]))
    s["r_cama_dpo"] = c["cama_dpo_share"] * decay
    s["r_cama_tdp"] = (1.0 - c["cama_dpo_share"]) * decay
    s["r_pent_dpo"] = 0.0
    s["r_pent_tdp"] = 0.0
    ddp = k[ix["r_dpo_ddp"]]
    s["r_dpo_ddp_h2s"] = c["strecker_share"] * ddp
    s["r_dpo_ddp"] = (1.0 - c["strecker_share"]) * ddp
    nf = (10.0 ** c["nf_log10_k_ph5"] * jnp.exp(-c["nf_ea"] / c["R_KJ"] * (1.0 / t_k - 1.0 / c["nf_t_anchor_k"]))
          * 10.0 ** (ph - c["ph_norm"]))
    s["r_nf_decay_oh_mgo"] = v["mgo_share_nf"] * nf
    s["r_nf_decay_oh"] = (1.0 - v["mgo_share_nf"]) * nf
    for key in ("r_tdp_mgo", "r_pent_mgo", "r_mgo_cys_strecker", "r_mgo_cys_h2s", "r_tdp_cys_trap", "r_odg_dpo",
                "ch_oligomer_mft", "ch_oligomer_fft"):
        s[key] = 0.0
    s["r_pent_amine"] = (_arrhenius_100c(10.0 ** v["log10_k_pent_amine"], c["entry_ea"], t_k, c)
                         * 10.0 ** (v["entry_ph_slope"] * (ph - 6.0)))
    ret = jnp.clip(v["amine_return"], 0.0, 1.0)
    dpo = c["pama_dpo_factor"] * (k[ix["r_arp_dpo"]] + k[ix["r_arp_dpo_thermal"]])
    tdp = c["pama_tdp_factor"] * (k[ix["r_arp_tdp"]] + k[ix["r_arp_tdp_thermal"]])
    s["r_pama_dpo"], s["r_pama_dpo_n"] = ret * dpo, (1.0 - ret) * dpo
    s["r_pama_tdp"], s["r_pama_tdp_n"] = ret * tdp, (1.0 - ret) * tdp
    s["r_dpo_aa"] = (_arrhenius_100c(10.0 ** v["log10_k_dpo_aa"], c["aa_ea"], t_k, c)
                     * 10.0 ** (v["aa_ph_slope"] * (ph - 6.0)))
    k_adduct = c["adduct_k_ref"] * jnp.exp(-(c["adduct_dh"] / c["R_KJ"]) * (1.0 / t_k - 1.0 / c["adduct_t_k"]))
    s["ch_thioether_release"] = k[ix["ch_thioether_mft"]] / k_adduct
    s["ch_thioether_release_fft"] = k[ix["ch_thioether_fft"]] / k_adduct
    s["r_amadori"] = c["amadori_factor"] * k[ix["r_schiff"]]
    keys = sorted(s)
    return k.at[jnp.array([ix[key] for key in keys])].set(jnp.stack([jnp.asarray(s[key], dtype=k.dtype)
                                                                      for key in keys]))


# --------------------------------------------------------------------------- the right-hand side and the solve
def mass_action(y, lane: Lane):
    yc = jnp.append(jnp.maximum(y, 0.0), 1.0)[jnp.asarray(lane.reactant_index)]
    p = jnp.asarray(lane.reactant_power)
    return jnp.prod(jnp.where(p == 2, yc * yc, jnp.where(p == 1, yc, 1.0)), axis=1)


class Pot:
    """One isothermal pot on a Lane. The pot's own numbers -- temperature, SID, buffer, initial state -- are ARGUMENTS
    of the compiled functions (``self.data``, ``self.y0``), not constants baked into them, so one compilation serves
    every pot of the same network: a fit compiles once per configuration, not once per pot."""

    def __init__(self, lane: Lane, y0: np.ndarray, t_k: float, label_ph: float):
        self.lane = lane
        self.solve_ph = make_ph_solver(lane)
        self.y0 = jnp.asarray(y0, dtype=jnp.float64)
        self.data = pot_data(lane, self.y0, t_k, label_ph)

    def ph_state(self, y, data=None):
        """(in-situ pH, I, cooled label pH, its I) of a state."""
        data = self.data if data is None else data
        conc = inventory(y, self.lane, data["buffer"])
        hot = self.solve_ph(data["t_k"], data["sid"], conc)
        cool = self.solve_ph(jnp.asarray(self.lane.const["label_t_k"]), data["sid"], conc)
        return hot[0], hot[1], cool[0], cool[1]

    def rhs(self, t, y, args):
        theta, data = args
        ph, ionic, ph_label, _ = self.ph_state(y, data)
        k = rate_constants(theta, data["t_k"], ph, ph_label, ionic, self.lane)
        return jnp.asarray(self.lane.stoichiometry) @ (k * mass_action(y, self.lane))

    def solve(self, theta, times: Sequence[float], *, data=None, y0=None, rtol=1e-10, atol=1e-16, adjoint=None,
              max_steps=200_000, solver=None, dt0=1e-7, pcoeff=0.3, icoeff=0.4):
        import diffrax

        ts = jnp.asarray(times, dtype=jnp.float64)
        return diffrax.diffeqsolve(
            diffrax.ODETerm(self.rhs), solver or diffrax.Kvaerno5(), t0=0.0, t1=ts[-1], dt0=dt0,
            y0=self.y0 if y0 is None else y0, args=(theta, self.data if data is None else data),
            saveat=diffrax.SaveAt(ts=ts),
            stepsize_controller=diffrax.PIDController(rtol=rtol, atol=atol, pcoeff=pcoeff, icoeff=icoeff),
            max_steps=max_steps, adjoint=adjoint or diffrax.RecursiveCheckpointAdjoint(), throw=True)

    def observables(self, ys, obs_index: Sequence[int], obs_rows: Sequence[int], data=None):
        """log10 of the observed species at the observed time rows, then the cooled pH at the last row."""
        vals = [jnp.log10(jnp.maximum(ys[r, i], 1e-300)) for i in obs_index for r in obs_rows]
        vals.append(self.ph_state(ys[obs_rows[-1]], data)[2])
        return jnp.stack(vals)


def pot_data(lane: Lane, y0, t_k: float, label_ph: float, buffer=None) -> Dict[str, Any]:
    """A pot's arguments: temperature, buffer groups, and the SID its declared label pH implies."""
    buffer = jnp.asarray(lane.group_offset) if buffer is None else jnp.asarray(buffer)
    conc0 = inventory(jnp.asarray(y0), lane, buffer)
    return {"t_k": jnp.asarray(float(t_k)), "buffer": buffer,
            "sid": strong_ion_difference(jnp.asarray(float(label_ph)), conc0, lane)}


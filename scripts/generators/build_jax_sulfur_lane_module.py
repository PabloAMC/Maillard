#!/usr/bin/env python
"""
build_jax_sulfur_lane_module.py -- the sulfur lane (``sulfur.FULL_REACTIONS``, 132 steps, 62 species, dynamic pH) as
an A2 contract module in JAX + diffrax (``src/kinetic_core_jax/sulfur_lane.py``; plan step A2 + H3; 2026-11-05):
does it reproduce the engine's sulfur lane at the shipped configuration, is its gradient right, and what do one solve
and one fit Jacobian cost against the engine at equal tolerance?

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml, which carries environment.yml's package set verbatim
plus JAX), with the ENGINE IN THE SAME PROCESS: every compared engine number is computed live (``integrate_sulfur``
on the same pot, the same node grid, the same parameters), back to back with the JAX run, so the parity is against
the code on disk and the timings are paired under one machine load. Two cores at most, niced:

    ./scripts/docker_maillard.sh run "set +u; conda activate maillard-jax; \
        XLA_FLAGS='--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2' OMP_NUM_THREADS=2 \
        nice -n 10 python scripts/generators/build_jax_sulfur_lane_module.py [--stage panel ...]"

Writes results/validation/jax_sulfur_lane_module.{json,md}. The helpers here (the panel pots, the engine's run on a
pot, the comparison and its tolerance) are what tests/scientific/test_jax_sulfur_parity.py imports.

THE PARITY TOLERANCE (from the solvers' tolerances, not from the observed difference). Both integrators run at the
engine's shipped setting, rtol 1e-8 / atol 1e-14 mmol/L (``engine._integrate_program``), and each keeps its LOCAL
error per step below rtol |y| + atol. For a dissipative mass-action network the global error is bounded by the sum of
the local errors (no growth), so after N accepted steps it is at most N (rtol |y| + atol) per solver. With N <= 1e4
steps per run (the engine's LSODA and the port's ESDIRK both take 1e2 - 1e4 on these pots), each side is within
1e-4 relative, the two together within 2e-4 relative = 8.7e-5 dex: the test bound is 1e-4 dex, on every species
whose level is above atol / rtol = 1e-6 mmol/L (where the relative term controls the step); below it the bound is
absolute, 2 N atol = 2e-10 mmol/L. The pH is -log10 of the H+ activity, a function of those levels: its bound is the
same 1e-4 (in pH units), the engine's own root tolerance (1e-9 pH) being negligible. The fixed point's sweep count
must be EQUAL (it is a discrete stopping rule; a different count is a different model, never inside a tolerance),
with the two recorded tie-breaks of ``Harness.matched`` where the engine's own discontinuities make the outcome a
matter of solver noise (a last residual at the fixed point's tolerance; a pot declared exactly at a table's jump).
At tight tolerance (rtol 1e-10 / atol 1e-16) the same argument gives the contract's 1e-6 dex (sec. 3), recorded here
on a subset.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUT = ROOT / "results" / "validation" / "jax_sulfur_lane_module"
FIT_TOL = (1e-8, 1e-14)          # engine._integrate_program's sulfur tolerance
TIGHT_TOL = (1e-10, 1e-16)
MAX_STEPS_BOUND = 1e4
DEX_TOL = 1e-4
PH_TOL = 1e-4
TIGHT_DEX_TOL = 1e-6
CELSIUS = 273.15
#: the random draws: the engine's own envelope draw of the sulfur lane (uncertainty.sulfur_joint_draw: the shipped
#: wave's Laplace covariance on the identified coordinates, flat across the admissible band on the others), seeds
DRAW_SEEDS = (11, 29, 47)
DRAW_POTS = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "mp_holdout_arp_glutathione_120C_60min_pH75_Zhou2023c")
TIGHT_POTS = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "thiamine_cys_glucose_120C_Bolton1994",
              "mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b")
#: Zhou's unbuffered ARP + glutathione pot: its pH drifts (ARP titrates), so the drift's coordinates matter there
GRADIENT_POT = "mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c"
GRADIENT_COORDINATES = ("arp_amine_pka", "acid_yield", "log10_k:k_arp_dpo", "ea:k_arp_dpo", "log10_k:k_nf_mft")
GRADIENT_SPECIES = ("MFT", "FFT", "NF", "H2S")
#: central-difference steps (the drift's acid yield ships at 7.5e-5, so its step is absolute and small)
FD_STEPS = {"arp_amine_pka": 1e-3, "acid_yield": 2e-5, "log10_k:k_arp_dpo": 1e-4, "ea:k_arp_dpo": 1e-2,
            "log10_k:k_nf_mft": 1e-4}
FIT_FD_STEPS = {"arp_amine_pka": 1e-2, "acid_yield": 2e-5, "log10_k:k_arp_dpo": 1e-3, "ea:k_arp_dpo": 1e-1,
                "log10_k:k_nf_mft": 1e-3}
#: synthetic pots for the modes no panel pot runs (not benchmarks)
PROFILE_POINTS = ((0.0, 25.0), (6.0, 145.0), (20.0, 145.0), (24.0, 100.0))
SYNTHETIC_INITIAL = {"PENT": 30.0, "Cys": 20.0, "THI": 2.0, "Glc": 20.0}


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def load() -> str:
    u = uptime()
    return u.split("load averages:")[-1].strip() if "load averages:" in u else u


# --------------------------------------------------------------------------- pots
def panel_pots() -> List[Dict[str, Any]]:
    """Every panel bundle the engine routes to the sulfur lane and answers (``panel.panel_bundles`` ->
    ``engine.declare_envelope``)."""
    from src.kinetic_core import engine, panel

    scored, _ = panel.panel_bundles()
    out = []
    for path, tag in scored:
        bench = panel.load_bundle(path)
        spec = panel.core_spec(bench)
        targets = panel.bundle_targets(bench)
        try:
            decl = engine.declare_envelope(spec, list(targets))
        except Exception:   # noqa: BLE001 - a refused declaration is not a sulfur-lane pot
            continue
        if engine.SULFUR not in engine._lanes_of(decl) or not decl.is_answerable:
            continue
        out.append({"id": str(bench.get("benchmark_id")), "panel": tag, "process": spec.process,
                    "initial": dict(decl.mapped_precursors), "targets": dict(targets)})
    return out


def synthetic_pots() -> List[Dict[str, Any]]:
    from src.kinetic_core import engine, ph_state

    phos = ph_state.BufferSpec(kind="phosphate", phosphate_mol_l=0.2)
    clamped = ph_state.BufferSpec(kind="clamped")
    staircase = engine.ThermalProgram(((5.0, 100.0), (15.0, 140.0)))
    return [
        {"id": "synthetic_profile_25_145_100", "process": engine.ProcessSpec(
            thermal=engine.ThermalProgram.from_profile(PROFILE_POINTS), ph=6.0, buffer=phos),
         "initial": dict(SYNTHETIC_INITIAL)},
        {"id": "synthetic_staircase_100_140", "process": engine.ProcessSpec(thermal=staircase, ph=5.5, buffer=phos),
         "initial": dict(SYNTHETIC_INITIAL)},
        {"id": "synthetic_clamped_pH6_130C", "process": engine.ProcessSpec(
            thermal=engine.ThermalProgram.isothermal(130.0, 30.0), ph=6.0, buffer=clamped),
         "initial": dict(SYNTHETIC_INITIAL)},
    ]


# --------------------------------------------------------------------------- the engine on a pot
def engine_run(lane, pot: Mapping[str, Any], parameters, drift, *, rtol: float, atol: float) -> Dict[str, Any]:
    """``engine._integrate_program``'s sulfur path on the pot's node grid (``integrate_sulfur`` per segment, the SID
    carried; a profile as one integration), at the given tolerance."""
    from src.kinetic_core import sulfur
    from src.kinetic_core_jax import sulfur_lane as sl

    process = pot["process"]
    state, _ = sl.engine_prelude(pot["initial"], process)
    spec = pot["jax"]["buffer_spec"]
    thermal = process.thermal
    t0 = time.perf_counter()
    states, ph, cool, sweeps, residual, sid = [], [], [], [], [], None
    if thermal.profile is not None:
        prof = (np.array([t for t, _ in thermal.profile]), np.array([c for _, c in thermal.profile]) + CELSIUS)
        segs = [(float(thermal.total_minutes), float(thermal.peak_temperature_c), prof)]
    else:
        segs = [(float(d), float(t), None) for d, t in thermal.segments]
    for (duration, t_c, prof), nodes in zip(segs, pot["jax"]["node_minutes"]):
        run = sulfur.integrate_sulfur(parameters, t_c + CELSIUS, state, nodes, ph=float(process.ph),
                                      ph_final=None if process.ph_final is None else float(process.ph_final),
                                      buffer_spec=spec, ph_drift=drift, rtol=rtol, atol=atol, sid_mmol_l=sid,
                                      temperature_profile=prof)
        if sid is None:
            sid = run.metadata.get("strong_ion_difference_mmol_L")
        states.append(run.concentrations)
        state = {k: float(run.concentrations[-1, i]) for i, k in enumerate(sulfur.SULFUR_STATE_KEYS)}
        md = run.metadata
        ph.append(None if run.ph_series is None else np.asarray(run.ph_series))
        cool.append(md.get("ph_final_cooled"))
        sweeps.append(md.get("ph_fixed_point_sweeps", 0))
        residual.append(md.get("ph_fixed_point_residual"))
    return {"states": np.stack(states), "ph": ph, "ph_final_cooled": cool, "sweeps": sweeps, "residual": residual,
            "sid": sid, "seconds": time.perf_counter() - t0}


def compare(jax_run, eng: Mapping[str, Any], *, rtol: float, atol: float) -> Dict[str, Any]:
    """Every species at every node of every segment: the max |dex| where either level is above atol / rtol, the max
    absolute difference below; the node pH and the final cooled pH; the sweep counts."""
    j = np.asarray(jax_run.ys)
    e = np.asarray(eng["states"])
    floor = atol / rtol
    above = np.maximum(j, e) >= floor
    dex = np.abs(np.log10(np.maximum(j, 1e-300)) - np.log10(np.maximum(e, 1e-300)))
    dex_above = np.where(above, dex, 0.0)
    worst = np.unravel_index(int(np.argmax(dex_above)), dex_above.shape)
    out = {"max_dex": float(dex_above.max()), "worst_species": None,
           "max_abs_below_floor": float(np.max(np.where(above, 0.0, np.abs(j - e)))),
           "floor_mmol_L": floor, "jax_sweeps": [int(s) for s in np.asarray(jax_run.sweeps)],
           "engine_sweeps": [int(s) for s in eng["sweeps"]]}
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS

    out["worst_species"] = SULFUR_STATE_KEYS[worst[2]] if out["max_dex"] > 0 else None
    jp = np.asarray(jax_run.ph)
    ph_diffs = [float(np.max(np.abs(jp[s] - eng["ph"][s]))) for s in range(len(eng["ph"])) if eng["ph"][s] is not None]
    out["max_abs_node_ph"] = max(ph_diffs) if ph_diffs else 0.0
    cool = [abs(float(np.asarray(jax_run.ph_cool)[s, -1]) - float(c)) for s, c in enumerate(eng["ph_final_cooled"])
            if c is not None]
    out["max_abs_cooled_ph"] = max(cool) if cool else 0.0
    return out


def table_jumps() -> Tuple[float, ...]:
    """The label pH values where a measured table JUMPS (not merely kinks): Zheng & Ho's mechanism boundary (the
    cysteine and glutathione sets switch to the pH-3 pair below it) and the thiamine share's floor."""
    from src.kinetic_core_jax import sulfur_lane as sl

    c = sl._term_constants()
    return (float(c["zheng_boundary"]), float(c["thi_share_floor"]))


def tie_tolerance(rtol: float) -> float:
    """How far apart the two sides' last fixed-point residuals may be by solver noise alone: the residual is a max of
    |fresh - node pH| over the nodes, each fresh pH within the pH bound of the other side's, so 2 x PH_TOL at the fit
    tolerance (scaled with rtol)."""
    return 2.0 * PH_TOL * rtol / FIT_TOL[0]


def within(c: Mapping[str, Any], *, rtol: float, atol: float, dex_tol: float = DEX_TOL,
           ph_tol: float = PH_TOL) -> List[str]:
    """The tolerance checks of ``compare``'s result (empty = pass)."""
    bad = []
    if c["jax_sweeps"] != c["engine_sweeps"]:
        bad.append(f"sweeps {c['jax_sweeps']} vs the engine's {c['engine_sweeps']}")
    if c["max_dex"] > dex_tol:
        bad.append(f"max dex {c['max_dex']:.3g} ({c['worst_species']}) > {dex_tol:g}")
    if c["max_abs_below_floor"] > 2.0 * MAX_STEPS_BOUND * atol:
        bad.append(f"below-floor abs {c['max_abs_below_floor']:.3g} > {2.0 * MAX_STEPS_BOUND * atol:g}")
    if c["max_abs_node_ph"] > ph_tol or c["max_abs_cooled_ph"] > ph_tol:
        bad.append(f"pH {c['max_abs_node_ph']:.3g} / cooled {c['max_abs_cooled_ph']:.3g} > {ph_tol:g}")
    return bad


def draw(seed: int):
    """One engine envelope draw of the sulfur lane: (parameters, PhDrift), as ``uncertainty`` composes it."""
    from src.kinetic_core import engine, uncertainty

    d = uncertainty.sulfur_joint_draw(np.random.default_rng(seed))
    frozen = engine.frozen_parameters(engine.SULFUR)
    blocks: Dict[str, Any] = {}
    for block, value in d["maillard_blocks"].items():
        if isinstance(value, dict):
            merged = dict(frozen[block])
            merged.update(value)
            blocks[block] = merged
        else:
            blocks[block] = value
    for block in ("log10_k_ref_at_145C", "lumped_formation_Ea_kJ_mol", "decay_Ea_kJ_mol",
                  "formation_Ea_by_route_kJ_mol", "oxygen", "oxygen_log10_k"):
        if block in frozen:
            blocks.setdefault(block, frozen[block])
    return engine.core_parameters(engine.SULFUR, frozen=blocks), d["ph_drift"] or engine.core_ph_drift(), d["coords"]


class Harness:
    """One lane, its compiled solves (cached per pot shape, mode and tolerance) and the pots."""

    def __init__(self, cfg=None, solver: str = "Kvaerno5", lane=None):
        from src.kinetic_core_jax import sulfur_lane as sl

        self.sl = sl
        self.lane = sl.build(cfg, solver=solver) if lane is None else lane
        self._compiled: Dict[Tuple, Any] = {}

    def prepare(self, pot: Dict[str, Any]) -> Dict[str, Any]:
        pot["jax"] = self.sl.make_pot(self.lane, pot["initial"], pot["process"])
        return pot

    def solve(self, pot, theta=None, *, rtol=FIT_TOL[0], atol=FIT_TOL[1], solver=None, force_sweeps=None):
        import jax

        jp = pot["jax"]
        key = (jp["mode"], np.asarray(jp["node_s"]).shape, np.asarray(jp["knot_s"]).shape, rtol, atol,
               type(solver).__name__ if solver is not None else None, force_sweeps)
        if key not in self._compiled:
            kw = {"rtol": rtol, "atol": atol, "force_sweeps": force_sweeps}
            if solver is not None:
                kw["solver"] = solver
            self._compiled[key] = self.sl.compiled(self.lane, jp, **kw)
        f = self._compiled[key]
        theta = self.lane.theta0() if theta is None else theta
        t0 = time.perf_counter()
        run = f(theta, self.sl.pot_arrays(jp))
        jax.block_until_ready(run.ys)
        return run, time.perf_counter() - t0

    def jump_nudge(self, pot, run, eng, drift=None):
        """A pot declared at a label pH that IS a jump of a measured pH table (``Module.breakpoints["ph_label"]``:
        Zheng & Ho's mechanism boundary at 5.0, ...) reads that table at t = 0 at the cooled pH of its initial state,
        which equals the declared label up to the root-finders' last digits: which side of the jump the first node
        reads is then noise, on both sides. If the two sides fall on opposite sides, returns the (S, N) nudge that puts
        the port's first-node table read at the engine's value (a difference below 1e-7 pH); None otherwise."""
        from src.kinetic_core import engine, ph_state

        label = float(pot["process"].ph)
        near = [b for b in table_jumps() if abs(label - b) < 1e-9]
        if not near or eng.get("sid") is None:
            return None
        y0 = pot["jax"]["y0"]
        state = {k: float(y0[i]) for i, k in enumerate(self.lane.species)}
        e_c = ph_state.ph_of_state(state, ph_state.LABEL_TEMPERATURE_K, float(eng["sid"]),
                                   drift if drift is not None else engine.core_ph_drift(), pot["jax"]["buffer_spec"])
        j_c = float(np.asarray(run.ph_cool)[0, 0])
        b = near[0]
        if (e_c < b) == (j_c < b) or abs(e_c - j_c) > 1e-7:
            return None
        nudge = np.zeros_like(np.asarray(pot["jax"]["label_nudge"]))
        nudge[0, 0] = e_c - j_c
        return nudge

    def matched(self, pot, run, eng, theta=None, *, drift=None, rtol=FIT_TOL[0], atol=FIT_TOL[1]):
        """``compare`` of a JAX run with the engine's, with the TWO exceptions the engine's discontinuities force,
        both recorded: (i) a first node ON a table's jump (``jump_nudge``) is read on the engine's side; (ii) when the
        sweep counts differ on a one-segment pot and the two last residuals straddle the fixed point's tolerance
        within ``tie_tolerance`` (solver noise alone decided one sweep more or less), the JAX run is redone at the
        engine's count. Either way it is then compared at the same tolerance (``jump_side_matched``,
        ``tie_broken``)."""
        nudge = self.jump_nudge(pot, run, eng, drift)
        jump = None
        if nudge is not None:
            jump = {"engine_first_node_cooled_ph": float(np.asarray(run.ph_cool)[0, 0] + nudge[0, 0]),
                    "port_first_node_cooled_ph": float(np.asarray(run.ph_cool)[0, 0]),
                    "unnudged_max_dex": compare(run, eng, rtol=rtol, atol=atol)["max_dex"]}
            pot = dict(pot, jax=dict(pot["jax"], label_nudge=nudge))
            run, _ = self.solve(pot, theta, rtol=rtol, atol=atol)
        c = compare(run, eng, rtol=rtol, atol=atol)
        c["jump_side_matched"] = jump
        c["tie_broken"] = False
        if c["jax_sweeps"] == c["engine_sweeps"] or len(c["engine_sweeps"]) != 1:
            return run, c
        n_e = int(c["engine_sweeps"][0])
        forced, _ = self.solve(pot, theta, rtol=rtol, atol=atol, force_sweeps=n_e)
        tol = self.lane.settings.tolerance
        r_j, r_e = float(np.asarray(forced.residual)[0]), float(eng["residual"][0])
        # the engine stopped first: at its count its residual is below the tolerance and the port's is not
        straddle = n_e < int(c["jax_sweeps"][0]) and r_e < tol <= r_j
        if straddle and (r_j - r_e) <= tie_tolerance(rtol):
            c = compare(forced, eng, rtol=rtol, atol=atol)
            c.update({"tie_broken": True, "unforced_jax_sweeps": [int(s) for s in np.asarray(run.sweeps)],
                      "residual_jax_at_engine_count": r_j, "residual_engine": r_e, "jump_side_matched": jump})
            return forced, c
        c.update({"residual_jax_at_engine_count": r_j, "residual_engine": r_e})
        return run, c


# --------------------------------------------------------------------------- stages
def stage_structure(h: Harness) -> Dict[str, Any]:
    from src.kinetic_core import sulfur
    from src.kinetic_core_jax import sulfur_structure as ss

    net = h.lane.network
    found = ss.classify()
    return {"reactions": len(net.reactions), "species": len(net.species),
            "reactions_equal_engine_order": tuple(net.reactions) == tuple(sulfur.FULL_REACTION_KEYS),
            "stoichiometry_equal": bool(np.array_equal(net.stoichiometry, sulfur.SULFUR_STOICHIOMETRY)),
            "coordinates": len(h.lane.module.params), "hand_ported": sorted(ss.HAND_PORTED), "off": list(ss.OFF),
            "engine_custom": sorted(found["custom"]), "problems": ss.structure_problems(),
            "linear_coordinates": [k for k in h.lane.module.params if k.startswith("k_ref:")],
            "titratable_groups": list(h.lane.balance.names)}


RATE_T = (25.0, 80.0, 100.0, 120.0, 145.0, 170.0, 200.0)
RATE_PH = (3.0, 4.4, 5.0, 6.2, 7.5, 9.0)
RATE_LABEL_SHIFT = (0.0, 0.45)
RATE_I = (0.0, 0.12, 0.9)


def rate_probe(h: Harness, parameters, drift) -> Dict[str, Any]:
    import jax.numpy as jnp
    from src.kinetic_core import sulfur
    from src.kinetic_core_jax.network import RateCtx

    theta = h.sl.theta_from(h.lane, parameters, drift)
    f = jax_rates(h)
    worst, n, zeros_bad = 0.0, 0, 0
    for t_c in RATE_T:
        for ph in RATE_PH:
            for shift in RATE_LABEL_SHIFT:
                for i in RATE_I:
                    e = sulfur.sulfur_rate_constants_at(parameters, t_c + CELSIUS, ph, ph_label=ph + shift,
                                                        ionic_strength_m=i)
                    e = np.array([e[r] for r in h.lane.network.reactions])
                    k = np.asarray(f(theta, jnp.asarray([t_c, ph, ph + shift, i])))
                    m = e != 0.0
                    zeros_bad += int(np.sum(k[~m] != 0.0))
                    worst = max(worst, float(np.max(np.abs(k[m] - e[m]) / np.abs(e[m]))))
                    n += len(e)
    return {"values": n, "worst_relative": worst, "nonzero_where_engine_zero": zeros_bad}


_RATES = {}


def jax_rates(h: Harness):
    import jax
    from src.kinetic_core_jax.network import RateCtx

    if id(h) not in _RATES:
        _RATES[id(h)] = jax.jit(lambda th, x: h.lane.network.rate_constants(
            th, RateCtx(t_c=x[0], t_k=x[0] + CELSIUS, ph=x[1], ph_label=x[2], ionic_m=x[3])))
    return _RATES[id(h)]


def stage_rates(h: Harness) -> Dict[str, Any]:
    from src.kinetic_core import engine

    out = {"shipped": rate_probe(h, engine.core_parameters(engine.SULFUR), engine.core_ph_drift())}
    for seed in DRAW_SEEDS[:2]:
        p, d, _ = draw(seed)
        out[f"draw_{seed}"] = rate_probe(h, p, d)
    return out


def stage_balance(h: Harness) -> Dict[str, Any]:
    """SID, in-situ (pH, I) and the cooled pH of states along a run, against ph_state."""
    import jax.numpy as jnp
    from src.kinetic_core import engine, ph_state

    pots = [h.prepare(p) for p in panel_pots() if p["id"] in DRAW_POTS]
    drift = engine.core_ph_drift()
    from src.kinetic_core_jax import ph_coupling

    bal = ph_coupling.solver_of(h.lane.network)     # the network's own solver (2026-11-06; was the lane's)
    th = h.lane.theta0()
    worst = {"sid": 0.0, "ph": 0.0, "ionic_rel": 0.0, "cool": 0.0}
    n = 0
    for pot in pots:
        run, _ = h.solve(pot)
        spec = pot["jax"]["buffer_spec"]
        buf = jnp.asarray(pot["jax"]["buffer"])
        y0 = pot["jax"]["y0"]
        init = {k: float(y0[i]) for i, k in enumerate(h.lane.species)}
        sid_e = ph_state.strong_ion_difference(float(pot["process"].ph), init, drift, spec)
        sid_j = float(bal.strong_ion_difference(float(pot["process"].ph), jnp.asarray(y0), th, buf))
        worst["sid"] = max(worst["sid"], abs(sid_j - sid_e) / max(1.0, abs(sid_e)))
        for row in np.asarray(run.ys[0])[::8]:
            m = {k: float(row[i]) for i, k in enumerate(h.lane.species)}
            for t_k in (298.15, 373.15, 418.15):
                pe, ie = ph_state.ph_and_ionic_strength_of_state(m, t_k, sid_e, drift, spec)
                pj, ij, cj = bal.state(jnp.asarray(row), t_k, sid_e, th, buf)
                ce = ph_state.ph_of_state(m, ph_state.LABEL_TEMPERATURE_K, sid_e, drift, spec)
                worst["ph"] = max(worst["ph"], abs(float(pj) - pe))
                worst["ionic_rel"] = max(worst["ionic_rel"], abs(float(ij) - ie) / max(ie, 1e-12))
                worst["cool"] = max(worst["cool"], abs(float(cj) - ce))
                n += 1
    return {"states": n, **worst, "engine_root_tolerance_ph": 1e-9}


def stage_panel(h: Harness, *, timing: bool = True) -> Dict[str, Any]:
    import diffrax
    from src.kinetic_core import engine

    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    rows = []
    for pot in [h.prepare(p) for p in panel_pots()]:
        run, t_compile = h.solve(pot)
        eng = engine_run(h.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1])
        run, c = h.matched(pot, run, eng)
        rec = {"id": pot["id"], "panel": pot["panel"], "initial": pot["initial"], **c,
               "pass": not within(c, rtol=FIT_TOL[0], atol=FIT_TOL[1]), "engine_s": eng["seconds"],
               "jax_steps": int(np.sum(np.asarray(run.num_steps))), "load": load()}
        # the engine's own path (the prelude included) on the final state
        lane_final, _ = engine._integrate_program(engine.SULFUR, params, dict(pot["initial"]), pot["process"])
        ef = np.array([float(lane_final.get(k, 0.0)) for k in h.lane.species])
        jf = np.asarray(run.ys)[-1, -1]
        m = np.maximum(ef, jf) >= FIT_TOL[1] / FIT_TOL[0]
        rec["final_vs_integrate_program_dex"] = float(np.max(np.abs(np.log10(np.maximum(jf[m], 1e-300))
                                                                    - np.log10(ef[m]))))
        if timing:
            js, es = [], []
            for _ in range(2):
                js.append(h.solve(pot)[1])
                es.append(engine_run(h.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1])["seconds"])
            run4, _ = h.solve(pot, solver=diffrax.Kvaerno4())
            k4 = [h.solve(pot, solver=diffrax.Kvaerno4())[1] for _ in range(2)]
            c4 = compare(run4, eng, rtol=FIT_TOL[0], atol=FIT_TOL[1])
            rec.update({"jax_k5_s": min(js), "engine_s": min(es + [eng["seconds"]]), "jax_k4_s": min(k4),
                        "jax_k4_steps": int(np.sum(np.asarray(run4.num_steps))), "k4_max_dex": c4["max_dex"],
                        "k4_pass": not within(c4, rtol=FIT_TOL[0], atol=FIT_TOL[1]), "first_call_s": t_compile})
        rows.append(rec)
        print(f"  {pot['id']}: {c['max_dex']:.2e} dex, sweeps {c['jax_sweeps']}/{c['engine_sweeps']}, "
              f"jax {rec.get('jax_k5_s', 0):.2f} s, engine {rec['engine_s']:.2f} s", flush=True)
    return {"tolerance": {"rtol": FIT_TOL[0], "atol": FIT_TOL[1], "dex": DEX_TOL, "ph": PH_TOL}, "pots": rows}


def stage_draws(h: Harness) -> Dict[str, Any]:
    pots = [h.prepare(p) for p in panel_pots() if p["id"] in DRAW_POTS]
    rows = []
    for seed in DRAW_SEEDS:
        params, drift, coords = draw(seed)
        theta = h.sl.theta_from(h.lane, params, drift)
        for pot in pots:
            run, _ = h.solve(pot, theta)
            eng = engine_run(h.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1])
            run, c = h.matched(pot, run, eng, theta, drift=drift)
            rows.append({"seed": seed, "id": pot["id"], "drawn_coordinates": len(coords), **c,
                         "pass": not within(c, rtol=FIT_TOL[0], atol=FIT_TOL[1])})
            print(f"  draw {seed} {pot['id']}: {c['max_dex']:.2e} dex", flush=True)
    return {"seeds": list(DRAW_SEEDS), "rows": rows}


def stage_modes(h: Harness) -> Dict[str, Any]:
    """The modes no panel pot runs: a profile, a staircase, a clamped pH; and the inline coupling on two panel pots."""
    from src.kinetic_core import config, engine

    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    rows = []
    for pot in [h.prepare(p) for p in synthetic_pots()]:
        run, _ = h.solve(pot)
        eng = engine_run(h.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1])
        run, c = h.matched(pot, run, eng)
        rows.append({"id": pot["id"], "coupling": "fixed_point_nodes", **c,
                     "pass": not within(c, rtol=FIT_TOL[0], atol=FIT_TOL[1])})
        print(f"  {pot['id']}: {c['max_dex']:.2e} dex", flush=True)
    cfg = config.Config.from_globals().replace({"ph.coupling": "inline"})
    inline = Harness(cfg)
    for pot in [inline.prepare(p) for p in panel_pots() if p["id"] in DRAW_POTS] + \
            [inline.prepare(p) for p in synthetic_pots() if p["id"] == "synthetic_profile_25_145_100"]:
        run, _ = inline.solve(pot)
        with config.applied(cfg):
            eng = engine_run(inline.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1])
        c = compare(run, eng, rtol=FIT_TOL[0], atol=FIT_TOL[1])
        c["jax_sweeps"], c["engine_sweeps"] = [0], [0]
        rows.append({"id": pot["id"], "coupling": "inline", **c,
                     "pass": not within(c, rtol=FIT_TOL[0], atol=FIT_TOL[1])})
        print(f"  inline {pot['id']}: {c['max_dex']:.2e} dex", flush=True)
    return {"rows": rows}


def stage_tight(h: Harness) -> Dict[str, Any]:
    from src.kinetic_core import engine

    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    rows = []
    for pot in [h.prepare(p) for p in panel_pots() if p["id"] in TIGHT_POTS]:
        run, t = h.solve(pot, rtol=TIGHT_TOL[0], atol=TIGHT_TOL[1])
        eng = engine_run(h.lane, pot, params, drift, rtol=TIGHT_TOL[0], atol=TIGHT_TOL[1])
        run, c = h.matched(pot, run, eng, rtol=TIGHT_TOL[0], atol=TIGHT_TOL[1])
        rows.append({"id": pot["id"], **c, "pass": not within(c, rtol=TIGHT_TOL[0], atol=TIGHT_TOL[1],
                                                              dex_tol=TIGHT_DEX_TOL, ph_tol=TIGHT_DEX_TOL),
                     "jax_s": t, "engine_s": eng["seconds"]})
        print(f"  tight {pot['id']}: {c['max_dex']:.2e} dex", flush=True)
    return {"tolerance": {"rtol": TIGHT_TOL[0], "atol": TIGHT_TOL[1], "dex": TIGHT_DEX_TOL}, "rows": rows}


def gradient_check(h: Harness, pot, coordinates=GRADIENT_COORDINATES, *, rtol=TIGHT_TOL[0], atol=TIGHT_TOL[1],
                   steps=None) -> Dict[str, Any]:
    """jax.grad (reverse mode: the fixed point as a scan) of one scalar of the run against the port's own central
    differences, and against forward mode."""
    import diffrax
    import jax
    import jax.numpy as jnp

    sl, lane = h.sl, h.lane
    theta = lane.theta0()
    arrays = sl.pot_arrays(pot["jax"])
    mode = pot["jax"]["mode"]
    idx = [lane.index[s] for s in GRADIENT_SPECIES]

    def scalar(x, adjoint, reverse):
        th = {sl.NAME: dict(theta[sl.NAME], **{c: x[i] for i, c in enumerate(coordinates)})}
        run = sl.integrate(lane, th, dict(arrays, mode=mode), adjoint=adjoint, reverse=reverse, rtol=rtol, atol=atol)
        ys = run.ys[-1]
        logs = jnp.log10(jnp.maximum(ys[-1, jnp.asarray(idx)], 1e-300))
        return jnp.sum(logs * jnp.asarray([1.0, 0.5, 0.25, 0.125])) + run.ph_cool[-1, -1]

    x0 = jnp.stack([theta[sl.NAME][c] for c in coordinates])
    t0 = time.perf_counter()
    grad = jax.jit(jax.grad(lambda x: scalar(x, diffrax.RecursiveCheckpointAdjoint(), True)))
    g = np.asarray(grad(x0))
    t_grad_first = time.perf_counter() - t0
    t0 = time.perf_counter()
    g = np.asarray(grad(x0))
    t_grad = time.perf_counter() - t0
    fwd = jax.jit(jax.jacfwd(lambda x: scalar(x, diffrax.ForwardMode(), False)))
    jf = np.asarray(fwd(x0))
    f = jax.jit(lambda x: scalar(x, diffrax.RecursiveCheckpointAdjoint(), False))
    steps = steps or {c: FD_STEPS[c] for c in coordinates}
    fd = []
    for i, c in enumerate(coordinates):
        e = np.zeros(len(coordinates))
        e[i] = steps[c]
        fd.append((float(f(x0 + e)) - float(f(x0 - e))) / (2.0 * steps[c]))
    fd = np.asarray(fd)
    scale = np.maximum(np.abs(fd), 1e-8)
    return {"pot": pot["id"], "coordinates": list(coordinates), "steps": [steps[c] for c in coordinates],
            "grad": g.tolist(), "jacfwd": jf.tolist(), "central_fd": fd.tolist(),
            "rel_grad_vs_fd": (np.abs(g - fd) / scale).tolist(), "rel_grad_vs_jacfwd": (np.abs(g - jf) / scale).tolist(),
            "finite": bool(np.all(np.isfinite(g))), "grad_first_call_s": t_grad_first, "grad_s": t_grad,
            "rtol": rtol, "atol": atol}


def stage_gradient(h: Harness) -> Dict[str, Any]:
    pot = next(h.prepare(p) for p in panel_pots() if p["id"] == GRADIENT_POT)
    return gradient_check(h, pot)


def stage_fit_jacobian(h: Harness) -> Dict[str, Any]:
    """One fit Jacobian (17 observables: log10 MFT / FFT / Cys / NF at four nodes and the cooled pH; the five gradient
    coordinates) at the fit tolerance: JAX forward mode against the engine's forward and central differences (the
    engine's runs timed, the same pot, the same tolerance)."""
    import dataclasses

    import jax
    from src.kinetic_core import engine

    pot = next(h.prepare(p) for p in panel_pots() if p["id"] == GRADIENT_POT)
    n = np.asarray(pot["jax"]["node_s"]).shape[1]
    rows = [n // 8, n // 3, (2 * n) // 3, n - 1]
    jac = h.sl.jacobian(h.lane, pot["jax"], GRADIENT_COORDINATES, GRADIENT_SPECIES, rows,
                        rtol=FIT_TOL[0], atol=FIT_TOL[1])
    theta, arrays = h.lane.theta0(), h.sl.pot_arrays(pot["jax"])
    t0 = time.perf_counter()
    value, jm = jac(theta, arrays)
    jax.block_until_ready(jm)
    first = time.perf_counter() - t0
    times = []
    for _ in range(2):
        t0 = time.perf_counter()
        value, jm = jac(theta, arrays)
        jax.block_until_ready(jm)
        times.append(time.perf_counter() - t0)
    solve_s = min(h.solve(pot)[1] for _ in range(2))
    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()

    def moved(name, hstep):
        kind, _, key = name.partition(":")
        p, d = dict(params), drift
        if kind == "log10_k":
            p[key] = dataclasses.replace(p[key], k_ref=float(p[key].k_ref) * 10.0 ** hstep)
        elif kind == "ea":
            p[key] = dataclasses.replace(p[key], ea_kj_mol=float(p[key].ea_kj_mol) + hstep)
        elif name == "acid_yield":
            d = dataclasses.replace(d, acid_yield=d.acid_yield + hstep)
        else:
            d = dataclasses.replace(d, arp_amine_pka=d.arp_amine_pka + hstep)
        return p, d

    def obs(eng):
        s = eng["states"][-1]
        vals = [math.log10(max(s[r, h.lane.index[sp]], 1e-300)) for sp in GRADIENT_SPECIES for r in rows]
        return np.array(vals + [eng["ph_final_cooled"][-1]])

    steps = dict(FIT_FD_STEPS)
    t0 = time.perf_counter()
    base = obs(engine_run(h.lane, pot, params, drift, rtol=FIT_TOL[0], atol=FIT_TOL[1]))
    fwd_cols = []
    for c in GRADIENT_COORDINATES:
        p, d = moved(c, steps[c])
        fwd_cols.append((obs(engine_run(h.lane, pot, p, d, rtol=FIT_TOL[0], atol=FIT_TOL[1])) - base) / steps[c])
    t_fwd = time.perf_counter() - t0
    t0 = time.perf_counter()
    cen_cols = []
    for c in GRADIENT_COORDINATES:
        up = obs(engine_run(h.lane, pot, *moved(c, steps[c]), rtol=FIT_TOL[0], atol=FIT_TOL[1]))
        dn = obs(engine_run(h.lane, pot, *moved(c, -steps[c]), rtol=FIT_TOL[0], atol=FIT_TOL[1]))
        cen_cols.append((up - dn) / (2.0 * steps[c]))
    t_cen = time.perf_counter() - t0
    jm = np.asarray(jm)
    cen = np.array(cen_cols).T
    fwd = np.array(fwd_cols).T
    col_scale = np.maximum(np.max(np.abs(cen), axis=0), 1e-12)
    return {"pot": pot["id"], "observables": len(value), "coordinates": list(GRADIENT_COORDINATES),
            "jax_jacobian_s": min(times), "jax_jacobian_first_call_s": first, "jax_solve_s": solve_s,
            "engine_forward_fd_s": t_fwd, "engine_central_fd_s": t_cen,
            "engine_runs_forward": 1 + len(GRADIENT_COORDINATES), "engine_runs_central": 2 * len(GRADIENT_COORDINATES),
            "max_col_rel_jax_vs_engine_central": float(np.max(np.max(np.abs(jm - cen), axis=0) / col_scale)),
            "max_col_rel_engine_forward_vs_central": float(np.max(np.max(np.abs(fwd - cen), axis=0) / col_scale)),
            "load": load()}


STAGES = ("structure", "rates", "balance", "panel", "draws", "modes", "tight", "gradient", "fit_jacobian")


# --------------------------------------------------------------------------- the record
def _fmt(x, p=2):
    return "n/a" if x is None else f"{x:.{p}e}" if isinstance(x, float) and (abs(x) < 1e-2 or abs(x) >= 1e4) and x != 0 \
        else f"{x:.{p}f}" if isinstance(x, float) else str(x)


def exceptions(rows) -> List[str]:
    """One line per row where ``Harness.matched`` used a tie-break."""
    out = []
    for p in rows:
        tag = f"{p['id']}" + (f" (draw {p['seed']})" if "seed" in p else "")
        j = p.get("jump_side_matched")
        if j:
            out.append(f"- {tag}: declared at a table jump; the first node's cooled pH is {j['port_first_node_cooled_ph']:.15f} "
                       f"in the port and {j['engine_first_node_cooled_ph']:.15f} in the engine, opposite sides of the "
                       f"jump (unmatched: {j['unnudged_max_dex']:.2e} dex); read on the engine's side.")
        if p.get("tie_broken"):
            out.append(f"- {tag}: unforced, the port ran {p['unforced_jax_sweeps'][0]} sweeps and the engine "
                       f"{p['engine_sweeps'][0]}; at sweep {p['engine_sweeps'][0]} the residuals are "
                       f"{p['residual_jax_at_engine_count']:.6f} (port) and {p['residual_engine']:.6f} (engine), "
                       "either side of the 1e-3 tolerance; rerun at the engine's count.")
    return out


def markdown(r: Mapping[str, Any]) -> str:
    import diffrax
    import jax

    s, panel = r["structure"], r["panel"]
    pots = panel["pots"]
    worst = max(pots, key=lambda p: p["max_dex"])
    L = [
        "# The sulfur lane as an A2 contract module (JAX + diffrax)",
        "",
        f"*Generated by `scripts/generators/build_jax_sulfur_lane_module.py` in the separate `maillard-jax` env (jax "
        f"{jax.__version__}, diffrax {diffrax.__version__}, float64, CPU), the engine in the same process (computed "
        "live on the same pots, node grids and parameters; paired timings). Plan step A2 + H3, 2026-11-05. "
        "Module: `src/kinetic_core_jax/sulfur_lane.py`; structure registry and drift check: "
        "`src/kinetic_core_jax/sulfur_structure.py`. Nothing shipped moves; `src/kinetic_core` does not import it.*",
        "",
        "## What was ported, and how",
        "",
        f"- ONE module, `sulfur`: {s['reactions']} reactions in `sulfur.FULL_REACTIONS` order (identical: "
        f"{s['reactions_equal_engine_order']}), {s['species']} species, the stoichiometry matrix equal to the engine's "
        f"({s['stoichiometry_equal']}), {s['coordinates']} coordinates (log10 k_ref or, for the "
        f"{len(s['linear_coordinates'])} keys the engine carries at exactly zero, a linear k_ref; Ea where a key has "
        "one; the drift's acid yield and ARP pK).",
        "- DERIVED from the engine at build time: species, elements, reactions, stoichiometry, mass-action orders, "
        "parameter values / Ea / T_ref / caps, the pH-factor kind of every step, every table the hand-ported terms "
        "read, the titratable groups and buffer groups, the pot prelude's constants and `integrate_sulfur`'s node "
        "and sweep defaults.",
        f"- HAND-PORTED ({len(s['hand_ported'])} reactions; equal to the set the engine computes outside the generic "
        f"rule, found numerically: {sorted(s['hand_ported']) == s['engine_custom']}): Zheng & Ho's cysteine and "
        "glutathione thermolysis, the thiamine H2S share, the thioether release, the Amadori pin, and "
        f"{len(s['off'])} steps held at zero because their routes are off at the shipped configuration. Plus the "
        "generic rate rule with the six pH factors, the charge balance and its (pH, I) root, and the pH coupling.",
        f"- The structure check (`sulfur_structure.structure_problems`, core env, CI unit tier): "
        f"{'clean' if not s['problems'] else s['problems']}.",
        f"- The charge balance's groups: {', '.join(s['titratable_groups'])}.",
        "",
        "## Rates and the charge balance (no integration)",
        "",
    ]
    for label, rp in r["rates"].items():
        L.append(f"- Rate constants, {label}: {rp['values']} values ((T, pH, label pH, I) grid x 132 reactions), worst "
                 f"relative difference {rp['worst_relative']:.1e}; non-zero where the engine is zero: "
                 f"{rp['nonzero_where_engine_zero']}.")
    b = r["balance"]
    L += [f"- Charge balance on {b['states']} (state, T) points of two panel runs: max |pH| difference "
          f"{b['ph']:.1e} (in situ) / {b['cool']:.1e} (cooled), I {b['ionic_rel']:.1e} relative, SID "
          f"{b['sid']:.1e} relative (the engine's bisection tolerance is 1e-9 pH).",
          "",
          "## Parity on every panel pot the sulfur lane answers (shipped configuration: pH coupling fixed_point_nodes)",
          "",
          f"Tolerance (derived from the solvers' tolerances, module docstring of the generator): rtol "
          f"{FIT_TOL[0]:g} / atol {FIT_TOL[1]:g} on both sides (the engine's shipped setting); |dex| <= {DEX_TOL:g} on "
          f"every species above atol / rtol = {FIT_TOL[1] / FIT_TOL[0]:g} mmol/L at every pH node (~81 per segment), "
          f"absolute <= {2 * MAX_STEPS_BOUND * FIT_TOL[1]:g} mmol/L below; pH (node and final cooled) <= {PH_TOL:g}; "
          "the fixed point's sweep count equal. One exception, forced by the engine's discrete stopping rule: where "
          "the engine stops at a sweep whose residual is below the fixed point's tolerance (1e-3 pH) and the port's "
          f"residual at that sweep is above it by less than {tie_tolerance(FIT_TOL[0]):g} pH (twice the pH bound: "
          "solver noise alone decided), the port is rerun at the engine's sweep count and compared at the same "
          "tolerance (marked * below, with both residuals). A second, for pots DECLARED at a table jump (label pH 5.0 is "
          "Zheng & Ho's mechanism boundary): the first node's cooled pH equals the declared label up to the "
          "root-finders' last digits, so which side of the jump it reads is noise on both sides; where they differ the "
          "port reads the engine's side (marked j).",
          "",
          f"**{sum(p['pass'] for p in pots)} / {len(pots)} pots pass.** Largest difference: {worst['max_dex']:.2e} dex "
          f"({worst['worst_species']}, {worst['id']}).",
          ""]
    L += exceptions(pots)
    L += ["",
          "| pot | sweeps (JAX / engine) | max dex | worst species | max pH diff (node / cooled) | final vs `_integrate_program` dex | JAX K5 s (steps) | JAX K4 s (steps; max dex) | engine s | load |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for p in pots:
        mark = ("*" if p.get("tie_broken") else "") + ("j" if p.get("jump_side_matched") else "")
        L.append(f"| {p['id']} | {p['jax_sweeps'][0]}{mark} / {p['engine_sweeps'][0]}"
                 f" | {p['max_dex']:.1e} | "
                 f"{p['worst_species']} | {p['max_abs_node_ph']:.1e} / {p['max_abs_cooled_ph']:.1e} | "
                 f"{p['final_vs_integrate_program_dex']:.1e} | {p.get('jax_k5_s', float('nan')):.2f} ({p['jax_steps']}) | "
                 f"{p.get('jax_k4_s', float('nan')):.2f} ({p.get('jax_k4_steps', '')}; {p.get('k4_max_dex', float('nan')):.1e}) | "
                 f"{p['engine_s']:.2f} | {p['load']} |")
    js = [p["jax_k5_s"] for p in pots if "jax_k5_s" in p]
    es = [p["engine_s"] for p in pots if "jax_k5_s" in p]
    k4 = [p["jax_k4_s"] for p in pots if "jax_k4_s" in p]
    if js:
        ratio = np.array(js) / np.array(es)
        L += ["", f"Per-pot solve time, JAX Kvaerno5 / engine LSODA at equal tolerance: median {np.median(ratio):.2f}x "
                  f"(range {ratio.min():.2f} - {ratio.max():.2f}); panel total {sum(js):.1f} s vs {sum(es):.1f} s; "
                  f"Kvaerno4 total {sum(k4):.1f} s. Paired back to back on a shared machine (load beside each row): "
                  "PROVISIONAL as absolute numbers, the ratio is the reading."]
    d = r["draws"]
    L += ["", "## Random parameter draws (the engine's own envelope draw, `uncertainty.sulfur_joint_draw`)", "",
          f"{sum(x['pass'] for x in d['rows'])} / {len(d['rows'])} (draw, pot) pairs pass at the same tolerance; "
          f"largest {max(x['max_dex'] for x in d['rows']):.2e} dex.", ""] + exceptions(d["rows"]) + ["",
          "| seed | pot | drawn coordinates | sweeps (JAX / engine) | max dex | worst species | max pH diff |",
          "|---|---|---|---|---|---|---|"]
    for x in d["rows"]:
        L.append(f"| {x['seed']} | {x['id']} | {x['drawn_coordinates']} | {x['jax_sweeps'][0]} / {x['engine_sweeps'][0]}"
                 f" | {x['max_dex']:.1e} | {x['worst_species']} | {max(x['max_abs_node_ph'], x['max_abs_cooled_ph']):.1e} |")
    m = r["modes"]
    L += ["", "## The other modes (synthetic pots, not benchmarks) and the inline coupling", "",
          "| pot | coupling | sweeps (JAX / engine) | max dex | worst species | max pH diff | pass |", "|---|---|---|---|---|---|---|"]
    for x in m["rows"]:
        L.append(f"| {x['id']} | {x['coupling']} | {x['jax_sweeps']} / {x['engine_sweeps']} | {x['max_dex']:.1e} | "
                 f"{x['worst_species']} | {max(x['max_abs_node_ph'], x['max_abs_cooled_ph']):.1e} | {x['pass']} |")
    t = r["tight"]
    L += ["", f"## Tight tolerance (rtol {TIGHT_TOL[0]:g} / atol {TIGHT_TOL[1]:g}; contract sec. 3's 1e-6 dex)", "",
          "| pot | max dex | worst species | max pH diff | JAX s (first call, compile included) | engine s | pass |", "|---|---|---|---|---|---|---|"]
    for x in t["rows"]:
        L.append(f"| {x['id']} | {x['max_dex']:.1e} | {x['worst_species']} | "
                 f"{max(x['max_abs_node_ph'], x['max_abs_cooled_ph']):.1e} | {x['jax_s']:.2f} | {x['engine_s']:.2f} | {x['pass']} |")
    g = r["gradient"]
    L += ["", "## Gradient (`jax.grad`, reverse mode through the fixed point)", "",
          f"Scalar: a weighted sum of log10 {' / '.join(GRADIENT_SPECIES)} at the end of {g['pot']} plus the final "
          "cooled pH, at "
          f"rtol {g['rtol']:g}. Finite: {g['finite']}.", "",
          "| coordinate | jax.grad | forward mode | central FD (step) | rel. grad vs FD | rel. grad vs forward |",
          "|---|---|---|---|---|---|"]
    for i, c in enumerate(g["coordinates"]):
        L.append(f"| {c} | {g['grad'][i]:.6g} | {g['jacfwd'][i]:.6g} | {g['central_fd'][i]:.6g} ({g['steps'][i]:g}) | "
                 f"{g['rel_grad_vs_fd'][i]:.1e} | {g['rel_grad_vs_jacfwd'][i]:.1e} |")
    L += ["", f"jax.grad: {g['grad_s']:.1f} s per call ({g['grad_first_call_s']:.0f} s with its compile)."]
    f = r["fit_jacobian"]
    L += ["", "## One fit Jacobian at fit tolerance (paired)", "",
          f"{f['pot']}, {f['observables']} observables x {len(f['coordinates'])} coordinates "
          f"({', '.join(f['coordinates'])}), rtol {FIT_TOL[0]:g}.", "",
          "| quantity | wall s |", "|---|---|",
          f"| JAX, one solve | {f['jax_solve_s']:.2f} |",
          f"| JAX, the exact Jacobian (forward mode) | {f['jax_jacobian_s']:.2f} (first call with compile "
          f"{f['jax_jacobian_first_call_s']:.0f}) |",
          f"| engine, forward differences ({f['engine_runs_forward']} runs) | {f['engine_forward_fd_s']:.2f} |",
          f"| engine, central differences ({f['engine_runs_central']} runs) | {f['engine_central_fd_s']:.2f} |",
          "", f"JAX exact / engine forward FD: {f['jax_jacobian_s'] / f['engine_forward_fd_s']:.2f}x; / central: "
              f"{f['jax_jacobian_s'] / f['engine_central_fd_s']:.2f}x. Column error of the exact Jacobian against the "
              f"engine's central FD at these (fit-sized) steps {f['max_col_rel_jax_vs_engine_central']:.1e}; the engine's "
              f"forward FD against its own central {f['max_col_rel_engine_forward_vs_central']:.1e}, so at these steps "
              "the differences disagree with each other by as much (truncation over the large steps, and a step can "
              "move the fixed point's sweep count); at small steps the exact gradient equals central differences to "
              f"{max(r['gradient']['rel_grad_vs_fd']):.1e} (gradient section). "
              f"Load: {f['load']}.",
          "", "## What remains (not in this port)", "",
          "- Every off-by-default rival (B58 / B59a / B59b / B61, the pentose MGO routes, the MGO Strecker, the 3-DX trap, "
          "the hexose C5 entry, the pentose + amine entry, B7a R-T, B62, B63, B51 off, the trunk pH terms and the trunk "
          "solutes in the balance) is refused (`UnsupportedConfiguration`) and stays NumPy-only.",
          "- The pH solver lives in `sulfur_lane.integrate`, not in `network.integrate`: `Module.titratables` cannot yet "
          "carry the inventory's dH / dCp / acid yield; the ONE charge-closure check (`network.charge_closure`) is "
          "still the typed stub. The fit driver (`fit.py`) and the density bridge call `network.integrate`, so they do "
          "not take the sulfur lane yet.",
          "- Composing this module with `trunk` / `fed_3deoxy` (one owner per species) is A2b's union network.",
          "", f"Machine load (the last stage run): start {r.get('uptime_start', '')}; end {r.get('uptime_end', '')}.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", nargs="*", choices=STAGES, help="run only these stages (merged into the record)")
    parser.add_argument("--render", action="store_true", help="rewrite the .md from the .json only")
    args = parser.parse_args(argv)
    import jax

    path = OUT.with_suffix(".json")
    if args.render:
        OUT.with_suffix(".md").write_text(markdown(json.loads(path.read_text())), encoding="utf-8")
        return 0
    report = json.loads(path.read_text()) if path.exists() and args.stage else {}
    report.update({"artifact": "jax_sulfur_lane_module", "generated_by": "scripts/generators/build_jax_sulfur_lane_module.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "uptime_start": uptime()})
    h = Harness()
    t_all = time.perf_counter()
    for stage in args.stage or STAGES:
        t0 = time.perf_counter()
        print(f"stage {stage} ...", flush=True)
        report[stage] = globals()[f"stage_{stage}"](h)
        report.setdefault("stage_seconds", {})[stage] = time.perf_counter() - t0
        path.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    report["uptime_end"] = uptime()
    report["wall_s"] = time.perf_counter() - t_all
    path.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    if all(s in report for s in STAGES):
        OUT.with_suffix(".md").write_text(markdown(report), encoding="utf-8")
    print(f"wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

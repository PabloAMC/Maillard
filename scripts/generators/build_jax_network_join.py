#!/usr/bin/env python
"""
build_jax_network_join.py -- the one composable JAX network with the sulfur lane in it (plan step A2 + H3, the
11-20 checkpoint "the one network runs with exact gradients, sulfur included"; 2026-11-06): the network integrator
owns the pH solve (``ph_coupling``: the merged titratable table with dH / dCp / the acid-yield coupling, the charge
balance and its implicit-function derivative, both couplings), and the sulfur lane runs through it -- as one module
and as three COMPOSED modules (``sulfur_lane.build_composed``: ``sulfur_trunk`` + ``sulfur_fed_3deoxy`` +
``sulfur``, declared competition refused) -- in ``network.integrate``, ``fit.py`` and ``densities.py``.

Runs in the SEPARATE ``maillard-jax`` env with the ENGINE IN THE SAME PROCESS (every compared engine number is
computed live; timings paired under one machine load). One core, niced, threads pinned:

    ./scripts/docker_maillard.sh run "set +u; conda activate maillard-jax; \
        XLA_FLAGS='--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1' OMP_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 nice -n 10 python scripts/generators/build_jax_network_join.py [--stage composed ...]"

Writes results/validation/jax_network_join.{json,md}. Stages (each merged into the record):

* ``structure`` -- the composed network equals the one-module lane (reactions, species, stoichiometry), declares
  every shared consumption, and its merged table passes the ONE charge-closure check; the movement it computes equals
  the engine's ``ph_state.centre_delta`` on every step; the charge each step moves at the label pH 6;
* ``rates`` -- the composed network's rate constants against the engine's ``sulfur_rate_constants_at`` (1e-12);
* ``composed`` -- five panel pots (glucose + cysteine with the trunk active, thiamine, an unbuffered ARP pot whose pH
  drifts, a buffered xylose pot): the composed network against the engine at the engine's shipped tolerance, under
  the sulfur parity test's own bound (1e-4 dex, ``build_jax_sulfur_lane_module.within``), and against the one-module
  lane; paired solve timings;
* ``gradient`` -- ``jax.grad`` through the composed network on five coordinates of all three modules, with the
  DEFAULT adjoint (no ``reverse`` flag: the fixed point's sweeps run in equinox's checkpointed while loop), against
  central differences and forward mode, at rtol 1e-10;
* ``fit_jacobian`` -- ``fit.Compiled`` on the composed network (the fit driver taking the sulfur lane): the residual
  equals ``network.integrate``'s, the exact Jacobian against ``fit.fd_check``'s central differences, its wall time
  against the engine's forward-difference Jacobian on the same pot (paired);
* ``densities`` -- ``densities.Predictor`` on the composed network: a batch of draws equals single solves.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"
for p in (ROOT, GEN):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import build_jax_sulfur_lane_module as g  # noqa: E402  (the sulfur parity harness: pots, engine run, tolerance)

OUT = ROOT / "results" / "validation" / "jax_network_join"
COMPOSED_POTS = ("mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7", "thiamine_cys_glucose_120C_Bolton1994",
                 "resconi_2023_pbma_beef_identity_benchmark", "mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c",
                 "mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b")
GRADIENT_POT = "thiamine_cys_glucose_120C_Bolton1994"
GRADIENT_SPECIES = ("MFT", "H2S", "HMF", "FFT")
FD_STEP = {"log10_k": 1e-4, "ea": 1e-2, "acid_yield": 2e-5, "arp_amine_pka": 1e-3}
STAGES = ("structure", "rates", "composed", "gradient", "fit_jacobian", "fit_column_diagnostic", "densities")


def gradient_coordinates(lane) -> List[tuple]:
    """Five coordinates over the three modules that move the gradient pot (glucose + cysteine + thiamine, no amino
    acid, unbuffered): the trunk's fructose -> acids step (its formic / acetic acid drift the pH through the charge
    balance) and the Ea of glucose -> 3-DG, the fed triangle's first constant, the thiamine route's constant and
    HMF + cysteine. (The Schiff step and the drift's acid yield are idle here: no glycine, no lumped acid.)"""
    from src.kinetic_core_jax import sulfur_structure as ss

    pk = ss.derive()["parameter_key"]
    fed = next(k for k in next(m for m in lane.modules if m.name == "sulfur_fed_3deoxy").params
               if k.startswith("log10_k:"))
    return [("sulfur_trunk", f"log10_k:{pk['r_fru_acids']}"), ("sulfur_trunk", f"ea:{pk['r_glc_tdg']}"),
            ("sulfur_fed_3deoxy", fed), ("sulfur", f"log10_k:{pk['r_thi_hmp']}"),
            ("sulfur", f"log10_k:{pk['r_hmf_cys']}")]


def step_of(coord) -> float:
    name = coord[1]
    return FD_STEP.get(name.partition(":")[0], FD_STEP.get(name))


class Join:
    """The one-module lane, the composed lane and their parity harnesses."""

    def __init__(self):
        from src.kinetic_core_jax import sulfur_lane as sl

        self.sl = sl
        self.one = g.Harness()
        self.composed = g.Harness(lane=sl.build_composed())
        self.pots = {p["id"]: p for p in g.panel_pots() if p["id"] in COMPOSED_POTS}
        missing = sorted(set(COMPOSED_POTS) - set(self.pots))
        if missing:
            raise SystemExit(f"panel pots not routed to the sulfur lane any more: {missing}")

    def pot(self, h, pid):
        return h.prepare(dict(self.pots[pid]))


# --------------------------------------------------------------------------- stages
def stage_structure(j: Join) -> Dict[str, Any]:
    from src.kinetic_core import ph_state, sulfur
    from src.kinetic_core_jax import network

    one, comp = j.one.lane.network, j.composed.lane.network
    closure = network.charge_closure(comp, ph=6.0)
    closure_one = network.charge_closure(one)
    engine_moved = {}
    for r in sulfur.FULL_REACTIONS:
        d = {k: v for k, v in ph_state.centre_delta(r.reactants, r.products).items() if v != 0}
        if d:
            engine_moved[r.key] = d
    same = set(engine_moved) == set(closure["moved"]) and all(
        all(abs(closure["moved"][r].get(k, 0.0) - v) < 1e-9 for k, v in d.items()) for r, d in engine_moved.items())
    sc = comp.charge.structural_weights().T @ comp.charge.centres.T
    centres_equal = all(
        all(abs(sc[comp.index[s], comp.charge.kinds.index(k)] - v) < 1e-12 for k, v in c.items())
        for s, c in ph_state.TITRATABLE_CENTRES.items() if s in comp.index)
    return {
        "modules": [{"name": m.name, "species": len(m.species), "reactions": len(m.reactions),
                     "coordinates": len(m.params), "requires": list(m.requires), "competes": list(m.competes),
                     "titratables": len(m.titratables), "ledger_entries": len(m.centre_ledger)}
                    for m in comp.modules],
        "reactions_equal": tuple(comp.reactions) == tuple(one.reactions),
        "species_equal": tuple(comp.species) == tuple(one.species),
        "stoichiometry_equal": bool(np.array_equal(comp.stoichiometry, one.stoichiometry)),
        "undeclared_competition": list(comp.undeclared_competition),
        "shared_species": len(comp.competition),
        "pools": list(comp.charge.names), "pot_pools": sorted(comp.charge.buffer_rows),
        "yield_couplings": [list(k) for k, _ in comp.charge.yields],
        "pka_coordinates": [list(k) for k, _ in comp.charge.pka_slots],
        "closure_composed": {k: closure[k] for k in ("closes", "problems", "kinds")},
        "closure_one_module": {k: closure_one[k] for k in ("closes", "problems")},
        "steps_moving_centres": len(closure["moved"]),
        "movement_equals_engine_centre_delta": bool(same),
        "species_centres_equal_TITRATABLE_CENTRES": bool(centres_equal),
        "charge_moved_at_label_ph6": closure["charge_moved"],
        "ph_mode": comp.ph_mode,
    }


def stage_rates(j: Join) -> Dict[str, Any]:
    from src.kinetic_core import engine

    h = j.composed
    out = {"shipped": g.rate_probe(h, engine.core_parameters(engine.SULFUR), engine.core_ph_drift())}
    p, d, _ = g.draw(g.DRAW_SEEDS[0])
    out[f"draw_{g.DRAW_SEEDS[0]}"] = g.rate_probe(h, p, d)
    return out


def stage_composed(j: Join) -> Dict[str, Any]:
    from src.kinetic_core import engine

    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    rows = []
    for pid in COMPOSED_POTS:
        pc, po = j.pot(j.composed, pid), j.pot(j.one, pid)
        run_c, first_c = j.composed.solve(pc)
        run_o, _ = j.one.solve(po)
        eng = g.engine_run(j.composed.lane, pc, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        run_c, c = j.composed.matched(pc, run_c, eng)
        run_o, co = j.one.matched(po, run_o, eng)
        yc, yo = np.asarray(run_c.ys), np.asarray(run_o.ys)
        tc = [j.composed.solve(pc)[1] for _ in range(2)]
        to = [j.one.solve(po)[1] for _ in range(2)]
        te = [g.engine_run(j.composed.lane, pc, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])["seconds"]
              for _ in range(2)]
        rec = {"id": pid, "initial": pc["initial"], **c, "pass": not g.within(c, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1]),
               "one_module_max_dex": co["max_dex"],
               "composed_vs_one_module_max_abs": float(np.max(np.abs(yc - yo))),
               "composed_vs_one_module_max_rel": float(np.max(np.abs(yc - yo) / np.maximum(np.abs(yo), 1e-300))),
               "composed_s": min(tc), "one_module_s": min(to), "engine_s": min(te + [eng["seconds"]]),
               "composed_first_call_s": first_c, "steps": int(np.sum(np.asarray(run_c.num_steps))),
               "load": g.load()}
        rows.append(rec)
        print(f"  {pid}: {c['max_dex']:.2e} dex (one module {co['max_dex']:.2e}), sweeps {c['jax_sweeps']}/"
              f"{c['engine_sweeps']}, composed {rec['composed_s']:.2f} s, engine {rec['engine_s']:.2f} s", flush=True)
    return {"tolerance": {"rtol": g.FIT_TOL[0], "atol": g.FIT_TOL[1], "dex": g.DEX_TOL, "ph": g.PH_TOL},
            "pots": rows}


def stage_gradient(j: Join) -> Dict[str, Any]:
    """jax.grad (default adjoint: the checkpointed sweep loop) of one scalar through ``network.integrate`` on the
    composed network, against central differences and forward mode, at rtol 1e-10."""
    import diffrax
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import network

    lane = j.composed.lane
    net = lane.network
    pot = j.pot(j.composed, GRADIENT_POT)["jax"]
    arrays, mode = network.pot_arrays(net, pot), network.pot_mode(net, pot)
    coords = gradient_coordinates(lane)
    theta0 = net.theta0()
    idx = jnp.asarray([net.index[s] for s in GRADIENT_SPECIES])
    weights = jnp.asarray([1.0, 0.5, 0.25, 0.125])
    rtol, atol = g.TIGHT_TOL

    def scalar(x, adjoint):
        th = {m: dict(v) for m, v in theta0.items()}
        for i, (m, k) in enumerate(coords):
            th[m][k] = x[i]
        run = network.integrate(net, th, arrays, mode=mode, adjoint=adjoint, rtol=rtol, atol=atol)
        logs = jnp.log10(jnp.maximum(run.ys[-1, idx], 1e-300))
        return jnp.sum(logs * weights) + run.ph.ph_cool[-1, -1]

    x0 = jnp.stack([theta0[m][k] for m, k in coords])
    t0 = time.perf_counter()
    grad = jax.jit(jax.grad(lambda x: scalar(x, diffrax.RecursiveCheckpointAdjoint())))
    gv = np.asarray(grad(x0))
    first = time.perf_counter() - t0
    t0 = time.perf_counter()
    gv = np.asarray(grad(x0))
    t_grad = time.perf_counter() - t0
    fwd = jax.jit(jax.jacfwd(lambda x: scalar(x, diffrax.ForwardMode())))
    t0 = time.perf_counter()
    jf = np.asarray(fwd(x0))
    t_fwd_first = time.perf_counter() - t0
    f = jax.jit(lambda x: scalar(x, diffrax.RecursiveCheckpointAdjoint()))
    steps = [step_of(c) for c in coords]
    fd = []
    for i in range(len(coords)):
        e = np.zeros(len(coords))
        e[i] = steps[i]
        fd.append((float(f(x0 + e)) - float(f(x0 - e))) / (2.0 * steps[i]))
    fd = np.asarray(fd)
    scale = np.maximum(np.abs(fd), 1e-8)
    return {"pot": GRADIENT_POT, "coordinates": [f"{m}.{k}" for m, k in coords], "steps": steps,
            "species": list(GRADIENT_SPECIES), "grad": gv.tolist(), "jacfwd": jf.tolist(), "central_fd": fd.tolist(),
            "rel_grad_vs_fd": (np.abs(gv - fd) / scale).tolist(),
            "rel_grad_vs_jacfwd": (np.abs(gv - jf) / scale).tolist(), "finite": bool(np.all(np.isfinite(gv))),
            "grad_s": t_grad, "grad_first_call_s": first, "jacfwd_first_call_s": t_fwd_first, "rtol": rtol,
            "atol": atol, "adjoint": "RecursiveCheckpointAdjoint (default; no reverse flag)", "load": g.load()}


def _fit_problem(j: Join):
    from src.kinetic_core_jax import fit, network

    lane = j.composed.lane
    net = lane.network
    pot = j.pot(j.composed, GRADIENT_POT)["jax"]
    n = np.asarray(pot["node_s"]).shape[1]
    node_rows = [n // 8, n // 3, (2 * n) // 3, n - 1]
    coords = gradient_coordinates(lane)
    run = network.integrate(net, net.theta0(), network.pot_arrays(net, pot), mode=network.pot_mode(net, pot),
                            rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
    ys = np.asarray(run.ys)
    obs = [fit.Observation(pot=0, row=r, observable=f"mmol_per_L:{s}",
                           value=math.log10(max(float(ys[r, net.index[s]]), 1e-300)) + 0.05, sigma=0.1)
           for s in GRADIENT_SPECIES for r in node_rows]
    base = net.theta0()
    x0 = np.array([float(base[m][k]) for m, k in coords])
    width = np.array([1.0 if k.startswith("log10_k") else 20.0 if k.startswith("ea") else 1e-3 for _, k in coords])
    lower = x0 - width
    lower[[i for i, (_, k) in enumerate(coords) if k == "acid_yield"]] = 0.0
    problem = fit.JaxProblem(net=net, keys=tuple(coords), kinds=tuple(fit.kind_of(k) or "acid_yield" for _, k in coords),
                             lower=lower, upper=x0 + width, pots=(pot,), observations=tuple(obs),
                             solve={"rtol": g.FIT_TOL[0], "atol": g.FIT_TOL[1]}, batch="loop")
    return problem, x0, ys, node_rows, coords


def stage_fit_jacobian(j: Join) -> Dict[str, Any]:
    import dataclasses

    from src.kinetic_core import engine
    from src.kinetic_core_jax import fit

    problem, x0, ys, node_rows, coords = _fit_problem(j)
    comp = fit.Compiled(problem)
    r, _ = comp.residual(x0)
    expect = np.array([(math.log10(max(float(ys[o.row, problem.net.index[o.observable.split(':')[1]]]), 1e-300))
                        - o.value) / o.sigma for o in problem.observations])
    t0 = time.perf_counter()
    _, jac, failed = comp.jacobian(x0)
    first = time.perf_counter() - t0
    times = []
    for _ in range(2):
        t0 = time.perf_counter()
        comp.jacobian(x0)
        times.append(time.perf_counter() - t0)
    check = fit.fd_check(problem, x0, compiled=comp)
    # the engine's forward-difference Jacobian on the same pot and observables (1 + p runs), paired
    pot = j.pot(j.composed, GRADIENT_POT)
    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()

    def moved(coord, h):
        _, name = coord
        kind, _, key = name.partition(":")
        p, d = dict(params), drift
        if kind == "log10_k":
            p[key] = dataclasses.replace(p[key], k_ref=float(p[key].k_ref) * 10.0 ** h)
        elif kind == "ea":
            p[key] = dataclasses.replace(p[key], ea_kj_mol=float(p[key].ea_kj_mol) + h)
        else:
            d = dataclasses.replace(d, acid_yield=d.acid_yield + h)
        return p, d

    def engine_obs(p, d):
        s = g.engine_run(j.composed.lane, pot, p, d, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])["states"][-1]
        return np.array([math.log10(max(s[r, problem.net.index[sp]], 1e-300)) for sp in GRADIENT_SPECIES
                         for r in node_rows])

    t0 = time.perf_counter()
    base = engine_obs(params, drift)
    cols = [(engine_obs(*moved(c, 10 * step_of(c))) - base) / (10 * step_of(c)) for c in coords]
    t_eng = time.perf_counter() - t0
    jac_dex = jac * np.array([o.sigma for o in problem.observations])[:, None]
    eng = np.array(cols).T
    col_scale = np.maximum(np.max(np.abs(eng), axis=0), 1e-12)
    t0 = time.perf_counter()
    comp.residual(x0)
    t_res = time.perf_counter() - t0
    return {"pot": GRADIENT_POT, "observables": len(problem.observations), "coordinates": [f"{m}.{k}" for m, k in coords],
            "residual_equals_network_integrate": float(np.max(np.abs(r - expect))), "failed": failed,
            "jax_jacobian_s": min(times), "jax_jacobian_first_call_s": first, "jax_residual_s": t_res,
            "fd_check_max_relative_column_error": check["max_relative_column_error"],
            "fd_check_per_column": check["per_column"], "engine_forward_fd_s": t_eng,
            "engine_runs": 1 + len(coords), "engine_steps": [10 * step_of(c) for c in coords],
            "max_col_rel_jax_vs_engine_forward_fd": float(np.max(np.max(np.abs(jac_dex - eng), axis=0) / col_scale)),
            "col_rel_jax_vs_engine_forward_fd": (np.max(np.abs(jac_dex - eng), axis=0) / col_scale).tolist(),
            "col_scale_dex_per_unit": np.max(np.abs(jac_dex), axis=0).tolist(),
            "load": g.load()}


def stage_fit_column_diagnostic(j: Join) -> Dict[str, Any]:
    """Why the worst column of the fit-Jacobian comparison disagrees with the engine's forward difference: the
    composed network and the engine at that coordinate's -h / 0 / +h (h the engine's step), with both sweep counts
    and the two sides' largest difference on the 16 observables (a sweep-count jump, or noise against a step)."""
    import dataclasses

    from src.kinetic_core import engine

    f = j_record().get("fit_jacobian", {})
    cols = f.get("col_rel_jax_vs_engine_forward_fd")
    _, _, _, node_rows, coords = _fit_problem(j)
    worst = int(np.argmax(cols)) if cols else 0
    module, name = coords[worst]
    kind, _, key = name.partition(":")
    h = 10 * step_of(coords[worst])
    pot = j.pot(j.composed, GRADIENT_POT)
    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    idx = [j.composed.lane.index[s] for s in GRADIENT_SPECIES]
    rows = []
    for x in (-h, 0.0, h):
        p = dict(params)
        if kind == "log10_k":
            p[key] = dataclasses.replace(p[key], k_ref=float(p[key].k_ref) * 10.0 ** x)
        else:
            p[key] = dataclasses.replace(p[key], ea_kj_mol=float(p[key].ea_kj_mol) + x)
        e = g.engine_run(j.composed.lane, pot, p, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        th = j.composed.lane.network.theta0()
        th[module][name] = th[module][name] + x
        run, _ = j.composed.solve(pot, th)
        ev = np.array([math.log10(e["states"][-1][r, i]) for i in idx for r in node_rows])
        jv = np.array([math.log10(float(np.asarray(run.ys)[-1, r, i])) for i in idx for r in node_rows])
        labels = [(s, r) for s in GRADIENT_SPECIES for r in node_rows]
        w = int(np.argmax(np.abs(ev - jv)))
        rows.append({"shift": x, "engine_sweeps": [int(s) for s in e["sweeps"]],
                     "jax_sweeps": [int(s) for s in np.asarray(run.sweeps)],
                     "max_abs_dex_jax_vs_engine": float(np.max(np.abs(ev - jv))),
                     "worst_observable": f"{labels[w][0]} at node {labels[w][1]}",
                     "worst_level_mmol_L": float(10.0 ** ev[w]), "engine_change_dex": None})
        rows[-1]["_ev"] = ev
    base = rows[1]["_ev"]
    for r in rows:
        r["engine_change_dex"] = float(np.max(np.abs(r.pop("_ev") - base)))
    return {"coordinate": f"{module}.{name}", "step": h, "rows": rows}


def j_record() -> Dict[str, Any]:
    path = OUT.with_suffix(".json")
    return json.loads(path.read_text()) if path.exists() else {}


def stage_densities(j: Join) -> Dict[str, Any]:
    import jax.numpy as jnp
    from src.kinetic_core_jax import densities, network

    problem, x0, _, _, coords = _fit_problem(j)
    net = problem.net
    rng = np.random.default_rng(5)
    draws = x0[None, :] + rng.normal(size=(3, len(x0))) * np.array([0.05 if k.startswith("log10_k") else 1.0
                                                                    if k.startswith("ea") else 1e-5 for _, k in coords])
    pred = densities.Predictor(net, problem.pots, problem.observations, coords, chunk=3, rtol=g.FIT_TOL[0],
                               atol=g.FIT_TOL[1])
    t0 = time.perf_counter()
    out, failed = pred(draws)
    t_batch = time.perf_counter() - t0
    pot = problem.pots[0]
    worst = 0.0
    for i, x in enumerate(draws):
        th = network.theta_from_vector(net, network.theta_vector(net, net.theta0()).at[
            jnp.asarray([net.coordinates.index(c) for c in coords])].set(jnp.asarray(x)))
        run = network.integrate(net, th, network.pot_arrays(net, pot), mode=network.pot_mode(net, pot),
                                rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        ys = np.asarray(run.ys)
        single = np.array([math.log10(max(float(ys[o.row, net.index[o.observable.split(':')[1]]]), 1e-300))
                           for o in problem.observations])
        worst = max(worst, float(np.max(np.abs(single - out[i]))))
    return {"draws": int(len(draws)), "failed": int(np.sum(failed)), "batch_vs_single_max_dex": worst,
            "batch_s_first_call": t_batch, "load": g.load()}


# --------------------------------------------------------------------------- the record
def markdown(r: Mapping[str, Any]) -> str:
    s = r.get("structure", {})
    L = ["# The one JAX network, sulfur included (A2 + H3, 2026-11-06)", "",
         f"*Generated by `scripts/generators/build_jax_network_join.py` in the separate `maillard-jax` env (jax "
         f"{r.get('jax_version', '?')}, float64, CPU, one core, niced), the engine in the same process (live, paired "
         "timings). Code: `src/kinetic_core_jax/ph_coupling.py` (the network's pH solve), `network.py` (v2.3 "
         "contract), `sulfur_lane.py` (one-module and composed forms), `fit.py` / `densities.py`. Nothing shipped "
         "moves; `src/kinetic_core` is not touched and does not import it.*", "",
         "## What moved where", "",
         "- The charge balance, its (pH, I) root with the implicit-function derivative, the two couplings "
         "(\"fixed_point_nodes\", shipped, and \"inline\"), the clamped pH, the node grids and the pot builder moved "
         "from `sulfur_lane.py` into `ph_coupling.py`, unchanged operation for operation; `network.integrate` calls it "
         "for any network with a `needs_ph` module.",
         "- `network.Titratable` carries what the engine's inventory needs: dH and dCp per step, pools fed by several "
         "species (`weights`), the acid-yield coupling (`yield_weights` x a coordinate), pK coordinates, pot-supplied "
         "(buffer) pools and the centres per pool. `Module.ph_model` declares the solvent and the coupling's numerics "
         "(one per network); `Module.centre_ledger` the steps that move centres.",
         "- `network.charge_closure` is the ONE check over the merged table (was a typed stub).",
         "- `fit.py` and `densities.py` take inline-pH pots (their mode is a compile key; their rows are every pH "
         "node, `network.segment_end_rows`).", ""]
    if s:
        L += ["## Structure and the charge closure", "",
              "| module | species | reactions | coordinates | requires | declared competition | pools | ledger entries |",
              "|---|---|---|---|---|---|---|---|"]
        for m in s["modules"]:
            L.append(f"| {m['name']} | {m['species']} | {m['reactions']} | {m['coordinates']} | "
                     f"{', '.join(m['requires']) or '-'} | {len(m['competes'])} species | {m['titratables']} | "
                     f"{m['ledger_entries']} |")
        L += ["",
              f"- Composed with `competition=\"refuse\"`: undeclared pairs {s['undeclared_competition'] or 'none'}; "
              f"{s['shared_species']} species consumed by a module that does not own them. Reactions / species / "
              f"stoichiometry equal to the one-module lane: {s['reactions_equal']} / {s['species_equal']} / "
              f"{s['stoichiometry_equal']}. pH mode {s['ph_mode']}.",
              f"- Merged table: {len(s['pools'])} pools ({', '.join(s['pools'])}); pot-supplied: "
              f"{', '.join(s['pot_pools'])}; yield coupling {s['yield_couplings']}; pK coordinate {s['pka_coordinates']}.",
              f"- Charge closure (composed): closes {s['closure_composed']['closes']}, problems "
              f"{s['closure_composed']['problems'] or 'none'}; one-module form: closes "
              f"{s['closure_one_module']['closes']}. {s['steps_moving_centres']} steps move centres; the movement "
              f"computed from the merged table equals the engine's `ph_state.centre_delta` on every step: "
              f"{s['movement_equals_engine_centre_delta']}; species centres equal `TITRATABLE_CENTRES`: "
              f"{s['species_centres_equal_TITRATABLE_CENTRES']}.", ""]
        cm = s["charge_moved_at_label_ph6"]["per_reaction"]
        top = sorted(cm.items(), key=lambda kv: -abs(kv[1]))[:8]
        L += [f"- Charge moved per event at label pH 6, 25 C, I = 0 (the protons a step puts into the balance), largest: "
              + ", ".join(f"{k} {v:+.3f}" for k, v in top) + ".", ""]
    if "rates" in r:
        for label, rp in r["rates"].items():
            L.append(f"- Rate constants of the composed network, {label}: {rp['values']} values, worst relative "
                     f"difference from the engine {rp['worst_relative']:.1e}; non-zero where the engine is zero: "
                     f"{rp['nonzero_where_engine_zero']}.")
        L.append("")
    if "composed" in r:
        c = r["composed"]
        pots = c["pots"]
        L += ["## The composed recipe on five panel pots (against the engine, live)", "",
              f"Bound: the sulfur parity test's own (rtol {g.FIT_TOL[0]:g} / atol {g.FIT_TOL[1]:g} both sides; "
              f"{g.DEX_TOL:g} dex above atol / rtol, pH {g.PH_TOL:g}, sweep counts equal). **"
              f"{sum(p['pass'] for p in pots)} / {len(pots)} pass**; largest {max(p['max_dex'] for p in pots):.2e} dex.",
              "No panel pot routes to the sulfur and the acrylamide lane together (25 sulfur-only, 5 acrylamide-only, "
              "1 acrylamide + lipid), so the acrylamide module has no place in this recipe; composing the declared-pH "
              "trunk / acrylamide modules with the inline-pH sulfur lane is refused (sec. 2.4, A2b).", ""]
        L += g.exceptions(pots)
        L += ["", "| pot | sweeps (composed / engine) | max dex vs engine | worst species | one-module max dex | "
                  "composed vs one-module (max abs, mmol/L) | composed s | one-module s | engine s | load |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for p in pots:
            mark = ("*" if p.get("tie_broken") else "") + ("j" if p.get("jump_side_matched") else "")
            L.append(f"| {p['id']} | {p['jax_sweeps'][0]}{mark} / {p['engine_sweeps'][0]} | {p['max_dex']:.1e} | "
                     f"{p['worst_species']} | {p['one_module_max_dex']:.1e} | {p['composed_vs_one_module_max_abs']:.1e} | "
                     f"{p['composed_s']:.2f} | {p['one_module_s']:.2f} | {p['engine_s']:.2f} | {p['load']} |")
        ratio = np.array([p["composed_s"] for p in pots]) / np.array([p["engine_s"] for p in pots])
        L += ["", f"Composed solve / engine at equal tolerance: median {np.median(ratio):.2f}x (range "
                  f"{ratio.min():.2f} - {ratio.max():.2f}); first call with compile "
                  f"{max(p['composed_first_call_s'] for p in pots):.0f} s at most. PROVISIONAL as absolute numbers "
                  "(one core of a shared machine; load beside each row); the ratio is the reading.", ""]
    if "gradient" in r:
        gr = r["gradient"]
        L += ["## `jax.grad` through the composed network", "",
              f"Scalar: a weighted sum of log10 {' / '.join(gr['species'])} at the end of {gr['pot']} plus the final "
              f"cooled pH, rtol {gr['rtol']:g}; adjoint: {gr['adjoint']}. Finite: {gr['finite']}.", "",
              "| coordinate | jax.grad | forward mode | central FD (step) | rel. grad vs FD | rel. grad vs forward |",
              "|---|---|---|---|---|---|"]
        for i, cname in enumerate(gr["coordinates"]):
            L.append(f"| {cname} | {gr['grad'][i]:.6g} | {gr['jacfwd'][i]:.6g} | {gr['central_fd'][i]:.6g} "
                     f"({gr['steps'][i]:g}) | {gr['rel_grad_vs_fd'][i]:.1e} | {gr['rel_grad_vs_jacfwd'][i]:.1e} |")
        L += ["", f"jax.grad: {gr['grad_s']:.1f} s per call ({gr['grad_first_call_s']:.0f} s with its compile); "
                  f"forward mode first call {gr['jacfwd_first_call_s']:.0f} s. Load: {gr['load']}.", ""]
    if "fit_jacobian" in r:
        f = r["fit_jacobian"]
        L += ["## The fit driver on the composed network", "",
              f"{f['pot']}: {f['observables']} log10 observables (four species at four pH nodes) x "
              f"{len(f['coordinates'])} coordinates ({', '.join(f['coordinates'])}), rtol {g.FIT_TOL[0]:g}.", "",
              f"- `fit.Compiled.residual` equals the residual of `network.integrate`'s states to "
              f"{f['residual_equals_network_integrate']:.1e}; failed pots {f['failed'] or 'none'}.",
              f"- The exact Jacobian against `fit.fd_check`'s central differences of the same residual (column-scaled "
              f"steps, at the fit tolerance, so solver noise / step bounds the agreement), per column: "
              + ", ".join(f"{x:.1e}" for x in f["fd_check_per_column"]) + ".",
              f"- Against the engine's own forward differences at steps {f['engine_steps']} (one-sided: truncation "
              "O(step), plus LSODA's noise / step), per column: "
              + ", ".join(f"{x:.1e}" for x in f.get("col_rel_jax_vs_engine_forward_fd", [])) + "; the columns' "
              "scale (max |d log10 / d coordinate|): "
              + ", ".join(f"{x:.2g}" for x in f.get("col_scale_dex_per_unit", [])) + ".", ""]
        dg = r.get("fit_column_diagnostic")
        if dg:
            L += [f"- The worst column, {dg['coordinate']} at step {dg['step']:g}, diagnosed (engine and composed "
                  "network at -step / 0 / +step):", "",
                  "| shift | sweeps (engine / JAX) | engine's largest change from 0 (dex) | JAX vs engine, largest (dex) "
                  "| on (level, mmol/L) |", "|---|---|---|---|---|"]
            for x in dg["rows"]:
                L.append(f"| {x['shift']:+g} | {x['engine_sweeps'][0]} / {x['jax_sweeps'][0]} | "
                         f"{x['engine_change_dex']:.1e} | {x['max_abs_dex_jax_vs_engine']:.1e} | "
                         f"{x.get('worst_observable', '')} ({x.get('worst_level_mmol_L', float('nan')):.1e}) |")
            counts = {(x["engine_sweeps"][0], x["jax_sweeps"][0]) for x in dg["rows"]}
            low = all(x.get("worst_level_mmol_L", 1.0) < g.FIT_TOL[1] / g.FIT_TOL[0] for x in dg["rows"])
            where = (" The largest differences sit on a level below atol / rtol = 1e-6 mmol/L, where the solvers "
                     "control only the absolute error (atol 1e-14): the observation set includes a species that "
                     "early in the run is at the solvers' atol, which a real fit would floor or drop.") if low else ""
            if len(counts) == 1:
                L += ["", "  Equal sweep counts at all three points, so no sweep-count jump: the step moves the "
                          "observables by about as much as the two solvers differ at rtol 1e-8 (inside the parity "
                          "bound), so a forward difference at this step reads solver noise." + where
                      + " The exact Jacobian is "
                          "free of it; at tight tolerance `jax.grad` equals central differences in this column "
                          "(previous section).", ""]
            else:
                L += ["", "  The sweep count CHANGES across the step: a jump of the discrete stopping rule, which "
                          "a difference sees and no derivative does (What remains).", ""]
        L += [
              "| quantity | wall s |", "|---|---|",
              f"| JAX residual (one solve) | {f['jax_residual_s']:.2f} |",
              f"| JAX exact Jacobian (forward mode) | {f['jax_jacobian_s']:.2f} (first call with compile "
              f"{f['jax_jacobian_first_call_s']:.0f}) |",
              f"| engine forward-difference Jacobian ({f['engine_runs']} runs) | {f['engine_forward_fd_s']:.2f} |",
              "", f"JAX exact / engine forward FD: {f['jax_jacobian_s'] / f['engine_forward_fd_s']:.2f}x. Load: "
                  f"{f['load']}.", ""]
    if "densities" in r:
        d = r["densities"]
        L += ["## The density bridge on the composed network", "",
              f"`densities.Predictor`: {d['draws']} draws of the five coordinates in one batch equal their single "
              f"solves to {d['batch_vs_single_max_dex']:.1e} dex; failed {d['failed']}.", ""]
    L += ["## What remains", "",
          "- **The sweep count is blind to its own jumps.** The fixed point stops on a discrete rule (residual < "
          "1e-3 pH after at least the minimum sweeps); every derivative here is the derivative AT a fixed sweep count. "
          "Where a perturbation changes the count, the trajectory jumps (by up to ~the tolerance in pH) and no "
          "derivative sees it; a central difference that straddles the change does. The run reports the count "
          "(`PhRun.sweeps`), so a fit can flag pots whose count moved between evaluations; a smooth surrogate (a "
          "fixed sweep count per fit, or the implicit-function derivative of the converged fixed point) is a model "
          "change and needs its own key.",
          "- **Inline-pH table jumps** (Zheng & Ho's boundary at label pH 5.0, the thiamine floor at 6.25) are crossed "
          "at a time that depends on theta; the derivative misses the crossing term (contract sec. 2.4, unchanged).",
          "- **The cooled pH is not yet a Network observable** (it needs SID, theta and the buffer, not one state): "
          "fits read it from `Run.ph`; `fit.Observation` reads species only.",
          "- **Two pH conventions in one network** (a declared-pH trunk / acrylamide module beside the inline-pH "
          "sulfur lane) are refused; merging them, the trunk solutes in the balance and one titratable table across "
          "lanes are A2b (contract sec. 6.3), judged by evidence, not identity.",
          "- **forcing.py** still refuses an inline-pH network (a forced pH there is the clamped mode; a forced species "
          "beside an inline pH is not ported).",
          "- **Declared competition** is refused in the composed sulfur recipe; `compose`'s default stays \"report\" "
          "for the trunk-lane / acrylamide modules, which predate the declaration.",
          "- **The pH-root performance step** (a cheaper residual; one root per step) is untouched.",
          "", f"Machine load: start {r.get('uptime_start', '')}; end {r.get('uptime_end', '')}.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", nargs="*", choices=STAGES, help="run only these stages (merged into the record)")
    parser.add_argument("--render", action="store_true", help="rewrite the .md from the .json only")
    args = parser.parse_args(argv)
    path = OUT.with_suffix(".json")
    if args.render:
        OUT.with_suffix(".md").write_text(markdown(json.loads(path.read_text())), encoding="utf-8")
        return 0
    import jax

    report = json.loads(path.read_text()) if path.exists() and args.stage else {}
    report.update({"artifact": "jax_network_join", "generated_by": "scripts/generators/build_jax_network_join.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "uptime_start": g.uptime()})
    j = Join()
    for stage in args.stage or STAGES:
        t0 = time.perf_counter()
        print(f"stage {stage} ...", flush=True)
        report[stage] = globals()[f"stage_{stage}"](j)
        report.setdefault("stage_seconds", {})[stage] = time.perf_counter() - t0
        path.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    report["uptime_end"] = g.uptime()
    path.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    OUT.with_suffix(".md").write_text(markdown(report), encoding="utf-8")
    print(f"wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

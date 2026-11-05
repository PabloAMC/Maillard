#!/usr/bin/env python
"""
build_jax_trunk_lane_module.py -- the trunk LANE (``network.TRUNK_REACTIONS``, 55 steps) as three A2 contract v2.1
modules in JAX + diffrax (``trunk.build(lane="trunk")``, ``trunk_lane.build``, ``trunk_lane.build_fed_3deoxy``;
plan step A2 + H3; 2026-11-05): do they reproduce the engine's trunk lane, are their exact Jacobians right, what do
one solve and its Jacobian cost, and what do the atol and Kvaerno4 levers buy?

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml), never in the core env:

    source <conda.sh>; conda activate maillard-jax
    XLA_FLAGS=--xla_cpu_multi_thread_eigen=false OMP_NUM_THREADS=1 nice -n 10 \
        python scripts/generators/build_jax_trunk_lane_module.py [--stage rates ...]

Every compared engine number comes from results/validation/jax_trunk_lane_module_reference.json (core env,
``build_jax_trunk_lane_reference.py``); the engine is called in this process only to BUILD the modules (their static
data) and, back to back with the JAX runs, for the paired timings. Writes results/validation/jax_trunk_lane_module.
{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

REFERENCE = ROOT / "results" / "validation" / "jax_trunk_lane_module_reference.json"
OUT = ROOT / "results" / "validation" / "jax_trunk_lane_module"
#: name -> (rtol, atol); "fit_tol_1e-8_atol_1e-10" is the atol lever the acrylamide port found (its finding 4)
TOLERANCES = {"tight_1e-10": (1e-10, 1e-16), "fit_tol_1e-8": (1e-8, 1e-14), "fit_tol_1e-8_atol_1e-10": (1e-8, 1e-10),
              "loose_1e-6": (1e-6, 1e-12)}
DEX_FLOOR = 1e-9
#: the named products: the trunk lane's panel targets and each block's end products
NAMED = ("HMF", "TDG", "DDG", "G", "GO", "MGO", "DMHF", "PZ", "DMP", "MPZ", "CML", "CEL", "MTAL", "MSH", "DMDS", "AP",
         "ATHP", "DGAL")
#: the structural part of a record's Config (what changes the compiled network); the rest is theta
STRUCTURAL = ("methionine_b22", "proline_b24", "exit_term")
TIMING_POTS = (
    ("isothermal hold, glucose only (18 min at 121 C)", "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021",
     "point"),
    ("fed dicarbonyls + glycine (pyrazine; 60 min at 120 C)", "synthetic_fed_glyoxal_methylglyoxal_glycine_120C",
     "point"),
    ("glucose + bound lysine (glycation; 40 min at 120 C)", "synthetic_glucose_bound_lysine_120C", "point"),
    ("fed 3-deoxyglucosone (B39 triangle; 60 min at 120 C)", "synthetic_fed_3deoxyglucosone_120C", "point"),
    ("staircase profile, every block", "synthetic_profile_trunk_lane_all_blocks_25_160_100", "point"),
)
#: (label, solver, rtol, atol) for the timing table: the default, the two levers alone and together, and rtol 1e-6
TIMING_SETTINGS = (("K5 rtol 1e-8 atol 1e-14", "Kvaerno5", 1e-8, 1e-14), ("K5 rtol 1e-8 atol 1e-10", "Kvaerno5", 1e-8, 1e-10),
                   ("K4 rtol 1e-8 atol 1e-14", "Kvaerno4", 1e-8, 1e-14), ("K4 rtol 1e-8 atol 1e-10", "Kvaerno4", 1e-8, 1e-10),
                   ("K4 rtol 1e-6 atol 1e-10", "Kvaerno4", 1e-6, 1e-10))


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


_NETS: Dict[Tuple, Any] = {}


def network_for(spec: Optional[Mapping[str, Any]] = None):
    """The composed trunk lane for a record's Config, cached on its STRUCTURAL part (which steps are zero, whether
    the exit term is on): the exit exponent and the aqueous-glyoxal values are theta, so they share the compile."""
    import build_jax_trunk_lane_reference as refgen
    from src.kinetic_core_jax import trunk_lane

    spec = dict(spec or {})
    key = tuple((k, spec.get(k)) for k in STRUCTURAL)
    if key not in _NETS:
        _NETS[key] = trunk_lane.compose_lane(refgen.config_for({k: v for k, v in key if v is not None},
                                                               engine_side=False))
    return _NETS[key]


def setup():
    from src.kinetic_core import network as engine_network

    reference = json.loads(REFERENCE.read_text())
    net = network_for()
    if list(net.reactions) != reference["reactions"] or tuple(net.reactions) != engine_network.TRUNK_REACTION_KEYS:
        raise AssertionError("reaction order differs from the engine's network.TRUNK_REACTIONS")
    return reference, net


def engine_columns_of_species(reference, net) -> List[int]:
    return [reference["engine_species"].index(s) for s in net.species]


def theta_of(rec: Mapping[str, Any], net):
    """The network theta for a record: k_ref / Ea per engine parameter (log10 or linear, as the module declares),
    the condition coordinates, and the exit exponent from the record's Config spec."""
    import jax.numpy as jnp

    theta = {m.name: {k: float(v) for k, v in m.params.items()} for m in net.modules}
    owner = {k: m for m, k in net.coordinates}
    for key, (k_ref, ea) in rec["parameters"].items():
        if f"log10_k:{key}" in owner:
            if not k_ref > 0.0:
                raise AssertionError(f"{key}: k_ref {k_ref} on a log10 coordinate")
            theta[owner[f"log10_k:{key}"]][f"log10_k:{key}"] = math.log10(k_ref)
        else:
            theta[owner[f"k_ref:{key}"]][f"k_ref:{key}"] = float(k_ref)
        theta[owner[f"ea:{key}"]][f"ea:{key}"] = float(ea)
    for key, value in rec["conditions"].items():
        theta[owner[key]][key] = float(value)
    spec = rec.get("config_spec") or {}
    if spec.get("exit_exponent") is not None and spec.get("exit_term") is not False:
        theta[owner["ph_exponent_exit:k_tdg_fa"]]["ph_exponent_exit:k_tdg_fa"] = float(spec["exit_exponent"])
    return {m: {k: jnp.asarray(v) for k, v in d.items()} for m, d in theta.items()}


def pot_of(rec: Mapping[str, Any], net, n_segments=None, segments_key="refined_segments") -> Dict[str, Any]:
    from src.kinetic_core_jax import trunk_lane

    lysp = rec["engine_initial"].get("LYSP") if "LYSP" not in rec["initial"] else None
    pot = trunk_lane.make_pot(net, rec["initial"], rec[segments_key], ph=rec["ph"],
                              water_activity=rec["water_activity"], n_segments=n_segments,
                              atmosphere=rec.get("atmosphere"), lysp=lysp)
    expect = np.zeros(net.n_species)
    for k, v in rec["engine_initial"].items():
        expect[net.index[k]] = float(v)
    if not np.array_equal(pot["y0"], expect):
        raise AssertionError(f"{rec['id']}: make_pot's charge differs from the engine's prelude")
    return pot


def compare(ys: np.ndarray, rec: Dict[str, Any], y0: np.ndarray, net, cols) -> Dict[str, Any]:
    eng = np.asarray(rec["states"])[:, cols]
    ys = ys[: len(eng)]
    scale = float(np.max(y0))
    above = eng > DEX_FLOOR * scale
    dex = np.where(above, np.abs(np.log10(np.where(above, np.maximum(ys, 1e-300), 1.0))
                                 - np.log10(np.where(above, eng, 1.0))), 0.0)
    r, i = np.unravel_index(int(np.argmax(dex)), dex.shape)
    below = np.where(~above, np.abs(ys - eng), 0.0) / scale
    named = [list(net.species).index(s) for s in NAMED]
    return {"max_dex": float(dex.max()), "worst_species": net.species[i], "worst_row": int(r),
            "max_dex_named_products": float(dex[:, named].max()),
            "max_abs_below_floor_share_of_y0": float(below.max())}


def group_of(rec) -> str:
    return {"synthetic": "synthetic isothermal", "synthetic_profile": "profile"}.get(rec["panel"], "panel")


# --------------------------------------------------------------------------- stages
def _rate_errors(k: np.ndarray, expect: np.ndarray):
    rel = np.abs(k - expect) / np.maximum(np.abs(expect), 1e-300)
    return np.where((k == 0.0) & (expect == 0.0), 0.0, rel)


def stage_rates(reference, net) -> Dict[str, Any]:
    """All 55 rate constants against trunk_conditions.apply + network.rate_constants_at on the (T, pH, a_w) grid."""
    import jax
    import jax.numpy as jnp
    import build_jax_trunk_lane_reference as refgen
    from src.kinetic_core_jax import network, trunk, trunk_lane

    probe = reference["rate_probe"]
    pts = [(ph, aw, t) for ph in probe["ph"] for aw in probe["aw"] for t in probe["temperatures_c"]]
    t_c = jnp.asarray([p[2] for p in pts])
    ph = jnp.asarray([trunk.nan_or(p[0]) for p in pts])
    aw = jnp.asarray([trunk.nan_or(p[1]) for p in pts])
    worst, worst_at, n, cases = 0.0, None, 0, {}
    for case in probe["cases"]:
        spec = case["config_spec"]
        cfg = refgen.config_for(spec, engine_side=False)
        if cfg.digest() != case["jax_config_digest"]:
            raise AssertionError(f"{case['label']}: the rebuilt Config differs from the reference's")
        net_c = network_for(spec)
        theta = theta_of(case, net_c)
        # the module built from the FULL Config carries the same values theta does for the coordinates only the
        # Config moves (the aqueous k_ama_g, the exit exponent; k_g_go is also moved by the probe's PERTURB)
        full = {m.name: m for m in trunk_lane.compose_lane(cfg).modules}
        for m, k in net_c.coordinates:
            if k in ("log10_k:k_ama_g", "ph_exponent_exit:k_tdg_fa") and "k_ama_g" not in refgen.PERTURB:
                if abs(float(full[m].params[k]) - float(theta[m][k])) > 1e-12:
                    raise AssertionError(f"{case['label']}: {k} from the Config {full[m].params[k]} vs the "
                                         f"engine's {float(theta[m][k])}")
        f = jax.jit(jax.vmap(lambda t, p, a: net_c.rate_constants(theta, network.RateCtx(
            t_c=t, t_k=t + network.CELSIUS, ph_label=p, aw=a))))
        k = np.asarray(f(t_c, ph, aw))
        rel = _rate_errors(k, np.asarray(case["values"]))
        n += rel.size
        i, j = np.unravel_index(int(np.argmax(rel)), rel.shape)
        cases[case["label"]] = {"worst_relative": float(rel[i, j]), "at": tuple(pts[i]) + (net_c.reactions[j],),
                                "n_values": int(rel.size), "n_nonzero": int(np.sum(np.asarray(case["values"]) != 0))}
        if rel[i, j] > worst:
            worst, worst_at = float(rel[i, j]), (case["label"],) + tuple(pts[i]) + (net_c.reactions[j],)
    out = {"n_rate_values": n, "n_grid_points": len(pts), "worst_relative": worst, "worst_at": worst_at,
           "cases": cases}
    print(f"  rates: {n} values, worst {worst:.1e} at {worst_at}")
    return out


def stage_trajectories(reference, net) -> Dict[str, Any]:
    import diffrax
    import jax
    from src.kinetic_core_jax import network

    cols = engine_columns_of_species(reference, net)
    records = reference["pots"]
    s_iso = max(len(r["refined_segments"]) for r in records if group_of(r) != "profile")
    rows = []
    compiled: Dict[Tuple, Any] = {}
    for solver_name in ("Kvaerno4", "Kvaerno5"):
        solver = getattr(diffrax, solver_name)()
        for tol_name, (rtol, atol) in TOLERANCES.items():
            for group in ("panel", "synthetic isothermal", "profile"):
                recs = [r for r in records if group_of(r) == group]
                worst = {"max_dex": 0.0}
                steps, below, secs, named, closure = 0, 0.0, 0.0, 0.0, 0.0
                for rec in recs:
                    net_r = network_for(rec["config_spec"])
                    key = (id(net_r), solver_name, tol_name)
                    if key not in compiled:
                        compiled[key] = jax.jit(lambda th, pot, net_r=net_r, solver=solver, rtol=rtol, atol=atol:
                                                network.integrate(net_r, th, pot, solver=solver, rtol=rtol,
                                                                  atol=atol))
                    run = compiled[key]
                    pot = pot_of(rec, net_r, None if group == "profile" else s_iso)
                    th = theta_of(rec, net_r)
                    out = run(th, pot)
                    out.ys.block_until_ready()
                    t0 = time.perf_counter()
                    out = run(th, pot)
                    out.ys.block_until_ready()
                    secs += time.perf_counter() - t0
                    ys = np.asarray(out.ys)
                    c = compare(ys, rec, pot["y0"], net_r, cols)
                    c["pot"] = f"{rec['id']}:{rec['variant']}"
                    steps += int(np.sum(out.num_steps))
                    below = max(below, c["max_abs_below_floor_share_of_y0"])
                    named = max(named, c["max_dex_named_products"])
                    e0 = net_r.element_matrix @ pot["y0"]
                    closure = max(closure, float(np.max(np.abs(net_r.element_matrix @ ys[-1] - e0)
                                                        / np.maximum(e0, 1e-300) * (e0 > 0))))
                    if c["max_dex"] >= worst["max_dex"]:
                        worst = c
                rows.append({"solver": solver_name, "tolerance": tol_name, "group": group, "n_pots": len(recs),
                             "max_dex": worst["max_dex"], "worst": worst, "steps_total": steps,
                             "max_dex_named_products": named, "max_rel_element_drift": closure,
                             "max_abs_below_floor_share_of_y0": below, "wall_s_all_pots_after_jit": secs})
                print(f"  {solver_name} {tol_name} {group}: max dex {worst['max_dex']:.2e} ({worst.get('pot')}, "
                      f"{worst.get('worst_species')}), named {named:.1e}, {steps} steps, {secs:.3f} s")
    return {"dex_floor_share_of_max_initial": DEX_FLOOR, "named_products": list(NAMED), "rows": rows}


RELATIVE_COLUMN_SHARE = 1e-4


def column_errors(a: np.ndarray, b: np.ndarray, share: float = RELATIVE_COLUMN_SHARE) -> Dict[str, Any]:
    """Per column: max |a - b| / max |b| where max |b| >= share x max |b| overall, else max |a - b| (absolute);
    plus the worst absolute difference over all columns."""
    flat = share * float(np.max(np.abs(b)))
    cols, absolute = [], []
    for c in range(b.shape[1]):
        scale = float(np.max(np.abs(b[:, c])))
        diff = float(np.max(np.abs(a[:, c] - b[:, c])))
        cols.append(diff / scale if scale >= flat else diff)
        absolute.append(diff)
    flat_cols = [c for c in range(b.shape[1]) if float(np.max(np.abs(b[:, c]))) < flat]
    worst = int(np.argmax(cols))
    return {"per_column": cols, "worst": float(cols[worst]), "worst_column": worst,
            "worst_absolute": float(max(absolute)), "flat_columns_absolute": flat_cols,
            "n_relative_columns": b.shape[1] - len(flat_cols)}


def stage_jacobians(reference, net) -> Dict[str, Any]:
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import network

    out = {}
    for pot_id, fd in reference["jacobian_central_fd"].items():
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
        net_r = network_for(fd["config_spec"])
        keys = [k for _, k in net_r.coordinates]
        obs_idx = [net_r.index[s] for s in fd["species"]]
        rows = fd["rows"]
        pot = pot_of(rec, net_r)

        def observe(ys, rows=rows, obs_idx=obs_idx):
            return jnp.log10(jnp.maximum(ys[jnp.asarray(rows)][:, jnp.asarray(obs_idx)].T.ravel(), 1e-300))

        theta = network.theta_vector(net_r, theta_of(fd, net_r))
        res, secs = {}, {}
        for tol_name in ("tight_1e-10", "fit_tol_1e-8", "fit_tol_1e-8_atol_1e-10", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            t0 = time.perf_counter()
            _, jac = network.jacobian(net_r, observe, rtol=rtol, atol=atol)(theta, pot)
            res[tol_name] = np.asarray(jac)
            secs[tol_name] = time.perf_counter() - t0
        f = jax.jit(lambda v: observe(network.integrate(net_r, network.theta_from_vector(net_r, v), pot,
                                                        rtol=1e-12, atol=1e-20).ys))
        h = fd["step"]
        linear = [i for i, k in enumerate(keys) if k.startswith("k_ref:")]
        #: a linear k_ref coordinate AT ZERO is differenced one-sidedly (second order: (-3 f0 + 4 f(h) - f(2h)) / 2h):
        #: a central difference's minus side is a NEGATIVE rate constant, which drives its products below zero, into
        #: the max(y, 0) of the mass-action law, so it straddles a kink the physical (k >= 0) derivative does not see
        one_sided = [i for i in linear if float(theta[i]) == 0.0]
        f0 = np.asarray(f(theta))
        jax_fd = []
        for c in range(len(theta)):
            e = jnp.zeros(len(theta)).at[c].set(h)
            if c in one_sided:
                jax_fd.append((-3.0 * f0 + 4.0 * np.asarray(f(theta + e)) - np.asarray(f(theta + 2 * e))) / (2 * h))
            else:
                jax_fd.append((np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h))
        jax_fd = np.array(jax_fd).T
        cols = [keys.index(c) for c in fd["columns"]]
        eng = np.asarray(fd["jacobian"])
        r = {
            "rows": rows, "species": fd["species"], "n_observables": int(eng.shape[0]),
            "n_coordinates_engine": len(cols), "n_coordinates": len(keys), "n_linear_k_ref": len(linear),
            "one_sided_columns": [keys[i] for i in one_sided],
            "forward_tight_vs_engine_central_fd": column_errors(res["tight_1e-10"][:, cols], eng),
            "forward_tight_vs_jax_central_fd": column_errors(res["tight_1e-10"], jax_fd),
            "jax_central_fd_vs_engine_central_fd": column_errors(jax_fd[:, cols], eng),
            "forward_fit_tol_vs_forward_tight": column_errors(res["fit_tol_1e-8"], res["tight_1e-10"]),
            "forward_fit_tol_atol_1e-10_vs_forward_tight": column_errors(res["fit_tol_1e-8_atol_1e-10"],
                                                                        res["tight_1e-10"]),
            "forward_loose_vs_forward_tight": column_errors(res["loose_1e-6"], res["tight_1e-10"]),
            "jacobian_wall_s_including_compile": secs,
        }
        for k in ("forward_tight_vs_engine_central_fd", "jax_central_fd_vs_engine_central_fd"):
            r[k]["worst_coordinate"] = fd["columns"][r[k]["worst_column"]]
        r["forward_tight_vs_jax_central_fd"]["worst_coordinate"] = keys[
            r["forward_tight_vs_jax_central_fd"]["worst_column"]]
        out[pot_id] = r
        print(f"  {pot_id}: forward vs engine FD {r['forward_tight_vs_engine_central_fd']['worst']:.2e} "
              f"({r['forward_tight_vs_engine_central_fd']['worst_coordinate']}), vs own FD "
              f"{r['forward_tight_vs_jax_central_fd']['worst']:.2e} "
              f"({r['forward_tight_vs_jax_central_fd']['worst_coordinate']})")
    return out


def stage_residuals(reference, net) -> Dict[str, Any]:
    """The panel's predictions (each scored row's log10 value, after the engine's post-integration operators:
    carried volatiles added, matrix binding applied) at fit tolerance against engine.predict's."""
    import diffrax
    import jax
    from src.kinetic_core_jax import network

    preds = {p["id"]: p for p in reference["panel_predictions"]}
    rows = []
    for solver_name in ("Kvaerno4", "Kvaerno5"):
        solver = getattr(diffrax, solver_name)()
        for tol_name in ("fit_tol_1e-8", "fit_tol_1e-8_atol_1e-10", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            run = jax.jit(lambda th, pot: network.integrate(net, th, pot, solver=solver, rtol=rtol, atol=atol).ys)
            worst, items = 0.0, []
            for rec in (r for r in reference["pots"] if r["variant"] == "point" and group_of(r) == "panel"):
                ys = np.asarray(run(theta_of(rec, net), pot_of(rec, net, segments_key="segments")))
                for row in preds[rec["id"]]["rows"]:
                    if row["predicted"] is None:
                        continue
                    name = f"ug_per_L:{row['species']}"
                    raw = float(net.observe(ys[-1], [name])[0])
                    mine = (raw + row["carried_added"]) * row["binding_remaining_fraction"]
                    d = abs(math.log10(mine) - math.log10(row["predicted"]))
                    items.append({"pot": rec["id"], "compound": row["compound"], "engine": row["predicted"],
                                  "jax": mine, "measured": row["measured"], "engine_residual_dex": row["residual_dex"],
                                  "abs_diff_dex": d})
                    worst = max(worst, d)
            rows.append({"solver": solver_name, "tolerance": tol_name, "n_rows": len(items), "max_abs_diff_dex": worst,
                         "rows": items})
            print(f"  residuals {solver_name} {tol_name}: {len(items)} rows, max |diff| {worst:.1e} dex")
    return {"rows": rows}


def stage_timing(reference, net) -> Dict[str, Any]:
    """Engine and JAX back to back on the same pots (the engine in this process for the timing only): one lane run
    (``engine._integrate_program``, its rtol 1e-8) and its forward-difference Jacobian in all module coordinates
    ((1 + n) runs, from the median run), against one JAX solve and one exact Jacobian, at each lever setting."""
    import diffrax
    import jax
    import jax.numpy as jnp
    import build_jax_trunk_lane_reference as refgen
    from src.kinetic_core import config, engine
    from src.kinetic_core_jax import network

    panel_by_id = {p["id"]: p for p in refgen.trunk_panel_pots()}
    obs_idx = jnp.asarray([net.index[s] for s in ("HMF", "MGO", "GO", "TDG")])
    n_coord = len(net.coordinates)
    rows = []
    for label, pot_id, variant in TIMING_POTS:
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == variant)
        if pot_id in panel_by_id:
            process = panel_by_id[pot_id]["spec"].process
        elif pot_id == refgen.PROFILE_ID:
            process = refgen.synthetic_process(refgen.PROFILE_POINTS, 5.5, 0.80)
        else:
            _, seg, ph, aw, _ = refgen.SYNTHETIC_ISOTHERMAL[pot_id]
            process = refgen.synthetic_process(seg, ph, aw)
        params = engine.core_parameters(engine.TRUNK)
        with config.applied(config.Config.from_globals()):
            engine._integrate_program(engine.TRUNK, params, dict(rec["initial"]), process)
            secs = []
            for _ in range(5 if len(rec["segments"]) > 1 else 10):
                t0 = time.perf_counter()
                engine._integrate_program(engine.TRUNK, params, dict(rec["initial"]), process)
                secs.append(time.perf_counter() - t0)
        med = float(np.median(secs))
        out = {"class": label, "pot": pot_id, "n_segments": len(rec["segments"]),
               "engine": {"lane_run_median_s": med, "forward_fd_jacobian_s": (1 + n_coord) * med,
                          "fd_runs": 1 + n_coord}, "jax": {}}
        pot = pot_of(rec, net, segments_key="segments")

        def observe(ys):
            return jnp.log10(jnp.maximum(ys[-1][obs_idx], 1e-300))

        theta = network.theta_vector(net, net.theta0())
        for name, solver_name, rtol, atol in TIMING_SETTINGS:
            solver = getattr(diffrax, solver_name)()
            solve = jax.jit(lambda v, p: network.integrate(net, network.theta_from_vector(net, v), p, solver=solver,
                                                           rtol=rtol, atol=atol))
            t0 = time.perf_counter()
            jax.block_until_ready(solve(theta, pot))
            compile_s = time.perf_counter() - t0
            js = []
            for _ in range(5):
                t0 = time.perf_counter()
                res = jax.block_until_ready(solve(theta, pot))
                js.append(time.perf_counter() - t0)
            jac_fn = network.jacobian(net, observe, solver=solver, rtol=rtol, atol=atol)
            t0 = time.perf_counter()
            jax.block_until_ready(jac_fn(theta, pot))
            jac_compile = time.perf_counter() - t0
            jj = []
            for _ in range(3):
                t0 = time.perf_counter()
                jax.block_until_ready(jac_fn(theta, pot))
                jj.append(time.perf_counter() - t0)
            out["jax"][name] = {"solve_median_s": float(np.median(js)), "solve_compile_s": compile_s,
                                "steps": int(np.sum(np.asarray(res.num_steps))),
                                "jacobian_median_s": float(np.median(jj)), "jacobian_compile_s": jac_compile,
                                "n_coordinates": n_coord}
        out["uptime"] = uptime()
        rows.append(out)
        d = out["jax"][TIMING_SETTINGS[0][0]]
        print(f"  timing {pot_id}: engine run {med:.2e} s, FD Jacobian {(1 + n_coord) * med:.2e} s; JAX solve "
              f"{d['solve_median_s']:.2e} s ({d['steps']} steps), Jacobian {d['jacobian_median_s']:.2e} s")
    return {"rows": rows, "settings": [s[0] for s in TIMING_SETTINGS]}


# --------------------------------------------------------------------------- the write-up
def write_markdown(report: Dict[str, Any]) -> str:
    rates, traj, jac, resid, timing = (report[k] for k in ("rates", "trajectories", "jacobians", "residuals",
                                                           "timing"))
    lines = [
        "# The trunk lane as A2 contract v2.1 modules in JAX + diffrax (plan step A2 + H3)",
        "",
        f"*Generated by `scripts/generators/build_jax_trunk_lane_module.py` in the separate `maillard-jax` env (jax "
        f"{report['jax_version']}, diffrax {report['diffrax_version']}, float64, CPU); the engine's numbers from "
        "`jax_trunk_lane_module_reference.json` (core env, `build_jax_trunk_lane_reference.py`). Nothing shipped "
        "moves; the port is not imported by `src/kinetic_core`. Contract: "
        "[A2_MODULE_CONTRACT.md](../../docs/guides/A2_MODULE_CONTRACT.md) v2.1. The trunk module and the acrylamide "
        "lane: [jax_trunk_acrylamide_module.md](jax_trunk_acrylamide_module.md).*",
        "",
        "## What was ported",
        "",
        f"- `src/kinetic_core_jax/trunk_lane.py`: the trunk lane's seven blocks as two modules composed after "
        f"`trunk.build(lane=\"trunk\")` -- `trunk_lane` (dicarbonyl 5, pyrazine 5, glycation 5, aqueous glyoxal 1, "
        f"methionine 6, proline 4: 26 steps, 23 species) and `fed_3deoxy` (3 steps, owns DGAL, also part of the "
        f"sulfur network). `compose_lane` gives `network.TRUNK_REACTIONS` in the engine's order: "
        f"{report['n_species']} species, {report['n_reactions']} reactions, {len(report['coordinates'])} "
        "coordinates.",
        "- Every switch the lane reads, ported or refused (finding 5). The trunk module's own terms (a_w on "
        "`k_schiff`, the Amadori pH term, the 3-DG exit term), ported at rate level by the previous port, are now "
        "checked on trajectories.",
        "",
        "## Rates (no integration)",
        "",
        f"- All 55 steps, {rates['n_rate_values']} values ({rates['n_grid_points']} (T, pH, a_w) points x 55 x "
        f"{len(rates['cases'])} cases: shipped; perturbed parameters, conditions, aqueous-glyoxal mapping and exit "
        "exponent (Config), B22 / B24 on; exit term off): worst relative difference "
        f"**{rates['worst_relative']:.1e}** (at {rates['worst_at']}).",
        "",
        "| case | values | non-zero | worst relative | at |", "|---|---|---|---|---|",
    ]
    for label, c in rates["cases"].items():
        lines.append(f"| {label} | {c['n_values']} | {c['n_nonzero']} | {c['worst_relative']:.1e} | {c['at']} |")
    lines += ["", "## Trajectories (state after every refined segment vs the engine at LSODA rtol 1e-12; dex where "
              f"the engine holds more than {traj['dex_floor_share_of_max_initial']:g} of the pot's largest initial "
              "charge)", "",
              "| solver | tolerance | pots | max dex, all species | worst (pot, species) | max dex, named products | steps (all pots) | wall s after JIT | max abs below floor / max y0 | element drift |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for r in traj["rows"]:
        w = r["worst"]
        lines.append(f"| {r['solver']} | {r['tolerance']} | {r['group']} ({r['n_pots']}) | {r['max_dex']:.1e} | "
                     f"{w.get('pot')} {w.get('worst_species')} | {r['max_dex_named_products']:.1e} | "
                     f"{r['steps_total']} | {r['wall_s_all_pots_after_jit']:.3f} | "
                     f"{r['max_abs_below_floor_share_of_y0']:.1e} | {r['max_rel_element_drift']:.1e} |")
    lines += ["", "Pots: the one panel bundle the engine routes to the trunk lane (Steinhagen 2021, glucose only, "
              "18 min at 121 C, pH 4.36) at its point and both band corners (a_w scale, Amadori pH exponent, 3-DG exit "
              "exponent, the disputed and HPD sinks, all at one end); four synthetic isothermal pots, one per block "
              "the panel pot leaves idle (fed dicarbonyls + glycine; glucose + bound lysine; fed 3-DG; glucose + "
              "methionine + proline with B22 / B24 on and the exit term off); and a synthetic staircase profile "
              "(25 -> 160 -> 100 C, 2 K steps, pH 5.5, a_w 0.80, every block charged) at the shipped point and with "
              f"the switches on. Named products: {', '.join(traj['named_products'])}.", "",
              "## Exact Jacobian (forward mode) vs central differences, at an interior point", "",
              "Worst column error = max over observables of |difference| / max |reference| in that column, on the "
              f"columns whose largest entry is at least {RELATIVE_COLUMN_SHARE:g} of the Jacobian's largest (\"rel\"); "
              "the worst absolute difference over ALL columns, in dex per coordinate unit, beside it (\"abs\"). "
              "The engine moves every log10 k_ref, Ea and condition coordinate through its arguments (the "
              "exit exponent through the Config); the linear k_ref coordinates of the zero-k steps are checked "
              "against the port's own differences only, ONE-SIDED at zero (a negative rate constant would push its "
              "products into the mass-action law's max(y, 0) clip, a kink the physical derivative does not see).", "",
              "| pot | observables x coordinates (engine / all; relative columns) | vs engine central FD: rel (worst coordinate) / abs | vs JAX central FD: rel (worst) / abs | JAX FD vs engine FD: rel / abs | rtol 1e-8 vs tight | rtol 1e-8 atol 1e-10 vs tight | rtol 1e-6 vs tight |",
              "|---|---|---|---|---|---|---|---|"]
    for pot_id, r in jac.items():
        e, o, x = (r[k] for k in ("forward_tight_vs_engine_central_fd", "forward_tight_vs_jax_central_fd",
                                  "jax_central_fd_vs_engine_central_fd"))
        lines.append(f"| {pot_id} | {r['n_observables']} x {r['n_coordinates_engine']} / {r['n_coordinates']} "
                     f"({e['n_relative_columns']}) | {e['worst']:.1e} ({e['worst_coordinate']}) / "
                     f"{e['worst_absolute']:.1e} | {o['worst']:.1e} ({o['worst_coordinate']}) / "
                     f"{o['worst_absolute']:.1e} | {x['worst']:.1e} / {x['worst_absolute']:.1e} | "
                     f"{r['forward_fit_tol_vs_forward_tight']['worst']:.1e} | "
                     f"{r['forward_fit_tol_atol_1e-10_vs_forward_tight']['worst']:.1e} | "
                     f"{r['forward_loose_vs_forward_tight']['worst']:.1e} |")
    lines += ["", "## Panel predictions at fit tolerance vs `engine.predict` (log10 of each scored row's value)", "",
              "| solver | tolerance | rows | max abs difference, dex |", "|---|---|---|---|"]
    for r in resid["rows"]:
        lines.append(f"| {r['solver']} | {r['tolerance']} | {r['n_rows']} | {r['max_abs_diff_dex']:.1e} |")
    lines += ["", "## Timing (CPU, one core, niced; shared machine: load beside each row; provisional)", "",
              f"The engine's forward-difference Jacobian is (1 + {len(report['coordinates'])}) lane runs, from the "
              "median run; the JAX Jacobian is one forward-mode solve over all coordinates. Seconds per call after "
              "compilation; steps of one solve.", ""]
    head = "| pot class | segments | engine run s | engine FD Jacobian s | " + " | ".join(
        f"{s}: solve s / steps / Jacobian s" for s in timing["settings"]) + " | compile s (solve / Jacobian, first setting) | load |"
    lines += [head, "|" + "---|" * (6 + len(timing["settings"]))]
    for r in timing["rows"]:
        e = r["engine"]
        load = r["uptime"].split("load averages:")[-1].strip() if "load" in r["uptime"] else r["uptime"]
        cells = [f"{r['jax'][s]['solve_median_s']:.1e} / {r['jax'][s]['steps']} / {r['jax'][s]['jacobian_median_s']:.1e}"
                 for s in timing["settings"]]
        first = r["jax"][timing["settings"][0]]
        lines.append(f"| {r['class']} | {r['n_segments']} | {e['lane_run_median_s']:.1e} | "
                     f"{e['forward_fd_jacobian_s']:.1e} | " + " | ".join(cells)
                     + f" | {first['solve_compile_s']:.1f} / {first['jacobian_compile_s']:.1f} | {load} |")
    lines += [""]
    lines += report.get("findings", [])
    return "\n".join(lines) + "\n"


STAGES = ("rates", "trajectories", "jacobians", "residuals", "timing")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", action="append", choices=STAGES)
    args = parser.parse_args(argv)
    import diffrax
    import jax

    reference, net = setup()
    stages = args.stage or list(STAGES)
    path = OUT.with_suffix(".json")
    report = json.loads(path.read_text()) if path.exists() and args.stage else {}
    report.update({"artifact": "jax_trunk_lane_module",
                   "generated_by": "scripts/generators/build_jax_trunk_lane_module.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "diffrax_version": diffrax.__version__,
                   "n_species": net.n_species, "n_reactions": len(net.reactions), "species": list(net.species),
                   "coordinates": [f"{m}.{k}" for m, k in net.coordinates], "reactions": list(net.reactions),
                   "modules": {m.name: {"reactions": len(m.reactions), "species": list(m.species),
                                        "coordinates": len(m.params), "notes": list(m.notes)} for m in net.modules},
                   "ph_mode": net.ph_mode, "config_digests": {m.name: m.config_digest for m in net.modules}})
    for stage in stages:
        print(f"[{stage}] {uptime()}")
        report[stage] = globals()[f"stage_{stage}"](reference, net)
        path.write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    report["findings"] = FINDINGS
    path.write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    if all(k in report for k in STAGES):
        OUT.with_suffix(".md").write_text(write_markdown(report), encoding="utf-8")
    print(f"wrote {path}")
    return 0


#: the write-up's closing section, edited by hand after reading the numbers (kept here so the .md regenerates)
FINDINGS: List[str] = [
    "## Findings",
    "",
    "(pending)",
    "",
]


if __name__ == "__main__":
    sys.exit(main())

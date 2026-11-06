#!/usr/bin/env python
"""
build_jax_trunk_acrylamide_module.py -- the trunk (``network.REACTIONS``) and the acrylamide block as modules of the
A2 contract v2, composed into the engine's acrylamide lane, in JAX + diffrax (plan step A2 + H3;
docs/guides/A2_MODULE_CONTRACT.md; 2026-11-05): do the modules reproduce the engine, are their exact Jacobians right,
what do one solve and its Jacobian cost?

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml), never in the core env:

    source <conda.sh>; conda activate maillard-jax
    XLA_FLAGS=--xla_cpu_multi_thread_eigen=false OMP_NUM_THREADS=1 nice -n 10 \
        python scripts/generators/build_jax_trunk_acrylamide_module.py

Every compared engine number comes from results/validation/jax_trunk_acrylamide_module_reference.json (core env,
``build_jax_trunk_acrylamide_reference.py``); the engine is called in this process only to BUILD the modules (their
static data) and, back to back with the JAX runs, for the paired timings. Writes
results/validation/jax_trunk_acrylamide_module.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

REFERENCE = ROOT / "results" / "validation" / "jax_trunk_acrylamide_module_reference.json"
OUT = ROOT / "results" / "validation" / "jax_trunk_acrylamide_module"
TOLERANCES = {"tight_1e-10": (1e-10, 1e-16), "fit_tol_1e-8": (1e-8, 1e-14), "loose_1e-6": (1e-6, 1e-12)}
#: a species is compared in dex where the engine holds more than this share of the pot's largest initial charge
DEX_FLOOR = 1e-9
#: the lane's scored products (what the panel's acrylamide-lane rows read)
OBSERVED = ("ACR", "HMF")
#: the observable each scored species is reported in (the engine's unit: ppb for acrylamide, ug/L otherwise)
OBSERVABLE_OF = {"ACR": "ppb:ACR", "HMF": "ug_per_L:HMF"}
#: 2026-11-06: the short-isothermal pot is the reference's SYNTHETIC numerics pot (the withdrawn SPI extrusion inputs).
TIMING_POTS = (("short isothermal (synthetic, 0.42 min at 130 C)", "synthetic_short_isothermal_acrylamide_130C_25s"),
               ("isothermal hold (30 min at 180 C)", "mp_holdout_glucose_asparagine_180C_30min_Chang2021"),
               ("dense isothermal (200 mM, 30 min at 180 C)", "mp_holdout_glucose_asparagine_180C_Ye2024"),
               ("staircase profile (all 42 steps)", "synthetic_profile_acrylamide_all_steps_25_180_100"))


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def setup():
    from src.kinetic_core_jax import acrylamide, network, trunk

    reference = json.loads(REFERENCE.read_text())
    net = network.compose([trunk.build(), acrylamide.build()])
    if list(net.reactions) != reference["reactions"]:
        raise AssertionError("reaction order differs from the engine's FULL_ACRYLAMIDE_REACTIONS")
    return reference, net


def engine_columns_of_species(reference, net) -> List[int]:
    """For each network species, its column in the reference's engine state."""
    return [reference["engine_species"].index(s) for s in net.species]


def pot_of(rec: Dict[str, Any], net, n_segments=None, segments_key="refined_segments") -> Dict[str, Any]:
    from src.kinetic_core_jax import acrylamide

    return acrylamide.make_pot(net, rec["initial"], rec[segments_key], ph=rec["ph"],
                               water_activity=rec["water_activity"], n_segments=n_segments)


def theta_of(theta_rec: Dict[str, Any], net):
    """The network theta for a reference record's ``theta`` (k_ref / Ea per engine parameter, condition coordinates)."""
    import jax.numpy as jnp

    theta = {m.name: {k: float(v) for k, v in m.params.items()} for m in net.modules}
    owner = {k: m for m, k in net.coordinates}
    for key, (k_ref, ea) in theta_rec["parameters"].items():
        theta[owner[f"log10_k:{key}"]][f"log10_k:{key}"] = math.log10(k_ref)
        theta[owner[f"ea:{key}"]][f"ea:{key}"] = ea
    for key, value in theta_rec["conditions"].items():
        theta[owner[key]][key] = float(value)
    return {m: {k: jnp.asarray(v) for k, v in d.items()} for m, d in theta.items()}


def compare(ys: np.ndarray, rec: Dict[str, Any], y0: np.ndarray, net, cols) -> Dict[str, Any]:
    """Worst dex over species above the floor, and the worst absolute difference (share of max y0) below it."""
    eng = np.asarray(rec["states"])[:, cols]
    ys = ys[: len(eng)]
    scale = float(np.max(y0))
    above = eng > DEX_FLOOR * scale
    dex = np.where(above, np.abs(np.log10(np.where(above, np.maximum(ys, 1e-300), 1.0))
                                 - np.log10(np.where(above, eng, 1.0))), 0.0)
    r, i = np.unravel_index(int(np.argmax(dex)), dex.shape)
    below = np.where(~above, np.abs(ys - eng), 0.0) / scale
    scored = [list(net.species).index(s) for s in OBSERVED]
    return {"max_dex": float(dex.max()), "worst_species": net.species[i], "worst_row": int(r),
            "max_dex_scored_products": float(dex[:, scored].max()),
            "max_abs_below_floor_share_of_y0": float(below.max())}


# --------------------------------------------------------------------------- stages
def _rate_errors(k: np.ndarray, expect: np.ndarray):
    rel = np.abs(k - expect) / np.maximum(np.abs(expect), 1e-300)
    return np.where((k == 0.0) & (expect == 0.0), 0.0, rel)


def stage_rates(reference, net) -> Dict[str, Any]:
    """The composed acrylamide network's rate vector against acrylamide_rate_constants_at(apply(...)) on the (T, pH,
    a_w) grid, and the trunk-lane module's against trunk_conditions.apply + network.rate_constants_at."""
    import jax
    import jax.numpy as jnp
    from src.kinetic_core import config
    from src.kinetic_core_jax import network, trunk

    def grid_points(probe):
        return [(ph, aw, t) for ph in probe["ph"] for aw in probe["aw"] for t in probe["temperatures_c"]]

    def run_grid(rate_fn, theta, probe):
        pts = grid_points(probe)
        t_c = jnp.asarray([p[2] for p in pts])
        ph = jnp.asarray([trunk.nan_or(p[0]) for p in pts])
        aw = jnp.asarray([trunk.nan_or(p[1]) for p in pts])
        f = jax.jit(jax.vmap(lambda t, p, a: rate_fn(theta, network.RateCtx(t_c=t, t_k=t + network.CELSIUS,
                                                                           ph_label=p, aw=a))))
        return np.asarray(f(t_c, ph, aw)), pts

    out: Dict[str, Any] = {}
    probe = reference["acrylamide_rate_probe"]
    worst, worst_at, n = 0.0, None, 0
    for case in probe["cases"]:
        k, pts = run_grid(net.rate_constants, theta_of(case["theta"], net), probe)
        rel = _rate_errors(k, np.asarray(case["k"]))
        n += rel.size
        i, j = np.unravel_index(int(np.argmax(rel)), rel.shape)
        if rel[i, j] > worst:
            worst, worst_at = float(rel[i, j]), (case["label"],) + tuple(pts[i]) + (net.reactions[j],)
    out["acrylamide_lane"] = {"n_rate_values": n, "worst_relative": worst, "worst_at": worst_at,
                              "n_grid_points": len(grid_points(probe)), "cases": [c["label"] for c in probe["cases"]]}
    probe = reference["trunk_lane_rate_probe"]
    base = config.Config.from_globals()
    worst, worst_at, n, digests = 0.0, None, 0, {}
    for case in probe["cases"]:
        exit_cfg = dict(base["trunk.three_deoxy_exit_ph"])
        changes: Dict[str, Any] = {"trunk.three_deoxy_exit_ph_term": case["exit_term"]}
        changes["trunk.three_deoxy_exit_ph"] = {k: (case["exit_ph"][k],) + tuple(v[1:]) for k, v in exit_cfg.items()}
        cfg = base.replace(changes)
        if cfg.digest() != case["config_digest"]:
            raise AssertionError(f"{case['label']}: the rebuilt Config differs from the reference's")
        mod = trunk.build(cfg, lane="trunk")
        theta = {k: float(v) for k, v in mod.params.items()}
        for key, (k_ref, ea) in case["parameters"].items():
            theta[f"log10_k:{key}"], theta[f"ea:{key}"] = math.log10(k_ref), ea
        theta.update(aw_scale=case["aw_scale"], ph_exponent=case["ph_exponent"])
        for key, e in case["exit_ph"].items():
            if case["exit_term"]:
                theta[f"ph_exponent_exit:{key}"] = e
        theta = {k: jnp.asarray(v) for k, v in theta.items()}
        k, pts = run_grid(mod.rate_constants, theta, probe)
        rel = _rate_errors(k, np.asarray(case["values"]))
        n += rel.size
        digests[case["label"]] = mod.config_digest
        i, j = np.unravel_index(int(np.argmax(rel)), rel.shape)
        if rel[i, j] > worst:
            worst, worst_at = float(rel[i, j]), (case["label"],) + tuple(pts[i]) + (mod.reactions[j],)
    out["trunk_lane_terms"] = {"n_rate_values": n, "worst_relative": worst, "worst_at": worst_at,
                               "cases": [c["label"] for c in probe["cases"]], "module_digests": digests}
    print(f"  rates: acrylamide lane {out['acrylamide_lane']['worst_relative']:.1e}, trunk-lane terms {worst:.1e}")
    return out


def stage_trajectories(reference, net) -> Dict[str, Any]:
    import diffrax
    import jax
    from src.kinetic_core_jax import network

    cols = engine_columns_of_species(reference, net)
    records = reference["pots"]
    panel_recs = [r for r in records if r["panel"] != "synthetic"]
    profile = [r for r in records if r["panel"] == "synthetic"]
    s_panel = max(len(r["refined_segments"]) for r in panel_recs)
    rows = []
    for solver_name, solver in (("Kvaerno4", diffrax.Kvaerno4()), ("Kvaerno5", diffrax.Kvaerno5())):
        for tol_name, (rtol, atol) in TOLERANCES.items():
            run = jax.jit(lambda th, pot: network.integrate(net, th, pot, solver=solver, rtol=rtol, atol=atol))
            for group, recs, s in (("panel", panel_recs, s_panel), ("profile", profile, None)):
                worst = {"max_dex": 0.0}
                steps, below, secs, scored, closure = 0, 0.0, 0.0, 0.0, 0.0
                for rec in recs:
                    pot = pot_of(rec, net, s)
                    th = theta_of(rec["theta"], net)
                    out = run(th, pot)
                    out.ys.block_until_ready()
                    t0 = time.perf_counter()
                    out = run(th, pot)
                    out.ys.block_until_ready()
                    secs += time.perf_counter() - t0
                    ys = np.asarray(out.ys)
                    c = compare(ys, rec, pot["y0"], net, cols)
                    c["pot"] = f"{rec['id']}:{rec['variant']}"
                    steps += int(np.sum(out.num_steps))
                    below = max(below, c["max_abs_below_floor_share_of_y0"])
                    scored = max(scored, c["max_dex_scored_products"])
                    e0 = net.element_matrix @ pot["y0"]
                    closure = max(closure, float(np.max(np.abs(net.element_matrix @ ys[-1] - e0)
                                                        / np.maximum(e0, 1e-300) * (e0 > 0))))
                    if c["max_dex"] >= worst["max_dex"]:
                        worst = c
                rows.append({"solver": solver_name, "tolerance": tol_name, "group": group, "n_pots": len(recs),
                             "max_dex": worst["max_dex"], "worst": worst, "steps_total": steps,
                             "max_dex_scored_products": scored, "max_rel_element_drift": closure,
                             "max_abs_below_floor_share_of_y0": below, "wall_s_all_pots_after_jit": secs})
                print(f"  {solver_name} {tol_name} {group}: max dex {worst['max_dex']:.2e} ({worst.get('pot')}, "
                      f"{worst.get('worst_species')}), scored {scored:.1e}, {steps} steps, {secs:.3f} s")
    return {"dex_floor_share_of_max_initial": DEX_FLOOR, "rows": rows}


#: a column is compared RELATIVELY where its largest entry is at least this share of the whole Jacobian's largest;
#: below it (|d log10 obs / d theta| ~ 1e-6, where a finite difference at rtol 1e-12 and step 1e-5 is noise-limited)
#: the column is compared ABSOLUTELY, and the worst absolute difference over every column is reported too
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

    obs_idx = [net.index[s] for s in reference["observable_species"]]
    keys = [k for _, k in net.coordinates]
    out = {}
    for pot_id, fd in reference["jacobian_central_fd"].items():
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
        rows = fd["rows"]
        pot = pot_of(rec, net)

        def observe(ys, rows=rows):
            return jnp.log10(jnp.maximum(ys[jnp.asarray(rows)][:, jnp.asarray(obs_idx)].T.ravel(), 1e-300))

        theta = network.theta_vector(net, theta_of(fd["theta"], net))
        res = {}
        for tol_name in ("tight_1e-10", "fit_tol_1e-8", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            _, jac = network.jacobian(net, observe, rtol=rtol, atol=atol)(theta, pot)
            res[tol_name] = np.asarray(jac)
        f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), pot,
                                                        rtol=1e-12, atol=1e-20).ys))
        h = reference["jacobian_central_fd"][pot_id]["step"]
        jax_fd = []
        for c in range(len(theta)):
            e = jnp.zeros(len(theta)).at[c].set(h)
            jax_fd.append((np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h))
        jax_fd = np.array(jax_fd).T
        cols = [keys.index(c) for c in fd["columns"]]
        eng = np.asarray(fd["jacobian"])
        r = {
            "rows": rows, "n_observables": int(eng.shape[0]), "n_coordinates": len(cols),
            "forward_tight_vs_engine_central_fd": column_errors(res["tight_1e-10"][:, cols], eng),
            "forward_tight_vs_jax_central_fd": column_errors(res["tight_1e-10"], jax_fd),
            "jax_central_fd_vs_engine_central_fd": column_errors(jax_fd[:, cols], eng),
            "forward_fit_tol_vs_forward_tight": column_errors(res["fit_tol_1e-8"], res["tight_1e-10"]),
            "forward_loose_vs_forward_tight": column_errors(res["loose_1e-6"], res["tight_1e-10"]),
            "max_abs_jacobian_per_column": dict(zip(keys, np.max(np.abs(res["tight_1e-10"]), axis=0).tolist())),
        }
        for k in ("forward_tight_vs_engine_central_fd", "forward_tight_vs_jax_central_fd",
                  "jax_central_fd_vs_engine_central_fd"):
            r[k]["worst_coordinate"] = keys[r[k]["worst_column"]] if k == "forward_tight_vs_jax_central_fd" \
                else fd["columns"][r[k]["worst_column"]]
        out[pot_id] = r
        print(f"  {pot_id}: forward vs engine FD {r['forward_tight_vs_engine_central_fd']['worst']:.2e} "
              f"({r['forward_tight_vs_engine_central_fd']['worst_coordinate']}), vs own FD "
              f"{r['forward_tight_vs_jax_central_fd']['worst']:.2e}")
    return out


def stage_residuals(reference, net) -> Dict[str, Any]:
    """The panel's predictions (the scored rows' log10 values) at fit tolerance against engine.predict's."""
    import diffrax
    import jax
    from src.kinetic_core_jax import network

    preds = {p["id"]: p for p in reference["panel_predictions"]}
    rows = []
    for solver_name in ("Kvaerno4", "Kvaerno5"):
        solver = getattr(diffrax, solver_name)()
        for tol_name in ("fit_tol_1e-8", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            run = jax.jit(lambda th, pot: network.integrate(net, th, pot, solver=solver, rtol=rtol, atol=atol).ys)
            worst, items = 0.0, []
            for rec in (r for r in reference["pots"] if r["variant"] == "point" and r["id"] in preds):
                ys = np.asarray(run(theta_of(rec["theta"], net), pot_of(rec, net, segments_key="segments")))
                for row in preds[rec["id"]]["rows"]:
                    if row["predicted"] is None or row["species"] not in OBSERVABLE_OF:
                        continue
                    mine = float(net.observe(ys[-1], [OBSERVABLE_OF[row["species"]]])[0])
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
    """Engine and JAX back to back on the same pots (the engine called in this process for the timing only):
    one lane run (``engine._integrate_program``, the engine's rtol 1e-8) and its forward-difference Jacobian in all
    module coordinates ((1 + n) runs, extrapolated from the median run), against one JAX solve and one exact
    Jacobian."""
    import jax
    import jax.numpy as jnp
    import build_jax_trunk_acrylamide_reference as refgen
    from src.kinetic_core import engine
    from src.kinetic_core_jax import network

    pots_by_id = {p["id"]: p for p in refgen.acrylamide_panel_pots()}
    params = engine.core_parameters(engine.ACRYLAMIDE)
    obs_idx = jnp.asarray([net.index[s] for s in OBSERVED])
    n_coord = len(net.coordinates)
    rows = []
    for label, pot_id in TIMING_POTS:
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
        if pot_id in pots_by_id:
            process = pots_by_id[pot_id]["spec"].process
        elif pot_id == refgen.SHORT["id"]:
            process = refgen.short_process()
        else:
            program = engine.ThermalProgram.from_profile(refgen.PROFILE["points"])
            process = engine.ProcessSpec(thermal=program, ph=refgen.PROFILE["ph"],
                                         water_activity=refgen.PROFILE["water_activity"])
        engine._integrate_program(engine.ACRYLAMIDE, params, dict(rec["initial"]), process)
        secs = []
        for _ in range(10):
            t0 = time.perf_counter()
            engine._integrate_program(engine.ACRYLAMIDE, params, dict(rec["initial"]), process)
            secs.append(time.perf_counter() - t0)
        med = float(np.median(secs))
        engine_row = {"lane_run_median_s": med, "forward_fd_jacobian_s": (1 + n_coord) * med,
                      "fd_runs": 1 + n_coord}
        pot = pot_of(rec, net, segments_key="segments")

        def observe(ys):
            return jnp.log10(jnp.maximum(ys[-1][obs_idx], 1e-300))

        theta = network.theta_vector(net, net.theta0())
        out = {"class": label, "pot": pot_id, "n_segments": len(rec["segments"]), "engine": engine_row, "jax": {}}
        for tol_name in ("fit_tol_1e-8", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            solve = jax.jit(lambda v, p: network.integrate(net, network.theta_from_vector(net, v), p,
                                                           rtol=rtol, atol=atol))
            t0 = time.perf_counter()
            jax.block_until_ready(solve(theta, pot))
            compile_s = time.perf_counter() - t0
            js = []
            for _ in range(5):
                t0 = time.perf_counter()
                res = jax.block_until_ready(solve(theta, pot))
                js.append(time.perf_counter() - t0)
            jac_fn = network.jacobian(net, observe, rtol=rtol, atol=atol)
            t0 = time.perf_counter()
            jax.block_until_ready(jac_fn(theta, pot))
            jac_compile = time.perf_counter() - t0
            jj = []
            for _ in range(3):
                t0 = time.perf_counter()
                jax.block_until_ready(jac_fn(theta, pot))
                jj.append(time.perf_counter() - t0)
            out["jax"][tol_name] = {"solve_median_s": float(np.median(js)), "solve_compile_s": compile_s,
                                    "steps": int(np.sum(np.asarray(res.num_steps))),
                                    "jacobian_median_s": float(np.median(jj)), "jacobian_compile_s": jac_compile,
                                    "n_coordinates": n_coord}
        out["uptime"] = uptime()
        rows.append(out)
        j8 = out["jax"]["fit_tol_1e-8"]
        print(f"  timing {pot_id}: engine run {med:.2e} s, FD Jacobian {engine_row['forward_fd_jacobian_s']:.2e} s; "
              f"JAX solve {j8['solve_median_s']:.2e} s, Jacobian {j8['jacobian_median_s']:.2e} s")
    return {"rows": rows}


# --------------------------------------------------------------------------- the write-up
def write_markdown(report: Dict[str, Any]) -> str:
    rates, traj, jac, resid, timing = (report[k] for k in ("rates", "trajectories", "jacobians", "residuals",
                                                           "timing"))
    ts = report["trunk_sets"]
    lines = [
        "# The trunk and the acrylamide lane as A2 contract v2 modules in JAX + diffrax (plan step A2 + H3)",
        "",
        f"*Generated by `scripts/generators/build_jax_trunk_acrylamide_module.py` in the separate `maillard-jax` env "
        f"(jax {report['jax_version']}, diffrax {report['diffrax_version']}, float64, CPU); the engine's numbers from "
        "`jax_trunk_acrylamide_module_reference.json` (core env, `build_jax_trunk_acrylamide_reference.py`). Nothing "
        "shipped moves; the port is not imported by `src/kinetic_core`. Contract: "
        "[A2_MODULE_CONTRACT.md](../../docs/guides/A2_MODULE_CONTRACT.md) v2.*",
        "",
        "## What was ported",
        "",
        f"- `src/kinetic_core_jax/trunk.py`: `network.REACTIONS` ({ts['network.REACTIONS']} steps: B1's fifteen "
        "Martins steps and B7's eleven furanic steps, `r_amadori` derived from `r_schiff`), the trunk the acrylamide "
        f"lane runs. `network.TRUNK_REACTIONS` ({ts['network.TRUNK_REACTIONS']}, the trunk LANE) is it plus seven "
        f"blocks (dicarbonyl, pyrazine, glycation, aqueous glyoxal, methionine, proline, fed 3-deoxy); sulfur's "
        f"{ts['sulfur.FULL_REACTIONS']} steps are it plus the fed 3-deoxy block plus the sulfur steps. Both start "
        "with these 26 steps in this order (checked by the reference generator).",
        f"- `src/kinetic_core_jax/acrylamide.py`: `ACRYLAMIDE_REACTIONS` ({ts['acrylamide.ACRYLAMIDE_REACTIONS']} "
        "steps, `a_sb_int1` derived from `a_asn_glc_sb`) and every shipped condition term of "
        "`acrylamide_conditions` (B14 / B15 a_w and pH factors, CF-039 Ea draws, CF-104 band outside the window), "
        "with CF-106's low-temperature refusal at the pot. Composed after the trunk: "
        f"{report['n_species']} species, {report['n_reactions']} reactions, {len(report['coordinates'])} coordinates.",
        "- `network.py`: pot scalars as arguments (`pot['ph_label']`, `pot['aw']`, read through `RateCtx` by "
        "modules that list them in `Module.pot_inputs`); pH mode \"prescribed\" (a declared pH, constant over the "
        "programme) now runs.",
        "- The trunk LANE's own declared terms (`trunk_conditions`: a_w on `k_schiff`, pH on the three Amadori-decay "
        "steps, the 3-DG exit term) are ported in `trunk.build(lane=\"trunk\")`, accepted at the rate level only "
        "(the trunk lane's trajectories need the seven blocks above).",
        "",
        "## Rates (no integration)",
        "",
        f"- Acrylamide lane, {rates['acrylamide_lane']['n_rate_values']} values ({rates['acrylamide_lane']['n_grid_points']} "
        f"(T, pH, a_w) points x {report['n_reactions']} reactions, shipped and perturbed): worst relative difference "
        f"**{rates['acrylamide_lane']['worst_relative']:.1e}** (at {rates['acrylamide_lane']['worst_at']}).",
        f"- Trunk-lane terms on the 26 trunk steps, {rates['trunk_lane_terms']['n_rate_values']} values (shipped, "
        f"perturbed through the Config, exit term off): worst **{rates['trunk_lane_terms']['worst_relative']:.1e}** "
        f"(at {rates['trunk_lane_terms']['worst_at']}).",
        "",
        "## Trajectories (state after every refined segment vs the engine at LSODA rtol 1e-12; dex where the engine "
        f"holds more than {traj['dex_floor_share_of_max_initial']:g} of the pot's largest initial charge)",
        "",
        "| solver | tolerance | pots | max dex, all species | worst (pot, species) | max dex, ACR and HMF | steps (all pots) | wall s after JIT | max abs below floor / max y0 | element drift |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in traj["rows"]:
        w = r["worst"]
        lines.append(f"| {r['solver']} | {r['tolerance']} | {r['group']} ({r['n_pots']}) | {r['max_dex']:.1e} | "
                     f"{w.get('pot')} {w.get('worst_species')} | {r['max_dex_scored_products']:.1e} | "
                     f"{r['steps_total']} | {r['wall_s_all_pots_after_jit']:.3f} | "
                     f"{r['max_abs_below_floor_share_of_y0']:.1e} | {r['max_rel_element_drift']:.1e} |")
    lines += ["", "Panel pots: the 6 bundles the engine routes to the acrylamide lane and answers, each at its point "
              "and at the low and high ends of every sampled condition band (18 records, one compile). Profile: a "
              "synthetic staircase (25 -> 180 -> 100 C, 2 K steps) with every competitor charged, at pH 5.0 and "
              "a_w 0.80, so all 42 steps and every condition term carry flux.", "",
              "## Exact Jacobian (forward mode) vs central differences, at an interior point", "",
              "Worst column error = max over observables of |difference| / max |reference| in that column, on the "
              f"columns whose largest entry is at least {RELATIVE_COLUMN_SHARE:g} of the Jacobian's largest (\"rel\"); "
              "the worst absolute difference over ALL columns, in dex per coordinate unit, beside it (\"abs\"). "
              "Observables: log10 ACR and HMF at the save rows. Every coordinate is moved in the engine through its "
              "own arguments (a parameter dict, the acrylamide draw); the engine's differences carry LSODA's noise "
              "at rtol 1e-12, which the JAX-FD-vs-engine-FD column measures.", "",
              "| pot | observables x coordinates (relative columns) | vs engine central FD: rel (worst coordinate) / abs | vs JAX central FD: rel / abs | JAX FD vs engine FD: rel / abs | fit tol vs tight | rtol 1e-6 vs tight |",
              "|---|---|---|---|---|---|---|"]
    for pot_id, r in jac.items():
        e = r["forward_tight_vs_engine_central_fd"]
        o = r["forward_tight_vs_jax_central_fd"]
        x = r["jax_central_fd_vs_engine_central_fd"]
        lines.append(f"| {pot_id} | {r['n_observables']} x {r['n_coordinates']} ({e['n_relative_columns']}) | "
                     f"{e['worst']:.1e} ({e['worst_coordinate']}) / {e['worst_absolute']:.1e} | "
                     f"{o['worst']:.1e} / {o['worst_absolute']:.1e} | {x['worst']:.1e} / {x['worst_absolute']:.1e} | "
                     f"{r['forward_fit_tol_vs_forward_tight']['worst']:.1e} | "
                     f"{r['forward_loose_vs_forward_tight']['worst']:.1e} |")
    lines += ["", "## Panel predictions at fit tolerance vs `engine.predict` (log10 of each scored row's value)", "",
              "| solver | tolerance | rows | max abs difference, dex |", "|---|---|---|---|"]
    for r in resid["rows"]:
        lines.append(f"| {r['solver']} | {r['tolerance']} | {r['n_rows']} | {r['max_abs_diff_dex']:.1e} |")
    lines += ["", "## Timing (CPU, one core, niced; shared machine: load beside each row; provisional)", "",
              f"The engine's forward-difference Jacobian is (1 + {len(report['coordinates'])}) lane runs, from the "
              "median run; the JAX Jacobian is one forward-mode solve over all coordinates (Kvaerno5).", "",
              "| pot class | segments | engine run s | engine FD Jacobian s | JAX solve s (rtol 1e-8) | steps | compile s | JAX Jacobian s (rtol 1e-8) | compile s | JAX solve s (1e-6) | JAX Jacobian s (1e-6) | load |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in timing["rows"]:
        e, j8, j6 = r["engine"], r["jax"]["fit_tol_1e-8"], r["jax"]["loose_1e-6"]
        load = r["uptime"].split("load averages:")[-1].strip() if "load" in r["uptime"] else r["uptime"]
        lines.append(f"| {r['class']} | {r['n_segments']} | {e['lane_run_median_s']:.1e} | "
                     f"{e['forward_fd_jacobian_s']:.1e} | {j8['solve_median_s']:.1e} | {j8['steps']} | "
                     f"{j8['solve_compile_s']:.1f} | {j8['jacobian_median_s']:.1e} | {j8['jacobian_compile_s']:.1f} | "
                     f"{j6['solve_median_s']:.1e} | {j6['jacobian_median_s']:.1e} | {load} |")
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
    report.update({"artifact": "jax_trunk_acrylamide_module",
                   "generated_by": "scripts/generators/build_jax_trunk_acrylamide_module.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "diffrax_version": diffrax.__version__,
                   "n_species": net.n_species, "n_reactions": len(net.reactions), "species": list(net.species),
                   "coordinates": [f"{m}.{k}" for m, k in net.coordinates], "reactions": list(net.reactions),
                   "ph_mode": net.ph_mode, "trunk_sets": reference["trunk_sets"],
                   "config_digests": {m.name: m.config_digest for m in net.modules}})
    for stage in stages:
        print(f"[{stage}] {uptime()}")
        report[stage] = globals()[f"stage_{stage}"](reference, net)
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
    "1. **Match: yes, at the contract's bar with margin.**",
    "   - Rates equal the engine's functions to ~1e-15 relative: 81 648 values on the composed acrylamide lane (every "
    "condition term on and off, a perturbed parameter set and draw), 75 816 on the trunk lane's declared terms.",
    "   - Every acrylamide panel pot, at its point and both band corners, matches to <= 2.2e-9 dex at rtol 1e-10 "
    "(Kvaerno5; Kvaerno4 1.4e-9); acrylamide and HMF to <= 1.9e-11. The staircase profile with all 42 steps: "
    "3.3e-9. The worst species is a trace pool (glucose formed from fructose in the Lin pot, cysteine in the profile).",
    "   - Panel predictions at the fits' tolerances equal `engine.predict`'s to <= 3e-9 dex at rtol 1e-8 and "
    "<= 6e-8 dex at rtol 1e-6, against the contract's 1e-4; so rtol 1e-6 is licensed for fits on this lane.",
    "   - C / N / S close to ~5e-15 relative on every pot.",
    "2. **Exact Jacobian: correct.** On 3 pots x 84 coordinates (every log10 k_ref and Ea, and the four condition "
    "coordinates), forward mode equals the port's own central differences to <= 1.3e-6 relative on the columns that "
    "matter and <= 1e-9 dex per unit everywhere. Against the engine's central differences it is <= 8e-5 relative / "
    "3.3e-7 absolute, and the JAX-FD-vs-engine-FD column shows the same numbers: that is the engine's LSODA noise "
    "(rtol 1e-12, step 1e-5), not the port. Unlike the lipid lane, the engine here is an ODE solve, so its differences "
    "are not noise-free; the contract's \"the port's own otherwise\" applies.",
    "3. **Kvaerno4 vs Kvaerno5: Kvaerno4 here.** On this network Kvaerno4 takes ~0.6x the steps at every tolerance "
    "(20 030 vs 33 368 on the panel at rtol 1e-8) at equal or better accuracy. The timing below is Kvaerno5 (the "
    "contract default); the composed union decides on its own pots.",
    "4. **atol drives the step count.** On the 30-min Chang pot at rtol 1e-8, Kvaerno5 takes 2248 steps at atol "
    "1e-14, 947 at 1e-10 and 410 at 1e-8, with acrylamide within 1.4e-10 / 5.9e-10 dex of the engine at the last "
    "two (Kvaerno4: 1328 / 709 / 462). Several pools stay at exactly zero or are depleted far below 1e-14. The "
    "smallest scored species here (acrylamide, ~1e-4 mmol/L at 10 ppb) licenses atol ~1e-10 for fits: a further "
    "2 - 3x on the solve and the Jacobian.",
    "5. **Cost (provisional, Kvaerno5, rtol 1e-8, atol 1e-14, one core, load ~4 - 5).**",
    "   - One solve: 0.6 - 0.8x the engine's time on the short isothermal pot and the staircase profile, ~3x slower on the "
    "isothermal 30-min holds (LSODA is very good at one long stiff segment).",
    "   - The Jacobian in all 84 coordinates: 2.8 - 2.9x faster than the engine's forward differences on the holds, "
    "11x on the profile, 15x on the short pot. With atol 1e-10 and Kvaerno4 (findings 3 - 4) the holds should gain "
    "another ~3x.",
    "   - Compile: 1 - 2.4 s per solve, 2.5 - 5.4 s per Jacobian, once per configuration and segment count.",
    "6. **What is and is not ported.**",
    "   - Ported for every value: the acrylamide lane's terms (formation and elimination a_w, both pH factors with "
    "the window hold, the CF-039 Ea draws, the CF-104 band outside the window), `panel.thermal_history` (the "
    "programme is pot data), and on the trunk lane `trunk.three_deoxy_exit_ph_term` (on / off) and "
    "`trunk.three_deoxy_exit_ph` for steps in this module.",
    "   - Refused (`UnsupportedConfiguration`): a programme peaking in (60, 120) C (CF-106, as the engine); a charged "
    "species the network does not carry; an exit term on a step outside `network.REACTIONS`; a non-air atmosphere "
    "on the trunk lane; a non-unit statistical factor on a trunk step (the two engine paths disagree there).",
    "   - Not ported: the trunk LANE's seven extra blocks (55 - 26 steps), its pyrazine / methionine / proline pH "
    "terms and the atmosphere factor (they act on those blocks), and the engine's convention of folding MET / PRO "
    "into Gly on the trunk lane. So the trunk lane's own trajectories are not yet accepted; its terms on these 26 "
    "steps are, at the rate level.",
    "   - Composition: the acrylamide module owns Cys, CBX, FRAG_N and FRAG_S, so `compose` refuses it with a "
    "sulfur module (as `engine.resolve_lanes` refuses the two lanes); declared Cys competition is A2b's.",
    "7. **Effort.** The two modules, the contract's pot scalars, the two generators and the tests took ~1.5 h of agent "
    "work against the 13 - 23 h estimated: no pH state, a pure Arrhenius rate law with two derived steps, and a "
    "reference path that is the engine's own per-segment functions. Remaining A2 + H3, revised:",
    "   - the trunk lane's seven blocks and their terms (dicarbonyl, pyrazine, glycation, glyoxal, methionine, "
    "proline, fed 3-deoxy; the Gly folding; the atmosphere factor): 3 - 6 h;",
    "   - sulfur (every remaining switch, the pH solver in `integrate`, profile mode): 12 - 20 h;",
    "   - network infrastructure (merged titratable table and charge closure, declared competition, forced inputs "
    "including the trajectory form of `prescribed`, the fit driver): 15 - 25 h;",
    "   - the pH-root performance step: 8 - 15 h.",
    "   - Total remaining ~40 - 65 h (p ~0.6 of landing inside). The pH-coupled sulfur lane and the merge carry most "
    "of the uncertainty; the pH-free ports have each come in far under estimate.",
    "",
]


if __name__ == "__main__":
    sys.exit(main())

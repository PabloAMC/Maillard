#!/usr/bin/env python
"""
build_jax_lipid_module.py -- the lipid lane as the first module of the A2 contract v2, in JAX + diffrax (plan step
A2 + H3; docs/guides/A2_MODULE_CONTRACT.md; 2026-11-05): does the module reproduce the engine, are its exact
Jacobians right, what do one solve, its Jacobian and a batch of draws cost?

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml), never in the core env:

    source <conda.sh>; conda activate maillard-jax
    XLA_FLAGS=--xla_cpu_multi_thread_eigen=false OMP_NUM_THREADS=1 nice -n 10 \
        python scripts/generators/build_jax_lipid_module.py

Every compared engine number comes from results/validation/jax_lipid_module_reference.json (core env,
``build_jax_lipid_module_reference.py``); the engine is called in this process only to BUILD the module (its static
data) and, back to back with the JAX runs, for the paired timings. Writes results/validation/jax_lipid_module.{json,md}.
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

REFERENCE = ROOT / "results" / "validation" / "jax_lipid_module_reference.json"
OUT = ROOT / "results" / "validation" / "jax_lipid_module"
TOLERANCES = {"tight_1e-10": (1e-10, 1e-16), "fit_tol_1e-8": (1e-8, 1e-14), "loose_1e-6": (1e-6, 1e-12)}
#: a species is compared in dex where the engine holds more than this share of the pot's initial hydroperoxide
DEX_FLOOR = 1e-9
#: the lane's scored products (what the panel's lipid rows read)
OBSERVED = ("HEXANAL", "PENTYLFURAN", "NONANAL", "DECADIENAL")


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def setup():
    from src.kinetic_core_jax import lipid, network

    reference = json.loads(REFERENCE.read_text())
    net = network.compose([lipid.build()])
    if list(net.species) != reference["species"]:
        raise AssertionError("species order differs from the reference")
    return reference, net


def pot_of(rec: Dict[str, Any], net, n_segments=None) -> Dict[str, Any]:
    from src.kinetic_core_jax import lipid, network

    return {"y0": lipid.initial_state(rec["carriers"], net.species),
            "segments": network.program_segments(rec["refined_segments"], n_segments)}


def theta_of(values: Dict[str, float]):
    import jax.numpy as jnp

    return {"lipid": {k: jnp.asarray(float(v)) for k, v in values.items()}}


def compare(ys: np.ndarray, rec: Dict[str, Any], y0: np.ndarray, species) -> Dict[str, Any]:
    """Worst dex over species above the floor, and the worst absolute difference (share of y0) below it."""
    eng = np.asarray(rec["states"])
    ys = ys[: len(eng)]
    scale = float(np.max(y0))
    above = eng > DEX_FLOOR * scale
    dex = np.where(above, np.abs(np.log10(np.where(above, np.maximum(ys, 1e-300), 1.0))
                                 - np.log10(np.where(above, eng, 1.0))), 0.0)
    r, i = np.unravel_index(int(np.argmax(dex)), dex.shape)
    below = np.where(~above, np.abs(ys - eng), 0.0) / scale
    scored = [list(species).index(s) for s in OBSERVED if s in species]
    return {"max_dex": float(dex.max()), "worst_species": species[i], "worst_row": int(r),
            "max_dex_scored_products": float(dex[:, scored].max()),
            "max_abs_below_floor_share_of_y0": float(below.max())}


# --------------------------------------------------------------------------- stages
def stage_rates(reference, net) -> Dict[str, Any]:
    """The module's rate vector against the engine's functions on the reference probe's (T, Q10, donor) grid."""
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import network

    probe = reference["rate_probe"]
    mod = net.modules[0]
    f = jax.jit(lambda th, t_c: mod.rate_constants(th, network.RateCtx(t_c=t_c, t_k=t_c + network.CELSIUS)))
    worst, worst_at, n = 0.0, None, 0
    for label, brk_key, k_brk_key in (("shipped", "breakage_shipped", "k_breakage_per_min"),
                                      ("perturbed", "breakage_perturbed", "k_breakage_perturbed_per_min")):
        brk = probe[brk_key]
        ship = probe["breakage_shipped"]
        for q in probe["q10"]:
            for d in probe["donors"]:
                th = {"q10": q, "log10_k_looh": 0.0, "donor": d,
                      "log10_k_decadienal": math.log10(brk["k_loss_per_min_at_ref"] / ship["k_loss_per_min_at_ref"]),
                      "ea_decadienal_kj_mol": brk["ea_loss_kj_mol"],
                      "share_hexanal_from_decadienal": brk["share_hexanal"]}
                th = {k: jnp.asarray(float(v)) for k, v in th.items()}
                for ti, t_c in enumerate(probe["temperatures_c"]):
                    k = np.asarray(f(th, jnp.asarray(t_c)))
                    expect = expected_rates(net.reactions, probe, repr(q), repr(d), ti, k_brk_key, brk)
                    rel = np.abs(k - expect) / np.maximum(np.abs(expect), 1e-300)
                    rel = np.where((k == 0.0) & (expect == 0.0), 0.0, rel)
                    n += len(k)
                    j = int(np.argmax(rel))
                    if rel[j] > worst:
                        worst, worst_at = float(rel[j]), (label, q, d, t_c, net.reactions[j])
    return {"n_rate_values": n, "worst_relative": worst, "worst_at": worst_at}


def expected_rates(reactions, probe, q_key, d_key, ti, k_brk_key, brk) -> np.ndarray:
    """The engine's rate of each module reaction, assembled from the probe's engine numbers."""
    k_looh = probe["k_looh_per_min"][q_key][ti]
    k_brk = probe[k_brk_key][ti] if probe["temperatures_c"][ti] >= probe["breakage_window_min_c"] else 0.0
    yields = probe["yields"][d_key]
    out = []
    for key in reactions:
        if key.startswith("lipid_decomp_"):
            out.append(k_looh)
        elif key.startswith("lipid_slate_"):
            rest = key[len("lipid_slate_"):]
            pool = next(p for p in yields["cells"] if rest.startswith(p + "_"))
            out.append(k_looh * yields["cells"][pool]["slate"][rest[len(pool) + 1:]])
        elif key.startswith("lipid_pentylfuran_"):
            out.append(k_looh * yields["cells"][key[len("lipid_pentylfuran_"):]]["pentylfuran"])
        elif key == "lipid_nonanal_LOOH_OL":
            out.append(k_looh * yields["nonanal_per_oleate_looh"])
        elif key == "lipid_decadienal_loss":
            out.append(k_brk)
        elif key == "lipid_decadienal_hexanal":
            out.append(k_brk * brk["share_hexanal"])
        else:
            raise KeyError(key)
    return np.asarray(out)


def stage_trajectories(reference, net) -> Dict[str, Any]:
    import diffrax
    import jax
    from src.kinetic_core_jax import network

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
                steps, below, secs, scored = 0, 0.0, 0.0, 0.0
                for rec in recs:
                    pot = pot_of(rec, net, s)
                    out = run(theta_of(rec["theta"]), pot)
                    out.ys.block_until_ready()
                    t0 = time.perf_counter()
                    out = run(theta_of(rec["theta"]), pot)
                    out.ys.block_until_ready()
                    secs += time.perf_counter() - t0
                    c = compare(np.asarray(out.ys), rec, pot["y0"], net.species)
                    c["pot"] = f"{rec['id']}:{rec['variant']}"
                    steps += int(np.sum(out.num_steps))
                    below = max(below, c["max_abs_below_floor_share_of_y0"])
                    scored = max(scored, c["max_dex_scored_products"])
                    if c["max_dex"] >= worst["max_dex"]:
                        worst = c
                rows.append({"solver": solver_name, "tolerance": tol_name, "group": group, "n_pots": len(recs),
                             "max_dex": worst["max_dex"], "worst": worst, "steps_total": steps,
                             "max_dex_scored_products": scored,
                             "max_abs_below_floor_share_of_y0": below, "wall_s_all_pots_after_jit": secs})
                print(f"  {solver_name} {tol_name} {group}: max dex {worst['max_dex']:.2e} ({worst.get('pot')}, "
                      f"{worst.get('worst_species')}), {steps} steps, {secs:.3f} s")
    return {"dex_floor_share_of_initial_looh": DEX_FLOOR, "rows": rows}


def stage_jacobians(reference, net) -> Dict[str, Any]:
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import network

    obs_idx = [net.index[s] for s in reference["observable_species"]]
    out = {}
    for pot_id, fd in reference["jacobian_central_fd"].items():
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
        rows = fd["rows"]
        pot = pot_of(rec, net)

        def observe(ys, rows=rows):
            return jnp.log10(jnp.maximum(ys[jnp.asarray(rows)][:, jnp.asarray(obs_idx)].T.ravel(), 1e-300))

        theta = jnp.asarray([fd["theta"][k] for _, k in net.coordinates])
        res = {}
        for tol_name in ("tight_1e-10", "fit_tol_1e-8", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            jac_fn = network.jacobian(net, observe, rtol=rtol, atol=atol)
            _, jac = jac_fn(theta, pot)
            res[tol_name] = np.asarray(jac)
        f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), pot,
                                                        rtol=1e-12, atol=1e-18).ys))
        jax_fd = []
        for c in range(len(theta)):
            h = 1e-5 * max(1.0, abs(float(theta[c])))
            e = jnp.zeros(len(theta)).at[c].set(h)
            jax_fd.append((np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h))
        jax_fd = np.array(jax_fd).T
        eng = np.asarray(fd["jacobian"])
        cols = engine_columns(net, fd)
        out[pot_id] = {
            "rows": rows, "theta": fd["theta"], "n_observables": int(eng.shape[0]),
            "engine_columns": fd["columns"],
            "forward_tight_vs_engine_central_fd": column_errors(res["tight_1e-10"][:, cols], eng),
            "forward_tight_vs_jax_central_fd": column_errors(res["tight_1e-10"], jax_fd),
            "jax_central_fd_vs_engine_central_fd": column_errors(jax_fd[:, cols], eng),
            "forward_fit_tol_vs_forward_tight": column_errors(res["fit_tol_1e-8"], res["tight_1e-10"]),
            "forward_loose_vs_forward_tight": column_errors(res["loose_1e-6"], res["tight_1e-10"]),
            "max_abs_jacobian_per_column": dict(zip([k for _, k in net.coordinates],
                                                    np.max(np.abs(res["tight_1e-10"]), axis=0).tolist())),
        }
        print(f"  {pot_id}: forward vs engine central FD {out[pot_id]['forward_tight_vs_engine_central_fd']['worst']:.2e}")
    return out


def engine_columns(net, fd: Dict[str, Any]) -> List[int]:
    """The network-coordinate indices of the reference Jacobian's columns (the coordinates the engine can move)."""
    keys = [k for _, k in net.coordinates]
    return [keys.index(c) for c in fd["columns"]]


def column_errors(a: np.ndarray, b: np.ndarray, flat: float = 1e-6) -> Dict[str, Any]:
    """Per column: max |a - b| / max |b| (absolute where the column's max |b| < ``flat``)."""
    cols = []
    for c in range(b.shape[1]):
        scale = float(np.max(np.abs(b[:, c])))
        diff = float(np.max(np.abs(a[:, c] - b[:, c])))
        cols.append(diff / scale if scale >= flat else diff)
    flat_cols = [c for c in range(b.shape[1]) if float(np.max(np.abs(b[:, c]))) < flat]
    return {"per_column": cols, "worst": float(max(cols)), "flat_columns_absolute": flat_cols}


def stage_timing(reference, net) -> Dict[str, Any]:
    """Engine and JAX back to back on the same pots (the engine called in this process for the timing only)."""
    import jax
    import jax.numpy as jnp
    import build_jax_lipid_module_reference as refgen
    from src.kinetic_core_jax import network

    obs_idx = [net.index[s] for s in reference["observable_species"]]
    pots_by_id = {p["id"]: p for p in refgen.lipid_panel_pots()}
    rows = []
    for pot_id in ("pea_isolate_uht_140C_Trikusuma2019", "external_validation_bi_2020_roasted_pea_hexanal",
                   "synthetic_profile_pea_soy_25_160_60"):
        rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
        engine_row = {}
        if pot_id in pots_by_id:
            engine_row = refgen.engine_timing(pots_by_id[pot_id], repeats=20)
        else:   # the synthetic profile: the lane function _run_lipid_lane calls, once per carrier
            from src.kinetic_core import engine as eng_mod
            from src.kinetic_core import lipid as lane
            from src.kinetic_core import parameters_lipid as pl

            segs = [tuple(s) for s in rec["segments"]]
            secs = []
            for _ in range(10):
                t0 = time.perf_counter()
                branch, comp = eng_mod.core_lipid_model()
                for c in rec["carriers"]:
                    lane.integrate_lipid(lane.charge_from_carrier(pl.LIPID_CARRIERS[c], comp), segs, branch)
                secs.append(time.perf_counter() - t0)
            med = float(np.median(secs))
            engine_row = {"pot": pot_id, "lane_run_median_s": med, "forward_fd_jacobian_s": 7 * med,
                          "fd_runs": 7, "note": "integrate_lipid per carrier; FD Jacobian = 7 runs"}
        pot = pot_of({**rec, "refined_segments": rec["segments"]}, net)
        rows_obs = [len(rec["segments"]) - 1]

        def observe(ys, r=rows_obs):
            return jnp.log10(jnp.maximum(ys[jnp.asarray(r)][:, jnp.asarray(obs_idx)].ravel(), 1e-300))

        theta = network.theta_vector(net, net.theta0())
        out = {"pot": pot_id, "n_segments": len(rec["segments"]), "engine": engine_row, "jax": {}}
        for tol_name in ("fit_tol_1e-8", "loose_1e-6"):
            rtol, atol = TOLERANCES[tol_name]
            solve = jax.jit(lambda v, p: network.integrate(net, network.theta_from_vector(net, v), p,
                                                           rtol=rtol, atol=atol).ys)
            t0 = time.perf_counter()
            solve(theta, pot).block_until_ready()
            compile_s = time.perf_counter() - t0
            secs = []
            for _ in range(10):
                t0 = time.perf_counter()
                solve(theta, pot).block_until_ready()
                secs.append(time.perf_counter() - t0)
            jac_fn = network.jacobian(net, observe, rtol=rtol, atol=atol)
            t0 = time.perf_counter()
            jax.block_until_ready(jac_fn(theta, pot))
            jac_compile = time.perf_counter() - t0
            jsecs = []
            for _ in range(10):
                t0 = time.perf_counter()
                jax.block_until_ready(jac_fn(theta, pot))
                jsecs.append(time.perf_counter() - t0)
            out["jax"][tol_name] = {"solve_median_s": float(np.median(secs)), "solve_compile_s": compile_s,
                                    "jacobian_median_s": float(np.median(jsecs)), "jacobian_compile_s": jac_compile}
        out["uptime"] = uptime()
        rows.append(out)
        print(f"  timing {pot_id}: engine {engine_row.get('lane_run_median_s', float('nan')):.2e} s, "
              f"JAX {out['jax']['fit_tol_1e-8']['solve_median_s']:.2e} s")
    return {"rows": rows}


def stage_batch(reference, net) -> Dict[str, Any]:
    """Draws in lockstep (``network.batched`` / ``network.grid``): cost per member vs batch size, and the batch equal
    to the single solves. CPU only; the same functions are the harness for the later CPU-vs-GPU benchmark."""
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import network

    rng = np.random.default_rng(20261105)
    mod = net.modules[0]
    rec = next(r for r in reference["pots"] if r["id"] == "external_validation_bi_2020_roasted_pea_hexanal"
               and r["variant"] == "point")
    pot = pot_of(rec, net)
    rtol, atol = TOLERANCES["fit_tol_1e-8"]

    def draws(n):
        lo, hi = mod.bands["q10"]
        base = {k: float(v) for k, v in mod.params.items()}
        out = []
        for _ in range(n):
            d = dict(base)
            d["q10"] = float(rng.uniform(lo, hi))
            d["log10_k_looh"] = float(rng.normal(0.0, 0.3))
            d["donor"] = float(rng.uniform(0.0, 0.5))
            d["log10_k_decadienal"] = float(rng.normal(0.0, 0.2))
            d["ea_decadienal_kj_mol"] = float(rng.uniform(20.0, 35.0))
            out.append({"lipid": d})
        return out

    single = jax.jit(lambda th, p: network.integrate(net, th, p, rtol=rtol, atol=atol).ys)
    sizes = []
    for n in (1, 16, 128, 1024):
        thetas = draws(n)
        batch = network.stack_theta(thetas)
        f = network.batched(net, theta_axis=0, pot_axis=None, rtol=rtol, atol=atol)
        t0 = time.perf_counter()
        ys = f(batch, pot)
        ys.block_until_ready()
        compile_s = time.perf_counter() - t0
        secs = []
        for _ in range(3):
            t0 = time.perf_counter()
            f(batch, pot).block_until_ready()
            secs.append(time.perf_counter() - t0)
        check = [0, n - 1]
        dex = 0.0
        for i in check:
            one = np.asarray(single(jax.tree_util.tree_map(jnp.asarray, thetas[i]), pot))
            b = np.asarray(ys[i])
            m = one > DEX_FLOOR * float(np.max(pot["y0"]))
            dex = max(dex, float(np.max(np.abs(np.log10(np.where(m, b, 1)) - np.log10(np.where(m, one, 1))))))
        med = float(np.median(secs))
        sizes.append({"draws": n, "wall_s": med, "per_draw_s": med / n, "compile_s": compile_s,
                      "batch_vs_single_max_dex": dex})
        print(f"  batch {n}: {med:.3f} s, {med / n * 1e3:.3f} ms per draw (batch vs single {dex:.1e} dex)")
    # every draw on every point pot: one compile for (pots x draws)
    recs = [r for r in reference["pots"] if r["variant"] == "point" and r["panel"] != "synthetic"]
    s = max(len(r["refined_segments"]) for r in recs)
    pots = [pot_of(r, net, s) for r in recs]
    pot_batch = {"y0": jnp.stack([jnp.asarray(p["y0"]) for p in pots]),
                 "segments": jnp.stack([jnp.asarray(p["segments"]) for p in pots])}
    thetas = network.stack_theta(draws(64))
    g = network.grid(net, rtol=rtol, atol=atol)
    t0 = time.perf_counter()
    g(thetas, pot_batch).block_until_ready()
    compile_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    out = g(thetas, pot_batch)
    out.block_until_ready()
    wall = time.perf_counter() - t0
    grid_row = {"pots": len(pots), "draws": 64, "shape": list(out.shape), "wall_s": wall, "compile_s": compile_s,
                "per_member_s": wall / (len(pots) * 64)}
    print(f"  grid {len(pots)} pots x 64 draws: {wall:.3f} s")
    return {"pot": rec["id"], "tolerance": "fit_tol_1e-8", "sizes": sizes, "grid": grid_row, "uptime": uptime()}


# --------------------------------------------------------------------------- the write-up
def write_markdown(report: Dict[str, Any]) -> str:
    rates, traj, jac, timing, batch = (report[k] for k in ("rates", "trajectories", "jacobians", "timing", "batch"))
    lines = [
        "# The lipid lane as an A2 contract v2 module in JAX + diffrax (plan step A2 + H3)",
        "",
        f"*Generated by `scripts/generators/build_jax_lipid_module.py` in the separate `maillard-jax` env "
        f"(jax {report['jax_version']}, diffrax {report['diffrax_version']}, float64, CPU); the engine's numbers from "
        "`jax_lipid_module_reference.json` (core env, `build_jax_lipid_module_reference.py`). Nothing shipped moves; "
        "the port is not imported by `src/kinetic_core`. Contract: "
        "[A2_MODULE_CONTRACT.md](../../docs/guides/A2_MODULE_CONTRACT.md) v2.*",
        "",
        "## What was ported",
        "",
        f"- `src/kinetic_core_jax/lipid.py`: the lane's {report['n_species']} species and {report['n_reactions']} "
        "mass-action reactions (decomposition of five hydroperoxide pools; one yield channel per (cell, product) at "
        "Frankel's frozen simplexes; 2-pentylfuran and nonanal on the 1981 ratios; decadienal breakage and its hexanal "
        "share above 100 C), built from the engine's own definitions; `src/kinetic_core_jax/network.py`: the "
        "lane-agnostic contract (Module, compose, integrate, batched / grid, jacobian).",
        f"- Free coordinates: {', '.join(report['coordinates'])}.",
        "- The engine solves the lane analytically per segment; the module is the network whose exact solution that "
        "arithmetic is, integrated by the same diffrax solve as every other module will be.",
        "",
        "## Rates (no integration)",
        "",
        f"{rates['n_rate_values']} rate values ((T, Q10, donor) grid x {report['n_reactions']} reactions, shipped and "
        f"perturbed breakage law): worst relative difference **{rates['worst_relative']:.1e}** "
        f"(at {rates['worst_at']}).",
        "",
        "## Trajectories (state after every refined segment vs the engine; dex where the engine holds more than "
        f"{traj['dex_floor_share_of_initial_looh']:g} of the pot's initial hydroperoxide)",
        "",
        "| solver | tolerance | pots | max dex, all species | worst (pot, species) | max dex, scored products | steps (all pots) | wall s after JIT | max abs below floor / y0 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in traj["rows"]:
        w = r["worst"]
        lines.append(f"| {r['solver']} | {r['tolerance']} | {r['group']} ({r['n_pots']}) | {r['max_dex']:.1e} | "
                     f"{w.get('pot')} {w.get('worst_species')} | {r.get('max_dex_scored_products', float('nan')):.1e} | "
                     f"{r['steps_total']} | "
                     f"{r['wall_s_all_pots_after_jit']:.3f} | {r['max_abs_below_floor_share_of_y0']:.1e} |")
    lines += ["", "Panel pots: the 8 bundles the engine routes to the lipid lane, each at its point and at both "
              "declared-band corners (24 records, one compile). Profile: a synthetic staircase (pea + soy isolate, "
              "25 -> 160 -> 60 C, 2 K steps) that crosses the decadienal window twice.", "",
              "## Exact Jacobian (forward mode) vs central differences, at an interior point",
              "", "Worst column error = max over observables of |difference| / max |reference| in that column "
              "(absolute on columns whose reference is below 1e-6).", "",
              "| pot | observables | vs engine central FD | vs JAX central FD | JAX FD vs engine FD | fit tol vs tight | rtol 1e-6 vs tight |",
              "|---|---|---|---|---|---|---|"]
    for pot_id, r in jac.items():
        lines.append(f"| {pot_id} | {r['n_observables']} | {r['forward_tight_vs_engine_central_fd']['worst']:.1e} | "
                     f"{r['forward_tight_vs_jax_central_fd']['worst']:.1e} | "
                     f"{r['jax_central_fd_vs_engine_central_fd']['worst']:.1e} | "
                     f"{r['forward_fit_tol_vs_forward_tight']['worst']:.1e} | "
                     f"{r['forward_loose_vs_forward_tight']['worst']:.1e} |")
    lines += ["", "## Timing (CPU, one core, niced; shared machine: load beside each row)", "",
              "| pot | segments | engine lane run s | engine forward-FD Jacobian s (7 lane runs) | JAX solve s (rtol 1e-8) | its compile s | JAX Jacobian s (rtol 1e-8) | its compile s | JAX Jacobian s (rtol 1e-6) | load |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for r in timing["rows"]:
        e, j8, j6 = r["engine"], r["jax"]["fit_tol_1e-8"], r["jax"]["loose_1e-6"]
        load = r["uptime"].split("load averages:")[-1].strip() if "load" in r["uptime"] else r["uptime"]
        lines.append(f"| {r['pot']} | {r['n_segments']} | {e['lane_run_median_s']:.1e} | "
                     f"{e['forward_fd_jacobian_s']:.1e} | {j8['solve_median_s']:.1e} | {j8['solve_compile_s']:.1f} | "
                     f"{j8['jacobian_median_s']:.1e} | {j8['jacobian_compile_s']:.1f} | "
                     f"{j6['jacobian_median_s']:.1e} | {load} |")
    lines += ["", "## Batches of draws in lockstep (`network.batched`, one pot, rtol 1e-8)", "",
              "| draws | wall s | per draw ms | compile s | batch vs single, max dex |", "|---|---|---|---|---|"]
    for r in batch["sizes"]:
        lines.append(f"| {r['draws']} | {r['wall_s']:.3f} | {r['per_draw_s'] * 1e3:.3f} | {r['compile_s']:.1f} | "
                     f"{r['batch_vs_single_max_dex']:.1e} |")
    g = batch["grid"]
    lines += ["", f"`network.grid`: {g['pots']} pots x {g['draws']} draws in one call, {g['wall_s']:.3f} s "
              f"({g['per_member_s'] * 1e3:.3f} ms per member), compile {g['compile_s']:.1f} s.", ""]
    lines += report.get("findings", [])
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", action="append", choices=("rates", "trajectories", "jacobians", "timing", "batch"))
    args = parser.parse_args(argv)
    import diffrax
    import jax

    reference, net = setup()
    stages = args.stage or ["rates", "trajectories", "jacobians", "timing", "batch"]
    path = OUT.with_suffix(".json")
    report = json.loads(path.read_text()) if path.exists() and args.stage else {}
    report.update({"artifact": "jax_lipid_module", "generated_by": "scripts/generators/build_jax_lipid_module.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "diffrax_version": diffrax.__version__,
                   "n_species": net.n_species, "n_reactions": len(net.reactions),
                   "coordinates": [k for _, k in net.coordinates], "reactions": list(net.reactions),
                   "config_digest": net.modules[0].config_digest})
    for stage in stages:
        print(f"[{stage}] {uptime()}")
        report[stage] = globals()[f"stage_{stage}"](reference, net)
    report["findings"] = FINDINGS
    path.write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    if all(k in report for k in ("rates", "trajectories", "jacobians", "timing", "batch")):
        OUT.with_suffix(".md").write_text(write_markdown(report), encoding="utf-8")
    print(f"wrote {path}")
    return 0


#: the write-up's closing section, edited by hand after reading the numbers (kept here so the .md regenerates)
FINDINGS: List[str] = [
    "## Findings",
    "",
    "1. **Match: yes, at the contract's bar with margin.**",
    "   - Rates equal the engine's functions to ~1e-15 relative.",
    "   - Every lipid panel pot, at its point and both corners, matches to <= 1.1e-9 dex at rtol 1e-10 (Kvaerno5); "
    "the scored products (hexanal, 2-pentylfuran, nonanal, decadienal) match to <= 4.5e-10. At the fits' rtol 1e-8 "
    "the scored products match to <= 4.3e-8, at rtol 1e-6 to <= 3.6e-6.",
    "   - The worst species is always a hydroperoxide pool depleted to ~1e-9 of its start. The engine's own "
    "split-vs-unsplit noise there is ~1.7e-11 relative: it computes the pool as state - consumed.",
    "2. **Exact Jacobian: correct.**",
    "   - Forward mode equals the engine's central differences to <= 5e-8 per column on four pots. The lane is "
    "analytic in the engine, so those differences carry no solver noise and this is a clean check.",
    "   - The check covers the five coordinates the engine can move. `log10_k_looh` multiplies a constant that is "
    "not a ConfigKey, so it is checked against the port's own central differences (<= 5.3e-8).",
    "   - Columns the physics zeroes (the decadienal law below 100 C) come out exactly 0.",
    "   - The rtol 1e-6 Jacobian is within 1.7e-5 of the tight one.",
    "3. **Kvaerno4 vs Kvaerno5: neither dominates.** On these pots Kvaerno5 takes ~4x fewer steps (973 vs 4216 on "
    "the panel at rtol 1e-8; 158 vs 794 on the 160 C / 30 min pot). Kvaerno4 stalls on pools driven far below atol. "
    "The bake-off found the opposite on the sulfur pot, so the contract makes the solver a declared per-network "
    "choice (`compose(..., solver=)`, default Kvaerno5).",
    "4. **Cost: the engine is faster for this lane alone, as a closed form must be.**",
    "   - One solve: 7x (UHT 140 C, 0.1 min) to ~170x (160 C, 30 min) slower.",
    "   - The Jacobian: ~3x to ~50x slower than seven engine runs.",
    "   - The staircase profile (119 segments): 10x per solve, 4x per Jacobian.",
    "   - Compile: ~2 s per solve, ~3 s per Jacobian, once per segment count.",
    "   - Hours are not at stake: the lane costs milliseconds either way. The port buys composition (the lane can now "
    "share species with sulfur, which M5 needs), exact gradients, and batches. A pH-free linear network has a cheaper "
    "exact path, a matrix exponential per segment, which this port does not take: the contract wants one integrator "
    "for the union.",
    "5. **Batches on CPU do not amortise.** The cost per draw is flat from 1 to 1024 draws (~10 - 12 ms on the "
    "160 C pot). vmap removes dispatch overhead only, and the arithmetic is serial on one core. A GPU could change "
    "this only with float64 throughput (A100 / H100 class, not T4 / L4); the harness above is ready for that run.",
    "6. **Effort.** The lipid port took about 1 h of agent work: the contract runtime (network.py), the module, the "
    "reference generator, the builder and the tests.",
    "   - The lane is the easiest one: pH-free, 22 linear reactions, no derived-rate chain, two switches.",
    "   - So this measures the per-lane fixed cost (reference generator + tests + builder ~0.5 h), not the "
    "chemistry cost.",
    "   - Remaining A2 + H3 estimates, given this:",
    "     - sulfur: every remaining switch, the pH solver moved into `integrate`, profile mode; 12 - 20 h;",
    "     - trunk: declared-pH scaling, trunk arms, methionine / proline, the a_w steps; 8 - 14 h;",
    "     - acrylamide: 5 - 9 h;",
    "     - network infrastructure: merged titratable table and charge closure, declared competition, forced "
    "inputs, the fit driver; 15 - 25 h;",
    "     - the pH-root performance step: 8 - 15 h.",
    "   - Total ~50 - 85 h, against the bake-off's 75 - 115 h. Confidence that the total lands inside that range: "
    "~0.6. The sulfur switches and the pH merge carry most of the uncertainty.",
    "",
]


if __name__ == "__main__":
    sys.exit(main())

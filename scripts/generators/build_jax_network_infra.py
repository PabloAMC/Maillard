#!/usr/bin/env python
"""
build_jax_network_infra.py -- the A2 network infrastructure in JAX + diffrax, checked (plan step A2 + H3,
docs/guides/A2_MODULE_CONTRACT.md v2.2; 2026-11-05): the forced-input driver (``src/kinetic_core_jax/forcing.py``), the
fit driver (``fit.py``), the declared-competition and titratable-table hooks (``network.compose``) and the bridge to
M1.1's predictive densities (``densities.py``).

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml), never in the core env:

    source <conda.sh>; conda activate maillard-jax
    XLA_FLAGS=--xla_cpu_multi_thread_eigen=false OMP_NUM_THREADS=1 nice -n 10 \
        python scripts/generators/build_jax_network_infra.py [--stage forced ...]

Every compared engine number comes from results/validation/jax_network_infra_reference.json (core env,
``build_jax_network_infra_reference.py``): the engine's ``fitting.fit`` (central-polished) on two SYNTHETIC toy
problems, an engine run with a pH + temperature staircase, and the engine's draws/s. Writes
results/validation/jax_network_infra.{json,md}.

Stages:

* ``forced``: a forced constant equals pH mode "prescribed"; a 'previous' pH + temperature course equals the engine's
  staircase run; a linear pH course and a linear temperature ramp against midpoint staircases (second-order
  convergence); a forced species against the network with its row removed, and the forced-exchange ledger's closure;
  the exact Jacobian through a forced linear course against central differences; the lipid window crossed by a ramp.
* ``fit``: the JAX fit driver on the trunk-lane and lipid toys against the engine's polished optimum (theta in Laplace
  sigmas, cost), its wall time against the engine's, vmap vs loop, the exact Jacobian against central differences at
  the optimum, and the contract's rtol 1e-6 residual check.
* ``draws``: ``densities.Predictor`` draws/s at 64 / 512 / 2048 draws on a trunk-lane toy pot (CPU, one core), against
  the engine's; the batch equals single solves; ``log_likelihoods`` equals ``gaussian_log_density``.
"""
from __future__ import annotations

import argparse
import dataclasses
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

REFERENCE = ROOT / "results" / "validation" / "jax_network_infra_reference.json"
OUT = ROOT / "results" / "validation" / "jax_network_infra"
STAGES = ("forced", "fit", "draws")
TIGHT = {"rtol": 1e-10, "atol": 1e-16}
DEX_FLOOR = 1e-9
DRAW_COUNTS = (64, 512, 2048)
CHUNK = 512


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def load_reference() -> Dict[str, Any]:
    return json.loads(REFERENCE.read_text())


def max_dex(a: np.ndarray, b: np.ndarray, scale: float) -> float:
    """Worst |log10 a - log10 b| over entries of b above DEX_FLOOR x scale."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    m = b > DEX_FLOOR * scale
    return float(np.max(np.abs(np.log10(np.maximum(a[m], 1e-300)) - np.log10(b[m])))) if m.any() else 0.0


# --------------------------------------------------------------------------- networks and problems
_NETS: Dict[str, Any] = {}


def trunk_net():
    from src.kinetic_core_jax import trunk_lane

    if "trunk" not in _NETS:
        _NETS["trunk"] = trunk_lane.compose_lane()
    return _NETS["trunk"]


def lipid_net():
    from src.kinetic_core_jax import lipid, network

    if "lipid" not in _NETS:
        _NETS["lipid"] = network.compose([lipid.build()])
    return _NETS["lipid"]


def toy_problem(toy: Mapping[str, Any], lane: str, batch: str = "vmap", **solve):
    """The reference's toy as a ``fit.JaxProblem`` (same pots, rows, coordinates, bounds; the rest shipped)."""
    from src.kinetic_core_jax import fit, lipid, network, trunk_lane

    if lane == "trunk":
        net = trunk_net()
        pots = []
        for rec in toy["pots"]:
            pot = trunk_lane.make_pot(net, rec["initial"], rec["refined_segments"], ph=rec["ph"],
                                      water_activity=rec["water_activity"])
            expect = np.zeros(net.n_species)
            for k, v in rec["engine_initial"].items():
                expect[net.index[k]] = float(v)
            if not np.array_equal(pot["y0"], expect):
                raise AssertionError(f"{rec['id']}: make_pot's charge differs from the engine's prelude")
            pots.append(pot)
    else:
        net = lipid_net()
        pots = [{"y0": lipid.initial_state(rec["carriers"], net.species),
                 "segments": network.program_segments(rec["refined_segments"])} for rec in toy["pots"]]
    keys = tuple((c["module"], c["coordinate"]) for c in toy["coordinates"])
    obs = tuple(fit.Observation(pot=r["pot"], row=r["row"], observable=r["observable"], value=r["value"],
                                sigma=r["sigma"]) for r in toy["rows"])
    return fit.JaxProblem(net=net, keys=keys, kinds=tuple(c["jax_kind"] for c in toy["coordinates"]),
                          lower=np.asarray(toy["lower"]), upper=np.asarray(toy["upper"]), pots=tuple(pots),
                          observations=obs, solve=solve or {"rtol": 1e-8, "atol": 1e-14}, batch=batch)


# --------------------------------------------------------------------------- forced inputs
def stage_forced(reference) -> Dict[str, Any]:
    import jax
    import jax.numpy as jnp
    from src.kinetic_core_jax import fit, forcing as fz, lipid, network, trunk_lane

    net = trunk_net()
    theta = net.theta0()
    out: Dict[str, Any] = {}

    def solve_plain(pot, **kw):
        return np.asarray(jax.jit(lambda th, p: network.integrate(net, th, p, **TIGHT, **kw).ys)(theta, pot))

    def solve_forced(F, fp, n=net, th=None):
        run = jax.jit(lambda t, p: fz.integrate_forced(n, t, p, F, **TIGHT))(th or theta, fp.pot)
        return np.asarray(run.ys), run

    # 1. a forced constant pH (linear course, constant values) equals pH mode "prescribed"
    segs = [(6.0, 120.0), (14.0, 120.0), (20.0, 140.0)]
    init = {"Glc": 100.0, "Gly": 50.0, "MGO": 5.0}
    pot = trunk_lane.make_pot(net, init, segs, ph=6.4, water_activity=0.9)
    plain = solve_plain(pot)
    F = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"),), ph_scale="label_25C")
    fp = fz.forced_pot(net, F, pot["y0"], {"t_c": fz.program_course(segs), "ph_label": ([0.0, 40.0], [6.4, 6.4])},
                       scalars={"aw": 0.9}, save_times=[6.0, 20.0, 40.0])
    forced, _ = solve_forced(F, fp)
    out["constant_equals_prescribed"] = {"max_dex": max_dex(forced[fp.rows], plain, 100.0),
                                         "max_abs_rel": float(np.max(np.abs(forced[fp.rows] - plain)) / 100.0)}
    # lipid: a forced constant temperature equals the segment mode
    lnet = lipid_net()
    carriers = {"pea_protein_isolate": {}}
    y0 = lipid.initial_state(carriers, lnet.species)
    lsegs = [(10.0, 160.0), (10.0, 160.0)]
    lpot = {"y0": y0, "segments": network.program_segments(lsegs)}
    lplain = np.asarray(network.integrate(lnet, lnet.theta0(), lpot, **TIGHT).ys)
    LF = fz.Forcing(inputs=("t_c",))
    lfp = fz.forced_pot(lnet, LF, y0, {"t_c": ([0.0, 20.0], [160.0, 160.0])}, save_times=[10.0, 20.0])
    lforced = np.asarray(fz.integrate_forced(lnet, lnet.theta0(), lfp.pot, LF, **TIGHT).ys)
    out["lipid_constant_temperature_equals_segments"] = {"max_dex": max_dex(lforced[lfp.rows], lplain,
                                                                            float(np.max(y0)))}

    # 2. a 'previous' pH + temperature course equals the engine's staircase run
    st = reference["forced_staircase"]
    knots = np.concatenate([[0.0], np.cumsum([d for d, _ in st["segments"]])])
    SF = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"), ("ph_label", "previous")),
                    ph_scale="label_25C")
    spot = trunk_lane.make_pot(net, st["initial"], [tuple(s) for s in st["segments"]], ph=st["ph"][0],
                               water_activity=st["water_activity"])
    expect = np.zeros(net.n_species)
    for k, v in st["engine_initial"].items():
        expect[net.index[k]] = float(v)
    assert np.array_equal(spot["y0"], expect), "staircase charge differs from the engine's prelude"
    sfp = fz.forced_pot(net, SF, spot["y0"], {"t_c": fz.program_course([tuple(s) for s in st["segments"]]),
                                              "ph_label": (knots, list(st["ph"]) + [st["ph"][-1]])},
                        scalars={"aw": st["water_activity"]}, save_times=knots[1:])
    sys_, _ = solve_forced(SF, sfp)
    cols = [st["engine_species"].index(s) for s in net.species]
    eng = np.asarray(st["states"])[:, cols]
    out["staircase_vs_engine"] = {"max_dex": max_dex(sys_[sfp.rows], eng, 100.0), "segments": len(st["segments"]),
                                  "ph": st["ph"], "temperatures_c": [t for _, t in st["segments"]],
                                  "engine_tolerance": st["tolerance"], "jax_tolerance": TIGHT}

    # 3. linear courses against midpoint staircases (second-order convergence)
    def convergence(n, F_lin, courses_lin, F_stair, make_stair, scalars, y0_, end, scale):
        lin = fz.forced_pot(n, F_lin, y0_, courses_lin, scalars=scalars, save_times=[end])
        ylin, _ = solve_forced(F_lin, lin, n, n.theta0())
        errs = []
        for n_sub in (2, 4, 8, 16, 32):
            sp = fz.forced_pot(n, F_stair, y0_, make_stair(n_sub), scalars=scalars, save_times=[end])
            ys, _ = solve_forced(F_stair, sp, n, n.theta0())
            errs.append({"n_sub": n_sub, "max_dex": max_dex(ys[sp.rows[-1]], ylin[lin.rows[-1]], scale)})
        ratios = [errs[i]["max_dex"] / max(errs[i + 1]["max_dex"], 1e-300) for i in range(len(errs) - 1)]
        return {"errors": errs, "ratios": ratios, "crossings": lin.crossings}

    ph_knots, ph_vals = [0.0, 50.0], [7.6, 5.4]
    PL = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"),), ph_scale="label_25C")
    PS = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"), ("ph_label", "previous")),
                    ph_scale="label_25C")
    tcourse = ([0.0, 50.0], [125.0, 125.0])
    out["linear_ph_vs_midpoint_staircase"] = convergence(
        net, PL, {"t_c": tcourse, "ph_label": (ph_knots, ph_vals)}, PS,
        lambda n_sub: {"t_c": tcourse, "ph_label": fz.staircase_course(ph_knots, ph_vals, n_sub)},
        {"aw": 0.9}, spot["y0"], 50.0, 100.0)
    t_knots, t_vals = [0.0, 8.0, 20.0, 26.0], [25.0, 160.0, 160.0, 100.0]
    TL = fz.Forcing(inputs=("t_c",))
    TS = fz.Forcing(inputs=("t_c",), interpolation=(("t_c", "previous"),))
    pot_t = trunk_lane.make_pot(net, {"Glc": 100.0, "Fru": 20.0, "Gly": 30.0, "MET": 5.0, "PRO": 5.0, "LYSP": 10.0},
                                [(26.0, 100.0)], ph=5.5, water_activity=0.8)
    out["linear_temperature_vs_midpoint_staircase"] = convergence(
        net, TL, {"t_c": (t_knots, t_vals)}, TS,
        lambda n_sub: {"t_c": fz.staircase_course(t_knots, t_vals, n_sub)},
        {"ph_label": 5.5, "aw": 0.8}, pot_t["y0"], 26.0, 100.0)
    out["lipid_ramp_across_window_vs_midpoint_staircase"] = convergence(
        lnet, LF, {"t_c": ([0.0, 20.0], [60.0, 160.0])}, fz.Forcing(inputs=("t_c",),
                                                                   interpolation=(("t_c", "previous"),)),
        lambda n_sub: {"t_c": fz.staircase_course([0.0, 20.0], [60.0, 160.0], n_sub)}, None, y0, 20.0,
        float(np.max(y0)))

    # 4. a forced species: constant equals the network without its row; the ledger closes
    gly = net.index["Gly"]
    stoich = np.array(net.stoichiometry)
    stoich[gly, :] = 0.0
    net_wo = dataclasses.replace(net, stoichiometry=stoich)
    pot4 = trunk_lane.make_pot(net, {"Glc": 100.0, "Gly": 50.0}, [(30.0, 130.0)], ph=6.8, water_activity=None)
    y_wo = np.asarray(jax.jit(lambda th, p: network.integrate(net_wo, th, p, **TIGHT).ys)(theta, pot4))
    GF = fz.Forcing(inputs=("t_c",), species=("Gly",))
    gfp = fz.forced_pot(net, GF, pot4["y0"], {"t_c": ([0.0, 30.0], [130.0, 130.0]), "Gly": ([0.0, 30.0], [50.0, 50.0])},
                        scalars={"ph_label": 6.8, "aw": float("nan")}, save_times=[30.0])
    yg, run_g = solve_forced(GF, gfp)
    clos_const = fz.closure(net, GF, gfp.pot["y0"], run_g)
    gfp2 = fz.forced_pot(net, GF, pot4["y0"], {"t_c": ([0.0, 30.0], [130.0, 130.0]),
                                               "Gly": ([0.0, 10.0, 30.0], [50.0, 42.0, 30.0])},
                         scalars={"ph_label": 6.8, "aw": float("nan")}, save_times=[10.0, 30.0])
    _, run_g2 = solve_forced(GF, gfp2)
    clos_lin = fz.closure(net, GF, gfp2.pot["y0"], run_g2)
    out["forced_species"] = {
        "constant_vs_row_removed_max_dex": max_dex(yg[gfp.rows], y_wo, 100.0),
        "closure_constant_max": float(np.max(np.abs(clos_const))),
        "closure_linear_course_max": float(np.max(np.abs(clos_lin))),
        "ledger_linear_course_mmol_per_l": np.asarray(run_g2.exchange)[:, 0].tolist(),
        "forced_values_at_rows": np.asarray(run_g2.ys)[:, gly].tolist()}

    # 5. the exact Jacobian through a forced linear pH course against central differences (fit.fd_check)
    keys = (("trunk", "log10_k:k_ama_tdg"), ("trunk", "ea:k_ama_tdg"), ("trunk", "ph_exponent"),
            ("trunk_lane", "pyrazine_ph_slope_above_knot"), ("trunk_lane", "pyrazine_ph_slope_below_knot"))
    lin = fz.forced_pot(net, PL, spot["y0"], {"t_c": tcourse, "ph_label": (ph_knots, ph_vals)}, scalars={"aw": 0.9},
                        save_times=[10.0, 25.0, 50.0])
    obs = tuple(fit.Observation(pot=0, row=r, observable=f"mmol_per_L:{s}", value=0.0, sigma=1.0)
                for r in lin.rows for s in ("TDG", "HMF", "MGO", "PZ", "MPZ"))
    x0 = np.asarray([float(net.theta0()[m][k]) for m, k in keys])
    prob = fit.JaxProblem(net=net, keys=keys, kinds=tuple(fit.kind_of(k) for _, k in keys),
                          lower=x0 - 1.0, upper=x0 + 1.0, pots=(lin.pot,),
                          observations=obs, forcing=PL, solve=TIGHT, batch="loop")
    out["jacobian_through_forced_course"] = {**fit.fd_check(prob, x0), "crossings": lin.crossings,
                                             "coordinates": [f"{m}.{k}" for m, k in keys]}
    return out


# --------------------------------------------------------------------------- the fit driver
def stage_fit(reference) -> Dict[str, Any]:
    from src.kinetic_core_jax import fit

    out: Dict[str, Any] = {}
    for lane in ("trunk", "lipid"):
        toy = reference[f"{lane}_toy"]
        eng = toy["polished"]
        theta_e = np.asarray(list(eng["theta"].values()))
        row: Dict[str, Any] = {"engine": {"theta": theta_e.tolist(), "cost": eng["cost"],
                                          "wall_seconds": eng["wall_seconds"],
                                          "residual_evaluations": eng["residual_evaluations"],
                                          "seconds_per_residual_evaluation": toy["seconds_per_residual_evaluation"],
                                          "status": eng["status"], "polish": eng.get("polish", {}).get("carried"),
                                          "forward_only_cost": toy["forward_only"]["cost"],
                                          "forward_only_theta": list(toy["forward_only"]["theta"].values()),
                                          "forward_only_wall_seconds": toy["forward_only"]["wall_seconds"]},
                               "keys": [f"{c['module']}.{c['coordinate']}" for c in toy["coordinates"]],
                               "truth": toy["truth"], "start": toy["start"]}
        for batch in ("vmap", "loop"):
            print(f"  fit {lane} ({batch}) {uptime()}")
            problem = toy_problem(toy, lane, batch=batch)
            comp = fit.Compiled(problem)
            t0 = time.perf_counter()
            result = fit.fit(problem, toy["start"], budget=200, compiled=comp, require_converged="raise")
            total = time.perf_counter() - t0
            # a second fit from the same start on the warm compile: the steady-state cost of a refit
            t0 = time.perf_counter()
            again = fit.fit(problem, toy["start"], budget=200, compiled=comp, require_converged="raise")
            warm = time.perf_counter() - t0
            lap = fit.laplace(problem, result)
            lap_scaled = fit.laplace(problem, result, scale_by_chi2=True)
            d = np.asarray(result.theta) - theta_e
            row[batch] = {**result.as_dict(problem.keys), "wall_seconds_with_compile": total,
                          "wall_seconds_warm_refit": warm, "warm_refit_cost": again.cost,
                          "sigma": lap["sigma"].tolist(), "sigma_chi2_scaled": lap_scaled["sigma"].tolist(),
                          "reduced_chi2": lap["reduced_chi2"],
                          "delta_theta_vs_engine": d.tolist(),
                          "max_delta_in_sigma": float(np.max(np.abs(d) / lap["sigma"])),
                          "max_delta_in_chi2_scaled_sigma": float(np.max(np.abs(d) / lap_scaled["sigma"])),
                          "cost_relative_difference": float((result.cost - eng["cost"]) / eng["cost"])}
            print(f"    cost {result.cost:.8g} vs engine {eng['cost']:.8g}; max |d theta| / sigma "
                  f"{row[batch]['max_delta_in_sigma']:.2e}; {total:.1f} s (warm {warm:.1f} s) vs engine "
                  f"{eng['wall_seconds']:.1f} s")
            if batch == "vmap":
                # the engine's residual definition at the engine's optimum, evaluated by JAX: the residual gap
                r_at_e = comp.residual(theta_e)[0]
                row["jax_residuals_at_engine_optimum_max_diff"] = float(
                    np.max(np.abs(r_at_e - np.asarray(toy["polished_residuals"]))))
                row["fd_check_at_optimum"] = fit.fd_check(dataclasses.replace(problem, solve=TIGHT), result.theta)
                row["rtol_check"] = fit.residual_tolerance_check(problem, result.theta)
        out[lane] = row
    return out


# --------------------------------------------------------------------------- the density bridge
def stage_draws(reference) -> Dict[str, Any]:
    import jax
    from src.kinetic_core import predictive_densities as pdm
    from src.kinetic_core_jax import densities, network

    toy = reference["trunk_toy"]
    problem = toy_problem(toy, "trunk")
    net = problem.net
    obs = [o for o in problem.observations if o.pot == 0]
    start = np.asarray(toy["start"])
    scale = np.asarray([0.1, 3.0, 0.1, 0.1, 0.05])
    rng = np.random.default_rng(20261107)
    out: Dict[str, Any] = {"pot": toy["pots"][0]["id"], "keys": [f"{m}.{k}" for m, k in problem.keys],
                           "n_rows": len(obs), "chunk": CHUNK, "rows": [], "engine": reference.get("engine_draws")}
    for rtol, atol in ((1e-8, 1e-14), (1e-6, 1e-12)):
        pred = densities.Predictor(net, [problem.pots[0]], obs, problem.keys, chunk=CHUNK, rtol=rtol, atol=atol)
        warm = start + rng.normal(size=(CHUNK, len(start))) * scale
        t0 = time.perf_counter()
        pred(warm)
        compile_s = time.perf_counter() - t0
        small = densities.Predictor(net, [problem.pots[0]], obs, problem.keys, chunk=64, rtol=rtol, atol=atol)
        small(warm[:64])
        for n in DRAW_COUNTS:
            draws = start + rng.normal(size=(n, len(start))) * scale
            p = small if n <= 64 else pred
            t0 = time.perf_counter()
            values, failed = p(draws)
            dt = time.perf_counter() - t0
            out["rows"].append({"rtol": rtol, "atol": atol, "draws": n, "chunk": 64 if n <= 64 else CHUNK,
                                "seconds": dt, "draws_per_second": n / dt, "failed": int(failed.sum()),
                                "compile_and_first_chunk_seconds": compile_s, "uptime": uptime()})
            print(f"  draws {n} rtol {rtol:g}: {n / dt:.1f} draws/s ({int(failed.sum())} failed)")
    # the batch equals single solves; the likelihoods equal gaussian_log_density row by row
    check = start + rng.normal(size=(8, len(start))) * scale
    pred = densities.Predictor(net, [problem.pots[0]], obs, problem.keys, chunk=8)
    values, _ = pred(check)
    base = network.theta_vector(net, net.theta0())
    where = {c: i for i, c in enumerate(net.coordinates)}
    single = jax.jit(lambda v, p: network.integrate(net, network.theta_from_vector(net, v), p, rtol=1e-8,
                                                    atol=1e-14).ys)
    worst = 0.0
    for x, v in zip(check, values):
        vec = np.asarray(base).copy()
        for (m, k), xi in zip(problem.keys, x):
            vec[where[(m, k)]] = xi
        ys = np.asarray(single(vec, problem.pots[0]))
        ref = np.asarray([math.log10(max(float(ys[o.row, net.index[o.observable.split(':')[1]]]), 1e-300))
                          for o in obs])
        worst = max(worst, float(np.max(np.abs(ref - v))))
    y = np.asarray([o.value for o in obs])
    cov = pdm.series_covariance([o.sigma for o in obs], "absolute", omega_dex=0.1)
    ll = densities.log_likelihoods(values, y, cov)
    direct = np.asarray([pdm.gaussian_log_density(y, v, cov) for v in values])
    out["batch_equals_single_solves_max_dex"] = worst
    out["log_likelihoods_vs_gaussian_log_density_max_abs"] = float(np.max(np.abs(ll - direct)))
    return out


# --------------------------------------------------------------------------- the write-up
def _fmt(x: float) -> str:
    return f"{x:.2g}" if x != 0 else "0"


def write_markdown(report: Dict[str, Any]) -> str:
    forced, fits, draws = report["forced"], report["fit"], report["draws"]
    L: List[str] = [
        "# The A2 network infrastructure in JAX + diffrax: forced inputs, the fit driver, the density bridge",
        "",
        f"*Generated by `scripts/generators/build_jax_network_infra.py` in the separate `maillard-jax` env (jax "
        f"{report['jax_version']}, diffrax {report['diffrax_version']}); engine numbers from "
        "`jax_network_infra_reference.json` (core env, `build_jax_network_infra_reference.py`). Plan step A2 + H3; "
        "contract v2.2 ([A2_MODULE_CONTRACT.md](../../docs/guides/A2_MODULE_CONTRACT.md)). No shipped number moves. "
        "The toy problems are SYNTHETIC: software equivalence only, no evidence about chemistry.*",
        "",
        "## 1. Forced inputs (`src/kinetic_core_jax/forcing.py`)",
        "",
        "| check | result |",
        "|---|---|",
        f"| forced constant pH == pH mode \"prescribed\" (trunk lane, 3 segments, rtol 1e-10) | "
        f"{_fmt(forced['constant_equals_prescribed']['max_dex'])} dex |",
        f"| forced constant temperature == segment mode (lipid, 160 C) | "
        f"{_fmt(forced['lipid_constant_temperature_equals_segments']['max_dex'])} dex |",
        f"| 'previous' pH + temperature course == the engine's staircase run ({forced['staircase_vs_engine']['segments']}"
        f" segments, pH {forced['staircase_vs_engine']['ph'][0]} -> {forced['staircase_vs_engine']['ph'][-1]}, "
        f"110 - 130 C; engine LSODA rtol 1e-12) | {_fmt(forced['staircase_vs_engine']['max_dex'])} dex |",
        f"| forced constant Gly == the network with Gly's row removed | "
        f"{_fmt(forced['forced_species']['constant_vs_row_removed_max_dex'])} dex |",
        f"| element closure with the forced-exchange ledger (constant / linear Gly course) | "
        f"{_fmt(forced['forced_species']['closure_constant_max'])} / "
        f"{_fmt(forced['forced_species']['closure_linear_course_max'])} (relative to E . y0) |",
        f"| exact Jacobian through a linear pH course (crossing B18's knot at pH 7) vs column-scaled central "
        f"differences, 5 coordinates x 15 rows | {_fmt(forced['jacobian_through_forced_course']['max_relative_column_error'])}"
        " (max column-relative) |",
        "",
        "Linear courses against midpoint staircases with n steps per interval (max dex at the end; a second-order "
        "rule should give ratios near 4):",
        "",
        "| course | n = 2 | 4 | 8 | 16 | 32 | ratios |",
        "|---|---|---|---|---|---|---|",
    ]
    for key, label in (("linear_ph_vs_midpoint_staircase", "pH 7.6 -> 5.4 over 50 min at 125 C (trunk lane)"),
                       ("linear_temperature_vs_midpoint_staircase", "T 25 -> 160 -> 160 -> 100 C (trunk lane)"),
                       ("lipid_ramp_across_window_vs_midpoint_staircase", "T 60 -> 160 C over 20 min (lipid, crosses "
                        "the 100 C window)")):
        c = forced[key]
        L.append(f"| {label} | " + " | ".join(_fmt(e["max_dex"]) for e in c["errors"]) + " | "
                 + ", ".join(f"{r:.2f}" for r in c["ratios"]) + " |")
    L += ["", "## 2. The fit driver (`src/kinetic_core_jax/fit.py`) against `fitting.fit` (central-polished)", "",
          "| toy | engine cost | JAX cost | rel. diff | max abs(d theta) / sigma (raw; chi2-scaled) | engine wall s "
          "(evaluations) | JAX wall s, vmap (compile; warm refit) | JAX loop | iterations (nfev / njev) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for lane in ("trunk", "lipid"):
        f = fits[lane]
        v, lp = f["vmap"], f["loop"]
        comp = sum(x for x in (v["compile_seconds"]["residual"] or 0.0, v["compile_seconds"]["jacobian"] or 0.0))
        L.append(f"| {lane} ({len(f['keys'])} coordinates) | {f['engine']['cost']:.6g} | {v['cost']:.6g} | "
                 f"{v['cost_relative_difference']:.1e} | {v['max_delta_in_sigma']:.1e}; "
                 f"{v['max_delta_in_chi2_scaled_sigma']:.1e} | {f['engine']['wall_seconds']:.1f} "
                 f"({f['engine']['residual_evaluations']}) | {v['wall_seconds_with_compile']:.1f} ({comp:.1f}; "
                 f"{v['wall_seconds_warm_refit']:.1f}) | {lp['wall_seconds_with_compile']:.1f} (warm "
                 f"{lp['wall_seconds_warm_refit']:.1f}) | {v['nfev']} / {v['njev']} |")
    L += ["", "At the optimum: exact Jacobian vs column-scaled central differences (rtol 1e-10), the residuals' move "
          "between rtol 1e-6 and 1e-8 (contract sec. 2.6), and the engine's forward-only fit (no polish):", ""]
    for lane in ("trunk", "lipid"):
        f = fits[lane]
        fo = np.asarray(f["engine"]["forward_only_theta"])
        d_fo = np.max(np.abs(fo - np.asarray(f["engine"]["theta"])) / np.asarray(f["vmap"]["sigma"]))
        L.append(f"- **{lane}**: FD check {_fmt(f['fd_check_at_optimum']['max_relative_column_error'])}; rtol move "
                 f"{_fmt(f['rtol_check']['max_move'])} dex ({'passes' if f['rtol_check']['passes'] else 'FAILS'} "
                 f"< 1e-3); JAX residuals at the engine's optimum differ from the engine's by "
                 f"{_fmt(f['jax_residuals_at_engine_optimum_max_diff'])} (sigma units); the engine's forward-only "
                 f"optimum is {d_fo:.2g} sigma from its polished one (cost {f['engine']['forward_only_cost']:.6g}, "
                 f"{f['engine']['forward_only_wall_seconds']:.1f} s).")
    eng = draws.get("engine") or {}
    L += ["", "## 3. The density bridge (`src/kinetic_core_jax/densities.py`): draws/s on a trunk-lane toy pot", "",
          f"Pot `{draws['pot']}` (55 steps, 3 segments), {draws['n_rows']} rows, {len(draws['keys'])} drawn "
          "coordinates; CPU, one core, niced; shared machine (load beside each row).", "",
          "| rtol | draws | chunk | seconds | draws/s | failed | load |", "|---|---|---|---|---|---|---|"]
    for r in draws["rows"]:
        load = r["uptime"].split("load averages:")[-1].strip()
        L.append(f"| {r['rtol']:g} | {r['draws']} | {r['chunk']} | {r['seconds']:.1f} | {r['draws_per_second']:.1f} | "
                 f"{r['failed']} | {load} |")
    for r in eng.get("rows", []):
        L.append(f"| engine LSODA {eng['tolerance']['rtol']:g} | {r['draws']} | - | {r['seconds']:.1f} | "
                 f"{r['draws_per_second']:.1f} | - | {eng['load_before'].split('load averages:')[-1].strip()} |")
    L += ["", f"The batch equals single solves to {_fmt(draws['batch_equals_single_solves_max_dex'])} dex; "
          f"`log_likelihoods` equals `gaussian_log_density` row by row to "
          f"{_fmt(draws['log_likelihoods_vs_gaussian_log_density_max_abs'])} nats.", "",
          "## 4. Findings", ""]
    L += [f"{i}. {t}" for i, t in enumerate(report.get("findings", []), 1)]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--stage", action="append", choices=STAGES)
    args = parser.parse_args(argv)
    import diffrax
    import jax

    reference = load_reference()
    stages = args.stage or list(STAGES)
    path = OUT.with_suffix(".json")
    report = json.loads(path.read_text()) if path.exists() and args.stage else {}
    report.update({"artifact": "jax_network_infra", "generated_by": "scripts/generators/build_jax_network_infra.py",
                   "env": "maillard-jax", "jax_version": jax.__version__, "diffrax_version": diffrax.__version__})
    for stage in stages:
        print(f"[{stage}] {uptime()}")
        report[stage] = globals()[f"stage_{stage}"](reference)
        path.write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    report["findings"] = FINDINGS
    path.write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    if all(k in report for k in STAGES):
        OUT.with_suffix(".md").write_text(write_markdown(report), encoding="utf-8")
    print(f"wrote {path}")
    return 0


#: the write-up's closing section, edited by hand after reading the numbers (kept here so the .md regenerates)
FINDINGS: List[str] = []


if __name__ == "__main__":
    raise SystemExit(main())

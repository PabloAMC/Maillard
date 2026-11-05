#!/usr/bin/env python
"""
profile_fit_evaluation.py -- plan step W1.3: where one fit evaluation's time goes (2026-11-05).

One evaluation of M1's objective (``fit_m1_supply_combinations.Objective``: R0 on B60's 52 rows; C1 on B60's rows
plus the route block), on the inline pH coupling, at the carried optima
(``results/validation/kinetic_core_m1_supply_combinations_carried_{R0,C1}.json``), timed and taken apart:

  --part scaling    wall time at MAILLARD_WORKERS = 1, 2, 4, 8 (min of ``--repeats``, interleaved), every pool
                    task's start / end / process, so per-pot times, stragglers, idle and pool overhead follow
  --part inside     one slow pot and one typical pot, single process: timer-wrapped engine functions (the
                    shares) and cProfile (call counts: RHS, Jacobian, integrator steps)
  --part jacobian   one forward-difference Jacobian of the objective (17 / 19 coordinates + centre) as the
                    fit runs it (``residuals_many``: points across workers, pots serial inside a worker)
  --part fixed      R0's evaluation on ``ph.coupling = "fixed_point_nodes"`` (B60's quick setting) at 4 workers
  --part percall    per-call costs on the RHS path (rate constants, pH solve and its numba kernel, mass action,
                    analytic Jacobian) at captured states of one pot
  --part memo       a bit-identical candidate MEASURED, not installed: the pure functions of temperature memoised
  --part render     the .md from the .json (with the findings and the ranked speedups declared below)

Nothing is changed in the engine (the memo candidate is installed only inside its own measurement): the timing wrappers sit around ``parallel.map_ordered``'s tasks and the engine's
module functions, and the residual vector of every timed evaluation is checked against the first one (bit-identical)
and its cost against the carried record. Each measurement records the 1-min load average before it; above
``LOAD_LIMIT`` it waits ``LOAD_WAIT_S`` and retries up to three times, then runs and labels the number.

Writes results/validation/profile_fit_evaluation_2026-11-05.{json,md}.
"""
from __future__ import annotations

import argparse
import contextlib
import cProfile
import json
import os
import pickle
import pstats
import statistics
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import fit_m1_supply_combinations as M1  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import config, parallel, ph_state, sulfur  # noqa: E402

GENERATED_BY = "scripts/generators/profile_fit_evaluation.py"
STEM = "profile_fit_evaluation_2026-11-05"
OUT_JSON = data_paths.VALIDATION_DIR / f"{STEM}.json"
OUT_MD = data_paths.VALIDATION_DIR / f"{STEM}.md"
LOAD_LIMIT = 4.0
LOAD_WAIT_S = 90
LOAD_RETRIES = 3
WORKER_COUNTS = (1, 2, 4, 8)

# ----------------------------------------------------------------------------------------------------------------
# The objective at a carried optimum
# ----------------------------------------------------------------------------------------------------------------


def objective_at(rival: str, numerics: Optional[M1.Numerics] = None):
    checks = M1._read(M1.CHECKS_OUT)
    car = M1._read(M1.carried_path(rival))
    numerics = numerics or M1.numerics_for_fits(checks)
    obj = M1.Objective(M1.RIVALS[rival], numerics, M1.fit_side_from_checks(checks), car.get("carried_noise") or {})
    theta = np.array([car["carried_theta"][k] for k in obj.keys], dtype=float)
    return obj, theta, float(car["carried_cost"])


# ----------------------------------------------------------------------------------------------------------------
# Load discipline
# ----------------------------------------------------------------------------------------------------------------


def load_gate(label: str) -> Dict[str, Any]:
    tries = []
    for attempt in range(LOAD_RETRIES + 1):
        one, five, fifteen = os.getloadavg()
        tries.append(round(one, 2))
        if one <= LOAD_LIMIT or attempt == LOAD_RETRIES:
            break
        print(f"[{label}] load {one:.2f} > {LOAD_LIMIT}: waiting {LOAD_WAIT_S} s", flush=True)
        time.sleep(LOAD_WAIT_S)
    one, five, fifteen = os.getloadavg()
    return {"load_1min": round(one, 2), "load_5min": round(five, 2), "load_15min": round(fifteen, 2),
            "load_tries": tries, "clean": one <= LOAD_LIMIT,
            "label": None if one <= LOAD_LIMIT else f"LOADED: 1-min load {one:.2f} > {LOAD_LIMIT} after retries"}


# ----------------------------------------------------------------------------------------------------------------
# The task timeline: every parallel.map_ordered task, timed in the process that runs it
# ----------------------------------------------------------------------------------------------------------------

_LOG: Dict[str, Any] = {"dir": None, "call": 0}
_ORIGINAL_MAP = parallel.map_ordered


def _write_line(record: Dict[str, Any]) -> None:
    d = _LOG["dir"]
    if d is None:
        return
    with open(os.path.join(d, f"{os.getpid()}.jsonl"), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")


def _label(phase: str, item: Any, i: int) -> str:
    if phase == "_integrate_system":
        return f"b60:{item[0]}"
    if phase == "budget_of":
        return f"flux:{item}"
    if phase == "_run_pot_safe":
        return f"route:{item.get('pot')}"
    return f"{phase}:{i}"


def _timed_map(fn, items, n_workers=None, *, config=None):
    items = list(items)
    _LOG["call"] += 1
    call = f"{os.getpid()}:{_LOG['call']}"
    phase = getattr(fn, "__name__", "fn")
    wrapped = []
    for i, it in enumerate(items):
        if isinstance(it, parallel.Configured):
            wrapped.append(parallel.Configured((i, it.item), it.config))
        else:
            wrapped.append((i, it))

    def timed(pair):
        i, it = pair
        t0, c0 = time.time(), time.process_time()
        out = fn(it)
        t1, c1 = time.time(), time.process_time()
        _write_line({"kind": "task", "call": call, "phase": phase, "i": i, "label": _label(phase, it, i),
                     "pid": os.getpid(), "t0": t0, "t1": t1, "cpu": c1 - c0})
        return out

    timed.__name__ = phase
    t0 = time.time()
    out = _ORIGINAL_MAP(timed, wrapped, n_workers, config=config)
    t1 = time.time()
    n = parallel.workers() if n_workers is None else max(1, int(n_workers))
    serial = n == 1 or len(items) < 2 or parallel._IN_WORKER
    _write_line({"kind": "call", "call": call, "phase": phase, "n_items": len(items), "workers": 1 if serial else
                 min(n, len(items)), "pid": os.getpid(), "in_worker": bool(parallel._IN_WORKER), "t0": t0, "t1": t1})
    return out


class Timeline:
    """Install the timed map for a block; read every process's records after it."""

    def __enter__(self):
        self.dir = tempfile.mkdtemp(prefix="w13_timeline_")
        _LOG["dir"] = self.dir
        parallel.map_ordered = _timed_map
        return self

    def __exit__(self, *exc):
        parallel.map_ordered = _ORIGINAL_MAP
        _LOG["dir"] = None
        recs = []
        for f in Path(self.dir).glob("*.jsonl"):
            recs += [json.loads(line) for line in f.read_text(encoding="utf-8").splitlines() if line]
            f.unlink()
        os.rmdir(self.dir)
        self.tasks = sorted([r for r in recs if r["kind"] == "task"], key=lambda r: r["t0"])
        self.calls = sorted([r for r in recs if r["kind"] == "call"], key=lambda r: r["t0"])
        return False


def lpt_makespan(durations: Sequence[float], workers: int) -> float:
    """Longest-processing-time-first list scheduling's makespan (the dynamic pool with the queue sorted)."""
    loads = [0.0] * workers
    for d in sorted(durations, reverse=True):
        j = min(range(workers), key=loads.__getitem__)
        loads[j] += d
    return max(loads)


def queue_makespan(durations: Sequence[float], workers: int) -> float:
    """The pool's own order (chunksize 1: the next task to the first free worker), no overhead."""
    loads = [0.0] * workers
    for d in durations:
        j = min(range(workers), key=loads.__getitem__)
        loads[j] += d
    return max(loads)


def analyse_calls(tl: Timeline, wall: float) -> Dict[str, Any]:
    """Per pool call (phase): items, workers, wall, busy sum, the slowest task, the pool's start / tail latency,
    per-worker busy fraction, and the makespans of the pool's order and of longest-first with no overhead."""
    by_call = defaultdict(list)
    for t in tl.tasks:
        by_call[t["call"]].append(t)
    out = []
    for c in tl.calls:
        if c["in_worker"]:
            continue
        ts = sorted(by_call.get(c["call"], []), key=lambda r: r["i"])
        if not ts:
            continue
        d = [t["t1"] - t["t0"] for t in ts]
        window = c["t1"] - c["t0"]
        busy_by_pid = defaultdict(float)
        for t in ts:
            busy_by_pid[t["pid"]] += t["t1"] - t["t0"]
        w = int(c["workers"])
        out.append({
            "phase": c["phase"], "n_items": c["n_items"], "workers": w, "wall_s": round(window, 4),
            "task_sum_s": round(sum(d), 4), "task_cpu_sum_s": round(sum(t.get("cpu", 0.0) for t in ts), 4),
            "slowest_task_s": round(max(d), 4),
            "slowest_task": ts[int(np.argmax(d))]["label"],
            "start_latency_s": round(min(t["t0"] for t in ts) - c["t0"], 4),
            "tail_latency_s": round(c["t1"] - max(t["t1"] for t in ts), 4),
            "busy_fraction_per_worker": sorted([round(b / window, 3) for b in busy_by_pid.values()], reverse=True),
            "mean_busy_fraction": round(sum(d) / (w * window), 3) if window > 0 else None,
            "makespan_pool_order_no_overhead_s": round(queue_makespan(d, w), 4),
            "makespan_longest_first_no_overhead_s": round(lpt_makespan(d, w), 4),
            "lower_bound_s": round(max(max(d), sum(d) / w), 4),
        })
    return {"calls": out, "wall_s": round(wall, 4),
            "pool_calls_wall_s": round(sum(c["wall_s"] for c in out), 4),
            "outside_pool_s": round(wall - sum(c["wall_s"] for c in out), 4)}


def task_table(tl: Timeline) -> List[Dict[str, Any]]:
    return [{"phase": t["phase"], "label": t["label"], "seconds": round(t["t1"] - t["t0"], 4),
             "cpu_s": round(t.get("cpu", float("nan")), 4)} for t in tl.tasks]


# ----------------------------------------------------------------------------------------------------------------
# --part scaling
# ----------------------------------------------------------------------------------------------------------------


def one_evaluation(obj, theta, workers: int):
    os.environ["MAILLARD_WORKERS"] = str(workers)
    with Timeline() as tl:
        t0 = time.time()
        r = obj.residuals(theta)
        wall = time.time() - t0
    return r, wall, tl


def part_scaling(rival: str, repeats: int) -> Dict[str, Any]:
    obj, theta, carried_cost = objective_at(rival)
    reference = None
    runs: Dict[int, List[Dict[str, Any]]] = {w: [] for w in WORKER_COUNTS}
    serial_tasks = None
    t0 = time.time()
    _, cold_wall, _ = one_evaluation(obj, theta, 1)   # warm-up: numba's cache load, first-call setup (untimed)
    cold = round(time.time() - t0, 4)
    for rep in range(repeats):
        for w in (WORKER_COUNTS if rep % 2 == 0 else tuple(reversed(WORKER_COUNTS))):
            gate = load_gate(f"{rival} W={w} rep {rep}")
            del M1.FAILURES[:]
            r, wall, tl = one_evaluation(obj, theta, w)
            if reference is None:
                reference = r
            identical = bool(np.array_equal(r, reference))
            rec = {"repeat": rep, "wall_s": round(wall, 4), "identical_to_first": identical,
                   "failed_pots": len(M1.FAILURES), **gate, "analysis": analyse_calls(tl, wall)}
            if w == 1 and (serial_tasks is None or wall < min(x["wall_s"] for x in runs[1])):
                serial_tasks = task_table(tl)
            runs[w].append(rec)
            print(f"{rival} W={w} rep {rep}: {wall:.2f} s (load {gate['load_1min']}), identical {identical}",
                  flush=True)
    cost = float(0.5 * reference @ reference)
    best = {w: min(runs[w], key=lambda x: x["wall_s"]) for w in WORKER_COUNTS}
    t1 = best[1]["wall_s"]
    durations = [t["seconds"] for t in serial_tasks]
    total = sum(durations)
    ranked = sorted(serial_tasks, key=lambda t: -t["seconds"])
    by_phase = defaultdict(list)
    for t in serial_tasks:
        by_phase[t["phase"]].append(t["seconds"])
    return {
        "rival": rival, "n_free": len(obj.keys), "cold_first_evaluation_s": cold, "n_residuals": int(len(reference)), "cost": cost,
        "carried_cost": carried_cost, "cost_matches_carried": bool(abs(cost - carried_cost) <= 1e-6 * carried_cost),
        "numerics": dict(M1.dataclasses.asdict(obj.numerics)),
        "wall_by_workers": {str(w): {"min_s": best[w]["wall_s"], "all_s": [x["wall_s"] for x in runs[w]],
                                     "speedup_vs_1": round(t1 / best[w]["wall_s"], 3),
                                     "efficiency": round(t1 / best[w]["wall_s"] / w, 3),
                                     "load_1min": [x["load_1min"] for x in runs[w]],
                                     "labels": [x["label"] for x in runs[w] if x["label"]]}
                            for w in WORKER_COUNTS},
        "best_run_analysis": {str(w): best[w]["analysis"] for w in WORKER_COUNTS},
        "all_identical": all(x["identical_to_first"] for w in WORKER_COUNTS for x in runs[w]),
        "failed_pots": sum(x["failed_pots"] for w in WORKER_COUNTS for x in runs[w]),
        "serial_pots": {
            "n_tasks": len(durations), "task_sum_s": round(total, 4),
            "per_phase": {p: {"n": len(v), "sum_s": round(sum(v), 4), "median_s": round(statistics.median(v), 4),
                              "max_s": round(max(v), 4)} for p, v in by_phase.items()},
            "quantiles_s": {q: round(float(np.quantile(durations, q / 100)), 4) for q in (10, 25, 50, 75, 90, 100)},
            "slowest_10": [dict(t, share=round(t["seconds"] / total, 4)) for t in ranked[:10]],
            "slowest_10_share": round(sum(t["seconds"] for t in ranked[:10]) / total, 4),
            "slowest_1_share": round(ranked[0]["seconds"] / total, 4),
            "all": serial_tasks,
        },
    }


# ----------------------------------------------------------------------------------------------------------------
# --part inside: one pot, single process
# ----------------------------------------------------------------------------------------------------------------

#: the engine functions timed by wrapper; the closures (rhs, jac, constants) are reached through cProfile
TIMED_FUNCTIONS = {
    "sulfur": ("sulfur_rate_constants_at", "rate_vector", "_rates_vectorised", "_rate_jacobian",
               "sulfur_reaction_rates", "reaction_flux_integrals", "rate_vector_at"),
    "ph_state": ("solve_ph_and_ionic_strength_compiled", "solve_ph_and_ionic_strength_warm",
                 "solve_ph_and_ionic_strength", "titratable_inventory", "ph_and_ionic_strength_of_state",
                 "ph_of_state", "strong_ion_difference"),
}


class FunctionTimers:
    """Wrap module functions with an exclusive timer (time inside a wrapped callee is not counted twice)."""

    def __enter__(self):
        self.saved, self.acc, self.n, self.stack = [], defaultdict(float), defaultdict(int), []
        mods = {"sulfur": sulfur, "ph_state": ph_state}
        for mname, names in TIMED_FUNCTIONS.items():
            mod = mods[mname]
            for name in names:
                if not hasattr(mod, name):
                    continue
                fn = getattr(mod, name)
                self.saved.append((mod, name, fn))
                setattr(mod, name, self._wrap(f"{mname}.{name}", fn))
        return self

    def _wrap(self, key, fn):
        acc, n, stack = self.acc, self.n, self.stack

        def wrapper(*a, **kw):
            t0 = time.perf_counter()
            stack.append(0.0)
            try:
                return fn(*a, **kw)
            finally:
                inner = stack.pop()
                dt = time.perf_counter() - t0
                acc[key] += dt - inner
                n[key] += 1
                if stack:
                    stack[-1] += dt
        return wrapper

    def __exit__(self, *exc):
        for mod, name, fn in self.saved:
            setattr(mod, name, fn)
        return False


def pot_runner(obj, theta, label: str) -> Callable[[], Any]:
    """A zero-argument call that integrates one pot exactly as the objective does (its Config installed)."""
    point = obj.point(theta)
    cfg = M1.config_of(point, obj.numerics, obj.fit_side)
    kind, name = label.split(":", 1)
    x = point.x()
    if kind == "b60":
        parameters, drift = M1.B23.build_parameters(x), M1.B23.unpack(x)[3]
        item = (name, M1.B23.SYSTEMS[name], parameters, drift, obj.numerics)
        return lambda: _under(cfg, lambda: M1._integrate_system(item))
    if kind == "flux":
        parameters, drift = M1.B23.build_parameters(x), M1.B23.unpack(x)[3]
        spec = M1.B23.SYSTEMS[name]
        return lambda: _under(cfg, lambda: sulfur.sulfur_flux_budget(
            parameters, float(spec["t_c"]) + M1.CELSIUS, spec["initial"], float(spec["minutes"]),
            ph=float(spec["ph"]), buffer_spec=spec.get("buffer"), ph_drift=drift,
            n_points=21 if obj.numerics.quick else 101))
    if kind == "route":
        (job,) = [j for p, j in M1.route_block()[1].items() if p == name]
        (task,) = M1.table_tasks(point, [job], obj.numerics)
        return lambda: _under(cfg, lambda: M1._run_pot_safe(task))
    raise ValueError(label)


def _under(cfg, fn):
    with config.applied(cfg):
        return fn()


def classify_profile(stats: pstats.Stats) -> Dict[str, Any]:
    """cProfile's call counts and inclusive times of the pieces that matter."""
    want = {"rhs": "sulfur.py", "jac": "sulfur.py", "constants": "sulfur.py", "state_ph": "sulfur.py",
            "_step_impl": "lsoda.py", "solve_ivp": "ivp.py", "_integrate_inline_ph": "sulfur.py",
            "_solve_on_time_nodes": "sulfur.py", "sulfur_rate_constants_at": "sulfur.py",
            "_rates_vectorised": "sulfur.py", "_rate_jacobian": "sulfur.py", "rate_vector": "sulfur.py",
            "titratable_inventory": "ph_state.py", "solve_ph_and_ionic_strength_compiled": "ph_state.py",
            "integrate_sulfur": "sulfur.py", "reaction_flux_integrals": "sulfur.py",
            "sulfur_reaction_rates": "sulfur.py", "rate_vector_at": "sulfur.py", "dense_output": "base.py",
            "_call_impl": "lsoda.py", "__call__": "common.py"}
    out = {}
    for (file, line, func), (cc, nc, tt, ct, callers) in stats.stats.items():
        if func in want and file.endswith(want[func]):
            key = func if func not in out else f"{func}@{line}"
            out[key] = {"calls": nc, "tottime_s": round(tt, 4), "cumtime_s": round(ct, 4),
                        "callers": {f"{c[2]}": {"calls": v[1], "cumtime_s": round(v[3], 4)}
                                    for c, v in callers.items()}}
    builtins = []
    for (file, line, func), (cc, nc, tt, ct, callers) in stats.stats.items():
        if file == "~":
            builtins.append({"func": func, "calls": nc, "tottime_s": round(tt, 4)})
    builtins.sort(key=lambda b: -b["tottime_s"])
    top = sorted(((k, v) for k, v in stats.stats.items()), key=lambda kv: -kv[1][2])[:25]
    return {"functions": out, "top_builtins_by_tottime": builtins[:12],
            "top_by_tottime": [{"func": f"{Path(k[0]).name}:{k[1]}:{k[2]}", "calls": v[1], "tottime_s": round(v[2], 4),
                                "cumtime_s": round(v[3], 4)} for k, v in top]}


def part_inside(rival: str, labels: Sequence[str], repeats: int) -> Dict[str, Any]:
    obj, theta, _ = objective_at(rival)
    os.environ["MAILLARD_WORKERS"] = "1"
    out = {}
    for label in labels:
        run = pot_runner(obj, theta, label)
        run()   # warm (numba cache, imports)
        gate = load_gate(f"inside {label}")
        walls = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            run()
            walls.append(time.perf_counter() - t0)
        with FunctionTimers() as ft:
            t0 = time.perf_counter()
            run()
            wrapped_wall = time.perf_counter() - t0
        prof = cProfile.Profile()
        t0 = time.perf_counter()
        prof.enable()
        result = run()
        prof.disable()
        prof_wall = time.perf_counter() - t0
        st = pstats.Stats(prof)
        t0 = time.perf_counter()
        blob = pickle.dumps(result, protocol=pickle.HIGHEST_PROTOCOL)
        pickle_s = time.perf_counter() - t0
        timers = {k: {"seconds": round(v, 5), "calls": ft.n[k], "share_of_wrapped_wall": round(v / wrapped_wall, 4)}
                  for k, v in sorted(ft.acc.items(), key=lambda kv: -kv[1])}
        out[label] = {"wall_s_min": round(min(walls), 4), "wall_s_all": [round(w, 4) for w in walls],
                      "wall_s_with_timers": round(wrapped_wall, 4), "wall_s_under_cprofile": round(prof_wall, 4),
                      "timers_exclusive": timers, "result_pickle_bytes": len(blob),
                      "result_pickle_s": round(pickle_s, 5), **gate, "cprofile": classify_profile(st)}
        print(f"inside {label}: {min(walls):.3f} s", flush=True)
    return out


# ----------------------------------------------------------------------------------------------------------------
# --part jacobian / --part fixed
# ----------------------------------------------------------------------------------------------------------------


def part_jacobian(rival: str, workers: int) -> Dict[str, Any]:
    obj, theta, _ = objective_at(rival)
    problem = obj.problem()
    os.environ["MAILLARD_WORKERS"] = str(workers)
    from src.kinetic_core import fitting

    gate = load_gate(f"jacobian {rival}")
    with Timeline() as tl:
        t0 = time.time()
        jac = fitting._forward_jacobian(problem, theta, problem.lower, problem.upper,
                                        lambda pts: fitting._evaluate_all(problem, pts, None))
        wall = time.time() - t0
    points = [t for t in tl.tasks if t["phase"] == "_evaluate_job"]
    (call,) = [c for c in tl.calls if c["phase"] == "_evaluate_job" and not c["in_worker"]]
    by_pid = defaultdict(list)
    for t in points:
        by_pid[t["pid"]].append(t)
    window = call["t1"] - call["t0"]
    per_worker = []
    for pid, ts in sorted(by_pid.items(), key=lambda kv: min(t["t0"] for t in kv[1])):
        busy = sum(t["t1"] - t["t0"] for t in ts)
        per_worker.append({"points": [t["i"] for t in sorted(ts, key=lambda r: r["t0"])], "busy_s": round(busy, 3),
                           "idle_s": round(window - busy, 3), "busy_fraction": round(busy / window, 3),
                           "first_start_s": round(min(t["t0"] for t in ts) - call["t0"], 3),
                           "last_end_s": round(max(t["t1"] for t in ts) - call["t0"], 3)})
    d = [t["t1"] - t["t0"] for t in sorted(points, key=lambda r: r["i"])]
    cpu = [t["cpu"] for t in sorted(points, key=lambda r: r["i"])]
    return {"rival": rival, "workers": workers, "n_points": len(points), "wall_s": round(wall, 3),
            "env_threads": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
            "point_cpu_seconds": [round(x, 3) for x in cpu],
            "cpu_over_wall": round(sum(cpu) / sum(d), 3),
            "pool_window_s": round(window, 3), "point_seconds": [round(x, 3) for x in d],
            "point_mean_s": round(float(np.mean(d)), 3), "point_cv": round(float(np.std(d) / np.mean(d)), 3),
            "rounds": int(np.ceil(len(points) / workers)), "per_worker": per_worker,
            "mean_busy_fraction": round(sum(d) / (workers * window), 3),
            "makespan_pool_order_no_overhead_s": round(queue_makespan(d, workers), 3),
            "makespan_longest_first_no_overhead_s": round(lpt_makespan(d, workers), 3),
            "lower_bound_s": round(sum(d) / workers, 3), "jacobian_finite": bool(np.all(np.isfinite(jac))), **gate}


def part_fixed(workers: int) -> Dict[str, Any]:
    obj, theta, _ = objective_at("R0", M1.fit_fixed_point())
    gate = load_gate("fixed point R0")
    del M1.FAILURES[:]
    r, wall, tl = one_evaluation(obj, theta, workers)
    tasks = task_table(tl)
    ranked = sorted(tasks, key=lambda t: -t["seconds"])
    return {"rival": "R0", "numerics": dict(M1.dataclasses.asdict(obj.numerics)), "workers": workers,
            "wall_s": round(wall, 4), "cost": float(0.5 * r @ r), "failed_pots": len(M1.FAILURES),
            "analysis": analyse_calls(tl, wall), "task_sum_s": round(sum(t["seconds"] for t in tasks), 3),
            "slowest_5": ranked[:5], **gate}


# ----------------------------------------------------------------------------------------------------------------
# --part percall: per-call costs on the RHS path at captured states of one pot
# ----------------------------------------------------------------------------------------------------------------


def part_percall(rival: str, label: str, n: int = 2000) -> Dict[str, Any]:
    """The arguments of the middle call of each engine function, captured during one run of ``label``, re-timed
    in isolation (timeit, min of 3): rate constants, rate_vector, the pH solve (wrapper + numba kernel) and the
    kernel alone, the mass-action RHS and the analytic Jacobian."""
    import timeit

    from src.kinetic_core import ph_kernel

    obj, theta, _ = objective_at(rival)
    os.environ["MAILLARD_WORKERS"] = "1"
    run = pot_runner(obj, theta, label)
    run()
    captured: Dict[str, list] = defaultdict(list)
    saved = (sulfur.sulfur_rate_constants_at, ph_state.solve_ph_and_ionic_strength_compiled)

    def capture(key, fn):
        def wrapper(*a, **kw):
            captured[key].append((a, kw))
            return fn(*a, **kw)
        return wrapper

    sulfur.sulfur_rate_constants_at = capture("rate_constants", saved[0])
    ph_state.solve_ph_and_ionic_strength_compiled = capture("ph_solve", saved[1])
    try:
        run()
    finally:
        sulfur.sulfur_rate_constants_at, ph_state.solve_ph_and_ionic_strength_compiled = saved
    gate = load_gate(f"percall {label}")
    cfg = M1.config_of(obj.point(theta), obj.numerics, obj.fit_side)

    def per_call(fn, number=n):
        return min(timeit.repeat(fn, number=number, repeat=3)) / number

    with config.applied(cfg):
        a, kw = captured["rate_constants"][len(captured["rate_constants"]) // 2]
        k = saved[0](*a, **kw)
        kv = sulfur.rate_vector(k)
        sa, skw = captured["ph_solve"][len(captured["ph_solve"]) // 2]
        temperature_k, sid, inventory, guess = sa[0], sa[1], sa[2], sa[3]
        ionic = sa[4] if len(sa) > 4 else 0.0
        groups = [(g, c) for g, c in inventory if c > 0.0]
        width = max((len(g.pkas) for g, _ in groups), default=1)
        pk = np.zeros((len(groups), width))
        npk = np.zeros(len(groups), dtype=np.int64)
        zmax, conc = np.zeros(len(groups)), np.zeros(len(groups))
        for row, (group, c) in enumerate(groups):
            cps = group.delta_cp_j_mol_k or (0.0,) * len(group.pkas)
            for i, (pka, dh, cp) in enumerate(zip(group.pkas, group.delta_h_kj_mol, cps)):
                pk[row, i] = ph_state.pk_at(pka, dh, cp, temperature_k)
            npk[row], zmax[row], conc[row] = len(group.pkas), float(group.z_max), float(c)
        kargs = (pk, npk, zmax, conc, float(sid), ph_state.pkw_at(temperature_k), ph_state.davies_a(temperature_k),
                 bool(ph_state.B51.ionic_strength), float(guess), float(ionic), 1e-10, ph_state.IONIC_STRENGTH_TOL_M,
                 ph_state.IONIC_STRENGTH_SWEEPS)
        y = np.full(sulfur.N_SULFUR_STATE, 0.1)
        us = {
            "rate_constants": per_call(lambda: saved[0](*a, **kw)),
            "rate_vector": per_call(lambda: sulfur.rate_vector(k)),
            "ph_solve_wrapper_and_kernel": per_call(lambda: saved[1](*sa, **skw)),
            "ph_kernel_alone": per_call(lambda: ph_kernel.solve(*kargs)),
            "mass_action_rhs": per_call(lambda: sulfur.SULFUR_STOICHIOMETRY @ sulfur._rates_vectorised(y, kv)),
            "analytic_jacobian": per_call(lambda: sulfur.SULFUR_STOICHIOMETRY @ sulfur._rate_jacobian(y, kv), 200),
        }
    us = {key: round(v * 1e6, 2) for key, v in us.items()}
    per_rhs = us["rate_constants"] + us["rate_vector"] + 2 * us["ph_solve_wrapper_and_kernel"] + us["mass_action_rhs"]
    return {"rival": rival, "pot": label, "microseconds_per_call": us, "ph_solves_per_rhs": 2,
            "microseconds_per_rhs_call": round(per_rhs, 1),
            "shares_of_one_rhs_call": {
                "rate_constants_python": round((us["rate_constants"] + us["rate_vector"]) / per_rhs, 3),
                "ph_solves": round(2 * us["ph_solve_wrapper_and_kernel"] / per_rhs, 3),
                "of_which_numba_kernel": round(2 * us["ph_kernel_alone"] / per_rhs, 3),
                "mass_action_rhs": round(us["mass_action_rhs"] / per_rhs, 3)},
            "ph_kernel_compiled": bool(ph_kernel.COMPILED),
            "captured_calls": {key: len(v) for key, v in captured.items()}, **gate}


# ----------------------------------------------------------------------------------------------------------------
# --part memo: a bit-identical candidate, measured (NOT installed): the pure functions of temperature memoised
# ----------------------------------------------------------------------------------------------------------------


def _same(a, b) -> bool:
    if hasattr(a, "concentrations"):
        return (np.array_equal(a.concentrations, b.concentrations) and np.array_equal(a.times_min, b.times_min)
                and _same(a.metadata, b.metadata))
    if isinstance(a, np.ndarray):
        return np.array_equal(a, b, equal_nan=True)
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a if k != "seconds")   # a wall time
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and a != a:
        return isinstance(b, float) and b != b
    return a == b


@contextlib.contextmanager
def temperature_memo():
    """davies_a(T), pkw_at(T), pk_at(..., T) and the temperature part of parameters_sulfur.pka_at, memoised (the
    same floats, computed once per temperature and B51 setting). The candidate of the ranked list; restored on exit."""
    import functools

    from src.kinetic_core import parameters_sulfur as ps

    saved = (ph_state.davies_a, ph_state.pkw_at, ph_state.pk_at, ps.pka_at)
    orig_pk, orig_pka = saved[2], saved[3]
    ph_state.davies_a = functools.lru_cache(maxsize=256)(saved[0])
    ph_state.pkw_at = functools.lru_cache(maxsize=256)(saved[1])

    @functools.lru_cache(maxsize=4096)
    def _pk(pk25, dh, dcp, t, heat_capacity):
        return orig_pk(pk25, dh, dcp, t)

    @functools.lru_cache(maxsize=4096)
    def _pka_t(name, t, b51):
        return orig_pka(name, t, 0.0)          # davies_term(0, T) is exactly 0.0

    ph_state.pk_at = lambda pk25, dh, dcp, t: _pk(pk25, dh, dcp, t, ph_state.B51.heat_capacity)
    ps.pka_at = lambda name, t, i=0.0: _pka_t(name, t, ph_state.B51) - ph_state.davies_term(i, t)
    try:
        yield
    finally:
        ph_state.davies_a, ph_state.pkw_at, ph_state.pk_at, ps.pka_at = saved


def part_memo(rival: str, labels: Sequence[str], repeats: int) -> Dict[str, Any]:
    obj, theta, _ = objective_at(rival)
    os.environ["MAILLARD_WORKERS"] = "1"
    gate = load_gate("memo")
    out: Dict[str, Any] = {"pots": {}}
    for label in labels:
        run = pot_runner(obj, theta, label)
        base = run()
        t_base = min(_wall(run) for _ in range(repeats))
        with temperature_memo():
            new = run()
            t_new = min(_wall(run) for _ in range(repeats))
        out["pots"][label] = {"before_s": round(t_base, 4), "after_s": round(t_new, 4),
                              "factor": round(t_base / t_new, 3), "bit_identical": _same(base, new)}
    r0 = obj.residuals(theta)
    t0 = time.perf_counter()
    obj.residuals(theta)
    t_eval = time.perf_counter() - t0
    with temperature_memo():
        r1 = obj.residuals(theta)
        t0 = time.perf_counter()
        obj.residuals(theta)
        t_eval_new = time.perf_counter() - t0
    out["evaluation_serial"] = {"before_s": round(t_eval, 4), "after_s": round(t_eval_new, 4),
                                "factor": round(t_eval / t_eval_new, 3),
                                "residuals_bit_identical": bool(np.array_equal(r0, r1))}
    return {**out, **gate}


def _wall(fn) -> float:
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


# ----------------------------------------------------------------------------------------------------------------
# --part render
# ----------------------------------------------------------------------------------------------------------------

#: the ranked speedups (W1.3's deliverable): estimated factor on a FIT's wall time (fits are ~90 % Jacobian), the
#: evidence in this record, agent effort in hours, and whether the numbers stay bit-identical. Estimates, not
#: measurements, except where the evidence says "measured".
SPEEDUPS: Tuple[Dict[str, str], ...] = (
    {"rank": "1", "change": "Compiled RHS path: the rate constants (T-stage once per pot, pH-stage per call), both pH "
     "solves and mass action in one numba (or JAX) kernel; scipy's LSODA called without the per-step Python loop",
     "factor": "5-10x per pot (est.)", "evidence": "rate constants + pH solves are ~80 % of every pot "
     "(inside_*: timers; percall: ~200 us of ~215 us per RHS call is Python); integrator + glue ~10 %",
     "effort_h": "16-30", "bit_identical": "no (same physics to ~1e-12 relative; trajectories to tolerance)"},
    {"rank": "2", "change": "Exact gradients (forward sensitivities: JAX + diffrax, or CVODES-style staggered) in place "
     "of the forward-difference Jacobian", "factor": "3-6x on Jacobians (est.), on top of 1",
     "evidence": "a Jacobian is 18 (R0) / 20 (C1) full evaluations and ~90 % of a fit's wall; a sensitivity solve "
     "shares the step and factorisation, ~2-5 solve-equivalents for 17-19 parameters; also removes FD noise",
     "effort_h": "40-80 (plan W3)", "bit_identical": "no (a better Jacobian; the optimiser path changes)"},
    {"rank": "3", "change": "Worker count matched to the Jacobian's points (R0: 6 for 18 points; C1: 5 for 20), and the "
     "two rivals' lanes run one after the other at that count", "factor": "1.38x R0 (W=6), 1.18x C1 (W=5), measured against W=4 alone; ~1.1-1.2x against M1's two concurrent 4-worker lanes (est.)",
     "evidence": "jacobian_R0_w{2,4,6,8}: 43.5 / 30.5 / 22.1 / 22.2 s; per-point time 4.8 / 6.3 / 7.3 / 8.3 s "
     "(serial 3.9 s): the M4's 4 P + 6 E cores give ~3.5-4 serial-equivalents, so W > 6 buys nothing",
     "effort_h": "0 (a setting)", "bit_identical": "yes"},
    {"rank": "4", "change": "Memoise the pure functions of temperature (pK(T), pKw(T), Davies A(T), pka_at's "
     "temperature part)", "factor": "1.13-1.15x (measured)", "evidence": "memo: bit-identical pots and residuals",
     "effort_h": "2-3", "bit_identical": "yes"},
    {"rank": "5", "change": "One pool per evaluation: the flux budget and the route pots scheduled with B60's systems, "
     "longest first (by the previous evaluation's times)", "factor": "1.15-1.4x per single evaluation; ~1.02x per fit",
     "evidence": "scaling_*: the lone flux budget (0.31 s at rtol 1e-9) is a second barrier, 20 % of R0's wall "
     "at W=4; longest-first makespans in best_run_analysis; single evaluations are ~10 % of a fit",
     "effort_h": "3-5 (in M1's runner: B23's generator is frozen)", "bit_identical": "yes"},
    {"rank": "6", "change": "Flatten a Jacobian's points x pots into one longest-first pool", "factor": "1.0-1.08x",
     "evidence": "jacobian_*: point-level busy 92 % (R0, W=4: last round half empty) and 99 % (C1); "
     "rank 3 recovers most of it for free", "effort_h": "4-6", "bit_identical": "yes"},
    {"rank": "7", "change": "Reuse the centre: least_squares has just evaluated f(x) when it asks for jac(x)",
     "factor": "1.0-1.06x", "evidence": "1 of 18 / 20 points; saves wall only when it removes a round",
     "effort_h": "1", "bit_identical": "yes (residuals); the evaluation count in the records changes"},
    {"rank": "8", "change": "Sparse analytic Jacobian (precomputed nonzeros of S dr/dy)", "factor": "~1.03x",
     "evidence": "percall: 99 us per Jacobian vs 13 us per RHS, but LSODA asks for it ~1 in 7 RHS calls in stiff pots "
     "and never in the typical (Adams-mode) pot; at n = 62 a sparse LU does not pay", "effort_h": "2",
     "bit_identical": "no (summation order)"},
    {"rank": "9", "change": "Fewer pH solves: the Jacobian call reuses its RHS call's pH (same t, y)", "factor": "~1.05x",
     "evidence": "13 % of solves in stiff pots; not bit-identical (the warm-started bisection lands elsewhere "
     "within its 1e-10 tolerance). The label-scale (25 C) solve is physics (B51) and stays; a pH state solved as an "
     "algebraic variable belongs to the rewrite", "effort_h": "1", "bit_identical": "no"},
    {"rank": "10", "change": "Persistent worker pool (MAILLARD_POOL=persistent)", "factor": "~1.01x",
     "evidence": "fork + start latency 10-24 ms and tail 3-8 ms per pool call against 1.7-30 s; result pickles "
     "7-117 kB in < 0.3 ms", "effort_h": "0 (exists)", "bit_identical": "yes"},
)


def fit_utilisation(d: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """M1's fit records against this record's serial evaluation times."""
    out = []
    for rival in ("R0", "C1"):
        serial = (d.get(f"scaling_{rival}") or {}).get("wall_by_workers", {}).get("1", {}).get("min_s")
        if not serial:
            continue
        for start in range(3):
            rec = M1._read(M1.fit_path(rival, start))
            if not rec.get("complete"):
                continue
            ev, sec, w = int(rec["evaluations"]), float(rec["seconds"]), int(rec.get("workers") or 1)
            out.append({"fit": f"{rival} s{start}", "evaluations": ev, "seconds": sec, "workers": w,
                        "seconds_per_evaluation": round(sec / ev, 2), "utilisation": round(ev * serial / (sec * w), 3)})
    return out


def _md_table(header: Sequence[str], rows: Sequence[Sequence[Any]]) -> List[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


#: what the measurements say (W1.3's answer), filled from this record's numbers
FINDINGS: Tuple[str, ...] = (
    "One evaluation: R0 3.92 s serial -> 2.28 / 1.67 / 1.18 s at 2 / 4 / 8 workers (efficiency 0.86 / 0.59 / 0.41); "
    "C1 9.86 s -> 5.42 / 3.22 / 2.92 s (0.91 / 0.77 / 0.42; C1's other W=4 repeats 3.88 - 3.96 s at load ~3).",
    "Stragglers are NOT what limits an evaluation. R0's 38 serial tasks: median 0.068 s, the slowest "
    "(b60:hofmann_pentose_pH5, 0.46 s) 12 % of the total, the slowest 10 60 %. At W=4 the systems phase is 90 % busy "
    "and within 4 % of a longest-first schedule (1.24 vs 1.20 s; lower bound 1.19 s). C1's route phase at W=4 is 89 % "
    "busy (longest first: 1.82 -> 1.66 s); only at W=8 does order matter (1.54 -> 1.19 s).",
    "What does limit it: (a) every task runs slower when several run at once -- R0's pot-time sum is 3.59 s serial, "
    "4.77 s at W=4 (+33 %), 5.50 s at W=8 (+53 %), with CPU time = wall time, i.e. slower cores (4 performance + 6 "
    "efficiency cores, background load 2 - 3), not waiting; (b) the one flux budget (cerny_ternary, 0.31 s: it "
    "re-integrates the pot at the engine's default rtol 1e-9, not the fit's 1e-5) runs alone after the systems, a "
    "second barrier worth 20 % of R0's wall at W=4 and 28 % at W=8.",
    "Pool overhead is negligible: fork + first task 9 - 24 ms, tail 3 - 8 ms per pool call (<= 2 % of an evaluation, "
    "< 0.1 % of a Jacobian); results pickle in < 0.3 ms (3 - 117 kB). A persistent pool would buy ~1 %.",
    "The finite-difference Jacobian is the fit's cost (~90 % of its wall: R0 s0's 161 evaluations are ~8 Jacobians "
    "of 18 points + ~17 single evaluations). Its points are equal-sized (CV 4 - 14 %) and the pool is 92 % (R0, 18 "
    "points on 4 workers: the 5th round half empty) to 99 % (C1, 20 points) busy. The wall follows the per-point slow-"
    "down: R0's point takes 3.9 s alone, 4.8 / 6.3 / 7.3 / 8.3 s with 2 / 4 / 6 / 8 running; the Jacobian takes 43.5 / "
    "30.5 / 22.1 / 22.2 s. C1: 74.8 s (W=4), 63.5 s (W=5, 4 full rounds), 54.3 s (W=8). The machine delivers ~3.5 - 4 "
    "serial-equivalents however many workers; OPENBLAS_NUM_THREADS=1 changes nothing (30.46 vs 30.46 s).",
    "The '~25 % busy' of M1's fits: their records give a serial-equivalent utilisation of 37 - 57 % (evaluations x the "
    "serial evaluation time / (seconds x 4 workers)). The deficit is the per-core slowdown -- R0's and C1's lanes ran "
    "together, 8 processes, each point ~2.1x slower than alone -- plus R0's half-empty last round; it is not workers "
    "idling behind a straggler.",
    "Inside a pot, single process (timer-wrapped functions; the same in the slow and the typical pot): ~52 % is "
    "Python building the rate constants (sulfur_rate_constants_at re-derives all 132 from the parameter registry at "
    "every RHS and Jacobian call, 127 us), ~30 % the two pH solves per call (in-situ and the 25 C label pH, 30 us "
    "each: 17 us the numba kernel, 13 us the Python setup), 4 - 6 % mass action, 2 - 3 % the analytic Jacobian, ~10 % "
    "LSODA's Fortran, solve_ivp's per-step Python and the closures. hofmann_pentose_pH5: 900 steps, 1443 RHS, 218 "
    "Jacobian calls; whitfield_nf_cys (typical): 116 steps, 231 RHS, 0 Jacobians (LSODA stays in Adams mode); the "
    "slowest route pot: 940 steps, 1729 RHS, 354 Jacobians.",
    "The fixed-point coupling at R0's optimum, W=4: 14.6 s against 1.67 s inline (8.8x slower); its flux budget alone "
    "is 3.0 s.",
    "Nothing was implemented: no change that is small and bit-identical reaches 1.5x. The best one, memoising the "
    "pure functions of temperature, is bit-identical and gives 1.12 - 1.14x per pot, 1.13x (R0) / 1.17x (C1) per "
    "serial evaluation. The large factors (compiled constants + pH, exact gradients) change numbers at the rounding "
    "level and belong to the rewrite (W3).",
)


def render(d: Mapping[str, Any]) -> str:
    L: List[str] = [f"# Where one fit evaluation's time goes (plan step W1.3, 2026-11-05)", "",
                    f"Generated by `{GENERATED_BY}` from `{OUT_JSON.name}`; do not edit by hand.", ""]
    m = d.get("machine", {})
    L += [f"Machine: {m.get('cpu')} ({m.get('performance_cores')} performance + {m.get('efficiency_cores')} "
          f"efficiency cores), Python {m.get('python')}, numpy {m.get('numpy')}, scipy {m.get('scipy')}. Objective: "
          "M1's (`fit_m1_supply_combinations.Objective`) at the carried optima, inline pH coupling, rtol 1e-5 "
          "(the fits' numerics). Every timed evaluation's residual vector is bit-identical to the first and its cost "
          "equals the carried record's. Load = the 1-min load average before the measurement (limit "
          f"{LOAD_LIMIT}; a measurement above it after {LOAD_RETRIES} waits is labelled).", ""]
    # 1. scaling
    L += ["## 1. One evaluation against the worker count", ""]
    rows = []
    for rival in ("R0", "C1"):
        s = d.get(f"scaling_{rival}")
        if not s:
            continue
        for w, v in s["wall_by_workers"].items():
            rows.append([rival, w, v["min_s"], ", ".join(map(str, v["all_s"])), v["speedup_vs_1"], v["efficiency"],
                         ", ".join(map(str, v["load_1min"])), "; ".join(v["labels"]) or "clean"])
    L += _md_table(["rival", "workers", "wall min (s)", "all repeats (s)", "speed-up", "efficiency", "load",
                    "label"], rows) + [""]
    for rival in ("R0", "C1"):
        s = d.get(f"scaling_{rival}")
        if not s:
            continue
        L += [f"**{rival}** ({s['n_free']} free, {s['n_residuals']} residuals, cost {s['cost']:.4f} = carried "
              f"{s['carried_cost']:.4f}; first, cold evaluation {s['cold_first_evaluation_s']} s, excluded). "
              "Per pool call of the fastest run at each worker count:", ""]
        rows = []
        for w, a in s["best_run_analysis"].items():
            for c in a["calls"]:
                rows.append([w, c["phase"], c["n_items"], c["workers"], c["wall_s"], c["task_sum_s"],
                             f"{c['slowest_task_s']} ({c['slowest_task']})",
                             c["start_latency_s"], c["tail_latency_s"], c["mean_busy_fraction"],
                             c["makespan_pool_order_no_overhead_s"], c["makespan_longest_first_no_overhead_s"]])
        L += _md_table(["W", "phase", "tasks", "workers", "wall (s)", "task sum (s)", "slowest (s)",
                        "start lat. (s)", "tail lat. (s)", "busy", "pool order, no overhead (s)",
                        "longest first (s)"], rows) + [""]
        sp = s["serial_pots"]
        L += [f"Serial per-pot times ({sp['n_tasks']} tasks, {sp['task_sum_s']} s): quantiles "
              + ", ".join(f"p{q} {v}" for q, v in sp["quantiles_s"].items())
              + f" s; the slowest pot is {sp['slowest_1_share']:.0%} of the total, the slowest 10 are "
              f"{sp['slowest_10_share']:.0%}:", ""]
        L += _md_table(["pot", "phase", "s", "share"], [[t["label"], t["phase"], t["seconds"], t["share"]]
                                                        for t in sp["slowest_10"]]) + [""]
    # 2. inside
    L += ["## 2. Inside one pot (single process)", "",
          "Exclusive time of the timer-wrapped engine functions (time in a wrapped callee is not counted twice; the "
          "rest is LSODA's Fortran, solve_ivp's per-step Python, the closures and the post-run audit), and cProfile's "
          "call counts (cProfile inflates Python-heavy code ~2.3x, so its times are not used for shares).", ""]
    for key in ("inside_R0", "inside_C1"):
        for label, p in (d.get(key) or {}).items():
            fns = p["cprofile"]["functions"]
            calls = {k: fns.get(k, {}).get("calls", 0) for k in ("rhs", "jac", "_step_impl", "constants")}
            share = {k: v["share_of_wrapped_wall"] for k, v in p["timers_exclusive"].items()}
            rc = share.get("sulfur.sulfur_rate_constants_at", 0) + share.get("sulfur.rate_vector", 0)
            ph = sum(v for k, v in share.items() if k.startswith("ph_state."))
            rhs = share.get("sulfur._rates_vectorised", 0)
            jac = share.get("sulfur._rate_jacobian", 0)
            obs = share.get("sulfur.reaction_flux_integrals", 0) + share.get("sulfur.sulfur_reaction_rates", 0) + \
                share.get("sulfur.rate_vector_at", 0)
            rest = max(0.0, 1.0 - rc - ph - rhs - jac - obs)
            L += [f"**{label}** ({key[-2:]}): {p['wall_s_min']} s (load {p['load_1min']}); {calls['_step_impl']} "
                  f"integrator steps, {calls['rhs']} RHS calls, {calls['jac']} analytic-Jacobian calls, "
                  f"{calls['constants']} rate-constant builds; result pickle {p['result_pickle_bytes']} B in "
                  f"{p['result_pickle_s'] * 1e3:.2f} ms.", ""]
            L += _md_table(["rate constants (Python)", "pH solves (2 per call)", "mass-action RHS",
                            "analytic Jacobian", "flux audit", "integrator + Python glue"],
                           [[f"{rc:.0%}", f"{ph:.0%}", f"{rhs:.0%}", f"{jac:.0%}", f"{obs:.0%}", f"{rest:.0%}"]]) + [""]
    pc = d.get("percall_R0")
    if pc:
        us = pc["microseconds_per_call"]
        L += [f"Per call at a captured state of {pc['pot']} (us): " + ", ".join(f"{k} {v}" for k, v in us.items())
              + f". One RHS call costs ~{pc['microseconds_per_rhs_call']} us: "
              + ", ".join(f"{k} {v:.0%}" for k, v in pc["shares_of_one_rhs_call"].items()) + ".", ""]
    # 3. jacobian
    L += ["## 3. One forward-difference Jacobian (points across workers, pots serial inside each)", ""]
    rows = []
    for key in sorted(k for k in d if k.startswith("jacobian_")):
        j = d[key]
        rows.append([j["rival"], j["workers"], (j.get("env_threads") or {}).get("OPENBLAS_NUM_THREADS") or "",
                     j["n_points"], j["rounds"], j["wall_s"], j["point_mean_s"], j.get("cpu_over_wall", ""),
                     j["mean_busy_fraction"], ", ".join(f"{w['idle_s']}" for w in j["per_worker"]),
                     j["makespan_longest_first_no_overhead_s"], j["load_1min"]])
    L += _md_table(["rival", "W", "OPENBLAS threads", "points", "rounds", "wall (s)", "mean point (s)", "CPU / wall",
                    "busy", "idle per worker (s)", "longest first (s)", "load"], rows) + [""]
    # 4. fixed point
    fp = d.get("fixed_point_R0")
    if fp:
        L += ["## 4. The fixed-point coupling, R0's optimum, 4 workers", "",
              f"{fp['numerics']['name']}: {fp['wall_s']} s (inline: "
              f"{d.get('scaling_R0', {}).get('wall_by_workers', {}).get('4', {}).get('min_s')} s); task sum "
              f"{fp['task_sum_s']} s; cost {fp['cost']:.4f}; load {fp['load_1min']}. Slowest: "
              + ", ".join(f"{t['label']} {t['seconds']} s" for t in fp["slowest_5"]) + ".", ""]
    # memo
    memos = [(r, d[f"memo_{r}"]) for r in ("R0", "C1") if d.get(f"memo_{r}")]
    if memos:
        L += ["## 5. A bit-identical candidate, measured and not installed: temperature memoisation", ""]
        L += _md_table(["rival", "pot", "before (s)", "after (s)", "factor", "bit-identical"],
                       [[r, k, v["before_s"], v["after_s"], v["factor"], v["bit_identical"]]
                        for r, mm in memos for k, v in mm["pots"].items()]
                       + [[r, "whole evaluation, serial", mm["evaluation_serial"]["before_s"],
                           mm["evaluation_serial"]["after_s"], mm["evaluation_serial"]["factor"],
                           mm["evaluation_serial"]["residuals_bit_identical"]] for r, mm in memos])
        L += ["", "Below the 1.5x bar, so not implemented.", ""]
    util = d.get("m1_fit_utilisation")
    if util:
        L += ["## 6. M1's fits: serial-equivalent utilisation", "",
              "evaluations x the serial evaluation time (section 1) / (fit seconds x workers); the R0 and C1 lanes "
              "ran at the same time, 4 workers each.", ""]
        L += _md_table(["fit", "evaluations", "seconds", "workers", "s per evaluation", "utilisation"],
                       [[u["fit"], u["evaluations"], u["seconds"], u["workers"], u["seconds_per_evaluation"],
                         u["utilisation"]] for u in util]) + [""]
    # findings + ranked list
    L += ["## Findings", ""] + [f"- {f}" for f in d.get("findings", [])] + [""]
    L += ["## Ranked speedups (estimated factor on a fit's wall time unless marked measured)", ""]
    L += _md_table(["#", "change", "factor", "evidence", "agent effort (h)", "bit-identical"],
                   [[s["rank"], s["change"], s["factor"], s["evidence"], s["effort_h"], s["bit_identical"]]
                    for s in SPEEDUPS]) + [""]
    return "\n".join(L)


# ----------------------------------------------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------------------------------------------


def machine() -> Dict[str, Any]:
    import platform
    import subprocess

    def sysctl(name):
        try:
            return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True).stdout.strip() or None
        except OSError:
            return None
    return {"cpu_count": os.cpu_count(), "platform": sys.platform, "machine": platform.machine(),
            "cpu": sysctl("machdep.cpu.brand_string"), "performance_cores": sysctl("hw.perflevel0.physicalcpu"),
            "efficiency_cores": sysctl("hw.perflevel1.physicalcpu"), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": __import__("scipy").__version__}


def update_json(key: str, value: Any) -> None:
    data = json.loads(OUT_JSON.read_text(encoding="utf-8")) if OUT_JSON.exists() else {}
    data.setdefault("artifact", STEM)
    data["generated_by"] = GENERATED_BY
    data.setdefault("machine", machine())
    data[key] = value
    OUT_JSON.write_text(json.dumps(data, indent=1, sort_keys=False) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--part", required=True, choices=("scaling", "inside", "jacobian", "fixed", "percall",
                                                         "memo", "render"))
    ap.add_argument("--rival", choices=("R0", "C1"), default="R0")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pots", default="", help="--part inside: comma-separated labels (b60:<sys>, flux:<sys>, "
                                                "route:<pot>)")
    ap.add_argument("--key", default="", help="the JSON key to store under (default from the part)")
    args = ap.parse_args(argv)
    if args.part == "scaling":
        update_json(args.key or f"scaling_{args.rival}", part_scaling(args.rival, args.repeats))
    elif args.part == "inside":
        update_json(args.key or f"inside_{args.rival}", part_inside(args.rival, args.pots.split(","), args.repeats))
    elif args.part == "jacobian":
        update_json(args.key or f"jacobian_{args.rival}_w{args.workers}", part_jacobian(args.rival, args.workers))
    elif args.part == "fixed":
        update_json(args.key or "fixed_point_R0", part_fixed(args.workers))
    elif args.part == "percall":
        update_json(args.key or f"percall_{args.rival}", part_percall(args.rival, args.pots or "b60:hofmann_pentose_pH5"))
    elif args.part == "memo":
        update_json(args.key or f"memo_{args.rival}", part_memo(args.rival, args.pots.split(","), args.repeats))
    elif args.part == "render":
        data = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        data["m1_fit_utilisation"] = fit_utilisation(data)
        data["findings"] = list(FINDINGS)
        data["ranked_speedups"] = list(SPEEDUPS)
        OUT_JSON.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
        OUT_MD.write_text(render(data), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

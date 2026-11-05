#!/usr/bin/env python
"""
build_jax_bakeoff_reference.py -- today's engine's numbers for the JAX bake-off (STATUS_AND_PLAN W4 step H2, pulled
forward; 2026-11-04). Runs in the CORE env (``./scripts/docker_maillard.sh run``), never in maillard-jax: every number
the JAX port is judged against comes from the package set that backs results/validation/.

One representative pot, Schieberle 2000 Table IV (ribose 100 + cysteine 33 mmol/L, 0.5 M phosphate, label pH 5.0,
100 C, 12 h; thiol_submodel), on the R3En joint fit's engine state and optimum (B60 + B58/B59b/B59a/B61, inline pH
coupling, trunk solutes in the balance; kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.json). Writes, for
scripts/generators/build_jax_bakeoff.py to read:

* the trajectory at the measured times and a dense grid, at the probe's tolerances (rtol 1e-8, atol 1e-14) and tight
  (rtol 1e-10, atol 1e-16), with the in-situ and the cooled (25 C label) pH;
* the rate constants ``sulfur_rate_constants_at`` returns on a grid of (pH, label pH, I) at two parameter points, and
  the charge balance's (pH, I) on states of the trajectory: the port's unit checks;
* the finite-difference Jacobian of the observables (log10 MFT / FFT / Cys / NF at the measured times, the final
  cooled pH) in the 16 free coordinates, as a fit forms it (forward, B8's absolute steps) and central, with its
  wall-clock;
* four stiff corners (every rate offset at its upper bound, label pH 8, cysteine x10, unbuffered).

Timings are taken on a shared machine: each carries ``uptime``'s load averages and is PROVISIONAL when they are high.

    python scripts/generators/build_jax_bakeoff_reference.py [--quick]
Writes results/validation/jax_bakeoff_sulfur_lane_reference.json.
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

from src import data_paths  # noqa: E402
from src.kinetic_core import config, engine, fit_rows, ph_state, sulfur  # noqa: E402
from src.kinetic_core import fitting  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402

FIT = "kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en"
ARM = "R3En"
TABLE, POT = "thiol_submodel", "schieberle2000_rib_cys_ph5_100C"
OBS_SPECIES = ("MFT", "FFT", "Cys", "NF")
OUT = data_paths.VALIDATION_DIR / "jax_bakeoff_sulfur_lane_reference.json"
CELSIUS = 273.15
#: the four stiff corners: name -> (job overrides, theta rule)
CORNERS = {
    "offsets_at_upper_bounds": ({}, "upper"),
    "label_ph_8": ({"label_ph": 8.0}, None),
    "cysteine_x10": ({"charge": {"PENT": 100.0, "Cys": 330.0}}, None),
    "unbuffered": ({"buffer": "none (water)"}, None),
}


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def arm_setup():
    arms = joint.ARMS[ARM]
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    joint.TDP_EXIT_STATIC = joint.tdp_exit_static(arms)
    joint.EA_STATIC = {}
    joint.STRECKER_SHARE = 0.042
    config.update({"ph.coupling": "inline"})
    return joint.coordinates(False, arms)


def base_job(times) -> Dict[str, Any]:
    pots, rows = fit_rows.load(TABLE)
    job = next(j for j in probe.plan(pots, rows) if j["pot"] == POT)
    measured = sorted({r.t_min for r in rows if r.pot == POT and r.t_min > 0})
    return dict(job, times=sorted({0.0, *measured, *times})), measured


def job_at(job, theta, keys, overrides=None) -> Dict[str, Any]:
    out = dict(job, **joint.trial(theta, keys), strecker_share=joint.STRECKER_SHARE)
    out.update(overrides or {})
    return out


def parameters_of(job):
    import dataclasses

    parameters = engine.core_parameters(engine.SULFUR)
    for key, k_ref in (job.get("k_ref") or {}).items():
        parameters[key] = dataclasses.replace(parameters[key], k_ref=float(k_ref))
    return parameters


def integrate(job, rtol=1e-8, atol=1e-14):
    """``probe.integrate_pot`` with the tolerances as arguments (no pH-stat pot here, so no clamp branch)."""
    spec = probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"])
    assert not spec.is_clamped
    state0 = dict(job["charge"])
    state0.setdefault("OX", engine.OX_AMBIENT_MMOL_L)
    state0.setdefault("OXR", engine.OX_RESERVOIR_DEFAULT_UNITS)
    return sulfur.integrate_sulfur(parameters_of(job), job["t_c"] + CELSIUS, state0, np.array(job["times"]),
                                   ph=job["label_ph"], buffer_spec=spec, ph_drift=probe._drift(job), rtol=rtol,
                                   atol=atol)


def run_record(job, rtol=1e-8, atol=1e-14, limit_s: int = 900) -> Dict[str, Any]:
    import signal

    def _alarm(signum, frame):
        raise probe.PotTimeout(f"over {limit_s} s")

    with config.applied(joint.job_config(job)):
        t0 = time.perf_counter()
        previous = signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(int(limit_s))
        try:
            run = integrate(job, rtol, atol)
        except Exception as exc:  # noqa: BLE001 - a corner may fail; that is a result
            return {"failed": f"{type(exc).__name__}: {exc}", "seconds": time.perf_counter() - t0}
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, previous)
        seconds = time.perf_counter() - t0
        drift = probe._drift(job)
        spec = probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"])
        sid = run.metadata["strong_ion_difference_mmol_L"]
        states = run.concentrations
        cooled = [ph_state.ph_of_state({k: float(row[sulfur.SULFUR_INDEX[k]]) for k in sulfur.SULFUR_STATE_KEYS},
                                       ph_state.LABEL_TEMPERATURE_K, sid, drift, spec) for row in states]
    return {"seconds": seconds, "rtol": rtol, "atol": atol, "times": list(map(float, job["times"])),
            "states": states.tolist(), "ph_in_situ": [float(v) for v in run.ph_series],
            "ph_cooled": [float(v) for v in cooled], "sid_mmol_L": float(sid),
            "worst_raw_negative": run.metadata["worst_raw_negative_before_clip"]}


def observables(rec, measured) -> np.ndarray:
    times = rec["times"]
    states = np.asarray(rec["states"])
    out = []
    for species in OBS_SPECIES:
        i = sulfur.SULFUR_INDEX[species]
        out += [math.log10(max(states[times.index(t), i], 1e-300)) for t in measured]
    out.append(rec["ph_cooled"][times.index(measured[-1])])
    return np.array(out)


def observable_names(measured) -> List[str]:
    return [f"log10 {s} @ {t:g} min" for s in OBS_SPECIES for t in measured] + [f"pH cooled @ {measured[-1]:g} min"]


def rate_constant_probe(job) -> Dict[str, Any]:
    grid = [(ph, phl, i) for ph in (3.6, 4.6, 5.0, 5.4, 6.3, 7.6, 9.2) for phl in (4.4, 4.99, 5.01, 6.4, 6.9, 8.3)
            for i in (0.0, 0.04, 0.95)]
    with config.applied(joint.job_config(job)):
        params = parameters_of(job)
        t_k = job["t_c"] + CELSIUS
        rows = [sulfur.rate_vector(sulfur.sulfur_rate_constants_at(params, t_k, ph, ph_label=phl,
                                                                   ionic_strength_m=i)).tolist()
                for ph, phl, i in grid]
    return {"grid_ph_phlabel_I": grid, "k": rows}


def charge_balance_probe(job, rec) -> Dict[str, Any]:
    drift = probe._drift(job)
    spec = probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"])
    out = []
    with config.applied(joint.job_config(job)):
        sid = rec["sid_mmol_L"]
        for row in rec["states"][::4]:
            state = {k: float(row[sulfur.SULFUR_INDEX[k]]) for k in sulfur.SULFUR_STATE_KEYS}
            hot = ph_state.ph_and_ionic_strength_of_state(state, job["t_c"] + CELSIUS, sid, drift, spec)
            cool = ph_state.ph_and_ionic_strength_of_state(state, ph_state.LABEL_TEMPERATURE_K, sid, drift, spec)
            out.append({"state": row, "hot": list(hot), "cool": list(cool)})
        initial = {k: float(rec["states"][0][sulfur.SULFUR_INDEX[k]]) for k in sulfur.SULFUR_STATE_KEYS}
        sid0 = ph_state.strong_ion_difference(job["label_ph"], initial, drift, spec)
    return {"sid_mmol_L": sid0, "points": out}


def fd_jacobian(job, theta, keys, kinds, lower, upper, measured, central: bool) -> Dict[str, Any]:
    t0 = time.perf_counter()
    f0 = observables(run_record(job_at(job, theta, keys)), measured)
    cols, runs = [], 1
    for j, kind in enumerate(kinds):
        h = fitting._step(kind)
        if central:
            up, dn = theta.copy(), theta.copy()
            up[j] += h
            dn[j] -= h
            fu = observables(run_record(job_at(job, up, keys)), measured)
            fd = observables(run_record(job_at(job, dn, keys)), measured)
            cols.append((fu - fd) / (2 * h))
            runs += 2
        else:   # fitting._forward_jacobian: inward at an upper bound
            h = h if theta[j] + h <= upper[j] else -h
            up = theta.copy()
            up[j] += h
            cols.append((observables(run_record(job_at(job, up, keys)), measured) - f0) / h)
            runs += 1
    return {"jacobian": np.array(cols).T.tolist(), "seconds": time.perf_counter() - t0, "engine_runs": runs,
            "uptime": uptime(), "steps": {k: fitting._step(kd) for k, kd in zip(keys, kinds)}}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--quick", action="store_true", help="the trajectory and the probes only (no Jacobians, corners)")
    parser.add_argument("--timing-only", action="store_true",
                        help="append one paired-timing record (5 single runs, the fit's forward-FD Jacobian) to the "
                             "existing reference file, for timings taken back to back with the JAX ones")
    args = parser.parse_args(argv)
    keys, kinds, lower, upper, _start = arm_setup()
    fit = json.loads((data_paths.VALIDATION_DIR / f"{FIT}.json").read_text())
    theta = np.array([fit["fit"]["theta"][k] for k in keys], dtype=float)
    if args.timing_only:
        report = json.loads(OUT.read_text())
        job0, measured = base_job([])
        job = job_at(job0, theta, keys)
        start = uptime()
        secs = [run_record(job)["seconds"] for _ in range(5)]
        fd = fd_jacobian(job0, theta, keys, kinds, lower, upper, measured, central=False)
        report.setdefault("paired_timing", []).append({
            "when": time.strftime("%Y-%m-%d %H:%M:%S"), "uptime_before": start, "single_run_s": secs,
            "fd_forward_jacobian_s": fd["seconds"], "fd_engine_runs": fd["engine_runs"], "uptime_after": uptime()})
        OUT.write_text(json.dumps(report, default=float) + "\n", encoding="utf-8")
        print(report["paired_timing"][-1])
        return 0
    dense = np.unique(np.concatenate([np.geomspace(1e-3, 720.0, 41), np.linspace(0.0, 720.0, 25)])).tolist()
    job0, measured = base_job(dense)
    job = job_at(job0, theta, keys)
    report: Dict[str, Any] = {
        "artifact": OUT.name, "generated_by": "scripts/generators/build_jax_bakeoff_reference.py", "env": "maillard (core)",
        "numpy": np.__version__, "fit": FIT, "arm": ARM, "pot": f"{TABLE}:{POT}", "keys": list(keys),
        "kinds": list(kinds), "lower": lower.tolist(), "upper": upper.tolist(), "theta": theta.tolist(),
        "measured_times": measured, "observables": observable_names(measured), "species": list(sulfur.SULFUR_STATE_KEYS),
        "reaction_keys": list(sulfur.FULL_REACTION_KEYS), "job": {k: v for k, v in job.items() if k != "times"},
        "uptime_start": uptime()}
    report["trajectory_probe_tol"] = run_record(job)
    report["trajectory_tight_tol"] = run_record(job, 1e-10, 1e-16)
    report["rate_constant_probe"] = {"theta_star": rate_constant_probe(job)}
    rng = np.random.default_rng(20261104)
    theta_b = np.clip(theta + rng.normal(0.0, 0.3, size=theta.size), lower, upper)
    report["rate_constant_probe"]["theta_perturbed"] = dict(rate_constant_probe(job_at(job0, theta_b, keys)),
                                                            theta=theta_b.tolist())
    report["charge_balance_probe"] = charge_balance_probe(job, report["trajectory_probe_tol"])
    if not args.quick:
        # one engine run, repeated, for the per-run cost at this load
        secs = [run_record(job)["seconds"] for _ in range(3)]
        report["single_run_seconds"] = {"runs": secs, "uptime": uptime()}
        report["fd_forward"] = fd_jacobian(job0, theta, keys, kinds, lower, upper, measured, central=False)
        report["fd_central"] = fd_jacobian(job0, theta, keys, kinds, lower, upper, measured, central=True)
        corners = {}
        for name, (over, rule) in CORNERS.items():
            th = theta.copy()
            if rule == "upper":
                th = np.where(np.array(kinds) == "log10_offset", upper, theta)
            cj = job_at(job0, th, keys, over)
            rec = run_record(cj)
            corners[name] = {"overrides": over, "theta": th.tolist(), "uptime": uptime(), **rec}
        report["corners"] = corners
    report["uptime_end"] = uptime()
    OUT.write_text(json.dumps(report, default=float) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

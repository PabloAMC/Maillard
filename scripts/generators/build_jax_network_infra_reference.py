#!/usr/bin/env python
"""
build_jax_network_infra_reference.py -- the engine's side of the A2 network infrastructure checks (plan step A2 + H3,
docs/guides/A2_MODULE_CONTRACT.md v2.2; 2026-11-05): the JAX fit driver's equivalence with ``fitting.fit``, the
forced-input driver's check against an engine run with the same pH / temperature staircase, and the engine's
draws/s for the predictive-density bridge. Runs in the CORE env (``./scripts/docker_maillard.sh run``), never in
maillard-jax; ``build_jax_network_infra.py`` (maillard-jax) reads it.

THE ENGINE IS MOVED ONLY THROUGH ITS OWN ARGUMENTS AND THE CONFIG (tasks/lessons.md, 2026-11-05): a trunk-lane rate
constant through the parameter dict ``integrate`` takes (``dataclasses.replace`` of a ``KineticParameter``), the
Amadori pH exponent through ``trunk_conditions.apply(ph_exponent=)``, the lipid breakage law through
``config.applied`` (``build_jax_lipid_module_reference.engine_at``) and the hydrogen-donor share through
``integrate_lipid(donor=)``. A staircase pH is a per-segment ``ProcessSpec`` (its declared pH), which is exactly the
engine's per-segment path (``trunk_conditions.apply`` then ``integrate``).

THE TOY PROBLEMS ARE SYNTHETIC (software equivalence only; no evidence about chemistry, no promotion):

* trunk lane (55 steps): three pots (glucose + glycine at 120 C / pH 6.0 and 140 C / pH 7.5 a_w 0.90; glyoxal +
  methylglyoxal + glycine at 120 C / pH 8.0), log10 levels of 3-DG, HMF, MGO, Fru and the pyrazines at 10 / 30 /
  100 % of the hold; five coordinates (log10 k and Ea of k_ama_tdg, log10 k of k_ddg_hmf and k_mgo_ak, the Amadori pH
  exponent); data = the engine at a declared truth (LSODA rtol 1e-12) plus a seeded N(0, 0.02 dex) draw, sigma
  0.05 dex; start = shipped.
* lipid lane: pea protein isolate at 140 / 160 / 180 C, log10 ug/L of hexanal, 2-pentylfuran, nonanal and
  2,4-decadienal; four coordinates (the breakage multiplier and Ea, the hexanal share, the donor suppression).

``fitting.fit`` runs each with the M1.1 central polish (2 calls) and ``require_converged="raise"``, serially
(``MAILLARD_WORKERS`` unset), at LSODA rtol 1e-8 / atol 1e-14 (trunk) and the analytic lipid lane.

    ./scripts/docker_maillard.sh run "python scripts/generators/build_jax_network_infra_reference.py"
Writes results/validation/jax_network_infra_reference.json.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

import build_jax_lipid_module_reference as lipref  # noqa: E402
import build_jax_trunk_lane_reference as refgen  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import config, engine, fitting, trunk_conditions  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.species import SPECIES_KEYS  # noqa: E402
from src.kinetic_core import species_lipid as sl  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "jax_network_infra_reference.json"
CELSIUS = 273.15
SEED = 20261105
NOISE_DEX = 0.02
SIGMA_DEX = 0.05
SAVE_FRACTIONS = (0.1, 0.3, 1.0)
FIT_TOL = {"rtol": 1e-8, "atol": 1e-14}
TIGHT = {"rtol": 1e-12, "atol": 1e-20}

# --------------------------------------------------------------------------- the trunk-lane toy
TRUNK_POTS = (
    {"id": "toy_glucose_glycine_120C_pH6", "initial": {"Glc": 100.0, "Gly": 50.0}, "hold": (60.0, 120.0),
     "ph": 6.0, "aw": None, "observed": ("TDG", "HMF", "MGO", "Fru")},
    {"id": "toy_glucose_glycine_140C_pH7.5_aw0.90", "initial": {"Glc": 100.0, "Gly": 50.0}, "hold": (30.0, 140.0),
     "ph": 7.5, "aw": 0.90, "observed": ("TDG", "HMF", "MGO", "Fru")},
    {"id": "toy_glyoxal_methylglyoxal_glycine_120C_pH8", "initial": {"GO": 20.0, "MGO": 20.0, "Gly": 40.0},
     "hold": (60.0, 120.0), "ph": 8.0, "aw": None, "observed": ("PZ", "MPZ", "DMP", "MGO")},
)
#: (module, coordinate, engine kind, JAX kind, lower, upper) -- bounds are offsets from the shipped value
TRUNK_COORDINATES = (
    ("trunk", "log10_k:k_ama_tdg", "log10k", "log10k", -1.5, 1.5),
    ("trunk", "ea:k_ama_tdg", "Ea", "Ea", -40.0, 40.0),
    ("trunk", "log10_k:k_ddg_hmf", "log10k", "log10k", -1.5, 1.5),
    ("trunk_lane", "log10_k:k_mgo_ak", "log10k", "log10k", -1.5, 1.5),
    ("trunk", "ph_exponent", "ph_slope", "ph_slope", None, None),       # bounds: trunk_conditions.PH_EXPONENT_BAND
)
TRUNK_TRUTH_SHIFT = (0.25, -8.0, -0.2, 0.3, None)
TRUNK_TRUTH_PH_EXPONENT = 0.80

# --------------------------------------------------------------------------- the lipid toy
LIPID_POTS = (
    {"id": "toy_pea_isolate_140C", "carriers": {"pea_protein_isolate": {}}, "hold": (30.0, 140.0)},
    {"id": "toy_pea_isolate_160C", "carriers": {"pea_protein_isolate": {}}, "hold": (20.0, 160.0)},
    {"id": "toy_pea_isolate_180C", "carriers": {"pea_protein_isolate": {}}, "hold": (15.0, 180.0)},
)
LIPID_OBSERVED = ("HEXANAL", "PENTYLFURAN", "NONANAL", "DECADIENAL")
#: (coordinate, engine kind, JAX kind, lower, upper) -- absolute bounds
LIPID_COORDINATES = (
    ("log10_k_decadienal", "log10k", "log10k", -1.5, 1.5),
    ("ea_decadienal_kj_mol", "Ea", "Ea", 5.0, 80.0),
    ("share_hexanal_from_decadienal", "acid_yield", "share", 0.0, 0.6),
    ("donor", "acid_yield", "share", 0.0, 0.9),
)
LIPID_TRUTH = {"log10_k_decadienal": 0.3, "ea_decadienal_kj_mol": 40.0, "share_hexanal_from_decadienal": 0.2,
               "donor": 0.3}

# --------------------------------------------------------------------------- the forced staircase
STAIRCASE = {"id": "forced_staircase_glucose_glycine_methylglyoxal", "initial": {"Glc": 100.0, "Gly": 50.0,
                                                                                 "MGO": 5.0},
             "segments": ((10.0, 110.0), (10.0, 120.0), (10.0, 130.0), (10.0, 130.0), (10.0, 120.0)),
             "ph": (7.6, 7.1, 6.5, 5.9, 5.4), "aw": 0.90}
DRAW_COUNTS = (64,)


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def _process(ph, aw, minutes=1.0, t_c=100.0):
    return refgen.synthetic_process((minutes, t_c), ph, aw)


# --------------------------------------------------------------------------- trunk-lane engine side
def trunk_shipped():
    cfg = refgen.config_for({}, engine_side=True)
    with config.applied(cfg):
        params = engine.core_parameters(engine.TRUNK)
    return cfg, params


def trunk_start(params) -> np.ndarray:
    out = []
    for module, coord, *_ in TRUNK_COORDINATES:
        kind, _, key = coord.partition(":")
        if kind == "log10_k":
            out.append(math.log10(float(params[key].k_ref)))
        elif kind == "ea":
            out.append(float(params[key].ea_kj_mol))
        else:
            out.append(float(trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT))
    return np.asarray(out)


def trunk_bounds(start: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    lo, hi = [], []
    for x, (_, coord, _, _, a, b) in zip(start, TRUNK_COORDINATES):
        if coord == "ph_exponent":
            lo.append(float(trunk_conditions.PH_EXPONENT_BAND[0]))
            hi.append(float(trunk_conditions.PH_EXPONENT_BAND[1]))
        else:
            lo.append(x + a)
            hi.append(x + b)
    return np.asarray(lo), np.asarray(hi)


def trunk_parameters_at(params, x: Sequence[float]):
    """The parameter dict and the condition terms at coordinate vector x (the engine's own arguments)."""
    out = dict(params)
    cond = refgen.full_conditions({})
    for value, (_, coord, *_rest) in zip(x, TRUNK_COORDINATES):
        kind, _, key = coord.partition(":")
        if kind == "log10_k":
            out[key] = dataclasses.replace(out[key], k_ref=10.0 ** float(value))
        elif kind == "ea":
            out[key] = dataclasses.replace(out[key], ea_kj_mol=float(value))
        else:
            cond["ph_exponent"] = float(value)
    return out, cond


def trunk_pot_record(pot) -> Dict[str, Any]:
    minutes, t_c = pot["hold"]
    refined, save_times = refgen.refine([[minutes, t_c]], SAVE_FRACTIONS)
    process = refgen.synthetic_process((minutes, t_c), pot["ph"], pot["aw"])
    return {"id": pot["id"], "initial": dict(pot["initial"]), "engine_initial": refgen.engine_initial(pot["initial"],
            process), "refined_segments": refined, "save_times_min": save_times, "ph": pot["ph"],
            "water_activity": pot["aw"], "observed": list(pot["observed"]), "process": process}


def trunk_states(rec, params, cond, cfg, tol) -> np.ndarray:
    return np.asarray(refgen.engine_states(params, rec["initial"], rec["process"], rec["refined_segments"], cond, cfg,
                                           tol=tol))


def trunk_observe(rec, states: np.ndarray) -> List[float]:
    out = []
    for row in range(len(rec["refined_segments"])):
        for s in rec["observed"]:
            out.append(math.log10(max(float(states[row, SPECIES_KEYS.index(s)]), 1e-300)))
    return out


def trunk_toy() -> Dict[str, Any]:
    cfg, params = trunk_shipped()
    start = trunk_start(params)
    lower, upper = trunk_bounds(start)
    truth = np.asarray([start[i] + d if d is not None else TRUNK_TRUTH_PH_EXPONENT
                        for i, d in enumerate(TRUNK_TRUTH_SHIFT)])
    recs = [trunk_pot_record(p) for p in TRUNK_POTS]
    p_truth, c_truth = trunk_parameters_at(params, truth)
    rng = np.random.default_rng(SEED)
    data, sigma, rows = [], [], []
    for pi, rec in enumerate(recs):
        clean = trunk_observe(rec, trunk_states(rec, p_truth, c_truth, cfg, TIGHT))
        noisy = np.asarray(clean) + rng.normal(0.0, NOISE_DEX, size=len(clean))
        k = 0
        for row in range(len(rec["refined_segments"])):
            for s in rec["observed"]:
                rows.append({"pot": pi, "row": row, "species": s, "observable": f"mmol_per_L:{s}",
                             "value": float(noisy[k]), "sigma": SIGMA_DEX, "clean": float(clean[k])})
                k += 1
        data += noisy.tolist()
        sigma += [SIGMA_DEX] * len(clean)
    y, sig = np.asarray(data), np.asarray(sigma)
    evals = {"n": 0, "seconds": 0.0}

    def residuals(x, _rows):
        t0 = time.perf_counter()
        p, c = trunk_parameters_at(params, x)
        pred = []
        for rec in recs:
            pred += trunk_observe(rec, trunk_states(rec, p, c, cfg, FIT_TOL))
        evals["n"] += 1
        evals["seconds"] += time.perf_counter() - t0
        return (np.asarray(pred) - y) / sig

    problem = fitting.Problem(keys=tuple(f"{m}.{c}" for m, c, *_ in TRUNK_COORDINATES),
                              kinds=tuple(k for _, _, k, *_ in TRUNK_COORDINATES), lower=lower, upper=upper,
                              rows=tuple(range(len(y))), residuals=residuals)
    out = _engine_fits(problem, start, evals)
    out.update({"pots": [{k: v for k, v in rec.items() if k != "process"} for rec in recs], "rows": rows,
                "coordinates": [{"module": m, "coordinate": c, "engine_kind": ek, "jax_kind": jk}
                                for m, c, ek, jk, *_ in TRUNK_COORDINATES],
                "lower": lower.tolist(), "upper": upper.tolist(), "start": start.tolist(), "truth": truth.tolist(),
                "fit_tolerance": FIT_TOL, "data_tolerance": TIGHT, "noise_dex": NOISE_DEX, "sigma_dex": SIGMA_DEX,
                "seed": SEED, "config_digest": cfg.digest()})
    return out


def _engine_fits(problem, start, evals) -> Dict[str, Any]:
    out: Dict[str, Any] = {"workers": int(os.environ.get("MAILLARD_WORKERS", "1"))}
    evals["n"], evals["seconds"] = 0, 0.0
    t0 = time.perf_counter()
    forward = fitting.fit(problem, start, budget=1200)
    out["forward_only"] = {**forward.as_dict(problem.keys), "wall_seconds": time.perf_counter() - t0,
                           "residual_evaluations": evals["n"], "seconds_in_residuals": evals["seconds"]}
    evals["n"], evals["seconds"] = 0, 0.0
    t0 = time.perf_counter()
    polished = fitting.fit(problem, start, budget=1200, central_polish=2, require_converged="raise")
    out["polished"] = {**polished.as_dict(problem.keys), "wall_seconds": time.perf_counter() - t0,
                       "residual_evaluations": evals["n"], "seconds_in_residuals": evals["seconds"]}
    out["seconds_per_residual_evaluation"] = evals["seconds"] / max(1, evals["n"])
    out["polished_residuals"] = np.asarray(polished.residuals).tolist()
    return out


# --------------------------------------------------------------------------- lipid engine side
def lipid_toy() -> Dict[str, Any]:
    shipped = lipref.theta_shipped()
    names = [c for c, *_ in LIPID_COORDINATES]
    start = np.asarray([shipped[c] for c in names])
    lower = np.asarray([lo for *_, lo, _ in LIPID_COORDINATES])
    upper = np.asarray([hi for *_, hi in LIPID_COORDINATES])
    recs = []
    for pot in LIPID_POTS:
        minutes, t_c = pot["hold"]
        refined, save_times = lipref.refine([(minutes, t_c)], SAVE_FRACTIONS)
        recs.append({"id": pot["id"], "carriers": pot["carriers"], "refined_segments": refined,
                     "save_times_min": save_times, "observed": list(LIPID_OBSERVED)})
    mw = dict(sl.MOLECULAR_WEIGHT_G_PER_MOL)

    def observe(rec, theta):
        states = np.asarray(lipref.engine_states(rec["carriers"], rec["refined_segments"], theta))
        return [math.log10(max(float(states[row, sl.LIPID_KEYS.index(s)]) * mw[s] * 1.0e3, 1e-300))
                for row in range(len(rec["refined_segments"])) for s in rec["observed"]]

    def theta_at(x):
        return dict(shipped, **{c: float(v) for c, v in zip(names, x)})

    truth = np.asarray([LIPID_TRUTH[c] for c in names])
    rng = np.random.default_rng(SEED + 1)
    rows, data = [], []
    for pi, rec in enumerate(recs):
        clean = observe(rec, theta_at(truth))
        noisy = np.asarray(clean) + rng.normal(0.0, NOISE_DEX, size=len(clean))
        k = 0
        for row in range(len(rec["refined_segments"])):
            for s in rec["observed"]:
                rows.append({"pot": pi, "row": row, "species": s, "observable": f"ug_per_L:{s}",
                             "value": float(noisy[k]), "sigma": SIGMA_DEX, "clean": float(clean[k])})
                k += 1
        data += noisy.tolist()
    y = np.asarray(data)
    evals = {"n": 0, "seconds": 0.0}

    def residuals(x, _rows):
        t0 = time.perf_counter()
        th = theta_at(x)
        pred = [v for rec in recs for v in observe(rec, th)]
        evals["n"] += 1
        evals["seconds"] += time.perf_counter() - t0
        return (np.asarray(pred) - y) / SIGMA_DEX

    problem = fitting.Problem(keys=tuple(f"lipid.{c}" for c in names),
                              kinds=tuple(k for _, k, *_ in LIPID_COORDINATES), lower=lower, upper=upper,
                              rows=tuple(range(len(y))), residuals=residuals)
    out = _engine_fits(problem, start, evals)
    out.update({"pots": recs, "rows": rows,
                "coordinates": [{"module": "lipid", "coordinate": c, "engine_kind": ek, "jax_kind": jk}
                                for c, ek, jk, *_ in LIPID_COORDINATES],
                "lower": lower.tolist(), "upper": upper.tolist(), "start": start.tolist(), "truth": truth.tolist(),
                "fixed": {k: v for k, v in shipped.items() if k not in names}, "noise_dex": NOISE_DEX,
                "sigma_dex": SIGMA_DEX, "seed": SEED + 1})
    return out


# --------------------------------------------------------------------------- forced staircase
def forced_staircase() -> Dict[str, Any]:
    """The engine's per-segment path with the pH and temperature changing per segment: what a 'previous' forced
    course must reproduce."""
    cfg, params = trunk_shipped()
    cond = refgen.full_conditions({})
    first = _process(STAIRCASE["ph"][0], STAIRCASE["aw"])
    state = refgen.engine_initial(STAIRCASE["initial"], first)
    initial = dict(state)
    rows = []
    with config.applied(cfg):
        for (minutes, t_c), ph in zip(STAIRCASE["segments"], STAIRCASE["ph"]):
            process = _process(ph, STAIRCASE["aw"], minutes, t_c)
            p, _ = trunk_conditions.apply(params, process, **refgen.apply_kwargs(cond))
            run = integrate(p, float(t_c) + CELSIUS, state, np.array([0.0, float(minutes)]), **TIGHT)
            state = {k: float(run.concentrations[-1, i]) for i, k in enumerate(SPECIES_KEYS)}
            rows.append([state[k] for k in SPECIES_KEYS])
    return {"id": STAIRCASE["id"], "initial": STAIRCASE["initial"], "engine_initial": initial,
            "segments": [list(s) for s in STAIRCASE["segments"]], "ph": list(STAIRCASE["ph"]),
            "water_activity": STAIRCASE["aw"], "states": rows, "engine_species": list(SPECIES_KEYS),
            "tolerance": TIGHT}


# --------------------------------------------------------------------------- engine draws/s
def engine_draw_timing(toy: Dict[str, Any]) -> Dict[str, Any]:
    """Draws/s of the engine on the first trunk toy pot (60 min at 120 C, three refined segments): one LSODA run per
    segment per draw at fit tolerance, serially, the draws a Gaussian around the shipped start (0.1 dex, 3 kJ/mol,
    0.05 in the exponent)."""
    cfg, params = trunk_shipped()
    rec = trunk_pot_record(TRUNK_POTS[0])
    start = np.asarray(toy["start"])
    scale = np.asarray([0.1, 3.0, 0.1, 0.1, 0.05])
    rng = np.random.default_rng(SEED + 2)
    out = {"pot": rec["id"], "tolerance": FIT_TOL, "load_before": uptime(), "rows": []}
    for n in DRAW_COUNTS:
        draws = start + rng.normal(size=(n, len(start))) * scale
        t0 = time.perf_counter()
        for x in draws:
            p, c = trunk_parameters_at(params, x)
            trunk_states(rec, p, c, cfg, FIT_TOL)
        dt = time.perf_counter() - t0
        out["rows"].append({"draws": n, "seconds": dt, "draws_per_second": n / dt})
    out["load_after"] = uptime()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", action="append", choices=("trunk", "lipid", "forced", "draws"))
    args = ap.parse_args(argv)
    stages = args.stage or ["trunk", "lipid", "forced", "draws"]
    report: Dict[str, Any] = json.loads(OUT.read_text()) if OUT.exists() and args.stage else {}
    report.update({"generator": "scripts/generators/build_jax_network_infra_reference.py", "env": "maillard (core)",
                   "numpy": np.__version__})
    t0 = time.perf_counter()
    if "trunk" in stages:
        report["trunk_toy"] = trunk_toy()
        print("trunk toy:", report["trunk_toy"]["polished"]["cost"], report["trunk_toy"]["polished"]["wall_seconds"])
    if "lipid" in stages:
        report["lipid_toy"] = lipid_toy()
        print("lipid toy:", report["lipid_toy"]["polished"]["cost"], report["lipid_toy"]["polished"]["wall_seconds"])
    if "forced" in stages:
        report["forced_staircase"] = forced_staircase()
    if "draws" in stages:
        report["engine_draws"] = engine_draw_timing(report["trunk_toy"])
        print("engine draws/s:", [r["draws_per_second"] for r in report["engine_draws"]["rows"]])
    report["seconds"] = time.perf_counter() - t0
    OUT.write_text(json.dumps(report, indent=None, sort_keys=True) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

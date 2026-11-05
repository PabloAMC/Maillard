#!/usr/bin/env python
"""
build_jax_lipid_module_reference.py -- today's engine's numbers for the JAX lipid module (plan step A2 + H3, the
first module under docs/guides/A2_MODULE_CONTRACT.md v2; 2026-11-05). Runs in the CORE env
(``./scripts/docker_maillard.sh run``), never in maillard-jax: every number the port is judged against comes from the
package set that backs results/validation/.

The pots: every panel bundle the engine routes to the lipid lane (``panel.panel_bundles`` -> ``engine.declare_envelope``
with the bundle's own targets; 8 on 2026-11-05: three trust-loop isolates, the Trikusuma UHT pot, the acrylamide SPI
extrusion pot's lipid half and four external hexanal pots), each at its point and at the engine's two declared-band
corners (Q10, lipid fraction and peroxide value all low / all high, ``engine._run_lipid_lane``), plus one SYNTHETIC
staircase profile (pea + soy isolate, 25 -> 160 C in 12 min, 20 min hold, cooled to 60 C; not a benchmark) that
crosses the decadienal window twice. For each, the state after every segment of the program refined at its save
times (``lipid.integrate_lipid`` on each prefix, summed over carriers as the engine does); the final state is checked
against ``engine._run_lipid_lane``'s own.

Also written: a rate probe (``k_looh_decomp_per_min``, ``DecadienalBreakage.k_loss_per_min``, the per-cell and pooled
slate yields, the 1981 ratios) on a (T, Q10, donor) grid and at a perturbed breakage law, and central-difference
Jacobians of log10 hexanal / 2-pentylfuran / nonanal / decadienal in the six module coordinates on four pots, with the
engine's wall-clock for one lane run and for the fits' forward-difference Jacobian.

    python scripts/generators/build_jax_lipid_module_reference.py
Writes results/validation/jax_lipid_module_reference.json.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import config, engine, panel  # noqa: E402
from src.kinetic_core import lipid as lane  # noqa: E402
from src.kinetic_core import parameters_lipid as pl  # noqa: E402
from src.kinetic_core import species_lipid as sl  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "jax_lipid_module_reference.json"
#: save times as fractions of a pot's total duration (each becomes a segment boundary)
SAVE_FRACTIONS = (0.01, 0.1, 0.3, 1.0)
OBS_SPECIES = ("HEXANAL", "PENTYLFURAN", "NONANAL", "DECADIENAL")
COORDINATES = ("q10", "log10_k_looh", "donor", "log10_k_decadienal", "ea_decadienal_kj_mol",
               "share_hexanal_from_decadienal")
#: the interior point the engine's central differences are taken at (donor > 0 so both sides are defined).
#: log10_k_looh stays 0: the 25 C anchor is a constant, not a ConfigKey, so the engine cannot move it without
#: patching a module global (test_engine_config forbids that); its column is checked against the port's own FD only.
THETA_STAR = {"q10": 2.8, "log10_k_looh": 0.0, "donor": 0.25, "log10_k_decadienal": 0.1,
              "ea_decadienal_kj_mol": 30.0, "share_hexanal_from_decadienal": 0.13}
#: coordinates the engine can move (by Config or integrate_lipid's own arguments)
ENGINE_COORDINATES = ("q10", "donor", "log10_k_decadienal", "ea_decadienal_kj_mol", "share_hexanal_from_decadienal")
JACOBIAN_POTS = ("pea_isolate_40C_PratapSingh2021", "pea_isolate_uht_140C_Trikusuma2019",
                 "external_validation_bi_2020_roasted_pea_hexanal", "synthetic_profile_pea_soy_25_160_60")
PROFILE = {"id": "synthetic_profile_pea_soy_25_160_60", "carriers": ("pea_protein_isolate", "soy_protein_isolate"),
           "points": ((0.0, 25.0), (12.0, 160.0), (32.0, 160.0), (40.0, 60.0))}


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def theta_shipped() -> Dict[str, float]:
    brk = lane.DECADIENAL_BREAKAGE
    return {"q10": float(pl.Q10_ASSUMPTION.default), "log10_k_looh": 0.0, "donor": 0.0, "log10_k_decadienal": 0.0,
            "ea_decadienal_kj_mol": float(brk.ea_loss_kj_mol),
            "share_hexanal_from_decadienal": float(brk.share_hexanal)}


@contextlib.contextmanager
def engine_at(theta: Mapping[str, float]):
    """The engine with the module's coordinates applied: the breakage law through the Config (the sanctioned route);
    Q10 and donor are integrate_lipid's own arguments. The 25 C anchor's multiplier is not an engine knob, so it must
    be 0 here."""
    if theta["log10_k_looh"] != 0.0:
        raise ValueError("the engine has no knob for log10_k_looh (the anchor is a constant, not a ConfigKey)")
    brk = config.Config.from_globals()["lipid.decadienal_breakage"]
    new = dataclasses.replace(
        brk, k_loss_per_min_at_ref=brk.k_loss_per_min_at_ref * 10.0 ** theta["log10_k_decadienal"],
        ea_loss_kj_mol=float(theta["ea_decadienal_kj_mol"]),
        share_hexanal=float(theta["share_hexanal_from_decadienal"]))
    with config.applied(config.Config.from_globals().replace({"lipid.decadienal_breakage": new})):
        yield


def refine(segments: Sequence[Tuple[float, float]], fractions=SAVE_FRACTIONS) -> Tuple[List[List[float]], List[float]]:
    """Split a single-hold program at the save fractions (a multi-segment program keeps its own boundaries)."""
    if len(segments) > 1:
        ends = np.cumsum([d for d, _ in segments]).tolist()
        return [[float(d), float(t)] for d, t in segments], ends
    (duration, temp), = segments
    cuts = [duration * f for f in fractions]
    pieces = [[b - a, float(temp)] for a, b in zip([0.0] + cuts[:-1], cuts)]
    return pieces, cuts


def engine_states(carriers: Mapping[str, Mapping[str, Any]], segments, theta: Mapping[str, float]) -> List[List[float]]:
    """The engine's state (lipid species order) after every segment, summed over carriers."""
    branch, composition = engine.core_lipid_model()
    rows = []
    with engine_at(theta):
        for i in range(len(segments)):
            total = np.zeros(len(sl.LIPID_KEYS))
            for key, over in carriers.items():
                charge = lane.charge_from_carrier(pl.LIPID_CARRIERS[key], composition,
                                                  lipid_fraction=over.get("lipid_fraction"),
                                                  peroxide_value_meq_per_kg=over.get("peroxide_value_meq_per_kg"))
                run = lane.integrate_lipid(charge, [tuple(s) for s in segments[: i + 1]], branch,
                                           q10=theta["q10"], donor=theta["donor"])
                total += np.array([run.state_mmol_per_l[k] for k in sl.LIPID_KEYS])
            rows.append(total.tolist())
    return rows


def lipid_panel_pots() -> List[Dict[str, Any]]:
    scored, _ = panel.panel_bundles()
    pots = []
    for path, tag in scored:
        bench = panel.load_bundle(path)
        spec = panel.core_spec(bench)
        try:
            decl = engine.declare_envelope(spec, list(panel.bundle_targets(bench)))
        except Exception:   # noqa: BLE001 - a refused declaration is not a lipid pot
            continue
        if engine.LIPID not in engine._lanes_of(decl):
            continue
        pots.append({"id": str(bench.get("benchmark_id")), "panel": tag, "spec": spec, "declaration": decl,
                     "carriers": [c for c in decl.lipid_carriers if c in pl.LIPID_CARRIERS],
                     "segments": [list(map(float, s)) for s in spec.process.thermal.segments]})
    return pots


def pot_record(pot_id, panel_tag, carriers, segments, theta, variant, lane_final=None) -> Dict[str, Any]:
    refined, save_times = refine(segments)
    states = engine_states(carriers, refined, theta)
    rec = {"id": pot_id, "panel": panel_tag, "variant": variant, "carriers": carriers, "segments": segments,
           "refined_segments": refined, "save_times_min": save_times, "theta": dict(theta), "states": states}
    if lane_final is not None:
        final = np.array([float(lane_final.get(k, 0.0)) for k in sl.LIPID_KEYS])
        mine = np.array(states[-1])
        # the engine's own split-vs-unsplit noise: a depleted pool is computed as state - consumed, which loses
        # ~1e-16 / exp(-k dt) of its relative precision (1.7e-11 on the SPI extrusion pot's upper corner)
        rel = float(np.max(np.abs(mine - final) / np.maximum(np.abs(final), 1e-300)))
        if rel > 1e-9:
            raise AssertionError(f"{pot_id} {variant}: prefix integration differs from _run_lipid_lane by {rel:.2e}")
        rec["engine_lane_final_state"] = final.tolist()
        rec["prefix_vs_lane_max_rel"] = rel
    return rec


def rate_probe() -> Dict[str, Any]:
    branch, composition = engine.core_lipid_model()
    brk = lane.DECADIENAL_BREAKAGE
    temps = [20.0, 25.0, 40.0, 60.0, 99.9, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 230.0]
    q10s = [pl.Q10_ASSUMPTION.lo, pl.Q10_ASSUMPTION.default, pl.Q10_ASSUMPTION.hi]
    donors = [0.0, 0.3, 0.9]
    perturbed = dataclasses.replace(brk, k_loss_per_min_at_ref=0.21, ea_loss_kj_mol=41.0, share_hexanal=0.2)
    cells = {}
    for d in donors:
        per = {}
        for pool in sl.LOOH_POOLS:
            unit = {k: (1.0 if k == pool else 0.0) for k in sl.LOOH_POOLS}
            comp = lane.LOOHComposition(f13_ct=unit["LOOH_13_ct"], f13_tt=unit["LOOH_13_tt"],
                                        f9_ct=unit["LOOH_9_ct"], f9_tt=unit["LOOH_9_tt"])
            y = lane.slate_yields(branch, comp, donor=d)
            per[pool] = {"slate": y, "pentylfuran": lane._b28_yields(y)["PENTYLFURAN"]}
        pooled = lane.slate_yields(branch, composition, donor=d)
        cells[repr(d)] = {"cells": per, "pooled": pooled,
                          "nonanal_per_oleate_looh": lane._b28_yields(pooled)["NONANAL_PER_OLEATE_LOOH"],
                          "pentylfuran_pooled": lane._b28_yields(pooled)["PENTYLFURAN"]}
    return {
        "temperatures_c": temps, "q10": q10s, "donors": donors,
        "k_looh_per_min": {repr(q): [pl.k_looh_decomp_per_min(t, q) for t in temps] for q in q10s},
        "breakage_window_min_c": brk.window_min_c,
        "k_breakage_per_min": [brk.k_loss_per_min(t) for t in temps],
        "breakage_shipped": dataclasses.asdict(brk),
        "k_breakage_perturbed_per_min": [perturbed.k_loss_per_min(t) for t in temps],
        "breakage_perturbed": dataclasses.asdict(perturbed),
        "yields": cells,
    }


def observables(states, rows) -> np.ndarray:
    idx = [sl.LIPID_KEYS.index(s) for s in OBS_SPECIES]
    return np.array([math.log10(max(states[r][i], 1e-300)) for i in idx for r in rows])


def central_jacobian(carriers, refined, theta, rows) -> Dict[str, Any]:
    cols, steps = [], {}
    t0 = time.perf_counter()
    for key in ENGINE_COORDINATES:
        h = 1e-5 * max(1.0, abs(theta[key]))
        up, dn = dict(theta), dict(theta)
        up[key] += h
        dn[key] -= h
        cols.append((observables(engine_states(carriers, refined, up), rows)
                     - observables(engine_states(carriers, refined, dn), rows)) / (2 * h))
        steps[key] = h
    return {"theta": dict(theta), "steps": steps, "columns": list(ENGINE_COORDINATES),
            "jacobian": np.array(cols).T.tolist(), "seconds": time.perf_counter() - t0}


def engine_timing(pot: Dict[str, Any], repeats: int = 20) -> Dict[str, Any]:
    """One lane run as the engine makes it (``_run_lipid_lane``, corners off), and a forward-difference Jacobian in
    seven such runs (the base point and one step per coordinate; the lane's own Q10 argument stands for each step, so
    no Config is installed inside the timing: that would time the install, not the lane)."""
    spec, decl = pot["spec"], pot["declaration"]
    engine._run_lipid_lane(spec, decl, corners=False)   # warm the fit-report cache
    secs = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        engine._run_lipid_lane(spec, decl, corners=False)
        secs.append(time.perf_counter() - t0)
    q10 = float(pl.Q10_ASSUMPTION.default)
    t0 = time.perf_counter()
    for _ in range(repeats):
        for step in range(1 + len(COORDINATES)):
            engine._run_lipid_lane(spec, decl, q10=q10 + 1e-3 * step, corners=False)
    fd = (time.perf_counter() - t0) / repeats
    return {"pot": pot["id"], "lane_run_s": secs, "lane_run_median_s": float(np.median(secs)),
            "forward_fd_jacobian_s": fd, "fd_runs": 1 + len(COORDINATES), "uptime": uptime()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.parse_args(argv)
    q10 = pl.Q10_ASSUMPTION
    shipped = theta_shipped()
    report: Dict[str, Any] = {
        "artifact": OUT.name, "generated_by": "scripts/generators/build_jax_lipid_module_reference.py",
        "env": "maillard (core)", "numpy": np.__version__, "species": list(sl.LIPID_KEYS),
        "coordinates": list(COORDINATES), "theta_shipped": shipped, "theta_star": THETA_STAR,
        "observable_species": list(OBS_SPECIES), "save_fractions": list(SAVE_FRACTIONS),
        "config_digest": config.Config.from_globals().digest(), "uptime_start": uptime()}
    pots: List[Dict[str, Any]] = []
    panel_pots = lipid_panel_pots()
    for pot in panel_pots:
        carriers = {c: {} for c in pot["carriers"]}
        point, meta, _ = engine._run_lipid_lane(pot["spec"], pot["declaration"], corners=True)
        pots.append(pot_record(pot["id"], pot["panel"], carriers, pot["segments"], shipped, "point", point))
        for variant, q, pick in (("lower", q10.lo, "lo"), ("upper", q10.hi, "hi")):
            over = {c: {"lipid_fraction": getattr(pl.LIPID_CARRIERS[c], f"lipid_{pick}"),
                        "peroxide_value_meq_per_kg": getattr(pl.LIPID_CARRIERS[c], f"pv_{pick}")}
                    for c in pot["carriers"]}
            final = meta[f"{variant}_corner_mmol_per_l"]
            pots.append(pot_record(pot["id"], pot["panel"], over, pot["segments"], dict(shipped, q10=q), variant,
                                   final))
    program = engine.ThermalProgram.from_profile(PROFILE["points"])
    profile_segments = [list(map(float, s)) for s in program.segments]
    pots.append(pot_record(PROFILE["id"], "synthetic", {c: {} for c in PROFILE["carriers"]}, profile_segments,
                           shipped, "point"))
    report["pots"] = pots
    report["rate_probe"] = rate_probe()
    jac = {}
    for pot_id in JACOBIAN_POTS:
        rec = next(p for p in pots if p["id"] == pot_id and p["variant"] == "point")
        rows = (list(range(len(rec["refined_segments"]))) if len(rec["segments"]) == 1
                else [len(rec["refined_segments"]) // 4, len(rec["refined_segments"]) // 2,
                      len(rec["refined_segments"]) - 1])
        jac[pot_id] = dict(central_jacobian(rec["carriers"], rec["refined_segments"], THETA_STAR, rows), rows=rows)
    report["jacobian_central_fd"] = jac
    report["timing"] = [engine_timing(p) for p in panel_pots if p["id"] in JACOBIAN_POTS[:3]]
    report["uptime_end"] = uptime()
    OUT.write_text(json.dumps(report, default=float) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(pots)} pot records, {len(panel_pots)} panel pots")
    return 0


if __name__ == "__main__":
    sys.exit(main())

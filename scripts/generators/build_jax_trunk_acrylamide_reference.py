#!/usr/bin/env python
"""
build_jax_trunk_acrylamide_reference.py -- today's engine's numbers for the JAX trunk and acrylamide modules (plan step
A2 + H3, docs/guides/A2_MODULE_CONTRACT.md v2; 2026-11-05). Runs in the CORE env (``./scripts/docker_maillard.sh
run``), never in maillard-jax: every number the port is judged against comes from the package set that backs
results/validation/.

The pots: every panel bundle the engine routes to the acrylamide lane and answers (``panel.panel_bundles`` ->
``engine.declare_envelope`` with the bundle's own targets; 6 on 2026-11-05: the SPI extrusion pot and five
maillard-path hold-outs), each at its point and at two corners of the lane's sampled declared bands (every condition
coordinate at its band's low / high end together: the formation a_w scale, the elimination deficit scale, both pH
exponents, and the CF-039 Ea draws of k_int1_acr and k_asn_glc, the latter at its 95 % CI), plus one SYNTHETIC
staircase profile (not a benchmark): glucose, fructose, asparagine, glycine and all four competitors, 25 -> 180 C in
8 min, 12 min hold, cooled to 100 C, at pH 5.0 and a_w 0.80, so every one of the 42 steps and every condition term
carries flux. The engine has no corner runs of its own for this lane (its envelope samples the bands); these corners
are the bands' ends.

States: after every segment of the programme refined at its save times, by the engine's own per-segment path
(``acrylamide_conditions.apply`` then ``integrate_acrylamide``, chained as ``engine._integrate_program`` chains them)
at LSODA rtol 1e-12 / atol 1e-20. The final state is checked against ``engine._integrate_program`` itself (its rtol
1e-8), and every engine species the port does not carry is checked constant.

Also written: a rate probe of the acrylamide lane (``acrylamide_rate_constants_at`` after ``apply``) on a (T, pH,
a_w) grid at the shipped and at a perturbed point; a rate probe of the trunk LANE's declared terms on the 26 steps of
``network.REACTIONS`` (``trunk_conditions.apply`` then ``network.rate_constants_at``), shipped, perturbed (the exit
exponent through the Config) and with the exit term off (Config); central-difference Jacobians of log10 acrylamide and
HMF in every module coordinate (84) on three pots, moved only through the engine's arguments (a parameter dict, the
``acrylamide_draw``); and the panel predictions (``engine.predict``) with their measured values.

    python scripts/generators/build_jax_trunk_acrylamide_reference.py
Writes results/validation/jax_trunk_acrylamide_module_reference.json.
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
from types import SimpleNamespace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import acrylamide as lane  # noqa: E402
from src.kinetic_core import acrylamide_conditions as ac  # noqa: E402
from src.kinetic_core import config, engine, network, panel, trunk_conditions  # noqa: E402
from src.kinetic_core.parameters_acrylamide import MEASURED_ACRYLAMIDE  # noqa: E402
from src.kinetic_core.species_acrylamide import ACRYLAMIDE_INDEX, ACRYLAMIDE_STATE_KEYS  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "jax_trunk_acrylamide_module_reference.json"
SAVE_FRACTIONS = (0.01, 0.1, 0.3, 1.0)
TIGHT = {"rtol": 1e-12, "atol": 1e-20}
OBS_SPECIES = ("ACR", "HMF")
CELSIUS = 273.15
#: the draw keys of engine._integrate_program's acrylamide branch -> the port's condition coordinates
DRAW_TO_COORDINATE = {"aw_scale": "aw_scale", "aw_elimination_scale": "aw_elimination_scale",
                      "ph_exponent_formation": "ph_exponent_formation",
                      "ph_exponent_elimination": "ph_exponent_elimination"}
PROFILE = {"id": "synthetic_profile_acrylamide_all_steps_25_180_100",
           "initial": {"Glc": 40.0, "Fru": 10.0, "Asn": 20.0, "Gly": 5.0, "Cys": 2.0, "Gln": 5.0, "Lys": 5.0,
                       "Ala": 5.0},
           "points": ((0.0, 25.0), (8.0, 180.0), (20.0, 180.0), (26.0, 100.0)), "ph": 5.0, "water_activity": 0.80}
JACOBIAN_POTS = ("mp_holdout_glucose_asparagine_180C_30min_Chang2021", "acrylamide_spi_extrusion_130C_ACSRef3",
                 PROFILE["id"])
#: the interior point of the Jacobians: the shipped parameters, the condition coordinates off-centre
THETA_STAR_CONDITIONS = {"aw_scale": 1.1, "aw_elimination_scale": 0.8, "ph_exponent_formation": 0.2,
                         "ph_exponent_elimination": 0.14}
FD_STEP = 1e-5


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def corner_draws() -> Dict[str, Dict[str, float]]:
    asn = MEASURED_ACRYLAMIDE["k_asn_glc"]
    lo, hi = {}, {}
    for key, band in (("aw_scale", ac.AW_SCALE_BAND), ("aw_elimination_scale", ac.AW_ELIMINATION_SCALE_BAND),
                      ("ph_exponent_formation", ac.PH_EXPONENT_FORMATION_BAND),
                      ("ph_exponent_elimination", ac.PH_EXPONENT_ELIMINATION_BAND),
                      ("ea_int1", ac.EA_INT1_BAND_KJ_MOL),
                      ("ea_asn_glc", (asn.ea_kj_mol - asn.ea_ci95_kj_mol, asn.ea_kj_mol + asn.ea_ci95_kj_mol))):
        lo[key], hi[key] = float(band[0]), float(band[1])
    return {"lower": lo, "upper": hi}


def shipped_conditions() -> Dict[str, float]:
    return {"aw_scale": float(ac.AW_MULTIPLIER), "aw_elimination_scale": 1.0,
            "ph_exponent_formation": float(ac.PH_EXPONENT_FORMATION),
            "ph_exponent_elimination": float(ac.PH_EXPONENT_ELIMINATION)}


def apply_kwargs(draw: Optional[Mapping[str, Optional[float]]]) -> Dict[str, Any]:
    """engine._integrate_program's acrylamide branch, verbatim in effect: draw -> acrylamide_conditions.apply kwargs."""
    d = draw or {}
    return dict(
        aw_scale=d.get("aw_scale"),
        aw_elimination_deficit_scale=1.0 if d.get("aw_elimination_scale") is None else float(d["aw_elimination_scale"]),
        ph_exponent_formation=(ac.PH_EXPONENT_FORMATION if d.get("ph_exponent_formation") is None
                               else float(d["ph_exponent_formation"])),
        ph_exponent_elimination=(ac.PH_EXPONENT_ELIMINATION if d.get("ph_exponent_elimination") is None
                                 else float(d["ph_exponent_elimination"])),
        ea_overrides={"k_int1_acr": d.get("ea_int1"), "k_asn_glc": d.get("ea_asn_glc")})


def module_parameter_keys() -> List[str]:
    return list(dict.fromkeys(r.parameter_key for r in lane.FULL_ACRYLAMIDE_REACTIONS if r.parameter_key))


def effective_theta(parameters: Mapping[str, Any], draw: Optional[Mapping[str, Optional[float]]]) -> Dict[str, Any]:
    """The port's coordinates for an engine (parameter dict, draw): k_ref and Ea per parameter after the draw's Ea
    overrides (k_ref held), and the four condition coordinates."""
    applied, _ = ac.apply(parameters, SimpleNamespace(ph=None, water_activity=None), **apply_kwargs(draw))
    params = {k: [float(applied[k].k_ref), float(applied[k].ea_kj_mol)] for k in module_parameter_keys()}
    cond = shipped_conditions()
    for k, c in DRAW_TO_COORDINATE.items():
        if draw and draw.get(k) is not None:
            cond[c] = float(draw[k])
    return {"parameters": params, "conditions": cond}


def engine_states(parameters, initial: Mapping[str, float], process, segments, draw, tol=TIGHT) -> List[List[float]]:
    """The engine's state (full acrylamide state order) after every segment: the per-segment path of
    engine._integrate_program's acrylamide branch, at ``tol``."""
    state = dict(initial)
    rows = []
    kwargs = apply_kwargs(draw)
    for duration, temp in segments:
        params, _ = ac.apply(parameters, process, **kwargs)
        run = lane.integrate_acrylamide(params, float(temp) + CELSIUS, state, np.array([0.0, float(duration)]),
                                        water_activity=process.water_activity, **tol)
        state = {k: float(run.concentrations[-1, i]) for i, k in enumerate(ACRYLAMIDE_INDEX)}
        rows.append([state[k] for k in ACRYLAMIDE_STATE_KEYS])
    return rows


def refine(segments: Sequence[Tuple[float, float]], fractions=SAVE_FRACTIONS):
    if len(segments) > 1:
        return [[float(d), float(t)] for d, t in segments], np.cumsum([d for d, _ in segments]).tolist()
    (duration, temp), = segments
    cuts = [duration * f for f in fractions]
    return [[b - a, float(temp)] for a, b in zip([0.0] + cuts[:-1], cuts)], cuts


def acrylamide_panel_pots() -> List[Dict[str, Any]]:
    scored, _ = panel.panel_bundles()
    pots = []
    for path, tag in scored:
        bench = panel.load_bundle(path)
        spec = panel.core_spec(bench)
        targets = panel.bundle_targets(bench)
        try:
            decl = engine.declare_envelope(spec, list(targets))
        except Exception:   # noqa: BLE001 - a refused declaration is not an acrylamide pot
            continue
        if engine.ACRYLAMIDE not in engine._lanes_of(decl) or not decl.is_answerable:
            continue
        pots.append({"id": str(bench.get("benchmark_id")), "panel": tag, "spec": spec, "declaration": decl,
                     "bench": bench, "targets": targets, "initial": dict(decl.mapped_precursors),
                     "segments": [list(map(float, s)) for s in spec.process.thermal.segments]})
    return pots


def pot_record(pot_id, panel_tag, initial, process, segments, parameters, draw, variant, check_engine=True):
    refined, save_times = refine(segments)
    states = engine_states(parameters, initial, process, refined, draw)
    rec = {"id": pot_id, "panel": panel_tag, "variant": variant, "initial": dict(initial),
           "ph": process.ph, "water_activity": process.water_activity, "segments": segments,
           "refined_segments": refined, "save_times_min": save_times, "draw": dict(draw or {}),
           "theta": effective_theta(parameters, draw), "states": states}
    y0 = np.array([float(initial.get(k, 0.0)) for k in ACRYLAMIDE_STATE_KEYS])
    final = np.array(states[-1])
    touched = {s for r in lane.FULL_ACRYLAMIDE_REACTIONS for s in list(r.reactants) + list(r.products)}
    idle = [i for i, k in enumerate(ACRYLAMIDE_STATE_KEYS) if k not in touched]
    drift = float(np.max(np.abs(final[idle] - y0[idle]))) if idle else 0.0
    if drift > 0.0:
        raise AssertionError(f"{pot_id}: an engine species no acrylamide-lane step touches moved by {drift:.2e}")
    rec["untouched_engine_species"] = [ACRYLAMIDE_STATE_KEYS[i] for i in idle]
    if check_engine:
        lane_final, _ = engine._integrate_program(engine.ACRYLAMIDE, parameters, dict(initial), process,
                                                  acrylamide_draw=draw)
        eng = np.array([float(lane_final.get(k, 0.0)) for k in ACRYLAMIDE_STATE_KEYS])
        floor = 1e-9 * float(np.max(y0))
        mask = eng > floor
        rel = float(np.max(np.abs(final[mask] - eng[mask]) / eng[mask]))
        if rel > 1e-4:
            raise AssertionError(f"{pot_id} {variant}: tight per-segment path differs from _integrate_program "
                                 f"(rtol 1e-8) by {rel:.2e}")
        rec["tight_vs_integrate_program_rtol_1e-8_max_rel"] = rel
    return rec


# --------------------------------------------------------------------------- rate probes
GRID_T = [25.0, 60.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 230.0]
GRID_PH = [None, 3.0, 4.0, 5.5, 6.0, 6.8, 7.3, 8.0, 9.0]
GRID_AW = [None, 0.2, 0.33, 0.34, 0.5, 0.8, 0.82, 0.86, 0.92, 0.98, 0.99, 1.0]
PERTURB = {"k_schiff": (1.7, 13.0), "k_tdg_fa": (0.6, -7.0), "k_asn_glc": (1.3, 9.0), "k_int1_acr": (0.8, -11.0),
           "k_acr_dp": (2.1, 5.0), "k_cys_glc": (0.5, 3.0), "k_ama_tdg": (1.2, 4.0)}
PERTURBED_DRAW = {"aw_scale": 1.25, "aw_elimination_scale": 0.6, "ph_exponent_formation": 0.15,
                  "ph_exponent_elimination": 0.13, "ea_int1": 120.0, "ea_asn_glc": 130.0}


def perturbed(parameters: Mapping[str, Any]) -> Dict[str, Any]:
    out = dict(parameters)
    for key, (k_mult, d_ea) in PERTURB.items():
        if key not in out:      # the trunk lane's dict has no acrylamide parameters
            continue
        p = out[key]
        out[key] = dataclasses.replace(p, k_ref=float(p.k_ref) * k_mult, ea_kj_mol=float(p.ea_kj_mol) + d_ea)
    return out


def acrylamide_rate_probe() -> Dict[str, Any]:
    shipped = engine.core_parameters(engine.ACRYLAMIDE)
    keys = [r.key for r in lane.FULL_ACRYLAMIDE_REACTIONS]
    out = {"temperatures_c": GRID_T, "ph": GRID_PH, "aw": GRID_AW, "reactions": keys, "cases": []}
    for label, params, draw in (("shipped", shipped, None), ("perturbed", perturbed(shipped), PERTURBED_DRAW)):
        values = []
        for ph in GRID_PH:
            for aw in GRID_AW:
                applied, _ = ac.apply(params, SimpleNamespace(ph=ph, water_activity=aw), **apply_kwargs(draw))
                for t in GRID_T:
                    k = lane.acrylamide_rate_constants_at(applied, t + CELSIUS)
                    values.append([k[r] for r in keys])
        out["cases"].append({"label": label, "theta": effective_theta(params, draw), "k": values})
    return out


def trunk_lane_rate_probe() -> Dict[str, Any]:
    """trunk_conditions.apply on the trunk lane's parameters, restricted to network.REACTIONS."""
    shipped = engine.core_parameters(engine.TRUNK)
    keys = [r.key for r in network.REACTIONS]
    base = config.Config.from_globals()
    exit_shipped = dict(base["trunk.three_deoxy_exit_ph"])
    (exit_key, (exit_e, exit_band, exit_src)), = exit_shipped.items()
    cases_def = (
        ("shipped", shipped, base, 1.0, trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT),
        ("perturbed", perturbed(shipped),
         base.replace({"trunk.three_deoxy_exit_ph": {exit_key: (0.70, exit_band, exit_src)}}), 0.6, 0.5),
        ("exit_term_off", shipped, base.replace({"trunk.three_deoxy_exit_ph_term": False}), 1.15, 0.8),
    )
    out = {"temperatures_c": GRID_T, "ph": GRID_PH, "aw": GRID_AW, "reactions": keys, "cases": []}
    for label, params, cfg, aw_scale, ph_exponent in cases_def:
        values = []
        with config.applied(cfg):
            for ph in GRID_PH:
                for aw in GRID_AW:
                    proc = SimpleNamespace(ph=ph, water_activity=aw, atmosphere=None)
                    applied, _ = trunk_conditions.apply(params, proc, aw_scale=aw_scale, ph_exponent=ph_exponent)
                    for t in GRID_T:
                        k = network.rate_constants_at(applied, t + CELSIUS)
                        values.append([k[r] for r in keys])
        out["cases"].append({
            "label": label, "values": values, "config_digest": cfg.digest(),
            "exit_term": bool(cfg["trunk.three_deoxy_exit_ph_term"]),
            "exit_ph": {k: v[0] for k, v in dict(cfg["trunk.three_deoxy_exit_ph"]).items()},
            "aw_scale": aw_scale, "ph_exponent": ph_exponent,
            "parameters": {k: [float(params[k].k_ref), float(params[k].ea_kj_mol)]
                           for k in dict.fromkeys(r.parameter_key for r in network.REACTIONS if r.parameter_key)}})
    return out


# --------------------------------------------------------------------------- Jacobians
def observables(states, rows) -> np.ndarray:
    idx = [ACRYLAMIDE_STATE_KEYS.index(s) for s in OBS_SPECIES]
    return np.array([math.log10(max(states[r][i], 1e-300)) for i in idx for r in rows])


def coordinate_list() -> List[str]:
    out = []
    for k in module_parameter_keys():
        out += [f"log10_k:{k}", f"ea:{k}"]
    return out + list(THETA_STAR_CONDITIONS)


def moved(parameters, draw, name: str, h: float):
    """(parameters, draw) with one port coordinate moved by h, through the engine's arguments only."""
    params, d = dict(parameters), dict(draw)
    kind, _, key = name.partition(":")
    if kind == "log10_k":
        params[key] = dataclasses.replace(params[key], k_ref=float(params[key].k_ref) * 10.0 ** h)
    elif kind == "ea":
        if key == "k_int1_acr" and d.get("ea_int1") is not None:
            d["ea_int1"] = float(d["ea_int1"]) + h
        elif key == "k_asn_glc" and d.get("ea_asn_glc") is not None:
            d["ea_asn_glc"] = float(d["ea_asn_glc"]) + h
        else:
            params[key] = dataclasses.replace(params[key], ea_kj_mol=float(params[key].ea_kj_mol) + h)
    else:
        d[name] = float(d[name]) + h
    return params, d


def central_jacobian(initial, process, refined, rows) -> Dict[str, Any]:
    parameters = engine.core_parameters(engine.ACRYLAMIDE)
    draw = dict(THETA_STAR_CONDITIONS)
    names = coordinate_list()
    cols, t0 = [], time.perf_counter()
    for name in names:
        (p_up, d_up), (p_dn, d_dn) = moved(parameters, draw, name, FD_STEP), moved(parameters, draw, name, -FD_STEP)
        cols.append((observables(engine_states(p_up, initial, process, refined, d_up), rows)
                     - observables(engine_states(p_dn, initial, process, refined, d_dn), rows)) / (2 * FD_STEP))
    return {"theta": effective_theta(parameters, draw), "columns": names, "step": FD_STEP, "rows": rows,
            "jacobian": np.array(cols).T.tolist(), "seconds": time.perf_counter() - t0}


# --------------------------------------------------------------------------- the panel's predictions
def panel_predictions(pots) -> List[Dict[str, Any]]:
    out = []
    for pot in pots:
        pred = engine.predict(pot["spec"], list(pot["targets"]))
        rows = []
        for compound, spec in pot["targets"].items():
            key = pot["declaration"].mapped_targets.get(compound)
            value = pred.concentrations_ug_per_l.get(compound)
            measured = panel.measured_value(pot["bench"], compound, spec)
            rows.append({"compound": compound, "species": key, "unit": spec.get("unit", "ppb"),
                         "predicted": value, "measured": measured,
                         "residual_dex": (math.log10(value / measured) if value and measured and value > 0
                                          and measured > 0 else None)})
        out.append({"id": pot["id"], "rows": rows})
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.parse_args(argv)
    t_start = time.perf_counter()
    parameters = engine.core_parameters(engine.ACRYLAMIDE)
    corners = corner_draws()
    report: Dict[str, Any] = {
        "artifact": OUT.name, "generated_by": "scripts/generators/build_jax_trunk_acrylamide_reference.py",
        "env": "maillard (core)", "numpy": np.__version__, "engine_species": list(ACRYLAMIDE_STATE_KEYS),
        "reactions": [r.key for r in lane.FULL_ACRYLAMIDE_REACTIONS], "observable_species": list(OBS_SPECIES),
        "save_fractions": list(SAVE_FRACTIONS), "tight_engine_tolerance": TIGHT, "corner_draws": corners,
        "config_digest": config.Config.from_globals().digest(), "uptime_start": uptime(),
        "trunk_sets": {"network.REACTIONS": len(network.REACTIONS),
                       "network.TRUNK_REACTIONS": len(network.TRUNK_REACTIONS),
                       "acrylamide.ACRYLAMIDE_REACTIONS": len(lane.ACRYLAMIDE_REACTIONS)}}
    from src.kinetic_core import sulfur

    report["trunk_sets"]["sulfur.FULL_REACTIONS"] = len(sulfur.FULL_REACTIONS)
    report["trunk_sets"]["sulfur.FULL_REACTIONS starts with REACTIONS + FED_3DEOXY_REACTIONS"] = bool(
        tuple(r.key for r in sulfur.FULL_REACTIONS[: len(network.REACTIONS) + len(network.FED_3DEOXY_REACTIONS)])
        == tuple(r.key for r in network.REACTIONS + network.FED_3DEOXY_REACTIONS))
    report["trunk_sets"]["network.TRUNK_REACTIONS starts with REACTIONS"] = bool(
        tuple(r.key for r in network.TRUNK_REACTIONS[: len(network.REACTIONS)])
        == tuple(r.key for r in network.REACTIONS))
    pots_out: List[Dict[str, Any]] = []
    panel_pots = acrylamide_panel_pots()
    for pot in panel_pots:
        proc = pot["spec"].process
        pots_out.append(pot_record(pot["id"], pot["panel"], pot["initial"], proc, pot["segments"], parameters,
                                   None, "point"))
        for variant, draw in corners.items():
            pots_out.append(pot_record(pot["id"], pot["panel"], pot["initial"], proc, pot["segments"], parameters,
                                       draw, variant))
        print(f"  {pot['id']}: 3 records, {time.perf_counter() - t_start:.0f} s")
    program = engine.ThermalProgram.from_profile(PROFILE["points"])
    profile_proc = engine.ProcessSpec(thermal=program, ph=PROFILE["ph"], water_activity=PROFILE["water_activity"])
    profile_segments = [list(map(float, s)) for s in program.segments]
    pots_out.append(pot_record(PROFILE["id"], "synthetic", PROFILE["initial"], profile_proc, profile_segments,
                               parameters, None, "point"))
    report["pots"] = pots_out
    print(f"  profile: {len(profile_segments)} segments, {time.perf_counter() - t_start:.0f} s")
    report["acrylamide_rate_probe"] = acrylamide_rate_probe()
    report["trunk_lane_rate_probe"] = trunk_lane_rate_probe()
    print(f"  rate probes, {time.perf_counter() - t_start:.0f} s")
    jac = {}
    by_id = {p["id"]: p for p in panel_pots}
    for pot_id in JACOBIAN_POTS:
        rec = next(p for p in pots_out if p["id"] == pot_id and p["variant"] == "point")
        n = len(rec["refined_segments"])
        rows = list(range(n)) if len(rec["segments"]) == 1 else [n // 4, n // 2, n - 1]
        proc = by_id[pot_id]["spec"].process if pot_id in by_id else profile_proc
        jac[pot_id] = central_jacobian(rec["initial"], proc, rec["refined_segments"], rows)
        print(f"  jacobian {pot_id}: {jac[pot_id]['seconds']:.0f} s")
    report["jacobian_central_fd"] = jac
    report["panel_predictions"] = panel_predictions(panel_pots)
    report["uptime_end"] = uptime()
    report["wall_s"] = time.perf_counter() - t_start
    OUT.write_text(json.dumps(report, default=float) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(pots_out)} pot records, {len(panel_pots)} panel pots, {report['wall_s']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())

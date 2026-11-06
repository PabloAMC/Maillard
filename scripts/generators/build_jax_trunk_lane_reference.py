#!/usr/bin/env python
"""
build_jax_trunk_lane_reference.py -- today's engine's numbers for the JAX trunk LANE (``network.TRUNK_REACTIONS``,
55 steps: the trunk module plus ``trunk_lane`` and ``fed_3deoxy``; plan step A2 + H3, docs/guides/A2_MODULE_CONTRACT.md
v2.1; 2026-11-05). Runs in the CORE env (``./scripts/docker_maillard.sh run``), never in maillard-jax.

THE ENGINE IS MOVED ONLY THROUGH ITS CONFIG AND ITS OWN ARGUMENTS (tasks/lessons.md, 2026-11-05): rate-phase keys
through ``config.applied(Config...replace(...))``; the build-phase methionine / proline values through
``core_parameters(frozen={"methionine": FROZEN_B22, "proline": FROZEN_B24})`` (the engine's own override blocks, which
compute exactly what the module would at ``shipped = True``); the declared bands through ``trunk_conditions.apply``'s
arguments (``aw_scale``, ``ph_exponent``, ``pyrazine_slopes``) and ``core_parameters(frozen={"disputed_sinks": ...})``.

The pots:

* every panel bundle the engine routes to the trunk lane and answers (``panel.panel_bundles`` ->
  ``engine.declare_envelope``; ONE on 2026-11-05, Steinhagen 2021's glucose-only autoclave pot), at its point and at
  two corners of the lane's declared bands (all at the band's low / high end together: the a_w scale, the Amadori pH
  exponent, the 3-DG exit exponent, ENV-B13's disputed sinks and ENV-B34's HPD sinks);
* four SYNTHETIC pots (not benchmarks), so every block carries flux: fed glyoxal + methylglyoxal + glycine (pyrazine),
  glucose + protein-bound lysine (glycation), fed 3-deoxyglucosone (the B39 triangle), and glucose + methionine +
  proline with the B22 / B24 values switched on (the engine's override blocks) and the exit term off (Config);
* one SYNTHETIC staircase profile, 25 -> 160 C in 8 min, 12 min hold, cooled to 100 C, at pH 5.5 and a_w 0.80, with
  glucose, fructose, glycine, methionine, proline and bound lysine charged, at the shipped point and with the switches
  on.

States: after every segment of the programme refined at its save times, by the engine's own per-segment path
(``trunk_conditions.apply`` then ``integrate``, chained as ``engine._integrate_program`` chains them, after its
trunk-lane prelude: Gly += PRO, Gly += MET; LYSP from the protein loading) at LSODA rtol 1e-12 / atol 1e-20. The final
state is checked against ``engine._integrate_program`` itself (its rtol 1e-8).

Also written: a rate probe of all 55 steps (``trunk_conditions.apply`` then ``network.rate_constants_at``) on a
(T, pH, a_w) grid, shipped, perturbed (parameters, the five condition coordinates, the aqueous-glyoxal mapping and
the exit exponent through the Config, B22 / B24 on) and with the exit term off; central-difference Jacobians of log10
observables in every log10 / Ea / condition coordinate on three pots; and the panel predictions (``engine.predict``)
with their post-integration operators (carried volatiles, matrix binding) and measured values.

    ./scripts/docker_maillard.sh run "python scripts/generators/build_jax_trunk_lane_reference.py"
Writes results/validation/jax_trunk_lane_module_reference.json.
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
from src.kinetic_core import config, engine, network, panel, trunk_conditions  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_dicarbonyl import (  # noqa: E402
    DISPUTED_SINK_BANDS, HPD_SINK_BANDS)
from src.kinetic_core.parameters_methionine import FROZEN_B22  # noqa: E402
from src.kinetic_core.parameters_proline import FROZEN_B24  # noqa: E402
from src.kinetic_core.parameters_pyrazine import PYRAZINE_PH_SLOPES  # noqa: E402
from src.kinetic_core.species import SPECIES_KEYS  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "jax_trunk_lane_module_reference.json"
SAVE_FRACTIONS = (0.01, 0.1, 0.3, 1.0)
TIGHT = {"rtol": 1e-12, "atol": 1e-20}
CELSIUS = 273.15
CONDITIONS = ("aw_scale", "ph_exponent", "pyrazine_ph_slope_above_knot", "pyrazine_ph_slope_below_knot")
SWITCHES_ON = {"methionine_b22": True, "proline_b24": True, "exit_term": False}
PROFILE_POINTS = ((0.0, 25.0), (8.0, 160.0), (20.0, 160.0), (26.0, 100.0))
PROFILE_INITIAL = {"Glc": 100.0, "Fru": 20.0, "Gly": 30.0, "MET": 5.0, "PRO": 5.0, "LYSP": 10.0}
#: id -> (initial, (minutes, T_C), pH, a_w, switches); synthetic, chosen so each block carries flux
SYNTHETIC_ISOTHERMAL = {
    "synthetic_fed_glyoxal_methylglyoxal_glycine_120C": ({"GO": 20.0, "MGO": 20.0, "Gly": 40.0}, (60.0, 120.0), 8.0,
                                                         None, {}),
    "synthetic_glucose_bound_lysine_120C": ({"Glc": 150.0, "LYSP": 30.0}, (40.0, 120.0), 6.8, 0.90, {}),
    "synthetic_fed_3deoxyglucosone_120C": ({"TDG": 0.2}, (60.0, 120.0), 5.0, None, {}),
    "synthetic_glucose_methionine_proline_140C_switches_on": ({"Glc": 50.0, "MET": 10.0, "PRO": 10.0},
                                                              (30.0, 140.0), 6.0, None, dict(SWITCHES_ON)),
}
PROFILE_ID = "synthetic_profile_trunk_lane_all_blocks_25_160_100"
JACOBIAN_POTS = {
    "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021": ("HMF", "TDG", "DDG", "G", "GO", "MGO"),
    "synthetic_fed_glyoxal_methylglyoxal_glycine_120C": ("PZ", "DMP", "MPZ", "AKG"),
    PROFILE_ID: ("HMF", "MGO", "GO", "PZ", "MPZ", "CML", "DGAL"),
}
#: the interior point of the Jacobians: shipped parameters, condition coordinates off-centre
THETA_STAR = {"aw_scale": 1.1, "ph_exponent": 0.6, "pyrazine_ph_slope_above_knot": 0.25,
              "pyrazine_ph_slope_below_knot": 0.5, "exit_exponent": 0.8}
FD_STEP = 1e-5


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


# --------------------------------------------------------------------------- Config and engine overrides
def config_for(spec: Mapping[str, Any], *, engine_side: bool) -> "config.Config":
    """The Config a record names: ``exit_term`` / ``exit_exponent`` / ``aqueous_k_ama_g_multiplier`` (rate phase), and,
    on the JAX side only, ``methionine_b22`` / ``proline_b24`` (build phase: the engine gets them as override blocks,
    ``frozen_for``). Shared with the build script, which rebuilds the same Config and checks its digest."""
    base = config.Config.from_globals()
    changes: Dict[str, Any] = {}
    if spec.get("exit_term") is not None:
        changes["trunk.three_deoxy_exit_ph_term"] = bool(spec["exit_term"])
    if spec.get("exit_exponent") is not None:
        changes["trunk.three_deoxy_exit_ph"] = {k: (float(spec["exit_exponent"]),) + tuple(v[1:])
                                                for k, v in dict(base["trunk.three_deoxy_exit_ph"]).items()}
    if spec.get("aqueous_k_ama_g_multiplier") is not None:
        m = dict(base["dicarbonyl.aqueous_glyoxal_parameters"])
        m["k_ama_g"] = dataclasses.replace(m["k_ama_g"],
                                           k_ref=float(m["k_ama_g"].k_ref) * float(spec["aqueous_k_ama_g_multiplier"]))
        changes["dicarbonyl.aqueous_glyoxal_parameters"] = m
    if not engine_side:
        if spec.get("methionine_b22"):
            changes["methionine.shipped_b22"] = True
        if spec.get("proline_b24"):
            changes["proline.shipped_b24"] = True
    return base.replace(changes) if changes else base


def frozen_for(spec: Mapping[str, Any], sinks: Optional[Mapping[str, Any]] = None) -> Optional[Dict[str, Any]]:
    out: Dict[str, Any] = {}
    if spec.get("methionine_b22"):
        out["methionine"] = dict(FROZEN_B22)
    if spec.get("proline_b24"):
        out["proline"] = dict(FROZEN_B24)
    if sinks:
        out["disputed_sinks"] = {k: dict(v) for k, v in sinks.items()}
    return out or None


def sink_corner(end: int) -> Dict[str, Dict[str, float]]:
    """ENV-B13's disputed sinks and ENV-B34's HPD sinks at their band's low (0) or high (1) end, as the envelope's
    ``disputed_sinks`` block states them."""
    shipped = engine.core_parameters(engine.TRUNK)
    out: Dict[str, Dict[str, float]] = {}
    for key, band in DISPUTED_SINK_BANDS.items():
        block = {}
        if band.get("log10_k_100C") is not None:
            block["log10_k_100C"] = float(band["log10_k_100C"][end])
        if band.get("ea_kj_mol") is not None:
            block["ea_kj_mol"] = float(band["ea_kj_mol"][end])
        out[key] = block
    for key, band in HPD_SINK_BANDS.items():
        p = shipped[key]
        rel = float(band["k_rel_hpd"]) * (-1.0 if end == 0 else 1.0)
        hpd = float(band["ea_hpd_kj_mol"]) * (-1.0 if end == 0 else 1.0)
        out[key] = {"log10_k_100C": math.log10(float(p.k_ref) * (1.0 + rel)), "ea_kj_mol": float(p.ea_kj_mol) + hpd}
    return out


def corner_specs() -> Dict[str, Dict[str, Any]]:
    tc = trunk_conditions
    (_, (_, exit_band, _)), = dict(config.Config.from_globals()["trunk.three_deoxy_exit_ph"]).items()
    out = {}
    for name, end in (("lower", 0), ("upper", 1)):
        out[name] = {"conditions": {"aw_scale": float(tc.AW_SCALE_BAND[end]),
                                    "ph_exponent": float(tc.PH_EXPONENT_BAND[end])},
                     "config": {"exit_exponent": float(exit_band[end])}, "sinks": sink_corner(end)}
    return out


def full_conditions(cond: Optional[Mapping[str, float]]) -> Dict[str, float]:
    out = {"aw_scale": 1.0, "ph_exponent": float(trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT),
           "pyrazine_ph_slope_above_knot": float(PYRAZINE_PH_SLOPES[0]),
           "pyrazine_ph_slope_below_knot": float(PYRAZINE_PH_SLOPES[1])}
    out.update({k: float(v) for k, v in (cond or {}).items()})
    return out


def apply_kwargs(cond: Mapping[str, float]) -> Dict[str, Any]:
    return {"aw_scale": cond["aw_scale"], "ph_exponent": cond["ph_exponent"],
            "pyrazine_slopes": (cond["pyrazine_ph_slope_above_knot"], cond["pyrazine_ph_slope_below_knot"])}


def parameter_table(parameters: Mapping[str, Any]) -> Dict[str, List[float]]:
    keys = dict.fromkeys(r.parameter_key for r in network.TRUNK_REACTIONS if r.parameter_key)
    return {k: [float(parameters[k].k_ref), float(parameters[k].ea_kj_mol)] for k in keys}


# --------------------------------------------------------------------------- the engine's states
def engine_initial(initial: Mapping[str, float], process) -> Dict[str, float]:
    """``engine._integrate_program``'s trunk-lane prelude, verbatim in effect (checked against it on every record)."""
    state = {k: float(v) for k, v in initial.items()}
    if state.get("PRO", 0.0) > 0.0:
        state["Gly"] = state.get("Gly", 0.0) + state["PRO"]
    if state.get("MET", 0.0) > 0.0:
        state["Gly"] = state.get("Gly", 0.0) + state["MET"]
    if "LYSP" not in state:
        from src.kinetic_core.matrix_sites import resolve
        from src.kinetic_core.parameters_glycation import available_fraction

        try:
            charged, _ = resolve(process)
        except Exception:  # noqa: BLE001
            charged = None
        if charged is not None and charged.amine > 0:
            state["LYSP"] = float(charged.amine) * available_fraction(charged.amine_band)
    return state


def engine_states(parameters, initial, process, segments, cond, cfg, tol=TIGHT) -> List[List[float]]:
    state = engine_initial(initial, process)
    rows = []
    with config.applied(cfg):
        for duration, temp in segments:
            params, _ = trunk_conditions.apply(parameters, process, **apply_kwargs(cond))
            run = integrate(params, float(temp) + CELSIUS, state, np.array([0.0, float(duration)]), **tol)
            state = {k: float(run.concentrations[-1, i]) for i, k in enumerate(SPECIES_KEYS)}
            rows.append([state[k] for k in SPECIES_KEYS])
    return rows


def refine(segments: Sequence[Tuple[float, float]], fractions=SAVE_FRACTIONS):
    if len(segments) > 1:
        return [[float(d), float(t)] for d, t in segments], np.cumsum([d for d, _ in segments]).tolist()
    (duration, temp), = segments
    cuts = [duration * f for f in fractions]
    return [[b - a, float(temp)] for a, b in zip([0.0] + cuts[:-1], cuts)], cuts


def pot_record(pot_id, panel_tag, initial, process, segments, variant, *, cond=None, cfg_spec=None, sinks=None):
    cfg_spec = dict(cfg_spec or {})
    cond = full_conditions(cond)
    cfg = config_for(cfg_spec, engine_side=True)
    with config.applied(cfg):
        parameters = engine.core_parameters(engine.TRUNK, frozen=frozen_for(cfg_spec, sinks))
    refined, save_times = refine(segments)
    states = engine_states(parameters, initial, process, refined, cond, cfg)
    rec = {"id": pot_id, "panel": panel_tag, "variant": variant, "initial": dict(initial),
           "engine_initial": engine_initial(initial, process), "ph": process.ph,
           "water_activity": process.water_activity, "atmosphere": getattr(process, "atmosphere", None),
           "segments": segments, "refined_segments": refined, "save_times_min": save_times,
           "config_spec": cfg_spec, "engine_config_digest": cfg.digest(),
           "jax_config_digest": config_for(cfg_spec, engine_side=False).digest(),
           "parameters": parameter_table(parameters), "conditions": cond, "states": states}
    # the final state against engine._integrate_program itself (the prelude included; the pyrazine slopes are the
    # shipped ones on every pot record, the only value that path takes)
    if cond["pyrazine_ph_slope_above_knot"] != PYRAZINE_PH_SLOPES[0]:
        raise AssertionError("pot records carry the shipped pyrazine slopes")
    with config.applied(cfg):
        lane_final, _ = engine._integrate_program(engine.TRUNK, parameters, dict(initial), process,
                                                  trunk_draw=(cond["aw_scale"], cond["ph_exponent"]))
    eng = np.array([float(lane_final.get(k, 0.0)) for k in SPECIES_KEYS])
    final = np.array(states[-1])
    floor = 1e-9 * float(max(rec["engine_initial"].values()))
    mask = eng > floor
    rel = float(np.max(np.abs(final[mask] - eng[mask]) / eng[mask]))
    if rel > 1e-4:
        raise AssertionError(f"{pot_id} {variant}: the tight per-segment path differs from _integrate_program "
                             f"(rtol 1e-8) by {rel:.2e}")
    rec["tight_vs_integrate_program_rtol_1e-8_max_rel"] = rel
    return rec


def trunk_panel_pots() -> List[Dict[str, Any]]:
    scored, _ = panel.panel_bundles()
    pots = []
    for path, tag in scored:
        bench = panel.load_bundle(path)
        spec = panel.core_spec(bench)
        targets = panel.bundle_targets(bench)
        try:
            decl = engine.declare_envelope(spec, list(targets))
        except Exception:   # noqa: BLE001 - a refused declaration is not a trunk-lane pot
            continue
        if engine.TRUNK not in engine._lanes_of(decl) or not decl.is_answerable:
            continue
        histories = {}
        for hist in ("isothermal", "as_recorded"):
            with config.applied(config.Config.from_globals().replace({"panel.thermal_history": hist})):
                histories[hist] = [list(map(float, s)) for s in panel.core_spec(bench).process.thermal.segments]
        pots.append({"id": str(bench.get("benchmark_id")), "panel": tag, "spec": spec, "declaration": decl,
                     "bench": bench, "targets": targets, "initial": dict(decl.mapped_precursors),
                     "segments": [list(map(float, s)) for s in spec.process.thermal.segments],
                     "thermal_history_segments": histories})
    return pots


def synthetic_process(segments_or_profile, ph, aw):
    if isinstance(segments_or_profile, tuple) and len(segments_or_profile) == 2 and not isinstance(
            segments_or_profile[0], tuple):
        minutes, t_c = segments_or_profile
        program = engine.ThermalProgram.isothermal(float(t_c), float(minutes))
    else:
        program = engine.ThermalProgram.from_profile(segments_or_profile)
    return engine.ProcessSpec(thermal=program, ph=ph, water_activity=aw)


# --------------------------------------------------------------------------- rate probe
GRID_T = [25.0, 60.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 230.0]
GRID_PH = [None, 3.0, 4.36, 5.5, 6.8, 7.0, 7.4, 8.0, 9.0]
GRID_AW = [None, 0.2, 0.33, 0.34, 0.5, 0.8, 0.85, 0.9, 0.97, 0.98, 0.99, 1.0]
PERTURB = {"k_schiff": (1.7, 13.0), "k_tdg_fa": (0.6, -7.0), "k_ama_tdg": (1.2, 4.0), "k_go_ak": (2.0, 9.0),
           "k_mgo_ak": (0.7, -6.0), "k_cond": (0.3, 5.0), "k_glyc": (1.4, -3.0), "k_flp_cml": (0.8, 6.0),
           "k_g_go": (1.3, 11.0), "k_glc_g": (0.9, -4.0), "k_go_met": (0.5, 2.0), "k_mgo_pro": (1.6, -5.0),
           "k_pyrl_ap": (0.4, 7.0), "k_ddg_dgal": (1.8, 3.0), "k_odg_da": (2.2, -9.0)}
PERTURBED_SPEC = {"methionine_b22": True, "proline_b24": True, "exit_exponent": 0.70, "aqueous_k_ama_g_multiplier": 1.5}
PERTURBED_CONDITIONS = {"aw_scale": 0.6, "ph_exponent": 0.5, "pyrazine_ph_slope_above_knot": 0.31,
                        "pyrazine_ph_slope_below_knot": 0.42}


def perturbed(parameters: Mapping[str, Any]) -> Dict[str, Any]:
    out = dict(parameters)
    for key, (k_mult, d_ea) in PERTURB.items():
        p = out[key]
        out[key] = dataclasses.replace(p, k_ref=float(p.k_ref) * k_mult, ea_kj_mol=float(p.ea_kj_mol) + d_ea)
    return out


def rate_probe() -> Dict[str, Any]:
    keys = list(network.TRUNK_REACTION_KEYS)
    cases_def = (("shipped", {}, None, {}, False),
                 ("perturbed", PERTURBED_SPEC, {"k_da_sink": {"log10_k_100C": -3.0}}, PERTURBED_CONDITIONS, True),
                 ("exit_term_off", {"exit_term": False}, None, {"aw_scale": 1.15, "ph_exponent": 0.8}, False))
    out = {"temperatures_c": GRID_T, "ph": GRID_PH, "aw": GRID_AW, "reactions": keys, "cases": []}
    for label, spec, sinks, cond_in, perturb in cases_def:
        cfg = config_for(spec, engine_side=True)
        cond = full_conditions(cond_in)
        with config.applied(cfg):
            params = engine.core_parameters(engine.TRUNK, frozen=frozen_for(spec, sinks))
            if perturb:
                params = perturbed(params)
            values = []
            for ph in GRID_PH:
                for aw in GRID_AW:
                    proc = SimpleNamespace(ph=ph, water_activity=aw, atmosphere=None)
                    applied, _ = trunk_conditions.apply(params, proc, **apply_kwargs(cond))
                    for t in GRID_T:
                        k = network.rate_constants_at(applied, t + CELSIUS)
                        values.append([k[r] for r in keys])
        out["cases"].append({"label": label, "config_spec": dict(spec), "engine_config_digest": cfg.digest(),
                             "jax_config_digest": config_for(spec, engine_side=False).digest(),
                             "parameters": parameter_table(params), "conditions": cond, "values": values})
    return out


# --------------------------------------------------------------------------- Jacobians
def observables(states, rows, species) -> np.ndarray:
    idx = [SPECIES_KEYS.index(s) for s in species]
    return np.array([math.log10(max(states[r][i], 1e-300)) for i in idx for r in rows])


def jacobian_coordinates(parameters: Mapping[str, Any]) -> List[str]:
    """Every coordinate the engine can move through its arguments: log10 k_ref and Ea of each parameter with a
    non-zero k_ref, the four condition coordinates and the exit exponent (Config). The LINEAR k_ref coordinates of
    the zero-k steps are checked against the port's own central differences only."""
    out = []
    for k, p in parameter_table(parameters).items():
        if p[0] > 0.0:
            out.append(f"log10_k:{k}")
        out.append(f"ea:{k}")
    return out + list(CONDITIONS) + ["ph_exponent_exit:k_tdg_fa"]


def moved(parameters, cond, cfg_spec, name: str, h: float):
    params, c, spec = dict(parameters), dict(cond), dict(cfg_spec)
    kind, _, key = name.partition(":")
    if kind == "log10_k":
        params[key] = dataclasses.replace(params[key], k_ref=float(params[key].k_ref) * 10.0 ** h)
    elif kind == "ea":
        params[key] = dataclasses.replace(params[key], ea_kj_mol=float(params[key].ea_kj_mol) + h)
    elif kind == "ph_exponent_exit":
        spec["exit_exponent"] = float(spec["exit_exponent"]) + h
    else:
        c[name] = float(c[name]) + h
    return params, c, spec


def central_jacobian(rec, process, rows, species) -> Dict[str, Any]:
    cond = full_conditions({k: v for k, v in THETA_STAR.items() if k in CONDITIONS})
    spec = dict(rec["config_spec"], exit_exponent=THETA_STAR["exit_exponent"])
    with config.applied(config_for(spec, engine_side=True)):
        parameters = engine.core_parameters(engine.TRUNK, frozen=frozen_for(spec))
    names = jacobian_coordinates(parameters)
    cols, t0 = [], time.perf_counter()
    for name in names:
        runs = []
        for h in (FD_STEP, -FD_STEP):
            p, c, s = moved(parameters, cond, spec, name, h)
            runs.append(observables(engine_states(p, rec["initial"], process, rec["refined_segments"], c,
                                                  config_for(s, engine_side=True)), rows, species))
        cols.append((runs[0] - runs[1]) / (2 * FD_STEP))
    return {"config_spec": spec, "parameters": parameter_table(parameters), "conditions": cond, "columns": names,
            "species": list(species), "step": FD_STEP, "rows": rows, "jacobian": np.array(cols).T.tolist(),
            "seconds": time.perf_counter() - t0}


# --------------------------------------------------------------------------- the panel's predictions
def panel_predictions(pots) -> List[Dict[str, Any]]:
    out = []
    for pot in pots:
        pred = engine.predict(pot["spec"], list(pot["targets"]))
        meta = pred.run_metadata
        rows = []
        for compound, spec in pot["targets"].items():
            key = pot["declaration"].mapped_targets.get(compound)
            value = pred.concentrations_ug_per_l.get(compound)
            measured = panel.measured_value(pot["bench"], compound, spec)
            binding = (meta.get("matrix_binding") or {}).get(compound)
            rows.append({"compound": compound, "species": key, "unit": spec.get("unit", "ppb"),
                         "predicted": value, "measured": measured,
                         "carried_added": float((meta.get("carried_volatiles") or {}).get(compound, 0.0)),
                         "binding_remaining_fraction": (float(binding["remaining_fraction"]) if binding else 1.0),
                         "residual_dex": (math.log10(value / measured) if value and measured and value > 0
                                          and measured > 0 else None)})
        out.append({"id": pot["id"], "rows": rows})
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.parse_args(argv)
    t_start = time.perf_counter()
    report: Dict[str, Any] = {
        "artifact": OUT.name, "generated_by": "scripts/generators/build_jax_trunk_lane_reference.py",
        "env": "maillard (core)", "numpy": np.__version__, "engine_species": list(SPECIES_KEYS),
        "reactions": list(network.TRUNK_REACTION_KEYS), "save_fractions": list(SAVE_FRACTIONS),
        "tight_engine_tolerance": TIGHT, "config_digest": config.Config.from_globals().digest(),
        "corners": corner_specs(), "uptime_start": uptime(),
        "blocks": {name: [r.key for r in getattr(network, name)] for name in (
            "REACTIONS", "DICARBONYL_REACTIONS", "PYRAZINE_REACTIONS", "GLYCATION_REACTIONS",
            "AQUEOUS_GLYOXAL_REACTIONS", "METHIONINE_REACTIONS", "PROLINE_REACTIONS", "FED_3DEOXY_REACTIONS")}}
    pots_out: List[Dict[str, Any]] = []
    panel_pots = trunk_panel_pots()
    for pot in panel_pots:
        proc = pot["spec"].process
        pots_out.append(pot_record(pot["id"], pot["panel"], pot["initial"], proc, pot["segments"], "point"))
        for variant, corner in report["corners"].items():
            pots_out.append(pot_record(pot["id"], pot["panel"], pot["initial"], proc, pot["segments"], variant,
                                       cond=corner["conditions"], cfg_spec=corner["config"], sinks=corner["sinks"]))
        print(f"  {pot['id']}: 3 records, {time.perf_counter() - t_start:.0f} s")
    report["panel_thermal_history"] = {p["id"]: p["thermal_history_segments"] for p in panel_pots}
    processes = {}
    for pot_id, (initial, seg, ph, aw, spec) in SYNTHETIC_ISOTHERMAL.items():
        proc = synthetic_process(seg, ph, aw)
        processes[pot_id] = proc
        pots_out.append(pot_record(pot_id, "synthetic", initial, proc, [list(map(float, s)) for s in
                                                                        proc.thermal.segments], "point",
                                   cfg_spec=spec))
        print(f"  {pot_id}, {time.perf_counter() - t_start:.0f} s")
    profile_proc = synthetic_process(PROFILE_POINTS, 5.5, 0.80)
    processes[PROFILE_ID] = profile_proc
    profile_segments = [list(map(float, s)) for s in profile_proc.thermal.segments]
    pots_out.append(pot_record(PROFILE_ID, "synthetic_profile", PROFILE_INITIAL, profile_proc, profile_segments,
                               "point"))
    pots_out.append(pot_record(PROFILE_ID, "synthetic_profile", PROFILE_INITIAL, profile_proc, profile_segments,
                               "switches_on", cfg_spec=SWITCHES_ON))
    print(f"  profile: {len(profile_segments)} segments x 2, {time.perf_counter() - t_start:.0f} s")
    report["pots"] = pots_out
    report["rate_probe"] = rate_probe()
    print(f"  rate probe, {time.perf_counter() - t_start:.0f} s")
    jac = {}
    by_id = {p["id"]: p for p in panel_pots}
    for pot_id, species in JACOBIAN_POTS.items():
        rec = next(p for p in pots_out if p["id"] == pot_id and p["variant"] == "point")
        n = len(rec["refined_segments"])
        rows = list(range(n)) if len(rec["segments"]) == 1 else [n // 4, n // 2, n - 1]
        proc = by_id[pot_id]["spec"].process if pot_id in by_id else processes[pot_id]
        jac[pot_id] = central_jacobian(rec, proc, rows, species)
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

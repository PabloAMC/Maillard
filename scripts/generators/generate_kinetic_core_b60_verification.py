#!/usr/bin/env python
"""
Wave B60 -- VERIFICATION OF THE pH-NODE GRID AND THE FLUX BUDGET (2026-10-28).

B60 is a kind-1 correction (``results/validation/kinetic_core_b60_prereg.md`` sec. 4): it ships on
verification, not on scores. This generator runs the pre-registered checks at B51's frozen vector:

  V1  every one of the 52 fit-row observables, ``ph_nodes`` 41 vs 81 on the new grid (and, recorded
      for the record, the fit's own "quick" setting and the old uniform grid);
  V2  three probe pots, the new grid at 41 against a uniform grid fine enough to resolve the transient;
  V3  the flux budget closes on the share-row system and the probe pots;
  V5  (amendment 1) every fit system's pH fixed point converges at the careful setting;
and chooses the optimiser's quick setting as the cheapest candidate whose 52 observables lie within
0.03 dex of the careful setting (amendment 1).

V4 (frozen waves reproduce under ``PH_NODE_GRID = "uniform"``, ``PH_FIXED_POINT = "fixed_sweeps"``) is the
unit-test suite. V2's reference is a 6 401-node uniform grid plus a 16x refinement of the new grid
(321 + 321 nodes); the 25 601-node uniform grid the prereg named needs ~8 converged sweeps per node once
the fixed point is iterated to convergence, which is impractical, and the deviation is recorded.

    python scripts/generators/generate_kinetic_core_b60_verification.py           # ~15 min
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b51_fit as B51  # noqa: E402  (installs B47's rows; B51's physics on)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, ph_state, sulfur  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

B23 = B51.B23
ph_state.B51 = ph_state.B51_ALL_ON

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b60_verification.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b60_verification.md"
B51_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b51_fit_report.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b60_prereg.md"

V1_TOL_DEX = 0.01
V2_TOL_DEX = 0.02
V3_TOL = 0.01
FLOOR = 1e-12
CLOSURE_SPECIES = ("PENT", "Cys", "H2S", "MFT", "FFT")

#: The probe pots of the prereg's sec. 1 (engine path, B51's shipped vector).
PROBE_POTS = {
    "hofmann_ribose_ph5_145C_20min": ({"L-Cysteine": 33.0, "D-Ribose": 100.0}, 145.0, 20.0, 5.0,
                                      BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True,
                                                 source="Hofmann 1998, 0.5 M phosphate")),
    "cao2017_xylose_90C_300min": ({"L-Cysteine": 400.0, "D-Xylose": 400.0}, 90.0, 300.0, 5.5,
                                  BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                                             source="Cao 2017, 0.2 M phosphate")),
    "zhang2025b_xylose_ph65_120C_50min": ({"L-Cysteine": 40.0, "D-Xylose": 40.0}, 120.0, 50.0, 6.5,
                                          BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                                                     source="Zhang 2025b, 0.2 M phosphate")),
}
PROBE_SPECIES = ("MFT", "FFT", "Cys")
BRUTE_FORCE_NODES = 6401
REFINED_NODES = 321
QUICK_TOL_DEX = 0.03
#: (ph_nodes, rtol) candidates for the optimiser's quick setting, cheapest first.
QUICK_CANDIDATES = ((9, 1e-5), (17, 1e-5), (17, 1e-6), (25, 1e-6), (41, 1e-6))


def _observables(x: np.ndarray, quick: bool, nodes: int, grid: str, rtol: float = None) -> Dict[str, float]:
    sulfur.PH_NODE_GRID = grid
    sulfur.PH_FIXED_POINT = "fixed_sweeps" if grid == "uniform" else "relaxed_converged"
    original = B23.integrate_sulfur
    if rtol is not None:
        def _with_rtol(*a, **k):
            k["rtol"] = rtol
            return original(*a, **k)
        B23.integrate_sulfur = _with_rtol
    saved = B23.PH_NODES_QUICK if quick else B23.PH_NODES_CAREFUL
    if quick:
        B23.PH_NODES_QUICK = nodes
    else:
        B23.PH_NODES_CAREFUL = nodes
    try:
        _f, _e, _d, drift = B23.unpack(x)
        predicted, _ = B23.observables(B23.build_parameters(x), quick, drift)
    finally:
        B23.integrate_sulfur = original
        if quick:
            B23.PH_NODES_QUICK = saved
        else:
            B23.PH_NODES_CAREFUL = saved
        sulfur.PH_NODE_GRID, sulfur.PH_FIXED_POINT = "log_and_linear", "relaxed_converged"
    return {k: float(v) for k, v in predicted.items()}


def _dex(a: float, b: float, kind: str) -> float:
    if kind == "ph_endpoint":  # a pH, not a concentration: report the difference in pH units / 10
        return abs(a - b) / 10.0
    if not (math.isfinite(a) and math.isfinite(b)):
        return float("inf")
    return abs(math.log10((abs(a) + FLOOR) / (abs(b) + FLOOR)))


def _probe(name: str, nodes: int, grid: str):
    prec, t_c, minutes, ph, buf = PROBE_POTS[name]
    captured = {}
    original = engine.integrate_sulfur

    def _capture(parameters, temperature_k, state, times, **kw):
        kw["ph_nodes"] = nodes
        captured["args"] = (parameters, temperature_k, dict(state), kw)
        return original(parameters, temperature_k, state, times, **kw)

    sulfur.PH_NODE_GRID = grid
    engine.integrate_sulfur = _capture
    try:
        spec = engine.FormulationSpec(name=name, precursors=prec, process=engine.ProcessSpec(
            thermal=engine.ThermalProgram.isothermal(t_c, minutes), ph=ph, buffer=buf))
        decl = engine.declare_envelope(spec, ["2-methyl-3-furanthiol"])
        state, _ = engine._integrate_program(engine.SULFUR, engine.core_parameters(engine.SULFUR),
                                             dict(decl.mapped_precursors), spec.process)
    finally:
        engine.integrate_sulfur = original
        sulfur.PH_NODE_GRID = "log_and_linear"
    return {k: float(state.get(k, 0.0)) for k in PROBE_SPECIES}, captured["args"]


def _closure(args, minutes: float) -> Dict[str, float]:
    parameters, temperature_k, state, kw = args
    kw = {k: v for k, v in kw.items() if k in ("ph", "ph_final", "buffer_spec", "ph_drift")}
    grid = np.unique(np.concatenate([np.linspace(0.0, minutes, 101), np.geomspace(minutes * 1e-7, minutes, 101)]))
    run = sulfur.integrate_sulfur(parameters, temperature_k, state, grid, **kw)
    flux = sulfur.reaction_flux_integrals(run)
    return sulfur.flux_closure(run, flux, CLOSURE_SPECIES)


def main() -> int:
    assert PREREG.exists(), "B60 is pre-registered; write the prereg before running the verification"
    t0 = time.time()
    x = frozen_vector(json.loads(B51_REPORT.read_text(encoding="utf-8")))
    kinds = {r["id"]: r["kind"] for r in B23.ACTIVE_FIT_ROWS}

    configs = {
        "new_41": (False, 41, "log_and_linear"), "new_81": (False, 81, "log_and_linear"),
        "old_41": (False, 41, "uniform"), "old_quick_9": (True, 9, "uniform"),
    }
    obs, wall = {}, {}
    for label, (quick, nodes, grid) in configs.items():
        t1 = time.time()
        obs[label] = _observables(x, quick, nodes, grid)
        wall[label] = round(time.time() - t1, 1)
        print(f"  V1 {label:12s} done ({time.time() - t0:6.0f}s)", flush=True)
    quick_choice = None
    quick_trials = []
    for nodes, rtol in QUICK_CANDIDATES:
        label = f"new_quick_{nodes}_rtol{rtol:g}"
        t1 = time.time()
        obs[label] = _observables(x, True, nodes, "log_and_linear", rtol=rtol)
        wall[label] = round(time.time() - t1, 1)
        d = max(_dex(obs[label][k], obs["new_81"][k], kinds[k]) for k in obs["new_81"])
        quick_trials.append({"ph_nodes": nodes, "rtol": rtol, "max_dex_vs_careful_81": d, "wall_seconds": wall[label]})
        print(f"  quick {label}: max {d:.4f} dex, {wall[label]} s", flush=True)
        if d <= QUICK_TOL_DEX:
            quick_choice = {"ph_nodes": nodes, "rtol": rtol, "max_dex_vs_careful_81": d}
            break

    def worst(a: str, b: str):
        d = {k: _dex(obs[a][k], obs[b][k], kinds[k]) for k in obs[b]}
        key = max(d, key=d.get)
        return {"max_dex": d[key], "worst_row": key, "rows_over_tol": sorted(k for k, v in d.items() if v > V1_TOL_DEX)}

    v1 = {"tolerance_dex": V1_TOL_DEX, "new_41_vs_new_81": worst("new_41", "new_81"),
          "old_41_vs_new_81": worst("old_41", "new_81"), "old_quick_9_vs_new_81": worst("old_quick_9", "new_81"),
          "wall_seconds": wall}
    v1["pass"] = v1["new_41_vs_new_81"]["max_dex"] <= V1_TOL_DEX

    v2, v3 = {}, {}
    for name, pot in PROBE_POTS.items():
        new, args = _probe(name, 41, "log_and_linear")
        refined, _ = _probe(name, REFINED_NODES, "log_and_linear")
        sulfur.PH_FIXED_POINT = "relaxed_converged"
        ref, _ = _probe(name, BRUTE_FORCE_NODES, "uniform")
        dev = {s: abs(math.log10((new[s] + FLOOR) / (refined[s] + FLOOR))) for s in PROBE_SPECIES}
        dev_uniform = {s: abs(math.log10((new[s] + FLOOR) / (ref[s] + FLOOR))) for s in PROBE_SPECIES}
        v2[name] = {"new_grid_41": new, f"new_grid_{REFINED_NODES}": refined, f"uniform_{BRUTE_FORCE_NODES}": ref,
                    "uniform_node_spacing_min": pot[2] / (BRUTE_FORCE_NODES - 1),
                    "abs_dex_vs_refined": dev, "abs_dex_vs_uniform": dev_uniform,
                    "pass": max(dev.values()) <= V2_TOL_DEX}
        v3[name] = _closure(args, pot[2])
        print(f"  V2/V3 {name} done ({time.time() - t0:6.0f}s)", flush=True)
    share_system = next(r["system"] for r in B23.ACTIVE_FIT_ROWS if r["kind"] == "share")
    spec = B23.SYSTEMS[share_system]
    _f, _e, _d, drift = B23.unpack(x)
    args = (B23.build_parameters(x), float(spec["t_c"]) + 273.15, spec["initial"],
            {"ph": float(spec["ph"]), "buffer_spec": spec.get("buffer"), "ph_drift": drift})
    v3[share_system] = _closure(args, float(spec["minutes"]))
    v3_pass = all(v <= V3_TOL for pot in v3.values() for v in pot.values())

    v5 = {}
    for name, spec in B23.SYSTEMS.items():
        run = sulfur.integrate_sulfur(B23.build_parameters(x), float(spec["t_c"]) + 273.15, spec["initial"],
                                      [0.0, float(spec["minutes"])], ph=float(spec["ph"]),
                                      buffer_spec=spec.get("buffer"), ph_drift=drift, ph_nodes=41,
                                      ph_iterations=4, rtol=1e-8, atol=1e-14)
        if "ph_fixed_point_residual" in run.metadata:
            v5[name] = {"residual_pH": run.metadata["ph_fixed_point_residual"],
                        "sweeps": run.metadata["ph_fixed_point_sweeps"]}
    v5_pass = all(r["residual_pH"] < 1e-3 for r in v5.values())

    verdict = {"V1": v1["pass"], "V2": all(p["pass"] for p in v2.values()), "V3": v3_pass, "V5": v5_pass,
               "quick_setting_found": quick_choice is not None}
    payload = {
        "artifact": "kinetic_core_b60_verification",
        "provenance": provenance.provenance_block(
            "kinetic_core_b60_verification", generated_by="scripts/generators/generate_kinetic_core_b60_verification.py",
            wave="B60", inputs=[B51_REPORT, ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/ph_state.py"]),
        "prereg": data_paths.rel(PREREG),
        "vector": "B51's frozen optimum (kinetic_core_b51_fit_report.json)",
        "V1_fit_panel_convergence": v1, "V2_probe_reference": v2,
        "V3_flux_closure": {"tolerance": V3_TOL, "relative_error": v3, "pass": v3_pass},
        "V5_fixed_point": {"tolerance_pH": 1e-3, "systems": v5, "pass": v5_pass},
        "quick_setting": {"tolerance_dex": QUICK_TOL_DEX, "trials": quick_trials, "choice": quick_choice},
        "observables": obs, "verdict": verdict, "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(payload, indent=2))
    lines = [
        "# Wave B60 verification — the pH-node grid and the flux budget",
        "",
        f"*Generated by `scripts/generators/generate_kinetic_core_b60_verification.py` at B51's frozen vector; "
        f"prereg `{data_paths.rel(PREREG)}` sec. 4.*",
        "",
        "| check | result | pass |",
        "| --- | --- | --- |",
        f"| V1 fit panel, new grid 41 vs 81 | max {v1['new_41_vs_new_81']['max_dex']:.4f} dex ({v1['new_41_vs_new_81']['worst_row']}) | {v1['pass']} |",
        "| quick setting for the fit (cheapest within 0.03 dex of careful 81) | "
        + "; ".join(f"{q['ph_nodes']} nodes rtol {q['rtol']:g}: {q['max_dex_vs_careful_81']:.3f} dex, {q['wall_seconds']} s" for q in quick_trials)
        + f" | {quick_choice is not None} |",
        f"| (record) old uniform 41 vs new 81 | max {v1['old_41_vs_new_81']['max_dex']:.3f} dex ({v1['old_41_vs_new_81']['worst_row']}); {len(v1['old_41_vs_new_81']['rows_over_tol'])} rows over tol | — |",
        f"| (record) old uniform quick 9 vs new 81 | max {v1['old_quick_9_vs_new_81']['max_dex']:.3f} dex ({v1['old_quick_9_vs_new_81']['worst_row']}); {len(v1['old_quick_9_vs_new_81']['rows_over_tol'])} rows over tol | — |",
    ]
    for name, r in v2.items():
        lines.append(f"| V2 {name}: 41 vs {REFINED_NODES} (uniform {BRUTE_FORCE_NODES}) | "
                     + ", ".join(f"{s} {r['abs_dex_vs_refined'][s]:.4f} ({r['abs_dex_vs_uniform'][s]:.3f})" for s in PROBE_SPECIES)
                     + f" dex | {r['pass']} |")
    for name, r in v3.items():
        lines.append(f"| V3 closure {name} | max {max(r.values()):.2e} ({max(r, key=r.get)}) | {max(r.values()) <= V3_TOL} |")
    worst5 = max(v5, key=lambda k: v5[k]["residual_pH"])
    lines.append(f"| V5 fixed point, {len(v5)} systems | worst {v5[worst5]['residual_pH']:.1e} pH ({worst5}, "
                 f"{v5[worst5]['sweeps']} sweeps); max sweeps {max(r['sweeps'] for r in v5.values())} | {v5_pass} |")
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(json.dumps(verdict), f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

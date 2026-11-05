#!/usr/bin/env python
"""
build_jax_bakeoff.py -- the JAX bake-off of the sulfur lane (STATUS_AND_PLAN W4 step H2, pulled forward; 2026-11-04):
can the lane's chemistry be expressed in JAX so that (a) trajectories match today's engine, (b) exact parameter
gradients are correct and much cheaper than the fits' finite-difference Jacobian, (c) the stiff, pH-coupled cases
run?

Runs in the SEPARATE ``maillard-jax`` env (environment-jax.yml), never in the core env:

    source <conda.sh>; conda activate maillard-jax
    XLA_FLAGS=--xla_cpu_multi_thread_eigen=false OMP_NUM_THREADS=1 nice -n 10 \
        python scripts/generators/build_jax_bakeoff.py [--stage NAME ...]

It reads the core engine's numbers from results/validation/jax_bakeoff_sulfur_lane_reference.json (written in the
core env by build_jax_bakeoff_reference.py) and the port from src/kinetic_core_jax/. It installs the same engine
configuration in this process only to READ the lane's static data (lane.build); every compared engine number comes
from the reference file. Writes results/validation/jax_bakeoff_sulfur_lane.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

REFERENCE = ROOT / "results" / "validation" / "jax_bakeoff_sulfur_lane_reference.json"
OUT = ROOT / "results" / "validation" / "jax_bakeoff_sulfur_lane"
CELSIUS = 273.15


def uptime() -> str:
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def setup_pot(reference: Dict[str, Any], overrides: Dict[str, Any] | None = None):
    """The engine configuration of the reference job installed in this process (to read the lane's static data),
    then the port's Lane and Pot. Returns (lane, pot, job)."""
    import build_submodel_probe as probe
    import fit_joint_deoxyosone_thiol as joint
    from src.kinetic_core import config, engine
    from src.kinetic_core.species_sulfur import initial_sulfur_state
    from src.kinetic_core_jax import lane as lane_mod
    from src.kinetic_core_jax import model

    arms = joint.ARMS[reference["arm"]]
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    joint.TDP_EXIT_STATIC = joint.tdp_exit_static(arms)
    joint.STRECKER_SHARE = 0.042
    config.update({"ph.coupling": "inline"})
    job = dict(reference["job"], **(overrides or {}))
    config.install(joint.job_config(job))
    spec = probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"])
    lane = lane_mod.build(engine.core_parameters(engine.SULFUR), reference["keys"], spec, probe._drift(job))
    state0 = dict(job["charge"])
    state0.setdefault("OX", engine.OX_AMBIENT_MMOL_L)
    state0.setdefault("OXR", engine.OX_RESERVOIR_DEFAULT_UNITS)
    pot = model.Pot(lane, initial_sulfur_state(state0), job["t_c"] + CELSIUS, job["label_ph"])
    return lane, pot, job


def timed(fn, *args, repeats: int = 3) -> Dict[str, Any]:
    """Compile (first call) then ``repeats`` timed calls; wall-clock seconds and the machine load around them."""
    import jax

    load_before = uptime()
    t0 = time.perf_counter()
    out = jax.block_until_ready(fn(*args))
    first = time.perf_counter() - t0
    secs = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = jax.block_until_ready(fn(*args))
        secs.append(time.perf_counter() - t0)
    return {"out": out, "first_call_s": first, "compile_s_est": first - min(secs), "runs_s": secs,
            "median_s": float(np.median(secs)), "uptime_before": load_before, "uptime_after": uptime()}


def dex_table(ys, ref_states, times, species, floor=1e-9) -> Dict[str, Any]:
    """|log10 JAX - log10 engine| per species at every time > 0 where either exceeds ``floor`` mmol/L."""
    ys, ref_states = np.asarray(ys), np.asarray(ref_states)
    per = {}
    for i, sp in enumerate(species):
        worst = (0.0, None)
        for r, t in enumerate(times):
            a, b = float(ys[r, i]), float(ref_states[r, i])
            if t <= 0.0 or max(a, b) < floor:
                continue
            d = abs(math.log10(max(a, 1e-300)) - math.log10(max(b, 1e-300)))
            if d > worst[0] or worst[1] is None:
                worst = (d, t)
        if worst[1] is not None:
            per[sp] = {"max_dex": worst[0], "at_min": worst[1]}
    return per


def column_errors(j_test, j_ref, keys, floor=1e-6) -> Dict[str, Any]:
    """Per coordinate: max |test - ref| over the observables, absolute and relative to the column's max |ref|."""
    j_test, j_ref = np.asarray(j_test), np.asarray(j_ref)
    out = {}
    for c, k in enumerate(keys):
        scale = float(np.max(np.abs(j_ref[:, c])))
        diff = float(np.max(np.abs(j_test[:, c] - j_ref[:, c])))
        out[k] = {"max_abs_ref": scale, "max_abs_diff": diff,
                  "relative": diff / scale if scale > floor else None}
    rel = [v["relative"] for v in out.values() if v["relative"] is not None]
    return {"per_coordinate": out, "worst_relative": max(rel) if rel else None,
            "worst_abs_on_flat_columns": max([v["max_abs_diff"] for v in out.values() if v["relative"] is None],
                                             default=None)}


#: 2026-11-04, the first full run (load averages ~17-85): diffrax.BacksolveAdjoint (Kvaerno5, PID rtol 1e-8) on the 17 x 16
#: Jacobian had not returned after ~23 min of one core, against 6.1 s for the discrete adjoint; it was killed and later
#: runs pass --skip-backsolve.
BACKSOLVE_NOTE = ("NOT COMPLETED: the continuous adjoint (BacksolveAdjoint, Kvaerno5, rtol 1e-8) of the 17 x 16 Jacobian "
                  "had not returned after ~23 min on one core (the discrete adjoint: ~6 s); killed. Not pursued: the "
                  "discrete adjoint is exact for the discretised solve and is diffrax's default.")
#: The judgement, written after reading the numbers below (2026-11-04); the tables are generated, this is not.
WHAT_WAS_HARD = """
1. **The pH coupling.** The engine's inline mode solves the charge balance at every right-hand-side call and
   warm-starts each solve from the previous call (a memo). A JAX right-hand side is a pure function, so every call
   solves cold: a 14-step bisection then 6 Newton steps on the coupled (pH, I), twice (100 C and the cooled 25 C
   label). Its derivative is a custom JVP (the implicit-function theorem on the 2 x 2 system, ~25 lines); reverse
   mode transposes it with no extra code. Davies's sqrt(I) is singular at I = 0, so Newton starts after one ionic
   strength sweep. Measured cost split (same solve, pH replaced by constants): **the root-find is ~75 % of the
   forward solve** (2.1 vs 0.5 ms per step). Two attempts to cut it did not: unrolling the loops (no gain), and
   carrying the four roots as extra ODE states to warm-start Newton (slower, 1.86 vs 1.18 s: the roots' own time
   derivative costs two more JVPs per call). Not tried: a cheaper residual (a precomputed pK table, fewer exps) or
   reusing one root across a step's stages.
2. **Discontinuities in the measured pH tables.** Zheng & Ho's cysteine thermolysis switches to its pH 3 pair below
   label pH 5.0, and this pot starts at a cooled pH of 4.9999999995 -- on the boundary (the 10 % jump in k_cys_h2s
   happens at t = 0+ in both codes); the thiamine H2S share has a floor at pH 6.25. Values port exactly (jnp.where);
   an exact derivative ignores the jump-crossing term that a finite difference across a boundary sees. No compared
   column crossed one here, but a pot whose pH walks through 5.0 or 6.25 would show it.
3. **The derived rates.** sulfur_rate_constants_at computes ~25 constants after its generic loop, order-dependent
   (the DDP law, then the Strecker split of the result, the entry from the Amadori steps ...). Ported one by one,
   specialised to the R3En configuration, with a refusal for every other switch. Worth it: 7e-15 relative on 252
   (pH, label pH, I) points x 132 reactions.
4. **The integrator.** diffrax has no BDF / variable-order method; Kvaerno5 (ESDIRK, 7 stages) takes ~670 steps
   here with ~25 % rejected, Kvaerno4 ~640 with almost none and is ~25 % faster; the engine's LSODA is cheaper per
   accuracy. Compile: ~12-16 s for the solve, ~20-30 s for its forward-mode Jacobian; since the pot's numbers are arguments,
   once per configuration, not per pot (the corners reran with no recompile).
5. **The continuous adjoint** (BacksolveAdjoint) did not return in ~23 min on this stiff problem. The discrete adjoint
   (RecursiveCheckpointAdjoint) is exact and works, but with 16 parameters it costs about twice the forward mode, and
   one scalar gradient costs as much as the full Jacobian (one backward pass dominates): forward mode is the right
   mode for this lane's fits; the adjoint matters only once parameters outnumber observables.
6. **Not exercised by this pot:** the five entry coordinates' derivatives are zero (no glycine), as are mgo_share_nf's
   (~3e-7); the gradient check covers the other 11. vmap over 8 parameter points of one pot gave only ~15 % over
   running them one after another: the cost is arithmetic, not dispatch overhead.
"""
RECOMMENDATION = """
**Build A2 + H3 in JAX + diffrax.** Confidence that it is the right call: ~0.7.

What the evidence says:
- (a) **Match: yes.** The port reproduces the engine to ~6e-9 dex on the fitted MFT / FFT rows at the fit's
  tolerance (the engine's own rtol 1e-8 vs 1e-10 difference is the same size), ~2e-11 at tight tolerance, and within
  ~5e-6 dex at every stiff corner, with the final pH to 1e-9. A port can be this faithful because the rate constants and
  the charge balance are checked against the engine BEFORE any integration.
- (b) **Gradients: correct, and more accurate than what the fits use today.** Forward mode, the discrete adjoint and
  central differences of the port agree (forward vs discrete adjoint 6e-15; vs central FD at h = 1e-4: <= 4e-6 relative on every
  column with |J| > 0.01, <= 1.5e-7 absolute everywhere). The port's central FD at the fit's steps equals the engine's
  central FD at the same steps to ~4e-7 relative on every column with |J| > 0.01, so the two compute the same
  function. The fit's own forward-FD
  Jacobian (B8 steps, 0.05 dex) is off by up to 7 % on some columns; the exact Jacobian removes that error.
- **Cost: cheaper, by less than hoped.** One forward solve is ~2.5x SLOWER than the engine's LSODA (1.12 vs
  0.44 s, paired; the pH root, item 1). The full 17 x 16 Jacobian by forward mode is ~2.5x cheaper than the fit's
  17-run FD Jacobian at equal tolerance (3.1 vs 7.8 s, paired), and ~6.5x cheaper with JAX at rtol 1e-6 (1.17 s),
  which exact derivatives allow: the FD Jacobian needs rtol 1e-8 only to keep solver noise below its step, and the
  rtol 1e-6 Jacobian differs from the 1e-8 one by <= 4e-6 relative. PROVISIONAL: shared machine (load ~7-8 on 10
  cores for the paired rounds, up to ~85 earlier); the paired table is the comparison to trust. The plan's "3 - 8x"
  holds only at the loose tolerance; at equal tolerance it is ~2.5x.
- (c) **Stiff, pH-coupled cases: yes.** Every corner (all rate offsets at their upper bounds, label pH 8, cysteine x10,
  unbuffered) ran, matched the engine and gave finite Jacobians; no failures. The only failure was the continuous
  adjoint, which is not needed.

Why JAX over the alternatives:
- **CasADi** (CVODES / IDAS: BDF, the pH as a true algebraic state in IDAS, forward and adjoint sensitivities, C
  codegen) would very likely be faster per forward solve on CPU (my estimate: 2-5x, p ~0.65; NOT measured here) and
  removes the root-find-in-the-RHS problem. But W4's posterior (H4: NumPyro HMC, simulation-based calibration) wants
  gradients inside JAX; a CasADi engine would need a hand-built bridge for every model, and its piecewise pH tables and
  derived rates are a second rewrite of the same chemistry in a less inspectable symbolic language. Choose CasADi
  only if the posterior moves off HMC (for example, Laplace + H0's offset wins H4's comparison).
- **A compiled CPU right-hand side** (numba / C + SUNDIALS) is the fastest forward, but exact sensitivities would
  need a hand-derived or generated Jacobian with respect to every parameter, which is the AD this bake-off shows JAX
  already gives. Not recommended.

Conditions the JAX build should meet (they are what this bake-off learned):
1. Rate constants and the charge balance checked against the engine on grids before any trajectory, per lane, at
   1e-12. A2a's "bit-identical" acceptance cannot hold across integrators (LSODA vs ESDIRK); replace it with "rates to
   1e-12, trajectories to 1e-6 dex at tight tolerance, panel residuals to 1e-4 dex".
2. Pots as arguments (one compile per configuration); structural switches static, fitted numbers in the pytree.
3. Forward mode for least squares and for HMC with ~20 parameters; the discrete adjoint only above ~50.
4. Fits at rtol 1e-6 with exact Jacobians, after a declared check that the panel's residuals move < 1e-3 dex.
5. The pH root-find is the performance budget: give it one focused step (a cheaper residual; one root per step
  reused across stages), with this pot's 2.1 vs 0.5 ms per step as the yardstick.

**Effort, hours of agent work (the JAX route):** the module contract in JAX form (Module emitting arrays, compose,
Config split into static structure and a dynamic pytree), 20-30 h; the trunk, acrylamide and lipid lanes plus every
remaining sulfur switch, each with its rate / charge-balance / trajectory equality tests against the engine, 30-45 h;
thermal profiles and the clamped (pH-stat) and prescribed pH modes, 6-10 h; the fit driver (least squares on
forward-mode Jacobians, per-pot failure handling: max_steps, since SIGALRM cannot interrupt XLA; one compile cache per
configuration), 10-15 h; the pH-root performance step, 8-15 h (outcome uncertain). **Total ~75-115 h**, about the
plan's ~2.5 weeks. NumPyro for H4 is H4's own budget.

Nothing blocks. The separate env works on this machine (CPU, float64); GPU was not tried (H2 says benchmark it before
any Colab run: a 62-species stiff solve is unlikely to gain on a GPU unless thousands of pots or draws are batched).
"""


def _f(x, spec=".2g") -> str:
    return "n/a" if x is None else format(float(x), spec)


def _load(text: str) -> str:
    return text.split("load averages:")[-1].strip() if text else "?"


def write_markdown(r: Dict[str, Any]) -> None:
    tm, g, m = r["timing"], r["gradients"], r["match"]
    uc = r["unit_checks"]
    L = ["# JAX bake-off: the sulfur lane in JAX + diffrax (plan W4 step H2, pulled forward)", "",
         f"*Generated by `scripts/generators/build_jax_bakeoff.py` in the separate `maillard-jax` env "
         f"(jax {r['env']['jax']}, diffrax {r['env']['diffrax']}, float64; {r['env']['cpu']}); the engine's numbers "
         f"from `{r['reference']}` (core env, `build_jax_bakeoff_reference.py`). Pot `{r['pot']}` (ribose 100 + "
         f"cysteine 33 mmol/L, 0.5 M phosphate, label pH 5.0, 100 C, 12 h) at the optimum of `{r['fit']}` (arm "
         f"{r['arm']}: B60 + B58 / B59b / B59a / B61, inline pH coupling, trunk solutes in the balance). Nothing "
         "shipped moves; the port is not imported by `src/kinetic_core`.*", "",
         "## What was ported", "",
         f"- `src/kinetic_core_jax/`: all {r['port']['species']} species and {r['port']['reactions']} reactions of the "
         "lane (stoichiometry and mass-action layout read from the engine), every rate law (Arrhenius with caps, the "
         "six pH-factor kinds with B51's Suleimenov / GKL pKa and Davies activities, Zheng & Ho's cysteine and "
         "glutathione tables, the thiamine share, the cysteine Amadori route, B59a's Strecker H2S share, B61's "
         "hydroxide decay with the MGO share, the pentose + amine entry, the acetic-acid cleavage, the thioether "
         "equilibrium, the Amadori pin), and the inline charge balance (in-situ pH and I at 100 C plus the cooled 25 C "
         f"label pH, over {len(r['port']['titratable_groups'])} titratable groups: "
         f"{', '.join(r['port']['titratable_groups'])}).",
         "- The 16 free coordinates of R3En are the parameter vector; one compiled function serves every pot of the "
         "configuration (temperature, buffer, SID and initial state are arguments).",
         "- **pH coupling: the root is solved inside the right-hand side and differentiated by the implicit-function "
         "theorem** (a custom JVP: du = -(dF/du)^-1 dF/dc, u = (pH, I)), not carried as a DAE state: diffrax has no "
         "DAE solver, a relaxation would change the model, and this is term for term the engine's inline physics. "
         "Forward solve: bisection (the residual is monotone in pH) then Newton on the coupled (pH, I).",
         "- Configuration scope: `lane.build` refuses any switch the port does not implement (it raises rather than "
         "compute another model).", "",
         "## Unit checks (port vs engine, no integration)", "",
         f"- Rate constants: {uc['rate_constants']['grid_points']} (pH, label pH, I) points x "
         f"{uc['rate_constants']['reactions']} reactions, at the optimum and at a perturbed point: worst relative "
         f"difference {_f(uc['rate_constants']['worst']['relative'])} ({uc['rate_constants']['worst']['reaction']}).",
         f"- Charge balance on {uc['charge_balance']['states']} states of the trajectory: max |pH difference| "
         f"{_f(uc['charge_balance']['max_abs_ph_diff'])} (the engine's bisection tolerance is 1e-9), max relative I "
         f"difference {_f(uc['charge_balance']['max_rel_I_diff'])}; SID {uc['charge_balance']['sid_mmol_L']['jax']:.10g} "
         f"vs {uc['charge_balance']['sid_mmol_L']['engine']:.10g} mmol/L.", "",
         "## Match (trajectories at the measured times and a 63-point grid; engine tight = LSODA rtol 1e-10, atol 1e-16)",
         "", f"The engine itself at the fit's tolerance (rtol 1e-8) differs from its tight run by "
         f"{_f(m['engine_fit_tol_vs_engine_tight_max_dex_MFT_FFT'])} dex on the fitted MFT / FFT rows.", "",
         "| JAX solver, tolerance | steps (rejected) | wall s (median of 3) | max dex MFT/FFT fitted rows vs engine tight | "
         "max dex, selected species | worst species above 1e-9 mmol/L | max abs pH diff (in situ / cooled) |",
         "|---|---|---|---|---|---|---|"]
    for name, x in m["runs"].items():
        sel = max(v["max_dex"] for v in x["max_dex_selected_species"].values())
        w = x["max_dex_any_species_above_1e-9"]
        L.append(f"| {name} | {x['steps']} ({x['rejected_steps']}) | {x['timing']['median_s']:.2f} | "
                 f"{_f(x['max_dex_fitted_rows_MFT_FFT'])} | {_f(sel)} | {w['species']} {_f(w['max_dex'])} at "
                 f"{w['at_min']:.3g} min | {_f(x['max_abs_ph_in_situ_diff'])} / {_f(x['max_abs_ph_cooled_diff'])} |")
    L += ["", "Selected species: MFT, FFT, Cys, NF, TTCA, H2S, DPO, PENT. The target was agreement well inside 0.01 dex.",
          "", "## Gradients (d observable / d coordinate; 17 observables: log10 MFT / FFT / Cys / NF at 30, 60, 360, "
          "720 min and the cooled pH at 720 min; 16 coordinates)", "",
          "Worst relative error per column = max over observables of |difference| / max |reference| in that column "
          "(columns whose reference is below 1e-6 everywhere are reported as an absolute difference).", "",
          "| comparison | worst column relative error, all columns | the same, columns with max abs J > 0.01 | "
          "worst absolute error on flat columns |", "|---|---|---|---|"]

    def big(c):
        vals = [v["relative"] for v in c["per_coordinate"].values() if v["max_abs_ref"] > 0.01]
        return max(vals) if vals else None

    for name, c in g["checks"].items():
        L.append(f"| {name.replace('_', ' ')} | {_f(c['worst_relative'])} | {_f(big(c))} | "
                 f"{_f(c['worst_abs_on_flat_columns'])} |")
    for name in ("discrete_adjoint_recursive_checkpoint", "continuous_adjoint_backsolve"):
        x = g.get(name, {})
        if "vs_forward_mode" in x:
            L.append(f"| {name.replace('_', ' ')} vs forward mode | {_f(x['vs_forward_mode']['worst_relative'])} | "
                     f"{_f(big(x['vs_forward_mode']))} | {_f(x['vs_forward_mode']['worst_abs_on_flat_columns'])} |")
        elif "failed" in x:
            L.append(f"| {name.replace('_', ' ')} | FAILED: {x['failed']} | | |")
    L += ["", "The 1e-3-level entries in the all-columns figure are the near-flat columns (offset_k_cys_actz, max |J| "
          "~5e-6; offset_k_h2s_loss, ~4e-4), where a finite difference's noise is the error. Forward mode vs the "
          "engine's central FD at the fit's 0.05-dex steps (~3e-3) is that FD's truncation error: the port's own "
          "central FD at the same steps equals the engine's to ~4e-7 on the columns that matter. The engine's forward "
          "FD -- what the fits use -- is 7 % off its own central FD on some columns."]
    L += ["", "## Timing (PROVISIONAL: taken while other fits loaded this 10-core machine; load averages beside each)", "",
          "| quantity | wall s | load (1 / 5 / 15 min) |", "|---|---|---|",
          f"| engine, one run (LSODA, rtol 1e-8) | {tm['engine_one_run_s']:.3f} | {_load(tm['engine_one_run_uptime'])} |",
          f"| engine, the fit's forward-FD Jacobian ({tm['engine_fd_jacobian_runs']} runs) | "
          f"{tm['engine_fd_jacobian_forward_s']:.2f} | {_load(tm['engine_fd_jacobian_uptime'])} |",
          f"| JAX, one forward solve after JIT (Kvaerno5, rtol 1e-8) | {tm['jax_forward_s']:.3f} | "
          f"{_load(g['forward_solve_observables_only']['uptime_after'])} |",
          f"| JAX, its compile | {tm['jax_forward_compile_s']:.1f} | |",
          f"| JAX, full Jacobian, forward mode | {tm['jax_forward_mode_jacobian_s']:.2f} | "
          f"{_load(g['forward_mode']['timing']['uptime_after'])} |",
          f"| JAX, its compile | {tm['jax_forward_mode_jacobian_compile_s']:.1f} | |",
          f"| JAX, one forward solve, rtol 1e-6 | {tm['jax_forward_loose_1e-6_s']:.3f} | |",
          f"| JAX, full Jacobian, forward mode, rtol 1e-6 (column error vs rtol 1e-8: "
          f"{_f(g['forward_mode_loose_1e-6']['vs_forward_mode_fit_tol']['worst_relative'])}) | "
          f"{tm['jax_forward_mode_jacobian_loose_1e-6_s']:.2f} | |",
          f"| JAX, full Jacobian, discrete adjoint (17 reverse passes) | {_f(tm['jax_discrete_adjoint_jacobian_s'], '.2f')} | |",
          f"| JAX, gradient of one scalar (discrete adjoint) | {_f(tm['jax_one_scalar_gradient_s'], '.2f')} | |", "",
          f"**Ratio, engine FD Jacobian / JAX forward-mode Jacobian: {tm['ratio_engine_fd_over_jax_forward_mode']:.1f}x at "
          f"equal tolerance, {tm['ratio_engine_fd_over_jax_forward_mode_loose_1e-6']:.1f}x with JAX at rtol 1e-6** "
          "(provisional; the timings above are from different moments of a varying load: the paired table below is "
          "the better comparison)."]
    if tm.get("ratio_engine_fd_over_jax_one_scalar_gradient"):
        L.append(f"Engine FD Jacobian / one JAX scalar gradient: {tm['ratio_engine_fd_over_jax_one_scalar_gradient']:.1f}x.")
    if tm.get("paired"):
        L += ["", "Paired timings (engine and JAX back to back, same load window):", "",
              "| when | engine one run s | engine FD Jacobian s | JAX forward s | JAX forward-mode Jacobian s | ratio | "
              "JAX Jacobian at rtol 1e-6 s | ratio | load after |", "|---|---|---|---|---|---|---|---|---|"]
        for p in tm["paired"]:
            L.append(f"| {p['when']} | {p['engine_one_run_s']:.3f} | {p['engine_fd_s']:.2f} | {p['jax_forward_s']:.3f} | "
                     f"{p['jax_jacobian_s']:.2f} | {p['engine_fd_s'] / p['jax_jacobian_s']:.1f}x | "
                     f"{p['jax_jacobian_loose_1e-6_s']:.2f} | {p['engine_fd_s'] / p['jax_jacobian_loose_1e-6_s']:.1f}x | "
                     f"{_load(p['uptime'])} |")
    L += ["", "## Robustness: stiff corners (Kvaerno5, rtol 1e-8; the engine's run of the same corner beside it)", "",
          "| corner | engine s | JAX s | max dex MFT/FFT/Cys/NF vs engine | worst species | final pH in situ JAX / engine | "
          "forward-mode Jacobian s, finite | failures |", "|---|---|---|---|---|---|---|---|"]
    for name, c in r["corners"].items():
        w = c.get("max_dex_any_species_above_1e-9", {})
        ph = c.get("final_ph_in_situ", {})
        fails = "; ".join(str(c[k])[:120] for k in ("engine_failed", "jax_failed", "jacobian_failed") if c.get(k)) or "none"
        L.append(f"| {name} | {c['engine_seconds']:.2f} | {_f(c.get('jax_seconds'), '.2f')} | "
                 f"{_f(c.get('max_dex_MFT_FFT_Cys_NF'))} | {w.get('species', '')} {_f(w.get('max_dex'))} | "
                 f"{_f(ph.get('jax'), '.4f')} / {_f(ph.get('engine'), '.4f')} | "
                 f"{_f(c.get('jax_forward_mode_jacobian_s'), '.2f')}, {c.get('jacobian_finite')} | {fails} |")
    L += ["", "## What was hard", "", WHAT_WAS_HARD.strip(), "", "## Recommendation for A2 + H3", "",
          RECOMMENDATION.strip(), ""]
    OUT.with_suffix(".md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT.with_suffix('.md')}")


def paired_timing() -> int:
    """One JAX timing paired with the engine's latest --timing-only record (same load window, back to back)."""
    import diffrax
    import jax
    import jax.numpy as jnp

    reference = json.loads(REFERENCE.read_text())
    engine_rec = reference["paired_timing"][-1]
    lane, pot, job = setup_pot(reference)
    species = reference["species"]
    obs_times = [0.0] + list(reference["measured_times"])
    obs_index = [species.index(s) for s in ("MFT", "FFT", "Cys", "NF")]
    obs_rows = list(range(1, len(obs_times)))
    theta = jnp.asarray(reference["theta"])

    def obs_fn(th, adjoint=None, rtol=1e-8):
        sol = pot.solve(th, obs_times, rtol=rtol, atol=rtol * 1e-6, adjoint=adjoint, solver=diffrax.Kvaerno5())
        return pot.observables(sol.ys, obs_index, obs_rows)

    fwd = timed(jax.jit(obs_fn), theta, repeats=5)
    jac = timed(jax.jit(jax.jacfwd(lambda th: obs_fn(th, diffrax.ForwardMode()))), theta, repeats=3)
    jac_loose = timed(jax.jit(jax.jacfwd(lambda th: obs_fn(th, diffrax.ForwardMode(), 1e-6))), theta, repeats=3)
    report = json.loads(OUT.with_suffix(".json").read_text())
    report["timing"].setdefault("paired", []).append({
        "when": engine_rec["when"], "engine_one_run_s": float(np.median(engine_rec["single_run_s"])),
        "engine_fd_s": engine_rec["fd_forward_jacobian_s"], "engine_uptime": engine_rec["uptime_after"],
        "jax_forward_s": fwd["median_s"], "jax_jacobian_s": jac["median_s"],
        "jax_jacobian_loose_1e-6_s": jac_loose["median_s"], "uptime": jac_loose["uptime_after"]})
    OUT.with_suffix(".json").write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    print(report["timing"]["paired"][-1])
    write_markdown(report)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--skip-backsolve", action="store_true")
    parser.add_argument("--markdown-only", action="store_true", help="rewrite the .md from the existing .json")
    parser.add_argument("--paired-timing", action="store_true",
                        help="time the JAX forward solve and forward-mode Jacobian right after the engine's latest "
                             "`build_jax_bakeoff_reference.py --timing-only` record, and add the pair to the .json")
    args = parser.parse_args(argv)
    if args.markdown_only:
        write_markdown(json.loads(OUT.with_suffix(".json").read_text()))
        return 0
    if args.paired_timing:
        return paired_timing()
    import diffrax
    import jax
    import jax.numpy as jnp

    from src.kinetic_core_jax import model

    reference = json.loads(REFERENCE.read_text())
    lane, pot, job = setup_pot(reference)
    keys = reference["keys"]
    theta = jnp.asarray(reference["theta"])
    species = reference["species"]
    measured = reference["measured_times"]
    obs_times = [0.0] + list(measured)
    obs_index = [species.index(s) for s in ("MFT", "FFT", "Cys", "NF")]
    obs_rows = list(range(1, len(obs_times)))
    report: Dict[str, Any] = {
        "artifact": OUT.name + ".json", "generated_by": "scripts/generators/build_jax_bakeoff.py",
        "env": {"name": "maillard-jax", "jax": jax.__version__, "diffrax": diffrax.__version__,
                "numpy": np.__version__, "x64": bool(jax.config.read("jax_enable_x64")),
                "XLA_FLAGS": os.environ.get("XLA_FLAGS"), "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
                "cpu": subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True,
                                      text=True).stdout.strip()},
        "reference": REFERENCE.name, "pot": reference["pot"], "fit": reference["fit"], "arm": reference["arm"],
        "keys": keys, "observables": reference["observables"], "uptime_start": uptime(),
        "port": {"species": lane.n_species, "reactions": len(lane.reactions), "titratable_groups": list(lane.group_names),
                 "notes": list(lane.notes)}}

    # ---------------------------------------------------------------- 0. unit checks: rate constants, charge balance
    worst_k = {"relative": 0.0}
    for name in ("theta_star", "theta_perturbed"):
        probe = reference["rate_constant_probe"][name]
        th = jnp.asarray(probe.get("theta", reference["theta"]))
        f = jax.jit(lambda ph, phl, i, th=th: model.rate_constants(th, pot.data["t_k"], ph, phl, i, lane))
        for (ph, phl, i), kref in zip(probe["grid_ph_phlabel_I"], probe["k"]):
            k, kref = np.asarray(f(ph, phl, i)), np.asarray(kref)
            rel = np.where((k == 0) & (kref == 0), 0.0, np.abs(k - kref) / np.maximum(np.abs(kref), 1e-300))
            j = int(np.argmax(rel))
            if rel[j] > worst_k["relative"]:
                worst_k = {"relative": float(rel[j]), "reaction": lane.reactions[j], "ph": ph, "ph_label": phl,
                           "I_M": i, "theta": name}
    cb = reference["charge_balance_probe"]
    ph_fn = jax.jit(pot.ph_state)
    worst_ph = 0.0
    worst_i = 0.0
    for pt in cb["points"]:
        hot_ph, hot_i, cool_ph, cool_i = map(float, ph_fn(jnp.asarray(pt["state"])))
        worst_ph = max(worst_ph, abs(hot_ph - pt["hot"][0]), abs(cool_ph - pt["cool"][0]))
        worst_i = max(worst_i, abs(hot_i / pt["hot"][1] - 1.0), abs(cool_i / pt["cool"][1] - 1.0))
    report["unit_checks"] = {
        "rate_constants": {"grid_points": 2 * len(reference["rate_constant_probe"]["theta_star"]["k"]),
                           "reactions": len(lane.reactions), "worst": worst_k},
        "charge_balance": {"states": len(cb["points"]), "max_abs_ph_diff": worst_ph, "max_rel_I_diff": worst_i,
                           "sid_mmol_L": {"jax": float(pot.data["sid"]), "engine": cb["sid_mmol_L"]}}}
    print("unit checks", report["unit_checks"])

    # ---------------------------------------------------------------- 1. match
    tight = reference["trajectory_tight_tol"]
    probe_tol = reference["trajectory_probe_tol"]
    times = tight["times"]
    solvers = {"Kvaerno5": diffrax.Kvaerno5(), "Kvaerno4": diffrax.Kvaerno4()}
    match: Dict[str, Any] = {}
    for sname, solver in solvers.items():
        for tol_name, (rtol, atol) in {"fit_tol_1e-8": (1e-8, 1e-14), "tight_1e-10": (1e-10, 1e-16),
                                       "loose_1e-6": (1e-6, 1e-12)}.items():
            def run(th, solver=solver, rtol=rtol, atol=atol):
                sol = pot.solve(th, times, rtol=rtol, atol=atol, solver=solver)
                return sol.ys, sol.stats["num_steps"], sol.stats["num_rejected_steps"]
            t = timed(jax.jit(run), theta)
            ys, steps, rejected = t.pop("out")
            ys = np.asarray(ys)
            ph = [list(map(float, ph_fn(jnp.asarray(y)))) for y in ys]
            vs_tight = dex_table(ys, tight["states"], times, species)
            obs = {sp: vs_tight[sp] for sp in ("MFT", "FFT", "Cys", "NF", "TTCA", "H2S", "DPO", "PENT")}
            obs_meas = max(abs(math.log10(ys[times.index(tm), species.index(sp)])
                               - math.log10(np.asarray(tight["states"])[times.index(tm), species.index(sp)]))
                           for sp in ("MFT", "FFT") for tm in measured)
            allmax = max(vs_tight.items(), key=lambda kv: kv[1]["max_dex"])
            match[f"{sname} {tol_name}"] = {
                "rtol": rtol, "atol": atol, "steps": int(steps), "rejected_steps": int(rejected),
                "timing": t, "max_dex_fitted_rows_MFT_FFT": obs_meas, "max_dex_selected_species": obs,
                "max_dex_any_species_above_1e-9": {"species": allmax[0], **allmax[1]},
                "max_abs_ph_in_situ_diff": float(max(abs(p[0] - q) for p, q in zip(ph, tight["ph_in_situ"]))),
                "max_abs_ph_cooled_diff": float(max(abs(p[2] - q) for p, q in zip(ph, tight["ph_cooled"]))),
                "vs_engine_fit_tol_max_dex_MFT_FFT": max(
                    abs(math.log10(ys[times.index(tm), species.index(sp)])
                        - math.log10(np.asarray(probe_tol["states"])[times.index(tm), species.index(sp)]))
                    for sp in ("MFT", "FFT") for tm in measured)}
            print(sname, tol_name, match[f"{sname} {tol_name}"]["max_dex_fitted_rows_MFT_FFT"], t["median_s"], int(steps))
    eng_self = max(abs(math.log10(np.asarray(probe_tol["states"])[times.index(tm), species.index(sp)])
                       - math.log10(np.asarray(tight["states"])[times.index(tm), species.index(sp)]))
                   for sp in ("MFT", "FFT") for tm in measured)
    report["match"] = {"engine_fit_tol_vs_engine_tight_max_dex_MFT_FFT": eng_self, "runs": match,
                       "engine_single_run_s": reference["single_run_seconds"]}

    # ---------------------------------------------------------------- 2. gradients
    solver = diffrax.Kvaerno5()

    def obs_fn(th, adjoint=None, rtol=1e-8, atol=1e-14, max_steps=200_000):
        sol = pot.solve(th, obs_times, rtol=rtol, atol=atol, adjoint=adjoint, solver=solver, max_steps=max_steps)
        return pot.observables(sol.ys, obs_index, obs_rows)

    grads: Dict[str, Any] = {}
    f_obs = jax.jit(lambda th: obs_fn(th))
    t_fwd = timed(f_obs, theta)
    o0 = np.asarray(t_fwd.pop("out"))
    grads["forward_solve_observables_only"] = t_fwd
    jac_fwd = jax.jit(jax.jacfwd(lambda th: obs_fn(th, adjoint=diffrax.ForwardMode())))
    t = timed(jac_fwd, theta)
    J_fwd = np.asarray(t.pop("out"))
    grads["forward_mode"] = {"timing": t}
    print("forward mode", t["median_s"])
    # the operating point exact gradients allow: a fit's FD Jacobian needs rtol 1e-8 to keep solver noise below its
    # 0.05-dex step; an exact derivative does not, so the forward-mode Jacobian at rtol 1e-6 is the fair rival too
    t = timed(jax.jit(lambda th: obs_fn(th, rtol=1e-6, atol=1e-12)), theta)
    o_loose = np.asarray(t.pop("out"))
    grads["forward_solve_loose_1e-6"] = dict(t, max_abs_obs_diff_vs_fit_tol=float(np.max(np.abs(o_loose - o0))))
    t = timed(jax.jit(jax.jacfwd(lambda th: obs_fn(th, adjoint=diffrax.ForwardMode(), rtol=1e-6, atol=1e-12))), theta)
    J_loose = np.asarray(t.pop("out"))
    grads["forward_mode_loose_1e-6"] = {"timing": t, "vs_forward_mode_fit_tol": column_errors(J_loose, J_fwd, keys)}
    print("forward mode loose", t["median_s"])
    jac_rev = jax.jit(jax.jacrev(lambda th: obs_fn(th, adjoint=diffrax.RecursiveCheckpointAdjoint(), max_steps=20_000)))
    try:
        t = timed(jac_rev, theta, repeats=2)
        J_rev = np.asarray(t.pop("out"))
        grads["discrete_adjoint_recursive_checkpoint"] = {"timing": t,
                                                          "vs_forward_mode": column_errors(J_rev, J_fwd, keys)}
        print("discrete adjoint", t["median_s"])
    except Exception as exc:  # noqa: BLE001
        grads["discrete_adjoint_recursive_checkpoint"] = {"failed": f"{type(exc).__name__}: {str(exc)[:400]}"}
    # one vector-Jacobian product: what an adjoint costs per scalar objective (a fit's cost function)
    w = jnp.ones(len(o0))
    vjp_fn = jax.jit(lambda th: jax.grad(lambda x: jnp.dot(w, obs_fn(
        x, adjoint=diffrax.RecursiveCheckpointAdjoint(), max_steps=20_000)))(th))
    try:
        t = timed(vjp_fn, theta, repeats=2)
        t.pop("out")
        grads["discrete_adjoint_one_gradient_of_a_scalar"] = {"timing": t}
    except Exception as exc:  # noqa: BLE001
        grads["discrete_adjoint_one_gradient_of_a_scalar"] = {"failed": f"{type(exc).__name__}: {str(exc)[:400]}"}
    if args.skip_backsolve:
        grads["continuous_adjoint_backsolve"] = {"failed": BACKSOLVE_NOTE}
    else:
        bs = diffrax.BacksolveAdjoint(solver=diffrax.Kvaerno5(),
                                      stepsize_controller=diffrax.PIDController(rtol=1e-8, atol=1e-14, pcoeff=0.3,
                                                                                icoeff=0.4))
        jac_bs = jax.jit(jax.jacrev(lambda th: obs_fn(th, adjoint=bs, max_steps=50_000)))
        try:
            t = timed(jac_bs, theta, repeats=1)
            J_bs = np.asarray(t.pop("out"))
            grads["continuous_adjoint_backsolve"] = {"timing": t, "vs_forward_mode": column_errors(J_bs, J_fwd, keys)}
            print("backsolve", t["median_s"])
        except Exception as exc:  # noqa: BLE001
            grads["continuous_adjoint_backsolve"] = {"failed": f"{type(exc).__name__}: {str(exc)[:400]}"}
    # checks: central FD of the JAX forward model (small step, tight tolerance), and the engine's own FD
    f_tight = jax.jit(lambda th: obs_fn(th, rtol=1e-11, atol=1e-17))
    jac_fwd_tight = jax.jit(jax.jacfwd(lambda th: obs_fn(th, adjoint=diffrax.ForwardMode(), rtol=1e-11, atol=1e-17)))
    J_fwd_tight = np.asarray(jac_fwd_tight(theta))
    h = 1e-4
    cols = []
    for c in range(len(keys)):
        e = jnp.zeros(len(keys)).at[c].set(h)
        cols.append((np.asarray(f_tight(theta + e)) - np.asarray(f_tight(theta - e))) / (2 * h))
    J_fd_jax = np.array(cols).T
    steps = np.array([reference["fd_central"]["steps"][k] for k in keys])
    cols = []
    for c in range(len(keys)):
        e = jnp.zeros(len(keys)).at[c].set(steps[c])
        cols.append((np.asarray(f_obs(theta + e)) - np.asarray(f_obs(theta - e))) / (2 * steps[c]))
    J_fd_jax_fitsteps = np.array(cols).T
    J_eng_c = np.asarray(reference["fd_central"]["jacobian"])
    J_eng_f = np.asarray(reference["fd_forward"]["jacobian"])
    grads["checks"] = {
        "forward_mode_tight_vs_jax_central_fd_h1e-4": column_errors(J_fwd_tight, J_fd_jax, keys),
        "forward_mode_fit_tol_vs_forward_mode_tight": column_errors(J_fwd, J_fwd_tight, keys),
        "jax_central_fd_at_fit_steps_vs_engine_central_fd_at_fit_steps": column_errors(J_fd_jax_fitsteps, J_eng_c, keys),
        "forward_mode_vs_engine_central_fd_at_fit_steps": column_errors(J_fwd, J_eng_c, keys),
        "engine_forward_fd_vs_engine_central_fd": column_errors(J_eng_f, J_eng_c, keys),
        "forward_mode_vs_engine_forward_fd": column_errors(J_fwd, J_eng_f, keys)}
    grads["jacobians"] = {"forward_mode": J_fwd.tolist(), "forward_mode_tight": J_fwd_tight.tolist(),
                          "jax_central_fd_h1e-4": J_fd_jax.tolist()}
    grads["observables_at_theta_star"] = {"jax": o0.tolist()}
    report["gradients"] = grads

    # ---------------------------------------------------------------- 3. timing summary
    eng_fd = reference["fd_forward"]
    report["timing"] = {
        "PROVISIONAL": "taken while other fits ran on this machine; see each entry's uptime",
        "engine_one_run_s": float(np.median(reference["single_run_seconds"]["runs"])),
        "engine_one_run_uptime": reference["single_run_seconds"]["uptime"],
        "engine_fd_jacobian_forward_s": eng_fd["seconds"], "engine_fd_jacobian_runs": eng_fd["engine_runs"],
        "engine_fd_jacobian_uptime": eng_fd["uptime"],
        "jax_forward_s": t_fwd["median_s"], "jax_forward_compile_s": t_fwd["compile_s_est"],
        "jax_forward_mode_jacobian_s": grads["forward_mode"]["timing"]["median_s"],
        "jax_forward_mode_jacobian_compile_s": grads["forward_mode"]["timing"]["compile_s_est"],
        "jax_forward_loose_1e-6_s": grads["forward_solve_loose_1e-6"]["median_s"],
        "jax_forward_mode_jacobian_loose_1e-6_s": grads["forward_mode_loose_1e-6"]["timing"]["median_s"],
        "jax_discrete_adjoint_jacobian_s": grads.get("discrete_adjoint_recursive_checkpoint", {}).get("timing", {}).get("median_s"),
        "jax_one_scalar_gradient_s": grads.get("discrete_adjoint_one_gradient_of_a_scalar", {}).get("timing", {}).get("median_s"),
    }
    tm = report["timing"]
    tm["ratio_engine_fd_over_jax_forward_mode"] = tm["engine_fd_jacobian_forward_s"] / tm["jax_forward_mode_jacobian_s"]
    tm["ratio_engine_fd_over_jax_forward_mode_loose_1e-6"] = (tm["engine_fd_jacobian_forward_s"]
                                                              / tm["jax_forward_mode_jacobian_loose_1e-6_s"])
    if tm["jax_one_scalar_gradient_s"]:
        tm["ratio_engine_fd_over_jax_one_scalar_gradient"] = tm["engine_fd_jacobian_forward_s"] / tm["jax_one_scalar_gradient_s"]
    print("timing", tm)

    # ---------------------------------------------------------------- 4. corners
    corners: Dict[str, Any] = {}
    run_corner = jax.jit(lambda th, data, y0: pot.solve(th, times, data=data, y0=y0, rtol=1e-8, atol=1e-14,
                                                        solver=solver).ys)
    jac_corner = jax.jit(jax.jacfwd(lambda th, data, y0: pot.observables(
        pot.solve(th, obs_times, data=data, y0=y0, rtol=1e-8, atol=1e-14, solver=solver,
                  adjoint=diffrax.ForwardMode()).ys, obs_index, obs_rows, data)))
    from src.kinetic_core import engine as _engine
    from src.kinetic_core.species_sulfur import initial_sulfur_state
    import build_submodel_probe as probe

    for name, rec in reference["corners"].items():
        over = rec["overrides"]
        cj = dict(job, **over)
        state0 = dict(cj["charge"])
        state0.setdefault("OX", _engine.OX_AMBIENT_MMOL_L)
        state0.setdefault("OXR", _engine.OX_RESERVOIR_DEFAULT_UNITS)
        y0 = jnp.asarray(initial_sulfur_state(state0))
        spec = probe.buffer_of({"buffer": cj["buffer"]}, cj["label_ph"])
        buffer = np.zeros(len(lane.group_names))
        for group, mmol in spec.groups():
            buffer[lane.group_names.index(f"buffer_{group.key}")] = mmol
        data = model.pot_data(lane, y0, cj["t_c"] + CELSIUS, cj["label_ph"], buffer)
        th = jnp.asarray(rec["theta"])
        out: Dict[str, Any] = {"overrides": over, "engine_seconds": rec["seconds"], "engine_failed": rec.get("failed"),
                               "sid": {"jax": float(data["sid"]), "engine": rec.get("sid_mmol_L")}}
        try:
            t = timed(run_corner, th, data, y0, repeats=1)
            ys = np.asarray(t.pop("out"))
            out["jax_seconds"] = t["median_s"]
            out["jax_first_call_s"] = t["first_call_s"]
            if "states" in rec:
                tab = dex_table(ys, rec["states"], times, species)
                out["max_dex_MFT_FFT_Cys_NF"] = max(tab[s]["max_dex"] for s in ("MFT", "FFT", "Cys", "NF"))
                allmax = max(tab.items(), key=lambda kv: kv[1]["max_dex"])
                out["max_dex_any_species_above_1e-9"] = {"species": allmax[0], **allmax[1]}
                ph_end = [float(v) for v in pot.ph_state(jnp.asarray(ys[-1]), data)]
                out["final_ph_in_situ"] = {"jax": ph_end[0], "engine": rec["ph_in_situ"][-1]}
                out["final_ph_cooled"] = {"jax": ph_end[2], "engine": rec["ph_cooled"][-1]}
        except Exception as exc:  # noqa: BLE001
            out["jax_failed"] = f"{type(exc).__name__}: {str(exc)[:400]}"
        try:
            t = timed(jac_corner, th, data, y0, repeats=1)
            Jc = np.asarray(t.pop("out"))
            out["jax_forward_mode_jacobian_s"] = t["median_s"]
            out["jacobian_finite"] = bool(np.all(np.isfinite(Jc)))
        except Exception as exc:  # noqa: BLE001
            out["jacobian_failed"] = f"{type(exc).__name__}: {str(exc)[:400]}"
        out["uptime"] = uptime()
        corners[name] = out
        print("corner", name, {k: v for k, v in out.items() if k != "overrides"})
    report["corners"] = corners
    report["uptime_end"] = uptime()
    OUT.with_suffix(".json").write_text(json.dumps(report, default=float, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.with_suffix('.json')}")
    write_markdown(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())

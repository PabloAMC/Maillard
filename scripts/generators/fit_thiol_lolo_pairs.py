#!/usr/bin/env python
"""
fit_thiol_lolo_pairs.py -- H0's leave-two-laboratories-out refits of the thiol lane, writing residuals only.

Pre-registration: results/validation/kinetic_core_between_lab_tau_h0_prereg.md (v3), secs. 3.1, 4.2, 4.3 and 7. For
each unordered pair of the five scorable laboratories, R3En is refitted without BOTH, under EX-056's settings
exactly (warm start from s042_r3en_fixed, s042_e_conv's fixed series noise, budget 1200, ftol 1e-4, up to 3
continuations, the Laplace on the training rows, the leak-prone rule and its reset arm). The record holds the
optimiser fields and the RAW residuals of both held-out laboratories and of the training rows; nothing else.

Blinding (prereg sec. 4.3): this script never computes, writes or prints a null residual, `training_terms`, a log
density, a coverage, a per-row interval, an error statistic or a paired difference. Its stdout carries optimiser
fields only (status, cost, evaluations, seconds). Convergence is read from `status`, `budget_exhausted`, `cost` and
`evaluations` only, and a pair fit is converged when the optimiser's status is > 0 (sec. 3.1). The scoring step is
scripts/generators/score_thiol_lolo_nested_tau.py, run once after every pair record exists.

The frozen EX-056 generator (fit_thiol_lolo.py) is imported, never modified.

    python scripts/generators/fit_thiol_lolo_pairs.py --list                       # the pairs, no engine run
    python scripts/generators/fit_thiol_lolo_pairs.py --guard-engine               # sec. 4.2 (b) and (c)
    python scripts/generators/fit_thiol_lolo_pairs.py --pair LAB_A LAB_B --freeze-commit HASH
    python scripts/generators/fit_thiol_lolo_pairs.py --all --freeze-commit HASH   # every unfinished pair
Writes results/validation/kinetic_core_thiol_lolo_pair_<A>__<B>.json (A < B), and for the guard
kinetic_core_thiol_lolo_pair_guard_<lab>.json.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows, fitting  # noqa: E402

import fit_deoxyosone_submodel as sm3  # noqa: E402
import fit_joint_deoxyosone_thiol as fj  # noqa: E402
import fit_thiol_lolo as lolo  # noqa: E402
from fit_thiol_lolo import (ARM, BUDGET, FTOL, FULL_FIT, MAX_CONTINUATIONS, NOISE_FROM,  # noqa: E402
                            REPRODUCTION_TOL, _fit_record, converge, engine_state, lab_of, load_fit_table)

STEM = "kinetic_core_thiol_lolo_pair"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_between_lab_tau_h0_prereg.md"
#: prereg sec. 1: the five scorable laboratories. Nestle has no paired held-out row in EX-056: it stays in every
#: training set and is never held out.
SCORABLE = ("dalian_polytechnic_chen_d", "dfa_hofmann_schieberle", "jiangnan_zhang_x", "ochanomizu_murata",
            "reading_mottram")
NEVER_HELD = ("nestle_blank",)
#: prereg sec. 4.2: the guard's checks
GUARD_REPREDICT = ("reading_mottram", "dfa_hofmann_schieberle")   # (b): warm folds re-predicted from saved theta
GUARD_FOLD = "reading_mottram"                                    # (c): one single-out fold through this runner
GUARD_RAW_TOL_B = 1e-6      # dex
GUARD_RAW_TOL_C = 0.01      # dex
#: keys a pair record must never carry (prereg sec. 4.3, unit-tested)
FORBIDDEN_KEYS = frozenset({"models", "log_density", "coverage", "per_row", "paired_log_density_difference",
                            "error", "null_raw"})
#: prereg sec. 10: at most 8 cores; 2 pair fits at a time, 4 pot workers each
PARALLEL_FITS = 2
WORKERS_PER_FIT = 4


def pairs() -> List[Tuple[str, str]]:
    return [tuple(sorted(p)) for p in itertools.combinations(SCORABLE, 2)]


def split(rows: Sequence[Any], held_labs: Sequence[str]) -> Tuple[tuple, Dict[str, tuple]]:
    """(training rows, {held lab: its rows}): training = every row whose laboratory is in none of ``held_labs``;
    row order kept, as fit_thiol_lolo.make_folds keeps it."""
    held_set = set(held_labs)
    train = tuple(r for r in rows if lab_of(r) not in held_set)
    held = {lab: tuple(r for r in rows if lab_of(r) == lab) for lab in held_labs}
    return train, held


def out_path(held_labs: Sequence[str], guard: bool = False) -> Path:
    if guard:
        return data_paths.VALIDATION_DIR / f"{STEM}_guard_{'__'.join(held_labs)}.json"
    return data_paths.VALIDATION_DIR / f"{STEM}_{'__'.join(sorted(held_labs))}.json"


def row_key(r) -> Dict[str, Any]:
    return {"ledger_id": r.ledger_id, "species": r.species, "t_min": float(r.t_min), "method": fj.method_of(r)}


def _write(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(provenance.dumps_compact(payload, default=float) + "\n", encoding="utf-8")


def _provenance(name: str) -> Dict[str, Any]:
    tables = fj.TABLES + (fj.ACTZ_TABLE, fj.ENTRY_TABLE, fj.ENTRY_PH_TABLE)
    return provenance.provenance_block(
        name, generated_by="scripts/generators/fit_thiol_lolo_pairs.py",
        inputs=[fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in tables]
               + [fit_rows.LABORATORIES, ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                  ROOT / "src/kinetic_core/fitting.py", ROOT / "scripts/generators/build_submodel_probe.py",
                  ROOT / "scripts/generators/fit_joint_deoxyosone_thiol.py", ROOT / "scripts/generators/fit_thiol_lolo.py",
                  PREREG,
                  data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{FULL_FIT}.json",
                  data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{NOISE_FROM}.json"])


def converged(fit_dict: Dict[str, Any]) -> bool:
    """Prereg sec. 3.1: the optimiser's status > 0 (scipy's ftol / xtol / gtol). Status 0 is a budget exhausted
    after its continuations, status -1 a failure; neither counts."""
    return int(fit_dict["status"]) > 0


def run_pair(held_labs: Sequence[str], budget: int = BUDGET, max_continuations: int = MAX_CONTINUATIONS,
             freeze_commit: str = "", guard: bool = False) -> Dict[str, Any]:
    """A copy of fit_thiol_lolo.run_fold up to (and excluding) its scoring block, for one or two held-out
    laboratories; then the raw residuals of every held-out row and every training row, per variant."""
    t0 = time.time()
    data, jobs, pots_all = load_fit_table()
    labs_present = {lab_of(r) for r in data}
    unknown = [lab for lab in held_labs if lab not in labs_present]
    if unknown:
        raise SystemExit(f"unknown laboratory {unknown}; laboratories: {sorted(labs_present)}")
    train, held = split(data, held_labs)
    assert not ({lab_of(r) for r in train} & set(held_labs))
    assert all({lab_of(r) for r in rows} == {lab} for lab, rows in held.items())
    out = out_path(held_labs, guard=guard)
    state: Dict[str, Any] = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {}
    if state.get("complete"):
        print(f"{'+'.join(held_labs)}: complete already ({out.name})")
        return state
    keys, kinds, lower, upper, start_declared = fj.coordinates(False, fj.ARMS[ARM])
    full = _fit_record(FULL_FIT)
    theta_full = np.array([full["fit"]["theta"][k] for k in keys], dtype=float)
    ident_full = {k: bool(full["laplace"]["identified"][k]) and not bool(full["laplace"]["at_bound"][k]) for k in keys}
    fixed_noise = dict(_fit_record(NOISE_FROM)["series_extra_noise"])
    missing = {fj.series_of(r) for r in data} - set(fixed_noise)
    if missing:
        raise SystemExit(f"{NOISE_FROM}: no noise for series {sorted(missing)[:5]}")
    fit_model = lolo._model(fixed_noise)
    raw_many = fj.make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    def make_problem(rows):
        return fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=rows,
                               residuals=lambda th, rs: fit_model.whiten(raw(th, rs), rs),
                               evaluate_many=lambda pts, rs: [fit_model.whiten(v, rs) for v in raw_many(pts, rs)],
                               group_of=sm3.paper_of, nuisance_dof=fit_model.nuisance_dof)

    problem = make_problem(train)
    if "reproduction" not in state:
        z = make_problem(data).evaluate(theta_full)
        cost = float(0.5 * z @ z)
        rel = abs(cost - float(full["fit"]["cost"])) / float(full["fit"]["cost"])
        state["reproduction"] = {"cost_at_full_optimum": cost, "record_cost": float(full["fit"]["cost"]),
                                 "relative_difference": rel, "tolerance": REPRODUCTION_TOL}
        print(f"reproduction cost {cost:.3f} vs record {full['fit']['cost']:.3f} (rel {rel:.2e})", flush=True)
        if rel > REPRODUCTION_TOL:
            raise SystemExit(f"the engine does not reproduce {FULL_FIT}'s cost ({cost:.3f} vs {full['fit']['cost']:.3f})")
    state.update({"artifact": out.stem, "labs": list(held_labs), "arm": ARM, "full_fit": FULL_FIT,
                  "noise_from": NOISE_FROM, "budget": budget, "ftol": FTOL, "max_continuations": max_continuations,
                  "freeze_commit": freeze_commit, "guard": guard,
                  "training_rows": len(train), "training_laboratories": sorted({lab_of(r) for r in train}),
                  "held_out_rows": {lab: len(rows) for lab, rows in held.items()}})

    def variant(name: str, start: np.ndarray) -> Dict[str, Any]:
        rec = state.setdefault(name, {"start": dict(zip(keys, map(float, start))), "history": []})

        def save(last):
            rec.update({"fit": last, "converged": converged(last)})
            _write(out, state)   # checkpoint after every segment: a restart resumes from the last end point

        converge(problem, start, budget, max_continuations, rec["history"], keys, save)
        if not rec.get("converged"):
            print(f"{name}: not converged (status {rec['fit']['status']}) after {len(rec['history'])} segments",
                  flush=True)
        return rec

    def laplace_of(name: str) -> None:
        rec = state[name]
        if "laplace" in rec:
            return
        theta = np.array([rec["fit"]["theta"][k] for k in keys])
        rec["laplace"] = fitting.laplace(problem, theta)
        _write(out, state)

    variant("warm", theta_full)
    laplace_of("warm")
    leak = [k for k in keys if ident_full[k]
            and (not state["warm"]["laplace"]["identified"][k] or state["warm"]["laplace"]["at_bound"][k])]
    state["leak_prone"] = leak
    if leak:
        reset_start = np.array([start_declared[i] if k in leak else theta_full[i] for i, k in enumerate(keys)])
        variant("reset", reset_start)
        laplace_of("reset")
    # ---- raw residuals only (prereg sec. 4.3); no score of any kind
    thetas = {"warm": state["warm"]["fit"]["theta"], **({"reset": state["reset"]["fit"]["theta"]} if leak else {})}
    state["row_keys"] = {"held": {lab: [row_key(r) for r in rows] for lab, rows in held.items()},
                         "train": [row_key(r) for r in train]}
    state["held_engine_raw"] = {v: {lab: [float(x) for x in raw(np.array([t[k] for k in keys]), rows)]
                                    for lab, rows in held.items()} for v, t in thetas.items()}
    state["train_engine_raw"] = {v: [float(x) for x in raw(np.array([t[k] for k in keys]), train)]
                                 for v, t in thetas.items()}
    state.update({"provenance": _provenance(out.stem), "failed_pot_runs": len(fj.FAILURES),
                  "complete": True, "converged_all": all(state[v]["converged"] for v in thetas),
                  "wall_seconds_this_run": round(time.time() - t0, 1)})
    _write(out, state)
    print(f"{'+'.join(held_labs)}: warm status {state['warm']['fit']['status']} "
          f"({sum(h['evaluations'] for h in state['warm']['history'])} ev); leak-prone {len(leak)}; "
          f"converged {state['converged_all']}; {state['wall_seconds_this_run']} s", flush=True)
    return state


# ----------------------------------------------------------------------------------------------------------------
# The engine guard (prereg sec. 4.2 (b) and (c)); (a) is run_pair's reproduction check
# ----------------------------------------------------------------------------------------------------------------

def _ex056(lab: str) -> Dict[str, Any]:
    return json.loads((data_paths.VALIDATION_DIR / f"{lolo.STEM}_{lab}.json").read_text(encoding="utf-8"))


def guard_engine() -> bool:
    """(b) re-predict two saved warm folds' held-out rows from their saved theta (to 1e-6 dex); (c) re-run one
    single-out fold through this runner (cost to a relative 1e-3, held-out engine_raw to 0.01 dex). Prints the
    largest differences and pass / fail only."""
    ok = True
    data, jobs, _ = load_fit_table()
    keys = fj.coordinates(False, fj.ARMS[ARM])[0]
    raw_many = fj.make_raw_many(jobs, keys)
    for lab in GUARD_REPREDICT:
        rec = _ex056(lab)
        _, held = split(data, (lab,))
        rows = held[lab]
        assert [row_key(r) for r in rows] == [{k: x[k] for k in ("ledger_id", "species", "t_min", "method")}
                                              for x in rec["rows"]], f"{lab}: held-out rows differ from EX-056's"
        theta = np.array([rec["warm"]["fit"]["theta"][k] for k in keys], dtype=float)
        new = np.asarray(raw_many([theta], rows)[0], dtype=float)
        old = np.array([x["engine_raw"]["warm"] for x in rec["rows"]], dtype=float)
        same_nan = bool(np.array_equal(np.isfinite(new), np.isfinite(old)))
        diff = float(np.nanmax(np.abs(new - old))) if np.isfinite(new).any() else 0.0
        passed = same_nan and diff <= GUARD_RAW_TOL_B
        ok &= passed
        print(f"guard (b) {lab}: max |diff| {diff:.2e} dex, same finite mask {same_nan} -> {'PASS' if passed else 'FAIL'}",
              flush=True)
    rec = _ex056(GUARD_FOLD)
    state = run_pair((GUARD_FOLD,), guard=True)
    cost_new, cost_old = float(state["warm"]["fit"]["cost"]), float(rec["warm"]["fit"]["cost"])
    rel = abs(cost_new - cost_old) / cost_old
    new = np.array(state["held_engine_raw"]["warm"][GUARD_FOLD], dtype=float)
    old = np.array([x["engine_raw"]["warm"] for x in rec["rows"]], dtype=float)
    same_nan = bool(np.array_equal(np.isfinite(new), np.isfinite(old)))
    diff = float(np.nanmax(np.abs(new - old))) if np.isfinite(new).any() else 0.0
    passed = rel <= REPRODUCTION_TOL and same_nan and diff <= GUARD_RAW_TOL_C
    ok &= passed
    print(f"guard (c) {GUARD_FOLD}: cost rel diff {rel:.2e}, max |raw diff| {diff:.2e} dex, same finite mask "
          f"{same_nan} -> {'PASS' if passed else 'FAIL'}", flush=True)
    return ok


# ----------------------------------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------------------------------

def _check_freeze(commit: str) -> str:
    if not commit:
        raise SystemExit("--freeze-commit is required: the prereg, the controls and both scripts are committed first "
                         "(prereg sec. 7.3)")
    full = subprocess.run(["git", "rev-parse", "--verify", f"{commit}^{{commit}}"], cwd=str(ROOT),
                          capture_output=True, text=True)
    if full.returncode != 0:
        raise SystemExit(f"--freeze-commit {commit!r} is not a commit")
    full_hash = full.stdout.strip()
    if subprocess.run(["git", "merge-base", "--is-ancestor", full_hash, "HEAD"], cwd=str(ROOT)).returncode != 0:
        raise SystemExit(f"--freeze-commit {full_hash} is not an ancestor of HEAD")
    frozen = ["scripts/generators/fit_thiol_lolo_pairs.py", "scripts/generators/score_thiol_lolo_nested_tau.py",
              str(PREREG.relative_to(ROOT))]
    changed = subprocess.run(["git", "diff", "--name-only", full_hash, "--"] + frozen, cwd=str(ROOT),
                             capture_output=True, text=True).stdout.split()
    if changed:
        raise SystemExit(f"frozen files changed since {full_hash}: {changed}")
    return full_hash


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--list", action="store_true", help="print the pairs and exit (no engine run)")
    g.add_argument("--guard-engine", action="store_true", help="prereg sec. 4.2 (b) and (c)")
    g.add_argument("--pair", nargs=2, metavar="LAB", help="run one pair fit (resumes from its checkpoint)")
    g.add_argument("--all", action="store_true", help="every unfinished pair, PARALLEL_FITS at a time")
    parser.add_argument("--freeze-commit", default="")
    parser.add_argument("--budget", type=int, default=BUDGET)
    parser.add_argument("--max-continuations", type=int, default=MAX_CONTINUATIONS)
    parser.add_argument("--log-dir", default="", help="--all: one log file per pair here (else inherited stdout)")
    args = parser.parse_args(argv)
    if args.list:
        data, _, _ = load_fit_table()
        for p in pairs():
            train, held = split(data, p)
            print(f"{p[0]} + {p[1]}: train {len(train)} rows; held {', '.join(f'{k} {len(v)}' for k, v in held.items())}")
        return 0
    if args.guard_engine:
        with engine_state():
            return 0 if guard_engine() else 1
    commit = _check_freeze(args.freeze_commit)
    if args.pair:
        p = tuple(sorted(args.pair))
        if p not in pairs():
            raise SystemExit(f"{p} is not a pair of scorable laboratories {SCORABLE}")
        with engine_state():
            run_pair(p, args.budget, args.max_continuations, freeze_commit=commit)
        return 0
    todo = [p for p in pairs()
            if not (out_path(p).exists() and json.loads(out_path(p).read_text(encoding="utf-8")).get("complete"))]
    env = {**os.environ, "MAILLARD_WORKERS": str(WORKERS_PER_FIT)}
    running: List[Tuple[Tuple[str, str], subprocess.Popen]] = []
    failed = []
    while todo or running:
        while todo and len(running) < PARALLEL_FITS:
            p = todo.pop(0)
            cmd = ["nice", "-n", "10", sys.executable, str(Path(__file__).resolve()), "--pair", *p,
                   "--freeze-commit", commit, "--budget", str(args.budget),
                   "--max-continuations", str(args.max_continuations)]
            log = open(Path(args.log_dir) / f"h0_pair_{p[0]}__{p[1]}.log", "a") if args.log_dir else None
            print(f"== {p[0]} + {p[1]} ({time.strftime('%H:%M:%S')})", flush=True)
            running.append((p, subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=log,
                                                stderr=subprocess.STDOUT if log else None)))
        time.sleep(5)
        for item in list(running):
            p, proc = item
            if proc.poll() is not None:
                running.remove(item)
                print(f"-- {p[0]} + {p[1]} exited {proc.returncode} ({time.strftime('%H:%M:%S')})", flush=True)
                if proc.returncode != 0:
                    failed.append(p)
    if failed:
        print(f"pairs that exited with an error: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

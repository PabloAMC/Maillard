#!/usr/bin/env python
"""
frozen_fit_drift_gate -- every frozen fit a shipped constant rests on still has its objective, or the drift is
registered (2026-11-07, PR #107).

WHY. B1, B18 and B2.3 drifted silently for weeks: the trunk and the pH physics changed under them, nothing
re-evaluated their objectives, and 68 shipped coordinates came to rest on objectives that no longer exist
(results/validation/frozen_fit_reproducibility.md). This gate makes such a drift fail the build the day it happens.

WHAT. For each fit (build_frozen_fit_reproducibility.FIT_ORDER), ONE evaluation of the fit's own residual function at
its frozen point, through the reproducibility generator's ``--probe`` (one process per fit: several wave generators
install module state at import; nothing is duplicated here). FAIL when

1. S = r . r moved beyond the fit's tolerance (|S_today - S_recorded| > tol x max(1, S_recorded), tol the fit's
   record's, else 1 %) or its row set changed (row ids where the report keeps rows, else the count), UNLESS
   docs/reference/frozen_fit_drift_allowlist.yml names the fit with its exact current S (to its stated precision), its
   row count, its class (BENIGN / UNRESOLVED / MATERIAL) and a CF id (a PR alone suffices only for BENIGN);
2. an allow-list entry is stale: its S no longer matches (the drift moved again), or the fit reproduces again;
3. an entry is malformed: an unknown fit or class, a CF id not in docs/reference/conflicts_register.yml, or a class
   and S that results/validation/frozen_fit_reproducibility.json does not give the fit (the class is judged there, by
   the audit's Jacobian and step; the gate cannot judge it cheaply, so it insists the audit ran at this S).

Slow fits (B2.3, B60: 22-35 s per evaluation) are evaluated only with ``--slow`` (the weekly scheduled CI job);
without it their entries get the structural and record checks only.

    python scripts/ci/frozen_fit_drift_gate.py              # the cheap fits
    python scripts/ci/frozen_fit_drift_gate.py --slow       # every fit
    python scripts/ci/frozen_fit_drift_gate.py --fit B18    # one fit (repeatable)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "scripts" / "generators" / "build_frozen_fit_reproducibility.py"
ALLOWLIST = ROOT / "docs" / "reference" / "frozen_fit_drift_allowlist.yml"
REGISTER = ROOT / "docs" / "reference" / "conflicts_register.yml"
RECORD = ROOT / "results" / "validation" / "frozen_fit_reproducibility.json"
SLOW_FITS = ("B2.3", "B60")
#: allow-list class -> the reproducibility record's class
CLASSES = {"BENIGN": "DRIFTED-BENIGN", "UNRESOLVED": "DRIFTED-UNRESOLVED", "MATERIAL": "DRIFTED-MATERIAL"}
#: classes whose entry must name a conflicts-register entry (a PR alone is not a record of the drift's cause)
NEEDS_CF = ("UNRESOLVED", "MATERIAL")
#: an entry must name S at least this tightly (the loosest fit tolerance is 1 %: an entry must not approach it)
MAX_PRECISION = 1e-3
CF_ID = re.compile(r"^CF-\d{3}$")


def _generator():
    """The reproducibility generator, for its fit list and probe prefix (cheap: the specs build lazily)."""
    spec = importlib.util.spec_from_file_location("frozen_fit_reproducibility_generator", GENERATOR)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _close(a: float, b: float, rel: float) -> bool:
    return abs(float(a) - float(b)) <= rel * max(1.0, abs(float(b)))


# ---------------------------------------------------------------------------------------------------------------------
# The judgement (pure: probes, allow-list, record, register ids in; problems out). Unit-tested on toy fits.
# ---------------------------------------------------------------------------------------------------------------------

def reproduces(p: Mapping[str, Any]) -> Tuple[bool, str]:
    """Does the probe reproduce its frozen record? (the generator's REPRODUCES rule; ids where the report keeps rows)"""
    s_ok = _close(p["S_today"], p["S_recorded"], float(p["rel_tol"]))
    if p.get("rows_then") is not None:
        then, now = set(p["rows_then"]), set(p["rows_now"])
        rows_ok = then == now
        rows = "" if rows_ok else (f"rows withdrawn {sorted(then - now)[:5]}, added {sorted(now - then)[:5]}")
    else:
        rows_ok = p.get("n_rows_then") is None or int(p["n_rows_then"]) == int(p["n_rows_now"])
        rows = "" if rows_ok else f"rows {p['n_rows_then']} -> {p['n_rows_now']}"
    change = (float(p["S_today"]) - float(p["S_recorded"])) / max(1.0, abs(float(p["S_recorded"])))
    what = f"S {float(p['S_recorded']):.6g} -> {float(p['S_today']):.6g} ({change:+.2%}, tol {float(p['rel_tol']):.0%})"
    return s_ok and rows_ok, what + (f"; {rows}" if rows else "")


def allowlist_problems(entries: Sequence[Mapping[str, Any]], fits: Sequence[str], record: Mapping[str, Any],
                       register_ids: Iterable[str]) -> List[str]:
    """Structure, register and record checks of the allow-list (no evaluation)."""
    problems: List[str] = []
    register_ids = set(register_ids)
    seen = set()
    rec_fits = (record or {}).get("fits") or {}
    for e in entries:
        fit = e.get("fit")
        tag = f"allow-list {fit}"
        if fit not in fits:
            problems.append(f"{tag}: unknown fit (not in the reproducibility generator's FIT_ORDER)")
            continue
        if fit in seen:
            problems.append(f"{tag}: duplicate entry")
        seen.add(fit)
        cls = e.get("class")
        if cls not in CLASSES:
            problems.append(f"{tag}: class {cls!r} not in {sorted(CLASSES)}")
        try:
            s, prec, n = float(e["S"]), float(e["S_rel_precision"]), int(e["n_rows"])
        except (KeyError, TypeError, ValueError):
            problems.append(f"{tag}: needs a numeric S, S_rel_precision and n_rows")
            continue
        if not (0.0 < prec <= MAX_PRECISION) or not math.isfinite(s):
            problems.append(f"{tag}: S_rel_precision {prec:g} must be in (0, {MAX_PRECISION:g}] and S finite")
        cf, pr = e.get("cf"), e.get("pr")
        if not cf and not pr:
            problems.append(f"{tag}: names neither a CF id nor a PR")
        elif cls in NEEDS_CF and not cf:
            problems.append(f"{tag}: a {cls} entry must name its conflicts-register entry (cf), not only a PR")
        held = e.get("held_in")
        if held is not None:
            held_rec = rec_fits.get(held)
            if held not in fits:
                problems.append(f"{tag}: held_in {held!r} is not a fit")
            elif held_rec is None or held_rec["verdict"]["class"] != "REPRODUCES":
                problems.append(f"{tag}: held_in {held} but the reproducibility record does not class {held} "
                                "REPRODUCES")
        if cf and (not CF_ID.match(str(cf)) or cf not in register_ids):
            problems.append(f"{tag}: {cf} is not an entry of docs/reference/conflicts_register.yml")
        if pr is not None and (not isinstance(pr, int) or pr <= 0):
            problems.append(f"{tag}: pr {pr!r} is not a PR number")
        rec = rec_fits.get(fit)
        if rec is None:
            problems.append(f"{tag}: the reproducibility record has no audit of {fit} (run "
                            f"build_frozen_fit_reproducibility.py --fit {fit})")
            continue
        rec_cls = rec["verdict"]["class"]
        if cls in CLASSES and rec_cls != CLASSES[cls]:
            problems.append(f"{tag}: class {cls} but the reproducibility record classes it {rec_cls}")
        if not _close(rec["S_today"], s, prec) or int(rec["n_rows_now"]) != n:
            problems.append(f"{tag}: the record's audit ran at S {rec['S_today']:.6f} / {rec['n_rows_now']} rows, not at "
                            f"the entry's {s:.6f} / {n}: its class is not judged at this S (re-audit: "
                            f"build_frozen_fit_reproducibility.py --fit {fit})")
    for fit, rec in rec_fits.items():
        if rec["verdict"]["class"] != "REPRODUCES" and fit not in seen:
            problems.append(f"{fit}: the reproducibility record classes it {rec['verdict']['class']} and no allow-list "
                            "entry names it")
    return problems


def probe_problems(probes: Mapping[str, Mapping[str, Any]], entries: Sequence[Mapping[str, Any]]) -> List[str]:
    """The live judgement of each evaluated fit against its allow-list entry (if any)."""
    by_fit = {e.get("fit"): e for e in entries}
    problems: List[str] = []
    for fit, p in probes.items():
        ok, what = reproduces(p)
        e = by_fit.get(fit)
        if ok and e is not None:
            problems.append(f"{fit}: reproduces again ({what}): its allow-list entry is stale, remove it (and close "
                            f"{e.get('cf') or 'its record'})")
        elif not ok and e is None:
            problems.append(f"{fit}: DRIFTED, {what}. Re-derive the fit, or register the drift (a CF entry), re-audit it "
                            f"(build_frozen_fit_reproducibility.py --fit {fit}) and allow-list it at S "
                            f"{float(p['S_today']):.6f}, {p['n_rows_now']} rows")
        elif not ok:
            try:
                s, prec, n = float(e["S"]), float(e["S_rel_precision"]), int(e["n_rows"])
            except (KeyError, TypeError, ValueError):
                continue  # reported by allowlist_problems
            if not _close(p["S_today"], s, prec) or int(p["n_rows_now"]) != n:
                problems.append(f"{fit}: STALE allow-list entry: S today {float(p['S_today']):.6f} / {p['n_rows_now']} "
                                f"rows, the entry names {s:.6f} / {n} (precision {prec:g}): the drift moved again; "
                                "re-audit and re-name it")
    return problems


# ---------------------------------------------------------------------------------------------------------------------
# Evaluation (one process per fit, at most ``jobs`` at a time, one thread each)
# ---------------------------------------------------------------------------------------------------------------------

def run_probe(fit: str, prefix: str) -> Dict[str, Any]:
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env.setdefault(var, "1")
    env["MAILLARD_WORKERS"] = "1"
    out = subprocess.run([sys.executable, str(GENERATOR), "--probe", fit], cwd=ROOT, env=env, capture_output=True,
                         text=True)
    for line in out.stdout.splitlines():
        if line.startswith(prefix):
            return json.loads(line[len(prefix):])
    raise RuntimeError(f"{fit}: the probe printed no result (exit {out.returncode}):\n{out.stderr[-2000:]}")


def evaluate(fits: Sequence[str], jobs: int = 2) -> Dict[str, Dict[str, Any]]:
    prefix = _generator().PROBE_PREFIX
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        results = list(pool.map(lambda f: run_probe(f, prefix), fits))
    return {r["fit"]: r for r in results}


def _load_inputs() -> Tuple[List[dict], dict, set]:
    entries = (yaml.safe_load(ALLOWLIST.read_text(encoding="utf-8")) or {}).get("entries") or []
    record = json.loads(RECORD.read_text(encoding="utf-8")) if RECORD.exists() else {}
    register = (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []
    return entries, record, {str(e.get("id")) for e in register}


def check(slow: bool = False, only: Optional[Sequence[str]] = None, jobs: int = 2) -> Tuple[List[str], Dict[str, Any]]:
    gen = _generator()
    fits = list(gen.FIT_ORDER)
    entries, record, register_ids = _load_inputs()
    problems = allowlist_problems(entries, fits, record, register_ids)
    unknown = [f for f in (only or []) if f not in fits]
    problems += [f"{f}: unknown fit" for f in unknown]
    chosen = [f for f in (only or fits) if f in fits]
    skipped = [] if slow or only else [f for f in chosen if f in SLOW_FITS]
    chosen = [f for f in chosen if f not in skipped]
    t0 = time.time()
    probes = evaluate(chosen, jobs)
    problems += probe_problems(probes, entries)
    by_fit = {e.get("fit"): e for e in entries}
    info = {"evaluated": chosen, "skipped": skipped, "seconds": time.time() - t0,
            "allowed": [f"{f} {by_fit[f]['class']} {by_fit[f].get('cf') or 'PR #' + str(by_fit[f].get('pr'))}"
                        for f in chosen if f in by_fit],
            "probes": probes}
    return problems, info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slow", action="store_true", help="also evaluate the slow fits (B2.3, B60)")
    ap.add_argument("--fit", action="append", help="evaluate only this fit (repeatable)")
    ap.add_argument("--jobs", type=int, default=2, help="fits evaluated at once (one process each; default 2)")
    a = ap.parse_args(argv)
    problems, info = check(slow=a.slow, only=a.fit, jobs=a.jobs)
    if problems:
        print("frozen_fit_drift_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    n_ok = len(info["evaluated"]) - len(info["allowed"])
    skipped = f"; skipped without --slow: {', '.join(info['skipped'])}" if info["skipped"] else ""
    print(f"frozen_fit_drift_gate: PASS ({len(info['evaluated'])} fits evaluated in {info['seconds']:.0f} s: {n_ok} "
          f"reproduce, allow-listed drifts: {', '.join(info['allowed']) or 'none'}{skipped})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

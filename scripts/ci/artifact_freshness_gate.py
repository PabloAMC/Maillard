#!/usr/bin/env python
"""Artifact freshness gate.

2026-09-03 (improvement backlog). ``results/validation/cutover_final_exam.json`` drifted from
the code silently (4/34 -> 3/34) until B2 regenerated it; nothing had compared a tracked
artifact with what the code produces. Three checks, all cheap enough for CI:

1. **Regenerate and compare** the artifacts whose generators run in under a minute
   (the scorecard, the directional scorecard, the fit-target record): the live payload must
   equal the tracked file once volatile keys (date, git head, wall time) are stripped.
2. **``--check`` runners**: the generators that already know how to compare themselves with
   their tracked output (wave manifest, model card, data/README, results/README, the two
   key registries) must all report current.
3. **Input hashes** for every ``results/validation/*.json`` that carries a ``provenance``
   block (the envelope included, which takes 35 minutes to regenerate): each recorded input
   must still hash the same. A stale input means the artifact describes a model that no
   longer exists.

Usage:
    python scripts/ci/artifact_freshness_gate.py
    python scripts/ci/artifact_freshness_gate.py --skip-regenerate   # checks 2 and 3 only
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_access, data_paths, provenance  # noqa: E402


def _say(line: str) -> None:
    """One progress line, flushed at once. 2026-11-12: without the flush, CI's pipe block-buffers stdout, so every line
    after the core panel's reached the log at the gate's exit, and the log read as if the directional scorecard alone
    took the whole tail (97 min on a slow runner) when it also covered the fit targets, the ranking and the --check
    runners."""
    print(line, flush=True)


def _machine() -> str:
    """The runner's CPU model and count: the gate's wall time scales with them (3.6 vs 9 min for the same core panel on
    two GitHub-hosted ubuntu runners, 2026-11-11/12), so the log must say which machine it ran on."""
    model = platform.processor() or platform.machine()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return (f"  [machine] {model}; {os.cpu_count()} CPUs; MAILLARD_WORKERS={os.environ.get('MAILLARD_WORKERS', 'unset')}")

CHECK_RUNNERS = (
    "scripts/generators/build_wave_manifest.py",
    "scripts/generators/generate_model_card.py",
    "scripts/generators/build_data_readme.py",
    "scripts/generators/build_results_readme.py",
    "scripts/generators/build_paper_registry.py",
    "scripts/generators/build_compound_registry.py",
)


def _live_scorecard() -> Dict:
    from src.kinetic_core import scoring

    return scoring.score_panel()


def _live_directional() -> Dict:
    from src.kinetic_core import directional

    return directional.score_panel()


def _live_fit_targets() -> Dict:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "fit_targets.json"
        subprocess.run(
            [sys.executable, "scripts/generators/generate_core_fit_targets.py", "--wave", "b60", "--output", str(out)],
            cwd=str(ROOT), check=True, capture_output=True, text=True, timeout=600,
        )
        return json.loads(out.read_text(encoding="utf-8"))


def _live_ranking() -> Dict:
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [sys.executable, "scripts/generators/generate_experiment_value_ranking.py", "--output-dir", tmp],
            cwd=str(ROOT), check=True, capture_output=True, text=True, timeout=600,
        )
        return json.loads((Path(tmp) / "experiment_value_ranking.json").read_text(encoding="utf-8"))


def _live_wishlist() -> Dict:
    from src.data_wishlist import build

    return build()


def _live_hypotheses() -> Dict:
    from src.network_hypotheses.report import build

    return build()


REGENERATE: Dict[str, Callable[[], Dict]] = {
    "results/validation/core_panel_scores.json": _live_scorecard,
    "results/validation/core_directional_scores.json": _live_directional,
    "results/validation/kinetic_core_b60_fit_targets.json": _live_fit_targets,
    "results/validation/experiment_value_ranking.json": _live_ranking,
    "results/validation/data_wishlist.json": _live_wishlist,
    "results/validation/network_hypotheses.json": _live_hypotheses,
}


#: the live payloads ``check_regenerate`` computed, by artifact path; ``check_runners`` hands the core panel's to the
#: model card (``src.model_card.LIVE_CORE_PANEL_ENV``) so the card does not score the panel a second time
LIVE: Dict[str, Dict] = {}
CORE_PANEL = "results/validation/core_panel_scores.json"


def check_regenerate(failures: List[str]) -> None:
    for relative, live_fn in REGENERATE.items():
        tracked_path = ROOT / relative
        if not tracked_path.exists():
            failures.append(f"{relative}: tracked artifact missing")
            continue
        tracked = data_access.load_json(tracked_path)
        start = time.monotonic()
        try:
            live = live_fn()
        except Exception as exc:  # noqa: BLE001 - report, then keep checking the rest
            failures.append(f"{relative}: regeneration failed: {exc}")
            _say(f"  [regenerate] {relative}: FAILED ({time.monotonic() - start:.0f} s)")
            continue
        LIVE[relative] = live
        diffs = provenance.payload_differences(tracked, live)
        if diffs:
            failures.append(f"{relative}: tracked != regenerated:\n      " + "\n      ".join(diffs))
            _say(f"  [regenerate] {relative}: DIFFERS ({time.monotonic() - start:.0f} s)")
        else:
            _say(f"  [regenerate] {relative}: matches ({time.monotonic() - start:.0f} s)")


def check_runners(failures: List[str]) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ)
        if CORE_PANEL in LIVE:
            from src.model_card import LIVE_CORE_PANEL_ENV

            panel = Path(tmp) / "core_panel_scores.live.json"
            panel.write_text(json.dumps(LIVE[CORE_PANEL]), encoding="utf-8")
            env[LIVE_CORE_PANEL_ENV] = str(panel)
        _run_checks(failures, env)


def _run_checks(failures: List[str], env: Dict[str, str]) -> None:
    for relative in CHECK_RUNNERS:
        start = time.monotonic()
        proc = subprocess.run(
            # 2026-11-07: 1800 s, not 600: the model card's --check re-scores the panel live, and since the S3 refusal
            # every balance-solved row also runs the carried pH-balance rivals (~10 min serial). 2026-11-12: unless the
            # gate scored it already this run (``LIVE``), in which case the card reads that payload.
            [sys.executable, relative, "--check"], cwd=str(ROOT), capture_output=True, text=True, timeout=1800,
            env=env,
        )
        if proc.returncode != 0:
            failures.append(f"{relative} --check: " + (proc.stderr.strip() or proc.stdout.strip())[-500:])
            _say(f"  [check] {relative}: FAILED ({time.monotonic() - start:.0f} s)")
        else:
            _say(f"  [check] {relative}: current ({time.monotonic() - start:.0f} s)")


#: 2026-10-29 (owner: "Yes, exempt them"). A wave's record -- its ship rule, fit report, verification, identifiability
#: or discrepancy study -- states what was true at the commit that wrote it, and that commit pins it. Re-running it
#: because an engine file it hashed has since changed (switched-off code included) can only re-confirm history, at
#: hours of compute: regenerating 11 such records cost about three hours on 2026-10-29. So a FROZEN record's inputs
#: are not re-hashed; it must carry the commit it was written at. The LIVE artifacts (scorecards, the envelope, flux
#: attribution, baselines, contrasts, the choice sweep, the wishlist, figures) are still re-hashed and, where cheap,
#: regenerated: they describe the model as it is.
#: 2026-10-30: the Phase 2 sub-model records (each fit arm, each no-refit probe, each frozen hold-out prediction) are
#: records of the same kind -- run at a commit, quoted by number in their pre-registration -- and are pinned likewise.
#: 2026-10-31: so are the structural feasibility records (feasibility_*.py: sampled constants within declared bands, each
#: run read in its pre-registration, e.g. dpo partition prereg secs. 17 - 22); a 600-draw run costs hours.
#: 2026-11-04: so are the scoring records of a pre-registered run-once analysis (score_thiol_lolo_nested_tau.py, H0 /
#: EX-059: its controls and its single scoring run, prereg sec. 7.5 "run once"); re-running them is forbidden, not
#: merely costly, so they are pinned by the commit they were written at.
#: 2026-11-05: the supply-structure signature screen (supply_structure_signature_screen.py) is a pre-design feasibility
#: record of the same kind as feasibility_*: plug-in runs read in its pre-registration (M1 + M2 prereg sec. 0B), pinned
#: by the commit that wrote it.
FROZEN_RECORD_GENERATOR = re.compile(
    r"scripts/generators/(generate_(kinetic_core_b\d+\w*|env_\w+)_(ship_rule|fit|verification|identifiability|"
    r"structural_discrepancy)|fit_\w+|build_submodel_probe|predict_ph_acid_holdout|diagnose_\w+|probe_\w+|"
    r"feasibility_\w+|score_thiol_lolo_nested_tau|rescore_\w+|supply_structure_signature_screen)\.py$")
#: 2026-11-05: frozen scoring records whose generator also writes LIVE files (score_contrast_benchmark.py writes the
#: items record, re-derived and checked by the scorer itself before every use, and the named M1 scores): the v1 score
#: is the benchmark's run-once "before" record (prereg v6), pinned by the commit it was written at. And the
#: rescore_* records (held-out re-scores, each run at a commit and quoted in its prereg) are pinned by the pattern above.
#: 2026-11-11: so is decision-relevance step 0 (build_decision_relevance_step0.py, run once under its pre-registration,
#: results/validation/decision_relevance_prereg.md sec. 5.0, frozen in 9e02ce3a; read in STATUS_AND_PLAN sec. 4.1): its inputs (the contrast and panel
#: scores) moved with PR #108's S3 refusal after it was run; re-running it would replace a frozen verdict.
FROZEN_RECORD_FILES = frozenset({"results/validation/kinetic_core_contrast_benchmark_score_v1.json",
                                 "results/validation/decision_relevance_step0.json"})


#: 2026-10-30 (owner: "I am not sure regenerating the figure is needed until we have the engine finished"; strategy
#: sec. 4: no regeneration of a structurally changing model before Phase 3). The figures are rebuilt ONCE, at Phase
#: 3's single regeneration; until then their record's staleness is REPORTED, not failed, so it stays visible.
DEFERRED_TO_PHASE3 = ("results/validation/figure_inputs.json",)


def is_frozen_record(relative: str, payload: Dict) -> bool:
    generated_by = str((payload.get("provenance") or {}).get("generated_by") or "")
    return relative not in REGENERATE and (relative in FROZEN_RECORD_FILES
                                           or bool(FROZEN_RECORD_GENERATOR.search(generated_by)))


def check_input_hashes(failures: List[str]) -> None:
    frozen = 0
    for path in sorted(data_paths.VALIDATION_DIR.glob("*.json")):
        try:
            payload = data_access.load_json(path)
        except Exception:  # noqa: BLE001 - not every JSON here is an object
            continue
        if not isinstance(payload, dict) or "provenance" not in payload:
            continue
        if is_frozen_record(data_paths.rel(path), payload):
            git = payload["provenance"].get("git") or {}
            if not (git.get("head") or git.get("commit")):
                failures.append(f"{data_paths.rel(path)}: a frozen record must carry the commit it was written at")
            frozen += 1
            continue
        stale = provenance.stale_inputs(payload)
        if stale and data_paths.rel(path) in DEFERRED_TO_PHASE3:
            print(f"  [deferred] {data_paths.rel(path)}: {len(stale)} stale inputs -- rebuilt at Phase 3 (owner, 2026-10-30)")
        elif stale:
            failures.append(f"{data_paths.rel(path)}: stale inputs:\n      " + "\n      ".join(stale))
        else:
            print(f"  [inputs] {data_paths.rel(path)}: {len(payload['provenance'].get('inputs', []))} inputs current")
    print(f"  [inputs] {frozen} frozen wave records pinned by their commit (not re-hashed)")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--skip-regenerate", action="store_true", help="skip check 1 (about a minute)")
    args = parser.parse_args(argv)
    failures: List[str] = []
    _say(_machine())
    if not args.skip_regenerate:
        check_regenerate(failures)
    check_runners(failures)
    check_input_hashes(failures)
    if failures:
        print("artifact_freshness_gate: FAIL", file=sys.stderr)
        for line in failures:
            print(f"  - {line}", file=sys.stderr)
        return 1
    print("artifact_freshness_gate: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

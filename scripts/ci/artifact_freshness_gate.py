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
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_access, data_paths, provenance  # noqa: E402

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


def check_regenerate(failures: List[str]) -> None:
    for relative, live_fn in REGENERATE.items():
        tracked_path = ROOT / relative
        if not tracked_path.exists():
            failures.append(f"{relative}: tracked artifact missing")
            continue
        tracked = data_access.load_json(tracked_path)
        try:
            live = live_fn()
        except Exception as exc:  # noqa: BLE001 - report, then keep checking the rest
            failures.append(f"{relative}: regeneration failed: {exc}")
            continue
        diffs = provenance.payload_differences(tracked, live)
        if diffs:
            failures.append(f"{relative}: tracked != regenerated:\n      " + "\n      ".join(diffs))
        else:
            print(f"  [regenerate] {relative}: matches")


def check_runners(failures: List[str]) -> None:
    for relative in CHECK_RUNNERS:
        proc = subprocess.run(
            [sys.executable, relative, "--check"], cwd=str(ROOT), capture_output=True, text=True, timeout=600
        )
        if proc.returncode != 0:
            failures.append(f"{relative} --check: " + (proc.stderr.strip() or proc.stdout.strip())[-500:])
        else:
            print(f"  [check] {relative}: current")


#: 2026-10-29 (owner: "Yes, exempt them"). A wave's record -- its ship rule, fit report, verification, identifiability
#: or discrepancy study -- states what was true at the commit that wrote it, and that commit pins it. Re-running it
#: because an engine file it hashed has since changed (switched-off code included) can only re-confirm history, at
#: hours of compute: regenerating 11 such records cost about three hours on 2026-10-29. So a FROZEN record's inputs
#: are not re-hashed; it must carry the commit it was written at. The LIVE artifacts (scorecards, the envelope, flux
#: attribution, baselines, contrasts, the choice sweep, the wishlist, figures) are still re-hashed and, where cheap,
#: regenerated: they describe the model as it is.
#: 2026-10-30: the Phase 2 sub-model records (each fit arm, each no-refit probe, each frozen hold-out prediction) are
#: records of the same kind -- run at a commit, quoted by number in their pre-registration -- and are pinned likewise.
FROZEN_RECORD_GENERATOR = re.compile(
    r"scripts/generators/(generate_(kinetic_core_b\d+\w*|env_\w+)_(ship_rule|fit|verification|identifiability|"
    r"structural_discrepancy)|fit_\w+|build_submodel_probe|predict_ph_acid_holdout|diagnose_\w+)\.py$")


#: 2026-10-30 (owner: "I am not sure regenerating the figure is needed until we have the engine finished"; strategy
#: sec. 4: no regeneration of a structurally changing model before Phase 3). The figures are rebuilt ONCE, at Phase
#: 3's single regeneration; until then their record's staleness is REPORTED, not failed, so it stays visible.
DEFERRED_TO_PHASE3 = ("results/validation/figure_inputs.json",)


def is_frozen_record(relative: str, payload: Dict) -> bool:
    generated_by = str((payload.get("provenance") or {}).get("generated_by") or "")
    return relative not in REGENERATE and bool(FROZEN_RECORD_GENERATOR.search(generated_by))


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

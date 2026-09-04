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
            [sys.executable, "scripts/generators/generate_core_fit_targets.py", "--wave", "b9", "--output", str(out)],
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


REGENERATE: Dict[str, Callable[[], Dict]] = {
    "results/validation/core_panel_scores.json": _live_scorecard,
    "results/validation/core_directional_scores.json": _live_directional,
    "results/validation/kinetic_core_b9_fit_targets.json": _live_fit_targets,
    "results/validation/experiment_value_ranking.json": _live_ranking,
    "results/validation/data_wishlist.json": _live_wishlist,
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


def check_input_hashes(failures: List[str]) -> None:
    for path in sorted(data_paths.VALIDATION_DIR.glob("*.json")):
        try:
            payload = data_access.load_json(path)
        except Exception:  # noqa: BLE001 - not every JSON here is an object
            continue
        if not isinstance(payload, dict) or "provenance" not in payload:
            continue
        stale = provenance.stale_inputs(payload)
        if stale:
            failures.append(f"{data_paths.rel(path)}: stale inputs:\n      " + "\n      ".join(stale))
        else:
            print(f"  [inputs] {data_paths.rel(path)}: {len(payload['provenance'].get('inputs', []))} inputs current")


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

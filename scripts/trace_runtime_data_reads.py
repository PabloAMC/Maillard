#!/usr/bin/env python3
"""Which curated data files does the front door actually read?

Runs every verb of ``scripts/maillard.py`` that turns a formulation into an answer on the
shipped examples, each in its own interpreter under an audit hook, and records every file
opened under ``data/`` or ``results/literature/``. A trace, not a static guess: until
2026-10-27 the provenance census called any record sitting in ``data/lit/`` "consumed at
runtime", and the trace showed the front door opens none of those files.

    ./scripts/docker_maillard.sh run "python scripts/trace_runtime_data_reads.py"
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[1]
MARKER = "TRACED_DATA_READS="

#: The verbs and the shipped inputs they are traced on. ``ui`` serves these same verbs and
#: ``wishlist`` reads a tracked artifact, so neither adds a data dependency of its own.
VERBS: Dict[str, List[str]] = {
    "compare": ["compare", "docs/examples/compare_ribose_vs_glucose.yml"],
    "predict": ["predict", "docs/examples/compare_ribose_vs_glucose.yml", "--system", "a"],
    "explain": ["explain", "2-methyl-3-furanthiol"],
    "rank": ["rank"],
    "score": ["score", "docs/examples/reading_2026_ladder.yml", "--no-write"],
    "calibrate": ["calibrate", "docs/examples/reading_2026_ladder.yml", "--lab", "trace", "--no-write"],
}

_CHILD = f"""
import atexit, os, runpy, sys
root = {str(ROOT)!r}
opened = set()
def hook(event, args):
    if event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
        path = os.fsdecode(args[0])
        rel = os.path.relpath(os.path.abspath(path), root)
        if rel.startswith(("data" + os.sep, os.path.join("results", "literature") + os.sep)):
            opened.add(rel.replace(os.sep, "/"))
sys.addaudithook(hook)
atexit.register(lambda: print({MARKER!r} + repr(sorted(opened)), file=sys.stderr))
sys.argv = ["scripts/maillard.py"] + sys.argv[1:]
runpy.run_path("scripts/maillard.py", run_name="__main__")
"""


def trace_verb(argv: List[str]) -> List[str]:
    """Data files one verb opens; raises if the verb fails, so a broken run never reads as clean."""
    done = subprocess.run(
        [sys.executable, "-c", _CHILD, *argv],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
    )
    if done.returncode != 0:
        raise RuntimeError(f"maillard {' '.join(argv)} exited {done.returncode}:\n{done.stderr[-2000:]}")
    lines = [line for line in done.stderr.splitlines() if line.startswith(MARKER)]
    if not lines:
        raise RuntimeError(f"maillard {' '.join(argv)} left no trace:\n{done.stderr[-2000:]}")
    return list(eval(lines[-1][len(MARKER):]))  # noqa: S307 - our own repr of a list of str


def trace_front_door() -> Dict[str, List[str]]:
    return {verb: trace_verb(argv) for verb, argv in VERBS.items()}


def main() -> int:
    traced = trace_front_door()
    union = sorted({path for paths in traced.values() for path in paths})
    print(json.dumps({"by_verb": traced, "all": union}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

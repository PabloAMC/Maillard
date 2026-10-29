#!/usr/bin/env python
"""
evidence_anchor_gate -- every rate constant's evidence anchor resolves to a file (2026-10-29, guard 4).

For every ``KineticParameter`` of the Maillard lanes (trunk, sulfur, acrylamide): every file token in its
``dossier_anchor`` resolves (``src/evidence_anchors.py``), and a constant whose ``evidence_class`` says MEASURED cites
at least one extraction dossier. Structural and offline; it does not judge whether the source supports the number
(the source-verification record and the reaction evidence table do that).

    python scripts/ci/evidence_anchor_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import evidence_anchors  # noqa: E402
from src.kinetic_core import engine  # noqa: E402


def check() -> tuple:
    problems, seen = [], {}
    for lane in engine.MAILLARD_LANES:
        for key, p in engine.core_parameters(lane).items():
            if key in seen or not hasattr(p, "dossier_anchor"):
                continue
            seen[key] = p
            tokens = evidence_anchors.anchor_tokens(p.dossier_anchor)
            if not tokens:
                problems.append(f"{key}: its dossier_anchor names no file")
            for token in tokens:
                if evidence_anchors.resolve(token) is None:
                    problems.append(f"{key}: anchor {token!r} does not resolve to a file")
            if p.evidence_class in evidence_anchors.MEASURED_CLASSES and not evidence_anchors.dossier_stems(p.dossier_anchor):
                problems.append(f"{key}: evidence_class {p.evidence_class} but no extraction dossier is cited")
    return problems, len(seen)


def main() -> int:
    problems, n = check()
    if problems:
        print("evidence_anchor_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    print(f"evidence_anchor_gate: PASS ({n} constants; every anchor resolves)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

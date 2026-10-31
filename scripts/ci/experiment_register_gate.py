#!/usr/bin/env python
"""
experiment_register_gate -- every experiment is in the register, and the register is well-formed (2026-10-31).

Owner: "are we logging the experiments we did, the hypothesis, results and conclusions, in case we want to revisit
them in the future?" The register is docs/reference/experiment_register.yml. This gate fails when
  * an entry lacks a field, repeats an id, or uses a kind / verdict outside the vocabulary;
  * an entry names an artifact that does not exist;
  * a fit, signature, flux-diagnosis or probe record in results/validation (COVERED below) appears in no entry --
    a run nobody logged is a result nobody can find again.
Stdlib + PyYAML only, like the other data-integrity gates.

    python scripts/ci/experiment_register_gate.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
REGISTER = ROOT / "docs" / "reference" / "experiment_register.yml"
FIELDS = ("id", "date", "kind", "question", "hypothesis", "design", "expected", "artifacts", "commit", "verdict",
          "conclusion", "follow_up")
KINDS = {"fit", "probe", "diagnosis", "literature", "review"}
VERDICTS = {"supported", "rejected", "partial", "inconclusive", "superseded", "pending"}
#: the records that must be registered (the Phase 2 experiment records; older waves are in the plan's sec. 7)
COVERED = ("kinetic_core_joint_deoxyosone_thiol_fit_*.md", "kinetic_core_ph_acid_submodel_fit_*.md",
           "kinetic_core_deoxyosone_submodel_fit_*.md", "kinetic_core_mechanism_signatures_*.md",
           "kinetic_core_thiol_flux_diagnosis_*.md", "kinetic_core_pentose_amine_entry_probe.md")


def main() -> int:
    problems = []
    entries = (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []
    seen, registered = set(), set()
    for e in entries:
        key = e.get("id", "?")
        missing = [f for f in FIELDS if f not in e]
        if missing:
            problems.append(f"{key}: missing fields {missing}")
        if key in seen:
            problems.append(f"{key}: duplicate id")
        seen.add(key)
        if e.get("kind") not in KINDS:
            problems.append(f"{key}: kind {e.get('kind')!r} not in {sorted(KINDS)}")
        if e.get("verdict") not in VERDICTS:
            problems.append(f"{key}: verdict {e.get('verdict')!r} not in {sorted(VERDICTS)}")
        if e.get("verdict") != "pending" and e.get("kind") in ("fit", "probe", "diagnosis") and not e.get("artifacts"):
            problems.append(f"{key}: a {e.get('kind')} with a verdict names no artifact")
        for a in e.get("artifacts") or []:
            if not (ROOT / a).exists():
                problems.append(f"{key}: artifact {a} does not exist")
            registered.add(a)
    validation = ROOT / "results" / "validation"
    for pattern in COVERED:
        for path in sorted(validation.glob(pattern)):
            rel = str(path.relative_to(ROOT))
            if rel not in registered:
                problems.append(f"{rel}: an experiment record in no register entry")
    if problems:
        print("experiment_register_gate: FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"experiment_register_gate: PASS ({len(entries)} entries, {len(registered)} artifacts registered)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

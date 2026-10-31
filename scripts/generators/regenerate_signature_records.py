#!/usr/bin/env python
"""
regenerate_signature_records.py -- re-score every mechanism / deoxyosone signature record whose engine inputs changed.

2026-10-31: an engine edit (even one that is bit-identical with its new switch off) moves the sha256 of
src/kinetic_core/sulfur.py or parameters_sulfur.py, and every signature record that lists them as inputs goes stale
for the artifact-freshness gate. The scorers are deterministic, so re-scoring reproduces the same numbers with the new
hashes. This walks results/validation, finds the signature and chance records whose recorded input hashes differ from
the files on disk, and re-runs their generator with the arguments recoverable from the record:

  kinetic_core_mechanism_signatures_<fit>.json       -> score_mechanism_signatures.py --fit <fit>
  kinetic_core_mechanism_signatures_chance_<arm>      -> score_mechanism_signatures.py --chance <draws> --arm <ARM> --seed <seed>
  kinetic_core_deoxyosone_signatures_<fit>.json       -> score_deoxyosone_signatures.py --fit <fit>
  kinetic_core_deoxyosone_signatures_chance_<arm>     -> score_deoxyosone_signatures.py --chance <draws> --arm <ARM> --seed <seed>

Feasibility records and other generators are reported, not re-run (their runs are long or take arguments the record
does not hold).

    python scripts/generators/regenerate_signature_records.py [--dry-run]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VALIDATION = ROOT / "results" / "validation"


def _stale(record: dict) -> list:
    out = []
    for item in (record.get("provenance") or {}).get("inputs") or []:
        path, sha = item.get("path"), item.get("sha256")
        if not path or not sha:
            continue
        f = ROOT / path
        if f.is_file() and hashlib.sha256(f.read_bytes()).hexdigest() != sha:
            out.append(path)
    return out


def _command(path: Path, record: dict):
    stem = path.stem
    for prefix, script in (("kinetic_core_mechanism_signatures_", "score_mechanism_signatures.py"),
                           ("kinetic_core_deoxyosone_signatures_", "score_deoxyosone_signatures.py")):
        if not stem.startswith(prefix):
            continue
        rest = stem[len(prefix):]
        gen = ["python", f"scripts/generators/{script}"]
        if rest.startswith("chance_"):
            return gen + ["--chance", str(record.get("draws", 200)), "--arm", str(record["arm"]),
                          "--seed", str(record.get("seed", 20261031))]
        return gen + ["--fit", rest]
    return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    todo, other = [], []
    # frozen wave records (fits, probes, hold-out predictions, diagnoses) are pinned by their commit in the freshness
    # gate, not re-hashed: they are not reported
    frozen = ("_fit_", "_probe", "_holdout_predictions_", "_flux_diagnosis_", "_ship_rule", "_verification",
              "_feasibility_", "figure_inputs")   # the last: deferred to Phase 3 by the gate itself
    for path in sorted(VALIDATION.glob("*.json")):
        if any(tag in path.name for tag in frozen):
            continue
        try:
            record = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        stale = _stale(record) if isinstance(record, dict) else []
        if not stale:
            continue
        cmd = _command(path, record)
        (todo if cmd else other).append((path, stale, cmd))
    for path, stale, cmd in todo:
        print(f"re-score {path.name} (stale: {', '.join(stale)})")
        if not args.dry_run:
            subprocess.run(cmd, cwd=ROOT, check=True)
    for path, stale, _ in other:
        print(f"NOT re-run (another generator): {path.name} (stale: {', '.join(stale)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

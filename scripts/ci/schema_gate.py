#!/usr/bin/env python
"""Benchmark schema gate.

Every JSON under ``data/benchmarks/`` (all tiers, recursively) must:

1. validate against ``data/schemas/benchmark.schema.json`` (jsonschema, draft 2020-12);
2. carry a ``benchmark_id`` equal to its file stem (the id IS the filename everywhere
   in the code base);
3. name its volatiles with spellings that resolve in ``data/keys/compounds.yml`` --
   a compound that resolves nowhere is a typo the matcher would have silently
   fuzzy-matched at ratio 0.75 before 2026-09-01;
4. if it lives under ``external_validation/`` (either level), declare
   ``evidence_class: external_validation_only`` (scripts/ci/holdout_guard.py checks
   the same thing from the calibration side; this is the schema side).

Until 2026-09-01 there was no schema: four top-level variants and nine distinct
``measured_volatiles`` entry shapes had accumulated, and the only validator in the
repository covered intake YAML, not benchmarks. Dependencies: pyyaml + jsonschema
(both pure Python; no science stack), plus ``src.compound_keys`` for check 3.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import jsonschema  # noqa: E402

from src import compound_keys, data_paths  # noqa: E402
from src.kinetic_core.panel import quantification_family  # noqa: E402


def main() -> int:
    schema = json.loads(data_paths.BENCHMARK_SCHEMA.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema)
    files = sorted(data_paths.BENCHMARKS_DIR.rglob("*.json"))
    violations: list[str] = []
    for path in files:
        rel = data_paths.rel(path)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except ValueError as exc:
            violations.append(f"{rel}: not valid JSON: {exc}")
            continue
        for err in sorted(validator.iter_errors(payload), key=lambda e: list(e.path)):
            where = "/".join(str(p) for p in err.path) or "<root>"
            violations.append(f"{rel}: {where}: {err.message[:160]}")
        if payload.get("benchmark_id") != path.stem:
            violations.append(f"{rel}: benchmark_id {payload.get('benchmark_id')!r} != file stem {path.stem!r}")
        for block in ("measured_volatiles", "reference_volatiles", "holdout_targets"):
            for name in (payload.get(block) or {}):
                if compound_keys.resolve(name) is None:
                    violations.append(f"{rel}: {block}/{name!r} resolves to no entry in {data_paths.rel(data_paths.COMPOUND_REGISTRY)}")
        if data_paths.EXTERNAL_VALIDATION_DIR in path.parents and payload.get("evidence_class") != "external_validation_only":
            violations.append(f"{rel}: hold-out bundle without evidence_class=external_validation_only")
        # 2026-09-03: every PANEL bundle says how its number was quantified -- a class the envelope
        # can place (anywhere in the bundle, as it searches) or an explicit top-level 'undeclared'
        # with a quantification_note. The headspace-band gating must never default silently.
        if path.parent in data_paths.PANEL_DIRS:
            family, _why = quantification_family(payload)
            if family == "undeclared" and payload.get("quantification_class") != "undeclared":
                violations.append(f"{rel}: panel bundle whose quantification class resolves to no family (declare it, or set quantification_class: 'undeclared' with a quantification_note)")
    if violations:
        print(f"schema_gate: FAIL -- {len(violations)} violation(s) across {len(files)} benchmark files", file=sys.stderr)
        for v in violations:
            print(f"  {v}", file=sys.stderr)
        return 1
    print(f"schema_gate: PASS ({len(files)} benchmark files validate against {data_paths.rel(data_paths.BENCHMARK_SCHEMA)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

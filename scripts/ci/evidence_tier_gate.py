#!/usr/bin/env python
"""
evidence_tier_gate -- every shipped coordinate has an evidence tier, and no tier rises unreviewed (plan sec. 7, S1).

Sec. 7 S1 (owner, 2026-11-06): "A CI gate fails on a shipped coordinate without a tier; the record is complete (the
gate green) before any S2 / S4 application." "Upgrades need an artifact and a second reader; downgrades need a
citation." The record is ``results/validation/incumbent_evidence_tiers.json`` (``build_evidence_tiers.py`` from
``src/evidence_tiers.py``). FAILS when:

1. a shipped coordinate (log10 k_ref and Ea of every constant the engine reads, enumerated by
   ``evidence_tiers.shipped_coordinate_ids``, the generator's own enumeration) has no row, or any row has no tier from
   ``evidence_tiers.TIERS`` (``unclassified`` included);
2. the record is stale: an input recorded in its provenance block no longer hashes the same (the freshness gate's
   test), or a register entry names no row;
3. a row is ``fitted_identified`` without its S1 profile (``identification.profile_artifact``: a file under
   results/validation/ that names the constant; every other unidentified coordinate profiled over its band);
4. a row's tier differs from the frozen baseline (``incumbent_evidence_tiers_baseline.json``) carried through its
   ``tier_changes`` (``docs/reference/evidence_tier_changes.yml``), or a change is unsupported: an UPGRADE (into
   measured / analogue / bounded / fitted_identified, ``evidence_tiers.is_upgrade``) without an ``artifact`` that
   exists and a ``second_reader``; any change without a ``citation``.

``--record-only`` skips the engine enumeration in 1 (it needs numpy): the dependency-light CI job runs that, the conda
job runs the full gate.

    python scripts/ci/evidence_tier_gate.py [--record-only]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, List, Mapping, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, evidence_tiers, provenance  # noqa: E402

RECORD = data_paths.VALIDATION_DIR / "incumbent_evidence_tiers.json"
BASELINE = data_paths.VALIDATION_DIR / "incumbent_evidence_tiers_baseline.json"
VALIDATION_DIR = data_paths.VALIDATION_DIR


def _validation_file(path: object, root: Path) -> Optional[Path]:
    """The named file when it exists under results/validation/ (relative to the repo root), else None."""
    if not isinstance(path, str) or not path.strip():
        return None
    p = (root / path).resolve()
    try:
        p.relative_to(VALIDATION_DIR.resolve())
    except ValueError:
        return None
    return p if p.is_file() else None


def check_tiered(rows: List[dict], shipped: Optional[Iterable[str]]) -> List[str]:
    """Failure mode 1: a shipped coordinate without a row; a row without a tier; a row twice."""
    problems = []
    ids = Counter(r.get("id") for r in rows)
    problems += [f"{i}: {n} rows" for i, n in sorted(ids.items(), key=lambda x: str(x[0])) if n > 1]
    for r in rows:
        if r.get("tier") not in evidence_tiers.TIERS:
            problems.append(f"{r.get('id')}: no tier (tier {r.get('tier')!r}; basis: {str(r.get('basis'))[:120]})")
    for cid in shipped or ():
        if cid not in ids:
            problems.append(f"{cid}: a shipped coordinate with no row in the record (regenerate it)")
    return problems


def check_fresh(payload: Mapping) -> List[str]:
    """Failure mode 2: the record no longer describes its inputs."""
    if not (payload.get("provenance") or {}).get("inputs"):
        return ["the record carries no provenance inputs; it cannot be checked for staleness"]
    problems = [f"stale input {s} (regenerate: scripts/generators/build_evidence_tiers.py)"
                for s in provenance.stale_inputs(payload)]
    problems += [f"{e.get('coordinate')}: a tier change names no row of the record"
                 for e in payload.get("tier_changes_without_a_row") or []]
    return problems


def check_profiles(rows: List[dict], root: Path = ROOT) -> List[str]:
    """Failure mode 3: fitted_identified without the S1 profile artifact."""
    problems = []
    for r in rows:
        if r.get("tier") != "fitted_identified":
            continue
        path = (r.get("identification") or {}).get("profile_artifact")
        f = _validation_file(path, root)
        if f is None:
            problems.append(f"{r.get('id')}: fitted_identified without its S1 profile (identification.profile_artifact "
                            f"{path!r} is not a file under results/validation/)")
        elif str(r.get("constant") or r.get("id")) not in f.read_text(encoding="utf-8", errors="replace"):
            problems.append(f"{r.get('id')}: its profile artifact {path} does not name {r.get('constant')}")
    return problems


def check_changes(rows: List[dict], baseline: Mapping[str, str], root: Path = ROOT) -> Tuple[List[str], Counter]:
    """Failure mode 4: a tier that differs from the baseline without a supported chain of changes."""
    problems: List[str] = []
    counts: Counter = Counter()
    for r in rows:
        rid, tier = r.get("id"), r.get("tier")
        current = baseline.get(rid)
        for e in r.get("tier_changes") or []:
            before, after = e.get("from"), e.get("to")
            where = f"{rid}: change {before} -> {after}"
            if before != current:
                problems.append(f"{where}: its 'from' is not the tier before it ({current})")
            if after not in evidence_tiers.TIERS or after == before:
                problems.append(f"{where}: 'to' is not a different tier")
            if not str(e.get("citation") or "").strip():
                problems.append(f"{where}: no citation")
            if evidence_tiers.is_upgrade(before, after):
                counts["upgrades"] += 1
                if not str(e.get("second_reader") or "").strip():
                    problems.append(f"{where}: an upgrade without a second_reader")
                if _validation_file(e.get("artifact"), root) is None:
                    problems.append(f"{where}: an upgrade without an artifact under results/validation/ "
                                    f"(artifact {e.get('artifact')!r})")
            else:
                counts["downgrades"] += 1
            current = after
        if tier != current and not (current is None and not evidence_tiers.is_upgrade(None, tier)):
            kind = "an upgrade" if evidence_tiers.is_upgrade(current, tier) else "a downgrade"
            problems.append(f"{rid}: tier {current} -> {tier} is {kind} with no entry in "
                            f"docs/reference/evidence_tier_changes.yml")
    return problems, counts


def check(record_path: Path = RECORD, baseline_path: Path = BASELINE, shipped: Optional[Iterable[str]] = None,
          root: Path = ROOT) -> Tuple[List[str], dict]:
    payload = json.loads(record_path.read_text(encoding="utf-8"))
    rows = payload.get("rows") or []
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))["tiers"]
    shipped = None if shipped is None else list(shipped)
    problems = check_tiered(rows, shipped) + check_fresh(payload) + check_profiles(rows, root)
    change_problems, changes = check_changes(rows, baseline, root)
    problems += change_problems
    stats = {"rows": len(rows), "shipped": None if shipped is None else len(shipped),
             "by_tier": Counter(r.get("tier") for r in rows),
             "inputs": len((payload.get("provenance") or {}).get("inputs") or []),
             "baseline": len(baseline), "upgrades": changes["upgrades"], "downgrades": changes["downgrades"]}
    return problems, stats


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--record-only", action="store_true",
                        help="skip the engine's enumeration of shipped coordinates (no numpy needed)")
    args = parser.parse_args(argv)
    shipped = None if args.record_only else evidence_tiers.shipped_coordinate_ids()
    problems, s = check(shipped=shipped)
    if problems:
        print("evidence_tier_gate: FAIL")
        for p in problems:
            print("  - " + p)
        return 1
    tiers = ", ".join(f"{t} {s['by_tier'][t]}" for t in evidence_tiers.TIERS if s["by_tier"].get(t))
    covered = ("engine enumeration skipped (--record-only)" if s["shipped"] is None
               else f"all {s['shipped']} shipped coordinates tiered")
    print(f"evidence_tier_gate: PASS ({covered}; {s['rows']} rows: {tiers}; {s['inputs']} inputs current; "
          f"{s['by_tier'].get('fitted_identified', 0)} fitted_identified, each with its profile; "
          f"{s['upgrades']} upgrades and {s['downgrades']} downgrades since the {s['baseline']}-row baseline, "
          f"each reviewed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

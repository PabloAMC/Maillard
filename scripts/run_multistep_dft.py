#!/usr/bin/env python
"""Run a multi-step (elementary-step decomposed) DFT campaign.

Targets are defined in ``data/lit/computational_gap_multistep_targets.json``.
Each target lists ordered elementary steps. For each step we:

  1. Build hydrated R/P endpoints (microsolvation: explicit waters as
     proton shuttles).
  2. Run DE-GSM (pyGSM) to obtain a TS guess.
  3. Validate xTB n_imag and TS mode-coordinate concentration.
  4. Compute the barrier via DFTRefiner.calculate_robust_barrier
     (xTB-Sella TS refine + r2SCAN-3c single-points, ddCOSMO water).
  5. Aggregate (default: max_step) into a target-level barrier and
     compare against literature.

Outputs are written under
``results/computational_gap_refinement/multistep/<target_id>/``.

Usage
-----
  python scripts/run_multistep_dft.py --target pe_amadori_decomposed --execute
  python scripts/run_multistep_dft.py --list
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

# Make the src/ package importable when running from the repo root
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.dft_refiner import DFTRefiner  # noqa: E402
from src.elementary_step_runner import (  # noqa: E402
    ElementaryStepRunner,
    find_target,
    load_multistep_targets,
)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target", type=str, help="multi-step target_id to run")
    p.add_argument("--list", action="store_true", help="list available multi-step targets and exit")
    p.add_argument("--execute", action="store_true", help="actually run the campaign (otherwise dry-run summary only)")
    p.add_argument("--registry", type=Path, default=None, help="alternate multi-step registry path")
    p.add_argument(
        "--output-root", type=Path,
        default=REPO_ROOT / "results" / "computational_gap_refinement" / "multistep",
        help="root directory for per-target outputs",
    )
    p.add_argument("--gsm-nodes", type=int, default=11)
    p.add_argument("--gsm-max-iters", type=int, default=80)
    p.add_argument("--gsm-timeout-s", type=float, default=1800.0)
    p.add_argument("--ts-mode-threshold", type=float, default=0.3)
    p.add_argument("--log-level", type=str, default="INFO")
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

    payload = load_multistep_targets(args.registry)

    if args.list or not args.target:
        print(f"Multi-step targets in {args.registry or 'default registry'}:")
        for t in payload.get("targets", []):
            steps = t.get("elementary_steps", [])
            lit = (t.get("literature") or {}).get("barrier_kcal_mol")
            print(f"  - {t['target_id']:<35s} steps={len(steps):d}  lit={lit}  agg={t.get('aggregation', 'max_step')}")
        return 0

    target = find_target(payload, args.target)

    if not args.execute:
        print(f"[dry-run] target={args.target} aggregation={target.get('aggregation')}")
        for s in target.get("elementary_steps", []):
            print(f"  step={s['step_id']:<25s} R={s['dry_reactant_xyz']}  P={s['dry_product_xyz']}")
            ms = s.get("microsolvation") or {}
            print(f"    microsolvation: {ms}")
        print("Pass --execute to run.")
        return 0

    output_dir = args.output_root / args.target
    output_dir.mkdir(parents=True, exist_ok=True)

    refiner = DFTRefiner()
    runner = ElementaryStepRunner(
        refiner,
        gsm_nodes=args.gsm_nodes,
        gsm_max_iters=args.gsm_max_iters,
        gsm_timeout_s=args.gsm_timeout_s,
        ts_mode_threshold=args.ts_mode_threshold,
    )
    result = runner.run_target(target, output_dir)

    summary_path = output_dir / "target_result.json"
    print(json.dumps(result.to_dict(), indent=2))
    print(f"\nSummary written to: {summary_path}")
    return 0 if result.status == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())

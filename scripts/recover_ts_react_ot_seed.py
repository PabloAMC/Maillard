"""scripts/recover_ts_react_ot_seed.py

Pilot wrapper: run React-OT over the CHON-eligible computational-gap
targets to generate one TS-guess geometry per target.

Writes per target:
  results/computational_gap_refinement/<target>_react_ot_seed.xyz
  results/computational_gap_refinement/<target>_react_ot_seed.json
and a manifest:
  results/computational_gap_refinement/react_ot_pilot_manifest.json

Strict scope:

* React-OT is treated as a *geometric seed generator only*. Energies are
  never propagated as runtime barrier authority.
* Sulfur- and phosphorus-containing targets are explicitly excluded
  (CHON applicability scope; see
  ``models/external/react_ot/provenance.json``).
* Mode inspection and downstream Sella DFT are *not* run from this
  script. The operator decides which seeds to escalate after reading the
  per-target JSON.
* The smoke gate ``scripts/run_react_ot_smoke.py`` must pass before this
  pilot is allowed to consume real targets.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.react_ot_backend import (  # noqa: E402
    predict_ts_from_reactant_product,
    probe_backend,
)


ELIGIBLE_TARGETS: Tuple[str, ...] = (
    "lysinoalanine_crosslink",
    "aa_ring_open_dicarbonyl",
    "pe_schiff_base",
    "asparagine_sugar_explicit_water_cluster",
)


def _read_xyz_pair(target: str) -> Tuple[str, str]:
    base = REPO_ROOT / "data" / "geometries" / "xtb_inputs" / target
    return (base / "reactant.xyz").read_text(), (base / "product.xyz").read_text()


def _resolve_checkpoint(arg: str | None) -> str | None:
    if arg:
        return arg
    default = REPO_ROOT / "models" / "external" / "react_ot" / "sb-pretrained.ckpt"
    if default.exists():
        return str(default)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        action="append",
        choices=ELIGIBLE_TARGETS,
        help=(
            "Eligible CHON target (repeat for multiple). Defaults to the "
            "full eligible set listed in this module."
        ),
    )
    parser.add_argument("--checkpoint", default=None)
    parser.add_argument("--solver", default="ode")
    parser.add_argument("--nfe", type=int, default=10)
    parser.add_argument("--n-candidates", type=int, default=1)
    parser.add_argument(
        "--out-dir",
        default=str(REPO_ROOT / "results" / "computational_gap_refinement"),
    )
    args = parser.parse_args()

    targets: List[str] = args.target or list(ELIGIBLE_TARGETS)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    available, reason = probe_backend(
        "deepprinciple_react_ot_sb_pretrained", "python:src.react_ot_backend"
    )
    if not available:
        manifest = {
            "status": "probe_failed",
            "reason": reason,
            "targets": targets,
        }
        (out_dir / "react_ot_pilot_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True)
        )
        print(f"React-OT probe failed: {reason}", file=sys.stderr)
        return 1

    checkpoint = _resolve_checkpoint(args.checkpoint)
    summaries: List[dict] = []
    for target in targets:
        try:
            reactant_xyz, product_xyz = _read_xyz_pair(target)
        except FileNotFoundError as exc:
            summaries.append(
                {
                    "target": target,
                    "status": "input_xyz_missing",
                    "error": str(exc),
                }
            )
            continue

        result = predict_ts_from_reactant_product(
            reactant_xyz=reactant_xyz,
            product_xyz=product_xyz,
            checkpoint=checkpoint,
            solver=args.solver,
            nfe=args.nfe,
            n_candidates=args.n_candidates,
        )
        ts_xyz = result.pop("xyz", None)
        ts_xyz_path: str | None = None
        if ts_xyz:
            ts_xyz_path_obj = out_dir / f"{target}_react_ot_seed.xyz"
            ts_xyz_path_obj.write_text(ts_xyz)
            ts_xyz_path = str(ts_xyz_path_obj.resolve())

        summary_payload = {
            "target": target,
            "checkpoint": checkpoint,
            "solver": args.solver,
            "nfe": args.nfe,
            "n_candidates": args.n_candidates,
            "result": result,
            "ts_xyz_path": ts_xyz_path,
        }
        (out_dir / f"{target}_react_ot_seed.json").write_text(
            json.dumps(summary_payload, indent=2, sort_keys=True)
        )
        summaries.append({"target": target, "status": result.get("status", "unknown")})

    manifest = {
        "status": "completed",
        "probe_reason": reason,
        "checkpoint": checkpoint,
        "targets": summaries,
    }
    (out_dir / "react_ot_pilot_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

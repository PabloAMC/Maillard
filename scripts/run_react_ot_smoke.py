"""scripts/run_react_ot_smoke.py

Smoke gate for the React-OT secondary env.

Runs from the main ``maillard`` env. Verifies in order:

  (1) the secondary conda env exists,
  (2) ``import reactot`` succeeds inside it,
  (3) the upstream checkpoint is reachable,
  (4) the inference script can be invoked end-to-end on one CHON pair
      and produces a TS XYZ.

Writes ``results/validation/react_ot_smoke.json`` with the result. This
gate must pass before ``scripts/recover_ts_react_ot_seed.py`` is allowed
to consume real targets. If this gate fails for ``torch_geometric`` or
device reasons, fall back to ``notebooks/react_ot_colab_gpu.ipynb`` on
Google Colab.

This script does *not* propagate React-OT energies anywhere. It only
proves the pipeline is wired correctly end-to-end.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.react_ot_backend import (  # noqa: E402  (after sys.path tweak)
    DEFAULT_CHECKPOINT_ENV,
    predict_ts_from_reactant_product,
    probe_backend,
)


def _read_xyz_pair(target_dir: Path) -> tuple[str, str]:
    return (
        (target_dir / "reactant.xyz").read_text(),
        (target_dir / "product.xyz").read_text(),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        default="lysinoalanine_crosslink",
        help=(
            "Eligible CHON target under data/geometries/xtb_inputs/<target>/. "
            "Defaults to lysinoalanine_crosslink (smallest CHON pair already "
            "available in repo)."
        ),
    )
    parser.add_argument(
        "--checkpoint",
        default=None,
        help=(
            "Override the React-OT checkpoint path. Defaults to "
            f"environment variable {DEFAULT_CHECKPOINT_ENV} or "
            "models/external/react_ot/sb-pretrained.ckpt."
        ),
    )
    parser.add_argument(
        "--out",
        default=str(REPO_ROOT / "results" / "validation" / "react_ot_smoke.json"),
    )
    parser.add_argument("--solver", default="ode")
    parser.add_argument("--nfe", type=int, default=10)
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload: dict = {
        "target": args.target,
        "checkpoint": args.checkpoint,
        "solver": args.solver,
        "nfe": args.nfe,
    }

    available, reason = probe_backend(
        "deepprinciple_react_ot_sb_pretrained", "python:src.react_ot_backend"
    )
    payload["probe_backend"] = {"available": available, "reason": reason}
    if not available:
        payload["status"] = "probe_failed"
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    target_dir = REPO_ROOT / "data" / "geometries" / "xtb_inputs" / args.target
    if not target_dir.exists():
        payload["status"] = "target_not_found"
        payload["error"] = f"missing {target_dir}"
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    reactant_xyz, product_xyz = _read_xyz_pair(target_dir)

    default_checkpoint = (
        REPO_ROOT / "models" / "external" / "react_ot" / "sb-pretrained.ckpt"
    )
    checkpoint_arg = args.checkpoint
    if checkpoint_arg is None and default_checkpoint.exists():
        checkpoint_arg = str(default_checkpoint)

    inference = predict_ts_from_reactant_product(
        reactant_xyz=reactant_xyz,
        product_xyz=product_xyz,
        checkpoint=checkpoint_arg,
        solver=args.solver,
        nfe=args.nfe,
    )
    payload["inference"] = {k: v for k, v in inference.items() if k != "xyz"}
    payload["status"] = inference.get("status", "unknown")
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

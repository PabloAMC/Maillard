"""scripts/_react_ot_inference.py

Inference entry-point invoked **inside the secondary `react_ot` conda env**.

This script imports `reactot` / `torch` / `torch_geometric` directly and is
intentionally not importable from the main `maillard` env. The main env
reaches it via `subprocess` from `src.react_ot_backend`.

Reads reactant + product XYZ from disk, runs React-OT to generate one TS
geometry, and writes:
  --out-xyz       : the generated TS XYZ string
  --out-summary   : a JSON with status, runtime metadata, and any error

Exit codes: 0 ok, 2 torch missing, 3 react-ot import failed,
4 inference failed, 5 unexpected result type.

NOTE: The exact upstream programmatic API of React-OT is verified at install
time (see `models/external/react_ot/provenance.json`). The import below is
kept lazy and best-effort; if the entry point moves between releases, fix it
here, not in the main env.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace


def _read_xyz(path: str) -> str:
    return Path(path).read_text()


def _write_summary(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def _resolve_pred_ts():
    """Locate React-OT's supported batch/file inference entry point.

    Upstream exposes `reactot.run_model.pred_ts` and threads solver / ODE
    settings through a loose `opt` namespace. We resolve that directly so the
    helper matches the published service and CLI path.
    """
    try:
        from reactot.run_model import pred_ts
    except Exception as exc:  # noqa: BLE001
        raise ImportError(
            "Could not import reactot.run_model.pred_ts"
        ) from exc
    return pred_ts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reactant", required=True)
    parser.add_argument("--product", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--out-xyz", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--solver", default="ode")
    parser.add_argument("--nfe", type=int, default=10)
    parser.add_argument("--n-candidates", type=int, default=1)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=72)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--diz", default="linear")
    parser.add_argument("--method", default="midpoint")
    parser.add_argument("--atol", type=float, default=1e-2)
    parser.add_argument("--rtol", type=float, default=1e-2)
    args = parser.parse_args()

    out_summary = Path(args.out_summary)
    summary = {
        "status": "unknown",
        "reactant": str(Path(args.reactant).resolve()),
        "product": str(Path(args.product).resolve()),
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "solver": args.solver,
        "nfe": args.nfe,
        "n_candidates": args.n_candidates,
        "device": args.device,
        "batch_size": args.batch_size,
        "method": args.method,
        "order": args.order,
    }

    try:
        import torch  # noqa: F401
    except ImportError as exc:
        summary["status"] = "torch_missing"
        summary["error"] = str(exc)
        _write_summary(out_summary, summary)
        return 2

    try:
        pred_ts = _resolve_pred_ts()
    except Exception as exc:  # noqa: BLE001
        summary["status"] = "react_ot_import_failed"
        summary["error"] = repr(exc)
        summary["traceback"] = traceback.format_exc()
        _write_summary(out_summary, summary)
        return 3

    started = time.time()
    out_dir = Path(args.out_xyz).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    opt = SimpleNamespace(
        rxyz=args.reactant,
        pxyz=args.product,
        output_path=str(out_dir),
        batch_size=args.batch_size,
        nfe=args.nfe,
        solver=args.solver,
        checkpoint_path=args.checkpoint,
        order=args.order,
        diz=args.diz,
        method=args.method,
        atol=args.atol,
        rtol=args.rtol,
    )
    expected_ts_path = out_dir / (Path(args.reactant).stem + "_ts.xyz")
    try:
        pred_ts(args.reactant, args.product, opt, str(out_dir))
    except Exception as exc:  # noqa: BLE001
        summary["status"] = "inference_failed"
        summary["error"] = repr(exc)
        summary["traceback"] = traceback.format_exc()
        summary["wall_seconds"] = time.time() - started
        _write_summary(out_summary, summary)
        return 4

    summary["wall_seconds"] = time.time() - started
    if not expected_ts_path.exists():
        summary["status"] = "unexpected_result_type"
        summary["error"] = (
            f"expected React-OT to write {expected_ts_path}, but it did not exist"
        )
        _write_summary(out_summary, summary)
        return 5

    ts_xyz = expected_ts_path.read_text()
    Path(args.out_xyz).write_text(ts_xyz)
    summary["status"] = "ok"
    summary["out_xyz"] = str(Path(args.out_xyz).resolve())
    summary["react_ot_output_dir"] = str(out_dir)
    rxn_path = out_dir / (Path(args.reactant).stem + "_rxn.xyz")
    if rxn_path.exists():
        summary["out_rxn_xyz"] = str(rxn_path)
    _write_summary(out_summary, summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())

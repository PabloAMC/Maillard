"""src/react_ot_backend.py

Plugin module wiring React-OT (deepprinciple/react-ot) into the repo's
MLP backend adapter contract. React-OT lives in a *secondary* conda env
named ``react_ot`` (see ``scripts/setup_react_ot_env.sh``); this module
runs inside the main ``maillard`` env and reaches the secondary env via
``subprocess`` so the main env's torch / mace stack stays untouched.

Trust posture (also recorded in
``models/external/react_ot/provenance.json``):

* React-OT is treated strictly as a *geometric* TS-seed generator.
* React-OT energies are never propagated as runtime barrier authority.
* Every promising seed must clear downstream Sella DFT + imaginary-mode
  validation before any barrier change is considered.

The TS-seed-benchmark contract used by
``src.ts_seed_benchmark_validator`` feeds adapters a single challenged
seed XYZ; React-OT instead requires reactant + product XYZ. The two
contracts are not equivalent, so :func:`prepare_ts_seed` raises
``NotImplementedError`` and the adapter reports
``available=False`` to the validator. The pilot wrapper
``scripts/recover_ts_react_ot_seed.py`` calls
:func:`predict_ts_from_reactant_product` directly.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


REACT_OT_ENV_NAME = os.environ.get("MAILLARD_REACT_OT_ENV", "react_ot")
CONDA_SH = os.environ.get("MAILLARD_CONDA_SH", "/opt/conda/etc/profile.d/conda.sh")
DEFAULT_INFERENCE_SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "_react_ot_inference.py"
)
DEFAULT_CHECKPOINT_ENV = "MAILLARD_REACT_OT_CHECKPOINT"


def _bash_in_secondary_env(command: str) -> subprocess.CompletedProcess:
    full = (
        f"set -eo pipefail; "
        f"source '{CONDA_SH}'; "
        f"conda activate '{REACT_OT_ENV_NAME}'; "
        f"{command}"
    )
    return subprocess.run(
        ["bash", "-lc", full],
        capture_output=True,
        text=True,
        check=False,
    )


def probe_backend(model_name: str, locator: Optional[str]) -> Tuple[bool, str]:
    """Check that the secondary env exists and React-OT imports.

    Returns ``(True, reason)`` only when both ``conda activate
    <secondary>`` succeeds and ``import reactot`` succeeds inside it.
    Otherwise returns ``(False, reason)`` with a diagnostic suitable for
    :class:`BackendAvailability`.
    """
    del model_name, locator  # adapter contract symmetry, not used here

    if not Path(CONDA_SH).exists() and not shutil.which("conda"):
        return False, (
            f"conda not available; expected activation script at {CONDA_SH}. "
            "Run scripts/setup_react_ot_env.sh inside the maillard_validation "
            "container."
        )

    proc = _bash_in_secondary_env(
        "python -c 'import reactot' >/dev/null 2>&1 && echo OK || echo MISSING"
    )
    if proc.returncode != 0:
        diag = (proc.stderr or proc.stdout).strip()
        return False, (
            f"could not activate secondary env '{REACT_OT_ENV_NAME}': {diag}"
        )
    if "OK" not in proc.stdout:
        return False, (
            f"react-ot not importable in secondary env '{REACT_OT_ENV_NAME}'; "
            "run scripts/setup_react_ot_env.sh"
        )
    return True, (
        f"react-ot importable in secondary env '{REACT_OT_ENV_NAME}'"
    )


def prepare_ts_seed(
    xyz_string: str,
    model_name: str,
    locator: Optional[str],
    fmax: float,
    max_steps: int,
) -> str:
    """Not supported: React-OT requires reactant and product geometries.

    The TS-seed-benchmark contract feeds a challenged single-structure
    seed; React-OT's contract is fundamentally different (it generates a
    TS from a reactant + product pair, not from a single seed). Use
    :func:`predict_ts_from_reactant_product` from the pilot script
    ``scripts/recover_ts_react_ot_seed.py`` instead.
    """
    del xyz_string, model_name, locator, fmax, max_steps
    raise NotImplementedError(
        "React-OT requires reactant and product XYZ inputs and does not "
        "implement the single-seed prepare_ts_seed contract. Use "
        "predict_ts_from_reactant_product via "
        "scripts/recover_ts_react_ot_seed.py."
    )


def predict_ts_from_reactant_product(
    reactant_xyz: str,
    product_xyz: str,
    *,
    checkpoint: Optional[str] = None,
    solver: str = "ode",
    nfe: int = 10,
    n_candidates: int = 1,
    device: str = "cpu",
    inference_script: Path = DEFAULT_INFERENCE_SCRIPT,
) -> Dict[str, Any]:
    """Run React-OT inference in the secondary env via subprocess.

    Returns a dict with at least ``status`` and, on success, ``xyz``,
    ``out_xyz`` and ``wall_seconds``. The TS-guess is treated as a
    geometric seed only; React-OT energies are never propagated as
    runtime barrier authority.
    """
    resolved_checkpoint = checkpoint or os.environ.get(DEFAULT_CHECKPOINT_ENV)
    if not resolved_checkpoint:
        return {
            "status": "missing_checkpoint",
            "error": (
                f"set --checkpoint or environment variable "
                f"{DEFAULT_CHECKPOINT_ENV}"
            ),
        }
    if not Path(resolved_checkpoint).exists():
        return {
            "status": "checkpoint_not_found",
            "checkpoint": str(resolved_checkpoint),
        }

    with tempfile.TemporaryDirectory(prefix="react_ot_") as tmp:
        tmp_path = Path(tmp)
        reactant_path = tmp_path / "reactant.xyz"
        product_path = tmp_path / "product.xyz"
        out_xyz = tmp_path / "ts_guess.xyz"
        out_summary = tmp_path / "ts_summary.json"

        reactant_path.write_text(reactant_xyz)
        product_path.write_text(product_xyz)

        command = (
            f"python {str(inference_script)!r} "
            f"--reactant {str(reactant_path)!r} "
            f"--product {str(product_path)!r} "
            f"--checkpoint {str(Path(resolved_checkpoint).resolve())!r} "
            f"--out-xyz {str(out_xyz)!r} "
            f"--out-summary {str(out_summary)!r} "
            f"--solver {solver} --nfe {int(nfe)} "
            f"--n-candidates {int(n_candidates)} "
            f"--device {device}"
        )
        proc = _bash_in_secondary_env(command)
        if out_summary.exists():
            summary: Dict[str, Any] = json.loads(out_summary.read_text())
        else:
            summary = {
                "status": "subprocess_failed",
                "stdout": proc.stdout,
                "stderr": proc.stderr,
                "returncode": proc.returncode,
            }
        if out_xyz.exists():
            summary["xyz"] = out_xyz.read_text()
        return summary


__all__ = [
    "DEFAULT_CHECKPOINT_ENV",
    "DEFAULT_INFERENCE_SCRIPT",
    "REACT_OT_ENV_NAME",
    "predict_ts_from_reactant_product",
    "prepare_ts_seed",
    "probe_backend",
]

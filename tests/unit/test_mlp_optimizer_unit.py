from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mlp_optimizer import (  # noqa: E402
    _resolve_backend_locator,
    _should_reject_bounded_preopt,
    compute_geometry_drift_metrics,
)


CYSTEINE_XYZ = """14
Cysteine
C -0.531 -0.103 0.093
C 0.916 -0.345 -0.341
O 1.309 -1.464 -0.627
O 1.731 0.722 -0.380
N -1.282 -1.328 -0.218
C -0.643 0.211 1.587
S -2.392 0.536 2.112
H -1.002 0.724 -0.457
H -2.227 -1.168 0.124
H -1.341 -1.472 -1.229
H 2.637 0.490 -0.655
H -0.088 -0.627 2.036
H -0.223 1.149 1.956
H -2.325 1.761 1.528
"""


def test_mace_omol_requires_explicit_backend_locator():
    try:
        _resolve_backend_locator("mace_omol", "large", None)
    except ImportError as exc:
        assert "explicit backend locator" in str(exc)
    else:
        raise AssertionError("MACE-OMol should require an explicit backend locator")


def test_compute_geometry_drift_metrics_reports_sulfur_local_delta():
    perturbed_xyz = CYSTEINE_XYZ.replace("S -2.392 0.536 2.112", "S -2.050 0.536 2.112")

    metrics = compute_geometry_drift_metrics(CYSTEINE_XYZ, perturbed_xyz)

    assert metrics["max_atom_displacement_angstrom"] is not None
    assert float(metrics["max_atom_displacement_angstrom"] or 0.0) > 0.1
    assert metrics["sulfur_local_rmsd_angstrom"] is not None
    assert metrics["sulfur_neighbor_max_delta_angstrom"] is not None
    assert metrics["sulfur_bond_max_delta_angstrom"] is not None
    assert metrics["sulfur_angle_max_delta_degrees"] is not None
    assert float(metrics["sulfur_neighbor_max_delta_angstrom"] or 0.0) > 0.05


def test_bounded_preopt_guard_rejects_large_local_sulfur_drift():
    metrics = {
        "max_atom_displacement_angstrom": 0.35,
        "sulfur_neighbor_max_delta_angstrom": 0.31,
    }

    assert _should_reject_bounded_preopt(
        metrics,
        max_atom_drift_threshold=0.8,
        sulfur_neighbor_delta_threshold=0.25,
    ) is True
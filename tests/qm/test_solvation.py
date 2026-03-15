import os
import pytest
import numpy as np
from src.solvation import SolvationEngine  # noqa: E402

# Simple water dimer as a test case
WATER_XYZ = """3
water
O 0.0 0.0 0.0
H 0.0 0.0 0.96
H 0.0 0.92 -0.25
"""

def test_solvation_engine_init():
    """Verify engine initializes and finds crest binary."""
    engine = SolvationEngine()
    assert engine.crest_bin is not None
    assert os.path.isabs(engine.crest_bin) or not ("/" in engine.crest_bin or "\\" in engine.crest_bin) or os.path.exists(engine.crest_bin)

def test_n_atoms_parsing():
    """Verify the internal atom count parser."""
    engine = SolvationEngine()
    assert engine._parse_n_atoms(WATER_XYZ) == 3

_CREST_AVAILABLE, _CREST_SKIP_REASON = SolvationEngine.probe_qcg_capability()

@pytest.mark.skipif(not _CREST_AVAILABLE, reason=_CREST_SKIP_REASON)
def test_explicit_solvation_run():
    """Run a real CREST/QCG call if the binary is present."""
    engine = SolvationEngine()
    
    # Add 1 water to a water molecule (3 atoms -> 6 atoms)
    res_xyz = engine.generate_solvated_cluster(WATER_XYZ, n_water=1, freeze_core=True)

    lines = res_xyz.strip().split('\n')
    n_atoms = int(lines[0])
    assert n_atoms == 6

    # Parse the actual coordinates for the first 3 atoms (the solute)
    solute_coords = np.array([
        [float(p) for p in line.split()[1:4]] 
        for line in lines[2:5]
    ])

    # The original expected coordinates
    expected_coords = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.96],
        [0.0, 0.92, -0.25]
    ])

    # Helper function to compute the pairwise distance matrix
    def get_distance_matrix(coords):
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        return np.linalg.norm(diff, axis=-1)

    # Calculate internal distances
    actual_dist = get_distance_matrix(solute_coords)
    expected_dist = get_distance_matrix(expected_coords)

    # Verify the internal structure hasn't changed (rigid body)
    np.testing.assert_allclose(actual_dist, expected_dist, atol=5e-2, 
                               err_msg="Solute internal geometry was altered despite freeze_core=True")

def test_error_handling_on_missing_bin():
    """Verify that a missing or invalid binary raises a RuntimeError."""
    engine = SolvationEngine(crest_bin="/path/to/nonexistent/crest")
    with pytest.raises(RuntimeError) as excinfo:
        engine.generate_solvated_cluster(WATER_XYZ, n_water=1, freeze_core=True)
    
    assert "Solvation error" in str(excinfo.value)
    assert "No such file or directory" in str(excinfo.value)
import pytest
import os
import numpy as np
from src.mlp_optimizer import MLPOptimizer

# Cysteine SMILES for drift test (sulfur species)
CYSTEINE_XYZ = """14
Cysteine
C   -0.531  -0.103   0.093
C    0.916  -0.345  -0.341
O    1.309  -1.464  -0.627
O    1.731   0.722  -0.380
N   -1.282  -1.328  -0.218
C   -0.643   0.211   1.587
S   -2.392   0.536   2.112
H   -1.002   0.724  -0.457
H   -2.227  -1.168   0.124
H   -1.341  -1.472  -1.229
H    2.637   0.490  -0.655
H   -0.088  -0.627   2.036
H   -0.223   1.149   1.956
H   -2.325   1.761   1.528
"""

def mace_not_found():
    try:
        import mace
        return False
    except ImportError:
        return True

@pytest.mark.skipif(mace_not_found(), reason="MACE not installed in current environment")
def test_mace_drift_cysteine():
    """
    Verify MACE optimization on a sulfur species.
    Checks if the drift is within acceptable bounds (per todo.md, it might drift).
    """
    optimizer = MLPOptimizer(device="cpu")
    opt_xyz = optimizer.optimize_geometry(CYSTEINE_XYZ, fmax=0.01)
    
    # Parse coordinates
    lines = opt_xyz.strip().split('\n')
    opt_coords = []
    for line in lines[2:]:
        parts = line.split()
        opt_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    # Original coords
    orig_lines = CYSTEINE_XYZ.strip().split('\n')
    orig_coords = []
    for line in orig_lines[2:]:
        parts = line.split()
        orig_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        
    # Calculate max displacement
    displacements = np.linalg.norm(np.array(opt_coords) - np.array(orig_coords), axis=1)
    max_drift = np.max(displacements)
    
    print(f"Max MACE drift for Cysteine: {max_drift:.3f} Å")
    # According to todo.md, sulfur species drift can be ~1.11 Å
    assert max_drift < 2.0 # Sanity check for "non-exploding" optimization

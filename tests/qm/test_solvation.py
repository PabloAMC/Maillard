import os
import shutil
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
    # If it's a relative path ("crest"), it's fine as a fallback to PATH lookup
    assert os.path.isabs(engine.crest_bin) or not ("/" in engine.crest_bin or "\\" in engine.crest_bin) or os.path.exists(engine.crest_bin)

def test_n_atoms_parsing():
    """Verify the internal atom count parser."""
    engine = SolvationEngine()
    assert engine._parse_n_atoms(WATER_XYZ) == 3

def _has_crest() -> bool:
    engine = SolvationEngine()
    # Check if the detected path is valid and executable
    return bool(engine.crest_bin and (os.path.exists(engine.crest_bin) or shutil.which(engine.crest_bin)))

def _has_xtb_iff() -> bool:
    """CREST's QCG mode requires xtb-IFF; without it CREST segfaults (exit -11)."""
    return shutil.which("xtb-IFF") is not None

_CREST_QCG_AVAILABLE = _has_crest() and _has_xtb_iff()
_CREST_QCG_SKIP_REASON = (
    "CREST binary not found" if not _has_crest()
    else "xtb-IFF not found (required by CREST -qcg mode)"
)

@pytest.mark.skipif(not _CREST_QCG_AVAILABLE, reason=_CREST_QCG_SKIP_REASON)
def test_explicit_solvation_run():
    """Run a real CREST/QCG call if the binary is present."""
    engine = SolvationEngine()
    # Add 1 water to a water molecule (3 atoms -> 6 atoms)
    res_xyz = engine.generate_solvated_cluster(WATER_XYZ, n_water=1, freeze_core=True)
    
    lines = res_xyz.strip().split('\n')
    n_atoms = int(lines[0])
    assert n_atoms == 6
    
    # Verify the first 3 atoms (solute) haven't moved significantly (freeze_core works)
    # The header line is line 0, comment is line 1
    solute_coords = []
    for line in lines[2:5]:
        parts = line.split()
        solute_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    expected_coords = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.96],
        [0.0, 0.92, -0.25]
    ]
    
    for i in range(3):
        np.testing.assert_allclose(solute_coords[i], expected_coords[i], atol=1e-3)

def test_error_handling_on_missing_bin():
    """Verify that a missing or invalid binary raises a RuntimeError."""
    engine = SolvationEngine(crest_bin="/path/to/nonexistent/crest")
    with pytest.raises(RuntimeError) as excinfo:
        engine.generate_solvated_cluster(WATER_XYZ, n_water=1, freeze_core=True)
    
    assert "Solvation error" in str(excinfo.value)
    assert "No such file or directory" in str(excinfo.value)

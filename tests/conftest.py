"""
Pytest configuration and shared fixtures for Phase 3 and Phase 12 tests.
"""

import pytest
import numpy as np


# ============================================================================
# PHASE 3 FIXTURES: DFT Refiner, Quasi-Harmonic, Barriers, IRC, Verification
# ============================================================================

@pytest.fixture
def dft_geometry_water():
    """Pre-optimized water molecule (XYZ format)."""
    # H2O: O at origin, H atoms ~0.96 Å away
    atoms = ['O', 'H', 'H']
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [-0.24, 0.93, 0.0]
    ])
    return atoms, coords


@pytest.fixture
def dft_geometry_formaldehyde():
    """Pre-optimized formaldehyde (CH2O)."""
    atoms = ['C', 'O', 'H', 'H']
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.20, 0.0, 0.0],
        [-0.50, 0.89, 0.0],
        [-0.50, -0.89, 0.0]
    ])
    return atoms, coords


@pytest.fixture
def dft_geometry_ethane():
    """Pre-optimized ethane (C2H6)."""
    atoms = ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.54, 0.0, 0.0],
        [-0.51, 0.89, 0.0],
        [-0.51, -0.44, 0.77],
        [-0.51, -0.44, -0.77],
        [2.05, 0.89, 0.0],
        [2.05, -0.44, 0.77],
        [2.05, -0.44, -0.77]
    ])
    return atoms, coords


@pytest.fixture
def literature_barriers_dict():
    """Dictionary of literature DFT barrier values (kcal/mol)."""
    return {
        "amadori_rearrangement": {"low": 18, "high": 24, "best": 20},
        "enolisation_1_2": {"low": 18, "high": 22, "best": 20},
        "enolisation_2_3": {"low": 22, "high": 28, "best": 25},
        "strecker_decarboxylation": {"low": 22, "high": 28, "best": 24},
        "fft_formation": {"low": 20, "high": 30, "best": 25},
        "retro_aldol": {"low": 20, "high": 26, "best": 23},
        "dha_elimination": {"low": 24, "high": 30, "best": 27},
        "trapping_hexanal": {"low": 15, "high": 22, "best": 18},
        "pyrazine_aminoketone": {"low": 25, "high": 35, "best": 28},
    }


@pytest.fixture
def xtb_barriers_dict():
    """Dictionary of xTB reference barrier values (kcal/mol)."""
    return {
        "amadori_rearrangement": 19,
        "enolisation_1_2": 21,
        "enolisation_2_3": 26,
        "strecker_decarboxylation": 25,
        "fft_formation": 27,
        "retro_aldol": 22,
        "dha_elimination": 28,
        "trapping_hexanal": 17,
        "pyrazine_aminoketone": 30,
    }


@pytest.fixture
def dft_ts_geometries():
    """Dictionary of pre-computed TS geometries (minimal representations)."""
    return {
        "amadori": {"atoms": ['C', 'N', 'O', 'H', 'H'], "energy": -100.5},
        "strecker": {"atoms": ['C', 'C', 'N', 'O', 'O', 'H'], "energy": -99.2},
        "pyrazine": {"atoms": ['C', 'C', 'C', 'N', 'N', 'H', 'H'], "energy": -98.8},
    }


@pytest.fixture
def temperature_list():
    """Common test temperatures (K)."""
    return {
        "room_temperature": 298,
        "maillard_low": 373,      # 100°C
        "maillard_mid": 423,      # 150°C (commonly used)
        "maillard_high": 473,     # 200°C
    }


# ============================================================================
# PHASE 12 FIXTURES: Cantera Integration and CLI
# ============================================================================

@pytest.fixture
def initial_state_simple():
    """Simple initial state for isothermal simulations."""
    return {
        "ribose": 0.1,
        "cysteine": 0.1,
    }


@pytest.fixture
def initial_state_complex():
    """Complex initial state with multiple precursors."""
    return {
        "ribose": 0.05,
        "cysteine": 0.05,
        "glucose": 0.05,
        "glycine": 0.05,
    }


@pytest.fixture
def formulations_list():
    """List of test formulations for sensory prediction."""
    return [
        {"name": "S-Maillard", "precursors": {"ribose": 0.1, "cysteine": 0.1}},
        {"name": "N-Maillard", "precursors": {"glucose": 0.1, "glycine": 0.1}},
        {"name": "Mixed", "precursors": {"ribose": 0.1, "cysteine": 0.05}},
    ]


@pytest.fixture
def temperature_profile_linear():
    """Linear temperature ramp function: RT → 150°C over 30 min."""
    def temp_func(t):
        """Return temperature (K) as function of time (s)."""
        t_ramp = 1800  # 30 minutes
        T_initial = 298  # RT (K)
        T_final = 423  # 150°C
        if t < t_ramp:
            return T_initial + (T_final - T_initial) * t / t_ramp
        else:
            return T_final
    return temp_func


@pytest.fixture
def cli_args_basic():
    """Basic CLI arguments for kinetic simulation."""
    return [
        "python", "scripts/run_cantera_kinetics.py",
        "--precursors", "ribose:0.1,cysteine:0.1",
        "--temperature", "150",
        "--output", "test_output.csv"
    ]


# ============================================================================
# UTILITY FIXTURES
# ============================================================================

@pytest.fixture
def tmp_geometry_file(tmp_path):
    """Temporary XYZ geometry file for testing."""
    xyz_path = tmp_path / "test.xyz"
    xyz_content = """3
Water molecule
O    0.0    0.0    0.0
H    0.96   0.0    0.0
H   -0.24   0.93   0.0
"""
    xyz_path.write_text(xyz_content)
    return xyz_path


@pytest.fixture
def mock_mechanism_dict():
    """Mock RMG mechanism as dictionary."""
    return {
        "species": [
            {"name": "ribose", "smiles": "OC[C@@H](O)[C@H](O)C=O"},
            {"name": "cysteine", "smiles": "N[C@H](CS)C(=O)O"},
            {"name": "furfural", "smiles": "O=Cc1ccco1"},
            {"name": "FFT", "smiles": "Cc1csc[nH]1"},
        ],
        "reactions": [
            {
                "reactants": ["ribose", "cysteine"],
                "products": ["furfural", "FFT"],
                "barrier_kcal_mol": 25.0
            }
        ]
    }


@pytest.fixture
def expected_sensory_profile():
    """Expected sensory profile dictionary."""
    return {
        "meaty": 0.0,
        "roasted": 0.0,
        "fruity": 0.0,
        "floral": 0.0,
        "beany": 0.0,
    }


# ============================================================================
# PYTEST MARKERS
# ============================================================================

def pytest_configure(config):
    """Register custom pytest markers."""
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "dft: marks Phase 3 DFT tests")
    config.addinivalue_line("markers", "cantera: marks Phase 12 Cantera tests")
    config.addinivalue_line("markers", "integration: marks end-to-end integration tests")

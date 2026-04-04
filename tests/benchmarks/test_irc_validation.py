"""Test suite for mounted Phase 3.4 IRC authority-lane fixtures."""

import numpy as np
import pytest

from src.authority_benchmark_data import load_irc_validation_cases
from src.dft_refiner import DFTRefiner, compute_irc

from tests.benchmarks._lane_policy import (
    HAS_IRC_FIXTURES,
    HAS_IRC_IMPLEMENTATION,
    IRC_SKIP_REASON,
)


pytestmark = [
    pytest.mark.slow,
    pytest.mark.optional_dft_authority_lane,
    pytest.mark.skipif(not (HAS_IRC_IMPLEMENTATION and HAS_IRC_FIXTURES), reason=IRC_SKIP_REASON),
]


IRC_CASES = load_irc_validation_cases() if (HAS_IRC_IMPLEMENTATION and HAS_IRC_FIXTURES) else []
IRC_CASES_BY_ID = {case["case_id"]: case for case in IRC_CASES}


def _coords_from_xyz(xyz_text: str) -> np.ndarray:
    lines = xyz_text.strip().splitlines()[2:]
    coords = []
    for line in lines:
        _, x_coord, y_coord, z_coord = line.split()
        coords.append([float(x_coord), float(y_coord), float(z_coord)])
    return np.array(coords)


def _install_fixture_backend(monkeypatch: pytest.MonkeyPatch, case_id: str) -> dict[str, object]:
    case = IRC_CASES_BY_ID[case_id]
    energy_by_xyz = {
        case["backward_endpoint_xyz"]: case["energies"][0],
        case["ts_xyz"]: case["energies"][1],
        case["forward_endpoint_xyz"]: case["energies"][2],
    }

    def fake_init(self, solvent_name="water", temp_k=423.15, geometry_backend="pyscf", **_kwargs):
        self.solvent_name = solvent_name
        self.temp_k = temp_k
        self.geometry_backend = geometry_backend

    def fake_generate_irc(self, ts_xyz, charge=0, spin=0, step_size=0.05):
        assert charge == 0
        assert spin == 0
        assert step_size > 0.0
        assert ts_xyz == case["ts_xyz"]
        return case["backward_endpoint_xyz"], case["forward_endpoint_xyz"]

    def fake_single_point(self, xyz_content, xc_method="wB97M-V", basis="def2-tzvp", charge=0, spin=0):
        assert charge == 0
        assert spin == 0
        assert xc_method
        assert basis
        return float(energy_by_xyz[xyz_content])

    monkeypatch.setattr(DFTRefiner, "__init__", fake_init)
    monkeypatch.setattr(DFTRefiner, "generate_irc", fake_generate_irc)
    monkeypatch.setattr(DFTRefiner, "single_point", fake_single_point)
    return case


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_path_continuity(case_id, monkeypatch):
    """Mounted IRC fixtures are wired through the stable compute_irc API."""
    case = _install_fixture_backend(monkeypatch, case_id)
    irc_path = compute_irc(case["ts_xyz"])

    assert irc_path["backward_endpoint"] == case["backward_endpoint_xyz"]
    assert irc_path["forward_endpoint"] == case["forward_endpoint_xyz"]
    for index in range(1, len(irc_path["energies"])):
        energy_jump = abs(irc_path["energies"][index] - irc_path["energies"][index - 1])
        assert energy_jump < 0.05


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_energy_profile_smooth(case_id, monkeypatch):
    """Mounted IRC energy traces remain smooth enough for a 3-point proxy path."""
    case = _install_fixture_backend(monkeypatch, case_id)
    irc_path = compute_irc(case["ts_xyz"])
    second_derivative = np.diff(np.array(irc_path["energies"]), n=2)

    assert np.all(np.abs(second_derivative) < 0.05)


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_energy_is_minimum_at_endpoints(case_id, monkeypatch):
    """Fixture endpoints remain lower in energy than the transition state."""
    case = _install_fixture_backend(monkeypatch, case_id)
    irc_path = compute_irc(case["ts_xyz"])
    energies = irc_path["energies"]

    assert energies[0] < energies[1]
    assert energies[2] < energies[1]


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_ts_at_maximum(case_id, monkeypatch):
    """The transition state remains at the energy maximum of the proxy path."""
    case = _install_fixture_backend(monkeypatch, case_id)
    irc_path = compute_irc(case["ts_xyz"])
    assert int(np.argmax(irc_path["energies"])) == len(irc_path["energies"]) // 2


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_endpoint_geometries_match_mounted_fixtures(case_id, monkeypatch):
    """compute_irc returns the mounted endpoint geometries without drift."""
    case = _install_fixture_backend(monkeypatch, case_id)
    irc_path = compute_irc(case["ts_xyz"])

    assert irc_path["backward_endpoint"] == case["backward_endpoint_xyz"]
    assert irc_path["forward_endpoint"] == case["forward_endpoint_xyz"]


@pytest.mark.parametrize("case_id", sorted(IRC_CASES_BY_ID))
def test_irc_barrier_consistency(case_id):
    """Mounted direct barriers remain consistent with the IRC energy gap."""
    case = IRC_CASES_BY_ID[case_id]
    barrier_from_irc = (case["energies"][1] - case["energies"][0]) * 627.509
    assert abs(barrier_from_irc - case["direct_barrier_kcal_mol"]) < 0.1


def test_irc_reaction_3_3c_strecker():
    """The mounted Strecker case remains present and annotated as CO2 elimination."""
    strecker_case = next(case for case in IRC_CASES if case["family"] == "strecker_decarboxylation")
    assert strecker_case["expected_bond_change"] == "co2_elimination"
    assert strecker_case["product_smiles"] == "CC(=O)CN.C(=O)=O"


def test_irc_reaction_geometry_change():
    """The mounted Amadori proxy shows a meaningful endpoint geometry displacement."""
    amadori_case = next(case for case in IRC_CASES if case["family"] == "amadori_rearrangement")
    reactant_coords = _coords_from_xyz(amadori_case["backward_endpoint_xyz"])
    product_coords = _coords_from_xyz(amadori_case["forward_endpoint_xyz"])
    displacement = np.linalg.norm(product_coords - reactant_coords, axis=1)

    assert float(np.mean(displacement)) > 0.5

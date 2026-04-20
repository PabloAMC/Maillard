"""Unit tests for src.microsolvation: pure-geometry primitives only.

Tests requiring xTB (relax_solvent, build_hydrated_endpoints) are exercised
in scientific tests where the docker container provides xtb.
"""
from __future__ import annotations

import numpy as np
import pytest

from src.microsolvation import (
    HOH_ANGLE_DEG,
    OH_BOND_A,
    check_no_clash,
    kabsch_align,
    make_water,
    molecular_normal,
    place_proton_shuttle_water,
    read_xyz,
    write_xyz,
)


def test_kabsch_align_identity():
    """Aligning a geometry to itself yields identity rotation and zero RMSD."""
    rng = np.random.default_rng(42)
    coords = rng.standard_normal((10, 3))
    rot, mc, tc = kabsch_align(coords, coords)
    aligned = (coords - mc) @ rot + tc
    assert np.allclose(rot, np.eye(3), atol=1e-6)
    assert np.allclose(aligned, coords, atol=1e-6)


def test_kabsch_align_recovers_known_rotation():
    rng = np.random.default_rng(1)
    coords = rng.standard_normal((20, 3))
    # Rotate by a known matrix and translate
    theta = np.pi / 6
    R_true = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                       [np.sin(theta), np.cos(theta), 0.0],
                       [0.0, 0.0, 1.0]])
    target = coords @ R_true.T + np.array([1.0, -2.0, 3.0])
    rot, mc, tc = kabsch_align(coords, target)
    aligned = (coords - mc) @ rot + tc
    rmsd = float(np.sqrt(((aligned - target) ** 2).sum() / len(coords)))
    assert rmsd < 1e-6


def test_kabsch_align_shape_mismatch_raises():
    with pytest.raises(ValueError):
        kabsch_align(np.zeros((3, 3)), np.zeros((4, 3)))


def test_make_water_oh_bond_and_angle():
    """Generated water has OH=0.97 Å and HOH≈104.5°."""
    o = np.array([0.0, 0.0, 0.0])
    donor = np.array([2.0, 0.0, 0.0])
    acceptor = np.array([-2.0, 0.0, 0.0])
    w = make_water(o, donor, acceptor)
    assert w.shape == (3, 3)
    o_pos, h1, h2 = w
    oh1 = float(np.linalg.norm(h1 - o_pos))
    oh2 = float(np.linalg.norm(h2 - o_pos))
    assert abs(oh1 - OH_BOND_A) < 1e-6
    assert abs(oh2 - OH_BOND_A) < 1e-6
    v1 = (h1 - o_pos) / oh1
    v2 = (h2 - o_pos) / oh2
    angle_deg = float(np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))))
    assert abs(angle_deg - HOH_ANGLE_DEG) < 1e-3


def test_molecular_normal_unit_vector():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                       [1.0, 1.0, 0.0]])
    n = molecular_normal(coords)
    assert abs(np.linalg.norm(n) - 1.0) < 1e-6
    # All points lie in the z=0 plane → normal is along ±z
    assert abs(abs(n[2]) - 1.0) < 1e-6


def test_place_proton_shuttle_water_no_overlap_with_endpoints():
    """The placed water O is offset from the donor-acceptor midpoint."""
    coords = np.array([
        [0.0, 0.0, 0.0],   # donor (idx 0)
        [2.5, 0.0, 0.0],   # acceptor (idx 1)
        [1.0, 0.0, -1.0],  # padding atom
        [1.0, 0.0, 1.0],
    ])
    w = place_proton_shuttle_water(coords, 0, 1, offset=1.6)
    assert w.shape == (3, 3)
    o_pos = w[0]
    d_to_donor = float(np.linalg.norm(o_pos - coords[0]))
    d_to_acceptor = float(np.linalg.norm(o_pos - coords[1]))
    # Water O should be roughly equidistant from donor and acceptor
    assert abs(d_to_donor - d_to_acceptor) < 1e-3
    # And further than 1.0 Å from each (no overlap)
    assert d_to_donor > 1.0


def test_check_no_clash_detects_overlap():
    syms = ["O", "H", "H", "O", "H", "H"]
    # Two waters with O atoms 0.5 Å apart → clash
    coords = np.array([
        [0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [-0.24, 0.94, 0.0],
        [0.5, 0.0, 0.0], [1.47, 0.0, 0.0], [0.26, 0.94, 0.0],
    ])
    rep = check_no_clash(syms, coords, n_solute=0)
    assert rep.has_clash
    assert rep.min_distance_a < 0.85


def test_check_no_clash_skips_intra_water_oh():
    """Intra-water O–H bonds (~0.97 Å) must not trigger a clash flag."""
    syms = ["C", "H", "H", "H", "H", "O", "H", "H"]
    # Solute (5 atoms) far from a water (3 atoms)
    coords = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0],
        [10.0, 0.0, 0.0], [10.97, 0.0, 0.0], [9.76, 0.94, 0.0],
    ])
    rep = check_no_clash(syms, coords, n_solute=5)
    # Min nonbonded distance is between solute atoms (1.0 Å) — not a clash
    assert not rep.has_clash
    assert rep.min_distance_a >= 1.0 - 1e-6


def test_xyz_round_trip(tmp_path):
    syms = ["C", "H", "O"]
    coords = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [-1.2, 0.0, 0.0]])
    p = tmp_path / "test.xyz"
    write_xyz(p, syms, coords, comment="test")
    syms2, coords2 = read_xyz(p)
    assert syms2 == syms
    assert np.allclose(coords2, coords, atol=1e-6)

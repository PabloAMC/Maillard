"""Tests for src.xyz_fragment_utils — fragment splitting and charge/spin assignment."""

import pytest
from src.xyz_fragment_utils import (
    split_xyz_into_fragments,
    assign_fragment_charge_spin,
    count_unpaired_electrons,
    _parse_xyz,
    _build_adjacency,
    _connected_components,
)


# ── Fixtures ───────────────────────────────────────────────────────────────

HEXANAL_SH_XYZ = """\
21
hexanal + SH radical
C      1.3551    0.0958   -0.1293
C      2.7921   -0.2862    0.2060
C      3.7519    0.8770   -0.0035
C      5.1962    0.5166    0.3553
C      6.1520    1.6688    0.1273
C      7.5563    1.2898    0.4729
O      7.9137    0.1672    0.7795
H      0.6811   -0.6897    0.2093
H      1.1168    0.1447   -1.1964
H      1.1855    1.0767    0.3267
H      3.1137   -1.1791   -0.3446
H      2.8410   -0.5878    1.2596
H      3.4387    1.7780    0.5362
H      3.6914    1.1572   -1.0626
H      5.5113   -0.3798   -0.1888
H      5.2529    0.2321    1.4139
H      5.8346    2.5614    0.6767
H      6.0906    1.9526   -0.9307
H      8.2798    2.1081    0.3592
H     10.6551    0.1528    0.6795
S     11.8677   -0.2072    0.1413
"""

WATER_XYZ = """\
3
water molecule
O      0.0000    0.0000    0.0000
H      0.7572    0.5860    0.0000
H     -0.7572    0.5860    0.0000
"""

# Two water molecules, far apart (no bond between them)
TWO_WATERS_XYZ = """\
6
two separate waters
O      0.0000    0.0000    0.0000
H      0.7572    0.5860    0.0000
H     -0.7572    0.5860    0.0000
O     10.0000    0.0000    0.0000
H     10.7572    0.5860    0.0000
H      9.2428    0.5860    0.0000
"""


# ── split_xyz_into_fragments ──────────────────────────────────────────────

def test_single_molecule_returns_one_fragment():
    frags = split_xyz_into_fragments(WATER_XYZ)
    assert len(frags) == 1


def test_two_waters_detected_as_two_fragments():
    frags = split_xyz_into_fragments(TWO_WATERS_XYZ)
    assert len(frags) == 2
    for frag in frags:
        syms, _ = _parse_xyz(frag)
        assert len(syms) == 3
        assert syms.count("O") == 1
        assert syms.count("H") == 2


def test_hexanal_sh_detected_as_two_fragments():
    frags = split_xyz_into_fragments(HEXANAL_SH_XYZ)
    assert len(frags) == 2
    sizes = sorted(len(_parse_xyz(f)[0]) for f in frags)
    assert sizes == [2, 19]  # SH (2 atoms) + hexanal (19 atoms)


def test_fragment_xyz_is_valid():
    """Each fragment XYZ should be parseable and have correct atom count header."""
    frags = split_xyz_into_fragments(HEXANAL_SH_XYZ)
    for frag in frags:
        lines = frag.strip().splitlines()
        n_claimed = int(lines[0].strip())
        syms, coords = _parse_xyz(frag)
        assert len(syms) == n_claimed
        assert len(coords) == n_claimed


# ── count_unpaired_electrons ──────────────────────────────────────────────

def test_sh_radical_is_doublet():
    assert count_unpaired_electrons(["S", "H"]) == 1  # 17 electrons → odd


def test_water_is_singlet():
    assert count_unpaired_electrons(["O", "H", "H"]) == 0  # 10 electrons → even


def test_hexanal_is_singlet():
    # C6H12O = 6*6 + 12*1 + 8 = 56 electrons → even
    syms = ["C"] * 6 + ["H"] * 12 + ["O"]
    assert count_unpaired_electrons(syms) == 0


# ── assign_fragment_charge_spin ───────────────────────────────────────────

def test_explicit_fragment_spins():
    frags = split_xyz_into_fragments(HEXANAL_SH_XYZ)
    result = assign_fragment_charge_spin(frags, total_charge=0, total_spin=1,
                                         manifest_fragment_spins=[0, 1])
    assert result == [(0, 0), (0, 1)]


def test_auto_detect_fragment_spins():
    frags = split_xyz_into_fragments(HEXANAL_SH_XYZ)
    result = assign_fragment_charge_spin(frags, total_charge=0, total_spin=1)
    spins = [s for _, s in result]
    assert sum(spins) == 1  # total spin = 1 (doublet)
    # The SH fragment (2 atoms, odd electrons) should get spin=1
    sizes = [len(_parse_xyz(f)[0]) for f in frags]
    sh_idx = sizes.index(2)
    assert result[sh_idx] == (0, 1)


def test_single_fragment_charge_spin():
    frags = split_xyz_into_fragments(WATER_XYZ)
    result = assign_fragment_charge_spin(frags, total_charge=0, total_spin=0)
    assert result == [(0, 0)]


# ── _build_adjacency / _connected_components ──────────────────────────────

def test_adjacency_water():
    syms, coords = _parse_xyz(WATER_XYZ)
    adj = _build_adjacency(syms, coords)
    assert len(adj) == 3
    assert 0 in adj[1] and 0 in adj[2]  # H's bonded to O


def test_components_two_waters():
    syms, coords = _parse_xyz(TWO_WATERS_XYZ)
    adj = _build_adjacency(syms, coords)
    comps = _connected_components(adj)
    assert len(comps) == 2
    assert sorted(comps[0]) == [0, 1, 2]
    assert sorted(comps[1]) == [3, 4, 5]

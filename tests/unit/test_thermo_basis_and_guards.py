"""
tests/unit/test_thermo_basis_and_guards.py

Audit remediation (2026-08-26), deferred-by-design subset:

* Task 3 — Joback entropy basis in `KineticsEngine.get_reaction_thermo`.
  The audit flagged `s298 = (H298 - G298) / 298.15` (a *formation*-basis
  entropy) being combined with the compound's own `int Cp/T dT`. These tests
  pin down what that combination actually means: both the enthalpy and the
  entropy carried per species are element-referenced at 298.15 K, so the
  element reference cancels exactly for an atom-balanced reaction and the
  reaction quantities are correct. They also pin the failure mode that
  survives — atom-unbalanced steps — and the new metadata that exposes it.

* Task 2 — physics guards in `HeadspaceModel`: the van't Hoff extrapolation
  clamp and the `protein_fraction = 1.0` sentinel guard.
"""

import math

import pytest
from rdkit import Chem

from src.headspace import (
    VANT_HOFF_MAX_TEMP_K,
    VANT_HOFF_MIN_TEMP_K,
    HeadspaceModel,
)
from src.kinetics import (
    ELEMENT_STANDARD_ENTROPY_J_MOL_K,
    KineticsEngine,
    reaction_element_balance,
)
from src.thermo import JobackEstimator

# NIST / CODATA standard molar entropies of the ideal gas at 298.15 K, J/(mol K).
NIST_ABSOLUTE_ENTROPY_J_MOL_K = {
    "O": 188.84,   # H2O (g)
    "N": 192.77,   # NH3 (g)
    "S": 205.81,   # H2S (g)
}


def _element_counts(smiles: str) -> dict:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    counts: dict = {}
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return counts


def _element_entropy_sum(smiles: str) -> float:
    return sum(
        n * ELEMENT_STANDARD_ENTROPY_J_MOL_K[symbol]
        for symbol, n in _element_counts(smiles).items()
    )


# ---------------------------------------------------------------------------
# Task 3 — the Joback entropy really is a formation-basis entropy
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", ["O", "N", "S"])
def test_joback_s298_is_a_formation_basis_entropy(smiles):
    """
    Known reaction check: the formation reaction of each override species.

    (H298 - G298) / 298.15 must reproduce dfS(298) = S(compound) - sum_e n_e S(element),
    using NIST absolute entropies. This is the fact the whole basis argument
    rests on, so it is asserted directly rather than inferred.
    """
    est = JobackEstimator.estimate(smiles)
    s_formation = (est["H298"] - est["G298"]) / 298.15
    expected = NIST_ABSOLUTE_ENTROPY_J_MOL_K[smiles] - _element_entropy_sum(smiles)
    assert s_formation == pytest.approx(expected, abs=0.8)
    # And it is emphatically NOT the absolute entropy.
    assert abs(s_formation - NIST_ABSOLUTE_ENTROPY_J_MOL_K[smiles]) > 100.0


@pytest.mark.parametrize(
    "reactants, products",
    [
        (["CC=O", "N"], ["CC=N", "O"]),   # Schiff base + water: C, H, N, O
        (["CC=O", "S"], ["CC(S)O"]),      # thiol addition: brings S in
        (["CC=O"], ["C=CO"]),             # unimolecular enolisation
    ],
)
@pytest.mark.parametrize("temperature_k", [298.15, 423.15, 573.15])
def test_balanced_reaction_thermo_is_basis_consistent(reactants, products, temperature_k):
    """
    For an atom-balanced reaction, `get_reaction_thermo` must agree with an
    independent calculation that carries the *absolute* entropy of every
    species (formation entropy + the standard entropies of its elements).

    Agreement proves the element reference cancels: enthalpy and entropy are on
    one consistent basis, and the agreement does not degrade with temperature —
    which is the specific concern the audit raised.
    """
    engine = KineticsEngine()
    result = engine.get_reaction_thermo(reactants, products, temperature_k)
    assert result["atom_balanced"] is True
    assert result["thermo_basis"] == "formation_298K_element_referenced"

    def absolute_entropy(smiles: str) -> float:
        est = JobackEstimator.estimate(smiles)
        s_formation_298 = (est["H298"] - est["G298"]) / 298.15
        s_absolute_298 = s_formation_298 + _element_entropy_sum(smiles)
        a, b, c, d = est["cp_coeffs"]

        def int_cp_over_t(t):
            return a * math.log(t) + b * t + 0.5 * c * t**2 + (1 / 3) * d * t**3

        return s_absolute_298 + int_cp_over_t(temperature_k) - int_cp_over_t(298.15)

    delta_s_absolute = sum(absolute_entropy(s) for s in products) - sum(
        absolute_entropy(s) for s in reactants
    )
    # The absolute basis and the formation basis must give the same reaction dS.
    assert result["delta_s_j_mol_k"] == pytest.approx(delta_s_absolute, abs=1e-6)

    delta_g_reference = result["delta_h_j_mol"] - temperature_k * delta_s_absolute
    assert result["delta_g_j_mol"] == pytest.approx(delta_g_reference, rel=1e-9, abs=1e-6)


def test_balanced_reaction_values_are_unchanged_by_the_refactor():
    """
    Regression pins captured from the pre-refactor implementation. The rewrite
    is a pure algebraic rearrangement for balanced reactions, so these must not
    move — this is what makes the change safe to ship enabled.
    """
    engine = KineticsEngine()
    schiff = engine.get_reaction_thermo(["CC=O", "N"], ["CC=N", "O"], 423.15)
    assert schiff["delta_g_kcal_mol"] == pytest.approx(2.0208098, abs=1e-6)
    assert schiff["delta_h_kcal_mol"] == pytest.approx(4.2130557, abs=1e-6)

    enolisation = engine.get_reaction_thermo(["CC=O"], ["C=CO"], 423.15)
    assert enolisation["delta_g_kcal_mol"] == pytest.approx(25.1857, abs=1e-4)


# ---------------------------------------------------------------------------
# Task 3 — the real exposure: atom-unbalanced steps
# ---------------------------------------------------------------------------

def test_reaction_element_balance_detects_atom_and_charge_imbalance():
    assert reaction_element_balance(["CC=O", "N"], ["CC=N", "O"])["balanced"] is True

    dropped_water = reaction_element_balance(["CC=O", "N"], ["CC=N"])
    assert dropped_water["balanced"] is False
    assert dropped_water["element_imbalance"] == {"H": -2, "O": -1}

    # The one unbalanced SMIRKS step in the current benchmark panel
    # (thiamine Additive_Thermal_Degradation) is charge-unbalanced too.
    thiamine = reaction_element_balance(
        ["Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"], ["S", "Cc1cccs1", "CC1=NCCS1"]
    )
    assert thiamine["balanced"] is False
    assert thiamine["charge_imbalance"] == -1

    assert reaction_element_balance(["not-a-smiles"], ["O"])["balanced"] is None


def test_unbalanced_reaction_thermo_is_labelled_not_silently_trusted():
    """
    A dropped water leaves the formation enthalpy of the missing fragment in
    dG (tens of kcal/mol). The number is still returned by default, but it must
    be labelled so a caller can tell it is not a property of the chemistry.
    """
    engine = KineticsEngine()
    result = engine.get_reaction_thermo(["CC=O", "N"], ["CC=N"], 423.15)
    assert result["atom_balanced"] is False
    assert result["thermo_basis"] == "element_reference_uncancelled"
    assert result["element_imbalance"] == {"H": -2, "O": -1}
    # Dominated by the element residual, far above the 30 kcal/mol gate.
    assert result["delta_g_kcal_mol"] > 40.0
    assert result["delta_g_kcal_mol_raw"] == pytest.approx(result["delta_g_kcal_mol"])


def test_neutralize_flag_makes_unbalanced_steps_gate_neutral(monkeypatch):
    """The opt-in behaviour, kept off by default (see the module header)."""
    monkeypatch.setattr("src.kinetics.NEUTRALIZE_UNBALANCED_THERMO_GATE", True)
    engine = KineticsEngine()
    result = engine.get_reaction_thermo(["CC=O", "N"], ["CC=N"], 423.15)
    assert result["delta_g_kcal_mol"] == 0.0
    assert result["is_spontaneous"] is True
    assert result["thermo_basis"] == "unavailable_unbalanced"
    # The contaminated value is retained for forensics.
    assert result["delta_g_kcal_mol_raw"] > 40.0

    feasible, dg = engine.is_reaction_feasible(["CC=O", "N"], ["CC=N"], threshold_kcal_mol=30.0)
    assert feasible is True
    assert dg == 0.0


def test_estimator_failure_is_gate_neutral_and_self_describing():
    engine = KineticsEngine()
    result = engine.get_reaction_thermo(["definitely not a smiles"], ["O"], 423.15)
    assert result["delta_g_kcal_mol"] == 0.0
    assert result["is_spontaneous"] is True
    assert result["thermo_basis"] == "unavailable_estimator_error"
    # The pre-audit fallback dropped these keys entirely.
    assert "delta_h_kcal_mol" in result and "delta_s_j_mol_k" in result


# ---------------------------------------------------------------------------
# Task 2 — headspace guards
# ---------------------------------------------------------------------------

def test_vant_hoff_extrapolation_is_clamped_above_the_boiling_point():
    """
    Audit example: hexanal at 453 K extrapolated to Kaw ~ 3.7, i.e. a soluble
    aldehyde predicted to sit almost entirely in the vapour phase.
    """
    model = HeadspaceModel()
    kaw_453 = model.get_kaw_at_temp("Hexanal", 453.15)
    kaw_clamp = model.get_kaw_at_temp("Hexanal", VANT_HOFF_MAX_TEMP_K)

    assert kaw_453 == pytest.approx(kaw_clamp)
    assert kaw_453 < 1.0
    # Everything at or above the clamp collapses onto the same value.
    assert model.get_kaw_at_temp("Hexanal", 600.0) == pytest.approx(kaw_clamp)
    # ... and the cold end is clamped too.
    assert model.get_kaw_at_temp("Hexanal", 200.0) == pytest.approx(
        model.get_kaw_at_temp("Hexanal", VANT_HOFF_MIN_TEMP_K)
    )


def test_vant_hoff_inside_the_window_is_untouched_and_still_monotonic():
    model = HeadspaceModel()
    kaw_298 = model.get_kaw_at_temp("Hexanal", 298.15)
    kaw_333 = model.get_kaw_at_temp("Hexanal", 333.15)
    # At the reference temperature the tabulated value is returned untouched.
    assert kaw_298 == pytest.approx(model.data["Hexanal"]["Kaw_25c"])
    assert kaw_333 > kaw_298

    # No Kaw <= 1 ceiling is imposed: sparingly soluble gases legitimately
    # exceed 1 near the boiling point.
    assert model.get_kaw_at_temp("Hydrogen Sulfide", VANT_HOFF_MAX_TEMP_K) > 1.0


def test_protein_fraction_sentinel_is_neutralised_in_predict_headspace(monkeypatch):
    """
    ReactionConditions.protein_fraction defaults to 1.0 ("unspecified").
    Passed through untouched it would suppress a carbonyl by ~1/(1 + 45*1.0).
    """
    import src.headspace as headspace_module

    model = HeadspaceModel()
    matrix = {"Methional": 1.0}

    warnings: list = []
    monkeypatch.setattr(
        headspace_module.logger,
        "warning",
        lambda msg, *args, **kwargs: warnings.append(msg % args if args else msg),
    )

    sentinel = model.predict_headspace(matrix, 25.0, protein_fraction=1.0)
    unspecified = model.predict_headspace(matrix, 25.0, protein_fraction=0.0)

    assert sentinel["Methional"] == pytest.approx(unspecified["Methional"])
    assert len(warnings) == 1
    assert "protein_fraction" in warnings[0]
    # The un-guarded value would have been ~46x smaller for a carbonyl.
    assert unspecified["Methional"] > 0.0


def test_real_protein_fractions_still_suppress_headspace():
    """The guard must not swallow genuine compositions."""
    model = HeadspaceModel()
    matrix = {"Methional": 1.0}

    unspecified = model.predict_headspace(matrix, 25.0, protein_fraction=0.0)
    real = model.predict_headspace(matrix, 25.0, protein_fraction=0.2)

    # Methional Kprot = 2.0 -> 1 / (1 + 2.0 * 0.2)
    assert real["Methional"] == pytest.approx(unspecified["Methional"] / 1.4)


def test_negative_fractions_are_clamped_to_zero():
    model = HeadspaceModel()
    matrix = {"Hexanal": 1.0}
    assert model.predict_headspace(matrix, 25.0, fat_fraction=-0.5)[
        "Hexanal"
    ] == pytest.approx(model.predict_headspace(matrix, 25.0, fat_fraction=0.0)["Hexanal"])

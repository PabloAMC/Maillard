import pytest

from src.barrier_constants import (  # noqa: E402
    DEFAULT_BARRIER,
    DEFAULT_REFERENCE_PREEXPONENTIAL,
    arrhenius_rate_constant,
    effective_barrier_from_rate_constant,
    get_arrhenius_params,
    get_barrier,
    get_reference_pre_exponential,
)

def test_get_barrier_exact_match():
    # Test valid exact matches and normalized variants
    assert get_barrier("schiff_condensation")[0] == 15.0
    assert get_barrier("Schiff Condensation")[0] == 15.0 # Normalizes to schiff_condensation
    assert get_barrier("1,2-enolisation")[0] == 28.0
    assert get_barrier("1_2-enolisation")[0] == 28.0
    assert get_barrier("cysteine_thermolysis")[0] == 30.0
    
    # Test invalid or unknown matches fall back
    assert get_barrier("unknown_reaction")[0] == DEFAULT_BARRIER
    assert get_barrier("")[0] == DEFAULT_BARRIER

def test_get_arrhenius_params_valid():
    # Test that valid families from the YAML load correctly
    # Note: tests depend on the actual contents of arrhenius_params.yml
    
    # Schiff base
    params = get_arrhenius_params("schiff_condensation")
    assert params is not None
    A, Ea_kcal, _, _ = params
    assert A == 1.5e11
    # Audit 2026-08-26: corrected to the thesis value (Table 5.2, Glu+Gly->DFG).
    assert abs(Ea_kcal - 23.18) < 0.1 # 97.0 kJ/mol / 4.184 = 23.18 kcal/mol
    
    # Amadori
    params = get_arrhenius_params("amadori_rearrangement")
    # A_value is 1.0e11 in YAML
    assert params is not None
    A, Ea_kcal, quality, _ = params
    assert A == 1.0e11
    assert quality == "literature_estimated"
    
    # Dehydration (fructose)
    params = get_arrhenius_params("dehydration")
    assert params is not None
    A, Ea_kcal, _, _ = params
    assert A == 3.09e11
    assert abs(Ea_kcal - 29.6) < 0.1 # 123.88 kJ/mol / 4.184

def test_get_arrhenius_params_invalid():
    # Unknown families should return None
    assert get_arrhenius_params("magic_reaction") is None
    assert get_arrhenius_params("") is None

def test_get_arrhenius_params_normalization():
    # Should normalize "schiff" to "schiff_condensation"
    assert get_arrhenius_params("schiff") == get_arrhenius_params("schiff_condensation")


def test_enolisation_uses_enolisation_arrhenius_contract_not_fructose_dehydration():
    params = get_arrhenius_params("1,2-enolisation")
    assert params is not None
    A, ea_kcal, quality, _ = params

    assert A == 1e12
    assert quality == "literature_estimated"
    assert abs(ea_kcal - (120.0 / 4.184)) < 0.1


def test_fast_enolisation_barrier_stays_close_to_arrhenius_reference():
    fast_barrier, _ = get_barrier("1,2-enolisation")
    params = get_arrhenius_params("1,2-enolisation")
    assert params is not None
    _, arrhenius_barrier, _, _ = params

    assert abs(fast_barrier - arrhenius_barrier) < 1.0


def test_reference_pre_exponential_uses_literature_then_fallback():
    assert get_reference_pre_exponential("schiff_condensation") == 1.5e11
    assert get_reference_pre_exponential("magic_reaction") == DEFAULT_REFERENCE_PREEXPONENTIAL


def test_arrhenius_round_trip_preserves_family_barrier():
    barrier = 13.62
    temperature_kelvin = 423.15
    family = "schiff_condensation"

    rate_constant = arrhenius_rate_constant(barrier, temperature_kelvin, family=family)
    recovered = effective_barrier_from_rate_constant(rate_constant, temperature_kelvin, family=family)

    assert recovered == pytest.approx(barrier, rel=1e-9)


# ── Wave T4 (2026-08-27) ──────────────────────────────────────────────────────


def test_the_thiamine_arrhenius_anchor_is_reachable_from_the_emitted_family():
    """T4.2. The literature A/Ea pair existed and NO emitted family could reach it.

    `data/lit/arrhenius_params.yml` has carried a `thiamine_degradation` entry all
    along, but the engine's thiamine steps emit `Additive_Thermal_Degradation`
    (`src/reaction_templates.py:1943-1966`), which canonicalises to
    `additive_thermal_degradation` -- a FAST_BARRIERS key in its own right since
    Wave G1 fix 8. `_arrhenius_yaml_key` therefore returned None and
    `src/cantera_export.py` silently fell back to the heuristic prefactor. Same
    defect class as the eight DEFAULT_BARRIER fallthroughs, one lane over.
    """
    from src.barrier_constants import _arrhenius_yaml_key

    assert _arrhenius_yaml_key("Additive_Thermal_Degradation") == "thiamine_degradation"

    params = get_arrhenius_params("Additive_Thermal_Degradation")
    assert params is not None, "the anchor is unreachable again"
    a_value, ea_kcal, quality, _ = params
    assert a_value == 1.0e9
    assert ea_kcal == pytest.approx(100.8 / 4.184, abs=0.01)

    # THE ANCHOR IS UNSOUND AND THE FIX DOES NOT LAUNDER IT. That YAML entry is a
    # dead DOI with an unanchored Ea; wiring it makes the Cantera export REPORT
    # "dead_doi" where it previously reported "heuristic" for a step whose
    # literature entry existed the whole time. If this ever silently becomes
    # "literature_estimated" without a re-anchoring wave, that is laundering.
    assert quality == "dead_doi"


def test_wiring_the_thiamine_anchor_left_the_screening_lane_untouched():
    """SCOPE GUARD for T4.2: FAST_BARRIERS is authoritative for every prediction.

    Per the authority statement in `src/barrier_constants.py`, FAST_BARRIERS drives
    `evaluate_benchmark_payload` and the whole recommend lane, while the YAML drives
    the Cantera export only. T4.2 touched the YAML key map and nothing else, so no
    shipped number may move.
    """
    assert get_barrier("Additive_Thermal_Degradation")[0] == 25.0
    assert get_barrier("thiamine_degradation")[0] == 25.0
    assert get_barrier("additive_degradation")[0] == 25.0


def test_the_mutarotation_anchor_is_still_unreachable_and_that_is_recorded():
    """T4.3, pinned so the ledger item cannot be quietly lost.

    The `Sugar_Ring_Opening` / mutarotation lane cannot fire: every sugar
    `src/precursor_resolver.py` produces is already open-chain, so the cyclic
    hemiacetal the rule matches never exists. This test does NOT assert the lane
    should stay dead -- it asserts that the CURRENT state is the documented one, so
    that implementing cyclic sugar forms (the owner option filed as [P] under Wave
    T4) makes this test fail loudly and get updated, rather than passing unnoticed.
    """
    from src.barrier_constants import _arrhenius_yaml_key

    assert _arrhenius_yaml_key("Sugar_Ring_Opening") == "mutarotation"

    from src.precursor_resolver import resolve_many

    for sugar in ("D-Glucose", "D-Ribose", "D-Xylose", "D-Fructose"):
        (species,) = resolve_many([sugar])
        assert "1" not in species.smiles, (
            f"{sugar} now resolves to a RING form ({species.smiles}). The "
            f"mutarotation lane may now be live -- re-read the Wave T4 note on "
            f"`_sugar_ring_opening` in src/reaction_templates.py and update the "
            f"[P] item in tasks/audit_remediation.md."
        )

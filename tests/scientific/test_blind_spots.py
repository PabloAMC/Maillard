import pytest
import sys
from pathlib import Path
from src.pipeline import MaillardPipeline
from src.conditions import ReactionConditions

def test_peptide_accessibility_blind_spot():
    """
    BLIND SPOT: The engine assumes 100% amino acid availability.
    In raw pea protein, amino acids are buried. 
    EXPECTATION: High hydrolysis should increase flavor score.
    CURRENT: Score remains constant regardless of 'hydrolysis' state.
    """
    designer = MaillardPipeline(target_tag="meaty")
    cond = ReactionConditions(temperature_celsius=100)
    
    # Formulation with a hypothetical 'degree_of_hydrolysis'
    low_hydrolysis = {
        "name": "IntactProtein",
        "sugars": ["ribose"],
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"ribose": 1.0, "lysine": 1.0, "cysteine": 0.5},
        "degree_of_hydrolysis": 0.1 # 10%
    }
    
    high_hydrolysis = {
        "name": "HydrolyzedProtein",
        "sugars": ["ribose"],
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"ribose": 1.0, "lysine": 1.0, "cysteine": 0.5},
        "degree_of_hydrolysis": 0.9 # 90%
    }
    
    res_low = designer.evaluate_single(low_hydrolysis, cond)
    res_high = designer.evaluate_single(high_hydrolysis, cond)
    
    assert res_high.target_score > res_low.target_score

@pytest.mark.xfail(
    strict=True,
    reason=(
        "Engine does not yet support matrix inhibition (volatile partitioning by "
        "fibre/starch). Verified still absent 2026-08-27: clear and matrix conditions "
        "both return meaty radar 0.0."
    ),
)
def test_matrix_inhibition_blind_spot():
    """
    BLIND SPOT: Volatiles are 'trapped' by fiber/starch.
    EXPECTATION: High fiber content should decrease sensory radar scores.
    CURRENT: Radar scores depend only on chemical concentration.

    CONVERTED 2026-08-27 (Wave J2, red-team finding: self-excusing skips). The body used to
    call ``pytest.xfail("Engine does not yet support matrix inhibition ...")`` on the line
    BEFORE the assertion. ``pytest.xfail()`` is imperative: it aborts the test at that line.
    So the assertion below never executed, and the test would have reported xfail forever --
    including after someone implemented matrix inhibition. A marker that excuses the exact
    property under test, evaluated before that property is examined, can never be retired by
    the code getting better.

    The marker is now declarative and ``strict=True``, so the assertion RUNS. If the engine
    gains volatile partitioning, this test XPASSes, strict mode turns that into a FAILURE,
    and whoever implemented it is forced to come here and promote the test to a real one.
    That is the intended lifecycle for a known-gap test.

    NOTE FOR WHOEVER IMPLEMENTS THIS: the comparison is currently DEGENERATE, not merely
    negative. Both branches return meaty radar 0.0 at these conditions, so ``0.0 < 0.0`` is
    False and the xfail is earned partly by the baseline being zero. The first assertion
    below states that explicitly. Implementing matrix inhibition is not enough to flip this
    test -- the formulation also needs a non-zero meaty baseline to inhibit.
    """
    designer = MaillardPipeline(target_tag="meaty")
    cond_clear = ReactionConditions(protein_fraction=1.0) # Pure solution
    cond_matrix = ReactionConditions(protein_fraction=0.1, matrix_fiber=0.5) # High bread/pea matrix
    
    form = {
        "name": "BaseForm",
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"lysine": 1.0, "cysteine": 0.5}
    }
    
    res_clear = designer.evaluate_single(form, cond_clear)
    res_matrix = designer.evaluate_single(form, cond_matrix)

    # Stated first so the failure message names the real blocker. Measured 2026-08-27: 0.0.
    assert res_clear.radar["meaty"][0] > 0.0, (
        f"degenerate baseline: the protein-free control scores meaty "
        f"{res_clear.radar['meaty'][0]}, so there is nothing for a matrix to inhibit"
    )
    assert res_matrix.radar["meaty"][0] < res_clear.radar["meaty"][0]

@pytest.mark.xfail(
    strict=True,
    reason=(
        "Engine does not yet support non-heme metal catalysis (Fe2+). Verified still "
        "absent 2026-08-27: with and without metal_catalyst='Fe2+' the roasted target "
        "score is bit-identical at 0.11253319637981546."
    ),
)
def test_metal_catalysis_blind_spot():
    """
    BLIND SPOT: Non-heme iron (common in pea) accelerates pyrazines.
    EXPECTATION: Presence of Iron should lower pyrazine barriers.
    CURRENT: Only temperature and pH are taken into account.

    CONVERTED AND GIVEN AN ASSERTION 2026-08-27 (Wave J2, red-team finding: self-excusing
    skips). This was the worst of the family. The body computed ``res1`` and ``res2`` and
    then ended on ``pytest.xfail("Engine does not yet support ...")`` as its FINAL
    statement -- there was no assertion in the function at all. It named an expectation in
    its docstring, evaluated nothing against it, and reported a tidy xfail. Nothing about
    the engine could ever have changed that outcome.

    The expectation stated above is now written as an assertion, under a declarative
    ``strict=True`` marker so it actually runs. Implementing Fe2+ catalysis makes this
    XPASS, which strict mode turns into a failure, forcing promotion to a real test.
    """
    designer = MaillardPipeline(target_tag="roasted") # Pyrazines
    cond_no_iron = ReactionConditions(temperature_celsius=120)
    cond_iron = ReactionConditions(temperature_celsius=120, metal_catalyst="Fe2+")
    
    form = {
        "name": "NuttyForm",
        "amino_acids": ["glycine", "alanine"],
        "sugars": ["glucose"]
    }
    
    res1 = designer.evaluate_single(form, cond_no_iron)
    res2 = designer.evaluate_single(form, cond_iron)

    # The assertion this test always claimed to make and never did. Iron should accelerate
    # pyrazine formation, so the roasted score with Fe2+ must exceed the score without it.
    assert res2.target_score > res1.target_score, (
        f"metal_catalyst='Fe2+' changed nothing: roasted score {res2.target_score} with "
        f"iron vs {res1.target_score} without. ReactionConditions accepts the field but no "
        f"barrier consumes it."
    )
    assert res2.target_score > res1.target_score

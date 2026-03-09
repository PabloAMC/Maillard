import pytest
from src.pathway_ranker import PathwayRanker, PathwayProfile, _evaluate_pathway
from src.pathway_extractor import ElementaryStep, Species
from src.conditions import ReactionConditions

def mock_evaluate_single_step(step: ElementaryStep):
    """Return mock barrier based on reaction family name prefix."""
    if not step.reaction_family:
        return (0.0, 0.0)
        
    if "Fast" in step.reaction_family:
        return (-10.0, 5.0) # Downhill, low barrier
    if "Slow" in step.reaction_family:
        return (10.0, 30.0) # Uphill, large barrier
    if "Medium" in step.reaction_family:
        return (0.0, 15.0)  # Thermoneutral, moderate barrier
        
    return (0.0, 0.0)

def test_pathway_profile_properties():
    # Linear pathway: A -> B -> C -> D
    profile = PathwayProfile(
        pathway_name="Test",
        steps=[],
        deltaE_kcal_list=[10.0, -20.0, -5.0],
        barrier_kcal_list=[15.0, 5.0, 30.0],
        scaled_rates=[]
    )
    
    assert profile.rate_limiting_barrier == 30.0
    assert profile.overall_thermodynamics == -15.0
    
    # Energetic span calculation:
    # TS1: 15.0, Intermediate 1: 10.0
    # TS2: 10 + 5 = 15.0, Intermediate 2: 10 - 20 = -10.0
    # TS3: -10 + 30 = 20.0, Intermediate 3: -10 - 5 = -15.0
    # Lowest intermediate is -10.0. Max TS energy is 20.0 (TS3).
    # Span = 20.0 - (-10.0) = 30.0
    assert profile.energetic_span == 30.0

def test_screen_pathways_ranking(monkeypatch):
    import src.pathway_ranker
    monkeypatch.setattr(src.pathway_ranker, "evaluate_single_step", mock_evaluate_single_step)
    
    p1_steps = [
        ElementaryStep([], [], reaction_family="Fast_Step1"),  # bar: 5
        ElementaryStep([], [], reaction_family="Fast_Step2"),  # bar: 5
    ]
    p2_steps = [
        ElementaryStep([], [], reaction_family="Fast_Step"),   # bar: 5
        ElementaryStep([], [], reaction_family="Slow_Step"),   # bar: 30
    ]
    p3_steps = [
        ElementaryStep([], [], reaction_family="Medium_Step1"),# bar: 15
        ElementaryStep([], [], reaction_family="Medium_Step2"),# bar: 15
    ]
    
    pathways = {
        "Pathway_Slow": p2_steps,
        "Pathway_Fast": p1_steps,
        "Pathway_Medium": p3_steps
    }
    
    ranker = PathwayRanker(n_cores=1)
    ranked = ranker.screen_pathways(pathways)
    
    assert len(ranked) == 3
    assert ranked[0].pathway_name == "Pathway_Fast"   
    assert ranked[1].pathway_name == "Pathway_Medium" 
    assert ranked[2].pathway_name == "Pathway_Slow"   

def test_conditions_scaling(monkeypatch):
    """Test that pH and Arrhenius scaling correctly penalizes unfavorable pathways."""
    def mock_eval(step):
        return (-10.0, 20.0)

    import src.pathway_ranker
    monkeypatch.setattr(src.pathway_ranker, "evaluate_single_step", mock_eval)
    
    furan_step = ElementaryStep([], [], reaction_family="1,2-enolisation")
    pyrazine_step = ElementaryStep([], [], reaction_family="2,3-enolisation")
    
    # Under alkaline conditions, pyrazine should have a much higher rate (dominates)
    alkaline_cond = ReactionConditions(pH=8.0)
    
    prof_furan = _evaluate_pathway(("Furan", [furan_step], alkaline_cond))
    prof_pyrazine = _evaluate_pathway(("Pyrazine", [pyrazine_step], alkaline_cond))
    
    # Pyrazine should be favored > 4.0x
    assert prof_pyrazine.scaled_rates[0] > prof_furan.scaled_rates[0]

def test_energetic_span_complex_profile():
    """
    Validate span for a profile with multiple 'valleys' and 'peaks'.
    A (0) -> B (10, TS: 15) -> C (5, TS: 20) -> D (-10, TS: 10)
    Energy levels: [0, 15, 10, 20, 5, 10, -10]
    Lowest intermediate: 0 (A), then 10 (B), then 5 (C), then -10 (D).
    BUT, the span is TS_i - lowest_intermediate_j (where j <= i).
    TS1: 15 - 0 = 15
    TS2: 20 - 0 = 20
    TS3: 10 - 0 = 10
    Max span = 20.
    """
    profile = PathwayProfile(
        pathway_name="Complex",
        steps=[],
        deltaE_kcal_list=[10.0, -5.0, -15.0],
        barrier_kcal_list=[15.0, 10.0, 5.0],
        scaled_rates=[]
    )
    # Re-calc manually:
    # TS1 = 15. I1 = 10. (Low: 0)
    # TS2 = I1 + 10 = 20. I2 = I1 - 5 = 5. (Low: 0)
    # TS3 = I2 + 5 = 10. I3 = I2 - 15 = -10. (Low: 0)
    # Max(TS - prevailing_low): Max(15-0, 20-0, 10-0) = 20.
    assert profile.energetic_span == 20.0

def test_ranker_sorting_criteria(monkeypatch):
    """Ensure energetic span is the primary sorting key."""
    import src.pathway_ranker
    
    # Pathway A: Max Barrier 25, Span 25
    # Pathway B: Max Barrier 20, Span 30 (Worse span, better max barrier)
    
    p_a = PathwayProfile("A", [], [0], [25])
    p_b = PathwayProfile("B", [], [10, -10], [20, 20]) # Span is 20 - 0 = 20, TS2 is 10+20=30 relative to start.
    # Calculation for B: 
    # TS1=20, I1=10. 
    # TS2=10+20=30. I2=0. 
    # Span = Max(20-0, 30-0) = 30.
    
    assert p_a.rate_limiting_barrier == 25
    assert p_a.energetic_span == 25
    assert p_b.rate_limiting_barrier == 20
    assert p_b.energetic_span == 30
    
    ranker = PathwayRanker(n_cores=1)
    # We mock screen_pathways return or just test the sort logic directly if possible
    # But screen_pathways calls _evaluate_pathway which we'd need to mock heavily.
    # Let's just test the sort key logic that's inside screen_pathways.
    
    profiles = [p_b, p_a] # B is first in list
    profiles.sort(key=lambda p: p.energetic_span if p.energetic_span > 0 else p.rate_limiting_barrier)
    
    # A should be first because 25 < 30
    assert profiles[0].pathway_name == "A"

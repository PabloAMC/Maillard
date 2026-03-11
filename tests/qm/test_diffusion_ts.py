from src.diffusion_ts import DiffusionTSEngine  # noqa: E402

def test_diffusion_engine_initialization():
    engine = DiffusionTSEngine()
    assert engine is not None
    # In mock mode, we check availability
    assert hasattr(engine, 'available')

def test_diffusion_engine_confidence():
    engine = DiffusionTSEngine()
    # Maillard-like (N and O) should be high confidence
    score = engine.get_confidence_score("CC(N)C(=O)O", "CC(N)C=O")
    assert score >= 0.85
    
    # Non-Maillard (just C) should be lower
    score = engine.get_confidence_score("CCCC", "CCC")
    assert score < 0.85

def test_diffusion_engine_predict_mock():
    engine = DiffusionTSEngine()
    # Mock returns None until real weights are integrated
    res = engine.predict_ts_geometry("sugar", "amino")
    assert res is None or isinstance(res, str)

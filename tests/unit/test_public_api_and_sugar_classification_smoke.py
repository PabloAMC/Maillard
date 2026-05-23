
def test_public_api_imports_smoke():
    try:
        from src import MaillardPipeline, FormulationResult, ReactionConditions, SmirksEngine, FormulationOptimizer
        assert all(
            symbol is not None
            for symbol in (MaillardPipeline, FormulationResult, ReactionConditions, SmirksEngine, FormulationOptimizer)
        )
    except ImportError as exc:
        raise AssertionError(f"Public API imports failed: {exc}") from exc
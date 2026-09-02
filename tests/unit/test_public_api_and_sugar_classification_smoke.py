
def test_public_api_imports_smoke():
    try:
        # 2026-09-03 (retirement step B5b): the public API is the kinetic core's front door.
        from src import EnvelopeDeclaration, FormulationSpec, ProcessSpec, ThermalProgram, predict
        assert all(
            symbol is not None
            for symbol in (EnvelopeDeclaration, FormulationSpec, ProcessSpec, ThermalProgram, predict)
        )
    except ImportError as exc:
        raise AssertionError(f"Public API imports failed: {exc}") from exc
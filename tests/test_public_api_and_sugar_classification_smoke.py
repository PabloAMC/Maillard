import sys
from pathlib import Path


# Add src to path
sys.path.append(str(Path.cwd()))


def test_public_api_imports_smoke():
    try:
        from src import MaillardPipeline, FormulationResult, ReactionConditions, SmirksEngine, FormulationOptimizer
        assert all(
            symbol is not None
            for symbol in (MaillardPipeline, FormulationResult, ReactionConditions, SmirksEngine, FormulationOptimizer)
        )
    except ImportError as exc:
        raise AssertionError(f"Public API imports failed: {exc}") from exc


def test_sugar_classification_smoke():
    from src.pathway_extractor import Species
    from src.sugar_classifier import is_hexose, is_pentose, is_sugar

    glucose = Species(label="Glucose", smiles="C(C1C(C(C(C(O1)O)O)O)O)O")
    ribose = Species(label="Ribose", smiles="C1C(C(C(C(O1)O)O)O)O")

    assert is_sugar(glucose) is True
    assert is_hexose(glucose) is True
    assert is_pentose(ribose) is True


if __name__ == "__main__":
    test_public_api_imports_smoke()
    test_sugar_classification_smoke()
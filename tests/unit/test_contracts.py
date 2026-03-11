from src.xtb_screener import XTBScreener  # noqa: E402
from src.dft_refiner import DFTRefiner  # noqa: E402
from src.skala_refiner import SkalaRefiner  # noqa: E402

def test_screener_interface_contract():
    """Ensure XTBScreener implements the required pipeline methods."""
    # Classes themselves or instances should have these methods
    required_methods = [
        "generate_3d_xyz",
        "optimize_species",
        "compute_reaction_energy"
    ]
    for method in required_methods:
        assert hasattr(XTBScreener, method), f"XTBScreener missing {method}"

def test_refiner_interface_contract():
    """Ensure quantum refiners implement the required refinement interface."""
    required_methods = [
        "single_point",
        "calculate_barrier"
    ]
    for cls in [DFTRefiner, SkalaRefiner]:
        for method in required_methods:
            assert hasattr(cls, method), f"{cls.__name__} missing {method}"

def test_result_objects_contract():
    """Ensure result objects from screeners have mandatory fields."""
    # This catches changes to the XTBResult or ThermoResult dataclasses
    from src.xtb_screener import XTBResult
    
    # XTBResult check (used in pipeline)
    assert "energy_hartree" in XTBResult.__dataclass_fields__
    assert "optimized_xyz" in XTBResult.__dataclass_fields__

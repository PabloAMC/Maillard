import pytest
from src.pre_processor import PreProcessor

def test_yeast_cleaning():
    processor = PreProcessor()
    ratios = {"hexanal": 1.0, "glucose": 0.5}
    
    cleaned = processor.apply(ratios, ["yeast_fermentation"])
    
    assert cleaned["hexanal"] == pytest.approx(0.2)
    assert cleaned["hexanol"] == pytest.approx(0.8)
    assert cleaned["glucose"] == pytest.approx(0.5)


def test_protease_hydrolysis():
    processor = PreProcessor()
    ratios = {"leucine": 0.1, "ribose": 1.0}
    
    hydrolyzed = processor.apply(ratios, ["protease_hydrolysis"])
    
    assert hydrolyzed["leucine"] == pytest.approx(0.2)
    assert hydrolyzed["ribose"] == pytest.approx(1.0)


import numpy as np
from src.kinetics import KineticsEngine  # noqa: E402

def test_rate_constant_sanity():
    """Verify Eyring-Polanyi rate constant scales correctly with barrier."""
    engine = KineticsEngine(temperature_k=423.15) # 150 C
    
    # A 25 kcal/mol barrier at 150C should be ~1.08 s^-1
    k_25 = engine.get_rate_constant(25.0)
    assert k_25 < 2.0
    assert k_25 > 0.5
    
    # Lower barrier should be significantly faster
    k_15 = engine.get_rate_constant(15.0)
    assert k_15 > k_25 * 1000 # Exponential scaling

def test_half_life():
    """Verify half-life logic sanity."""
    engine = KineticsEngine(temperature_k=373.15) # 100 C
    k = engine.get_rate_constant(20.0)
    assert k > 0
    
    t_half = np.log(2) / k
    assert t_half > 0

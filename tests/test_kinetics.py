import unittest
from src.kinetics import KineticsEngine

class TestKinetics(unittest.TestCase):
    def test_rate_constant_sanity(self):
        engine = KineticsEngine(temperature_k=423.15) # 150 C
        
        # A 25 kcal/mol barrier at 150C is ~1.08 s^-1
        k_25 = engine.get_rate_constant(25.0)
        self.assertLess(k_25, 2.0)
        
        # A 15 kcal/mol barrier should be very fast
        k_15 = engine.get_rate_constant(15.0)
        self.assertGreater(k_15, k_25)
        
    def test_half_life(self):
        engine = KineticsEngine(temperature_k=373.15) # 100 C
        k = engine.get_rate_constant(20.0)
        # Check that we get a positive real number
        self.assertTrue(k > 0)

if __name__ == "__main__":
    unittest.main()

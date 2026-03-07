import unittest
import os
import sys
import numpy as np

# Add project root to path
sys.path.append(os.getcwd())

from src.dft_refiner import DFTRefiner

class TestIRC(unittest.TestCase):
    def setUp(self):
        self.refiner = DFTRefiner()
        # Use fast settings for test
        self.refiner.opt_basis = 'sto-3g'
        self.refiner.opt_method = 'hf'

    def test_irc_stub_simple(self):
        """
        Verify that generate_irc can run and identify an imaginary mode.
        We'll use a very simple hydrogen molecule at an 'activated' distance 
        or just a dummy system to check the plumbing.
        Actually, let's use water as a dummy 'TS' and expect it to fail 
        (no imaginary frequency) to test the error handling.
        """
        water_xyz = "3\n\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0"
        
        print("\nTesting IRC error handling (stable structure)...")
        with self.assertRaises(ValueError) as cm:
            self.refiner.generate_irc(water_xyz)
        self.assertIn("No imaginary frequency found", str(cm.exception))

if __name__ == "__main__":
    unittest.main()

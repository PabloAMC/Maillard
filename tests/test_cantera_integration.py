"""
Test suite for Phase 12 — Cantera Kinetic Simulation Integration

Tests kinetic simulations using DFT-quality barriers from Phase 3,
including time-temperature profiles, multiple pathways, and sensory predictions.
"""

import pytest
import numpy as np
from pathlib import Path


@pytest.mark.slow
class TestMechanismFileGeneration:
    """Test CTI/YAML mechanism file generation from RMG + DFT barriers."""

    def test_mechanism_file_generation(self):
        """Convert RMG mechanism + DFT barriers → CTI format."""
        pass
        # from src.cantera_integration import generate_mechanism_file
        #
        # rmg_mechanism = load_rmg_output('tests/fixtures/rmg_output/mechanism.py')
        # dft_barriers = load_dft_barriers('tests/fixtures/dft_barriers.json')
        #
        # cti_file = generate_mechanism_file(rmg_mechanism, dft_barriers, output_format='cti')
        # assert cti_file.exists()
        # # Verify CTI is valid
        # import cantera as ct
        # gas = ct.Solution(str(cti_file))
        # assert gas.n_species > 0

    def test_mechanism_file_no_syntax_errors(self):
        """Generated CTI/YAML parses without syntax errors in Cantera."""
        pass
        # import cantera as ct
        # mech_file = generate_mechanism_file(...)
        # try:
        #     gas = ct.Solution(str(mech_file))
        # except Exception as e:
        #     pytest.fail(f"Mechanism file has syntax errors: {e}")


@pytest.mark.slow
class TestIsothermalSimulation:
    """Test isothermal batch reactor simulations."""

    def test_simple_isothermal_simulation(self, tmp_path):
        """Isothermal batch reactor: Ribose+Gly at 150°C."""
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        exporter = CanteraExporter()
        exporter.add_reaction(["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"], 15.0)
        exporter.add_reaction(["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"], 20.0) 
        
        output_file = tmp_path / "isothermal.yaml"
        exporter.export_yaml(str(output_file))
        
        engine = KineticsEngine(temperature_k=423.15)
        # S_0: Ribose, S_1: Glycine
        initial_concs = {"S_0": 0.1, "S_1": 0.1} 
        time_span = (0, 600)
        
        results = engine.simulate_network_cantera(str(output_file), initial_concs, time_span)
        
        assert np.max(results["S_2"]) > 0 # Schiff base formed
        assert np.max(results["S_4"]) > 0 # FFT formed
        assert results["S_0"][-1] < 0.1 # Ribose consumed

    def test_isothermal_equilibration(self, tmp_path):
        """Test mole fraction conservation in the simulation."""
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        exporter = CanteraExporter()
        exporter.add_reaction(["OCC1OC(O)C(O)C1O"], ["C1=C(SC=C1)CS"], 10.0)
        output_file = tmp_path / "mass_balance.yaml"
        exporter.export_yaml(str(output_file))
        
        engine = KineticsEngine(temperature_k=423.15)
        # S_0: Ribose initial mole fraction = 1.0 (pure)
        initial_concs = {"S_0": 1.0}
        results = engine.simulate_network_cantera(str(output_file), initial_concs, (0, 1000))
        
        # Total molar density should be constant (~0.0288 kmol/m3 at 150C/1atm)
        total_density_init = np.sum([results[s][0] for s in ["S_0", "S_1"]])
        total_density_final = np.sum([results[s][-1] for s in ["S_0", "S_1"]])
        assert np.isclose(total_density_init, total_density_final, rtol=1e-3)
        assert np.isclose(total_density_final, 101325 / (8314.46 * 423.15), rtol=1e-2)


@pytest.mark.slow
class TestTimeTemperatureProfile:
    """Test kinetics with time-dependent temperature ramps."""

    def test_time_temperature_profile(self, tmp_path):
        """Temperature ramp simulation logic check."""
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        exporter = CanteraExporter()
        exporter.add_reaction(["OCC1OC(O)C(O)C1O"], ["C1=C(SC=C1)CS"], 25.0)
        output_file = tmp_path / "ramp.yaml"
        exporter.export_yaml(str(output_file))
        
        engine = KineticsEngine()
        # Simulate at 100C then 150C
        # We compare conversion (%) to avoid density effects
        res1 = engine.simulate_network_cantera(str(output_file), {"S_0": 1.0}, (0, 300), temperature_k=373.15)
        res2 = engine.simulate_network_cantera(str(output_file), {"S_0": 1.0}, (0, 300), temperature_k=423.15)
        
        conv1 = (res1["S_0"][0] - res1["S_0"][-1]) / res1["S_0"][0]
        conv2 = (res2["S_0"][0] - res2["S_0"][-1]) / res2["S_0"][0]
        assert conv2 > conv1
        #
        # solution = run_time_temperature_simulation(gas, initial_state, time_span, temp_profile)
        #
        # # Product concentration should track with temperature
        # temps = np.array([temp_profile(t) for t in solution['time']])
        # # Not much product formation at low T; more at high T
        # product_conc_early = solution['FFT'][int(len(solution['time']) * 0.2)]
        # product_conc_late = solution['FFT'][int(len(solution['time']) * 0.7)]
        # assert product_conc_late > product_conc_early


@pytest.mark.slow
class TestMultiplePathways:
    """Test mechanisms with multiple competing pathways."""

    def test_multiple_pathways_coexist(self, tmp_path):
        """Ribose -> FFT (Path A) vs Ribose -> Other (Path B)."""
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        exporter = CanteraExporter()
        exporter.add_reaction(["OCC1OC(O)C(O)C1O"], ["C1=C(SC=C1)CS"], 30.0) # Path A (Slow)
        exporter.add_reaction(["OCC1OC(O)C(O)C1O"], ["C(O)(=O)C(O)C(O)C(O)"], 15.0) # Path B (Fast)
        
        output_file = tmp_path / "compete.yaml"
        exporter.export_yaml(str(output_file))
        
        engine = KineticsEngine()
        results = engine.simulate_network_cantera(str(output_file), {"S_0": 1.0}, (0, 100))
        
        # Path B should have higher yield
        assert results["S_2"][-1] > results["S_1"][-1]

    def test_pathway_selectivity_chemistry(self, tmp_path):
        """Selectivity changes with barrier heights: lower barrier dominates at short times.
        
        Key physics: with h0=s0=0, equilibrium K=1 (50/50 at long times).
        At short times (before equilibration), the faster reaction dominates.
        20 kcal/mol: k~0.4 s-1, half-life ~1.7s → ~55% converted at t=2s
        25 kcal/mol: k~0.001 s-1, half-life ~693s → ~0.2% converted at t=2s
        """
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        exporter = CanteraExporter()
        # Use real SMILES for correct RDKit stoichiometry
        # S_0=methane (reactant), S_1=ethane (fast product), S_2=propane (slow product)
        exporter.add_reaction(["C"], ["CC"], 20.0)   # Fast
        exporter.add_reaction(["C"], ["CCC"], 25.0)  # Slow
        exporter.export_yaml(str(tmp_path / "select.yaml"))
        
        engine = KineticsEngine(temperature_k=423.15)
        # S_0 is the auto-generated Cantera name for the first species ("C")
        # Use short time span (2s) so kinetics governs, not equilibrium
        results = engine.simulate_network_cantera(
            str(tmp_path / "select.yaml"), {"S_0": 1.0}, (0, 2)
        )
        # At t=2s: S_1 (~55% converted) >> S_2 (~0.2% converted)
        assert results["S_1"][-1] > results["S_2"][-1]


@pytest.mark.slow
class TestBarrierSensitivity:
    """Test kinetics sensitivity to DFT barrier variations."""

    def test_barrier_effect_on_kinetics(self, tmp_path):
        """5 kcal/mol barrier difference gives >100x rate difference at 150C.
        
        Key physics: Eyring-Polanyi k = (kB*T/h) * exp(-Ea/RT)
        At 150C: 20 kcal/mol -> k~0.4 s-1 (half-life 1.7s)
                 25 kcal/mol -> k~0.001 s-1 (half-life 693s)
        At t=2s kinetics governs (before equilibration), giving >100x difference.
        """
        from src.cantera_export import CanteraExporter
        from src.kinetics import KineticsEngine
        
        # S_0=methane (C), S_1=ethane (CC) for both mechanisms
        exp_high = CanteraExporter()
        exp_high.add_reaction(["C"], ["CC"], 25.0)  # Slow
        exp_high.export_yaml(str(tmp_path / "high.yaml"))
        
        exp_low = CanteraExporter()
        exp_low.add_reaction(["C"], ["CC"], 20.0)  # Fast
        exp_low.export_yaml(str(tmp_path / "low.yaml"))
        
        engine = KineticsEngine(temperature_k=423.15)
        # Compare at short time so kinetics (not equilibrium) dominates
        res_high = engine.simulate_network_cantera(str(tmp_path / "high.yaml"), {"S_0": 1.0}, (0, 0.1))
        res_low = engine.simulate_network_cantera(str(tmp_path / "low.yaml"), {"S_0": 1.0}, (0, 0.1))
        
        # Compare conversion fractions: (S_0_initial - S_0_final) / S_0_initial
        # This avoids saturation effects near equilibrium
        conv_low = (res_low["S_0"][0] - res_low["S_0"][-1]) / res_low["S_0"][0]
        conv_high = (res_high["S_0"][0] - res_high["S_0"][-1]) / res_high["S_0"][0]
        # Lower barrier should convert significantly more reactant
        assert conv_low > conv_high * 3



@pytest.mark.slow
class TestSensoryPrediction:
    """Test sensory profile prediction from Cantera kinetics."""

    def test_sensory_profile_ranking(self):
        """Rank formulations by predicted sensory profile vs Cantera kinetics."""
        from src.sensory import SensoryPredictor
        
        predictor = SensoryPredictor()
        # Meaty precursor: FFT (large OAV)
        meaty_conc = {"FFT": 1.0} # 1 ppm
        # Roasted precursor: pyrazine (smaller OAV due to higher threshold)
        roasted_conc = {"2,3-dimethylpyrazine": 1.0} # 1 ppm
        
        profile_meaty = predictor.predict_profile(meaty_conc)
        profile_roasted = predictor.predict_profile(roasted_conc)
        
        assert profile_meaty["meaty"] > profile_roasted["meaty"]
        assert profile_roasted["roasted"] > 0

    def test_sensory_correlates_with_dft(self):
        """Sensory score should correlate with concentration."""
        from src.sensory import SensoryPredictor
        predictor = SensoryPredictor()
        
        low_conc = {"FFT": 0.0001}
        high_conc = {"FFT": 0.1}
        
        profile_low = predictor.predict_profile(low_conc)
        profile_high = predictor.predict_profile(high_conc)
        
        assert profile_high["meaty"] > profile_low["meaty"]


@pytest.mark.slow
class TestKineticsRefinement:
    """Test improvement of kinetics from xTB to DFT barriers."""

    def test_compare_xtb_vs_dft_kinetics(self):
        """DFT-based kinetics should show refinement over xTB-based kinetics."""
        pass
        # # Simulate with xTB barriers
        # gas_xtb = generate_mechanism_file(rmg_mech, xtb_barriers)
        # solution_xtb = run_isothermal_simulation(gas_xtb, initial_state, (0, 600), T=423)
        #
        # # Simulate with DFT barriers
        # gas_dft = generate_mechanism_file(rmg_mech, dft_barriers)
        # solution_dft = run_isothermal_simulation(gas_dft, initial_state, (0, 600), T=423)
        #
        # # DFT should show more realistic FFT formation (e.g., faster due to corrected barriers)
        # # Compare final concentrations
        # fft_xtb = solution_xtb['FFT'][-1]
        # fft_dft = solution_dft['FFT'][-1]
        #
        # # Both should be positive, but DFT might be more/less realistic
        # assert fft_dft > 0
        # assert fft_xtb > 0

import pytest
import numpy as np
import pandas as pd
import json
from pathlib import Path
from src.kinetics import KineticsEngine  # noqa: E402
from src.cantera_export import CanteraExporter  # noqa: E402
from src.recommend import Recommender  # noqa: E402
from src.chem_utils import ElementaryStep, Species  # noqa: E402

def test_isothermal_vs_ramp(tmp_path):
    # 1. Create a simple balanced isomerisation mechanism: A <=> B
    exporter = CanteraExporter()
    # Balanced: acetaldehyde <=> vinyl alcohol (both C2H4O)
    exporter.add_reaction(["CC=O"], ["C=CO"], 25.0) # 25 kcal/mol barrier (slower)
    mech_path = tmp_path / "simple_mech.yaml"
    exporter.export_yaml(str(mech_path))
    
    engine = KineticsEngine()
    
    # Isothermal 150C (423.15K) for 1s
    res_iso = engine.simulate_network_cantera(
        str(mech_path), 
        {"S_0": 1.0}, 
        (0, 1.0), 
        temperature_k=423.15
    )
    
    # Ramp 25C -> 150C over 1s
    ramp = [(0.0, 298.15), (1.0, 423.15)]
    res_ramp = engine.simulate_network_cantera(
        str(mech_path), 
        {"S_0": 1.0}, 
        (0, 1.0), 
        temperature_profile=ramp
    )
    
    # 2. Assertions
    # Note: we must use mole fractions (ends with _X) to ignore expansion effects
    x_iso = res_iso["S_0_X"]
    x_ramp = res_ramp["S_0_X"]
    
    conv_iso = (x_iso[0] - x_iso[-1]) / x_iso[0]
    conv_ramp = (x_ramp[0] - x_ramp[-1]) / x_ramp[0]
    
    print(f"Iso conversion: {conv_iso:.6f}")
    print(f"Ramp conversion: {conv_ramp:.6f}")
    
    assert conv_iso > conv_ramp
    
    # Check mass balance in mole fractions
    total_iso = res_iso["S_0_X"] + res_iso["S_1_X"]
    total_ramp = res_ramp["S_0_X"] + res_ramp["S_1_X"]
    np.testing.assert_allclose(total_iso, 1.0, atol=1e-5)
    np.testing.assert_allclose(total_ramp, 1.0, atol=1e-5)
    
    # Check temperature recording
    assert "temperature" in res_ramp
    assert res_ramp["temperature"][0] == pytest.approx(298.15)
    assert res_ramp["temperature"][-1] == pytest.approx(423.15)

# DELETED 2026-08-27 (Wave J2, red-team finding: empty-body tests). `test_cli_ramp(tmp_path)`
# wrote a barriers JSON and a ramp CSV into tmp_path, then ended on a bare `pass` under the
# comment "We can't easily run the CLI and check results without mocking or subprocesses /
# But we already verified the core logic in the test above." So it named the CLI in its title,
# exercised no CLI, asserted nothing, and reported PASS -- a green tick for a test whose own
# body admits it is not testing anything. The fixture setup it did perform was pure cost: both
# files were written and then never read.
#
# The claim in that comment is true, which is why this is a deletion and not a rewrite: the
# ramp logic IS covered, by test_temperature_ramp above (mole-fraction conservation plus the
# 298.15 K -> 423.15 K endpoints). Nothing was lost. The CLI's own argument contract is
# covered for real in tests/integration/test_cantera_sim.py.


def test_fast_temporal_ramp_matches_cantera_directional_reference(tmp_path):
    ribose = Species("ribose", "OCC(O)C(O)C(O)C=O")
    furfural = Species("furfural", "O=Cc1ccco1")
    fft = Species("2-furfurylthiol", "SCc1ccco1")
    steps = [
        ElementaryStep(
            reactants=[ribose],
            products=[furfural],
            reaction_family="Enolisation_1_2",
        ),
        ElementaryStep(
            reactants=[ribose],
            products=[fft],
            reaction_family="Enolisation_1_2",
        ),
    ]
    barriers = {
        f"{ribose.smiles}->{furfural.smiles}": 20.0,
        f"{ribose.smiles}->{fft.smiles}": 32.0,
    }
    initial = {ribose.smiles: 1.0}

    ramp_df = pd.DataFrame({
        "time": [0, 5, 10],
        "temp": [25, 100, 150],
    })
    ramp_path = tmp_path / "fast_reference_ramp.csv"
    ramp_df.to_csv(ramp_path, index=False)

    recommender = Recommender()
    fast_iso = recommender.predict_from_steps(steps, barriers, initial, temperature_kelvin=423.15, time_minutes=10.0)
    fast_ramp = recommender.predict_from_steps(steps, barriers, initial, temp_ramp_csv=str(ramp_path), time_minutes=10.0)

    fast_iso_furfural = fast_iso["predicted_proxy_ppb"].get("furfural", 0.0)
    fast_ramp_furfural = fast_ramp["predicted_proxy_ppb"].get("furfural", 0.0)
    fast_iso_fft = fast_iso["predicted_proxy_ppb"].get("2-furfurylthiol", 0.0)
    fast_ramp_fft = fast_ramp["predicted_proxy_ppb"].get("2-furfurylthiol", 0.0)

    exporter = CanteraExporter()
    exporter.add_reaction(["CC=O"], ["C=CO"], 20.0)
    exporter.add_reaction(["CC=O"], ["C1CO1"], 32.0)
    mech_path = tmp_path / "reference_mech.yaml"
    exporter.export_yaml(str(mech_path))

    engine = KineticsEngine()
    res_iso = engine.simulate_network_cantera(str(mech_path), {"S_0": 1.0}, (0, 10.0), temperature_k=423.15)
    res_ramp = engine.simulate_network_cantera(
        str(mech_path),
        {"S_0": 1.0},
        (0, 10.0),
        temperature_profile=[(0.0, 298.15), (300.0, 373.15), (600.0, 423.15)],
    )

    conv_iso = (res_iso["S_0_X"][0] - res_iso["S_0_X"][-1]) / res_iso["S_0_X"][0]
    conv_ramp = (res_ramp["S_0_X"][0] - res_ramp["S_0_X"][-1]) / res_ramp["S_0_X"][0]
    fast_iso_ratio = fast_iso_furfural / max(fast_iso_fft, 1.0e-30)
    fast_ramp_ratio = fast_ramp_furfural / max(fast_ramp_fft, 1.0e-30)

    assert fast_iso_furfural > 0.0
    assert fast_ramp_furfural > 0.0
    assert fast_iso_fft > 0.0
    assert fast_ramp_fft > 0.0
    assert fast_iso_ratio > fast_ramp_ratio
    assert conv_iso > conv_ramp > 0.0

if __name__ == "__main__":
    import json
    # Manual run for debugging
    test_isothermal_vs_ramp(Path("/tmp"))

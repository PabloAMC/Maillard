import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from scripts.run_cantera_kinetics import run_simulation

precursors = {"ribose": 0.1, "cysteine": 0.05, "leucine": 0.05}
db_path = "results/maillard_results.db"

out_prefix = "debug_fft"

run_simulation(
    barriers_json=db_path,
    precursors=precursors,
    temp_c=150.0,
    time_sec=3600.0,
    from_smirks=True,
    track_species=["2-furfurylthiol", "H2S"],
    output_prefix=out_prefix,
    verbose_reactions=False,
    no_gating=True,
    pH=6.0,
    solvent="water"
)

import json
import pathlib
import re

# Paths
ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCHMARK_MD = ROOT / "data" / "benchmarks" / "maillard_validation_benchmarks.md"
DEST_DIR = ROOT / "data" / "benchmarks"

def extract_hofmann_1998():
    """Extract Hofmann & Schieberle (1998) quantitative benchmark."""
    # Data from section 1.3
    bench = {
        "benchmark_id": "cys_ribose_140C_Hofmann1998",
        "source_doi": "10.1021/jf9705983",
        "precursors": {
            "L-Cysteine": {"concentration_mM": 10.0},
            "D-Ribose": {"concentration_mM": 10.0}
        },
        "conditions": {
            "temp_C": 140,
            "ph": 5.0,
            "water_activity": 0.98,
            "time_min": 30
        },
        "measured_volatiles": {
            # 0.03 mol% of 10mM = 0.003 mM = 3 uM = 3000 ppb (approx)
            # Actually MFT is MW 114. 3 uM * 114 = 342 ug/L = 342 ppb
            "2-methyl-3-furanthiol": {"conc_ppb": 342, "uncertainty_pct": 20},
            "2-furfurylthiol": {"conc_ppb": 200, "uncertainty_pct": 25}
        }
    }
    return bench

def extract_pratap_singh_2021_pea():
    """Extract Pratap-Singh et al. (2021) pea isolate baseline."""
    # Data from section 3.1
    bench = {
        "benchmark_id": "pea_isolate_40C_PratapSingh2021",
        "source_doi": "10.3390/molecules26134104",
        "precursors": {
            "Pea Protein Isolate": {"concentration_mM": 1000.0} # Representing matrix
        },
        "conditions": {
            "temp_C": 40,
            "ph": 6.0,
            "water_activity": 0.95,
            "time_min": 10
        },
        "measured_volatiles": {
            "hexanal": {"conc_ppb": 260, "uncertainty_pct": 14},
            "2-pentylfuran": {"conc_ppb": 638, "uncertainty_pct": 8},
            "hexanol": {"conc_ppb": 80, "uncertainty_pct": 19}
        },
        "protein_type": "pea_iso",
        "denaturation_state": 0.0 # Raw isolate
    }
    return bench

def main():
    benchmarks = [
        extract_hofmann_1998(),
        extract_pratap_singh_2021_pea()
    ]
    
    for b in benchmarks:
        outfile = DEST_DIR / f"{b['benchmark_id']}.json"
        print(f"Writing {outfile}...")
        with open(outfile, "w") as f:
            json.dump(b, f, indent=2)

if __name__ == "__main__":
    main()

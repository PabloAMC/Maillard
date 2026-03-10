"""
scripts/mlp_benchmark.py

Phase B.3: Benchmark MACE-OFF against existing xTB pathway endpoints.
"""

import os
import json
from pathlib import Path
from src.mlp_barrier import MLPBarrier

def run_benchmark():
    input_dir = Path("data/geometries/xtb_inputs")
    results_file = Path("results/mlp_benchmark.json")
    results_file.parent.mkdir(parents=True, exist_ok=True)
    
    mlp = MLPBarrier(model="medium", device="cpu")
    
    benchmark_results = {}
    
    for pathway_dir in input_dir.iterdir():
        if not pathway_dir.is_dir():
            continue
            
        reactant_file = pathway_dir / "reactant.xyz"
        product_file = pathway_dir / "product.xyz"
        
        if not (reactant_file.exists() and product_file.exists()):
            continue
            
        print(f"Processing pathway: {pathway_dir.name}...")
        
        reactant_xyz = reactant_file.read_text()
        product_xyz = product_file.read_text()
        
        try:
            e_reac = mlp.get_energy(reactant_xyz)
            e_prod = mlp.get_energy(product_xyz)
            barrier_est = mlp.estimate_barrier(reactant_xyz, product_xyz)
            
            benchmark_results[pathway_dir.name] = {
                "reactant_energy_ev": e_reac,
                "product_energy_ev": e_prod,
                "estimated_barrier_kcal": barrier_est
            }
        except Exception as e:
            print(f"  Error processing {pathway_dir.name}: {e}")
            
    with open(results_file, "w") as f:
        json.dump(benchmark_results, f, indent=2)
        
    print(f"\nBenchmark complete. Results saved to {results_file}")

if __name__ == "__main__":
    run_benchmark()

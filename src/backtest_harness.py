import os
import json
import sys
from pathlib import Path

# Ensure we can import from src
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.benchmark_validation import get_benchmark_files, evaluate_benchmark, summarize_evaluation

def run_backtest():
    print("Starting Computational Gap Closure Backtest...")
    files = get_benchmark_files()
    
    masked_offsets = {
        "schiff_condensation": 3.58,
        "amadori": 3.58,
        "enol": 3.58,
        "strecker": 3.58,
        "cys": 3.58
    }
    
    # Run Baseline (Unmasked)
    print("Running baseline (unmasked) benchmarks...")
    if "BARRIER_OFFSETS" in os.environ:
        del os.environ["BARRIER_OFFSETS"]
        
    baseline_results = {}
    for f in files:
        try:
            eval_res = evaluate_benchmark(f)
            summary = summarize_evaluation(eval_res)
            if isinstance(summary, tuple):
                summary = summary[0]
            baseline_results[f.name] = summary
        except Exception as e:
            print(f"Baseline error {f}: {e}")

    # Run Masked (Computational Surrogate Only)
    print("Running masked (computational only) benchmarks...")
    os.environ["BARRIER_OFFSETS"] = json.dumps(masked_offsets)
    
    masked_results = {}
    for f in files:
        try:
            eval_res = evaluate_benchmark(f)
            summary = summarize_evaluation(eval_res)
            if isinstance(summary, tuple):
                summary = summary[0]
            masked_results[f.name] = summary
        except Exception as e:
            print(f"Masked error {f}: {e}")
            
    # Remove env var
    if "BARRIER_OFFSETS" in os.environ:
        del os.environ["BARRIER_OFFSETS"]

    # Scorecard Evaluation
    scorecard = {
        "generated_at": "2026-03-27",
        "description": "Backtest scorecard demonstrating recovery of ranking/direction when empirical anchors are masked with +15 kJ/mol computational uncertainty.",
        "lanes": []
    }
    
    total_evals = 0
    preserved_ranking = 0
    preserved_ratio = 0
    
    for f in files:
        name = f.name
        if name not in baseline_results or name not in masked_results:
            continue
            
        b_res = baseline_results[name]
        m_res = masked_results[name]
        
        total_evals += 1
        r1 = getattr(b_res, "ranking_status", "pass")
        r2 = getattr(m_res, "ranking_status", "fail")
        if r1 == r2:
            preserved_ranking += 1
            
        s1 = getattr(b_res, "scale_status", "pass")
        s2 = getattr(m_res, "scale_status", "fail")
        if s1 == s2:
            preserved_ratio += 1 

    metrics = {
        "total_benchmarks_evaluated": total_evals,
        "ranking_preservation_rate": preserved_ranking / max(total_evals, 1),
        "ratio_band_recovery_rate": preserved_ratio / max(total_evals, 1),
        "recommendation_stability": "high" if (preserved_ranking / max(total_evals, 1)) > 0.8 else "low"
    }
    
    scorecard["overall_metrics"] = metrics
    
    out_dir = Path("results/validation")
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "computational_gap_closure_scorecard.json"
    md_path = out_dir / "computational_gap_closure_scorecard.md"
    
    with open(json_path, "w") as f:
        json.dump(scorecard, f, indent=2)
        
    with open(md_path, "w") as f:
        f.write("# Computational Gap Closure Scorecard\n\n")
        f.write(f"- **Ranking Preservation Rate:** {metrics['ranking_preservation_rate']:.1%}\n")
        f.write(f"- **Ratio-Band Recovery Rate:** {metrics['ratio_band_recovery_rate']:.1%}\n")
        f.write(f"- **Recommendation Stability:** {metrics['recommendation_stability'].upper()}\n\n")
        f.write("This artifact proves whether computational closure lanes earn their promotion by actually backtesting them against masked literature anchors.\n")
        
    print(f"Backtest complete. Preserved {preserved_ranking}/{total_evals} benchmarks.")
    print(f"Scorecard saved to {md_path}")

if __name__ == "__main__":
    run_backtest()


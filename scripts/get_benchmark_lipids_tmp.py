from src.benchmark_validation import evaluate_benchmark
import json

def get_values(path):
    eval_res = evaluate_benchmark(path)
    preds = eval_res.predicted_ppb
    print(f"--- {path} ---")
    print(f"Hexanal: {preds.get('Hexanal')}")
    print(f"Nonanal: {preds.get('Nonanal')}")
    print()

get_values("data/benchmarks/pea_isolate_ribose_cysteine_100C_45min_Internal2026.json")
get_values("data/benchmarks/soy_isolate_ribose_cysteine_100C_45min_Internal2026.json")

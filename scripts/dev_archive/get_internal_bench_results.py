
import json
from pathlib import Path
from src.benchmark_validation import evaluate_benchmark

bench_file = "data/benchmarks/pea_isolate_ribose_cysteine_100C_45min_Internal2026.json"
eval_res = evaluate_benchmark(bench_file)

print(json.dumps(eval_res.predicted_ppb, indent=2))

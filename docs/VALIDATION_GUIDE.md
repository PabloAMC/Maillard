# Scientific Validation Guide

## 1. Validation Contract

The project currently uses a layered benchmark contract rather than a single release number.

- **Coverage**: all measured compounds in a supported benchmark must resolve to predicted outputs.
- **Ranking**: Pearson R is only reported when at least 3 compounds match.
- **Scale**: the strict free-amino-acid gate requires full coverage plus a max measured/predicted ratio <= 1.5x.
- **Support status**: benchmarks that cannot run through the current execution path are reported as `unsupported`, not silently skipped.

The current benchmark summary artifact is generated with:

```bash
python scripts/generate_benchmark_summary.py
```

This writes `results/validation/benchmark_summary.md` and `results/validation/benchmark_summary.json` and is the default diagnostic surface for supported vs unsupported cases, ranking quality, and remaining scale gaps.

The explicit strict gate for free-amino-acid benchmarks is:

```bash
MAILLARD_STRICT_BENCHMARKS=1 python -m pytest tests/scientific/test_benchmarks.py -q
```

## 2. Current Validated Envelope

As of the current Docker-validated benchmark summary:

- `cys_glucose_150C_Farmer1999` is `pass` and `strict-ready`.
- `cys_ribose_140C_Hofmann1998` remains `scale-gap`.
- `cys_ribose_150C_Mottram1994` remains `scale-gap`.
- `pea_isolate_40C_PratapSingh2021` remains `unsupported` because it is matrix-only and cannot yet run through the free-precursor path.

This means the framework is now honest about its validated free-amino-acid envelope, but its absolute concentration projection is still not a generally validated headspace model.

## 3. How We Validate

### A. Literature Benchmarks

Benchmark execution is centralized in `src/benchmark_validation.py` and reused by:

- `tests/scientific/test_benchmarks.py`
- `tests/scientific/test_free_aa_quantitative_regression.py`
- `scripts/compare_sim_to_lit.py`
- `scripts/generate_benchmark_summary.py`

### B. Sensory Fidelity

Sensory scoring remains downstream of concentration prediction. Stevens' law and masking logic are useful for ranking perceived aroma, but they do not replace benchmark validation of predicted ppb.

### C. Safety Conservatism

Safety scoring is still treated conservatively and independently from the benchmark gate so that flavor optimization does not hide formation-risk regressions.

## 4. Blind Spots That Still Matter

- **Matrix-only systems**: plant-isolate benchmarks still require a dedicated matrix-aware execution path.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.

## 5. Recommended Verification Workflow

1. Run `python scripts/generate_benchmark_summary.py` to inspect the current validated envelope.
2. Run `python -m pytest tests/scientific/test_free_aa_quantitative_regression.py -q` to guard the currently calibrated free-amino-acid ratios.
3. Run `MAILLARD_STRICT_BENCHMARKS=1 python -m pytest tests/scientific/test_benchmarks.py -q` when you want an honest go/no-go signal for the strict free-AA gate.
4. Treat any unsupported matrix benchmark as a scope gap, not as a passing scientific result.

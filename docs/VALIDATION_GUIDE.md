# Scientific Validation Guide

## 1. Validation Contract

The project currently uses a layered benchmark contract rather than a single release number.

- **Coverage**: all measured compounds in a supported benchmark must resolve to predicted outputs.
- **Ranking**: Pearson R is only reported when at least 3 compounds match.
- **Scale**: the strict free-amino-acid gate requires full coverage plus a max measured/predicted ratio <= 1.5x.
- **Support status**: benchmarks that cannot run through the current execution path are reported as `unsupported`, not silently skipped.

The current benchmark summary artifact is generated with the Docker-validated helper:

```bash
./scripts/docker_maillard.sh summary
```

This writes `results/validation/benchmark_summary.md` and `results/validation/benchmark_summary.json` and is the default diagnostic surface for supported vs unsupported cases, ranking quality, and remaining scale gaps.

The benchmark index artifact is generated with:

```bash
./scripts/docker_maillard.sh index
```

This writes `results/validation/benchmark_index.md` and `results/validation/benchmark_index.json` so execution-path limits such as matrix-only benchmarks are explicit.

The explicit strict gate for free-amino-acid benchmarks is:

```bash
MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q
```

## 2. Current Validated Envelope

As of the current Docker-validated benchmark summary:

- `cys_glucose_150C_Farmer1999` is `pass` and `strict-ready`.
- `cys_ribose_140C_Hofmann1998` is `partial-pass` and `strict-ready`.
- `cys_ribose_150C_Mottram1994` is `pass` and `strict-ready`.
- `pea_isolate_40C_PratapSingh2021` remains `unsupported` because it is matrix-only and cannot yet run through the free-precursor path.

This means the framework now has three strict-ready free-amino-acid benchmarks and one explicit matrix-only scope gap (`pea isolate`). The absolute concentration projection is therefore materially better calibrated for the current free-precursor envelope, but it is still not a generally validated matrix headspace model.

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

### D. Execution Lanes

The validated execution contract now uses named Docker lanes instead of ad hoc command sequences:

- `./scripts/docker_maillard.sh core`: unit and integration correctness checks.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index generation plus scientific FAST regressions.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend validation.
- `./scripts/docker_maillard.sh hofmann`: reproducible diagnostic trace for the calibrated Hofmann sulfur benchmark.

## 4. Blind Spots That Still Matter

- **Matrix-only systems**: plant-isolate benchmarks still require a dedicated matrix-aware execution path.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.

## 5. Recommended Verification Workflow

1. For macOS/OrbStack or Docker Desktop, bring up the validated Linux environment with `./scripts/docker_maillard.sh up` and create or refresh the env with `./scripts/docker_maillard.sh bootstrap`.
2. Run `./scripts/docker_maillard.sh summary` to inspect the current validated envelope.
3. Run `./scripts/docker_maillard.sh index` to inspect benchmark tiering and execution-path scope, especially for matrix-only benchmarks.
4. Run `./scripts/docker_maillard.sh pytest tests/scientific/test_free_aa_quantitative_regression.py -q` to guard the currently calibrated free-amino-acid ratios.
5. Run `./scripts/docker_maillard.sh pytest tests/` for the full Docker-validated suite.
6. Run `MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q` when you want an honest go/no-go signal for the strict free-AA gate.
7. Treat any unsupported matrix benchmark as a scope gap, not as a passing scientific result.

## 6. Expected Skips In The Docker Lane

- `tests/benchmarks/` is intentionally skip-heavy today. Those tests are Phase 3 placeholders and HPC-oriented literature checks, not part of the current release gate.
- Capability-gated QM tests should skip only when the backend is genuinely unavailable or unusable in the active Docker environment.
- A path-based skip for a binary that is actually present in `PATH` is a test bug, not a valid environment gate.

## 7. Named Docker Lanes

- `./scripts/docker_maillard.sh core`: unit and integration correctness gate.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index plus scientific regression lane.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend lane.
- `./scripts/docker_maillard.sh hofmann`: diagnostic trace for the calibrated Hofmann sulfur benchmark.

These lanes keep infrastructure validation separate from unresolved scientific scope gaps such as matrix-only precursor systems.

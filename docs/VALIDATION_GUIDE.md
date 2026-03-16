# Scientific Validation Guide

## 1. Validation Contract

The project currently uses a layered benchmark contract rather than a single release number.

- The single source of truth for the benchmark contract now lives in `src/validation_contract.py`.
- **Coverage**: all measured compounds in a supported benchmark must resolve to predicted outputs.
- **Ranking**: Pearson R is only reported when at least 3 compounds match.
- **Scale**: the strict free-amino-acid gate requires full coverage plus a max measured/predicted ratio <= 1.5x.
- **Support status**: benchmarks that cannot run through the current execution path are reported as `unsupported`, not silently skipped.

The contract now separates three meanings that had started to drift together in docs and tests:

- **Directional validity**: preserves the reported ordering or trend without claiming absolute concentration fidelity.
- **Quantitative replication**: meets full-coverage plus thresholded Pearson / ratio criteria.
- **Formulation utility**: remains useful for recipe ranking or intervention comparison even when a benchmark is not strict-ready.

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

For reproducible target-level scientific inspection, use:

```bash
./scripts/docker_maillard.sh targets data/benchmarks/cys_ribose_140C_Hofmann1998.json
./scripts/docker_maillard.sh targets data/benchmarks/cys_glucose_150C_Farmer1999.json competing
```

This exposes the current target snapshot (`ppb`, `span`, `depth`, weighted flux) without relying on brittle inline Python or shell quoting.

For the aggregated target artifact with headspace observability metadata, use:

```bash
./scripts/docker_maillard.sh targets-report
```

This writes `results/validation/benchmark_targets.md` and `results/validation/benchmark_targets.json`, including per-target headspace class and Henry-law metadata when available.

`matrix_only` benchmarks are intentionally excluded from this target artifact. They remain executable through the summary and index artefacts, but they do not run through the free-precursor FAST target-ranking path, so exposing them as ordinary target snapshots would be scientifically misleading.

The explicit strict gate for free-amino-acid benchmarks is:

```bash
MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q
```

## 2. Current Validated Envelope

As of the current Docker-validated benchmark summary:

- `cys_glucose_150C_Farmer1999` is `pass` and `strict-ready`.
- `cys_ribose_140C_Hofmann1998` is `partial-pass` and `strict-ready`.
- `cys_ribose_150C_Mottram1994` is `pass` and `strict-ready`.
- `pea_isolate_40C_PratapSingh2021` is now executable through a dedicated `matrix_only` intake path, with full coverage, but it remains outside the strict release gate.
- `soy_isolate_40C_PratapSingh2021` is now executable through the same dedicated `matrix_only` intake path, with full coverage, but it also remains outside the strict release gate.

The headspace module now also carries a conservative matrix fallback for `pea_iso` and `soy_iso` when no explicit fat/protein fractions are supplied. That fallback is intentionally tied to the same retention estimates already used in `src/matrix_correction.py`, and those retention estimates now relax with `denaturation_state` in the same way across matrix correction, output projection, and sensory scoring.

The headspace layer now also carries a narrow pH-release correction for plant matrices: in pea/soy systems, acid-sensitive aldehydes and furans are allowed to release more strongly under acidic conditions, anchored so the existing Pratap-Singh pH ~6 baselines remain stable while the Pouvreau pea-isolate family is reproduced as an acidic-vs-less-acidic trend rather than as a forced absolute benchmark with incomplete metadata.

The matrix-accessibility layer now uses explicit native/desaturated endpoints for reactive lysine and cysteine instead of implicitly relaxing all the way to free-amino-acid behavior. In practice, pea-isolate cysteine remains strongly suppressed even after denaturation, matching the repo's literature note that free SH is near undetectable in commercial PPI.

For soy, the current envelope is still conservative, but it is no longer a free-floating placeholder: the repo now anchors it explicitly to the internal plant-protein literature synthesis that describes soy glycinin/β-conglycinin as compact globulins with buried reactive groups, soy as lysine-rich but sulfur-poor relative to meat-like systems, and extrusion / protein-polysaccharide conjugation as the mechanisms that partially reopen accessibility while still trapping volatiles. That is enough to justify the current relative contract `free > soy_iso > pea_iso`, but not enough yet to promote soy into an absolute quantitative accessibility benchmark.

For `matrix_only` today, the contract is narrow and explicit:

- The benchmark must declare a non-`free` `protein_type` with a registered matrix model.
- Validation is an intake/headspace execution check, not a general precursor-resolved chemistry claim.
- It can appear as `supported` in summary/index outputs while still remaining `strict_ready = False`.
- It is deliberately omitted from `targets-report` until there is a benchmark-facing target-ranking model for matrix systems.
- The new `protein_type` headspace fallback is a conservative bridge for unspecified matrix fractions, not a replacement for explicit pea/soy composition calibration.
- Denaturation-aware retention is now internally consistent across `src/matrix_correction.py`, `src/headspace.py`, `src/recommend.py`, and `src/sensory.py`, but its numeric endpoints are still literature-estimate placeholders rather than benchmark-fitted matrix physics.
- Denaturation-aware amino-acid accessibility is now explicit too, and canonical lysine/cysteine precursor keys are recognized inside the recommendation path so matrix corrections actually apply to those species.
- That accessibility contract is now frozen end-to-end by `tests/integration/test_matrix_accessibility_recommendation.py`: sulfur-rich meaty recommendations are penalized in the order `free > soy_iso > pea_iso`, and denatured pea recovers signal without collapsing back to free-amino-acid behavior.
- The unit layer now freezes the same relative hierarchy in `tests/unit/test_matrix_correction.py`, including the stricter expectation that soy and pea concentrates stay more buried than their corresponding isolates across native, midpoint, and denatured states.
- pH-dependent plant-matrix headspace validation is now covered as a relative trend check against the Pouvreau pea-isolate family; it is intentionally narrower than the absolute Pratap-Singh intake benchmarks.

This means the framework now has three strict-ready free-amino-acid benchmarks and two executable matrix-headspace intake benchmarks (`pea isolate` and `soy isolate`). The absolute concentration projection is materially better calibrated for the current free-precursor envelope, but the matrix headspace lane should still be treated as an executable intake model rather than a broadly calibrated release gate.

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
- `./scripts/docker_maillard.sh targets ...`: reproducible target snapshot for a specific benchmark and target family (`desirable`, `competing`, `toxic`; aliases `off_flavour` and `off-flavour`).
- `./scripts/docker_maillard.sh targets-report`: reproducible aggregate target artifact with headspace observability metadata.

## 4. Blind Spots That Still Matter

- **Matrix-only systems**: the first plant-isolate benchmark path is now executable, but broader matrix calibration beyond the pea-isolate intake case is still missing.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.

## 5. Recommended Verification Workflow

1. For macOS/OrbStack or Docker Desktop, bring up the validated Linux environment with `./scripts/docker_maillard.sh up` and create or refresh the env with `./scripts/docker_maillard.sh bootstrap`.
2. Run `./scripts/docker_maillard.sh summary` to inspect the current validated envelope.
3. Run `./scripts/docker_maillard.sh index` to inspect benchmark tiering and execution-path scope, especially for matrix-only benchmarks.
4. Run `./scripts/docker_maillard.sh pytest tests/scientific/test_free_aa_quantitative_regression.py -q` to guard the currently calibrated free-amino-acid ratios.
5. Run `./scripts/docker_maillard.sh pytest tests/` for the full Docker-validated suite.
6. Run `MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q` when you want an honest go/no-go signal for the strict free-AA gate.
7. Run `./scripts/docker_maillard.sh targets data/benchmarks/<benchmark>.json` when you need the current target snapshot for branch-level analysis.
8. Run `./scripts/docker_maillard.sh targets-report` when you need the aggregate target artifact regenerated before review or comparison.
9. Treat `matrix_only` benchmarks as executable intake checks unless they are explicitly promoted into the strict gate.

## 6. Expected Skips In The Docker Lane

- `tests/benchmarks/` is intentionally skip-heavy today. Those tests are Phase 3 placeholders and HPC-oriented literature checks, not part of the current release gate.
- Capability-gated QM tests should skip only when the backend is genuinely unavailable or unusable in the active Docker environment.
- A path-based skip for a binary that is actually present in `PATH` is a test bug, not a valid environment gate.

## 7. Named Docker Lanes

- `./scripts/docker_maillard.sh core`: unit and integration correctness gate.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index plus scientific regression lane.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend lane.
- `./scripts/docker_maillard.sh hofmann`: diagnostic trace for the calibrated Hofmann sulfur benchmark.
- `./scripts/docker_maillard.sh targets ...`: benchmark target snapshot lane for scientific inspection.
- `./scripts/docker_maillard.sh targets-report`: aggregate target artifact lane for scientific inspection and review.

These lanes keep infrastructure validation separate from unresolved scientific scope gaps such as matrix-only precursor systems.

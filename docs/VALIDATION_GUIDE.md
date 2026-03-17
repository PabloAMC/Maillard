# Scientific Validation Guide

If you are new to the repository, read [guides/SCIENTIFIC_RELIABILITY.md](guides/SCIENTIFIC_RELIABILITY.md) first. This document is the deeper contract and methodology reference.

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

The benchmark-comparison policy is now explicit too:

- `free_precursor` benchmarks compare the FAST observable path, i.e. the benchmark-facing `predicted_ppb` output from the recommender.
- Cantera is currently a diagnostic reference lane for mechanism export, temporal debugging, and chemistry investigations; it is not the authoritative comparator for strict benchmark pass/fail.
- `matrix_only` benchmarks compare the dedicated intake/headspace path and therefore remain outside the strict gate and outside FAST target snapshots.

The absolute projection contract is explicit as well:

- The current projection strategy is `precursor_limited_observable_v1` in `src/recommend.py`.
- The limiting precursor pool is converted from mM to M, multiplied by a conservative severity-dependent volatile yield fraction, allocated across selected terminal outputs, and then converted to ppb on an aqueous mass-equivalent basis before matrix and headspace penalties are applied.
- Those strategy constants and assumptions are emitted in `projection_context` and in each target's `projection_metadata`.

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
The index also records which engine is authoritative for each benchmark and whether Cantera is only diagnostic or not currently authoritative.

For a single graphical snapshot of reliability and current limitations, use:

```bash
./scripts/docker_maillard.sh validation-figures
```

This writes `results/validation/validation_overview.png`, `results/validation/validation_overview.md`, and `results/validation/validation_overview.json`.

For a reproducible thermodynamic-gating audit, use:

```bash
./scripts/docker_maillard.sh thermo-gating
```

This writes `results/validation/thermodynamic_gating_audit.md` and `results/validation/thermodynamic_gating_audit.json` and reports whether the gated variant materially improves benchmark error. Until that audit says otherwise, thermodynamic gating remains diagnostic-only for benchmark pass/fail.
The benchmark metadata and generated summary/index artefacts now also expose the current thermodynamic-gating policy so the `auto` benchmark path resolves through an explicit contract rather than an implicit default.

For the current matrix-specific validation surfaces, use:

```bash
./scripts/docker_maillard.sh matrix-deltas
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh matrix-branch-deltas main
./scripts/docker_maillard.sh coverage-gaps
```

These artifacts make matrix progress and matrix limitations explicit without pretending the matrix lane is already strict-ready.

For marker-separated validation lanes, use:

```bash
./scripts/docker_maillard.sh scientific-fast
./scripts/docker_maillard.sh kinetics-validation
```

The fast lane covers FAST regression and benchmark contract checks, while the kinetics lane covers the slower Cantera reference cases, including temperature-ramp validation.
Benchmarks with fewer than 3 matched compounds now report `pass-no-ranking` when full coverage and scale pass, which keeps the summary aligned with the contract that Pearson is not meaningful for 1-2 point systems.

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

The headspace layer now also carries the explicit observable calibration for the Pratap-Singh plant-matrix intake family: pea remains the reference yield baseline, soy-vs-pea release gaps for the tracked off-flavour markers are modeled in src/headspace.py, and the narrower pH-release correction still preserves the Pratap-Singh pH ~6 baselines while reproducing the Pouvreau acidic-vs-less-acidic pea trend.

The matrix-accessibility layer now uses explicit native/desaturated endpoints for reactive lysine and cysteine instead of implicitly relaxing all the way to free-amino-acid behavior. The code now also exposes explicit literature windows for the two isolate families it treats as quantitatively anchored: pea-isolate lysine is kept in the approximate 0.30-0.45 reactive envelope while pea cysteine remains near-zero to low-single-digit accessibility, and soy-isolate lysine/cysteine stay in a broader but still conservative 0.36-0.50 / 0.10-0.24 envelope derived from the repo's soy literature synthesis.

That calibration is now connected to a temperature/time/pH-aware denaturation heuristic in `src/matrix_correction.py`. In the formulation and optimizer paths, `denaturation_state` is no longer forced to `0.5` by default: if the user does not specify it explicitly, the engine infers an effective matrix state from the reaction conditions and exposes that `effective_denaturation_state` in the formulation result. Explicit overrides still win when the user or a benchmark needs to pin the state manually.

The user-facing explainability surface is also broader now:

- benchmark target artefacts expose proxy-vs-observable projection semantics directly (`Proxy ppb`, `Observable ppb`, `Obs/Proxy`, matrix factor, headspace factor, volatile class)
- formulation-level explainability can be generated with `./scripts/docker_maillard.sh explain-formulation <name>`
- domain-of-applicability warnings can be generated with `./scripts/docker_maillard.sh validated-envelope`

Finally, the previous bulk matrix fallback has been narrowed into a compound-class-aware overlay in the recommendation/headspace/sensory paths. Aldehydes and furans remain more strongly retained than sulfur volatiles and alcohols, which makes the default pea/soy fallback more interpretable without claiming benchmark-fitted explicit binding chemistry for every compound.

For soy, the current envelope remains conservative rather than benchmark-fitted, but it is no longer a free-floating placeholder: the repo now anchors it explicitly to the internal plant-protein literature synthesis that describes soy glycinin/β-conglycinin as compact globulins with buried reactive groups, soy as lysine-rich but sulfur-poor relative to meat-like systems, and extrusion / protein-polysaccharide conjugation as the mechanisms that partially reopen accessibility while still trapping volatiles. That is enough to justify both the relative contract `free > soy_iso > pea_iso` and the narrower calibrated window frozen in unit tests, while still keeping soy outside any strict absolute accessibility benchmark.

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
- `./scripts/docker_maillard.sh validated-envelope`: reproducible validated-envelope / domain-of-applicability artifact.
- `./scripts/docker_maillard.sh explain-formulation <name>`: reproducible formulation explainability artifact for a formulation in `data/formulation_grid.yml`.

## 4. Blind Spots That Still Matter

- **Matrix-only systems**: the first plant-isolate benchmark path is now executable and the default matrix-state handling is less manual, but broader matrix calibration beyond the current pea/soy intake family is still missing.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.
- **User-facing explainability**: the engine now computes a more physical matrix state internally, but that matrix metadata is not yet surfaced broadly in benchmark summary/report artifacts.
- **Default CLI integration**: explainability and validated-envelope artifacts exist, but the default user commands do not yet surface those warnings inline during every formulation run.

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

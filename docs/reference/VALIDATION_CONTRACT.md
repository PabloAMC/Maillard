# Scientific Validation Guide

If you are new to the repository, read [architecture.md](../architecture.md) first. This document is the deeper contract and methodology reference.

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

> **Provenance audit (2026-08-26).** An adversarial citation audit found that the
> sources behind three panel benchmarks are dead or resolve to unrelated papers:
> `cys_ribose_150C_Mottram1994`, `cys_glucose_150C_Farmer1999`, and
> `thiamine_cys_ribose_100C_Hofmann1996`. They are quarantined under
> `data/benchmarks/quarantined/` (excluded from the panel) pending human source
> verification — see the README there and `tasks/audit_remediation.md`. Two of
> them were previously listed below as strict-ready; that status is withdrawn.
> The same audit found the panel headline mixed literature-measured rows with
> internal synthetic comparators; coverage is now reported split by signal
> origin (see §3E).

As of the current benchmark summary:

- **Strict-ready count is now 0/16 (2026-08-27).** The panel has no strict-ready benchmark
  at all. The sequence, all on this branch and all recorded in `tasks/audit_remediation.md`:
  6/16 → 4/16 when the projection retune cost `spi_hvp_xylose_120C_PMC9905368` and
  `wheat_gluten_hvp_xylose_120C_PMC9905368` their status, then 4/16 → 0/16 when the Wave G1
  chemistry rebuild removed the fabricated MFT shortcut and the fabricated lipid radical
  chain. `cys_ribose_140C_Hofmann1998`, the last strict-ready free-amino-acid benchmark, now
  under-predicts MFT by 5.58x and FFT by 3.26x against a 1.45x contract; the safety-marker
  benchmarks fail on the separately documented unit collisions and the SME desaturation.
  Every one of these tolerances is UNCHANGED — nothing was widened to recover the count.
  The Hofmann-only refit that followed
  (`results/validation/sulfur_barrier_refit_hofmann.md`) established that no barrier value in
  any defensible range recovers it either: the residual is a volatile-budget ALLOCATION
  deficit, not a barrier deficit.
- `pea_isolate_40C_PratapSingh2021` is executable through a dedicated `matrix_only` intake path, with full coverage, but it remains outside the strict release gate. Note: its marker yields are back-fitted to this same measurement, so this benchmark validates the intake plumbing, not predictive accuracy.
- `soy_isolate_40C_PratapSingh2021` is executable through the same dedicated `matrix_only` intake path, with full coverage, with the same caveat.

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

This means the framework currently has NO strict-ready benchmark (2026-08-27; see §2) and two executable matrix-headspace intake benchmarks (`pea isolate` and `soy isolate`). The matrix headspace lane should be treated as an executable intake model rather than a broadly calibrated release gate.

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

### E. External Hold-out Methodology (added 2026-08-26)

The external hold-out is the only surface that tests the frozen model on data it
has never seen. Its rules, stated explicitly:

- **Exclusion.** Hold-out bundles live in `data/benchmarks/external_validation/`
  and carry `evidence_class = external_validation_only`. The benchmark loader's
  non-recursive glob makes the directory structurally invisible to every fitting
  and panel-evaluation path, and `calibrate_from_intake` raises on the evidence
  class as a belt-and-braces guard. No hold-out value appears in any fitted
  parameter.
- **Prior tier.** Hold-out envelopes use the `uncalibrated` matrix prior tier
  (matrix_headspace ln-sigma 2.86, ~±110x at 90% CI) because their process states
  have no calibrated registry entry. This is much wider than the in-panel tier,
  so **hold-out coverage and in-panel coverage are not comparable numbers**. The
  sigma is sized from data, not judgment (2026-08-26): leave-lane-out transfer
  error over the in-panel matrix anchors — hold-out untouched — gives RMS
  ln-sigma 2.86 with a 90% interval of [1.98, 5.48] (n=6, see
  `results/validation/matrix_sigma_residual_derivation.md`). The originally
  shipped 2.0 sat at the lower edge of that interval; the tier now uses the
  residual-derived point estimate.
- **Extrapolation vs re-scoring.** Bundles whose executable conditions are
  copied from an in-panel benchmark re-score that anchor's prediction against a
  different paper's measurement — a reproducibility comparison, not an
  extrapolation test. Only bundles at genuinely new process states (HME
  extrusion, roasting) test transfer, and the honest statement is that the model
  largely fails those.
- **Current numbers (2026-08-27 regeneration).** Nominal coverage **5/8**
  (62.5%): **3/3** on the re-scoring bundles (`bi_2020_raw_pea`,
  `liu_2023_ppi_offnote_baseline`) and **2/5** on the genuine extrapolations
  (`bi_2020_roasted_pea`, `li_2026_spi_wg_hme_control`). Median fold error
  **32.79x** (median |log10| 1.516 dex); worst miss 2474x on roasted pea. Read
  the fold error, not the coverage: the move from the previously reported 3/8 to
  5/8 came entirely from raising the uncalibrated sigma 2.0 → 2.86, which widens
  the interval without changing a single point prediction — the median fold
  error was 33.84x before and is 32.79x now.

### F. Reading The Coverage Split

`results/validation/prediction_uncertainty.md` now reports coverage split by
signal origin:

- **External literature** rows are validation evidence.
- **Internal synthetic** rows compare the model against its own frozen output
  (drift detection); they carry zero validation weight and exist so refactors
  cannot silently change predictions.
- **Not evaluable** rows have degenerate (near-zero-width) envelopes: the Monte
  Carlo perturbs nothing on their path, so pass/fail is meaningless for them.
- Coverage is only interpretable next to median CI width: a 90% interval
  spanning several orders of magnitude makes coverage cheap. Report both,
  always.

**Current in-panel numbers (2026-08-27 regeneration, n=200, seed 0).** 16 benchmarks,
41 matched rows, 13 of them in the Monte Carlo panel; aggregate 90% CI coverage
**28/41 (68.3%)**. Split:

| Signal origin | Inside 90% CI | Not evaluable | Median CI width |
| --- | ---: | ---: | ---: |
| External literature | **2/11 (18.2%)** | 4 | 3.00 dex |
| Internal synthetic | 18/18 (100%) | 8 | 3.19 dex |

Read the first row with the third column: the intervals GREW (2.0 → 3.0 dex) and coverage
still FELL (7/11 → 2/11). That is not a widening-of-error-bars story — the point predictions
moved further from the measurements than even a ~1000x interval allows, because the
2026-08-27 Wave G1 rebuild removed the fabricated chemistry the sulfur branch had been
standing on. Out-of-CI cells: 13/41. Strict-ready: 0/16 (see §2).

## 4. How To Read The Validation Figures

The repository now relies on a small set of figures that each have a distinct purpose.

- `validation_overview.png` is the first-pass trust surface. It is now a single-panel parity figure for the authoritative benchmark points, with formatted study references instead of raw benchmark ids, so it can be embedded directly in the README without multi-panel cropping.
- `validated_envelope.png` is the boundary figure. It shows how much of the benchmark set is strict-ready, how much is only directional, and which caveats still define the current edge of reliable use.
- `*_comparison.png` files are benchmark cards. They show parity, absolute yields, and the benchmark summary for one literature system at a time.

Read them in that order. The overview tells you whether the repo has a real quantitative proof surface. The envelope tells you where that proof stops. The benchmark card tells you why a single literature case passed, failed, or remained directional.

Correlation is informative but not sufficient on its own. The contract should always be interpreted through four signals together:

- coverage
- Pearson $R$ only when at least three compounds match
- MAE in ppb
- max ratio as the release-facing scale tolerance

This is why a benchmark may be scientifically acceptable as `pass-no-ranking` while still lacking a meaningful correlation coefficient.

## 5. Blind Spots That Still Matter

- **Citation provenance**: literature anchors were ingested by LLM-assisted deep-research campaigns and are not yet human-verified. An automated sweep (2026-08-26) found ~20% of registry DOIs unresolvable, plus a harder-to-detect class of live DOIs pointing at the wrong paper. Kinetic parameters and benchmark values should be treated as provisional pending manual source verification; `data/lit/computational_priors.json` (`dft_kinetic_priors`) is the provenance template the rest of the registries should be held to.
- **Matrix-only systems**: the first plant-isolate benchmark path is now executable and the default matrix-state handling is less manual, but broader matrix calibration beyond the current pea/soy intake family is still missing.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.
- **User-facing explainability**: the engine now computes a more physical matrix state internally, but that matrix metadata is not yet surfaced broadly in benchmark summary/report artifacts.
- **Default CLI integration**: explainability and validated-envelope artifacts exist, but the default user commands do not yet surface those warnings inline during every formulation run.

## 6. Recommended Verification Workflow

1. For macOS/OrbStack or Docker Desktop, bring up the validated Linux environment with `./scripts/docker_maillard.sh up` and create or refresh the env with `./scripts/docker_maillard.sh bootstrap`.
2. Run `./scripts/docker_maillard.sh summary` to inspect the current validated envelope.
3. Run `./scripts/docker_maillard.sh index` to inspect benchmark tiering and execution-path scope, especially for matrix-only benchmarks.
4. Run `./scripts/docker_maillard.sh pytest tests/scientific/test_free_aa_quantitative_regression.py -q` to guard the currently calibrated free-amino-acid ratios.
5. Run `./scripts/docker_maillard.sh pytest tests/` for the full Docker-validated suite.
6. Run `MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q` when you want an honest go/no-go signal for the strict free-AA gate.
7. Run `./scripts/docker_maillard.sh targets data/benchmarks/<benchmark>.json` when you need the current target snapshot for branch-level analysis.
8. Run `./scripts/docker_maillard.sh targets-report` when you need the aggregate target artifact regenerated before review or comparison.
9. Treat `matrix_only` benchmarks as executable intake checks unless they are explicitly promoted into the strict gate.

## 7. Expected Skips In The Docker Lane

- `tests/benchmarks/` is intentionally skip-heavy today. Those tests are Phase 3 placeholders and HPC-oriented literature checks, not part of the current release gate.
- Capability-gated QM tests should skip only when the backend is genuinely unavailable or unusable in the active Docker environment.
- A path-based skip for a binary that is actually present in `PATH` is a test bug, not a valid environment gate.

## 8. Named Docker Lanes

- `./scripts/docker_maillard.sh core`: unit and integration correctness gate.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index plus scientific regression lane.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend lane.
- `./scripts/docker_maillard.sh hofmann`: diagnostic trace for the calibrated Hofmann sulfur benchmark.
- `./scripts/docker_maillard.sh targets ...`: benchmark target snapshot lane for scientific inspection.
- `./scripts/docker_maillard.sh targets-report`: aggregate target artifact lane for scientific inspection and review.

These lanes keep infrastructure validation separate from unresolved scientific scope gaps such as matrix-only precursor systems.

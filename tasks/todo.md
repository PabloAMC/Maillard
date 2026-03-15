# Maillard Literature Replication Plan

## Active Matrix-Only Contract — 2026-03-15

- [x] Auditar y fijar el contrato explícito de benchmarks `matrix_only` en `src/benchmark_validation.py`.
- [x] Añadir tests de aceptación para `pea_isolate_40C_PratapSingh2021` que cubran soporte, artefactos y exclusión deliberada del snapshot de targets.
- [x] Aclarar en la documentación qué significa `matrix_only`, qué sí valida y qué queda fuera del strict gate.
- [x] Revalidar en Docker con `./scripts/docker_maillard.sh` dentro del entorno `maillard`.
- [x] Documentar el resultado y actualizar el review section.

## Active Headspace Calibration — 2026-03-15

- [x] Alinear `src/headspace.py` con los perfiles de retención pea/soy cuando no haya fracciones explícitas de matriz.
- [x] Añadir tests unitarios que congelen el fallback de headspace para `pea_iso` y `soy_iso` sin doble conteo.
- [x] Validar el subset de headspace, matrix_correction y scientific lane en Docker.
- [x] Documentar el alcance de esta calibración mínima y qué sigue pendiente.

## Short Priority — 2026-03-15

- [x] Corregir scoring por concentración (`src/recommend.py` / `src/inverse_design.py`) — validado en Docker (`tests/unit/test_safety_and_flux.py::test_concentration_aware_ranking` y `tests/integration/test_recommendation_engine.py::test_concentration_boltzmann_scoring`)
- [x] Añadir manejo robusto cuando CREST falla (capability probe / pytest marker) — validado en Docker (`tests/qm/test_solvation.py::test_explicit_solvation_run` y `tests/qm/test_explicit_solvation_integration.py::test_dft_refiner_explicit_solvation`)
- [x] Documentar y automatizar el flujo Docker reproducible (`scripts/docker_maillard.sh`, `README.md`, `Installation.md`, `docs/VALIDATION_GUIDE.md`)
- [x] Exponer una inspección reproducible de targets de benchmark sin `python -c` frágil (`./scripts/docker_maillard.sh targets ...`, con tipos `desirable`/`competing`/`toxic` y alias legibles)
- [x] Auditar los `skip` del suite y corregir guards obsoletos de backend (`tests/qm/test_dft_solvation_integration.py`, `tests/qm/test_sella_ts.py`)
- [x] Reprobar tests clave en Docker y cerrar Tier 0/benchmark lane — suite completa verde en Docker (`327 passed, 41 skipped, 4 xfailed`)
- [x] Documentar cambios y lecciones en `tasks/lessons.md`
- [x] Generar y validar resumen de benchmarks (`scripts/generate_benchmark_summary.py`) — validado en Docker

## Next Scientific Priority — 2026-03-15

- [x] Cerrar el `scale-gap` de `cys_ribose_140C_Hofmann1998` sin reabrir Farmer ni Mottram. Docker-validado con `./scripts/docker_maillard.sh summary` (`max ratio 1.442`, `strict-ready: yes`).
- [x] Diseñar una lane Docker explícita para `core`, `scientific` y `qm-heavy`, con skips justificados por capacidad real.
- [x] Convertir la inspección puntual de Hofmann/Farmer/Mottram en un comando reproducible de snapshot de targets dentro del wrapper Docker.
- [x] Exportar snapshots agregados de targets de benchmark a `results/validation/benchmark_targets.{md,json}` y regenerarlos desde la lane `scientific`.
- [x] Definir el primer camino ejecutable para benchmarks matrix-only (`pea_isolate_40C_PratapSingh2021`) sin falsear soporte científico.

## Validation Contract

- [ ] Define the project-wide meaning of literature replication.
- [ ] Lock benchmark tiers and acceptance thresholds in one place.
- [ ] Separate directional validity, quantitative replication, and formulation utility in docs and tests.

### Acceptance thresholds

- [ ] Free amino acid PRIMARY systems: Pearson R >= 0.85 when at least 3 matched compounds are available.
- [ ] Free amino acid PRIMARY systems: major volatile ratio <= 1.5x once strict mode is enabled.
- [ ] Plant matrix PRIMARY systems: major volatile ratio <= 2.0x once headspace and matrix corrections are calibrated.
- [x] Benchmarks with fewer than 3 matched compounds: do not report Pearson R; rely on coverage, ratio, and MAE only.

## Priority Reset — 2026-03-14

### Immediate Utility Priorities

1. Surface a single benchmark summary artifact so users can see supported benchmarks, ranking quality, and remaining scale gaps without reading test logs.
2. Turn the free-amino-acid benchmark gate from soft guidance into an explicit strict mode once the summary is the default diagnostic surface.
3. Defer matrix and broader chemistry work until the tool communicates its validated envelope clearly.

- [x] Confirm the real bottleneck with Docker-only benchmark runs in the `maillard` conda env.
- [x] Confirm that the main free-amino-acid systems are no longer blocked by species coverage.
- [x] Treat absolute concentration projection from `src/recommend.py` through `src/inverse_design.py` as the current P0 blocker.
- [x] Freeze additional sulfur-chemistry expansion unless a PRIMARY benchmark loses coverage again.
- [ ] Defer plant-matrix calibration work until the free-amino-acid projection is quantitatively credible.
- [x] Write a benchmark summary artifact that separates ranking success from absolute-scale failure.
- [x] Make the free-amino-acid strict gate consume the same summary criteria instead of a loose per-compound tolerance.
- [x] Confirm via Docker sweeps that no single local sulfur-barrier tweak closes Hofmann, Mottram, and Farmer simultaneously.

## Test Stability Gate — 2026-03-14

Rationale: yes, but not as an undifferentiated "make all pytest green" task. The useful next move is to restore a trustworthy test contract in layers, so Phase B-C work is not driven by a noisy or internally inconsistent suite.

### Tier 0: Re-establish the core deterministic contract

- [ ] Fix the barrier constant contract drift so `tests/unit/test_arrhenius_params.py` and the calibrated FAST barrier table agree on the intended canonical value for `1,2-enolisation`.
- [ ] Fix the headspace protein-sequestration invariant so `tests/unit/test_headspace.py` reflects the current partitioning model and the implementation matches the documented physics.
- [x] Fix concentration-aware ranking so `tests/unit/test_safety_and_flux.py` and `tests/integration/test_recommendation_engine.py` once again respond monotonically to precursor loading.
- [x] Remove or gate stray debug output in the recommendation path while stabilizing the concentration-sensitive tests.

### Tier 1: Fix mechanism serialization and simulation validity

- [ ] Make exported Cantera species and reaction equations parser-safe for names such as `bis(2-methyl-3-furyl) disulfide` so the current parser failures disappear from `tests/integration/test_fft_bottleneck.py` and `tests/integration/test_regression.py`.
- [ ] Re-run the Cantera simulation suite after the serialization fix and determine whether the remaining failures are expectation drift or real kinetics regressions.
- [ ] Reconcile `tests/integration/test_cantera_sim.py` with the current reversible chemistry and barrier model instead of letting outdated magnitude assertions mask real regressions.
- [ ] Restore the ribose plus cysteine plus leucine literature gate so Strecker output (`3-Methylbutanal`) is visible again in `tests/integration/test_barrier_calibration.py`.

### Tier 2: Separate environment failures from code failures

- [x] Treat the CREST/QCG segmentation faults in `tests/qm/test_solvation.py` and `tests/qm/test_explicit_solvation_integration.py` as an environment-stability issue first; verify whether they reproduce inside the required Docker `maillard` environment.
- [x] If CREST is unstable only on the host, add an explicit capability check or pytest marker so unsupported host setups skip cleanly instead of failing the suite.
- [x] Document which test subsets are required to pass in Docker before Phase B-E chemistry work continues.

### Tier 3: Lock the suite into useful execution lanes

- [x] Define a small "core correctness" lane that must pass before any new benchmark-calibration work.
- [x] Define a "scientific validation" lane for FAST plus benchmark tests in Docker.
- [x] Define a "heavy kinetics/QM" lane for Cantera and external-tool tests, with explicit environment prerequisites.
- [x] Update docs and local commands so contributors stop interpreting host-only external-tool crashes as chemistry regressions.

### Acceptance criteria for this gate

- [ ] All Tier 0 failures are green locally and in Docker.
- [ ] The Cantera parser failures are gone before any further benchmark recalibration.
- [ ] The remaining Cantera assertions are either fixed or intentionally re-baselined with scientific justification.
- [ ] External-tool tests fail only for real code defects, not for missing or unstable binaries on unsupported environments.

## Phase A: Benchmark Infrastructure [Low]

- [x] Create initial benchmark JSON files.
- [x] Add shared benchmark execution core in `src/benchmark_validation.py`.
- [x] Refactor `tests/scientific/test_benchmarks.py` to use the shared core.
- [x] Refactor `scripts/compare_sim_to_lit.py` to use the shared core.
- [x] Make benchmark execution work inside the Docker conda env.
- [x] Prevent misleading Pearson R reporting for 2-point comparisons.
- [x] Tighten fuzzy matching so trivial substrings do not count as species matches.
- [ ] Add benchmark metadata fields for tier and benchmark family directly in each JSON file.
- [x] Add a benchmark index in code so tests and reports can filter by PRIMARY, SECONDARY, matrix, and safety.

## Phase B: Free Amino Acid Replication Gate [Medium]

- [x] Make `cys_ribose_140C_Hofmann1998` quantitatively credible for MFT and FFT.
- [x] Make `cys_ribose_150C_Mottram1994` reproduce sulfur ordering without absurd absolute overprediction.
- [x] Make `cys_glucose_150C_Farmer1999` discriminate ribose vs glucose cleanly.
- [ ] Audit the scaling path from `src/recommend.py` to `src/inverse_design.py` so `predicted_ppb` is physically interpretable, not just a heuristic weight.
- [ ] Add species alias tables for benchmark matching where literature and internal names diverge.
- [ ] Enable strict ratio assertions behind the benchmark flag once the first three PRIMARY systems are stable.

## Phase C: FAST Physics and Quantitative Scaling [High]

- [ ] Replace arbitrary ppb scaling in `src/recommend.py` with a documented concentration projection strategy.
- [ ] Review the relationship between barrier, rate constant, and output concentration in `src/conditions.py`, `src/recommend.py`, and `scripts/run_cantera_kinetics.py`.
- [ ] Add temperature and time monotonicity tests for benchmarked compounds.
- [ ] Add ratio-sensitivity tests for precursor loading, especially cysteine and ribose.
- [ ] Decide whether benchmark comparison should use FAST output, Cantera output, or both, and document the intended role of each.
- [ ] Implement thermodynamic gating where it materially changes benchmark error rather than as a roadmap placeholder.
- [x] Add explicit observability/headspace groundwork in `src/recommend.py` so the projection layer can distinguish physically low-headspace species before the full benchmark-facing projection redesign is enabled.

## Phase D: Plant Matrix Replication [Very High]

- [x] Promote `pea_isolate_40C_PratapSingh2021` from xfail to executable benchmark.
- [x] Add matrix-aware precursor handling or a dedicated matrix benchmark pathway that does not rely on free-precursor resolution.
- [ ] Calibrate `src/headspace.py` against pea and soy headspace literature.
- [ ] Calibrate `src/matrix_correction.py` against reactive lysine and cysteine accessibility literature.
- [ ] Add pH-dependent headspace validation using the Pouvreau benchmark family.
- [ ] Keep plant-matrix benchmarks outside the strict gate until coverage and matrix physics are both credible.

## Phase E: Safety and Temporal Validation [Medium]

- [ ] Validate `src/safety.py` against acrylamide formation and elimination literature, not just monotonic formation.
- [ ] Add an explicit non-monotonic acrylamide benchmark.
- [ ] Validate temperature-ramp predictions against the temporal FAST suite and Cantera reference cases.
- [ ] Separate fast scientific regression tests from slower kinetics validation tests with pytest markers.

## Phase F: Regression Gate and Reporting [Low]

- [x] Add a benchmark summary report that writes coverage, MAE, ratio failures, and supported vs unsupported benchmarks.
- [x] Expose strict benchmark mode in documented local commands.
- [x] Update `docs/VALIDATION_GUIDE.md` so it reflects the actual gate instead of aspirational release checks.
- [x] Update `README.md` and `docs/architecture.md` to distinguish current validated capability from long-term SOTA goals.

## Research Roadmap After Replication Gate

- [ ] Peptide-bound reactivity in SmirksEngine.
- [ ] Oligosaccharide and maltodextrin precursor handling.
- [ ] Broad-spectrum phytochemical scavenging.
- [ ] Flavor-texture coupling through DHA and rheology.
- [ ] High-shear extrusion dynamics.
- [ ] Large-scale DFT and delta-ML scaling.

## Current Status Notes

- [x] Shared benchmark infrastructure is now in place.
- [x] The Docker conda env named `maillard` is usable for validation runs.
- [x] The benchmark script runs in Docker and produces reports and plots.
- [x] Spurious lipid oxidation contamination of free-amino-acid benchmarks was removed.
- [ ] Free-amino-acid benchmarks still show severe absolute scaling error, but the projection no longer wastes volatile budget on stoichiometric coproducts.
- [ ] The current FAST output is still a proxy signal, not a validated concentration model.
- [x] The glucose benchmark chemistry coverage gap is closed; Farmer now has full species coverage.
- [x] The Mottram disulfide now receives non-zero projected ppb after aligning the curated target SMILES with the generated species identity.
- [x] Pentose and hexose sulfurization are now separated in the FAST barrier model, and pyrazine condensation is less overexpressed in acidic cysteine systems.
- [x] A previous Pearson R = 1.0000 report was misleading because it came from a 2-point benchmark plus permissive name matching; the validation layer now guards against that.
- [x] Current Docker measurement confirms `cys_ribose_140C_Hofmann1998` has 100% coverage, but still misses absolute scale badly: MFT ratio ~653x, FFT ratio ~253x.
- [x] Current Docker measurement confirms `cys_ribose_150C_Mottram1994` has 100% coverage and useful ranking (`Pearson R ~= 0.87`), but absolute ratios remain far off: MFT ~185x, disulfide ~327x, furfural ~853x.
- [x] Current Docker measurement confirms `cys_glucose_150C_Farmer1999` has 100% coverage and strong ranking (`Pearson R ~= 0.99`), but absolute ratios remain far off: MFT ~80x, furfural ~313x, pyrazine ~69x.
- [x] Current Docker measurement confirms `pea_isolate_40C_PratapSingh2021` now runs through a dedicated `matrix_only` intake path with 100% coverage while remaining outside the strict release gate.
- [x] `predicted_ppb` now excludes exogenous precursors and other non-output species; the interface is once again semantically aligned with projected volatile outputs.
- [x] Projection now allocates the volatile budget across curated, non-exogenous, budget-relevant outputs instead of diffusing it across the full tracked network.
- [x] The dominant program risk is no longer global quantitative compression of `predicted_ppb`; the remaining benchmark gaps are now compound-family-specific.
- [x] The benchmark path and the diagnostic path now both propagate `time_minutes` into FAST prediction, so branch calibration is temporally consistent.
- [x] Residual branch calibration has been localized and corrected at the family level instead of by ad hoc benchmark-specific patches.
- [x] Current benchmark-summary status is now explicit: Farmer, Hofmann, and Mottram are `strict-ready`; the pea-isolate case is executable as `matrix_only` but not `strict-ready`.
- [x] Docker sweeps confirmed that the remaining free-AA blocker is not a missing single sulfur barrier entry but the translation from FAST activity into observed concentration/headspace.
- [x] `src/recommend.py` now carries Henry-constant lookup and observability classification helpers, but the benchmark-facing budget still intentionally stays on the stable pre-headspace allocation until the full projection redesign is ready.
- [x] The benchmark-facing budget now applies a conservative post-projection penalty to explicitly low-headspace targets (for example HMF) without redistributing that budget into the validated free-amino-acid benchmark compounds.
- [ ] A full host `python -m pytest tests/` run is currently not a trustworthy release gate: 13 failures span core logic drift, Cantera serialization breakage, expectation drift in kinetics tests, and CREST host instability.
- [ ] The highest-value stabilization work is to fix Tier 0 logic and Tier 1 serialization before returning to deeper Phase B-E benchmark calibration.

## Review Section

- Verified in Docker with `conda activate maillard` on Python 3.12.
- Headspace calibration result: `src/headspace.py` now applies a conservative `protein_type` fallback for `pea_iso` and `soy_iso` only when explicit matrix fractions are absent, keeping explicit fat/protein-fraction physics unchanged and validated by Docker (`27 passed` on the focused unit subset, `66 passed, 3 xfailed` on `./scripts/docker_maillard.sh scientific`).
- Matrix-only contract review result: `pea_isolate_40C_PratapSingh2021` remains executable and supported in summary/index, remains outside the strict gate, and is now explicitly excluded from FAST target snapshots and `targets-report` by contract.
- Docker reproducibility is now standardized around `./scripts/docker_maillard.sh`; `status`, `summary`, `index`, targeted `pytest`, and full `pytest tests/` all ran successfully in the validated container.
- Target-level benchmark introspection is now standardized too: `./scripts/docker_maillard.sh targets data/benchmarks/<benchmark>.json [desirable|competing|toxic]` replaces brittle ad hoc inline Python for scientific inspection, while still accepting off-flavour aliases.
- The skip audit is now explicit: placeholder-heavy `tests/benchmarks/` remain intentionally out of the release gate, while QM skips must be capability-based rather than path-based.
- The Hofmann sulfur benchmark is now back inside the strict-ready envelope through a projection-layer calibration in `src/recommend.py`, validated via `./scripts/docker_maillard.sh summary` and `./scripts/docker_maillard.sh scientific`.
- The next scientific expansion target is no longer free-amino-acid sulfur calibration but extending the executable envelope beyond `free_precursor`, starting with matrix-aware benchmark handling.
- Current scientific harness status: `tests/unit/test_budget_projection.py` plus the key scientific suite now pass in Docker (`8 passed, 1 xfailed`).
- Re-prioritization outcome from the latest Docker benchmark readout:
  - Coverage for the free-amino-acid PRIMARY systems is no longer the gating issue.
  - Farmer and Mottram now have enough ordering fidelity that more chemistry expansion is unlikely to be the highest-value next move.
  - The largest blocker is that `predicted_ppb` is still a proxy projection with severe magnitude compression, so benchmark work should center on projection semantics and benchmark-facing scaling before more matrix or SOTA work.
- Projection refactor result verified in Docker:
  - Focused regression suite now passes at `14 passed, 1 xfailed`.
  - `cys_ribose_140C_Hofmann1998`: MFT improved from ~653x error to ~2.5x; FFT improved from ~253x to ~3.2x.
  - `cys_ribose_150C_Mottram1994`: furfural improved from ~853x to ~1.45x and MFT from ~185x to ~1.74x, while disulfide remains strongly underpredicted and is now clearly a branch-specific chemistry/selectivity problem.
  - `cys_glucose_150C_Farmer1999`: furfural improved from ~313x to ~1.25x, while glucose-derived MFT and pyrazine remain underpredicted and now look like branch-specific pathway calibration gaps rather than a global projection collapse.
- Temporal consistency and local branch calibration result verified in Docker:
  - Focused regression suite now passes at `18 passed, 1 xfailed`.
  - `cys_ribose_140C_Hofmann1998`: MFT is now ~2.52x and FFT ~1.58x.
  - `cys_ribose_150C_Mottram1994`: MFT is now ~1.12x, furfural ~1.72x, and bis(2-methyl-3-furyl) disulfide ~1.07x with Pearson `~0.9994`.
  - `cys_glucose_150C_Farmer1999`: MFT is now ~1.76x, furfural ~1.23x, and pyrazine ~1.23x with Pearson `~0.9999`.
  - The remaining unsupported PRIMARY-adjacent benchmark is still the matrix-only pea-isolate case, which is blocked by execution-path scope rather than free-AA selectivity.
  - Local sulfur-family refinement in `src/barrier_constants.py` now improves the last quantitative free-AA outliers without reopening the global projection problem: Hofmann reaches MFT ~1.99x and FFT ~1.52x, Mottram holds MFT ~1.30x / disulfide ~1.03x / furfural ~1.78x, and Farmer reaches MFT ~1.45x / furfural ~1.22x / pyrazine ~1.22x.
  - Added `tests/scientific/test_free_aa_quantitative_regression.py` so these benchmark-family ratios are enforced directly instead of relying on ad hoc manual inspection.
  - Added `scripts/generate_benchmark_summary.py` plus summary helpers in `src/benchmark_validation.py`; Docker now writes `results/validation/benchmark_summary.md` and `results/validation/benchmark_summary.json`, making supported benchmarks, ranking quality, and remaining scale gaps visible in one place.
  - The opt-in `MAILLARD_STRICT_BENCHMARKS=1` path in `tests/scientific/test_benchmarks.py` now uses the same centralized summary thresholds, so strict failures report benchmark-level blocking issues instead of a disconnected per-compound tolerance.
  - Current stable Docker regression state after the latest review is `14 passed, 1 xfailed` for `tests/unit/test_budget_projection.py`, `tests/scientific/test_free_aa_quantitative_regression.py`, `tests/scientific/test_benchmark_summary.py`, and `tests/scientific/test_benchmarks.py`.
  - Short Docker sweeps on local sulfur-family barriers improved individual benchmarks but could not make Hofmann, Mottram, and Farmer all pass together; the next P0 is the concentration/headspace projection layer rather than another local barrier tweak.
  - A conservative B-C groundwork step is now in place: the recommender can classify low-headspace species from Henry constants and preserve multi-role target metadata without disturbing the current validated benchmark budget.
  - The next conservative step is now in place too: low-headspace targets are lightly suppressed at the projection output stage, so the benchmark-target artifact reflects headspace observability rather than raw liquid-phase proxy mass for near-nonvolatile species.
  - Latest conservative projection step: the low-headspace suppression now uses the existing HeadspaceModel temperature dependence inside Docker validation, but remains capped so free-amino-acid benchmark calibration stays stable (`15 passed` for budget/targets tests and `7 passed` for the benchmark scientific subset).
  - Projection semantics are now explicit in `src/recommend.py`: FAST first produces a proxy output budget, then a separate observable-output projection applies fallback matrix retention and relative headspace suppression; Docker suite remains green at `337 passed, 41 skipped, 4 xfailed`.
  - The first matrix-only benchmark lane is now executable too: `pea_isolate_40C_PratapSingh2021` runs through a dedicated oxidation/headspace intake model with 100% coverage, `Pearson 1.000`, `max ratio 1.002`, and remains outside the strict gate by design.
- Root causes identified so far:
  - Free-amino-acid benchmarks were previously contaminated by lipid oxidation products injected without lipid inputs.
  - Benchmark matching was too permissive and could overstate model performance.
  - FAST output scaling remains architecturally inconsistent with experimental headspace concentrations.
  - The FAST projection was also diluting endpoint concentrations by allocating volatile budget to bookkeeping coproducts such as water, CO2, H2, sulfur, and formaldehyde.
  - After filtering the budget to curated targets and budget-relevant volatile endpoints, Farmer and Mottram roughly doubled their predicted ppb scale, so the remaining gap is calibration/selectivity rather than coproduct dilution.
  - A further projection bug was that `predicted_ppb` still mixed real volatile outputs with exogenous precursors and arbitrary tracked species; constraining it to curated, non-exogenous endpoints removed that interface-level distortion and closed the large global compression error for the main benchmark compounds.
  - A second structural bug was that the main benchmark path did not propagate formulation `time_minutes` into FAST prediction, so branch diagnostics and benchmark evaluation were calibrating different temporal models.
  - Once temporal semantics were unified, the residual free-AA benchmark errors became cleanly localisable to three families: pentose thiol addition, hexose thiol addition, and furyl disulfide / pyrazine-forming condensation.
  - The Mottram disulfide zero was traced to a curated-target identity mismatch, not to missing chemistry; fixing the target SMILES restored a non-zero disulfide projection.
  - A single generic sulfurization barrier was flattening ribose-system selectivity. Splitting pentose vs hexose thiol addition and raising aminoketone condensation improved Farmer and Mottram simultaneously while keeping Hofmann closer to the expected MFT/FFT balance.
  - The current failing full-suite run is not one bug. It decomposes into four clusters: Tier 0 logic-contract drift (`arrhenius`, `headspace`, concentration-aware ranking), Tier 1 Cantera serialization and simulation semantics, Tier 1 pathway-selectivity regression in the leucine Strecker branch, and Tier 2 host-side CREST instability.
  - Because those clusters cut across the same projection and mechanism-export layers that Phase B-C depend on, stabilizing the targeted test gate first is the more efficient path than continuing benchmark calibration on top of a partially broken suite.
  - Latest Docker stabilization result:
    - The priority subset `tests/scientific/test_free_aa_quantitative_regression.py`, `tests/unit/test_safety_and_flux.py::test_concentration_aware_ranking`, `tests/integration/test_recommendation_engine.py::test_concentration_boltzmann_scoring`, `tests/qm/test_solvation.py::test_explicit_solvation_run`, and `tests/qm/test_explicit_solvation_integration.py::test_dft_refiner_explicit_solvation` now passes at `5 passed`.
    - Host-only CREST crashes are now separated from real code defects by a runtime QCG capability probe rather than a bare binary-presence check.
    - `DFTRefiner` no longer treats optional TS-optimizer import failures such as JAX/AVX runtime incompatibilities as hard blockers for explicit-solvation smoke tests.

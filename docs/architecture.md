# Maillard Reaction Computational Framework

## Purpose

Maillard exists to help alternative-protein scientists narrow the wet-lab search space before committing to GC-MS or process work.

The repository is designed to answer questions like these:

- which precursor combinations are directionally most promising for meaty sulfur chemistry
- which changes improve desirable aroma without worsening off-flavour or safety risk
- which predictions are benchmark-backed and which are only heuristic or matrix-directional

The framework is useful when it improves prioritization. It is not useful if it creates confident-looking but weakly anchored quantitative claims.

## Current Validated Scope

The most important distinction in this repository is between:

- the current validated operating surface
- the broader architecture the codebase is capable of growing into

Those are not the same thing.

### What Is Shareable Today

- free-precursor screening inside the strict benchmark envelope
- scientist-facing run reports with confidence and provenance
- side-by-side comparison artifacts for named formulations
- campaign-level shareable packages with per-run bundles and review context
- matrix intake and directional ranking for pea and soy with explicit caveats

### What Is Not Yet Closed Scientifically

- broad quantitative matrix prediction across plant-protein systems
- strict-gate meaty-positive promotion for pea or soy matrices
- benchmark-closed transfer across real process states such as extrusion-heavy histories
- internal experiment-ingestion loops that automatically convert wet-lab data into new benchmark and calibration surfaces

## Current Trust Surface

### High trust

- free-precursor comparative screening
- benchmark-aware ranking inside the strict-ready literature envelope
- safety-aware comparisons that stay close to the validated benchmark neighborhood

### Moderate trust

- pea and soy matrix comparisons when used as directional prioritization
- matrix explainability and off-flavour triage
- deciding what to test next in the wet lab

### Low trust

- new protein families without nearby benchmarks
- intact-protein or peptide-bound systems
- extrusion-heavy process claims without dedicated evidence
- absolute matrix concentration claims beyond the validated envelope

For the review surface that should accompany external sharing, see [guides/SCIENTIFIC_RELIABILITY.md](guides/SCIENTIFIC_RELIABILITY.md), [reference/SCIENTIFIC_REFERENCE.md](reference/SCIENTIFIC_REFERENCE.md), and [guides/SHARING_RESULTS.md](guides/SHARING_RESULTS.md).

## Architecture Layers

The repository uses a layered architecture so benchmark-backed screening, explainability, and future higher-fidelity physics can coexist without pretending they are all equally validated.

### 1. Reaction Enumeration

Main module: `src/smirks_engine.py`

Role:

- generate deterministic, atom-balanced reaction candidates from precursor sets
- keep the chemistry surface explicit and inspectable

Why it matters:

- if the reaction graph is not chemically coherent, no later scoring or DFT work is trustworthy

### 2. FAST Screening And Ranking

Main modules:

- `src/recommend.py`
- `src/inverse_design.py`
- `src/barrier_constants.py`

Role:

- produce fast laptop-feasible rankings
- propagate literature-calibrated barrier assumptions into a formulation score
- support forward prediction, inverse design, and scientist-facing summaries

Current state:

- this is the main validated daily-use surface inside the free-precursor benchmark envelope

### 3. Matrix And Headspace Translation

Main modules:

- `src/matrix_correction.py`
- `src/headspace.py`
- `src/matrix_calibration_registry.py`

Role:

- translate beaker chemistry into plant-matrix observability
- model accessibility, denaturation, retention, and release effects

Current state:

- useful for directional pea/soy work
- not yet benchmark-closed for release-grade quantitative claims

### 4. Safety And Decision Support

Main modules:

- `src/safety.py`
- `src/sensory.py`
- `src/usability_reports.py`
- `src/reporting.py`

Role:

- map chemistry into scientist-facing decisions
- expose confidence, calibration diagnostics, and reportable artifacts
- keep the trust boundary visible in every result surface

Current state:

- this layer is mature enough for sharing and review
- its quality still depends on the scientific support level of the underlying prediction mode

### 5. Validation And Benchmark Surfaces

Main modules:

- `src/benchmark_validation.py`
- `src/validation_contract.py`

Main artifacts:

- `results/validation/benchmark_summary.md`
- `results/validation/validated_envelope.md`
- `results/validation/validation_overview.png`
- `results/validation/matrix_benchmark_assertions.md`

Role:

- define what counts as benchmark-backed evidence
- keep matrix evidence, internal reproducibility harnesses, and strict-ready benchmarks separate

### 6. Future Higher-Fidelity Physics

Main modules:

- `src/mlp_barrier.py`
- `src/diffusion_ts.py`
- `src/dft_refiner.py`
- `src/skala_refiner.py`

Role:

- improve barrier quality and structural fidelity for selected high-value steps
- support future expansion of the chemistry surface and higher-fidelity rate estimation

Current state:

- architecturally important
- not the main reason the tool is or is not scientist-shareable today

## What Makes The Tool Useful In Practice

The repository is most useful when it does four things well:

- ranks candidate formulations inside a known trust boundary
- exposes why a result looks the way it does
- shows confidence and calibration posture explicitly
- produces artifacts that can be reviewed without re-running the whole CLI by hand

That is why the reporting and validation surfaces matter as much as the chemistry internals.

## What Still Blocks State-Of-The-Art Status

Scientifically, the main blockers are:

- no primary quantitative pea or soy meaty-positive matrix benchmark with desirable and adverse targets in the same run
- no benchmark-closed time-series matrix data for the target PPI and SPI systems
- incomplete process-state calibration for real matrix release behavior

Architecturally, the main blockers are:

- no full experiment-ingestion path for internal wet-lab data
- no external-team API layer beyond scripts and report artifacts
- no automatic loop from new experiments to benchmark promotion review

## Review Sequence For External Sharing

If the goal is to share repository state with scientists, reviewers, or collaborators, the shortest honest review sequence is:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh campaign data/campaigns/shareable_meaty_screen.yml
```

That sequence gives:

- benchmark status
- the validated boundary
- the current reliability visuals
- matrix readiness posture
- a reproducible campaign package with provenance

## Immediate Scientific Contract Gap

The highest-value missing experiment is now explicit and documented in [protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md](protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md).

That protocol is the bridge between the current directional matrix surface and the next benchmark-backed scientific upgrade.

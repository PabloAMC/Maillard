# Maillard — Architecture & Design

## Purpose

Maillard helps alternative-protein scientists narrow the wet-lab search space before committing
to GC-MS runs or process trials. It answers:

- Which precursor combinations are most promising for meaty sulfur chemistry?
- Which changes improve desirable aroma without worsening off-flavour or safety risk?
- Which predictions are benchmark-backed and which are only directional?

The framework is useful when it improves prioritisation. It is not useful if it produces
confident-looking but weakly anchored quantitative claims.

---

## Design Principles

**Observable-first.** Every model tier must eventually face an observable. Claims that cannot
be reduced to a measurable headspace concentration, benchmark ratio, or kinetic observable are
kept explicitly labelled as exploratory.

**Tiered confidence, not a single model.** The system treats free-precursor FAST kinetics,
matrix physics, and selective QM as separate, independently validated layers — not a single
monolithic pipeline. Each layer has its own trust tier; results are never silently promoted
across tiers.

**Speed vs. rigour as a dial, not a binary.** Laptop-speed FAST-mode screening (< 1 s) and
multi-hour DFT refinement coexist in the same codebase. Scientists pick the appropriate level
for each decision; the tool makes the tradeoffs visible.

**Three operating regimes**

| Regime | Timescale | Primary use | Trust level |
|--------|-----------|-------------|-------------|
| FAST screening | < 1 s | Formulation ranking, what-if questions | Calibrated heuristics from literature |
| Matrix physics | seconds | Headspace / retention / pea-soy translation | Directional; benchmark-partial |
| Selective QM | hours–days | Barrier refinement on specific rate-limiting steps | xtb_derived → DFT → literature surrogate |

---

## Trust Surface

### High trust — use freely
- Free-precursor comparative screening inside the benchmark envelope
- Benchmark-aware ranking (literature barrier constants; see `src/barrier_constants.py`)
- Safety-aware comparisons close to the validated benchmark neighbourhood

### Moderate trust — useful for prioritisation, verify before deciding
- Pea and soy matrix directional comparisons
- Matrix explainability and off-flavour triage
- Deciding what to test next in the wet lab

### Low trust — exploratory only
- New protein families without nearby benchmarks
- Intact-protein or peptide-bound systems
- Extrusion-heavy process claims without dedicated evidence
- Absolute matrix concentration claims beyond the validated envelope

> **Confidence tiers in code**: priors carry `bounded_calibration` / `transferred_literature` /
> `surrogate_family` / `xtb_derived` labels. These propagate through the pipeline and surface
> in every report. Never silently upgrade a tier.

---

## Architecture Layers

```mermaid
graph TD
    subgraph "1. Input & Accessibility"
        A["Precursors (Sugars, AAs, Lipids)"] --> C["Matrix Correction"]
        B["Process Conditions (pH, T, t, Matrix)"] --> C
    end

    subgraph "2. Reactive Core"
        C --> | "Accessible Molarity" | D["SMIRKS Rule Engine"]
        D --> | "Reaction Network" | E["Thermodynamic Gating (Joback)"]
        E --> | "Feasible Paths" | F["Cantera ODE Solver"]
        G["Literature Kinetics"] --> F
    end

    subgraph "3. Observability Projection"
        F --> | "Aqueous Moles" | I["Projection Module"]
        I --> | "Volatilization (Henry's Law)" | J["Headspace Calibration"]
        I --> | "Surface Adsorption" | K["Matrix Retention"]
    end

    subgraph "4. Decision Layer"
        J & K --> L["Validation Surface"]
        L --> M["Bayesian Optimizer"]
        M --> | "Formulation Candidate" | A
    end
```

### Layer 1 — Reaction Enumeration

**Module:** `src/smirks_engine.py`

Generates deterministic, atom-balanced reaction candidates from precursor sets. The chemistry
surface is kept explicit and inspectable — if the reaction graph is not chemically coherent,
no later scoring or QM work is trustworthy.

### Layer 2 — FAST Screening & Ranking

**Modules:** `src/recommend.py`, `src/pipeline.py`, `src/barrier_constants.py`

Produces laptop-feasible rankings in < 1 s using literature-calibrated barrier constants
(Schiff base 15 kcal/mol → enolisation 28 kcal/mol; anchored to Hofmann/Martins/Nursten
references). This is the main validated daily-use surface.

### Layer 3 — Matrix & Headspace Translation

**Modules:** `src/matrix_correction.py`, `src/headspace.py`, `src/matrix_calibration_registry.py`

Translates beaker chemistry into plant-matrix observability: accessibility, denaturation,
retention, and headspace release. Useful for directional pea/soy work; not yet
benchmark-closed for release-grade quantitative claims.

**Calibration single-application rule:** `HeadspaceModel.get_matrix_benchmark_headspace_factor()`
already applies the matrix observable factor — never multiply `calibration_observable_factor`
again downstream.

### Layer 4 — Safety & Decision Support

**Modules:** `src/safety.py`, `src/sensory.py`, `src/usability_reports.py`, `src/reporting.py`

Maps chemistry into scientist-facing decisions. Exposes confidence, calibration diagnostics,
and reportable artifacts. The trust boundary is visible in every result surface.

### Layer 5 — Validation & Benchmark Surfaces

**Modules:** `src/benchmark_validation.py`, `src/validation_contract.py`

**Key artifacts:** `results/validation/benchmark_summary.md`, `results/validation/validated_envelope.md`

Defines what counts as benchmark-backed evidence. Keeps matrix evidence, internal
reproducibility harnesses, and strict-ready benchmarks separate.
Observable-first governance: never promote a target without a justifying artifact in
`results/validation/`.

### Layer 6 — Selective QM (xTB → DFT)

**Modules:** `src/mlp_barrier.py`, `src/dft_refiner.py`, `src/skala_refiner.py`

Refines barrier quality for selected high-value steps. Not the main reason the tool is
scientist-shareable today, but architecturally critical for the next tier of confidence.

---

## xTB: Role & Known Limitations

GFN2-xTB (Tier 1) identifies which pathways are kinetically viable out of thousands generated
by RMG. **It is a pathfinder, not a barrier authority.** Final barriers come from DFT
(r2SCAN/def2-svp + ddCOSMO water via PySCF/Sella) or explicit literature surrogates.

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Implicit solvation — no explicit water proton shuttles | Proton-transfer barriers overestimated 10–25 kcal/mol | Refine Amadori/enolisation steps with 1–2 explicit H₂O in DFT |
| Zwitterionic intermediates (amino acids at pH 5–8) | Schiff base energetics slightly distorted | DFT required for initial nucleophilic-attack steps |
| β-elimination (DHA pathway) | Barriers off by 5–15 kcal/mol | DHA elimination flagged for Tier 2 refinement |
| Open-shell sulfur radicals | Less accurate than closed-shell | MLP-first for radical-trapping lanes |

See the [Computational Gap Runbook](guides/COMPUTATIONAL_GAP_RUNBOOK.md) for the current
selective-QM queue and copy-paste commands.

---

## What Still Blocks State-of-the-Art Status

**Scientifically:**
- No primary quantitative pea or soy meaty-positive matrix benchmark with desirable and adverse targets in the same run
- No benchmark-closed time-series matrix data for target PPI/SPI systems
- Incomplete process-state calibration for real matrix release behaviour

**Architecturally:**
- No full experiment-ingestion path for internal wet-lab data
- No external-team API layer beyond scripts and report artifacts

The highest-value missing experiment is documented in
[protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md](protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md).

---

## Review Sequence for External Sharing

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh campaign data/campaigns/shareable_meaty_screen.yml
```

For the full validation contract see
[reference/VALIDATION_CONTRACT.md](reference/VALIDATION_CONTRACT.md).

# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Recommended](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Maillard is a computational screening framework designed for food scientists who want to predict and optimize meat-like Maillard chemistry in plant-based matrices (pea, soy, mycoprotein) *before* running expensive wet-lab campaigns.

Our philosophy: **A model is useful when it separates known bounds from structural gaps.** By combining determinisic kinetic ODEs with matrix-aware retention and headspace physics, Maillard helps answer: *Which formulation and process changes are most worth testing next if the goal is meaty aroma under plant-matrix constraints?*

For a deeper dive into the physics and framework rationale, see [docs/PHILOSOPHY.md](docs/PHILOSOPHY.md).

## What This Does (Vs. Existing Tools)

Unlike general databases (like NIST WebBook) that only store static spectra, or raw chemical kinetics solvers (like RMG) that assume free chemistry environments, Maillard specifically translates:
1. **Ingredients + Process** (precursors, protein matrix type, physical state, pH, temp, time)
2. **Into actionable, benchmark-calibrated aroma outputs** (meaty thiols vs. beany off-notes, retention effects, and safety limits).

```mermaid
graph LR
    A[Inputs: Sugars, AAs, Matrix, T, pH] -->|Maillard Engine\n(ODEs + Matrix Correction)| B[Predictions: Aroma Profile, Off-Notes]
    B -.-> C[Decision Dashboard:\nConfidence & Interventions]
    style A fill:#e6f3ff,stroke:#4a90e2,stroke-width:2px;
    style B fill:#e6ffea,stroke:#2ecc71,stroke-width:2px;
    style C fill:#fff3e6,stroke:#f39c12,stroke-width:2px;
```
For a detailed diagram of the pipeline's SMIRKS, thermodynamic gating, and retention layers, read [docs/architecture.md](docs/architecture.md).

## 🚀 Quick Start

This project requires conda or Docker.

**Install and Boot:**
```bash
# Bring up the validated Linux container
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
```

**Run a Screening Prediction:**
```bash
python scripts/run_pipeline.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --time-minutes 45 \
  --protein-type pea_iso \
  --target meaty \
  --minimize beany \
  --report \
  --output-dir results/first_run
```

**Optimize a Formulation:**
```bash
python scripts/optimize_formulation.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --target-tag meaty \
  --minimize-tag beany \
  --protein-type pea_iso \
  --n-iterations 50
```

## 📊 Example Output & Scientific Validation

Maillard intentionally distances itself from arbitrary "sensory scores". The `--report` output yields highly actionable **Decision-Risk Dashboards**:
- **Relative Intervention Results:** e.g., "Increasing pH shifted meaty thiols +20% but triggered acrylamide bounds."
- **Confidence Overlays:** Explicit bounding telling you whether the prediction relies on a known experimental anchor or is a purely directional extrapolation.

### What is Validated Today?
We track 16 distinct reaction families internally. Quantitative accuracy relies strictly on extracted literature benchmarks where exact headspace (ppb) measurements were taken under bounded conditions.

A snapshot of our current quantitative proof surface:
![Validation Overview](docs/assets/validation_overview.png)

And a separate snapshot of family coverage, which shows which chemistry lanes are already wired into the runtime even when they still have zero direct quantitative parity points:
![Family Coverage](docs/assets/family_coverage.png)

This README plot is intentionally benchmark-centric: it only shows quantitative measured-versus-predicted parity for literature systems with matched numeric compounds. It does not add every active computational-gap parametrization target as a plotted point.

The current no-wet-lab parametrization queue spans these formal selective-QM targets plus prep-only companion lanes:
- Family 11: `hexanal_radical_quench`
- Family 12: `lysinoalanine_crosslink`
- Family 13: `quinone_cys_michael`
- Family 13 safety lane: `asparagine_sugar_explicit_water_cluster`
- Family 14: `aa_ring_open_dicarbonyl`
- Family 15: `pe_schiff_base`, `pe_amadori` (now formal DFT candidates after managed xTB pair-gate success)
- Family 16: `melanoidin_radical_trapping`

In other words: the README parity plot is the quantitative proof surface, while the computational-gap queue above is the selective parametrization surface. Families can be active in the second surface without yet appearing as new points in the first.

### Where the next experiments matter most

The S20–S22 trust loop turns the panel residuals into ranked, bookable wet-lab requests. The heatmap below shows the per-(benchmark × compound) value-of-information score: yellow → red is "go run this experiment first"; cells marked `*` are predictions whose measured value sits **outside** the 90 % Monte-Carlo envelope (barriers + matrix/Henry/retention multipliers swept jointly).

![Experiment-value gap heatmap](results/validation/gap_heatmap.png)

Underlying artifacts (all regenerated in Docker):

- Monte-Carlo envelope per matched compound: [results/validation/prediction_uncertainty.md](results/validation/prediction_uncertainty.md)
- Ranked experiment requests (DOE template + intake YAML per top candidate): [results/validation/experiment_value_ranking.md](results/validation/experiment_value_ranking.md), [results/validation/experiment_requests/index.md](results/validation/experiment_requests/index.md)
- Leave-one-benchmark-out leverage on panel residual: [results/validation/loo_leverage.md](results/validation/loo_leverage.md)

Check [docs/VALIDATION_GUIDE.md](docs/VALIDATION_GUIDE.md) to explore the most current accuracy snapshot, execution coverage gaps, and explicit proof contracts.
To explore all machine-readable reports, benchmarking artifacts, and matrix coverage metrics generated by the framework, dive into the `results/validation/` directory.

## Documentation Hub
- **Installation and basic usage:** [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md)
- **Copy-paste runbook for the official computational-gap QM queue plus proxy lanes:** [docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md](docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md)
- **Literature gaps and recommended experiment requests:** [results/validation/literature_learning_loop.md](results/validation/literature_learning_loop.md)
- **Understanding terminology (FAST mode, bounds, etc):** [docs/guides/GLOSSARY.md](docs/guides/GLOSSARY.md)
- **Detailed scientific positioning and constraints:** [docs/PHILOSOPHY.md](docs/PHILOSOPHY.md)
- **Validation rules and benchmark testing:** [docs/VALIDATION_GUIDE.md](docs/VALIDATION_GUIDE.md)
- **Complete framework architecture:** [docs/architecture.md](docs/architecture.md)

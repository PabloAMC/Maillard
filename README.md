# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Calibration: 39/48 inside 90% CI](https://img.shields.io/badge/calibration-39%2F48%20inside%2090%25%20CI-brightgreen.svg)](results/validation/prediction_uncertainty.md)
[![16 reaction families](https://img.shields.io/badge/chemistry-16%20reaction%20families-blueviolet.svg)](data/lit/family_ingestion_plan.json)

**Maillard** is a computational screening framework that predicts which combinations of
sugars, amino acids, temperatures, pH, and protein matrices will produce **meat-like aroma**
in plant-based formulations — and which will produce **off-notes** or **safety concerns** —
*before* you run a single wet-lab experiment.

Plant-based meat smells "beany" because the Maillard reactions that create real meat aroma
behave very differently inside a soy or pea protein matrix than in a free aqueous system.
Maillard simulates 16 families of reaction chemistry, corrects for how the plant matrix traps
or releases each volatile, and tells you exactly where its predictions are trustworthy, where
they are directional, and **what experiment would improve them most**.

> **Who is this for?** Alternative protein scientists who want to triage formulations before
> burning GC-MS time, and computational chemists who want a transparent, benchmarked platform
> for matrix-aware Maillard chemistry.

---

## How it works

```mermaid
graph LR
    subgraph Inputs
        A["🧪 Precursors<br/>(sugars, amino acids, lipids)"]
        B["⚙️ Process conditions<br/>(pH, T, time, matrix)"]
    end

    subgraph Engine
        C["SMIRKS rule engine<br/>16 reaction families"]
        D["Kinetic solver<br/>Arrhenius + literature barriers"]
    end

    subgraph Physics
        E["Matrix correction<br/>binding, denaturation, SH access"]
        F["Headspace model<br/>Henry's law, retention"]
    end

    subgraph Output
        G["Per-compound ppb<br/>with 90% CI + ODT"]
        H["🔬 Ranked experiment<br/>requests (VoI)"]
    end

    A & B --> C --> D --> E --> F --> G --> H
```

**Inputs:** Choose your protein source (pea, soy, wheat gluten, mycoprotein, or 10 others),
precursors, and process conditions. **Engine:** The SMIRKS rule engine generates
atom-balanced reaction candidates across 16 chemistry families; a kinetic solver ranks them
using literature-calibrated Arrhenius parameters. **Physics:** Matrix corrections account for
protein-volatile binding, denaturation, and accessibility; a headspace model applies Henry's
law partitioning and retention penalties. **Output:** Per-compound predicted concentrations
(ppb) with 90% confidence intervals, odour-activity ratios, and a ranked list of which
wet-lab experiment would reduce uncertainty the most.

---

## What the output looks like

Before installing anything, here is what a typical prediction report contains:

| Compound              | Predicted ppb | 90% CI      | ODT (ppb) | Above ODT?    | Confidence                 |
| --------------------- | ------------- | ----------- | --------- | ------------- | -------------------------- |
| 2-methyl-3-furanthiol | 12.4          | 3.1 – 49.7 | 0.007     | ✅ ~1770×    | `bounded_calibration`    |
| 2-furfurylthiol       | 4.8           | 0.9 – 25.3 | 0.01      | ✅ ~480×     | `bounded_calibration`    |
| hexanal               | 284           | 71 – 1136  | 4.5       | ⚠️ off-note | `transferred_literature` |
| methional             | 0.8           | 0.1 – 6.4  | 0.2       | ✅ ~4×       | `surrogate_family`       |

> **How to read this:** *Predicted ppb* is the observable headspace concentration after matrix
> and retention corrections. *90% CI* is the Monte Carlo uncertainty envelope — narrow
> whiskers mean the prediction is well-anchored, wide ones mean the model needs more data.
> *ODT* is the odour detection threshold — compounds far above ODT dominate the sensory
> profile. *Confidence* shows the provenance tier of the underlying barrier constants.

The full report also includes a compound confidence overlay chart, an intervention waterfall,
and provenance metadata linking every number to its literature source.

![Compound confidence overlay](docs/assets/report_compound_confidence_overlay.png)

---

## How well calibrated is it?

We publish four orthogonal evidence surfaces rather than a single accuracy number.

### Headline: **39 out of 48** matched literature data points fall within the model's 90% confidence interval

<table>
<tr>
<td width="33%" valign="top"><a href="docs/assets/validation_overview.png"><img src="docs/assets/validation_overview.png" alt="Parity plot"/></a><br/><sub><b>Parity.</b> Predicted vs. measured ppb across 16 benchmarks.</sub></td>
<td width="33%" valign="top"><a href="docs/assets/family_coverage.png"><img src="docs/assets/family_coverage.png" alt="Family coverage"/></a><br/><sub><b>Coverage.</b> Which of the 16 families are wired and evidence-backed.</sub></td>
<td width="33%" valign="top"><a href="results/validation/gap_heatmap.png"><img src="results/validation/gap_heatmap.png" alt="Gap heatmap"/></a><br/><sub><b>Gaps.</b> Which experiments would close the largest blind spots.</sub></td>
</tr>
</table>

| Surface                       | Question                                                                 | Status                                                       |
| ----------------------------- | ------------------------------------------------------------------------ | ------------------------------------------------------------ |
| **Parity**              | On matched systems, how close is predicted ppb to measured?              | 16 benchmarks · 48 matched rows                             |
| **External hold-out**   | On systems excluded from calibration, does the frozen model still cover? | 4 bundles ·**0/8 inside 90% CI** · median 36× error |
| **Coverage**            | Which reaction families are wired and calibrated?                        | 16/16 wired · 7 with DFT anchors                            |
| **Experiment priority** | Where would the next experiment improve confidence the most?             | 9/48 cells outside 90% CI — all queued                      |

> **On the external hold-out (0/8):** These 8 data points are intentionally excluded from
> calibration. The large error (median 36×) quantifies how much structured plant-protein
> trapping is *not yet* captured by free-precursor kinetics alone. The gap heatmap converts
> each miss into a ranked, bookable wet-lab request.

### When to trust the predictions

| Trust level                                  | Use for                                            | Example                                                         |
| -------------------------------------------- | -------------------------------------------------- | --------------------------------------------------------------- |
| **High** — use freely                 | Ranking formulations within the benchmark envelope | Free-precursor sulfur chemistry (Cys + ribose, 80–150°C)      |
| **Moderate** — verify before deciding | Directional prioritisation                         | Pea vs soy matrix comparisons; choosing what to test next       |
| **Low** — exploratory only            | Hypothesis generation                              | New protein sources without nearby benchmarks; extrusion claims |

Full validation methodology: [VALIDATION_CONTRACT.md](docs/reference/VALIDATION_CONTRACT.md).
Regenerate all evidence artifacts: `./scripts/docker_maillard.sh summary`.

---

## The 16 reaction families

Maillard covers 16 families of chemistry, each independently calibrated and wired:

| #  | Family                        | Key compounds                     | Role                   | Status                   |
| -- | ----------------------------- | --------------------------------- | ---------------------- | ------------------------ |
| 01 | Amino acid–sugar core        | MFT, FFT, methional, pyrazines    | Core meaty aroma       | ✅ Benchmarked           |
| 02 | Lipid oxidation & crosstalk   | Hexanal, 2-pentylfuran, nonanal   | Off-note & competition | ✅ Benchmarked           |
| 03 | Thiamine degradation          | MFT (via thiamine), thiazoles     | Sulfur support         | 📐 Literature-calibrated |
| 04 | Nucleotide & ribose support   | IMP, GMP, umami precursors        | Umami/kokumi support   | 📐 Literature-calibrated |
| 05 | Glutathione & peptide support | GSH-derived sulfur volatiles      | Sulfur boost           | 📐 Literature-calibrated |
| 06 | Alternative protein matrices  | Matrix-specific modifiers         | Scope extension        | 🔄 Directional           |
| 07 | Carbonyl donor hierarchy      | Sugar reactivity ranking          | Formulation variable   | ✅ Benchmarked           |
| 08 | Off-notes & suppression       | Dicarbonyl traps, inhibitors      | Guardrail              | ✅ Benchmarked           |
| 09 | Caramelization & pyrolysis    | HMF, furfural                     | Severity flag          | 🔄 Directional           |
| 10 | Fermentation pretreatment     | Free AA/nucleotide enrichment     | Upstream modifier      | 📐 Literature-calibrated |
| 11 | Lipid–Maillard crosstalk     | MFT quenching by aldehydes        | Competition surface    | 🧮 DFT-queued            |
| 12 | Protein damage markers        | CML, CEL, furosine, lysinoalanine | Safety guardrail       | 🧮 DFT-queued            |
| 13 | Polyphenol–amino capping     | Quinone–thiol trapping           | Precursor sink         | 🧮 DFT-queued            |
| 14 | Ascorbic acid Maillard        | Dicarbonyl source, crosslinks     | Upstream modifier      | 🧮 DFT-queued            |
| 15 | PE stealth sugar sink         | Phospholipid glycation            | Sugar depletion        | 🧮 DFT-queued            |
| 16 | Melanoidin polymerization     | Thiol scavenging                  | Trapping burden        | 🧮 DFT-queued            |

✅ = benchmark-validated &emsp; 📐 = literature-calibrated priors &emsp;
🔄 = directional only &emsp; 🧮 = computational closure in progress (xTB → DFT queue)

Full literature basis: [docs/slr_benchmark_evaluation.md](docs/slr_benchmark_evaluation.md) (76 peer-reviewed papers evaluated against an 8-criterion quality checklist).

---

## Guiding experiments: what to measure next

Maillard doesn't just predict — it tells you which experiment would improve its predictions
the most, using a **Value of Information (VoI)** framework that ranks compounds by:
uncertainty × sensory impact.

![Gap heatmap](results/validation/gap_heatmap.png)

**How to read the heatmap:** Each cell is a benchmark × compound combination. Bright cells
are high-value: the model is uncertain *and* the compound matters for sensory quality.
Dark cells are already well-calibrated.

### The highest-value missing experiment

The single experiment that would close the most gaps is a **quantitative PPI/SPI meaty-positive
benchmark** — pea and soy protein isolates with ribose + cysteine, measuring both desirable
sulfur targets and adverse off-flavour markers in the same GC-MS run:

- Full protocol: [PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md](docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md)
- Ranked experiment requests: [experiment_value_ranking.md](results/validation/experiment_value_ranking.md)

---

## Getting started

### Prerequisites

- [Docker](https://www.docker.com/products/docker-desktop/) (recommended) or Python 3.12 with conda
- Git

### Install and boot

```bash
git clone https://github.com/PabloAMC/Maillard.git
cd Maillard
./scripts/docker_maillard.sh up && ./scripts/docker_maillard.sh bootstrap
```

### Step 1 — Predict & screen

Run a simulation for a custom formulation:

```bash
./scripts/docker_maillard.sh run "python scripts/run_pipeline.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
  --ph 5.5 --temp 105 --time-minutes 45 \
  --protein-type pea_iso \
  --target meaty --minimize beany \
  --report --output-dir results/first_run"
```

Open `results/first_run/report.md` for the full analysis.

### Step 2 — Optimise & compare

Auto-search the ingredient space over 50 Bayesian iterations:

```bash
./scripts/docker_maillard.sh run "python scripts/optimize_formulation.py \
  --sugars ribose,glucose --amino-acids cysteine,leucine \
  --target-tag meaty --minimize-tag beany \
  --protein-type pea_iso --n-iterations 50"
```

Or run the built-in quickstart comparison:

```bash
./scripts/docker_maillard.sh quickstart
```

![Comparison intervention waterfall](docs/assets/report_comparison_intervention_waterfall.png)

### Step 3 — Ingest lab results & close the loop

When your GC-MS results come back, feed them back to shrink the uncertainty whiskers:

```bash
./scripts/docker_maillard.sh ingest \
  --file path/to/my_results.csv \
  --protein-type pea_iso \
  --process-state extrusion_structured \
  --temp-c 105 --ph 5.5 \
  --precursor cysteine=15.0 \
  --precursor glucose=30.0 \
  --confirm
```

The engine updates the benchmark database and regenerates validation artifacts.

Full quickstart guide: [QUICKSTART.md](docs/guides/QUICKSTART.md). Glossary for non-computational scientists: [GLOSSARY.md](docs/guides/GLOSSARY.md).

---

<details>
<summary><b>Computational methods under the hood</b></summary>

### Reaction enumeration

The SMIRKS rule engine (`src/smirks_engine.py`) generates deterministic, atom-balanced reaction
candidates from precursor sets. The chemistry surface is explicit and inspectable — if the
reaction graph is not coherent, no later scoring is trustworthy.

### Kinetics

Literature-calibrated Arrhenius parameters (Schiff base 15 kcal/mol → enolisation 28 kcal/mol;
anchored to Hofmann, Martins, Nursten) drive laptop-speed FAST screening in < 1 second.
A Cantera ODE solver is available for temporal profiles and mechanism debugging.

### Matrix physics

Accessibility corrections model how protein denaturation, sulfhydryl availability, and pH
alter precursor reactivity in structured plant matrices. Volatile retention factors
(non-covalent binding to denatured proteins) and Henry's law headspace partitioning translate
beaker-chemistry yields into what a sensory panel would actually perceive.

### Uncertainty quantification

Monte Carlo propagation across the full benchmark panel produces per-compound 90% confidence
intervals. Each prediction carries a confidence tier label (`bounded_calibration`,
`transferred_literature`, `surrogate_family`, `xtb_derived`) that propagates through the
pipeline and surfaces in every report.

### Selective quantum chemistry (xTB → DFT)

For high-value rate-limiting steps, GFN2-xTB identifies kinetically viable pathways (pathfinder,
not a barrier authority). Final barriers come from DFT (r2SCAN-3c/def2-svp + ddCOSMO implicit
water via PySCF/Sella). See the [Computational Gap Runbook](docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md)
for the current queue and copy-paste commands.

### Key dependencies

RDKit (reaction enumeration), Cantera (ODE kinetics), PySCF + Sella (DFT refinement),
xtb (pathfinding), MACE (ML potentials), NumPy/SciPy (numerics), Matplotlib (figures).

### Supported protein matrices

14 protein sources are registered with endogenous composition profiles and matrix correction
factors: pea isolate, pea concentrate, soy isolate, soy concentrate, wheat gluten, mycoprotein,
chickpea, lentil, faba bean, oat, potato, sunflower, spirulina, and yeast extract.
Pea and soy have the strongest literature backing; others are directional.

</details>

---

## Where to look next

| If you are a…                                          | Start with                                                                                                                                 |
| ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **Food scientist** — first run                   | [QUICKSTART.md](docs/guides/QUICKSTART.md)                                                                                                    |
| **Scientist** — understanding the output         | [GLOSSARY.md](docs/guides/GLOSSARY.md)                                                                                                        |
| **Reviewer** — auditing what is verified         | [VALIDATION_CONTRACT.md](docs/reference/VALIDATION_CONTRACT.md) → [results/validation/](results/validation/)                                    |
| **Experimentalist** — closing the gaps           | [PPI_SPI protocol](docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md) → [experiment ranking](results/validation/experiment_value_ranking.md) |
| **Maintainer** — extending the chemistry         | [architecture.md](docs/architecture.md) → [SMIRKS_SYSTEM.md](docs/reference/SMIRKS_SYSTEM.md)                                                   |
| **QM operator** — running the DFT queue          | [COMPUTATIONAL_GAP_RUNBOOK.md](docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md)                                                                      |
| **Literature curator** — ingestion & calibration | [data/lit/README.md](data/lit/README.md)                                                                                                      |
| **New contributor**                               | [CONTRIBUTING.md](CONTRIBUTING.md)                                                                                                            |

---

## Citation

If you use Maillard in your research, please cite:

```
Moreno Casares, P. A. (2026). Maillard: Computational screening for meat-like
Maillard chemistry in plant-based protein matrices. GitHub repository.
https://github.com/PabloAMC/Maillard
```

## License

[Apache 2.0](LICENSE)

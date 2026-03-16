# Maillard Reactant Framework

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Required](https://img.shields.io/badge/Docker-Required-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Maillard** is a pure-Python, high-fidelity chemical discovery engine designed for the next generation of plant-based foods. It explores the high-dimensional chemical space of the Maillard reaction to help you design flavor systems that are indistinguishable from animal meat.

- 🌿 **Plant-Based Focus**: Mechanistic modeling of precursor accessibility and volatile retention in soy, pea, and fungal matrices.
- ⚡ **Multi-Tiered Screening**: A hierarchical funneling strategy ranging from combinatorial rule-based graph traversal to high-precision DFT refinement.
- 🔬 **Scientific Rigor**: Domain-specific physics for pH-dependent kinetics, water activity effects, and heme-mediated lipid/Maillard synergy.

---

## 🎯 Mission

To empower food scientists to rationally design precursor combinations (pea, soy, sugars, fats) that maximize meaty volatiles (MFT, pyrazines) while minimizing off-flavors (beany hexanal) and toxic by-products (HMF, acrylamide).

### 🌟 Highlights
- **Hybrid SmirksEngine**: Automated discovery of complex reaction cascades (Amadori, Strecker, caramelization) with strict stoichiometric mass conservation to ensure physical validity of derived kinetics.
- **Formulation Inverse Design**: Multi-objective Pareto optimization across formulation grids to identify precursor combinations that satisfy specific sensory targets.
- **Bayesian Formulation Optimization**: Active learning-driven search of continuous multi-dimensional spaces (precursor ratios, pH, T(t)) facilitated by `optuna`.
- **Sensory & Safety Radar**: Psychophysical mapping via Stevens' Power Law to correlate chemical flux with human perception, while simultaneously decoupling flavor from toxic marker formation (Acrylamide, HMF, AGEs).
- **Matrix Physics Engine**: Quantitative modeling of "precursor accessibility" and "headspace partitioning" (Henry's Law constants corrected for protein/fat binding) to account for the unique challenges of plant-based isolates and concentrates.

---

## 🧬 The Challenge

Plant-based proteins (pea, soy, etc.) lack the native precursor matrix of animal meat (Ribose, Cysteine, Heme). This leads to:
1.  **Aroma Gap**: Insufficient production of key meat odorants like 2-methyl-3-furanthiol (MFT).
2.  **Off-Flavors**: Dominance of beany/grassy notes from lipid oxidation.
3.  **Competition**: Stoichiometric competition from the Dehydroalanine (DHA) pathway consuming critical amino acids.

Exploring this space in the wet-lab is **combinatorially explosive**, **slow**, and **expensive**.

## 🧪 Scientific Case Studies & Reports

We provide detailed scientific reports that illustrate how to use the framework for real-world Alt-Protein formulation challenges. Each case study compares predicted outcomes against peer-reviewed literature.

- 🍔 **[Premium Roast Pea Protein Patty](docs/use_cases/pea_protein_report.md)**: Mitigating beany off-flavors by optimizing the Thiol:Aldehyde ratio through pH and sulfur supplementation.
- 🥜 **[Alkali-Induced Roasted Nutty Profile](docs/use_cases/roasted_nutty_report.md)**: Controlling pyrazine regioselectivity in alkaline plant-based beverage formulations.
- ☢️ **[Toxicity-Flavor Decoupling](docs/use_cases/toxicity_decoupling_report.md)**: Advanced Pareto optimization to resolve the temperature-dependent conflict between Maillard aroma and Acrylamide safety limits.
- 📑 **[Report Template](docs/use_cases/REPORT_TEMPLATE.md)**: Guidelines for contributing new scientific validations.

---

## 🛠️ How It Works

Maillard uses a funnel strategy: generate broadly, then refine precisely. Most users only need the top two tiers.

| Tier | What it does | Speed | When you need it |
|---|---|---|---|
| **Tier 0** | Generates reaction networks (SMIRKS + templates) | Seconds | **Always** — this is the core engine |
| **FAST** | Concentration-aware kinetic ranking + sensory prediction | Seconds | **Always** — ranks pathways by flavor impact |
| **ML (MACE-OFF24)** | Near-DFT activation barriers via machine learning | Minutes | When you need accurate barrier energies without HPC |
| **xTB / DFT** | Semi-empirical or full quantum chemistry | Hours | Research-grade validation of specific bottlenecks |

> [!TIP]
> **For most formulation work**, Tier 0 + FAST is all you need. The Bayesian optimizer uses these tiers internally and runs entirely on a laptop.

### Key Capabilities

- **Reaction Discovery**: Deterministic enumeration of Maillard cascades, Strecker degradations, and radical-mediated lipid oxidation via a domain-specific SMIRKS library.
- **Physicochemical Modeling**: Smooth sigmoid transitions for pH-dependent reactivity and non-isothermal Arrhenius time-integrals for complex temperature ramps.
- **Headspace & Partitioning**: Advanced calculation of the observable volatilome using Henry's Law corrected by matrix-specific protein and lipid binding coefficients.
- **Safety & Toxicity Decoupling**: Rigorous modeling of toxic byproduct kinetics (Acrylamide, CML, HMF) to provide clear safety boundaries for high-temperature processing.
- **Sensory Psychophysics**: Mapping chemical concentrations to human olfactory thresholds and intensities using power-law scaling across meaty, roasted, and vegetable-like dimensions. 

## 🗺️ Framework Capabilities Roadmap

The following scientific capabilities define the framework's scope and modular growth path. For real-time validation status against specific literature benchmarks, see the **[Scientific Validation Guide](docs/VALIDATION_GUIDE.md)**.

### Core Chemical Foundation
- **Reaction SmirksEngine**: Deterministic rule-based mechanism generation with atom-by-atom balancing.
- **Radical Lipid Oxidation**: Mechanistic propagation and β-scission modeling to predict "beany" off-flavor generation from PUFAs.
- **Temporal FAST Solver**: Rapid kinetic ranking using Boltzmann-weighted flux and Arrhenius temperature integration.
- **Formulation Design Engine**: Integrated grid search and Bayesian optimization (Optuna) for Pareto-optimal design.

### Advanced Food Matrix Physics
- **Matrix-Aware Kinetics**: Accessibility scaling for protein isolates (Pea, Soy) and concentrates to account for reactive site burial.
- **Dynamic Headspace Projection**: Multi-phase partitioning model considering native fat/protein binding constants.
- **Safety Conservatism**: Kinetic models for Acrylamide (Knol 2009) and HMF formation during processing.
- **Heme & Ion Catalysis**: Specialized kinetic modifiers for metal-promoted pyrazine and sulfur pathways.

### Future Research Frontiers
- **Flavor-Texture Coupling**: Integrating the DHA pathway's rheological feedback with flavor precursor consumption.
- **Mechanochemistry of Extrusion**: Modeling the impact of high-shear environments on precursor release and reaction barriers.
- **Phytochemical Sequestration**: Modeling how plant-derived polyphenols scavenge reactive osones to modulate flavor profiles.

## ⚙️ Installation & Setup

Maillard requires Python 3.12 and several specialized chemical informatics binaries (RDKit, Cantera, xTB). 

### 🍎 macOS (Apple Silicon M1/M2/M3) - Recommended
Because chemical binaries are often x86-specific, we recommend using **Docker** (via OrbStack or Docker Desktop) for a consistent environment:

1.  **Start the container**:
    ```bash
    ./scripts/docker_maillard.sh up
    ```
2.  **Bootstrap the environment** (one-time setup for dependencies):
    ```bash
    ./scripts/docker_maillard.sh bootstrap
    ```

### 🐧 Linux & 🪟 Windows (WSL2)
Use **Miniforge** to manage native dependencies:
```bash
conda env create -f environment.yml
conda activate maillard
# See Installation.md for mandatory patches for 'xtbiff' and 'pytorch'.
```

> [!NOTE]
> For a full list of manual patches and external binary requirements, see the **[Detailed Installation Guide](Installation.md)**.

---

## 🐳 Working with Docker (macOS/Workflow)

The `./scripts/docker_maillard.sh` script is your primary tool for day-to-day work. It ensures you are always running against the validated environment.

### Daily Lifecycle
- **Start**: `./scripts/docker_maillard.sh up`
- **Interactive Shell**: `./scripts/docker_maillard.sh shell` (activates `maillard` conda env automatically).
- **Stop**: `./scripts/docker_maillard.sh down`

### Verification Lanes
Run specific test lanes to verify your changes:
- `core`: Unit + integration correctness.
- `scientific`: FAST regressions and benchmark generation.
- `qm-heavy`: Quantum chemistry and external-backend validation.

Example:
```bash
./scripts/docker_maillard.sh core
```

---

## 🔬 Scientific Validation & Results

We use **Test-Driven Science** to monitor our correlation with literature. For a deep dive into our methodology, see the **[Scientific Validation Guide](docs/VALIDATION_GUIDE.md)**.

### Generating Reports
- **Benchmark Summary**: `./scripts/docker_maillard.sh summary` — Writes `results/validation/benchmark_summary.md`.
- **Benchmark Index**: `./scripts/docker_maillard.sh index` — Maps execution paths and matrix support status.
- **Targets Report**: `./scripts/docker_maillard.sh targets-report` — Aggregates all predicted volatiles.

### Inspecting Specific Benchmarks
View predicted targets for a specific study:
```bash
./scripts/docker_maillard.sh targets data/benchmarks/cys_glucose_150C_Farmer1999.json
```


### Execution Gates
To enforce the high-fidelity scientific gate during development:
```bash
MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py
```
This requires Pearson `>= 0.85` and error ratios `<= 1.5x` for all PRIMARY systems.


---

## 🛠️ Usage

### 1. Python API Quickstart (For Food Scientists)
Maillard is designed to be easily scriptable in Jupyter Notebooks or standard Python workflows. Here is how you run a Bayesian formulation optimization to find the perfect mix of ingredients:

```python
from src.bayesian_optimizer import FormulationOptimizer

# 1. Define your goal: Maximize "meaty" notes, mask "beany" off-flavors
optimizer = FormulationOptimizer(
    target_tag="meaty", 
    minimize_tag="beany", 
    risk_aversion=1.5 # Penalize Acrylamide/HMF by 1.5x
)

# 2. Define your available ingredients (e.g., from a pea protein matrix)
sugars = ["ribose", "glucose"]
amino_acids = ["cysteine", "leucine"]
lipids = ["hexanal"] # Source of the beany off-flavor

# 3. Optimize! Optuna will search the concentration, pH, and temp space
study = optimizer.optimize(
    fixed_sugars=sugars,
    fixed_amino_acids=amino_acids,
    fixed_lipids=lipids,
    n_trials=25
)

best = study.best_trial
print(f"Best Target Score: {best.user_attrs['target_score']:.2f}")
print(f"Optimal pH: {best.params['ph']:.2f}")
print(f"Optimal Temp: {best.params['temp']:.1f} °C")
```

### 2. Command Line Interface (CLI)
Identify the precursors and tags available in your current database:
```bash
python scripts/run_pipeline.py --list-precursors
python scripts/run_pipeline.py --list-tags
```

### 3. Forward Mode: Predict Aroma
Predict the volatiles and sensory profile produced by a formulation.
```bash
python scripts/run_pipeline.py \
    --sugars ribose:0.5 \
    --amino-acids cysteine:0.2,leucine:0.1 \
    --ph 5.5 \
    --temp 105 \
    --protein-type pea_iso
```
*Note: Add `--xtb` for rigorous (but slow) structural optimization.*

### 4. Bayesian Optimizer
Search the continuous space to find the Pareto-optimal formulation for flavor vs. safety.
```bash
python scripts/optimize_formulation.py \
    --sugars ribose,glucose \
    --amino-acids cysteine,leucine \
    --target-tag meaty \
    --minimize-tag beany \
    --n-iterations 50
```

### 5. Kinetics & Validation: Simulated vs Experimental Yields
Run rigorous ODE-based microkinetic simulations (supporting temperature ramps) and validate against experimental benchmarks.
```bash
# Run simulation using the structured results database
python scripts/run_cantera_kinetics.py \
    --precursors ribose:0.1,glycine:0.1 \
    --temp-ramp data/temp_profiles/isothermal_150.csv \
    --input results/maillard_results.db \
    --predict-sensory

# Validate framework against literature benchmarks
python scripts/compare_sim_to_lit.py
```
*Use `--export mech.yaml` to save the Cantera mechanism for external use.*

---

## 📂 What You Get Back

Every evaluation (whether from the CLI, Python API, or Bayesian optimizer) returns:

| Output | What it tells you |
|---|---|
| **Sensory Radar** | Multi-axis flavor profile (meaty, roasted, beany, malty, earthy) scaled by Stevens' power law |
| **Target Score** | How well this formulation hits your desired flavor tag (e.g., "meaty") |
| **Safety Score** | Penalty from predicted toxic marker formation (Acrylamide, CML, HMF) |
| **Flagged Toxics** | Specific compounds flagged as safety risks for this formulation |
| **Off-Flavour Risk** | Predicted intensity of undesirable notes (e.g., beany/grassy) |
| **Lipid Trapping %** | How effectively your amino acids sequester reactive aldehydes (like hexanal) |
| **Lysine Budget** | % of lysine consumed by the competing Dehydroalanine (DHA) pathway |

The Bayesian optimizer additionally tracks the full optimization trajectory so you can inspect how it converged on the optimal formulation.

<details>
<summary><strong>🔬 Advanced: DFT Barrier Refinement (for computational chemists)</strong></summary>

If you need research-grade activation barriers, the framework supports a full quantum chemistry pipeline:

```bash
# Generate 3D geometries for reactants/products
python scripts/generate_mapped_geometries.py

# Run xTB transition state search
python scripts/run_tier2_dft.py --reaction strecker

# With IRC validation
python scripts/run_tier2_dft.py --reaction strecker --irc
```

This uses the `r2SCAN-3c // wB97M-V` composite protocol. Requires `pyscf`, `geometric`, and optionally `CREST` for explicit solvation. See `src/skala_refiner.py` for details.

</details>

---

## 🧩 Architecture: A Codebase Tour

If you are new to the project, here is how the core modules plug together to build the simulation:

### 1. Generative Chemistry (`src/`)
- **`smirks_engine.py`**: The heart of the network generator. Applies reaction SMIRKS templates to discover thousands of possible Maillard pathways while enforcing strict stoichiometric mass conservation.
- **`conditions.py`**: Defines the physical `ReactionConditions` (pH, temperature, water activity). Enforces physical kinetics using smooth sigmoid transitions instead of hard cutoffs.

### 2. Physical & Quantum Chemistry (`src/`)
- **`results_db.py`**: SQLite caching layer so we never compute the same activation barrier twice.
- **`mlp_barrier.py`**: Leverages the **MACE-OFF24** ML potential for near-DFT accurate barrier approximations in milliseconds.
- **`xtb_screener.py`** & **`skala_refiner.py`**: The heavy-duty quantum chemistry layers. Used to resolve unknown pathway bottlenecks via Semi-Empirical (GFN2-xTB) or full DFT (r2SCAN-3c) calculations.

### 3. Food Science & Sensory Prediction (`src/`)
- **`recommend.py`**: The `FAST` kinetic solver. Ranks active flavor pathways by resolving rate-limiting bottlenecks against reactant concentrations.
- **`headspace.py`**: Corrects liquid concentrations into air-phase (headspace) concentrations. Accounts for matrix effects (hydrophobic flavors getting trapped in plant fats/proteins).
- **`sensory.py`**: Translates chemical concentrations into human perception using Stevens' Power Law, generating multi-dimensional flavor radar charts.

### 4. Formulation Design (`src/` & `scripts/`)
- **`inverse_design.py`**: Evaluates static grids of formulations, applying Pareto-ranking to balance desired flavor profiles against safety risks (like Acrylamide or HMF).
- **`bayesian_optimizer.py`** / `scripts/optimize_formulation.py`: Uses `optuna` to actively search the continuous multi-dimensional space (varying pH, time, temperatures, and exact ingredient ratios) to find the absolute mathematically optimal formulation.
- **`scripts/run_cantera_kinetics.py`**: Exports the discovered network to Cantera for rigorous, time-dependent ODE temperature-ramp simulations.

## ⚖️ License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

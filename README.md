# Maillard Reactant Framework

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Maillard** is a pure-Python, high-fidelity chemical discovery engine designed for the next generation of plant-based foods. It explores the high-dimensional chemical space of the Maillard reaction to help you design flavor systems that are indistinguishable from animal meat.

- 🌿 **Plant-Based Focus**: Tailored rules for soy, pea, and fungal protein precursors.
- ⚡ **Multi-Tiered Screening**: Balances generative breadth with DFT precision.
- 🔬 **Scientific Rigor**: Includes native support for pH, water activity, and heme catalysis.

---

## 🎯 Mission

To empower food scientists to rationally design precursor combinations (pea, soy, sugars, fats) that maximize meaty volatiles (MFT, pyrazines) while minimizing off-flavors (beany hexanal) and toxic by-products (HMF, acrylamide).

### 🌟 Highlights
- **Hybrid SmirksEngine**: Automated discovery of thousands of pathways with strict mass conservation.
- **Formulation Inverse Design**: Don't just predict; search the formulation grid for the optimal precursor matrix to hit a target sensory profile.
- **Bayesian Formulation Optimization**: Actively learn and search a continuous space of precursor concentrations, pH, and temperature using `optuna`.
- **Sensory & Safety Radar**: Target flavor profiles (e.g., meaty) using Stevens' psychophysical power-law scaling while strictly penalizing toxic markers (Acrylamide, HMF) via Pareto ranking.
- **PBMA Metrics**: Native calculation of "Lysine Budgets" and "Lipid Trapping Efficiency" to account for competition in complex food matrices.

---

## 🧬 The Challenge

Plant-based proteins (pea, soy, etc.) lack the native precursor matrix of animal meat (Ribose, Cysteine, Heme). This leads to:
1.  **Aroma Gap**: Insufficient production of key meat odorants like 2-methyl-3-furanthiol (MFT).
2.  **Off-Flavors**: Dominance of beany/grassy notes from lipid oxidation.
3.  **Competition**: Stoichiometric competition from the Dehydroalanine (DHA) pathway consuming critical amino acids.

Exploring this space in the wet-lab is **combinatorially explosive**, **slow**, and **expensive**.

## 🧪 Scientific Case Studies & Reports

We provide detailed scientific reports that illustrate how to use the framework for real-world Alt-Protein formulation challenges. Each case study compares predicted outcomes against peer-reviewed literature.

- 🍔 **[Premium Roast Pea Protein Patty](docs/use_cases/pea_protein_report.md)**: Strategies for masking beany off-flavors while maximizing meaty thiols (MFT/FFT).
- 🥜 **[Alkali-Induced Roasted Nutty Profile](docs/use_cases/roasted_nutty_report.md)**: Optimizing pyrazine formation for plant-based milks and beverages.
- ☢️ **[Toxicity-Flavor Decoupling](docs/use_cases/toxicity_decoupling_report.md)**: Balancing high-heat searing flavor with Acrylamide/HMF safety limits.
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
- **Reaction Discovery**: Automated enumeration of Maillard, Strecker, Amadori, retro-aldol, and thiol pathways with strict mass conservation.
- **pH & Temperature Physics**: Smooth sigmoid kinetics model how pH shifts favor different pathways (acidic → furans/thiols, alkaline → pyrazines).
- **Headspace Partitioning**: Converts liquid-phase concentrations to what the consumer actually smells, accounting for fat/protein binding in plant matrices.
- **Safety Scoring**: Automatically flags and penalizes toxic marker formation (Acrylamide, CML, CEL, HMF) using Pareto ranking.
- **Sensory Radar**: Stevens' power-law psychophysical model generates multi-axis flavor profiles (meaty, roasted, beany, malty, earthy).

## 🚀 Completed Milestones & Roadmap

### 🟢 Phase 1: Complexity & Matrix (Done)
- [x] **Radical Lipid Oxidation**: Modeling the *generation* of beany off-flavors (hexanal) from PUFAs via `src/lipid_oxidation.py`.
- [x] **Matrix Correction Layer**: Scaling reactivity for protein-bound precursors via `src/matrix_correction.py`.
- [x] **Enzymatic Pre-Processing**: Simulating fermentation/hydrolysis clean-up via `src/pre_processor.py`.

### 🟡 Phase 2: Advanced Physics (Current Priority)
- [ ] **Temporal FAST Mode**: Integrating temperature ramps into the Boltzmann scoring loop.
- [ ] **Heme/Iron Catalysis**: Explicit kinetic models for transition-metal promoted oxidation.
- [ ] **Rheological Feedback**: Linking DHA cross-linking to physical texture changes.

## 🚀 Installation

### 1. Recommended Setup (Conda / Mamba)
The Maillard framework relies on complex scientific binaries (`CREST`, `xTB`, `PySCF`) which are most reliably managed via Conda.

**For Linux & Windows (WSL2):**
```bash
# Create the unified environment
conda env create -f environment.yml

# Activate it
conda activate maillard
```

**For macOS (Apple Silicon: M1/M2/M3):**
Do not attempt to run this natively on macOS. Instead, use Docker for a seamless Linux environment:
```bash
# Start a new container (First time)
docker run --platform linux/amd64 -it -v "$(pwd):/workspace" -w /workspace condaforge/miniforge3

# Inside the container, set up environment:
# conda create -n maillard python=3.12 -y && conda activate maillard
```

**Returning to Work (macOS):**
```bash
docker start -ai maillard_container
conda activate maillard
```

> [!NOTE]
> For detailed per-OS setup (Windows WSL2, Linux native) and chemistry library patching (xtbiff), please refer to the [Installation Guide](Installation.md).


### 2. Verify Scientific Dependencies
Ensure that the QM engines are correctly detected by the framework:
```bash
# Check if binaries are in your PATH
which crest
which xtb

# Run core validation tests
python -m pytest tests/qm/test_solvation.py
```

### 3. Install Skala (Tier 2 DFT)
If you intend to run Tier 2 DFT refinement, install Microsoft Skala:
```bash
pip install git+https://github.com/microsoft/skala.git
```

---

## 🔬 Scientific Accuracy & Monitoring

The framework uses a **Test-Driven Science** approach. We maintain specific tests in `tests/scientific/` that monitor our correlation with literature and document known gaps. For a detailed breakdown of how we verify our predictions against curated literature benchmarks, see the **[Scientific Validation Guide](docs/VALIDATION_GUIDE.md)** and the **[Literature Benchmark Reference](data/benchmarks/maillard_validation_benchmarks.md)**.

### 🚩 Known Blind Spots (Tracked)
We proactively document and test for current engine limitations to prevent over-confidence in edge cases:
- **Heme Optimization**: While supported in the CLI, explicit leghemoglobin-specific kinetics are still being refined.
- **Supplier Variability**: Batch-level differences between isolate suppliers (PURIS vs Roquette) are not yet modeled.
- **Metal Catalysis**: General iron/copper synergistic effects on pyrazine formation are currently heuristic.

*Run these baseline tests with:* `python -m pytest tests/scientific/test_blind_spots.py`

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
List the supported sugars, amino acids, and lipids in the framework:
```bash
python scripts/run_pipeline.py --list-precursors
```

### 3. Forward Mode: Predict Aroma
Predict the volatiles produced by a specific precursor combination.
```bash
python scripts/run_pipeline.py \
    --sugars ribose \
    --amino-acids cysteine:2.0,leucine:1.0 \
    --catalyst heme \
    --ph 6.5 \
    --temp 160
```
*Note: Use `--xtb` to run rigorous structural optimizations (slow).*

### 4. Bayesian Formulation Optimization
Instead of evaluating a static grid, dynamically search the continuous parameter space (concentrations, pH, temp) to find the absolute Pareto-optimal formulation for flavor vs safety.
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

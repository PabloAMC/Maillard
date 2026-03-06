# Maillard Reaction Computational Framework

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Accelerating wet-lab experimentation for plant-based meat alternatives by computationally exploring the chemical pathways of the Maillard reaction.

## 🎯 Mission

Enable food chemists and formulation scientists to rationally design precursor combinations that maximize desirable "meaty" volatiles while minimizing off-flavors and toxic by-products.

## 🧬 The Challenge

Plant-based proteins (pea, soy, etc.) lack the native precursor matrix of animal meat (Ribose, Cysteine, Heme). This leads to:
1.  **Aroma Gap**: Insufficient production of key meat odorants like 2-methyl-3-furanthiol (MFT).
2.  **Off-Flavors**: Dominance of beany/grassy notes from lipid oxidation.
3.  **Competition**: Stoichiometric competition from the Dehydroalanine (DHA) pathway consuming critical amino acids.

Exploring this space in the wet-lab is **combinatorially explosive**, **slow**, and **expensive**.

## 🛠️ Computational Architecture

This framework employs a multi-tiered approach to bridge the gap between screening breadth and physical accuracy:

| Tier | Method | Scope | Purpose |
|---|---|---|---|
| **Tier 0** | **Generative SMIRKS Engine** | Thousands of pathways | Automated discovery & mass conservation |
| **Tier 1** | **Semi-Empirical (GFN2-xTB)** | Hundreds of steps | Rapid energetic screening & pathway ranking |
| **Tier 2** | **DFT (r2SCAN-3c / wB97M-V)** | Decisive bottlenecks | High-precision activation barriers |

### Key Technological Pillars
- **Hybrid SmirksEngine**: A rule-based generative engine combining SMARTS functional group transforms with expert-parameterized templates (Strecker, Amadori, etc.) to ensure strict mass conservation.
- **Physical Parametrization**: Native support for **pH**, **Temperature**, and **Water Activity ($a_w$)** affects pathway bifurcations and kinetics.
- **Boltzmann scoring**: Concentration-dependent ranking of formulations based on active flux predictions rather than binary presence/absence.

## 🚀 Getting Started

### Prerequisites
- Conda or Mamba
- RDKit, xTB (via conda-forge)
- PySCF, Skala (for Tier 2)

### Installation
```bash
conda env create -f environment.yml
conda activate maillard
python scripts/check_env.py
```

### Basic Usage

#### Forward Mode (Predict Volatiles)
Predict the aroma profile of a specific formulation:
```bash
python scripts/run_pipeline.py --sugars ribose --amino-acids cysteine,leucine --catalyst heme
```

#### Inverse Design Mode (Optimize Formulation)
Search the `formulation_grid.yml` for the best precursor combination to hit a target sensory profile:
```bash
python scripts/run_pipeline.py --target meaty --minimize beany
```

## 📊 Current Status (Phase 8.E Complete)

Internal milestone reached: The core generative engine is fully atom-balanced and calibrated against literature baselines.
- [x] **Lipid synergy** (alkylthiazole formation) modeled.
- [x] **Concentration sensitivity** (Boltzmann scoring) implemented.
- [x] **Literature Validation Gate** passed (benchmarked against Ribose/Cys/Leu systems).
- [ ] **Phase 3.3 (Ongoing)**: High-precision DFT refinement for rate-limiting bottlenecks.

## 📂 Project Structure

- `src/`: Core logic (engine, recommender, screeners).
- `data/`: Species databases, sensory tags, and formulation grids.
- `docs/`: Technical architecture, pathway analysis, and theoretical limitations.
- `tests/`: Comprehensive unit and integration test suite (ensuring mass conservation & chemistry rules).

## ⚖️ License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

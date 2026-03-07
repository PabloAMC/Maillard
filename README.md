# Maillard Reactant Framework 🍳

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Maillard** is a high-fidelity chemical discovery engine designed for the next generation of plant-based foods. It explores the high-dimensional chemical space of the Maillard reaction to help you design flavor systems that are indistinguishable from animal meat.

- 🌿 **Plant-Based Focus**: Tailored rules for soy, pea, and fungal protein precursors.
- ⚡ **Multi-Tiered Screening**: Balances generative breadth with DFT precision.
- 🔬 **Scientific Rigor**: Includes native support for pH, water activity, and heme catalysis.

---

## 🎯 Mission

To empower food scientists to rationally design precursor combinations that maximize meaty volatiles (MFT, pyrazines) while minimizing off-flavors and toxic by-products like HMF.

### 🌟 Highlights
- **Hybrid SmirksEngine**: Automated discovery of thousands of pathways with strict mass conservation.
- **Formulation Inverse Design**: Don't just predict; search the formulation grid for the optimal precursor matrix to hit a target sensory profile.
- **PBMA Metrics**: Native calculation of "Lysine Budgets" and "Lipid Trapping Efficiency" to account for competition in complex food matrices.

---

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

## 🚀 Installation

### 1. Prerequisites (Conda)
Maillard requires RDKit and xTB for screening. We recommend using a dedicated environment:

```bash
conda env create -f environment.yml
conda activate maillard
```

### 2. Verify Scientific Dependencies
Ensure that the QM engines are correctly mapped:
```bash
# Check xTB and RDKit installation
python scripts/check_env.py
```

### 3. Install Skala (Tier 2 DFT)
If you intend to run Tier 2 DFT refinement, install Microsoft Skala:
```bash
pip install git+https://github.com/microsoft/skala.git
```

---

## 🛠️ Usage

### 1. Identify Available Precursors
List the supported sugars, amino acids, and lipids in the framework:
```bash
python scripts/run_pipeline.py --list-precursors
```

### 2. Forward Mode: Predict Aroma
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

### 3. Inverse Design Mode: Optimize Formulation
Generate the optimal precursor matrix to maximize a specific flavor profile.
```bash
python scripts/run_pipeline.py --target meaty --minimize beany
```

### 4. Advanced: Tier 2 DFT Refinement
High-precision activation barriers using the `r2SCAN-3c // wB97M-V` composite protocol.

> [!IMPORTANT]
> **Prerequisite:** You must generate the starting geometries and xTB transition state guesses before running the final DFT refinement.

**Step A: Generate Atom-Mapped Geometries**
This builds the initial 3D structures for both reactant and product for all 8 target reactions.
```bash
python scripts/generate_mapped_geometries.py
```

**Step B: Generate xTB Transition State Guess**
Navigate to the reaction directory and run the path search.
```bash
cd data/geometries/xtb_inputs/strecker
./run_xtb.sh  # Output: xtbpath_ts.xyz
```

**Step C: Run Tier 2 DFT Refinement**
Execute the orchestration script to refine the xTB guess into a high-level DFT barrier.
```bash
# Return to root
cd ../../../..
python scripts/run_tier2_dft.py --reaction strecker

# Use --fast for rapid testing (HF/STO-3G in Vacuum)
python scripts/run_tier2_dft.py --reaction strecker --fast

# Use --irc to perform automated Phase 3.4 validation (recommended)
python scripts/run_tier2_dft.py --reaction strecker --irc
```

> [!TIP]
> **Performance:** The refiner dynamically detects and utilizes all available CPU cores (via `os.cpu_count()`) to maximize throughput. On Apple M-series chips, this typically yields a 10x speedup by engaging all performance cores automatically.

---

## 📂 Output

Results are displayed in a formatted terminal table and include:

- **Predicted Compound**: The target volatile species (e.g., thiols, pyrazines).
- **Formation Tag**: Classification of the compound (Desirable, Toxic, etc.).
- **Barrier**: The rate-limiting energetic span for the pathway.
- **Sensory Character**: Qualitative descriptor (e.g., "savory, meaty").
- **Toxicity/Risk**: Metadata regarding known by-product risks (e.g., HMF, LAL).

### PBMA Metrics
The pipeline also calculates proprietary formulation indices:
- **Lipid Trapping Efficiency**: % of reactive aldehydes successfully sequestered by amino acid traps.
- **Lysine Budget (DHA)**: % of available lysine consumed by the competing Dehydroalanine pathway.

---

## 🧩 Architecture

- **`src/smirks_engine.py`**: Rule-based reaction network generator.
- **`src/xtb_screener.py`**: GFN2-xTB energy evaluation and NEB optimization.
- **`src/skala_refiner.py`**: High-precision DFT refinement using the Skala XC functional.
- **`src/recommend.py`**: Sensory scoring and pathway bottleneck identification.
- **`src/inverse_design.py`**: Automated search over formulation matrices.

---

## 📊 Project Roadmap (Phase 8.E Complete)

Internal milestone reached: The core generative engine is fully atom-balanced and calibrated against literature baselines.
- [x] **Lipid synergy** (alkylthiazole formation) modeled.
- [x] **Concentration sensitivity** (Boltzmann scoring) implemented.
- [x] **Literature Validation Gate** passed (benchmarked against Ribose/Cys/Leu systems).
- [ ] **Phase 3.3 (Ongoing)**: High-precision DFT refinement for rate-limiting bottlenecks.

## ⚖️ License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

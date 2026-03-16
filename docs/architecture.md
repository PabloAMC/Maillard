# Maillard Reaction Computational Framework
## Problem Statement, Proposed Solution, and Initial Architecture

> **Mission**: Accelerate wet-lab experimentation for plant-based meat alternatives by computationally exploring the chemical pathways of the Maillard reaction, enabling GFI scientists and food chemists to rationally design precursor combinations that maximise desirable "meaty" volatiles while minimising off-flavours and toxic by-products.

---

## 1. The Problem

### 1.1 Scientific Background

The Maillard reaction is a cascade of hundreds of interconnected, competing, non-enzymatic browning reactions between reducing sugars and amino groups. In animal meat, the native biological matrix—rich in ribose (from ATP breakdown), sulfur amino acids (cysteine, methionine), glutathione, and thiamine—spontaneously produces the signature meaty volatiles (2-methyl-3-furanthiol, 2-furfurylthiol, pyrazines, Strecker aldehydes) upon cooking.

Plant-based proteins (pea, soy, faba bean, wheat gluten) **lack** this native precursor matrix. They are deficient in:
- Sulfur-containing amino acids → no spontaneous S-Maillard pathway
- Open-chain pentose sugars (ribose/xylose) → slower, less reactive cascades
- Catalytic heme iron → absence of the pyrazine-boosting oxidation catalyst

Instead, they generate undesirable volatiles (hexanal, nonanal, 1-octen-3-ol, 2-pentylfuran) from PUFA oxidation, and compete for lysine through the Dehydroalanine (DHA) pathway, which simultaneously builds texture (lysinoalanine cross-links) but starves the Maillard network of nitrogen donors.

### 1.2 The Challenge for Formulation Scientists

Wet-lab experiments exploring precursor combinations are:
- **Combinatorially explosive**: dozens of amino acids × sugars × pH × temperature × water activity
- **Slow and expensive**: GC-MS volatilome characterisation per condition takes days
- **Difficult to interpret**: side pathways (DHA, AGE formation, PUFA oxidation) confound results

A computational tool that can **screen and rank candidate precursor systems** before committing to lab work would provide substantial leverage.

---

## 2. Target Outputs, Success Criteria, and Usefulness

The tool should predict, for a given input formulation (amino acid/peptide composition, reducing sugar, pH, T, aw):

| Output | Description | Status |
|--------|-------------|:---:|
| **Off-flavour risk** | Mechanistic enzymatic cleaning + trapping metric | Supported |
| **Competing pathway load** | DHA / lysinoalanine formation consuming lysine | Supported |
| **Toxicity flags** | AGE (CML, CEL) and HAA (PhIP, MeIQx) risk | Supported |
| **Concentration Sensitivity** | Bayesian Optimization of formulated precursors | Supported |
| **Reaction Physics** | Radical propagation & termination cycle | Supported |


**Headline success criterion**: Reducing the number of wet-lab conditions by ≥10× while maintaining experimental hit rate.

### 2.1 What "Useful" Concretely Looks Like

A formulation scientist at a pea-protein company wants to hit a *savory, meaty* profile without a *beany* off-note. They have:
- Access to 5 amino acid supplements, 2 reducing sugars, a pH dial, and a temperature profile.
- A GC-MS instrument that takes 3 days and £1,500 per sample.

The tool is **useful** if it can:

1. **Correctly predict winners and losers** within a set of formulations they'd test anyway.
   - *Example*: Given {ribose, cysteine, pH 5} vs {glucose, glycine, pH 7}, the tool correctly predicts that the former gives higher FFT (sulfur, roasted) and the latter gives more pyrazines (nutty, roasted).
   - *This is verified through the scientific validation suite.*

2. **Reveal non-obvious trade-offs** that they would not have intuited.
   - *Example*: "Adding hexanal (lipid) to a cysteine+leucine formulation actually increases alkylthiazoles via Strecker catalysis, not just masking."
   - *This leverages the Lipid-Maillard synergy pathways.*

3. **Give concentration guidance** beyond binary present/absent.
   - *Example*: "You need ≥0.3% cysteine to shift from pyrazine-dominant to FFT-dominant at pH 5."
   - *This requires precursor concentration resolution and Arrhenius-weighted flux models.*

The tool is **not useful** if it merely confirms what any experienced food chemist already knows from memory.

### 2.2 Current Validated Envelope (2026-03-15)

The architecture above describes the intended computational stack, but the currently validated operating envelope is narrower and now documented explicitly.

- The canonical validation path is Docker-first, using the `maillard` conda environment through `./scripts/docker_maillard.sh`.
- Three PRIMARY free-amino-acid benchmarks are currently strict-ready: `cys_glucose_150C_Farmer1999`, `cys_ribose_140C_Hofmann1998`, and `cys_ribose_150C_Mottram1994`.
- Two PRIMARY matrix benchmarks are now executable through a dedicated intake path rather than a free-precursor route: `pea_isolate_40C_PratapSingh2021` and `soy_isolate_40C_PratapSingh2021` run as `matrix_only`, with full benchmark coverage, but are not yet part of the strict release gate.

This distinction matters: the framework is already useful for ranking and benchmarking within the validated free-precursor envelope, and it now has a small executable matrix-headspace intake family, but it is not yet a generally validated matrix-headspace replication engine.

### 2.3 Risk Mitigation

> [!CAUTION]
> The biggest risk is building a confident-looking tool that doesn't correlate with experiment. A tool that produces wrong rankings is worse than no tool at all, because it erodes trust and wastes wet-lab resources.

**Primary risk**: The scoring function (`score = Σ max(0, 40 − barrier)`) produces near-identical scores for all formulations containing the same reaction families. A scientist who tests the "top" formulation and sees no differentiation from the second-best will immediately distrust the tool.

**Mitigations**:

1.  **Rigorous Literature Validation Gate**: The framework's core chemistry is governed by a "Hard Gate" standard. Before new reaction families or physical models are integrated into the primary simulation, they must demonstrate the ability to reproduce measured outcomes from high-fidelity, published model systems (e.g., the datasets of Hofmann, Mottram, and Farmer). This prevents the tool from drifting into "internally consistent" but experimentally invalid predictions.
2.  **Thermodynamic Boltzmann Scoring**: To provide the necessary sensitivity for formulation differentiation, the framework implements a Boltzmann-weighted flux model: `Σ [conc] ⋅ exp(−ΔG‡/RT)`. This ensures that small differences in activation barriers (ΔG‡) result in exponential changes in predicted concentration, allowing researchers to distinguish between "good" and "great" ingredient ratios.
3.  **Automated Scientific Regression**: Every architectural update or chemical rule modification is automatically re-run against the library of primary literature benchmarks. This continuous validation loop ensures that improvements in one flavor branch (e.g., pyrazines) do not unintentionally regress the accuracy of another (e.g., meaty thiols).

---

## 3. Proposed Computational Architecture

The framework has three tiers of increasing physical fidelity and computational cost, mirroring the tiered approach used in computational drug discovery.

```
┌──────────────────────────────────────────────────────────────┐
│  TIER 0: Pathway Enumeration (Rule-Based)                    │
│  (seconds per query)                                         │
│  • SmirksEngine: Hybrid SMIRKS + Parametric Templates        │
│  • Stoichiometrically balanced graph traversal               │
│  • Radical Autoxidation & Enzymatic Pre-processing           │
│  Output: Enumerated reaction graph, balanced chemical network │
└────────────────────┬─────────────────────────────────────────┘
                     │ Complete network with family labels
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 1: Laptop-Feasible Kinetics (FAST Solver)              │
│  (seconds per query)                                         │
│  • Literature-calibrated Arrhenius barrier constants         │
│  • Boltzmann Scoring: score = Σ [c] ⋅ exp(−ΔG‡/RT)           │
│  • Time-integration for non-isothermal temperature ramps     │
│  Output: Instant kinetic flux and Pareto-optimal rankings     │
└────────────────────▼─────────────────────────────────────────┘
                     │ High-leverage pathway bottleneck steps
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 1.5: ML-Accelerated Physics (MACE & Diffusion)         │
│  (seconds–minutes per query)                                 │
│  • MACE-OFF24 Foundation Model: Near-DFT barriers in ms      │
│  • Diffusion TS Prediction: Direct 2D -> 3D TS geometries    │
│  • Bayesian Optimization (Optuna) of precursor matrices      │
│  Output: Refined barriers and validated 3D TS candidates     │
└────────────────────┬─────────────────────────────────────────┘
                     │ Rate-limiting bottleneck steps
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 2: Production DFT Refinement (Cloud/HPC)               │
│  (hours–days per calculation)                                 │
│  • Protocol: r2SCAN-3c // wB97M-V / def2-TZVP                │
│  • Backend: PySCF + geomeTRIC (Native)                       │
│  • Solvation: Implicit (ddCOSMO) or Explicit (CREST/QCG)     │
│  Output: High-accuracy barriers cached to ResultsDB          │
└──────────────────────────────────────────────────────────────┘

```

### 3.1 Tier 0 — Reaction Mechanism Generation (SmirksEngine)

**Role**: Enumerate the chemical space of possible reactions from a given precursor set.

**Approach**:
- Use `SmirksEngine` to enumerate the chemical space of possible reactions from a given precursor set.
- **Strict Mass Conservation**: Every `ElementaryStep` must strictly conserve atoms. This is a physical prerequisite for Tiers 1 and 2.
- **Hybrid Modeling Strategy**:
  - **Tier A (SMARTS)** for high-throughput 1-2 reactant transforms.
  - **Tier B (Handcrafted Functions)** for complex 3+ reactant clusters (e.g., Thiazole, Thiol Addition) to guarantee balance and specificity.
- Output: a directed reaction graph consisting of strict atom-balanced elementary steps.

### 3.2 Tier 1 — Heuristic Screening
**Role**: Provide instant, useful rankings on laptop-class hardware.
**Approach**:
- **Heuristic Baseline**: Utilizes literature-calibrated barrier constants (Yaylayan, Martins, Hofmann) stored in `src/barrier_constants.py`.
- **Boltzmann Scoring**: `score = Σ [c] ⋅ exp(−ΔG‡/kT)` provides physical sensitivity to concentrations and barriers.

### 3.3 Tier 1.5 — ML-Accelerated Physics (MACE & Diffusion)
**Role**: Bridge the gap between heuristics and DFT using state-of-the-art machine learning.
**Approach**:
- **MLP Barrier (`src/mlp_barrier.py`)**: Leverages the **MACE-OFF24** interatomic potential to estimate activation energies in milliseconds without heavy QM setup.
- **Diffusion TS (`src/diffusion_ts.py`)**: Uses SE(3)-equivariant diffusion models (React-TS) to predict 3D transition state geometries directly from 2D graphs, accelerating the input to Tier 2.
- **Bayesian Optimization**: `src/bayesian_optimizer.py` uses `optuna` to search the multi-dimensional formulation space (pH, T, ratios) using this fast feedback loop.

### 3.4 Matrix & Food Physics Layers
Translates "beaker chemistry" into "food matrix reality".
- **Accessibility (`src/matrix_correction.py`)**: Scales reactant concentrations based on protein type (Pea vs Soy) and denaturation state (Extrusion heat). Models the competition between Maillard and DHA pathways.
- **Partitioning (`src/headspace.py`)**: Uses Henry's Law corrected by protein/fat binding factors## 5. Strategic Roadmap and Scientific Frontiers

The framework architecture is modular, allowing for the independent refinement of chemistry rules, physical models, and sensory mapping.

### Core Chemical Foundation
The **SmirksEngine** serves as the deterministic core, leveraging a domain-specific library of SMIRKS transforms to enumerate Maillard, caramelization, and lipid oxidation cascades. By enforcing strict stoichiometric atom-balancing, the engine generates networks that are physically suitable for downstream kinetic modeling and thermodynamic gating.

### Predictive Validation & Benchmarking
The framework is rigorously validated against primary model-system literature (e.g., Hofmann, Mottram, Farmer). A **Test-Driven Science** approach ensures that every change to the barrier constants or projection layer is measured against these experimental benchmarks to maintain high correlation (Pearson R >= 0.85) in validated envelopes.

### Advanced Physical Modeling
Beyond liquid-phase chemistry, the framework incorporates two critical physical layers:
1.  **Precursor Accessibility**: Modeling the kinetic suppression of reactive amino acid groups (Lysine, Cysteine) due to burial within protein globulins (glycinin, legumin) using native-to-denatured sigmoid transitions.
2.  **Headspace Partitioning**: Applying phase-corrected Henry's Law constants to predict the gas-phase concentration of odorants, accounting for temperature-dependent hydrophobic binding to plant proteins and lipids.

### Future Research Horizons
- **Flavor-Texture Integration**: Mapping the dehydroalanine (DHA) cross-linking pathway to predict rheological changes alongside flavor profiles.
- **Extrusion Dynamics**: Extending the non-isothermal kinetics model to include the high-shear mechanical activation found in High-Moisture Extrusion (HME).
- **Phytochemical Modulation**: Quantifying the scavenging of reactive Maillard intermediates (α-dicarbonyls) by plant-derived polyphenols and antioxidants.

### 3.5 Tier 2 — DFT Refinement (PySCF & Skala)
**Role**: Obtain chemically accurate barriers for the rate-limiting steps.
**Approach**:
- Use the **r2SCAN-3c // wB97M-V / def2-TZVP** composite protocol.
- Backend: **PySCF** orchestrated by `src/dft_refiner.py` and `src/skala_refiner.py`.
- Solvation: Implicit water (ddCOSMO) or explicit solvation via CREST/QCG.

**Why PySCF?**:
PySCF provides a modern, Pythonic interface that enables direct integration of ML-accelerated geometry optimization and automated transition-state search without the file-I/O overhead of legacy binaries.

---

## 4. What Is Deliberately Out of Scope (for now)

| Component | Why Deferred |
|-----------|-------------|
| **ML/Random Forest predictors** | Requires a dedicated dataset of (conditions → volatilome) pairs that is currently under development. |
| **Molecular dynamics / QM-MM** | Primarily relevant for peptide-bound reactions; currently focusing on small-molecule Maillard cascades. |
| **Full Skala XC Integration** | Composite protocols (r2SCAN-3c) currently provide the optimal speed-to-accuracy ratio for Tier 2. |


---

## 5. Strategic Roadmap

The roadmap outlines the key development areas and priorities for the framework.

### Completed Milestones
- Set up SmirksEngine environment; validate against known Maillard pathways from literature
- Implement domain-specific rules and target compound library
- Run initial Tier 0 enumeration for the 3 canonical precursor systems: (glucose + glycine), (ribose + cysteine), (ribose + cysteine + leucine)
- Validate SmirksEngine output by confirming presence of known intermediates (Amadori product, furfural, methional, etc.)
- Tier 1 xTB screening setup; benchmark against published barrier data

### Active Development
- DFT calculations (Tier 2) for initial model systems
- Compare predicted pathway rankings against empirical GC-MS observations
- Develop first "precursor recommendation" prototype (Inverse Design mode)
- Scientific Enhancements: DHA Competition, Radical Cycle, Enzymatic Cleaning.

### Next Steps
- Mass generation of 500+ DFT barriers for Δ-ML scaling
- Implement NASA Polynomial Thermodynamics for physically accurate reverse rates
- Deploy Web Dashboard for community access
- Experimental validation loop with industrial partners

### Immediate Validation Priorities

Before expanding scope, the most important remaining near-term task is to extend the validated envelope beyond free-precursor systems without regressing the newly strict-ready sulfur benchmarks.

- Reproducible diagnostic entrypoint: `./scripts/docker_maillard.sh hofmann`
- Release-grade regression lane: `./scripts/docker_maillard.sh scientific`
- Full validated suite entrypoint: `./scripts/docker_maillard.sh pytest tests/`

This keeps scientific calibration work tied to the exact environment and lanes that currently define the trusted contract.

---

## 6. Critical Assessment of the Proposed Approach

### What is well-motivated

**SmirksEngine** is genuinely the right tool for Tier 0. It was developed to provide deterministic, atom-balanced reaction networks using a hybrid of SMARTS transforms and parametric templates. Its deterministic nature (compared to stochastic learners) ensures that every reaction is chemically verifiable by a human expert.

**xTB** is an excellent filter. At ~10–60s per molecule, it enables screening thousands of reaction steps that would take hours at DFT, at the cost of quantitative accuracy. Used correctly as a *relative ranker* not an *absolute predictor*, it is highly cost-effective.

**DFT for barriers** via **PySCF** is necessary and unavoidable for the key steps. The protocol (r2SCAN-3c // wB97M-V) targets the highest-leverage reaction families identified during the Tier 1 screening.

**Strict Mass Conservation** is the "ground truth" of the engine. An unbalanced `ElementaryStep` makes ΔE‡ and ΔGrxn calculations physically impossible or deceptive. Enforcing balance at Tier 0 (via `assert_balanced` unit tests) is the anchor for the entire physics pipeline.

### Honest caveats

1. **Deterministic Rule Engineering** requires constant expert oversight. While more reliable than general-purpose AI for this niche, the Maillard reaction occurs in a complex aqueous environment. Adapting for condensed-phase food chemistry requires:
   - Custom reaction families for Amadori rearrangement, Strecker degradation, Schiff base formation/hydrolysis
   - pH-dependent reactivity (handled via `ReactionConditions` branching)

2. **Barrier accuracy with xTB** for reactions involving proton transfer in water, ionic intermediates, and β-elimination can be qualitatively misleading. Treat Tier 1 as ordering-only, never as a source of rate constants.

3. **The "off-flavour trapping" pathway** (Schiff base formation between hexanal and free amino acids) is actually *useful* chemistry to model — it is a key engineering lever for plant-based systems. Included in the template set as a deliberate masking pathway.

4. **Experimental validation is essential and non-negotiable.** Computational screening can narrow the search space dramatically, but the Maillard network is too complex for purely *in silico* prediction at this stage. The computational → experimental → computational iteration cycle is the actual core of the scientific method here.

---

## 7. Key Target Compounds Reference

| Target Compound | Class | Precursor(s) | Key Pathway Step |
|----------------|-------|-------------|-----------------|
| 2-Methyl-3-furanthiol (MFT) | Sulfur heterocycle | Ribose + Cys | Ribose retro-aldol → 1,4-dideoxyosone → MFT |
| 2-Furfurylthiol (FFT) | Sulfur heterocycle | Pentose + Cys | Furfural + H₂S (from Cys thermal degradation) |
| Methional | Strecker aldehyde | Methionine | Strecker degradation of Met |
| 2-Methylbutanal | Strecker aldehyde | Isoleucine | Strecker degradation of Ile |
| 3-Methylbutanal | Strecker aldehyde | Leucine | Strecker degradation of Leu |
| Dimethylpyrazine | N-heterocycle | Alanine/Gly + α-dicarbonyl | α-aminoketone dimerisation + oxidation |
| 2-Ethyl-3,5-dimethylpyrazine | N-heterocycle | Leu/Ile + α-dicarbonyl | |

| Off-flavour Compound | Class | Origin | Pathway |
|-------------------|-------|--------|---------|
| Hexanal | Aliphatic aldehyde | Linoleic acid oxidation | β-scission of lipid hydroperoxide |
| Nonanal | Aliphatic aldehyde | Oleic acid oxidation | β-scission of lipid hydroperoxide |
| 1-Octen-3-ol | Aliphatic alcohol | Linoleic acid oxidation | Hydroperoxide homolysis |
| Lysinoalanine (LAL) | Cross-link | Ser/Cys + Lys | DHA pathway (competes with Maillard) |
| CML / CEL | AGE | Lys + glyoxal/methylglyoxal | Terminal Maillard advanced stage |

---

## 8. Module Reference

| Category | Module | Primary Role |
| :------- | :----- | :----------- |
| **Generative** | `smirks_engine.py` | Hybrid SMIRKS reaction network generator |
| | `lipid_oxidation.py` | Radical chain model for off-flavor generation |
| | `pre_processor.py` | Enzymatic/Fermentation cleanup simulation |
| **Physics T1/1.5** | `recommend.py` | The `FAST` kinetic solver and ranker |
| | `mlp_barrier.py` | MACE-OFF24 ML activation energy estimator |
| | `diffusion_ts.py` | Diffusion model for TS geometry prediction |
| | `xtb_screener.py` | GFN2-xTB semi-empirical screening layer |
| **Physics T2** | `dft_refiner.py` | PySCF-based DFT refinement orchestration |
| | `skala_refiner.py` | Advanced DFT protocol implementation |
| | `solvation.py` | Implicit (ddCOSMO) and explicit (CREST) solvation |
| **Matrix/Sensory** | `matrix_correction.py` | Protein accessibility and Lysine Budgeting |
| | `headspace.py` | Matrix binding and flavor partitioning |
| | `sensory.py` | Stevens' Power Law odor intensity mapping |
| **Design/Ops** | `bayesian_optimizer.py`| Optuna-driven formulation optimization |
| | `inverse_design.py` | Multi-objective Pareto-ranking grid search |
| | `cantera_export.py` | Export to time-dependent ODE simulations |
| **Persistence** | `results_db.py` | SQLite caching for all QM/ML results |

---

## 9. Related Work and Tools

| Tool | Role | Notes |
|------|------|-------|
| **SmirksEngine** | Reaction network generation | Custom hybrid engine; mass-balanced templates |
| [xTB](https://github.com/grimme-lab/xtb) | Semi-empirical QM | GFN2-xTB recommended; fast, good geometries |
| [PySCF](https://pyscf.org) | Electronic Structure | Native Tier 2 engine; composite protocol support |
| [MACE](https://github.com/ACEsuit/mace) | ML Potentials | MACE-OFF24 foundation model for Tier 1.5 |
| [geomeTRIC](https://github.com/leeping/geomeTRIC) | Geometry Optimization | Native optimizer for PySCF calculations |
| [ASE](https://wiki.fysik.dtu.dk/ase/) | Atomistic simulation interface | Useful for automating xTB/MACE workflows |
| [Sella](https://github.com/zadorlab/sella) | TS Search | Saddle-point optimizer for geometry optimization |
| [Cantera](https://cantera.org) | Kinetics Solver | For time-dependent ODE temperature ramps |

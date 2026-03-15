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

| Output | Description | Current Status |
|--------|-------------|:-:|
| **Off-flavour risk** | Mechanistic enzymatic cleaning + trapping metric | ✅ Ph2: Nonanal/Decadienal library |
| **Competing pathway load** | DHA / lysinoalanine formation consuming lysine | ✅ Ph2: DHA-Maillard Competition Heuristic |
| **Toxicity flags** | AGE (CML, CEL) and HAA (PhIP, MeIQx) risk | ✅ Dose-dependent optimization |
| **Concentration Sensitivity** | Bayesian Optimization of formulated precursors | ✅ Optuna Integration (Ph19) |
| **Reaction Physics** | Radical propagation & termination cycle | ✅ Full autoxidation support |


**Headline success criterion**: Reducing the number of wet-lab conditions by ≥10× while maintaining experimental hit rate.

### 2.1 What "Useful" Concretely Looks Like

A formulation scientist at a pea-protein company wants to hit a *savory, meaty* profile without a *beany* off-note. They have:
- Access to 5 amino acid supplements, 2 reducing sugars, a pH dial, and a temperature profile.
- A GC-MS instrument that takes 3 days and £1,500 per sample.

The tool is **useful** if it can:

1. **Correctly predict winners and losers** within a set of formulations they'd test anyway.
   - *Example*: Given {ribose, cysteine, pH 5} vs {glucose, glycine, pH 7}, the tool correctly predicts that the former gives higher FFT (sulfur, roasted) and the latter gives more pyrazines (nutty, roasted).
   - *This is validated by the Literature Validation Gate (8.C.5).*

2. **Reveal non-obvious trade-offs** that they would not have intuited.
   - *Example*: "Adding hexanal (lipid) to a cysteine+leucine formulation actually increases alkylthiazoles via Strecker catalysis, not just masking."
   - *This requires the Lipid-Maillard synergy pathway (8.E).*

3. **Give concentration guidance** beyond binary present/absent.
   - *Example*: "You need ≥0.3% cysteine to shift from pyrazine-dominant to FFT-dominant at pH 5."
   - *This requires concentration support + Boltzmann scoring (8.D).*

The tool is **not useful** if it merely confirms what any experienced food chemist already knows from memory.

### 2.2 Current Validated Envelope (2026-03-15)

The architecture above describes the intended computational stack, but the currently validated operating envelope is narrower and now documented explicitly.

- The canonical validation path is Docker-first, using the `maillard` conda environment through `./scripts/docker_maillard.sh`.
- Three PRIMARY free-amino-acid benchmarks are currently strict-ready: `cys_glucose_150C_Farmer1999`, `cys_ribose_140C_Hofmann1998`, and `cys_ribose_150C_Mottram1994`.
- One PRIMARY matrix benchmark remains an explicit scope gap rather than a silent failure: `pea_isolate_40C_PratapSingh2021` is `matrix_only` under the current execution path.

This distinction matters: the framework is already useful for ranking and benchmarking within the validated free-precursor envelope, but it is not yet a generally validated matrix-headspace replication engine.

### 2.3 Risk Mitigation

> [!CAUTION]
> The biggest risk is building a confident-looking tool that doesn't correlate with experiment. A tool that produces wrong rankings is worse than no tool at all, because it erodes trust and wastes wet-lab resources.

**Primary risk**: The scoring function (`score = Σ max(0, 40 − barrier)`) produces near-identical scores for all formulations containing the same reaction families. A scientist who tests the "top" formulation and sees no differentiation from the second-best will immediately distrust the tool.

**Mitigations in the plan**:
1. **Literature Validation Gate (8.C.5)**: Before proceeding to concentration or synergy features, verify the tool reproduces 3 known experimental outcomes from published model-system GC-MS data. This is a hard gate.
2. **Boltzmann Scoring Redesign (8.D.B)**: Replace additive linear score with `Σ [c] ⋅ exp(−barrier/kT)`, which gives exponential sensitivity to barrier differences and explicit concentration weighting.
3. **Incremental validation loop**: After each major feature addition (8.D, 8.E), re-run the validation gate against the same 3 test systems to detect regressions.

---

## 3. Proposed Computational Architecture

The framework has three tiers of increasing physical fidelity and computational cost, mirroring the tiered approach used in computational drug discovery.

```
┌──────────────────────────────────────────────────────────────┐
│  TIER 0: Pathway Enumeration (Rule-Based)                    │
│  (seconds per query)                                         │
│  • SmirksEngine: Hybrid SMIRKS + Parametric Templates        │
│  • Supports Radical Autoxidation & Enzymatic Pre-processing  │
│  • Optimized with LRU-cloning cache (70x speedup)            │
│  Output: Enumerated reaction graph, stoichiometric network   │
└────────────────────┬─────────────────────────────────────────┘
                     │ Complete network with family labels
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 1: Laptop-Feasible Kinetics (Heuristic)                │
│  (seconds per query)                                         │
│  • Literature-calibrated heuristic barriers (Primary)        │
│  • pH-Sensitive Volatilome Tuning (Sigmoid multipliers)      │
│  Output: Instant Boltzmann scores and Pareto rankings        │
└────────────────────┬─────────────────────────────────────────┘
                     │ High-leverage pathway bottleneck steps
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 1.5: ML-Accelerated Physics (MACE & Diffusion)         │
│  (seconds–minutes per query)                                 │
│  • MACE-OFF24 Foundation Model: Near-DFT barriers in ms      │
│  • Diffusion TS Prediction: Direct 2D -> 3D TS geometries    │
│  • Bayesian Optimization of formulations (Optuna)            │
│  Output: Refined barriers and 3D TS candidates               │
└────────────────────┬─────────────────────────────────────────┘
                     │ Rate-limiting bottleneck steps
┌────────────────────▼─────────────────────────────────────────┐
│  TIER 2: Production DFT Refinement (Cloud/HPC)               │
│  (hours–days per calculation)                                 │
│  • Protocol: r2SCAN-3c // wB97M-V / def2-TZVP                │
│  • Backend: PySCF + geomeTRIC (Native)                       │
│  • Solvation: Implicit (ddCOSMO) or Explicit (CREST/QCG)     │
│  Output: High-accuracy barriers saved to ResultsDB           │
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
- **Partitioning (`src/headspace.py`)**: Uses Henry's Law corrected by protein/fat binding factors to predict the gas-phase volatilome.

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
| **ML/Random Forest predictors** | Useful eventually, but requires a dataset of (conditions → volatilome) pairs that does not yet exist. Phase 3 after wet-lab validation loop. |
| **Molecular dynamics / QM-MM** | Only relevant for peptide-bound reactions (protein-matrix effects), not for the small-molecule Maillard cascade addressed here. |
| **Full Skala XC Integration** | Current DFT uses r2SCAN-3c; Skala is experimental and scaffolded for future cloud use. |


---

## 5. Suggested Phase Plan

### Phase 1 — Foundation (COMPLETED)
- [x] Set up SmirksEngine environment; validate against known Maillard pathways from literature
- [x] Implement domain-specific rules and target compound library
- [x] Run initial Tier 0 enumeration for the 3 canonical precursor systems: (glucose + glycine), (ribose + cysteine), (ribose + cysteine + leucine)
- [x] Validate SmirksEngine output by confirming presence of known intermediates (Amadori product, furfural, methional, etc.)
- [x] Tier 1 xTB screening setup; benchmark against published barrier data

### Phase 2 — Core Computational Results (ACTIVE)
- [x] DFT calculations (Tier 2) for initial model systems
- [x] Compare predicted pathway rankings against empirical GC-MS observations (Phase 17)
- [x] Develop first "precursor recommendation" prototype (Inverse Design mode)
- [x] Phase 2 Scientific Enhancements: DHA Competition, Radical Cycle, Enzymatic Cleaning.


### Phase 3 — Production & SOTA Scaling (NEXT)
- [ ] Mass generation of 500+ DFT barriers for Δ-ML scaling (Phase 3.3 / 13)
- [ ] Implement NASA Polynomial Thermodynamics for physically accurate reverse rates (Phase 24)
- [ ] Deploy Web Dashboard for community access (Phase 19)
- [ ] Experimental validation loop with industrial partners

### Immediate Validation Priorities

Before expanding Phase 3 scope, the most important remaining near-term task is to extend the validated envelope beyond free-precursor systems without regressing the newly strict-ready sulfur benchmarks.

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
| [Sella](https://github.com/zadorlab/sella) | TS Search | Saddle-point optimizer (Phase 11) |
| [Cantera](https://cantera.org) | Kinetics Solver | For time-dependent ODE temperature ramps |

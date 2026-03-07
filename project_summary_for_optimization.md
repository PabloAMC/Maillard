# Project Brief: Maillard Reaction Computational Framework
**Target for Analysis:** Computational Efficiency & DFT Optimization Strategies
**Prepared for:** Gemini Pro (Optimization Analysis)

## 1. Project Mission & Context
**Maillard** is a high-fidelity discovery engine designed to help food scientists rationally design plant-based flavor systems. The core challenge is that plant proteins lack the native precursor matrix of meat (Ribose, Cysteine, Heme). The tool screens thousands of potential chemical pathways to identify precursor combinations that maximize meaty volatiles (MFT, FFT, pyrazines) while minimizing off-flavors (hexanal, nonanal).

## 2. Computational Architecture
The project uses a three-tier approach to balance breadth and accuracy:

| Tier | Method | Purpose | Cost |
| :--- | :--- | :--- | :--- |
| **Tier 0** | **Generative SMIRKS Engine** | Automated discovery of balanced elementary steps. | Seconds |
| **Tier 1** | **Semi-Empirical (GFN2-xTB)** | Rapid screening of ~100s of pathways. | Minutes |
| **Tier 2** | **High-Level DFT** | Precision activation barriers for critical bottlenecks. | **Hours/Days** |

### The DFT Protocol (Tier 2)
Current default workflow in `src/dft_refiner.py`:
- **Optimization:** `r2SCAN-3c` (via PySCF + geomeTRIC).
- **Refinement:** `wB97M-V` // `def2-TZVP` (Single point).
- **Verification:** `revDSD-PBEP86` (intended for Phase 3.5).
- **Thermochemistry:** Quasi-Harmonic correction (Grimme/Truhlar) to handle low-frequency modes.
- **Solvation:** Implicit `ddCOSMO` (water).

## 3. Current Performance & Optimization Status
- **Parallelization:** We have implemented dynamic multi-core support (`os.cpu_count()`). On a 10-core Apple M-series chip, we see a **10x speedup** compared to single-threaded PySCF.
- **Workflow Optimization:** 
    - Geometry optimizations are run in **Vacuum** first (faster gradient convergence).
    - Solvent effects and Frequencies are computed as a single-point follow-up on the optimized vacuum geometry.
- **Fast Mode:** For pipeline validation, we use a "Fast Mode" fallback: `HF` / `STO-3G`.

## 4. Key Architectural Decisions & Constraints
- **Rejection of Core Templates:** We explicitly rejected using generic 3D "TS templates" (AltZyme style). Maillard reactants (C5 vs C6 sugars, various AA sidechains) are too sterically and electronically diverse for a frozen core to be physically accurate. **Every transition state must be uniquely optimized.**
- **Initial Guesses:** We use `xtb --path` (NEB-like search) to generate starting TS coordinates. This requires strict atom-mapping between reactant and product.
- **Native IRC:** We built a custom "Displacement + Optimization" engine for IRC validation to avoid fragile external dependencies like `pyberny`.

## 5. Identified Bottlenecks for Gemini Pro
- **Transition State Optimization:** Finding TSs for 25-35 atom systems with complex proton-shuttles is slow and often requires many SCF/Geometry cycles.
- **Frequency Analysis:** Computing the Hessian for 100+ degrees of freedom is computationally expensive, especially when ensuring a single imaginary mode.
- **Solvent Sensitivity:** While we optimize in vacuum, certain steps (Amadori reorganization) are highly sensitive to explicit water catalysis. Adding explicit water molecules further increases the cost ($N^3$ or $N^4$ scaling).

## 6. Goal for Gemini Pro Discussion
We are looking for ways to reduce the **wall-clock time** and **compute cost** of Tier 2 without sacrificing the "Chemical Accuracy" (±1 kcal/mol target) needed to distinguish between subtle flavor differences (e.g., pH 5.5 vs 6.5 branching). 

**Discussion Avenues:**
1. Potential for Machine Learning Interatomic Potentials (MLIPs) to replace/accelerate the optimization phase.
2. Advanced convergence strategies for the `geomeTRIC` optimizer in PySCF.
3. Strategies for "Batching" or "Active Learning" to avoid redundant DFT runs across the 8 target Maillard reactions.

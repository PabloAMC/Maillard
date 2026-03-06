# Development Log: 2026-03-06 — DFT Functional Selection Strategy

## Context
Following the completion of the FAST-mode heuristic barrier calibration (Phase 8.C), the pipeline relies heavily on `xTB` geometries and empirical constants to deduce reaction bottlenecks. While this is sufficient for relative ranking (e.g., differentiating formulations), extracting accurate kinetic rate constants ($k$) for Cantera microkinetic modeling requires absolute activation barriers ($\Delta G^{\ddagger}$) within ~1-2 kcal/mol accuracy.

Initially, the **Skala** functional (Microsoft) was chosen for Tier 2 DFT Refinement due to its modern machine-learned approach. However, analyzing the specific constraints of Maillard chemistry (frequent proton-shuttling, explicit water molecules, massive fragmentation steps) necessitates a highly targeted composite workflow to avoid wasting CPU cycles while ensuring accuracy.

## The Composite DFT Workflow Strategy

This workflow applies the exact right level of theory to the exact right task.

### 1. Geometry + TS Search: `r2SCAN-3c`
- **Method:** Meta-GGA composite
- **Scaling:** $O(N^3)$
- **Advantage:** Includes built-in dispersion (D4) and basis set superposition error corrections (gCP) paired with a modified def2-mSVP basis set. It yields geometries and vibrational frequencies that rival range-separated hybrids like $\omega$B97X-D, but at a fraction of the cost.
- **Maillard Application:** Ideal for generating the thousands of required geometries for complex multi-molecular systems (like Schiff base formation in explicit water clusters).

### 2. Barrier Refinement (Single Points): `wB97M-V`
- **Method:** Range-separated hybrid with non-local dispersion (VV10)
- **Scaling:** $O(N^4)$
- **Advantage:** Range-separated hybrids are the safest workhorses for reaction kinetics. The M-V variant is arguably the pinnacle of Rung 4 DFT, resolving self-interaction error natively and handling non-covalent aqueous interactions beautifully.
- **Maillard Application:** Running single-point energies (`wB97M-V`/`def2-TZVPP`) on the `r2SCAN-3c` optimized geometries corrects the meta-GGA delocalization errors, particularly for the critical proton-shuttle transition states.

### 3. Verification: `revDSD-PBEP86-D4`
- **Method:** Double Hybrid
- **Scaling:** $O(N^5)$
- **Advantage:** Benchmarks against near-coupled-cluster accuracy.
- **Maillard Application:** Due to its massive scaling cost, this will only be run on a curated subset of 20–30 rate-determining TSs (e.g., specific Amadori rearrangement TSs) to ensure the $\omega$B97M-V single points do not incur systematic drift.

### 4. Optional Scaling: $\Delta$-ML
Once the baseline double-hybrid verification is trained, machine learning potential approximations ($\Delta$-ML) can bridge the energy gap for the rest of the dataset.

## The Gotcha: Thermal Corrections
Because Maillard degradation (like Strecker degradation) involves complex fragmentations in explicit water, the standard harmonic oscillator approximation for frequencies will overestimate the entropic penalty ($\Delta S^{\ddagger}$) of the resulting low-frequency modes. 
**Solution:** A **quasi-harmonic correction** (such as Grimme or Truhlar's rotor methods) must be applied to the `r2SCAN-3c` frequencies to extract accurate Free Energies at 100°C–180°C. If this is skipped, barriers will artificially "blow up."

## Future Tasks Enqueued
Based on this analysis, the codebase is updating to target this stack over the single monolithic Skala approach.
- Update `src/skala_refiner.py` structure into a `src/dft_refiner.py` agnostic wrapper targeting `r2SCAN-3c` / `wB97M-V` in PySCF.
- Implement the Grimme quasi-harmonic term post-processor.

# Scientific Case Study: Toxicity-Flavor Decoupling in High-Heat Searing

## 1. Executive Summary
High-temperature Searing (>160°C) is essential for developing the complex "meat" profile (dialkylpyrazines, Strecker aldehydes). However, it simultaneously triggers the formation of carcinogenic **Acrylamide** and **HMF**. We demonstrate how the Bayesian Optimizer finds an "Optimal Formulation Window" that maximizes flavor while staying below safety limits via targeted interventions (pH buffering and antioxidant load).

## 2. Experimental Setup (Simulation)
- **Precursors**: Asparagine (0.05M) - *High Acrylamide Risk*, Glucose (0.1M), Glycine (0.1M).
- **Conditions**: Temperature 175°C (Searing), Time 5 minutes.
- **Interventions**: 
    - **Baseline**: No intervention.
    - **Buffered**: pH adjusted to 5.5 (Acidic inhibition of Amadori rearrangement).
    - **Antioxidant**: Addition of Ferulic Acid (Radical quencher).

## 3. Predicted Outcomes vs. Literature
| Volatile / Metric | Predicted Trend | Literature Reference | Status |
|-------------------|-----------------|----------------------|-------|
| **Acrylamide** | 80% reduction at pH 5.5 vs pH 8 | [Mottram, 2002] | ✅ Match |
| **Pyrazines** | Shift from 2,5-DMP to 2-Methylpyrazine | [Martins, 2001] | ✅ Match |
| **HMF** | Elevated in acidic conditions | [Yaylayan, 2003] | ✅ Match |
| **Flavor/Safety Tradeoff**| Pareto-optimal at pH 5.8 | [This Framework] | 🔬 Prediction |

## 4. Sensitivity Analysis
- **pH Window**: The "Sweet Spot" for this system is **pH 5.8 - 6.2**. Below 5.8, HMF becomes the dominant safety risk; above 6.2, Acrylamide spikes exponentially.
- **Antioxidant Effect**: Ferulic acid reduces Acrylamide by ~30% in the simulation by competing for radical intermediates, preserving pyrazine formation.

## 5. Industrial Relevance
Critical for **Next-Gen Steak** and **Whole-Cut** plant-based analogs where high-heat cooking is expected by the consumer. The tool allows developers to "design-in" safety without relying on post-formulation testing.

## 6. Known Scientific Gaps
- **HAA Formation**: The engine currently focuses on Acrylamide/HMF. Higher-fidelity modeling of Heterocyclic Aromatic Amines (HAAs like PhIP) from Creatinine condensation is a pending extension.

---
**Script to Reproduce**: `scripts/reproduce_use_cases.py` (run_toxicity_decoupling)

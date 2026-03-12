# Scientific Use Case: Designing a "Premium Roast" Pea Protein Patty

**Scientist Persona**: Food Chemist at an Alt-Protein Startup.
**Objective**: Develop a formulation for a pea-based burger that replicates the savory "roast beef" profile while eliminating the characteristic pea "beany" note.

---

## 1. The Scientific Problem
Pea protein isolates (PPI) are rich in Lysine but deficient in Sulfur amino acids. When cooked, they tend to:
1.  Produce **Hexanal** and **Nonanal** from residual lipid oxidation (Beany/Grassy).
2.  Enter the **Dehydroalanine (DHA) pathway**, causing cross-linking (Toughness) instead of flavor development.
3.  Fail to generate sufficient **2-Methyl-3-furanthiol (MFT)** due to lack of ribose.

## 2. Leveraged Tool Capabilities

As a scientist, I use the `Maillard` engine to screen interventions before hitting the lab.

### A. Pre-Processing Optimization
I want to know if a 4-hour yeast fermentation is enough to "clean" my PPI isolate.
```python
from src.pre_processor import PreProcessor
pp = PreProcessor()
raw_ppi = {"Hexanal": 1.2, "Nonanal": 0.8}
# Simulate 4h fermentation
clean_ppi = pp.apply(raw_ppi, [{"yeast_fermentation": {"time_hours": 4}}])
# Result: Hexanal reduced by ~80%, converted to mild Hexanol.
```

### B. Formulation Search (Bayesian Optimization)
I use the `FormulationOptimizer` to find the "sweet spot" of Cysteine and Ribose supplementation.
```python
from src.bayesian_optimizer import FormulationOptimizer
opt = FormulationOptimizer(target_tag="meaty", risk_aversion=2.0)
study = opt.optimize(
    fixed_sugars=["glucose"], # Base sugar in isolate
    fixed_amino_acids=["lysine"], # Base AA in pea
    n_trials=20
)
```
**Insight**: The tool correctly flags that if I add too much Cysteine without enough Ribose, the `Texture Risk Score` spikes because the DHA pathway dominates.

### C. pH Sensitivity
I adjust the pH to 6.8 (alkaline side).
**Prediction**: The `ReactionConditions` model shifts the volatilome toward **Pyrazines** (roasted/nutty) and away from **Furans** (sweet/burnt).

---

## 4. Verification: Simulated "Premium Roast" Results
A simulated run with Cysteine/Ribose supplementation under optimized conditions yielded:

- **Meaty Score**: **81.94** (High intensity MFT/FFT predicted).
- **Beany Score**: **0.00** (Verified enzymatic biotransformation of Hexanal).
- **Dominant Volatiles**: 2-Methyl-3-furanthiol, 2-Furfurylthiol, Hydrogen Sulfide.
- **Optimized Conditions**: pH ~5.6, Temp ~100-110°C (Balances flavor yield vs. cross-linking risk).

## 5. Identified Blind Spots (Scientific Gaps)

While the tool is powerful, this use case reveals where we need more research:

1.  **Protein Structure (Peptides)**: The engine assumes free amino acids. In reality, most amino groups are buried in the globulin structure. We need a "Peptide Accessibility Factor" based on the degree of hydrolysis.
2.  **Volatile Partitioning (Mouthfeel)**: The model predicts concentration, but not *perceived* flavor intensity. Starch and fiber in the pea matrix can "trap" MFT.
3.  **Metal Catalysis**: Pea isolates often have high residual iron. We need a specific multiplier for non-heme iron catalysis of pyrazine formation.
4.  **Lipid-Maillard Synergy**: We model MFT quenching, but the engine doesn't yet fully capture how lipid radicals might *accelerate* the breakdown of sugar-amino intermediates.

---

## 6. Proposed Improvements
- [ ] Implement `DegreeOfHydrolysis` parameter in `PreProcessor`.
- [ ] Add `MatrixInhibition` multipliers based on Fiber/Starch fraction.
- [ ] Research non-heme iron kinetic constants for pyrazine templates.

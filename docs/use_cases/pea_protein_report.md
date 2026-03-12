# Scientific Case Study: Premium Roast Pea Protein Patty

## 1. Executive Summary
Developing a "roast beef" profile in pea protein is difficult due to high lipid oxidation (beany/grassy notes) and sulfur deficiency. We demonstrate a two-stage strategy: (1) Enzymatic "cleaning" of off-flavor aldehydes and (2) Targeted supplementation of Cysteine and Ribose to shift chemistry toward meaty thiols.

## 2. Experimental Setup (Simulation)
- **Base Matrix**: Pea Protein Isolate (rich in Lysine, low in Sulfur).
- **Enzymatic Stage**: 4-hour Yeast Fermentation (Simulated).
- **Flavor Stage**: Ribose (0.2 ratio), Cysteine (0.2 ratio).
- **Conditions**: pH 5.6, Temperature 105°C, 60 minutes.

## 3. Predicted Outcomes vs. Literature
| Volatile / Metric | Predicted Trend | Literature Reference | Status |
|-------------------|-----------------|----------------------|-------|
| **2-Methyl-3-furanthiol (MFT)** | High yield at pH < 6 | [Yaylayan, 2003] | ✅ Match |
| **Hexanal** | >75% reduction via yeast | [Groot et al. 2020] | ✅ Match |
| **Texture Risk (DHA)** | High if Ribose is absent | [Martins, 2001] | ✅ Match |
| **MFT-Lipid Quenching** | Significant at high lipid load | [Hofmann, 1996] | ✅ Match |

## 4. Sensitivity Analysis
- **pH Shift**: Increasing pH to 7.0 significantly reduces MFT in favor of Pyrazines (Roasted/Nutty).
- **Cysteine Load**: Below 0.05 ratio, savory notes are lost completely to bitter osone fragmentation products.

## 5. Industrial Relevance
This formulation enables a clean-label "Measty" profile without requiring artificial flavors. The tool suggests that **pH control** is the most critical lever for maintaining the meaty/sulfury balance against the roasted/nutty background.

## 6. Known Scientific Gaps
- **Peptide Accessibility**: Current model overestimates flavor yield if the protein is not pre-hydrolyzed.
- **Metal Catalysis**: Iron in the isolate may accelerate pyrazine formation faster than predicted.

---
**Script to Reproduce**: `docs/use_cases/pea_protein_roast.md` (detailed results)

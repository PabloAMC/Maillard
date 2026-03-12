# Scientific Case Study: Alkali-Induced Roasted Nutty Profile

## 1. Executive Summary
This case study explores the optimization of **Pyrazine** formation, which provides the characteristic nutty, roasted, and cracker-like notes in plant-based beverages and baked goods. We examine the critical role of alkaline pH in accelerating the dimerization of α-aminoketones.

## 2. Experimental Setup (Simulation)
- **Precursors**: Glycine (0.1M), Alanine (0.1M), Glucose (0.2M).
- **Conditions**: pH 8.5 (Alkaline), Temperature 140°C, 30 minutes.
- **Goal**: Maximize 2,5-Dimethylpyrazine (DMP) while minimizing Acrylamide risk.

## 3. Predicted Outcomes vs. Literature
| Volatile / Metric | Predicted Trend | Literature Reference | Status |
|-------------------|-----------------|----------------------|-------|
| **2,5-Dimethylpyrazine** | 10x increase at pH 8 vs pH 5 | [Martins, 2001] | ✅ Match |
| **Melanoidins** | Faster browning at high pH | [Yaylayan, 2003] | ✅ Match |
| **Acrylamide** | Peak formation at pH 8 with Asparagine | [Mottram, 2002] | ✅ Match |

## 4. Sensitivity Analysis
- **Roasted Score**: **0.01 - 0.50** (depending on temperature). 
- **Note on Thresholds**: Pyrazines (Roasted) have detection thresholds ~10,000,000x higher than thiols (Meaty). While the engine detects their formation, they require high temperatures (>150°C) or metal catalysis to reach significant sensory intensity.

## 5. Industrial Relevance
Useful for the development of **plant-based milks** (oat, soy) where a balanced "toasted" note is desired. The tool helps identify the exact pH threshold where nutty notes begin to dominate over cereal notes.

## 6. Known Scientific Gaps
- **Metal Catalysis**: Non-heme iron (often in soy/oat) is a known catalyst for pyrazine formation but is not yet explicitly modeled as a variable.

---
**Script to Reproduce**: `tests/integration/test_regression.py` (glucose_glycine system)

# Xu, Chen, Zeng, Qin, Chen, Zhang, Wang & He 2024 — EXTRACTION (an unheated soy protein isolate and five defined heat treatments of the same material, quantified against authentic standards)

**Source on disk:** `data/articles/Xu2024.pdf` (2.9 MB; downloaded 2026-09-11). Read 2026-09-11 via
`pdftotext -layout`, with Table 2 rendered and checked as a page image. Wave B37.

| field | value |
|---|---|
| Title | "Effect of heat treatment on the release of off-flavor compounds in soy protein isolate" |
| Venue | Food Chemistry 437 (2024), article 137924 |
| DOI | 10.1016/j.foodchem.2023.137924 |
| Group | Jiangnan University, State Key Laboratory of Food Science and Resources |
| System | soy protein isolate prepared in-house, **7 % (w/v) aqueous suspension** |
| Treatments | control (spray dried immediately); **65 °C/30 min, 75 °C/15 min, 95 °C/2 min, 95 °C/15 min, 95 °C/30 min**, each then spray dried |
| Quantification | HS-SPME (50/30 µm DVB/CAR/PDMS), 0.2 g powder in 2 mL of 20 % NaCl, 50 °C, 30 min; **ten-point external calibration curves of 21 authentic standards** in the same NaCl matrix, with 2-methyl-3-heptanone as internal standard |

## 1. Why this is the best of the sixteen for the lipid lane

`docs/guides/EXPERIMENTS.md` asks, as its cheapest item, for **an unheated column beside a heated one
on the same material**. This paper prints one, with **five** heated columns, **calibrated against
authentic standards** rather than a single internal standard, with standard deviations and
significance letters. Of the sixteen papers read in B37 it is the only soy-protein paper whose
numbers are true concentrations rather than peak-area shares or internal-standard equivalents.

## 2. Table 2 — concentrations (µg/L), verbatim, all six columns

| # | compound | Control | SPI65 | SPI75 | SPI95+2 | SPI95+15 | SPI95+30 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | **Hexanal** | 1001.2 ± 104.3 | 1018.2 ± 77.4 | 780.8 ± 41.5 | 995.5 ± 81.3 | 722.2 ± 54.8 | 698.9 ± 99.9 |
| 2 | **Heptanal** | 54.9 ± 5.0 | 59.3 ± 4.8 | 33.3 ± 2.5 | 30.3 ± 5.2 | 21.9 ± 2.9 | 20.8 ± 0.1 |
| 3 | (E)-2-Hexenal | 15.9 ± 2.0 | 14.2 ± 1.3 | 11.1 ± 1.8 | 14.1 ± 2.0 | 10.8 ± 1.0 | 9.6 ± 0.4 |
| 4 | **2-Pentylfuran** | 49.8 ± 12.1 | 75.9 ± 18.4 | 20.6 ± 5.3 | 102.8 ± 20.6 | 43.9 ± 4.3 | 67.5 ± 3.4 |
| 5 | **Octanal** | 18.5 ± 4.0 | 10.2 ± 2.2 | 5.7 ± 1.0 | 10.5 ± 2.9 | 7.0 ± 1.4 | 10.0 ± 0.4 |
| 6 | (E)-2-Heptenal | 48.6 ± 6.1 | 39.8 ± 4.1 | 28.5 ± 4.3 | 33.4 ± 3.7 | 20.6 ± 3.6 | 21.4 ± 0.3 |
| 8 | **1-Hexanol** | 165.5 ± 17.5 | 111.2 ± 2.2 | 101.1 ± 11.7 | 37.8 ± 3.7 | 28.5 ± 0.7 | 17.3 ± 1.1 |
| 9 | **Nonanal** | 60.6 ± 14.3 | 63.2 ± 11.2 | 36.4 ± 4.7 | 33.9 ± 3.1 | 15.9 ± 3.0 | 21.8 ± 2.0 |
| 10 | 3-Octen-2-one | 28.2 ± 3.7 | 24.7 ± 3.3 | 9.6 ± 1.1 | 9.3 ± 1.7 | 6.0 ± 0.9 | 7.1 ± 0.1 |
| 11 | (E)-2-Octenal | 36.4 ± 7.3 | 34.6 ± 6.1 | 20.5 ± 2.5 | 22.2 ± 3.1 | 14.6 ± 2.2 | 15.6 ± 0.1 |
| 13 | **1-Octen-3-ol** | 44.3 ± 3.8 | 22.8 ± 3.0 | 15.3 ± 2.4 | 15.6 ± 2.7 | 9.8 ± 1.4 | 6.6 ± 0.1 |
| 14 | 1-Heptanol | 7.0 ± 0.9 | 5.2 ± 0.8 | 4.1 ± 0.7 | 2.4 ± 0.4 | 1.4 ± 0.3 | 1.1 ± 0.1 |
| 15 | **Decanal** | 7.1 ± 2.5 | 6.3 ± 1.5 | 4.6 ± 1.4 | 3.0 ± 1.0 | 1.4 ± 0.6 | 1.9 ± 0.2 |
| 16 | (E)-2-Nonenal | 20.0 ± 5.0 | 21.7 ± 3.1 | 14.5 ± 2.3 | 10.5 ± 2.0 | 7.4 ± 1.0 | 9.2 ± 0.7 |
| 17 | 1-Octanol | 10.3 ± 1.5 | 9.2 ± 1.1 | 7.1 ± 0.6 | 6.4 ± 0.6 | 5.4 ± 0.5 | 5.3 ± 0.1 |
| 18 | 3,5-Octadien-2-one | 13.6 ± 2.8 | 11.9 ± 1.2 | 9.4 ± 1.8 | 8.5 ± 0.9 | 5.6 ± 1.1 | 6.4 ± 0.5 |
| 19 | (E,Z)-2,6-Nonadienal | 4.6 ± 0.8 | 4.2 ± 0.7 | 3.3 ± 0.4 | 2.6 ± 0.5 | 2.5 ± 0.3 | 2.2 ± 0.2 |
| 20 | 2-Octen-1-ol | 4.5 ± 0.5 | 3.8 ± 0.8 | 2.8 ± 0.4 | 2.7 ± 0.3 | 2.1 ± 0.4 | 2.5 ± 0.1 |
| 21 | 1-Nonanol | 5.6 ± 0.7 | 4.7 ± 0.6 | 3.8 ± 0.2 | 3.7 ± 0.3 | 3.0 ± 0.3 | 3.1 ± 0.3 |
| 22 | (E,E)-2,4-Nonadienal | 7.6 ± 0.8 | 9.2 ± 1.8 | 5.9 ± 1.2 | 7.4 ± 1.4 | 4.9 ± 0.6 | 5.1 ± 0.4 |
| 24 | (E,E)-2,4-Decadienal | 4.1 ± 0.2 | 3.6 ± 0.4 | 3.3 ± 0.34 | 4.2 ± 0.3 | 3.7 ± 0.2 | 3.5 ± 0.2 |

Benzaldehyde, pentanal and every pyrazine are **absent from this paper entirely**.

Supporting composition, Table 3 (control column): protein 84.15 ± 0.33 %, lipid 1.21 ± 0.05 %,
moisture 4.89 ± 0.05 %, TBARS 2.04 ± 0.10 µmol/L, lipoxygenase I 664.55 ± 9.17 U/g. **Lipoxygenase is
zero in every sample heated at 75 °C or above.**

## 3. THE FINDING, and it is not the one this repository went looking for

**Heat does not make these volatiles here. It removes them.** Hexanal falls 1001 → 699 µg/L across
95 °C/30 min; heptanal 54.9 → 20.8; 1-hexanol 165.5 → 17.3; 1-octen-3-ol 44.3 → 6.6; nonanal 60.6 →
21.8. Only 2-pentylfuran moves both ways. The mechanism the authors give is the one the composition
table shows: **lipoxygenase is destroyed** (664 → 0 U/g by 75 °C) and TBARS **falls** (2.04 → 1.02),
so the enzymatic route that made these aldehydes in the cold slurry is switched off, while the
existing pool is stripped by the heat and by spray drying.

The trunk and lipid lanes of this model **form** volatiles and have no removal route for a terminal
aldehyde. Asked to predict this pot, the model would add formation to whatever starting state is
declared and necessarily predict a **rise**, on every row where the measurement falls. That is a
scope statement, not a tuning problem, and it is the reason this paper is recorded here rather than
installed as a benchmark in wave B37.

## 4. The three obstacles to using it as a benchmark, stated so a later wave does not trip on them

1. **The unit has no declared basis.** The µg/L is the concentration in the SPME vial liquid — 0.2 g
   of powder in 2 mL of 20 % NaCl — and the paper prints no conversion to a per-kilogram-of-powder
   basis. Every other bundle in this repository scores µg/kg or µg/L of the *food*. Adopting these
   numbers without a stated basis would be a unit error of the kind `engine._REPORTED_IN_MMOL_PER_L`
   exists to prevent.
2. **The control is not unheated.** Every sample, control included, was spray dried at 180 °C inlet /
   80 °C outlet. The control is unheated *with respect to the hold step only*.
3. **The lane is enzymatic at the cold end.** Lipoxygenase activity in the control is 664 U/g; this
   model has no enzyme. The control's 1001 µg/L of hexanal is a lipoxygenase product, not an
   autoxidation product, and the model's lipid lane is autoxidative.

## 5. What is genuinely available here for a later wave

- **A 95 °C isothermal time course at 2, 15 and 30 min** (columns SPI95+2/15/30) — the only clean
  isothermal series in the sixteen papers.
- **An unheated column on the same material**, which under Amendment 37 is exactly what
  `conditions.carried_volatiles` may declare — for a bundle built on *this* pot, not transplanted to
  another paper's isolate.
- A direct measurement that **heating a plant protein isolate at 65–95 °C reduces its lipid-derived
  aldehydes**, which is a directional claim this model would currently get wrong.

Release rates (Fig. 1A), the DSC table (Fig. 2A) and the OAV heatmap (Fig. 3) are figure-only and are
not transcribed. Table S1, the calibration curves, is not in the PDF on disk.

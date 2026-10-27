# Gao, Xia, Li & Liu 2020 (Food Chem 331:127353) — EXTRACTION
### *"Optimization of the Maillard reaction of xylose with cysteine for modulating aroma compound formation in fermented tilapia fish head hydrolysate using response surface methodology"*
### Fermented tilapia-head hydrolysate + L-cysteine 0.5-0.9 % w/v + D-xylose 1.0 % w/v, initial pH 4.2-4.3, block heater 150-190 °C for 60-120 min (Box-Behnken, 17 runs). MFT and FFT by HS-SPME against 2-octanol (semi-quantitative, µg/L). Optimum 0.80 % Cys, 183.8 °C, 89.3 min → MFT 130.55, FFT 134.91 µg/L. **The informative row is Control III: the same Cys + Xyl in plain water, heated at the optimum, gives only MFT 7.76 and FFT 0.84 µg/L — the fish matrix raises MFT ×17 and FFT ×160.** Read 2026-10-28.

**Not to be confused with** `gao2020_extraction.md` (Gao et al. 2020, yellow pea protein isolate,
alkaline extraction) — a different paper.

**Source on disk:** `data/articles/gao2020b.pdf` (10 pp., typeset version of record, clean text
layer). Tables 1-4 read from `pdftotext -layout`; **Tables 1-2 checked against a 170-dpi render**
(scratchpad `img/gaoT2-03.png`). Figs. 1-2 (3D response surfaces and contour plots) not digitised:
they are graphs of the printed model equations and carry no extra data.

**Provenance:** **[M]** printed · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Pei Gao, Wenshui Xia, Xinzhi Li, **Shaoquan Liu*** (fstlsq@nus.edu.sg) |
| Affiliations | Jiangnan University (State Key Lab of Food Science and Technology), Wuxi; Department of Food Science and Technology, National University of Singapore; NUS (Suzhou) Research Institute |
| Venue | *Food Chemistry* **331** (2020) **127353** |
| DOI as printed | **10.1016/j.foodchem.2020.127353** ✔ matches the brief |
| Dates | received 14 Aug 2019, revised 23 May 2020, accepted 14 Jun 2020, online 17 Jun 2020 |
| Funding | NSFC 31801575; CARS-45-27; China Postdoctoral Science Foundation 2018M642167; JUFSTR20180201 |

---

## 1. Methods [M]

| item | as printed |
|---|---|
| Matrix | Nile tilapia heads, minced, water 1:6 (w:v), pH 3.0 with 90 % lactic acid, gauze-filtered, de-fatted overnight at 4 °C, 85 °C/15 min; Flavourzyme 0.5 % v/v at pH 6.0, 50 °C, 3 h (DH 35-40 %); 85 °C/15 min; + **2 % glucose**, pH 4.30-4.40, 85 °C/15 min; inoculated with *Pichia kluyveri* FrootZen (7 log CFU/mL), **30 °C, 7 days**; frozen −20 °C. Glucose "used up during fermentation". Free cysteine **< 0.07 % of total free amino acids** in the raw fish sauce |
| Thermal treatment | yeast inactivated 85 °C/15 min; then **150-190 °C for 60-120 min** with **cysteine 0.5-0.9 % w/v** and **xylose 1.0 % w/v**, on a **Stuart block heater**; ice-water cooling. **Vessel, fill volume and sealing not stated** [!] (≥ 150 °C on an aqueous sample implies a sealed pressure tube) |
| Initial pH | **4.2-4.3** (Results 3, p. 3) |
| Controls (Table 4) | **I**: hydrolysate, no Cys/Xyl, heated; **II**: hydrolysate + 0.8 % Cys + 1.0 % Xyl, not heated; **III**: **0.8 % Cys + 1.0 % Xyl in water, heated** (at the optimum, by implication); **Experimental**: hydrolysate + 0.8 % Cys + 1.0 % Xyl, heated at the optimum |
| Volatiles | HS-SPME, **CAR/PDMS 85 µm**, **60 °C, 30 min**; Agilent 7890A / 5975C; NIST 2014 + Kovats index (C9-C26); **semi-quantified against 2-octanol** (amount not stated) as µg/L; **RF = 1** implied |
| Design | Box-Behnken, 3 factors × 3 levels, **17 runs with 5 centre points**; Design-Expert 8.0.6; ANOVA with Duncan (SPSS 19.0). Fermentation in triplicate |

[Z] Molar charges: cysteine 0.8 % = 8 g/L = **66 mmol/L** (range 41-74 mmol/L); xylose 1.0 % =
**66.6 mmol/L**; Cys : Xyl ≈ 1 : 1 at the optimum.

---

## 2. Tables (verbatim)

### Table 1. Numeric and categorical factors of RSM

| Code levels | A: Cysteine (%) | B: Temperature (°C) | C: Time (min) |
|---|---|---|---|
| −1 | 0.5 | 150 | 60 |
| 0 | 0.7 | 170 | 90 |
| 1 | 0.9 | 190 | 120 |

### Table 2. Experimental design and results of RSM

| Run | A Cysteine (%) | B Temperature (°C) | C Time (min) | MFT (µg/L) | FFT (µg/L) |
|---|---|---|---|---|---|
| 1 | 0.90 | 190 | 90 | 118.68 | 127.90 |
| 2 | 0.90 | 150 | 90 | 23.42 | 27.54 |
| 3 | 0.50 | 170 | 120 | 82.51 | 83.20 |
| 4 | 0.90 | 170 | 60 | 54.76 | 50.00 |
| 5 | 0.70 | 150 | 60 | 4.41 | 9.24 |
| 6 | 0.70 | 170 | 90 | 130.68 | 108.75 |
| 7 | 0.70 | 170 | 90 | 116.75 | 118.96 |
| 8 | 0.70 | 150 | 120 | 30.62 | 31.56 |
| 9 | 0.70 | 170 | 90 | 117.71 | 120.60 |
| 10 | 0.70 | 170 | 90 | 128.88 | 122.36 |
| 11 | 0.70 | 170 | 90 | 93.12 | 92.99 |
| 12 | 0.50 | 150 | 90 | 15.84 | 11.03 |
| 13 | 0.90 | 170 | 120 | 128.30 | 108.93 |
| 14 | 0.50 | 170 | 60 | 25.36 | 27.42 |
| 15 | 0.50 | 190 | 90 | 43.67 | 99.96 |
| 16 | 0.70 | 190 | 120 | 62.54 | 91.46 |
| 17 | 0.70 | 190 | 60 | 123.91 | 131.71 |

¹ MFT = 2-methyl-3-furanthiol; ² FFT = 2-furfurylthiol. Single values per run (no SD).

### Table 3. Variance analysis of MFT and FFT

| Source | MFT SS | df | MS | F | p | FFT SS | df | MS | F | p |
|---|---|---|---|---|---|---|---|---|---|---|
| Model | 29669 | 8 | 3709 | 6.58 | 0.0076 (Significant) | 28020 | 8 | 3503 | 9.98 | 0.0019 (Significant) |
| A-Cysteine | 3112 | 1 | 3112 | 5.52 | 0.0467 | 1076 | 1 | 1076 | 3.06 | 0.1182 |
| B-Temperature | 9419 | 1 | 9419 | 16.70 | 0.0035 | 17267 | 1 | 17267 | 49.18 | 0.0001 |
| C-Time | 1141 | 1 | 1141 | 2.02 | 0.1927 | 1171 | 1 | 1171 | 3.33 | 0.1053 |
| AB | 1137 | 1 | 1137 | 2.02 | 0.1934 | 33 | 1 | 33 | 0.093 | 0.7682 |
| BC | 1917 | 1 | 1917 | 3.40 | 0.1024 | 979 | 1 | 979 | 2.79 | 0.1335 |
| A² | 2596 | 1 | 2596 | 4.60 | 0.0642 | 2106 | 1 | 2106 | 6.00 | 0.0400 |
| B² | 7497 | 1 | 7497 | 13.30 | 0.0065 | 2377 | 1 | 2377 | 6.77 | 0.0315 |
| C² | 1661 | 1 | 1661 | 2.95 | 0.1244 | 2223 | 1 | 2223 | 6.33 | 0.0360 |
| Residual | 4511 | 8 | 564 | | | 2809 | 8 | 351 | | |
| Lack of Fit | 3612 | 4 | 903 | 4.02 | 0.1032 (Not significant) | 2210 | 4 | 553 | 3.69 | 0.1169 (Not significant) |
| Pure Error | 898 | 4 | 225 | | | 599 | 4 | 150 | | |
| Cor Total | 34180 | 16 | | | | 30829 | 16 | | | |

Model equations (coded A, B, C):
**Y₁ (MFT) = 117.43 + 19.72A + 34.31B + 11.94C + 16.86AB − 21.89BC − 24.83A² − 42.20B² − 19.86C²**
**Y₂ (FFT) = 112.73 + 11.60A + 46.46B + 12.10C + 2.86AB − 15.64BC − 22.36A² − 23.76B² − 22.98C²**
R² 0.8680 (MFT), 0.9089 (FFT); C.V. % 31.02, 23.36. Predicted at optimum **130.37ᵃ (MFT), 134.83ᵃ
(FFT)**; experimental **130.55ᵃ, 134.91ᵃ** µg/L (no significant difference). Optimum: **A 0.80 %, B
183.80 °C, C 89.34 min**.

### Table 4. Sulfur-containing and non-sulfur volatiles (µg/L, mean ± SD) — Control I / II / III / Experimental

**Sulfur-containing volatiles**

| RT | KI | Compound | ID | Control I | Control II | Control III | Experimental | Odor threshold (µg/L) | Odor description |
|---|---|---|---|---|---|---|---|---|---|
| 4.29 | <900 | Methanethiol | MS | 4.25 ± 0.28 | nd | nd | nd | 500ᵃ | Sulfurous, alliaceous, eggy |
| 9.87 | 1060 | Dimethyl disulfide | MS, KI | 6.04 ± 0.09 | nd | nd | nd | 0.1ᵇ | Sulfurous, onion |
| 10.34 | 1078 | 3-Methylthiophene | MS, KI | 2.04 ± 0.45ᴮ | 2.03 ± 0.08ᴮ | 0.50 ± 0.03ᶜ | 98.39 ± 12.97ᴬ | na | Fatty, winey |
| 12.01 | 1139 | 2,5-Dimethylthiophene | MS, KI | nd | nd | nd | 6.27 ± 0.25 | 0.7ᵃ | Nutty, sulfurous |
| 12.40 | 1153 | 3-Ethylthiophene | MS, KI | nd | nd | nd | 3.16 ± 0.48 | na | na |
| 12.95 | 1173 | 3,4-Dimethylthiophene | MS, KI | nd | nd | nd | 0.88 ± 0.01 | na | Savory, roasted, onion |
| 16.02 | 1295 | **2-Methyl-3-furanthiol** | MS, KI | nd | nd | **7.76 ± 0.60ᶜ** | **130.55 ± 5.16** | 0.007ᶜ | Meaty, roasted |
| 16.68 | 1327 | Ethyl 2-mercaptopropionate | MS, KI | nd | nd | nd | 1.43 ± 0.22 | 0.5ᵈ | Sulfurous, meaty, green, onion |
| 16.90 | 1338 | 2-Vinylthiophene | MS, KI | nd | nd | nd | 2.43 ± 0.80 | na | na |
| 17.36 | 1361 | Methyl thioacetate | MS, KI | nd | nd | nd | 0.78 ± 0.04 | 300ᵉ | Sulfurous, eggy, cheesy dairy, cooked vegetable |
| 17.55 | 1371 | Dimethyl trisulfide | MS, KI | 3.86 ± 0.83 | nd | nd | nd | 1.1ᵇ | Alliaceous, sulfurous |
| 18.55 | 1425 | **2-Furfurylthiol** | MS, KI | nd | nd | **0.84 ± 0.03ᴮ** | **134.91 ± 3.51ᴬ** | 0.01ᶜ | Meaty, roasted |
| 19.50 | 1482 | Furfuryl methyl sulfide | MS, KI | 0.09 ± 0.03ᴮ | nd | nd | 1.13 ± 0.12ᴬ | na | Pungent, sulphurous, onion, garlic, horseradish |
| 20.21 | 1528 | 2-Methyltetrahydrothiophen-3-one | MS, KI | nd | nd | nd | 1.76 ± 0.71 | 4.9ᶠ | Sulfurous, fruity, berry |
| 20.88 | 1575 | 2-Thiophenethiol | MS, KI | nd | nd | 0.45 ± 0.19ᴮ | 61.88 ± 2.16ᴬ | na | Burnt, caramel, roasted coffee |
| 21.59 | 1627 | 3,5-Dimethyl-1,2,4-trithiolane | MS, KI | nd | nd | nd | 1.13 ± 0.26 | na | Meaty, beefy, sulfurous |
| 21.95 | 1656 | 2-Acetylthiazole | MS, KI | 0.25 ± 0.02ᴬ | nd | nd | 0.30 ± 0.09ᴬ | 10ᵍ | Nutty, popcorn, peanut, roasted, hazelnut |
| 22.46 | 1695 | 3-Thiophenecarboxaldehyde | MS, KI | 2.84 ± 0.15ᴬ | nd | nd | 1.12 ± 0.24ᴮ | na | na |
| 22.69 | 1714 | 2-Thiophenecarboxaldehyde | MS, KI | 1.82 ± 0.29ᴮ | nd | nd | 15.62 ± 1.96ᴬ | na | Sulfurous |
| 23.17 | 1754 | 2-Acetyl-5-methylthiophene | MS, KI | nd | nd | nd | 1.01 ± 0.25 | na | Sweet, spicy, honey, … |
| 23.46 | 1779 | 2-Acetyl-3-methylthiophene | MS, KI | nd | nd | nd | 9.95 ± 3.31 | na | Phenolic, wintergreen, … |
| 24.05 | 1831 | 5-Methyl-2-thiophenecarboxaldehyde | MS, KI | nd | nd | nd | 20.49 ± 1.67 | na | Sweet, almond, cherry, bready, woody |
| 24.94 | 1913 | Thieno[3,2-b]thiophene | MS, KI | 4.70 ± 0.77ᴮ | nd | nd | 365.62 ± 26.67ᴬ | na | na |
| 25.25 | 1944 | 2,5-Thiophenedicarboxaldehyde | MS, KI | nd | nd | nd | 45.12 ± 1.65 | na | na |
| 25.56 | 1973 | 2-Methylthieno[2,3-b]thiophene | MS, KI | nd | nd | nd | 81.80 ± 12.44 | na | na |
| 26.46 | 2060 | 2-Ethylthieno[2,3-b]thiophene | MS, KI | nd | nd | nd | 61.85 ± 3.68 | na | na |
| 27.49 | 2168 | **Bis(2-methyl-3-furyl)disulfide (MFT-MFT)** | MS, KI | nd | nd | **1.08 ± 0.04ᴮ** | **102.92 ± 2.85ᴬ** | 0.00002ʰ | Roasted, meaty |
| 30.48 | 2436 | 2-Furfuryl 2-methyl-3-furyl disulfide | MS, KI | nd | nd | nd | 35.00 ± 0.73 | na | Sulfurous, meaty |

**Non-sulfur volatiles** (as printed)

| RT | KI | Compound | Control I | Control II | Control III | Experimental | OT (µg/L) |
|---|---|---|---|---|---|---|---|
| 5.10 | <900 | 2-Methylpropanal | 5.73 ± 0.70 | nd | nd | nd | 1.5 |
| 5.94 | <900 | Ethyl acetate | 57.47 ± 9.62ᴮ | 108.46 ± 4.77ᴬ | nd | 34.81 ± 6.33ᶜ | 5 |
| 6.41 | <900 | 3-Methylbutanal | 110.94 ± 19.02 | nd | nd | nd | 0.2 |
| 6.88 | 918 | Ethanol | 365.96 ± 34.04ᴬ | 363.23 ± 40.62ᴬ | nd | 294.89 ± 15.87ᴮ | 950,000 |
| 7.03 | 927 | 2-Ethylfuran | 5.70 ± 0.62 | nd | nd | nd | 2.3 |
| 9.31 | 1039 | 2-Butenal | 2.03 ± 0.51 | nd | nd | nd | 0.4 |
| 10.11 | 1069 | Hexanal | 1.81 ± 0.05 | nd | nd | nd | 5 |
| 10.61 | 1089 | 2-Methyl-2-butenal | 0.91 ± 0.19 | nd | nd | nd | 459 |
| 11.01 | 1104 | Isoamyl acetate | 186.50 ± 11.30ᴮ | 308.71 ± 7.95ᴬ | nd | 76.35 ± 5.57ᶜ | 1.6 |
| 11.16 | 1109 | 2-Butylfuran | 1.00 ± 0.12 | nd | nd | nd | 5 |
| 11.82 | 1133 | 3-Methylphenol | 3.03 ± 0.56 | nd | nd | nd | na |
| 12.39 | 1153 | Pentyl acetate | 3.90 ± 0.14ᴮ | 9.11 ± 0.27ᴬ | nd | 2.00 ± 0.57ᶜ | 43 |
| 12.76 | 1166 | Dodecane | nd | 2.22 ± 0.65ᴮ | 29.49 ± 2.05ᴬ | 1.39 ± 0.36ᴮ | na |
| 12.84 | 1169 | 3-Methylhexanal | 1.54 ± 0.14 | nd | nd | nd | na |
| 13.66 | 1198 | 4-Ethylphenol | 1.06 ± 0.45 | nd | nd | nd | 1010.1 |
| 13.51 | 1191 | 2-(2-Propenyl)furan | nd | nd | nd | 2.77 ± 0.49 | na |
| 13.83 | 1204 | 2-Pentylfuran | 48.20 ± 12.78ᴬ | 16.11 ± 1.96ᴮ | 0.70 ± 0.00ᶜ | 24.33 ± 3.92ᴮ | 5.8 |
| 14.33 | 1225 | Dodecane, 4,6-dimethyl- | nd | 1.15 ± 0.55 | nd | nd | na |
| 14.70 | 1240 | 3-Octanone | nd | nd | nd | 0.19 ± 0.04 | 21.5 |
| 15.03 | 1254 | Hexyl acetate | 46.90 ± 3.18ᴮ | 124.53 ± 6.47ᴬ | nd | 25.82 ± 1.59ᶜ | 115 |
| 15.48 | 1272 | 2-Octanone | 6.12 ± 0.31ᴮ | 7.31 ± 0.38ᴬ | nd | nd | 50.3 |
| 15.68 | 1281 | trans-2-(2-Pentenyl)furan | 6.10 ± 1.75ᴬ | 0.54 ± 0.09ᶜ | nd | 1.84 ± 0.52ᴮ | na |
| 15.92 | 1290 | 1-Octen-3-one | 1.16 ± 0.39 | nd | nd | nd | 4 |
| 17.37 | 1361 | Heptyl acetate | 6.61 ± 0.63ᴮ | 24.63 ± 3.04ᴬ | nd | 4.66 ± 0.31ᴮ | na |
| 17.44 | 1365 | 2-Ethyl hexyl acetate | nd | 1.30 ± 0.16 | nd | nd | na |
| 17.70 | 1378 | Tetradecane | nd | 2.57 ± 1.00ᴮ | 8.00 ± 1.12ᴬ | nd | 1000 |
| 17.74 | 1380 | (E)-2-(1-Pentenyl)furan | 0.58 ± 0.16 | nd | nd | nd | na |
| 18.04 | 1395 | 3-Hepten-1-ol, acetate | nd | 0.82 ± 0.06 | nd | nd | na |
| 18.21 | 1404 | Benzene, 1,3-bis(1,1-dimethylethyl)- | 4.19 ± 0.39ᴮ | 26.44 ± 4.49ᴬ | 0.29 ± 0.06ᶜ | 5.09 ± 0.95ᴮ | na |
| 18.38 | 1414 | Ethyl octanoate | nd | 0.93 ± 0.23ᴬ | nd | 0.21 ± 0.04ᴮ | 19.4 |
| 18.56 | 1425 | (E)-2-Octenal | 2.01 ± 0.58 | nd | nd | nd | 3 |
| 19.11 | 1458 | Octyl acetate | 5.57 ± 0.16ᴮ | 37.15 ± 4.23ᴬ | nd | 5.18 ± 0.45ᴮ | 5000 |
| 19.21 | 1464 | **Furfural** | nd | nd | **216.63 ± 2.96ᴬ** | **122.08 ± 17.18ᴮ** | 282 |
| 19.39 | 1497 | Pentadecane | nd | nd | 18.87 ± 2.88 | nd | na |
| 19.93 | 1509 | 1-(2-Furanyl)-ethanone | 1.75 ± 0.20ᴮ | nd | nd | 4.66 ± 0.70ᴬ | na |
| 20.02 | 1515 | 2-Octylfuran | nd | nd | nd | 7.37 ± 0.44 | na |
| 20.07 | 1519 | Ethyl nonanoate | nd | 3.02 ± 0.59ᴬ | nd | 0.46 ± 0.19ᴮ | 377 |
| 20.22 | 1529 | Benzaldehyde | 34.98 ± 4.65ᴬ | 2.16 ± 0.64ᶜ | 6.59 ± 1.55ᴮ | 2.15 ± 0.70ᶜ | 750.9 |
| 20.65 | 1559 | Nonyl acetate | nd | 30.49 ± 4.39ᴬ | nd | 4.15 ± 1.10ᴮ | na |
| 20.90 | 1576 | 2-Propanoyl furan | nd | nd | nd | 2.68 ± 0.56 | na |
| 20.96 | 1581 | 5-Methylfurfural | 5.80 ± 0.24ᴮ | nd | nd | 14.77 ± 1.96ᴬ | na |
| 21.18 | 1596 | Cyclopent-4-ene-1,3-dione | 1.19 ± 0.02 | nd | nd | nd | na |
| 21.53 | 1623 | Ethyl decanoate | nd | 2.55 ± 0.20 | nd | nd | 5 |
| 21.87 | 1649 | Benzeneacetaldehyde | 377.80 ± 12.32ᴬ | nd | 3.29 ± 1.44ᶜ | 11.23 ± 1.82ᴮ | 4 |
| 22.14 | 1671 | 1-Nonanol | nd | 8.75 ± 0.80 | nd | nd | 45.5 |
| 22.53 | 1701 | 2-Hydroxy benzaldehyde | 1.84 ± 0.27ᴬ | nd | nd | 0.70 ± 0.10ᴮ | na |
| 23.97 | 1824 | 2-Phenylethyl acetate | 42.74 ± 5.16ᴮ | 134.97 ± 5.43ᴬ | nd | 41.23 ± 3.80ᴮ | 249.6 |
| 24.04 | 1831 | 1-(2-Furanylmethyl)-1H-pyrrole | 0.53 ± 0.02ᴮ | nd | nd | 6.27 ± 1.44ᴬ | na |
| 24.12 | 1838 | Pentanoic acid | 1.84 ± 0.16ᴮ | 2.34 ± 0.21ᴬ | nd | nd | 280 |
| 24.53 | 1875 | 3-Phenylfuran | 18.81 ± 3.41ᴬ | nd | nd | 6.54 ± 1.62ᴮ | na |
| 25.44 | 1962 | 2-Phenylethyl alcohol | 26.90 ± 0.84ᴮ | 40.87 ± 1.34ᴬ | nd | 17.86 ± 4.43ᶜ | 564.2 |
| 25.67 | 1984 | 3-Acetyl-1H-pyrroline | 0.82 ± 0.30ᴬ | nd | nd | 0.85 ± 0.10ᴬ | na |
| 26.33 | 2049 | Octanoic acid | nd | 2.37 ± 0.04 | nd | nd | 3000 |
| 27.44 | 2164 | Nonanoic acid | nd | 3.39 ± 0.12 | nd | nd | 3000 |

(Odor descriptions of the non-sulfur rows omitted; threshold sources are footnote letters a-m as
printed: Zhu 2018, Giri 2010, Münch & Schieberle 1998 [MFT 0.007, FFT 0.01 µg/L], Meng 2012, Zhao 2017,
Pino & Queris 2011, Yu & Chen 2010, Dreher 2003, Niu 2019, Combet 2006, Song 2018, Javed 2018.)
Footnotes: Control I = no Cys/Xyl, heated; II = + 0.8 % Cys + 1.0 % Xyl, unheated; III = 0.8 % Cys +
1.0 % Xyl **dissolved in water**, heated; Experimental = hydrolysate + 0.8 % Cys + 1.0 % Xyl, heated.
Letters A-C: Duncan groups within a row (p < 0.05).

### Text claims [Q]
- "Only one set of the optimized conditions": 0.80 %, 183.80 °C, 89.34 min.
- Single-factor pre-screen (data not shown): maximum MRPs at 170 °C, 90 min, 0.7 % Cys.
- MFT/FFT "only detected in the heat-treated sample in the presence of cysteine and xylose".
- MFT-MFT "not considered as the response" because it may form by oxidation in storage or in the GC.
- DMDS/DMTS/MeSH only in Control I (methionine Strecker) — absent whenever cysteine is added.
- The authors write that xylose and cysteine form TTCA by "Amadori rearrangement" [!: TTCA is the
  thiazolidine, not the Amadori product] and cite Cao 2017 that TTCA forms "even at room temperature".

---

## 3. Derived quantities [Z]

**Centre-point replication (0.7 %, 170 °C, 90 min, n = 5):** MFT 117.4 ± 15.0 µg/L (pure-error SD √225),
CV 13 %; FFT 112.7 ± 12.2 µg/L, CV 11 %.

**Time behaviour at 0.7 % Cys (fish matrix):**

| T | MFT 60 → 120 min | FFT 60 → 120 min |
|---|---|---|
| 150 °C | 4.41 → 30.62 (×6.9, still forming) | 9.24 → 31.56 (×3.4) |
| 190 °C | 123.91 → 62.54 (÷2.0, net loss) | 131.71 → 91.46 (÷1.4) |

So by 190 °C the thiols pass a maximum before 60 min and decay with an apparent first-order net
loss of **~1.1 × 10⁻² min⁻¹ (MFT)** and **~6 × 10⁻³ min⁻¹ (FFT)** between 60 and 120 min (from
two points, assuming formation has stopped; a floor on the true loss rate). At 150 °C MFT is
still rising at 120 min.

**Molar yields at the optimum (RF = 1 against 2-octanol, so order-of-magnitude only):**

| species | Experimental µM | per mol Cys | Control III µM | per mol Cys |
|---|---|---|---|---|
| MFT (114.17) | 1.14 | 1.7 × 10⁻⁵ | 0.068 | **1.0 × 10⁻⁶** |
| FFT (114.17) | 1.18 | 1.8 × 10⁻⁵ | 0.0074 | **1.1 × 10⁻⁷** |
| MFT-MFT (226.3) | 0.455 (= 0.91 µM MFT-eq) | — | 0.0048 | — |
| furfural (96.08) | 1.27 | 1.9 × 10⁻⁵ per mol Xyl | 2.25 | 3.4 × 10⁻⁵ per mol Xyl |

MFT : FFT (molar): **0.97** (fish matrix) vs **9.2** (water, Control III). Matrix factor: MFT ×17, FFT
×160, MFT-MFT ×95, 2-thiophenethiol ×138, 3-methylthiophene ×197; furfural ×0.56.
(Headspace partitioning differs between a protein hydrolysate and water, so part of each factor is
SPME matrix effect, but a ×160 FFT effect against a ×0.56 furfural effect cannot be partitioning
alone.)

**Control III vs de Roos 2005:** MFT 1.0 × 10⁻⁶ mol/mol Cys at ~184 °C, pH ~4.2, 66 mM, ~90 min,
vs de Roos's 3.2 × 10⁻⁴ (pH 5, boiling SDE, ~4 mol/kg). Different extraction (HS-SPME standing
headspace vs SDE stripping), RF = 1 here, closed vs open; not directly comparable — but the
~300× gap is consistent with heavy thiol loss in a closed pot at 184 °C (cf. the 190 °C decay above).

---

## 4. What this can and cannot anchor in the model

**Can (within-study, qualitative to semi-quantitative):**
1. **A large matrix effect on Cys + Xyl thiol yields** (MFT ×17, FFT ×160 in fermented fish
   hydrolysate vs water, same Cys/Xyl, same heating). Direction and size are a hold-out for any
   protein-hydrolysate / free-amino-acid co-reactant effect (cf. de Roos 2005's glycine general-acid
   catalysis of the Cys-Amadori route; Cao 2017's glycine acceleration of TTCA drain; Zhang 2025b's
   lysine). The hydrolysate supplies a large excess of free amino acids, which would push the
   cysteine route through the Amadori compound and make the Strecker / dicarbonyl → H₂S chemistry
   richer. **Not attributable to a single mechanism here.**
2. **Thiol maximum and decay at high temperature**: at 190 °C MFT and FFT are past their maxima by
   60 min (net loss ~0.6-1.1 × 10⁻² min⁻¹); at 150 °C still rising at 120 min. A shape constraint for
   the MFT/FFT formation-loss balance at 150-190 °C.
3. **MFT : FFT ≈ 9 in water at pH ~4.2 and ~184 °C** (Control III; RF = 1 for two isomers) vs ≈ 1 in
   the matrix.

**Cannot:**
- **No absolute concentrations** (2-octanol IS amount unstated, RF = 1, no thiol calibration).
- No intermediates (TTCA, Cys-Amadori, dicarbonyls, H₂S), no free Cys/Xyl, no final pH.
- Vessel and headspace not stated; heating on a block heater at 150-190 °C (pressure vessel
  assumed).
- The RSM model omits the AC interaction (8 model df instead of 9 for a full quadratic BBD) and has
  C.V. 23-31 %; its surfaces are not kinetic information.
- Fish matrix composition (free amino acids, peptides, residual lipids, yeast metabolites) not
  quantified in this paper (in Gao 2019, not on disk).

---

## 5. Flags

1. **"Amadori rearrangement … to form TTCA"** is a mechanistic error in the paper; TTCA is the
   thiazolidine (de Roos 2005, Cao 2017).
2. Text says "The value of the C.V. % greater than 10 % means the model is reproducible" — the
   opposite of the usual criterion; C.V. 31 % (MFT) is high.
3. Control III's heating conditions are not stated explicitly (assumed = optimum).
4. Semi-quantification with RF = 1 and an unstated IS amount: µg/L are relative responses.
5. No heated-water blank for the hydrolysate's own sulfur (Control I covers it partly).
6. Registry: thieno[3,2-b]thiophene (365 µg/L, the largest S-volatile) and the
   thieno[2,3-b]thiophenes are MS + KI identifications only.

---

## 6. Relations to existing dossiers

| dossier | relation |
|---|---|
| `gao2020_extraction.md` | **Different paper** (pea protein isolate). Do not merge. |
| `deroos2005_extraction.md` | de Roos predicts that other amino acids catalyse the Cys + sugar → Cys-Amadori route and that Cys reacts with other amino acids' Amadori compounds — a candidate explanation for the ×17-160 matrix factor. |
| `cao2017_extraction.md`, `hou2017_extraction.md` | Cited by the paper for TTCA as the intermediate and for glycine accelerating TTCA → meaty compounds. |
| `zhang2025b_extraction.md` | Lysine (0-16 × Cys) modulates Cys + Xyl thiols at 120 °C, pH 6.5 — the controlled version of an amino-acid co-reactant effect. |
| `xu2011_extraction.md` | Cited (FFT favoured at lower pH; thiophenes from cysteine + lipid carbonyls); lipid suppresses, hydrolysate enhances. |

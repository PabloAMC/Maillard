# Zhang C, Tan, He, Hu, Li & Xie 2025 — EXTRACTION
### *"Effect of lysine on the cysteine-xylose Maillard reaction to form flavor compounds"*
### Cysteine + xylose 40 + 40 mmol/L in 0.2 mol/L phosphate pH 6.5, 120 °C, 50 min, with lysine at 0-16 × cysteine: 36 volatiles (MFT, FFT, bis(2-methyl-3-furyl) disulfide included) at 50 min, plus 0-50 min time courses of TTCA, Cys-Amadori, Lys-Amadori, free cysteine and free xylose. Read 2026-10-27.

**Source on disk:** `data/articles/zhang2025b.pdf` (publisher's PDF, 12 pp.). Text from `pdftotext -layout`;
Table 2 checked cell by cell against 150 dpi page images (pp. 4-5); Figure 2 (p. 7) read from 300 dpi
crops. The Supplementary file (Figs. S1-S2: structures, heat map) is **not on disk** and was not read.

**Not the same group as `zhang2025_extraction.md`.** That file is Zhang **H**., ..., Ho (Jiangnan/Rutgers,
di-xylulosyl-lysine ARP, Food Chem 479, 143828). This paper is Zhang **C**., ..., Xie (Beijing Technology
and Business University), the Xie group of Hou 2017 (`hou2017_extraction.md`), Cao 2017 and Zhao 2019a.

**Provenance:** **[M]** printed number · **[Q]** qualitative statement · **[F]** read off a figure by me
(read error stated) · **[Z]** derived by me from the paper's numbers, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Chenping Zhang, Jia Tan, Jinxin He, Qingqing Hu, Juan Li, Jianchun Xie* |
| Title | "Effect of lysine on the cysteine-xylose Maillard reaction to form flavor compounds" |
| Venue | *Food Chemistry* **469** (2025) **142529** |
| DOI | 10.1016/j.foodchem.2024.142529 |
| Affiliation | School of Light Industry Science and Engineering, School of Food Science and Health, Beijing Technology and Business University (BTBU), Beijing 100048 |
| Dates | received 7 Aug 2024; revised 17 Nov 2024; accepted 14 Dec 2024; online 24 Dec 2024 |
| Data availability | "Data will be made available on request." |

---

## 1. THE ONE-PARAGRAPH ANSWER

0.2 mmol cysteine + 0.2 mmol xylose in 5 mL 0.2 mol/L phosphate pH 6.5 (**40 + 40 mmol/L** [Z]) with
lysine at 0, 1, 2, 4, 6, 8 or 16 × cysteine, sealed 15 mL pressure vials, **120 °C, 50 min**. In the
cysteine-xylose blank, **2-methyl-3-furanthiol (MFT) is 22.63 ng/mL, 2-furfurylthiol (FFT) 0.29 ng/mL
and bis(2-methyl-3-furyl) disulfide 18.25 ng/mL** at 50 min (SPME-GC-MS, internal standard; Table 2).
Lysine at 1 × cysteine raises MFT to 32.56 and the disulfide to 51.98 ng/mL; more lysine lowers both
(MFT 0.58 ng/mL at 16 ×, disulfide not detected at ≥ 8 ×). FFT falls monotonically and is not
detected at ≥ 6 ×. Pyrazines rise ~20-fold (84 → 1660 ng/mL). In the blank, **TTCA accumulates to
~18 mmol/L and Cys-Amadori to ~1.6 mmol/L** by 40-50 min (Fig. 2), i.e. ~45 % of the cysteine
charged sits as TTCA [Z]. With lysine, both peak at 10-15 min and are nearly gone by 50 min. Apparent
first-order loss of free cysteine and xylose (fits printed on Fig. 2g-h): **blank k_Cys = 0.059,
k_Xyl = 0.073 min⁻¹**; at 16 × lysine 0.328 and 0.547 min⁻¹.

---

## 2. METHODS [M] (pp. 2-3)

| item | as printed |
|---|---|
| Reactants | L-cysteine 99 %, L-lysine 99 % (**free base or HCl not stated** [!]), D-xylose 99 % (J&K) |
| Model reaction | 15 mL pressure glass vial; **0.2 mmol cysteine, 0.2 mmol xylose**, lysine at **0, 1, 2, 4, 6, 8 or 16 ×** the cysteine amount, **5 mL phosphate buffer 0.2 mol/L, pH 6.5** [Z: 40 mmol/L each; lysine 40-640 mmol/L] |
| Heating | **120 °C, total 50 min**, Parallel Synthesis Poly-Block 4 (HEL). Volatiles at the end only; other analyses on vials "taken at intervals"; ice-water quench. **No come-up time stated** [!] |
| pH after mixing / after heating | **not reported** [!]. 640 mmol/L lysine (16 ×) in 0.2 mol/L buffer may move the start pH if added as free base |
| Replicates | n = 3; Duncan's test p < 0.05; Pearson correlation in SPSS 24 |
| SPME | 1 µL internal standard (**1,2-dichlorobenzene 73 µg/mL in methanol**) + 2.0 g NaCl; equilibrate **50 °C, 15 min**; extract **20 min**, "50/35 μm carboxen/polydimethylsiloxane/divinylbenzene" fiber (as printed) |
| GC-MS | Agilent 7890A/5975C, DB-Wax and HP-5MS, 70 eV, m/z 50-450 |
| Quantification | internal-standard ratio; no response factors or calibration curves printed for volatiles. Table 2 footnote: **"concentrations were the sums detected from both the DB-WAX and HP-5 MS capillary columns"** [!] (meaning unclear: summed across two separate injections would double-count) |
| GC-O | AEDA by split-ratio steps, HP-5MS column, 3 panellists, run on the 1 × lysine system only |
| Intermediates | TTCA and Cys-Amadori prepared per Hou 2017 (> 98 %), ε-Lys-Amadori from Fmoc-Lys-OH (> 95 %); quantified by LC-MS/MS MRM, external standards (Table 1). α-Lys-Amadori approximated with the ε-isomer curve; Fig. 2 Lys-Amadori is the α + ε sum |
| Cysteine, xylose | HPLC-ELSD, XBridge Amide, external standard. **Both calibration curves are printed identical: y = 263.28 x − 58.733, R² = 0.99** (x mmol/L) [!] (one is evidently a copy) |
| Spectrophotometry | A294, A420 (UV2300), ΔE (Hunter UltraScan PRO) |

**Table 1** (LC-MS/MS MRM; x = concentration mmol/L, y = peak area), verbatim:

| Compound | Quantitative ion (m/z) | Cone voltage (V) | Collision energy (eV) | Standard concentrations (mmol/L) | Calibration curve |
|---|---|---|---|---|---|
| TTCA | 254 > 122 | 80 | 10 | 0.1-20 | y = 2418 x, R² = 0.993 |
| Cys-Amadori | 254 > 218 | 80 | 8 | 0.1-3 | y = 6703 x, R² = 0.991 |
| ε-Lys-Amadori | 279 > 225 | 110 | 15 | 0.1-10 | y = 5000 x, R² = 0.994 |
| α-Lys-Amadori | 279 > 84 | 110 | 10 | – | – (ε-isomer curve used) |

---

## 3. RESULTS

### 3.1 Table 2 (pp. 4-5): volatiles at 120 °C, 50 min, ng/mL, mean ± SD (n = 3), Duncan letters within a row [M]

"–" = not detected. Columns: lysine at 0 (blank), 1, 2, 4, 6, 8, 16 × cysteine. RI on DB-Wax / HP-5MS.
Odor and FD from AEDA on the 1 × system. ID: MS, RI, O, STD (authentic standard).

| No. | RI Wax | RI HP-5 | Compound | blank | ×1 | ×2 | ×4 | ×6 | ×8 | ×16 | Odor | FD | ID |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| | | | **Sulfur-containing** | | | | | | | | | | |
| 1 | 1024 | 662 | thiophene | 3.80 ± 0.58b | 6.42 ± 0.73a | 3.07 ± 0.44bc | 2.48 ± 0.44c | 0.88 ± 0.15d | 0.58 ± 0.15d | 0.29 ± 0.02e | meaty | 2 | MS, RI, O, STD |
| 2 | 1090 | 767 | 2-methylthiophene | 2.34 ± 0.29e | 13.72 ± 2.63a | 7.15 ± 1.02b | 4.53 ± 0.58c | 3.50 ± 0.44d | 1.02 ± 0.15f | – | sulfury | 8 | MS, RI, O, STD |
| 3 | 1168 | 1114 | 1,2,3-trithiolane | 10.07 ± 1.90cd | 32.85 ± 5.26a | 40.44 ± 5.11a | 22.34 ± 2.92b | 13.43 ± 1.61c | 9.2 ± 1.17d | 4.82 ± 0.58e | garlic | 8 | MS, RI, O |
| 4 | 1231 | 781 | thiazole | 30.08 ± 5.40a | 16.94 ± 3.07b | 11.53 ± 2.04c | 7.59 ± 1.31d | 6.13 ± 1.17d | 4.23 ± 0.73e | 0.88 ± 0.15f | meaty | 2 | MS, RI, O, STD |
| 5 | 1268 | 807 | 2-methylthiazole | 34.46 ± 3.80a | 23.07 ± 2.48b | 13.43 ± 1.46c | 9.64 ± 1.02d | 5.11 ± 0.58e | 4.09 ± 0.44f | 3.34 ± 0.29g | roasted | 2 | MS, RI, O, STD |
| 6 | 1280 | 813 | 3-mercapto-2-butanone | 0.58 ± 0.15a | 0.45 ± 0.15a | 0.44 ± 0.15a | 0.29 ± 0.02b | 0.26 ± 0.02b | 0.15 ± 0.01c | – | meaty | 16 | MS, RI, O, STD |
| 7 | 1298 | 1322 | [1,2,3,4]tetrathiane | 5.40 ± 0.73d | 11.24 ± 1.61b | 21.46 ± 2.77a | 12.26 ± 1.46b | 8.18 ± 1.31c | 4.67 ± 0.88d | 2.19 ± 0.29e | sulfury | 4 | MS, RI, O |
| 8 | 1328 | 865 | **2-methyl-3-furanthiol** | **22.63 ± 4.09b** | **32.56 ± 3.94a** | **21.75 ± 4.09b** | **9.49 ± 1.46c** | **5.26 ± 0.73d** | **3.21 ± 0.44e** | **0.58 ± 0.04f** | meaty | 256 | MS, RI, O, STD |
| 9 | 1336 | 902 | 3-mercapto-2-pentanone | 0.29 ± 0.02a | 0.15 ± 0.01b | 0.12 ± 0.02b | 0.10 ± 0.01b | – | – | – | meaty | 256 | MS, RI, O, STD |
| 10 | 1376 | 897 | 2,3-dihydro-5-methythiophene | 4.67 ± 0.58a | 2.63 ± 0.44b | 1.9 ± 0.29c | 0.88 ± 0.15d | – | – | – | meaty | 2 | MS, RI, O |
| 11 | 1422 | 909 | **2-furfurylthiol** | **0.29 ± 0.03a** | **0.25 ± 0.02a** | **0.16 ± 0.01b** | **0.15 ± 0.02b** | – | – | – | meaty | 256 | MS, RI, O, STD |
| 12 | 1511 | 968 | 3-thiophenethiol | 17.37 ± 2.48a | 4.67 ± 0.73b | 1.75 ± 0.15c | 1.31 ± 0.29c | 0.35 ± 0.01d | 0.29 ± 0.02e | – | meaty | 256 | MS, RI, O, STD |
| 13 | 1525 | 987 | dihydro-2-methy-3(2H)-thiophenone | 4.67 ± 0.58e | 28.03 ± 4.23b | 37.67 ± 4.53a | 19.86 ± 2.34c | 10.95 ± 1.31d | 5.11 ± 0.58e | 3.21 ± 0.44f | meaty | 2 | MS, RI, O, STD |
| 14 | 1600 | 1129 | 3,5-dimethyl-1,2,4-trithiolane | 0.89 ± 0.06e | 3.50 ± 0.44d | 5.26 ± 0.73c | 6.42 ± 0.88c | 10.37 ± 1.31b | 12.12 ± 1.46ab | 13.72 ± 1.75a | meaty | 8 | MS, RI, O |
| 15 | 1630 | 1020 | 2-acetylthiazole | 114.13 ± 15.98d | 241.15 ± 33.75a | 215.00 ± 30.12a | 186.46 ± 26.11b | 170.15 ± 23.82bc | 140.91 ± 19.73c | 107.75 ± 15.10d | roasted, popcorn | 2 | MS, RI, O, STD |
| 16 | 1686 | 990 | 3-thiophenecarboxaldehyde | 44.97 ± 8.03a | 16.79 ± 1.90b | 9.78 ± 1.75c | 3.36 ± 0.44d | 1.61 ± 0.29e | 1.02 ± 0.15f | 0.29 ± 0.03g | sulfury | 16 | MS, RI, O, STD |
| 17 | 1711 | 1086 | 2-thiophenemethenethiol | 12.56 ± 1.61b | 24.67 ± 3.21a | 12.85 ± 2.04b | 11.10 ± 1.90b | 7.74 ± 1.31c | 2.48 ± 0.44d | 1.17 ± 0.15e | meaty | 16 | MS, RI, O, STD |
| 18 | 1721 | 1001 | 2-thiophenecarboxaldehyde | 28.91 ± 4.09a | 17.96 ± 2.92b | 14.31 ± 2.48b | 9.78 ± 1.17c | 7.15 ± 0.88d | 4.67 ± 0.58e | 2.04 ± 0.44f | sulfury | 32 | MS, RI, O, STD |
| 19 | 1744 | 1115 | 5-methyl-2-thiophenecarboxaldehyde | 9.34 ± 1.46e | 35.62 ± 4.96a | 37.38 ± 4.09a | 24.53 ± 4.38b | 16.06 ± 2.92c | 13.43 ± 1.61d | 10.80 ± 1.9e | almond | 64 | MS, RI, O, STD |
| 20 | 1776 | 1120 | 3-methyl-2-thiophenecarboxaldehyde | 71.69 ± 12.85b | 92.86 ± 13.87ab | 109.06 ± 13.14a | 74.02 ± 8.91b | 43.07 ± 7.30c | 20.88 ± 3.80d | 9.78 ± 1.17e | almond | 16 | MS, RI, O, STD |
| 21 | 1801 | 1336 | 2,5-thiophenedicarboxaldehyde | 3.36 ± 0.44b | 3.50 ± 0.58b | 7.45 ± 0.88a | 4.09 ± 0.44b | 1.75 ± 0.15c | – | – | sulfury | 1 | MS, RI, O, STD |
| 22 | 1812 | 1180 | 1-(3-thienyl)-1-propanone | 7.30 ± 1.02a | 3.36 ± 0.58b | 2.77 ± 0.29bc | 2.34 ± 0.29c | 1.31 ± 0.15d | 0.88 ± 0.15e | 0.44 ± 0.02f | onion | 2 | MS, RI, O |
| 23 | 1821 | 1203 | thieno[2,3-b]thiophene | 23.07 ± 4.09a | 15.62 ± 2.34b | 12.70 ± 1.46b | 7.88 ± 0.88c | 1.75 ± 0.15d | 1.17 ± 0.15e | – | meaty | 4 | MS, RI, O, STD |
| 24 | 1833 | 1209 | thieno[3,2-b]thiophene | 15.48 ± 2.63c | 20.29 ± 2.63b | 32.56 ± 4.82a | 28.32 ± 4.23ab | 21.17 ± 3.80b | 12.85 ± 2.34c | 6.72 ± 1.02d | meaty | 4 | MS, RI, O, STD |
| 25 | 2145 | 1525 | **bis(2-methyl-3-furyl) disulfide** | **18.25 ± 2.77b** | **51.98 ± 6.72a** | **16.94 ± 2.77b** | **4.67 ± 0.73c** | **2.48 ± 0.44c** | – | – | meaty | 256 | MS, RI, O, STD |
| | | | Subtotal (S) | 486.60 | 700.28 | 636.93 | 453.89 | 338.66 | 242.96 | 168.02 | | | |
| | | | **Nitrogen-containing** | | | | | | | | | | |
| 26 | 1213 | 739 | pyrazine | 4.67 ± 0.88e | 36.35 ± 4.38d | 56.65 ± 7.15c | 61.76 ± 8.03bc | 70.81 ± 8.76b | 76.80 ± 9.64ab | 93.73 ± 10.68a | roasted | 1 | MS, RI, O |
| 27 | 1315 | 915 | 2,5-dimethylpyrazine | 13.14 ± 1.9g | 27.59 ± 4.96f | 52.12 ± 9.34e | 80.15 ± 10.45d | 149.8 ± 27.01c | 232.14 ± 41.76b | 496.11 ± 89.35a | roasted | 8 | MS, RI, O, STD |
| 28 | 1334 | 921 | 2,6-dimethylpyrazine | 16.50 ± 1.90d | 40.00 ± 5.55c | 38.98 ± 5.40c | 91.98 ± 14.31b | 118.7 ± 16.64ab | 133.59 ± 19.27a | 141.91 ± 21.61a | roasted | 8 | MS, RI, O, STD |
| 29 | 1339 | 927 | 2,3-dimethylpyrazine | 14.60 ± 2.34e | 20.44 ± 2.19d | 26.13 ± 2.92c | 29.93 ± 4.23c | 38.54 ± 4.53b | 40.00 ± 4.38b | 65.41 ± 9.20a | roasted | 8 | MS, RI, O, STD |
| 30 | 1391 | 991 | 2-ethyl-6-methylpyrazine | 20.29 ± 2.77e | 155.78 ± 25.70d | 230.39 ± 36.79c | 325.29 ± 35.77b | 389.09 ± 62.20ab | 431.72 ± 69.06a | 474.65 ± 52.27a | roasted | 16 | MS, RI, O, STD |
| 31 | 1411 | 1001 | 2,3,5-trimethylpyrazine | – | 8.18 ± 1.02d | 11.10 ± 1.31c | 11.83 ± 1.90c | 22.92 ± 2.77b | 27.16 ± 3.21b | 46.28 ± 7.45a | roasted | 16 | MS, RI, O, STD |
| 32 | 1457 | 1089 | 3-ethyl-2,5-dimethylpyrazine | 8.91 ± 1.02f | 24.24 ± 2.92e | 29.49 ± 3.50de | 34.75 ± 4.82d | 52.85 ± 6.42c | 133.88 ± 16.79b | 308.94 ± 43.51a | roasted | 32 | MS, RI, O |
| 33 | 1501 | 1128 | 3,5-diethyl-2-methylpyrazine | 6.13 ± 1.17c | 7.01 ± 1.31c | 11.83 ± 1.46b | 13.72 ± 2.04b | 15.77 ± 2.92b | 27.30 ± 3.50a | 32.70 ± 3.80a | roasted | 32 | MS, RI, O |
| | | | Subtotal (N) | 84.24 | 319.59 | 456.69 | 649.41 | 858.48 | 1102.59 | 1659.73 | | | |
| | | | **Oxygen-containing heterocyclic** | | | | | | | | | | |
| 34 | 879 | 606 | 2-methylfuran | 4.38 ± 3.07e | 18.25 ± 1.46a | 15.18 ± 1.17b | 11.83 ± 0.88c | 8.61 ± 0.73d | 8.32 ± 0.73d | 3.21 ± 0.29e | – | – | MS, RI |
| 35 | 1246 | 954 | 2-(methoxymethyl)-furan | 1.61 ± 0.29e | 36.94 ± 2.92a | 36.21 ± 2.91a | 29.35 ± 2.34b | 24.67 ± 1.90c | 21.32 ± 1.75c | 8.76 ± 0.73d | – | – | MS, RI |
| 36 | 1666 | 857 | 2-furanmethanol | 4.53 ± 0.88c | 9.34 ± 0.73a | 6.86 ± 0.58b | 3.36 ± 0.29cd | 2.92 ± 0.29d | 1.75 ± 0.15e | 0.58 ± 0.05f | bready | 2 | MS, RI, O, STD |
| | | | Subtotal (O) | 10.52 | 64.53 | 58.25 | 44.54 | 36.2 | 31.39 | 12.55 | | | |
| | | | **Total** | 581.36 | 1084.40 | 1151.87 | 1147.84 | 1233.34 | 1376.94 | 1840.30 | | | |

Subtotals and totals re-add exactly from the rows (blank, 1 × and 16 × checked) [Z]. No furfural is
reported in any system [!] (a cysteine-xylose pot at 120 °C normally gives it; the SPME/GC window or the
compound list excluded it). Footnote: odor/FD by AEDA on HP-5MS; letters are within-row.

**Molar scale [Z]** (MW: MFT 114.17, FFT 114.17, bis(2-methyl-3-furyl) disulfide 226.31; per mol cysteine
charged, 40 mmol/L):

| system | MFT µmol/L | MFT mol/mol Cys | FFT nmol/L | FFT mol/mol Cys | disulfide µmol/L | MFT + 2 × disulfide, mol/mol Cys |
|---|---|---|---|---|---|---|
| blank | 0.198 | 5.0 × 10⁻⁶ | 2.5 | 6.4 × 10⁻⁸ | 0.081 | 9.0 × 10⁻⁶ |
| 1 × Lys | 0.285 | 7.1 × 10⁻⁶ | 2.2 | 5.5 × 10⁻⁸ | 0.230 | 1.9 × 10⁻⁵ |
| 16 × Lys | 0.0051 | 1.3 × 10⁻⁷ | n.d. | – | n.d. | 1.3 × 10⁻⁷ |

MFT : FFT ≈ **78 : 1** in the blank (mass ratio 22.63/0.29), unusually MFT-heavy for a pentose pot.

### 3.2 Figure 2a-c (p. 7): early-stage intermediates, mmol/L [F]

Read from the 300 dpi crop; read error **±0.05 mmol/L (Cys-Amadori, Lys-Amadori) and ±0.3 mmol/L
(TTCA)**; error bars are small and not transcribed. Sampling at 0, 5, 10, 15, 20, 25, 30, 40, 50 min.

**Cys-Amadori (N-(1-deoxy-D-xylulos-1-yl)-L-cysteine):**

| t (min) | blank | ×1 | ×2 | ×4 | ×6 | ×8 | ×16 |
|---|---|---|---|---|---|---|---|
| 5 | ~0.00 | 0.10 | 0.05 | 0.03 | 0.02 | 0.01 | 0.00 |
| 10 | 0.07 | 0.85 | 0.60 | 0.43 | 0.35 | 0.25 | 0.07 |
| 15 | 1.05 | **1.35** | 1.10 | 0.78 | 0.64 | 0.51 | 0.33 |
| 20 | 1.38 | 1.19 | 0.84 | 0.65 | 0.52 | 0.41 | 0.21 |
| 25 | 1.52 | 0.95 | 0.60 | 0.49 | 0.35 | 0.29 | 0.11 |
| 30 | 1.60 | 0.62 | 0.35 | 0.27 | 0.20 | 0.11 | 0.04 |
| 40 | 1.61 | 0.22 | 0.12 | 0.08 | 0.06 | 0.03 | 0.01 |
| 50 | 1.61 | 0.08 | 0.03 | 0.01 | ~0 | ~0 | ~0 |

**TTCA (2-threityl-thiazolidine-4-carboxylic acid):**

| t (min) | blank | ×1 | ×2 | ×4 | ×6 | ×8 | ×16 |
|---|---|---|---|---|---|---|---|
| 5 | 7.4 | 7.6 | 7.3 | 6.9 | 6.2 | 5.4 | 3.7 |
| 10 | 11.2 | **13.1** | 11.6 | 10.5 | 9.3 | 8.0 | 5.2 |
| 15 | 13.5 | 10.6 | 8.4 | 7.0 | 5.6 | 4.7 | 2.4 |
| 20 | 15.1 | 7.1 | 5.6 | 4.3 | 3.8 | 3.1 | 1.4 |
| 25 | 16.4 | 4.1 | 3.2 | 2.3 | 2.0 | 1.6 | 0.6 |
| 30 | 17.7 | 1.6 | 1.3 | 1.0 | 0.8 | 0.5 | 0.3 |
| 40 | 18.1 | 0.3 | 0.3 | 0.2 | ~0.2 | ~0.1 | 0.1 |
| 50 | 18.2 | ~0 | ~0 | ~0 | ~0 | ~0 | ~0 |

**Lys-Amadori (α + ε, ε calibration):**

| t (min) | ×1 | ×2 | ×4 | ×6 | ×8 | ×16 |
|---|---|---|---|---|---|---|
| 5 | 0.55 | 0.88 | 1.58 | 2.07 | 2.78 | **5.12** |
| 10 | 0.68 | 1.15 | 2.15 | 2.02 | 1.85 | 3.02 |
| 15 | 0.55 | 0.78 | 1.50 | 1.35 | 1.27 | 1.87 |
| 20 | 0.32 | 0.45 | 0.80 | 0.85 | 0.88 | 1.05 |
| 25 | 0.18 | 0.25 | 0.33 | 0.35 | 0.55 | 0.62 |
| 30 | 0.10 | 0.15 | 0.20 | 0.25 | 0.33 | 0.38 |
| 40 | 0.05 | 0.07 | 0.08 | 0.10 | 0.15 | 0.17 |
| 50 | ~0 | ~0 | ~0 | ~0 | ~0 | ~0 |

Text anchors [M]: blank TTCA and Cys-Amadori "increased rapidly to 30 min and then climbed up slowly";
with lysine Cys-Amadori peaks at **15 min**, TTCA at **10 min**; Lys-Amadori peak shifts from **10 to
5 min** above 4 × lysine. TTCA ↔ Cys-Amadori r = 0.97; both vs Lys-Amadori r < −0.90 (sign implied).

**[!] Text vs figure.** The text says Cys-Amadori and TTCA "were less formed as lysine reacted with
xylose against the cysteine", but at 1 × lysine both **exceed** the blank at 10-15 min (TTCA 13.1 vs
11.2 at 10 min; Cys-Amadori 1.35 vs 1.05 at 15 min). Lysine lowers the late pools, not the early peaks
at low dose.

### 3.3 Figure 2d-h (p. 7): rate fits as printed on the figure [M]

Zero-order fits through the origin for A294, A420 and ΔE (y in AU or ΔE units, x in min); first-order
fits of ln(C/C₀) for xylose and cysteine (y = −k x + b). All printed on the panels; assignment of lines
to series by colour, and checked against the text ratios (below).

| lysine | k294 (AU/min) | R² | k420 (AU/min) | R² | kΔE (/min) | R² | k_Xyl (min⁻¹) | b_Xyl | R² | k_Cys (min⁻¹) | b_Cys | R² |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| blank | 0.0581 | 0.952 | 0.004 | 0.989 | 0.016 | 0.981 | **0.073** | −0.324 | 0.969 | **0.059** | −0.760 | 0.976 |
| ×1 | 0.485 | 0.984 | 0.030 | 0.993 | 0.156 | 0.965 | 0.169 | −0.337 | 0.962 | 0.085 | −0.736 | 0.900 |
| ×2 | 0.742 | 0.981 | 0.048 | 0.989 | 0.225 | 0.977 | 0.259 | −0.307 | 0.984 | 0.097 | −0.747 | 0.900 |
| ×4 | 1.102 | 0.983 | 0.082 | 0.989 | 0.364 | 0.972 | 0.330 | −0.303 | 0.981 | 0.120 | −0.792 | 0.928 |
| ×6 | 1.313 | 0.995 | 0.130 | 0.983 | 0.506 | 0.970 | 0.370 | −0.305 | 0.981 | 0.157 | −0.790 | 0.919 |
| ×8 | 1.719 | 0.998 | 0.180 | 0.995 | 0.698 | 0.996 | 0.465 | −0.310 | 0.982 | 0.237 | −0.632 | 0.943 |
| ×16 | 2.987 | 0.908 | 0.270 | 0.998 | 0.878 | 0.998 | 0.547 | −0.336 | 0.981 | 0.328 | −0.900 | 0.975 |

Text [M] (p. 8): at 1 × lysine k294, k420, kΔE, k_Xyl, k_Cys are "8.4, 7.5, 9.8, 2.3, and 1.4 times" the
blank; at 16 × "51.5, 67.5, 54.9, 7.5, and 5.6 times". These reproduce from the slopes above (0.485/0.0581
= 8.3; 0.030/0.004 = 7.5; 0.156/0.016 = 9.75; 0.169/0.073 = 2.3; 0.085/0.059 = 1.44; 2.987/0.0581 =
51.4; 0.547/0.073 = 7.5; 0.328/0.059 = 5.6) [Z], so the colour assignment is right.

**[!] Fit caveats.** (i) The intercepts are far from 0 although the t = 0 point is plotted at
ln(C/C₀) = 0: about **−0.76 for cysteine (≈ 53 % gone in the first ~5 min) and −0.32 for xylose
(≈ 27 %)** [Z]. The printed k describes the 5-50 min phase after a fast initial drop (TTCA
formation). (ii) The 16 × cysteine and xylose fits stop at 15-25 min (concentrations below detection).
(iii) The A294 16 × points curve over (R² 0.908); the zero-order slope is not a rate there. (iv) With
lysine present, cysteine and xylose loss include TTCA, Cys-Amadori and Lys-Amadori formation; these are
total-loss constants, not single-step constants. (v) The cysteine and xylose calibrations are printed
identical (section 2).

### 3.4 Mass balance in the blank [Z]

At 50 min in the blank: free cysteine e^(−0.059·50 − 0.760) = **2.4 %** of 40 mmol/L (≈ 1.0 mmol/L);
free xylose e^(−0.073·50 − 0.324) = **1.9 %** (≈ 0.8 mmol/L). TTCA 18.2 + Cys-Amadori 1.6 = 19.8 mmol/L
hold about half the cysteine and half the xylose. **~19 mmol/L of each reactant is not in any measured
pool**, while A420 is only 0.2 AU and total volatiles 0.58 µg/mL. Candidates (none measured): cystine,
Schiff base / open-chain TTCA that the MRM transition does not catch, and H₂S. The fitted k_Cys and
k_Xyl in the blank are therefore **not** decay to products the model would recognise; they are
disappearance from the HPLC-ELSD peak.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** The sulfur lane carries `Cys`, `PENT` (pentose), `TTCA` with two fates
(`r_ttca_cys` TTCA → Cys + PENT; `r_ttca_deg`), and MFT/FFT/disulfide products
(`src/kinetic_core/sulfur.py` ~line 577; `species_sulfur.py` line 160). TTCA is a **fed** species
(Kang 2026 at 10 mmol/L); there is **no Cys + PENT → TTCA forward step** in `sulfur.py` (checked by
grep, 2026-10-27). Lysine as such does not appear in the sulfur lane. Nothing below was run against
the model.

**(a) Usable as a held-out TEST pot (not a fit target):** cysteine + xylose 40 + 40 mmol/L, 0.2 mol/L
phosphate pH 6.5 (start), 120 °C, 50 min, sealed:
- MFT **22.63 ± 4.09 ng/mL**, FFT **0.29 ± 0.03 ng/mL**, bis(2-methyl-3-furyl) disulfide **18.25 ±
  2.77 ng/mL**, 3-mercapto-2-pentanone 0.29, 3-mercapto-2-butanone 0.58, thiophene 3.80, 2-acetylthiazole
  114.13 ng/mL (Table 2 blank). Semi-quantitative (internal standard, response factor presumably 1,
  "sum of both columns"); treat as order-of-magnitude, ±× 3 at least.
- Directional test: lysine at 1 × cysteine **raises** MFT (+44 %) and the MFT disulfide (×2.8) while
  lowering FFT; ≥ 4 × lysine collapses all thiols. A model with a generic competing amine could be asked
  only for the direction.

**(b) Usable as a SHAPE / partition check:** in a cysteine-xylose pot at 120 °C, pH 6.5, **~45 % of the
cysteine is held as TTCA (18 mmol/L) after 30-50 min and only ~4 % as Cys-Amadori** (Fig. 2). This is
measured, and it is the sulfur-availability question the engine cannot answer, because it does not form
TTCA from Cys + PENT. Any engine pot of cysteine + pentose at pH 6.5 with all cysteine free is
overstating the sulfur available to H₂S/MFT routes (by up to ~2× here). The blank TTCA rising to a
plateau, not decaying, over 50 min at 120 °C also bounds net TTCA loss in a pot that keeps re-forming it;
compare `zhai2021_extraction.md` (10 mmol/L fed TTCA at pH 6, 120 °C) before using it.

**(c) Rate constants:** k_Cys 0.059 and k_Xyl 0.073 min⁻¹ (blank, 120 °C) are **apparent total-loss
constants after a fast initial drop**, with half of each reactant going to unmeasured pools (section
3.4). They cannot pin `k_cys_h2s` (Zheng & Ho thermolysis) or a condensation constant. At most: the
early drop (−0.76 in ln units within ~5 min) says Cys + Xyl condensation at 40 mmol/L, 120 °C is
**fast, t½ of order minutes**, consistent with the authors' TTCA reaching 7.4 mmol/L at 5 min.

**(d) Not usable:** A294/A420/ΔE zero-order slopes (no model observable); Pearson correlations and the
proposed pathways of Fig. 4 (no labelling; the authors say isotope work is needed); GC-O FD values.

**Caveats, collected.**
1. Start pH only; no pH after lysine addition or after heating. Lysine form (free base/HCl) unstated.
2. No come-up time; 15 mL pressure vials in a heating block, 5 mL fill.
3. Volatile quantification is internal-standard semi-quant with an unexplained two-column sum.
4. Cysteine and xylose calibrations printed identical; the mass balance leaves ~half of each reactant
   unaccounted.
5. Text claims for Cys-Amadori/TTCA ("less formed" with lysine) contradict Fig. 2 at 1 × lysine.
6. No furfural reported; MFT : FFT ≈ 78 : 1, which is unusual for a pentose + cysteine pot and should
   be checked against other Xie-group pots (Hou 2017, Cao 2017) before weighting FFT.
7. Supplementary (Figs. S1-S2) not read; they are a structure figure and a heat map, no new numbers
   expected.

# Fan et al. 2018 — EXTRACTION
### *"Aroma Compounds in Chicken Broths of Beijing Youji and Commercial Broilers"*
### Whole ~1 kg carcass stewed in 1500 mL water, crock pot, ~100 °C, 3 h, atmospheric pressure. SAFE + GC-MS(SIM) against one internal standard with SAFE recovery factors, ng/g broth. **Beijing-You (BYJ) broth: MFT 1.07, FFT 2.71, methional 158, MFT disulfide 1.62, hexanal 107 ng/g. Commercial broiler (CB): MFT 1.48, methional 118, hexanal 141 ng/g; FFT not quantified.** Plus breast-muscle ribose, glucose, thiamine, cysteine, methionine and fatty acids. Read 2026-10-27.

**Full read of `data/articles/fan2018.pdf`** (37 pp.: the ACS "Just Accepted" manuscript, not the
typeset version; tables at the end, pp. 26-34). Text layer clean; all five tables transcribed from
`pdftotext -layout` and checked by recomputing every OAV in Table 5 from Tables 4 and 5 (one
inconsistency, flag 1). Figure 1 (aroma-profile spider plots) not digitised.

**Provenance:** **[M]** printed number · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.
1 ng/g broth = 1 µg/kg.

---

## 0. IDENTITY — matches

| field | value |
|---|---|
| Authors | Mengdie Fan, Qunfei Xiao, Jianchun Xie*, Jie Cheng, Baoguo Sun, Wenbin Du, Yaxin Wang, Tianze Wang |
| Title | "Aroma Compounds in Chicken Broths of Beijing Youji and Commercial Broilers" |
| Venue | *J. Agric. Food Chem.* — the PDF is the Just Accepted manuscript ("Publication Date (Web): 08 Sep 2018"), volume/pages not printed. Crossref 2026-10-27: **66 (39), 10242-10251**, issue date 2018-10-03 |
| DOI | 10.1021/acs.jafc.8b03297 (printed on p. 1) |
| Affiliation | Beijing Technology and Business University (BTBU, Xie group); CAAS Institute of Quality Standard and Testing Technology |
| Companions | Same group as Cao 2017, Hou 2017, Zhao 2019, Du 2023/2024 (`cao2017_extraction.md`, `du2023_extraction.md`, `du2024_extraction.md`). Cited as top priority in `sohail2022_extraction.md` |

---

## 1. THE ONE-PARAGRAPH ANSWER

A stewed chicken broth contains, in absolute terms, **MFT ≈ 1.1-1.5 µg/kg, FFT ≈ 2.7 µg/kg (BYJ
only), bis(2-methyl-3-furyl) disulfide ≈ 1.6 µg/kg (BYJ only), 2-methyl-3-(methylthio)furan 0.5
µg/kg, methional ≈ 120-160 µg/kg, 2-acetylthiazole 3.7-6.3 µg/kg, hexanal ≈ 110-140 µg/kg and
(E,E)-2,4-decadienal ≈ 70-105 µg/kg**. On thresholds in water these give the highest OAVs of the
broth: methional ~3000-4000, MFT ~2700 (BYJ), (E,E)-2,4-decadienal ~1400-2100, MFT disulfide 810,
FFT 543. Recombination of the OAV ≥ 1 odorants reproduced both broths (4.6-4.7 of 5), and omission
of every one except decanal was significant. The slower-growing BYJ hen has more ribose (1.37 vs
1.09 mg/100 g), glucose, thiamine (0.19 vs 0.12 mg/100 g), cysteine (4.26 vs 3.56 mg/g, acid
hydrolysate) and PUFA, and its broth has more sulfur odorants by number; but **CB had the higher
MFT** (1.48 vs 1.07 ng/g), so the precursor argument does not carry to MFT.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Birds | frozen eviscerated carcasses from a Beijing supermarket; **BYJ hens fed 120 d; CB fed 42 d**; "ca. 1000 g"; −20 °C storage |
| Broth | "The whole chicken carcass was stewed in **1500 mL of water** in a crock pot for **3 h at ca. 100 °C** and atmosphere pressure." No salt, no lid status, no final broth mass, no pH [!] |
| Replication | six carcasses per breed -> six broths; concentrates pooled in pairs -> **three** samples for quantification (so n = 3 means 3 pooled pairs) |
| Isolation | SAFE; distillate extracted 3 × 1500 mL CH2Cl2; Na2SO4; Vigreux to 0.4 mL |
| Screening | AEDA/GC-O; quantification only for log2FD ≥ 7 (BYJ) or ≥ 3 (CB) with an authentic standard |
| Quantification | GC-MS SIM on DB-5MS; **one internal standard (1,2-dichlorobenzene)** with **compound-specific calibration curves** (Table 1) and **SAFE recovery factors** measured in an artificial broth (odourless broth residue + water + standards at ~native levels, triplicate). Concentration = calibration × extraction yield / RF. **Not stable-isotope dilution.** |
| Thresholds | orthonasal in water, copied from van Gemert 2011 (Table 5) |
| Meat analyses | breast muscle (M. pectoralis major): moisture, Kjeldahl protein, Soxhlet fat; ribose/glucose by HPLC-ELSD (recovery 96.3 / 99.6 %, not corrected); thiamine HPLC-DAD 267 nm (recovery 98.1 %); FAMEs (C19:0 IS); amino acids after **6 M HCl, 110 °C, 24 h** hydrolysis (total, not free) |

---

## 3. TABLES (verbatim)

### Table 1. Standards, SIM ions, calibration and SAFE recovery factors (RF)

Calibration: x = peak-area ratio to 1,2-dichlorobenzene, y = concentration (mg/g) in the SAFE
concentrate. Standard solutions in CH2Cl2, mg/L. "–" = not involved in CB.

| RI | compound | standards (mg/L) | ions (m/z) | calibration | R² | RF BYJ | RF CB |
|---|---|---|---|---|---|---|---|
| 683 | 3-hydroxy-2-butanone | 56.25, 225, 450, 900, 1800, 3600, 7200 | 88, 73, 55, 45 | y = 0.0451x − 0.0007 | 0.9989 | 0.71 | 0.76 |
| 696 | dimethyl disulfide | 0.14, 0.28, 1.09, 4.36, 8.72, 17.44, 34.88 | 94, 79, 61, 45 | y = 0.0009x − 0.00002 | 0.9996 | 0.76 | 0.73 |
| 805 | hexanal | 3, 6, 25, 102, 205, 410, 820 | 100, 82, 72, 56 | y = 0.0535x − 0.0046 | 0.9980 | 0.89 | 0.85 |
| 817 | 2-methylpyrazine | 0.12, 0.23, 0.89, 3.54, 7.08, 14.16, 28.32 | 94, 79, 67, 53 | y = 0.0216x − 0.0002 | 0.9993 | 0.84 | 0.87 |
| 827 | 4-methylthiazole | 0.91, 1.82, 7.25, 29, 58, 116, 232 | 99, 84, 71 | y = 0.0359x − 0.0039 | 0.9990 | 0.81 | – |
| 829 | furfural | 1.55, 3.1, 6.2, 12.4, 24.8, 49.6 | 96, 95, 67 | y = 0.013x − 0.0001 | 0.9995 | – | 0.83 |
| 852 | 1-pentanol | 0.74, 1.47, 5.85, 23.4, 46.8, 93.6, 187.2 | 80, 70, 58 | y = 0.1832x − 0.0018 | 0.9993 | 0.82 | – |
| 865 | furfuryl alcohol | 0.37, 0.74, 2.95, 11.8, 23.6, 47.2, 94.4 | 98, 81, 70, 69 | y = 0.0482x − 0.0008 | 0.9992 | 0.91 | – |
| 876 | **2-methyl-3-furanthiol** | 0.15, 0.30, 1.19, 4.76, 9.52, 19.04, 38.08 | 114, 113, 85, 71 | y = 0.0013x + 0.0002 | 0.9994 | **0.94** | **0.89** |
| 898 | heptanal | 0.46, 0.93, 1.85, 3.70, 7.40, 14.80, 29.60 | 114, 96, 86, 70 | y = 0.017x − 0.00009 | 0.9892 | 0.85 | 0.91 |
| 907 | **3-(methylthio)propanal** | 11.70, 46.80, 93.60, 187.20, 374.40, 748.80 | 104, 76, 69, 61 | y = 0.4401x − 0.0006 | 0.9985 | 0.79 | 0.82 |
| 911 | **2-furfurylthiol** | 0.08, 0.31, 1.23, 2.46, 4.92, 9.84, 19.68 | 114, 96, 85, 81 | y = 0.0027x + 0.0003 | 0.9945 | 0.91 | – |
| 922 | benzaldehyde | 0.16, 0.31, 1.23, 4.90, 9.80, 19.60, 39.20 | 106, 84, 77, 73 | y = 0.0021x − 0.0002 | 0.9997 | 0.82 | 0.78 |
| 953 | (E)-2-heptenal | 0.72, 1.43, 5.69, 22.75, 45.50, 91.00, 182.00 | 112, 97, 83, 70 | y = 0.0150x − 0.0012 | 0.9926 | 0.83 | – |
| 980 | 1-octen-3-ol | 4.69, 18.75, 37.50, 75.00, 150.00, 300.00 | 127, 110, 99, 85 | y = 0.1097x − 0.0024 | 0.9950 | 0.75 | 0.72 |
| 988 | 2-pentylfuran | 0.45, 0.89, 1.78, 3.56, 7.13, 14.25, 28.5 | 138, 109, 95, 81 | y = 0.0032x − 0.00005 | 0.9966 | 0.87 | 0.93 |
| 997 | 2-methyl-3-(methylthio)furan | 0.04, 0.13, 0.25, 0.50, 0.99, 1.98, 3.96 | 128, 113, 99, 85 | y = 0.0010x + 0.00004 | 0.9987 | 0.94 | – |
| 1004 | octanal | 1.28, 2.55, 5.10, 10.20, 20.40, 40.80, 81.60 | 128, 110, 100, 84 | y = 0.0588x − 0.0015 | 0.9980 | 0.91 | 0.92 |
| 1019 | 2-acetylthiazole | 0.30, 0.60, 1.19, 2.38, 4.76, 9.53, 19.05 | 127, 112, 99, 85 | y = 0.0057x + 0.0007 | 0.9969 | 0.84 | 0.81 |
| 1039 | phenylacetaldehyde | 0.13, 0.25, 0.97, 3.87, 7.74, 15.48, 30.96 | 120, 102, 91, 85 | y = 0.0019x + 0.0002 | 0.9995 | 0.94 | 0.89 |
| 1051 | (E)-2-octenal | 0.70, 1.40, 5.57, 22.25, 44.50, 89.00, 178.00 | 125, 108, 97, 83 | y = 0.0161x + 0.00002 | 0.9930 | 0.97 | 0.95 |
| 1056 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 0.79, 3.16, 6.32, 25.25, 50.50, 101.00, 202.00 | 128, 85, 72, 67 | y = 0.2939x + 0.0031 | 0.9993 | 0.82 | 0.77 |
| 1062 | 1-octanol | 0.15, 0.29, 1.13, 4.50, 9.00, 18.00, 36.00 | 129, 112, 97, 84 | y = 0.0038x + 0.0003 | 0.9983 | 0.78 | 0.87 |
| 1077 | 3-octen-2-one | 0.39, 0.78, 3.12, 12.45, 24.90, 49.80, 99.60 | 126, 111, 97, 83 | y = 0.0421x − 0.0018 | 0.9947 | 0.86 | – |
| 1082 | 2-ethyl-3,5-dimethylpyrazine | 0.12, 0.48, 0.96, 1.92, 3.84, 7.68, 15.36 | 136, 135, 121, 108 | y = 0.0014x + 0.00004 | 0.9993 | 0.95 | 0.91 |
| 1098 | nonanal | 1.16, 2.32, 9.25, 37, 74, 148, 296 | 142, 124, 114, 98 | y = 0.1299x − 0.0005 | 0.9928 | 0.94 | 0.92 |
| 1145 | (E,E)-2,4-octadienal | 0.04, 0.07, 0.26, 1.02, 2.04, 4.08, 8.16 | 124, 95, 91, 81 | y = 0.0271x − 0.0006 | 0.9941 | 0.91 | 0.89 |
| 1152 | (E,Z)-2,6-nonadienal | 0.08, 0.15, 0.59, 2.33, 4.65, 9.30, 18.60 | 138, 120, 109, 94 | y = 0.0010x − 0.00006 | 0.9981 | 0.95 | 0.9 |
| 1159 | (E)-2-nonenal | 0.13, 0.26, 1.02, 4.05, 8.10, 16.20, 32.40 | 139, 111, 96, 83 | y = 0.0429x − 0.0004 | 0.9957 | 0.92 | 0.91 |
| 1176 | 2-pentylthiophene | 0.02, 0.04, 0.13, 0.50, 0.99, 1.98, 3.96, 7.92 | 154, 97, 134, 111 | y = 0.0007x − 0.000003 | 0.9997 | 0.84 | – |
| 1192 | 2-pentylpyridine | 0.02, 0.06, 0.11, 0.43, 0.85, 1.70, 3.40 | 148, 134, 120, 93 | y = 0.0005x + 0.00005 | 0.9986 | 0.95 | – |
| 1205 | decanal | 0.10, 0.39, 0.77, 1.53, 3.05, 6.10, 12.20 | 156, 128, 112, 95 | y = 0.0012x + 0.00002 | 0.9973 | 0.85 | 0.87 |
| 1212 | (E,E)-2,4-nonadienal | 0.68, 1.35, 5.38, 21.50, 43.00, 86.00, 172.00 | 138, 109, 95, 81 | y = 0.0189x − 0.0007 | 0.9937 | 0.93 | 0.96 |
| 1247 | (E)-2-decenal | 0.38, 1.50, 3.00, 12, 24, 48, 96 | 154, 136, 121, 110 | y = 0.0227x − 0.00008 | 0.9998 | 0.97 | 0.93 |
| 1275 | 4-methyl-5-thiazoleethanol | 0.04, 0.13, 0.26, 1.02, 2.03, 4.06, 8.12 | 143, 125, 112, 98 | y = 0.0424x − 0.0092 | 0.9969 | 0.94 | – |
| 1292 | 2-undecanone | 0.10, 0.20, 0.78, 3.12, 6.24, 12.48, 24.96 | 170, 155, 127, 112 | y = 0.0023x + 0.000006 | 0.9988 | 0.82 | – |
| 1307 | undecanal | 0.41, 0.81, 3.22, 12.85, 25.7, 51.4, 102.8 | 170, 126, 110, 96 | y = 0.0141x − 0.0005 | 0.9974 | 0.94 | – |
| 1314 | (E,E)-2,4-decadienal | 1.47, 2.94, 11.75, 47, 94, 188, 376 | 152, 123, 109, 95 | y = 0.1762x − 0.0045 | 0.9962 | 0.95 | 0.97 |
| 1362 | (E)-2-undecenal | 0.8, 1.6, 6.4, 25.6, 51.2, 102.4, 204.8 | 168, 126, 111, 97 | y = 0.2406x − 0.0025 | 0.9915 | 0.89 | 0.91 |
| 1378 | γ-nonanolactone | 0.80, 1.45, 5.90, 23.3, 46.5, 93, 186 | 155, 128, 114, 100 | y = 0.0016x + 0.00009 | 0.9992 | 0.87 | 0.83 |
| 1460 | 1-dodecanol | 0.18, 0.35, 1.37, 5.45, 10.9, 21.8, 43.6 | 168, 140, 125, 111 | y = 0.0342x − 0.0008 | 0.9917 | 0.81 | – |
| 1491 | β-ionone | 0.02, 0.04, 0.14, 0.53, 1.05, 2.1, 4.2 | 192, 177, 135, 121 | y = 0.0008x + 0.00007 | 0.9973 | 0.83 | 0.85 |
| 1540 | **bis(2-methyl-3-furyl)disulfide** | 0.21, 0.41, 1.63, 6.5, 13, 26, 52 | 258, 226, 183, 162 | y = 0.0106x − 0.0009 | 0.9956 | 0.87 | – |
| 1597 | 1-tridecanol | 0.24, 0.93, 1.85, 7.4, 14.8, 29.6, 59.2 | 182, 139, 125, 111 | y = 0.0618x + 0.000006 | 0.9985 | 0.86 | 0.88 |
| 1002 | 1,2-dichlorobenzene (IS) | 25 | 150, 146, 131, 111 | internal standard | | 0.95 | 0.95 |

### Table 2. Breast muscle composition (mean ± SD, n = 6; letters = Duncan p < 0.05)

| item | BJY | CB |
|---|---|---|
| moisture (%) | 73.16 ± 0.35 a | 74.44 ± 0.80 a |
| crude fat (%) | 1.07 ± 0.02 a | 0.91 ± 0.01 b |
| crude protein (%) | 22.33 ± 0.21 b | 23.53 ± 0.16 a |
| **ribose (mg/100 g)** | **1.37 ± 0.05 a** | **1.09 ± 0.04 b** |
| glucose (mg/100 g) | 141.73 ± 2.26 a | 120.55 ± 2.91 b |
| **thiamine (mg/100 g)** | **0.19 ± 0.02 a** | **0.12 ± 0.01 b** |

### Table 3. Breast muscle fatty acids (% of total FA) and amino acids (mg/g meat, acid hydrolysate), mean ± SD, n = 6

| fatty acid | BJY | CB | | amino acid | BYJ | CB |
|---|---|---|---|---|---|---|
| SFA | 32.99 ± 0.34 b | 34.48 ± 0.19 a | | Aspartic acid | 18.81 ± 0.02 b | 20.09 ± 0.77 a |
| C14:0 | 0.43 ± 0.02 a | 0.49 ± 0.01 a | | Threonine | 8.55 ± 0.14 b | 9.66 ± 0.57 a |
| C16:0 | 20.65 ± 0.25 b | 21.77 ± 0.13 a | | Serine | 8.01 ± 0.45 b | 9.03 ± 0.37 a |
| C18:0 | 11.25 ± 0.04 a | 10.96 ± 0.02 b | | Glutamic acid | 30.09 ± 0.13 b | 32.18 ± 0.19 a |
| C20:0 | 0.17 ± 0.01 a | 0.16 ± 0.01 a | | Glycine | 8.69 ± 0.72 b | 9.92 ± 0.43 a |
| C22:0 | 0.19 ± 0.01 a | 0.16 ± 0.01 b | | Alanine | 12.13 ± 0.23 b | 12.88 ± 0.71 a |
| C24:0 | 0.30 ± 0.01 b | 0.94 ± 0.01 a | | **Cysteine** | **4.26 ± 0.11 a** | **3.56 ± 0.12 b** |
| MUFA | 28.42 ± 0.26 b | 31.80 ± 0.52 a | | Valine | 12.93 ± 1.02 b | 13.76 ± 1.12 a |
| C16:1 | 1.59 ± 0.04 b | 2.09 ± 0.12 a | | **Methionine** | **5.69 ± 0.34 a** | **5.41 ± 0.38 b** |
| C18:1 | 26.21 ± 0.18 b | 28.93 ± 0.36 a | | Isoleucine | 9.80 ± 0.11 b | 11.11 ± 0.17 a |
| C20:1 | 0.23 ± 0.02 b | 0.34 ± 0.02 a | | Leucine | 19.79 ± 0.90 b | 21.37 ± 1.25 a |
| C24:1 | 0.39 ± 0.02 a | 0.44 ± 0.02 a | | Tyrptophan [sic] | 6.57 ± 0.07 a | 6.86 ± 0.13 a |
| PUFA | 39.80 ± 0.39 a | 34.88 ± 0.60 b | | Phenylalanine | 9.86 ± 0.71 b | 11.00 ± 0.65 a |
| C18:2 | 27.52 ± 0.23 a | 24.06 ± 0.18 b | | Histidine | 7.61 ± 0.34 a | 7.81 ± 0.75 a |
| C18:3 | 0.95 ± 0.01 b | 1.25 ± 0.01 a | | Lysine | 13.24 ± 0.21 b | 14.16 ± 0.03 a |
| C20:2 | 0.65 ± 0.01 b | 0.72 ± 0.02 a | | Arginine | 13.22 ± 0.48 b | 13.90 ± 0.69 a |
| C20:3 | 1.18 ± 0.05 b | 4.35 ± 0.13 a | | Proline | 7.46 ± 0.02 b | 7.87 ± 0.28 a |
| C20:4 | 4.60 ± 0.02 a | 1.96 ± 0.12 b | | Total | 196.71 ± 6.00 b | 210.57 ± 10.61 a |
| C20:5 | 0.35 ± 0.01 a | 0.13 ± 0.01 b | | | | |
| C22:4 | 2.32 ± 0.03 a | 1.09 ± 0.09 b | | | | |
| C22:5 | 0.81 ± 0.01 a | 0.59 ± 0.01 b | | | | |
| C22:6 | 1.42 ± 0.02 a | 0.73 ± 0.03 b | | | | |

### Table 4. AEDA/GC-O odorants and quantified concentrations (DB-5MS RI; ng/g broth, mean ± SD of n = 3 pooled concentrates; "–" in log2FD = not detected by GC-O; "–" in concentration = not selected for quantification)

| RI | compound | odour | log2FD BYJ | log2FD CB | BYJ (ng/g) | CB (ng/g) | ID |
|---|---|---|---|---|---|---|---|
| 549 | unknown | spicy, sweet | 1 | – | – | – | – |
| 656 | methanethiol | sulfur, gasoline | 8 | 4 | – | – | RI/Odor (t) |
| 683 | 3-hydroxy-2-butanone | buttery | 8 | 4 | 1592.27 ± 1.28 a | 1078.49 ± 1.42 b | S/RI/Odor/MS |
| 696 | dimethyl disulfide | onion, cabbage | 10 | 3 | 0.75 ± 0.03 a | 0.62 ± 0.02 b | S/RI/Odor/MS |
| 716 | 1-penten-3-ol | pungent, butter | 5 | – | – | – | S/RI/Odor/MS |
| 727 | pentanal | grass, pungent | 5 | – | – | – | S/RI/Odor/MS |
| 754 | thiazole | meaty | 2 | – | – | – | S/RI/Odor |
| 768 | dihydro-2-methyl-3(2H)-furanone | caramel | 5 | – | – | – | RI/Odor (t) |
| 777 | 2-methylthiophene | sulfury | 0 | – | – | – | S/RI/Odor/MS |
| 805 | **hexanal** | grass, green | 10 | 10 | **107.39 ± 2.42 b** | **141.16 ± 3.32 a** | S/RI/Odor/MS |
| 817 | 2-methylpyrazine | roasted, nutty | 7 | 11 | 1.96 ± 0.10 b | 2.29 ± 0.14 a | S/RI/Odor/MS |
| 827 | 4-methylthiazole | roasted meat | 9 | – | 3.91 ± 0.20 | – | S/RI/Odor |
| 829 | **furfural** | caramel | 3 | 3 | – | **4.11 ± 0.32** | S/RI/Odor/MS |
| 852 | 1-pentanol | green, fruity | 7 | – | 46.52 ± 1.87 | – | S/RI/Odor/MS |
| 865 | furfuryl alcohol | caramel | 12 | 0 | 13.03 ± 0.29 | – | S/RI/Odor/MS |
| 876 | **2-methyl-3-furanthiol** | meaty | 11 | 12 | **1.07 ± 0.09 b** | **1.48 ± 0.18 a** | S/RI/Odor |
| 898 | heptanal | fatty, green | 7 | 5 | 6.55 ± 0.42 a | 4.89 ± 0.42 b | S/RI/Odor/MS |
| 907 | **3-(methylthio)propanal** | cooked potato | 12 | 10 | **158.15 ± 2.98 a** | **117.99 ± 0.25 b** | S/RI/Odor |
| 911 | **2-furfurylthiol** | coffee, roasted | 12 | 0 | **2.71 ± 0.22** | – | S/RI/Odor |
| 913 | 2,5-dimethylpyrazine | roasted, popcorn | 6 | – | – | – | S/RI/Odor/MS |
| 922 | benzaldehyde | almond | 8 | 3 | 5.50 ± 0.21 a | 2.35 ± 0.12 b | S/RI/Odor/MS |
| 953 | (E)-2-heptenal | green, fatty | 10 | – | 16.11 ± 0.26 | – | S/RI/Odor/MS |
| 980 | 1-octen-3-ol | mushroom | 11 | 6 | 107.69 ± 0.19 a | 92.13 ± 1.30 b | S/RI/Odor/MS |
| 984 | 3-octanol | herbal, mushroom | 0 | – | – | – | S/RI/Odor/MS |
| 988 | 2-pentylfuran | earthy, green | 9 | 5 | 4.99 ± 0.16 a | 3.54 ± 0.12 b | S/RI/Odor/MS |
| 997 | 2-methyl-3-(methylthio)furan | meaty | 10 | – | 0.51 ± 0.04 | – | S/RI/Odor/MS |
| 1004 | octanal | citrus, green | 10 | 3 | 10.02 ± 0.88 a | 3.68 ± 0.36 b | S/RI/Odor/MS |
| 1011 | (E,E)-2,4-heptadienal | green, fatty | 3 | 1 | – | – | S/RI/Odor/MS |
| 1019 | **2-acetylthiazole** | roasted, sesame | 7 | 3 | **6.32 ± 0.41 a** | **3.69 ± 0.18 b** | S/RI/Odor/MS |
| 1033 | limonene | citrus | 3 | 0 | – | – | S/RI/Odor/MS |
| 1039 | phenylacetaldehyde | sweet, flowery | 10 | 8 | 1.35 ± 0.07 a | 1.09 ± 0.09 b | S/RI/Odor/MS |
| 1045 | 2-acetylpyrrole | licorice, nutty | 4 | 2 | – | – | S/RI/Odor/MS |
| 1051 | (E)-2-octenal | green, fatty | 12 | 12 | 13.97 ± 0.07 a | 14.47 ± 0.11 a | S/RI/Odor/MS |
| 1056 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | caramel | 7 | 12 | 6.82 ± 0.35 b | 13.65 ± 0.87 a | S/RI/Odor/MS |
| 1062 | 1-octanol | fatty, citrus | 7 | 11 | 1.28 ± 0.09 b | 3.87 ± 0.32 a | S/RI/Odor/MS |
| 1077 | 3-octen-2-one | mushroom, earthy | 9 | – | 4.81 ± 0.30 | – | S/RI/Odor/MS |
| 1082 | 2-ethyl-3,5-dimethylpyrazine | earthy, roasted | 11 | 12 | 1.63 ± 0.02 b | 2.32 ± 0.19 a | S/RI/Odor |
| 1098 | nonanal | green, fatty | 9 | 10 | 25.45 ± 3.26 a | 26.99 ± 1.05 a | S/RI/Odor/MS |
| 1145 | (E,E)-2,4-octadienal | green, fatty | 12 | 9 | 2.59 ± 0.14 a | 2.13 ± 0.05 b | S/RI/Odor/MS |
| 1152 | (E,Z)-2,6-nonadienal | cucumber | 7 | 4 | 0.47 ± 0.02 a | 0.37 ± 0.02 b | S/RI/Odor/MS |
| 1159 | (E)-2-nonenal | green, fatty | 8 | 11 | 2.9 ± 0.19 b | 7.62 ± 0.75 a | S/RI/Odor/MS |
| 1176 | 2-pentylthiophene | fruity | 7 | – | 1.18 ± 0.05 | – | S/RI/Odor/MS |
| 1192 | 2-pentylpyridine | burnt, fatty | 7 | – | 0.22 ± 0.03 | – | S/RI/Odor/MS |
| 1205 | decanal | green, fatty | 8 | 11 | 1.94 ± 0.10 b | 2.63 ± 0.20 a | S/RI/Odor/MS |
| 1212 | (E,E)-2,4-nonadienal | fatty, green | 10 | 9 | 7.52 ± 0.27 a | 5.72 ± 0.52 b | S/RI/Odor/MS |
| 1223 | 2,3,5-trimethylpyrazine | earthy, musty | 3 | – | – | – | S/RI/Odor/MS |
| 1230 | α-terpineol | herbal, floral | 7 | 4 | – | – | RI/Odor/MS (t) |
| 1247 | (E)-2-decenal | fatty, tallowy | 8 | 12 | 11.33 ± 0.25 b | 21.58 ± 1.58 a | S/RI/Odor/MS |
| 1250 | unknown | pine, green | 2 | – | – | – | – |
| 1261 | phenylethyl alcohol | rose | 3 | – | – | – | S/RI/Odor/MS |
| 1275 | 4-methyl-5-thiazoleethanol | meaty | 7 | – | 0.18 ± 0.01 | – | S/RI/Odor/MS |
| 1278 | (E,Z)-2,4-decadienal | fatty | 5 | 4 | – | – | RI/Odor/MS (t) |
| 1281 | unknown | herbal | 0 | – | – | – | – |
| 1292 | 2-undecanone | fruity, fatty | 8 | – | 1.09 ± 0.03 | – | S/RI/Odor/MS |
| 1302 | 1-nonen-3-ol | green, fatty | 2 | – | – | – | S/RI/Odor/MS |
| 1307 | undecanal | fatty, tallowy | 9 | – | 4.93 ± 0.17 | – | S/RI/Odor/MS |
| 1314 | **(E,E)-2,4-decadienal** | deep-fried, fatty | 12 | 11 | **104.22 ± 2.82 a** | **71.37 ± 5.16 b** | S/RI/Odor/MS |
| 1320 | unknown | earthy | 2 | – | – | – | – |
| 1341 | hexyl hexanoate | sweet, fruity | 4 | – | – | – | S/RI/Odor/MS |
| 1362 | (E)-2-undecenal | fatty, metallic | 11 | 5 | 34.97 ± 0.11 a | 14.78 ± 1.15 b | S/RI/Odor/MS |
| 1378 | γ-nonanolactone | coconut | 11 | 6 | 3.89 ± 0.12 a | 3.25 ± 0.19 b | S/RI/Odor/MS |
| 1391 | naphthalene | camphor | 4 | – | – | – | RI/Odor/MS (t) |
| 1405 | vanillin | vanilla-like | 2 | – | – | – | S/RI/Odor/MS |
| 1450 | (E)-2-dodecenal | fatty, pungent | 5 | – | – | – | S/RI/Odor/MS |
| 1460 | 1-dodecanol | wax, sweet | 7 | – | 1.22 ± 0.28 | – | S/RI/Odor/MS |
| 1467 | γ-decalactone | peach | 4 | – | – | – | S/RI/Odor/MS |
| 1491 | β-ionone | violet | 7 | 5 | 0.13 ± 0.01 a | 0.08 ± 0.01 b | S/RI/Odor/MS |
| 1514 | unknown | spicy | 1 | – | – | – | – |
| 1520 | tridecanal | fatty, sweet | 3 | 4 | – | – | RI/Odor/MS (t) |
| 1540 | **bis(2-methyl-3-furyl)disulfide** | meat broth | 7 | – | **1.62 ± 0.02** | – | S/RI/Odor/MS |
| 1560 | (E,E)-2,4-dodecadienal | fatty, roasted | 2 | – | – | – | RI/Odor/MS (t) |
| 1589 | tetradecanal | fatty, floral | 2 | – | – | – | RI/Odor/MS (t) |
| 1597 | 1-tridecanol | tallowy | 7 | 7 | 0.76 ± 0.04 b | 1.74 ± 0.14 a | S/RI/Odor/MS |
| 1626 | unknown | musty | – | 2 | – | – | – |
| 1696 | 2-pentadecanone | fatty, sweet | 4 | – | – | – | RI/Odor/MS (t) |
| 1795 | 4-hydroxy-5-ethyl-2-methyl-3(2H)-furanone | caramel | 1 | – | – | – | RI/Odor/MS (t) |
| 1864 | hexadecanal | fatty | 3 | – | – | – | RI/Odor/MS (t) |

ID: S = all parameters agree with the authentic standard; RI; MS (NIST 15); Odor; (t) tentative.

### Table 5. Recombination and omission: odorants with OAV ≥ 1 (OAV = conc / threshold in water, van Gemert 2011; impact: ** p ≤ 0.01, * p ≤ 0.05, – n.s.)

| compound | OAV BYJ | OAV CB | threshold (µg/kg) | impact | [Z] OAV recomputed BYJ / CB |
|---|---|---|---|---|---|
| 3-hydroxy-2-butanone | 2 | 1 | 800 | ** | 2.0 / 1.3 |
| dimethyl disulfide | 12 | 9 | 0.06 | ** | 12.5 / 10.3 |
| hexanal | 23 | 30 | 4.5 | ** | 23.9 / 31.4 |
| **2-methyl-3-furanthiol** | **2686** | **14522** | **0.0004** | ** | 2675 / **3700** [!] |
| heptanal | 2 | 1 | 3 | ** | 2.2 / 1.6 |
| 3-(methylthio)propanal | 3953 | 2949 | 0.04 | ** | 3954 / 2950 |
| 2-furfurylthiol | 543 | – | 0.005 | ** | 542 / – |
| (E)-2-heptenal | 1 | – | 13 | ** | 1.2 / – |
| 1-octen-3-ol | 53 | 46 | 2 | ** | 53.8 / 46.1 |
| 2-methyl-3-(methylthio)furan | 102 | – | 0.005 | ** | 102 / – |
| octanal | 14 | 5 | 0.7 | ** | 14.3 / 5.3 |
| (E)-2-octenal | 5 | 5 | 3 | * | 4.7 / 4.8 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 1 | 3 | 5 | ** | 1.4 / 2.7 |
| 2-ethyl-3,5-dimethylpyrazine | 1 | 2 | 1 | ** | 1.6 / 2.3 |
| nonanal | 25 | 26 | 1 | ** | 25.5 / 27.0 |
| (E,Z)-2,6-nonadienal | 23 | 18 | 0.02 | ** | 23.5 / 18.5 |
| (E)-2-nonenal | 36 | 95 | 0.08 | * | 36 / 95 |
| decanal | 19 | 2 | 0.1 | – | 19.4 / **26.3** [!] |
| (E,E)-2,4-nonadienal | 125 | 95 | 0.06 | ** | 125 / 95 |
| (E,E)-2,4-decadienal | 2084 | 1427 | 0.05 | ** | 2084 / 1427 |
| (E)-2-undecenal | 11 | 4 | 3 | ** | 11.7 / 4.9 |
| bis(2-methyl-3-furyl)disulfide | 810 | – | 0.002 | ** | 810 / – |

[Z] molar view of the sulfur targets (MW: MFT/FFT 114.17; MFT-SS 226.3; methional 104.17;
hexanal 100.16; 2-acetylthiazole 127.16):

| compound | BYJ nmol/kg | CB nmol/kg |
|---|---|---|
| MFT | 9.4 | 13.0 |
| MFT disulfide (as MFT equivalents) | 7.2 (14.3) | – |
| MFT + 2 × MFT-SS (total MFT-derived thiol) | 23.7 | ≥ 13.0 |
| FFT | 23.7 | – |
| 2-methyl-3-(methylthio)furan | 4.0 | – |
| methional | 1518 | 1133 |
| 2-acetylthiazole | 49.7 | 29.0 |
| hexanal | 1072 | 1409 |

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor (external real-food target levels, not benchmark pots):**
- **Chicken broth MFT ≈ 1.1-1.5 µg/kg (9-13 nmol/kg); FFT ≈ 2.7 µg/kg (BYJ); methional ≈ 120-160
  µg/kg; hexanal ≈ 110-140 µg/kg; 2-acetylthiazole ≈ 4-6 µg/kg.** With Guth & Grosch 1994 (beef juice:
  FFT 0.5, methional 13, hexanal 72 µg/kg, MFT not found) these bracket a "cooked meat" target band
  for an analogue: **FFT 0.5-3 µg/kg, methional 13-160 µg/kg, hexanal 70-140 µg/kg, MFT ~1-1.5 µg/kg
  (chicken only)**. Factor-of-~10 spread on methional between the two meats/methods.
- **MFT : FFT ≈ 0.4 (BYJ, mass)**, and MFT disulfide holds ~60 % of MFT-derived sulfur in BYJ
  ([Z] 14.3 / 23.7 nmol/kg as MFT equivalents). If the engine predicts MFT and MFT-SS for a broth-
  like pot, the disulfide fraction is a within-food ratio worth checking (with the caveat that SAFE
  work-up can oxidise MFT).
- **Precursor inventory of the raw meat** (breast): ribose 1.09-1.37 mg/100 g (≈ 73-91 µmol/kg),
  glucose 121-142 mg/100 g (≈ 6.7-7.9 mmol/kg), thiamine 0.12-0.19 mg/100 g (≈ 4.5-7.2 µmol/kg at
  265.4 g/mol; [!] salt form not stated), cysteine 3.56-4.26 mg/g **total (acid
  hydrolysate), not free**. These are the right order of inputs for a "chicken meat" initial
  condition, except cysteine, which must not be read as free cysteine.

**Cannot anchor:**
- **No kinetics.** One cook (3 h, ~100 °C), no time series, no pH, no water loss, no broth mass.
- **The pot is not a closed mass-action system:** a whole carcass (skin, fat, bone) in 1.5 L water,
  crock pot, open or lidded (not stated); the meat analyses are breast only. Broth concentrations
  cannot be converted to yields per precursor.
- **Quantification is not SIDA.** One internal standard (1,2-dichlorobenzene) with external
  calibration curves in CH2Cl2 and recovery factors measured by spiking odourless broth residue.
  For reactive thiols (MFT, FFT) the recovery spike was at ~native level and RFs are 0.89-0.94, but
  the method does not track losses during the 3 h cook or any thiol-disulfide exchange during SAFE.
  Treat the thiol numbers as ±2× (my judgement, ~70 % confidence that the true broth MFT lies in
  0.5-3 µg/kg).
- **BYJ vs CB is not a precursor experiment.** Breeds differ in age (120 vs 42 d), sex (hens vs
  unstated), fat, PUFA, ribose, thiamine and cysteine simultaneously; and CB has **more** MFT despite
  less ribose, thiamine and cysteine. No causal reading.

**Flags:**
1. [!] **CB MFT OAV 14522 is inconsistent with its own concentration.** 1.48 ng/g / 0.0004 µg/kg =
   3700. 14522 × 0.0004 = 5.8 ng/g, which is not printed anywhere. Either the CB MFT concentration in
   Table 4 or the OAV in Table 5 is wrong. Use the concentration (1.48), which is the measured
   quantity; the Abstract and Discussion do not quote either number.
2. [!] **CB decanal OAV printed 2; recomputed 26** (2.63 / 0.1). Likely a typo; immaterial for the
   engine.
3. [!] CB FFT is log2FD 0 (FD 1) and not quantified, yet the CB recombinate is said to reproduce the
   broth; FFT is simply absent from the CB model.
4. [!] The MFT threshold used (0.0004 µg/kg in water) is at the low end of the literature range;
   OAVs scale with it. The concentrations do not depend on it.
5. [!] Tables 1 and 4 say "BYJ" and the text says "BJY"; same breed.
6. The PDF is the Just Accepted manuscript; the typeset version (66, 10242-10251) may differ in
   minor edits. Numbers here are from the manuscript.

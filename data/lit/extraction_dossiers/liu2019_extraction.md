# Liu et al. 2019 — EXTRACTION
### *"Characterization of Key Aroma Compounds in Beijing Roasted Duck by Gas Chromatography-Olfactometry-Mass Spectrometry, Odor Activity Values and Aroma Recombination Experiments"*
### Restaurant Beijing roast ducks (4 brands), hung-oven ("Gua Lu") roasting **60 min at 250 °C** (Dadong 80 min); breast skin and breast muscle minced. HS-SPME (CAR/PDMS, 55 °C) + GC-MS-SIM against one internal standard (2-methyl-3-heptanone) with **matrix-matched calibration curves in solvent-stripped duck skin**, ng/g. **FFT 18–158, methional 27–215, DMTS 12.6–117, hexanal 11–1116 ng/g. MFT and DMDS not reported.** Read 2026-10-28.

**Full read of `data/articles/liu2019.pdf`** (37 pp.: the ACS "Just Accepted" manuscript, not the typeset
version; tables pp. 27-33). Text layer clean; Tables 1-5 transcribed from `pdftotext -layout`, and every
OAV in Table 4 recomputed for the sulfur compounds and hexanal (all consistent). Figures 1 (heat map) and 2
(spider plot, recombinate vs duck) not digitised (no numbers needed).

**Provenance:** **[M]** printed number · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.
1 ng/g = 1 µg/kg.

---

## 0. IDENTITY — matches

| field | value |
|---|---|
| Authors | Huan Liu, Zhenyu Wang (equal), Dequan Zhang*, Qingwu Shen, Teng Pan, Teng Hui, Jianrong Ma |
| Title | as above (first page matches the requested paper) |
| Venue | *J. Agric. Food Chem.* — PDF is the Just Accepted manuscript ("Publication Date (Web): 02 May 2019"). Crossref 2026-10-28: **67 (20), 5847-5856** |
| DOI | **10.1021/acs.jafc.9b01564** (printed on p. 1) |
| Affiliation | Institute of Food Science and Technology, CAAS, Beijing; Hunan Agricultural University |
| Companions | Chen, Song & Ma 2010 (Flavour Fragr. J., aroma-active compounds of Beijing roast duck; not on disk). Cites Fan 2018 (`fan2018_extraction.md`) and Zhao 2017 (`zhao2017_extraction.md`) |

---

## 1. THE ONE-PARAGRAPH ANSWER

In Beijing roast duck the three sulfur odorants quantified are **2-furfurylthiol 18-158 ng/g, methional
27-215 ng/g and dimethyl trisulfide 12.6-117 ng/g**; hexanal spans **11-1116 ng/g** (breast skin low,
breast muscle high; Quanjude muscle an outlier at 1116). **2-Methyl-3-furanthiol, its disulfide and
dimethyl disulfide do not appear anywhere in the paper** (not in the 42 identified compounds). On water
thresholds of 0.005 µg/kg FFT and DMTS dominate the OAV sum (FFT 44-68 %, DMTS 27-41 %). The FFT levels are
**7-60× Fan 2018's chicken broth and 36-316× Guth & Grosch's beef juice**; the method (HS-SPME of minced
meat at 55 °C, calibration in a solvent-stripped, freeze-dried duck-skin matrix with one non-isotopic
internal standard) is the weakest quantification class in this set, and the printed calibration table has
copy errors (sec. 4). Treat these as order-of-magnitude levels.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Meat | Beijing duck (GI-protected breed), 38 d, 1.40 ± 0.09 kg; 4 restaurants (Quanjude, Bianyifang, Dongxinglou, Dadong); 3 batches each, 2 ducks drawn per batch of 20 -> **6 ducks per brand** |
| Cooking | "All Beijing roasted ducks were heated for **60 min at 250 °C**, except the Dadong roasted duck which was heated for **80 min**." Discussion: Gua Lu roasting, duck hung in an oven at **200-250 °C**. Oven (air) temperature; no product temperature, no pre-treatment (scalding, maltose glaze, inflation) described [!] |
| Samples | breast skin and breast muscle, trimmed of connective tissue and visible fat, minced, vacuum film, liquid-N2 snap-frozen, −80 °C |
| Extraction | **HS-SPME**: 4.0 g mince in 20 mL vial + 1.5 µL 2-methyl-3-heptanone (1.68 µg/µL in methanol; [Z] 2.52 µg = 630 ng/g); 10 min incubation then 45 min extraction at **55 °C**, CAR/PDMS 75 µm; desorb 250 °C 3 min |
| GC | DB-WAX 30 m; 1:1 split to MS and ODP; EI 70 eV full scan for ID, **SIM for quantitation** |
| Screening quantitation | "semiquantitated based on 2-methyl-3-heptanone as an internal standard" (response factor 1 implied) |
| Final quantitation | compounds with OAV > 1 re-quantified "using calibration curves constructed from the odorless matrix of roasted duck": Quanjude breast **skin** stripped with diethyl ether/pentane (2:1:1 m:m:m), 8 h shaking ≥ 3×, freeze-dried 20 h at −60 °C; standards "1 mg/L to 1 g/L" with IS added to the matrix, extracted and run as samples. **Not stable-isotope dilution; one IS; matrix = skin, applied to muscle too** [!] |
| Replicates | n = 6 ducks per brand (mean ± SD) |
| Thresholds | in water, from literature and "flavor base" software (Table 4 footnotes) |
| Sensory | 9 panelists; recombinate in odourless duck skin + water (78 % water); omission by triangle test; similarity 4.5/5 |

---

## 3. TABLES (verbatim)

### Table 1. Identification of aroma compounds (LRI on DB-WAX; ID: MS, LRI, S = authentic standard, O = olfactometry)

| compound | LRI lit | LRI calc | ID | | compound | LRI lit | LRI calc | ID |
|---|---|---|---|---|---|---|---|---|
| pentanal | 964 | 958 | MS, LRI, S, O | | acetic acid | 1401 | 1397 | MS, LRI |
| hexanal | 1064 | 1060 | MS, LRI, S, O | | **methional** | 1426 | 1412 | MS, LRI, S, O |
| heptanal | 1163 | 1160 | MS, LRI, S, O | | furfural | 1421 | 1417 | MS, LRI |
| (E)-2-hexenal | 1196 | 1188 | MS, LRI, S, O | | 1-octene-3-ol | 1430 | 1425 | MS, LRI, S, O |
| 2-pentylfuran | 1215 | 1211 | MS, LRI, S | | 2-acetylfuran | 1461 | 1457 | MS, LRI |
| styrene | 1236 | 1226 | MS, LRI | | pyrrole | 1470 | 1468 | MS, LRI |
| 2-methylpyrazine | 1238 | 1231 | MS, LRI | | benzaldehyde | 1476 | 1473 | MS, LRI, S |
| 2,4,6-trimethyl-1,3,5-triazine | – (f) | 1249 | MS | | propanoic acid | 1490 | 1488 | MS, LRI |
| octanal | 1273 | 1265 | MS, LRI, S, O | | (E)-2-nonenal | 1509 | 1500 | MS, LRI, S, O |
| 1-pentanol | 1274 | 1275 | MS, LRI, S, O | | 5-methylfuran | 1522 | 1523 | MS, LRI |
| 2,5-dimethylpyrazine | 1281 | 1287 | MS, LRI, S | | butanoic acid | 1589 | 1577 | MS, LRI |
| (E)-2-heptenal | 1291 | 1292 | MS, LRI, S, O | | phenyl acetaldehyde | 1606 | 1589 | MS, LRI, S, O |
| 2-ethylpyrazine | 1302 | 1297 | MS, LRI | | (E)-2-decenal | 1616 | 1609 | MS, LRI, S, O |
| 2,3-octanedione | 1325 | 1301 | MS, LRI | | isovaleric acid | 1633 | 1621 | MS, LRI |
| hexenyl acetate | 1313 | 1308 | MS, LRI | | 2(5H)-furanone | 1702 | 1687 | MS, LRI, S, O |
| **dimethyl trisulfide** | 1354 | 1340 | MS, LRI, S, O | | 2-undecenal | 1755 | 1719 | MS, LRI |
| 2-ethyl-6-methylpyrazine | 1353 | 1353 | MS, LRI | | (E,E)-2,4-decadienal | 1763 | 1771 | MS, LRI, S, O |
| 2-nonanone | – | 1363 | MS, S | | hexanoic acid | 1802 | 1790 | MS, LRI |
| nonanal | 1369 | 1368 | MS, LRI, S, O | | guaiacol | 1814 | 1800 | MS, LRI, S, O |
| **2-furfurylthiol** | 1402 | 1395 | MS, LRI, S, O | | maltol | 1936 | 1905 | MS, LRI |
| (E)-2-octenal | 1396 | 1395 | MS, LRI, S, O | | phenol | 1956 | 1940 | MS, LRI |

(f) not identified (no literature LRI). [!] "5-methylfuran" as printed (not a valid name; probably 5-methylfurfural or 2-methylfuran).

### Table 2. Standards, SIM ions, calibration equations ("x is the peak area relative to that of the internal standard 2-methyl-3-heptanone, y is the concentration (mg/g) in the samples relative to that of the internal standard")

| compound | ions (m/z) | calibration equation (as printed) | R² | odour (GC-O) |
|---|---|---|---|---|
| pentanal | 58,57,56,55 | y = 0.00002x + 2.41 | 0.9994 | fermented, fruity, nutty |
| hexanal | 56,57,72,82 | **y = 0.000002 + 28.70** [!] | 0.9999 | grassy, green, fatty |
| heptanal | 70,55,57,71 | **y = 0.00004 + 2.10** [!] | 0.9999 | fatty, green |
| octanal | 84,56,57,55 | y = 0.00003x + 123.65 | 0.9965 | fatty, green |
| nonanal | 57,56,98,70 | **y = 0.00008x – 358.28** | **0.9923** | fatty, fresh |
| phenyl acetaldehyde | 91,120,92,65 | y = 0.0001x – 25.04 | 0.9951 | green, sweet |
| (E)-2-hexenal | 69,83,55,57 | y = 0.0001x + 4.30 | 0.9998 | unpleasant, fatty |
| (E)-2-heptenal | 83,55,57,56 | y = 0.00006x + 10.74 | 0.9985 | unpleasant, fatty |
| (E)-2-octenal | 70,55,83,57 | **y = 0.00008x – 358.25** | **0.9923** | unpleasant, fatty |
| (E)-2-nonenal | 70,83,55,69 | y = 0.00006x + 3.91 | 0.9999 | unpleasant, fatty |
| (E)-2-decenal | 70,55,83,57 | y = 0.00007x + 1.32 | 0.9988 | unpleasant, fatty |
| (E,E)-2,4-decadienal | 81,152,67,83 | y = 0.00005x + 174.49 | 0.9993 | fatty, meat |
| **dimethyl trisulfide** | 126,79,128,80 | y = 0.00003x – 6.47 | 0.9998 | sulfurous, cooked, meaty |
| **2-furfurylthiol** | 81,114,53,82 | y = 0.00002x – 5.76 | 0.9984 | sulfury, roasted, smoky |
| **methional** | 104,76,61,57 | y = 0.00009x – 3.93 | 0.9913 | potato, creamy |
| 1-octene-3-ol | 57,72,85,55 | y = 0.00003x –19.88 | 0.9937 | green, mushroom |
| 2-pentylfuran | 81,138,82,53 | y = 0.000002x + 3.64 | 0.9968 | vegetable, fruity, metallic |
| guaiacol | 109,124,81,53 | y = 0.00002x – 6.88 | 0.9961 | smoke, vanilla, woody |

### Table 3. Quantitation (ng/g; mean ± SD, n = 6)

| compound | Dongxinglou skin | Dongxinglou muscle | Quanjude skin | Quanjude muscle | Bianyifang skin | Bianyifang muscle | Dadong skin | Dadong muscle |
|---|---|---|---|---|---|---|---|---|
| pentanal | 6.54±1.10 | 54.20±10.14 | 207.11±43.35 | 520.90±94.64 | 53.17±17.47 | 67.88±14.14 | 2.44±0.17 | 52.33±9.66 |
| **hexanal** | **13.79±0.29** | **100.99±10.82** | **113.73±15.77** | **1116.47±102.74** | **49.53±9.06** | **110.27±29.89** | **10.87±0.58** | **78.09±9.69** |
| heptanal | 59.04±11.90 | 106.02±13.39 | 245.78±53.52 | 1151.97±205.83 | 164.50±20.16 | 167.65±61.91 | 42.41±4.75 | 82.03±8.14 |
| octanal | 55.79±2.87 | 125.93±10.44 | 111.36±20.80 | 382.27±97.31 | 99.58±12.05 | 130.40±25.12 | 54.95±2.66 | 85.43±12.74 |
| nonanal | 127.45±28.71 | 1149.65±99.26 | 238.54±53.88 | 1441.81±495.30 | 313.74±88.89 | 1207.32±188.59 | 159.71±43.63 | 715.49±151.51 |
| furfural | 4.75±0.28 | 4.44±0.70 | 14.82±1.04 | 7.84±2.51 | 6.30±0.83 | 3.94±0.91 | 4.36±1.25 | 1.63±0.18 |
| benzaldehyde | 5.12±1.08 | 50.74±3.94 | 12.82±3.30 | 84.56±14.63 | 14.03±3.54 | 58.23±11.28 | 8.69±3.42 | 47.89±8.19 |
| phenyl acetaldehyde | 136.96±14.09 | 60.96±12.16 | 40.16±12.87 | 48.30±24.96 | 108.26±24.93 | 48.74±19.27 | 138.62±33.69 | 43.47±8.56 |
| 5-methylfuran | 1.09±0.22 | 1.07±0.43 | 2.37±0.87 | 1.08±0.54 | 1.99±0.36 | 0.51±0.19 | 2.19±0.92 | 0 |
| (E)-2-hexenal | 44.14±19.86 | 35.22±4.11 | 195.4±38.40 | 306.92±66.34 | 107.86±29.36 | 45.63±18.86 | 47.03±9.96 | 6.13±0.84 |
| (E)-2-heptenal | 59.05±21.22 | 45.92±6.74 | 145.68±52.85 | 167.35±22.16 | 133.10±28.52 | 47.18±6.55 | 104.89±27.36 | 39.25±4.25 |
| (E)-2-octenal | 73.05±11.21 | 36.97±3.46 | 78.74±14.57 | 427.05±143.72 | 73.86±23.26 | 254.64±48.20 | 41.22±6.62 | 157.61±18.87 |
| (E)-2-nonenal | 18.08±2.31 | 53.71±5.98 | 23.83±3.50 | 51.92±17.95 | 27.51±10.38 | 54.91±7.44 | 14.18±1.80 | 53.70±4.16 |
| (E)-2-decenal | 29.99±3.51 | 45.55±2.73 | 12.29±3.59 | 30.02±11.60 | 35.35±12.23 | 68.03±16.67 | 21.26±4.43 | 67.60±8.71 |
| 2-undecenal | 0.88±0.18 | 1.29±0.19 | 0 | 0 | 0.43±0.06 | 1.60±0.41 | 0.84±0.55 | 3.20±0.39 |
| (E,E)-2,4-decadienal | 91.29±12.56 | 82.14±6.43 | 69.84±3.09 | 60.41±2.85 | 75.39±2.94 | 77.97±7.91 | 60.93±18.30 | 107.14±5.62 |
| **methional** | **77.78±12.35** | **30.26±3.90** | **80.01±20.18** | **139.69±44.55** | **215.22±32.93** | **71.20±25.26** | **110.87±17.06** | **27.22±3.99** |
| **dimethyl trisulfide** | **12.58±0.73** | **45.55±6.43** | **16.43±4.01** | **117.35±10.72** | **23.56±5.29** | **66.81±21.91** | **17.33±2.38** | **23.08±4.81** |
| **2-furfurylthiol** | **27.22±1.87** | **112.96±14.66** | **18.10±6.17** | **157.40±43.06** | **32.58±6.11** | **158.03±40.61** | **30.83±4.72** | **56.55±7.12** |
| 2,3-octanedione | 3.22±0.32 | 0.72±0.15 | 5.18±1.68 | 238.66±48.43 | 1.94±0.35 | 42.08±21.44 | 1.11±0.13 | 35.84±8.75 |
| 2-nonanone | 0.59±0.08 | 0.85±0.07 | 0.01±0.00 | 0.01±0.00 | 0.01±0.00 | 0.01±0.00 | 1.23±0.28 | 0.40±0.10 |
| 2(5H)-furanone | 4.03±0.56 | 5.08±0.86 | 9.99±2.24 | 4.15±0.50 | 6.24±1.98 | 4.12±0.23 | 4.02±0.94 | 4.29±0.48 |
| 1-pentanol | 0.92±0.10 | 1.21±0.39 | 3.64±1.20 | 4.07±1.86 | 1.69±0.27 | 2.53±0.77 | 0 | 0 |
| 1-octene-3-ol | 5.91±1.43 | 5.86±2.52 | 16.36±4.41 | 8.20±3.21 | 5.74±1.32 | 3.91±1.91 | 12.45±4.76 | 7.54±2.37 |
| 2-pentylfuran | 2.87±0.37 | 8.47±0.60 | 5.65±1.39 | 17.93±3.04 | 5.30±1.21 | 8.89±2.09 | 1.94±0.15 | 4.52±0.83 |
| 2-acetylfuran | 5.11±0.82 | 3.60±1.75 | 9.48±3.16 | 7.63±1.25 | 7.06±1.19 | 8.16±2.38 | 6.32±2.43 | 3.69±0.39 |
| 2-methylpyrazine | 23.02±3.74 | 9.14±0.70 | 36.31±13.20 | 24.04±4.28 | 49.90±4.46 | 19.18±7.25 | 42.73±16.08 | 8.00±0.81 |
| 2,5-dimethylpyrazine | 6.33±1.29 | 1.72±0.30 | 5.38±1.49 | 5.62±0.78 | 12.63±1.97 | 4.26±1.65 | 8.36±1.29 | 1.73±0.12 |
| 2-ethylpyrazine | 0 | 59.22±8.98 | 2.63±0.78 | 0 | 3.94±0.60 | 0 | 3.86±1.22 | 1.41±0.29 |
| 2,4,6-trimethyl-1,3,5-triazine | 0 | 0 | 0 | 0 | 3.00±0.81 | 0.80±0.22 | 0 | 0 |
| 2-ethyl-6-methylpyrazine | 2.27±0.49 | 1.02±0.08 | 2.53±0.86 | 2.89±1.12 | 4.21±0.83 | 1.60±0.68 | 4.43±1.40 | 0.50±0.14 |
| pyrrole | 1.15±0.06 | 0.86±0.08 | 4.67±1.47 | 0 | 5.27±1.11 | 0 | 3.78±1.65 | 0 |
| guaiacol | 18.72±3.00 | 22.71±11.13 | 60.96±104.6 | 19.78±5.69 | 10.11±3.70 | 9.19±3.05 | 4.75±1.22 | 0.45±0.31 |
| maltol | 8.05±1.32 | 2.78±1.00 | 24.91±4.58 | 0 | 6.60±0.50 | 1.14±0.19 | 2.92±0.65 | 1.93±0.40 |
| phenol | 12.67±0.76 | 10.37±4.76 | 10.82±2.55 | 5.20±3.07 | 1.37±0.49 | 0.70±0.29 | 2.61±1.29 | 2.71±0.69 |
| hexenyl acetate | 19.42±1.52 | 0 | 24.50±8.48 | 0 | 19.32±4.46 | 0 | 25.97±11.43 | 0 |
| styrene | 1.33±0.14 | 2.88±0.44 | 2.42±1.15 | 0 | 3.97±0.89 | 4.48±2.42 | 1.65±0.22 | 2.21±0.23 |
| acetic acid | 0 | 0 | 44.66±5.56 | 0 | 0 | 0 | 0 | 0 |
| propanoic acid | 1.69±0.20 | 1.41±0.27 | 6.28±1.90 | 3.90±1.72 | 6.93±0.67 | 1.42±0.53 | 0.84±0.17 | 0.84±0.15 |
| butanoic acid | 2.15±0.40 | 1.62±0.28 | 7.52±1.80 | 2.45±0.84 | 8.62±4.01 | 1.83±0.20 | 2.17±0.46 | 1.82±0.26 |
| isovaleric acid | 1.38±0.49 | 0 | 5.34±1.10 | 11.07±2.91 | 2.08±1.82 | 0 | 0.72±0.14 | 0.87±0.11 |
| hexanoic acid | 3.79±0.66 | 4.79±0.80 | 8.57±2.15 | 25.72±4.05 | 4.81±1.26 | 5.66±0.76 | 2.68±0.31 | 6.92±1.15 |

"Muscle" = the paper's "breast" column (breast muscle). [!] Quanjude skin guaiacol "60.96±104.6" as printed (SD > mean).

### Table 4. Compounds with OAV > 1 (OAV = conc / threshold in water, µg/L)

| compound | threshold | DXL skin | DXL muscle | QJD skin | QJD muscle | BYF skin | BYF muscle | DD skin | DD muscle |
|---|---|---|---|---|---|---|---|---|---|
| pentanal | 22 (b) | <1 | 2 | 9 | 24 | 2 | 3 | <1 | 2 |
| hexanal | 4.5 (c) | 3 | 22 | 25. [sic] | 248 | 11 | 25 | 2 | 17 |
| heptanal | 3 (d) | 20 | 35 | 82 | 384 | 55 | 56 | 14 | 27 |
| octanal | 0.7 (e) | 80 | 180 | 159 | 546 | 142 | 186 | 79 | 122 |
| nonanal | 8 (f) | 16 | 144 | 30 | 180 | 39 | 151 | 20 | 89 |
| phenyl acetaldehyde | 4 (d) | 34 | 15 | 10 | 12 | 27 | 12 | 35 | 11 |
| (E)-2-hexenal | 40 (a) | 1 | <1 | 5 | 8 | 3 | 1 | 1 | <1 |
| (E)-2-heptenal | 13 (b) | 5 | 4 | 11 | 13 | 10 | 4 | 8 | 3 |
| (E)-2-octenal | 4 (g) | 18 | 9 | 20 | 107 | 18 | 64 | 10 | 39 |
| (E)-2-nonenal | 0.25 (e) | 72 | 215 | 95 | 208 | 110 | 220 | 57 | 215 |
| (E)-2-decenal | 0.4 (h) | 75 | 114 | 31 | 75 | 88 | 170 | 53 | 169 |
| (E,E)-2,4-decadienal | 0.2 (e) | 456 | 411 | 349 | 302 | 377 | 390 | 305 | 536 |
| **methional** | **0.2 (e)** | 389 | 151 | 400 | 698 | 1076 | 356 | 554 | 136 |
| **dimethyl trisulfide** | **0.005 (d)** | 2515 | 9110 | 3286 | 23470 | 4711 | 13362 | 3466 | 4616 |
| **2-furfurylthiol** | **0.005 (d)** | 5444 | 22593 | 3620 | 31480 | 6517 | 31606 | 6166 | 11311 |
| 1-octene-3-ol | 1 (f) | 6 | 6 | 16 | 8 | 6 | 4 | 12 | 8 |
| 2-pentylfuran | 6 (d) | <1 | 1 | 1 | 3 | <1 | 1 | <1 | <1 |
| guaiacol | 3 (d) | 6 | 8 | 20 | 7 | 3 | 3 | 2 | <1 |

Threshold sources: (a) literature; (b) Zhu 2016; (c) van Gemert 2003; (d) "flavor base"; (e) Kerler & Grosch 1997; (f) Czerny 2008 [printed "2009"]; (g) Greger & Schieberle 2007; (h) Christlbauer & Schieberle 2011.
[Z] Recompute checks: FFT BYF muscle 158.03/0.005 = 31606 ✓; QJD muscle 157.40/0.005 = 31480 ✓; DMTS QJD muscle 117.35/0.005 = 23470 ✓; methional BYF skin 215.22/0.2 = 1076 ✓; hexanal QJD muscle 1116.47/4.5 = 248 ✓. The table is internally consistent with Table 3.

### Table 5. Contribution rate (% of summed OAV; mean ± SD)

| compound | DXL skin | DXL muscle | QJD skin | QJD muscle | BYF skin | BYF muscle | DD skin | DD muscle |
|---|---|---|---|---|---|---|---|---|
| pentanal | 0 | 0.01±0.00 | 0.13±0.05 | 0.04±0.01 | 0.02±0.01 | 0.01±0.00 | 0 | 0.01±0.00 |
| hexanal | 0.03±0.00 | 0.07±0.01 | 0.32±0.08 | 0.45±0.12 | 0.08±0.02 | 0.05±0.01 | 0.02±0.00 | 0.10±0.02 |
| heptanal | 0.22±0.04 | 0.11±0.02 | 1.01±0.09 | 0.67±0.10 | 0.42±0.06 | 0.12±0.02 | 0.13±0.01 | 0.16±0.02 |
| octanal | 0.87±0.04 | 0.55±0.07 | 1.98±0.23 | 0.94±0.14 | 1.10±0.18 | 0.41±0.05 | 0.74±0.10 | 0.71±0.12 |
| nonanal | 0.18±0.05 | 0.44±0.04 | 0.40±0.18 | 0.31±0.09 | 0.29±0.05 | 0.34±0.07 | 0.19±0.07 | 0.52±0.10 |
| phenyl acetaldehyde | 0.37±0.03 | 0.05±0.01 | 0.12±0.03 | 0.02±0.01 | 0.20±0.04 | 0.03±0.00 | 0.32±0.07 | 0.06±0.01 |
| (E)-2-hexenal | 0.01±0.01 | 0 | 0.06±0.01 | 0.01±0.00 | 0.02±0.00 | 0 | 0.01±0.00 | 0 |
| (E)-2-heptenal | 0.05±0.02 | 0.01±0.00 | 0.13±0.02 | 0.02±0.00 | 0.08±0.00 | 0.01±0.00 | 0.08±0.03 | 0.02±0.00 |
| (E)-2-octenal | 0.20±0.02 | 0.03±0.00 | 0.25±0.05 | 0.18±0.04 | 0.14±0.02 | 0.14±0.02 | 0.10±0.02 | 0.23±0.05 |
| (E)-2-nonenal | 0.79±0.08 | 0.66±0.14 | 1.24±0.38 | 0.35±0.07 | 0.81±0.18 | 0.49±0.12 | 0.54±0.11 | 1.26±0.19 |
| (E)-2-decenal | 0.82±0.08 | 0.35±0.03 | 0.38±0.12 | 0.13±0.05 | 0.67±0.19 | 0.37±0.06 | 0.50±0.11 | 0.99±0.16 |
| (E,E)-2,4-decadienal | 4.99±0.59 | 1.25±0.12 | 4.52±1.31 | 0.53±0.08 | 2.92±0.49 | 0.87±0.18 | 2.86±0.89 | 3.14±0.40 |
| methional | 4.25±0.61 | 0.46±0.06 | 4.90±0.44 | 1.20±0.24 | 8.24±1.19 | 0.76±0.12 | 5.24±1.17 | 0.79±0.11 |
| dimethyl trisulfide | 27.53±0.81 | 27.56±0.54 | 40.51±5.76 | 41.42±6.54 | 35.62±4.28 | 28.34±1.62 | 32.10±0.78 | 26.50±2.00 |
| 2-furfurylthiol | 59.54±1.39 | 68.41±0.55 | 43.54±6.77 | 53.69±6.45 | 49.28±3.60 | 68.05±1.38 | 57.03±2.13 | 65.44±1.32 |
| 1-octene-3-ol | 0.06±0.02 | 0.02±0.01 | 0.21±0.05 | 0.01±0.01 | 0.04±0.01 | 0.01±0.00 | 0.12±0.06 | 0.04±0.01 |
| 2-pentylfuran | 0.01±0.00 | 0 | 0.01±0.00 | 0.01±0.00 | 0.01±0.00 | 0 | 0 | 0 |
| guaiacol | 0.07±0.01 | 0.02±0.01 | 0.26±0.43 | 0.01±0.00 | 0.03±0.01 | 0.01±0.00 | 0.01±0.00 | 0 |

### [Z] Molar view of the engine's targets (MW: FFT 114.17; methional 104.17; DMTS 126.26; hexanal 100.16)

| compound | range (ng/g) | range (nmol/kg) | median of 8 cells (ng/g) |
|---|---|---|---|
| 2-furfurylthiol | 18.1 – 158.0 | 159 – 1384 | 44.6 |
| methional | 27.2 – 215.2 | 261 – 2066 | 79 |
| dimethyl trisulfide | 12.6 – 117.4 | 100 – 929 | 23.3 |
| hexanal | 10.9 – 1116 | 109 – 11 147 | 89.5 |
| MFT, MFT disulfide, DMDS | not reported | – | – |

Skin vs muscle: FFT, DMTS and hexanal are higher in muscle in all four brands; methional is higher in skin in
three of four (Quanjude the exception).

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor (external real-food levels, weak class):**
- **A roasted-poultry upper band** for FFT (18-158 µg/kg), methional (27-215 µg/kg) and DMTS (13-117 µg/kg),
  and a wide hexanal band (11-1116 µg/kg). Use as "roast duck, HS-SPME + matrix-matched IS calibration,
  order-of-magnitude", not alongside SIDA numbers as equals.
- **DMTS at tens of µg/kg in a dry-roasted meat** is a direct statement that the methionine-side
  polysulfide exists at levels comparable to methional in roast poultry. The engine refuses DMTS
  (`METHIONINE_NO_DMTS_REASON`); this paper gives a real-food level, not a mechanism.
- **Skin vs muscle contrast** (higher FFT/DMTS/hexanal in muscle, higher methional mostly in skin) as a
  qualitative within-food pattern.

**Cannot anchor:**
- **No kinetics, no precursors** (no ribose, cysteine, thiamine, fat composition, pH, water loss).
- **No defined thermal history:** oven air 250 °C (200-250 °C in the Discussion), 60 or 80 min, restaurant
  practice; product/skin temperature unknown. Not a temperature-defined pot.
- **MFT absence is not informative:** HS-SPME on CAR/PDMS with full-scan ID and no MFT standard purchased
  (the chemicals list has no MFT); absence reflects the method, not the duck.

**Caveats / flags:**
1. [!] **Quantification class is the weakest here.** HS-SPME at 55 °C of minced meat; one non-isotopic IS
   (2-methyl-3-heptanone); calibration in solvent-extracted, freeze-dried Quanjude **skin**, applied also to
   **muscle** (different fat/water, hence different headspace partition). No recovery or matrix-effect
   check reported. My judgement: true values within a factor 3-5 at ~60 % confidence. FFT may be
   overestimated: it co-elutes on DB-WAX with (E)-2-octenal (both LRI 1395, Table 1), so its SIM ions
   (81, 114, 53, 82) carry interference risk from an aldehyde present at 37-427 ng/g.
2. [!] **Table 2 has copy errors:** hexanal and heptanal equations lack "x"; **nonanal and (E)-2-octenal
   have the same equation and R² (y = 0.00008x − 358.28/−358.25, 0.9923)**; the footnote says y is "mg/g"
   while Table 3 reports ng/g; large intercepts (+174.49 decadienal, +123.65 octanal) relative to many
   sample levels. The calibration cannot be re-applied from the paper.
3. [!] **FFT levels are far above other meats:** 36-316× Guth & Grosch beef juice (0.5 µg/kg, SIDA),
   7-60× Fan 2018 chicken broth (2.7 µg/kg), 3-60× Du 2020 grilled mutton (2.6-5.9 µg/kg). Roasting a skin
   does push FFT up (furfural + H2S on a dry surface), but the method spread is at least as large as the
   real difference. Do not set a target from these alone.
4. [!] Thresholds for FFT and DMTS (0.005 µg/kg) are from "flavor base" software; OAVs and contribution
   rates scale with them. The concentrations do not depend on them.
5. [!] Table 3's "5-methylfuran" is not a valid name; Quanjude skin guaiacol SD exceeds the mean.
6. PDF is the Just Accepted manuscript; typeset version JAFC 67(20) 5847-5856 may differ in minor edits.

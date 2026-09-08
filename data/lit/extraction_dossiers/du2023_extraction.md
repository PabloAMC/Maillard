# Du et al. 2023 — EXTRACTION (cysteine + glucose + one C9 aldehyde of 0, 1 or 2 double bonds, pH 6.5 phosphate, 150 C / 120 min; HS-SPME-GC-MS with one internal standard, µg/L; CAMOLA 13C6-glucose labelling)
### The paper that puts 2-pentylthiophene beside 2-butylpyridine from the same 2,4-nonadienal pot, with the 2-enal and the alkanal as built-in negative controls; H2S and NH3 come from cysteine's own Strecker degradation.

**Source on disk:** `data/articles/Du2023.pdf` (owner's download, 2026-09-08). Read from the `pypdf` text
layer; Table 1 (88 compounds x 5 systems) and Table 2 (isotopomer distributions) came through with the
dashes preserved, so column placement could be checked against the printed class sums (all thiazole
sums reproduce; the pyridine sum for system B does not, flag 3). Fig. 1 (a)-(h) and Fig. 2 are drawn
ChemDraw schemes: **FIGURE-ONLY**; the text describes them in one clause each and only those clauses are
used. Fig. 3 (heatmap, PLS-DA) and Fig. S1 are not used. Table S1 (aldehyde-derived volatiles: acids,
alcohols, ketones) is in the Supplement, not available.

## 0. Identity

| field | value |
|---|---|
| Title | "The number and position of unsaturated bonds in aliphatic aldehydes affect the cysteine-glucose Maillard reaction: Formation mechanism and comparison of volatile compounds" |
| Authors | Wenbin Du, Yutang Wang, Qinghong Yan, Shuang Bai, Yatao Huang, Long Li, Yuwen Mu, Ashbala Shakoor, Bei Fan, Fengzhong Wang (CAAS Institute of Food Science and Technology, Beijing) |
| Venue | Food Research International 173 (2023) 113337 |
| DOI | 10.1016/j.foodres.2023.113337 |
| Naming | System M = cysteine + glucose (control); A = + nonanal; B = + (E)-2-nonenal; C = + (E,E)-2,4-nonadienal; D = + (E,Z)-2,6-nonadienal. Table 1 compound numbers (No. 1-88) are used below. CAMOLA = carbohydrate module labelling (1 : 1 [13C6]/[12C6]-glucose). |
| Companions | Zhao et al. 2019 (Food Chem 270 and 274: the same group's CAMOLA on glutathione/glucose + fat and cysteine/xylose/glycine); Wang et al. 2020 (Food Chem 305: initial intermediates of glucose + glutathione + aldehydes); Mottram 1998 (the 2-alkylthiophene mechanism they cite); Adams et al. 2011 (amino-acid-catalysed 2-alkylfuran from 2-alkenals) |

## 1. Why it matters

Programme 7 names 2-alkylthiophenes and alkylthiazoles as the sulfur half of the lipid-Maillard cross
products. This paper heats cysteine + glucose (the engine's Cys and Glc) with one C9 aldehyde at a time
and reports every sulfur volatile with a concentration. It gives (i) **2-pentylthiophene 238.65 ± 23.37
µg/L from the 2,4-dienal system only** (none from nonanal, 2-nonenal, 2,6-nonadienal or the control),
with 92 % of its molecules carrying no glucose carbon; (ii) the same pot's **2-butylpyridine 191.24 ±
1.49 µg/L** (the C9 analogue of 2-pentylpyridine from 2,4-decadienal), with the ammonia from cysteine's
Strecker degradation; (iii) 2-hexylthiophene (10.77 / 2.58 / 6.38 µg/L from nonanal / 2-nonenal /
2,4-nonadienal, all unlabelled: the extra carbon is not from glucose); (iv) 2-heptylthiophene 471.07
µg/L from nonanal, 16 % doubly labelled (aldol with acetaldehyde first); (v) the suppression of the
control's sulfur volatiles by every aldehyde (thiophene sum 2067 -> 537-1181 µg/L; thieno[3,2-b]thiophene
1038 -> 0-420) and the general statement that dienals inhibit more than the enal or the alkanal; and
(vi) the thiazole slate per system. **It does not detect 2-pentyl-4-methylthiazole or
2-hexyl-4-methylthiazole** (the roadmap's named targets) in any system: the thiazoles it finds are
short-chain (thiazole, 2-ethyl-, 2,5-dimethyl-, 2,4,5-trimethyl-, 4,5-dimethyl-, 2-acetylthiazole).

## 2. Methods as they matter to a model

- **Charge (section 2.2):** L-cysteine 1 mmol + D-glucose 1 mmol + aldehyde 0.05 mmol in 5 mL of 0.2 M
  phosphate buffer pH 6.50, 15 mL screw-cap glass vial. So Cys 200 mM, Glc 200 mM, aldehyde 10 mM;
  aldehyde : Cys = 1 : 20. Control M without aldehyde. Aldehydes 99 % (J&K): nonanal, (E)-2-nonenal,
  (E,E)-2,4-nonadienal, (E,Z)-2,6-nonadienal. Aldehyde solubility in the buffer is not discussed.
- **Heating:** 150 C oil bath, magnetic stirring, 120 min; one time point.
- **HS-SPME:** internal standard "1 µL of 1,2-dichlorobenzene (100 µg, in 1 mL of methanol)" added to
  the 5 mL (read: 1 µL of a 100 µg/mL solution = 0.1 µg IS = 20 µg/L; the sentence is ambiguous, flag
  5); CAR/PDMS/DVB 50/35 µm; equilibrate 50 C 30 min; expose 50 C 30 min; desorb 250 C 3 min splitless.
- **GC-MS:** Agilent 7890B / 5977A; DB-Wax 30 m x 0.25 mm x 0.25 µm; He 1 mL/min; 40 C (2 min) ->
  210 C at 3.5 C/min -> 240 C at 10 C/min; EI 70 eV; m/z 40-450. GC-O on the same method (Sniffer 9000,
  three trained panelists).
- **Quantification:** c_x = S_x · c_i / S_i, one internal standard, **response factor 1 for every
  compound, no calibration curves, no recovery**. Units µg/L of the 5 mL reaction. Three replicates;
  ANOVA + Duncan, p < 0.05. Numbers are therefore semi-quantitative levels, comparable within a compound
  across systems and only loosely across compounds.
- **Identification:** RI (C5-C30 on DB-Wax) + NIST 15 MS; "S" = agrees with an authentic standard; "O"
  = odour detected at GC-O. Authentic standards were run for 2-pentylthiophene, 2-hexylthiophene,
  2-pentylpyridine, 2-methyl- to 2-butylthiophene, 2-methylthiazole, 2-ethylthiazole,
  2,4,5-trimethylthiazole, 4,5-dimethylthiazole and others (section 2.1). 2-Heptylthiophene,
  2-butylpyridine, thiazole and 2-acetylthiazole are RI/MS(/O) without a standard.
- **CAMOLA (section 2.6):** same systems with 0.5 mmol [13C6]- + 0.5 mmol [12C6]-glucose; isotopomer
  fractions from M+ to M+n after correcting for natural 13C, 33S, 34S and for M-1. A compound made from
  an intact glucose skeleton shows 1 : 1 unlabelled : fully labelled; one made only from the aldehyde
  (or from cysteine) is 99 % unlabelled; partial labelling means a glucose fragment was incorporated.

## 3. Tables re-typed

### Table 1. "Volatile flavor compounds derived from the Maillard reaction or the interaction of the Maillard reaction with aliphatic aldehydes." Quantity µg/L, mean ± SD (n = 3); "–" = not detected; letters compare systems within a row; RI on DB-Wax; ID = identification methods.

Conditions for every column: Cys 1 mmol + Glc 1 mmol (+ aldehyde 0.05 mmol) in 5 mL pH 6.5 phosphate,
150 C, 120 min. M control; A nonanal; B (E)-2-nonenal; C (E,E)-2,4-nonadienal; D (E,Z)-2,6-nonadienal.

**Furans**

| No. | RI | compound | odour | M | A | B | C | D | ID |
|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 1 | 955 | 2-ethylfuran | – | 4.44 ± 0.25 | – | – | – | – | RI/MS/S |
| 2 | 1032 | 2-propylfuran | – | – | – | – | 4.35 ± 1.03 | – | RI/MS/S |
| 3 | 1070 | 2-vinylfuran | – | 0.06 ± 0.01 | – | – | – | – | RI/MS/S |
| 4 | 1120 | 2-butylfuran | – | – | – | – | 8.45 ± 2.09 | – | RI/MS/S |
| 5 | 1152 | 2-butyltetrahydrofuran | – | 36.77 ± 6.02 | – | – | – | – | RI/MS |
| 6 | 1230 | **2-pentylfuran** | green bean | – | – | 366.97 ± 69.66 a | 7.79 ± 0.40 b | – | RI/MS/O/S |
| 7 | 1416 | (E)-2-(2-pentenyl)furan | coffee | – | – | – | 7.19 ± 1.40 b | 133.01 ± 4.54 a | RI/MS/O |
| 8 | 1493 | 1-(2-furanyl)ethanone | – | 49.13 ± 16.55 | 31.81 ± 1.67 | 44.41 ± 24.81 | 37.86 ± 1.24 | – | RI/MS/S |
| 9 | 1519 | furfural | caramel | – | – | – | 26.71 ± 1.87 | – | RI/MS/O/S |
| 10 | 1654 | 2-furanmethanol | burnt | – | – | 31.07 ± 0.71 | – | – | RI/MS/O/S |
| 11 | 1750 | 2-hexanoylfuran | – | – | – | – | 106.58 ± 6.00 | – | RI/MS/S |
| 12 | 1776 | 1(3H)-isobenzofuranone | herb | – | – | – | 186.34 ± 0.00 | – | RI/MS/O |
| 13 | 1900 | 5-butyldihydro-2(3H)-furanone | caramel | – | 12.77 ± 3.27 | – | – | – | RI/MS/O/S |
| 14 | 2015 | dihydro-5-pentyl-2(3H)-furanone | caramel | – | 198.32 ± 12.39 a | 54.29 ± 6.59 bc | 59.83 ± 3.91 bc | – | RI/MS/O |
| 15 | 2026 | 2,5-dimethyl-4-hydroxy-3(2H)-furanone | caramel | 9.71 ± 3.84 ab | 3.84 ± 0.82 c | 5.16 ± 1.43 ab | – | – | RI/MS/O/S |
| 16 | 2317 | 3-formylfuran-2-carboxylic acid | – | 11.93 ± 1.67 | – | – | – | – | RI/MS |
| | | furan sum | | 112.03 | 246.74 | 501.91 | 445.10 | 133.01 | |

**Thiophenes**

| No. | RI | compound | odour | M | A | B | C | D | ID |
|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 17 | 1021 | thiophene | garlic | 7.67 ± 0.86 | – | – | – | – | RI/MS/O |
| 18 | 1084 | 2-methylthiophene | sulfur | 12.15 ± 11.46 a | 0.25 ± 0.01 b | – | – | – | RI/MS/O/S |
| 19 | 1142 | 2,4-dimethylthiophene | – | 12.34 ± 5.68 | – | – | – | – | RI/MS |
| 20 | 1156 | 2,3-dihydro-5-methylthiophene | – | 12.11 ± 1.87 | – | – | – | – | RI/MS |
| 21 | 1157 | 2-ethylthiophene | meaty | 19.22 ± 0.82 a | – | – | – | 7.34 ± 0.52 b | RI/MS/O/S |
| 22 | 1197 | 2,3-dimethylthiophene | – | 16.06 ± 3.86 | – | – | – | – | RI/MS |
| 23 | 1234 | 2-propylthiophene | meaty | – | – | – | 14.54 ± 5.99 b | 105.33 ± 14.82 a | RI/MS/O/S |
| 24 | 1250 | 3,4-dimethylthiophene | – | – | – | – | – | 9.73 ± 2.03 | RI/MS |
| 25 | 1265 | 2-(1-methylethyl)thiophene | – | 0.06 ± 0.01 | – | – | – | – | RI/MS |
| 26 | 1770 | 1-(2-thienyl)ethanone (2-acetylthiophene) | – | 50.99 ± 6.53 a | 19.42 ± 0.62 d | 138.32 ± 0.15 bc | – | 32.97 ± 5.35 bc | RI/MS |
| 27 | 1340 | 2-butylthiophene | fruity | – | – | – | 43.21 ± 2.96 | – | RI/MS/O/S |
| 28 | 1452 | **2-pentylthiophene** | fruit, meaty | – | – | – | **238.65 ± 23.37** | – | RI/MS/O/S |
| 29 | 1514 | dihydro-2-methyl-3(2H)-thiophenone | sulfur | 8.47 ± 0.52 | – | – | – | – | RI/MS/O |
| 30 | 1546 | dihydro-3(2H)-thiophenone | burnt garlic | 6.72 ± 1.22 | – | – | – | – | RI/MS/O |
| 31 | 1554 | 2-thiophenethiol | burnt | 18.63 ± 3.51 | – | – | – | – | RI/MS/O/S |
| 32 | 1590 | **2-hexylthiophene** | meaty | – | 10.77 ± 0.40 a | 2.58 ± 0.31 c | 6.38 ± 0.15 b | – | RI/MS/O/S |
| 33 | 1634 | 2-(methylthio)thiophene | – | 15.84 ± 1.00 | – | – | – | – | RI/MS |
| 34 | 1676 | 2-thiophenecarboxaldehyde | meaty | 16.96 ± 0.63 a | 3.04 ± 0.35 b | – | – | – | RI/MS/O/S |
| 35 | 1683 | 2-ethyl-5-propylthiophene | garlic | – | – | – | 25.70 ± 0.28 | – | RI/MS/O/S |
| 36 | 1698 | **2-heptylthiophene** | meaty | – | 471.07 ± 17.11 | – | – | – | RI/MS |
| 37 | 1745 | 3,4-diethylthiophene | – | 13.40 ± 0.33 | – | – | – | – | RI/MS |
| 38 | 1785 | 5-methyl-2-thiophenecarboxaldehyde | – | 53.09 ± 2.68 de | 44.31 ± 13.28 de | 123.87 ± 7.98 c | 113.61 ± 12.54 a | 188.58 ± 16.62 b | RI/MS/S |
| 39 | 1795 | 3-methyl-2-thiophenecarboxaldehyde | herb | 201.89 ± 18.07 | – | – | – | – | RI/MS/O/S |
| 40 | 1812 | 2-methoxy-5-methylthiophene | – | 10.26 ± 3.13 | – | – | – | – | RI/MS |
| 41 | 1821 | 1-(2-thienyl)-1-propanone | cream | 15.19 ± 4.27 | – | – | – | – | RI/MS/O |
| 42 | 1851 | 2,5-diethylthiophene | – | 8.28 ± 2.38 | – | – | – | 2.60 ± 0.47 | RI/MS/S |
| 43 | 1856 | thieno[3,2-b]thiophene | meaty | 1037.62 ± 81.63 a | 419.56 ± 10.85 b | 147.00 ± 4.92 c | – | 109.35 ± 0.55 d | RI/MS/O/S |
| 44 | 1864 | 5-methylbenzo[b]thiophene | – | 10.03 ± 2.35 | – | – | – | – | RI/MS |
| 45 | 1894 | 2,5-thiophenedicarboxaldehyde | – | 32.23 ± 0.04 a | – | 14.10 ± 1.85 cd | 14.86 ± 1.17 cd | 17.40 ± 1.30 b | RI/MS/S |
| 46 | 1917 | 2,5-dipropylthiophene | – | – | – | – | 99.69 ± 11.74 | – | RI/MS/S |
| 47 | 1922 | 2-methylthieno[2,3-b]thiophene | meaty | 375.69 ± 13.83 a | 211.14 ± 6.55 b | 110.91 ± 3.42 c | – | 105.02 ± 0.37 d | RI/MS/O |
| 48 | 1932 | 2-thiophenemethanol | – | 4.63 ± 0.94 | – | – | – | – | RI/MS |
| 49 | 1940 | 5-methyl-2-thiophenecarboxylic acid | – | 5.77 ± 1.47 | – | – | – | – | RI/MS |
| 50 | 1954 | thieno[2,3-b]thiophene | – | 29.62 ± 5.90 | – | – | – | – | RI/MS/S |
| 51 | 1998 | 2-ethyl-5-pentylthiophene | – | – | – | – | 37.78 ± 4.56 | – | RI/MS |
| 52 | 2002 | 1-(2-thienyl)-1-hexanone | – | – | – | – | 42.64 ± 2.02 | – | RI/MS |
| 53 | 2037 | 2-butyl-5-ethylthiophene | – | – | – | – | 46.20 ± 0.96 | – | RI/MS/S |
| 54 | 2103 | 2,5-dibutylthiophene | – | 28.05 ± 1.88 b | – | – | 133.10 ± 4.57 a | – | RI/MS/S |
| 55 | 2219 | 5-acetyl-3-thiopheneacetic acid | – | 26.31 ± 3.64 | – | – | – | – | RI/MS |
| 56 | 2434 | 2-(3-thienylthio)thiophene | – | 18.01 ± 5.33 | – | – | – | – | RI/MS |
| | | thiophene sum | | 2067.32 | 1180.69 | 536.78 | 816.36 | 582.60 | |

**Thiazoles**

| No. | RI | compound | odour | M | A | B | C | D | ID |
|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 57 | 1237 | 2-methylthiazole | cabbage | 0.95 ± 0.06 | – | – | – | – | RI/MS/O/S |
| 58 | 1246 | isothiazole | – | – | – | – | 3.04 ± 0.09 | – | RI/MS |
| 59 | 1250 | 2,5-dimethylthiazole | – | – | – | 4.90 ± 0.66 | – | – | RI/MS |
| 60 | 1265 | thiazole | meaty | – | 1.31 ± 0.14 | – | 82.67 ± 3.19 | – | RI/MS/O |
| 61 | 1303 | 4,5-dihydro-2-methylthiazole | – | 23.42 ± 0.00 | – | – | – | – | RI/MS |
| 62 | 1304 | 2-ethylthiazole | meaty | – | – | – | 93.35 ± 3.79 | – | RI/MS/O/S |
| 63 | 1375 | 2,4,5-trimethylthiazole | meaty | 31.73 ± 2.94 a | 9.51 ± 2.93 d | 20.80 ± 2.84 bc | 14.33 ± 4.33 bc | – | RI/MS/O/S |
| 64 | 1401 | 4,5-dimethylisothiazole | roast, smoke | – | – | – | – | 1.64 ± 0.34 | RI/MS/O |
| 65 | 1410 | 5-ethyl-2-methylthiazole | – | – | – | 7.28 ± 0.59 | – | – | RI/MS |
| 66 | 1440 | 4-methyl-2-(1-methylethyl)thiazole | – | 6.64 ± 0.45 | – | – | – | – | RI/MS |
| 67 | 1451 | 4,5-dimethylthiazole | – | – | – | – | – | 24.09 ± 2.53 | RI/MS/S |
| 68 | 1512 | 5-ethenyl-4-methylthiazole | – | – | – | – | – | 3.77 ± 0.98 | RI/MS |
| 69 | 1633 | 2-acetylthiazole | roast | – | – | – | – | 8.87 ± 1.12 | RI/MS/O |
| 70 | 2299 | 4-methyl-5-thiazoleethanol | meaty | – | 2.75 ± 0.76 | – | – | – | RI/MS/O/S |
| | | thiazole sum | | 62.74 | 13.58 | 32.98 | 193.40 | 38.37 | |

(Column sums recomputed from the rows: M 62.74, A 13.57, B 32.98, C 193.39, D 38.37: placement confirmed.)

**Pyridines**

| No. | RI | compound | odour | M | A | B | C | D | ID |
|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 71 | 1378 | 2,4,6-trimethylpyridine | – | – | – | – | 133.76 ± 1.24 | – | RI/MS |
| 72 | 1469 | **2-butylpyridine** | – | – | – | 52.32 ± 0.03 c | **191.24 ± 1.49 a** | 12.93 ± 2.81 b | RI/MS |
| 73 | 1527 | **2-pentylpyridine** | burnt | – | – | 44.15 ± 9.51 | – | – | RI/MS/O/S |
| 74 | 2072 | 2-(2-phenylethyl)pyridine | – | – | – | 22.27 ± 1.21 | – | – | RI/MS |
| | | pyridine sum (as printed) | | – | – | 68.75 | 325.00 | 12.93 | |

(Row entries for B add to 118.74, not the printed 68.75; C and D sums reproduce. See flag 3.)

**Naphthalenes, benzenes, others**

| No. | RI | compound | odour | M | A | B | C | D | ID |
|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| 75 | 1720 | naphthalene | tar | – | – | – | – | 14.43 ± 2.45 | RI/MS/O |
| 76 | 1875 | 1-methylnaphthalene | – | – | – | 13.22 ± 1.05 | – | – | RI/MS |
| 77 | 1831 | 2-methylnaphthalene | – | – | – | – | 35.12 ± 4.00 a | 15.01 ± 0.49 b | RI/MS |
| | | naphthalene sum | | – | – | 63.22 (printed; rows give 13.22) | 35.12 | 29.45 (printed; rows give 29.44) | |
| 78 | 1188 | propylbenzene | – | – | – | – | 7.00 ± 1.75 | – | RI/MS |
| 79 | 1508 | benzaldehyde | almond | – | – | – | 23.67 ± 1.55 a | 15.52 ± 4.20 b | RI/MS/O/S |
| 80 | 1658 | 2-ethylbenzaldehyde | – | – | – | – | – | 21.11 ± 2.60 | RI/MS |
| | | benzene sum | | – | – | – | 30.67 | 36.63 | |
| 81 | 1290 | 1-(methylthio)butane | – | 5.78 ± 0.28 | – | – | – | – | RI/MS |
| 82 | 1315 | 2-methyl-3-furanthiol | meat broth | 12.03 ± 8.80 | 11.36 ± 1.63 | 8.25 ± 0.25 | – | – | RI/MS/O/S |
| 83 | 1422 | 2-furfurylthiol | coffee | 6.18 ± 2.25 | – | – | – | – | RI/MS/O/S |
| 84 | 1584 | anti-3,5-dimethyl-1,2,4-trithiolane | garlic | 190.60 ± 7.15 | – | – | – | – | RI/MS/O |
| 85 | 1842 | syn-3,5-dimethyl-1,2,4-trithiolane | garlic | 108.09 ± 9.49 | 104.28 ± 5.27 | 106.13 ± 5.69 | 154.54 ± 8.78 | 99.37 ± 3.18 | RI/MS/O |
| 86 | 1730 | 1,2,3-trithiane | sulfur | 15.42 ± 0.85 cd | 15.58 ± 2.77 cd | 70.34 ± 5.55 a | – | 28.89 ± 2.07 b | RI/MS/O |
| 87 | 1919 | 1,3,5-trithiane | – | – | – | 4.13 ± 0.26 b | – | 8.15 ± 0.63 a | RI/MS |
| 88 | 1959 | 1-(1H-pyrrol-2-yl)ethanone | – | – | – | 23.82 ± 1.39 c | 28.33 ± 1.57 b | 68.41 ± 3.56 a | RI/MS |
| | | others sum | | 338.10 | 131.22 | 212.67 | 182.86 | 204.82 | |
| | | **Total** | | 2580.19 | 1572.23 | 1334.34 | 1074.91 | 958.81 | |

Number of compounds detected per system (text): M 46, A 18, B 23, C 32, D 23 (Table 1 only; Table S1
adds 19 / 27 / 11 / 11 aldehyde-derived volatiles for A / B / C / D).

### Table 2. "Isotopic distribution pattern of the identified aroma compounds that had 13C-labeled isotopomers detected in the CAMOLA experiment." Percent of molecules with 0, 1, 2, ... 13C atoms (M+ = molecular ion; S = system; "< 1" omitted below).

| No. | compound | M+ | system | isotopomer fractions (%) | reading |
|---:|---|---:|---|---|---|
| 2 | 2-propylfuran | 110 | A; C | 0: 99 | aldehyde only |
| 4 | 2-butylfuran | 124 | C | 0: 99 | aldehyde only |
| 6 | 2-pentylfuran | 138 | B; C | 0: 99 | aldehyde only |
| 7 | (E)-2-(2-pentenyl)furan | 136 | C; D | 0: 99 | aldehyde only |
| 8 | 1-(2-furanyl)ethanone | 110 | A / B / C / D | 0: 51 / 50 / 50 / 49; 6: 49 / 50 / 50 / 51 | intact glucose |
| 9 | furfural | 96 | C | 0: 99 | aldehyde-derived (authors) |
| 10 | 2-furanmethanol | 98 | B | 0: 38; 1: 13; 4: 14; 5: 35 | mixed |
| 11 | 2-hexanoylfuran | 166 | C | 0: 99 | aldehyde only |
| 12 | 1(3H)-isobenzofuranone | 148 | C | 0: 74; 1: 9; 2: 2; 4: 13 | mostly aldehyde |
| 13 | 5-butyldihydro-2(3H)-furanone | 142 | A | 0: 51; 7: 49 | glucose C6 + C1 |
| 14 | dihydro-5-pentyl-2(3H)-furanone | 156 | A; B; C | A 0: 99; B 0: 94, 1: 5, 2: 1; C 0: 99 | aldehyde |
| 15 | furaneol | 128 | A; B | 0: 54 / 53; 6: 46 / 47 | intact glucose |
| 18 | 2-methylthiophene | 98 | A | 0: 50; 5: 50 | glucose C5 (3-deoxypentosone) |
| 21 | 2-ethylthiophene | 112 | D | 0: 99 | aldehyde + H2S |
| 23 | 2-propylthiophene | 130 | C; D | C 0: 99; D 0: 59, 1: 8, 6: 28, 7: 5 | aldehyde (C); mixed (D) |
| 34 | 2-thiophenecarboxaldehyde | 112 | D (first entry); A | D 0: 25, 1: 32, 4: 16, 5: 27; A 0: 99, 1: 1 | mixed (D); aldehyde (A). Listed twice, as printed |
| 26 | 1-(2-thienyl)ethanone | 126 | A / B / C / D | A 0: 51, 4: 26, 6: 22; B 0: 59, 4: 13, 6: 28; C 0: 47, 4: 13, 6: 40; D 0: 47 (rest not printed) | mixed. Listed under system C although Table 1 has no C entry; as printed |
| 27 | 2-butylthiophene | 140 | C | 0: 99 | aldehyde + H2S |
| 28 | **2-pentylthiophene** | 154 | C | **0: 92; 1: 8** | aldehyde + H2S; 8 % carries one glucose carbon |
| 32 | **2-hexylthiophene** | 168 | A; B; C | 0: 99 in all three | aldehyde + one non-glucose carbon + H2S |
| 35 | 2-ethyl-5-propylthiophene | 154 | C | 0: 99 | aldehyde |
| 36 | **2-heptylthiophene** | 182 | A | 0: 84; 2: 16 | nonanal + acetaldehyde (aldol to 2-undecenal) + H2S; the C2 unit is 16 % glucose-derived |
| 38 | 5-methyl-2-thiophenecarboxaldehyde | 126 | A / B / C / D | 0: 36 / 30 / 30 / 31; 1: 4-5; 2: 3-4; 3: 29-33; 4: 24-26; 5: 2-3; 6: 1 | glucose C3/C4 fragment + aldehyde |
| 42 | 2,5-diethylthiophene | 140 | D | 0: 99 | aldehyde |
| 43 | thieno[3,2-b]thiophene | 140 | A; B; D | A, B 0: 73, 1: 7, 2: 20; D 0: 99 | acetaldehyde + crotonaldehyde + H2S |
| 45 | 2,5-thiophenedicarboxaldehyde | 140 | B; C; D | B 0: 50, 4: 39, 6: 11; C 0: 99, 1: 1; D 0: 50, 4: 40, 6: 10 | glucose C4 + acetaldehyde (B, D); aldehyde only (C) |
| 46 | 2,5-dipropylthiophene | 168 | C | 0: 85; 1: 10; 2: 5 | mostly aldehyde |
| 47 | 2-methylthieno[2,3-b]thiophene | 154 | A / B / D | 0: 37 / 31 / 35; 1: 15 / 10 / 7; 2: 37 / 38 / 38; 3: 20 / 21 / 20 | crotonaldehyde + C1/C2 + H2S |
| 51 | 2-ethyl-5-pentylthiophene | 182 | C | 0: 61; 2: 39 | 2-pentylthiophene + C2 (authors: "with ethanol") |
| 52 | 1-(2-thienyl)-1-hexanone | 182 | C | 0: 56; 4: 44 | glucose C4 + hexyl from aldehyde |
| 53 | 2-butyl-5-ethylthiophene | 168 | C | 0: 89; 2: 11 | 2-butylthiophene + C2 |
| 54 | 2,5-dibutylthiophene | 196 | C | 0: 43; 1: 6; 2: 5; 3: 3; 4: 38; 5: 3; 6: 2 | half from a glucose C4 |
| 82 | 2-methyl-3-furanthiol | 114 | A; B | 0: 52 / 57; 5: 48 / 43 | glucose C5 (1-deoxypentosone + H2S) |

**Pyridines are absent from Table 2**: no labelling data for 2-butylpyridine or 2-pentylpyridine.

### Mechanisms (Fig. 1, Fig. 2: FIGURE-ONLY; the text's clauses)

- Fig. 1(a)(b): "2-alkylthiophene could be formed from the reaction of a 2-enal or 2,4-dienal with
  H2S (Mottram, 1998)"; the alkylthiophenes Nos. 21, 23, 27, 28, 32, 34, 35, 42, 46 "were produced
  from aliphatic aldehydes through cleavage, condensation with H2S, cyclization, or dehydration".
  Which carbons become which ring atoms is not stated in the text.
- Fig. 1(c): 2-alkylfuran from the (E)-2-alkenal (Adams et al. 2011, amino-acid-catalysed).
- Fig. 1(d): alkylthiophenes from thiophenes (alkylation; the "with ethanol under alkaline conditions"
  remark for Nos. 51 and 53).
- Fig. 1(e)-(h), Fig. 2: Nos. 38, 42, 43, 47 and Nos. 12, 26, 34, 38, 45 from glucose fragments +
  acetaldehyde / crotonaldehyde / formaldehyde + H2S.
- 2-Heptylthiophene (No. 36): "derived from undecenal, and the retro-aldol condensation of nonanal and
  acetaldehyde continued to react with H2S" (i.e. nonanal + acetaldehyde -> 2-undecenal -> + H2S).
- H2S and NH3 source: "the Strecker degradation of cysteine generates H2S, NH3, acetaldehyde,
  formaldehyde, and mercaptoacetaldehyde"; aldehydes "competed to react with cysteine, which weakened
  Strecker degradation to release H2S".
- No pyridine mechanism is given anywhere in the paper.

## 4. Routes and numbers the repository can use

Conditions for every row: Cys 200 mM + Glc 200 mM + aldehyde 10 mM, pH 6.5 phosphate 0.2 M, 150 C,
120 min, HS-SPME, one IS with response factor 1, n = 3, µg/L of the 5 mL pot.

| route | reactant -> product | mechanism as drawn / stated | measured numbers (units, conditions) | evidence class |
|---|---|---|---|---|
| DU-24D-THIO | **(E,E)-2,4-nonadienal + H2S -> 2-pentylthiophene** | Fig. 1(a)(b), FIGURE-ONLY; text: 2,4-dienal + H2S, cyclisation, dehydration (Mottram 1998). Carbon bookkeeping (ours): C9 in, C9 out; S bridges dienal C1 and C4; pentyl = C5-C9; no oxidation needed (C9H14O + H2S -> C9H14S + H2O) | 238.65 ± 23.37 µg/L (C); not detected in M, A, B, D; CAMOLA 92 % unlabelled, 8 % one 13C | level_only; within_study_ratio; label-supported carbon source |
| DU-24D-PYR | **(E,E)-2,4-nonadienal + NH3 (from Cys) -> 2-butylpyridine** (the C9 analogue of 2,4-decadienal -> 2-pentylpyridine) | none given; by analogy the Zhou 2000 / Buttery route | 191.24 ± 1.49 µg/L (C); 52.32 ± 0.03 (B, from 2-nonenal); 12.93 ± 2.81 (D, from 2,6-nonadienal); none in M, A; RI/MS only, no standard | level_only; within_study_ratio (dienal : enal : 2,6-dienal = 15 : 4 : 1) |
| DU-2EN-PYR | (E)-2-nonenal -> 2-pentylpyridine (C10 from a C9 enal: one extra carbon) | none given | 44.15 ± 9.51 µg/L (B) with authentic standard; not detected in A, C, D, M | level_only; unexplained (flag 4) |
| DU-ENAL-NOT-THIO | (E)-2-nonenal + H2S -> 2-pentylthiophene | the text allows "a 2-enal or 2,4-dienal"; the data do not | **not detected** in B (would need an oxidation: C9H16O + H2S -> C9H14S + H2O + H2) | null result (level_only) |
| DU-HEX-THIO | C9 aldehyde + C1 (not from glucose) + H2S -> 2-hexylthiophene | not drawn; 99 % unlabelled in A, B, C, so the tenth carbon is from cysteine's Strecker formaldehyde or from the aldehyde pool | 10.77 ± 0.40 (A), 2.58 ± 0.31 (B), 6.38 ± 0.15 (C) µg/L; none in D, M | level_only; label-constrained |
| DU-HEPT-THIO | nonanal + acetaldehyde -> 2-undecenal; + H2S -> 2-heptylthiophene | stated (text); the C2 unit is 84 % non-glucose (cysteine Strecker acetaldehyde), 16 % glucose | 471.07 ± 17.11 µg/L (A only); the largest single aldehyde-derived sulfur volatile in the paper | level_only; label-supported |
| DU-24D-FURAN | 2-alkenal / dienal -> 2-alkylfuran (Adams 2011) | Fig. 1(c) | 2-pentylfuran 366.97 ± 69.66 (B, from 2-nonenal), 7.79 ± 0.40 (C); (E)-2-(2-pentenyl)furan 133.01 ± 4.54 (D), 7.19 (C); 2-butylfuran 8.45 (C); 2-propylfuran 4.35 (C); all 99 % unlabelled | level_only |
| DU-SUPPRESS | any C9 aldehyde suppresses the Cys/Glc sulfur slate | competition of the aldehyde with glucose for cysteine (Schiff base, thiazolidine); less Strecker H2S | thiophene sum M 2067 -> A 1181, B 537, C 816, D 583 µg/L; thieno[3,2-b]thiophene 1038 -> 420 / 147 / 0 / 109; 2-methyl-3-furanthiol 12.0 -> 11.4 / 8.3 / 0 / 0; 2-furfurylthiol 6.2 -> 0 in all; total volatiles 2580 -> 1572 / 1334 / 1075 / 959 | within_study_ratio |
| DU-THIAZOLE | thiazoles per system | not drawn; "degradation products of aliphatic aldehydes involved" | sums M 62.74, A 13.58, B 32.98, C 193.40, D 38.37 µg/L; C's 193 is thiazole 82.67 + 2-ethylthiazole 93.35 + 2,4,5-trimethyl 14.33 + isothiazole 3.04; **no 2-alkyl-4-methylthiazole with a chain > C2 in any system** | level_only |
| DU-25-DIALKYL | 2,4-nonadienal -> 2,5-dipropyl-, 2,5-dibutyl-, 2-ethyl-5-propyl-, 2-butyl-5-ethyl-, 2-ethyl-5-pentylthiophene | Fig. 1(d) alkylation of a thiophene; the 2,5-dibutyl is 38 % glucose-C4 | 99.69, 133.10, 25.70, 46.20, 37.78 µg/L (all C only, except 2,5-dibutyl 28.05 in M) | level_only |

Within-study ratios worth registering (system C, one pot):
- 2-pentylthiophene : 2-butylpyridine : 2-hexylthiophene : 2-pentylfuran = 238.65 : 191.24 : 6.38 :
  7.79 µg/L. With response factor 1 and SPME these are order-of-magnitude. On a molar basis (MW 154.3,
  135.2, 168.3, 138.2; ours) 1.55 : 1.41 : 0.038 : 0.056 µmol/L, i.e. the S and the N heterocycle from
  the same dienal are formed in about equal amount, and each is ~1.5 x 10^-4 of the 10 mM dienal charged
  (ours; SPME headspace levels, not a mass balance).
- 2-butylpyridine across aldehydes: C 191 : B 52 : D 13 : A 0 : M 0.
- 2-pentylthiophene across aldehydes: C 239 : B 0 : D 0 : A 0 : M 0. The enal and the alkanal are the
  negative controls for the dienal + H2S rule inside one study.

## 5. Rule sketches (repository suggestions, not the paper's)

Registry check (`data/keys/compounds.yml`): `hydrogen_sulfide` (alias h2s), `2_methylthiophene`,
`2_pentylfuran`, `hexanal`, `nonanal`, `2_pentyl_4_methylthiazole` (SMILES `CCCCCC1=NC(C)=CS1`),
`2_hexyl_4_methylthiazole`, `4_5_dihydro_2_methylthiazole` exist. **2-pentylthiophene, 2-hexylthiophene,
2-heptylthiophene, 2-butylpyridine and 2-pentylpyridine have no key.** Species: `Cys`, `H2S`,
`DECADIENAL` exist in `data/species/structures.yml`; no 2,4-nonadienal, no 2-nonenal, no NH3.

**S1. 2,4-alkadienal + H2S -> 2-alkylthiophene (net; DU-24D-THIO).** Reactant R-CH=CH-CH=CH-CHO.
Change: S bonds to C1 and C4; the C1-C4 chain plus S is the ring; C5 onward is the 2-alkyl; loses one
H2O; no oxidation. Atom map is ours from the carbon count (the paper's scheme is not readable in text);
no label experiment fixes it, unlike the pyridine.
- positive: 2,4-nonadienal `CCCC/C=C/C=C/C=O` + `S` -> `CCCCCc1cccs1` (2-pentylthiophene) + `O` (238.65 µg/L, system C)
- second positive (the engine's dienal; proposed by homology, not measured here): `CCCCC/C=C/C=C/C=O` (DECADIENAL) + `S` -> `CCCCCCc1cccs1` (2-hexylthiophene)
- negative: (E)-2-nonenal `CCCCCC/C=C/C=O` + `S` -> must NOT give 2-pentylthiophene (system B: none); nonanal `CCCCCCCCC=O` + `S` -> no fire; (E,Z)-2,6-nonadienal `CC/C=C\CC/C=C/C=O` + `S` -> no 2-pentylthiophene (system D: none; its 2-propylthiophene 105 µg/L and (E)-2-(2-pentenyl)furan 133 come from the non-conjugated diene by a different cut).
- required substructure: `O=CH-CH=CH-CH=CH-` conjugated; the 2,6-isomer must fail.

**S2. 2,4-alkadienal + NH3 -> 2-alkylpyridine (DU-24D-PYR; same rule as `zamora2020_extraction.md` S1
and `zhou2000_extraction.md` S1).**
- positive here: `CCCC/C=C/C=C/C=O` + `N` -> `CCCCc1ccccn1` (2-butylpyridine; 191.24 µg/L, C)
- negative here: nonanal + `N` -> no fire (A: no pyridine at all); the 2-nonenal and 2,6-nonadienal
  pots did give 2-butylpyridine (52, 13 µg/L), so the enal is not a clean negative for the N ring the
  way it is for the S ring; the rule's negative control should be the alkanal.

**S3. Cysteine -> H2S + NH3 + acetaldehyde (+ formaldehyde, mercaptoacetaldehyde) (the source of both
heteroatoms).** Already in the rule file as R08 (cysteine thermolysis); this paper is a second anchor
for "NH3 from cysteine" feeding a pyridine rule: `NC(CS)C(=O)O` -> `S` + `N` + `CC=O` (+ CO2). Negative:
`NCC(=O)O` (Gly) -> no H2S.

**S4. Nonanal + acetaldehyde -> 2-undecenal; 2-undecenal + H2S -> 2-heptylthiophene (net; DU-HEPT-THIO).**
- positive: `CCCCCCCCC=O` + `CC=O` -> `CCCCCCCC/C=C/C=O` (UNDECENAL_2E exists in literature_structures) ; `CCCCCCCC/C=C/C=O` + `S` -> `CCCCCCCc1cccs1` + `O` + H2 (an oxidation is needed for the enal, see S1 negative; the authors do not address this).
- Because the enal route needs an oxidation that the dienal route does not, and system B shows no
  2-pentylthiophene from 2-nonenal, S4's second step should be `status: proposed` and carry the
  contradiction as a note.

**S5. Alkylthiazoles.** No rule can be anchored here for `2_pentyl_4_methylthiazole` or
`2_hexyl_4_methylthiazole`: neither was detected. The thiazoles that did appear with the dienal
(thiazole 82.67, 2-ethylthiazole 93.35 µg/L) are C2-C5 bodies made from glucose/cysteine fragments
(authors: "degradation products of aliphatic aldehydes involved", not drawn). A 2-alkyl-4-methylthiazole
would need alkanal + mercaptoacetone (or 1-mercapto-2-propanone) + NH3; that is R27's mercaptoketone
plus an alkanal, which would be `proposed` with no source in this batch: `CCCCCC=O` (hexanal) +
`CC(=O)CS` + `N` -> `CCCCCc1nc(C)cs1` + 2 H2O + H2. Negative: hexanal + `S` + `N` without the
mercaptoketone -> no fire.

## 6. Flags

1. **Semi-quantitative:** one internal standard (1,2-dichlorobenzene), response factor 1 for all 88
   compounds, HS-SPME on a CAR/PDMS/DVB fibre, no calibration, no recovery. µg/L values are
   headspace-weighted levels; use as presence/absence and within-compound ratios across systems.
   Evidence class level_only throughout.
2. **Single time point (120 min), single temperature (150 C), single aldehyde loading (10 mM, 1 : 20
   to cysteine).** No rate can be extracted.
3. **Printed sums that do not match their rows:** pyridines in B (printed 68.75; rows 118.74) and
   naphthalenes in B (printed 63.22; rows 13.22). Thiazole, thiophene and furan sums reproduce. Table 2
   lists 1-(2-thienyl)ethanone for system C although Table 1 has no C entry, lists 2-propylfuran for A
   although Table 1 has no A entry, and prints No. 34 twice.
4. **2-Pentylpyridine (C10) from (E)-2-nonenal (C9) and not from 2,4-nonadienal** is unexplained by the
   paper and by any dienal + NH3 rule; the identity is standard-confirmed. A C1 addition (formaldehyde
   from cysteine Strecker) before ring closure is the obvious candidate; no label data (pyridines are
   absent from Table 2).
5. **Internal-standard amount is ambiguous** ("1 µL of 1,2-dichlorobenzene (100 µg, in 1 mL of
   methanol)"): 0.1 µg or 100 µg. A 1000-fold uncertainty in the absolute scale of every µg/L value;
   ratios are unaffected.
6. **No 2-alkyl-4-methylthiazole** (the roadmap's 2-pentyl- and 2-hexyl-4-methylthiazole) in any
   system, with a standard library that did not include them either (section 2.1 lists 2-methyl-,
   2-ethyl-, 2,4,5-trimethyl-, 4,5-dimethylthiazole). Absence is therefore "not identified", not
   "not formed".
7. **Fig. 1 and Fig. 2 are unreadable in the text layer**; the atom maps in §5 are ours from carbon
   counts and the CAMOLA fractions, not the authors' arrows.
8. **The 2-enal + H2S -> 2-alkylthiophene statement in the text is contradicted by the data for the
   C9 enal** (no 2-pentylthiophene in B) while 2-hexylthiophene (C10) and 2-heptylthiophene (C11) do
   appear from the C9 alkanal/enal: the alkyl thiophenes from saturated/mono-unsaturated aldehydes go
   through a chain extension (aldol with a C1/C2 Strecker aldehyde) first. A rule "enal + H2S ->
   thiophene of the same carbon number" is not supported here.
9. **Aldehyde solubility:** 10 mM of a C9 aldehyde in 5 mL water at 150 C in a sealed vial is a
   two-phase system; the authors do not comment.
10. **Ammonia is never measured or named as a product**; the pyridines are the only evidence that NH3
    from cysteine reached the dienal. The paper also never names 2-butylpyridine's mechanism.
11. **Table S1** (acids, alcohols, ketones, residual aldehydes per system) was not available; the
    remaining aldehyde after 120 min is unknown, so even an order-of-magnitude conversion of dienal to
    thiophene (ours, ~1.5 x 10^-4 on the headspace number) has no mass-balance check.

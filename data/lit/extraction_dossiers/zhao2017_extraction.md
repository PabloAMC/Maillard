# Zhao et al. 2017 — EXTRACTION
### *"Volatile flavor constituents in the pork broth of black-pig"*
### Yunan black-pig hindquarter, 200 g lean cubes + 200 mL water, 3-neck flask under reflux, **~100 °C for 3 h** (oil bath). SAFE + GC-MS/MS (SRM) against one internal standard (1,2-dichlorobenzene) with external calibration in CH2Cl2, **no recovery correction**, µg/kg **wet meat**. **MFT 3.39, methional 6.65, DMTS 0.89, hexanal 26.7 µg/kg; DMDS 0.36 (semi-quant, response factor 1); FFT detected by GC-O (log2FD 17) but below S/N 10, not quantified.** Plus meat fatty acids and acid-hydrolysate amino acids. Read 2026-10-28.

**Full read of `data/articles/zhao2017.pdf`** (10 pp., typeset). Text layer clean apart from the Elsevier
font substitutions ("lg" = µg, "c-decalactone" = γ-decalactone, "6" = ≤, "P" = ≥, "" = ×). Tables 1-4
transcribed from `pdftotext -layout`. Figure 1 (spider web of log2FD ratios) not digitised.
Supplementary Tables S1-S3 (common white-pig comparison, Wang et al. 2015) not on disk.

**Provenance:** **[M]** printed number · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value |
|---|---|
| Authors | Jian Zhao, Meng Wang, Jianchun Xie*, Mengyao Zhao, Li Hou, Jingjing Liang, Shi Wang, Jie Cheng |
| Title | "Volatile flavor constituents in the pork broth of black-pig" (first page matches) |
| Venue | *Food Chemistry* **226 (2017) 51–60** (printed; Crossref agrees, issue date July 2017) |
| DOI | **10.1016/j.foodchem.2017.01.011** (printed) |
| Dates | received 2 Aug 2016; revised 17 Nov 2016; accepted 3 Jan 2017; online 4 Jan 2017 |
| Affiliation | Beijing Technology and Business University (Xie group); CAAS |
| Companions | Wang et al. 2015 (Food Sci. (China) 36:209, common white-pig broth; not on disk); Fan 2018 and Du 2020 use the same lab method with recovery factors (`fan2018_extraction.md`, `du2020_extraction.md`) |

---

## 1. THE ONE-PARAGRAPH ANSWER

A 3-h reflux pork broth from lean black-pig hindquarter contains, per kg of meat stewed, **MFT 3.39 µg/kg,
methional 6.65 µg/kg, dimethyl trisulfide 0.89 µg/kg, hexanal 26.7 µg/kg, 2-acetylthiazole 2.2 µg/kg**,
quantified by GC-MS/MS-SRM on the SAFE concentrate. **FFT was strongly odour-active (log2FD 17) but below
S/N 10 in SRM and was not quantified; DMDS appears only in the full-scan table at 0.36 µg/kg on a response
factor of 1.** By AEDA the top odorants are MFT and methional (log2FD 18) and γ-decalactone (18), then FFT
(17). The meat is lean (1.81 % fat; MUFA 58 %, PUFA 6.4 %). No recovery correction is applied, so the thiol
values are likely low (sec. 4).

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Meat | hindquarters of 8 Yunan black pigs (Huainan × Duroc; two slaughter days), Wu-Mart Beijing; skin, visible fat and connective tissue removed; cubes ~0.5 cm³; "stored in a refrigerator at 20 °C for 72 h maximum" [!] (sign lost; −20 °C) |
| Composition | crude protein **22.85 %** (wet; AOAC 1990); crude fat **1.81 %** (wet; Soxhlet, GB/T 5009.6-2003) |
| Broth | "**Two hundred grams of the meat and 200 mL of water** were placed in a 1000 mL 3-neck flask, which was fitted with a reflux condenser and a mechanical stirrer. The meat was stewed at **ca. 100 °C for 3 h** by an oil bath." Three replicates. No salt, no pH, no final broth mass [!] |
| Isolation | meat residue filtered off; broth extracted 3 × 200 mL CH2Cl2; SAFE at 40 °C (10⁻⁴-10⁻⁵ Pa); Na2SO4; Vigreux to ~1 mL, N2 to **0.35 mL** |
| Full-scan GC-MS | DB-Wax and DB-5MS; "quantity ... approximately calculated by its peak area to that of 1,2-dichlorobenzene (internal standard, **200 µg/mL** in dichloromethane) using a **calibration factor of 1**", converted to µg/kg wet meat by the concentrate yield (Table 3) |
| GC-O | AEDA on DB-5 (Table 4 footnote says HP-5), 1:2 dilution series, three sniffers |
| Targeted quantitation | Trace 1310 / TSQ 8000 GC-MS/MS, **SRM**, TG-5MS; odorants confirmed by standards and S/N > 10; calibration = relative peak area vs 1,2-dichlorobenzene in CH2Cl2 (Table 1); converted to µg/kg wet meat "according to the yield of the SAFE concentrate relative to the meat used". **No recovery factor, no isotope standard** |
| Meat analyses | FAME by GC-MS (C13:0 IS); amino acids after **6 N HCl, 110 °C, 24 h** (total, not free), ion-exchange/ninhydrin, mg/g dry meat |
| Replicates | n = 3 broths |

---

## 3. TABLES (verbatim)

### Table 1. GC-MS/MS-SRM quantitation parameters (x = area ratio to 1,2-dichlorobenzene; y = µg/mL; standards in CH2Cl2)

| compound | RT (min) | SRM, parent > daughter (CE, eV) | calibration | R² | standards (µg/mL) |
|---|---|---|---|---|---|
| pentanal | 3.63 | 58 > 29 (10) | y = 9.7724 x | 0.9948 | 1.0, 5, 10, 25, 50 |
| 3-hydroxy-2-butanone | 3.84 | 73 > 45 (5); 88 > 45 (5) | y = 17.727 x | 0.9996 | 5, 10, 25, 50, 500 |
| **hexanal** | 5.56 | 56 > 41 (10); 72 > 57 (10) | y = 7.6846 x | 0.9959 | 1.0, 5, 10, 25, 50 |
| (E)-2-hexenal | 6.94 | 69 > 41 (10); 83 > 55 (5) | y = 7.4210 x | 0.9985 | 0.1, 0.5, 1.0, 2.0, 5 |
| **2-methyl-3-furanthiol** | 7.35 | 85 > 45 (5); 114 > 85 (5) | y = 1.8746 x | 0.9962 | 0.1, 0.5, 1, 5, 10 |
| heptanal | 8.61 | 55 > 29 (10); 70 > 55 (5) | y = 2.4138 x | 0.9975 | 0.1, 0.5, 1, 5, 10 |
| **3-(methylthio)propanal** | 8.75 | 76 > 48 (5); 104 > 48 (10) | y = 3.6742 x | 0.9994 | 0.1, 0.5, 1, 5, 10 |
| (E)-2-heptenal | 11.23 | 55 > 29 (10); 83 > 55 (5) | y = 3.8603 x | 0.9988 | 0.1, 0.5, 1, 5, 10 |
| **dimethyl trisulfide** | 11.81 | 79 > 64 (15); 126 > 79 (15) | y = 1.2189 x | 0.9979 | 0.1, 0.5, 1.0, 2.0, 5 |
| 1-octen-3-ol | 12.74 | 57 > 29 (10); 72 > 43 (5) | y = 4.0858 x | 0.9965 | 1.0, 5, 10, 25, 50 |
| 2-thiophenecarboxaldehyde | 13.88 | 82.9 > 39 (10); 111 > 83 (10) | y = 20.721 x | 0.9985 | 0.1, 0.5, 1.0, 2.0, 5 |
| 2-ethyl-3-methylpyrazine | 14.22 | 94 > 67 (10); 121 > 94 (10) | y = 2.7390 x | 0.9989 | 0.1, 0.5, 1.0, 2.0, 5 |
| 2-acetylthiazole | 15.29 | 99 > 58 (20); 127 > 99 (5) | y = 8.4361 x | 0.9990 | 0.1, 0.5, 1, 5, 10 |
| phenylacetaldehyde | 17.45 | 91 > 65 (10); 120 > 91 (10) | y = 1.3696 x | 0.9981 | 0.1, 0.5, 1.0, 2.0, 5 |
| nonanal | 24.30 | 57 > 29 (10); 70 > 55 (5) | y = 8.9647 x | 0.9940 | 0.1, 0.5, 1, 5, 10 |
| (E)-2-nonenal | 31.85 | 83 > 55 (5) | y = 4.3039 x | 0.9918 | 0.1, 0.5, 1.0, 2.0, 5 |
| (E)-2-decenal | 39.72 | 55 > 29 (10); 70 > 41 (10) | y = 15.937 x | 0.9902 | 0.1, 0.5, 1.0, 2.0, 5 |
| 4-methyl-5-thiazoleethanol | 39.94 | 112 > 85 (5); 143 > 112 (10) | y = 1.5534 x | 0.9991 | 0.1, 0.5, 1.0, 5, 10 |
| (E,E)-2,4-decadienal | 40.45 | 67 > 41 (10); 81 > 53 (15) | y = 2.6367 x | 0.9917 | 0.1, 0.5, 1.0, 2.0, 5 |
| 1,2-dichlorobenzene | 16.39 | 111 > 75 (10); 146 > 111 (15) | internal standard | | **0.5** |

### Table 2. Meat of the black-pig hindquarter: fatty acids and amino acids (mean ± SD, n = 3)

| fatty acid | mg/g dry meat (as printed) | % total | | amino acid | mg/g dry meat | % total |
|---|---|---|---|---|---|---|
| C14:0 | 3.05 ± 0.26 | 5.83 ± 0.17 | | Valine | 35.15 ± 0.24 | 4.86 ± 0.02 |
| C16:0 | 11.35 ± 0.10 | 21.69 ± 1.07 | | Isoleucine | 32.91 ± 0.16 | 4.55 ± 0.03 |
| C18:0 | 4.12 ± 0.37 | 7.87 ± 0.25 | | Leucine | 63.06 ± 0.37 | 8.73 ± 0.04 |
| C20:0 | 0.13 ± 0.01 | 0.25 ± 0.00 | | Phenylalanine | 32.60 ± 0.63 | 4.51 ± 0.04 |
| C16:1 | 3.12 ± 0.02 | 5.96 ± 0.34 | | Histidine | 32.37 ± 0.32 | 4.48 ± 0.04 |
| C18:1 | 22.56 ± 1.86 | 43.11 ± 1.07 | | Tyrosine | 27.58 ± 0.99 | 3.82 ± 0.11 |
| C20:1 | 4.63 ± 0.05 | 8.85 ± 0.44 | | Tryptophan | 5.49 ± 0.13 | 0.76 ± 0.02 |
| C22:1 | 0.01 ± 0.00 | 0.02 ± 0.00 | | Arginine | 49.87 ± 0.39 | 6.90 ± 0.03 |
| C14:2 | 0.02 ± 0.00 | 0.04 ± 0.00 | | **Methionine** | **17.98 ± 0.14** | 2.49 ± 0.02 |
| C18:2 | 2.89 ± 0.30 | 5.52 ± 0.25 | | Glutamic acid | 130.63 ± 1.84 | 18.07 ± 0.07 |
| C20:2 | 0.14 ± 0.02 | 0.27 ± 0.02 | | Aspartic acid | 74.65 ± 0.94 | 10.33 ± 0.05 |
| C20:3 | 0.01 ± 0.01 | 0.02 ± 0.01 | | Threonine | 36.33 ± 0.30 | 5.03 ± 0.07 |
| C20:4 | 0.20 ± 0.02 | 0.38 ± 0.01 | | Serine | 32.00 ± 0.27 | 4.43 ± 0.02 |
| C22:3 | 0.03 ± 0.01 | 0.06 ± 0.01 | | Glycine | 34.50 ± 0.14 | 4.77 ± 0.03 |
| C22:4 | 0.06 ± 0.00 | 0.11 ± 0.01 | | Alanine | 44.48 ± 0.59 | 6.15 ± 0.02 |
| C22:6 | 0.01 ± 0.00 | 0.02 ± 0.01 | | Lysine | 38.50 ± 0.21 | 5.33 ± 0.04 |
| SFA | 18.65 | 35.64 | | **Cysteine** | **4.97 ± 0.14** | 0.69 ± 0.01 |
| MUFA | 30.32 | 57.94 | | Proline | 29.63 ± 0.90 | 4.10 ± 0.10 |
| PUFA | 3.36 | 6.42 | | Umami (Glu + Asp) | 205.28 | 28.40 |
| Total | 52.33 | 100 | | Total | 722.71 | 100 |

[!] The text calls the fatty-acid total "52.33 mg/g meat" while the header says "mg/g dry meat"; with 1.81 %
fat on a wet basis, 52.33 mg/g is plausible only on a dry basis (~28 % dry matter -> ~65 mg fat/g dry; FAs
~80 % of fat). Read as dry basis.

### Table 3. Full-scan GC-MS of the SAFE extract (µg/kg wet meat, mean ± SD, n = 3, **semi-quantitative: 1,2-dichlorobenzene, calibration factor 1**; "–" undetected; DB-Wax / DB-5 / averaged)

| compound | RI Wax | RI DB-5 | DB-Wax | DB-5 | averaged | ID |
|---|---|---|---|---|---|---|
| *Sulfur compounds* | | | | | | |
| **Dimethyl disulfide** | 1073 | 737 | 0.23 ± 0.03 | 0.49 ± 0.02 | **0.36 ± 0.02** | MS/RI/S |
| 3-(Methylthio)propanal | 1424 | 906 | 1.16 ± 0.03 | 1.12 ± 0.04 | 1.14 ± 0.04 | MS/RI/S |
| 2-(Methylthio)ethanol | 1523 | 823 | 0.31 ± 0.01 | 0.65 ± 0.05 | 0.48 ± 0.03 | MS/RI |
| 2-Acetylthiazole | 1650 | 1040 | 1.78 ± 0.03 | 4.19 ± 0.37 | 2.98 ± 0.20 | MS/RI/S |
| 2-Thiophenecarboxaldehyde | 1679 | – | 0.65 ± 0.03 | – | 0.33 ± 0.02 | MS/RI/S |
| Benzothiazole | 1961 | – | 0.69 ± 0.01 | – | 0.34 ± 0.01 | MS/RI/S |
| 4-Methyl-5-thiazoleethanol | 2275 | 1266 | 6.41 ± 0.97 | 5.69 ± 0.27 | 6.05 ± 0.62 | MS/RI/S |
| Subtotal | | | 11.23 | 12.14 | 11.69 | |
| *Heterocyclic compounds* | | | | | | |
| Furfural | – | 818 | – | 0.51 ± 0.25 | 0.26 ± 0.12 | MS/RI/S |
| Pyridine | 1178 | 752 | 1.69 ± 0.17 | 2.12 ± 0.22 | 1.91 ± 0.20 | MS/RI |
| Pyrazine | 1206 | – | 0.51 ± 0.06 | – | 0.26 ± 0.03 | MS/RI |
| 2-Pentylfuran | 1234 | 993 | 1.21 ± 0.01 | 2.71 ± 0.08 | 1.96 ± 0.05 | MS/RI/S |
| 2-Ethyl-3-methylpyrazine | 1356 | 981 | 8.42 ± 0.20 | 3.64 ± 0.28 | 6.03 ± 0.24 | MS/RI/S |
| Furfuryl alcohol | 1635 | 864 | 1.10 ± 0.17 | 0.29 ± 0.01 | 0.69 ± 0.09 | MS/RI/S |
| 1-Methyl-2-pyrrolidinone | 1660 | 1033 | 1.73 ± 0.13 | 2.71 ± 0.03 | 2.22 ± 0.08 | MS/RI |
| Subtotal | | | 14.66 | 11.98 | 13.32 | |
| *Aldehydes* | | | | | | |
| Pentanal | 1002 | 693 | 3.11 ± 0.07 | 5.01 ± 0.09 | 4.06 ± 0.08 | MS/RI/S |
| **Hexanal** | 1098 | 789 | 22.65 ± 1.27 | 43.28 ± 1.20 | 32.97 ± 1.24 | MS/RI/S |
| 2-Methyl-2-butenal | 1103 | 738 | 3.27 ± 0.09 | 0.54 ± 0.04 | 1.90 ± 0.06 | MS/RI/S |
| Heptanal | 1192 | 881 | 0.52 ± 0.05 | 1.34 ± 0.05 | 0.93 ± 0.05 | MS/RI/S |
| (E)-2-hexenal | 1206 | 856 | 0.15 ± 0.01 | 0.23 ± 0.01 | 0.19 ± 0.01 | MS/RI/S |
| (E)-2-heptenal | 1316 | 970 | 1.30 ± 0.08 | 0.69 ± 0.07 | 0.99 ± 0.08 | MS/RI/S |
| Nonanal | 1388 | 1090 | 6.28 ± 0.21 | 11.66 ± 0.39 | 8.97 ± 0.30 | MS/RI/S |
| (E)-2-octenal | 1416 | 1065 | 1.04 ± 0.08 | 1.05 ± 0.02 | 1.05 ± 0.05 | MS/RI/S |
| Benzaldehyde | 1530 | 939 | 7.16 ± 0.13 | 15.93 ± 0.68 | 11.55 ± 0.40 | MS/RI/S |
| (E)-2-nonenal | – | 1157 | – | 0.64 ± 0.03 | 0.32 ± 0.02 | MS/RI/S |
| Phenylacetaldehyde | 1648 | 1031 | 0.93 ± 0.04 | 2.04 ± 0.01 | 1.48 ± 0.02 | MS,RI/S |
| (E)-2-decenal | – | 1247 | – | 0.49 ± 0.03 | 0.24 ± 0.02 | MS/RI/S |
| (E,E)-2,4-nonadienal | – | 1204 | – | 0.11 ± 0.00 | 0.05 ± 0.00 | MS/RI/S |
| (E,E)-2,4-decadienal | 1781 | 1306 | 0.76 ± 0.14 | 1.46 ± 0.05 | 1.11 ± 0.10 | MS/RI/S |
| Tetradecanal | – | 1609 | – | 0.62 ± 0.04 | 0.31 ± 0.02 | MS/RI |
| Hexadecanal | 2133 | 1791 | 6.90 ± 0.49 | 12.46 ± 0.52 | 9.68 ± 0.50 | MS/RI |
| Octadecanal | 2380 | 1998 | 8.62 ± 0.45 | 5.57 ± 0.43 | 7.09 ± 0.44 | MS/RI |
| Subtotal | | | 62.69 | 103.12 | 82.91 | |
| *Ketones* | | | | | | |
| 2-Heptanone | – | 892 | – | 1.12 ± 0.04 | 0.56 ± 0.02 | MS/RI |
| 2,3-Octanedione | 1333 | 983 | 0.38 ± 0.01 | 1.03 ± 0.02 | 0.71 ± 0.01 | MS/RI |
| 3-Methyl-2-cyclopenten-1-one | 1498 | 976 | 1.88 ± 1.29 | 0.77 ± 0.05 | 1.33 ± 0.67 | MS/RI |
| 2-Pentadecanone | 2008 | 1688 | 1.85 ± 0.03 | 0.89 ± 0.05 | 1.37 ± 0.04 | MS/RI |
| Subtotal | | | 4.11 | 3.81 | 3.96 | |
| *Alcohols* | | | | | | |
| 1-Butanol | 1155 | – | 1.23 ± 0.09 | – | 0.62 ± 0.04 | MS/RI |
| 1-Penten-3-ol | 1172 | – | 0.27 ± 0.01 | – | 0.14 ± 0.00 | MS/RI/S |
| 1-Pentanol | 1250 | 761 | 1.67 ± 0.08 | 8.62 ± 1.44 | 5.14 ± 0.76 | MS/RI/S |
| 3-Methyl-3-buten-1-ol | 1254 | 733 | 0.20 ± 0.01 | 0.80 ± 0.03 | 0.50 ± 0.02 | MS/RI |
| 3-Methyl-2-buten-1-ol | 1325 | 767 | 0.64 ± 0.05 | 1.20 ± 0.03 | 0.92 ± 0.04 | MS/RI |
| 1-Hexanol | 1358 | – | 0.81 ± 0.38 | – | 0.40 ± 0.19 | MS/RI/S |
| 1-Octen-3-ol | 1451 | 960 | 4.39 ± 0.11 | 8.04 ± 0.14 | 6.22 ± 0.12 | MS/RI/S |
| 1-Heptanol | – | 970 | – | 0.64 ± 0.03 | 0.32 ± 0.02 | MS/RI/S |
| 2-Ethyl-1-hexanol | 1493 | – | 0.90 ± 0.01 | – | 0.45 ± 0.00 | MS/RI |
| 1-Octanol | 1553 | 1072 | 2.62 ± 0.03 | 5.84 ± 0.20 | 4.23 ± 0.12 | MS/RI/S |
| 2,3-Butanediol | 1568 | 778 | 3.02 ± 0.04 | 10.85 ± 0.38 | 6.93 ± 0.21 | MS/RI |
| 1-Phenethyl alcohol | 1821 | 1039 | 1.07 ± 0.01 | 3.00 ± 0.08 | 2.04 ± 0.05 | MS/RI |
| Benzyl alcohol | 1848 | 1045 | 1.30 ± 0.08 | 3.45 ± 0.47 | 2.37 ± 0.27 | MS/RI |
| Phenylethyl alcohol | 1883 | 1112 | 1.95 ± 0.03 | 0.61 ± 0.00 | 1.28 ± 0.02 | MS/RI/S |
| Pentadecanol | – | 1765 | – | 1.64 ± 0.12 | 0.82 ± 0.06 | MS/RI |
| Subtotal | | | 20.07 | 44.69 | 32.38 | |
| *Acids* | | | | | | |
| Acetic acid | 1431 | – | 4.10 ± 0.05 | – | 2.05 ± 0.03 | MS/RI/S |
| Propanoic acid | 1539 | – | 2.27 ± 1.82 | – | 1.14 ± 0.91 | MS/RI |
| Butanoic acid | 1635 | – | 13.99 ± 7.60 | – | 6.99 ± 3.80 | MS/RI |
| Pentanoic acid | 1732 | – | 11.59 ± 9.41 | – | 5.80 ± 4.71 | MS/RI |
| Hexanoic acid | 1857 | 1036 | 5.60 ± 0.70 | 0.57 ± 0.00 | 3.08 ± 0.35 | MS/RI |
| Heptanoic acid | 1952 | 1071 | 0.90 ± 0.02 | 0.81 ± 0.06 | 0.86 ± 0.04 | MS/RI |
| Octanoic acid | 2049 | – | 2.39 ± 0.11 | – | 1.19 ± 0.05 | MS/RI |
| Nonanoic acid | 2153 | – | 2.55 ± 0.07 | – | 1.27 ± 0.04 | MS/RI |
| Decanoic acid | – | 1376 | – | 1.24 ± 0.00 | 0.62 ± 0.00 | MS/RI |
| Dodecanoic acid | – | 1574 | – | 3.08 ± 0.04 | 1.54 ± 0.02 | MS/RI |
| Tetradecanoic acid | – | 1769 | – | 9.34 ± 0.66 | 4.67 ± 0.33 | MS/RI |
| Pentadecanoic acid | – | 1856 | – | 3.08 ± 0.11 | 1.54 ± 0.06 | MS/RI |
| (Z)-9-hexadecenoic acid | – | 1948 | – | 3.21 ± 0.16 | 1.60 ± 0.08 | MS/RI |
| (Z)-11-hexadecenoic acid | – | 1958 | – | 7.42 ± 0.19 | 3.71 ± 0.09 | MS/RI |
| Hexadecanoic acid | – | 1952 | – | 33.87 ± 0.32 | 16.93 ± 0.16 | MS/RI |
| Heptadecanoic acid | – | 2068 | – | 10.19 ± 0.45 | 5.10 ± 0.23 | MS/RI |
| (E)-9-octadecenoic acid | – | 2125 | – | 42.54 ± 0.84 | 21.27 ± 0.42 | MS/RI |
| Octadecanoic acid | – | 2148 | – | 32.79 ± 0.53 | 16.40 ± 0.27 | MS/RI |
| Subtotal | | | 43.39 | 148.14 | 95.77 | |
| *Esters* | | | | | | |
| Acetic acid, butyl ester | 1089 | – | 0.13 ± 0.00 | – | 0.07 ± 0.00 | MS/RI |
| Formic acid, octyl ester | 1544 | – | 0.71 ± 0.06 | – | 0.35 ± 0.03 | MS/RI |
| γ-Butyrolactone | 1610 | 908 | 0.15 ± 0.00 | 2.01 ± 0.01 | 1.08 ± 0.01 | MS/RI |
| γ-Hexalactone | 1693 | 1063 | 1.26 ± 0.00 | 0.41 ± 0.04 | 0.84 ± 0.02 | MS/RI |
| δ-Hexalactone | 1800 | 1068 | 0.99 ± 0.05 | 1.45 ± 0.08 | 1.22 ± 0.07 | MS/RI |
| Hexadecanoic acid, methyl ester | 2196 | 1923 | 1.28 ± 0.02 | 3.82 ± 0.23 | 2.55 ± 0.13 | MS/RI |
| 15-Octadecenoic acid, methyl ester | 2379 | – | 3.11 ± 0.07 | – | 1.56 ± 0.03 | MS |
| Subtotal | | | 7.63 | 7.69 | 7.66 | |
| *Aliphatic hydrocarbons* | | | | | | |
| 2,3,5-Trimethyl-hexane | – | 803 | – | 3.29 ± 0.15 | 1.64 ± 0.07 | MS/RI |
| 2,4-Dimethyl-heptane | 823 | 844 | 0.49 ± 0.21 | 0.29 ± 0.01 | 0.39 ± 0.11 | MS/RI |
| 2,4-Dimethyl-1-heptene | – | 851 | – | 0.71 ± 0.06 | 0.36 ± 0.03 | MS/RI |
| 4-Methyl-octane | – | 863 | – | 4.03 ± 0.29 | 2.01 ± 0.15 | MS/RI |
| Nonane | – | 896 | – | 1.14 ± 0.05 | 0.57 ± 0.02 | MS/RI |
| Decane | – | 993 | – | 3.45 ± 0.26 | 1.73 ± 0.13 | MS/RI |
| 3,7-Dimethyl-decane | 1085 | 1100 | 2.08 ± 0.04 | 3.72 ± 0.21 | 2.90 ± 0.12 | MS/RI |
| Dodecane | 1189 | 1203 | 0.19 ± 0.01 | 6.93 ± 0.15 | 3.56 ± 0.08 | MS/RI |
| Tridecane | 1301 | 1293 | 1.06 ± 0.07 | 2.26 ± 0.07 | 1.66 ± 0.07 | MS/RI |
| Tetradecane | 1399 | 1397 | 0.23 ± 0.00 | 0.84 ± 0.04 | 0.53 ± 0.02 | MS/RI |
| Pentadecane | 1495 | 1498 | 0.33 ± 0.01 | 0.82 ± 0.01 | 0.57 ± 0.01 | MS/RI |
| Hexadecane | – | 1605 | – | 0.80 ± 0.02 | 0.40 ± 0.01 | MS/RI |
| 1-Heptadecene | 1759 | 1685 | 0.23 ± 0.00 | 1.48 ± 0.07 | 0.85 ± 0.04 | MS/RI |
| 2,6,10,14-Tetramethyl-pentadecane | – | 1689 | – | 0.76 ± 0.03 | 0.38 ± 0.01 | MS/RI |
| Octadecane | – | 1792 | – | 0.83 ± 0.03 | 0.42 ± 0.01 | MS/RI |
| Nonadecane | – | 1896 | – | 2.18 ± 0.09 | 1.09 ± 0.05 | MS,RI |
| Eicosane | – | 1998 | – | 3.83 ± 0.07 | 1.92 ± 0.04 | MS/RI |
| Docosane | 2193 | 2199 | 7.09 ± 0.12 | 19.91 ± 0.90 | 13.50 ± 0.51 | MS/RI |
| Tricosane | 2287 | 2289 | 11.42 ± 0.00 | 32.30 ± 0.56 | 21.86 ± 0.28 | MS/RI |
| Tetracosane | 2405 | 2397 | 12.74 ± 0.29 | 29.94 ± 0.71 | 21.34 ± 0.50 | MS/RI |
| Subtotal | | | 35.86 | 119.51 | 77.68 | |
| *Others* | | | | | | |
| p-Xylene | 1118 | 875 | 0.36 ± 0.01 | 1.47 ± 0.08 | 0.91 ± 0.04 | MS/RI |
| Ethylbenzene | 1132 | – | 0.40 ± 0.01 | – | 0.20 ± 0.01 | MS/RI |
| D-limonene | 1170 | 1014 | 1.45 ± 0.00 | 3.09 ± 0.09 | 2.27 ± 0.05 | MS/RI |
| o-Xylene | 1176 | – | 0.42 ± 0.01 | – | 0.21 ± 0.01 | MS/RI |
| Styrene | – | 895 | – | 0.34 ± 0.11 | 0.17 ± 0.06 | MS/RI |
| 1,2,3,4-Tetramethyl-benzene | – | 1144 | – | 0.62 ± 0.01 | 0.31 ± 0.00 | MS/RI |
| Butylated hydroxytoluene | 1902 | 1484 | 0.88 ± 0.07 | 0.54 ± 0.19 | 0.71 ± 0.13 | MS/RI |
| p-Cresol | 2071 | 1093 | 1.17 ± 0.05 | 2.59 ± 0.13 | 1.88 ± 0.09 | MS/RI |
| 2,4-Bis(1,1-dimethylethyl)-phenol | 2248 | 1487 | 1.07 ± 0.05 | 2.67 ± 0.09 | 1.87 ± 0.07 | MS/RI |
| Subtotal | | | 5.75 | 11.32 | 8.53 | |
| **Total** | | | 205.39 | 462.40 | 333.90 | |

[!] Where a compound is undetected on one column, the "averaged" column halves the other column's value
(e.g. 2-thiophenecarboxaldehyde 0.65 -> 0.33): it averages with zero.

### Table 4. GC-O and GC-MS/MS-SRM quantitation of the odour-active compounds (RI and FD on HP-5; µg/kg meat, mean ± SD, n = 3; "/" = not quantitated, no standard or S/N ≤ 10)

| compound | RI | odour | amount (µg/kg meat) | log2FD | ID |
|---|---|---|---|---|---|
| *Sulfur compounds* | | | | | |
| Methanthiol [sic] | 455 | sulfur, gasoline | / | 5 | RI/odor |
| **2-Methyl-3-furanthiol** | 875 | cooked meat | **3.393 ± 0.181** | **18** | RI/S/odor |
| **3-(Methylthio)propanal** | 912 | cooked potato | **6.649 ± 0.332** | **18** | RI/S/odor/MS |
| **2-Furfurylthiol** | 916 | sesame | **/** | **17** | RI/S/odor |
| **Dimethyl trisulfide** | 984 | sulfur, meaty | **0.888 ± 0.034** | 6 | RI/S/odor |
| 2-Mercaptothiophene | 986 | sesame | / | 6 | RI/S/odor |
| 2-Thiophenecarboxaldehyde | 1005 | sulfur | 0.453 ± 0.024 | 2 | RI/S/odor/MS |
| 2-Acetylthiazole | 1023 | sesame, roasted | 2.172 ± 0.017 | 9 | RI/S/odor/MS |
| 4-Methyl-5-thiazoleethanol | 1299 | sesame | 2.398 ± 0.092 | 4 | RI/S/odor/MS |
| *Heterocyclic compounds* | | | | | |
| 2-Acetyl-1-pyrroline | 950 | nutty, roasted | / | 8 | RI/S/odor |
| 2-Ethyl-3-methylpyrazine | 1023 | roasted | 0.195 ± 0.004 | 9 | RI/S/odor/MS |
| *Aldehydes* | | | | | |
| Pentanal | 687 | green | 14.295 ± 0.025 | 7 | RI/S/odor/MS |
| **Hexanal** | 794 | green | **26.654 ± 0.092** | 6 | RI/S/odor/MS |
| (E)-2-hexenal | 835 | green, fatty | 0.663 ± 0.042 | 8 | RI/S/odor/MS |
| Heptanal | 886 | sweet melon | 4.758 ± 0.102 | 5 | RI/S/odor/MS |
| (E)-2-heptenal | 960 | fatty | 2.667 ± 0.023 | 1 | RI/S/odor/MS |
| Phenylacetaldehyde | 1049 | rose | 0.884 ± 0.008 | 4 | RI/S/odor/MS |
| Nonanal | 1104 | sweet melon | 3.198 ± 0.064 | 15 | RI/S/odor/MS |
| (E)-2-nonenal | 1150 | cucumber | 0.317 ± 0.011 | 15 | RI/S/odor/MS |
| (E,E)-2,4-nonadienal | 1218 | fatty | / | 6 | RI/S/odor/MS |
| (E)-2-decenal | 1261 | fatty | 0.680 ± 0.018 | 5 | RI/S/odor/MS |
| (E,E)-2,4-decadienal | 1319 | fatty | 0.345 ± 0.014 | 15 | RI/S/odor/MS |
| *Ketones* | | | | | |
| 3-Hydroxy-2-butanone | 705 | yogurt | 560.500 ± 12.150 | 11 | RI/S/odor |
| 4-Hydroxy-2,5-dimethyl-3(2H)-furanone | 1055 | caramel | / | 13 | RI/S/odor |
| *Alcohols* | | | | | |
| 1-Octen-3-ol | 982 | mushroom | 11.737 ± 0.352 | 13 | RI/S/odor/MS |
| *Esters* | | | | | |
| γ-decalactone | 1430 | peach | / | 18 | RI/S/odor |
| *Others* | | | | | |
| Unknown | 1174 | fatty | / | 6 | |
| 2,4-Bis(1,1-dimethylethyl)-phenol | 1407 | woody | / | 15 | RI/odor/MS |

### [Z] Molar view, per kg wet meat stewed (MW: MFT 114.17; methional 104.17; DMTS 126.26; DMDS 94.20; hexanal 100.16)

| compound | µg/kg meat | nmol/kg meat | class |
|---|---|---|---|
| MFT | 3.39 | 29.7 | SRM + IS calibration, no RF |
| methional | 6.65 | 63.8 | same |
| DMTS | 0.888 | 7.0 | same |
| hexanal | 26.65 | 266 | same |
| DMDS | 0.36 | 3.8 | full scan, response factor 1 (semi-quant) |
| FFT | not quantified (< S/N 10 in SRM) | – | – |
| MFT : methional (molar) | 0.47 | | |
| methional : DMTS (molar) | 9.1 | | |

[Z] Per kg broth: the broth mass is not given. With 200 g meat in 200 mL water under reflux (no
evaporation) the filtered broth is plausibly 200-300 g, so µg/kg broth ≈ (0.7-1.0) × the per-meat
numbers. Order of magnitude unchanged.

[Z] Meat precursors on a wet basis (assuming ~28 % dry matter, consistent with 22.85 % protein and 1.81 %
fat; my assumption): total cysteine ~1.4 mg/g wet, total methionine ~5.0 mg/g wet. Both are acid
hydrolysates (protein-bound + free), not free amino acids.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor (external real-food levels, IS-calibrated without recovery correction):**
- **Pork broth: MFT 3.4, methional 6.6, DMTS 0.9, hexanal 27 µg/kg meat.** Joins Guth 1994 (beef juice) and
  Fan 2018 (chicken broth) as a third boiled-meat matrix, and is the only one of the three with a
  calibrated DMTS.
- **MFT ≈ 0.5 × methional (molar)** in a water-boiled lean pork: higher than Fan 2018's chicken (MFT :
  methional ≈ 0.006-0.011 molar) by 40-80×; the difference is 3× in MFT and 20× in methional. Treat as
  between-lab/matrix spread, not a pork-specific ratio.
- **DMTS present at ~1 µg/kg in a water-stewed meat with no added sulfur** (methional : DMTS ≈ 9 molar).
  A real-food counterpoint for the methionine lane, which refuses DMTS (`METHIONINE_NO_DMTS_REASON`): meat
  has cysteine-derived H2S, so this does not contradict the refusal reason for sugar + methionine pots.

**Cannot anchor:**
- **No kinetics** (one time point), **no pH, no sugar/ribose/thiamine/IMP data**, broth mass not given.
- **FFT: no number** despite log2FD 17; S/N ≤ 10 in SRM. Do not read as zero.
- **Amino acids are acid hydrolysates on a dry basis**, not free precursors.

**Caveats / flags:**
1. [!] **No recovery correction.** The same lab later (Fan 2018, Du 2020) measured SAFE recovery factors of
   0.77-0.94 for these compounds; for MFT the in-process loss (oxidation to the disulfide during 3 h
   reflux and work-up) is not covered by any RF. Treat MFT/DMTS as ±2× (judgement, ~65 %), more likely low
   than high.
2. [!] **IS concentration inconsistent:** 200 µg/mL in sec. 2.7 vs 0.5 µg/mL in Table 1. The SRM
   calibration (Table 1) is the one used for Table 4.
3. [!] **Full-scan vs SRM disagree** for the same compounds: methional 1.14 (Table 3, RF 1) vs 6.65 µg/kg
   (Table 4, calibrated); hexanal 32.97 vs 26.65; 2-acetylthiazole 2.98 vs 2.17. Use Table 4.
4. [!] **log2FD values of 17-18** imply dilutions of ~1:130 000-260 000, far beyond the "1:1024 and so on"
   series described; possibly log2FD of a stock rather than of the 0.35-mL concentrate. Rankings only.
5. [!] MFT was not seen in full-scan GC-MS (absent from Table 3) but quantified by SRM; the SRM transitions
   (85>45, 114>85) are specific enough for that to be credible.
6. [!] 2-Ethyl-3-methylpyrazine and 2-acetylthiazole share HP-5 RI 1023 in Table 4.
7. [!] "Stored in a refrigerator at 20 °C" (minus sign lost); fatty-acid unit "mg/g dry meat" conflicts with
   the text's "mg/g meat" (sec. 3, Table 2 note).

# Du, Ma, Li, Bai, Huang, Cui, Accoroni, Fan & Wang 2025 — EXTRACTION
### *"Effects of unsaturated C18 fatty acids on 'glucose-glutathione' Maillard reaction: Comparison and formation pathways of initial stage and meaty flavor compounds"*
### Glutathione + glucose 200 + 200 mmol/L in 0.2 M phosphate pH 6.50, ± 10 mmol/L-equivalent oleic, linoleic or linolenic acid, 145 °C: glutathione and glutathione-Amadori at 10-180 min, 65 volatiles at 120 min (MFT included, no FFT), [13C6]/[12C6]-glucose CAMOLA on five sulfur volatiles. Read 2026-10-27.

**Source on disk:** `data/articles/du2025.pdf` (publisher's PDF, 10 pp.). Text from `pdftotext -layout`;
Table 1 checked cell by cell against a 200 dpi page image (p. 5); Fig. 2B-C (p. 4) read from a 250 dpi
crop. Supplementary material (Fig. S1 lipid scheme, Fig. S2 HSQC) **not on disk**, not read.

**Provenance:** **[M]** printed number · **[Q]** qualitative statement · **[F]** read off a figure by me
(read error stated) · **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Wenbin Du, Qianli Ma, Yang Li, Shuang Bai, Yatao Huang, Weiye Cui, Cecilia Accoroni, Bei Fan*, Fengzhong Wang* |
| Title | "Effects of unsaturated C18 fatty acids on "glucose-glutathione" Maillard reaction: Comparison and formation pathways of initial stage and meaty flavor compounds" |
| Venue | *Food Research International* **201** (2025) **115645** |
| DOI | 10.1016/j.foodres.2024.115645 |
| Affiliation | Institute of Food Science and Technology, Chinese Academy of Agricultural Sciences, Beijing (+ Lanzhou, BTBU, CangZhou CAAS, INTA Argentina) |
| Dates | received 15 Oct 2024; revised 30 Nov 2024; accepted 28 Dec 2024; online 30 Dec 2024 |
| Data availability | "No data was used for the research described in the article." [!] (boilerplate) |

Same first author as Du et al. 2023 (FRI 173; 2-alkenals on cysteine-glucose), see `du2023_extraction.md`.

---

## 1. THE ONE-PARAGRAPH ANSWER

1 mmol glutathione (GSH) + 1 mmol glucose in 5 mL 0.2 M phosphate pH 6.50 (**200 + 200 mmol/L** [Z]),
alone (GG) or with 0.05 mmol oleic (GGO), linoleic (GGL) or linolenic acid (GGN) (**1 : 20 fatty acid :
GSH**, 10 mmol/L-equivalent, two-phase), heated at **145 °C**. Free GSH in all four systems falls
~4-fold per 10 min between 10 and 30 min (**apparent first-order 0.12-0.22 min⁻¹** [Z]) and is
essentially gone by 60 min; the GSH-Amadori compound is already at its highest at the first sample
(10 min) and decays as fast (Fig. 2). Fatty acids leave **more** GSH at 10 min (oleic 51, linoleic 49,
linolenic 44 vs 36 mg/L). At 120 min, **MFT is 135 µg/L in GG and 161-175 µg/L with any fatty acid
(+19 to +30 %)**, while summed thiophenes fall by half (2716 → 1028-1400 µg/L) and summed "sulfides"
other than MFT fall by 55-74 % [Z]. CAMOLA: in GG, 2-methylthiophene is 100 % glucose-carbon; with fatty
acids only 38-48 % is, i.e. **52-62 % of 2-methylthiophene carbon comes from the lipid** [Z]. The
absolute GSH scale is unusable (36 mg/L at 10 min is 0.06 % of the 61 g/L charged; section 3.1).

---

## 2. METHODS [M] (pp. 2-4)

| item | as printed |
|---|---|
| Reactants | oleic, linoleic, linolenic acid (99 %), L-glutathione (99 %), [12C6]-D-glucose (99 %), [13C6]-D-glucose (99 %), J&K |
| Systems | **GG**: GSH 1 mmol + glucose 1 mmol; **GGO / GGL / GGN**: + 0.05 mmol oleic / linoleic / linolenic acid. "dissolved in 5 ml phosphate buffer" (**0.2 M, pH 6.50**) [Z: 200 mmol/L GSH and glucose; fatty acid 10 mmol/L-equivalent, not soluble] |
| Heating | **145 °C, 120 min** for volatiles; **145 °C for 10, 20, 40, 60, 80, 120, 180 min** for GSH and GSH-Amadori (Fig. 2 x-axis instead shows 10, 20, 30, 40, 50, 60, 120, 180) [!] |
| Vessel, headspace, heating device, come-up, quench | **not stated** [!] (oxygen supply to the lipid is unknown) |
| pH after mixing / after heating | **not reported** [!]. 200 mmol/L GSH (an acid with three ionisable carboxyl/ammonium groups) in 0.2 M buffer will very likely pull the start pH well below 6.50 [Z, inference, not measured] |
| GSH-Amadori standard | GSH 0.2 mol + glucose 0.2 mol in phosphate buffer, 65 °C 1 h then 85 °C 6 h; AG 50W-X4 (H⁺), eluted water / 0.15 M NH₃; freeze-dried, ≥ 95 % (HPLC-UV); [M+H]⁺ 470.1739 |
| GSH, GSH-Amadori assay | UPLC-MS/MS (Xevo TQ-XS), MRM 308/179 (GSH), 470/386 (Amadori), BEH C18; external standard. **Dilution of the reaction mixture before injection not stated** [!] |
| Volatiles | SPME CAR/PDMS/DVB, sample at 50 °C 20 min then fibre 30 min; IS **1 µL of 1,2-dichlorobenzene 100 µg/mL** (methanol); GC-MS DB-Wax, 7890B/5977A; Sniffer 9000 olfactometer. Response factor 1 implied, not stated |
| CAMOLA | [13C6]-glucose 0.5 mmol + [12C6]-glucose 0.5 mmol in place of 1 mmol glucose; isotopomer proportions from M⁺ to M⁺+n, corrected for 1.10 % 13C, 0.76 % 33S, 4.20 % 34S |
| Replicates | "All the data were repeated three times" |

---

## 3. RESULTS

### 3.1 Fig. 2B-C (p. 4): glutathione and GSH-Amadori, mg/L [F]

Bars, mean ± error bar; read error **±1 mg/L (GSH > 10), ±0.2 mg/L (GSH < 5), ±0.2 mg/L (Amadori)**.
Letters are within a time point, as printed.

**Glutathione (mg/L):**

| t (min) | GG | GGO | GGL | GGN |
|---|---|---|---|---|
| 10 | 36.0 ± 2 (c) | 50.9 ± 1.2 (a) | 48.8 ± 3 (a) | 43.5 ± 1 (b) |
| 20 | 11.3 ± 1 (ab) | 12.3 ± 0.5 (a) | 10.4 ± 0.5 (b) | 11.1 ± 0.4 (b) |
| 30 | 2.9 (b) | 3.8 (a) | 1.2 (c) | 1.8 (c) |
| 40 | 0.4 (c) | 1.1 (a) | 0.6 (b) | 1.0 (a) |
| 50 | 0.3 (b) | 0.6 (a) | 0.2 (c) | 0.7 (a) |
| 60 | ~0.1 | ~0.3 | ~0.1 | ~0 |
| 120 | ~0.1 | ~0.1 | ~0.1 | ~0.1 |
| 180 | ~0 | ~0 | ~0 | ~0 |

**GSH-Amadori (N-(1-deoxy-fructos-1-yl)-glutathione, mg/L):**

| t (min) | GG | GGO | GGL | GGN |
|---|---|---|---|---|
| 10 | 8.35 (c) | 15.1 ± 1.9 (a) | 11.2 ± 1.2 (b) | 9.8 ± 1.1 (b) |
| 20 | 2.4 ± 0.8 (a) | 2.1 ± 0.7 (a) | 1.7 ± 0.4 (a) | 1.1 ± 0.3 (ab) |
| 30 | 0.45 (a) | 0.55 (a) | 0.65 (a) | 0.5 (a) |
| 40 | ~0.05 (c) | 0.1 (b) | 0.45 (a) | 0.1 (b) |
| 50 | ~0.02 (c) | ~0.02 (c) | 0.2 (a) | 0.15 (b) |
| ≥ 60 | ~0 | ~0 | ~0 | ~0 |

**Apparent first-order GSH loss, 145 °C [Z]** (k = ln(C₁/C₂)/10 min):

| interval | GG | GGO | GGL | GGN |
|---|---|---|---|---|
| 10 → 20 min | 0.116 | 0.142 | 0.155 | 0.137 min⁻¹ |
| 20 → 30 min | 0.136 | 0.117 | 0.216 | 0.182 min⁻¹ |

**[!] The absolute scale cannot be right as a whole-pot concentration.** 200 mmol/L GSH = 61.5 g/L
charged; 36 mg/L at 10 min would be **0.06 %** remaining, which requires k ≥ 0.74 min⁻¹ over the first
10 min, five times the 10-30 min slope [Z]. Either the samples were diluted before UPLC-MS/MS by an
unstated factor (~10³), or the number is not a pot concentration. **Use only ratios between systems and
the decay shape.** GSH-Amadori is likewise 0.018 mmol/L at 10 min against 200 mmol/L charged.

**[!] Text vs figure.** (i) Text: GSH-Amadori "first increase and then decrease"; the figure has no
point before 10 min and only decreases. (ii) Text: GSH is consumed "the fastest" in GG and fatty acids
inhibit consumption "oleic > linoleic > linolenic". True at 10 min only; at 30 min GGL (1.2) and GGN
(1.8) are below GG (2.9). (iii) Methods list sampling at 80 min; the figure has 30 and 50 min instead.

### 3.2 Table 1 (pp. 5-6): volatiles at 145 °C, 120 min, µg/L, mean ± SD, letters as printed [M]

Blank = not detected. OD = odour at the sniffing port; ID: RI, MS (NIST 15), O, S (authentic standard).

| No. | RI | Compound | GG | GGO | GGL | GGN | OD | ID |
|---|---|---|---|---|---|---|---|---|
| | | **thiophenes** | | | | | | |
| 1 | 1021 | thiophene | 45.75 ± 1.74 a | 27.32 ± 3.37 d | 35.92 ± 2.32 bc | 38.73 ± 3.27 bc | meaty | RI/MS/O/S |
| 2 | 1078 | 2-methylthiophene | 213.32 ± 21.56 a | 166.69 ± 4.89 bc | 168.04 ± 6.92 bc | 113.86 ± 1.82 d | meaty | RI/MS/O/S |
| 3 | 1192 | 2,5-dimethylthiophene | 167.04 ± 4.12 ab | 159.71 ± 7.63 ab | 67.55 ± 12.00 cd | 72.96 ± 8.03 cd | meaty | RI/MS/O/S |
| 4 | 1210 | 2,3-dimethylthiophene | 105.81 ± 10.16 | | | | meaty | RI/MS/O/S |
| 5 | 1259 | 2-methoxy-5-methylthiophene | 98.14 ± 2.49b | 38.54 ± 0.47c | 146.95 ± 8.27 a | | meaty | RI/MS/O/S |
| 6 | 1491 | 2-hexylthiophene | | "21 ± 1.93 ± 1.80c" [!] | 105.19 ± 3.47 ab | 106.04 ± 7.30 ab | meaty | RI/MS/O/S |
| 7 | 1496 | 3-(methylthio)-thiophene | 185.33 ± 8.11 a | | 105.86 ± 1.56b | 42.76 ± 12.18c | meaty | RI/MS/O |
| 8 | 1502 | 2-methyl-5-(methylthio)-thiophene | 130.43 ± 5.89 | | | | – | RI/MS |
| 9 | 1652 | 2-heptylthiophene | | 272.19 ± 22.65 ab | 153.97 ± 12.79c | 289.05 ± 22.20 ab | meaty | RI/MS/O/S |
| 10 | 1677 | 2-thiophenecarboxaldehyde | 89.77 ± 2.67 | 95.70 ± 1.70 | | | meaty | RI/MS/O/S |
| 11 | 1721 | 5-methyl-2(5H)-thiophenone | 112.81 ± 1.20 | | | 11.34 ± 1.87 | – | RI/MS/S |
| 12 | 1785 | 5-methyl-2-thiophenecarboxaldehyde | 270.25 ± 12.42 ab | 246.75 ± 27.18 ab | 134.91 ± 13.80c | 104.61 ± 0.27 d | – | RI/MS/S |
| 13 | 1789 | 1-(3-thienyl)-ethanone | 19.82 ± 1.82b | 27.74 ± 4.80 a | | 16.47 ± 0.78c | – | RI/MS/S |
| 14 | 1796 | 2-acetyl-5-methylthiophene | 27.75 ± 7.36b | 14.11 ± 0.46c | 37.80 ± 3.11 a | | meaty | RI/MS/O/S |
| 15 | 1817 | 1-(2-thienyl)-1-propanone | 167.08 ± 13.36 a | 29.46 ± 5.21cd | 53.60 ± 0.76b | 26.15 ± 1.48 cd | – | RI/MS/S |
| 16 | 1822 | benzo[b]thiophene-2-ol | 130.93 ± 8.61 ab | 127.84 ± 6.20 ab | 71.03 ± 6.17c | 56.85 ± 6.30 d | – | RI/MS |
| 17 | 1840 | thieno[3,2-b]thiophene | 325.56 ± 29.16 a | 80.70 ± 4.79c | 172.87 ± 23.90b | 58.14 ± 1.57d | meaty | RI/MS/O/S |
| 18 | 1918 | 2-methylthieno[2,3-b]thiophene | 381.62 ± 33.87 a | 95.70 ± 9.79b | 79.22 ± 1.80c | 58.10 ± 1.57 d | meaty | RI/MS/O/S |
| 19 | 1930 | 2-thiophenecarboxylic acid | 8.73 ± 0.31 | | | | – | RI/MS/S |
| 20 | 1951 | 3-(vinylthio)thiophene | 67.94 ± 9.24 a | 17.21 ± 0.49 d | 26.03 ± 1.56b | 20.74 ± 0.47c | – | RI/MS |
| 21 | "19,668" [!] | 2-ethylthieno[2,3-b]thiophene | 167.52 ± 17.00 | | | 12.28 ± 0.08 | – | RI/MS/S |
| | | sum | 2715.60 | 1399.66 | 1358.94 | 1028.08 | | |
| | | **furans** | | | | | | |
| 22 | 1230 | 2-pentylfuran | | 41.77 ± 4.22c | 108.87 ± 8.58b | 225.30 ± 11.49 a | sweet | RI/MS/O/S |
| 23 | 1416 | trans-2-(2-pentenyl)furan | | 21.33 ± 4.85c | 158.71 ± 2.46b | 203.25 ± 30.36 a | – | RI/MS |
| 24 | 1425 | 2-heptylfuran | | 119.55 ± 14.35 ab | 99.23 ± 0.62 ab | 54.08 ± 2.38c | sweet | RI/MS/O/S |
| 25 | 1713 | 5-methyl-2(3H)-furanone | 13.95 ± 2.13 d | 15.71 ± 0.16 bc | 18.63 ± 0.31 a | 14.40 ± 1.15 bc | sweet | RI/MS/O/S |
| 26 | 1425 | 3-furaldehyde | 39.94 ± 6.84 | 35.44 ± 0.09 | | | –- | RI/MS |
| 27 | 1502 | 1-(2-furanyl)-ethanone | 469.00 ± 29.45 ab | 449.95 ± 35.09cd | 413.46 ± 27.22 cd | 486.77 ± 12.08 ab | – | RI/MS/S |
| 28 | 1538 | 2-octylfuran | | 201.30 ± 14.26 a | 99.10 ± 5.87b | 62.80 ± 11.36c | – | RI/MS/S |
| 29 | 1567 | 5-methyl-2-furancarboxaldehyde | 670.24 ± 97.38 a | 100.31 ± 4.67 d | 125.77 ± 16.60c | 518.30 ± 30.53b | – | RI/MS/S |
| 30 | 2020 | dihydro-5-pentyl-2(3H)-furanone | | 147.18 ± 2.88 a | 50.65 ± 5.66b | 23.00 ± 0.46c | sweet | RI/MS/O/S |
| 31 | 2450 | 5-hydroxymethylfurfural | 35.63 ± 4.58c | 117.62 ± 10.29 ab | 126.73 ± 10.40 ab | | – | RI/MS/S |
| | | sum | 1228.76 | 1250.16 | 1201.15 | 1587.9 | | |
| | | **others (sulfides)** | | | | | | |
| 32 | 1314 | **2-methyl-3-furanthiol** | **135.01 ± 14.02 d** | **160.79 ± 9.34 abc** | **175.26 ± 7.78 abc** | **170.37 ± 15.89 abc** | meaty | RI/MS/O/S |
| 33 | 1428 | 2-methyl-5-(methylthio)-furan | 116.57 ± 4.67 | | | | – | RI/MS/S |
| 34 | 1583 | anti-3,5-dimethyl-1,2,4-trithiolane | 70.69 ± 9.21b | 127.43 ± 12.96 a | 46.45 ± 8.06 cd | 37.81 ± 6.73 cd | – | RI/MS/S |
| 35 | 1842 | syn-3,5-dimethyl-1,2,4-trithiolane | 108.51 ± 5.29 a | 10.51 ± 3.04 d | 58.07 ± 3.56 bc | 52.56 ± 5.57 bc | – | RI/MS/S |
| 36 | 1728 | 1,2,3-trithiolane | 36.49 ± 13.39cd | 51.31 ± 9.21 ab | 61.63 ± 10.21 ab | 25.27 ± 3.18 cd | meaty | RI/MS/O/S |
| 37 | 1751 | 1,2,5-trithiolane | 481.18 ± 11.31 a | 57.49 ± 7.80 bc | 58.07 ± 3.56 bc | 32.96 ± 4.07 d | – | RI/MS/S |
| 38 | 1795 | 4-mercaptophenol | 19.15 ± 0.47 a | 7.33 ± 0.16 d | 13.16 ± 0.31c | 15.25 ± 0.31b | meaty | RI/MS/O |
| 39 | 1821 | 2-acetylthiazole | 86.11 ± 13.51 ab | 35.53 ± 1.07c | 90.70 ± 5.57 ab | | brunt | RI/MS/O/S |
| 40 | 1852 | 1(3H)-isobenzofuranthione | 76.57 ± 1.56 | | | 58.67 ± 3.11 | – | RI/MS |
| 41 | 1910 | 1,3-benzenedithiol | 135.19 ± 9.60 a | 27.97 ± 1.56cd | 51.76 ± 9.34b | 28.69 ± 4.16 cd | – | RI/MS |
| 42 | 1949 | benzothiazole | 115.67 ± 6.22b | 182.70 ± 7.64 a | | 81.96 ± 4.67c | – | RI/MS/S |
| 43 | 1917 | 1,2,4-trithiolane | 54.30 ± 13.82 a | | 16.93 ± 4.78 bc | 10.07 ± 2.32 bc | – | RI/MS |
| 44 | 1929 | 1,3,5-trithiane | 35.83 ± 1.69 | 100.68 ± 4.20 | | | – | RI/MS |
| | | sum | 1471.27 | 600.95 | 572.03 | 513.61 | | |
| | | **lipid degradation produces** | | | | | | |
| 45 | 1096 | hexanal | | 5.64 ± 0.12c | 63.83 ± 1.56 a | 6.70 ± 0.31b | green | RI/MS/O/S |
| 46 | 1299 | octanal | | 135.84 ± 20.63 a | 23.15 ± 5.76 bc | 25.00 ± 1.27 bc | fat | RI/MS/O/S |
| 47 | 1390 | nonanal | | 216.69 ± 42.32 a | 54.10 ± 3.11b | 26.81 ± 1.23c | fat | RI/MS/O/S |
| 48 | 1420 | 2-propyl-2-heptenal | | 15.48 ± 1.91b | 109.18 ± 3.11 a | 9.34 ± 0.62c | – | RI/MS |
| 49 | 1544 | (E)-2-nonenal | | | 72.79 ± 11.39 | 38.42 ± 1.56 | fat | RI/MS/O/S |
| 50 | 1804 | (E,E)-2,4-decadienal | | 233.72 ± 20.67 a | 114.65 ± 47.28b | 17.32 ± 2.16c | fat | RI/MS/O/S |
| 51 | 1810 | (E,Z)-2,4-decadienal | | | 30.00 ± 0.93 | 10.01 ± 0.16 | fat | RI/MS/O/S |
| 52 | 1462 | 1-heptanol | | 18.74 ± 1.82c | 147.62 ± 7.31 a | 53.44 ± 4.35b | sweet | RI/MS/O/S |
| 53 | 1288 | 1-hydroxy2-propanone | | 21.01 ± 0.31 | | 68.36 ± 6.62 | – | RI/MS |
| 54 | 1315 | 1-octen-3-one | | 22.25 ± 1.56 | | 269.44 ± 11.92 | mushroom | RI/MS/O/S |
| 55 | 1494 | 2-decanone | | 149.72 ± 4.09 a | 117.79 ± 10.47b | 19.46 ± 7.37c | – | RI/MS/S |
| 56 | 2035 | 2-pentadecanone | | 8.13 ± 0.94 | | 15.92 ± 0.26 | – | RI/MS/S |
| 57 | 1462 | acetic acid | | 73.70 ± 5.51 | | 80.86 ± 3.51 | acid | RI/MS/O/S |
| 58 | 1635 | 2-propenoic acid | | | 125.65 ± 13.93 | 23.51 ± 3.85 | – | RI/MS/S |
| 59 | 1555 | hexanoic acid | | 337.82 ± 22.34 ab | 309.46 ± 9.13 ab | 133.97 ± 2.41c | – | RI/MS/O/S |
| 60 | 1962 | heptanoic acid | | 715.72 ± 39.40 | 640.12 ± 59.82 | 503.22 ± 34.45 | – | RI/MS/S |
| 61 | 2066 | octanoic acid | | | 320.92 ± 22.37 | 993.63 ± 108.72 | acid | RI/MS/O/S |
| 62 | 2174 | nonanoic acid | | 391.66 ± 32.11 a | 269.63 ± 8.58b | 74.93 ± 10.28c | – | RI/MS/S |
| 63 | 2280 | decanoic acid | | 42.44 ± 3.01 bc | 181.09 ± 16.19 a | 39.00 ± 3.11 bc | – | RI/MS/S |
| 64 | 2433 | benzoic acid | | 26.56 ± 0.16 | 26.39 ± 1.24 | | – | RI/MS/S |
| 65 | 2868 | hexadecanoic acid | | 288.56 ± 33.68b | 155.68 ± 11.76c | 514.13 ± 3.32 a | – | RI/MS/S |
| | | sum | | 2703.68 | 2762.05 | 2923.47 | | |
| | | **total** | 10831.26 | 9005.22 | 9026.29 | 9182.65 | | |

Footnotes as printed: RI against C5-C30 n-alkanes; system codes; OD "–" not detected; ID RI, MS (NIST 15),
O odour, S authentic standard. Compound counts in the text (37 / 52 / 49 / 53) match the table for GG and
GGO [Z].

**Sum checks [Z].** GG thiophenes, furans, sulfides; all four furan sums; GGL and GGN thiophene and
sulfide sums; all three lipid sums re-add exactly. Exceptions [!]:
- **GGO thiophene sum 1399.66 omits 2-hexylthiophene** (row 6, garbled cell "21 ± 1.93 ± 1.80c");
  with 21 it would be 1420.66.
- **GGO sulfide sum 600.95 omits MFT** (761.74 − 160.79 = 600.95).
- **"total" is not the sum of the class sums.** GG 10831.26 = exactly 2 × 5415.63 (the three class sums);
  GGO/GGL/GGN totals exceed their class sums by 3051 / 3132 / 3130 µg/L. The totals are unusable.

**Molar MFT [Z]** (MW 114.17): GG 1.18 µmol/L = **5.9 × 10⁻⁶ mol/mol GSH**; GGO 1.41, GGL 1.54, GGN 1.49
µmol/L (7.0-7.7 × 10⁻⁶ mol/mol). No FFT, no MFT disulfide reported in any system.

### 3.3 Table 2 (p. 7): CAMOLA isotopomer distributions, % of molecular-ion cluster [M]

[13C6]-glucose : [12C6]-glucose 1 : 1. Columns: number of 13C atoms in the isotopomer.

| No. | Compound | M⁺ (m/z) | System | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|---|---|
| a, b | 2-methylthiophene | 98 | GG | 48 | 1 | <1 | <1 | 1 | 50 | |
| | | | GGO | 75 | <1 | <1 | <1 | 1 | 24 | |
| | | | GGL | 80 | <1 | <1 | 1 | <1 | 19 | |
| | | | GGN | 81 | <1 | <1 | <1 | <1 | 19 | |
| c | 2-thiophenecarboxaldehyde | 112 | GG | 50 | <1 | <1 | <1 | <1 | 50 | |
| | | | GGO | 49 | <1 | 1 | <1 | <1 | 50 | |
| d | 4-mercaptophenol | 126 | GG | 50 | <1 | <1 | <1 | <1 | <1 | 50 |
| | | | GGO | 51 | <1 | <1 | <1 | <1 | 1 | 49 |
| | | | GGL | 49 | 1 | <1 | <1 | <1 | <1 | 50 |
| | | | GGN | 50 | <1 | <1 | <1 | <1 | 1 | 49 |
| e | 2-acetylthiazole | 127 | GG | 49 | <1 | <1 | 50 | <1 | 1 | <1 |
| | | | GGO | 50 | <1 | <1 | 49 | <1 | <1 | <1 |
| | | | GGL | 50 | 1 | <1 | 50 | <1 | <1 | <1 |
| f | thieno[3,2-b]thiophene | 140 | GG | 60 | <1 | 4 | <1 | 2 | 6 | 28 |
| | | | GGO | 59 | <1 | 5 | <1 | 2 | 5 | 29 |
| | | | GGL | 59 | <1 | 4 | <1 | 2 | 6 | 29 |
| | | | GGN | 61 | <1 | 5 | <1 | 2 | 7 | 25 |

**Derived carbon origin [Z]** (with a 1 : 1 label mix, the glucose-derived fraction of a compound that
takes an intact glucose fragment is ≈ 2 × the fully labelled isotopomer; assumes equal reactivity of 12C
and 13C glucose):

| compound | GG glucose-derived | GGO | GGL | GGN | other carbon source |
|---|---|---|---|---|---|
| 2-methylthiophene (C5 fragment) | ~100 % | **48 %** | **38 %** | **38 %** | the rest unlabelled: lipid (authors: 2-pentenal + H₂S) |
| 2-thiophenecarboxaldehyde (C5) | ~100 % | ~100 % | n.d. | n.d. | none |
| 4-mercaptophenol (C6) | ~100 % | ~98 % | ~100 % | ~98 % | none (lipid does not contribute) |
| 2-acetylthiazole (C3 from glucose) | C3 ~100 % | ~100 % | ~100 % | n.d. | C2 of the ring from GSH cysteine |
| thieno[3,2-b]thiophene (C6) | 56 % intact C6 | 58 % | 58 % | 50 % | ~30 % unlabelled already in GG, so not lipid; GSH-derived carbon or a non-glucose C2 + C4 route |

Applied to Table 1 [Z]: lipid-derived 2-methylthiophene ≈ 87 (GGO), 104 (GGL), 71 (GGN) µg/L; glucose-
derived 2-methylthiophene falls from 213 (GG) to 80, 64, 43 µg/L, a 62-80 % cut.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** `GSH` is "a sulfide source only" (`species_sulfur.py` line 367): one step,
`r_gsh_h2s` GSH → H₂S + fragments, constant `k_gsh_h2s` from Zheng & Ho 1994 Table III
(`parameters_sulfur.py` lines 369-380; 80-110 °C, pH 3-9). Its amine, glycation and peptide chemistry are
not modelled. The lipid lane makes oxidation aldehydes; there is no lipid → thiophene step. Nothing below
was run against the model.

**(1) Glutathione consumption in a sugar pot is far faster than its thermolysis (the key number).**
The model's own Table III constants, extrapolated by the code's refit Arrhenius pairs to 145 °C, give
roughly **4.5 × 10⁻³ min⁻¹ at pH 7 and 1.3 × 10⁻³ min⁻¹ at pH 5** [Z: from k(110 °C) 3.89 × 10⁻⁴ and
7.46 × 10⁻⁵ min⁻¹, Ea ≈ 93 and 109 kJ/mol from the 80/110 °C points; an extrapolation 35 K beyond the
data]. Free GSH here disappears at **0.12-0.22 min⁻¹** at 145 °C with equimolar glucose, i.e. **25-170 ×
faster**. Most GSH therefore leaves through glucose-dependent routes (Amadori/condensation), not through
direct H₂S release. **Consequence for the model:** in a GSH + reducing-sugar pot, a GSH species with only
`r_gsh_h2s` will keep GSH (and its sulfur) far too long, and will route all of it to H₂S. This is a
**qualitative-to-semi-quantitative contradiction**, not a fit target: the absolute scale of Fig. 2 is
broken (section 3.1), the pH is unknown, and the first sample is at 10 min. The ratio
k(loss)/k(thermolysis) ≥ ~25 is robust to the scale problem because it uses only the decay shape.

**(2) Lipid × sulfur interaction, quantified.**
- **MFT goes up, not down, with free C18:1/2/3 acids**: 135 → 161-175 µg/L (+19 to +30 %; letters d vs
  abc). Direction test for any lipid → sulfur coupling. The authors' "meaty flavour decreased" applies to
  summed thiophenes (−48 to −62 %) and the non-MFT sulfides (GG 1336 → 601 / 397 / 343 µg/L, −55 to −74 % [Z]), not to MFT.
- **Carbon origin**: with fatty acid present, **52-62 % of 2-methylthiophene carbon is non-glucose**
  (lipid). Only thiophene-type products incorporate lipid carbon; 4-mercaptophenol, 2-acetylthiazole and
  2-thiophenecarboxaldehyde do not. If the model ever adds an aldehyde + H₂S → alkylthiophene sink for the
  lipid lane, these give the size of the lipid-carbon share at 1 : 20 fatty acid : GSH.
- **H₂S competition**: summed sulfur volatiles fall when lipid is present (sections 3.2). Consistent
  with lipid-derived carbonyls competing for H₂S, but not separable from the higher residual GSH at 10 min.

**(3) Lipid lane yields (weak).** At 145 °C, 120 min, aqueous, 10 mmol/L-equivalent fatty acid,
hexanal from linoleic acid is **63.83 µg/L = 0.64 µmol/L ≈ 6 × 10⁻⁵ mol/mol** [Z]; nonanal from oleic
216.69 µg/L ≈ 1.5 × 10⁻⁴ mol/mol; 2,4-decadienal (E,E) from linoleic 114.65 µg/L ≈ 7.5 × 10⁻⁵ mol/mol.
**Do not use these as precursor-specific yields** [!]: GGO (oleic) carries the most (E,E)-2,4-decadienal
(233.72), a linoleate product; GGN (linolenic) carries the most 2-pentylfuran (225.30) and 1-octen-3-one
(269.44), both n-6 products; hexadecanoic acid (C16) appears in all lipid systems up to 514 µg/L. The
"99 %" acids are evidently not pure, or the assignments are wrong. The oxygen supply (vessel, headspace)
is not stated, so no oxidation rate can be inferred.

**(4) Not usable:** absolute GSH and GSH-Amadori concentrations (scale); the "total" row; RI 19,668 and
the garbled GGO 2-hexylthiophene cell; the heat map and PCA (Fig. 4); the proposed pathways (Fig. 3,
beyond what Table 2 shows).

**Caveats, collected.**
1. pH 6.50 is the buffer label; 200 mmol/L GSH in 0.2 M buffer almost certainly acidifies the pot; no
   measured pH.
2. Vessel, headspace, heating, come-up and quench unstated.
3. Fatty acids at 0.05 mmol in 5 mL water form a separate phase; effective concentration unknown.
4. GSH/Amadori absolute scale off by ~10³ against the charge; sampling times in Methods and Fig. 2
   disagree.
5. Volatiles semi-quantitative (internal standard, RF 1 implied). Two printed sums omit rows; totals are
   inconsistent.
6. Lipid product patterns contradict the stated precursor (decadienal from oleic, pentylfuran from
   linolenic).
7. No FFT (glucose pot) and no MFT disulfide reported.

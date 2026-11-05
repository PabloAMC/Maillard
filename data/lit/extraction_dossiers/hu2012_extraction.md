# Hu, Wang, Xie & Liu 2012 (Chin. J. Chromatogr. 30:298) — EXTRACTION (file stem `hu2012`)
### 胡斌, 王昇, 谢复炜, 刘惠民, "高效液相色谱-蒸发光散射检测法在烟用香精料液中水溶性糖及山梨醇分析中的应用" / *"Application of HPLC-ELSD in determination of water-soluble sugars and sorbitol in tobacco flavourings and casings"*
### A method paper: rhamnose, **xylose**, fructose, sorbitol, glucose, sucrose, maltose by **HPLC-ELSD** (Prevail Carbohydrate ES 250 × 4.6 mm, MeCN / water 75:25, 1 mL/min, 5 µL; drift tube 79 °C, N2 2.0 L/min) with **double-logarithmic calibrations**: **xylose ln A = 8.44 + 2.00 ln c (c in g/L), 0.064-1.278 g/L, r² 0.9984; LOD (S/N = 3) 19 mg/L = 127 µM** (all seven: 12-26 mg/L). No Maillard chemistry; its value is as an **independent measurement-model floor and response law for ELSD of sugars**.

**Read 2026-11-05 (first reader) from `data/articles/hu2012.pdf`** (6 pp., publisher PDF, clean Chinese / English text
layer; tables in English). All numbers from `pdftotext -raw`, Tables 2 and 5 checked for column alignment. Figures
(chromatograms; ln C1 and C2 vs mobile phase) not digitised — Table 5 carries their numbers.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **HU Bin*** (binhu74@hotmail.com), WANG Sheng, XIE Fuwei, LIU Huimin — Zhengzhou Tobacco Research Institute of China National Tobacco Corporation |
| Venue | *Chinese Journal of Chromatography* (色谱) **30**(3), 298-303, March 2012 |
| DOI | **10.3724/SP.J.1123.2011.11014** (printed) — matches the brief |
| Received | 8 Nov 2011 |

**Laboratory.** Analytical-chemistry group of a tobacco institute; no Maillard kinetics. A laboratory id is not needed
unless a fit row ever uses it (proposed `ztri_hu` if so).

## 1. What this paper answers for the model

| model need | what Hu 2012 gives | tag |
|---|---|---|
| ELSD LOD for xylose | **19 mg/L (127 µM) at S/N = 3**, 5 µL injection → **~95 ng on column** [Z] | [M] |
| ELSD response law | **power law, A ∝ m^C2 with C2 = 2.00 (xylose, 75:25)**; 1.94 at 80:20 and 70:30; 1.67-2.39 across the seven analytes; C2 falls with retention time [Q] | [M] |
| calibrated range | xylose 0.064-1.278 g/L (0.43-8.5 mM) at 75:25; 0.128-1.278 g/L at 80:20; 0.032-1.278 g/L at 70:30 | [M] |
| precision | 6 replicate analyses of a casing sample: RSD 1.1-5.1 % | [M] |
| recovery | spiked fructose / sorbitol / glucose / sucrose / maltose 88-109 % (xylose not spiked) | [M] |

## 2. Conditions [M]

Agilent 1200 + Alltech 2000 ELSD; Prevail Carbohydrate ES (polymer-based) 250 × 4.6 mm, 5 µm, guard 7.5 × 4.6 mm;
MeCN / water 75:25 (v/v), 1.0 mL/min, 5 µL injection; drift tube 79 °C, N2 2.0 L/min. Standards (Sigma, > 99 %) in a mixed
stock ~1.2 g/L, serially diluted. Samples: 0.2 g in 20 mL water, ultrasonic 40 min, centrifuged, 0.5 mL through Isolute
ENV+ SPE (50 mg), eluted with 1 mL water into 2 mL. LOD: "signal-to-noise ratio of 3" [V, translated] by serial dilution.

## 3. Data [M]

### Table 2 (75:25): double-log calibration, Y = peak area, X = g/L; LOD in mg/L

| analyte | ln Y = a + b ln X | range (g/L) | r² | LOD (mg/L) | LOD (µM) [Z] |
|---|---|---|---|---|---|
| rhamnose | 8.94 + 2.17 | 0.065-1.033 | 0.9988 | 26 | 158 (MW 164.2) |
| **xylose** | **8.44 + 2.00** | **0.064-1.278** | **0.9984** | **19** | **127** (MW 150.1) |
| fructose | 9.49 + 1.87 | 0.063-1.001 | 0.9997 | 12 | 67 |
| sorbitol | 9.43 + 1.91 | 0.062-0.987 | 0.9997 | 20 | 110 |
| glucose | 8.97 + 1.83 | 0.062-1.231 | 0.9999 | 18 | 100 |
| sucrose | 8.84 + 1.80 | 0.061-1.222 | 0.9983 | 22 | 64 |
| maltose | 8.25 + 1.76 | 0.062-1.239 | 0.9988 | 23 | 67 (anhydrous, 342.3) |

### Table 5 (xylose rows): retention time / ln C1 / C2 / r² / range

80:20: 9.74 min / 7.41 / 1.94 / 0.9986 / 0.128-1.278 g/L · 75:25: 7.26 / 8.44 / 2.00 / 0.9984 / 0.064-1.278 ·
70:30: 7.12 / 8.71 / 1.94 / 0.9961 / 0.032-1.278.

### Arithmetic [Z]

1. With C2 ≈ 2, a straight-line (A vs c) calibration through this range is badly biased: a least-squares line through
   A ∝ c² on seven standards spanning 0.064-1.278 g/L has a large negative intercept (−0.22 to −0.29 of the top
   response) and, inverted, returns ~0.17-0.22 g/L for the 0.064 g/L standard (×3 high) and +11 % at the top — the
   signature of linear ELSD calibrations such as `chen2026`'s (Thr intercept −5.23 × 10⁵).
2. The LOD in moles is ~0.1 mM for monosaccharides: **a floor for any HPLC-ELSD sugar / sugar-derivative series of this
   method class**, about three orders above `chen2026`'s printed ELSD LODs (0.2-0.4 µM) and the same order as
   `deng2022`'s methionine ELSD LOD (89.7 µM, log-log calibration, exponent 1.47).

## 4. What it can and cannot test

**Can:** set an independent, printed ELSD detection floor (~0.1 mM, ~100 ng on column) and response exponent (C2 ≈ 2 for
xylose) for the measurement model of any ELSD-quantified sugar, ARP or amino acid in the corpus, and flag records whose
printed ELSD LODs are implausibly low or whose calibrations are linear.

**Cannot:** say anything about Maillard kinetics; transfer exactly to other instruments (C1 and C2 depend on mobile
phase, retention and detector design, as the paper shows); give xylose recovery (xylose was not spiked).

## 5. Conflicts and agreements with the repo's records

1. **[!] CONTRADICTS `chen2026`'s ELSD LODs** (Thr-ARP 0.199 µM, Thr 0.421 µM): ~300-600× below this paper's sugar LODs (~200-450× below `deng2022`'s methionine LOD),
   on comparable injection volumes (8 vs 5 µL).
2. **AGREES with `deng2022`'s ELSD floor** (methionine 89.7 µM, log-log calibration).
3. Other ELSD records on disk (`cao2017`, `cui2022`, `xia2022`, `zhang2025b`, ... — 20 dossiers mention ELSD) were not
   re-checked against this floor today; a sweep is the obvious follow-up.

## 6. Proposed ledger entry

None: a method paper with no reaction. If the measurement-model work wants a record, an `other` entry in
`aroma_other.yml` carrying the Table 2 xylose row and LOD (first reader only, second reader pending) would do.

## 7. Supplementary material

None.

## Second reader (2026-11-05)

**Independent second reading** by `agent-second-reader-2026-11-05`, PDF read before this dossier. Source:
`data/articles/hu2012.pdf` (Chin. J. Chromatogr. 30:298), Table 2 and Table 5 (bilingual captions). Printed numbers.

| item | location | first reader | second reader | verdict |
|---|---|---|---|---|
| xylose LOD | Table 2, sec. 2.2.1 | 19 mg/L, S/N = 3 | **19 mg/L, S/N = 3** (by serial dilution of the standards); 19 / 150.13 = 127 µM | AGREE |
| calibration form | Table 2 footnote | ln A = 8.44 + 2.00 ln c, c in g/L | **lnY = 8.44 + 2.00 lnX, Y = peak area, X = mass concentration in g/L** (double-log) | AGREE |
| linear range, r² | Table 2 | 0.064-1.278 g/L, 0.9984 | **0.064-1.278 g/L, r² 0.9984** | AGREE |
| conditions | sec. 1.2, abstract | Prevail Carbohydrate ES 250 × 4.6 mm 5 µm, MeCN / water 75:25, 1.0 mL/min, 5 µL, drift tube 79 °C, N2 2.0 L/min | same; Agilent 1200 + Alltech 2000 ELSD | AGREE |
| other mobile phases | Table 5 | 80:20: 7.41 + 1.94, 0.128-1.278; 70:30: 8.71 + 1.94, 0.032-1.278 | same | AGREE |

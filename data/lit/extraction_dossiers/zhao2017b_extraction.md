# Zhao et al. 2017 (Food Science, in Chinese) — EXTRACTION (L-cysteine + D-xylose 2 + 2 mmol in 5 mL 0.2 M sodium phosphate pH 5.5 = 0.40 mol/L each, ± 0.5 mmol (0.10 mol/L) hexanal / (E)-2-heptenal / (E,E)-2,4-decadienal, sealed 15 mL vial, 90 °C, 0.5-5 h; free Cys, free Xyl, TTCA, Cys-Amadori by HPLC-ELSD, pH, A294, A420 time-resolved; volatiles at 3 h as raw peak areas)

### A SECOND, SAME-LAB, SAME-RECIPE RECORD OF THE CAO 2017 BLANK POT, WITH A DIFFERENT TTCA / Cys-Amadori ASSAY (ELSD, not LC-MS SIM). It agrees with Cao on the 0.5-1 h TTCA plateau (0.355 vs 0.356 mol/L) and on free Cys / Xyl to ± 0.01 mol/L, but its Cys-Amadori is 30-40 % lower, its late TTCA 0.04-0.05 mol/L higher, and its A420 ~80× smaller. Not a duplicate of Cao's data; a sibling run. No MFT/FFT time course (3 h peak areas only).

**Source on disk:** `data/articles/zhao2017b.pdf` (7 pp., journal pp. 63-69). Born-digital, Chinese, with
English abstract and bilingual figure/table captions. Prose and methods from `pdftotext -layout`; **Fig. 2a-g
digitised from 400-dpi crops** (scratchpad `img/zhA4-3.png` = 2a, `zhB-4.png` = 2b, `zhC-4.png` = 2c,
`zhD-4.png` = 2d, `zhE-4.png` = 2e, `zhFG-4.png` = 2f-g); Table 1 cross-checked on a 300-dpi crop
(`zhT1-5.png`). Translations are mine.

## 0. Identity

| field | value |
|---|---|
| Title (Chinese) | 脂肪醛对"半胱氨酸-木糖"美拉德反应进程及肉香物质形成的影响 |
| Title (English, as printed) | "Effect of Aliphatic Aldehydes on Maillard Reaction between Cysteine and Xylose and Meaty Flavor Formation" ✔ |
| Authors | ZHAO Jian, ZHAO Mengyao, CAO Changchun, ZHEN Dawei, WANG Jiaxin, LI Li, **XIE Jianchun*** (corresponding, xjchun@th.btbu.edu.cn) |
| Affiliations | Beijing Advanced Innovation Center for Food Nutrition and Human Health / Beijing Laboratory for Food Quality and Safety / Beijing Key Laboratory of Flavor Chemistry, BTBU; Zhen: Taiyuan University of Technology |
| Venue | 食品科学 Food Science **38**(20), 63-69 (2017) ✔; received 2017-03-23 |
| DOI as printed | **10.7506/spkx1002-6630-201720010** in the citation block (p. 64) ✔ matches the brief. The line under the keywords on the same page prints **"DOI:10.7506/spkx1002-6630-2017200010"** (one extra 0) — a typesetting slip; the citation-block form is the one consistent with the article number 1002-6630(2017)20-0063-07. |
| Funding | NSFC 31371838, 31671895 (same grants as Cao 2017); Beijing NSF 6172004; 2017YFD0400106 |
| Same lab | `cao2017_extraction.md` (Cao is third author here), `hou2017_extraction.md`, `zhao2019_extraction.md` |

## 1. Why it matters

Same chemistry, same lab, same charge, buffer, vessel, reactor block and sampling times as Cao 2017's
blank arm — so the blank arm here is an **independent re-run of the TTCA formation pot** (four species, six
times), measured with a different TTCA/Cys-Amadori assay. It also adds three aldehyde arms that show
lipid aldehydes (0.10 mol/L) take cysteine away from xylose: heptenal lowers 0.5 h TTCA from 0.355 to 0.255
mol/L and free Cys from 0.081 to 0.033 mol/L.

## 2. Methods as they matter to a model (translated, §1.3)

- **Charge (§1.3.1).** "In a 15 mL pressure-resistant sealed tube, add 5 mL buffer (NaOH added dropwise to
  NaH2PO4 solution, **0.2 mol/L, pH 5.5**), L-cysteine **2 mmol (0.242 g)**, D-(+)-xylose **2 mmol (0.300 g)**.
  Then add **0.5 mmol** of each of the three aldehydes separately: (E)-2-heptenal (0.061 g), hexanal (0.050 g),
  (E,E)-2,4-decadienal (0.076 g); react at **90 °C for 0.5, 1, 2, 3, 4 and 5 h**." → **0.40 mol/L Cys, 0.40 mol/L
  Xyl, 0.10 mol/L aldehyde** (5 mL basis; solid volume ignored). "A blank experiment was run at the same time.
  Each system in triplicate, mean reported."
- **Vessel.** 15 mL high-temperature pressure tube (Synthware, USA) in a Parallel Synthesis Poly-BLOCK 4 reactor
  (HEL, UK); 5 mL liquid, ~10 mL headspace, atmosphere not stated (air). Stirring not stated (Cao: stirred).
  Quench: "rapid cooling in running water".
- **pH (§1.3.2).** "0.5 mL of reaction solution was diluted to 3 mL and pH measured at room temperature
  (25 °C)" — i.e. **pH of a 6-fold dilution**, not of the pot (Flag 3).
- **A294, A420.** 1 cm cell, diluted into range, back-multiplied by the dilution factor; n = 3, SD < 5 %.
- **HPLC-ELSD (§1.3.3).** Waters XBridge Amide 4.6 × 150 mm, 3.5 µm; 0.01 mol/L ammonium formate (pH 6) :
  ACN = 25 : 75, isocratic, 1 mL/min, 10 min, 45 °C; SEDEX 75 ELSD, drift tube 40 °C, nebuliser air 3.5 bar,
  gain 5; 0.45 µm filter; 2 µL. **Standards of Cys-Amadori and TTCA prepared after ref. 13** (Cao et al.
  2016, Food Science 37(2):63-71). External-standard calibrations (y peak area, x mol/L, all R² > 0.99):
  - TTCA (2-threityl-thiazolidine-4-carboxylic acid) **y = 49 119x − 4 516**
  - Cys-Amadori **y = 24 024.1x − 143.9**
  - xylose **y = 61 719x − 587.1**
  - cysteine **y = 38 046x − 311.57**
  n = 3, SD < 5 %. Calibration ranges NOT stated.
- **SPME-GC-MS (§1.3.4), 3 h samples only.** Reaction mixture in a 15 mL vial, 30 °C, 10 min equilibration + 30
  min headspace, 75 µm CAR/PDMS fibre; HP-5MS and DB-Wax; **no internal standard**; amounts as TIC peak area
  (×10⁵), mean of 3 replicates on both columns, SD < 15 %.

## 3. Tables and figures re-typed

### Fig. 2 (图2, "Changes in the concentrations of xylose (a), cysteine (b), pH (c), 2-threityl-4-carboxythiazolidine (d) and Cys-Amadori (e), and absorbance values at 294 nm (f) and 420 nm (g) in the reaction mixtures with or without aliphatic aldehydes as a function of reaction time") — digitised [D]

Legend: 空白 blank (∇), 反-2-庚烯醛 (E)-2-heptenal (□), 己醛 hexanal (△), 反,反-2,4-癸二烯醛 (E,E)-2,4-decadienal (○).
t = 0: Cys = Xyl = 0.40 mol/L, TTCA = Cys-Amadori = 0, pH 5.5 in all arms (plotted).

Axis calibrations (400-dpi crops): 2a y 0 → px 506, 0.4 → px 106 (1000 px per mol/L); 2b 0 → 457, 0.4 → 56
(1002 px/mol/L); 2c pH 3.5 → 437, 5.5 → 47 (195 px/pH); 2d 0 → 453, 0.4 → 46 (1018 px/mol/L); 2e 0 → 453,
0.12 → 44 (3408 px/mol/L); 2f 0 → 433, 20 → 27; 2g 0 → 968, 3 → 616. **Read-off uncertainty:** ± 0.004 mol/L
(2a, 2b), ± 0.03 pH (2c), ± 0.005 mol/L (2d), ± 0.002 mol/L (2e), ± 0.3 (2f), ± 0.05 (2g). Where markers
overlap (2a, 2b, 2c at 0.5-1 h) the series assignment is by line continuity and carries ± 0.006 mol/L / ± 0.05 pH.

| species | arm | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|---|
| **Xyl** (mol/L) | blank | 0.074 | 0.069 | 0.067 | 0.065 | 0.064 | 0.060 |
| | heptenal | 0.094 | 0.091 | 0.089 | 0.090 | 0.085 | 0.079 |
| | hexanal | 0.083 | 0.078 | 0.078 | 0.078 | 0.075 | 0.075 |
| | decadienal | 0.081 | 0.080 | 0.079 | 0.073 | 0.072 | 0.066 |
| **Cys** (mol/L) | blank | 0.081 | 0.080 | 0.078 | 0.070 | 0.065 | 0.062 |
| | heptenal | 0.033 | 0.033 | 0.033 | 0.030 | 0.029 | 0.028 |
| | hexanal | 0.055 | 0.052 | 0.050 | 0.046 | 0.044 | 0.043 |
| | decadienal | 0.049 | 0.047 | 0.045 | 0.042 | 0.041 | 0.038 |
| **pH** (6× diluted, 25 °C) | blank | 4.51 | 4.47 | 4.43 | 4.28 | 4.24 | 4.10 |
| | heptenal | 4.45 | 4.40 | 4.25 | 4.05 | 3.91 | 3.82 |
| | hexanal | 4.47 | 4.45 | 4.37 | 4.32 | 4.29 | 4.12 |
| | decadienal | 4.47 | 4.45 | 4.23 | 4.13 | 3.95 | 3.80 |
| **TTCA** (mol/L) | blank | 0.355 | 0.333 | 0.327 | 0.304 | 0.296 | 0.290 |
| | heptenal | 0.255 | 0.241 | 0.233 | 0.226 | 0.217 | 0.211 |
| | hexanal | 0.297 | 0.299 | 0.296 | 0.274 | 0.266 | 0.261 |
| | decadienal | 0.298 | 0.286 | 0.282 | 0.260 | 0.249 | 0.242 |
| **Cys-Amadori** (mol/L) | blank | 0 (n.d.) | 0.037 | 0.058 | 0.070 | 0.079 | 0.092 |
| | heptenal | 0 | 0.037 | 0.058 | 0.072 | 0.080 | 0.092 |
| | hexanal | 0 | 0.031 | 0.049 | 0.055 | 0.065 | 0.075 |
| | decadienal | 0 | 0.039 | 0.071 | 0.092 | 0.098 | 0.110 |
| **A294** (back-diluted) | blank | ~0.3 | 0.7 | 1.3 | 2.1 | 3.9 | 5.1 |
| | decadienal | 1.0 | 2.8 | 7.0 | 11.0 | 13.0 | 17.8 |
| **A420** (back-diluted) | blank | ~0.03 | 0.12 | 0.24 | 0.40 | 0.72 | 0.94 |
| | heptenal | 0.07 | 0.33 | 0.69 | 0.97 | 1.47 | 1.74 |
| | hexanal | ~0.01 | 0.03 | 0.19 | 0.27 | 0.44 | 0.68 |
| | decadienal | 0.17 | 0.57 | 1.37 | 2.03 | 2.46 | 3.26 |

Text (§2.1.1): Cys-Amadori "was detected only after 0.5 h" (2e plots 0 at 0.5 h); TTCA "rises rapidly within
0.5 h, then declines, showing that conversion to Cys-Amadori exceeds formation".

### Table 1 (表1 添加、不添加3种脂肪醛反应3 h产物检测到的挥发性物质, "Volatile compounds found in the 3 h reaction mixtures with or without three aliphatic aldehydes by GC-MS") — Maillard-derived block only

Peak area ×10⁵ (TIC, mean of HP-5MS and DB-Wax; no internal standard); "—" not detected.

| compound | RI HP-5MS | RI DB-Wax | blank | hexanal | heptenal | decadienal |
|---|---|---|---|---|---|---|
| 2-methylfuran | 624 | 908 | 134.0 | — | — | — |
| thiophene | 683 | 1018 | 4.7 | — | — | — |
| 2-methylthiophene | — | 1087 | 9.3 | — | — | — |
| 3-methylthiophene | 776 | — | 3.9 | — | — | — |
| **2-methyl-3-furanthiol (MFT)** | 864 | 1303 | **1 133.0** | **13.7** | — | — |
| **2-furfurylthiol (FFT)** | 910 | 1425 | **33.9** | — | — | — |
| 2-thiophenethiol | 964 | — | 10.4 | 0.6 | — | — |
| bis(2-methyl-3-furyl) disulfide | 1521 | 2155 | 68.0 | 3.5 | — | — |
| bis(2-furfuryl) disulfide | — | 2003 | 10.7 | — | — | — |
| subtotal | | | 1 407.8 | 17.8 | — | — |
| **total, all classes** | | | 1 407.8 | 421.4 | 15 298.4 | 15 389.5 |

The aldehyde-derived block (alkylthiophenes, alkanethiols, alkylfurans; e.g. 2-hexylthiophene 3 700 and
2-pentylfuran 6 320 with decadienal; 2-propylfuran 12 815 with heptenal) is in the PDF and not re-typed here
(no engine species).

## 4. Arithmetic (mine) [Z]

**1. Apparent TTCA ⇌ Cys + Xyl quotient, blank.** K_app = [TTCA]/([Cys][Xyl]) = **59.2 / 60.3 / 62.6 / 66.8 /
71.2 / 78.0 L/mol** at 0.5 / 1 / 2 / 3 / 4 / 5 h. Early values match Cao 2017's ≈ 60 L/mol; here the quotient
drifts up by 30 % as pH falls (4.5 → 4.1), where Cao's stayed at 58-65. Read-error propagation ± 8 %.

**2. Cysteine pool closure, blank.** Cys + TTCA + Cys-Amadori = **0.436 / 0.451 / 0.463 / 0.444 / 0.440 / 0.444
mol/L** vs 0.400 charged (+9 to +16 %, same over-closure as Cao's +6 to +17 %, although the assays differ).
Xylose pool: 0.429-0.452. The pool does not measurably drain over 0.5-5 h (within ± 0.02), so the
cysteine-destroying channels at 90 °C, pH 4.1-4.5 are **≤ ~2 × 10⁻⁴ min⁻¹** pooled (≤ 0.02/0.4 = 5 % over 270 min)
— about 2× tighter than Cao’s ≤ 4 × 10⁻⁴.

**3. Cys-Amadori / TTCA, blank:** 0.11 / 0.18 / 0.23 / 0.27 / **0.32** at 1-5 h (Cao: 0.56 at 5 h). Net
TTCA → Cys-Amadori floor: 0.055 mol/L gained over 1-5 h against mean TTCA ≈ 0.31 → **~7 × 10⁻⁴ min⁻¹**
(Cao: ~1.0 × 10⁻³).

**4. Aldehyde arms take cysteine.** Heptenal, 0.5 h: Cys + TTCA = 0.288 mol/L vs blank 0.436, a deficit of
**0.148 mol/L against 0.10 mol/L heptenal added** — larger than 1:1, even allowing the blank's +9 % over-closure
(deficit ≈ 0.11-0.15). Hexanal deficit 0.084, decadienal 0.089 (≤ the 0.10 aldehyde). The aldehydes therefore
bind ~1 equivalent of cysteine within 30 min (thiazolidine / hemithioacetal / Michael adduct; not measured).

## 5. Comparison with Cao 2017 (the brief's question)

Identical recipe: 2 + 2 mmol in 5 mL 0.2 M phosphate pH 5.5, 15 mL Synthware tube, HEL Poly-BLOCK 4, 90 °C,
0.5/1/2/3/4/5 h, n = 3; **identical Cys and Xyl calibration equations** (38046x − 311.57; 61719x − 587.1) —
the same instrument calibration. Different: TTCA and Cys-Amadori by **HPLC-ELSD here** (TTCA y = 49119x −
4516) vs **LC-MS SIM in Cao** (y = 72.398x + 0.613); pH of a 6× dilution here vs pH meter (unspecified
dilution) in Cao; fibre 75 µm CAR/PDMS here vs CAR/PDMS/DVB in Cao; no internal standard here.

Blank arm, this paper minus Cao (Cao values from `cao2017_extraction.md`):

| species | 0.5 h | 1 h | 2 h | 3 h | 4 h | 5 h |
|---|---|---|---|---|---|---|
| Cys (Zhao / Cao) | 0.081 / 0.076 | 0.080 / 0.068 | 0.078 / 0.066 | 0.070 / 0.064 | 0.065 / 0.062 | 0.062 / 0.056 |
| Xyl | 0.074 / 0.080 | 0.069 / 0.077 | 0.067 / 0.073 | 0.065 / 0.071 | 0.064 / 0.070 | 0.060 / 0.069 |
| TTCA | 0.355 / 0.356 | 0.333 / 0.338 | 0.327 / 0.290 | 0.304 / 0.268 | 0.296 / 0.254 | 0.290 / 0.236 |
| Cys-Amadori | 0 / 0.006 | 0.037 / 0.063 | 0.058 / 0.090 | 0.070 / 0.106 | 0.079 / 0.109 | 0.092 / 0.133 |
| pH | 4.51 / 4.67 | 4.47 / 4.56 | 4.43 / 4.49 | 4.28 / 4.33 | 4.24 / 4.26 | 4.10 / 4.15 |
| A420 at 5 h | **0.94** / ≈ 78 | | | | | |

**Verdict: not a duplicate of Cao's data set; a sibling run (or the same samples re-assayed) from the same lab.**
Cys and Xyl agree within 0.012 mol/L (≈ 2-3× read error — consistent with either the same samples or a
replicate run); TTCA agrees at 0.5-1 h and diverges by 0.04-0.05 mol/L after 2 h; Cys-Amadori is 30-40 % lower
throughout; **A420 at 5 h differs ~80×** (0.94 vs ≈ 78) and A294 ~20× (5 vs ≈ 100), which a re-assay of the same
samples with the same spectrophotometric method could not produce — one of the two papers mis-states the
dilution back-correction, or the samples differ. My probabilities: separate physical runs 70 %, same samples
re-assayed 25 %, other 5 %. **For weighting: do not count the two blank arms as two independent pots** —
same lab, recipe, calibration and (likely) operator; treat them as replicates whose spread (TTCA ± 0.025,
Cys-Amadori ± 0.02 mol/L at 3-5 h) is an honest estimate of the between-run + between-assay error.

## 6. What the model can take

| step / observable | quantity | value | unit | conditions | grade |
|---|---|---|---|---|---|
| Cys + Xyl ⇌ TTCA | blank time courses of Cys, Xyl, TTCA, Cys-Amadori | Fig. 2 table | mol/L | 0.4 M each, 90 °C, pH 5.5 → 4.1, 0.5-5 h | **F** (digitised; pair with Cao as replicate) |
| TTCA ⇌ Cys + Xyl | K_app | 59-78 (≈ 60 at ≤ 1 h) | L/mol | 90 °C, pH 4.1-4.5 | **Z** |
| all cysteine-destroying channels | first-order ceiling | ≤ ~2 × 10⁻⁴ | min⁻¹ | 90 °C, pH 4.1-4.5 | **Z** (bound) |
| TTCA → Cys-Amadori | net floor | ~7 × 10⁻⁴ | min⁻¹ | 90 °C, 1-5 h | **Z** |
| aldehyde + Cys | ~1 eq. Cys bound by 0.10 M aldehyde in 30 min | — | — | 90 °C | **Z**, qualitative check for any lipid-aldehyde/thiol trap |
| pH trajectory | Fig. 2c | — | pH | **6× diluted, 25 °C** | **!** not the pot pH |
| MFT, FFT, disulfides | 3 h peak areas, no IS | Table 1 | area ×10⁵ | 3 h only | **!** level-only; no time course, no calibration |
| A294, A420 | Fig. 2f-g | — | AU | | **!** conflicts with Cao by ~80× |

## 7. Flags

1. **DOI printed twice, differently** (…-201720010 in the citation block; …-2017200010 under the keywords).
2. **"TTCA"**: the Chinese text calls it 2-木糖化噻唑烷-4-羧酸 (literally "2-xylosylated thiazolidine-4-carboxylic acid"); the
   English caption calls it "2-threityl-4-carboxythiazolidine" — the same compound as Cao's TTCA.
3. **pH is measured on a 6-fold dilution at 25 °C**; for a 0.2 M phosphate / 0.4 M cysteine pot this shifts pH
   by tenths. Do not score the engine's pot pH against it without a dilution correction.
4. **Pool over-closure (+9-16 %)** under a different assay from Cao's suggests the TTCA standard purity or
   response is biased high in both labs' standards (both made after the same lab's ref. 13 / de Roos route).
5. **No calibration ranges.** The 0.40 mol/L t = 0 values are the charge, not measurements.
6. **No MFT/FFT kinetics.** Volatiles are 3 h raw TIC areas, not comparable to Cao's relative areas
   (blank MFT 1 133 × 10⁵ here vs relative area 0.38 in Cao — different normalisation, not convertible).
7. **The aldehydes kill MFT and FFT completely at 3 h** (heptenal, decadienal: not detected; hexanal: MFT 1.2 %
   of blank) while they raise A294/A420 — a qualitative constraint on any lipid-aldehyde × sulfur coupling.
8. Registry: paper absent from `data/keys/papers.yml`; TTCA and Cys-Amadori absent from `compounds.yml`
   (as noted in the Cao dossier).

## Use in the model

- **TTCA formation / equilibrium (the missing Cys + PENT ⇌ TTCA step).** The blank arm is a second fit row
  (or hold-out) for a reversible condensation at 90 °C; it confirms Cao's K_app ≈ 60 L/mol at ≤ 1 h and
  shows the quotient rising with falling pH (59 → 78), which a speciation-aware K should reproduce. Treat
  Cao + Zhao blank arms as **one pot with two replicates**, not two pots.
- **Cysteine-sink ceiling at 90 °C, pH ~4.3: ≤ ~2 × 10⁻⁴ min⁻¹** (pool conserved for 4.5 h). Together with
  `li2011_extraction.md` (glucose, pH 7, K_cys(90 °C) ≈ 4 × 10⁻⁴ min⁻¹ by extrapolation) this is a direct
  check against any trunk that consumes cysteine on a seconds-to-minutes timescale.
- **No MFT/FFT time course** — nothing here for MFT/FFT formation or loss kinetics.
- **Lipid-aldehyde arms**: a qualitative check (≈ 1 eq. cysteine trapped by 0.1 M enal within 30 min; MFT/FFT
  suppressed to n.d.) for the lipid lane's coupling to sulfur; not a fit row (no adduct measured).

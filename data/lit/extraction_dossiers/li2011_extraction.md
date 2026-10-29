# Li, Zeng, Chen & Fan 2011 — EXTRACTION (glucose + L-cysteine, nominally 0.2 mol/L each, in pH 7 phosphate buffer, 100 / 110 / 120 / 130 °C in an oil bath, 10-120 min; free glucose by HPLC-RI, free cysteine by HPLC after DNFB derivatisation, A420, pH; first-order loss constants for both substrates, zero-order constants for A420 and pH, Arrhenius Ea for all four)

### A MEASURED FREE-CYSTEINE LOSS CONSTANT FOR A GLUCOSE + CYSTEINE POT AT FOUR TEMPERATURES: k = 0.0011 / 0.0032 / 0.0081 / 0.0203 min⁻¹ at 100 / 110 / 120 / 130 °C, so the cysteine half-life is 10.5 h / 3.6 h / 86 min / 34 min, and Ea = 121 kJ/mol. That is the timescale check the trunk needs: at 100 °C, pH 7, 0.2 mol/L, free cysteine lasts HOURS, not seconds.

**Source on disk:** `data/articles/li2011.pdf` (3 pp., journal pp. 73-75). The text layer is unusable (broken
CJK font encoding: `pdftotext` returns symbol soup). Everything below was read from 150-dpi page renders
(scratchpad `img/li-1.png`, `li-2.png`, `li-3.png`) with Read; Table 1 was re-checked on a 300-dpi crop
(`img/liT1-3.png`), Figs. 1-2 on a 300-dpi crop (`img/liF12-2.png`). The paper is in Chinese with an English
title/abstract; translations below are mine and marked as such.

## 0. Identity

| field | value |
|---|---|
| Title (Chinese, as printed) | 葡萄糖/L-半胱氨酸美拉德反应动力学的研究 |
| Title (English, as printed) | "Kinetics of Maillard reaction in glucose/L-cysteine model systems" ✔ |
| Authors | 李伶俐 LI Ling-li, 曾茂茂 ZENG Mao-mao, 陈洁* CHEN Jie (corresponding), 范柳萍 FAN Liu-ping ✔ |
| Affiliation | State Key Laboratory of Food Science and Technology, School of Food Science and Technology, Jiangnan University, Wuxi 214122 |
| Venue | 食品工业科技 (Science and Technology of Food Industry) **32**(10), 73-75 (2011) ✔ |
| Article no. as printed | 文章编号 1002-0306(2011)10-0073-03; CLC TS201.2; document code A |
| DOI as printed | **none printed**. Do not invent one. |
| Received | 2010-09-13 |
| Paper type | Short kinetic study; one pH (7), four temperatures, one concentration; all constants are the authors' SAS regressions |

## 1. Why it matters

The engine's trunk consumes free cysteine 10³-10⁴× faster than measured pots. This paper is a direct
measured counterweight: **glucose + cysteine at 0.2 mol/L, pH 7, 100 °C — free cysteine falls by only ~12 % in
2 h** (0.137 → 0.120 mol/L, Table 1 fit; Fig. 2 points agree). Even at 130 °C the half-life is 34 min. It also
gives an Arrhenius slope for that loss (121 kJ/mol) that I reproduce from the four printed k (arithmetic 2).

## 2. Methods as they matter to a model (translated, §1.2)

- **Charge.** "Glucose and cysteine were dissolved in pH 7 phosphate buffer so that **both concentrations were
  0.2 mol/L**" (§1.2.1, my translation). **Buffer concentration NOT stated**; salts listed: NaH2PO4, Na2HPO4
  (§1.1). Volume NOT stated.
- **Vessel.** 反应钢柱, "reaction steel column(s)" (Jiangsu Hanbang Technology) heated in a 恒温油浴锅
  thermostatted oil bath (Wuxi Hujing) — i.e. sealed steel tubes; fill and headspace NOT stated. Atmosphere NOT
  stated (air presumed).
- **T, t.** **100, 110, 120, 130 °C**; **10, 20, 30, 60, 90, 120 min**; quench "immediately in an ice bath".
  Heat-up time of the steel tubes NOT stated (no come-up correction).
- **Glucose.** HPLC (Waters 1525) after an amino-SPE clean-up cartridge; Thermo ASP-2 Hypersil amino column
  250 × 4.6 mm, 5 µm, 30 °C, ACN : water 70 : 30, 1 mL/min, **refractive-index detector** (示差检测器). Method
  after Wilson & Work (ref. 6).
- **Cysteine.** HPLC (Shimadzu) **after derivatisation, method of Lü Yingguo et al. (ref. 7: 2,4-dinitrofluoro-
  benzene pre-column derivatisation)**. Note: DNFB also derivatises secondary amines, so whether the
  glucose-cysteine thiazolidine is read as "cysteine" is unknown (Flag 3).
- **pH.** Mettler-Toledo SevenEasy meter; measurement temperature NOT stated.
- **Colour.** A420 on a UV2800 (Unico) after dilution; "a larger absorbance means more brown product"; values
  are reported back-calculated (they exceed 18).
- **Replicates.** Error bars on Figs. 1-4; n NOT stated. Constants by SAS regression on Eq. (1)
  A = A0·exp(−Kt) (first order) and Eq. (2) A = A0 + Kt (zero order); Arrhenius Eq. (3) K = K0·exp(−Ea/RT),
  R = 8.314 J/(mol·K).

## 3. Tables re-typed

### Table 1 (表1 不同反应温度条件下各反应动力学参数, "Kinetic parameters of each reaction at different reaction temperatures")

Units, inferred (not printed): A0 in mol/L for glucose and cysteine, absorbance units for A420, pH units for pH;
**K in min⁻¹ (first order) for glucose and cysteine, A420·min⁻¹ and pH·min⁻¹ (zero order)** — the min⁻¹ reading is
confirmed because A0·exp(−K·120) reproduces the 120-min points of Figs. 1-2 (e.g. cysteine 130 °C:
0.1739·e^(−0.0203·120) = 0.0152 mol/L; Fig. 2 last point ≈ 0.015-0.020).

| T (°C) | Glc A0 | Glc K | Glc R² | Cys A0 | Cys K | Cys R² | A420 A0 | A420 K | A420 R² | pH A0 | pH K | pH R² |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.1194 | 0.0013 | 0.9452 | 0.1373 | 0.0011 | 0.9398 | −0.018 | 0.0028 | 0.9626 | – | – | – |
| 110 | 0.1206 | 0.0032 | 0.9232 | 0.1323 | 0.0032 | 0.9840 | −0.2601 | 0.0208 | 0.9786 | 6.3104 | −0.0014 | 0.9380 |
| 120 | 0.1232 | 0.0073 | 0.9660 | 0.1486 | 0.0081 | 0.9814 | −0.5005 | 0.0554 | 0.9968 | 6.4295 | −0.0074 | 0.9842 |
| 130 | 0.1236 | 0.0151 | 0.9904 | 0.1739 | 0.0203 | 0.9803 | −1.1886 | 0.1613 | 0.9964 | 6.3349 | −0.0142 | 0.9449 |

"–" as printed at 100 °C for pH (text: "at 100 °C the pH change was not obvious").

### Table 2 (表2 美拉德反应体系各反应活化能, "Activation energies of each reaction in the Maillard system")

| variable | Ea (kJ/mol) | R² |
|---|---|---|
| glucose concentration | 87.85 | 0.9890 |
| cysteine concentration | 120.96 | 0.9997 |
| A420 | 164.77 | 0.9751 |
| final pH (终pH) | 149.17 | 0.9468 |

### Text values (§2.2)

A420 at 120 min: **0.3465** (100 °C), **2.12** (110), **5.99** (120), **18.34** (130 °C).

### Fig. 2 (free cysteine vs time) — digitised [D]

Calibration on the 300-dpi crop: y = 0 at px 975, y = 0.16 mol/L at px 677 (1862 px per mol/L); read-off
uncertainty **± 0.005 mol/L** (marker size + error bars). **The x positions of the markers do not match the
sampling times in §1.2.1**: the last markers sit at ≈ 72, 110 and 150 min on the printed axis (ticks 10-130 in
steps of 20), not at 60, 90, 120. Table 1 reproduces the points only with t = 10, 20, 30, 60, 90, 120 min, so
I assign the markers to the method's times (Flag 2).

| T (°C) | 10 min | 20 | 30 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|
| 100 | 0.136 | 0.138 | 0.135 | 0.132 | 0.128 | 0.122 |
| 110 | 0.132 | 0.127 | 0.124 | 0.112 | 0.099 | 0.096 |
| 120 | 0.143 | 0.130 | 0.119 | 0.089 | 0.074 | 0.063 |
| 130 | 0.139 | 0.114 | 0.097 | 0.069 | 0.034 | ~0.018 |

Fig. 1 (glucose) not digitised point by point; it starts at 0.11-0.125 mol/L at 10 min in all four series and ends
at ≈ 0.100 / 0.085 / 0.054 / 0.022 mol/L (± 0.005) at 120 min for 100 / 110 / 120 / 130 °C.

## 4. Arithmetic (mine) [Z]

**1. Cysteine half-life, t½ = ln 2 / K.**

| T (°C) | K_cys (min⁻¹) | K_cys (s⁻¹) | **t½ cysteine** | K_glc (min⁻¹) | t½ glucose |
|---|---|---|---|---|---|
| 100 | 0.0011 | 1.83 × 10⁻⁵ | **630 min = 10.5 h** | 0.0013 | 533 min = 8.9 h |
| 110 | 0.0032 | 5.33 × 10⁻⁵ | **217 min = 3.6 h** | 0.0032 | 217 min |
| 120 | 0.0081 | 1.35 × 10⁻⁴ | **85.6 min** | 0.0073 | 95.0 min |
| 130 | 0.0203 | 3.38 × 10⁻⁴ | **34.1 min** | 0.0151 | 45.9 min |

**2. Arrhenius re-fit of the printed K (least squares on ln K vs 1/T, T in K).**
- cysteine: **Ea = 121.05 kJ/mol**, ln A = 32.22 (A in min⁻¹) — reproduces the printed 120.96. ✔
- A420: 164.90 kJ/mol (printed 164.77) ✔; pH (110-130 °C only): 149.28 (printed 149.17) ✔.
- **glucose: 102.4 kJ/mol, printed 87.85 ✘.** The four printed glucose K do not give 87.85 kJ/mol by any
  pairing I tried (end-points 100/130 °C: 102.2). The printed glucose Ea is internally inconsistent with Table 1
  (Flag 1); prefer the re-fit, or prefer Table 1 K directly.

**3. Extrapolation to 90 °C (cysteine), using my fit:** K(90 °C) = 0.0011·exp[−(121050/8.314)(1/363.15 −
1/373.15)] = **3.8 × 10⁻⁴ min⁻¹, t½ ≈ 31 h**. For comparison the pooled cysteine-destroying ceiling read from
Cao 2017 (`cao2017_extraction.md`, cysteine + xylose 90 °C, pH 4.2-4.6) is ≤ ~4 × 10⁻⁴ min⁻¹ — the same order.
(Different sugar and pH; agreement is order-of-magnitude, not a validation.)

**4. Apparent second-order constant**, if the loss is written k₂[Glc][Cys]: k₂ ≈ K_cys / [Glc]_mean. At
100 °C with [Glc] ≈ 0.11 mol/L: **k₂ ≈ 1.7 × 10⁻⁴ L mol⁻¹ s⁻¹**; at 130 °C with [Glc] ≈ 0.06: ≈ 5.6 × 10⁻³
L mol⁻¹ s⁻¹. Lumped over every cysteine sink (Strecker, thermolysis, thiazolidine if not read as Cys), pH
drifting 6.3 → 4.6 at 130 °C.

**5. The missing 30-40 % at t = 10 min.** Charged 0.2 mol/L each; fitted A0 is 0.119-0.124 (glucose) and
0.132-0.174 (cysteine), and the 10-min points are 0.11-0.14 at every temperature including 100 °C, where the
subsequent loss is only 0.1 % per min. So **~0.06-0.08 mol/L of each substrate is "gone" before the first
sample at all four temperatures**, about equally for glucose (0.076-0.081) and cysteine (0.026-0.068).
Two readings: (a) a fast, temperature-insensitive 1:1 sequestration of glucose + cysteine (the
2-polyhydroxyalkyl-thiazolidine-4-carboxylic acid, which forms within minutes, cf. Cao 2017 TTCA 89 % by 30
min with xylose) that neither assay reads as free substrate; (b) the nominal 0.2 mol/L was not what was in
the tube. The paper does not discuss it. My weight: (a) 55 %, (b) 30 %, other (assay bias) 15 %.

## 5. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | grade |
|---|---|---|---|---|---|
| free cysteine loss, glucose pot | first-order K | 0.0011 / 0.0032 / 0.0081 / 0.0203 | min⁻¹ | 100/110/120/130 °C, pH 7 (drifting to 6.3-4.6), 0.2 M nominal each | **M** (authors' fit, Table 1) |
| free cysteine loss | Ea | 120.96 (re-fit 121.05) | kJ/mol | 100-130 °C | **M** ✔ reproduced |
| free cysteine | t½ | 630 / 217 / 86 / 34 | min | as above | **Z** |
| free glucose loss | first-order K | 0.0013 / 0.0032 / 0.0073 / 0.0151 | min⁻¹ | as above | **M** |
| free glucose loss | Ea | 87.85 printed; **102.4 from Table 1** | kJ/mol | 100-130 °C | **!** inconsistent |
| A420 | zero-order K | 0.0028 / 0.0208 / 0.0554 / 0.1613 | A·min⁻¹ | 110-130 °C; 100 °C weak | **M** (dilution-corrected A) |
| A420 | Ea | 164.77 | kJ/mol | | **M** ✔ |
| pH | zero-order K | −0.0014 / −0.0074 / −0.0142 | pH·min⁻¹ | 110/120/130 °C; 100 °C "not obvious" | **M** |
| pH | Ea | 149.17 | kJ/mol | 110-130 °C | **M** ✔ |
| free cysteine time course | 6 points × 4 T | Fig. 2 table | mol/L | | **F** (digitised ± 0.005) |
| early substrate deficit | ~0.06-0.08 | mol/L | all T, by 10 min | **!** unexplained (arith. 5) |

## 6. Flags

1. **Printed glucose Ea (87.85) does not follow from the printed glucose K (102.4 kJ/mol).** Cysteine, A420 and
   pH Ea do reproduce. Use the K.
2. **Fig. 1-4 x-axis positions do not match the stated sampling times**; Table 1 is consistent with the stated
   times (10-120 min). Digitised values are assigned to the stated times.
3. **"Free cysteine" may or may not include the thiazolidine.** DNFB derivatisation reads primary and secondary
   amines; the thiazolidine ring NH is secondary. If the thiazolidine is read as cysteine, K_cys is the loss of
   (Cys + thiazolidine), i.e. even slower true free-Cys turnover is not implied but the sink is lumped.
   Unresolvable from the paper.
4. **Buffer concentration, volume, fill, atmosphere, replicate number and pH-measurement temperature are not
   reported.** pH A0 is 6.31-6.43, not 7: the buffer did not hold pH 7 even at t ≈ 0, and pH fell to ~4.6 by
   120 min at 130 °C (Fig. 4). The constants are for a **drifting** pH.
5. **Steel tubes in an oil bath, no come-up time.** At 10-min first sampling the come-up (likely 1-3 min) biases
   the 130 °C constant low by a few per cent at most.
6. **The English abstract says "L-cystine"** once (typo for cysteine; the Chinese and every table say 半胱氨酸).
7. The introduction promises "different initial pH" conditions; only pH 7 is reported.
8. Registry: paper absent from `data/keys/papers.yml` (not checked against `compounds.yml`).

## Use in the model

- **Timescale check on the trunk's cysteine consumption (the purpose of this read).** Score the engine's
  free-cysteine trajectory for glucose + cysteine 0.2 M, pH 7, at 100-130 °C against Table 1 / Fig. 2: the
  measured t½ is **10.5 h at 100 °C** and **34 min at 130 °C**. An engine that consumes cysteine in seconds to
  minutes at 100 °C is off by the 10³-10⁴ the brief describes; this pot makes that a measured statement. Use
  K_cys(T) as a **check / hold-out**, or the Fig. 2 points as a **fit row** (24 points, ± 0.005 mol/L), but
  decide first how to treat the 10-min deficit (arith. 5): score the *slope* (relative loss from the 10-min
  value), not the absolute level, unless the engine carries a fast thiazolidine pre-equilibrium.
- **Ea_cys = 121 kJ/mol** is usable as a measured lumped temperature dependence of free-cysteine loss in a
  hexose pot at pH ~6-7 (100-130 °C).
- **Not usable:** the glucose Ea as printed; anything about pH dependence (one pH); any product identity
  (no volatiles, no H2S).

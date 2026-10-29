# Deng et al. 2021 — EXTRACTION (purified methionine-glucose Amadori compound, MG-ARP 0.2 mol/L, ± DL-Met 0.1 mol/L, water with initial pH 5.5 or 7.5 (NaOH), sealed vials, 100 / 120 / 130 °C, 0-240 min; plus 0.1 M phosphate pH 7.5 at 100 °C; MG-ARP and Met by HPLC-ELSD, 3-DG / 1-DG / glucosone / GO / MGO by OPD-HPLC, A420, pH (buffer arm only))

### FED-ARP DECAY AT THREE TEMPERATURES AND TWO pH, WITH THE AMINO ACID RELEASED AND FIVE DICARBONYLS — all in figures except the Met-yield table. No rate constant or Ea is printed. From the digitised Fig. 3 (mine): **initial (0-30 min) first-order MG-ARP loss k = 0.0048 / 0.016 / 0.032 min⁻¹ at pH 5.5 and 0.012 / 0.037 / 0.066 min⁻¹ at pH 7.5 (100 / 120 / 130 °C) → Ea ≈ 78 kJ/mol (pH 5.5) and ≈ 72 kJ/mol (pH 7.5)**; pH 7.5 is 2.1-2.4× faster than pH 5.5. **Added Met (0.1 M) does NOT change the ARP decay rate** (k within ± 15 %, inside read error) but raises MGO (to ~1.5× at 100 °C buffered), lowers C6 dicarbonyls and consumes more Met. **3-DG dominates (peak ~14 mM ≈ 7 mol % of ARP at pH 5.5; ≤ 3.5 mM at pH 7.5)**; glucosone up to ~10.6 mM at pH 7.5 / 100 °C; MGO ≤ 1.7 mM (≤ 2.1 with Met); GO ≤ 0.085 mM. **No methional** (that is the sister paper Deng 2022).

**Source on disk:** `data/articles/deng2021.pdf` (11 pp., publisher PDF from pubs.acs.org, clean text
layer). Prose and Table 1 from `pdftotext -layout`. **Fig. 3 digitised** from 300 dpi crops (scratchpad
`lit5/den_f3a-05.png`, `den_f3b-05.png`, read error ± 3 mmol/L); **Fig. 6** from a 200 dpi crop
(`den_f6-08.png`); **Fig. 5** read only for peaks from a 110 dpi page render (`den-07.png`, coarse:
± 5 % of each panel's full scale). Fits: scratchpad `lit5/deng_fit.py` (maillard env). Supporting
Information (S1 fragmentation, S2 NMR) **not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Methionine on the Thermal Degradation of N-(1-Deoxy-D-fructos-1-yl)-methionine Affecting Browning Formation" |
| Authors | Shibin Deng, Heping Cui, Khizar Hayat, Shahzad Hussain, Muhammad Usman Tahir, Yun Zhai, Qiang Zhang, **Xiaoming Zhang*** (Jiangnan), **Chi-Tang Ho*** (Rutgers) |
| Venue | J. Agric. Food Chem. **69** (17), 5167-5177, 2021. Received 7 Apr 2021, accepted 13 Apr 2021, published 23 Apr 2021 |
| DOI as printed | **10.1021/acs.jafc.1c02023** ✔ matches the brief. **Identity verified.** |
| Funding | National Key R&D Program of China 2017YFD0400105 (+ Fujian plan S202011498008) |
| **Relation to `deng2022_extraction.md`** | **Deng 2022 (Food Chem. 382:132500) is the sister paper**, same group, same MG-ARP preparation (Met + Glc 0.2 M, pH 7.5, vacuum dehydration 95 °C, Dowex purification), same 100/120/130 °C, 30-240 min, **unbuffered initial pH 7.5**. Deng 2022 compares MG-ARP with free Met + Glc and measures **methional and pyrazines** (HS-SPME, semi-quantitative); **this paper (2021) measures the ARP itself, released Met, the C6 and C2/C3 dicarbonyls and colour, ± added Met, at pH 5.5 and 7.5**, and has **no volatiles**. Deng 2022's dossier lists this paper as "not on disk"; it now is. The two do not share any plotted series, so no cross-check of numbers is possible except the ARP-loss timescale (both: most ARP gone within 60 min at 120-130 °C). |
| Naming | MG-ARP = MA (figures) = N-(1-deoxy-D-fructos-1-yl)-methionine (C11H21O7NS, 311 Da; four isomers by NMR); MAM = MG-ARP + Met; 3-DG, 1-DG (1-deoxyglucosone), glucosone, GO (glyoxal), MGO (methylglyoxal) |

## 1. Why it matters

The brief: ARP degradation rate constants vs temperature (Ea), pH, buffer, dicarbonyl yields (3-DG,
GO, MGO), methional, and whether added Met changes rates. For the engine: a fed hexose-ARP decay with its
amine release and deoxyosone split is the trunk's ARP → deoxyosone step on a sulfur amino acid; and the
"free amino acid does not accelerate ARP decay" result bears on whether ARP degradation should carry an
amine-catalysis term.

## 2. Methods as they matter to a model

- **ARP:** prepared by Cui et al.'s vacuum-dehydration method (Met + Glc 0.2 M each, pH 7.5, 25 mbar,
  95 °C, 20 min), Dowex 50WX4 (H⁺) purification, NH4OH elution, freeze-dried; **96 % purity**
  (UPLC-Q-TOF), identity by MS/MS and ¹H/¹³C NMR (S2).
- **Model reactions:** MG-ARP **0.2 mol/L** ± Met **0.1 mol/L** in **deionised water**, pH adjusted to
  **5.5 ± 0.1 or 7.5 ± 0.1 with 6 M NaOH**, "sealed vials" (volume, headspace, atmosphere **not
  stated**; continuous stirring), **100 ± 1, 120 ± 1, 130 ± 1 °C**, 0, 30, 60, 120, 180, 240 min, ice
  quench. **Buffer arm:** same in **0.1 M phosphate pH 7.5, 100 °C only** (the text calls these
  "unbuffered models" — a typo, Flag 1). Met-alone control 0.1 M (no colour).
- **pH during heating:** measured **only in the buffer arm** (Fig. 6a); the unbuffered arms' pH drift
  is **not reported**, although the authors note pH drops more from higher initial pH.
- **ARP and Met:** HPLC-ELSD, XBridge Amide, external calibration in log-log form (ARP lg y = 1.3434 lg x
  + 5.6697, R² 0.9968, LOD 0.18 mM, LOQ 0.34 mM; Met lg y = 1.1571 lg x + 5.6220, R² 0.9897, LOD 0.067,
  LOQ 0.27 mM); samples diluted 5×.
- **Dicarbonyls:** OPD (0.5 % in 0.1 M phosphate pH 7.0 + 11 mM DTPA), 4 h dark, RT; HPLC-DAD 315 nm.
  External standards for **3-DG, GO, MGO**; **1-DG and glucosone semi-quantified on the 3-DG curve**.
  GO LOD 0.0086, MGO 0.0042, 3-DG 0.0015 mM.
- **Met yield** (Table 1) = measured Met / ARP reacted (for MAM, the added 0.1 M is excluded).
- Triplicates, mean ± SD.

## 3. Table 1 (verbatim) — Met yield (%) = Met released and still present / MG-ARP reacted

| model | t (min) | 100 °C pH 5.5 | 100 °C pH 7.5 | 120 °C pH 5.5 | 120 °C pH 7.5 | 130 °C pH 5.5 | 130 °C pH 7.5 |
|---|---|---|---|---|---|---|---|
| MG-ARP | 30 | 67.06 ± 9.02 | 86.73 ± 7.4 | 90.64 ± 1.56 | 66.26 ± 2.20 | 58.40 ± 3.30 | 51.86 ± 1.79 |
| | 60 | 76.85 ± 7.92 | 50.90 ± 3.84 | 72.46 ± 1.36 | 67.12 ± 2.96 | 60.05 ± 2.67 | 51.20 ± 2.02 |
| | 120 | 88.89 ± 7.62 | 49.15 ± 2.17 | 65.56 ± 1.57 | 42.46 ± 1.74 | 50.21 ± 1.64 | 49.95 ± 2.36 |
| | 180 | 88.67 ± 3.99 | 47.14 ± 2.86 | 51.24 ± 2.17 | 36.06 ± 2.03 | 36.76 ± 2.64 | 43.25 ± 2.09 |
| | 240 | 82.60 ± 4.22 | 49.16 ± 2.25 | 37.68 ± 1.61 | 34.33 ± 1.11 | 31.64 ± 1.23 | 32.51 ± 1.52 |
| MG-ARP/Met | 30 | 60.02 ± 4.94 | 51.00 ± 3.89 | 35.71 ± 3.16 | 33.13 ± 2.25 | 50.57 ± 3.42 | 59.80 ± 2.52 |
| | 60 | 47.30 ± 7.24 | 48.24 ± 1.85 | 34.93 ± 2.52 | 33.71 ± 2.35 | 34.83 ± 1.85 | 27.49 ± 2.16 |
| | 120 | 54.29 ± 3.48 | 51.54 ± 2.69 | 30.67 ± 1.05 | 33.01 ± 2.18 | 24.50 ± 1.11 | 13.99 ± 1.64 |
| | 180 | 59.80 ± 5.11 | 51.64 ± 1.48 | 23.62 ± 2.73 | 30.28 ± 2.43 | 15.77 ± 2.07 | 10.14 ± 1.88 |
| | 240 | 74.97 ± 4.96 | 40.26 ± 1.61 | 21.57 ± 1.92 | 26.65 ± 1.21 | 15.20 ± 1.84 | 3.40 ± 0.90 |

(Significance letters omitted.)

## 4. Figures digitised

### Fig. 3 — MG-ARP (mmol/L), unbuffered, start 200 [D], ± 3 mmol/L

| arm | 30 | 60 | 120 | 180 | 240 min |
|---|---|---|---|---|---|
| ARP, 100 °C, pH 5.5 | 173 | 157 | 135 | 123 | 111 |
| ARP, 100 °C, pH 7.5 | 141 | 93 | 57 | 43 | 31.5 |
| ARP, 120 °C, pH 5.5 | 124 | 68 | 34 | 20.5 | 14 |
| ARP, 120 °C, pH 7.5 | 67 | 39 | 22 | 14 | 10 |
| ARP, 130 °C, pH 5.5 | 77 | 32 | 12 | 8 | 6 |
| ARP, 130 °C, pH 7.5 | 28 | 12 | 5 | 4 | 3.5 |
| ARP + Met, 100 °C, pH 5.5 | 178 | 156 | 139 | 125 | 115 |
| ARP + Met, 100 °C, pH 7.5 | 140 | 86 | 59 | 45 | 38 |
| ARP + Met, 120 °C, pH 5.5 | 126 | 77 | 38 | 24 | 19 |
| ARP + Met, 120 °C, pH 7.5 | 72 | 34 | 21 | 16 | 10 |
| ARP + Met, 130 °C, pH 5.5 | 69 | 32 | 10 | 6 | 4 |
| ARP + Met, 130 °C, pH 7.5 | 21 | 10 | 4 | 3 | 3 |

### Fig. 6 — 0.1 M phosphate, initial pH 7.5, 100 °C [D] (MA = ARP alone; MAM = ARP + Met)

| t (min) | 0 | 30 | 60 | 120 | 180 | 240 |
|---|---|---|---|---|---|---|
| pH (both) | 7.5 | 7.3 | 6.9 | 6.2 | 5.95 | 5.6 |
| ARP, MA / MAM (mM) | 200 | 105 / 110 | 57 / 71 | 44 / 43 | 33 / 36 | 30 / 31 |
| Met, MA (mM) | 0 | 94 | 126 | 117 | 118 | 102 |
| Met, MAM (mM; 100 added) | 100 | 172 | 213 | 200 | 182 | 166 |
| MGO, MA / MAM (mM) | 0 | 0.095 / 0.10 | **0.155 / 0.227** | 0.142 / 0.250 | 0.140 / 0.236 | 0.133 / 0.206 |
| GO, MA / MAM (mM) | 0.005 | 0.07 / 0.065 | 0.078 / 0.08 | 0.017 / 0.024 | 0.012 / 0.019 | 0.011 / 0.014 |
| 3-DG, MA / MAM (mM) | 0.15 | 1.5 / 1.4 | 1.43 / 0.87 | 0.94 / 0.94 | 0.83 / 0.56 | 0.83 / 0.54 |
| glucosone, MA / MAM (mM) | 0.7 | 1.23 / 1.07 | 0.73 / 0.61 | 0.42 / 0.28 | 0.30 / 0.22 | 0.20 / 0.13 |
| 1-DG, MA / MAM (mM, right axis) | 0.028 | 0.064 / 0.057 | 0.041 / 0.039 | 0.037 / 0.024 | 0.026 / 0.010 | 0.027 / 0.006 |
| A420, MA / MAM | 0.01 | 0.17 / 0.18 | 0.26 / 0.39 | 0.32 / 0.48 | 0.50 / 0.63 | 0.56 / 0.73 |

### Fig. 5 — dicarbonyl peaks, unbuffered [D, coarse]

| species | ARP alone | ARP + Met |
|---|---|---|
| 3-DG | **pH 5.5: ~14 mM** (100 °C still 14 at 240 min; 120 °C peak at 30-60, 5.3 at 240; 130 °C peak at 30, ~1.5 at 240). pH 7.5: ≤ 3.5 mM | pH 5.5: ~13.5 mM peaks, same shape; pH 7.5 ≤ 2.5 |
| 1-DG | pH 5.5: 1.55 (120 °C, 30 min), 1.25 (130 °C, 30), 1.05 (100 °C, 120); pH 7.5 ≤ 0.25 | pH 5.5: ≤ 1.3; pH 7.5 ≤ 0.2 |
| glucosone | **pH 7.5 / 100 °C: 10.6 mM at 60 min**; others ≤ 3 | pH 7.5 / 100 °C: 7.4 at 120 min; others ≤ 2.8 |
| GO | ≤ 0.085 mM (130 °C pH 7.5, 60 min) | ≤ 0.073 mM (120 °C pH 5.5) |
| MGO | pH 5.5: 1.7 mM (130 °C, 60 min), 1.2 (120 °C, 180); pH 7.5 ≤ 0.7 | pH 5.5: **2.05** (130 °C, 60), 1.85 (120 °C, 180); pH 7.5 ≤ 1.28 |

## 5. Arithmetic (mine)

1. **ARP decay, apparent first order, initial 0-30 min (unbuffered):**

   | | 100 °C | 120 °C | 130 °C | Ea |
   |---|---|---|---|---|
   | pH 5.5, ARP | 0.0048 | 0.0159 | 0.0318 min⁻¹ | **78 kJ/mol** |
   | pH 7.5, ARP | 0.0117 | 0.0365 | 0.0655 min⁻¹ | **72 kJ/mol** |
   | pH 5.5, ARP + Met | 0.0039 | 0.0154 | 0.0355 | 91 |
   | pH 7.5, ARP + Met | 0.0119 | 0.0341 | 0.0751 | 75 |

   Read error on the 30-min point (± 3 mM) moves the 130 °C pH 7.5 k by ± 0.004 and each Ea by
   ± ~8 kJ/mol. Three temperatures only. **Use Ea ≈ 75 ± 10 kJ/mol.** The full 0-240 min ln-fit is not
   first order (rapid early loss, slow tail — consistent with the unbuffered pH falling), giving
   meaningless Ea of 27-32 at pH 7.5; do not use it.
2. **pH effect:** k(7.5)/k(5.5) = 2.4, 2.3, 2.1 at 100, 120, 130 °C → **~0.17 dex per pH unit** (initial
   pH; true in-situ pH unknown).
3. **Buffer effect:** at 100 °C, initial pH 7.5, the 0.1 M phosphate arm loses ARP at 0.020-0.022 min⁻¹
   (0-30 min) vs 0.012 unbuffered — **~1.8× faster in phosphate** although its pH falls only to 7.3 in that
   window. Either the unbuffered pot acidifies fast within 30 min, or phosphate catalyses ARP
   enolisation (known for glycine ARP in the literature; not tested here).
4. **Added Met does not change ARP decay**: MA vs MAM initial k differ by −19 … +15 % with no consistent
   sign, inside read error; same conclusion in phosphate (Fig. 6b). The authors' conclusion stands.
5. **Met balance, 130 °C pH 7.5, 240 min:** ARP alone: 196.5 mM reacted, 32.5 % → 64 mM Met left, **133 mM
   consumed**. ARP + Met: total available 196.5 + 100; yield 3.4 % → ~107 mM Met left, **~189 mM consumed**.
   Extra Met is consumed roughly additively (Strecker/dicarbonyl reactions), not catalytically.
6. **Amine release is fast:** buffer arm, 60 min: 143 mM ARP lost, 126 mM Met present → **88 %** of the
   reacted ARP returned free Met (the rest consumed or not yet released).
7. **Dicarbonyl yields (mol % of ARP charged):** 3-DG ~7 % (pH 5.5), ≤ 1.8 % (pH 7.5); glucosone ≤ 5.3 %
   (pH 7.5, 100 °C — oxidative, so headspace O2 matters and is not reported); 1-DG ≤ 0.8 %; MGO ≤ 0.85 %
   (≤ 1.0 % with Met); GO ≤ 0.04 %. All are transient levels (net of formation and loss), i.e. lower bounds
   on flux.

## 6. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| MG-ARP → products | initial first-order k | 0.0048 / 0.016 / 0.032 | min⁻¹ | 0.2 M, water, initial pH 5.5, 100 / 120 / 130 °C | Fig. 3a [D]; arith. 1 | derived_from_digitised |
| same | initial first-order k | 0.012 / 0.037 / 0.066 | min⁻¹ | initial pH 7.5 | Fig. 3a [D] | derived_from_digitised |
| same | Ea | ≈ 72-78 (± 10) | kJ/mol | 100-130 °C, 3 T | arith. 1 | derived |
| same, phosphate | initial k | 0.020-0.022 | min⁻¹ | 0.1 M phosphate pH 7.5, 100 °C | Fig. 6b [D] | derived_from_digitised |
| effect of + 0.1 M Met on ARP decay | ratio | 1.0 ± 0.15 | — | all arms | Fig. 3, 6b; arith. 4 | derived |
| Met released | yield | Table 1 | % | all arms | Table 1 | printed |
| 3-DG / 1-DG / glucosone / GO / MGO | peak levels | sec. 4 tables | mM | all arms | Figs. 5, 6 [D] | digitised (1-DG, glucosone semi-quantitative) |
| pH trajectory | pH vs t | 7.5 → 5.6 in 240 min | — | phosphate 0.1 M, 100 °C | Fig. 6a [D] | digitised |
| methional, volatiles, H2S | — | **NOT MEASURED** (see Deng 2022) | — | — | — | — |

## 7. What is NOT reported

- No printed rate constant or Ea for anything.
- No pH trajectory in the unbuffered arms (where all the temperature data are).
- No methional, methanethiol, DMDS or any volatile; no H2S.
- Vial volume, fill, headspace gas — glucosone and GO formation are oxidative, so this matters.
- No mass balance for the ARP carbon (dicarbonyls account for ≤ 10 % at their peaks).

## 8. Flags

1. Methods call the phosphate arm "unbuffered models" (p. 5170) — evidently a typo for "buffered".
2. 1-DG and glucosone are quantified on the 3-DG calibration (no standards) — semi-quantitative.
3. ARP calibration is log-log (lg y = 1.34 lg x + …), so relative error grows at low concentration; the
   < 10 mM tail points are less reliable than the error bars suggest.
4. Fig. 3a/b markers at 120/130 °C pH 5.5 and 120 °C pH 7.5 overlap at 30-60 min; series assignment
   checked against line styles, ± 1 series possible at 120 °C pH 7.5 vs 130 °C pH 5.5 at 60 min.
5. The authors invoke retro-aldol of 3-DG → MGO "facilitated by Met's thioether"; it is an interpretation
   of net levels, not a measured flux.

## 9. Use in the model

- **ARP decay barrier:** a measured ~75 kJ/mol (100-130 °C) for a hexose-amino-acid ARP's total loss,
  pH 5.5-7.5, well below the carbohydrate-class median of 103.7 used for B55/B61 single-temperature
  constants. If the engine's ARP → deoxyosone steps (`r_arp_dpo` and its TDP sibling on the sulfur lane;
  the trunk's ARP step) sit on a higher shared barrier, this is a candidate anchor — with the caveat that
  it is a glucose ARP of methionine, unbuffered (in-situ pH unknown), and total loss, not a single exit.
- **pH law:** ~0.17 dex per pH unit (5.5 → 7.5), far shallower than a specific-base (1 dex/pH) law; ARP
  enolisation is not hydroxide-catalysed in this range. Do not transfer B61's hydroxide law to ARPs.
- **No amine catalysis of ARP decay:** 0.1 M free Met leaves the ARP loss rate unchanged → an ARP decay
  step should not carry a free-amino-acid factor. Free amino acid acts downstream (dicarbonyl
  consumption, MGO formation, browning).
- **Deoxyosone split:** pH 5.5 → 3-DG (1,2-enolisation) dominant at ~7 mol %; pH 7.5 → 3-DG falls ~4×,
  1-DG stays low (consumed fast), glucosone rises (oxidative) — qualitatively the engine's acid-favoured
  1,2- vs base-favoured 2,3-enolisation, here on a hexose ARP.
- **Phosphate catalysis:** the ~1.8× faster decay in 0.1 M phosphate at 100 °C, pH 7.5 is a (weak, one
  point) data point for the phosphate-catalysis hypothesis named in `SULFUR_FLUX_MAP.md` sec. 5 as a
  candidate low-pH MFT source (Hofmann's 0.5 M phosphate arms).
- **Methional / Met chain:** nothing here; use Deng 2022.

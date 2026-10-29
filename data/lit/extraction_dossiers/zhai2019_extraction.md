# Zhai et al. 2019 (JAFC) — EXTRACTION (D-xylose 165 + L-cysteine 82.6 mmol/L in water, pH 7.4 set with NaOH, 70 / 90 / 110 / 130 °C for 0-160 min: TTCA and Cys-ARP TIME-RESOLVED; purified TTCA 5 mmol/L ± Cys 5 mmol/L at 120 °C, pH 7.4: 1-DX, 3-DX, glyoxal 10-60 min; a two-step "cysteine-indicator" browning test)

### THE PARENT OF THE ZHAI TTCA SERIES, AND ITS ONLY TEMPERATURE LADDER OF TTCA *FORMATION*. Four temperatures × nine times of TTCA and the cysteine Amadori compound from a free xylose + cysteine pot — the engine has no step that forms TTCA, so this is a hold-out it cannot yet be scored on. The same paper also makes the mechanism claim that added cysteine "completely restor[es] TTCA to sugar and amino acid", which runs against mass action (Flag 1).

**Source on disk:** `data/articles/zhai2019.pdf` (9 pp., typeset JAFC version, downloaded via McMaster
Univ. 23 Oct 2019). Prose read from `pdftotext` (raw and `-layout`); **Table 1 is an image** and was
re-typed from a 300 dpi render; **Figs. 1A, 1B, 3A and 3B were digitised from 400-500 dpi renders**
(scratchpad `img/zhai19_f1-3.png`, `zhai19_f3A-5.png`, `zhai19_f3b-5.png`). Figs. 2, 4 and 5 are
chromatograms / schemes with no quantitative axis beyond detector mV. Every digitised value is
tagged **[D]** with its read error. The Supporting Information (TIC and MS/MS of TTCA and ARP) is
**not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Interaction of Added L-Cysteine with 2-Threityl-thiazolidine-4-carboxylic Acid Derived from the Xylose−Cysteine System Affecting Its Maillard Browning" |
| Authors | Yun Zhai, Heping Cui, Khizar Hayat, Shahzad Hussain, Muhammad Usman Tahir, Jingyang Yu, Chengsheng Jia, **Xiaoming Zhang** (corresponding, xmzhang@jiangnan.edu.cn), **Chi-Tang Ho** (corresponding, ho@aesop.rutgers.edu) |
| Affiliations | State Key Laboratory of Food Science and Technology, Jiangnan University (Wuxi); King Saud University (Riyadh); Rutgers University |
| Venue | J. Agric. Food Chem. **2019, 67 (31), 8632-8640**. Received 12 Jul 2019, revised 15 Jul 2019, accepted 16 Jul 2019, published 16 Jul 2019 |
| DOI as printed | **10.1021/acs.jafc.9b04374** ✔ matches the brief |
| Funding | NSFC 31671826; Jiangnan Postgraduate JNKY19_009; National Key R&D 2017YFD0400105; JUFSTR20180204; King Saud RG-1440-020 |
| Abbreviations as defined | **TTCA** 2-threityl-thiazolidine-4-carboxylic acid; **ARP** Amadori rearrangement product (here = **Cys-Amadori**, N-(1-deoxy-D-xylulos-1-yl)-L-cysteine — Table 1 labels it "Cys-Amadori"); MRI Maillard reaction intermediate; 1-DX 1-deoxypentose-2-ulose; 3-DX 3-deoxypentose-2-ulose; GO glyoxal |
| Children on disk | `zhai2020_extraction.md` (spray-drying TTCA/ARP ratio), `zhai2021_extraction.md` (TTCA ± extra Xyl, the kinetic fit), `zhai2023b_extraction.md` = `zhai2023jafc_extraction.md` (TTCA ± extra Cys, 90-150 °C), `zhai2023c_extraction.md` (dehydration), `kang2026_extraction.md` |

## 1. Why it matters

**(a) A formation ladder the engine cannot represent yet.** Fig. 3A measures TTCA and the
Cys-Amadori ("ARP") against time at **70, 90, 110 and 130 °C** from a free Xyl + Cys charge. TTCA is
already **14.5 mmol/L at t = 0** [D] (formed during mixing and NaOH adjustment at room temperature),
reaches a plateau of **54-63 mmol/L** at 70-90 °C (65-77 % of the cysteine), peaks at **58.5 mmol/L
by 20 min at 110 °C** then drifts down, and at **130 °C** peaks at ~46 mmol/L by 20 min and falls to
**6.9 mmol/L by 160 min**. The Cys-Amadori stays at **≤ 5.6 mmol/L** throughout. `sulfur.py` has no
Cys + PENT → TTCA step and no Cys-Amadori species, so none of this can be scored today; with a
reversible condensation (see `cao2017_extraction.md`, K_app ≈ 60 L/mol at 90 °C, pH ~4.5) it would
be the natural **four-temperature hold-out**.

**(b) The dicarbonyl split from fed TTCA.** Fig. 3B (TTCA 5 mmol/L, 120 °C, pH 7.4) gives
**3-DX > GO > 1-DX at every time, with 3-DX / 1-DX = 2.1-2.6** [D]. The engine's direct TTCA
degradation `r_ttca_deg` produces **only DPO (1-deoxyosone)**; 3-DX (`TDP`) can reach a TTCA pot only
through `r_ttca_cys` → PENT → `r_pent_tdp` / `r_pent_caramel`, both **acid-catalysed** and slow at
pH 7.4. **Prediction to test (not run here):** the engine will put DPO ≥ TDP in a fed-TTCA pot at
pH 7.4, against a measured TDP/DPO ≈ 2.5. The authors' route is TTCA → Cys-Amadori → 1-DX / 3-DX.

**(c) Added cysteine suppresses all three dicarbonyls** (at 60 min: 3-DX to 20 %, GO to 23 %, 1-DX to
45 % of the TTCA-alone values [D]) — the same perturbation as Zhai 2023b, at one temperature, with
numbers. This is the origin of the "cysteine blocks TTCA" line that `zhai2023b_extraction.md` reads
as evidence that `r_ttca_cys` is reversible.

### How it relates to the companion dossiers

| | **this paper (2019)** | `zhai2021_extraction.md` | `zhai2023b_extraction.md` | `cao2017_extraction.md` |
|---|---|---|---|---|
| start | **free Xyl 165 + Cys 82.6 mM** (2 : 1); and TTCA 5 mM ± Cys 5 mM | purified TTCA 10 mM ± Xyl 10 mM | purified TTCA 10 mM ± Cys 0.1-10 mol/mol | free Cys + Xyl 0.4 M each |
| pH | 7.4 (NaOH), unbuffered | 5.5-8 | 5.5-8.5 | 5.5 phosphate, drifts to 4.15 |
| T | **70 / 90 / 110 / 130 °C**; 120 °C | 100 / 120 / 140 | 90-150 (A420), 100/120/140 | 90 |
| TTCA measured | **yes** (formation, 4 T) | yes (decay) | no | yes (formation) |
| free Cys / Xyl measured | **no** (chromatograms only) | yes | no | **yes** |
| dicarbonyls | 1-DX, 3-DX, GO at 120 °C | 4 at 100 °C | 4 + ARP at 100-140 °C | none |

## 2. Methods as they matter to a model

- **Formation pot.** "The reaction system of Xyl and Cys was mixed in a **2:1 molar ratio** and
  dissolved in **80 mL** of deionized water. The concentration of final solution was **34.8 mg/mL**.
  The solution pH was adjusted to **pH 7.4** with NaOH (6 mol/L)." → Cys **82.6 mmol/L**, Xyl
  **165.2 mmol/L** (arithmetic 1). **Unbuffered; pH never re-measured.**
- **Temperature ladder.** "heated at temperatures of **70 °C, 90 °C, 110 °C, and 130 °C** for 0, 20,
  40, 60, 80, 120, 140, and 160 min" (the figure also has a **100 min** point — Flag 4); ice-quenched,
  stored at 4 °C. **The vessel for this ladder is not stated** (110 and 130 °C in water require a
  sealed pressure vessel; the headspace is unknown — Flag 5).
- **Two-step "indicator" test.** Six 90 °C aliquots for 20-70 min, then **Cys added at "1% of the
  mixture"** (basis not stated — Flag 6), pH re-adjusted to 7.4, then **120 °C for 110 min in an oil
  bath**; A420 of the final MRPs, stepwise samples **diluted 4-fold**, controls undiluted (Flag 7).
- **TTCA–Cys model.** Cys : TTCA = 1 : 1, **Cys 5 mmol/L** (so TTCA 5 mmol/L), **10 mL** deionized
  water, pH 7.4 (NaOH 6 M), **10 min at room temperature**, then **120 °C oil bath, 0-60 min**; ice
  quench. Control: TTCA alone, same.
- **TTCA/ARP purification.** Rotary evaporation 40 °C, ethanol reconstitution, Dowex 50WX4 (H⁺),
  elution water then 0.2 M NH4OH (TTC test for xylose), semi-prep RP-HPLC XBridge Amide; identity by
  UPLC-ESI-Q-TOF MS and ¹³C NMR (Table 1). **Purity by TIC peak area: TTCA 95 %, ARP 92 %.**
- **TTCA/ARP quantitation.** HPLC-ELSD, XBridge Amide 4.6 × 150 mm 3.5 µm, A = 10 mM ammonium
  formate pH 6, B = ACN, gradient 95 % B 0-3 min → 75 % 5 min → 65 % 15 min → 95 % 17 min, 1 mL/min,
  10 µL, 25 °C; **external standards** of the purified compounds. Retention: Xyl 5.65, Cys 6.80,
  **TTCA 8.80, ARP 10.54 min**. **Calibration equations and ranges not printed.**
- **Dicarbonyls.** OPD 2 % in 50 mL HEPES (1 mol/L, pH 7) + DTPA 0.216 g; **hexanedione internal
  standard** (960 ppm, 40 µL) + 2 mL sample + 2 mL reagent; **12 h** in the dark; HPLC-DAD SunFire
  C18 5 µm 4.6 × 150 mm, 0.1 % formic acid / methanol gradient, 0.8 mL/min, 35 °C, 20 µL.
- **Statistics.** "All the experiments in the research process were performed three times";
  mean ± SD; p < 0.05.

## 3. Tables re-typed

### Table 1 (verbatim). Derivative Fragmentation in the LC-MS/MS Spectrum and ¹³C NMR Data of the Two Prepared Intermediates

| Compound | ¹³C NMR / ppm (D₂O) | Ions in MS/MS spectra of [M+H]⁺ (m/z) |
|---|---|---|
| TTCA | 174.40 (s) (C1), 75.07 (d) (C6), 74.42 (d) (C5), 72.92 (d) (C7), 69.77 (d) (C2), 67.48 (d) (C4), 64.83 (t) (C8), 35.60 (t) (C3) | 254.03 [M+H]⁺; 236.02 [M+H−H₂O]⁺; 218.02 [M+H−2H₂O]⁺; 200.02 [M+H−3H₂O]⁺; 188.02 [M+H−2H₂O−CH₂O]⁺ |
| Cys-Amadori, structure 1 (S-cyclised) | 175.00 (s) (C₁), 82.79 (s) (C₅), 76.38 (d) (C₆), 73.20 (d) (C₇), 66.62 (t) (C₃), 58.35 (t) (C₈), 55.51 (t) (C₄), 28.04 (t) (C₃) | 254.03 [M+H]⁺; 236.02 [M+H−H₂O]⁺; 218.02 [M+H−2H₂O]⁺; 200.02 [M+H−3H₂O]⁺; 188.02 [M+H−2H₂O−CH₂O]⁺ (one entry for both structures) |
| Cys-Amadori, structure 2 (pyranose) | 174.60.00 (s) (C₁) [sic], 82.17 (d) (C₄), 78.01 (d) (C₆), 72.93 (d) (C₅), 66.52 (d) (C₇), 62.78 (t) (C₈), 51.40 (d) (C₂), 26.83 (t) (C₃) | (as above) |

Transcription notes: structure 1 lists **C₃ twice** (66.62 t and 28.04 t) and has no C₂; structure
2's first shift is printed "174.60.00". Both are reproduced as printed (Flag 8). The text states the
open-chain ARP was not found by NMR; two six-membered-ring isomers were.

### Digitised Fig. 3A — TTCA and Cys-ARP vs time, Xyl 165 + Cys 82.6 mmol/L, pH 7.4 initial [D]

Axis 0-70 mmol/L; read error **± 0.7 mmol/L (TTCA), ± 0.3 mmol/L (ARP)**. No error bars drawn.

**TTCA (mmol/L)**

| T | 0 | 20 | 40 | 60 | 80 | 100 | 120 | 140 | 160 min |
|---|---|---|---|---|---|---|---|---|---|
| 70 °C | 14.5 | 41.1 | 50.5 | 55.6 | 54.0 | 54.1 | 53.1 | 54.0 | 54.3 |
| 90 °C | 14.5 | 48.3 | 59.2 | 57.7 | 59.0 | 58.3 | 59.1 | 61.9 | 63.2 |
| 110 °C | 14.5 | 58.5 | 55.9 | 53.9 | 52.0 | 52.8 | 48.7 | 42.6 | 39.0 |
| 130 °C | 14.5 | 45.7 | 42.5 | 34.5 | 32.2 | 21.2 | 15.8 | 11.3 | 6.9 |

**Cys-ARP (mmol/L)**

| T | 0 | 20 | 40 | 60 | 80 | 100 | 120 | 140 | 160 min |
|---|---|---|---|---|---|---|---|---|---|
| 70 °C | 0 | 0.2 | 0.5 | 1.0 | 1.3 | 1.6 | 1.9 | 2.1 | 2.8 |
| 90 °C | 0 | 0.9 | 1.4 | 1.7 | 1.9 | 2.6 | 3.0 | 3.4 | 3.7 |
| 110 °C | 0 | 2.3 | 3.7 | 5.6 | 5.3 | 4.8 | 4.6 | 4.2 | 3.8 |
| 130 °C | 0 | 4.3 | 4.8 | 3.9 | 3.1 | 1.6 | 1.2 | 0.0 | 0.0 |

**Fig. 1B** is the 90 °C series re-plotted (identical within read error; it carries small error bars:
TTCA 14.4, 48.0, 59.0, 57.6, 58.8, 58.0, 58.9, 61.5, 63.0 mmol/L).

### Digitised Fig. 3B — α-dicarbonyls from TTCA (5 mmol/L) and TTCA–Cys (5 + 5 mmol/L), 120 °C, pH 7.4 initial [D]

Axis 0-300 µmol/L; read error **± 5 µmol/L**. No 0 min point plotted.

| species (system) | 10 | 20 | 30 | 40 | 50 | 60 min |
|---|---|---|---|---|---|---|
| 3-DX (TTCA) | 89 | 164 | 199 | 205 | 228 | 249 |
| GO (TTCA) | 55 | 95 | 112 | 121 | 138 | 149 |
| 1-DX (TTCA) | 42 | 68 | 76 | 89 | 91 | 98 |
| 3-DX (TTCA+Cys) | 21 | 25 | 30 | 33 | 49 | 51 |
| 1-DX (TTCA+Cys) | 17 | 22 | 31 | 34 | 41 | 44 |
| GO (TTCA+Cys) | 10 | 12 | 18 | 25 | 30 | 34 |

### Digitised Fig. 1A — A420 of the final MRPs (two-step test) [D], read error ± 0.005

| first-step time at 90 °C (min) | 20 | 30 | 40 | 50 | 60 | 70 |
|---|---|---|---|---|---|---|
| with added Cys | 0.617 | 0.608 | **0.430** | 0.511 | 0.594 | 0.683 |
| without added Cys | 0.425 | 0.454 | 0.516 | 0.537 | 0.576 | 0.620 |

### Values printed in the text

| quantity | value | where |
|---|---|---|
| optimal MRI (TTCA) formation | **pH 7.4, 90 °C, 40 min** | Results / Conclusions |
| TTCA and ARP yields from Fig. 2 peak areas | **62 % and 5 %** | Results, "Method validation" |
| "consumption of Cys was up to **62%** at the maximum formation point of TTCA" | — | "Proposed Mechanism" |
| purity (TIC peak area) | TTCA **95 %**, ARP **92 %** | Results |
| [M+H]⁺ | **254.03** (C₈H₁₅NO₆S, 253.03 Da) | Results |

### Arithmetic (mine)

**1. Charge.** 34.8 g/L ÷ (2 × 150.13 + 121.16 g per 2 Xyl + 1 Cys) = **0.0826 mol/L Cys, 0.165 mol/L
Xyl**. The later Zhai papers prepare TTCA from "an equimolar solution of Cys and Xyl (0.0827 mol/L)"
citing this paper — same cysteine concentration, **different xylose** (Flag 3).

**2. TTCA yield vs the printed "62 %".** Plateau TTCA at 90 °C, 40 min = 59.2 mmol/L = **72 %** of the
cysteine; the printed 62 % is a peak-area percentage of Fig. 2 (ELSD response is non-linear) and is
not a molar yield.

**3. Apparent condensation quotient, lower bound.** Free Cys and Xyl are not reported, so take them
by balance (no other cysteine sink): at 90 °C, 40 min, Cys ≈ 82.6 − 59.2 − 1.4 = 22.0, Xyl ≈ 104.6
mmol/L → [TTCA]/([Cys][Xyl]) ≥ **26 L/mol**; at 90 °C, 160 min ≥ **41 L/mol**; at 70 °C, 160 min ≥
**20 L/mol**. **Lower bounds**: any cysteine lost elsewhere lowers free Cys and raises the true quotient.
The final pH is unknown. Cao 2017's measured-free-Cys value at 90 °C and pH ~4.5 is **≈ 60 L/mol**.
Class **derived_assumption**.

**4. Net TTCA loss once formation has stalled.** 130 °C, 20 → 160 min: ln(45.7/6.9)/140 =
**0.0135 min⁻¹**; 110 °C, 100 → 160 min: ln(52.8/39.0)/60 = **0.0050 min⁻¹**. These are net
(formation from residual Cys + Xyl is still running) and carry an unknown, falling pH. For scale,
Zhai 2021's purified TTCA 10 mM at pH 7 has first-order equivalents 0.0065 (120 °C) and 0.0081
min⁻¹ (140 °C).

**5. 3-DX / 1-DX from TTCA, 120 °C:** 2.1 / 2.4 / 2.6 / 2.3 / 2.5 / 2.5 at 10-60 min. Sum of the three
dicarbonyls at 60 min = 496 µmol/L = **9.9 % of the 5 mmol/L TTCA charge**.

**6. Cysteine suppression at 60 min:** 3-DX 51/249 = **0.20**, GO 34/149 = **0.23**, 1-DX 44/98 =
**0.45**. The 1-deoxy branch is suppressed least.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | order as fitted | source | evidence class |
|---|---|---|---|---|---|---|---|
| TTCA formation from free Xyl + Cys | TTCA vs t, 4 T | Fig. 3A table | mmol/L | Xyl 165 + Cys 82.6 mM, pH 7.4 initial, unbuffered, 70-130 °C, 0-160 min | **nothing fitted** | Fig. 3A [D] | **digitised time series** (± 0.7) — **hold-out candidate once a formation step exists** |
| Cys-ARP in the same pots | vs t, 4 T | Fig. 3A table | mmol/L | same | — | Fig. 3A [D] | digitised (± 0.3) — **no engine species** |
| TTCA present before heating | t = 0 | **14.5** | mmol/L | RT mixing + pH adjustment | — | Fig. 3A/1B [D] | digitised |
| TTCA ⇌ Cys + Xyl | quotient lower bound | **≥ 20-41** | L/mol | 70-90 °C, pH ≤ 7.4 | equilibrium | arithmetic 3 | derived_assumption |
| net TTCA loss in the formation pot | first-order | 0.0135 (130 °C); 0.0050 (110 °C, late) | min⁻¹ | as above | 1st, mine | arithmetic 4 | derived_from_digitised |
| 3-DX, 1-DX, GO from fed TTCA | vs t | Fig. 3B table | µmol/L | TTCA 5 mM, 120 °C, pH 7.4 | — | Fig. 3B [D] | **digitised time series** (± 5) — **fit/hold-out candidate on the TTCA → TDP/DPO split** |
| 3-DX / 1-DX from fed TTCA | ratio | **2.1-2.6** | — | 120 °C, pH 7.4 | — | arithmetic 5 | within_study_ratio |
| same with Cys 5 mM added | vs t | Fig. 3B table | µmol/L | 120 °C, pH 7.4 | — | Fig. 3B [D] | digitised |
| Cys suppression of dicarbonyls | ratio at 60 min | 0.20 / 0.45 / 0.23 (3-DX / 1-DX / GO) | — | 120 °C | — | arithmetic 6 | within_study_ratio |
| TTCA fate with added Cys | "TTCA disappeared completely" at 10 min, 120 °C | — | — | 5 + 5 mM | — | Fig. 4B (chromatogram) | **threshold** — see Flag 1 |
| A420, two-step test | Fig. 1A table | absorbance | — | 90 °C then 120 °C 110 min | — | Fig. 1A [D] | level_only (dilution basis inconsistent, Flag 7) |

## 5. What this can and cannot anchor in the model

**Can:**
1. **A four-temperature hold-out for TTCA formation** (Fig. 3A), conditional on adding
   Cys + PENT ⇌ TTCA (and, to score the "ARP" curve, a Cys-Amadori species). Pair with Cao 2017 as
   the fit row: Cao fixes K and k at 90 °C / pH ~4.5 with measured free Cys and Xyl; this paper tests
   70-130 °C at pH 7.4.
2. **A fit or hold-out row on the products of fed TTCA at 120 °C, pH 7.4** (Fig. 3B): six times × three
   dicarbonyls at 5 mmol/L TTCA — the same pot type as Kang 2026 and Zhai 2021, at a temperature
   inside both. The 3-DX / 1-DX ≈ 2.5 ratio directly tests whether `r_ttca_deg` → DPO alone is the
   right topology.
3. **The cysteine-perturbation arm** (Fig. 3B dashed), as a quantitative companion to Zhai 2023b's
   qualitative delay.

**Cannot:**
- **No free cysteine or free xylose is quantified** in any pot (chromatograms only), so the formation
  ladder cannot alone fix an equilibrium constant (only the lower bound of arithmetic 3).
- **No TTCA time course in the 120 °C TTCA / TTCA–Cys pots** (only the 0 and 10 min chromatograms), so
  `k_ttca_cys` / `k_ttca_deg` get no direct row from Fig. 3B — only their products.
- **pH is set once and never measured** after heating; the vessel and headspace of the formation
  ladder are not stated.
- A420 values are on inconsistent dilution bases.

## 6. Flags

1. **The central mechanism claim contradicts mass action.** The paper concludes that "the added Cys
   reacted with TTCA to prevent it from further MR by **completely restoring it to sugar and amino
   acid** under heat treatment", based on Fig. 4B: after 10 min at 120 °C with 5 mM Cys, "the peak of
   TTCA disappeared completely, but the peak of Cys still remained, and a new peak of Xyl appeared".
   For TTCA ⇌ Cys + Xyl, **adding Cys shifts the equilibrium toward TTCA, not away from it**; complete
   loss of TTCA within 10 min with excess cysteine present is the opposite of what the condensation
   equilibrium predicts, and the TTCA-alone control at the same time *retains* most of its TTCA
   (Fig. 4D). Possible readings: (i) a co-elution / ELSD artefact; (ii) a cysteine-driven route
   that consumes TTCA without regenerating it (e.g. trans-thiazolidination or thiol addition to the
   open-chain Schiff base); (iii) the claim is wrong. **Do not encode "Cys converts TTCA back to
   Cys + Xyl" as a forward rate law.** The quantitative result that survives is Fig. 3B (downstream
   dicarbonyls suppressed). Note that `zhai2023b_extraction.md` leans on the same claim from the
   2023 paper; the mass-action reading (Cys addition *stabilises* TTCA, which delays everything
   downstream) explains the 2023 delays without needing this mechanism.
2. **"ARP" here is the cysteine Amadori compound**, not the engine's `ARP` (alanine Amadori, Zhou
   2023). Table 1 names it "Cys-Amadori". Never map this paper's ARP curve onto `ARP`.
3. **Xyl : Cys ratio.** This paper prepares the intermediate at **Xyl : Cys = 2 : 1** (Cys 82.6 mM).
   Zhai 2023b and 2023d describe the same preparation as "an equimolar solution of Cys and Xyl
   (0.0827 mol/L)", citing this paper. One of the two descriptions is wrong; the excess xylose
   matters for any equilibrium reading (arithmetic 3 uses 2 : 1, as printed here).
4. **Time grid.** Methods list 0, 20, 40, 60, 80, 120, 140, 160 min; Figs. 1B and 3A also plot
   **100 min**.
5. **Vessel not stated** for the 70-130 °C ladder (the 120 °C steps use an oil bath; no vessel or fill
   is named anywhere). 110 and 130 °C aqueous runs imply a sealed vessel of unknown headspace.
6. **"Indicator (Cys) was then added (1% of the mixture)"** — mass or molar basis not stated.
7. **A420 dilution bases differ:** stepwise samples diluted 4-fold, controls undiluted, and Fig. 1A
   does not say which basis it plots.
8. **Table 1 typos:** Cys-Amadori structure 1 assigns C₃ twice; structure 2's carboxyl reads
   "174.60.00".
9. **Calibration not printed** for TTCA, ARP (ELSD) or the dicarbonyls; ELSD response is non-linear
   in concentration (Zhai 2023d prints log-log ELSD calibrations for the same method), so the
   mmol/L axis depends on an unseen fit.
10. **Dicarbonyl quantitation basis conflicts across the series:** here an *internal* standard
    (hexanedione) is used; Zhai 2023b and 2023d use *external* standards for the "same" method.
11. **Registry gaps:** TTCA, Cys-Amadori, cysteine, xylose, 1-DX, 3-DX absent from
    `data/keys/compounds.yml`; `glyoxal` present. Paper not in `papers.yml`.
12. **What to request:** free Cys and Xyl concentrations behind Figs. 3A and 4; the final pH of every
    pot; TTCA remaining against time in the 120 °C TTCA and TTCA–Cys pots; ELSD calibration curves;
    the vessel and fill for the temperature ladder.

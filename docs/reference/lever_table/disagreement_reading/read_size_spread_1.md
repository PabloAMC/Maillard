# Why the effect sizes differ: the largest same-direction spreads, part 1

*Reader's note, 2026-11-12. Literature and contrast files only; no engine output; no repository file changed.
Sign convention as in the cell file: effect = log10(arm at the lower lever value, or the "none" arm) minus log10(the other arm), per step. A negative number therefore means "the higher pH, the larger dose or the added compound gives more". Every cell value was recomputed from the keyed arms (the recomputations are given below).*

## What I read

| file | how |
|---|---|
| DISAGREEMENT_READING.md, method_read_cysteine_dose.md | in full |
| method_read_thiols.md | in part (Part A lines on FFT; Parts B and C in full) |
| method_read_interlab.md | not read |
| contrast files | cerny2007, li2021ph, mottram2002, li2022heat, ding2023, hwang1995, martin2001, xu2007, zhou2023, kumazawa2003, grosch1992, feng2022: every contrast behind the eight cells, arm by arm (values, censoring, limits, conditions) |
| mottram2002.pdf | text layer, in full |
| cerny2007.pdf | text layer, in full |
| li2021.pdf | text layer pp. 1-7 (methods, pH, xylose, amino acids, volatiles); p. 4 as a page image (Fig. 1, final pH). Supplement (docx): Table S2 furfural rows only |
| li2022b.pdf | Table 4A/4B furfural rows and results text pp. 8-9. Methods and pH taken from method_read_cysteine_dose.md |
| ding2023.pdf | Table 2 furfural and furfuryl mercaptan rows and the discussion pp. 12-13. Methods taken from method_read_cysteine_dose.md |
| hwang1995.pdf | text layer pp. 1-3; p. 3 as a page image (Table 1, because the text layer drops the ctrl column) |
| martin2001.pdf | text layer pp. 1-4 (methods, Tables 1-4, results) |
| xu2007.pdf | Table 1 (pH rows, FFT, furfural) and results text p. 333 |
| Zhou2023b.pdf | Table 1 FFT row and results text pp. 2475-2476. **Fig. 2 (final pH) not viewed**: the final-pH figures are taken from DISAGREEMENT_READING.md |
| kumazawa2003.pdf | methods, Tables 2-3 and results text in full; introduction and GC-O part skimmed |
| grosch1992.pdf | experimental section and Table II |
| feng2022.pdf | methods (ARP preparation, model reactions, HS-SPME) and Tables 1-2 FFT rows. Supplement not read |

**A finding that recurs: how a censored or non-monotone ladder becomes one number.** I recomputed every cell value from the keyed arms. Four mechanisms in the table, not in the papers, set most of the extreme values:

1. **A censored arm enters at its limit.** For example, Mottram's not-detected arm enters as 0.1, and Xu's as 0.07. The effect is then set by the detection limit.
2. **Two censored arms give a ratio of two limits.** In Grosch, both arms are "<" values (<1.5 and <1.7), and the table still prints -0.05.
3. **A zero-dose arm is dropped and only the longest monotone run is kept.** This is the rule found earlier for the cysteine cells. Here it also drops the not-detected arms of Cerny's pH ladder.
4. **When SDs are printed, the slope is an inverse-variance weighted regression.** This reproduces zhou2023's -0.11 exactly (the unweighted value is -0.20). It is my reconstruction from the arms; I did not read the generator.

---

## Cell 1. Furfural x pH (spread 1.39 dex): cerny2007 -0.12; li2021ph -0.19 to -0.90; mottram2002 -2.93

**Designs.**

| | cerny2007 (Firmenich) | li2021ph (NUS) | mottram2002 (Reading) |
|---|---|---|---|
| system | xylose 150 µmol + cysteine 50 + thiamin 50 µmol in 1 mL **0.5 M K-phosphate** | pork-trimming hydrolysate (~300 mg free AA / 10 mL; cysteine ~3 mg / 10 mL) + xylose 0 - 3.5 g/100 mL, **unbuffered**, pH set with lactic acid | cysteine + ribose 5-phosphate (or ribose, or IMP), 0.3 mmol each in 6 mL **0.2 M pyrophosphate** |
| heat | 145 °C, 20 min | 90 / 95 / 100 °C, 60 min | 140 °C, 30 min |
| pH arms | 4.0, 5.0, 5.5, 6.0, 7.0 (buffer pH) | initial 5.5 vs 4.5 | initial 5.6 vs 4.2 (final 5.5 vs 3.9) |
| measurement | HS-SPME, TIC area, no IS | HS-SPME GC-FID area | Tenax headspace, area vs 1,2-dichlorobenzene (RF 1) |

**Recomputation.**
- **Mottram -2.93.** This is log10(1260 / 0.1) / 1.4. The pH 5.6 Rib-PO4 arm is "−" (below ~0.1 ng/mmol, Table 1 footnote), and it enters at the limit. The value is a **lower bound set by the detection limit**, not a measured ratio.
  - The two sister pairs are keyed with censor "trace" (limit 0.5): ribose 2290 against <0.5, giving ≤ -2.62; IMP 27 against <0.5, giving ≤ -0.88.
  - They do not appear in the cell. I cannot tell from the files I may read whether the table drops "trace" arms.
- **Cerny -0.12.** This is log10(158 / 208), the 4.0 → 5.0 pair only.
  - Table 2 prints 208, 158, 165, **0, 0** at pH 4.0, 5.0, 5.5, 6.0, 7.0. The paper says furfural formed "only at acidic pH, and no traceable amounts" at 7.0 (p. 1553).
  - The two zeros are dropped. Of the two length-2 runs (4 → 5 falling, 5 → 5.5 rising), the first is kept.
  - The keyed slope is therefore the **flat part** of a response that collapses between pH 5.5 and 6.0, where ≥ 2 dex is lost in half a unit (bounded only by detection).
- **li2021ph.** The values reproduce exactly from Table S2. For example, 90 °C, 0.5 % xylose: log10(113.36 / 16.25) = 0.84. Every arm is detected, and each contrast is a two-arm pair.

**Explanations.**
- **(e) Censoring and segment choice set both extremes.**
  - Mottram's magnitude is a detection-limit bound. With the ribose pair it would be ≥ 2.6. With a lower LOD it would be larger.
  - Cerny's magnitude comes from the plateau below his own cliff. If both labs' ladders were scored over the same step (≈ pH 5.5 → 6), both would read "detected → not detected", i.e. ≥ 2 dex per unit.
  - Confidence: **high** that these two numbers are artefacts of censored arms. Medium that the underlying responses agree (one LOD-bounded pair and one ladder).
- **(a) The response is a switch, not a line.**
  - Cerny: flat from 4 to 5.5 (-0.12, then +0.02), then gone by 6.0.
  - Mottram: furfural absent or trace at final pH 5.3-5.7 in every buffered arm except phthalate-ribose (260), and 1260-4990 at final pH 3.9.
  - li2021ph spans final pH ≈ 4.4 → 5.4. That is below Cerny's cliff but well above Mottram's acid arm, and it gives the intermediate -0.5 to -0.9.
  - Medium confidence. The switch point probably depends on buffer and temperature, see (b).
- **(c) li2021ph's pH step is smaller than labelled.**
  - Fig. 1 (page image, read by eye): the pH 5.5 group ends at ≈ 5.4 (0 % xylose) down to 4.86 (A100-35). The pH 4.5 group ends at ≈ 4.45 down to 4.35.
  - The final-pH gap is therefore ≈ 0.95 units at 0 % xylose and ≈ 0.5 at 3.5 % xylose and 100 °C.
  - Per unit of final pH, li's xylose contrasts become ≈ -1.0 (e.g. 100 °C 3.5 %: -0.51 / 0.51; 90 °C 0.5 %: -0.84 / ≈ 0.85). This moves li toward Mottram and away from Cerny. It does not close the gap. Medium-high confidence on the numbers (read off a bar chart).
- **(b) Other factors set differently.**
  1. **Xylose.** Inside li2021ph the effect at 0 % xylose is -0.19 to -0.49 and with xylose -0.5 to -0.9. At 0 % xylose the heated arms (2.09, 1.98, 2.36) are **below** the unheated control (4.13), so these contrasts measure loss of furfural already present, not formation.
  2. **Buffer species.** Mottram's ribose arms at pH 5.6 give trace in phosphate and 260 in phthalate (Table 1). The authors note that furfural "decreased in the presence of buffer". A phosphate-pH contrast is therefore partly a buffer contrast.
  3. **A sink.** Cerny's pot holds cysteine and thiamin. Furfural is the precursor of FFT, which falls in parallel there (431 → 364 → 185 → 0). Mottram's Rib-PO4 arm forms its carbonyls by dephosphorylation to the furanone, not via 3-deoxypentosone (p. 4083), so furfural is suppressed at pH 5.6 for a reason specific to that sugar.
  - Medium confidence for each, low on how much each one contributes.
- **(d) A possible lever-linked offset in li2021ph.** The **unheated** controls already differ by pH: 4.13 at pH 5.5 against 9.41 at pH 4.5 (Table S2), i.e. 0.36 dex.
  - The 0 % xylose heated contrasts (0.19-0.49) are of the same size, so they may carry no reaction signal.
  - Low confidence on the cause (formation during pH adjustment and storage, or recovery); the paper does not discuss it.

**What would confirm or refute.**
- Score all three on a matched final-pH step (e.g. 5.5 → 4.5 final) with censored arms treated as bounds. This predicts consistent ≥ 1 dex per unit near 5.5-6 and smaller slopes below 5.
- Refuted if a buffered ladder through 4.5-6.5 at one temperature shows a constant slope.

**Keying.** No transcription error.
- cerny2007-1 matches Table 2. li2021ph-151 to -165 match Table S2. mottram2002-6/7/8 match Table 1, including the final pH values in the labels.
- **Flag (table, not keying):** mottram2002-6 is an LOD-bounded value and should be scored as a bound. cerny2007-1 loses its two not-detected arms.

---

## Cell 2. Furfural x precursor dose (cysteine) (spread 1.12): li2022heat +0.02, -0.09; ding2023 -0.53

**Designs.**
- **li2022heat.** Pork-trimming extract (NPE) or hydrolysate (EPH) + **xylose 3.5 g/100 mL (233 mM)** + cysteine 0 / 0.5 / 1.0 g/100 mL (0 / 41 / 83 mM). Initial pH 5.5, unbuffered (final ≈ 5.15 NPE, 4.75 EPH, per the cysteine note). 100 °C, 60 min. HS-SPME, area against 2-octanol.
- **ding2023.** Fish-waste hydrolysate + **xylose 0.25 % (16.7 mM)** + cysteine 0 / 0.25 / 0.5 / 0.75 / 1.0 % (0 - 83 mM). Initial pH 7.0, unbuffered. 110 °C, 60 min. HS-SPME, area against 1,2-dichlorobenzene.
- **Cys : xylose (mol/mol).** li 0 / 0.18 / 0.35. ding 0 / 1.2 / 2.5 / 3.7 / 5.0.

**Recomputation.** The zero arm is dropped in all three ladders.

| ladder | arms | step kept | keyed effect |
|---|---|---|---|
| li N35 (Table 4A) | 97.79, 41.87, 44.20 | 0.5 → 1.0 | log10(44.20 / 41.87) = +0.02 |
| li E35 (Table 4B) | 235.14, 101.39, 82.69 | 0.5 → 1.0 | -0.09 |
| ding (Table 2 row 2) | 121.7, 32.33, 9.54, 10.94, 3.30 | 0.25 → 0.5 (first of three length-2 runs) | log10(9.54 / 32.33) = -0.53 |

**Explanations.**
- **(e) + (a): the dropped first dose and the segment chosen.**
  - The **first dose** is similar in both labs: li 0 → 0.5 % gives -0.37 dex in both NPE (97.79 → 41.87) and EPH (235.14 → 101.39); ding 0 → 0.25 % gives -0.58 (121.7 → 32.33).
  - After that the labs part. li is flat (Cys still ≤ 0.35 per xylose). ding keeps falling about 0.5 dex per doubling (0.25 → 1.0 %: log10(3.30 / 32.33) / 2 = -0.50) once cysteine exceeds xylose.
  - The 25-fold spread is mostly the table comparing li's **flat second step** with ding's **steep second step**. Confidence: **medium-high**.
- **(b) The dose relative to the sugar.**
  - In ding the cysteine : xylose ratio passes 1 at the first dose. FFT ("furfuryl mercaptan") appears at 0.75 and 1.0 % (209.8 and 111.6 ng/g, Table 2) while furfural falls to 3.3, which fits furfural being consumed by H2S once sulfur is in excess. The paper itself attributes the fall to cysteine competing with glycine for xylose (p. 13).
  - In li, xylose stays in 3-6-fold excess. FFT appears only in E35-10 at 8.13 ng/mL, far too little to account for the 0.37 dex furfural drop, so the first-dose drop there is some other diversion that the paper does not name.
  - The pH also differs (7.0 against 5.5). Furfural formation is favoured at low pH and H2S release at higher pH.
  - Confidence: medium on the ratio explanation; low on the mechanism.
- **(d) Not the sign; possibly a little of the size.** Both are HS-SPME semi-quantification, with total load falling with cysteine in both (cysteine note), so any fibre-competition bias runs the same way in both.

**What would confirm or refute.**
- Keep the 0 → first-dose step as its own contrast, and band by Cys : sugar. The prediction is agreement (about -0.4 to -0.6 dex) at the first step, and a continued fall only where Cys : sugar > 1.
- Refuted if a ladder at Cys : sugar < 0.5 keeps falling steeply, or one at > 1 goes flat.

**Keying.** li2022heat-81 and -86 and ding2023-1 match Tables 4A, 4B and 2. No error.
- **Flag (table):** the zero-dose step, the most comparable number in both studies, is lost. ding2023's kept pair is one of three tied runs.

---

## Cells 3, 6, 7, 8. Pyrazines x amino acid identity, "none" vs Ile / Asn / Gln: hwang1995 against martin2001

| cell | spread | hwang1995 (none → X) | martin2001 (G → G+X) |
|---|---|---|---|
| 2-methylpyrazine, none vs Ile | 1.05 | 146.77 → 148.95: -0.01 | 6 (a) → 10 (a): -0.22 |
| pyrazine, none vs Ile | 0.90 | 6.71 → 8.43: -0.10 | 1 (a) → 6 (ab): -0.78 |
| 2-methylpyrazine, none vs Asn | 0.85 | 146.77 → 226.53: -0.19 | 6 (a) → 130 (c): -1.34 |
| 2,3-dimethylpyrazine, none vs Gln | 0.80 | 45.68 → 60.66: -0.12 | 1 (a) → 6 (b): -0.78 |

The letters are martin's Table 4 multiple-range groups. Hwang prints no SD or statistics.

**Designs.**
- **hwang1995.** Glucose + [15N]glycine + test amino acid, **equimolar**, in 20 g wheat starch. pH 7, freeze-dried, re-humidified to 12-14 % water, **180 °C for 1 h**. Purge-and-trap at 80 °C, toluene-d8 IS. µg/g glucose. **"ctrl contained glycine only"** (Table 1 footnote).
- **martin2001.** Potato slices leached of sugars and amino acids, then infused with glucose (G), or glucose + one amino acid at **potato-proportional, unequal doses** (Table 2). Deep-fried 2 min (oil 180 → 160 °C). Tenax headspace, area against 1,2-dichlorobenzene. Amounts "to two significant figures", CV < 25 %.
- **Molar doses in martin (uptake per 100 g of slices).** Glucose 258.7 mg = 1.44 mmol. Asn 210 mg = 1.59 mmol (1.1 per glucose). Gln 82.3 mg = 0.56 mmol (0.39). Ile 15 mg = 0.11 mmol (0.08).

**Explanations.**
- **(b) The reference arm is not the same pot. This is the main reason, as in the arginine conflict cell.**
  - Hwang's "none" already holds one amine per glucose. A second amino acid at most doubles the amine, and in this dry, hour-long roast pyrazine formation is not nitrogen-starved: ctrl total 599 µg/g.
  - Martin's "none" is glucose with essentially no amine (the blank W shows only traces of 2- and 3-methylbutanal, p. 3887). The infused amino acid is the nitrogen source.
  - Hence martin's large effects where the dose is substantial: Asn at 1.1 per glucose gives -1.34, the largest. Hwang's equimolar Asn gives only -0.19.
  - Confidence: **high** for the Asn and Gln cells.
- **(d)/(e) The isoleucine cells sit on martin's floor.**
  - Martin's isoleucine dose is 0.08 mol per mol glucose. Its effects rest on 6 → 10 and 1 → 6 peak-area units, printed as integers.
  - A printed "1" spans 0.5 - 1.5, i.e. ±0.18 dex from rounding alone. The paper's own groups do not separate G from GI (2-MP: a, a; pyrazine: a, ab).
  - The size gap in the two Ile cells is therefore within martin's measurement floor, not a laboratory difference. Confidence: medium-high.
- **(b, secondary) Extent of reaction.**
  - Hwang: 1 h at 180 °C in a low-water starch matrix. Martin: 2 min of frying with the oil cooling to 160 °C.
  - Martin's pots are far from completion and limited by nitrogen. Hwang's are long-cooked, where amino acids also compete: the paper reports that arginine depresses glycine's yield and lysine raises it.
  - Medium confidence; untested.
- **(a)** Not applicable (categorical).

**What would confirm or refute.**
- Martin's own pairs that share an amine background, if any were keyed (e.g. GA against the mixed system X), should give hwang-sized effects.
- Hwang's design with glycine removed from the control (glucose + X against glucose) should give martin-sized effects.
- Refuted if a glycine-containing reference still gives ≥ 0.8 dex for Asn.

**Keying.** hwang1995-1, -2, -4 match Table 1, checked against the page image because the text layer drops the ctrl column. martin2001-37, -45, -66 match Table 4. No error.
- **Flag (table):** martin's Ile arms are not significantly different from G in the source. Hwang has no replicate statistics.

---

## Cell 4. FFT x pH (pH 6 - 7.5) (spread 1.01): xu2007 -1.09, -0.46; zhou2023 -0.11; kumazawa2003 -1.18, -1.14

**Designs.**

| | xu2007 (CAU) | zhou2023 (Ho network) | kumazawa2003 (Ogawa) |
|---|---|---|---|
| what is measured | FFT **formed** from ribose + cysteine (0.05 M each) | FFT **formed** from Xyl-Ala ARP + cysteine (20 mM each) | **survival** of 1 ppm FFT spiked into buffer (no sugar, no amino acid), canned without deoxidation |
| buffer | 0.2 M phosphate | none (NaOH), water | McIlvaine citrate-phosphate |
| pH arms (initial / label) | 5.6, 7.0, 8.0 | 6, 7, 8 | Table 2: 6.0 vs 6.8. Table 3: 3.0 - 7.0 |
| final pH | control 5.5 / 6.2 / 6.8; SC-CO2 5.5 / 6.0 / 6.3 (Table 1) | ≈ 3.3 / 3.4 / 5.0 (Fig. 2, via DISAGREEMENT_READING) | buffered, label only |
| heat | 140 °C, 1 h | 120 °C, 1 h | 121 °C 10 min / 123 °C 20 min |
| measurement | HS-SPME, TIC, RF 1 | HS-SPME, external calibration in water | DCM extraction, GC-FID, residual % |

**Recomputation.**
- **xu2007 control, -1.09.** Arms 22.09, ND, ND. The two ND arms enter at the 0.07 limit, and the regression over 5.6 / 7 / 8 gives exactly -1.09. **The value is LOD-set.**
- **xu2007 SC-CO2, -0.46.** Unweighted regression of 25.13, 14.57, 1.79 over 5.6 / 7 / 8.
- **zhou2023, -0.11.** Arms 813.65 ± 13.37, 757.965 ± 13.03, 325.22 ± 14.61. The inverse-variance weighted regression reproduces -0.114. The unweighted value is -0.20. The weights favour the two low-pH arms, where nothing moves (Duncan letters c, c), and down-weight the pH 8 arm, where FFT falls.
- **kumazawa2003, -1.18 / -1.14.** log10(5.5 / 48.3) / 0.8 and log10(5.3 / 43.4) / 0.8.

**Explanations.**
- **(c) + (a): where the reaction pH actually sits on a steep curve. Main.**
  - Kumazawa's own pH ladder (Table 3, relative FFT after heating) is **flat from 3 to 5** (0.46, 0.44, 0.41) and **steep above 5.4** (0.37 → 0.21 → 0.05 → < 0.01 at 5.4, 6.0, 6.4, 7.0, i.e. -0.87 per unit from 5.4 to 6.4).
  - zhou2023's pots ended at pH 3.3 - 5.0: initial 6 and 7 both end near 3.3 - 3.4, so its 6 → 7 step is no step (-0.03, not significant). That puts it on the flat part. Its one real move (7 → 8, final ≈ 3.4 → 5.0) gives -0.37 per initial unit and ≈ -0.23 per final unit.
  - xu2007's final pH 5.5 → 6.8 and Kumazawa's buffered 6.0 → 6.8 sit on the cliff. Per unit of **final** pH, xu's SC-CO2 arm is log10(1.79 / 25.13) / 0.8 ≈ -1.4, close to Kumazawa's -1.1 to -1.2.
  - Confidence: **medium-high.** zhou's final pH comes from an earlier reading of a plotted figure, not checked here.
- **(e) Table artefacts on both ends.** xu2007's -1.09 is set by the 0.07 limit. zhou2023's -0.11 is halved by SD weighting. Confidence high on the arithmetic.
- **(b) Different processes under one lever.**
  - Kumazawa measures **thermal-oxidative loss** of pre-formed FFT. The loss runs mainly to difurfuryl disulfide, which rises 0.07 → 0.17 from pH 5.4 to 7.0 (Table 3), "Fenton-type" with air in the can.
  - xu and zhou measure **net formation**. The paper text of zhou2023 attributes the fall to less furfural at higher pH (p. 2475), and xu's control furfural also goes 3.19 → 1.74 → ND (Table 1).
  - The two processes point the same way, which is why the signs agree. Their sizes need not match.
  - Medium confidence. Kumazawa arguably belongs in its own "stability" sub-cell.
- **(d) A measurement bias that moves with the lever, in xu2007 only.** The thiols note found browned-matrix thiol binding of 0.3 - 1 dex in xu's undiluted, browned control arm; browning rises with pH. This would steepen xu's control slope; the less-browned SC-CO2 arm is shallower (-0.46). Low-medium confidence.

**What would confirm or refute.**
- Rescore on final pH, with censored arms as bounds. The prediction: slopes cluster near -1 per unit wherever the final pH spans 5.5 - 7, and near 0 below 5.
- Refuted if a buffered formation study held at pH 3.5 - 5 shows ≥ 0.5 dex per unit.

**Keying.** xu2007-3/-4 match Table 1 (FFT row and final pH rows). zhou2023-2 matches Table 1. kumazawa2003-1/-2/-5 match Tables 2-3. No error.
- **Flags (table):** xu2007-3 is LOD-set. zhou2023-2 is down-weighted toward its non-moving arms. kumazawa2003 is a spiked-stability study, not a Maillard formation contrast.

---

## Cell 5. FFT x sulfur source, none vs glutathione (spread 0.96): grosch1992 -0.05; feng2022 -0.63, -0.33

**Designs.**
- **grosch1992.** Thiamine 30 µmol + ribose 670 µmol, with or without GSH 230 µmol, in 500 mL 0.2 M pyrophosphate at pH 5.7, refluxed 2 h (100 °C). Isotope dilution (d-FFT).
- **feng2022.** **Ribose-glutathione ARP** (RG-ARP, made from ribose + GSH), 20 mM, with or without 20 mM free GSH. Water, initial pH 7.0, 100 or 120 °C. HS-SPME with authentic standards.
- **Concentrations.** Grosch's pot is about 15-40-fold more dilute: ribose 1.3 mM, GSH 0.46 mM.

**Recomputation.**
- **grosch1992-8.** Both arms are **censored**: thiamine + ribose < 1.5 µg; thiamine + GSH + ribose < 1.7 µg (Table II, nos. 3 and 7). The -0.05 is log10(1.5 / 1.7), **a ratio of two reporting limits**.
- **feng2022.** log10(10.46 / 44.63) = -0.63 (120 °C, 120 min). log10(5.05 / 10.92) = -0.33 (100 °C).

**Explanations.**
- **(e) The Grosch value carries no information.** Neither its size nor its sign is measured. It should be "keyed but uninformative", like the all-ND pairs elsewhere. Confidence: **high**. Without it the cell has one laboratory.
- **(b) The reference arms differ anyway.**
  - Feng's "none" is itself a glutathione Amadori product, so sulfur is already present. Its FFT is 10.46 µg/L without added GSH.
  - Grosch's "none" has sulfur only from thiamine.
  - Feng's own effect also grows with temperature and time: at 120 °C, -0.37, -0.57, -0.61, -0.63 at 30, 60, 90, 120 min; at 100 °C, -0.33. Grosch ran at 100 °C and far more dilute, where Feng's trend predicts a smaller effect.
  - Medium confidence (an inference from Feng's trend; Grosch cannot test it).

**What would confirm or refute.** Not testable from Grosch's data. A Grosch-type dilute pot with a lower LOD, or Feng's design at Grosch's concentrations, would decide it.

**Keying.** grosch1992-8 matches Table II (limits 1.5 and 1.7). feng2022-6/-8 match Tables 1-2. No transcription error.
- **Flag (table):** an effect is computed from two censored arms. This should be blocked generally.

---

## Summary

| cell | spread (dex) | explanation | confidence |
|---|---|---|---|
| Furfural x pH (4.5-6) | 1.39 | e (Mottram LOD bound; Cerny's zeros dropped, plateau kept), a (switch near pH 5.5-6), b (xylose, buffer species, sugar form), c (li's final-pH gap ≈ 0.5-0.95, steepens li), d (li unheated controls differ 0.36 dex) | medium-high for e + a; low-medium for the rest |
| Furfural x cysteine dose | 1.12 | e (zero arm dropped; ties), a (flat vs steep second step), b (Cys : xylose 0.18-0.35 vs 1.2-5; pH 5.5 vs 7) | medium-high (mechanism: low-medium) |
| 2-Methylpyrazine x none vs Ile | 1.05 | b (glycine-containing vs amine-free reference), d (martin 6 → 10 units, not significant; Ile dose 0.08 per glucose) | medium-high |
| FFT x pH (6-7.5) | 1.01 | c + a (zhou's final pH 3.3-5.0 on the flat part; xu and Kumazawa on the cliff), e (xu LOD-set; zhou SD-weighted), b (formation vs spiked-FFT survival), d (xu browning binding) | medium-high |
| FFT x none vs glutathione | 0.96 | e (Grosch: ratio of two "<" limits), b (Feng's reference is a GSH-ARP; T and dilution) | high (e); medium (b) |
| Pyrazine x none vs Ile | 0.90 | b, d (martin 1 → 6 units, groups a / ab) | medium-high |
| 2-Methylpyrazine x none vs Asn | 0.85 | b (reference arm; Asn 1.1 per glucose in martin) | high |
| 2,3-Dimethylpyrazine x none vs Gln | 0.80 | b (d: martin "1" at the floor) | medium-high |

No transcription error in any of the 12 sources. The flags are table-level: LOD substitution, two-censored-arm ratios, dropped zero and not-detected arms, SD weighting, and pooling a stability study with formation studies.

**Overall pattern.**
- The biggest size spreads are mostly made by the table, not the laboratories. A detection limit sets the effect (Mottram, Xu control), two limits are divided (Grosch), a not-detected or zero arm is dropped so that a plateau segment is kept (Cerny, li2022heat), or SD weighting mutes the arm that moves (zhou2023).
- Where the spread is real, it again comes from what the reference arm and the lever meant. A "none" arm that already holds the limiting reagent (glycine, a GSH-ARP) gives small effects. A dose far below or far above the sugar on a molar basis decides whether the slope has saturated.
- pH effects on furfural and FFT behave like a switch near pH 5.5-6.5. A laboratory's slope depends on whether its reaction pH, the final one rather than the label, straddles that switch.

# Leksrisompong, Barbano, Foegeding, Gerard & Drake 2010 — COMPLETE TRANSCRIPTION
### The corpus's only dataset with sensory thresholds AND headspace partition coefficients measured on the SAME eight matrices.

**Full extraction of every number in `data/articles/leksrisompong2010.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

Tables 3 and 5 re-read from a **400 dpi render** to confirm the `E⁵` superscript in the
Table 5 header — the single most consequential character in the paper.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **Pattarin Leksrisompong¹, David M. Barbano², Allen E. Foegeding¹, Patrick Gerard³, MaryAnne Drake¹** |
| Title | **"THE ROLES OF FAT AND PH ON THE DETECTION THRESHOLDS AND PARTITION COEFFICIENTS OF THREE COMPOUNDS: DIACETYL, δ-DECALACTONE AND FURANEOL"** |
| Venue | ***Journal of Sensory Studies* 25 (2010) 347–370** |
| DOI | **10.1111/j.1745-459X.2009.00264.x** |
| Accepted | **24 August 2009** (published in the 2010 volume — so "joss 2009/2010" both appear in the wild; the volume year is **2010**) |
| Affiliations | ¹ Dept. of Food, Bioprocessing and Nutrition Sciences, **Southeast Dairy Foods Research Center, North Carolina State University**, Raleigh NC · ² Dept. of Food Science, **Northeast Dairy Foods Research Center, Cornell University** · ³ Dept. of Applied Economics and Statistics, **Clemson University** |
| Funding | **Dairy Management, Inc.** and **Danisco Singapore Pte. Ltd.**; paper **FSR 09-40** |
| PDF character | 24 pages, clean digital text layer. **All six tables fully legible; no cell unreadable.** |

**This is the paper Hong & Kim 2020 cites for both its Z statistic and its sniffing
protocol** — the two are methodologically the same lineage (Drake lab, ASTM E679-91,
sure/not-sure modification, delta-method SE, pairwise Z). §5.3 uses that link to settle a
defect in Hong.

**Provenance codes:** **[M]** measured and printed · **[C]** cited · **[Z]** derived by this
wave, never printed.

---

## 1. THE ONE-PARAGRAPH ANSWER

**All 24 best-estimate thresholds (3 compounds × 8 matrices, Table 3, each with an SE), all
24 headspace/matrix partition coefficients on the same 8 matrices (Table 5, each with an
SD), all 18 protein-binding peak areas (Table 2) and all 48 pairwise significance calls
(Tables 4 and 6) are printed, legible and complete.** Because the BETs and the K values were
measured on **the same eight matrices with the same panel and the same emulsions**, this is
the only place in the corpus where the question *"does headspace partition predict the
matrix threshold shift?"* can be answered on paired data. **The answer, computed here for
the first time, is: only for the strongly hydrophobic compound. Across the eight matrices
the correlation between log K_HS/matrix and log BET is r = −0.91 for δ-decalactone
(log P 3.4, correct sign), but +0.60 for diacetyl (log P −2.0) and +0.43 for furaneol
(log P 1.4) — both the WRONG sign** (§6). Used as a predictor of the oil-vs-water threshold
shift, the partition model **under-predicts by 30.6× for diacetyl, over-predicts by 7.2× for
δ-decalactone, and errs by 2.4× for furaneol** (§6.2). The paper asserts the opposite —
*"A high partition coefficient value ideally correlates to a low threshold value ... and we
observed these trends in our study"* — and its own Table 3/Table 5 refute it for two of
three compounds (§9, hazard L-1). One header character, `E⁵`, controls the absolute scale of
every K and is resolved in §4.1.

---

## 2. SYSTEM COMPOSITION AND METHOD — applies to every number below

### 2.1 The eight matrices (Table 1, p. 351)

| matrix as named | composition | pH |
|---|---|---|
| **Water control** | deionized water | NA |
| **0 % control** | **1 % w/v fresh curd calcium caseinate (Cas)** in deionized water | **5.5 and 7.0** |
| **10 % w/o emulsion** | **10 % v/v soybean oil** + **1 % w/v Cas** | **5.5 and 7.0** |
| **20 % w/o emulsion** | **20 % v/v soybean oil** + **1 % w/v Cas** | **5.5 and 7.0** |
| **Oil control** | **soybean oil (Wesson, ConAgra)** | NA |

⚠️ Table 1's heading says *"10 % w/o emulsion"* / *"20 % w/o emulsion"* — **w/o = water-in-oil**
— but the body text and every other table say **oil-in-water**: *"Full-fat cream is an o/w
emulsion..."*, *"10 % v/v and 20 % v/v oil in water emulsions"* (p. 350), *"emulsion
containing 10 % (v/v) fat"* (T3/T5 footnotes). **Table 1's `w/o` is a typo for `o/w`.**

**⚠️ THE CENTRAL CONFOUND, and the paper states it three times.** *"the same amount of
emulsifying agent was used with all the treatments"* — **1 % w/v caseinate in every
matrix.** Therefore *"in 10 and 20 % fat emulsions at pH 7.0, a portion of the 1 % protein
was removed by interacting at the interface and the protein concentration in the aqueous
phase of the emulsions was lower than the initial 1 % protein in the 0 % fat control."*
And *"20 % fat emulsions had more protein adsorbed on the droplet interface compared with
10 % fat emulsions and had a lower residual concentration of emulsifier in the water phase."*
⇒ **Fat content and aqueous-phase protein content move in OPPOSITE directions across the
0/10/20 % series. The design cannot separate a fat effect from a protein effect.** The paper
concedes: *"the precise impact of proteins on the detection threshold or partition
coefficients of diacetyl and furaneol in the 10 and 20 % fat emulsions is unknown even
though proteins were found to have an effect."* **This must ride on every emulsion row.**

### 2.2 Emulsion preparation (p. 350, after Meynier et al. 2005 [C])

Cas 1 % w/v heated to **50 °C in a scraped-surface steam kettle for 30 min**; soybean oil
heated to **50 °C**; coarse dispersion blended **1 min at high speed** (KB7207 Krups) with
hot oil added slowly; fed immediately into a **Panda 2K homogenizer (Niro) at 80 MPa
(72 MPa first stage, 8 MPa second stage)**; **chilled to 5 °C overnight**; pH adjusted at
5 °C with **lactic acid**; volatiles added **as propylene-glycol stock, after the emulsion
was made and cooled to 5 °C**, or to water/oil at **25 °C**, **on the day of the test**.
**All emulsions manufactured in duplicate.**

### 2.3 Measured particle size (p. 351–352, Mastersizer 2000) [M]

| | **d[4,3] volume mean** | **d[0.9]** |
|---|---:|---:|
| **10 % fat** | **0.59 µm** | **1.10 µm** |
| **20 % fat** | **1.58 µm** | **1.64 µm** |

*"pH did not have an effect on the droplet size (d[4.3]) and emulsion distribution (d[0.9])
(results not shown)."* Increasing fat **increased** droplet size and flocculation
(P < 0.05), attributed to the lower free-emulsifier concentration at 20 % fat.
⚠️ Note the direction: **20 % fat has 2.7× the mean droplet diameter of 10 % fat**, so
interfacial area does NOT scale with fat volume — another reason 10 % and 20 % are not a
clean fat series.

### 2.4 Threshold determination (pp. 352–353)

| variable | value as printed |
|---|---|
| Method | **ASTM E679-91 ascending forced-choice method of limits** (ASTM 1991) |
| Series | **serially diluted by a factor of 3**; **seven ascending series tested each time**; **15 mL** per **56 mL plastic soufflé cup**, lidded (Solo) |
| Blanks | **two appropriate blanks per set** (deodorized water, oil, 0 % control or emulsion) |
| **Incubation in the cups** | **2–3 h at room temperature (22 °C)**; *"Preliminary studies confirmed that this was a sufficient time to achieve equilibrium for each compound."* |
| Serving temperature | **22 °C** |
| Panel | **n = 25**; **at least three practice orthonasal threshold tests** before the study |
| Sniffing instruction | open the cup **from the side**, **briefly sniff the headspace without completely removing the lid**; **rest 1 min between each set of three**; **sniff their sleeve** to clear the nasal passage |
| Response | choose the odd sample of three **and give a certainty judgment (sure/not sure)** |
| Individual BET | geometric mean of the last incorrect and first correct concentration; **if "not sure" on a correct choice, that concentration was multiplied by 1.41** [C] Lawless et al. 2000 |
| Group BET | **geometric mean of the individual BETs**; *"the threshold values for each individual were log transformed and then the geometric mean was computed. The antilog then produced the group estimate."* |
| Replication | **"Threshold testing was conducted in duplicate for each compound on different days."** |
| Software | Compusense Five v4.6 |
| SE | **delta method** [C] Sen & Singer 1993 |
| Statistic | **Z = (BET2 − BET1)/(√[SE1² + SE2²])**, standard normal under H₀ |
| Stated rationale for pairwise Z | *"This pairwise comparison approach was adopted to reduce the type II error rate which, if too large in an initial study like this, may result in failing to find differences that can be subsequently investigated more thoroughly."* |

### 2.5 Static headspace partition determination (pp. 354–355)

| variable | value as printed |
|---|---|
| Sample | **10 g of liquid matrix** in a **40 mL amber screw-top vial** (28 × 98 mm, Supelco), PTFE/silicone septum ⇒ **V_HS/V_matrix ≈ 3** |
| **Dosing** | **diacetyl — 10 ppm in ALL matrices**; **δ-decalactone — 100 ppm in water, 0 % fat pH 5.5 and 7; 500 ppm in all emulsions; 5 000 ppm in oil**; **furaneol — 2 000 ppm for water and oil; 3 000 ppm for all other matrices** |
| ⚠️ | **δ-decalactone and furaneol were dosed at DIFFERENT concentrations in different matrices** (5–50× for δ-decalactone). Linearity of the response over that range is **not demonstrated anywhere in the paper.** |
| PG stock concentrations | diacetyl **2 380 ppm**, δ-decalactone **11 600 ppm**, furaneol **12 500 ppm** |
| Equilibration | **sonicated 30 min at 25 °C**, then **incubated 40 °C for 30 min** |
| **Measurement temperature** | **40 °C** ⚠️ — **the thresholds were measured at 22 °C. The two datasets are 18 °C apart.** |
| Injection | **1 mL gastight syringe (Hamilton), splitless**, directly onto **HP 6890 / 5973 inert GC-MS**, ZB-5ms 30 m × 0.25 mm × 0.25 µm |
| SIM ions | **diacetyl m/z 86, δ-decalactone m/z 99, furaneol m/z 128** |
| Ovens | diacetyl 40 °C 1 min → 10 °C/min → 100 → 30 °C/min → 200 (**10.33 min**); δ-decalactone 40 °C 1 min → 10 → 200 → 30 → 250 (**18.67 min**); furaneol 40 °C 1 min → 10 → 100 → 30 → 250 (**12 min**) |
| Injector | **250 °C, 0.049 MPa** |
| ⚠️ **Calibration** | *"A standard curve prepared in **methanol** (δ-decalactone and furaneol) or **diethyl ether** (diacetyl) was generated ... **One microliter of the solution** ... was **directly injected onto the GC**."* **Liquid-phase standards used to quantify a 1 mL gas-phase injection.** Same class of calibration transfer as Meynier 2002 §2 — see `Meynier2002_extraction.md` §7.3 and §7 below. |
| Replication | **all injections in duplicate**, on the duplicate emulsion manufactures |
| Equations | `K_HS/matrix = C_HS / C_matrix` **(1)**, with `C_matrix = C′_matrix − C_HS·(V_HS/V_matrix)` **(2)** [C] Jung & Ebeler 2003 — i.e. a **mass-balance-corrected** liquid concentration |
| Statistics | SAS 9.1, **ANOVA + Tukey HSD, α = 0.05** |

### 2.6 The protein-binding side experiment (pp. 352–353)

Each compound added to **1, 2 and 3 % (w/v) caseinate in deionized water**, **10 mL in
40 mL vials**, same headspace procedure, **SIM peak areas compared** (not converted to
concentrations). Final in-vial concentrations: **diacetyl 10 ppm, δ-decalactone 100 ppm,
furaneol 1 000 ppm.** *"This experiment was replicated two times and all injections were
done in duplicate."* The high furaneol level was needed *"likely because of its thermal
degradation under normal gas chromatography conditions"* [C] Shu et al. 1985; Sanz et al.
1994. **The GC oven program used was the δ-decalactone one (18.67 min) for all three
compounds** — i.e. **a different program from the one used for diacetyl and furaneol in the
partition experiment**, so these peak areas are **not** comparable to Table 5's.

---

## 3. TABLE 3 — THE 24 BEST ESTIMATE THRESHOLDS. **Anchor: Table 3, p. 358.** [M]

Title as printed: *"BEST ESTIMATE THRESHOLDS (BET) OF DIACETYL, δ-DECALACTONE AND FURANEOL
IN DIFFERENT MATRICES AT 22C"*. Footnote: *"BET = Best Estimate Threshold values (PPB = in
parts per billion); S.E., standard errors of BET values."*
**Units as printed: PPB. S.E. is a LINEAR standard error in ppb** (verified in §5.3).

| matrix | **Diacetyl BET (ppb)** | S.E. | **δ-Decalactone BET (ppb)** | S.E. | **Furaneol BET (ppb)** | S.E. |
|---|---:|---:|---:|---:|---:|---:|
| **Water** | **6.0** | 2.21 | **66.0** | 15.2 | **22.3** | 6.99 |
| **0 % fat pH 7.0** | **40.8** | 17.6 | **43.5** | 17.4 | **90.8** | 41.0 |
| **0 % fat pH 5.5** | **44.9** | 20.1 | **35.8** | 12.0 | **46.4** | 15.8 |
| **10 % fat pH 7.0** | **5.6** | 1.66 | **546** | 176 | **56.4** | 25.2 |
| **10 % fat pH 5.5** | **9.2** | 2.62 | **294** | 103 | **67.8** | 23.4 |
| **20 % fat pH 7.0** | **21.8** | 7.43 | **113** | 42.8 | **66.9** | 20.3 |
| **20 % fat pH 5.5** | **8.6** | 2.90 | **329** | 123 | **134** | 60.0 |
| **Oil** | **99.5** | 38.0 | **1,550** | 365 | **27.4** | 9.71 |

**Relative SE [Z]: 23.6 % to 45.2 % (mean 33.6 %) — entirely plausible for a 25-panellist
BET and internally consistent across all 24 cells.** Contrast Hong 2020's ± column, whose
implied linear RSD ranges from 0.0024 % to 120 % (see `hong2020_extraction.md` §4.1).

### 3.1 Matrix/water threshold ratios, computed by this wave [Z]

| matrix | **diacetyl** (log P −2.0) | **δ-decalactone** (log P 3.4) | **furaneol** (log P 1.4) |
|---|---:|---:|---:|
| 0 % fat pH 7.0 | **6.80×** | 0.66× | **4.07×** |
| 0 % fat pH 5.5 | **7.48×** | 0.54× | 2.08× |
| 10 % fat pH 7.0 | **0.93×** | **8.27×** | 2.53× |
| 10 % fat pH 5.5 | 1.53× | **4.45×** | 3.04× |
| 20 % fat pH 7.0 | 3.63× | 1.71× | 3.00× |
| 20 % fat pH 5.5 | 1.43× | **4.98×** | **6.01×** |
| **Oil** | **16.6×** | **23.5×** | **1.23×** |

**⚠️ Three of the twenty-one matrix/water ratios are BELOW 1** — δ-decalactone is **easier**
to smell in 1 % caseinate than in water (0.54–0.66×) and diacetyl is marginally easier in
10 % fat at pH 7 (0.93×). **Inversions are not rare.** Hong 2020 has one in ten; this paper
has three in twenty-one. **Any matrix model that can only increase thresholds is refuted by
both papers.**

**⚠️ The whole ratio range here is 0.54× to 23.5× — a 44× span, and NOT the 29–2 035× of
Hong 2020's soybean paste.** The difference is the matrix: **1 % protein plus ≤20 % oil is a
mild matrix; a whole-soybean paste is not.** These two datasets together bracket the
problem: a lightly loaded model emulsion moves thresholds by singles-to-tens, and a real
food paste by hundreds-to-thousands. **`k2_matrix_and_thresholds.md` §D.4.1's lookup-table
recommendation is reinforced, and the lookup key must include matrix loading, not just
matrix identity.**

---

## 4. TABLE 5 — THE 24 HEADSPACE PARTITION COEFFICIENTS. **Anchor: Table 5, p. 361.** [M]

Title as printed at 400 dpi: *"HEADSPACE PARTITION COEFFICIENT (K**_HS/MATRIX_** **E⁵**)
RESULTS OF DIACETYL, δ-DECALACTONE AND FURANEOL IN DIFFERENT MATRICES AT 40C"*.
Column headers: `K_HS/matrix E⁵` and `Stdev.` **The compound header reads "Furaneol™"**
(trademark symbol; Firmenich mark). **All at 40 °C.**

| matrix | **Diacetyl** | Stdev | **δ-Decalactone** | Stdev | **Furaneol** | Stdev |
|---|---:|---:|---:|---:|---:|---:|
| **Water** | **2.60** | 0.28 | **1.70** | 0.73 | **0.031** | 0.0073 |
| **0 % Cas pH 7.0** | **1.50** | 0.33 | **2.20** | 0.57 | **0.082** | 0.038 |
| **0 % Cas pH 5.5** | **2.70** | 0.45 | **1.60** | 0.28 | **0.022** | 0.0056 |
| **10 % fat pH 7.0** | **1.50** | 0.21 | **0.23** | 0.053 | **0.025** | 0.0076 |
| **10 % fat pH 5.5** | **1.60** | 0.60 | **0.26** | 0.091 | **0.090** | 0.0097 |
| **20 % fat pH 7.0** | **1.70** | 0.21 | **0.15** | 0.028 | **0.029** | 0.011 |
| **20 % fat pH 5.5** | **1.60** | 0.41 | **0.10** | 0.046 | **0.072** | 0.029 |
| **Oil** | **4.80** | 1.47 | **0.01** | 0.0019 | **0.060** | 0.014 |

### 4.1 ⚠️ RESOLVING `E⁵` — the single character that sets the absolute scale

The header is `K_HS/matrix E⁵`. Two readings: the printed number **is** K × 10⁵ (so
K = value × 10⁻⁵), or the printed number **is** K × 10⁻⁵ (so K = value × 10⁵). **The second
is impossible** — a dimensionless air/liquid partition coefficient of 2.6 × 10⁵ would mean
diacetyl is 260 000× more concentrated in air than in water. **Therefore:**

> **K_HS/matrix = (printed value) × 10⁻⁵, dimensionless.**
> Diacetyl in water at 40 °C = **2.60 × 10⁻⁵**. δ-Decalactone in oil = **1.0 × 10⁻⁷**.

**Corroborated by the paper's own cited comparisons** (both of which the paper claims agree
with it — see hazard L-2):

| quantity | this study | cited comparison | source | **ratio [Z]** |
|---|---:|---:|---|---:|
| K, diacetyl, **water** | **2.60 × 10⁻⁵** | **4.5 × 10⁻⁴** | Salvador et al. 1994, via p. 361 | **17.3× higher than this study** |
| K, diacetyl, **oil** | **4.80 × 10⁻⁵** | **6.3 × 10⁻⁴** | " | **13.1×** |
| K, δ-decalactone, **water** | **1.70 × 10⁻⁵** | **1.03 × 10⁻⁴** (pH 5.2) | Guyot et al. 1996, via p. 363 | **6.1×** |

**⇒ The absolute scale of Table 5 is 6–17× below every value it is compared to.** Same
direction, same order of magnitude, and the **same class of cause** as the 6.24× systematic
offset independently derived for Meynier 2002 (`Meynier2002_extraction.md` §7.3): both
papers **calibrate a gas-phase headspace injection against liquid-phase standards**.
**Two independent static-headspace studies, two different labs, two different decades, both
land ~6–17× below the reference values for the same reason.** That is now a *pattern*, not
an anomaly, and it should be recorded as a standing caveat on any static-headspace-derived
partition constant the repo ingests: **absolute values carry a factor-of-10 downward bias
risk; within-study ratios do not.**

### 4.2 K_water/K_matrix suppression factors at 40 °C [Z]

| matrix | **diacetyl** | **δ-decalactone** | **furaneol** |
|---|---:|---:|---:|
| 0 % Cas pH 7.0 | **1.73×** (suppressed) | 0.77× (enhanced) | **0.38×** (enhanced) |
| 0 % Cas pH 5.5 | 0.96× | 1.06× | **1.41×** |
| 10 % fat pH 7.0 | 1.73× | **7.39×** | 1.24× |
| 10 % fat pH 5.5 | 1.63× | **6.54×** | **0.34×** |
| 20 % fat pH 7.0 | 1.53× | **11.3×** | 1.07× |
| 20 % fat pH 5.5 | 1.63× | **17.0×** | **0.43×** |
| **Oil** | **0.54× (ENHANCED — oil makes diacetyl MORE volatile)** | **170×** | **0.52× (enhanced)** |

**Two of three compounds are MORE volatile over oil than over water** — diacetyl by 1.85×
and furaneol by 1.94× — because both are hydrophilic enough that soybean oil is a poor
solvent for them. **Only the log P 3.4 lactone behaves the way a "fat retains aroma" model
expects, and it does so by 170×.** ⚠️ **A repo term of the form "fat suppresses headspace"
is flatly wrong for hydrophilic and amphiphilic odourants and right by two orders of
magnitude for lipophilic ones. There is no middle setting.**

---

## 5. TABLES 4 AND 6 — the 48 pairwise significance calls. **Anchors: Table 4 p. 358, Table 6 p. 362.** [M]

`Sig. = P < 0.05`; `N.S. = P > 0.05`.

### 5.1 Table 4 — effect of FAT (A) and pH (B) on **BET**

| comparison | **Diacetyl** | **δ-Decalactone** | **Furaneol** |
|---|:--:|:--:|:--:|
| **A.** Water vs Oil | **Sig.** | **Sig.** | N.S. |
| Water vs 0 % pH 7 | N.S. | N.S. | N.S. |
| 0 % pH 7 vs 10 % pH 7 | **Sig.** | **Sig.** | N.S. |
| 10 % pH 7 vs 20 % pH 7 | **Sig.** | **Sig.** | N.S. |
| 0 % pH 7 vs 20 % pH 7 | N.S. | N.S. | N.S. |
| Water vs 0 % pH 5.5 | N.S. | N.S. | N.S. |
| 0 % pH 5.5 vs 10 % pH 5.5 | N.S. | **Sig.** | N.S. |
| 10 % pH 5.5 vs 20 % pH 5.5 | N.S. | N.S. | N.S. |
| 0 % pH 5.5 vs 20 % pH 5.5 | N.S. | **Sig.** | N.S. |
| **B.** 0 % pH 7 vs 0 % pH 5.5 | N.S. | N.S. | N.S. |
| 10 % pH 7 vs 10 % pH 5.5 | N.S. | N.S. | N.S. |
| 20 % pH 7 vs 20 % pH 5.5 | N.S. | N.S. | N.S. |

**⇒ pH (7.0 vs 5.5) has NO significant effect on ANY threshold, at ANY fat level. 0/12
significant.** And **furaneol's threshold is unaffected by everything — 0/9 significant,
including water vs pure oil.**

### 5.2 Table 6 — effect of FAT (A) and pH (B) on **partition coefficient**

| comparison | **Diacetyl** | **δ-Decalactone** | **Furaneol** |
|---|:--:|:--:|:--:|
| **A.** Water vs Oil | **Sig.** | **Sig.** | N.S. |
| Water vs 0 % pH 7 | **Sig.** | N.S. | **Sig.** |
| 0 % pH 7 vs 10 % pH 7 | N.S. | **Sig.** | **Sig.** |
| 10 % pH 7 vs 20 % pH 7 | N.S. | N.S. | N.S. |
| 0 % pH 7 vs 20 % pH 7 | N.S. | **Sig.** | **Sig.** |
| Water vs 0 % pH 5.5 | N.S. | N.S. | N.S. |
| 0 % pH 5.5 vs 10 % pH 5.5 | N.S. | **Sig.** | **Sig.** |
| 10 % pH 5.5 vs 20 % pH 5.5 | N.S. | N.S. | N.S. |
| 0 % pH 5.5 vs 20 % pH 5.5 | N.S. | **Sig.** | **Sig.** |
| **B.** 0 % pH 7 vs 0 % pH 5.5 | **Sig.** | N.S. | **Sig.** |
| 10 % pH 7 vs 10 % pH 5.5 | N.S. | N.S. | **Sig.** |
| 20 % pH 7 vs 20 % pH 5.5 | N.S. | N.S. | **Sig.** |

**⇒ pH affects the PARTITION coefficient in 4/12 comparisons (all three furaneol rows and
one diacetyl row) but NEVER the threshold (0/12).** That dissociation is itself a finding:
**a statistically detectable change in headspace concentration produced no detectable change
in perception.** It is the same conclusion as `k2_matrix_and_thresholds.md` §B.5 (Andriot and
Barallat-Pérez), reached here on a third, independent design.

**⚠️ Sixteen of the twenty-four pairwise calls disagree between Table 4 and Table 6** —
i.e. the sensory and instrumental experiments, on identical matrices, reach opposite
significance conclusions two-thirds of the time.

### 5.3 The Z statistic reproduces exactly — which settles a defect in Hong 2020 [Z]

Applying the printed formula `Z = (BET2 − BET1)/√(SE1² + SE2²)` to Table 3:

| comparison | ΔBET | √(SE²+SE²) | **Z [Z]** | |Z| > 1.96? | **Table 4 says** |
|---|---:|---:|---:|:--:|:--:|
| diacetyl, water vs oil | 93.5 | 38.06 | **2.457** | yes | **Sig.** ✅ |
| diacetyl, water vs 0 % pH 7 | 34.8 | 17.74 | **1.962** | marginal | N.S. ⚠️ borderline |
| diacetyl, 0 % pH 7 vs 10 % pH 7 | −35.2 | 17.68 | **−1.991** | yes | **Sig.** ✅ |
| diacetyl, 10 % vs 20 % pH 7 | 16.2 | 7.613 | **2.128** | yes | **Sig.** ✅ |
| diacetyl, 0 % vs 20 % pH 7 | −19.0 | 19.10 | **−0.995** | no | N.S. ✅ |

**The formula reproduces the printed significance calls.** The SEs are therefore **linear
ppb SEs**, and the Z is a genuine standard normal deviate of order 1–3.

**⇒ CONSEQUENCE FOR HONG 2020.** Hong & Kim (2020) cite *this* paper for *this* formula, but
apply it with a **log-scale dispersion in the denominator and a linear difference in the
numerator**, producing "Z-values" up to 24,942 (see `hong2020_extraction.md` §4.2, hazard
H-1). **The defect is now diagnosed, not merely suspected: Hong inherited a correct
statistic from Leksrisompong and broke it by changing the scale of only one of its two
inputs.** This is a clean, citable resolution of an open item in that dossier.

---

## 6. THE PAIRED ANALYSIS — does headspace partition predict the matrix threshold shift? [Z]

**The paper never does this.** It has the data to. Both tables cover **the same eight
matrices**, so `log K_HS/matrix` and `log BET` can be correlated cell by cell.

### 6.1 Correlation across the eight matrices, per compound

Theory: more compound in the headspace ⇒ lower detection threshold ⇒ **r should be NEGATIVE.**

| compound | **log P** | **r(log K, log BET), n = 8** | slope | **sign** |
|---|---:|---:|---:|:--|
| **δ-decalactone** | **+3.4** | **−0.907** | −0.66 | **✅ CORRECT** |
| **diacetyl** | **−2.0** | **+0.600** | +1.55 | **❌ WRONG SIGN** |
| **furaneol** | **+1.4** | **+0.427** | +0.44 | **❌ WRONG SIGN** |

**One of three compounds behaves as partition theory requires. Two behave backwards.**
For diacetyl the wrong sign is driven by the two extreme matrices: **oil has the HIGHEST
partition coefficient (4.80 × 10⁻⁵, 1.85× water) AND the HIGHEST threshold (99.5 ppb,
16.6× water).** More diacetyl in the headspace, and it is *harder* to smell. There is no
partition model that produces that.

### 6.2 The partition model as a threshold predictor — quantified error

Take the water threshold as given and predict each matrix threshold as
`BET_pred = BET_water × (K_water / K_matrix)`. Full 21-cell result:

**Diacetyl** (log P −2.0)
| matrix | obs BET | partition-predicted | **obs/pred** |
|---|---:|---:|---:|
| 0 % pH 7 | 40.8 | 10.4 | **3.92×** |
| 0 % pH 5.5 | 44.9 | 5.8 | **7.77×** |
| 10 % pH 7 | 5.6 | 10.4 | 0.54× |
| 10 % pH 5.5 | 9.2 | 9.8 | 0.94× |
| 20 % pH 7 | 21.8 | 9.2 | 2.38× |
| 20 % pH 5.5 | 8.6 | 9.8 | 0.88× |
| **Oil** | **99.5** | **3.3** | **30.6×** |

**δ-Decalactone** (log P +3.4)
| matrix | obs BET | partition-predicted | **obs/pred** |
|---|---:|---:|---:|
| 0 % pH 7 | 43.5 | 51.0 | 0.85× |
| 0 % pH 5.5 | 35.8 | 70.1 | 0.51× |
| 10 % pH 7 | 546 | 487.8 | **1.12×** |
| 10 % pH 5.5 | 294 | 431.5 | 0.68× |
| 20 % pH 7 | 113 | 748.0 | **0.15×** |
| 20 % pH 5.5 | 329 | 1 122.0 | 0.29× |
| **Oil** | **1 550** | **11 220** | **0.14× (7.2× over-prediction)** |

**Furaneol** (log P +1.4)
| matrix | obs BET | partition-predicted | **obs/pred** |
|---|---:|---:|---:|
| 0 % pH 7 | 90.8 | 8.4 | **10.8×** |
| 0 % pH 5.5 | 46.4 | 31.4 | 1.48× |
| 10 % pH 7 | 56.4 | 27.7 | 2.04× |
| 10 % pH 5.5 | 67.8 | 7.7 | **8.83×** |
| 20 % pH 7 | 66.9 | 23.8 | 2.81× |
| 20 % pH 5.5 | 134 | 9.6 | **14.0×** |
| **Oil** | **27.4** | **11.5** | **2.38×** |

**Summary — the headline number of this dossier:**

| compound | **oil/water partition-model error** | direction |
|---|---:|---|
| **diacetyl** | **30.6×** | model **under**-predicts the threshold shift |
| **δ-decalactone** | **7.2×** | model **over**-predicts |
| **furaneol** | **2.4×** | model under-predicts |
| **across all 21 matrix cells** | **0.14× to 30.6× — a 219× span, in both directions** | — |

**⇒ On the only same-sample paired dataset in the corpus, a partition-based prediction of
matrix threshold shift errs by up to 30× and gets the direction wrong for two of three
compounds.** This is the strongest available quantitative extension of
`k2_matrix_and_thresholds.md` §B.4/§B.5 and §D.4.5 ("stop deriving matrix odour activity
from `f_free`"): **the earlier conclusion was reached by comparing a binding ceiling to
cross-study threshold shifts; here it is reached directly, on matched samples, with the
partition term measured rather than modelled. The conclusion survives and hardens.**

⚠️ **One legitimate defence of the partition model, recorded for fairness:** the K values
were measured at **40 °C** and the thresholds at **22 °C** (§2.4, §2.5). An 18 °C gap
changes K by roughly 2–3× for these compounds (cf. `Meynier2002_extraction.md` §7.2). **That
can explain a factor of 2–3. It cannot explain 30.6×, and it cannot explain a sign
reversal**, because the temperature correction is monotone and compound-specific in
magnitude but never in sign.

---

## 7. TABLE 2 — the protein-binding side experiment. **Anchor: Table 2, p. 357.** [M]

Title as printed: *"PEAK AREAS (AREA COUNTS) OF DIACETYL, δ-DECALACTONE AND FURANEOL IN 1,
2, AND 3 % PROTEIN (W/V) SOLUTIONS AT PH 7.0 AND PH 5.5 AT 40C"*. Footnote: *"a,b means in a
row not followed by a common letter are statistically different (P < 0.05)."*

| pH | compound | **1 % Cas** | **2 % Cas** | **3 % Cas** | **3 %/1 % [Z]** |
|---:|---|---:|---:|---:|---:|
| **7.0** | **Diacetyl** | **174,000ᵃ** | 138,000ᵃ·ᵇ | **137,000ᵇ** | **0.787× — headspace DECREASES** |
| 7.0 | δ-Decalactone | 12,800ᵃ | 13,600ᵃ | 14,200ᵃ | 1.109× (N.S.) |
| **7.0** | **Furaneol** | **21,600ᵇ** | 26,500ᵃ | **36,600ᵃ** | **1.694× — headspace INCREASES** |
| 5.5 | Diacetyl | 384,000ᵃ | 409,000ᵃ | 405,000ᵃ | 1.055× (N.S.) |
| 5.5 | δ-Decalactone | 21,800ᵃ | 18,300ᵃ | 11,600ᵃ | 0.532× (N.S. despite a 1.9× fall) |
| 5.5 | Furaneol | 59,900ᵃ | 95,600ᵃ | 46,300ᵃ | 0.773× (N.S., **non-monotone**) |

**Findings as the paper states them (p. 352):** *"Static headspace analysis confirmed that
diacetyl headspace concentration decreased as protein in solution increased from 1 to 3 % at
pH 7.0. The opposite prevailed for furaneol. Proteins did not interact with δ-decalactone at
either pH (P > 0.05) nor had any interaction at pH 5.5 for all compounds (P > 0.05)."*
Mechanism offered: *"At acidified pH, proteins are more aggregated compared with neutral pH,
and this could lead to less binding sites available ... Thus, at neutral or alkali pH,
casein binds more carbonyl, alcohol and ester volatile compounds than at acidic pH."*
For diacetyl specifically, quoting **Fares et al. 1998** [C]: *"diacetyl reacted with the
primary amino groups of protein, the terminal amino groups and the ε-lysyl residues ...
Strong bonding occurred via **covalent bonding** between amino groups and the carbonyl group
on the diacetyl whereas weak bonding occurred via hydrogen bonding."*

**⚠️ THREE THINGS THIS TABLE SAYS THAT THE REPO SHOULD RECORD.**
1. **The maximum measured protein effect on headspace, over a 3× protein range at 10–30 g/L,
   is 1.27× (diacetyl, suppression) and 1.69× (furaneol, ENHANCEMENT).** A fourth
   independent confirmation of `k2_matrix_and_thresholds.md` §B.4's single-digit ceiling —
   after Andriot (1.25–3.7× at 40 g/L), Meynier (1.07–1.39× at ~34 g/L, this wave) and
   Barallat-Pérez.
2. **Protein can push a volatile OUT of solution as well as hold it in.** Furaneol's
   headspace rises 1.69× with added caseinate. `f_free` formalisms that can only reduce free
   concentration cannot represent this. *"proteins drove furaneol to the headspace"* (p. 364).
3. **The pH dependence has the sign the covalent-adduct literature predicts** — binding at
   pH 7, none at pH 5.5 — which is the same direction as
   `anantharamkrishnan2020_extraction.md` (this wave): **carbonyl–lysine adduct formation
   requires the free base amine and is suppressed at acid pH.** Two independent papers in
   this wave, on two different proteins and two different carbonyls, agree.

⚠️ **Peak areas are not concentrations and were acquired under a different oven program from
Table 5** (§2.6). **Do not convert them to partition coefficients or binding constants.**
The 2.2× and 2.8× jumps in *absolute* area between pH 7 and pH 5.5 for diacetyl and furaneol
are un-normalised instrument response, not a pH effect on volatility.

---

## 8. THE CITED COMPARISON VALUES — all [C]

| quantity | value | source as cited |
|---|---:|---|
| diacetyl orthonasal threshold, **water** | **5 ppb** | Druaux & Voilley 1997, via p. 357 |
| diacetyl orthonasal threshold, **oil** | **50 ppb** | " |
| diacetyl recognition threshold, water | **5 ppb**, skewed distribution; a **hypersensitive group perceived below 0.625 ppb** | Lawless et al. 1994 |
| diacetyl orthonasal, water vs **sunflower oil** | **no difference reported** | van Gemert 2003 |
| **δ-decalactone**, water | **100 ppb** | van Gemert 2003 |
| δ-decalactone **taste** thresholds, water vs **butter oil** | **140 ppb vs 1,400 ppb (10×)** | Siek et al. 1969 |
| **furaneol**, water | **1–10 ppb** | Larsen & Poll 1992 |
| furaneol, **sunflower oil** | **25–50 ppb** | Preininger & Grosch 1994 |
| furaneol, published range | **1–10 ppb to 1,700 ppb** | van Gemert 2003 |
| furaneol threshold **lower at pH 4.5 than pH 7** | direction only | Buttery, Takeoka & Ling 1995 |
| δ-decalactone K in water at **pH 5.2** | **1.03 × 10⁻⁴** | Guyot et al. 1996 |
| diacetyl K, **sunflower oil / water** | **6.3 × 10⁻⁴ / 4.5 × 10⁻⁴** | Salvador et al. 1994 |
| pKa of the δ-hydroxy decanoic acid form | **4.46** | Kearney et al. 1993 |
| max real-time diacetyl exposure, butter-popcorn workers | **525 ppm** | Martyny et al. 2008 |
| log P values used | diacetyl **−2.0**, δ-decalactone **3.4** (both Guyot et al. 1996); furaneol **1.4** (Relkin et al. 2004) | — |

**⚠️ The furaneol water threshold literature spans 1 to 1 700 ppb — 1 700×** — and this
study's 22.3 ppb sits in the middle. Same class of aqueous-threshold irreproducibility as
`hong2020_extraction.md` §5 (butyric acid, 5 540×). **Two papers in this wave independently
establish that cross-study aqueous orthonasal thresholds are reproducible to no better than
~3 orders of magnitude.**

**⚠️ A pH prediction the paper made and could not test.** *"Kearney et al. (1993) reported
that the pKa of the carboxyl group on the δ-hydroxy decanoic acid form was 4.46. Below that
pKa, a large proportion of the compound is protonated and forms a lactone ring whereas above
the pKa, a higher amount of the compound will be in the hydroxyacid form. **If thresholds
were determined below the pKa, 4.46, differences in BET values may be observed.**"*
The tested range (5.5–7.0) is entirely above the pKa. **This is a specific, falsifiable,
untested prediction about lactone behaviour below pH 4.5** — worth recording as an open
question for any acidic matrix the repo models.

---

## 9. NAMED LAUNDERING HAZARDS

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| **L-1** | **"A high partition coefficient value ideally correlates to a low threshold value (as more compounds are presented in the headspace to be detected by human nose), and we observed these trends in our study."** | **Observed for ONE of three compounds.** Over the eight matrices, r(log K, log BET) = **−0.91** for δ-decalactone but **+0.60** for diacetyl and **+0.43** for furaneol — both the wrong sign. Diacetyl in oil has the **highest K and the highest threshold simultaneously.** The next sentence does narrow the claim to δ-decalactone, but the general statement is what a citing paper will lift. | p. 364–365 vs T3/T5 |
| **L-2** | **"Our results agreed with previous published values: Salvador et al. (1994) also found diacetyl in sunflower oil to have higher K_HS/matrix (6.3E-4) compared with water (4.5E-4)."** | Only the **ordering** agrees. The **magnitudes differ by 13–17×** (this study: 4.80 × 10⁻⁵ and 2.60 × 10⁻⁵). "Agreed with previous published values" is not a defensible description of a 17× disagreement. | p. 361 |
| **L-3** | **Table 5's `E⁵` header** | Must be read as **K = printed × 10⁻⁵**. A reader who takes `E⁵` at face value gets partition coefficients of order 10⁵ — wrong by **10 orders of magnitude**. | T5 header, p. 361 |
| **L-4** | **Table 1's "10 % w/o emulsion" / "20 % w/o emulsion"** | **w/o = water-in-oil.** Every other statement in the paper says **oil-in-water**. Typo, but it inverts the phase topology of the matrix. | T1, p. 351 |
| **L-5** | The **0/10/20 % fat series** presented as a fat effect | **Aqueous protein concentration falls as fat rises** (constant 1 % total emulsifier), and **droplet size rises 2.7× from 10 % to 20 %**. Three variables move together; the paper concedes the effect "is unknown". Do not ingest any of these rows as a pure fat dependence. | pp. 352, 358 |
| **L-6** | **Table 2 peak areas** across pH | Acquired **under a different GC oven program** from Table 5 and **not normalised to any internal standard**. The 2.2–2.8× area jump between pH 7 and pH 5.5 is instrument response, not chemistry. Only the **within-row** 1→3 % trends are interpretable. | pp. 352–353, T2 |
| **L-7** | Thresholds at **22 °C** and partition coefficients at **40 °C** presented as a matched pair throughout | An 18 °C gap. Worth ~2–3× on K. The paper never mentions the mismatch. It does not explain the §6 failures but it must be carried on any joint use of the two tables. | §2.4 vs §2.5 |
| **L-8** | **δ-decalactone and furaneol dosed at 5–50× different concentrations in different matrices** for the partition measurement | Linearity over that range is **never demonstrated**. Table 5's cross-matrix comparisons for those two compounds assume it. | p. 354 |

---

## 10. CONSOLIDATED NEW-PARAMETER TABLE

**Common conditions:** BETs — orthonasal 3-AFC ASTM E679-91, factor-3 series, 7 series,
n = 25 panellists, **22 °C**, 15 mL in 56 mL lidded cups, 2–3 h equilibration, duplicate on
different days, group BET = geometric mean, delta-method SE.
**Partition — static headspace GC-MS SIM, 10 g in a 40 mL vial (V_HS/V_matrix ≈ 3), 30 min
sonication at 25 °C then 30 min at 40 °C, mass-balance-corrected (eq. 2), duplicate
injections on duplicate emulsions. K reported as printed × 10⁻⁵, dimensionless.**

| # | parameter | value | units | matrix / conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1–8 | **BET diacetyl**, 8 matrices | **6.0 / 40.8 / 44.9 / 5.6 / 9.2 / 21.8 / 8.6 / 99.5** (SE 2.21 / 17.6 / 20.1 / 1.66 / 2.62 / 7.43 / 2.90 / 38.0) | ppb | water / 0 % pH7 / 0 % pH5.5 / 10 % pH7 / 10 % pH5.5 / 20 % pH7 / 20 % pH5.5 / oil, all 22 °C | M | T3 p.358 |
| 9–16 | **BET δ-decalactone**, same 8 | **66.0 / 43.5 / 35.8 / 546 / 294 / 113 / 329 / 1,550** (SE 15.2 / 17.4 / 12.0 / 176 / 103 / 42.8 / 123 / 365) | ppb | " | M | T3 |
| 17–24 | **BET furaneol**, same 8 | **22.3 / 90.8 / 46.4 / 56.4 / 67.8 / 66.9 / 134 / 27.4** (SE 6.99 / 41.0 / 15.8 / 25.2 / 23.4 / 20.3 / 60.0 / 9.71) | ppb | " | M | T3 |
| 25–32 | **K_HS/matrix diacetyl**, same 8 | **2.60 / 1.50 / 2.70 / 1.50 / 1.60 / 1.70 / 1.60 / 4.80** ×10⁻⁵ (SD 0.28 / 0.33 / 0.45 / 0.21 / 0.60 / 0.21 / 0.41 / 1.47 ×10⁻⁵) | dimensionless | same 8 matrices, **40 °C** | M | T5 p.361 |
| 33–40 | **K_HS/matrix δ-decalactone**, same 8 | **1.70 / 2.20 / 1.60 / 0.23 / 0.26 / 0.15 / 0.10 / 0.01** ×10⁻⁵ (SD 0.73 / 0.57 / 0.28 / 0.053 / 0.091 / 0.028 / 0.046 / 0.0019 ×10⁻⁵) | dimensionless | " | M | T5 |
| 41–48 | **K_HS/matrix furaneol**, same 8 | **0.031 / 0.082 / 0.022 / 0.025 / 0.090 / 0.029 / 0.072 / 0.060** ×10⁻⁵ (SD 0.0073 / 0.038 / 0.0056 / 0.0076 / 0.0097 / 0.011 / 0.029 / 0.014 ×10⁻⁵) | dimensionless | " | M | T5 |
| 49 | **oil/water BET ratio**, diacetyl / δ-decalactone / furaneol | **16.6 / 23.5 / 1.23** | × | 22 °C | Z | §3.1 |
| 50 | **K_water/K_oil**, same order | **0.54 / 170 / 0.52** | × | 40 °C | Z | §4.2 |
| 51 | **partition-model error on the oil/water threshold shift** | **30.6× (diacetyl, under) / 7.2× (δ-decalactone, over) / 2.4× (furaneol, under)** | × | paired, same samples | Z | §6.2 |
| 52 | **r(log K_HS/matrix, log BET) over 8 matrices** | **+0.600 (diacetyl) / −0.907 (δ-decalactone) / +0.427 (furaneol)** | r | " | Z | §6.1 |
| 53 | **full 21-cell partition-prediction error range** | **0.14× to 30.6×** | × | " | Z | §6.2 |
| 54 | **max caseinate effect on headspace, 10→30 g/L** | **0.787× (diacetyl, pH 7) to 1.694× (furaneol, pH 7)** | × | 40 °C, pH 7.0 | Z from M | §7 |
| 55 | caseinate effect at pH 5.5 | **all N.S.** (0.53–1.06×) | × | 40 °C, pH 5.5 | M | T2 |
| 56 | emulsion droplet size | **d[4,3] 0.59 µm @10 % fat; 1.58 µm @20 % fat**; d[0.9] 1.10 / 1.64 µm | µm | 80 MPa homogenisation, pH-independent | M | p.351 |
| 57 | **absolute-scale offset vs cited K values** | **6.1× to 17.3× low** | × | vs Salvador 1994 and Guyot 1996 | Z | §4.1 |
| 58 | pH effect on BET | **0 of 12 pairwise comparisons significant** | — | 5.5 vs 7.0, all fat levels | M | T4B |
| 59 | pH effect on K | **4 of 12 significant** (all 3 furaneol + 1 diacetyl) | — | " | M | T6B |
| 60 | log P values used | diacetyl **−2.0**, δ-decalactone **+3.4**, furaneol **+1.4** | — | cited | C | pp.350, 359, 364 |

---

## 11. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Leksrisompong 2010 is a **new source**; a declaration amendment must
> be approved before any wave fits any row. This dossier does not edit the declaration.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **Table 3, 24 BETs** | 24 | **HOLD-OUT** | Same-panel matrix thresholds with real linear SEs. Together with Hong 2020 these are the only two same-panel matrix threshold sets in the corpus. Both should stay clean. |
| **Table 5, 24 K values — RATIOS** | 21 | **FIT-ELIGIBLE (ratios only)** | The 6–17× absolute-scale offset (§4.1) cancels in a within-study ratio. `K_water/K_matrix` is defensible; raw K is not. Same ruling as `Meynier2002_extraction.md` §12. |
| **Table 5, 24 K values — ABSOLUTES** | 24 | **REJECT for calibration; INGEST only with `absolute_scale_suspect: true, offset_vs_literature: 6-17x_low`** | Two independent static-headspace papers in this wave land 6–17× low for the same methodological reason. |
| **§6.1 correlations and §6.2 error table** | 3 + 21 | **HOLD-OUT — the wave's primary falsification test for `f_free`/partition→perception** | Proposed guard: any repo path that derives a matrix odour threshold from a partition or binding term must be checked against **r = +0.60 for diacetyl**. A model that cannot produce a *positive* correlation for a hydrophilic carbonyl is refuted by this dataset. **Do not fit; use to reject.** |
| **§5.2 pH dissociation (4/12 K vs 0/12 BET)** | 1 | **HOLD-OUT — diagnostic** | Independent third confirmation of `k2_matrix_and_thresholds.md` §B.5. |
| **Table 2, 18 peak areas** | 18 | **DIRECTION ONLY; NOT FITTABLE** | Un-normalised areas, different oven program (L-6). The usable content is the three within-row trends and the pH-7-vs-5.5 contrast. |
| **Table 4 / Table 6, 48 significance calls** | 48 | **METADATA** | Useful as acceptance gates ("the model must not predict a significant pH effect on any BET"), not as parameters. |
| **§2.3 droplet sizes** | 4 | **FIT-ELIGIBLE as matrix metadata** | Clean, measured, replicated. Needed if the repo ever models interfacial area. |
| **All 12 emulsion rows (10 % and 20 % fat)** | 12 | **QUARANTINE from any pure-fat-dependence fit** | Fat, aqueous protein and droplet size are confounded by design (L-5), stated by the authors. |
| **§8 cited threshold and K values** | 15 | **REJECT as parameters; RETAIN the furaneol 1–1,700 ppb spread as a calibration fact** | Second-hand throughout. |

---

## 12. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Guyot, C., Bonnafont, C., Lesschaeve, I., Issanchou, S., Voilley, A. & Spinnler, H. E.
   (1996)** — the source of **both** the log P values and a **δ-decalactone partition
   coefficient in water at pH 5.2 (1.03 × 10⁻⁴)** that is 6.1× above this paper's. It also
   reportedly contains **partition coefficients correlated with sensory odour intensity**
   (cited on p. 349), i.e. **a second paired instrumental/sensory dataset.** Highest-value
   retrieval from this paper.
2. **Salvador, D., Bakker, J., Langley, G. J., Potjewijd, R., Martin, A. & Elmore, J. S.
   (1994)** — the source of the diacetyl **4.5 × 10⁻⁴ (water) / 6.3 × 10⁻⁴ (sunflower oil)**
   partition pair. Needed to adjudicate the 13–17× absolute-scale question in §4.1, which
   now bears on two papers in this wave.
3. **Fares, K., Landy, P., Guilard, R. & Voilley, A. (1998)** — the **covalent** diacetyl /
   sodium-caseinate binding study: two site classes, covalent to amino groups plus hydrogen
   bonding, and the glycosylated-caseinate control that abolished binding. **This is a
   quantitative covalent flavour–protein binding paper on a food protein** — precisely the
   class of data `anantharamkrishnan2020_extraction.md` (this wave) turns out **not** to
   supply.
4. **Jung, D.-M. & Ebeler, S. E. (2003)** — source of eq. (2), the mass-balance correction.
   Worth having so the repo's own headspace mass balance can be traced to a primary source.
5. **Adhikari, K., Hein, K. A., Elmore, J. R., Heymann, H. & Willott, A. M. (2006)**,
   *J. Sensory Studies* **21**, 626–643 — **retronasal thresholds of diacetyl, HEXANAL and
   δ-decalactone in skim milk at 7 °C, with measured compound–compound interactions**
   (cited p. 349: *"The presence of one other compound either suppressed or increased the
   retronasal threshold value of the other compound"*). **A hexanal threshold in a real dairy
   matrix, plus the only mixture-interaction threshold data the corpus has any pointer to.**

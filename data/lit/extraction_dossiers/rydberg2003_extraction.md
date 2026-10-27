# Rydberg et al. 2003 — EXTRACTION
### *"Investigations of Factors That Influence the Acrylamide Content of Heated Foodstuffs"*
### Homogenised potato, oven-heated 25 min. Net acrylamide vs initial (cold) pH peaks near **pH 7.7 (180 °C) and 8.0 (160 °C)**. It is ~40–50 % lower at the native pH 5.7–5.8 and falls again by pH 9.7. Acidification with HCl cuts formation and speeds the loss of spiked acrylamide (47 % lost at pH 3 vs 22 % at pH 5.7, 180 °C). Also: time course at 200 °C, temperature series, amino-acid, sugar, antioxidant, citrate and fish additions. All single samples. Read 2026-10-27.

**Full read of `data/articles/rydberg2003.pdf`** (7 pages, born-digital text plus five 600-dpi
bilevel figure images). Tables were taken from the text layer. **Figures 1–5 were digitised from
the embedded images**: Fig. 4 by automated marker detection (blob centroids, calibrated on the
axis ticks), Figs. 1, 2, 3, 5 by reading pixel positions against the tick marks. The Fig. 4
digitisation was validated by refitting the authors' quadratic: R² 0.769 reproduced exactly at
180 °C, and 0.763 vs printed 0.755 at 160 °C.

**Provenance:** **[M]** printed number · **[F]** read off a figure by me (read error stated) ·
**[Q]** qualitative statement · **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Per Rydberg, Sune Eriksson, Eden Tareke, Patrik Karlsson, Lars Ehrenberg, Margareta Törnqvist* |
| Title | "Investigations of Factors That Influence the Acrylamide Content of Heated Foodstuffs" |
| Venue | *J. Agric. Food Chem.* **2003**, 51 (24), **7012–7018** |
| DOI | 10.1021/jf034649+ |
| Affiliation | Dept. of Environmental Chemistry, Stockholm University; AnalyCen Nordic AB, Lidköping |
| Dates | received 18 June 2003, revised 12 Sept 2003, accepted 15 Sept 2003, web 18 Oct 2003 |

---

## 1. METHODS [M]

| item | as printed |
|---|---|
| Matrix | Peeled Swedish potatoes, ~500 g from 5–6 tubers per sample, chopped; 400 g homogenised (Waring blender, 3 + 3 min, warms 22 → ~40 °C); additions mixed in. Cultivars: **Frieslander batch 1** (amino acids, sugars), **Frieslander batch 2** (HCl), **Ukama** (antioxidants, peroxides, trisodium citrate), **Rosor** (citric acid) |
| Sample form | homogenate poured into ice-cube trays, ~15 g cubes frozen at −20 °C; **3 cubes (~45 g) per sample**, heated from frozen |
| Oven program (default) | GC oven (HP 5790A): 70 °C 0.1 min → **40 °C/min** → **180 °C, held 25 min**. Other temperatures (120–180 °C) for Figs. 4–5 |
| Time course (Fig. 1) | 30 g commercial frozen pre-fried potato strips, 70 °C → 40 °C/min → **200 °C**, held 1 s, 4, 8, 12, 16, 20, 24 or 60 min |
| pH | Orion 290A meter, Ross 8102SC electrode, calibrated at 22 °C. Measured on the **upper layer of the homogenate incl. foam, before freezing and heating**. The pH of these samples "was measured at ambient temperature prior to heating" [!]. No hot or post-heating pH |
| pH adjustment | raised with **NaOH** (Fig. 4); lowered with **HCl** (Table 2, Fig. 5), citric acid (Table 3) or ascorbic acid |
| Endogenous composition (Frieslander batch 1, wet basis) | free **Asn 17.0 mM**, **glucose 21 mM**, **fructose 2.8 mM** |
| Reporting basis | µg/kg of **original (unheated) sample weight**, i.e. corrected for weight loss (50–70 %); "as consumed" is 2–3× higher |
| Analysis | LC-MS/MS, ¹³C3-acrylamide internal standard; LOD ≈ 10 µg/kg (50 µg/kg in the microwave test); recovery ~100 ± 7.5 %; CV ≈ 5 % |
| Replication | **single samples** for every addition. Control CV 11.9 % (Table 1) and 10.1 % (Table 4); contrast error 13.7 % (Table 1, n = 3 controls), 11.3 % (Table 4, n = 4), pooled 12.5 %. Authors: contrasts > ~25 % unlikely to be measurement error |

---

## 2. pH DEPENDENCE

### 2.1 Fig. 4 — NaOH series, 180 and 160 °C, 25 min [F]

pH raised by NaOH; cultivar/batch not stated for this figure. [Z] The three 180 °C points at pH
5.73–5.81 (541, 430, 499) coincide with Table 1's controls I–III (542 @ 5.74, 427 @ 5.82, 497 @
5.82), so the 180 °C series is probably Frieslander batch 1. Read error: **±0.02 pH, ±5 µg/kg**.

| 180 °C: pH | acrylamide (µg/kg) | | 160 °C: pH | acrylamide (µg/kg) |
|---|---|---|---|---|
| 5.73 | 541 | | 5.82 | 318 |
| 5.81 | 430 | | 6.33 | 294 |
| 5.81 | 499 | | 6.89 | 449 |
| 6.00 | 639 | | 8.10 | 549 |
| 6.34 | 566 | | 8.58 | 589 |
| 6.34 | 649 | | 9.16 | 412 |
| 6.90 | 636 | | 9.66 | 376 |
| 6.90 | 828 | | | |
| 8.11 | 844 | | | |
| 8.58 | 738 | | | |
| 9.17 | 597 | | | |
| 9.67 | 460 | | | |

**Printed quadratic fits [M]:** 180 °C: y = −86.5x² + 1330x − 4300, R² = 0.769; 160 °C:
y = −58.7x² + 934x − 3230, R² = 0.755 (y in µg/kg, x = pH).

**The maximum [Z]:**

| | vertex pH | value at vertex | value at pH 5.7–5.8 | ratio max / native |
|---|---|---|---|---|
| 180 °C, printed fit | **7.69** | 812 | 471 (pH 5.7) | 1.7 |
| 180 °C, my refit of [F] points | 7.69 | 800 | — | — |
| 160 °C, printed fit | 7.96 | 485 | 213 (pH 5.8) | 2.3 |
| 160 °C, my refit of [F] points | **8.01** | 537 | 253 (pH 5.8) | 2.1 |

**[!] The printed 160 °C coefficients are rounded too coarsely to reproduce the drawn curve:**
they give 485 at the vertex, while the drawn curve and the refit peak at ~537. Use the refit or the
points, not the printed 160 °C equation. The abstract and text say the maximum is "around pH 8".
The data place it at **7.7 ± ~0.3 (180 °C) and 8.0 ± ~0.3 (160 °C)** [Z; the curvature is set by
only 3–4 points above pH 8]. The fall from the maximum to pH 9.7 is to ~57 % (180 °C) and ~70 %
(160 °C) of the peak [Z, from the points].

### 2.2 Table 2 and Fig. 5 — HCl series, Frieslander batch 2, 25 min

**Table 2 [M]** (180 °C, 25 min, single samples, corrected for weight loss; spike 1000 µg/kg
acrylamide added before heating):

| measured pH | acrylamide (µg/kg) | acrylamide after prior addition of 1000 µg/kg (µg/kg) | loss of added acrylamide (%) |
|---|---|---|---|
| 3.0 (HCl) | 229 | 759 | 47 |
| 4.8 (HCl) | 356 | 1050 | 31 |
| 5.7 (unadjusted) | 600 | 1380 | 22 |

[Z] check: loss = (no-spike + 1000 − spiked)/1000 → 47.0, 30.6, 22.0 %. Matches.

**Fig. 5 [F]** (pH of original 5.72, acidified with HCl, 25 min; read error ±5 µg/kg; the 120 °C
points are indistinguishable from 0):

| T (°C) | pH 5.72 | pH 4.78 | pH 2.96 |
|---|---|---|---|
| 120 | ~3 | ~0 | ~0 |
| 140 | 21 | 13 | 6 |
| 160 | 154 | 99 | 18 |
| 180 | 601 | 356 | 232 |

The 180 °C column equals Table 2 (600 / 356 / 229), so Fig. 5 and Table 2 are the same samples.
**[!] Batch effect:** at 160 °C, native pH, Fig. 5 (batch 2) reads 154 µg/kg, while Fig. 4 reads
318 at pH 5.82 (probably batch 1). Same nominal conditions, a factor of ~2 between batches.

### 2.3 Table 3 — citric acid, cv. Rosor, 180 °C, 25 min [M]

| addition | amount added (mmol/kg; % wet weight) | measured pH | acrylamide (µg/kg) |
|---|---|---|---|
| none (controls) | — | 5.75 | 275 |
| citric acid | 6.5 (0.125) | 5.61 | 225 |
| citric acid | 26 (0.500) | 5.24 | 180 |
| citric acid | 104 (2.00) | 4.54 | 149 |

Text [Q]: the effects of citric acid resembled those of HCl acidification.

### 2.4 Ascorbic acid, microwave [M]

30 g potato with 0.30, 1.7 or 8.3 wt % ascorbic acid, microwave 3 min at 750 W: acrylamide
**3100 → 420 → < 50 µg/kg** (LOD of that test), at pH **5.19, 4.10, 3.35**. Authors [Q]: the drop
is not explained by pH alone; they propose pH plus water binding.

### 2.5 Apparent pH slopes of NET content [Z]

These are slopes of the final content, d ln C / d pH, not of rate constants. Formation and
elimination are confounded, and the pH is the cold initial pH.

| source | T | pH span | d ln C / d pH |
|---|---|---|---|
| Fig. 4 printed fit | 180 °C | 5.7 → 7.0 | 0.38 |
| Fig. 4 refit | 160 °C | 5.8 → 7.0 | 0.53 |
| Table 2 / Fig. 5 | 180 °C | 2.96 → 5.72 | 0.35 |
| Table 2 / Fig. 5 | 180 °C | 4.78 → 5.72 | 0.56 |
| Fig. 5 | 160 °C | 2.96 → 5.72 | 0.78 (pH 2.96 point is 18 ± 5, so ±0.1) |
| Fig. 5 | 160 °C | 4.78 → 5.72 | 0.47 |
| Table 3 (citric) | 180 °C | 4.54 → 5.75 | 0.51 |

---

## 3. OTHER FACTOR STUDIES

### 3.1 Time at 200 °C — Fig. 1 [F] (commercial potato strips, 30 g)

Squares = corrected for weight loss (per original weight); circles = uncorrected; triangles =
weight loss (right axis). Read error **±70 µg/kg, ±1.5 % weight loss**. The x-axis is categorical
(the 24 → 60 min gap is broken).

| time at 200 °C | corrected (µg/kg) | uncorrected (µg/kg) | weight loss (%) |
|---|---|---|---|
| 1 s | ~0 | ~0 | 4 |
| 4 min | ~60 | ~120 | 26 |
| 8 min | 375 | 710 | 41 |
| 12 min | 830 | 1560 | 46 |
| 16 min | 1460 | 3770 | 61 |
| 20 min | **2630** | **7700** | 66 |
| 24 min | 2190 | 5620 | 62 |
| 60 min | 620 | 1950 | 70 |

[Z] check: corrected ≈ uncorrected × (1 − loss) holds within read error (20 min: 7700 × 0.34 =
2620). Maximum near 20 min, then net loss: −17 % by 24 min and −76 % by 60 min. Text [Q]: content
first increases "exponentially", then decreases as degradation predominates.

### 3.2 Temperature, 15 min — Fig. 2 [F] (French fries; data from Tareke et al. 2002, ref 2)

Read error ±30 µg/kg, ±1.5 %.

| oven T (°C) | corrected (µg/kg) | uncorrected (µg/kg) | weight loss (%) |
|---|---|---|---|
| control (unheated) | ~160 | ~160 | ~1 |
| 100 | ~160 | ~190 | 15 |
| 120 | ~205 | ~250 | 21 |
| 140 | ~290 | ~420 | 30 |
| 160 | ~520 | ~840 | 42 |
| 180 | ~1050 | ~2010 | 50 |
| 200 | ~1650 | ~3530 | 55 |
| 220 | ~2320 | ~5130 | ~55 |

**[!]** The commercial strips are pre-fried. The unheated control at ~160 µg/kg is the factory
load, not a blank.

### 3.3 Fish (cod) admixture — Fig. 3 [F]

Grated potato + minced cod, 0–60 g fish in a 60 g patty (0–100 %). Frying pan (conditions per ref 2,
not restated). Inset: oven (default 180 °C, 25 min program). Read error ±2 µg/kg (pan), ±15 µg/kg
(oven inset).

| fish (%) | pan, uncorrected | pan, corrected | oven, uncorrected | oven, corrected |
|---|---|---|---|---|
| 0 | 188 | 146 | 757 | 516 |
| 16.7 | 126 | 100 | — | — |
| 33.3 | 44 | 36 | — | — |
| 50 | 27 | 22 | 101 | 75 |
| 66.7 | 15 | 12 | — | — |
| 83.3 | 4 | 3 | — | — |
| 100 | ~2 | ~2 | ~18 | ~18 |

Text [M]: −70 % at equal potato and fish versus a purely additive (dilution) expectation. [Z]
pan: 27 vs ~95 expected → −71 %; oven: 101 vs ~388 → −74 %.

### 3.4 Table 1 — amino acids and sugars, Frieslander batch 1, 180 °C, 25 min [M]

Single samples; corrected for weight loss (66.6 ± 1.6 %).

| addition | amount (mmol/kg) | measured pH | acrylamide (µg/kg) | change (%) |
|---|---|---|---|---|
| none (control I) | — | 5.82 | 427 | |
| none (control II) | — | 5.74 | 542 | |
| none (control III) | — | 5.82 | 497 | |
| mean (± SD) controls | | | 488 (± 58.0) | |
| glutamine | 35 | 5.79 | 115 | −76 |
| | 140 | 5.72 | 37 | −92 |
| glycine | 35 | 5.78 | 146 | −70 |
| | 140 | 5.79 | 44 | −91 |
| 4-aminobutyric acid | 35 | 5.81 | 262 | −46 |
| | 140 | 5.86 | 67 | −86 |
| lysine (acetate salt) | 35 | 5.83 | 276 | −43 |
| | 140 | 5.99 | 59 | −88 |
| sodium glutamate | 35 | 5.98 | 463 | −5.3 |
| | 140 | 6.16 | 226 | −54 |
| alanine | 35 | 5.73 | 422 | −14 |
| | 140 | 5.76 | 243 | −50 |
| asparagine | 35 | 5.71 | 1550 | +220 |
| | 140 | 5.61 | 1885 | +290 |
| fructose | 35 | 5.77 | 2720 | +460 |
| | 140 | 5.81 | 2060 | +320 |
| glucose | 35 | 5.76 | 1280 | +160 |
| | 140 | 5.72 | 691 | +41 |
| glucose-6-phosphate | 35 | 5.56 | 1180 | +140 |
| glyceraldehyde | 15 | 5.77 | 765 | +56 |
| | 45 | 5.72 | 1080 | +120 |

Text [M]: mean of the five non-Asn amino acids at 35 mM, −42 ± 29 % (SD, n = 5); at 140 mM, −75 ±
19 %. Asn additions raise free Asn ×3 and ×9 over the endogenous 17 mM. Glucose ×3 and ×8 over
21 mM; fructose ×14 and ×50 over 2.8 mM.

### 3.5 Table 4 — antioxidants, peroxides, trisodium citrate, cv. Ukama, 180 °C, 25 min [M]

| addition | amount (mmol/kg) | measured pH | wt loss (%) | acrylamide (µg/kg) | change (%) |
|---|---|---|---|---|---|
| none (control I) | — | 5.70 | 68.8 | 407 | |
| none (control II) | — | 5.78 | 70.3 | 477 | |
| none (control III) | — | 5.92 | 67.3 | 396 | |
| none (control IV) | — | 5.89 | 67.7 | 476 | |
| mean ± SD controls | | | | 439 (± 43.5) | (± 10) |
| sodium ascorbate | 1.51 | 5.76 | 67.3 | 494 | +12 |
| | 4.54 | 5.81 | 64.6 | 429 | −2.3 |
| | 13.6 | 6.05 | 66.4 | 379 | −14 |
| | 40.9 | 5.88 | 65.6 | 422 | −3.9 |
| | 123 | 6.01 | 61.5 | 239 | −46 |
| ascorbyl palmitate | 1.51 | 5.86 | 65.1 | 569 | +30 |
| | 4.54 | 5.8 | 65.5 | 499 | +14 |
| | 13.6 | 5.57 | 65.1 | 491 | +12 |
| benzoyl peroxide | 1.03 | 5.73 | 68.4 | 640 | +46 |
| | 4.13 | 5.73 | 67.8 | 501 | +14 |
| hydrogen peroxide | 6.0 | 5.77 | 62.3 | 453 | +3.2 |
| | 18 | 5.6 | 67.3 | 510 | +16 |
| | 54 | 5.74 | 64.8 | 452 | +3.0 |
| | 162 | 5.66 | 68.3 | 470 | +7.1 |
| trisodium citrate | 4.54 | 5.96 | 65.1 | 442 | +0.7 |
| | 13.6 | 6.06 | 64.1 | 537 | +22 |
| | 40.9 | 6.29 | 62.2 | 283 | −36 |
| | 123 | 6.53 | 60.4 | 208 | −53 |

Authors [Q]: the −46 % (ascorbate) and −53 % (citrate) at 123 mM are attributed to water binding
(lower weight loss, less browning), not antioxidant action. Note trisodium citrate **raises** pH
(to 6.53) yet lowers acrylamide, the opposite sign to the NaOH series. This is a matrix/water effect
that confounds any pH reading of these rows.

### 3.6 Table 5 — authors' summary [Q]

Higher T: +. Longer time at high T: + then −. pH: maximum ~8. Monosaccharides: + (fructose >
glucose ≈ G6P), saturating. Asn: +, saturating. Other amino acids: −. Antioxidants/radical
initiators: (+)?, minor. Water-binding additives: −. Lean fish: strong −.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context (read, not run):** `acrylamide_conditions.py` applies a declared initial-pH
factor 10^(e·(pH − 6.8)) to `k_asn_glc` (e = 0.235 dec/pH) and `k_acr_dp` (e = 0.149), from De
Vleeschouwer 2006 phosphate model systems. The factor is clamped to the window pH 4–8. The backlog
(`tasks/data_restructure_plan.md` §7) names Rydberg's pH-8 maximum as "the deciding test" for
acrylamide speciation.

**Can anchor (as held-out directional/shape tests, not as fit targets for rate constants):**
1. **A net-content maximum in initial pH between ~7.4 and ~8.3**, at 160 and 180 °C, 25 min, in
   a real potato matrix, with a decline of 30–45 % from the peak to pH 9.7. A monotone power law in
   [H⁺] cannot itself produce a maximum, and the model's clamp at pH 8 makes its pH > 8 response
   flat. Whether the full lane (formation × elimination × Asn speciation) produces a maximum is a
   question for a scored run, not answered here.
2. **Acidification direction and size on the low side:** net content falls ×2.6 from pH 5.7 to
   3.0 (180 °C) and ×8.5 at 160 °C. Net d ln C/d pH ≈ 0.35–0.56 at 180 °C (§2.5). These are the
   same order as De Vleeschouwer's ln-slopes (formation 0.54, potato 0.43), but are net quantities.
3. **Elimination direction [!] — a potential contradiction:** spiked acrylamide is lost **faster at
   LOW pH** (47 % at pH 3.0 vs 22 % at pH 5.7, 180 °C, 25 min). The model's elimination factor
   **rises with pH** (+0.149 dec/pH). These are different systems (potato matrix with net
   formation + spike vs phosphate kinetics), and spike "loss" also includes reaction with matrix
   nucleophiles. The sign disagreement should be recorded, not resolved by tuning.
4. **Time course** at 200 °C (peak ~20 min, −76 % by 60 min) and **temperature** response (Fig. 2,
   Fig. 5), as shape checks for formation/elimination balance in a potato-strip or homogenate
   matrix.
5. **Additions** (Table 1): Asn and sugar saturation, and the strong suppression by other amino
   acids at 35–140 mM. These are directional claims for a competitive-amine term. Endogenous Asn
   17.0 mM, glucose 21 mM and fructose 2.8 mM make the matrix composition declarable.

**Cannot anchor:**
- **Rate constants or Ea.** Single samples, net contents, an oven ramp at 40 °C/min with frozen
  45 g cubes (internal T unknown), and simultaneous drying (60–70 % weight loss). Nothing here
  separates formation from elimination except Table 2's single spike experiment.
- **In-situ pH.** All pH values are cold, pre-heating readings of the homogenate surface layer. The
  pH during heating and drying is not known. Any mapping to the model's in-situ pH needs the model's
  own pH state, and is a declared assumption.
- **Absolute levels across batches.** A factor ~2 between Frieslander batches at nominally equal
  conditions (Fig. 4 vs Fig. 5, 160 °C) and 275–600 µg/kg control spread across cultivars. Score
  ratios within a series, not absolute values.
- **pH > 8 with buffering composition.** NaOH addition raises ionic strength and Na⁺ and changes
  water binding (cf. trisodium citrate's −53 % at higher pH). The high-pH decline may be partly a
  matrix effect.
- **Fig. 2 and Fig. 3 pan conditions.** They come from Tareke 2002 (ref 2) and are not restated here.

**Caveats.** All figure values are my digitisation (errors as stated per figure). Fig. 4 is the
most reliable (automated, R² reproduced). The Fig. 4 cultivar is inferred, not printed.
`claeys2005_extraction.md` quotes Rydberg's lysine effect as −57 % with "pH not stated". Table 1
prints **−43 % (35 mM) and −88 % (140 mM)** and gives pH for every row. The discrepancy is either in
Claeys's citation or in that dossier, and is not resolved here.

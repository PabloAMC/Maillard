# Troise, Fiore, Wiltafsky & Fogliano 2015 — EXTRACTION (one-hydrolysate SIDA LC-MS/MS for lysine, furosine, CML and CEL; validation numbers; expanded soybean autoclaved at 110 C for 0-60 min)
### The analytical method behind Nguyen 2016 and the only soy-protein CML/CEL/furosine time course on disk — a methods paper, no rate constants, the profile itself is figure-only.

**Source on disk:** `data/articles/troise2015.pdf` (owner's download, 2026-09-08; Elsevier accepted manuscript,
28 pages, Food Chemistry, PII S0308-8146(15)00707-4, DOI 10.1016/j.foodchem.2015.04.137; the version of
record is Food Chem. 188 (2015) 357-364 — NOT verified against this file, which carries no volume or page
numbers). Read from the `pdftotext` text layer in the scratchpad; Tables 1-3 are clean and re-typed in
full below. Figures 1-3 are images without a text layer (PDF pp. 25-27 hold only the words "Figure-1",
"Figure-2", "Figure-3"): the soybean kinetic profile of Fig. 3 is **FIGURE-ONLY**; the numbers quoted from
it in the Results text are transcribed in section 3. Repo status before this dossier: not cited anywhere
(grep of `data/lit`, `src/`, `results/`); `nguyen2016_extraction.md` (written the same evening) cites the
method through Troise et al. 2014.

## 0. Identity

| field | value |
|---|---|
| Title | "Quantification of Nε-(2-Furoylmethyl)-L-lysine (furosine), Nε-(Carboxymethyl)-L-lysine (CML), Nε-(Carboxyethyl)-L-lysine (CEL) and Total Lysine through Stable Isotope Dilution Assay and Tandem Mass Spectrometry" |
| Authors | Antonio Dario Troise (Wageningen / Napoli Federico II, corresponding), Alberto Fiore (Abertay), Markus Wiltafsky (Evonik Industries), Vincenzo Fogliano (Wageningen) |
| Venue | Food Chemistry, accepted 29 Apr 2015 (received 14 Jan 2015, revised 28 Apr 2015); section "Analytical methods" |
| Type | analytical-method paper with one industrial kinetic demonstration (soybean feed) |
| Naming | "AP" = Amadori product; furosine is its acid-hydrolysis marker; "MRPs" = Maillard reaction (end) products; concentrations in mg per 100 g protein (lysine in g per 100 g protein) |

## 1. Why it matters

Two things. First, it is the method Nguyen 2016 used (through Troise 2014) and the template for what an
isotope-dilution CML/CEL/lysine measurement looks like: the repository's evidence-class rules for AGE rows
(the `cml_cel_commercial_pbma_Foods2023` benchmark records `quantification_class:
isotope_dilution_lcmsms`) can point at this paper's Table 2 for the LOD, LOQ, linearity, RSD and recovery
that "isotope dilution" implies. Second, section 3.4 gives the only protein-bound-lysine time course in a
SOY matrix on disk — lysine loss, furosine rise and fall, CML rise, CEL rise and fall over 60 min at
110 C in an autoclave — which is the closest thing in the corpus to Programme 7's pea/soy isolate being
glycated. Its usefulness is bounded by what is printed: six numbers from the curves, no rate constants,
no model, no replicate table for the kinetic run.

## 2. Methods as they matter to a model

- **Soybean sample and heat treatment (section 2.2.1).** One batch of quartered raw soybeans (Rieder
  Asamhof, Kissing) processed at the Amandus Kahl hydrothermal plant: conditioned to 80 C in 45 s;
  hydrothermal belt cooker 72 C entry, 3 min at 70 C; expanded at 117 C in an annular-gap expander;
  10 min in a drying wagon; dried with 65 C air 10 min, cooled 10 min to **12 % moisture**. Then
  **autoclaved at 110 C and 1470 mbar for 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50 and 60 min** (Zirbus
  HST 6x9x12). So: a whole-seed feed (~35-40 % protein, ~20 % oil, sucrose + raffinose-family sugars — none
  of this measured or printed), 12 % moisture, saturated steam at 110 C (a_w near 1 at the surface,
  seed interior at its 12 % moisture a_w, not measured), no buffer, no pH. Twelve time points; number of
  replicates for the kinetic run not stated separately (2.6 says all analyses in quadruplicate; Table 3's
  caption says 8 replicates for the food survey).
- **Hydrolysis (2.3).** 100 mg sample + 4 mL 6 N HCl in a PTFE-capped flask; N2-saturated 15 min at 2 bar;
  110 C for 20 h; PVDF 0.22 µm filtration; 400 µL dried under N2; reconstituted in 370 µL water + 10 µL
  of each internal standard (**d4-lysine, d2-CML, d4-CEL**) to a final **200 ng per mg sample**; Oasis HLB
  1 cc SPE; 5 µL injected. **No sodium borohydride reduction** — tried and abandoned (3.4: "After several
  preliminary measurements it was decided to avoid the reduction"), because reduction converts
  fructosyl-lysine to hexitol-lysine and so suppresses furosine; the N2 saturation is the guard against
  AP -> CML conversion during hydrolysis. Consequence: CML may carry an analytical contribution from AP in
  AP-rich samples (the soybean run at 30 min has AP at its maximum).
- **LC (2.4).** Kinetex C18 core-shell 2.6 µm, 2.1 x 100 mm; A = 5 mM perfluoropentanoic acid in water,
  B = 5 mM PFPA in acetonitrile; 200 µL/min; gradient (min / % B): 0/10, 2/10, 5/70, 7/70, 9/90, 10/90,
  12/10, 15/10. Retention: CML and d2-CML 7.11 min, lysine and d4-lysine 7.23, CEL and d4-CEL 7.36,
  furosine 7.91 (shift < 0.5 min over a batch).
- **MS (2.4, Table 1).** API 3000 triple quadrupole, ESI+, 5.0 kV, 350 C, dwell 100 ms, CAD 45, curtain 5;
  MRM. Quantifier / qualifier: CML 205 -> 84.1 / 205 -> 130.2; CEL 219.1 -> 84.1 (both listed as 84.1 in
  the text; Table 1 gives 130.3 and 84.0); furosine 255.1 -> 130.2 / 255.1 -> 84; lysine 147.2 -> 130.2 /
  147.2 -> 84.1. Internal standards: d2-CML 207 -> 144.1 (quant) / 207 -> 84; d4-CEL 223 -> 134.1 /
  223 -> 84 (Table 1: 88.0); d4-lysine 151.2 -> 134.1 / 151.1 -> 88. **Furosine has no labelled analogue and
  is quantified against d4-lysine** — justified by an infusion test in which furosine (m/z 255) and
  d4-lysine (151.2) gave similar intensities, both ~15 % above d2-CML and d4-CEL.
- **Calibration (2.5, 3.3).** Linear curves of analyte / internal-standard area ratio against
  concentration, standards in water, internal standards at 200 ng/mL; linearity 5-1000 ng/mL (CML, CEL,
  lysine) and 9-1000 ng/mL (furosine); r^2 > 0.99; intraday x3 and interday x3 days, RSD < 8 %; carry-over
  checked after each calibration point; recovery from the internal-standard intensity in every matrix
  (Table 2). The paper reports the analyte / IS ratio matched to the curve in each sample, with the IS
  intensity in sample vs standard agreeing within 10 % RSD.
- **Units and conversions.** CML, CEL, furosine in mg per 100 g protein; lysine in g per 100 g protein
  (Table 3's UHT-milk lysine "4.71 ± 0.22 mg/100 g of protein" in the text is a unit slip for g). Protein
  content of the soybean feed is not printed; per-protein numbers therefore cannot be turned into per-kg
  seed. Molar conversions per g protein: lysine (146.19 g/mol) 1 g/100 g = 68.4 µmol/g; furosine (254.28)
  1 mg/100 g = 0.0393 µmol/g; CML (204.22) 1 mg/100 g = 0.0490 µmol/g; CEL (218.25) 1 mg/100 g =
  0.0458 µmol/g. No furosine -> AP factor is applied in this paper (Nguyen 2016 uses 3.1, Berk 2021 uses
  2.2).

## 3. Tables re-typed

### Table 1. "Mass spectrometry set up" (all values as printed)

| compound | [M+H]+ | fragments | CE (V) | DP (V) |
|---|---|---|---|---|
| CML | 205 | 84 | 29 | 30 |
| | | 130.2 | 27 | 30 |
| d2-CML | 207 | 84 | 30 | 20 |
| | | 144 | 21 | 20 |
| | | 130 | 17 | 20 |
| Furosine | 255.1 | 130 | 18 | 21 |
| | | 84.4 | 28 | 21 |
| Lys | 147.2 | 130.2 | 16 | 30 |
| | | 84.1 | 24 | 30 |
| d4-Lys | 151.3 | 134.1 | 15 | 30 |
| | | 88.2 | 26 | 30 |
| CEL | 219.2 | 130.3 | 20 | 30 |
| | | 84.0 | 28 | 30 |
| d4-CEL | 223 | 134.1 | 18 | 25 |
| | | 88.0 | 30 | 25 |

### Table 2. "Analytical performances for the four analytes and their respective internal standards"

| compound | LOD | LOQ | RSD (%) | linearity range | r^2 | recovery (%) |
|---|---|---|---|---|---|---|
| CML | 0.5 ppb | 5 ppb | 7 | 5-1000 ng/mL | > 0.99 | 91.1 ± 8.4 |
| CEL | 1 ppb | 5 ppb | 5 | 5-1000 ng/mL | > 0.99 | 84.2 ± 7.4 |
| Lysine | 0.5 ppb | 5 ppb | 5 | 5-1000 ng/mL | > 0.99 | 88.0 ± 6.9 |
| Furosine | 3 ppb | 9 ppb | 8 | 9-1000 ng/mL | > 0.99 | 88.0 ± 6.9 |

ppb = ng/mL in the injected solution. Furosine's recovery is the d4-lysine figure repeated (it has no
labelled standard). Text (3.3): 0.1 ppb gave no signal; S/N > 3 at the LODs; recoveries "91.1 ± 8.4,
84.2 ± 7.4, 88.0 ± 6.9 for d2-CML, d4-CEL and d4-lysine" across all matrices.

### Table 3. "MRPs concentration after 8 replicates in different samples" — CML, CEL, furosine in mg / 100 g protein; lysine in g / 100 g protein; "Age Database" rows = TU Dresden AGE database ranges quoted for comparison

| food | CML | CEL | furosine | lysine (g / 100 g protein) |
|---|---|---|---|---|
| Infant formula 1 | 8.22 ± 0.31 | 0.71 ± 0.02 | 471.91 ± 22.31 | 9.89 ± 0.88 |
| Infant formula 2 | 10.4 ± 0.52 | 0.85 ± 0.06 | 542.53 ± 11.91 | 12.24 ± 0.91 |
| Infant formula 3 | 10.9 ± 1.03 | 1.10 ± 0.05 | 574.5 ± 44.12 | 13.12 ± 0.78 |
| Infant formula 4 | 14.81 ± 0.92 | 1.31 ± 0.11 | 639.4 ± 21.11 | 10.28 ± 1.01 |
| AGE database (infant formula) | 0.6-40.5 | / | up to 1819 | / |
| Low lactose milk | 1.28 ± 0.11 | 0.28 ± 0.01 | 12.32 ± 0.31 | 5.21 ± 0.30 |
| AGE database (milk) | 1.4 | / | / | / |
| Lab-scale UHT milk | 18.41 ± 0.93 | 1.12 ± 0.02 | 14.41 ± 1.02 | 4.71 ± 0.22 |
| AGE database (UHT milk) | 0.9-8.3 | / | 12.4-220.0 | / |
| Biscuits | 43.75 ± 2.02 | 46.25 ± 3.01 | 10.01 ± 0.61 | 5.01 ± 0.04 |
| Bread slices | 27.15 ± 0.61 | 10.91 ± 0.01 | 98.55 ± 4.61 | 5.81 ± 0.04 |
| AGE database (bakery) | 2.6-45.1 | / | / | / |

Biscuit and bread protein contents 6 % and 8 % (text). Bread: 20 min at 200 C (text, in the furosine
comparison with Capuano 2008). The infant-formula lysine values (9.9-13.1 g / 100 g protein) exceed the
lysine content of any milk protein (~8 g / 100 g) — the paper does not comment; take Table 3's lysine
column as method output, not composition.

### The soybean kinetic profile (Fig. 3, FIGURE-ONLY) — every number printed in the text (section 3.4)

Expanded soybeans, 110 C autoclave, 1470 mbar, 12 % initial moisture; mg per 100 g protein unless stated.

| time (min) | lysine (g / 100 g protein) | furosine | CML | CEL | note |
|---|---|---|---|---|---|
| 0 | 3.45 ± 0.12 | 24.24 ± 1.74 | 9.94 ± 0.74 | 0.98 ± 0.04 | "initial concentration" (after expansion at 117 C — not raw) |
| 30 | — | **108.01 ± 8.97** (maximum) | — | — | "After 30 minutes the concentration of furosine reached the highest values" |
| 45 | — | — | — | **2.41 ± 0.24** (maximum) | "CEL reached the maximum concentration after 45 minutes ... then it decreased" |
| 55 | — | 60.58 ± 3.75 | — | — | "rapidly decreased up to 60.58 ± 3.75 mg/100 g of protein after 55 min" |
| 60 | 2.60 ± 0.08 | — | "> 76" | — | "at the end of the thermal treatment its concentration was higher than 76 mg/100 g of protein"; lysine loss "around 23 %" (3.45 -> 2.60 is 24.6 %) |

Shape statements: lysine "degradation ... was constant throughout the thermal treatment" (i.e. roughly
linear over 60 min); furosine rise to 30 min then fall; CML rising "according to the reaction mechanism the
degradation of the Amadori products was followed by the increase of CML"; CEL rise to 45 min then fall
("probably due to degradation processes or to the blockage of methylglyoxal by other compounds").

Molar form (per g protein): lysine 236 -> 178 µmol/g (−58 µmol/g); furosine 0.95 -> 4.25 (30 min) ->
2.38 µmol/g (55 min); CML 0.49 -> > 3.7 µmol/g; CEL 0.045 -> 0.110 µmol/g. So of the ~58 µmol/g lysine
lost in 60 min, CML accounts for ~3.2 µmol/g (~6 %), CEL for < 0.1 µmol/g, and furosine-measurable AP at
its 30-min peak for 3.3 µmol/g as furosine (x 2.2-3.1 as AP: 7-10 µmol/g, 12-18 % of the lysine loss).
The balance is unmeasured (other AP fates, crosslinks, hydrolysis-resistant adducts). This paragraph is my
arithmetic on printed numbers, not the authors'.

## 4. Kinetic numbers the repository can use

Registry keys: CML -> `cml`; CEL -> `cel`; furosine -> `furosine`; lysine (protein-bound) -> not in
registry (`reactive_lysine` is a marker set). No glyoxal, methylglyoxal or 3-deoxyglucosone measured.

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| protein-bound lysine loss (all routes) | fraction lost in 60 min | 0.246 (printed "around 23 %") | — | expanded soybean, 110 C saturated steam, 12 % moisture | 3.4 | level_only |
| protein-bound lysine loss | apparent first-order k = ln(3.45/2.60)/60 | 4.7e-3 | min^-1 | 110 C, as above; two points; authors describe the decline as linear (zero order: 0.97 µmol g^-1 min^-1) | derived by me from 3.4 | within_study_ratio (derived, two-point) |
| AP (as furosine) rise | furosine 24.2 -> 108.0 mg / 100 g protein in 30 min | +3.3 µmol/g protein (x 2.2-3.1 for AP) | — | 110 C | 3.4 | level_only |
| AP (as furosine) fall | 108.0 -> 60.6 in 25 min: ln(108.0/60.6)/25 | 2.3e-2 | min^-1 (apparent, net of continuing formation) | 110 C, 30 -> 55 min | derived by me | within_study_ratio (derived, two-point, lower bound on the AP loss constant) |
| CML formation | 9.94 -> > 76 mg / 100 g protein in 60 min | > 3.2 µmol/g protein | — | 110 C | 3.4 | level_only (end value is a bound) |
| CEL formation | 0.98 -> 2.41 (45 min), then falls | +0.065 µmol/g protein | — | 110 C | 3.4 | level_only |
| CEL / CML at the CEL maximum | ≈ 2.41 / (between 10 and 76) | 0.03-0.2 | — | 110 C, 45 min | 3.4 | within_study_ratio (bounded) |
| furosine max / lysine(0) | 4.25 / 236 | 1.8 % of lysine sites as furosine (4-6 % as AP) | — | 110 C, 30 min | derived by me | within_study_ratio |
| whole time course (12 points, 4 analytes) | — | mg / 100 g protein vs min | 110 C | Fig. 3 | figure_only |
| Method: LOD / LOQ / RSD / recovery | Table 2 | see section 3 | ng/mL; % | — | Table 2 | measured (analytical) |

Nothing here is a rate constant of a named step. The two derived constants above are two-point apparent
rates from a whole seed and should not enter a parameter file; they are the right size to check a
Programme 7 simulation against (a 25 % lysine loss in 60 min at 110 C in a moist protein solid).

## 5. Flags

1. **Figure-only kinetics.** Fig. 3 has twelve time points per analyte; only the six numbers above are
   printed. The 60-min CML value is a bound ("higher than 76"). No replicate spread is printed for the
   kinetic run other than the ± on those six numbers.
2. **The lysine value is low for soy.** 3.45 g / 100 g protein at t = 0 against ~6.0-6.4 g / 100 g protein
   in soybean meal composition tables; the authors themselves say that protein-to-HCl ratio can cause
   "underestimation of lysine content" and that the hydrolysis was tuned for furosine release, not lysine
   release. Use the fractional loss (25 %), not the absolute lysine, and expect the CML / CEL per-protein
   numbers to be on a consistent but possibly shifted scale.
3. **The t = 0 sample is already processed** (expanded at 117 C, hydrothermally cooked): furosine 24 and
   CML 9.9 mg / 100 g protein before the autoclave. A model starting from a raw isolate should start lower.
4. **No NaBH4 reduction**: CML may include some AP converted during hydrolysis, most at the 30-min AP peak.
   The paper's own Table 3 comparison with literature is its argument that the effect is small.
5. **Furosine is quantified against d4-lysine** (no labelled furosine); the 15 % response-factor argument is
   from a single infusion at 10 ppm. Furosine is the analyte with the largest RSD (8 %) and LOD (3 ppb).
6. **Matrix, moisture, pressure**: 1470 mbar autoclave = saturated steam at ~110-111 C; the seed starts at
   12 % moisture and will take up water; a_w in the seed is not measured. "110 C" is the chamber, not a
   guaranteed seed-interior temperature during the first minutes.
7. **Protein content of the feed not printed** — per-protein numbers cannot be turned into per-kg seed,
   and sugar content is not measured (the reducing-sugar side of the reaction is invisible).
8. **Replicate count ambiguity**: 2.6 says quadruplicate; Table 3's caption says 8 replicates.
9. **Accepted manuscript**, not the version of record; table values may have been edited in production.
10. **Small text/table inconsistencies** in the MS transitions (CEL qualifier listed as 84.1 in the text and
    130.3 / 84.0 in Table 1; d4-CEL confirmation 84 in text, 88.0 in Table 1) — analytical detail only.

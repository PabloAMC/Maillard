# Thomas, Mercier, Tournayre, Martin & Berdagué 2014 — EXTRACTION (thiamine vs cysteine + sugar as sources of 2-methyl-3-furanthiol and its disulfides in model cooked ham, 69 °C)

**Source on disk:** `data/articles/thomas2014.pdf` (the publisher's PDF, 7 pages, with a text layer). Read and
checked by eye on 2026-10-09: every number below was read from the page image and, where the page has
text, cross-checked against `pdftotext -layout`. Fig. 2, 3 and 4 have no data table; their values were
read from the embedded figure images (extracted at native 300 ppi and zoomed) and are labelled
"read from graph, approx.". Written for the thiamine-vs-cysteine question raised by the composition sweep
(thiamine predicted to matter for MFT at 100 °C/20 min but not at 140 °C/5 min), not for the core fit.

| field | value |
|---|---|
| Title | "Identification and origin of odorous sulfur compounds in cooked ham" |
| Authors | Caroline Thomas, Frédéric Mercier, Pascal Tournayre, Jean-Luc Martin, Jean-Louis Berdagué (INRA UR 370 QuaPA, Saint-Genès-Champanelle; IFIP, Maisons-Alfort) |
| Venue | Food Chemistry 2014, 155, 207-213 |
| DOI | 10.1016/j.foodchem.2014.01.029 |

## 1. Methods

**Two experiments.** (i) Identification of sulfur volatiles in one commercial cooked ham (Table 1, Fig. 1):
SPME-GC×GC-TOFMS and dynamic headspace (DHS) GC-MS with 8-way and single-port olfactometry. (ii) The
precursor experiment (Fig. 2, 3, 4), which is the part that matters here.

**Matrix and additions (§2.1, p. 207-208).** Model hams of 300 g from pig *M. semimembranosus*. 1 kg chopped
raw meat mixed with 100 g brine (18 g curing salt containing 0.6 % nitrite, in 82 g water). Three dose
series, amounts "added to 100 g of brine":

| series | precursor doses (mg per 100 g brine) | co-added |
|---|---|---|
| T | thiamine 0, 8, 80, 800, 8000 | none |
| C + F | cysteine 0, 100, 1000, 3000, 10,000 | fructose 3 g per 100 g brine |
| C + X | cysteine 0, 100, 1000, 3000, 10,000 | xylose 3 g per 100 g brine |

Conversion printed in the Fig. 2, 3 and 4 captions: mg per 100 g brine × 0.91 = mg added per kg of ham before
cooking. So (derived here, × 0.91): thiamine 7.3, 73, 728, 7280 mg/kg ham; cysteine 91, 910, 2730,
9100 mg/kg ham; fructose or xylose 2730 mg/kg ham. Reagent form: thiamine "purity 99%, Ref.: T4625"
(Sigma); the salt form is not printed in the paper. Cysteine purity 97 % (Sigma W326305).

**Cooking (§2.1, p. 207).** Closed glass jars under vacuum; "40–69 °C at 0.1 °C/min⁻¹ and 120 min hold at
69 °C". Ramp duration derived here: (69 − 40) / 0.1 = 290 min, so about 410 min total above 40 °C. pH:
not printed. Water activity, moisture: not printed.

**Native thiamine (§2.6 and p. 212).** AFNOR NF EN 14122 (acid hydrolysis, enzymatic dephosphorylation,
HPLC, thiochrome fluorimetry), so total thiamine. Raw pork 9.5 mg/kg; 8.6 mg/kg after brining; 7.5 mg/kg
after cooking (control ham). Only the 0, 8 and 80 mg/100 g brine hams were assayed.

**Volatile measurement (§2.4, p. 208).** 7 g of minced ham in a Pyrex cartridge, DHS at 30 °C for 60 min,
helium 40 ml/min, Tenax TA trap; GC-MS (RTX5-MS 60 m) quadrupole. For the precursor trials "the sulphur
odourants of interest were semi-quantified by measuring their peak areas from specific ions acquired in
single ion monitoring mode". Ions (Fig. 2 caption): MFT m/z 114, 2-methyl-3-(methyldithio)furan m/z 160,
bis(2-methyl-3-furyl) disulfide m/z 113, 2-methylthiophene m/z 97. **No internal standard is mentioned and
no calibration**: values are arbitrary units (a.u.) of chromatographic area. "Analyses of model hams were
performed in triplicate" (p. 208); whether three hams or three injections of one ham is not stated.
Fig. 2 and 4 plot median with minimum and maximum; no statistical test is reported for any volatile. The
only p-value in the paper is for thiamine loss on cooking (p < 0.05, Fig. 3).

**GC-O (§2.5).** On the commercial ham only (2 sessions × 8 sniffers on DHS-GC-MS/8O; 2 assessors on the
single-port and heart-cut instrument). Not done on the precursor hams.

2-Furfurylthiol (2-furylmethanethiol) is **not reported anywhere in the paper** (not in Table 1, not in
Fig. 2 or 4).

## 2. Findings that matter

### 2a. Fig. 2 (p. 211): sulfur odorants vs precursor dose, model hams

All values **read from graph, approx.**, in the printed unit "a.u." with the axis multiplier "x 10000"
printed on each panel (so a reading of 40 means 40 × 10⁴ a.u.). Median (min–max where the bar is
visible). The 0-dose point is drawn as a common starting point for all three series. x-axis is a broken
log scale; thiamine points sit at 8, 80, 800, 8000 and cysteine points at 100, 1000, 3000, 10,000.

| compound (ion) | control (0) | T 8 | T 80 | T 800 | T 8000 | C+F and C+X, 100 to 10,000 |
|---|---|---|---|---|---|---|
| 2-methyl-3-furanthiol (m/z 114) | ~0.3 | ~0.8 | ~1.0 | ~1.8 (max ~8.5) | ~40 (~37–52) | ~0.3–0.5 at every dose, not distinguishable from control |
| 2-methyl-3-(methyldithio)furan = methyl 2-methyl-3-furyl disulfide (m/z 160) | ~0 | ~20 (~1–39) | ~11 (max ~32) | ~27 (~5–92) | ~188 (~101–242) | ~0; at most ~3–4 at 1000 |
| bis(2-methyl-3-furyl) disulfide (m/z 113) | ~0 | ~0.3 | ~0.5 | ~0.8 (max ~3.5) | ~27.5 (~12–53) | ~0 |
| 2-methylthiophene (m/z 97) | ~125 | ~620 (~530–700) | ~620 (~600–680) | ~945 (max ~1070) | ~2570 | C+F ~125 at 100, ~10 at 1000 and 3000, ~0 at 10,000; C+X ~60 at 1000, ~0 at 3000 and 10,000 |

Text (p. 211-212), paraphrased: thiamine "at the highest doses clearly induced" all four; MFT and the
symmetric disulfide rise clearly "only from significant enrichments of thiamine (higher than 80 mg of
thiamine per 100 g of brine)"; cysteine "induced no significant production" of MFT, the methyl disulfide
or the symmetric disulfide whichever sugar was co-added, and lowered 2-methylthiophene at very high doses.
Authors' explanation: "The cooking of ham at 69 °C for several hours seems insufficient to induce Maillard
reactions forming compounds from cysteine." "Significant" here is not backed by a test.

Fold changes, **derived here** from the graph readings: MFT at T 8000 vs control ≈ 40 / 0.3 ≈ 130×, but
the control sits on the axis, so the ratio is only order-of-magnitude (≥ 50×). At T 8 (which, derived
here, raises total thiamine from 8.6 to about 8.6 + 7.3 = 15.9 mg/kg, ≈ 1.85×; Fig. 3 reads ~16 before
cooking) MFT moves from ~0.3 to ~0.8 × 10⁴ a.u., inside the scatter of the higher points; the paper does
not claim it as an effect.

### 2b. Fig. 4 (p. 212): other compounds (read from graph, approx., × 10⁴ a.u.)

| compound (ion) | control | T 8 | T 80 | T 800 | T 8000 | C+F / C+X |
|---|---|---|---|---|---|---|
| 2-methyl-5-(methylthio)furan (m/z 114) | ~0 | ~7 | ~12 | ~72 | ~480 (~515 max) | ~0 |
| 2-methyl-4,5-dihydrothiophene (m/z 97) | ~0 | ~2 | ~3 | ~14 | ~82 | ~0 |
| dimethyl sulfide (m/z 62) | ~1.9 | ~4.8 | ~5.5 | ~4.5 | ~4.0 | C+F ~1.9, ~2.7, ~1.5, ~0.9 at 100/1000/3000/10,000; C+X ~1.8, ~1.7, ~1.0, ~2.5 |

### 2c. Fig. 3 (p. 212): thiamine in the hams (read from graph, approx.)

| thiamine added (mg/100 g brine) | before cooking | after cooking |
|---|---|---|
| 0 | ~8.7 | ~7.6 |
| 8 | ~15.8 | ~14.8 |
| 80 | ~81 (SE bar ~69–93) | ~57 (SE bar ~46–69) |

Units: the y-axis reads "Concentration of thiamine (mg/kg of ham)" but the caption says "(mg/100 g)". The
text numbers (9.5, 8.6, 7.5 mg per kg) match the axis, so mg/kg is taken here as correct. Text: a loss of
"about 30% (p < 0.05)" only at the 80 mg dose (derived here from the readings: (81 − 57)/81 ≈ 30 %); at
lower doses assay precision was insufficient. Authors: thiamine "is only partially consumed during
cooking".

### 2d. Table 1 (p. 209-210): commercial ham, for context only

Relative peak area, % of total sulfur-compound area (SPME-GC×GC-MS / DHS-GC-MS), and GC-O mean
intensity (1–5): MFT tr / 0.5, intensity 4.2; 2-methyl-3-(methyldithio)furan tr / 3.8, 4;
bis(2-methyl-3-furyl) disulfide tr / 0.1, 4 (smelt only on the single-port devices); dimethyl trisulfide
1.3 / 8.1, 4; 2-methylthiophene 0.1 / 16.2, 2.6. "tr" = peak area < 0.1 %. Relative, not absolute.

### 2e. Internal inconsistencies noticed

- p. 211: "From 10 mg of added thiamine, production of both compounds was observed"; the methods and the
  plotted points give 8 mg, not 10 mg.
- p. 211 text gives 2-methylthiophene a mean intensity "greater than 3.5/5"; Table 1 prints 2.6.
- Fig. 3 caption unit (mg/100 g) vs axis and text (mg/kg), above.

## 3. What it means for the model

At 69 °C core (slow ramp, 2 h hold), in a real pork matrix with nitrite:

- **Thiamine is the MFT source; cysteine (+ fructose or xylose) is not.** Cysteine up to 9100 mg/kg ham with
  2730 mg/kg reducing sugar gave no visible MFT, methyl 2-methyl-3-furyl disulfide or bis(2-methyl-3-furyl)
  disulfide. The sweep's claim is about 100 °C; this paper sits 31 °C lower. It is consistent in
  direction with the sweep's trend (thiamine weighs more against cysteine as temperature falls), and
  stronger than it: at 69 °C cysteine contributes nothing measurable, not "as much as thiamine".
- **The thiamine effect is only clear at supraphysiological doses.** A clear MFT and bis-disulfide rise
  starts above 80 mg/100 g brine (≈ 73 mg/kg added, ≈ 8.5× the native 8.6 mg/kg, derived: 72.8 / 8.6) and is large at 8000
  (≈ 7280 mg/kg, ≈ 850× native, derived: 7280 / 8.6). The methyl disulfide and 2-methylthiophene respond
  from the lowest dose. Near native levels (the "restore thiamine" case), the MFT change is not resolved.
- Thiamine is only partly consumed (about 30 % loss at the 80 mg dose), so at 69 °C the pool is not
  exhausted, unlike the beef pan-cooking in Lombardi-Boccia 2005 where cooked thiamine was not detected.

**Strength of evidence: weak-to-moderate, qualitative.** Semi-quantitative SIM areas in arbitrary units, no
internal standard, no calibration, no statistical test on volatiles; "triplicate" of unclear unit; one
muscle, one cooking programme, no pH. Good for the sign (thiamine yes, cysteine no, at 69 °C), not for any
rate constant or yield.

## What it does not give

- No absolute concentrations (no ng/g, no mM) of any volatile; no internal standard.
- No 2-furfurylthiol data.
- No temperature series, so nothing direct at 100 °C or 140 °C.
- No cysteine-alone series (cysteine always with 3 g sugar per 100 g brine) and no thiamine + cysteine
  combination.
- No pH, no free cysteine, ribose or IMP of the meat.
- No GC-O on the precursor hams; odour impact is only for the commercial ham.
- No thiamine assay for the 800 and 8000 mg doses.

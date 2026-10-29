# US 9,005,689 B2, Nestec, 2015 — EXTRACTION (multi-step vs single-step Maillard: glucose + glutamic acid + cysteine, each 0.4 mol/L, 0.5 mol/L phosphate pH 7, 100 °C; separate pre-reaction of cysteine at pH 9 (H2S) or pH 3 (cysteamine), 100 °C, 6 h; FFT, 2-acetyl-2-thiazoline and 5-acetyl-2,3-dihydro-1,4-thiazine yields per mol cysteine; pH series for H2S, cysteamine, methylglyoxal and 2,3-butanedione as bar charts)

### SULFUR-ROUTE pH SELECTIVITY, STATED AS A PATENT CLAIM: H2S from cysteine is favoured at alkaline pH (≈4× more at pH 9 than pH 3 after 5 h at 100 °C), cysteamine at acidic pH; splitting the reaction into pH-optimised steps raises FFT 10-fold (2.3 → 23 μg/mol Cys, 2 h) and 35-fold at 30 min, 2-AT ≈2× and ADHT ≈4×. BUT no analytical method, no replicate count and no error bars are printed; yields are "per mol cysteine" charged, at the μg/mol level (ppb conversion).

**Source on disk:** `data/articles/US9005689.pdf` (front page, four drawing sheets, specification cols 1-8). Text read with
`pdftotext -layout` (prose reliable; drawing sheets garbled in the text layer). Drawing sheets 1 and 3 (PDF pp. 2, 4) viewed
as 80-dpi renders to read axis titles, units and legends. **Figures are described, not digitised**; the only numbers taken
from figures are those printed as bar labels (Figs. 2 and 8), which the text repeats.

**Status of the evidence.** A patent, not peer reviewed; examples are written to support the multi-step claim. Inventors
include Hofmann and Schieberle, but **no analytical method (extraction, GC, standards, internal standards) is described
anywhere**, so the quantification basis is unknown. No replicates, no SD.

## 0. Identity

| field | value |
|---|---|
| Title | "Process for the preparation of a flavor active composition" |
| Patent number | **US 9,005,689 B2** (4 claims, 4 drawing sheets) |
| Assignee | Nestec S.A., Vevey (CH) |
| Inventors | Tomas Davidek (Correvon, CH); Imre Blank (Savigny, CH); Thomas Hofmann (Neufahrn, DE); Peter Schieberle (Zolling, DE) |
| Appl. No. | 13/390,177; PCT/EP2010/059683, PCT filed **Jul. 7, 2010**; §371 date Feb. 13, 2012 |
| PCT publication | WO2011/018280, Feb. 17, 2011 |
| Priority | EP 09167771(.6), **Aug. 13, 2009** |
| Prior publication | US 2012/0141643 A1, Jun. 7, 2012 |
| Date of patent | **Apr. 14, 2015** (term extended 367 days) |
| Cited literature | Wang et al. 2008 JAFC 56:7405 (HDMF via methylglyoxal); Kerler et al. 1997 Food Rev Int 13:553; Adams & De Kimpe 2006 Chem Rev 106:2299 |

## 1. Why it matters

Three things the sulfur lane can use as **directional** constraints, all at 100 °C in phosphate:
1. **H2S release from cysteine rises with pH** (pH 3 < 5 < 7 < 9, at 2 h and 5 h) — the text puts the pH 9 / pH 3 ratio at
   "about 4x" after 5 h.
2. **Cysteamine (cysteine decarboxylation) falls with pH** (pH 3 > 5 > 7 > 9 after 6 h) — the opposite selectivity.
3. **Sugar fragments (methylglyoxal, 2,3-butanedione) from glucose + glutamic acid peak near neutral pH** (30 min).
These are the pH levers behind the H2S-route vs cysteamine-route split for FFT vs 2-AT / ADHT.

## 2. Methods as they matter to a model

General ranges (col. 3–4; claims): steps at 60–180 °C (preferred 90–130 °C), 1 min–12 h (preferred 0.5–2 h), pH 2–9; water
5–99 % (preferred 60–90 %); catalysts "disodium hydrogen phosphate or citric acid". Claim 1 narrows each step to
**90–130 °C, 0.5–2 h, pH 2–9** — note that the examples' 6 h cysteine steps fall outside the claimed 0.5–2 h window.

### Example 4 — FFT (cols. 5–6, verbatim conditions)
- **Single step:** *"heating of all three precursors (each 0.4 mol/L) in phosphate buffer (0.5 mol/L; pH 7) at 100° C. for
  2 hours"* — glucose, glutamic acid, cysteine.
- **Multi step:** *"First, glucose was reacted with glutamic acid (each 0.8 mol/L) in phosphate buffer (0.5 mol/L; pH 7) at
  100° C. for 30 min (step 1: conversion of glucose to 2-furaldehyde). In another reaction a solution of cysteine (0.8 mol/L)
  was heated in phosphate buffer (0.5 mol/L. pH 9) at 100° C. for 6 hours (step 2: conversion of cysteine to hydrogen
  Sulphide). And finally, the reaction mixtures obtained from step 1 and step 2 were mixed together (ratio 1:1) and after
  adjusting pH value to pH 7 heated at 100° C. for 2 hours (step 3: generation of FFT)"*. After 1:1 mixing every precursor is
  nominally 0.4 mol/L, the same charge as the single step (my arithmetic).
- Fig. 3 legend repeats it: *"Step 1: Mix 1 - Glc/Glu, pH 7, 100°C / 30 min; Step 2: Mix 2 - Cys, pH 9, 100°C / 6 h; Step 3:
  Mix 1 + Mix 2, pH7, 100°C / 120 min"* vs *"Single step approach Step 1: Mix 1 - Glc/Glu/cys, pH 7, 100°C / 120 min"*.

### Example 5 — 2-AT and ADHT (col. 6)
- **Single step:** glucose + glutamic acid + cysteine, each 0.4 mol/L, pH 7 (phosphate 0.5 mol/L), **100 °C, 30 min**.
- **Multi step:** *"(i) reaction of glucose with glutamic acid (each 0.8 mol/L) in phosphate buffer (0.5 mol/L. pH 7) at
  100° C. for 30 min (step 1: generation of Sugar fragments); (ii) reaction of cysteine (0.8 mol/L) in phosphate buffer
  (0.5 mol/L; pH 3) at 100° C. for 6 hours (step 2: generation of cysteamine); and (iii) reaction of mixtures obtained in
  step 1 and 2 (mixed in ratio 1:1) at 100° C. for 30 min at pH 7 (step 3: generation of odorants)."*

### Example 6 — roasted-beef flavouring (col. 6), conditions only
- Multi step: xylose + glycine (each 1.5 mol/L) in **pyrophosphate 0.2 mol/L, pH 5.5**, 100 °C for 60 min (flavouring B) or
  120 min (flavouring C); then cysteine (1.5 mol/L) added, pH re-adjusted to 5.5, 100 °C, 120 min.
- Single step: xylose, glycine, cysteine (each 1.5 mol/L), same buffer and pH, 100 °C for 120 and 240 min (flavouring A).
- Sensory only, no numbers. ⚠️ The text contradicts itself: *"Flavouring B and C gave pleasant well balanced overall
  flavour ... On the other hand, the flavouring A and B developed rather unpleasant flavour"* — B is in both groups
  (probably "A" refers to the two single-step flavourings).

Examples 1–3 (roasty/popcorn, biscuit, chocolate) are recipes in words: no concentrations, temperatures or yields.

**Not printed anywhere:** analytical method and standards; replicate number; final or in-process pH; reaction volume and
vessel; headspace/atmosphere; come-up time; whether step-2 mixtures were analysed before mixing.

## 3. Every printed number

| quantity | single step | multi step | ratio (printed) | conditions | anchor |
|---|---|---|---|---|---|
| **FFT, 2 h** | **2.3 μg/mol Cys** | **23 μg/mol Cys** | "10 fold" | Ex. 4; step 3 = 2 h | col. 5 l. 17–19, col. 5 l. 43–44; Fig. 2 bar labels |
| **FFT, 30 min** | **0.4 μg/mol Cys** | **14 μg/mol Cys** | "35 fold" | Ex. 4; step 3 = 30 min | col. 5 l. 56–60; Fig. 2 |
| **2-acetyl-2-thiazoline (2-AT), 30 min** | **17 μg/mol Cys** | **37 μg/mol Cys** | "about twice" | Ex. 5 | col. 6 l. 31–43; Fig. 8 |
| **ADHT, 30 min** | **29 μg/mol Cys** | **114 μg/mol Cys** | "about 4 times" | Ex. 5 | col. 6 l. 31–44; Fig. 8 |
| **H2S, pH 9 vs pH 3** | — | — | "about 4x more ... after 5 hours at 100° C." | "glucose/cysteine", Fig. 1 | col. 5 l. 25–27 |

Arithmetic `[D]`: 23/2.3 = 10.0 ✔; 14/0.4 = 35 ✔; 37/17 = 2.18; 114/29 = 3.93 ✔. At 0.4 mol/L cysteine, 2.3 μg/mol is
0.9 μg/L FFT (single step, 2 h) and 23 μg/mol is 9.2 μg/L (multi step) — molar conversions of ≈2.0 × 10⁻⁸ and 2.0 × 10⁻⁷
(MW FFT 114.2). ADHT multi step 114 μg/mol ≈ 46 μg/L.

## 4. Figures (described, not digitised)

| fig. | content | axes / units | conditions | what it shows |
|---|---|---|---|---|
| 1 | H2S vs pH, "Glc/Cys system" | **H2S (mmol/mol Cys)**, 0–100; pH 3, 5, 7, 9 | 100 °C, 2 h (black) and 5 h (hatched) | monotonic rise with pH at both times; pH 9 / 5 h bar is the tallest, near the top of the axis. **Note the unit is mmol/mol — a percent-level H2S yield, five to six orders above the FFT yields.** |
| 2 | FFT, single vs multi step | FFT (μg/mol Cys), 0–40 | 100 °C 30 min, 2 h | bars labelled 0.4, 2.3, 14, 23 (= §3) |
| 3 | sensory spider plots | intensity 0–2.5; rubber, pungent, sulfurous, sweaty, spicy, meat-like, chicken-like, roasty, earthy, popcorn-like, caramel-like, sweet | Ex. 4 recipes (legend quoted in §2) | multi step higher in roasty, sweet, popcorn, caramel (text) |
| 4 | scheme: cysteamine + methylglyoxal → 2-AT | — | — | structural |
| 5 | scheme: cysteamine + 2,3-butanedione → ADHT | — | — | structural |
| 6 | methylglyoxal ("MeGlx") and 2,3-butanedione vs pH, "Glc/Glu system" | **c (mmol/mol Glc)**, 0–6; pH 3, 5, 7, 9 | "100°C, 30" (min) | MeGlx maximal at pH 7, 2,3-butanedione maximal at pH 5; both near zero at pH 3 |
| 7 | cysteamine vs pH, "Glc/Cys system" | **cysteamine (mmol/mol Cys)**, 0–400; pH 3, 5, 7, 9 | 100 °C, 6 h | monotonic fall with pH; the pH 3 bar is a large fraction of the axis (tens of percent of cysteine) |
| 8 | 2-AT and ADHT, single vs multi step | c (μg/mol Cys), 0–140 | Ex. 5 | bars labelled 17, 29, 37, 114 (= §3) |
| 9 | schematic of single vs multi step | — | — | text garbled in the text layer; not needed |

⚠️ **Figs. 1 and 7 are labelled "Glc/Cys system"** (glucose + cysteine), but the step-2 reactions they justify are
**cysteine alone** in buffer (Ex. 4: pH 9; Ex. 5: pH 3). The pH dependence was therefore measured with glucose present; the
step-2 yields of H2S and cysteamine actually fed into step 3 are not reported.

## 5. What this can and cannot anchor

- **Can:** the *direction* of pH selectivity — H2S ↑ with pH, cysteamine ↓ with pH, sugar C2/C3/C4 fragments maximal at
  pH 5–7 — at 100 °C in 0.5 mol/L phosphate; and the order of magnitude of FFT / 2-AT / ADHT yields per mol cysteine in a
  0.4 mol/L glucose/Glu/Cys pot (10⁻⁸–10⁻⁷ mol/mol). The 30 min → 2 h single-step FFT rise (0.4 → 2.3 μg/mol) is a two-point
  time course.
- **Cannot:** any kinetic constant (two times at most, no replicates, no method); absolute H2S or cysteamine values (figure
  only, and measured in Glc/Cys, not in the step-2 pots); hot pH (0.5 mol/L phosphate set at room temperature; nothing
  measured during heating).

## 6. Flags

1. **No analytical method** — quantification basis unknown (`quant_class: unknown`).
2. **No replicates, no SD.**
3. **Figs. 1 and 7 measured in Glc/Cys**, while the multi-step step 2 is cysteine alone.
4. **Example 6 text contradiction** (flavouring B both pleasant and unpleasant).
5. **6 h step-2 heating lies outside claim 1's 0.5–2 h window** — the examples are not within the claimed process.
6. Yields are per mol cysteine **charged**; in the multi-step arm cysteine was pre-heated 6 h, so the denominator is not the
   cysteine present at step 3.

## Full text (read 2026-10-29)

Read in full on 2026-10-29 from `data/articles/US9005689.pdf` (new dossier; sections 0–6 above are the reading). Every
number printed in the specification is in §3; the drawing sheets were read for axes, units and legends only.

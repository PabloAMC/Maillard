# Brewer & Vega 1995 — EXTRACTION (the same six aldehydes as the gelatin paper, dosed into raw lean ground beef mixed 1:1 with water, then COOKED to 70 C internal and sniffed at 45 C; triangle tests, ASTM E-1432-91, 10 panellists, group geometric means in ppm)

### THE COMPANION PAPER, AND THE ONE THAT MUST NOT BE TREATED LIKE ITS PARTNER: the dose is added BEFORE a 70 C cook, so these six numbers are quantities weighed into raw meat, not concentrations present at the moment of perception — and the paper's own molar column is wrong by a factor of 1 000 on four rows and 10 000 on two.

**Source on disk:** `data/articles/brewer1995.pdf` (4 pp., Journal of Food Science **60** (3) 1995,
592-595).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/brewer1995.txt`). The scan is a two-column journal page and the text layer
**interleaves the two columns**, so the prose must be read by following sentences rather than lines;
this was done and the whole article is accounted for. **Table 1 — the only table of thresholds —
came through clean** and every value in it reproduces what `k2_matrix_and_thresholds.md` sec. A.1
already carries. **Tables 2 to 6 (the descriptive-analysis means) are badly damaged**: descriptor
labels drop out, superscript letters merge into the numbers, and several "0" cells are rendered as
`0"`, `00:`, `,"i`, `$` or `;I:`. They are transcribed below with every uncertain cell marked
`[?]`, and no uncertain value is carried into section 4. Figure 1 (sensory descriptors for beef
containing added pentanal) is an image: **the pentanal descriptor profile is figure-only** — pentanal
is the one compound of the six with no descriptor table. There is no supplementary material. Repo
status before this dossier: Brewer 1995 is named in `src/kinetic_core/matrix_oav.py` in
`SEALED_OR_REFUSED_MATRICES["cooked_beef"]` (as declaration D.6 Module 7 **HOLD-OUT**, reclassified
`dose_added_pre_cook`), is excluded by name from the fit in
`src/kinetic_core/parameters_matrix.py`'s `ALPHA_BETA_UNSATURATION_OBSERVATIONS`, and is summarised
in `k2_matrix_and_thresholds.md` sec. A.1 — but has **no extraction dossier**.

## 0. Identity

| field | value |
|---|---|
| Title | "Detectable Odor Thresholds of Selected Lipid Oxidation Compounds in a Meat Model System" |
| Authors | M. Susan Brewer (corresponding) and Juan D. Vega — Division of Foods & Nutrition, University of Illinois at Urbana-Champaign, Urbana, IL 61801 |
| Venue | Journal of Food Science **60** (3), 1995, pp. 592-595. MS received 21 August 1994; revised 8 December 1994; accepted 15 January 1995 |
| DOI / article ID | none printed on the scan |
| Matrix | **lean beef top round, thawed, ground through a 0.50 cm plate, mixed 1:1 with distilled water**, 15 g sealed into a 25 mL amber vial, **cooked in a circulating water bath to 70 C internal temperature**, then **held at 45 C** for the sniff |
| Quantity measured | **DOT**, defined here as the concentration at which a panellist gives the correct triangle-test response **50 % of the time ABOVE CHANCE**, fitted per panellist by non-linear regression, then combined as a **group geometric mean**. ASTM E-1432-91 |
| Naming | the abstract and Methods both write "**t-2-heptanal**" where Table 1, the Results and the Conclusions all write "**heptanal**"; there is no 2-heptenal anywhere else in the paper and the compound was bought as an alkanal. Read as **heptanal**; the "t-2-" is a typographical error (Flags 5) |
| Companion paper | **Vega & Brewer 1994**, the same two authors reversed, the same six compounds, in a 3 % gelatin gel — `vega1994_extraction.md` (this batch). Cited here as ref. "Vega and Brewer, 1994" and compared against directly in the first sentence of Results |
| Second-hand content | Wick et al. 1967 thresholds for methional, phenylacetaldehyde and nonanal, described here as "**in beef**" (Vega 1994 calls the same three numbers "meat slurries") |

## 1. Why it matters

`src/kinetic_core/matrix_oav.py` refuses on every protein pot because `MATRIX_THRESHOLDS` holds
entries for only two matrices. **A meat model at ~100 g protein/L is exactly the kind of pot that
would fill the gap, and this paper is the reason the layer declines to use it.** The module's
`SEALED_OR_REFUSED_MATRICES["cooked_beef"]` entry states the reason in one sentence — the numbers
are "doses added to RAW beef before a 70 C cook, not concentrations present at the moment of
perception". **This dossier confirms that reading from the primary Methods paragraph and finds two
further defects the code does not yet record.**

| repository object | file | what this paper does to it |
|---|---|---|
| `SEALED_OR_REFUSED_MATRICES["cooked_beef"]` | `src/kinetic_core/matrix_oav.py` | this paper IS that entry; the `dose_added_pre_cook` reclassification is confirmed verbatim against the Methods |
| `MATRIX_THRESHOLDS` (matrix `cooked_beef`) | same file | **stays empty.** Six values exist and are re-typed in section 3; none is a threshold in this repository's sense |
| `ALPHA_BETA_UNSATURATION_OBSERVATIONS` | `src/kinetic_core/parameters_matrix.py` line ~490 | the beef t-2-hexenal/hexanal contrast (**1.34x**, mine) is the third observation of the unsaturation penalty and is **excluded from the fit** by the comment above that table; this dossier gives the number the exclusion applies to |
| `MATRIX_LOADING` | same file | there is **no** `cooked_beef` entry, and this paper supports building one only in part: the 1:1 meat:water dilution is stated, but **no proximate analysis of the beef is reported at all** (protein, fat, moisture and pH are every one of them unstated) |

**How it differs from the gelatin paper — the six differences that matter, in order of size.**

1. **A cook happens between dosing and sniffing.** Vega 1994 dosed a cooled gel and never heated it
   again; Brewer 1995 doses raw ground beef, seals it, and takes it to **70 C internal**, then holds
   it at **45 C**. Aldehydes are lost to the headspace of a sealed vial, to reaction with the meat's
   own protein, and to thermal chemistry, and **none of that is measured**. The paper's own
   discussion offers protein binding as the explanation for its high numbers, which is the same as
   saying that the number weighed in is not the number smelled.
2. **The matrix contains lipid, and the lipid is oxidising.** Gelatin was chosen in 1994 precisely
   because it is lipid-free; beef is a lipid-oxidation substrate that generates **these very six
   compounds** on its own. The blank is beef + water with nothing added, so the panel is detecting an
   **increment on an unmeasured background** of the same molecules. Neither the background level nor
   its variation over the six months of testing is reported.
3. **The method changed.** 1994: one sample against one control, ascending, **linear** fit, **75 %
   detection, uncorrected for chance**, 16 panellists. 1995: **triangle test**, ascending, **non-
   linear (PROC NLIN) fit per panellist**, **50 % correct ABOVE CHANCE**, ASTM E-1432-91, 10
   panellists, group **geometric** mean. These are different quantities. The 1995 criterion is the
   more defensible of the two — it is chance-corrected and it is a recognised standard — but it is
   **not the same criterion**, so the beef/gelatin ratios this paper computes in its own first
   paragraph of Results are cross-method even though the two studies share both authors.
4. **The panel changed.** 16 panellists (2 men, 14 women, mean age 23) became **10 panellists, all
   women, mean age 25**, trained in ten 1-h sessions rather than two.
5. **One temperature instead of four.** Everything here is at **45 C** serving temperature. The
   gelatin ladder's temperature axis has no counterpart.
6. **It adds something the gelatin paper does not have**: a full **descriptive analysis** on a
   0-5 category scale, for each of the six compounds, at three treatments (raw above threshold,
   cooked below threshold, cooked above threshold), with SEMs and significance letters — Tables 2-6
   and Fig. 1. That is a genuine, usable, *qualitative* result even though the thresholds are not
   usable quantitatively.

## 2. Methods as they matter to a model

- **The meat.** "Fresh lean beef **top round (35 kg)**" from the University of Illinois Meat Science
  Laboratory, double-wrapped in freezer paper, **frozen at -18 C in 454 g aliquots**. Thawed **4 C
  for 4 h**, ground through a **0.50 cm plate**, and **mixed 1:1 with distilled water** — either
  plain (the blank) or with the aldehyde stock solution — to predetermined concentrations. **No
  proximate composition is given**: not protein, not fat, not moisture, not pH, not the storage time
  before use. The ~100 g protein/L figure in `k2_matrix_and_thresholds.md` sec. A.1 is arithmetic on
  a conventional lean-beef protein content halved by the 1:1 dilution, **not** a number this paper
  prints.
- **The stocks.** GC-grade pentanal, hexanal, t-2-hexenal, "t-2-heptanal" (= heptanal, Flags 5),
  t-2-octenal and t,t-2,4-decadienal from Aldrich. Stock solutions in **distilled water at 22 C, at
  parts-per-million**, nitrogen-flushed, capped, stored at **4 C in amber glass bottles**. Same
  supplier and same stock-making protocol as the 1994 gelatin paper.
- **The cook — the load-bearing paragraph.** "Samples (15 g) were placed in **25-mL amber vials and
  sealed**. Samples were **cooked in a circulating water bath to 70 C (internal temperature)** and
  **maintained at 45 C** in a water bath for odor analysis. Samples were cooked **15-20 min before
  odor evaluation**." So: the compound is weighed into raw meat, the vial is closed, the whole thing
  is taken to 70 C, cooled to 45 C, and opened for sniffing 15-20 min later. **A 15 g sample in a
  25 mL vial leaves roughly 10 mL of headspace that equilibrates during the cook and is opened to the
  room at the sniff.** The concentration in the meat at the moment of perception is therefore lower
  than the dose by an unmeasured factor with at least three contributions (headspace partition at
  70 C, covalent and hydrophobic binding to the meat protein, thermal loss). **The paper measures
  none of them and makes no correction.**
- **Panel.** A **10-member experienced panel**, aged 18-35 (mean 25), **all women, all non-smokers**,
  trained in **ten 1-h sessions** with preliminary triangle testing to set the concentration ranges.
  Familiarisation used the added compounds in meat over **10 ppb to 100 ppm**, plus rancid oils, pork
  and beef as reference oxidised foods. The descriptor list was generated by group consensus from a
  composite of the individual panellists' terms. Practice concentrations for the descriptors:
  pentanal 50 ppm, hexanal 220 ppm, t-2-hexenal 120 ppm, heptanal 100 ppm, t-2-octenal 150 ppm,
  t,t-2,4-decadienal 110 ppm — **identical to the 1994 paper's training doses**.
- **Threshold test.** Samples presented **in increasing order of concentration**. Judges were given
  **6 sets of 3 samples** and asked to pick the odd sample **on odour alone**, pausing 10-20 s
  between samples. **10 replicates of each compound at each concentration, spread over 6 months.**
  Room at 22 C and **60 % relative humidity**, fans for positive air pressure.
- **How the threshold is computed, exactly.** "**PROC NLIN** (SAS, 1993) was used to analyze triangle
  test threshold data for **each aldehyde and panelist individually** by fitting nonlinear regression
  models using the least squares method (**ASTM, 1992** = E-1432-91). ... the percent of correct
  responses (**above chance**) was calculated at each concentration for each compound. ... The
  concentration at which panelists gave the correct response **50 % of the time (above chance)** was
  considered to be the detectable odor threshold **for a specific panelist**. **Geometric means** were
  calculated for group detectable odor thresholds." So Table 1's headline column is a **group
  geometric mean of ten individually fitted thresholds**, and the "lowest"/"highest" columns are the
  extreme individual panellists — a genuine dispersion measure that the gelatin paper does not
  provide.
- **`n/geometric mean = 180`** (Table 1 footnote a). The paper does not say how 180 is composed. Ten
  panellists times ten replicates is 100; ten panellists times eighteen concentration presentations
  is 180, as is ten replicates times eighteen. **Not asserted** (Flags 6).
- **Descriptive analysis.** A **5-point category scale: 0 = absent, 1 = slight, 2 = mild,
  3 = moderate, 4 = high, 5 = extremely high.** Three treatments per compound: **cooked at a
  BELOW-threshold concentration = 55 ppb; cooked at an ABOVE-threshold concentration = 55 ppm; and
  raw at 55 ppm**. Note that **55 ppb and 55 ppm are the same two doses for every compound** — they
  are not scaled to each compound's own threshold, so "above threshold" is 7x the DOT for t-2-hexenal
  and 239x for heptanal (mine, Flags 4). Three replications; LS means by GLM; means separated at
  p < 0.05 by probability of difference; **n/LS mean = 30**.
- **What is never measured.** No GC. No headspace analysis. No concentration verification of any
  kind — this paper has no instrumental analysis at all, unlike its 1994 companion, which at least
  ran headspace GC on the gel. No TBARS or peroxide value on the beef. No pH. No temperature other
  than the 45 C serving temperature.

## 3. Tables re-typed

### Table 1. "Detection odor threshold group geometric means for selected aldehydes in lean ground beef"

Column headings exactly as printed: `Geometric means` spanning `Threshold (M)` and `Threshold
(ppm)`, then `Lowest threshold value (ppm)` and `Highest threshold value (ppm)`. Footnote a:
"**n/geometric mean = 180**."

| | Threshold (M) | Threshold (ppm) | Lowest threshold value (ppm) | Highest threshold value (ppm) |
|---|---|---:|---:|---:|
| Pentanal | 3.0 x 10^-8 | **2.67** | 1.41 | 5.55 |
| Hexanal | 5.8 x 10^-8 | **5.87** | 2.37 | 36.81 |
| Heptanal | 0.2 x 10^-9 | **0.23** | 0.004 | 0.484 |
| t-2-Hexenal | 8.0 x 10^-8 | **7.87** | 2.04 | 10.62 |
| t-2-Octenal | 3.3 x 10^-8 | **4.20** | 1.12 | 26.30 |
| t,t-2,4-Decadienal | 0.3 x 10^-9 | **0.47** | 0.013 | 3.08 |

(The exponents in the (M) column are printed with OCR damage — `1O-8`, `IO-9`, `IO-*`, `10-g` — but
each is unambiguous from its neighbours and from the arithmetic in Flags 1: four rows read `x 10^-8`
and the heptanal and decadienal rows read `x 10^-9`. **The whole column is wrong regardless of how
the exponents are read** — see Flags 1.)

### Tables 2-6. Descriptive analysis, "Effect of cooking treatment on sensory characteristics of ground beef containing added <compound>"

Common headings and footnotes on all five tables: columns `Descriptor | Raw, above | Cook, below |
Cook, above | SEM`; footnote z "Category scale: 0 = absent, 5 = extremely high"; footnote y "n/LS
mean = 30"; footnote x "**Concentration of compound added = 55 ppm**" (this applies to both the
"Raw, above" and the "Cook, above" columns); footnote w "**Concentration of compound added =
55 ppb**" (the "Cook, below" column); footnote v "SEM = standard error of the mean"; footnote a,b,c
"Means in a row with different superscript letters are different (p < 0.05)".

**Superscript letters are merged into the digits by the scan in most cells.** Where a letter is
legible it is given; where it is not, it is omitted rather than guessed. `[?]` marks a value whose
digits themselves are uncertain; those values are not used anywhere else in this dossier.

**Table 2 — hexanal**

| Descriptor | Raw, above (55 ppm) | Cook, below (55 ppb) | Cook, above (55 ppm) | SEM |
|---|---:|---:|---:|---:|
| Putrid | 1.20 a | 0 c | 0.77 b | 0.11 |
| Sour | 1.33 b | 2.17 a | 1.87 a | 0.09 |
| Sweaty | 1.63 a | 0 | 2.00 a | 0.10 |
| Rancid | 3.03 a | 0 | 3.07 a | 0.16 |
| Animal | 1.10 b | 3.10 a | 1.50 b | 0.16 |
| Blood | 1.90 | 0 b | 0 | 0 |
| Fatty | 1.40 b | 3.00 a | 1.60 | 0.14 |
| Oily | 0.67 c | 3.13 a | 1.80 | 0.14 |
| Meaty | 1.00 c | 3.13 a | 1.67 b | 0.14 |
| Raw meat | 3.40 a | 0 | 0 b | 0.18 |
| Fishy | 1.70 a | 0 c | 1.27 b | 0.10 |
| Painty | 4.20 a | 0 [?] | 3.43 | 0.22 |
| Herbal | 3.43 a | 0 [?] | 3.43 a | 0.18 |

**Table 3 — heptanal**

| Descriptor | Raw, above (55 ppm) | Cook, below (55 ppb) | Cook, above (55 ppm) | SEM |
|---|---:|---:|---:|---:|
| Putrid | 1.34 a | [?] | 0.66 b | 0.10 |
| Sour | 1.20 c | [?] (the cell is unreadable: `*.*",“a`) | 1.77 b | 0.11 |
| Sweaty | 1.93 a | [?] | 2.00 a | 0.10 |
| Rancid | 3.47 a | [?] | 3.27 a | **0.8 [?]** (printed as `0.8`; every other SEM in the five tables is 0.03-0.22, so a lost digit is likely — not asserted) |
| Animal | 0.37 b | 3.77 a | 0.66 b | 0.19 |
| Blood | 0.61 a | 0 b | 0 | 0.07 |
| Fatty | 0.90 c | 3.90 a | 1.10 b | 0.17 |
| Oily | 0.90 c | 3.80 a | 1.40 b | 0.15 |
| Meaty | 0.90 c | 3.80 a | 1.40 b | 0.16 |
| Raw meat | 4.37 a | [?] | 0 b | 0.21 |
| Fishy | 1.80 a | [?] | 1.53 a | 0.12 |
| Painty | 4.23 a | [?] | 3.73 b | 0.22 |
| Herbal | 4.00 a | [?] | 4.00 a | 0.20 |

**Table 4 — t-2-hexenal**

| Descriptor | Raw, above (55 ppm) | Cook, below (55 ppb) | Cook, above (55 ppm) | SEM |
|---|---:|---:|---:|---:|
| Putrid (label lost in the scan; recovered from the row order shared by all five tables) | 1.44 a | 0 c | 0.88 b [?] | 0.12 |
| Sour | 1.20 | 2.30 a | 1.90 a | 0.08 |
| Sweaty | 1.40 b | [?] | 2.30 a | 0.12 |
| Rancid | 3.77 a | 0 | 3.53 a | 0.19 |
| Animal | 1.14 b | 2.53 a | 1.27 b | 0.13 |
| Blood | 0.11 | 0 | 0 | 0.03 |
| Fatty | 0.71 c | 2.60 a | 1.50 b | 0.14 |
| Oily | 1.03 c | 3.07 a | 1.60 b | 0.14 |
| Meaty | 0.83 c | 3.07 a | 1.?? b [?] | 0.14 |
| Raw meat | 3.87 a | 0 | 0 | 0.20 |
| Fishy | 1.97 | 0 c | 1.?? a [?] | 0.11 |
| Painty | 2.87 a | 0 c | 2.43 b | 0.20 |
| Herbal | 4.00 | 0 | 4.00 | 0.20 |

**Table 5 — t-2-octenal**

| Descriptor | Raw, above (55 ppm) | Cook, below (55 ppb) | Cook, above (55 ppm) | SEM |
|---|---:|---:|---:|---:|
| Putrid | 1.20 a | 0 c | 0.77 b | 0.11 |
| Sour | 1.13 c | 3.00 a | 1.60 b | 0.11 |
| Sweaty | 1.63 b | [?] | 2.00 a | 0.10 |
| Rancid | 3.87 a | 0 | 3.87 a | 0.21 |
| Animal | 0.54 c | 3.73 a | 1.0? b [?] | 0.20 |
| Blood | 0.14 a | 0 b | 0 | 0.04 |
| Fatty | 0.21 c | 4.20 a | 0.81 b | 0.21 |
| Oily | 1.16 b | 4.26 a | 1.13 b | 0.17 |
| Meaty | 1.13 b | 4.27 a | 1.00 b | 0.17 |
| Raw meat | 4.00 a | 0 b | 0 | 0.20 |
| Fishy | 1.40 a | 0 c | 1.2? b [?] | 0.10 |
| Painty | 3.93 a | [?] | 4.00 a | 0.21 |
| Herbal | 4.00 a | [?] | 4.00 a | 0.20 |

**Table 6 — t,t-2,4-decadienal**

| Descriptor | Raw, above (55 ppm) | Cook, below (55 ppb) | Cook, above (55 ppm) | SEM |
|---|---:|---:|---:|---:|
| Putrid | 1.04 a | 0 c | 0.?? b [?] | 0.10 |
| Sour | 1.?? c [?] | 3.00 a | ?? b [?] | 0.10 |
| Sweaty | 1.30 b | [?] | 2.23 a | 0.12 |
| Rancid | 2.87 | [?] | 2.80 a | 0.19 |
| Animal | 0.67 c | 3.67 a | 1.34 b | 0.18 |
| Blood | 1.11 | 0 b | 0 b | 0.09 |
| Fatty | 1.70 b | 3.80 a | 1.70 b | 0.17 |
| Oily | 0.90 c | 3.90 a | 1.80 b | 0.16 |
| Meaty | 0.60 c | 3.90 a | 1.80 b | 0.17 |
| Raw meat | 4.43 a | 0 | 0 | 0.22 |
| Fishy | 2.10 a | 0 c | 2.27 a | 0.15 |
| Painty | 3.00 a | 0 c | 4.00 b | 0.22 |
| Herbal | 4.60 a | 0 c | 2.20 b | 0.20 |

**Pentanal has no descriptor table.** Its profile is Fig. 1 ("Sensory descriptors for beef containing
added pentanal. Scale: 0 = none, 5 = intense", three series, thirteen descriptors) and is therefore
**figure_only**.

### Numbers printed in the running text

| quantity | value | where | whose measurement |
|---|---|---|---|
| DOT range in cooked beef | 0.23 ppm (heptanal) to 7.87 ppm (t-2-hexenal) | Abstract, Results, Conclusions | this paper |
| rank order of cooked-beef DOT | heptanal > t,t-2,4-decadienal > pentanal > t-2-octenal > hexanal > t-2-hexenal (printed as an *increasing* list) | Conclusions | this paper |
| beef vs gelatin, claimed | "**three orders of magnitude** for pentanal, hexanal, t-2-hexenal, t-2-octenal, and **two orders of magnitude** for heptanal and t,t-2,4-decadienal" | Results, first sentence | this paper — **and it is refuted by its own tables**, Flags 2 |
| "in beef": methional | 6.1 ppm | Results | **Wick et al. 1967**, quoted |
| "in beef": phenylacetaldehyde | 0.94 ppm | Results | Wick et al. 1967, quoted |
| "in beef": nonanal | 7.6 ppm | Results | Wick et al. 1967, quoted |
| descriptive-analysis doses | 55 ppb (below threshold) and 55 ppm (above threshold) | Methods and every table footnote | this paper |
| predominant-descriptor cutoff | "more than 3.0, scale = 1-5" | Results | this paper (note the scale is defined as 0-5 in Methods and as 1-5 here) |

### Arithmetic on the printed thresholds (all mine)

**1. Beef / gelatin, compound by compound.** Against the companion paper's 22 C gelatin column
(`vega1994_extraction.md` Table 1): pentanal 2670/41 = **65x**; hexanal 5870/58 = **101x**; heptanal
230/79 = **2.9x**; t-2-hexenal 7870/109 = **72x**; t-2-octenal 4200/109 = **39x**; decadienal
470/64 = **7.3x**. **The whole set spans 2.9x to 101x — between half an order and two orders of
magnitude.** Choosing a different gelatin temperature does not rescue the paper's claim: against the
37 C gelatin column the ratios are 79 / 173 / 3.7 / 100 / 40 / 5.3, and against 60 C they are
121 / 154 / 4.6 / 131 / 52 / 7.3. **No choice of comparison temperature produces "three orders of
magnitude" for any compound**, and the two compounds the paper says differ by "two orders" (heptanal,
decadienal) are in fact the two that differ *least*, by 2.9x and 7.3x.

**2. Beef / water.** Against the Guadagni values the companion paper quotes: pentanal 2670/12 =
**223x**; hexanal 5870/4.5 = **1 304x**; heptanal 230/3 = **77x**; t-2-hexenal 7870/3 = **2 623x**;
t-2-octenal 4200/3 = **1 400x**; decadienal 470/0.07 = **6 714x**. These reproduce
`k2_matrix_and_thresholds.md` sec. A.1 exactly. They are the largest matrix/water ratios in the
corpus **and the least interpretable**, because the numerator is a pre-cook dose and the denominator
is a 1960s forced-choice water threshold.

**3. The unsaturation contrast in beef.** t-2-hexenal / hexanal = 7.87/5.87 = **1.34x**, against
**1.88x** in the gel (same authors, same compounds) and **6.88/1.39 = 4.95x** in Meynier's skim milk
by headspace partition. Expressed the way `parameters_matrix.py` expresses it — beef/water on
t-2-hexenal over beef/water on hexanal — it is 2623/1304 = **2.01x**, which is the "2.0x beef
observation" the code's comment says is excluded from the fit. **Confirmed: the number is 2.0x and
the exclusion is on provenance, not on the value being an outlier** — 2.0x sits neatly between the
gelatin 2.81x and nothing else, and would have *strengthened* the fit had it been eligible.

**4. Inter-panellist spread, which is the real headline of Table 1.** Highest/lowest per compound:
pentanal 5.55/1.41 = **3.9x**; hexanal 36.81/2.37 = **15.5x**; heptanal 0.484/0.004 = **121x**;
t-2-hexenal 10.62/2.04 = **5.2x**; t-2-octenal 26.30/1.12 = **23.5x**; decadienal 3.08/0.013 =
**237x**. **Two compounds show more than two orders of magnitude between the least and the most
sensitive of ten panellists.** The heptanal group geometric mean of 0.23 ppm sits on a spread from
4 ppb to 484 ppb. Any use of these six numbers must carry these ranges; the group mean alone is
misleading by a wide margin, and — unlike the 1994 paper, which reports no threshold dispersion at
all — **this paper gives you the honest picture and it is not reassuring.**

**5. The descriptive doses relative to each compound's own threshold.** The "above threshold" dose is
55 ppm for every compound, so it is 55/2.67 = **21x** the DOT for pentanal, 55/5.87 = **9.4x** for
hexanal, 55/0.23 = **239x** for heptanal, 55/7.87 = **7.0x** for t-2-hexenal, 55/4.20 = **13x** for
t-2-octenal and 55/0.47 = **117x** for decadienal. The "below threshold" dose is 55 ppb, i.e.
55/2670 = **0.021x** the DOT for pentanal down to 55/230 = **0.24x** for heptanal. **The five
descriptor tables therefore compare compounds at wildly different multiples of their own
thresholds**, which limits any cross-compound reading of Tables 2-6 (Flags 4).

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Same as for the companion paper:
`hexanal` and `heptanal` are keyed; `e_2_octenal` is keyed and covers t-2-octenal if the naming is
confirmed; **`pentanal`, `t_2_hexenal` and `tt_2_4_decadienal` are absent**. Every row below shares:
**lean beef top round, thawed, ground through a 0.50 cm plate, mixed 1:1 with distilled water,
composition unstated, 15 g sealed in a 25 mL amber vial, dosed RAW, cooked in a water bath to 70 C
internal, held and sniffed at 45 C, 15-20 min after cooking; triangle test, ascending, ASTM
E-1432-91, per-panellist non-linear fit to 50 % correct above chance, group geometric mean of 10
panellists (all women, 18-35 y, non-smokers), 10 replicates per concentration over 6 months, room
22 C / 60 % RH.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| pentanal, dose added to raw beef giving 50 %-above-chance detection after cooking | 2.67 (range 1.41-5.55) | ppm (mg/kg) | 45 C serving, cooked to 70 C after dosing | Table 1 p. 592 | **threshold** — but **`dose_added_pre_cook`**: NOT a concentration at perception (Flags 3) |
| hexanal, same | 5.87 (range 2.37-36.81) | ppm | as above | Table 1 | **threshold**, `dose_added_pre_cook` |
| heptanal, same | 0.23 (range 0.004-0.484) | ppm | as above | Table 1 | **threshold**, `dose_added_pre_cook` |
| t-2-hexenal, same | 7.87 (range 2.04-10.62) | ppm | as above | Table 1 | **threshold**, `dose_added_pre_cook` |
| t-2-octenal, same | 4.20 (range 1.12-26.30) | ppm | as above | Table 1 | **threshold**, `dose_added_pre_cook` |
| t,t-2,4-decadienal, same | 0.47 (range 0.013-3.08) | ppm | as above | Table 1 | **threshold**, `dose_added_pre_cook` |
| inter-panellist spread (highest/lowest) | 3.9 / 15.5 / 121 / 5.2 / 23.5 / 237 | x | pentanal / hexanal / heptanal / t-2-hexenal / t-2-octenal / decadienal, 10 panellists | derived from Table 1 (mine) | within_study_ratio — **the most reliable quantity in the paper** |
| unsaturation contrast in beef (t-2-hexenal / hexanal, raw ratio) | 1.34 | x | 45 C | derived (mine) | within_study_ratio |
| unsaturation penalty, beef, in the code's ratio-of-ratios form | 2.01 | x | 45 C; both legs divide by a Guadagni water value | derived (mine) | within_study_ratio — **HOLD-OUT; excluded from the `ALPHA_BETA_UNSATURATION_OBSERVATIONS` fit by name** |
| beef/gelatin ratio at the companion paper's 22 C | 65 / 101 / 2.9 / 72 / 39 / 7.3 | x | same authors, DIFFERENT method and panel | derived (mine) | within_study_ratio — **cross-method despite the shared authorship** (Flags 2) |
| beef/water ratio | 223 / 1 304 / 77 / 2 623 / 1 400 / 6 714 | x | vs Guadagni 1963/72 | derived (mine) | within_study_ratio — cross-study, cross-method, **and the numerator is a pre-cook dose** |
| the "Threshold (M)" molar column | 3.0e-8 / 5.8e-8 / 0.2e-9 / 8.0e-8 / 3.3e-8 / 0.3e-9 | M | — | Table 1 | **REFUSE. Wrong by 10^3 on four rows and 10^4 on two** (Flags 1) |
| beef 1:1 with water | 1:1 | w/w | ground beef : distilled water | Methods p. 592 | level_only (**the beef's own composition is not stated**) |
| cook endpoint | 70 | C internal | circulating water bath, sealed 25 mL vial | Methods | level_only |
| serving temperature | 45 | C | held in a water bath | Methods | level_only |
| panel size | 10 | panellists | all women, 18-35 y, 10 replicates over 6 months | Methods | level_only |
| descriptor intensities, five compounds x three treatments x thirteen descriptors, with SEM and p<0.05 letters | see section 3 | 0-5 category scale | 55 ppb cooked / 55 ppm cooked / 55 ppm raw | Tables 2-6 pp. 593-594 | level_only (**and heavily OCR-damaged**; the uncertain cells are marked and not used) |
| pentanal descriptor profile | — | 0-5 scale | three treatments | Fig. 1 | **figure_only** |
| Wick 1967 "in beef": methional / phenylacetaldehyde / nonanal | 6.1 / 0.94 / 7.6 | ppm | matrix labelled "beef" here and "meat slurries" by Vega 1994 | Results p. 593 | **threshold** — second-hand, and the matrix label is contested (Flags 7) |

### What could and could not be done with these

**(a) They cannot enter `MATRIX_THRESHOLDS`, and the reason is chemical, not bureaucratic.** The
module's contract is that a `ThresholdRecord` holds "the compound's concentration in the stated
matrix at which the panel detects it". These six numbers are the mass of aldehyde stirred into raw
meat that, after being sealed in a vial, taken to 70 C and cooled to 45 C, ends up detectable. The
gap between the two is the very quantity — protein binding of aldehydes at cooking temperature —
that `src/kinetic_core/matrix_sites.py` exists to model, and **this paper's own Discussion says so**:
it attributes its high numbers to aldehydes reacting "with side chains of lysine, arginine,
methionine, and tyrosine to form Schiff bases", to hydrophobic interaction increasing as the protein
unfolds during heating, and to bound aldehydes behaving "as non-volatile compounds which reduce their
concentrations in the headspace". Ingesting the dose as a threshold would silently absorb an
unmeasured binding term into the threshold table, in the one place the layer is built to keep them
apart.

**(b) But the paper is a qualitative witness for the binding layer.** Its three cited mechanisms
(Solms 1973; Arai 1970, 1980; Sydow 1975) are the same literature the `BINDING_CLASSES` brackets rest
on, and its observation that "hexanal and t-2-hexenal, which exhibited higher DOT values in our
study, have been shown to bind more readily to protein than do the other aldehydes evaluated" is a
sign-level corroboration of the `unsaturated_aldehyde_amine` bracket sitting **above** the
`saturated_aldehyde_amine` bracket (5.3-7.9e-5 vs 6.0e-6-2.5e-5 M^-1 s^-1). It is corroboration of a
sign, not of a magnitude, and it is not independent of the sources already cited.

**(c) The one number here the repository could newly use is the inter-panellist spread.** Ten
panellists on the same six compounds in the same matrix on the same day give 3.9x to 237x between
the least and most sensitive. `matrix_oav.py` carries `HS_SPME_SAME_SAMPLE_DISPERSION` (10-23x) and
`K_AW_UNCERTAINTY_DECADES` (±0.5) as its two reliability bands and has **no panel-dispersion band at
all**. This paper measures one, in a protein matrix, and it is larger than either.

**(d) What cannot be transported.** Nothing from this paper reaches a plant-protein pot. Beef is
myofibrillar protein with heme iron and a live lipid-oxidation background; the composition is not
even reported. The matrix has **no** entry in `data/species/protein_matrices.yml` and nothing here
would let one be written — no thiol assay, no lysine assay, no protein assay.

## 5. Flags

1. **The "Threshold (M)" column is arithmetically wrong, and I reproduced the error independently.**
   Converting the printed ppm values with the compounds' molar masses: pentanal 2.67 mg/L / 86.13 =
   **3.1e-5 M** against a printed 3.0e-8; hexanal 5.87/100.16 = **5.86e-5** against 5.8e-8; heptanal
   0.23/114.19 = **2.0e-6** against 0.2e-9 (= 2e-10); t-2-hexenal 7.87/98.14 = **8.02e-5** against
   8.0e-8; t-2-octenal 4.20/126.20 = **3.33e-5** against 3.3e-8; decadienal 0.47/152.23 = **3.09e-6**
   against 0.3e-9 (= 3e-10). **The mantissas are right in all six rows and the exponents are wrong by
   10^3 in four and 10^4 in two.** The most likely mechanism is a mg/L-to-g/L slip compounded, on the
   two sub-ppm rows, by an additional decimal shift when the mantissa was normalised to a leading
   "0.". `k2_matrix_and_thresholds.md` already carries this warning; it is now verified from the
   printed table. **Refuse the (M) column outright and derive molarity from the ppm column if it is
   ever needed.**
2. **The paper's headline comparison against its own companion is refuted by its own tables, and this
   is a named laundering hazard.** The first sentence of Results claims "three orders of magnitude"
   against the gelatin system for four compounds and "two orders" for the other two. The true ratios
   are **2.9x to 101x** against the 22 C gelatin column, and no gelatin temperature gives a
   different answer (§3 arithmetic 1). A reader who quotes the Results sentence rather than the
   tables will overstate the matrix effect by a factor of 10 to 300. **The claim also has the
   direction of the two extremes backwards**: heptanal and decadienal, which the paper calls the
   *smaller* effect at "two orders", are the two smallest at 2.9x and 7.3x, while the "three orders"
   compounds top out at 101x.
3. **`dose_added_pre_cook` — the reclassification is confirmed and should be considered the paper's
   defining property.** The dose goes into raw meat; a 70 C cook happens; the sniff is at 45 C.
   Nothing is measured after the dose. There are at least three unquantified sinks (headspace of the
   sealed vial at 70 C, covalent and hydrophobic binding to meat protein, thermal chemistry) and the
   paper's own Discussion argues that the second is large. **The six numbers are upper bounds on the
   perceptual threshold, of unknown tightness, and they are not comparable to any threshold in any
   other paper in this corpus.**
4. **The descriptive analysis uses one dose pair for all six compounds, so Tables 2-6 are not
   comparable across compounds.** 55 ppm is **7.0x t-2-hexenal's** threshold
   but **239x heptanal's**, and 55 ppb is 0.021x pentanal's threshold but 0.24x heptanal's. So the
   "above threshold" column is a mild suprathreshold dose for some compounds and a massive one for
   others, and the "below threshold" column is nearly detectable for heptanal and utterly
   undetectable for pentanal. Read each table down its own column, never across tables.
5. **Compound naming is inconsistent within the paper.** The Abstract and the Materials paragraph
   both list "**t-2-heptanal**"; Table 1, the Results and the Conclusions all say "**heptanal**".
   There is no 2-heptenal elsewhere in the paper, the rank orders in Results and Conclusions treat
   the compound as the third saturated aldehyde in the C5-C6-C7 series, and the 1994 companion used
   heptanal. **Read as heptanal.** Also note that the 1994 companion's own Materials paragraph omits
   hexanal from its purchase list while measuring it — the pair of papers is casual about naming.
6. **`n/geometric mean = 180` is unexplained.** Ten panellists x ten replicates is 100. The paper
   says judges got "6 sets of 3 samples" and performed "10 replicates of samples containing each
   compound at each concentration over 6 mos", but never states the number of concentrations per
   compound. 180 is consistent with 10 panellists x 18 presentations or with 10 replicates x 18, and
   **neither is asserted**. Without the concentration ladder, the per-panellist fits cannot be
   reconstructed or re-analysed.
7. **The Wick 1967 numbers are second-hand here and are labelled differently by the two companion
   papers.** Brewer 1995 calls them thresholds "**in beef**"; Vega 1994 calls the identical three
   numbers "in **meat slurries**". Neither paper measured them. Attributing them to Brewer or to Vega
   would be laundering, and the matrix itself is unresolved between the two labels. Carry them as
   Wick 1967 or not at all.
8. **The blank is not blank.** Beef oxidises to give pentanal, hexanal, heptanal, t-2-hexenal,
   t-2-octenal and t,t-2,4-decadienal — the paper's own Introduction says exactly this, citing
   St. Angelo 1987 for their presence "in parts per million" in reheated beef. The panel is
   therefore detecting an increment on a background of the same molecules, over six months of
   testing, and **the background is never measured**. If the background is at the ppm level as the
   Introduction says, it is comparable to the thresholds being measured, which would inflate them.
9. **No instrumental analysis whatsoever.** Unlike the 1994 companion (which at least ran headspace
   GC on the gel), this paper has no GC, no headspace measurement, no TBARS, no peroxide value, no
   proximate analysis and no pH. Every number in it is sensory.
10. **What this paper does not contain**: any composition of the beef (protein, fat, moisture, pH,
    ash, heme, iron); any concentration measured after cooking; any temperature series; any
    same-method water or gelatin arm; any partition coefficient; any binding constant; any
    thiol or amine assay; any lipid-oxidation index; any supplementary material.
11. **What to request from the authors or a follow-up**: (i) the concentration ladder used for each
    compound and the raw triangle-test counts, which would let the per-panellist fits be redone and
    would resolve the `n = 180`; (ii) a headspace or extraction measurement of the six compounds in
    the cooked vial at 45 C, which would convert the six doses into six actual thresholds and is the
    single experiment that would make this paper usable; (iii) a proximate analysis of the beef;
    (iv) the descriptor tables in a legible form.
12. **Registry gaps against `data/keys/compounds.yml`**: identical to the companion paper —
    `hexanal`, `heptanal` and (as `e_2_octenal`) t-2-octenal are keyed; **`pentanal`, `t_2_hexenal`
    and `tt_2_4_decadienal` are not**. In addition, the compounds this paper's Discussion names as
    the binding partners on the protein side (**lysine, arginine, methionine, tyrosine**) are not
    registry ids either, and only the lysine epsilon-amine pool has a density on file in
    `data/species/protein_matrices.yml`; **arginine is counted for beta-lactoglobulin (3 per
    monomer) but is not charged as a binding pool by `matrix_sites.py`**, which is a real modelling
    gap this paper points at.

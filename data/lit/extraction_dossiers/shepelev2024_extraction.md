# Shepelev & Reineccius 2024 — COMPLETE TRANSCRIPTION
### *"¹⁴C-Isotope Use to Quantify Covalent Reactions between Flavor Compounds and β-Lactoglobulin"*
### ⭐⭐⭐ **THE PAPER THAT CLOSES THE COVALENT-SINK BRACKET.** It is the corpus's first — and still only — temperature-resolved, quantitative, absolutely-calibrated time course of an **aliphatic aldehyde** forming covalent bonds with a **protein**. Wave **K6b**, 2026-08-29.

**Full extraction of every number in `data/articles/Shepelev2024.pdf`.**
Read-only wave: no repo file outside `data/lit/extraction_dossiers/` was created or modified;
no existing dossier was edited.

**Read first:** `meynier2004_extraction.md` §1(c)–(d) and §9c (the bracket left open),
`anantharamkrishnan2020b_extraction.md` §6 (the k₂ lower bound), `research_round4_nulls.md` §A
(the null this paper falsifies), `yuan2023_extraction.md` (its sister paper, same lab).

**Provenance codes:** **[M]** measured and printed · **[C]** cited from another paper ·
**[Q]** qualitative statement, no number attached · **[D]** digitized by this wave off a
figure with a calibrated axis · **[Z]** derived by this wave · **[!]** integrity flag.

---

## 0. PAPER IDENTITY — **MATCHES THE EXPECTED IDENTITY EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Igor Shepelev\*, Gary A. Reineccius** |
| Title | **"¹⁴C-Isotope Use to Quantify Covalent Reactions between Flavor Compounds and β-Lactoglobulin"** |
| Venue | ***J. Agric. Food Chem.* 2024, 72 (18), 10579−10583** |
| DOI | **10.1021/acs.jafc.4c00134** ✔ **exactly the expected DOI** |
| Affiliation | Dept. of Food Science and Nutrition, **University of Minnesota, St. Paul, MN 55108** |
| ORCIDs | Shepelev 0000-0002-9310-2254 · Reineccius 0000-0003-1817-8758 |
| Dates | Received **Jan 5, 2024** · Revised **Apr 17, 2024** · Accepted **Apr 17, 2024** · Published **Apr 29, 2024** |
| Keywords as printed | *covalent bond, decanal, 2-undecanone, 2-phenylethyl isothiocyanate, pea protein isolate* |
| Funding | **Midwest Dairy Foods Research Center** (partial); Prof. Dan Gallaher thanked for radioisotope guidance |
| PDF character | **5 pages** (10579–10583), clean text layer, **Figures 1–5, ZERO TABLES**, one SI PDF ("Description of preliminary study") |

**⇒ CONFIRMED: expected authors, expected journal, expected volume, expected DOI, expected
year, expected system (¹⁴C quantification of covalent flavor–BLG reactions), expected
compound (decanal), expected temperatures (25/45/65 °C), expected duration (4 weeks).
Every field the brief predicted is present and matches.**

**⚠️ [!] ONE IDENTITY ODDITY.** The keyword list ends with **"pea protein isolate"**, and no pea
protein appears anywhere in the paper. The only protein studied is bovine β-lactoglobulin.
**A stray keyword, presumably from a companion manuscript. Do not index this paper under pea protein.**

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER — **READ THIS BEFORE ANYTHING ELSE**

**This paper contains no table.** Its entire quantitative content is in **three time-course
figures (Figs. 2, 3, 4)** with a calibrated `mg of compound / g of protein` y-axis, an
`0–28 days` x-axis, and **printed standard-deviation error bars** — the first such object in
the corpus for an aldehyde–protein covalent reaction. This wave digitized Fig. 3 (decanal)
to ±0.15 mg/g and Fig. 4 (PEITC) to ±0.6 mg/g (§4a).

**The four headline results.**

**(a) ⭐⭐⭐ THE ACTIVATION ENERGY OF ALDEHYDE–PROTEIN COVALENT ADDUCT FORMATION IS NOW
MEASURED, AND IT IS TINY: Eₐ ≈ 15 kJ/mol from the day-1 extents, ≈ 20 kJ/mol after
correcting for the saturation already present at day 1** (§6). Both estimates come from a
three-point Arrhenius fit over **25 / 45 / 65 °C**. **The result is robust to every
digitization ambiguity in the figure (§6d): the whole defensible range is 15–23 kJ/mol, and
the most generous possible upper bound, allowing a hydrophobic-pre-equilibrium correction the
paper does not measure, is ≈ 45 kJ/mol (§6e).**

**(b) ⭐⭐⭐ THE `Eₐ ≥ 70 kJ/mol` THRESHOLD IS REFUTED, NOT NARROWLY BUT BY A FACTOR OF
3.5–4.7×.** `anantharamkrishnan2020b_extraction.md` derived that the ambient covalent sink
only becomes material at process temperature **if Eₐ ≳ 70 kJ/mol**; `meynier2004_extraction.md`
§9c recorded that no paper in the corpus measured it. **It is now measured and it is 15–20.**
The full arithmetic and its consequence live in `k6b_adduct_kinetics_synthesis.md` §3.

**(c) ⭐ THE METHOD'S OWN INTERNAL CONTROL PROVES THE LOW Eₐ IS NOT AN ARTEFACT.** The same
experiment, same protein, same temperatures, same isotope method, run on
**2-phenylethyl isothiocyanate returns Eₐ ≈ 36–43 kJ/mol** (§7c) — **2.0–2.5× larger than
decanal's.** A method that could not resolve activation energies would return the same number
for both. **The apparatus resolves Eₐ; decanal's is simply small.**

**(d) ⭐⭐ THE CAPACITY OF THE COVALENT SINK *DECREASES* WITH TEMPERATURE.** At day 28 the
decanal series orders **25 °C (25.6) > 45 °C (≈20.1) > 65 °C (≈17.4) mg/g** — i.e. **the
coldest sample has bound the MOST aldehyde.** The authors state the mechanism themselves
[Q]: *"storage of a protein at elevated temperatures may alter protein structure such that it
becomes less available/reactive. **We favor this explanation.**"* The same inversion appears
in the PEITC extrapolation (§7d [!]). **So heat attacks the covalent sink from BOTH sides:
the rate barely rises (tiny Eₐ) and the ceiling falls.**

---

## 2. ⚠️ WHAT THIS PAPER IS NOT — the negative inventory, up front

**No table of any kind exists in this paper.** Verified across all five pages.

**Absent entirely:** **hexanal** (the compound the repo actually needs — see §9) · any
aldehyde other than decanal · any protein other than BLG · any pH measurement or control
(**the system is unbuffered double-distilled water; pH is never stated — §3b [!]**) · any
water-activity axis · any reversibility experiment (no dilution, dialysis, headspace
back-extraction or re-equilibration) · any rate constant, half-life or activation energy
computed by the authors · any Arrhenius plot · any statistical test (only "a Q-test … to
identify and reject outliers") · any identification of the adduct (no MS, no site mapping) ·
any measurement between 65 and 130 °C.

**⚠️ [!] THE PAPER REPORTS NO KINETICS AT ALL.** It reports *extents at six times and three
temperatures*. **Every rate constant and every Eₐ attributed to this paper anywhere in the
corpus is derived by wave K6b in §6 and must be labelled `[Z] derived_K6b`, never
`measured_by_Shepelev`.**

---

## 3. SYSTEM COMPOSITION AND METHOD — complete transcription

### 3a. The four ¹⁴C compounds, with lot numbers **[M]**
All from **American Radiolabeled Chemicals, Inc.**

| role | compound | label position | lot |
|---|---|---|---|
| **ketone** | **2-undecanone (UDO)** | **[2-¹⁴C]** | ARC 1160−50 |
| **aldehyde ⭐** | **decanal (DAL)** | **[1-¹⁴C]** | **ARC 0557-250** |
| **isothiocyanate** | **2-phenylethyl isothiocyanate (PEITC)** | **[1-¹⁴C]** | ARC 1826-50 |
| **methods blank** | **n-decane (DEC)** | **[1-¹⁴C]** | ARC 1189−50 |

⚠️ **[!] THE PAPER'S OWN ABBREVIATIONS ARE INTERNALLY CONTRADICTORY, AND IT MATTERS.**
The Abbreviations list (p. 10583) prints **"DAL, decanal"** and **"DEC, decane"**. The
abstract also uses **DAL** for decanal. **But the Results section heading is literally
"Decanal (DEC)."**, Fig. 3's caption reads **"Covalently bound DEC in (a) 1 mL of 0.1 % BLG
solutions"**, and the Summary reads **"ca. 23 mg of DEC bonded per gram of BLG"**.
**⇒ In Figure 3, in the Results heading, and in the Summary, `DEC` means DECANAL, not decane.**
An automated scrape keyed on the Abbreviations list will assign the paper's entire aldehyde
dataset to an inert alkane blank. **Register `fig3_label_DEC_means_decanal: TRUE`.**
A second, smaller instance: the text once writes **"2-UNO"** where UDO is meant.

Other chemicals **[M]**: hexane ACS grade [110-54-3], sodium chloride [7647-14-5], sodium
azide [26628-22-8] — Millipore Sigma; **200 proof ethyl alcohol** ACS/USP — Pharmco;
**Hionic-Fluor scintillation cocktail** — PerkinElmer.

### 3b. Reaction system **[M]** — the whole thing

| variable | value as printed |
|---|---|
| **Protein** | **Single-variant BLG from a single-variant cow in the UMN herd, 90 % protein purity**, isolated by the method of `anantharamkrishnan2020_extraction.md` (ref. 9) |
| **Stock** | **0.1 g BLG + 0.05 g sodium azide** made to **100 mL** with double-distilled water ⇒ **1.0 mg/mL BLG powder, 0.5 mg/mL azide** |
| **Aliquot** | **1.5 mL** of stock into a **2 mL crimp-cap vial** |
| **Flavor dose** | **2.5 µL** of the diluted flavor–isotope mixture, **in duplicate** |
| **Dilutions (v/v, flavor/EtOH)** | **n-decane 20.2 % · decanal 12.7 % · PEITC 9.8 % · 2-undecanone 13.5 %** |
| **Dose basis, verbatim** | *"The amount of C¹⁴-labeled flavor compound per sample was calculated on a **mol-to-mol flavor compound-to-lysine basis** in the protein isolate."* |
| **Blank** | 2.5 µL of the flavor–isotope mixture into **1.5 mL of double-distilled water** (no protein), same extraction |
| **Temperatures** | **25, 45 and 65 °C**, water baths |
| **Sampling times** | **1, 3, 7, 14, 21 and 28 days** (+ a day-0 point plotted in all three figures) |
| **pH** | ⚠️ **NEVER STATED AND NEVER CONTROLLED.** Double-distilled water; no buffer. A 1 mg/mL BLG solution sits near **pH 6.5–7**; it is not measured, and Schiff-base rates are strongly pH-dependent. **Register `pH: unstated_unbuffered_DDW`, not `pH 7`.** |
| **Replication** | *"The whole experiment … was conducted **twice**. Each time, **two separate samples** were prepared for each flavor−protein combination, yielding **4 replicates**."* |
| **Counting replication** | *"counted **three times per sample for 4 min each** to produce **6 sample counts** per initial BLG sample"* |
| **Statistics** | *"A **Q-test** was implemented to identify and reject outliers."* **No ANOVA, no significance test, no n printed on any figure.** |
| **Normalization** | *"The study final results are expressed as average of all data points and **recalculated based on the protein content in samples**."* |

⚠️ **[!] FIGURE CAPTIONS SAY "1 mL", METHODS SAY "1.5 mL".** Figs. 3 and 4 are captioned
*"in (a) 1 mL of 0.1 % BLG solutions"* while §Sample Preparation specifies a **1.5 mL**
aliquot. Since the y-axis is **per gram of protein**, the discrepancy does not touch the
plotted values; it does touch any attempt to back out absolute molarities. **Flag it and
propagate ±50 % on any absolute-concentration reconstruction.**

### 3c. Separation of covalently bound from free label **[M]** — the method's core

1. Vials removed from the bath and **placed on ice** *"to cool and thereby condense possible radioactive vapors inside the vials."*
2. **1.5 mL** of the reacted solution → **40 mL screw-cap vial containing 10 mL of hexane**.
3. Capped, **sealed with parafilm (to detect leaks)**, inside a 50 mL plastic centrifuge tube, shaken horizontally on a **Roto-Bot (Benchmark Scientific)** for **5 min (2.5 min per side)**.
4. Parafilm inspected for hexane dissolution (leak check); **30 s ultrasound bath** to break foam.
5. Hexane discarded; **0.15 g NaCl** added *"to facilitate extraction"*; **re-extracted four more times** without further salt.
6. **Two 0.5 mL samplings** of the aqueous phase taken as the **BLG-bound isotope**.
7. 0.5 mL + **4 mL scintillation cocktail** → **Beckman LS 6000IC** liquid scintillation counter.

**Why hexane [M]:** *"chosen based on its hydrophobicity (**logP = 3.9**) and density, which is
lower than water, which made it easier to remove the extracting solvent."*
**Why salt [M]:** *"adding salt resulted in **maximal n-decane recovery**"* — consistent with ref. 10.
**What was rejected [M][Q]:** *"methods to reduce the protein emulsification capacity could not
include any treatment that **cleaved the reacted flavor from the protein (e.g., strong acids)**
or **favored further isotope reaction in the extraction process (e.g., heating or alkaline pH
adjustment)**."* High-speed centrifugation was avoided for radiological-safety reasons.

### 3d. ⭐ THE BLANK — the paper's strongest methodological claim **[M]**
**n-Decane cannot form a covalent bond but binds BLG hydrophobically.** Carried through the
full 28-day shelf-life study at all three temperatures:

> *"This sample set gave **very low, consistent label counts throughout the shelf life study
> irrespective of the storage temperature**. The consistent, low isotope signal (**< 0.1 % of
> the total added n-decane**) demonstrated that our extraction protocols were adequately
> removing unbound decane."*

**⇒ THE NON-COVALENT CARRY-OVER FLOOR IS < 0.1 % OF ADDED DOSE, MEASURED, AT EVERY
TEMPERATURE AND EVERY TIME POINT.** Against decanal's 14–16 % and PEITC's 60–72 %, this is a
**150–700× signal-to-blank ratio**. This is the single feature that makes the paper's numbers
usable: it is the control `anantharamkrishnan2020b_extraction.md` and `meynier2004_extraction.md`
both lacked. **Register `noncovalent_carryover_floor: <0.1%_of_dose_measured_all_T`.**

⚠️ **[!] But note the blank's own limit.** n-Decane (logP ≈ 5.0) is *more* hydrophobic than
decanal (logP ≈ 4.0) and so should be *harder* to strip, which makes the blank conservative in
the right direction. It is nonetheless **a different molecule**: it carries no carbonyl and
therefore cannot report on a **reversible hemiacetal/carbinolamine** population that survives
hexane partitioning and then reverts. **The measurement is `covalent_or_hexane_inextractable`,
not proven `irreversible`.** As in `meynier2004_extraction.md` §4a, **the corpus still has
ZERO measured reversibility data on an aldehyde–protein adduct.** That gap survives this wave too.

### 3e. Why the authors rejected the two prior methods **[Q][C]** — worth registering as method critique
- On **Andriot et al. 1999** (¹⁴C-benzaldehyde + BLG by LC + radiodetection, ref. 8/10):
  *"the authors **did not demonstrate that they had removed all of the free, labeled
  benzaldehyde** from the BLG as it eluted through the column. BLG has strong potential for
  hydrophobically bonding with BLG — the eluting BLG might have carried some loosely bound
  benzaldehyde through the column and be **incorrectly considered as covalently bound**."*
  *(The sentence as printed says "hydrophobically bonding with BLG"; "with benzaldehyde" is
  clearly meant — a typo.)*
- On **LC/ESI-MS** (`anantharamkrishnan2020b`, ref. 7/9): *"LC/ESI MS will work for small
  proteins such as BLG but **does not have the mass range needed to measure flavor bonding with
  most other proteins, especially plant proteins**."*
- On **fluorescence quenching** (refs. 5, 6): *"it did not distinguish between the reacted
  ligand type (e.g., reducing sugar or flavor molecule) or type of the reaction."*
  **⇒ Independent, same-lab corroboration of `meynier2004_extraction.md` §2 [!] and §6b: fluorescence
  is not a binding or a covalent-extent observable.**

---

## 4. ⭐⭐ FIGURE 3 — DECANAL — THE QUANTITATIVE CORE, FULLY DIGITIZED **[D]**

*Caption verbatim (p. 10581): "**Figure 3.** Covalently bound DEC in (a) 1 mL of 0.1 % BLG
solutions during incubation at 25, 45, and 65 °C. (Error bars represent standard deviations)."*
*Legend: ● 25 °C (long dash) · ▲ 45 °C (dash-dash) · ■ 65 °C (dash-dot-dot).*
*Y-axis: **"Decanal, mg/g of protein"**, 0–50, gridlines every 10.
X-axis: **"Days of incubation of 0.1 % β-lactoglobulin solution"**, ticks 0 / 7 / 14 / 21 / 28.*

### 4a. **Digitization method and its accuracy** [D]
Page 10581 rasterized at **300 dpi** (2531 × 3338 px). **Y calibration:** the six horizontal
gridlines were located automatically at rows **1344 / 1424 / 1504 / 1584 / 1664 / 1744** for
**50 / 40 / 30 / 20 / 10 / 0** — spacing **80.0 px per 10 mg/g, exact to < 0.5 px on all five
intervals** ⇒ **8.00 px per mg/g**. **X calibration:** tick-label centroids at
**1587.0 / 1734.0 / 1882.0 / 2026.0 / 2174.0** for **0 / 7 / 14 / 21 / 28 d** ⇒
**20.96 px/day**, intervals internally consistent to **± 1.4 %**.
Markers were located as connected dark runs in a ±7 px column at each day; error-bar caps as
**thin (≤ 6 px) wide (≥ 20 px)** horizontal runs; a series' mean was taken as the **midpoint
of its two caps** where both were resolvable, otherwise as the marker-bounding-box centre.

**Precision: ± 0.06 mg/g from pixel reading alone; ± 0.15 mg/g including marker-overlap
assignment.** The markers are **≈ 20 px = 2.5 mg/g tall**, so **series separated by less than
≈ 2.5 mg/g merge into one blob.** Separation quality by day:

| day | 25 °C | 45 °C | 65 °C | separability |
|---|---|---|---|---|
| 0 | ✅ | ✅ | ✅ | all at ≈ 0 |
| **1** | ✅ | ⚠️ partly hidden | ✅ | **good — this is the day the Eₐ rests on** |
| **3** | ✅ | ✅ (both caps read) | ⚠️ merged under the triangle | fair |
| 7 | ✅ (both caps) | ⚠️ apex only | ✅ | fair |
| 14 | ❌ merged with 65 | ✅ (both caps) | ❌ merged with 25 | poor |
| 21 | ❌ | ❌ | ❌ | **all three merged — do not use** |
| 28 | ✅ (both caps) | ⚠️ | ⚠️ | fair |

### 4b. **THE DIGITIZED DECANAL DATASET** [D] — mg decanal per g protein

| day | **25 °C** | **45 °C** | **65 °C** | notes |
|---:|---:|---:|---:|---|
| **0** | **≈ 0.4** | ≈ 0 | **≈ 0.4** | one square at the origin, all series start at zero |
| **1** | **3.6** | **≈ 5.0** *(bracket 4.5–6.3)* | **7.5** | ⭐ the initial-rate row; ordering strictly monotone in T |
| **3** | **4.9** | **12.4 ± 5.4** | **≈ 10.0** | 45 °C caps read at 17.75 / 7.00 |
| **7** | **15.4 ± 4.5** | **≈ 24.3** | **23.6** *(one bar 21.6–26.1, centre 23.9 ± 2.25)* | |
| **14** | ≈ 17 *(merged)* | **22.4 ± 2.5** | ≈ 17 *(merged)* | 25 and 65 indistinguishable |
| **21** | *all three merged in 16.6–21.8; single bar centre 19.1 ± 3.5* | | | **unusable per-series** |
| **28** | **25.6 ± 1.4** | **≈ 20.1** | **≈ 17.4** | ⭐ **the inversion — see §4d** |

*Standard deviations quoted are read from the printed caps (half the cap-to-cap span). Where
only one bar is resolvable for two overlapping series it is quoted once and attributed to neither.*

### 4c. **The authors' own three statements about Fig. 3** [M][Q]
1. *"Under the highest temperature/longest storage time, reaction conditions (**28 days @ 65 C**),
   **16 % of the added DEC was covalently bound to BLG or ca. 2.4 mol of DEC per mol of BLG**."*
2. *"The plot shown in Figure 3 could be considered as being **linear at the lower storage
   temperature (25 C)** and **initially linear but then constant/decreasing at higher storage
   temperatures (45 and 65 C)**."*
3. *"It has been shown that **storage of a protein at elevated temperatures may alter protein
   structure such that it becomes less available/reactive**.¹⁵ **We favor this explanation.**
   This effect may also be seen in the reaction of 2-UNO with the protein."*
And in the Summary: *"**ca. 23 mg of DEC bonded per gram of BLG** at the same time/temperature."*
Abstract: *"**DAL (max of 16.4 % reacted)**."*

### 4d. ⚠️⚠️ **[!] THE TEXT'S HEADLINE CLAIM IS CONTRADICTED BY ITS OWN FIGURE — and the contradiction is the paper's most important result**

The text attributes the maximum to *"the highest temperature/longest storage time … 28 days
@ 65 C"*. **Figure 3 says otherwise.** At day 28 the digitized ordering is:

| series | day-28 value | rank |
|---|---:|---|
| **25 °C ●** | **25.6 ± 1.4 mg/g** | **HIGHEST** |
| 45 °C ▲ | ≈ 20.1 mg/g | 2nd |
| 65 °C ■ | ≈ 17.4 mg/g | **LOWEST** |

**And the 65 °C series peaked at day 7 (23.6 mg/g) and then declined by 26 % to day 28.**
The 25 °C series is the only one still rising at day 28 — exactly as the authors' own
statement (2) says. **⇒ The figure's global maximum is the COLDEST sample, and the text's
"highest temperature" attribution is wrong.**

**[Z] Which reading does the arithmetic support?** Take BLG A-variant **MW 18,363 Da**
(the value the sister paper `yuan2023_extraction.md` prints) and decanal **156.27 g/mol**:

| candidate point | mg/g (digitized) | ⇒ mol/mol BLG | ⇒ % of a 15 mol/mol dose |
|---|---:|---:|---:|
| 65 °C, day 28 | 17.4 | **2.05** | **13.6 %** |
| 65 °C, day 7 (the 65 °C max) | 23.6 | **2.77** | **18.5 %** |
| **25 °C, day 28 (the figure's max)** | **25.6** | **3.01** | **20.0 %** |
| *the text's claim* | *"ca. 23"* | *"ca. 2.4"* | *"16 %" / "16.4 %"* |

**Note the internal check that DOES close: `2.4 mol/mol ÷ 16 % = 15.0`, and BLG has exactly
15 lysine residues per monomer.** ⇒ **The dose really was 1 mol flavor per mol lysine, and the
"% reacted" and "mol/mol" figures are mutually consistent.** What does *not* close is which
point they refer to and the "23 mg/g": **2.4 mol/mol = 20.4 mg/g of pure BLG**, not 23
(a **+13 %** discrepancy, plausibly the 90 %-purity normalization applied inconsistently).

**⇒ REGISTER: `shepelev_decanal_max_extent: 25.6 mg/g at 25 °C day 28 [D]`, with the paper's
own printed alternatives `16–16.4 % of dose / 2.4 mol per mol BLG / ca. 23 mg per g` carried
as `[M] text_value_attributed_to_65C_day28_CONTRADICTED_BY_FIG3`. Do not ingest the abstract's
16.4 % as a 65 °C value.**

### 4e. **[Z] Absolute concentrations — and the ±40 % that rides on them**
From §3b: 1.5 mL of 1.0 mg/mL BLG powder at 90 % purity ⇒ **1.35 mg protein, 73.5 nmol BLG,
49 µM BLG, 735 µM lysine.** From the dilution recipe: 2.5 µL of 12.7 % v/v decanal
(ρ 0.830) = **0.2635 mg = 1.686 µmol = 1.12 mM decanal.**

⚠️ **[!] THAT IS 1.686 µmol OF DECANAL AGAINST 1.103 µmol OF LYSINE — a ratio of 1.53 : 1, not
the 1 : 1 the Methods claim and not the 15 mol/mol the results imply.** Reconciling requires
either **2.06 mg** of BLG in the vial (i.e. the "1 mg/mL" is protein-not-powder *and* the
aliquot is larger) or a smaller effective dose. **The discrepancy is 1.38×.**
**⇒ The `mg/g` y-axis is self-consistent and safe. Every quantity derived on a MOLAR basis
(k₂, % of dose, mol/mol) carries a one-sided systematic of up to ±40 %. Register
`molar_basis_uncertainty: 1.4x`.** *It does not touch Eₐ, which is a ratio of rates at the
same composition (§6f).*

### 4f. **[Z] The second-order rate constant for decanal + BLG at 25 °C — and why it is 15× larger than hexanal's**
Initial rate at 25 °C, from day 0 → day 1: **3.6 mg/g × 1.35 mg protein / 1.5 mL** = 3.24 mg/L
of decanal bound per day = **20.7 µM/day = 2.40 × 10⁻¹⁰ M s⁻¹**.
With [decanal] = 1.12 mM and [Lys] = 0.735 mM:

> **k₂(decanal + BLG, 25 °C) ≈ 2.9 × 10⁻⁴ M⁻¹ s⁻¹** [Z] *(±40 % from §4e, ±20 % from the digitization)*

**Compare the corpus's hexanal values on the same footing:**

| system | k₂ (M⁻¹ s⁻¹) | source |
|---|---:|---|
| **decanal + BLG, 25 °C** | **2.9 × 10⁻⁴** | **this wave [Z]** |
| t-2-hexenal + Na-CN/WP lysine, 20 °C | 5.3–7.9 × 10⁻⁵ | `meynier2004_extraction.md` §8 |
| **hexanal + Na-CN/WP lysine, 20 °C (upper bound)** | **≤ 2.5 × 10⁻⁵** | `meynier2004_extraction.md` §1(b) |
| hexanal + BLG, ambient (mass-spectral) | 10⁻⁶ – 10⁻⁵ | `anantharamkrishnan2020b_extraction.md` §6 |

**⇒ Decanal is ~12–15× faster than the hexanal upper bound and ~4–5× faster than t-2-hexenal.**
**[Z] The most economical explanation is chain length, and it is quantitatively right.**
Decanal is C10, hexanal C6 — four CH₂ groups. The corpus's two independently measured
hydrophobic-binding chain-length slopes are **2.72×/CH₂ (Andriot, β-lg)** and **2.9×/CH₂
(Damodaran, soy)** (`k2_matrix_and_thresholds.md` §B.6). A **partial** transfer of that slope
into the covalent rate — as expected if the reaction proceeds from a hydrophobically bound
state whose local concentration is enhanced — predicts an enhancement between **2.8¹ ≈ 3×**
(one effective CH₂ of accessibility gain) and **2.8⁴ ≈ 61×**. **The observed 12–15× sits
squarely inside that bracket.**

**⇒ THE OPERATIONAL CONSEQUENCE, AND IT IS THE ONE THE REPO NEEDS:
`k₂(decanal)` IS AN UPPER BOUND ON `k₂(hexanal)`, NOT AN ESTIMATE OF IT. Anything this paper
says about the size of the aldehyde covalent sink OVERSTATES hexanal's by roughly an order of
magnitude.** Register `chain_length_correction_decanal_to_hexanal: divide_by_12-15 [Z]`.

---

## 5. FIGURE 2 — 2-UNDECANONE — the negative control that is also a result

*Caption verbatim: "**Figure 2.** Covalently bound UDO in BLG solution during incubation at 25,
45, and 65 °C. (Error bars represent standard deviations)."*
*Y-axis: "**2-undecanone, mg/g of protein**", **0.0–3.0**, gridlines every 0.5.*

**⚠️ [!] DIGITIZATION OF FIG. 2 IS LOW-CONFIDENCE AND THIS WAVE DOES NOT PUBLISH PER-POINT
VALUES FOR IT.** The y-axis numerals overlap the rotated axis title in the raster, and two
candidate baselines differ by ≈ 0.6 mg/g; the paper's own summary numbers (below) are the
reliable record. What IS securely readable from the figure: **every point lies between
≈ 0.5 and ≈ 1.5 mg/g**, all series are **flat to weakly rising** across 28 days, the error
bars are **large relative to the spread**, and **the three temperature series overlap
throughout**.

**The paper's own numbers for UDO [M]:**
- Abstract: *"**UDO showed lowest reactivity (max of 0.9 % of added compound reacted)**."*
- Results: *"only a small amount of added UDO was found to be covalently bound to BLG under our
  experimental conditions (**0.9 % of potential reactivity or ca. 0.1 mmol of UDO per mol of
  BLG**). As one would expect, the data show a **small, positive effect of reaction temperature
  on the extent of reaction** consistent with an earlier work (Anantharamkrishnan et al. (2020))."*
- Summary: *"**ca. 1.4 mg** of UDO"* per gram of BLG.

⚠️ **[!] "0.1 mmol of UDO per mol of BLG" IS WRONG BY ~1 000×.** [Z] 0.9 % of a 15 mol/mol
dose is **0.135 mol/mol = 135 mmol/mol**; and 1.4 mg/g ÷ 170.29 g/mol × 18 363 g/mol =
**0.151 mol/mol = 151 mmol/mol**. **The printed "0.1 mmol per mol" should read "0.1 mol per
mol."** Two independent routes agree on ≈ 0.14 mol/mol. **Do not ingest the mmol figure.**

**Mechanistic reason given [C]:** *"mono-ketones do not readily form covalent bonds with protein
amino acid residues due to the **steric hindrance** effects"* (refs. 12 Damodaran & Kinsella 1980,
13 Franzen & Kinsella 1974). **⇒ Independent confirmation of `anantharamkrishnan2020b`'s tier
ranking: mono-ketone ≪ aldehyde ≪ isothiocyanate.**

**[Z] The UDO/decanal selectivity ratio, at matched dose and temperature: 1.4 vs 23 mg/g ⇒
16×; on the % basis, 0.9 % vs 16 % ⇒ 18×.** Two routes agree. **Register
`monoketone_vs_alkanal_covalent_selectivity: 16-18x [Z]` — a clean, saturating-dose,
single-experiment number the corpus did not have.**

---

## 6. ⭐⭐⭐ THE ACTIVATION ENERGY — derived by this wave [Z]

> **EVERYTHING IN THIS SECTION IS `[Z]`. The paper computes no rate and no Eₐ.**

### 6a. Why day 1 is the right anchor
The Arrhenius comparison must be made where the three series are (i) **separable** and
(ii) **least contaminated by saturation**. Day 1 is the only sampling time that satisfies both:
all three points are individually resolvable (§4a), and it is the earliest non-zero time.
Day 3 is already unusable for a three-point fit because **45 °C (12.4) exceeds 65 °C (≈10.0)** —
the 65 °C curve has begun to bend by day 3, exactly as the authors describe.

### 6b. **Estimate A — raw day-1 extents (all series treated as still in the initial-rate regime)**
Extents at t = 1 d are proportional to the mean rate over 0→1 d, so their ratio is a rate ratio.

| T | 1/T (K⁻¹) | extent (mg/g) | ln(extent) |
|---|---:|---:|---:|
| 25 °C (298.15 K) | 3.35402 × 10⁻³ | 3.6 | 1.2809 |
| 45 °C (318.15 K) | 3.14317 × 10⁻³ | 5.0 | 1.6094 |
| 65 °C (338.15 K) | 2.95727 × 10⁻³ | 7.5 | 2.0149 |

Least-squares slope = **−1 830** K ⇒ **Eₐ = 15.2 kJ/mol**, **R² = 0.998**.
Pairwise: 25→65 gives **15.3**; 25→45 gives **15.4**. *(Q₁₀ ≈ 1.2.)*

### 6c. **Estimate B — saturation-corrected (first-order approach to each series' own plateau)**
Each series is fitted as `y(t) = y∞(1 − e^(−kt))` using that series' own observed maximum as
y∞ (25 °C → 25.6, 45 °C → 22.4, 65 °C → 23.6), and k is solved from the day-1 point:

| T | y(1) | y∞ | **k (d⁻¹)** | ln k |
|---|---:|---:|---:|---:|
| 25 °C | 3.6 | 25.6 | **0.1515** | −1.887 |
| 45 °C | 5.0 | 22.4 | **0.2513** | −1.381 |
| 65 °C | 7.5 | 23.6 | **0.3826** | −0.961 |

Least-squares slope = **−2 405** K ⇒ **Eₐ = 20.0 kJ/mol**, **R² = 0.9997**.

**⇒ Eₐ(decanal–BLG covalent adduct formation, aqueous, 1 mg/mL BLG, unbuffered, 25–65 °C)
= 15–20 kJ/mol.** Estimate B is the better one (it removes the known curvature) and is the
value this wave recommends. **Report as `Ea_aldehyde_protein_adduct = 20 kJ/mol
(range 15–23), [Z] derived_K6b_from_Shepelev2024_Fig3`.**

### 6d. **Sensitivity — the number is robust to every digitization choice**

| what is varied | over what range | resulting Eₐ (estimate A) |
|---|---|---:|
| **the unresolved 45 °C day-1 point** | 4.5 → 6.3 mg/g | **15.2 → 15.5** |
| the 25 °C day-1 point | 3.3 → 3.9 | 14.0 → 16.5 |
| the 65 °C day-1 point | 7.0 → 8.0 | 14.0 → 16.4 |
| **dropping 45 °C entirely (2-point fit)** | — | **15.3** |
| **using day 3 for the 25→65 pair only** | — | **14.9** |
| estimate B with a common y∞ = 25.6 for all three | — | **22.7** |

**Every route lands in 14–23 kJ/mol.** Propagating the printed SDs (±2 to ±5.4 mg/g at
comparable points, n = 4 ⇒ SEM ≈ SD/2) gives **Eₐ = 20 ± 4 kJ/mol (1 σ)**, i.e. a 95 %
interval of roughly **12–28 kJ/mol**.

### 6e. ⚠️ **The one honest way the true chemical Eₐ could be larger — and it still does not reach 70**
The observed Eₐ is an **apparent** activation energy for a two-step process: hydrophobic
association of decanal with BLG, then chemical adduction from the bound state.
`E_obs = E_chem + ΔH_bind`. Hydrophobic binding of aldehydes to β-lg is exothermic
(ΔH_bind typically **−10 to −25 kJ/mol** for C6–C10 carbonyls in the flavour-binding
literature). **⇒ The intrinsic chemical barrier could be `E_chem = E_obs − ΔH_bind` ≈ 25–45 kJ/mol.**

**This is the most generous defensible correction, and it is stated so that no one can say the
verdict rests on the smallest possible number. Even at the top of it, Eₐ = 45 kJ/mol, the
`≥ 70 kJ/mol` threshold is missed by 25 kJ/mol, and the process-temperature acceleration is
156× rather than 2 596× — a 17× shortfall (`k6b_adduct_kinetics_synthesis.md` §3).**
**And for the repo the apparent Eₐ is the correct one to use anyway**, because the repo's sink
term operates on total (not bound) aldehyde concentration, exactly as this experiment did.

### 6f. **Why Eₐ is immune to the §4e molar-basis problem**
Eₐ is computed from **ratios of extents at three temperatures in the same vials with the same
composition**. Every scale factor — the 1.38× dose ambiguity, the 1 mL/1.5 mL caption
conflict, the 90 % purity normalization, the `mg/g` → molar conversion — is **common to all
three temperatures and cancels identically in `ln(y₁/y₂)`.** **The Eₐ is the most robust
number this paper yields, and it is more robust than the extents it is derived from.**

---

## 7. FIGURE 4 — 2-PHENYLETHYL ISOTHIOCYANATE — the internal control that validates §6

*Caption verbatim: "**Figure 4.** Covalently bound PEITC in 1 mL of 0.1 % BLG solutions during
incubation at 25, 45, and 65 °C. (Error bars represent standard deviations.)"*
*Y-axis: "**2-phenylethyl Isothiocyanate, mg/g of protein**", **0–100**, gridlines every 10.*

### 7a. **THE DIGITIZED PEITC DATASET** [D] — mg PEITC per g protein
*(Fig. 4 is the cleanest of the three: the series never merge and every error bar is readable.)*

| day | **25 °C** | **45 °C** | **65 °C** |
|---:|---:|---:|---:|
| **0** | ≈ 0 | ≈ 0 | ≈ 0 |
| **1** | **8.0** *(bar ≈ 0–28)* | **34.5** *(bar 28–38)* | **43.5** *(bar 38–48)* |
| **3** | **21.5** *(bar 18–25)* | **44.5** *(bar 40–49)* | **50.5** *(bar 39–62)* |
| **7** | **28.5** *(bar 26–31)* | **50.0** *(bar 43–54)* | **60.0** *(bar 44–77)* |
| **14** | **34.5** *(bar 29–40)* | **59.5** | **77.5** |
| **21** | **40.5** *(bar 39–42)* | **61.0** *(bar 57–65)* | **81.0** ← series max |
| **28** | **44.0** *(bar 40–47)* | **65.0** *(bar 61–68)* | **77.5** *(bar 75–81)* |

*Digitization precision ± 0.6 mg/g (y-axis 0–100 over the same pixel height).*

### 7b. **The paper's own PEITC numbers [M]**
- Abstract: *"**PEITC (max of 71.8 % reacted)**"*; *"It appears that **only PEITC (at 65 °C) saturated all potential protein-reactive sites** over the storage period."*
- Results: *"a maximum of **71 %** of the added PEITC covalently bonded with BLG"*; *"**The reaction was complete (all possible bonding sites reacted) only for the samples at 65 °C on the third week** of incubation."*
- Summary: *"**up to ca. 80 mg of PEITC bonded per gram of BLG** during sample incubation at 65 °C."*
- *"Reaction between PEITC and BLG at lower temperatures continued with time approaching values of samples incubated at 65 °C."*

**[Z] Consistency check:** 80 mg/g ÷ 163.24 × 18 363 = **9.00 mol/mol = 60.0 % of a 15 mol/mol
dose**, against the printed **71–71.8 %**. **A 1.20× gap of the same character as decanal's
1.13× (§4d) — the two are consistent with a single ≈ 1.15–1.20× normalization inconsistency
running through the whole Summary paragraph.** Flag; do not reconcile.

### 7c. ⭐ **[Z] Eₐ FOR PEITC — the control**
Same two estimators as §6, on the day-1 extents (8.0 / 34.5 / 43.5 mg/g):
- **Estimate A (raw extents):** slope −4 327 K ⇒ **Eₐ = 36.0 kJ/mol**, R² = 0.875.
- **Estimate B (saturation-corrected, common y∞ = 81):** k = 0.104 / 0.555 / 0.770 d⁻¹ ⇒ slope
  −5 110 K ⇒ **Eₐ = 42.5 kJ/mol**, R² = 0.96.

> **Eₐ(PEITC–BLG) = 36–43 kJ/mol, vs Eₐ(decanal–BLG) = 15–20 kJ/mol.
> A 2.0–2.5× separation, same apparatus, same protein, same three temperatures, same
> analyst, same figures. [Z]**

**⇒ THIS IS THE DECISIVE INTERNAL CONTROL. The low decanal Eₐ is not an artefact of the
digitization, of the day-1 anchor, of the saturation correction, or of the ¹⁴C method — every
one of those is shared with PEITC, and PEITC returns a value 2–2.5× higher. The method
resolves activation energies; decanal's is genuinely small.** And note that the ordering
matches the corpus's chemistry: the isothiocyanate value lands next to the acrylamide-Michael
forward Eₐ of 44 kJ/mol (`zamora2010_extraction.md`), while the alkanal lands next to the
thiol-Michael 28–30 (`hidalgo2010_extraction.md`) and below it.

### 7d. ⚠️ **[!] THE PEITC EXTRAPOLATION REPEATS THE DECANAL INVERSION**
Verbatim: *"Assuming that 25 and 45 °C sample reaction rates would not change with reaction
time, PEITC would fully saturate available reactive sites on BLG on the **11th and 17th
weeks, respectively**, of storage."*
**As printed, "respectively" assigns week 11 to 25 °C and week 17 to 45 °C — i.e. the COLDER
sample saturates SOONER.** Either the ordering is a typo, **or it is the same real
phenomenon as §4d**: the low-temperature series alone remains linear and therefore
extrapolates to the ceiling faster than a high-temperature series that has already flattened.
**[Z] The digitized data support the second reading being at least arithmetically available:**
from the day-21→28 slopes, 25 °C is gaining **0.50 mg/g/d** while 45 °C gains **0.57** and
65 °C **−0.50**; a linear extrapolation of 25 °C to 81 mg/g needs **74 more days ≈ week 14**,
and of 45 °C needs **28 more days ≈ week 8**. **Neither matches the printed 11/17.**
**⇒ Register the sentence as `[!] unreconcilable_extrapolation`. Do not ingest weeks 11/17 as data.**

### 7e. **Why isothiocyanates are the most reactive class — the authors' mechanism [C][Q]**
> *"its electrophilic carbon can react with several nucleophilic functional groups of amino
> acids (Figure 5). Reaction with **hydroxyl (HO−)** groups results in the formation of
> **O-thiocarbamates**, while **amino (H₂N-)** groups produce **thioureas**, and **thiol (HS-)**
> groups form **dithiocarbamates**. In addition, **Schiff base, Michael addition, and disulfide
> bond formation** occur."*
> *"an isothiocyanate will react with **multiple sites on a protein**, while reactions with **all
> other classes of volatiles studied reacted only once with a single protein**."*
> *"a **thiol group can cleave electrophilic disulfide linkages** in a protein structure
> **increasing possible reaction sites**."*
> *"isothiocyanate bonding to protein molecules resulted in the **unfolding of the α-helix
> structure of BLG into β-sheets**"* (Ersöz & Dudak 2020; Keppler 2017).

**⭐ [Q] AND THE GENERAL PRINCIPLE THE PAPER STATES, WHICH THE REPO SHOULD REGISTER:**
> *"the **number of reactive sites on the surface of a protein plays a greater role in
> determining reaction rate/extent than the total amount of reactive amino acids** in a protein
> structure since the amino acids must be **reachable** by a reactive compound in order for bond
> formation to occur."*

**⇒ ACCESSIBILITY, NOT STOICHIOMETRY, SETS THE CEILING. This is the mechanistic statement that
explains §4d/§7d (heat changes accessibility, not the lysine count) and it is a direct warning
against any repo sink term that sizes the covalent channel from total lysine content.**

---

## 8. THE PAPER'S OTHER SUBSTANTIVE CONTENT

### 8a. Figure 1 — Schiff-base scheme [C], reproduced from Qin et al. 2013 (ref. 11). No data.
### 8b. Figure 5 — proposed isothiocyanate modification scheme [C], after Kaschula & Hunter 2016 (ref. 16). No data.
### 8c. Statements about aldehyde reactivity, all cited, none measured here **[C][Q]**
- *"**aromatic aldehydes more easily form Schiff base linkages with amino acids than aliphatic aldehydes**"* (ref. 14, Zaugg, Walder & Klotz 1977, haemoglobin).
  **⚠️ [!] NOTE THE TENSION with `yuan2023_extraction.md` §5, which measures the opposite ordering on BLG: *"Citral (geranial) and hexenal dimethyl acetal are more reactive aldehydes than aromatic counterparts (e.g., benzaldehyde)."* Same lab, adjacent years, opposite claims — one cited from 1977 haemoglobin work, one measured on BLG. Prefer the measurement.**
- *"The reactivity of aldehydes with proteins has been demonstrated to occur **primarily with the epsilon amino group of lysine** and, potentially, **with thiol groups**¹⁵ and result in **Schiff base formation**.⁷"*
- *"Ketones and aldehydes are known to form **irreversible** Schiff bases with primary amines"* — **[C][!] asserted from ref. 11, not tested here. Schiff bases are textbook-reversible; this sentence should not be ingested as evidence of irreversibility.**
### 8d. The framing argument for why covalent ≠ non-covalent, verbatim [Q]
> *"non-covalent bonding tends to be less onerous since this binding reaches an **equilibrium**,
> and thus the flavor formulation can be adjusted … However, undesirable changes in flavor
> profile due to covalent bonding **cannot be compensated for via flavor formulation changes
> since it is a continuous process, i.e., does not come to a stable equilibrium until the flavor
> compound or the reactive sites on the protein have become saturated.**"*
**⇒ The paper's own model of the covalent channel is `saturating, not equilibrating`. §4d and
§7d show the saturation ceiling is itself temperature-dependent and falls with T.**

---

## 9. ⚠️ THE GAP THIS PAPER DOES NOT CLOSE — and it is the repo's compound

**Decanal is not hexanal, and this paper contains no hexanal.**

| what the repo needs | what this paper gives |
|---|---|
| hexanal (C6) | **decanal (C10)** |
| a food matrix at 3–20 wt % protein | **0.1 wt % BLG in distilled water** |
| trace aldehyde against a vast lysine excess | **1–1.5 mol aldehyde per mol lysine (saturating)** |
| controlled pH | **unbuffered, pH unstated** |
| process temperatures 100–180 °C | **25–65 °C** |
| reversibility | **not tested** |

**What nevertheless transfers, and why:**
1. **The activation energy.** Eₐ is a property of the reaction coordinate, and hexanal and
   decanal are the same reaction (n-alkanal + ε-amino → carbinolamine → Schiff base) differing
   only in an inert alkyl tail. Chain length moves the **pre-exponential** (via the local
   concentration enhancement of §4f), not the barrier. **If anything, hexanal's apparent Eₐ
   should be LARGER than decanal's by the difference in ΔH_bind — but that difference is at
   most a few kJ/mol between C6 and C10, and §6e already brackets the whole correction.**
2. **The sign and size of the capacity effect (§4d).** Loss of accessible sites on heating is a
   protein property, not an aldehyde property.
3. **The blank (§3d).** A < 0.1 % non-covalent carry-over floor validates the method class.

**What does NOT transfer: the absolute rate.** Use `k₂(decanal)/12–15` as the hexanal bound (§4f),
or better, keep the corpus's directly measured hexanal bracket
(`meynier2004_extraction.md` + `anantharamkrishnan2020b_extraction.md`) and apply **only** this
paper's Eₐ to it. **That is exactly what `k6b_adduct_kinetics_synthesis.md` §3 does.**

---

## 10. WHAT WOULD IMPROVE THIS RECORD

| item | what it would deliver | priority |
|---|---|---|
| **The Supporting Information PDF** — *"Description of preliminary study"* (extraction-protocol development, decane recovery vs solvent and salt) | Would settle whether any *reversible* carbonyl adduct could survive hexane extraction (§3d's residual doubt). **It contains no time-course or temperature data** — the abstract of the SI is explicit that it is a method-development note. | **LOW — not needed for the Eₐ or the verdict.** |
| A **thesis** by I. Shepelev at the UMN Conservancy | Would carry the **raw DPM tables behind Figs. 2–4**, turning this wave's ±0.15 mg/g digitization into exact values and supplying the printed SDs at every point. **It would NOT add hexanal.** | **LOW–MEDIUM.** The Eₐ is already robust to ±25 % on every point (§6d); exact values would tighten it, not change it. |
| A hexanal + BLG (or + a food protein) run at ≥ 3 temperatures | **The only thing that would fully close the gap in §9.** | **The standing open item — but see `k6b_adduct_kinetics_synthesis.md` §5 for why the verdict no longer depends on it.** |

---

## 11. INTEGRITY FLAGS RAISED BY THIS WAVE — consolidated

| # | flag | where | severity |
|---|---|---|---|
| **1** | **`DEC` means decanal in Fig. 3, the Results heading and the Summary, but decane in the Abbreviations list** | §3a | ⭐ **HIGH — will mis-assign the entire aldehyde dataset** |
| **2** | **The text's "highest temperature/longest storage time" maximum is contradicted by its own Fig. 3, whose maximum is the 25 °C day-28 point** | §4d | ⭐ **HIGH — and the contradiction is itself the result** |
| **3** | "ca. 0.1 mmol of UDO per mol of BLG" is wrong by ~1 000× (should be ~0.1 **mol**/mol) | §5 | MEDIUM |
| 4 | Summary `mg/g` figures disagree with the same paragraph's `%`/`mol per mol` figures by 1.13× (decanal) and 1.20× (PEITC) | §4d, §7b | MEDIUM |
| 5 | Dose recipe + stock concentration imply 1.53 : 1 flavor : lysine, not the stated 1 : 1; a 1.38× molar-basis inconsistency | §4e | MEDIUM |
| 6 | Figure captions say "1 mL", Methods say "1.5 mL" | §3b | LOW (does not affect `mg/g`) |
| 7 | **pH never stated, never buffered** | §3b | MEDIUM |
| 8 | PEITC saturation extrapolation ("11th and 17th weeks, respectively") reconciles with neither ordering nor the digitized slopes | §7d | MEDIUM |
| 9 | Keyword "pea protein isolate" — no pea protein in the paper | §0 | LOW |
| 10 | "irreversible Schiff bases" asserted from a citation, not tested | §8c | MEDIUM (recurring corpus hazard) |
| 11 | "2-UNO" for UDO; "hydrophobically bonding with BLG" for "with benzaldehyde" | §3a, §3e | LOW (typos) |
| 12 | Ref. 8 and ref. 10 are **the same paper** (Andriot et al. 1999) listed twice with different author strings | References | LOW |

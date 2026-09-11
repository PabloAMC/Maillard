# Schieberle & Hofmann — EXTRACTION (fed 1-pyrroline + hydroxyacetone or methylglyoxal, 2 mmol/L, pH 3-9, 100 C / 30 min; plus the first time course of the 2-acetylpyrrolidine -> AP air oxidation at 25 C)

### The book-chapter companion of Hofmann & Schieberle 1998b: it repeats that paper's ATHP pH ladder and its aqueous/dry AP pair number for number, adds one thing the journal paper does not have (a three-point time course of the final oxidation), and states in words that hydroxyacetone gives ONLY the tetrahydropyridine while methylglyoxal gives ONLY the pyrroline.

**Source on disk:** `data/articles/schieberle2005.pdf` (7 pp., owner's download, 2026-09-08). The text layer is
an OCR layer with heavy glyph damage in prose ("Hofinann", "2-acetyl-l-pyrroline", "Skecker", "*9" for
">99") and it silently DROPPED the pH 9.0 row of Table 2. **Every table below was therefore re-read cell
by cell from 200-dpi rasters of the printed pages 211, 212, 213 and 214** (`scratchpad/img/sch05b-3.png`,
`sch05-4.png`, `sch05-5.png`, `sch05-6.png`); the numbers below are the raster's, not the OCR's. Figures
2-5 are structures and a chromatogram, FIGURE-ONLY, and nothing was read off them.

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanistic Studies on the Formation of the Cracker-like Aroma Compounds 2-Acetyltetrahydropyridine and 2-Acetyl-1-pyrroline by Maillard-type Reactions" |
| Authors | Peter Schieberle (Institute of Food Chemistry, TU München) and Thomas Hofmann (German Research Centre for Food Chemistry, Lichtenbergstrasse 4, 85748 Garching) |
| Venue | Book chapter, running head *The Maillard Reaction in Foods and Medicine* / *Flavour Chemistry*, printed pages **209-215**. The book's title page, editors, publisher and year are NOT in this PDF (flag 1); the file name says 2005 |
| Dating evidence inside the text | reference 13 is a 1996 RSC volume, and references 10 and 11 (the two Hofmann & Schieberle JAFC papers of 1998) are cited as "**submitted**". The chapter was therefore written in 1996-97 and predates the journal paper it duplicates |
| Naming | AP = 2-acetyl-1-pyrroline; ATHP = 2-acetyltetrahydropyridine; APD (also "ADP" in the same table) = 2-acetylpyrrolidine; 2-oxopropanal = methylglyoxal; hydroxy-2-propanone = hydroxyacetone = acetol; the ring-opened ATHP intermediate is called **4,5-dioxoheptylamine** here (`hofmann1998b_extraction.md` reports it as 5,6-dioxoheptylamine — flag 2) |
| Companions on disk | `hofmann1998b_extraction.md` (Hofmann & Schieberle, JAFC 1998, 46, 2270-2277) is the same chemistry with more tables; `blank2003_extraction.md` quotes its yields |

## 1. Why it matters

`results/validation/kinetic_core_b24_prereg.md` section 6 refuses 2-acetyl-1-pyrroline from proline. The
fed acylation of 1-pyrroline by methylglyoxal fits (+0.18 and +0.30 dex); the chain from proline does not
(-1.29, +0.13, +1.16 dex over the methylglyoxal ladder), because the arm carries neither the competing
2-acetyltetrahydropyridine branch through hydroxyacetone nor any loss of 1-pyrroline. The section names
what a next wave needs: hydroxyacetone as a species, the tetrahydropyridine as the competing product, a
pyrroline loss. This chapter is the second reading of the only paper that feeds those precursors. It
matters in three ways and only three:

1. It **confirms** `hofmann1998b_extraction.md` Table 4 to the digit (Table 2 here), so the ATHP branch
   constant and its pH ladder are two independent printings of the same experiment, not one.
2. It states the **branch exclusivity** in words — hydroxyacetone gives only ATHP, methylglyoxal gives
   only AP — which is the structural claim a two-branch network would encode.
3. It prints, uniquely, a **rate for the last step of the AP chain** (2-acetylpyrrolidine oxidising to AP
   in air), which tells the modeller that this step is not the slow one and should not be written as a
   free constant.

It does NOT contain a 1-pyrroline loss (section 4).

## 2. Methods as they matter to a model

- **ATHP pot (Table 2):** 1-pyrroline **10 µmol** + hydroxy-2-propanone **10 µmol** in **5 mL of 0.5 M
  phosphate** -> **2 + 2 mmol/L**; **30 min at 100 C**; pH 3.0 / 5.0 / 7.0 / 9.0. Vessel not stated.
  Identical charge, buffer, time and temperature to `hofmann1998b_extraction.md` Table 4.
- **HOP text statement:** the synthesised intermediate 2-(1-hydroxy-2-oxo-1-propyl)-pyrrolidine heated
  **30 min at 100 C in phosphate buffer** gave "very high yields of ATHP (23.4 % at pH 7.0 and 35.0 % at
  pH 9)". The concentration is not given here (Hofmann's Table 5 gives 1 µmol in 5 mL = 0.2 mmol/L).
- **AP pot (Table 4):** aqueous arm — 1-pyrroline **10 µmol (690 µg)** + 2-oxopropanal **10 µmol** in
  **5 mL of 0.5 M phosphate pH 7.0** (-> **2 + 2 mmol/L**), **30 min at 100 C**. Non-aqueous arm — the
  same 10 + 10 µmol mixed with **silica gel (3 g containing 20 mg KH2PO4 in 300 µL water)**, **5 min at
  180 C**. (Hofmann's Table 7 expt 4 says 2.7 g silica and 300 µL of 0.1 M phosphate; 20 mg KH2PO4 in
  300 µL is ~0.49 mol/L, five times that — flag 3.)
- **Oxidation pot (Table 3):** 2-acetylpyrrolidine trifluoroacetate **45 µg, 0.4 µmol** in **1 mL of tap
  water** (-> **0.4 mmol/L**) containing **dicyclohexylamine 0.4 µmol**, standing at **25 C** in a brown
  glass vial, open to air ("in the presence of oxygen"). Sampled at 5, 30 and 120 min (the column head
  says only "Reaction time"; the text's "within 2 h" fixes the unit as minutes).
- **Labelling pot (Table 1):** **L-proline 1 mmol + D-glucose or [U-13C]-D-glucose 2 mmol** with silica
  gel (3 g) containing **300 µL of 0.1 M phosphate pH 7.0**, **10 min at 160 C**. Note the ratio: proline
  1 : glucose 2. `hofmann1998b_extraction.md` records the same experiment as proline 2 mmol + glucose
  1 mmol, the inverse (flag 4).
- **Quantification:** volatiles recovered by **sublimation in vacuo**; **stable isotope dilution assays**
  with the deuterium-labelled analogues as internal standards. Isotopomer ratios by MS. No response
  factors beyond the labelled internal standards; no replicate count, no error bar anywhere in the
  chapter; no LOD.
- **Yield basis:** "mol %" is molar on the 10 µmol limiting reagent. 10 µmol ATHP = 1251.7 µg, 10 µmol
  AP = 1111.4 µg, 0.4 µmol AP = 44.46 µg. Recomputed below; the printed percentages check out.
- **1-Pyrroline provenance:** "prepared as described previously" (ref. 8, Schieberle 1989). Purity,
  trimer content and recovery are not stated here either.

## 3. Tables re-typed

### Table 1. "Main isotopomers (represented by their molecular ions) determined in 2-acetylpyrroline and 2-acetyltetrahydropyridine generated from proline in the presence of either unlabelled or [U-13C]-labelled glucose (Glc)"

| odorant | unlabelled Glc, m/z (%) | labelled Glc, m/z (%) |
|---|---|---|
| 2-Acetyl-1-pyrroline | 111 (93.9) | 113 (76.8) |
| 2-Acetyltetrahydropyridine | 125 (88.3) | 128 (89.2) |

Only the main isotopomer is printed. `hofmann1998b_extraction.md` Table 6 prints the full distribution
for AP and shows a **second, 19.1 % M+3 isotopomer (m/z 114)** that this chapter omits; the text here
says only "mainly two carbon atoms in the AP stem from glucose". For ATHP the chapter's claim is
stronger and is the useful one: "only the isotopomer with three labelled carbons (m/z 128 vs. m/z 125)
was formed, suggesting **only one reaction pathway in ATHP formation**".

### Table 2. "Influence of pH on the formation of 2-acetyltetrahydropyridine (ATHP) from 1-pyrroline and hydroxy-2-propanone" (10 + 10 µmol in 5 mL of 0.5 M phosphate = 2 + 2 mmol/L; 30 min, 100 C)

| pH | ATHP µg | printed mol % | recomputed mol % (µg / 1251.7) (mine) |
|---:|---:|---:|---:|
| 3.0 | <0.1 | – | <0.008 |
| 5.0 | 0.9 | 0.1 | 0.072 |
| 7.0 | 10.8 | 0.9 | 0.863 |
| 9.0 | 38.4 | 3.1 | 3.07 |

**Identical, digit for digit, to `hofmann1998b_extraction.md` Table 4.** The OCR text layer of this PDF
lost the pH 9.0 row; the raster has it.

### Table 3. "Time course of the formation of 2-acetyl-1-pyrroline (AP) from 2-acetylpyrrolidine (APD)" (0.4 µmol in 1 mL of tap water + 0.4 µmol dicyclohexylamine = 0.4 mmol/L, 25 C, brown glass vial, air)

| reaction time (min) | AP (µg) | printed conversion of ADP (%) | recomputed mol % (µg / 111.14 / 0.4 µmol) (mine) | first-order k from the printed conversion (min^-1) (mine) |
|---:|---:|---:|---:|---:|
| 5 | 12.0 | 26 | 27.0 | 0.060 |
| 30 | 32.5 | 72 | 73.1 | 0.042 |
| 120 | 44.4 | >99 | 99.9 | 0.038 (taking >99 as 99) |

The printed conversions are the molar yields on 0.4 µmol, so the assay closes: the pyrrolidine goes to
AP essentially quantitatively. The three k values fall by a third over the run, so the step is
approximately, not exactly, first order — consistent with a slow oxygen supply into an unstirred vial.
Half-life at 25 C ~12-18 min. **This table has no counterpart in `hofmann1998b_extraction.md`.**

### Table 4. "Influence of the presence of water on the amounts of 2-acetyl-1-pyrroline (AP) generated from 1-pyrroline and 2-oxopropanal"

| reaction system | AP (µg) | recomputed mol % (µg / 1111.4) (mine) |
|---|---:|---:|
| aqueous (10 + 10 µmol, 5 mL of 0.5 M phosphate pH 7.0, 30 min, 100 C) | 58.5 | 5.26 |
| non-aqueous (10 + 10 µmol on 3 g silica gel with 20 mg KH2PO4 in 300 µL water, 5 min, 180 C) | 1.1 | 0.099 |

Ratio 53x; the text says "50-fold". These are `hofmann1998b_extraction.md` Table 7 experiments 2 and 4.
**Experiments 1 (five-fold methylglyoxal, 28.7 mol %) and 3 (five-fold 1-pyrroline, 0.33 mol % of the
methylglyoxal) are NOT in this chapter** — and experiment 3 is exactly the row that would size a
1-pyrroline loss.

### Statements printed in prose (no table)

- Summary: "Reacting 1-pyrroline with hydroxy-2-propanone yielded high amounts of **only the ATHP**,
  whereas the reaction with 2-oxopropanal gave **only AP**."
- Summary: "Synthesized 2-(1-hydroxy-2-oxo-1-propyl)-pyrrolidine was shown to be the key precursor of
  ATHP"; heating it 30 min at 100 C in phosphate gave **23.4 % at pH 7.0 and 35.0 % at pH 9**.
- "the thermal degradation of the amino acid **ornithine**, when reacted with 2-oxopropanal, is **more
  effective in generating AP than proline**" (quoting ref. 7, Schieberle 1995 — not measured here).
- "However, **ornithine did not yield ATHP**" (quoting ref. 9).
- Ring enlargement: 2-methyl-1-pyrroline + hydroxy-2-propanone gave 2-acetyl-3-methyltetrahydropyridine
  as "the main reaction product" — **no yield is printed here** (Hofmann's text gives 3 mol %, "data not
  shown").
- Mechanism (Figure 5): hydrated 2-oxopropanal attacks C-2 of 1-pyrroline -> 2-(1,2-dioxo-1-propyl)-
  pyrrolidine -> **O2** oxidation -> hydration -> rearrangement to 2-acetyl-2-pyrrolidinic acid ->
  **-CO2** -> 2-acetylpyrrolidine -> **O2** -> AP. Two oxygen steps on one chain.
- Six "popcorn-like" odour regions were found by AEDA in the proline-glucose extract, AP the highest FD
  factor. No FD factors are printed.

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`, 75 ids): `2_acetyl_1_pyrroline` is present (added by the B24 wave).
**2-acetyltetrahydropyridine, 1-pyrroline, hydroxyacetone, methylglyoxal, 2-acetylpyrrolidine, HOP,
2-methyl-1-pyrroline and proline all have no row.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| 1-pyrroline + hydroxyacetone -> ATHP | <0.008 / 0.072 / 0.863 / 3.07 | mol % of the 10 µmol charge | pH 3 / 5 / 7 / 9; 2 + 2 mmol/L; 0.5 M phosphate; 100 C; 30 min | Table 2 (recomputed; printed – / 0.1 / 0.9 / 3.1) | fed_intermediate_yield, pH ladder |
| ATHP pH shape | k(9)/k(7) = 3.6; k(7)/k(5) = 12; k(5)/k(3) > 9 | – | same pot | Table 2 (mine) | within_study_ratio |
| apparent second-order constant, 1-pyrroline + hydroxyacetone -> ATHP (mine) | ~1.4e-4 (pH 7); ~5.1e-4 (pH 9); ~1.2e-5 (pH 5) | L/(mmol·min) | Y / (30 min x 2 mmol/L); assumes no loss of either reactant and a bilinear law | derived from Table 2 | derived_assumption, lower bound |
| HOP -> ATHP | 23.4 (pH 7.0); 35.0 (pH 9) | mol % | 30 min, 100 C, phosphate; concentration not stated here | text, p. 212 | fed_intermediate_yield |
| 1-pyrroline + methylglyoxal -> AP, aqueous | 58.5 µg = 5.26 mol % | mol % of the 10 µmol charge | 2 + 2 mmol/L, pH 7.0, 0.5 M phosphate, 100 C, 30 min | Table 4 | fed_intermediate_yield |
| 1-pyrroline + methylglyoxal -> AP, dry | 1.1 µg = 0.099 mol % | mol % | 10 + 10 µmol on 3 g silica, 180 C, 5 min | Table 4 | fed_intermediate_yield |
| AP: aqueous vs dry | 53x (authors: "50-fold") | – | Table 4 rows | Table 4 (mine) | within_study_ratio |
| **2-acetylpyrrolidine -> AP (air oxidation)** | **0.060 / 0.042 / 0.038** | **min^-1 (first order)** | **0.4 mmol/L in tap water + 0.4 mmol/L dicyclohexylamine, 25 C, air, brown glass; from the 5 / 30 / 120 min points** | **Table 3 (k values mine; conversions printed)** | **measured_rate** |
| 2-acetylpyrrolidine -> AP, conversion | 26 / 72 / >99 | % | 5 / 30 / 120 min, 25 C | Table 3 | measured, printed |
| ATHP label pattern from [U-13C]Glc | 89.2 % at M+3 and no other isotopomer printed: three sugar carbons, one route | – | proline 1 mmol + glucose 2 mmol, silica, 160 C, 10 min | Table 1 | measured isotopomer distribution |
| AP label pattern from [U-13C]Glc | 76.8 % at M+2 | – | same | Table 1 | measured isotopomer distribution (partial; the M+3 isotopomer is omitted here, see `hofmann1998b_extraction.md` Table 6) |
| branch exclusivity | hydroxyacetone -> only ATHP; methylglyoxal -> only AP | – | the two fed pots above | Summary, p. 209 | qualitative statement, no number |
| ornithine vs proline for AP | ornithine "more effective"; ornithine gives no ATHP | – | not measured in this chapter | text, p. 209 (refs. 7, 9) | level_only, secondary quotation |

### Can the tetrahydropyridine branch and a 1-pyrroline loss now be written?

**The tetrahydropyridine branch: yes, and the numbers are these.** A step
`1-pyrroline + hydroxyacetone -> ATHP` is carried by a four-point pH ladder at 2 + 2 mmol/L, 100 C,
30 min, printed twice and identically — Table 2 here and `hofmann1998b_extraction.md` Table 4
(<0.008 / 0.072 / 0.863 / 3.07 mol % at pH 3 / 5 / 7 / 9). As a bilinear lower bound that is
~1.4e-4 L/(mmol·min) at pH 7 and ~5.1e-4 at pH 9 (mine, derived_assumption). The step's second half is
independently pinned by `hofmann1998b_extraction.md` Table 5 (HOP -> ATHP, 1.7 / 11.4 / 23.4 / 35.0 mol %
over the same pH ladder), restated here for pH 7 and 9 in the text; comparing the two ladders says the
pH dependence sits mostly BEFORE the HOP intermediate (the addition), not after it, so a single-step
condensation carrying the pH term is the right shape. Whether ATHP or AP wins is `hofmann1998b`
Table 9's within-study ratio (AP : ATHP = 0.16 / 0.51 / 12.8 as methylglyoxal : proline goes 0.01 /
0.1 / 1), and this chapter supplies the mechanism behind it in words: the two dicarbonyls are not
interchangeable, each gives one product only.

**Under which assumption.** Three, and they must be declared, because nothing on disk measures them:
(a) the **barrier of the ATHP step is not measured** — every aqueous pot in both papers is 100 C and
30 min, and the only second temperature is the dry 180 C run, so an activation energy has to be
declared (glycine's, or the AP step's); (b) the **hydroxyacetone supply is not measured** — the
mechanism says the Strecker degradation of an amino acid on methylglyoxal makes hydroxyacetone and
that proline makes 1-pyrroline the same way, but neither paper prints a hydroxyacetone yield from any
pot, so the branch is only as good as the trunk's guess at it, which is the same conditionality
`PYRAZINE_SUPPLY_CAVEAT` already carries for the dicarbonyls; (c) **0.5 M phosphate** throughout, and
phosphate catalyses both steps.

**A 1-pyrroline loss: no — not from a number in either paper.** Nothing here or in
`hofmann1998b_extraction.md` measures 1-pyrroline: no time course, no recovery, no purity statement, no
mass balance. This chapter is strictly weaker than the journal paper on the point, because it prints
only the aqueous/dry pair and **omits Hofmann's Table 7 experiment 3** — the five-fold-excess pyrroline
run (0.33 mol % of the methylglyoxal, sixteen times less per methylglyoxal than the 1:1 run) that is the
sole quantitative hint that 1-pyrroline consumes itself. That single row remains the only handle, and it
is an inference from a yield, not a measurement of a loss: a sink fitted on it would be one free
parameter on one row, with the trimerisation of 1-pyrroline as the named but unmeasured mechanism.

**What this chapter does settle about the AP chain.** Table 3 gives the last step, 2-acetylpyrrolidine
-> AP, a measured first-order constant of ~0.04-0.06 min^-1 **at 25 C**, i.e. a half-life of a quarter
of an hour at room temperature. A step that fast at 25 C is not rate-limiting at 100 C and should be
written as instantaneous (or folded into the acylation constant) rather than fitted; the B24 acylation
constant 2.74e-03 L/(mmol·min) at 100 C already absorbs it. The chain's slow step is the acylation, and
its measured pair of rows is what B24 found did fit.

## 5. Flags

1. **Bibliographic identity is incomplete.** The PDF carries no title page: no editors, no publisher, no
   year, no ISBN, no DOI. All that is printed is the running head (*The Maillard Reaction in Foods and
   Medicine* / *Flavour Chemistry*) and pages 209-215. The file name's "2005" is not evidenced by the
   document; the internal evidence (a 1996 reference, and the two 1998 JAFC papers cited as
   "submitted") dates the writing to 1996-97. Cite it as a chapter with pages and note the uncertainty.
2. **This chapter is not an independent measurement.** Tables 1, 2 and 4 and the HOP percentages are the
   same experiments, with the same charges and the same digits, as `hofmann1998b_extraction.md` Tables
   6/3, 4, 7 and 5. Fitting both as separate rows would double-count. Only **Table 3 is new**.
3. **The dry-run buffer disagrees between the two papers.** Here: 3 g silica with 20 mg KH2PO4 in 300 µL
   water (~0.49 mol/L, and monobasic, so acidic). There: 2.7 g silica with 300 µL of 0.1 M phosphate at
   pH 7. Same yield (1.1 µg) is reported for both, so one description is wrong; the dry arm's pH is not
   trustworthy.
4. **The labelling pot's proline : glucose ratio is inverted between the two papers** (1 : 2 here,
   2 : 1 there) for what is presented as the same experiment giving the same isotopomer percentages.
   The atom bookkeeping is unaffected; the pot composition is not usable.
5. **Table 3's mass and molarity do not both hold.** "2-Acetylpyrrolidine trifluoroacetate (45 µg,
   0.4 µmol)": 0.4 µmol of the trifluoroacetate salt (MW ~227) is 91 µg, while 45 µg is 0.4 µmol of the
   **free base** (MW 113.2). The conversions close on 0.4 µmol (44.4 µg AP = 99.9 %), so read the charge
   as 0.4 µmol and ignore the 45 µg.
6. **Table 3 is not isothermal with anything else in the corpus.** It is 25 C, in tap water (undefined
   ionic strength, undefined pH — dicyclohexylamine is a base), open to air in a single vial, n = 1, no
   error bar. It fixes an order of magnitude, not a constant, and it has no barrier.
7. **Oxygen is a reagent and is never measured.** Figure 5 puts two O2 steps on the AP chain and Table 3
   is explicitly an air oxidation. Nothing in either paper quantifies dissolved oxygen; in a sealed
   extruder barrel this chain could be oxygen-limited and both papers' yields would not transfer.
8. **The AP isotopomer table is truncated here** (m/z 114, 19.1 %, is absent), so a reader of this
   chapter alone would conclude AP has a single route. Use `hofmann1998b_extraction.md` Table 6.
9. **No 1-pyrroline measurement of any kind**, and no replicate counts, error bars, LODs or recoveries
   anywhere in the chapter.
10. **Registry gaps** against `data/keys/compounds.yml`: 2-acetyltetrahydropyridine (the competing
    product a next wave must carry), 1-pyrroline, hydroxyacetone, methylglyoxal (a kinetic-core species
    key only), 2-acetylpyrrolidine, 2-(1-hydroxy-2-oxo-1-propyl)-pyrrolidine, 2-methyl-1-pyrroline,
    2-acetyl-3-methyltetrahydropyridine and proline all have no id. `2_acetyl_1_pyrroline` is the only
    one of this chapter's molecules that is registered.
11. **To request from the authors:** a 1-pyrroline time course or recovery in the fed pot; a
    hydroxyacetone yield from a proline + methylglyoxal pot; either branch at a second temperature.

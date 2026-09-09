# Selke, Frankel & Neff 1978 — EXTRACTION (pure methyl OLEATE hydroperoxides, 27/23/23/27 % 8-/9-/10-/11-OOH, thermolysed neat in a GC injector port at 210 C; 19 peaks by GC-MS; a full four-isomer cleavage-product accounting)

### THE NONANAL NUMBER THE ENGINE REFUSES EXISTS AND IS PRINTED HERE: **nonanal = 15 % of the total volatile peak area** (32.6 % of the non-ester volatiles) from thermally decomposed methyl oleate hydroperoxides — the same laboratory, the same injector-port thermolysis technique and the same first author group as Frankel & Gardner 1989, run 32 C hotter and on the oleate substrate Frankel 1989 never fed.

**Source on disk:** `data/articles/selke1978.pdf` (3 pp., *Lipids* 13 (7), 511-513, 1978; a
"COMMUNICATIONS" short paper). The `pdftotext -layout` text layer is an OCR layer whose running
prose is glyph-spaced ("h y d r o p e r o x i d e s") and whose **Table I is materially corrupt**
(it renders peak 7 as "2-Octanol" for *1-Octanol*, peak 10 as "20", peak 11 as "12", methyl
nonanoate as 2.5 for *1.5*, methyl 9-oxononanoate as 25 for *15*, methyl 10-oxodecanoate as 22 for
*12*, and the triolein t-2-undecenal column as 22/22.1 for *11/12.1*). **Every number below was
re-read off 300-dpi rasters of printed pages 512 and 513**
(`scratchpad/img/selke-t1.png`, `scratchpad/img/selke-t2.png`) and the corrected table closes on all
four of its own printed subtotals (45.9, 47.7, 90.8, 100), which the OCR version does not. Figure 1
(the gas chromatogram) is **FIGURE-ONLY**; no number is read off it. Figure 2 is referenced by
Table I footnote *a* ("See Figure 2 for numbered peaks") but the paper prints only one figure — see
flag 8.

## 0. Identity

| field | value |
|---|---|
| Title | "Thermal Decomposition of Methyl Oleate Hydroperoxides and Identification of Volatile Components by Gas Chromatography-Mass Spectrometry" |
| Authors | **E. Selke, E. N. Frankel and W. E. Neff** — Northern Regional Research Center, Federal Research, Science and Education Administration, U.S. Department of Agriculture, Peoria, Illinois 61604 |
| Venue | *Lipids* **13** (7), 511-513 (1978). Received January 30, 1978. Printed as a COMMUNICATION, not a full paper |
| **DOI** | **NO DOI IS PRINTED IN THE PDF.** There is no DOI, no CrossRef stamp, no copyright line and no "Available online" line anywhere in the three pages; the file is a 2006 scan (`Creator: 0M6127010.TIF`, `Producer: PageGenie PDFGenerator`). Cite by volume/page only |
| Naming | "Rel %" = relative percent of the **total peak area of the chromatogram**, all peaks including the methyl esters and the unassigned minor peaks (it sums to 45.9 + 47.7 + 6.4 = **100.0**). "Normalized" = the same value renormalised to the eleven **non-ester** volatiles only (their 45.9 rescaled to 100). "Yield %" in Table II is **not** a yield — see flag 1. "8-, 9-, 10-, 11-OOH" are hydroperoxide positions on the oleate chain; oleate has one double bond, so there is no ct/tt geometry axis here as there is in the linoleate papers |
| Its own forward reference | The last sentence promises a full paper "that will deal also with the decomposition of linoleate and linolenate hydroperoxides" — **that full paper is `data/articles/frankel1981.pdf`**, also on disk, also in this wave |
| Repo status before this dossier | **Not cited anywhere in `src/kinetic_core/`.** `species_lipid.LOOH_OL` says its branch fractions are "MEASURED BY NOTHING IN THE FIT CORPUS"; that statement was true of the corpus as constituted and is what this dossier changes |

## 1. Why it matters

**The engine refuses `nonanal` by name, and this paper is the missing measurement.**
`src/kinetic_core/lipid.py` line 554 raises a refusal whose text is:

> "REFUSED. This matrix carries an OLEATE hydroperoxide pool (… mmol/L), and the oleate -> nonanal
> branch fraction is measured NOWHERE in the fit corpus. Frankel 1989 fed linoleate only, and
> nonanal appears in no table, figure or sentence of it."

Both halves of that sentence are correct about Frankel 1989 and the first half is **no longer
correct about the corpus**. Selke 1978 feeds *pure methyl oleate hydroperoxides* — the exact parent
of `species_lipid.LOOH_OL`, right down to its lumped composition, since `LOOH_OL` is documented as
the "oleate hydroperoxide pool (8-/9-/10-/11-OOH, lumped)" and this paper's feed is 27 % 8-, 23 %
9-, 23 % 10-, 27 % 11-OOH — and prints nonanal at **15 % of total volatile peak area**, the single
largest non-ester peak in the chromatogram. `species_lipid.NONANAL`'s docstring says nonanal "has
exactly one incoming edge in this network, from `LOOH_OL`, whose branch fraction is unmeasured".
That edge now has a number attached to it.

**But the number is not the same kind of number as `FRANKEL_ZERO_ADDITIVE`, and the difference is
the whole finding.** Frankel 1989's shares are percent of a **six-peak sum**; Selke's Rel % is
percent of the **whole chromatogram** (19 identified peaks plus a 6.4 % unassigned remainder), and
his "Normalized" column is percent of the **eleven non-ester** peaks. Three different denominators
across two papers from the same bench. Dropping Selke's 15 straight into a slate normalised the
Frankel way would be a category error; section 4 states the three candidate forms explicitly and
says which one a `branch_fraction` row can honestly carry.

**Relation to `frankel1989_extraction.md`, stated explicitly as required.**

| axis | Selke 1978 (this paper) | Frankel & Gardner 1989 |
|---|---|---|
| laboratory | Northern Regional Research Center, USDA, Peoria IL | **the same** Northern Regional Research Center, USDA, Peoria IL |
| shared author | **E. N. Frankel** | **E. N. Frankel** (with H. W. Gardner) |
| method | injector-port thermolysis, GC-MS ("reaction chromatography") | **the same** injector-port thermolysis, GC |
| substrate | methyl **oleate** hydroperoxides, 8/9/10/11-OOH | methyl **linoleate** hydroperoxides, 9/13-OOH |
| **injector temperature** | **210 C** | **180 C** — a 30 C gap, uncalibrated in either paper |
| sample form | **neat**, 3.8 µL | **hexane solution**, 1 µL of a 200 µL solution |
| trapping | GC program started at **25 C** ("instead of -60 C") | column head **cryo-trapped at -65 C** |
| hydroperoxide source | autoxidation in O2 at **40 C** to PV **1051**, column partition chromatography | autoxidation with O2 at **40 C**, silicic acid chromatography (plus HPLC and lipoxygenase preparations) |
| quantification | Rel % of **total** peak area; no internal standard named | % of a **six-peak** sum; internal standard methyl hexanoate |
| replication | **none stated** — one injection, no RSD, no duplicates | duplicate GC analyses, RSD ±3.9 to ±4.8 % |
| rate, Ea, absolute yield | **none** | **none** |

So: **same laboratory, same method, non-overlapping substrates, and only partly overlapping
products.** They are complementary, not redundant, and there is exactly one product in common to
cross-check on — see the next paragraph.

**The one cross-check that exists, and what it says.** `ME_9_OXONONANOATE` is in both slates.
Frankel 1989 makes it from the **linoleate** 9-hydroperoxide (both homolytic pathway B and the Hock
route) at 13 / 4.3 / 26 % of his six-peak sums; Selke 1978 makes it from the **oleate** 9- and
10-hydroperoxides at 15 % of total peak area, the joint-largest ester peak. **The two are not the
same quantity and cannot be compared numerically** (different denominators, different substrates,
different injector temperatures). What they establish *qualitatively* is a real and unmodelled
structural fact: **methyl 9-oxononanoate has two parents in any matrix that carries both oleate and
linoleate lipid** — `LOOH_OL` and the linoleate 9-hydroperoxide pool — while
`src/kinetic_core/species_lipid.py` gives it only the linoleate edge. Any real food fat is mostly
oleate. This is a mis-specification flag for the lane, independent of the nonanal question, and it
is the concrete answer to "do the product distributions agree?": **they cannot be made to agree,
because no printed percentage in the two papers shares a denominator.**

**And `2-pentylfuran` is not here.** This paper does not name it, which is structurally correct —
2-pentylfuran is a linoleate product and this feed is pure oleate. The alkylfuran refusal is
untouched by this paper; `frankel1981.pdf` is where to look.

## 2. Methods as they matter to a model

- **The pot.** There is no pot. This is **reaction chromatography**: a neat sample of hydroperoxides
  is injected into a hot injector port, decomposes there, and the volatiles are swept onto the
  column. There is **no solvent, no reaction time, no reactor volume, no conversion, no
  concentration at any moment and no atmosphere other than the GC carrier gas.**
- **The substrate and how it was made.** "The hydroperoxides from methyl oleate **autoxidized in O2
  at 40 C to a peroxide value of 1051** were purified by **column partition chromatography** (6). The
  hydroperoxides (**checked for purity by thin layer chromatography**) were analyzed for isomeric
  composition by GC-MS (7): **27 % 8-, 23 % 9-, 23 % 10-, and 27 % 11-OOH isomers.**"
- **Exactly which isomers were isolated: none of them individually.** Unlike Frankel 1989, which ran
  three separate preparations (mixed ct/tt, HPLC-separated tt, and a pure lipoxygenase ct-13),
  **Selke 1978 ran ONE pool containing all four positional isomers at once** and disentangled them
  arithmetically afterwards, by matching each product to the side of the chain it must have come
  from (Table II). There is no isomer-resolved experiment in this paper. **Every per-isomer number in
  Table II is an assignment of a measured whole-pool peak, not a separate measurement.**
- **The thermolysis.** "A **neat sample of oleate hydroperoxides (3.8 µL)** was injected into the
  same GC-MS system used previously (1). The GC parameters were similar except **the injector port
  temperature was 210 C**, and temperature programming was **initiated at 25 C instead of -60 C**."
  So: **210 C injector, no cryo-trap**, in contrast to Frankel 1989's 180 C injector with a -65 C
  cryo-trap. The two most volatile classes in an oleate slate (heptane, octane) are exactly the ones
  a 25 C start handles worse than a -60 C start; see flag 5.
- **Identification.** "Identifications of volatile compounds were based on **mass spectra matched
  manually and by computer with those of reference compounds** and were **confirmed by GC-retention
  data**." Two products carry footnote *b*: "**Tentative identification based on GC elution and MS
  without reference compound**" — methyl 10-oxo-8-decenoate and methyl 11-oxo-9-undecenoate.
- **Quantification.** "volatiles were separated, identified, and their **relative proportion
  estimated by direct analysis of micro samples**." That is the entire statement of method.
  **No internal standard is named, no response factors are given, no replicates are reported and no
  error bar of any kind appears in the paper.** The Rel % column sums to 100.0 over the whole
  chromatogram, which is what fixes its denominator.
- **The comparison arm.** The Triolein columns of Table I are **not measured in this paper**; they
  are carried over from ref (1) = Selke, Rohwedder & Dutton, *JAOCS* **54**:62 (1977), triolein
  heated in air at **192 C**. They are `[C]` cited throughout.
- **The mechanism, as the authors state it.** "The well-recognized mechanism of **carbon-carbon
  scission on either side of the alkoxy radical intermediate** produced from hydroperoxides (8) was
  checked by matching the concentration of cleavage products expected from each part of the oleate
  hydroperoxide isomers (Table II)." Side **A** is the cleavage on the carboxyl side, side **B** on
  the methyl side (the labels are printed over the structure in Table II).
- **One extra assumption the authors add, which the model should know about.** "**We assumed further
  that 1-enols are produced from the reaction of hydroxy radicals with 1-olefins to form the
  corresponding saturated aldehydes by tautomerism.** For example, decanal would be formed as
  follows: CH3(CH2)7CH=CH· + ·OH → CH3(CH2)7CH=CH-OH → CH3(CH2)7CH2CHO". **The saturated aldehydes
  in this slate — decanal, nonanal, octanal — are therefore attributed to a route that needs a
  hydroxyl radical and a tautomerisation, not to plain beta-scission.** Nonanal, the compound the
  repository wants, is one of these. This is the authors' hypothesis, not a measurement.
- **What the paper claims it proves.** "These results clearly indicate that **oleate hydroperoxides
  are the major precursors of volatiles produced from triolein even at 192 C**. However, these data
  are **insufficient to prove that hydroperoxidation is the only route** by which such products
  form."
- **What the paper says about kinetics — the sentence that matters most for the lane's Q10
  assumption.** Introduction: hydroperoxides "are **rapidly decomposed at temperatures exceeding
  100 C** (2), producing **90 % polymeric and 10 % volatile materials** (3) … If hydroperoxides are
  formed as intermediates at high temperatures, **there is no information about their finite
  existence and about the kinetics that control their decomposition.**" Both the 90/10 split and the
  ">100 C" claim are **cited to 1960/1961 sources not on disk**, and the kinetic statement is an
  explicit declaration of absence. This paper adds **no rate and no activation energy**, exactly as
  Frankel 1989 adds none.

## 3. Tables re-typed

Marks: `[M]` measured in this paper, `[C]` cited from another paper, `[F]` fitted/derived by the
authors from their own measured values.

### TABLE I (p. 512). "Comparison of Volatiles from Decomposed Oleate Hydroperoxides and Heated Triolein"

Footnotes as printed: *a* "See Figure 2 for numbered peaks." *b* "Tentative identification based on
GC elution and MS without reference compound."

| Peak no.*a* | Compound | Oleate-hydroperoxides Rel % | Oleate-hydroperoxides Normalized | Triolein (1) Rel % | Triolein (1) Normalized |
|---:|---|---:|---:|---:|---:|
| 1 | Heptane | 4.4 `[M]` | 9.6 `[F]` | 8.6 `[C]` | 9.5 `[C]` |
| 2 | Octane | 2.7 `[M]` | 5.9 `[F]` | 9.7 `[C]` | 10.8 `[C]` |
| 3 | Heptanal | 0.5 `[M]` | 1.1 `[F]` | 5.1 `[C]` | 5.6 `[C]` |
| 4 | 1-Heptanol | 0.4 `[M]` | 0.9 `[F]` | 1.6 `[C]` | 1.8 `[C]` |
| 5 | Octanal | 11 `[M]` | 23.9 `[F]` | 8.5 `[C]` | 9.4 `[C]` |
| 7 | 1-Octanol | 0.4 `[M]` | 0.9 `[F]` | 2.5 `[C]` | 2.4 `[C]` |
| **8** | **Nonanal** | **15** `[M]` | **32.6** `[F]` | **22** `[C]` | **24.3** `[C]` |
| 10 | 2-Nonenal | 0.5 `[M]` | 1.1 `[F]` | 2.0 `[C]` | 2.2 `[C]` |
| 11 | Decanal | 3.9 `[M]` | 8.5 `[F]` | 2.8 `[C]` | 3.1 `[C]` |
| 13 | c/t-2-Decenal | 5.4 `[M]` | 11.8 `[F]` | 17 `[C]` | 18.8 `[C]` |
| 14 | t-2-Undecenal | 1.7 `[M]` | 3.7 `[F]` | 11 `[C]` | 12.1 `[C]` |
| | **(non-ester subtotal, printed)** | **45.9** | **100** | **90.8** | **100** |
| 6 | Me heptanoate | 1.5 `[M]` | — | — | — |
| 9 | Me octanoate | 5.0 `[M]` | — | — | — |
| 12 | Me nonanoate | 1.5 `[M]` | — | — | — |
| 15 | Me 8-oxooctanoate | 3.5 `[M]` | — | — | — |
| 16 | Me 9-oxononanoate | 15 `[M]` | — | — | — |
| 17 | Me 10-oxodecanoate | 12 `[M]` | — | — | — |
| 18 | Me 10-oxo-8-decenoate*b* | 3.4 `[M]` | — | — | — |
| 19 | Me 11-oxo-9-undecenoate*b* | 5.8 `[M]` | — | — | — |
| | **(ester subtotal, printed)** | **47.7** | — | — | — |
| | Other minor peaks | 6.4 `[M]` | — | — | — |

**Arithmetic checks (mine).** Non-esters sum to 4.4+2.7+0.5+0.4+11+0.4+15+0.5+3.9+5.4+1.7 = **45.9**,
matching the printed subtotal exactly. Esters sum to 1.5+5.0+1.5+3.5+15+12+3.4+5.8 = **47.7**,
matching exactly. 45.9 + 47.7 = **93.6**, which is the figure the text prints ("These compounds,
together with the methyl ester fragments, represent **93.6 % of the relative total peak area of
Figure 1**"), and 93.6 + 6.4 = **100.0**. The Normalized column reproduces as x/45.9 to the printed
digit (15/45.9 = 32.68 → 32.6; 11/45.9 = 23.97 → 23.9; 4.4/45.9 = 9.59 → 9.6). The Triolein Rel %
column sums to **90.8** and its Normalized column reproduces as x/90.8 (22/90.8 = 24.2 → 24.3;
11/90.8 = 12.1). **The corrected table closes on all four printed subtotals; the OCR text layer's
version closes on none of them, which is how the OCR errors were caught.**

### TABLE II (p. 513). "Decomposition of Methyl Oleate Hydroperoxide Isomers"

Column head: "Cleavage products" | "Yield %". The parenthetical percentages are the same Rel %
values as Table I, re-quoted; the Yield % column is their arithmetic sum. Footnote *a*, printed
once at the foot: "**Because nonanal and 9-oxononanoate arise from both the 9- and
10-hydroperoxides, the concentrations of these compounds were divided by assuming that the same
amount of saturated aldehydes would be produced from the 9- and from the 8-hydroperoxides.**"

Each structure is printed with the scission points marked **B** (methyl side) and **A** (carboxyl
side) and the carbon bearing the O· labelled with the isomer number.

**8-hydroperoxide** — CH3(CH2)7-CH=CH ⟦B⟧ CH(8)(O·) ⟦A⟧ (CH2)6COOMe

| side | cleavage products | Yield % |
|---|---|---:|
| A | 2-Undecenal (1.7 %) + Me heptanoate (1.5 %) | 3.2 `[F]` |
| B | Decanal (3.9 %) + Me 8-oxooctanoate (3.5 %) | 7.4 `[F]` |
| | **subtotal** | **10.6** |

**9-hydroperoxide** — CH3(CH2)6-CH=CH ⟦B⟧ CH(9)(O·) ⟦A⟧ (CH2)7COOMe

| side | cleavage products | Yield % |
|---|---|---:|
| A | 2-Decenal (5.4 %) + Me octanoate (5.0 %) | 10.4 `[F]` |
| B | **Nonanal (4.0 %)*a*** + Me 9-oxononanoate (4.0 %)*a* | 8.0 `[F]` |
| | **subtotal** | **18.4** |

**10-hydroperoxide** — CH3(CH2)7 ⟦B⟧ CH(10)(O·) ⟦A⟧ CH=CH-(CH2)6COOMe

| side | cleavage products | Yield % |
|---|---|---:|
| A | **Nonanal (11 %)*a*** + Me 9-oxononanoate (11 %)*a* | 22.0 `[F]` |
| B | Octane (2.7 %) + Me 10-oxo-8-decenoate (3.4 %) + 1-octanol (0.4 %) | 6.5 `[F]` |
| | **subtotal** | **28.5** |

**11-hydroperoxide** — CH3(CH2)6 ⟦B⟧ CH(11)(O·) ⟦A⟧ CH=CH-(CH2)7COOMe

| side | cleavage products | Yield % |
|---|---|---:|
| A | Octanal (11 %) + Me 10-oxodecanoate (12 %) | 23.0 `[F]` |
| B | Heptane (4.4 %) + Me 11-oxo-9-undecenoate (5.8 %) + 1-heptanol (0.4 %) | 10.6 `[F]` |
| | **subtotal** | **33.6** |

### Statements printed in the running text (p. 512-513), transcribed

- "Eleven of these peaks are due to the same compounds previously identified from heated triolein
  (Table I)." `[M]`
- "**These compounds, together with the methyl ester fragments, represent 93.6 % of the relative
  total peak area of Figure 1.** The remaining peaks are due to minor components which were too
  small to identify reliably." `[M]`
- "**The yields of cleavage products arising from each side of the hydroperoxide isomers were in
  remarkably good agreement.**" (the authors' own reading of Table II) `[F]`
- "This mechanism accounts for all the volatiles listed in Table I **except for heptanal,
  2-nonenal, and methyl nonanoate**." `[M]`
- "Although the isomeric composition of the oleate hydroperoxides was symmetrical (e.g., 8- =
  11-OOH, 9- = 10-OOH), the **total relative concentration of volatiles from the 10- and
  11-hydroperoxides was 61 %** and that from the **8- and 9-hydroperoxides was 29 %**." `[F]`
- "This divergence is apparently related to the position of the hydroperoxide group and allylic
  unsaturation. **Cleavages yielding octanal, nonanal, and decanal seemed favored kinetically in the
  oleate hydroperoxides.**" `[F]` — this is the paper's only use of the word "kinetically" and it is
  a statement about product ratios, not about a rate.
- "The implication of these results will be discussed in the full paper that will deal also with the
  decomposition of linoleate and linolenate hydroperoxides." (→ Frankel, Neff & Selke 1981)

### Derived numbers (mine, arithmetic on the printed tables — NOT the paper's)

- **Table II's four subtotals sum to 10.6 + 18.4 + 28.5 + 33.6 = 91.1** (mine). Table I's identified
  peaks total 93.6, and the three products the mechanism does not explain are heptanal 0.5 +
  2-nonenal 0.5 + Me nonanoate 1.5 = **2.5** (mine). 91.1 + 2.5 = **93.6** exactly. The accounting
  is internally airtight.
- **The text's "61 %" does not reproduce.** 28.5 + 33.6 = **62.1** (mine), not 61; the 8-/9- pair
  reproduces exactly (10.6 + 18.4 = **29.0** vs the printed 29 %). Recorded as a printed
  inconsistency, not corrected. See flag 4.
- **The A/B pairings the mechanism implies, as ratios** (mine), all from Table II's own
  parentheticals: 2-undecenal/Me heptanoate = 1.7/1.5 = **1.13**; decanal/Me 8-oxooctanoate =
  3.9/3.5 = **1.11**; 2-decenal/Me octanoate = 5.4/5.0 = **1.08**; octanal/Me 10-oxodecanoate =
  11/12 = **0.92**. All four fall within 1.13/0.92 = **1.23x** of each other. That is what the
  authors mean by "remarkably good agreement", and it is a much tighter closure than the **8.5x**
  span Frankel 1989's linoleate pairings showed (`frankel1989_extraction.md` section 3). Two of the
  four Selke pairings involve a C10-C11 aldehyde against a C7-C11 ester — chain lengths far closer
  than Frankel's C5 alkane against a C14 dienoate — so the contrast is at least partly a detector-
  commensurability artefact and should not be read as "the oleate mechanism is cleaner".
- **Nonanal expressed three ways** (mine, from the printed 15): **15.0 %** of total chromatogram
  peak area; **32.6 %** of non-ester volatiles (printed); **16.0 %** of identified peak area
  (15/93.6). Section 4 says which of these can be used as what.
- **Nonanal + its scission partner methyl 9-oxononanoate = 15 + 15 = 30 % of total peak area**
  (mine) — i.e. **the nonanal channel and its ester half together are just under a third of the
  entire volatile chromatogram**, and the strict 1:1 stoichiometry the mechanism demands is met
  exactly (15/15 = 1.00, mine), the only exact pairing in either Frankel-lab paper.

## 4. Numbers the repository can use

**All rows share these conditions unless stated: neat methyl oleate hydroperoxides (27 % 8-, 23 %
9-, 23 % 10-, 27 % 11-OOH; from autoxidation in O2 at 40 C to PV 1051), 3.8 µL injected, thermolysed
in a GC injector port at 210 C, GC program from 25 C, GC-MS identification, quantified as relative
peak area with no response factors, no internal standard, no replicates and no stated error.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **oleate LOOH → nonanal, share of the whole chromatogram** | **15** | % of total volatile peak area (denominator = 100.0, all peaks incl. esters and 6.4 % unassigned) | 210 C injector-port thermolysis, whole-pool oleate LOOH | Table I, peak 8, p. 512 | **branch_fraction** — the primary candidate for the refused `LOOH_OL → NONANAL` edge |
| **oleate LOOH → nonanal, share of the non-ester volatiles** | **32.6** | % of the eleven non-ester volatiles (denominator = 45.9) | as above | Table I, peak 8 "Normalized", p. 512 | **within_study_ratio** — a *different* branch fraction with a *different* denominator; do not mix with the row above |
| oleate LOOH → nonanal, share of identified peaks | **16.0** (mine, 15/93.6) | % of identified peak area | as above | derived from Table I | derived_assumption — my arithmetic, not printed |
| nonanal attributed to the 9-hydroperoxide | 4.0 | % of total peak area | as above | Table II, 9-OOH side B, p. 513 | **derived_assumption** — footnote *a* says this split was **assumed**, not measured (nonanal from 9-OOH was set equal to the saturated aldehyde from 8-OOH); only the total 15 is measured |
| nonanal attributed to the 10-hydroperoxide | 11 | % of total peak area | as above | Table II, 10-OOH side A, p. 513 | **derived_assumption** — the remainder of the same assumed split |
| **full oleate product slate, share of total peak area** | heptane 4.4 · octane 2.7 · heptanal 0.5 · 1-heptanol 0.4 · octanal 11 · 1-octanol 0.4 · **nonanal 15** · 2-nonenal 0.5 · decanal 3.9 · c/t-2-decenal 5.4 · t-2-undecenal 1.7 · Me heptanoate 1.5 · Me octanoate 5.0 · Me nonanoate 1.5 · Me 8-oxooctanoate 3.5 · **Me 9-oxononanoate 15** · Me 10-oxodecanoate 12 · Me 10-oxo-8-decenoate 3.4 · Me 11-oxo-9-undecenoate 5.8 · other minor 6.4 | % of total peak area, summing to 100.0 | as above | Table I, p. 512 | **branch_fraction** (the whole distribution) — the oleate analogue of `FRANKEL_ZERO_ADDITIVE`, but on a whole-chromatogram denominator |
| octanal share | 11 (23.9 normalized) | % | as above | Table I, peak 5 | branch_fraction |
| decanal share | 3.9 (8.5 normalized) | % | as above | Table I, peak 11 | branch_fraction |
| c/t-2-decenal share | 5.4 (11.8 normalized) | % | as above | Table I, peak 13 | branch_fraction |
| t-2-undecenal share | 1.7 (3.7 normalized) | % | as above | Table I, peak 14 | branch_fraction |
| 2-nonenal share | 0.5 (1.1 normalized) | % | as above | Table I, peak 10 | branch_fraction — **the compound Frankel 1989 names in his introduction and never measures**, measured here but from the wrong parent (oleate, not linoleate) |
| methyl 9-oxononanoate share | 15 | % of total peak area | as above | Table I, peak 16 | branch_fraction — **the only product shared with Frankel 1989's slate; see flag 3** |
| nonanal : methyl 9-oxononanoate | **1.00** (mine, 15/15) | — | as above | derived from Table I | within_study_ratio — the exact 1:1 scission pairing |
| the four A/B pairing ratios | 1.13 / 1.11 / 1.08 / 0.92 (mine) | — | as above | derived from Table II | within_study_ratio |
| products of the 10- + 11-hydroperoxides | **62.1** (mine, from Table II) / **61** (printed in text) | % of total peak area | as above | Table II subtotals; text p. 513 | within_study_ratio — **the two disagree, see flag 4** |
| products of the 8- + 9-hydroperoxides | **29.0** (mine) / **29** (printed) | % of total peak area | as above | Table II subtotals; text p. 513 | within_study_ratio — these agree |
| identified fraction of the chromatogram | 93.6 | % of total peak area | as above | text p. 512 | measured_bound — **6.4 % of the chromatogram is unassigned**, the direct analogue of `LIPID_FRAG_C` |
| unexplained-by-mechanism fraction | 2.5 (mine: heptanal 0.5 + 2-nonenal 0.5 + Me nonanoate 1.5) | % of total peak area | as above | text p. 513 + Table I | measured_bound |
| feed isomer composition | 27 / 23 / 23 / 27 | % 8- / 9- / 10- / 11-OOH | the pool that was thermolysed | Experimental, p. 512 | measured_ratio — **directly parameterises `LOOH_OL`'s lumped composition** |
| peroxide value of the feed | 1051 | (meq/kg, unit not printed) | methyl oleate autoxidised in O2 at 40 C | Experimental, p. 512 | level_only — the unit is not stated in the paper |
| volatile vs polymeric split of decomposing hydroperoxide | **90 % polymeric, 10 % volatile** | % of decomposed hydroperoxide | ">100 C", generic heated fat | Introduction, p. 511, **cited to ref (3)** = Evans, Flavor Chemistry Symposium, Campbell Soup Co., 1961 | **`[C]` cited, not measured here** — and it is precisely the factor that would convert every peak-area share above into an absolute yield per hydroperoxide. **It is second-hand, from a 1961 conference volume not on disk. Do not use it to close the mass balance.** |
| triolein comparison slate | heptane 8.6 · octane 9.7 · heptanal 5.1 · 1-heptanol 1.6 · octanal 8.5 · 1-octanol 2.5 · **nonanal 22** · 2-nonenal 2.0 · decanal 2.8 · c/t-2-decenal 17 · t-2-undecenal 11 (sum 90.8) | % of that paper's total peak area | **triolein heated in air at 192 C** | Table I, right columns, `[C]` from Selke, Rohwedder & Dutton, *JAOCS* 54:62 (1977) | **`[C]` cited from a paper NOT on disk** — a genuine second-system nonanal number (22 %, 24.3 % normalized) but not measured here |
| **rate constant at any temperature** | **NOT PRESENT** | — | — | — | — |
| **activation energy** | **NOT PRESENT** | — | — | — | — |
| **absolute yield (mol/mol, mass, mmol/L)** | **NOT PRESENT** | — | — | — | — |
| **2-pentylfuran** | **NOT PRESENT** — not named anywhere | — | — | — | structurally correct: pure oleate feed |
| **hexanal, pentane, 2,4-decadienal, methyl 13-oxo-tridecadienoate** | **NOT PRESENT** | — | — | — | structurally correct: all four are linoleate products |

### **WHICH KIND OF QUANTITY IS "15 %"? — the question the brief demands be answered**

**It is a relative GC peak area, expressed as a percent of the total peak area of the whole
chromatogram.** It is **not** a fraction of the hydroperoxide consumed, and it is **not** a percent
of total volatiles in the sense of a mass or molar balance. Specifically:

1. **NOT a fraction of the hydroperoxide consumed.** Nothing in this paper measures how much
   hydroperoxide decomposed, how much survived, or how much went to non-volatile products. The
   paper's own introduction quotes a *cited* 90 % polymeric / 10 % volatile split for heated fats,
   which — if it applied here, and there is no evidence it does — would mean the entire 100 % of
   this chromatogram is a share of the **10 %** that volatilised. **A branch fraction per
   hydroperoxide consumed cannot be built from this paper without importing that 1961 number.**
2. **It IS a share of total volatiles, in the peak-area sense, because the denominator closes to
   100.** 45.9 + 47.7 + 6.4 = 100.0. That is stronger than Frankel 1989's six-peak sum, which is a
   share of a *selected* subset. **Selke's denominator includes an explicit 6.4 % unassigned
   remainder**, which is the honest thing and is what `LIPID_FRAG_C` exists to represent.
3. **It is an AREA share, not a molar or mass share.** No response factors, no internal standard, no
   calibration curve. Over a slate running from C7 alkane to a C11 oxo-ester on a 1978 GC-MS with a
   25 C start, area is not moles. **Do not convert 15 % to mmol/L.**
4. **The "Normalized" column (32.6) is a fourth thing again** — a share of the eleven non-ester
   volatiles only, which excludes 47.7 % of the chromatogram. It exists to make the oleate slate
   comparable to the triolein slate of ref (1), which reported only non-esters. **Its only correct
   use is that comparison.**
5. **Table II's "Yield %" is none of the above.** It is the arithmetic sum of Table I Rel % values
   for the products assigned to one cleavage side. Calling it a yield is the paper's word, not a
   measurement; the column head is misleading and section 5 flags it.

**Recommended form for a `branch_fraction` row.** `oleate_looh_to_nonanal = 0.15`, denominator
= *total volatile peak area including unassigned*, method = `injector_port_thermolysis_210C`,
`response_factor_corrected = False`, `n = 1`, `error = none stated`. If the lane instead wants a
Frankel-commensurable share, the honest construction is **not** to renormalise Selke onto a
six-peak sum — the two slates share only one product — but to carry the oleate slate as its own
distribution with its own denominator, and to make the denominator a first-class field.

## 5. Flags

1. **Table II's column head says "Yield %" and the column is not a yield.** Every entry is the sum
   of two or three Table I relative peak-area percentages. There is no conversion, no response
   factor, no mass balance and no measurement of hydroperoxide consumed anywhere in the paper. A
   reader (or a scraper) taking "Yield %" at face value would import a 22.0 % *molar yield of
   nonanal from the 10-hydroperoxide* that does not exist.
2. **n = 1, no replicates, no error bar, no internal standard.** The paper reports a single neat
   3.8 µL injection. Frankel 1989 at least prints an RSD of ±3.9-4.8 % on duplicate GC analyses;
   **this paper prints nothing.** Every number in it is a single-shot value from a 1978 GC-MS. Any
   uncertainty band placed on the 15 % is invented.
3. **The one product shared with Frankel 1989 cannot be numerically cross-checked.** Methyl
   9-oxononanoate is 15 % here (of a whole chromatogram, from oleate 9- and 10-OOH at 210 C) and
   13 / 4.3 / 26 % there (of a six-peak sum, from linoleate 9-OOH at 180 C). Different denominator,
   different substrate, different temperature, different trapping. **The papers agree
   mechanistically and are numerically incomparable.** Do not report a "consistency check" between
   them. The one transferable finding is structural: in any oleate+linoleate matrix, methyl
   9-oxononanoate has two parents and the repository models one.
4. **A printed internal inconsistency.** The text says volatiles from the 10- and 11-hydroperoxides
   were "61 %"; Table II's own subtotals give 28.5 + 33.6 = 62.1. The companion figure (29 % for
   8-/9-) reproduces exactly. Recorded as printed. It is one point out of 62 and changes nothing
   about nonanal, but it means the text's summary statistics were not machine-checked against the
   table.
5. **210 C is not 180 C, and the trapping differs too.** Selke ran the injector 30 C hotter than
   Frankel 1989 and started the column at **25 C instead of -60 C**. The paper says so explicitly.
   The most volatile peaks in an oleate slate (heptane b.p. 98 C, octane b.p. 126 C) are the ones a
   25 C start recovers worst, so **the light end of this slate is plausibly under-counted relative
   to Frankel's cryo-trapped runs**, which inflates every heavier share including nonanal's by an
   unknown amount. There is no way to correct for this from the paper.
6. **The saturated aldehydes — nonanal included — are attributed to a hypothesised route.** The
   authors' 1-enol/hydroxyl-radical tautomerisation mechanism ("We assumed further that…") is how
   decanal, nonanal and octanal are explained. If that route needs an ·OH population, then the
   nonanal share depends on radical chemistry that a neat 210 C injector port supplies and a food
   matrix may not. **A 15 % branch fraction transferred to 150 C dough is transferred across a
   mechanism the authors themselves flagged as an assumption.**
7. **Two products are tentatively identified.** Methyl 10-oxo-8-decenoate (3.4 %) and methyl
   11-oxo-9-undecenoate (5.8 %) carry footnote *b*, "without reference compound" — together 9.2 % of
   the chromatogram. Nonanal is **not** among these; it was confirmed against a reference compound
   and by retention data.
8. **A dangling figure reference.** Table I footnote *a* reads "See Figure 2 for numbered peaks",
   but the communication prints only **Figure 1** (the chromatogram, whose caption is "Gas
   chromatogram of volatile compounds from thermally decomposed methyl oleate hydroperoxides"). The
   peak numbers used in Table I are the ones annotated on Figure 1. Harmless, but it means the peak
   numbering has no independent printed key.
9. **The 90/10 volatile/polymeric split is second-hand and 17 years older than the paper.** It is
   cited to a 1961 Campbell Soup Co. symposium volume, for "heated fats" generically, at
   temperatures ">100 C", with no substrate specified. It is the only bridge in existence between
   this paper's shares and an absolute yield, and it is not strong enough to carry that weight. **Do
   not use it.**
10. **The triolein arm is not this paper's data.** The right-hand half of Table I is cited from
    *JAOCS* 54:62 (1977), a paper **not on disk**. Its nonanal value (22 % rel, 24.3 % normalized,
    triolein heated in air at **192 C**) is a genuinely independent second system and would be worth
    retrieving — but as printed here it is a citation, and treating it as a replicate of the 15 %
    would be double-counting a number nobody in this corpus has verified.
11. **The whole pool was decomposed at once; the per-isomer numbers are bookkeeping.** Unlike
    Frankel 1989's three physically separated preparations, this paper ran one mixed pool. Table
    II's per-isomer assignment is an inference from product identity, and the nonanal split in
    particular (4.0 / 11) rests on footnote *a*'s explicit assumption. **`LOOH_OL`'s lumped
    treatment is therefore the *right* level of resolution for this source — a per-isomer oleate
    branch model would be fitting the authors' assumption, not their data.**
12. **No DOI exists to cite.** See section 0. Anything in the repository that requires a DOI field
    for this reference must record it as absent rather than reconstructing one.
13. **No aqueous phase, no amine, no pH, no water activity, no matrix.** Neat hydroperoxide in a
    hot metal port. This paper contributes nothing to the aldehyde-lysine channel, nothing to
    matrix retention, and nothing to any pH-dependent term.
14. **What to request.** (i) Selke, Rohwedder & Dutton, *JAOCS* **54**:62 (1977) — the triolein-at-
    192 C paper whose nonanal number is quoted here, which would give a second, independent oleate-
    system nonanal share in a real triglyceride rather than a methyl ester. (ii) Frankel, Evans,
    McConnell & Jones, *JAOCS* **38**:134 (1961), ref (6), the purification method — only if the
    isomer composition needs auditing. (iii) The same experiment with response factors and a mass
    balance, which is the same request `frankel1989_extraction.md` already makes and which no paper
    in this corpus satisfies for either substrate.

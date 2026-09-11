# Frankel & Gardner 1989 — EXTRACTION (methyl linoleate hydroperoxides thermolysed in a GC injector port at 180 C, with alpha-tocopherol or 1,4-cyclohexadiene as hydrogen donors; three hydroperoxide preparations, six volatile products, 26 distribution columns)

### The whole product slate of the repository's lipid lane, printed as relative percent of a six-peak sum at one temperature, with the total peak area beside it — and with no rate, no yield, no barrier, no nonanal and no 2-pentylfuran anywhere in the paper.

**Source on disk:** `data/articles/frankel1989.pdf` (6 pp., Lipids 24 (7), 603-608, 1989). The
`pdftotext -layout` text layer is an OCR layer with spaced-out glyphs in the justified prose
("t h e r m a l") and mangled Greek ("a-tocopherol", "~-tocopherol"), but **Tables 1 and 2 came
through with every cell and every column head intact and were verified against 190-dpi rasters of
printed pages 604 and 606** (`scratchpad/img/fr89t1b-2.png`, `fr89t1c-2.png`, `fr89t2b-4.png`);
every number below matches the raster. Figures 1-4 (the four regression plots) and Scheme 1 (the
cleavage mechanism) are **FIGURE-ONLY**; their captions carry printed statistics, which are
transcribed as statistics, and no data point is read off any plot.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of a-Tocopherol on the Volatile Thermal Decomposition Products of Methyl Linoleate Hydroperoxides" |
| Authors | **E. N. Frankel** (corresponding) and **H. W. Gardner** — two authors, not "et al." — Northern Regional Research Center, Agricultural Research Service, USDA, Peoria IL |
| Venue | Lipids 24 (7), 603-608 (1989). Received February 10 1989; revision accepted April 28 1989 |
| Naming | "ct" / "cis,trans" and "tt" / "trans,trans" are hydroperoxide **double-bond geometries**; "9-" and "13-" are the **hydroperoxide positions**. "Me octanoate" = methyl octanoate; "Me 9-oxononanoate" and "Me 13-oxo-9,11-tridecadienoate" are the methyl esters. "Relative percent" always means percent of the sum of the six major volatiles in that one column |
| Where the repository already cites it | `src/kinetic_core/parameters_lipid.py` (`FRANKEL_ZERO_ADDITIVE`, `FRANKEL_SYSTEM_GEOMETRY`, `FRANKEL_ZERO_ADDITIVE_TOTAL_AREA`, `FRANKEL_ANCHOR`); `src/kinetic_core/species_lipid.py` (`LIPID_PRODUCTS`, `FRANKEL_SLATE`, `POSITION_PRODUCTS`, `CLEAVAGE_MECHANISM`, `NONANAL`, `NAMED_UNQUANTIFIED_COPRODUCTS`); `results/validation/kinetic_core_b6_prereg.md`; `docs/reference/FIT_HOLDOUT_DECLARATION.md` D.6 Module 5 |
| Documented until now only through | `k3_final_parameter_inventory.md` sec. A.5 and sec. C.9 — the `FRANKEL_DOSSIER` string points at those, and there has been no standalone dossier |

## 1. Why it matters

**This paper is the entire product slate of `src/kinetic_core/lipid.py`.** Six species — `PENTANE`,
`HEXANAL`, `ME_OCTANOATE`, `DECADIENAL`, `ME_9_OXONONANOATE`, `ME_13_OXO_TRIDECADIENOATE` — exist in
`species_lipid.LIPID_PRODUCTS` because they are the six rows of these two tables, in this order, and
`FRANKEL_SLATE` is literally that tuple. The three zero-additive columns are the lane's eighteen fit
rows and appear as `FRANKEL_ZERO_ADDITIVE` in `parameters_lipid.py`; the alpha-tocopherol columns are
a declared hold-out; the absence of nonanal is a declared negative-test hold-out. The
position-to-product map `POSITION_PRODUCTS` and the homolytic/heterolytic labels in
`CLEAVAGE_MECHANISM` are this paper's introduction and Scheme 1. Nothing else in the corpus supplies
a measured branch distribution for lipid hydroperoxide decomposition.

**And the paper is equally the source of the lane's biggest declared gap.** `parameters_lipid.py`
opens with "The branch DISTRIBUTION is measured. The absolute RATE is not", and cites `k3` sec. C.9
on this paper: "NOT a yield source ... no absolute yield and no Ea exist in it (one temperature,
180 C)". This dossier confirms that verdict against the full text and makes it sharper: the
experiment is **reaction chromatography** — a hexane solution injected into a 180 C injector port,
with the products cryo-trapped on the column head — so there is no reaction time, no reactor volume,
no conversion and no concentration at any moment. A rate cannot be extracted from it even in
principle, and the paper's own last paragraph says so: "We also did not analyze the products of
decomposition of hydroperoxides at lower temperatures than 180 C ... additional studies are needed
with methods permitting the analyses of volatile decomposition products at or below physiological
conditions." The lane's absolute rate therefore comes from Schroen & Berton-Carabin 2022 at 25 C
with a declared Q10 band, which is `kinetic_core_b6_prereg.md`'s stated assumption and remains one.

The hold-out disclosure in `kinetic_core_b6_prereg.md` section 0 is a fact about this paper's
layout, and reading the full text confirms it: **the zero-additive column and the tocopherol columns
are the same table rows**, and the abstract states the hold-out result in prose. There is no way to
read one without the other, which is why B6 scored that hold-out `seen_diagnostic`.

## 2. Methods as they matter to a model

- **Three hydroperoxide preparations, each a different pool composition.**
  1. **Mixed cis,trans / trans,trans 9- and 13-hydroperoxides**: "made by the **autoxidation of pure
     methyl linoleate with oxygen at 40 C** and purified by silicic acid chromatography (16)". This
     is the composition "a real, autoxidised food lipid most nearly resembles", as
     `parameters_lipid.py` puts it.
  2. **trans,trans 9- and 13-hydroperoxides**: separated from that autoxidation mixture by
     semi-preparative reversed-phase HPLC, C-18 column "25.0 X 2.14 mm, 5 microns" (as printed; the
     dimensions are those of a semi-preparative column in cm, see flag 7), 70:30
     acetonitrile:water v/v at 3.0 mL/min, UV 235 nm. "**The cis,trans isomers were eluted between
     40 and 44 min, and the trans,trans isomers between 44 and 52 min.**" Concentrated on a rotary
     evaporator until cloudy, extracted with hexane, the last water removed by azeotropic vacuum
     distillation with absolute acetone.
  3. **Pure methyl cis,trans-13-hydroperoxide**: "prepared by **soy lipoxygenase oxidation of
     linoleic acid**, followed by silicic acid column chromatography (17) and esterification of the
     purified hydroperoxide with **diazomethane**."
- **Additives.** alpha-Tocopherol (Eastman Kodak) and 1,4-cyclohexadiene (Chemical Samples Co.),
  both used as **hydrogen donors, not as antioxidants** — the Discussion is explicit: "a-tocopherol
  and 1,4-cyclohexadiene were used in this study **in relatively large concentrations as hydrogen
  donors to suppress alkoxyl radicals and not as antioxidants, which suppress peroxyl radicals**."
  And: "Relatively high concentrations of a-tocopherol were required in this study to show
  sufficient changes ... Therefore, our results are only relevant to clarifying the mechanism of
  hydroperoxide decomposition and not the antioxidant effect of a-tocopherol."
- **The reaction: thermolysis in a GC injector port ("reaction chromatography").** "a hexane
  solution of hydroperoxides was subjected to thermolysis in the **injector port of a gas
  chromatograph held at 180 C**. The volatile decomposition products generated were **trapped in a
  capillary column, cooled at -65 C**, and separated and identified by gas chromatography."
  The authors' own caution: "this technique varies according to the conditions used for
  decomposition of hydroperoxides, and **more variability in the data can be expected than in
  standard GC analyses**."
- **The one loading the paper prints.** "In a typical run, a **one microliter injection** was made
  from a **200 microliter hexane solution containing 13.1 mg methyl linoleate hydroperoxides, 2.4 mg
  a-tocopherol (15 wt % a-tocopherol), and 1.17 mg methyl hexanoate as internal standard**." So one
  injection carries about **65.5 µg of hydroperoxide and 5.85 µg of methyl hexanoate** (mine, 1/200
  of each). The "15 wt %" is on the **hydroperoxide-plus-tocopherol** basis: 2.4/(13.1 + 2.4) =
  15.5 %, against 2.4/13.1 = 18.3 % on hydroperoxide alone (mine). The paper states the basis
  nowhere else, and every other wt % in the two tables is unaccompanied by its own loading.
- **Chromatography.** Perkin-Elmer Sigma 300, capillary **DB-5, 60 m x 0.315 mm, 1 micron film**,
  cooled to -65 C with liquid nitrogen; initial hold 5 min, programmed to 260 C at 5 C/min, final
  hold 20 min. Identities "confirmed by EI mass spectrometry (10, 18)".
- **Quantification: relative peak area against an internal standard, and nothing absolute.** "Total
  volatiles were calculated as **percent of peak areas of volatiles relative to the peak area of the
  internal standard methyl hexanoate**. **The relative standard deviations of duplicate GC analyses
  ranged between 4-5 %.**" The table footnotes give the same figure as "**+/- 3.9 and +/- 4.8 %**".
  **There is no response factor for any product**, so a percent-of-six-peaks is a percent of *areas*,
  not of moles or of mass — for a slate that runs from a C5 alkane to a C14 dienoate, that is a
  material limitation and it is what `species_lipid.PENTANE`'s note already records.
- **The mechanism the repository takes from the introduction and Scheme 1** (attributed by the
  authors to their refs 3-10, all pre-1989, and **not** derived from the tocopherol arms):
  - **Homolytic beta-scission** of the alkoxyl radical, **pathway A**: from the 13-hydroperoxide,
    **pentane + methyl 13-oxo-9,11-tridecadienoate**; from the 9-hydroperoxide, **methyl octanoate +
    2,4-decadienal**. Pathway A is the favoured homolytic route because pathway B "is energetically
    less favorable because the heat of formation and the related bond dissociation energy required
    for the formation of a **vinyl radical** is larger".
  - **Homolytic pathway B**: from the 13-hydroperoxide, **hexanal**; from the 9-hydroperoxide,
    **methyl 9-oxononanoate**.
  - **Heterolytic (Hock) cleavage**, "between the carbon bearing the hydroperoxide group and the
    allylic double bond": from the 13-hydroperoxide, **hexanal + 12-oxo-10-dodecenoic acid**; from
    the 9-hydroperoxide, **2-nonenal + 9-oxononanoic acid**.
  - Scheme 1's substituents: **R = (CH2)7COOCH3, R' = (CH2)4CH3.**
  - So hexanal and methyl 9-oxononanoate are reachable by two routes and the other four by one.
    That is exactly `species_lipid.CLEAVAGE_MECHANISM` ("both" for those two, "homolytic" for the
    other four), and this dossier confirms the mapping against the source.
  - **The Hock partners 2-nonenal and methyl 12-oxo-10-dodecenoate are named in the introduction and
    measured in no table.** `species_lipid.NAMED_UNQUANTIFIED_COPRODUCTS` is the correct place for
    them.
- **The paper's own summary statistic (Figure 4).** "the sum of volatiles from the heterolytic
  and/or homolytic pathway B (**hexanal and methyl 9-oxononanoate**), was divided by the sum of
  volatiles from homolytic pathway A (**pentane, methyl octanoate**). This product ratio increased
  with increasing a-tocopherol concentration, to a greater extent with the trans,trans hydroperoxide
  isomers than with the mixture of cis,trans and trans,trans isomers (Fig. 4). **No increase was
  observed with the pure cis,trans isomers.**" Figure 4 is figure-only; the ratio can be recomputed
  from Tables 1 and 2, and is, in section 3.

## 3. Tables re-typed

Both tables use the same footnotes, printed identically: *a* "Individual volatiles are reported as
percent relative to the sum of six major volatiles observed. Relative standard deviations from
duplicate GC analyses ranged between +/-3.9 and +/-4.8 %." *b* "Relative to that of Me hexanoate
used as internal standard."

### Table 1. "Effect of a-Tocopherol and 1,4-Cyclohexadiene on the Volatile Thermal Decomposition Products of Mixed Isomers of Methyl Linoleate Hydroperoxides, cis,trans and trans,trans 9- and 13-Hydroperoxides"

Twelve columns: eight alpha-tocopherol levels and four 1,4-cyclohexadiene levels, both in wt %.
Column "0" is the same control for both blocks (it is printed once).

| Major volatiles | \_ | a-Toc 0 | 11 | 15 | 26 | 30 | 33 | 39 | 58 | 1,4-CHD 31 | 62 | 95 | 97 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pentane | | 16 | 10 | 10 | 9.7 | 7.8 | 8.8 | 6.7 | 7.1 | 9.0 | 9 | 9.5 | 12 |
| Hexanal | | 11 | 17 | 16 | 13 | 16 | 22 | 16 | 18 | 17 | 17 | 18 | 16 |
| Me Octanoate | | 17 | 13 | 12 | 12 | 11 | 12 | 10 | 10 | 13 | 13 | 13 | 3.2 |
| 2,4-Decadienal | | 23 | 23 | 22 | 24 | 23 | 23 | 24 | 20 | 21 | 20 | 20 | 25 |
| Me 9-oxononanoate | | 13 | 18 | 17 | 15 | 20 | 20 | 17 | 23 | 18 | 20 | 20 | 24 |
| Me 13-oxo-9,11-tridecadienoate | | 20 | 19 | 24 | 26 | 22 | 14 | 26 | 22 | 22 | 21 | 20 | 20 |
| **Total peak areas** (footnote b) | | **16** | **15** | **7.7** | **9.2** | **12** | **6.9** | **7.8** | **6.5** | **5.5** | **7.2** | **5.5** | **5.0** |

Column sums of the six shares (mine): 100, 100, 101, 99.7, 99.8, 99.8, 99.7, 100.1, 100, 100, 100.5,
100.2. Every column closes to within about one point, as a normalised distribution should.

### Table 2. "Effect of a-Tocopherol on the Volatile Thermal Decomposition Products of Different Isomers of Methyl Linoleate Hydroperoxides, Relative Percent"

Fourteen columns in two blocks: **cis,trans-13-hydroperoxide** (eight alpha-tocopherol levels) and
**trans,trans 9 + 13-hydroperoxides** (six levels), both in wt %.

| Major volatiles | ct-13: 0 | 9.4 | 11 | 21 | 35 | 38 | 39 | 48 | tt 9+13: 0 | 10 | 11 | 18 | 31 | 52 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pentane | 21 | 23 | 22 | 20 | 22 | 20 | 22 | 24 | 1.5 | 2.4 | 1.9 | 2.3 | 1.1 | 5.6 |
| Hexanal | 20 | 15 | 17 | 18 | 13 | 14 | 17 | 22 | 13 | 12 | 17 | 16 | 13 | 26 |
| Me octanoate | 5.0 | 5.3 | 3.6 | 4.1 | 3.4 | 6.0 | 3.1 | 5.4 | 13 | 13 | 15 | 9.8 | 8.3 | 4.6 |
| 2,4-Decadienal | 3.5 | 4.7 | 2.1 | 4.4 | 18 | 6.9 | 1.9 | 5.7 | 30 | 32 | 28 | 27 | 29 | 10 |
| Me 9-oxononanoate | 4.3 | 4.4 | 3.2 | 3.9 | 2.8 | 5.4 | 2.9 | 3.7 | 26 | 26 | 30 | 30 | 30 | 47 |
| Me 13-oxo-9,11-tridecadienoate | 46 | 48 | 52 | 50 | 40 | 48 | 53 | 39 | 16 | 15 | 8.7 | 15 | 18 | 6.3 |
| **Total peak areas** (footnote b) | **2.2** | **0.77** | **0.71** | **0.93** | **0.54** | **0.34** | **0.35** | **0.33** | **2.3** | **1.5** | **1.3** | **1.0** | **0.99** | **0.77** |

Column sums of the six shares (mine): 99.8, 100.4, 99.9, 100.4, 99.2, 100.3, 99.9, 99.8; 99.5,
100.4, 100.6, 100.1, 99.4, 99.5.

### The figure captions, transcribed as statistics (the plots themselves are figure-only)

| figure | system | printed statistics: probability of a zero slope (R2) for pentane, hexanal, methyl octanoate, 9-oxononanoate |
|---|---|---|
| 1 | mixed ct and tt 9- and 13-, vs wt % a-tocopherol | <0.01, 0.03, <0.01, <0.01 (R2: 0.64, 0.47, 0.79, 0.70) |
| 2 | mixed, vs wt % 1,4-cyclohexadiene | <0.16, 0.02, 0.05, <0.01 (R2: 0.36, 0.68, 0.57, 0.87) |
| 3 | trans,trans 9 and 13-, vs wt % a-tocopherol | 0.13, 0.05, <0.01, <0.01 (R2: 0.34, 0.49, 0.82, 0.72) |
| — (in the text, p. 604) | **pure cis,trans-13**, vs wt % a-tocopherol, "plots not shown" | **0.05, 0.01, 0.00, 0.02 (R2: 0.49, 0.75, 0.95, 0.65)**, described as "**plots not significantly different from zero**" |
| 4 | (hexanal + Me 9-oxononanoate) / (pentane + Me octanoate) vs wt % a-tocopherol, three systems overlaid | no statistics printed |

### Derived numbers (mine, arithmetic on the tables above)

**Figure 4's ratio, recomputed from the printed shares** — (hexanal + Me 9-oxononanoate) /
(pentane + Me octanoate), at zero additive and at the highest tocopherol of each system:

| system | at 0 wt % | at the highest a-tocopherol | change |
|---|---:|---:|---|
| mixed ct/tt 9+13 (Table 1) | (11 + 13)/(16 + 17) = **0.73** | 58 wt %: (18 + 23)/(7.1 + 10) = **2.40** | x3.3 |
| pure ct-13 (Table 2) | (20 + 4.3)/(21 + 5.0) = **0.93** | 48 wt %: (22 + 3.7)/(24 + 5.4) = **0.87** | x0.94 (no change) |
| tt 9+13 (Table 2) | (13 + 26)/(1.5 + 13) = **2.69** | 52 wt %: (26 + 47)/(5.6 + 4.6) = **7.16** | x2.7 |

This reproduces the Discussion's claim exactly: the ratio rises with tocopherol in the mixed and the
trans,trans systems and does not rise in the pure cis,trans one.

**Total volatiles, the paper's other hold-out claim** (total peak area relative to methyl hexanoate,
zero additive over highest additive): mixed **16 -> 6.5** with 58 % tocopherol (0.41x) and
**16 -> 5.0** with 97 % cyclohexadiene (0.31x); pure ct-13 **2.2 -> 0.33** (0.15x); tt 9+13
**2.3 -> 0.77** (0.33x). "This decrease in total volatiles was observed for all the hydroperoxide
samples tested."

**The 1 : 1 pairings the mechanism implies, and how far the data are from them.** Pathway A on the
13-hydroperoxide makes pentane and methyl 13-oxo-9,11-tridecadienoate as the two halves of one
scission; pathway A on the 9-hydroperoxide makes methyl octanoate and 2,4-decadienal as the two
halves of another. Across the three zero-additive columns:

| pairing | mixed | pure ct-13 | tt 9+13 | span |
|---|---:|---:|---:|---:|
| pentane / Me 13-oxo-tridecadienoate | 16/20 = 0.80 | 21/46 = 0.46 | 1.5/16 = 0.094 | **8.5x** |
| Me octanoate / 2,4-decadienal | 17/23 = 0.74 | 5.0/3.5 = 1.43 | 13/30 = 0.43 | **3.3x** |

The 8.5x is the number `kinetic_core_b6_prereg.md` pre-registered as expectation F-1, and this
dossier confirms it against the source. Either the GC areas of a C5 alkane and a C14 dienoate are
not commensurable after cryo-trapping at -65 C, or the pairing is not 1 : 1 in this experiment; the
paper's own Discussion supplies a reason for the second reading — "These two conjugated dienals,
which are expected products of homolytic pathway A in Scheme 1, **may be more reactive and less
stable than their saturated counterparts under the conditions used in this study**."

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`):** hexanal -> `hexanal`; **nonanal -> `nonanal`
(present in the registry, absent from this paper — that pairing is the point of the negative test)**;
2-pentylfuran -> `2_pentylfuran` (present in the registry, absent from this paper).
**Not in the registry:** pentane, methyl octanoate, 2,4-decadienal, methyl 9-oxononanoate, methyl
13-oxo-9,11-tridecadienoate, methyl hexanoate (the internal standard), 2-nonenal, methyl
12-oxo-10-dodecenoate, alpha-tocopherol, 1,4-cyclohexadiene, methyl linoleate and its hydroperoxides
— i.e. **five of the six products the lipid lane models have no compound id**, and the two named but
unquantified Hock partners have none either.

All rows share these conditions: **thermolysis in a GC injector port at 180 C**, hexane solution,
products cryo-trapped at -65 C on a DB-5 capillary, quantified as peak area relative to methyl
hexanoate, duplicate analyses with RSD 3.9-4.8 %, one temperature, no reaction time, no reactor
volume, no conversion.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **branch distribution, mixed ct/tt 9- and 13-hydroperoxides, no additive** | pentane 16 · hexanal 11 · Me octanoate 17 · 2,4-decadienal 23 · Me 9-oxononanoate 13 · Me 13-oxo-9,11-tridecadienoate 20 | % of the six-peak area sum | 180 C, autoxidised methyl linoleate hydroperoxide | Table 1 col. "0", p. 604 | **within_study_ratio** (the lane's FIT row) |
| **branch distribution, pure cis,trans-13-hydroperoxide, no additive** | 21 · 20 · 5.0 · 3.5 · 4.3 · 46 | % of sum | 180 C, soy-lipoxygenase preparation | Table 2 ct-13 col. "0", p. 606 | **within_study_ratio** (FIT row) |
| **branch distribution, trans,trans 9 + 13-hydroperoxides, no additive** | 1.5 · 13 · 13 · 30 · 26 · 16 | % of sum | 180 C, HPLC-separated tt fraction | Table 2 tt col. "0", p. 606 | **within_study_ratio** (FIT row) |
| total peak area at zero additive, three systems | 16 · 2.2 · 2.3 | relative to methyl hexanoate | as above; **the three preparations were injected at different loadings** | Tables 1, 2 footnote b | level_only — **not comparable across systems** |
| the 23 additive columns (8 + 4 tocopherol/cyclohexadiene in Table 1; 7 + 5 tocopherol in Table 2) | as re-typed in section 3 | % of sum, and relative total area | 180 C, additive wt % as printed | Tables 1, 2 | within_study_ratio — **HOLD-OUT under `FIT_HOLDOUT_DECLARATION.md` D.6 Module 5**; see flag 1 |
| hexanal's measured range | **11 to 20 %** across the three zero-additive columns; **11 to 26 %** across all 26 columns (mine, checked cell by cell) | % of the six-peak sum | 180 C | Tables 1, 2 | within_study_ratio |
| (hexanal + Me 9-oxononanoate) / (pentane + Me octanoate) | 0.73 / 0.93 / 2.69 at zero additive; rising x3.3 and x2.7 with tocopherol in the mixed and tt systems, flat (x0.94) in the pure ct-13 system | — | 180 C | derived by me from Tables 1-2; the paper plots it as Figure 4 | within_study_ratio (derived) |
| total-volatile suppression by a hydrogen donor | 0.41x (mixed, 58 % toc), 0.31x (mixed, 97 % CHD), 0.15x (ct-13, 48 % toc), 0.33x (tt, 52 % toc) | — | 180 C | derived by me from the footnote-b rows | within_study_ratio (derived, hold-out side) |
| pentane / Me 13-oxo-tridecadienoate, and Me octanoate / 2,4-decadienal | 0.80 / 0.46 / 0.094 and 0.74 / 1.43 / 0.43 | — | 180 C, zero additive | derived by me | within_study_ratio (derived) — the falsified 1 : 1 pairings |
| statistical significance of the tocopherol effect | p(zero slope) and R2 as tabulated in section 3 | — | four products x four systems | Fig. 1-3 captions and text p. 604 | figure_only (the plots); the statistics are printed |
| loading of a typical run | 13.1 mg hydroperoxide + 2.4 mg a-tocopherol + 1.17 mg methyl hexanoate in 200 µL hexane; 1 µL injected | mg per 200 µL; = 65.5 / 12.0 / 5.85 µg injected (mine) | 180 C injector | p. 604 | level_only |
| **rate constant at any temperature** | **NOT PRESENT** | — | — | — | — |
| **activation energy** | **NOT PRESENT** | — | — | — | — |
| **absolute yield, mass yield or mol %** | **NOT PRESENT** | — | — | — | — |
| **nonanal** | **NOT PRESENT — no table, no figure, no sentence** | — | — | — | the declared negative-test hold-out |
| **2-pentylfuran** | **NOT PRESENT** | — | — | — | — |
| Scheme 1, Figures 1-4 | mechanism and four regression plots | — | — | — | figure_only |

### The six products the repository's slate takes, and where each comes from

`species_lipid.LIPID_PRODUCTS` takes exactly the six rows of Tables 1 and 2, in the printed order,
and `parameters_lipid.FRANKEL_ZERO_ADDITIVE` takes exactly the three "0" columns. Verified cell by
cell against the source:

| repo key | Frankel's row label | mixed "0" | pure ct-13 "0" | tt 9+13 "0" | position (`POSITION_PRODUCTS`) | mechanism (`CLEAVAGE_MECHANISM`) |
|---|---|---:|---:|---:|---|---|
| `PENTANE` | Pentane | 16.0 | 21.0 | 1.5 | 13 | homolytic (pathway A) |
| `HEXANAL` | Hexanal | 11.0 | 20.0 | 13.0 | 13 | **both** (Hock and homolytic B) |
| `ME_OCTANOATE` | Me Octanoate | 17.0 | 5.0 | 13.0 | 9 | homolytic (pathway A) |
| `DECADIENAL` | 2,4-Decadienal | 23.0 | 3.5 | 30.0 | 9 | homolytic (pathway A) |
| `ME_9_OXONONANOATE` | Me 9-oxononanoate | 13.0 | 4.3 | 26.0 | 9 | **both** |
| `ME_13_OXO_TRIDECADIENOATE` | Me 13-oxo-9,11-tridecadienoate | 20.0 | 46.0 | 16.0 | 13 | homolytic (pathway A) |
| (total peak area, carried but not fitted) | Total peak areas | 16.0 | 2.2 | 2.3 | — | — |

**The branch fractions in the source are these percentages and nothing else.** They are percent of a
six-peak area sum within one column — footnote a says so — which means (i) they are shares, not
yields, and cannot be turned into a yield per hydroperoxide without a mass balance the paper does not
provide; (ii) they are area shares, not molar shares, because no response factor is given. That is
why `parameters_lipid.PROHIBITED_DERIVATIONS` refuses "a mass yield of hexanal per hydroperoxide" and
why `LIPID_FRAG_C` exists to absorb the unclosed remainder. Nothing in the full text contradicts
either decision.

### Do the tocopherol arms change the branch fractions?

**Yes in two of the three systems, no in the third, and the paper's whole argument is that
asymmetry.** The direction is uniform where it exists: **hexanal and methyl 9-oxononanoate rise,
pentane and methyl octanoate fall, and the total falls in every system.** In the **mixed** system
hexanal goes 11 -> 22 and methyl 9-oxononanoate 13 -> 23 while pentane goes 16 -> 6.7 and methyl
octanoate 17 -> 10 (Table 1, zero to the highest tocopherol levels); 1,4-cyclohexadiene, a pure
hydrogen donor with no antioxidant chemistry, does the same thing (hexanal 11 -> 16-18, methyl
9-oxononanoate 13 -> 18-24, methyl octanoate 17 -> 3.2 at 97 wt %), which is the paper's control for
the mechanism. In the **trans,trans** system hexanal goes 13 -> 26 and methyl 9-oxononanoate
26 -> 47 while methyl octanoate goes 13 -> 4.6 and 2,4-decadienal 30 -> 10. In the **pure
cis,trans-13** system the four tracked products are flat: "**Linear regression plots of these four
volatiles vs the concentration of a-tocopherol gave plots not significantly different from zero**".
Two products never trend anywhere: "With increasing a-tocopherol concentration, **no trends were
noted in the relative amounts of methyl 13-oxotridecadienoate and 2,4-decadienal**."

The authors' reading, which is the mechanism `species_lipid.CLEAVAGE_MECHANISM` encodes: hydrogen
donors "block the homolytic pathways A and B without affecting heterolytic pathways", so the
products reachable only homolytically lose share to the two that are also reachable by the Hock
route, and the total falls. The trans,trans isomers are the ones that do the shifting — "the
trans,trans isomers of 9- and 13-hydroperoxides of methyl linoleate were previously found to be more
susceptible to heterolytic decomposition by acid than the corresponding cis,trans isomers (22)".

**A caveat the tables carry and the prose does not.** In the trans,trans block, pentane is 1.5 at
zero additive and **5.6** at 52 wt % tocopherol — it rises, against the stated trend. Figure 3's
caption gives that regression a zero-slope probability of **0.13**, the only non-significant one in
the four systems, so the paper's claim for pentane in the trans,trans system rests on a
non-significant fit and the single highest-tocopherol column runs the other way. Any hold-out score
that treats "pentane falls" as a prediction in the tt system should note this.

### What the source does NOT contain, and what the repository does instead

1. **A rate at cooking temperature — or at any temperature.** One temperature (180 C), no reaction
   time, no reactor volume, no conversion, no concentration-time data, no Arrhenius plot, no
   activation energy. Reaction chromatography cannot supply one: the sample is decomposed inside an
   injector port during an injection. The authors close the paper by asking for exactly this
   experiment. **The lipid lane's absolute rate is therefore a declared assumption**, anchored at
   Schroen & Berton-Carabin 2022's `k4 = 6e-3 h^-1` at 25 C with a Q10 band of [2, 3] carried as a
   user-visible parameter (`parameters_lipid.K_LOOH_DECOMP_ANCHOR`, `Q10_ASSUMPTION`,
   `kinetic_core_b6_prereg.md` section 1), and every lipid prediction returns
   `in_envelope_extrapolated`. This dossier does not change that and confirms it cannot be closed
   from this paper.
2. **Nonanal.** The word does not appear in any table, figure, caption, scheme or sentence. Nonanal
   is the C9 fragment of the **oleate** double bond and this paper fed pure methyl **linoleate**
   hydroperoxides, so its absence is structural rather than a detection failure. That is the
   declared negative-test hold-out and the reason `species_lipid.NONANAL` has exactly one incoming
   edge, from `LOOH_OL`, whose branch fraction is `None`. The shipped FAST-lane value
   `nonanal 0.15` in `data/lit/lipid_oxidation_calibration.json` is refuted **structurally**, not by
   a number.
3. **2-Pentylfuran.** Not in the slate, not named anywhere in the paper. The alkylfuran route has no
   measured branch fraction here or elsewhere in the corpus, which is why
   `parameters_lipid.PROHIBITED_DERIVATIONS` refuses one and why 2-pentylfuran is not a species
   despite having a registry id.
4. **Propanal.** Also absent, and for a reason worth recording: propanal is a linolenate product and
   Frankel fed linoleate only. `PROHIBITED_DERIVATIONS` refuses importing one.
5. **Any absolute yield.** Shares and a relative total area, nothing else. The `hexanal 0.37` in the
   shipped calibration file sits **above** the paper's entire measured range (11-26 % over all 26
   columns), which is `SHIPPED_VALUES_REFUTED`'s claim, and this dossier confirms the range against
   the raster.
6. **The Hock partners.** 2-Nonenal and methyl 12-oxo-10-dodecenoate are named in the introduction
   and quantified nowhere, so no share exists for them.
7. **Any aqueous, emulsified or protein system.** The medium is hexane, injected into a hot metal
   port. There is no water, no pH, no water activity, no matrix and no amine — which also means this
   paper contributes nothing to the aldehyde-lysine channel the lipid lane carries as a bounded,
   inert ceiling.

## 5. Flags

1. **The FIT and HOLD-OUT columns are printed in the same table rows, and the abstract states the
   hold-out result in prose.** `kinetic_core_b6_prereg.md` section 0 discloses this and scores the
   tocopherol hold-out `seen_diagnostic`; reading the full text confirms there is no way to extract
   a zero-additive column without seeing its neighbours. **This dossier is itself a disclosure
   event**: it re-types all 26 columns, because the house rule requires every printed table to be
   re-typed. Any future wave using this file must treat the tocopherol and cyclohexadiene columns as
   seen, and the B6 firewall test (which asserts no hold-out number appears in
   `src/kinetic_core/*lipid*.py` or in the B6 fit report) still governs the code — this file is not
   code and is not covered by that grep, so the discipline has to be carried by the reader.
2. **Percent of six peaks is not percent of products.** The six were chosen as "the six major
   volatiles observed"; the paper says other volatiles exist ("Other radical reactions lead to the
   formation of relatively minor amounts of additional volatiles") and quantifies none of them, and
   two of the six scissions' partners (2-nonenal, methyl 12-oxo-10-dodecenoate) are named and never
   measured. **The slate does not close.** `LIPID_FRAG_C` is the correct handling and its size is
   the honest measure of the ignorance.
3. **No response factors.** A percent of *areas* on an FID for a slate running from C5H12 to a C14
   dienoate is not a percent of moles or of mass. The pentane / methyl 13-oxo-tridecadienoate ratio
   spanning 8.5x across three preparations of the same molecule (section 3) is the visible symptom.
   Do not convert any Frankel share to mmol/L.
4. **The three systems' total peak areas are not comparable** (16, 2.2, 2.3 at zero additive), and
   `parameters_lipid.FRANKEL_ZERO_ADDITIVE_TOTAL_AREA` correctly carries them without fitting them:
   the three preparations were injected at different loadings and the paper prints only one run's
   loading. Within one system the totals ARE comparable across additive levels, which is what makes
   the suppression claim meaningful.
5. **Every additive concentration is a wt % with no stated basis except one.** Only the "typical run"
   sentence ties a wt % to milligrams, and it does so at 15 %; the other 22 additive columns
   (11, 26, 30, 33, 39, 58, 31, 62, 95, 97, 9.4, 21, 35, 38, 39, 48, 10, 11, 18, 31, 52) have no
   loading. **A 97 wt % 1,4-cyclohexadiene column is a different chemical system**, not a perturbed
   version of the control, and the methyl octanoate collapse to 3.2 in that column should not be
   read as a dose response.
6. **One trend runs backwards.** trans,trans pentane rises 1.5 -> 5.6 with tocopherol while the
   paper claims it falls; Figure 3 gives that regression p = 0.13, the only non-significant slope
   among the sixteen. Section 4 states the consequence for hold-out scoring.
7. **Printed dimensions that cannot be right as printed.** The semi-preparative HPLC column is given
   as "25.0 X 2.14 mm" at 3.0 mL/min — a 25 mm x 2.14 mm column cannot carry preparative loads at
   that flow, and the numbers are the standard 25.0 cm x 2.14 cm semi-preparative format. Recorded as
   printed; it affects no number in this dossier.
8. **The RSD is quoted twice with two values.** Methods: "ranged between 4-5 %"; both table
   footnotes: "ranged between +/-3.9 and +/-4.8 %". Take the footnote range; it is on duplicate GC
   analyses only, and it does not include preparation-to-preparation variation, of which there is
   none measured (one preparation per system).
9. **The mechanism assignment is citation, not measurement.** The homolytic/heterolytic product map
   in the introduction and Scheme 1 is attributed to refs 3-10, all pre-1989, and none of them is on
   disk. `species_lipid.CLEAVAGE_MECHANISM`'s docstring already says the assignment comes from the
   introduction rather than from the tocopherol arms; that is the right provenance, and this dossier
   adds that the introduction's own provenance is second-hand.
10. **A second consumption route for hexanal is named and not modelled.** Introduction: "Hexanal may
    form either after rearrangement of the 9-hydroperoxide to the 13-hydroperoxide (13), or **after
    oxidative decomposition of 2,4-decadienal (14)**" (Schieberle & Grosch 1981).
    `parameters_lipid.DECADIENAL_TO_HEXANAL_GAP` already carries this as declared and not modelled;
    if that route is live, hexanal and 2,4-decadienal are not independent branches and the simplex
    the lane fits is mis-specified. Nothing in this paper tests it.
11. **Registry gaps against `data/keys/compounds.yml`:** five of the six modelled products
    (pentane, methyl octanoate, 2,4-decadienal, methyl 9-oxononanoate, methyl
    13-oxo-9,11-tridecadienoate) have no compound id, nor do methyl hexanoate, 2-nonenal, methyl
    12-oxo-10-dodecenoate, alpha-tocopherol, 1,4-cyclohexadiene, methyl linoleate or linoleic acid.
    The registry does hold `hexanal`, `nonanal` and `2_pentylfuran` — the last two being precisely
    the compounds this paper does **not** contain, which is a legible way for a coverage report to
    go wrong.
12. **What to request.** (i) The one experiment the authors themselves ask for: linoleate
    hydroperoxide decomposition products at temperatures below 180 C, by a method that gives a
    time course — that is the single retrieval that would replace the lane's Q10 assumption with a
    measurement. (ii) A response-factor-corrected or mass-balanced version of the same slate, which
    would let the shares become yields and let `LIPID_FRAG_C` shrink to a measured remainder.
    (iii) An oleate-hydroperoxide slate, which is the only thing that could ever give the nonanal
    branch fraction the engine currently refuses.

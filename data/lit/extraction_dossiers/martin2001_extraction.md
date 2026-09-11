# Martin & Ames 2001 — EXTRACTION (real potato slices stripped of sugars and amino acids by water/ethanol soaking, re-infused with chosen sugars and amino acids, deep-fried in palmolein 2 min from 180 C falling to 160 C; four Strecker aldehydes AND ten pyrazines measured by dynamic-headspace Tenax GC-MS in the SAME pots, six single-amino-acid pots plus a six-amino-acid competition pot)

### The one paper in the corpus that measures Strecker aldehydes and pyrazines side by side on a common glucose pool with one amino acid changed at a time — and its most valuable content is not the amino-acid identity ranking (which is a cross-compound peak-area comparison with no response factors, and whose four single-amino-acid pyrazine totals are all one statistical group) but the COMPETITION result: put six amino acids on one glucose pool and leucine's aldehyde goes UP 1.24x while isoleucine's, phenylalanine's and methionine's fall 2.4x, 2.7x and 5.2x, and total pyrazines fall to 36.5 % of the sum of the singles.

**Source on disk:** `data/articles/martin2001.pdf` (8 pp., J. Agric. Food Chem. 2001, 49 (8), 3885-3892).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/martin2001.txt`, 578 lines). **All six tables came through and are re-typed in
full below.** Tables 3, 4 and 5 are wide sparse grids in which blank cells collapse to whitespace,
so I recovered every cell by matching character offsets against the header row rather than by eye;
the assignment is independently verified against four statements the paper makes about its own
numbers (section 3, "Column assignment verified"). The only true damage is one wrapped cell:
Table 5's `Z` column for the parent pyrazine (4.5b) is pushed onto its own line. Figures 1 and 2
are images. **Figure 2 — "Total relative yield of pyrazines in model systems" — is the ONLY place
the paper's headline yield comparison appears, and it is figure-only; but it can be reconstructed
arithmetically from Table 4's printed totals and Table 2's printed uptakes, and I do that below.**
There is no supplementary material; the fuller account is F. L. Martin's Reading PhD thesis
(reference 17), not on disk. Repo status before this dossier: **Martin & Ames 2001 has no
extraction dossier and is not cited by any file in `src/kinetic_core/`.**

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Strecker Aldehydes and Pyrazines in a Fried Potato Model System" |
| Authors | Fiona L. Martin, Jennifer M. Ames (corresponding) — School of Food Biosciences, The University of Reading, P.O. Box 226, Whiteknights, Reading RG6 6AP, United Kingdom |
| Venue | J. Agric. Food Chem. 2001, 49 (8), 3885-3892. Received 7 March 2001, revised and accepted 25 May 2001, web 11 July 2001 |
| DOI / article ID | 10.1021/jf010310g (printed as `JF010310G`) |
| Funding | BBSRC (U.K.) and United Biscuits (U.K.), via a studentship to F.L.M. Potatoes, palmolein and the tuber composition data supplied by United Biscuits |
| Naming | **"relative amount" (RA) = GC peak area divided by the internal standard's**, two significant figures. **"relative yield" = relative amount divided by the moles of amino acid infused.** Model-system codes are single letters (Table 2): `W` blank (steeped in water only), `U` untreated (never stripped, never infused), `A` `Gn` `L` `I` `P` `M` single amino acid, `G` glucose alone, `GA` `GGn` `GL` `GI` `GP` `GM` glucose + one amino acid, `X` glucose + all six, `F` fructose alone, `FA` fructose + asparagine, `GFA` both sugars + asparagine, `Y` both sugars + all six amino acids, `Z` = Y + threonine |
| Lineage | the stripped-slice model is developed from Khanbari & Thompson 1993 (ref 14), a fry-colour method; the Strecker-aldehyde reactivity comparison is against Hofmann, Münch & Schieberle 2000 (ref 32); pyrazine yield comparisons against Koehler 1969 (ref 39) and Chun & Ho 1997 (ref 40) |
| Companions on disk | none — this is the only fried-potato model-system paper in the corpus |

## 1. Why it matters

Two pre-registrations stalled this week on the same missing object.

`results/validation/kinetic_core_b19_prereg_draft.md` section 5 refuses a per-amino-acid Strecker
rate fit and names what it would need instead: "an identity layer on the one Strecker step already
fitted: glycine's two constants from the pyrazine step (`FROZEN_B18`) anchor the rate and its
barrier; each other amino acid enters as a **partition ratio of the same dicarbonyl pool**, FIT on
the within-study ratios". It lists the ratio sources it has — Balagiannis' yield fractions
(Ile : Leu 1.6), Kocadagli 2021's aldehyde ratios, Amrani-Hemaimi 1995's isotope fractions — and
says the pre-registration "is a separate document to write when the three ratio sources have been
re-read for the ratio rows".

`results/validation/kinetic_core_b22_prereg.md` section 6 records the wave that did run and was
refused: methional could not be reproduced as free dicarbonyl times methionine at any single
identity ratio, because "no single ratio serves the two pots". Its stated next structure is a
methionine Amadori compound with its own first-order decomposition.

**Martin & Ames 2001 is a fourth ratio source, and it is the only one that measures the Strecker
aldehydes and the pyrazines in the same pots.** That matters specifically because the B19 layer
proposes to anchor every amino acid's Strecker rate to *glycine's constants from the pyrazine
step*: the two channels are being tied together, and this paper is the only place in the corpus
where both channels are read off the same dicarbonyl pool with one amino acid swapped at a time.
Six pots (`GA`, `GGn`, `GL`, `GI`, `GP`, `GM`) share an identical glucose loading (258.7 mg per
100 g of slices) and an identical thermal history, and each carries one amino acid at a printed
molar loading; a seventh (`X`) puts all six on that same pool at once.

What it delivers, in order of usefulness:

1. **A competition experiment, fully printed.** `X` against the singles is a direct test of the
   "partition of one pool" structure B19 wants to write. It fails additivity in an informative
   direction (section 3, arithmetic 4): leucine's aldehyde is **superadditive** (1.24x its single-pot
   value) while the other three aldehydes fall 2.4-5.2x and total pyrazines fall to 36.5 %. **A
   partition layer fitted on single-amino-acid pots will not reproduce a mixed pot**, and this paper
   says so with numbers rather than in principle.
2. **A methionine result that bears directly on the refused wave.** In `GM`, dimethyl disulfide
   (330) is **3.0x higher than methional** (110), and in the methionine-only pot `M` the factor is
   **16.7x**. Any methionine step that routes the amino acid to methional as its principal fate is
   wrong on carbon by at least a factor of four in a frying matrix. The paper's own reading is that
   methionine is oxidised to the sulfoxide first (ref 38) and that methional is a minority branch.
3. **A pyrazine identity ratio across six amino acids at constant glucose** — with two heavy
   qualifications that the dossier's section 4 states plainly: the four small pots `GL`, `GI`, `GP`,
   `GM` are all **one statistical group (superscript "a")** on the total row, and the glucose-only
   blank `G` supplies 27-50 % of their totals.
4. **A within-study Ile : Leu aldehyde ratio on isomers**, which is the one cross-compound
   comparison in this paper that survives the missing response factors, and is directly comparable
   with the Balagiannis 1.6 the B19 draft already holds.

What this paper does NOT give: **any rate constant, any activation energy, any time course, any
concentration in mass or molar units, any pH, any water activity.** One frying time (2 min), one
thermal profile, one matrix. Every number in it is a headspace peak-area ratio. It cannot supply a
`k`, and nothing in it belongs in a `MEASURED_*` registry as a rate.

## 2. Methods as they matter to a model

- **The matrix is real potato, stripped and re-loaded.** *Solanum tuberosum* cv. Saturna (a chipping
  cultivar), specific gravity 1.086-1.096, supplied by United Biscuits. Washed, peeled, sliced to
  **1.4 ± 0.1 mm**, rinsed 2 min in cold water, held in water at room temperature for at most 5 h.
- **Stripping.** 100 g of slices into 500 mL water at **60 C, stirred, 5 min**; transferred to
  500 mL **50 % ethanol at 45 C, stirred, 15 min**; rinsed with **four 1 L aliquots** of cold water;
  stored in water **overnight at 4 C**. Glucose removal checked with Clinistix reagent strips,
  detection limit ~0.5 mM ~ 0.1 mg/mL — glucose was **below** it. Fructose and the amino acids were
  *assumed* removed on a solubility argument, not measured (Flags 2). "The potato cell wall network
  and starch were not removed."
- **Re-infusion.** 100 g of slices into **200 mL of steeping solution at 40 C for 10 min**, then
  drained, blotted and **immediately fried**. Uptake determined by capillary electrophoresis on the
  steeping solution before and after use, by difference, with the printed formula
  `uptake per 100 g = (C V1) - (area2 x C V2 / area1)`. **Triplicate; average standard deviation
  < 15 %.** "The presence of more than one component in the steeping solution did not significantly
  affect uptake by the slices" (Figure 1, for glucose).
- **Frying.** 3 L of palmolein preheated to **180 C for 20 min** in a Magimix Sélection professional
  deep-fat fryer. **25 g of slices** added, and *the fryer was switched off at that moment and stayed
  off for the whole 2 min* — a deliberate choice to make the thermal history reproducible. The
  printed profile: **180 C falling to ~165 C over the first 30 s, then a further 5 C to 160 C
  between 30 and 150 s.** Fresh palmolein each day, maximum 30 batches per day. Drained and cooled
  on paper towels, sealed in metalized laminate chip bags.
- **Volatile isolation is a dynamic headspace purge from a slurry, not an extraction of the chip.**
  **10 g of crushed chips + 40 mL of water** in a 250 mL flask at **37 C**; **nitrogen at
  40 mL/min for 1 h** through a Tenax TA trap (155 mm x 3 mm i.d., 85 mg). The internal standard,
  **1,2-dichlorobenzene, 1 uL of a 130.6 ug/mL methanol solution = 0.131 ug**, is injected **onto
  the trap**, not into the sample. **So the internal standard corrects GC-MS response and
  desorption, and corrects nothing about purge efficiency, matrix retention or partition.** Every
  "relative amount" in this paper is a *release-weighted* number, and the weighting differs between
  a C5 aldehyde, an aromatic aldehyde, methional and a pyrazine (Flags 1).
- **GC-MS.** HP 5890 II + HP 5972 MSD; CP-SIL8, 60 m x 0.25 mm i.d., 0.25 um film; helium
  1.5 mL/min. First 0.5 m cooled in solid CO2 for 4 min with the oven off; trap desorbed at 280 C
  into that cryofocus for 5 min; oven then 40 C (2 min), 4 C/min to 200 C, 10 C/min to 250 C,
  15 min hold. EI 70 eV, source 165-175 C, scan m/z 32-450, 1.82 scans/s. Identification by mass
  spectrum **and** linear retention index against C6-C22 n-alkanes; RI values printed in Tables 3
  and 4 against literature values, and they agree to within 4 units on every compound that carries
  both.
- **Replication and statistics.** All systems in triplicate. **CV < 25 %** on the relative amounts.
  One-way ANOVA plus multiple-range tests (Statgraphics Plus 4.1); **means carrying different
  superscript letters within a row differ at P < 0.05**. The letters are the paper's only error
  information — there is no standard deviation printed anywhere in Tables 3-6.
- **What is not measured at all.** No pH. No water activity or moisture, before or after frying. No
  sugar or amino acid measurement *in the fried chip* (only uptake into the raw slice). No
  dicarbonyl measurement — glyoxal, methylglyoxal and 3-deoxyglucosone are discussed as mechanism
  and never quantified. No time series: **one frying time.** No oil-uptake measurement, though the
  paper attributes benzaldehyde partly to linoleic acid in the palmolein.

## 3. Tables re-typed

### Table 1 (p. 3885). "Concentrations of Sugars and Amino Acids in a Potato Cultivar Used for Chipping (Saturna) (United Biscuits, Personal Communication)"

**This is a literature/communication value, not measured in this study.**

| sugar | concn (g/100 g) |
|---|---|
| glucose | 0.1 |
| fructose | 0.08 |
| sucrose | 1.07 |

| amino acid | concn (mg/100 g) | amino acid | concn (mg/100 g) |
|---|---|---|---|
| Ala | 4.7 | Lys | 4.7 |
| Arg | 16.4 | Met | 4.7 |
| Asn | 93.9 | Phe | 4.7 |
| Asp | 4.7 | Pro | 4.7 |
| Gln | 28.2 | Ser | 4.7 |
| Glu | 9.4 | Thr | 18.8 |
| Gly | 0 | Trp | 0 |
| His | 7 | Tyr | 7 |
| Ile | 7 | Val | 9.4 |
| Leu | 4.7 | | |

### Table 2 (p. 3886). "Composition and Codes of the Model Systems"

Footnote a: "Initial concentration of components in the steeping solution." Footnote b:
"Concentration of components in potato slices. Average of three experiments. Average standard
deviation < 15 %." Footnote c: "Blank. Steeped in water only."

| component | initial concn^a (mg/100 mL) | uptake^b (mg/100 g) | W^c | A | Gn | L | I | P | M | G | GA | GGn | GL | GI | GP | GM | X | F | FA | GFA | Y | Z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| glucose | 500 | **258.7** | | | | | | | | x | x | x | x | x | x | x | x | | | x | x | x |
| fructose | 500 | **218.2** | | | | | | | | | | | | | | | | x | x | x | x | x |
| asparagine | 475 | **210** | | x | | | | | | | x | | | | | | x | | x | x | x | x |
| glutamine | 150 | **82.3** | | | x | | | | | | | x | | | | | x | | | | x | x |
| leucine | 25 | **11.6** | | | | x | | | | | | | x | | | | x | | | | x | x |
| isoleucine | 35 | **15** | | | | | x | | | | | | | x | | | x | | | | x | x |
| phenylalanine | 25 | **11.9** | | | | | | x | | | | | | | x | | x | | | | x | x |
| methionine | 25 | **13.7** | | | | | | | x | | | | | | | x | x | | | | x | x |
| threonine | 100 | **46.9** | | | | | | | | | | | | | | | | | | | | x |

### Table 3 (p. 3888). "Relative Amounts^a,^b of Strecker Aldehydes and Their Selected Degradation Products Formed in Model Potato Chips"

Footnote a: "Amounts of components are quoted in relative GC peak area units to two significant
figures... Figures quoted are the means of triplicate analyses. **CV < 25 %.**" Footnote b: "Means
with different superscript letters within a row are significantly different (P < 0.05)." Footnote c:
"Calculated RI values for identified components." Footnote d: "Linear retention indices obtained for
authentic compounds analyzed on the same GC column or from the literature (19, 20)." Footnote e:
"See Table 2 for explanation of codes other than U. Model system U was prepared by frying potato
slices without removal of or infusion with sugars and amino acids."

Empty cells are empty in the printed table (compound not reported for that system).

| compound | RI exptl^c | RI lit.^d | W | U | L | I | P | M | G | GL | GI | GP | GM | X |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-methylbutanal | 671 | 655 | 2a | 940c | 770b | 3a | 4a | 4a | 30a | **4500d** | — | 30a | 33a | **5600e** |
| 2-methylbutanal | 680 | 665 | 2a | 840c | 9a | 550b | 3a | 2a | 15a | 55a | **4600e** | 31a | 18a | **1900d** |
| phenylacetaldehyde | 1063 | 1066 | — | 89b | — | — | 69b | — | 2a | 2a | 4a | **490d** | 2a | **180c** |
| benzaldehyde | 981 | 983 | 13ab | 26bc | 6a | 4a | 14ab | 5a | 7a | 8a | 16abc | **85d** | 7a | 30c |
| methional | — | 924 | — | 6a | — | — | — | 6a | — | — | — | — | **110c** | 21b |
| dimethyl sulfide | — | — | — | 6 | — | — | — | 3 | — | — | — | — | 4 | 1 |
| dimethyl disulfide | 753 | 744 | 6a | 290cd | 6a | 5a | 7a | 100b | 19a | 8a | 7a | 17a | **330d** | 240c |
| dimethyl trisulfide | 989 | 990 | — | 8 | — | — | — | 7 | 1 | — | — | — | 7 | 6 |

### Table 4 (p. 3888). "Relative Amounts^a,^b of Selected Pyrazines Formed in Model Potato Chips"

Same four footnotes as Table 3. **There is no W column in this table.**

| pyrazine | RI exptl^c | RI lit.^d | U | A | Gn | G | GA | GGn | GL | GI | GP | GM | X |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pyrazine | — | — | 6ab | — | 1a | 1a | 57d | 38c | 5ab | 6ab | 11b | 6ab | 38c |
| methylpyrazine | 836 | 833 | 97b | 3a | 3a | 6a | 130c | 120bc | 13a | 10a | 26a | 17a | 91b |
| 2,5(6)-dimethylpyrazine | 924 | 925 | 80d | — | — | 6a | 32c | 22abc | 7a | 6a | 11ab | 13ab | 25bc |
| ethylpyrazine | 927 | 930 | 45b | — | — | 1a | 91c | 43b | 4a | 5a | 8a | 4a | 82c |
| 2,3-dimethylpyrazine | 930 | 932 | 15d | — | — | 1a | 19e | 6b | 1a | 1a | 1a | 1a | 10c |
| vinylpyrazine | 946 | 948 | 9b | — | — | — | 62e | 22c | 2ab | 2ab | 3ab | 2ab | 33d |
| 2-ethyl-6-methylpyrazine | 1009 | 1010 | 29b | — | — | — | 8a | 6a | — | — | — | — | 9a |
| 2-ethyl-3(5)-methylpyrazine | 1014 | 1016 | 27d | — | — | 1a | 15c | 5ab | — | — | — | 3a | 12bc |
| 2-vinyl-6-methylpyrazine | 1033 | 1034 | 11 | — | — | — | 12 | 10 | — | — | — | — | 14 |
| 3-ethyl-2,5-dimethylpyrazine | 1086 | 1086 | 40b | — | — | — | 3a | — | — | — | — | — | 3a |
| **total** | | | **359bc** | **3a** | **4a** | **16a** | **429c** | **272b** | **32a** | **30a** | **60a** | **46a** | **317b** |

### Table 5 (p. 3890). "Percentage RA Values for Pyrazines^a,^b"

Footnote a: "Percentage of the sum of the relative amount (RA) for all monitored pyrazines. Figures
quoted are the means of triplicate analyses. CV < 25 %." Footnotes b and c as before.

| component | U | G | GA | GGn | GL | GI | GP | GM | X | F | FA | GFA | Y | Z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pyrazine | 1.8a | 8.4c | 13de | 14ef | 16f | 12d | 18g | 13de | 12d | 1.4a | 0.6a | 5.4b | 4.8b | 4.5b |
| methylpyrazine | 27a | 34bcd | 30ab | 44gh | 40ef | 37de | 43fg | 36d | 29a | 37de | 32bc | 47h | 35cd | 33bcd |
| 2,5(6)-dimethylpyrazine | 22cd | 37g | 7.3a | 7.5a | 23d | 23de | 18bc | 29f | 8.0a | 55i | 45h | 17b | 26def | 28ef |
| ethylpyrazine | 13de | 8.3b | 21g | 16f | 12de | 18fg | 13e | 9.4bc | 26h | 3.6a | 4.7a | 10bcd | 12de | 12cde |
| 2,3-dimethylpyrazine | 4.1de | 2.8bcd | 4.5e | 1.8ab | 1.9ab | 2.2abcd | 2.3bcd | 3.1bcde | 3.2bcde | 1.8ab | 1.2a | 2.6bcd | 3.8cde | 2.1abc |
| vinylpyrazine | 2.4a | 1.5a | 14f | 7.7d | 6.6cd | 7.7cd | 5.9bc | 4.7b | 11e | — | 1.2a | 7.0cd | 4.6b | 4.8b |
| 2-ethyl-6-methylpyrazine | 8.0f | — | 1.9a | 2.0ab | — | — | — | — | 2.9d | — | 2.6c | 2.2b | 3.1d | 4.4e |
| 2-ethyl-3(5)-methylpyrazine | 7.3ef | 7.6f | 3.6abcde | 1.8abc | — | — | — | 5.5cdef | 3.6abcd | — | 7.3ef | 6.0def | 4.6bcdef | 6.5def |
| 2-vinyl-6-methylpyrazine | 3.2cd | — | 2.7bc | 3.7de | — | — | — | — | 4.3f | — | 2.0a | 2.3ab | 4.2ef | 3.8e |
| 3-ethyl-2,5-dimethylpyrazine | 11f | — | 0.4a | — | — | — | — | — | 0.9ab | 2.3d | 3.5e | 0.8ab | 1.1b | 1.8c |

### Table 6 (p. 3890). "Percentage RA Values for Strecker Aldehydes^a,^b"

Footnotes as for Table 5.

| compound | U | X | Y |
|---|---|---|---|
| 3-methylbutanal | 50 | 73 | 67 |
| 2-methylbutanal | 45 | 25 | 30 |
| phenylacetaldehyde | 4.7b | 2.3a | 2.8a |
| methional | 0.3 | 0.3 | 0.5 |

### Column assignment verified

Tables 3-5 are sparse and the blanks carry no placeholder, so I checked the recovered grid against
four independent statements the paper makes about its own numbers. All four reproduce:

1. **Table 4's `total` row is the exact sum of its column** for every one of the eleven systems
   (U 359, A 3, Gn 4, G 16, GA 429, GGn 272, GL 32, GI 30, GP 60, GM 46, X 317 — checked
   individually, mine). A mis-assigned cell would break at least one of these.
2. **Table 5's columns sum to ~100 %** (GI 99.9, GL 99.5, GP 100.2, GM 100.7 — mine).
3. **Table 4 divided by its own totals reproduces Table 5** to within two-significant-figure
   rounding on 40 of the 41 cells I checked. The exception is Flags 6.
4. **The text's own ratio claims come out right.** "This ratio [methylpyrazine : 2,5(6)-dimethyl-
   pyrazine] is < 3 for GL, GI, GP, and GM but > 4 for GA and GGn": I get 1.74, 1.61, 2.39, 1.24 and
   4.11, 5.87 (mine). "3-methylbutanal : 2-methylbutanal is 1.1, 3.0 and 2.2, respectively, in U, X
   and Y": Table 6 gives 1.11, 2.92, 2.23 (mine). "3-ethyl-2,5-dimethylpyrazine ... in U (11.1 %)"
   against the printed 11f, and "~10-fold lower in Y" against Y = 1.1 (mine).

### Numbers printed in the running text (everything else is in the tables or figure-only)

| quantity | value | where |
|---|---|---|
| glucose in stripped slices | below the Clinistix limit, **~0.5 mM ~ 0.1 mg/mL** | Results, p. 3886 |
| model loading vs the real tuber | "concentrations that were proportional to levels in a potato cultivar used for chipping ..., although the absolute concentrations were **~2-fold higher**" | Results, p. 3886 |
| thermal profile of the fry | **180 -> ~165 C over the first 30 s, then a further 5 C to 160 C between 30 and 150 s** | Results, p. 3887 |
| residual amino acids in the blank W | 2- and 3-methylbutanal in W were "at least **275-385-fold lower** than in the models infused with these amino acids" | Results, p. 3887 |
| **glucose effect on Strecker aldehydes** | "The relative amount of each Strecker aldehyde in the systems infused only with amino acid was **~12 %** of that in the corresponding system containing glucose; that is, the relative rate of reaction with and without glucose was **independent of the amino acid**" | Discussion, p. 3888 (my check: 5.5-17.1 %, Flags 4) |
| **3-methylbutanal in the competition pot** | "The relative amount of 3-methylbutanal in X was **~120 %** of that in GL" | Discussion, p. 3888 |
| the other three aldehydes in X | "relative amounts of 2-methylbutanal, phenylacetaldehyde, and methional were **60-80 % lower** than those in GI, GP, and GM" | Discussion, p. 3888 (my check: 59 %, 63 %, **81 %**) |
| **dimethyl disulfide vs methional** | "relative amounts being higher than those of methional by factors of **17, 3, and 11**, respectively, in M, GM, and X" | Discussion, p. 3888 (my check: 16.7, 3.0, 11.4) |
| **total pyrazines in the competition pot** | "the total relative amount of pyrazines in X being **38 %** of the sum of the values in the single amino acid-glucose systems" | Discussion, p. 3889 (my check: 36.5 %) |
| which pyrazine escapes the competition | "With the exception of **3-ethyl-2,5-dimethylpyrazine**, the relative amount of each individual pyrazine in X was lower than the sum of its relative amounts in the single amino acid-glucose systems" | Discussion, p. 3889 (my check: 3 vs 3, exactly 1.00) |
| **which amino acid gives the most pyrazine per mole** | "**GP** gave a total relative yield of pyrazines that was significantly higher than those for the other models containing glucose as the only sugar" | Discussion, p. 3888 — **read off Figure 2** |
| glutamine vs asparagine | "The total relative yield of pyrazines in **GGn was higher than that in GA, but not significantly so**" | Discussion, p. 3888 — **Figure 2** |
| the sugar was limiting | "The significantly higher total relative yield for pyrazines in **GFA compared to both GA and FA** (Figure 2) confirms that sugar was limiting pyrazine formation in both"; and "The total relative yield of pyrazines in **Y was significantly higher than in X** (Figure 2), confirming that the level of reducing sugar was limiting pyrazine formation in the latter" | Discussion, pp. 3889-3890 — **Figure 2** |
| 3-ethyl-2,5-dimethylpyrazine and threonine | percentage relative amount rose from **1.1 % in Y to 1.8 % in Z**, a significant increase, but still far from U's **11.1 %** | Discussion, p. 3890 |
| fructose vs glucose, the mechanism claim | 2,5(6)-dimethylpyrazine, 2-ethyl-6-methylpyrazine and 3-ethyl-2,5-dimethylpyrazine significantly **higher** in F and FA than in G and GA; pyrazine, ethylpyrazine, vinylpyrazine and 2-vinyl-6-methylpyrazine significantly **lower** — read as **more pyruvaldehyde from fructose than from glucose** | Discussion, p. 3890 |

**Figure-only in this paper (not typed as numbers).** **Figure 1**, glucose uptake with and without
amino acids — superseded by Table 2's printed uptakes, and the paper states "No significant
difference in uptake of glucose was observed among the systems." **Figure 2, "Total relative yield
of pyrazines in model systems" with significance letters — the only place the relative-yield
comparison the abstract and conclusion rest on is drawn, and the only place the fructose systems F,
FA, GFA, Y and Z have any absolute quantity at all.** Everything in Table 5 for those five systems
is a percentage composition.

### Arithmetic on the printed numbers (all mine)

**1. The molar loadings, which the paper never prints in moles.** From Table 2's uptakes
(mg per 100 g of slices), with standard molecular weights:

| component | uptake (mg/100 g) | MW (g/mol) | **mmol/kg of slices** |
|---|---|---|---|
| glucose | 258.7 | 180.16 | **14.36** |
| fructose | 218.2 | 180.16 | **12.11** |
| asparagine | 210 | 132.12 (anhydrous) | **15.90** (13.99 if the monohydrate, Flags 5) |
| glutamine | 82.3 | 146.15 | **5.63** |
| leucine | 11.6 | 131.17 | **0.884** |
| isoleucine | 15 | 131.17 | **1.144** |
| phenylalanine | 11.9 | 165.19 | **0.720** |
| methionine | 13.7 | 149.21 | **0.918** |
| threonine | 46.9 | 119.12 | **3.937** |

**The pots are not equimolar and the imbalance is not uniform.** Glucose : amino acid is **0.90 in
GA and 2.55 in GGn, but 16.2, 12.6, 19.9 and 15.6 in GL, GI, GP and GM.** So comparing GA and GGn
against the four small pots as an *identity* ratio confounds identity with stoichiometry by more
than an order of magnitude in sugar excess; comparing GL, GI, GP and GM against **each other** does
not (their glucose excesses agree to within 1.6x). This is the single most important structural fact
for anyone taking ratios out of this paper.

**2. Figure 2 reconstructed from printed numbers.** Total relative yield = Table 4's total divided
by the amino acid's mmol/kg:

| system | total RA (Table 4) | amino acid (mmol/kg) | **relative yield per mmol/kg** | normalised to Leu = 1.00 |
|---|---|---|---|---|
| GA (asparagine) | 429 | 15.90 | **27.0** | 0.75 |
| GGn (glutamine) | 272 | 5.63 | **48.3** | 1.33 |
| GL (leucine) | 32 | 0.884 | **36.2** | 1.00 |
| GI (isoleucine) | 30 | 1.144 | **26.2** | 0.72 |
| GP (phenylalanine) | 60 | 0.720 | **83.3** | **2.30** |
| GM (methionine) | 46 | 0.918 | **50.1** | 1.38 |
| X (all six) | 317 | 25.20 (sum) | **12.6** | 0.35 |

This reproduces both claims the paper reads off Figure 2 — **GP is the highest** and **GGn is above
GA (48.3 vs 27.0, 1.79x)** — from printed numbers only. It is a derivation, not a figure read, and
it is the most useful single object this paper offers the identity layer.

**3. The same, background-corrected, and it moves by up to 2x.** The glucose-only pot `G` produces
16 total pyrazine units with no infused amino acid — i.e. **27-50 % of the totals of GL, GI, GP and
GM.** Subtracting it: GL 16, GI 14, GP 44, GM 30, GA 413, GGn 256; yields 18.1, 12.2, 61.1, 32.7,
26.0, 45.5; normalised to Leu = 1.00: **Ile 0.68, Phe 3.38, Met 1.81, Asn 1.43, Gln 2.51.** Against
the uncorrected 0.72, 2.30, 1.38, 0.75, 1.33. **Asparagine moves by 1.9x and phenylalanine by
1.5x.** The paper does not subtract the blank and does not discuss it. Any ratio taken from this
paper carries that factor-of-two ambiguity, and the dossier records it rather than choosing.

**4. The competition experiment, compound by compound.** X against the corresponding single pot:

| compound | single pot | X | **X / single** |
|---|---|---|---|
| 3-methylbutanal | GL 4500 | 5600 | **1.24** (superadditive) |
| 2-methylbutanal | GI 4600 | 1900 | 0.41 |
| phenylacetaldehyde | GP 490 | 180 | 0.37 |
| methional | GM 110 | 21 | **0.19** |
| dimethyl disulfide | GM 330 | 240 | 0.73 |
| **total pyrazines** | sum of GA+GGn+GL+GI+GP+GM = 869 | 317 | **0.365** |
| parent pyrazine | sum 123 | 38 | 0.31 |
| methylpyrazine | sum 316 | 91 | 0.29 |
| 2,5(6)-dimethylpyrazine | sum 91 | 25 | 0.27 |
| ethylpyrazine | sum 155 | 82 | 0.53 |
| 2,3-dimethylpyrazine | sum 29 | 10 | 0.34 |
| vinylpyrazine | sum 93 | 33 | 0.35 |
| 2-ethyl-6-methylpyrazine | sum 14 | 9 | 0.64 |
| 2-ethyl-3(5)-methylpyrazine | sum 23 | 12 | 0.52 |
| 2-vinyl-6-methylpyrazine | sum 22 | 14 | 0.64 |
| 3-ethyl-2,5-dimethylpyrazine | sum 3 | 3 | **1.00** |

**A single shared pool, partitioned linearly, predicts 1.00 in every row of this table.** The
aldehyde rows span 0.19 to 1.24 and the pyrazine rows 0.27 to 1.00. That is the quantitative content
of the B19 draft's "partition ratio of the same dicarbonyl pool" being tested and coming out wrong
by a factor of five across amino acids — with the direction informative: **the fastest reactant
(leucine) gains at the expense of the rest**, which is what competition for a limiting pool does when
the rate constants differ, and what a linear partition cannot represent.

**5. The Ile : Leu aldehyde ratio, the one cross-compound ratio that survives.** 2-methylbutanal and
3-methylbutanal are structural isomers (C5H10O, MW 86.13) with near-identical EI response, Tenax
retention and volatility, so their peak-area ratio is a defensible mole ratio in a way that
methional : 3-methylbutanal is not. Relative amounts: **GI 4600 : GL 4500 = 1.02**. Per mole of
amino acid infused: 4600/1.144 = 4021 against 4500/0.884 = 5090, i.e. **Ile : Leu = 0.79**, or
**Leu : Ile = 1.27**. In the competition pot X the ratio inverts: 1900 : 5600 = **0.34** at equal
molar loading of the two amino acids. The B19 draft holds Balagiannis' **Ile : Leu = 1.6**; this
paper's single-pot value is **0.79** and its mixed-pot value **0.34**. Three numbers spanning 4.7x
for the same nominal ratio (Flags 8).

**6. The glucose effect, and it is not "independent of the amino acid".** Amino-acid-only pot
divided by the same amino acid with glucose: 3-methylbutanal L/GL = 770/4500 = **17.1 %**;
2-methylbutanal I/GI = 550/4600 = **12.0 %**; phenylacetaldehyde P/GP = 69/490 = **14.1 %**;
methional M/GM = 6/110 = **5.5 %**. Mean 12.2 %, which is the paper's "~12 %", but the spread is
**3.1x**, and the claim in the same sentence that "the relative rate of reaction with and without
glucose was independent of the amino acid" is not supported by these four numbers. Methionine is the
outlier in both directions — least helped by glucose in methional, and its dimethyl-disulfide branch
is helped only 3.3x (M 100 -> GM 330).

**7. The model loading against the real tuber (mine).** Uptake divided by Table 1: Asn 2.24x,
Gln 2.92x, Leu 2.47x, Ile 2.14x, Phe 2.53x, Met 2.91x, Thr 2.49x, glucose 2.59x, fructose 2.73x.
The paper's "~2-fold higher" is really **2.1-2.9x**, a 1.4x spread — so the *ratios* between
precursors are preserved to within 40 %, which is what matters if this pot is used as a proxy for a
real chip.

**8. Total aldehyde is nearly conserved while total pyrazine collapses (mine, and weakly held).**
Summing the four Strecker aldehydes: singles 4500 + 4600 + 490 + 110 = 9700 against X's
5600 + 1900 + 180 + 21 = 7701, i.e. **79 %**, against **36.5 %** for total pyrazines. This is a sum
over different compounds and so carries the unmeasured response factors; take it as a direction, not
a number. The direction is chemically sensible: pyrazine needs two amino-acid-derived nitrogen
fragments to condense while a Strecker aldehyde needs one, so a nitrogen-fragment shortage should
hurt pyrazines quadratically and aldehydes linearly.

## 4. Kinetic numbers the repository can use

**There are none.** This paper contains no rate constant, no activation energy, no reaction order
and no time course. It has one heating time and one thermal profile. What it supplies is
**identity ratios**, and section 4's job is to say which of them the paper actually supports.

**The governing distinction.** Every number is a GC peak area divided by an internal standard
spiked onto the trap after a 1 h nitrogen purge from a 37 C slurry. Therefore:

- **Same compound, different pot -> the response factor and the purge efficiency cancel exactly.**
  These ratios are `within_study_ratio` and are usable.
- **Different compounds, same pot -> nothing cancels.** No response factors were determined and no
  authentic-standard calibration was run for quantification (authentic standards were used only for
  retention indices). These are `peak_area_only` and must never be read as mole ratios. **This
  disqualifies the paper's headline amino-acid ranking of Strecker aldehydes**, and it is why the
  abstract's "leucine gave the highest relative amount ... of its Strecker aldehyde" cannot be
  checked — and in fact conflicts with Table 3, where 2-methylbutanal in GI (4600) exceeds
  3-methylbutanal in GL (4500) (Flags 7).
- **The one exception** is 2- vs 3-methylbutanal, isomers whose response factors are close enough to
  make the comparison defensible with the caveat stated.

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Keyed and measured here:
`3_methylbutanal`, `2_methylbutanal`, `phenylacetaldehyde`, `benzaldehyde`, `methional`,
`dimethyl_disulfide`, `dimethyl_trisulfide`, `methylpyrazine`, `2_ethylpyrazine`,
`2_3_dimethylpyrazine`, and the aggregate `pyrazines`. Gaps in Flags 11.

All rows below share: cv. Saturna slices 1.4 mm, stripped by water/ethanol soaking, re-infused at
40 C for 10 min, deep-fried in palmolein **2 min from 180 C falling to 160 C**, volatiles purged
from a 10 g-in-40 mL slurry at 37 C with N2 for 1 h onto Tenax, triplicate, CV < 25 %, one time
point.

| quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **pyrazine identity ratio across six amino acids at one glucose loading**, per mole of amino acid, normalised to leucine | Asn **0.75**, Gln **1.33**, Leu **1.00**, Ile **0.72**, Phe **2.30**, Met **1.38** | dimensionless | glucose 14.36 mmol/kg in every pot; amino acid as in the table above | none — a single-time-point yield ratio, not a rate | derived from Table 4 totals + Table 2 uptakes (mine) | **within_study_ratio**, but see the three caveats immediately below |
| the same, with the glucose-only blank subtracted | Asn 1.43, Gln 2.51, Leu 1.00, Ile 0.68, Phe **3.38**, Met 1.81 | dimensionless | as above, minus G = 16 | — | derived (mine) | **within_study_ratio** — the alternative reading; the paper chooses neither |
| **per-compound pyrazine identity ratio** (methylpyrazine, the best-populated row), per mole, normalised to leucine | Asn 0.56, Gln 1.45, Leu 1.00, Ile 0.59, Phe **2.46**, Met 1.26 | dimensionless | as above | — | derived from Table 4 row 2 (mine) | **within_study_ratio** — the cleanest kind here: one compound, six pots, response factor cancels exactly |
| the same for the parent pyrazine | Asn 0.63, Gln 1.19, Leu 1.00, Ile 0.93, Phe **2.71**, Met 1.16 | dimensionless | as above | — | derived from Table 4 row 1 (mine) | within_study_ratio |
| the same for 2,5(6)-dimethylpyrazine | Asn 0.25, Gln 0.49, Leu 1.00, Ile 0.66, Phe **1.93**, Met **1.79** | dimensionless | as above | — | derived (mine) | within_study_ratio (**a co-elution of the 2,5- and 2,6- isomers**, Flags 11) |
| **glucose multiplier on a Strecker aldehyde** (with glucose / without) | Leu **5.8x**, Ile **8.4x**, Phe **7.1x**, Met **18.3x** | dimensionless | GL/L, GI/I, GP/P, GM/M | — | Table 3 (mine) | **within_study_ratio** — same compound, two pots |
| **competition penalty**: six amino acids on one glucose pool, vs the same amino acid alone | 3-methylbutanal **1.24**, 2-methylbutanal **0.41**, phenylacetaldehyde **0.37**, methional **0.19**, dimethyl disulfide 0.73 | dimensionless | X vs GL/GI/GP/GM, identical glucose | — | Table 3 (mine, and the paper prints "~120 %" and "60-80 % lower") | **within_study_ratio** — the highest-value rows in this dossier |
| **competition penalty, pyrazines** | total **0.365**; per compound 0.27 to 1.00 | dimensionless | X vs the sum of the six single pots | — | Table 4 (mine); the paper prints 38 % | within_study_ratio |
| **Ile : Leu aldehyde ratio**, per mole of amino acid | **0.79** (single pots); **0.34** (competition pot X) | dimensionless | GI vs GL; within X | — | Table 3 (mine) | within_study_ratio — **the one cross-compound ratio defensible here, because the two aldehydes are isomers** |
| 3-methylbutanal : 2-methylbutanal percentage ratio | **1.1 (U), 3.0 (X), 2.2 (Y)** | dimensionless | untreated chip; glucose + six amino acids; both sugars + six amino acids | — | Table 6, printed in the text | within_study_ratio (isomers) |
| **dimethyl disulfide : methional** | **16.7 (M), 3.0 (GM), 11.4 (X)** | dimensionless | methionine alone; + glucose; + glucose and five other amino acids | — | Table 3; the paper prints 17, 3, 11 | **peak_area_only** for the absolute ratio (two different compounds), but the **change across pots (5.6x from M to GM) is a within_study_ratio** and is the load-bearing part |
| pyrazine percentage composition, all fourteen systems | 10 compounds x 14 systems | % of total pyrazine RA | as above | — | Table 5 | **peak_area_only** (a composition across different compounds) — usable as a fingerprint, not as mole fractions |
| Strecker aldehyde percentage composition | 4 compounds x 3 systems | % of total aldehyde RA | U, X, Y | — | Table 6 | peak_area_only |
| all absolute levels | 4500, 429, 110 ... | relative GC peak area units | as above | — | Tables 3, 4 | **peak_area_only** — there is no mass or molar concentration of any volatile anywhere in this paper |
| precursor loadings in the fried slice | glucose 14.36, fructose 12.11, Asn 15.90, Gln 5.63, Leu 0.884, Ile 1.144, Phe 0.720, Met 0.918, Thr 3.937 | mmol/kg of slices | uptake by CE difference, triplicate, SD < 15 % | — | Table 2 (mg/100 g printed; the molar conversion is mine) | **measured level** (a real analytical measurement, the only one in the paper) |
| thermal profile | 180 C -> ~165 C at 30 s -> 160 C at 150 s, 2 min total | C vs s | 25 g of slices into 3 L of palmolein, fryer switched off | — | Results p. 3887 | **level_only** (printed as prose, from thesis ref 17) |
| tuber composition, cv. Saturna | Table 1 | g/100 g and mg/100 g | — | — | Table 1 | **level_only** — a personal communication, not measured here |
| total relative yield of pyrazines for the fructose systems F, FA, GFA, Y, Z | — | — | — | — | Figure 2 | **figure_only** |
| glucose uptake with and without amino acids | — | — | — | — | Figure 1 | **figure_only** (and superseded by Table 2) |

### Which ratios does the paper support, and which does it not?

**Supported, and usable as ratio rows.**

1. **The competition rows** (X against the singles). Same compound, two pots, identical glucose
   loading and thermal history, and the paper prints the headline versions itself. These are the
   rows a partition layer should actually be tested against, because they are the only measurement in
   the corpus of what happens when several amino acids compete for one dicarbonyl pool in a real
   matrix. They say a linear partition is wrong by up to 5x, and in which direction.
2. **The glucose-multiplier rows** (with sugar over without). Same compound, two pots. They bound how
   much of each Strecker aldehyde is sugar-derived in this matrix: 82-95 %.
3. **The per-compound pyrazine identity ratios across GL, GI, GP and GM.** Same compound, four pots
   at the same glucose excess (12.6-19.9x, within 1.6x of each other). This is the paper's cleanest
   identity comparison **and it is not statistically resolved** — every one of GL, GI, GP and GM
   carries superscript "a" on Table 4's total row, i.e. **the four are not significantly different
   from one another at P < 0.05**. Read them as a two-level partition (Asn and Gln high; Leu, Ile,
   Phe, Met low and mutually indistinguishable), not as four numbers.
4. **The Ile : Leu aldehyde ratio**, with the isomer argument stated.
5. **The dimethyl-disulfide branch on methionine**, as a change across pots. This is the row that
   speaks to the refused methionine wave.

**Not supported, and figure-only or response-factor-blocked.**

6. **The amino-acid ranking of Strecker aldehydes** — "leucine gave the highest relative amount and
   relative yield of its Strecker aldehyde", "phenylalanine gave the highest total relative yield of
   pyrazines" for the *aldehyde* half. Comparing 3-methylbutanal against phenylacetaldehyde and
   methional requires response factors that were never measured. The paper's own ranking is also
   internally contradicted for the "relative amount" half (Flags 7).
7. **Every absolute yield in the fructose systems** (F, FA, GFA, Y, Z). Table 4 has no fructose
   columns; those systems appear only as percentage compositions in Table 5 and as bars in
   Figure 2. **The fructose-vs-glucose comparison the abstract and conclusion rest on is therefore a
   composition comparison, not a yield comparison**, except for the significance statements read off
   Figure 2.
8. **"Phenylalanine gave the highest total relative yield of pyrazines"** as printed is a Figure 2
   read — but section 3 arithmetic 2 reconstructs it from Table 4 and Table 2, so the *conclusion*
   is available from printed numbers even though the paper's own presentation of it is not.
9. **Anything about sugar identity in the mixed pots.** X versus Y differ in both fructose presence
   and total sugar loading (14.36 vs 26.47 mmol/kg), and the paper itself says "the level of reducing
   sugar was limiting pyrazine formation" in X. **The Y-vs-X comparison confounds sugar identity with
   sugar amount and cannot be used for either alone.** Only GA vs FA holds the amount roughly
   constant (14.36 vs 12.11 mmol/kg, a 1.19x difference), and even that pair is only available as
   percentages.

### What this changes for the two stalled pre-registrations

For `kinetic_core_b19_prereg_draft.md`: this paper adds a fourth ratio source with a genuinely
different matrix (real fried potato rather than an aqueous or dry model), and it adds the object the
draft's structure most needs and does not have — **a measurement of what a shared pool does when
several amino acids draw on it at once.** It also supplies a second, independent Ile : Leu value
(0.79 against Balagiannis' 1.6), which is a useful disagreement to carry rather than average.
What it cannot supply is any barrier: there is one temperature profile, so no amino acid's Ea can be
separated from any other's, which is consistent with the draft's decision to declare every barrier
equal to glycine's.

For `kinetic_core_b22_prereg.md`: the methional rows here point the same way the wave's own verdict
did. Methional is a **minority** methionine product in a frying matrix — three to seventeen times
below dimethyl disulfide in the same chromatogram — and it is the aldehyde most suppressed by
competition (5.2x, the largest penalty of the four). A structure that treats methional as
methionine's principal sink will mis-state the branch even before it gets the rate wrong. Nothing
here supports or refutes the proposed methionine-Amadori route directly, because no Amadori compound
was fed or measured.

## 5. Flags

1. **Every number in this paper is a dynamic-headspace peak area, and the internal standard does not
   correct for release.** 1,2-dichlorobenzene is spiked **onto the Tenax trap**, downstream of the
   purge, so it corrects desorption and GC-MS response and nothing about how much of each compound
   left the 37 C slurry in an hour of nitrogen. Methional, phenylacetaldehyde and 3-methylbutanal
   have very different volatilities and matrix affinities. **Cross-compound comparisons are blocked;
   same-compound cross-pot comparisons are fine**, because the matrix is nominally identical between
   pots.
2. **The stripping was verified only for glucose.** Fructose and all the amino acids were "assumed"
   removed on a solubility argument. The blank W then produces 2- and 3-methylbutanal, which the
   paper reads, correctly, as proof that leucine and isoleucine were not fully removed. **The
   glucose-only pot G produces 16 total pyrazine units and 30 units of 3-methylbutanal with no
   infused amino acid at all** — a background that is 27-50 % of the totals in GL, GI, GP and GM,
   and which the paper never subtracts. Section 3 arithmetic 3 shows the identity ratios move by up
   to 1.9x depending on whether it is subtracted.
3. **The four small single-amino-acid pyrazine totals are one statistical group.** GL 32a, GI 30a,
   GP 60a, GM 46a all carry "a" on Table 4's total row: **not significantly different at P < 0.05.**
   Only GA (429c) and GGn (272b) separate. The per-mole yields do separate GP (Figure 2, and my
   reconstruction), but the raw totals do not. **Do not fit four distinct partition ratios to Leu,
   Ile, Phe and Met from this paper.**
4. **"Independent of the amino acid" is over-stated.** The glucose effect ranges from 5.5 % (Met) to
   17.1 % (Leu), a 3.1x spread, around the quoted ~12 % mean (section 3, arithmetic 6).
5. **The asparagine molar loading is ambiguous by 14 %.** Asparagine is commonly supplied as the
   monohydrate (MW 150.13); the paper says only "Asparagine ... minimum purity = 99 % ... from
   Sigma". At 132.12 the uptake is 15.90 mmol/kg, at 150.13 it is 13.99. **Every asparagine ratio in
   this dossier carries that 1.14x ambiguity**, which is smaller than the blank-subtraction ambiguity
   but should be stated. It also shifts the glucose : asparagine molar ratio in GA between 0.90 and
   1.03.
6. **Tables 4 and 5 disagree for one cell.** Table 4's GI column gives the parent pyrazine as 6 out
   of a total of 30, i.e. 20 %; Table 5 prints 12 % for the same cell. The other 40 cells I checked
   reconcile to within two-significant-figure rounding, and Table 5's GI column sums to 99.9 %, so
   the likelier error is Table 4's "6ab" (12 % of 30 would be 3.6). A second, smaller disagreement is
   3-ethyl-2,5-dimethylpyrazine in GA: Table 4's 3 out of 429 is 0.70 %, against Table 5's 0.4 %.
7. **The abstract's aldehyde ranking is contradicted by Table 3.** "For the single amino acid-glucose
   systems, leucine gave the highest relative amount and relative yield of its Strecker aldehyde" —
   but 2-methylbutanal in GI is 4600 against 3-methylbutanal in GL at 4500. The **relative yield**
   half of the claim does survive (5090 against 4021 per mmol/kg, mine, because less leucine was
   infused); the **relative amount** half does not. Quote only the yield version.
8. **Three values of Ile : Leu now exist in the corpus and they span 4.7x**: Balagiannis 1.6 (held by
   the B19 draft), this paper's single-pot 0.79, this paper's competition-pot 0.34. They are not
   measuring the same thing — the competition value is a different experiment — but the 1.6 against
   0.79 is a genuine two-laboratory disagreement on the same nominal ratio and should be carried as
   a spread, not averaged.
9. **No pH, no water activity, no moisture, and no measurement of anything inside the fried chip.**
   Uptake is measured into the raw slice by difference on the steeping solution; what survives the
   fry is never measured. A fried 1.4 mm slice loses most of its water and takes up oil, and neither
   is quantified. **Nothing in this paper can be put on a water-activity or pH axis.**
10. **Two small internal inconsistencies in the methods.** The frying time is given as 2 min but the
    thermal profile runs to 150 s. And "the mixture containing ... 25 g of slices" is fried in 3 L of
    oil while volatiles are isolated from 10 g of *crushed chips*, with no statement of the yield of
    chips per gram of slices — so no mass basis connects a relative amount to a precursor loading.
    That, more than the response factors, is why nothing here can become an absolute yield.
11. **Registry gaps against `data/keys/compounds.yml`.** Present and measured here:
    `3_methylbutanal`, `2_methylbutanal`, `phenylacetaldehyde`, `benzaldehyde`, `methional`,
    `dimethyl_disulfide`, `dimethyl_trisulfide`, `methylpyrazine`, `2_ethylpyrazine`,
    `2_3_dimethylpyrazine`. **Absent: the parent (unsubstituted) pyrazine as a molecule** — the id
    `pyrazines` is a `compound_class` whose aliases include the bare word "pyrazine", and its own
    `identity_note` says that bare word "means the family", so **the paper's largest single pyrazine
    signal in GA (57 units) has no molecular id and would silently resolve to a class**; **vinyl-
    pyrazine**; **2-ethyl-6-methylpyrazine**; **2-vinyl-6-methylpyrazine**; **2-ethyl-3(5)-methyl-
    pyrazine**; **dimethyl sulfide**. Two further identity problems: (i) the paper reports
    **"2,5(6)-dimethylpyrazine"**, an unresolved co-elution, while the registry carries
    `2_5_dimethylpyrazine` and `2_6_dimethylpyrazine` as separate ids — the paper's number can be
    mapped to **neither**; (ii) the paper's **3-ethyl-2,5-dimethylpyrazine** (RI 1086, named in the
    introduction as an important chip aroma compound) is **not** the registry's
    `2_ethyl_3_5_dimethylpyrazine` (SMILES `CCc1ncc(C)nc1C`, CAS 13360-64-0) — they are different
    ring isomers, and the paper's compound has no id. Similarly "2-ethyl-3(5)-methylpyrazine" is
    printed with an explicit ambiguity between two isomers. **Six of the ten pyrazines in Table 4
    cannot currently be keyed.**
12. **What this paper does NOT contain**: any rate constant, activation energy, reaction order or
    time course; any concentration in mass or molar units for any volatile; any dicarbonyl
    measurement; any pH or water activity; any response-factor calibration; any standard deviation
    (only significance letters); any absolute yield for the fructose systems; any measurement of the
    precursors remaining after frying; any replicate-level data.
13. **What to request from the authors**: (i) response factors, or authentic-standard calibration
    curves, for 3-methylbutanal, 2-methylbutanal, phenylacetaldehyde, methional and the ten
    pyrazines — this single item would convert the whole paper from `peak_area_only` to a real
    identity-ratio source; (ii) the numeric values behind Figure 2, in particular for F, FA, GFA, Y
    and Z, which have no printed totals at all; (iii) the resolution of the Table 4 / Table 5
    disagreement on GI's parent pyrazine; (iv) the chip mass recovered per gram of slices, and the
    moisture and oil content of the fried chips, so a relative amount can be put on a mass basis;
    (v) whether the asparagine was the anhydrous form or the monohydrate; (vi) F. L. Martin's Reading
    thesis (2001), reference 17, which is cited for the thermal profile, the CE linearity and "more
    details concerning the development of the model system" and is the only route to the underlying
    numbers.

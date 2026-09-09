# Hwang, Hartman & Ho 1995b — EXTRACTION (the same nine pots as `hwang1995_extraction.md`: glucose + 15N-glycine + one tested amino acid, 20 g wheat starch, 12-14 % moisture, pH 7, 180 C / 1 h; 25 pyridines, pyrroles, oxazoles, amines and benzonitrile)

### The second half of one experiment, not a second experiment: the same nine reaction beds analysed for the non-pyrazine nitrogen heterocycles, with the same figure-borne glycine-versus-competitor split — but it prints three isotope-derived nitrogen partitions in prose, which the pyrazine paper does not.

**Source on disk:** `data/articles/hwang1995b.pdf` (5 pp., owner's download, 2026-09-08). J. Agric. Food
Chem. 1995, 43, 2917-2921. The text layer is an OCR layer that **dropped two whole columns of Table 1**
(Ctrl and Gln), so Table 1 was re-read cell by cell from a 200-dpi raster and a 320-dpi crop of printed
page 2919 (`scratchpad/img/h95b-3.png`, `h95bzoom-3.png`). **Seven of the nine column sums reproduce the
printed totals exactly; the Lys and Asp columns do not** (flag 2). Figures 1-4 (bar charts with
percentage labels) are FIGURE-ONLY; their labels are recorded in section 3 under an explicit warning and
are not typed as numbers in section 4.

## 0. Identity

| field | value |
|---|---|
| Title | "Relative Reactivities of Amino Acids in the Formation of Pyridines, Pyrroles, and Oxazoles" |
| Authors | Hui-Ing Hwang, Thomas G. Hartman, Chi-Tang Ho (Rutgers, New Brunswick NJ) |
| Venue | J. Agric. Food Chem. 1995, 43, 2917-2921; received 28 February 1995, accepted 30 June 1995 |
| Article id | JF950123G (no DOI printed) |
| Relationship to the companion | The authors state it plainly: "In an earlier paper (Hwang et al., 1995), we reported a total of 56 pyrazines in the reaction systems containing 15N labeled glycine and 8 other amino acids. **In this paper, we further report another 25 nitrogen-containing reaction products**". Same materials, same charges, same drying, same 180 C / 1 h, same purge-and-trap, same internal standard, same 15N arithmetic. **Same pots.** |
| Naming | Ctrl = glycine only; Gln / Lys / Asn / Phe / Glu / Asp / Ile / Arg = labelled glycine PLUS that amino acid. "Ref" in the figures = "Ctrl" in the table |
| Companion on disk | `hwang1995_extraction.md` |

## 1. Why it matters

The same refusals as its companion: `results/validation/kinetic_core_b19_prereg_draft.md` section 5 and
`results/validation/kinetic_core_b22_prereg.md` section 6 leave an identity-ratio layer on glycine's
fitted Strecker step (`FROZEN_B18` in `src/kinetic_core/parameters_pyrazine.py`) as the only honest
structure, fitted on within-study ratios of one amino acid against another in the same pot.

This paper adds three things its companion does not:

1. **Three isotope-derived nitrogen partitions printed in the running text**, not in a figure: >90 % of
   the nitrogen in 1-methylpyrrole and 1-methylpyrrole-2-carboxaldehyde comes from glycine; >80 % of
   the nitrogen in 4,5-dimethyloxazole comes from glycine; ~96 % of the nitrogen in benzonitrile comes
   from phenylalanine. Those are printed within_study_ratios of exactly the kind the layer wants — but
   they are per-compound, on compounds the repository does not carry.
2. A **second, independent ordering** of the same nine pots, on different products, which either
   corroborates the pyrazine ordering or does not. It does not: aspartic acid and isoleucine lead here,
   lysine and arginine led and trailed there (section 4).
3. Direct evidence that the nitrogen partition is **product-dependent within a single pot** — glycine
   dominates the N-methyl compounds, the competitor dominates its own side-chain-derived products — so
   a single scalar identity ratio per amino acid is an approximation the source itself contradicts.

It contains no rate, no barrier and no time course.

## 2. Methods as they matter to a model

Identical to `hwang1995_extraction.md` section 2, and re-verified against this paper's own text:

- **Charge (verbatim):** "Twenty grams of wheat starch and an equal amount (**2.66 µmol of each**) of
  glucose, L-glycine-α-amine-15N, and tested amino acid ... were mixed with **150 mL of deionized water
  and adjusted to pH 7**." Ctrl carries glycine only, hence half the total amine. The same
  µmol/mmol inconsistency applies (flag 1 of the companion, repeated as flag 4 here).
- **Drying / rehydration:** freeze-dried, then rehydrated over water in a desiccator to **12-14 %
  moisture** (the moisture figure is quoted here without the AOAC citation the companion gives).
- **Heating:** "transferred into a reaction vessel and heated at **180 °C for 1 h**". One point.
- **Isolation:** 2 g of heated sample, **1 µL of 1.001 mg/mL deuterated toluene** as internal standard,
  SIS solid-sample purge-and-trap, **nitrogen 40 mL/min, 80 °C, 1 h**, Tenax TA + Carbotrap.
- **Quantification:** GC-MS after Hwang et al. (1993); linear retention indices against C5-C25
  n-paraffins; NIST library or literature spectra. **No response factors, no calibration, no LOD, no
  replicate count, no error bars.**
- **Unit conflict:** Table 1's head reads "yield (**mg**/g of glucose)" while the y-axes of Figures 1-3
  read "Yield (**µg**/g glucose)". The figures are right and the table head is a typo: the pyrrole sums
  computed from Table 1 (Asp 226.0, Lys 198.4, Ctrl 58.4) reproduce the bar heights of Figure 2 on a
  0-250 µg/g axis, and the pyridine sums (Asp 63.5, Gln 49.0, Ctrl 27.5) reproduce Figure 1 on a 0-75
  axis. **Read Table 1 as µg per g of glucose**, the same unit as the companion (flag 3).
- **The 15N bookkeeping.** These products carry **one** ring nitrogen, not two, so the arithmetic is
  simpler than the pyrazines': each compound has W1 (14N, from the tested amino acid) and W2 (15N, from
  labelled glycine), solved from the M and M+1 abundances of the labelled and unlabelled runs, and

      % contribution of tested amino acid = [W1/(W1 + W2)] x 100 %
      % contribution of labelled glycine  = [W2/(W1 + W2)] x 100 %

  There is no half-weighted middle term, so this paper's percentages are a **cleaner partition** than
  the companion's — but they still assign every non-labelled nitrogen, α-amino and side chain alike, to
  the tested amino acid. No supplementary material is offered here (the companion's 13-page supplement
  covers the pyrazines only).

## 3. Tables re-typed

### Table 1. "Pyridines, Pyrroles, Oxazoles, and Other Nitrogen-Containing Compounds Identified in the Reaction of Glucose, Glycine-α-amine-15N, and Tested Amino Acids"

Unit as printed in the table head: "mg/g of glucose". **Read as µg/g of glucose** (section 2). "–" = not
observed. Footnote: "Ctrl, glycine only; Gln, labeled glycine and glutamine; Lys, labeled glycine and
lysine; Asn, labeled glycine and asparagine; Phe, labeled glycine and phenylalanine; Glu, labeled
glycine and glutamic acid; Asp, labeled glycine and aspartic acid; Ile, labeled glycine and isoleucine;
Arg, labeled glycine and arginine."

| compound | Ctrl | Gln | Lys | Asn | Phe | Glu | Asp | Ile | Arg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **pyridines** | | | | | | | | | |
| pyridine | – | – | – | 6.3 | – | – | 8.5 | – | – |
| 4-methylpyridine | – | – | – | – | – | – | 14.5 | – | – |
| 2-methylpyridine | – | – | – | – | – | – | 4.6 | – | – |
| 2-ethylpyridine | – | – | – | – | – | – | 8.3 | 2.2 | – |
| 2,5-dimethylpyridine | – | – | – | – | – | – | 1.4 | – | – |
| 3-ethyl-2,6-dimethylpyridine | 2.5 | 11.6 | – | 5.1 | 4.5 | 6.5 | – | – | – |
| 3-butylpyridine | 1.7 | 4.3 | – | – | – | 9.6 | 1.7 | – | – |
| 2-acetylpyridine | 23.3 | 31.5 | 33.3 | 14.3 | 10.2 | 32.2 | 18.0 | 5.5 | 6.0 |
| 2-propionylpyridine | – | 1.6 | – | 0.8 | – | – | 6.5 | – | – |
| 2,6-diphenylpyridine | – | – | – | – | 11.8 | – | – | – | – |
| **pyrroles** | | | | | | | | | |
| pyrrole | – | – | – | 9.3 | – | 6.1 | 10.0 | 5.0 | 4.8 |
| 1-methylpyrrole | 9.4 | 30.2 | 4.6 | 26.7 | 58.9 | 59.0 | 131.5 | 26.8 | 11.8 |
| 2,5-dimethylpyrrole | 1.2 | – | – | 4.4 | – | 3.9 | 4.0 | 4.7 | 1.8 |
| tetramethylpyrrole | – | – | – | – | – | 2.2 | – | – | – |
| 2-acetylpyrrole | 3.8 | 17.0 | 42.1 | 9.2 | – | 15.8 | 8.7 | – | 16.8 |
| 1-methyl-2-acetylpyrrole | 44.0 | 29.5 | 102.7 | 27.0 | 69.1 | 54.1 | 52.0 | 25.0 | 27.8 |
| 1-methylpyrrole-2-carboxaldehyde | – | 6.4 | – | – | 10.0 | – | – | 7.5 | 6.0 |
| 1-ethylpyrrole-2-carboxaldehyde | – | 20.1 | 49.0 | 6.1 | 14.2 | 14.5 | 9.6 | 18.6 | 8.3 |
| 1-(2-furanylmethyl)pyrrole | – | 8.2 | – | 6.9 | 6.1 | 9.2 | 10.2 | – | 6.4 |
| **oxazoles** | | | | | | | | | |
| 4,5-dimethyloxazole | – | 9.8 | 17.0 | – | – | 13.9 | 14.0 | – | 7.7 |
| trimethyloxazole | – | – | – | – | – | – | 9.0 | – | – |
| 2-acetyl-4,5-dimethyloxazole | – | – | – | – | – | – | 10.1 | – | – |
| **other nitrogen-containing compounds** | | | | | | | | | |
| benzonitrile | – | – | – | – | 2.4 | – | – | – | – |
| N-(2-methylbutylidene)-2-methylbutylamine | – | – | – | – | – | – | – | 69.7 | – |
| bis(2-methylbutyl)amine | – | – | – | – | – | – | – | 65.9 | – |
| tris(2-methylbutyl)amine | – | – | – | – | – | – | – | 83.0 | – |
| **totals (as printed)** | **85.9** | **170.2** | **248.57** | **116.1** | **187.2** | **227.0** | **322.4** | **313.9** | **97.4** |

**Arithmetic check (mine).** Ctrl 85.9 ✓, Gln 170.2 ✓, Asn 116.1 ✓, Phe 187.2 ✓, Glu 227.0 ✓, Ile 313.9 ✓,
Arg 97.4 ✓ — seven columns reproduce exactly. **Lys sums to 248.7 against a printed 248.57** (and it is
the only total printed to two decimals), and **Asp sums to 322.6 against a printed 322.4.** See flag 2.

**Class subtotals (mine), the quantities the figures plot:**

| pot | pyridines | pyrroles | oxazoles | other | total |
|---|---:|---:|---:|---:|---:|
| Ctrl | 27.5 | 58.4 | 0 | 0 | 85.9 |
| Gln | 49.0 | 111.4 | 9.8 | 0 | 170.2 |
| Lys | 33.3 | 198.4 | 17.0 | 0 | 248.7 |
| Asn | 26.5 | 89.6 | 0 | 0 | 116.1 |
| Phe | 26.5 | 158.3 | 0 | 2.4 | 187.2 |
| Glu | 48.3 | 164.8 | 13.9 | 0 | 227.0 |
| Asp | 63.5 | 226.0 | 33.1 | 0 | 322.6 |
| Ile | 7.7 | 87.6 | 0 | 218.6 | 313.9 |
| Arg | 6.0 | 83.7 | 7.7 | 0 | 97.4 |

### Percentage labels drawn on Figures 1 and 2 — FIGURE_ONLY, NOT to be fitted

Read from a 200-dpi raster of printed page 2919 (`scratchpad/img/h95b-3.png`) and cross-checked against
the class subtotals above (each bar height matches its subtotal, which is how the pot-to-label mapping
was fixed). Caption: "The numbers on the tops of the columns show the percent contributions of each
tested amino acid."

| pot | Figure 1 (pyridines) | Figure 2 (pyrroles) |
|---|---:|---:|
| Gln | 54 % | 18 % |
| Glu | 33 % | 22 % |
| Asn | 56 % | 32 % |
| Asp | 62 % | 30 % |
| Lys | 35 % | 34 % |
| Arg | 50 % | 26 % |
| Phe | 46 % | 17 % |
| Ile | 21 % | 17 % |

Figure 3 (oxazoles, y-axis 0-40) and Figure 4 (all nitrogen-containing compounds, y-axis 0-3000) carry
the same kind of labels and were not rastered; nothing here depends on them.

### Statements printed in prose — these are NOT figure-borne

- Pyrroles: "If we further examine the relative contributions of the tested amino acids to the formation
  of 1-methylpyrrole and 1-methylpyrrole-2-carboxaldehyde, **more than 90 % of the nitrogen atoms were
  from glycine**." (The paper's reason: glycine, having no side chain, forms the N-methyl compounds
  directly; the others need an extra cleavage.)
- Oxazoles: "**more than 80 % of the nitrogen atoms in the 4,5-dimethyloxazole ring come from glycine**",
  taken by the authors as evidence for the direct glycine + diacetyl route over the Strecker route.
- Benzonitrile: "**About 96 % of the nitrogen atoms in benzonitrile were contributed from
  phenylalanine.** This implies that benzonitrile is mainly the direct degradation product of
  phenylalanine."
- Pyrroles, overall: "The overall contributions of the tested amino acids were **below 35 %** in
  generating pyrroles when they competed with glycine ... **glycine is superior to other amino acids in
  the production of pyrroles**."
- Pyridines: "Aspartic acid had the highest contribution and isoleucine the lowest contribution in the
  formation of pyridines. Aspartic acid also generated the highest quantity of pyridines, while arginine
  had the lowest quantity of pyridine in the presence of glycine."
- Across all nitrogen compounds (Figure 4): "**Glutamic acid was the lowest contributor, while
  asparagine was the highest contributor**, to flavor formation among the tested amino acids in the
  presence of labeled glycine"; glycine's own yield is highest in the lysine pot and lowest in the
  arginine pot — "lysine acts as a synergist ... the arginine could depress specifically pyrroles,
  oxazoles, and pyridines at the expense of other products."
- The three 2-methylbutylamines appear **only** in the isoleucine pot and are attributed to
  2-methylbutanal (isoleucine's Strecker aldehyde) condensing with an amino group; 2,6-diphenylpyridine
  and benzonitrile appear only in the phenylalanine pot.
- Mechanism claimed for pyrroles: amino acid + 3-deoxyhexosone via Strecker degradation, then
  dehydration and ring closure (Kato & Fujimaki 1968); or furans + amines (Rizzi 1974). For pyridines:
  condensation of aldehydes/ketones with ammonia — which the authors then argue against, since the
  amide-bearing amino acids did not give more pyridines than their acids.

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`, 75 ids): **not one of the 25 compounds in Table 1 has an id.** No
pyridine, no pyrrole, no oxazole, no amine, no benzonitrile is registered. `2_acetylfuran` and
`furfural` are the nearest registered furan-side neighbours; `2_acetyl_1_pyrroline` is a pyrroline, not
one of these pyrroles.

**The pot for every row below:** glucose + 15N-glycine + one tested amino acid, equimolar (2.66 µmol
each as printed, flag 4), on 20 g of wheat starch rehydrated to 12-14 % moisture, the pre-drying
solution at pH 7, **180 C for 1 h**, one time point, no replicates stated, Ctrl carrying half the amine.

### The identity ratios this paper actually supports

| ratio | value | unit | what it is a ratio OF | conditions | source location | evidence class |
|---|---:|---|---|---|---|---|
| glycine's share of the nitrogen in **1-methylpyrrole** and **1-methylpyrrole-2-carboxaldehyde** | >90 (glycine); <10 (tested AA) | % of ring N | an **isotope partition** in the same pot | as above | text, p. 2918 | within_study_ratio (printed in prose) |
| glycine's share of the nitrogen in **4,5-dimethyloxazole** | >80 | % of ring N | isotope partition | " | text, p. 2920 | within_study_ratio (printed in prose) |
| phenylalanine's share of the nitrogen in **benzonitrile** | ~96 | % of ring N | isotope partition | " | text, p. 2920 | within_study_ratio (printed in prose) |
| tested amino acid's share of **pyrrole** nitrogen, all eight | <35 | % | isotope partition, ceiling only | " | text, p. 2918 | within_study_ratio (a bound, printed) |
| tested amino acid's share of **pyridine** and **pyrrole** nitrogen, pot by pot | not typed — see §3, Figures 1 and 2 | % | the per-pot partition the identity layer wants | " | Figures 1-4 labels | **figure_only** |
| total N-heterocycle yield, Gly+X pot ÷ Gly-only pot | Gln 1.981; Lys 2.894; Asn 1.352; Phe 2.179; Glu 2.643; Asp 3.753; Ile 3.654; Arg 1.134 | – | a **yield** (µg/g glucose), pot against pot | " | Table 1 totals (mine) | within_study_ratio |
| pyrrole subtotal, Gly+X ÷ Gly-only | Gln 1.91; Lys 3.40; Asn 1.53; Phe 2.71; Glu 2.82; Asp 3.87; Ile 1.50; Arg 1.43 | – | yield, by class | " | Table 1 (mine) | within_study_ratio |
| pyridine subtotal, Gly+X ÷ Gly-only | Gln 1.78; Lys 1.21; Asn 0.96; Phe 0.96; Glu 1.76; Asp 2.31; Ile 0.28; Arg 0.22 | – | yield, by class | " | Table 1 (mine) | within_study_ratio |
| ordering of pot totals | Asp > Ile > Lys > Glu > Phe > Gln > Asn > Arg > Ctrl | – | ordinal | " | Table 1 totals | within_study_ratio (ordinal) |
| individual compound yields, all 25 x 9 cells | see Table 1 | µg per g of glucose (table head says mg — flag 3) | absolute amount, internal-standard-normalised | " | Table 1 | level_only (no response factors) |
| any rate, barrier or time course | none | – | – | – | – | not reported |

### What this adds to, and takes away from, an identity ratio layer

**The reference amino acid is glycine**, in all nine pots, as the layer requires. On the printed pot
totals **the span is 1.13 to 3.75** (Arg to Asp), and on the pyrrole class alone 1.43 to 3.87 — narrower
than the companion's pyrazine span of 0.84 to 4.34, and **ordered differently**: aspartic acid and
isoleucine lead here, lysine and phenylalanine led there, and arginine, which was the only pot below the
glycine control for pyrazines (0.84), is above it for the N-heterocycles (1.13). Two product classes
from the same nine beds do not rank the amino acids the same way.

That disagreement is the useful finding, and the three prose percentages explain it. In one pot the
nitrogen partition is **product-specific**: glycine takes >90 % of the N-methylpyrrole nitrogen and
>80 % of the dimethyloxazole nitrogen while the competitor takes 21-62 % of the pyridine nitrogen and
its own side chain shows up as its own products (isoleucine's three 2-methylbutylamines, 218.6 of its
313.9 total; phenylalanine's benzonitrile at ~96 % phenylalanine nitrogen and its 2,6-diphenylpyridine).
A layer that gives each amino acid **one** scalar ratio against glycine, applied to every product of a
shared dicarbonyl pool, is contradicted by this paper's own numbers: which product forms depends on
which amine is available, not only on how much. If such a layer is fitted anyway, that is a declared
approximation and this paper is the citation for its size.

**And, as with the companion, the medium is wrong for `FROZEN_B18`.** 180 C at 12-14 % moisture with no
buffer, against fed dicarbonyls in water at 100-120 C and initial pH 8. There is no temperature ladder,
no time course, no moisture ladder, and the reacting pH is unknown, so nothing here can be placed on
`PYRAZINE_PH_STEPS`' axis or given a barrier.

## 5. Flags

1. **This is not an independent study.** Same pots, same charges, same day's chemistry as
   `hwang1995_extraction.md`; the authors say so. Fitting both papers' pot ratios as separate rows
   double-counts one experiment. They are two analyses of nine beds.
2. **Two printed totals do not equal their columns.** Lys sums to 248.7 against a printed 248.57 (the
   only total carrying two decimals in the table), and Asp to 322.6 against 322.4. The other seven close
   exactly. Either a cell is misprinted or the totals were computed from unrounded values; use the
   column sums, and do not read either total as a checksum.
3. **The unit in the table head is wrong.** "mg/g of glucose" in Table 1 against "µg/g glucose" on the
   y-axes of Figures 1-3. The subtotals reproduce the bar heights on the µg axes, so µg/g is the unit.
   Taken literally the table head would make these minor heterocycles a thousand times more abundant
   than the pyrazines of the companion paper, which is absurd.
4. **The charge is internally inconsistent by three orders of magnitude**, exactly as in the companion:
   2.66 µmol of glucose is 0.479 mg in 20 g of starch, which cannot support the reported yields;
   2.66 mmol can. The paper prints µmol. No absolute concentration derived from this paper is safe; the
   µg-per-g-of-glucose basis and every ratio within it are unaffected.
5. **No response factors, no calibration, no LOD, no replicates, no error bars**, one deuterated toluene
   internal standard for 25 analytes. Comparing tris(2-methylbutyl)amine with pyridine across a row is
   comparing responses, not moles.
6. **The percentages that matter most are figure-borne.** The per-pot pyridine and pyrrole partitions
   exist only as labels on Figures 1-3 (recorded in §3, not typed in §4). **To request from the
   authors or ACS:** the W1/W2 abundance table behind Figures 1-4. Unlike the companion, this paper
   offers no supplementary material at all, so there is nothing to order — it must be asked for.
7. **The 15N label is on glycine's α-amine only**, and every nitrogen that is not that label counts as
   the tested amino acid's, side chains included. For lysine, arginine, asparagine and glutamine that is
   two nitrogens against glycine's one.
8. **One temperature, one time.** 180 C, 1 h, single point, no vessel description. No rate and no
   barrier can come from this paper, and formation cannot be separated from subsequent loss.
9. **The reacting pH is unknown** (pH 7 was set in water before freeze-drying; the bed reacted at
   12-14 % moisture with no buffer).
10. **Registry gap, total:** none of the 25 compounds has a `compounds.yml` id. If any of this paper's
    products is ever to be answered for, the pyrroles (1-methyl-2-acetylpyrrole is the largest single
    number in the table at 131.5 in the Asp pot for 1-methylpyrrole) and 2-acetylpyridine (present in
    all nine pots, the one compound with a complete row) would be the first to register.

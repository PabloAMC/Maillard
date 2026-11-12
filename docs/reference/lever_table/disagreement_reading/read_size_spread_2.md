# Size disagreements, part 2: why the laboratories' effect sizes differ about 5x

Reader: one agent, 2026-11-12. Literature and keyed contrast files only; no engine output. Effect sign convention
checked against the keyed arms: for categorical cells the effect is log10(none / added amino acid), so negative means
"adding the amino acid raises the odorant". Time slopes are dex per doubling of time, temperature slopes dex per 10 C,
pH slopes dex per pH unit. Every keyed number quoted below was recomputed from the keyed arms and matched the cell
table.

What was read:

| source (PDF) | read |
|---|---|
| hwang1995 | in full (text layer + Table 1 as a page image, p. 181) |
| martin2001 | in full (text layer, Tables 2-5) |
| park2022 (park2022mealworm) | in part: Methods, Table 2 pyrazine rows and sums, pyrazine discussion |
| wei2019 (wei2019flax) | in part: Methods, Table 1 (all class rows), nitrogen-compound discussion |
| zhou2023c | in part: Methods, Fig. 1 and Fig. 5 / Table 1 as page images (pp. 3, 7) |
| shibamoto1977 | in part: pp. 144-147 as page images (Table I conditions, Table II, discussion) |
| vanlancker2012 | in part: Methods, Table 10 and its discussion |
| schirle1992 | in full text layer (Methods, Tables 1-2); Figs. 3-5 not viewed |
| chan1994b (file chan2005.pdf) | in full (text + Figs. 1-3 as page images) |
| hofmann2000 | in part: pH-related passages only; Fig. 4 not viewed |
| zamora2006, zamora2013, zamora2015 | not opened; keyed (digitised) values taken as is |

---

## 1. Pyrazines x amino acid identity (none vs glutamine): spread 0.78 dex

| lab | effect | the "none" pot | the added arm |
|---|---|---|---|
| hwang1995 | -0.20 (598.6 vs 958.3, Table 1 totals) | glucose + 15N-glycine, equimolar (2.66 "umol" each, as printed) on 20 g wheat starch, pH adjusted to 7, freeze-dried, 12-14 % moisture, 180 C 1 h | + glutamine, equimolar to glucose |
| martin2001 | -1.23 (G 16 vs GGn 272, Table 4) | potato slices leached of sugars AND amino acids (hot water + 50 % ethanol), infused with glucose only, fried 2 min | + glutamine 82.3 mg/100 g (0.56 mmol, about 0.4 mol per mol glucose taken up, 258.7 mg = 1.44 mmol) |
| park2022mealworm | -0.24 (0.729 vs 1.280, Table 2 sum) | 2 g mealworm protein + 0.1 g xylose in 2 mL water, 125 C 60 min | + 0.1 g glutamine (0.68 mmol, about 1 mol per mol xylose) |

**Explanation: (b), the reference arm's base mixture; secondary (d)/(e), a near-floor reference value in martin2001.**

- The two small effects both have an amine already in the reference pot: glycine at 1:1 to glucose (hwang1995) or
  2 g of protein against 0.1 g of sugar (park2022mealworm; sugar-limited, amine in large excess). Adding glutamine is
  a second amine competing for the same dicarbonyls: +0.20 to +0.24 dex.
- martin2001's reference has had its amino acids removed. The authors' own arms show the amine is the bottleneck:
  glutamine without glucose (Gn) gives total pyrazines 4, glucose without amino acid (G) gives 16, both together 272
  (Table 4). The residual 16 in G is attributed by the authors to incomplete leaching (2- and 3-methylbutanal found in
  the water blank, p. 3887). Adding the only real amine source to an amine-starved pot is a different lever from
  adding a second amine to an amine-replete pot.
- The dose relative to sugar does not explain it: martin2001's glutamine dose is the LOWEST of the three relative to
  sugar (about 0.4 vs 1 mol/mol) yet gives the largest effect.
- Floor: martin2001 prints relative peak-area units to two significant figures (CV < 25 %); G values of 16 (total) are
  readable, so the -1.23 is not a rounding artefact (rounding moves it by < 0.05 dex).
- Time (2 min frying vs 60 min) and matrix (potato vs starch vs insect protein) also differ, but the base-mixture
  account predicts the direction and size without them.

Would confirm: a martin2001-type reference that keeps the potato's own amino acids (their untreated U chips, total
359), or hwang1995's design without glycine; a study adding glutamine to a glucose pot with and without a background
amine. Would refute: an amine-free reference giving only about +0.2 dex.

**Confidence: high** for (b). Keying: hwang1995 ctrl column (lost from the text layer) checked on the page image:
598.63, 146.77, 45.68, 10.40 as keyed. martin2001 and park2022mealworm as keyed. Side note, no effect on ratios:
hwang1995's "2.66 umol" of glucose against 20 g starch is implausibly small (probably mmol); keyed "as printed".

## 2. Pyrazines x time (T 100-130 C): spread 0.77 dex

| lab | slope (dex/doubling) | window | system |
|---|---|---|---|
| wei2019flax | +0.98 (110 C), +0.36 (120 C) | 90 -> 150 min | flaxseed protein hydrolysate 2 g + xylose 1 g + cysteine, initial pH 7.5 (presumed), HS-SPME |
| zhou2023c | +0.11 | GSH pre-heated 30 -> 120 min, then 60 min with ARP at 120 C | Ala-xylose Amadori product + glutathione, pH 7.5 initial |
| shibamoto1977 | 0.00 (open), +0.11 (open, propyl gallate), +0.28 (closed), +0.22 (closed CuCl2), +0.23 (closed NaOH) | 2 h -> 18 h, 100 C | 1 M glucose + 8 M NH4OH |

**Explanation: (a) the window sits at a different stage of the reaction; (b) open vs closed vessel; (e) one
contrast's lever is not heating time; (d) unreplicated single SPME values in wei2019flax.**

- (a) wei2019flax itself shows the spread inside one laboratory: the same 90 -> 150 min window gives +0.97, +0.35 and
  +0.12 dex/doubling at 110, 120 and 130 C (Table 1: 0.010 / 0.028 / 0.052; 0.041 / 0.067 / 0.074;
  0.084 / 0.102 / 0.103 ug/g). At 110 C the window is the early rise (thiophenes and thiazoles not yet detected at
  110 C, p. 91); at 130 C it is near a plateau. shibamoto1977's 2 -> 18 h window at 100 C in 8 M ammonia is long and
  late, so its per-doubling slope is small. The 130 C row of wei2019flax (+0.12) is not in this cell but agrees with
  shibamoto1977 and zhou2023c.
- (b) In shibamoto1977, the open flask (atmospheric pressure, condenser) gives no time effect (0.940 -> 0.951 %
  yield, Exps. 13/14), the sealed ampoule gives +0.28 (0.681 -> 5.289 %, Exps. 15/16; Table II, p. 147). The authors
  attribute the large 18 h closed yields partly to measurement: "the use of a closed system, which does not give
  emulsions and thus facilitates quantitative recovery during methylene chloride extraction" (p. 145). So the open
  arm's flat slope is part chemistry (my inference: ammonia escapes an open flask at 100 C) and part a recovery bias
  that moves with time (d).
- (e) zhou2023c's keyed "time" contrast (zhou2023c-21, Fig. 5a, digitised 25 / 34 / 42 / 49 / 58 ug/L) is the time
  GSH was pre-heated alone before a fixed 60 min reaction with ARP (Methods, p. 14313; Fig. 5 caption). It is a
  precursor pre-treatment, not reaction time. By coincidence the true heating-time series (Fig. 1, sum of the three
  printed pyrazines in ARP+GSH: 17.4 at 30 min, 25.2 at 60, 21.2 at 90, 30.8 at 120) also gives about +0.12 per
  doubling over 30 -> 120 min (non-monotone), so the size conclusion does not change, but the row is in the wrong cell.
- (d) wei2019flax reports single values without SD. Three columns that the keying treats as the same nominal pot
  (0.15 g cysteine, pH 7.5, 120 C, 2 h: cysteine series 0.035, pH series 0.051, thermal series 0.067 ug/g) differ by
  0.28 dex if they are indeed the same pot (the paper does not print the heating conditions of the first two series).
  Spread over only 0.74 doublings, 0.1 - 0.15 dex per arm becomes about +-0.2 - 0.3 dex on the slope. This widens but
  does not create the 110 C value.

Would confirm (a): any lab's time series at two temperatures showing the per-doubling slope fall as temperature
rises; a covariate for extent of reaction (e.g. browning A420 or sugar consumed) would collapse the slopes.

**Confidence: medium - high** (a and b are visible inside single papers; the open-vessel mechanism is my inference).
Keying: values in all three sources verified (wei2019flax Table 1 text layer, shibamoto1977 Table II page image,
zhou2023c Figs. 1 and 5 page images). **Flag: zhou2023c-21 is keyed as lever class "time" but the lever is GSH
pre-degradation time.**

## 3. Furfural x temperature (time > 120 min): spread 0.75 dex

| lab | slope (dex/10 C) | temperatures | system and method |
|---|---|---|---|
| vanlancker2012 | +0.03 (glycine), +0.24 (Gly-Gly) | 150 -> 180 C (100 and 130 C not detected) | glucose + glycine 2 mmol or diglycine 1 mmol, unbuffered initial pH 8, 2 h; HS-SPME, peak area (Table 10) |
| schirle1992 | +0.75 | 100 -> 120 C, amounts at 10 h | glucose 0.5 mol + cysteine 0.1 mol in 400 mL, stirred Parr reactor; solvent extraction, AED carbon channel vs octane IS |

**Explanation: (a) different temperature segments; (e) a two-arm slope with the non-detects dropped; (b) the amine
(cysteine vs glycine vs diglycine); (d) possible SPME competition at 180 C. Low-medium on mechanism.**

- (a) schirle1992's own rate constants flatten with temperature: zero-order k = 1.692e3, 5.030e4, 6.829e5 UNA/h at
  100, 120, 150 C (Table 1) give +0.74 dex/10 C over 100 -> 120 and +0.38 over 120 -> 150; its Arrhenius fit
  (Table 2, Ea 35.2 kcal/mol, r = 0.94) is curved. Continuing that trend to 150 -> 180 C gives roughly 0.2 - 0.3,
  which brackets vanlancker2012's diglycine value. At 150 - 180 C for 2 h furfural is also consumed (vanlancker2012's
  glycine arm at 180 C shows pyrazines and 5-methylfurfural rising 1 - 1.5 dex while furfural barely moves: 20.6 ->
  25.4). An aside on schirle1992: the 150 C extract was concentrated 150x against 300x at 100 and 120 C (Methods);
  if the areas were not corrected, the true k at 150 C is twice the printed value and the curvature is smaller. The
  keying already excludes 150 C from the amount contrasts.
- (e) vanlancker2012 has nd at 100 and 130 C and 20.6 at 150 C. The keyed slope uses only 150 -> 180 C, so it drops
  the steep nd -> detected rise between 130 and 150 C. That rise cannot be sized without a detection limit, and none
  is printed.
- (b) Inside vanlancker2012 the amine form changes the slope by 0.21 dex (glycine +0.03, diglycine +0.24). The
  authors chose the comparison because diglycine favours furfurals at high temperature (p. 4706). schirle1992 uses
  cysteine, a weak browning amino acid at a 5:1 sugar excess.
- (d) At 180 C the glycine pot has a much larger headspace load than the diglycine pot (for example trimethylpyrazine
  405.9 vs 9.9, 2,5-dimethylpyrazine 210.6 vs 54.2). Fibre competition would depress furfural in the glycine arm
  more, and that arm has the flatter slope. This fits the documented 0.1 - 0.45 dex lever-linked fibre bias but is
  not tested.

Would confirm: schirle1992-type kinetics at 150 - 180 C, or vanlancker2012 at 110 - 130 C with a detection limit;
an SPME dilution check of the 180 C glycine pot.

**Confidence: medium** for (a) + (e); low for the split between (b) and (d). Keying: Table 10 values and schirle1992
Table 1 constants as keyed. **Flag: vanlancker2012-226/227 are two-arm slopes; the nd arms at 100 and 130 C were
dropped (no log of zero).**

## 4. 2,3-Dimethylpyrazine x amino acid identity (none vs asparagine): spread 0.72 dex

hwang1995 -0.24 (45.68 vs 79.64) vs martin2001 -1.28 (G 1 vs GA 19, Table 4).

**Explanation: (b), as in cell 1**: hwang1995's reference holds glycine; martin2001's is glucose on amino-acid-leached
potato. Asparagine is dosed at about 1.1 mol per mol glucose in martin2001 (210 mg = 1.59 mmol against 1.44 mmol)
and 1:1 in hwang1995, so dose does not explain the gap. **Secondary (d)/(e):** martin2001's G value is 1 relative
unit, the smallest printed unit, so rounding alone spans about 0.5 - 1.5 (effect -1.10 to -1.58). The value is
imprecise but still about 1 dex more than hwang1995.

**Confidence: high** for (b). Keying correct (hwang1995 checked on the page image).

## 5. Phenylacetaldehyde x pH (pH 6 - 7.5): spread 0.67 dex

| lab | slope (dex/pH) | design |
|---|---|---|
| hofmann2000 (x5) | -0.07 to -0.09 | Phe + one alpha-dicarbonyl or glucose, 1:1; keyed slope is the pH 5 -> 9 run of a 3/5/7/9 ladder |
| zamora2006 | -0.07 | Phe + 4,5-epoxy-2-decenal, 50 mM phosphate pH 6 -> 8, 37 C |
| zamora2013 | -0.11, -0.19 | Phe + 4-oxo-2-hexenal or -nonenal, 0.3 M phosphate pH 6 -> 8, 180 C |
| zamora2015 | -0.07 | Phe + 2,4-decadienal, phosphate, 140 C |
| chan1994b | -0.41 | glucose 0.5 mol + Met, Phe, Pro AND Leu, 0.075 mol each, in 400 mL 0.1 M phosphate, 95 C, 60 min; AED vs 4-methylthiazole |

**Explanation: (b) a co-reactant present only in chan1994b consumes phenylacetaldehyde faster at higher pH; (d) the
digitised pH 8 point is uncertain but not decisive; (c) possible, not testable.**

- (b) chan1994b's mixture contains leucine, so isovaleraldehyde is formed alongside phenylacetaldehyde. The authors
  explain the pH 6 maximum by the pH-dependent aldol condensation of the two aldehydes to 5-methyl-2-phenyl-2-hexenal,
  supported by more hexenal at pH 7 and 8 (pp. 133, Figs. 1-3). The figures show phenylacetaldehyde FALLING with time
  at every pH (pH 6: about 1.0e-3 -> 0.68e-3 mol from 20 to 60 min), and at pH 7 the unheated zero-time extract is
  already about 3e-3 and drops to about 0.35e-3 by 20 min while the hexenal bar rises to about 2e-3 (Fig. 2). So the
  60 min amount measures what is left after consumption, and consumption is what depends on pH. In the authors'
  ANCOVA, pH had no significant effect on phenylacetaldehyde's rate of formation (p. 134). All the other laboratories
  use phenylalanine as the only amino acid, so there is no second Strecker aldehyde.
- (d) The pH 8 value (1.05e-4, re-digitised; the two readers had 1.3e-4 and 2.5 - 3e-4) sits at the bottom of a
  linear 0 - 3e-3 axis. Over 6 -> 8 the slope ranges from -0.22 to -0.41 depending on that reading. The pH 6 -> 7
  step alone (6.8e-4 -> 1.8e-4, both well resolved) is -0.58 dex, so the steepness does not depend on the pH 8 point.
- (c) chan1994b holds 1.25 M glucose and 0.75 M amino acids in 0.1 M phosphate. Drift during heating is plausible;
  the authors measured pH after cooling but do not report it. zamora2013 uses 0.3 M buffer. Not testable from the
  papers.
- (a) is minor: hofmann2000's slope spans 5 -> 9 across its pH 5 maximum (local 7 -> 9 for 2-oxopropanal is about
  -0.13), and temperatures range from 37 to 180 C without changing the cluster at -0.07 to -0.19. In hofmann2000 the
  fall with pH is a shift of phenylalanine toward phenylacetic acid (aldehyde:acid 3:1 at pH 3, 1:5 at pH 9;
  abstract): one more consumption route, but a weaker one.

Would confirm: chan1994b's design without leucine, or with 5-methyl-2-phenyl-2-hexenal counted back as phenylacetaldehyde
equivalents (the pH slope should shrink toward -0.1).

**Confidence: medium** (the mechanism is the authors' own, with their hexenal data; not tested by removing leucine).
Keying: chan1994b Figs. 1-3 checked visually; values consistent with the plots. The `pdf:` field points to
chan2005.pdf, which is this paper (a file-name mismatch, not a value error). **Flag: zamora2015-110 keys the 6 -> 6.5
piece (3.9 -> 3.6) of a flat, non-monotone ladder (3.9 / 3.6 / 4.1 / 3.7 at pH 6 / 6.5 / 7 / 8). The whole-ladder
slope is about -0.01/unit, so the keyed -0.07 is scatter selected by the longest-monotone-run rule. All runs have
length 2 and the first was taken.** zamora2006/2013 values are digitised and were not re-checked against their PDFs.

## 6. 2-Methylpyrazine x amino acid identity (none vs glutamine): spread 0.67 dex

hwang1995 -0.28 (146.77 vs 279.85) vs martin2001 -1.30 (G 6 vs GGn 120). **Explanation (b), as in cells 1 and 4**:
the glycine-containing reference against the amine-leached reference. martin2001's G value of 6 is readable (rounding
< 0.04 dex). **Confidence: high.** Keying correct.

## 7. 2-Ethyl-3,6-dimethylpyrazine x amino acid identity (none vs arginine): spread 0.66 dex

| lab | effect | reference / arginine arm |
|---|---|---|
| hwang1995 | -0.11 (10.40 vs 13.25; printed as 3-ethyl-2,5-dimethylpyrazine, the same compound) | glucose + glycine; all arms adjusted to pH 7 before freeze-drying; arginine LOWERS total pyrazines (598.6 -> 500.4) |
| park2022mealworm | -0.48 (0.142 vs 0.432, Table 2 no. 11) | protein + xylose; 0.1 g arginine added to 2 mL water with no pH adjustment reported; arginine RAISES total pyrazines 0.729 -> 3.911 (+0.73 dex) |

**Explanation: (b) the base mixture, plus (c) arginine changing pH in park2022mealworm.**

- (b) As in cell 1: an amine-replete reference with glycine (hwang1995) against protein + xylose (park2022mealworm).
  park2022mealworm's glutamine arm raises this pyrazine only +0.16 dex, so the base mixture alone does not give its
  +0.48 for arginine.
- (c) park2022mealworm adds 0.1 g arginine (0.57 mmol, about 0.29 M in 2 mL) and reports no pH measurement or
  adjustment (Methods, p. 1648). If it was the free base (the form is not stated), the pot became alkaline, which
  speeds pyrazine formation. hwang1995 removed that route by adjusting every arm to pH 7, and the authors describe
  arginine as an inhibitor of glycine's pyrazine formation (p. 182). park2022mealworm's threonine arm (+1.03 dex)
  shows that other routes also matter there (threonine is itself a pyrazine source).

Would confirm: the pH of park2022mealworm's MPArg-RF before and after heating; an arginine arm adjusted to the
reference pH. This is the same pair and the same reading as the earlier trimethylpyrazine x arginine conflict.

**Confidence: medium.** Keying correct (both checked against the printed tables).

## 8. 2,5-Dimethylpyrazine x time (T 100-130 C): spread 0.66 dex

| lab | slope (dex/doubling) | window |
|---|---|---|
| zhou2023c | +0.77 | ARP + GSH, 120 C, 30 -> 120 min (Fig. 1: 0.226 / 3.75 / 5.39 / 7.97) |
| shibamoto1977 | +0.03, +0.04 (open), +0.32, +0.23, +0.22 (closed) | glucose + 8 M NH4OH, 100 C, 2 h -> 18 h (Table II: open 6.3 -> 8.1 mg; closed 4.0 -> 39.9 mg) |

**Explanation: (a) zhou2023c's slope is dominated by a lag-phase first step; (b) open vs closed vessel in
shibamoto1977; (d) the shibamoto1977 recovery remark (as in cell 2).**

- (a) In zhou2023c, 30 -> 60 min is +1.22 dex (0.226 -> 3.75). Over 60 -> 120 min the slope is +0.33 dex/doubling,
  the same as shibamoto1977's closed, additive-free arm (+0.32). 2,5-Dimethylpyrazine is absent from ARP alone at
  every time (Fig. 1) and appears only after GSH breaks down. Pyrazine and methylpyrazine in the same pots are already
  near their level at 30 min (8.92 -> 12.3 and 8.26 -> 10.5 by 120 min). So the steep first step is a lag in
  supplying this compound's precursors (mechanism: my reading of Figs. 1 and 3, low-medium). The 0.226 value may also
  sit near the method's LOQ; the LOQ list in the Methods cannot be assigned compound by compound.
- (b)/(d) shibamoto1977's open arms give almost no time effect (+0.03, +0.04) and its closed arms give +0.22 to
  +0.32. The paper's comment on emulsions and extraction recovery (p. 145) applies.

Would confirm: zhou2023c's 2,5-dimethylpyrazine at intermediate times (10 - 45 min). Restricting the zhou2023c slope
to the post-lag window should put both labs at 0.2 - 0.35.

**Confidence: medium - high** (the numbers are direct; the lag mechanism is inferred). Keying: zhou2023c Fig. 1
values verified on the page image; shibamoto1977 Table II verified on the page image.

---

## Summary

| cell | spread | explanation | confidence |
|---|---|---|---|
| Pyrazines x amino acid identity (none vs Gln) | 0.78 | b (amine-leached vs amine-replete reference); d/e minor | high |
| Pyrazines x time (T 100-130 C) | 0.77 | a (stage of reaction; T x time inside wei2019flax), b (open vs closed), e (zhou2023c lever is GSH pre-heating), d (unreplicated SPME; extraction recovery) | medium - high |
| Furfural x temperature (> 120 min) | 0.75 | a (100-120 vs 150-180 C; schirle1992's own k flattens), e (two-arm slope, nd arms dropped), b/d (amine form; possible fibre competition) | medium (b/d low) |
| 2,3-Dimethylpyrazine x amino acid identity (none vs Asn) | 0.72 | b; d/e (martin2001 reference = 1 unit) | high |
| Phenylacetaldehyde x pH (6-7.5) | 0.67 | b (leucine -> isovaleraldehyde aldol consumption in chan1994b only), d (pH 8 digitisation, not decisive), c possible; e (zamora2015 piece is scatter) | medium |
| 2-Methylpyrazine x amino acid identity (none vs Gln) | 0.67 | b | high |
| 2-Ethyl-3,6-dimethylpyrazine x amino acid identity (none vs Arg) | 0.66 | b + c (arginine, no pH adjustment in park2022mealworm) | medium |
| 2,5-Dimethylpyrazine x time (T 100-130 C) | 0.66 | a (zhou2023c 30 -> 60 min lag step; 60-120 matches shibamoto1977), b/d (open vessel; recovery) | medium - high |

Keying flags: no transcription errors found in the values checked. Structural flags: zhou2023c-21 (pre-treatment
time keyed as heating time); vanlancker2012-226/227 (two-arm slopes with nd arms dropped); zamora2015-110 (keyed piece
of a flat non-monotone ladder is scatter); chan1994b `pdf:` names chan2005.pdf (same paper). zamora2006/2013/2015 and
hofmann2000 Fig. 4 were not re-checked against their PDFs.

Overall pattern:
1. None of the 5x size gaps needs measurement noise. They come from what the reference state or window was: an
   amine-starved against an amine-replete reference pot (four amino-acid cells), where a time or temperature window
   sat on the reaction curve (three cells), and a co-reactant that consumes the odorant (phenylacetaldehyde).
2. Single papers reproduce most of the spread themselves: wei2019flax 0.97 / 0.35 / 0.12 at 110 / 120 / 130 C over
   the same window; shibamoto1977 open 0.00 against closed 0.28; schirle1992 0.74 against 0.38 per 10 C;
   vanlancker2012 glycine 0.03 against diglycine 0.24.
3. For the table this argues for three things: record the reference arm's base mixture (amine present or not) for
   categorical cells; record the extent of reaction (or the vessel, open or closed) alongside time and temperature
   slopes; check the lever class (pre-treatment time is not heating time). The longest-monotone-run rule matters
   here too (zamora2015-110, vanlancker2012 nd arms).

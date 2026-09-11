# Schieberle 1989 — EXTRACTION (ACS book chapter: 2-acetyl-1-pyrroline in wheat bread crust; aroma extract dilution analysis of a chemically leavened vs a yeast-leavened crust; 13-C labelling of proline and glucose; ground baker's yeast + sugar mixtures and dilute aqueous proline + phosphorylated sugar-degradation products, all quantified by stable isotope dilution assay; four tables of FD-factors, isotope ratios and absolute microgram yields; NO rate constant, NO temperature series, NO barrier)

### The origin paper for 2-acetyl-1-pyrroline as a bread-crust impact compound: it prints absolute microgram yields per pot and one bread-level number (34 ug/kg flour), it fixes the two-carbon acetyl group to the sugar and the ring to proline by labelling, and it shows the pyrroline + pyruvaldehyde couple giving 72 % of the volatile fraction — but every experiment is at a single temperature, so nothing here is a rate.

**Source on disk:** `data/articles/schieberle1989.pdf` (8 pp., ACS Symposium Series 409,
*Thermal Generation of Aromas*, Parliment, McGorrin and Ho (eds.), chapter 25, pp. 268-275).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/schieberle1989.txt`, 564 lines). The PDF is an ACS scan with an OCR text
layer stamped "Downloaded by EAST CAROLINA UNIV on November 3, 2016"; the body text is
letter-spaced by the OCR ("f o r m a t i o n") but readable, and **the digits in all four tables
came through clean**. Tables I-IV are re-typed in full below. **The OCR mangles the Greek mu
throughout**: it renders ug/kg as "pg/kg", "/ug/kg" and "j i g", and the header of Table III as
"(pg)". Every such instance is read as **microgram** — this is forced by the internal arithmetic
(section 4) and by the isotope-dilution method, and it is recorded as an OCR correction, not a
guess. Two further OCR artefacts: Table IV No. 4 prints "«C0.1" for "<0.1", and Table IV No. 5
prints the footnote marker `c)` run into the value as "Dihydroxyacetone phosphate )13 / 13.6".
There are no figures in this chapter and no supplementary material (1989). Repo status before this
dossier: the chapter is cited inside `chan1994b_extraction.md` (Chan & Reineccius name it as the
source for 2-acetyl-1-pyrroline) but has **no extraction dossier**.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of 2-Acetyl-1-pyrroline and Other Important Flavor Compounds in Wheat Bread Crust" |
| Author | Peter Schieberle — Deutsche Forschungsanstalt für Lebensmittelchemie, Lichtenbergstrasse 4, D-8046 Garching, Federal Republic of Germany (**sole author**) |
| Book | *Thermal Generation of Aromas*, T. H. Parliment, R. J. McGorrin, C.-T. Ho (eds.), ACS Symposium Series 409, American Chemical Society, Washington DC, 1989 |
| Chapter / pages | Chapter 25, pp. 268-275 |
| DOI | **10.1021/bk-1989-0409.ch025** (printed in the page furniture) |
| Publication date | 3 October 1989; "RECEIVED May 8, 1989" |
| Series code printed | `0097-6156/89/0409-0268$06.00/0` |
| Naming | **Acp = 2-acetyl-1-pyrroline** throughout (the author's own abbreviation); "pyruvaldehyde" = methylglyoxal = 2-oxopropanal; DHA = dihydroxyacetone; "dihydroxyacetone phosphate" = DHAP; CL-WB = chemically leavened wheat bread; YL-WB = yeast-leavened wheat bread; FD-factor = flavour dilution factor from aroma extract dilution analysis; "1-pyrroline" = the Strecker degradation product of proline |
| Companion chapter in the same book | "see Chapter by Schieberle and Grosch, this book" — cited twice, for the meaning of FD-factors and for the crust-vs-crumb Acp comparison; **not on disk** |
| Related dossiers on disk | `schieberle1998_extraction.md`, `schieberle2000_extraction.md`, `schieberle2005_extraction.md` (all later Schieberle work); `chan1994b_extraction.md` (cites this chapter for Acp); `hofmann1998b_extraction.md` (the methylglyoxal / dicarbonyl routes) |
| Who quotes it | Chan & Reineccius 1994 (`chan1994b_extraction.md`, reference list: "Schieberle, P. 1989 ... pp. 268-275", cited as the source for the proline route to Acp) |

## 1. Why it matters

The repository has **`2_acetyl_1_pyrroline` in `data/keys/compounds.yml`** (id at line 191) but no
dossier carrying a primary measurement of it. The only other on-disk source that touches the
compound is `chan1994b_extraction.md`, and what that chapter prints for Acp is a single averaged
activation energy (14.4 kcal/mol = 60.2 kJ/mol) with **no rate constant, no yield and no level** —
and it prints it from a four-amino-acid pot whose authors warn that secondary products in such a
pot are "very system dependent". So the proline route reaches the repository at present as one
barrier and nothing else.

**This chapter supplies the other half: absolute yields and a mechanism, but no kinetics.** It is
the paper that establishes (a) that Acp is the highest-FD-factor odorant of yeast-leavened wheat
bread crust and essentially absent from a chemically leavened one, (b) by 13-C labelling, that the
pyrroline ring comes from proline **minus its carboxyl carbon** and the acetyl group's **two**
carbons come from the sugar, and (c) that the immediate precursor pair is **1-pyrroline +
pyruvaldehyde (methylglyoxal)**, which together give 72 % of the volatile fraction where proline +
glucose gives none at all under the same boiling conditions.

That last point is the one that bears directly on how a model would have to be wired. The trunk
reaches methylglyoxal through the Amadori/deoxyosone supply steps (`martins2005b_extraction.md`,
`MARTINS_M4`: Amadori -> methylglyoxal, Ea 125.0 +/- 4.7 kJ/mol) and the pyrazine lane consumes
it with a measured second-order Strecker constant (`parameters_pyrazine.py`, `FROZEN_B18`). **This
chapter says the Acp route draws on the same methylglyoxal pool**, through a Strecker step on
proline that yields 1-pyrroline rather than an aldehyde, and then a condensation. If Acp is ever
added to the network it should be hung off the existing dicarbonyl pool and a proline Strecker
step, not off a direct proline + glucose reaction — the chapter's Table IV shows the direct
reaction giving **< 0.1 ug** with glucose, fructose or sucrose.

What this chapter does **not** give: any rate constant, any activation energy, any time series,
any second temperature for any experiment, and any concentration expressed per unit volume or mass
of reaction mixture (every yield is an absolute microgram figure per pot). It cannot fix a
magnitude for `k_strecker` or for any Acp-forming step, and it cannot be put on an Arrhenius axis.

## 2. Methods as they matter to a model

- **Breads.** The chemically leavened model bread: "A dough prepared from 500 g wheat flour (type
  550), 12 g glucono-delta-lactone, 10 g salt and ca. 270 ml of tap water was baked for **30 min at
  220 °C**" (Table I footnote a). The yeast-leavened comparator was prepared per reference 4
  (Schieberle & Grosch 1985) and its dough formula is **not printed here**. Crust volatiles were
  isolated "immediately after baking by extraction with dichloromethane and sublimation in vacuo".
- **Aroma extract dilution analysis.** Applied to the **acid-free** crust extract; 31 odorants
  detected; identification by MS/EI, MS/CI and retention data on two columns of different polarity
  against reference compounds; aroma quality assessed. Retention index measured on a **30 m x
  0.32 mm fused silica Supelcowax 10** (Supelchem, Germany) — a polar wax column, so the RI values
  in Table I are polar-column indices and are not comparable with DB-5 indices.
- **Quantification of Acp: stable isotope dilution assay** throughout (Table III footnote b,
  Table IV footnote b), on extracts obtained by **simultaneous distillation/extraction for 2 h**;
  for the one sub-0.1 ug entry of Table III, by diethyl ether extraction plus sublimation in vacuo
  instead. This is the strongest quantification method in this five-paper set — contrast Chan
  1994b's uncalibrated GC-AED response. **Absolute levels here are transportable in a way Chan's
  are not.**
- **Labelling experiments (Table II).** "4 mM of proline, 4 mM of glucose and 3 g of silica gel
  (10 % H2O) were heated for **30 min at 170 °C**." Isotope distribution of the Acp determined by
  **MS/CI** on extracts from simultaneous distillation/extraction of the reaction mixtures
  suspended in water. Experiment A: 1-13C-proline (carboxyl carbon labelled) + unlabelled glucose.
  Experiment B: unlabelled proline + U-13C-glucose. C: synthetic Acp reference.
- **Yeast experiments (Table III).** "30 g of yeast and 23 g of silica gel", **ground in a mortar
  for 10 min**, plus the stated sugar additive. Heating conditions for the ground mixtures are
  **not printed** except where footnote d says "Boil 2 hours". Controls: No. 8 omits the grinding;
  No. 10 adds 10 g ammonium sulfate to precipitate the yeast protein.
- **Dilute aqueous model (Table IV).** "A mixture of **2 mM (230 mg) of proline** and 2 mM of
  various sugar or phosphate ester mixtures dissolved in **400 ml of distilled water** was
  continuously steam distilled and extracted according to (7)." Phosphates used as their sodium or
  potassium salts. No. 10: 4 mM fructose 1,6-diphosphate pre-incubated 20 min at 25 °C with 200
  units of aldolase before heating.
- **THE UNIT "mM" IN THIS CHAPTER IS MILLIMOLES, NOT MILLIMOLAR.** The Table IV footnote settles
  it: "2 mM (230 mg) of proline" against proline's molar mass 115.13 g/mol gives 2.00 **mmol**
  (230.3 mg). In 400 mL that is 5 mmol/L. The same reading applies to the "4 mM" of Table II
  (proline and glucose on silica gel, where a molarity would be meaningless) and, less certainly,
  to the three back-flush experiments below. This is recorded as Flag 2 and every derived
  concentration in section 4 carries it.
- **The three back-flush experiments (running text, p. 273, no table).** All in 100 mL of
  **0.1 mol/L phosphate buffer, pH 7.0, boiled at back-flush for 2 hours**, volatiles isolated by
  ether extraction:
  1. proline (2 mM) + dihydroxyacetone (1 mM) — Acp = **0.1 % of the volatile fraction**;
  2. proline (2 mM) + pyruvaldehyde (0.1 mM) — Acp = **0.3 % of the volatile fraction**;
  3. **1-pyrroline (2 mM) + pyruvaldehyde (0.1 mM)** — Acp = **72 % of the volatile fraction and
     1140 ug**.
  Only the third has an absolute mass printed. "% of the volatile fraction" is a **composition
  share of an ether extract**, not a yield on a precursor, and the first two carry no absolute
  mass at all, so 0.1 % and 0.3 % cannot be converted (Flag 3).
- **Free proline determined by amino acid analysis**: flour 12 mg/kg, rising during fermentation to
  32 mg/kg in the dough; the baker's yeast used "more than 200 mg/kg yeast".
- **No temperature is varied anywhere in the chapter.** Baking 220 °C; labelling 170 °C; the
  aqueous models at the boil. Three isolated single-temperature points in three different matrices.

## 3. Tables re-typed

### Table I (p. 270). "Comparison of the Important Neutral/Basic Volatile Crust Flavor Compounds of a Chemically Leavened Wheat Bread (CL-WB)^a) With Those of a Yeast-leavened Wheat Bread (YL-WB)^b)"

RI per footnote c); FD-factors as printed.

| No. | Compound | Odor Description | RI | FD-factor CL-WB | FD-factor YL-WB |
|---:|---|---|---:|---:|---:|
| 1 | 2(E)-Nonenal | green, tallowy | 1508 | 512 | 256 |
| 2 | 2(E),4(E)-Decadienal | fatty, waxy | 1778 | 256 | 32 |
| 3 | 3-Methylbutanal | malty | 915 | 256 | 128 |
| 4 | 1-Octen-3-one | mushroom-like | 1281 | 128 | 32 |
| 5 | 2(Z)-Nonenal | green, fatty | 1480 | 128 | 64 |
| 6 | 2(E),4(E)-Nonadienal | fatty, waxy | 1669 | 128 | 16 |
| 7 | unknown | metallic | n.d. | 128 | <1 |
| 8 | Diacetyl | buttery | 968 | 64 | 64 |
| 9 | 2(E)-Octenal | fatty, nutty | 1393 | 64 | 4 |
| 10 | Unknown | boiled apple | 1786 | 64 | 16 |
| 11 | Hexanal | green | 1064 | 32 | 8 |
| 12 | 4(Z)-Heptenal | biscuit-like | 1214 | 32 | 32 |
| 13 | 1,5(Z)-Octadien-3-one ^d) | green, geranium-like | 1353 | 32 | 64 |
| 14 | 2(E),6(Z)-Nonadienal | cucumber-like | 1557 | 32 | 32 |
| 15 | Phenylacetaldehyde | honey-like | 1600 | 32 | 32 |
| 16 | **2-Acetyl-1-pyrroline** | roasty | 1299 | **4** | **512** |

Footnotes as printed: a) "A dough prepared from 500 g wheat flour (type 550), 12 g
glucono-delta-lactone, 10 g salt and ca. 270 ml of tap water was baked for 30 min at 220 °C."
b) "The dough was prepared according to (4)." c) "The retention index (RI) was determined on a
30 m x 0.32 mm fused silica column (Supelcowax 10; Supelchem; Germany)." d) "Structure not
established by MS."

**OCR note.** Row 13 reads "lr 5(Ζ)-0ctadien-3-one >" with the footnote marker run in; read as
1,5(Z)-octadien-3-one. Row 4 reads "l-0cten-3-one" (letter l and digit 0 for 1 and O). The table
is sorted by descending CL-WB FD-factor, which is why Acp sits last.

### Table II (p. 271). "Carbon Isotope Ratio in the Acp Formed from 1-13C-Proline/Unlabeled Glucose (A) or Unlabeled Proline/U-13C-Glucose (B)^a)"

| m/z | A | B | C (synthetic Acp) |
|---:|---:|---:|---:|
| 112 ^c) | 91.9 | 0.8 | 92.4 |
| 113 | 6.9 | 3.6 | 6.2 |
| 114 | 1.2 | 74.8 | 1.4 |
| 115 | <0.4 | 20.1 | <0.4 |
| 116 | <0.4 | 1.5 | <0.4 |

Column header as printed: "Isotope Distribution (%)^b)". Footnotes: a) "4 mM of proline, 4 mM of
glucose and 3 g of silica gel (10 % H2O) were heated for 30 min at 170 °C." b) "Values were
determined by mass chromatography (MS/CI)." c) "m/z 112 is the M+1-ion of unlabeled
2-acetyl-1-pyrroline."

The author's reading, printed in the text: column A "agreed with the data for synthetic Acp (C)",
and "because no upward shift of the M+1-ion (m/z 112 -> m/z 113) was observed, it can be concluded
that the carbon atom of the carboxyl group from proline is absent in the Acp formed". Column B's
major isotope at m/z 114 "shows a shift of 2", so "two carbon atoms in the Acp originate from
glucose"; and in the MS/EI (**data not shown**) the acetyl fragment moved from m/z 43 to m/z 45,
"indicating that both carbon atoms are derived from glucose".

### Table III (p. 272). "Formation of 2-Acetyl-1-Pyrroline from Ground Baker's Yeast Cells and Sugar Mixtures"

Header as printed: "Additive to 30 g of yeast and 23 g of silica gel ^a)" and "2-acetyl-1-pyrroline
^b) (ug)" — the unit is OCR'd "(pg)" and is read as **microgram** (see Source note and Flag 1).

| No. | Additive | 2-acetyl-1-pyrroline (ug) |
|---:|---|---:|
| 1 | none | <0.1 ^c) |
| 2 | none ^d) | 1.5 |
| 3 | 1 g Sucrose | 6.3 |
| 4 | 10 g Sucrose | 10.3 |
| 5 | 30 g Sucrose | 20.3 |
| 6 | 10 g Fructose | **26.9** |
| 7 | 10 g Glucose | 5.1 |
| 8 | 10 g Sucrose ^e) | 0.6 |
| 9 | 10 g Sorbitol | 1.2 |
| 10 | 10 g Fructose ^f) | 1.4 |

Footnotes as printed: a) "The mixture was ground in a mortar for 10 min." b) "The compound was
determined by an isotope dilution assay (2) in an extract which was obtained by simultaneous
distillation/extraction of the reaction mixture for 2 h." c) "The compound was determined by an
isotope dilution assay in an extract obtained by extraction of the reaction mixture with diethyl
ether and sublimation in vacuo (4)." d) "Boil 2 hours." e) "The grinding procedure was omitted."
f) "10 g of ammonium sulfate were added for protein precipitation."

**No heating temperature or time is printed for rows 1 and 3-10.** Only row 2 carries "Boil
2 hours" (footnote d), and it is attached to the *no-additive* row. Whether rows 3-10 were boiled
on the same schedule is not stated (Flag 4).

### Table IV (p. 274). "Formation of 2-Acetyl-1-Pyrroline (Acp) by Reaction of Proline and Sugars or Phosphorylated Sugar Degradation Products in Dilute Aqueous Solution^a)"

Header as printed: "Additive (2 mM)" and "Acp ^b) (ug)".

| No. | Additive | Acp (ug) |
|---:|---|---:|
| 1 | none | <0.1 |
| 2 | Glucose | <0.1 |
| 3 | Fructose | <0.1 |
| 4 | Sucrose | <0.1 |
| 5 | Dihydroxyacetone phosphate ^c) | **13.6** |
| 6 | 3-Phosphoglyceraldehyde | 1.0 |
| 7 | Phosphoenolpyruvate ^c) | 0.3 |
| 8 | D,L-alpha-glycerophosphate ^c) | <0.1 |
| 9 | Fructose 1,6-diphosphate ^c) | <0.1 |
| 10 | Fructose 1,6-diphosphate ^c) ^d) | **11.2** |

Footnotes as printed: a) "A mixture of 2 mM (230 mg) of proline and 2 mM of various sugar or
phosphate ester mixtures dissolved in 400 ml of distilled water was continuously steam distilled
and extracted according to (7)." b) "The Acp was determined by a stable isotope dilution assay."
c) "The compounds were used as their sodium or potassium salts." d) "4 mM of fructose
1,6-diphosphate were incubated for 20 min at 25 °C with 200 units of aldolase prior to heating."

**OCR note.** Row 4 prints "«C0.1"; read as <0.1. Row 5 prints the footnote marker run into the
value as ")13 / 13.6"; the value is 13.6 and the marker is c).

### Numbers in the running text (nothing else in this chapter is tabulated)

| quantity | value | where |
|---|---|---|
| free proline in the flour | 12 mg/kg | p. 271, from ref. 8 |
| free proline in the dough after fermentation | 32 mg/kg | p. 271, from ref. 8 |
| Acp produced by that bread | **34 ug/kg** | p. 271 |
| Acp when yeast replaced by a commercial leavening agent | **9.6 ug/kg flour** (down from 34) | p. 272 |
| proline addition levels tested | 120, 200, 500, 2000, 10000 mg/kg flour | p. 271 |
| marginal Acp gain from added proline | "the increase ranged only between **4.2 and 6.4 ug per 100 mg proline added**" | p. 271 |
| free proline in the baker's yeast used | "more than **200 mg/kg yeast**" | p. 273 |
| Acp from proline (2 mM) + DHA (1 mM), 100 mL 0.1 mol/L phosphate pH 7.0, 2 h back-flush | **0.1 % of the volatile fraction** (no absolute mass) | p. 273 |
| Acp from proline (2 mM) + pyruvaldehyde (0.1 mM), same conditions | **0.3 % of the volatile fraction** (no absolute mass) | p. 273 |
| Acp from **1-pyrroline (2 mM) + pyruvaldehyde (0.1 mM)**, same conditions | **72 % of the volatile fraction; 1140 ug** | p. 273 |
| odorants detected in the crust extract by AEDA | 31 | p. 269 |
| yeast low-molecular-weight fraction | compounds of MW < 1000, boiled and continuously extracted, "produced substantial amounts of Acp" | p. 273 (no number) |

### Arithmetic on the printed numbers (all mine)

Molar mass of 2-acetyl-1-pyrroline (C6H9NO) = 111.14 g/mol; proline 115.13 g/mol.

**1. The 1-pyrroline + pyruvaldehyde yield.** 1140 ug / 111.14 g/mol = **10.3 umol** of Acp.
Reading the charges as absolute millimoles (Flag 2): pyruvaldehyde 0.1 mmol = 100 umol, so the
yield on the **limiting dicarbonyl is 10.3 %**; on 1-pyrroline (2 mmol) it is 0.51 %. This is the
single most useful derived number in the chapter: it says the condensation of a Strecker-derived
1-pyrroline with methylglyoxal is an efficient reaction, not a trace one, which is what makes a
4-ug/kg-scale odorant reachable from a dicarbonyl pool at millimolar levels.

**2. The bread's conversion of proline to Acp.** 34 ug/kg / 111.14 = 0.306 umol/kg Acp against
32 mg/kg / 115.13 = 278 umol/kg free proline: **0.110 % of the dough's free proline appears as
Acp**.

**3. The marginal conversion of ADDED proline is ~20x worse.** 4.2 to 6.4 ug Acp per 100 mg
proline is 0.0378 to 0.0576 umol per 868.6 umol, i.e. **0.0044 % to 0.0066 %** — against the
0.110 % that the endogenous 32 mg/kg achieves. The ratio is **17x to 25x**. This is the
quantitative content of the author's conclusion that "another more effective way exists to form
Acp than the reaction between proline and sucrose or glucose": the yeast-borne route converts its
proline roughly twenty times more efficiently than bulk proline stirred into the flour. **The
chapter states the conclusion but does not compute this ratio; it is mine.**

**4. The yeast's share of the crust's Acp.** (34 - 9.6) / 34 = **72 %** of the Acp in this crust is
attributable to the yeast. (The coincidence with the 72 % volatile-fraction share in the
1-pyrroline experiment is arithmetically unrelated; do not conflate them.)

**5. Table IV's aqueous yields are minute.** No. 5 (DHAP), the best of them: 13.6 ug / 111.14 =
0.122 umol against 2 mmol proline = **0.0061 %**. So the dilute aqueous route is roughly 18x less
efficient than the bread and ~1700x less efficient than the pyrroline + pyruvaldehyde couple. The
ordering DHAP (13.6) > FDP+aldolase (11.2) >> 3-phosphoglyceraldehyde (1.0) > PEP (0.3) > all
sugars and glycerophosphate and untreated FDP (<0.1) is the chapter's mechanistic result.

**6. Sucrose dose-response in Table III is strongly sub-linear.** 1 g -> 6.3 ug; 10 g -> 10.3 ug;
30 g -> 20.3 ug. A 30x rise in sucrose gives a 3.2x rise in Acp; log-log slope = **0.34** (mine).
The author says only that "this increase was not proportional to the amount of sucrose added".
Fructose at 10 g (26.9 ug) beats sucrose at 30 g (20.3 ug) and beats glucose at 10 g (5.1 ug) by
**5.3x**.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Present: `2_acetyl_1_pyrroline`
(line 191), `3_methylbutanal`, `phenylacetaldehyde`, `2_3_butanedione` (alias "diacetyl"),
`hexanal`, `e_2_octenal`. **Absent: 2(E)-nonenal, 2(Z)-nonenal, 2(E),4(E)-decadienal,
2(E),4(E)-nonadienal, 2(E),6(Z)-nonadienal, 4(Z)-heptenal, 1-octen-3-one (only the -3-**ol** is
keyed, line 116 — a different compound), 1,5(Z)-octadien-3-one, methylglyoxal / pyruvaldehyde,
1-pyrroline, proline, glucose, fructose, sucrose, sorbitol, dihydroxyacetone, DHAP,
3-phosphoglyceraldehyde, phosphoenolpyruvate, glycerophosphate, fructose 1,6-diphosphate.**

Note against `chan1994b_extraction.md` Flag 11, which recorded that a `2_acetyl_1_pyrroline` key
would be needed before any proline row could be written: **the key now exists.**

| step | quantity | value | unit | conditions | order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| 1-pyrroline + pyruvaldehyde -> Acp | absolute Acp yield | **1140** | ug per pot | 1-pyrroline 2 mM + pyruvaldehyde 0.1 mM in 100 mL 0.1 mol/L phosphate, pH 7.0, back-flush boil 2 h, ether extraction | **none — single time, single temperature** | text p. 273 | level_only (isotope dilution; the strongest number in the chapter) |
| " | Acp share of the volatile fraction | 72 | % | as above | — | text p. 273 | level_only (a composition share, not a yield) |
| " | yield on the limiting dicarbonyl | **10.3** | % | as above, reading mM as mmol (Flag 2) | — | derived (mine) | derived_assumption |
| proline + pyruvaldehyde -> Acp | Acp share of the volatile fraction | 0.3 | % | proline 2 mM + pyruvaldehyde 0.1 mM, otherwise as above | — | text p. 273 | level_only (**no absolute mass printed — not convertible**) |
| proline + dihydroxyacetone -> Acp | Acp share of the volatile fraction | 0.1 | % | proline 2 mM + DHA 1 mM, otherwise as above | — | text p. 273 | level_only (no absolute mass) |
| proline + DHAP, dilute aqueous | absolute Acp yield | 13.6 | ug per pot | proline 2 mM + DHAP 2 mM in 400 mL water, continuous steam distillation/extraction | — | Table IV No. 5 | level_only |
| proline + fructose 1,6-diphosphate **after aldolase** | absolute Acp yield | 11.2 | ug per pot | 4 mM FDP + 200 U aldolase, 20 min at 25 °C, then heated with 2 mM proline | — | Table IV No. 10 | level_only |
| proline + 3-phosphoglyceraldehyde | absolute Acp yield | 1.0 | ug per pot | as Table IV | — | Table IV No. 6 | level_only |
| proline + phosphoenolpyruvate | absolute Acp yield | 0.3 | ug per pot | as Table IV | — | Table IV No. 7 | level_only |
| **proline + glucose / fructose / sucrose, dilute aqueous** | absolute Acp yield | **< 0.1** | ug per pot | as Table IV, at the boil | — | Table IV Nos. 2, 3, 4 | **null_result — the direct sugar route gives nothing under these conditions** |
| proline + glycerophosphate; untreated FDP; proline alone | absolute Acp yield | < 0.1 | ug per pot | as Table IV | — | Table IV Nos. 1, 8, 9 | null_result |
| ground yeast + fructose (10 g) | absolute Acp yield | **26.9** | ug per pot | 30 g yeast + 23 g silica gel ground 10 min; heating schedule not printed | — | Table III No. 6 | level_only (Flag 4: conditions incomplete) |
| ground yeast + sucrose (1 / 10 / 30 g) | absolute Acp yield | 6.3 / 10.3 / 20.3 | ug per pot | as above | — | Table III Nos. 3, 4, 5 | level_only |
| ground yeast + glucose (10 g) | absolute Acp yield | 5.1 | ug per pot | as above | — | Table III No. 7 | level_only |
| ground yeast, no additive, boiled 2 h | absolute Acp yield | 1.5 | ug per pot | as above | — | Table III No. 2 | level_only |
| ground yeast, no additive, not boiled | absolute Acp yield | < 0.1 | ug per pot | as above | — | Table III No. 1 | null_result |
| **grinding is required** | Acp with vs without grinding, 10 g sucrose | 10.3 vs **0.6** | ug per pot | Table III No. 4 vs No. 8 | — | Table III | within_study_ratio (**17x**, mine) |
| **yeast protein must stay in solution** | Acp with vs without ammonium sulfate, 10 g fructose | 26.9 vs **1.4** | ug per pot | Table III No. 6 vs No. 10 | — | Table III | within_study_ratio (**19x**, mine) |
| sugar alcohol control | sorbitol (10 g) | 1.2 | ug per pot | Table III No. 9 | — | Table III | null_result (against fructose 26.9) |
| bread crust, yeast-leavened | Acp | **34** | ug/kg | dough with 32 mg/kg free proline; 500 g flour type 550, baked per Table I footnotes | — | text p. 271 | measured_level (**the one food-matrix number in the chapter**) |
| bread crust, chemically leavened | Acp | **9.6** | ug/kg flour | yeast replaced by a commercial leavening agent | — | text p. 272 | measured_level |
| conversion of dough free proline to Acp | — | **0.110** | % (mol/mol) | as above | — | derived (mine) | derived_assumption |
| marginal conversion of ADDED proline to Acp | — | **0.0044 to 0.0066** | % (mol/mol) | proline added at 120-10000 mg/kg flour | — | derived (mine) from "4.2 to 6.4 ug per 100 mg" | derived_assumption |
| yeast's share of crust Acp | (34 - 9.6)/34 | **72** | % | this bread pair | — | derived (mine) | within_study_ratio |
| **carbon origin of Acp** | acetyl group = **2 carbons from the sugar**; pyrroline ring from proline **without** its carboxyl carbon | — | — | proline + glucose, 4 mmol each, 3 g silica gel 10 % H2O, 30 min at 170 °C | — | Table II + text pp. 270-271 | **structural — the mechanism result, and the most durable content of the chapter** |
| crust odorant ranking, yeast-leavened | Acp FD-factor **512**, the joint highest with 2(E)-nonenal | — | AEDA, Supelcowax 10 | — | Table I No. 16 | ordinal_only (an FD-factor is a dilution rank, **not** a concentration and **not** an odour threshold) |
| crust odorant ranking, chemically leavened | Acp FD-factor **4** | — | as above | — | Table I No. 16 | ordinal_only |
| any rate constant, any activation energy, any time series | — | **not printed** | — | — | — | whole chapter | — |

### Can anything here be put on the trunk's axes? Step by step.

**(a) No, for kinetics.** There is no k and no Ea, and no experiment is run at two temperatures, so
nothing in this chapter can enter a rate registry or be transported to `T_REF`. The chapter's Acp
barrier does not exist; the only Acp barrier in the corpus remains Chan 1994b's averaged
14.4 kcal/mol = 60.2 kJ/mol, which this chapter neither supports nor contradicts.

**(b) Yes, for network topology, and this is the transportable result.** The pairing of Table IV
(proline + glucose/fructose/sucrose -> **< 0.1 ug**) against the back-flush experiment
(1-pyrroline + pyruvaldehyde -> **72 % of the volatile fraction, 1140 ug**) is a clean within-study
contrast on the same analytical method. It licenses one structural statement for any future Acp
lane: **Acp is not formed from proline and a sugar directly at the boil; it is formed from the
Strecker product of proline and a small dicarbonyl.** The repository already carries a
methylglyoxal pool with a measured supply barrier (`martins2005b_extraction.md`) and a measured
second-order Strecker consumption (`FROZEN_B18`), so the hook exists.

**(c) Partly, for a level check.** The 34 ug/kg crust figure is a real isotope-dilution measurement
in a real food and could serve as a benchmark row **if** a bread lane ever existed — but the crust
is not a defined matrix in this repository (no water activity, no pH, no crust temperature profile
is printed here; only "30 min at 220 °C" oven), and the yeast contributes 72 % of it through a
biochemical route the trunk does not model at all. Carry it as context, not as a target.

**(d) No, for the Table I odorants.** FD-factors are dilution ranks from one AEDA on one extract.
They are not concentrations, they are not thresholds, and they are not comparable between the two
breads except as the ordering the author draws. Six of the sixteen compounds are lipid oxidation
products, which are outside the Maillard network entirely.

## 5. Flags

1. **The microgram sign is destroyed throughout the OCR.** "34 pg/kg", "9.6 /ug/kg", "1140 jig",
   "4.2 and 6.4 ug", and the Table III column header "(pg)" all denote **micrograms**. The reading
   is forced, not guessed: a picogram-per-kilogram bread level is four orders below what a 1989
   stable isotope dilution assay could reach, and the arithmetic in section 3 closes only in
   micrograms (a 1140 pg yield could not be 72 % of an ether-extracted volatile fraction). If the
   page images are ever consulted, verify this once and delete this flag.
2. **"mM" in this chapter means MILLIMOLES, not millimolar.** Table IV's own footnote proves it:
   "2 mM (230 mg) of proline" is 2.00 mmol at 115.13 g/mol. Table II's "4 mM of proline, 4 mM of
   glucose and 3 g of silica gel" must be the same (a molarity is meaningless on silica gel).
   **The three back-flush experiments in the running text are the uncertain case**: "Proline (2 mM)
   and dihydroxyacetone (DHA; 1 mM) were combined in 100 ml of 0.1 mol/l phosphate buffer" reads
   naturally as absolute mmol (giving 20 and 10 mmol/L) but could be millimolar. **Every derived
   percentage yield in section 3 and section 4 depends on this reading and must be re-derived if it
   is wrong.** The 1140 ug absolute mass does not depend on it.
3. **Two of the three back-flush experiments print no absolute mass.** "0.1 %" and "0.3 % of the
   volatile fraction" are shares of an unquantified ether extract. They establish the ordering
   (DHA < pyruvaldehyde << 1-pyrroline + pyruvaldehyde) and nothing more; they cannot be turned
   into yields, and the 720x and 240x ratios one might compute from the percentages are ratios of
   *composition shares of different extracts*, not of amounts.
4. **Table III's heating conditions are incomplete.** Only row 2 carries "Boil 2 hours"
   (footnote d), and it is attached to the no-additive control. Rows 3-10 print no temperature and
   no time. The Acp values of rows 3-10 are therefore levels from an **undescribed** thermal
   treatment. Request the schedule before any of Table III is used quantitatively.
5. **The yeast route is enzymatic, and the chapter says so.** The grinding control (10.3 vs 0.6 ug)
   and the ammonium sulfate control (26.9 vs 1.4 ug) are read by the author as showing that
   "glycolytic enzymes liberated by the grinding process are involved in the formation of
   precursors for Acp", and the active additives in Table IV are glycolytic intermediates (DHAP,
   3-phosphoglyceraldehyde) rather than sugars. **The dominant Acp route in bread crust is
   therefore biochemical up to the dicarbonyl and only then thermal.** A purely thermal Maillard
   model cannot reproduce the 34 ug/kg figure and should not be asked to.
6. **The two comparator breads differ in more than the leavening.** The chemically leavened dough
   formula is printed (Table I footnote a); the yeast-leavened one is not (footnote b defers to
   reference 4, Schieberle & Grosch 1985, **not on disk**). Flour type, salt, water and bake are
   only known for one of the pair, so the FD-factor comparison and the 34-vs-9.6 comparison are
   not fully specified controls.
7. **Fructose beats glucose by 5.3x in the yeast system (26.9 vs 5.1 ug, Table III Nos. 6 and 7)**
   and the author notes glucose "had the opposite effect". No mechanism is offered. If a sugar
   identity term is ever fitted anywhere in the repository, this is a data point of the right sign
   but the wrong matrix (a ground-yeast slurry with active enzymes, not a defined pot).
8. **One MS result is "data not shown"** — the MS/EI of the Acp from U-13C-glucose, on which the
   claim that *both* acetyl carbons come from the sugar rests. Table II carries only the MS/CI
   molecular-ion distribution, which fixes the total number of labelled carbons at two but not
   their position. The position claim is one un-shown spectrum.
9. **No pH, no water activity, no temperature series, no replication and no dispersion anywhere.**
   Not one value in the four tables carries a standard deviation, an n, or a replicate count. The
   dilute aqueous model is unbuffered distilled water (Table IV) while the back-flush experiments
   are at pH 7.0 phosphate; the two are not comparable to each other.
10. **What this chapter does not contain**: any rate constant; any activation energy; any time
    course; any second temperature; any concentration per unit volume or mass of reaction mixture;
    any odour threshold (the FD-factors are dilution ranks and the thresholds behind them are in
    the Schieberle & Grosch companion chapter, not here); any identification of the "unknown"
    odorants Nos. 7 and 10 of Table I; any supplementary material.
11. **What to request or fetch next**: (i) the **Schieberle & Grosch chapter in the same ACS volume
    409**, cited twice here, for the FD-factor definition and the crust-vs-crumb Acp comparison;
    (ii) Schieberle & Grosch 1985, Z. Lebensm. Unters. Forsch. 180:474-478 (ref. 4), for the
    yeast-leavened dough formula that Table I footnote b depends on; (iii) Schieberle 1988, Getreide
    Mehl Brot 44:334-335 (ref. 8), the source of the 12 -> 32 mg/kg proline figures and of the
    claim that Acp forms "only from proline"; (iv) Hodge, Mills & Fisher 1972, Cereal Science Today
    17:34-40 (ref. 11), for the pyruvaldehyde-catalysed Strecker degradation of proline to
    1-pyrroline — the step this chapter assumes rather than measures; (v) the heating schedule for
    Table III rows 3-10.
12. **Registry gaps against `data/keys/compounds.yml`**: `2_acetyl_1_pyrroline` is present (so
    `chan1994b_extraction.md` Flag 11 is now stale). **Absent and needed before any row of Tables I,
    III or IV could be written: methylglyoxal (pyruvaldehyde), 1-pyrroline, proline, and the five
    lipid-derived crust odorants of Table I (2(E)-nonenal, 2(Z)-nonenal, 2(E),4(E)-decadienal,
    2(E),4(E)-nonadienal, 2(E),6(Z)-nonadienal, 4(Z)-heptenal).** Note the trap at line 116: the
    registry keys `1_octen_3_ol`, which is **not** Table I's No. 4 `1-octen-3-one`.

# Modelling the Maillard reaction: an introduction

*For a reader who knows what the Maillard reaction is and nothing else. Seven short sections, each
led by a figure: the chemistry, what has been measured, how a model is built from it, how well this
one does, the one problem that stops it, what would settle it, and what to read next. Updated
2026-09-07; every figure is generated from the repository's own records, so they move with the
evidence. The step-by-step reaction trees are in the [appendix](REACTION_TREES.md).*

## 1. The chemistry, as the field draws it

![The Maillard reaction as the field draws it, and how well each part has been measured](../assets/thiol_sink/14_field_scheme.png)

Heat a reducing sugar with an amino acid and they condense to an Amadori compound. It decays along
two branches: one to the 3-deoxyosone, parent of HMF and furfural; the other to the 1-deoxyosone,
parent of the caramel furanones. Both branches also fragment into small dicarbonyls (glyoxal,
methylglyoxal, diacetyl), which strip amino acids of their nitrogen in the Strecker degradation and
give the amino-acid-specific aldehydes and, from those, pyrazines. Everything reactive ends in
melanoidins, the brown polymers. Cysteine adds a sulfur branch: it releases hydrogen sulfide, which
adds to the furanones and to furfural to make the two thiols that smell of cooked meat, MFT and FFT.
Asparagine has its own branch to acrylamide. Fats oxidise alongside, and their aldehydes cross into
the Maillard chemistry. This scheme has been stable since Hodge drew it in 1953; nobody disputes it.

The colours are the point. Dark blue steps have published rate constants at several temperatures;
green steps a rate or a yield at one temperature; grey steps are known as mechanisms from
isotope-labelling work but have no rate; the dashed red step has no published rate at all. Read that
way, the sugar and acrylamide branches are quantified, thiol formation is measured once, and thiol
removal, the step every meaty-aroma prediction ends on, has never been measured in a cooking pot.

## 2. Where the measurements sit

![Where the quantitative measurements sit](../assets/thiol_sink/15_field_coverage.png)

Each bar is one study and the temperatures it measured at. Almost every quantitative study sits
between 100 and 145 °C, in water; below 100 °C there are only milk-powder and storage studies, above
150 °C only dry glasses and roasting matrices. The thiol branch is measured at single temperatures by
different laboratories, so no laboratory has a temperature series of its own except two five-rung
grids published as figures. And the matrix changes the answer: the same dicarbonyls come out in the
opposite order in a sugar glass at 180 °C and in water at 120 °C.

## 3. How a kinetic model is built from that

![How a kinetic model works](../assets/thiol_sink/16_how_a_kinetic_model_works.png)

A kinetic model is the scheme above written as a list of steps, each with a rate constant, and
integrated over a cooking programme. Two ideas carry everything. First, steps in a row give curves in
time: an aroma compound is formed by one step and removed by the next, so it rises, peaks and falls,
and its level at the end of a cook is formation minus removal (left panel). Second, a rate constant
rises with temperature at a pace set by the step's activation energy: measured at 145 °C, a step with
a barrier of 60 kJ/mol is a fifth as fast at 100 °C, one with 160 kJ/mol a hundredth (right panel).
A constant measured at one temperature therefore says nothing about another unless the barrier is
measured too, and most of the thiol branch's constants come with no barrier.

![What the repository is made of](../assets/thiol_sink/09_repository_flow.png)

In this repository the papers are read into dossiers with every table re-typed, and three kinds of
evidence are kept apart on purpose: rate constants, which build the model; benchmark pots, end-of-cook
measurements the model is scored against and never tuned on; and directional claims ("this thiol
falls as pH rises") scored right or wrong. Every re-calibration is written down before it runs, with
the tests it must pass, and its outcome is recorded whether it passed or not. Of nearly 300 papers
registered, 34 supply a constant; the rest are the examiner, not the author. To run a prediction,
see the [QUICKSTART](QUICKSTART.md).

## 4. How well this model does

![The reaction paths the model carries, coloured by how well each is predicted](../assets/thiol_sink/00_map.png)

![Path by path scorecard](../assets/thiol_sink/08_path_scorecard.png)

Browning is predicted within a factor of 1.5 on a study the model never saw; HMF within 2 to 12; the
caramel furanone is 50 to 270 times off, and the small dicarbonyls come out in the wrong order in
water because their constants came from a sugar glass. Acrylamide reproduces the Leuven series it was
built from and other laboratories' pots within a median factor of 9. The meaty thiols land within 2 to
7 in the one laboratory, one temperature and one pH their constants came from, and are 10 to 500
times off everywhere else. Glucose plus cysteine has no route at all, so the model says "unknown"
rather than guessing.

## 5. The one problem that matters most: where the thiols go

![The model's own reference pot held at 100 C](../assets/thiol_sink/01_hofmann_pot_100C.png)

The reference pot, ribose and cysteine in buffer, the very system the constants came from, held at
100 °C for twelve hours. Blue is measured: the thiols rise forty-fold and are still rising at the end.
Green is the model as shipped: it peaks after an hour and then loses what it made. Orange and purple
are attempts to re-tune the removal on this series, with every existing bound respected and with one
bound relaxed; neither reaches the measurement, and both broke the 145 °C data the model was tuned
on. At the other end, a Beijing pot at 140 °C declines gently after an hour where the model loses a
thousandfold in two and a half, and at 168 °C the model predicts nothing at all where a pot measured
MFT halving between twenty and sixty minutes.

![A pot the model was never tuned on](../assets/thiol_sink/04_yiltirak_ladder.png)

Why it matters beyond the thiols: a Reading pot at 100 and 110 °C that no version of the model has
seen. As shipped, the model is 480 times too high for FFT; the re-tuned model, asked only to match the
100 °C series above, brings that to 18 without seeing this data. Most of what looked like laboratories
disagreeing with each other is the same removal problem. One removal step with one activation energy
cannot serve 100 °C and 145 °C at once; the removal needs a different structure, either reversible
(the disulfides give the thiol back) or saturating (it runs out of the partner it consumes). Two
smaller faults are known: the ring intermediate that carries the sulfur opens about ten times too
fast, and the formation steps have no pH dependence where the pots show a strong one.

## 6. What would settle it

**One experiment.** The reference pot (ribose 100 mmol/L, cysteine 33 mmol/L, 0.5 mol/L phosphate,
pH 5) in 20 mL vials with 5 mL of liquid, the volumes written down. Two temperatures, 100 °C for 0.5
to 12 hours and 140 °C for 5 to 120 minutes, one vial per time point, three replicates. Alongside,
the same buffer with MFT alone and FFT alone at 1 mg/L on the same grid, so removal is measured with
nothing forming. Measure the thiols by stable-isotope dilution and, in the same run, their
disulfides; the sugar and cysteine remaining; the pH. One arm at 100 °C under nitrogen. About 130
vials, two weeks of GC-MS. The shipped model predicts the fed thiol decays to zero and the pot peaks
after an hour; the alternative is a plateau with its disulfide and a pot still rising at twelve hours.
Either result decides the next version of the removal step.

**Without a laboratory.** Two published grids exist only as figures (Beijing 100 to 140 °C; a
five-temperature ladder from a second Chinese group); their numbers, obtainable from the authors,
would give removal its first tuning rows from a third laboratory. And the reversible or saturating
removal can be built and pre-registered against the existing series in a day.

## 7. Ten papers to read first

1. **Hodge 1953**, J. Agric. Food Chem. 1:928. The scheme everyone still draws.
2. **Martins & van Boekel 2005**, Food Chem. 90:257. The glucose-glycine cascade with rate constants at three temperatures; the template for every multiresponse model since.
3. **Hofmann & Schieberle 1998**, JAFC 46:235. Each intermediate fed on its own, the yield of MFT and FFT from each; the quantitative backbone of meaty-aroma chemistry.
4. **Cerny & Davidek 2003**, JAFC 51:2714. Labelled ribose with cysteine: which carbons become which sulfur compound.
5. **Whitfield & Mottram 1999 and 2001**, JAFC 47:1626 and 49:816. Norfuraneol fed with cysteine or H2S at pH 4.5 and 6.5: the pH switch in thiol formation.
6. **De Vleeschouwer et al. 2009**, Food Chem. 114:116. Acrylamide formation and elimination by multiresponse modelling; the best-parameterised branch in the field.
7. **Kocadağlı & Gökmen 2016**, JAFC 64:6333. Dicarbonyls and HMF in a glucose glass at 160 to 200 °C; the dry-side counterpart of the Wageningen work.
8. **Schieberle, Hofmann & Münch 2000**, ACS Symp. Ser. 756 ch. 10. The one published time series of the thiols at 100 °C.
9. **Hofmann & Schieberle 2000**, JAFC 48:4301. The Amadori compound's oxidative Strecker route: the only measured oxygen effect.
10. **Yiltirak et al. 2026**, Food Res. Int. The most recent thiol ladder, with the vessel and buffer stated, from an independent laboratory.

Nine of the ten are on disk in `data/articles` (Hodge 1953 is cited for its scheme); six have an
extraction dossier in `data/lit/extraction_dossiers`, and Cerny 2003, Whitfield 1999 and
De Vleeschouwer 2009 are read into the model's constants directly.

## Words used here

| word | meaning |
|---|---|
| pot | one published cooking experiment: stated ingredients, buffer, temperature, time, and a measured result |
| held out | a pot the model was never tuned on; every score in this document is on held-out pots |
| rate constant | how fast one step runs at a stated temperature |
| activation energy | the number that says how fast a rate constant rises with temperature |
| fed-intermediate experiment | one intermediate heated on its own; the cleanest source of a rate constant or a yield |
| fold error | measured over predicted, or the reverse, whichever is larger; 1 is perfect, 3 is the working threshold |
| directional claim | a published statement of direction, scored independently of absolute levels |
| removal, sink | the steps that take an aroma compound out of the pot once it has formed |
| multiresponse model | a kinetic model fitted to all measured species of one experiment at once |
| dossier | the repository's re-typed record of one paper |

## Sources

| figure or number | measurement | repository record |
|---|---|---|
| sections 1, 2 | the papers the repository has read | `data/lit/extraction_dossiers/`, `scripts/generators/build_thiol_sink_figures.py` |
| section 4 | the benchmark scorecard and the directional panel | `results/validation/core_panel_scores.json`, `core_directional_scores.json` |
| section 5, 100 °C pot | Schieberle, Hofmann & Münch 2000, Table IV | `schieberle2000_extraction.md` |
| section 5, 140 and 168 °C | Wang et al. 2022 (Flavour Fragr. J.); Liu et al. 2023 (LWT) | `wang2022_extraction.md`, `liu2023b_extraction.md` |
| section 5, Reading ladder | Yiltirak et al. 2026, Food Res. Int. | `yiltirak2026_extraction.md` |
| the re-tuning attempts | the repository's pre-registrations and their outcomes | `results/validation/kinetic_core_b*_prereg.md`, `scripts/generators/WAVES.md` |

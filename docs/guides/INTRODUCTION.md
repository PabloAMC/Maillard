# Modelling the Maillard reaction: an introduction

*For a reader who knows what the Maillard reaction is and nothing else. Each section is one figure and
a few sentences. The figures are made from the repository's own records, so they change when the
evidence changes. Updated 2026-09-07. Step-by-step reaction trees: [appendix](REACTION_TREES.md).*

The chemistry is known. The rates are measured in one temperature window. The model works where the
rates were measured and fails where they were not. Each section below shows one part of that.

## 1. The chemistry

![The Maillard reaction as the field draws it, and how well each part has been measured](../assets/thiol_sink/14_field_scheme.png)

*The map of the reaction has not changed since 1953. The colours show how well each step has been
measured: only the sugar and acrylamide branches have rates at several temperatures.*

A sugar and an amino acid join to form an Amadori compound. It breaks down along two branches. One
gives HMF and furfural. The other gives the caramel furanones. Both also break into small fragments
(glyoxal, methylglyoxal, diacetyl). The fragments react with amino acids and make the aldehydes that
smell of each amino acid, and from those, pyrazines. Everything ends in melanoidins, the brown
polymers. Cysteine adds a sulfur branch: it releases hydrogen sulfide, which joins the furanones and
furfural to make the two thiols that smell of cooked meat, MFT and FFT. Asparagine has its own branch
to acrylamide. Fats oxidise alongside and their aldehydes join in.

## 2. Where the measurements are

![Where the quantitative measurements sit](../assets/thiol_sink/15_field_coverage.png)

*Almost every rate was measured between 100 and 145 °C in water. Below that there are only powders
and storage tests; above it only dry glasses and roasting. The thiol branch has no temperature series
from a single laboratory except two published as figures.*

The matrix matters as much as the temperature. The same small fragments come out in the opposite
order in a sugar glass at 180 °C and in water at 120 °C.

## 3. How a kinetic model works

![How a kinetic model works](../assets/thiol_sink/16_how_a_kinetic_model_works.png)

*Left: a compound is formed by one step and removed by the next, so it rises, peaks and falls. Right:
a rate measured at 145 °C may be five times slower at 100 °C, or a hundred times; only a second
measurement tells which.*

The model is the map in section 1 written as a list of steps, each with a rate, run forward in time
for a given recipe and cooking programme. Most of the thiol branch's rates were measured at one
temperature, so away from it they are guesses.

## 4. What the repository is made of

![What the repository is made of](../assets/thiol_sink/09_repository_flow.png)

*Papers are read into dossiers with every table re-typed. Three kinds of evidence are kept apart:
rates that build the model, benchmark pots it is scored on but never tuned on, and statements of
direction it is scored on. Every re-calibration is written down before it runs, pass or fail.*

![The papers this model rests on most](../assets/thiol_sink/17_papers_by_weight.png)

*Counted from the repository's own records: which papers the model was tuned on and which it is
scored against. One paper, Hofmann and Schieberle 1998, carries most of the sulfur branch.*

## 5. How well the model does

![The reaction paths the model carries, coloured by how well each is predicted](../assets/thiol_sink/00_map.png)

*Green: predicted within a factor of about 1.5 on pots the model never saw. Amber: right in shape, and
within a factor of 3 only in the laboratory the rates came from. Red: wrong by more than tenfold.
Grey dashed: no route in the model.*

![Path by path scorecard](../assets/thiol_sink/08_path_scorecard.png)

*The same five paths with what we have, how each does, and what is missing. Browning is predicted;
the small fragments come out in the wrong order in water; the meaty thiols are right only in one
laboratory at one temperature and pH; glucose with cysteine has no route at all.*

## 6. The one problem that matters most: where the thiols go

![The model's own reference pot held at 100 C](../assets/thiol_sink/01_hofmann_pot_100C.png)

*The reference pot, held at 100 °C for twelve hours. The pot keeps making thiols the whole time.
Every version of the model stops after an hour and then loses what it made. Re-tuning the removal
step on this series (orange, purple) does not reach the measurement and breaks the 145 °C data.*

![A pot the model was never tuned on](../assets/thiol_sink/04_yiltirak_ladder.png)

*A pot from another laboratory that no version of the model has seen. Weakening the removal step
alone cuts the error from 480 times to 18. Most of what looked like laboratories disagreeing is the
same removal problem.*

The model removes thiols far faster than any pot does, at 100 °C and at 140 °C alike. One removal
step with one temperature dependence cannot fit both, so the step needs a different form: either the
disulfides give the thiol back, or the removal stops when it runs out of the partner it needs. Two
smaller faults are also known: the ring intermediate that carries the sulfur opens about ten times
too fast, and the formation steps have no pH dependence where the pots show a strong one.

## 7. What is needed next

**One experiment.** The reference pot (ribose 100 mmol/L, cysteine 33 mmol/L, 0.5 mol/L phosphate,
pH 5) in 20 mL vials with 5 mL of liquid, volumes written down. Two temperatures, 100 °C for 0.5 to
12 hours and 140 °C for 5 to 120 minutes, one vial per time point, three replicates. Beside it, the
same buffer with MFT alone and FFT alone at 1 mg/L on the same grid, so removal is measured with
nothing forming. Measure the thiols by stable-isotope dilution and their disulfides in the same run,
plus the sugar and cysteine left and the pH. One arm at 100 °C under nitrogen. About 130 vials and two
weeks of GC-MS. The model as shipped predicts that the fed thiol decays to zero and that the pot peaks
after an hour. If instead the thiol levels off with its disulfide and the pot keeps rising, the
removal step is reversible. Either result decides the next version.

**Without a laboratory.** Two published data sets exist only as figures (a Beijing grid at 100 to
140 °C and a five-temperature ladder from another Chinese group). Their numbers, from the authors,
would give the removal step its first data from a third laboratory. And the reversible or saturating
removal step can be built and tested against the existing series in a day.

## Words used here

| word | meaning |
|---|---|
| pot | one published cooking experiment: stated ingredients, buffer, temperature, time, and a measured result |
| held out | a pot the model was never tuned on; every score here is on held-out pots |
| rate constant | how fast one step runs at a stated temperature |
| activation energy | the number that says how fast a rate rises with temperature |
| fed-intermediate experiment | one intermediate heated on its own; the cleanest source of a rate or a yield |
| fold error | measured over predicted, or the reverse, whichever is larger; 1 is perfect, 3 is the working threshold |
| removal, sink | the steps that take an aroma compound out of the pot once it has formed |
| dossier | the repository's re-typed record of one paper |

## Sources

| figure | measurement | repository record |
|---|---|---|
| sections 1, 2, 4 | the papers the repository has read | `data/lit/extraction_dossiers/`, `data/keys/papers.yml`, `scripts/generators/build_thiol_sink_figures.py` |
| section 5 | the benchmark scorecard and the directional panel | `results/validation/core_panel_scores.json`, `core_directional_scores.json` |
| section 6, 100 °C pot | Schieberle, Hofmann & Münch 2000, Table IV | `schieberle2000_extraction.md` |
| section 6, second pot | Yiltirak et al. 2026, Food Res. Int. | `yiltirak2026_extraction.md` |
| the re-tuning attempts | the repository's pre-registrations and outcomes | `results/validation/kinetic_core_b*_prereg.md`, `scripts/generators/WAVES.md` |
| the trees | the model's reaction lists | `scripts/generators/build_reaction_tree.py`, [REACTION_TREES.md](REACTION_TREES.md) |

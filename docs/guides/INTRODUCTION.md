# Modelling the Maillard reaction: an introduction

*For a reader who knows what the Maillard reaction is and nothing else. Each section is one or two
figures and a few sentences. The figures are made from the repository's own records, so they change
when the evidence changes. Updated 2026-09-08. Step-by-step reaction trees: [appendix](REACTION_TREES.md).
Every paper used: [sources](SOURCES.md).*

The chemistry is known. The rates are measured in one temperature window, mostly by one laboratory
per path. The model works where the rates were measured and fails where they were not.

1. [The chemistry](#1-the-chemistry)
2. [Where the measurements are, and what they can support](#2-where-the-measurements-are-and-what-they-can-support)
3. [How a kinetic model works](#3-how-a-kinetic-model-works)
4. [What the repository is made of](#4-what-the-repository-is-made-of)
5. [How well the model does](#5-how-well-the-model-does)
6. [How the model got here](#6-how-the-model-got-here)
7. [The one problem that matters most](#7-the-one-problem-that-matters-most-where-the-thiols-go)
8. [What is needed next](#8-what-is-needed-next)

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

## 2. Where the measurements are, and what they can support

![Where the quantitative measurements sit](../assets/thiol_sink/15_field_coverage.png)

*Almost every rate was measured between 100 and 145 °C in water. Below that there are only powders
and storage tests; above it only dry glasses and roasting. The thiol branch has no temperature series
from a single laboratory except two published as figures.*

![The same two thiols, five laboratories](../assets/thiol_sink/18_lab_spread.png)

*Cysteine with a pentose, or with its first product, heated in five laboratories: the two meaty
thiols come out a thousandfold apart. Recipes, vessels and calibration all differ, and no pot has
been repeated by a second laboratory. No model can be more accurate than this spread until it is
explained.*

![What the papers behind the test panel state](../assets/thiol_sink/19_what_papers_report.png)

*What each paper on the test panel tells us. Temperature, time and pH are always given. The vessel,
the headspace and the atmosphere, which set how much oxygen the pot sees, are missing from a third
to three quarters of them. Fitting harder cannot recover what was not written down.*

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
scored against. One paper, Hofmann and Schieberle 1998, holds up most of the sulfur branch.*

![One measurement's path through the repository](../assets/thiol_sink/20_one_number.png)

*One value followed from the paper to the figure that shows it. At every step there is a file, and
the file says where the value came from and what it may be used for. This one is marked "never used
for tuning", so the model's miss on it is a real test.*

## 5. How well the model does

![The same map as section 1, coloured by how well this model predicts each part](../assets/thiol_sink/00_map.png)

*The map from section 1 again. Boxes: how far the model's prediction is from the measurement, on pots it
was never tuned on. Arrows: whether the model knows the step's rate at several temperatures, at one, or
only as a range. Brown colour is right. HMF is close. The furanones, the thiols, acrylamide and the fat
aldehydes are off by more than tenfold on other laboratories' pots. Glucose with cysteine has no route.*

![Predicted against measured, every pot on the panel](../assets/thiol_sink/21_parity.png)

*Every scored pot. Filled points are from laboratories the rates did not come from. The misses
cluster by path: the thiols from other laboratories all sit too high, the fat aldehydes all too low,
acrylamide both ways. That is what missing steps look like, not noise.*

![Path by path scorecard](../assets/thiol_sink/08_path_scorecard.png)

*The same five paths with what we have, how each does, and what is missing. Browning is predicted;
the small fragments come out in the wrong order in water; the meaty thiols are right only in one
laboratory at one temperature and pH; glucose with cysteine has no route at all.*

![Why some rows get no number](../assets/thiol_sink/22_no_number.png)

*The rows the panel asks for and the model declines to answer. It says "no route" or "no such
compound" rather than guessing. Most are the thiols from glucose, and compounds from fat that the
model does not make.*

## 6. How the model got here

Every re-calibration was written down before it ran, with the test it had to pass. Most failed.
This is the list, in order, in plain words.

| what was tried | why | what happened | kept |
|---|---|---|---|
| the sugar path fitted to one glucose and glycine study at three temperatures | the only study with every step measured | browning on pots it never saw came out within a factor of about 1.5 | yes |
| the pentose and cysteine path fitted to one laboratory's fed-intermediate experiments at 145 °C, with that laboratory's pH series | the only step-by-step data for the thiols | fits that laboratory; the pH trend of the last step comes out backwards | yes |
| the asparagine path fitted to one laboratory's rate constants | the only rate constants for acrylamide | right shape; levels off by tenfold in other laboratories | yes |
| the fat path built from one product slate, with the rate assumed | no measured rate exists at cooking temperatures | flagged as an assumption on every answer | yes, flagged |
| HMF and the caramel furanone added to the sugar path | the panel measures them | HMF close; the furanone off by fiftyfold | yes |
| the thiol removal steps given measured temperature dependence | the removal step decides the level | fits the reference laboratory | yes |
| end-of-cook levels taken out of the thiol fit; only rates, yields and ratios tuned | levels are what the model should predict, not what it is tuned on | glucose and fructose thiols fell to zero: they had no step-level support, and the model now says so instead of guessing | yes, the current thiol parameters |
| the temperature dependence split into two barriers, with a second laboratory's temperature ladder in the fit | one barrier for every formation step is too crude | neither barrier could be pinned by the data | no |
| oxygen made an input, with the vessel's headspace as a reservoir | the laboratories' vessels differ | the oxygen steps could not be pinned; no gap between laboratories closed | no; kept inert |
| water activity and pH terms on the sugar path, declared from measurements rather than fitted | the panel moves both | the water claims split; the pH claims agree | yes |
| the small dicarbonyls added from constants measured in a dry glass | the panel measures them | they come out in the wrong order in water | yes, flagged |
| water activity and pH terms on the acrylamide path, declared from one laboratory's series | the panel moves both | the pH claims agree | yes |
| thiol removal re-tuned on the twelve-hour 100 °C series and the intermediate's own decay | the model loses thiols too fast at low temperature | peaks at six hours where the pot keeps rising; breaks the 145 °C pots | no |

## 7. The one problem that matters most: where the thiols go

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
smaller faults are also known: the ring intermediate that holds the sulfur opens about ten times
too fast, and the formation steps have no pH dependence where the pots show a strong one.

## 8. What is needed next

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

**What the tool can do today.** It compares two recipes, predicts one, explains where a compound
comes from, scores your own measurements, calibrates itself to your laboratory from those
measurements without touching the shipped model, and lists what to measure next. It refuses questions the
evidence cannot answer instead of guessing, and when it refuses it says whether the chemistry has no
route or only no rate: a small layer of cited reaction rules, run over the model's species, lists the
steps the literature draws that the model does not have
([network_hypotheses.md](../../results/validation/network_hypotheses.md)). For the thiols, the same
list with the literature's numbers is [the sink table](../validation/thiol_sink_candidates.md). The
table of what the tool answers and what it refuses is at the top of the [quick start](QUICKSTART.md).

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

Every paper the model uses, alphabetically, with what was taken from each: [SOURCES.md](SOURCES.md),
made from the repository's records. Where each figure's numbers come from:

| figure | measurement | repository record |
|---|---|---|
| sections 1, 2, 4: the scheme, the coverage map, the papers by weight | the papers the repository has read | `data/lit/extraction_dossiers/`, `data/keys/papers.yml`, `scripts/generators/build_thiol_sink_figures.py` |
| section 2: the five laboratories | Hofmann & Schieberle 1998, Bolton 1994, Zhou 2023, Kang 2026, Yiltirak 2026 | the benchmark panel and the sulfur fit's row table |
| section 2: what papers state | the benchmark files' vessel and calibration fields | `data/benchmarks/`, `results/validation/core_panel_scores.json` |
| section 4: one number's path | Yiltirak 2026, Table S3 | `data/benchmarks/external_validation/maillard_path/`, `docs/validation/directional_claims_panel.yml` |
| section 5: the map, the parity plot, the scorecard, the rows with no number | the benchmark scorecard and the directional panel | `results/validation/core_panel_scores.json`, `core_directional_scores.json` |
| section 6: the history | the pre-registrations and their outcomes | `results/validation/kinetic_core_b*_prereg.md`, `kinetic_core_b*_ship_rule.md`, `scripts/generators/WAVES.md` |
| section 7, 100 °C pot | Schieberle, Hofmann & Münch 2000, Table IV | `schieberle2000_extraction.md` |
| section 7, second pot | Yiltirak et al. 2026, Food Res. Int. | `yiltirak2026_extraction.md` |
| the trees | the model's reaction lists | `scripts/generators/build_reaction_tree.py`, [REACTION_TREES.md](REACTION_TREES.md) |
| section 8: what the rules propose | the cited reaction rules and the species' structures | `data/lit/reaction_rules.yml`, `data/species/structures.yml`, `results/validation/network_hypotheses_prereg.md` |

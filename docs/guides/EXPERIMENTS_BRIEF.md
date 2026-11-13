# Would a few experiments make this model predictive?

*A decision brief, first written 2026-10-09 and updated 2026-11-12 (the update is the next section; text it
supersedes is marked). The model's purpose is to help alternative-protein scientists
build meat aroma into plant-based products: which precursors and which cooking programme give the meaty
thiols in a pea or soy matrix, without the off-notes and the acrylamide. Every number below comes from the
model's own records (paths given). The probabilities are judgement and are labelled as such. Protocols
are in [EXPERIMENTS.md](EXPERIMENTS.md).*

## The goal, where we are, and what changed since October (2026-11-12)

**The goal has not changed: predict the aroma a formulation and a cooking programme will give, in anyone's
kitchen.** A kitchen is a set of conditions (pH, the temperature the sample really saw, water, the precursors the
matrix brings), and chemistry maps conditions to products whoever does the cooking. The model cannot do this
today. This section says how far it is and what would move it.

**Four milestones, each with a test that exists in the repository.**

| | milestone | the test | today |
|---|---|---|---|
| 1 | Know the direction and rough size of each big lever, per odorant | a lever-table cell with three laboratories whose interval excludes zero ([lever_table.md](../../results/validation/lever_table.md)) | **4 of 75** cells: leaving out cysteine lowers 2-methyl-3-furanthiol about 6x; xylose gives about 3x the methional of glucose; raising the pH lowers 2-furfurylthiol and furfural (both marked fragile by the table) |
| 2 | The kinetic model beats that table on contrasts it never saw | the display rule of [STATUS_AND_PLAN](STATUS_AND_PLAN.md) sec. 0, step 5, on held-out laboratories | **cannot be run**: 55 cells now have held-out contrasts from two laboratories, but the model can be run on none of their 396 contrasts (no target for the odorant, a hydrolysate or food with no composition, or no printed pH; [step 5a](../../results/validation/step5a_display_check.md)). The model is paused until new data arrive ([record](../reference/lever_table/STEP5_RECORD.md)). On the contrasts that can be scored, the model loses to "no change" (median error 0.30 against 0.15 log units, 19 pairs; [contrast_scores.md](../../results/validation/contrast_scores.md)) |
| 3 | Ratios within threefold in a new laboratory after a few calibration pots | `maillard calibrate`, then `score` on that laboratory's other pots | one stand-in case (the Yiltirak ladder, below): 3 of 4 rows; not yet tried on a real partner |
| 4 | Absolute levels within threefold in a new plant matrix, no calibration | the out-of-sample panel | **10 of 23** rows outside the meaty thiols; absolute thiol levels are refused for a new laboratory (they differ between laboratories by about one log unit) |

**What the month between the two versions showed.**

- **The literature was keyed as a table of levers, and it answers almost nothing.** 70 sources, 1,643
  within-study contrasts, 31 laboratories, read by two blind readers: four cells of 75 are answered. The cause is the
  shape of the data, as this brief said in October: liquid model systems, one laboratory per question, peak
  areas, several things changed at once. **Four sources on disk have a heated plant- or fungal-protein contrast** (flaxseed, goji seed meal, two mushroom hydrolysates); none is a pea or soy isolate.
- **Laboratories disagree by more than most levers move things.** How much of that is measurement (isotope
  dilution against peak areas), how much is unrecorded conditions and how much is chemistry nobody has
  described is not known yet; that split is the next analysis, and it decides how far milestone 3 can go.
- **The model's pH error is in its structure, not its constants.** No choice of constants inside their bands
  repairs it ([decision_relevance_step0_addendum.md](../../results/validation/decision_relevance_step0_addendum.md)).
  Outputs that rest on the falsified pH balance are now refused.
- **A designed study was sized before asking anyone to run it.** 36 pots in pea isolate (below, experiment L).
  The check says what it will and will not buy: it halves the wrong-direction rate for a new kitchen where the
  literature has no plant-protein data (0.50 to about 0.31), it does not remove the spread between
  laboratories (a floor near 0.17 remains), and it is the only way milestone 2 can be tested at all
  ([lever_design_preposterior.md](../../results/validation/lever_design_preposterior.md)).

**What this changes in the advice below.** The October ranking stands, with one experiment added at the top
(L). The October sentence "the model will not become the design tool" described the evidence then and still
describes it now; it is a statement of where the model is, not a decision to stop wanting it. The path is the
four milestones, and each experiment below is marked with the milestone it moves.

## The short answer (2026-10-09)

**Would these experiments help design better plant-based products? Through their data, yes. Through the
model, only modestly.**

- **Building the product is an experimental job, and the method exists.** Sensomics (isotope-dilution
  quantification, odour activity, recombination and omission tests) defines the meat target. The same
  analysis on the plant product gives the gap, and a designed precursor study on the real product and
  process closes it. The beef target is largely published; the plant side is thin (see "The experimental
  route" below). This route does not need the kinetic model.
- **The model will not become the design tool.** After the experiments below, its realistic role is
  narrow: ranking precursor and cooking choices for ribose–cysteine reaction flavours in one plant
  protein, and explaining why a result changes when the process changes. The chance that it then changes
  a real formulation decision a product team would not reach faster by cooking and tasting: **~15–25 %**
  (judgement).
- **Two of the experiments are worth doing for the product, model or no model.** Odour thresholds in pea
  protein, and hot binding on pea protein, answer questions every plant-based flavour team has. Nobody
  has published either. The chance their results are directly useful to formulators: **~70 %**.
- **Stop scraping literature.** On 2026-10-09 four targeted searches against the 438 registered papers
  filled one gap out of about fifteen (`tasks/data_restructure_plan.md`, "The gap search"). The marginal
  paper no longer changes what the model can predict. *(2026-11-12: the lever table tested this from the
  other side. Keying 54 sources as contrasts answered one cell of 75. The one exception worth a search is a
  paper from a laboratory the corpus does not hold yet, because laboratories are what the cells lack;
  17 such papers from 13 new laboratories, keyed the same day with four supplements, answered three more cells, two of them weakly.)*

## Where the model stands

*Numbers as of 2026-11-12; the October figures follow in brackets. The counts fell because rows were refused,
not because predictions got better: absolute thiol levels for a new laboratory, outputs that rest on the
falsified pH balance, and three rows the source data cannot support.*

Out of sample, 10 of 23 scored rows land within threefold of the measurement, with a median miss of 4.4×
[11 of 44, median 9×, worst 9,630×]. Directional claims ("more X gives more Y") are right 10 times out of 19
[16 of 31], which is a coin flip (`results/validation/core_panel_scores.md`, `core_directional_scores.md`).

| lane | what it is for in a plant-based product | rows | within 3× | October |
|---|---|---|---|---|
| sulfur | the meaty thiols, the aroma you are trying to build | 0: absolute levels refused; ratios within one laboratory are given with an interval | – | 2 of 19, median miss 29× |
| lipid | the beany and green off-notes of pea and soy (hexanal and friends) | 7 | 4 | 4 of 8 |
| acrylamide | the safety constraint on hot processing | 10 | 4 | 4 of 12 |
| trunk | the sugar chemistry underneath everything | 6, all from one pot | 2 | 2 of 6 |

The lane that matters most for the product, the meaty thiols, is the worst: 29 of the 39 thiol rows removed
by the refusal were misses.

## Why more papers have not helped

The failures are not mainly a shortage of numbers. Five measured facts say so.

1. **Widening the uncertainty on every parameter barely helps.** With every prior uncapped, a nominal
   90 % interval still covers only 21 % of rows, up from 19 %. The model is wrong in a way no choice of
   rate constants fixes. The problem is missing or misshapen chemistry, which is to say structure, not
   parameters (`tasks/data_restructure_plan.md`, the prior-width sweep).
2. **The literature almost never measures what sets temperature dependence.** The meaty-thiol fit reads
   35 distinct systems, and **none of them appears at more than one temperature**. Most are pots held at
   145 °C for 20 minutes. Both of the lane's activation energies are therefore unidentified, and its
   held-out misses grow with distance from 145 °C (Spearman ρ = −0.82). An extruder, a pan and a retort
   all run at different temperatures, so a plant-based product formulated with this model is
   formulated by extrapolation.
3. **The plant matrix is the gap the literature cannot close.** Every plant-protein paper in the corpus
   computes odour activity with a threshold measured in **water**; no paired water-versus-pea threshold
   exists anywhere in reach. The only binding constants for plant protein are at 37 °C. On the one hot
   soy pot where it can be checked, the model moves hexanal the wrong way and by four orders of magnitude
   too little: heat and acid *release* bound aldehyde, and the model has no channel for that.
4. **Most of the rest is single endpoints, often relative.** Of 47 fitted coordinates, 27 are weakly
   identified or not identified at all (`kinetic_core_b38_identifiability.md`). There are 301 PDFs on
   disk but only four curated time series, all from one laboratory. Almost every constant comes from one
   group, and where a second group has been read, one constant was off by 466×.
5. **Some inputs are not chemistry at all.** On the fat lane the dominant uncertainty is the oil's
   peroxide value, which papers rarely report. One paper's "formed" volatiles were 36–42 % present
   before heating; correcting for the unheated blank turned two 30× misses into 2×.

Reading another hundred papers of the same kind adds rows of the same shape, and the identifiability audit
says it is the shape that fails: "it is the design that pins, not the number of rows."

## What the literature does not contain (searched 2026-10-09)

Four targeted searches against the ~440 papers already registered filled one gap out of about fifteen.
The gaps that decide the model's accuracy are empty; reading cannot close them. Each paper named here has
a dossier under `data/lit/extraction_dossiers/`, read by eye against its PDF.

| gap | verdict | closest evidence |
|---|---|---|
| Meaty-thiol formation at ≥ 2 temperatures, absolute, in a defined system | empty | one four-temperature emulsion series with one time point each (`yiltirak2026`) |
| Loss rate of free meaty thiols in buffer at ≥ 2 temperatures | empty | mechanism only: copper/iron thiol oxidation (`kreitman2016`); protein turning disulfides back into thiols at 100 °C (`mottram1996`) |
| Free cysteine against time in a heated sugar–cysteine pot | empty | copper rate law at 25–37 °C (`ehrenberg1989`) |
| Thiamine to meaty thiol at ≥ 2 temperatures | empty | total thiamine loss only, Ea 23–121 kJ/mol by matrix (`mauri1992`, `ramaswamy1990`, `brehm2020`) |
| Odour thresholds in plant protein vs water, same panel | effectively empty | one soymilk thesis with a different panel per test (Norton 2003, no dossier) |
| Flavour binding to plant protein at ≥ 50 °C in water | empty | 21 °C constants for five proteins (`snel2023`, `snel2023b`); heat pre-treatment read at 37 °C (`heng2005`) |
| Hexanal formation rate in wet plant protein | empty | whey at 50–90 °C, figure only (`jansson2020`) |
| Same-sample quantified comparison of beef and a plant analogue | empty | sniffing panels only (`thong2024`) |
| Designed precursor study inside an analogue, isotope-dilution measured | empty | dose screens with peak areas (`ma2026`, `milani2024`); patents, presence only (`us9943096`) |
| Furfural loss in water | filled | `almhofer2022` |

## Three testable hypotheses (judgement, not measurement)

| hypothesis | for | against | test | P(holds) |
|---|---|---|---|---|
| **Most cysteine is lost in the first minutes of heating**, so the model over-produces meaty thiols | one added fast cysteine sink takes meaty-thiol rows within 3× from 2 to 7 of 19 and weakens the temperature trend (ρ −0.80 → −0.38) | it needs ~90 % of cysteine gone in 5 min at 100 °C, pH 5.5; copper reaches that only with an unmeasured Ea of ~90–110 kJ/mol | free cysteine against time, 0–30 min at 100 °C, ± chelator: one day | ~0.3 for fast loss of any cause; ~0.1–0.15 that copper is the cause |
| **Thiamine pays off when braising, not searing** | direction matches two meat studies (`thomas2014`, `madruga1997`) | thiamine breakdown at meat pH has Ea ~115 kJ/mol against the engine's 64; a xylose + thiamine flavour at 100–120 °C shows no free meaty thiols (`us20250221436a1`) | experiment 1b below | ~0.3 |
| **Plant protein can protect meaty aroma**, not only trap it | egg albumin turned up to 45 % of a meaty disulfide back into its thiol at 100 °C (`mottram1996`) | no plant protein tested; no rates | fed disulfide with pea and soy protein at 100 °C, isotope dilution | open |

## The evidence that the right experiment does help

The repo has two cases where the missing data shape existed. Both worked.

- **Fed-intermediate time courses pin what they touch.** Feeding pure 3-deoxyglucosone at one temperature
  and following it over time let three new steps be added and fitted, every one to a tenth of a decade.
  The only fit on the branch that pinned every constant it touched was built this way.
- **A four-temperature ladder from one laboratory.** Using Yiltirak 2026 (ribose and cysteine in an
  oil-in-water emulsion) as a stand-in laboratory, fitting at 100 and 120 °C and predicting 110 and
  130 °C, cut the held-out median miss from **115× to 1.2×**, with 3 of 4 rows within threefold
  (`results/validation/calibration_prereg.md`, T5). Most of that gain was a per-lab response factor, and
  the 130 °C thiol still trended the wrong way. A calibration fixes a level but not the sign of a trend,
  and that limit is exactly what experiment 1 below probes.

## The experiments worth paying for, in order

Step-by-step protocol cards for experiments 0, 1, 1b, 2 and 3 are at the top of [EXPERIMENTS.md](EXPERIMENTS.md#protocol-cards).
Experiment L was added on 2026-11-12; its design is [FACTORIAL_LEVER_STUDY.md](FACTORIAL_LEVER_STUDY.md). The
milestone each experiment moves is in brackets. The order is judgement: L comes first because it is the
smallest, is useful with or without the model, and is the only one that can test the model on new laboratories.

| # | experiment | size | what it buys | why it matters for plant-based products |
|---|---|---|---|---|
| L | **The big levers in pea protein** ([FACTORIAL_LEVER_STUDY.md](FACTORIAL_LEVER_STUDY.md)): sugar (ribose or glucose) × sulfur source (cysteine or thiamine) × pH (5.5 or 6.5) in 15 % pea isolate, 140 °C for 20 min, isotope-dilution GC, in two laboratories from one shipped lot (one laboratory is acceptable) | 36 pots in all | The first heated plant-protein contrasts for the sugar, sulfur-source and pH cells of every odorant family, and the held-out laboratories the model's first real test needs [milestones 1 and 2]. It will not remove the spread between laboratories; record the conditions listed in the next section so that spread can be explained. | These are the three choices every reaction-flavour formulation makes first, and today each is answered only in buffer or in meat. |
| 0 | **Residual cysteine, one day**: free cysteine against time in a ribose–cysteine pot, 0–30 min at 100 °C, pH 5–5.5, ± chelator in excess | 42 samples (two arms), one day | Whether fast early cysteine loss explains the meaty-thiol over-prediction (first hypothesis above) | Decides the model's main structural question before anyone pays for experiment 1 |
| 1 | **Where the thiols go** ([EXPERIMENTS.md §1](EXPERIMENTS.md)): ribose + cysteine at pH 5, 100 °C and 140 °C, against time, with fed-thiol, nitrogen, chelator and diketone arms; isotope-dilution GC | ~130 vials, ~2 weeks of GC | Settles where the 29× over-prediction of meaty thiols comes from: formation, the thiol sink, or precursor loss. A copper-catalysed loss of cysteine is one candidate, now unlikely at 33 mM cysteine (Ehrenberg 1989). Gives the lane its first activation energies. | Ribose + cysteine is the standard reaction flavour behind plant-based meat aroma. Without this, the model cannot say which cooking programme gets the most meat aroma out of it. |
| 1b | **Thiamine arm in the same pots** ([pre-registration §10](../../results/validation/cultivated_tissue_invariance_prereg.md)): thiamine vs cysteine at both temperatures | +24 vials | Tests the one prediction a flavour chemist would not make unaided: thiamine pays off at braising temperatures and hardly at all when searing. | Thiamine and yeast extract are cheap, common reaction-flavour ingredients. If the claim holds, it is a direct formulation rule; if not, it is one fewer false lead. |
| 2 | **Binding on pea protein, hot** ([EXPERIMENTS.md §4](EXPERIMENTS.md)): hexanal and an alkenal at 40, 70 and 90 °C, pH 4.5 and 7, headspace *and* total on the same aliquot | one design, not costed in the guide | The reversible protein binding the model lacks, and the temperature dependence of hexanal formation in a wet protein. | Off-note control is half the problem in pea and soy, and this is where the model is out by four orders of magnitude and has the sign wrong. |
| 3 | **Odour thresholds in pea protein** ([EXPERIMENTS.md §3](EXPERIMENTS.md)): 6–8 compounds, including a meaty thiol and hexanal, water vs 3 % pea dispersion on the same panel | a sensory panel, not a GC | Lifts the refusal of every matrix-corrected odour activity. | Whether a predicted concentration is smelled at all in the product. Only a panel can measure it. |
| 4 | **Fed-intermediate ladder** ([EXPERIMENTS.md §2](EXPERIMENTS.md)) | one design, several pots | Pins seven sliding constants and returns eight refused rows. | Precision after correctness: do it only if experiment 1 says the structure is right. |
| — | **No laboratory needed**: ask authors for raw data, and add unheated blanks to any partner's next run | email | Seven refused rows become answerable. | Nearly free. |

Experiments 1 and 1b share one bench and one extraction method. Experiment 2 needs headspace and total
extraction on a protein suspension, and experiment 3 needs a sensory panel; these may be different
partners. The Reading group behind the Yiltirak ladder already heats ribose–cysteine systems across
temperatures and quantifies these thiols, which makes it the natural first call for 1 and 1b.

## What to record, so that the data add up

Most of the 1,643 keyed contrasts could not be combined, and each item below is a reason found while keying
them. A study that records these is worth several that do not, to this model and to anyone else's.

| record | why (what its absence cost here) |
|---|---|
| **Absolute amounts**, by isotope dilution or a stated calibration; never peak areas alone | peak areas from two laboratories cannot be compared, and a response factor per compound is then one more unknown per laboratory |
| **An unheated blank** of every matrix | in one paper 36–42 % of the "formed" volatiles were there before heating; correcting for it turned two 30× misses into 2× |
| **pH before heating and after cooling**, in every vessel, and **the buffer by name and concentration** (or "none") | pH drifts during heating wherever the pot is not strongly buffered (Zhou 2023's pots ended near pH 3.3 from an initial 6 or 7), and the model's largest structural error is its pH response; the two thiol x pH conflicts between laboratories in the lever table dissolve once an unbuffered initial pH is read as what it is |
| **What the reference arm holds**: the full base mixture of the "none" or control pot | in 9 of the 16 widest disagreements in size, the "none" arm was a different pot in each study (glycine in one, no amine in another, 2 g of protein in a third) |
| **Every dose, including zero, in one unit**, with the sugar's amount beside it | cysteine both feeds and starves a pathway depending on its ratio to the sugar; two conflicts between laboratories were opposite halves of rise-then-fall dose curves |
| **The temperature the sample saw** (a logged thermocouple in a dummy vessel), not the oven setting | a 1 g sample in a 150 °C oven for 10 min is not at 150 °C for 10 min, and at a typical activation energy 30 °C is roughly a tenfold change in rate |
| **Water before and after** (mass loss, or water activity) and whether the vessel was sealed | dry and wet arms move 2-furfurylthiol in opposite directions for pentoses and hexoses (Hofmann 1998) |
| **What the matrix brings**: free sugars, free amino acids, and the lot; for a hydrolysate or food, its composition | a flour charges precursors nobody added; a laboratory difference must not be a lot difference; 148 of the 396 held-out contrasts of step 5a could not be put to the engine because their hydrolysate's or food's composition is not printed |
| **One lever at a time, or a full factorial** | 183 keyed contrasts changed several things at once and feed no cell |
| **A detection limit wherever a compound is "not detected"** | the thiamine cells are unusable because the no-thiamine arm is "not detected" with no limit |
| **Vessel volume and fill** | headspace sets how much of a volatile thiol is in the liquid that is extracted |
| **The numbers as a table**, in the supplement | values read off a figure carry a reading error that two readers here resolved by hand |

A few **calibration pots** make a new laboratory comparable with the rest: a handful of defined pots that
every laboratory runs on its own method (which pots is not fixed yet; experiment L's two check pots, the matrix with nothing added and the
buffer without the matrix, are the first version). `maillard calibrate` turns a laboratory's own measurements into its
response factors [milestone 3].

## What to expect (judgement, not measurement)

- Experiment 1 settles the meaty-thiol structure one way or the other: **~90 %**. The design separates
  the hypotheses whatever the outcome.
- Meaty-thiol predictions then land within threefold in buffered pots like those tested, between 100 and
  140 °C: **~60 %**.
- Those predictions carry over to a pea or soy product within threefold *without* experiments 2 and 3:
  **~25 %**. *With* them: **~40 %**. The matrix is the larger unknown, and nobody has measured it.
- More than half of the whole held-out panel lands within threefold: **~20–25 %**. Fat-lane inputs and
  the trunk's single pot remain.
- The thiamine claim survives experiment 1b: **~30 %**.

So the realistic product is **a model that is right, and knows it is right, for ribose–cysteine reaction
flavours over the cooking range, with a measured correction for one plant protein**, plus a clear yes or
no on its central structural question. It would not be a general Maillard predictor. It would, though,
be enough to rank precursor and process choices for meaty aroma in a pea-based product, which is the
decision the model exists to help with.

## The experimental route, and how far the literature has taken it

| step | what it does | where the literature stands (dossiers on file) |
|---|---|---|
| 1. Define the meat target | isotope-dilution quantification, odour activity, recombination | largely done for stewed and boiled beef: the meaty thiols, methional, furaneol, fat-derived aldehydes and, for stewed beef, the tallowy 12-methyltridecanal (`christlbauer2011`; Kerscher & Grosch 1998). No recombinate of grilled or roasted beef has been confirmed by omission tests |
| 2. Measure the plant product | the same analysis on the cooked analogue | pea isolate (Utz 2022) and one soy extrudate (`wang2026c`: 41-odorant recombinate, **no sulfur odorant found**, hexanal the key off-note). Commercial burgers vs beef only by sniffing (`thong2024`: meatiness tracks the sulfur set) |
| 3. The gap | which target odorants are missing, which off-notes are in excess | never done on the same samples with quantification: **the first missing study** |
| 4. Matrix test | spike the missing odorants into the cooked product | not done; experiments 2 and 3 above are its lab version |
| 5. Precursor design | designed formulation on the real product and process, measured on the key odorants and by a panel | dose screens only (`ma2026`: xylose:cysteine ~3:2 at 0.8 % in extruded soy; `milani2024`: thiamine raises meat odour), peak areas or sensory alone. **No response-surface study inside an analogue measured by isotope dilution** |

Several odorants the targets need are ones the model cannot produce at all: 12-methyltridecanal, dimethyl
trisulfide, most of the fat-derived 2-alkenals and dienals, 1-octen-3-one, sotolon and most pyrazines.
That is a second reason the model is a helper on step 5, not the route.

## When to stop

Write the stopping rule down before the data arrive. **If experiment 1 shows the structure is right (the
model with the missing step added) and the meaty-thiol miss on the new pots is still above about 5×, stop
extending the engine,** and publish it as a reasoned map of what is and is not known, which it already
is. In every case: stop scraping literature for this model. The evidence above says the marginal paper
no longer changes what it can predict.

*2026-11-12.* The rule stands as written, and it is a rule about effort, not about the goal. "Stop extending
the engine" means no further fitting on the literature on disk, which the lever table has since shown cannot
decide these questions. Work on the engine resumes when designed data arrive (experiment L, experiment 1, a
partner's time courses), judged by milestones 2 and 3 above.

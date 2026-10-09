# Would a few experiments make this model predictive?

*A two-page decision brief, 2026-10-09. The model's purpose is to help alternative-protein scientists
build meat aroma into plant-based products: which precursors and which cooking programme give the meaty
thiols in a pea or soy matrix, without the off-notes and the acrylamide. Every number below comes from the
model's own records (paths given). The probabilities are judgement and are labelled as such. Protocols
are in [EXPERIMENTS.md](EXPERIMENTS.md).*

## The short answer

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
  paper no longer changes what the model can predict.

## Where the model stands

Out of sample, 11 of 44 scored rows land within threefold of the measurement. The median miss is 9×
(about one decade) and the worst is 9,630×. Directional claims ("more X gives more Y") are right 16 times
out of 31 once pH and water activity are set aside, which is close to a coin flip
(`results/validation/core_panel_scores.md`, `core_directional_scores.md`).

| lane | what it is for in a plant-based product | rows | within 3× | median miss |
|---|---|---|---|---|
| sulfur | the meaty thiols, the aroma you are trying to build | 19 | 2 | 29× |
| lipid | the beany and green off-notes of pea and soy (hexanal and friends) | 8 | 4 | 2.8× |
| acrylamide | the safety constraint on hot processing | 12 | 4 | 7.7× |
| trunk | the sugar chemistry underneath everything | 6, all from one pot | 2 | 8.0× |

The lane that matters most for the product, the meaty thiols, is the worst.

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

| # | experiment | size | what it buys | why it matters for plant-based products |
|---|---|---|---|---|
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

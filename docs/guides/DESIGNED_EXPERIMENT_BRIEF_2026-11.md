# A designed experiment for a partner laboratory: brief (2026-11)

*Draft for Pablo to review. Nothing has been sent. Written for a food-chemistry PI (Wageningen FQD, or TUM /
Leibniz-LSB). It merges protocol experiments E3 and E1 of
[DESIGNED_TIMECOURSE_PROTOCOL_DRAFT.md](DESIGNED_TIMECOURSE_PROTOCOL_DRAFT.md) with the rows of
[ENGINE_VS_LITERATURE.md](../reference/pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3.1 that can share their pots.
Numbering: "E1 - E20" below are the model's mechanism errors (atlas sec. 3); the experiments are Parts A and B.*

## Why we need new data, not more papers

We build an open kinetic model of the Maillard reaction ([INTRODUCTION.md](INTRODUCTION.md)). A review of the
literature against the model found 20 mechanism errors (E1 - E20). Fixing them is the easy part. Testing the
fixes is hard:

- **The literature that showed an error cannot also test its fix.** So we froze a split before scoring anything
  ([contrast_benchmark/README.md](../reference/contrast_benchmark/README.md)): a laboratory whose papers were used
  to find an error is "development" for it; only other laboratories can test it. Even under the more permissive
  version of this rule, **8 of the 20 errors have fewer than two independent held-out laboratories**: E1, E2, E4,
  E8, E10, E15, E16, E17. No published dataset can test their fixes.
- **Some regimes were never measured** (atlas sec. 1.4): the pH *during* heating (almost every published pH is the
  starting pH); air vs inert headspace at 100 °C and above (no reaction family has it); the pentose
  deoxyosone split measured with an in-situ trap; 1:1 labelled pentose at 100 - 140 °C.

A new laboratory, a pre-registered design, and model predictions frozen and sent to you **before** the data
exist give clean test data. The design below covers five of the eight untestable errors (E1, E2, E8, E10, E17),
four thinly tested ones (E5, E6, E11, E12), and all four unmeasured regimes.

## The design: two pot families, about 105 - 110 vials

Sealed glass pressure vials, one vial per time point, the same stated fill and headspace in every vial, sample
temperature logged, ice quench. **pH measured after cooling in every vial** (and in situ if you can). Air
headspace unless stated.

**Part A: the meaty-thiol pot, with a labelled MFT spike** (protocol E3). Core: D-ribose 30 mM + L-cysteine 30 mM in
0.2 M sodium phosphate, pH 5.5. Times: 100 °C 0 / 15 / 30 / 60 / 120 / 240 min; 120 °C 0 / 5 / 10 / 20 / 40 /
80; 140 °C 0 / 2 / 5 / 10 / 20 / 40 (subsets: early, middle, late).

| arm | what changes | vials | what it decides |
|---|---|---|---|
| A1 | 1:1 [13C5]- / [12C5]-ribose (CAMOLA), 100 / 120 / 140 °C, 5 times | 15 | which carbon routes make MFT over time and how each speeds up with temperature (E11) |
| A2 | unlabelled ribose, [2H3]-MFT spiked at t = 0; 100 and 140 °C, air vs argon, 3 times | 12 | how fast MFT is lost in a working pot, and whether oxygen drives the loss (E12) |
| A3 | as A2 (air), spike added after cooling, before work-up; 100 / 140 °C, 2 times | 4 | loss in the pot vs loss in the work-up |
| A4 | as A2 (air, 140 °C), cysteine 0 or 100 mM, 3 times | 6 | whether cysteine protects MFT (E12) |
| A5 | as A2 (air, 140 °C), pH 4.5 or 6.5, 3 times | 6 | the pH law of furfural / FFT supply; the model's is 10 - 70x too steep (E8) |
| A6 | fully [13C5]-ribose + unlabelled norfuraneol (30 mM each) + Cys; 120 / 140 °C, 2 times | 4 | how much MFT comes via free norfuraneol at cooking temperature (E6) |
| A7 | 1:1 [13C6]-glucose 30 mM + Cys 30 mM, 120 °C, pH 5.5 / 6.5, 3 times | 6 | whether glucose makes MFT through an intact C5 chain; the model has no such route (E1, its top-ranked error) |

**Part B: the pentose hub ladder, in-situ vs ex-situ OPD** (protocol E1). Core: D-xylose 50 mM + glycine 50 mM, 100
°C, 0.2 M phosphate; times 15 / 45 / 120 / 240 min. "In-situ": o-phenylenediamine (OPD, in excess) added before
heating, so it traps every dicarbonyl as it forms (a flux). "Ex-situ": OPD added after cooling, as usual (a
standing level). The difference tests whether published deoxyosone levels are partly a trapping artefact.

| arm | what changes | vials | what it decides |
|---|---|---|---|
| B1 | pH 5.5 / 6.5 / 7.5, in-situ and ex-situ OPD, 4 times | 24 | the split of the pentose into 1-DX vs 3-DX and fragments, and its pH law (E2); amine return from the Amadori compound (E17) |
| B2 | pH 6.5 in water (NaOH, no buffer), 0.05 M phosphate, or 0.2 M citrate; ex-situ, 3 times | 9 | whether phosphate catalyses the Amadori steps, and how it saturates (E5) |
| B3 | pH 6.5 + cysteine 5 or 25 mM; ex-situ 3 times each, in-situ at 25 mM, 3 times | 9 | whether cysteine traps 3-DX, and how much |
| B4 | no glycine: xylose 10 + cysteine 10 mM, pH 4.5 / 6.5 in water and pH 6.5 in phosphate, 3 times | 9 | dilute cysteine fate: TTCA formation and release, and its pH law (E10) |

Plus 3 - 5 replicate vials at the centre points (A1 120 °C, B1 pH 6.5 both OPD modes). **If the budget is tight
(~55 vials):** A1, A2, B1 and B3's in-situ arm. Phosphate barely buffers at pH 4.5 - 5.5, so those pots will
drift; that is why the after-cooling pH matters (the model runs on the measured pH). Optional extras: 160 °C in A1;
a fed xylose-glycine Amadori compound in B1 (+6 vials), if you can make it.

Not covered, because they need a different pot: the Strecker / pyrazine errors (E3, E4, E13 - E15), lipid (E7),
acrylamide (E9) and norfuraneol decay (E16). A fed-dicarbonyl pot could be a second round.

## Analytes and methods (suggested; your methods take precedence)

| analyte | where | method |
|---|---|---|
| MFT, FFT, bis(2-methyl-3-furyl) disulfide, 2-mercapto-3-pentanone | A | SIDA GC-MS |
| MFT isotopologues M+0 to M+5 | A1, A6, A7 | GC-MS (CAMOLA) |
| furfural, norfuraneol, furfuryl alcohol | A, B | HPLC-DAD or GC-MS |
| quinoxalines of 1-DX, 3-DX, glyoxal, methylglyoxal | A5, B | OPD, HPLC-DAD (~315 nm) or LC-MS, with standards |
| ribose / xylose / glucose | A, B | HPAEC-PAD |
| free glycine, the Xyl-Gly Amadori compound | B | LC-MS/MS or amino-acid analyser |
| free cysteine, TTCA, the Cys Amadori compound | A, B3, B4 | LC-MS/MS with synthesised standards |
| H2S (optional) | A2 | 34S-SIDA after in-vial derivatisation (Li, Schieberle & Steinhaus 2017) |
| pH after cooling | every vial | glass electrode, 25 °C |

**One point to agree first.** The [2H3]-MFT spike (M+3) and the C2 + C3 isotopologues of a 1:1 labelled pot (M+2,
M+3) overlap on a unit-resolution GC-MS. Hence spike and CAMOLA run in separate vials (A1 vs A2), unless you have
high-resolution MS (the 13C3 / 2H3 gap is 8.8 mDa). The spike, the work-up standard and the CAMOLA isotopologues
must all carry distinguishable labels.

## What the result decides

- **Part A alone** fixes the most likely missing thiol element (the MFT loss law) and the route mix over time.
  Our judgement: P(it changes a structural choice in the model) ~0.65.
- **Part B** settles the most upstream open question (the pentose split as a rate, not a level) and tests the
  OPD method; P ~0.6.
- Together: P(at least two of the hub / thiol questions resolved) ~0.7, against ~0.3 from the literature alone.

Every vial's prediction is frozen and sent to you before the data arrive. Each error's pass / fail rule is written
down beforehand. A fix counts only if it beats a typical-yield baseline on your data, out of sample. Negative
results count equally.

## What the partner gets

- **Co-authorship** on the data paper and on the model-validation paper (the blind test).
- **The model**: open source ([github.com/PabloAMC/Maillard](https://github.com/PabloAMC/Maillard)); we run it for
  your own recipes and conditions on request.
- **The pathway atlas**: which Maillard reactions happen under which conditions, cell by cell, with laboratory and
  evidence type (`docs/reference/pathway_atlas/`).
- **The analysis**: we write the pre-registration and do the kinetic modelling. You keep your data and decide
  when to publish.

Open points: which arms suit your analytical set-up (Part B fits multiresponse kinetics with HPLC; Part A needs
thiol SIDA and CAMOLA, so the parts could go to two groups); the labelled standards ([13C5]-ribose, [13C6]-glucose,
[2H3]-MFT, TTCA, the Amadori compounds) and their cost; air vs argon handling.

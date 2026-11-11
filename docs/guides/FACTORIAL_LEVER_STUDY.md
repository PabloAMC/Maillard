# One factorial lever study in a plant-protein matrix, run in one or two laboratories (plan P3, STATUS_AND_PLAN sec. 4.1 and sec. 0)

*Draft, 2026-11-11. Design and sizing only; nothing has been sent. Sending a request to any laboratory is the
owner's action (STATUS_AND_PLAN sec. 0, owner row). Written for food-chemistry groups with SIDA GC-MS; redesigned the same day after the preposterior
and predictive checks.*

## Why this study, and why now

The lever table (`results/validation/lever_table.md`) keyed the literature's within-study contrasts: 54 sources,
946 usable contrasts, 20 laboratories. It answers **one** formulator question at the three-laboratory standard
(leaving out cysteine lowers MFT, about 6x). Six cells read "direction unclear"; every other cell has fewer than
three laboratories. The experiment map (`results/validation/experiment_map.md`) shows two binding gaps, and they are
the same in every odorant family:

1. **Laboratories, not replicates.** Most cells have one or two laboratories; a third is what turns a listed contrast
   into an answer.
2. **Matrix.** No keyed cell has a heated plant-protein contrast. Every published lever is measured in buffer, meat or
   a process flavour.

One designed study in a pea-isolate matrix gives the first plant-protein data to many cells at once, in the matrix
users cook in. It is also the only route by which the kinetic
engine can earn a cell: step 5 needs held-out contrasts from two new laboratories, and no thiol x pH, sugar or
temperature cell has two on disk. One laboratory alone would not open step 5 (preposterior check below), though for the prediction itself one
laboratory does about as well as two or three (predictive check below).

## The design: one 2^3 factorial, 36 pots, in two laboratories (or one)

*Redesigned 2026-11-11 after the preposterior check below. The 2^4 single-laboratory draft is kept in git history.
What the study must buy is laboratories, not replicates.*

**Base pot.** 15 % w/w commercial pea protein isolate in water, 1 kg/L, sealed glass vials, one vial per pot, the
same fill and headspace in every vial. One product and one lot for every laboratory, shipped from one batch (state
product, lot and protein %), so that a laboratory difference is not a lot difference. Heat 140 C, 20 min (sample
temperature logged), ice quench. pH measured before heating and after cooling in every vial.

| factor | low | high | the cells it feeds |
|---|---|---|---|
| A. sugar | glucose 30 mM | ribose 30 mM | sugar identity x every family |
| B. sulfur source | cysteine 10 mM | thiamine 10 mM | sulfur source (thiamine connects to cysteine, which no keyed study does) |
| C. pH (set before heating) | 5.5 | 6.5 | pH x every family, 6 - 7.5 band (midpoint 6.0) |

Yeast extract is dropped: the lever table has no yeast-extract class, so it would feed no cell (preposterior check,
point 4). It returns only if the owner wants it and a PREREG amendment adds the class before any data are read.

**Per laboratory**, randomised run order:

| | two laboratories (default) | one laboratory (acceptable) | three laboratories (not worth the coordination) |
|---|---|---|---|
| factorial pots | 8 cells x 2 pots | 8 cells x 4 pots | 8 cells x 1 pot |
| matrix background (base pot, no sugar, no sulfur source) | 1 | 2 | 2 |
| buffer only (A-low / B-low / C-low without pea isolate: the matrix contrast) | 1 | 2 | 2 |
| pots per laboratory | 18 | 36 | 12 |
| total | 36 | 36 | 36 |

The number of laboratories is chosen by the predictive check below, not by the 'answered' count: for the
prediction a new kitchen gets, one, two and three laboratories are nearly equal. Two is the default because one
laboratory cannot show any between-laboratory spread (its intervals are the worst calibrated) and step 5 needs two.

**Analytes** (SIDA where a labelled standard exists; otherwise one internal standard and the same method in every
pot, so within-study ratios cancel the response factor; each laboratory may use its own method, which is part of
what a new laboratory measures):
- meaty thiols: MFT, FFT, bis(2-methyl-3-furyl) disulfide, H2S (if the group can trap it);
- Strecker aldehydes and methional: 2- and 3-methylbutanal, phenylacetaldehyde, methional;
- pyrazines: 2,5-dimethylpyrazine, 2,3,5-trimethylpyrazine, 2-ethyl-3,5-dimethylpyrazine, total pyrazines;
- off-notes: hexanal, nonanal, 2-pentylfuran, 1-octen-3-ol (the isolate's lipid background);
- furanones: DMHF, norfuraneol, furfural.
- Optional, if available: GC-O / AEDA on the 4 corner pots that differ in A and B (the perception layer).

**Keying (PREREG amendment A4).** Each laboratory gives ONE contrast per lever and odorant: the mean log10 value of
the factorial pots at the lever's high level against those at its low level (4 against 4, or 8 against 8 with two pots
per cell). Interactions are estimated within each laboratory and reported with the
study; they do not enter the table.

## Sizing

Per-pot reproducibility of a log value is taken as 0.10 - 0.15 dex (PREREG sec. 5 uses 20 % CV, n = 2). One
laboratory's main effect is a difference of two means of 4 pots (three-laboratory arm) or 8 pots (two-laboratory
arm): standard error 0.07 - 0.11 dex or 0.05 - 0.075 dex. A 2x effect (0.30 dex) is therefore seen by each
laboratory on its own, and the measurement is not what limits the answer. What limits it is how much the effect
varies between laboratories, and only more laboratories measure that. With one pot per cell, a laboratory's own
error estimate rests on the duplicated checks and the three-factor interaction; the table does not need it (sec. 5
uses the declared measurement variance).

## What each result does to the tool

From the preposterior check (matrix shift 0.3 dex, A4 keying; means of 200 simulated outcomes):

| | today | two laboratories | three laboratories |
|---|---|---|---|
| lever-table cells newly answered | - | 11 (8 with a literature estimate today) | 36 (8 with a literature estimate; the rest depend on the prior, read as optimistic) |
| experiment-map family x lever cells newly answered | 1 answered | 7 | 18 |
| cells where step 5 can test the engine (two TEST laboratories) | 0 | 57 | 57 |
| plant-protein matrix in the keyed record | none | every cell the study touches | every cell the study touches |

Every prediction (the lever table's and the engine's) is frozen and sent before the data exist, with its pass rule.
A cell is "answered" only by the PREREG sec. 5 rule.

## Preposterior check: does the study unblock anything? (2026-11-11)

Before anything is sent, `scripts/generators/build_lever_design_preposterior.py` simulated each candidate design 200
times. Each run draws true effects from the lever table's own posterior predictive (cells with no keyed estimate draw
from the PREREG prior), adds a laboratory offset, a 0.1 dex within-laboratory residual and SIDA / GC-MS measurement
error, and optionally a pea-matrix shift of 0.3 dex shared by all new laboratories. It then keys the synthetic
contrasts and refits with the table's own model and reading rule. Full table:
[lever_design_preposterior.md](../../results/validation/lever_design_preposterior.md). Matrix shift 0.3 dex,
main-effect keying:

| design (36 pots unless stated) | cells newly answered, mean [90 %] | of which with a literature estimate today | experiment-map cells newly answered | cells with 2 TEST laboratories (step 5) |
|---|---|---|---|---|
| the draft above: one laboratory, 2^4 x 2 | 4.1 [3, 5] | 4.1 | 3.1 | 0 -> 2 |
| two laboratories, 2^3 x 2 each, no yeast extract | 11.2 [8, 14] | 7.7 | 6.7 | 0 -> 57 |
| three laboratories, 2^3 x 1 each, no yeast extract | 35.8 [30, 41] | 8.2 | 18.5 | 0 -> 57 |
| three laboratories, sugar x sulfur source 2 x 3, x 2 | 30.7 [25, 36] | 6.0 | 8.9 | 0 -> 57 |
| reference: the draft in two laboratories (72 pots) | 11.2 [8, 14] | 7.8 | 6.8 | 0 -> 57 |

What it shows:

1. **As drafted, the study barely unblocks anything.** One laboratory can only finish cells that already have two;
   there are about four. It gives step 5 (the engine against the lever table on held-out contrasts from two or more
   TEST laboratories) two cells. Its strength (power for 0.15 dex effects and interactions) is the wrong currency: the
   binding constraint is laboratories, and replicates inside one laboratory do not buy them.
2. **The same pots in two laboratories do about 3x better on every count, and give step 5 its 57 cells.** Running the
   whole draft in a second laboratory (72 pots) adds nothing over two laboratories at 36: the yeast-extract factor and
   the second replicate are spent on cells the table cannot read.
3. **A third laboratory mostly answers cells no paper has keyed** (34 against 11; the literature-backed count stays
   about 8). Those cells would be answered by the three new laboratories alone, and how many depends on the prior
   for effect sizes (mu ~ N(0, 0.5) dex). Read the 34 as optimistic; the 8 and the step-5 count barely depend on it.
4. **Yeast extract feeds no cell.** PREREG sec. 4 has no yeast-extract lever class. Keep it only after an amendment
   adds one (categorical, reference "none") before any data are read.
5. **Keying a factorial pair by pair would overstate the answer** (per-pair keying: up to +3 cells). The categorical
   network has one heterogeneity scale per pair and no laboratory random effect, so a laboratory's 2^(k-1) contrasts
   for one lever would count as independent draws. Key a designed study as one contrast per laboratory per lever (the
   mean over the other factors), and state that in the pre-registration. The one cell the literature answers today
   (MFT, no sulfur source vs cysteine; 7 contrasts, 4 from one laboratory) survives this: refitted with one contrast
   per laboratory, it stays answered in all 8 subsamples (estimate -0.81 to -1.21 dex).

**Redesign (done the same day, owner-approved):** the design section above is now the 2^3 (sugar x sulfur source x pH)
with checks, keyed by PREREG amendment A4. The laboratory count was then revised by the predictive check below: two by
default, one acceptable, three not worth it.

## Predictive check: does the number of laboratories matter? (2026-11-11)

The 'answered' count above is partly circular: it counts how often new laboratories clear our own three-laboratory
rule. The non-circular question is how much better the table predicts **a new kitchen cooking in pea isolate**. The
same simulated worlds were scored for the literature alone and for the same 36 pots in one (2^3 x 4), two (2^3 x 2) or
three (2^3 x 1) laboratories ([the last section of the record](../../results/validation/lever_design_preposterior.md);
200 worlds per setting). The 57 cells the 2^3 touches, matrix shift 0.3 dex:

| | literature alone | one laboratory | two laboratories | three laboratories | floor: the pea mean known exactly |
|---|---|---|---|---|---|
| sign wrong for a new kitchen, all cells | 0.44 | 0.30 | 0.29 | 0.29 | 0.16 |
| the same, 12 cells with a literature estimate | 0.22 | 0.21 | 0.22 | 0.22 | 0.14 |
| the same, 45 cells with none | 0.50 | 0.33 | 0.31 | 0.30 | 0.17 |
| 90 % interval coverage, cells with none (nominal 0.90) | 0.88 | 0.57 | 0.62 | 0.64 | - |
| confident calls (P(same sign) >= 0.8) that are wrong, cells with none | - | 0.30 | 0.26 | 0.23 | - |

With the design laboratories sharing our protocol and lot (offsets half the literature's), sign errors are 0.28 for
all three designs. Coverage then FALLS with more laboratories (0.61, 0.53, 0.49): laboratories that agree with each
other make the fit understate the spread a new kitchen has.

What it shows:

1. **The study's value is pea-matrix data where the literature has none.** It cuts the sign error on those 45 cells
   from a coin flip (0.50) to about 0.31. On the 12 cells the literature already estimates, it changes almost nothing
   (0.22 -> 0.21).
2. **How the 36 pots are split across laboratories barely matters for prediction** (0.30 / 0.29 / 0.29). A third
   laboratory does not earn its coordination cost.
3. **The binding problem is calibration, not data.** With 1 - 3 laboratories in a cell, the table's 90 % new-laboratory
   interval covers about 60 %, and a third of 'confident' signs are wrong, against the 20 % the 0.8 threshold promises.
   Even the literature-only cells cover 0.83. Whether the fit or the simulator is at fault is being checked; until it
   is, no cell this study adds should be read as answered on its interval alone.
4. **The floor is high.** Even a predictor that knew the pea-matrix mean exactly would get the sign wrong for 16 % of
   new kitchens, because laboratories scatter that much. No study design removes that; only conditioning on what
   differs between kitchens would.

## What the partner gets

Co-authorship on the data and on the validation; the open tool and the lever table run on their own recipes; the
analysis and the pre-registration done by us. They keep their data and decide when to publish.

## Open points for the owner

- Which one or two groups. Each must be a separate laboratory (not two groups sharing one instrument or one method
  lead). Preferably each is absent from `docs/reference/contrast_benchmark/split.yml`: PREREG sec. 2 makes an absent
  laboratory TEST, while a laboratory already in the split keeps its partition, and a DEVELOPMENT one gives step 5
  nothing. The preposterior counts assume all new laboratories are TEST.
- Whether yeast extract matters enough to bring back (it needs a PREREG amendment adding the class, and then it
  needs three laboratories of its own).
- Whether to add the pan-fried arm (180 C surface) at the two A x B corners: the cooking-method cell has no contrast.
- Whether the blind-test framing already used once (predictions frozen and sent before the data) applies here too.

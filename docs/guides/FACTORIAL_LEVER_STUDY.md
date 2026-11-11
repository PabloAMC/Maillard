# One factorial lever study in a plant-protein matrix (plan P3, STATUS_AND_PLAN sec. 4.1 and sec. 0)

*Draft, 2026-11-11. Design and sizing only; nothing has been sent. Sending a request to any laboratory is the
owner's action (STATUS_AND_PLAN sec. 0, owner row). Written for a food-chemistry group with SIDA GC-MS.*

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

A single designed study in a pea-isolate matrix adds one laboratory to many cells at once, in the matrix users cook
in. It is also the only route by which the kinetic engine can earn a cell (step 5 needs held-out contrasts from new
laboratories, and no thiol x pH, sugar or temperature cell has two on disk).

## The design: 2^4 full factorial, two replicate pots, plus four checks (36 pots)

**Base pot.** 15 % w/w commercial pea protein isolate (state the product, lot and its protein %) in water, 1 kg/L,
sealed glass vials, one vial per pot, the same fill and headspace in every vial. Heat 140 C, 20 min (sample
temperature logged), ice quench. pH measured before heating and after cooling in every vial.

| factor | low | high | the cells it feeds |
|---|---|---|---|
| A. sugar | glucose 30 mM | ribose 30 mM | sugar identity x every family |
| B. sulfur source | cysteine 10 mM | thiamine 10 mM | sulfur source (thiamine connects to cysteine, which no keyed study does) |
| C. yeast extract | 0 % | 1 % w/w (state product) | yeast extract x every family (no contrast on disk) |
| D. pH (set before heating) | 5.5 | 6.5 | pH x every family, in the 4.5-6 and 6-7.5 bands |

16 cells x 2 replicate pots = **32 pots**, randomised run order, blocked by day if needed. **Four checks:** base pot
with no sugar and no sulfur source (the matrix's own background), duplicated; the A-low/B-low/C-low/D-low pot
without pea isolate (buffer only), duplicated: the matrix contrast the map lacks.

**Analytes** (SIDA where a labelled standard exists; otherwise one internal standard and the same method in every
pot, so within-study ratios cancel the response factor):
- meaty thiols: MFT, FFT, bis(2-methyl-3-furyl) disulfide, H2S (if the group can trap it);
- Strecker aldehydes and methional: 2- and 3-methylbutanal, phenylacetaldehyde, methional;
- pyrazines: 2,5-dimethylpyrazine, 2,3,5-trimethylpyrazine, 2-ethyl-3,5-dimethylpyrazine, total pyrazines;
- off-notes: hexanal, nonanal, 2-pentylfuran, 1-octen-3-ol (the isolate's lipid background);
- furanones: DMHF, norfuraneol, furfural.
- Optional, if available: GC-O / AEDA on the 4 corner pots that differ in A and C (the perception layer).

## Sizing

Per-pot reproducibility of a within-study log ratio is taken as 0.10 - 0.15 dex (20 - 30 % CV per pot; PREREG sec. 5
uses 20 %, n = 2). In a 2^4 design with 2 replicates, each main effect is a difference of two means of 16 pots:
standard error 0.035 - 0.053 dex. A main effect of **0.15 dex (1.4x) is detected at alpha 0.05 (two-sided) with power
about 0.85 - 0.99** across that range; in a full factorial every two-factor interaction (e.g. ribose x thiamine,
yeast extract x pH) has the same standard error, so it is detected at the same size. The study is therefore not
limited by replicates; what it adds is one laboratory and one matrix to each cell below. Adding a third replicate
does not change what the table can show; a second laboratory running the same 36 pots would.

## What each result does to the tool

| map cell (family x lever) | today | after this study |
|---|---|---|
| thiols x sugar identity | 2 laboratories, listed | 3: can be answered (ribose vs glucose) |
| thiols x pH | 4 laboratories, direction unclear | 5, and the first plant-protein matrix in the cell |
| thiols x sulfur source (thiamine) | keyed, not connected to cysteine | connected: thiamine vs cysteine estimated |
| any family x yeast extract | no contrast on disk | 1 laboratory: listed (the map's experiment line still asks for 2 more) |
| Strecker / pyrazines / furanones x sugar, pH | 1 - 3 laboratories | one more each |
| off-notes x any lever | thin | the isolate's own lipid background, measured |
| engine, step 5 | no held-out new-laboratory contrasts | 16+ new-laboratory contrasts per odorant (TEST by PREREG sec. 2) |

Every prediction (the lever table's and the engine's) is frozen and sent before the data exist, with its pass rule.
A cell is "answered" only by the PREREG sec. 5 rule; the study cannot answer a cell on its own, only add its
laboratory.

## What the partner gets

Co-authorship on the data and on the validation; the open tool and the lever table run on their own recipes; the
analysis and the pre-registration done by us. They keep their data and decide when to publish.

## Open points for the owner

- Which group, and whether the full 36 pots or the 2^(4-1) half fraction (8 cells x 2 pots = 16 pots plus the checks;
  resolution IV: main effects clear of two-factor interactions, which are then aliased in pairs) fits their budget.
- Whether to add the pan-fried arm (180 C surface) at the two A x C corners: the cooking-method cell has no contrast.
- Whether the blind-test framing already used once (predictions frozen and sent before the data) applies here too.

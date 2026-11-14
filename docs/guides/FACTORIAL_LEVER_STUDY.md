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
prediction a new kitchen gets, one, two and three laboratories are nearly equal, and one is the best calibrated. Two
is the default only because step 5 (the engine against held-out contrasts from two TEST laboratories) needs two.

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

From the preposterior check, rerun on the 2026-11-14 table (99 sources, 44 laboratories; matrix shift 0.3 dex, A4
keying; means of 400 simulated outcomes; the 2026-11-12 numbers are in the sections below):

| | today | two laboratories | three laboratories |
|---|---|---|---|
| lever-table cells newly answered | - | 16.5 (12.4 with a literature estimate today) | 28.4 (13.5 with a literature estimate; the rest depend on the prior, read as optimistic) |
| experiment-map family x lever cells newly answered | 8 answered | 8.7 | 14.2 |
| cells where step 5 can test the engine (two TEST laboratories) | 9 | 57 | 57 |
| plant-protein matrix in the keyed record | none | every cell the study touches | every cell the study touches |

Every prediction (the lever table's and the engine's) is frozen and sent before the data exist, with its pass rule.
A cell is "answered" only by the PREREG sec. 5 rule.

## Preposterior check: does the study unblock anything? (2026-11-11, rerun on the fixed table 2026-11-12)

Before anything is sent, `scripts/generators/build_lever_design_preposterior.py` simulated each candidate design 200
times. Each run draws true effects from the lever table's own posterior predictive (cells with no keyed estimate draw
from the PREREG prior), adds a laboratory offset, a 0.1 dex within-laboratory residual and SIDA / GC-MS measurement
error, and optionally a pea-matrix shift of 0.3 dex shared by all new laboratories. It then keys the synthetic
contrasts and refits with the table's own model and reading rule. Full table:
[lever_design_preposterior.md](../../results/validation/lever_design_preposterior.md). Matrix shift 0.3 dex,
main-effect keying:

| design (36 pots unless stated) | cells newly answered, mean [90 %] | of which with a literature estimate today | experiment-map cells newly answered | cells with 2 TEST laboratories (step 5) |
|---|---|---|---|---|
| the draft above: one laboratory, 2^4 x 2 | 3.5 [2, 5] | 3.5 | 2.5 | 0 -> 2 |
| two laboratories, 2^3 x 2 each, no yeast extract | 9.7 [6, 13] | 6.5 | 6.0 | 0 -> 57 |
| three laboratories, 2^3 x 1 each, no yeast extract | 30.7 [25, 36] | 7.0 | 17.8 | 0 -> 57 |
| three laboratories, sugar x sulfur source 2 x 3, x 2 | 27.2 [21, 33] | 5.0 | 8.7 | 0 -> 57 |
| reference: the draft in two laboratories (72 pots) | 9.8 [7, 13] | 6.6 | 6.2 | 0 -> 57 |

What it shows:

1. **As drafted, the study barely unblocks anything.** One laboratory can only finish cells that already have two;
   there are about three to four. It gives step 5 (the engine against the lever table on held-out contrasts from two or more
   TEST laboratories) two cells. Its strength (power for 0.15 dex effects and interactions) is the wrong currency: the
   binding constraint is laboratories, and replicates inside one laboratory do not buy them.
2. **The same pots in two laboratories do about 2 - 3x better on every count, and give step 5 its 57 cells.** Running the
   whole draft in a second laboratory (72 pots) adds nothing over two laboratories at 36: the yeast-extract factor and
   the second replicate are spent on cells the table cannot read.
3. **A third laboratory mostly answers cells no paper has keyed** (31 against 10; the literature-backed count stays
   about 7). Those cells would be answered by the three new laboratories alone, and how many depends on the prior
   for effect sizes (mu ~ N(0, 0.5) dex). Read the 31 as optimistic; the 7 and the step-5 count barely depend on it.
4. **Yeast extract feeds no cell.** PREREG sec. 4 has no yeast-extract lever class. Keep it only after an amendment
   adds one (categorical, reference "none") before any data are read.
5. **Keying a factorial pair by pair would overstate the answer** (per-pair keying: up to +5 cells). The categorical
   network has one heterogeneity scale per pair and no laboratory random effect, so a laboratory's 2^(k-1) contrasts
   for one lever would count as independent draws. Key a designed study as one contrast per laboratory per lever (the
   mean over the other factors), and state that in the pre-registration. The one cell the literature answers today
   (MFT, no sulfur source vs cysteine; 7 contrasts, 4 from one laboratory) survives this: refitted with one contrast
   per laboratory, it stays answered in all 8 subsamples (estimate -0.81 to -1.21 dex).

**Redesign (done the same day, owner-approved):** the design section above is now the 2^3 (sugar x sulfur source x pH)
with checks, keyed by PREREG amendment A4. The laboratory count was then revised by the predictive check below: two by
default, one acceptable, three not worth it.

## Predictive check: does the number of laboratories matter? (2026-11-11, rerun 2026-11-12)

The 'answered' count above is partly circular: it counts how often new laboratories clear our own three-laboratory
rule. The non-circular question is how much better the table predicts **a new kitchen cooking in pea isolate**. The
same simulated worlds were scored for the literature alone and for the same 36 pots in one (2^3 x 4), two (2^3 x 2) or
three (2^3 x 1) laboratories ([the last section of the record](../../results/validation/lever_design_preposterior.md);
200 worlds per setting). The 57 cells the 2^3 touches, matrix shift 0.3 dex:

| | literature alone | one laboratory | two laboratories | three laboratories | floor: the pea mean known exactly |
|---|---|---|---|---|---|
| sign wrong for a new kitchen, all cells | 0.45 | 0.31 | 0.30 | 0.29 | 0.17 |
| the same, 12 cells with a literature estimate | 0.27 | 0.25 | 0.24 | 0.24 | 0.18 |
| the same, 45 cells with none | 0.50 | 0.33 | 0.31 | 0.30 | 0.17 |
| 90 % interval coverage, cells with none (nominal 0.90) | 0.88 | 0.86 | 0.81 | 0.77 | - |
| the same, design laboratories sharing our protocol (half the scatter) | 0.88 | 0.89 | 0.78 | 0.67 | - |

The first run of this check found the table's intervals covering only about 60 % in sparse cells. The cause was a bug
in the lever table: the scale prior on its log-spaced grid was not the declared half-normal (it put 61 % of its mass
below 0.1 dex instead of 16 %). Fixed in PR #124 (readings unchanged, intervals wider); the numbers above are the rerun.

What it shows:

1. **The study's value is pea-matrix data where the literature has none.** It cuts the sign error on those 45 cells
   from a coin flip (0.50) to about 0.31. On the 12 cells the literature already estimates, it helps a little
   (0.27 -> 0.24).
2. **How the 36 pots are split across laboratories barely matters for accuracy** (0.31 / 0.30 / 0.29).
3. **More laboratories make the intervals worse calibrated, not better.** Laboratories that share our protocol and lot
   agree with each other, so the fit narrows its interval below the spread a new kitchen has: 0.89 coverage with one
   laboratory, 0.67 with three. One laboratory is the best calibrated.
4. **The floor is high.** Even a predictor that knew the pea-matrix mean exactly would get the sign wrong for about
   17 % of new kitchens, because laboratories scatter that much. No study design removes that.

**So the laboratory count is decided by step 5, not by the table.** For the table's own predictions, one laboratory
(2^3 x 4) is as accurate as two or three and the best calibrated. Two laboratories are worth the coordination only
because step 5 (the engine judged on held-out contrasts from two TEST laboratories) cannot run on one. The design
section keeps two as the default for that reason; if step 5 is not wanted from this round, one laboratory is the
better and cheaper design.

## Changes after steps 5a - 5b (2026-11-13)

The disagreement audit ([DISAGREEMENT_READING.md](../reference/lever_table/DISAGREEMENT_READING.md)) and the step-5a
coverage check ([step5a_display_check.md](../../results/validation/step5a_display_check.md)) found why published
contrasts fail to add up; the design already avoids most of it, and four points are added:

- **pH, buffered or not.** The pea isolate is the only buffer; no buffer salt is added (state "none" in the report).
  pH is measured before heating and after cooling in every vial (as above), and each laboratory's pH contrast is
  keyed with its final pH beside the initial one, so that the table can read it as the lever it really was. If the
  two pH levels end less than 0.3 units apart in a laboratory, that laboratory's pH contrast is reported but read
  "pH not held".
- **Reference arms.** Every factor's two levels are both compounds, and the matrix-background pot is the same in
  every laboratory, so no "none" arm means a different pot in different studies (the cause of 9 of the 16 widest
  disagreements in size).
- **Detection limits.** A compound not detected is reported with that run's limit, never as zero or blank.
- **The isolate's composition.** The engine cannot run a pot whose matrix composition is not on file (148 of step
  5a's 396 held-out contrasts), and no dossier on disk gives free sugars or free amino acids for a pea isolate. One
  laboratory measures them (and total cysteine and methionine, and fat) once on the shipped lot; without them the
  study still serves the lever table, but the engine cannot be judged on it.

## Changes after steps 13 - 15 (2026-11-14; STATUS_AND_PLAN step 16)

Three literature steps tried to fill the map's cells one laboratory short before asking anyone to cook. They settle
what this study must carry:

- **The literature cannot fill this study's levers.** A fixed 200-paper search for pH cells (step 13) found no
  qualifying paper; one for cysteine added vs none (step 14) found at most two. Only time and temperature series were
  fillable (step 15: seven papers, two cells answered). The study's three factors (sugar, sulfur source, pH) are
  exactly the ones the literature does not supply, so none of its pots is redundant with a search.
- **Add one buffered-vs-unbuffered pair (A6's decisive measurement).** PREREG amendment A6-R1 (does buffering change a
  pH cell?) cannot be scored on the literature: no laboratory runs both control states at matched pH
  (`docs/reference/decisive_measurements.yml`, `ph_control_buffered_vs_unbuffered`). At the glucose + cysteine corner,
  each laboratory adds the two pH levels in 0.1 M phosphate (buffered) beside the factorial's own matrix-only pots
  (adjusted, unbuffered), 2 pots each: **4 pots per laboratory** (22 per laboratory and 44 in all with two laboratories), final pH measured in every vial. Within one
  laboratory and one lot, the buffered and unbuffered slopes differ only in control state, which is the comparison
  the literature confounds with laboratory. It also gives the buffer cells (phosphate vs none, one laboratory today)
  a designed contrast.
- **Time and temperature stay out of the factorial.** Step 15's cells now read "direction unclear" on 3 - 5
  laboratories (FFT and MFT x time at 100 - 130 C among them): more laboratories of the same kind do not settle them,
  and a fourth factor halves the pots per cell. An optional time pair (10 vs 40 min at the glucose + cysteine corner,
  2 pots each) is listed under the open points.
- **New laboratories must be new to the keyed record.** The lever table now holds 44 laboratories; a candidate must be
  absent from `docs/reference/contrast_benchmark/split.yml` and from `docs/reference/lever_table/new_labs.yml`, or it
  adds no laboratory to any cell.

### Sizing rerun on the 2026-11-14 table (step 16a)

`build_lever_design_preposterior.py` was rerun unchanged on the table after steps 13 - 15 (99 sources, 2,581
contrasts, 44 laboratories), with 400 simulated worlds per setting instead of 200
([record](../../results/validation/lever_design_preposterior.md)). Matrix shift 0.3 dex, main-effect keying:

| design (36 pots unless stated) | cells newly answered, mean [90 %] | of which with a literature estimate today | experiment-map cells newly answered | cells with 2 TEST laboratories (step 5) |
|---|---|---|---|---|
| one laboratory, 2^4 x 2 (the first draft) | 4.1 [3, 5] | 4.1 | 2.0 | 9 -> 24 |
| two laboratories, 2^3 x 2 each (this design) | 16.5 [13, 20] | 12.4 | 8.7 | 9 -> 57 |
| three laboratories, 2^3 x 1 each | 28.4 [23, 33] | 13.5 | 14.2 | 9 -> 57 |
| three laboratories, sugar x sulfur source 2 x 3, x 2 | 31.3 [26, 37] | 13.9 | 8.0 | 3 -> 57 |
| reference: the first draft in two laboratories (72 pots) | 16.6 [13, 20] | 12.5 | 8.8 | 9 -> 57 |

Predictive check, the 57 cells the 2^3 touches (now 29 with a literature estimate and 28 without, against 12 and 45
on 11-12):

| | literature alone | one laboratory | two laboratories | three laboratories | floor: the pea mean known exactly |
|---|---|---|---|---|---|
| sign wrong for a new kitchen, all cells | 0.38 | 0.28 | 0.28 | 0.28 | 0.17 |
| the same, 29 cells with a literature estimate | 0.26 | 0.25 | 0.24 | 0.24 | 0.18 |
| the same, 28 cells with none | 0.50 | 0.32 | 0.32 | 0.31 | 0.17 |
| 90 % interval coverage, cells with none (nominal 0.90) | 0.88 | 0.79 | 0.75 | 0.73 | - |
| the same, design laboratories sharing our protocol (half the scatter) | 0.88 | 0.83 | 0.72 | 0.65 | - |

What changed and what did not:

1. **The two-laboratory design now answers about 16 cells, not 10.** The literature steps put more cells one or two
   laboratories short, so two new laboratories finish more of them; 12 of the 16 have a literature estimate today,
   so the count leans less on the prior than it did.
2. **The literature alone is better than it was** (sign error over the 57 cells 0.45 -> 0.38), because 17 more of them
   now have a literature estimate. The study's gain is still where the literature has nothing: 0.50 -> 0.32.
3. **The laboratory split still barely matters for accuracy** (0.28 for one, two or three), and more laboratories
   sharing our protocol still under-cover (0.83 / 0.72 / 0.65 with half the scatter). The PREREG for this study must
   therefore score coverage on laboratories outside the protocol, or widen the design laboratories' scale, before
   its intervals are shown to a user (premortem item for 16c).
4. **One laboratory now gives step 5 24 cells** (9 already have two TEST laboratories from steps 13 - 15). Two
   laboratories still give all 57, so the default stays two; one remains acceptable.
5. The buffered pair (4 pots per laboratory) is outside the simulation: it serves A6-R1 and the buffer cells, not
   the counts above.

## Premortem (2026-11-14; step 16c)

An independent reviewer imagined the study had failed by 2027-12 and listed twelve causes, ranked by likelihood
times impact. Everything below is decided before anything is sent; the rest of this document is history wherever it
disagrees.

**Gates before sending** (new sub-steps, STATUS_AND_PLAN step 16):

1. **16e, engine dry run.** The case for two laboratories rests on step 5. But "two TEST laboratories" does not make a
   cell predictable: in step 5a the engine had no target for 2- and 3-methylbutanal, phenylacetaldehyde or
   1-octen-3-ol, refused some arms, and lacked matrix composition for others. Run the engine on the 8 corners with a
   placeholder pea composition and count the cells it predicts. If too few cells remain to judge it, the default
   becomes one laboratory, 2^3 x 4, which the predictive check found as accurate and better calibrated.
2. **16f, PREREG amendment A7**, frozen before sending:
   - **Censored pots in A4 keying.** Glucose arms will often leave MFT, FFT and the disulfide below the limit.
     - A censored pot enters at its limit, through a censored-likelihood mean.
     - An arm with more than half its pots censored is keyed, but enters no cell.
     - The count is reported.
   - **pH band by the nominal level.** This needs the move to 6.0 / 7.0 (below).
     - The "pH not held" rule (the levels end less than 0.3 units apart) moves from this document into the PREREG.
   - **Corner sign check.** Where the four corner contrasts of a lever disagree in sign, the cell is flagged
     "interaction-dominated" and the corner contrasts are reported beside the A4 mean.
   - **Shared-protocol scale.** The two design laboratories share one protocol and one lot, and the sizing showed
     their intervals under-covering (0.72 against 0.90). They therefore get a shared protocol-cluster effect.
     - A pea-only estimate is reported beside the pooled one.
     - Its interval is never shown to a user as a pea answer until a laboratory outside the protocol has scored it.
   - **The buffered pair's own rule** (below). It is separate from A6's confirmation set: A6 needs at least 3
     laboratories, and two design laboratories can never meet that alone.
3. **16g, rerun the sizing on each laboratory's real analyte panel.** The simulation assumed all 19 analytes in both
   laboratories with no censoring. Some candidates have no thiol SIDA, and H2S is optional. Rerun once the two
   laboratories have confirmed their panels.

**Design changes:**

- **pH levels move to 6.0 and 7.0** (from 5.5 and 6.5), for two reasons:
  - Pea isolate aggregates near its isoelectric point (about 4.5 - 5), so at 5.5 the pH factor would also change
    protein state and volatile binding.
  - A 5.5 / 6.5 midpoint sits on the 6.0 band edge of the lever table (`PH_BANDS`, lower bound inclusive).
  6.0 / 7.0 stays inside the 6 - 7.5 band and near the isolate's native pH. Adjust with HCl or NaOH after a stated
  hydration time, and read the pH at 25 C.
- **Heating protocol named:**
  - pressure-rated crimp vials with a stated fill;
  - one heating medium per laboratory, stated;
  - the 20 min hold timed from a thermocouple dummy vial reaching 138 C;
  - every heating batch holds all 8 cells.
- **Sample handling:**
  - The labelled standard goes in through the septum straight after the ice quench.
  - The time from quench to extraction is fixed, and so is homogenisation.
  - Every vial is analysed within one stated window, or all are frozen at -80 C the same way.
  - A pooled QC injection runs every 6 samples.
- **Core and optional analytes.**
  - The core panel is MFT and FFT (SIDA, or a recovery shown in this matrix), methional, 2- and 3-methylbutanal,
    phenylacetaldehyde, 2,5-dimethylpyrazine, hexanal, furfural and DMHF.
  - Everything else, H2S included, is optional and declared in advance.
  - Labelled MFT and FFT standards are offered to a laboratory without them.
- **The no-isolate pot (glucose + cysteine without pea isolate, the matrix contrast) is made up in 0.1 M phosphate at
  pH 6.0, 2 pots in every design.** As drafted it was in water, which drifts, and the two-laboratory design had one
  pot. The matrix-background pot is also 2 pots, set to pH 6.0.
- **The buffered pair is re-scoped.**
  - 15 % pea isolate buffers strongly itself, and phosphate catalyses the reaction more at higher pH (its HPO4(2-)
    fraction rises with pH, pKa2 7.2).
  - The pair therefore measures "phosphate-buffered vs matrix-buffered pH slope". That is the literature's real
    contrast, with catalysis and control state together. It does not separate the two.
  - At 2 pots per arm, the difference of slopes has a standard error of about 0.1 - 0.15 dex pooled over both
    laboratories.
  - Its rule is an equivalence margin: |difference of slopes| < 0.2 dex per pH unit reads "same", the 90 % interval
    entirely outside it reads "different", and anything else reads "unresolved".
  - A water model system (glucose + cysteine, unbuffered vs phosphate) would test the literature's drift confound
    directly. It is listed under the open points (4 pots), not added.
- **Independence confirmed before the protocol goes out.** The criteria:
  - its own method lead, instrument and standards source;
  - no co-author from a laboratory already in the record;
  - (`new_labs.yml`'s collaboration rule).
  The ZHAW (Steinhaus co-authorship) and Ohio State (shared university) candidates are confirmed against this, and
  a third laboratory is kept as backup.

**The brief, made honest:**
- the effort is weeks of method work, not "a few days of instrument time";
- "about 12 questions" newly answered with a literature estimate, not 16;
- the core and optional analyte split;
- an embargo: a partner's data are keyed on a private branch until the partner publishes;
- the composition analysis is listed as part of the ask.

Pot count (two laboratories): 2^3 x 2 = 16, matrix background 2, no-isolate 2, buffered pair 4: **24 per laboratory
(18 before), 48 in all.**

## Gates 16e and 16f (2026-11-14): one laboratory, A7 frozen

**16e, the engine dry run** ([p3_engine_dry_run.md](../../results/validation/p3_engine_dry_run.md),
`scripts/generators/build_p3_engine_dry_run.py`). The shipped engine was run on the 8 corners through its own front
door, with A4 main effects. In 15 % pea isolate (protein sites only, since no dossier gives its free sugars or amino
acids) it predicts **3 of the 57** lever x odorant cells in full and 10 in part. With the pH held it predicts 6 in
full; with the added precursors in water, none. Why:
- 7 of the 19 analytes have no engine target: 2- and 3-methylbutanal, phenylacetaldehyde, 1-octen-3-ol, two
  alkylpyrazines and the pyrazine total.
- The lipid and Strecker lanes refuse the isolate arm.
- The thiol lane refuses most ratios.

Some of the ratios it does return are implausible (DMHF -10 dex for ribose against glucose, norfuraneol +29 dex).
This is an engine follow-up, not a design input.

**Decision: one laboratory, 2^3 x 4.** Step 5 was the only reason for two laboratories (the predictive check
above), and it cannot run on 3 cells. A7-R4 also counts two laboratories sharing one protocol and one lot as one
cluster toward the three-laboratory rule. So the two-laboratory "answered" counts (16.5) were overstated, and a
second laboratory under the same protocol buys little.

For the table's own prediction, one laboratory is as accurate as two:
- sign error 0.28 over the 57 cells, 0.32 on the 28 with no literature;
- the best calibrated (coverage 0.79, or 0.83 with half the scatter).

It newly answers about 4 cells (the one-laboratory row of the preposterior table; the 2^3 x 4 was not run for that
count): the cells already two laboratories short.

The better second laboratory is a later, independent replication under its own protocol. That counts as a new
laboratory under R4, and it is what R4 needs before any interval is shown as a pea answer.

**Pots in the one laboratory:**

| | pots |
|---|---|
| factorial, 8 cells x 4 | 32 |
| matrix background (pH 6.0) | 2 |
| no isolate (glucose + cysteine in 0.1 M phosphate, pH 6.0) | 2 |
| buffered pair (0.1 M phosphate at the glucose + cysteine corner, pH 6.0 and 7.0, 4 each; its matrix-only side is the factorial's corner) | 8 |
| **total** | **44** |

**16f, PREREG amendment A7** is frozen in [PREREG.md](../reference/lever_table/PREREG.md). It was revised once,
before freezing, after an independent review of its first draft. Its rules:
- R1, censored pots;
- R2, the pH band at the nominal levels and the 0.3-unit holding rule;
- R3, the corner sign check, at about 0.1 % false flags;
- R4, the protocol cluster;
- R5, the buffered pair.

R5's power is stated now. With one laboratory, the pair can show a large difference (P about 0.64 at 0.4 dex per pH
unit) but is not expected to show sameness (P about 0.28 at D = 0).

**Still open:**
- the R1 - R5 code in `build_lever_table.py`, written and tested on synthetic data before any P3 value is read;
- 16g, the sizing rerun on the laboratory's confirmed analyte panel, with R4 in the simulation.

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
- Whether to add the optional time pair (10 vs 40 min at the glucose + cysteine corner, 2 pots per laboratory).

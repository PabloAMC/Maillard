# A designed Maillard flavour study in pea protein: one-page brief for partner laboratories

*Draft, 2026-11-14 (STATUS_AND_PLAN sec. 0, step 16). For a food-chemistry group with quantitative GC-MS. The full
design, its sizing and its checks are in [FACTORIAL_LEVER_STUDY.md](FACTORIAL_LEVER_STUDY.md). Nothing has been sent.*

## The problem

Formulators of plant-based meat ask simple questions: if I switch glucose for ribose, cysteine for thiamine, or raise
the pH by one unit, what happens to the meaty thiols, the Strecker aldehydes, the pyrazines and the off-notes? We
keyed the published literature as a table of such levers. Two blind readers read every paper: 99 sources, 2,581
within-study contrasts, 44 laboratories ([lever table](../../results/validation/lever_table.md)). It answers only
**8 of 75** of these questions ([experiment map](../../results/validation/experiment_map.md)), for three reasons:
- almost every question has been studied by one or two laboratories;
- no keyed study cooks in a pea or soy isolate;
- laboratories disagree by more than most levers move.

Searching the literature further does not fix this for sugar, sulfur source or pH (we tried, with fixed,
pre-registered searches). Only new laboratories cooking in the users' matrix can.

## What we ask

One small factorial, the same in each of two independent laboratories (one laboratory with more replicates is the
fallback, decided before we write to you).

- **Base:** 15 % w/w pea protein isolate in water. Heated in pressure-rated sealed vials at 140 C for 20 min, timed
  from a thermocouple vial. We ship the isolate from one lot to both laboratories.
- **Factors**, 8 cells x 2 vials:
  - sugar: glucose vs ribose, 30 mM;
  - sulfur source: cysteine vs thiamine, 10 mM;
  - pH set before heating: 6.0 vs 7.0.
- **Checks:**
  - a matrix-background pair;
  - a no-isolate pair in phosphate buffer;
  - one phosphate-buffered pair at both pH levels: does a buffer salt change the pH effect seen in the matrix alone?
    We found no published laboratory that tested this in one study at matched pH.
- **24 vials per laboratory.** Measure pH before heating and after cooling in every vial.
- **Composition:** one analysis of the shipped lot (free sugars, free amino acids, total cysteine and methionine,
  fat). We can arrange it with a contract laboratory if you prefer.
- **Analytes:**
  - **Core panel:** 2-methyl-3-furanthiol and 2-furfurylthiol (stable isotope dilution, or a recovery shown in this
    matrix; we can supply labelled standards), methional, 2- and 3-methylbutanal, phenylacetaldehyde,
    2,5-dimethylpyrazine, hexanal, furfural, DMHF.
  - **Optional, declared in advance:** disulfides, other alkylpyrazines, nonanal, 2-pentylfuran, 1-octen-3-ol,
    norfuraneol, H2S.
- **Handling:** the internal standard goes in straight after the ice quench; the time to extraction is fixed. A
  compound not detected is reported with that run's detection limit.
- **Your own method is welcome** within these points. That is part of what a new laboratory contributes.
- **Honest effort:** weeks of method set-up and a few days of instrument time, mostly for the thiols.

## What you get

- Co-authorship on the data and on the validation.
- The analysis and the pre-registration, done by us. Every prediction (the lever table's and our kinetic model's)
  is frozen and published before your data exist.
- The open tool, run on your own recipes.
- You keep your data: they are keyed on a private branch and made public only when you publish.

## Why it is worth it

In our simulations (`results/validation/lever_design_preposterior.md`, rerun 2026-11-14), two laboratories cut the
chance that the table gets an effect's direction wrong for a new kitchen:
- from one in two to about one in three, on the 28 questions where the literature has nothing in any matrix;
- from 0.38 to 0.28, over all 57 questions the study touches.

About 12 questions that already have a literature estimate would become answered. More would be answered among the
rest, but that count depends on our prior and we do not quote it.

Some error no study can remove: laboratories scatter enough that even the true pea-matrix mean gets the sign wrong
for about one new kitchen in six.

The vials also give the first held-out test of a mechanistic Maillard model on plant-protein data.

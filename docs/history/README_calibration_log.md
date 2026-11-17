# README calibration log, 2026-09-10 to 2026-09-11

*Moved verbatim from the README's "How well calibrated is it?" section on 2026-11-17; only the relative link
paths were changed so they resolve from here. It is the dated record of how the panel's denominator moved in
September 2026. The current numbers are in the [README](../../README.md#how-well-calibrated-is-it).*

**Read the denominator, 2026-09-10.** It moved twice this week, in opposite directions, and the
rate below is worthless without both halves.

It first went from 39 rows to 46, because seven questions stopped being refused: a branch fraction
measured in 1978 and 1981 was read into the model, and a missing unit conversion was found. The
numerator did not move, so **the pass rate FELL on changes that made the model strictly more
capable.** That is what happens when a model answers more questions, and it is why a bare rate is a
poor headline.

It then went from 46 back to **39**, and this direction deserves more suspicion, because it flatters
the rate. Seven rows were refused on the finding that **four of the panel's pots were never cooked**:
their "40 °C for 10 minutes" is the headspace incubation of the measurement itself, and each
bundle's own vessel record says so. Asking a model of what heat FORMS to reproduce what a raw flour
ARRIVED WITH is not a test of its chemistry, and the misses were 3357x to 42301x. Refusing them
raises within-3x from 7/46 to 7/39 **with no prediction improving at all**. Three things are offered
against the obvious objection, and the full argument is in
[`kinetic_core_b31_prereg.md`](../../results/validation/kinetic_core_b31_prereg.md): the criterion looks
only at the pot's conditions and never at an error; it needs two independent declarations to agree;
and **every lipid miss in a pot that WAS cooked still stands**, including this panel's largest, an
alkylfuran at 366x. The control is Bi 2020, one laboratory measuring one compound in raw and in
roasted pea: raw missed by 3717x and is refused, roasted misses by 3.7x and is scored.

**A fourth, the same day, and it came from asking whether the third was a one-off.** It was not.
Cross-referencing all 48 benchmark bundles against the 245 PDFs on disk found **two more** whose
papers print compounds the bundle never scored — and, in passing, one bundle that contradicted itself
about whether its own source was on disk, and one stale note left behind by the fix a few hours
earlier. Eight more measurements went on the panel.

Five of the eight came back as **refusals**, not misses, and that is the useful part. Three are
compounds whose waves were refused and whose steps ship inert, so the model says so by name instead of
returning a number. Two more exposed something worse: a pot that charges only a protein isolate was
**answering 0.0** for any Maillard compound — an absence of a prediction reported as one — because the
isolate is a lipid carrier and charges nothing on the other networks. This model's own record had
named that failure two waves earlier and the code kept doing it. It now refuses, and a test asserts
that **no scored row anywhere has a prediction of exactly zero**.

**A third move, 2026-09-11, and it is the one worth reading.** Five rows were added to the panel —
not measured, not modelled, just *noticed*. A hold-out bundle scored **one** of the six species its
own paper measures, and the paper has been on this disk since 7 September under a different author's
name than the bundle gives it. Scoring the other five moves the rate to **9 of 44** and makes the
median fold error **worse**, 9.3x to 10.6x, because three of the five are 32x, 35x and 63x. Nothing
about the model changed; the panel simply stopped ignoring measurements it already had.

What those five say is the sharpest result on this page. In an amine-free pot the model's
3-deoxyglucosone is right to **11 %** and its methylglyoxal to **28 %** — and then the very next
species, 3,4-dideoxyglucosone, is **32x low**. The hydroxymethylfurfural deficit is not the entry to
the sugar chemistry; it is one step, `k_tdg_ddg` — and that is one of four constants whose interval
this model declined to widen *because two laboratories agree on it*. Both of them measured it dry.
Agreement is not accuracy when both laboratories share a matrix.

The numerator did move, once: **three rows entered the 3x band together** when a beverage's own
unheated column was declared as what the pot started with. Between 36 % and 42 % of every level that
paper prints was in the beverage before any heat, and the model had been charged for the raw
material. None of the three is a fit row.

**Do not read those three as chemistry.** A declared starting level is a measurement handed to the
model, so the fold error on the total grades the model partly on the number it was given. Split the
declaration out of both sides and the picture changes:

| compound | fold on the total | the declared part is | fold on what the cook FORMED |
|---|---:|---:|---:|
| hexanal | 2.2x | **93.5 %** of the prediction | **19.7x** |
| 2-pentylfuran | 2.5x | **92.2 %** | **20.6x** |
| nonanal | 1.5x | 52.8 % | **2.1x** |

Two of the three are inside the band on about six per cent of their own answer, and on the part the
chemistry is responsible for the lipid lane still under-predicts by about twentyfold — which is where
it already was. What improved is the scoring, not the model. Only nonanal is a chemistry result.
Every scored row now carries this split (`declared_share_of_prediction`, `fold_error_formed_only`)
so the total can never be quoted without it. The median fold error falls from 31x to 9.3x across all
of this, and most of that fall is the refusals, not the fix.

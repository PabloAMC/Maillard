# Pre-registration: wave B28, the two compounds the lipid lane refuses (written 2026-09-09, before anything ran)

## 1. Why

The scoring panel refuses seventeen rows outright. Seven of them are two compounds, and the lane's
own registry says exactly why in `PROHIBITED_DERIVATIONS`:

> **oleate → nonanal branch fraction.** "Measured nowhere in the fit corpus. Frankel 1989 fed
> linoleate only. The retired screening lane's shipped 0.15 had no source. Requests for absolute
> nonanal in an oleate-bearing matrix are REFUSED."
>
> **linoleate → 2-pentylfuran branch fraction.** "Not in Frankel's six-product slate and measured
> nowhere else in the corpus. The retired screening lane's shipped 0.08 had no source."

The lane was built so that these refusals could be lifted by a number rather than rewritten:
`NONANAL` is already a species, `LOOH_OL` is already the only pool with a nonanal edge, and that
edge's branch fraction is `None`. Two papers fetched on 2026-09-09 print both numbers, from the
**same laboratory and the same injector-port thermolysis method** as the Frankel & Gardner 1989
slate the lane is fitted on.

## 2. What enters, and in what form

| what | value | source | form |
|---|---:|---|---|
| oleate → nonanal, autoxidised pool | 15 % of total volatile peak area | Frankel 1981 Table II | **`[C]`, republished from Selke 1978** |
| oleate → nonanal, photosensitized pool | 10 % | Frankel 1981 Table II | `[M]`, the only independent determination |
| linoleate → 2-pentylfuran, autoxidised | 2.4 % | Frankel 1981 Table III | `[M]` |
| linoleate → 2-pentylfuran, photosensitized | 0.6 % | Frankel 1981 Table III | `[M]` |
| **2-pentylfuran / hexanal, same chromatogram** | **0.160** and **0.0353** | derived from Table III | **the form to ship** |

**Two traps this wave must not fall into, both named now.**

**(a) The two 15 % figures are ONE measurement, not two agreeing ones.** Frankel 1981's autoxidised
oleate column is footnoted "Data from ref. 21", which is Selke 1978. Anyone comparing the two papers
and finding 15 % twice is reading the same number twice. The independent second determination is the
photosensitized 10 %, on a different hydroperoxide distribution, and the honest statement of the
oleate nonanal share is therefore **one measurement at 15 % and one at 10 %, not a replicate pair.**

**(b) The 1981 and 1989 denominators are different.** Frankel 1989's shares are fractions of **six
measured peaks**; Frankel 1981's are fractions of a **whole ~20-peak chromatogram**, including 9.9 to
12 % explicitly unidentified. A 1981 share is therefore systematically smaller than a 1989 share for
the same product, by the ratio of the denominators, and the two must never be pooled as printed. For
2-pentylfuran the wave ships the **ratio to hexanal**, which is denominator-free and hangs the
alkylfuran on a node the lane already models.

## 2b. A problem found while building, before any result existed

Written the same evening, after reading the lane's integrator and before running anything. It is
recorded here rather than discovered later because it changes what this wave can honestly claim.

**Frankel 1981 prints peak-area shares, not yields.** No internal standard, no response factors, no
replicates, no stated error. The lane converts its LINOLEATE distribution into moles per
hydroperoxide by anchoring it to Schroen's separately measured hexanal yield. **There is no such
anchor for oleate anywhere in the corpus.** So the oleate column gives nonanal's SHARE of the oleate
slate and cannot give moles of nonanal per mole of oleate hydroperoxide.

The consequence, stated before the run: **an absolute nonanal prediction still requires a declared
assumption**, namely that the named-product molar yield per oleate hydroperoxide equals the measured
one per linoleate hydroperoxide. That is not a measurement. This wave therefore ships it as a
DECLARED ASSUMPTION with a band, propagated into the interval on every nonanal answer and carrying a
mandatory warning, on the precedent of the `soy_paste_hong` protein loading — and the refusal changes
from "no branch fraction exists" to "the share is measured, the molar anchor is assumed". Those are
different states and the layer must say which one it is in.

**Prediction 1 of section 5 is therefore revised down, before the run, from 80 % to 55 %**, and the
revision is recorded rather than made silently. The share being measured does not by itself make an
absolute answer legitimate, and if the panel's rows demand absolutes then fewer of them lift than the
first draft of this pre-registration assumed. 2-pentylfuran is unaffected: it ships as a ratio to
hexanal, which is denominator-free and needs no anchor.

## 3. What this moves

`NONANAL`'s branch fraction stops being `None`, so requests for absolute nonanal in an
oleate-bearing matrix stop being refused. `2-pentylfuran` becomes a species with one incoming edge
from the linoleate pool. Nothing else in the lane changes: no rate constant, no barrier, no
temperature dependence, and the 180 °C Frankel 1989 fit that sets the lane's constants is untouched.

## 4. What counts as success, declared before the run

- **T1 arithmetic.** Both new slates reproduce their printed columns; the nonanal edge is a number;
  2-pentylfuran has exactly one incoming edge; the unidentified remainder of each 1981 column is
  routed to `LIPID_FRAG_C` exactly as the lane routes Frankel 1989's unclosed carbon.
- **T2 the refusals, decisive.** The panel's refused-row count falls, and every row that changes
  changes from REFUSED to answered. No row may go the other way, and no currently answered row may
  move by more than 0.05 dex.
- **T3 the cross-laboratory check, and it is free.** Renormalise Frankel 1981's linoleate column onto
  Frankel 1989's six-product slate and compare product by product. Same laboratory, eight years
  apart, 210 °C against 180 °C, and one product (the C13 oxo-ester) that 1981 could not identify for
  want of an authentic reference. Report every fold. This is the first time the lane's own fit source
  has had any external check at all.
- **T4 the newly answered rows, reported.** Do they land within 3x of measurement? Reported, not
  decisive: a branch fraction from a neat hydroperoxide at 210 °C predicting a real food is a long
  transfer and this wave does not pretend otherwise.
- **T5 nothing else moves.** The kinetic panel's sulfur, trunk and acrylamide rows bit-for-bit
  unchanged; the matrix layer untouched.

Ship rule: **SHIP if T1, T2 and T5 hold.** T3 and T4 are reported.

## 5. Predictions, before the run

1. T2 holds and at least four refused rows become answered. **80 %.**
2. T3's folds are worse than 2x on at least two of the five comparable products. **60 %.** The
   denominators differ, the temperatures differ by 30 °C, and 1981 reports no replicates and no error
   at all. If they agree better than that, it is a genuinely reassuring result for a lane that has
   never been checked.
3. T4 is bad: the newly answered rows miss by more than a decade. **70 %.** These are neat
   hydroperoxides pyrolysed in an injector port, and the panel's rows are foods. Lifting a refusal is
   not the same as being right, and the wave claims only the first.
4. Nobody will later mistake the 15 % for two agreeing measurements, because trap (a) is written
   into the parameter's own note and not only here. **90 %.**

## 6. Outcome (2026-09-09, run the same evening) — SHIP, with one half of it withdrawn during the run

**Verdict SHIP.** T1, T2 and T5 held; T3 and T4 are reported. Three refused rows lifted, not seven,
and the four that did not lift are the more instructive half.

| | before | after |
|---|---:|---:|
| refused panel rows | 25 | **22** |
| nonanal, Trikusuma pea beverage | REFUSED | **7.37 against 24 measured, 3.3x** |
| nonanal, Li 2026 soy/wheat extrudate | REFUSED | **23.9 against 72.7, 3.0x** |
| nonanal, Liu 2023 pea isolate | REFUSED | **0.109 against 0.802, 7.3x** |
| 2-pentylfuran, four rows | REFUSED | **still REFUSED, for a different reason** |

**Prediction 3 was wrong, and wrong in the good direction.** It said at 70 % that the newly answered
rows would miss by more than a decade, because these are neat hydroperoxides pyrolysed in an injector
port and the panel's rows are foods. Nonanal lands at 3.0x, 3.3x and 7.3x on three matrices from
three laboratories. That is better than most of what this model does, on a share measured in 1978 and
an anchor this wave declared rather than measured. It should not be over-read: three rows, and the
anchor's band spans a factor of five.

**What was withdrawn, and why it is the useful part.** The first run lifted 2-pentylfuran too, and
answered it **six to nine orders of magnitude** below measurement. The branch fraction is not the
problem; it is measured and it is right. The problem is that on the matrix-only path **the hexanal
these rows are scored against does not come from the lipid lane at all** — the lane's own hexanal in
the same pot is about 1e5 smaller — so an alkylfuran hung off the lane's hexanal is nearly zero. By
this layer's own rule a degenerate value is the absence of a prediction dressed as one, so the
refusal was restored with a sharper reason that now names what would lift it: a lipid charge these
matrices can actually integrate, which is the same gap the hexanal rows already carry.

**The ship rule was missing the test that catches this, and now has it.** A rule that counts refusals
falling would have called the first run a success. T2 now also requires that every row lifted out of
REFUSED be answered within three decades of its measurement. A lift into a near-zero is a regression
in honesty, not a gain in coverage.

**T3, the first external check this lane's fit source has ever had.** Frankel 1981 against Frankel
1989, renormalised onto the five products both quantify, agrees **within 1.6x on every one**. It is
not a two-point Arrhenius: 210 °C neat against 180 °C in hexane, a 25 °C column start against a
−65 °C cryotrap, packed against capillary. Temperature and light-end loss are confounded and these
two papers cannot separate them. Prediction 2 said at 60 % that at least two products would be worse
than 2x; none was.

**Both traps held.** The two 15 % nonanal figures are recorded in the species note as one measurement
republished, and the two slates are never pooled as printed.

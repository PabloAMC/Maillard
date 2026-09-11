# Pre-registration: wave B34, the five observables already on disk (written 2026-09-11, before anything was changed)

## 1. Why, and how it was found

B33 was registered yesterday as a *prediction*: that the trunk's three amine-free sugar entries are
about ten times too slow in water, on the evidence that scaling them ×10 moves four of five
hydroxymethylfurfural rows inside 3×. **That prediction rested on watching one observable.** Asked
whether the unlocks really needed new reading, a search of the 245 PDFs already on disk found
`data/articles/Leitzen2021.pdf`, open access, downloaded 2026-09-07 — and it is **the same paper the
`..._Steinhagen2021` hold-out bundle already cites**, DOI `10.3390/ph14111121`, Pharmaceuticals 2021,
14, 1121. The bundle names an author who is not on it (the authors are Leitzen, Vogel, Steffens, Zapf,
Müller and Brandl; "Steinhagen" appears nowhere in the paper) and its buffer note says "SOURCE PAPER
NOT ON DISK".

The paper's Table 4 measures **six species the trunk carries** in that exact pot — 10 % (w/v) glucose
in water, 121 °C, 18 min, amine-free. The bundle scores **one** of them.

## 2. What the other five say, measured before this wave is written

At the shipped vector, nothing altered:

| species | measured µg/mL | model µg/mL | fold |
|---|---:|---:|---:|
| 3-deoxyglucosone | 52.2 ± 4.0 | 46.98 | **1.11×** |
| methylglyoxal | 2.6 ± 0.2 | 2.03 | **1.28×** |
| 5-HMF | 17.4 ± 3.9 | 1.459 | 11.9× |
| 3,4-dideoxyglucosone | 55.5 ± 1.7 | 1.711 | **32.4×** |
| glyoxal | 5.6 ± 1.3 | 0.161 | 34.8× |
| glucosone | 7.5 ± 1.4 | 0.119 | 63.0× |

**B33 IS REFUTED.** The amine-free entry to 3-deoxyglucosone is right to 11 %, and the route to
methylglyoxal to 28 %. Scaling the entries ×10, as B33 forecast, takes 3-DG from 1.11× to 8.89× and
methylglyoxal from 1.28× to 7.76× in order to buy HMF — it destroys two observables I had no data for
in order to fix the one I could see. The probe behind B33 was not wrong about what it measured; it was
blind, and this paper is what it was blind to.

## 3. What is built

**No constant moves and nothing is fitted.** These are end-of-cook levels in a declared hold-out
bundle: under the standing rule they VALIDATE and may never fit. The wave is three changes to a
benchmark file and its record.

1. **The citation is corrected** to Leitzen et al., and the "SOURCE PAPER NOT ON DISK" note is
   replaced with the path to the PDF that has been on disk since 7 September. The `benchmark_id` is
   NOT renamed: it is cited by artifacts, tests and figures across the branch, and a wrong ID that
   says so in its own provenance is safer than a silent re-key. The mistaken name is kept in the file
   as the audit record.
2. **Five hold-out targets are added** from Table 4, scheme A, the 121 °C row, each with the printed
   value, its SD and a verbatim quote, in the µg/mL → ppb rescaling the existing HMF row already uses.
3. **The record is written**: B33 retracted in place, and the deficit localised.

## 4. What counts as success, declared before the run

- **T1 the existing row does not move.** HMF stays at 17400 ppb and 11.9×, bit for bit.
- **T2 the five new rows score at the folds in section 2**, within rounding.
- **T3 nothing else on the panel moves.** No other bundle's predictions change at all.
- **T4 the headline is reported in both directions.** Adding five rows of which two are inside 3×
  RAISES the pass rate by arithmetic (7/39 → 9/44) and WORSENS the median fold (three of the five are
  32× to 63×). Both must be stated where the rate is stated, or the wave is not honest.

Ship rule: **SHIP if T1, T2 and T3 hold.** T4 is a reporting obligation, not a test.

## 5. Predictions, before the run

1. T1 and T3 hold. **97 %.** Adding declared targets to one bundle cannot touch another, and the HMF
   row is untouched — but this is exactly the assumption B31's T1 was written to catch, so it is tested.
2. The panel's median fold error rises above 12× (from 9.31×). **80 %.** Three of five new rows are
   32×, 35× and 63×, against 39 existing rows with a 9.31× median.
3. **The deficit is `k_tdg_ddg`, not the entry.** 3-DG is right to 11 % and the very next species,
   3,4-dideoxyglucosone, is 32× low, so the 3-DG limb's rate-determining step is where HMF is lost.
   That constant is one of the four ENV-B13 declined to band **because two laboratories agree on it to
   1.5×** — and both determinations are in DRY matrices (a glucose melt and a hazelnut). Agreement is
   not accuracy when both laboratories share a matrix. **75 %** that an aqueous determination of
   3-DG → 3,4-DGE lands more than 10× above the shipped value.
4. Glucosone at 63× low is a SECOND and separate finding: in an amine-free pot the only route to it is
   `k_glc_g`, which B32 measured carrying no flux at all (k_ref 5.34e-8). **This is now a gap with a
   measured size rather than an inference.** 90 % that it survives review as stated.

---

# Outcome (2026-09-11)

## Ship: T1, T2 and T3 all hold

| test | result | pass |
|---|---|---|
| T1 the existing row does not move | 5-HMF 17400 ppb, 11.93×, bit for bit | **yes** |
| T2 the five new rows score as measured | 1.11 / 1.28 / 32.43 / 34.83 / 62.99 | **yes** |
| T3 nothing else moves | **not one pre-existing row's prediction changed**, 39 compared | **yes** |

**SHIPPED.**

## T4, the reporting obligation, in both directions

Panel **7/39 → 9/44** within 3× (out of sample 6/38 → 8/43). The rate rises from 17.9 % to 20.5 %
*by arithmetic*: two of the five rows added are inside the band. And the median fold error **worsens,
9.31× → 10.62×**, because three of the five are 32×, 35× and 63×. Both numbers are in the README next
to each other. Nothing about the model changed in this wave — not one constant moved.

## A unit bug found by the wave, and it is the third time

The first score of the five new rows read **3-deoxyglucosone 180 144×** and **methylglyoxal 92 306×**.
Neither is a chemistry result: 46 980 µg/L was being reported as 0.2898, exactly a factor of the
molar mass. `MOLECULAR_WEIGHT_G_PER_MOL` had no entry for either species, and the reporting loop's
`else` branch **silently returns mmol/L for anything it cannot weigh**.

This is the same failure mode as B28's 2-pentylfuran — where it was diagnosed for a day as a routing
problem and written into five documents before anyone checked — and it was still live three days
later. It is fixed structurally, not locally:

* `TDG`, `MGO` and `ODG` get their molar masses (all three are ordinary molecules the trunk has
  carried since B1/B7/B13; nothing had ever asked for them as targets).
* The `else` is gone. `engine._REPORTED_IN_MMOL_PER_L` **names** the sixteen accounting pools that
  are not molecules — a mole of "melanoidin nitrogen" is a mole of atoms — and anything else with no
  molar mass now **raises**, with a message that says what to do.
* A test walks every name reachable through `TARGET_ALIASES` and asserts it can be weighed or is a
  declared pool. That test would have failed on the day B28 shipped.

## What the six observables say about the trunk

Two findings, on one pot, with no fit anywhere.

**1. The HMF deficit is one step, and it is not the entry.** 3-deoxyglucosone lands at **1.11×** and
methylglyoxal at **1.28×** — the amine-free entries are right. The very next species, 3,4-dideoxy-
glucosone, is **32× low**, and 5-HMF downstream of it is 11.9× low. So the loss is in
`k_tdg_ddg`, the step this repository's own note calls "THE RATE-DETERMINING STEP OF THE 3-DG LIMB".

That constant is one of the four ENV-B13 declined to band **because two laboratories agree on it to
1.5×** — and both determinations are in DRY matrices, a freeze-dried glucose melt and a hazelnut.
**Agreement is not accuracy when both laboratories share a matrix**, and this is the first aqueous
measurement to reach it. Prediction 3 stands at 75 %: an aqueous determination will land more than
10× above the shipped value.

**2. Glucosone is 63× low and glyoxal 35× low, and that is a structural absence with a measured
size.** In an amine-free pot the only route to glucosone is `k_glc_g`, which B32 measured carrying no
flux at all (k_ref 5.34e-8 — six decades of scaling changed nothing). Until today that was an
inference from a rate constant; it is now a number: the pot makes 7.5 µg/mL and the model makes 0.12.

## Predictions, scored

1. T1 and T3 hold (97 %) — **right**, and worth having tested.
2. The median rises above 12× (80 %) — **wrong**. It rose, to 10.62×, not past 12×.
3. The deficit is `k_tdg_ddg`, not the entry (75 %) — **supported**, decisively, by 1.11× against 32×.
4. Glucosone at 63× is a second, separate finding (90 %) — **stands**.

## A fourth finding the envelope made visible, recorded not built

The five new rows went through the Monte-Carlo envelope and **three of them came back with intervals
of about 1e-6 decades** — 3-deoxyglucosone, 3,4-dideoxyglucosone and glucosone are published as
*exact*. That is precisely the defect ENV-B13 was written to fix for the five hydroxymethylfurfural
rows, and it is still true of the 3-DG limb and the amine-free entries, which have no `CorePrior` row
at all.

It matters most on the wave's best row: 3-deoxyglucosone lands at **1.11× on fold error and OUTSIDE
its own interval**, because the interval has no width. A model that is right to 11 % and says it is
certain is making a stronger claim than its evidence supports.

**The fix is available and is not invented**: `k_glc_tdg` already carries a declared 95 % HPD that is
58 % of its estimate (flag `hpd_58_percent_of_estimate`), and `k_tdg_ddg` carries one too. Turning
those into prior rows is an ENV wave of the ENV-B13 shape, with the same method note — measure the
sampler's noise floor from two seeds and compare against the observed maximum. It is named here
rather than bolted onto a wave that was about something else.

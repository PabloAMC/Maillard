# Wave B16 pre-registration — the thiol sinks against a 100 C time series and measured TTCA decay

*Written 2026-09-07, after the owner downloaded Schieberle, Hofmann & Münch 2000 (10.1021/bk-2000-0756.ch010)
and Zhai et al. 2021 (10.1021/acs.jafc.1c03727) and before any B16 number. Base wave: B9. Generator
`generate_kinetic_core_b16_fit.py`; ship rule `generate_kinetic_core_b16_ship_rule.py`; licence
FIT_HOLDOUT_DECLARATION.md Amendment 26.*

## 1. The finding that motivates the wave

Scored on the shipped B9 lane before this wave (SCH-T-01, RIB-T-01/02 in `core_directional_scores.json`):

- Schieberle 2000 Table IV — the sulfur fit's own Hofmann 1998 pot (ribose 100 mM + cysteine 33 mM, 0.5 M
  phosphate pH 5.0) held at 100 C — has MFT rising 4.5 -> 13.8 -> 156 -> 179 ug and FFT 2.0 -> 3.1 -> 110
  -> 132 ug over 30 / 60 / 360 / 720 min (stable isotope dilution). The shipped lane predicts MFT 66 -> 130
  -> 36 -> 16 ug/L: a PEAK at 60 min and a fourfold LOSS by 12 h. **DISAGREE.**
- Liu 2023 (LWT 182:114874), a cysteine-rich ribose pot at 168 C: MFT falls by half between 20 and 60 min,
  FFT flat. The shipped lane predicts both below 1 ng/L by 60 min. **Not evaluable (zero).**

Both say the same thing: the lane's thiol removal is far too strong away from the 121-145 C window its
sink constants were fitted in. `Ea_decay_thiol_sink` sits on its 102 kJ/mol ceiling in every wave since
B8 (the slice profile: "wants to go UP"); the B12/B13 wishlist names `k_thiol_decay`, `k_dimer_decay`
and this barrier as the coordinates no primary evidence identifies. Schieberle's series is the first
primary evidence on the sinks' temperature dependence in the fit's own reference frame.

## 2. What is FIT (the owner's rule: within-study ratios and measured rates fit; levels validate)

Four new systems: the Hofmann pentose pot at 100 C for 30 / 60 / 360 / 720 min (same initial charge,
buffer and pH as `hofmann_pentose_pH5`). Seven cross-system ratio rows, sigma_log 0.10 (the source
prints no SD; 0.10 is the B10 fold sigma):

| row | target | note |
|---|---:|---|
| MFT 60 / 30 min | 3.07 | |
| MFT 360 / 30 | 34.7 | |
| MFT 720 / 30 | 39.8 | still rising at 12 h |
| FFT 60 / 30 | 1.55 | |
| FFT 360 / 30 | 55.0 | |
| FFT 720 / 30 | 66.0 | |
| MFT 145 C-20 min / 100 C-360 min | 1/13 = 0.077 | the text's "factor of 13"; links the series to the objective's 145 C level |

Three new systems for Zhai 2021: TTCA 10 mmol/L in water, pH 7.0 (NaOH), no buffer, 100 / 120 / 140 C,
60 min. Three concentration rows on TTCA remaining (from the zero-order fits c = c0 - k t, sigma_log
0.05): 8.71 / 6.07 / 4.46 mmol/L. These are measured rates on a species the lane carries (the
first temperature series on `k_ttca_deg`).

Objective: B9's 54 rows + 10 = 64 rows; free set B9's 23 (no new coordinate). **Two variants, both
run, both reported:** (a) `b16` with every band as B9 (the thiol-sink ceiling 102 kJ/mol, the owner's
2026-09-04 decision, kept); (b) `b16_lift` with the thiol-sink ceiling at 160 kJ/mol, run as
INFORMATION on what the series asks for. Only (a) can ship under this prereg; if (b) is what the data
want, that is reported to the owner as the decision it is.

## 3. What is VALIDATED

- Schieberle levels (ug/L at 100 C): validation, not read.
- SCH-T-01 (the monotone rise through 720 min), RIB-T-01/02 (Liu 2023 at 168 C), the Yiltirak ladder
  (100-130 C levels), Bolton, the Hofmann Table-1 pH-5 rows, Kang 140 C: scored after the fit is frozen.

## 4. Pre-registered tests and the ship rule

- **T1 shape:** SCH-T-01 AGREES under the frozen B16 parameters, and all seven fold rows land within
  0.3 dex. Falsifier: MFT still peaks before 720 min.
- **T2 in-sample discipline:** no B9 row's |residual| grows by more than 0.3 dex (the 121-145 C window
  must not be traded for the 100 C series). Falsifier: any row over 0.5 dex.
- **T3 Yiltirak 100 C / 4 h and 110 C / 2 h levels:** fold errors not worse than B9's (today MFT 9x /
  FFT 480x at 100 C).
- **T4 TTCA rows:** all three within 0.1 dex.
- **T5 Laplace:** `Ea_decay_thiol_sink` identified (finite sigma) OR off its bound.
- **T6 (recorded, not gating):** RIB-T-01/02 become evaluable (non-zero at 168 C) and their status.

**Ship rule.** B16 (a) ships as the engine's sulfur report if T1, T2 and T4 hold. If T1 fails with the
ceiling kept and passes only in (b), the wave does NOT ship and the finding is put to the owner: the
100 C series cannot be reproduced under the Gigl 2021 ceiling.

## 5. Forecasts

P(T1 passes with the ceiling kept) 0.25. P(T1 passes with the ceiling lifted) 0.60. P(T2 holds) 0.6.
P(B16 (a) ships) 0.20. P(RIB-T-01 becomes evaluable) 0.5.

## 6. OUTCOME (2026-09-07, after the run) — DO NOT SHIP, and a diagnosis

Two starts per variant (`kinetic_core_b16_members`, `kinetic_core_b16_lift_members`), reports, Laplace on
`b16`, `kinetic_core_b16_ship_rule.md`.

| | ceiling kept (102) | ceiling lifted (160) |
|---|---|---|
| cost (64 rows; B9 alone 18.7) | 931 (budget-exhausted, both starts) | 802 (converged, both starts) |
| thiol-sink Ea | 102.0 (on the ceiling) | 160.0 (on the NEW ceiling) |
| T1: MFT at 100 C, 30/60/360/720 min (measured 45 / 138 / 1560 / 1790 ug/L) | 22 / 55 / 109 / 96 -- peaks at 6 h | 24 / 75 / 358 / 329 -- peaks at 6 h |
| T1 worst fold row | 1.50 dex | (fails) |
| T2 in-sample worst | +1.90 dex (`fed_ribose_h2s_MFT`) | +0.57 dex |
| T3 Yiltirak 100 / 110 C levels, median fold | 198x -> **14x** | (not gating) |
| T4 TTCA rows | 1.25 dex (clipped) | fail |
| T6 Liu 2023 at 168 C | evaluable; MFT 12.2 / 8.5 / 5.9 ug/L DECREASING (agrees); FFT 18.6 -> 8.0 (not flat) | same shape |
| Laplace (b16) | 15 of 23 identified; thiol-sink barrier on its bound | — |

**Ruling (sec. 4):** T1, T2 and T4 all fail on `b16`; T1 fails on `b16_lift` as well, so the wave does not
ship and there is no ceiling decision to put to the owner in the form the prereg anticipated -- lifting
the ceiling is NOT enough. The engine keeps reading B9.

**What the run says.**
1. The 100 C series cannot be reproduced by moving the thiol-sink barrier alone. At 160 kJ/mol the sink
   at 100 C is 2.4 dex slower than at 145 C and MFT still turns over at 6 h: what stops the rise is not
   only the sink but the FORMATION drying up -- the pot's pentose / intermediate supply is spent in the
   model by 6 h, where the real pot keeps producing for 12 h. The formation side (the route's barrier and
   the pentose-consuming competition) is as much the problem as the sink.
2. Both variants buy the 100 C rows with the fed rows (+1.9 / +0.6 dex on the ribose + H2S pot): a single
   set of sink constants cannot serve 100 C and the 145 C fed pots. The structure needs a sink whose
   effective order or temperature dependence differs from a first-order Arrhenius step -- e.g. a
   reversible dimerisation (the disulfide pool re-releasing thiol) or a sink that scales with the
   melanoidin / carbonyl pool rather than with time.
3. **The between-lab misses are a low-temperature sink problem.** Weakening the sinks at 100-110 C moved
   the Yiltirak 100 C / 110 C median from 198x to 14x without any lab-specific term. Since B10 and B11
   this is the first change that closes most of that gap; it fails only because it breaks the 145 C rows.
4. The TTCA rows fail because the core CONSUMES TTCA 18x faster than Zhai measured at 120-140 C (the
   residual is clipped). Probe on the shipped B9 lane, TTCA 10 mM at pH 7 after 60 min: 1.58 mM left at
   100 C (measured 8.71), 0.05 at 120 C (measured 6.07), 0.00 at 140 C (measured 4.46), with free cysteine
   at 2.2-2.9 mM and the pentose consumed: the RING-OPENING step `r_ttca_cys` (TTCA -> Cys + pentose,
   fitted on Kang 2026's free-cysteine readings) runs about ten times too fast, and the released pentose
   is then eaten by the sugar trunk. The fitted k_ttca_deg is not the issue; k_ttca_cys is.
5. Liu 2023's 168 C decline (RIB-T-01) becomes evaluable and AGREES under both variants; FFT is not flat.

**Forecast scoring.** P(T1 with ceiling kept) 0.25 -> failed; P(T1 lifted) 0.60 -> failed (the surprise);
P(ships) 0.20 -> did not.

**Next (plan W7):** a sink-structure wave -- reversible thiol dimerisation with a temperature-dependent
equilibrium, and the pentose supply at 100 C -- pre-registered against the same ten rows plus the fed
pots. Until then the panel keeps SCH-T-01 as a miss and the Yiltirak ladder as the 14x-vs-198x lever.

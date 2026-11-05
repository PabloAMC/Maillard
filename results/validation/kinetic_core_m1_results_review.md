# M1 (EX-061): independent results review

- Date: 2026-11-05. Reviewer: an independent agent with no part in the design, the build or the run. Read-only;
  every number recomputed from the fit, carried and comparison records (scripts in the session scratchpad).
- Verdict on the run: no blocking problem. The frozen rule (prereg v3.3) was applied exactly; all 24 unit x reading
  values reproduce to 0.01; tau 0.964 is H0's tau_ABS(all); both rivals "compared" (R0 starts within 0.00065 nats,
  C1 within 0.0124); amendment A1 applied correctly (C1's "compared" depends on it; benign: the clean start 0 agrees
  within 0.012 nats).

## The verdict "not separated / conflict" stands, but it is artefact-driven

The tie rests on one number: R0's linearised RES_Cys lppd under the alternative floor, -175.4. R0's cysteine moves
with one identified coordinate (`k_cys_thermal`); the exact predictive is a 1-D integral. Scanned on the engine
(delta -1.2 ... +1.2 dex, Laplace sd 0.279), the exact R0 value is about -746, and d_Cys(alt) is about +105 (against
C1's MC) to +250 (against C1's linearised). Evaluated exactly, every reservoir setting favours C1 and the frozen rule
would read "C1 favoured" under all six readings. The rule cannot be changed after the fact: the verdict stands as
"not separated", recorded here as driven by a linearisation artefact, not by the data.

## The measurement model needs three fixes before M2 (MAJOR)

1. **Floors at measured values zero a rival's width.** The primary floors were data values (MFT 45 ug/L and FFT
   20 ug/L = the 30-min measurements; Cys 0.80 and PENT 1.22 mM = the 50-min measurements). R0's MFT (3/3 points),
   Cys (8/8) and PENT (8/8) fall below them, so R0 became a flat line at a measured value with zero width; the
   alternative floor (1 ug/L; 0.22 mM) restores R0's own prediction. The whole flip is on R0's side; C1's reservoir
   is never floored. The stable numbers are the plug-in ones: Cys +375 / +205, PENT +229 / +230, MFT +261 / +236
   (isothermal), FFT -54 / -164 (primary / alternative floor).
2. **H0's between-laboratory tau was applied to a C/C0 series.** The reservoir rows plot ln(C/C0) with C0 the 40 mM
   charge (`route_submodel.yml` value_basis): a calibration offset cancels. tau ~1 dex hid R0's main failure; at tau
   0.03 every reservoir number favours C1 (primary Cys +889 / +505, PENT +404 / +223; alternative Cys +850 / +653,
   PENT +1102 / +921). Mandated by the prereg, so not a misapplication; a design error.
3. **Linearised densities on floor sensitivities.** The MC check ran only on the primary setting; on C1 the MC SE
   reached 1.00 (effective sample size ~1), so it said nothing where it mattered.

## What the data show

| species | R0 (y - f, dex, 5 -> 50 min or 30 / 360 / 720 min) | C1 |
|---|---|---|
| Cys | +1.68 ... +0.38 (99 % gone before 5 min) | -0.21 ... -1.34 (decay ~5x too slow: 0.014 vs 0.072 /min) |
| PENT | +3.7 ... +15, then the integrator floor | -0.24 ... -1.08 (~3.4x too slow: 0.018 vs 0.060 /min) |
| MFT | +0.68 / +1.73 / +2.07 | +2.04 / +1.08 / +1.02 |
| FFT | -0.05 / -0.16 / +0.05 (essentially exact) | +1.94 / +0.64 / +0.20 (2 dex low at 30 min) |

- Supply is real: free cysteine and pentose persist on a ~10 - 15 min timescale at 120 C, where R0 loses ~99 % of
  both within 5 min. C1's reservoir is right in kind and wrong in rate (3 - 5x too slow).
- The MFT 360/30 ratio: R0 0.48, C1 2.50, measured 1.54: both miss by ~1 dex in opposite directions.
- FFT favours R0 for real. C1 is 36 nats worse than R0 on B60's own 52 rows (58.0 vs 21.9).
- Compensating structure: the rivals use shared coordinates very differently (`k_fur_decay` +0.18 vs -9.98;
  `k_osone_decay` -1.35 vs -9.98; `k_arp_dpo` -10 vs -2.75). C1's `k_fur_decay` and `k_osone_decay` sit 0.02 dex above
  their -10 bound: in effect switched off.
- Bound-active coordinates: R0's `Ea_decay_thiol_sink` profile is flat (Delta cost 0.04 - 0.59 nats); C1's binds (profile
  sd 2.3 kJ/mol), but relaxing it to 112 / 122 kJ/mol keeps the D1 signs: "conditional on the Gigl ceiling" applies
  to C1's D1 numbers, not to the sign pattern or the verdict. C1's `ea_form` binds at 40 kJ/mol (the route's
  temperature law is under strain, Cao vs Zhai); the reservoir prediction is insensitive to it.

## Summary for the owner

M1 ran as frozen and the verdict "not separated" follows the rule, but the tie rests on one linearisation artefact:
evaluated exactly, C1 wins the reservoir block. What the data show: supply through a reservoir is real (R0 has none and
loses ~99 % of cysteine and pentose in 5 min), but C1's reservoir is 3 - 5x too slow, it costs 36 nats on B60's rows,
and it loses the FFT shape that R0 fits almost exactly. M1 did not establish an adequate supply model, did not
attribute the D1 result to supply, and did not decide the loss law. Before M2: no between-laboratory tau on C/C0 data,
floors that are not data points, exact or MC densities on every floor sensitivity, and a reservoir whose release rate
is free to be 3 - 5x faster than C1's.

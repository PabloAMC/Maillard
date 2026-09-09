# Wave B14 pre-registration — water activity on the acrylamide lane, declared flat inside its measured window

*Written 2026-09-07, after the owner downloaded De Vleeschouwer 2008 (10.1021/jf8006294) and before
the directional panel was re-scored. Module: `src/kinetic_core/acrylamide_conditions.py`; licence:
FIT_HOLDOUT_DECLARATION.md Amendment 24; dossier `devleeschouwer2008_extraction.md`.*

## 1. What the paper holds for the lane

The lane's shipped constants (`k_int1_acr` 3.57e-3 /min, Ea 159.2; `k_asn_glc` 1.70 /M/min, Ea 117.5;
both at a_w 0.92) are the a_w 0.92 column of this paper's Table 2. The same constants are fitted at
a_w 0.88 / 0.96 / 0.99 and the authors find, at 95 % HPD, no significant change in formation or
elimination; the point estimates of k_Fref still span 1.45-3.57 x 1e-3 /min (0.41-1.0 of the shipped
value). Nothing below a_w 0.88 is measured.

## 2. What changes (declared, no fit)

- Inside a_w 0.88-0.99 the acrylamide-forming step carries a multiplier of exactly 1.0 (flat) with an
  envelope band (0.41, 1.39) sampled uniformly (`acrylamide.aw_multiplier`).
- Outside the window: no term (the value is recorded, changes no rate, and the run says so); a
  comparison whose arms straddle the boundary is REFUSED with the window named.
- At a_w None the term is inert. Every existing panel row, claim and fit report reproduces exactly.

## 3. Pre-registered expectations on the directional panel

- AW-05 (new; `fit_adjacent`: declared from the same finding): evaluable and AGREE (flat by
  construction). It does not enter the independent headline.
- AW-02 (extrusion, a_w 0.3 / 0.6 / 0.9): stays NOT EVALUABLE (two arms below the window).
- Every other claim unchanged; headline 18/30 unchanged; all-claims 20/42 -> 21/43.

## 4. What would falsify the declared term

A measurement of acrylamide in the same system at two water activities inside 0.88-0.99 differing by
more than the band (a factor above 1.39 or below 0.41). Table 3 (potato matrix) is the nearest
existing check: k_Fref 0.83-2.51 x 1e-3 across the window, again overlapping at 95 % HPD.

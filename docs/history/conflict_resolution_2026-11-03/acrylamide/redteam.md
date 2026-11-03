# Red team: acrylamide cluster (CF-039, CF-104), 2026-11-03

This review is read-only. I ran the acrylamide lane (`core_parameters("acrylamide")` + `integrate_acrylamide`, no
envelope) through `./scripts/docker_maillard.sh run "nice -n 10 python ..."`. My probes are in the session scratchpad
under `redteam_acr/`:
- `probe_ea.py`: the lane's own emergent Ea on Gökmen's estimator, and the output ratios at the proposed band edges.
- `probe2.py`: the printed-CI edges, and the emergent Ea with k_int1_acr at 95.

The page reads were: Gökmen 2012 p. 169 (Methods) and Fig. 2b, cropped at 300 dpi; and the De Vleeschouwer
2008 / 2008b / 2009b text layers. The code reads were: `parameters_acrylamide.py` (T_REF_A_K, MEASURED_ACRYLAMIDE),
`acrylamide_conditions.py`, `uncertainty.py` (`_b3_priors`, `_declared_band_priors`), the B3 fit report,
`core_panel_scores.md` and `core_directional_scores.json`.

## Claim 1 (Proposal A): a 95-170 kJ/mol formation-Ea band below 140 C, anchored at 160 C

**The anchor: UPHELD (P 0.95).** `T_REF_A_K = 433.15` (160 C) is the lane's T_ref. Both measured constants are printed
there: k_int1_acr has 3.57e-3 /min with Ea 159.2 ± 29.5, and k_asn_glc has Ea 117.5 ± 25.2. Neither measured Ea is
drawn by the envelope. `_b3_priors` draws only B3's identified k_acr_dp / Ea_acr_dp. **The printed CIs of the two
formation constants are therefore unused today.** The proposal does not notice this.

**The widening factors: WEAKENED (P 0.9 that they are overstated).** The values ×2.7 / ×7.9 / ×26 are the ratio of
the band's two edges (170 vs 95), and even that is slightly mis-computed: exactly, it is 2.74 / 8.32 / 28.5. Relative
to the shipped centre (159.2), the upward widening is smaller. I confirmed it with the engine on acrylamide output
(20 min, 10 mM, k_int1_acr at 95):

| temperature | 100 C | 120 C | 130 C | 140 C |
|---|---|---|---|---|
| upward widening | ×17.6 | ×6.1 | ×3.8 | ×2.4 |

Above T_ref the response is not a monotone "narrowing": at 180 C the output is ×0.58, and at 200 C it is ×2.16.

**"No current score row moves (all 140-200 C)": REFUTED (P 0.97).** Several rows and claims sit below 140 C:
- `acrylamide_spi_extrusion_130C_ACSRef3` is a PRIMARY panel row at 130 C, and it is also an in-sample B3 fit row.
  Its envelope is [0.0022, 0.109] ppb against an observed 150 ppb.
- The directional claims TEMP-01 (an arm at 130 C, verdict disagree), TEMP-02 (130 C vs 25 C) and ACR-02 (an arm at
  130 C) all have arms below 140 C.

The intervals of these would widen by about ×3.8 upward. No verdict flips: 0.109 × 3.8 = 0.4 ppb, which is still far
below 150. The sentence is false all the same.

**Selective application: REFUTED as designed (P 0.85).** Applying the band "below 140 C" while anchoring it at 160 C
creates a ×2.4 discontinuity at 140 C. It also widens only the region where the lane under-predicts. The SPI row is
×4,000-9,600 low, and TEMP-01 has the model rising where the data fall. Restricting a band to where the misfit lives
is exactly what the misfit protocol forbids. The evidence (Gökmen, Knol, the DV 2008 k_F columns 84.6 / 89.9) is
temperature-agnostic, so the band must be too.

**The upper edge 170 truncates the measured CI.** DV 2009's 95 % upper bound is 188.7. Class-A fits reach 182.2
(DV 2007), and Claeys Cys reaches 206.

**The lower edge 95: UPHELD, with caveats (P 0.6).** The 95 is the depletion-corrected value, so it does not
double-count the yield-per-time artefact. It is a composite (lumped) estimand placed on one step. That is a category
mismatch, but in effect it is conservative: the lane's emergent early-time Ea is about the sum of its sequential
steps.

**New finding: the like-for-like comparison is far worse than the band.** I applied Gökmen's estimator to the
lane's own output in a Gökmen-like pot (100 mM Asn = Glc, 60-min holds at 90 / 120 / 150 C, and 10 min at 180 C):

| lane setting | Ea 90-120 C | Ea 120-150 C | Ea 150-180 C | 60-min yield vs Gökmen, 90 C | 60-min yield vs Gökmen, 120 C |
|---|---|---|---|---|---|
| shipped | 296 | 154 | 96 | ×2000 low | ×15 low |
| k_int1_acr at 95 | 232 | 89 | — | ×60 low | ×2.4 low |

Gökmen's own 90-120 C value on this estimator is 103 (87-131). The B3 report already carries the same miss:
`knol2005_EaF` has observed 94.4 against predicted 189.5. **Below 120 C the lane is unanchored by orders of
magnitude, and no Ea band covers that.** Below 120 C, refusing is better than a band (P 0.8). Note also that 120 C is
the measured constants' own `temperature_range_c` floor; 140 C is not a floor anywhere in the code.

## Claim 2 (the CF-039 diagnosis): Gökmen's k is yield / hold time

**UPHELD (P 0.85).** The Methods sample 5, 10, 20, 30 and 60 min at 90, 120, 150 and 180 C. The Fig. 2b caption says
60 min at every temperature. The 10-min hold at 180 C is the cluster's inference, but it is the natural one: Fig. 2a's
180 C curve is flat after 10 min, so the 0-10 min slope is yield / 10.

The same zero-order reading follows from any initial-slope estimator, so the diagnosis does not depend on the exact
hold. The 90 C bar reads 2-3 px at 300 dpi, i.e. 3.0e-6 ± 1.5e-6 mmol (1 px ≈ 1.5e-6 mmol). The cluster's
×0.5-1.5 read error is right.

## Claim 3 (Proposal B): flat AW_MULTIPLIER, a moisture window, and the "residual 1.01"

**Keeping AW_MULTIPLIER = 1.0: UPHELD (P 0.75).**

**The De Vleeschouwer basis: the cluster's P(wet) 0.6 should RISE to about 0.8.** DV 2009b (p. 120) says that
"due to differences in initial moisture content" the Asn–sucrose system has the lower initial Asn concentration.
Moisture can dilute a concentration only on a whole-sample (or per-water) basis, and per-water is already excluded at
7.6 M. The confound is sucrose's molar mass (342 vs 180 g/mol), which would also lower Asn per dry kg. DV 2008b
(p. 6 col. 2) contrasts the "much higher initial precursor concentrations" of the dry systems, which is consistent
with a whole-sample basis.

**"Mass action predicts 0.904 → 0.430 (residual 1.01)": WEAKENED (P 0.5).** I reproduce the likelihood ratio of
about 11. However, Bassama's k_F and k_E come from one acrylamide-only fit, in which they are strongly co-determined
(peak height ≈ k_F/k_E, peak time ≈ 1/(k_F + k_E)). The same step's k_E rises ×2.4 with no explanation. Calling the
formation half "explained" and the elimination half "unexplained" from one correlated pair is not consistent. The
ratio CI also ignores the covariance.

**Supporting evidence the cluster missed.** DV 2008's own wettest column, a_w 0.99 (34 % water), already has
k_F ×0.41 relative to 0.92. The HPDs separate: 1.03-1.87 vs 2.26-4.64. This is the source of AW_SCALE_BAND's 0.41.
Formation slowing at the wet edge is therefore seen in two laboratories, which supports a lower band edge, not a
centre change.

**The options as written: REFUTED for "flag out-of-window" (P 0.9).** The engine has no moisture input; only
`process.water_activity` exists. Plantain at a_w 0.972 (70 % water) and DV's powders at a_w 0.96-0.99 (≤ 34 % water)
are indistinguishable on a_w, so no a_w proxy can draw the line.

Worse, out of window `aw_multiplier()` returns 1.0 and ignores the draw. "Flag out-of-window" through `in_window`
would remove the 0.41-1.39 band, which narrows the high-moisture envelope: the opposite of the intent. The same trap
already applies to a_w 1.0 pots, for example PH-ACR-01.

**Plant-based matrices.** High-moisture extrudates and patties (55-70 % water) have a_w of about 0.97-0.99, so today
they are in window with the 0.41-1.39 band. Acrylamide in a patty forms in the crust (< 20 % water), not the bulk, so
a bulk-moisture flag would mis-flag them.

No shipped row declares moisture. The SPI row declares a_w 0.35, TEMP-01 0.5 and ACR-02 0.4. **Under either option,
no shipped number moves today.**

## Claim 4: does either proposal tune a band to absorb a misfit?

- **A as written: YES in form (P 0.7).** Its temperature restriction coincides with the misfit region, and its
  "no row moves" premise is false. Made temperature-agnostic, it is a second-laboratory band in the ENV-B13 style,
  which is allowed.
- **B: NO (P 0.75).** The 0.28 edge comes from an external laboratory and is corroborated by DV's own 0.99 column. It
  does not touch any scoring row.

## Verdicts

**Proposal A: APPLY-CHANGED.**
1. Refuse acrylamide-lane answers for 60 C < T < 120 C. Refusal is at `declare_envelope` / `_acrylamide_warnings`,
   using `temperature_range_c` (120-200 C, already declared). Keep the current near-zero at ≤ 60 C, so TEMP-02's
   25 C control stays evaluable.
2. Do not use the "below 140 C" band. Instead, in `uncertainty._declared_band_priors`, add
   `acrylamide.k_int1_acr.ea_kj_mol`:
   - uniform over [95, 189] kJ/mol: 95 is Gökmen depletion-corrected, Knol k4 and the DV 2008 outliers; 189 is the
     DV 2009 CI upper bound;
   - k(160 C) held at 3.57e-3 /min;
   - applied at all temperatures, with reason `SECOND_LABORATORY_BAND`.

   Plumb it through `acrylamide_draw["ea_int1_acr"]` with `dataclasses.replace(parameters["k_int1_acr"], ea_kj_mol=...)`
   in engine.py before `acrylamide_conditions.apply`.
3. Recommended in the same change: draw k_asn_glc's Ea as normal(117.5, 25.2/1.96).

The effects:
- No central value moves.
- Envelopes widen at 130 C by up to ×3.8 upward and ×0.54 downward. At 180 C they widen to ×0.58-1.18, and at 200 C
  to ×0.10-2.16. This moves every acrylamide envelope, including the Chang 2021 and Lin 2022 rows at 180 C.
- Arms in 60-120 C become refusals. None exist on the panel today.

Because shipped intervals move, ask the owner first.

**Proposal B: APPLY-CHANGED.**
1. Keep AW_MULTIPLIER = 1.0 and the band (0.41, 1.39).
2. Add an optional `moisture_wb` to the process spec. When it is > 0.35 (DV's calibrated maximum, 34.2 %), draw the
   formation band (0.28, 1.39), and do not route this through `in_window`. Warn that the a_w window is calibrated on
   ≤ 34 % water, and that bulk moisture is not crust moisture.
3. Add the whole-sample-basis declaration, as the cluster proposed.
4. Separately, fix the trap where out-of-window a_w drops the band.

Effect: no shipped number moves (no row passes moisture). For high-moisture plant-based inputs that declare moisture,
the envelope's low edge falls ×0.68.

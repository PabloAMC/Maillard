# Pre-registration: B55 (an activation-energy prior for the single-temperature sulfur constants) and B56 (a per-lane structural-discrepancy term) — written 2026-10-28, before any coverage number under either change existed

Two envelope waves on one branch, each behind its own switch, each with its own ship rule, so that
whatever moves is attributable to one of them. Neither moves a fitted centre, a point prediction or
any fit report. Both change only what `src/kinetic_core/uncertainty.py` draws.

## 1. Why

The shipped envelope's nominal 90 % interval covers **19 of 48** literature rows and **17 of 43** out of
sample (`core_prediction_uncertainty.json`, n = 200, seed 0, as committed at `2cab8e5`). Its own
docstring says why: it samples **parameter** uncertainty only, and "what is left is model-structure
error". Two pieces of that are named here.

1. **B55.** 31 of the 43 sulfur-lane rate constants are fitted (or carried) as log10 k at **145 °C** and
   given a barrier nobody measured: 30 share the lumped formation barrier **64.08 kJ/mol**, which the
   B2.3 fit reported **UNIDENTIFIED** (`kinetic_core_b2_3_fit_report.md`) and B9 froze, and one
   (`k_mmft`) has no barrier at all and is held at its 145 °C value (policy 2). The envelope samples
   none of the 31 barriers (`b8.lumped_formation_Ea_kJ_mol`: "frozen in the sulfur fit, source
   uncertainty unrecorded"). So every prediction away from 145 °C treats a temperature slope nobody
   measured as known. This is the green "rate known at one temperature only" class of
   `scripts/generators/build_reaction_tree.py`.
2. **B56.** Even a perfectly propagated parameter posterior is centred on a model that is wrong in
   structure. The residuals the fits themselves leave on the rows they read are the only
   in-sample measurement of that error the corpus has, and the envelope adds none of it.

## 2. B55 — what changes

**Targets (computed, not listed by hand).** Every key of the shipped sulfur vector
(`kinetic_core_b9_fit_report.json` → `log10_k_ref_at_145C`, 43 keys) whose barrier is neither a
measurement (`parameters_sulfur.MEASURED_EA_OVERRIDES`: five keys) nor a fitted decay-family barrier
(`DECAY_FAMILIES`: seven keys, already drawn flat across their capped bands since 2026-09-04). That is
**31 keys**: 14 on the `sugar_trunk` formation route, 16 on `thiol_assembly`, and `k_mmft` (no route;
a thiol condensation, so the sulfur class).

**The prior: the measured activation energies of the same reaction class in this corpus.** Two
classes, keyed by the route table the code already carries (`FORMATION_ROUTE_OF`, B10):

- **carbohydrate class** (target route `sugar_trunk`): every MEASURED barrier in the registries for a
  step whose reactant is a sugar, an Amadori/Schiff product or a sugar-derived carbonyl and carries no
  sulfur. Martins & van Boekel 2005 Table 2 (ten: `k_glc_fru`, `k_fru_glc`, `k_schiff`, `k_ama_tdg`,
  `k_ama_odg`, `k_ama_mgo`, `k_fru_acids`, `k_odg_aa`, `k_tdg_fa`, `k_tdg_mel`); Kocadagli & Gokmen
  2016 Table 2 (seven: `k_fru_int`, `k_fru_odg`, `k_glc_g`, `k_glc_tdg`, `k_int_hmf`, `k_odg_da`,
  `k_tdg_mgo`); De Vleeschouwer 2009 Part I Table 3 (two: `k_asn_glc`, `k_int1_acr`); Zhang 2026 k16
  (Cys-Amadori → α-dicarbonyl). **n = 20**, median 103.7, range 30–237 kJ/mol.
- **sulfur class** (target route `thiol_assembly`, and `k_mmft`): every MEASURED barrier for a step with a
  sulfur-bearing reactant. Zheng & Ho 1994 (cysteine → H2S, glutathione → H2S); Zhang 2026 (k17 thiol →
  disulfide; cysteine → Cys-Amadori); Kang 2026 and Zhai 2021 (free-cysteine depletion, two
  laboratories); Chan & Reineccius 1994 Table I (methional at pH 6 / 7 / 8, 2-acetylthiophene); Gigl 2021
  (thiol covalent capture, central value); Hamzalioglu & Gokmen 2018 (HMF + cysteine); De Vleeschouwer
  2009 (acrylamide + cysteine, cysteine sink, cysteine + glucose). **n = 15**, median 100.9, range
  29.7–133 kJ/mol.

Every value is READ from the registry that already carries it (parameter key, or the named module
constant); the four Chan values, which live only in a source string, are transcribed once and a unit
test asserts the string still prints them. Values whose registry barrier is 0.0 (`k_da_sink`,
`k_go_sink`, `k_hmf_self`: temperature held, not measured) are excluded, and so is anything
`derived_from_fit_data`, `bounded_from_a_timescale_bracket` or `structural_constant`.

**The draw.** Each of the 31 keys gets its own prior row `b55.<key>.ea_kj_mol` and its own stream
(ENV-M1), drawn **independently** from the class's **empirical distribution** — the inverse empirical
CDF with linear interpolation at a uniform quantile, so the draw never leaves the measured range and no
shape is assumed. log10 k at 145 °C stays where the shipped wave (and its Laplace draw) puts it; the
engine receives the drawn barrier through a new override block `single_temperature_ea_kj_mol` applied
after `with_fitted_sulfur`, which refuses a key carrying a measured barrier. Because every sulfur constant
is quoted at T_ref = 418.15 K, **k(145 °C) is unchanged in every draw by construction**.

Declared simplifications: draws are independent across keys and independent of the k(145 °C) draw (the
fit's k was obtained with the barrier fixed at 64.08, so the joint is ignored, which overstates width on
rows near the fit's own temperatures); the prefactor cap (12 decades) is not applied, because the class
range is itself a measured statement about barriers and its extreme (237 vs a 64.08 centre) is
173 kJ/mol ≈ 21.6 decades of prefactor — this is recorded, not hidden. `k_mmft` gaining a barrier in the
envelope does not touch policy 2: the point prediction still holds it at 145 °C; the interval stops
asserting that a slope nobody measured is zero.

Switch: `uncertainty.core_priors(single_temperature_ea_prior=...)`, module default
`B55_SINGLE_TEMPERATURE_EA_PRIOR = False` (the shipped behaviour) during the run.

## 3. B56 — what changes

**The term.** Each Monte-Carlo sample of a row on lane L is multiplied by 10^(σ_L · z), z ~ N(0, 1),
one z per (draw, row) from a stream keyed on (seed, benchmark id, compound) — so the B56 arm and the
B55 + B56 arm use the same z, and the point prediction is untouched.

**The estimator — a declared model-error assumption, not a measurement.** For each lane, from the
residuals of the fits the shipped engine reads for that lane, on the rows those fits read, and never
from a panel row:

  σ_L² = Σ_f c_f · Σ_{i∈f} r_i² / Σ_f n_used,f,   with c_f = n_f / (n_f − p_f)

where r_i is the log10(predicted / observed) residual at the frozen optimum of a two-sided level, yield,
ratio, rate or share row (pH-unit rows, kJ/mol rows and one-sided ceiling/floor rows excluded), n_f the
fit's total objective rows and p_f its free coordinates. c_f is the degrees-of-freedom correction:
residuals at an optimum are shrunk by the p_f parameters the fit spent on them, E[Σ e_i²] = σ²(n − p),
so c_f is exactly the parameter-uncertainty accounting at the fit rows. **The parameter share is
therefore NOT subtracted a second time**: the envelope already samples it, and a fitted model's
in-sample residuals already sit *below* the noise by the leverage term, so "RMS minus the parameter
share" would remove it twice and floor more lanes at zero than the evidence says.

| lane | residual source (shipped, engine-read) | n_f, p_f |
|---|---|---|
| sulfur | `kinetic_core_b9_laplace_covariance.json` `residuals_at_optimum` × each row's `sigma_log` (B2.3/B8 row table) | 54, 23 |
| trunk | B1 Martins 2005 fit (`variant_A.per_response_fit`: per (species, T) group, r² ≈ n_g·(1.4826·log10 median fold error)² — a declared robust approximation, the report prints no row residuals); B7 `fit.residual_log10_fold`; B18 `residual_by_row_dex`; B41 `rows[].after_dex` | 354/8; 3/1; 10/6; 12/5 |
| acrylamide | `kinetic_core_b3_fit_report.json` `rows[]` of kind rate / ppb: log10(fold_error) | 30, 11 |
| lipid | `kinetic_core_b6_fit_report.json` `fit_table`: log10(predicted / measured relative %) | 18, 12 |

The σ_L are written by `scripts/generators/generate_kinetic_core_b56_structural_discrepancy.py` to
`results/validation/kinetic_core_b56_structural_discrepancy.json` BEFORE the B56 envelope arms run, and
the generator asserts that no residual it reads carries a panel benchmark id.

Switch: `propagate_panel(structural_discrepancy=...)`, module default
`B56_STRUCTURAL_DISCREPANCY = False` during the run.

## 4. Arms (n = 200, seed 0, the union panel, the same code)

| arm | B55 | B56 |
|---|---|---|
| A baseline | off | off |
| B | on | off |
| C | off | on |
| D | on | on |

## 5. Ship rules

Coverage is counted on a **fixed denominator**: every out-of-sample literature row with an envelope
formed, a not-evaluable (zero-width) interval counting as a miss, so a term that makes a degenerate
row evaluable cannot shrink the denominator. The artifact's own definition is reported beside it.

**B55** (B vs A). SHIP iff all three hold:
- **T1** — invariance at 145 °C: every row of an isothermal 145 °C pot, and every row on a lane other
  than sulfur, has p5 / p50 / p95 equal to arm A within a relative 1e-9.
- **T2** — widening away from 145 °C: over the sulfur rows at T ≠ 145 °C, the median change in 90 %
  width is ≥ +0.10 dex and no such row narrows by more than 0.05 dex.
- **T3** — no coverage loss: out-of-sample coverage (fixed denominator) in B ≥ A − 1 row.

**B56** (D vs B if B55 ships, else C vs A). SHIP iff
- **T4** — out-of-sample 90 % coverage (fixed denominator) lies in **[70 %, 95 %]**,
- with the median out-of-sample interval width REPORTED (no threshold; a wide interval that covers is a
  weak claim and the width is the price), and per-lane coverage and width reported.

If B55 ships, it ships whether or not B56 does; B56 is judged on top of whatever B55 decided.

## 6. Predictions, before the run

1. B55 T1 passes. **95 %** (an isothermal 145 °C pot sees k_ref exactly; the only failure is a bug).
2. B55 T2 passes. **75 %.** The Yiltirak 100–130 °C rows are 15–45 °C from T_ref; a barrier spread of
   ~±40 kJ/mol is ~±1.3 decades at 100 °C per constant, but many products are ratios of constants on
   the same route, which partly cancels.
3. B55 T3 passes. **80 %.** The class medians (~101–104) sit well above 64.08, so the off-145 medians
   move, direction not claimed; a wider interval rarely loses a hit.
4. B55 adds 0 to 3 out-of-sample hits. **70 %.**
5. σ_sulfur in [0.25, 0.6] dex, σ_acrylamide in [0.4, 0.9], σ_trunk in [0.05, 0.25],
   σ_lipid in [0.1, 0.5]. **60 %** that all four land in their ranges.
6. **B56 T4 passes (out-of-sample coverage ≥ 70 %). 15 %.** The out-of-sample errors on the sulfur
   lane are 1.7–1.9 dex at the median; in-sample residuals of a least-squares fit are, by construction,
   the smallest errors the model makes. Expected out-of-sample coverage under D: 45–60 %.
7. If T4 fails it fails LOW (under-coverage), not high. **95 %.**

What a failure of T4 would mean, said now: the model's structural error is not visible on the rows it
was fitted to, so an honest discrepancy term needs out-of-sample evidence the rule forbids reading — the
next step would be a held-out calibration split declared in `FIT_HOLDOUT_DECLARATION.md`, not a bigger
in-sample σ.

## 7. Verdict

(filled in after the run, below this line, never above it)

### Verdict (2026-10-28; `kinetic_core_b55_b56_ship_rule.{json,md}`, arms in `_b55_b56_arms/`)

The baseline arm reproduces the committed envelope exactly (19/48, 0.9793 dex; out of sample 17/43,
or 17/44 on the fixed denominator).

**B55: SHIP.** T1 pass (40 rows checked — every 145 °C row and every non-sulfur row — bit-identical);
T2 pass (median Δwidth **+0.70 dex** over the ten sulfur rows away from 145 °C, smallest +0.013 at
Resconi's 150 °C); T3 pass (out of sample **17 → 20 of 44**: Yiltirak 100 °C FFT, 110 °C MFT and Bolton
120 °C MFT move inside). Predictions 1–4 held (prediction 4, 0–3 hits, at its edge). The class medians
(~101–104 kJ/mol) sit above the shipped 64.08, so the off-145 medians moved — downward for FFT at
100–120 °C (171 → 61 µg/L at 100 °C against a measured 1.3), i.e. toward the measurement, not claimed in
advance. Shipped: `B55_SINGLE_TEMPERATURE_EA_PRIOR = True`; envelope 22/48 (1.14 dex), 20/43 out of sample.

**B56: DO NOT SHIP.** σ_lane (in sample, dof-corrected): sulfur 0.491, trunk 0.135, acrylamide 1.188,
lipid 0.120 dex. Prediction 5 failed on acrylamide (1.19 above the declared [0.4, 0.9]; B3's rate rows
miss by 0.95 dex RMS in sample). T4 on arm D (B55 + B56): **26/44 = 59 %** against [70 %, 95 %], median
width 2.37 dex (B: 45 %, 1.25 dex). Per lane, D vs B: acrylamide 8/12 vs 4/12 (width 3.98 vs 0.91),
sulfur 9/18 vs 7/18, lipid 7/8 both, trunk 2/6 both. Predictions 6 (15 % pass) and 7 (fails low) held.
The term stays in the code behind `B56_STRUCTURAL_DISCREPANCY = False`. What it says: on the sulfur and
trunk lanes the out-of-sample misses are 1–3 dex and the in-sample residuals 0.1–0.5 dex, so the
structural error is invisible on the rows the model was fitted to; an honest discrepancy term needs a
declared held-out calibration split, not a larger in-sample σ.

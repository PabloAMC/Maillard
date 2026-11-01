# CF-077: the statistical factor on the mixed aminoketone condensation -- PRE-REGISTRATION

**Written 2026-11-01 on `cf077-pyrazine-cross-pair-factor` (off `main` @ `19cf422f`), BEFORE the fix is applied to
the engine and BEFORE any panel artifact is regenerated.** Register entry: `docs/reference/conflicts_register.yml`
CF-077 (`code_mismatch`); found by the pyrazine-rivals red team (EX-052, `docs/history/pyrazine_rivals_2026-10-31/red_team.md`
sec. 2).

What was run before this note was written (disclosed, so the reader knows what the predictions below could lean on):
a scratch probe that multiplied the cross-pair constant by 2 through a monkeypatch of `integrate.rate_constants_at`,
on (a) the bare condensation rates, (b) two glucose + glycine pots through `predict`, and (c) the B18 fit generator's
`predictions()` at the frozen vector on today's trunk. Section 4 quotes them. No panel scorecard, envelope or ranking
was computed under the fix.

## 1. The defect

`src/kinetic_core/network.py` writes the three condensations on one `k_cond`, in plain mass action:

| reaction | stoichiometry | rate today |
|---|---|---|
| `r_akg_pz` | 2 AKG -> PZ | k_cond [AKG]^2 |
| `r_akm_dmp` | 2 AKM -> DMP | k_cond [AKM]^2 |
| `r_ak_mpz` | AKG + AKM -> MPZ | k_cond [AKG][AKM] |

so rate_MPZ = sqrt(rate_PZ rate_DMP) and the mixing ratio R = MPZ^2 / (PZ DMP) = 1 at every instant. The design
(`parameters_pyrazine.py` docstring, the `r_ak_mpz` note, `kinetic_core_b18_fit_report.md` "mixed_route") is random
pairing of the two pools with one encounter constant: with N_A and N_B molecules there are N_A^2 / 2 like pairs and
N_A N_B unlike pairs, so rate_MPZ = 2 sqrt(rate_PZ rate_DMP) and R = 4. The cross pair lacks its factor 2.

Reproduced in the env (probe (a)): R = 1.000 at AKG : AKM = 0.5 : 0.5 and 0.3 : 0.7; 4.000 with the factor.
Through `predict` (probe (b)), glucose + glycine 100 + 100 mM: R = 1.016 (120 C, 60 min, pH 6.8) and 1.021
(95 C, 120 min, pH 9); the excess over 1 is the time heterogeneity of the pools (Cauchy-Schwarz: a time-integrated
R is >= the instantaneous one).

## 2. The fix and the alternative

**Chosen: `r_ak_mpz` runs at 2 k_cond; the two homo pairs keep k_cond.** The alternative, k_cond / 2 on the homo
pairs with the cross pair at k_cond, gives the same partition (only the ratio of the constants matters for R). Why
the first:

1. `k_cond` keeps its meaning: it is the rate constant of the homo pairs' product formation, which is what B18's
   pre-registration, its `k_cond_sensitivity` diagnostic (x0.1, x10 on the Zhou pots) and the ship rule describe.
   Halving it on the homo pairs would change the constant every one of those records refers to.
2. The fitted pots are single-pool (Zhou 2024 fed glyoxal OR methylglyoxal + alanine): they run only homo pairs, so
   with the chosen fix the B18 fit's Zhou rows are bit-identical (no second-order lag change either). Under the
   alternative they would move by the aminoketone lag at k_cond / 2 (tiny, ~1e-3 dex by the x0.1 sensitivity, but
   not zero), and the fit report's frozen predictions would no longer be reproducible.
3. `K_COND` stays one declared parameter with one value; the factor is a structural (combinatorial) constant, not a
   second declared rate, so it is written as a named constant on the reaction, not as a new parameter key.

## 3. What should move, and by how much (fast-condensation limit)

At k_cond = 1e4 L/(mmol min) the condensation is fast (B18 diagnostics: aminoketone ~1e-5 mM), so every
aminoketone formed ends in a pyrazine and the Strecker fluxes F_G (to AKG) and F_M (to AKM) are untouched. Nitrogen
conservation then fixes the TOTAL pyrazine moles; the fix only re-partitions them. With p = F_G / (F_G + F_M):

- after the fix, PZ : MPZ : DMP = p^2 : 2pq : q^2 (binomial, q = 1 - p);
- today, PZ : MPZ : DMP = r^2 : r : 1 with 2 r^2 + (1 - rho) r - 2 rho = 0, rho = p / q.

Predicted ratio new / old per species (computed from these two formulas, not from a run):

| p (AKG share of flux) | PZ | DMP | MPZ |
|---:|---:|---:|---:|
| 0.001 | 0.251 | 1.000 | 1.003 |
| 0.01 | 0.265 | 1.000 | 1.029 |
| 0.1 | 0.380 | 0.980 | 1.221 |
| 0.25 | 0.538 | 0.913 | 1.401 |
| 0.5 | 0.750 | 0.750 | 1.500 |
| 0.75 | 0.913 | 0.538 | 1.401 |
| 0.9 | 0.980 | 0.380 | 1.221 |

So: **MPZ rises by a factor between 1 and 1.5 (largest at balanced pools); the MINOR pool's homo pyrazine falls by up
to 4x; the major pool's homo pyrazine falls by at most 25 %; total pyrazine moles are unchanged.** "Methylpyrazine
2x low" holds for MPZ relative to sqrt(PZ DMP); in absolute terms the shift is at most +0.18 dex on MPZ and between
-0.60 and 0 dex on PZ and DMP. In a single-pool pot nothing moves.

Expected on the shipped surfaces (P = my probability that the statement holds):

- `core_panel_scores`: no scored number moves (its pyrazine rows are refusals on the lipid lane). P 0.95.
- `core_directional_scores`: PH-06 and TEMP-04 (2,5-DMP, trunk) keep `agree`; their values move by <= 0.6 dex.
  P 0.9.
- `tests/unit/test_kinetic_core_b18.py`: the order PZ > MPZ > DMP in the glucose + glycine pot at 120 C, pH 6.8
  holds (probe (b): 298 > 167 > 23 ug/L after, against 316 > 124 > 48 before). P 0.95.
- Every other live artifact that runs a trunk pot with both pools present moves its PZ / DMP / MPZ inside the table's
  envelope; anything that moves outside it, or any non-pyrazine number that moves, is a defect in the fix. P 0.95.
- Sulfur, acrylamide and lipid lanes: bit-identical (the pyrazine steps are trunk-only). P 0.99.

## 4. Is the B18 fit affected?

The B18 fit has six free coordinates and ten rows (`generate_kinetic_core_b18_fit.py`):

- **The two Strecker constants and barriers** are set by Zhou 2024's six rates, each on ONE fed dicarbonyl: one
  aminoketone pool, so no cross pair, so no change. Probe (c) confirms it to every printed digit: the six Zhou
  predictions at the frozen vector are identical with and without the factor (e.g. zhou_PZ_rate_120C 0.16004,
  zhou_DMP_rate_100C 0.0041261 umol/(L min) in both).
- **The two pH slopes** are set by Leahy 1989's four within-study ratios (PZ and MPZ, k(pH 7) / k(pH 9) and
  k(pH 5) / k(pH 9)), simulated as glucose + glycine pots, i.e. BOTH pools. These rows can move. On TODAY's trunk
  (probe (c)) they do: PZ ratios x1.72 / x1.76, MPZ ratios x0.87 / x0.80 at the frozen vector. But today's trunk is
  not the trunk the fit ran on, and the same probe shows the frozen vector already misses the PZ rows by ~1.4 dex
  on today's trunk (residuals 9.6 and 11.3 sigma, against +0.125 and -0.007 dex in the fit report): later waves
  (B21's aqueous glyoxal at least) changed these pots. That drift is pre-existing and is not CF-077's (follow-up).
- The question that decides whether the FROZEN B18 values depend on the bug is therefore the fit-time trunk. There
  the Leahy pots ran almost on aminoacetone alone (B18 ship rule T3: shares PZ 6e-6 %, MPZ 0.023 %, DMP 99.98 %),
  i.e. p ~ 1e-4, the minor-pool limit of the table, where the fix divides PZ by 4 at every pH and leaves MPZ and DMP
  unchanged. A RATIO across pH is then invariant to O(p).

**Prediction:** at the fit commit (`34da5b51`), re-running the generator's `predictions()` at the frozen vector with
the factor changes each Leahy ratio row by < 0.01 dex and each Zhou row by 0, and the cost (2.644) by < 0.01. If so
the frozen B18 values are not conditional on the bug and no refit is needed. P 0.9. **Decision rule, fixed now:** if
any Leahy row moves by > 0.05 dex (a third of its sigma 0.15) or the cost by > 0.25, the frozen slopes are flagged as
conditional on CF-077 and a B18 re-fit is proposed to the owner as its own decision; this PR does not refit.

## 5. What this PR does not touch

The declared one-sided MPZ band and the pH-slope-below-7 band (EX-052 follow-ups) are separate decisions. The frozen
wave records (`kinetic_core_b18_fit_report.*`, `kinetic_core_b18_ship_rule.*`) state what was true at their commit and
are not regenerated (artifact-freshness policy of 2026-10-29); their "2 sqrt" text was the design, which the engine
now meets.

---

## 6. OUTCOME (appended after the runs; sections 1 - 5 are unchanged since commit `55d20af6`)

### 6.1 The B18 fit at its own commit: the prediction FAILED, the decision rule fired

Re-run at `34da5b51` (a detached worktree; the generator's own `predictions()` at the report's frozen vector; the
factor by the same monkeypatch as the probe). Factor 1 reproduces the report to every printed digit.

| row | report | factor 1 | factor 2 | shift (dex) |
|---|---:|---:|---:|---:|
| zhou_PZ_rate_100C / 110C / 120C | 0.029544 / 0.070337 / 0.16004 | same | same | 0 |
| zhou_DMP_rate_100C / 110C / 120C | 0.0041261 / 0.0094551 / 0.020565 | same | same | 0 |
| leahy_PZ_k_ph7_over_ph9 | 0.4989 | 0.4989 | 0.5624 | +0.052 |
| leahy_PZ_k_ph5_over_ph9 | 0.025699 | 0.025699 | 0.083836 | +0.514 |
| leahy_MPZ_k_ph7_over_ph9 | 0.35665 | 0.35665 | 0.37585 | +0.023 |
| leahy_MPZ_k_ph5_over_ph9 | 0.0023924 | 0.0023924 | 0.0042644 | +0.251 |

Cost 2.644 -> 15.364. Why sec. 4's argument failed: it assumed the fast-condensation limit. In the Leahy pots the
aminoketone flux is tiny (2,5-DMP reaches 1.6e-4 mM in 2 h at pH 9 and 3.3e-8 mM at pH 5), and at a quasi-steady
pool the aminoketone lifetime is 1 / sqrt(2 k_cond F): ~6 min at pH 9 but ~400 min at pH 5, longer than the 120-min
window. So the condensation is not fast in the acid arm, the products grow as k [A][B] on accumulating pools, and
the factor reaches the cross pair (MPZ x1.9 at pH 5, x1.07 at pH 9) and the minor pool unevenly across pH. The same caveat applies to sec. 3's table: it is
the fast limit, and in slow (dilute, acid, cool) pots MPZ can rise by up to x2, not x1.5.

By the rule fixed in sec. 4 (any Leahy row > 0.05 dex, cost > 0.25), **the frozen B18 pH slopes are conditional on
CF-077.** The Strecker constants are not (Zhou's rows are bit-identical). This PR does not refit. Sizing, for the owner's
decision only (the generator's `fit_member` at the fit commit with the factor, starts 0 / 1 / 2, all at cost 3.439):

| coordinate | frozen | re-fit with the factor | shift | Laplace sigma |
|---|---:|---:|---:|---:|
| log10_k_go_ak_100C | -6.5415 | -6.5807 | -0.039 | 0.094 |
| log10_k_mgo_ak_100C | -7.5297 | -7.5682 | -0.038 | 0.088 |
| ea_go_ak / ea_mgo_ak (kJ/mol) | 103.1 / 114.9 | same (upper bounds) | 0 | 18 |
| ph_slope_above_7 | 0.1974 | 0.2143 | +0.017 | 0.047 |
| ph_slope_below_7 | 0.5795 | 0.6907 | +0.111 | 0.055 |

The two log10 k shifts are the pH-term transfer from Zhou's pH 8 to the stored pH 6.8 (0.017 x 1 + 0.111 x 0.2 =
0.039), not a change at pH 8. On TODAY's trunk the frozen vector already misses Leahy's PZ rows by ~1.4 dex with or
without the factor (sec. 4), so a B18 re-fit now is a larger question than CF-077 and is the owner's.

### 6.2 The live artifacts after the fix

The artifact-freshness gate re-ran every live artifact it regenerates: `core_panel_scores`, the B60 fit targets, the
experiment-value ranking, the data wishlist, network hypotheses, the model card and the wave manifest all match
(bit-identical). Only `core_directional_scores` moved, in its two trunk 2,5-DMP claims, and both keep `agree`:

| claim | arm | before (ug/L) | after | shift (dex) |
|---|---|---:|---:|---:|
| PH-06 (DMP rises with pH) | pH 4 | 3.107e-4 | 2.581e-4 | -0.081 |
| PH-06 | pH 7 | 6.920 | 6.393 | -0.034 |
| PH-06 | pH 9 | 21.84 | 21.77 | -0.001 |
| TEMP-04 (A > B) | T4 profile, 160 C max | 1.596e-3 | 1.590e-3 | -0.002 |
| TEMP-04 | T1 profile, 145 C max | 4.37e-6 | 4.37e-6 | 0.000 |

The envelope (`core_prediction_uncertainty`) carries pyrazines only as refusals and does not hash `network.py`; it is
not regenerated (35 min) and cannot move. `figure_inputs.json` is stale on `network.py` but its rebuild is deferred
to Phase 3 (owner, 2026-10-30), as the gate reports. Sulfur, acrylamide and lipid lanes: bit-identical (every scorecard
row matches). Front-door pots (the probe of sec. 1, not shipped artifacts), glucose + glycine 100 + 100 mM:

| pot | PZ | DMP | MPZ | R |
|---|---|---|---|---|
| 120 C, 60 min, pH 6.8 | 316.0 -> 297.7 (-0.026) | 47.6 -> 23.2 (-0.311) | 123.6 -> 167.0 (+0.131) | 1.016 -> 4.030 |
| 95 C, 120 min, pH 9 | 2.25 -> 1.28 (-0.246) | 13.88 -> 12.51 (-0.045) | 5.65 -> 8.08 (+0.156) | 1.021 -> 4.082 |

Every pre-registered expectation of sec. 3 held except the B18 one (sec. 6.1).

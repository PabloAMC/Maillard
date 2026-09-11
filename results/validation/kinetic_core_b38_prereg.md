# Pre-registration: wave B38, the identifiability audit (written 2026-09-11, BEFORE any Jacobian was computed)

## 1. Why

The B37 methodology note said the one analysis this repository has not done is an identifiability
analysis: **which of the fitted constants do the fit data actually pin, and which are set by their
declared bands, by collinearity with a neighbour, or by nothing at all.** Pieces exist — the sulfur
fit carries a Laplace covariance and slice profiles, the acrylamide and trunk fits carry Gauss-Newton
intervals — but no artifact reads them together, none names the *directions* the data cannot see,
and none crosses the answer with the envelope's priors. The owner asked for it. This wave does it,
and it moves no constant.

## 2. What is audited

Every fit whose ship rule says SHIP, at the optimum that shipped, using the fit's own residual
vector:

| fit | lane | free coordinates | residual rows | optimum read from |
|---|---|---:|---:|---|
| B8 | sulfur (on top of B2.4-half; B2.3 is its superseded parent and is not audited twice) | 23 free of 48 | 62 | `kinetic_core_b8_laplace_covariance.json` |
| B3 | acrylamide | 11 | 30 | `kinetic_core_b3_fit_report.json` |
| B18 | dicarbonyl → Amadori-ketone sinks | 6 | 10 | `kinetic_core_b18_fit_report.json` |
| B20 | glycation | 5 | 10 | `kinetic_core_b20_fit_report.json` |
| B21 | aqueous glucosone → glyoxal | 2 | 6 | `kinetic_core_b21_fit_report.json` |
| B1 | trunk | its report's own standard errors are read, not recomputed | | `kinetic_core_b1_fit_report.json` |

B8 is also probed on its **25 frozen** coordinates, so a frozen-but-sensitive constant is named.

## 3. Method, fixed before running

At the shipped optimum x*, a central-difference Jacobian J of the sigma-weighted residuals
(step 10⁻³ in log10 k, 0.5 kJ/mol in a barrier, 10⁻³ otherwise); FIM = JᵀJ; χ²_red = 2·cost/dof;
Σ = pinv(FIM)·χ²_red. Per coordinate: **marginal σ** = √Σᵢᵢ (everything else free), **conditional σ**
= √(χ²_red / FIMᵢᵢ) (everything else fixed), and whether x* sits on a declared bound. Sloppy
directions are the eigenvectors of the correlation-normalised FIM with eigenvalue below 10⁻³ of the
largest, reported as their two heaviest loadings. Pairwise collinearity is |corr| from Σ.

Verdict per coordinate, in this order: **AT_BOUND** (within 10⁻⁶ of the band); **PINNED** if the 95 %
half-width 1.96·σ_marg is ≤ 0.5 dex for a rate, ≤ 30 kJ/mol for a barrier, ≤ 0.15 for a yield or
0.5 for a pKa; **WEAK** if ≤ 1.5 dex / 90 kJ/mol / 0.5 / 1.5; otherwise **UNIDENTIFIED**. A coordinate
is **collinear-not-insensitive** when its conditional σ passes the PINNED threshold and its marginal
does not.

Cross with the envelope: for each prior in `uncertainty.core_priors()` whose key names a fitted
coordinate, the prior's own half-width against 1.96·σ_marg: **prior-dominated** if the prior is
narrower, **data-dominated** otherwise.

## 4. Predictions

- **P1.** Across the shipped fits, **fewer than 40 %** of free coordinates come out PINNED.
- **P2.** B3 has **at least three** near-null directions, and its three competitor channels
  (`k_gln_glc`, `k_ala_glc`, `k_acr_ala`) are insensitive in the strong sense: conditional σ above
  1 dex, so no amount of collinearity accounting rescues them. The fit report's own halfwidths of
  10⁵–10⁶ dex already say so; the audit says whether the failure is insensitivity or collinearity.
- **P3.** In B8, `Ea_decay_thiol_sink` is AT_BOUND at its 102 kJ/mol ceiling, matching the existing
  profile's "bound_limited" call, and at least one further barrier is at a bound.
- **P4.** Collinearity, not insensitivity, is the dominant failure: **at least half** of the
  coordinates that are not PINNED have a conditional σ that *would* pass the PINNED threshold.
- **P5.** Of the envelope priors attached to fitted coordinates, **the majority are prior-dominated**:
  the declared band is narrower than what the data say. That is the honest reading of "the envelope
  covers 16 of 44": it is mostly the priors' width, not the fits'.
- **P6.** B18's two barriers sit in declared bands about 2.5 and 3.2 kJ/mol wide; they come out
  PINNED by the band, with a data σ wider than the band.
- **P7.** Nothing moves: no constant, no artifact the freshness gate watches except the new one.

## 5. What is built

`scripts/generators/generate_kinetic_core_b38_identifiability.py` writing
`results/validation/kinetic_core_b38_identifiability.{json,md}`, a test, a results-README entry, a
WAVES row, and a declaration amendment. The deliverable is the table and the list of named debts
re-ranked by whether a refit on the existing data could even see them.

## 6. Outcome (written 2026-09-11, after the run; artifact `kinetic_core_b38_identifiability.{json,md}`)

**Four predictions held, three were refuted, and the refutations say more than the holds.**

### The table in one line

**47 free coordinates across the five shipped fits: 12 pinned by the data (25.5 %), 6 sitting on a
declared bound, 29 weak or unidentified — and of those 29, only 9 are collinear. Twenty are simply
not seen.** The fits are starved, not tangled.

| fit | rows | free | χ²_red | pinned | at bound | weak | unidentified | of which insensitive |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B8 sulfur | 62 | 23 | 1.03 | 4 | 3 | 11 | 5 | 3 |
| B3 acrylamide | 30 | 11 | 5.55 | 1 | 1 | 1 | 8 | 8 |
| B18 dicarbonyl sinks | 10 | 6 | 56.9 | 0 | 2 | 3 | 1 | 2 |
| B20 glycation | 10 | 5 | 3.44 | **5** | 0 | 0 | 0 | 0 |
| B21 glucosone → glyoxal | 6 | 2 | 0.89 | **2** | 0 | 0 | 0 | 0 |

The two fits that are fully pinned are the two built on **fed or within-study rows on a small
network** (B20 on Nguyen's casein rates, B21 on Hamzalioglu's aqueous route). The two that are
mostly blind are the two that scored **end-of-cook levels through a large network**.

### P1 — HELD. 25.5 % pinned, against a prediction of "fewer than 40 %".

### P2 — REFUTED on its first clause, held on its second, and the reason is a method fact worth keeping.
B3 has **one** null direction, not three: `k_int1_mel` against `Ea_int1_mel` with correlation −1.00
— the rate-versus-barrier degeneracy of a step measured at effectively one temperature. The three
competitor channels are indeed dead to the data (conditional half-widths of 3×10⁴ to 4×10⁵ dex), but
they do not show up as null *directions*, because the eigen-analysis runs on the
correlation-normalised FIM, which puts 1 on every diagonal and so **cannot see insensitivity, only
collinearity**. Insensitivity lives in the conditional half-width. I predicted them in the wrong
column. Eight of B3's eleven coordinates are unidentified and all eight are insensitive: **the
acrylamide fit's thirty rows determine two numbers, `k_acr_dp` (±0.44 dex) and, weakly, its barrier
(±46 kJ/mol).** `Ea_competitor_sugar` sits on its 20 kJ/mol floor.

### P3 — REFUTED on its second clause. `Ea_decay_thiol_sink` is on its 102 kJ/mol ceiling as
predicted, but no second barrier is on a bound: `Ea_decay_carbonyl_sink` is at 167 inside (20, 250)
and unidentified (±166 kJ/mol). What *is* on a bound — and this was not predicted — is **both dimer
rates**, `k_dimer_mft` and `k_dimer_fft`, at their +0.5 ceiling. So the sulfur data push against
three bounds at once around the thiol sink: the barrier and both dimerisation rates. That is the
identifiability form of B27's finding that the disulfide deficit is a rate at its ceiling, and it
means **a refit on the present data cannot move the thiol sink; only the thiol-against-time
experiment (`EXPERIMENTS.md` #1) can.**

### P4 — REFUTED. 9 of 29, not "at least half". Collinearity explains under a third of what the
data fail to pin; the rest is rows that never touch the coordinate. In B8 the three insensitive
constants are `k_dimer_decay`, `k_thiol_decay` and `k_arp_dpo`; in B3 all eight; in B18 both
barriers. Where collinearity *does* act it is in B8's furanic-and-thiol cluster (`k_tdp_fur`,
`k_glc_ha`, `k_fur_decay`, `k_thiolate_loss`, `k_fur_fft`, `k_osone_decay`, `k_ttca_deg`,
`k_arp_tdp`), each of which the data would pin to ±0.2–0.5 dex if its neighbours were held.

### P5 — REFUTED, and correctly so. 29 priors sit on fitted coordinates; **21 are data-dominated,
8 prior-dominated.** The B8 priors are the fit's own search bands, ±5.25 dex, far wider than the
data on every coordinate the data see. Prior-dominated are exactly the coordinates the data are
blind to: the two dead B8 decays, `Ea_int1_mel`, `Ea_competitor_sugar`, both B18 barriers and both
B18 pH slopes. **The envelope is therefore not hiding a prior where a fit exists; it is honest
width where the data exist and declared width where they do not.** The prediction was wrong in the
direction that is good for the repository.

### P6 — HELD. B18's two barriers sit on the ceilings of bands 2.5 and 3.2 kJ/mol wide while the
data would allow ±290 kJ/mol: the band is the number. And B18's χ²_red of 56.9 says the fit itself is
poor — recorded here, not acted on.

### P7 — HELD. One artifact written; nothing else moved.

### Also found

- **Fourteen of B8's 25 frozen coordinates would be data-visible if freed** (conditional half-width
  inside the pinned threshold), including `k_mft_decay`, `k_h2s_loss` and the lumped formation
  barrier. Conditional visibility is not marginal identifiability, so this is a candidate list for a
  future B8 re-merge, not a verdict.
- **The B1 trunk report stores no per-parameter standard errors** (its script computes them; the
  report does not carry them). A gap, recorded.

## 7. The named debts, re-ranked by whether a refit could see them

| debt | what the audit says | so the next move is |
|---|---|---|
| the thiol sink (B25/B27) | barrier and both dimer rates on their ceilings; 3 sulfur decays dead to the data | **an experiment**, not a refit: thiol against time at two temperatures |
| the 3-deoxy pool lifetime (B37) | not a fitted coordinate today; the constants are Kocadagli bands | a **B21-style fed fit on Mittelmaier's pot** — small network, fed rows, the design that pins |
| the acrylamide competitor channels | insensitive: the competition rows never touch them | leave frozen; a refit is pointless without competitor-varying data |
| the lipid Q10 → barrier (B37) | not fitted at all; the lane has no fit | a declared-band install with a frozen before/after pair, as ENV-B34 did |
| the extrusion residence time (B37) | changes a FIT row's condition, and that row's fit sees two numbers | a re-run of B3 with 40 s declared, compared frozen-to-frozen |
| B18's dicarbonyl-sink barriers | pure band, χ²_red 57 | a premise check on the fit before any refit |

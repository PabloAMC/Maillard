# Wave B60 pre-registration — a converged pH trajectory and a flux budget that closes

*Written 2026-10-28 before any B60 fit, score or cost existed. The only B60 numbers computed before this
file are the convergence probes in §1 (three pots at B51's shipped vector). Kind 1 under
`scripts/generators/WAVES.md` ("How a change is accepted"): a NUMERICAL correction, accepted on
verification (§4); its refit is kind 2 and ships with it. Switch: `sulfur.PH_NODE_GRID`
("uniform" reproduces every wave through B58). Generators: `scripts/generators/generate_kinetic_core_b60_fit.py`,
`scripts/generators/generate_kinetic_core_b60_verification.py`.*

## 1. Why

The sulfur lane's dynamic pH (B2.2) is a fixed point solved at `ph_nodes` = 41 times EQUALLY SPACED over
the run; between nodes the rate vector is interpolated linearly in time. A flux attribution of five
reference pots (2026-10-28, `build_sulfur_flux_attribution.py`) found that in a cysteine + pentose pot the
in-situ pH falls most of its way in the first ~0.1 % of the run: Hofmann ribose pH 5 at 145 °C goes 5.27 →
4.47 by 0.07 min of 20; Cao's xylose pot at 90 °C goes 5.28 → 3.39 by 0.09 min of 300. The first uniform node
sits at duration / 40 (0.5 and 7.5 min), so through the whole transient every pH-dependent constant is
evaluated at the pH of the starting state. The solution is not converged:

| pot (B51 vector) | output | 41 uniform | 1601 uniform | 41 + 41 log-and-linear | 81 + 81 log-and-linear |
| --- | --- | --- | --- | --- | --- |
| Hofmann ribose pH 5, 145 °C, 20 min | MFT, mM | 1.47e-5 | 5.77e-4 | 6.02e-4 | 5.99e-4 |
| | FFT, mM | 2.25e-5 | 1.98e-3 | 2.08e-3 | 2.07e-3 |
| | free Cys, mM | 0.0044 | 0.967 | 1.02 | 1.01 |
| Cao xylose, 90 °C, 300 min | FFT, mM | 0.0433 | 0.199 | 0.253 | 0.250 |
| | free Cys, mM | 96 | 211 | 216 | 214 |
| Zhang 2025b xylose pH 6.5, 120 °C, 50 min | MFT, mM | 3.93e-6 | 7.33e-6 | 7.77e-6 | 7.78e-6 |

Every fit since B2.2 integrated its systems on the 41-point uniform grid, so the fitted constants have
absorbed a discretisation error of up to ~40x on Hofmann's MFT.

A second, smaller defect sits in `sulfur.sulfur_flux_budget`, which scores the one active `share` row
(`cerny_ternary_thiamine_share`) and every `branch_shares` diagnosis: it integrates the reaction rates on a
LINEAR output grid (101 points in the fit), which misses a transient lasting seconds, and it recomputes
the constants from the pH series alone, without B51's ionic strength or label-scale pH, so its rates are
not the rates the trajectory saw. A budget built that way does not close: in Cao's pot it books 2 575 mM of
3-deoxypentosone made from 400 mM of xylose.

## 2. What changes

- **`sulfur.PH_NODE_GRID = "log_and_linear"`**: the dynamic pH is solved at the union of the 41 uniform
  nodes and 41 nodes geometrically spaced from 1e-5 x duration to the end; the rate vector is interpolated
  between bracketing nodes by bisection (`sulfur.rate_vector_at`, shared by the ODE right-hand side and the
  flux audit). The prescribed and clamped pH paths are unchanged.
- **`SulfurRun.rate_nodes`**: the run keeps the node rate vectors it was integrated with.
- **`sulfur_flux_budget`** integrates on a log-dense output grid with the run's own node constants, and
  reports the closure of the cysteine, pentose and thiol balances it implies.
- Nothing in the network, the parameters, the data or the objective changes.

## 3. Candidates

- **(A) B60 refit.** B51's 52-row objective with the new grid; same free set, bands, sigmas, weighting,
  optimiser, budget and two-start protocol (start 0 = B51's optimum, start 1 = its perturbation).
- **(B) B51's frozen vector on the new grid** (what the correction alone does).
- **Baseline:** B51 as shipped.

## 4. Acceptance — verification, not scores (kind 1)

B60 ships when all four hold; the scores of (A) and (B) are recorded as diagnosis, and the rows the
correction makes worse are listed on the model card (`KNOWN_EXPOSED_ASSUMPTIONS`).

- **V1, convergence on the fit panel:** at B51's vector, every one of the 52 fit-row observables changes by
  at most 0.01 dex (2.3 %) between `ph_nodes` 41 and 81 on the new grid.
- **V2, an independent reference:** on the three probe pots of §1, the new grid at 41 agrees within 0.02
  dex on MFT, FFT and free cysteine with a brute-force uniform grid fine enough to resolve the transient
  (25 601 nodes); where the brute-force grid itself is unconverged, that is recorded instead.
- **V3, the flux budget closes:** on the share-row system and the five reference pots, initial − final
  equals the net integrated flux within 1 % for PENT, Cys, H2S, MFT and FFT.
- **V4, frozen waves reproduce:** every frozen-wave test passes with `PH_NODE_GRID = "uniform"`.

(A) is the shipped vector if its cost is no higher than (B)'s on the same objective (the optimiser did its
job); otherwise (B).

## 5. Expectations

- V1–V4 should pass by construction of the grid; probability of shipping **about 0.9** (the risk is a fit
  system with a transient faster than 1e-5 x duration, or an LSODA slowdown that makes the refit budget bind).
- At B51's vector the Hofmann pH 5 thiols rise 40–90x, which removes most of the "~100x low" those rows
  showed under B51; B51's diagnosis of that gap (the thiolate-acting thiol sink) is therefore at least
  partly a numerical artefact, and the refit will move the trunk and sink constants that compensated it.
- What B60 does NOT fix: the trunk's timescale. In every cysteine + pentose pot the lane destroys the
  cysteine in seconds (Strecker step `r_dpo_ddp` after the cysteine-catalysed enolisation), against Cao 2017's
  ~0 % cysteine destroyed in 5 h at 90 °C and Zhang 2025b's ~48 % in 50 min at 120 °C, pH 6.5 (§1 of
  `docs/guides/SULFUR_FLUX_MAP.md`). That is a model-structure and data problem, B59, kind 3.

## 5a. Amendment 1 (2026-10-28, while building V1, before any B60 fit or score existed)

Measuring V1 at the optimiser's own "quick" setting (9 nodes, 2 pH sweeps, rtol 1e-5) exposed a second
numerical defect: **the pH fixed point is not converged, and the code does not check it.** The dynamic pH
is a Picard iteration (sweep: integrate with the node pH, recompute the node pH from the states); it stops
after `ph_iterations` sweeps whatever the residual. In the Kang TTCA pots (B51 vector, 41 nodes) it
OSCILLATES between two branches — final pH 5.52 / 6.55 / 5.56 / 6.37 after sweeps 2 / 3 / 4 / 5, residual
~1 pH unit — and settles only after ~20 sweeps (residual 1.4e-3; MFT 3–5x the 4-sweep value). The fit
stopped at 2 sweeps, the careful scorers at 4. The residual was recorded in the run's metadata
(`ph_fixed_point_residual`) and read by nobody.

Added to §2: the sweep is under-relaxed (node pH ← node pH + ω·(fresh − node pH), ω declared in
`sulfur.PH_FIXED_POINT_RELAXATION`), iterated until the residual is below 1e-4 pH or a declared cap, and a
run that ends above 1e-2 pH is flagged in its metadata (`ph_fixed_point_converged`). The fit's quick setting
is chosen, before the fit, as the cheapest (nodes, sweeps, rtol) whose 52 observables lie within 0.03 dex of
the careful setting at B51's vector; the choice and its V1 numbers are recorded in the verdict.

Added to §4: **V5, the fixed point converges:** at B51's vector every one of the fit's systems ends with a
residual below 1e-3 pH at the careful setting.

## 6. Verdict

*(Written after the run, below this line.)*

**SHIP (A), the B60 refit.** Written 2026-10-28, after the fit (`kinetic_core_b60_fit_report.json`, two
starts, `kinetic_core_b60_members/`) and the verification (`kinetic_core_b60_verification.{json,md}`).

**Verification (§4 and Amendment 1), all pass:**

| check | threshold | measured |
| --- | --- | --- |
| V1, 41 vs 81 nodes on the 52 fit rows, B51's vector | ≤ 0.01 dex | 0.0010 dex |
| V2, probe pots vs an independent fine reference | ≤ 0.02 dex | ≤ 0.005 dex vs a 321-node refined grid |
| V3, flux closure (PENT, Cys, H2S, MFT, FFT) | ≤ 1 % | ≤ 0.43 % |
| V4, frozen-wave tests on the legacy numerics | pass | 129 pass (`test_kinetic_core_b2{,_1,_2,_4}.py`), after deviation 4 |
| V5, fixed-point residual over the fit's 37 systems | < 1e-3 pH | worst 9.6e-4 pH, at most 10 sweeps |

Before the correction, 17 of the 52 fit observables at B51's vector were off their converged value by
more than 0.03 dex, and the worst by 0.833 dex.

**Deviations from the text above, all recorded in the generators:**

1. *V2's reference.* The 25 601-node uniform grid named in §4 was replaced by a 6 401-node uniform grid
   plus a 321-node refined grid on the new layout, because the 25 601 run had been launched on the
   pre-amendment numerics and was stopped. The uniform 6 401 grid is itself unconverged on Cao's MFT,
   0.145 dex off the refined grid, which §4 allows ("where the brute-force grid itself is unconverged, that
   is recorded instead"). The convergence argument therefore rests on V1 together with the refined grid,
   not on a brute-force uniform grid.
2. *Fixed-point tolerance.* Amendment 1 says the iteration runs "until the residual is below 1e-4 pH". The
   code uses 1e-3 (`sulfur.PH_FIXED_POINT_TOLERANCE`), which is V5's threshold. Below about 1e-3 the
   residual reaches a floor of 3e-5 to 6e-4 pH set by the ODE tolerances, so 1e-4 is not always
   reachable; the reason is in the constant's comment.
3. *Quick setting.* 17 pH nodes at rtol 1e-5 (0.0076 dex from the careful 81-node setting). 9 nodes
   missed the 0.03 dex tolerance, at 0.0315. After the fit, the generator was changed to install this
   setting at import rather than inside `main()`, so the Laplace and the profile differentiate the same
   objective the optimiser minimised. The Laplace reproduces the optimum cost of 27.1320.
4. *A third legacy switch.* §2's closed flux budget also changes the share rows (`branch_shares` of
   `sulfur_flux_budget`). With only the grid and fixed point pinned to their legacy values, B2.4's
   frozen-control cost moved by 4.6e-5 relative, past the test's 1e-5. `sulfur.FLUX_BUDGET` now names the
   two budgets, `"closed"` and the legacy `"linear_rederived"`, and the frozen fixtures pin all three
   switches. The B60 fit, the verification and the flux attribution all ran on `"closed"`.
5. *A Laplace correction found on the way.* The covariance is a pseudo-inverse on JᵀJ's eigenbasis that
   drops near-null directions. A dropped direction got zero variance where it should get unbounded
   variance, so three constants no row responds to were labelled IDENTIFIED with σ between 1e-9 and
   0.02. `k_glc_ha` has an exactly zero Jacobian column; the other two are `k_thiol_decay` and
   `k_dimer_decay`. B51's covariance has the same defect. The envelope then pinned those constants at
   the optimum instead of drawing them across their bands. `generate_kinetic_core_b8_laplace.py` now
   records each coordinate's weight in the dropped subspace and writes σ = null (unidentified) above
   1e-6. The B60 covariance was regenerated with it.

**Ship rule:** (A) 27.132 ≤ (B) 33.156 on the same objective (52 rows, quick setting), so (A) ships.
Start 1, from the perturbed vector, reached 27.671, and start 0 is the reported optimum. The active
bounds are `k_dimer_mft` (upper) and `ph_acid_yield_per_sink_event` (lower). `Ea_decay_thiol_sink` came
off its 102 kJ/mol ceiling to 101.3.

**Diagnosis: the scores were flattered by the unconverged numerics.** As shipped, B51 reported a cost of
20.41. That number came from the 2-sweep oscillating fixed point. The same vector integrated correctly
costs 33.16, and the refit recovers only part of that, to 27.13. The level sum of squares over the 52
rows is 40.8 for B51 as shipped, 66.3 for (B) and 54.3 for (A). These are the rows that get worse (A
against B51 as shipped, by more than 0.5 σ); they go on the model card:

| row | B51 shipped | (B) | (A) |
| --- | --- | --- | --- |
| `kang_120C_MFT` | +1.19 | +3.01 | +2.51 |
| `kang_100C_MFT` | +0.76 | +2.47 | +1.75 |
| `zhai_MFT_fold_120_over_100` | +1.16 | +1.43 | +2.03 |
| `kang_100C_MFT_over_FFT` | +0.91 | +1.98 | +1.56 |
| `kang_120C_FFT` | +0.27 | +1.69 | +0.86 |
| `hofmann_ribose_FUR_insitu` | +0.04 | +0.96 | +0.62 |
| `zhou_pH7_FUR_arp_alone` | −0.85 | −0.84 | −1.38 |
| `kang_100C_FFT` | +0.08 | +0.99 | +0.59 |

Three rows improve by more than 0.5 σ: `cerny_ternary_thiamine_share` (+1.30 → −0.34),
`whitfield_nf_h2s_MFT` (−1.49 → −0.81) and `hofmann_ribose_NF_insitu` (+1.23 → +0.69).

**Reading.** The rows that get worse are almost all thiol over-predictions in the cysteine + pentose pots
at 100–120 °C (Kang, Zhai's temperature fold). The old numerics' oscillating pH held these down by
accident. On a converged trajectory the lane makes too much MFT and FFT in exactly the pots where
`docs/guides/SULFUR_FLUX_MAP.md` finds the trunk 10³–10⁴× too fast: cysteine gone in seconds where Cao
and Zhang measure hours. The refit cannot repair this, because no constant in the free set slows the
trunk; the parameter moves are all ≤ 0.7 dex (`k_ttca_deg` −0.68, `k_dimer_fft` −0.55, `k_thiol_decay`
+0.41, `Ea_decay_carbonyl_sink` −6.7 kJ/mol). So the correction exposes the trunk-timescale defect that
B59 addresses rather than creating a new one. B51's §7 diagnosis, which blamed the thiolate-acting thiol
sink for the Hofmann pH 5 misses, is superseded; see `tasks/lessons.md`, 2026-10-28.

The 41 → 81 node agreement settles the numerics. It says nothing about the trunk's structure.

# Premortem: inference track H0-H4 + wave V (as of 2026-12-31)

What the premortem rests on, from EX-056:
- 6 thiol labs, 5 of them scorable. Jiangnan holds 341 of the 383 paired held-out rows.
- R3En has 16 coordinates, not 20 - 40.
- In-sample absolute-class tau is 0 - 0.17 dex; the held-out error is ~1 dex.
- Folds move constants by 3 - 4.5 dex, which rule (c) flags as a warning sign.

## Causes

**1. H0 is circular and has no second split. P 0.7.**
- Tuning tau until the LOLO residuals cover 0.9 makes coverage 0.9 by construction.
- 5 labs give a coverage SE of at least 0.13, since rows within a lab are correlated.
- Leaks bias tau low: R3En and the series noise were both chosen on all labs, and the warm starts come from the full optimum.
- The Jiangnan fold (67 training rows) biases tau high. Pooled by row, it is the whole answer.
- *Early sign:* tau changes by more than 2x when one lab is dropped.
- *Fix:*
  - Nested LOLO: calibrate on the inner folds, score the outer lab.
  - Weight by lab, not by row.
  - One tau per absolute vs relative class.
  - Declare tau a lower bound.
  - The real second split is a different lane's LOLO, or labs newly keyed by P1, frozen before computing.
  - With 10 - 20 labs: about 20 repeated grouped 2/3 : 1/3 lab partitions, stratified by method class, reporting the spread across partitions.

**2. Calibrated but useless. P 0.6.**
- An error of ~1 dex needs a 90 % half-width of ~1.6 dex. That is wider than the nulls' intervals (1.11 dex wide at coverage 0.53 - 0.63).
- *Fix:* accept on a proper score (paired log density plus interval score at matched coverage).
- Pre-declare: if absolute tau > 0.8 dex, ship relative within-lab predictions only.

**3. Lab effects on k cannot be identified. P 0.6.**
- Most constants are informed by 1 - 2 labs, so the hyper-SD posteriors are just their priors.
- A lab shift in ln k is nearly collinear with an observation offset when the product is far from depletion.
- *Fix:* by default, lab effects go on observation levels (kept analytically marginalised) and on bounded latent conditions: heat-up lag, pH drift, buffer. A lab effect on k needs 3 or more labs.

**4. The hierarchy hides wrong structure. P 0.5.** In sample, the constants already act as lab effects, which is why tau collapses to 0.
- *Detection:*
  - (a) lab effects that regress on recorded pH, T or buffer point to a missing mechanism;
  - (b) lab effects cannot change within-lab frozen signature verdicts; if they do, they are absorbing structure;
  - (c) a tau_k above 0.3 dex on a physical constant is structural.
- *Fix:* H2's negative control fits a wrong structure, with lab effects on k, to data simulated on EX-056's lopsided design, and must detect it.

**5. The JAX port stalls. P 0.55 that H3 overruns 2x.**
- pH is a nested solve: brentq inside an ionic-strength fixed point, with a warm-start memo. Node coupling is an outer Picard loop around whole integrations, so only inline (B64) mode can be ported.
- The clips and caps give kinked gradients.
- diffrax has no BDF and no DAEs. `vmap` runs in lockstep with the stiffest pot, and compiling takes minutes.
- *Recommendation:*
  - Spike the charge balance for 1 day first.
  - Compile A2's arrays into CasADi/IDAS with pH as an algebraic state (an index-1 DAE). That gives exact forward and adjoint sensitivities in float64 on CPU, with no nested root. `numbalsoda` is the alternative for raw speed.
  - Keep JAX for the forced-input layers.
  - At 16 - 30 dims, run SMC or ensemble samplers on the compiled engine, or GP emulators with history matching. History matching carries an explicit discrepancy term, so structural error is named rather than hidden.

**6. NUTS pathologies. P 0.6.**
- Sloppy spectra (Hessian eigenvalues spread over more than 6 decades) collapse the step size and hit max tree depth.
- Warmup reaches regions where the solver hits `max_steps` (NaNs, divergences).
- A funnel with 5 groups.
- Gradient noise when the solver tolerance is loose.
- *Fix:*
  - Initialise from Laplace or Pathfinder, with a dense mass matrix from the Laplace covariance.
  - rtol at or below 1e-8; bounds through a logit, never a clamp.
  - Count solver failures.
  - A pre-declared sampler-health stop rule.
- P(zero divergences) ~0.3.

**7. GPU / Colab. P 0.6 that the GPU gains less than 3x over the 8-core CPU.**
- FP64 TFLOPS, corrected: T4 ~0.25 (1/32), L4 ~0.47 (1/64), A100 9.7 (19.5 on tensor cores), H100 26 (PCIe) or 34 (SXM), roughly double on tensor cores. The G4 (RTX PRO 6000) is 1/64.
- Small stiff solves are bound by launches and lockstep, not FLOPs. The GPU pays only at about 1e3 or more solves per batch.
- Colab facts:
  - Pro: $9.99 for 100 compute units (CU). Pro+: $49.99 for 500 - 600 CU.
  - A100: ~5.4 CU/h (40 GB) and ~7.5 CU/h (80 GB), measured in March 2026. H100: no published rate.
  - Premium GPUs are "subject to availability", so the A100 is not reliable.
  - 24 h background execution is Pro+ only.
- Money is trivial (a 24 h A100 run is about $13); the risk is disconnects.
- *Fix:* benchmark the forced-input layer under `vmap` on CPU, L4 and A100 first.

**8. The speed-ups are in the wrong order. P 0.5.**
- No profiling step.
- "10 - 100x" holds only for one gradient in about 100 dims. Gauss-Newton at p = 16 needs a full Jacobian: about 3 - 8x, mostly from compiling.
- *Fix:* the order should be profile, then pool, then compiled RHS + sensitivities on CPU, then a gated JAX / GPU step.

**9. Dependencies and opportunity cost. P 0.5.**
- H0 only post-processes the EX-056 JSONs, so it does not need A3 phase 2.
- The binding constraint is TEST data.
- Adding jax can bump the numpy / scipy pins and drift frozen evidence.
- *Fix:* drop the H0 dependency; a separate `maillard-jax` env with a gate that the core lock and `core-scores` are unchanged; cap H at about 20 h before the structure fixes land.

**10. Wave V overruns and biases. P 0.65.**
- "Every cited number" is thousands of numbers (521 rows, 562 cells, 297 contrasts), many read off figures, so not 1 day.
- Re-reading only what disagrees with the model is asymmetric scrutiny.
- TEST corrections made after the 11-07 scoring are forking paths.
- Cluster agents editing the register will collide.
- *Fix:*
  - Verifiers are blind to engine residuals.
  - Full re-read of the TEST and high-leverage rows; the rest by stratified sample. The atlas needs n ~150 to tell 9 % from 4 %.
  - Errata go in with a new split hash before scoring.
  - Agents write patches to the scratchpad; one writer applies them serially, after the in-flight batch.

**11. The V3 stop rule. P 0.5.**
- "Two empty hops" means paying for hop 2: about 1e3 papers at hop 1, about 3e4 at hop 2.
- "Discriminating" needs the full text.
- OpenAlex reference lists before 1990 are patchy.
- *Fix:* best-first search per mechanism. Stop when 10 triaged abstracts yield fewer than 1 discriminating candidate, or at a cap of 30 abstracts and 5 full texts. Old papers' references come from the PDFs on disk.

## Calibration

| claim | draft | mine | note |
|---|---|---|---|
| H0 by 11-07 | 0.6 | 0.7 without A3; 0.35 with it | |
| H2 passes SBC | 0.5 | 0.35 | needs ~100+ SBC fits: hours, not minutes |
| H3 | 0.35 | 0.2 | A2a not started |
| converged posterior by mid-December | 0.3 | 0.15 | |
| intervals widen > 2x | 0.6 | 0.85, conditional on calibration | H0 alone does it, so it does not test H |
| hierarchy beats the nulls | 0.15 | 0.1 - 0.15 | agree |

## Fold into the plan

1. H0: nested LOLO, lab-weighted, tau a lower bound, no A3 dependency.
2. Second split: a cross-lane LOLO or P1's labs; repeated grouped partitions once there are 10 or more labs.
3. Accept on a proper score against the nulls, not on coverage alone.
4. Pre-declared: absolute tau > 0.8 dex means relative predictions only.
5. Lab effects on observation levels and latent conditions; on k only with 3 or more labs.
6. Masking tests: covariate regression, within-lab signature invariance, the lopsided negative control.
7. Profile first, then a 1-day pH-DAE spike (CasADi/IDAS vs JAX) before any port.
8. Posterior at 16 - 30 dims by SMC, ensemble or emulator history matching; NUTS last.
9. A separate JAX env with a core-lock gate; benchmark the GPU on the forced-input layer.
10. V: blind verifiers, influence-ranked plus sampled, errata before scoring, one writer, a marginal-yield stop.

## Verdict

| stage | verdict | change |
|---|---|---|
| H0 | change | nested; cross-lane check |
| H1 | change | lab effects on observation / conditions; masking tests |
| H2 | keep | as the tool bake-off and the negative control |
| H3 | change | CasADi / compiled path, gated on the pH spike and a lane passing its structure tests |
| H4 | drop for 2026 | |
| V | change | rescoped, ~2 days |

Sources: [McCormick, Colab GPUs (Mar 2026)](https://mccormickml.com/2024/04/23/colab-gpus-features-and-pricing/); [Colab pricing overview](https://aitoolsatlas.ai/tools/google-colab/pricing).

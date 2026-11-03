# Red team: Munich sulfur cluster (CF-016, CF-065), 2026-11-03

This review is read-only. I re-ran the cluster's engine assembly (`generate_kinetic_core_b60_fit` frozen vector + fitted
pH drift) through `./scripts/docker_maillard.sh run "nice -n 10 python ..."`. My probes are in the session scratchpad
under `redteam_munich/`:
- `rt_probe.py`: the confound split, the grids and the time course.
- `rt_sens.py`: pH, the Cys -> H2S scale, volume and a thiazolidine proxy.
- `rt_loss2.py`, `rt_slow.py`, `rt_furdecay.py`: sinks switched on and off.

I read the PDFs through their text layer: munch1997 Table 5 and Methods; hofmann1998 Tables 3 and 5 and p. 239; and
hofmann1995, 1995b, 1998b, 1998c, schieberle1996 / 1998 / 2000 and munch1998 for furfural.

## Claim 1: the engine over-predicts Munch T5's Cys-vs-H2S FFT contrast by +0.71 dex, through fast FFT formation then ~99 % loss

**The number: UPHELD (P 0.97).** I reproduce R_FUR = 0.234, 0.83 mol % in the H2S arm and 0.19 mol % in the Cys arm.

**The unexplained 2.8 vs 0.83 mol %: EXPLAINED.** The generator assembly and `core_parameters(SULFUR)` agree on every
float parameter (0 differences). The whole gap is the pH mode:

| parameters | pH mode | mol % | R_FUR |
|---|---|---|---|
| generator or core | fitted drift (5.0 -> 5.34 in situ) | 0.83 | 0.234 |
| generator or core | clamped 5.0 | 2.82 | 0.099 |

The FFT sinks act through the thiolate (`ch_thiolate_loss_fft`) and formation through HS-, so 0.34 pH units are worth
0.53 dex on the absolute value and 0.37 dex on R. The results' one-line attribution is therefore right, but it was
asserted, not shown, because the probe changed the parameter source and the pH mode together.

**Grid convergence: the cluster's own flux number is not converged.** In the H2S pot, the `r_fur_fft_hs` integral and
the closure are:

| grid points | r_fur_fft_hs integral (mM) | closure |
|---|---|---|
| 201 | 64.6 | 70 % |
| 2001 | 20.5 | 5.6 % |
| 20001 | 19.4 | 0.1 % |

The 20.5 mM quoted in the register is more than the 20 mM of furfural charged, which is impossible. The converged
figures are 19.4 mM, plus 0.14 mM via `r_fur_fft`. That is 97.7 % of the furfural, and 99.1 % of the FFT formed is
then destroyed. The final concentrations do not depend on the grid. "Within minutes" understates the speed: 88 % of
the furfural is already FFT at 15 s.

**Robustness of the +0.71 (shipped engine, one change at a time):**

| change | R_FUR | dex vs 0.046 |
|---|---|---|
| 100 mL (CF-102; Munch's Methods line covers Tables 4-6) | 0.194 | +0.62 |
| pH clamped 4.8 / 5.2 / 5.6 | 0.074 / 0.137 / 0.302 | +0.21 / +0.47 / +0.82 |
| k_cys_h2s x0.5 / x0.2 / x0.1 | 0.151 / 0.073 / 0.039 | +0.52 / +0.20 / -0.07 |
| irreversible furfural-Cys adduct: 50 % / 80 % of furfural bound at t0 | 0.155 / 0.087 | +0.53 / +0.28 |
| k_fur_decay x3 / x6 (H2S arm still 0.80 / 0.77 mol %) | 0.023 / 0.008 | -0.30 / -0.77 |
| k_fur_decay x0.3 | 1.08 | +1.37 |
| FFT loss (k_fft_decay, k_thiolate_loss) x0.3 | 0.046 | 0 (but the H2S arm goes to 23 mol %) |
| FFT loss x0 | 0.026 | -0.25 |

- **Zheng escape.** The register's "~5x below Zheng" is not enough: 5x leaves +0.20 dex, and closing the gap needs
  about 8x. A ±20 kJ/mol error in the activation energy over the 35 C extrapolation is worth only ~1.7x. P 0.15
  (the cluster said 0.3).
- **Thiazolidine escape.** The engine carries HMF + Cys (Hamzalioglu, 23 M^-1 day^-1 at 50 C) but no furfural + Cys
  adduct. An HMF-like rate is <= ~3 % of `k_fur_decay` even with a generous extrapolation, so it is negligible.
  Only a fast equilibrium binding >= 80 % of the furfural helps, and if bound furfural is protected from decay the
  sign flips, because R rises steeply with furfural lifetime. P 0.1 that this explains most of the gap.

**The inference: REFUTED as stated.** The claim is that the contrast "counts against product loss" and that the data
"favour slow formation with a surviving product". I give P 0.2 that the contrast discriminates on FFT survival at all.

- In the engine, R is governed by the furfural lifetime in the Cys arm, where the H2S trickle competes with
  `r_fur_decay` (0.38 1/min).
  - R falls roughly as k_fur_decay^-2: x3 gives 0.023.
  - Hofmann's absolute value is untouched (0.80 %). So a factor-2 change in `k_fur_decay` closes the contrast while
    keeping fast formation and ~99 % FFT loss.
  - `k_fur_decay` is fitted, has no literature value, and shares the carbonyl-family 143.5 kJ/mol that B62 already
    flags against Yaghmur's measured 46.5. The misfit protocol's "fitted constant far from its measured analogue"
    sign points here first.
- The cluster's own world fails worse with the engine's furfural sink.
  - Slow formation (k_fur_fft x1/3000) with zero FFT loss gives R = 0.0055 with the engine's `k_fur_decay`
    (-0.92 dex) and only 0.020-0.023 with stable furfural.
  - The prereg's donor band 0.028-0.075 is a ratio of H2S integrals. It is the right null only if the co-reactant's
    time course is flat, and furfural's is not.
  - So H2c's "in band on the furfural side" is conditional on furfural surviving 20 min at 145 C. That is the very
    quantity decider (ii) would measure.
- What the contrast does constrain is the joint (furfural lifetime x Cys -> H2S rate x FFT survival).
- Independent evidence on FFT survival does exist.
  - Kumazawa 2003 measures 89 / 80 % survival at pH 5.0 / 5.4, 121 C, 10 min. Extrapolated first order with the
    engine's 101 kJ/mol, that is 7-26 % survival at 145 C / 20 min, against 0.85 % in the engine's fed pot.
  - The weak time dependence at pH 6.0 (48.3 -> 43.4 % for +10 min) suggests a saturating, oxidant-limited sink.
  - So there may be a case that the engine's FFT loss is too fast, but it comes from Kumazawa, not from Munch T5.

## Claim 2: H2S order n = 1.04 [0.90, 1.18]; engine passes at -0.03 dex

- **The values: UPHELD (P 0.97).** The text layer reads 130.1 / 550.8 ug, matching the cluster's eye read; the
  arithmetic and the ±10 % band check out.
- **Integrity nit.** Expt 2's printed 0.12 mol % does not follow from its 130.1 ug: 1.140 umol / 1 mmol = 0.114 %.
  The printed % implies about 137 ug, a ratio of 4.02 and n = 1.00. That is immaterial to the verdict, but Munch
  Table 5 is not internally clean.
- **The reading: WEAKENED (P 0.25 that it supports cluster 3's rate-constant transfer).** Every structure I ran gives
  a dose ratio of 3.65-4.1 against 4.23 observed:
  - shipped (H2S-limited, ~complete conversion then loss) 3.95;
  - no FFT loss 3.91;
  - no furfural decay 4.09;
  - slow formation 3.71-3.88;
  - k_fur_decay x3 3.81.
- The ratio shows only that the FFT yield scales with the H2S charge. That rules out a saturating or threshold YIELD
  at 5-20 mM. It does not show low-conversion second-order kinetics, which is what carrying Hofmann's constant to
  Zhai's pot needs.

## Claim 3: CF-016, "no third print" and "67.5 is not a misprint"

- **No third print: UPHELD (P 0.95).**
  - The cluster's list omits hofmann1998c (Z Lebensm Unters Forsch A 207:229, Table 4). It prints FD factors only;
    furfural in aqueous ribose/Cys column I is "~1".
  - hofmann1995's text layer (61k characters) has no furfural string at all.
  - None of the on-disk Munich-group PDFs prints a furfural level for this pot beyond the two known prints.
- **67.5 is not a misprint: UPHELD (P(single-cell slip) 0.04).** The p. 239 sentence is explicitly about that table:
  "As shown in Table 5, nearly 10 times more furan-2-aldehyde was generated in the ribose compared to the glucose
  system".
  - The caveat: the same authors wrote the text from the same data. The sentence excludes a typesetting slip in the
    cell, not an upstream calculation or internal-standard error that would feed both (the "column high" branch).
  - A pentose making far more furfural than a hexose is also the chemical prior, so the text adds little beyond
    excluding the slip.
- The remedy (keep 67.5, sigma_log 0.5; carry 5.2 as a choices-register alternative; do not average) is sound.

## Claim 4: the two Table 5 ratio rows (count once under CF-101)

- **(a) Dose ratio: add only as a VALIDATION row; do not make it a fit row.**
  - It has near-zero discriminating power (claim 2).
  - Its denominator is expt 3, which is already the in-fit `hofmann1998_t3_fft_c5_h2s`. As a fit row it would weight
    expt 3's error twice, which defeats "count once".
  - The sigma of 0.06 dex assumes Hofmann's ±10 % replicate bound for Munch-only cells whose replication is
    unstated, alongside the 0.12 / 0.114 % inconsistency. Use sigma >= 0.1 dex.
- **(b) Cys-vs-H2S ratio: ADD as a HELD-OUT DIAGNOSTIC.** At sigma 0.1 the engine fails by about 7 sigma. Its
  description must be neutral: it constrains furfural lifetime x Cys -> H2S rate x FFT survival jointly, and is not a
  test of product loss.

## REQUIRED changes to register_edits.yml (before the register's writer applies them)

1. **CF-065, test (4).**
   - Replace "r_fur_fft_hs integral 20.5 mM, closure 5.6 % on a 2001-point grid" with "r_fur_fft_hs + r_fur_fft
     19.5 mM (97.7 % of the furfural; closure 0.1 % on a 20001-point grid; 2001 points over-count by 6 %)".
   - Replace "within minutes" with "within ~15 s".
   - Add: "0.83 vs 2.8 mol % is the fitted pH drift (5.0 -> 5.34) vs a pH-5.0 clamp; the parameter sets are
     identical".
2. **CF-065, test (4), last sentence.** Replace it with: "The contrast is governed in the engine by furfural's
   lifetime in the Cys arm: k_fur_decay x3 gives R 0.023 with Hofmann's 0.48 % unchanged; slow formation with a
   surviving product gives 0.0055 under the engine's furfural sink. It constrains furfural lifetime x k_cys_h2s x FFT
   survival jointly and does not discriminate product loss. Closing it via k_cys_h2s needs ~8x (P 0.15); a
   furfural-Cys adduct at HMF-like rates is negligible."
3. **CF-065, test (1).** Replace "so cluster 3's second-order transfer holds on a row outside the conflict" with "the
   yield scales with the H2S charge. Rate-limited and H2S-limited (complete conversion then loss) structures both
   predict 3.7-4.1, so this excludes a saturating yield, not the low-conversion assumption behind cluster 3's
   transfer".
4. **CF-065, remedy.**
   - Restore FFT "product loss in Hofmann's hold" to 0.25. If it is lowered at all, do so on Kumazawa's
     extrapolated survival, not on Munch T5.
   - Hold the MFT "Cys-relieved depression" at 0.25 until the engine's own R_NF = 0.14 has been decomposed. The
     engine already carries `r_mgo_cys_h2s` (Strecker H2S release catalysed by NF-derived MGO) and `ch_cys_ox`, so
     "more than Zheng thermolysis" is not "more than an H2S donor". I did not decompose it.
   - Dose-ratio row: "validation only, sigma >= 0.1".
5. **CF-065, explanation 4.** Delete "(the shipped engine is this world); 0.046 is observed". The ~0.2 prediction
   belongs to the engine's fitted furfural sink, not to the product-loss world as such.
6. **CF-016, test (Munich cluster).** Add hofmann1998c Table 4 (FD factors only) to the list of non-prints. No status
   change.
7. **results.md table row "munch1997 Table 5 ... no pH, no volume".** Munch's Methods give 100 mL and 0.5 M phosphate
   for Tables 4-6; the pH is unstated, and expt 3 is pH 5.0 via Hofmann T3. Not a register change, but correct it if
   results.md is revised.

**Decider list (CF-065).** Keep (ii), residual furfural at 20 min, now as the decider of THIS contrast. Add the direct
product-loss test: authentic FFT at 20 mM heated 20 min at 145 C, pH 5, in 0.5 M phosphate, with and without H2S,
measured by SIDA.

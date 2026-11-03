# Results: Munich sulfur cluster (CF-016, CF-065), 2026-11-03

Pre-registration: `prereg.md` (this folder). Scripts and outputs (not archived): session scratchpad
`conflicts/munich_sulfur/scripts/munich_t5_probe_b60.py`, `munich_closure.py`, `munich_t5_probe_b60.json`.
Engine runs: the shipped B60 parameter set, assembled by `generate_kinetic_core_b60_fit` from its frozen vector, with
the fitted pH drift. That assembly reproduces the B60 report's `fed_furfural_h2s_FFT` residual (0.83 mol% against
0.48, +0.24 dex). The same parameters run through `core_parameters` with `ph_drift=None` give 2.8 mol%. That
difference is a pH-drift effect; it is not chased here.

## Page reads (eye, 110 dpi renders)

| source | table / page | what it prints for the ribose 10 + Cys 3.3 mmol pot | furfural level? |
|---|---|---|---|
| munch1997 | Table 4, p. 1341 | FFT only (22.9 / 12.1 / 1.2 ug at pH 3 / 5 / 7) | no |
| munch1997 | Table 5, p. 1341 | fed furfural 1 mmol: + Cys 1 mmol 25.3 ug (0.02 %); + H2S 0.25 mmol 130.1 ug (0.12 %); + H2S 1 mmol 550.8 ug (0.48 %); no pH, no volume | (fed, not formed) |
| munch1998 | all tables (dossier + text layer) | yeast extracts and Strecker pots; furfural not analysed | no |
| hofmann1995 | AEDA paper (dossier: no concentration at all) | FD factors only | no |
| schieberle1996 | Table 3, p. 178 | 12 odorants, ug/L; no furfural row | no |
| schieberle2000 | Table III (text layer) | reprint of the 1996 / 1998 aqueous column; no furfural row | no |
| hofmann1998 | Table 5, p. 238 | FA 67.5, FFT 12.1, NF 54530.0, MFT 19.8 ug (100 mL, pH 5.0) | 67.5 (confirmed) |
| schieberle1998 | Table IV, p. 325 | aqueous column II: furan-2-aldehyde 52 ug/100 mmol ribose, OAV "<1" | 52 (confirmed) |

## CF-016

- **H1a (third print): FALSIFIED.** No on-disk source prints a third furfural level for this pot (P 0.95 that none
  exists in the corpus; the 5 % allows for an OCR miss in hofmann1995, which I did not render).
- **H1b (a 10x slip in Hofmann's ribose cell): REJECTED.** Hofmann's own text on p. 239 (left column) reads "nearly
  10 times more furan-2-aldehyde was generated in the ribose compared to the glucose system". The table gives
  67.5 / 6.9 = 9.8. A slip in the ribose cell would make that ratio ~1, so the slip is ruled out (P 0.03). A scale
  error over the whole column is not tested by this.
- **H1c (observation).** Except for furfural, the Schieberle 1998 aqueous column reprints Schieberle 1996 Table 3
  cell for cell (FFT 121, 2-acetyl-2-thiazoline 7, thiazine 424, 3-mercapto-2-pentanone 599, MFT 198 ug/L). Hofmann
  1998 Table 5 does the same for FFT, MFT and NF (54530 ug/100 mL = 545300 ug/L). **So each furfural number is the
  only print in its paper family.** They are two separate determinations, probably from different heatings and
  work-ups. They are not two readings of one measurement. Hofmann's method is stated: 5-methylfurfural added to the
  reaction mixture after cooling, before ether extraction, with an MS response factor of 0.95. Because the internal
  standard goes in before extraction, the result is corrected for recovery. Schieberle 1998 states no furfural method.
  The SIDA standards it names are for the thiols. Its other aqueous furanaldehyde (5-methylfurfural in the rhamnose
  pot) is also tiny (< 10 ug/L, against 147 000 dry). That fits a work-up that under-recovers aldehydes, but no other
  print exists to test it.
- **H1d (FFT closure): NON-DISCRIMINATING, as pre-registered.** Scaling the slow second-order rate of Munch T5 expt 1
  (furfural + Cys) to the pot's Cys 33 mM, a time-mean free furfural of ~2.9 mM (~28 mg / 100 mL) would be needed to
  make 12.1 ug FFT. Even with all of the Cys counted as H2S (the expt 3 rate), ~0.13 mM is needed. The printed finals
  are 0.0070 mM (Hofmann) and 0.00054 mM (Schieberle), so the shortfall is 19-410x and 250-5300x. Under either print,
  free furfural at the printed level cannot be the main FFT source unless the furfural pool turns over fast (that is
  the engine's picture; see CF-065). Closure cannot choose between the prints.

**Revised P (cluster 4's values in brackets):**
- E1, Schieberle under-recovered: 0.50 (0.45).
- E2, a slip: now only on Schieberle's side (52 printed per pot rather than per 100 mmol), 0.20 (0.25, both sides).
- Hofmann column high: 0.12 (0.15).
- Run-to-run (separate heatings): 0.13 (0.15).
- Other: 0.05.

None reaches 0.6, so the entry **stays open (P 0.85 that it remains unresolved on corpus evidence)**. P(Hofmann's
67.5 is within 2x of the truth) ~0.7; P(Schieberle's 5.2 is) ~0.2.

**Remedy change (recommended):** do NOT move `hofmann_ribose_FUR_insitu` to the geometric mean (18.7 ug). With this
much belief on one print, averaging puts the target where little probability lies, and OPERATING_PROCEDURE sec. 4
forbids averaging conflicting laboratories. Keep the row (67.5 ug, sigma_log 0.5: Schieberle's print already sits at
2.2 sigma). Record the 5.2 ug alternative in the choices register with weight 0.2, to run as an envelope scenario
when the envelope next runs. Expected effect on the shipped scorecard: **none**. At the shipped parameters the row is
+0.29 dex (FUR 131 ug / 100 mL, 0.58 sigma). For comparison, the previously pending geometric-mean remedy (sigma
0.6) would raise it to +0.85 dex (1.4 sigma), +1.7 to the row's sum r^2 before any refit.

## CF-065

- **H2a (order in H2S): LINEAR, supported.** n = log(550.8 / 130.1) / log 4 = **1.04 [0.90, 1.18]** (+/-10 % per
  cell). A fixed zero-order H2S sink fitted to the two points comes to 0.018 mmol (1.8 % of the high dose), which is
  negligible. So the second-order law that cluster 3 used to carry Hofmann's constant to Zhai's 10 mM holds within the
  Munich data, on a row outside the conflict. Saturation or a threshold in sulfide cannot explain the gap, and any
  order > 1 would make it worse (Zhai has less sulfide). Shipped engine: ratio 3.95 against 4.23 observed (-0.03 dex);
  3.90 on the 100 mL reading (CF-102).
- **H2b (fraction of furfural ending as FFT):** 0.114 % (0.25 mmol H2S) and 0.482 % (1 mmol) per furfural, or 0.456 %
  and 0.482 % per H2S charged. Both are below 1 %, so if the product survives, neither reagent is used up by FFT
  formation. The yield is then limited by the rate, not by the supply.
- **H2c (cysteine only as an H2S donor): FALSIFIED on the NF side, upheld on the furfural side.**
  - The engine's measured Cys -> H2S constant (Zheng, extrapolated to 145 C) gives R_pred = 0.028-0.037 when the fed
    H2S is conserved, and up to 0.056-0.075 when it is lost quickly (pH 5.0-5.34).
  - Observed R_FUR = 25.3 / 550.8 = **0.046** (in band).
  - Observed R_NF = 50.8 / 211.2 = **0.24**, 0.5-0.9 dex above the band. Both cells are in Hofmann 1998 Table 10
    (1.0 mmol each, 50 mL, pH 5.0). The prereg's "T4" names the table that reprints NF/H2S.
  - R_FUR / R_NF = 0.19, outside [0.5, 2]. Cysteine does more than release H2S on the NF -> MFT route: it may act as a
    reductant (Hofmann's Fig. 1 needs one) or protect MFT from loss.
- **The shipped engine against the Cys contrast (a new structural misfit).**
  - The engine gives R_FUR = 0.234, **+0.71 dex**: absolute 0.19 mol % against 0.022, +0.94 dex. On the 100 mL
    reading it is 0.194 (+0.62 dex).
  - It gives R_NF = 0.14, -0.24 dex.
  - Cause: in the engine's fed furfural + H2S pot, about all 20 mM of furfural becomes FFT within minutes
    (r_fur_fft_hs integral 20.5 mM, closure 5.6 % on a 2001-point grid), and about 99 % of that FFT is then destroyed
    (ch_fft_decay + ch_thiolate_loss_fft). That leaves 0.83 mol %. The fitted engine is therefore the "product loss in
    Hofmann's hold" world for FFT. That world predicts R_FUR ~0.2, because FFT formed slowly from a Cys trickle
    survives better than FFT formed early from fed H2S. The data give 0.046, the value expected for slow formation
    with a surviving product.
  - **So the Cys contrast is a test outside the conflict, and it counts against product loss as the FFT half of the
    CF-065 residual.** The way out is a true Cys -> H2S rate at 145 C about 5x below Zheng's extrapolation. That is
    possible: the parameter is flagged as an upper bound and is extrapolated 35 C. P 0.3.
  - (On a 201-point grid the flux audit does not close (70 %). Do not read flux integrals from coarse grids on
    fast-turnover pots.)

**Revised residual allocation:**

| half | quantification (Zhai RF = 1) | NH3 chemistry | product loss / Cys-relieved depression | other |
|---|---|---|---|---|
| FFT | 0.50 | 0.25 | 0.15 (was 0.25) | 0.10 |
| MFT | 0.40 | 0.20 | 0.30 (was 0.25) | 0.10 |

The MFT figure rises because Cys relieves something on the NF route that free H2S does not. The status stays
**open**, explained in part (conditions ~1.1 dex, now with the rate law checked on an outside row).

## Recommended data / engine changes (none applied)

1. Ledger: munch1997 Table 5 expts 1 and 2 enter as two within-study ratio rows against expt 3 (count once: expt 3 is
   hofmann1998_t3_fft_c5_h2s, CF-101). sigma ~0.06 dex comes from the +/-10 % per cell. Both rows are robust to
   CF-102's volume question (engine ratios move <= 0.09 dex between 50 and 100 mL).
   - (a) `munch1997_t5_fft_h2s_dose_ratio` (130.1 / 550.8). The shipped engine passes it (-0.03 dex). Fit or
     validate: the scorecard does not change.
   - (b) `munch1997_t5_fft_cys_over_h2s` (25.3 / 550.8). The shipped engine fails it by +0.71 dex (~12 sigma at
     0.06). Enter it as a HELD-OUT DIAGNOSTIC row, not a fit row, until the FFT-loss structure has been re-examined
     under the misfit protocol (mechanism step). Effect: one new failing diagnostic; the shipped fit is unchanged.
2. Do not change any constant on this evidence: the test points at structure (FFT turnover in the fed pot), and
   tuning k_fft_decay or k_cys_h2s to absorb it is forbidden. Next test (mechanism step): residual free furfural or
   H2S at 20 min in a fed furfural + H2S pot. The engine predicts ~0 for both; slow formation predicts > 90 %
   furfural left. Zhai Fig. 6a's furfural loss (80 % in 10 min at 120 C) points the engine's way, but that pot holds
   NH4+.
3. CF-016: keep the fit row; add the choices-register alternative (above).

## Expectations vs outcome

- n ~1: confirmed.
- R_pred in 0.02-0.3: confirmed (0.028-0.075).
- CF-016 stays open: confirmed.
- Not expected: the shipped engine's +0.71 dex failure on the furfural Cys contrast, and the opposite sign on the NF
  contrast.

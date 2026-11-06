# M0 checks batch 2: CK-23, CK-34, CK-42, CK-42b (MF-15 a) and CK-56

> **Misfit protocol, cheapest checks on seen data** ([OPERATING_PROCEDURE.md](../../docs/guides/OPERATING_PROCEDURE.md)
> sec. 4), for the misfit inventory (draft PR #82, branch `misfit-inventory`, `docs/reference/misfit_inventory.md`,
> secs. 3 and 5). 2026-11-06, from `origin/main` at `fdea5cf0`. Batch 1 (draft PRs #84, #86) and the argon Strecker
> checks (PRs #81, #86, #87) are cited, not repeated.
>
> **Rules kept.** The selection and every pass rule below were written and committed BEFORE any number of the check was
> computed (commit "pass rules before computing" on this branch). Every outcome is recorded, failures included.
> **Disclosure:** while scoping, the sources' own printed or dossier-derived numbers were read (for example the Wolf
> 1977a dossier's k at named conditions, CF-013's engine sink at 110 C, CF-065's dissolved-sulfide range), and rough
> mental arithmetic on them shaped which checks looked discriminating. No engine run and no scripted number of any
> check existed when the rules were committed.
> Constants come only from cited measurements, or are declared assumptions with their band. No constant is tuned to
> absorb a misfit. A check ships nothing: it nominates or retires candidates for a sec. 7 case. Nothing under `data/`,
> `src/`, the registers, cases, preregs or the inventory was edited. Engine use is read-only (`predict`-level runs and
> scripted parameter variants in [scripts/diagnostics/m0_checks_batch2.py](../../scripts/diagnostics/m0_checks_batch2.py)),
> one core, niced, threads pinned to 1.
>
> **Blinding.** Not read: the B59b hold-outs (the US 10,076,126 MFT shape is not used; see CK-23), van Seeventer 2001
> Table 1 (a declared FIT_HOLDOUT hold-out), hofmann2000b Figs. 2-3, the E7 / E13 / E14 TEST sets, the WUR samples.
> Not computed under any variant: the pots that open cases name as decisive (listed per check).

## 1. Selection, and why the others are out

**Rule.** The next highest-ranked checks of the inventory (sec. 5 ranking, then the unranked runnable checks) that are
(i) not owned by a case, (ii) not done, and (iii) whose spent data do not overlap the decisive data of any open case:
the M1b prereg (branch `m1b-prereg`; part (i)'s controls C-Mey, C-Hou, C-Zhou-ARP from CF-084), the M1c case (branch
`m1c-case`, with amendment A1), the M2 route design (`kinetic_core_m1m2_supply_routes_refit_prereg.md`, main), the M4,
M5 and M6 cases (main) and the M3b case (branch `m3b-case`). The M1 case (main) is not on the open list (EX-061 ran its
decisive comparison); its D-rows are still checked and named where a check touches them.

### Selected

| check | MF | rank | why selected | data spent | overlap test (iii) |
|---|---|---|---|---|---|
| **CK-23** | 01 (a) | 8 | highest-ranked eligible thiol-sink check (G3) | the US 10,076,126 Methods statement that the [2H3]-MFT / [2H2]-FFT standards were heated with the pot (CF-013 analysis; its MFT shape NOT used); Schieberle 2000 Table IV MFT (SEEN); Grosch 1992 Table VIII (SEEN, as an inferred in-pot λ); the B60 in-sample rows listed below | Schieberle T4 is M1's D1 (spent by EX-061), M2's D1 **guard** units (not claim units) and M1c's development signature (d): SEEN everywhere, decisive in no open case. M2's decisive pots (Cerny pot D, the Hofmann 145 C ribose + Cys and NF + Cys MFT pots of Rt-H, the Zhai 2023 JAFC 140 C pots) are **not computed** |
| **CK-34** | 13 (a) | 9 | G9 / G10; CF-065's proposed operator | Hofmann 1998 T3 / T4 / T6 / T7 fed-sulfide rows (B60 in-sample); Munch 1997 Table 5 (SEEN) | Hofmann T10's NF + Cys cell is M2's Rt-H denominator: **excluded** (so the T10 cysteine-over-H2S NF contrast is not computed). M4's D-rows use Tables 1, 2, 5, 9 (not touched). Whitfield 1999 (M1 D6) excluded |
| **CK-42** | 15 (a) | 10 | the Th x cys / Th x thi levers; engine-free | Grosch 1992 Tables II, III, VI, VIII (SEEN); Zheng & Ho 1994's measured k_cys_h2s (a shipped constant) | no open case names Grosch as decisive. **Flag:** Grosch's three atlas cells are v2 TEST(E18) contrasts F27#15-17; they were already scored (W1.1 decision set) and engine-run (batch 1, CK-10), so this engine-free check re-spends them at no new loss. Any E18 case must list them as SEEN |
| **CK-42b** | 15 (a) | new (batch 1) | the thiamine route's missing cysteine term, found by CK-10 in batch 1 | as CK-42, plus Arnold 1971's text (SEEN, CF-089) | as CK-42 |
| **CK-56** | 27 (a) | unranked (batch J) | engine-free; no case owns glycation's barrier | Wolf 1977a secs. 3.2-3.3 (the 105 C order experiment's printed k; the Case 1 / 2 response surfaces; SEEN); Wolf 1977 (companion) pass 1 as a side reading | Wolf 1977a is M4 sec. 2.7 development evidence, not decisive |

### Excluded

| check | rank | reason |
|---|---|---|
| CK-06, CK-01 + CK-50, CK-02 | 1, 2, 4 | done (batch 1, PR #86) |
| CK-12 | 3 | spends Hou 2026 Fig. 3, Zhou 2023 Fig. 2 and the Meynier pH statement: M1b part (i)'s controls C-Hou, C-Zhou-ARP, C-Mey (CF-084) |
| CK-43 (+ CK-44) | 5 | done / parked (PRs #81, #86, #87) |
| CK-46 | 6 | spends Piloty 1979's pH-flat amine consumption: M6 D3 |
| **CK-26 (+ CK-27)** | 7 | **its datum is M2's.** The Zhai 2023 JAFC 16-comparison T-order includes the 140 C pots, whose MFT shapes are M2's decisive units T-140 / T-140c (prereg sec. 5); running per-route barriers on them computes M2's temperature-law claim (E11) on its decisive pots. Cerny pot D is M2's Rt-D. CK-27 is M2's own sec. 5 requirement ("the T verdict requires the Gigl-ceiling profile to agree"). Without the 140 C cells and pot D, what is left (the 100 / 120 C cells) is B60 in-sample and cannot discriminate. **Proposed inventory correction:** mark CK-26 and CK-27 "owned by M2" |
| CK-13 | -- | Hou 2026 (M1b C-Hou) |
| CK-15 | -- | Hofmann 1998 T2 / T5: M4 D1 / D2; needs CK-12's state |
| CK-22 | -- | Xu 2007 is keyed held-out for M4 item 7 |
| CK-28 | -- | Zhou 2023 Table 1 / Zhang 2024 shares are M1's D3 / D4 (part C, never run); held back while M1's C3 is unrun |
| CK-29, CK-37 | -- | Elmore 2002 / Yiltirak 2026: M5 D5; US 8,920,862 a B59b hold-out; Xu 2011 E12 TEST |
| CK-45 | -- | Hofmann's fed MGO + Phe rate is M6 D4's quantity |
| CK-52 | -- | Leahy's ratios are M6 D6 |
| CK-16, -17, -19, -25, -33, -39, -40, -41, -51, -54 | -- | owned by a case (inventory) |
| CK-18, -21, -35, -47, -48, -57 - -59, -62 | -- | lower-ranked, or waiting on a batch-1 owner decision (CK-21 on MF-28; CK-57 - -59 on CK-04's pot corrections; CK-62 on CK-05's) |

---

## 2. The checks

Each check: pass rule (written before computing), constants with citations, outcome, data spent, nominated / retired
candidates, proposed CF text. The generic rule (inventory header): a check **EXPLAINS** if, with cited constants only,
the prediction falls within the unit's own noise band of the datum without moving in-sample rows out of theirs;
**sufficient only if \<constant\>** if it needs a fitted or unmeasured constant; **REFUTES** a candidate when the
candidate's best cited-constant prediction misses the datum by more than twice its noise band on the decisive side.

**In-sample "moves out" (all engine checks).** A B60 row moves out if its whitened residual (B60's sigma_log) grows in
magnitude by more than 1 and ends above 2 in magnitude.

### CK-23: the MFT sink bounded by the co-heated labelled standards (MF-01 (a))

**Pass rule (written before computing).**
- *Bound.* US 10,076,126 (Example 1, as quoted in CF-013's analysis and the us10076126 dossier sec. 2): the [2H3]-MFT
  and [2H2]-FFT standards were weighed into the vials before 0.5 - 4 h at 110 C, and the 4 h vial was quantified. No
  standard amount, recovery or LOD is printed. So the only bound is: the standard survived 240 min at a fraction
  f >= f_min, where f_min is the dynamic range a SIM quantitation tolerates. **Declared band (an assumption, not a
  measurement): f_min in [1e-6, 1e-2].** Then lambda_max(110 C) = ln(1 / f_min) / 240 min, i.e. 0.019 - 0.058 /min.
  The patent's MFT levels (its shape, a B59b hold-out) are not used.
- *Engine side.* The engine's first-order MFT loss at 110 C is k_mft_decay(110 C) plus the thiolate share of
  `k_thiolate_loss` (arithmetic on the shipped constants; no run on the patent pot). FFT likewise with `k_fft_decay`.
- *Verdict on the sink:* REFUTED as too fast at 110 C if the engine's loss exceeds lambda_max at the loosest f_min
  (1e-6); consistent if below lambda_max at the tightest (1e-2); inconclusive between.
- *Forward (MF-01).* Schieberle 2000 Table IV pot (ribose 100 + Cys 33 mM, 0.5 M phosphate pH 5, 100 C; the M1.0
  runner and settings), supply as shipped, with `k_mft_decay` scaled so its 110 C value equals lambda_max at
  f_min = 1e-6 and at 1e-2 (Ea held at B60's 101.25), and the limit k_mft_decay = 0. Units: MFT 360 / 30 and
  720 / 30 (dex). **Noise band 0.06 dex per ratio** (Table IV prints no SD; the same laboratory's SIDA triplicates
  "differ by <= 10 %" (Hofmann 1998), two points in quadrature). In-sample rows: the B60 rows of the 145 C Hofmann
  fed pots (ribose, 3-DP, furfural, NF + H2S; C2 + C3 at pH 3 / 5 / 7; thiamine; MGO + H2S), Zhou 2023 pH 7 (MFT,
  FFT, ratio, share, ACTZ) and the Kang / Zhai 100 / 120 C rows. **Not computed:** Hofmann T10 NF + Cys
  (`fed_nf_cys`, M2 Rt-H), the Hofmann 145 C ribose + Cys pot (Rt-H numerator), Whitfield 1999 (M1 D6).
- EXPLAINS MF-01 if at a sink inside the bound both ratios fall within 0.06 dex of the measured and no in-sample row
  moves out. **REFUTES candidate (a) as sufficient** if even k_mft_decay = 0 leaves 360 / 30 more than 0.12 dex
  below the measured. A cross-check, not a pass condition: Grosch 1992 Table VIII's inferred in-pot loss at 100 C
  (0.034 - 0.044 /min, assumption-laden) against the engine's k_mft_decay(100 C) and against the bound transported
  to 100 C under Ea in Gigl 2021's measured range (7 - 102 kJ/mol; central 60.2).

**Outcome.** *(pending)*

### CK-34: an H2S vessel partition on the fed-sulfide rows (MF-13 (a))

**Pass rule (written before computing).**
- *Operator (a diagnostic, not an engine change).* The engine charges every fed sulfide as fully dissolved. Under fast
  gas-liquid equilibrium in a closed vessel only phi = 1 / (1 + Kaw(T) alpha0 Vg / Vl) of the sulfide is in the
  liquid, where alpha0 is the neutral H2S fraction of the dissolved sulfide. Applied by scaling the k_ref of every
  reaction with H2S as a reactant by phi (r_fur_fft, r_fur_fft_hs, r_nf_mft, r_nf_mp3p, r_nf_shet, r_ddp_mft,
  r_ddp_mft_hs, r_mgo_mp, r_dmhf_h2s, ch_redox_mp3p, r_h2s_loss; sensitivity: r_h2s_loss unscaled).
- *Constants.* Kaw(25 C) 0.403 and delta_H_sol -17.5 kJ/mol (`henry_constants.yml`, Sander 5.0.0), van 't Hoff to
  the pot temperature; pK1(T) from Suleimenov & Seward 1997 (`parameters_sulfur.h2s_pk1_suleimenov_seward`) at the
  run's in-situ pH. Vg / Vl = 3 for the Hofmann fed pots (50 mL in the 200 mL Roth Type II autoclave,
  hofmann1998_reconciliation dossier); Munch 1997 (same autoclave) 3 (50 mL, P 0.7) and 1 (100 mL, the Methods'
  generic sentence), both reported.
- *Units.* (i) Munch T5 cysteine-over-H2S FFT contrast (expt 1 / expt 3, 1 mmol each in the furfural pot): measured
  25.3 / 550.8 = 0.046. (ii) Munch T5 H2S dose ratio (expt 3 / expt 2, 1 vs 0.25 mmol): 4.23. **Noise band 0.10 dex
  each** (CF-065's declared sigma >= 0.1; the same lab's <= 10 % per value gives 0.06). In-sample rows: the B60
  fed-sulfide rows (`fed_ribose_h2s_FFT`, `fed_ribose_h2s_MFT`, `fed_tdp_h2s_FFT`, `fed_furfural_h2s_FFT`,
  `fed_nf_h2s_MFT`, `fed_mgo_h2s_1to1_MP`, `fed_mgo_h2s_1to2_MP`). **Not computed:** `fed_nf_cys` (M2 Rt-H), so
  Hofmann T10's NF cysteine-over-H2S contrast is not scored here; Whitfield 1999 (M1 D6).
- EXPLAINS MF-13 if, at shipped constants with the operator, the contrast falls within 0.10 dex of 0.046, the dose
  ratio within 0.10 dex of 4.23, and no in-sample row moves out. "Sufficient only if" if it needs k_fur_fft_hs /
  k_nf_mft rescaled. REFUTES candidate (a) if the best prediction (over the Vg / Vl readings and the r_h2s_loss
  sensitivity) misses the contrast by more than 0.20 dex on the high side.
- *Reported, not a pass condition:* the factor on `k_fur_fft_hs` (with `k_fur_fft`) and on `k_nf_mft` that would
  restore `fed_furfural_h2s_FFT` and `fed_nf_h2s_MFT` under the operator (the bias CF-065 asks for; a diagnostic
  bisection, not a proposed constant).

**Outcome.** *(pending)*

### CK-42: Grosch's cysteine-dose optimum, engine-free (MF-15 (a))

**Pass rule (written before computing).** The candidate: "Cys -> H2S at the measured cysteine rate, plus a
Cys-dose-dependent MFT capture", tested on Grosch 1992 (thiamine 60 uM, ribose 1.34 mM, 0.2 M pyrophosphate pH 5.7,
2 h reflux, SIDA after boiling).
- *Rising limb (H2S part).* The sulfide cysteine releases in 2 h at 100 C at Zheng & Ho 1994's measured rate
  (`k_cys_h2s`, Ea band 123 - 135 kJ/mol) is compared with the exogenous sulfide of Table VI, which raised
  thiamine-MFT only x3 (15 +- 7 vs 4.7 +- 1.2 ug) at 26.8 mmol H2S flushed through 500 mL; the dissolved part of that
  dose is unstated, declared band 0.1 - 100 % of the charge. Table II no. 4 (thiamine + Cys, **no ribose**) gives x7
  (33 +- 3.7). The H2S part EXPLAINS the rise only if the cysteine-derived sulfide reaches the Table VI arm's
  dissolved sulfide at some point of the bands; it is REFUTED if the rate needed exceeds the measured by > 10x
  (1 dex) at every point of the bands.
- *Falling limb (capture part).* Table III (thiamine 30 umol): MFT 19 ug at 325 umol Cys falls to 2.0 at 2720 umol.
  Under constant formation and a first-order in-pot loss, MFT ~ F / (lambda0 + k_c c^n). The first-order capture
  (n = 1) sized by the 325 -> 2720 fall predicts lambda(0.65 mM Cys); Grosch Table VIII's inferred lambda with Cys
  0.65 mM is 0.034 (0.026 - 0.041) /min and without Cys 0.044 (0.035 - 0.053). n = 1 is REFUTED if its predicted
  lambda(0.65 mM) lies above the Table VIII band's upper end (0.041); n = 2 is recorded the same way. The high-dose
  log slope d ln MFT / d ln c over 680 -> 2720 umol is printed for both thiamine series.
- Noise: Grosch's printed SDs (Table II: 0.4 - 3.7 ug, n = 2 - 4; ~0.05 - 0.11 dex); Table III prints none (the
  Table II SDs transfer).

**Outcome.** *(pending)*

### CK-42b: the form of the thiamine route's cysteine term (MF-15 (a), from batch 1's CK-10)

**Pass rule (written before computing).**
- *Ratio or concentration.* A stoichiometric term (cysteine acting through a Cys : thiamine ratio, e.g. an adduct with
  a thiamine-derived intermediate) moves the dose optimum with the thiamine charge: at thiamine 6 umol (Table III
  second row) the optimum should sit at 0.2x the 30 umol series' Cys (34 - 65 umol). A concentration term (a
  reductant, catalyst or medium effect) keeps it at 170 - 340 umol. Discriminates if the 6 umol series' optimum,
  bracketed by its printed points (170 / 340 / 680 / 1360 umol), lies outside one law's range by more than one grid
  step.
- *A total-thiol concentration law, predicted out of construction.* Build the law from Table III (thiamine 30 umol,
  ribose 670 umol; log-log interpolation in thiol concentration) and predict Table II no. 7 (thiamine + GSH 0.46 mM +
  ribose; 23 +- 1.2 ug, n 3) and no. 8 (thiamine + GSH 0.46 + Cys 0.65 mM + ribose, total thiol 1.11 mM; 9.0 +- 2.8,
  n 4). PASSES if each prediction is within 2 x sqrt((SD / mean / ln 10)^2 + 0.05^2) dex of the measured (0.05 dex
  the interpolation's declared read error); FAILS otherwise. No. 6 (thiamine + Cys 0.65 mM + ribose, 19 +- 2.5)
  against Table III's 325 umol cell (19) is printed as a same-design check.
- *pH artefact.* Cysteine is charged as Cys.HCl: the acid load at 2720 umol in 500 mL of 0.2 M pyrophosphate (pKa
  6.6 for the third proton) shifts the pH by delta-pH (arithmetic). Arnold 1971 (Fig. 2): thiamine loss at ~100 C is
  <= 2 % below pH 6 and 14 - 17 % at pH 7, i.e. at most ~1 dex per pH unit. The pH artefact is RETIRED as the cause
  of the x9.5 fall if |delta-pH| x 1 dex/pH < 0.5 dex (it can then give at most x3 of the x9.5).

**Outcome.** *(pending)*

### CK-56: Wolf 1977a's barrier without the plateau-prone end point (MF-27 (a))

**Pass rule (written before computing).** Candidate (a) of MF-27: the 80 C / 360 min end points are plateau-capped, so
k(80 C) is biased low and Wolf's two-temperature E (32.9 kcal/mol = 137.7 kJ/mol) is inflated. A plateau-free E uses two
phase-one rates only: the printed first-order slopes of the 105 C order experiment (0.024 and 0.026 /min; 18 min; pH 7,
a_w 0.65, 4 % glucose) and the response surface's k at 125 C (1.5 min end points, 36 - 45 % loss, inside phase one) at
the same composition (Case 1 and Case 2; glucose coded +1.68 (3.68 %, the design's edge) and +2.0 (4 %)).
- E(105 -> 125) = R ln(k125 / k105) / (1 / 378.15 - 1 / 398.15).
- Noise: k125 carries the surface's residual SD (k_100: 0.003 /min, ~0.1 dex; E: 2.0 kcal/mol, Table 4), propagated;
  k105 its two printed values. **The band is the larger of Wolf's residual SD of E (8.4 kJ/mol) and the propagated SD.**
- EXPLAINS MF-27 through (a) if the plateau-free E lies within one band of the engine's 96.8 kJ/mol. REFUTES (a) if it
  exceeds 96.8 by more than two bands. Side reading (not a pass condition): the companion paper's 125 C / 2 min pass 1
  (9.6 -> 5.0 mg/g) as a second k125.

**Outcome.** *(pending)*

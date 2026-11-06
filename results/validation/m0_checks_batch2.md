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

## Outcomes

| check | misfit | outcome in one line |
|---|---|---|
| CK-23 | MF-01 (a) | Shipped MFT sink **refuted** as too fast at 110 C (1.9 - 5.6x over the co-heated-standard bound). No pot-independent sink inside the bound fits Schieberle and keeps B60's rows. Candidate (a) survives only as a **cysteine-slowed loss** (M1's C2) |
| CK-34 | MF-13 (a) | The H2S partition (phi 0.14) moves the Munch contrast only +0.71 -> +0.51 dex: **refuted** as the explanation. The fitted furfural constants absorb it (bias +0.84 dex); `k_nf_mft`'s bias is +0.23 |
| CK-42 | MF-15 (a) | Cysteine acting through H2S **refuted** (needs >= 17x Zheng & Ho's rate). A first-order cysteine capture is refuted; second order survives at the band edge |
| CK-42b | MF-15 (a) | The thiamine route's cysteine term is a **total-thiol concentration law**, not a ratio law: it predicts Grosch's GSH pots out of construction. The pH artefact is retired |
| CK-56 | MF-27 (a) | Plateau-inflated barrier **refuted**: the plateau-free E(105 -> 125 C) is 154 - 196 kJ/mol, and the apparent barrier rises with temperature |

No misfit is explained. Three candidates are retired (CK-34 (a), CK-42's H2S route, CK-56 (a)) and one is refuted as
built (the shipped MFT sink at 110 C). Two candidates are nominated with a form (CK-23: C2; CK-42b: g(c_SH)).

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

**Outcome: the shipped MFT sink is REFUTED as too fast at 110 C (by 1.9 - 5.6x). MF-01 is NOT explained by a
pot-independent sink inside the bound, and candidate (a) is NOT refuted as sufficient. It survives only as a sink that
is slower in cysteine-rich pots (M1's C2 form), or with a barrier outside Gigl's measured range.**

*The bound.* A trace MFT / FFT "standard" alone in 0.2 M phosphate pH 5.5 at 110 C (the patent's conditions, no
precursors; in-situ pH 5.64) loses MFT at **0.107 /min** in the shipped engine. Of that, `k_mft_decay` is 0.095 and
the rest (thiolate loss) 0.012. It survives 240 min at 7e-12 of its charge. FFT loses 0.022 /min.

| f_min (declared) | lambda_max(110 C) | engine MFT / bound | engine FFT / bound | bound at 100 C, shipped Ea | bound at 100 C, Gigl Ea 7.2 / 60.2 / 102.4 |
|---|---|---|---|---|---|
| 1e-6 | 0.058 /min | **1.86** | 0.38 | 0.025 | 0.054 / 0.035 / 0.024 |
| 1e-2 | 0.019 /min | **5.58** | 1.15 | 0.008 | 0.018 / 0.012 / 0.008 |

- **MFT sink: refuted as too fast** at both corners. The 2026-10-29 analysis (CF-013's test line) said ">= 5x". That
  figure leaned on the patent's 4 h level, i.e. its B59b-held-out shape. On the standard's survival alone it is
  >= 1.9x (loose) to >= 5.6x (tight).
- FFT sink: consistent at the loose corner, inconclusive at the tight one.

*Grosch cross-check.* The engine's MFT loss at 100 C, pH 5.7 (standard alone) is 0.048 /min. Grosch's inferred in-pot
loss is 0.044 (0.035 - 0.053) /min without cysteine and 0.034 (0.026 - 0.041) /min with 0.65 mM. So the engine
matches Grosch's cysteine-free pot. The patent bound, transported to 100 C at the shipped Ea (0.025 / 0.008), sits
below both Grosch values. One pot-independent sink fits Grosch and the patent only if Ea <= 62.6 kJ/mol at the loose
corner, and never at the tight corner.

*Forward on Schieberle 2000 (100 C, 33 mM Cys; M1.0's runner reproduces -1.094 / -1.421 at shipped constants):*

| sink | k_mft_decay at 100 C / 145 C (/min) | MFT 360/30 resid (dex) | MFT 720/30 resid (dex) | MFT peak (min) | in-sample rows moved out |
|---|---|---|---|---|---|
| shipped | 0.041 / 1.36 | -1.09 | -1.42 | 158 | -- |
| bound, f_min 1e-6 (Ea held) | 0.020 / 0.65 | -0.83 | -1.17 | 190 | none (but `kang_120C_MFT` 2.51 -> 3.15) |
| bound, f_min 1e-2 (Ea held) | 0.0031 / 0.10 | -0.23 | -0.38 | 391 | 5: `fed_nf_h2s_MFT` 0.64 -> 3.01, `fed_c2c3_MFT_pH3` -> 3.64, `kang_120C_MFT` -> 4.11, `kang_100C_MFT_over_FFT` -> 2.75, `zhai_MFT_fold_120_over_100` -> 3.90 |
| k_mft_decay = 0 | 0 / 0 | **-0.03** | **+0.06** | >= 720 | 7 (the five above, plus `fed_c2c3_MFT`, `fed_thiamine_MFT`) |
| bound 1e-6 with Gigl's central Ea 60.2 | 0.028 / 0.22 | -0.96 | -1.30 | 173 | 3 |

- No sink inside the bound fits Schieberle (both ratios within 0.06 dex) while keeping the in-sample rows. The
  zero-sink limit comes closest: 360/30 is inside the band, 720/30 is 0.004 dex outside it but inside two bands. It
  moves seven B60 rows out, all of them cysteine-free 145 C fed pots or the fed-TTCA ladder.
- FFT's ratios are unchanged (+0.36 / +0.27) in every variant, since FFT has its own sink.
- So MF-01 is **not explained**. Candidate (a) is **not refuted as sufficient**: under the shipped supply, the late
  rise appears once the 100 C loss falls to <~0.003 /min.

*What the pattern says.*
- The rows that need the fast sink are cysteine-free: H2S-, C2 + C3- and thiamine-fed pots at 145 C, with
  k_mft_decay ~1.4 /min.
- Grosch's dilute pot (0 - 0.65 mM Cys, 100 C) sits at ~0.04 /min, slightly slower with cysteine.
- Schieberle (33 mM Cys, 100 C) needs <~0.003 /min under the shipped supply.
- The patent (250 mM Cys, 110 C) bounds it at <= 0.019 - 0.058 /min.
- A temperature-only reading (0.003 at 100 C, 1.36 at 145 C) needs an apparent Ea of ~176 kJ/mol. That is outside
  Gigl's measured 7 - 102 range, above the engine's 102 ceiling, and against Grosch's 100 C value. So it is not
  admissible with cited constants.
- What remains is a **loss that free cysteine slows** (M1's combination C2, EX-050's form), which ordered by cysteine
  fits all four. This is conditional on the shipped (instant) supply. M1 found a real reservoir (C1, 3 - 5x too slow as
  built). With a slower supply the sink needed at 100 C changes, so the C2 reading is a nomination, not a size.

**Nominated:** for M1c / the joint refit, M1's C2 (a cysteine-slowed thiol loss), judged jointly with the supply
structure. **Retired:** a frozen pot-independent `k_mft_decay` as the sole answer (it cannot meet the patent bound, the
Schieberle shape and the 145 C fed rows at one Ea in Gigl's range). **Decisive data** (not on disk):
- CF-013's designed protocol: [2H3]-MFT charged before heating in sister vials;
- at two cysteine levels (e.g. 0 and >= 30 mM, same pentose, 100 and 120 C), which separates C2 from a
  temperature-only law;
- atlas E12 (labelled MFT spiked mid-cook).

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

**Outcome: candidate (a) REFUTED as the explanation of the cysteine-over-H2S contrast.** The partition is real
arithmetic, but it is absorbed by the fitted furfural constants and leaves the contrast where it was.

- *The operator.* At 145 C, Kaw(H2S) = 2.18 and pK1 = 6.49 (Suleimenov). At the in-situ pH of 5.34 the neutral
  fraction is 0.933. So **phi = 0.141** at Vg / Vl 3 and 0.330 at Vg / Vl 1 (CF-065's estimate was 15 - 24 %).
- *Baseline reproduced.* With no partition the contrast is 0.234 (+0.71 dex) and the dose ratio 3.95 (-0.03 dex),
  as CF-065 records. The guard rows reproduce B60's residuals to 7e-11.

| reading | contrast Cys / H2S (measured 0.046) | residual (dex) | dose ratio (measured 4.23) | residual (dex) |
|---|---|---|---|---|
| no partition, 50 mL | 0.234 | +0.71 | 3.95 | -0.03 |
| partition, 50 mL (Vg / Vl 3) | 0.157 | **+0.53** | 3.70 | -0.06 |
| partition, 50 mL, r_h2s_loss unscaled | 0.148 | **+0.51** | 3.71 | -0.06 |
| partition, 100 mL (Vg / Vl 1) | 0.154 | +0.52 | 3.74 | -0.05 |
| partition, 100 mL, r_h2s_loss unscaled | 0.147 | +0.51 | 3.75 | -0.05 |
| partition, 50 mL, furfural constants restored (+0.84 dex) | 0.237 | +0.71 | 3.94 | -0.03 |

- The best reading misses the contrast by +0.51 dex, against the 0.20 dex refutation line. The dose ratio stays
  within its band in every reading.
- *Why it cannot work.* Furfural is fully consumed within 20 min in both arms (residual furfural <= 0.04 % of the
  charge). The H2S arm stays saturated: FFT moves -0.02 dex under phi. The cysteine arm loses only ~0.19 dex, because
  its furfural goes mostly to `r_fur_decay` either way.
- *The bias.* Restoring `fed_furfural_h2s_FFT` under the operator needs `k_fur_fft` and `k_fur_fft_hs` x10^0.84,
  i.e. ~1 / phi: the shipped constants already carry the partition implicitly. With them restored the contrast
  returns to +0.71. Restoring `fed_nf_h2s_MFT` needs `k_nf_mft` x10^0.23 only. That row is not saturated, but it
  sits near its formation-sink optimum: the scan's maximum, 0.382 mol %, is barely above the row's B60 value of
  0.367. So the "0.6 - 0.8 dex bias" CF-065 estimated holds for the furfural constants (0.84), not for `k_nf_mft`
  (0.23).
- *In-sample rows.* None moves out by the rule. Three shift by more than 1 sigma: `fed_ribose_h2s_MFT` -0.83 ->
  -1.84, and `fed_mgo_h2s_1to1_MP` / `_1to2_MP` +0.06 / -0.06 -> -1.72 / -1.84. The operator is not free: it would
  have to be refitted with `k_mgo_mp` and the ribose + H2S route.

**Retired:** an H2S vessel partition as the cause of MF-13's contrast (it remains a correct conditions term whose size
the fitted constants absorb; carry it in any refit that mixes vessels). **Nominated:** (b) the furfural lifetime in
the cysteine arm (`k_fur_decay`, clipped on its upper bound at B60) and (c) FFT survival in the hold. **Decisive
data:** CF-065's deciders (ii), residual free furfural and H2S at 20 min in Hofmann's fed furfural + H2S pot (the
engine says ~0 for both), and (iii), authentic FFT heated alone. Neither is on disk.

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

**Outcome: the H2S half is REFUTED; the capture half survives only as a term of order >= 2 in cysteine (no cited
constant).**

*Rising limb.* Cysteine 0.65 mM releases **2.0 - 3.1 uM** sulfide in 2 h at 100 C at Zheng & Ho's rate (Ea 135 /
133 / 123 kJ/mol; 0.31 - 0.47 % of the charge). Table VI's arm charged 26.8 mmol H2S over 500 mL (53.6 mM if all
dissolved). Even if only 0.1 % of it stayed in solution (54 uM), the cysteine-derived sulfide would need **>= 17x**
the measured rate to match it (174x at 1 %, 1700x at 10 %). That exogenous dose gave x3.2 (15 +- 7 vs 4.7 +- 1.2 ug),
while cysteine gave **x7.0** (33 +- 3.7, Table II no. 4, with no ribose). So cysteine does more for thiamine-MFT
than >= 17x its own sulfide release does. A dicarbonyl-assisted (Strecker-type) release in no. 4 is capped by the
thiamine charge itself (60 uM, the only carbonyl source there): even complete conversion only matches the most
pessimistic Table VI arm, which gave x3, not x7. **Refuted at every point of the bands (needed factor >= 17 > 10).**

*Falling limb.* Under constant formation, MFT ~ F / (lambda0 + k_c c^n), with the 325 -> 2720 umol fall of x9.5
(2 SD band x5.3 - x16.9, Table II's relative SDs transferred):

| law | central (lambda0 0.044, fall 9.5) | across the bands (lambda0 0.035 - 0.053; fall 5.3 - 16.9) |
|---|---|---|
| n = 1 (first order in Cys) | **cannot** make the fall: its ceiling is the cysteine ratio, x8.37 | refuted everywhere (lambda(0.65 mM) >= 0.084 /min where feasible, against Table VIII's <= 0.041) |
| n = 2 | lambda(0.65 mM) 0.050 /min: above 0.041, refuted at the central values | survives at lambda0 0.035 for every fall in the band (lambda(0.65 mM) 0.037 - 0.040) |

The high-dose log slopes are -1.16 (thiamine 30 umol, 680 -> 2720 umol) and -1.68 (thiamine 6 umol, 680 -> 1360).
Table VIII's own reading points the other way from any additive capture: the in-pot loss with 0.65 mM cysteine
(0.034 /min) is *lower* than without (0.044), i.e. at low dose cysteine protects MFT (as CF-065's Munich cluster
found on the NF route: "cysteine does more than release H2S").

**Nominated** (for an MF-15 case): (i) a non-sulfide cysteine term on the thiamine route at low dose (CK-42b gives
its form); (ii) a high-dose depression of order >= 2 in cysteine, or a cysteine inhibition of thiamine degradation
(Dwivedi & Arnold 1971, qualitative and mainly at pH 7). **Retired:** cysteine acting on thiamine-MFT through its H2S
release. **Sufficient only if** a second-order depression constant k_c ~ 0.006 - 0.015 mM^-2 min^-1 at 100 C (no
source).

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

**Outcome: a concentration law, not a ratio law; a total-thiol concentration law PASSES out of construction; the pH
artefact is RETIRED.**

- *Ratio or concentration.* At thiamine 6 umol the printed series is 2.1 / 2.8 / 1.6 / 0.5 ug at 170 / 340 / 680 /
  1360 umol Cys: **the optimum sits at 340 umol**, as at 30 umol (within 10 % of the maximum at 170 - 325 umol). The
  ratio law puts it at 34 - 65 umol, below the first printed point. There the series would already be falling from
  170 to 340 umol (Cys : thiamine 28 -> 57; on the 30 umol curve that stretch, 850 -> 1700 umol, falls about
  0.35 dex). It rises instead (+0.12 dex). The ratio law misses the optimum by more than two grid steps and is
  **refuted**; the concentration law (170 - 340 umol) holds.
- *Total-thiol law* (nodes from Table III, thiamine 30 umol; the 325 / 340 umol near-duplicates 19 / 14 ug merged
  to 16.5 at 0.665 mM, declared):

  | Table II row | total thiol (mM) | predicted (ug) | measured (ug) | residual (dex) | tolerance (dex) | |
  |---|---|---|---|---|---|---|
  | no. 7, thiamine + GSH + ribose | 0.46 | 18.0 | 23 +- 1.2 | -0.106 | 0.110 | pass (marginal) |
  | no. 8, thiamine + GSH + Cys + ribose | 1.11 | 11.5 | 9.0 +- 2.8 | +0.108 | 0.288 | pass |
  | no. 6, thiamine + Cys + ribose (same design as Table III's 325 umol cell) | 0.65 | 16.6 | 19 +- 2.5 | -0.058 | 0.152 | pass |

  A thiol counts as a thiol: GSH stands in for cysteine at equal concentration (no. 7 slightly above, so GSH may be a
  little more effective per mol). Two thiols together fall on the same curve (no. 8). The shipped engine gives GSH
  only `r_gsh_h2s` (inventory MF-15 (c)), and so predicts no such effect.
- *pH artefact.* Cys.HCl at 2720 umol adds 5.44 mM HCl to 0.2 M pyrophosphate: pH 5.70 -> 5.57 (**delta-pH -0.13**,
  Henderson-Hasselbalch on pKa3 6.6, declared). At Arnold's at most ~1 dex per pH unit this moves thiamine-MFT by
  <= 0.13 dex of the 0.98 dex fall. **Retired.**

**Nominated candidate** (for a sec. 7 case on MF-15): the thiamine -> MFT steps carry a factor g(c_SH), a function of
the total free-thiol concentration only (not the thiol : thiamine ratio). It rises x2.7 - 3.5 from 0 to 0.17 - 0.34
mM, is flat to ~0.65 mM, and falls about as c^-1.2 to c^-1.7 above ~0.7 mM. GSH counts as a thiol. With CK-42: not
through H2S. Mechanistic readings to separate (both [P], Grosch's "SH-group replacement" equilibrium): a thiol
reductant on a thiamine-derived disulfide or oxofuranone, versus thiol-catalysed ring opening. **Decisive data:** none
on disk (Grosch is spent). CF-089's designed series, plus a non-amine thiol arm (e.g. N-acetylcysteine or
2-mercaptoethanol) at the optimum concentration and thiamine measured: this separates a thiol role from an amine
role.

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

**Outcome: candidate (a) REFUTED.** Removing the plateau-prone end point raises Wolf's barrier instead of lowering it.

| surface | glucose (coded) | surface E (kJ/mol) | k125 (/min) | E(105 -> 125), measured k105 0.024 - 0.026 (kJ/mol) | band (kJ/mol) | E(80 -> 105) with the surface's k80 (reported) |
|---|---|---|---|---|---|---|
| Case 1 | +1.68 (3.68 %) | 135.7 | 0.306 | **156.7** (154.3 - 159.3) | 14.5 | 120.8 |
| Case 1 | +2.0 (4 %) | 131.8 | 0.302 | **156.0** (153.6 - 158.6) | 14.1 | 114.6 |
| Case 2 | +1.68 | 143.9 | 0.377 | **169.8** | 14.2 | 125.5 |
| Case 2 | +2.0 | 152.4 | 0.572 | **196.0** | 13.0 | 121.5 |

- In every reading the plateau-free E exceeds the engine's 96.8 by 59 - 99 kJ/mol, i.e. by more than four bands.
  The side reading (the companion's 125 C pass 1, k >= 0.33 /min) gives 161 kJ/mol.
- **The sign is the opposite of (a).** The 80 -> 105 C segment (115 - 125 kJ/mol) is *lower* than the 105 -> 125 C
  segment (156 - 196). A plateau-capped 80 C point would make the low segment *higher*, so the curvature is not the
  plateau's. The apparent barrier rises with temperature.
- Bias directions checked. A come-up lag in the 1.5 min retort-pouch run at 125 C would make k125 low, so E is
  underestimated. An onset of the plateau inside 1.5 min at 125 C (36 - 58 % loss, phase one 30 - 60 %) works the same
  way. Both push E further from 96.8.
- The response surface reads k(105 C) 0.035 - 0.050 /min, against the order experiment's measured 0.024 - 0.026. That
  is 0.15 - 0.30 dex high, beside a residual SD of ~0.1 dex.

**Nominated:** (c) the water-phase glucose (dissolution, the glass, moisture shifts) changing with temperature. It
carries the observed signature: an apparent E that rises with temperature, i.e. non-Arrhenius in a low-moisture solid.
Also (b), a genuine barrier difference between free glycine (Martins 2005, the engine's source) and protein-bound
lysine. (b) alone predicts one E across 80 - 125 C, which the two segments contradict, unless the 80 C point is
itself biased. **Retired:** (a). **Decisive data:** CF-100's soy or pea time course at >= 3 temperatures with moisture
printed. A glass-transition temperature (or a_w / moisture at each temperature) would separate (c) from (b).

---

## 3. What this batch spends

Datasets read here, now SEEN for any later case. A case that later names one of them as decisive must list it under
its item 2.

| dataset | spent on | already SEEN? |
|---|---|---|
| US 10,076,126 Example 1 Methods (standards before heating; no amounts) | CK-23's bound | yes (CF-013). Its MFT **shape stays a B59b hold-out, not used** (CF-013's register text that prints it was displayed while reading the register) |
| Schieberle 2000 Table IV, MFT and FFT ratios | CK-23 forward, under five sink variants | yes (M1 D1, EX-061; M2 D1 guards; M1c signature (d)) |
| Grosch 1992 Table VIII (the inferred in-pot loss) | CK-23 cross-check, CK-42 | yes (dossier, EX-050) |
| Grosch 1992 Tables II, III, VI | CK-42, CK-42b | yes (W1.1 decision set; batch 1 CK-10). **v2 TEST(E18) F27#15-17:** re-spent here, engine-free |
| Hofmann 1998 T3 / T4 / T6 / T7 fed-sulfide rows; Munch 1997 Table 5 | CK-34 | yes (B60 fit rows; CF-065 Munich cluster) |
| B60 in-sample rows (Hofmann fed pots, Kang / Zhai 100 / 120 C, van Seeventer conversion, Feng) | the guards of CK-23 and CK-34 | yes (fit rows) |
| Wolf 1977a secs. 3.2 - 3.3; Wolf 1977 pass 1 | CK-56 | yes (CF-096, CF-100) |
| Zheng & Ho 1994's k_cys_h2s; Kaw(H2S); Suleimenov pK1; Gigl 2021 barriers; Arnold 1971 text | constants | yes |

Engine outputs were computed under variants for:
- the Schieberle pot (the k_mft_decay variants);
- the Munch / Hofmann furfural pots (the partition);
- the guard systems listed in the script (`GUARD_SYSTEMS_CK23`, `GUARD_SYSTEMS_CK34`).

They were **not** computed for:
- `fed_nf_cys` or the 145 C ribose + Cys pot (M2 Rt-H);
- Cerny pot D;
- the Zhai 2023 JAFC 140 C pots;
- Whitfield 1999;
- the Zhou 2023 / Zhang 2024 pots;
- the US 10,076,126 pot itself.

## 4. Proposed CF texts and corrections (not applied)

1. **CF-013, tests, add** (`m0_checks_batch2`, CK-23):
   > The co-heated [2H3]-MFT standard bounds the in-pot MFT loss at 110 C to lambda <= ln(1 / f_min) / 240 min =
   > 0.019 - 0.058 /min for a declared f_min of 1e-2 - 1e-6. The patent prints no standard amount, recovery or LOD.
   > The shipped engine loses a lone standard at 0.107 /min (`k_mft_decay` 0.095), 1.9 - 5.6x too fast. The earlier
   > ">= 5x" used the patent's 4 h level (the B59b-held-out shape) and holds only at the tight corner. FFT
   > (0.022 /min) is within the bound.
   >
   > On Schieberle 2000 no pot-independent sink inside the bound fits the shape and keeps B60's 145 C fed rows. The
   > zero-sink limit gives 360/30 -0.03 and 720/30 +0.06 dex, but moves 7 rows out. The fast-sink rows are
   > cysteine-free and the slow-sink pots cysteine-rich (Grosch 0 / 0.65 mM: 0.044 / 0.034 /min; Schieberle 33 mM:
   > <~0.003 under the shipped supply; patent 250 mM: <= 0.019 - 0.058), which points to M1's C2 (a loss free
   > cysteine slows). A temperature-only law needs Ea ~176 kJ/mol, outside Gigl's measured 7 - 102.

2. **CF-065, tests, add** (CK-34):
   > H2S vessel partition (Kaw(145 C) 2.18, pK1 6.49, in-situ pH 5.34; phi 0.141 at Vg/Vl 3, 0.330 at 1) applied to
   > every H2S-consuming step: Munch's cysteine-over-H2S FFT contrast 0.234 -> 0.148 - 0.157 (+0.51 to +0.53 dex
   > against 0.046), the dose ratio 3.70 - 3.75 (passes). Furfural is fully consumed in both arms. Restoring
   > `fed_furfural_h2s_FFT` needs `k_fur_fft` / `k_fur_fft_hs` x10^0.84 (~1 / phi), which returns the contrast to
   > +0.71. `k_nf_mft`'s bias is +0.23 dex (not 0.6 - 0.8). Remedy, amended: the partition is a correct conditions
   > term absorbed by the fitted furfural constants and does not explain the contrast. The contrast's candidates
   > are the furfural lifetime and FFT survival (deciders (ii), (iii)).

3. **CF-089, tests and remedy, add** (CK-42, CK-42b):
   > Engine-free on Grosch Tables II, III, VI, VIII:
   > - Cysteine's x7 on thiamine-MFT is not its H2S release. 0.65 mM Cys releases 2 - 3 uM sulfide in 2 h at Zheng &
   >   Ho's rate, >= 17x less than the Table VI arm, which gave only x3.2. Table II no. 4 has no ribose.
   > - The dose optimum sits at 170 - 340 umol cysteine at both 30 and 6 umol thiamine: a concentration law, not a
   >   Cys : thiamine ratio law.
   > - A total-thiol concentration law built on Table III predicts the GSH pots (no. 7 -0.11, no. 8 +0.11 dex) within
   >   their SDs.
   > - A first-order cysteine capture cannot make the high-dose fall. The pH shift from Cys.HCl is -0.13.
   >
   > Remedy: model the cysteine effect on the thiamine route as a factor g(total free thiol), with GSH counted, sized
   > by a designed series with a non-amine thiol arm and thiamine measured.

4. **CF-100, tests, add** (CK-56):
   > Explanation 1 (plateau-capped 80 C end points) refuted. The plateau-free pair, the measured 105 C phase-one k
   > (0.024 / 0.026) and the surface's 125 C k, gives E 154 - 196 kJ/mol, above the published 137.7 and far above the
   > engine's 96.8. The 80 -> 105 C segment is lower (115 - 125), so the apparent barrier rises with temperature, the
   > signature of explanation 3 (water-phase glucose / glass).

5. **Inventory corrections** (for the inventory's owner):
   - Mark CK-26 and CK-27 "owned by M2": the Zhai 140 C cells are T-140 / T-140c, pot D is Rt-D, and CK-27 is M2's
     sec. 5 Gigl-ceiling requirement.
   - Mark CK-28 "held while M1's C3 is unrun" (M1 D3 / D4).
   - Record CK-23 / 34 / 42 / 56 and the new CK-42b with this record's path.
   - In the *spends* column of CK-42, add "v2 TEST(E18) F27#15-17, already scored".

## 5. Owner-level calls this record raises

1. Whether **M1's C2** (a cysteine-slowed thiol loss) enters M1c's or the joint refit's rival set. CK-23 finds it the
   only sink reading consistent with the patent bound, Grosch, Schieberle and the 145 C fed rows under cited
   barriers.
2. The **Grosch E18 TEST status**: the v2 split lists F27#15-17 as E18 TEST, yet they were scored in W1.1 and run
   under the engine in batch 1. They should be declared spent in the split (a `docs/` edit, not made here).

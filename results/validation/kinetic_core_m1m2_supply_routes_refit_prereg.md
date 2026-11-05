# M1 + M2 pre-registration: thiol supply and MFT route shares, judged as combinations in one refit (v1 DRAFT)

**Status: v1 DRAFT, not frozen; design review pending.** Nothing in this file has been reviewed. No build, fit, probe
or score of any rival below exists. Freezing needs the design review (sec. 11 (a)), the answers, and a freeze commit.

*Plan steps M1 and M2 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, rows M1 and M2; sec. 7 items
1 - 9; the hard dependency "M1 and M2 in one refit", sec. 6). Predecessors: M1 (EX-061),
[kinetic_core_m1_supply_loss_combinations_prereg.md](kinetic_core_m1_supply_loss_combinations_prereg.md) (v3.3, frozen
at `68247583a7`), its [decisive comparison](kinetic_core_m1_decisive_comparison.md), its
[amendment A1](kinetic_core_m1_amendment_A1.md) and its [independent results review](kinetic_core_m1_results_review.md);
the E6 + E11 proposal [kinetic_core_mft_route_barriers_prereg.md](kinetic_core_mft_route_barriers_prereg.md) (2026-11-01,
never frozen) and its stage 0 (`origin/mft-route-barriers-2026-11-01`, ff041ee2, unmerged). Conflict: CF-107. Case:
[M1_thiol_loss_oxidant_ph_drift.md](../../docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md). There is no M2 case
file; secs. 1 - 6 of this file are M2's case items 1 - 5.*

- Written 2026-11-05 on branch `m2-refit-prereg` from main at 46d538c4. No engine was run to write it; the only numbers
  computed here are arithmetic on printed or keyed values, each marked "(arithmetic)".
- Engine state: main at 46d538c4, shipped calibration B60, inline pH coupling as in M1. The registers
  (`experiment_register.yml`, `choices_register.yml`, `conflicts_register.yml`) and `STATUS_AND_PLAN.md` are not edited
  here (another session holds them); the entries this design needs are listed in sec. 16.3 for that session.

## 0. Seen by the author

**The expectations in sec. 15 are not blind.** Read while writing:

- **M1, all of it:** the frozen design; R0's and C1's fit records and carried optima (R0 cost 21.883 on B60's 52 rows;
  C1 186.75 = 58.00 on the 52 rows + 128.75 on the route block; C1's `k_fur_decay` -9.98 and `k_osone_decay` -9.98
  against R0's +0.18 / -1.35; `Ea_decay_thiol_sink` at 102.0 in both; C1's `ea_form` at 40.0, `k_ttca_cys` -1.23); the
  decisive comparison (every unit x reading x form); the results review's exact R0 Cys (about -746 under the alternative
  floor) and its residual table (R0: Cys +1.68 ... +0.38, PENT +3.7 ... +15, MFT +0.68 / +1.73 / +2.07, FFT -0.05 /
  -0.16 / +0.05; C1: Cys -0.21 ... -1.34, PENT -0.24 ... -1.08, MFT +2.04 / +1.08 / +1.02, FFT +1.94 / +0.64 / +0.20
  dex); the measured free-reactant decay constants (Cys 0.072, PENT 0.060 /min at 120 C) against C1's (0.014, 0.018).
- **The E6 + E11 stage 0** (B60's constants, no refit): Cerny pot D fed-NF share 0.981 - 0.991 in every member;
  Givaudan ribose / NF 0.0021 - 0.0088; Hofmann 145 C NF + Cys / ribose + Cys 11.5 (measured 5.13); Cerny pot C
  MP3P / MFT 0.075 / 0.30 / 1.09 (shipped / A100 / A133); Wang 2012 pot C xylose share 0.52 - 0.70; the guard
  violations (Kang 100 C FFT and FUR, Kang 120 C FUR, the Zhai folds, `zhou_pH7_MFT_over_FFT`, the Zhou dimer).
- **Earlier diagnostics on pot D:** EX-023 / 024 / 030 / 040 (pot D 0.98 - 0.99 in every joint arm; 0.805 with the NF
  step at 133 kJ/mol; 0.070 at 133 kJ/mol with the DDP route raised 3 decades, a diagnostic, not a fit).
- **Measured values, no engine number:** Cerny & Davidek 2003 Tables 1 - 3 (pot D MFT 7 % unlabelled, MP3P 96 %,
  3-mercapto-2-butanone 94 %, FFT 0 %); de Roos 2005 Table I (MFT 43.5 / 36 / < 0.05 / 0.9 mg per mol Cys from ribose
  / xylose / glucose / NF); Hofmann & Schieberle 1998 Tables 1 - 4, 8, 10; the Zhai 2023 Food Chem Table 1 cells and
  the Zhai 2023 JAFC Fig. 3 cells **at 140 C** (MFT 5.907 +- 0.085, FFT 11.439 +- 0.265 ug/L at 120 min, keyed in
  `thiol_submodel.yml`; FIT_HOLDOUT_DECLARATION already records the 140 / 120 MFT ratio as seen at design); Zhang 2025b
  Fig. 2 (TTCA 7.4 ... 18.2 mmol/L at 5 - 50 min; the printed ELSD curve); Cao 2017 Figs. 4 - 5 and its calibration
  table; Zhai 2021's LOD / LOQ table and its TTCA-alone zero-order fits; Mikami 2017's fed-NF decay (0.032 /min,
  pH 8, 100 C).
- **Arithmetic seen on the decisive reservoir pot (Zhang 2025b):** the ratio [TTCA] / ([Cys][Xyl]) climbs from about
  20 to about 18 600 L/mol over 5 - 50 min, against Cao 2017's apparent 60 L/mol at 90 C (an assay-reader agent's
  arithmetic, 2026-11-05). This informs no band (sec. 3.3); it is listed because it shaped the S2 design.
- Not read: any engine number on the Zhai 140 C pots, on pot D's MP3P share, on the Hofmann in-situ FFT / MFT ratio,
  or on any rival of this file.

## 1. The claim (sec. 7 item 1)

Two mechanisms, judged together because they share rows and their errors stack (CF-107; plan sec. 6):

- **Supply (M1, E10 / E12).** Free cysteine and pentose are fed from a TTCA / cysteine-Amadori reservoir whose
  exchange is slow and pH-dependent: TTCA forms on a minutes timescale, its ring opens slowly at pH 6 - 7 and faster
  in acid, so the free reactants fall on a ~10 - 20 min timescale at 120 C while TTCA holds the charge. Changed steps:
  `r_ttca_form`, `r_ttca_cys` (an exchange no longer pinned to K = 0.060 at every pH), the route steps `r_cama_*`
  (B58 / B59b) on, `r_pent_dpo` / `r_pent_tdp` retired.
- **Route shares (M2, E6 + E11).** MFT forms mostly NOT through free norfuraneol (NF) at <= 100 C (labelling: 7 %),
  through the 1,4-dideoxyosone (DDP) + H2S route instead, and the two routes carry different barriers inside measured
  classes, so the mix moves with temperature. Changed steps: `r_dpo_nf`, `r_nf_mft`, `r_ddp_mft`, `r_ddp_mft_hs`
  (rate and barrier), and on the FFT side `r_tdp_fur`, `ch_fur_decay` (barrier).

The loss law stays C1's (C2 specified, not run, M1 sec. 2.3; the W6 data need stands). The disulfide side stays the
registered C3 follow-up. Nothing here touches the pH balance (M1b), the hexose entry (M4) or the fragment route.

## 2. The evidence, by kind (sec. 7 item 2)

Ranked by the plan's fixed order. "Against" names a source that contradicts the claim. DEV / TEST from
`split.yml` (v1) and `split_v2_per_change.yml` (v2); "hold-out" from FIT_HOLDOUT_DECLARATION.

| rank | kind | source (lab) | what it shows | for / against | status |
|---|---|---|---|---|---|
| 1 | labelling | Cerny & Davidek 2003 pot D (Nestle; 95 C, pH 5, 4 h; NF + [13C5]ribose + Cys) | MFT 7 % from fed NF, 93 % from ribose; MP3P 96 % from NF | for E6 | v2 DEV for E6 / E11; hold-out, **pot D seen at design** (2026-10-30). Decisive (sec. 5) |
| 1 | labelling | Cerny & Davidek 2003 pot C (ribose + [13C5]ribose + Cys) | MP3P not detected | for E6 (NF is not formed in quantity, or not converted) | hold-out, blind companion; reported guard (sec. 5.2) |
| 1 | labelling | Cerny & Briffod 2007; Cerny 2007b (Firmenich; 145 C, [13C5]xylose, pH 4 - 7) | NF marker absent at pH 4 - 5.5, 145 C (F23#10); thiamine / xylose shares | for E6 at 145 C; the ternary share is a B60 row | calibration input; `cerny_ternary_thiamine_share` in-sample |
| 1 | labelling | Wang 2012 pot C (BTBU; 132 C) | > 99 % of MFT non-xylose in a thiamine pot | for E6 (the shipped NF route makes 0.68 of it) | v2 TEST for E6 (F23#5); BTBU independence lost in M1; side effects only |
| 1 | labelling | Liu 2023b (CAAS, 168 C); Du 2023 (CAAS, 150 C) | intact-C5 vs fragment shares | ~0 power on E6 (both MFT routes are intact C5) | v2 TEST (E6 / E11); keyed held-out tables; side effects only |
| 2 | fed intermediates | de Roos 2005 Table I (Givaudan; Laarhoven, unpublished; boiling, 1.5 h, 37 % water) | ribose gives 48x the MFT of NF per mol Cys | for E6 | v2 DEV; ledger `simulable: false`, `condition_uncertain`. Decisive with sensitivities (sec. 5) |
| 2 | fed intermediates | Hofmann & Schieberle 1998 (DFA; 145 C, 20 min, pH 5) Tables 3, 4, 10 | fed NF + Cys gives 5.13x the MFT (ppb) of ribose + Cys; NF + H2S 0.19 mol %; ribose + H2S 0.013 mol % | **against at 145 C** (NF a good precursor there): the T-signature E11 needs | fed rows FIT (B60); ribose + Cys MFT removed from B60 (validation bundle). Decisive ratio (sec. 5) |
| 2 | fed intermediates | Whitfield 1999 (Reading; fed NF, 140 C) | NF + Cys / H2S make MFT | against-ish at 140 C | B60 rows (in-sample); v2 declared held out for E6 / E11, already in R0's objective |
| 2 | fed intermediates | Zhai 2023 Food Chem / JAFC (Jiangnan; fed TTCA, TTCA + Cys, 100 / 120 / 140 C, pH 7 water) | MFT, FFT from the reservoir; 140 C pots | the route's T-law | 100 / 120 C in-sample (B60 folds, Kang = reprint, CF-004 semi-quantitative); **every Zhai >= 140 C pot hold-out**, 140 / 120 MFT ratio seen at design. Decisive (sec. 5) |
| 2 | fed intermediates | Zhai 2023 JAFC Fig. 6 (NF + sulfide, furfural + sulfide, 120 C, 20 min) | NF + sulfide -> MFT ~34 % | conflicts with Hofmann's 0.19 % (CF-065: sulfide charged as dissolved) | not in any objective; not used (CF-065 unresolved); recorded |
| 2 | fed intermediates | Mikami 2017 (fed NF alone, pH 8, 100 C) | NF decays first-order, 0.032 /min (t1/2 ~22 min), ~700x the engine's `k_nf_decay` | the rival R-NF's band (sec. 3.4) | not a fit row; a band source |
| 4 | measured intermediates | Cao 2017 Figs. 4 - 5 (BTBU; 90 C, pH 5.5 -> 4.2, 400 mM each) | TTCA + CAMA hold 90 - 92 % of the charge 0.5 - 5 h; free Cys and Xyl plateau at 15 - 20 % (they do not vanish; sec. 13 item 3); apparent K 58.5 - 64.6 L/mol (dossier arithmetic, not printed) | for a reservoir; sets K at acid pH | C1 / S2 objective (route block) |
| 4 | measured intermediates | Zhang 2025b Fig. 2 (BTBU; 120 C, pH 6.5, 40 mM each) | free Cys and Xyl fall at 0.072 / 0.060 /min (keyed points; the printed fits do not reproduce them); TTCA 7.4 -> 18 mmol/L | for a reservoir; free reactants fall much further than K = 0.060 allows | out of every objective; decisive (sec. 5) |
| 4 | measured intermediates | Zhai 2019, 2021, 2023 (TTCA, CAMA, Cys, Xyl decay, 70 - 130 C; 140 C hold-out) | TTCA alone at pH 5.5 loses 0.027 / 0.065 / 0.081 mmol/(L min) at 100 / 120 / 140 C from ~10 mmol/L | slow net TTCA loss: release slow or K large at pH 5.5 | route block (in-sample); 140 C pots hold-out |
| 7 | fitted time courses | Schieberle 2000 Table IV (DFA; 100 C, pH 5) | MFT 360 / 30 rises 1.54x, FFT shape | R0 fits FFT almost exactly; C1 lost it | out of every objective; decisive D1 (M1's) |
| -- | whole-model fits | B10 (barriers not identified), R3En, EX-024 / 030 / 035 / 036 / 041, the E6 + E11 stage 0 | the single barrier cannot move the route mix; pot D fails in every joint arm | weak (compensated states, CF-107) | recorded |

**The class evidence for the barriers (bands in sec. 3.3):** NF formation, Lau 2003 81.5 +- 3.2 and Pandit 2006 93.0 +-
6.4 kJ/mol (whey gel 116 - 131 C; potato) - appearance barriers, not a step measurement (the thiol prereg sec. 13 says
so); sulfur-volatile formation, Chan & Reineccius 1994 (81 - 137) and the B55 measured sulfur class (n 15, median
100.9, range 29.7 - 133.0, a `source_line` string not recomputed here: UNVERIFIED); Cys -> H2S, Zheng & Ho 1994 Table V
123 - 135 (measured, `k_cys_h2s` 133.0, unchanged); furfural formation, Schirle-Keller & Reineccius 1992, 147.2 (band
117 - 189, B62); furfural loss with Cys, Yaghmur 2005, 46.5 (band 46.5 - 115, B62). **No step-level barrier is
measured for NF + H2S -> MFT, DDP + H2S -> MFT, DPO -> DDP, thiamine -> MFT or TTCA -> MFT.** No barrier in this
file is computed (`assert_no_dft_sulfur` runs at import, `parameters_sulfur.py` l. 3346).

**Sources left out, with the reason:** Zhai 2023 JAFC Fig. 6 (CF-065 open); van den Ouweland & Peer 1975 (identifies
no MFT, CF-092); Grosch 1992 (E18 TEST, M1 v3.1); Yiltirak 2026 and Hofmann 1998 pH 3 / 7 (declared hold-outs:
side effects); Bornhorst 2017b's NF barriers (pH 8.4 - 9.5); Shu 1988 (area %, no standard).

## 3. The rivals (sec. 7 item 3)

### 3.1 The design: two axes, six fitted combinations

Errors stack, so the supply and the route structure are judged in combination. Supply levels: **S0** instant supply
(R0's), **S1** M1's C1 reservoir (K held at 0.060, pH- and T-flat), **S2** the reservoir with a free, pH-dependent
exchange (sec. 3.3). Route levels: **Q0** the shipped single barrier with the DDP route frozen ~4 decades weak, **Q1**
the E6 + E11 claim (sec. 3.3), **Q2** the rival explanation of the same route evidence: one NF route, with fed NF lost
to its own decay before it can make MFT (sec. 3.4).

| rival | supply | route | role | free coordinates |
|---|---|---|---|---|
| **R0** | S0 | Q0 | incumbent, refit (M1's R0, comparable to M1) | 17 |
| **C1** | S1 | Q0 | M1's C1, refit (comparable to M1) | 19 |
| **C1F** | S2 | Q0 | supply alone | 21 |
| **R0-Q1** | S0 | Q1 | routes with instant supply | 24 |
| **M2** | S2 | Q1 | **the claim** | 28 |
| **C1F-NF** | S2 | Q2 | the rival explanation of the route evidence | 23 |

S1 x Q1 and S0 x Q2 are not fitted: S1 is carried only for comparability with M1, and Q2 is a rival to Q1 at the
supply level the claim uses. Every rival runs on M1's numerics and fit-side settings: `ph.coupling = "inline"`; the
dimer's encounter-limit ceiling (M1 sec. 2.5, log10 9.47 at 145 C); the `k_mft_decay` weak prior (M1 sec. 2.6, centre
-1.391 at 100 C, sd 1 dex); `Ea_decay_thiol_sink` free in Gigl's band (7 - 102); every other switch at M1's settings
(`ph.b51` on, `sulfur.flux_budget` closed, B59a / B61 / B62 / B63 / B7a / `tdp_own_decay` / the pentose-amine entry
off); `panel.thermal_history` isothermal for every fit.

### 3.2 The base coordinates (every rival; M1 sec. 2.5)

B60's 16 identified coordinates (`k_tdp_fur`, `k_nf_mft`, `k_nf_mp3p`, `k_mgo_mp`, `k_fft_decay`, `k_fur_decay`,
`k_pent_caramel`, `k_cys_thermal`, `k_thiolate_loss`, `ph_arp_secondary_ammonium_pKa`, `k_osone_decay`, `k_ttca_deg`,
`Ea_decay_thiol_sink`, `k_dimer_mft`, `k_arp_dpo`, `k_arp_tdp`) plus `k_mft_decay` under its prior: 17, B60's bounds
except the dimer ceiling. Held at B60: the 7 unidentified and the 24 other frozen keys, except those a route or supply
level frees below. **The constants M1 held at values fitted under instant supply** (`k_dpo_ddp`, `k_ddp_mft`,
`k_dpo_nf`, `k_dpo_c2c3`, `Ea_lumped_formation`, ...; M1 sec. 2.5, review M4): Q1 frees `k_dpo_nf` and the DDP pair;
`k_dpo_ddp`, `k_dpo_c2c3` and `Ea_lumped_formation` stay held in every rival (sec. 17 question 4).

### 3.3 The changed levels, their coordinates, bands, sources and starts

**S1 (C1, as M1 sec. 2.2):** + `k_ttca_cys` (log10 k at 145 C; bounds B60 -0.579 +- 3 dex, the stage-1 convention;
Ea on `Ea_lumped_formation`, held) and `ea_form_kj_mol` (40 - 135, the `sugar_trunk` class band); K held at 0.060 L/mmol
(Cao 2017), pH- and T-flat; `k_form`, the de Roos CAMA-decay table, `ea_decay` 85.7 and `dpo_share` 0.30 declared as in
M1.

**S2 (C1F, M2, C1F-NF): the exchange made free.** S1 plus two coordinates and one new key (`sulfur.ttca_exchange`, not
built; field bands in `config._BANDS` before any fit):

| coordinate | meaning | band | source | starts 0 / 1 / 2 |
|---|---|---|---|---|
| `k_ttca_cys` | ring opening TTCA -> Cys + PENT at pH_ref, log10 k(145 C) | B60 -0.579 +- 3 dex (as S1) | stage-1 convention. **Check, not a bound:** Cao 2017's dossier floor k >= 7e-3 /min at 90 C (equilibrium within 30 min, dossier arithmetic, not printed), i.e. log10 k(145 C) >= -0.94 at 64.08 kJ/mol | C1's carried -1.23 / B8 perturbation / R0-carried route at B60 -0.579 |
| `log10 K_ref` | TTCA equilibrium at pH_ref, L/mmol | log10 0.045 - log10 0.075 | Cao 2017 blank, Figs. 4 - 5: apparent K 58.5 / 64.6 / 60.2 / 59.0 / 58.5 / 61.1 L/mol at 0.5 - 5 h (dossier arithmetic on digitised points; the paper prints no K); 90 C, pH 4.67 -> 4.15 after cooling | 0.060 / perturbed / 0.060 |
| `m_ring` | pH slope of the ring opening: k_r(pH) = k_r,ref x 10^(-m (pH - pH_ref)); formation k_f = K_ref k_r,ref, pH-flat; so K(pH) = K_ref x 10^(m (pH - pH_ref)) | 0 - 1 per pH unit | **declared mechanism, UNVERIFIED for thiazolidines:** one proton on the rate-limiting step (specific acid catalysis of the ring opening, or equivalently the cysteine amine's protonation on the formation side) gives a slope of at most 1; 0 is B58's pH-flat declaration. pH_ref 4.4 on the 25 C label scale (the middle of Cao's measured 4.15 - 4.67). A measured thiazolidine ring-chain pH law is a literature need (sec. 16.2) | 0 (= C1 exactly) / perturbed / 1 |

- **Why this form (lesson 4 of the brief, the review's "free the release rate").** In C1, K = 0.060 holds free cysteine
  near 19 mmol/L in Zhang 2025b's 40 + 40 mM pot at any `k_ttca_cys` (arithmetic: 0.06 c^2 + c - 40 = 0), so freeing the
  release rate alone cannot let free cysteine fall at the measured rate. The two measured anchors of the exchange sit at
  different pH: Cao 2017 at pH ~4.4 (apparent K ~0.060 for 5 h) and Zhai 2021's TTCA-alone pots at pH 5.5 (a slow net
  loss, 0.0027 / 0.0065 / 0.0087 /min first-order equivalents at 100 / 120 / 140 C, arithmetic on its printed zero-order
  fits). An exchange that slows with pH reconciles them without either lab's number being wrong.
- **The band comes from Cao 2017 and the mechanism, not from the scored pot** (`tasks/lessons.md`, 2026-10-31: a
  rival's bands come from data other than its scored tests). The brief asked for a band "wide enough to include" Zhang
  2025b's timescale. The check, recorded as seen (sec. 0) and used to set nothing: at m = 1, pH 6.5, K reaches ~7.5
  L/mmol and the equilibrium free cysteine in Zhang's pot ~2 mmol/L (arithmetic), so the band does not exclude the
  measured fall. If a reviewer finds the slope ceiling unsupported, widening it on Zhang's numbers would make the
  reservoir block partly in-sample: that is an owner-level decision (sec. 16.1 item 1), not an author edit.
- The route block (in every S1 / S2 objective) spans pH 4.2 - 8 (Cao, Zhai 2019 / 2021 / 2023, Zhao 2017b), so
  `m_ring` is identified in-sample if at all; Zhang 2025b (pH 6.5) stays out.
- **Limits the S2 form leans on, to be checked per pot (`tasks/lessons.md`, 2026-11-01, limit arguments):** the build
  prints, for every route-block and decisive pot at C1F's start and at each carried optimum, k_r x window and
  k_f x (Cys + Xyl) x window. Neither "fast equilibrium" nor "irreversible formation" is assumed anywhere.

**Q1 (R0-Q1, M2): the E6 + E11 claim.** + 7 coordinates through the B62 hook (`FurfuralBarriers.barriers_kj_mol` takes
any fitted key, read by `with_fitted_sulfur` via `measured_ea_overrides()`; a ship adds a named switch). Each re-barriered
step keeps its fitted log10 k at 145 C, so a barrier moves only the off-145 C rates.

| coordinate | meaning | band | source | starts 0 / 1 / 2 |
|---|---|---|---|---|
| `k_dpo_nf` | DPO -> NF, log10 k(145 C) (B60 frozen +0.434) | -10 to 0.5 (`FITTED_SULFUR_BOUNDS_LOG10K`, numerical) | none measured | B60 / perturbed / B60 |
| `d_ddp` | one log10 offset on `k_ddp_mft` and `k_ddp_mft_hs` together (B60 frozen -6.541 / -7.693) | -2.31 to +7.04 (each key inside -10 to 0.5; arithmetic) | none measured; the identifiable product (round-3 review); EX-040's diagnostic pass sat at +3 | 0 / perturbed / +3 |
| `Ea_dpo_nf` | barrier of DPO -> NF | 75.1 - 105.8 kJ/mol | Lau 2003 81.5 +- 3.2 and Pandit 2006 93.0 +- 6.4, their +-2 sd envelope (`lau2003_extraction.md`, `pandit2006_extraction.md`); appearance barriers, declared transfer | 83.8 (inverse-variance mean) / perturbed / 93.0 |
| `Ea_nf_mft` | barrier of NF + H2S -> MFT | 40 - 133 | sulfur-volatile class: Chan & Reineccius 1994; B55 class (range 29.7 - 133.0); the thiol prereg's sec. 13 ceiling 133 | 64.08 / perturbed / 100.9 (class median) |
| `Ea_ddp_mft` | barrier of DDP + H2S -> MFT (and `_hs`, tied) | 40 - 133 | as `Ea_nf_mft` | 64.08 / perturbed / 100.9 |
| `Ea_tdp_fur` | barrier of TDP -> furfural | 117 - 189 | Schirle-Keller & Reineccius 1992 Table 2, 147.2 (refit 156); B62's band | 147.2 / perturbed / 156 |
| `Ea_fur_decay` | barrier of furfural loss (leaves the carbonyl-sink family's shared 143.54) | 46.5 - 115 | Yaghmur 2005 46.5 (with Cys); B62's band | 46.5 / perturbed / 80 |

- Held: `k_dpo_ddp` (+0.042) and its barrier (64.08): no measurement exists, and freeing it beside `d_ddp` adds a
  non-identified direction (the earlier proposal's rule e). The fragment route (`r_ha_mp_mft`), thiamine (`k_thi_hmp`,
  `k_hmp_mft`) and `k_cys_h2s` (133.0, measured) unchanged.
- **Why not the B10 route map:** `FORMATION_ROUTE_OF` puts `k_dpo_nf` and `k_dpo_ddp` both in `sugar_trunk` and
  `k_nf_mft` and `k_ddp_mft` both in `thiol_assembly`, so two route barriers give both MFT routes one composite barrier
  and the NF share stays T-invariant (the earlier proposal, sec. 0). Per-step barriers are the only form that can move
  the mix.
- The composite gap NF route minus DDP route spans about -82 to +135 kJ/mol across the bands (arithmetic: (75.1 + 40)
  - (64.08 + 133) to (105.8 + 133) - (64.08 + 40)), so either route can win at either end of 95 - 145 C: the band does
  not encode E6's answer.

**Q2 (C1F-NF): one NF route, fed NF lost to its own decay.** + 2 coordinates; the DDP pair and every Q0 barrier as R0:

| coordinate | meaning | band | source | starts 0 / 1 / 2 |
|---|---|---|---|---|
| `k_nf_decay` | NF -> fragments, log10 k at **100 C** (reparameterised so the measured bound applies directly) | -10 to log10 0.032 /min | Mikami 2017 Fig. 6B: 0.032 /min (+- ~10 %), 1.7 mM NF alone, 0.2 M phosphate **pH 8**, 100 C. Base-catalysed (E16, M1b), so an upper bound for the pH 5 - 7 pots | B60's value / perturbed / log10 0.032 - 0.5 |
| `Ea_nf_decay` | its own barrier (leaves the carbonyl-sink family) | 40 - 135 | `FORMATION_ROUTE_OF` puts `k_nf_decay` in `sugar_trunk`; that class band (`FORMATION_EA_BOUNDS_BY_ROUTE`); no NF-decay barrier is measured at pH <= 8 | 135 (nearest to B60's 143.54) / perturbed / 80 |

- **What Q2 predicts, from its mechanism:** fed NF is spent early, so pot D's MFT comes from ribose-made NF later. But
  then MP3P, which in the engine forms only from NF (`k_nf_mp3p`), must also come mostly from ribose-made NF: pot D's
  MP3P share falls with the MFT share. Measured, MP3P is 96 % from fed NF while MFT is 7 %. Q1 keeps the NF pool fed-
  dominated (MP3P high) and moves MFT to DDP (MFT share low). **Pot D's MP3P share is therefore the unit that separates
  Q1 from Q2** (sec. 5.1). A second separation: Hofmann's 145 C fed NF + Cys pot needs NF to survive 20 min at 145 C.

### 3.4 Start vectors and order

- Wave 1: R0, C1, C1F. Start 0 = M1's carried optimum for R0 and C1 (`kinetic_core_m1_supply_combinations_carried_*`;
  the M2 objective's cost at it is the V-reference); C1F's start 0 = C1's carried optimum with `K_ref` 0.060 and
  `m_ring` 0, i.e. C1 exactly. Start 1 = B8's perturbation of start 0 (seed printed). Start 2 as in the tables (R0:
  B60's start-1 endpoint, as M1; C1: R0's carried optimum with the route coordinates at their starts).
- Wave 2, after wave 1's carried optima: R0-Q1 from R0's carried optimum, M2 and C1F-NF from C1F's, the new coordinates
  at the tables' starts. A Q1 start 0 is NOT the parent rival exactly (`Ea_dpo_nf` and `Ea_tdp_fur` start inside their
  bands, not at 64.08); V3 checks the reproduction with the bands disabled.

## 4. Objective, fit and numerics

### 4.1 Objective

As M1 sec. 4.1, unchanged, so M1 and M2 are comparable:

| block | rows | rivals | measurement model |
|---|---|---|---|
| B60's 52 rows | B60's `row_ids` | all | B60's sigmas, `weight_tag` half, pH endpoint weight 0.7, sigma_pH 0.5 |
| route block | the 296 verified rows in 34 pots of `route_submodel.yml`, minus the Zhang 2025b pots and the five Zhai 140 C pots (M1 sec. 4.1) | S1, S2 rivals | read error / printed SD + s_series (M1's estimate - re-estimate protocol) |
| prior | the `k_mft_decay` term | all | M1 sec. 2.6 |

Not in any objective: every decisive unit (sec. 5), every side-effect row (sec. 10), Grosch 1992, Zhai 2023 JAFC Fig. 6.
**What identifies the route coordinates in-sample** (stated so the reader knows the predictions are out of sample, not
unconstrained): the fed NF, ribose and TDP pots at 145 C (Hofmann), Whitfield's fed NF at 140 C, the in-situ NF level
(`hofmann_ribose_NF_insitu`), the Kang / Zhai 100 / 120 C rows, Zhou 2023 at 120 C and the Cerny ternary share. No
labelling row is in the objective. If a Q1 coordinate is not identified by these (Laplace rule), it is held at the
optimum for the predictive (sec. 7) and its profile is printed; the decisive comparison then reads a sharp prediction
at a held value, and the result says so.

### 4.2 Optimiser, starts, convergence, gradient accuracy

- `fitting.fit`, `LEAST_SQUARES_SETTINGS` (trf, x_scale jac, ftol 1e-6, xtol 1e-8), the module's absolute FD steps;
  budget 1 200 evaluations per start; status > 0 counts; one continuation (+1 200) from a status-0 endpoint.
- **Gradient accuracy (the bake-off's finding: the forward-difference Jacobian up to 7 % off its own central difference
  on some columns).** Pre-declared: every start runs forward differences to convergence, then a **central-difference
  polish**: continued from its endpoint with the central-difference Jacobian at the same absolute steps until status
  > 0 (<= 300 evaluations; status 0: once more; still 0, the forward optimum is carried, flagged). The polished point is
  the start's optimum. Laplace, profiles' Hessians and every Jacobian in sec. 7 are central. **Stability record:** at
  each carried optimum, the cost change and each decisive unit's prediction change from the forward endpoint to the
  polished one are printed beside the unit's tie band (sec. 7.3); a change larger than 0.1 band flags the unit
  "FD-sensitive" for the results review. (UNVERIFIED that `fitting.fit` exposes a central-difference option at the
  module's absolute steps; if not, the build adds it as a numerics option with a unit test that it reproduces
  `scipy.optimize.approx_fprime`'s central form on a toy objective.)
- **Comparison rule (M1 sec. 4.3, carried):** a rival is compared iff two starts end status > 0 within 1.92 nats and
  agree on every decisive prediction within max(the unit's measurement sd, the Laplace predictive sd), read without y;
  otherwise the verdict is computed per qualifying start and must agree, else "not compared, flagged". Fallbacks as
  M1 (R0's lowest converged start; R0-shipped (b)).
- **Series noise** for S1 / S2 rivals: M1's estimate at start, re-estimate at the optimum, one continuation (<= 300 + 300).
- **Run conditions (carried from M1 and H0's amendment A2):** every pot under its wall-clock alarm; free memory checked
  at launch; two fits in parallel only if both fit in RAM; never SIGSTOP a fit. **A record with `failed_pot_runs` > 0
  is re-run with the same settings, and amendment A1's rule applies:** it is accepted only if every failure was inside
  the optimiser's segments, the cost recomputed from the optimum's residuals matches (relative 1e-9), and a
  same-settings re-run reproduced it exactly (cost to 1e-12, evaluations, failure count and pots). `failure_accepted`
  in the generator, as M1.
- **G1 (M1.0)** at each carried optimum: oxidant use and the phi = 1 rise, beside B60's.

### 4.3 Verification before any fit

- **V1** (M1's V1 and V6): every new key off: the M2 generator reproduces M1's carried costs (R0 21.88294, C1 186.75261)
  at their vectors, inline, fit rtol, to relative 1e-9.
- **V2** S2 at `K_ref` 0.060, `m_ring` 0 reproduces C1's carried cost and every residual exactly.
- **V3** Q1 with its bands disabled, barriers at 64.08, `d_ddp` 0, `k_dpo_nf` at B60: reproduces the parent; Q2 with
  `Ea_nf_decay` 143.54 and `k_nf_decay` at B60's: reproduces C1F.
- **V4** the pot D tracer extended to MP3P (the same fed-fraction weighting on `r_nf_mp3p`), checked against an explicit
  two-isotopologue integration (NF duplicated as NF_fed and NF_rib, each with every NF reaction) on pot D at R0's and
  C1F's start vectors: both shares agree to 1e-4 absolute. The tracer's premise (every NF sink first-order in NF) is
  checked in code against `sulfur.FULL_REACTIONS`.
- **V5** inline vs `fixed_point_nodes` (81 nodes, rtol 1e-8) on every decisive pot at R0's and C1F's starts, as M1's V7:
  each decisive prediction within 0.1 of the unit's measurement sd.
- **V6** smoothness (M1's V8) on `m_ring`, `K_ref`, `d_ddp`, `Ea_nf_mft`, `Ea_ddp_mft`, `k_nf_decay` at C1F's start: the
  inline FD at +-0.025 / 0.05 dex (+-0.25 / 0.5 kJ/mol) keeps its sign and stays within a factor 2 of the fixed point's;
  if it fails, every fit runs on the fixed point (M1's fallback; sec. 14's basis changes, the 8.8x of W1.3 sec. 4).
- **V7** the predictive-density machinery (sec. 7.2) on a case with a known answer: R0's RES_Cys under the M1
  alternative floor by 1-D quadrature on `k_cys_thermal` must reproduce the M1 review's scan (about -746); and the IS
  estimator must match quadrature on that case within its own MC SE.
- **V8** the scorer: the decisive units at B60 (R0-shipped (b)) reproduce stage 0's values for pot D (0.981), Givaudan
  (0.00261 at 2 M) and Hofmann (11.5) to 1e-3 relative, on the same paths.
- **V9** every decisive pot's limit numbers (sec. 3.3) and the maximum MFT on the panel (the dimer cap, M1 V5).

### 4.4 Feasibility probe (M1 sec. 3, carried, per new level)

For C1F, M2 and C1F-NF: 256 Sobol draws over the new coordinates' bands, everything else at the parent's start; each
draw scored on the rival's own claim units (sec. 5.1) by the plug-in density; a local search (<= 200 evaluations) from
the best draw. **One rule:** a rival is killed only if its best value is not above R0-shipped (b)'s on those units;
**stacking guard:** a would-be kill is confirmed by one fit (start 0) before it stands. Passing is screening, never
evidence. Artifact `kinetic_core_m1m2_feasibility_probe.{md,json}`.

## 5. The decisive data (sec. 7 item 4)

Every unit is out of every objective. A **unit** is the independent piece the measurement model scores.

| block | unit | data | rivals predict differently because | lab; v1; v2 | seen |
|---|---|---|---|---|---|
| **S** reservoir | S-Cys, S-Xyl | Zhang 2025b Fig. 2g - h blank, 8 points each, 5 - 50 min (M1's) | S0: gone in seconds (floored); S1: held ~19 mM by K; S2: falls with the exchange | BTBU; DEV; E12 TEST lab (F28#25), not these rows | at start (M1) |
| S | S-TTCA (reservoir rivals only) | Zhang 2025b Fig. 2a - b blank TTCA, 8 points | S1's K caps TTCA; S2 lets it rise to the plateau | as above | read (sec. 0) |
| **D1** thiol shape | D1-MFT, D1-FFT | Schieberle 2000 Table IV, 360/30 and 720/360 ratios, bivariate per species (M1's) | supply timing and the FFT branch | DFA; DEV; E12 DEV | at start (M1) |
| **Rt** route shares | Rt-D | Cerny 2003 pot D, MFT fed-NF share 0.07 | Q0: NF route ~0.98; Q1: DDP can carry MFT; Q2: fed NF spent | Nestle; F23#2 DEV; DEV | **seen at design** |
| Rt | Rt-DMP | Cerny 2003 pot D, MP3P fed-NF share 0.96 | Q2 lowers it with Rt-D; Q0 / Q1 keep it high | as above | read (sec. 0) |
| Rt | Rt-G | de Roos 2005 Table I, MFT ribose / NF per mol Cys = 48 | as Rt-D, at ~2 M and 100 C | Givaudan; DEV; DEV | stage 0 (0.0026) |
| Rt | Rt-H | Hofmann 1998, 145 C, pH 5: MFT(ribose + Cys) / MFT(NF + Cys), ppb, panel path; measured 1 / 5.13 | the T-signature: NF a good precursor at 145 C; Q2 must keep NF alive 20 min at 145 C | DFA; DEV; DEV | stage 0 (11.5) |
| **T** temperature (E11) | T-MFT, T-FFT | Zhai 2023 Food Chem fed TTCA, pH 7, MFT and FFT at 120 min, 140 C / 120 C | per-route barriers vs the single 64.08; supply's T-law | Jiangnan (`ho_network`); v1 UNVERIFIED; DEV | 140 / 120 MFT seen at design; cells read |
| T | T-MFTc | Zhai 2023 JAFC fed TTCA + Cys, MFT at 120 min, 140 C / 120 C | as T-MFT, with free cysteine present | as above | cells read |

- **Hofmann's ratio has an in-sample numerator** (`fed_nf_cys_MFT` is a B60 row; ribose + Cys MFT is a removed row,
  out of every objective). Every rival carries the same row, so the comparison is fair; the unit's information is the
  out-of-sample ribose arm. Its 145 C pot has the Hofmann autoclave come-up (sec. 6.3).
- **Givaudan** is scored on the earlier proposal's reconstruction (Cys : carbonyl : lactate 1 : 1 : 1.67, 37 % water
  read as ~2 M, pH 5 clamped, 100 C, 90 min) although the ledger marks it `simulable: false` and `condition_uncertain`;
  the 0.3 M reading is a required-agreement sensitivity (sec. 7.3) and the owner decides whether it stays decisive
  (sec. 16.1 item 2).
- **Pot D's concentration** is nominal (500 mg total in a 2 mL vial, read as 0.5 mL: NF 145 + [13C5]ribose 145 + Cys
  96.6 mmol/L); x0.5 and x2 are required-agreement sensitivities (bimolecular steps make the share
  concentration-dependent).
- **Moved out of the decisive set since the earlier proposal:** its D4 (F23#13, one lab's unlabelled in-house estimate,
  contradicted by Cerny 2007 at the same pH) and its held-out H1 / H2 (side effects, sec. 10).

### 5.1 Contrasts and their claim units (pre-declared from the mechanisms, sec. 3)

| contrast | pairs (a vs b) | claim units (count toward "a favoured") | guard units (none may favour b) |
|---|---|---|---|
| **supply** S2 vs S0 | C1F vs R0; M2 vs R0-Q1 | S-Cys, S-Xyl | D1-MFT, D1-FFT, every Rt and T unit |
| supply, reported | C1 vs R0 (M1 comparability); C1F vs C1 | S-Cys, S-Xyl (+ S-TTCA for C1F vs C1) | as above |
| **route** Q1 vs Q0 | M2 vs C1F; R0-Q1 vs R0 | Rt-D, Rt-G, Rt-H | Rt-DMP, D1, S, T |
| **temperature law** (E11) | M2 vs C1F; R0-Q1 vs R0 | T-MFT, T-FFT, T-MFTc | D1, S |
| **route vs rival explanation** Q1 vs Q2 | M2 vs C1F-NF | Rt-DMP, Rt-H, T-MFT, T-MFTc | Rt-D, Rt-G, D1, S |

Where the rivals predict the same thing on every claim unit (e.g. S0 and S2 on pot D's MP3P), the unit is a guard, and
the record says the contrast cannot be decided there.

### 5.2 Reported, not decisive

Cerny pot C MP3P / MFT (blind companion; "not detected", no detection limit printed): reported against the thiol
prereg's declared 0.1, never a unit. Hofmann's in-situ FFT / MFT at 145 C (both removed rows): reported. The reservoir
record of M1 (R0's structural zeros on TTCA / CAMA; Cao Fig. 4).

## 6. Tolerances: the measurement model per series (sec. 7 item 5)

### 6.1 Error terms per series

**omega** is H0's within-laboratory engine error (omega(all) 0.390 dex, reset set; EX-059). It is a model-discrepancy
term, not the measurement's noise, so (sec. 7 item 5: "the measurement's own reported noise") it enters every unit that
compares different pots or runs **as a required-agreement sensitivity, not the primary**; it never enters a ratio inside
one injection. **tau** (H0's tau_ABS(all) 0.964 dex) applies only to absolute levels compared against constants from
other laboratories; **no decisive unit is such a level**, so no unit carries tau (lesson 2).

| unit | value basis | measurement noise (primary) | omega | tau | floor (primary; sensitivities) |
|---|---|---|---|---|---|
| S-Cys, S-Xyl | ln(C / C0), C0 "the concentration of original reaction solution" (caption): read as measured by the same assay, so a multiplicative calibration error cancels (lesson 2) | per point: Cao 2017's same-method statement, "RSD were all less than 10 %" (ms p. 10; the method is Cao's, Zhang p. 3), sd log10(1.10) = 0.041 dex, in quadrature with the read error (0.026 / 0.035 dex). Zhang's Fig. 2 caption states n = 3, mean +- SD, but no bars are visible on panels g - h at 400 dpi | sensitivity, per point | none | **0.025 / 0.020 mmol/L**: Zhai 2021's LOD for Cys / Xyl by HPLC-ELSD (3.025 / 3.02 mg/L, p. 10650; a same-method-class transfer: Zhang prints no LOD). Sensitivities: Zhang's printed curve's zero-signal point 0.223 mmol/L (58.733 / 263.28, arithmetic on p. 3's curve); Zhai 2021's LOQ 0.090 / 0.070. Every one is below the blank's lowest point (0.80 / 1.22) |
| S-TTCA | absolute mmol/L, LC-MS/MS external standard (Table 1; standards 0.1 - 20 mmol/L) | Cao 2017's RSD < 10 % (same lab, same analyte class) with the read error (+-0.3 mmol/L) | sensitivity | **none, declared:** the reservoir rivals' TTCA comes from K and k measured in the same laboratory's assay (Cao 2017, BTBU), so the comparison is within-lab; a between-lab offset is a sensitivity at tau_ABS(all) | Zhai 2021 TTCA LOD 0.0087 mmol/L (2.21 mg/L); sensitivity: the lowest standard 0.1 |
| D1-MFT, D1-FFT | log ratios within one series (scale cancels) | M1's: Hofmann 1998 Table 1 fn b ("triplicates ... not more than 10 %"), 0.041 dex per point, 0.059 - 0.061 per ratio, bivariate (covariance from the shared 360 min point) | sensitivity, per point | none | **1 ug/L**: Table IV's own "< 0.1 ug" per 100 mL pot for co-analysed compounds (no LOD printed; SIDA inferred, not stated). Sensitivity: no floor (the integrator's value) |
| Rt-D, Rt-DMP | log10 odds of the share, one injection (HS-SPME GC-MS, natural-abundance corrected, 98 % enrichment) | Cerny 2003 prints none ("duplicate", p. 2715). Declared transfer: Cerny & Briffod 2007's "standard deviation did not exceed 23 %" (p. 1553, TIC areas, same first author) on each isotopologue area, independent: sd(log10 odds) = sqrt(2) log10(1.23) = **0.127 dex**; conservative (co-eluting isotopologues of one injection share most of their error) | none (one injection) | none | none needed (a share); a rival whose total MFT or MP3P in pot D falls below 1e-6 of the charge is flagged |
| Rt-G | log10 ratio of two pots of one study | none printed; declared: the same 23 % per value, 0.127 dex per ratio | sensitivity | none | 0.05 mg MFT per mol Cys, the table's own "< 0.05" (glucose arm) |
| Rt-H | log10 ratio of two pots of one study, SIDA, ppb | Hofmann 1998's 10 % per value: 0.059 dex | sensitivity | none | 0.1 ug per pot (2 ug/L): the same laboratory's reporting floor (Schieberle 2000 "< 0.1 ug"); Hofmann 1998 prints no "<" entry |
| T-MFT, T-FFT | log10 ratio 140 / 120 C within one study (calibration curve shared, Zhai 2023 JAFC Table S2) | the printed SDs (Food Chem Table 1, n = 3): MFT 5.907 +- 0.085 / 1.388 +- 0.256, FFT 11.439 +- 0.265 / 4.107 +- 0.137: **0.080 / 0.018 dex** per ratio (arithmetic) | sensitivity | none | MFT 0.28 ug/L (the printed curve's blank-equivalent, 0.0011 / 0.0039, arithmetic); FFT 0.01 ug/L; sensitivity: no floor |
| T-MFTc | as T-MFT (JAFC Fig. 3, read errors only) | declared: the TTCA arm's relative SD at the same 120 min, 120 C cell (0.256 / 1.388 = 18 %) per value: 0.104 dex per ratio | sensitivity | none | as T-MFT |

**Why S-Cys / S-Xyl carry no level term, and the check.** If C0 were the nominal charge instead of a measured value, a
calibration factor would not cancel. Required-agreement sensitivity: each series re-expressed relative to its own first
printed point (5 min), ln(C(t) / C(5)), which cancels any multiplicative calibration factor whatever C0 is (the 5 min
point's noise enters as a shared covariance). R0's failure survives that form (floored at every point: a flat line
against a 3.2 ln-unit fall).

**A floored prediction** is scored at the floor (the assay would read the floor). Floors are never data points of the
scored series. The count of floored points per unit, rival and form is printed.

### 6.2 Noise sensitivities (required agreement)

Each unit's measurement sd at 0.5x and 2x; omega added where the table says; the alternative floors; S's first-point
form; Rt-D's concentration x0.5 / x2; Rt-G at 0.3 M. Never used to choose a value; used only in the sign agreement of
sec. 7.3.

### 6.3 Thermal readings (M1 sec. 6, carried)

Every fit isothermal. Decisive units scored, no refit, under every registered reading (CR-03, CR-04, B64's tau 4 and
6.7): the Hofmann autoclave pots (Schieberle's 100 C pot for D1; Hofmann's 145 C, 20 min pots for Rt-H: a 20 min hold
after a 20 -> 145 C come-up of the same order) take all six readings. Zhang 2025b, Cerny pot D (heating block, 2 mL
vials), Givaudan (boiling) and the Zhai pots (no come-up stated) are isothermal under every reading (B64's convention, "no
invented time constants"); the Hofmann-rate stand-in ramp is printed for them outside the agreement requirement (M1
v3.3). A reading applies to every pot of one evaluation pairwise by name.

## 7. Predictive densities and the tie band

### 7.1 The predictive (lesson 3: no linearised density decides)

For unit u, rival k: lppd_uk = log of the integral of N(y_u; f_uk(theta), Sigma_u) over the rival's parameter
distribution: the Laplace Gaussian at the carried optimum on the identified coordinates (central Hessian); bound-active
coordinates by their one-sided profile density (M1 sec. 7); unidentified coordinates held at the optimum. f is computed
by the engine at every draw (exact in f; only the parameter posterior is approximated). Floors as sec. 6.1.

**Primary estimator: importance sampling with an unit-tuned proposal.** The posterior predictive misses y by many sd for
some rivals (M1: R0's Cys), where draws from the posterior never reach y and plain MC has an effective sample size near 1
(M1's MC SE 1.00). So the proposal is a multivariate Student-t (5 degrees of freedom, heavier-tailed than the target)
centred at theta*_u = argmax [log N(y_u; f_uk(theta), Sigma_u) + log posterior(theta)] (found on the unit's pots only,
<= 300 pot-evaluations), with the inverse central Hessian there as its scale. Weights posterior x likelihood / proposal;
lppd = log mean weight. The proposal's draws are reused across the noise sensitivities (the likelihood is re-weighted
analytically); floors and thermal readings need new pot runs at the same draws.

**Also printed, every unit:** the plug-in density (theta at the optimum; required agreement, sec. 7.3, so a sign bought
by parameter vagueness cannot decide); the linearised density (diagnostic only; flagged where it differs from the IS
value by more than the tie band); plain posterior MC (64 draws, diagnostic).

**Where IS fails** (ESS below its requirement after S = 4 096 draws): the unit is computed by Gauss-Hermite quadrature
on the one or two identified coordinates that carry > 90 % of the unit's linearised predictive variance (the
review's method), the rest held; if two do not carry 90 %, the unit is "density not resolved" and counts as a tie.

### 7.2 Effective sample size, derived from the tie band

The tie band b_u (sec. 7.3) is estimated from R = 4 000 replicate redraws, so its own relative precision is about
1 / sqrt(2R) = 1.1 %. MC error adds in quadrature to the replicate spread; it is negligible when it widens the band by
less than the band's own precision: sqrt(1 + (SE_MC / sd_rep)^2) - 1 <= 1 / sqrt(2R), i.e. **SE_MC(d_u) <= 0.15 sd_rep(d_u)**
at R = 4 000 (arithmetic). With SE_MC(lppd) ~ sqrt(1 / ESS - 1 / S), each rival's ESS must reach 2 / (0.15 sd_rep(d_u))^2
(the factor 2: two rivals' errors in d_u). S starts at 256 and doubles until both rivals meet it (cap 4 096, then sec.
7.1's quadrature). ESS and SE_MC are printed per unit, rival and form; no fixed draw count is used.

### 7.3 The tie band and when a unit favours a rival

- **d_u = lppd_u(a) - lppd_u(b).** Its replicate spread sd_rep(d_u): redraw y_u ~ N(y_u, Sigma_u) R times and recompute
  d_u at the same predictive distributions (the data's own noise moves the comparison by this much).
- **b_u = z_0.95 x sd_rep(d_u)**, the half-width of d_u's central 90 % replicate interval. The 90 % level is the plan's
  interval convention (H0 ships 90 % intervals), not a choice made here.
- **A unit favours a** iff d_u > b_u under the primary (IS) density AND d_u has the same sign under every
  required-agreement form: the plug-in density; each noise sensitivity, floor alternative and concentration reading of
  sec. 6.2; each thermal reading (sec. 6.3); the forward-difference optimum (sec. 4.2). Favours b symmetrically.
  Otherwise a tie. (In M1 a sign flip on one non-resolved form made the deciding unit a tie; here every form that can
  flip a sign is resolved by sec. 7.1 - 7.2 before it is read.)

## 8. The decision rule, the outcomes, what ships

### 8.1 Blocks and contrasts

- **A block favours a over b** iff more than half of the contrast's claim units in that block favour a and no unit of
  the block favours b; **conflict** if units favour both; otherwise **not separated**. With two claim units this is M1's
  rule (both must favour); with three, two of three.
- **A contrast credits a** iff, in EVERY pair listed for it (sec. 5.1), its claim units' blocks favour a, no guard unit
  favours b, and the in-sample guard (8.2) holds. If the pairs disagree (e.g. S2 wins at Q0 but not at Q1), the
  contrast is **"conflict: the errors stack"** and nothing on that axis ships.
- **Evidence rank breaks no tie and overrides no win; it limits what a loss can be excused by** (8.2).

### 8.2 The in-sample guard: B60's rows and the route block (lesson 5)

The incumbent's own rows are in every objective, so a rival cannot "win" on them (more coordinates fit more), but it can
lose them, as C1 lost 36 nats on B60's 52 rows. Per row group (paper x species, e.g. Kang FFT, Hofmann fed FFT,
Kumazawa, Zhou pH 7): Delta c_g = cost_g(a) - cost_g(b) at the carried optima, and its replicate spread from the rows'
declared sigmas (sd = sqrt(sum over the group's rows of (r_a - r_b)^2), arithmetic for Gaussian rows); a group is
**lost** if Delta c_g > z_0.95 x that sd. A credited rival with a lost group:
- ships only if every lost group is named in the compensating-error table (sec. 10) with its candidate error, AND no
  lost group's evidence kind ranks above the best-ranked claim unit the rival won on (losing fed-intermediate rows,
  rank 2, cannot be excused by winning measured intermediates, rank 4);
- otherwise the contrast is **"credited, held: an unexplained in-sample loss"** (sec. 7 item 7: look further before
  shipping), and nothing ships on that axis until the results review answers it.
- The FFT case the review names: D1-FFT (a decisive guard unit) and the FFT row groups. A supply rival that wins S but
  has D1-FFT favouring R0 fails the supply contrast outright (a guard unit); a loss confined to in-sample FFT groups
  follows the two bullets above.

### 8.3 Outcomes and what each would mean

| supply contrast | route contrast | Q1 vs Q2 | meaning | what ships |
|---|---|---|---|---|
| S2 credited | Q1 credited | Q1 credited | the reservoir with a pH-dependent exchange and the DDP route with per-route barriers explain the reservoir, the labelling and the T-shape together | **M2's structure and constants** (sec. 8.4) |
| S2 credited | Q1 credited | not separated / Q2 | supply settled; the route evidence equally explained by fed-NF loss | the S2 structure at C1F's constants; E6 stays open, Q1 not shipped; E16's NF-decay pH law (M1b) becomes the deciding experiment |
| S2 credited | not separated / conflict | -- | supply settled, routes not decided by the data on disk | S2 at C1F's constants; E6 / E11 to the W6 request (Cerny pot D at 120 / 145 C; Cerny's raw 145 C isotopologues) |
| not separated / conflict | Q1 credited | Q1 credited | routes settled, supply not separable | Q1 at R0-Q1's constants (instant supply kept); CF-107 open |
| S0 credited | -- | -- | the measured reservoir is better explained by instant supply: contradicts EX-060 and M1's data; investigate before anything | nothing; a results-review item |
| not separated | not separated | -- | the data on disk cannot decide either step under this measurement model | nothing; CF-107 and E6 / E11 open, recorded as not evidence against S2 or Q1 |
| any, with T favouring Q0 | -- | -- | the route shares may be right at 95 - 100 C for a reason the T-law contradicts | Q1 not shipped (T is a guard of the route contrast); E11 open |

E11 is credited separately only if the temperature-law contrast credits Q1; otherwise a shipped Q1 is recorded as "route
shares at 95 - 145 C, T-law not attributed".

**Comparability with M1** is reported, not decided: C1 vs R0 on S and D1 under M1's rule and under this measurement
model side by side (what the review's three fixes change), and C1F vs C1.

### 8.4 What "ships" means here

- The winning cell's switches (`sulfur.cys_amadori_route`, the new `sulfur.ttca_exchange`, a named route-barrier
  switch replacing the in-memory B62 hook) and its carried constants replace B60 for the sulfur lane, after the
  side-effect record (sec. 10) is complete, every review objection answered, and the coupling flip shipped first or
  declared part of the ship (M1 sec. 2.5).
- **Owner sign-off before merge** (sec. 16.1 item 3): every pentose MFT and FFT prediction moves.
- The constants stay provisional until the plan's joint refit (JR); the loss law (C1 vs C2) stays undecided.

## 9. Compensating structure: bound-active and switched-off coordinates (lesson 6)

- **Bound-active** (`active_mask` at the carried optimum): profiled at three points inside the bound (0.25 / 0.5 / 1.0
  dex; 1 / 3 / 10 kJ/mol), others re-optimised (<= 300 + 300 evaluations per point, central polish). Each decisive unit's
  d_u is recomputed at each profile point against the other rival's carried optimum. If any unit's side changes, every
  verdict that uses the unit is labelled **"conditional on <coordinate> at its <measured | numerical> bound"**, and a
  verdict conditional on a NUMERICAL bound does not ship (a measured class bound does: the class is the physics).
- **Switched off** (a rate at its lower bound -10, "the step does not happen", as C1's `k_fur_decay` and
  `k_osone_decay`): one extra profile point with the coordinate fixed at the value it takes in the other rival's
  carried optimum (or B60's, if that one is also off), others re-optimised. A side change labels the verdict
  "conditional on <step> switched off", and that verdict does not ship.
- **Shared-coordinate divergence record:** every coordinate shared by a compared pair whose optima differ by more than
  the root-sum-square of their Laplace sds, with the blocks whose cost moves most along it.
- The Gigl-ceiling rule of M1 (sec. 2.6) is this rule applied to `Ea_decay_thiol_sink` (a measured bound): D1 is
  labelled "conditional on the Gigl ceiling" if the sides change, and the no-refit 112 / 122 kJ/mol sensitivity is
  printed as in M1.

## 10. Side effects (sec. 7 item 7)

Before (R0-shipped (b), B60 on inline) and after (every compared rival), every row, as M1 sec. 8: the numerics row
(inline vs fixed point); the full held-out contrast benchmark (`score_contrast_benchmark.py --score` on a named Config
and vector; coverage first); the in-sample rows per group and the route residuals by pH; G1; the E6 / E11 TEST cells the
scorer can answer (CAAS: F23#6 / F24#6 Liu 2023b, F23#8 / F24#5 / F26#24 Du 2023; BTBU's F23#5 / F23#7 / F24#4 / F25#23
tagged "laboratory independence lost in M1"); the earlier proposal's held-out tests H1 (Liu 2023b fragment share) and H2
(Wang 2012 pot C xylose share, with the thiamine flux beside it); Cerny pot C MP3P / MFT; the declared hold-outs
(Yiltirak 2026, Hofmann 1998 pH 3 / 7, the Zhai 140 C TTCA / Cys / Xyl decay rows, the Zhou 2023 pH 6 / 8 columns); M1's
disulfide record. Never summed, never used to choose.

**Candidate compensating errors, named before the run** (each loss names one or says none is known; each named one
gets a CF entry at the decision):

| expected loss | why | candidate compensating error in the incumbent |
|---|---|---|
| Kang 100 / 120 C FFT and FUR, `kang_100C_MFT_over_FFT`, `zhai_FFT_fold_120_over_100` | the furfural barriers move off 64.08 (stage 0: -1.6 dex on Kang 100 C FFT at B60's constants) | the single 64.08 barrier on `k_tdp_fur` standing in for furfural's measured 147; `k_fur_decay` on the carbonyl-sink family's 143.5 against Yaghmur's 46.5 |
| `zhai_MFT_fold_120_over_100`, Kang MFT rows | the route's barriers and supply both move | one barrier for every route (E11); the instant cysteine-catalysed enolisation (`r_pent_dpo`) standing in for the route (M1) |
| `zhou_pH7_MFT`, `_FFT`, `_MFT_over_FFT`, `_FUR_arp_alone` | 120 C, off the 145 C reference | E5 (no buffer catalysis; an unbuffered ARP + Cys pot on constants from 0.5 M phosphate); E2 (no pH on the trunk) |
| the dimer rows (`zhou_pH7_dimer_over_MFT`, `zhang_fig1_*`) | MFT supply and timing move | the declared ambient OX (1.0 mM); the irreversible dimer (C3) |
| Hofmann fed FFT rows (`fed_tdp_h2s_FFT`, `fed_furfural_h2s_FFT`, `fed_ribose_h2s_FFT`) | furfural barriers | CF-065 (fed sulfide charged as dissolved H2S, biasing `k_fur_fft_hs` low by 0.6 - 0.8 dex) |
| `fed_nf_h2s_MFT`, `fed_nf_cys_MFT`, Whitfield NF rows (C1F-NF) | `k_nf_decay` free | CF-065; the NF + H2S step's unmeasured barrier |
| `cerny_ternary_thiamine_share`, `fed_thiamine_MFT` | the NF / DDP mix changes the non-thiamine MFT | the thiamine route fitted against an NF-route-only network |
| the 145 C fed + H2S rows | the `k_mft_decay` prior | the frozen fast sink (CF-107); the isothermal reading (CR-03) |
| Kumazawa FFT retention | none expected | -- (a loss is a build error until shown otherwise) |
| any row, by up to the numerics row's size | inline vs fixed point | none |

**Spending (tags):** de Roos 2005 Fig. 4 (E10 TEST table) is a constant of every S1 / S2 rival, already spent in M1;
BTBU and DFA as in M1 sec. 8; this design adds Nestle (pot D), Givaudan and Jiangnan's 140 C pots as choosing data for
E6 / E11 (all three are DEV laboratories for E6 / E11 in v2, so no TEST contrast loses independence).

## 11. The adversarial review (sec. 7 item 8)

- **Who:** fresh general-purpose agents with no prior context, never a fork of the author's session, never the author.
- (a) **Design review**, before any build: this file, M1's prereg, review and decisive comparison, the stage 0, the
  data files and dossiers named here. Its brief: the strongest case that this design cannot decide, or decides by
  construction; specifically (i) whether S2's slope band is supported or invented, (ii) whether any band is set from a
  scored unit, (iii) whether any required-agreement form makes a tie inevitable, (iv) whether the Q1 band encodes E6's
  answer, (v) the IS estimator and the ESS rule, (vi) sec. 17's open questions. Then v2, a re-review, freeze.
- (b) **Results review**, before the decision: fit records, probe, decisive and side-effect tables, the build diff and
  the data; not sec. 15 or the author's verdict paragraph. Its task: the strongest case AGAINST shipping the favoured
  cell and against the rule's application. Its brief names the constants every rival holds at values fitted under
  instant supply (sec. 3.2).
- Every objection answered in the record (sec. 18) or the change waits.

## 12. Order of operations, blinding, run-once

0. Prerequisites: this design frozen (sec. 11 (a)); H0's omega and tau on disk (they are, EX-059); V1 - V9.
1. Build, off by default: `sulfur.ttca_exchange` with its bands; the Q1 / Q2 coordinates through the B62 hook in the
   generator; the pot D tracer extended to MP3P; the decisive scorer (IS, quadrature fallback, tie bands, the forms of
   sec. 7.3); the central polish; generator `scripts/generators/fit_m1m2_supply_routes.py`, modelled on
   `fit_m1_supply_combinations.py` (its `Objective`, `failure_accepted`, run-once and freeze checks reused). Unit tests.
2. Feasibility probe (sec. 4.4).
3. Wave 1 fits (R0, C1, C1F), then wave 2 (R0-Q1, M2, C1F-NF); polish; `failed_pot_runs` read on every record;
   comparison rule; series-noise continuation; Laplace; profiles (sec. 9); G1.
4. **Only after every fit has ended:** one scorer run computes every decisive number
   (`kinetic_core_m1m2_decisive_comparison.{md,json}`). No decisive number on an intermediate member.
5. Side effects. 6. Results review and answers. 7. Decision; EX entry; CF entries (sec. 10, sec. 16.3).

**Run-once and freeze checks (M1's and H0's discipline):** `--score` refuses unless the prereg, the generator,
`fitting.py`, `fit_m1_supply_combinations.py` and the scorer are unchanged since the freeze commit (an ancestor of HEAD)
and every fit record carries that commit; a second `--score` refuses unless `--amendment REASON`, which writes
`..._amendment_<k>` with the reason. Amendments are numbered (A1, A2, ...), written before the number they affect, and
say which fields were read.

## 13. Differences from the plan, M1 and the E6 + E11 proposal

1. **"Labelling in 3 labs (Cerny, Givaudan, Hofmann)"** (plan row M2): only Cerny's is labelling. Givaudan (de Roos 2005,
   Laarhoven) and Hofmann are fed-intermediate comparisons (rank 2). Cerny 2003 is Nestle (Lausanne / Amiens), not
   Givaudan or Firmenich.
2. **"Zhai TTCA"** is read as the fed-TTCA MFT / FFT at 140 vs 120 C (the held-out T-law test). The 100 / 120 C cells are
   in-sample (B60's folds; Kang is a reprint). UNVERIFIED that this is what the plan meant.
3. **"In Cao 2017 free cysteine and xylose vanish within 30 min"** (plan row M1, M1.0) is not what the figures show: they
   plateau at 15 - 20 % of the charge (Cys 0.0752 -> 0.056 mol/L, Xyl 0.080 -> 0.069 over 0.5 - 5 h); TTCA + CAMA hold
   90 - 92 %. The reservoir finding stands; the wording needs correcting where it lives (sec. 16.3).
4. **Supply has three levels, not two.** M1 judged S1 only; the review showed S1 is right in kind and wrong in rate.
5. **The measurement model changes from M1's** (the review's three fixes): floors from assays (sec. 6.1), no tau on any
   decisive unit, IS densities on every form, and a tie band from the data's replicate noise (M1 read a tie only from
   sign disagreement). M1's own rule is re-applied side by side for comparability (sec. 8.3).
6. **The earlier E6 + E11 proposal** (2026-11-01) had fixed pass bands (pot D [0.03, 0.15], Givaudan >= 10, Hofmann
   within 3x, D5 <= 0.1, "worse by > 0.3 dex"); this design replaces every one by likelihood comparisons with tie bands
   from the data (sec. 7 item 5: no round-number thresholds). Its barriers were declared; here they are fitted inside
   the same classes (the plan's "bounded by measured classes, fitted together"). Its `k_nf_mft` profile {64, 100, 133}
   becomes a fitted coordinate in [40, 133]. Its held-out H1 / H2 move to side effects.
7. **Pot D's MP3P share** is added as a unit (not in the plan): it is what separates Q1 from the rival Q2.
8. **Q2 (fed NF lost to its own decay)** is the plan's "at least one other explanation" (sec. 7 item 3) for E6; it
   overlaps M1b's E16 (base-catalysed NF decay), which owns the pH law. Q2 here is pH-flat with Mikami's pH 8 rate as
   its ceiling.
9. **The objective is M1's** (B60's 52 + the route block), not an M2-specific one: comparability with M1 and with B60.
10. **Gradient accuracy:** fits polished on a central-difference Jacobian (sec. 4.2); M1 used forward differences
    throughout and central only in Laplace.

## 14. Compute estimate (current engine, <= 8 cores)

**Measured basis** (W1.3, `profile_fit_evaluation_2026-11-05.md` on `origin/w13-profile`, Apple M4, inline, rtol 1e-5):
one serial evaluation R0 3.92 s (52 rows), C1 9.86 s (+ 33 route pots); the machine delivers ~3.5 - 4 serial
evaluation-equivalents however many workers (4 performance + 6 efficiency cores); the FD Jacobian is ~90 % of a fit's
wall; ~80 % of a pot is the Python rate-constant rebuild (52 %) and the two pH solves (30 %); M1's fits took 146 - 280
evaluations (R0, 17 coordinates) and 179 - 226 (C1, 19). So machine-wide wall = serial seconds / ~3.75, whether two fits
share the machine or one runs alone (two 4-worker lanes each run ~2x slower). New coordinates do not change the per-
evaluation time; S2 / Q1 / Q2 add no pots. Assumed iterations: 10 - 15 Jacobians per start.

| step | count | serial h | wall h (/3.75) |
|---|---|---|---|
| build + tests (work) | key, hook, tracer, scorer with IS, polish, generator | -- | **10 - 14 h work** |
| V1 - V9 | reproductions; V5 at 81 nodes on ~10 pots x 2 vectors; V6 6 coordinates x 4 steps x 2 numerics; V7 quadrature | ~4 | ~1 + 2 h work |
| probe | 3 rivals x 256 draws on their claim pots + local search | ~1.5 | ~0.4 |
| fits, forward | 3 starts x (R0 18 x 12 x 3.9 s; C1 20 x 12 x 9.9; C1F 22 x 12 x 9.9; R0-Q1 25 x 15 x 3.9; M2 29 x 15 x 9.9; C1F-NF 24 x 13 x 9.9) | ~12 | ~3.3 |
| central polish | 2 central Jacobians per start (2n evaluations each) | ~3.5 | ~0.9 |
| continuations | ~40 % of starts, half a start each | ~2.5 | ~0.7 |
| series-noise continuation | 4 reservoir rivals x <= 300 evaluations x 9.9 s | ~3.3 | ~0.9 |
| Laplace (central) | 6 x (2n + 1) | ~0.6 | ~0.2 |
| profiles (sec. 9) | expected 2 - 3 bound-active or switched-off per rival (M1: 2 - 4); 3 - 4 points x ~150 - 300 evaluations | ~14 (typical) | ~4 |
| decisive densities | 6 rivals x 12 units: proposal search + IS (ESS-driven, typically 256 - 1 024 draws) on the unit's pots, x 6 thermal readings for D1 / Rt-H; plug-in, linearised, MC; tie bands (cheap) | ~12 - 18 | ~3 - 5 |
| side effects | the benchmark x 7 configurations; in-sample groups; numerics row | ~5 | ~1.5 |
| reviews | design review + re-review + results review + answers | -- | **5 - 7 h work** |
| **total** | | **~60 h serial typical** | **~16 h unattended wall typical (~1 - 1.5 calendar days with the serial steps); worst ~45 h** (every start at its budget and continuation, 4 profiled coordinates per rival, IS at its 4 096 cap on half the units) |

Plus **~17 - 23 h of author work** (build, V-steps, reviews). Wave 2 waits for wave 1's carried optima, so the critical
path is wave 1 (~1.5 h wall) -> wave 2 (~2.5 h) -> profiles (~4 h) -> scoring (~3 - 5 h). Every job `nice`d; at most 8
worker processes (two 4-worker lanes, or one 6-worker lane, W1.3's best count for an 18 - 20 point Jacobian); memory
checked at launch (M1 sec. 4.3). The plan's row ("2 days + ~1 day on >= 8 cores") is consistent. W3's compiled engine
would cut every compute line by the W1.3 estimate (5 - 10x per pot) but is not assumed.

## 15. Expectations (written now; not blind, sec. 0)

| event | P |
|---|---|
| V5 (inline vs fixed point) and V6 (smoothness) pass on the new coordinates | 0.7 / 0.55 |
| every rival compared (sec. 4.2) | 0.55 |
| at C1F's optimum `m_ring` > 0.5 | 0.6 |
| supply contrast: S2 credited / not separated or conflict / S0 credited | 0.45 / 0.5 / 0.05 |
| D1-FFT favours the S0 rival in at least one supply pair (S2 then not credited) | 0.4 |
| route contrast: Q1 credited / not separated or conflict / Q0 credited | 0.25 / 0.65 / 0.1 |
| Rt-D favours Q1 in both pairs | 0.45 |
| Rt-G favours Q1 in both pairs | 0.35 |
| `d_ddp` identified at M2's optimum (Laplace rule) | 0.4 |
| Q1 vs Q2: Q1 credited / not separated / Q2 credited | 0.35 / 0.55 / 0.1 |
| the temperature-law contrast credits Q1 (E11) | 0.2 |
| at least one verdict "conditional on a numerical bound" | 0.5 |
| at least one unit "density not resolved" | 0.2 |
| outcome row 1 of sec. 8.3 (M2 ships) | 0.12 |
| outcome "S2 ships, routes open" (rows 2 - 3) | 0.3 |
| nothing ships | 0.45 |
| a held-out benchmark row gets worse by more than its noise (favoured cell vs R0-shipped (b)) | 0.85 |
| M1 + M2 verdicts by 11-20 (plan sec. 8) | 0.4 |

## 16. Owner-level decisions and registrations

### 16.1 Decisions for the owner

1. **S2's slope band.** Default: m in [0, 1] from a one-proton mechanism, Cao 2017's K, and no use of the scored pot
   (the 2026-10-31 lesson). The brief asked for a band "wide enough to include" Zhang 2025b's timescale. It does include
   it at m = 1 (sec. 3.3), but if review finds the mechanism ceiling unsupported, the only way to keep the band wide
   enough is to set it on Zhang: then the reservoir block is no longer out of sample and supply can be credited only on
   Cao and the route block. Which?
2. **Givaudan as a decisive unit** although its ledger row is `simulable: false`, `condition_uncertain` (37 % water,
   lactate, ~2 M, outside M5's water-activity terms). Default: decisive, with the 0.3 M reading as a required-agreement
   form. Alternative: a reported guard, leaving Rt with two claim units for Q1 vs Q0 (both must favour).
3. **Shipping moves shipped predictions** (every pentose MFT / FFT; `merge-authority` rule): owner sign-off before
   merge, after the results review.
4. **Overlap with M1b.** Q2 uses a pH-flat NF decay with Mikami's pH 8 rate as a ceiling; M1b owns E16's pH law. If M1b
   runs first, Q2 takes M1b's law instead. Default: M2 runs Q2 as written; M1b's verdict supersedes it.
5. **The run waits for this design's review**, not for B64 (the thermal readings stay scored across, as in M1).

### 16.2 Literature needs (to the owner's download queue; lists stay out of the repo)

- A measured pH dependence of 2-substituted thiazolidine-4-carboxylic acid ring opening or formation (the S2 slope).
- Any measured NF decay rate at pH 5 - 6.5 and >= 95 C (Q2's band; E16).
- Cerny's raw 145 C isotopologue data, or pot D at 120 / 145 C (already a W6 request).

### 16.3 Entries for the register session (not written here)

- EX: this design's ID at freeze (next free, checked against origin/main and open branches).
- CF candidates: (a) **Cao 2017 vs Zhang 2025b (same laboratory): one K cannot hold both** (apparent 60 L/mol at 90 C,
  pH ~4.4 vs a ratio climbing to ~18 600 L/mol at 120 C, pH 6.5; about half of Zhang's charge unaccounted, CF-063 /
  CF-076); (b) Zhai 2021's printed TTCA-degradation Ea 80.99 kJ/mol is not reproduced by its own three rate constants
  (~35.5, arithmetic); (c) the plan's M1.0 wording on Cao 2017 (sec. 13 item 3); (d) Zhang 2025b's high-lysine arms plot
  free cysteine (~0.03 - 0.06 mmol/L) below its printed ELSD curve's zero-signal point (0.223), so the curve is not the
  assay's effective floor.
- CR: the floor choices of sec. 6.1 and the noise transfers (Cao's RSD to Zhang; Cerny 2007's 23 % to Cerny 2003 and
  de Roos), each with its sensitivities.

## 17. Open design questions (not settled by the author)

1. **Does anything in the objective identify the DDP route?** No labelling row is in-sample. If `d_ddp` is unidentified
   at every Q1 optimum, it is held at its optimum and the Rt units score a sharp prediction at an arbitrary value; the
   alternative is to move one labelling unit into the objective (e.g. Cerny 2007's 145 C shares), which removes it from
   the decisive set. Which is better is a design-review question.
2. **The S2 slope's split.** The slope is put on the ring opening (formation pH-flat). Putting it on formation (amine
   protonation; ring opening pH-flat) gives the same K(pH) but faster equilibration at high pH, which changes Zhang's
   free-reactant timescale. One coordinate for the split, or a no-refit sensitivity? (pH_ref 4.4 on the label scale is a
   second choice to check: Cao's in-situ pH at 90 C is not measured.)
3. **Is Q2 a fair rival or a strawman?** With its barrier free in [40, 135], NF decay can be fast at 95 C and moderate at
   145 C. Whether Mikami's pH 8 value is a fair ceiling for pH 5 pots is E16's question (M1b).
4. **`k_dpo_ddp`, `k_dpo_c2c3`, `Ea_lumped_formation` held** at values fitted under instant supply (sec. 3.2). Freeing
   them adds non-identified directions; holding them leaves M1's review item M4 open. Is a no-refit sensitivity at their
   B60 Laplace +-1 sd enough?
5. **The Zhai fed-TTCA pots are unbuffered pH 7 water.** The B51 balance moves their pH; S2's `m_ring` and the route's
   pH-flat steps meet that drift. A fixed-label-pH sensitivity (M1's) is printed; should it be required-agreement for
   the T units?
6. **Hofmann's ratio carries an in-sample numerator** (`fed_nf_cys_MFT`). Removing that row from every objective would
   make Rt-H fully out of sample but weakens `k_nf_mft`'s identification. Keep (default) or remove?
7. **The noise transfers** (Cao's 10 % RSD to Zhang 2025b; Cerny 2007's 23 % to Cerny 2003 and de Roos 2005) are the
   largest declared numbers in the measurement model. Is there a better within-laboratory statement on disk for any of
   them?
8. **The ESS rule** assumes the band's estimation precision is the right yardstick for MC error (sec. 7.2). A reviewer
   may prefer a rule on the probability that MC error alone flips the unit's side.
9. **Required-agreement forms vs power.** Every added form can only turn a "favours" into a tie. With 6 thermal readings,
   two noise scales, omega, floors and the FD scheme, a design review should estimate the chance that a true 1-dex
   improvement on Rt-D is read as a tie (a power check on synthetic data, ~1 h of compute at R0's speed), before freezing.
10. **The in-sample group test** (sec. 8.2) uses the rows' declared sigmas, which are B60's fit sigmas, not measurement
   noise; for groups whose sigma is a weight, the "lost" test is an approximation. Should the group test use H0's omega
   instead?
11. **Is S-TTCA a within-lab comparison?** Its K and k come from the same laboratory (Cao 2017, BTBU), but `k_form`,
    `ea_decay` and the route block's Zhai rows are other laboratories'. If a reviewer reads it as absolute across
    laboratories, tau applies and S-TTCA becomes uninformative (it is the one supply unit where S1 and S2 differ most).

## 18. Results

*(appended after the run; nothing above this line is edited after the first number exists)*

## 19. Revision log

- v1 DRAFT, 2026-11-05: first version, written from M1's prereg (v3.3), decisive comparison, amendment A1 and results
  review; the E6 + E11 proposal and its stage 0; H0 (EX-059); W1.3's profile; the dossiers and keyed rows named in
  secs. 2 and 6. Not reviewed.

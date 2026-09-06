# Wave B10 pre-registration — the temperature structure of the sulfur lane (programme step R2(c))

*Written 2026-09-06 before any B10 number existed. Owner's instruction: "take the most sensible
decision long term." This wave was chosen over the oxygen wave (now B11) for the reasons in sec. 1.
Generator to be written: `scripts/generators/generate_kinetic_core_b10_fit.py`. Programme:
`tasks/data_restructure_plan.md`, "Reaction-modelling programme".*

## 1. Why this wave runs first

1. **Every prediction away from 145 C rides on one frozen number.** The sulfur lane carries a
   single `lumped_formation_Ea_kJ_mol` = 64.1 on 37 formation and miscellaneous steps, frozen
   since B9; only the two sink barriers are free (the thiol sink at its 102 ceiling). The
   2026-09-04 diagnosis measured the consequence: the lane is unbiased at 145 C and runs
   ~+2 dex at 100-130 C.
2. **The oxygen probe (B11 prereg, sec. 2) showed the oxidant channels are inert at trace thiol.**
   The oxygen wave needs a new sink structure and two unidentified constants; this wave needs
   no new constant, only a split of one that exists, and has data.
3. **The corpus holds four thiol temperature ladders no fit reads, and they agree on the shape:
   non-monotone.** Kang 2026 (100/120/140 C, TTCA + Cys, 120 min): flat then steep (x1.12 then
   x4.26 for MFT). Meng 2017 (80/95/120 C, fermented soy sauce, 5 and 20 min, n = 4): rising with
   a FALLING apparent barrier (122.8 -> 25.0 kJ/mol for MFT at 5 min). Wang 2026 (85-125 C, fed
   Cys-Amadori + Glu-Amadori, five rungs, digitised): MFT peaks at 115 C and falls, FFT peaks at
   95 C and collapses. Yiltirak 2026 (100/110/120/130 C under a time-compensated cook, SIDA,
   buffered, n = 3): MFT falls 4.0x, FFT rises 1.27x. A single Arrhenius formation barrier cannot
   produce a peak; a formation barrier below a sink barrier can. The structure is what the data
   ask for.
4. **The directional panel carries no sulfur temperature claim at all** (its temperature claims are
   acrylamide, HMF and pyrazines), so the lane's temperature behaviour has never been scored as a
   direction. B10 adds the claims BEFORE it fits, so the pre-wave baseline is on record.

## 2. The data, and which role each takes under the fit / validate rule

Owner's rule (2026-09-03): rate constants, activation energies, fed-intermediate yields, conversions
and WITHIN-STUDY RATIOS fit; end-to-end levels in full precursor systems validate.

| source | what it gives | role in B10 | why |
|---|---|---|---|
| Yiltirak 2026 buffer arm (dossier `yiltirak2026_extraction.md`) | 4 T-t rungs, MFT + FFT, SIDA, stated charge and buffer, one lab | **FIT: 6 within-study folds** (consecutive rungs, MFT and FFT); **VALIDATION: the 8 levels stay on the hold-out panel** and become `in_core_fit`-adjacent | the corpus's only clean thiol ladder; a fold cancels the lab factor and is primary evidence under the rule; the level is not |
| Kang 2026 / Zhai 2023 (in the objective) | 100/120 C MFT, FFT, furfural levels and folds; cys conversions 100/120 | unchanged (fit) | already declared; the 140 C rung stays the B2.x hold-out |
| Feng 2022 (in the objective) | ARP conversion 100/120 C | unchanged (fit) | |
| Wang 2026 (dossier) | five-rung MFT/FFT shape, digitised; pH and time of the series NOT stated; Glu-Amadori co-charged (no core species); SI Table S2 not on disk | **VALIDATION, shape only:** MFT peak in (105, 125) C, FFT peak in (85, 105) C, as directional claims with `evaluable` set by whether the engine can charge the pot | too many unstated inputs for a fit row; the SHAPE is robust to them |
| Meng 2017 (dossier) | 80/95/120 C, 5 and 20 min, MFT + FFT, n = 4, soy-sauce matrix; 80/95 in an open cylinder, 120 in an autoclave | **VALIDATION, ordering only** (`evaluable: false` until a chargeable proxy exists; recorded) | no stated precursor charge; vessel discontinuity |
| Chan & Reineccius 1994 (dossier) | 75-115 C ladders of methional, DMDS, 2-acetylthiophene, six Ea in (81, 137) kJ/mol, non-linear legs | **PRIOR band for the thiol-assembly barrier** (81-137) | none of the three is a core species; the class is the right one |
| Zamora 2013 (dossier) | Ea ladder for carbonyl-amine Strecker products, 27.6-78.0 kJ/mol | prior information only, recorded | downstream products; the dossier itself forbids substituting them into a sink |
| Hofmann 2002 brew 80 C, van Seeventer 50 C | thiol loss at low temperature | unchanged hold-outs (B2.x) | the sink barrier's only low-temperature tests |
| Hofmann 1998 dry-heat 180 C series | one point, confounded (water, T, time) | stays ordinal, unchanged | the reconciliation dossier says the T/water series does not exist |

**What is NOT done.** No level enters the fit. No proxy charge is invented for Meng. Wang's digitised
bars are not fit rows. The Kang 140 C rung stays a hold-out.

## 3. The structural change

**(a) One barrier becomes two, by route.** The 37 keys on the lumped barrier split into
`Ea_sugar_trunk` (pentose / Amadori / hexose entries and the deoxyosone branchings: `k_pent_*`,
`k_arp_*_th`, `k_glc_*`, `k_dpo_*`, `k_tdp_fur`, `k_ttca_*`) and `Ea_thiol_assembly` (every step
that joins a sulfur nucleophile to a carbonyl: `k_ddp_mft*`, `k_nf_mft`, `k_nf_mp3p`, `k_fur_fft*`,
`k_mgo_mp`, `k_ha_mp_mft`, `k_hmp_*`, `k_thi_*`, `k_cys_actz`, plus `k_h2s_loss`). The measured
overrides (`k_cys_thermal` 55.1, `k_dimer_*` 122.2, `k_arp_dpo/tdp` 85.7) and the two sink families
are untouched; the residual decay keys not in a family (`k_dimer_decay`) join `thiol_sink`. The
exact key table is printed by the generator and pinned by its test.

**(b) Both barriers are FREE, with declared bands narrowed by the prefactor rule** (2026-09-04:
holding k(145 C) fixed, 8.0 kJ/mol per decade of prefactor; 12 decades). `Ea_sugar_trunk`: centre
85.7 (Zhang 2026's measured Amadori enolisation), band (40, 135). `Ea_thiol_assembly`: centre
100 (Chan 1994's class), band (55, 145). Both sink barriers keep their B9 bands.

**(c) The ambient oxidant is charged consistently** (B11 prereg finding 2.1): the engine charges
`OX` at the fit's ambient 1.0 mmol/L when no vessel is declared, and from the vessel block when
one is. Effect on every panel row is below 1 % at trace thiol (probe), so this is a consistency
fix, not a modelling change, and the B9 numbers are reproduced to 1e-6 by a unit test when the
ambient charge is set to zero.

**Free set: 23 + 2 = 25** (`lumped_formation_Ea` was frozen in B9; both route barriers are free
in B10). Objective: 54 + 6 = 60 rows.

## 4. The six new rows

Systems: the Yiltirak buffer arm at (100 C, 240 min), (110 C, 120 min), (120 C, 60 min),
(130 C, 30 min); 25 mM ribose + 25 mM cysteine; 0.5 M potassium phosphate, pH 5.5 bench; vessel
3 mL in 20 mL under air, tap water (recorded, not modelled in B10). Rows: `kind =
cross_system_ratio`, consecutive rungs, targets from Table S3 (MFT 3.29/6.88, 2.4/3.29, 1.71/2.4;
FFT 1.46/1.28, 1.68/1.46, 1.62/1.68), `sigma_log` 0.10 (SIDA, n = 3; the printed SDs give
0.03-0.07 dex per rung, and 0.10 leaves room for the unreported come-up time). Anchors quote the
table row verbatim per the dossier.

## 5. Pre-registered tests, the falsifier, and the ship rule

Scored after the fit is frozen; baselines are today's numbers.

- **T1 (structure): the Laplace at the optimum identifies both route barriers** (finite sigma, not
  rank-deficient). Falsifier: either lies in the null space. That outcome is reported as "the
  corpus's temperature contrast does not identify a split" and the wave does NOT ship the split.
- **T2 (in-sample discipline):** no B9 row's |residual| grows by more than 0.3 dex.
- **T3 (Yiltirak levels, now fit-adjacent):** the median fold error of the 8 Yiltirak rows falls
  from 100x to below 10x. This is NOT an out-of-sample claim and is labelled so.
- **T4 (hold-out shapes, strictly out of sample):** Wang 2026 MFT peak inside (105, 125) C and FFT
  peak inside (85, 105) C on the directional panel, if the pot is chargeable; Kang 140 C rung
  direction (MFT 120 -> 140 rises); Hofmann 2002 brew FFT loss at 80 C not worse than today.
- **T5 (leave-Yiltirak-out):** the same fit WITHOUT the six fold rows, reported alongside. What
  the difference shows is what Yiltirak alone taught; if the two optima agree within the Laplace
  sigma the temperature structure was already in the objective and the fold rows cost nothing.
- **T6 (directional):** the sulfur temperature claims added in sec. 6 move from their pre-wave
  score (recorded before the fit) toward agreement; the panel headline does not fall.

**Ship rule.** B10 ships if T1, T2 and T4's Kang rung hold and T3 improves at all. If T1 fails, the
two barriers are re-merged, the consistent oxidant charge ships alone as B10, and the finding is
recorded: the next data item is a two-temperature SIDA measurement in one buffered pot (R7).

## 6. Directional claims added before the fit

`YIL-01` MFT falls monotonically across the compensated ladder (100/4 h -> 130/0.5 h);
`YIL-02` FFT at 120 C / 1 h exceeds FFT at 100 C / 4 h (1.68 vs 1.28, Tukey c vs a);
`WANG-01` MFT peaks at 115 C over 85-125 C; `WANG-02` FFT peaks at 95 C;
`MENG-01` MFT rises 80 -> 95 -> 120 C at 5 min and at 20 min (recorded `evaluable: false`: no
chargeable soy-sauce pot). Fit status: YIL-* become `fit_adjacent` when B10 fits the folds and are
labelled so; WANG-* and MENG-* stay independent. Their pre-wave score is the baseline T6 reads.

**Pre-wave baseline, recorded 2026-09-06 after the claims were added and before any fit:** `YIL-01`
DISAGREE (the shipped lane predicts MFT 79.5 -> 121 -> 163 -> 197 ug/L, rising, against a measured
fall 6.88 -> 1.71); `YIL-02` DISAGREE (FFT 1220 vs 2310 ug/L, the opposite ordering); `WANG-01/02`
and `MENG-01` not evaluable. Panel headline 17/26 -> 17/28; temperature axis 5/7 -> 5/9.

## 7. Forecasts

P(T1: both barriers identified) 0.55. P(B10 ships under the rule) 0.50. P(Yiltirak median below
10x, fit-adjacent) 0.45. P(Wang's two peaks land, out of sample) 0.30. P(a Yiltirak level within
3x) 0.10.

## 8. What is read and what is not

The fit reads its row table (B9's, plus the six fold rows quoted from the Yiltirak dossier) and
nothing under `data/benchmarks/`. The Yiltirak bundles' `holdout_targets` are the same numbers as
the dossier's Table S3 transcription; the fold rows are declared as reading the dossier, and the
fit-target index lists the four Yiltirak bundles as fit-adjacent (`in_core_fit` on the fold basis)
so the scorecard's out-of-sample count excludes them. The hold-out guard runs on the B10 generator
like any other.

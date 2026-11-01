# Cluster 5 verdicts: CF-076 and CF-011 (2026-11-01)

Pre-registration: `prereg.md` (secs. 1-6 before any run; Addendum A, sec. 7, written while the main run was in
progress and before any result was read). Numbers: `results_c5.json` (600 + 80 draws per rival), `search_c5.json`
(differential evolution, maxiter 30), `label_check.json`, `profile_pka.json` (post hoc), `cf011_probe.json`. Code:
`cys_exchange_c5.py` (EX-054's module copied; E0 and controls reproduce EX-054's T1-T7 predictions bit for bit,
`smoke_c5.py`).

## 1. Battery health

- **Controls pass:** NC1 fails T4 and T8b (feasible share 0). NC2 fails T5, T7 and T8a (0). The battery can be read.
- **Power:** T1-T6, T8a and T8b are scored. **T7 (Cao, the only test EX-054 had from a second lab) is dropped:** its best
  non-control pair reaches 0.79 (E3|E5), under the 0.8 rule. The Xie laboratory is now represented by **T8, Zhang
  2025b's blank**, which EX-054 never listed. T8a/T8b power is E0|E4c 0.87 / 0.96 and E0|E5 0.84 / 0.85. Inside the
  slow family the power is low on every single test: E2c|E4c 0.48-0.59, E4c vs E5 / E6 / E1c 0.56-0.76.
- **T3 fails every rival and NC2** (feasible share 0 for all). By the pre-registered clause it is read as a failure of
  the shared skeleton, not used to rank. **S1 (all scored tests minus T3) is therefore the governing set.**

## 2. Search results (best max |z|; > 3 rejects, 2-3 inconclusive)

| rival | primary (T1-T8) | **S1 (no T3)** | S3 core T4-T8 | joint log-marginal S1 (nats) | verdict |
|---|---|---|---|---|---|
| E0 engine | 33.3 (sampled) | 33.3 | 33.3 | -3086 | **REJECTED, now in two labs:** T4 (Jiangnan) and T8b (Xie, feasible share 0) |
| E3 E0 + Cys sink | 33.3 (sampled) | 33.3 | 33.3 | -3159 | **REJECTED** (same as E0) |
| E1c catalysed, Kallen-corrected | 4.18 | 4.02 | 4.04 | 14.3 | **REJECTED**: T8a/T8b ~4.0, T5 3.9 |
| E2c in-situ pH + Kallen k_r(pH) | 3.49 (T3) | **1.55** | 1.39 | 49.3 | feasible, but 3.4 nats behind E4c -> **INCONCLUSIVE** under the rule |
| E4c E2c + catalysis | 3.38 (T3) | **1.45** | 1.43 | **52.7** (best) | **ADMITTED, PROVISIONAL** under the rule (caveats below) |
| E5 secondary-amine catalysis | 5.19 | 5.19 | 4.26 | -4.8 | **REJECTED**: T5, T6, T8a, T2 ~5 |
| E6 storage operator | 5.10 | 5.10 | 4.10 | 3.9 | **REJECTED**, with the storage dose at its band top (f 0.98): T5, T8a, T2 ~5 |
| S4 K released (plain / E5) | 2.60 / 2.64 | -- | 2.11 / 2.16 | -- | hypothesis generator only: a slow exchange with K 0.26-26 L/mmol gets close, but not inside 2 |

**Caveats on the admission:**
1. **E2c and E4c are one family.** E2c|E4c power is ~0.5 on every test. The catalysis term is not identified: c6 is
   2.2 at the S1 point and 0.02 at the S3 point. What the data support is **in-situ pH with Kallen's K(pH) and k_r(pH)
   shape**. Neither amine nor phosphate catalysis is needed.
2. **pKa corner (post hoc profile, `profile_pka.json`).** Both searches put the TTCA ring-N pKa at the band floor,
   5.6. With it fixed: 5.9 still feasible (E2c 1.77, with s_lo 0.63, Kallen's central shape); 6.24 (Kallen's TC value)
   gives E2c 2.43 and E4c 2.40, inconclusive but not rejected. A 2-polyhydroxyalkyl substituent plausibly lowers the
   pKa below TC's, but TTCA's own pKa has never been measured.
3. **Latent pH.** Zhai 2019's pH trajectory is unmeasured; it is a banded latent (6.0-7.4 falling to 5.0-6.5). Hou
   2026 (7 -> 5.05 in 20 min, unbuffered, same laboratory network) makes the band realistic, not proven.
4. **Isotope check at face value (S5a, sensitivity only).** Zhai 2021 Fig. 7a, read here for the first time, rejects
   the E4c S1 point at z 6.3: it over-predicts 13C5-TTCA by 0.4-0.9 dex at pH 7-8. The E2c S1 point is near the edge
   (z 2.6), and no sampled draw of any rival passes (minimum z 2.2-2.9). The face-value reading is doubtful because the
   Methods quantify 13C5-TTCA as a difference between two ~10 mM pots. The resolution-limited bound (S5b, <= 0.3 mM) is
   passed by both feasible points (0.18 / 0.04 mM).
5. **Laboratories.** The Jiangnan laboratory supplies T1-T6, and the Xie laboratory supplies T8 only (one pot, one
   panel). The cross-lab check that the protocol requires is T8, and it is passed.

**Expectations vs outcome:**
- 6 of 9 came out as stated. Controls failed as designed (0.95). T8b was powered (0.9). E1c was rejected (0.85).
  E5 was rejected (0.6). E6 was rejected (0.75). The CF-011 storage bias was < 30 % (0.9).
- 3 of 9 were wrong:
  - "T8 fails every rival": it does not (I gave 0.65). In-situ pH reaches Zhang's K through K(pH).
  - "E2c rejected": it was not (I gave 0.85).
  - "No rival admitted": one was (I gave 0.8).
- The corrected k_r(pH) sign is what changed E2's fate. EX-054's acid limb (a > 0) was the wrong sign on the N-H
  analogue.

## 3. Task 1: corrected bands (summary; full text in prereg sec. 2 and `kallen1971_correction.md`)

- **TC k_r at 25 C**, from Kallen 1971b's k_f and Kallen 1971's K: 5.4e-5 /min at pH 5.6 (t1/2 ~9 d), 1.4e-5 at 4.5.
  It rises with pH up to 5.6 and falls above 6.6. The dossier's 1.5-5e-4 /s was the relaxation, 120-270x too fast.
- **k0(100 C) band** stays [1e-5, 1] /min on a rewritten basis: TC Arrhenius-transferred (1.4e-3-0.18 at pH 5.6) times
  phi in [0.03, 5]. Phi, the TC-to-TTCA transfer factor, is not identified by any data on disk.
- **c_aa:** 0.01-3 /M at pH 6, rising about x8 per pH unit. None below pH ~5.5. Saturating at x3.5-5.
- **c_P:** 90-800 /M HPO4 2-, under the same gate and ceiling.
- **k_r(pH):** s_lo 0-0.85, s_hi 0-0.8.
- **E5:** c_N 0.5-500 /M, from morpholine's 55-70 /M at 10 C, transfer x0.01-8. k_open(120 C) 0.08-0.16 /min, from
  Zhang's lysine arms.
- **E6:** storage dose <= D0 = 1.225e-3 L/mmol, from Zhai 2019's t = 0 point.
- **Transfer assumptions:**
  - TTCA releases by the same C-N cleavage as TC (N-H ring);
  - one Ea from 25 to 140 C;
  - ionic strength ignored;
  - formaldehyde to aldose only through phi;
  - catalysis ratios are temperature-independent.
- **Analogue tension, stated:** at Ea 64 and phi 1, TC itself releases ~1e-2 /min at 100 C, 10x above the dilute
  ceiling. The feasible family needs k_r45(100 C) ~2-8e-4 /min, i.e. phi ~0.1 or Ea near 40.

## 4. CF-011 probe 1 (forced-input version; no engine)

Apparent k_form = [dCAMA + k_cd integral CAMA] / integral TTCA, at both edges of k_cd:

| pot | T | medium, pH | k_app (/min) |
|---|---|---|---|
| Cao 2017 | 90 C | 0.2 M phosphate, 5.5 -> 4.2 | 3.0-3.2e-3 (30-120 min); 8.8e-4 over 0-30 min |
| Zhai 2019 | 70 / 90 / 110 / 130 C | water, label 7.4 | 4.1e-4 / 4.0-5.8e-4 / 2.1-4.3e-3 / 0.9-1.8e-2. Own Arrhenius Ea 69 kJ/mol, scatter +-0.3 dex |
| Zhang 2025b | 120 C | 0.2 M phosphate, 6.5 | 0.9-2.0e-2 (10-25 min); 2.3e-3-1.6e-2 from the CAMA plateau |

- **Ratio Xie / Zhai line:** x3.1-3.3 at 90 C and x1.6-3.6 at 120 C. A roughly constant factor within Zhai's scatter.
  This contradicts the "barrier pinned at the 40 kJ floor" part of the earlier reading: Zhai's own forced-input Ea is 69.
- **The two Xie pots bracket Zhai in pH (4.2-4.7 and ~6.3-6.5 vs ~5-6) and in concentration (400 and 40 mM vs 83 mM),
  and both are faster.** No monotone pH or concentration law gives that. That leaves phosphate (0.2 M in both Xie
  pots) or a laboratory / quantification effect (Jiangnan ELSD at the low end, 0.2-5.6 mM CAMA; calibration
  unprinted).
- **Storage operator:** inflates Zhai's TTCA by <= 5-15 % (<= 0.07 dex of the x3).
- **Guan 2012 cuts both ways.** It is a third laboratory showing phosphate speeds Amadori-compound chemistry at
  110-130 C. But it measured decay, which would lower CAMA. Phosphate must therefore speed formation more than decay.
- **Unexplained, noted:** Cao's first 30 min run at Zhai's rate (8.8e-4 vs 9.8e-4); the excess appears after Cao's pH
  has fallen below ~4.7.

**CF-011 verdict:** **EXPLAINED in part** (pH and concentration excluded, storage negligible, barrier confound not
needed). P(phosphate catalysis of formation) 0.5, P(Jiangnan ELSD quantification of low CAMA) 0.3, other
(laboratory, heat-up) 0.2. **Deciding data:**
- buffer arm C: one pot, water vs 0.05 / 0.2 M phosphate at fixed pH, CAMA by LC-MS;
- or, cheaper, Zhai's ELSD calibration for CAMA (author request).

## 5. CF-076 verdict

**EXPLAINED, provisional (P 0.45): in-situ pH with Kallen-shaped K(pH) and k_r(pH).**
- **The pH difference does the work.** The unbuffered concentrated pots (Zhai 2019, pH ~6-7 falling) and Zhang's
  phosphate pot (pH ~6.5) sit at higher K and faster release than Cao (pH 4.2-4.7) and the drifting dilute pots
  (~4.5-5.5).
- **No catalysis term is needed.**
- **Rejected outright:**
  - the engine exchange, now in two laboratories (T4 and T8b);
  - Kallen-corrected catalysed exchange;
  - secondary-amine nucleophilic catalysis;
  - the storage operator.
- **Why "provisional":**
  - it is feasible only from pKa_T <= ~5.9; at 6.24 it is inconclusive (2.4);
  - E2c vs E4c is not separable;
  - the Zhai 2019 pH is latent;
  - the face-value isotope reading disfavours the E4c point;
  - T3 is a skeleton failure that nobody passes.
- **Probabilities:**
  - P(in-situ pH / K(pH) is the main cause) 0.45;
  - P(slow exchange with a much larger, pH-dependent K, the S4 direction) 0.2;
  - P(quench / storage) 0.05;
  - P(catalysis, any kind) 0.1;
  - P(something unlisted, e.g. an unmeasured Cys + Xyl adduct pool, CF-063) 0.2.

**Task 4: what decides it:**
1. **TTCA ring-N pKa (titration, 25 and 90 C):** the cheapest datum. Below ~5.9 the pH family is feasible; at >= 6.2 it
   falls to inconclusive.
2. **Author requests:**
   - **Zhai 2021 Fig. 7a (Jiangnan):** the direct m/z 259 / 254 area ratio of TTCA in the 13C5-xylose pots. If the
     plotted 3-20 uM are real, the dilute re-formation sits 0.4-0.9 dex below the admitted point at pH 7-8.
   - **Zhai 2019 (Jiangnan):** pH of the cooled samples; storage time and temperature before HPLC; any re-injection.
   - **Zhang 2025b (Xie):** final pH.
3. **Designed, extending E4.** The core pot is TTCA at 10 vs 100 mM, 100 C, with a 13C5-xylose spike and a logged or
   pH-stat pH. Three contrasts separate the rivals:

   | arm | pH | condition | decides |
   |---|---|---|---|
   | buffer | 6.5 | +-0.2 M phosphate | catalysis |
   | buffer | 4.5 | +-0.2 M acetate | catalysis |
   | amine | -- | +-50 mM morpholine vs N-methylmorpholine | E5 |
   | handling | -- | immediate 20x cold dilution vs 4 C storage for 0 / 4 / 24 h | E6 |

   The pH-stat contrast (4.5 vs 6.5 at fixed concentration) decides E2c directly.

**For the model:**
- The engine's k_ttca_cys (0.0286 /min at 100 C, pH-flat) is rejected in both laboratories.
- The feasible family puts release at k_r(100 C, pH 4.5) ~2-8e-4 /min, with K and k_r following Kallen's pH shape.
- That is a direction and a declared band for the A2 writer, not an estimate. No engine change before a red team and
  item 1 above.

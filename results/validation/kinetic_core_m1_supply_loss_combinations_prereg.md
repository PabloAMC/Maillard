# M1 pre-registration: thiol supply, judged as R0 against C1; the loss law specified, not run (v3.2)

*Plan step M1 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, row M1; sec. 7 items 6 - 9). Case:
[M1_thiol_loss_oxidant_ph_drift.md](../../docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md) (items 1 - 5; not
edited here). Conflict: CF-107 ([conflicts_register.yml](../../docs/reference/conflicts_register.yml)). M1.0: EX-060,
[kinetic_core_m1_preflight_checks.md](kinetic_core_m1_preflight_checks.md).*

- v1 and v2 written 2026-11-04; v3, v3.1 and v3.2 the same day after the design review and its re-review (sec. 16),
  all BEFORE any M1 build, fit or score. No engine was run with a changed switch or constant to write any version.
  Nothing above the "Results" line (sec. 15) is edited after the first M1 number exists; changes before that are
  versioned and logged in sec. 16.
- Engine state today: main at dcbf8348, shipped calibration B60, configuration digest `40b746476905786c` (M1.0).

## 0. Seen by the author

Read while writing this file. **The expectations in sec. 13 are not blind.**

- **M1.0 (EX-060), all of it:** G1 at B60 (oxidant use <= 0.009 %, phi = 1 rise <= 0.08 %); R0's Schieberle 2000
  ratios (MFT -1.09 / -0.33 / -1.42 dex on 360/30, 720/360, 720/30; FFT +0.36 / -0.09 / +0.27; MFT peak 2.6 h); R0's
  D2 tables (Cao 2017, Zhang 2025b, Zhai 2019).
- **B9 / B16 numbers in the records:** B16 prereg and ship rule (Schieberle fold residuals at both ceilings; T3
  Yiltirak median 198x -> 14x / 42x; T4 TTCA rows; the 6 h turnover). B16's addendum gives B9-engine shapes on
  **Wang 2022 (F28#25)** and its T6 gives B9-engine values on **Liu 2023 at 168 C (F28#26)**: both are v2 E12 TEST
  contrasts. They are B9 numbers, not B60, and they are recorded here as seen.
- **The route at B60's constants (= C1's start vector), seen at start:** B59b Amendment 3 (Schieberle MFT over 30 min
  6.1 / 84 / 73 at 60 / 360 / 720 min; Nestle 110 C shape; Zhang 2025b TTCA +0.15 to -0.37 dex; Hofmann pH 7 MFT -1.6
  dex; Zhai 2019 Cys-Amadori over-made 2 - 5x; Kang fed-TTCA MFT down 260x); the route sub-model no-refit probe
  ([kinetic_core_route_submodel_probe_with_b58_with_b59b_with_b59a_with_b61.md](kinetic_core_route_submodel_probe_with_b58_with_b59b_with_b59a_with_b61.md),
  per-data-set table only: Cao RMS 0.03 - 0.31 dex; Zhang 2025b TTCA -0.10 / CAMA +0.66 / Cys +0.55 / PENT +0.44 dex).
  **The design review's probe (v3):** at C1's start the Schieberle residuals are MFT 360/30 +0.49 dex and FFT
  360/30 +1.58 dex (R0: -1.09 / +0.36), so D1's two species point opposite ways at the start. **The re-review (v3.2):**
  R0's free pentose is predicted 0.00e+00 (the integrator floor, atol 1e-12 - 1e-14) against ~60 mM measured.
- **The stage-1 route calibration**
  ([kinetic_core_cysteine_route_calibration_prereg.md](kinetic_core_cysteine_route_calibration_prereg.md)): its
  constants, its checks (Zhang 2025b, Zhai 2019 Fig. 3B, Zhao 2017, pH), and its flavour-level diagnosis (Schieberle
  12.7 / 345 / 171 over 30 min).
- **R3En** ([kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.md](kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.md)):
  residuals by data set, including Schieberle MFT / FFT, Wang 2026b MFTD / FFTD and Zhou 2022 MFTD (after a free
  scale), the declared hold-outs (Yiltirak, Hofmann T1 / T2, Zhai 140 C) and the score-only rows.
- **B60:** its prereg, its fit report (`kinetic_core_b60_fit_report.json`: start 0 cost 27.132 in 246 evaluations and
  7601.5 s, start 1 cost 27.671 in 298 evaluations and 8600.2 s, both status 2; active bounds `k_dimer_mft` upper and
  `ph_acid_yield_per_sink_event` lower) and `residual_by_row` for all 52 rows (e.g. `zhou_pH7_dimer_over_MFT` -2.71,
  `zhang_fig1_cys_dimer_over_MFT` -1.99, `kang_120C_MFT` +2.51 sigma); B58's ship rule; B27's G1 / G2 as quoted in
  the case.
- **B64 and the thermal readings (v3.2):** `kinetic_core_b64_prereg.md` (the readings, their probabilities, the probe
  at B60's constants, its independent review's item 3: later waves run on every reading); `choices_register.yml`
  CR-03 / CR-04 (their alternatives and sensitivities); the C5 entry of `tasks/data_restructure_plan.md` (the inline
  coupling: within 0.002 dex of the node fixed point on five probe pots; one objective evaluation 27 s -> 7.5 s).
- Register entries EX-024, EX-027, EX-050 (and its README), EX-055 (as quoted in EX-050), EX-060; CF-055; CF-107.
- **Measured values only (no engine number):** Zhou 2023 Table 1 at pH 6 / 7 / 8 (MFT, FFT, MFTD, FUR, with SDs);
  Whitfield 1999 Table 1 row 48 (bis(2-methyl-3-furyl) disulfide 2 / 3 ug Cys arm, 7 / 2 ug H2S arm) and the dossier's
  7.5 % / 35 % shares; Grosch 1992 Tables II and VIII with their SDs (grosch1992 dossier; seen, used for nothing in
  M1 since v3.1); Hofmann 1998's error footnote (triplicates within 10 %); Guth 1995's untrapped arm (k_net <= 4e-4
  /min); Zhang 2025b Fig. 2g - h's lowest printed points and its printed ELSD calibration curve; the F28 atlas cells
  #1, #2, #9 - #12 and F27's Grosch cells; the de Roos 2005 keyed-table header (no value).
- **v3 / v3.2 additions:** the design review's and the re-review's findings; `route_submodel.yml`'s pot list
  (temperatures, species, row counts); `split_v2_per_change.yml`'s E6 / E7 / E11 / E12 / E18 lists; the
  FIT_HOLDOUT_DECLARATION rows on the Zhou 2023 and Zhang 2024 dimer fractions.
- Not read: any engine number on the Zhou 2023 pH 6 / 8 disulfide shares, on Zhang 2024's free-MFT ratio, on
  Whitfield's disulfides, or on Grosch 1992 Table VIII.

## 1. Question

**Does the cysteine-pentose route's reservoir explain the thiol data the shipped structure cannot?** M1 gives one
verdict, **supply: R0 (the shipped structure) against C1 (the TTCA / cysteine-Amadori route)**, with three possible
outcomes (C1 favoured; R0 favoured; not separated / conflict), on two pre-declared blocks of decisive data (sec. 5):
the reservoir (where supply is credited) and the thiol shape (which C1 must not lose).

- **The loss law (C1 vs C2: does free cysteine slow the thiol loss?)** is not decided in M1: C2 is **specified, not
  run** (sec. 2.3). No out-of-objective development datum discriminates the two; in-sample cysteine contrasts might
  identify K_c (review M3), which a later step can test. The question is registered as a data need for plan W6 (sec.
  14).
- **C3 (the disulfide side: Guth's release and an osone oxidant)** is withdrawn from M1 as a registered follow-up
  (sec. 2.4).

## 2. The rivals, exactly

### 2.1 Config keys

| rival | switches ON (Config key) | new step(s) | removed step(s) | fitted coordinates | in M1 |
|---|---|---|---|---|---|
| R0 | none (the shipped configuration) | -- | -- | 17 (sec. 4.2) | fitted, compared |
| C1 | `sulfur.cys_amadori_route` (`enabled = True`; it also turns on `r_ttca_form` by detailed balance with `sulfur.ttca_formation`'s K) | `r_ttca_form`, `r_cama_form`, `r_cama_dpo`, `r_cama_tdp` (exist, B58 / B59b) | `r_pent_dpo`, `r_pent_tdp` (retired by the route) | 19 | fitted, compared |
| C2 | C1 + `sulfur.cys_suppressed_thiol_loss` (a NEW key; not built in M1) | none: a factor on `ch_mft_decay` and `ch_fft_decay` | -- | 20 | **specified, not run** |

**Every M1 computation runs on `ph.coupling = "inline"`** (sec. 2.5, numerics). **R0 is scored in three forms.**
*R0-shipped (a)* is B60 on the node fixed point, the frozen shipped record, printed only. *R0-shipped (b)* is B60 on
inline, like for like; wherever a shipped vector enters a computation (the probe's reference) it is (b). *R0* is the
shipped structure refitted on its own objective (B60's 52 rows) on inline with the two fit-side changes C1 also gets:
the dimer ceiling (sec. 2.5) and the `k_mft_decay` prior (sec. 2.6). The verdict uses R0 (or its declared fallback,
sec. 4.3). Refitting R0 changes no switch; if R0 is favoured, its refit constants are recorded, not shipped (M2's refit
owns the calibration).

### 2.2 Constants of the route (C1; C2 would share them)

| constant | step | status | value, band, source |
|---|---|---|---|
| K, TTCA equilibrium | `r_ttca_form` = K x `k_ttca_cys` | **declared, measured** | 0.060 L/mmol; band 0.045 - 0.075 (Cao 2017, `TtcaFormation`). Held at 0.060. pH- and T-flat (B58's declared assumption). Measured on the Cao 2017 points, which is why Cao Fig. 4 is not decisive (sec. 5) |
| **`k_ttca_cys`**, ring opening TTCA -> Cys + PENT | `r_ttca_cys` | **FITTED** (in the structure, not frozen; CF-107) | start B60 (log10 -0.579 at 145 C); bounds +-3 decades about the start (stage-1 convention). Its Ea stays on `Ea_lumped_formation` (64.08 kJ/mol, B60, held). Check, not a bound: the Cao dossier's floor >= 7e-3 /min at 90 C |
| `k_form`, TTCA -> Cys-Amadori at 90 C | `r_cama_form` | **declared** (anchored on Cao 2017's blank, B59b V2) | 3.29e-3 /min; probe band +-0.44 dex (the anchor vs stage 1's 1.2e-3) |
| **`ea_form_kj_mol`** | `r_cama_form` | **FITTED** | start 100.9 (Zhang 2026, glucose system: a transfer B59b Amendment 3 found too steep for pentose); bounds 40 - 135 kJ/mol (`sugar_trunk` class band, as stage 1). Needs a field band in `config._BANDS["sulfur.cys_amadori_route"]` before any fit (the key refuses an unbanded field) |
| Cys-Amadori decay k(pH) at ~100 C | `r_cama_dpo`, `r_cama_tdp` | **declared, measured** | de Roos 2005 Fig. 4 table; above pH 6 `extend_last_slope` (CR-14; stage 1 found the two arms identical). A keyed E10 TEST table: spent, sec. 11 item 3 |
| `ea_decay_kj_mol` | same | **declared** (transfer) | 85.7 (Zhang 2026); probe band 40 - 135 |
| `dpo_share` | same | **declared, measured** | 0.30 (Zhai 2019; alternatives 0.28 / 0.32, CR-12); `split_law = "none"` |

**The route's temperature and pH dependence (declared; review MINOR).** Of the route's five steps, one barrier is
fitted (`ea_form`), one is a transfer (`ea_decay` 85.7), `k_ttca_cys` rides on the lumped formation barrier (64.08,
held), and K is T-flat and pH-flat. Without the Zhai >= 140 C pots (sec. 4.1) the route block spans 70 - 130 C, and
B60's 145 C rows exercise the route by extrapolation through these declarations. The pH law is the de Roos table on
the Cys-Amadori decay only; K and `k_ttca_cys` are pH-flat across the route block's pH 5.5 - 8 pots (58 % of the
route rows are at pH >= 7), while the decisive data sit at pH 5.0 and 6.5. Both are reported, not tuned (secs. 7, 8).

### 2.3 C2, specified, not run

| constant | status | form, value, band, source |
|---|---|---|
| `K_c`, the cysteine level that halves the loss | **would be FITTED** | `ch_mft_decay` and `ch_fft_decay` run at k / (1 + [Cys]_free / K_c). [Cys]_free is the free-cysteine species, not TTCA- or CAMA-bound. log10 K_c (mmol/L) start 1.5, bounds 0 - 3 (1 - 1000 mM). No measurement pins it at process temperature: van Seeventer 2001 (+150 mM Cys) is a 50 C storage cell and cannot set a processing law (atlas rule) |

The factor would not act on `ch_thiolate_loss_*` (the Kumazawa-anchored step) or on the dimer; the `k_mft_decay` prior
(sec. 2.6) would then be a prior on lambda at zero free cysteine. Form: EX-050's cysteine factor without the
electrophile pool E (sec. 11 item 7). **Why it is not run (owner, v3.1):** no out-of-objective development datum
discriminates C1 from C2 (sec. 5.1), and a C2 fit costs compute for no M1 verdict. **What that costs M1 (review M3):**
B60's rows include within-study cysteine contrasts (the Zhou 2023 ARP vs ARP + Cys pots at pH 7, Zhang 2024 Fig. 1's
cysteine vs cystine arms, Whitfield's Cys vs H2S arms) that could pin K_c in-sample, and a C2 (or the case's R2) could
then change D1. So a D1 gain of C1 over R0 cannot be credited to supply (sec. 7). The key is not built in M1; the
specification stands for the W6 data need (sec. 14) and for M2, which may start from either loss law.

### 2.4 C3, withdrawn: the registered follow-up "C3: the disulfide side"

v2's C3 (Guth 1995's release declared from the trapped arm, plus a first-order-in-thiol oxidation by the
deoxypentosone pool) is withdrawn. The design review (B3, B4, M6, M7) showed it was built to fail and that its
nesting claim was false:

- Guth's 0.0165 /min is the **trapped** arm (4-vinylpyridine, pH probably 7 - 8). The untrapped arm conserves MFT-MFT
  (k_net <= 4e-4 /min). The declared 0.57 one-thiol branch makes the reservoir a net sink.
- A first-order-in-thiol law with an MFTD product contradicts its own mechanism: with cysteine at mM and MFT at uM, the
  second thiol is cysteine, and the product is MFT-S-S-Cys.
- Its pool was model-internal (TDP 1.6 -> 40 mM in the review's probe), and D6 was set by the declared ambient OX
  (1.0 mM).

**The follow-up needs, before it is pre-registered:**
1. Guth's **untrapped** arm as the release's data, with hydrolysis to RSH + RSOH and recondensation;
2. an oxidation product split by thiol abundance (MFT-S-S-Cys where cysteine dominates), or a law second order in MFT;
3. the measured DPO / TDP series (`deoxyosone_submodel.yml`, Zhai 2023 Fig. 2) in its objective, so the pool is
   measured, not model-internal;
4. D6 (Whitfield's fed-norfuraneol arms) dropped, or the ambient OX carried as a declared sensitivity.

v2's D4 - D6 (the disulfide shares, the cystine-arm free-MFT ratio, Whitfield's arms) move to the side-effect record
(sec. 8), where they are recorded for R0 and C1 and become the follow-up's decisive data. The B17 release
(`k_dimer_release`) stays inert in every rival.

### 2.5 Held fixed, freed, the dimer ceiling, and the numerics (R0 and C1)

- **Free, from B60's start:** the 16 coordinates B60's Laplace marks identified (B59b Amendment 2's rule):
  `k_tdp_fur`, `k_nf_mft`, `k_nf_mp3p`, `k_mgo_mp`, `k_fft_decay`, `k_fur_decay`, `k_pent_caramel`, `k_cys_thermal`,
  `k_thiolate_loss`, `ph_arp_secondary_ammonium_pKa`, `k_osone_decay`, `k_ttca_deg`, `Ea_decay_thiol_sink`,
  `k_dimer_mft`, `k_arp_dpo`, `k_arp_tdp`. B60's bounds, except the dimer's upper bound (below).
- **Held at B60's values:** the 7 unidentified at B60 (`k_glc_ha`, `k_dimer_fft`, `ph_acid_yield_per_sink_event`,
  `k_fur_fft`, `k_dimer_decay`, `k_thiol_decay`, `Ea_decay_carbonyl_sink`) and the 24 other frozen keys of B60's
  report, except `k_mft_decay` (now under the sec. 2.6 prior, both rivals) and, in C1, `k_ttca_cys` (fitted).
  `k_pent_dpo`, `k_pent_tdp` drop out with their steps in C1. **C1 therefore carries 31 constants fitted under R0's
  instant supply**, among them `k_dpo_ddp`, `k_ddp_mft`, `k_dpo_nf`, `k_dpo_c2c3` and `Ea_lumped_formation` (review
  M4; named again in secs. 9 and 14).
- **Switches at their shipped defaults, except the coupling:** `ph.node_grid` and `ph.fixed_point` unused under inline;
  `ph.b51` all on; `ph.trunk_solutes_in_balance` False and `ph.trunk_terms` False (the pH balance is M1b's);
  `sulfur.flux_budget` closed; `panel.thermal_history` isothermal for every fit (sec. 6, the thermal readings); B59a
  (`sulfur.strecker_sulfur`), B61, B62, B63, B7a's split law, `sulfur.tdp_own_decay`, the pentose-amine entry and every
  other sub-model switch off. (R3En ran with B59a and B61 on; M1 does not.)

**The numerics (owner question 2026-11-04; re-review addendum B).** `ph.coupling = "inline"` (the charge balance
solved at every right-hand-side call, warm-started Brent) for every fit, probe, Laplace, profile and score, in R0, C1
and R0-shipped (b) alike. It is a numerics change (kind 1): recorded within 0.002 dex of the converged node fixed
point on five probe pots and ~3.6x faster per evaluation (data_restructure_plan C5; EX-056 and H0 fit R3En on it).
- **Tolerances:** fits, probe, Laplace and profiles at integrator rtol 1e-5; decisive data and side effects scored at
  rtol 1e-8; the engine's atol (printed); the Brent solve at `ph_state`'s default xtol 1e-10 pH, ionic strength to
  `IONIC_STRENGTH_TOL_M`; the compiled pH kernel if present (printed). The 17 / 81 node settings no longer apply.
- **Checked before any fit** (V7, V8, sec. 4.4): agreement with the converged fixed point on the decisive pots and on
  B60's 52 rows, and smoothness of the cost (most of B64's inline members ended status -9).
- **Shipping:** if C1 ships, the coupling flip either ships first as its own kind-1 change (V7 supplies its diff) or is
  declared part of M1's ship. Never silently.

**The dimer ceiling (review B4).** `k_dimer_mft` ends B60 on its upper bound, log10 0.5 at 145 C
(`FITTED_SULFUR_BOUNDS_LOG10K`). That bound is NUMERICAL, and its stated reason ("a half-life of well under a second")
holds for a first-order step, not for this one: `ch_dimer_mft` is 2 MFT + OX -> MFTD, rate k theta [MFT]^2 [OX], so at
k = 10^0.5 (L/mmol)^2/min, 1 uM MFT and 1 mM OX the per-MFT rate is ~6e-3 /min, slow. The bound therefore holds the
disulfide below the measured shares. In both rivals it is replaced by a physical ceiling:

- The rate-limiting step of thiol oxidation is bimolecular (thiolate + oxidant); the second thiol's attack is fast.
  The true rate is at most k_enc theta [MFT][OX], with k_enc the encounter limit. The engine's form equals
  k theta [MFT]^2 [OX], so k <= k_enc / [MFT] (theta cancels).
- k_enc = 8RT / (3 eta): at 145 C (eta of water ~0.19 mPa s) 4.9e10 L/(mol s) = 2.9e9 L/(mmol min).
- [MFT] is capped at 1 mmol/L, a declared value >= 10x above any MFT the objective's pots reach (UNVERIFIED; V5 prints
  the panel's maximum at B60, and if it exceeds 0.1 mmol/L the cap becomes 10x that maximum).
- **Ceiling: log10 k_dimer_mft(145 C) <= 9.47** (L/mmol)^2/min. Above 145 C the Arrhenius factor of the measured
  122.2 kJ/mol outgrows k_enc, so the build takes the minimum over the objective's pot temperatures and prints it
  (one number; at 145 C it is 9.47).
- The ceiling is valid but loose (re-review MINOR 7): it removes a numerical bound, it does not constrain the step; the
  in-sample share rows decide where `k_dimer_mft` goes. Integration at large k is stiff: V5 times it, and sec. 12
  carries a stiffness allowance.
- The same ceiling replaces the numerical bound of `k_dimer_fft`, which is held at B60 in both rivals, so the change
  is inert for it. The lower bound (1e-10, "the step does not happen") is unchanged.
- This is a fit bound, not a model constant: the shipped engine is unchanged by it.

### 2.6 `k_mft_decay`: a declared weak prior, both rivals

v2 banded `k_mft_decay` on Grosch 1992 Table VIII. The review (M1, M2) showed that table is held-out TEST data (Table
VIII's reflux values are Table II nos. 1 and 4: E18's F27#15 and the keyed `heldout_grosch1992_tab2.yml`, CF-089), that
the band is an inference from a thiamine pot (another matrix) with a chosen width, that its centre is C2's reading, and
that the band and v3's D2b came from one table. Since v3.1, M1 uses Grosch for nothing.

- **No development-only source of the free-thiol loss rate at process temperature exists on disk.** Searched: the
  evidence map (2.4: "SINGLE plus a bound"), the sulfur ledger and the thiol dossiers. Grosch 1992 is held out; the
  Nestec patent (US 10,076,126; lambda(110 C) <= 0.02 - 0.05 /min) is not measurement-grade (P 0.15) and is a 0.25 M
  cysteine pot, so it would build C2's term into R0's and C1's prior; Liu 2023b (MFT halves 20 -> 60 min at 168 C) is
  E12's F28#26 and a keyed held-out table; Hofmann 2002's 80 C brew is a hold-out row; every other thiol-stability
  source is at 22 - 50 C (Blank 2002, Charles-Bernard 2005, van Seeventer 2001, Dreher 2003, Zhou 2025, Sun 2019),
  which the atlas rule does not transfer to processing.
- **The prior:** Gaussian on log10 lambda_mft(100 C), centre **log10 0.0406 /min = -1.391** (the shipped value, B60),
  sd 1.0 dex. The centre is the incumbent's, so the prior moves no rival away from R0's value; the same term enters R0
  and C1. **Physical basis:** the loss is covalent capture by the pot's electrophiles plus oxidation; no measurement
  sets its rate in these pots, so the prior states only the order of magnitude (2 sd: 4e-4 to 4 /min). It is not
  blind: the author has seen Grosch's inferred 0.03 - 0.05 /min (sec. 0). The hard bounds stay
  `FITTED_SULFUR_BOUNDS_LOG10K`'s; their upper end (log10 0.5 at 145 C, a first-order step for which the module's
  numerical argument holds: t1/2 ~14 s against holds of >= 20 min) lies 0.36 dex above B60's value, so lambda(100 C)
  cannot exceed ~0.094 /min at Ea 101.25; an optimum on that bound is profiled (sec. 7).
- **The identifiability cost (corrected in v3.2, re-review B1d).** At 1 dex the prior constrains nothing the data
  constrain, and the sink can trade against the supply in-sample (EX-024's pattern). v3 claimed that a poorly fixed
  sink would make a rival "predict D1 less sharply and score lower". **That is wrong where both rivals miss D1 by many
  sigma:** there a wider predictive density scores HIGHER, so parameter width would reward the looser rival. v3.2
  therefore (i) holds unidentified coordinates at the optimum, (ii) uses the linearised predictive (sec. 7), in which
  lambda's width enters only through its Laplace covariance, (iii) prints the plug-in density beside it, and (iv)
  decides by the sign pattern of the units (sec. 7), so a magnitude bought by vagueness cannot decide. lambda's Laplace
  sd is reported per rival.
- **The fast 145 C sink is kept, stated plainly.** With lambda(100 C) centred at 0.0406 /min and
  `Ea_decay_thiol_sink` at 101.25 kJ/mol (0.75 below its 102 ceiling), lambda(145 C) = 33.5 x lambda(100 C) ~ 1.4 /min
  (t1/2 ~30 s). That is B60's fast 145 C sink, the anchor CF-107 names. The prior does not remove it and the ceiling
  bounds the ratio at 34.4. M1 does not test it; the 145 C rows of the objective decide where it goes.
- **`Ea_decay_thiol_sink` stays FREE** inside Gigl's measured band (7 - 102 kJ/mol), as in B60. Reason: it is shared
  by the whole thiol-sink family (`k_mft_decay`, `k_fft_decay`, `k_thiol_decay`, `k_thiolate_loss`). Banding it needs a
  centre, and the only measured one (Gigl's 60.2, a coffee-brew NMR transfer) would cut lambda(145) / lambda(100) from
  33.5 to 8.1 for all four sinks at once: a structural change carried by a transfer, applied to both rivals, that the
  verdict would then partly measure. **If it binds on its 102 ceiling in either rival** (`active_mask`), it is
  profiled (sec. 7), D1 is labelled **"conditional on the Gigl ceiling"** (review M4), and a no-refit D1 sensitivity is
  printed with the ceiling relaxed (Ea at 112 and 122 kJ/mol, every other coordinate at the optimum, k at 145 C fixed).
- Not used: EX-050 / EX-055's Cys-conditional bound, for the reasons above.

## 3. Feasibility before any fit (no fit)

The owner's rule: test the structures on the mechanism signatures before any constant is fitted. v3 replaces v2's
two signatures and their thresholds (review M10) with **one rule on the decisive data**.

**Probe.** For C1, 256 Sobol draws over its bands, every other constant at B60: `k_ttca_cys` (+-3 dex),
`ea_form_kj_mol` (40 - 135), K (0.045 - 0.075), `k_form` (+-0.44 dex), `ea_decay_kj_mol` (40 - 135),
`Ea_decay_thiol_sink` (7 - 102), and lambda(100 C) of `k_mft_decay` from the prior centre -2 dex up to the fit's
reachable maximum at that draw's Ea (k(145 C) <= 10^0.5; ~0.094 /min at Ea 101.25; re-review MINOR 3).

Each draw is scored on the decisive units (sec. 5: D1 and the reservoir block) by the plug-in density of sec. 7
(predictions below an assay floor scored at the floor), isothermal reading, inline, fit tolerance. R0-shipped (b) is
scored the same way. From the best draw a local search on that same log-likelihood (least squares, <= 200 evaluations,
inside the same bands) is run, so a thin region is not missed.

**The one rule.** C1 is **killed** only if its best value (the best draw or the local search's end, whichever is
higher) is not higher than R0-shipped (b)'s on the decisive units. No other threshold.

**The guard against stacking.** The probe holds B60's 16 fitted coordinates at values fitted under R0's instant supply,
so the errors R0's constants absorb can make a correct structure look worse than R0 (CF-107's stacking). The probe
therefore cannot remove C1 on its own: if the rule would kill it, C1 is fitted once (start 0 only, sec. 4.3), and the
kill stands only if the same rule fails at that fit's optimum. Passing is screening, never evidence in the verdict.
v2's S1 (the Schieberle rise) and S2 (the Cao reservoir) values are printed at the best draw as descriptions, with no
pass mark.

Artifact: `results/validation/kinetic_core_m1_feasibility_probe.{md,json}`, generator
`scripts/generators/m1_combination_feasibility.py` (to write). Seed fixed and printed.

## 4. The fit rows and the fit

### 4.1 Objective

| block | rows | rivals | measurement model |
|---|---|---|---|
| **B60's 52 rows** (R0's own objective) | exactly B60's `row_ids` | R0, C1 | B60's: its sigmas, `weight_tag` half, pH endpoint weight E 0.7 decades per pH unit, sigma_pH 0.5 |
| **route block** | the verified rows of [route_submodel.yml](../../data/lit/fit_rows/route_submodel.yml): **296 rows in 34 pots** (Cao 2017 Figs. 4 - 5; Zhai 2019 Fig. 3A, 70 - 130 C; Zhai 2021 Figs. 3a (120 C), 4b, 6; Zhao 2017b; Zhai 2023 Food Chem and JAFC at 100 / 120 C). Out: the two Zhang 2025b pots (31 rows; the reservoir block and its record) and **the five Zhai 140 C pots** (30 rows: `zhai2021_ttca_ph55_140c_fig3a`, `zhai2023foodchem_{ttca,ttcaxyl}_ph7_140c_fig3a`, `zhai2023jafc_{ttca,ttcacys}_ph7_140c_fig5a`), a declared hold-out (FIT_HOLDOUT_DECLARATION: "every Zhai >= 140 C pot"). The table has no pH rows | C1 | per row sigma^2 = read_error^2 (or the printed SD in dex) + s_series^2, s_series per data set x species by `fitting.estimate_series_noise`: estimated at C1's start 0, held for its starts, re-estimated at its optimum (sec. 4.3) |
| **prior** | one term: (log10 lambda_mft(100 C) - log10 0.0406) / 1.0 | R0, C1 | sec. 2.6 |

**Why not R3En's 521 rows.** R3En has Schieberle 2000 (D1) and Wang 2026b / Zhou 2022 disulfides as fit rows, runs
B59a and B61 on, and lacks the fed-pot disulfide rows and Kang rows of R0's objective. B60's rows keep the comparison
with R0 like-for-like; the route block adds the rows that measure the route's own species, without which no fit row
measures a timescale (B58) and `k_ttca_cys` is not identified. R0 is not fitted on the route block: half its species
(TTCA, CAMA) are structural zeros for R0, and its pentose steps (`k_pent_*`) are frozen (sec. 5).

**Before the build, three checks on the rows:**
- *Duplicates.* `lab_ids.yml` lists `kang2026_SI` as a copy of `zhai2023foodchem`. B60's Kang rows
  (`kang_*_cys_conversion`, `kang_free_cys_yield_ceiling`) may be the same data as the route block's Zhai 2023 Food
  Chem 120 C pots. If they are, the route block drops those pots (B60's rows stay, so R0's objective remains a subset).
- *Declared hold-outs inside B60's objective (re-review MINOR 10).* `zhou_pH7_dimer_over_MFT` and
  `zhang_fig1_cys_dimer_over_MFT` measure the Zhou 2023 / Zhang 2024 dimer fractions that FIT_HOLDOUT_DECLARATION
  lists as HOLD-OUT ("holding both out keeps that agreement as an out-of-sample fact"), yet they are rows of B60's
  objective. M1 keeps them (removing them would change R0's objective, and both rivals carry them equally), registers
  the contradiction (sec. 11 item 22), and tags them "declared hold-out, in-sample since B60" in sec. 8.
- *Weight.* The route block has 296 rows against 52. The series noise keeps a misfitting series from dominating, but
  CF-107 found 77 % of B16's cost on TTCA rows. Each optimum reports the share of its cost from each block. A record,
  not a gate.

Not in any objective: Schieberle 2000 (D1); Zhang 2025b; every row of the side-effect record (sec. 8). Grosch 1992 is
in no objective and no decisive set.

### 4.2 Free coordinates

| | B60's 16 | `k_mft_decay` (prior) | `k_ttca_cys` | `ea_form` | K_c | total |
|---|---|---|---|---|---|---|
| R0 | yes | yes | -- (held at B60) | -- | -- | 17 |
| C1 | yes | yes | yes | yes | -- | 19 |
| C2 (specified, not run) | yes | yes | yes | yes | yes | 20 |

### 4.3 Optimiser, starts, convergence

- `src/kinetic_core/fitting.py`: `fit`, `LEAST_SQUARES_SETTINGS` (trf, x_scale jac, ftol 1e-6, xtol 1e-8), the
  module's absolute finite-difference steps (never a relative `diff_step`). Budget 1 200 evaluations per start. Cost =
  the module's 0.5 x sum of squared weighted residuals (nats).
- **Three starts per rival.** R0's starts and C1's starts 0 and 1 run in parallel, one R0 fit and one C1 fit at a
  time, each with `MAILLARD_WORKERS=4` (8 cores in all, the owner's cap of ~8 of 10). C1's start 2 waits for R0's
  carried optimum:

  | rival | start 0 | start 1 | start 2 |
  |---|---|---|---|
  | R0 | B60's optimum (ceiling and prior added; its inline cost is V6's reference) | B8's perturbation of start 0 (seed printed) | B60's start-1 endpoint (the other B60 optimum) |
  | C1 | B60's vector, new coordinates at their sec. 2 starts | B8's perturbation of start 0 | R0's carried optimum, new coordinates at their starts |

- **Run conditions (`tasks/lessons.md` 2026-11-04, "Never pause a fit that runs under a wall-clock limit"; on the h0
  branch, not yet on this one).** Every pot integration runs under a wall-clock
  alarm, so a paused or memory-starved integration enters the objective as a failed pot run. Before each launch, free
  memory is checked (`vm.swapusage`, `top -o mem`) and the two parallel fits are launched only if both fit in RAM
  (V4's peak memory x 2); otherwise they run one at a time. A running fit is never paused (no SIGSTOP); concurrency
  changes only at launch. **Before any fit record is used, its `failed_pot_runs` is read: a record with
  `failed_pot_runs` > 0 is re-run with the same settings** and the failed record is kept, labelled.
- **A start counts only with `least_squares` status > 0.** A start that ends at status 0 (budget) is continued once
  from its endpoint with another 1 200 evaluations; still status 0, it is "not converged" and not used.
- **The comparison rule (re-review M2; replaces `two_starts`' rule, sec. 11 item 18).** A rival is **compared** iff
  (a) at least two of its starts end with status > 0 and costs within **1.92 nats** of the lowest (Delta chi^2 3.84,
  one coordinate's 95 % profile interval; B60's two starts differ by 0.54 nats), and (b) their predictions on every
  decisive point (sec. 5) agree within max(sigma_decision of that point, the Laplace predictive sd at the better
  optimum), compared without reading the data y, under the isothermal reading at the scoring rtol (v3.3). If (a) holds and (b) fails, the verdict (sec. 7) is computed at each
  qualifying start after all fits have ended, and the rival is compared only if those verdicts agree. Otherwise it is
  **"not compared, flagged"**, with every start's cost, vector and decisive predictions listed.
- **Fallbacks, declared now.** If R0 is not compared, the verdict uses R0's lowest converged start, flagged; if no R0
  start converged, R0-shipped (b), flagged "different dimer bound and `k_mft_decay` treatment" (review B4). If C1 is
  not compared there is no verdict and no switch changes.
- **Series noise (review M11; re-review MINOR 5).** For C1, after the comparison rule, s_series is re-estimated at the
  carried optimum and the fit is continued once from it under the new noise (<= 300 evaluations). Status 0: continued
  once more (+300); still 0, the pre-re-estimation optimum is carried, flagged, with both noise sets recorded.
  Otherwise the continued optimum is C1's carried optimum. (R0's objective has no series noise.)
- **Laplace** at each carried optimum (the module's central differences and identification rule). Coordinates on a
  bound are profiled (sec. 7), never given a Laplace sd.
- **M1.0's G1 at each optimum (re-review MINOR 6):** oxidant use and the phi = 1 rise, re-run at R0's and C1's carried
  optima and printed beside M1.0's B60 values.
- R0 and C1 are not nested (the route retires two steps), so there is no nesting check.

### 4.4 Verification before any fit

- V1: the route key off: every existing test passes and B60's 52-row cost reproduces 27.132 on the node fixed point
  (quick setting), the frozen reference.
- V2: C1's switch on at B60's vector reproduces the Cao residuals of a reference **recomputed on main with M1's switch
  set** (B59a and B61 off; re-review MINOR 9), not the route sub-model probe's (which ran with them on).
- V4: one evaluation of each objective (R0's 52 rows; C1's 52 rows + route block) timed on inline, with its peak
  memory at 4 workers. Run after the plan's W1.3 profiling, itself on inline (sec. 10).
- V5: the dimer ceiling: one evaluation each at log10 `k_dimer_mft` = 4 and at the ceiling integrates without
  failure; its time is recorded; the panel's maximum MFT at B60 is printed (the 1 mM cap). A failed integration is a
  build error to fix before any fit; a slow one sets the stiffness allowance (sec. 12), and the ceiling stays.
- V6: R0 at B60's vector with the ceiling and the prior on inline: its cost is printed as the new reference (R0's
  start-0 cost); the prior term is 0 at B60's lambda.
- **V7 (inline vs the converged fixed point; re-review addendum B2).** Inline (scoring rtol 1e-8) against
  fixed_point_nodes at 81 nodes, rtol 1e-8, at B60's vector and at C1's start vector (route on), on: every decisive pot
  (Schieberle 100 C 12 h; Zhang 2025b 120 C) and the reservoir record's Cao Figs. 4 - 5 pots (none of them among
  B64's five probe pots); and B60's 52 rows (B64 item 5's missing full-objective diff), with the per-row differences
  printed. Tolerance: every decisive point within 0.1 sigma_decision (~0.004 dex per D1 point, 0.006 dex per ratio; 0.1
  x the row sigma for the reservoir points). A failure is a build error, investigated before any fit.
- **V8 (smoothness; re-review addendum B6).** On five coordinates (`k_mft_decay`, `Ea_decay_thiol_sink`, `k_dimer_mft`,
  `k_ttca_cys`, `ea_form`), the objective cost at steps of +-0.025 and +-0.05 dex (+-0.25 / +-0.5 kJ/mol for a barrier)
  about C1's start, on inline at fit tolerance and on the fixed point. Fails if any inline finite difference changes
  sign between the two step sizes where the fixed point's does not, or differs from the fixed point's by more than a
  factor 2. **If V8 fails, every fit runs on fixed_point_nodes at the quick setting** (B60's numerics), the scoring
  stays as declared, and sec. 12 reverts to v3.1's basis.

(V3, C2's factor check, is not run: C2 is not built in M1.)

## 5. The decisive data

Every row below was checked against [split.yml](../../docs/reference/contrast_benchmark/split.yml) (v1),
[split_v2_per_change.yml](../../docs/reference/contrast_benchmark/split_v2_per_change.yml) (E10 and E12, the errors M1
addresses; E6 / E7 / E11 / E18, whose laboratories M1 uses) and
[FIT_HOLDOUT_DECLARATION.md](../../docs/reference/FIT_HOLDOUT_DECLARATION.md). Both blocks are out of every
objective.

**A unit** is the independent piece the measurement model scores (sec. 6): one paper x species series (its points
share one level offset), or for D1 one species' pair of ratios (they share the 360 min point).

| block | data (record) | units | lab: v1; v2 | seen | other status |
|---|---|---|---|---|---|
| **D1, the thiol shape** | Schieberle 2000 Table IV (100 C, 0.5 M phosphate pH 5.0, ribose 100 + Cys 33 mM, laboratory autoclave): MFT and FFT, log ratios 360/30 and 720/360 (720/30 is their sum, reported only; the 60 min point is not used, as in v1) | 2 (MFT, FFT) | DFA: DEV; E10 not DEV, no TEST contrast; E12 DEV | at start (sec. 0) | calibration input, never a v2 TEST contrast. CF-055: the paper's text contradicts its own Table IV (x13 vs 7.9); M1 reads the table |
| **the reservoir** | Zhang 2025b Fig. 2g - h (120 C, 0.2 M phosphate pH 6.5, Cys 40 + Xyl 40 mM, 5 - 50 min, sealed vials, ice-water quench): free cysteine and free xylose, 8 points each | 2 (Cys, PENT) | BTBU: DEV in v1; v2: an E12 TEST laboratory (F28#25; these rows are not TEST); declared held-out lab for E6 / E11 (sec. 11 item 17) | at start (Cys +0.55, PENT +0.44 dex) | its Table 2 volatiles are a frozen hold-out bundle (Amendment 45): NOT read. CF-063: Cys and Xyl print one identical ELSD calibration curve |

**Moved out of the decisive set (re-review M1): Cao 2017 Fig. 4** (v3's D2-in). K = 0.060 is measured on those same
points (sec. 2.2) and Fig. 4 is in C1's objective, while R0 cannot fit it (its `k_pent_*` are frozen); the -2
"optimism correction" of v3 did not repair that and is dropped. Cao Fig. 4 joins the reservoir record.

**The reservoir record (reported, never summed):**
- R0's **structural zeros** (re-review B2 of v2): a structural zero is a species the rival's network cannot form from
  the pot's charge, decided from the network, not from the number. R0 has no step that forms TTCA or CAMA from cysteine
  and pentose, so its TTCA and CAMA in Cao 2017 Fig. 5, Zhai 2019 Fig. 3A (four pots) and Zhang 2025b Fig. 2a - b are
  structural zeros, reported as **"R0 cannot represent the measured reservoir (EX-060)"**, never summed (never -inf).
- C1's log densities on those rows (in-sample for Cao and Zhai; out-of-sample for Zhang 2025b, seen at start), and
  both rivals' on Cao Fig. 4 (in-sample for C1).

**Assay floors (re-review B1a).** A predicted level below the assay's floor is scored at the floor, not at the
integrator's value; the count of floored points per series is printed. Floors, per series, the lowest printed point
(no paper prints a detection limit for these): D1 MFT 4.5 ug and FFT 2.0 ug per 100 mL (Table IV's 30 min values;
the table's "< 0.1 ug" for other compounds is the alternative); the reservoir's Cys 0.80 and PENT 1.22 mmol/L (Fig. 2g
- h at 50 min; the shared ELSD curve's zero crossing, 0.22 mmol/L, is the alternative). The alternatives are printed as
sensitivities. R0's pentose (0.00e+00 at the integrator floor) is expected to be floored at most points.

Not decisive, recorded only: Yiltirak 2026 vs Schieberle (case D7: not separable on disk); Kumazawa 2003 (v1 TEST; a
constraint already in every objective); the pH controls (M1b's); everything in sec. 8.

### 5.1 What each rival predicts, from its mechanism

Directions are read from the structures and from the start-vector numbers of sec. 0, not computed for v3.

| block | R0 | C1 |
|---|---|---|
| D1 | MFT 360/30 short by 1.09 dex; falls 6 -> 12 h; FFT +0.36 (M1.0, at B60) | at its start, MFT +0.49 and FFT +1.58 dex: the reservoir over-feeds FFT. The fit does not see D1 |
| reservoir | no reservoir: free cysteine held up only by its own losses; pentose gone in seconds (M1.0), so floored | the reservoir holds cysteine and pentose for hours (start: Cys +0.55, PENT +0.44 dex) |

**What M1 cannot decide.**
- **Supply from D1.** D1 is out of every objective; v3.1's reason why it cannot separate C1 from C2 ("lambda absorbs
  C2's factor") was wrong, because nothing is fitted to D1 (re-review M3). The right statement: a loss law that free
  cysteine slows (C2, or the case's R2) is not run, and B60's in-sample cysteine contrasts might identify it; such a
  rival could move D1 as much as the route does. **A D1 gain is therefore "consistent with supply, not attributed".**
  Supply is credited on the reservoir block only.
- **The loss law.** No out-of-objective development datum discriminates C1 from C2 (Grosch 1992 Table VIII, the one
  candidate, is held-out E18 TEST data and was dropped unread in v3.1). C1 is carried as the simpler form, C2 is
  recorded as not excluded, and the question goes to W6 as a data need (sec. 14).

## 6. Tolerances, numerics of the scoring, thermal readings

**D1: one decision sigma (review B1; re-review M5).**
- Schieberle 2000 Table IV prints no replicate statement and no method statement for the table; the series is "data
  from reference 16" / "Hofmann and Schieberle, unpublished". **That it is SIDA is inferred** (the same laboratory's
  method for this pot; the dossier's reading), **not stated**. The same laboratory and pot (ribose 10 + cysteine 3.3
  mmol in 100 mL 0.5 M phosphate pH 5.0, laboratory autoclave) give one measured within-lab replicate statement:
  **Hofmann & Schieberle 1998 (10.1021/jf9705983), Table 1 footnote b, "mean values of triplicates and differed by not
  more than 10 %"**, applying to every table (hofmann1998_reconciliation.md sec. 2). The transfer to the 2000 series is
  declared; it is used as a noise floor. The paper's text is not reliable on its own numbers (CF-055: "by a factor of
  13" against Table IV's 7.9), so M1 reads only Table IV.
- **Per point: sigma_p = log10(1.10) = 0.041 dex**, taken as the sd of a printed value. This is conservative twice (10 %
  is a maximum difference among three, which bounds the sd near 6 %; and a mean of three divides that by sqrt 3) and
  offsets, without measuring it, the run-to-run variation a within-series statement omits.
- **Per ratio:** sqrt(2 x (sigma_p^2 + read^2)) = **0.059 - 0.061 dex** (read error 0.001 - 0.011 dex per point,
  thiol_submodel.yml). The two ratios of one species share the 360 min point: covariance -(sigma_p^2 + read_360^2),
  correlation about -0.5. Each species is one bivariate Gaussian unit.
- **Sensitivities, never used to choose:** (a) the read error alone (0.001 - 0.011 dex); (b) the s042 series noise,
  MFT 0.80 / FFT 1.16 dex (a misfit property, M1.0).
- Not used: Grosch 1992 Table II's SIDA SDs (DFA, n = 2 - 4, about 5 - 30 %): that table is held out.

| block | noise as the source states it | sigma used |
|---|---|---|
| D1 | none printed; Hofmann 1998's triplicates within 10 % (above) | 0.059 - 0.061 dex per ratio, bivariate per species |
| reservoir | Zhang 2025b: n = 3, SD; points read from a raster ln(C/C0) plot (read error 0.026 dex Cys, 0.035 dex PENT); identical Cys / Xyl calibration curves (CF-063) | per point: the row's SD (in dex, where printed) and read error in quadrature; per series: a level offset N(0, tau_ABS(all)^2) integrated out (below). CF-063 is flagged: the two species may share a calibration error, so the 2 units may be fewer in effect |

**Across laboratories (H0, review M11; re-review MINOR 2).** The reservoir's levels are absolute concentrations
against constants shared across laboratories, so each paper x species series carries a level offset
N(0, tau_ABS(all)^2), integrated out in closed form (H0's form). **tau_ABS(all) is read from EX-059 (H0); if it is not
on disk when the decisive data are scored, the scoring waits; it never runs at tau = 0.** **Declared transfer:** H0
estimates tau_ABS from thiol and volatile residuals (sida / calibrated groups); applying it to HPLC-ELSD reactant
levels assumes a reactant assay's between-laboratory level spread is no larger than an odorant assay's. Sensitivities
at 0.5x and 2x tau_ABS(all) are printed. D1's ratios sit inside one series and carry no tau: the scale cancels.

**Numerics of the scoring.** Inline at rtol 1e-8 (sec. 2.5), checked by V7.

**Thermal readings (re-review B3, BLOCKING; supersedes v3.2-draft's dead-time sensitivity, M5).** The B64 prereg's
review (item 3) records that later waves run on every reading, isothermal included, and CR-03 / CR-04 are registered,
open alternatives. Schieberle's pot is the Hofmann autoclave, run "by raising the temperature within 20 min from 20 to
145 C" (Hofmann 1995, 1997). Its 30 min point may get far less than nominal heat, which moves both rivals' 360/30
ratios by much more than 0.06 dex, with an unknown sign of the R0 - C1 difference. **Every fit runs on isothermal**
(the shipped panel setting). **The decisive units are scored, no refit, under every registered reading:**

| reading | source | the Schieberle pot (100 C) | the reservoir pot (Zhang 2025b, 120 C) |
|---|---|---|---|
| isothermal | CR-03 `isothermal` (every wave through B63) | 100 C from t = 0 | 120 C from t = 0 |
| ramp only, linear | CR-03 `ramp_only` + CR-04 `linear_no_cooling` | 20 -> 100 C at the Hofmann rate (125 K / 20 min = 6.25 K/min: 12.8 min), then 100 C; sampling times counted from the start of heating | isothermal (v3.3: no invented come-up; the stand-in ramp is a printed sensitivity only) |
| ramp + hold | CR-03 `ramp_plus_hold` | the same come-up, then the stated time at 100 C | isothermal (v3.3) |
| ramp only, Newtonian | CR-04 `newtonian_no_cooling`, at B64's two members tau 4 and 6.7 min | exponential approach from 20 C, times from the start of heating | isothermal (v3.3) |
| ramp only, linear with cooling | CR-04 `linear_with_cooling` | reading 2 plus a cool-down of the ramp's duration after each sampling time | not applicable: an ice-water quench is stated |

That is 6 readings on the Schieberle pot. **The reservoir pot (Zhang 2025b, sealed vials) is isothermal under every
reading (v3.3, final check):** its come-up is not stated, and B64's convention records such pots as isothermal, "no
invented time constants". The Hofmann-rate stand-in ramp (20 -> 120 C in 16 min) is printed as a sensitivity OUTSIDE
the agreement requirement, so an invented ramp cannot veto the only block that credits supply. A reading applies to both pots of the same
evaluation pairwise by name (isothermal with isothermal, ...); the reservoir pot takes isothermal where its column says
not applicable. **The verdict counts only if it is the same under every reading;** otherwise it is "not separated",
with the reading(s) that disagree named (sec. 7).

## 7. The comparison and the decision rule

**One measurement model.** Each decisive unit u has data y_u (log ratios or log levels, dex), covariance Sigma_u from
sec. 6 (including the integrated level offset where it applies), and for rival k the prediction f_uk(theta_hat) at the
carried optimum (floored where sec. 5 says).

**The predictive density (re-review B1).**
- **Primary: the linearised predictive.** lppd_uk = log N(y_u; f_uk(theta_hat), Sigma_u + J_uk Sigma_theta J_uk^T).
  J_uk is the Jacobian of the unit's predictions on the identified free coordinates (central differences, the module's
  steps; a floored prediction has zero rows); Sigma_theta is the Laplace covariance. A coordinate on a bound
  (`active_mask`) enters with the variance of its one-sided profile density (below), uncorrelated (declared
  approximation). Unidentified free coordinates are **held at the optimum**; held constants at B60.
- **Printed beside it:** the plug-in density log N(y_u; f_uk(theta_hat), Sigma_u); and the Monte Carlo check, S = 64
  draws from the same distribution (identified from the Laplace, reflected at bounds; bound-active from the profile;
  unidentified held), with its MC SE per unit (delta method). A unit whose MC and linearised values differ by more than
  2 MC SE is flagged "non-linear predictive" for the results review.
- **Sensitivity, never used to choose:** the MC version with unidentified coordinates drawn log-uniform across their
  bounds (a prior-predictive draw across numerical bounds, which favours looser rivals: re-review B1b).
- **Profiles (re-review MINOR 4).** A bound-active coordinate is fixed at three values inside its bound (0.25, 0.5 and
  1.0 dex for a log rate constant; 1, 3 and 10 kJ/mol for a barrier), the other free coordinates re-optimised from the
  optimum with <= 300 evaluations per point; status 0: continued once (+300); still 0, the point is used and flagged
  "not converged" (its cost is then an upper bound, so the profile density there a lower bound). The profile cost c(v)
  gives a one-sided density exp(-(c(v) - c_min)), interpolated linearly in log density and truncated at the bound.

**Blocks, differences and SE (re-review B2).** For each block b with n_b units: d_u = lppd_u(C1) - lppd_u(R0);
Delta L_b = sum_u d_u; SE_b = sqrt(n_b / (n_b - 1) x sum_u (d_u - mean d)^2). With n_b = 2, SE_b = |d_1 - d_2|, so
"Delta L_b > SE_b" holds iff **both** units favour C1, and "Delta L_b < -SE_b" iff both favour R0: the rule reads as a
sign pattern, and no single unit's magnitude can decide it. No unit is credited across blocks, and v3's -2 correction
is dropped (no decisive datum is in either objective).

**The decision rule.**
1. The verdict is computed only between compared rivals (sec. 4.3), on the two blocks of sec. 5, under each thermal
   reading of sec. 6.
2. **C1 favoured** iff Delta L_reservoir > SE_reservoir (both reservoir units favour C1) AND Delta L_D1 >= -SE_D1 (C1 is
   not worse on the thiol shape: not both D1 units favour R0).
3. **R0 favoured** iff Delta L_reservoir < -SE_reservoir (both reservoir units favour R0).
   **Tie band (v3.3, final check):** with two units per block the rule is a sign test, so a unit "favours" a rival
   only if the sign of its d_u is the same under the linearised and the plug-in densities and under the floor and
   tau x0.5 / x2 sensitivities (all printed). Otherwise the unit is a tie. A tie in the reservoir block gives "not
   separated"; a tie in the D1 block counts as "not worse".
4. **Otherwise "not separated / conflict":** CF-107 stays open, no switch changes, and the result is **not evidence
   against C1**.
5. The verdict counts only if it is the same under every thermal reading; otherwise "not separated", naming the
   readings. If `Ea_decay_thiol_sink` binds in either rival, D1 is labelled "conditional on the Gigl ceiling" (sec.
   2.6).
6. **Credit:** supply is credited on the reservoir block only. The verdict text states that the reservoir block is
   one pot from one laboratory whose Cys and Xyl print one ELSD calibration curve (CF-063): "both units" is close to
   one datum, so any supply credit rests on Zhang 2025b alone. A D1 gain of C1 is reported as "consistent with
   supply, not attributed" (C2 / R2 not run, sec. 5.1).
7. **What ships:** C1's structure (the route switch on) ships if C1 is favoured, C1 was compared, the side-effect
   record (sec. 8) is complete with a candidate compensating error registered for each loss, no review objection (sec.
   9) is left unanswered, and the coupling flip is shipped first or declared part of the ship (sec. 2.5). If R0 is
   favoured, no switch changes and CF-107 stays open. The loss law is not decided either way.

**Also printed, per unit, per block, per reading and per sigma scale (decision, read error, s042):** both rivals'
linearised, plug-in and MC lppd, so a rival that wins a block while losing a unit is visible; the reservoir record
(sec. 5); and the sensitivities, each computed once at the carried optimum with no refit and never used to choose: the
floors' alternatives; tau_ABS at 0.5x and 2x; the reservoir pot at its label pH held fixed instead of the B51 balance's
drift (the pH exposure, M1b's error); K at 0.045 and 0.075 and `ea_decay` at 40 and 135 kJ/mol (the route's
T-dependence) on the decisive units and on B60's 145 C rows; the Gigl-ceiling relaxation (sec. 2.6); the
unidentified-wide MC.

**What "ships" means here.** M1 admits one supply structure for the sulfur lane. Its constants stay provisional until
M2's refit, because the plan joins M1 and M2 in one refit (same rows; CF-107).

## 8. Side effects

Every form is scored on inline. Before (R0-shipped (b), B60 on inline) and after (R0 refit and C1, if compared), every
row:

- **The numerics row:** B60 on inline against B60 on the node fixed point, over every item below, so no before / after
  change is partly numerics.
- **The full held-out contrast benchmark**, with `scripts/generators/score_contrast_benchmark.py --score`. It scores
  the shipped engine only today; the build adds an option to score a named Config and vector (UNVERIFIED that no
  other change is needed). Coverage first, then every answerable item: the 65 v2 TEST atlas cells and the 15 keyed
  held-out tables (Grosch 1992 Table II among them, scored as TEST after the verdict like every keyed table, never
  used to choose). The frozen before-record stays
  [kinetic_core_contrast_benchmark_score_v1.md](kinetic_core_contrast_benchmark_score_v1.md) (B60, fixed point).
- **The in-sample rows:** B60's 52 rows and the route block, per row, at each optimum vs B60; **the route residuals by
  pH** (5.5, 6, 7, 7.4, 8; re-review MINOR 8: 58 % of route rows are at pH >= 7).
- **M1.0's G1** at each optimum (sec. 4.3).
- **E10 and E12 TEST contrasts,** read here and only here: F04#26 (Mizutani 1960; E10), F28#8 (Ehrenberg 1989), F28#20
  (Xu 2011), F28#25 (Wang 2022), F28#26 (Liu 2023b) (E12), where the scorer can answer them.
- **Declared hold-outs of every fit:** Yiltirak 2026 (100 - 130 C ladder); Hofmann 1998 T1 / T2 pH 3 / 7; every Zhai
  >= 140 C pot, including the five route-table pots removed in sec. 4.1 (30 rows); the Zhou 2023 pH 6 / 8 columns.
  Not unseen (sec. 0); recorded, never used to choose.
- **The disulfide record (v2's D4 - D6, moved here; the C3 follow-up's data):** the Zhou 2023 pH 7 share and the Zhang
  2024 Fig. 1 share ratio (both in B60's objective, both declared hold-outs, sec. 4.1); the Zhou 2023 pH 6 / 8 shares
  (held out, above); Zhou 2022 MFTD / MFT at four times, **not independent of the in-sample Zhou 2023 pH 7 pot it
  replicates** (review M3); Wang 2026b MFTD / MFT and FFTD / FFT (60 - 150 min); Whitfield 1999's Cys arm and NF + H2S
  arm, **each arm counted once** (review M3); Zhang 2024 Fig. 1's free-MFT ratio, cystine arm over cysteine arm (x1.5
  measured). Per rival, with each source's printed noise (printed SD, or half the drawn bar plus read error); never
  summed, never used to choose.

**Tags (what M1 spends; sec. 11):**
- de Roos 2005 Fig. 4 (keyed E10 TEST table): a constant source of C1: "spent" in C1's rows.
- **Lost laboratory independence** ("same laboratory as M1's fitting or choosing data"): BTBU (M1 fits Cao 2017 and
  Zhao 2017b and chooses on Zhang 2025b): E12's F28#25 (Wang 2022); M2's E6 / E11 contrasts F23#5 (Wang 2012), F23#7,
  F24#4, F25#23 (Zhao 2019b); E18's F27#5 (Wang 2012); E7's keyed `heldout_zhao2017b_tab1.yml` (re-review NOTE). DFA
  (M1 chooses on Schieberle 2000): E18's F27#15 - F27#17 (Grosch 1992; re-review NOTE).
- Reading (Whitfield 1999, declared held out for E6 / E11 in split v2): already in R0's objective (B60's three
  Whitfield MFT rows); M1 adds none: "pre-existing, R0's objective".
- `zhou_pH7_dimer_over_MFT`, `zhang_fig1_cys_dimer_over_MFT`: "declared hold-out, in-sample since B60".
- Every sensitivity of sec. 7: printed, never used.

**Candidate compensating errors, named before the run.** A row that gets worse is expected where R0 relied on one of
these. Each loss in the record names one, or says none is known; each named one gets a CF entry at the decision.

| expected loss | why | candidate compensating error in R0 |
|---|---|---|
| Kang fed-TTCA rows (`kang_100C_MFT`, `kang_120C_MFT`, FFT, `MFT_over_FFT`; B60 +0.6 to +2.5 sigma) flip low | the route retires `r_pent_dpo`, which carried 9.4 of ~10 mM of their flux (B59b Amendment 3) | the instant cysteine-catalysed enolisation (`r_pent_dpo` / `r_pent_tdp`) standing in for the route; `k_nf_mft` ~x21 high (Cerny pot D, E6 / E11: M2) |
| `zhai_MFT_fold_120_over_100`, `zhai_FFT_fold_120_over_100` | the supply's barrier moves to `ea_form` | one barrier for every route (E11, M2); the hexose transfer of 100.9 kJ/mol |
| the 145 C fed + H2S rows (`fed_ribose_h2s_*`, `fed_tdp_h2s_FFT`) | the prior lets `k_mft_decay` move; the fast 145 C sink is no longer frozen | the frozen fast sink (CF-107); the H2S headspace bias (CF-065); the isothermal reading of the Hofmann pots (CR-03) |
| free MFT in oxidant-bearing pots (Zhang 2024's cystine arm; the Zhou 2023 pots), both rivals | the dimer ceiling lets `k_dimer_mft` rise above the numerical bound that held it | the declared ambient OX (1.0 mM, no measurement) and the irreversible dimer (no release: the C3 follow-up's subject) |
| `cerny_ternary_thiamine_share`, `fed_thiamine_MFT`; benchmark F27#15 (Grosch thiamine / ribose / Cys) | cysteine is sheltered in TTCA, so less is free to compete | the thiamine route (`k_thi_hmp`, `k_hmp_mft`) fitted under instant supply |
| `zhou_pH7_FUR_arp_alone`, `hofmann_ribose_FUR_insitu` | the pentose's fate changes | `k_pent_caramel` and `k_fur_decay` fitted against an instant pentose sink |
| the Zhou 2023 final-pH rows | TTCA formation and CAMA titration move the charge balance | the B51 balance's acid allocation (CF-084, M1b) |
| Kumazawa FFT retention | none expected: no route species in an FFT-only buffer | -- (a loss here is a build error until shown otherwise) |
| Yiltirak 100 - 130 C levels | the sink's T-dependence moves | laboratory and work-up offsets (CF-068, rival R4) |
| any row, by up to the numerics row's size | inline vs fixed point | none: numerics, read off the numerics row |

## 9. The adversarial review

- **Who:** an independent agent (a fresh general-purpose subagent with no prior context, launched by the session
  that runs M1; never a fork of the author's context and never the author).
- **Reviews.** (a) **Design**, before any build or fit: done on v2 (findings B1 - B5, M1 - M12, MINOR, NOTE). (a')
  **Re-review** of v3.1: done (verdict REVISE; findings B1 - B3, M1 - M5, MINOR 1 - 10, NOTE); both answered in sec.
  16.1. The review files are filed by the session that ran them. (a'') **A check of v3.2** by a fresh agent on the
  answers to the re-review, before the build. (b) **Results**, before the decision: the fit reports, the probe, the
  decisive-data and side-effect tables, the build diff and the data. It does not see the author's verdict paragraph or
  sec. 13. Its task: the strongest case AGAINST shipping the favoured rival and against the rule's application (the
  sign pattern per block, the readings, the floors, the profiles, the linearisation). **Its brief names the 31 constants
  C1 holds at values fitted under instant supply** (sec. 2.5; among them `k_dpo_ddp`, `k_ddp_mft`, `k_dpo_nf`,
  `k_dpo_c2c3`, `Ea_lumped_formation`) and asks whether any verdict rests on them.
- No reviewer sees sec. 13. Every objection is answered in the record (sec. 15) or the change waits.

## 10. Order of operations (blinding)

0. **Prerequisites, before any fit is scheduled:** the plan's W1.3 profiling of one evaluation on inline (Jacobian, pH
   solve, integrator) and V4's timing, so sec. 12 is re-estimated from measurement; H0's tau_ABS(all) (EX-059) on disk
   before the probe and the decisive data are scored (both wait for it).
1. Design review (done) -> v3 -> v3.1; re-review (done) -> v3.2; the v3.2 check (sec. 9 (a'')); v4 if it changes
   anything.
2. Build the dimer ceiling (a fit bound), the scorer option, the thermal-reading scoring and the two generators, off by
   default; V1, V2, V4 - V8.
3. Feasibility probe (sec. 3), and the confirmation fit if the rule would kill C1.
4. Fits, R0 and C1 in parallel under sec. 4.3's run conditions (C1's start 2 after R0's optimum); `failed_pot_runs`
   read on every record; comparison rule; C1's series-noise continuation; Laplace; profiles; G1.
5. Only after both fits have ended: one generator run scores both rivals on the decisive data under every reading
   (`scripts/generators/score_m1_decisive_data.py`, to write; artifact
   `kinetic_core_m1_decisive_comparison.{md,json}`). No decisive number is computed on an intermediate member.
6. Side effects (sec. 8).
7. Results review; answers.
8. Decision and its basis in a new EX entry (next free id, checked against origin/main and open branches); CF entries
   for the compensating errors named in sec. 8; the C3 follow-up and the W6 loss-law data need registered in the
   plan's backlog.

## 11. Differences from the case and the plan

1. **D2 does not separate C1 from C2 on mechanism.** The plan says the cysteine and sugar time courses separate C1
   from C2. C2's term acts on the thiols, not on cysteine or pentose, so on D2 the two predict the same thing. Zhang
   2025b's cysteine and pentose are the supply verdict's reservoir block.
2. **BTBU is a v2 E12 TEST laboratory.** The case marks Cao 2017, Zhang 2024 and Zhang 2025b "DEV (BTBU)"; true in v1.
   In v2, E12's test labs include `btbu_beijing` (F28#25, Wang 2022). Choosing on BTBU rows does not read F28#25, but
   F28#25 is then not laboratory-independent of the choice: tagged (sec. 8). (M1.0's "no TEST contrast" statement is
   about rows and holds.)
3. **de Roos 2005 Fig. 4 is a keyed E10 TEST table** (`heldout_deroos2005_fig4.yml`, keyed 2026-11-01: "TEST data
   for engine error E10 ... never a fit input"). C1's Cys-Amadori decay table (B59b Amendment 4) is read from that
   figure, and M1 addresses E10. It is used as a declared constant, which spends it as a test of E10: tagged "spent".
   No other measured source of this decay is on disk (Davidek 2002 is a glycine-Amadori analogue in another buffer).
   E10's remaining held-out contrast is F04#26.
4. **Zhou 2023 pH 6 / 8 columns are a declared HOLD-OUT** in FIT_HOLDOUT_DECLARATION.md while both frozen splits mark
   Zhou 2023 DEV. Decided 2026-11-04 (v2): the stricter file governs; the columns are read only in the side-effect
   record. With C3 withdrawn no rival is judged on the share's pH slope, so nothing is lost.
5. **C3 is withdrawn** (sec. 2.4), and with it the case's R1 replacement, open questions 1 and 5 (the oxidant step,
   the release stoichiometry) and v2's items 5, 6 and 9. They pass to the registered follow-up.
6. **M1 is the supply verdict only; C2 is specified, not run (owner, v3.1).** The plan expects M1 to choose the loss
   law. No out-of-objective development datum discriminates C1 from C2; in-sample cysteine contrasts might (re-review
   M3), which is why a D1 gain is not attributed to supply. C1 is carried as the simpler form; the loss-law question
   is a W6 data need (sec. 14). M2's refit may revisit C2 if W6 data arrive.
7. **C2's form is EX-050's cysteine factor without the electrophile pool E** (the case quotes k_E E / (1 + Cys / K_c)).
   E has no declared identity or measurement, and EX-050 / EX-055 judged the E form flexible, not identifiable (P 0.1).
8. **`k_mft_decay` has a declared weak prior** (sec. 2.6), not v2's Grosch band and not EX-050's Cys-conditional bound.
9. **R0 is refitted** (sec. 2.1) so that both rivals have the same dimer ceiling, prior and numerics.
10. **Out of M1, as the plan says:** headspace O2 (B11; case part B; open questions 3 and 6 lapse; Ehrenberg 1989 is
    not used); the pH balance (case part D, moved to M1b; held at B60); B27's source. Rival R4 (laboratory offsets) is
    not a combination; it enters only through H0's tau.
11. **D1's noise:** the case uses the measurement model's per-series noise. M1.0 found the s042 noise is a misfit
    property, so M1 decides on Hofmann 1998's measured replicate noise and reports the s042 and read-error scales as
    sensitivities (sec. 6).
12. **The fit objective** is B60's 52 rows plus the route block (minus the Zhai 140 C pots), not R3En's 521 rows.
13. **"Ships" is provisional until M2** (sec. 7): the plan's hard dependency "M1 and M2 in one refit".
14. **The case's D1 row** says R0 "falls after ~1 h"; at B60 it peaks at 2.6 h (M1.0). Immaterial to the design.
15. **The dimer's numerical bound is replaced by a physical ceiling** (sec. 2.5), a fit-side change in both rivals.
16. **Grosch 1992 is not used (v3.1).** Table VIII's reflux values are Table II nos. 1 and 4, E18's TEST contrast
    F27#15 and the keyed `heldout_grosch1992_tab2.yml` (CF-089). Grosch Table II / VIII stays unread by M1's choice and
    fit; the side-effect benchmark scores it as TEST after the verdict. Its laboratory (DFA) is the choosing
    laboratory of D1, so F27#15 - F27#17 lose laboratory independence (item 23).
17. **BTBU and Reading are M2's declared held-out laboratories.** Split v2 lists `btbu_beijing`, `reading_mottram_parker`
    and the study `whitfield1999` as `prereg_declared_heldout_not_dev` for E6 and E11 (M2's errors). M1 fits BTBU rows
    (the route block's Cao 2017 and Zhao 2017b) and chooses on BTBU rows (the reservoir block), and M1 and M2 are one
    refit, so M2's BTBU TEST contrasts (F23#5 Wang 2012, F23#7, F24#4, F25#23 Zhao 2019b) lose laboratory
    independence: registered here, tagged in sec. 8. No row of those contrasts is read. The alternative, dropping BTBU,
    is not taken: without Cao, Zhao and Zhang 2025b the route has no timescale (B58) and there is no reservoir block.
    Whitfield 1999 is already in R0's objective (B60); M1 adds no Whitfield row.
18. **The comparison rule departs from `fitting.two_starts`** (re-review M2). `two_starts` asks costs within 1 % of the
    cost and every identified constant within its Laplace sigma; a 1 % rule has a scale that depends on the rival's
    cost, and B60's own starts (0.54 nats apart) fail it. M1 uses 1.92 nats plus agreement of the decisive predictions
    (sec. 4.3). `two_starts` is not changed.
19. **Every decisive unit is scored under every registered thermal reading** (re-review B3; CR-03, CR-04, the B64
    prereg's review item 3), and the verdict counts only if it agrees across them (sec. 6). Fits stay isothermal, the
    shipped setting, until B64 decides CR-03.
20. **The numerics change to `ph.coupling = "inline"`** for every M1 computation (sec. 2.5), a kind-1 change checked by
    V7 / V8, with a fixed-point fallback for the fits if V8 fails, and a declared shipping rule.
21. **Cao 2017 Fig. 4 is not decisive** (re-review M1): K is measured on it and it is in C1's objective. The supply
    verdict's reservoir is Zhang 2025b only.
22. **Two declared hold-outs are rows of B60's objective** (re-review MINOR 10): `zhou_pH7_dimer_over_MFT` and
    `zhang_fig1_cys_dimer_over_MFT` (FIT_HOLDOUT_DECLARATION, the Zhou 2023 / Zhang 2024 dimer fractions). Kept in both
    rivals' objectives; registered and tagged.
23. **Lost laboratory independence beyond E6 / E11 / E12** (re-review NOTE): E18's F27#5 (Wang 2012, BTBU) and
    F27#15 - F27#17 (Grosch 1992, DFA), and E7's keyed `heldout_zhao2017b_tab1.yml` (BTBU). Tagged in sec. 8.

## 12. Compute estimate

**Measured basis.** B60, the same 52 rows, node fixed point, quick setting: 7601.5 s / 246 evaluations = 30.9 s per
evaluation (start 1: 28.9 s; `kinetic_core_b60_fit_report.json`). The inline coupling is recorded at ~3.6x faster per
evaluation (data_restructure_plan C5: one objective evaluation 27 s -> 7.5 s), so **R0's objective is estimated at
7.5 - 10 s per evaluation on inline** (30.9 / 3.6 = 8.6 s). **C1's objective adds 34 route pots: 10 - 15 s**
(UNVERIFIED). The W1.3 profiling and V4 confirm both before any fit is scheduled; if V8 fails and the fits revert to
the fixed point, v3.1's estimate applies (~25 - 50 h typical). Machine: 16 GB RAM, now available in full (the owner
closed the VM that held 10 GB); every job `nice`d; R0's and C1's fits run side by side, each with
`MAILLARD_WORKERS=4` (8 of 10 cores), under the run conditions of sec. 4.3. The per-evaluation times above assume
that two parallel fits do not slow each other (UNVERIFIED; V4 times one evaluation with the other fit running).

| step | count | estimate (unattended unless stated) |
|---|---|---|
| build + tests (author work) | the dimer ceiling, a scorer option, the thermal-reading scoring, two generators | 5 - 8 h of work |
| W1.3 on inline + V1, V2, V4 - V8 | one evaluation profiled; B60 cost; anchors; V7 on 3 pots + 52 rows at 81 nodes rtol 1e-8 at two vectors; V8 on 5 coordinates x 4 steps x 2 numerics | ~2 h of work + ~1 - 1.5 h |
| feasibility probe | 256 draws x 2 pots, inline, + a local search | 0.1 - 0.3 h; + one start (~1 - 1.5 h) if the rule would kill C1 |
| R0 fits | 3 starts x 250 - 300 evaluations x 7.5 - 10 s | 1.6 - 2.5 h |
| C1 fits | 3 starts x 250 - 330 evaluations x 10 - 15 s | 2.1 - 4.1 h |
| budget continuations | ~40 % of starts, ~half a start each | 1 - 2 h |
| C1's series-noise continuation | <= 300 (+300) evaluations | 0.8 - 2.5 h |
| Laplace | 2 x (2n + 1) evaluations | ~0.2 h |
| profiles | per coordinate 3 points x <= 600 evaluations x 7.5 - 15 s = 1.9 - 7.5 h; expected 0 - 2 per rival | typical ~4 h; worst ~20 h |
| decisive data | linearised (Jacobians on the decisive pots), plug-in and MC (64 draws) for 3 forms (R0, R0-shipped (b), C1) x 2 decisive pots + Cao, under 6 / 5 thermal readings, + the sensitivities | 1.5 - 3 h |
| G1 at each optimum | 2 runs | ~0.1 h |
| side effects | the numerics row; the contrast benchmark x 3 inline configurations; in-sample rows; route by pH | 1 - 3 h (UNVERIFIED: the scorer's run time) |
| stiffness allowance (V5) | 0 - 50 % on the fit lines if V5 finds slow integration at high `k_dimer_mft` | 0 - 4 h |
| reviews | v3.2 check + results review + answers | 4 - 5 h of work |
| **total (machine-hours)** | | **~11 - 15 h of work; ~12 - 25 h of compute typical (~18 h); worst ~75 h (every start at its cap and continuation, four profiled coordinates), ~100 h with the full stiffness allowance** |
| **wall clock, R0 and C1 in parallel** | the fit block (R0 fits ‖ C1's starts 0 - 1; then C1's start 2, its continuations and the noise continuation; the profiles of the two rivals in parallel) | **~9 - 18 h typical (~13 h); worst ~50 h (~65 h with the stiffness allowance); about 1.5 - 2 calendar days** |

The parallel saving is in the fit and profile block (R0's ~2 h runs under C1's ~3 - 4 h), so the wall clock is the C1
path plus the serial steps (probe, V-steps, scoring, side effects). v3.1 (fixed point, sequential, without the readings
and the profiles' new budget): ~25 - 50 h typical, ~110 h worst. Memory: one evaluation with 4 workers is measured in
V4; two fits run in parallel only if twice that fits in free memory at launch; otherwise sequentially (the
machine-hours row then is the wall clock).

## 13. Expectations (written now; not blind, sec. 0)

| event | P |
|---|---|
| V7 passes (inline within 0.1 sigma_decision on every decisive point) | 0.7 |
| V8 passes (inline cost smooth; fits stay on inline) | 0.5 |
| C1 passes the feasibility rule without a confirmation fit | 0.65 |
| R0 / C1 compared (sec. 4.3) | 0.8 / 0.7 |
| verdict: C1 favoured / R0 favoured / not separated or conflict (incl. disagreement across readings) / no verdict (C1 not compared) | 0.22 / 0.08 / 0.40 / 0.30 |
| the isothermal-reading verdict differs from at least one other reading's | 0.4 |
| isothermal: both reservoir units favour C1 | 0.6 |
| isothermal: C1 not worse on D1 (not both D1 units favour R0) | 0.5 |
| at least one unit flagged "non-linear predictive" | 0.3 |
| lambda_mft(100 C) at C1's optimum within 0.5 dex of 0.0406 /min | 0.6 |
| `Ea_decay_thiol_sink` on its 102 ceiling in at least one rival (D1 then "conditional on the Gigl ceiling") | 0.5 |
| `k_dimer_mft` above log10 0.5 in both refit rivals | 0.85 |
| C1's 52-row cost below R0's refit cost | 0.3 |
| at least one held-out benchmark row gets worse by more than its noise (favoured rival vs R0-shipped (b)) | 0.85 |
| Zhang 2025b TTCA within stage 1's sigma (0.20 dex) at C1's optimum / CAMA (reservoir record) | 0.60 / 0.30 |
| a W6 dataset that decides C1 vs C2 exists by M2's refit | 0.10 |
| M1 reaches a decision by 11-20 (plan sec. 8) | 0.45 |

## 14. What this does not do

**The loss law, registered as a data need for plan W6.** The question: C1 vs C2, does free cysteine slow the thiol
loss at process temperature? The measurement that would decide it (from the case's open question 2 and the reviews'
M4 / M5 / M3): **a fed-MFT (and FFT) loss rate at process temperature (100 - 145 C, pH 5 - 6), with and without free
cysteine (0 and at least two levels between ~1 and ~50 mM), in one study, from a development laboratory**, with the
thiol quantified by SIDA (labelled standard added after heating, so the in-pot loss is not compensated) over a time
course, free cysteine measured in the same pots, and the vessel's heating curve logged (so CR-03 / CR-04 do not
apply). Best in a reacted ribose + cysteine background, so the pot's electrophiles are present; a buffer-only arm
separates the sink from the matrix. C2 predicts the loss falling as 1 / (1 + [Cys] / K_c); C1 predicts no cysteine
dependence. Nothing on disk meets it: Grosch 1992 is held out and its Table III confounds formation with loss; van
Seeventer 2001 is 50 C storage; the Nestec patent is not measurement-grade. Before W6, the in-sample cysteine contrasts
of B60's rows (sec. 2.3) can be checked for whether they identify K_c at all.

**The constants M1 does not refit (review M4).** C1 holds 31 constants at values fitted under R0's instant supply,
among them `k_dpo_ddp`, `k_ddp_mft`, `k_dpo_nf`, `k_dpo_c2c3` and `Ea_lumped_formation` (sec. 2.5). Any verdict is
conditional on them; M2's joint refit is where they move.

Also not done: no disulfide side (the C3 follow-up, sec. 2.4); no oxidant from headspace O2; no pH-balance change; no
MFT route shares (M2); no posterior (W4); no thermal-reading decision (B64 owns CR-03 / CR-04; M1 only scores across
them). It does not test the Yiltirak vs Schieberle level gap (not separable on disk). It does not settle the Cao vs Zhai
Amadori-formation tension (phosphate or laboratory; stage 1); `ea_form` absorbs part of it, which the results review is
asked to weigh. It does not test the fast 145 C sink (sec. 2.6).

## 15. Results

*(appended after the run; nothing above this line is edited after the first M1 number exists)*

## 16. Revision log

- v1, 2026-11-04: first version.
- v2, 2026-11-04 (before the design review and any build): sec. 11 item 4 decided; the Zhou 2023 pH 6 / 8 columns
  stay held-out and D3 leaves the comparison.
- v3, 2026-11-04 (after the design review, before any build): C3 withdrawn to a registered follow-up; two verdicts on
  pre-declared subsets, R0's structural zeros reported apart; one decision sigma for D1 (Hofmann 1998); the dimer's
  numerical bound replaced by a physical ceiling in every rival, and R0 refitted under it; `k_mft_decay` on a weak
  prior, not on Grosch; Grosch (E18) and BTBU (E6 / E11) spending registered; the Zhai 140 C pots out of the
  objective; the one-SE decision rule with posterior-predictive SE and profiles; the comparison rule; one feasibility
  rule with a stacking guard; tau from H0 or wait; series noise re-estimated; compute from B60's measurement. Owner's
  decisions of 2026-11-04 (in the session that commissioned v3): items 1 - 12 of that brief.
- v3.1 (2026-11-04, before the re-review and any build): D2b dropped; C2 specified, not run; M1 = the supply verdict.
  Consequences: Grosch 1992 used for nothing (no E18 spend); the loss-law question registered as a W6 data need (sec.
  14); C2's key not built; compute and expectations revised (secs. 12, 13). Owner's decisions relayed by the
  coordinator, 2026-11-04.
- v3.2 (2026-11-04, after the re-review, before any build): assay floors instead of the integrator floor; unidentified
  coordinates held; the linearised predictive as primary, plug-in and MC printed; sec. 2.6 corrected; two blocks
  (D1, the reservoir = Zhang 2025b) and three outcomes, supply credited on the reservoir only, the -2 dropped; Cao Fig.
  4 out of the decisive set; sec. 5.1's reason corrected; the comparison rule at 1.92 nats plus prediction agreement;
  the Gigl-ceiling label and sensitivity; the SIDA wording and CF-055; every decisive unit scored under every
  registered thermal reading (re-review B3, superseding the dead-time sensitivity); `ph.coupling = "inline"` with V7 /
  V8, the R0-shipped forms (a) / (b), the numerics side-effect row and the shipping rule; MINOR 1 - 10; the NOTE's
  prior centre and independence tags; compute on inline, with R0's and C1's fits in parallel under the 2026-11-04 run conditions (failed_pot_runs
  read before any record is used) (secs. 4.3, 12, 13). Owner's decisions relayed by the coordinator, 2026-11-04 (items
  1 - 8, the addendum A, B, and the parallel-run item).

### 16.1 The reviews' findings and the answers

**The design review (on v2), as amended by v3.1 / v3.2.**

| finding | what the review found | change (section) or why not |
|---|---|---|
| B1 | sec. 7 summed D1 without saying which sigma; MFT and FFT point opposite ways at C1's start | one decision sigma from measured replicate noise (Hofmann 1998, 0.059 - 0.061 dex per ratio, bivariate per species); read error and s042 printed as sensitivities (6); the start values recorded as seen (0) |
| B2 | non-discriminating data (D2p's 31 points) summed into L would decide by compensation; R0's TTCA / CAMA terms are log(0) | blocks scored as sign patterns (7); the reservoir block reduced to Zhang 2025b's free Cys and pentose; structural zeros defined from the network and reported as "R0 cannot represent the measured reservoir (EX-060)", never summed (5) |
| B3 | C3 built to fail: the trapped release arm, a net-sink branch, a mechanism-contradicting law, false nesting | C3 withdrawn; the follow-up's four requirements registered (2.4) |
| B4 | `k_dimer_mft` sits on a numerical bound that holds R0 / C1 / C2 below the shares while C3's `k_ox` floats | the bound replaced in both rivals by the encounter-limit ceiling, log10 9.47 at 145 C (2.5); R0 refitted under it (2.1, 4.3); V5 |
| B5 | D3 remnants (sec. 5 D4 row, sec. 7, sec. 13, the title "v1"); the five Zhai 140 C pots (30 rows) in the objective though declared held out | every D3 / D4 - D6 reference rewritten (5 - 8, 13); title updated; the five pots named and removed (4.1) and listed as hold-outs (8) |
| M1 | Grosch Table VIII = Table II nos. 1 and 4 = E18's F27#15 and a keyed held-out table; BTBU declared held out for E11; M1 and M2 are one refit | Grosch used for nothing (2.6; v3.1 drops D2b); the lost laboratory independence registered (11 items 17, 23) and tagged (8) |
| M2 | the `k_mft_decay` band is an inference from another matrix, with a chosen sd and C2's centre, from D2b's table, and with Ea near 102 it keeps the fast 145 C sink | weak prior centred on the incumbent, physical basis and identifiability cost stated; the fast 145 C sink stated plainly; `Ea_decay_thiol_sink` freed inside Gigl's band, with the reason (2.6) |
| M3 | Zhou 2022 replicates an in-sample pot; Whitfield counted twice | moot for the verdict (those data left the decisive set with C3); in the side-effect record Zhou 2022 is marked not independent and each Whitfield arm is counted once (8) |
| M4 | D2b spans [-0.22, 0] for C2, so it can only support C2; sigma 0.10 - 0.13, not 0.07; computation unspecified | v3 defined it; v3.1 drops it: it cannot move the loss verdict and is held-out E18 data (5.1, 11 item 16) |
| M5 | free Cys is near-flat in Schieberle, so D1 cannot separate C1 from C2 | superseded by the re-review's M3 (below): the loss law is not run, and a D1 gain is not attributed to supply |
| M6 | C3's pool is model-internal (TDP 1.6 -> 40 mM) | moot (C3 withdrawn); the follow-up must put the measured DPO / TDP in its objective (2.4) |
| M7 | D6 is set by the declared ambient OX 1.0 mM | moot (C3 withdrawn); the follow-up drops D6 or carries OX as a sensitivity (2.4) |
| M8 | the tie rule is intransitive; one sd chosen; the Laplace sd is invalid at active bounds | replaced in v3.2 by the block rule with three outcomes; bound-active coordinates profiled (7) |
| M9 | the convergence rule is ambiguous | replaced in v3.2 by the 1.92-nat rule with prediction agreement (4.3) |
| M10 | the feasibility screen has hidden thresholds and an undefined kill rule; S2 is circular | one kill rule (not better than R0-shipped (b) on the decisive units); a confirmation fit guards against stacking; S1 / S2 printed only (3) |
| M11 | tau = 0 while EX-059 is pending; s_series estimated at the start | tau_ABS(all) from H0 or the scoring waits; s_series estimated at C1's start and re-estimated at its optimum with a continuation (4.3, 6) |
| M12 | C1's start values on D1 and D2p were seen | marked "seen at start" (0, 5); the claim weakened to out-of-objective, not blind (5, 13) |
| MINOR: compute | ~40 - 60 h typical, not 20 - 25 (B60 ran ~31 s / evaluation) | measured basis stated; on inline 12 - 25 machine-hours typical (~18 h), worst ~75 - 100 h; ~9 - 18 h wall clock with R0 and C1 in parallel; W1.3 and V4 as prerequisites (10, 12) |
| MINOR: route T-dependence | the route's temperature law was not stated | the route's barriers and T-flat K declared, the 140 C extrapolation named, K and `ea_decay` band-end sensitivities (2.2, 7) |
| MINOR: pH exposure | which pH the rivals' pH-flat pieces meet was not stated | stated (2.2); a fixed-label-pH sensitivity (7); route residuals by pH (8) |
| MINOR: tie-break order | the order was asserted, not defined | no tie-break remains: the three-outcome rule needs none; R0 17 / C1 19 coordinates recorded (2.1) |
| MINOR items on C3 | -- | moot (C3 withdrawn) |
| NOTE | -- | not itemised in the brief v3 was written from; no change made for it |

**The re-review (on v3.1; verdict REVISE).**

| finding | what the re-review found | v3.2's change (section) or why not |
|---|---|---|
| B1a | R0's pentose is predicted at the integrator floor (0.00e+00; atol 1e-12 - 1e-14) against ~60 mM: a 14 - 16 dex miss set by atol | predictions below an assay floor scored at the floor (the lowest printed point per series, stated; alternatives as sensitivities); floored counts printed (5) |
| B1b | unidentified coordinates drawn log-uniform across numerical bounds: a prior-predictive draw favouring loose rivals | held at the optimum; the bounds-wide draw a sensitivity only (7) |
| B1c | at sigma 0.06 with both rivals > 5 sigma off, a vague prediction beats a sharp one by ~230 nats; MC noise 7 - 23 nats per unit at S = 64 | the linearised predictive is primary (no MC noise), plug-in and MC (with MC SE) printed, non-linear units flagged; the block rule reads signs, so magnitudes cannot decide (7) |
| B1d | sec. 2.6 claims the opposite | corrected (2.6) |
| B2 | "Delta L - 2 > SE" is a one-sample t on 6 non-exchangeable units; with one dominant unit SE = Delta L and C1 can never win; a conflict is labelled "R0 favoured" | two blocks, three outcomes, supply credited on the reservoir only, conflict reported as "not separated" and not as evidence against C1; with two units per block the rule is a sign pattern (7) |
| B3 (addendum, BLOCKING) | thermal history: later waves run on every reading; Schieberle's 30 min point may get far less heat; the R0 - C1 sign is unknown | fits isothermal; every decisive unit scored under every registered reading (CR-03, CR-04, B64's taus); the verdict counts only if it agrees across them (6, 7, 11 item 19) |
| M1 | D2-in used twice (K from the same points; in C1's objective); R0 cannot fit it; the -2 is not an optimism correction | Cao Fig. 4 moved to the reservoir record; the -2 dropped (5, 7, 11 item 21) |
| M2 | the 1 %-of-cost rule has a rival-dependent scale; B60's own starts (0.54 nats) would fail it | 1.92 nats plus agreement of decisive predictions without reading y; verdicts compared at each start if predictions disagree; R0 fallback its lowest converged start; departure declared (4.3, 11 item 18) |
| M3 | a D1 win cannot be credited to supply while C2 / R2 are not run; "lambda absorbs" is wrong (D1 is out of every objective); B60's rows hold cysteine contrasts that could pin K_c | D1 gain "consistent with supply, not attributed"; sec. 5.1's reason corrected; the in-sample contrasts named (2.3, 5.1, 7, 14) |
| M4 | with Ea capped at 102, D1 is tied to Gigl's ceiling; C1 holds 31 constants fitted under instant supply | "conditional on the Gigl ceiling" label and a no-refit relaxed-ceiling sensitivity (2.6, 7); the held constants named in secs. 2.5, 9 (the results-review brief) and 14 |
| M5 | the sigma transfer is fine as a floor, but SIDA is inferred for Table IV; CF-055; the 30 min point is exposed to come-up | wording fixed and CF-055 cited (6); the dead-time sensitivity is superseded by the thermal readings (B3) |
| MINOR 1 | sec. 5 calls BTBU "E12 neither"; E12's test labs include it (F28#25) | fixed (5, 11 item 2) |
| MINOR 2 | tau_ABS transferred from thiol residuals to HPLC reactant levels, undeclared | declared, with 0.5x / 2x sensitivities (6) |
| MINOR 3 | the probe's lambda band exceeds the fit's reachable range (cap ~0.094 /min at Ea 101.25) | the probe draws up to the reachable maximum at each draw's Ea (3) |
| MINOR 4 | profiles: 100 evaluations too few; no status-0 rule; worst case ~20 h | <= 300 per point, one +300 continuation, an unconverged point flagged as an upper bound on cost; worst case ~20 h carried (7, 12) |
| MINOR 5 | no status-0 rule for the series-noise continuation | one +300 continuation, else the pre-re-estimation optimum, flagged (4.3) |
| MINOR 6 | re-run M1.0's G1 at each optimum | added (4.3, 8) |
| MINOR 7 | the dimer ceiling is valid but loose; no compute allowance for stiffness | stated as loose (2.5); a stiffness allowance from V5 (12) |
| MINOR 8 | print route residuals by pH (58 % of route rows at pH >= 7) | added (2.2, 8) |
| MINOR 9 | V2's reference came from a probe with B59a / B61 on | the reference recomputed on main with M1's switch set (4.4) |
| MINOR 10 | `zhou_pH7_dimer_over_MFT` and `zhang_fig1_cys_dimer_over_MFT` are declared hold-outs in B60's objective | registered and tagged; kept in both objectives, with the reason (4.1, 8, 11 item 22) |
| NOTE | the prior centre is exactly log10 0.0406 = -1.391; tag the lost independence for E18 (F27#5, F27#15 - 17) and E7 (`heldout_zhao2017b_tab1`) | centre corrected (2.6, 4.1); tags added (8, 11 item 23) |
| owner question (numerics) | run on `ph.coupling = "inline"` | adopted for every computation, with V7 (agreement on the decisive pots and the 52 rows), V8 (smoothness, fixed-point fallback), R0-shipped (a) / (b), the numerics side-effect row and the shipping rule (2.5, 4.4, 8, 11 item 20, 12) |
- v3.3, 2026-11-04 (final check: FREEZE with two edits before the first decisive number): the reservoir pot is
  isothermal under every thermal reading (the stand-in ramp a sensitivity outside the agreement requirement); a tie band
  on the two-unit sign rule; the CF-063 single-datum note in the verdict text; sec. 13's no-verdict odds made
  consistent with P(C1 compared) 0.7; sec. 4.3's cross-start prediction check reads the isothermal reading at the
  scoring rtol; if V8 fails, V7's diff is printed at the carried optima.

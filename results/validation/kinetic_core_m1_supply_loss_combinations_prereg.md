# M1 pre-registration: thiol supply, judged as R0 against C1; the loss law specified, not run (v3.1)

*Plan step M1 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, row M1; sec. 7 items 6 - 9). Case:
[M1_thiol_loss_oxidant_ph_drift.md](../../docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md) (items 1 - 5; not
edited here). Conflict: CF-107 ([conflicts_register.yml](../../docs/reference/conflicts_register.yml)). M1.0: EX-060,
[kinetic_core_m1_preflight_checks.md](kinetic_core_m1_preflight_checks.md).*

- v1 and v2 written 2026-11-04, v3 and v3.1 the same day after the design review (sec. 16), all BEFORE any M1 build,
  fit or score. No engine was run with a changed switch or constant to write any version. Nothing above the "Results"
  line (sec. 15) is edited after the first M1 number exists; changes before that are versioned and logged in sec. 16.
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
  360/30 +1.58 dex (R0: -1.09 / +0.36), so D1's two species point opposite ways at the start.
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
- Register entries EX-024, EX-027, EX-050 (and its README), EX-055 (as quoted in EX-050), EX-060; CF-107.
- **Measured values only (no engine number):** Zhou 2023 Table 1 at pH 6 / 7 / 8 (MFT, FFT, MFTD, FUR, with SDs);
  Whitfield 1999 Table 1 row 48 (bis(2-methyl-3-furyl) disulfide 2 / 3 ug Cys arm, 7 / 2 ug H2S arm) and the dossier's
  7.5 % / 35 % shares; Grosch 1992 Tables II and VIII with their SDs (grosch1992 dossier; seen, used for nothing in
  M1 since v3.1); Hofmann 1998's error footnote (triplicates within 10 %); Guth 1995's untrapped arm (k_net <= 4e-4
  /min); the F28 atlas cells #1, #2, #9 - #12 and F27's Grosch cells; the de Roos 2005 keyed-table header (no value).
- **v3 additions:** the design review's findings (B1 - B5, M1 - M12, MINOR, NOTE); `route_submodel.yml`'s pot list
  (temperatures, species, row counts); `split_v2_per_change.yml`'s E6 / E11 / E12 / E18 lists.
- Not read: any engine number on the Zhou 2023 pH 6 / 8 disulfide shares, on Zhang 2024's free-MFT ratio, on
  Whitfield's disulfides, or on Grosch 1992 Table VIII.

## 1. Question

**Does the cysteine-pentose route's reservoir explain the thiol data the shipped structure cannot?** M1 gives one
verdict, **supply: R0 (the shipped structure) against C1 (the TTCA / cysteine-Amadori route)**, judged together with
the thiol sink because the errors stack (CF-107), on a pre-declared decisive subset (sec. 5).

- **The loss law (C1 vs C2: does free cysteine slow the thiol loss?)** cannot be decided on disk: no development
  source has a datum that discriminates the two (sec. 5.1). C2 is **specified, not run** (sec. 2.3), and the question
  is registered as a data need for plan W6 (sec. 14).
- **C3 (the disulfide side: Guth's release and an osone oxidant)** is withdrawn from M1 as a registered follow-up
  (sec. 2.4).

## 2. The rivals, exactly

### 2.1 Config keys

| rival | switches ON (Config key) | new step(s) | removed step(s) | fitted coordinates | in M1 |
|---|---|---|---|---|---|
| R0 | none (the shipped configuration) | -- | -- | 17 (sec. 4.2) | fitted, compared |
| C1 | `sulfur.cys_amadori_route` (`enabled = True`; it also turns on `r_ttca_form` by detailed balance with `sulfur.ttca_formation`'s K) | `r_ttca_form`, `r_cama_form`, `r_cama_dpo`, `r_cama_tdp` (exist, B58 / B59b) | `r_pent_dpo`, `r_pent_tdp` (retired by the route) | 19 | fitted, compared |
| C2 | C1 + `sulfur.cys_suppressed_thiol_loss` (a NEW key; not built in M1) | none: a factor on `ch_mft_decay` and `ch_fft_decay` | -- | 20 | **specified, not run** |

**R0 is scored in two forms.** *R0-shipped* is B60 as shipped (the incumbent, printed everywhere, used to choose
nothing except in the declared fallback of sec. 4.3). *R0* is the shipped structure refitted on its own objective
(B60's 52 rows) with the two fit-side changes C1 also gets: the dimer ceiling (sec. 2.5) and the `k_mft_decay` prior
(sec. 2.6). The verdict uses R0, so it compares structures under the same bounds and priors (review B4). Refitting R0
changes no switch; if R0 is favoured, its refit constants are recorded, not shipped (M2's refit owns the calibration).

### 2.2 Constants of the route (C1; C2 would share them)

| constant | step | status | value, band, source |
|---|---|---|---|
| K, TTCA equilibrium | `r_ttca_form` = K x `k_ttca_cys` | **declared, measured** | 0.060 L/mmol; band 0.045 - 0.075 (Cao 2017, `TtcaFormation`). Held at 0.060. pH- and T-flat (B58's declared assumption) |
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
the Cys-Amadori decay only; K and `k_ttca_cys` are pH-flat across the route block's pH 5.5 - 8 pots (Zhai 2021), while
the supply verdict's data sit at pH 5.0 - 6.5. Both are reported, not tuned (sec. 7, sensitivities).

### 2.3 C2, specified, not run

| constant | status | form, value, band, source |
|---|---|---|
| `K_c`, the cysteine level that halves the loss | **would be FITTED** | `ch_mft_decay` and `ch_fft_decay` run at k / (1 + [Cys]_free / K_c). [Cys]_free is the free-cysteine species, not TTCA- or CAMA-bound. log10 K_c (mmol/L) start 1.5, bounds 0 - 3 (1 - 1000 mM). No measurement pins it at process temperature: van Seeventer 2001 (+150 mM Cys) is a 50 C storage cell and cannot set a processing law (atlas rule) |

The factor would not act on `ch_thiolate_loss_*` (the Kumazawa-anchored step) or on the dimer; the `k_mft_decay` prior
(sec. 2.6) would then be a prior on lambda at zero free cysteine. Form: EX-050's cysteine factor without the
electrophile pool E (sec. 11 item 7). **Why it is not run (v3.1):** the loss verdict has no decisive unit on disk
(sec. 5.1), so a C2 fit would decide nothing in M1 and costs ~12 - 25 h of compute. The key is not built in M1. The
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

### 2.5 Held fixed, freed, and the dimer ceiling (R0 and C1)

- **Free, from B60's start:** the 16 coordinates B60's Laplace marks identified (B59b Amendment 2's rule):
  `k_tdp_fur`, `k_nf_mft`, `k_nf_mp3p`, `k_mgo_mp`, `k_fft_decay`, `k_fur_decay`, `k_pent_caramel`, `k_cys_thermal`,
  `k_thiolate_loss`, `ph_arp_secondary_ammonium_pKa`, `k_osone_decay`, `k_ttca_deg`, `Ea_decay_thiol_sink`,
  `k_dimer_mft`, `k_arp_dpo`, `k_arp_tdp`. B60's bounds, except the dimer's upper bound (below).
- **Held at B60's values:** the 7 unidentified at B60 (`k_glc_ha`, `k_dimer_fft`, `ph_acid_yield_per_sink_event`,
  `k_fur_fft`, `k_dimer_decay`, `k_thiol_decay`, `Ea_decay_carbonyl_sink`) and the 24 other frozen keys of B60's
  report, except `k_mft_decay` (now under the sec. 2.6 prior, both rivals) and, in C1, `k_ttca_cys` (fitted).
  `k_pent_dpo`, `k_pent_tdp` drop out with their steps in C1.
- **Switches at their shipped defaults:** `ph.coupling` fixed_point_nodes; `ph.node_grid` log_and_linear;
  `ph.fixed_point` relaxed_converged; `ph.b51` all on; `ph.trunk_solutes_in_balance` False and `ph.trunk_terms`
  False (the pH balance is M1b's); `sulfur.flux_budget` closed; `panel.thermal_history` isothermal; B59a
  (`sulfur.strecker_sulfur`), B61, B62, B63, B7a's split law, `sulfur.tdp_own_decay`, the pentose-amine entry and
  every other sub-model switch off. (R3En ran with B59a and B61 on; M1 does not.)
- **The fit's numerics:** B60's quick setting (17 pH nodes, rtol 1e-5). Scoring of decisive data and side effects:
  the careful setting (81 nodes).

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
- **The prior:** Gaussian on log10 lambda_mft(100 C), centre log10 0.041 /min (= -1.39, the shipped value, B60), sd
  1.0 dex. The centre is the incumbent's, so the prior moves no rival away from R0's value; the same term enters R0 and
  C1. **Physical basis:** the loss is covalent capture by the pot's electrophiles plus oxidation; no measurement sets
  its rate in these pots, so the prior states only the order of magnitude (2 sd: 4e-4 to 4 /min, a thiol pool turning
  over in seconds to days at 100 C). It is not blind: the author has seen Grosch's inferred 0.03 - 0.05 /min (sec. 0).
  The hard bounds stay `FITTED_SULFUR_BOUNDS_LOG10K`'s; their upper end (log10 0.5 at 145 C, a first-order step for
  which the module's numerical argument holds: t1/2 ~14 s against holds of >= 20 min) lies 0.36 dex above B60's value,
  and an optimum on it is profiled (sec. 7).
- **The identifiability cost.** At 1 dex the prior constrains nothing the data constrain. The sink can then trade
  against the supply in-sample (EX-024's pattern). M1 pays this in two places: lambda's Laplace sd is reported per
  rival, and the posterior-predictive scoring (sec. 7) carries lambda's width into each rival's predictive density on
  D1, so a rival whose D1 rests on a poorly fixed sink predicts D1 less sharply and scores lower. (It would also let
  C2's K_c trade against lambda in any pot whose free cysteine is near-constant, review M5: one more reason the loss
  law needs new data, sec. 14.)
- **The fast 145 C sink is kept, stated plainly.** With lambda(100 C) centred at 0.041 /min and
  `Ea_decay_thiol_sink` at 101.25 kJ/mol (0.75 below its 102 ceiling), lambda(145 C) = 33.5 x lambda(100 C) ~ 1.4 /min
  (t1/2 ~30 s). That is B60's fast 145 C sink, the anchor CF-107 names. The prior does not remove it and the ceiling
  bounds the ratio at 34.4. M1 does not test it; the 145 C rows of the objective decide where it goes.
- **`Ea_decay_thiol_sink` stays FREE** inside Gigl's measured band (7 - 102 kJ/mol), as in B60. Reason: it is shared
  by the whole thiol-sink family (`k_mft_decay`, `k_fft_decay`, `k_thiol_decay`, `k_thiolate_loss`). Banding it needs a
  centre, and the only measured one (Gigl's 60.2, a coffee-brew NMR transfer) would cut lambda(145) / lambda(100) from
  33.5 to 8.1 for all four sinks at once: a structural change carried by a transfer, applied to both rivals, that the
  supply verdict would then partly measure. If it ends on its 102 ceiling (`active_mask`), it is profiled (sec. 7), and
  the record states that the fast 145 C sink is then set by the ceiling and lambda(100 C), not by a measurement.
- Not used: EX-050 / EX-055's Cys-conditional bound, for the reasons above.

## 3. Feasibility before any fit (no fit)

The owner's rule: test the structures on the mechanism signatures before any constant is fitted. v3 replaces v2's
two signatures and their thresholds (review M10) with **one rule on the supply verdict's own data**.

**Probe.** For C1, 256 Sobol draws over its bands, every other constant at B60: `k_ttca_cys` (+-3 dex),
`ea_form_kj_mol` (40 - 135), K (0.045 - 0.075), `k_form` (+-0.44 dex), `ea_decay_kj_mol` (40 - 135), lambda(100 C) of
`k_mft_decay` (prior centre +-2 dex), `Ea_decay_thiol_sink` (7 - 102).

Each draw is scored on the supply subset (sec. 5: D1, D2-in, D2p; free cysteine and free pentose only) with the
decision measurement model of sec. 6 (point prediction, no posterior averaging), at the quick setting. R0-shipped is
scored the same way. From the best draw a local search on that same log-likelihood (least squares, <= 200 evaluations,
inside the same bands) is run, so a thin region is not missed.

**The one rule.** C1 is **killed** only if its best value (the best draw or the local search's end, whichever is
higher) is not higher than R0-shipped's log-likelihood on the supply subset. No other threshold.

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
| **route block** | the verified rows of [route_submodel.yml](../../data/lit/fit_rows/route_submodel.yml): **296 rows in 34 pots** (Cao 2017 Figs. 4 - 5; Zhai 2019 Fig. 3A, 70 - 130 C; Zhai 2021 Figs. 3a (120 C), 4b, 6; Zhao 2017b; Zhai 2023 Food Chem and JAFC at 100 / 120 C). Out: the two Zhang 2025b pots (31 rows; D2p) and **the five Zhai 140 C pots** (30 rows: `zhai2021_ttca_ph55_140c_fig3a`, `zhai2023foodchem_{ttca,ttcaxyl}_ph7_140c_fig3a`, `zhai2023jafc_{ttca,ttcacys}_ph7_140c_fig5a`), a declared hold-out (FIT_HOLDOUT_DECLARATION: "every Zhai >= 140 C pot"). The table has no pH rows | C1 | per row sigma^2 = read_error^2 (or the printed SD in dex) + s_series^2, s_series per data set x species by `fitting.estimate_series_noise`: estimated at C1's start 0, held for its starts, re-estimated at its optimum (sec. 4.3) |
| **prior** | one term: (log10 lambda_mft(100 C) - log10 0.041) / 1.0 | R0, C1 | sec. 2.6 |

**Why not R3En's 521 rows.** R3En has Schieberle 2000 (D1) and Wang 2026b / Zhou 2022 disulfides as fit rows, runs
B59a and B61 on, and lacks the fed-pot disulfide rows and Kang rows of R0's objective. B60's rows keep the comparison
with R0 like-for-like; the route block adds the rows that measure the route's own species, without which no fit row
measures a timescale (B58) and `k_ttca_cys` is not identified. R0 is not fitted on the route block: half its species
(TTCA, CAMA) are structural zeros for R0 (sec. 5).

**Before the build, two checks on the rows:**
- *Duplicates.* `lab_ids.yml` lists `kang2026_SI` as a copy of `zhai2023foodchem`. B60's Kang rows
  (`kang_*_cys_conversion`, `kang_free_cys_yield_ceiling`) may be the same data as the route block's Zhai 2023 Food
  Chem 120 C pots. If they are, the route block drops those pots (B60's rows stay, so R0's objective remains a subset).
- *Weight.* The route block has 296 rows against 52. The series noise keeps a misfitting series from dominating, but
  CF-107 found 77 % of B16's cost on TTCA rows. Each optimum reports the share of its cost from each block. A record,
  not a gate.

Not in any objective: Schieberle 2000 (D1); Zhang 2025b (D2p and its TTCA / CAMA); every row of the side-effect record
(sec. 8). Grosch 1992 is in no objective and no decisive set.

### 4.2 Free coordinates

| | B60's 16 | `k_mft_decay` (prior) | `k_ttca_cys` | `ea_form` | K_c | total |
|---|---|---|---|---|---|---|
| R0 | yes | yes | -- (held at B60) | -- | -- | 17 |
| C1 | yes | yes | yes | yes | -- | 19 |
| C2 (specified, not run) | yes | yes | yes | yes | yes | 20 |

### 4.3 Optimiser, starts, convergence

- `src/kinetic_core/fitting.py`: `fit`, `LEAST_SQUARES_SETTINGS` (trf, x_scale jac, ftol 1e-6, xtol 1e-8), the
  module's absolute finite-difference steps (never a relative `diff_step`). Budget 1 200 evaluations per start.
- **Three starts per rival**, run in the order R0, C1 (C1's start 2 needs R0's optimum):

  | rival | start 0 | start 1 | start 2 |
  |---|---|---|---|
  | R0 | B60's optimum (ceiling and prior added; V6) | B8's perturbation of start 0 (seed printed) | B60's start-1 endpoint (cost 27.671), the other B60 optimum |
  | C1 | B60's vector, new coordinates at their sec. 2 starts | B8's perturbation of start 0 | R0's carried optimum, new coordinates at their starts |

- **A start counts only with `least_squares` status > 0.** A start that ends at status 0 (budget) is continued once
  from its endpoint with another 1 200 evaluations; still status 0, it is "not converged" and not used.
- **The comparison rule (review M9).** A rival is **compared** only if its lowest-cost converged start agrees with at
  least one other converged start under `two_starts`' rule: costs within `COST_AGREEMENT` (1 %), and every identified
  constant within its Laplace sigma at the better optimum. A coordinate on a bound (`least_squares` `active_mask`)
  agrees when both starts put it on the same bound; unidentified coordinates are listed, not compared. Otherwise the
  rival is **"not compared, flagged"**, and the record lists every start's cost and vector. (B60's own two starts,
  27.132 and 27.671, differ by 2.0 % and would fail this rule; sec. 13 prices that.)
- **Fallback, declared now.** If R0 is not compared, the verdict uses R0-shipped, flagged "different dimer bound and
  `k_mft_decay` treatment" (the asymmetry of review B4). If C1 is not compared there is no verdict and no switch
  changes.
- **Series noise (review M11).** For C1, after the comparison rule, s_series is re-estimated at the carried optimum
  and the fit is continued once from it under the new noise (<= 300 evaluations). The continued optimum is C1's
  carried optimum; both noise sets and both costs are recorded. (R0's objective has no series noise.)
- **Laplace** at each carried optimum (the module's central differences and identification rule). Coordinates on a
  bound are profiled (sec. 7), never given a Laplace sd.
- R0 and C1 are not nested (the route retires two steps), so there is no nesting check.

### 4.4 Verification before any fit

- V1: the route key off: every existing test passes and B60's 52-row cost reproduces (27.132 at the quick setting).
- V2: C1's switch on at B60's vector reproduces the route sub-model probe's Cao residuals (the B59b anchor).
- V4: one evaluation of each objective (R0's 52 rows; C1's 52 rows + route block), timed, with its peak memory at 4
  workers. Run after the plan's W1.3 profiling (sec. 10).
- V5: the dimer ceiling: one evaluation each at log10 `k_dimer_mft` = 4 and at the ceiling integrates without
  failure; its time is recorded; the panel's maximum MFT at B60 is printed (the 1 mM cap). A failed integration is a
  build error to fix before any fit; a slow one is a compute cost, recorded, and the ceiling stays.
- V6: R0 at B60's vector with the ceiling and the prior reproduces 27.132 plus a prior term of 0 (the ceiling is
  inert at an interior point; B60's lambda is the prior's centre).

(V3, C2's factor check, is not run: C2 is not built in M1.)

## 5. The decisive data

Every row below was checked against [split.yml](../../docs/reference/contrast_benchmark/split.yml) (v1),
[split_v2_per_change.yml](../../docs/reference/contrast_benchmark/split_v2_per_change.yml) (E10 and E12, the errors M1
addresses; E6 / E11, M2's, because M1 and M2 are one refit) and
[FIT_HOLDOUT_DECLARATION.md](../../docs/reference/FIT_HOLDOUT_DECLARATION.md). "In" = in a rival's objective; "out" =
in no objective. **No datum enters Delta L that the review showed does not discriminate R0 from C1.**

**A unit** is the independent piece the measurement model scores (sec. 6): one paper x species series (its points
share one level offset), or for D1 one species' pair of ratios (they share the 360 min point).

| # | data (record) | in / out | units | lab: v1; v2 | seen | other status |
|---|---|---|---|---|---|---|
| D1 | Schieberle 2000 Table IV (100 C, 0.5 M phosphate pH 5.0, ribose 100 + Cys 33 mM): MFT and FFT, log ratios 360/30 and 720/360 (720/30 is their sum, reported only; the 60 min point is not used, as in v1) | out, both rivals | 2 (MFT, FFT) | DFA: DEV; E10 not DEV, no TEST contrast; E12 DEV | at start (sec. 0) | calibration input, never a v2 TEST contrast |
| D2-in | Cao 2017 Fig. 4 (90 C, 0.2 M phosphate pH 5.5, Cys 400 + Xyl 400 mM, 0.5 - 5 h): free cysteine and free xylose, 6 points each | in for C1 (route block); out for R0 | 2 (Cys, PENT) | BTBU: DEV; E10 DEV; E12 neither DEV nor TEST; E6 / E11 declared held-out lab (sec. 11 item 17) | Cao RMS at start | -- |
| D2p | Zhang 2025b Fig. 2g - h (120 C, 0.2 M phosphate pH 6.5, Cys 40 + Xyl 40 mM, 10 - 50 min): free cysteine and free xylose, 8 points each | out, both rivals | 2 (Cys, PENT) | BTBU: as D2-in | at start (Cys +0.55, PENT +0.44 dex) | its Table 2 volatiles are a frozen hold-out bundle (Amendment 45): NOT read |

**The decisive subset (supply, R0 vs C1): D1 + D2-in + D2p**, free cysteine and free pentose only. 6 units, 32 points.

**R0's structural zeros (review B2).** A structural zero is a species the rival's network cannot form from the pot's
charge, decided from the network, not from the number. R0 has no step that forms TTCA or CAMA from cysteine and
pentose, so its TTCA and CAMA in Cao 2017 Fig. 5, Zhai 2019 Fig. 3A (four pots) and Zhang 2025b Fig. 2 are structural
zeros. They are reported as **"R0 cannot represent the measured reservoir (EX-060)"**, never summed (never -inf). C1's
log densities on the same rows are printed beside them as the reservoir record (in-sample for Cao and Zhai;
out-of-sample for Zhang 2025b, seen at start). A representable species that a rival predicts below the integrator's
absolute tolerance (R0's pentose is expected there: gone in seconds, M1.0) is scored at that tolerance and flagged.

Not decisive, recorded only: Yiltirak 2026 vs Schieberle (case D7: not separable on disk); Kumazawa 2003 (v1 TEST; a
constraint already in every objective); the pH controls (M1b's); everything in sec. 8.

### 5.1 What each rival predicts, from its mechanism

Directions are read from the structures and from the start-vector numbers of sec. 0, not computed for v3.

| # | R0 | C1 |
|---|---|---|
| D1 | MFT 360/30 short by 1.09 dex; falls 6 -> 12 h; FFT +0.36 (M1.0, at B60) | at its start, MFT +0.49 and FFT +1.58 dex: the reservoir over-feeds FFT. The fit does not see D1 |
| D2-in / D2p | no reservoir: free cysteine held up only by its own losses; pentose gone in seconds (M1.0) | the reservoir holds cysteine and pentose for hours (start: Cys +0.55, PENT +0.44 dex on Zhang 2025b) |

**What the data cannot decide: the loss law.** C1 vs C2 has no decisive unit on disk. D1 cannot separate them,
because free cysteine is near-flat in Schieberle's pot, so C2's factor is a constant that lambda absorbs (review M5).
D2-in / D2p cannot, because C2's factor does not act on cysteine or pentose (sec. 11 item 1). The one candidate, Grosch
1992 Table VIII (v3's D2b), was one-sided (at most ~0.3 log-likelihood units for C2, one unit, no SE) and is held-out
E18 TEST data, so v3.1 drops it unread. The loss verdict is **"cannot decide on disk"**: C1 is carried as the simpler
form, C2 is recorded as not excluded, and the question goes to W6 as a data need (sec. 14).

## 6. Tolerances from each measurement's own noise

**D1: one decision sigma (review B1).**
- Schieberle 2000 Table IV prints no replicate statement; the series is "Hofmann and Schieberle, unpublished"
  (dossier). The same laboratory, the same pot (ribose 10 + cysteine 3.3 mmol in 100 mL 0.5 M phosphate pH 5.0,
  laboratory autoclave) and the same SIDA ([2H3]-MFT, [2H2]-FFT) give one measured within-lab replicate statement:
  **Hofmann & Schieberle 1998 (10.1021/jf9705983), Table 1 footnote b, "mean values of triplicates and differed by not
  more than 10 %"**, the paper's only error statement, applying to every table (hofmann1998_reconciliation.md sec. 2).
  The transfer from the 1998 tables to the 2000 time series is declared.
- **Per point: sigma_p = log10(1.10) = 0.041 dex**, taken as the sd of a printed value. This is conservative twice (10 %
  is a maximum difference among three, which bounds the sd near 6 %; and a mean of three divides that by sqrt 3) and
  offsets, without measuring it, the run-to-run variation a within-series statement omits.
- **Per ratio:** sqrt(2 x (sigma_p^2 + read^2)) = **0.059 - 0.061 dex** (read error 0.001 - 0.011 dex per point,
  thiol_submodel.yml). The two ratios of one species share the 360 min point: covariance -(sigma_p^2 + read_360^2),
  correlation about -0.5. Each species is one bivariate Gaussian unit.
- **Sensitivities, never used to choose:** (a) the read error alone (0.001 - 0.011 dex); (b) the s042 series noise,
  MFT 0.80 / FFT 1.16 dex (a misfit property, M1.0). The verdict is printed at (a) and (b) beside the decision sigma.
- Not used: Grosch 1992 Table II's SIDA SDs (DFA, n = 2 - 4, about 5 - 30 %): that table is held out.

| # | noise as the source states it | sigma used |
|---|---|---|
| D1 | none printed; Hofmann 1998's triplicates within 10 % (above) | 0.059 - 0.061 dex per ratio, bivariate per species |
| D2-in, D2p | Cao: printed SD (n = 3) where drawn, plus read error; Zhang 2025b: n = 3, SD; identical Cys / Xyl calibration curves (CF-063) | per point: the row's SD (in dex) and read error in quadrature; per series: a level offset N(0, tau_ABS(all)^2) integrated out (below). CF-063 is flagged: the two species of one pot may share a calibration error, so the 2 units per pot may be fewer in effect |

**Across laboratories (H0, review M11).** Terms whose scale does not cancel (D2-in and D2p levels: absolute
concentrations against constants shared across laboratories) carry a level offset per paper x species,
N(0, tau_ABS(all)^2), integrated out in closed form (H0's form). **tau_ABS(all) is read from EX-059 (H0). If it is not
on disk when the decisive data are scored, the run waits; it never runs at tau = 0.** D1's ratios sit inside one
series and carry no tau: the scale cancels.

## 7. The comparison and the decision rule

**One measurement model.** Each decisive unit u has data y_u (log ratios or log levels, dex), covariance Sigma_u from
sec. 6 (including the integrated level offset where it applies), and for rival k the prediction f_uk(theta).

**Posterior-predictive averaging.** For each rival, S = 64 draws theta_s:
- identified free coordinates from the Laplace covariance at the carried optimum (draws outside a bound reflected);
- a coordinate on a bound (`active_mask`) from its **profile**: the coordinate fixed at three values inside the bound
  (0.25, 0.5 and 1.0 dex for a log rate constant; 1, 3 and 10 kJ/mol for a barrier), the other free coordinates
  re-optimised from the optimum (<= 100 evaluations each, status recorded); the profile cost c(v) gives a one-sided
  density exp(-(c(v) - c_min)), interpolated linearly in log density and truncated at the bound; the other coordinates
  are then drawn from the Laplace with that coordinate removed (an approximation, declared);
- unidentified free coordinates across their bounds (log-uniform for rate constants), the envelope convention after
  B60's Laplace correction;
- held constants at B60 (they are part of the declared structure).

lppd_uk = log((1/S) sum_s N(y_u; f_uk(theta_s), Sigma_u)). Careful numerical setting throughout.

**The difference and its SE.** On the subset's n = 6 units: d_u = lppd_u(C1) - lppd_u(R0); Delta L = sum_u d_u;
**SE = sqrt(n / (n - 1) x sum_u (d_u - mean d)^2)**. D2-in is in C1's objective and not in R0's, so C1's Delta L is
reduced by 2 (one log-likelihood unit for each of C1's two coordinates beyond R0's, AIC's correction, a lower bound on
in-sample optimism). The verdict is also printed on the out-of-sample units alone (D1, D2p; n = 4); if the two
disagree, the record says the verdict rests on in-sample data and that goes to the results review as an objection.

**The decision rule (review M8; plan sec. 7 applied).**
1. The verdict is computed on the decisive subset (sec. 5), only between compared rivals (sec. 4.3).
2. **The favoured rival is the simplest one whose L is within one SE of the best.** Simplicity is the fitted-coordinate
   count: R0 (17) < C1 (19). So **C1 is favoured iff Delta L(C1 - R0) - 2 > SE; otherwise R0.**
3. A dominance flag (a record, not a threshold): if one unit carries more than half of |Delta L|, the record says so
   and the verdict goes to the results review as resting on one unit. Expected here, where R0's pentose sits at the
   integrator's floor (sec. 5).
4. **What ships:** C1's structure (the route switch on) ships if the verdict favours it, C1 was compared, the
   side-effect record (sec. 8) is complete with a candidate compensating error registered for each loss, and no review
   objection (sec. 9) is left unanswered. If R0 is favoured, no switch changes and CF-107 stays open. The loss law is
   not decided either way (sec. 5.1).

**Also printed, per unit and per sigma scale (decision, read error, s042):** both rivals' lppd, so a rival that wins
the sum while losing a unit is visible; the reservoir record (sec. 5); and the sensitivities, each computed once at
the carried optimum with no refit and never used to choose: D2-in / D2p at the label pH held fixed instead of the B51
balance's drift (the pH exposure, M1b's error); K at 0.045 and 0.075 and `ea_decay` at 40 and 135 kJ/mol (the route's
T-dependence) on the supply subset and on B60's 145 C rows.

**What "ships" means here.** M1 admits one supply structure for the sulfur lane. Its constants stay provisional until
M2's refit, because the plan joins M1 and M2 in one refit (same rows; CF-107).

## 8. Side effects

Before (R0-shipped = B60) and after (R0 refit and C1, if compared), every row:

- **The full held-out contrast benchmark**, with `scripts/generators/score_contrast_benchmark.py --score`. It scores
  the shipped engine only today; the build adds an option to score a named Config and vector (UNVERIFIED that no
  other change is needed). Coverage first, then every answerable item: the 65 v2 TEST atlas cells and the 15 keyed
  held-out tables (Grosch 1992 Table II among them, scored as TEST after the verdict like every keyed table, never
  used to choose). R0-shipped's record is the existing
  [kinetic_core_contrast_benchmark_score_v1.md](kinetic_core_contrast_benchmark_score_v1.md).
- **The in-sample rows:** B60's 52 rows and the route block, per row, at each optimum vs B60.
- **E10 and E12 TEST contrasts,** read here and only here: F04#26 (Mizutani 1960; E10), F28#8 (Ehrenberg 1989), F28#20
  (Xu 2011), F28#25 (Wang 2022), F28#26 (Liu 2023b) (E12), where the scorer can answer them.
- **Declared hold-outs of every fit:** Yiltirak 2026 (100 - 130 C ladder); Hofmann 1998 T1 / T2 pH 3 / 7; every Zhai
  >= 140 C pot, including the five route-table pots removed in sec. 4.1 (30 rows); the Zhou 2023 pH 6 / 8 columns.
  Not unseen (sec. 0); recorded, never used to choose.
- **The disulfide record (v2's D4 - D6, moved here; the C3 follow-up's data):** the Zhou 2023 pH 7 share and the Zhang
  2024 Fig. 1 share ratio (both in B60's objective); the Zhou 2023 pH 6 / 8 shares (held out, above); Zhou 2022
  MFTD / MFT at four times, **not independent of the in-sample Zhou 2023 pH 7 pot it replicates** (review M3); Wang
  2026b MFTD / MFT and FFTD / FFT (60 - 150 min); Whitfield 1999's Cys arm and NF + H2S arm, **each arm counted once**
  (review M3: v2 counted the Cys arm in two rows); Zhang 2024 Fig. 1's free-MFT ratio, cystine arm over cysteine arm
  (x1.5 measured). Per rival, with sec. 6's kind of sigmas (printed SD or half the drawn bar plus read error); never
  summed, never used to choose.

**Tags (what M1 spends; sec. 11):**
- de Roos 2005 Fig. 4 (keyed E10 TEST table): a constant source of C1: "spent" in C1's rows.
- BTBU contrasts: E12's F28#25 (Wang 2022) and M2's E6 / E11 contrasts F23#5 (Wang 2012), F23#7, F24#4, F25#23
  (Zhao 2019b): "same laboratory as M1's fitting and choosing data".
- Reading (Whitfield 1999, declared held out for E6 / E11 in split v2): already in R0's objective (B60's three
  Whitfield MFT rows); M1 adds none: "pre-existing, R0's objective".
- The pH-fixed and route-barrier sensitivities (sec. 7): printed, never used.

**Candidate compensating errors, named before the run.** A row that gets worse is expected where R0 relied on one of
these. Each loss in the record names one, or says none is known; each named one gets a CF entry at the decision.

| expected loss | why | candidate compensating error in R0 |
|---|---|---|
| Kang fed-TTCA rows (`kang_100C_MFT`, `kang_120C_MFT`, FFT, `MFT_over_FFT`; B60 +0.6 to +2.5 sigma) flip low | the route retires `r_pent_dpo`, which carried 9.4 of ~10 mM of their flux (B59b Amendment 3) | the instant cysteine-catalysed enolisation (`r_pent_dpo` / `r_pent_tdp`) standing in for the route; `k_nf_mft` ~x21 high (Cerny pot D, E6 / E11: M2) |
| `zhai_MFT_fold_120_over_100`, `zhai_FFT_fold_120_over_100` | the supply's barrier moves to `ea_form` | one barrier for every route (E11, M2); the hexose transfer of 100.9 kJ/mol |
| the 145 C fed + H2S rows (`fed_ribose_h2s_*`, `fed_tdp_h2s_FFT`) | the prior lets `k_mft_decay` move; the fast 145 C sink is no longer frozen | the frozen fast sink (CF-107); the H2S headspace bias (CF-065) |
| free MFT in oxidant-bearing pots (Zhang 2024's cystine arm; the Zhou 2023 pots), both rivals | the dimer ceiling lets `k_dimer_mft` rise above the numerical bound that held it | the declared ambient OX (1.0 mM, no measurement) and the irreversible dimer (no release: the C3 follow-up's subject) |
| `cerny_ternary_thiamine_share`, `fed_thiamine_MFT`; benchmark F27#15 (Grosch thiamine / ribose / Cys) | cysteine is sheltered in TTCA, so less is free to compete | the thiamine route (`k_thi_hmp`, `k_hmp_mft`) fitted under instant supply |
| `zhou_pH7_FUR_arp_alone`, `hofmann_ribose_FUR_insitu` | the pentose's fate changes | `k_pent_caramel` and `k_fur_decay` fitted against an instant pentose sink |
| the Zhou 2023 final-pH rows | TTCA formation and CAMA titration move the charge balance | the B51 balance's acid allocation (CF-084, M1b) |
| Kumazawa FFT retention | none expected: no route species in an FFT-only buffer | -- (a loss here is a build error until shown otherwise) |
| Yiltirak 100 - 130 C levels | the sink's T-dependence moves | laboratory and work-up offsets (CF-068, rival R4) |

## 9. The adversarial review

- **Who:** an independent agent (a fresh general-purpose subagent with no prior context, launched by the session
  that runs M1; never a fork of the author's context and never the author).
- **Three reviews.** (a) **Design**, before any build or fit: done on v2 (findings B1 - B5, M1 - M12, MINOR, NOTE;
  answered in sec. 16.1; the review file `kinetic_core_m1_supply_loss_combinations_prereg_review.md` is filed by the
  session that ran it). (a') **A re-review of v3.1**, by a fresh agent, on the changed sections (1, 2.3, 2.5, 2.6,
  3, 5 - 7, 14) and the answers table, before the build: whether each answer closes its finding and whether v3 / v3.1
  opened a new one. (b) **Results**, before the decision: the fit reports, the probe, the decisive-data and side-effect
  tables, the build diff and the data. It does not see the author's verdict paragraph or sec. 13. Its task: the
  strongest case AGAINST shipping the favoured rival, and against the rule's application (the dominance flag, the
  in-sample D2-in, the profiles).
- No reviewer sees sec. 13. Every objection is answered in the record (sec. 15) or the change waits.

## 10. Order of operations (blinding)

0. **Prerequisites, before any fit is scheduled:** the plan's W1.3 profiling of one evaluation (Jacobian, pH solve,
   integrator) and V4's timing, so sec. 12 is re-estimated from measurement; H0's tau_ABS(all) (EX-059) on disk before
   the probe and the decisive data are scored (both wait for it).
1. Design review (done) -> v3 -> v3.1; the re-review (sec. 9 (a')); v4 if it changes anything.
2. Build the dimer ceiling (a fit bound), the scorer option and the two generators, off by default; V1, V2, V4 - V6.
3. Feasibility probe (sec. 3), and the confirmation fit if the rule would kill C1.
4. Fits, sequentially R0 then C1 (sec. 4); comparison rule; C1's series-noise continuation; Laplace; profiles.
5. Only after both fits have ended: one generator run scores both rivals on the decisive data
   (`scripts/generators/score_m1_decisive_data.py`, to write; artifact
   `kinetic_core_m1_decisive_comparison.{md,json}`). No decisive number is computed on an intermediate member.
6. Side effects (sec. 8).
7. Results review; answers.
8. Decision and its basis in a new EX entry (next free id, checked against origin/main and open branches); CF entries
   for the compensating errors named in sec. 8; the C3 follow-up and the W6 loss-law data need registered in the
   plan's backlog.

## 11. Differences from the case and the plan

1. **D2 does not separate C1 from C2 on mechanism.** The plan says the cysteine and sugar time courses separate C1
   from C2. C2's term acts on the thiols, not on cysteine or pentose, so on D2 the two predict the same thing. D2 is
   the supply verdict's data.
2. **BTBU is a v2 E12 TEST laboratory.** The case marks Cao 2017, Zhang 2024 and Zhang 2025b "DEV (BTBU)"; true in v1.
   In v2, E12's TEST contrasts include F28#25 (Wang 2022), whose laboratory is `btbu_beijing`. Using BTBU rows to
   choose does not read F28#25, but F28#25 is then not laboratory-independent of the choice: tagged (sec. 8). (M1.0's
   "no TEST contrast" statement is about rows and holds.)
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
6. **M1 is the supply verdict only; C2 is specified, not run (decided 2026-11-04, v3.1).** The plan expects M1 to
   choose the loss law. On disk no development datum discriminates C1 from C2 (sec. 5.1), so the loss verdict is
   "cannot decide on disk", C1 is carried as the simpler form, and C2 is not fitted (saving ~12 - 25 h of compute).
   The loss-law question is registered as a W6 data need (sec. 14). M2's refit starts from C1 (or R0) and may revisit
   C2 if W6 data arrive.
7. **C2's form is EX-050's cysteine factor without the electrophile pool E** (the case quotes k_E E / (1 + Cys / K_c)).
   E has no declared identity or measurement, and EX-050 / EX-055 judged the E form flexible, not identifiable (P 0.1).
8. **`k_mft_decay` has a declared weak prior** (sec. 2.6), not v2's Grosch band and not EX-050's Cys-conditional bound.
9. **R0 is refitted** (sec. 2.1) so that both rivals have the same dimer ceiling and the same `k_mft_decay` prior.
10. **Out of M1, as the plan says:** headspace O2 (B11; case part B; open questions 3 and 6 lapse; Ehrenberg 1989 is
    not used); the pH balance (case part D, moved to M1b; held at B60); B27's source. Rival R4 (laboratory offsets) is
    not a combination; it enters only through H0's tau.
11. **D1's noise:** the case uses the measurement model's per-series noise. M1.0 found the s042 noise is a misfit
    property, so v3 decides on Hofmann 1998's measured replicate noise and reports the s042 and read-error scales as
    sensitivities (sec. 6).
12. **The fit objective** is B60's 52 rows plus the route block (minus the Zhai 140 C pots), not R3En's 521 rows.
13. **"Ships" is provisional until M2** (sec. 7): the plan's hard dependency "M1 and M2 in one refit".
14. **The case's D1 row** says R0 "falls after ~1 h"; at B60 it peaks at 2.6 h (M1.0). Immaterial to the design.
15. **The dimer's numerical bound is replaced by a physical ceiling** (sec. 2.5), a fit-side change in both rivals.
16. **Grosch 1992 is not used (v3.1).** Table VIII's reflux values are Table II nos. 1 and 4, E18's TEST contrast
    F27#15 and the keyed `heldout_grosch1992_tab2.yml` (CF-089). v3 stopped banding `k_mft_decay` on it; v3.1 drops
    v3's D2b, which could not move the loss verdict. Grosch Table II / VIII stays unread by M1's choice and fit, so
    nothing of E18 is spent; the side-effect benchmark scores it as TEST after the verdict, as for every keyed table.
17. **BTBU and Reading are M2's declared held-out laboratories.** Split v2 lists `btbu_beijing`, `reading_mottram_parker`
    and the study `whitfield1999` as `prereg_declared_heldout_not_dev` for E6 and E11 (M2's errors). M1 fits BTBU rows
    (the route block's Cao 2017 and Zhao 2017b) and chooses on BTBU rows (D2-in, D2p), and M1 and M2 are one refit, so
    M2's BTBU TEST contrasts (F23#5 Wang 2012, F23#7, F24#4, F25#23 Zhao 2019b) lose laboratory independence:
    registered here, tagged in sec. 8. No row of those contrasts is read. The alternative, dropping BTBU, is not taken:
    without Cao, Zhao and Zhang 2025b the route has no timescale (B58) and R0 has no representable supply data besides
    D1. Whitfield 1999 is already in R0's objective (B60); M1 adds no Whitfield row.

## 12. Compute estimate

**Measured basis:** B60, the same 52 rows, quick setting: start 0 7601.5 s / 246 evaluations = 30.9 s per evaluation;
start 1 8600.2 s / 298 = 28.9 s (`kinetic_core_b60_fit_report.json`). v2's "~13 s with 4 workers" was not a
measurement and is withdrawn. C1's objective adds 34 route pots; its evaluation is estimated at **40 - 60 s**
(UNVERIFIED; W1.3 and V4 measure it before any fit is scheduled). Machine: 16 GB RAM, ~6 GB usable; every job `nice`d,
4 workers, one fit at a time.

| step | count | estimate (unattended unless stated) |
|---|---|---|
| build + tests (author work) | the dimer ceiling, a scorer option, two generators | 4 - 7 h of work |
| W1.3 profiling + V1, V2, V4 - V6 | one evaluation profiled; B60 cost, anchors, timings | ~2 h of work + ~0.5 h |
| feasibility probe | 256 draws x 3 pots, quick, + a local search | 0.3 - 0.5 h; + one start (~3 - 5 h) if the rule would kill C1 |
| R0 fits | 3 starts x 250 - 300 evaluations x 29 - 31 s | 6 - 8 h |
| C1 fits | 3 starts x 250 - 330 evaluations x 40 - 60 s | 8 - 16 h |
| budget continuations | ~40 % of starts, ~half a start each | 3 - 6 h |
| C1's series-noise continuation | <= 300 evaluations | 3.5 - 5 h |
| Laplace | 2 x (2n + 1) evaluations | ~1 h |
| profiles | expected 0 - 2 coordinates per rival, 3 points x <= 100 evaluations | 0 - 7 h (typical ~3.5) |
| decisive data | 64 draws x 3 rival forms (R0, R0-shipped, C1) x 3 pots, careful setting, + sensitivities | 1 - 1.5 h |
| side effects | the contrast benchmark x 2 new configurations, in-sample rows | 1.5 - 3 h (UNVERIFIED: the scorer's run time) |
| reviews | v3.1 re-review + results review + answers | 4 - 5 h of work |
| **total** | | **~10 - 14 h of work; ~25 - 50 h unattended compute typical (~36 h); worst ~110 h if most starts hit the 1 200 cap and need their continuation; about 3 - 4 calendar days** |

v3's estimate with C2 was ~45 - 75 h typical (~60 h), ~170 h worst. Memory: one evaluation with 4 workers is expected
to fit in ~6 GB; recorded in V4. If not, 3 workers (~+30 % wall time).

## 13. Expectations (written now; not blind, sec. 0)

| event | P |
|---|---|
| C1 passes the feasibility rule without a confirmation fit | 0.65 |
| R0 / C1 compared (the lowest converged start agrees with another) | 0.65 / 0.6 |
| verdict: C1 favoured / R0 favoured (incl. the R0-shipped fallback) / no verdict (C1 not compared) | 0.30 / 0.30 / 0.40 |
| the verdict on the out-of-sample units alone agrees with the full subset's | 0.8 |
| the dominance flag fires | 0.6 |
| D1 at C1's optimum: MFT 360/30 within 0.3 dex of measured / FFT 360/30 within 0.3 dex | 0.35 / 0.25 |
| D1 at C1's optimum: MFT rises from 6 to 12 h (720/360 > 1) | 0.45 |
| lambda_mft(100 C) at C1's optimum within 0.5 dex of 0.041 /min | 0.6 |
| `Ea_decay_thiol_sink` on its 102 ceiling in at least one rival | 0.5 |
| `k_dimer_mft` above log10 0.5 in both refit rivals | 0.85 |
| C1's 52-row cost below R0's refit cost | 0.3 |
| at least one held-out benchmark row gets worse by more than its noise (favoured rival vs R0-shipped) | 0.85 |
| D2p: TTCA within stage 1's sigma (0.20 dex) at C1's optimum / CAMA | 0.60 / 0.30 |
| a W6 dataset that decides C1 vs C2 exists by M2's refit | 0.10 |
| M1 reaches a decision by 11-20 (plan sec. 8) | 0.35 |

## 14. What this does not do

**The loss law, registered as a data need for plan W6.** The question: C1 vs C2, does free cysteine slow the thiol
loss at process temperature? The measurement that would decide it (from the case's open question 2 and the review's
M4 / M5): **a fed-MFT (and FFT) loss rate at process temperature (100 - 145 C, pH 5 - 6), with and without free
cysteine (0 and at least two levels between ~1 and ~50 mM), in one study, from a development laboratory**, with the
thiol quantified by SIDA (labelled standard added after heating, so the in-pot loss is not compensated) over a time
course, and free cysteine measured in the same pots. Best in a reacted ribose + cysteine background, so the pot's
electrophiles are present; a buffer-only arm separates the sink from the matrix. C2 predicts the loss falling as
1 / (1 + [Cys] / K_c); C1 predicts no cysteine dependence. Nothing on disk meets it: Grosch 1992 is held out and its
Table III confounds formation with loss; van Seeventer 2001 is 50 C storage; the Nestec patent is not
measurement-grade.

Also not done: no disulfide side (the C3 follow-up, sec. 2.4); no oxidant from headspace O2; no pH-balance change; no
MFT route shares (M2); no posterior (W4). It does not test the Yiltirak vs Schieberle level gap (not separable on
disk). It does not settle the Cao vs Zhai Amadori-formation tension (phosphate or laboratory; stage 1); `ea_form`
absorbs part of it, which the results review is asked to weigh. It does not test the fast 145 C sink (sec. 2.6).

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
  rule with a stacking guard; tau from H0 or wait; series noise re-estimated; compute from B60's measurement. Every
  finding's answer is in sec. 16.1. Owner's decisions of 2026-11-04 (in the session that commissioned v3): items 1 -
  12 of that brief.
- v3.1 (2026-11-04, before the re-review and any build): D2b dropped; C2 specified, not run; M1 = the supply verdict.
  Consequences: Grosch 1992 used for nothing (no E18 spend); the loss-law question registered as a W6 data need (sec.
  14); C2's key not built; compute and expectations revised (secs. 12, 13). Owner's decisions relayed by the
  coordinator, 2026-11-04.

### 16.1 The design review's findings and the answers (v3, as amended by v3.1)

| finding | what the review found | change (section) or why not |
|---|---|---|
| B1 | sec. 7 summed D1 without saying which sigma; MFT and FFT point opposite ways at C1's start | one decision sigma from measured replicate noise (Hofmann 1998, 0.059 - 0.061 dex per ratio, bivariate per species); read error and s042 printed as sensitivities (6); the start values recorded as seen (0) |
| B2 | non-discriminating data (D2p's 31 points) summed into L would decide by compensation; R0's TTCA / CAMA terms are log(0) | one verdict on a pre-declared subset; D2p reduced to free Cys and pentose (16 points); structural zeros defined from the network and reported as "R0 cannot represent the measured reservoir (EX-060)", never summed; the per-unit SE stops one-sided sums from deciding (5, 7) |
| B3 | C3 built to fail: the trapped release arm, a net-sink branch, a mechanism-contradicting law, false nesting | C3 withdrawn; the follow-up's four requirements registered (2.4) |
| B4 | `k_dimer_mft` sits on a numerical bound that holds R0 / C1 / C2 below the shares while C3's `k_ox` floats | the bound replaced in both rivals by the encounter-limit ceiling, log10 9.47 at 145 C (2.5); R0 refitted under it (2.1, 4.3); V5 |
| B5 | D3 remnants (sec. 5 D4 row, sec. 7, sec. 13, the title "v1"); the five Zhai 140 C pots (30 rows) in the objective though declared held out | every D3 / D4 - D6 reference rewritten (5 - 8, 13); title updated; the five pots named and removed (4.1) and listed as hold-outs (8) |
| M1 | Grosch Table VIII = Table II nos. 1 and 4 = E18's F27#15 and a keyed held-out table; BTBU declared held out for E11; M1 and M2 are one refit | Grosch used for nothing (2.6; v3.1 drops D2b); BTBU's lost independence registered (11 item 17) and tagged (8) |
| M2 | the `k_mft_decay` band is an inference from another matrix, with a chosen sd and C2's centre, from D2b's table, and with Ea near 102 it keeps the fast 145 C sink | weak prior centred on the incumbent, physical basis and identifiability cost stated; the fast 145 C sink stated plainly; `Ea_decay_thiol_sink` freed inside Gigl's band, with the reason (2.6) |
| M3 | Zhou 2022 replicates an in-sample pot; Whitfield counted twice | moot for the verdict (those data left the decisive set with C3); in the side-effect record Zhou 2022 is marked not independent and each Whitfield arm is counted once (8) |
| M4 | D2b spans [-0.22, 0] for C2, so it can only support C2; sigma 0.10 - 0.13, not 0.07; computation unspecified | v3 defined it (sigma 0.149, one-sided Delta L <= ~0.28); v3.1 drops it: it cannot move the loss verdict and is held-out E18 data (5.1, 11 item 16) |
| M5 | free Cys is near-flat in Schieberle, so D1 cannot separate C1 from C2 | accepted: with D2b gone the loss verdict has no unit; C2 is specified, not run; the deciding measurement is registered for W6 (2.3, 5.1, 14) |
| M6 | C3's pool is model-internal (TDP 1.6 -> 40 mM) | moot (C3 withdrawn); the follow-up must put the measured DPO / TDP in its objective (2.4) |
| M7 | D6 is set by the declared ambient OX 1.0 mM | moot (C3 withdrawn); the follow-up drops D6 or carries OX as a sensitivity (2.4) |
| M8 | the tie rule is intransitive; one sd chosen; the Laplace sd is invalid at active bounds | "simplest within one SE of the best", SE from per-unit Delta L with posterior-predictive averaging, bound-active coordinates profiled (7) |
| M9 | the convergence rule is ambiguous | one rule: compared iff the lowest converged start agrees with another under `two_starts`' rule, else "not compared, flagged"; fallbacks declared; three starts per rival (4.3) |
| M10 | the feasibility screen has hidden thresholds and an undefined kill rule; S2 is circular | one kill rule (not better than R0 on the supply data under the decision sigma); a confirmation fit guards against stacking; S1 / S2 printed only (3) |
| M11 | tau = 0 while EX-059 is pending; s_series estimated at the start | tau_ABS(all) from H0 or the scoring waits; s_series estimated at C1's start and re-estimated at its optimum with a continuation (4.3, 6) |
| M12 | C1's start values on D1 and D2p were seen | marked "seen at start" (0, 5); the claim weakened to out-of-objective, not blind (5, 13) |
| MINOR: compute | ~40 - 60 h typical, not 20 - 25 (B60 ran ~31 s / evaluation) | measured basis stated; with C2 dropped, 25 - 50 h typical (~36 h), worst ~110 h; W1.3 and V4 as prerequisites (10, 12) |
| MINOR: route T-dependence | the route's temperature law was not stated | the route's barriers and T-flat K declared, the 140 C extrapolation named, K and `ea_decay` band-end sensitivities (2.2, 7) |
| MINOR: pH exposure | which pH the rivals' pH-flat pieces meet was not stated | the route's pH-flat pieces across pH 5.5 - 8 and the supply data's pH 5.0 - 6.5 stated; a fixed-label-pH sensitivity (2.2, 7) |
| MINOR: tie-break order | the order was asserted, not defined | the order is the fitted-coordinate count, R0 17 < C1 19 (C2 20, not run) (2.1, 7) |
| MINOR items on C3 | -- | moot (C3 withdrawn) |
| NOTE | -- | not itemised in the brief v3 was written from; no change made for it |

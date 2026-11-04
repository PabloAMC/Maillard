# M1 pre-registration: thiol supply and loss, judged as three combinations (v1)

*Plan step M1 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, row M1; sec. 7 items 6 - 9). Case:
[M1_thiol_loss_oxidant_ph_drift.md](../../docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md) (items 1 - 5; not
edited here). Conflict: CF-107 ([conflicts_register.yml](../../docs/reference/conflicts_register.yml)). M1.0: EX-060,
[kinetic_core_m1_preflight_checks.md](kinetic_core_m1_preflight_checks.md).*

- v1 written 2026-11-04, BEFORE any M1 build, fit or score. No engine was run with a changed switch or constant to
  write it. Nothing above the "Results" line (sec. 15) is edited after the first M1 number exists; changes before that
  are versioned (v2, ...) and logged in sec. 16.
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
- **The route at B60's constants (= C1's start vector):** B59b Amendment 3 (Schieberle MFT over 30 min 6.1 / 84 / 73
  at 60 / 360 / 720 min; Nestle 110 C shape; Zhang 2025b TTCA +0.15 to -0.37 dex; Hofmann pH 7 MFT -1.6 dex; Zhai
  2019 Cys-Amadori over-made 2 - 5x; Kang fed-TTCA MFT down 260x); the route sub-model no-refit probe
  ([kinetic_core_route_submodel_probe_with_b58_with_b59b_with_b59a_with_b61.md](kinetic_core_route_submodel_probe_with_b58_with_b59b_with_b59a_with_b61.md),
  per-data-set table only: Cao RMS 0.03 - 0.31 dex; Zhang 2025b TTCA -0.10 / CAMA +0.66 / Cys +0.55 / PENT +0.44 dex).
- **The stage-1 route calibration**
  ([kinetic_core_cysteine_route_calibration_prereg.md](kinetic_core_cysteine_route_calibration_prereg.md)): its
  constants, its checks (Zhang 2025b, Zhai 2019 Fig. 3B, Zhao 2017, pH), and its flavour-level diagnosis (Schieberle
  12.7 / 345 / 171 over 30 min).
- **R3En** ([kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.md](kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.md)):
  residuals by data set, including Schieberle MFT / FFT, Wang 2026b MFTD / FFTD and Zhou 2022 MFTD (after a free
  scale), the declared hold-outs (Yiltirak, Hofmann T1 / T2, Zhai 140 C) and the score-only rows.
- **B60:** its prereg and `residual_by_row` for all 52 rows (e.g. `zhou_pH7_dimer_over_MFT` -2.71,
  `zhang_fig1_cys_dimer_over_MFT` -1.99, `kang_120C_MFT` +2.51 sigma); B58's ship rule; B27's G1 / G2 as quoted in
  the case.
- Register entries EX-024, EX-027, EX-050 (and its README), EX-055 (as quoted in EX-050), EX-060; CF-107.
- **Measured values only (no engine number):** Zhou 2023 Table 1 at pH 6 / 7 / 8 (MFT, FFT, MFTD, FUR, with SDs);
  Whitfield 1999 Table 1 row 48 (bis(2-methyl-3-furyl) disulfide 2 / 3 ug Cys arm, 7 / 2 ug H2S arm) and the dossier's
  7.5 % / 35 % shares; the F28 atlas cells #1, #2, #9 - #12; the de Roos 2005 keyed-table header (no value).
- Not read: any engine number on the Zhou 2023 pH 6 / 8 disulfide shares, on Zhang 2024's free-MFT ratio, on
  Whitfield's disulfides, or on Grosch 1992 Table VIII.

## 1. Question

Which combination of thiol supply and thiol loss do the mechanism-direct data favour, judged together because the
errors stack (CF-107): the shipped engine R0, or C1, C2 or C3? All three combinations are fitted and compared at
once. No combination is judged as "C1 plus one addition".

## 2. The rivals, exactly

### 2.1 Config keys

| rival | switches ON (Config key) | new step(s) | removed step(s) |
|---|---|---|---|
| R0 | none (the shipped configuration, B60 vector) | -- | -- |
| C1 | `sulfur.cys_amadori_route` (`enabled = True`; it also turns on `r_ttca_form` by detailed balance with `sulfur.ttca_formation`'s K) | `r_ttca_form`, `r_cama_form`, `r_cama_dpo`, `r_cama_tdp` (exist, B58 / B59b) | `r_pent_dpo`, `r_pent_tdp` (retired by the route) |
| C2 | C1 + **`sulfur.cys_suppressed_thiol_loss`** (NEW key, to build; off by default and bit-identical when off) | none: a factor on `ch_mft_decay` and `ch_fft_decay` | -- |
| C3 | C2 + **`sulfur.disulfide_release`** (NEW key, to build; replaces the inert `k_dimer_release` table entry) + **`sulfur.osone_thiol_oxidation`** (NEW key, to build) | `ch_dimer_release_mft_2`, `_mft_1`, `_fft_2`, `_fft_1` (the existing release steps, split by yield, sec. 2.4); `ch_osone_ox_mft`, `ch_osone_ox_fft` | the B17 single-constant release (`k_dimer_release`, inert today) |

Names of the three new keys are proposals; the build may rename them, and the rename is logged in sec. 16.

### 2.2 Constants every combination shares (C1, C2, C3)

| constant | step | status | value, band, source |
|---|---|---|---|
| K, TTCA equilibrium | `r_ttca_form` = K x `k_ttca_cys` | **declared, measured** | 0.060 L/mmol; band 0.045 - 0.075 (Cao 2017, `TtcaFormation`). Held at 0.060. pH- and T-flat (B58's declared assumption) |
| **`k_ttca_cys`**, ring opening TTCA -> Cys + PENT | `r_ttca_cys` | **FITTED** (in the structure, not frozen; CF-107) | start B60 (log10 -0.579 at 145 C); bounds +-3 decades about the start (stage-1 convention). Its Ea stays on `Ea_lumped_formation` (64.08 kJ/mol, B60, held). Check, not a bound: the Cao dossier's floor >= 7e-3 /min at 90 C |
| `k_form`, TTCA -> Cys-Amadori at 90 C | `r_cama_form` | **declared** (anchored on Cao 2017's blank, B59b V2) | 3.29e-3 /min; probe band +-0.44 dex (the anchor vs stage 1's 1.2e-3) |
| **`ea_form_kj_mol`** | `r_cama_form` | **FITTED** | start 100.9 (Zhang 2026, glucose system: a transfer B59b Amendment 3 found too steep for pentose); bounds 40 - 135 kJ/mol (`sugar_trunk` class band, as stage 1). Needs a field band in `config._BANDS["sulfur.cys_amadori_route"]` before any fit (the key refuses an unbanded field) |
| Cys-Amadori decay k(pH) at ~100 C | `r_cama_dpo`, `r_cama_tdp` | **declared, measured** | de Roos 2005 Fig. 4 table; above pH 6 `extend_last_slope` (CR-14; stage 1 found the two arms identical). See sec. 11, item 3: this table is a keyed E10 TEST table |
| `ea_decay_kj_mol` | same | **declared** (transfer) | 85.7 (Zhang 2026); probe band 40 - 135 |
| `dpo_share` | same | **declared, measured** | 0.30 (Zhai 2019; alternatives 0.28 / 0.32, CR-12); `split_law = "none"` |
| **`k_mft_decay`** | `ch_mft_decay` | **BANDED** (sec. 2.6) | Gaussian prior on log10 lambda(100 C): centre log10 0.044 /min, sd 0.30 dex (Grosch 1992 Table VIII). Moves from the frozen set into the free set |

### 2.3 C2's added constant

| constant | status | form, value, band, source |
|---|---|---|
| `K_c`, the cysteine level that halves the loss | **FITTED** | `ch_mft_decay` and `ch_fft_decay` run at k / (1 + [Cys]_free / K_c). [Cys]_free is the free-cysteine species, not TTCA- or CAMA-bound. log10 K_c (mmol/L) start 1.5, bounds 0 - 3 (1 - 1000 mM). No measurement pins it at process temperature: van Seeventer 2001 (+150 mM Cys) is a 50 C storage cell and cannot set a processing law (atlas rule); Grosch 1992 Table VIII is a check (D2b), not a fit row |

The factor does NOT act on `ch_thiolate_loss_*` (the Kumazawa-anchored step) or on the dimer. `k_mft_decay`'s band
(sec. 2.6) is then lambda at zero free cysteine, which is what Grosch's thiamine arm measures.

### 2.4 C3's added constants

**Release (Guth 1995, declared, not fitted).**

| constant | status | value, band, source |
|---|---|---|
| `k_rel_mft`, MFTD lost at 100 C | **declared, measured** | 0.0165 /min (86 % gone in 2 h, 4-vinylpyridine trap). Probe band 0.0134 - 0.0210 /min: the author's assumption of +-6 points on the 86 % (Guth prints no replicate), UNVERIFIED as a noise statement |
| yield, MFT per MFTD lost | **declared, measured** | 1.43 (62 % of two thiols back over 86 % lost; D2O shows hydrolysis RSSR + H2O -> RSH + RSOH). Built as two integer steps: fraction 0.43 -> 2 MFT, 0.57 -> 1 MFT + FRAG (the sulfenic acid leaves) |
| `k_rel_fft` at 100 C | **declared, a lower bound** | 0.034 /min (>= 98 % gone). Probe band 0.033 - 0.10 /min; the upper end is an assumption |
| yield, FFT per FFTD lost | **declared, measured** | 1.88 (92 % of two back): fraction 0.88 -> 2 FFT, 0.12 -> 1 FFT + FRAG |
| `Ea_rel` | **declared** (transfer) | 122.2 kJ/mol, the corpus's only measured thiol-to-disulfide barrier (Zhang 2026 k17, Cys-Gly), already on `ch_dimer_*`, as B17 variant b: the equilibrium is then T-flat. Probe band 60 - 140 kJ/mol, an assumption (no disulfide-hydrolysis barrier on disk, UNVERIFIED) |
| pH law of the release | **declared** | none (pH-flat): Guth does not state the pH. A base-catalysed hydrolysis is the named alternative; it is not run |

**Oxidant: the deoxypentosones as a POOL (not B27's co-product source).**

| item | status | declaration |
|---|---|---|
| rate law | **declared** | r = k_ox x theta(pH, T) x [RSH] x ([DPO] + [TDP]), for RSH = MFT and FFT. First order in thiol, first order in the osone pool. Each event consumes 2 RSH and 1 osone and makes 1 RSSR; the reduced osone (the hydroxy-ketone of Whitfield 1999 Fig. 6) leaves as FRAG_C 5 |
| why first order in thiol | -- | the rate-limiting step is the thiolate's addition to the carbonyl (one thiol, one dicarbonyl); the second thiol's attack is fast. The engine's `ch_dimer_*` is second order in a trace thiol, which is why it cannot reach the share at any oxidant level the pot holds (G1). The engine computes rates by mass action from stoichiometry (`network.Reaction.order`), so the build needs a rate-order field; it reports whether this changes any existing step (it must not) |
| pH law | **declared** | theta = the thiolate fraction at reaction temperature, the factor `ch_dimer_*` and `ch_thiolate_loss_*` already use (`sulfur.REACTION_PH_FACTOR` "thiolate"), each thiol's own pKa. The pool's own pH dependence comes from the network. A pH-flat law is computed once at C3's optimum, no refit, as a recorded sensitivity; it is never used to choose |
| `k_ox` at 145 C | **FITTED** | one constant shared by MFT and FFT. Start: the value at which the new channel's MFT flux equals `ch_dimer_mft`'s in the Zhou 2023 pH 7 pot at B60 (data-free: "as large as the existing channel"). Bounds +-4 decades about the start |
| `Ea_ox` | **declared** (transfer) | 122.2 kJ/mol (Zhang 2026 k17, Cys-Gly thiol to disulfide: the same bond formed). Probe band 60 - 140 |
| `ch_dimer_*` with OX | unchanged | kept: OX carries the ambient 1.0 mM and the cystine of Zhang 2024's GCys arm (thiol-disulfide exchange is a different mechanism); zero in the fed mol % rows as today |

C3 therefore adds one fitted coordinate (`k_ox`) and no fitted release constant.

### 2.5 Held fixed at B60 (every rival)

- **Free, from B60's start:** the 16 coordinates B60's Laplace marks identified (B59b Amendment 2's rule):
  `k_tdp_fur`, `k_nf_mft`, `k_nf_mp3p`, `k_mgo_mp`, `k_fft_decay`, `k_fur_decay`, `k_pent_caramel`, `k_cys_thermal`,
  `k_thiolate_loss`, `ph_arp_secondary_ammonium_pKa`, `k_osone_decay`, `k_ttca_deg`, `Ea_decay_thiol_sink`,
  `k_dimer_mft`, `k_arp_dpo`, `k_arp_tdp`. B60's bounds.
- **Held at B60's values:** the 7 unidentified at B60 (`k_glc_ha`, `k_dimer_fft`, `ph_acid_yield_per_sink_event`,
  `k_fur_fft`, `k_dimer_decay`, `k_thiol_decay`, `Ea_decay_carbonyl_sink`) and the 24 other frozen keys of B60's
  report, except `k_mft_decay` (now banded) and `k_ttca_cys` (now fitted). `k_pent_dpo`, `k_pent_tdp` drop out with
  their steps.
- **Switches at their shipped defaults:** `ph.coupling` fixed_point_nodes; `ph.node_grid` log_and_linear;
  `ph.fixed_point` relaxed_converged; `ph.b51` all on; `ph.trunk_solutes_in_balance` False and `ph.trunk_terms`
  False (the pH balance is M1b's); `sulfur.flux_budget` closed; `panel.thermal_history` isothermal; B59a
  (`sulfur.strecker_sulfur`), B61, B62, B63, B7a's split law, `sulfur.tdp_own_decay`, the pentose-amine entry and
  every other sub-model switch off. (R3En ran with B59a and B61 on; M1 does not.)
- **The fit's numerics:** B60's quick setting (17 pH nodes, rtol 1e-5). Scoring of decisive data and side effects:
  the careful setting (81 nodes).

### 2.6 `k_mft_decay`: banded, not freed. The reason

- **Frozen** is ruled out by the plan and by CF-107: the frozen fast sink is the anchor of the degenerate optimum.
- **Free** is ruled out because a free sink can absorb a supply error. EX-024 shows the pattern (the budget-stopped
  E arm found MFT loss ~8x slower; converged, it went back). And a free sink shares its effect on D1 with C2's term,
  which removes the separation C1 vs C2 needs (case, open question 2).
- **Banded on a measurement at the engine's own size:** Grosch 1992 Table VIII (DFA, DEV in v1 and in v2 E12; atlas
  F28#12): in-pot loss at 100 C, pH 5.7, 0.044 (0.035 - 0.053) /min in the thiamine pot, 0.034 (0.030 - 0.038) with
  cysteine. Centre: the thiamine arm (no added cysteine, matching C2's lambda at zero cysteine). Width 0.30 dex,
  wider than the printed interval for two one-sided caveats the atlas cell names: it is a lower bound (equal formation
  under SDE and reflux assumed), and it may be a thiamine-specific electrophile sink (jhoo2002). B60's own lambda(100
  C) is 0.041 /min, inside the band, so the band does not pull B60 away from its value at 100 C; the 145 C value can
  still move through `Ea_decay_thiol_sink`.
- Not used: EX-050 / EX-055's Cys-conditional bound (lambda(110 C) <= 0.02 - 0.03 /min at Cys >= 10 mM), because it
  rests on a patent's four digitised single-vial points (not measurement-grade, P 0.15) and because a Cys-conditional
  band in C1 would build C2's term into C1.

## 3. Mechanism-signature feasibility first (no fit)

The owner's rule: test the structures on the mechanism signatures before any constant is fitted.

**Signatures.**
- **S1, the Schieberle rise to 12 h** (Schieberle 2000 Table IV; 100 C, 0.5 M phosphate pH 5, ribose 100 + Cys 33
  mM): MFT(720) / MFT(360) > 1 (the pot still rises) AND log10 MFT(360) / MFT(30) within 0.80 dex of the measured
  +1.54. The 0.80 dex is the s042 MFT series noise, the most lenient noise on file; a feasibility test asks "at all",
  so it takes the lenient one. R0 at B60: 0.54 and +0.45 (fails both).
- **S2, the Cao reservoir** (Cao 2017 blank, 90 C, 5 h): (free Cys + TTCA + CAMA) / Cys charge and (free Xyl + TTCA +
  CAMA) / Xyl charge each within 0.20 dex of the measured 1.06 and 1.09, and free xylose at 5 h within 0.20 dex of
  0.172. The 0.20 dex is the common sigma stage 1 estimated on these route rows. R0: 0.67 and 0.00 (fails).

**Probe.** Per structure, 256 Sobol draws over its bands, the other constants at B60:

| structure | coordinates drawn (bands) |
|---|---|
| C1 | `k_ttca_cys` (+-3 dex), `ea_form_kj_mol` (40 - 135), K (0.045 - 0.075), `k_form` (+-0.44 dex), `ea_decay_kj_mol` (40 - 135), lambda(100 C) of `k_mft_decay` (0.005 - 0.10 /min), `Ea_decay_thiol_sink` (7 - 102) |
| C2 | C1's + log10 K_c (0 - 3) |
| C3 | C2's + `k_rel_mft`, `k_rel_fft`, `Ea_rel`, `Ea_ox` (sec. 2.4 bands) + log10 `k_ox` (+-4 dex) |

A structure **passes** if one draw meets S1 and S2 together. If no draw does, a local search on the signature
residuals alone (least squares from the best draw, <= 200 evaluations, inside the same bands) is run before the
verdict, so a thin feasible region is not missed. A structure that fails both is **out, without fitting**, and the
record says which signature it cannot meet. R0 is not probed: its failure is M1.0's, and it stays in the comparison
as the incumbent.

The probe tests the structures on their new and banded constants only; the 16 B60 coordinates stay at B60. C1 is
nearly nested in C2 (K_c -> 1000 mM) and C2 in C3 (`k_ox` -> its floor; the release is then nearly inert, because
the disulfide made by OX is tiny), so C2 and C3 pass wherever C1 passes, approximately. The probe mostly tests C1
against R0's failure.

Artifact: `results/validation/kinetic_core_m1_feasibility_probe.{md,json}`, generator
`scripts/generators/m1_combination_feasibility.py` (to write). Seed fixed and printed.

## 4. The fit rows and the fit

### 4.1 Objective (identical for C1, C2, C3)

| block | rows | measurement model |
|---|---|---|
| **B60's 52 rows** (R0's own objective) | exactly B60's `row_ids` | B60's: its sigmas, `weight_tag` half, pH endpoint weight E 0.7 decades per pH unit, sigma_pH 0.5 |
| **route block** | every verified row of [route_submodel.yml](../../data/lit/fit_rows/route_submodel.yml) (Cao 2017 Figs. 4 - 5; Zhai 2019 Fig. 3A; Zhai 2021 Figs. 3a, 4b, 6; Zhao 2017b; Zhai 2023 Food Chem / JAFC TTCA pots), except: the Zhang 2025b pots (`zhang2025b_blank`, `zhang2025b_blank_fig2gh`: kept out as D2's predictive part) and every pH row | per row sigma^2 = read_error^2 (or the printed SD in dex) + s_series^2, with s_series per data set x species estimated ONCE by `fitting.estimate_series_noise` at C1's start vector and then fixed for every combination and start |
| **prior** | one term: (log10 lambda_mft(100 C) - log10 0.044) / 0.30 | sec. 2.6 |

**Why not R3En's 521 rows.** R3En has Schieberle 2000 (D1) and Wang 2026b / Zhou 2022 disulfides (D4) as fit rows,
runs B59a and B61 on, and lacks the fed-pot disulfide rows and Kang rows of R0's objective. B60's rows keep the
comparison with R0 like-for-like; the route block adds the rows that measure the route's own species, without which
no fit row measures a timescale (B58) and `k_ttca_cys` is not identified.

**Before the build, two checks on the rows:**
- *Duplicates.* `lab_ids.yml` lists `kang2026_SI` as a copy of `zhai2023foodchem`. B60's Kang rows
  (`kang_*_cys_conversion`, `kang_free_cys_yield_ceiling`) may be the same data as the route block's Zhai 2023 Food
  Chem pots. If they are, the route block drops those pots (B60's rows stay, so R0's objective remains a subset).
- *Weight.* The route block has ~300 rows against 52. The fixed series noise keeps a misfitting series from
  dominating, but CF-107 found 77 % of B16's cost on TTCA rows. Each optimum reports the share of its cost from each
  block. This is a record, not a gate.

Not in the objective (scored only as decisive data, sec. 5): Schieberle 2000; Zhang 2025b's route species; Zhou 2023
pH 6 and 8; Zhou 2022 and Wang 2026b disulfides; Whitfield 1999's disulfides; Grosch 1992 Table VIII; Zhang 2024's
free-MFT ratio between arms.

### 4.2 Free coordinates

| | B60's 16 | `k_mft_decay` (banded) | `k_ttca_cys` | `ea_form` | K_c | `k_ox` | total |
|---|---|---|---|---|---|---|---|
| C1 | yes | yes | yes | yes | -- | -- | 19 |
| C2 | yes | yes | yes | yes | yes | -- | 20 |
| C3 | yes | yes | yes | yes | yes | yes | 21 |

### 4.3 Optimiser, budget, starts, convergence

- `src/kinetic_core/fitting.py`: `fit` / `two_starts`, `LEAST_SQUARES_SETTINGS` (trf, x_scale jac, ftol 1e-6, xtol
  1e-8), the module's absolute finite-difference steps (never a relative `diff_step`).
- **Budget:** 1 200 evaluations per start.
- **Starts per combination:** (0) B60's vector with the new coordinates at their declared starts; (1) B8's
  perturbation protocol on that vector (as B60's start 1); (2) for C2 and C3 only, the next-simpler combination's
  best optimum with the new coordinate at its start (K_c at 1000 mM; `k_ox` at its floor).
- **A start counts only with `least_squares` status > 0.** A start that ends at status 0 (budget) is continued once
  from its endpoint with another 1 200 evaluations. Still status 0: reported "not converged" and not used.
- **A combination is converged** when at least two starts reach status > 0 and `two_starts`' agreement holds (costs
  within the module's `COST_AGREEMENT`, every identified constant within its Laplace sigma). Disagreeing starts:
  both reported, the lower carried, flagged. A combination with no converged start is not compared (plan sec. 7 item
  6: only converged fits count); the record says so.
- **Laplace** at each carried optimum (the module's central differences and identification rule). Bounds hit are
  reported as such.
- **Nesting check (diagnostic):** C2's cost must not exceed C1's, nor C3's C2's, by more than the near-nesting gap
  measured at start (2). If it does, the higher fit is not at its optimum, and is continued once.

### 4.4 Verification before any fit

- V1: all new keys off: every existing test passes and B60's 52-row cost reproduces (27.132 at the quick setting).
- V2: C1's switch on at B60's vector reproduces the route sub-model probe's Cao residuals (the B59b anchor).
- V3: C2's factor is exactly 1 at [Cys]_free = 0; C3's release, run on Guth's pot with a thiol sink standing in for
  the trap, reproduces 86 % of MFTD gone and 1.43 thiols per MFTD lost at 100 C / 2 h; C3's steps conserve C and S
  (closed flux budget, <= 1 %).
- V4: peak memory of one evaluation with 4 workers, measured and recorded (sec. 12).

## 5. The decisive data per comparison

Every row below was checked against [split.yml](../../docs/reference/contrast_benchmark/split.yml) (v1) and
[split_v2_per_change.yml](../../docs/reference/contrast_benchmark/split_v2_per_change.yml) (E10 and E12, the errors
M1 addresses) and against [FIT_HOLDOUT_DECLARATION.md](../../docs/reference/FIT_HOLDOUT_DECLARATION.md). "In" = in the
objective of every C (in-sample); "out" = never fitted (predictive).

| # | comparison | data (record) | in / out | lab: v1; v2 E10 / E12 | other status |
|---|---|---|---|---|---|
| D1 | R0 vs C1 - C3; weakly C1 vs C2 | Schieberle 2000 Table IV: MFT and FFT ratios 360/30 and 720/360 (720/30 is their sum, reported only) | out | DFA: DEV; E10 not DEV, no TEST contrast; E12 DEV | calibration input (never a v2 TEST contrast) |
| D2 | requirement on all C; R0 fails it (M1.0) | Cao 2017 and Zhai 2019 Fig. 3A route species | in | BTBU, ho_network: DEV; E10 DEV; E12 neither DEV nor TEST (see sec. 11 item 2) | -- |
| D2p | same, predictive | Zhang 2025b Fig. 2 route species (Cys, Xyl, TTCA, CAMA, 10 - 50 min) | out | BTBU: DEV; as D2 | its Table 2 volatiles are a frozen hold-out bundle (Amendment 45): NOT read |
| D2b | C1 vs C2 | Grosch 1992 Table VIII: lambda(+Cys) / lambda(thiamine) = 0.034 / 0.044 | out | DFA: DEV; E12 DEV (F28#12 is a cited DEV cell) | not the keyed held-out Table 2 |
| D3 (NOT USED, v2: held-out, sec. 11 item 4) | C3 vs C1 / C2 / R0 (pH law of the oxidation) | Zhou 2023 Table 1: log10 of the MFT disulfide share (thiol equivalents) at pH 7 minus pH 6, and pH 8 minus pH 7 | out | ho_network: DEV; E12: excluded as calibration input (F28#4) | **pH 6 / 8 columns are a declared HOLD-OUT in FIT_HOLDOUT_DECLARATION.md**: see sec. 11 item 4 |
| D4 | C3 vs C1 / C2 / R0 (share level) | out: Zhou 2023 pH 6 and 8 shares; Zhou 2022 MFTD / MFT (4 times); Wang 2026b MFTD / MFT and FFTD / FFT (60 - 150 min); Whitfield 1999 Cys arm MFT disulfide / MFT. In (reported apart): Zhou 2023 pH 7, Zhang 2024 Fig. 1 cysteine arm | mixed | ho_network, Reading: DEV; Wang 2026b (Dalian): in no split, development by use | Zhang 2024 Fig. 1: FIT in the declaration; Figs. 2d - f (115 C) are a HOLD-OUT and are not used |
| D5 | C3's release vs C1 / C2 / R0 | Zhang 2024 Fig. 1: free MFT, cystine arm over cysteine arm (x1.5 measured). The share ratio (x6.3) is in-sample, reported apart | out (ratio of free MFT) | BTBU: DEV; E12: excluded as calibration input (F28#5) | figure read-off only |
| D6 | C3's oxidant reaches fed pots | Whitfield 1999: bis(MFT) disulfide / MFT in the NF + H2S arm over the NF + Cys arm, and each arm's level | out | Reading: DEV; E12 DEV | FIT source in the declaration (its three MFT rows are in B60's objective; the disulfide columns are not) |

Not decisive, recorded only: Yiltirak 2026 vs Schieberle (case D7: not separable on disk); Kumazawa 2003 (v1 TEST; a
constraint already in every objective); the pH controls (M1b's).

### 5.1 What each rival predicts, from its mechanism

Directions are read from the structures, not computed.

| # | R0 (B60) | C1 | C2 | C3 |
|---|---|---|---|---|
| D1 | MFT 360/30 short by 1.09 dex; falls 6 -> 12 h (M1.0) | the reservoir feeds DDP for hours: 360/30 near the measured or above; 720/360 near 1, falling if the sink outruns the drain | as C1, loss slowed where free Cys is high: 720/360 >= C1's | as C2; release adds back at most ~1.1x (bounded by the share); the osone oxidation is small at pH 5 (thiolate) |
| D2 / D2p | no reservoir; pentose gone in seconds (M1.0) | the reservoir holds for hours | = C1 (its term does not touch Cys or PENT) | = C2 (the oxidation consumes a negligible part of the osone pool) |
| D2b | ratio 1 (pot-independent) | ratio 1 | 1 / (1 + [Cys]_free / K_c) < 1 | as C2 |
| D3 | share rises with pH (thiolate law on a second-order dimer) | as R0 | as R0 | thiolate x the osone pool: rises unless the pool falls with pH. Does not predict the measured flat, non-monotone slope by construction |
| D4 | about a decade short (B60 -2.71 / -1.99 sigma in-sample) | as R0, scaled by its free-MFT level (second order in thiol) | as C1 | higher: first order in thiol on a sub-mM pool |
| D5 | free MFT DOWN with cystine (irreversible dimer) | as R0 | as R0 | closer to the measured rise: the release returns thiol from the 54 % held as disulfide |
| D6 | no disulfide if OX is zero in Whitfield's pots, from ambient OX if not (UNVERIFIED which; the build reports it) | as R0 | as R0 | little: a fed-norfuraneol pot holds almost no DPO / TDP, so the pool is near zero. A pass needs an oxidant the NF pot makes (sec. 11 item 6) |

What the data cannot decide: C1 vs C2 rests on D1's 720/360 ratio and one ratio in one paper (D2b, whose interval
includes no effect); expected weak. C3's two parts (release, oxidant) are judged on different data (D5; D3, D4, D6);
a split verdict (one part favoured, the other not) is reported as such and does not ship the unfavoured part.

## 6. Tolerances from each measurement's own noise

| # | noise as the source states it | sigma used |
|---|---|---|
| D1 | none printed; read error 0.001 - 0.011 dex per point (thiol_submodel.yml) | reported at TWO scales: (a) the read error in quadrature (the only stated noise); (b) the s042 series noise, MFT 0.80 / FFT 1.16 dex (a misfit property, M1.0). The sign of each comparison is read at both; its strength is printed at both |
| D2, D2p | Cao: printed SD (n = 3) where drawn, plus read error; Zhang 2025b: n = 3, SD; identical Cys / Xyl calibration curves (CF-063) | the row's SD or read error in dex |
| D2b | the atlas cell's intervals (0.035 - 0.053; 0.030 - 0.038 /min); their type UNVERIFIED (dossier "cluster 4 sec. 2") | half-width of each interval as one sd, in quadrature on the log ratio (~0.07 dex) |
| D3, D4 Zhou 2023 | mean +- SD, n = 3 per cell (Table 1) | SD propagated to the log share; ~0.02 - 0.03 dex per share, ~0.04 dex per pH difference |
| D4 Zhou 2022, Wang 2026b | printed SD per point | propagated to log(MFTD / MFT) |
| D4, D6 Whitfield 1999 | duplicates, both printed | half the duplicate range, in dex, per arm |
| D5 Zhang 2024 | triplicate error bars, figure read-off | half the drawn bar plus the read error |

**Across laboratories (H0).** Every term is within one laboratory, but the constants are shared across them. Terms
whose scale does not cancel (D2p levels, D4 and D6 share levels: a ratio of two species calibrated apart) add the
between-laboratory level spread: sigma^2 = s^2 + 2 tau_ABS(all)^2 for a two-species ratio, s^2 + tau_ABS(all)^2 for a
level. **EX-059, pending; this prereg reads its tau_ABS(all) when it exists, else the measurement noise alone,
flagged** (anti-conservative: the shares' cross-species response factors then count as exact). Ratios inside one
series or one study (D1, D2b, D3, D5's free-MFT ratio, D6's arm ratio) carry no tau: the scale cancels.

## 7. The comparison and the decision rule

**One measurement model.** Each decisive datum i is a log ratio or a log level in dex. For rival k, r_ik = predicted
minus measured, and the term is log N(r_ik; 0, sigma_i^2) with sigma_i from sec. 6. A series of time points from one
pot shares one level offset N(0, tau_ABS^2), integrated out in closed form (H0's form); with EX-059 pending, tau = 0
and the series' points count as independent (flagged). Every rival is scored on the same data with the same sigmas,
at the careful numerical setting.

**What is summed.** The predictive log-likelihood L_k = the sum over the "out" data of sec. 5 (D1, D2p, D2b, D3, the
"out" part of D4, D5, D6). In-sample terms (D2, the "in" part of D4, D5's share ratio, the 52 rows) are reported per
rival with the AIC correction for added coordinates, and are not summed into L_k. R0 is scored as shipped (B60
vector, default configuration). Its D2 failure is M1.0's structural finding.

**Uncertainty of a difference.** For each pair, the sd of Delta L across 64 draws from each fit's Laplace
covariance (draws outside a bound are reflected). This is the comparison's own noise, not a chosen number.

**The decision rule (sec. 7 of the plan, applied).**
1. Score R0, C1, C2 and C3 on the predictive decisive data under the one measurement model above; every combination
   that converged is scored, together.
2. The favoured rival has the highest L. A margin smaller than one sd of its Delta L against a rival is a tie, broken
   toward fewer fitted coordinates (R0 < C1 < C2 < C3).
3. It ships if it is favoured over every rival including R0, its fit converged, the side-effect record (sec. 8) is
   complete with a candidate compensating error registered for each loss, and no review objection (sec. 9) is left
   unanswered.

**Also printed, per D and per sigma scale:** each rival's terms, so a rival that wins the sum while losing one decisive
datum is visible. If the favoured rival loses a decisive datum of a higher evidence rank (plan sec. 7 item 2:
measured intermediates D2 / D2p above one-factor-in-one-study D2b / D3 / D5 / D6 above fitted time courses D1) by more
than that datum's propagated sd, that loss is put to the results reviewer as an objection to be answered. There is no
pass mark.

**What "ships" means here.** M1's decision admits one combination as the sulfur lane's supply-and-loss structure. Its
constants stay provisional until M2's refit, because the plan joins M1 and M2 in one refit (same rows; CF-107). If R0
is favoured, no switch changes and CF-107 stays open.

## 8. Side effects

Before (R0 = B60) and after (each converged combination), every row:

- **The full held-out contrast benchmark**, with the existing scorer `scripts/generators/score_contrast_benchmark.py
  --score`. It scores the shipped engine only today; the build adds an option to score a named Config and vector
  (UNVERIFIED that no other change is needed). Coverage first, then every answerable item: the 65 v2 TEST atlas cells
  and the 15 keyed held-out tables. R0's record is the existing
  [kinetic_core_contrast_benchmark_score_v1.md](kinetic_core_contrast_benchmark_score_v1.md).
- **The in-sample rows:** B60's 52 rows, per row, at each optimum vs B60.
- **E10 and E12 TEST contrasts,** read here and only here: F04#26 (Mizutani 1960; E10), F28#8 (Ehrenberg 1989), F28#20
  (Xu 2011), F28#25 (Wang 2022), F28#26 (Liu 2023b) (E12), where the scorer can answer them.
- **Declared hold-outs of every fit:** Yiltirak 2026 (100 - 130 C ladder), Hofmann 1998 T1 / T2 pH 3 / 7, every Zhai
  140 C pot. Not unseen (sec. 0); recorded, not used to choose.
- **Tags:** de Roos 2005 Fig. 4 (keyed E10 TEST table) is a constant source of C1 - C3 and is tagged "spent" in their
  rows (sec. 11 item 3). F28#25 is tagged "same laboratory as the choosing data" (sec. 11 item 2). The pH-flat C3
  sensitivity is printed, never used.

**Candidate compensating errors, named before the run.** A row that gets worse is expected where R0 relied on one of
these. Each loss in the record names one, or says none is known; each named one gets a CF entry at the decision.

| expected loss | why | candidate compensating error in R0 |
|---|---|---|
| Kang fed-TTCA rows (`kang_100C_MFT`, `kang_120C_MFT`, FFT, `MFT_over_FFT`; B60 +0.6 to +2.5 sigma) flip low | the route retires `r_pent_dpo`, which carried 9.4 of ~10 mM of their flux (B59b Amendment 3) | the instant cysteine-catalysed enolisation (`r_pent_dpo` / `r_pent_tdp`) standing in for the route; `k_nf_mft` ~x21 high (Cerny pot D, E6 / E11: M2) |
| `zhai_MFT_fold_120_over_100`, `zhai_FFT_fold_120_over_100` | the supply's barrier moves to `ea_form` | one barrier for every route (E11, M2); the hexose transfer of 100.9 kJ/mol |
| the 145 C fed + H2S rows (`fed_ribose_h2s_*`, `fed_tdp_h2s_FFT`) | `k_mft_decay` can no longer be fast at 145 C for free; C3's osone oxidation acts in fed TDP pots | the frozen fast sink (CF-107); the H2S headspace bias (CF-065) |
| `cerny_ternary_thiamine_share`, `fed_thiamine_MFT`; benchmark F27#15 (Grosch thiamine / ribose / Cys) | cysteine is sheltered in TTCA, so less is free to compete | the thiamine route (`k_thi_hmp`, `k_hmp_mft`) fitted under instant supply |
| `zhou_pH7_FUR_arp_alone`, `hofmann_ribose_FUR_insitu` | the pentose's fate changes | `k_pent_caramel` and `k_fur_decay` fitted against an instant pentose sink |
| the Zhou 2023 final-pH rows | TTCA formation and CAMA titration move the charge balance | the B51 balance's acid allocation (CF-084, M1b) |
| Kumazawa FFT retention | none expected: no osone in an FFT-only buffer; the release acts on FFTD | -- (a loss here is a build error until shown otherwise) |
| Yiltirak 100 - 130 C levels | the sink's T-dependence moves | laboratory and work-up offsets (CF-068, rival R4) |

## 9. The adversarial review

- **Who:** an independent agent (a fresh general-purpose subagent with no prior context, launched by the session
  that runs M1; never a fork of the author's context and never the author).
- **Two reviews.** (a) **Design**, before any build or fit: on secs. 1 - 8, 10 - 12 of this file, the case, CF-107,
  M1.0 and the cited dossiers. It does not see sec. 13 (expectations). Its task: missed direct measurements,
  undeclared transfers, rows that are not development rows, a structure built to fail or to pass. (b) **Results**,
  before the decision: the fit reports, the probe, the decisive-data and side-effect tables, the build diff and the
  data. It does not see the author's verdict paragraph or sec. 13. Its task: the strongest case AGAINST shipping the
  favoured rival, and against the rule's application.
- Every objection is answered in the record (sec. 15) or the change waits. The reviews are files beside this one
  (`..._prereg_review.md`, `..._results_review.md`).

## 10. Order of operations (blinding)

1. Design review; answer; freeze (v2 if changed).
2. Build the three keys and the scorer option, off by default; V1 - V4.
3. Feasibility probe (sec. 3). Structures that fail are out.
4. Fits of the surviving structures, sequentially (sec. 4); Laplace.
5. Only after every fit has ended: one generator run scores every rival on the decisive data
   (`scripts/generators/score_m1_decisive_data.py`, to write; artifact
   `kinetic_core_m1_decisive_comparison.{md,json}`). No decisive number is computed on an intermediate member.
6. Side effects (sec. 8).
7. Results review; answers.
8. Decision and its basis in a new EX entry (next free id, checked against origin/main and open branches); CF entries
   for the compensating errors named in sec. 8.

## 11. Differences from the case and the plan

1. **D2 does not separate C1 from C2 on mechanism.** The plan says the cysteine and sugar time courses separate C1
   from C2. C2's term acts on the thiols, not on cysteine or pentose, so on D2 the two predict the same thing. D2 is a
   requirement both must meet (S2; in-sample), and the input C2's law reads. C1 vs C2 is judged on D1's 720/360 and on
   Grosch 1992 Table VIII (added here as D2b). Expected weak (sec. 5.1).
2. **BTBU is a v2 E12 TEST laboratory.** The case marks Cao 2017, Zhang 2024 and Zhang 2025b "DEV (BTBU)"; true in v1.
   In v2, E12's TEST contrasts include F28#25 (Wang 2022), whose laboratory is `btbu_beijing`. Using BTBU rows to
   choose does not read F28#25, but F28#25 is then not laboratory-independent of the choice. It is read in sec. 8,
   tagged. (M1.0's "no TEST contrast" statement is about rows and holds.)
3. **de Roos 2005 Fig. 4 is a keyed E10 TEST table** (`heldout_deroos2005_fig4.yml`, keyed 2026-11-01: "TEST data
   for engine error E10 ... never a fit input"). C1's Cys-Amadori decay table (B59b Amendment 4, 2026-10-29) is read
   from that figure, and M1 addresses E10. It is used here as a declared constant, which spends it as a test of E10:
   registered here, tagged "spent" in sec. 8. No other measured source of this decay is on disk (Davidek 2002 is a
   glycine-Amadori analogue in another buffer). E10's remaining held-out contrast is F04#26.
4. **Zhou 2023 pH 6 / 8 columns are a declared HOLD-OUT** in FIT_HOLDOUT_DECLARATION.md ("the strongest new hold-out
   available") while both frozen splits mark Zhou 2023 DEV and a calibration input. The task and the plan name D3 as
   decisive. This prereg keeps them OUT of every objective (predictive), but choosing among rivals on them spends them.
   After M1 their status is "seen at design". The two governance files disagree.
   **Decided 2026-11-04 (session, under the owner's standing order to act; v2):** the stricter file governs. The
   pH 6 / 8 columns stay held-out: they are read only in the side-effect record (sec. 8), never to choose a rival. D3
   (the share's pH slope) is therefore not available to the comparison; C3 is judged on D4 - D6. Reason: a declared
   hold-out is spent once and cannot be recovered, while C3 is the least likely rival (sec. 13, P 0.20) and its
   other decisive data remain. The design review may challenge this.
5. **C3's oxidation is first order in thiol,** not the engine's second-order dimer form, and needs a rate-order field
   (sec. 2.4). The case left the oxidant step open (open question 1). The reason is a mechanism (rate-limiting thiolate
   addition) and G1's arithmetic; but the author has seen that the share is a decade short (sec. 0), so this choice is
   not blind.
6. **C3 as planned predicts little disulfide in Whitfield's fed-norfuraneol pots** (no pentose, so almost no DPO /
   TDP), against ~35 % of MFT measured. D6 therefore tests the pool's IDENTITY. An alternative pool that includes the
   dicarbonyls a norfuraneol pot makes (methylglyoxal, 2,3-pentanedione) is named here as a rival for after M1. It is
   not built in M1, because the plan specifies the deoxypentosones.
7. **C2's form is EX-050's cysteine factor without the electrophile pool E** (the case quotes k_E E / (1 + Cys / K_c)).
   E has no declared identity or measurement, and EX-050 / EX-055 judged the E form flexible, not identifiable (P 0.1).
8. **`k_mft_decay` is banded on Grosch 1992 Table VIII,** not on EX-050's Cys-conditional bound (sec. 2.6).
9. **Release stoichiometry (case open question 5):** Guth's measured 1.43 MFT and 1.88 FFT per disulfide lost, built as
   split integer steps; separate MFT and FFT constants instead of the one shared `k_dimer_release`.
10. **Out of M1, as the plan says:** headspace O2 (B11; case part B; open questions 3 and 6 lapse; Ehrenberg 1989 is
    not used); the pH balance (case part D, moved to M1b; held at B60); B27's source. The case's R1 (parts A - D) is
    replaced by C3. Rival R4 (laboratory offsets) is not a combination; it enters only through H0's tau.
11. **D1's noise:** the case says it uses the measurement model's per-series noise. This prereg reports both the read
    error and the s042 noise, because M1.0 found the s042 noise is a misfit property.
12. **The fit objective** is B60's 52 rows plus the route block, not R3En's 521 rows (sec. 4.1).
13. **"Ships" is provisional until M2** (sec. 7): the plan's hard dependency "M1 and M2 in one refit".
14. **The case's D1 row** says R0 "falls after ~1 h"; at B60 it peaks at 2.6 h (M1.0). Immaterial to the design.

## 12. Compute estimate

Machine: 16 GB RAM, ~6 GB usable; one fit evaluation ~13 s with 4 workers on B60's objective. A long refit is
running now, so every M1 job is `nice`d and runs alone with 4 workers (never two fits at once).

| step | count | estimate (unattended unless stated) |
|---|---|---|
| build + tests (author work) | 3 keys, a rate-order field, a scorer option, two generators | 6 - 9 h of work |
| V1 - V4 | B60 cost, anchors, unit tests | ~0.5 h |
| feasibility probe | 256 draws x 3 structures x 2 pots, quick setting, + local searches | 0.5 - 1.5 h |
| one fit evaluation, M1 objective | 52 rows' systems + ~35 route pots | ~20 - 25 s (UNVERIFIED; measured in V4) |
| one start | ~15 iterations x (19 - 21 + 1) evaluations ~ 300 - 330 evaluations | ~2 h typical; cap 1 200 evaluations ~8 h |
| fits | 3 combinations x 2 - 3 starts (8 starts), ~40 % needing one continuation | ~16 - 20 h typical; ~40 h if most hit the cap |
| Laplace | 3 x ~42 evaluations | ~1 h |
| decisive data | 4 rivals, careful setting | < 0.5 h |
| side effects | the contrast benchmark x 3 configurations, in-sample rows | ~1 - 2 h (UNVERIFIED: the scorer's run time) |
| reviews | 2 independent reviews + answers | 3 - 4 h of work |
| **total** | | **~10 - 14 h of work; ~20 - 25 h unattended compute typical (~45 h worst); about 3 calendar days** |

Memory: one evaluation with 4 workers is expected to fit in the ~6 GB usable; per-worker peak is UNVERIFIED and
recorded in V4. If it does not fit, 3 workers (~+30 % wall time).

## 13. Expectations (written now; not blind, sec. 0)

| event | P |
|---|---|
| C1 passes the feasibility probe (S1 and S2 in one draw) | 0.75 |
| C2 passes, given C1 passes / C3 passes, given C1 passes | 0.95 / 0.90 |
| every surviving combination converges within budget + one continuation | 0.6 |
| R0 favoured over every combination on the predictive data | 0.10 |
| C1 favoured (outright or by the tie rule) | 0.40 |
| C2 favoured | 0.30 |
| C3 favoured | 0.20 |
| D1: every converged C improves MFT 360/30 on R0's -1.09 dex | 0.80 |
| D1: some C makes MFT rise from 6 to 12 h (720/360 > 1) | 0.45 |
| D2b separates C1 from C2 by more than its own sd | 0.15 |
| D3: some rival reproduces both pH differences of the share within their noise | 0.15 |
| D4: C3 brings the Zhou pH 7 share (in-sample) from -2.71 sigma to within +-1 sigma | 0.35 |
| D5: C3's free-MFT ratio cystine / cysteine >= 1 | 0.30 |
| D6: C3 makes at least a tenth of Whitfield's MFT disulfide share | 0.15 |
| C3's two parts get a split verdict (release favoured, oxidant not, or the reverse) | 0.35 |
| the favoured combination's lambda_mft(100 C) stays within one prior sd of 0.044 /min | 0.65 |
| the favoured combination's 52-row cost is below B60's 27.13 | 0.30 |
| at least one held-out benchmark row gets worse by more than its noise | 0.85 |
| D2p: TTCA within stage 1's sigma (0.20 dex) at the favoured optimum / CAMA | 0.60 / 0.30 |
| M1 reaches a decision by 11-20 (plan sec. 8) | 0.40 |

## 14. What this does not do

No oxidant from headspace O2; no pH-balance change; no MFT route shares (M2); no posterior (W4). It does not test the
Yiltirak vs Schieberle level gap (not separable on disk). It does not settle the Cao vs Zhai Amadori-formation tension
(phosphate or laboratory; stage 1); `ea_form` absorbs part of it, which the results review is asked to weigh.

## 15. Results

*(appended after the run; nothing above this line is edited after the first M1 number exists)*

## 16. Revision log

- v1, 2026-11-04: first version.
- v2, 2026-11-04 (before the design review and any build): sec. 11 item 4 decided; the Zhou 2023 pH 6 / 8 columns
  stay held-out and D3 leaves the comparison.

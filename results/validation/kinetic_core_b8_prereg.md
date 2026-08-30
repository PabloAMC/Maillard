# Build Wave B8 — pre-registration

**Written and saved to disk BEFORE any B8 fit member ran, before any scorer was
re-run, and before any B8 number existed.** Scope: Amendments 16 and 17 of
`docs/reference/FIT_HOLDOUT_DECLARATION.md`. B8 is the FINAL PARAMETER WAVE of
the kinetic-core rebuild.

Sources of record read by this wave, all read-only:
`data/lit/extraction_dossiers/k6a_sulfur_ladders_synthesis.md`,
`k6b_adduct_kinetics_synthesis.md`, `zhang2026_extraction.md`,
`gigl2021_extraction.md`, `zhai2023foodchem_extraction.md`,
`feng2022_extraction.md`, `results/validation/prefactor_audit.md`,
`results/validation/kinetic_core_b2_2_diagnosis.md`,
`results/validation/cutover_final_exam.md`.

**FIREWALL, restated.** `data/benchmarks/external_validation/` is opened by
`scripts/generators/generate_cutover_final_exam.py` and by nothing else in this
wave. No B8 fit-side file names a hold-out row or carries a hold-out literal;
`tests/unit/test_kinetic_core_b8.py` asserts it by literal grep and by a
SYSTEMS-walk over every condition the fit integrates.

---

## §0. EXPOSURE DISCLOSURE — read this before any score

This wave was instructed to read `cutover_final_exam.md` (which prints the
MEASURED value of all 40 exam points) and the K6a/K6b syntheses (which print
every value of the three hold-outs B8 adds). Under the Amendment 9 clause 1 /
Amendment 10 clause 1 / Amendment 14 clause 1 precedent:

1. **The three replacement hold-outs enter as `seen_diagnostic` and may NEVER
   gate.** They are scored, printed, and pre-registered with quantitative
   predictions below — but the gating denominator does not include them.
2. **The two retired rows are therefore NOT replaced in the gating total.** The
   gating denominator FALLS by 2. Both bases are printed side by side in every
   B8 artifact so the retirement cannot read as score laundering.
3. The structural mitigation is that **no B8 parameter is fitted to any of the
   three new hold-outs**, and that the four new barriers B8 ingests are
   MEASUREMENTS the optimiser cannot move.

---

## §1. WHAT CHANGES — the complete list, declared before it is measured

### 1a. The sulfur T-structure (Amendment 17 clauses 4 and 5)

| # | change | from | to | basis |
|---|---|---|---|---|
| T1 | `Ea_decay_thiol_sink` search band | (20, 250), landed 248.0 / 216.1 / 212.9 | **(7, 102)**, prior centre **60.2** | Gigl 2021 covalent capture, `R ln 67.2 / (1/279 − 1/333)`; 3-point fit 58.5, R² 0.886; K6a §5.2. The 248 is REFUTED: it predicts k(333)/k(279) = 3.4e7 against a measured 67.2 |
| T2 | `k_dimer_mft`, `k_dimer_fft` | **no Ea at all** (held at 145 °C, policy 2) | **Ea = 122.2 kJ/mol, MEASURED, immovable** | Zhang 2026 k17 (Cys-Gly → di-(Cys-Gly) oxidation), R² 0.971, 130–150 °C, legs 87.1/159.7 |
| T3 | `k_arp_dpo`, `k_arp_tdp` | lumped formation Ea (64.08) | **Ea = 85.7 kJ/mol, MEASURED, immovable** | Zhang 2026 k16 (Cys-Amadori → α-dicarbonyl), R² **1.000**, legs 86.7/84.7 — the flat-leg member of a family whose other four members are plateau artefacts |
| T4 | Zhang k15 (Cys + sugar → Cys-Amadori, 100.9) | — | **NO SHIPPED SITE** | The sulfur network has no `Cys + sugar → ARP` edge; Zhou 2023 charges the ARP directly. Recorded as a declared prior with no site rather than forced onto a step that is not it |
| T5 | π-stacking reservoir (ΔH° = −19.5 kJ/mol, reversible, < 5 min) | — | **DECLARED MODEL-FORM GAP, not implemented** | Gigl §6e. It is a reversible negative-ΔH sequestration term the module does not have. Adding a species after reading a hold-out would be a refit in disguise |

**Why T1's band is a PRIOR band and not a bound trade.** (7, 102) is the
defensible range Gigl's six derivation routes span, not a numerical convenience.
The band's ceiling being *below* B2.3's landed 216.1 is the point: the previous
value is refuted by measurement, and a fit that presses the new ceiling is
reporting a residual conflict rather than a barrier.

**Why T2 is not the prohibited derivation.** Policy 2 (`parameters_sulfur.py`,
B2.1) forbids drawing ONE Arrhenius line through TWO of four named channels
measured in different papers at different temperatures. Zhang k17 is a SINGLE
paper's SINGLE step measured at THREE temperatures. It is the same class of
object as Kang/Zhai's 55.1 on `k_cys_thermal`, which `MEASURED_EA_OVERRIDES`
already carries. `k_mmft`, `k_thioether`, `k_oligomer` and `k_protein_ss` keep
no activation energy and the prohibition on pairing them stands.

### 1b. Provenance corrections (Amendment 17 clauses 1 and 2)

* Every `kang_*` anchor string re-points to **Zhai et al., *Food Chem.* 404:134420
  (online 28 Sep 2022)** as the primary source, with Kang 2026 recorded as a
  re-publication of the same 102 cells (101/102 identical, including a
  reproduced arithmetic error).
* **The Tier A ±15 % label is WITHDRAWN.** The primary source prints
  `Wi = f′ × (Ai·ms)/(As·V)` with `f′ = 1`: single-IS semi-quant, n = 1 pot.
  The re-band is Tier B, × ÷ 2–3 on magnitude.
  **DECLARED CONSEQUENCE: this costs no value change.** The six shipped Kang
  level rows already carry `sigma_log = 0.4` = log10(2.51), which is inside the
  × ÷ 2–3 interval. The re-band is therefore a LABEL correction, and saying so
  is more honest than widening a sigma to manufacture a visible edit.
* `kang_switch_on_MFT` and `kang_switch_on_FFT` are **RETIRED** from the panel:
  the target is a hold-time artefact of a single 120-min slice through curves
  that peak and turn over (four axes, three of them inside the source lab).
* 2-methylthiophene at 120 °C is `not detected`, not 1.128 (copy-paste in both
  Zhai papers; Kang's 0.000 is the correction). **No shipped row uses it.**

### 1c. Four new declared FIT rows (Amendment 17 clause 5)

All four are RATIO- or CONVERSION-scale, i.e. legal for a semi-quant source.
None is a level. None touches a 140 °C rung.

| id | kind | target | source |
|---|---|---|---|
| `feng_ARP_conversion_100C` | conversion of `ARP` | **0.48** (σ_log 0.30) | Feng 2022 Fig. 2a, RG-ARP 20 mmol/L, pH 7 unbuffered, 120 min, 100 °C. Digitised endpoint — coarse, hence the wide sigma |
| `feng_ARP_conversion_120C` | conversion of `ARP` | **0.935** (σ_log 0.15) | Feng 2022 Fig. 2a, same pot, 120 °C: 93.50 % consumed. Printed in the text |
| `zhai_MFT_fold_120_over_100` | cross-system ratio | **1.12** (σ_log 0.15) | Zhai *Food Chem.* Table 1, TTCA arm, 120 min: MFT 1.237 → 1.388 µg/L |
| `zhai_FFT_fold_120_over_100` | cross-system ratio | **1.10** (σ_log 0.15) | same, FFT 3.734 → 4.107 µg/L |

The two Zhai rows are the response-factor-free form of information the two
Kang level pairs already carry at σ_log 0.4 each; entering the RATIO at a tight
sigma and leaving the LEVELS at their semi-quant width is exactly what
"shapes and ratios, never absolutes" means operationally.

The Feng rows are new information: nothing in the objective previously measured
ARP consumption at two temperatures, and T3 gives ARP consumption a measured
barrier for the first time.

**Objective size: 58 → 62 declared FIT rows.**

### 1d. Prefactor fixes (Amendment 16 clause 1)

* `cysteine_thermolysis` in `data/lit/arrhenius_params.yml` (the Cantera lane):
  **A = 1e14 → 1.931e12 s⁻¹**, and **Ea 130.4 → 133.0 kJ/mol** with it. The Ea
  moves because a prefactor cannot be transplanted across a different barrier:
  1.931e12 is the refit of Zheng & Ho 1994 Table I's **pH 5.0 column**, whose
  paired Ea is 133.0. Keeping 130.4 while adopting the paired A would introduce
  a fresh ~2× rate error at 150 °C in the name of fixing a 51.8× one. The
  shipped 130.4 was the pH 3–9 MEAN, which is not a matched pair with any
  single column's prefactor.
  This also CLOSES the cross-lane flag: `k_cys_h2s` in
  `parameters_sulfur.py` already ships (133.0, 1.158e14 min⁻¹ = 1.93e12 s⁻¹).
  The two lanes now agree to 0.1 %.
* **Nothing else moves.** Per the Amendment 16 ruling, the six Kocadagli-lane
  flags stay diagnostics-only (the shipped values come from the source's
  simultaneous global ODE fit, not from a transcribed k(T) table; the artifact
  says so itself); `k_thioether`'s 2.3e5× is by design (Ea only transferred);
  and `schiff_condensation` / `k_schiff` stay flagged because NO source k(T)
  table exists for either, so the 15.9× disagreement is undecidable and
  choosing a side would be inventing a number.

### 1e. Covalent-sink retirement (Amendment 17 clause 6)

The state of B4's covalent ceiling changes from
*"unmeasured, could matter at process temperature if Ea ≥ 70 kJ/mol"*
to *"MEASURED NEGLIGIBLE: Ea = 15–23 kJ/mol"*, on three independent legs:

1. **Rate.** Shepelev 2024 Fig. 3, ¹⁴C on β-lactoglobulin, 25/45/65 °C:
   Ea = 15.2 (raw) / 20.0 (saturation-corrected), envelope 14–23, corroborated
   at 10.0 by Hamzalioglu's HMF + lysine in a different lab by a different
   method. The channel removes **0.005 %–0.21 %** of a hexanal dose in any real
   thermal process, against the 1 %–91 % the ≥ 70 assumption implied.
2. **Capacity falls with temperature.** Day-28 binding is 25.6 > 20.1 > 17.4
   mg/g at 25 / 45 / 65 °C; the 65 °C series peaks at day 7 and declines 26 %.
3. **The equilibrium unwinds on heating.** Zamora 2010's matched forward/reverse
   pair (44 fwd, 52 rev) gives ΔH ≈ −8 kJ/mol, K_eq falling 3.0× from 25 to
   180 °C, demonstrated by a 15× recovery of "vanished" acrylamide.

**DECLARED, BEFORE THE EXAM RUNS: this changes no number.** The channel was
already INERT BY RULING and contributed exactly 0.0 to every point prediction.
The retirement replaces an *unmeasured* zero with a *measured* zero and rewrites
what the reports say about it. See prediction X-5.

`k_protein_ss` — the thiol/protein-**disulfide** exchange channel — keeps its
bounded carry unchanged. It is a different reaction (disulfide exchange, not
carbonyl–amine addition), a different lane, and its flux is identically zero in
every scored row.

### 1f. Milk thresholds (Amendment 17 clause 6)

`?/kg = µg/kg`, closed by arithmetic: seven of seven `tian2020` Table 1
concentrations reproduce `tian2019`'s Y6 column digit for digit, and that
column is headed `Concentration (μg/kg)` in plain type.
**The `milk_tian` seal is lifted; NO VALUE CHANGES**, because no milk threshold
was ever tabulated in `src/` — the seal blocked the lookup, and what is edited
is the seal's reason string and the three caveats that must now travel with the
rows (`same_matrix: FALSE`; the reference thresholds are unsourced and
mislabelled "in air"; the printed SDs are dilution-series steps, not threshold
uncertainty).

### 1g. `safety.py`'s four uncited pairs (Amendment 16 clause 2)

`safety_dicarbonyl_pool_{formation,elimination}` and
`safety_furosine_{formation,elimination}` are LABELLED `no_verifiable_source`
per the Wave T3 convention — a module-level provenance record plus a
`RuntimeWarning` at first use naming the defect, what depends on it, and the
fact that **no value is substituted or rescaled.** Four numbers, zero moves.

---

## §2. THE FIT — declared before it runs

* **Weighting: W-HALF only**, `E = 0.70` decades per pH unit
  (`σ_ph = 0.35 / 0.70 = 0.50`), the value B2.4's own exam-blind and
  panel-blind criterion shipped. B8 does not re-open the weighting question.
* **Starts: 2, seeded.** Start 0 is B2.4-half's ensemble-best vector (the
  incumbent, so a B8 member worse than its predecessor is visibly worse);
  start 1 draws uniformly inside a declared local neighbourhood of it
  (± 2 decades on a log10 k, ± 40 kJ/mol on a barrier — clipped to the new
  bands), seed `20260830 + 1000`.
* **Free set: 23 of 48.** B2.4's 20 (clauses R1–R4) plus three coordinates the
  T-structure moves out from under:
  `k_dimer_mft` (gains a barrier where it had none),
  `k_arp_dpo` and `k_arp_tdp` (barrier 64.08 → 85.7).
  Freezing a k_ref whose Ea just changed would report the old rate at the new
  barrier, which is neither the old model nor the new one.
  `k_dimer_fft` and the other 19 are already in B2.4's free set.
* **Budget:** 600 evaluations per member, same exhaustion rule as B2.4 (a
  budget-exhausted member returns best-so-far with status −9 and is excluded
  from any spread statistic). ~2.35 s/evaluation measured on this container.
* **Comparability.** B8's objective has four more rows than B2.4's, so total
  cost is NOT comparable. The comparator reported for every member is
  `sum_r2_level_shared_55`: the sum of squared residuals over exactly the 55
  non-pH rows B2.4 scored, at their unchanged σ_log.

**PROMOTION, DECLARED BLIND.** B8's frozen fit report becomes the engine's
FIRST parameter candidate, ahead of B2.3's. This is declared here,
unconditionally, before any B8 score exists, and it is NOT contingent on any
scorecard: what B8 carries is four measured barriers and the removal of a
refuted one, and shipping a refuted 248 because its successor scores worse
would be preferring a known defect. If B8 scores worse, it ships anyway and the
regression is reported as the price of the measurement.

---

## §3. PRE-REGISTERED PREDICTIONS

### The fit

| # | claim | confidence |
|---|---|---|
| **F-1** | **`Ea_decay_thiol_sink` lands ON its new 102 ceiling** (within 2 kJ/mol). B2.2/B2.3/B2.4 all wanted 213–248; a band whose ceiling is 2.1× below the incumbent will most likely be pressed. **This is a prediction that the wave's own headline parameter comes back bound-limited again, and it is stated first for that reason.** | **60 %** |
| **F-2** | `sum_r2_level_shared_55` **WORSENS** against B2.4-half's 15.123 — four barriers are now immovable measurements and the thiol band is 2.5× narrower. Predicted range **15.5 – 25.0**. | 70 % (that it worsens); 55 % (that it lands in the range) |
| **F-3** | The lumped formation Ea stays inside **40 – 100 kJ/mol** (it was 64.08). | 75 % |
| **F-4** | The two Feng ARP conversion rows are reproduced within **3×** at 100 °C and **2×** at 120 °C. The 120 °C row is the harder one: 93.5 % consumed is close to exhaustion and the model must not overshoot past it. | 55 % |
| **F-5** | Giving the dimer channels a 122.2 kJ/mol barrier makes them run **much slower below 145 °C** than the previously-frozen 145 °C value, so the thiol pool at Kumazawa's 121 °C grows unless something else takes it. Predict **at least one of `k_dimer_fft`, `k_thiolate_loss` moves ≥ 0.5 decades** from its B2.4-half value. | 70 % |
| **F-6** | Still **≤ 5 of 48** parameters carry a finite Gauss-Newton interval. Nothing in B8 adds identifiability; it removes three fitted degrees of freedom's worth of freedom and adds four rows. | 80 % |

### Blind re-sit (a) — the full hold-out panel

| # | claim | confidence |
|---|---|---|
| **H-1** | On the **OLD basis** (the same 34-row gating denominator B2.4 used, with the two retired rows scored and counted), B8's gating total is within **±2** of B2.4-half's. | 65 % |
| **H-2** | Had they been scored, **both retired rows would still have FAILED.** They are reported scored-but-not-counted so that the retirement is visibly not a rescue. | 85 % |
| **H-3** | ★ **`hofmann2002_brew_80C_FFT` and `vanseeventer_50C_MFT_per_day` are AT RISK, and this is the direct price of T1.** A thiol-sink barrier falling from 216 to ≤ 102 with `k_ref` anchored at 145 °C makes the sink **relatively faster when cold** — which is the exact mechanism that destroyed both rows in B2 and that B2.1 fixed by giving the sinks a barrier at all. Predict **at least one regresses**. | 60 % |
| **H-4** | `kang_140C_FFT` **stays failing** (it was 0.0487×). | 80 % |
| **H-5** | The three `zhou_120C_dimer_share_pH*` rows **move**, in either direction, by ≥ 1.5× — they are the only rows that directly score the dimer branch T2 re-barriers. | 70 % |

### Blind re-sit (b) — the both-ways cutover exam

| # | claim | confidence |
|---|---|---|
| **X-1** | ★★ **THE HEADLINE. The Yiltirak family improves.** It is 0/8 at a family median of **193.5×**, OVER-predicting at every rung and worst at the long holds, and `cutover_final_exam.md` names the time axis. K6a supplies the missing mechanism three times over (Zhai's FFT peaking at 80 min at 140 °C and falling 2.5×; Wang's collapse over 105→125 °C while furfural saturates; Zhang's supply saturation above 140 °C): **a core whose sink barely runs cannot produce a peak, so it must over-predict every long hold.** T1 is the named fix. Predict the family median improves by **≥ 3×**, to **≤ 65×**. | **55 %** |
| **X-2** | Yiltirak enters the band on **at least 1 of 8** points. | 35 % |
| **X-3** | The exam's core paired median improves against **42.23×**. | 60 % |
| **X-4** | The `sulfur_hofmann1998_145C` family median moves by **< 1.5×** in either direction. Its rows sit at 145 °C, which is the `k_ref` reference temperature, where a barrier change is nearly inert by construction. | 70 % |
| **X-5** | ★ **The `acrylamide_180C`, `furan_browning_glc_alanine` and `matrix_path_lipid` families are BIT-IDENTICAL.** No acrylamide, furanic or lipid parameter is touched, and 1e's covalent retirement is documentation over a term that already contributed exactly 0.0. **If any of these three moves, the covalent retirement was not the no-op this wave claims it is, and that is a finding, not a rounding error.** | 90 % |
| **X-6** | The exam's 6 declined points stay declined and the answered count stays **34**. | 90 % |

### Blind re-sit (c) — the three replacement hold-outs

All three are `seen_diagnostic` and none gates. Predictions are made anyway,
because a hold-out that carries no prediction teaches nothing.

| # | row | claim | confidence |
|---|---|---|---|
| **N-1** | `wang2026_MFT_peak_and_fall_125_over_115` — measured **0.30×** (MFT peaks at 115 °C and falls) | The model **predicts a RISE (> 1.0×) and FAILS.** Its thiol formation is monotone in T once the sink is anchored at 145 °C, and 115→125 °C is below every turnover the network can express. | 70 % |
| **N-2** | `wang2026_FFT_peak_and_fall_115_over_105` — measured **0.15×** | Same direction, same reason; **FAILS**, but by less than N-1 because FFT's precursor (furfural) is itself saturating in the measurement and the model does carry a furfural sink. | 70 % |
| **N-3** | `zhai_13C_exogenous_carbon_threshold` — measured 0 % → 19–21 % exogenous sugar carbon between 100 and 120 °C | **The module cannot answer it in kind.** The three compounds isotope-traced are thiophenes, all declared out of scope, and MFT/FFT were never traced. Scored on the only in-scope proxy — the co-substrate increment on the module's own retro-aldol product — and reported as a SCOPE GAP rather than a pass or a fail. | 80 % that it is reported as a scope gap |
| **N-4** | `ames2001_excess_Ea_class_split` — free thiols positive, thiazoles negative, on the 120→150 °C leg | The **sign split reproduces** (MFT and FFT excess Ea over the bulk positive, 2-acetylthiazole negative) **in the module's aqueous pot**, but the row is DIAGNOSTIC because Ames measures a low-moisture extrudate and K6a §5.1 measures a ~2× aqueous→matrix gap on exactly these two barriers with no transfer term in the model. | 45 % on the split |

---

## §4. WHAT WOULD FALSIFY THE WHOLE WAVE

1. **If X-5 fails** — if any acrylamide, furanic or lipid exam point moves —
   then the covalent retirement touched a live number and this wave's central
   claim that it is documentation is false.
2. **If F-1 holds AND X-1 fails** — a barrier pressed against 102 that does not
   buy the Yiltirak turnover means the turnover is not a barrier phenomenon at
   all, and K6a's diagnosis (§9.4) is wrong about the mechanism even though it
   is right about the 248 being refuted.
3. **If H-3 realises on BOTH rows** — the low-temperature rows dying is the
   price T1 was always going to risk. Two of two would mean the corpus cannot
   hold a single thiol-sink barrier at any value, which is the two-channel
   model-form falsification Gigl §7c predicts and T5 declines to implement.

---

## §5. WHAT THIS WAVE WILL NOT DO

* It will **not** refit B1 (the trunk), B3, B4, B6 or B7.
* It will **not** implement the π-stacking reservoir (T5). It is declared as a
  model-form gap with its measured enthalpy on the record.
* It will **not** move Kang/Zhai's 55.1 to Zhai's primary-source 54.7. They are
  the same experiment read two ways (a digitised figure vs a printed sentence)
  and the 0.7 % gap is a validation of the digitisation, not new information.
  The provenance is corrected; the number is not.
* It will **not** adjudicate the contradictions K6a and K6b report. They are
  carried forward into the B8 reports as contradictions.
* It will **not** commit.

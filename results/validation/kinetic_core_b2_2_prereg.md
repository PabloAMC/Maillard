# Build Wave B2.2 — PRE-REGISTRATION

**Written 2026-08-29, BEFORE the fit was run and BEFORE any hold-out row was
scored.** Nothing below was edited after a number arrived. The outcome column
is filled in by `kinetic_core_b2_2_holdout_report.md`, which scores every claim
here whether it held or not.

Predecessors: `kinetic_core_b2_1_fit_report.{json,md}`,
`kinetic_core_b2_1_holdout_report.md` (**15 / 33 gating**, and **15 / 27** on
B2's own row set), `kinetic_core_b2_1_diagnosis.md`,
`cutover_final_exam.md` (Yiltirak **0 / 8**, median **290.5×**; Chang: one
number for two arms).

---

## 0. What was changed, and what it cost

| change | new free parameters | where they are identified |
|---|---:|---|
| **Task 1** — two named decay-barrier families (`thiol_sink`, `carbonyl_sink`) replacing the shared lumped-formation Ea on seven of the nine `*_decay` lumps | **2** | Kumazawa 121 °C formation-free survival grid + Hofmann 145 °C pH-5 anchors + Kang's declared 100/120 °C rungs (thiol); Kang's furfural ladder + Hofmann T5 in-situ furfural/norfuraneol (carbonyl) |
| **Task 2** — dynamic pH from a charge balance against the tracked acid pool and the declared buffer | **2** (`acid_yield_per_sink_event`, `arp_secondary_ammonium_pKa`) | Amendment 7's three declared Zhou drift endpoints, and nothing else |
| **total** | **4** | 58 declared FIT rows against 48 free parameters |

`k_dimer_decay` and `k_h2s_loss` deliberately keep B2.1's setting: nothing in
the corpus measures a disulfide or a sulfide sink at two temperatures.

### Two structural side-effects that are NOT new parameters, and are declared here because both move numbers

1. **Two pH scales.** Every pH label in this corpus is a cooled-sample reading.
   B2/B2.1 used one scale and therefore evaluated temperature-corrected pKa
   against a bench pH. B2.2 uses the label only to set the conserved strong-ion
   difference at 25 °C, then solves the in-situ pH at reaction temperature.
   **Measured consequence, before the fit: an unbuffered pot labelled pH 7 runs
   at an in-situ pH of ~5.6–5.9, not 7.0.** Every Zhou, Kang and drift system
   moves by ~1.1 units. Buffered systems move by <0.2.
2. **`r_ttca_deg` now carries one acid equivalent.** TTCA is a carboxylic acid;
   a degradation route that deletes its carboxylate while leaving the NaOH that
   titrated it manufactures strong base. Measured doing exactly that in this
   wave's probe: Kang's pot drove to **pH 11.4** against a measured 4.9. With
   the carboxyl carried it lands at **4.62**.

---

## 1. PRE-REGISTERED PREDICTIONS — the fit

| # | claim | confidence |
|---|---|---|
| F-1 | The fit converges with **reduced χ² between 1.5 and 6.0**, i.e. WORSE than B2.1's 2.01, because four of the new rows/constants are fought over by systems that previously could not disagree. | 70 % |
| F-2 | The **three Zhou drift endpoints are reproduced to within 0.5 pH units each**. A standalone scan at the frozen B2.1 rate constants already reached (3.28, 3.45, 4.80) against (3.22, 3.42, 5.07), so the drift model *can* do it; the risk is that the joint fit trades it away. | 75 % |
| F-3 | **The two drift constants are CORRELATED and at most one is individually identified.** The same scan showed a ridge: (yield 0.8, pKa 11.0) and (yield 1.0, pKa 10.0) score within 6 % of each other. | 80 % |
| F-4 | The fitted **Amadori secondary-ammonium pKa lands ABOVE 9**, which is 1–2 units higher than the ~7–8 a literature N-glycosylated secondary amine would give. If it does, that is reported as a **contradiction between the calibrated value and the chemistry**, not as a measurement — the fit has an incentive to push the pKa up simply to make the Amadori electrically invisible. | 70 % |
| F-5 | **`Ea_decay_thiol_sink` lands BELOW the lumped formation Ea.** This is the direction B2.1's diagnosis asks for and the direction Kumazawa-at-121 vs Hofmann-at-145 should supply. | 55 % — genuinely uncertain, because Kang's near-flat 100→120 °C thiol ladder pulls the other way |
| F-6 | **`Ea_decay_carbonyl_sink` lands between 30 and 90 kJ/mol**, near the 32.8–43.5 the Kang dossier reads off the furfural ladder itself. | 50 % |
| F-7 | **The buffered systems do NOT stay near-constant.** The wave brief expects "near-constant pH (Hofmann 0.5 M phosphate)". They will not be, and the reason is arithmetic: 0.5 M phosphate at pH 5.0 sits between its pKa1 (2.15) and pKa2 (7.20) and has a van Slyke capacity of only **≈15.7 mmol/L per pH unit**, against a modelled acid load of **≈73 mmol/L** in the Hofmann pentose pot and **≈227 mmol/L** in the Cerny ternary. At acid yield 0.9 the probe puts Hofmann at a cooled endpoint of **3.67** and Cerny at **3.36**. **This is pre-registered as a CONTRADICTION between the brief's expectation and the model's own buffer arithmetic**, and it is the single largest risk in this wave. | 85 % that the drift exceeds 0.5 units at acid yield ≥ 0.5 |
| F-8 | Because of F-7 the fit is expected to **push the acid yield DOWN**, below the ~0.9 the Zhou anchors alone want, trading the drift anchors against the buffered volatile rows. | 60 % |

**Falsifier for the wave as a whole:** if F-2 fails *and* F-7 holds, the drift
model is fitting one set of systems by breaking another and should be reported
as not yet usable, not as shipped.

---

## 2. PRE-REGISTERED PREDICTIONS — blind re-sit (a): Kang 140 °C

The gating rows are `kang_140C_MFT`, `kang_140C_FFT`, `kang_switch_on_MFT`,
`kang_switch_on_FFT`. B2.1 scored **0 / 4**, under-predicting by 1 500× (MFT)
and 6.6 million× (FFT).

| # | claim | confidence |
|---|---|---|
| K-1 | **The 140 °C under-prediction improves by at least 100× on FFT.** B2.1's FFT collapse was 6.6 million-fold and its *mechanism* is named: the sinks share the formation barrier and eat a precursor-exhausted pot. Separating the barrier is a direct attack on that mechanism. | 60 % |
| K-2 | **The switch-on does NOT emerge and the block still scores 0 / 4.** Separating the decay barrier can stop the model destroying its own product; it cannot manufacture the *acceleration* the measurement shows, because B2.1's diagnosis §2b found the second half of the failure elsewhere — the sulfide pool peaks at 120 °C and FALLS at 140 °C, since the Nedvidek step burns cysteine before thermolysis can turn it into sulfide. Nothing in B2.2 touches that. | **70 %** |
| K-3 | If any of the four passes, it is **`kang_140C_FFT` before `kang_140C_MFT`** — FFT's measured 120→140 leg is the shallower of the two (2.78× vs 4.26×). | 60 % conditional |
| K-4 | The 140 °C **cysteine-conversion** row (diagnostic, passed at 1.38× in B2.1) **stays passing**: it is governed by Kang's own measured 55.1 kJ/mol, which no change in this wave may move. | 90 % |

---

## 3. PRE-REGISTERED PREDICTIONS — blind re-sit (b): the full B2.1 hold-out panel

Baseline: **15 / 33 gating**, 6 / 9 diagnostic. Target: zero regressions.

| # | claim | confidence |
|---|---|---|
| H-1 | **Net gating score lands between 13 and 20 of 33.** | 70 % |
| H-2 | **There WILL be regressions — between 1 and 5 rows.** "Zero regressions" is the target, not the expectation: the two-pH-scale correction moves every unbuffered system's operating pH by ~1.1 units and the buffered systems' by up to 1.3, and no re-fit recovers that for free. Pre-registering a non-zero number is the honest position. | 65 % |
| H-3 | **The Hofmann pH-3 / pH-7 aqueous rows are the most likely regressions**, because they are the rows whose *nominal* pH is now converted to a different in-situ pH AND whose buffer now drifts. | 55 % conditional on H-2 |
| H-4 | **`zhou_dimer_share_pH_invariance` (187× in B2.1, the wave's one material regression) IMPROVES.** B2.1's own §4 named the reconciliation: Zhou's three runs converge in pH, so their *effective* pH span is far smaller than their labels. B2.2 is the first version that can actually compute that convergence — the probe gives in-situ starts of 4.2 / 5.6 / 6.8 collapsing to ~3.3 / 3.4 / 4.5. | **75 %** |
| H-5 | The three `zhou_120C_dimer_share_pH{6,7,8}` rows do **not** all survive; satisfying the invariance row and the three level rows simultaneously was already flagged as an unresolved conflict between Kumazawa (buffered, strongly pH-dependent) and Zhou (unbuffered, share invariant). Expect **2 of 3** to pass. | 55 % |
| H-6 | **`zhou_pH6_MFT` (diagnostic, 0.058× in B2.1) improves by ≥ 3×.** The declaration's §5 decision 1 scores it diagnostic *only* "until the model carries a pH-trajectory state"; this is that state, and the pH-6 and pH-7 runs converging is exactly what the row needs. | 65 % |
| H-7 | **`vanseeventer_50C_MFT_per_day` and `hofmann2002_brew_80C_FFT` — the two low-temperature rows B2.1 FIXED by giving the sinks a barrier — are at risk of being un-fixed** if `Ea_decay_thiol_sink` lands low. This is the direct cost of F-5 and it is stated in advance. | 45 % that at least one regresses |
| H-8 | `cerny2007_branch_responsiveness` **still fails**. B2.1's §3 proved the ceiling is arithmetic (the model's 1× xylose share of 0.412 caps responsiveness at 2.43× against a measured 3.067×) and nothing in B2.2 changes a branch share. | **90 %** |

---

## 4. PRE-REGISTERED PREDICTIONS — blind re-sit (c): the cutover final exam

Re-run of `generate_cutover_final_exam.py`, pure re-scoring, no scorer change.
Baseline: core **5 / 23** answered points inside 3.0×, median **24.93×**;
Yiltirak family **0 / 8**, median **290.5×**; Chang two arms → one number.

| # | claim | confidence |
|---|---|---|
| X-1 | **The Yiltirak ladder improves substantially — median fold error falls below 100×** (from 290.5×). Its named defect is the time axis: the core over-predicts at all four rungs because it accumulates thiol over holds of 30 min to 4 h. A thiol sink that no longer shares the formation barrier is precisely the lever on a long low-temperature hold. | 55 % |
| X-2 | **Yiltirak still scores 0 / 8.** Getting from 290× to inside 3× is a two-order-of-magnitude move and the band is unforgiving. | **70 %** |
| X-3 | **The Chang two-arm defect is NOT fixed and the core still returns one number for both arms.** The arms differ by 1 % acetic acid vs deionised water, and they are on the **acrylamide** lane, which carries no pH term at all. B2.2's pH state lives in the sulfur lane only. This is pre-registered as an EXPECTED NON-FIX so it cannot later be presented as a surprise. | **95 %** |
| X-4 | The `sulfur_hofmann1998_145C` family (the core's best result anywhere, 4 / 10) **loses at least one row**, as the direct cost of F-7's buffered drift. | 55 % |
| X-5 | **The overall paired median does not beat the old lane** (core 24.93× vs old 12.65×). Two barriers and a pH state are not the missing lipid lane. | 85 % |
| X-6 | The 17 declensions stay 17. Nothing in B2.2 adds a species to the acrylamide or matrix lanes. | 95 % |

---

## 5. What would make this wave a failure, stated in advance

1. Gating score below 13 / 33 (a net regression against B2.1).
2. The Zhou drift anchors missed by more than 1.0 pH unit on any of the three —
   the drift model would then be a parameter that fits nothing.
3. Either decay barrier landing **on a bound** (20 or 250 kJ/mol), which would
   mean the family is not identified and the parameter should be withdrawn
   rather than reported.
4. Any measured parameter moving. `MEASURED_EA_OVERRIDES` (Kang's 55.1 kJ/mol),
   Zheng & Ho's four-pH thermolysis set, Stack's equilibrium, and every
   textbook pKa in `ph_state.py` are immovable and unit-tested as such.

## 6. Firewall disclosure, stated before the scores

`kang2026_SI_extraction.md` prints the 140 °C column (MFT 5.907, FFT 11.439,
sulfur subtotal 60.400) and this wave was directed to read it; the same dossier
is where the declared 100/120 °C FIT rows live. `cutover_final_exam.md` prints
every Yiltirak measured value and this wave was directed to read that too.
**Those hold-out values were therefore SEEN before this fit ran.** They enter no
row, no bound and no initialisation; `tests/unit/test_kinetic_core_b2_2.py`
enforces that by a literal grep over the fit-side files and by a walk over the
fit script's own `SYSTEMS` dictionary. The frozen bundles under
`data/benchmarks/external_validation/` were opened by nothing except the exam
generator's own pre-existing scoring path.

**One new exposure this wave creates and declares:** the fit now integrates the
Zhou pH-6 and pH-8 *systems*, which are held-out **volatile** columns. Amendment
7 licenses this for the **endpoint pH only**. The mechanical guarantee is a
per-row assertion: every fit row whose system name starts with `zhou_drift_`
must have `kind == "ph_endpoint"`. No volatile level from those columns can
reach the objective.

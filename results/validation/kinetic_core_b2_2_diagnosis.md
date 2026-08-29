# Build Wave B2.2 — what the numbers said, AFTER the scorecards were frozen

**Everything in this file was written after looking at the hold-out and exam
scores. NO PARAMETER WAS CHANGED, before or after.** The frozen fit is
`kinetic_core_b2_2_fit_report.json`, the scorecard is
`kinetic_core_b2_2_holdout_report.md`, the exam is `cutover_final_exam.md`, and
the pre-registration is `kinetic_core_b2_2_prereg.md` — written and committed
to disk before the fit ran. This file exists so that a post-hoc explanation is
filed as a post-hoc explanation rather than folded back into a report where it
would read like a prediction.

The diagnostics below were run with `scratch/b22_kang140_diag.py` and
`scratch/b22_buffer_sensitivity.py` against the frozen parameters.

---

## 1. The headline

| | B2.1 | **B2.2** |
|---|---:|---:|
| fit reduced chi-square | 2.01 | **1.80** |
| fit rows / free parameters | 51 / 44 | **58 / 48** |
| gating hold-out | **15 / 33** | **12 / 33** |
| diagnostic hold-out | 6 / 9 | 6 / 9 |
| rows FIXED / REGRESSED vs predecessor | 11 / 0 | **6 / 9** |
| Kang 140 °C block (gating) | 0 / 4 | **1 / 4** |
| Kang 140 °C FFT fold error | 1.5 × 10⁻⁷ | **0.0706** |
| exam, core paired median fold | 24.93× | **10.65×** |
| exam, old-lane paired median | 12.65× | 12.65× |
| exam, Hofmann 145 °C family median | 7.09× | **4.16×** |
| exam, Yiltirak family | 0 / 8, median 290.5× | 0 / 8, median **262.3×** |
| Zhou drift endpoints (3 declared FIT anchors) | not modelled | **3.29 / 3.50 / 5.07** vs 3.22 / 3.42 / 5.07 |

**The wave traded three gating hold-out rows for the thing it was sent to get.**
The Kang 140 °C collapse — B2.1's worst single failure, six and a half orders of
magnitude on FFT — is essentially gone, and the exam's paired median improved
2.3× to the point where the core now beats the old lane on the only
apples-to-apples number in that report. Against that, nine rows regressed and
the gating total fell from 15 to 12. **Neither half of that should be reported
without the other.**

---

## 2. Pre-registration scorecard — every claim, held or not

### The fit

| # | claim | outcome |
|---|---|---|
| F-1 | χ²ᵣ between 1.5 and 6.0, and **worse** than B2.1's 2.01 | **HALF-FALSIFIED.** 1.80 is inside the band but it is *better*, not worse. |
| F-2 | the three Zhou drift endpoints reproduced within 0.5 pH units | **HELD**, with room to spare: worst error **0.08 units** (pH-7 anchor). |
| F-3 | the two drift constants are correlated; at most one individually identified | **HELD.** Neither carries a finite interval; **3 of 48** parameters are identified at all. |
| F-4 | the Amadori ammonium pKa lands above 9 and must be reported as a contradiction | **HELD.** 10.77 — 3 units above what an N-glycosylated secondary amine should give. See §5. |
| F-5 | `Ea_decay_thiol_sink` lands **below** the formation barrier | **FALSIFIED, and this is the wave's most informative negative.** It landed at **248.0** against a formation barrier of **76.2**, i.e. pressed against its 250 ceiling. |
| F-6 | `Ea_decay_carbonyl_sink` between 30 and 90 kJ/mol | **FALSIFIED.** 150.2. |
| F-7 | the buffered systems do **not** stay near-constant | **HELD.** Hofmann's declared 0.5 M phosphate drifts 5.00 → **3.59** (cooled), because its van Slyke capacity at pH 5 is only ≈15.7 mmol L⁻¹ per pH unit against a modelled 86.7 mmol L⁻¹ acid load. |
| F-8 | the fit pushes the acid yield down, below the ~0.9 the Zhou anchors alone want | **FALSIFIED.** It went to **0.968**, i.e. the fit chose to accept the buffered drift rather than suppress the acid. |

### Blind re-sit (a) — Kang 140 °C

| # | claim | outcome |
|---|---|---|
| K-1 | the 140 °C under-prediction improves by ≥ 100× on FFT | **HELD, by 4.6 × 10⁵×** (1.52e-07 → 0.0706). |
| K-2 | the switch-on does not emerge; the block still scores 0 / 4 | **HALF-FALSIFIED.** `kang_140C_MFT` PASSES at 2.15×, so the block is 1 / 4 — but both `kang_switch_on_*` rows still fail, so the *switch-on itself* did not emerge, which was the substance of the claim. **And the pass should not be counted — see §3.** |
| K-3 | if anything passes it is FFT before MFT | **FALSIFIED.** MFT passed, FFT did not. |
| K-4 | the 140 °C cysteine-conversion row stays passing | **HELD** (1.24×), governed by Kang's own immovable 55.1 kJ/mol. |

### Blind re-sit (b) — the full B2.1 hold-out panel

| # | claim | outcome |
|---|---|---|
| H-1 | gating between 13 and 20 of 33 | **FALSIFIED.** 12. |
| H-2 | between 1 and 5 regressions | **FALSIFIED.** Nine. |
| H-3 | the Hofmann pH-3 / pH-7 rows are the most likely regressions | **HELD for pH 7** (both MFT and FFT regressed); **falsified for pH 3**, where FFT was FIXED (24.8× → 1.47×). |
| H-4 | `zhou_dimer_share_pH_invariance` improves | **HELD.** 187× → **72.4×**. Still fails; the direction is right and the mechanism is the one B2.1 named. |
| H-5 | 2 of the 3 `zhou_120C_dimer_share_pH*` rows pass | **FALSIFIED.** **0 of 3.** Two regressed. |
| H-6 | `zhou_pH6_MFT` improves by ≥ 3× | **HELD, by 32×** (0.0576× → 1.82×) and it now PASSES. This is the row the declaration's §5 decision 1 scored diagnostic-only "until the model carries a pH-trajectory state". |
| H-7 | the two low-temperature rows B2.1 fixed are at risk | **NOT REALISED.** `hofmann2002_brew_80C_FFT` (2.96×) and `vanseeventer_50C_MFT_per_day` (1.66×) both still pass. |
| H-8 | `cerny2007_branch_responsiveness` still fails | **HELD** (0.323×) — the ceiling is arithmetic, as B2.1 proved. |

### Blind re-sit (c) — the cutover final exam

| # | claim | outcome |
|---|---|---|
| X-1 | Yiltirak median falls below 100× | **FALSIFIED.** 290.5× → 262.3×, and the *worst* point got worse (515× → 1475×). |
| X-2 | Yiltirak still 0 / 8 | **HELD.** |
| X-3 | the Chang two-arm defect is not fixed | **HELD.** Still 4041 ppb for both arms; the acrylamide lane has no pH term and B2.2's pH state is sulfur-lane only. |
| X-4 | the Hofmann 145 °C family loses at least one row | **FALSIFIED.** It held 4 / 10 and its median *improved* from 7.09× to 4.16×. |
| X-5 | the core does not beat the old lane on the paired median | **FALSIFIED.** Core **10.65×** vs old **12.65×**. First time the core wins that comparison. |
| X-6 | the 17 declensions stay 17 | **HELD.** |

⚠️ **A scorer wording defect, reported not fixed** (the exam re-run was mandated
as pure re-scoring): `generate_cutover_final_exam.py` marks the pre-registered
claim "median 10×–100×, **and NOT better than the old lane**" as **HELD**,
because its check only tests the numeric band. Its own detail column then says
"the core is BETTER on the paired subset". The compound claim is half
falsified and the label is wrong. Several narrative paragraphs in that report
("median 290x"; "the sulfur CONSUMPTION channels carry no activation energy at
all") are also hard-coded from the B5 run and are now stale.

---

## 3. THE KANG 140 °C PASS SHOULD NOT BE COUNTED, AND HERE IS WHY

B2.1's diagnosis §2a localised its 140 °C collapse exactly: the model's *peak*
FFT was 1.28× the measured endpoint and its *endpoint* was 1.5 × 10⁻⁷ of it,
because sinks sharing the formation barrier kept running in a
precursor-exhausted pot. **That mechanism is gone.** Peak-versus-endpoint in
B2.2's own 140 °C run:

| | end / peak, B2.1 | **end / peak, B2.2** |
|---|---:|---:|
| MFT | ~0.002 | **0.893** |
| FFT | ~1.2 × 10⁻⁷ | **0.534** |

So Task 1 did the mechanical thing it was asked to do. **But it did it for a
reason opposite to the one the diagnosis predicted.** The thiol-sink barrier did
not go *below* the formation barrier; it went to 248 against 76.2. What the fit
actually bought is a sink whose *absolute* rate is tiny across the whole 100–145
°C window — it has to be, to survive Kumazawa's 121 °C survival grid and Kang's
100/120 °C rungs with a 248 kJ/mol slope — rather than a sink that decelerates
relative to formation. The product survives at 140 °C because the sink barely
runs at all below its 145 °C reference, not because the ratio of barriers is
right.

**And the pot's pH trajectory is wrong at exactly the rungs that were fitted.**
Model-only, from the frozen parameters:

| | in-situ pH start → end | cooled endpoint | measured (Kang p. 3242) |
|---|---|---:|---:|
| 100 °C | 5.59 → 9.69 | **11.42** | 4.9 |
| 120 °C | 5.31 → 9.40 | **11.42** | 4.9 |
| 140 °C | 5.05 → 5.43 | **5.42** | 4.9 |

(The 4.9 comparator is Kang's own measured pH drift, `kang2026_SI_extraction.md`
§5 quoting main text p. 3242. It is declared in **no** column of
`FIT_HOLDOUT_DECLARATION.md` — Amendment 7 declares only Zhou's three endpoints
— so it is used here as an **out-of-sample diagnostic**, never as a scored row.
It carries the same caveat the dossier attaches to Table S5: the temperature of
that pH series is `100_or_120_C_UNRESOLVED`, so it brackets the ladder rather
than pinning one rung.)

The two FIT rungs run alkaline and the hold-out rung does not. **A gating
hold-out passed by a model whose process state is six pH units wrong at the two
conditions it was fitted on is not evidence, and it is recorded here as not
evidence** — the same treatment B2.1's own exam gave the Chang acrylamide
"pass". The honest statement is that the Kang block went from *catastrophically*
wrong to *plausibly* wrong, which is worth having, and that one of its four rows
landing inside the band is a coincidence of a broken trajectory.

### The cause: a charge-bookkeeping hole, named and sized

The charge balance conserves the **strong-ion difference**, which is correct —
no reaction makes or destroys sodium. But it does **not** conserve *solute*
charge: when a charged solute is consumed and its products are uncharged
fragments, the sodium that titrated it is left with nothing to balance and the
pot drifts alkaline for a bookkeeping reason.

This wave found and fixed one instance before freezing — `r_ttca_deg` now
carries TTCA's carboxyl into the acid pool, which took Kang's 120 °C pot from a
predicted pH 11.4 to 4.6 in the pre-freeze probe. **It did not fix the others,
and the frozen fit then routed TTCA down the route that was still broken.** The
surviving instances are the cysteine sinks: `r_cys_thermal`, `r_cys_h2s`,
`r_dpo_ddp` and `r_cys_actz` all delete cysteine's α-carboxylate without
depositing an acid equivalent.

**NAMED FIX FOR THE NEXT WAVE, not applied here:** every step that consumes an
amino acid or amino-acid conjugate must carry its carboxylate into `ACID`, and a
conservation test should assert that the *total titratable acid + solute
carboxylate* is non-decreasing along any trajectory. **It was deliberately not
applied after the scores were seen** — a stoichiometry change made after reading
a hold-out scorecard is a refit in disguise.

---

## 4. THE FROZEN BUNDLES DO NOT RECORD THEIR BUFFER, AND THAT NOW COSTS AN ORDER OF MAGNITUDE

A schema probe over the 21 frozen external-validation bundles (**condition KEYS
only — no value was read**) returns exactly `{ph, temp_C, time_min,
water_activity}`. There is **no buffer field**. The engine's declared default
for an undeclared buffer is *unbuffered plus an extrapolation warning*, and that
warning does travel into the exam report. But the exam's Yiltirak bundles are
literally named `..._ribose_cysteine_buffer_...`, and they are being integrated
as if they were water.

Sized on the fit panel's own system definition — no bundle opened, no measured
value read:

| system | unbuffered (engine default) | 0.5 M phosphate (Hofmann's stated buffer) | clamped (B2.1) |
|---|---|---|---|
| pentose/Cys, 145 °C, 20 min | pH → 2.76; MFT 360 µg/L, FFT 91 | pH → 3.59; MFT 214, FFT 106 | pH 5.00; MFT 18, FFT 52 |
| Yiltirak-shaped, 100 °C, 4 h | pH → 3.13; MFT 595, **FFT 3895** | pH → 4.30; MFT 39, **FFT 341** | pH 5.00; MFT 27, FFT 180 |

**An 11× swing in predicted FFT on the buffer field alone, on the exact system
shape whose exam row is the worst point in the whole report (1475×).** The
Yiltirak blow-up is therefore no longer purely a chemistry failure; a large part
of it is now a *data-schema gap*. **Recommended, and not done here because the
exam was mandated as pure re-scoring: add a declared `buffer` block to the
frozen bundle schema.**

---

## 5. Contradictions found — reported, not adjudicated

1. **The brief expects buffered systems to hold their pH; the arithmetic says
   they cannot.** 0.5 M phosphate at pH 5.0 sits between its pKa₁ (2.148) and
   pKa₂ (7.198) and has a van Slyke capacity of ≈15.7 mmol L⁻¹ per pH unit. The
   modelled acid load is 86.7 mmol L⁻¹ in the Hofmann pentose pot and 269 in the
   Cerny ternary. The buffered systems drift 1.4–1.7 units. This was
   pre-registered as F-7 before the fit ran. **Only the McIlvaine
   citrate/phosphate systems (Kumazawa) hold, and they hold exactly** — their
   capacity is 56–147 mmol L⁻¹ per pH unit and they contain nothing that makes
   an acid.
2. **The calibrated Amadori ammonium pKa (10.77) is not chemistry.** An
   N-glycosylated secondary amine should sit near 7–8. The fit has an incentive
   to push it high because that makes the Amadori electrically invisible and
   lets the strong-ion difference be set by cysteine alone. **It is a fitted
   nuisance parameter that happens to live in a pKa slot, and it must never be
   quoted as an acid-base constant.**
3. **Amendment 7 dates the Zhou drift endpoints "120 °C / 2 h".** The dossier
   (`zhou2023_extraction.md` §3, quoting the Figure 2 caption verbatim) says
   **120 °C for 60 min**, which is also the condition of the Table 1 volatiles
   the endpoints qualify. This wave used **60 min**, the dossier's value. The
   amendment's "2 h" appears to be a transcription slip; it is reported rather
   than silently corrected in the declaration.
4. **B2.1's diagnosis §4 says Zhou's pH convergence is something "the module's
   pH-trajectory state already carries".** It did not: B2.1 carried a
   *prescribed* linear interpolation between two measured endpoints, which
   cannot be run where nobody published the endpoint. That sentence overstated
   what B2.1 had. B2.2 is the first version for which it is true.
5. **Only 3 of 48 parameters are individually identified**, on 58 rows with 10
   degrees of freedom. Neither decay barrier and neither drift constant is among
   them. `Ea_decay_thiol_sink` at 248.0 clears its 250 ceiling by 2 kJ/mol — it
   is a **bound-limited estimate, not a barrier**, and
   `tests/unit/test_kinetic_core_b2_2.py` pins it as such so it cannot be
   lifted into another module as a measurement.

---

## 6. Where the nine regressions are, and what they have in common

`zhou_pH8_MFTD` · `cerny2003_NF_share_ceiling_95C` ·
`cerny2003_intact_skeleton_share_95C` · `route_mix_moves_with_temperature` ·
`hofmann_ribose_pH7_MFT` · `hofmann_ribose_pH7_FFT` · `meynier_FFT_pH_fold` ·
`zhou_120C_dimer_share_pH7` · `zhou_120C_dimer_share_pH8`

Seven of the nine are **pH-response or route-mix rows**, and both structural
changes push on exactly those:

* the two-scale correction moves every unbuffered pot's operating pH down by
  ~1.1 units and every buffered pot's endpoint down by 1.4–1.7, so every row
  scored against a *nominal* pH label is now scored against a different
  operating point;
* a 248 kJ/mol thiol sink and a 150 kJ/mol carbonyl sink against a 76 kJ/mol
  formation lane invert the route mix's temperature response — which is
  precisely what `route_mix_moves_with_temperature` (a diagnostic row) and the
  two Cerny 95 °C share rows measure.

The dimer rows are the tension B2.1 already recorded as unresolved: Kumazawa's
buffered system says the disulfide branch is strongly pH-dependent, Zhou's
unbuffered system says the dimer *share* is not. B2.2 improved the invariance
row (187× → 72.4×) at the cost of two of the three level rows. **Both cannot
currently be satisfied and this wave did not adjudicate it either.**

---

## 7. Process note — the fit was run twice, and why

The first two attempts at this fit were lost to the container's 8 GB memory
ceiling after ~2 h each: repeated stiff `solve_ivp` calls accumulate memory that
`gc.collect()` does not recover, and at ~3 × 10⁵ integrations the process
reached 4.6 GB and was killed before any report was written. The shipped run
therefore differs from the first attempt in three declared ways — 9 pH nodes
instead of 17 in the optimiser's quick mode, `ftol` 1e-6 instead of 1e-8, and
**two random starts instead of three**. The seed is fixed, so starts 0 and 1
reproduce exactly; start 0 converged to cost 16.59 and start 1 to **8.99**,
which is the frozen point. A third start was budgeted and not run. **The
multistart claim in the fit report is therefore two starts, not three, and the
report says so.**

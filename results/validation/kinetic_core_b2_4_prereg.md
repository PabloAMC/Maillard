# Build Wave B2.4 — PRE-REGISTRATION

**Written and saved to disk BEFORE any fit ran, BEFORE the hold-out panel was
scored under any weighting, and BEFORE the cutover exam was regenerated under
any weighting.** Every claim below is scored in
`kinetic_core_b2_4_fit_report.md` and `kinetic_core_b2_4_holdout_report.md`
whether it held or not, and the count of each is printed.

Mandate: `results/validation/d1_exam_panel_reconciliation.md` §8, "What a B2.4
should change", items 1–4, and its drafted pre-registration text. This file is
that draft, adapted to what the container can actually run, with every
deviation from the draft stated in §6.

---

## 0. What this wave changes

**No stoichiometry, no species, no reaction, no benchmark, no engine input, no
hold-out row, no pass band.** Four things only:

1. **A DECLARED WEIGHTING.** A new fit generator,
   `scripts/generators/generate_kinetic_core_b2_4_fit.py`, carries a single
   declared scalar `PH_ENDPOINT_WEIGHT` with a written basis, applied to the
   three `zhou_final_pH_*` rows, and the fit is run at **three** declared
   values of it. The B2.3 generator is imported unchanged and is not edited:
   its systems, its 58 rows, its observables and its bounds are reused verbatim
   so that nothing but the weighting and the free-parameter set differs.
2. **AN ENSEMBLE, NOT A POINT.** Six seeded starts per weighting, over a
   declared free-parameter subset (§2). The ensemble spread is reported; the
   best member is reported as *a member*, not as *the* optimum.
3. **TWO SCORER CONDITIONS.** The exam artifact gains the **geometric-mean
   fold** beside the paired median; the panel artifact gains **median |log₁₀
   fold|** beside the gating count; and the **four shared rows** are declared as
   shared in both artifacts.
4. **THE AMINE-FATE PROBE**, rebuilt against the current species set or the
   basis rewritten to claim only what is demonstrable (§5).

### The exchange rate, and its basis

The B2.3 objective has an **undeclared** pH-unit-vs-log-fold exchange rate. Its
three `zhou_final_pH_*` rows are scored in pH units at σ = 0.25 while its 55
level and ratio rows are scored in log₁₀ folds at σ_log = 0.20–0.60. D1 §3
measured the emergent rate: one pH unit of endpoint error costs ~9× what a 3×
level miss costs. Nobody chose that; it fell out of two σ picked independently
in different units.

B2.4 makes it a declaration. Define

> **E — the number of decades of level error that ONE pH unit of endpoint
> error is declared to be worth** —

and implement it as `sigma_ph = SIGMA_LOG_REFERENCE / E`, with
`SIGMA_LOG_REFERENCE = 0.35`, the modal σ_log of the corpus's level rows (every
Hofmann 1998 level anchor carries exactly 0.35). Nothing else in the objective
changes. Three declared values, each with a basis:

| tag | E (decades per pH unit) | σ_ph implied | basis |
|---|---:|---:|---|
| **W-SHIPPED** | **1.40** | 0.250 | The accidental rate B2.3 shipped, made explicit. `0.35 / 0.25 = 1.40`. This value reproduces B2.3's objective **exactly** and is the control: it is what "declaring the weighting" costs when the declared value is the one already in force. |
| **W-HALF** | **0.70** | 0.500 | Half the shipped rate, i.e. the pH endpoints down-weighted by **4× in cost** — precisely the intervention D1's W-3 names. Chosen as the smallest intervention large enough to be visible against a 10-dof objective, and chosen for that reason and not for any property of its outcome. |
| **W-MEASURED** | **0.28** | 1.250 | The rate the corpus measures **on the fit side**. Kumazawa 2003's six DECLARED FIT retention rungs move 2-furfurylthiol from 0.995 to 0.110 survival across pH 3.0–6.4: 0.956 decades over 3.4 pH units, i.e. **0.281 decades per pH unit**. This prices a pH miss at the level consequence the fit corpus's own pH-response measurement says a pH miss has. **See Correction C1 below — this value replaces a first draft derived from hold-out rows.** |

E is a **declared exchange rate, not a measurement uncertainty.** σ_ph = 1.25 pH
units at W-MEASURED does not claim Zhou's endpoints are uncertain by that much
(the dossier's digitisation error is ±0.06). It claims that missing an endpoint
by 1.25 pH units should cost the same as missing a level by 0.35 decades,
because that is the level consequence such a miss carries. Anyone who disagrees
with the exchange rate now has a number to disagree with, which is the whole
point.

> **CORRECTION C1 — a firewall violation in this pre-registration's own first
> draft, caught by this wave's own test, before any fit at the affected
> weighting ran.**
>
> W-MEASURED was first written as **E = 0.32**, justified by D1 §4's
> observation that Hofmann's pH ladder moves 2-furfurylthiol 1.28 decades
> across 4 pH units. **That number is computed from HOLD-OUT rows** — Hofmann's
> pH-3 and pH-7 aqueous points are hold-out in this repo and their difference
> is precisely what the exam and the panel score. A fit weighting derived from
> it would have put a hold-out value into the objective, which is the one thing
> the firewall exists to prevent.
>
> `tests/unit/test_kinetic_core_b2_4.py::test_no_holdout_literal_appears_on_the_b2_4_fit_side`
> found it, as the literal `1.28` on the fit side, when the test was first run.
> The basis is re-derived from the fit side only, from Kumazawa 2003's six
> declared FIT retention rungs — which is in any case the better anchor, being
> the corpus's only measurement of a thiol's pH response with formation held
> out of the picture. It gives **0.281**, so the two bases happen to agree
> closely; that is a coincidence and it is **not** the reason for the change.
>
> **Order of operations, disclosed:** the correction was made while the
> W-SHIPPED ensemble was still running and **before any W-MEASURED member
> existed on disk**. No result at the affected weighting was seen before or
> after. The original 0.32 is recorded here so the correction is auditable
> rather than invisible.

### The cross-weighting comparator

**Total cost is NOT comparable across the three weightings** — they are three
different objectives. The quantity that IS comparable, and which this wave
reports for every member of every ensemble, is

> **Σr²_level** — the sum of squared residuals over the **55 non-pH rows only**,
> at their **unchanged** σ_log.

Every claim below that compares weightings compares Σr²_level, never total cost.

---

## 1. THE ROUTE, AND WHY IT IS NOT THE DRAFT'S

D1's draft asked for six starts on the full 48-parameter vector at each of three
weightings. **That is not runnable here and the arithmetic is on the record
before the wave starts.** B2.3's own `multistart_trace` reports 4731 s and
2280 s for its two starts: ~1.2 h per start on the full vector. Eighteen such
starts is ~21 h of container time, and the container's 7 GiB ceiling already
killed two B2.3 attempts.

D1's own fallback is therefore the route taken, and it is taken **for the
declared reason above and not because of anything seen in a result**: freeze
the parameters the corpus does not identify, and fit a declared subset with six
starts.

---

## 2. THE FREE SET — declared by rule, before any fit

**20 of 48 free, 28 frozen at their B2.3 values.** The rule, in four clauses:

| clause | what | count | why |
|---|---|---:|---|
| **R1** | Every parameter the B2.3 report's own Gauss-Newton intervals call **individually identified** | 11 | These are the only constants the corpus determines. `k_tdp_fur`, `k_nf_mft`, `k_nf_mp3p`, `k_mgo_mp`, `k_glc_ha`, `k_dimer_fft`, `k_fft_decay`, `k_fur_decay`, `k_pent_caramel`, `k_cys_thermal`, `k_thiolate_loss`. |
| **R2** | The two pH-drift constants `ph_acid_yield_per_sink_event`, `ph_arp_secondary_ammonium_pKa` | 2 | They are the **direct price-takers of the declared weighting**. Freezing them would make the weighting experiment vacuous — there would be nothing for a re-priced pH row to buy. |
| **R3** | The five rate constants D1 §3 records as having moved **≥3 decades** between B2.2 and B2.3 and not already in R1: `k_fur_fft`, `k_osone_decay`, `k_dimer_decay`, `k_ttca_deg`, `k_thiol_decay` | 5 | These are the **topology switches** D1 identifies as the regression's mechanism. Freezing them at B2.3 would hard-code the very basin under test, and would make D1's W-3 (`k_fur_fft` falls back below −2.0) unfalsifiable by construction. |
| **R4** | The two decay barriers `Ea_decay_thiol_sink`, `Ea_decay_carbonyl_sink` | 2 | D1 §3 shows one traded **onto** its 250 kJ/mol ceiling and the other **off** it in the refit. A bound trade is exactly what a basin change looks like. |

The 28 frozen constants are held at B2.3's shipped values, **unchanged, and not
re-derived**. This wave therefore measures basin choice *within* the subspace
the corpus can speak about, and says nothing about the 28 it cannot. That is a
narrower claim than the draft's and it is stated as such.

### The six starts

- **Start 0 is B2.3's own vector**, exactly, on all 20 free coordinates. Every
  ensemble therefore contains the incumbent, so a member that is worse than
  B2.3 is visibly worse and cannot hide.
- **Starts 1–5 are seeded random draws** (`numpy.default_rng(20260829 + 1000·i)`),
  **stratified into a LOCAL and a GLOBAL arm** — see Amendment A2 below.
- Optimiser settings are **identical to B2.3's**: `trf`, `x_scale="jac"`,
  `diff_step=3e-2`, `ftol=1e-6`, `xtol=1e-8`, `max_nfev=250`, quick pH mode
  (9 nodes / 2 iterations). Changing the optimiser between a fit and its
  re-weighting would confound the comparison this wave exists to make.
- A member that terminates on `max_nfev` (scipy `status = 0`) or on the
  evaluation budget of Amendment A1 is **reported and excluded from the spread
  statistic**, and the count of such members is printed. A budget exhaustion is
  not a basin.

### AMENDMENTS A1 and A2 — made BEFORE any member finished, and why

**Both amendments were written, and this file re-saved, while zero ensemble
members had completed and zero member results existed on disk.** They are
forced by a measurement of the container, not by anything seen in a result. The
measurement, which is the whole justification and is reproducible:

> One residual evaluation of this objective costs **1.68 s** with the box to
> itself. With five processes running it costs **8.64 s** — a 5.1× slowdown on
> five processes, i.e. the container delivers **about one core of real
> throughput** for this workload whatever `nproc` reports. Parallelism buys
> nothing here. B2.3's own full-vector fit implies ~3150 evaluations for one
> start; eighteen members at that size is ~11.7 h of serial compute.

> **A1 — A HARD EVALUATION BUDGET.** Each member is given **500 residual
> evaluations**, enforced inside the objective. `max_nfev` does not bound this:
> scipy counts trial points in `nfev` but the finite-difference Jacobian's
> evaluations are not counted there, and with 20 free parameters the Jacobian is
> where essentially all the cost is. On exhaustion the member returns its
> **best-so-far** vector with `status = -9`, is labelled `budget_exhausted`, and
> is excluded from the spread statistic by the rule already declared above for
> `max_nfev` terminations. **This wave expects most GLOBAL-arm members to
> exhaust the budget** and says so here rather than discovering it later.

> **A2 — STRATIFIED STARTS.** Starts 1–2 are the **LOCAL arm**: uniform within
> ±2.0 decades of B2.3 on each free rate constant, ±40 kJ/mol on each decay
> barrier, ±0.15 on `acid_yield`, ±1.0 on the pKa, clipped to the declared
> bounds. Starts 3–5 are the **GLOBAL arm**: uniform over the full declared
> bounds, exactly as originally written. Both spreads are reported separately
> and the pooled spread is what feeds the shipping criterion.
>
> The reason is that the two arms answer different questions and the original
> design could only afford to answer neither. A full-bound draw starts at cost
> 1200–2500 against an optimum near 8; under A1's budget such a member reports
> where it *got to*, not where the basin *is*, which makes a pooled spread a
> statistic about the budget. The local arm asks "is there more than one basin
> within two decades of the incumbent" and can actually converge inside the
> budget; the global arm keeps the original, harder question and reports
> honestly that it is budget-limited. **Damage, stated plainly:** the local arm
> measures curvature and near-basin structure, not the global landscape. W-1 is
> therefore scored on the **pooled converged** members and, separately, on the
> local arm alone, and both are reported.

---

## 3. THE SHIPPING CRITERION — declared before any fit ran, and EXAM-BLIND

The shipping value is chosen by the following, computed from the **fit side
only**. Neither the exam nor the hold-out panel enters it, in any form, at any
stage. Both are computed and reported under **all three** weightings.

> **Gate G1 — the calibration must survive.** At the ensemble's best-cost
> member, all three `zhou_final_pH_*` predictions are within **1.0 pH unit** of
> their anchors. A weighting that fails G1 has stopped being a calibration of
> the drift constants and is disqualified whatever its spread.
>
> **Gate G2 — the ensemble must be an ensemble.** At least **4 of 6** members
> converged (scipy `status ∈ {1, 2, 3}`).
>
> **Score S — the basin spread.** `S = log10(max cost / min cost)` over the
> **converged** members of that weighting's ensemble, each evaluated under its
> own objective.
>
> **Ship** the qualifying weighting with the **smallest S**. Ties within 0.05
> in S break toward the **largest E** — keep the pH endpoints weighted as
> heavily as the data will bear. If no weighting passes both gates, ship the
> smallest-S weighting among those passing G1 and say so in the report.

Rationale, stated in advance: D1's finding is not "B2.3 chose a bad weighting",
it is "**the objective was selecting a basin and nobody declared the rate that
made it do so**". The defect is irreproducibility. The criterion therefore
selects on reproducibility, gated on the calibration still being a calibration —
and explicitly not on the out-of-sample score, which is what a hold-out is for.

---

## 4. THE CLAIMS

D1's W-1 … W-7 are carried. W-1 and W-2 are **restated** for the route in §2
(the draft's thresholds were set for a 48-parameter, fully-random ensemble and
would not mean the same thing over a 20-parameter subset containing the
incumbent); the restatement is here, before the fit, and the draft's original
wording is quoted beside it so the change is auditable.

| # | claim | falsified if |
|---|---|---|
| **W-1** | *Draft: "at the shipped weighting, six starts produce a spread in final cost of more than 2.0."* **Restated:** at **W-SHIPPED**, the six starts produce a spread `S = log10(max/min cost)` over converged members of **more than 0.30** (i.e. the worst converged member costs >2× the best). B2.3's two starts gave 10.38 / 8.18 (S = 0.104) on the full vector; B2.2's gave 16.59 / 8.99 (S = 0.266). The two-start protocol was selecting a basin, and six starts over the identified subset should show at least as much scatter as two starts showed over the full one. | the six converged members span less than 0.30 in S |
| **W-2** | Across the six-member ensemble at **W-SHIPPED**, the exam's **geometric-mean fold** spans a range of at least **1.5×** (max/min over members). If it does, "50.13×" was never a property of the model. | the range is below 1.5× |
| **W-3** | At **W-HALF** and at **W-MEASURED**, the ensemble-best `k_fur_fft` is **below −2.0** (B2.3: −0.576; B2.2: −4.257) **and** the Hofmann pH3→pH7 FFT slope is **below +3.0 decades** (B2.3: +4.84; measured: +1.28). | either fails at either weighting |
| **W-4** | At every weighting tried, the **charge-conservation fix stays free**: this wave predicts, in advance, that its predecessor's conservation fix is not the problem. Scored as the D1 §1 `S_only` cell — the exam's geometric-mean fold at B2.2's parameters under B2.3's stoichiometry — staying within 10% of **26.8×**. | it moves by more than 10% |
| **W-5** | **THE WAVE'S PREDICTION THAT ITS OWN REMEDY IS INSUFFICIENT.** The Hofmann pH-3/pH-7 hold-out block does **not** reach 3/4 at any weighting, and the model's formation-vs-pH ladder stays **humped** (a maximum at label pH 4.0–5.0) against the measured **monotone-falling** ladder, at every weighting. Re-weighting is not a fix for a shape defect. | the block reaches 3/4 at any weighting, or any weighting's ladder is monotone falling from pH 3 to pH 7 on both FFT and MFT |
| **W-6** | `fed_c2c3_MFT_pH3` remains the **largest single residual** in the fit at every weighting, because the route it scores still carries no pH factor. | it is not the largest at some weighting |
| **W-7** | The exam and the hold-out panel, once both report **median \|log₁₀ fold\|**, agree in **SIGN** on every weighting tried (both improve relative to W-SHIPPED, or both worsen). If they ever disagree in sign, the four shared rows are masking something and the wave says so. | they disagree in sign at any weighting |
| **W-8** | *(new, this wave)* **Σr²_level improves monotonically as E falls.** W-SHIPPED ≥ W-HALF ≥ W-MEASURED on the 55 non-pH rows at unchanged σ_log, at the ensemble-best member. This is the mechanical prediction of D1-A: the pH endpoints were buying level accuracy, so stop paying for them and the levels come back. | the ordering is not monotone |
| **W-9** | *(new, this wave)* **The three zhou_final_pH residuals grow monotonically as E falls**, and **W-MEASURED fails gate G1** (some endpoint misses by more than 1.0 pH unit). The measured exchange rate is predicted here, in advance, to be too loose to keep the drift calibration — which would make W-HALF the shipped value. | W-MEASURED passes G1, or the residuals are not monotone in E |
| **W-10** | *(new, this wave)* **The shipping criterion and the exam disagree.** The weighting the pre-declared criterion selects is **not** the weighting with the best exam geometric-mean fold. Recorded in advance because a criterion that always agrees with the score it is meant to be blind to is not blind. | the criterion selects the best-exam weighting |

### What would make B2.4 a failure

- **W-1 falsified** — the objective was not multi-modal over the identified
  subset and D1's basin diagnosis is wrong at least there; or
- **the ensemble's exam spread being narrow while the between-wave difference
  stays large**, which would mean the B2.2→B2.3 regression is a real physics
  change and D1 mis-attributed it; or
- **W-8 falsified** — the levels do not come back when the pH endpoints are
  cheapened, which would mean D1-A named the wrong mechanism and the objective's
  currency mismatch is not what re-routed the network.

---

## 5. THE AMINE-FATE PROBE

D1 §7 found that `scratch/b23_encoding_probe.py`, cited by
`ph_state.AMINE_FATE_BASIS` and by the B2.3 pre-registration §0 as the evidence
for the wave's "strongest single assumption", **cannot run on the shipped
tree**: it raises `KeyError: 'AMN'`, and two of its three encodings collapse
onto the same code path because there is no amine pool to patch.

**Declared before doing it:** the probe is rebuilt against the **current**
species set. The liberated amino nitrogen is not a species, but it is
**derivable**: the amine centres in play are exactly those on `Cys`, `ARP` and
`TTCA`, so the amount released at time *t* is
`(initial amine total) − (Cys + ARP + TTCA at t)`, computable from the
concentration vector alone with no new species and no new parameter. The rebuilt
probe injects that quantity as ammonium (pK_a 9.25) to realise the "carry both
centres" encoding, and compares it against the shipped encoding at B2.2's own
frozen parameters, on Amendment 7's three declared FIT anchors.

**Pre-registered expectation:** the rebuilt "carry both" encoding leaves Zhou's
three cooled endpoints **above pH 4.5** on all three anchors (against measured
3.22 / 3.42 / 5.07), reproducing the *direction* of the B2.3 prereg's third row
(5.04 / 5.11 / 5.89) but **not necessarily its digits** — the original ran
against a species that no longer exists and the reconstruction is not the same
object. **If the rebuild fails to reproduce the direction, `AMINE_FATE_BASIS` is
rewritten to claim only the arithmetic argument (40 mmol/L of surviving
ammonium at pK_a 9.25 cannot reach pH 3.2) and the probe citation is struck.**
Either way the citation is corrected to point at a probe that runs.

---

## 6. DEVIATIONS FROM D1'S DRAFT, DECLARED IN ADVANCE

| draft | this wave | why |
|---|---|---|
| six starts on all 48 parameters | six starts on a **declared 20-parameter subset**, 28 frozen at B2.3 | §1: 18 full-vector starts is ~21 h of container time against a 7 GiB ceiling that has already killed two attempts. This is the draft's own stated fallback. |
| "2–3 declared weightings" | **3**: 1.40, 0.70, 0.28 | the control is worth a slot: without it the comparison has no anchor. |
| W-1 threshold "spread > 2.0 in final cost" | **S = log10(max/min) > 0.30** over converged members | a 20-parameter ensemble containing the incumbent cannot be compared to a 48-parameter one on an absolute cost gap; a ratio is the scale-free statement of the same claim. |
| draft is silent on how the shipping value is chosen | §3, a two-gate + spread criterion, exam-blind | the mandate requires the criterion to be written before the results are seen and not to be the exam score. |
| draft is silent on cross-weighting comparison | **Σr²_level** on the 55 non-pH rows at unchanged σ_log | total cost across three different objectives is not a comparison. |

---

## 7. FIREWALL

Unchanged. `data/benchmarks/external_validation/` is opened by the exam
generator and by nothing else in this wave; the fit generator, the ensemble
driver, the amine probe and the analysis scripts are literal-grepped against
that path and against every hold-out row id before they run. No hold-out value
is read outside a scorer. The B2.3 fit generator, the B2.3 hold-out generator
and the cutover exam generator are the **existing** scorers and are reused; any
edit to them in this wave is additive (a new optional output basename, a new
reported statistic) and is listed in the fit report.

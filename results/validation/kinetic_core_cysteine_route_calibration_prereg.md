# Cysteine-pentose route calibration (stage 1) — pre-registration

*Written 2026-10-29, before any number of this calibration existed. Owner, 2026-10-29: "Let's do this modeling
correctly"; "get clarity before expensive runs". Context: `tasks/data_restructure_plan.md` §7, "THE SULFUR LANE,
REVISED UNDERSTANDING"; `kinetic_core_b59b_prereg.md` Amendments 2-4. Generator:
`scripts/generators/generate_kinetic_core_cysteine_route_calibration.py`. Nothing in the shipped engine changes:
the route stays off and the outputs feed stage 2 (the joint flavour-level refit, pre-registered separately).*

## 1. Question

Can de Roos's route, with the ring opening and cysteine's own loss, describe the measured time courses of ITS OWN
SPECIES — TTCA, the cysteine Amadori compound, free cysteine, free pentose — in three laboratories, at 70-140 °C
and starting pH 5.5-7.4? If so, with which constants? Flavour levels are NOT data here. They sit three or more
steps downstream, behind constants the old structure tuned; using them would let compensation back in.

## 2. What is calibrated, what is held

**Free (8), starting values and bounds.** Bounds: ±3 decades around the start for rates, the same width as the
codebase's identification threshold for a log10 k; barriers inside the sugar-trunk class band, 40-135 kJ/mol
(`FORMATION_EA_BOUNDS_BY_ROUTE["sugar_trunk"]`, inherited); K inside its measured span.

| coordinate | start | bounds | currently from |
|---|---|---|---|
| log10 k_ttca_cys (TTCA → Cys + pentose, at 145 °C) | B60 (0.264 1/min) | ±3 | fitted on Kang free Cys under the old structure |
| Ea of the ring opening, shared by `k_ttca_deg` | 64.08 (the lumped formation Ea) | 40-135 | lumped, never fitted on its own |
| log10 k_ttca_deg (TTCA → DPO without free Cys) | B60 | ±3 | fitted (B2.1) for Kang's ceiling |
| log10 K (TTCA equilibrium, L/mmol) | 0.060 | [0.019, 0.075] | Cao 0.060 (band 0.045-0.075); Zhai 2019 plateau ≥ 0.019 |
| log10 k_form (TTCA → Cys-Amadori, at 90 °C) | 3.29e-3 | ±3 | V2 anchor on Cao's pot |
| Ea of formation | 100.9 | 40-135 | Zhang 2026, glucose system (transfer) |
| Ea of the Cys-Amadori decay | 85.7 | 40-135 | Zhang 2026, glucose system (transfer) |
| log10 k_cys_thermal (Cys → fragments, at 145 °C) | B60 | ±3 | fitted; its Ea 55.1 is Kang's, measured, held |

**Held.** The decay's pH profile is de Roos 2005 Fig. 4 at ~100 °C, pH 3-6, implemented as printed (Amendment 4
of B59b); above pH 6 there are two arms (sec. 4). The DPO share is 0.30 (Zhai 2019; CR-12). Every other constant
of the sulfur lane stays at B60's isothermal optimum, including the pH drift (`ph_acid_yield_per_sink_event`).
This calibration is therefore CONDITIONAL on B60's downstream chemistry, declared. Stage 2 re-estimates everything
jointly, with these rows in its objective. Numerics: the inline pH coupling, rtol 1e-7.

## 3. Data (the fit rows), each pot as its source describes it

| source (lab) | pot | species, times | rows | read error (floor) |
|---|---|---|---|---|
| Cao 2017 (Xie, BTBU) | Cys + Xyl 0.40 M, 0.2 M phosphate pH 5.5, 90 °C | TTCA, Cys-Amadori (Fig. 5), free Cys, free Xyl (Fig. 4) at 0.5-5 h | 24 | ±5 mM (Fig. 5), ±4 mM (Fig. 4) |
| Zhai 2019 Fig. 3A (Jiangnan) | Xyl 165 + Cys 82.6 mM, water, pH 7.4 (NaOH), unbuffered; started from the measured t = 0 state (TTCA 14.5 mM) | TTCA, Cys-Amadori at 20-160 min, 70 / 90 / 110 / 130 °C | 64 | ±0.7 / ±0.3 mM |
| Zhai 2021 (Jiangnan) | TTCA 10 mM, water, pH 7 | TTCA left at 60 min from the zero-order fits, 100 / 120 / 140 °C (B16's rows) | 3 | — |
| Kang 2026 SI Fig. S3 (= Zhai 2022, Jiangnan) | TTCA 10 mM, water, pH 7 | free Cys at 20-120 min, 100 / 120 / 140 °C | 18 | ±0.03 mM |
| Kang 2026 SI Fig. S4 | free Cys 10 mM, water, pH 7 (read as a free-cysteine control, the dossier's 85 % reading) | free Cys at 20-120 min, 100 / 120 °C (the 140 °C rung stays a hold-out, as in B60) | 12 | ±0.03 mM |

All pots carry the fit systems' ambient dissolved oxygen. **Measurement model** (the owner's C3 design: one error
scale per quantification class). Every row here is a calibrated HPLC / LC-MS concentration. Residuals are
log10((pred + e) / (meas + e)), where e is the figure's read error (the resolution of the reading, so that a
measured 0.0 is not a log of zero); for the rows taken from fitted lines, e = 0. One common sigma is estimated from
the residuals at the optimum. It does not move the optimum and is used only for the Laplace step. Points within a
digitised series share calibration and are not independent. That is handled by the leave-one-source-out check
(sec. 4), not by a weighting choice.

## 4. Arms, checks and rules

**Arms.** CR-14: the decay above pH 6 (a) extended at de Roos's measured 5.5 → 6 slope, (b) held flat at its pH 6
value. Both are fitted. One is preferred only if its cost is lower by more than 0.5 × 3.84 × the reduced
chi-square: the chi-square(1) 95 % threshold on these direct-species rows, the same threshold the envelope uses.
Otherwise both are carried forward.

**Robustness.** Leave-one-source-out refits (Cao, Zhai 2019, Zhai 2021, Kang S3, Kang S4) on the preferred arm,
or on both if neither is preferred. A constant that moves by more than 2 sigma when one source is dropped is
recorded as **source tension**, with the source named. It is not tuned away.

**Identification.** The Laplace covariance at each optimum (sigma² (JᵀJ)⁻¹ from the final Jacobian), with the
codebase's rule (sigma within 3 decades for a log10 k, 60 kJ/mol for a barrier; not pinned at a bound). Bounds hit
are reported as such.

**Verification (before the fit).** V1: route off, every existing test and the B60 cost (27.08) unchanged. V2:
the decay table reproduces de Roos's Fig. 4 constants at its eight pH values. V3: at the start vector, the
objective's rows are computed and listed (the "no refit" baseline).

**Pre-declared checks, not fitted (recorded at each optimum).**
- Zhang 2025b blank (Xie lab; Cys + Xyl 40 mM, 0.2 M phosphate pH 6.5, 120 °C): TTCA 11.2 → 18.2 and Cys-Amadori
  0.07 → 1.61 mM at 10-50 min.
- Zhai 2019 Fig. 3B (TTCA 5 mM, 120 °C, pH 7.4): 1-DX and 3-DX at 10-60 min, and their ratio.
- pH: Cao's blank trajectory (5.49 → 4.67 at 0.5 h → 4.15 at 5 h, room-temperature readings) and Kang's final pH
  (TTCA 10 mM, 120 °C, 2 h: initial 5.5 / 7 / 8 → 5.1 / 4.9 / 4.7). The pH drift is held at B60, so a miss here
  is a stage-2 item, not a stage-1 failure.
- Zhao 2017: the cysteine pool (Cys + TTCA + Cys-Amadori) conserved within ~2 × 10⁻⁴ 1/min at 90 °C.
- Flavour-level diagnosis at the stage-1 constants (B60 downstream): the 52-row cost; Schieberle's 100 °C shape
  and his 100 °C / 145 °C factor under each of the five B64 heating readings; Nestlé's 110 °C shape; the Hofmann pH
  3 / 5 / 7 sweep.

**Expectations, written now.** The ring opening comes out roughly 10× slower at 100 °C than B60's (B16's and
Amendment 3's diagnosis). The formation barrier falls below 100.9, towards 55-80. K sits below Cao's 0.060. The
decay's barrier is weakly identified. Kang's free-cysteine rise-and-fall is reproduced within ×2, with
`k_ttca_deg` carrying most of TTCA's loss. The Kang final pH stays too high.

**What goes forward (the rule).** Stage 1's constants go to stage 2 as its start, and its rows join stage 2's
objective, when: every fit has converged or reached 1 200 evaluations (reported as which); V1-V3 pass; and an
independent reviewer finds no error in sections 1-4. Checks and costs do not gate. If a source is in tension, both
the full and the leave-out constants go forward as alternatives.

## 5. Results

*(appended after the run; nothing above this line is edited after that)*

### Run log 1 (2026-10-29): first fits withdrawn, false convergence

The first fit of each arm reported convergence (extend 20.34 → 6.57, flat 19.92 → 8.27). A restart of each arm from
the other's optimum then went lower in both (extend 6.39; flat 4.45), so neither had reached an optimum, and the
arm preference reversed. Cause: the generator left `least_squares` at its default finite-difference step (~1e-8
relative), far below the integrator's noise at rtol 1e-7, so the Jacobian was mostly noise. Every fit in the
codebase uses `diff_step = 3e-2` (B2-B8); this one now does too. The four results are kept, unused, in
`kinetic_core_cysteine_route_calibration/withdrawn_no_diff_step/`. **Convergence criterion, added now, before
the reruns:** each arm is fitted from two starts, the pre-registered start and the lowest-cost point found so far
for that arm. The arm counts as converged when the two costs agree within 1 % and every constant within its
Laplace sigma. If they do not, both are reported and the lower is carried, flagged.

### Run log 2 (2026-10-29): converged; arms identical; K on its bound

With `diff_step = 3e-2`, all four fits (two arms × two starts) reach the same optimum: cost 20.34 → 2.33, every
constant within a fraction of its sigma. Both arms are converged by the criterion of run log 1. Arm costs 2.327
(extend) and 2.326 (flat): neither is preferred, and both are carried. They coincide numerically, because these
pots spend little time above pH 6. RMS by source: Cao 0.26, Kang S3 0.22, Kang S4 0.06, Zhai 2019 0.18,
Zhai 2021 0.22 dex; common sigma 0.20 dex. Constants (extend arm; flat within 0.01):

| constant | value ± sigma | before |
|---|---|---|
| k_ttca_cys at 145 °C | 10^-1.96 ± 0.16 (0.011 1/min; 2.5e-3 at 100 °C) | 0.264 (B60; 0.029 at 100 °C) |
| Ea ring (opening and k_ttca_deg) | 42.8 ± 10.0 | 64.1 (lumped) |
| k_ttca_deg at 145 °C | 10^-1.91 ± 0.13 | 10^-2.21 (B60) |
| K | 0.075, **on the upper bound** | 0.060 (Cao) |
| k_form at 90 °C | 10^-2.92 ± 0.05 (1.2e-3) | 3.29e-3 (V2 anchor on Cao alone) |
| Ea formation | 42.5 ± 8.6 | 100.9 (glucose transfer) |
| Ea decay | 123.5 ± 15.4 | 85.7 (glucose transfer) |
| k_cys_thermal at 145 °C | 10^-2.32 ± 0.55 | 10^-2.12 (B60) |

Systematic residuals: Cao's Cys-Amadori −0.33 to −0.57 dex at 1-5 h, while Zhai 2019's at 90 °C is +0.14 to +0.27.
Every Xie-lab pot (Cao; Zhang 2025b) is in 0.2 M phosphate and every Jiangnan pot is in water, so phosphate
catalysis of the Amadori step (a condition; de Roos calls the step general-acid-catalysed) and a between-lab
calibration difference are confounded in the corpus. Recorded as an open hypothesis, not fitted. The first sample
of each Kang S3 and Zhai 2019 series is off by 0.2-0.5 dex, consistent with a few minutes of heat-up in oil-bath
vessels (not modelled here). Zhai 2019 130 °C Cys-Amadori at 140-160 min: +0.6 (measured 0, read-error floor).

**K on its bound (declared before running).** The upper bound, 0.075, is the top of Cao's measured band and was a
choice. A sensitivity fit with K's upper bound at 0.60 (8× Cao's value) is run once on the extend arm. If K moves
far above Cao's band, the Jiangnan water pots want a different equilibrium than Cao measured in phosphate; that is
recorded as part of the confound above, and the shipped value stays Cao's measured band (a measurement over an
inference).

### Run log 3 (2026-10-29): leave-one-source-out, the K sensitivity, the checks

**Leave-one-source-out** (both arms; they agree to 0.1 sigma). Moves beyond 2 sigma: without Zhai 2019,
k_form(90 °C) +5.1 sigma and the decay barrier −5.4 sigma (**source tension, Cao vs Zhai 2019**: the formation
rate is ~5× higher in the Xie lab's phosphate pots than in Jiangnan's water, and the decay barrier rests on Zhai's
110-130 °C series alone); without Cao, k_form −1.9 sigma (the same tension from the other side). Without Kang S4
or without Zhai 2021, k_cys_thermal falls to its floor (−4 to −5 sigma; sigma 0.55 dex): weakly identified, held
by those two sources. Without Kang S3, the ring opening +2.1 sigma. Nothing else moves beyond 2 sigma.

**K sensitivity** (upper bound 0.60): K = 0.150, cost 2.33 → 2.17, a fall larger than the 0.5 × 3.84 × chi²_red
threshold (0.08); k_cys_thermal again to its floor. As declared: the Jiangnan water pots want a higher K than Cao
measured in phosphate (part of the same confound). The carried value stays Cao's band (0.075, its top).

**Checks at the full optimum (extend arm; flat identical within 0.01 dex):**
- Zhang 2025b blank (Xie lab, within-lab): TTCA −0.23 to −0.37 dex; Cys-Amadori −0.9 to −1.1 at 20-50 min.
  The same direction as Cao: both Xie pots, both in phosphate, form the Amadori compound faster than the fit.
- Zhai 2019 Fig. 3B (TTCA 5 mM, 120 °C): 3-DX 0 to +0.27 dex; **1-DX −0.46 falling to −2.09 dex**; predicted
  3-DX/1-DX 594 vs measured 2.5. The 1-deoxyosone (the MFT precursor) is destroyed far too fast downstream.
- Zhao 2017: cysteine pool lost at 2.1e-3 1/min at 90 °C, against ≤ 2e-4 measured (10× too fast).
- pH: Cao's trajectory 0.5-0.9 units too acidic. Kang's final pH 4.0 / 7.5 / 8.9 from 5.5 / 7 / 8 against
  5.1 / 4.9 / 4.7: the drift cannot acidify an unbuffered pot. B60's acid yield is 7.5e-5 per sink event, i.e. the
  model makes almost no organic acid, and in alkaline pots the opened ring releases cysteine's strongly basic amine.
- Flavour level (downstream at B60): the 52 rows cost 69.4 (B60 27.1; route on at B60's constants 121.5). Worst:
  Zhai's 120/100 °C MFT and FFT folds (+5.3, +5.1 sigma), Kang's 120 °C FFT. Schieberle's MFT at 60 / 360 / 720
  min over 30: 12.7 / 345 / 171 (measured 3.1 / 35 / 40). The 100 °C / 145 °C factor by heating reading:
  isothermal 0.96, ramp + hold 0.75, Newtonian τ 4 4.9, τ 6.7 28.5, linear ramp 235 (measured 7.9-13). Nestlé
  0.070 / 0.46 / 1.26 (measured 0.008 / 0.067 / 0.31). Hofmann MFT at pH 3 / 5 / 7: −0.40 / −0.45 / −1.95 dex
  (B60 −1.35 / −0.32 / −3.69).

**Reading (for stage 2, not acted on here).** (1) One downstream constant likely explains the 1-DX and
cysteine-pool misses: `k_dpo_ddp`, the cysteine Strecker step on the 1-deoxyosone, fitted under the old structure,
consumes both. Zhai's dicarbonyl series and Zhao's pool are direct measurements of exactly this and belong in
stage 2's fit. (2) The pH of unbuffered pots needs a mechanism, not a refit: measured organic-acid yields of sugar
degradation, with Cao's, Kang's and Zhou's pH trajectories as data. (3) The Amadori formation's phosphate-or-lab
confound stays open; both the full constants and the without-Zhai-2019 constants (Cao-like) go forward as
alternatives, by the rule of sec. 4.

### Independent review (2026-10-29): PASS WITH CHANGES, and the corrections

Reviewer verified all 121 rows and the check data against the dossiers (no transcription, unit or threshold
error), the convergence by run log 1's criterion, and the leave-out sigma counts. No error in sections 1-4.
Corrections to run logs 2-3, as found:

1. **V1 and V2, now recorded.** V1: route off, the B60 isothermal cost is 27.0841 (unchanged). V2: the decay table
   reproduces de Roos's eight Fig. 4 constants (unit test, relative 1e-9). V3: `baseline_start_vector.json`.
2. **The cysteine-pool miss is stage 1's own `k_ttca_deg`, not `k_dpo_ddp`** (reading 1 withdrawn). Pool loss
   2.06e-3 /min; with `k_ttca_deg` / 1000 it falls to 1.5e-4 (inside Zhao's limit); with `k_dpo_ddp` / 1000 it stays
   at 1.3e-3. `k_ttca_deg` was fitted at pH 7 and 100-140 °C and carried down to Cao's 90 °C, pH 4.2-4.7 pot on the
   shared 42.8 kJ/mol barrier; Cao's fitted TTCA and free-cysteine residuals already drift with time. For the
   1-deoxyosone, `k_dpo_ddp` / 1000 raises the 60-min prediction only from 0.79 to 15 µM (measured 98): about 1.3
   of the 2.1 dex.
3. **The heat-up reading had the wrong sign.** Heat-up delays a series and predicts positive early residuals;
   Zhai's first samples are negative (TTCA −0.23 / −0.17 / −0.14), as are Cao's at 30 min and Zhang's at 10 min. Only
   Kang S3 fits heat-up. The pattern is fast equilibration in the formation pots against slow cysteine release in
   the TTCA-fed pots; the fitted ring opening at 90 °C (1.7e-3 /min) is 4× below the Cao dossier's floor of
   ≥ 7e-3 /min.
4. **"~5× faster" was a shift in sigmas, not a fold.** Without Zhai 2019, k_form moves 10^-2.919 → 10^-2.648,
   ×1.87. The 3-6× is the residual gap at the joint fit.
5. **The bound flags were wrong** (absolute tolerance 1e-6). The generator now uses B8's rule (on a bound within
   0.1 % of the band width; unidentified if on a bound with sigma above a tenth of the band), re-applied to every
   stored fit (`--reflag`). Full fit: only K is on its bound. Without Zhai 2019: the ring, formation and decay
   barriers and K are all on bounds, so the "−5.4 sigma" was the decay barrier LOSING identification, not source
   tension; Zhai 2019 alone identifies the three barriers. Without Cao: the ring barrier on its floor. In the full
   fit the ring and formation barriers sit 0.3 sigma above the 40 kJ/mol floor, where a symmetric Laplace sigma
   misdescribes them; and the sigmas ignore correlation within a digitised series, so the 2-sigma tests are
   anti-conservative.
6. **"The drift cannot acidify an unbuffered pot" is false.** Zhai 2019 (unbuffered) falls from 7.4 to 4.9
   cooled at 90 °C and 3.8 at 130 °C; Kang's pH 5.5 pot overshoots to 4.0. Kang's pH 7 and 8 pots fail to fall. The
   acid comes from TTCA formation and amine-centre speciation, so the pH item belongs to the centre ledger, not
   only to acid yields.
7. **3-DX at 10 min is −0.24, not "0 to +0.27".**
8. **Design gaps (sections 1-4 are not wrong, but incomplete):** (a) Zhai 2021 Fig. 4b measures TTCA degrading
   faster at pH 7-8 than at 5.5-6 at 100 °C, and `k_ttca_deg` / `k_ttca_cys` carry no pH law (audit item 5). (b)
   Kallen 1971 measures the thiazolidine K rising ~20× from pH 5.6 to 7.5 (formaldehyde analogue; the shape
   transfers, not the magnitude): pH is a third confound, with phosphate and lab, for K and for k_form, which is
   pH-independent in the engine. Concentration also differs (0.4 M vs 0.08 M). (c) The formation barrier (42.5
   kJ/mol) is far below the measured analogue (100.9), near its floor, and outside the pre-registered expectation
   of 55-80; the decay barrier (123.5) rests on Zhai 2019 alone, whose pH was never measured, while the model's pH
   in Cao's pot is 0.5-0.9 units too acidic.

Two further findings of mine, made while the review ran: Davidek 2002 measures phosphate catalysis of an Amadori
compound's degradation in one lab (0.1 M phosphate at 90 °C: ×3 / ×3.5 / ×2 / ×1.3 at pH 5 / 6 / 7 / 8, dihydrogen
phosphate as the catalyst, after Potman & van Wijk 1989). So de Roos's decay, measured in 0.5 M phosphate and applied
here to water pots, is a buffer transfer that stage 1 did not declare. And Cao's glycine arms do not show strong
general-acid catalysis of the Amadori formation (the Amadori compound at 1 h falls with glycine).

## 6. Verdict (stage 1)

**The procedure passes; the constants do not go forward as estimates.** By sec. 4's rule the stage-1 constants would
go to stage 2. But the review shows the apparent Cao-Zhai conflict and the K tension are plausibly the chemistry of
pH and buffer, which the corpus partly measures (Kallen's K(pH) shape, Zhai 2021's TTCA loss vs pH, Davidek's
phosphate catalysis, de Roos's decay vs pH), and which stage 1 omitted. Fitting two labs' data with those effects
missing would build compensations into stage 2. Stage 1's constants are kept only as the start vector of a
redesign, **stage 1b**, pre-registered separately: the route with pH laws for K and the TTCA ring kinetics,
phosphate catalysis of the Amadori steps, de Roos's decay corrected from 0.5 M phosphate to the pot's buffer, the
pH trajectory (Cao's measured series) as rows, and Zhai 2021 Figs. 4b and 5 digitised as data.

### Run log 4 (2026-10-29): the shared fitting module reproduces stage 1

The generator now fits through `src/kinetic_core/fitting.py` (absolute per-kind finite-difference steps, B8's
central-difference Laplace and identification rule, leave-one-group-out with a predictive check). Extend arm, from
the pre-registered start and from the stored optimum: cost 20.34 → 2.327 (95 evaluations) and 2.327 → 2.327 (46);
reduced chi-square 0.0412; every constant within 0.2 sigma of the stage-1 optimum and of each other (e.g. ring
opening 10^-1.999 / 10^-1.983 vs 10^-1.964; formation barrier 43.5 / 43.0 vs 42.5; decay barrier 125.3 / 124.2 vs
123.5; K on its bound in all three). Converged by run log 1's criterion. Files:
`kinetic_core_cysteine_route_calibration/extend_last_slope_fitting_module.json` and `..._second_start_fitting_module.json`.

### Run log 5 (2026-10-29): a condition error in three fit rows

Zhai 2021's three TTCA rows (zero-order fits at 100 / 120 / 140 °C) were entered as pH 7 pots, following the
dossier's heading. The paper ran them at TTCA's natural pH of 5.5, unadjusted (text layer; the newly digitised
Figs. 3a and 4b, `zhai2021_extraction.md` sec. 5). At pH 7, TTCA at 100 °C falls to 0.89 mmol/L by 120 min against
the pH 5.5 line's 7.08: the rows understated the pH 7 decay about tenfold. Every stage-1 constant touched by those
rows (the ring opening and `k_ttca_deg`, their shared barrier) is affected; this is a further reason the stage-1
constants are not forwarded (sec. 6). The dossier is corrected. Stage 1b uses the digitised Fig. 4b series instead
(TTCA alone and + xylose, pH 5.5 / 6 / 7 / 8, 100 °C, 7 times each) and Fig. 3a at pH 5.5 for 120 and 140 °C. B16,
a closed wave, used the same three rows as pH 7; recorded here, not re-run.

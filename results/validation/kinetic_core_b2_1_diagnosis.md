# Build Wave B2.1 — what the numbers said, AFTER the scorecard was frozen

**Everything in this file was written after looking at the hold-out scores.
NO PARAMETER WAS CHANGED, before or after.** The frozen fit is
`results/validation/kinetic_core_b2_1_fit_report.json` and the scorecard is
`results/validation/kinetic_core_b2_1_holdout_report.md`; both are unchanged by
anything here. This document exists so that a post-hoc explanation is filed as
a post-hoc explanation, in its own file, rather than being folded back into the
pre-registered report where it would read like a prediction.

The two diagnostics below were run with `scratch/b21_kang_diag.py` and
`scratch/b21_branch_diagnosis.py` against the frozen parameters.

---

## 1. The headline

| | B2 | B2.1 |
|---|---:|---:|
| gating, on **B2's own 27 rows** | **7 / 27** | **15 / 27** |
| gating, on the full B2.1 exam | — | **15 / 33** |
| diagnostic | 1 / 7 | 6 / 9 |
| rows FIXED | — | **11** |
| rows REGRESSED | — | **0** |
| fit reduced chi-square | 10.80 | **2.01** |

Eleven rows moved from FAIL to PASS and **none moved the other way**. The six
gating rows B2.1 added — Kang's 140 °C block and Sun 2019's pH-9 block — all
FAILED, and one of them failed in a way that is worth more than the eleven
passes.

---

## 2. Kang 140 °C: the switch-on did NOT emerge, and the failure splits cleanly in two

The declaration pre-registered a 3.8× / 2.5× **under**-prediction for a
single-Ea branch. B2.1 under-predicts by **1500× (MFT)** and **6.6 million×
(FFT)** — far worse than a single-Ea model would have done. That is not a
subtle miss and it needs a mechanism, not a shrug. There are two, and they are
independent.

### 2a. The endpoint collapses because the pot is precursor-exhausted long before 120 min

Peak versus endpoint in the model's own 140 °C run:

| | model **peak** | model **at 120 min** | measured at 120 min |
|---|---:|---:|---:|
| MFT | 1.79e-05 mmol/L | 3.40e-08 | 5.17e-05 |
| FFT | 1.28e-04 mmol/L | 1.52e-11 | 1.00e-04 |

**The model's PEAK FFT is 1.28× the measured 120-minute value and its peak MFT
is 0.35× it.** The formation side is close to right; what is wrong is that
everything then gets destroyed. TTCA is fully consumed, free cysteine is down
to 2e-4 mmol/L and the sulfide to 1e-6 mmol/L well before the measurement is
taken, so the thiols spend the back half of the hold decaying into an
unreplenished pool.

**This is a cost of the very change that fixed the low-temperature rows.**
B2.1 gave the unassigned decay lumps the lumped formation activation energy,
because pinning them to 145 °C and evaluating them unchanged at 50 °C was a
strong claim with no evidence behind it — and that single change is most of why
Hofmann's 80 °C brew went from 18× to 2.9×, why Cerny's 95 °C ceiling now
passes, and why the route mix now moves with temperature. The same change makes
the sinks 4× faster at 140 °C than at 120 °C, which is fine while precursor is
still arriving and fatal once it is not.

**The honest statement is that B2.1 traded a low-temperature defect for a
high-temperature one, and the trade was worth making — eleven rows for one
block — but it was a trade, not a free win.** What it localises precisely: the
`*_decay` lumps should not share the FORMATION barrier. They need their own,
and Kang's ladder is the first dataset in the corpus that could identify one,
because it is the first to measure the same pot at three temperatures.

### 2b. The sulfide supply does not rise with temperature in this pot, so the gating mechanism never fires

The switch-on mechanism this wave relied on was that the thiol-forming flux is
`[carbon skeleton] × [sulfide]` with the sulfide gated behind Zheng & Ho's
measured 133 kJ/mol — the largest barrier in the module — so the thiol lane
should accelerate late. **It does not, and the model says why:**

| | 100 °C | 120 °C | 140 °C |
|---|---:|---:|---:|
| peak free cysteine, mmol/L | 1.23 | 1.03 | **0.064** |
| peak sulfide, mmol/L | 0.0050 | 0.0092 | **0.0075** |

The sulfide pool **peaks at 120 °C and falls at 140 °C**, because its source —
free cysteine — collapses. And it collapses through the coupling this module
has always had: the Nedvidek step `r_dpo_ddp` consumes cysteine **as a
reagent** for the carbon skeleton, and it carries the lumped 119 kJ/mol
barrier, so at 140 °C the skeleton lane burns the cysteine before the
133 kJ/mol thermolysis lane can turn it into sulfide. The two factors of the
two-factor law are coupled through a shared depleting reactant, exactly as
`PH_TERM_PROVENANCE` says — and here that coupling works against the model.

**⇒ A structural switch-on requires either a cysteine supply that outruns both
sinks at 140 °C (Kang's TTCA release, whose barrier is fitted and landed too
low) or a separate, steeper barrier on the thiol-forming step itself. Neither
is available from a declared FIT row. This is reported as an open defect with a
named next measurement, not fixed after the fact.**

---

## 3. Branch responsiveness: the ceiling is arithmetic, not a missing nonlinearity

`cerny2007_branch_responsiveness` measures a 2× precursor change moving the
xylose share of MFT from 15 % to 46 %, i.e. 3.067×. B2 predicted 1.49×; B2.1
predicts **1.13×**, which is worse. B2 diagnosed "a magnitude error in the
relative reaction ORDERS". That is right, and it can now be made exact.

The model's own shares are **0.412 at 1× and 0.466 at 2×**. Because a share
cannot exceed 1, the largest responsiveness this model could produce is
`1 / 0.412 = 2.43×` — **below the measured 3.067× no matter what any rate
constant does.** Reproducing the row REQUIRES a 1× xylose share at or below
0.326, and the model's is 0.412. The companion row says the same thing
directly: `cerny2007_1x_xylose_share` is a FAIL at 2.74× (predicted 0.41
against a measured 0.15) while `cerny2007_2x_xylose_share` PASSES at 1.01×.

**So the row is not failing through a missing nonlinearity. It is failing
through the 1× LEVEL**: the sugar route is too strong at low loading, which
caps the achievable response by simple arithmetic.

A removal scan over every shared pool that could plausibly damp the response
confirms that no coupling is responsible (each candidate switched off one at a
time, in the diagnostic only, never in the model):

| coupling removed | responsiveness |
|---|---:|
| *(frozen model)* | **1.133×** |
| sulfide sink `k_h2s_loss` | 1.123× |
| cysteine thermal sink `k_cys_thermal` | 1.125× |
| thiazole drain `k_cys_actz` | 1.132× |
| **Nedvidek cysteine-reagent step `k_dpo_ddp`** | **1.305×** |
| osone decay | 1.085× |
| MFT decay / dimerisation / thiolate loss | 1.133× (unchanged) |

Only the Nedvidek step moves it at all, and only from 1.13× to 1.31× — the
shared-cysteine coupling damps the response slightly, which is a real finding,
but it is nowhere near the gap. **No ad-hoc nonlinearity was added, and none
would have helped.**

---

## 4. One row got materially worse, and it is not hidden

`zhou_dimer_share_pH_invariance` went from 17.8× to **187×**. Zhou measures the
dimer's share of MFT-equivalents as nearly pH-invariant (8.6 / 6.5 / 9.6 % over
pH 6–8, a 1.48× spread) while [MFT] swings 3×. B2.1 puts a **thiolate** pH
factor on the dimerisation channels, because thiol oxidation to a disulfide
proceeds through the thiolate and because Kumazawa measures the disulfide
growing monotonically with pH in the same tubes as the thiol collapses. That
made the three individual dimer-share rows much better — pH 7 went 0.21× →
0.75× and pH 8 went 0.079× → 1.12×, both now PASSING — at the cost of the
invariance row, which now has too much spread across the three.

Both cannot currently be satisfied, and the tension is real rather than
numerical: Kumazawa's buffered system says the disulfide branch is strongly
pH-dependent, and Zhou's unbuffered system says the dimer SHARE is not. The
most likely reconciliation is that Zhou's three runs all drift to within 0.2 pH
units of each other by the end of heating, so their *effective* pH span is far
smaller than their labels suggest — which the module's pH-trajectory state
already carries, evidently not strongly enough. **Recorded as an unresolved
conflict between two FIT/hold-out sources, not adjudicated.**

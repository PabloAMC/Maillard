# Pre-registration: wave B39, the fed 3-deoxyglucosone fit (written 2026-09-11, BEFORE any fit was run)

## 1. Why

Three independent measurements now say the model's 3-deoxyglucosone limb is wrong in the same
direction: 3,4-dideoxyglucosone 32× low from glucose at 121 °C (B34, Leitzen 2021), 7–10× low in
ratio from glucose at 90–110 °C (B36, Zhang 2021), and 5× low **and three times too early** from
pure 3-deoxyglucosone at 120 °C (B37, Mittelmaier 2011). B37's probe refuted its own prediction:
the model does not so much under-make the intermediate as **destroy the whole 3-deoxy pool far too
fast**, and B38 said the constants involved have never been fitted at all — they are Kocadagli's
glass values carried with printed bands. B38 also said what pins a constant: **a fed pot on a small
network** (B20, B21 fully pinned) rather than an end-of-cook level through a large one. This wave
is that design on Mittelmaier's pot.

## 2. What the source gives, and what the network lacks

Mittelmaier et al. 2011 charge **~200 µM of one pure compound** in a glucose-free peritoneal-dialysis
model at pH 5, 0.5 mL in a sealed vial, **120 °C**, sampled at 0/10/20/30/60/90/120 min. Printed
numbers (`mittelmaier2010_extraction.md` §3):

| fed | observed | value |
|---|---|---|
| 3-DG | 3,4-DGE maximum | **26.7 µM at 30 min** |
| 3-DG | share of the other diastereomer (3-DGal) at 60 min | **26 %** |
| 3-DGal | 3,4-DGE maximum | **46.2 µM at 20 min** |
| 3-DGal | share of the other diastereomer (3-DG) at 60 min | **48 %** |
| 3,4-DGE | 3-DG maximum | **26.9 µM at 30 min** |
| 3,4-DGE | 3-DGal maximum | **37.9 µM at 20 min** |

And a mechanism, proved by the fed experiments: **3-DG ⇌ 3,4-DGE ⇌ 3-DGal**, reversible
dehydration and hydration. The trunk today has `r_tdg_ddg` one way, `r_ddg_hmf` as the only exit
from 3,4-DGE, and **no 3-deoxygalactosone at all.** At 120 °C the model's 3-DG leaves through three
exits — formic acid (Martins, 0.057 /min), methylglyoxal (Kocadagli, 0.010 /min) and 3,4-DGE
(0.007 /min) — with a 9.5-minute half-life, and its 3,4-DGE leaves to HMF at 0.119 /min on a
barrier the authors fixed to zero.

Second source, already on disk: Zhang et al. 2021 (`zhang2020_extraction.md`), 0.3 M glucose in
water at 90/95/100/105/110 °C for 6 h, printing 3-DG and 3,4-DDG slopes whose **within-study ratio**
is unit-free (their absolute levels do not mass-balance and are not used). Under the owner's rule
fed yields and within-study ratios FIT.

## 3. The structural change, declared before the fit

Three reactions and one species enter the trunk, inert until fitted:

| reaction | what | barrier |
|---|---|---|
| `r_ddg_tdg` | 3,4-DGE + H₂O → 3-DG (the reverse hydration the paper proves) | declared equal to the forward step's 36.9 kJ/mol; flagged `barrier_declared_equal_to_forward` |
| `r_ddg_dgal` | 3,4-DGE + H₂O → 3-DGal (the C4 epimer) | same declaration |
| `r_dgal_ddg` | 3-DGal → 3,4-DGE | same declaration |
| species `DGAL` | 3-deoxygalactosone, C6, intermediate, molar mass 162.14 | measured in the fit corpus (this paper) |

The declaration of one shared barrier is the honest minimum: the paper measures at one temperature,
so it cannot give a barrier, and inventing four would be four fabricated numbers. Zhang's five
temperatures then test the transfer.

## 4. The fit

**Coordinates (5), all log10 k at the trunk's 100 °C reference:** `k_tdg_ddg`, `k_ddg_tdg`,
`k_ddg_dgal`, `k_dgal_ddg`, `k_ddg_hmf`. The last is freed because its shipped value is a 160–200 °C
glass rate with a barrier fixed to zero, i.e. a bracket, and B38 classed such constants as pure band.
**Not freed:** `k_tdg_fa` and `k_tdg_mgo` (measured on their own systems and, for formic acid, a
B1 fit row) and `k_glc_tdg` (the glucose entry, which the fed pot does not touch).

**Rows (12):** the six Mittelmaier maxima and shares (σ = 0.15 dex on a maximum; 0.5 on a share,
because the paper does not say whether "percentage of the other diastereomer" is of the pair or of
the total, and the pair reading is taken), the **two times of maximum** (30 and 20 min, σ = 0.25 dex,
since the sampling grid brackets them between 20 and 60 / 10 and 30), and Zhang's four
temperature-resolved ratios at 95–110 °C (σ = 0.2 dex; the 90 °C row's 3,4-DDG slope is the
smallest printed and is reported, not decisive). Fed pots: 0.2 mM of the fed compound, water, pH 5,
120 °C, integrated to 120 min. Zhang pots: 300 mM glucose, water, pH 6.5, 6 h.

**Bounds:** ±2 dex around each start (the shipped value for the two existing constants; for the
three new ones, the forward step's shipped value, since the paper's forward and reverse maxima are
equal to within 1 %). Two starts, as B21.

## 5. Predictions

- **P1.** The fit reaches the four maxima within 0.3 dex and both times of maximum inside their
  brackets, with χ²_red below 3.
- **P2.** `k_tdg_ddg` rises by **0.5–1.2 dex** — less than the 32× (1.5 dex) that B34's Leitzen row
  alone would ask for, because the reverse step and the epimer take part of the load.
- **P3.** `k_ddg_hmf` **falls** by at least 0.5 dex: the 3,4-DGE pool must live longer than a
  6-minute half-life for the paper's maxima to exist at 20–30 min.
- **P4.** On Leitzen 2021, the hold-out this fit never reads, 3,4-DGE's fold error improves from
  32× to **under 10×**, 3-DG stays within 3×, and HMF does not worsen beyond its current 12×.
- **P5.** At least three of the five coordinates come out PINNED under B38's thresholds, and the
  reverse pair (`k_ddg_tdg`, `k_dgal_ddg`) is the likeliest to be collinear.
- **P6.** The panel headline does not fall: no row now within 3× leaves the band.

## 6. Ship rule

SHIP if P1, P4 and P6 hold. P2, P3 and P5 are reported. If P4 fails while P1 holds, the fit
describes the fed pot and not the glucose pot, and the wave records that as its finding and
**does not ship**: the point of a hold-out is to be allowed to say no. Both sides of every
comparison are frozen artifacts under `results/validation/_b39_baseline/`.

If it ships: the five values become frozen literals in `parameters_dicarbonyl.py` asserted equal to
the report by a test; the envelope gains a `b39.` prior row per coordinate from the fit's own
Laplace σ (the B18 pattern), and the ENV-B34 printed band on `k_tdg_ddg` is retired in favour of the
data-derived width, with the reason recorded.

## 7. Outcome (written 2026-09-11, after the run) — DO NOT SHIP, by the rule as written

**The fit is good and the rule says no.** Artifacts: `kinetic_core_b39_fit_report.{json,md}`,
`kinetic_core_b39_ship_rule.{json,md}`, the frozen pair under `_b39_baseline/` (the AFTER panel is
the candidate scored in a process that installed it for itself only; nothing was installed in the
tree).

| prediction | result |
|---|---|
| P1 rows fit | **REFUTED on one clause.** All four maxima within 0.09 dex, all four Zhang ratios within 0.05 dex, χ²_red 0.58 — but the fed-3-DG pot's 3,4-DGE peaks at **12.5 min**, outside the 20–60 min bracket the paper's sampling grid puts around its 30-minute maximum. |
| P2 `k_tdg_ddg` up 0.5–1.2 dex | HELD: +0.57 dex |
| P3 `k_ddg_hmf` down ≥ 0.5 dex | HELD: −0.60 dex |
| P4 Leitzen hold-out | **REFUTED by 2 %.** 3,4-DGE 32.4× → **8.1×** (the rule asked for under 10×), 3-DG 1.11× → 1.23×, but HMF 11.93× → 12.21×, and the rule said HMF must not worsen. |
| P5 ≥ 3 pinned | HELD: **all five pinned**, σ 0.08–0.21 dex, no collinear pair. The fed design pins what the level-scored fits could not, as B38 said it would. |
| P6 no row leaves the band | HELD — and two rows **enter** it: within-3× would go 10/45 → 12/45, out-of-sample 9/44 → 11/44. |

**What the failure says.** With the five triangle constants free, the fit cannot delay the peak,
because the peak's timing is set by how fast 3-DG *leaves* — and 3-DG's dominant exit, the
formic-acid step, is Martins' `k_tdg_fa` at 0.056 /min (120 °C), a B1 fit row that this wave did
not free. That constant was measured at **pH 6.8**; the fed pot is at **pH 5**. And the same
laboratory prints what pH does to it: Martins & van Boekel 2003, Part II, Table 3 — the step
3-DG → formic acid is **14× slower at pH 5.5 than at 6.8 at 100 °C, and 7× slower at 120 °C**;
the other 3-DG exit (3-DG → fragments) is 6.6× slower at 100 °C. The trunk's pH term (B12, from the
same table) scales only the three Amadori-decay steps. **The 3-DG exits have never had one.** That is
wave B40, pre-registered before anything else moves, and B39's five values are its starting point,
not its result.

**Not installed:** `SHIPPED_B39` stays False, the three new steps stay at k = 0, and the two existing
constants keep their literals. The candidate values are in the report and nowhere else.

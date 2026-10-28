# Wave B63 pre-registration — single-temperature barriers at their measured class's median

*Written 2026-10-28, before any B63 number existed. Pablo asked "Can you fix it?" about the incoherence below.
Accepted under the frozen rule (`scripts/generators/WAVES.md`, "Evidence decides, scores diagnose"): on the
evidence audit in sec. 3, with V1 as its mechanical gate. T1-T4 are recorded as diagnosis. An independent
reviewer reads the audit before the engine switches. Generators: `generate_kinetic_core_b63_fit.py`,
`generate_kinetic_core_b63_ship_rule.py`. Switch: `parameters_sulfur.SINGLE_TEMPERATURE_BARRIERS`. Baseline:
what ships when the fit runs (B61 if accepted, else B60). It runs BEFORE B62, whose baseline then includes it.
North-star targets served: 1 (measured where the data allow) and 2 (the point and the interval come from one
model).*

## 1. Why

31 sulfur constants are known at 145 °C only. 30 of them share the lumped route barrier, 64.08 kJ/mol, and
`k_mmft` has none. The fit never identified that barrier. B2.3 left it unidentified, and B10's attempt to split
it by route returned Laplace σ 74 and 52 kJ/mol. Since B55 (2026-10-28) the ENVELOPE draws each of these
barriers from the empirical distribution of the corpus's measured barriers of its class: carbohydrate, n = 20,
median 103.7, range 30-237; sulfur, n = 15, median 100.9, range 30-133. The POINT kept 64.08. Off 145 °C the
point therefore sits in its own prior's lower tail. The independent review of 2026-10-28 found Yiltirak's and
Zhou's FFT points 1.0-1.4 dex above their envelope medians: the tool reports two different models.

The 30 constants carry most of the flux. `r_pent_tdp`, `r_pent_dpo` and `r_dpo_ddp` alone are 74 % of the mean
reference-pot flux (`SULFUR_FLUX_MAP.md` sec. 5c).

## 2. What changes

- Each route-sharing single-temperature constant takes its class's median barrier: `sugar_trunk` → 103.7,
  `thiol_assembly` → 100.9 kJ/mol. Its 145 °C rate is held where the fit put it.
- Untouched: measured overrides (they take precedence, so B62's furfural barriers win when B62 ships), the two
  decay families' fitted barriers, `k_mmft` (no barrier, policy 2), and every constant not on the lumped barrier.
- The envelope's B55 rows keep their class distributions; their centre becomes the class median, the point.
- **Candidates:** (A) refit (kind 2): the shipped objective, same free set and protocol, start 0 = the shipped
  optimum; (B) the shipped vector with the switch on (recorded).

## 3. Evidence audit

| audit item | finding |
| --- | --- |
| 1. established for the regime | The class distributions are the corpus's measured barriers for reactions of the same class (sugar-side steps with no sulfur reactant; steps joining a sulfur nucleophile to a carbonyl), aqueous, 75-180 °C, read from the registries (`ea_class_prior.py`). They are a class prior, not step measurements: the strongest evidence available for steps nobody has measured, and the evidence B55 already accepted for the envelope. |
| 2. carries flux | The 30 include the three steps that carry 74 % of the mean reference-pot flux. |
| 3. constants measured, bands carried | Medians of measured distributions; the bands are the class ranges, already drawn by the envelope. No new free parameter: the 145 °C constants stay the fit's coordinates. |
| 4. verified in the engine | **V1**: with the switch on, every route-sharing target carries its class median and an unchanged k(145 °C), and no other barrier moves. Checked when the generator was written: PASS, exactly 30 constants moved. |
| 5. no contradiction with a direct measurement | The step-specific measurements found on 2026-10-28 bracket the medians: the Amadori analogue's loss 75 ± 10 (Deng 2021), norfuraneol formation 81.5 and 93 (Lau 2003, Pandit 2006), furfural formation 147 (Schirle-Keller 1992), TTCA loss over 100-140 °C (Zhai 2021, Kang 2026, in the fit). None sits at 64; the median of the four is ~87 and the range 75-147. No measurement contradicts ~100 more than it contradicts 64. Declared caveat: the class median is not a step measurement, so each constant's own measurement, when it comes, overrides it. |

**Audit verdict, before any result: grounded as a class prior. It replaces an unidentified number with the
median of measured ones and makes the point and the interval one model.**

## 4. Recorded diagnosis (not gating)

T1 Yiltirak's FFT rows; T2 the hold-out sulfur rows; T3 the Hofmann ribose + xylose pH 5 rows; T4 the
directional panel; the fit cost of each arm; Kang's furfural 120/100 °C ratio; the typical-yield baseline and the
grounding share after the wave. A new coherence check: the median over off-145 °C sulfur rows of |log10(point /
envelope median)|, before and after.

**Expectations (mine, before any number):** the off-145 °C trunk slows (×0.25 at 100 °C, ×0.5 at 120 °C for the
sugar-side steps at fixed k(145)). The refit will raise some 145 °C constants to defend the 100-120 °C fit rows
(Kang, Zhai, Kumazawa 121 °C), so the net change at 100 °C will be smaller than ×4. Yiltirak FFT (234x high at
100 °C) should fall; MFT likewise, and it is already within 1.5-7x, so it may go low. The point-vs-median gap
should shrink to < 0.3 dex on the FFT rows.

## 5. Verdict

*(appended after the rule runs; nothing above this line is edited after that)*

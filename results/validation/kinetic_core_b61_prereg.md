# Wave B61 pre-registration — norfuraneol's hydroxide-catalysed decay, measured (Mikami 2017)

*Written 2026-10-28, before any B61 number existed. Approved by Pablo the same day ("Go ahead with the Mikami
norfuraneol decay wave"). Kind 3 under `scripts/generators/WAVES.md` ("How a change is accepted"): a new
reaction, so it is accepted on pre-registered hold-out tests. Generators:
`scripts/generators/generate_kinetic_core_b61_fit.py`, `scripts/generators/generate_kinetic_core_b61_ship_rule.py`.
Switch: `parameters_sulfur.NF_HYDROXIDE_DECAY` (`enabled`, off until a SHIP verdict). Baseline: what ships (B60).*

## 1. Why

Norfuraneol (NF, 4-hydroxy-5-methyl-3(2H)-furanone) is MFT's only carbon source in the shipped network
(`docs/guides/SULFUR_FLUX_MAP.md` sec. 6: the DDP route carries nothing). The engine loses NF by `r_nf_decay`,
first order, **no pH factor**, log10 k(145 °C) = −2.18 min⁻¹ on the carbonyl-sink barrier (143.5 kJ/mol). That
value is not a measurement: B2.3 fitted it (B2 −3.59, B2.1 −8.63, B2.2 −4.56, B2.3 −2.18) and every wave since
has frozen it. At 100 °C it gives 4.5 × 10⁻⁵ min⁻¹ (t½ ≈ 250 h).

**Mikami et al. 2017** (`data/lit/extraction_dossiers/mikami2017_extraction.md`; Food Sci. Technol. Res. 23:283,
10.3136/fstr.23.283) heats pure NF (1.7–1.8 mM) alone in 0.2 M Na/K phosphate at pH 8.0, 100 °C, and follows it
by HPLC. It is the only NF-alone decay in the corpus:

- pH 8, 100 °C: first order, **k = 0.032 min⁻¹** (t½ 22 min; ln-linear over 0 / 0.5 / 1 h, pointwise 0.034 and
  0.032; figure-read, ±10 %). That is **~710×** the engine's value at the same temperature.
- pH 7, 100 °C (xylose + lysine pot): NF falls while still being formed, so its decay is **≥ 1.8 × 10⁻³ min⁻¹**.
- pH 6, 100 °C: no net loss to 2 h, with the maximum at ~4 h. Decay is much slower there.
- NF fragments to diacetyl (≥ 5.5 mol %) and methylglyoxal (≥ 0.6 mol %), transient lower bounds; most NF carbon
  goes to colourless polymers and pigments.

Mikami's observations are what **specific-base (hydroxide) catalysis** predicts, one decade per pH unit: anchored
at 0.032 min⁻¹ at pH 8, the law gives 3.2 × 10⁻³ at pH 7 (≥ 1.8 × 10⁻³ ✓) and 3.2 × 10⁻⁴ at pH 6 (~4 % loss in
2 h, inside "no net loss while forming" ✓). 4-Hydroxy-3(2H)-furanones are enols. Base-catalysed ring opening
through the enolate is the textbook route, so the mechanism claim is ordinary. The engine's `"base"` factor is
not this: it is the cysteine free-amine fraction (amine catalysis of 2,3-enolisation).

Where it should matter: every pot at pH ≥ 7. The blind Zhou 2023c hold-outs (Ala-xylose ARP + glutathione,
unbuffered, initial pH 7.5, 120 °C) over-predict MFT 14×, 27×, 29× and 161× at 30, 60, 90 and 120 min. The measured
MFT plateaus and then collapses (59 → 64 → 57 → 8 µg/L), which a fast NF sink would produce. Where it should hurt:
Hofmann's pH 7 arm (MFT already 4 900× low from H2S starvation) and Zhang 2025b at pH 6.5 (MFT 44× low).

## 2. What changes

- **New reaction `r_nf_decay_oh`: NF → FRAG_C**, first order, pH factor `"hydroxide"` = 10^(pH − 5), normalised
  to 1 at pH 5 at the same temperature, as every engine pH factor is. It runs in parallel with the existing
  pH-independent `r_nf_decay`, which is unchanged and stays frozen. `ph_factor` gains the kind `"hydroxide"`.
  No existing reaction is retagged.
- **Its rate constant is MEASURED, not fitted:** k_OH(100 °C, pH 8) = 0.032 min⁻¹ (Mikami), so at the engine's
  pH 5 normalisation log10 k_OH(100 °C) = log10(0.032) − 3 = **−4.495**.
- **Declared assumption, its barrier:** never measured. Following B55's convention for a single-temperature
  constant, the point value is the **median of the corpus's measured carbohydrate-class barriers, 103.7 kJ/mol**
  (`src/kinetic_core/ea_class_prior.py`, n = 20, range 30–237). The constant is anchored at **100 °C**, the
  measurement temperature, so log10 k_OH(145 °C, pH 5) = −4.495 + 103.7 / (R ln10) × (1/373.15 − 1/418.15) =
  **−2.933**. The envelope draws this barrier from the carbohydrate class's empirical distribution with
  **k(100 °C) held**, the mirror of B55's k(145 °C) hold. That changes no point prediction.
- **Declared assumption, the pH law:** one decade per pH unit across the whole pH range, measured at one pH. At
  pH 5 the new term is 0.18× the neutral `r_nf_decay` at 145 °C and 0.7× at 100 °C (the barriers differ). Each pH
  unit lower divides it by 10, so it is negligible at pH ≤ 4 and a modest addition at pH 5–5.5. At pH 7 it is
  18× the neutral term at 145 °C.
- **Products:** FRAG_C, like the existing decay. Mikami's diacetyl and methylglyoxal are lower bounds on transients.
  Routing ≥ 0.6 % to methylglyoxal (which `r_mgo_mp` and `r_cys_actz` draw on) is a separate structural question and
  NOT part of this wave.
- **pH used:** the in-situ pH trajectory the engine already integrates (B60's converged fixed point), as for every
  pH factor.
- **Candidates, fixed order of preference:** (A) **refit**: B60's 52-row objective with the step on, the same
  23-coordinate free set, bands, sigmas, optimiser, budget, quick setting and two-start protocol (start 0 = B60's
  optimum); (B) **B60's frozen vector with the step on**. Neither touches k_OH, its barrier or the neutral
  `k_nf_decay`.

## 3. Ship rule (`generate_kinetic_core_b61_ship_rule.py`), each arm scored in a child process against B60 on the same numerics

- **T1 (the measurement is what the engine does, both must hold):** NF fed alone at 1.75 mM in 0.2 M phosphate,
  pH 8.0, 100 °C, 60 min: the engine's effective first-order NF loss constant over 0–60 min is within 0.1 dex of
  0.032 min⁻¹. The same pot at pH 7.0 gives ≥ 1.8 × 10⁻³ min⁻¹. This verifies the implementation. A failure is a
  bug, not a verdict, and is fixed before T2–T5 are read.
- **T2 (the hypothesis):** the four Zhou 2023c MFT rows (blind, Amendment 45): median |log10 error| falls by
  **≥ 0.3 dex** against B60.
- **T3:** the Maillard-path hold-out sulfur rows: within-3x not lower than B60's by more than 2, and median
  |log10 error| not worse by more than max(0.1 dex, one bootstrap standard error of B60's median), as in B58's T2.
- **T4:** the Hofmann ribose pH 5 + xylose pH 5 hold-out rows (MFT, FFT): median |log10 error| not worse than
  B60's by more than the same tolerance.
- **T5:** directional panel, independent claims agreeing: not fewer than B60's minus 3.
- **SHIP (A)** if (A) passes T1–T5, else **(B)** if (B) does, else **DO NOT SHIP** (the step stays off and the
  measurement is recorded as a conflict).
- **Recorded, not gating:** (i) the Hofmann pH 3 / 5 / 7 sweep shape (MFT and FFT ratios pH 3 : 5 : 7); (ii) Zhang
  2025b MFT and Yiltirak's four MFT rows; (iii) the fit cost of each arm against B60's 27.132, and whether (A)'s
  rise exceeds the 3.59 admissible-slice threshold (`kinetic_core_b60_profile.json`); (iv) T2 re-evaluated with
  the barrier at the carbohydrate class's 10th and 90th percentiles (B60's vector, step on), to show how much rides
  on the unmeasured barrier.

## 4. Expectations

Expected moves: Zhou 2023c MFT falls, by most at 90–120 min where the in-situ pH has had longest above 7. It
cannot fall far once the unbuffered pot acidifies (Zhou 2023's own pH 7 pot ends at pH 3.42), so the 161× row
may stay high. Hofmann pH 7 MFT falls further (a worse row, but already 3.7 dex off for a different reason).
Zhang 2025b MFT falls (worse). The pH ≤ 5.5 rows (Hofmann pH 3/5, Yiltirak) move < 0.2 dex. In the refit, k_nf_mft
and k_dpo_nf should rise to defend the pH 7–8 fit rows (Zhou 2023's ARP + Cys pots at pH 7 and 8), and that would
lift the pH 5 rows somewhat.

Probabilities (mine, before any number): T1 ≥ 0.95 (arithmetic). T2 ≈ 0.65. T3 ≈ 0.6 (the rows it hurts already
miss, so the within-3x count is robust, but the median can move). T4 ≈ 0.75. T5 ≈ 0.8. **SHIP (A or B) ≈ 0.4.**

## 4a. Amendment 1 (2026-10-28, while writing T1's probe, before any B61 fit or score existed)

Mikami's "pH 8.0" is the phosphate buffer's label. The engine converts every buffer to its in-situ pH at
reaction temperature, and for Mikami's pot (0.2 M phosphate, label 8.0, 100 °C, 1.75 mM NF) it computes
**8.082**, constant over the run. Anchored at the label, the engine reproduced 0.0387 min⁻¹ for the measured pot:
+0.08 dex, inside T1's 0.1 dex, but a systematic offset rather than noise. The anchor is therefore the engine's
own in-situ pH for the measured pot: log10 k_OH(100 °C, pH 5) = log10(0.032) − 3.082 = **−4.577**, and
log10 k_OH(145 °C, pH 5) = **−3.015**. Every other number in sec. 2 stands. The pH 7 check gives 3.9 × 10⁻³ min⁻¹
(≥ 1.8 × 10⁻³ ✓). The probe that found this ran the engine only on Mikami's NF-alone pot, never on a panel row.

## 5. Verdict

*(appended after the ship rule runs; nothing above this line is edited after that)*

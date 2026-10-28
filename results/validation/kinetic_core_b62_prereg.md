# Wave B62 pre-registration — furfural's own barriers: formation steep, loss shallow

*Written 2026-10-28, before any B62 number existed. Approved by Pablo the same day ("Go ahead with the furfural
barrier swap wave"), and re-scoped before writing (see sec. 1, "What the first proposal got wrong"). Kind 3 under
`scripts/generators/WAVES.md`: a structural change to which barrier two steps carry, accepted on pre-registered
hold-outs. Generators: `scripts/generators/generate_kinetic_core_b62_fit.py`,
`scripts/generators/generate_kinetic_core_b62_ship_rule.py`. Switch: `parameters_sulfur.FURFURAL_BARRIERS`
(`enabled`, off until a SHIP verdict). Baseline: what ships when the fit runs (B61 if its verdict is SHIP, else
B60), read from the engine. North-star targets served (`tasks/data_restructure_plan.md` sec. 7): 1 (measured
where the data allow) and 4 (timescales and temperature slopes, not only endpoints).*

## 1. Why

**The model's furfural runs the wrong way with temperature.** Furfural forms on the lumped sugar-trunk barrier
(64 kJ/mol) and is lost through `r_fur_decay` on the carbonyl-sink family's fitted 143.5 kJ/mol. A loss barrier
above the formation barrier makes furfural fall as temperature rises. Two measurements say the opposite:

- Kang 2026 (purified TTCA; fit rows): the model's furfural FALLS 2.5x from 100 to 120 °C; measured, it RISES 1.7x.
- Yiltirak 2026 (ribose + Cys, 0.5 M phosphate, pH 5.5; hold-outs): FFT, 95-99 % furfural-derived, is 234x high
  at 100 °C / 4 h, falling to 59x at 130 °C / 30 min. MFT on the same pots is within 1.5-7x.

`docs/guides/SULFUR_FLUX_MAP.md` sec. 5 names this. The flux attribution shows that every scored pot's furfural is
100 % `r_tdp_fur`, and 51-95 % of it leaves by `r_fur_decay` (`results/validation/sulfur_flux_attribution.md`).

**The measured barriers for the operative steps:**

- **Loss.** Yaghmur, Aserin, Abbas & Garti 2005 (`yaghmur2005_extraction.md`): furfural disappearance with
  cysteine, 0.5 M phosphate pH 5.0, 40-70 °C, first order, **Ea 46.50 ± 1.0 kJ/mol** (water arm). It is the one
  measurement of furfural loss in a cysteine system. Over 98 % of that flux is not FFT, which is what
  `r_fur_decay` represents. Caveats: 40-70 °C against our 90-145 °C, a pseudo-first-order constant at ~0.6 M
  cysteine, and a lump. A consistency bridge, recorded and not used: Yaghmur's k_obs at 65 °C
  (2.1 × 10⁻³ min⁻¹) against the engine's fitted k(145 °C) implies ~77 kJ/mol at Yaghmur's cysteine, and ~115 at
  a 40 mM pot if the sink is first order in cysteine. Every reading puts the loss barrier below 143.5.
- **Formation.** Schirle-Keller & Reineccius 1992 (`schirle1992_extraction.md`; glucose 0.5 M + cysteine 0.1 M,
  water, unbuffered, 80-150 °C): net furfural accumulation, **Table 2 Ea 35.18 kcal/mol = 147 kJ/mol** (the
  dossier's refit of Table 1 gives 156, or 171-189 under another reading of the 150 °C rate). Raw: furfural rises
  ~30x from 100 to 120 °C, where 64 kJ/mol gives 2.9x. Caveats: a hexose rather than a pentose, net accumulation
  rather than a formation step, and arbitrary area units (so a shape, not a level).

**What the first proposal got wrong.** It cited Garrett & Dvorchik 1969 and Danon et al. 2014 for formation
(~128 kJ/mol) and Garrett's furfural solvolysis (89) for loss. Both measure amine-free, acid-catalysed chemistry:
~3 × 10⁻⁷ min⁻¹ for pentose → furfural and ~5 × 10⁻⁸ min⁻¹ for furfural solvolysis at pH 5 and 145 °C. Neither
carries flux in any scored pot (`tasks/lessons.md`, 2026-10-28). They are not used.

**Why only furfural's steps move.** Norfuraneol formation shares the same lumped 64 kJ/mol. Two measurements put it
at 81.5 ± 3.2 (Lau 2003) and 93.0 ± 6.4 kJ/mol (Pandit 2006), ribose + amines at 100-130 °C
(`lau2003_extraction.md`, `pandit2006_extraction.md`). Raising the shared barrier to 147 would push norfuraneol
55-65 kJ/mol past both. So the split is furfural-only. MFT's route is untouched by construction.

## 2. What changes

- **`k_tdp_fur`** (TDP → furfural; the furfural-forming step in every pentose pot) takes a measured barrier,
  **147 kJ/mol** (Schirle-Keller Table 2), instead of the lumped sugar-trunk 64.
- **`k_fur_decay`** (furfural → fragments) leaves the carbonyl-sink family and takes a measured barrier,
  **46.5 kJ/mol** (Yaghmur 2005, water arm). The family's fitted barrier still governs `k_osone_decay` and
  `k_nf_decay`.
- Both enter through `parameters_sulfur.measured_ea_overrides()`, so the fit cannot move them, exactly like
  `MEASURED_EA_OVERRIDES`. Their 145 °C constants stay free in the refit.
- Nothing else changes: no new species, no new step, no pH factor.
- **Candidates, fixed order of preference:** (A) **refit**: the shipped wave's 52-row objective with the switch
  on, the same free set, bands, sigmas, optimiser, budget, quick setting and two-start protocol (start 0 = the
  shipped optimum); (B) **the shipped vector with the switch on**, no refit.

## 3. Ship rule (`generate_kinetic_core_b62_ship_rule.py`), each arm in a child process against the shipped baseline

- **T1 (the hypothesis):** Yiltirak's four FFT hold-out rows (100 / 110 / 120 / 130 °C): median |log10 error| falls
  by **≥ 0.3 dex** against the baseline.
- **T2:** the Maillard-path hold-out sulfur rows: within-3x not lower than the baseline's by more than 2, and the
  median |log10 error| not worse by more than max(0.1 dex, one bootstrap standard error of the baseline median).
  The same rule as B58 and B61.
- **T3:** the Hofmann ribose + xylose pH 5 hold-out rows: median |log10 error| not worse by more than the same
  tolerance.
- **T4:** directional panel, independent claims agreeing: not fewer than the baseline's minus 3.
- **SHIP (A)** if (A) passes T1-T4, else **(B)** if (B) does, else **DO NOT SHIP**.
- **Recorded, not gating:**
  - (i) Kang's furfural 120/100 °C ratio (fit rows), against the measured 1.7x rise;
  - (ii) the Zhou 2023c furfural rows (blind; 25-45x low at 120 °C);
  - (iii) Resconi's furfural (150 °C);
  - (iv) Yiltirak's FFT temperature slope (log FFT ratio 130/100 °C), model against measured;
  - (v) each arm's fit cost against the baseline's;
  - (vi) the typical-yield baseline (`build_null_model_baselines.py`) re-run on the winning arm.
  - Schirle-Keller's own 30x shape cannot be scored: a glucose pot's furfural runs on `k_glc_fur`, which this
    wave does not touch, and the hexose entry is unidentified. Recorded as untested.

## 4. Expectations

Moves expected:

- Furfural at 100-110 °C falls (faster loss, slower formation at low T) and stays near the baseline at 145 °C.
- Yiltirak FFT falls most at 100 °C: the right direction, but it starts 59-234x high, so T1's 0.3 dex is reachable.
- Kang 100 °C MFT/FFT moves AWAY from its target. FFT falls, the measured ratio is 0.33 and the model's is
  already 0.97 (B60 slice).
- That is the fit's main resistance, as in the slice of the family barrier (+4 cost at 90 kJ/mol). In the refit,
  k_tdp_fur and k_fur_fft should rise to defend Kang's 100 °C FFT. That could undo part of the low-T fall, which
  is why (B) is also scored.
- Zhou's furfural (already too low at 120 °C) gets lower: a worse recorded row.

Probabilities (mine, before any number): T1 ≈ 0.6. T2 ≈ 0.6. T3 ≈ 0.75. T4 ≈ 0.8. **SHIP (A or B) ≈ 0.35.**

## 5. Verdict

*(appended after the ship rule runs; nothing above this line is edited after that)*

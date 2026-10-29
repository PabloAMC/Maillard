# Wave B59b pre-registration — de Roos's route: the TTCA reservoir and the cysteine Amadori compound

*Written 2026-10-29, before any number with the switch on existed. Approved as the second half of B59 (owner,
2026-10-28; D2: on the literature alone; D3: the shared trunk). Accepted under the frozen rule (`WAVES.md`,
"Evidence decides, scores diagnose") on the audit in sec. 3 and the mechanical rule in sec. 4. Switch:
`parameters_sulfur.CYS_AMADORI_ROUTE`. New species `CAMA`; new steps `r_cama_form`, `r_cama_dpo`, `r_cama_tdp`;
retired when on: `r_pent_dpo`, `r_pent_tdp`. North-star targets: 1 (measured constants carry the trunk) and 4
(timescales), and C4 (it REMOVES two fitted directions, `k_pent_dpo` and `k_pent_tdp`, and adds one anchored
constant).*

## 1. Why

The flux map (`SULFUR_FLUX_MAP.md` sec. 4) found the trunk running 10³-10⁴× too fast: the model's cysteine acts as
an instant catalyst on free pentose (`r_pent_dpo`, `r_pent_tdp`) and is destroyed in seconds, where measured pots
hold it for hours. On 2026-10-28/29 five laboratories' data bear on it: Cao 2017 and Zhao 2017 (90 °C: the cysteine
pool Cys + TTCA + Cys-Amadori conserved at 0.44-0.47 M over 5 h; all cysteine-destroying channels together
≤ ~2×10⁻⁴ min⁻¹), Li 2011 (glucose + Cys: cysteine first order, t½ 10.5 h at 100 °C to 34 min at 130 °C, Ea 121
kJ/mol), Schieberle 2000 (MFT still rising at 12 h, 100 °C), and US 10,076,126 (MFT lagging then accelerating over
4 h at 110 °C, which the shipped model front-loads 50×). The chemistry that explains all five is de Roos 2005's:
cysteine + pentose ⇌ thiazolidine (TTCA, fast, K measured), the open Schiff base → cysteine Amadori compound
(slow), and its enolisation → deoxyosones with the cysteine released.

## 2. What changes (when the switch is on)

- `r_ttca_form` on, by detailed balance with Cao's K = 0.060 L/mmol (band 0.045-0.075) against the fitted
  `r_ttca_cys` (B58's step).
- `r_cama_form`: TTCA → CAMA, first order in TTCA. Rate at 90 °C ANCHORED on Cao 2017's blank pot (sec. 3 item 4),
  barrier 100.9 kJ/mol (Zhang 2026, Cys + glucose → Cys-Amadori; a declared hexose → pentose transfer).
- `r_cama_dpo` / `r_cama_tdp`: CAMA → DPO / TDP + free Cys. Total rate de Roos 2005 (4.5×10⁻³ min⁻¹ at pH 5,
  1.4×10⁻² at pH 6, ~100 °C: slope 0.493 decades per pH unit), barrier 85.7 kJ/mol (Zhang 2026); DPO share 0.30
  (Zhai 2019: 3-DX / 1-DX = 2.1-2.6, i.e. 0.28-0.32, both ends carried as alternatives, register CR-12).
- `r_pent_dpo` / `r_pent_tdp` set to zero (retired; their constants stay in the free set and become flat, which the
  profile will show). The amine-independent pentose steps stay.
- B59a-r (the Strecker sulfur) is run with and without, because its own verdict is pending.

## 3. Evidence audit

| item | finding |
| --- | --- |
| 1. established for the regime | de Roos 1992/2005 (cysteine + pentose in water/phosphate at 90-100 °C), Cao 2017 and Zhao 2017 (same, 90 °C, all four species measured), Zhai 2019 (70-130 °C) show the route and its intermediates in exactly these pots. |
| 2. carries flux | It replaces the steps carrying 74 % of the mean reference-pot flux (`r_pent_dpo`, `r_pent_tdp`, with `r_dpo_ddp` downstream). |
| 3. constants measured, bands carried | K (Cao), decay rate and pH slope (de Roos), both barriers (Zhang 2026), the split (Zhai 2019): measured. The one unmeasured constant, the formation rate, is anchored on its source pot, not fitted globally. Declared transfers: the barriers are glucose-system values; de Roos's decay is at ~100 °C in 0.5 M phosphate and is applied at the in-situ pH; CAMA is titrated as the calibrated ARP analogue. No new free coordinate in the global fit; two are removed. |
| 4. verified in the engine | V1: off, every existing test and number unchanged (checked: 83 tests). V2: the anchor — the engine, in its own description of Cao's blank pot, reproduces the Cys-Amadori series (0.006 → 0.133 mol/L, 0.5-5 h) with one constant; the TTCA and free-cysteine series and the conserved cysteine pool are then CHECKS, not fitted (each within the figure-read error, ± ~10 %, or recorded as missed). V3: de Roos's decay reproduced: a CAMA-only pot at pH 5 and 6, ~100 °C, decays at 4.5×10⁻³ / 1.4×10⁻² min⁻¹. |
| 5. no contradiction with a direct measurement | Kang 2026's fed-TTCA ladder (free cysteine rising then falling, ceiling 16.3 mol %) is fitted through `r_ttca_cys` / `r_ttca_deg`, which stay; Zhang 2025b's ~48 % cysteine loss in 50 min at 120 °C, pH 6.5 is a hold-out; Li 2011's glucose system is a hexose check, not a constraint. Any tension found is recorded, not tuned away. |

## 4. Mechanical ship rule (written now) and recorded diagnosis

**Ship when ALL hold:** V1-V3 pass; an independent reviewer finds no error in sec. 1-3; the refit is run on every
B64 reading (isothermal, ramp + hold, ramp only linear, Newtonian τ 4 and 6.7 min), with and without B59a-r, each
to the optimiser's convergence or 1 200 evaluations, and reported as converged or not. Scores do not gate.

**Pre-declared hold-outs (diagnosis):** US 10,076,126's MFT shape at 110 °C (0.008 / 0.067 / 0.31 / 1 at 0.5 / 1 / 2 /
4 h; the shipped model 0.42 / 0.80 / 1.05 / 1), Schieberle 2000's 100 °C series and its 100 °C-6 h / 145 °C factor
(7.9-13×), Li 2011's cysteine k(T) slopes (as a hexose check), Zhang 2025b's cysteine loss, the Hofmann pH 3 / 5 / 7
shape, every out-of-sample row, the contrast diagnostic and the typical-yield baseline. Expectation (mine, before any
number): the Nestlé shape gains a lag; the 100 °C factor rises toward 8-13 only if the released cysteine lives for
hours; the pH 3 arm moves only if the route's pH law, not the retired one's, governs it.

**Amendment 1 (2026-10-29, before any B59b hold-out number existed; only candidate B's cost, 121.5, had been
computed).** The diagnosis is one generator for B59a-r and B59b, `generate_kinetic_core_b59_holdout_diagnosis.py`.
(i) Zhang 2025b's blank is used as its full TTCA and Cys-Amadori series at 10-50 min (120 °C, pH 6.5, 40 + 40 mM:
TTCA 11.2 → 18.2, Cys-Amadori 0.07 → 1.61 mmol/L), not only its cysteine loss, whose fitted k the dossier flags as
disappearance from an HPLC peak with ~19 mmol/L unaccounted. (ii) US 10,076,126 does not print its buffer; the
sibling patent's 0.2 M phosphate is declared. (iii) Li 2011 is NOT evaluated: the sulfur lane carries no hexose
species, so the declared hexose check has no pot; recorded, not substituted.

**Amendment 2 (2026-10-29, after three all-free refits and one hold-out diagnosis; before any refit on the corrected
free set).** *What went wrong.* The first refits freed all 23 coordinates. The frozen rule says "Kind 2 refits ship
with the change they follow, refitting only coordinates the data identify" (`WAVES.md`); sec. 4 did not restate it,
and I did not apply it. The all-free optimum (isothermal, cost 35.99, converged) moved coordinates that B60's Laplace
leaves unidentified: the carbonyl-sink barrier from 143.5 to 20.0 kJ/mol, `k_thiol_decay` by +5.4 dex, `k_dimer_fft`
by −6.4 dex, `k_fur_fft` by +2.1 dex. Its hold-outs collapsed: Schieberle factor 0.08 (B60 0.37; measured 7.9-13);
Nestlé MFT at 30 min ÷ 4 h 6.0 (measured 0.008); Zhang 2025b TTCA −1.4 dex at 50 min. *Correction.* The refit frees
only the 16 coordinates B60's Laplace marks identified, by its own rule
(`kinetic_core_b60_laplace_covariance.json`). The other 7 (`k_glc_ha`, `k_dimer_fft`,
`ph_acid_yield_per_sink_event`, `k_fur_fft`, `k_dimer_decay`, `k_thiol_decay`, `Ea_decay_carbonyl_sink`) stay at
the start vector's values. The same set serves every reading; identification is a property of the rows, and the
readings change only the 145 °C and 130 °C profiles (declared). The all-free members and their diagnosis are kept,
unused, in `withdrawn_all_free_set/`. *Arms.* B59a-r's two shares run ON TOP of the route
(`--strecker-share`). B59a-r on the old trunk is deferred: it is needed only if this wave fails. Corrected runs
per reading: route on, without the Strecker sulfur and at each share; 15 runs.

**Amendment 3 (2026-10-29): the corrected refits were stopped before any result; the audit in sec. 3 fails item 5.**
A probe at B60's constants with the route on and NO refit (scratch, recorded here) found:
- *The route does what it was adopted for.* Schieberle's 100 °C MFT over 30 min: 6.1 / 84 / 73 at 60 / 360 / 720 min
  (measured 3.1 / 35 / 40; route off 2.4 / 2.8 / 1.5). Nestlé's 110 °C MFT at 0.5 / 1 / 2 h over 4 h: 0.056 / 0.22 /
  0.57 (measured 0.008 / 0.067 / 0.31; off 0.43 / 0.81 / 1.05). Zhang 2025b's TTCA: +0.15 to −0.37 dex. Hofmann pH 7 MFT:
  −1.6 dex (off −3.7).
- *A direct, in-regime measurement missed in sec. 3 contradicts the formation step.* Zhai 2019 Fig. 3A (Xyl 165 +
  Cys 82.6 mM, water, pH 7.4 start, 70-130 °C, TTCA and Cys-Amadori over 0-160 min). Predicted Cys-Amadori at 160 min:
  3.1 / 20 / 50 / 38 mM, measured 2.8 / 3.7 / 3.8 / 0.0 at 70 / 90 / 110 / 130 °C. Per-temperature estimates in the
  engine (formation and decay multipliers fitted to the digitised series): formation must be 2-5× slower at 90-130 °C
  under either pH reading. The decay matches de Roos's law extrapolated to pH 7.4 if the pot's pH stays near 7.4
  (clamped). It must be 4-30× faster at 110-130 °C if the pH falls as the engine computes (to 3.7-4.4). Zhai never
  measured it. The transferred glucose barrier (100.9 kJ/mol) is too steep for pentose.
- *The fed-TTCA fit rows (Kang 2026 / Zhai 2022) made their MFT through the retired step.* With the route off,
  `r_pent_dpo` carries 9.4 of ~10 mM of flux in `kang_ttca_100`. With it on, MFT falls 260× and free cysteine peaks at
  63 % of the charge (ceiling 12 %). TTCA drains in 2 h, where Zhai 2021 measured 87 % left at 100 °C after 60 min;
  this is B16's diagnosis that `r_ttca_cys` (fitted on Kang's free cysteine) opens the ring ~10× too fast.
Consequence: a refit on B60's 52 rows would ask the optimiser to compensate for a mis-transferred barrier and a
mis-fitted ring opening. B59b as written does not ship (item 5 fails); its redesign constrains the route on the data
that measure its own species (Cao 2017, Zhai 2019 Figs. 3A-B, Zhai 2021, Kang's free cysteine, de Roos 2005) before
any flavour-level refit.

**Amendment 4 (2026-10-29): an implementation error in the decay's pH law, corrected (kind 1).** Sec. 2 cites
de Roos's decay at pH 5 and 6 and its slope. Fig. 4 measures the whole profile at pH 3-6 (% left after 2 h,
boiling 0.5 M phosphate): U-shaped, slowest at pH 5. Read as first-order constants: 5.6 / 5.3 / 5.1 / 4.8 / 4.5 /
5.1 / 7.0 / 14.3 × 10⁻³ 1/min at pH 3 / 3.5 / 4 / 4.5 / 5 / 5.25 / 5.5 / 6. The engine extended the 5 → 6 slope
downward, which made the decay 3.5× too slow at pH 4, where the acidifying Cao, Zhai and Kang pots spend most of
their time. `CysAmadoriRoute.decay_ph_table` now carries the measured profile, interpolated log-linearly. It is held
at the pH 3 value below 3. Above 6 it is either extended at the measured 5.5 → 6 slope or held flat (register
CR-14; 7× apart at pH 7.4). V3 still passes (4.50 / 14.3 × 10⁻³ at pH 5 / 6). V2 re-run: k_form 2.73 → 3.29 ×
10⁻³ 1/min at 90 °C; Cao's Cys-Amadori within ±0.09 dex at 1-5 h; TTCA −0.08 to −0.26 dex; free Cys within 0.05 dex;
final in-situ pH 3.94 (Cao measured 4.15, at room temperature). The route stays off; this wave's constants go to the
route calibration (`kinetic_core_cysteine_route_calibration_prereg.md`).

## 5. Verdict

*(appended after the rule runs; nothing above this line is edited after that)*

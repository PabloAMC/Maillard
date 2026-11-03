# Results — acrylamide cluster (CF-039, CF-104), 2026-11-03

Pre-registration: `prereg.md` (this folder). All arithmetic in `scratchpad/conflicts/acrylamide/tests.awk` →
`tests.out` (awk; no engine run was needed: every test is on printed or read-off values). Page reads done for the
values the verdicts rest on: Gökmen 2012 Fig. 2b (600 dpi) and Fig. 3 (300 dpi), Bassama 2011 Table 2 (200 dpi) and
water contents (text p. 1453), De Vleeschouwer 2008 Tables 1-3 (text layer, clean) and Figs. 2 / 4 (110-200 dpi).

## CF-039 — formation Ea 84 (Gökmen 2012) vs 154-182 (Leuven, Bassama)

**T2a (integrity) — Gökmen's "k" is a yield per hold time, not a rate constant.** Fig. 2b (ASN–GLC, 0.01 mmol each,
60 min) read 3.0e-6 / 4.08e-5 / 1.29e-4 / 1.58e-4 mmol at 90 / 120 / 150 / 180 °C; the text gives 1.62e-4 mmol "within
10 min" at 180 °C. ln[(AA / C0) / t_hold] with t_hold = 60 min at 90-150 °C and **10 min at 180 °C** gives −12.21 /
−9.60 / −8.44 / −6.43 against the plotted −12.14 / −9.59 / −8.45 / −6.41 (Fig. 3), and a slope of 84.5 kJ/mol against
the printed 83.94. Premise holds (P 0.9). Two consequences: the holds are mixed (60 vs 10 min); with 60 min at 180 °C
too the same estimator returns 60.9 kJ/mol. The paper says glucose and Asn are essentially gone after 5 min at 180 °C,
and the 150 °C 60-min yield is 80 % of the 180 °C plateau, so the 150 and 180 points are depletion-capped.

**T2b (outside the conflict) — the same estimator on Bassama's time courses**, whose rate-constant Ea_F is 168-178 on
the same data (Fig. 3 points, dossier sec. 4c): 140-180 °C, hold ~45-60 min → 66-69 kJ/mol; hold ~10 min → 78-133;
a_w 0.430 (peaks before 10 min) → 43 or negative. **The estimator can halve a 170 kJ/mol Ea.** Supports H2 in
principle.

**T2c — how much of Gökmen's gap the estimator explains.** Depletion-corrected (Y = Ymax(1 − e^{−kt}), Ymax =
1.62e-4): k = 3.1e-4 / 4.8e-3 / 2.65e-2 min⁻¹ at 90 / 120 / 150 °C → Ea 95.0 (90-150), 108.5 (90-120); the
low-conversion pair 90-120 uncorrected 103.3, with the 90 °C bar read ×1.5 / ×0.5: 87-131. **The estimator explains
84 → ~95-110, not 84 → 160** (rule: ≤ 115 = explains at most part).

**T1 — medium / water state (H1).** Single-response lumped Ea_F (class A), 10 rows outside Gökmen: Claeys aqueous
168.3; De Vleeschouwer 2007 powders a_w 0.34-0.92: 153.6-182.2; Bassama plantain a_w 0.43-0.97: 168.5-177.8. Range
153.6-182.2 (29 kJ/mol), mean 168.2, slope 17 kJ/mol per unit a_w. Aqueous, powder and food give the same Ea within
one estimand. **H1 falsified as the cause of an 80 kJ/mol gap** (rule: spread < 30 → falsified; P that medium explains
more than half the gap ≤ 0.1). Note the low values also come from an aqueous system: Knol 2005 k4 94.4 / k3 102
(0.2 M phosphate), so "dry silica vs aqueous" cannot be the axis.

**T3 — heat model (H3).** A 1.5 / 3 min come-up lag on Gökmen's holds moves the Ea 84.5 → 86.3 / 88.5. Small tubes:
H3 < 5 kJ/mol here (rule: ≥ 20 to matter). **T7 (outside, food scale)**: Gökmen 2006 whole-strip acrylamide at 9 min
(10 / 418 / 1422 ng/g) gives 104 kJ/mol between 170-190 °C referenced to oil T and 86 referenced to the measured
surface T (291 vs 197 for 150-170 °C). The sprague2011 lesson holds for foods (20-90 kJ/mol from the thermal
reference alone) but not for Gökmen 2012's 0.5 g tubes.

**What remains (H4).** After T2c the low cluster is 94-110 kJ/mol (Gökmen 2012 low-conversion 95-109; Knol k4 94.4,
k3 102; Gökmen & Şenyuva 2006 k′1 93.8, cited in Palazoğlu 2008 Eq. 7; De Vleeschouwer 2008 potato k_INT 105-126 and
2009 I k_INTg 117.5) against the high cluster 146-182 (lumped single-response k_F: Claeys, De Vleeschouwer 2007,
Bassama; and multiresponse k_F(Int1 → AA) 136-172 in most Leuven columns, with outliers 84.6 / 89.9). The split runs
along **which step the estimator weights** (condensation / early yield vs the slow Int1 → AA step and lumped fits
with co-fitted elimination), not along medium. P(estimand/step split) 0.5; P(genuine laboratory or method
difference) 0.3; P(Arrhenius curvature between 90-120 and 120-200 °C) 0.2. Decider: one laboratory, one pot,
acrylamide AND asparagine time courses at ≥ 5 temperatures spanning 100-180 °C, fitted both ways.

**Verdict CF-039: open, narrowed.** The 84 is explained down to ~100 by its estimator (P 0.75); the remaining
~100 vs ~165 is not a medium effect (P 0.9) and not a heat-model effect for these pots (P 0.85).

**Lane.** The shipped constants are unchanged (k_asn_glc 117.5 ± 25.2, k_int1_acr 159.2 ± 29.5; the lane already
contains both clusters, one per step). **No formation-Ea band is drawn in the envelope** (`uncertainty.py` draws
B3's fitted Ea and k_acr_dp, not the measured formation Ea). Recommendation (P 0.65 it is worth doing): below
140 °C, the lane's fit window floor, declare the formation-Ea ambiguity as a band 95-170 kJ/mol anchored at 160 °C,
or refuse below 120 °C. Expected effect: central predictions unchanged; the envelope on formation rate widens by
×2.7 at 140 °C, ×7.9 at 120 °C, ×26 at 100 °C (upwards: a lower Ea means more acrylamide below T_ref), and
narrows nothing above 160 °C (×0.41 at 180, ×0.18 at 200 on the low-Ea side). This matters for plant-based
surfaces held at 100-140 °C; it does not move any current score row (all at 140-200 °C).

## CF-104 — Bassama's ×8 k_F at low a_w vs the lane's flat formation multiplier

**Basis check (integrity).** De Vleeschouwer 2008 Fig. 2 (basic system, a_w 0.96, 29.6 % water): Asn0 ≈ 2650 mM,
Glc0 ≈ 2300 mM; per wet kg 2.25, per dry kg 3.20, per L of water 7.6 mol. Per-water is excluded; wet (whole-sample)
vs dry basis is not decided (Glc0 fits wet, Asn0 sits between; P(wet / whole-sample) 0.6). The lane works in mmol/L
of the pot, i.e. on the whole-sample basis, with mass action; T4 does not need DV's basis.

**T4 — concentration basis.** Water 70 / 41 / 8 % wb (a_w 0.972 / 0.904 / 0.430). Sugars ~600 × Asn, so Bassama's
first-order k_F = k2·[reducing sugar]; [sugar] per wet kg ∝ (1 − w) = 1 : 1.97 : 3.07; per dry kg 1 : 1 : 1; per kg
water 1 : 3.36 : 26.8. Observed k_F 1 : 5.10 : 8.06.
| comparison | observed | wet (lane) | residual | dry | per water |
|---|---|---|---|---|---|
| 0.904 / 0.972 | 5.10 (CI 3.7-7.0) | 1.97 | **2.60** | 5.10 | 1.52 |
| 0.430 / 0.904 | 1.58 (CI 1.05-2.38) | 1.56 | **1.01** | 1.58 | 0.20 |
On the lane's own basis the two drier arms agree with a flat multiplier (residual 1.01; likelihood ratio ~11 for
wet-basis mass action over a per-dry flat rate). The wet arm is 2.6 × slower (CI 1.9-3.6) than any basis predicts.

**T5 (outside the conflict) — moisture, not a_w, bounds the calibration.** De Vleeschouwer 2008 lumped formation
k_INT·k_F/(k_F + k_M): basic system 4.4 / 4.9 / 8.7 / 5.4 × 10⁻³ at 10.4 / 14.5 / 29.6 / 34.2 % water; potato-based
10.3 / 10.2 / 13.7 / 14.4 × 10⁻³ at 18.5-25.0 % water: flat within ×2 and ×1.4. The lane's a_w window 0.34-0.99
therefore means ≤ 34 % water in those powders; plantain at a_w 0.972 holds 70 %. Bassama's agreeing arms (8 and
41 %) are at or near that range; the disagreeing arm is twice outside it. H6 holds (rule met).

**T6 — elimination.** k_E 1 : 6.07 : 14.4. Lane first-order (a_w shape 1.0 / 1.0 / 0.70): residual ×6.1 / ×20.6.
Lane bimolecular reading (matrix amines ∝ (1 − w), × the a_w shape): residual ×3.1 / ×6.7; 0.430 / 0.904 ratio 2.38
(CI 1.55-3.66) vs 1.09 predicted. **The elimination conflict is real on every reading** (rule: > ×3 on both).

**Verdict CF-104: open, narrowed.** Formation: the "×8" is mostly the lane's own mass action (wet basis) plus a
high-moisture arm outside the calibrated moisture range (P 0.6 that this is the whole formation story; the
remaining risk is DV's unresolved basis and a_w-at-25 °C). Elimination: a real, unexplained conflict in one food
matrix (P 0.65 real rather than a single-response artefact).

**Lane.** (1) Keep `AW_MULTIPLIER = 1.0` (P 0.75). (2) State the window in moisture as well as a_w: the 0.34-0.99
window is calibrated on ≤ 34 % water; above ~40 % wb either flag out-of-window or extend the band's lower edge from
0.41 to 0.28 (the 95 % lower bound of 1/2.6). Expected effect: no change to any central prediction; for high-moisture
plant-based matrices (meat analogues 55-70 % water) the formation envelope's low edge drops ×0.68. (3) Declare that
the multiplier presumes pot concentrations per whole sample; a caller passing per-dry-matter values loses the (1 − w)
effect (×3 between 70 and 8 % water). (4) Elimination: no change now; record that the only food-matrix series rises
dry-side (×3-7 on the most favourable reading) while the shape falls; if it transferred, dry-side acrylamide would be
over-predicted ×3-7. Decider unchanged (Asn + AA time courses in one food at ≥ 3 moistures, a_w at temperature).

## What a red team should attack
1. DV's concentration basis (wet vs dry; Asn0 2650 mM sits between 2.25 and 3.20). If dry, T5 still holds but the
   statement "DV and Bassama agree on the same basis" weakens to "Bassama's two dry arms agree with the lane's mass
   action".
2. T2a's reading of the 90 °C bar (2-4 px at 600 dpi); the 90-120 Ea is 87-131 across that error.
3. T2b uses dossier-digitised Bassama points (not re-read here); a ±0.05 mg/kg error moves the 140-180 slopes ±10.
4. Whether Knol's k4 and Leuven's k_INT are the "same step" as Gökmen's yield rate (H4's classification).
5. a_w at 25 °C vs at 140-200 °C in Bassama's sealed cell; moisture loss to the cell headspace.

## Side findings (for the single writer)
- `gokmen2012_fig2b_fig3_acrylamide_temperature_ea` (acrylamide.yml) and `gokmen2012_extraction.md` sec. 3.3: the
  "k" behind 83.94 is (AA / C0) / t_hold with t_hold 60 min (90-150 °C) and 10 min (180 °C) — record the estimator
  and flag `not_rate_constants`, `mixed_hold_times`, `precursor_depleted_high_t`.
- `palazoglu2008` gives a second Hacettepe formation Ea, 93.8 kJ/mol (Gökmen & Şenyuva 2006, cited), now in CF-039's
  sources; the primary (10.1080/02652030500482355) is not on disk and would settle what "k′1" is.

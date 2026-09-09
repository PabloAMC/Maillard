# Wave B15 pre-registration — pH and the dry-side water activity on the acrylamide lane

*Written 2026-09-07 after the owner downloaded De Vleeschouwer 2006 (10.1021/jf0611264) and 2007
(10.1021/bp060389f) and before the directional panel was re-scored. Module:
`src/kinetic_core/acrylamide_conditions.py`; licence: FIT_HOLDOUT_DECLARATION.md Amendment 25;
dossiers `devleeschouwer2006_extraction.md`, `devleeschouwer2007_extraction.md`.*

## 1. What the two papers hold for the lane

- **2006 (pH):** 0.1 M Asn + Glc in 0.05 M phosphate at initial pH 4 / 6 / 8, 120-200 C, the same
  formation/elimination fit as the lane's sources. ln k is linear in pH: 0.5414 +/- 0.106 (formation)
  and 0.3442 +/- 0.013 (elimination) per pH unit in natural-log units, i.e. **0.235 and 0.149 decades
  per pH unit**; the potato matrix gives 0.187 / 0.148. Both barriers rise as pH falls (recorded, not
  modelled). Constants are indexed to INITIAL pH; the pot drifts by up to 2 units while heating.
- **2007 (a_w 0.34-0.92):** formation k_F "varies only slightly" (0.71-1.09 of the 0.92 value);
  elimination k_E has a minimum at a_w 0.82 (0.33 of the 0.92 value; SEs 20-90 %); the Maillard
  competition peaks there. Together with 2008 (0.88-0.99, flat) the lane's a_w axis is measured over
  0.34-0.99.

## 2. What changes (declared, no fit)

- The B14 formation window widens from 0.88-0.99 to **0.34-0.99** (still flat, band 0.41-1.39).
- An ELIMINATION a_w multiplier on `k_acr_dp`: 0.76 / 0.60 / 0.33 / 0.37 at a_w 0.34 / 0.59 / 0.73 /
  0.82, joining 1.0 at 0.88 (the 2008 flat window; 2007's own 0.88 point, 0.66 +/- 0.33, is within one SE
  of 1); the envelope scales the deficit by 0-1.2.
- A pH factor 10^(0.235 (pH - 6.8)) on `k_asn_glc` (band 0.114-0.281) and 10^(0.149 (pH - 6.8)) on
  `k_acr_dp` (band 0.116-0.155), measured window pH 4-8, HELD at the window edge outside it with a
  warning; exactly 1 at the lane's reference pH 6.8.
- Engine: the acrylamide lane leaves `NO_PH_TERM_LANES`; pH comparisons inside 4-8 are answered,
  comparisons that leave the window are refused with the window named; the declaration replaces the
  old "NO pH term" line. Every panel row at a_w None and pH 6.8 reproduces exactly.

## 3. Pre-registered expectations on the directional panel

- PH-ACR-01 (new; `fit_adjacent`, declared from the same paper): evaluable and AGREE (formation slope
  exceeds elimination slope, so net acrylamide at 20 min rises with pH). Not in the independent headline.
- AW-02 (extrusion, a_w 0.3 / 0.6 / 0.9): the 0.3 arm is BELOW the window floor 0.34, so it stays NOT
  EVALUABLE. Recorded: had it been 0.34, the lane would now predict a monotone RISE with a_w (formation
  flat, elimination slowest at 0.82) against the claim's peak at 0.6.
- TEMP-01 / TEMP-02 (a_w 0.5, pH 6.0): both arms now carry the same elimination factor (a_w 0.5 ->
  x0.65) and the same pH factor, so directions are unchanged; levels move.
- Headline 18/30 unchanged; all claims 21/43 -> 22/44 (PH-ACR-01), plus DIC-01 from the dicarbonyl
  paper if the ranking scorer accepts it.

## 4. What would falsify the declared terms

A pH series in an aqueous or low-moisture Asn-Glc system whose net-acrylamide ordering between pH 4
and 8 at a short hold reverses; or an a_w series between 0.34 and 0.92 whose acrylamide at fixed
temperature and time is not highest near a_w 0.82.

## 5. OUTCOME (2026-09-07, after the re-score)

- PH-ACR-01: AGREE (fit_adjacent, wiring check). AW-05: still AGREE (the elimination shape joins the
  2008 flat window at 0.88). AW-02: still NOT EVALUABLE (0.30 below the 0.34 floor).
- Alongside, the two scorer additions of the same day: the UNSTATED-INPUT SWEEP (WANG-01/02 charged
  as TTCA 197 mM in 0.2 M phosphate over the paired ladder, scored at pH 5 / 7 / 9): **WANG-02 (FFT
  peak) AGREES at every pH** and enters the independent headline; **WANG-01 (MFT peak) is NOT
  EVALUABLE** because its verdict flips with the unstated pH (agree at 7, disagree at 5 and 9). And the
  RANKING scorer: DIC-01 (Zhang 2020 dicarbonyl ordering) DISAGREES (see `kinetic_core_b13_prereg.md`).
- Headline 18/30 -> **19/32**; all claims 21/43 -> 23/46; panel 75 -> 78 claims.

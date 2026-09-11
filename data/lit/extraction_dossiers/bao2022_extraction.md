# Bao et al. 2022 — EXTRACTION (hexanal from neat linoleic acid at seven temperatures — uncalibrated peak areas, REFUSED as a rate source; kept as an onset ordering)

**Source on disk:** `data/articles/bao2022.pdf` (downloaded 2026-09-11 at this repository's request).
Read 2026-09-11 via `pdftotext -layout`, with both halves of the split tables verified against
rendered page images. Wave B45.

| field | value |
|---|---|
| Title | "Detailed temperature-dependent study of linoleic acid oxidative decomposition into volatile compounds in the heating process" |
| Authors | Y. Bao, J. Du, C. Xu, M. Wang, B. Wang, L. Xiao, K. Cheng, L. Dong (Dalian Polytechnic University; CFAA Beijing) |
| Venue | Journal of Food Processing and Preservation **46** (2022) e16445 |
| DOI | 10.1111/jfpp.16445 |
| System | **20 µl of neat linoleic acid** (analytical standard) in a 20 ml glass vessel inside a Markes M-CTE250 micro-chamber. **No emulsion, no water, no buffer, no protein, no matrix.** |
| Temperatures | **30, 60, 90, 120, 150, 180, 210 °C, 30 min each**; plus time series at 60 °C (printed) and 90 °C (figure only) at 2, 4, 6, 8, 10 h |
| Units | **"Peak area (×10⁶)"** — see §2 |
| Data availability | "Research data are not shared." |

## 1. Why this paper was fetched

The reading list wanted a **linoleic acid temperature series** that could adjudicate between the
bulk-oil hexanal barrier (114–122 kJ/mol) and the moist-system one (61–65 kJ/mol). Seven
temperatures over 180 K, hexanal at every one of them, triplicates with standard deviations: on the
face of it the densest temperature series on disk.

## 2. Why the numbers cannot be a rate — the units first

Both tables carry the printed column header **"Peak area (×10⁶)"**. There is **no internal standard,
no calibration curve and no recovery correction** anywhere in §2.1–2.4. Authentic standards are used
for **identity only**, verbatim: "identification accuracy was determined by separating relevant
standard compounds through GC–MS analysis under the same conditions"; the table footnote defines STD
as "comparison with a standard compound".

So every number below is an **uncalibrated GC–MS peak area of a purged headspace trapped on Tenax
TA.** Under this repository's rule a peak area is not a concentration and is not scored as one.

**And the within-study-ratio escape does not open here.** A ratio of two peak areas of the *same*
compound on the *same* instrument would ordinarily survive that rule. It fails in this design for a
reason specific to it: after each hold, "the air valve was opened with a flow rate of 100 ml/min for
30 min" to strip the chamber onto the trap. A peak area therefore conflates formation with
**evaporative partition, purge efficiency and trap breakthrough — all of which are themselves
temperature-dependent** — with no mass balance and no measurement of the linoleic acid remaining.
A hotter cell delivers more of whatever it holds to the trap. The ratio is not a ratio of amounts.

## 3. Table 1 — hexanal and the lipid-lane companions, 30 min at each temperature

Peak area ×10⁶, mean ± SD (n = 3); superscripts are the printed one-way ANOVA groupings (p < .05);
"nd" = not detected.

| compound | 30 °C | 60 °C | 90 °C | 120 °C | 150 °C | 180 °C | 210 °C |
|---|---|---|---|---|---|---|---|
| **Hexanal** | **83.76 ± 2.96 c** | **54.26 ± 7.42 d** | **364.26 ± 28.75 b** | **581.74 ± 82.97 a** | **520.7 ± 65.4 a** | **568.38 ± 30.25 a** | **594.42 ± 43.7 a** |
| Pentanal | 19.4 ± 2.71 b | 20.44 ± 4.04 b | 52.19 ± 10.27 a | 52.19 ± 10.27 b | 46.01 ± 8.20 a | 40.61 ± 6.20 a | 138.47 ± 2.12 a |
| Heptanal | nd | nd | 5.11 ± 0.71 d | 17.48 ± 3.14 c | 33.29 ± 6.25 b | 86.45 ± 6.28 a | 93.1 ± 6.56 a |
| Octanal | nd | nd | 0.81 ± 0.01 e | 1.55 ± 0.17 d | 2.25 ± 0.22 c | 6.8 ± 0.67 b | 12.96 ± 0.04 a |
| Nonanal | nd | nd | nd | 4.88 ± 1.65 b | 3.35 ± 0.17 c | 4.62 ± 0.18 b | 6.5 ± 0.25 a |
| (E,E)-2,4-Decadienal | 0.78 ± 0.11 d | 22.93 ± 5.21 d | 112.78 ± 28.31 c | 84.98 ± 8.90 c | 350.92 ± 31.68 a | 218.25 ± 18.68 b | 199.91 ± 36.79 b |
| 1-Hexanol | 0.58 ± 0.10 e | 1.73 ± 1.09 d | 3.31 ± 0.31 c | 3.99 ± 1.32 c | 5.61 ± 0.08 b | 12.33 ± 1.09 a | 13.69 ± 1.78 a |
| 2-Pentylfuran | 15.47 ± 2.88 c | 10.59 ± 1.42 c | 10.59 ± 1.42 c | 78.8 ± 2.25 b | 798.95 ± 143.74 a | 768.65 ± 83.46 a | 698.74 ± 6.90 a |
| (Z)-2-Heptenal | 66.54 ± 7.24 c | 144.19 ± 18.94 c | 312.34 ± 23.03 b | 656.12 ± 79.05 a | 697.29 ± 111.23 a | 686.39 ± 20.78 a | 659.03 ± 46.74 a |
| (E)-2-Octenal | 23.02 ± 2.49 d | 83.08 ± 16.45 d | 315.01 ± 7.98 c | 508.44 ± 74.00 b | 426.76 ± 77.66 b | 191.38 ± 11.32 d | 630.69 ± 51.52 a |
| 1-Octen-3-ol | 12.11 ± 2.39 c | 58.74 ± 7.21 b | 116.05 ± 9.44 a | 152 ± 36.80 a | 118.1 ± 34.08 a | 120.24 ± 13.67 a | 120.04 ± 14.49 a |
| Hexanoic acid | 4.17 ± 0.93 f | 54.87 ± 7.51 e | 128.21 ± 7.05 d | 151.98 ± 18.31 d | 310.08 ± 146.24 c | 458.64 ± 60.97 b | 571.38 ± 92.95 a |
| Furfural | 0.74 ± 0.16 e | 1.2 ± 0.28 e | 3.41 ± 0.12 d | 10.75 ± 2.33 c | 23.59 ± 2.88 b | 46.86 ± 3.22 a | 54.96 ± 3.88 a |

(The full 42-compound table is in the paper; the rows above are the ones that map onto this model's
lipid lane or its furan/acid exits. The remaining 29 rows carry no species this model tracks.)

**Printed anomalies, flagged and not corrected:** pentanal at 90 and 120 °C is the identical string
"52.19 ± 10.27" with different ANOVA letters; 2-pentylfuran at 60 and 90 °C is identically
"10.59 ± 1.42 c"; 2-butyl-furan at 60 °C reads "1.27 ± 0.21 c" and at 90 °C "1.27 ± 0.21 e".
1-Butanol's retention index is 596 in Table 1 and 569 in Table 2. These read as typesetting
duplications in the source.

## 4. Table 2 — hexanal against time at 60 °C (peak area ×10⁶)

| 2 h | 4 h | 6 h | 8 h | 10 h |
|---:|---:|---:|---:|---:|
| 484.18 ± 52.17 c | 470.85 ± 36.42 c | 463.95 ± 25.46 c | 547.12 ± 28.66 b | 640.79 ± 40.29 a |

The 90 °C time series exists only as a colour heat-map (Fig. 6b) — verbatim, "data obtained during
heating at 90°C are not shown". **Figure only.**

## 5. Why no barrier follows — five reasons beyond the units

1. **One time point per temperature** in the temperature series (30 min). A single endpoint gives a
   rate only under an assumption of linearity, and the data refute linearity: hexanal above 120 °C is
   flat (581.74 → 520.7 → 568.38 → 594.42, all ANOVA group "a", i.e. **not significantly different**).
   A plateau is not a kinetically limited regime.
2. **The low end runs the wrong way.** Hexanal at 30 °C is **83.76 ± 2.96** and at 60 °C is
   **54.26 ± 7.42** — a *decrease* with rising temperature, and the paper's own letters (c vs d) call
   it significant. No Arrhenius form accommodates that, and the paper offers no explanation.
3. **Time series at only two temperatures, one of them unavailable.** 60 °C is printed; 90 °C is
   figure only. Two temperatures is the bare minimum for a barrier; one is not enough.
4. **The two 60 °C runs disagree by ~9×** — 54.26 at 30 min (Table 1) against 484.18 at 2 h
   (Table 2), from different runs of the same nominal system with no stated normalisation between
   them.
5. **No substrate depletion, no extent of conversion, no headspace gas volume, and no stated
   atmosphere.** The only gas statement is the 100 ml/min purge; oxygen concentration is never given
   and the chamber is not said to be purged before heating. Nothing to normalise a rate against, and
   for an oxidation that is disqualifying on its own.

**No hydroperoxide measurement of any kind** — no peroxide value, no conjugated dienes, no acid
value, no p-anisidine, no TBARS, no residual linoleic acid. Hydroperoxides appear only as
mechanistic discussion citing Frankel 1980/1982, with no number attached.

## 6. Kinetics

**None.** No rate constant, no activation energy, no reaction order, no half-life, no Arrhenius plot,
no kinetic model; the words "activation energy", "rate constant", "Arrhenius" and "kinetics" do not
appear. The statistics are one-way ANOVA (p < .05), PCA (PC1 29.5 %, PC2 29.1 %) and PLS-DA with
VIP > 1. Mechanism is deferred: "In the following research, we conducted an in-depth study on the
mechanism of thermal oxidation of linoleic acid in combination with quantum chemical calculation
methods." **That future work is a computational study and is out of scope here under the no-computed-
numbers rule, whatever it reports.**

## 7. What the repository can bank — an onset ordering, not a rate

The paper's "forming temperature" (Fig. 2, a donut chart) is **never formally defined**. Verified
against Table 1, it is operationally the **lowest of the seven sampled temperatures at which the
compound is detected at all** — a detection-limit-bounded first appearance on a coarse 30 K grid, not
a thermodynamic or kinetic onset. Recorded on those terms:

- **Hexanal's forming temperature is ≤ 30 °C** — it is already present at the lowest temperature
  sampled, so the true onset is not resolved.
- Heptanal, octanal, 2-heptanone and (E)-2-octenoic acid first appear at **90 °C**; nonanal,
  2-octanone, 6-undecanone and nonanoic acid at **120 °C**; propanoic acid, 2-undecenal,
  5-methyl-2-furancarboxaldehyde and dihydro-5-propyl-2(3H)-furanone at **210 °C**. No new compound
  appears at 180 °C.
- General claim, verbatim: "most volatiles with shorter carbon chains were generated at lower
  temperatures, while volatiles with longer carbon chains were generated at higher temperatures."

For this model the useful structural statement is that **hexanal from linoleic acid saturates**: it
rises steeply from 60 to 120 °C and then does not move significantly across the next 90 K. Whatever
limits it above 120 °C is not the forward step's barrier. That is consistent with the lipid lane's
existing hydroperoxide-pool bottleneck and is recorded as a direction.

## 8. Verdict

**REFUSED as a rate or barrier source** on the units alone, and independently on the single time
point, the inverted 30→60 °C step, and the undefined atmosphere. It cannot adjudicate the
114–122 vs 61–65 kJ/mol dispute. Kept for the onset ordering in §7 and the saturation direction.

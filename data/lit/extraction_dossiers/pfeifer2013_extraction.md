# Pfeifer, Haase & Kroh 2013 (J. Agric. Food Chem. 61:3090-3096) — EXTRACTION (file stem `pfeifer2013`)
**DOI:** 10.1021/jf302959k (printed on every page footer; Crossref 2026-11-05: same title, authors Yvonne V. Pfeifer, Paul T. Haase, Lothar W. Kroh, JAFC 61(12) 3090-3096, online 2013-03-18) — matches the brief.
### *"Reactivity of Thermally Treated α-Dicarbonyl Compounds"*
### 3-Deoxyglucosone (3-DG, "3-DH" in the paper) or methylglyoxal (MGO), "0.1 mmol" in 10 mL, alone or with L-alanine (1 mmol), water or phosphate buffer, **initial pH 5 or 8** (3 N HCl / NaOH), sealed ampoules in a thermoblock, 100 or 130 °C, up to 300 min; OPD trapping, quinoxalines by HPLC-DAD. **Only 130 °C is shown, MGO and 3-DG ALONE only (Figs. 1-4); no MGO + Ala data, no 100 °C data, no replicates.** MGO alone, pH 5, 130 °C, aerated: 82 / 78 / 48 / 46 / 20 % left at 30 / 60 / 120 / 180 / 300 min (k1 ≈ 0.0053 min⁻¹ over 0-300 min, noisy); pH 8: 13 % left at 30 min. **The same Berlin laboratory lineage as Kanzler / Bork (sec. 0): not a second laboratory.**

**Read 2026-11-05 (first reader) from `data/articles/pfeifer2013.pdf`** (7 pp., ACS publisher PDF, clean text layer).
Prose from `pdftotext -layout`. **Figs. 1, 3 and 4 digitised** from a 300 dpi render of p. 4 (Figs. 3, 4) and a 110 dpi
render of p. 3 (Fig. 1): marker centres located with a dark-pixel run finder (scratchpad `m4tools/runs.py`, crops
`m4img/pf_f3.png`, `m4img/pf_f4.png`) against the frame lines (35 mg/L = top frame, 0 = axis). **Read error ±0.3 mg/L
(±1 % of C0) for Figs. 3-4; ±1 mg/L for Fig. 1.** Arithmetic: scratchpad `m4tools/m4_arith.py` (maillard env, `-I`). No SI.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Yvonne V. Pfeifer, **Paul T. Haase** (corresponding, paul.haase@tu-berlin.de), Lothar W. Kroh |
| Affiliation | Institute of Food Technology and Food Chemistry, Technical University of Berlin, **Gustav-Meyer-Allee 25, 13355 Berlin** |
| Venue | *J. Agric. Food Chem.* **61** (2013) 3090-3096 |
| DOI | 10.1021/jf302959k |
| Dates | received 11 Jul 2012; revised 19 Feb 2013; accepted 23 Feb 2013; published 23 Feb 2013 (ASAP) |

**Laboratory: the same group lineage as Kanzler 2026 / 2021 and Bork 2022, by the record.**
- Same institute and the same street address as `bork2022` (Bork, **Haase**, Rohn & Kanzler 2022: TU Berlin, Dept. of
  Food Chemistry and Analysis, Gustav-Meyer-Allee 25). Haase is the corresponding author here and a co-author there.
- Kanzler 2021 (`kanzler2021`) takes its OPD / HPLC-DAD MGO method from its ref. 29, **Kanzler, Schestkowa, Haase &
  Kroh 2017** (JAFC 65:8957), and cites Kanzler, Haase, Schestkowa & Kroh 2016 and Kanzler & Haase 2020 [C]: Kanzler
  published with both Haase and Kroh in this department before the Rohn-chair papers.
- Same protocol family: dicarbonyl in water, pH adjusted once with HCl / NaOH, sealed ampoules in a Behr heating block,
  OPD trapping, quinoxalines by HPLC-DAD.

**Verdict (first reader): one laboratory with `tu_berlin_kanzler` (proposed in `bork2022`) for leave-one-laboratory-out.**
The chair changed (Kroh → Rohn) but people, address, assay and protocol carry over; `data/keys/laboratories.yml` defines
a laboratory as "shared methods, instruments and calibrations". Counting Pfeifer as a second laboratory would not be
conservative. Owner call; the evidence above is from the two PDFs' first pages and reference lists.

## 1. What this paper answers for the model

| model need | what Pfeifer 2013 gives | tag |
|---|---|---|
| MGO loss alone, pH 5, 130 °C | Fig. 3: 31.1 / 25.6 / 24.2 / 14.95 / 14.2 / 6.3 mg/L at 0 / 30 / 60 / 120 / 180 / 300 min; text: "20% of the starting value" at 300 min | [F] / [M] |
| MGO loss alone, pH 8, 130 °C | Fig. 4: 31.0 / 4.1 / 1.15 / 1.15 / 0.78 / ~0.45 mg/L | [F] |
| MGO alone at a second temperature | **none shown** (100 °C is in the methods only) | — |
| MGO + amino acid | **no numbers**: "degradation of proceeded faster" with Ala [V]; A420 0.8 vs 1.3 (×100 dilution) after 300 min at pH 5 / 130 °C without / with Ala [M] | [Q] / [M] |
| 3-DG loss alone, pH 5, 130 °C | Fig. 1: 58 / 33 / 12 / 2.5 / ~0.5 / 0 mg/L at 0 / 30 / 60 / 120 / 180 / 300 min; text "nearly 50%" decomposed at 30 min, not detectable after 180 min | [F] / [M] |
| products | MGO pH 5: 3-deoxytetrosulose (≤ 0.7 mg/L, Fig. 3), diacetyl (≤ 0.02 mg/L), formaldehyde and acetaldehyde (not quantified), an MGO dimer "almost up to 50%" of starting MGO (LC/MS, HPTLC; data not shown); pH 8: acetic acid (6.4 / 6.0 mg/L at 30 / 60 min, 0 after; absent under deaeration [Q]), diacetyl ≤ 0.87 mg/L; 3-DG pH 5: HMF (max 28.5 mg/L at 30 min), 3,4-dideoxyglucoson-3-ene (≤ 4 mg/L), traces of MGO and diacetyl | [F] / [Q] |
| glycine, cysteine, NF, MFT | **none** | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Dicarbonyl pots | MGO or 3-DG, **0.1 mmol** in **10 mL** water or phosphate buffer (pH 5) [Z: 10 mM nominal], ± L-Ala 1 mmol (100 mM, 10× excess) |
| Glucose pots | glucose 1 mmol ± Ala 1 mmol in 10 mL (100 mM each); no glucose data shown |
| pH | adjusted at the start to 5 (3 N HCl) or 8 (3 N NaOH); **label pH; not re-measured**; at pH 8 the authors state the pH "decreased" toward the pH 5 models as acids formed [Q] |
| Buffer | "H2O or phosphate buffer (pH 5)" — [!] **which one Figs. 1-4 used is not stated**; phosphate molarity not stated |
| Vessel / heating | sealed ampoules (fill and ampoule volume not stated) in a thermoblock (Behr behrotest ET2), 100 ± 1 or 130 ± 1 °C, up to 300 min |
| Atmosphere | Fig. 3 "under aerated conditions" [V]; a deaerated arm exists (water degassed with He, samples Ar-flushed) but no deaerated MGO curve is shown, only the statement that acetic acid is absent there |
| Sampling | 500 µL sample + 500 µL OPD 0.05 M, 3 h at 25 °C, then HPLC-DAD (320 nm) and GC/MS (acetylated) |
| Replicates | **not stated**; no error bars on any figure |

## 3. Assays [M]

- **Dicarbonyls:** OPD quinoxalines, HPLC-DAD at 320 nm, quantified against "standard solutions containing amounts of
  pure, authentic reference compounds" (2-methylquinoxaline for MGO). **No LOD / LOQ, no recovery, no calibration
  range.**
- Acetylated quinoxalines by GC/MS (SIM, diphenylquinoxaline internal standard); acetic acid as decyl ester by GC-FID
  (chlorosuccinic acid internal standard): "standard deviations <20 mmol/mol methylglyoxal", CV < 3.5 %, R² > 0.97 [M]
  — the only precision statement in the paper, for acetic acid only.
- Colour: A420 (×100 dilution for the MGO pots).
- **Replicate noise for MGO: not reported.** The Fig. 3 series is non-monotone in rate (pointwise k1 0.0065, 0.0019,
  0.0080, 0.0009, 0.0068 min⁻¹ [Z]), which suggests ≥ ±5 % of C0 noise per point.

## 4. Data

### 4.1 Fig. 3 — MGO alone, pH 5 (label), 130 °C, aerated [F ±0.3 mg/L]

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| MGO (mg/L) | 31.1 | 25.6 | 24.2 | 14.95 | 14.2 | 6.3 |
| fraction left [Z] | 1.000 | 0.823 | 0.778 | 0.481 | 0.457 | 0.203 |
| 3-deoxytetrosulose (mg/L, right axis) | 0 | 0.70 | 0.50 | 0.22 | 0.12 | 0 |
| diacetyl (mg/L, right axis) | 0 | ≤ 0.01 | ≤ 0.01 | ≤ 0.01 | ≤ 0.01 | ≤ 0.01 |

[!] **The time axis is categorical** (0, 30, 60, 120, 180, 300 equally spaced): slopes read off the plot are wrong;
use the tick labels.

[!] **Concentration scale.** 0.1 mmol MGO in 10 mL is 10 mM = 721 mg/L; the plotted start is 31.1 mg/L (0.43 mM,
4.3 % of nominal). Fig. 1's 3-DG starts at 58 mg/L vs 1621 mg/L nominal (3.6 %). The 1 + 1 OPD dilution explains a
factor 2 only. Either the plotted values are of a further-diluted sample, or the charge differs from the methods.
**Fractions within a series are unaffected; absolute levels and any second-order constant are not usable.**

### 4.2 Fig. 4 — MGO alone, pH 8 (label), 130 °C [F ±0.3 mg/L]

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| MGO (mg/L) | 31.0 | 4.1 | 1.15 | 1.15 | 0.78 | ~0.45 |
| acetic acid (mg/L) | 0 | 6.4 | 6.0 | 0 | 0 | 0 |
| diacetyl (mg/L, right axis) | 0 | 0.33 | 0.82 | 0.87 | 0.65 | 0.38 |

[!] Acetic acid falling to zero after 60 min in a sealed pot is unexplained (volatility in sampling, or the GC-FID
calibration range).

### 4.3 Fig. 1 — 3-DG alone, pH 5, 130 °C [F ±1 mg/L]

3-DG 58 / 33 / 12 / 2.5 / ~0.5 / 0; HMF 0 / 28.5 / 21 / 16 / 13 / 9; 3,4-DDH 0 / 4.0 / 4.0 / 1.0 / ~0.5 / ~0 mg/L at
0 / 30 / 60 / 120 / 180 / 300 min. 3-DG pointwise k1 0.019 → 0.034 → 0.026 min⁻¹ [Z] (accelerating, not first order).
Fig. 2 (3-DG, pH 8): 3-DG ≈ 0 throughout on the plotted scale; 2-deoxyribonic acid 220 → 20 mg/L (present at t = 0 from
oxidative cleavage at room temperature [Q]); pyruvate to 3.2 mg/L at 180 min. Not used further.

### 4.4 Arithmetic [Z] (`m4_arith.py`)

1. **MGO alone, pH 5, 130 °C (Fig. 3).** First order, cumulative k1 = 0.0065 / 0.0042 / 0.0061 / 0.0044 / 0.0053 min⁻¹
   (0-30 … 0-300): roughly constant, **k1 ≈ 0.0053 min⁻¹ (range 0.0042-0.0065)**. Second order at the nominal 10 mM:
   cumulative k2 rises 7.2 × 10⁻⁴ → 1.31 × 10⁻³ mM⁻¹ min⁻¹ (×1.8 drift); at the plotted 0.43 mM, 0.011-0.030. **First
   order describes this series better than second order**; with no replicates, neither is identified.
2. **pH 8:** k1(0-30) = 0.067 min⁻¹, ×13 faster than pH 5 [Z]; afterwards a tail at 1-4 % (≈ 0.01 dex per pH unit is
   meaningless here: label pH only, and it drifts toward 5).

## 5. What it can and cannot test in the model

**Can:**
- **The engine's MGO sink at 130 °C.** `k_mgo_mel` (B1 fit: log10 k(100 °C) = −1.643, Ea 20.04 kJ/mol at its lower
  bound; `results/validation/core_prediction_uncertainty.md` rows `b1.k_mgo_mel.*`) gives 0.0368 min⁻¹ at 130 °C [Z].
  Pfeifer's MGO alone at pH 5 / 130 °C: 0.0053 min⁻¹ — **the engine is ×7 faster than the measured amine-free loss
  at 130 °C** (the engine term is amine-free too).
- A pH contrast for MGO alone (pH 5 vs 8 at 130 °C, ×13 at 30 min), label pH only.
- A 3-DG lifetime at 130 °C / pH 5 (t½ ≈ 30 min, accelerating) and its HMF yield (≈ 60 mol % at the 30 min peak on the
  plotted scale [Z]: 28.5 / 126 vs 58 / 162 mg/mmol).

**Cannot:**
- **A temperature law for MGO alone**: one temperature only. The 100 °C arm in the methods has no data in the paper.
- MGO + amine kinetics (no data), glycine, cysteine, NF, MFT.
- Absolute concentrations (sec. 4.1 scale flag); replicate noise (none reported).

## 6. Conflicts and agreements with the repo's records

1. **[!] Against Kanzler 2026 / Bork 2022 (same lineage), on order and size.** Kanzler 2026: MGO alone 100 mM, 100 °C,
   second order, k2 ≈ 1.85 × 10⁻⁴ mM⁻¹ min⁻¹. `bork2022`: MGO alone 23 mM, 130 °C, pointwise k2 2.84 → 1.84 × 10⁻³,
   pointwise k1 0.036 → 0.0087 min⁻¹ (decelerating). Pfeifer, 130 °C, 10 mM nominal: k1 ≈ 0.0053 min⁻¹, not
   decelerating; Bork's second-order law at 10 mM would give an initial pseudo-first-order 0.018-0.028 min⁻¹, ×3-5 above
   Pfeifer's. Cross-paper Ea (second-order reading, Pfeifer at the nominal 10 mM): from Kanzler 2026's 100 °C value
   39-82 kJ/mol; from Kanzler 2021's ([Z] in `kanzler2021`) −9 to 33 kJ/mol. **These spans say only that the three
   papers' MGO-alone numbers do not share one rate law within this lineage**; none of them fixes the MGO-alone Ea.
2. **[!] Against the engine's `k_mgo_mel` (CF-122 on branch `registers-2026-11-05`):** at 130 °C the fitted term is
   ×7 above Pfeifer's MGO-alone loss. This does not by itself contradict Ea 20 (no temperature series here); it says the
   fitted level at 130 °C exceeds the amine-free measurement, as it does at 80 °C against Kanzler 2026.
3. Agrees with Kanzler 2021 / 2026 qualitatively: MGO alone at pH 5 forms aldol oligomers (the dimer, m/z 145 / 146.06),
   and fragments to formaldehyde and acetaldehyde [Q].
4. [!] Internal: concentration scale (sec. 4.1); categorical time axis; buffer of the plotted series unstated.

## 7. Provenance tiers (first reader only; second reader pending)

- Fig. 3 / Fig. 4 MGO fractions: `direct_measurement`, figure-digitised, **confidence low** (single series, no
  replicates, no LOD, buffer unstated, absolute scale inconsistent with the methods).
- Every k and every Ea in secs. 4.4 and 6 is [Z], tier `derived`; not for write-back.
- No ledger entry proposed: one temperature, no replicates; it can enter a band only as a sensitivity.

## M4 role (case 8660c005)

**Role fixed before reading:** band source for MGO's self-loss law, "a second laboratory beside Kanzler (lab boundary to
be confirmed)"; measurement that would count: "MGO alone at >= 2 temperatures"; P ~0.4.

- **What it supplies.** Coordinate: the MGO loss law (`k_mgo_mel` and any MGO sink on the one network). One value: MGO
  alone, pH 5 (label, unbuffered or phosphate, unstated), 130 °C, aerated, 10 mM nominal: **k1 ≈ 0.0053 min⁻¹
  (0.0042-0.0065), first order fits better than second**, 0-300 min. pH 8: k1(0-30) ≈ 0.067 min⁻¹. Transfer: amine-free
  and thiol-free; one temperature; level only. It can check a band's level at 130 °C (the case pins the level at
  100 °C to B60's 0.0228 min⁻¹ and widens Ea to 95-121; at Ea 95-121 that pinned level predicts 0.22-0.42 min⁻¹ at
  130 °C [Z], ×42-78 above Pfeifer's amine-free 0.0053). That gap is the amine-free vs amine-carrying distinction the
  case already names, not a new conflict — but it shows the pinned level is an amine-pot level.
- **The named measurement is NOT met:** MGO alone at one temperature only.
- **Second laboratory: NO.** Same department lineage, address, people (Haase) and assay as Kanzler / Bork (sec. 0). It
  widens no laboratory count.
- **Decisive criterion:** none met, none expected (no thiol, no hexose MFT).

# Hayashi & Namiki 1986 (Agric. Biol. Chem. 50:1965-1970) — EXTRACTION (file stem `hayashi1986`)
**DOI:** 10.1080/00021369.1986.10867692 (printed on the Taylor & Francis cover page; Crossref 2026-11-05: same title, authors Tateki Hayashi, Mitsuo Namiki, Agric. Biol. Chem. 50(8) 1965-1970, issued 1986-08) — matches the brief.
### *"Role of Sugar Fragmentation in an Early Stage Browning of Amino-carbonyl Reaction of Sugar with Amino Acid"*
### D-Glucose + **β-alanine, 2 M each, in water, initial pH 3.51 / 6.42 / 9.29** (adjusted once with NaOH / HCl, unbuffered), **boiling water bath**, 0-120 min: C2 and C3 sugar fragments ("C2 imine", "C3 imine", measured by a modified OPD method with capillary GC), A420 and the Amadori product (HPLC-RI). **C3 fragment: ~0 at pH 3.51, peak 0.013 mM at pH 6.42, peak 0.117 mM at pH 9.29; initial rate ×15 from neutral to alkaline (0.41 dex per pH unit); C2 ×4.6; Amadori ×1.5.** Plus browning rates of glucose, xylose, MGO, glyoxal, glycolaldehyde, glyceraldehyde and Amadori / osone intermediates with β-Ala (Tables I, II).

**Read 2026-11-05 (first reader) from `data/articles/hayashi1986.pdf`** (7 pp.: one T&F cover page plus the six
article pages, a scan with an OCR text layer). Prose and Tables I-IV from `pdftotext -raw` (OCR noise in symbols;
numbers checked against the page image). **Fig. 1A digitised by eye** from a 400 dpi render of p. 3 (scratchpad
`m4img/hayashi1986_400-3.png`, crop `m4img/hay_f1a.png`, ~300 dpi effective): **read error ±0.02 mM (C2 panel),
±0.003 mM (C3 panel), ±3 min** in time (sampling times are not printed). Figs. 1B, 1C and 2 not digitised (their initial rates are printed in Table III). Arithmetic: scratchpad
`m4tools/m4_arith.py`.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Tateki Hayashi, Mitsuo Namiki |
| Affiliation | Department of Food Science & Technology, Faculty of Agriculture, Nagoya University, Nagoya 464 |
| Venue | *Agric. Biol. Chem.* **50**(8) (1986) 1965-1970 |
| DOI | 10.1080/00021369.1986.10867692 |
| Dates | received 23 Dec 1985 |

**Laboratory.** Namiki group, Nagoya — the same laboratory as `hayashi1980` (glyoxal dialkylimine from glucose +
alkylamines in ethanol). The C3 assay is described in Hayashi, Mase & Namiki 1986 (ref. 8, Agric. Biol. Chem. 50:1959),
not on disk.

## 1. What this paper answers for the model

| model need | what Hayashi 1986 gives | tag |
|---|---|---|
| C2 / C3 fragment supply vs pH from glucose + amine | Fig. 1A time courses at three initial pH (sec. 4.1); Table III initial rates (neutral, alkaline) | [F] / [M] |
| acid limb | "negligible in acidic mixtures" [V]; the pH 3.51 series sits on the baseline in both panels | [Q] / [F] |
| Amadori product formation vs pH | Fig. 1C; Table III 9.8 (neutral) vs 14.8 mM/min (alkaline) | [F] / [M] |
| browning ability of fragments vs sugars | Table I (sec. 4.2): MGO 654×, glycolaldehyde 2109×, glyceraldehyde 1967×, glyoxal 122× glucose (relative 1/min to reach O.D. 5.0) | [M] |
| browning of intermediates | Table II: glucosone + β-Ala 23.7×, 3-DG + β-Ala 137.5×, Amadori 0.82× | [M] |
| melanoidin N/C vs pH | Table IV | [M] |
| identity of the fragments | assumed to be glycolaldehyde, glyoxal, MGO, glyceraldehyde "or their imine derivatives" [Q]; measured as "C2 / C3 imine" after OPD | [Q] |
| α-amino acid, thiol, MFT | **none**: the amine is β-alanine | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Fig. 1 pots | glucose + β-alanine **2 M each** in water, "heated in a boiling water bath" (≈ 100 °C; not printed) |
| pH | adjusted "at the starting point" with NaOH or HCl: **3.51 / 6.42 / 9.29**; label pH, unbuffered, not re-measured [!] |
| Fig. 2 / Table I pots | sugar or carbonyl + β-Ala **1 M each**, water, 80 °C (fragments, xylose) and 95 °C (xylose, glucose, fructose); pH not stated |
| Table II pots | intermediates 1 M each (± β-Ala), boiling water bath, 20 min |
| Table IV | glucose + β-Ala 2 M, initial pH 2.3 / 6.5 / 9.2, heated 190 / 53 / 15-25 min to matched browning; glycolaldehyde + β-Ala 1 M, pH 5.9, 30 min |
| Vessel / headspace / replicates | **not stated** (Table IV gives ± for N/C only) |

## 3. Assays [M]

- **C2 / C3 products:** "a modification of the o-phenylenediamine method" with capillary GC (details in ref. 8). No
  calibration, LOD or recovery printed. Units mM.
- Amadori product: HPLC (JEOL TWINCLE, Shodex RI SE-II, LiChrosorb NH2, acetonitrile-water 70:30, 2 mL/min).
- Browning: O.D. at 420 nm. Melanoidin: Sephadex G-15 fraction, elemental N/C.
- **Replicate noise: not reported** (no error bars).

## 4. Data

### 4.1 Fig. 1A — C2 and C3 products (mM), glucose + β-Ala 2 M each, boiling water [F]

| initial pH | product | ~3 min | ~8 min | ~15 min | ~20 min | ~30 min | ~45 min | ~62 min |
|---|---|---|---|---|---|---|---|---|
| 9.29 (squares) | C2 imine | 0.39 | 0.73 (peak) | 0.64 | 0.60 | 0.52 | — | 0.17 |
| 6.42 (triangles) | C2 imine | 0.03 | 0.11 | — | 0.26 | 0.29 (peak) | 0.20 | 0.085 |
| 3.51 (circles) | C2 imine | ~0 | — | — | ~0 | — | ~0 | ~0 |
| 9.29 | C3 product | 0.013 | 0.092 | 0.111 | — | 0.117 (peak) | — | 0.014 |
| 6.42 | C3 product | 0.005 | 0.010 | — | 0.009 | 0.013 (peak) | 0.010 | 0.006 |
| 3.51 | C3 product | ~0 (≤ 0.002, read floor) | | | | | | ~0 |

Sampling times are read from the marker positions (±3 min; the two panels place the last point at 55-60 min); "—" = no marker there.

### 4.2 Tables I-IV (printed values)

**Table I** (sugar or carbonyl + β-Ala, 1 M each, 20 min; browning activity = 1 / time to O.D. 5.0, min⁻¹; relative to
glucose): glucose 95 °C 0.019 (1); fructose 95 °C 0.014 (0.74); xylose 95 °C 0.166 (8.74); xylose 80 °C 0.037 (bridge,
no relative value printed); MGO 80 °C 2.77 (654.3); glyceraldehyde 80 °C 8.33 (1967); glyoxal 80 °C 0.515 (121.6);
glycolaldehyde 80 °C 8.93 (2109). The relative values for the 80 °C rows are scaled through xylose (×8.74 / 0.037) [Z,
checked: 2.77 / 0.037 × 8.74 = 654.3].

**Table II** (1 M each, boiling water bath, 20 min; min⁻¹, relative): glucose + β-Ala 0.028 (1); glucosylamine 0.138
(4.93); Amadori product 0.023 (0.821); Amadori + β-Ala 0.072 (2.57); glucosone + β-Ala 0.666 (23.7); 3-DG + β-Ala 3.85
(137.5).

**Table III** (initial rates from Fig. 1; browning O.D./min, others mM/min):

| condition | browning | Amadori product | C2 imine | C3 imine |
|---|---|---|---|---|
| neutral (6.42) | 0.243* | 9.8 | 0.017 | 0.001 |
| alkaline (9.29) | 0.989** | 14.8 | 0.079 | 0.015 |
| alkaline / neutral | 4.0 | 1.5 | 4.6 | 15 |

\* average over 30 min, ** over 15 min [M].

**Table IV** (N/C of melanoidin): glucose 2 M, pH 2.3, 190 min, O.D. 21: 0.083 ± 0.006; pH 6.5, 53 min, 21: 0.111 ± 0.009;
pH 9.2, 15 / 25 min, 31 / 37: 0.132 ± 0.008; glycolaldehyde 1 M, pH 5.9, 30 min, 183: 0.15.

### 4.3 Arithmetic [Z]

1. **pH slope of the fragment supply (initial rates, 6.42 → 9.29):** C3 ×15 = **0.41 dex per pH unit**; C2 ×4.6 = 0.23;
   Amadori ×1.5 = 0.06; browning ×4.1 = 0.21.
2. **Peak standing levels:** C3 0.117 / 0.013 mM → 0.33 dex per pH unit (6.42 → 9.29). Acid limb (3.51 → 6.42): ≥ 0.28
   dex per pH unit, a lower bound set by my read floor (≤ 0.002 mM), not by an assay LOD.
3. Molar yields at peak: C3 5.9 × 10⁻⁵ mol/mol glucose (pH 9.29), 6.5 × 10⁻⁶ (pH 6.42); C2 3.6 × 10⁻⁴ and 1.4 × 10⁻⁴.
4. Fragments vs Amadori at pH 6.42: initial C3 rate 0.001 vs Amadori 9.8 mM/min, ~10⁻⁴ of the Amadori flux.

## 5. What it can and cannot test in the model

**Can:**
- **The sign and rough steepness of the C2 / C3 fragment supply vs pH** from glucose + an amine near 100 °C: negligible in
  acid, rising ~0.2-0.4 dex per pH unit from neutral to alkaline (initial-rate basis, label pH).
- The relative browning power of fragments (Table I) — no engine observable.

**Cannot:**
- An α-amino-acid or thiol pot: β-alanine is a β-amine (no Strecker degradation with α-dicarbonyls, different pKa); no
  cysteine, no H2S, no MFT.
- Held pH: unbuffered 2 M solutions drift (acids form); only the initial pH is known.
- Absolute fragment rates as constants: the "C3 imine" is a standing level of a reactive intermediate measured through
  an uncalibrated derivatisation; no LOD; 2 M concentrations (~20× Hofmann's sugar).
- pH 3.5 → 6.4 resolution: one acid point, at the floor.

## 6. Conflicts and agreements with the repo's records

1. Agrees in sign with Hofmann 1998 Table 8 (fed HA + MP → MFT rising with pH) and with `yu2017` (MGO from glucose +
   α-amino acids rising 0.5-0.8 dex per pH unit from pH 5 to 7 at 120 °C): fragment chemistry grows with pH.
2. Agrees with `hayashi1980` (same laboratory) that a C2 fragment forms early, before most browning.
3. [!] Fig. 1A shows the alkaline C3 peak at ~30 min while Table III's alkaline browning rate is averaged over 15 min;
   rates and peaks use different windows (author's choice; noted, not a conflict).

## 7. Provenance tiers (first reader only; second reader pending)

- Table I-IV values: `direct_measurement` (printed), confidence medium for Tables III / IV, low-medium for Tables I / II
  (time-to-O.D. ratios across temperatures).
- Fig. 1A values: `direct_measurement`, figure-digitised by eye from an OCR scan, confidence **low** (no replicates, no
  calibration printed, acid series at the floor).
- pH slopes in sec. 4.3: [Z], `derived`.

## M4 role (case 8660c005)

**Role fixed before reading:** band source for R3's fragment-supply pH law ("C2 / C3 from glucose + beta-Ala negligible
in acid, rising with pH"); measurement that would count: "fragment yields at >= 2 pH"; P ~0.5.

- **The named measurement is met:** C2 and C3 at three initial pH in one study.
- **What it supplies.** Coordinate: R3's fragment supply pH law (`k_glc_ha` vs pH; sec. 5.2 band "0 - 1 dex per pH unit,
  rising", "Hayashi 1986 (candidate, by amendment); until then numerical").
  - Numbers: C3 initial rate **0.41 dex per pH unit** (6.42 → 9.29); C3 peak level 0.33 dex/pH; C2 0.23 dex/pH; acid
    limb (3.51 → 6.42) **≥ 0.28 dex/pH** (censored). All inside the band's 0-1 dex/pH and rising, as the band assumes.
  - Temperature / pH: ≈ 100 °C, initial pH 3.51 / 6.42 / 9.29, unbuffered.
  - Transfer: weak. β-alanine (not an α-amino acid, not cysteine), 2 M each, water, label pH only; standing levels of
    a derivatised intermediate, not a formation constant. Hofmann's D2 limb (pH 5 → 7, 145 °C, 0.5 M phosphate, Cys
    33 mM) sits between Hayashi's acid and neutral points, where Hayashi has only the censored lower bound.
  - Use: it supports a rising law of a few tenths of a dex per pH unit; it does not narrow the band's 0-1 range to a
    measured value on the 5 → 7 limb. As a band source it must not be scored or screened in M4 (it is not).
- **Decisive criterion:** none met (no thiol, no hexose MFT, no NF).

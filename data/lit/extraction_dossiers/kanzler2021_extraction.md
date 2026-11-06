# Kanzler, Wustrack & Rohn 2021 (J. Agric. Food Chem. 69:11960-11970) — EXTRACTION (file stem `kanzler2021`)
**DOI:** 10.1021/acs.jafc.1c04594 (printed on every page footer; Crossref 2026-11-05: same title, authors Clemens Kanzler, Felix Wustrack, Sascha Rohn, JAFC 69(40) 11960-11970, online 2021-09-30) — matches the brief.
### *"High-Resolution Mass Spectrometry Analysis of Melanoidins and Their Precursors Formed in a Model Study of the Maillard Reaction of Methylglyoxal with L-Alanine or L-Lysine"*
### Methylglyoxal (MGO) **100 mM** alone or 1:1 with L-Ala or L-Lys (100 mM), water, **pH adjusted to 5.0 ± 0.1** (HCl / NaOH, unbuffered), 2.5 mL sealed glass ampoules, heating block **100 ± 1 °C, 0-300 min** (6 times), triplicates; MGO by OPD / HPLC-DAD. **MGO alone: 36.5 / 21.9 / 11.6 / 9.0 / 7.4 % left at 30 / 60 / 120 / 180 / 300 min — second order, k2 ≈ 5.9 × 10⁻⁴ mM⁻¹ min⁻¹, ×3.2 the k2 of Kanzler 2026's nominally identical pot (same laboratory).** With Ala or Lys: 12.9 / 10.2 % left at 30 min, complete conversion by 120 min. No temperature series.

**Read 2026-11-05 (first reader) from `data/articles/kanzler2021.pdf`** (11 pp., ACS open-access PDF, clean text
layer). Prose and Table 2 from `pdftotext -raw`. **Fig. 1B digitised** from a 400 dpi render of p. 3 (scratchpad
`m4img/kanzler2021_400-03.png`, crop `m4img/k21_f1b.png`): marker centres located with a dark-pixel run finder
(`m4tools/runs.py`) against the gridlines (4.105 px per mol %). **Read error ±1 mol %** (black series); the blue / red
series overlap and are read ±2 mol %. Arithmetic: scratchpad `m4tools/m4_arith.py` (maillard env, `-I`). SI (MS² spectra,
HRMS tables S-1 to S-3, Kendrick plots) not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Clemens Kanzler** (corresponding, clemens.kanzler@tu-berlin.de), Felix Wustrack, Sascha Rohn |
| Affiliation | Technische Universität Berlin, Institute of Food Technology and Food Chemistry, Dept. of Food Chemistry and Analysis, 13355 Berlin |
| Venue | *J. Agric. Food Chem.* **69** (2021) 11960-11970 |
| DOI | 10.1021/acs.jafc.1c04594 |
| Dates | received 29 Jul 2021; revised 20 Sep 2021; accepted 20 Sep 2021; published 30 Sep 2021 |

**Laboratory.** Kanzler group, TU Berlin — **the same laboratory as `kanzler2026`, `bork2022`, `bork2022b`, `bork2023`,
`bork2024`, `bork2024b`** (proposed id `tu_berlin_kanzler`). Its MGO method is cited to Kanzler, Schestkowa, **Haase &
Kroh** 2017 (ref. 29), and it cites Pfeifer, Haase & Kroh 2013 (ref. 27) for MGO aldol dimers [C]: the lineage that
`pfeifer2013` sec. 0 merges with this group. Same nominal protocol as Kanzler 2026 (MGO 40 % Merck, 100 mM, pH 5.0 ± 0.1,
2.5 mL ampoule, Behr block, 100 ± 1 °C, OPD, ProntoSil Phenyl column, 318 nm).

## 1. What this paper answers for the model

| model need | what Kanzler 2021 gives | tag |
|---|---|---|
| MGO loss alone, pH 5, 100 °C, time course | Fig. 1B: 100 / 36.5 / 21.9 / 11.6 / 9.0 / 7.4 mol % at 0 / 30 / 60 / 120 / 180 / 300 min; Table 2: Δc(30 min) = 63.5 ± 4.6 mM; text "7.5 mol %" left at 300 min | [F] / [M] |
| order | cumulative k2 5.8 / 5.9 / 6.4 / 5.6 / 4.2 × 10⁻⁴ mM⁻¹ min⁻¹ (constant to 180 min); cumulative k1 0.034 → 0.0087 min⁻¹ (decelerating) — **second order fits, first order does not** | [Z] |
| MGO + Ala / Lys (100 mM each) | Table 2: Δc(30 min) = 87.1 ± 1.8 (Ala), 89.8 ± 0.9 (Lys) mM → 12.9 / 10.2 % left; Fig. 1B ~5 % at 60 min, ~0 from 120 min; text "complete conversion" at 120 min | [M] / [F] |
| amino acid consumed | **not measured** | — |
| temperature series | **none** (100 °C only) | — |
| colour | Table 2: ΔA420 (30 min) 1.9 ± 0.1 / 19.1 ± 0.9 / 31.7 ± 2.2 (alone / Ala / Lys); Fig. 1A to 300 min; text: 4.5× (Ala) and 6.6× (Lys) the MGO-alone colour at 300 min | [M] |
| products | ESI(−)-Orbitrap HRMS of the 300 min pots: MGO aldol oligomers (dimer m/z 143.0350, trimer 215.0561, tetramer 287.0774, ± H2O), building blocks formaldehyde, acetaldehyde, acetol, and aminoacetone with amino acids (Tables 3, 4; Kendrick and van Krevelen analysis) | [M] |
| SEC (MGO + Ala, 300 min) | coloured fractions at 11.5 / 14.5 / 16.5 min (≥ 800 kDa, ~100 kDa); most material colourless at the permeation limit (RID) | [M] / [Q] |
| glycine, cysteine, NF, MFT | **none** | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (Table 1) | MGO 100 mmol/L; MGO 100 + Ala 100; MGO 100 + Lys 100; MGO / diacetyl 50 + 50 (± Ala 100; HRMS only) |
| pH | adjusted "to (5.0 ± 0.1)" with HCl or NaOH before incubation; **"measured after 300 min" but the value is not printed** in the main text [!]; label pH, unbuffered |
| Vessel / heating | **2.5 mL** sealed in a glass ampoule (ampoule volume not stated), heating block (behr Labor-Technik) at (100 ± 1) °C, 0, 30, 60, 120, 180, 300 min; cooling not described |
| Headspace / O2 | **not stated** (no degassing, no inert gas) |
| Replicates | "Every sample was prepared in triplicate", mean ± SD |

## 3. Assays [M]

- **MGO:** 25 µL sample + 475 µL water [the PDF prints "475 mL", a typo] + 500 µL OPD (50 mM, water / methanol 1:1),
  overnight dark RT; HPLC-DAD 318 nm, ProntoSil 60-5 Phenyl. **Calibration, LOD and LOQ not printed** (method cited to
  Kanzler 2017, ref. 29).
- Colour A420 × dilution; SEC (PL aquagel-OH MIXED-H, pullulan 10 / 100 / 800 kDa); ESI(±)-Orbitrap HRMS by direct
  infusion at 1 mg/mL.
- **Replicate noise:** MGO alone SD ±4.6 mM at 30 min (Table 2; the Fig. 1B bar spans about ±4.7 mol %), ±3 at 60 min,
  ±1.5 later [F]; with amino acids SD ±0.9-1.8 mM at 30 min.

## 4. Data

### 4.1 Fig. 1B — MGO remaining (mol % of 100 mM), pH 5 (label), 100 °C [F ±1; SD]

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| MGO alone | 100 | **36.5** (4.6; Table 2) | 21.9 (~3) | 11.6 (~1.5) | 9.0 (~1.5) | 7.4 (~1) |
| MGO + Ala | 100 | **12.9** (1.8; Table 2) | ~5 | ≤ 1 | ≤ 1 | ≤ 1 |
| MGO + Lys | 100 | **10.2** (0.9; Table 2) | ~4 | ≤ 1 | ≤ 1 | ≤ 1 |

Text: "fell below 10 mol %" at 30 min in both amino-acid pots [M] — [!] Table 2 gives 12.9 % for Ala, above 10 %; the
figure agrees with the table.

### 4.2 Table 2 (verbatim values), 30 min

| system | ΔA | Δc (mmol/L) | ΔA/Δc (1/mol/L) |
|---|---|---|---|
| methylglyoxal | 1.9 ± 0.1 | 63.5 ± 4.6 | 30 ± 4 |
| methylglyoxal / alanine | 19.1 ± 0.9 | 87.1 ± 1.8 | 219 ± 15 |
| methylglyoxal / lysine | 31.7 ± 2.2 | 89.8 ± 0.9 | 353 ± 28 |

### 4.3 Arithmetic [Z] (`m4_arith.py`)

1. **MGO alone, 100 mM, pH 5, 100 °C.** Cumulative k2 = 5.80 / 5.94 / 6.35 / 5.62 / 4.17 × 10⁻⁴ mM⁻¹ min⁻¹ (0-30 …
   0-300); **k2 ≈ 5.9 × 10⁻⁴ (0-180 mean), 4.8-7.1 × 10⁻⁴ across Table 2's ±1 SD at 30 min**. Cumulative k1 0.0336 →
   0.0087 min⁻¹. The 300 min point falls below the second-order line (k2 4.2 × 10⁻⁴): a slower tail.
2. **Amino acids, k1(0-30):** Ala 0.068, Lys 0.076 min⁻¹ vs MGO alone 0.034: ×2.0 (Ala), ×2.3 (Lys).
3. **Against Kanzler 2026 (same laboratory, same nominal pot):** MGO alone left at 30 min 36.5 % vs 65 % (SD 10);
   k2 5.9 vs 1.85 × 10⁻⁴ (**×3.2, 0.51 dex**). MGO + Ala at 30 min 12.9 % vs 32 % (×2.5 in remaining); MGO + Lys 10.2 vs
   11 % (agree).
4. **Cross-paper MGO-alone Ea with `bork2022`** (23 mM, 130 °C, pointwise k2 2.84 → 1.84 × 10⁻³): from this paper's
   100 °C k2, **Ea ≈ 47-66 kJ/mol**; from Kanzler 2026's, 96-114 kJ/mol. Same laboratory throughout; the choice of the
   100 °C paper moves the estimate by ~50 kJ/mol.

## 5. What it can and cannot test in the model

**Can:**
- The size and order of the amine-free MGO loss at pH 5 / 100 °C, with six time points (Kanzler 2026 has two).
  Second order identified within this series (k2 constant to 180 min within ±7 %).
- The engine's `k_mgo_mel` (first order, amine-free, fitted 0.0228 min⁻¹ at 100 °C): sits inside this paper's
  MGO-alone pseudo-first-order range (0.034 → 0.0087 min⁻¹ across 0-300 min at 100 mM), i.e. the fitted level at 100 °C
  is not contradicted by this pot. The order is (second, not first).
- Amine acceleration of MGO loss with Ala and Lys at 30 min (×2.0-2.3), with Lys and Ala equal in MGO loss but Lys ×1.7
  in colour.

**Cannot:**
- **Any temperature dependence** (100 °C only); so it cannot test `k_mgo_mel`'s Ea directly.
- Amino-acid consumption (not measured), glycine, cysteine, NF, MFT.
- pH drift: measured at 300 min but not reported; O2 unknown.

## 6. Conflicts and agreements with the repo's records

1. **[!] CONFLICTS with `kanzler2026` (same laboratory, same nominal protocol).** MGO alone at 30 min: 36.5 ± 4.6 % (this
   paper) vs 65 ± 10 % (Kanzler 2026 Fig. 1D); at 120 min 11.6 vs 31 ± 6.5 %. k2 ×3.2. Outside both papers' SDs. MGO + Ala
   at 30 min 12.9 ± 1.8 vs 32 ± 1.5 %. MGO + Lys agrees (10.2 vs 11 %). Possible causes (none stated by either paper):
   different MGO stock batches (40 wt % MGO solutions contain oligomers and varying free MGO), calibration, or fill /
   headspace. **Within one laboratory, the between-paper spread of the MGO-alone k2 is 0.5 dex.** Proposed CF text
   (not registered; `docs/reference/` out of scope): "MGO alone, 100 mM, pH 5, 100 C, TU Berlin: Kanzler 2021 k2 5.9e-4
   vs Kanzler 2026 k2 1.85e-4 /mM/min (x3.2), same protocol".
2. **[!] Bears on CF-122 (branch `registers-2026-11-05`).** CF-122 quotes "~95 - 113 [kJ/mol] only as a cross-paper
   estimate with Bork 2022" for MGO alone; with this paper's 100 °C value the same estimate is 47-66 kJ/mol. The
   amino-acid-pot Ea of 95-121 kJ/mol (Kanzler 2026 Fig. 2D, one incubation set) is untouched. Both readings remain far
   above the engine's 20 kJ/mol.
3. Agrees with `kanzler2026` and `pfeifer2013` on products (aldol oligomers, formaldehyde / acetaldehyde / acetol
   building blocks) and with Kanzler 2026 that amino acids accelerate MGO loss.
4. [!] Minor: "475 mL" water in the OPD step (µL meant); "below 10 mol %" vs Table 2's 12.9 % for Ala.

## 7. Provenance tiers (first reader only; second reader pending)

- Fig. 1B / Table 2 MGO values: `direct_measurement` (Table 2 printed; Fig. 1B figure-digitised), external calibration
  (method cited), confidence **medium** (triplicates with SDs; calibration range and LOD not printed; between-paper
  conflict, sec. 6.1).
- k2, k1, ratios and Ea in sec. 4.3: [Z], tier `derived`; not for write-back.
- Proposed ledger entry (trunk.yml; NOT written): `kanzler2021_fig1b_methylglyoxal_alone_and_with_ala_lys_ph5_100c`,
  lab "TU Berlin (Kanzler / Rohn)", MGO 100 mM ± Ala / Lys 100 mM, water, pH 5.0 ± 0.1 initial, 2.5 mL ampoule, 100 °C,
  t = 0-300 min, species [methylglyoxal], kind time_series, location "Fig. 1B p. 11962; Table 2", quant_class calibrated,
  flags [ph_start_only, unbuffered, o2_not_stated, no_lod, conflicts_kanzler2026_same_lab].

## M4 role (case 8660c005)

**Role fixed before reading:** band source (MGO + Lys conversion; same laboratory as Kanzler 2026, so it widens no
laboratory count); measurement that would count: "MGO conversion vs time or T"; P ~0.3.

- **What it supplies.** Coordinate: the MGO loss law (sec. 5.2 band "Ea 95-121 kJ/mol, level at 100 C pinned to B60's
  0.0228 /min; order 1 (Kanzler's apparent order ~2 is a sensitivity)").
  - MGO conversion vs time: **yes** (six points, triplicates), at 100 °C, label pH 5, unbuffered, 100 mM. MGO alone
    second order, k2 ≈ 5.9 × 10⁻⁴ mM⁻¹ min⁻¹; MGO + Lys 10.2 % left at 30 min (k1(0-30) 0.076 min⁻¹), MGO + Ala 12.9 %
    (0.068 min⁻¹). Transfer: no thiol, no glycine, 100 mM (the trunk's pots are ~1000× more dilute, as CF-122 notes);
    it informs the order (second, now identified in one series) and the 100 °C level, not Ea.
  - vs T: **no**.
- **Unexpected, and material to the band (not decisive):** the same laboratory's nominally identical MGO-alone pot loses
  MGO ×3.2 faster (k2) than Kanzler 2026's. The MGO-alone half of the band's basis ("~95 - 113 kJ/mol for MGO alone
  against Bork 2022") becomes **47-66 kJ/mol** if this paper's 100 °C point is used instead. The amino-acid-pot basis
  (95-121 kJ/mol, Kanzler 2026 Fig. 2D) is unaffected. The band's required-agreement sensitivity at B60's Ea 20 stays
  necessary; a sensitivity near 50 kJ/mol now has a measured basis too. The within-laboratory between-paper spread
  (0.5 dex in k2) is the size of omega (0.419 dex) for this coordinate.
- **Decisive criterion:** none met, none expected (no thiol, no hexose MFT).

# Kreitman, Danilewicz, Jeffery & Elias 2016 (Part 2) — EXTRACTION (Fe(III) alone and Fe(III) + Cu(II) oxidation of H2S, cysteine and hexanethiols in model wine)

**Source on disk:** `data/articles/kreitman2016b.pdf`, the ACS ASAP version (9 pages lettered A–I, with volume and
pages printed as "XXXX, XXX, XXX−XXX"). The text layer was extracted with `pdftotext -layout`, and pages C–E
(Figures 4 and 5 and the results text) were read and checked by eye from the page images (Figs 4A and 5B again at
300 dpi) on 2026-10-09. Page letters below are the printed ones. Companion dossiers:
`kreitman2016_extraction.md` (Part 1, Cu alone) and `ehrenberg1989_extraction.md` (the rate law).

| field | value |
|---|---|
| Title | "Reaction Mechanisms of Metals with Hydrogen Sulfide and Thiols in Model Wine. Part 2: Iron- and Copper-Catalyzed Oxidation" |
| Authors | G. Y. Kreitman, J. C. Danilewicz, D. W. Jeffery, R. J. Elias |
| Venue | J. Agric. Food Chem. 2016 (received 5 Feb 2016, revised 29 Apr 2016, accepted 30 Apr 2016, p. A; volume and pages 64:4105 as supplied, not printed in this PDF) |
| DOI | 10.1021/acs.jafc.6b00642 (printed) |

## 1. Methods (pp. B–C)

- Model wine as in Part 1: tartaric acid 5 g/L, ethanol 12 % v/v, **pH 3.6**. Air-saturated; 60 mL BOD bottles
  overfilled and stoppered (**no headspace**), dark, **"ambient temperature"** (no number printed). Triplicates,
  sacrificial bottles.
- Thiols 300 µM: H2S, Cys, 6SH, 3SH. **Fe(III) as FeCl3·6H2O.** Treatments:
  - Fe alone: Fe(III) 200 µM.
  - Fe + Cu: Fe(III) 200 µM + Cu(II) 50 µM (H2S, 6SH, 3SH); for **Cys halved to Fe(III) 100 µM + Cu(II) 25 µM**
    because at 200/50 "Cys reacted rapidly and was completely consumed within 5 min (data not shown)" (pp. D–E).
  - Thiol + H2S: H2S 50 or 100 µM with Fe(III) 200 + Cu(II) 50 µM; Cys + H2S 50 µM also at 100/25.
- Fe(III) speciation from the Fe(III)–tartrate absorbance at 336 nm (Fe(II)–tartrate does not absorb). Thiols by
  Ellman's; mixed systems by bimane HPLC-MS/MS; O2 by oxidots; AC by DNPH-HPLC.
- No chelator treatment; no temperature or pH series.

## 2. Findings that matter

### 2a. Fe(III) alone, 200 µM (p. D, Fig. 4)

| thiol | consumed | Fe(III) reduced | O2 / µM | O2 : thiol (as printed) |
|---|---|---|---|---|
| H2S | 262 µM after 144 h | up to ~66 % within 96 h | 135 | 1:1.9 raw; ~1:1.5 after subtracting ~66 µM H2S held by the Fe(II) |
| Cys | **192 µM after 193 h**; "incomplete after 200 h" (p. E) | max ~17 % within 24 h, then a steady state | 49 | ~1:3.2 (157 µM Cys after subtracting 35 µM for the 17.5 % Fe(II) left) |
| 6SH, 3SH | < 15 µM | minimal | minimal | not calculable |

- No initial uptake with Fe(III), unlike Cu(II)'s immediate 2 equiv (p. C).
- Thiols do not displace tartrate from Fe(III); Cys "can presumably compete for Fe" through its carboxylate, and
  6SH/3SH cannot, which is the stated reason they barely react (p. D).
- AC 15–30 µM in the Cys and H2S systems (data not shown).
- Fig. 4A, Cys, read from graph, approx.: ~302 µM at 0 → ~295 at 2 h → ~275 at 8 h → ~242 at ~24 h → ~221 at
  ~48 h → ~185 at ~96 h → ~110 at ~193 h. Average ≈ 1 µM/h (derived here: (302 − 110)/193).

### 2b. Fe(III) + Cu(II) (pp. D–F, Figs 5, 6)

| thiol | metals / µM | Fe(III) reduction | thiol course | O2 / µM | O2 : thiol (as printed) | AC / µM |
|---|---|---|---|---|---|---|
| Cys 300 | Fe 200 + Cu 50 | — | "completely consumed within 5 min" (data not shown) | — | — | — |
| Cys 296 | **Fe 100 + Cu 25** | almost fully within 5 min (< 5 % of 336 nm absorbance left) | ~150 µM in the rapid phase, **complete within 7 h** | 110 | **1:2.6** (284 µM after subtracting 12 µM for the ~12 % Fe left reduced) | 60 (O2:AC ≈ 2:1) |
| 6SH 273 | Fe 200 + Cu 50 | initial ~40 %, settling at ~25 % Fe(II) | fully oxidised within 7 h | 106 | ~1:2.1 (223 µM) | 146 (≈ 1:1) |
| 3SH | Fe 200 + Cu 50 | as 6SH | incomplete at 150 h; 267 µM consumed | 82 | ~1:2.6 (217 µM) | — |
| H2S | Fe 200 + Cu 50 | near complete within 30 min | sharp ~135 µM drop; 308 µM consumed by ~48 h | 160 | ~1:1.6 (262 µM after subtracting 46 µM for 92 µM Fe(II) at 120 h) | 100 (≈ 1.6:1) |

- Mechanism (p. E, Fig. 6): Cu(II) is reduced by the thiol to the Cu(I)–SR complex; **Fe(III) oxidises Cu(I)
  rapidly**; Fe(II)–tartrate is reoxidised by O2 ("known to be fast"); "with copper alone, overall thiol oxidation
  is dependent upon the rate of reaction of O2 with the Cu(I) complex; however, when iron is present, the reaction
  rate is dependent upon the oxidation rate of the Fe(II)–tartrate complex". The Cu(I)–SR aggregate "reacts more
  slowly with O2 than with Fe(III)" (p. F).
- The authors' accounting of the Cys 100/25 rapid phase (p. E): 25 µM Cu(I) complex formed at once, Cu recycled
  "a further 3 times" to reduce ~100 µM Fe(III), then 25 µM Cys to cystine and 25 µM bound to Cu(I): "150 µM
  Cys would be consumed when all Fe(III) and Cu(II) were reduced", with no O2 consumed yet.
- Fig. 5B, Cys, read from graph, approx.: ~296 µM → ~150 at ~5 min → ~127 at 0.5 h → ~83 at 2 h → 0 at 7 h.

### 2c. Thiol + H2S mixtures (pp. F–G, Figs 7, 8)

- At least 60 % of free H2S was removed within 5 min in all four mixtures; by 24 h there was virtually none.
  Thiols kept oxidising afterwards without a copper precipitate.
- Cys + H2S: high metals (Fe 200 + Cu 50) oxidised all H2S and Cys within 2 h; low metals (100 + 25) needed 24 h.
  Total Cys + H2S consumed 302 and 326 µM; O2 132 and 138 µM; **~1:2.3 O2 : (Cys + H2S)** under both; AC 150 µM
  (high) and 81 µM (low).
- Mixed organic polysulfanes with up to five S atoms (6SH; similarly 3SH); S5-bimane from H2S alone (p. G).
- Context (pp. B, F, G–H): wine Fe is ~10-fold over Cu ("5−10-fold higher"); wine H2S 0.3–1 µM, GSH up to 40 µM,
  Cys and analogues 20 µM; Zn, Al, Mn average 0.54, 0.41, 0.97 mg/L.

## 3. What it means for the model

- **Rate constants: none**, as in Part 1: one pH (3.6), one unstated ambient temperature, no headspace, metals at
  25–200 µM.
- **The Fe–Cu interaction is the main new fact.** Fe alone oxidises Cys slowly (≈ 1 µM/h at 200 µM Fe, pH 3.6;
  derived here from Fig. 4A). Cu with Fe is far faster than either: Cu is reduced by the thiol, Fe(III) reoxidises
  Cu(I), and O2 reoxidises Fe(II). At pH 3.6 the rate is set by Fe(II)–tartrate oxidation. Ehrenberg 1989 found
  the opposite sign at pH 7.2 (0.1–1 µM Fe cuts the 1 µM Cu rate to about 40 %, and 10 µM Fe raises it again).
  **The sign of the Fe–Cu interaction depends on pH, ligand (tartrate) and the Fe:Cu ratio.** Neither paper
  covers a pH-5 phosphate pot.
- **Caution on "complete within 5 min"** (derived here): Fe(III) 200 µM can take 200 µM Cys by one-electron
  reduction, and Cu(II) 50 µM takes 100 µM (2 per Cu, Part 1). Together that is 300 µM, the whole charge. So the
  5-min disappearance at 200/50 is consistent with **stoichiometric metal reduction**, not demonstrated catalytic
  turnover. The authors' 150 µM accounting at 100/25 is the same arithmetic. The catalytic rate at 100/25 is
  the second phase: ~150 µM in ≤ 7 h, ≥ ≈ 0.36 µM/min (derived here, a lower bound from one 7-h point), i.e.
  ≥ ≈ 0.014 thiols per Cu per minute. That is roughly ten times Part 1's Cu-alone turnover (≈ 1–2×10⁻³).
- **Stoichiometry for any implementation**: O2:Cys 1:3.2 (Fe alone), 1:2.6 (Fe + Cu), 1:4.5 (Cu alone, Part 1),
  1:2.3 for Cys + H2S. The lower ratios include ethanol oxidation via Fenton (AC formed), which a buffer pot without
  ethanol lacks.
- **Live engine values** (keys and arithmetic in `ehrenberg1989_extraction.md` sec. 3): `k_cys_thermal` log10
  k(145 °C) = −2.066359667088082 /min, Ea 55.1 kJ/mol; `k_cys_h2s` Zheng & Ho; `k_cys_ox` = 0 (inert). The engine
  carries no Fe or Cu species, so neither the Cu channel nor the Fe–Cu cycle can be represented. A metal-catalysed
  Cys channel would need a catalyst-concentration input, and, if Fe matters, a second catalyst whose effect is
  not even of fixed sign across the two papers.
- **For EXPERIMENTS.md arm D**: the papers give no chelator data in wine. Ehrenberg's EDTA result (1:1 no effect,
  10:1 abolished; Fe on reaction 2 only) is the only chelator evidence in the three.

## What it does not give

- No rate constants, orders, Michaelis constants or temperature dependence; "ambient" has no number.
- One pH (3.6) in tartrate/ethanol; Fe speciation is tied to tartrate, so the Fe results do not carry to
  phosphate.
- No Cu concentration below 25 µM and no Fe series; no chelator arm.
- The Cys result at Fe 200 + Cu 50 is "data not shown".

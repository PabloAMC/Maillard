# Elfsson, Wallin, Eksborg, Rudaeus, Ros & Ehrsson 1998/1999 — EXTRACTION
### *"Stability of 5-aminolevulinic acid in aqueous solution"*
### 5-Aminolevulinic acid (ALA, H2N-CH2-CO-CH2-CH2-COOH, an alpha-aminoketone) disappearing by SECOND-ORDER self-condensation (to the dihydropyrazine CHPY, and some porphobilinogen) in buffered water: pH 2.35 / 4.81 / 7.42-7.53, 50 °C (pH series) and 37-85 °C (pH 7.53, Ea), 0.031-0.074 M, ± EDTA. Ion-pair RP-HPLC of the parent only. **The only measured alpha-aminoketone self-condensation rate constant on disk; an analogue for aminoacetone + aminoacetone, not a pyrazine measurement.**

**Read 2026-10-31 from `data/articles/elfsson1999.pdf`** (5 pp., born-digital, clean text layer). Prose and numbers
from `pdftotext -layout`; Figs. 1-2 (a chromatogram, one 1/[ALA] vs t plot) not rendered: no number is taken from
them. Every derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | B. Elfsson, I. Wallin, S. Eksborg, K. Rudaeus, A.M. Ros, H. Ehrsson (Karolinska Pharmacy / Karolinska Hospital, Stockholm) |
| Venue | *Eur. J. Pharm. Sci.* **7**, 87-91. **The PDF prints "(1998)"** in the running header; received 31 Oct 1997, accepted 13 Jan 1998 [!] |
| DOI | PII S0928-0987(98)00009-8 printed; Crossref 2026-10-31: **10.1016/S0928-0987(98)00009-8**, same title, vol. 7, pp. 87-91, issued **1999** (issue date). The stem's "1999" follows Crossref; the article says 1998 |

## 1. What this paper answers for the model

Context: the engine's `k_cond` (`src/kinetic_core/parameters_pyrazine.py`) is **declared fast**, 1.0e4 L/(mmol·min)
= 1e7 M⁻¹ min⁻¹, no pH term, no Ea, shared by the three pairings; its justification is Jousse 2002 calling the step
"fast" and Zhou 2024's pyrazines rising linearly from t = 0 (figure only).

| need | what Elfsson gives | tag |
|---|---|---|
| condensation rate law | **second order in ALA**: −d[ALA]/dt = k2[ALA]²; t1/2 vs 1/[ALA]0 linear over 0.52-1.24 % (0.031-0.074 M), y = 0.219x − 0.589, r² 0.9645, intercept CI −2.347 to 1.169 (includes 0) | [M] |
| k at one T, pH | 50 °C: **pH 2.35 no loss in 37 days**; **pH 4.81 k2 = 0.069 M⁻¹ h⁻¹** (SE 6.3e-3, n 9, [ALA]0 0.057 M, t1/2 257 h); **pH 7.42 k2 = 5.53 M⁻¹ h⁻¹** (SE 0.58, n 7, [ALA]0 0.060 M, t1/2 3.0 h) | [M] |
| pH law | k(7.42)/k(4.81) = **80** → **0.73 dex per pH unit** over 4.8-7.4 [Z]. One free amine required (pKa2 8.05 at 25 °C) predicts 330 (or ~200 if the amine pKa falls to ~7.4 at 50 °C, an assumed −0.025 /K) [Z]; two free amines ~1e5 [Z]. Observed is BELOW the one-free-amine law. **Buffer confounded with pH** (0.9 M acetate vs 1.15 M phosphate) [!] | [M] / [Z] |
| Ea | **43.7 kJ/mol** (pH 7.53, 0.06 M, 37-85 °C, n 7, r² 0.9927; slope CI −2.507 to −2.056); "about 1.5 times" per 10 °C [Q]. Printed line y = −2.281x + 4.217 reproduces k2(50 °C) = 0.086 M⁻¹ min⁻¹ = 5.2 M⁻¹ h⁻¹ if x = 1000/T and k2 in M⁻¹ s⁻¹ (printed 5.53 at pH 7.42) [Z] | [M] / [Z] |
| k at process temperatures (extrapolated beyond 85 °C [!]) | pH 7.5: **0.76 M⁻¹ min⁻¹ at 100 °C, 1.6 at 120 °C, 2.2 at 130 °C** (ALA consumption) [Z]; the pair constant (rate of DHP formation = k2/2 [ALA]²) is half | [Z] |
| products | none measured. Cited: CHPY (2,5-bis(β-carboxyethyl)dihydropyrazine) is the major product and porphobilinogen (a pyrrole) the minor, ratio ~1:10 "largely independent" of conditions (Franck & Stratmann 1981); CHPY oxidises to the pyrazine CPY; 4-7 % of ALA may exist as dihydropyrazine at physiological conditions (Jaffe 1990) | [C] |
| catalysis | EDTA 1 % no effect (k2 4.83 vs 5.53 M⁻¹ h⁻¹, P = 0.3) — no metal catalysis | [M] |
| 13C / 15N, peptides, mixing ratio | none | — |

**Answers.**

1. **A measured alpha-aminoketone self-condensation is ~1e7 times slower than the engine's declared `k_cond`** [Z]:
   ~0.4 M⁻¹ min⁻¹ (pair constant, pH 7.5, 100 °C, extrapolated with Ea 43.7) vs 1e7 M⁻¹ min⁻¹. ALA's reacting unit
   (H2N-CH2-CO-CH2R, amine pKa 8.05) is the same as aminoacetone's (R = H); the remote carboxylate is the main
   difference. Aminoacetone's alpha-ammonium pKa is lower, **6.8** (`dutra2001_extraction.md`, written in the same
   ingestion; a lead until two-reader checked), so at pH 7.5 its free-amine fraction is ~0.83 vs ALA's ~0.2-0.5:
   a factor ~2-4 faster per free amine, nowhere near 1e7 [Z]. My guess that aminoacetone self-condenses within 10x of
   ALA under the same conditions: ~0.6 (no aminoacetone condensation rate on disk; Dutra measures its O2 oxidation,
   k2 0.160 M⁻¹ s⁻¹ at 37 °C, a competing sink). Aminoacetaldehyde (an aldehyde, the engine's AKG) is not covered.
2. **Consequence for the engine's steady-state assumption** [Z]: at ALA's rate, an aminoacetone pool fed at Zhou
   2024's DMP-equivalent flux (2 x 0.0035-0.023 µmol L⁻¹ min⁻¹, 100-120 °C) relaxes to steady state in
   τ = 1/(2 sqrt(flux · k2)) ≈ 115 h (100 °C, 0.007 µM/min, k2 0.76) to 31 h (120 °C, 0.046 µM/min, k2 1.6), so DMP would rise quadratically with a long lag over a 2 h run, not linearly.
   If Zhou's linear rise is real, either aminoacetone condenses >= ~1e3x faster than ALA, or DMP does not form by
   free-aminoacetone dimerisation. Both are live; this paper cannot decide. It does make "declared fast, any value
   above about 1e3 gives the same yields" an assumption with a measured analogue against it.
3. **Pyrrole vs pyrazine** [C]: with no other carbonyl present, the alpha-aminoketone gives DHP over pyrrole ~10:1
   (Franck & Stratmann, cited). Adams 2008's pyrrole 34 dominance (Ala + MGO) needs excess MGO
   (`adams2008_extraction.md`); a slow self-condensation is what lets an MGO-consuming pyrrole sink win. Consistent
   with EX-052 S5b's direction (pyrrole / DMP rises with MGO), not a test of it.
4. **pH law of the condensation itself**: ~0.7 dex / pH between 4.8 and 7.4 [Z]. The engine puts all of its pH
   dependence on the two Strecker steps; if condensation is rate-relevant, part of Leahy's pH slope belongs to it.
   Single lab, two pH points, buffer-confounded: weight low.

## 2. Conditions [M]

| item | value |
|---|---|
| Substrate | ALA·HCl (Sigma, 98 %); 1 % = 60 mM; concentration series 0.52-1.24 % (0.031-0.074 M) at pH 7.42 |
| Buffers | pH 2.35: 0.3 M phosphate; pH 4.8: 0.9 M acetate; pH 7.4-7.5: **1.15 M phosphate**; pH set with 1 M NaOH |
| pH (label vs measured) | measured at 25 °C before and after: drift < 0.03 (pH 4.81), < 0.06 (7.42), < 0.1 (7.53, T series) |
| Temperature | 50 °C (pH, concentration, EDTA); 37-85 °C (pH 7.53, seven temperatures) |
| Vessel / O2 | thermostatted water bath; vessel, headspace, deaeration not stated [!] |
| pKa | potentiometric, 25.0 °C, 0.1 M NaCl, n 15: pKa1 3.90 (CI 3.76-4.03), pKa2 8.05 (7.88-8.30); slight degradation during titration admitted [Q] |
| Analysis | aliquot diluted 1:50 in 0.05 M H3PO4 (quench), C18 SPE clean-up, ion-pair RP-HPLC (octanesulfonate), UV 263 nm; linear 0.4-2.2 mM; degradation products shown not to interfere |
| Kinetic fit | 1/[ALA] vs t; k2 for the first ~50 % of loss in the T series |

## 3. Quantitative results [M]

| condition | k2 | n | [ALA]0 | t1/2 |
|---|---|---|---|---|
| pH 2.35, 50 °C | no degradation in 37 days | | 1 % | |
| pH 4.81, 50 °C | 0.069 M⁻¹ h⁻¹ (SE 6.3e-3) | 9 | 0.057 M | 257 h |
| pH 7.42, 50 °C | 5.53 M⁻¹ h⁻¹ (SE 0.58) | 7 | 0.060 M | 3.0 h |
| pH 7.41, 50 °C, + 1 % EDTA | 4.83 M⁻¹ h⁻¹ (SE 0.30) | 8 | 0.062 M | |
| pH 7.53, 37-85 °C | Ea 43.7 kJ/mol; log k2 = −2.281 (1000/T) + 4.217 (k2 in M⁻¹ s⁻¹ by my check) | 7 | 0.06 M | |
| concentration series, pH 7.42, 50 °C | t1/2 (h) = 0.219 (1/[ALA]0) − 0.589, r² 0.9645 | 5 | 0.031-0.074 M | |

Predictions printed [M]: shelf-life t90 of 1 % ALA at pH 7.53, 20 °C = 1.9 h; of 10 % = 10 min.

## 4. Mechanism evidence and its kind

- **Kinetic order** (kind 6, within study: concentration varied) — second order in ALA, consistent with a bimolecular
  condensation as the rate-determining loss. Narrow range (a factor 2.4), chosen to limit pH drift [Q].
- **pH** (kind 6): protonated amine is unreactive (pH 2.35 stable) [M]; the authors' argument that the free amine
  attacks a neighbour's ketone is [Q].
- **Products**: not measured here; product identities [C] (Franck & Stratmann 1981; Butler & George 1992: CHPY under
  anaerobic moderate pH, oxidised to CPY; pseudo-porphobilinogen in strong alkali).
- Fit to parent decay only (kind 7 in the ranking); no mass balance. Lab: Karolinska (pharmaceutics), independent of
  every food-chemistry group in the corpus.

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| `k_cond` "declared fast" (B18) | **CONTRADICTS (analogue)** | the one measured alpha-aminoketone self-condensation is ~1e7x slower; its pH and Ea are measurable, not zero |
| EX-052 "one shared condensation constant probably wrong" | SUPPORTS (indirectly) | if condensation is slow, pair constants and pool lifetimes matter, and AKG (aldehyde) vs AKM (ketone) need not share a constant (P1x's premise) |
| P1x / P1y vs P2b | **does not distinguish** | ketone-ketone only; nothing on the cross pair or on an extra route |
| triage 5 / evidence map 2.7 "aminoketone time course: NONE" | LEAVES UNCHANGED for aminoacetone; gives an analogue | first aminoketone lifetime number on disk |
| pyrrole sink (S5, Adams) | SUPPORTS the direction [C] | DHP dominates without excess dicarbonyl; a slow DHP route leaves room for an MGO-driven pyrrole sink |
| CF-077 | LEAVES UNCHANGED | |

## 6. Disagreements found

1. **Engine `k_cond` = 1e7 M⁻¹ min⁻¹ (declared) vs ALA ~0.4 M⁻¹ min⁻¹ (pair constant, pH 7.5, 100 °C,
   extrapolated)** — 7 decades. Candidate causes: different substrate (ALA vs aminoacetone / aminoacetaldehyde;
   aldehyde partner much faster), different pH (7.5 vs 8 in Zhou's pots: ~+0.4 dex at 0.73 dex/pH), extrapolation
   beyond 85 °C, and the engine's step being a lumped "net" step whose speed was never measured.
2. **Year**: article header 1998, Crossref 1999 (issue date). Cosmetic, but the stem says 1999.
3. **Printed Arrhenius line vs printed k2**: the line gives 5.2 M⁻¹ h⁻¹ at 50 °C, pH 7.53, vs 5.53 printed at pH 7.42
   (the higher pH should be faster). Within the SE (0.58); not a real conflict.

## 7. Ledger-ready measurements

1. k2 (pH 4.81, 7.42; 50 °C) and the pH 2.35 null — a two-point pH law for an alpha-aminoketone condensation,
   tagged `analogue_ALA`, buffer-confounded.
2. Ea 43.7 kJ/mol with the printed regression (37-85 °C, pH 7.53) — the only aminoketone-condensation Ea on disk.
3. Second-order evidence (t1/2 vs 1/[ALA]0) — supports the bimolecular form the engine uses.

Lead (no DOI checked): a measured aminoacetone self-condensation rate (to 2,5-dimethyl-dihydropyrazine) would
replace this analogue. Dutra 2001 (Bechara group, `dutra2001_extraction.md`) measures aminoacetone OXIDATION only;
the same group's aminoacetone papers are the place to look.

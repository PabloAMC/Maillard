# Saunders & Jervis 1966 — EXTRACTION (buffer salts as a THIRD reactant in glucose / invert sugar + glycine browning: 0.25 M sodium phosphate pH 7.5, sodium citrate 0.25-0.5 M at pH 3.0 and 8.0, phosphate / citrate / acetate at pH 3.0; sealed jars 2 h at ~121 °C, one reflux pair at 100 °C; browning (T490), carbonyls (A285), CO2, pH change)
DOI: 10.1002/jsfa.2740170604 (not printed on the paper; Crossref returns "The rôle of buffer salts in non-enzymic browning", Saunders, Jervis, J. Sci. Food Agric. 17(6) 245-249, June 1966; checked 2026-11-01)

**Read 2026-11-01 from `data/articles/saunders1966.pdf`** (5 pp., J. Sci. Fd Agric. 17:245-249; an Acrobat Paper
Capture scan with an OCR layer). All five pages read from the page images; **Figs. 1-5 read from 200 dpi renders**
(scratchpad `ing6/saund-3.png`, `saund-4.png` and crops). Read error on the figures: ±2 % transmittance, ±0.02
absorbance, ±8 mg CO2 / 100 ml, ±0.1 pH unit. Table II is printed and was read from the image.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown · **[Z]**
derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Title | "The rôle of buffer salts in non-enzymic browning" |
| Authors | J. Saunders, F. Jervis — Cerebos Group Research Limited, West Mills, Colinton, Edinburgh 13 |
| Venue | J. Sci. Fd Agric. 1966, **17** (June), 245-249; received 12 July 1965, amended 15 October 1965 |
| Paper type | factorial mixture designs (two of three components varied at constant total, or equimolar ladders), single end-point after heating; no rate constants ("no attempt has been made to calculate reaction constants" [Q]) |
| Lab | Cerebos (Edinburgh). **Not in `docs/reference/contrast_benchmark/lab_ids.yml`**: a new laboratory, no chain to any listed group |

## 1. What this paper answers for the model

The question (E5 in `docs/reference/pathway_atlas/ENGINE_VS_LITERATURE.md`; CF-011): buffer anion vs rate at equal pH,
concentration dependence (saturation?), temperature; is phosphate special; is the effect confounded with pH or
metals?

| what the model needs | what Saunders & Jervis give | tag |
|---|---|---|
| buffer vs no buffer, same start pH, amine present | 0.25 M glucose + 0.25 M glycine, start pH 7.5, 2 h ~121 °C: **with 0.25 M phosphate T490 ≈ 41 %, without ≈ 89 %** (100-fold dilution) → A490 0.387 vs 0.051, **×7.6 (range ×6-10)**; A285 0.50 vs 0.07 (×7); CO2 ≈ 130 vs ≈ 25 mg/100 ml (×5) | [F] / [Z] |
| the same at 100 °C, dilute | 0.05 M glucose + 0.05 M glycine, start pH 8.05, reflux: **0.2 M phosphate ×8-13 on A490 at every time 2-11.5 h** (2 h: T 85 vs 98.5 %; 11.5 h: 46 vs 92.5 %, 10-fold dilution) | [F] / [Z] |
| buffer vs no buffer, acid (pH 3.0), amine present | 0.5 M invert sugar + 0.5 M glycine (Table II): A490 ×**2.0** (0.25 M phosphate), ×**1.8** (0.25 M citrate), ×**1.8** (0.75 M acetate); A285 0.497 → 0.895 / 1.080 / 0.866 | [M] / [Z] |
| amine-free | 0.5 M invert sugar, pH 3.0 (Table II): A490 ×**11** (phosphate 0.25 M), ×**12** (citrate 0.25 M), ×**12** (acetate 0.75 M) over no buffer; no CO2 without glycine except citrate (2.8 mg/100 ml) | [M] / [Z] |
| is phosphate special? | **No, at pH 3**: phosphate ≈ citrate per mole, acetate needs ~3× the molarity for the same effect (Table II). Phosphate / citrate mixtures at 0.25 M total, pH 7.5: "not strictly equal for a given molarity, they were almost exactly additive" [V, System E, no figure or numbers printed] | [M] / [Q] |
| phosphate vs citrate per mole near neutral | only across two different systems: 0.25 M phosphate pH 7.5 ×7.6 (direct control) vs 0.25 M citrate pH 8.0 ~×3-4.6 (control INFERRED from the constant-total arms, assumes rate ∝ [Glc][Gly]) → phosphate ~×2 more per mole near neutral | [Z] (transfer-laden) |
| concentration law (saturation?) | **no buffer ladder at fixed reactants.** Citrate at pH 3: 0.25 M ×1.8 (Table II, invert sugar 0.5 M + Gly 0.5 M) vs 0.5 M ×3.7 (Fig. 2, glucose 0.25 + Gly 0.25 M) — different sugar and loadings, so only "no sign of saturation to 0.5 M". Equimolar Glc = Gly ladder at FIXED 0.25 M phosphate (Fig. 1): A490 rises **linearly** in c (0.066 at 0.05 M → 0.68 at 0.4 M), not as c² | [F] / [Z] |
| temperature | two temperatures (~121 °C jars, 100 °C reflux) but **in different systems** (0.25 vs 0.05 M reactants; pH 7.5 vs 8.05); the buffer factor is ×7.6 and ×8-13 respectively. Not a temperature series | [F] / [Z] |
| pH confound | addressed by the authors. In the jars the unbuffered pot falls MORE (pH change −2.4 vs −2.0, Fig. 1). In the reflux pair the PHOSPHATE pot is the more acid one for the first 8 h (−0.65 vs −0.05 at 2 h, Fig. 5) yet ×11 browner: "the difference was highly significant from the outset" [V]. So at pH 7.5-8 the effect is not a pH artefact | [F] / [Q] |
| metal confound | **not addressed**: "Analar" salts, glycine "ordinary reagent grade", lacquered tinplate caps, no chelator, no metal analysis | [M] / [!] |
| ionic strength / Na⁺ | "sodium and chloride ions in the quantities used had little effect" [V]; the control experiment and its quantities are not shown | [Q] |
| O2 | N2 headspace in the earliest jars "found to make little difference" [Q], then air | [Q] |

**Answers.**

1. **Buffer catalysis is general, not phosphate-specific** (confidence ~0.8 at pH 3, ~0.6 near neutral). At pH 3
   citrate and phosphate are equal per mole, acetate ~1/3 per mole, for both the amine-free and the glycine pot
   [M]; near neutral phosphate may be ~2× citrate per mole [Z, cross-system]. The authors' mechanism is general
   acid-base catalysis of the Amadori rearrangement and its decomposition [P], citing Rosen et al. and Reynolds [C].
2. **The effect is LARGE at pH 7.5-8 (×6-13 on browning), smaller at pH 3 with glycine (×2) and largest amine-free at
   pH 3 (×11-12)** [M/F/Z]. Nothing here fades at pH 8; see sec. 5.
3. **No saturation is visible** up to 0.5 M citrate (pH 3) [Z, weak, cross-design]; there is no single buffer ladder.
4. **pH is not the explanation at pH 7.5-8** (the reflux pair, sec. 1 row "pH confound"). **Metals are not
   excluded** [!].

## 2. Conditions

- **Jars** [M]: 40 ml solution in a tall fish-paste jar, an open tube of 2 ml 2.5 N NaOH inside to trap CO2,
  rubber-ring lacquered tinplate cap; 2 h at 15 psi steam (**"about 121°"** [V]), cooled in water. Headspace air
  (N2 only in the earliest runs).
- **Systems** (Table I) [M]: A sodium phosphate 0.25 M, glucose + glycine (A1 constant total 0.5 M; A2 equimolar
  0-0.4 M each; A3 no phosphate, 0.25 + 0.25 M), **initial pH 7.5**. B sodium citrate 0.5 M, same layout, **pH 3.0**.
  C citrate + glucose varied (0.25 M glycine), **pH 8.0**. D citrate + glycine varied (0.25 M glucose), **pH 8.0**.
  E 0.25 M total phosphate / citrate mixtures, 0.25 M glucose + 0.25 M glycine, pH 7.5 (no figure). F / G (Table II)
  0.5 M invert sugar ± 0.5 M glycine with phosphate 0.25, citrate 0.25 or acetate 0.75 M, **initial pH 3.0**.
- **Reflux** [M]: 0.05 M glucose + 0.05 M glycine ± 0.2 M sodium phosphate, initial pH 8.05, boiling under reflux
  (~100 °C) up to ~11.5 h.
- **pH**: initial label pH; final pH measured (Table II; pH-change panels); temperature of measurement not stated
  (presumably room temperature) [!].
- **Analytics** [M]: transmittance at 490 mμ (Unicam SP500, 1 cm; 100-fold dilution unless stated; Fig. 2 and
  Table II system F undiluted; Fig. 5 1 + 9), absorbance at 285 mμ (1000-fold; Fig. 2 1 + 399); CO2 by acid
  titration of the NaOH trap with a blank jar. Replicates and errors: **not stated**.
- **Invert sugar** by acid inversion of sucrose + neutralisation. Water activity: dilute aqueous.

## 3. Quantitative results

**Table II** (initial pH 3.0, 0.5 M invert sugar, jars) [M]:

| system | buffer | molar | final pH | A285 | T490 % | CO2 mg/100 ml |
|---|---|---:|---:|---:|---:|---:|
| F (no Gly) | none | — | 2.40 | 0.079 | 95.0* | 0.0 |
| F | phosphate | 0.25 | 2.58 | 0.282 | 56.3* | 0.0 |
| F | citrate | 0.25 | 2.67 | 0.386 | 55.0* | 2.8 |
| F | acetate | 0.75 | 2.40 | 0.333 | 53.9* | † |
| G (0.5 M Gly) | none | — | 3.12 | 0.497 | 87.2 | 34.7 |
| G | phosphate | 0.25 | 3.02 | 0.895 | 76.2 | 42.9 |
| G | citrate | 0.25 | 2.85 | 1.080 | 78.2 | 41.3 |
| G | acetate | 0.75 | 2.99 | 0.866 | 78.2 | † |

\* undiluted; † acetic acid volatilises. **[Z]** A490 = −log10(T): F 0.022 / 0.249 / 0.260 / 0.268 (×11.2 / 11.7 /
12.0); G (100-fold dilution) 0.060 / 0.118 / 0.107 / 0.107 (×2.0 / 1.8 / 1.8).

**Fig. 1** (0.25 M phosphate, pH 7.5, jars) [F]:
- equimolar Glc = Gly, 0.05 → 0.40 M: T490 86 / 73 / 62 / 51 / 41 / 33 / 27 / 21 % (A490 0.07 → 0.68, linear);
  A285 linear to 0.80; CO2 linear to ~205 mg/100 ml; pH change −0.4 → −2.4.
- constant total 0.5 M: maximum near Glc 0.25-0.3 / Gly 0.2-0.25 (T ≈ 40 %, A285 ≈ 0.50, CO2 ≈ 130).
  **Glucose 0.5 M + phosphate, NO glycine: T ≈ 77 %, A285 ≈ 0.30, CO2 ≈ 35 mg/100 ml, pH −1.65.**
- □ (A3, no phosphate, 0.25 + 0.25 M): T ≈ 89 %, A285 ≈ 0.07, CO2 ≈ 25, pH ≈ −2.4.

**Fig. 2** (0.5 M citrate, pH 3.0; T undiluted, A285 at 1 + 399) [F]: equimolar ladder T 85 / 73 / 60 / 41 / 20 / 6 %
at 0.05-0.30 M (super-linear in c); □ (no citrate, 0.25 + 0.25 M) T ≈ 65 %, A285 ≈ 0.09 vs ≈ 0.34 with citrate →
**×3.7-3.8** [Z]. Constant-total maximum T ≈ 12 % at Glc ~0.2 / Gly ~0.3.

**Fig. 3** (0.25 M glycine, citrate vs glucose, pH 8.0) [F]: Glc 0.5 / citrate 0: T ≈ 80 %; maximum T ≈ 55 % at Glc
~0.3-0.35 / citrate ~0.15-0.2; equimolar Glc = citrate 0.25 M: T ≈ 59 %; □ (citrate + glucose, no glycine) T ≈ 90 %,
pH change ≈ −1.95.

**Fig. 4** (0.25 M glucose, citrate vs glycine, pH 8.0) [F]: citrate 0.5 / Gly 0: T ≈ 92 %; Gly 0.5 / citrate 0:
T ≈ 71 %, pH change −3.3; maximum T ≈ 51-53 % near Gly 0.35 / citrate 0.15; equimolar arm pH change flat at −2.4,
T 89 → 40 % (0.05 → 0.4 M), A490 roughly linear in c; □ (no glucose) T 100 %, pH 0.

**Fig. 5** (reflux, 0.05 + 0.05 M, pH 8.05, T at 1 + 9) [F]:

| h | T no buffer | T 0.2 M phosphate | ΔpH none | ΔpH phosphate | A490 ratio [Z] |
|---:|---:|---:|---:|---:|---:|
| 2 | 98.5 | 85 | −0.05 | −0.65 | ×11 (×7-20 on the read error) |
| 4 | 95.5 | 69 | −0.35 | −0.95 | ×8 |
| 6 | 96 | 58 | −0.85 | −1.30 | ×13 |
| 8.4 | 94 | 50.5 | −1.45 | −1.40 | ×11 |
| 11.5 | 92.5 | 46 | −2.0 | −1.45 | ×10 |

## 4. Mechanism evidence and its kind

- **Within-study factor** (kind 6): buffer present / absent at the same start pH, several buffers at one pH, mixtures
  at constant total. End-point browning, A285 "carbonyls" and CO2 only: **no intermediate is measured, no fit**.
- The general acid-base catalysis of the Amadori rearrangement and its decomposition is **[P]** (and [C] Rosen et al.
  1957-58, Reynolds 1959). The claim "buffer salts accelerate browning by the same mechanism" rests on additivity of
  phosphate and citrate (System E, unshown) [Q].
- CO2 from glucose + phosphate with no amine (Fig. 1, ~35 mg/100 ml ≈ 8 mM [Z]) shows an amine-free decarboxylating
  route at pH 7.5 [F]; the authors doubt CO2 comes from the organic buffer anions [Q].
- Lab: Cerebos, Edinburgh — one industrial group, one paper.

## 5. Bearing on open questions

- **E5 (no buffer-species catalysis in the engine): SUPPORTS that buffer catalysis is real and large** (×6-13 at pH
  7.5-8, ×2 with glycine at pH 3, ×11-12 amine-free at pH 3), a new laboratory. **CONTRADICTS two parts of the
  proposed fix**: (a) "fading at pH 8" — browning with 0.2-0.25 M phosphate is ×6-13 at start pH 7.5-8.05; the fade
  rests on Davidek 2002's ARP decay at pH 8 (×1.27, F03), a different response, so the fade may be step-specific, not
  general; (b) a phosphate-specific factor — citrate and acetate act too (per mole: citrate ≈ phosphate at pH 3,
  acetate ~1/3). Says nothing on saturation beyond "none seen to 0.5 M citrate".
- **CF-011 (Xie phosphate pots form Cys-Amadori 3-6× faster than Jiangnan water pots): LEAVES UNCHANGED**, with mild
  support for the plausibility of "buffer catalysis" at that magnitude (0.2-0.25 M buffer → ×2-13 on a glucose-glycine
  response). It adds a testable consequence: if the explanation is general buffer catalysis, a citrate pot at the
  same pH should also run fast; the CF-011 decider could include a citrate arm. Response here is browning, not
  TTCA → Cys-Amadori.
- **F14 amine-free degradation ("phosphate REQUIRED at pH 7, 1 lab")**: SUPPORTS buffer-promoted amine-free sugar
  degradation at a second lab (×11-12 at pH 3), but shows it is **not phosphate-specific** at pH 3. At pH 7.5-8 there
  is amine-free browning with phosphate (T 77 %) and with citrate (T ≈ 90 %), but no buffer-free amine-free control.
- **Evidence map 1.2** ("sugar loss, browning: many labs"): adds one lab, one temperature pair in different systems.
  `LITERATURE_SEARCH_2026-10_VARIABLES.md` row 10 ("citrate vs phosphate at 121 °C and ~100 °C") **overstates it**:
  the 100 °C reflux has phosphate vs none only (sec. 6).

## 6. Disagreements found

1. **Citrate's activity.** Saunders: 0.25 M citrate ≈ 0.25 M phosphate at pH 3 / ~121 °C, and ×3-4.6 over no buffer at
   pH 8 [M/Z]. Bell 1997 (`bell1997_extraction.md`): citrate 0.02-0.5 M at pH 7.0 / 25 °C indistinguishable from zero.
   Davies 2005 (`davies2005_extraction.md`): 0.5 M citrate gives no amine-free browning at pH 5.5 / 50 °C. Rizzi 2004:
   citrate ~1/2 of phosphate at pH 7.3 / ~100 °C. Candidate causes: a high activation energy for carboxylate catalysis
   (activity appears only ≥ 100 °C), pH (pH 3 and 8 vs 5.5-7), response (amine-free browning vs Gly loss).
   Bell's own dossier already flags the Rizzi tension (sec. "Contradiction with Rizzi 2004"); Saunders makes it two
   high-temperature labs against one 25 °C lab.
2. **"Fades at pH 8"** (F03 atlas summary; E5 fix): Saunders ×6-13 on browning at pH 7.5-8.05 vs Davidek ×1.27 on
   DFG decay at pH 8 (pH-stat, 90 °C, 0.1 M phosphate). Causes: response (browning integrates all steps after the
   Amadori product; DFG decay is one step), buffer concentration (0.2-0.25 vs 0.1 M), pH control (drift vs pH-stat).
3. **`LITERATURE_SEARCH_2026-10_VARIABLES.md` row 10** says this paper compares citrate vs phosphate at two
   temperatures; it does not (citrate is only in the jars; the reflux is phosphate vs none, in another system).
4. **[!] Figure legends.** The note under Fig. 3 says □ "indicates value where phosphate was present", while the Fig. 1
   caption says phosphate is present "except where otherwise indicated". Every □ matches the third row of its
   system in Table I (A3 no phosphate; B3 no citrate; C3 no glycine; D3 no glucose — the last is T 100 %, "no reaction
   was of course possible in the absence of sugar"), so I read □ as the omission control. The printed note is wrong
   or garbled.

## 7. Ledger-ready measurements

1. **Table II, all eight rows** (printed): buffer identity at pH 3.0, ± glycine; the cleanest "is phosphate special"
   evidence. Key as paired arms (buffer vs none) with A490 derived from T (state the dilutions).
2. **Fig. 5 reflux pair** (two readers): phosphate vs none at 100 °C with a pH trajectory — the only arm here with
   time resolution and with the buffered pot more acid than the control.
3. **Fig. 1 □ vs the equimolar line at 0.25 M** (two readers): phosphate vs none at pH 7.5, ~121 °C, with A285 and
   CO2.
4. Fig. 2 □ vs line at 0.25 M (citrate 0.5 M vs none, pH 3): a concentration point beside Table II's 0.25 M citrate
   (different sugar; key with that caveat).

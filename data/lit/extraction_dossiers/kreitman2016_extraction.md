# Kreitman, Danilewicz, Jeffery & Elias 2016 (Part 1) — EXTRACTION (Cu(II)-mediated oxidation of H2S, cysteine and hexanethiols in model wine)

**Source on disk:** `data/articles/kreitman2016.pdf`, the ACS "Just Accepted" manuscript (44 PDF pages: a cover
sheet followed by manuscript pages 1–43; the printed manuscript page number is the PDF page minus one; figures on
manuscript pp. 29–42). This is not the typeset version, so journal page numbers do not appear and are not
quoted. The text layer was extracted with `pdftotext -layout`, and the results pages (ms pp. 10–12) and Figures
3, 6 and 8 (ms pp. 31, 36, 38) were read and checked by eye from the page images on 2026-10-09. Written for the
missing-metal-catalysis question in `docs/guides/EXPERIMENTS.md` sec. 1. The companion dossiers are
`ehrenberg1989_extraction.md` (the rate law) and `kreitman2016b_extraction.md` (Part 2, Fe + Cu).

| field | value |
|---|---|
| Title | "Reaction Mechanisms of Metals with Hydrogen Sulfide and Thiols in Model Wine. Part 1: Copper Catalyzed Oxidation." |
| Authors | G. Y. Kreitman, J. C. Danilewicz, D. W. Jeffery, R. J. Elias (Penn State; Canterbury, Kent; Adelaide) |
| Venue | J. Agric. Food Chem. 2016; on the cover sheet "Just Accepted Manuscript", published on the web 30 Apr 2016 (volume and pages 64:4095 as supplied, not printed in this PDF) |
| DOI | 10.1021/acs.jafc.6b00641 (printed) |

## 1. Methods (ms pp. 5–10)

- **Model wine**: tartaric acid 5 g/L, ethanol 12 % v/v, adjusted to **pH 3.6** with NaOH. No added iron (iron is
  Part 2).
- **Thiols, each 300 µM** in air-saturated model wine (1 L): H2S (from NaSH), L-cysteine (Cys),
  6-sulfanylhexan-1-ol (6SH, a primary thiol), 3-sulfanylhexan-1-ol (3SH, a secondary thiol). **Cu(II) as
  CuSO4: 50 µM** (H2S, Cys, 6SH) or **100 µM** (3SH). Mixed system: the Methods give H2S 100 µM + Cys 400 µM +
  Cu(II) 100 µM (ms p. 5) while the Results give H2S 75 µM + Cys 468 µM (ms p. 12); Fig. 6's caption says
  "~100 µM" and "~400 µM". The three disagree as printed.
- Transferred at once to 60 mL glass BOD bottles, overfilled and stoppered (**no headspace**), dark, **"ambient
  temperature"** (no number printed). Triplicates, one sacrificial bottle per time point.
- Starting O2 "approximately 7 mg/L (~220 µM)" (ms p. 6), by PSt3 oxidots and a NomaSense meter.
- Thiols by Ellman's reagent (DTNB); H2S and Cys in the mixed system by monobromobimane derivatisation and
  HPLC-MS/MS; Cu(I) by bathocuproinedisulfonate (BCDA); total Cu after 0.45 µm filtration by ICP-OES; EPR of
  Cu(II) (0.5 mM Cu + 1.5 mM thiol, 100 K); acetaldehyde (AC) as its DNPH derivative.
- Anaerobic disulfide test: 6SH 600 µM in 3 mL, argon-sparged, Cu(II) 50/100/200 µM, 5 min.
- No chelator treatment arm. DTPA appears only in the bimane derivatisation buffer at pH 9.5 (analytical), and
  BCDA only to dissolve the Cu(I) complex.

## 2. Findings that matter

### 2a. Fast stoichiometric step: Cu(II) + 2 RSH (ms pp. 10–12, Fig. 3, Fig. 4A)

| system | Cu(II) / µM | immediate thiol uptake | rest consumed |
|---|---|---|---|
| H2S 300 µM | 50 | ~1.4 mol equiv (72 µM) | "fully consumed within 72 h" |
| Cys 300 µM | 50 | **101 µM** (~2 equiv) | "fully consumed within 48 h" |
| 6SH 300 µM | 50 | 121 µM (~2 equiv) | within 48 h |
| 3SH 300 µM | 100 | 2 equiv (210 µM) after 2 hours | "not fully reacted after 168 h" |
| H2S 75 + Cys 468 µM (Results) | 100 | 53 µM H2S + 135 µM Cys within 5 min = 189 µM, "~2:1" | H2S gone within 40 min; Cys after 48 h |

- EPR: Cu(II) "immediately reduced to Cu(I)" by Cys, 6SH and H2S; 3SH complete after 2 h (ms p. 11).
- **Thiol:Cu stoichiometry**: abstract, "~1.4:1 H2S:Cu and ~2:1 thiol:Cu complexes". Under argon, Cu(II) 50 / 100 /
  200 µM gave 19.7 ± 3.6 / 43.4 ± 3.1 / 98.2 ± 3.6 µM 6SH disulfide, i.e. **0.5 mol disulfide per mol Cu(II)**
  in 5 min (ms pp. 12, 16): one thiol oxidised, one left bound to Cu(I). The dried 6SH–Cu(I) aggregate dissolved
  in BCDA released 1.17 ± 0.02 mM Cu(I) and 1.17 ± 0.13 mM 6SH, **~1:1 Cu(I):thiolate in the aggregate** (ms p. 12).
- Proposed mechanism (ms p. 13, Fig. 10): Cu(II)(SR)2 → Cu(I) intermediate; two associate and form RSSR bound to
  Cu(I), without free thiyl radicals (4-MeC and DMPO traps did not change disulfide formation, ms p. 12); the
  Cu(I)–SR aggregates (fine white/yellow precipitate, removable at 0.45 µm from 5 to 45 min for 6SH).
- Fig. 3, Cys trace, read from graph, approx.: 300 µM at t = 0 → ~197 at ~5 min → ~190 at 0.5 h → ~169 at 2 h →
  ~137 at ~10 h → ~55 at ~24 h → ~0 at ~48 h.
- Fig. 6, Cys in the mixed system, read from graph, approx.: ~467 → ~332 at ~5 min → ~311 at ~0.7 h → ~274 at 2 h
  → ~167 at ~10 h → ~10 at ~25 h → 0 at ~48 h. H2S ~75 → ~20 at ~5 min → ~0 by ~0.7 h.

### 2b. Slow catalytic phase: O2 stoichiometry and products (ms pp. 12–13, 18–20, Figs 7–9)

| system | O2 consumed / µM | O2 : thiol (as printed) | AC / µM | O2 : AC (as printed) |
|---|---|---|---|---|
| H2S (284 µM reacted) | 175 ± 9 | ~1:1.6 | 79 ± 2 | 2.2:1 |
| Cys (299 µM reacted) | 66 ± 6 | **~1:4.5** | 26 ± 0.3 | 2.5:1 |
| 6SH | 76 ± 6 | — | 52 ± 4 | 1.5:1 |
| 6SH 240 µM + Cu 50 µM, 262 h | 69 ± 8.0 | ~1:3.3 (231 ± 2.5 µM reacted, 116 ± 2.7 µM disulfide) | — | — |
| 3SH (~74 µM reacted beyond the Cu-bound 100) | 23 ± 1 (28 in the ratio calc.) | 1:2.6 | 13 ± 0.8 | 1.8:1 |
| Cys + H2S | 117 ± 5.2 | — | 54 ± 3 | 2.1:1 |

- "Minimal O2 uptake (<5 µM in all treatments) ... within the first 30 min" (ms p. 13): the fast step is
  anaerobic, a Cu(II) reduction.
- **Disulfide was essentially the sole product** for 6SH (1:0.5 RSH:RSSR, ms p. 17). Expected limits: 1:4 O2:thiol
  if H2O2 is reduced two-electron by the Cu(I) complex (Figs 13–14); 1:3 if H2O2 goes through the Fenton route and
  oxidises ethanol to AC (Fig. 15). Cys at 1:4.5 sits at or beyond the 1:4 end.
- No hydroperoxyl radicals (4-MeC not consumed, no catechol–thiol adducts, ms p. 18). O2 is proposed to be
  reduced two-electron to H2O2 by adjacent Cu(I) ions in the aggregate (Fig. 13).
- Cu partition: H2S, ca. 60 % of Cu filterable within 5 min and up to 24 h, ~90 % after 72 h (green-black
  precipitate); 6SH, essentially all Cu(I) complex retained at 0.45 µm from 5 to 45 min, with Cu partly released
  later (ms p. 11, Fig. 5).

### 2c. Wine-context numbers (ms pp. 3–4, Introduction)

Cys + N-acetylcysteine + homocysteine ca. 20 µM in white wines; GSH ca. 40 µM (Sauvignon blanc); Cu fining dose
3–6 µM; H2S ca. 300 nM.

## 3. What it means for the model

- **Rate constants: none.** The paper prints stoichiometries and time courses only, at one unstated ambient
  temperature, one pH (3.6), with no headspace and with Cu at 50–100 µM, 200–400× the 0.25 µM of
  EXPERIMENTS.md.
- **What it adds to Ehrenberg 1989**: the Cu(II) → Cu(I)–thiolate step is **instant and stoichiometric even at
  pH 3.6** (2 RSH per Cu(II), half to disulfide), in the presence of a large tartrate excess. Low pH does not
  stop the metal from complexing the thiol. What is slow at pH 3.6 is the **catalytic turnover**, the
  reoxidation of the Cu(I)–SR aggregate by O2. For the EXPERIMENTS.md pot the stoichiometric step is
  negligible: 0.25 µM Cu(II) removes 0.5 µM of 33 mM cysteine (derived here). Only turnover matters.
- **Turnover estimate** (derived here, from the Fig. 3 graph reading, approx.): Cys falls from ~169 µM at 2 h to
  ~0 at ~48 h, ≈ 3.7 µM/h (169/46), or ≈ 5.9 µM/h over 10–24 h ((137 − 55)/14). With 50 µM Cu that is ≈ 1–2×10⁻³
  thiols per Cu per minute at pH 3.6 and ambient temperature. Ehrenberg's 37 °C, pH 7.2 Cu turnover at 30 mM
  cysteine is 240 µM min⁻¹ per µM Cu (Table 3). The ratio, ~10⁵, mixes pH, temperature, tartrate, ethanol,
  aggregation at 50 µM Cu and a closed O2 supply, so it is an order-of-magnitude bracket, not a pH law. It does
  say that **a pH-5 Cu channel cannot be assumed to run at its pH-7 rate**.
- **Stoichiometry for any implementation**: O2:Cys ≈ 1:4.5 here (Cu only), against 1:3.2 with Fe alone and 1:2.6
  with Fe + Cu (Part 2). The product is cystine (disulfide), with little H2S-type or sulfur-oxyanion chemistry.
  This matches Ehrenberg's eqn (3) (4 RSH per O2) better than eqn (1) (2 per O2).
- **Live engine values** (see `ehrenberg1989_extraction.md` sec. 3 for keys and arithmetic): `k_cys_thermal`
  log10 k(145 °C) = −2.066359667088082 /min with Ea 55.1 kJ/mol; `k_cys_h2s` from Zheng & Ho; `k_cys_ox` = 0
  (inert, b9 has no oxygen block). None has a metal term. This paper confirms the missing channel's **kind**
  (Cu-dependent, O2-consuming, cystine-forming, thiolate-complex-mediated) but gives no number to set it with.
- **H2S**: Cu(II) binds and oxidises H2S faster than Cys in a mixture (53 of 75 µM H2S vs 135 of 468 µM Cys in
  5 min; H2S gone within 40 min). If a metal channel is added for Cys, the sulfide pool should not be exempt.

## What it does not give

- No rate constants, half-lives fitted by the authors, reaction orders or temperature dependence. "Ambient
  temperature" has no number.
- One pH (3.6) in a tartrate/ethanol matrix; no phosphate, no pH series.
- No iron (Part 2), no chelator arm, no Cu concentration series in air (the Cu series is anaerobic, 6SH only).
- Cu at 50–100 µM, so Cu(I)–thiolate aggregation and precipitation, which govern the slow phase here, may not
  occur at 0.25 µM.

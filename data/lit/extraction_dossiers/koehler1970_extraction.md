# Koehler & Odell 1970 — EXTRACTION
DOI: 10.1021/jf60171a041 (printed: none; Crossref 2026-11-01 resolves to "Factors affecting the formation of pyrazine compounds in sugar-amine reactions", J. Agric. Food Chem. 18, 895-898, 1970, authors Koehler / Odell)
### *"Factors Affecting the Formation of Pyrazine Compounds in Sugar-Amine Reactions"*
### Glucose + asparagine 0.1 mol each in 200 mL diethylene glycol + 20 mL water (~0.45 M each), 120 °C, 24 h (the "normal" pot), varied one factor at a time: temperature 80-150 °C (Fig. 1), time to 72 h (Fig. 2), Asn : glucose ratio (Fig. 3), H2SO4 / NaOH equivalents (Fig. 4); N source (Table I), sugar (Table II), fed C2 / C3 / C4 carbonyls (Table III). **MPZ and "dimethylpyrazine" only in most pots; PZ is reported only in Table III. So R is computable only for five fed-carbonyl / glucosamine pots.**

**Read 2026-11-01 from `data/articles/koehler1970.pdf`** (4 pp., JAFC 18(5):895-898; Acrobat Capture scan with an OCR text
layer). Prose and Tables I-III from `pdftotext -layout`, checked against a 170 dpi render; **Figs. 1-4 read off the
170 dpi render** (read error ~±15 µmol on the 0-500 axes, ±50 µmol on Fig. 4's 0-2000 axis, ±0.2 on the DMP/MP ratio
axis). Every derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Philip E. Koehler, George V. Odell (Biochemistry, Oklahoma State University; Koehler later Univ. Georgia) |
| Venue | J. Agric. Food Chem. 18 (5), 895-898, 1970; received 25 Mar 1970 |
| Laboratory | Oklahoma State (same as koehler1969). **Not in `lab_ids.yml`** |

## 1. What this paper answers for the model

| need | what Koehler 1970 gives | tag |
|---|---|---|
| PZ / MPZ / DMP per pot | **Table III only** (µmol / 220 mL): glycerol-Asn 2 / 7 / 2; acetaldehyde-Asn 29 / 51 / 291; glyoxal-Asn 527 / 91 / ND; glyoxal-NH4OH 253 / 16 / ND; glucosamine (neutralised) 107 / 252 / 27. Tables I-II and Figs. 1-4: MPZ and DMP only | [M] |
| R | glucosamine **22**; glycerol-Asn 12 (2-7 µmol, near LOQ); **acetaldehyde-Asn 0.31**; glyoxal-Asn: DMP ND, so R >= 157 if the LOD is 0.1 µmol (>= 16 at 1 µmol); glyoxal-NH4OH >= 10 (0.1 µmol LOD) | [Z] |
| response factors | FID, "comparison of peak areas of samples with those of standards" [V], three injections each; recovery of PZ, MPZ and 2,5-DMP (2.5 mmol each) through one falling-film pass: 98 % [M]. LOD not printed; "less than 0.1 µmol" [V] is used as "essentially none" below 100 °C | [M] / [!] |
| 2,5- vs 2,6-DMP | "dimethylpyrazine" on packed Carbowax 20M at 150 °C isothermal; identity of some peaks by IR. 2,5 / 2,6 not reported separately [!] | [!] |
| pH dependence | **Fig. 4** (acid / base equivalents in glycol-water, pH not measurable [Q]): ~0.1 mol H2SO4 equivalents -> "practically zero"; ~0.1 mol NaOH -> MPZ x10, DMP x5 [V/M in text]; read: MPZ ~120 -> ~1100, DMP ~400 -> ~1950 µmol [F] | [M] / [F] |
| temperature dependence | **Fig. 1** (24 h end points, 100-140 °C): "essentially no" pyrazines below 100 °C; DMP x10 and MPZ x4 from 100 to 140 °C [M, text]; DMP/MP ratio peaks ~4.5 at ~110-115 °C [M], ~2.8 at 105 and 135 °C [F] | [M] / [F] |
| nitrogen source | **Table I** (glucose + Asn / Gly / Ala / Lys / Asp / Asp + NaOH / NH4OH / NH4Cl, 120 °C 24 h) | [M] |
| separates P1x / P1y from P2b? | **No.** The glyoxal-Asn pot (MPZ without detectable DMP) bounds R from below only if the LOD is known; it is not | [Z] |

**Answers.**

1. **MPZ forms from glyoxal (a pure C2 source) with asparagine: PZ 527, MPZ 91 µmol, DMP ND** [M]. Under pairing,
   MPZ needs an aminoacetone-type C3 partner, which would also give DMP: at R = 4 the predicted DMP is 3.9 µmol, ~40x the
   stated sub-0.1 µmol "essentially none" level [Z]. So either R >= ~16-160 (P1x / P1y with a large selectivity), or MPZ
   here does not come from A x M pairing (a P2b-type route: e.g. a C1 addition to the pyrazine dihydropyrazine, or a C3
   unit built from glyoxal + Asn carbon [P]). The authors read it as "some 3-carbon units must be formed" [Q]. Not
   decisive (no LOD; no label), but it is the cleanest "MPZ without its DMP sibling" row on disk.
2. **Acetaldehyde-Asn gives R 0.31** (DMP 291 µmol dominant) [Z]: below every pairing rival's floor 4h with h >= 0.1.
   It needs extra DMP routes or strong pool segregation; another low-R row beside CF-087.
3. **Glucose + NH4OH: MPZ 1415, DMP 56 µmol (DMP/MP 0.04)** [M], the same MPZ-dominated pattern as UC Davis's glucose +
   NH3 (shibamoto1976 / 1977). **NH4Cl: ND** [M] (disagrees with koehler1969's 59 µmol).
4. **Base raises MPZ more than DMP** (x10 vs x5) [M]: the pyrazine distribution is pH-dependent in this system.
   Without PZ, R's pH trend cannot be computed.

## 2. Conditions [M]

| item | value |
|---|---|
| "Normal" pot | D-glucose 0.1 mol + L-asparagine 0.1 mol in 200 mL diethylene glycol + 20 mL water [~0.45 M each, Z], 120 ± 1 °C, 24 h, reflux |
| Variations | 80-150 °C; 0-72 h; Asn : glucose 1:3, 1:1, 3:1; 0-0.08+ equivalents H2SO4 or NaOH; other N sources and sugars at 0.1 mol; Table III: 0.1 mol each of a carbonyl and Asn or NH4OH |
| pH | not measurable in the glycol medium [Q]; acid / base given as equivalents |
| Isolation | two falling-film passes (steam jacket), liquid-N2 traps, 5 x 5 mL CH2Cl2, concentrated to 5 mL |
| GC | PE 801 FID; 20 ft x 1/4 in glass, 15 % Carbowax 20M, 150 °C isothermal; external standards same day; 3 injections each |
| Optimised run | glycol without water, 0.1 mol glucose + Asn + NaOH, 130 °C, 48 h: 190 mg DMP (1.8 %), 75 mg MPZ (0.8 %) [M] |

## 3. Quantitative results [M]

**Table I** (µmol / 220 mL; 0.1 mol glucose + N source, 120 °C, 24 h):

| N source | MPZ | DMP | DMP/MP |
|---|---:|---:|---:|
| Asn | 110 | 402 | 3.7 |
| Gly | ND | 14 | — |
| Ala | 23 | 82 | 3.6 |
| Lys | 11 | 20 | 1.8 |
| Asp | 16 | 81 | 5.1 |
| Asp + NaOH | 92 | 282 | 3.1 |
| NH4OH | 1415 | 56 | 0.04 |
| NH4Cl | ND | ND | — |

**Table II** (0.1 mol Asn + sugar): glucose 110 / 402 (3.7); fructose 137 / 1142 (8.4); sucrose 85 / 166 (2.0);
arabinose 94 / 48 (0.5).

**Table III** (0.1 mol each, µmol: PZ / MPZ / DMP / tetramethyl): propionaldehyde-Asn ND / ND / ND / ND; glycerol-Asn
2 / 7 / 2 / ND; acetaldehyde-Asn 29 / 51 / 291 / ND; acetaldehyde-NH4OH 399 / ND / ND / ND; glyoxal-Asn 527 / 91 / ND / ND;
glyoxal-NH4OH 253 / 16 / ND / ND; 2,3-butanedione-Asn ND / ND / ND / 2662; 2,3-butanedione-NH4OH ND / ND / ND / 640;
hydroxyacetone-Asn ND / ND / 9725 / ND; hydroxyacetone-NH4OH ND / ND / 2208 / ND; glucosamine 107 / 252 / 27 / ND.

**Figures [F]** (normal pot unless stated):
- Fig. 1 (24 h): MPZ ~5 (100 °C), ~30 (107), ~65 (112), ~100 (117), ~205 µmol (137 °C); DMP ~0 (100), ~120 (110),
  ~290 (115), ~400 (120), ~550 µmol (137 °C). The printed "ten-fold" (DMP) and "four-fold" (MPZ) 100 -> 140 °C are
  text values [M].
- Fig. 2 (120 °C): MPZ ~120 µmol at 24 h, ~150 at 72 h; DMP ~400 at ~24 h, plateau; DMP/MP ratio rises to ~3 by ~9 h
  and stays [M, text]; MPZ is the larger product only before ~3 h [M, text].
- Fig. 3: 1:3 Asn : glucose -> both near zero on the plot (text: MPZ x1/10, i.e. ~11 µmol; DMP x1/125, ~3 µmol [Z]);
  3:1 -> MPZ ~85, DMP ~345 µmol (text: MPZ -25 %, DMP "very little effect").
- Fig. 4: acid 0.04-0.08 eq -> both < ~30 µmol; base 0.08+ eq -> MPZ ~1100, DMP ~1950 µmol.

## 4. Mechanism evidence and its kind

- **Fed intermediates (kind 2), one lab:** glyoxal -> PZ (+ some MPZ); hydroxyacetone -> only DMP; 2,3-butanedione ->
  only tetramethylpyrazine; acetaldehyde + NH4OH -> only PZ, + Asn -> PZ, MPZ, DMP. Products from the fed carbonyls
  match the C2 / C3 / C4 dicarbonyl-or-hydroxycarbonyl + N -> dihydropyrazine -> pyrazine scheme [Q]; no labels, so the
  carbon of the "extra" products (MPZ from glyoxal; DMP from acetaldehyde + Asn) is not traced.
- **Within-study factors (kind 6):** T, time, reactant ratio, acid / base, N source, sugar.
- The glucosylamine -> tetrahydroxybutylpyrazine route at low T is **proposed** [P].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 / CF-077 | CONTRADICTS P0 (glucosamine R 22) | one calibrated pot |
| P1 (R <= 4h) | CONTRADICTS (glucosamine 22; glyoxal-Asn >= 16 at any LOD <= 1 µmol) | unresolved 2,5 / 2,6 lowers R, so conservative |
| P1x / P1y vs P2b | **LEAVES UNSEPARATED**; the glyoxal-Asn row is a lead for P2b | MPZ without detectable DMP from a C2 feed: P1x / P1y need R >= 16-160 there; P2b needs no C3 pool. Label test: [13C2]glyoxal + Asn, read MPZ M+2 vs M+3 / M+4 |
| low-R rows (CF-087) | bears on it | acetaldehyde-Asn R 0.31: forbidden by every pairing rival unless h < 0.08 or an extra DMP route exists (DMP 291 µmol from a C2 aldehyde implies one: aldol to a C3 / C4 unit [Q]) |
| E4 (pyrazine pH law) | SUPPORTS direction only | acid -> ~0, base -> x5-10, glycol medium, no pH value; supply and step not separated |
| E14 (aminoketone condensation, extra routes) | SUPPORTS "extra routes exist" (direction) | fed C2 gives MPZ; fed C2 aldehyde gives DMP |
| evidence map 2.7 "MPZ extra route: NONE" | LEAVES UNCHANGED (kind 2 lead, not kind 1) | no label |

## 6. Disagreements found

1. **Within lab vs koehler1969: Asn-glucose DMP/MP 3.7 here vs ~0.5 there**, same T, time, solvent; 10x concentration
   here. Candidate causes: concentration (the A : M supply ratio and any second-order step), two vs one falling-film
   pass, the packed-column DMP peak.
2. **Within lab vs koehler1969: NH4Cl-glucose ND here vs 59 µmol (mostly PZ) there.**
3. **Glycine: MPZ ND, DMP 14 µmol here** vs wong1988 glycine-glucose-NaOH (MPZ 0.78, DMP25 6.9, TrMP 28.6 % of identified
   compounds) and bemisyoung1993 (glycine pH 4-12, MPZ 0.2-2.2 %). Directionally the same (glycine is a poor MPZ
   source, DMP and trialkyl dominate); method (no base here, glycol medium) differs.
4. **Temperature optimum:** pyrazine yield rises to 150 °C here (above that "highly variable") vs a maximum at 120 °C in
   shibamoto1976 (2 h, NH3). Different N source, time (24 vs 2 h) and closure.

## 7. Ledger-ready measurements

1. Table III, glyoxal-Asn and glyoxal-NH4OH (PZ, MPZ, DMP ND) — MPZ without DMP from a C2 feed; key with
   `dmp: below_unstated_lod`.
2. Table III, acetaldehyde-Asn and glucosamine — the only pots with all three species (R 0.31, 22).
3. Table I (N source) and Table II (sugar) MPZ / DMP µmol — within-study ratios.
4. Fig. 4 (acid / base, text factors x10 / x5) — a direction for E4 in a low-water medium; figure-only beyond the text.

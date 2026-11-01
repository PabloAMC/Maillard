# Shibamoto & Bernhard 1976 — EXTRACTION
DOI: 10.1021/jf60206a045 (printed: none; Crossref 2026-11-01 resolves to "Effect of time, temperature, and reactant ratio on pyrazine formation in model systems", J. Agric. Food Chem. 24, 847-852, 1976, authors Shibamoto / Bernhard)
### *"Effect of Time, Temperature, and Reactant Ratio on Pyrazine Formation in Model Systems"*
### D-glucose 1 M (18 g) + ammonium hydroxide 0.1-15 M in water (100 mL), open reflux or N2-flushed sealed ampoules; 44 experiments: NH3 ladder at 100 °C / 2 h (1-20), temperature -5 to 160 °C (21-29, 40-41), time 15 min to 48 h (30-39), ± NaOH (42-44). Table II: weight % of PZ, MPZ, 2,5-DMP, 2,6-DMP and six more pyrazines per experiment, mg, and total yield % of glucose. **No amino acid, no pH measured, no labels. The method paper behind shibamoto1977.**

**Read 2026-11-01 from `data/articles/shibamoto1976.pdf`** (6 pp.; JAFC 24(4):847-852; the last page also carries the start
of Galetto & Hoffman on alkylthiophenes). Acrobat Capture scan with OCR. Prose and Table I from `pdftotext -layout`;
**Table II (printed rotated, no usable text layer) read cell by cell from a 300 dpi render rotated 90°**, columns 1-4
(PZ, MPZ, 2,5-DMP, 2,6-DMP, % and mg) and the yield column. Every derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Takayuki Shibamoto, Richard A. Bernhard (Food Science and Technology, UC Davis) |
| Venue | J. Agric. Food Chem. 24 (4), 847-852, 1976; received 20 Oct 1975, accepted 27 Feb 1976 |
| Laboratory | UC Davis: `ucdavis_shibamoto` (same lab as shibamoto1977, wong1988, bemisyoung1993) |
| Relation to shibamoto1977 | 1977's ref 3 ("J. Agric. Food Chem., in press") is the method source; this paper's method matches it (CH2Cl2 x 4, packed Carbowax FID, PZ / MPZ recovery > 98 %), so it is very probably ref 3 [Z]. Different experiments: the 8 M / 100 °C / 2 h rows here (exps 17, 33) are not 1977's control rows |

## 1. What this paper answers for the model

| need | what Shibamoto 1976 gives | tag |
|---|---|---|
| PZ / MPZ / 2,5-DMP / 2,6-DMP per pot | **yes, 44 experiments**, weight % of total pyrazines + mg | [M] |
| R (molar = weight x 0.978) | **8 M NH3, 100 °C, 2 h: 357** (exps 17 / 33, 8 replicates); NH3 ladder 0.2-15 M: 313-1200; 60-160 °C: 357-3209; 30 min-48 h: 357-485; + NaOH: 432-534. **25 °C / -5 °C, 30 days: 0.73 / 0.31.** 50 °C: 73. Median over rows >= 60 °C and NH3 >= 0.2 M: 484 | [Z] |
| response factors | FID relative RFs "using the methods described by Ettre and Zlatkis" [C]; PZ and MPZ recovered > 98 % through the extraction [M]; DMP recovery not tested [!]; emulsions above 4 h (open) lower recovery [Q] | [M] / [!] |
| 2,5- vs 2,6-DMP | reported separately (peaks 3, 4) on a 6 m packed Carbowax 20M + Igepal column; identities by MS on a 100 m glass capillary. 2,6/2,5 = 0.6-4.9 [Z] | [M] |
| pH dependence | none measured; NH3 concentration and 0.4 g NaOH act as base [Q]; NH3 escapes on reflux (8 M at 25 °C, <= 4.3 M soluble at 100 °C) [M] | [M] / [Q] |
| temperature dependence | yield: < 0.05 % below 70 °C, maximum 0.945 % at 120 °C (closed), lower at 140 / 160 °C [M]. R: ~0.3-0.7 at -5 / 25 °C (30 d), 73 at 50 °C, 360-3200 at 60-160 °C [Z] | [M] / [Z] |
| nitrogen source | NH3 only | — |
| separates P1x / P1y from P2b? | **No** (no labels; the supply ladders change basicity too) | — |

**Answers.**

1. **Agrees with shibamoto1977's R ~300-600** [Z]: 8 M / 100 °C / 2 h open gives R 357 here vs 351-483 for 1977's four
   controls (437 for the open 2 h control); NaOH pots 432-534 here vs 388-534 there. Two papers, one lab, ~0.1 dex.
2. **R is a strong function of temperature: ~0.3-0.7 at -5 / 25 °C (30 days) and ~360-500 at 100 °C** [Z]. Read as an
   Arrhenius contrast of R between 25 and 100 °C this is an apparent ΔEa of ~78 kJ/mol [Z]. Caveats: the cold pots
   yielded 0.005-0.007 % of glucose and the authors say these data "should be considered with caution" [V]; the times
   differ (30 d vs 2 h). Under P1x a T-dependent u (pair-specific Ea) is allowed; under P2b an extra MPZ route with a
   higher Ea is allowed. The contrast does not separate them, but it bounds either: the extra-route share or the
   selectivity must rise ~3 dex between 25 and 100 °C.
3. **R vs NH3 concentration (100 °C, 2 h): ~550-1200 at 0.2-4 M, ~310-440 at 5-15 M** [Z]. Under any pairing rival
   with fixed constants R does not depend on how much aminocarbonyl is made; here R falls ~0.3 dex as NH3 rises, with
   basicity changing at the same time [Q]. Not decisive.
4. **Closed systems were N2-flushed**: ampoules cooled 30 min in ice, N2 bubbled 5 min, flame sealed [M]. Open systems
   refluxed in air.

## 2. Conditions [M]

| item | value |
|---|---|
| Reactants | 18 g D-glucose (1 M) + ammonium hydroxide to 0.1-15 M, water to 100 mL |
| Vessels | open: 300 mL round-bottomed flask + Friedrich condenser open to air, oil bath; closed: 100 mL Kjeldahl flask as heavy-wall ampoule, N2 bubbled 5 min, flame sealed, oven (used above 100 °C) |
| Temperature | solution temperature; 100 °C standard; -5 to 160 °C series |
| Time | 2 h standard; 15 min to 48 h; 30 days at 25 / -5 °C |
| Base | 0.4 g NaOH (0.1 M) in exps 42-44 |
| Work-up | 4 x 50 mL CH2Cl2; Na2SO4 12 h; rotary flash evaporation; residue weighed; GC |
| GC | Varian 1200 FID, 6 m x 2.5 mm packed 5 % Carbowax 20M + 0.25 % Igepal CO-880, 70-150 °C at 0.5 °C/min; Infotronics integrator. ID: HP 5700A with 110 m glass capillary Carbowax 20M; Finnigan 3200 GC-MS |
| Excluded | seven minor pyrazines (< 0.1 % of total area) ignored for quantitation; peaks 12-16 (< 1 %) excluded [M] |
| Replicates | 1-8 per experiment (Table I); no SD |

## 3. Quantitative results (Table II, weight % of total pyrazines; yield = weight of pyrazines / 18 g glucose x 100) [M]

R [Z] = molar MPZ² / (PZ x 2,5-DMP). Rows in the order of Table I.

| exp | condition | PZ | MPZ | 2,5-DMP | 2,6-DMP | yield % | R |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | NH3 0.1 M, 100 °C, 2 h | 6.28 | 84.60 | 1.69 | 6.97 | 0.002 | 660 |
| 2 | 0.2 M | 4.45 | 86.18 | 1.36 | 6.62 | 0.091 | 1200 |
| 3 | 0.3 M | 4.66 | 84.16 | 1.78 | 6.86 | 0.212 | 835 |
| 4 | 0.4 M | 4.42 | 82.27 | 2.11 | 7.85 | 0.315 | 710 |
| 5 | 0.5 M | 4.40 | 81.74 | 2.17 | 8.26 | 0.389 | 684 |
| 6 | 0.6 M | 4.20 | 82.31 | 2.26 | 8.38 | 0.417 | 698 |
| 7 | 0.7 M | 3.55 | 81.99 | 2.99 | 8.26 | 0.467 | 619 |
| 8 | 0.8 M | 3.02 | 80.93 | 3.34 | 9.02 | 0.484 | 635 |
| 9 | 0.9 M | 2.10 | 76.73 | 5.02 | 11.33 | 0.483 | 546 |
| 10 | 1.0 M | 1.69 | 76.32 | 5.46 | 11.59 | 0.501 | 617 |
| 11 | 2.0 M | 1.50 | 74.32 | 5.47 | 12.84 | 0.527 | 658 |
| 12 | 3.0 M | 1.80 | 74.09 | 5.41 | 12.58 | 0.684 | 551 |
| 13 | 4.0 M | 2.15 | 76.51 | 5.20 | 11.21 | 0.763 | 512 |
| 14 | 5.0 M | 3.61 | 76.95 | 4.98 | 10.87 | 0.833 | 322 |
| 15 | 6.0 M | 3.65 | 78.08 | 4.81 | 9.97 | 0.887 | 340 |
| 16 | 7.0 M | 3.71 | 79.33 | 4.27 | 9.18 | 0.918 | 389 |
| 17 | 8.0 M (8 reps) | 3.81 | 81.65 | 4.79 | 8.99 | 0.945 | 357 |
| 18 | 9.0 M | 3.86 | 80.43 | 3.75 | 8.14 | 0.940 | 437 |
| 19 | 10 M | 4.00 | 80.22 | 4.25 | 7.91 | 0.934 | 370 |
| 20 | 15 M | 4.10 | 79.49 | 4.82 | 7.89 | 0.939 | 313 |
| 21 | 50 °C, 2 h | 2.39 | 62.28 | 21.86 | 1.35 | 0.001 | 73 |
| 22 | 60 °C | 5.10 | 87.16 | 2.90 | 3.60 | 0.027 | 502 |
| 23 | 70 °C | 0.64 | 92.34 | 4.06 | 2.38 | 0.044 | 3209 |
| 24 | 80 °C | 1.96 | 90.91 | 2.60 | 3.48 | 0.155 | 1586 |
| 25 | 90 °C | 3.44 | 82.31 | 4.45 | 8.17 | 0.384 | 433 |
| 26 | 100 °C (3 reps) | 3.44 | 82.31 | 4.45 | 8.17 | 0.722 | 433 |
| 27 | 120 °C closed | 3.81 | 81.65 | 4.79 | 8.99 | 0.945 | 357 |
| 28 | 140 °C closed | 4.44 | 82.60 | 2.88 | 7.59 | 0.362 | 522 |
| 29 | 160 °C closed | 3.72 | 80.78 | 2.81 | 7.67 | 0.154 | 611 |
| 30 | 100 °C, 15 min | 1.78 | 91.05 | 3.30 | 2.34 | 0.013 | 1380 |
| 31 | 30 min | 3.44 | 82.56 | 5.42 | 4.51 | 0.164 | 358 |
| 32 | 1 h | 3.52 | 82.23 | 5.25 | 6.62 | 0.434 | 358 |
| 33 | 2 h (8 reps) | 3.81 | 81.65 | 4.79 | 8.99 | 0.945 | 357 |
| 34 | 3 h | 2.81 | 77.89 | 4.36 | 10.70 | 1.162 | 484 |
| 35 | 4 h | 3.18 | 77.62 | 3.82 | 10.01 | 1.253 | 485 |
| 36 | 6 h | 4.00 | 78.63 | 3.36 | 9.31 | 0.746 | 450 |
| 37 | 18 h | 3.79 | 78.20 | 3.27 | 9.60 | 0.681 | 483 |
| 38 | 30 h | 3.03 | 74.64 | 4.27 | 12.38 | 0.724 | 421 |
| 39 | 48 h | 2.66 | 73.74 | 4.46 | 12.36 | 0.567 | 448 |
| 40 | 25 °C, 30 days | 63.94 | 23.02 | 11.41 | 1.63 | 0.007 | 0.73 |
| 41 | -5 °C, 30 days | 71.89 | 15.34 | 10.46 | 1.76 | 0.005 | 0.31 |
| 42 | 0.1 M NH3 + NaOH | 3.11 | 79.16 | 4.10 | 9.00 | 0.056 | 481 |
| 43 | 8 M + NaOH, 2 h | 3.75 | 80.91 | 3.95 | 10.10 | 1.048 | 432 |
| 44 | 8 M + NaOH, 18 h closed | 2.73 | 76.93 | 3.97 | 12.16 | 5.469 | 534 |

Text values [M]: best closed yield 5.47 % (100 °C, 18 h, NaOH); best open 1.25 % (100 °C, 4 h); MPZ 74.06 % (exp 12)
to 86.19 % (exp 2) [the table prints 74.09 and 86.18; small text-table mismatch]; 2-ethylpyrazine ~0.344 % of total;
-5 °C PZ 71.87 %, MPZ 15.34 % [M].

**Flags.**
1. **Exps 17, 27 and 33 are identical in every cell read** (% and mg of peaks 1-5, % of peak 6, yield 0.945; peaks
   7-10 not cross-checked) [!]. 17 and 33 are
   the same condition (8 M, 100 °C, 2 h, 8 replicates; 33 is the time series' 2 h point), so their identity is a reuse.
   **27 is 120 °C in a closed ampoule**: its row is very probably a copy of 17, so the "maximum at 120 °C" rests on a
   duplicated row and the 120 °C condition has no independent data.
2. **Exps 25 (90 °C) and 26 (100 °C) have identical % in peaks 1-6 but different mg and yields** (0.384 vs 0.722) [!]:
   a copied distribution or a coincidence; R at 90 °C is not independent.
3. **Table I's NH3 column prints "0.8" for exps 21-41** [!], yet exp 33 (listed as 0.8 M) equals exp 17 (8.0 M) and
   differs from exp 8 (0.8 M, 100 °C, 2 h: R 635, yield 0.484), and the Fig. 1 caption is 8 M / 100 °C / 2 h. So the
   T and time series are most likely 8 M and "0.8" a misprint. Read conservatively: concentration uncertain for 21-41.
4. pH never measured; NH3 escapes from open vessels (8 M cannot stay dissolved at 100 °C) [M].

## 4. Mechanism evidence and its kind

- **Within-study factors (kind 6), one lab:** NH3 concentration, temperature, time, NaOH. Distribution nearly
  invariant at >= 60 °C (MPZ 74-92 %) [Q]; changes sharply below 50 °C (PZ dominant).
- No labels, no fed intermediates. The amine-carbonyl condensation scheme and the low-T "lower activation energy"
  reading for PZ are interpretations [Q].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 / P1 | CONTRADICTS (same as shibamoto1977, replicated within lab) | R 310-3200 at >= 60 °C |
| P1x / P1y vs P2b | LEAVES UNSEPARATED | no labels. The 25 -> 100 °C rise of R (~3 dex) constrains both (sec. 1, answer 2) |
| CF-087 / rivals' floor R >= 4h, h >= 0.1 | the -5 / 25 °C rows (R 0.3-0.7) sit at or below the floor 0.4 | low yields; read with care |
| evidence map 2.7 "MPZ extra route: NONE" | LEAVES UNCHANGED | no kind 1-4 evidence |
| E4 (pH law) | LEAVES UNCHANGED | no pH |
| shibamoto1977 dossier, O2 (sec. 1, item 4 there: "air-filled ampoules") | CONTRADICTS that reading if 1977 used this paper's method | here closed pots are N2-flushed before sealing |

## 6. Disagreements found

1. **With the shibamoto1977 dossier, O2 in closed vessels**: that dossier says every pot had air ("air-filled
   ampoules"); this paper (the probable method source, 1977's ref 3) bubbles N2 for 5 min before sealing. 1977's own
   text says only "sealed ampoules". Consequence: 1977's closed pots were probably N2-flushed (not strictly anoxic), so
   its "O2 not needed" claim has a closed-vessel arm with reduced O2. Which reading is right needs 1977's ref 3 confirmed.
2. **Same lab, glucose + NH4OH: R ~360-530 here vs 14.6 in wong1988** (1 M NaOH, 120 °C, 4 h, capillary columns,
   continuous extraction). See wong1988 sec. 6.1 for candidate causes; this paper's own 120 °C row is a duplicate (flag 1),
   so temperature cannot be excluded from here.
3. **Two labs, glucose + inorganic N: R ~360-530 here (8 M NH3) vs ~0.07 for NH4Cl-glucose (koehler1969, glycol-water,
   120 °C)**; Koehler 1970's glucose + NH4OH (MPZ 1415, DMP 56 µmol; no PZ) matches the MPZ dominance here.
4. **Temperature of maximum yield: 120 °C here (2 h) vs rising to 150 °C in koehler1970 (24 h, Asn)** — and the 120 °C
   value here is a duplicated row (flag 1).

## 7. Ledger-ready measurements

1. Table II, exps 1-20 (NH3 ladder) and 30-39 (time series): PZ, MPZ, 2,5-DMP, 2,6-DMP weight % + yield — key with
   `quantification: weight_percent_rf_corrected_recovery_checked_PZ_MPZ`; exps 17 / 33 once only.
2. Exps 40-41 vs 33 (-5 / 25 °C 30 d vs 100 °C 2 h): the temperature contrast of R; low-yield caveat.
3. Exps 42-44 (NaOH) — base factor on yield (0.002 -> 0.056 % at 0.1 M NH3).
4. Do NOT key exps 25-27 as independent temperature points until a second reader confirms the duplications (flags 1-2).

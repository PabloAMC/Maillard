# Shibamoto & Bernhard 1977 — EXTRACTION
### *"Investigation of Pyrazine Formation Pathways in Glucose-Ammonia Model Systems"*
### D-glucose 1 M + ammonium hydroxide 8 M (100 mL), 100 °C solution temperature, 2 or 18 h, open (reflux, air) or closed (sealed ampoule), with propyl gallate, BHT, O2 bubbling, CuCl2, ZnCl2 or NaOH added; nine pyrazines quantified as weight % of total pyrazines and total yield as weight % of glucose (Table II, 16 experiments). **No labels, no fed intermediates, no amino acid.** Methylpyrazine is 77-81 % of all pyrazines in every pot.

**Read 2026-10-31 from `data/articles/shibamoto1977.pdf`** (11 pp., Agric. Biol. Chem. 41(1):143-153; publisher scan
with an OCR text layer). Prose from `pdftotext -layout`; **Tables I and II read from 140 / 220 dpi renders** (Table II
is printed rotated; the 220 dpi render was rotated and read cell by cell). Figs. 1-2 (proposed fragments, glucose
cycle) and Table III (reaction types) are schemes, not data. Every derived number below was computed in the
`maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Takayuki Shibamoto, Richard A. Bernhard (UC Davis, Dept. of Food Science and Technology) |
| Venue | *Agric. Biol. Chem.* **41** (1), 143-153, 1977; received 9 Aug 1976 |
| DOI | not printed. Crossref 2026-10-31: **10.1271/bbb1961.41.143** resolves to this title, vol. 41, pp. 143-153, 1977 |
| Methods | "described previously" (ref 3: Shibamoto & Bernhard, J. Agric. Food Chem., "in press"; not on disk). Recovery check printed here: pyrazine and 2-methylpyrazine (0.25 g each in 100 mL water, 100 °C, 2 h) recovered at "over 98 %" [M] |

## 1. What this paper answers for the model

Context: EX-052 (`docs/history/pyrazine_rivals_2026-10-31/`). The mixing ratio R = MPZ² / (PZ x DMP) is <= 4 for one
shared condensation constant with the statistical factor (P1), = 1 in the engine as coded (CF-077), any value >= 4h
for pair-specific constants (P1x) or pair-specific DHP fate (P1y), and any value for an extra MPZ route (P2b). The
prereg names this paper as the source of P2b's "dehydration, no oxidant" variant.

| need | what Shibamoto 1977 gives | tag |
|---|---|---|
| pyrazine ratios for R | **yes, all three species, 16 pots.** Weight % of total pyrazines: PZ 2.73-4.57, MPZ 76.93-81.13, 2,5-DMP 2.97-4.75, 2,6-DMP 8.16-12.16 | [M] |
| R (molar; MW factor 80.09 x 108.14 / 94.11² = 0.978) | **R(2,5-DMP) = 322-631** (control pot 13: 437); with 2,5 + 2,6 summed as "DMP", R = 83-138 | [Z] |
| response factors | PZ and MPZ recovered > 98 % by the same work-up; 2,5-/2,6-DMP recovery not tested; detector and calibration in ref 3 (not on disk). Weight basis, so the molar conversion is the MW factor only | [M] / [!] |
| 13C / 15N labels | none | — |
| aminoketone condensation k, pH law, Ea | none (no time series beyond 2 vs 18 h; one temperature) | — |
| pyrrole vs pyrazine | none measured | — |
| peptide vs free amino acid | not applicable (NH3, no amino acid) | — |
| O2 requirement of aromatisation | O2 bubbling: 89.4 % of the open control's yield (2 h); propyl gallate 32.2 % (2 h) but 98.7 % at 18 h (control pairing in doubt, sec. 3); BHT 91.6 % (18 h, closed). Authors: dihydropyrazine aromatisation "does not appear to require oxygen" [Q] | [M] / [Q] |
| ring-forming route for MPZ | proposed: alpha-amino-alpha-hydroxyacetaldehyde + aminoacetone, aromatised by loss of water (hydroxy-dihydropyrazine dehydration), not by oxidation | **[P]** |

**Answers.**

1. **This is the most extreme mixing-bound violation on disk.** R = 322-631 with 2,5-DMP, R >= 83 even if every
   2,6-DMP were counted as an aminoacetone dimer. Leahy's lysine pots (R 7-70) and the Ho rows (R up to ~14) are far
   below it. Spread across 16 pots is small (CV of R ~20 %), so it is a property of the system, not of one run.
   Confidence the violation is real (not a response-factor artefact) ~0.85: PZ and MPZ recoveries were checked at
   > 98 %, and a DMP recovery deficit large enough to fake R ~ 400 vs 4 would need DMP at ~1 % recovery while its
   lower homologues are at 98 %.
2. **It does NOT separate P1x / P1y from P2b.** All three accommodate any R >= 4. Under P1x, R = h (1 + u)² / u
   needs **u >= ~320** (aldehyde-vs-ketone pairing selectivity) in these pots [Z]; a chemistry-fixed u that large would
   then predict R/h >= 320 in every pot, so the Leahy and Ho rows would need h ~ 0.02-0.2 [Z]. That is a sharpened
   prediction for P1x, not a rejection. P2b with e = extra / pairing MPZ flux needs e >= ~9 (R = 4h(1+e)², h <= 1) [Z].
   The authors' MPZ mechanism (a hydroxylated C2 aminocarbonyl + aminoacetone, dehydration) is a P2b-type route, but
   it is **proposed only** [P]: no intermediate was isolated ("we have been unsuccessful", p. 151 [Q]).
3. **Transfer to the engine is weak.** 8 M NH3 (pH ~11-12, not measured [!]), 1 M glucose, no amino acid, no Strecker
   step. Aminocarbonyls here come from amination of sugar fragments (glycolaldehyde, hydroxyacetone, pyruvaldehyde,
   glyceraldehyde) by NH3 [P], so the C2 pool feeding MPZ need not be the C2 pool feeding PZ. That is exactly the
   two-pool heterogeneity that lets R exceed any single-pool bound, and it is untested here.
4. **The O2 test is weak** [!]: every pot had air (open reflux or air-filled ampoules); "no oxygen needed" was never
   tested by deaeration. Excess O2 did not raise the yield (89.4 %), which says O2 is not limiting at 100 °C in 8 M
   NH3, not that aromatisation is non-oxidative.

## 2. Conditions [M]

| item | value |
|---|---|
| Sugar / amine | D-glucose 1 M (18 g) + ammonium hydroxide 8 M, made up to 100 mL with distilled water |
| pH | not stated or measured [!]; 8 M aqueous NH3 is strongly basic [Z] |
| Buffer | none |
| Temperature / time | 100 °C solution temperature; 2 h or 18 h |
| Vessel / O2 | "open" = at atmospheric pressure (reflux condenser implied); "closed" = sealed ampoules. Exp. 4: O2 bubbled 5 min, then a ~500 mL O2 balloon on the condenser |
| Additives (Table I) | propyl gallate 0.1 g (Exp. 1-2); BHT 0.067 g (3); CuCl2 0.17 g (5-7); ZnCl2 0.13 / 0.68 g (8 / 9); NaOH 0.4 g (10-12); controls 13-16 (none) |
| Replicates | 1-8 per experiment (Table I); no SD printed |
| Analysis | methylene chloride extraction, GC quantification per ref 3 (not on disk); yields as weight % |

## 3. Quantitative results (Table II, p. 147) [M]

Weight % of total pyrazines for A pyrazine, B 2-methyl, C 2,5-dimethyl, D 2,6-dimethyl; last column total yield %
(weight of total pyrazines / 18 g glucose x 100). E-I (ethyl, 2,3-dimethyl, ethylmethyls, trimethyl) are each <= 2.1 %
and omitted here; the mg columns are printed too. R computed [Z].

| Exp | system, time, additive | A % | B % | C % | D % | total yield % | R (C only) | R (C + D) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | open 2 h, propyl gallate | 3.79 | 79.18 | 4.73 | 8.99 | 0.303 | 342 | 118 |
| 2 | open 18 h, propyl gallate | 3.24 | 78.81 | 2.97 | 11.71 | 0.672 | 631 | 128 |
| 3 | closed 18 h, BHT | 4.57 | 77.33 | 3.97 | 11.41 | 4.846 | 322 | 83 |
| 4 | open 2 h, O2 | 3.86 | 79.49 | 4.72 | 9.27 | 0.840 | 339 | 114 |
| 5 | open 2 h, CuCl2 | 4.05 | 80.86 | 3.16 | 9.09 | 0.607 | 500 | 129 |
| 6 | closed 2 h, CuCl2 | 3.91 | 80.02 | 3.28 | 9.29 | 0.795 | 488 | 127 |
| 7 | closed 18 h, CuCl2 | 3.61 | 79.52 | 3.63 | 11.11 | 3.867 | 472 | 116 |
| 8 | open 2 h, ZnCl2 0.13 g | 3.76 | 80.18 | 4.03 | 9.81 | 0.683 | 415 | 121 |
| 9 | open 2 h, ZnCl2 0.68 g | 3.45 | 81.01 | 3.75 | 10.10 | 0.630 | 496 | 134 |
| 10 | open 2 h, NaOH | 3.75 | 80.91 | 3.95 | 10.10 | 1.047 | 432 | 122 |
| 11 | closed 2 h, NaOH | 3.82 | 79.28 | 4.15 | 9.13 | 1.047 | 388 | 121 |
| 12 | closed 18 h, NaOH | 2.73 | 76.93 | 3.97 | 12.16 | 5.469 | 534 | 131 |
| 13 | open 2 h, control | 3.86 | 80.43 | 3.75 | 8.16 | 0.940 | 437 | 138 |
| 14 | open 18 h, control | 3.80 | 80.51 | 4.75 | 8.81 | 0.951 | 351 | 123 |
| 15 | closed 2 h, control | 3.79 | 78.20 | 3.27 | 9.60 | 0.681 | 483 | 123 |
| 16 | closed 18 h, control | 3.71 | 81.13 | 4.19 | 8.67 | 5.289 | 414 | 135 |

Authors' ratios to controls (text, p. 145) [M]: propyl gallate 32.2 % (2 h) and 98.7 % (18 h); BHT 91.6 %; O2 89.4 %;
CuCl2 64.6 / 83.6 / 73.1 %; ZnCl2 72.7 / 67.0 %; NaOH 111.4 / 110.1 / 103.4 %. Exp. 1 / 13: 0.303 / 0.940 = 32.2 %
(checked [Z]). **[!] Controls 14 and 15 look swapped between Table I and the text** [Z]: the text pairs Exp. 2 (open,
18 h) with 15 (Table I: closed, 2 h), and Exp. 6 and 11 (closed, 2 h) with 14 (Table I: open, 18 h). Every printed
percentage reproduces from Table II with those pairings (0.672 / 0.681 = 98.7 %; 0.795 / 0.951 = 83.6 %; 1.047 / 0.951
= 110.1 %), so either Table I's labels for 14 / 15 or the text's pairings are wrong. With Table I's labels the
matched-control ratios are 70.7 % (propyl gallate, 18 h), 116.7 % (CuCl2, closed 2 h) and 153.7 % (NaOH, closed 2 h).
The "antioxidant effect diminished greatly" at 18 h depends on which reading is right. R is unaffected. The 18 h closed pots' high totals (3.9-5.5 %) are attributed partly to better extraction (no emulsion) [Q].
Other numbers: MPZ / PZ molar ratio 14-24 (control 13: 17.7) [Z].

## 4. Mechanism evidence and its kind

- **Within-study factor (kind 6), one lab (UC Davis):** antioxidants, O2, metal ions, base. O2 not limiting; propyl
  gallate slows formation at 2 h only. The distribution pattern is nearly invariant to every additive (MPZ 77-81 %).
- **No labelling, no fed intermediate, no trapping, no isolated intermediate** (attempts failed, p. 151 [Q]).
- Ten alpha-aminocarbonyl fragments and the "glucose cycle" (Figs. 1-2) are **proposed** [P], built backwards from the
  pyrazines found. The dehydration-of-hydroxy-DHP route is argued from the absence of piperazines and the presence
  of ethyl groups [Q].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 (engine, R = 1) | CONTRADICTS (third system) | R >= 83 under any DMP reading |
| CF-077 (missing cross-pair factor) | LEAVES UNCHANGED | the factor 2 is a code-vs-docstring fact; with or without it, R = 1 or 4 is far below 322 |
| EX-052 P1 (one constant + factor 2) | CONTRADICTS (in NH3 systems) | R/4 = 80-160; a third laboratory beside Reineccius and the Ho network, but a different system |
| P1x vs P1y vs P2b | **LEAVES UNSEPARATED** | each admits any R >= 4; P1x would need u >= ~320 here; P2b e >= ~9 |
| "shared condensation constant probably wrong" (Leahy bound) | SUPPORTS (direction) | same channel, stronger violation; NH3 instead of amino acid, so it is not a replicate of Leahy |
| evidence map 2.7 "methylpyrazine extra route: NONE" | LEAVES UNCHANGED | no kind 1-4 evidence; a fitted-yield signature only |
| triage item 5 (pyrazines, pyrrole sink) | LEAVES UNCHANGED | no pyrroles measured |

## 6. Disagreements found

1. **R magnitude across systems** [Z]: Shibamoto R(2,5) ~ 440 (glucose + 8 M NH3, 100 °C) vs Leahy R ~ 7-70 (lysine +
   glucose / ribose, 75-95 °C) vs Van Lancker 2010 R ~ 0.3-3.5 (lysine + glucose, 130 °C, peak areas;
   `vanlancker2010_extraction.md`). Candidate causes: the N source (NH3 aminates sugar fragments directly, with no
   Strecker step, so the C2 and C3 aminocarbonyl pools come from different fragments); pH; quantification (weight %
   with recovery checked for PZ / MPZ vs peak areas with theoretical recoveries); 2,5 / 2,6 co-elution in the Ghent
   data (the 2,6-DMP that Shibamoto resolves is 2-3x the 2,5 isomer here).
2. **Self-contradiction, controls 14 / 15** (sec. 3): the text's control pairings match Table II only if Table I's
   labels for Exp. 14 (open, 18 h) and 15 (closed, 2 h) are swapped. Cause: a labelling error in one place; which one
   cannot be decided from the paper.
3. **2,6-DMP > 2,5-DMP** (ratio 1.8-4.0) [Z]. 2,6-DMP is not an aminoacetone dimer (it needs aminoacetone + a
   2-aminopropanal-type partner, or aldehyde addition); in Low 2007 ~45 % of it is formaldehyde addition. Any
   R computed from a co-eluting "2,5(6)-DMP" peak (Van Lancker, Scalone) understates the aminoacetone-dimer R.

## 7. Ledger-ready measurements

1. Table II, all 16 rows: PZ, MPZ, 2,5-DMP, 2,6-DMP weight % + total yield % (100 °C, glucose 1 M / NH3 8 M) — the
   cleanest R > 4 rows on disk with all three species; key with `quantification: weight_percent_recovery_checked_PZ_MPZ`.
2. Exp. 13 vs 4 (O2 bubbled vs open air, 2 h): 0.940 vs 0.840 % — a within-study O2 factor (air vs O2-enriched, not
   vs inert).
3. Exp. 1 vs 13 and Exp. 2 vs its 18 h open control: propyl gallate at 2 h vs 18 h. Key only after the 14 / 15
   label question (sec. 6.2) is settled by a second reader.

# Bemis-Young, Huang & Bernhard 1993 — EXTRACTION
DOI: 10.1016/0308-8146(93)90009-5 (printed as the PII "0308-8146(93)90009-5" in the PDF metadata; Crossref 2026-11-01 resolves to "Effect of pH on pyrazine formation in glucose-glycine model systems", Food Chemistry 46, 383-387, 1993, authors Bemis-Young / Huang / Bernhard)
### *"Effect of pH on pyrazine formation in glucose-glycine model systems"*
### D-glucose 0.5 M + glycine 1.5 M in water, unbuffered, initial pH 1.00 / 2.34 / 4.00 / 5.97 / 7.00 / 9.00 / 9.64 / 12.0 (NaOH or HCl), sealed pressure bottles, 120 °C oven, 5 h; continuous CH2Cl2 extraction at pH 10, capillary GC-FID with response factors. Table 2: relative weight % of 19 pyrazines per pH; Table 3: 13 non-pyrazines. **Relative shares only (no absolute amounts); end pH measured but not printed.**

**Read 2026-11-01 from `data/articles/bemis-young1993.pdf`** (5 pp., Food Chem. 46:383-387; Acrobat 3 scan with OCR).
Prose from `pdftotext -layout`; **Table 2 read from a 220 dpi render** (the text layer drops the pH 12.0 column);
Table 3 from the text layer, spot-checked on the render. No figure carries data (Fig. 1 is a proposed scheme). Every
derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Gwen L. Bemis-Young, Jeane Huang, Richard A. Bernhard (Food Science and Technology, UC Davis) |
| Venue | Food Chemistry 46, 383-387, 1993; received 17 Apr 1992, accepted 1 Jun 1992 |
| Laboratory | UC Davis (Bernhard): `ucdavis_shibamoto` in lab_ids.yml (same lab as shibamoto1976 / 1977, wong1988). The atlas string would need adding |

## 1. What this paper answers for the model

| need | what Bemis-Young 1993 gives | tag |
|---|---|---|
| PZ / MPZ / 2,5-DMP / 2,6-DMP per pot | relative weight % of identified pyrazines, 8 pH values, triplicate pots x triplicate injections, means only (no SD) | [M] |
| R (molar; weight-basis factor 0.978) | **pH 5.97: 0.020; pH 7.00: 0.087; pH 9.00: 0.54; pH 9.64: 3.3.** pH 1.00 / 2.34 / 12.0: PZ ND (and MPZ ND at 1.00 / 2.34), not computable; pH 4.00: PZ ND, MPZ 0.69 | [Z] |
| response factors | FID relative RFs and weight % by McNair & Bonelli [C]; identities by Kovats index + MS + co-injection for PZ, MPZ, 2,5- and 2,6-DMP (Table 1) | [M] |
| 2,5- vs 2,6-DMP | **resolved** (DB-Wax capillary): 2,6/2,5 = 0.03 (pH 5.97) to 0.27 (pH 9.64), 0.04 at pH 12 | [M] / [Z] |
| pH dependence | **yes: R rises 2.2 dex from pH 5.97 to 9.64** [Z]; MPZ share x2.6 (0.75 -> 1.93 %); PZ share falls (1.31 -> 0.06 %); 2,5-DMP flat (18-21 %) from pH 4 to 9.64, 53 % at pH 12 | [M] / [Z] |
| absolute pH slope (E4) | none: shares only. CH2Cl2-soluble extract 1.56-5.01 % of glucose (lowest at pH 5.97, highest at 12) [M], not pyrazine-specific | [M] |
| temperature | none (120 °C oven) | — |
| nitrogen source | glycine only (Strecker aldehyde = formaldehyde) | [M] |
| separates P1x / P1y from P2b? | **No.** See sec. 5 | [Z] |

**Answers.**

1. **Glycine-glucose R is far below 4 at pH 6-9 (0.02-0.54) and reaches 3.3 only at pH 9.64.** Every pairing rival
   has R >= 4h, so these rows need h << 0.1 or routes that make PZ or 2,5-DMP without the A x M pairing (extra PZ / DMP
   routes lower R; only extra MPZ routes raise it). Second laboratory (after the Ho network's glycine row, CF-087) with
   R << 0.4 in a glycine pot; the method here is a calibrated capillary FID with RFs, not near-LOQ SPME.
2. **R rises 160x with pH (5.97 -> 9.64)** [Z]. Under P1 / P1x / P1y, with pH-independent constants and fates, R moves
   only through h. A strong monotone pH trend needs pH-dependent selectivity (P1x with a pH-dependent u), pH-dependent
   DHP fate (P1y), or extra routes whose share depends on pH. It does not pick among them.
3. **2,3,5-Trimethylpyrazine is the largest pyrazine at every pH (19-47 %)**, and 2-ethyl-3,5-dimethylpyrazine is
   71 % / 55 % at pH 1.00 / 2.34 [M]. In a glycine pot the Strecker aldehyde is formaldehyde; trialkyl dominance fits
   aldehyde addition to dihydropyrazines (Reading's [2-13C]glycine labels, low2007 [C]), which consumes DHPs and moves
   the PZ : MPZ : DMP ratios away from pairing values.

## 2. Conditions [M]

| item | value |
|---|---|
| Reactants | D-glucose 0.05 mol + glycine 0.15 mol in 100 mL (0.5 M + 1.5 M) |
| pH | initial, set with 2 N NaOH or 10 % HCl: 1.00, 2.34 (glycine pKa1), 4.00, 5.97 (pI), 7.00, 9.00, 9.64 (pKa2), 12.0. **Unbuffered by design** [Q]. pH measured after heating but **not reported** [!] |
| Vessel / O2 | Wheaton pressure bottles, silicone gaskets, porcelain stoppers; air headspace (no deaeration stated) |
| T / time | 120 °C oven temperature (solution T not stated), 5 h |
| Work-up | cooled, end pH measured, adjusted to pH 10.0, continuous liquid-liquid extraction 16 h with 300 mL CH2Cl2, Na2SO4, rotary evaporation 40 °C / 3.3 kPa, residue weighed |
| GC | Varian 6000 FID, 30 m DB-WAX, 70 °C 8 min, 2 °C/min to 160 °C, hold; split 1:100. GC-MS VG Trio-2, 70 eV |
| Replicates | triplicate samples, triplicate analyses; means only |

## 3. Quantitative results (Table 2, mean relative weight % of identified pyrazines) [M]

| pH | PZ | MPZ | 2,5-DMP | 2,6-DMP | 2,3-DMP | 2,3,5-TrMP | 2-Et-3,5-DMP | TeMP | R(2,5) [Z] | R(2,5+2,6) [Z] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | ND | ND | 3.73 | 1.03 | ND | 18.72 | 71.36 | 2.13 | — | — |
| 2.34 | ND | ND | 3.03 | ND | ND | 21.52 | 55.43 | 1.19 | — | — |
| 4.00 | ND | 0.69 | 22.42 | 0.21 | ND | 28.66 | 5.31 | 7.29 | — (PZ ND) | — |
| 5.97 | 1.31 | 0.75 | 20.65 | 0.60 | ND | 38.53 | 2.45 | 10.41 | 0.020 | 0.020 |
| 7.00 | 0.79 | 1.14 | 18.43 | 0.65 | 0.63 | 37.06 | 3.89 | 13.17 | 0.087 | 0.084 |
| 9.00 | 0.48 | 2.22 | 18.58 | 3.84 | 1.94 | 39.55 | 4.07 | 12.03 | 0.54 | 0.45 |
| 9.64 | 0.06 | 1.93 | 18.55 | 5.03 | 2.84 | 47.42 | 4.70 | 10.23 | 3.3 | 2.6 |
| 12.0 | ND | 0.20 | 53.44 | 2.24 | 0.99 | 31.18 | 1.76 | 2.71 | — | — |

Other rows (2,3-diethyl, 2,3-diethyl-5-methyl, quinoxalines ...) as printed; 2,3-diethyl-5-methylpyrazine is 16.65 /
24.13 % at pH 2.34 / 4.00 [M]. PZ at pH 9.64 (0.06 %) is the smallest number entering R; a ±0.03 point error moves R
by x0.7-2 [Z]. Table 3 (non-pyrazines, relative weight %): acetol 63.21 % (pH 1) falling to 0.56 % (pH 12);
2-methylimidazole 83.93 % at pH 12; pyrrole 0.08-4.61 % at pH >= 7 [M]. Extract yield 1.56 % (pH 5.97) to 5.01 %
(pH 12) of glucose [M]. pH 5.97 pots stayed colourless [Q].

## 4. Mechanism evidence and its kind

- **Within-study factor (kind 6), one lab (UC Davis):** pH. Shares only.
- No labels, no fed intermediates. The 2,3-diethyl-5-methylpyrazine route (2-aminopropanal + 2-hydroxy-3-amino-4-
  hexanone from 2-hydroxybutanal + acetaldehyde) is **proposed** [P] (Fig. 1).
- [!] The pH 1-2.34 product set (2-ethyl-3,5-dimethylpyrazine 55-71 %, PZ and MPZ ND) is unusual; with relative shares
  and a small total, a few large peaks dominate.

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 / P1 (R <= 4h) | LEAVES UNCHANGED as a bound test (no row above 4); P0's R = 1 is crossed in both directions (0.02 -> 3.3) | R < 4 rows do not test the upper bound |
| **CF-087 (Ho glycine row, R < 0.19)** | **SUPPORTS the row being real chemistry, not a near-LOQ artefact** (second lab) | glycine-glucose R 0.02-0.09 at pH 6-7 with calibrated capillary FID. Weakens CF-087's "low-level data" explanation; leaves "segregated formation" and adds "extra 2,5-DMP / PZ routes in glycine pots" |
| the rivals' floor R >= 4h, h >= 0.1 | CONTRADICTS as a general floor | pH 5.97 / 7.00 need h ~ 0.005-0.02, or extra PZ / DMP routes |
| P1x / P1y vs P2b | LEAVES UNSEPARATED | no labels; a pH trend in R is allowed to all three if their parameters may depend on pH |
| E4 (pyrazine pH law; where the below-7 slope sits) | bears on it weakly | relative shares only: MPZ share rises x2.6 from pH 6 to 9.6 while 2,5-DMP share is flat, so MPZ and DMP do not share one pH law in this pot; no absolute slope |
| evidence map 2.7 (aldehyde addition) | SUPPORTS (direction) | trialkyl pyrazines dominate a glycine (HCHO) pot at every pH |

## 6. Disagreements found

1. **Same lab, same N source, alkaline: R(glycine) 3.3 at pH 9.64 / 0.54 at pH 9.00 here vs 5.4 (± 0.31 dex) in
   wong1988** (glucose 1 M + glycine 8 M + NaOH 1 M, 120 °C, 4 h; initial pH not printed, end pH ~3 units lower).
   Within ~1-2 sigma at the alkaline end; consistent direction (R grows with base).
2. **Glycine MPZ share: 0.2-2.2 % here vs koehler1970 ND (glycine-glucose, glycol, 120 °C 24 h)** — same direction
   (glycine is a poor MPZ source), different medium.
3. **pH direction of R differs between N sources across studies** [Z]: here R rises with pH (glycine). In wong1988 the
   high-pH NH4OH pot gives R 14.6 while the ammonium-salt pots (lower pH: 0.8 mol NH4+ with 0.1 mol NaOH) give 350-560,
   i.e. R falls with pH for NH3 / NH4+. Candidate causes: different N chemistry (glycine Strecker gives HCHO, and with it
   aldehyde addition; NH3 has none), pH-dependent fates, or work-up.

## 7. Ledger-ready measurements

1. Table 2, pH 5.97 / 7.00 / 9.00 / 9.64: PZ, MPZ, 2,5-DMP, 2,6-DMP relative weight % (RF-corrected, capillary) — a
   within-study pH ladder of R in a glycine pot; key as `relative_weight_percent_rf_corrected`, initial pH only.
2. Table 2, pH 1-12 TrMP and 2-ethyl-3,5-DMP shares — aldehyde-addition products vs pH (direction).
3. Extract yields 1.56-5.01 % vs pH — not pyrazine-specific; record only as context.

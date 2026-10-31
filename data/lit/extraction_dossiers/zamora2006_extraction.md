# Zamora, Navarro, Gallardo & Hidalgo 2006 — EXTRACTION
### *"Chemical Conversion of α-Amino Acids into α-Keto Acids by 4,5-Epoxy-2-decenal"*
### 4,5-Epoxy-2-decenal (ED, 0-31 µmol) + phenylalanine (25 µmol) in 500 µL acetonitrile / 50 mM buffer (2:1), 37 or 60 °C, up to 24 h: phenylpyruvic acid (PPA, the α-keto acid, by PFBHA / BSTFA GC-MS) and phenylacetaldehyde (PA, direct GC-MS). PPA > PA at 37 °C (~2.5x at pH 2-4); pH optimum 3 for both; both linear in ED; time courses at 37 and 60 °C; PPA degrades at 60 °C, and fed PPA + Phe gives benzaldehyde and PA. **All numbers figure-only; Seville; yields ~0.5-1.7 % of Phe.**

**Read 2026-10-31 from `data/articles/zamora2006.pdf`** (5 pp., JAFC 54:6101-6105; born-digital, clean text layer).
Figs. 1-6 read by me from 200 dpi renders (scratchpad `z06-2.png` - `z06-4.png`): read error ~ +-0.3 µmol/mmol
(Figs. 1-5), ~ +-1 nmol/µmol (Fig. 6), ~ +-0.3 h. No tables of results (Table 1 is retention indices and spectra).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Evidence map 2.10 (lipid carbonyls as Strecker donors, SINGLE Seville; zamora2006 named as queued) and the
Strecker layer's acid / keto-acid branch (the prereg's omega route: "e.g. oxidative decarboxylation of the
transamination keto acid; mechanism not asserted").

| model need | what Zamora 2006 gives | tag |
|---|---|---|
| a transamination branch (amino acid -> α-keto acid) competing with Strecker decarboxylation | **yes, measured:** PPA from ED + Phe, identified against the authentic compound; ~2.5x PA at pH 2-4, 37 °C, 21 h | [F] / [Q] |
| PA yield from a lipid donor | **~5 µmol/mmol Phe (0.5 %)** at 25 µmol ED + 25 µmol Phe, pH 3, 37 °C, 21-24 h; ~8-8.5 (0.8 %) at 60 °C | [F] |
| time course | Fig. 4 (37 °C) and Fig. 5 (60 °C), 0-24 h, PA and PPA | [F] |
| pH law | Fig. 1, pH 2-8, 37 °C, 21 h: both maximal at pH 2-4; PPA falls ~4x from pH 3 to 7-8, PA only ~1.7x | [F] |
| order in the donor | linear in ED, 0-31 µmol (r > 0.995 at 37 °C; > 0.996 at 60 °C); PPA / PA slope ratio 3.4 at 37 °C | [M] (r, ratio) / [F] |
| keto acid as an aldehyde source | fed PPA (0.5 µmol) + Phe (25 µmol), pH 3, 60 °C: BZA and PA rise linearly, ~2.4-2.7 % of the PPA each in 24 h | [F] -> [Z] |
| Strecker ACID (phenylacetic acid) | **not measured** | [M] (absence) |
| O2 | not controlled or stated | [M] |

**Answers.**
1. **A keto-acid exit of the amino acid exists with a lipid carbonyl donor, larger than the aldehyde exit at 37 °C**
   [F]. This is a second chemical route from the imine (no decarboxylation), the "pathway b" of Scheme 1 [P for the
   mechanism]. For the sugar Strecker layer it is analogue evidence only: it shows the transamination product can
   be the MAJOR amino-acid-derived carbonyl in a Strecker-type system, which makes a keto-acid intermediate a
   chemically credible carrier of Hofmann's oxygen-dependent acid route (the prereg does not assert one). It does
   **not** test that route: no PAA, no O2 variation, a different donor class.
2. **PPA -> PA conversion is slow** at 60 °C, pH 3: ~2.4 % of fed PPA to PA (and ~2.7 % to BZA) in 24 h [Z]. In the
   ED pots at 60 °C, PPA falls by ~64 % between 4 h and 21 h [F] — far more than the fed-PPA conversion explains, so
   most of the PPA lost goes elsewhere (not printed; the authors mention oxidation to BZA, oxalic acid not found).
3. **Yields are small:** PA 0.5 % (37 °C) / 0.8 % (60 °C); PPA up to 1.4 % (37 °C, 21-24 h) and 1.7 % (60 °C, 4 h,
   then falling) [F]. These are ~1.5 dex below the same laboratory's ED + Phe PA yield in hidalgo2005 (sec. 6).

## 2. Conditions [M]

- **ED** prepared from 2,4-decadienal with 3-chloroperoxybenzoic acid in chloroform, column-purified
  (hexane / acetone 95:5), checked by 1H / 13C NMR and GC-MS.
- **Pots:** 0-31 µmol ED + 25 µmol L-Phe in 500 µL acetonitrile-buffer (2:1): **Phe 50 mmol/L, ED 0-62 mmol/L**;
  buffers 50 mM sodium citrate (pH 2.15-6) or 50 mM sodium phosphate (pH 6-8) (Fig. 1 caption); 37 or 60 °C. Time
  courses (Figs. 4-5): 25 µmol ED + 25 µmol Phe, pH 3 citrate. Vessel, headspace and atmosphere not stated. The
  "pH" is that of the aqueous buffer before mixing with acetonitrile [Z: not stated otherwise].
- **Fed PPA:** 0.5 µmol PPA + 25 µmol Phe in 500 µL ACN-citrate (2:1), pH 3.0, 60 °C. [!] Fig. 6's caption says
  "0.5 µmol of 4,5-epoxy-2-decenal"; the Methods say phenylpyruvic acid. The Methods are the likelier truth.
- **PA and BZA:** 75 µL + 125 µL ACN-water + 25 µL IS (3-(Z)-nonenol), direct GC-MS (HP-5MS, injector 250 °C);
  calibrations 2-60 nmol PA, 5-60 nmol BZA, five levels, r > 0.99, CV < 5 %.
- **PPA:** 100 µL + 200 µL ACN-water + 50 µL IS (stearic acid), O-(2,3,4,5,6-pentafluorobenzyl)hydroxylamine 2 h at
  37 °C, dried, BSTFA 30 min at 60 °C, GC-MS; calibration 15-400 nmol, r > 0.99, CV < 5 %.
- **Replicates:** not stated; no error bars on any figure.

## 3. Quantitative results (all figure-only) [F]

**Fig. 1 — pH, 37 °C, 21 h** (µmol / mmol Phe; ED amount for this series not stated):

| pH | 2 | 3 | 4 | 5 | 6 (cit / phos) | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| PA | 4.9 | 5.1 | 4.5 | 4.2 | 3.6 / 3.8 | 3.4 | 2.9 |
| PPA | 12.1 | 13.5 | 12.8 | 8.7 | 6.0 / 5.5 | 3.6 | 3.6 |

Text [Q]: PPA ~2.5x PA at pH 2-4; "preference ... reduced at a higher pH".

**Fig. 2 (37 °C) and Fig. 3 (60 °C) — ED dose, pH 3, 21 h** (µmol / mmol Phe):

| ED mmol | 0 | 0.0013 | 0.0025 | 0.0125 | 0.025 | 0.031 |
|---|---:|---:|---:|---:|---:|---:|
| PA 37 °C | 1.9 | 1.9 | 2.0 | 3.4 | 5.1 | 6.0 |
| PPA 37 °C | 3.1 | 2.8 | 3.4 | 7.8 | 15.0 | 16.4 |
| PA 60 °C | 2.9 | 3.3 | 3.2 | 6.0 | 8.0 | 9.3 |
| PPA 60 °C | 2.2 | 2.4 | 2.6 | 3.9 | 6.0 | 6.7 |

Printed: linear, r > 0.995 (37 °C), r > 0.996 (60 °C); slope ratio PPA / PA = 3.4 at 37 °C [M]. [!] Non-zero
intercepts at 0 ED (PA 1.9-2.9, PPA 2.2-3.1 µmol/mmol): a Phe-only blank (or reagent / injector) signal of ~0.2-0.3 %.

**Fig. 4 — time course, 37 °C, pH 3, 25 + 25 µmol** (µmol / mmol Phe):

| h | 0 | 0.5 | 1 | 2 | 3 | 6 | 14 | 18 | 21 | 24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PA | ~2.4 | 2.9 | 3.2 | 3.6 | 4.0 | 4.4 | 4.9 | 5.1 | 5.2 | 5.4 |
| PPA | ~2.4 | — | — | — | — | 4.9 | 10.8 | 12.2 | 13.5 | 14.3 |

Printed slopes [M]: PPA linear over the whole period, slope 0.530; PA 0.541 over the first 3 h, then 0.075 (units:
µmol / mmol Phe per h [Z: by the axes]).

**Fig. 5 — time course, 60 °C, pH 3, 25 + 25 µmol:** PA 3.0 (0 h), 6.3 (1 h), 7.0 (2 h), 7.3 (3 h), 7.6 (4 h), 8.0
(21 h), 8.5 (24 h); PPA 2.0, 7.5, 10.3, 15.6, **17.4 (4 h)**, 6.2 (21 h), 6.5 (24 h). PPA linear for 3 h (r = 0.993)
[M], then falls.

**Fig. 6 — fed PPA (0.5 µmol) + Phe (25 µmol), pH 3, 60 °C** (nmol / µmol PPA): PA 18.5 (0 h), 23 (3 h), 38.5
(21 h), 42 (24 h); BZA 23.5, 28.5, 43.5, 50.5. Linear, r > 0.992 [M]. [!] ~2 % of the PPA already reads as PA and
BZA at t = 0 (impurity, or conversion in the 250 °C injector).

## 4. Mechanism evidence and its kind

Laboratory: **Instituto de la Grasa (CSIC), Seville — Zamora / Hidalgo.**
- **Fed intermediate (kind 2):** ED fed to Phe (PA and PPA measured); PPA fed to Phe (PA, BZA measured).
- **Within-study factors (kind 6):** pH (2-8), donor dose, temperature (37 / 60 °C), time.
- **Proposed [P]:** Scheme 1, pathway a (decarboxylation -> PA + unsaturated hydroxylamine -> 2-pentylpyridine) vs
  pathway b (no decarboxylation -> PPA + the same amine). No labelling, no amine or CO2 balance; 2-pentylpyridine
  not quantified here. Whether the fed-PPA PA comes from PPA decarboxylation or from PPA acting as a Strecker donor on
  Phe is not distinguished (same product) [!].
- [!] Possible artefact not addressed: PA measured by direct injection at 250 °C in samples holding 2-3x more PPA;
  any injector decarboxylation of PPA would read as PA (Fig. 6's t = 0 PA suggests some).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| Evidence map 2.10 (lipid donors, SINGLE Seville) | **LEAVES SINGLE** | same laboratory; adds the keto-acid exit |
| Forced-input Strecker layer: the omega (O2) acid route, mechanism not asserted (prereg sec. 1.1) | **CONTEXT, weak** | a transamination keto acid can dominate the amino-acid carbonyls in a Strecker-type system; no PAA, no O2, lipid donor: not a test |
| Forced-input layer SINGLE-LAB (EX-049 / EX-051) | **LEAVES UNCHANGED** | no sugar donor |
| hidalgo2013 (α-keto acids as Strecker donors; same lab) | **SUPPORTS the loop** lipid -> keto acid -> aldehyde, weakly | fed PPA + Phe gives PA at ~0.1 % / h at 60 °C |
| CF-073 | **LEAVES UNCHANGED** | — |
| CF-074 / CF-075 | **LEAVES UNCHANGED** | Hofmann-internal |
| Lipid lane (alt-protein matrices) | **CONTEXT** | at pH 3-8 the lipid-donor PA yield in this system is 0.3-0.5 % per Phe in 21 h at 37 °C |

## 6. Disagreements found

1. **Within Seville, ED + Phe -> PA at 37 °C, 21 h, 50 mmol/L each: ~5.1 µmol/mmol here vs 174 +- 18 in hidalgo2005
   Table 2 (34x, 1.53 dex).** Conditions: here ACN-50 mM citrate 2:1, pH 3, GC-MS; there ACN-water 2:1, unbuffered pH
   ~5, GC-FID. The pH series here (PA 4.2 at pH 5) rules out pH as the cause. Candidates: citrate buffer (catalysis or
   inhibition), ED batch / purity, detector / calibration, vessel and headspace, blank treatment. Zamora 2015
   (140 °C, sand, pH 3, air: 12.3 % from ED; `zamora2015_extraction.md`) sits nearer hidalgo2005. **New; a candidate
   for the conflicts register** (within-lab, same donor, same nominal temperature and time).
2. **pH optimum of ED + Phe -> PA: pH 2-3 here (nearly flat, 4.9 -> 2.9 over pH 2-8) vs a Gaussian centred at pH
   6.9 in hidalgo2004** (aqueous 0.3 M citrate / phosphate / borate, overnight, 37 °C, GC-MS area ratios;
   `hidalgo2004_extraction.md`). Same laboratory. Candidate causes: solvent (2:1 acetonitrile shifts the effective pH
   and solubility), buffer strength (50 mM vs 0.3 M), and that hidalgo2004 is area ratios. New.
3. `hidalgo2015_extraction.md` row (Zamora 2006b: PPA higher than PA overnight at 37 °C) and
   `hidalgo2013_extraction.md` (companion: amino acid -> α-keto acid by ED) **agree** with the PDF.
4. Internal [!]: Fig. 6 caption (ED) vs Methods (PPA) for the 0.5 µmol reagent.

## 7. Ledger-ready measurements

- **Figs. 4-5** (PA and PPA vs time, 37 / 60 °C, pH 3): the only lipid-donor Strecker / keto-acid time courses on disk;
  key two readers. PA's fast-then-slow shape (first 1-3 h) is a kinetic signature for any lipid-donor arm.
- **Fig. 1** (pH 2-8, both products): within-study pH factor; carry the hidalgo2004 conflict.
- **Figs. 2-3** (dose, linear): order in the donor (first order in ED at fixed Phe) and the blank intercepts.
- **Fig. 6** (fed PPA -> PA, BZA, 60 °C): a keto-acid -> aldehyde conversion rate (~0.1 % of PPA per h to PA).

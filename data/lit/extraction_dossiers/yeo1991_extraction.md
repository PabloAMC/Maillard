# Yeo & Shibamoto 1991 (J. Agric. Food Chem. 39:370-373) — EXTRACTION (file stem `yeo1991`)
### *"Microwave-Induced Volatiles of the Maillard Model System under Different pH Conditions"*
### L-cysteine + D-glucose, **0.05 + 0.05 mol in 50 mL water (1.0 M each)**, unbuffered, adjusted to **pH 2 / 5 / 7 / 9** with 6 N HCl / NaOH, **microwave 700 W "high", 15 min, temperature not measured**, heated until a "brown mass" formed. 22 volatiles by GC-FID as **peak-area ratio to undecane** (means of duplicates), no response factors, no recovery. **2-methyl-3-furanthiol (MFT) and norfuraneol (NF) are not reported at any pH.** 2-furfurylthiol (FFT) area ratio **0.18 (pH 5), 0.22 (pH 7), not detected at pH 2 and 9.**
DOI: 10.1021/jf00002a029 (not printed in the PDF; Crossref 2026-11-05: same title, authors Helen C. H. Yeo and Takayuki Shibamoto, JAFC 39(2) 370-373, 1991)

**Read 2026-11-05 (first reader) from `data/articles/yeo1991.pdf`** (4 pp., ACS scan with an OCR text layer). Prose from
`pdftotext -layout`; **Table I was checked cell by cell against a 150 dpi render of p. 371** (scratchpad
`m4pdf/yeo1991/p-2.png`); the OCR's "008" (2-methylpyridine, pH 7) and "0 .71" are 0.08 and 0.71 on the render. Fig. 1 was
read off the same render (marker centres against the 0-700 axis, 0.344 px per unit): **read error ±10 area units**
(±0.05 A420 on the right axis). Fig. 2 is a mechanism scheme. No supplementary material.

**Reservation (M4, case `docs/reference/cases/M4_glucose_mft_protein_amine.md` on branch `m4-case-refresh`, commit
8660c005, sec. 2.9).** The case fixed this paper's role before it was read: a candidate second laboratory for D2 (the
glucose + Cys MFT pH shape). Its thiol values below carry the mark
**RESERVED: unseen by the engine (M4); no engine run on these pots until M4's design is frozen** (abbreviated
**[RESERVED M4]** after this point). The role is not met (sec. 8), so release is the coordinator's call. No generator,
probe or prediction was run on any pot of this paper.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited by the authors · **[Z]** derived by me, never printed ·
**[!]** flag. Provenance: **first reader only, second reader pending.** Nothing here is in a ledger, fit table or
benchmark.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Helen C. H. Yeo and Takayuki Shibamoto*** (byline p. 370) |
| Affiliation | Department of Environmental Toxicology, University of California, Davis, CA 95616 |
| Venue | *J. Agric. Food Chem.* **39** (2) (1991) **370-373** |
| DOI | **10.1021/jf00002a029** (not printed; the page carries only the code 0021-8561/91/1439-0370; Crossref record matches title, authors, volume, pages) |
| Dates | received 10 Apr 1990; revised 6 Jul 1990; accepted 30 Jul 1990 |

**Laboratory.** UC Davis, Shibamoto group: `ucdavis_shibamoto` in the dossiers (`shibamoto1976_extraction.md`,
`shibamoto1977_extraction.md`). Not in `data/keys/laboratories.yml` as a group (no ledger row uses it). **[!] UC Davis
is a TEST laboratory for E13 and E14 (pyrazines) in split v2** (`split_v2_per_change.yml`, `test_labs`). This paper's
pyrazine rows (pH 9 only) are therefore not development evidence for M6; they are recorded below only because the
brief asks for every quantitative result.

## 1. What this paper answers for the model

| model need | what Yeo 1991 gives | tag |
|---|---|---|
| Glucose + Cys MFT at >= 3 pH (D2) | **none**: MFT is not among the 22 compounds of Table I, at any pH | [M] |
| Glucose + Cys NF | **none** (the furanone listed is 4-hydroxy-2,5-dimethyl-3(2H)-furanone, HDMF, not NF) | [M] |
| Glucose + Cys FFT vs pH | area ratio 0.18 / 0.22 at pH 5 / 7; nd at pH 2 and 9 (no LOD) | [M] |
| Total volatiles vs pH | Fig. 1: ~38 / 90 / 154 / 610 area units at pH 2 / 5 / 7 / 9 | [F] |
| Temperature, time course, rate constants | **none**: microwave heating, temperature not measured, one time | — |

## 2. Conditions [M] (p. 370)

| item | as printed | derived [Z] |
|---|---|---|
| Charge | "L-Cysteine (0.05 mol) and D-glucose (0.05 mol)" dissolved in ~30 mL deionized water, pH adjusted, made up to **50 mL** | **1.0 M Cys + 1.0 M glucose** (6.06 g + 9.01 g in 50 mL) |
| pH and scale | adjusted to **2, 5, 7, 9** with 6 N NaOH or 6 N HCl before heating; **label pH, room temperature; no buffer; pH after heating not reported** | at 1.0 M Cys and no buffer the pH will drift; the four arms are label pH only [!] |
| Heating | microwave oven (Quasar MQ 7796 AW, **700 W**, "high" setting), **15 min**, irradiation interrupted every 4 min and the sample rotated 90° | **temperature not measured or stated** [!]; 15 min "coincided with the onset of sample browning", charring after it; the product is a "brown mass", i.e. the water largely boiled off [Z] |
| Vessel / headspace | **not stated** (open or closed unknown) | — |
| Work-up | brown mass dissolved in 100 mL water, **adjusted to pH 8** (6 N NaOH), continuous liquid-liquid extraction with 50 mL CH2Cl2 for 6 h, Na2SO4 12 h, Vigreux to ~1 mL, N2 to 1 mL, **undecane 0.5 mL** [sic] added as GC internal standard, N2 to 0.5 mL | the "0.5 mL" of undecane into a 1 mL extract is as printed; probably a mass or a dilute solution was meant [!] |
| Replicates | "All data obtained in this study were collected in duplicate" | no per-value spread |
| Uniformity check | identical samples at different oven positions gave no significant variation in total yield [Q] | — |

## 3. The assay [M] (p. 370)

- GC-FID (HP 5890), 60 m × 0.25 mm DB-Wax, 60 °C (4 min) → 4 °C/min → 180 °C (30 min), split 1:45. GC-MS (VG Trio II),
  70 eV. Identification by mass spectra and Kovats indices against authentic compounds.
- **Quantity reported: "GC peak area of product / GC peak area of internal standard"**, mean of duplicates (Table I
  footnote a). **No response factors, no calibration, no recovery, no LOD.** "Not detected" (footnote b) has no stated
  limit. Values are therefore within-compound relative only across the pH arms; between compounds they are not molar.
- Browning: A420 (HP 8452A).

## 4. Data

### 4.1 Table I (p. 371): GC peak area ratio to undecane, means of duplicates [M]

Thiol rows: **[RESERVED M4]**.

| compound (as printed) | KI (DB-Wax) | pH 2 | pH 5 | pH 7 | pH 9 |
|---|---:|---:|---:|---:|---:|
| 4,5-dimethyloxazole | 1152 | nd | nd | nd | 13.44 |
| trimethyloxazole | 1206 | nd | 0.27 | 1.08 | 1.68 |
| 2-methylpyridine | 1234 | nd | 0.11 | 0.08 | 0.54 |
| thiazole | 1242 | 0.22 | 0.12 | 0.23 | 1.26 |
| 2,5-dimethylthiazole | 1354 | nd | nd | nd | 5.06 |
| 2-methylthiazolidine | 1415 | 1.85 | 0.56 | 0.69 | nd |
| 2-methylpyrazine [E13/E14 TEST lab] | 1263 | nd | nd | nd | 4.26 |
| 2,5-dimethylpyrazine [E13/E14 TEST lab] | 1316 | nd | nd | nd | 5.40 |
| 2,6-dimethylpyrazine [E13/E14 TEST lab] | 1319 | nd | nd | nd | 3.84 |
| 2-ethylpyrazine [E13/E14 TEST lab] | 1323 | nd | nd | nd | 3.81 |
| 2,3-dimethylpyrazine [E13/E14 TEST lab] | 1335 | nd | nd | nd | 4.94 |
| 2-ethyl-6-methylpyrazine [E13/E14 TEST lab] | 1363 | nd | nd | nd | 3.27 |
| trimethylpyrazine [E13/E14 TEST lab] | 1391 | nd | nd | nd | 10.75 |
| tetramethylpyrazine [E13/E14 TEST lab] | 1457 | nd | nd | nd | 2.92 |
| **furfurylmercaptan (FFT) [RESERVED M4]** | 1402 | **nd** | **0.18** | **0.22** | **nd** |
| 2-acetylfuran | 1483 | 5.69 | 31.19 | 41.49 | 5.46 |
| 2-furanmethanol | 1573 | nd | 0.11 | 0.26 | nd |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone (HDMF) | 2039 | 3.39 | 1.04 | 2.81 | 8.12 |
| **2-thiophenethiol [RESERVED M4]** | 1530 | 0.75 | 0.71 | 0.70 | nd |
| **3-thiophenethiol [RESERVED M4]** | 1555 | 0.04 | 0.32 | 0.47 | nd |
| 2-acetylpyrrole | 1952 | 2.79 | 3.60 | 5.64 | 3.88 |
| 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | 2266 | 0.92 | 4.31 | 13.01 | 76.10 |

**Absent from the table at every pH: 2-methyl-3-furanthiol, norfuraneol (4-hydroxy-5-methyl-3(2H)-furanone), furfural,
1-mercapto-2-propanone.** "Not detected" is footnote b; compounds never listed were either not detected or not
reported, and the paper does not say which [!].

### 4.2 Fig. 1 (p. 371): total GC area of volatiles and A420 vs pH [F, ±10 area units, ±0.05 A420]

| pH | 2 | 5 | 7 | 9 |
|---|---:|---:|---:|---:|
| total GC peak area | ~38 | ~90 | ~154 | ~610 |
| browning A420 | ~1.83 | ~0.00 | ~0.03 | ~0.95 |

Text: total area at pH 9 "12-fold greater" than at pH 2 [M]. **[!] The figure gives ~16× (610 / 38)**; the pH 2 marker
sits ~13 px above zero, so the ratio is sensitive to the read (12× would need ~51 units at pH 2, outside my read error).
The pH 2 sample was the darkest, pH 5 "very light yellow" [Q].

### 4.3 Arithmetic [Z]

- FFT pH 7 / pH 5 = 0.22 / 0.18 = **1.22 (+0.087 dex)**, within-study, same assay **[RESERVED M4]**. Duplicate spread
  is not printed, so the 1.22 is not distinguishable from 1.
- Authors' reading [Q]: C4 / C5 heterocycles (2-acetylpyrrole, 3-thiophenethiol, FFT, 2-acetylfuran) are higher at pH 5
  and 7 than at pH 2 and 9; C2-unit products (pyrazines, oxazoles) appear at pH 9 by base-catalysed retro-aldol
  fragmentation [C, Shibamoto & Bernhard 1977].

## 5. What it can and cannot test in the model

**Cannot:**
- Anything on MFT or NF (not reported). No D1, no D2.
- Any absolute level, rate constant or temperature law: microwave heating with no temperature, an evaporating charge,
  no calibration.
- The pH axis is label pH of an unbuffered 1 M Cys solution heated to a "brown mass"; the in-situ pH is unknown.

**Can (qualitative only):**
- FFT forms from glucose + Cys at label pH 5 and 7 and is below detection (no stated limit) at pH 2 and 9 under this
  heating; HDMF and 2-acetylfuran have a non-monotone pH pattern. Engine status: no FFT run exists on these pots.
- 2-methylthiazolidine (Cys + acetaldehyde) falls with pH and vanishes at pH 9 [M]: a Cys-carbonyl adduct signature,
  not an engine species.

## 6. Conflicts and agreements with existing records

- **Agreement in shape with Hofmann & Schieberle 1998 Table 2** (`hofmann1998_reconciliation.md` sec. 2.2; DFA, SIDA,
  145 °C, 0.5 M phosphate): glucose FFT 0.7 / 2.8 / 0.6 µg at pH 3 / 5 / 7, a maximum at pH 5. Yeo's FFT is nd at pH 2,
  present at 5 and 7, nd at 9. Not a test: different heating, no temperature, area ratios, and Yeo's pH 7 > pH 5 is
  within unknown duplicate noise. [RESERVED M4]
- **Partial agreement with `zhang1991b_extraction.md`** (Rutgers, glucose + Cys 50 / 50 mM, label pH 7.5, 180 °C): no MFT
  or FFT detected there; Yeo detects FFT but not MFT. Neither paper reports a detection limit.
- **[!] Internal:** the text's "12-fold" vs the figure's ~16-fold (sec. 4.2); "0.5 mL" of undecane added to a 1 mL
  extract (sec. 2).
- Compound keys: FFT (`2_furfurylthiol`) resolves; nothing was keyed (no data rows).

## 7. Provenance tiers (first reader only)

All Table I cells: `direct_measurement`, quant class **relative (area ratio to an internal standard, no response
factor)**, confidence low (duplicates without spread, temperature unknown). Fig. 1: `figure_digitised`, low. Every [Z]
value: `derived`. **Proposed ledger entries: none** (no MFT, no temperature, no calibration).

## 8. M4 role (case 8660c005)

**Role fixed before reading** (case sec. 2.9, candidate table): "decisive for D2 (a second laboratory) only if it is a
glucose + Cys pot with MFT quantified at >= 3 pH; otherwise evidence". The measurement that would count: "MFT at >= 3
pH, one design, stated scale".

| element of the criterion | what the paper has | verdict |
|---|---|---|
| a glucose + Cys pot | yes: 1.0 M glucose + 1.0 M Cys, water, one design | met |
| at >= 3 pH | yes: label pH 2 / 5 / 7 / 9 (unbuffered, adjusted with 6 N HCl / NaOH, room-temperature scale implied, not stated) | met |
| MFT quantified at each | **no: MFT is absent from Table I at all four pH** (and from the text) | **not met** |

**Verdict: NOT MET.** The paper is evidence only: FFT area ratio 0.18 (pH 5) and 0.22 (pH 7), nd at pH 2 and 9, under
microwave heating of unstated temperature. D2 keeps one laboratory (DFA, Hofmann 1998 Table 2).

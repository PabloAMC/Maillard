# Low, Koutsidis, Parker, Elmore, Dodson & Mottram 2006 — EXTRACTION
### *"Effect of Citric Acid and Glycine Addition on Acrylamide and Flavor in a Potato Model System"*
### Potato cakes (drum-dried flake : water or additive solution 1 : 1.3), baked at 180 °C in an impingement oven 10-60 min: control, + 1.3 % (and 0.39 %) glycine, + 0.39 % citric acid, + both (dry-weight basis). Acrylamide (SIDA, 13C3), Strecker aldehydes and diones (headspace, relative to IS), and **24 alkylpyrazines in µmol/kg dry weight by SDE with calibration and recovery correction** (Table 4). **The concentrations behind Low 2007's [2-13C]glycine isotopomer shares — but unsubstituted pyrazine was NOT quantified, so R cannot be computed.**

**Read 2026-10-31 from `data/articles/low2006.pdf`** (8 pp., JAFC 54:5976-5983; born-digital, clean text layer).
Tables 1-5 from `pdftotext -layout` (row alignment checked against the column letters of each row). Figs. 1-5
(acrylamide curves, schemes) not used for numbers. Derived numbers computed by hand, checked in the env where they
enter a claim.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Mei Yin Low, Georgios Koutsidis, Jane K. Parker, J. Stephen Elmore, Andrew T. Dodson, Donald S. Mottram (Reading) |
| Venue | *J. Agric. Food Chem.* **2006**, 54 (16), 5976-5983; published on web 6 Jul 2006 |
| DOI | **10.1021/jf060328x** (printed; Crossref 2026-10-31: same title, vol. 54, pp. 5976-5983) |
| Companion | Low, Parker & Mottram 2007, JAFC 55:4087 (`low2007_extraction.md`): same potato flake system with [2-13C]glycine at 1.3 % |

## 1. What this paper answers for the model

| need | what Low 2006 gives | tag |
|---|---|---|
| pyrazine ratios for R | **No.** "Unsubstituted pyrazine was found ... but it could not be isolated in sufficient quantities by SDE for quantification" (p. 5978) | [M] / [Q] |
| MPZ and 2,5-DMP (µmol/kg dry wt) | control: MPZ 80.6 / 385 / 395 and 2,5-DMP 77.7 / 327 / 301 at 15 / 30 / 60 min; + 1.3 % Gly: MPZ 30.2 / 96.8 / 182, 2,5-DMP 121 / 277 / 367 | [M] |
| MPZ / 2,5-DMP (under pairing = 2 A / M) | control 1.04 / 1.18 / 1.31; glycine 0.25 / 0.35 / 0.50; citric 1.42 / 1.29; Gly + citric 0.68 / 0.85 [Z]. Glycine cuts the C2 aminocarbonyl pool ~3-4x (MPZ falls, 2,5-DMP unchanged) | [Z] |
| response factors | SDE with 10-point calibration per compound (authentic standards; closest isomer otherwise), recovery by spiking; same-mass pyrazines assumed equal recovery [M]. Calibrated molar data — the best-quantified MPZ / DMP pair in this set | [M] |
| 13C / 15N | none here; Low 2007 gives the isotopomer shares for the 1.3 % Gly pot at 15 / 30 min | [C] |
| condensation k, pH, Ea | none; 180 °C only; three times | — |
| pyrrole vs pyrazine | not measured | — |
| peptide vs free amino acid | not applicable (free glycine added) | — |

**Answers.**

1. **No R.** The paper cannot test the mixing bound and cannot separate P1x / P1y from P2b.
2. **Absolute amounts for Low 2007's label shares** [Z]: at 30 min, + 1.3 % Gly, labelled 2,6-DMP ≈ 0.476 x 117 ≈ 56
   µmol/kg (SL + DL, `low2007_extraction.md` Table 1); labelled MPZ ≈ 0 x 96.8 (0 % at 30 min) and ≤ 0.048 x 30.2 ≈
   1.5 µmol/kg at 15 min. So glycine-formaldehyde addition makes tens of µmol/kg of 2,6-DMP from the
   AKG + AKM dihydropyrazine but ≤ ~1.5 µmol/kg of MPZ from the AKG + AKG dihydropyrazine. If both DHPs took HCHO
   with similar efficiency per molecule, the AKG-AKG DHP flux would be ≲ 0.03-0.1 of the AKG-AKM DHP flux, i.e.
   A/M ≲ 0.06-0.2 under pairing — and pairing then predicts MPZ / 2,5-DMP = 2 A / M ≲ 0.1-0.4, against 0.35 observed
   at 30 min. Consistent with pairing alone, but only under that equal-efficiency assumption (confidence ~0.4).
3. **Glycine lowers MPZ but not 2,5-DMP** [M]: under pairing this needs the AKG pool to fall ~3-4x; Low 2007's labels
   show glycine raises 2,3-butanedione (C1 from glycine) and the authors propose glycine chain-elongates MGO (Fig. 4)
   [P]. Without PZ we cannot check the predicted ~10x fall in PZ (PZ ∝ A²).

## 2. Conditions [M]

| item | value |
|---|---|
| Matrix | drum-dried potato flake (Maris Piper : Pentland Dell 1:1, McCain), deionised water or additive solution 1 : 1.3 w/w; cakes 3 mm x 73 mm, ~18 g; 15 cakes per batch pooled and milled |
| Precursors (Table 1, mmol/kg dry wt; set 1 / set 2) | total reducing sugars 26.9 / 33.0; asparagine 38.0 / 27.0; glycine (native) 0.42 / 0.54; total amino acids 106.4 / 83.0; AA / TRS 3.95 / 2.52. Added glycine 52.0 (0.39 %) or 173.2 (1.3 %); AA / TRS 10.38 / 7.77 with 1.3 % |
| pH (raw dough, 25 °C, slurry) | control 5.78 / 6.00; 1.3 % Gly 5.99 / 6.01; 0.39 % Gly 5.98; citric 5.15 / 5.16; Gly + citric 5.15 / 5.14 (Table 5). Not measured after baking |
| Temperature / time | 180 °C impingement oven; set 1: 10-60 min (12 points, acrylamide only, single measurements); set 2: 15 / 30 / 60 min, triplicate analyses of a pooled powder |
| Water / a_w | not measured; dough ~56 % added water at the start [Z]; end moisture not stated [!] |
| O2 | open oven, air |
| Analysis | acrylamide: dibromo derivative, GC-MS SIM, [13C3]acrylamide IS. Volatiles: dynamic headspace onto Tenax (Strecker aldehydes, diones; areas relative to 1,2-dichlorobenzene, RF 1) and SDE (pyrazines, calibrated, recovery-corrected); ZB-Wax 60 m; CV < 30 % |
| Set 2 design | no 0.39 % Gly arm; at 60 min only control and 1.3 % Gly |

## 3. Quantitative results [M]

### 3.1 Table 4, alkylpyrazines (µmol/kg dry wt; means of 3; letters omitted)

| compound | ID | 15 C | 15 Gly | 15 CA | 15 G+CA | 30 C | 30 Gly | 30 CA | 30 G+CA | 60 C | 60 Gly |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| methyl- | A | 80.6 | 30.2 | 121 | 74.4 | 385 | 96.8 | 264 | 230 | 395 | 182 |
| 2,5-dimethyl- | A | 77.7 | 121 | 85.5 | 109 | 327 | 277 | 204 | 271 | 301 | 367 |
| 2,6-dimethyl- | A | 32.0 | 35.1 | 40.2 | 50.5 | 139 | 117 | 95.9 | 148 | 139 | 194 |
| ethyl- | A | 51.3 | 13.1 | 65.5 | 38.2 | 224 | 53.5 | 179 | 127 | 252 | 87.2 |
| 2,3-dimethyl- | A | 10.9 | 29.4 | 11.7 | 30.5 | 57.9 | 124 | 35.8 | 105 | 61.7 | 211 |
| 2-ethyl-6-methyl- | A | 33.0 | 16.2 | 38.1 | 31.7 | 148 | 53.1 | 122 | 99.5 | 156 | 71.5 |
| 2-ethyl-5-methyl- | A | 25.6 | 27.5 | 27.0 | 30.8 | 105 | 69.5 | 80.4 | 85.4 | 107 | 80.4 |
| trimethyl- | A | 31.48 | 308 | 28.6 | 140 | 138 | 719 | 87.0 | 366 | 135 | 861 |
| 2-ethyl-3-methyl- | A | 21.0 | 20.2 | 22.0 | 23.5 | 88.9 | 61.8 | 68.3 | 70.7 | 91.0 | 82.2 |
| 2-ethyl-3,6-dimethyl- | A | 29.2 | 37.3 | 27.5 | 36.8 | 100 | 78.7 | 78.3 | 90.2 | 100 | 93.4 |
| 2-ethyl-3,5-dimethyl- | A | 6.40 | 37.2 | 6.38 | 19.6 | 32.9 | 104 | 25.8 | 69.1 | 34.8 | 125 |
| tetramethyl- | A | 0.35 | 11.8 | 0.31 | 3.80 | 2.46 | 48.8 | 1.62 | 17.4 | 2.62 | 60.5 |
| total alkylpyrazines | | 411 | 696 | 486 | 600 | 1800 | 1830 | 1290 | 1720 | 1830 | 2440 |

Minor rows (propyl-, diethyl-, isobutyl-, methylpropyl-, methylisobutyl-, diethylmethyl-; 0.06-11 µmol/kg) as
printed in Table 4. ID: A authentic standard; B library + LRI, closest-isomer calibration; C library only.

### 3.2 Table 3, Strecker aldehydes and diones (headspace, area relative to IS = 100)

2,3-butanedione: control 18 / 52 / 77, + Gly 42 / 71 / 57 (15 / 30 / 60 min); 2,3-pentanedione: control 20 / 30 / 26,
+ Gly 9 / 12 / 12. Strecker aldehydes fall 30-70 % with glycine (e.g. 2-methylbutanal 9700 → 4200 at 30 min).

### 3.3 Table 2, acrylamide (µg/kg dry wt; set 2 triplicate means)

control 2110 / 9480 / 5900; + 1.3 % Gly 409 / 2340 / 2790; + 0.39 % CA 1940 / 7750 / ND; Gly 0.39 % + CA 1150 / 4960 /
ND (15 / 30 / 60 min).

## 4. Mechanism evidence and its kind

- **Within-study factors** (kind 6), Reading: glycine (added amine) and citric acid (pH 6.0 → 5.15). Calibrated SDE
  concentrations, triplicates of pooled material.
- "x + x" vs "x + x + y" pathways (Figs. 2-3) and the glycine chain-elongation of MGO / hydroxyacetone (Figs. 4-5) are
  [P] here and tested by labels only in Low 2007.

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 mixing test | LEAVES UNCHANGED | no PZ |
| P1x / P1y vs P2b | **does not distinguish** | |
| P2a_aa (glycine HCHO adds to the AKG-AKG DHP → MPZ) | CONTRADICTS a large share (with Low 2007) | labelled MPZ ≤ ~1.5 µmol/kg vs labelled 2,6-DMP ~56 µmol/kg at the same time |
| CF-077 | LEAVES UNCHANGED | |
| evidence map 2.7 (Reading K1 row) | SUPPORTS, gives the denominators | Low 2007's shares become absolute amounts |
| acrylamide lane: glycine / citric acid mitigation | adds within-study factors | Table 2 |

## 6. Disagreements found

1. **Low 2007 Table 1 "change" column vs Low 2006 Table 4 letters** (the 2007 paper quotes this one): 44 of the 46
   filled cells (24 rows x 15 / 30 min, minus the two blank dione cells at 30 min) agree with this paper's letters. Two do not: **2-methyl-6-isobutylpyrazine at 15 min** ("nc" in 2007; 1.95 c vs 0.77 a,
   a significant fall, here) and **2,5-diethylpyrazine at 15 min** ("nc" in 2007; 0.63 a vs 1.08 c, a significant
   rise, here). Both tentatively identified compounds (B / C); cause: transcription in 2007. Neither enters the
   pyrazine lane.
2. **Low 2007's dossier says glycine lowers 2,3-pentanedione by 50-60 % "quoted from Low 2006"**: confirmed (−55 %,
   −60 %, −54 %). It does NOT say what happens to 2,3-butanedione at 60 min: it falls with glycine (77 → 57), while it
   rises at 15 and 30 min. The "glycine raises butanedione" claim holds only up to 30 min.

## 7. Ledger-ready measurements

1. Table 4 MPZ, 2,5-DMP, 2,6-DMP, 2,3-DMP, TMP, TeMP (µmol/kg, calibrated) for control vs 1.3 % Gly at 15 / 30 min —
   pairs with Low 2007's isotopomer shares to give labelled amounts.
2. MPZ / 2,5-DMP (control vs glycine): a within-study shift of the C2 / C3 aminocarbonyl balance.
3. Table 2 acrylamide (set 2), for the acrylamide lane's glycine and pH factors.

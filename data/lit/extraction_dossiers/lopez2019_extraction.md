# Lopez, Kerley, Jenkinson, Luckett & Munafo 2019 (J. Agric. Food Chem. 67:11444) — EXTRACTION (file stem `lopez2019`)
### *"Odorants from the Thermal Treatment of Hydrolyzed Mushroom Protein and Cysteine Enhance Saltiness Perception"*
### Enzymatic hydrolysate of freeze-dried button mushroom (eHMP), **± cysteine 1.00 g per 200 mL (~41 mM)**, **refluxed 4 h then baked open in a glass dish at 160 °C for 45 min** (to ~50 mL); 16 odorants by SIDA (GC-MS), reported **per kg of mushroom, single values, no SD, no LOD**. With cysteine: **2-furfurylthiol (FFT) 22, 3-sulfanylpentan-2-one 8, 1-(2-furyl)ethanethiol 4.1, 2-methyl-3-(methyldithio)furan 1.4 µg/kg**; none detected without cysteine. **2-Methyl-3-furanthiol (MFT) was neither identified nor quantified.** pH never stated for the hydrolysate pots. Not simulable by the engine.

**Read 2026-11-05 (first reader) from `data/articles/lopez2019.pdf`** (10 pp., ACS publisher PDF, clean text layer).
Prose and tables from `pdftotext -raw`, checked against the printed table layout. No figure carries a number the model
needs (Figs. 2-6 are sensory; Fig. 1 mass spectra; Fig. 7 a proposed pathway); none digitised. No supplementary
material on disk (none is announced in the paper).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

**Every thiol value in this dossier is unseen by the engine; no engine run on these pots.**

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Jordan Lopez, Trent Kerley, Lindsay Jenkinson, Curtis R. Luckett, **John P. Munafo, Jr.** (corresponding, jmunafo@utk.edu) |
| Affiliation | Department of Food Science, University of Tennessee, Knoxville, TN 37996, USA |
| Venue | *J. Agric. Food Chem.* **2019**, 67(41), 11444-11453 |
| DOI | **10.1021/acs.jafc.9b04153** (printed on every page) — matches the brief |
| Dates | received 3 Jul 2019; revised 12 Sep 2019; accepted 13 Sep 2019; published 8 Oct 2019 |

**File name.** The PDF is on disk as `lopez2019.pdf`, and the local download queue (`data/articles/DOWNLOAD_QUEUE.md`,
untracked) already lists it as `lopez2019.pdf` with the first-author correction; no rename or alias is needed (the stem
matches this dossier). **Laboratory:** Munafo group, U. Tennessee Knoxville — **not in `data/keys/laboratories.yml`**;
proposed id `utk_munafo`. A laboratory new to the repo.

## 1. What this paper answers for the model

| model need | what Lopez 2019 gives | tag |
|---|---|---|
| absolute FFT / MFT from a cysteine pot, second laboratory | **FFT 22 µg per kg mushroom** in eHMP + cys (one value, no SD); MFT **absent** from the odorant list | [M] |
| other thiols | 3-sulfanylpentan-2-one 8, 1-(2-furyl)ethanethiol 4.1 µg/kg (+ the disulfide 2-methyl-3-(methyldithio)furan 1.4 µg/kg) | [M] |
| cysteine dependence | all four sulfur odorants **not detected** without cysteine (Table 2 "0") | [M] |
| sugar charge of the matrix | IC-PAD, mg per g dried mushroom: mannose 118.0 ± 4.9, xylose 47.0 ± 2.6, glucose 17.1 ± 0.9, maltose 10.2 ± 1.1, fructose 6.0 ± 0.3, sucrose 2.3 ± 0.1 (n = 3) | [M] |
| defined model pot | mannose 100 mM + cysteine 33 mM, phosphate pH 5, 140 °C, 1 h, Parr reactor: FFT and 1-(2-furyl)ethanethiol **detected** (GC-O / GC-MS), **not quantified** | [M] / [Q] |
| saltiness (out of scope) | the 7-odorant aroma model raised perceived saltiness of low-sodium broth (2-AFC d′ 0.56, p 0.0025; n = 87) | [M] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Matrix (eHMP) | fresh immature *Agaricus bisporus* (cv. Delta), **lyophilised**, milled; **30 g powder + 1.5 mL *A. oryzae* protease + 540 mL deionised water**, 50 °C, 4 h, shaken; vacuum-filtered; 500 mL filtrate, held 12 h at 4 °C. Amino-acid profile **not measured** (Kimatu et al. cited for the hydrolysis pattern [C]) |
| Cysteine charge | **cysteine 1.00 g in 200 mL eHMP** (8.25 mmol [Z] → **~41 mM** at the start; ~165 mM nominal after concentration to ~50 mL, ignoring losses [Z]). Control: eHMP without cysteine. Third pot: cysteine 1.00 g in 200 mL water |
| pH | **not stated for any eHMP or cysteine pot** — not set, not measured [!]. The mannose / cysteine model is in "a phosphate buffer (pH 5, 100 mL)" (strength not stated; label) |
| Temperature / time / vessel | (1) **reflux 4 h with stirring** in a 500 mL round-bottom flask (~100 °C [Z]; open to air through the condenser); then (2) poured into a **46 × 31 cm glass dish, baked in a kitchen convection oven at 160 °C for 45 min**, to ~50 mL "viscous liquid"; stored overnight at 0 °C. Model: sealed Parr pressure reactor, 140 °C, 1 h, stirred |
| Headspace / O2 / metals | eHMP pots **fully air-exposed** (reflux, open oven dish); no chelator; the model reactor's headspace not stated |
| Replicates | **not stated for the odorant work**; sugar analysis n = 3 |

## 3. Assays [M]

- **Isolation:** 50 mL treated sample + 50 mL water; labelled standards (**50 µg each**) added before extraction;
  diethyl ether 2 × 100 mL; SAFE at 40 °C, 10⁻³ Pa; Vigreux to ~2 mL, N2 to ~200 µL; pentadecane-d32 as a
  concentration check.
- **SIDA by GC-MS** (Agilent 5973, EI 70 eV, **full scan m/z 50-450**, FFAP column); RF from fixed mixtures of labelled
  and unlabelled standard. Labelled standards quantified by GC-FID calibration.

| no. | analyte | internal standard | m/z (analyte / IS) | RF |
|---|---|---|---|---|
| 9 | **2-furfurylthiol (FFT)** | **d-9 (labelled FFT, purchased, Aromalab)** | 114 / 116 | 1.40 |
| 8 | 1-(2-furyl)ethanethiol | **d3-8, synthesised in-house** (furfural-d3 → MeMgBr → Lawesson's reagent; assay by GC-FID against unlabelled 1-(2-furyl)ethanol) | 95 / 98 | 1.37 |
| 6 | 3-sulfanylpentan-2-one | **d-9 (labelled FFT) — not an isotopologue** ("Isotopomers were not available for 6 and 17" [V]) | 75 / 116 | 1.42 |
| 17 | 2-methyl-3-(methyldithio)furan | **d-9 — not an isotopologue** | 160 / 116 | 0.28 |
| 11 | 3-(methylsulfanyl)propanal (methional) | d-11 | 104 / 107 | 0.82 |

- **No LOD / LOQ** for any analyte; "0" in Table 2 means not detected, with no floor given [!].
- **No replicate noise**: Table 2 prints single values, no SD, no n [!].
- [!] **Basis ambiguous.** Table 2's footnote says "per kg of raw mushrooms"; the methods say the calculation used "the
  amount of dried mushrooms present in the sample". Fresh vs dry differs ~×10 for button mushroom (~90 % water;
  handbook value, not from the paper).
- [!] The synthesised d3-8's purity / isotopic enrichment and its absolute amount rest on a GC-FID calibration against a
  different compound (the alcohol), so its absolute scale carries an unquantified error.

## 4. Data

### 4.1 Table 2 — SIDA concentrations (µg/kg mushroom, basis per sec. 3) and OAVs [M]

| no. | odorant | threshold in water (µg/kg) | eHMP | OAV | eHMP + cys | OAV | factor |
|---|---|---|---|---|---|---|---|
| 4 | 1-octen-3-one | 0.016 | 21.7 | 1356 | 21.5 | 1342 | 0.99 |
| 6 | **3-sulfanylpentan-2-one** | 0.19 | 0 | < 1 | **8** | 42 | — |
| 8 | **1-(2-furyl)ethanethiol** | 0.053 | 0 | < 1 | **4.1** | 78 | — |
| 9 | **2-furfurylthiol (FFT)** | 0.036 | 0 | < 1 | **22** | 610 | — |
| 11 | methional | 0.43 | 1.5 | 4 | 1.4 | 3 | 0.93 |
| 12 | 2-ethyl-3,5-dimethylpyrazine | 0.08 | 0.1 | 1 | 0.1 | 1 | 1.00 |
| 16 | butanoic acid | 2400 | 12.5 | < 1 | 10.5 | < 1 | 0.84 |
| 17 | 2-methyl-3-(methyldithio)furan | 0.4 | 0 | < 1 | 1.4 | 3 | — |
| 18 | 2-acetylthiazole | 10 | 12.4 | 1 | 0.5 | < 1 | 0.04 |
| 19 | phenylacetaldehyde | 5.4 | 30.4 | 6 | 4.1 | < 1 | 0.13 |
| 20 | 3-methylbutanoic acid | 490 | 348.4 | < 1 | 324.2 | < 1 | 0.93 |
| 28 | HDMF (furaneol) | 40 | 629.6 | 16 | 338 | 9 | 0.54 |
| 29 | p-cresol | 3.9 | 0.7 | < 1 | 4.9 | 1 | 7.00 |
| 31 | sotolon | 0.49 | 7.7 | 16 | 9.9 | 20 | 1.29 |
| 34 | indole | 11 | 30.6 | 3 | 83.5 | 8 | 2.73 |
| 35 | phenylacetic acid | 6100 | 334.9 | < 1 | 458.2 | < 1 | 1.37 |

**Unseen by the engine; no engine run on these pots** (rows 6, 8, 9, 17 and every other value).
[!] OAVs divide a per-kg-mushroom concentration by a threshold in water: the OAVs are not concentrations in any
tasted or smelled medium.

### 4.2 Table 1 — cAEDA, sulfur odorants [M]

FD factors (eHMP / eHMP + cys): FFT – / 64; 3-sulfanylpentan-2-one – / 16; 1-(2-furyl)ethanethiol – / 16;
2-methyl-3-(methyldithio)furan – / 16; methional 256 / 256; dimethyl trisulfide 4 / –; 2-acetylthiazole 64 / 16;
2-acetyl-2-thiazoline 16 / 16. **2-Methyl-3-furanthiol does not appear** among the 36 odorants (FFAP and DB-5 RIs
listed for all 36).

### 4.3 Arithmetic [Z]

1. **Cysteine charge:** 1.00 g / 121.16 g/mol = 8.25 mmol in 200 mL → 41 mM. 200 mL is 40 % of the 500 mL hydrolysate
   from 30 g powder → ~12 g dried mushroom per pot (assuming the filtrate carries the solubles evenly).
2. **FFT as a concentration (order of magnitude only):** 22 µg/kg × 12 g dry (if "per kg" is dry) = 0.26 µg, or × ~120
   g fresh-equivalent (if fresh, ~90 % water) = 2.6 µg, in ~50 mL of concentrate → **~5-50 µg/L (0.05-0.5 µM)**; yield
   on cysteine **~3 × 10⁻⁵ to 3 × 10⁻⁴ %**. The factor-10 ambiguity is the basis question of sec. 3.
3. Sugars per pot (if the 12 g dry estimate holds): mannose ~1.4 g (~7.9 mmol), xylose ~0.56 g (~3.8 mmol), glucose
   ~0.2 g: cysteine 8.25 mmol is of the same order as the reducing sugar.

## 5. What it can and cannot test in the model

**Simulable by the engine: NO.** The pot is a whole-mushroom enzymatic hydrolysate of unmeasured amino-acid and peptide
composition, with mannose (not a sugar the engine's trunk carries as such) as the main sugar, pH never stated, heated by
4 h reflux open to air and then an open-dish oven step at 160 °C that evaporates ~75 % of the water (a moving
water activity the engine does not model). The one defined pot (mannose + cysteine, phosphate pH 5, 140 °C, 1 h) is
**qualitative only**.

**Second-laboratory ABSOLUTE thiol evidence for H0 (whose refusal of absolute thiol levels for a new lab rests on a
~1 dex between-laboratory tau): NO, not even as a non-simulable matrix check.** Reasons: (i) the concentration basis is
ambiguous by ~×10 (fresh vs dried mushroom; sec. 3), (ii) single values with no SD, no n and no LOD, so no
within-laboratory noise to set against a between-laboratory one, (iii) the matrix volume and water loss are
approximate (~50 mL), (iv) MFT, the thiol most of the sulfur lane's absolute evidence concerns, is not measured, and
(v) 3-sulfanylpentan-2-one and the disulfide are quantified against labelled FFT, not their own isotopologues. **FFT is
measured with its true isotopologue (d-9)**, so its value is the least doubtful of the four; as a check it can only say
that FFT forms in the µg/kg range from cysteine + a mushroom hydrolysate under air, which no record disputes.

**Can:** qualitative facts — cysteine is required for FFT, 1-(2-furyl)ethanethiol, 3-sulfanylpentan-2-one and the
methyldithio-furan in this matrix; cysteine lowers HDMF (×0.54), phenylacetaldehyde (×0.13) and 2-acetylthiazole
(×0.04), consistent with cysteine competing for dicarbonyls and sugar-derived carbonyls [Z, direction only]; mannose +
cysteine at pH 5 makes FFT and 1-(2-furyl)ethanethiol.

## 6. Conflicts and agreements with the repo's records

1. **No conflict**: no ledger or fit record covers mushroom hydrolysates, mannose + cysteine, or these four sulfur
   odorants from this laboratory.
2. **Agrees in direction** with the corpus that cysteine suppresses HDMF / furaneol-type and Strecker carbonyls in
   mixed systems (the authors' "driven off, degraded, or possibly precursors" [Q]); not quantified against any record.
3. [!] Internal: concentration basis (sec. 3); OAV bases (sec. 4.1); the cysteine-alone pot was only smelled (egg-like),
   never quantified.

## 7. Proposed ledger entries — none (qualitative check only)

Not simulable, and not usable as absolute second-laboratory evidence (sec. 5). If a presence record is wanted:
`provenance_tier: direct_measurement`, `quant_class: sida` for FFT and 1-(2-furyl)ethanethiol (true isotopologues),
`calibrated` for 3-sulfanylpentan-2-one and the disulfide (surrogate IS), confidence **low** (basis ambiguity, no
replicates), `simulable: false`, flag `unseen_by_engine`. First reader only; second reader pending.

## 8. Supplementary material

None on disk and none announced (the paper has no Supporting Information section). Nothing to read.

## 9. Open questions it bears on

| item looked for | result |
|---|---|
| conditions: matrix, cysteine charge, pH and its scale, temperature, time, vessel | **FOUND except pH**: eHMP from 30 g dry mushroom / 540 mL; cysteine 1.00 g / 200 mL (~41 mM); reflux 4 h then 160 °C open dish 45 min; **pH NOT FOUND** for the eHMP pots (model pot: phosphate "pH 5", label) |
| SIDA-quantified FFT, with IS, LOD, noise | **FOUND (value, IS); NOT FOUND (LOD, noise)** — FFT 22 µg/kg, d-FFT, m/z 114 / 116, RF 1.40 |
| SIDA-quantified MFT | **NOT FOUND** — MFT not identified, not quantified |
| other thiols | **FOUND** — 1-(2-furyl)ethanethiol 4.1 (d3 isotopologue, synthesised); 3-sulfanylpentan-2-one 8 (surrogate d-FFT); disulfide 2-methyl-3-(methyldithio)furan 1.4 µg/kg (surrogate d-FFT) |
| simulable by the engine | **NO** |
| second-laboratory absolute thiol evidence for H0 | **NO** (basis ambiguity, no replicate noise, no MFT); H0's refusal of absolute thiol levels for a new laboratory is neither eased nor tested |

**Unseen by the engine; no engine run on these pots.**

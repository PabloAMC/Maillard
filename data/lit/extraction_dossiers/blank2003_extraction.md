# Blank 2003 — EXTRACTION (glucose + proline and fructosyl-proline, 100 mmol/L each, pH 6 / 7 / 8, reflux 1-4 h; AP and ATHP by isotope dilution)
### The only pH-resolved, isotope-dilution series for 2-acetyl-1-pyrroline and 2-acetyltetrahydropyridine from a sugar + proline pot in the corpus.

**Source on disk:** `data/articles/blank2003.pdf` (owner's download, 2026-09-08). Read-only extraction
from the pypdf text layer in the scratchpad; Tables 2-5 are clean and re-typed in full below, Table 1
(HPAEC waveforms) is instrument detail and omitted. The four odorant time courses (Figure 5) and the
Fru-Pro time courses (Figure 6) are FIGURE-ONLY; only the values the authors state in the text or in
Table 5 (the 2 h point) are carried. Repo status before this dossier: roadmap §5b lists 2-acetyl-1-pyrroline
and 2-acetyltetrahydropyridine with "data on disk: none".

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Odorants in Maillard Model Systems Based on L-Proline as Affected by pH" |
| Authors | Imre Blank, Stéphanie Devaud, Walter Matthey-Doret, Fabien Robert (Nestlé Research Center, Lausanne) |
| Venue | J. Agric. Food Chem. 2003, 51, 3643-3650 |
| DOI | 10.1021/jf034077t |
| Naming | AP = 2-acetyl-1-pyrroline; ATHP = 6-acetyl-1,2,3,4-tetrahydropyridine = "2-acetyltetrahydropyridine" (two tautomers; quantified on the first-eluting one); HDMF = 4-hydroxy-2,5-dimethyl-3(2H)-furanone = furaneol (repo `hdmf`); AcOH = acetic acid; Fru-Pro = N-(1-deoxy-D-fructos-1-yl)-L-proline, the Amadori compound of glucose and proline; Glc/Pro = glucose + proline |
| Companions | ref 16 = Hofmann & Schieberle 1998 (`hofmann1998b_extraction.md`, the fed 1-pyrroline yields this paper quotes); ref 14 = Schieberle 1995 (popcorn SIDA; 38 µg ATHP per mmol proline at pH 7, 2 h, boiling) |

## 1. Why it matters

The repository has no step for the proline odorants (roadmap §5b, §5c Programme 6). This paper gives,
from one laboratory with stable-isotope dilution, the AP and ATHP levels reached by glucose + proline
and by the Amadori compound at three pH values and three times, in a buffered aqueous pot at 100 mmol/L
each. It is an END-OF-COOK validation set for whatever proline step the next wave builds (the fed
1-pyrroline yields to fit against are in Hofmann 1998b), and it fixes the ORDER OF MAGNITUDE: AP is
0.002-0.003 mol % of the sugar and ATHP 0.01-0.04 mol % after 2 h, ten to a hundred times below HDMF
and four orders below acetic acid. It also carries a measured glucose-consumption time series (Table 3)
and the Amadori compound's formation and decay, which are trunk quantities for a proline pot.

## 2. Methods as they matter to a model

- **Odorant pots:** 2 mmol glucose + 2 mmol proline, or 2 mmol Fru-Pro, in 20 mL 0.2 M phosphate ->
  **100 mmol/L each** (the authors: "same initial concentration of the precursors (0.1 mol/L)").
  pH adjusted to 6, 7 or 8. **Refluxed** (atmospheric boiling; the temperature is not printed, ~100 C)
  for **1, 2 or 4 h**. Open reflux with condenser, so the pH drift is real (Table 2) and volatiles are
  retained. Each sample at least in duplicate; each injected twice.
- **Quantification: stable isotope dilution (IDA)**, standards added to the cooled mixture:
  [13C1]acetic acid 440 µg; [13C2]HDMF 36.5 µg; [2H2-4]AP 2.385 µg (0.05 mL of 47.7 µg/mL in propylene
  glycol, purity 95 %, tautomers 93 % AP + 7 % 2-acetyl-2-pyrroline); [2H2-5]ATHP 14.45 µg (0.05 mL of
  289.0 µg/mL, purity 85 %, tautomers 60 % 6-acetyl-2,3,4,5-tetrahydropyridine + 40 % ATHP). Acidic
  (pH 3) and basic (pH 10) fractions, Et2O 2 x 5 mL, dried, concentrated to ~1 mL. GC-MS EI: AcOH and
  HDMF on DB-Wax (m/z 60/61; 128/130); AP and ATHP on HP-PONA (m/z 111 vs 113-115; 125 vs 127-130),
  **ATHP quantified on the first-eluting tautomer**.
- **Precision:** "coefficient of variation lower than 20 % for concentrations (µg/mmol) above 300 (AcOH),
  5.0 (HDMF), 1.0 (ATHP) and 0.3 (AP)". All Table 5 AP and ATHP entries are above those floors.
- **Reporting unit:** µg per mmol precursor (glucose or Fru-Pro) in Figure 5; **mg/L** in Table 5.
  Conversion: 2 mmol precursor in 20 mL, so X µg/mmol = 2X µg / 0.020 L = 0.1 X mg/L, i.e.
  **µg/mmol = 10 x mg/L**. mol % of precursor = (µg/mmol) / (10 x MW). MW: AP 111.14, ATHP 125.17,
  HDMF 128.13, AcOH 60.05. Checks: 4.5 mg/L ATHP = 45 µg/mmol, the value the text quotes against
  Schieberle's 38; 325 mg/L AcOH = 3.25 mg/mmol = the text's "3.3 mg/mmol".
- **Nonvolatile pots (Table 3, Figure 6):** 5 mmol glucose + 5 mmol proline in 50 mL 0.2 M phosphate
  (**100 mmol/L each**), reflux up to 7 h, pH 6 and 7; or 5 mmol Fru-Pro (1.39 g) in 50 mL (100 mmol/L).
  HPAEC with pulsed electrochemical detection, calibration curves, duplicates, injected twice.
- **Sensory:** OAVs from thresholds in water (AcOH 50, HDMF 0.06, ATHP 0.0016, AP 0.0001 mg/L),
  thresholds not pH-adjusted.

## 3. Tables re-typed

### Table 2. "Evolution of the pH in the Maillard Reaction Samples" (phosphate 0.2 mol/L; measured on a cooled aliquot)

| time (h) | pH 6 Glc/Pro | pH 6 Fru-Pro | pH 7 Glc/Pro | pH 7 Fru-Pro | pH 8 Glc/Pro | pH 8 Fru-Pro |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 6.00 | 6.00 | 7.00 | 7.00 | 8.00 | 8.00 |
| 1 | 6.15 | 5.80 | 7.04 | 6.70 | 8.04 | 7.30 |
| 2 | 6.16 | 5.42 | 7.00 | 6.30 | 8.02 | 7.02 |
| 4 | 6.13 | 5.05 | 7.00 | 6.20 | 8.02 | 6.97 |

The Glc/Pro pots hold their pH; the Fru-Pro pots fall by about one unit, driven by acetic acid.

### Table 3. "Consumption of Glc in the Presence of Pro" (100 mmol/L each, 0.2 M phosphate, reflux; HPAEC; mean +/- as printed)

| time (min) | pH 6 concn (mmol/50 mL) | pH 7 concn (mmol/50 mL) | pH 6 loss (mol %) | pH 7 loss (mol %) |
|---:|---:|---:|---:|---:|
| 0 | 5.00 | 5.00 | 0 | 0 |
| 30 | 4.95 +/- 0.10 | 4.47 +/- 0.10 | 1 | 11 |
| 60 | 4.93 +/- 0.06 | 4.75 +/- 0.33 | 1 | 5 |
| 90 | 4.59 +/- 0.10 | 4.49 +/- 0.10 | 8 | 10 |
| 120 | 4.74 +/- 0.09 | 4.44 +/- 0.24 | 5 | 11 |
| 180 | 4.66 +/- 0.05 | 4.02 +/- 0.15 | 7 | 20 |
| 240 | 4.65 +/- 0.05 | 3.68 +/- 0.40 | 7 | 26 |
| 300 | 4.33 +/- 0.11 | 3.45 +/- 0.05 | 13 | 31 |
| 360 | 4.44 +/- 0.08 | 3.37 +/- 0.06 | 11 | 33 |
| 420 | 4.44 +/- 0.13 | 3.05 +/- 0.07 | 11 | 39 |

(mmol/50 mL = 20 x mmol/L: 5.00 = 100 mmol/L; 3.05 = 61 mmol/L.) The pH 7 series is noisy at 30-120 min
(11, 5, 10, 11 %) and clean afterwards.

### Table 4. "Absolute and Relative Yields (mol %) of AcOH" (abs = per precursor charged; rel = per precursor consumed)

| time (h) | Glc/Pro pH 6 abs | rel | Glc/Pro pH 7 abs | rel | Fru-Pro pH 6 abs | rel | Fru-Pro pH 7 abs | rel |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | nd | — | 2.5 | 50 | 16.6 | 28 | 51.6 | 54 |
| 2 | nd | — | 5.4 | 50 | 17.8 | 25 | 54.6 | 55 |
| 4 | 1.0 | 13 | 9.7 | 37 | 24.5 | 30 | 58.7 | 59 |

### Table 5. "OAV of AcOH, HDMF, ATHP, and AP" — concentrations after 2 h, in mg/L, with the conversions this repo uses

Printed columns: concn (mg/L) and OAV per odorant. Added here: µg/mmol (= 10 x mg/L) and mol % of
precursor (= µg/mmol / (10 x MW)).

| pH | system | AcOH mg/L (OAV) | HDMF mg/L (OAV) | ATHP mg/L (OAV) | AP mg/L (OAV) |
|---:|---|---:|---:|---:|---:|
| 6 | Glc/Pro | 41 (<1) | <0.05 (<1) | 1.3 (810) | 0.25 (2500) |
| 6 | Fru-Pro | 648 (13) | 4 (67) | 0.2 (125) | 0.03 (300) |
| 7 | Glc/Pro | 325 (6.5) | 1.6 (27) | 4.5 (2810) | 0.36 (3600) |
| 7 | Fru-Pro | 3274 (65) | 20 (330) | 0.5 (310) | 0.04 (400) |
| 8 | Glc/Pro | 2310 (46) | 26 (430) | 1.1 (690) | 0.33 (3300) |
| 8 | Fru-Pro | 3810 (76) | 23 (380) | 1.3 (810) | 0.22 (2200) |

Converted (2 h, 100 mmol/L precursor):

| pH | system | AP µg/mmol | AP mol % | AP µmol/L | ATHP µg/mmol | ATHP mol % | ATHP µmol/L | HDMF µg/mmol | HDMF mol % | AcOH mg/mmol | AcOH mol % |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | Glc/Pro | 2.5 | 0.0022 | 2.2 | 13 | 0.0104 | 10.4 | <0.5 | <0.0004 | 0.41 | 0.68 |
| 6 | Fru-Pro | 0.3 | 0.00027 | 0.27 | 2 | 0.0016 | 1.6 | 40 | 0.031 | 6.48 | 10.8 |
| 7 | Glc/Pro | 3.6 | 0.0032 | 3.2 | 45 | 0.0360 | 36 | 16 | 0.012 | 3.25 | 5.4 |
| 7 | Fru-Pro | 0.4 | 0.00036 | 0.36 | 5 | 0.0040 | 4.0 | 200 | 0.156 | 32.7 | 54.5 |
| 8 | Glc/Pro | 3.3 | 0.0030 | 3.0 | 11 | 0.0088 | 8.8 | 260 | 0.203 | 23.1 | 38.5 |
| 8 | Fru-Pro | 2.2 | 0.0020 | 2.0 | 13 | 0.0104 | 10.4 | 230 | 0.180 | 38.1 | 63.4 |

Arithmetic: AP pH 7 Glc/Pro 0.36 mg/L x 10 = 3.6 µg/mmol; 3.6 / (10 x 111.14) = 0.0032 mol %;
0.36 / 111.14 = 3.2 µmol/L. ATHP 4.5 x 10 = 45; 45 / 1251.7 = 0.036 mol %. These reproduce the text's
"up to 0.004 mol %" (AP), "up to 0.04 mol %" (ATHP), "up to 0.2 mol %" (HDMF), "up to 60 mol %" (AcOH).
Table 4 cross-check: Glc/Pro pH 7 at 2 h 5.4 mol % = Table 5's 5.4; Fru-Pro pH 7 54.6 = 54.5. The two
pH 6 arms DISAGREE between tables (Table 4: Glc/Pro nd, Fru-Pro 17.8 mol %; Table 5: 0.68 and 10.8
mol %) — flag 5.

### Figure 5 (FIGURE-ONLY time courses, 1 / 2 / 4 h) — what the text states

- AcOH at 2 h: "about 400 to 38 000 µg/mmol"; Glc/Pro pH 7: 3.3 mg/mmol (2 h), 5.8 mg/mmol (4 h);
  Glc/Pro pH 8: "almost 20 mg/mmol already after 1 h"; more AcOH from Fru-Pro at every pH.
- HDMF at 2 h: "0 to about 250 µg/mmol"; favoured from Fru-Pro at pH 6 and 7; at pH 8 both systems
  comparable with similar kinetic curves; Fru-Pro jumps at pH 7, Glc/Pro needs pH 8.
- ATHP at 2 h: "2 to 45 µg/mmol"; higher from Glc/Pro than Fru-Pro at pH 6 and 7, similar at pH 8;
  **favoured at pH 7; "at pH 8 more ATHP is decomposed than formed"**.
- AP at 2 h: "about 0.5-4 µg/mmol"; favoured from Glc/Pro at all pH; similar in both systems at pH 8;
  Fru-Pro gives its highest AP at pH 8.

### Figure 6 (FIGURE-ONLY) — what the text states

- Fru-Pro formed from Glc/Pro: up to **1.3 mg per mmol glucose (~0.5 mol %; 1.3 / 277.27 = 0.47 %)**.
  pH 6: continuous rise, rapid within 2 h, plateau after ~4 h. pH 7: faster rise, **maximum at 1 h**,
  then decay to **half the maximum by 7 h**.
- Fed Fru-Pro (100 mmol/L): **< 10 % left after 1 h at pH 7; ~40 % left after 1 h at pH 6**.

## 4. Numbers and steps the repository can use

Registry keys: `hdmf` exists; 2-acetyl-1-pyrroline, 2-acetyltetrahydropyridine, acetic acid, proline,
Fru-Pro and 1-pyrroline are **not in registry** (`data/keys/compounds.yml`).

| quantity or step | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| AP level, Glc + Pro | 2.5 / 3.6 / 3.3 (0.0022 / 0.0032 / 0.0030 mol %) | µg/mmol Glc | pH 6 / 7 / 8, 100 mmol/L each, 0.2 M phosphate, reflux 2 h, IDA | Table 5 (mg/L x 10) | level_only (end of a 2 h cook; validation) |
| AP level, Fru-Pro | 0.3 / 0.4 / 2.2 (0.00027 / 0.00036 / 0.0020 mol %) | µg/mmol | same, 100 mmol/L Fru-Pro | Table 5 | level_only |
| ATHP level, Glc + Pro | 13 / 45 / 11 (0.0104 / 0.0360 / 0.0088 mol %) | µg/mmol Glc | pH 6 / 7 / 8, 2 h | Table 5 | level_only |
| ATHP level, Fru-Pro | 2 / 5 / 13 (0.0016 / 0.0040 / 0.0104 mol %) | µg/mmol | same | Table 5 | level_only |
| HDMF level, Glc + Pro | <0.5 / 16 / 260 (<0.0004 / 0.012 / 0.203 mol %) | µg/mmol Glc | pH 6 / 7 / 8, 2 h | Table 5 | level_only (a furaneol pH ladder with proline as the amine; B7 comparator) |
| HDMF level, Fru-Pro | 40 / 200 / 230 (0.031 / 0.156 / 0.180 mol %) | µg/mmol | same | Table 5 | level_only |
| AcOH yield, Glc + Pro pH 7 | 2.5 / 5.4 / 9.7 abs; 50 / 50 / 37 of glucose consumed | mol % | 1 / 2 / 4 h | Table 4 | level_only time series |
| AcOH yield, Fru-Pro pH 7 | 51.6 / 54.6 / 58.7 abs | mol % | 1 / 2 / 4 h | Table 4 | level_only (Amadori -> acetic acid conversion, ~55 mol %) |
| glucose loss with proline | 11 % (pH 6), 39 % (pH 7) at 420 min; full series above | mol % | 100 mmol/L each, reflux | Table 3 | level_only time series (measured) |
| apparent first-order glucose loss (mine, from Table 3 endpoints) | ~2.8e-4 (pH 6), ~1.2e-3 (pH 7) | 1/min | ln(5/4.44)/420; ln(5/3.05)/420 | derived | derived from level_only; orientation only |
| Fru-Pro maximum from Glc + Pro | ~1.3 mg/mmol Glc (~0.5 mol %) | — | pH 6 plateau after ~4 h; pH 7 maximum at 1 h | text on Figure 6A | figure_only (text-quoted) |
| fed Fru-Pro remaining at 1 h | < 10 % (pH 7); ~40 % (pH 6) | — | 100 mmol/L, reflux | text on Figure 6B | figure_only (text-quoted) |
| implied Fru-Pro decay (mine) | k >= 0.038 (pH 7); ~0.015 (pH 6); t1/2 <= 18 min, ~45 min | 1/min | from the two text-quoted fractions | derived | derived from figure_only; do not fit |
| ATHP pH ratios, Glc + Pro, 2 h | pH 7 / pH 6 = 3.5; pH 8 / pH 7 = 0.24 | — | same pot, same method | Table 5 | within_study_ratio |
| AP pH ratios, Glc + Pro, 2 h | pH 7 / pH 6 = 1.4; pH 8 / pH 7 = 0.92 | — | same | Table 5 | within_study_ratio (AP nearly pH-flat where ATHP peaks at 7) |
| ATHP : AP (molar), Glc + Pro, 2 h | 4.6 / 11.1 / 3.0 | — | pH 6 / 7 / 8 | Table 5 | within_study_ratio |
| Glc/Pro : Fru-Pro, 2 h | ATHP 6.5 / 9.0 / 0.85; AP 8.3 / 9.0 / 1.5 | — | pH 6 / 7 / 8 | Table 5 | within_study_ratio (the Amadori route is NOT the route to the N-heterocycles at pH <= 7) |
| HDMF : ATHP : AP (molar), Glc + Pro pH 8, 2 h | 68 : 3.0 : 1 | — | | Table 5 | within_study_ratio |
| literature cross-check | ATHP 38 µg/mmol Pro (Schieberle 1995) vs 45 here | µg/mmol | pH 7, 2 h, boiling | text | level_only, two laboratories agree within 20 % |
| fed-intermediate yields quoted from ref 16 | ~1 mol % ATHP (1-pyrroline + hydroxyacetone); ~5 mol % AP (1-pyrroline + 2-oxopropanal) | mol % | equimolar, pH 7, 30 min, boiling | Discussion | see `hofmann1998b_extraction.md` (0.9 % and 5.3 % there) |

**What the numbers say for a proline step.** (i) AP and ATHP are minor products of side reactions
limited by 1-pyrroline, which the authors could not quantify ("reliable quantitative data of IX under
Maillard reaction conditions are missing"). (ii) The two odorants respond to pH differently: ATHP
peaks at 7 and is lost at 8 ("more decomposed than formed"), AP is almost flat 6-8 from Glc/Pro. A
single shared pH law would miss this. (iii) The Amadori compound is the route to acetic acid and HDMF,
not to AP/ATHP; at pH 8 the distinction disappears because Glc/Pro forms and degrades Fru-Pro fast.
(iv) Roughly half of the glucose consumed at pH 7 becomes acetic acid (Table 4 rel 50 %), which is a
carbon-balance constraint the trunk's 2,3-enolisation branch could be checked against for a proline pot.

## 5. Flags

1. **Time courses are figure-only.** The dossier carries 2 h levels (Table 5) and the text-stated
   points; the 1 h and 4 h AP/ATHP values are not transcribed. Do not fit rates to this paper.
2. **Reflux, temperature unstated** (~100 C for 0.2 M phosphate); open system under a condenser.
3. **IDA on one tautomer of ATHP**; the labelled ATHP standard was 85 % pure, a tautomer and isotopomer
   mixture, calibrated against N-acetylpiperidine; the d-AP standard 95 % pure with 7 % of the
   2-acetyl-2-pyrroline tautomer. Levels inherit those calibrations.
4. **Reporting unit in the odorant table is mg/L**, not µg/mmol; the conversion factor 10 is exact for
   2 mmol in 20 mL and is the one the text's own quoted values reproduce.
5. **Table 4 and Table 5 disagree for the pH 6 arms**: acetic acid from Fru-Pro at 2 h is 17.8 mol % in
   Table 4 and 10.8 mol % (648 mg/L) in Table 5; Glc/Pro pH 6 is "nd" at 2 h in Table 4 and 41 mg/L
   (0.68 mol %; 410 µg/mmol, just above the 300 µg/mmol CV floor) in Table 5. The pH 7
   arms agree exactly. Treat the pH 6 acetic acid numbers as uncertain by a factor ~2.
6. **Fru-Pro pots drift a full pH unit** (Table 2); their "pH 6/7/8" labels are initial values only.
   Glc/Pro pots hold.
7. **1-pyrroline is never measured**; the proline pot is precursor-limited by an unquantified
   intermediate, so end-of-cook AP/ATHP levels constrain the product of (pyrroline supply) x (fed
   step), not the fed step alone (which Hofmann 1998b measures).
8. The odor thresholds behind the OAVs were determined in water without pH adjustment; AP and ATHP
   are bases whose headspace partition is pH-dependent.
9. Glucose consumption (Table 3) is with proline, a secondary amine that forms a different Amadori
   compound from the trunk's glycine; use it as a comparator, not a fit row, for the trunk's glucose loss.

# Alaiz, Hidalgo & Zamora 1999 — EXTRACTION
### *"Effect of pH and Temperature on Comparative Antioxidant Activity of Nonenzymatically Browned Proteins Produced by Reaction with Oxidized Lipids and Carbohydrates"*
### BSA (1 mg/mL) heated 24 h at 25–120 °C (pH 7) or at 80 °C and pH 4 / 7 / 10 (acetate / phosphate / borate, 0.3 M) alone (COBSA) or with 10 mM ribose (RIBSA), methyl-linoleate hydroperoxides (HPBSA) or methyl-linoleate secondary oxidation products (SPBSA); desalted, freeze-dried, added at **100 ppm to soybean oil, 10 g in a Petri dish, 60 °C in air, 264 h**; TBARS (nmol/g oil, figures only) and a BHT-referenced **protection index PI** (Tables 1–2), mean ± SD of three assays. **Every modified protein is antioxidant; ribose-browned BSA's PI rises with preparation temperature (6 at 25 °C → 85 at 120 °C) and peaks at pH 7 (36 / 66 / 53 at pH 4 / 7 / 10).** In TBARS terms the ribose effect is **−0.07 to −0.20 dex** against heated control BSA [Z]. A bulk-oil antioxidant readout, no aldehyde measured. Read 2026-11-15.

**Full read of `data/articles/alaiz1999.pdf`** (5 pp., born-digital ACS PDF, JAFC 47 (2), 748–752). Prose
and Tables 1–3 from `pdftotext -layout` (clean; "±" is rendered "(" and "=" as ")" in the text layer).
Tables checked against a 130 dpi render of p. 750–751. TBARS end points were read from **400 dpi crops**
of Figs. 4, 6 and 7 (p. 750) and are labelled as figure readings; Figs. 1–3 and 5 from the 130 dpi render.

**Provenance:** **[M]** measured, printed · **[Q]** qualitative · **[Z]** derived by me (incl. values
read off figures) · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Manuel Alaiz, Francisco J. Hidalgo, Rosario Zamora* |
| Affiliation | Instituto de la Grasa, CSIC, Seville |
| Title | "Effect of pH and Temperature on Comparative Antioxidant Activity of Nonenzymatically Browned Proteins Produced by Reaction with Oxidized Lipids and Carbohydrates" |
| Venue | *J. Agric. Food Chem.* **1999, 47 (2), 748–752** |
| DOI | **10.1021/jf980733r** (printed p. 748 footer, "CCC: $18.00"; Crossref 2026-11-15: same title, Alaiz / Hidalgo / Zamora, vol. 47, 748–752; matches the expected DOI) |
| Dates | received 8 Jul 1998; revised 5 Nov 1998; accepted 6 Nov 1998; published on Web 14 Jan 1999 |
| Funding | CICYT OLI96-2124, ALI97-0358; Junta de Andalucía AGR 0135 |
| Lab | Seville (Zamora / Hidalgo) — **the same group as `zamora2005`, `hidalgo2006`, the `zamora20xx` / `hidalgo20xx` series and `hidalgo1999` (the companion paper: same protein preparations)** [!] |
| Registry | not in `data/keys/papers.yml` before this dossier (2026-11-15) |

---

## 1. THE ONE-PARAGRAPH ANSWER

BSA modified by ribose, by linoleate hydroperoxides or by linoleate secondary oxidation products (the
proteins characterised in `hidalgo1999`) all slowed TBARS formation in soybean oil at 60 °C. The size of
the effect depends on how the protein was browned. At 25–50 °C the oxidised-lipid-modified proteins were
the better antioxidants; at 80–120 °C ribose-modified BSA was (PI 65.9 at 80 °C, 84.6 at 120 °C, where
PI = 100 is "as good as 100 ppm BHT"). PI correlates linearly with preparation temperature for RIBSA
(r = 0.977), HPBSA (0.965) and SPBSA (0.975). With pH of preparation at 80 °C, RIBSA peaks at **pH 7**
(PI 35.6 / 65.9 / 52.8 at pH 4 / 7 / 10), while HPBSA and SPBSA peak at pH 10. Across all preparations
PI tracks **lysine loss** best (RIBSA r = −0.984; Table 3), not fluorescence. The authors' reading:
Maillard and oxidised-lipid/protein browning give analogous, complementary antioxidant products, the
lipid route dominating at low temperature and the carbohydrate route at high temperature.

---

## 2. METHODS [M] (pp. 748–749)

| item | as printed |
|---|---|
| Protein preparation | BSA **2.5 mg in 2.5 mL** (1 mg/mL) of 0.3 M sodium acetate pH 4.0, sodium phosphate pH 7.4, or sodium borate pH 10.0; untreated (COBSA) or with **10 mM** ribose (RI), methyl-linoleate hydroperoxides (HP) or secondary products of methyl-linoleate oxidation (SP); **24 h at 25, 37, 50, 80 or 120 °C**; desalted (PD-10), freeze-dried. Details in `hidalgo1999` |
| Oil | soybean oil from the Instituto de la Grasa pilot plant, "stripped ... with no antioxidant" |
| Dose | proteins and BHT dissolved in the oil at **100 ppm** |
| Oxidation | **10 g oil** in 90 × 20 mm Petri dishes, **264 h, air, dark, 60 °C**; TBARS periodically (Kosugi et al. 1989), nmol TBARS/g oil |
| PI | **PI = 100 − 100·(TBARS_sample − TBARS_BHT)/(TBARS_oil − TBARS_BHT)** at the end of incubation; 100 = as BHT, 0 = no protection, < 0 = pro-oxidant |
| Stats | mean of three experiments (figures: mean ± SD of three assays); ANOVA + Student–Newman–Keuls, p < 0.05 |

---

## 3. RESULTS

### 3.1 Table 1 — PI vs preparation temperature, pH 7 [M] (mean ± SD; SNK letters within a row)

| T (°C) | COBSA | RIBSA | HPBSA | SPBSA |
|---|---|---|---|---|
| 25 | 1.33 ± 4.16 a | 6.31 ± 6.42 a | 16.78 ± 3.79 b | 25.75 ± 3.44 c |
| 37 | 4.03 ± 5.60 a | 8.96 ± 5.05 a | 27.20 ± 3.19 b | 32.73 ± 3.54 b |
| 50 | 6.86 ± 4.33 a | 29.38 ± 3.04 b | 40.30 ± 2.03 c | 47.83 ± 2.73 d |
| 80 | 7.84 ± 3.21 a | **65.94 ± 3.39 b** | 44.34 ± 1.61 c | 53.11 ± 1.95 d |
| 120 | 5.77 ± 5.51 a | **84.59 ± 4.48 b** | 62.69 ± 5.16 c | 72.86 ± 4.87 d |

(The text gives COBSA at 25 °C as 1.31; Table 1 prints 1.33.)

### 3.2 Table 2 — PI vs preparation pH, 80 °C [M]

| pH | COBSA | RIBSA | HPBSA | SPBSA |
|---|---|---|---|---|
| 4 | 9.43 ± 3.40 a | 35.59 ± 4.21 b | 46.26 ± 2.51 c | 57.53 ± 2.05 d |
| 7 | 7.84 ± 3.21 a | **65.94 ± 3.39 b** | 44.34 ± 1.61 c | 53.11 ± 1.95 d |
| 10 | 26.75 ± 4.51 a | 52.78 ± 3.40 b | 70.34 ± 5.33 c | 83.57 ± 3.60 d |

Printed also: SPBSA (80 °C, pH 10) PI is not different from RIBSA (120 °C, pH 7), p = 0.774. At 80 °C /
pH 7 "HPBSA reduced TBARS production ... by 44%" and RIBSA "by 66%", **both relative to the TBARS
decrease given by BHT** (i.e. these are the PIs, not TBARS ratios).

### 3.3 Table 3 — correlations over all pH/temperature preparations [M]

| correlation | RIBSA | HPBSA | SPBSA |
|---|---|---|---|
| PI / colour difference | r = 0.853, p = 0.015 | 0.784, 0.037 | 0.539, 0.212 |
| PI / fluorescence | 0.614, 0.143 | 0.368, 0.417 | 0.068, 0.885 |
| PI / lysine recovered | **−0.984, p < 0.0001** | −0.824, 0.023 | −0.787, 0.036 |
| PI / arginine recovered | −0.911, 0.0043 | −0.870, 0.011 | −0.742, 0.056 |

### 3.4 TBARS end points at 264 h — **figure only** [Z] (read from the 400 dpi crops, ±~50 nmol/g)

| figure (prep.) | oil | COBSA | RIBSA | HPBSA | SPBSA | BHT |
|---|---:|---:|---:|---:|---:|---:|
| Fig. 4 (80 °C, pH 7) | ~2000 | ~1920 | ~1210 | ~1470 | ~1350 | ~780 |
| Fig. 6 (80 °C, pH 4) | ~3250 | ~3080 | ~2600 | ~2300 | ~2200 | ~1400 |
| Fig. 7 (80 °C, pH 10) | ~3250 | ~2750 | ~2250 | ~1900 | ~1700 | ~1400 |

The oil-alone end point differs between figures (~2000 in Fig. 4, ~3250 in Figs. 1, 2, 6, 7): each
figure is a separate oxidation run. Cross-check: with r_BHT = 780/2000 = 0.39, RIBSA PI 65.94 predicts
TBARS ≈ 1200, matching the ~1210 read.

---

## 4. DERIVED [Z] — TBARS ratios in dex

From the PI definition, TBARS_x / TBARS_oil = 1 − (PI/100)·(1 − r_BHT), r_BHT = TBARS_BHT/TBARS_oil read
from the same figure (0.39 Fig. 4; 0.43 Figs. 6, 7). Printed PIs, figure-read r_BHT:

| preparation | RIBSA / oil | RIBSA / COBSA | SPBSA / COBSA | HPBSA / COBSA |
|---|---:|---:|---:|---:|
| 80 °C, pH 4 | −0.098 | **−0.074** | −0.148 | −0.108 |
| 80 °C, pH 7 | −0.223 | **−0.202** | −0.149 | −0.116 |
| 80 °C, pH 10 | −0.155 | **−0.083** | −0.209 | −0.150 |

Temperature series at pH 7 (r_BHT bracketed 0.35–0.50, since only Fig. 4 was read at 400 dpi):
RIBSA / COBSA = −0.01 (25 °C), −0.01 (37 °C), −0.05 to −0.07 (50 °C), −0.16 to −0.22 (80 °C),
**−0.23 to −0.33 (120 °C)** dex.

So: (i) ribose browning of a protein lowers TBARS in oil by at most ~0.2–0.3 dex at a 100 ppm dose,
after 11 days at 60 °C; (ii) the **pH of the browning** moves that effect by ~0.12 dex between pH 4 and
pH 7 (−0.07 vs −0.20) and ~0.12 dex between pH 7 and pH 10; (iii) the heated control protein itself is
near neutral (−0.02 dex at pH 4 and 7) except at pH 10 (−0.07). SDs on PI are 2–6 PI units, i.e. roughly
±0.01–0.03 dex on these ratios.

---

## 5. USE IN THIS REPO

**Question it was downloaded for (L1 scoping, `docs/reference/cases/L1_lipid_lever_coupling_scoping.md`):**
do Maillard reaction products act as antioxidants (or pro-oxidants) on lipid oxidation, and does pH
or sugar change that, enough to move hexanal / nonanal / 2-pentylfuran / dienals in the P3 pot?

- **Sign: antioxidant.** A Maillard-browned protein lowers lipid oxidation (TBARS) in a bulk oil; it
  never turned pro-oxidant here (all PI > 0). The effect scales with the extent of browning (lysine loss,
  r = −0.98), which supports "more Maillard → less lipid oxidation" as the MRP channel.
- **Magnitude: small.** The largest ribose effect is −0.2 to −0.3 dex on TBARS, after 264 h at 60 °C,
  in a bulk oil with 100 ppm of a pre-browned protein. TBARS is an aggregate (malondialdehyde-like
  carbonyls), not hexanal, and the readout is at a fixed late time in an autocatalytic run, which
  amplifies small rate differences.
- **pH lever.** RIBSA prepared at pH 7 protects more than at pH 4 or 10 (−0.20 vs −0.07 / −0.08 dex
  against COBSA). pH 6 vs 7 is not covered; the three pHs use three different buffers (acetate,
  phosphate, borate; borate complexes ribose), so pH is confounded with buffer [!]. Direction, if it
  carried over: pH 7 browning gives somewhat more antioxidant MRP than more acid pH, i.e. **hexanal
  slightly lower at pH 7** through this channel. Not a basis for a magnitude at pH 6 vs 7.
- **Sugar lever.** Only ribose; no glucose arm. Says nothing about glucose vs ribose, beyond "a more
  reactive sugar browns more and so protects more" (via the lysine-loss correlation).
- **Sulfur lever.** Not tested.
- **Distance from the P3 pot:** the MRP is made separately (24 h, 1 mg/mL BSA, 10 mM ribose) and tested
  in a different medium (bulk soybean oil, 60 °C, 11 days). In P3 the MRP forms in the same 20 min,
  140 °C aqueous pot in which the lipid oxidises, from 15 % pea isolate, so both formation time and
  contact time are far shorter. Not usable for a declared band.
- **Same-lab note [!].** Same Seville group as `zamora2005`, `hidalgo2006`, `hidalgo1999` and the
  zamora/hidalgo series. Under a leave-one-laboratory-out split these all fall in one fold. None of them
  is a designated hold-out (the hold-outs `farmer1992`, `elmore2002` are Reading).

---

## 6. FLAGS

1. **TBARS only in figures.** §3.4 values are my readings; §4 uses printed PIs plus a figure-read
   r_BHT. r_BHT uncertainty (±0.05) moves the RIBSA/COBSA ratios by ≤ 0.03 dex.
2. **Fig. 3 time axis** (50 °C preparations) runs only to 168 h, while Methods say every oil was
   oxidised for 264 h and the PIs are "at the end of the incubation period". Either the 50 °C run was
   shorter or the axis is truncated; PIs at 50 °C may not be at 264 h.
3. **Separate oil runs per figure**: oil-alone TBARS at the end point varies ~1.6× between figures; ratios
   are only meaningful within a figure.
4. **pH confounded with buffer** (acetate / phosphate / borate) and the "pH 7" arm is phosphate pH 7.4.
5. **COBSA 25 °C PI**: 1.31 in text, 1.33 in Table 1.
6. **Dose.** 100 ppm of protein in oil; the modified protein's lysine content is ~0.2–0.9 mmol/g
   (`hidalgo1999` Tables 1, 4), so this is ≤ 0.09 mmol amine-equivalent per kg oil.
7. **Registry:** not in `data/keys/papers.yml` before this dossier.

**Relation to existing dossiers:** `hidalgo1999_extraction.md` (companion: the same protein preparations,
browning and amino acid losses), `hidalgo2006_extraction.md` (same group; amines/PE as antioxidants in
oil, Rancimat), `zamora2005_extraction.md` (same group, review), `morales2001_extraction.md` (MRP radical
scavenging vs sugar).
